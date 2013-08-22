import os, sys
import re

filtered_vcf_file = sys.argv[1]
training_call_set_outprefix = sys.argv[2]
top_percent_good_snps = int(sys.argv[3])

training_snps_outfile = training_call_set_outprefix + ".snps.train_1.vcf"
training_indels_outfile = training_call_set_outprefix + ".indels.train_1.vcf"
fTRAIN_SNP = open(training_snps_outfile, 'w')
fTRAIN_INDEL = open(training_indels_outfile, 'w')

variants_dict = {}
snps_call_score_dict = {}
indels_call_score_dict = {}
num_variants, num_snps, num_indels = 0, 0, 0
headers = ""
fFVCF = open(filtered_vcf_file, 'r')
for line in fFVCF :
	if re.match("^#", line.strip()) :
		if line.startswith("#CHROM") :
			individuals = re.split('\t', line.strip())[9:]
			print len(individuals)
		fTRAIN_SNP.write(line)
		fTRAIN_INDEL.write(line)
	else :
		tmp_line = re.split('\t', line.strip())
		chrom_id = tmp_line[0]
		chrom_pos = tmp_line[1]
		ref_base = tmp_line[3]
		alt_base = tmp_line[4]
		variants_dict[chrom_id+":"+chrom_pos] = line.strip()
		qual = float(tmp_line[5])
		if re.findall("DP=\d+", tmp_line[7]) :
			dp = re.findall("DP=\d+", tmp_line[7])[0]
			if int(re.split('=', dp)[1]) >= 40 and int(re.split('=', dp)[1]) <= 60 :
				if len(re.split(',', alt_base)) == 1 :
					########
					# SNPs #
					########
					if len(ref_base) == len(alt_base) and len(ref_base) == 1 :
						num_het, num_homo_nonRef, num_homo_Ref = 0, 0, 0
						tmp_num_individual = len(individuals)
						for i in range(9, len(tmp_line)) :
							if re.split(':', tmp_line[i])[0] == "1/1" and int(re.split(':', tmp_line[i])[2]) >= 4 :
								num_homo_nonRef += 1
							#if re.split(':', tmp_line[i])[0] == "0/0" and int(re.split(':', tmp_line[i])[2]) >= 4 :
							#	num_homo_Ref += 1
							#if re.split(':', tmp_line[i])[0] == "./." :
							#	tmp_num_individual -= 1
							if re.split(':', tmp_line[i])[0] == "0/1" :
								num_het += 1
						if num_homo_nonRef >= 1 and num_het <= 1 :
							#fTRAIN_SNP.write(line)
							if not snps_call_score_dict.has_key(qual) :
								snps_call_score_dict[qual] = [chrom_id+":"+chrom_pos]
							else :
								snps_call_score_dict[qual].append(chrom_id+":"+chrom_pos)
							num_snps += 1
					##########
					# INDELs #
					##########
					elif len(ref_base) > len(alt_base) or len(ref_base) < len(alt_base) :
						num_het, num_homo_nonRef = 0, 0
						for i in range(9, len(tmp_line)) :
							if re.split(':', tmp_line[i])[0] == "1/1" and int(re.split(':', tmp_line[i])[2]) >= 4 :
								num_homo_nonRef += 1
							if re.split(':', tmp_line[i])[0] == "0/1" :
								num_het += 1
						if num_homo_nonRef >= 1 and num_het <= 1 :
							fTRAIN_INDEL.write(line)
							if not indels_call_score_dict.has_key(qual) :
								indels_call_score_dict[qual] = [chrom_id+":"+chrom_pos]
							else :
								indels_call_score_dict[qual].append(chrom_id+":"+chrom_pos)
							num_indels += 1
		num_variants += 1
fTRAIN_INDEL.close()
print num_variants
print num_snps, num_indels

num_high_scoring_calls = 0.0
high_number_calls = []
for score in sorted(snps_call_score_dict.iterkeys(), reverse=True) :
	#if num_high_scoring_calls == 0.0 :
	#	print snps_call_score_dict[score]
	for i in range(len(snps_call_score_dict[score])) :
		num_high_scoring_calls += 1
		if num_high_scoring_calls/num_snps <= (top_percent_good_snps/100.0) :
			high_number_calls.append(snps_call_score_dict[score][i])

print len(high_number_calls)

tmp_high_number_calls = high_number_calls
fFVCF.seek(0)
for line in fFVCF :
	if re.match("^#", line.strip()) :
		continue
	else :
		tmp_line = re.split('\t', line.strip())
		chrom_id = tmp_line[0]
		chrom_pos = tmp_line[1]
		
		if chrom_id+":"+chrom_pos in tmp_high_number_calls :
			fTRAIN_SNP.write(line)
			tmp_high_number_calls.pop(tmp_high_number_calls.index(chrom_id+':'+chrom_pos))
fFVCF.close()
fTRAIN_SNP.close()
