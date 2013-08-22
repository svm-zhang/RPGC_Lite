import re
from argparse import ArgumentParser
from os import makedirs
from os.path import join, exists, dirname, realpath
from sys import exit, stdout, stderr
from datetime import datetime

def get_training_calls(vcf_file, outprefix_training_file, top_percent) :
	training_snps_outfile = outprefix_training_file + ".snp.vcf"
	training_indels_outfile = outprefix_training_file + ".indels.vcf"
	fTRAIN_SNP = open(training_snps_outfile, 'w')
	fTRAIN_INDEL = open(training_indels_outfile, 'w')

	variants_dict = {}
	snps_call_score_dict = {}
	indels_call_score_dict = {}
	num_variants, num_snps, num_indels = 0, 0, 0
	headers = ""
	stdout.write("Parsing the VCF file: %s\n" %(vcf_file))
	fVCF = open(vcf_file, 'r')
	for line in fVCF :
		if line.startswith("#") :
			if line.startswith("#CHROM") :
				individuals = re.split('\t', line.strip())[9:]
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
							if not snps_call_score_dict.has_key(qual) :
								snps_call_score_dict[qual] = [chrom_id+":"+chrom_pos]
							else :
								snps_call_score_dict[qual].append(chrom_id+":"+chrom_pos)
							num_snps += 1
						##########
						# INDELs #
						##########
						elif len(ref_base) > len(alt_base) or len(ref_base) < len(alt_base) :
							fTRAIN_INDEL.write(line)
							if not indels_call_score_dict.has_key(qual) :
								indels_call_score_dict[qual] = [chrom_id+":"+chrom_pos]
							else :
								indels_call_score_dict[qual].append(chrom_id+":"+chrom_pos)
							num_indels += 1
			num_variants += 1
	fTRAIN_INDEL.close()
	stdout.write("%d variants in %s\n" %(num_variants, vcf_file))
	stdout.write("\t%d SNPs parsed out\n" %(num_snps))
	stdout.write("\t%d INDELS parsed out\n" %(num_indels))

	stdout.write("Output the top %d%% of the variants\n" %(top_percent))
	num_high_scoring_calls = 0.0
	high_number_calls = []
	for score in sorted(snps_call_score_dict.iterkeys(), reverse=True) :
		for i in range(len(snps_call_score_dict[score])) :
			num_high_scoring_calls += 1
			if num_high_scoring_calls/num_snps <= (top_percent/100.0) :
				high_number_calls.append(snps_call_score_dict[score][i])

	tmp_high_number_calls = high_number_calls
	fVCF.seek(0)
	for line in fVCF :
		if not line.startswith("#") :
			tmp_line = re.split('\t', line.strip())
			chrom_id = tmp_line[0]
			chrom_pos = tmp_line[1]
			
			if chrom_id+":"+chrom_pos in tmp_high_number_calls :
				fTRAIN_SNP.write(line)
				tmp_high_number_calls.pop(tmp_high_number_calls.index(chrom_id+':'+chrom_pos))
	fVCF.close()
	fTRAIN_SNP.close()

def time_stamper() :
	"""
	generate a dash-separated time stamped string
	"""
	return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

def make_dir_if_needed(dir) :
	""" Make a single directory if it does not exist. """
	if not exists(dir) :
		os.makedirs(dir)

if __name__ == "__main__" :
	parser = ArgumentParser(description="Get the top X percent of the highest-scoring calls and use it as the input to GATK VQSR")
	parser.add_argument("-vcf", metavar="FILE", required=True, dest="vcf_file", help="input VCF file obtained from variant_calls_filter.py")
	parser.add_argument("-out", metavar="PREFIX", required=True, dest="outprefix_training_file", help="output prefix of VCF file used for training")
	parser.add_argument("-percent", metavar="INT", dest="top_percent", type=int, default=10,  help="user specified top X percent of the hightest-scoring calls")

	args = parser.parse_args()

	if not exists(realpath(args.vcf_file)) :
		stderr.write(time_stamper() + " Error: cannot find the file you provided: %s" %(realpath(args.vcf_file)))
		exit()

	make_dir_if_needed(dirname(realpath(args.outprefix_training_file)))

	get_training_calls(realpath(args.vcf_file), realpath(args.outprefix_training_file), args.top_percent)
