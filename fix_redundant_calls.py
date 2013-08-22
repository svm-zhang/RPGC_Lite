import os,sys
import re

final_vqsr_call = sys.argv[1]
final_raw_call = sys.argv[2]
out_fixed_file = sys.argv[3]
fOUT = open(out_fixed_file, 'w')

variants_dict = {}
num_redundant = 0
num_kept = 0
fFINAL_VQSR = open(final_vqsr_call, 'r')
for line in fFINAL_VQSR :
	if line.startswith('#') :
		if line.startswith("#CHROM") :
			individuals = re.split('\t', line.strip())[9:]
			print len(individuals), individuals
	else :
		tmp_line = re.split('\t', line.strip())
		ref_base = tmp_line[3]
		alt_base = tmp_line[4]
		chrom_id = tmp_line[0]
		chrom_pos = tmp_line[1]
		qual = tmp_line[5]
		filter = tmp_line[6]					# PASS or FILTERED by VQSR #
		# fix sites having different types of calls: redundant calls #
		if not variants_dict.has_key(chrom_id+':'+chrom_pos) :
			variants_dict[chrom_id+':'+chrom_pos] = line.strip()
		else :
			num_redundant += 1
			same_site_diff_call = re.split('\t', variants_dict[chrom_id+':'+chrom_pos])
			print chrom_id, chrom_pos
			tmp_qual = same_site_diff_call[5]
			tmp_filter = same_site_diff_call[6]
			tmp_alt_base = same_site_diff_call[4]
			if (tmp_filter != "PASS" and filter != "PASS") or (filter == "PASS" and tmp_filter == "PASS") :		# if two different call both passed the VQSR or both not, we remove it from the final call set #	
				variants_dict.pop(chrom_id+':'+chrom_pos)
				if filter == "PASS" :
					print chrom_id, chrom_pos+"both pass"
				else :
					print chrom_id, chrom_pos+"both filtered"
			elif filter == "PASS" and tmp_filter != filter :
				print chrom_id, chrom_pos + "second kept"
				variants_dict.pop(chrom_id+':'+chrom_pos)
				num_kept += 1
				#variants_dict[chrom_id+':'+chrom_pos] = line.strip()
			elif tmp_filter == "PASS" and tmp_filter != filter :
				variants_dict.pop(chrom_id+':'+chrom_pos)
				print chrom_id, chrom_pos + "first kept"
				num_kept += 1
print num_redundant, num_kept
fFINAL_VQSR.close()

fFINAL_RAW = open(final_raw_call, 'r')
for line in fFINAL_RAW :
	if line.startswith("#") :
		fOUT.write(line)
	else :
		tmp_line = re.split('\t', line.strip())
		chrom_id = tmp_line[0]
		chrom_pos = tmp_line[1]
		if variants_dict.has_key(chrom_id+':'+chrom_pos) :
			fOUT.write(line)
fOUT.close()
fFINAL_RAW.close()
