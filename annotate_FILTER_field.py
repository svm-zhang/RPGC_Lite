# The script is used to tag "PASS" or "Filtered" to a given vcf file based on another vcf with such informations #
# Genotyping individuals at all sites where some of them were filtered during the VQSR #
# So the program expects two vcf file: (1) vcf file from genotyping individuals at all raw sites;
#									   (2) annotate sites who passed the VQSR with "PASS", and ones who didn't with "FILTERED"

import os, sys
import re

genotyping_at_all_sites_vcf_file = sys.argv[1]
vcf_with_sites_passed_VQSR = sys.argv[2]
out_tagged_genotyping_vcf = sys.argv[3]

sites_passed_VQSR_dict = {}
all_sites_dict = {}
num_passed_sites,  num_all_sites = 0, 0
fPASSED_VQSE_VCF = open(vcf_with_sites_passed_VQSR, 'r')
for line in fPASSED_VQSE_VCF :
	if not line.startswith('#') :
		tmp_line = re.split('\t', line.strip())
		chrom_id = tmp_line[0]
		chrom_pos = tmp_line[1]
		passed = tmp_line[6]
		if passed == "PASS" :
			sites_passed_VQSR_dict[chrom_id+':'+chrom_pos] = 1
			if not all_sites_dict.has_key(chrom_id+':'+chrom_pos) :
				all_sites_dict[chrom_id+':'+chrom_pos] = 1
			else :
				print chrom_id+':'+chrom_pos
			num_passed_sites += 1
			num_all_sites += 1
		else :
			if not all_sites_dict.has_key(chrom_id+':'+chrom_pos) :
				all_sites_dict[chrom_id+':'+chrom_pos] = 1
			else :
				print chrom_id+':'+chrom_pos
			num_all_sites += 1
fPASSED_VQSE_VCF.close()
print num_passed_sites, num_all_sites
print len(all_sites_dict)

fOUT = open(out_tagged_genotyping_vcf, 'w')

num = 0
fGT_VCF = open(genotyping_at_all_sites_vcf_file, 'r')
for line in fGT_VCF :
	if not line.startswith('#') :
		tmp_line = re.split('\t', line.strip())
		chrom_id = tmp_line[0]
		chrom_pos = tmp_line[1]
		if sites_passed_VQSR_dict.has_key(chrom_id+':'+chrom_pos) :
			fOUT.write('\t'.join(tmp_line[0:6])+"\tPASS\t"+'\t'.join(tmp_line[7:])+'\n')
		else :
			fOUT.write('\t'.join(tmp_line[0:6])+"\tFILTERED\t"+'\t'.join(tmp_line[7:])+'\n')
		if not all_sites_dict.has_key(chrom_id+':'+chrom_pos) :
			print chrom_id, chrom_pos
		else :
			all_sites_dict.pop(chrom_id+':'+chrom_pos)
		num += 1
	else :
		fOUT.write(line)
print num

for key, value in all_sites_dict.iteritems() :
	print key
