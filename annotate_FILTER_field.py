# The script is used to tag "PASS" or "FILTERED" to a given vcf file based on another vcf with such informations #
# Genotyping individuals at all sites where some of them were filtered during the VQSR #
# So the program expects two vcf file: (1) vcf file with its FILTER field unannotated, e.g. ".";
#									   (2) annotate sites who passed GATK VQSR with "PASS", and ones who didn't with "FILTERED"

import re
from argparse import ArgumentParser
from os.path import exists, dirname, realpath
from os import makedirs
from sys import exit, stdout, stderr
from datetime import datetime

def annotate_VCF(raw_vcf, vqsr_vcf, out_annotated_vcf) :

	sites_passed_VQSR_dict = {}
	all_sites_dict = {}
	num_passed_sites,  num_all_sites = 0, 0
	fVQSR_VCF = open(vqsr_vcf, 'r')
	for line in fVQSR_VCF :
		if not line.startswith('#') :
			tmp_line = re.split('\t', line.strip())
			chrom_id = tmp_line[0]
			chrom_pos = tmp_line[1]
			passed = tmp_line[6]
			if passed == "PASS" :
				sites_passed_VQSR_dict[chrom_id+':'+chrom_pos] = 1
				num_passed_sites += 1
				num_all_sites += 1
				"""
				debug purpose
				"""
				#if not all_sites_dict.has_key(chrom_id+':'+chrom_pos) :
				#	all_sites_dict[chrom_id+':'+chrom_pos] = 1
				#else :
				#	stdout.write(time_stamper() + " Warning: passed sites where different types of variants were identified: %s\n" %(chrom_id+':'+chrom_pos))
			else :
				num_all_sites += 1
				"""
				debug purpose
				"""
				#if not all_sites_dict.has_key(chrom_id+':'+chrom_pos) :
				#	all_sites_dict[chrom_id+':'+chrom_pos] = 1
				#else :
				#	stdout.write(time_stamper() + " Warning: filtered sites where different types of variants were identified: %s\n" %(chrom_id+':'+chrom_pos))
	fVQSR_VCF.close()
	stdout.write("%d sites passed GATK VQSR in %s\n" %(num_passed_sites, vqsr_vcf))
	stdout.write("%d sites filtered by GATK VQSR in %s\n" %(num_all_sites-num_passed_sites, vqsr_vcf)
	stdout.write("%d sites contained in %s\n" %(num_all_sites, vqsr_vcf))

	fOUT = open(out_annotated_vcf, 'w')

	num_filtered = 0
	fRAW_VCF = open(raw_vcf, 'r')
	for line in fRAW_VCF :
		if not line.startswith('#') :
			tmp_line = re.split('\t', line.strip())
			chrom_id = tmp_line[0]
			chrom_pos = tmp_line[1]
			if sites_passed_VQSR_dict.has_key(chrom_id+':'+chrom_pos) :
				fOUT.write('\t'.join(tmp_line[0:6])+"\tPASS\t"+'\t'.join(tmp_line[7:])+'\n')
			else :
				fOUT.write('\t'.join(tmp_line[0:6])+"\tFILTERED\t"+'\t'.join(tmp_line[7:])+'\n')
				num_filtered += 1
			"""
			debug purpose
			"""
			#if not all_sites_dict.has_key(chrom_id+':'+chrom_pos) :
			#	print chrom_id, chrom_pos
			#else :
			#	all_sites_dict.pop(chrom_id+':'+chrom_pos)
		else :
			fOUT.write(line)
	fRAW_VCF.close()
	fOUT.close()

def check_file_existence(*files) :
	""" check file existence """
	for file in files :
		if not exists(file) :
			stderr.write(time_stamper() + " Error: cannot find the file %s" %(file))
			exit()

def make_dir_if_needed(dir) :
	""" create directory if necessary """
	if not exists(dir) :
		makedirs(dir)

def time_stamper() :
	"""
	generate dash-separated time stamped string
	"""
	return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

if __name__ == "__main__" :
	parser = ArgumentParser(description="")
	parser.add_argument("-raw_vcf", metavar="FILE", dest="raw_vcf", required=True, help="VCF file with its \"FILTER\" field unannotated")
	parser.add_argument("-vqsr_vcf", metavar="FILE", dest="vqsr_vcf", required=True, help="VCF file obtained after GATK VQSR, with its \"FILTER\" field annotated with either \"PASS\" or SOMETHING ELSE")
	parser.add_argument("-out", metavar="FILE", dest="out_annotated_vcf", required=True, help="output VCF file with its \"FILTER\" field annotated")

	args = parser.parse_args()

	# check the existence of input files #
	check_file_existence(args.raw_vcf, args.vqsr_vcf)

	# set up the output, creat directory if necessary #
	make_dir_if_needed(dirname(realpath(args.out_annotated_vcf)))

	annotate_VCF(args.raw_vcf, args.vqsr_vcf, realpath(args.out_annotated_vcf))
