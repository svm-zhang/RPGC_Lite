import os,sys
import re
from argparse import ArgumentParser
from sys import exit, stdout, stderr
from os import makedirs
from os.path import exists, realpath, dirname

def process_VCF(input_vcf, out_dup_file, out_vcf = None) :
	"""
	identify duplicates and remove them if -fixoff is not specified
	"""

	fVCF_OUT = None
	if out_vcf is not None :
		fVCF_OUT = open(out_vcf, 'w')
	fDUP_OUT = open(out_dup_file, 'w')

	variants_dict = {}
	variants_list = []
	num_redundant, num_kept = 0, 0
	fINVCF = open(input_vcf, 'r')
	for line in fINVCF :
		if line.startswith('#') :
			if line.startswith("#CHROM") :
				individuals = re.split('\t', line.strip())[9:]
				stdout.write("%d individuals included in the VCF file: %s\n" %(len(individuals), input_vcf))
			if fVCF_OUT :
				fVCF_OUT.write(line)
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
				variants_list.append(chrom_id+':'+chrom_pos)
			else :
				num_redundant += 1
				same_site_diff_call = re.split('\t', variants_dict[chrom_id+':'+chrom_pos])
				tmp_qual = same_site_diff_call[5]
				tmp_filter = same_site_diff_call[6]
				tmp_alt_base = same_site_diff_call[4]
				fDUP_OUT.write("%s\n%s\n" %(variants_dict[chrom_id+':'+chrom_pos], line.strip()))
				if (tmp_filter != "PASS" and filter != "PASS") or (filter == "PASS" and tmp_filter == "PASS") :		# if two different call both passed the VQSR or both not, we remove it from the final call set #	
					variants_dict.pop(chrom_id+':'+chrom_pos)
					variants_list.remove(chrom_id+':'+chrom_pos)
					if filter == "PASS" :
						stdout.write(chrom_id+" "+chrom_pos+" both pass\n")
					else :
						stdout.write(chrom_id+" "+chrom_pos+" both filtered\n")
				elif filter == "PASS" and tmp_filter != filter :
					stdout.write(chrom_id+" "+chrom_pos + " second kept\n")
					variants_dict[chrom_id+':'+chrom_pos] = line.strip()
					num_kept += 1
				elif tmp_filter == "PASS" and tmp_filter != filter :
					stdout.write(chrom_id+" "+chrom_pos + " first kept\n")
					num_kept += 1
	stdout.write("%d\t%d\n" %(num_redundant, num_kept))

	if fVCF_OUT :
		for i in range(len(variants_list)) :
			fVCF_OUT.write("%s\n" %(variants_dict[variants_list[i]]))
		fVCF_OUT.close()
	fINVCF.close()

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
	parser = ArgumentParser(description="Identify and remove variants sites where two types of variants (SNP, INDEL) were identified.")
	parser.add_argument("-vcf", metavar="FILE", dest="vcf_file", required=True, help="input VCF file")
	parser.add_argument("-duplicates", metavar="FILE", dest="out_dup_file", required=True, help="output plain text with duplicated items identified, Whether or not -fixoff is specified, the program will spit out this file")
	parser.add_argument("-out", metavar="FILE", dest="out_file", help="output VCF file with duplicates removed. Not required if -fixoff is specified")
	parser.add_argument("-fixoff", action="store_true", dest="fix_or_not", help="turn off removing identified duplicates from input VCF file. No output VCF will be expected if this is specified")

	args = parser.parse_args()

	# file path sanity check #
	if not exists(args.vcf_file) :
		stderr.write(time_stamper() + " Error: cannot find the VCF file you provided: %s" %(args.vcf_file))
		exit()

	make_dir_if_needed(dirname(realpath(args.out_dup_file)))
	out_dup_file = realpath(args.out_dup_file)

	fix = -1
	if args.fix_or_not :
		fix = 0
	else :
		fix = 1

	if fix :
		if not args.out_file :
			stderr.write(time_stamper() + " Error: you need to provide a path to the output VCF file with duplicate removed, since you specified -fix\n")
			exit()
		else :
			make_dir_if_needed(dirname(realpath(args.out_file)))
			out_file = realpath(args.out_file)
			process_VCF(args.vcf_file, out_dup_file, out_file)
	else :
		process_VCF(args.vcf_file, out_dup_file)
