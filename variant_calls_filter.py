#Notes: change the way "-filter" argument was defined in next version

from sys import exit, stdout, stderr
from os.path import exists, dirname, realpath
from os import makedirs
from argparse import ArgumentParser
from datetime import datetime
import re

def apply_filter(raw_vcf_file, filtered_vcf_file, parsed_filter) :
	"""
	Parsing the raw VCF file and apply the pre-defined filter
	"""

	min_homo_alt, min_FORMAT_DP, max_het = parsed_filter

	fOUT = open(filtered_vcf_file, 'w')

	individuals = []
	num_hets, num_homo_ref, num_homo_alt = 0, 0, 0
	tot_var_sites, num_var_sites_left = 0, 0
	fVCF = open(raw_vcf_file, 'r')
	for line in fVCF :
		if re.match("^#", line.strip()) :
			if line.startswith("#CHROM") :
				individuals = re.split('\t', line.strip())[9:]
			fOUT.write(line)
		else :
			tot_var_sites += 1
			tmp_line = re.split('\t', line.strip())
			for i in range(9, len(tmp_line)) :
				gt_per_indv = re.split(':', tmp_line[i])[0]
				if gt_per_indv != "./." :
					cov_per_individual = int(re.split(':', tmp_line[i])[2])
					if gt_per_indv == "0/1" :
						num_hets += 1
					elif gt_per_indv == "1/1" and cov_per_individual >= int(min_FORMAT_DP) :
						num_homo_alt += 1
			if num_hets < int(max_het) and num_homo_alt >= int(min_homo_alt) :
				fOUT.write(line)
				num_var_sites_left += 1
			num_homo_alt, num_hets = 0, 0
			if tot_var_sites % 10000 == 0 :
				stdout.write("%d variants being processed\n" %(((tot_var_sites / 10000)) * 10000))
	stdout.write("%d variants being processed\n" %(((tot_var_sites % 10000))))
	fVCF.close()
	fOUT.close()
	stdout.write("%d variants left in the filtered VCF: %s\n" %(num_var_sites_left, realpath(filtered_vcf_file)))

def time_stamper() :
	"""
	generate a dash-separated time stamped string
	"""
	return datetime.now().strftime("%Y-%d-%m-%H-%M-%S")

if __name__ == "__main__" :
	parser = ArgumentParser(description="Apply pre-defined filter to raw VCF file.")
	parser.add_argument("-raw_vcf", metavar="FILE", dest="raw_vcf", required=True, help="raw VCF file to be filtered")
	parser.add_argument("-filtered_vcf", metavar="FILE", dest="filtered_vcf", required=True, help="filtered VCF file")
	parser.add_argument("-filter", metavar="A:B:C", dest="filter", help="sites (1) where at least A individuals that are homozygous for the non-reference base with a minimum of BX coverage, and (2) where at most C individuals is genotyped as heterozygous. Example: 1:4:1")

	args = parser.parse_args()

	parsed_filter = re.split(':', args.filter)
	if len(parsed_filter) != 3 :
		stderr.write(time_stamper() + "Error: incorrect argument value for -filter. Please use -h for details\n")
		exit()


	if not exists(args.raw_vcf) :
		stderr.write(time_stamper() + " Error: cannot find your raw VCF file %s\n" %(args.raw_vcf))
		exit()

	outdir = dirname(realpath(args.filtered_vcf))
	if not exists(outdir) :
		makedirs(outdir)
	out_filtered_vcf = realpath(args.filtered_vcf)

	apply_filter(args.raw_vcf, out_filtered_vcf, parsed_filter)
