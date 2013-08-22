import re
from argparse import ArgumentParser
from sys import exit, stdout, stderr
from os import makedirs
from os.path import exists, dirname, realpath
from datetime import datetime

def process_VCF(raw_vcf, targets_file, fixed_vcf = None) :

	fOUT = None
	if fixed_vcf is not None :
		fOUT = open(fixed_vcf, 'w')

	fTARGET = open(targets_file, 'w')

	num_targets = 0
	variants_dict = {}
	stdout.write("Parsing the VCF file and identifying potential targets\n")
	fRAW_VCF = open(raw_vcf, 'r')
	for line in fRAW_VCF :
		if line.startswith('#') :
			if line.startswith("#CHROM") :
				individuals = re.split('\t', line.strip())[9:]
			if fOUT :
				fOUT.write(line)
		else :
			tmp_line = re.split('\t', line.strip())
			ref_base = tmp_line[3]
			alt_base = tmp_line[4]
			chrom_id = tmp_line[0]
			chrom_pos = tmp_line[1]
			qual = tmp_line[5]
			filter = tmp_line[6]					# PASS or FILTERED by VQSR #
			if len(re.split(',', alt_base)) > 1 :
				num_refAlleles = 0
				for i in range(9, len(tmp_line)) :
					tmp_genotype = re.split(':', tmp_line[i])[0]
					if tmp_genotype == "0/1" or tmp_genotype == "0/0" or tmp_genotype == "0/2" :
						num_refAlleles += 1
				if num_refAlleles/float(len(individuals)) <= 0.05 :
					fTARGET.write(line)
					num_targets += 1
					ref_base = re.split(',', alt_base)[0]
					alt_base = re.split(',', alt_base)[1]
					variants_dict[chrom_id+':'+chrom_pos] = "%s\t%s\t%s\t%s" %('\t'.join(tmp_line[0:3]), ref_base, alt_base, '\t'.join(tmp_line[5:9]))
					if fOUT :
						fOUT.write("%s" %(variants_dict[chrom_id+':'+chrom_pos]))
					for i in range(9, len(tmp_line)) :
						if re.split(':', tmp_line[i])[0] == "1/1" :
							variants_dict[chrom_id+':'+chrom_pos] += "\t%s" %(re.sub("1/1", "0/0", tmp_line[i]))
							if fOUT :
								fOUT.write("\t%s" %(re.sub("1/1", "0/0", tmp_line[i])))
						elif re.split(':', tmp_line[i])[0] == "2/2" :
							variants_dict[chrom_id+':'+chrom_pos] += "\t%s" %(re.sub("2/2", "1/1", tmp_line[i]))
							if fOUT :
								fOUT.write("\t%s" %(re.sub("2/2", "1/1", tmp_line[i])))
						elif re.split(':', tmp_line[i])[0] == "1/2" :
							variants_dict[chrom_id+':'+chrom_pos] += "\t%s" %(re.sub("1/2", "0/1", tmp_line[i]))
							if fOUT :
								fOUT.write("\t%s" %(re.sub("1/2", "0/1", tmp_line[i])))
						else :
							variants_dict[chrom_id+':'+chrom_pos] += "\t./."
							if fOUT :
								fOUT.write("\t./.")
					if fOUT :
						fOUT.write('\n')
				else :
					if fOUT :
						fOUT.write(line)
			else :
				if fOUT :
					fOUT.write(line)
	fRAW_VCF.close()
	fTARGET.close()

	if fOUT :
		stdout.write("\t%d targets were targeted and fixed\n" %(num_targets))
		fOUT.close()
	else:
		stdout.write("\t%d targets were targeted\n" %(num_targets))

def timestamper() :
	"""
	generate dash-separated time stamped string
	"""
	return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

def make_dir_if_needed(dir) :
	""" make directory if necessary """
	if not exists(dir) :
		os.makedirs(dir)

if __name__ == "__main__" :
	parser = ArgumentParser(description="Fix variant sites with three alleles called originally by GATK UG, but only two of them actually show up in the genotypes")
	parser.add_argument("-raw_vcf", metavar="FILE", dest="raw_vcf", required=True, help="input VCF file to be fixed")
	parser.add_argument("-fixed_vcf", metavar="FILE", dest="fixed_vcf", help="output fixed VCF file. Not required if -fixoff is specified")
	parser.add_argument("-targets", metavar="FILE", dest="targets_file", required=True, help="output file with a set of targets identified from the input VCF file")
	parser.add_argument("-fixoff", action="store_true", dest="fix_or_not", help="turn off fixing identified target sites")

	args = parser.parse_args()

	if not exists(args.raw_vcf) :
		stderr.write(timestamper() + " Error: cannot find the VCF file you provided: %s\n" %(args.raw_vcf))
		exit()

	make_dir_if_needed(dirname(realpath(args.fixed_vcf)))

	if not args.fix_or_not :
		if not args.fixed_vcf :
			stderr.write(timestamper() + " Error: you need to provide a path to the output fixed VCF file\n")
			exit()
		else :
			make_dir_if_needed(dirname(realpath(args.fixed_vcf)))
			process_VCF(args.raw_vcf, realpath(args.targets_file), realpath(args.fixed_vcf))
	else :
		process_VCF(args.raw_vcf, realpath(args.targets_file))
