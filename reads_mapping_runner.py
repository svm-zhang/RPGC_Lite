# The script is used to map pair-end (fragment libs) reads (call bwa) from each RIL line against the assembled scaffolds sequences #
# for current version, only bwa is supported #

import os, sys
import re
import getopt

def usage() :
	print "python rils_mapping.py [options] ..."

def file_existence_check(infile) :
	if not os.path.exists(infile) :
		print "Error: cannot find the file you specify: %s" %(infile)
		sys.exit()

def parse_reads_list(reads_list_file) :
	reads_file_dict = {}
	reads_list_handle = open(reads_list_file, 'r')
	for line in reads_list_handle :
		if len(re.split(' +', line.strip())) != 3 :
			print "Error: wrong format for the reads list file"
			print "Example: [line_id]\t[*_1.fastq]\t[*_2.fastq]"
			sys.exit(1)
	reads_list_handle.seek(0)
	line_id = ""
	num_line = 0
	for line in reads_list_handle :
		if not line.startswith('#') and line.strip() != "" :
			line_id = re.split(' +', line.strip())[0]
			fastq_1 = re.split(' +', line.strip())[1]
			fastq_2 = re.split(' +', line.strip())[2]
			reads_file_dict[line_id] = "%s\t%s" %(fastq_1, fastq_2)
			num_line += 1
	print num_line
	reads_list_handle.close()
	return reads_file_dict

def main() :
	try :
		opts, args = getopt.getopt(sys.argv[1:], "hf:l:o:t:", ["help", "--refDB_prefix", "--reads_list", "--outdir", "--threads"])
	except getopt.GetoptError, err :
		print str(err)
		usag()
		sys.exit(2)

	reads_list_file, outdir, refDB_prefix = "", "", ""
	num_threads = 8
	for o, a in opts :
		if o in ("-h", "--help") :
			usage()
			sys.exit(2)
		elif o in ("-f", "--refDB_prefix") :
			refDB_prefix = a
		elif o in ("-l", "--reads_list") :
			reads_list_file = a
		elif o in ("-o", "--outdir") :
			outdir = a
		elif o in ("-t", "--threads") :
			num_threads = int(a)
		else :
			print "Error: %s option cannot be recognized" %(o)
			print "Please try python rils_mapping.py -h or --help for more command information"
			sys.exit()

	# parse the reads list and at same time check the file format #
	reads_file_dict = parse_reads_list(reads_list_file)

	# map each of the RIL lines against the assembled genomes #
	for ril_line in reads_file_dict.iterkeys() :
		fastq_1 = re.split('\t', reads_file_dict[ril_line])[0]
		fastq_2 = re.split('\t', reads_file_dict[ril_line])[1]
		file_existence_check(fastq_1)
		file_existence_check(fastq_2)
		ril_outdir = os.path.join(outdir, "%s" %(ril_line))
		if not os.path.exists(ril_outdir) :
			os.makedirs(ril_outdir)
		sai_1 = os.path.join(ril_outdir, "%s_1.sai" %(ril_line))
		sai_2 = os.path.join(ril_outdir, "%s_2.sai" %(ril_line))
		sam_out = os.path.join(ril_outdir, "%s.sam" %(ril_line))
		bwa_aln_cmd = "bwa aln -t %d %s" %(num_threads, refDB_prefix)
		os.system("%s %s > %s" %(bwa_aln_cmd, fastq_1, sai_1))
		os.system("%s %s > %s" %(bwa_aln_cmd, fastq_2, sai_2))
		bwa_sampe_cmd = "bwa sampe %s %s %s %s %s > %s" %(refDB_prefix, sai_1, sai_2, fastq_1, fastq_2, sam_out)
		print bwa_sampe_cmd
		os.system(bwa_sampe_cmd)

if __name__ == "__main__" :
	main()
