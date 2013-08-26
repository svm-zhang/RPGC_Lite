# Note: 1. only BWA supported
#	2. only BWA default command line setting supported
#	3. for each dataset, only one pair of reads and one single-end reads file allowed

import re, shlex
from argparse import ArgumentParser
from os import makedirs, devnull, system
from os.path import exists, join, dirname, realpath
from sys import exit, stdout, stderr, argv
from datetime import datetime
from subprocess import Popen

def get_data(config_file) :
	""" get all the fastq files from the configuration file """
	stdout.write(timestamper() + " Getting data from the configuration file\n")
	pe_dict, se_dict = {}, {}
	fCONFIG = open(config_file, 'r')
	line_id = ""
	ndataset = 0
	for line in fCONFIG :
		if not line.startswith('#') and line.strip() != "" :
			tmp_line = re.split("\s+", line.strip())
			line_id = tmp_line[0]
			""" in case of single end reads available """
			if len(tmp_line) == 4 :
				se_dict[line_id] = tmp_line[3]
			elif len(tmp_line) > 4 :	# in case where data is provided more than it should be# 
				stderr.write(timestamper() + " [Config Error] only one pair of reads file and one single-end file can be provided in each line\n")
				stderr.write(timestamper() + " [Config Error] %s\n" %(line.strip()))
				exit()
			fastq_1 = tmp_line[1]
			fastq_2 = tmp_line[2]
			pe_dict[line_id] = "%s|%s" %(fastq_1, fastq_2)
			ndataset += 1
	stdout.write(timestamper() + " %d datasets parsed\n" %(ndataset))
	fCONFIG.close()
	return pe_dict, se_dict


def timestamper() :
	""" generate a dash-separated time stamped string """
	return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

def check_files_existence(*files) :
	""" check file existence """
	for file in files :
		if not exists(file) :
			stderr.write(timestamper() + " [IO Error] Cannot find the file you provided: %s\n" %(file))
			exit()

def make_dirs_if_needed(*dirs) :
	""" make any directories when necessary """
	for dir in dirs :
		if not exists(dir) :
			makedirs(dir)

def call_bwa(bwa_cmd, out = None, err = None) :
	""" call BWA """
	if err is None :
		fDEVNULL = open(devnull, 'wb')
		proc = Popen(shlex.split(bwa_cmd), stderr=fDEVNULL, stdout=file(out, 'w'), shell=False)
	else :
		proc = Popen(shlex.split(bwa_cmd), stderr=file(err, 'w'), stdout=file(out, 'w'), shell=False)
	proc_stderr = proc.communicate()[1]
	exit()
	if proc.returncode != 0 :
		stderr.write(timestamper() + " [BWA Error] something wrong when running %s\n" %(bwa_cmd))
		stderr.write(timestamper() + " [BWA Error] %s\n" %(proc_stderr))
		exit()

if __name__ == "__main__" :
	parser = ArgumentParser(description="")
	parser.add_argument("-config", metavar="FILE", dest="config_file", required=True, help="Specify a configuration file with a list of fastq files to be processed")
	parser.add_argument("-db", metavar="PREFIX", dest="refDB_prefix", required=True, help="Specify the DB prefix you want to run BWA against. Currently, no sanity check on the existence of the DB you provided")
	parser.add_argument("-outdir", metavar="DIR", dest="outdir", required=True, help="Specify the output directory")
	parser.add_argument("-t", metavar="INT", dest="num_threads", type=int, default=4, help="Specify the number of threads you want to run BWA with")
	parser.add_argument("-mapper_log", metavar="FILE", dest="mapper_log", help="Specify a log for recording mapper output")
	args = parser.parse_args()

	check_files_existence(args.config_file)
	make_dirs_if_needed(args.outdir)

	pe_dict, se_dict = get_data(args.config_file)

	#fMAPPER_LOG = None
	#if args.mapper_log is not None :
	#	if exists(args.mapper_log) :
	#		system("rm %s" %(args.mapper_log))
	#	fMAPPER_LOG = open(args.mapper_log, 'w')
	num_threads = args.num_threads
	refDB_prefix = args.refDB_prefix
	for each_dataset in pe_dict.iterkeys() :
		"""
		pair-end reads mapping
		"""
		stdout.write(timestamper() + " running BWA for pair-end reads from dataset %s\n" %(each_dataset))
		tmp_data = re.split('\|', pe_dict[each_dataset])
		fastq_1 = realpath(tmp_data[0])
		fastq_2 = realpath(tmp_data[1])
		check_files_existence(fastq_1, fastq_2)
		outdir_per_dataset = join(args.outdir, "%s" %(each_dataset))
		make_dirs_if_needed(outdir_per_dataset)
		sai_1 = join(outdir_per_dataset, "%s_1.sai" %(each_dataset))
		sai_2 = join(outdir_per_dataset, "%s_2.sai" %(each_dataset))
		pe_sam = join(outdir_per_dataset, "%s.pe.sam" %(each_dataset))
		bwa_aln_cmd = "/N/soft/mason/bwa/bwa-0.6.2-gcc-4.4.6/bwa aln -t %d %s" %(num_threads, refDB_prefix)
		call_bwa(bwa_aln_cmd + " %s" %(fastq_1), sai_1)
		call_bwa(bwa_aln_cmd + " %s" %(fastq_2), sai_2)
		bwa_sampe_cmd = "/N/soft/mason/bwa/bwa-0.6.2-gcc-4.4.6/bwa sampe %s %s %s %s %s" %(refDB_prefix, sai_1, sai_2, fastq_1, fastq_2)
		call_bwa(bwa_sampe_cmd, pe_sam, args.mapper_log)

		"""
		single-end reads mapping
		"""
		if se_dict.has_key(each_dataset) :
			stdout.write(timestamper() + " running BWA for single-end reads from dataset %s\n" %(each_dataset))
			se_fastq = realpath(se_fastq[each_dataset])
			check_files_existence(se_fastq)
			se_sai = join(outdir_per_dataset, "%s_se.sai" %(each_dataset))
			se_sam = join(outdir_per_dataset, "%s.pe.sam" %(each_dataset))
			call_bwa(bwa_aln_cmd + " %s" %(se_fastq), se_sai)
			bwa_samse_cmd = "/N/soft/mason/bwa/bwa-0.6.2-gcc-4.4.6/bwa samse %s %s %s" %(refDB_prefix, se_sai, se_fastq)
			call_bwa(bwa_samse_cmd, se_sam, args.mapper_log)
	stdout.write(timestamper() + "finish\n")
