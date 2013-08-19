# Author: Simo Zhang
# Email: simozhan@indiana.edu
# Usage: python post_mapping_pipeline.py -h #
# version: 1.01

#Updates (8/8/2013):
#		1. add option "-picard_log" to record log from Picard STDERR
#		2. add logging module with customized filter to keep track records at different logging level
#		3. add Call_Picard_SubProgram(cmd, logger_std, logger_prog)

# Future work:
#				1. Complete -sam and -RG options
#				2. Hope to add in break-and-continue function so that no need to re-run entire pipeline everytime something come up at later steps
#				3. Perhaps add a functio to spit out an example of config file
#				4. Perhaps change the format of config file

import os, sys
import re
import argparse
import logging
import shlex, subprocess
from subprocess import Popen, PIPE
from datetime import datetime


class Logging_Level_Filter(logging.Filter) :
	def __init__(self, level) :
		self.level = level

	def filter(self, record) :
		return (record.levelno == self.level)

# call Picard CleanSam #
def Clean_SAM_Files(raw_samfile, prog, rgid, logger_std, logger_prog) :
	cleaned_samfile = '.'.join(re.split('\.', raw_samfile)[0:3]) + ".cleaned.sam"
	picard_clean_sam_cmd = "java -Xmx2g -jar %s INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=STRICT" %(os.path.join(prog, "CleanSam.jar"), raw_samfile, cleaned_samfile)
	logger_std.info("[Run Picard] - [%s INPUT]: %s" %(rgid, raw_samfile))
	logger_std.info("[Run Picard] - [%s CleanSam]: %s" %(rgid, picard_clean_sam_cmd))
	Call_Picard_SubProgram(picard_clean_sam_cmd, "CleanSam", logger_std, logger_prog)
	return cleaned_samfile

# call Picard SamFormatConverter to convert sam to bam #
def SAM_to_BAM(cleaned_samfile, prog, rgid, logger_std, logger_prog) :
	cleaned_bamfile = '.'.join(re.split('\.', cleaned_samfile)[0:3]) + ".cleaned.bam"
	picard_samtobam_cmd = "java -Xmx2g -jar %s INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=LENIENT" %(os.path.join(prog, "SamFormatConverter.jar"), cleaned_samfile, cleaned_bamfile)
	logger_std.info("[Run Picard] - [%s INPUT]: %s" %(rgid, cleaned_samfile))
	logger_std.info("[Run Picard] - [%s SamFormatConverter]: %s " %(rgid, picard_samtobam_cmd))
	Call_Picard_SubProgram(picard_samtobam_cmd, "SamFormatConverter", logger_std, logger_prog)

	logger_std.info("[Run Bash] - [%s rm]: %s " %(rgid, "rm %s" %(cleaned_samfile)))
	proc = Popen(["rm", "%s" %(cleaned_samfile)], stderr=PIPE, stdout=PIPE)
	proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 or re.search("(Error|error)", proc_stderr) or re.search("(Error|error)", proc_stdout) :
		logger_std.error("[Pipeline Checker]: Failure to remove %s" %(cleaned_samfile))
		sys.exit(1)
	return cleaned_bamfile

# call Picard MergeSamFiles to merge bam files #
def Merge_Bams(list_bamfiles, outfile, prog, logger_std, logger_prog) :
	input_samfiles = ""
	for each_file in list_bamfiles :
		logger_std.info("[Merging INPUT]: %s" %(each_file))
		if input_samfiles == "" :
			input_samfiles += "INPUT=%s" %(each_file)
		else :
			input_samfiles += " INPUT=%s" %(each_file)
	picard_mergesam_cmd = "java -Xmx2g -jar %s %s OUTPUT=%s SORT_ORDER=coordinate USE_THREADING=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true" %(os.path.join(prog, "MergeSamFiles.jar"), input_samfiles, outfile)
	logger_std.info("[Run Picard] - [MergeSamFiles]: %s" %(picard_mergesam_cmd))
	Call_Picard_SubProgram(picard_mergesam_cmd, "MergeSamFiles", logger_std, logger_prog)

# call Picard AddOrReplaceReadGroups to put in @RG tag #
def AddRGs(bamfile, rg_def, outfile, prog, rgid, logger_std, logger_prog) :
	picard_addRG_cmd = "java -Xmx2g -jar %s INPUT=%s OUTPUT=%s %s SORT_ORDER=coordinate CREATE_INDEX=true" %(os.path.join(prog, "AddOrReplaceReadGroups.jar"), bamfile, outfile, re.sub('\:', '=', ' '.join(re.split(r'\\t', rg_def)[1:])))
	logger_std.info("[Run Picard] - [%s INPUT]: %s" %(rgid, bamfile))
	logger_std.info("[Run Picard] - [%s AddOrReplaceReadGroups]: %s" %(rgid, picard_addRG_cmd))
	Call_Picard_SubProgram(picard_addRG_cmd, "AddOrReplaceReadGroups", logger_std, logger_prog)

# get rid of the duplicates #
def Mark_Duplicates(cleaned_bamfile, prog, rgid, logger_std, logger_prog) :
	sorted_bamfile = '.'.join(re.split('\.', cleaned_bamfile)[0:3]) + ".cleaned.sort.bam"
	picard_sort_cmd = "java -Xmx2g -jar %s I=%s O=%s SO=coordinate VALIDATION_STRINGENCY=LENIENT" %(os.path.join(prog, "SortSam.jar"), cleaned_bamfile, sorted_bamfile)
	logger_std.info("[Run Picard] - [%s INPUT]: %s" %(rgid, cleaned_bamfile))
	logger_std.info("[Run Picard] - [%s SortSam]: %s" %(rgid, picard_sort_cmd))
	Call_Picard_SubProgram(picard_sort_cmd, "SortSam", logger_std, logger_prog)
	dup_log_file = '.'.join(re.split('\.', cleaned_bamfile)[0:3]) + ".dup.report"
	rmdup_bamfile = '.'.join(re.split('\.', sorted_bamfile)[0:3]) + ".cleaned.sorted.rmdup.bam"
	picard_markdup_cmd = "java -Xmx2g -jar %s REMOVE_DUPLICATES=true I=%s O=%s VALIDATION_STRINGENCY=LENIENT M=%s" %(os.path.join(prog, "MarkDuplicates.jar"), sorted_bamfile, rmdup_bamfile, dup_log_file)
	logger_std.info("[Run Picard] - [%s INPUT]: %s" %(rgid, sorted_bamfile))
	logger_std.info("[Run Picard] - [%s MarkDuplicates]: %s" %(rgid, picard_markdup_cmd))
	Call_Picard_SubProgram(picard_markdup_cmd, "MarkDuplicates", logger_std, logger_prog)
	return rmdup_bamfile

def Call_Picard_SubProgram(cmd, subprogram, logger_std, logger_prog) :
	proc = Popen(shlex.split(cmd), stderr=PIPE, stdout=PIPE)
	proc_stdout, proc_stderr = proc.communicate()
	if proc.returncode != 0 or re.search("(Error|error)", proc_stderr) or re.search("(Error|error)", proc_stdout) :
		logger_std.error("[Pipeline Checker]: Picard %s terminates unexpectedly" %(subprogram))
		sys.exit(1)
	else :
		if proc_stderr != "" :
			logger_prog.info("%s\n" %(proc_stderr))

# check whether or not file and directory exist #
def Path_Checker(path, type, logger_std) :
	if type == 'f' :
		tmp_path = re.split(' ', path)
		for i in range(len(tmp_path)) :
			if not os.path.exists(tmp_path[i]) :
				logger_std.error("[Path Checker]: cannnot reach the file: %s" %(tmp_path[i]))
				sys.exit(1)
			else :
				logger_std.debug("[Path Checker]: %s ... [check]" %(tmp_path[i]))
	elif type == 'd' :
		if not os.path.exists(path) :
			os.makedirs(path)
			logger_std.debug("[Path Checker]: Creating directory: %s\n" %(path))

# check RG definition #
def RG_Def_Checker(rg_def, logger_std) :
	tmp_def = re.split(r'\\t', rg_def)
	num_component = 0
	for i in range(1, len(tmp_def)) :
		if len(re.split(':', tmp_def[i])) == 2 :
			tmp_component = re.split(':', tmp_def[i])
			if tmp_component[0] == "" or not tmp_component[0] in ["ID", "SM", "LB", "PL", "PU"] :
				logger_std.error("[@RG Definition Checker]: missing @RG key %s") %(tmp_def[i])
				logger_std.error("[@RG Definition Checker]: @RG\tID:cs1\tSM:cs1\tLB:FRAG\tPL:Illumina\tPU:RUN")
				sys.exit(1)
			else :
				if tmp_component[1] == "" :
					logger_std.error("[@RG Definition Checker]: missing value for @RG_%s %s") %(tmp_component[0], tmp_def[i])
					sys.exit(1)
				else :
					num_component += 1
		else :
			logger_std.error("[@RG Definition Checker]: incorrect @RG format %s") %(tmp_def[i])
			logger_std.error("[@RG Definition Checker]: @RG\tID:cs1\tSM:cs1\tLB:FRAG\tPL:Illumina\tPU:RUN")
			sys.exit(1)
	if num_component != 5 :
		logger_std.error("[@RG Definition Checker]: missing @RG keys in the definition %s" %(rg_def))
		sys.exit(1)

# run Picard pipeline #
def Run_PICARD_Pipe(config_dict, merge_off, outprefix, prog, logger_std, logger_prog) :
	bamfiles_per_RGID = []
	bamfiles_all = []
	for rg_def, samfiles in config_dict.iteritems() :
		rg_id = ""
		tmp_def = re.split(r'\\t', rg_def)
		for i in range(1, len(tmp_def)) :
			tmp_component = re.split(':', tmp_def[i])
			if tmp_component[0] == "ID" :
				rgid = tmp_component[1]
		logger_std.info("[Run Picard] - [%s BEGIN]" %(rgid))
		for each_raw_samfile in re.split('\s+', samfiles) :
			each_cleaned_samfile = Clean_SAM_Files(each_raw_samfile, prog, rgid, logger_std, logger_prog)
			each_cleaned_bamfile = SAM_to_BAM(each_cleaned_samfile, prog, rgid, logger_std, logger_prog)
			each_rmdup_bamfile = Mark_Duplicates(each_cleaned_bamfile, prog, rgid, logger_std, logger_prog)
			bamfiles_per_RGID.append(each_rmdup_bamfile)
		processed_bamfile = os.path.join(each_rmdup_bamfile.strip(re.split('\/', each_rmdup_bamfile)[-1]), rgid+".processed.bam")
		if len(bamfiles_per_RGID) > 1 :
			sample_bamfile = os.path.join(each_rmdup_bamfile.strip(re.split('\/', each_rmdup_bamfile)[-1]), rgid+".merged.bam")
			logger_std.info("[Merging BAM files of the same @RG]")
			#for bamfile in bamfiles_per_RGID :
			#	logger_std.info("[Run Picard] - [%s INPUT]: %s" %(rgid, bamfile))
			Merge_Bams(bamfiles_per_RGID, sample_bamfile, prog, logger_std, logger_prog)
			AddRGs(sample_bamfile, rg_def, processed_bamfile, prog, rgid, logger_std, logger_prog)
			bamfiles_all.append(processed_bamfile)
		else :
			#logger_std.info("[Run Picard] - [%s INPUT]: %s" %(rgid, each_rmdup_bamfile))
			AddRGs(each_rmdup_bamfile, rg_def, processed_bamfile, prog, rgid, logger_std, logger_prog)
			bamfiles_all.append(processed_bamfile)
		logger_std.info("[RUn Picard] - [%s Finish]\n" %(rgid))
		bamfiles_per_RGID = []
	if len(bamfiles_all) > 1 :
		if not merge_off :
			logger_std.info("[Merging all BAM files of different @RG]")
			#for bamfile in bamfiles_all :
			Merge_Bams(bamfiles_all, outprefix+".all.processed.bam", prog, logger_std, logger_prog)

# parse config file and read in each to-be-processed SAM file #
def Config_Parser(config_file, logger_std) :
	line_id, rg_id, rg_lb, rg_pl, rg_pu, rg_sm = "", "", "", "", "", ""			# line id is defined by users, and should be consistent with the each sam filename #
	num_diff_RGs = 0				# denote how many RG tags defined in the file #
	raw_samfiles = []
	config_dict = {}
	fCon = open(config_file, 'r')
	for line in fCon :
		if not line.startswith('#') :
			tmp_line = re.split('\t', line.strip())
			rg_def = tmp_line[0]
			RG_Def_Checker(rg_def, logger_std)
			if len(tmp_line) > 2 :
				raw_samfiles = ' '.join(tmp_line[1:])
			else :
				raw_samfiles = tmp_line[1]
			Path_Checker(raw_samfiles, 'f', logger_std)
			if not config_dict.has_key(rg_def) :
				config_dict[rg_def] = raw_samfiles
			else :
				config_dict[rg_def] += ' ' + raw_samfiles
				logger_std.warning("[Config File Checker]: Repeated @RG definition %s in the config file" %(rg_def))
				for i in range(len(re.split(' ', config_dict[rg_def]))) :
					logger_std.warning("[Config File Checker]: %s" %(re.split(' ', config_dict[rg_def])[i]))
	fCon.close()
	return config_dict

# get SAM files from -sam option #
def Get_SAMs_from_CMD(samfiles, logger_std) :
	'''
	fix this function
	'''
	contig_dict = {}

	return contig_dict

# initiate loggers #
def Logger_Init(picard_log) :
	# logger_std records everything coming out from this script #
	logger_std = logging.getLogger("post_mapping_pipeline.py")
	logger_std.setLevel(logging.DEBUG)
	formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %I:%m:%S %p")

	'''
	normal events and warnings go to STDOUT
	error messages go to STDERR
	'''

	# add a info-level stream handler #
	stdout_handler = logging.StreamHandler(sys.stdout)
	stdout_handler.setFormatter(formatter)
	stdout_handler.addFilter(Logging_Level_Filter(logging.INFO))
	logger_std.addHandler(stdout_handler)

	# add a warning-level stream handler #
	warning_handler = logging.StreamHandler(sys.stdout)
	warning_handler.setFormatter(formatter)
	warning_handler.addFilter(Logging_Level_Filter(logging.WARNING))
	logger_std.addHandler(warning_handler)

	# add a error-level stream handler #
	stderr_handler = logging.StreamHandler(sys.stderr)
	stderr_handler.setFormatter(formatter)
	stderr_handler.addFilter(Logging_Level_Filter(logging.ERROR))
	logger_std.addHandler(stderr_handler)

	# logger_prog to record everything comming out from Picard #
	logger_prog = logging.getLogger("Picard")
	logger_prog.setLevel(logging.DEBUG)

	# add a info-level file handler #
	prog_stderr_handler = logging.FileHandler(picard_log, mode='a')
	prog_stderr_handler.setFormatter(formatter)
	prog_stderr_handler.addFilter(Logging_Level_Filter(logging.INFO))
	logger_prog.addHandler(prog_stderr_handler)

	return logger_std, logger_prog

if __name__ == "__main__" :

	# parse the options and arguments #
	parser = argparse.ArgumentParser(description="Calling PICARD to post-process mapping results in sam files. Java 1.6+ by default is assumed to be installed")
	parser.add_argument("-config", metavar="FILE", dest="config_file", help="a config file providing paths to sam files and @RG tag information. Note: either -config or -s can be specified at a time. In cases where they are specified together, only sam files in the config file will be processed")
	parser.add_argument("-sam", metavar="FILE", nargs='*', dest="samfiles", help="a list of sam files separated by space. Note: either -config or -s can be specified at a time. [NOT READY]")
	parser.add_argument("-RG", metavar="STR", nargs='*', dest="rgs", help="@RG definitions, e.g. @RG\tID:foo\tSM:foo\tLB:FRAG\tPL:Illumina\tPU:RUN. Note: if multiple sam files provided, then the same number of @RG definitions provided. If only one @RG defition provided, then the @RG will be added in each of the sam files. [NOT READY]")
	parser.add_argument("-picard", metavar="DIR", dest="prog", help="the full path to the directory where PICARD is installed. Note: if the path to PICARD is not under your system path, please specify")
	parser.add_argument("-mergeOFF", action='store_true', dest="merge_off", help="toggle off the mergeSamFiles function")
	parser.add_argument("-p", metavar="FILE", dest="out_prefix",  help="the prefix of the output file. Note: if -mergeOFF is provided, no output will be created no matter what you provide after -p")
	parser.add_argument("-picard_log", metavar="FILE", dest="picard_log", required=True, help="log generated from Picard STDERR")

	args = parser.parse_args()

	# initiate logger_std #
	if os.path.exists(args.picard_log) :
		os.remove(args.picard_log)
	logger_std, logger_prog = Logger_Init(args.picard_log)

	out_prefix, outdir = "", ""
	if not args.merge_off :
		if not args.out_prefix :
			logger_std.error("[Command Line Checker]: no prefix for output file provided\n")
			sys.exit(1)
		else :
			if re.search('\/', args.out_prefix) :
				outdir = args.out_prefix.strip(re.split('\/', args.out_prefix)[-1])
				Path_Checker(outdir, 'd', logger_std)
				out_prefix = args.out_prefix
			else :
				outdir = os.getcwd()
				out_prefix = os.path.join(outdir, args.out_prefix)

	# check path to Picard #
	if not os.path.exists(args.prog) :
		logger_std.error("[Command Line Checker]: Path to Picard -%s- does not exist" %(args.prog))
		sys.exit(1)

	contig_dict = {}
	if args.config_file and (args.samfiles or not args.samfiles) :
		if args.samfiles :
			logger_std.warning("[Command Line Checker]: all SAM file provided after -s will be ignored")
			args.samfiles = []
		Path_Checker(args.config_file, 'f', logger_std)
		config_dict = Config_Parser(args.config_file, logger_std)
		Run_PICARD_Pipe(config_dict, args.merge_off, out_prefix, args.prog, logger_std, logger_prog)
	else :
		if not args.rgs :
			logger_std.error("[Command Line Checker]: no @RG definition found for your input SAM files\n")
			sys.exit(1)
		else :
			if len(args.rgs) > 1 :
				if len(args.rgs) != len(args.samfiles, args.rgs) :
					logger_std.error("[Command Line Checker]: number of @RG definition is more than the number of SAM files provided\n")
					sys.exit(1)
				else :
					config_dict = Get_SAMs_from_CMD(args.samfiles, logger_std)
			else :
				if len(args.rgs) < len(args.samfiles) and len(args.rgs) == 1 :
					logger_std.warning("[Command Line Checker]: one @RG definition for all %d SAM files\n" %(len(args.samfiles)))
					RG_Def_Checker(args.rgs[0], logger_std)
					contig_dict[args.rgs[0]] = args.samfiles
