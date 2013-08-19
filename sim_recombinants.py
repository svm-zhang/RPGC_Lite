# The script is used to simulate recombinants for a given diploid genome #
# Currently, it can only support one cross between maternal and paternal chromosome #
# The program works in the following way:
#	1. for each pair, choose either maternal or paternal chromosome (based on a uniform distribution) on which a recombination break point will be created #
#	2. a break point will be created uniformly and interchange chromosomal segments between the maternal and paternal chromosome #
#	3. repeat step 1 and 2 until all the chromosome pairs being processed #

import re
import random
import math
from argparse import ArgumentParser
from os import makedirs
from os.path import exists, join, dirname, realpath
from sys import exit, stdout, stderr

def read_genome(seq_file, which_parent) :
	"""
	read each of the parental genomes
	"""
	fSEQ = open(seq_file, 'r')
	chr_id, chr_seq = "", ""
	nchr = 0
	chr_dict = {}
	for line in fSEQ :
		if line.startswith('>') :
			if chr_seq != "" :
				chr_dict[chr_id] = chr_seq
				nchr += 1
				chr_seq = ""
			chr_id = re.split(" ", line.strip())[0][1:]
			header = ">" + chr_id + "|%s" %(which_parent)
		else :
			chr_seq += line.strip()
	# last chromosome sequence #
	if chr_seq != "" :
		chr_dict[chr_id] = chr_seq
		nchr += 1
	fSEQ.close()

	return chr_dict, nchr

def outputter(new_chr_seq, m_or_p, chr_id, breakpoints, out_handle) :
	tmp_breakpoints = ""
	for breakpoint in breakpoints :
		if tmp_breakpoints == "" :
			tmp_breakpoints += str(breakpoint)
		else :
			tmp_breakpoints += "_%s" %(str(breakpoint))
	if m_or_p == 1 :
		out_handle.write(">%s|maternal_%s\n" %(chr_id, tmp_breakpoints))
		out_handle.write("%s\n" %(make_fasta_format(new_chr_seq)))
	else :
		out_handle.write(">%s|paternal_%s\n" %(chr_id, tmp_breakpoints))
		out_handle.write("%s\n" %(make_fasta_format(new_chr_seq)))
	out_handle.flush()

def make_fasta_format(seq) :
	i = 0
	len_seq = len(seq)
	formatted_seq = ""
	for i in range(len_seq) :
		if i % 100 == 0 and i != 0 :
			formatted_seq += "\n"
		formatted_seq += seq[i]
	return formatted_seq

def make_gametes(chr_seq_1, chr_seq_2, num_recomb_per_chr) :
	# chr_seq_1 could be either maternal or paternal copy, so does chr_seq_2. It depends on how the function being called #
	# breakpoints will be always made on the chr_seq_1 #
	breakpoints = get_breakpoints(num_recomb_per_chr, len(chr_seq_1), chr_seq_1)
	new_chr_seq = ""
	i = 0
	while i <= len(breakpoints) :
		if i%2 == 0 :
			if i == 0 :
				new_chr_seq += chr_seq_1[0:breakpoints[i]]
			elif i == len(breakpoints) :
				new_chr_seq += chr_seq_1[breakpoints[i-1] : len(chr_seq_1)]
			else :
				new_chr_seq += chr_seq_1[breakpoints[i-1] : breakpoints[i]]
		else :
			if i == len(breakpoints) :
				new_chr_seq += chr_seq_2[breakpoints[i-1] : len(chr_seq_2)]
			else :
				new_chr_seq += chr_seq_2[breakpoints[i-1] : breakpoints[i]]
		i += 1
	new_chr_seq = new_chr_seq.upper()

	"""
		DEBUG purpose
	"""
	# Checking if the recombined sequences are in the places as expected #
	#if new_chr_seq[breakpoints[0]-100:breakpoints[0]] == chr_seq_1[breakpoints[0]-100:breakpoints[0]].upper() :
	#	print "sequence to the first breakpoint is right"
	#else :
	#	print "sequence to the first breakpoint is not right"
	#	print new_chr_seq[breakpoints[0]-5:breakpoints[0]]
	#	print chr_seq_1[breakpoints[0]-5:breakpoints[0]]
	#	sys.exit()
	#if new_chr_seq[breakpoints[0]:breakpoints[1]] == chr_seq_2[breakpoints[0]:breakpoints[1]].upper() :
	#	print "sequence from the frist to the second breakpoint is right"
	#else :
	#	print "sequence from the first to the second breakpoint is not right"
	#	print new_chr_seq[breakpoints[1]-5:breakpoints[1]]
	#	print chr_seq_2[breakpoints[1]-5:breakpoints[1]]
	#	sys.exit()
	#if new_chr_seq[breakpoints[1]+1:len(new_chr_seq)] == chr_seq_1[breakpoints[1]+1:len(new_chr_seq)].upper() :
	#	print "sequence from the second breakpoint to the end of the sequence is right"
	#else :
	#	print "sequence from the second breakpoint to the end of the sequence is right"
	#	print new_chr_seq[breakpoints[1]+1:len(new_chr_seq)]
	#	print chr_seq_1[breakpoints[1]+1:len(new_chr_seq)]
	return new_chr_seq, breakpoints

# get the breakpoints #
def get_breakpoints(num_recomb_per_chr, len_chr_seq, chr_seq) :
	breakpoints = []
	len_to_ends, len_in_between = 500, 1000								# length from breakpoints to the ends, and length in between the breakpoints #
	which_pos, pre_which_pos = 0, 0
	for j in range(num_recomb_per_chr) :
		while True :
			which_pos = random.randint(0, len_chr_seq - 1)				# the position of the breakpoint #
			if len_chr_seq - which_pos >= len_to_ends and which_pos >= len_to_ends and chr_seq[which_pos] != 'N' :				# the breakpoints cannot be inside of a block of Ns #
				break
		if pre_which_pos == 0 :
			breakpoints.append(which_pos)
			pre_which_pos = which_pos
		else :
			if math.fabs(which_pos - pre_which_pos) >= len_in_between :
				breakpoints.append(which_pos)
	breakpoints.sort()
	return breakpoints

def sim_recombinants(cross_design, num_recomb_per_chr, mchr_dict, pchr_dict, popsize, outprefix) :
	# simulate recombinants #
	for i in range(popsize) :
		mHap_outfile = outprefix+"%s_m.fasta" %(i+1)
		pHap_outfile = outprefix+"%s_p.fasta" %(i+1)
		mHapOut = open(mHap_outfile, 'w')
		pHapOut = open(pHap_outfile, 'w')

		for chr_id, mchr_seq in mchr_dict.iteritems() :
			# 1. get the corresopnding pair of chromosome, e.g. chr01_m and chr01_p #
			pchr_seq = pchr_dict[chr_id]
			# 2. for each of the pair, we need to simulate the meiosis to generate maternal and paternal gametes, from which to obtain a recombinant for the current pair #
			m_or_p = random.randint(1,2)					# randomly choose to either maternal (1) or paternal (2) copy of the current chromosome for generate one of the two gametes #
			if m_or_p == 1 :
				one_gamete, breakpoints_gamete = make_gametes(mchr_seq, pchr_seq, num_recomb_per_chr)			# breakpoints will be made on the maternal copy of the current chromosome #
			else :
				one_gamete, breakpoints_gamete = make_gametes(pchr_seq, mchr_seq, num_recomb_per_chr)			# breakpoints will be made on the paternal copy of the current chromosome #
			outputter(one_gamete, m_or_p, chr_id, breakpoints_gamete, mHapOut)
			if cross_design == "RILs" :
				outputter(one_gamete, m_or_p, chr_id, breakpoints_gamete, pHapOut)
			else :
				m_or_p = random.randint(1,2)					# randomly choose to either maternal (1) or paternal (2) copy of the current chromosome for generate the second gamete #
				if m_or_p == 1 :
					second_gamete, breakpoints_gamete = make_gametes(mchr_seq, pchr_seq, num_recomb_per_chr)
				else :
					second_gamete, breakpoints_gamete = make_gametes(pchr_seq, mchr_seq, num_recomb_per_chr)
				outputter(second_gamete, m_or_p, chr_id, breakpoints_gamete, pHapOut)

		pHapOut.close()
		mHapOut.close()
		stdout.write("recombinant %d is simulated\n" %(i+1))

if __name__ == "__main__" :
	parser = ArgumentParser(description="Simulating a population of recombinants with user-defined number of cross-over events per chromosome, e.g. F2s, RILs")
	parser.add_argument("-cross_design", metavar="STR", dest="cross_design", required=True, choices=["RILs", "F2"], help="Specify the type of crossing design: F2s, RILs")
	parser.add_argument("-p1", metavar="FILE", dest="p1", required=True, help="Specify sequence file (in fasta) of one of the parents")
	parser.add_argument("-p2", metavar="FILE", dest="p2", required=True, help="Specify sequence file (in fasta) of the other parent")
	parser.add_argument("-popsize", metavar="INT", dest="popsize", type=int, default=96, help="Specify the number of recombinants you'd like in your population. DEFAULT: 96")
	parser.add_argument("-ncross", metavar="INT", dest="num_recomb_per_chr", type=int, default=5, help="Specify the number of crossing-over events per chrmosome. DEFAULT: 5")
	parser.add_argument("-outprefix", metavar="STR", dest="outprefix", help="Specify the output prefix")

	args = parser.parse_args()

	# setup the output path #
	abs_outprefix = realpath(args.outprefix)
	print abs_outprefix
	outdir = dirname(abs_outprefix)
	if not exists(outdir) :
		makedirs(outdir)

	nchr_p1 = 0
	stdout.write("reading one parent: %s\n" %(args.p1))
	mchr_dict, nchr_p1 = read_genome(args.p1, "p1")
	stdout.write("\t%d chromosomes parsed\n" %(nchr_p1))

	nchr_p2 = 0
	stdout.write("reading the other parent: %s\n" %(args.p2))
	pchr_dict, nchr_p2 = read_genome(args.p2, "p2")
	stdout.write("\t%d chromosomes parsed\n" %(nchr_p2))

	if nchr_p1 != nchr_p2 :
		stderr.write("Error: two parents contain different numbers of chromosomes")
		exit()

	sim_recombinants(args.cross_design, args.num_recomb_per_chr, mchr_dict, pchr_dict, args.popsize, abs_outprefix)
