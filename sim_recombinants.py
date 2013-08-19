# The script is used to simulate recombinants for a given diploid genome #
# Currently, it can only support one cross between maternal and paternal chromosome #
# The program works in the following way:
#	1. for each pair, choose either maternal or paternal chromosome (based on a uniform distribution) on which a recombination break point will be created #
#	2. a break point will be created uniformly and interchange chromosomal segments between the maternal and paternal chromosome #
#	3. repeat step 1 and 2 until all the chromosome pairs being processed #

import os, sys
import glob
import re
import random
import string
import math

def read_genome(haploid_genome_file, which_copy) :
	haploid_handle = open(haploid_genome_file, 'r')
	chr_id, chr_seq = "", ""
	num_chr = 0
	haploid_chr_dict = {}
	for line in haploid_handle :
		if line.startswith('>') :
			if chr_seq != "" :
				haploid_chr_dict[chr_id] = chr_seq
				num_chr += 1
				chr_seq = ""
			chr_id = re.split(" ", line.strip())[0][1:]
			header = ">" + chr_id + "|%s" %(which_copy)
		else :
			chr_seq += line.strip()
	# last chromosome sequence #
	if chr_seq != "" :
		haploid_chr_dict[chr_id] = chr_seq
		num_chr += 1
	haploid_handle.close()

	return haploid_chr_dict, num_chr

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
	print breakpoints
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

	# Checking if the recombined sequences are in the places as expected #
	if new_chr_seq[breakpoints[0]-100:breakpoints[0]] == chr_seq_1[breakpoints[0]-100:breakpoints[0]].upper() :
		print "sequence to the first breakpoint is right"
	else :
		print "sequence to the first breakpoint is not right"
		print new_chr_seq[breakpoints[0]-5:breakpoints[0]]
		print chr_seq_1[breakpoints[0]-5:breakpoints[0]]
		sys.exit()
	if new_chr_seq[breakpoints[0]:breakpoints[1]] == chr_seq_2[breakpoints[0]:breakpoints[1]].upper() :
		print "sequence from the frist to the second breakpoint is right"
	else :
		print "sequence from the first to the second breakpoint is not right"
		print new_chr_seq[breakpoints[1]-5:breakpoints[1]]
		print chr_seq_2[breakpoints[1]-5:breakpoints[1]]
		sys.exit()
	if new_chr_seq[breakpoints[1]+1:len(new_chr_seq)] == chr_seq_1[breakpoints[1]+1:len(new_chr_seq)].upper() :
		print "sequence from the second breakpoint to the end of the sequence is right"
	else :
		print "sequence from the second breakpoint to the end of the sequence is right"
		#print new_chr_seq[breakpoints[1]+1:len(new_chr_seq)]
		#print chr_seq_1[breakpoints[1]+1:len(new_chr_seq)]
	return new_chr_seq, breakpoints

# get the breakpoints #
def get_breakpoints(num_recomb_per_chr, len_chr_seq, chr_seq) :
	breakpoints = []
	len_to_ends, len_in_between = 500, 1000								# length from breakpoints to the ends, and length in between the breakpoints #
	which_pos, pre_which_pos = 0, 0
	for j in range(num_recomb_per_chr) :
		while True :
			which_pos = random.randint(0, len_chr_seq - 1)				# the position of the breakpoint #
			#if len_chr_seq - which_pos >= len_to_ends and which_pos >= len_to_ends and not (re.search('N+', chr_seq[which_pos:which_pos+10]) or re.search('N+', chr_seq[which_pos-9:which_pos+1])) :
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

def sim_recombinants(cross_design, num_recomb_per_chr, mchr_dict, pchr_dict, num_recombinants, out_prefix) :
	# simulate recombinants #
	for i in range(num_recombinants) :
		out_m_haploid = out_prefix+"%s_m.fasta" %(i+1)
		out_p_haploid = out_prefix+"%s_p.fasta" %(i+1)
		out_m_haploid_handle = open(out_m_haploid, 'w')
		out_p_haploid_handle = open(out_p_haploid, 'w')

		for chr_id, mchr_seq in mchr_dict.iteritems() :
			# 1. get the corresopnding pair of chromosome, e.g. chr01_m and chr01_p #
			pchr_seq = pchr_dict[chr_id]
			# 2. for each of the pair, we need to simulate the meiosis to generate maternal and paternal gametes, from which to obtain a recombinant for the current pair #
			m_or_p = random.randint(1,2)					# randomly choose to either maternal (1) or paternal (2) copy of the current chromosome for generate one of the two gametes #
			if m_or_p == 1 :
				print "Making one of the two gametes from the maternal copy of the current chromosome %s" %(chr_id)
				one_gamete, breakpoints_gamete = make_gametes(mchr_seq, pchr_seq, num_recomb_per_chr)			# breakpoints will be made on the maternal copy of the current chromosome #
			else :
				print "Making one of the two gametes from the paternal copy of the current chromosome %s" %(chr_id)
				one_gamete, breakpoints_gamete = make_gametes(pchr_seq, mchr_seq, num_recomb_per_chr)			# breakpoints will be made on the paternal copy of the current chromosome #
			sys.stdout.write("Output the first gamete of the current chromosome to file ...")
			outputter(one_gamete, m_or_p, chr_id, breakpoints_gamete, out_m_haploid_handle)
			print " [Done]"
			if cross_design == "RILs" :
				sys.stdout.write("Output the first gamete of the current chromosome to file ...")
				outputter(one_gamete, m_or_p, chr_id, breakpoints_gamete, out_p_haploid_handle)
				print " [Done]"
			else :
				m_or_p = random.randint(1,2)					# randomly choose to either maternal (1) or paternal (2) copy of the current chromosome for generate the second gamete #
				if m_or_p == 1 :
					print "Making the second gametes from the maternal copy of the current chromosome %s" %(chr_id)
					second_gamete, breakpoints_gamete = make_gametes(mchr_seq, pchr_seq, num_recomb_per_chr)
				else :
					print "Making the second gametes from the paternal copy of the current chromosome %s" %(chr_id)
					second_gamete, breakpoints_gamete = make_gametes(pchr_seq, mchr_seq, num_recomb_per_chr)
				sys.stdout.write("Output the second gamete of the current chromosome to file ...")
				outputter(second_gamete, m_or_p, chr_id, breakpoints_gamete, out_p_haploid_handle)
				print " [Done]"

		out_p_haploid_handle.close()
		out_m_haploid_handle.close()
		print "recombinant %d simulation is finished" %(i+1)

if __name__ == "__main__" :
	cross_design = sys.argv[1]
	num_recomb_per_chr = int(sys.argv[2])			# this is the number of crossovers per chromosome of the last generation #
	mhaploid = sys.argv[3]
	phaploid = sys.argv[4]
	num_recombinants = int(sys.argv[5])
	out_prefix = sys.argv[6]

	# setup the output path #
	if re.search("\/", out_prefix) :
		tmp_dir = out_prefix.rstrip(os.path.basename(out_prefix))
		if not os.path.exists(tmp_dir) :
			os.makedirs(tmp_dir)
	else :
		# the output file with the defined prefix name will be created under the current working directory, if no other output path being specified #
		print os.getcwd()
		out_prefix = os.path.join(os.getcwd(), out_prefix)

	num_mchr = 0
	sys.stdout.write("reading maternal copy ... ")
	mchr_dict, num_mchr = read_genome(mhaploid, "maternal")
	sys.stdout.write("%d chromosomes are included in the maternal copy of the genome" %(num_mchr))
	print " ... [Done]"

	num_pchr = 0
	sys.stdout.write("reading paternal copy ... ")
	pchr_dict, num_pchr = read_genome(phaploid, "paternal")
	sys.stdout.write("%d chromosomes are included in the paternal copy of the genome" %(num_pchr))
	print " ... [Done]"

	if num_pchr != num_mchr :
		print "Error: you have different number of chromosomes between your maternal (num_mchr) and paternal (num_pchr) copies of the genome"
		sys.exit()

	sim_recombinants(cross_design, num_recomb_per_chr, mchr_dict, pchr_dict, num_recombinants, out_prefix)

