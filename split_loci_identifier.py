# This script will push all the progress report to stderr, and results to stdout

import re, shlex
from argparse import ArgumentParser
from os import makedirs, stat, system, devnull
from os.path import join, exists, dirname, realpath
from sys import exit, stdout, stderr, argv
from datetime import datetime
from subprocess import Popen

def identify_splited_loci(alignment_file, assembled_genome_file, all_sample_bam_file, out_dir, splitted_sample_bam_dir, variants_dict, samples, smt_genotype_dict, assembled_genome_dict, min_aln_len, min_idn) :
	"""
	get a set of pairs of loci to be candidates of falsely splitted loci, given an alignment of an assembly against itself
	"""

	stderr.write(timestamper() + " Identifying splitted loci candidates\n")
	"""
	setting up sub-directories for pileups
	1. pileups around markers
	2. pileups of a pair of splitted loci candidate
	"""
	pileups_dir = join(out_dir, "pileups")
	markers_wide_pileups_dir = join(pileups_dir, "markers_wide")
	loci_wide_pileups_dir = join(pileups_dir, "loci_wide")
	make_dirs_if_needed(markers_wide_pileups_dir, loci_wide_pileups_dir)

	# setting up sub-directory for sequences #
	# sequences of pairs of splitted loci candidates with(out) markers #
	seqs_dir = join(out_dir, "seqs")
	make_dirs_if_needed(seqs_dir)
	splitted_loci_with_markers_seq_outfile = join(seqs_dir, "splitted_loci_with_markers.fasta")
	splitted_loci_without_markers_seq_outfile = join(seqs_dir, "splitted_loci_without_markers.fasta")
	fSEQ_LOCI_MARKER = open(splitted_loci_with_markers_seq_outfile, 'w')
	fSEQ_LOCI_NO_MARKER = open(splitted_loci_without_markers_seq_outfile, 'w')

	# setting up subdirectory for summary #
	report_dir = join(out_dir, "summary")
	make_dirs_if_needed(report_dir)
	proposed_genotypes_outfile = join(report_dir, "non_tandem_splitted_loci.proposed_genotypes")
	fProPose = open(proposed_genotypes_outfile, 'w')

	num_potentials = 0
	potential_dup_loci = []					# a list of pairs of loci that are potentially falsely splitted loci #
	all_splitted_loci_dict = {}
	cigar_dict = {}							# a dictionary to record the alignment detail (CIGAR string) between each pair of loci #
	variants_in_loci_candidates_db = {}
	num_localized_variants, num_loci_with_no_localized_variants = 0, 0
	num_pairs_predicted_falsely_splitted_loci, num_pairs_predicted_falsely_non_tandem_splitted_loci = 0, 0
	num_pairs_loci_predicted_correctly_splitted_loci = 0
	loci_without_markers_localized = []
	predicted_falsely_splitted_loci = []
	num_loci_failed_cov_filter = 0
	total_num_loci, nloci_pass_aln_filter = 0, 0
	fALIGN = open(alignment_file, 'r')
	for line in fALIGN :
		tmp_line = re.split(' ', line.strip())
		t_strand = tmp_line[8]
		q_strand = tmp_line[4]
		t_start = int(tmp_line[6]) + 1
		t_end = int(tmp_line[7]) + 1
		q_start = int(tmp_line[2]) + 1
		q_end = int(tmp_line[3]) + 1
		total_num_loci += 1
		if q_start > t_end or t_start > q_end :
			"""
			get the length and the identity of the alignment using CIGAR
			"""
			cigar = tmp_line[10:]
			i, len_alignment, num_matches = 0, 0, 0
			while i < len(cigar) :
				if cigar[i] in ['M', 'D', 'I'] :
					if cigar[i] == 'M' :
						match_flag = 1
				else :
					if match_flag == 1 :
						num_matches += int(cigar[i])
						match_flag = 0
					len_alignment += int(cigar[i])
				i += 1
			identity_alignment = float(num_matches) / len_alignment
			if len_alignment >= min_aln_len and identity_alignment >= min_idn :		# apply alignment length and identity filters #
				t_id = tmp_line[5]
				q_id = tmp_line[1]
				one_locus = "%s:%d-%d" %(t_id, t_start, t_end)
				if q_strand == '-' :
					the_other_locus = "%s:%d-%d" %(q_id, q_end, q_start)
				else :
					the_other_locus = "%s:%d-%d" %(q_id, q_start, q_end)
				if not all_splitted_loci_dict.has_key(the_other_locus+'|'+one_locus) :
					all_splitted_loci_dict[one_locus+'|'+the_other_locus] = 1
					nloci_pass_aln_filter += 1
					stdout.write(">Candidate pair: %s %s\n" %(one_locus, the_other_locus))

					proposed_genotype_at_variant_site = {}
					tmp_loc_1 = re.sub(':', '-', re.sub('_', '', one_locus))
					out_pileups_loc_1 = join(loci_wide_pileups_dir, tmp_loc_1+".pileups")
					tmp_loc_2 = re.sub(':', '-', re.sub('_', '', the_other_locus))
					out_pileups_loc_2 = join(loci_wide_pileups_dir, tmp_loc_2+".pileups")
					"""
					apply coverage filter
					"""
					pass_or_filter = apply_cov_filter(one_locus, the_other_locus, assembled_genome_file, all_sample_bam_file, out_pileups_loc_1, out_pileups_loc_2)

					if pass_or_filter :
						"""
						localizing markers in each pair of loci
						"""
						variants_loc_1 = find_markers_in_locus(t_id, t_start, t_end, t_strand, variants_dict, smt_genotype_dict)
						variants_loc_2 = find_markers_in_locus(q_id, q_start, q_end, q_strand, variants_dict, smt_genotype_dict)

						"""
						match markers in each pair of loci
						"""
						variant_pairs = []
						if len(variants_loc_1) > 0 :
							one_or_the_other = 0
							variant_pairs += match_markers_in_loci(t_id, t_start, t_end, t_strand, q_id, q_start, q_end, q_strand, one_or_the_other, variants_loc_1, cigar, variants_dict, assembled_genome_dict)
						if len(variants_loc_2) > 0 :
							one_or_the_other = 1
							variant_pairs += match_markers_in_loci(q_id, q_start, q_end, q_strand, t_id, t_start, t_end, t_strand, one_or_the_other, variants_loc_2, cigar, variants_dict, assembled_genome_dict)

						"""
						creating a database of variants that are identified in each pair of loci candidates #
						"""
						for i in range(len(variant_pairs)) :
							tmp_pair = re.split('\|', variant_pairs[i])
							for j in range(len(tmp_pair)) :
								tmp_info = tmp_pair[j]
								tmp_id = re.split(':', tmp_info)[0]
								tmp_pos = int(re.split(':', tmp_info)[1])
								if not variants_in_loci_candidates_db.has_key(tmp_id) :
									variants_in_loci_candidates_db[tmp_id] = [tmp_pos]
								else :
									variants_in_loci_candidates_db[tmp_id].append(tmp_pos)

						num_localized_variants += len(variants_loc_1) + len(variants_loc_2)
						if len(variants_loc_1) > 0 or len(variants_loc_2) > 0 :
							#stdout.write("\t%d pairs of markers matched\n" %(num_localized_variants))
							fSEQ_LOCI_MARKER.write(">%s_%d_%d_1\n%s\n" %(t_id, t_start, t_end, assembled_genome_dict[t_id][t_start-1:t_end]))
							if q_strand == '+' :
								fSEQ_LOCI_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, assembled_genome_dict[q_id][q_start-1:q_end]))
							else :
								fSEQ_LOCI_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, rc(assembled_genome_dict[q_id][q_end-1:q_start])))
							loci_candidate_out_subdir = join(markers_wide_pileups_dir, "%s_%d_%d-%s_%d_%d" %(t_id, t_start, t_end, q_id, q_start, q_end))
							make_dirs_if_needed(loci_candidate_out_subdir)
							proposed_genotype_at_variant_site = investigate_genotypes(variant_pairs, variants_dict, samples, loci_candidate_out_subdir, splitted_sample_bam_dir, assembled_genome_file)
							if proposed_genotype_at_variant_site != None :
								#print "\tproposed genotype:", proposed_genotype_at_variant_site
								predicted_falsely_splitted_loci.append("%s|%s" %(one_locus, the_other_locus))
								num_pairs_predicted_falsely_splitted_loci += 1
								if t_id != q_id :
									num_pairs_predicted_falsely_non_tandem_splitted_loci += 1
									for variant in proposed_genotype_at_variant_site.iterkeys() :
										tmp_variant = re.split(':', variant)
										tmp_id = tmp_variant[0]
										tmp_pos = tmp_variant[1]
										fProPose.write("%s\t%s\t%s\n" %(tmp_id, tmp_pos, '\t'.join(proposed_genotype_at_variant_site[variant])))
							else :
								num_pairs_loci_predicted_correctly_splitted_loci += 1
						else :
							loci_without_markers_localized.append("%s|%s" %(one_locus, the_other_locus))
							fSEQ_LOCI_NO_MARKER.write(">%s_%d_%d_1\n%s\n" %(t_id, t_start, t_end, assembled_genome_dict[t_id][t_start-1:t_end]))
							if q_strand == '+' :
								fSEQ_LOCI_NO_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, assembled_genome_dict[q_id][q_start-1:q_end]))
							else :
								fSEQ_LOCI_NO_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, rc(assembled_genome_dict[q_id][q_end-1:q_start])))
							num_loci_with_no_localized_variants += 1
					else :
						num_loci_failed_cov_filter += 1
						system("rm %s" %(out_pileups_loc_1))
						system("rm %s" %(out_pileups_loc_2))
	fALIGN.close()

	# for each loci without localized markers, we try to use the loci with markers to help #
	# for example, one copy in the reference, but more than two copies in the assembly	   #
	#print "use loci with markers to help others who don't have"
	#for i in range(len(loci_without_markers_localized)) :
	#	print loci_without_markers_localized[i]
	#	tmp_pair = re.split('\|', loci_without_markers_localized[i])
	#	for j in range(len(tmp_pair)) :
	#		tmp_id = re.split(':', tmp_pair[j])[0]
	#		tmp_start = int(re.split('-', re.split(':', tmp_pair[j])[1])[0])
	#		tmp_end = int(re.split('-', re.split(':', tmp_pair[j])[1])[1])
	#		print tmp_id, tmp_start, tmp_end
	#		if variants_in_loci_candidates_db.has_key(tmp_id) :
	#			if variants_in_loci_candidates_db[tmp_id] >= tmp_start and variants_in_loci_candidates_db[tmp_id] <= tmp_end :
	#				print "found", tmp_pair[j], variants_in_loci_candidates_db[tmp_id]

	stdout.write("\n" + timestamper() + " Predicted splitted loci:\n"), predicted_falsely_splitted_loci
	for loci in predicted_falsely_splitted_loci :
		stdout.write(timestamper() + " %s\n" %(loci))
	report_file = join(report_dir, "splitted_loci_identification.summary")
	fReport = open(report_file, 'w')
	fReport.write("total number of loci in the alignment: %d\n" %(total_num_loci))
	fReport.write("number of loci alignment having at least 1000 bp long and 90%% identity: %d\n" %(nloci_pass_aln_filter))
	fReport.write("number of loci alignment who passed the length and identity filter failing to pass the coverage filter: %d\n" %(num_loci_failed_cov_filter))
	fReport.write("number of loci alignment who passed all three filters having no localized variants: %d\n" %(num_loci_with_no_localized_variants))
	fReport.write("number of pairs of loci being identified as falsely splitted loci: %d\n" %(num_pairs_predicted_falsely_splitted_loci))
	fReport.write("among these predicted splitted loci, %d of them are non tandem ones, while %d are tandems\n" %(num_pairs_predicted_falsely_non_tandem_splitted_loci, num_pairs_predicted_falsely_splitted_loci-num_pairs_predicted_falsely_non_tandem_splitted_loci))
	fReport.write("total number of pairs of markers being localized in falsely splitted loci: %d\n" %(num_localized_variants))

	return potential_dup_loci, cigar_dict

def generate_pileups(region, assembled_genome_file, bam_file, individual_sample_or_all, pileups_outfile) :
	if individual_sample_or_all == "all" :
		pileups_cmd = "samtools mpileup -d 9999 -f %s -r %s %s" %(assembled_genome_file, region, bam_file)
	else :
		pileups_cmd = "samtools mpileup -f %s -r %s %s" %(assembled_genome_file, region, bam_file)
	fDEVNULL = open(devnull, 'wb')
	proc = Popen(shlex.split(pileups_cmd), stdout=file(pileups_outfile, 'w'), stderr=fDEVNULL, shell=False)
	proc.communicate()
	if proc.returncode != 0 :
		stderr.write(timestamper() + " [SAMtools Error] : something wrong when generate pileup file %s\n" %(pileups_outfile))
		exit()

def compute_median_cov(pileups_file) :
	""" return the median cov of a given pileups file """
	covs = []
	median_cov = 0.0
	fPILEUP = open(pileups_file, 'r')
	for line in fPILEUP :
		tmp_line = re.split('\t', line.strip())
		covs.append(int(tmp_line[3]))
	fPILEUP.close()
	if len(covs) > 0 :
		sort_covs = sorted(covs)
		if len(sort_covs) % 2 == 0 :
			median_cov = (sort_covs[len(sort_covs)/2-1] + sort_covs[len(sort_covs)/2]) / 2.0
		else :
			median_cov = sort_covs[len(sort_covs)/2]
	return median_cov

def apply_cov_filter(one_locus, the_other_locus, assembled_genome_file, all_sample_bam_file, out_pileups_loc_1, out_pileups_loc_2) :
	"""
	apply coverage filter to a given pair of split loci candidate identified from aligning assembled genome against itself
	exclude pairs that fail to pass the filter in the future analysis
	"""
	generate_pileups(one_locus, assembled_genome_file, all_sample_bam_file, "all", out_pileups_loc_1)
	generate_pileups(the_other_locus, assembled_genome_file, all_sample_bam_file, "all", out_pileups_loc_2)
	if stat(out_pileups_loc_1).st_size == 0 and os.stat(out_pileups_loc_2).st_size == 0 :
		stdout.write("\tLoci coverage: N/A N/A\n" %(out_pileups_loc_1, out_pileups_loc_2))
		system("rm %s" %(out_pileups_loc_1))
		system("rm %s" %(out_pileups_loc_2))
		return 0
	else :
		median_cov_loc_1 = compute_median_cov(out_pileups_loc_1)
		median_cov_loc_2 = compute_median_cov(out_pileups_loc_2)
		stdout.write("\tLoci coverage: %s %s\n" %(median_cov_loc_1, median_cov_loc_2))
		""" hard choice, fix me in the future """
		if median_cov_loc_1 + median_cov_loc_2 < 600 :
			stdout.write("\tCoverage filter: PASS\n")
			return 1
		else :
			stdout.write("\tCoverage filter: FAIL\n")
			return 0

def find_markers_in_locus(locus_id, locus_start, locus_end, locus_strand, variants_dict, smt_genoytpe_dict) :
	"""
	localize variants in the potentially splitted loci
	find markers within each of the identified potentially splitted loci
	for each pair of markers, define an interval whose pileup will be generated
	"""
	stdout.write("\t%s|%d|%d|%s" %(locus_id, locus_start, locus_end, locus_strand))
	tmp_variants_in_locus = []
	for variant, info in variants_dict.iteritems() :
		chrom_id = re.split(':', variant)[0]
		chrom_pos = int(re.split(':', variant)[1])
		if locus_strand == '+' :
			if chrom_id == locus_id and chrom_pos > locus_start + 5 and chrom_pos <= locus_end - 5 :
				tmp_variants_in_locus.append(variant)
		else :
			if chrom_id == locus_id and chrom_pos > locus_end + 5 and chrom_pos <= locus_start - 5 :
				tmp_variants_in_locus.append(variant)

	variants_in_locus = []
	for i in range(len(tmp_variants_in_locus)) :
		if smt_genoytpe_dict.has_key(tmp_variants_in_locus[i]) :
			variants_in_locus.append(tmp_variants_in_locus[i])
	if len(variants_in_locus) == 0 :
		stdout.write("\tno markers found in the locus\n")
	else :
		stdout.write("\t %s\n" %(" ".join(variants_in_locus)))
	return variants_in_locus

def match_markers_in_loci(loc_1_id, loc_1_start, one_locus_end, strand_loc_1, loc_2_id, the_other_locus_start, the_other_locus_end, strand_loc_2, one_or_the_other, variants_in_locus, cigar, variants_dict, assembled_genome_dict) :
	variant_pairs = []
	for i in range(len(variants_in_locus)) :
		tmp_chrom_id = re.split(':', variants_in_locus[i])[0]
		variant_pos_loc_1 = int(re.split(':', variants_in_locus[i])[1])
		stdout.write("\t--Marker: %s %s\n" %(variants_in_locus[i], '|'.join(re.split('\t', variants_dict[variants_in_locus[i]])[0:2])))
		if strand_loc_1 == '+' :
			relative_pos_loc_1 = variant_pos_loc_1 - loc_1_start
		else :
			relative_pos_loc_1 = loc_1_start - variant_pos_loc_1 - 1
		stdout.write("\t\trelative position in one locus: %d\n" %(relative_pos_loc_1))

		"""
		analyze the CIGAR string to parse the alignment
		"""
		tmp_relative_pos_in_the_other_locus = relative_pos_loc_1
		num_match, num_del, num_ins = 0, 0, 0
		match_flag, insertion_flag, deletion_flag, last_flag = 0, 0, 0, ''
		for j in range(len(cigar)) :
			if cigar[j] == 'M' :
				match_flag = 1
			elif cigar[j] == 'I' :
				insertion_flag = 1
			elif cigar[j] == 'D' :
				deletion_flag = 1
			else :
				if match_flag == 1 :
					num_match += int(cigar[j])
					match_flag = 0
					last_flag = 'M'
				elif deletion_flag == 1 :
					num_del += int(cigar[j])
					deletion_flag = 0
					last_flag = 'D'
					if one_or_the_other == 0 :
						tmp_relative_pos_in_the_other_locus -= int(cigar[j])
					else :
						tmp_relative_pos_in_the_other_locus += int(cigar[j])
				elif insertion_flag == 1 :
					num_ins += int(cigar[j])
					insertion_flag = 0
					last_flag = 'I'
					if one_or_the_other == 0 :
						tmp_relative_pos_in_the_other_locus += int(cigar[j])
					else :
						tmp_relative_pos_in_the_other_locus -= int(cigar[j])
			if one_or_the_other == 0 :
				if num_match+num_del >= relative_pos_loc_1 :
					if last_flag == 'M' :
						relative_pos_in_the_other_locus = relative_pos_loc_1 - num_del + num_ins
					elif last_flag == 'D' :
						relative_pos_in_the_other_locus = -1
					break
			else :
				if num_match + num_ins >= relative_pos_loc_1 :
					if last_flag == "M" :
						relative_pos_in_the_other_locus = relative_pos_loc_1 + num_del - num_ins
					elif last_flag == 'I' :
						relative_pos_in_the_other_locus = -1
					break
		stdout.write("\t\trelative position in the other locus: %d\n" %(relative_pos_in_the_other_locus))
		if relative_pos_in_the_other_locus == -1 :
			stdout.write("\tcannot find the relative position in the other locus: %s\n" %(variants_in_locus[i]))
		else :
			if one_or_the_other == 0 :
				if strand_loc_2 == '+' :
					variant_pos_the_other_locus = the_other_locus_start + relative_pos_in_the_other_locus
					if not assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-1] in re.split('\t', variants_dict[variants_in_locus[i]])[0:2] :
						return []
				else :
					variant_pos_the_other_locus = the_other_locus_start - relative_pos_in_the_other_locus - 1
					if not rc(assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-1]) in re.split('\t', variants_dict[variants_in_locus[i]])[0:2] :
						return []
			else :
				variant_pos_the_other_locus = the_other_locus_start + relative_pos_in_the_other_locus
			#print "\t%s:%s_%s:%s" %(loc_1_id, variant_pos_loc_1, loc_2_id, variant_pos_the_other_locus)
			if one_or_the_other == 0 :
				variant_pairs.append("%s:%s:%s|%s:%s:%s" %(loc_1_id, variant_pos_loc_1, strand_loc_1, loc_2_id, variant_pos_the_other_locus, strand_loc_2))
				stdout.write("\t\t%s\t\t%s %s %s\n" %(loc_1_id+':'+str(variant_pos_loc_1), assembled_genome_dict[loc_1_id][variant_pos_loc_1-11:variant_pos_loc_1-1], assembled_genome_dict[loc_1_id][variant_pos_loc_1-1], assembled_genome_dict[loc_1_id][variant_pos_loc_1:variant_pos_loc_1+10]))
				if strand_loc_2 == '+' :
					stdout.write("\t\t%s\t\t%s %s %s\n" %(loc_2_id+':'+str(variant_pos_the_other_locus), assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-11:variant_pos_the_other_locus-1], assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-1], assembled_genome_dict[loc_2_id][variant_pos_the_other_locus:variant_pos_the_other_locus+10]))
				else :
					stdout.write("\t\t%s\t\t%s %s %s\n" %(loc_2_id+':'+str(variant_pos_the_other_locus), rc(assembled_genome_dict[loc_2_id][variant_pos_the_other_locus:variant_pos_the_other_locus+10]), rc((assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-1])), rc(assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-11:variant_pos_the_other_locus-1])))
			else :
				variant_pairs.append("%s:%s:%s|%s:%s:%s" %(loc_2_id, variant_pos_the_other_locus, strand_loc_2, loc_1_id, variant_pos_loc_1, strand_loc_1))
				if strand_loc_1 == '-' :
					stdout.write("\t\t%s %s %s %s\n" %(loc_2_id+":"+str(variant_pos_the_other_locus), rc(assembled_genome_dict[loc_2_id][variant_pos_the_other_locus:variant_pos_the_other_locus+5]), rc(assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-1]), rc(assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-6:variant_pos_the_other_locus-1])))
				else:
					stdout.write("\t\t%s %s %s %s\n" %(loc_2_id+":"+str(variant_pos_the_other_locus), assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-6:variant_pos_the_other_locus-1], assembled_genome_dict[loc_2_id][variant_pos_the_other_locus-1], assembled_genome_dict[loc_2_id][variant_pos_the_other_locus:variant_pos_the_other_locus+5]))
				stdout.write("\t\t%s %s %s %s\n" %(loc_1_id+":"+str(variant_pos_loc_1), assembled_genome_dict[loc_1_id][variant_pos_loc_1-6:variant_pos_loc_1-1], assembled_genome_dict[loc_1_id][variant_pos_loc_1-1], assembled_genome_dict[loc_1_id][variant_pos_loc_1:variant_pos_loc_1+5]))
	return variant_pairs

def investigate_genotypes(variant_pairs, variants_dict, samples, loci_candidate_out_subdir, splitted_sample_bam_dir, assembled_genome_file) :
	"""
	use genotypes, coverage information to determine potential falsely split loci
	"""
	num_variant_pair_with_expected_gt_pattern_between_loci = 0
	proposed_genotype_at_variant_site = {}
	for i in range(len(variant_pairs)) :
		stdout.write("\t--Analyzing %s\n" %(variant_pairs[i]))
		one_locus = re.split('\|', variant_pairs[i])[0]
		loc_1_id = re.split(':', one_locus)[0]
		marker_pos_loc_1 = int(re.split(':', one_locus)[1])
		strand_loc_1 = re.split(':', one_locus)[2]
		the_other_locus = re.split('\|', variant_pairs[i])[1]
		loc_2_id = re.split(':', the_other_locus)[0]
		marker_pos_loc_2 = int(re.split(':', the_other_locus)[1])
		strand_loc_2 = re.split(':', the_other_locus)[2]
		tmp_genotypes = []
		variant_in_one_or_the_other = ""
		if variants_dict.has_key(loc_1_id+":"+str(marker_pos_loc_1)) :
			variant_in_one_or_the_other = "one"
			tmp_genoytpes = re.split('\t', variants_dict[loc_1_id+":"+str(marker_pos_loc_1)])[2:]
			rallele_loc_1 = re.split('\t', variants_dict[loc_1_id+":"+str(marker_pos_loc_1)])[0]
			alt_allele_loc_1 = re.split('\t', variants_dict[loc_1_id+":"+str(marker_pos_loc_1)])[1]
			proposed_genotype_at_variant_site[loc_1_id+':'+str(marker_pos_loc_1)] = []
		elif variants_dict.has_key(loc_2_id+":"+str(marker_pos_loc_2)) :
			variant_in_one_or_the_other = "other"
			tmp_genoytpes = re.split('\t', variants_dict[loc_2_id+":"+str(marker_pos_loc_2)])[2:]
			rallele_loc_2 = re.split('\t', variants_dict[loc_2_id+":"+str(marker_pos_loc_2)])[0]
			alt_allele_loc_2 = re.split('\t', variants_dict[loc_2_id+":"+str(marker_pos_loc_2)])[1]
			proposed_genotype_at_variant_site[loc_2_id+':'+str(marker_pos_loc_2)] = []
		variant_pairs_out_subdir = join(loci_candidate_out_subdir, "%s_%d-%s_%d" %(loc_1_id, marker_pos_loc_1, loc_2_id, marker_pos_loc_2))
		make_dirs_if_needed(variant_pairs_out_subdir)
		num_individual_with_expected_genotype_pattern = 0
		tmp_num_individual = len(samples)
		for j in range(len(samples)) :
			individual_name = re.sub("rils", 'r', samples[j])
			in_bam_file = join(splitted_sample_bam_dir, individual_name+".bam")
			if variant_in_one_or_the_other == "one" :
				individual_out_pileups = join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(loc_2_id, marker_pos_loc_2))
				generate_pileups("%s:%d-%d" %(loc_2_id, marker_pos_loc_2-25, marker_pos_loc_2+25), assembled_genome_file, in_bam_file, "individual", individual_out_pileups)
				genotype_loc_2, cov_loc_2, type_loc_2 = analyze_pileups(marker_pos_loc_2, individual_out_pileups)
				if strand_loc_2 == "-" :
					if genotype_loc_2 != "./." and genotype_loc_2 != "NA/NA" :
						genotype_loc_2 = rc(genotype_loc_2)
				individual_out_pileups = join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(loc_1_id, marker_pos_loc_1))
				genotype_loc_1, cov_loc_1, type_loc_1 = parse_individual_genotype(tmp_genoytpes[j], rallele_loc_1, alt_allele_loc_1)
				if genotype_loc_1 == genotype_loc_2 and cov_loc_1 + cov_loc_2 <= 7 and cov_loc_1 + cov_loc_2 > 0 :
					stdout.write("\t\t %s %s %s %s %s agree\n" %(samples[j], genotype_loc_1, cov_loc_1, genotype_loc_2, cov_loc_2))
					num_individual_with_expected_genotype_pattern += 1
					if type_loc_1 == "ref" :
						proposed_genotype_at_variant_site[loc_1_id+':'+str(marker_pos_loc_1)].append("0/0:%d" %(cov_loc_2+cov_loc_1))
					else :
						proposed_genotype_at_variant_site[loc_1_id+':'+str(marker_pos_loc_1)].append("1/1:%d" %(cov_loc_2+cov_loc_1))
				elif genotype_loc_2 == "./." and not genotype_loc_1 in ["./.", "NA/NA"] and cov_loc_2 == 0 and cov_loc_1 <= 8 and cov_loc_1 > 0 :
					stdout.write("\t\t %s %s %s %s %s agree\n" %(samples[j], genotype_loc_1, cov_loc_1, genotype_loc_2, cov_loc_2))
					num_individual_with_expected_genotype_pattern += 1
					if type_loc_1 == "ref" :
						proposed_genotype_at_variant_site[loc_1_id+':'+str(marker_pos_loc_1)].append("0/0:%d" %(cov_loc_1))
					else :
						proposed_genotype_at_variant_site[loc_1_id+':'+str(marker_pos_loc_1)].append("1/1:%d" %(cov_loc_1))
				elif genotype_loc_1 == "./." and not genotype_loc_2 in ["./.", "NA/NA"] and cov_loc_1 == 0 and cov_loc_2 <= 8 and cov_loc_2 > 0 :
					stdout.write("\t\t%s %s %s %s %s agree\n" %(samples[j],  genotype_loc_1, cov_loc_1, genotype_loc_2, cov_loc_2))
					num_individual_with_expected_genotype_pattern += 1
					if type_loc_1 == "ref" :
						proposed_genotype_at_variant_site[loc_1_id+':'+str(marker_pos_loc_1)].append("0/0:%d" %(cov_loc_2))
					else :
						proposed_genotype_at_variant_site[loc_1_id+':'+str(marker_pos_loc_1)].append("1/1:%d" %(cov_loc_2))
				elif genotype_loc_2 == "NA/NA" or genotype_loc_1 == "NA/NA" or (genotype_loc_1 == "./." and genotype_loc_2 == "./."):
					stdout.write("\t\t %s %s %s %s %s N/A\n" %(samples[j], genotype_loc_1, cov_loc_1, genotype_loc_2, cov_loc_2))
					tmp_num_individual -= 1
					proposed_genotype_at_variant_site[loc_1_id+':'+str(marker_pos_loc_1)].append("./.:0")
				else :
					stdout.write("\t\t %s %s %s %s %s disagree\n" %(samples[j], genotype_loc_1, cov_loc_1, genotype_loc_2, cov_loc_2))
					proposed_genotype_at_variant_site[loc_1_id+':'+str(marker_pos_loc_1)].append("./.:0")
			else :
				individual_out_pileups = join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(loc_1_id, marker_pos_loc_1))
				generate_pileups("%s:%d-%d" %(loc_1_id, marker_pos_loc_1-25, marker_pos_loc_1+25), assembled_genome_file, in_bam_file, "individual", individual_out_pileups)
				genotype_loc_1, cov_loc_1, type_loc_1 = analyze_pileups(marker_pos_loc_1, individual_out_pileups)
				if strand_loc_2 == "-" :
					if genotype_loc_1 != "./." and genotype_loc_1 != "NA/NA" :
						genotype_loc_1 = rc(genotype_loc_1)
				individual_out_pileups = join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(loc_2_id, marker_pos_loc_2))
				genotype_loc_2, cov_loc_2, type_loc_2 = parse_individual_genotype(tmp_genoytpes[j], rallele_loc_2, alt_allele_loc_2)
				if genotype_loc_1 == genotype_loc_2 and cov_loc_1 + cov_loc_2 <= 7 and cov_loc_1 + cov_loc_2 > 0 :
					stdout.write("\t\t %s %s %s %s %s agree\n" %(samples[j], genotype_loc_2, cov_loc_2, genotype_loc_1, cov_loc_1))
					num_individual_with_expected_genotype_pattern += 1
					if type_loc_2 == "ref" :
						proposed_genotype_at_variant_site[loc_2_id+':'+str(marker_pos_loc_2)].append("0/0:%d" %(cov_loc_2+cov_loc_1))
					else :
						proposed_genotype_at_variant_site[loc_2_id+':'+str(marker_pos_loc_2)].append("1/1:%d" %(cov_loc_2+cov_loc_1))
				elif genotype_loc_2 == "./." and (not genotype_loc_1 in ["./.", "NA/NA"]) and cov_loc_2 == 0 and cov_loc_1 <= 8 and cov_loc_1 > 0 :
					stdout.write("\t\t %s %s %s %s %s agree\n" %(samples[j], genotype_loc_2, cov_loc_2, genotype_loc_1, cov_loc_1))
					num_individual_with_expected_genotype_pattern += 1
					if type_loc_2 == "ref" :
						proposed_genotype_at_variant_site[loc_2_id+':'+str(marker_pos_loc_2)].append("0/0:%d" %(cov_loc_1))
					else :
						proposed_genotype_at_variant_site[loc_2_id+':'+str(marker_pos_loc_2)].append("1/1:%d" %(cov_loc_1))
				elif genotype_loc_1 == "./." and (not genotype_loc_2 in ["./.", "NA/NA"]) and cov_loc_1 == 0 and cov_loc_2 <= 8 and cov_loc_2 > 0 :
					stdout.write("\t\t %s %s %s %s %s agree\n" %(samples[j], genotype_loc_2, cov_loc_2, genotype_loc_1, cov_loc_1))
					num_individual_with_expected_genotype_pattern += 1
					if type_loc_2 == "ref" :
						proposed_genotype_at_variant_site[loc_2_id+':'+str(marker_pos_loc_2)].append("0/0:%d" %(cov_loc_2))
					else :
						proposed_genotype_at_variant_site[loc_2_id+':'+str(marker_pos_loc_2)].append("1/1:%d" %(cov_loc_2))
				elif genotype_loc_2 == "NA/NA" or genotype_loc_1 == "NA/NA" or (genotype_loc_1 == "./." and genotype_loc_2 == "./."):
					stdout.write("\t\t %s %s %s %s %s N/A\n" %(samples[j], genotype_loc_2, cov_loc_2, genotype_loc_1, cov_loc_1))
					tmp_num_individual -= 1
					proposed_genotype_at_variant_site[loc_2_id+':'+str(marker_pos_loc_2)].append("./.:0")
				else :
					stdout.write("\t\t %s %s %s %s %s disagree\n" %(samples[j], genotype_loc_2, cov_loc_2, genotype_loc_1, cov_loc_1))
					proposed_genotype_at_variant_site[loc_2_id+':'+str(marker_pos_loc_2)].append("./.:0")
		stdout.write("\tReport: %d individuals having expected genotypes out of %d\n" %(num_individual_with_expected_genotype_pattern, tmp_num_individual))
		if tmp_num_individual > len(samples)/2.0 :
			if num_individual_with_expected_genotype_pattern/float(tmp_num_individual) > 0.7 :
				num_variant_pair_with_expected_gt_pattern_between_loci += 1
	if len(variant_pairs) == 1 :
		if num_variant_pair_with_expected_gt_pattern_between_loci == 1 :
			return proposed_genotype_at_variant_site
	else :
		if num_variant_pair_with_expected_gt_pattern_between_loci > len(variant_pairs)-num_variant_pair_with_expected_gt_pattern_between_loci :
			return proposed_genotype_at_variant_site

def parse_individual_genotype(individual_genotype, ref_allele, alt_allele) :
	"""
	given an individual genotype as defined in the INFO field in the VCF file, the function spit out the genotype and coverage at the site #
	"""
	tmp_genotype = re.split(':', individual_genotype)[0]
	genotype = "./."
	cov, total_cov = 0, 0
	type = "NA"
	if tmp_genotype != "./." :
		total_cov = int(re.split(':', individual_genotype)[2])
		num_ref_allele = int(re.split(',', re.split(':', individual_genotype)[1])[0])
		num_alt_allele = int(re.split(',', re.split(':', individual_genotype)[1])[1])
		if total_cov >= 2 :
			if num_ref_allele > num_alt_allele :
				if tmp_genotype == "0/0" :
					genotype = ref_allele + ref_allele
					cov = num_ref_allele
					type = "ref"
				else :
					genotype = "NA/NA"
					cov = -1
			elif num_alt_allele > num_ref_allele :
				if tmp_genotype == "1/1" :
					genotype = alt_allele + alt_allele
					cov = num_alt_allele
					type = "alt"
				else :
					genotype = "NA/NA"
					cov = -1
			elif num_ref_allele == num_alt_allele and num_ref_allele == 1 :
				genotype = "NA/NA"
				cov = -1
			elif num_ref_allele == num_alt_allele and num_ref_allele > 1 :
				if tmp_genotype == "0/0" :
					genotype = ref_allele+ref_allele
					cov = num_ref_allele
					type = "ref"
				elif tmp_genotype == "1/1" :
					genotype = alt_allele+alt_allele
					cov = num_alt_allele
					type = "alt"
				elif tmp_genotype == "0/1" :
					genotype = "NA/NA"
					cov = -1
		else :
			#if num_ref_allele > num_alt_allele and tmp_genotype == "0/0" :
			#	genotype = ref_allele + ref_allele
			#	cov = num_ref_allele
			#	type = "ref"
			#elif num_alt_allele > num_ref_allele and tmp_genotype == "1/1" :
			#	genotype = alt_allele + alt_allele
			#	cov = num_alt_allele
			#	type = "alt"
			#else :
			genotype = "NA/NA"
			cov = -1
	return genotype, cov, type

def analyze_pileups(marker_pos, individual_out_pileups) :
	"""
	given a marker position and a pileup file, spit out the genotype at the site #
	"""
	genotype = "./."
	cov_at_marker_pos = 0
	type = "NA"
	num_ref_allele, num_alt_allele = 0, 0
	fPILEUP = open(individual_out_pileups, 'r')
	for line in fPILEUP :
		tmp_line = re.split('\t', line.strip())
		if marker_pos == int(tmp_line[1]) and int(tmp_line[3]) > 0 :
			pileups = tmp_line[4]
			total_cov = float(tmp_line[3])
			i = 0
			alt_allele_dict = {}
			while i < len(pileups) :
				if pileups[i] in ['.', ','] :
					num_ref_allele += 1
					i += 1
				elif pileups[i] == '$' :
					i += 1
					continue
				elif pileups[i] == '^' :
					i += 2
				elif pileups[i] in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'] :
					if not alt_allele_dict.has_key(pileups[i].upper()) :
						alt_allele_dict[pileups[i].upper()] = 1
					else :
						alt_allele_dict[pileups[i].upper()] += 1
					i += 1
				elif pileups[i] in ['+', '-'] :
					len_indel = int(pileups[i+1])
					indel = ""
					for j in range(len_indel) :
						indel += pileups[i+2+j]
					if not alt_allele_dict.has_key(indel.upper()) :
						alt_allele_dict[indel.upper()] = 1
					else :
						alt_allele_dict[indel.upper()] += 1
					i += 2+len_indel
			alt_allele = ""
			for tmp_alt_allele, count in alt_allele_dict.iteritems() :
				if num_alt_allele == 0 :
					num_alt_allele = count
					alt_allele = tmp_alt_allele
				else :
					if num_alt_allele < count :
						alt_allele = tmp_alt_allele
						num_alt_allele = count
			if total_cov >= 2 :
				if num_ref_allele > num_alt_allele and num_ref_allele/total_cov >= 0.8 :
					genotype = tmp_line[2]+tmp_line[2]
					cov_at_marker_pos = num_ref_allele
					type = "ref"
				elif num_alt_allele > num_ref_allele and num_alt_allele/total_cov >= 0.8 :
					genotype = alt_allele + alt_allele
					cov_at_marker_pos = num_alt_allele
					type = "alt"
				else :
					genotype = "NA/NA"
					cov_at_marker_pos = -1
			else :
				genotype = "NA/NA"
				cov_at_marker_pos = -1
	fPILEUP.close()
	return genotype, cov_at_marker_pos, type

def rc(seq) :
	""" return reversed and complemented sequence of the give sequence """
	rc_seq = ""
	i = len(seq) - 1
	while i >= 0 :
		if seq[i] == 'A' :
			rc_seq += 'T'
		elif seq[i] == 'T' :
			rc_seq += 'A'
		elif seq[i] == 'C' :
			rc_seq += 'G'
		elif seq[i] == 'G' :
			rc_seq += 'C'
		elif seq[i] == 'N' :
			rc_seq += 'N'
		i -= 1
	return rc_seq

def get_genotypes(gatk_vcf_file, smt_vcf_file) :
	""" get the genotypes from a given VCF file """
	stderr.write(timestamper() +  " Getting genotypes from GATK VCf file %s\n" %(gatk_vcf_file))
	samples = []
	variants_dict = {}
	num_variants = 0
	fGATK = open(gatk_vcf_file, 'r')
	for line in fGATK :
		if line.startswith('#') :
			if line.startswith("#CHROM") :
				samples = re.split('\t', line.strip())[9:]
		else :
			tmp_line = re.split('\t', line.strip())
			chrom_id = tmp_line[0]
			chrom_pos = tmp_line[1]
			ref_base = tmp_line[3]
			alt_base = tmp_line[4]
			if re.findall("DP=\d{1,}", tmp_line[7]) :
				dp = re.findall("DP=\d{1,}", tmp_line[7])[0]
				tmp_dp = 0
				num_het = 0
				tmp_num_individual = len(samples)
				genotypes = tmp_line[9:]
				for i in range(len(genotypes)) :
					tmp_genotype_per_individual = re.split(':', genotypes[i])[0]
					if tmp_genotype_per_individual == "./." :
						tmp_num_individual -= 1
					else :
						if tmp_genotype_per_individual == "0/1" or tmp_genotype_per_individual == "1/2" or tmp_genotype_per_individual == "0/2" :
							num_het += 1
						tmp_dp_per_individual = int(re.split(':', genotypes[i])[2])
						tmp_dp += tmp_dp_per_individual
				""" put this hard filter to command line option, fix me in the future """
				if num_het <= 4 :
					tmp_genotypes = '\t'.join(tmp_line[9:])
					num_variants += 1
					if not variants_dict.has_key(chrom_id+':'+chrom_pos) :
						variants_dict[chrom_id+':'+chrom_pos] = "%s\t%s\t%s" %(ref_base, alt_base, tmp_genotypes)
					else :				# there should not have any redundant calls #
						stderr.write(timestamper() + " redundant calls found: %s %s\n" %(chrom_id, chrom_pos))
						stderr.write(timestamper() + " please run fix_redundant_calls.py first to solve the problem\n")
						exit()
	fGATK.close()
	stderr.write(timestamper() + " genotypes of %d individuals at %d variant sites parsed\n" %(len(samples), num_variants))

	# read in genotypes by samtools #
	stderr.write(timestamper() +  " Getting genotypes from SAMtools VCf file %s\n" %(smt_vcf_file))
	fSMT = open(smt_vcf_file, 'r')
	num_variants = 0
	smt_genoytpe_dict = {}
	for line in fSMT :
		if not line.startswith('#') :
			tmp_line = re.split('\t', line.strip())
			chrom_id = tmp_line[0]
			variant_pos = tmp_line[1]
			ref_allele = tmp_line[3]
			alt_allele = tmp_line[4]
			smt_genoytpe_dict[chrom_id+':'+variant_pos] = "%s\t%s" %(ref_allele, alt_allele)
			num_variants += 1
	fSMT.close()
	stderr.write(timestamper() + " genotypes of %d individuals at %d variant sites parsed\n" %(len(samples), num_variants))
	return variants_dict, samples, smt_genoytpe_dict

def read_assembly(assembled_genome_file) :
	""" read the assembled genome """
	stderr.write(timestamper() + " Reading the assembled genome %s\n" %(assembled_genome_file))
	header, seq = "", ""
	num_seq = 0
	assembled_genome_dict = {}
	fREF = open(assembled_genome_file, 'r')
	for line in fREF :
		if line.startswith('>') :
			if seq != "" :
				assembled_genome_dict[header] = seq
				num_seq += 1
				seq = ""
			header = line.strip()[1:]
		else :
			seq += line.strip().upper()
	if seq != "" :
		assembled_genome_dict[header] = seq
		num_seq += 1
	fREF.close()
	stderr.write(timestamper() + " %d sequences parsed\n" %(num_seq))
	return assembled_genome_dict

def timestamper() :
	""" generate a dash-separated time stamped string """
	return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

def check_files_existence(*files) :
	""" check file existence """
	for file in files :
		if not exists(file) :
			stderr.write(timestamper() + " [IO Error]: Cannot find the file you provided %s\n" %(file))
			exit()

def check_dirs_existence(*dirs) :
	""" check input directory existence, not making any upon non-existence """
	for dir in dirs :
		if not exists(dir) :
			stderr.write(timestamper() + " [IO Error]: Cannot find the input directory you provided %s\n" %(dir))
			exit()

def make_dirs_if_needed(*dirs) :
	""" make any directories if necessary """
	for dir in dirs :
		if not exists(dir) :
			makedirs(dir)

def main(alignment_file, assembled_genome_file, gatk_vcf_file, smt_vcf_file, all_sample_bam_file, splitted_sample_bam_dir, out_dir, min_aln_len, min_idn) :
	assembled_genome_dict = read_assembly(assembled_genome_file)
	variants_dict, samples, smt_genotype_dict = get_genotypes(gatk_vcf_file, smt_vcf_file)
	split_loci_candidates, cigar_dict = identify_splited_loci(alignment_file, assembled_genome_file, all_sample_bam_file, out_dir, splitted_sample_bam_dir, variants_dict, samples, smt_genotype_dict, assembled_genome_dict, min_aln_len, min_idn/100.0)

if __name__ == "__main__" :
	parser = ArgumentParser(description="Identify alleles/loci that are erroneously split by assemblers")
	parser.add_argument("-assembly", metavar="FILE", dest="assembled_genome_file", required=True, help="Specify assembled genome in fasta")
	parser.add_argument("-aln", metavar="FILE", dest="alignment_file", required=True, help="Specify the file with alignment results of aligning assembled genome against itself. Note: lastz alignment results in CIGAR format is supported only in current version")
	parser.add_argument("-gatk_vcf", metavar="FILE", dest="gatk_vcf_file", required=True, help="Specify the VCF file with genotypes of the mapping population obtained from GATK")
	parser.add_argument("-bam_file", metavar="FILE", dest="all_sample_bam_file", required=True, help="Specify the BAM file of individuals of the entire mapping population")
	parser.add_argument("-bam_dir", metavar="DIR", dest="splitted_sample_bam_dir", required=True, help="Specify the directory with individual BAM files (not merged as a single BAM file like what -bam_file requires)")
	parser.add_argument("-smt_vcf", metavar="FILE", dest="smt_vcf_file", required=True, help="Specify the VCF file with genotypes of the mapping population obtained from SAMtools")
	parser.add_argument("-out_dir", metavar="DIR", dest="out_dir", required=True, help="Specifying the output directory under which multiple subfolders will be created")
	parser.add_argument("-min_len", metavar="INT", dest="min_aln_len", type=int, default=1000, help="Specify the minimum alignment length of a pair of loci/alleles to be considered. DEFAULT: 1000")
	parser.add_argument("-min_idn", metavar="INT", dest="min_idn", type=int, default=90, help="Specify the minimum alignment identity (e.g. 90 mean 90%% identity) of a pair of loci/alleles to be considered. DEFAULT: 90")

	args = parser.parse_args()

	check_files_existence(args.alignment_file, args.assembled_genome_file, args.gatk_vcf_file, args.smt_vcf_file)
	check_dirs_existence(args.splitted_sample_bam_dir)
	make_dirs_if_needed(args.out_dir)

	"""
	print out command line settings
	"""
	stderr.write("\npython %s" %(argv[0]))
	for k, v in vars(args).iteritems() :
		stderr.write(" -%s %s" %(k, v))
	stderr.write("\n\n")
	stderr.write("Settings:\n")
	stderr.write("\tassembly: \t\t%s\n" %(args.assembled_genome_file))
	stderr.write("\talignment: \t\t%s\n" %(args.alignment_file))
	stderr.write("\tGATK VCF: \t\t%s\n" %(args.gatk_vcf_file))
	stderr.write("\tSAMtools VCF: \t\t%s\n" %(args.smt_vcf_file))
	stderr.write("\tBAM file: \t\t%s\n" %(args.all_sample_bam_file))
	stderr.write("\tBAM dir: \t\t%s\n" %(args.splitted_sample_bam_dir))
	stderr.write("\toutput dir: \t\t%s\n" %(args.out_dir))
	stderr.write("\tmin alignment length: \t%d\n" %(args.min_aln_len))
	stderr.write("\tmin alignment identity: %d\n\n" %(args.min_idn))

	main(args.alignment_file, args.assembled_genome_file, args.gatk_vcf_file, args.smt_vcf_file, args.all_sample_bam_file, args.splitted_sample_bam_dir, args.out_dir, args.min_aln_len, args.min_idn)
