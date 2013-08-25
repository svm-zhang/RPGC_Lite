# The script is used to identify falsely splitted loci of a given assembly using alignment, genotype information, coverage, and insert size #

import os, sys
import re
import argparse

# get a set of pairs of loci to be candidates of falsely splitted loci, given an alignment of an assembly against itself #
#def Get_Splitted_loci_Candidates(alignment_file, assembled_genome_file, all_sample_bam_file, loci_wide_pileups_dir, markers_wide_pileups_dir, splitted_sample_bam_dir, variants_dict, samples, assembled_genome_dict) :
def Get_Splitted_loci_Candidates(alignment_file, assembled_genome_file, all_sample_bam_file, out_dir, splitted_sample_bam_dir, variants_dict, samples, smt_genotype_dict, assembled_genome_dict) :

	# setting up sub-directories for pileups		  #
	# 1. pileups around markers						  #
	# 2. pileups of a pair of splitted loci candidate #
	pileups_dir = os.path.join(out_dir, "pileups")
	markers_wide_pileups_dir = os.path.join(pileups_dir, "markers_wide")
	loci_wide_pileups_dir = os.path.join(pileups_dir, "loci_wide")
	if not os.path.exists(markers_wide_pileups_dir) :
		os.makedirs(markers_wide_pileups_dir)
	if not os.path.exists(loci_wide_pileups_dir) :
		os.makedirs(loci_wide_pileups_dir)

	# setting up sub-directory for sequences #
	# sequences of pairs of splitted loci candidates with(out) markers #
	seqs_dir = os.path.join(out_dir, "seqs")
	if not os.path.exists(seqs_dir) :
		os.makedirs(seqs_dir)
	splitted_loci_with_markers_seq_outfile = os.path.join(seqs_dir, "splitted_loci_with_markers.fasta")
	splitted_loci_without_markers_seq_outfile = os.path.join(seqs_dir, "splitted_loci_without_markers.fasta")
	fSEQ_LOCI_MARKER = open(splitted_loci_with_markers_seq_outfile, 'w')
	fSEQ_LOCI_NO_MARKER = open(splitted_loci_without_markers_seq_outfile, 'w')

	# setting up subdirectory for summary #
	summary_dir = os.path.join(out_dir, "summary")
	if not os.path.exists(summary_dir) :
		os.makedirs(summary_dir)
	proposed_genotypes_outfile = os.path.join(summary_dir, "non_tandem_splitted_loci.proposed_genotypes")
	fProPose = open(proposed_genotypes_outfile, 'w')

	threshold_len_alignment, threshold_identity_alignment = 1000, 0.9
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
	total_num_loci, num_loci_passed_len_idn_filter = 0, 0
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
			################################################################
			# get the length and the identity of the alignment using CIGAR #
			################################################################
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
			if len_alignment >= threshold_len_alignment and identity_alignment >= threshold_identity_alignment :		# apply alignment length and identity filters #
				t_id = tmp_line[5]
				q_id = tmp_line[1]
				one_locus = "%s:%d-%d" %(t_id, t_start, t_end)
				if q_strand == '-' :
					the_other_locus = "%s:%d-%d" %(q_id, q_end, q_start)
				else :
					the_other_locus = "%s:%d-%d" %(q_id, q_start, q_end)
				if not all_splitted_loci_dict.has_key(the_other_locus+'|'+one_locus) :
					all_splitted_loci_dict[one_locus+'|'+the_other_locus] = 1
					num_loci_passed_len_idn_filter += 1
					print ">", one_locus, the_other_locus
					#########################
					# apply coverage filter #
					#########################
					proposed_genotype_at_variant_site = {}
					tmp_one_locus = re.sub(':', '-', re.sub('_', '', one_locus))
					one_locus_pileups_out = os.path.join(loci_wide_pileups_dir, tmp_one_locus+".pileups")
					tmp_the_other_locus = re.sub(':', '-', re.sub('_', '', the_other_locus))
					the_other_locus_pileups_out = os.path.join(loci_wide_pileups_dir, tmp_the_other_locus+".pileups")
					pass_or_filter = Applying_Cov_Filter(one_locus, the_other_locus, assembled_genome_file, all_sample_bam_file, one_locus_pileups_out, the_other_locus_pileups_out)
					if pass_or_filter :
						###########################################
						# Localizing markers in each pair of loci #
						###########################################
						variants_in_one_locus = Localize_Markers_in_Locus(t_id, t_start, t_end, t_strand, variants_dict, smt_genotype_dict)
						variants_in_the_other_locus = Localize_Markers_in_Locus(q_id, q_start, q_end, q_strand, variants_dict, smt_genotype_dict)

						#############################################
						# Find each pair of markers in each of loci #
						#############################################
						variant_pairs = []
						if len(variants_in_one_locus) > 0 :
							one_or_the_other = 0
							variant_pairs += Get_Variant_pair_in_Splitted_Loci_Candidates(t_id, t_start, t_end, t_strand, q_id, q_start, q_end, q_strand, one_or_the_other, variants_in_one_locus, cigar, variants_dict, assembled_genome_dict)
						if len(variants_in_the_other_locus) > 0 :
							one_or_the_other = 1
							variant_pairs += Get_Variant_pair_in_Splitted_Loci_Candidates(q_id, q_start, q_end, q_strand, t_id, t_start, t_end, t_strand, one_or_the_other, variants_in_the_other_locus, cigar, variants_dict, assembled_genome_dict)

						##############################################
						# Creating a database of variants that are   #
						# identified in each pair of loci candidates #
						##############################################
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

						num_localized_variants += len(variants_in_one_locus) + len(variants_in_the_other_locus)
						if len(variants_in_one_locus) > 0 or len(variants_in_the_other_locus) > 0 :
							fSEQ_LOCI_MARKER.write(">%s_%d_%d_1\n%s\n" %(t_id, t_start, t_end, assembled_genome_dict[t_id][t_start-1:t_end]))
							if q_strand == '+' :
								fSEQ_LOCI_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, assembled_genome_dict[q_id][q_start-1:q_end]))
							else :
								fSEQ_LOCI_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, RC(assembled_genome_dict[q_id][q_end-1:q_start])))
							loci_candidate_out_subdir = os.path.join(markers_wide_pileups_dir, "%s_%d_%d-%s_%d_%d" %(t_id, t_start, t_end, q_id, q_start, q_end))
							if not os.path.exists(loci_candidate_out_subdir) :
								os.makedirs(loci_candidate_out_subdir)
							proposed_genotype_at_variant_site = Genotypes_Investigation(variant_pairs, variants_dict, samples, loci_candidate_out_subdir, splitted_sample_bam_dir, assembled_genome_file)
							if proposed_genotype_at_variant_site != None :
								print "\tproposed genotype:", proposed_genotype_at_variant_site
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
							print "\tnumber of pairs of variants can be localized:", num_localized_variants
							#print "\tdatabase", variants_in_loci_candidates_db
						else :
							loci_without_markers_localized.append("%s|%s" %(one_locus, the_other_locus))
							fSEQ_LOCI_NO_MARKER.write(">%s_%d_%d_1\n%s\n" %(t_id, t_start, t_end, assembled_genome_dict[t_id][t_start-1:t_end]))
							if q_strand == '+' :
								fSEQ_LOCI_NO_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, assembled_genome_dict[q_id][q_start-1:q_end]))
							else :
								fSEQ_LOCI_NO_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, RC(assembled_genome_dict[q_id][q_end-1:q_start])))
							num_loci_with_no_localized_variants += 1
					else :
						num_loci_failed_cov_filter += 1
						print "\tfailed the coverage filter", one_locus, the_other_locus
						os.system("rm %s" %(one_locus_pileups_out))
						os.system("rm %s" %(the_other_locus_pileups_out))
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

	print "#########################################################"
	print "predicted splitted loci", predicted_falsely_splitted_loci
	print "variants database", variants_in_loci_candidates_db
	out_summary_file = os.path.join(summary_dir, "splitted_loci_identification.summary")
	fSummary = open(out_summary_file, 'w')
	fSummary.write("total number of loci in the alignment: %d\n" %(total_num_loci))
	fSummary.write("number of loci alignment having at least 1000 bp long and 90%% identity: %d\n" %(num_loci_passed_len_idn_filter))
	fSummary.write("number of loci alignment who passed the length and identity filter failing to pass the coverage filter: %d\n" %(num_loci_failed_cov_filter))
	fSummary.write("number of loci alignment who passed all three filters having no localized variants: %d\n" %(num_loci_with_no_localized_variants))
	fSummary.write("number of pairs of loci being identified as falsely splitted loci: %d\n" %(num_pairs_predicted_falsely_splitted_loci))
	fSummary.write("among these predicted splitted loci, %d of them are non tandem ones, while %d are tandems\n" %(num_pairs_predicted_falsely_non_tandem_splitted_loci, num_pairs_predicted_falsely_splitted_loci-num_pairs_predicted_falsely_non_tandem_splitted_loci))
	fSummary.write("total number of pairs of markers being localized in falsely splitted loci: %d\n" %(num_localized_variants))

	# if verbose, I should output this to stdout #
	#print "total number of loci in the alignment:", total_num_loci
	#print "number of loci alignment having at least 1000 bp long and 90% identity:", num_loci_passed_len_idn_filter
	#print "number of loci alignment who passed the length and identity filter failing to pass the coverage filter:", num_loci_failed_cov_filter
	#print "number of loci alignment who passed all three filters having no localized variants:", num_loci_with_no_localized_variants
	#print "number of pairs of loci being identified as falsely splitted loci: %d" %(num_pairs_predicted_falsely_splitted_loci)
	#print "among these predicted splitted loci, %d of them are non tandem ones, while %d are tandems" %(num_pairs_predicted_falsely_non_tandem_splitted_loci, num_pairs_predicted_falsely_splitted_loci-num_pairs_predicted_falsely_non_tandem_splitted_loci)
	#print "total number of markers being localized in falsely splitted loci:", num_localized_variants
	return potential_dup_loci, cigar_dict

def Pileups_Generator(region, assembled_genome_file, bam_file, individual_sample_or_all, pileups_outfile) :
	if individual_sample_or_all == "all" :
		pileups_cmd = "samtools mpileup -d 9999 -f %s -r %s %s > %s" %(assembled_genome_file, region, bam_file, pileups_outfile)
	else :
		pileups_cmd = "samtools mpileup -f %s -r %s %s > %s" %(assembled_genome_file, region, bam_file, pileups_outfile)
	#print pileups_cmd
	os.system(pileups_cmd)

# Given a pileup file of a region, the function gives you the median coverage #
def Compute_Regional_Median_Coverage(pileups_file) :
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

def Applying_Cov_Filter(one_locus, the_other_locus, assembled_genome_file, all_sample_bam_file, one_locus_pileups_out, the_other_locus_pileups_out) :
	Pileups_Generator(one_locus, assembled_genome_file, all_sample_bam_file, "all", one_locus_pileups_out)
	Pileups_Generator(the_other_locus, assembled_genome_file, all_sample_bam_file, "all", the_other_locus_pileups_out)
	if os.stat(one_locus_pileups_out).st_size == 0 and os.stat(the_other_locus_pileups_out).st_size == 0 :
		print "no coverage in both loci", one_locus_pileups_out, the_other_locus_pileups_out
		os.system("rm %s" %(one_locus_pileups_out))
		os.system("rm %s" %(the_other_locus_pileups_out))
		return 0
	else :
		median_cov_one_locus = Compute_Regional_Median_Coverage(one_locus_pileups_out)
		median_cov_the_other_locus = Compute_Regional_Median_Coverage(the_other_locus_pileups_out)
		print "\t", median_cov_one_locus, median_cov_the_other_locus
		if median_cov_one_locus + median_cov_the_other_locus < 600 :
			return 1
		else :
			return 0

# get the variants information for the given vcf file #
######################
# parse the VCF file #
######################
def VCF_Parser(vcf_file) :
	print "parsing the vcf file ..."
	samples = []
	variants_dict = {}
	num_variants = 0
	fVCF = open(vcf_file, 'r')
	for line in fVCF :
		if line.startswith('#') :
			if line.startswith("#CHROM") :
				samples = re.split('\t', line.strip())[9:]
				print "\t%d" %(len(samples)), "parsed"
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
				#if tmp_num_individual/float(len(samples)) >= 0.85 :
				#if num_het/float(tmp_num_individual) <= 0.05 :
				if num_het <= 4 :
					tmp_genotypes = '\t'.join(tmp_line[9:])
					num_variants += 1
					if not variants_dict.has_key(chrom_id+':'+chrom_pos) :
						variants_dict[chrom_id+':'+chrom_pos] = "%s\t%s\t%s" %(ref_base, alt_base, tmp_genotypes)
					else :				# there should not have any redundant calls #
						print "redundant calls:", chrom_id, chrom_pos
						sys.exit(1)
	print "\t%d" %(num_variants), "variants parsed"

	# read in genotypes by samtools #
	smt_genotypes_file = "/N/dc/projects/simozhan/rpgc/c_elegan/genotyping/assembly213/genotypes/multi_sm_1_70.smt.genotypes.final.vcf"
	fSMT = open(smt_genotypes_file, 'r')
	smt_genoytpe_dict = {}
	for line in fSMT :
		if not line.startswith('#') :
			tmp_line = re.split('\t', line.strip())
			chrom_id = tmp_line[0]
			variant_pos = tmp_line[1]
			ref_allele = tmp_line[3]
			alt_allele = tmp_line[4]
			smt_genoytpe_dict[chrom_id+':'+variant_pos] = "%s\t%s" %(ref_allele, alt_allele)
	fSMT.close()
	return variants_dict, samples, smt_genoytpe_dict

######################################################
# Localize variants in the potentially splitted loci #
######################################################
# find markers within each of the identified potentially splitted loci #
# for each pair of markers, define an interval whose pileup will be generated #
#def Localize_Markers_in_Locus(split_loci_candidates, cigar_dict, variants_dict, assembled_genome_dict) :
def Localize_Markers_in_Locus(locus_id, locus_start, locus_end, locus_strand, variants_dict, smt_genoytpe_dict) :
	sys.stdout.write("\t%s|%d|%d|%s" %(locus_id, locus_start, locus_end, locus_strand))
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
		print "\tcannot find markers in the locus"
	else :
		print "\t", variants_in_locus
	return variants_in_locus

def Get_Variant_pair_in_Splitted_Loci_Candidates(one_locus_id, one_locus_start, one_locus_end, one_locus_strand, the_other_locus_id, the_other_locus_start, the_other_locus_end, the_other_locus_strand, one_or_the_other, variants_in_locus, cigar, variants_dict, assembled_genome_dict) :
	variant_pairs = []
	for i in range(len(variants_in_locus)) :
		tmp_chrom_id = re.split(':', variants_in_locus[i])[0]
		variant_pos_one_locus = int(re.split(':', variants_in_locus[i])[1])
		print "\tvariant:", variants_in_locus[i], '|'.join(re.split('\t', variants_dict[variants_in_locus[i]])[0:2])
		if one_locus_strand == '+' :
			relative_pos_in_one_locus = variant_pos_one_locus - one_locus_start
		else :
			relative_pos_in_one_locus = one_locus_start - variant_pos_one_locus - 1
		print "\trelative in the first locus:", relative_pos_in_one_locus
		# analyzing the CIGAR string to get the details of alignment #
		tmp_relative_pos_in_the_other_locus = relative_pos_in_one_locus
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
						#num_match += int(cigar[j])
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
				if num_match+num_del >= relative_pos_in_one_locus :
					if last_flag == 'M' :
						relative_pos_in_the_other_locus = relative_pos_in_one_locus - num_del + num_ins
					elif last_flag == 'D' :
						relative_pos_in_the_other_locus = -1
					break
			else :
				if num_match + num_ins >= relative_pos_in_one_locus :
					if last_flag == "M" :
						relative_pos_in_the_other_locus = relative_pos_in_one_locus + num_del - num_ins
					elif last_flag == 'I' :
						relative_pos_in_the_other_locus = -1
					break
		print "\trelative position in the other locus:", relative_pos_in_the_other_locus, tmp_relative_pos_in_the_other_locus
		if relative_pos_in_the_other_locus == -1 :
			print "\tcannot find the relative position in the other locus:", variants_in_locus[i]
		else :
			if one_or_the_other == 0 :
				if the_other_locus_strand == '+' :
					variant_pos_the_other_locus = the_other_locus_start + relative_pos_in_the_other_locus
					if not assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1] in re.split('\t', variants_dict[variants_in_locus[i]])[0:2] :
						return []
				else :
					variant_pos_the_other_locus = the_other_locus_start - relative_pos_in_the_other_locus - 1
					if not RC(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1]) in re.split('\t', variants_dict[variants_in_locus[i]])[0:2] :
						return []
			else :
				variant_pos_the_other_locus = the_other_locus_start + relative_pos_in_the_other_locus
			print "\t%s:%s_%s:%s" %(one_locus_id, variant_pos_one_locus, the_other_locus_id, variant_pos_the_other_locus)
			if one_or_the_other == 0 :
				variant_pairs.append("%s:%s:%s|%s:%s:%s" %(one_locus_id, variant_pos_one_locus, one_locus_strand, the_other_locus_id, variant_pos_the_other_locus, the_other_locus_strand))
				print "\t>%s" %(one_locus_id+':'+str(variant_pos_one_locus)), "\t\t", assembled_genome_dict[one_locus_id][variant_pos_one_locus-11:variant_pos_one_locus-1], assembled_genome_dict[one_locus_id][variant_pos_one_locus-1], assembled_genome_dict[one_locus_id][variant_pos_one_locus:variant_pos_one_locus+10]
				if the_other_locus_strand == '+' :
					print "\t>%s" %(the_other_locus_id+':'+str(variant_pos_the_other_locus)), "\t\t", assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-11:variant_pos_the_other_locus-1], assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1], assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus:variant_pos_the_other_locus+10]
				else :
					print "\t>%s" %(the_other_locus_id+':'+str(variant_pos_the_other_locus)), "\t\t", RC(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus:variant_pos_the_other_locus+10]), RC((assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1])), RC(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-11:variant_pos_the_other_locus-1])
			else :
				variant_pairs.append("%s:%s:%s|%s:%s:%s" %(the_other_locus_id, variant_pos_the_other_locus, the_other_locus_strand, one_locus_id, variant_pos_one_locus, one_locus_strand))
				if one_locus_strand == '-' :
					print "\t>%s" %(the_other_locus_id+":"+str(variant_pos_the_other_locus)), RC(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus:variant_pos_the_other_locus+5]), RC(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1]), RC(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-6:variant_pos_the_other_locus-1])
				else:
					print "\t>%s" %(the_other_locus_id+":"+str(variant_pos_the_other_locus)), assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-6:variant_pos_the_other_locus-1], assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1], assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus:variant_pos_the_other_locus+5]
				print "\t>%s" %(one_locus_id+":"+str(variant_pos_one_locus)), assembled_genome_dict[one_locus_id][variant_pos_one_locus-6:variant_pos_one_locus-1], assembled_genome_dict[one_locus_id][variant_pos_one_locus-1], assembled_genome_dict[one_locus_id][variant_pos_one_locus:variant_pos_one_locus+5]
	return variant_pairs

def Genotypes_Investigation(variant_pairs, variants_dict, samples, loci_candidate_out_subdir, splitted_sample_bam_dir, assembled_genome_file) :
	num_variant_pair_with_expected_gt_pattern_between_loci = 0
	proposed_genotype_at_variant_site = {}
	for i in range(len(variant_pairs)) :
		print "\t\t***", variant_pairs[i]
		one_locus = re.split('\|', variant_pairs[i])[0]
		one_locus_id = re.split(':', one_locus)[0]
		marker_pos_one_locus = int(re.split(':', one_locus)[1])
		one_locus_strand = re.split(':', one_locus)[2]
		the_other_locus = re.split('\|', variant_pairs[i])[1]
		the_other_locus_id = re.split(':', the_other_locus)[0]
		marker_pos_the_other_locus = int(re.split(':', the_other_locus)[1])
		the_other_locus_strand = re.split(':', the_other_locus)[2]
		tmp_genotypes = []
		variant_in_one_or_the_other = ""
		if variants_dict.has_key(one_locus_id+":"+str(marker_pos_one_locus)) :
			variant_in_one_or_the_other = "one"
			tmp_genoytpes = re.split('\t', variants_dict[one_locus_id+":"+str(marker_pos_one_locus)])[2:]
			ref_allele_one_locus = re.split('\t', variants_dict[one_locus_id+":"+str(marker_pos_one_locus)])[0]
			alt_allele_one_locus = re.split('\t', variants_dict[one_locus_id+":"+str(marker_pos_one_locus)])[1]
			print ref_allele_one_locus, alt_allele_one_locus, tmp_genoytpes
			proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)] = []
		elif variants_dict.has_key(the_other_locus_id+":"+str(marker_pos_the_other_locus)) :
			variant_in_one_or_the_other = "other"
			tmp_genoytpes = re.split('\t', variants_dict[the_other_locus_id+":"+str(marker_pos_the_other_locus)])[2:]
			ref_allele_the_other_locus = re.split('\t', variants_dict[the_other_locus_id+":"+str(marker_pos_the_other_locus)])[0]
			alt_allele_the_other_locus = re.split('\t', variants_dict[the_other_locus_id+":"+str(marker_pos_the_other_locus)])[1]
			print ref_allele_the_other_locus, alt_allele_the_other_locus, tmp_genoytpes
			proposed_genotype_at_variant_site[the_other_locus_id+':'+str(marker_pos_the_other_locus)] = []
		variant_pairs_out_subdir = os.path.join(loci_candidate_out_subdir, "%s_%d-%s_%d" %(one_locus_id, marker_pos_one_locus, the_other_locus_id, marker_pos_the_other_locus))
		if not os.path.exists(variant_pairs_out_subdir) :
			os.makedirs(variant_pairs_out_subdir)
		num_individual_with_expected_genotype_pattern = 0
		tmp_num_individual = len(samples)
		for j in range(len(samples)) :
			individual_name = re.sub("rils", 'r', samples[j])
			in_bam_file = os.path.join(splitted_sample_bam_dir, individual_name+".bam")
			if variant_in_one_or_the_other == "one" :
				individual_out_pileups = os.path.join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(the_other_locus_id, marker_pos_the_other_locus))
				Pileups_Generator("%s:%d-%d" %(the_other_locus_id, marker_pos_the_other_locus-25, marker_pos_the_other_locus+25), assembled_genome_file, in_bam_file, "individual", individual_out_pileups)
				genotype_the_other_locus, cov_the_other_locus, type_the_other_locus = Pileups_Analyzer(marker_pos_the_other_locus, individual_out_pileups)
				if the_other_locus_strand == "-" :
					if genotype_the_other_locus != "./." and genotype_the_other_locus != "NA/NA" :
						genotype_the_other_locus = RC(genotype_the_other_locus)
				individual_out_pileups = os.path.join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(one_locus_id, marker_pos_one_locus))
				genotype_one_locus, cov_one_locus, type_one_locus = Individual_Genotype_Analyzer(tmp_genoytpes[j], ref_allele_one_locus, alt_allele_one_locus)
				if genotype_one_locus == genotype_the_other_locus and cov_one_locus + cov_the_other_locus <= 7 and cov_one_locus + cov_the_other_locus > 0 :
					print "\t\t***",samples[j], genotype_one_locus, cov_one_locus, genotype_the_other_locus, cov_the_other_locus, ":D"
					num_individual_with_expected_genotype_pattern += 1
					if type_one_locus == "ref" :
						proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("0/0:%d" %(cov_the_other_locus+cov_one_locus))
					else :
						proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("1/1:%d" %(cov_the_other_locus+cov_one_locus))
				elif genotype_the_other_locus == "./." and not genotype_one_locus in ["./.", "NA/NA"] and cov_the_other_locus == 0 and cov_one_locus <= 8 and cov_one_locus > 0 :
					print "\t\t***",samples[j], genotype_one_locus, cov_one_locus, genotype_the_other_locus, cov_the_other_locus, ":D"
					num_individual_with_expected_genotype_pattern += 1
					if type_one_locus == "ref" :
						proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("0/0:%d" %(cov_one_locus))
					else :
						proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("1/1:%d" %(cov_one_locus))
				elif genotype_one_locus == "./." and not genotype_the_other_locus in ["./.", "NA/NA"] and cov_one_locus == 0 and cov_the_other_locus <= 8 and cov_the_other_locus > 0 :
					print "\t\t***",samples[j],  genotype_one_locus, cov_one_locus, genotype_the_other_locus, cov_the_other_locus, ":D"
					num_individual_with_expected_genotype_pattern += 1
					if type_one_locus == "ref" :
						proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("0/0:%d" %(cov_the_other_locus))
					else :
						proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("1/1:%d" %(cov_the_other_locus))
				elif genotype_the_other_locus == "NA/NA" or genotype_one_locus == "NA/NA" or (genotype_one_locus == "./." and genotype_the_other_locus == "./."):
					print "\t\t***",samples[j], genotype_one_locus, cov_one_locus, genotype_the_other_locus, cov_the_other_locus, "what the hell"
					tmp_num_individual -= 1
					proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("./.:0")
				else :
					print "\t\t***disagree",samples[j], genotype_one_locus, cov_one_locus, genotype_the_other_locus, cov_the_other_locus
					proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("./.:0")
			else :
				individual_out_pileups = os.path.join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(one_locus_id, marker_pos_one_locus))
				Pileups_Generator("%s:%d-%d" %(one_locus_id, marker_pos_one_locus-25, marker_pos_one_locus+25), assembled_genome_file, in_bam_file, "individual", individual_out_pileups)
				genotype_one_locus, cov_one_locus, type_one_locus = Pileups_Analyzer(marker_pos_one_locus, individual_out_pileups)
				if the_other_locus_strand == "-" :
					if genotype_one_locus != "./." and genotype_one_locus != "NA/NA" :
						genotype_one_locus = RC(genotype_one_locus)
				individual_out_pileups = os.path.join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(the_other_locus_id, marker_pos_the_other_locus))
				genotype_the_other_locus, cov_the_other_locus, type_the_other_locus = Individual_Genotype_Analyzer(tmp_genoytpes[j], ref_allele_the_other_locus, alt_allele_the_other_locus)
				if genotype_one_locus == genotype_the_other_locus and cov_one_locus + cov_the_other_locus <= 7 and cov_one_locus + cov_the_other_locus > 0 :
					print "\t\t***",samples[j], genotype_the_other_locus, cov_the_other_locus, genotype_one_locus, cov_one_locus, ":D"
					num_individual_with_expected_genotype_pattern += 1
					if type_the_other_locus == "ref" :
						proposed_genotype_at_variant_site[the_other_locus_id+':'+str(marker_pos_the_other_locus)].append("0/0:%d" %(cov_the_other_locus+cov_one_locus))
					else :
						proposed_genotype_at_variant_site[the_other_locus_id+':'+str(marker_pos_the_other_locus)].append("1/1:%d" %(cov_the_other_locus+cov_one_locus))
				elif genotype_the_other_locus == "./." and (not genotype_one_locus in ["./.", "NA/NA"]) and cov_the_other_locus == 0 and cov_one_locus <= 8 and cov_one_locus > 0 :
					print "\t\t***",samples[j], genotype_the_other_locus, cov_the_other_locus, genotype_one_locus, cov_one_locus, ":D"
					num_individual_with_expected_genotype_pattern += 1
					if type_the_other_locus == "ref" :
						proposed_genotype_at_variant_site[the_other_locus_id+':'+str(marker_pos_the_other_locus)].append("0/0:%d" %(cov_one_locus))
					else :
						proposed_genotype_at_variant_site[the_other_locus_id+':'+str(marker_pos_the_other_locus)].append("1/1:%d" %(cov_one_locus))
				elif genotype_one_locus == "./." and (not genotype_the_other_locus in ["./.", "NA/NA"]) and cov_one_locus == 0 and cov_the_other_locus <= 8 and cov_the_other_locus > 0 :
					print "\t\t***",samples[j], genotype_the_other_locus, cov_the_other_locus, genotype_one_locus, cov_one_locus, ":D"
					num_individual_with_expected_genotype_pattern += 1
					if type_the_other_locus == "ref" :
						proposed_genotype_at_variant_site[the_other_locus_id+':'+str(marker_pos_the_other_locus)].append("0/0:%d" %(cov_the_other_locus))
					else :
						proposed_genotype_at_variant_site[the_other_locus_id+':'+str(marker_pos_the_other_locus)].append("1/1:%d" %(cov_the_other_locus))
				elif genotype_the_other_locus == "NA/NA" or genotype_one_locus == "NA/NA" or (genotype_one_locus == "./." and genotype_the_other_locus == "./."):
					print "\t\t***",samples[j], genotype_the_other_locus, cov_the_other_locus, genotype_one_locus, cov_one_locus, "what the hell"
					tmp_num_individual -= 1
					proposed_genotype_at_variant_site[the_other_locus_id+':'+str(marker_pos_the_other_locus)].append("./.:0")
				else :
					print "\t\t***disagree",samples[j], genotype_the_other_locus, cov_the_other_locus, genotype_one_locus, cov_one_locus
					proposed_genotype_at_variant_site[the_other_locus_id+':'+str(marker_pos_the_other_locus)].append("./.:0")
		print "\tnumber of samples having expected genotypes at the current pair of markers:", num_individual_with_expected_genotype_pattern, tmp_num_individual
		if tmp_num_individual > len(samples)/2.0 :
			if num_individual_with_expected_genotype_pattern/float(tmp_num_individual) > 0.7 :
				num_variant_pair_with_expected_gt_pattern_between_loci += 1
	if len(variant_pairs) == 1 :
		if num_variant_pair_with_expected_gt_pattern_between_loci == 1 :
			#print "\tproposed genotype:", proposed_genotype_at_variant_site
			return proposed_genotype_at_variant_site
	else :
		if num_variant_pair_with_expected_gt_pattern_between_loci > len(variant_pairs)-num_variant_pair_with_expected_gt_pattern_between_loci :
			#print "\tproposed genotype:", proposed_genotype_at_variant_site
			return proposed_genotype_at_variant_site

# given an individual genotype as defined in the INFO field in the VCF file, the function spit out the genotype and coverage at the site #
def Individual_Genotype_Analyzer(individual_genotype, ref_allele, alt_allele) :
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

# given a marker position and a pileup file, spit out the genotype at the site #
def Pileups_Analyzer(marker_pos, individual_out_pileups) :
	genotype = "./."
	cov_at_marker_pos = 0
	type = "NA"
	num_ref_allele, num_alt_allele = 0, 0
	fPILEUP = open(individual_out_pileups, 'r')
	for line in fPILEUP :
		tmp_line = re.split('\t', line.strip())
		if marker_pos == int(tmp_line[1]) :
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
				#if num_alt_allele > num_ref_allele :
				#	genotype = alt_allele + alt_allele
				#	cov_at_marker_pos = num_alt_allele
				#	type = "ref"
				#elif num_ref_allele > num_alt_allele :
				#	genotype = tmp_line[2] + tmp_line[2]
				#	cov_at_marker_pos = num_ref_allele
				#	type = "alt"
				#else :
				genotype = "NA/NA"
				cov_at_marker_pos = -1
	fPILEUP.close()
	return genotype, cov_at_marker_pos, type

# reverse and complementary of the given sequence #
def RC(seq) :
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

# the assembly that you'd like to identify falsely splitted loci #
def Read_Ref_Seq(assembled_genome_file) :
	header, seq = "", ""
	assembled_genome_dict = {}
	fREF = open(assembled_genome_file, 'r')
	for line in fREF :
		if line.startswith('>') :
			if seq != "" :
				assembled_genome_dict[header] = seq
				seq = ""
			header = line.strip()[1:]
		else :
			seq += line.strip().upper()
	if seq != "" :
		assembled_genome_dict[header] = seq
		seq = ""
	fREF.close()
	return assembled_genome_dict

def Main(alignment_file, assembled_genome_file, vcf_file, all_sample_bam_file, splitted_sample_bam_dir, out_dir) :
	variants_dict, samples, smt_genotype_dict = VCF_Parser(vcf_file)
	assembled_genome_dict = Read_Ref_Seq(assembled_genome_file)
	split_loci_candidates, cigar_dict = Get_Splitted_loci_Candidates(alignment_file, assembled_genome_file, all_sample_bam_file, out_dir, splitted_sample_bam_dir, variants_dict, samples, smt_genotype_dict, assembled_genome_dict)

if __name__ == "__main__" :
	assembled_genome_file = sys.argv[1]
	alignment_file = sys.argv[2]
	vcf_file = sys.argv[3]
	splitted_sample_bam_dir = sys.argv[4]					# directory of each sample's bam file #
	all_sample_bam_file = sys.argv[5]						# single bam file with all the sample included #
	out_dir = sys.argv[6]									# the base output directory under which multiple sub-directory will be created #

	###############################
	# check the file availability #
	###############################
	if not os.path.exists(alignment_file) :
		print "Error: cannot find your alignment file %s, please check your path" %(alignment_file)
		sys.exit(1)
	if not os.path.exists(assembled_genome_file) :
		print "Error: cannot find your assembly sequence file, please check your path" %(assembled_genome_file)
		sys.exit(1)
	if not os.path.exists(vcf_file) :
		print "Error: cannot find  your genotypes file (vcf), please check your path" %(vcf_file)
		sys.exit(1)
	if not os.path.exists(splitted_sample_bam_dir) :
		print "Error: cannot find your directory %s with each sample's bam file" %(splitted_sample_bam_dir)
		sys.exit(1)

	Main(alignment_file, assembled_genome_file, vcf_file, all_sample_bam_file, splitted_sample_bam_dir, out_dir)
