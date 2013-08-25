import re
from argparse import ArgumentParser
from os import makedirs, stat, system
from os.path import join, exists, dirname, realpath
from sys import exit, stdout, stderr
from datetime import datetime

def identify_splited_loci(alignment_file, assembled_genome_file, all_sample_bam_file, out_dir, splitted_sample_bam_dir, variants_dict, samples, smt_genotype_dict, assembled_genome_dict, min_aln_len, min_idn) :
	"""
	get a set of pairs of loci to be candidates of falsely splitted loci, given an alignment of an assembly against itself
	"""

	stdout.write(timestamper() + " Identifying splitted loci candidates\n")
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
	summary_dir = join(out_dir, "summary")
	make_dirs_if_needed(summary_dir)
	proposed_genotypes_outfile = join(summary_dir, "non_tandem_splitted_loci.proposed_genotypes")
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
					num_loci_passed_len_idn_filter += 1
					stdout.write("> %s %s\n" %(one_locus, the_other_locus))
					"""
					apply coverage filter
					"""
					proposed_genotype_at_variant_site = {}
					tmp_one_locus = re.sub(':', '-', re.sub('_', '', one_locus))
					one_locus_pileups_out = join(loci_wide_pileups_dir, tmp_one_locus+".pileups")
					tmp_the_other_locus = re.sub(':', '-', re.sub('_', '', the_other_locus))
					the_other_locus_pileups_out = join(loci_wide_pileups_dir, tmp_the_other_locus+".pileups")
					pass_or_filter = apply_cov_filter(one_locus, the_other_locus, assembled_genome_file, all_sample_bam_file, one_locus_pileups_out, the_other_locus_pileups_out)
					if pass_or_filter :
						"""
						localizing markers in each pair of loci
						"""
						variants_in_one_locus = localize_markers_in_loci(t_id, t_start, t_end, t_strand, variants_dict, smt_genotype_dict)
						variants_in_the_other_locus = localize_markers_in_loci(q_id, q_start, q_end, q_strand, variants_dict, smt_genotype_dict)

						"""
						find each pair of markers in each of loci
						"""
						variant_pairs = []
						if len(variants_in_one_locus) > 0 :
							one_or_the_other = 0
							variant_pairs += get_variant_pair_in_splitted_loci_Candidates(t_id, t_start, t_end, t_strand, q_id, q_start, q_end, q_strand, one_or_the_other, variants_in_one_locus, cigar, variants_dict, assembled_genome_dict)
						if len(variants_in_the_other_locus) > 0 :
							one_or_the_other = 1
							variant_pairs += get_variant_pair_in_splitted_loci_Candidates(q_id, q_start, q_end, q_strand, t_id, t_start, t_end, t_strand, one_or_the_other, variants_in_the_other_locus, cigar, variants_dict, assembled_genome_dict)

						"""
						 Creating a database of variants that are identified in each pair of loci candidates #
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

						num_localized_variants += len(variants_in_one_locus) + len(variants_in_the_other_locus)
						if len(variants_in_one_locus) > 0 or len(variants_in_the_other_locus) > 0 :
							fSEQ_LOCI_MARKER.write(">%s_%d_%d_1\n%s\n" %(t_id, t_start, t_end, assembled_genome_dict[t_id][t_start-1:t_end]))
							if q_strand == '+' :
								fSEQ_LOCI_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, assembled_genome_dict[q_id][q_start-1:q_end]))
							else :
								fSEQ_LOCI_MARKER.write(">%s_%d_%d_2\n%s\n" %(q_id, q_start, q_end, rc(assembled_genome_dict[q_id][q_end-1:q_start])))
							loci_candidate_out_subdir = join(markers_wide_pileups_dir, "%s_%d_%d-%s_%d_%d" %(t_id, t_start, t_end, q_id, q_start, q_end))
							if not exists(loci_candidate_out_subdir) :
								make_dirs_if_needed(loci_candidate_out_subdir)
							proposed_genotype_at_variant_site = investigate_genotypes(variant_pairs, variants_dict, samples, loci_candidate_out_subdir, splitted_sample_bam_dir, assembled_genome_file)
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
							stdout.write("\tnumber of pairs of variants can be localized: %d\n" %(num_localized_variants))
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
						stdout.write("\tfailed the coverage filter: %s %s\n" %(one_locus, the_other_locus))
						system("rm %s" %(one_locus_pileups_out))
						system("rm %s" %(the_other_locus_pileups_out))
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

	print ""
	print "predicted splitted loci", predicted_falsely_splitted_loci
	print "variants database", variants_in_loci_candidates_db
	out_summary_file = join(summary_dir, "splitted_loci_identification.summary")
	fSummary = open(out_summary_file, 'w')
	fSummary.write("total number of loci in the alignment: %d\n" %(total_num_loci))
	fSummary.write("number of loci alignment having at least 1000 bp long and 90%% identity: %d\n" %(num_loci_passed_len_idn_filter))
	fSummary.write("number of loci alignment who passed the length and identity filter failing to pass the coverage filter: %d\n" %(num_loci_failed_cov_filter))
	fSummary.write("number of loci alignment who passed all three filters having no localized variants: %d\n" %(num_loci_with_no_localized_variants))
	fSummary.write("number of pairs of loci being identified as falsely splitted loci: %d\n" %(num_pairs_predicted_falsely_splitted_loci))
	fSummary.write("among these predicted splitted loci, %d of them are non tandem ones, while %d are tandems\n" %(num_pairs_predicted_falsely_non_tandem_splitted_loci, num_pairs_predicted_falsely_splitted_loci-num_pairs_predicted_falsely_non_tandem_splitted_loci))
	fSummary.write("total number of pairs of markers being localized in falsely splitted loci: %d\n" %(num_localized_variants))

	return potential_dup_loci, cigar_dict

def generate_pileups(region, assembled_genome_file, bam_file, individual_sample_or_all, pileups_outfile) :
	if individual_sample_or_all == "all" :
		pileups_cmd = "samtools mpileup -d 9999 -f %s -r %s %s > %s" %(assembled_genome_file, region, bam_file, pileups_outfile)
	else :
		pileups_cmd = "samtools mpileup -f %s -r %s %s > %s" %(assembled_genome_file, region, bam_file, pileups_outfile)
	system(pileups_cmd)

# Given a pileup file of a region, the function gives you the median coverage #
def get_regional_median_coverage(pileups_file) :
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

def apply_cov_filter(one_locus, the_other_locus, assembled_genome_file, all_sample_bam_file, one_locus_pileups_out, the_other_locus_pileups_out) :
	generate_pileups(one_locus, assembled_genome_file, all_sample_bam_file, "all", one_locus_pileups_out)
	generate_pileups(the_other_locus, assembled_genome_file, all_sample_bam_file, "all", the_other_locus_pileups_out)
	if stat(one_locus_pileups_out).st_size == 0 and os.stat(the_other_locus_pileups_out).st_size == 0 :
		stdout.write("no pileups generated for both loci: %s %s\n" %(one_locus_pileups_out, the_other_locus_pileups_out))
		system("rm %s" %(one_locus_pileups_out))
		system("rm %s" %(the_other_locus_pileups_out))
		return 0
	else :
		median_cov_one_locus = get_regional_median_coverage(one_locus_pileups_out)
		median_cov_the_other_locus = get_regional_median_coverage(the_other_locus_pileups_out)
		stdout.write("\t %s %s\n" %(median_cov_one_locus, median_cov_the_other_locus))
		""" hard choice, fix in the future """
		if median_cov_one_locus + median_cov_the_other_locus < 600 :
			return 1
		else :
			return 0

def get_genotypes(gatk_vcf_file, smt_genotypes_file) :
	""" get the genotypes from a given VCF file """
	stdout.write(timestamper() +  " Getting genotypes from VCf file %s\n" %(gatk_vcf_file))
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
				if num_het <= 4 :
					tmp_genotypes = '\t'.join(tmp_line[9:])
					num_variants += 1
					if not variants_dict.has_key(chrom_id+':'+chrom_pos) :
						variants_dict[chrom_id+':'+chrom_pos] = "%s\t%s\t%s" %(ref_base, alt_base, tmp_genotypes)
					else :				# there should not have any redundant calls #
						stderr.write(timestamper() + " redundant calls found: %s %s\n" %(chrom_id, chrom_pos))
						stderr.write(timestamper() + " please run fix_redundant_calls.py first to solve the problem\n")
						exit(1)
	fGATK.close()
	stdout.write(timestamper() + " genotypes of %d individuals at %d variant sites parsed\n" %(len(samples), num_variants))

	# read in genotypes by samtools #
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

def localize_markers_in_loci(locus_id, locus_start, locus_end, locus_strand, variants_dict, smt_genoytpe_dict) :
	"""
	Localize variants in the potentially splitted loci #
	find markers within each of the identified potentially splitted loci #
	for each pair of markers, define an interval whose pileup will be generated #
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
		stdout.write("\tcannot find markers in the locus\n")
	else :
		print "\t", variants_in_locus
	return variants_in_locus

def get_variant_pair_in_splitted_loci_Candidates(one_locus_id, one_locus_start, one_locus_end, one_locus_strand, the_other_locus_id, the_other_locus_start, the_other_locus_end, the_other_locus_strand, one_or_the_other, variants_in_locus, cigar, variants_dict, assembled_genome_dict) :
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
					if not rc(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1]) in re.split('\t', variants_dict[variants_in_locus[i]])[0:2] :
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
					print "\t>%s" %(the_other_locus_id+':'+str(variant_pos_the_other_locus)), "\t\t", rc(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus:variant_pos_the_other_locus+10]), rc((assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1])), rc(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-11:variant_pos_the_other_locus-1])
			else :
				variant_pairs.append("%s:%s:%s|%s:%s:%s" %(the_other_locus_id, variant_pos_the_other_locus, the_other_locus_strand, one_locus_id, variant_pos_one_locus, one_locus_strand))
				if one_locus_strand == '-' :
					print "\t>%s" %(the_other_locus_id+":"+str(variant_pos_the_other_locus)), rc(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus:variant_pos_the_other_locus+5]), rc(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1]), rc(assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-6:variant_pos_the_other_locus-1])
				else:
					print "\t>%s" %(the_other_locus_id+":"+str(variant_pos_the_other_locus)), assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-6:variant_pos_the_other_locus-1], assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus-1], assembled_genome_dict[the_other_locus_id][variant_pos_the_other_locus:variant_pos_the_other_locus+5]
				print "\t>%s" %(one_locus_id+":"+str(variant_pos_one_locus)), assembled_genome_dict[one_locus_id][variant_pos_one_locus-6:variant_pos_one_locus-1], assembled_genome_dict[one_locus_id][variant_pos_one_locus-1], assembled_genome_dict[one_locus_id][variant_pos_one_locus:variant_pos_one_locus+5]
	return variant_pairs

def investigate_genotypes(variant_pairs, variants_dict, samples, loci_candidate_out_subdir, splitted_sample_bam_dir, assembled_genome_file) :
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
		variant_pairs_out_subdir = join(loci_candidate_out_subdir, "%s_%d-%s_%d" %(one_locus_id, marker_pos_one_locus, the_other_locus_id, marker_pos_the_other_locus))
		if not exists(variant_pairs_out_subdir) :
			makedirs(variant_pairs_out_subdir)
		num_individual_with_expected_genotype_pattern = 0
		tmp_num_individual = len(samples)
		for j in range(len(samples)) :
			individual_name = re.sub("rils", 'r', samples[j])
			in_bam_file = join(splitted_sample_bam_dir, individual_name+".bam")
			if variant_in_one_or_the_other == "one" :
				individual_out_pileups = join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(the_other_locus_id, marker_pos_the_other_locus))
				generate_pileups("%s:%d-%d" %(the_other_locus_id, marker_pos_the_other_locus-25, marker_pos_the_other_locus+25), assembled_genome_file, in_bam_file, "individual", individual_out_pileups)
				genotype_the_other_locus, cov_the_other_locus, type_the_other_locus = analyze_pileups(marker_pos_the_other_locus, individual_out_pileups)
				if the_other_locus_strand == "-" :
					if genotype_the_other_locus != "./." and genotype_the_other_locus != "NA/NA" :
						genotype_the_other_locus = rc(genotype_the_other_locus)
				individual_out_pileups = join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(one_locus_id, marker_pos_one_locus))
				genotype_one_locus, cov_one_locus, type_one_locus = analyze_individual_genotype(tmp_genoytpes[j], ref_allele_one_locus, alt_allele_one_locus)
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
					print "\t\t***",samples[j], genotype_one_locus, cov_one_locus, genotype_the_other_locus, cov_the_other_locus
					tmp_num_individual -= 1
					proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("./.:0")
				else :
					print "\t\t***disagree",samples[j], genotype_one_locus, cov_one_locus, genotype_the_other_locus, cov_the_other_locus
					proposed_genotype_at_variant_site[one_locus_id+':'+str(marker_pos_one_locus)].append("./.:0")
			else :
				individual_out_pileups = join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(one_locus_id, marker_pos_one_locus))
				generate_pileups("%s:%d-%d" %(one_locus_id, marker_pos_one_locus-25, marker_pos_one_locus+25), assembled_genome_file, in_bam_file, "individual", individual_out_pileups)
				genotype_one_locus, cov_one_locus, type_one_locus = analyze_pileups(marker_pos_one_locus, individual_out_pileups)
				if the_other_locus_strand == "-" :
					if genotype_one_locus != "./." and genotype_one_locus != "NA/NA" :
						genotype_one_locus = rc(genotype_one_locus)
				individual_out_pileups = join(variant_pairs_out_subdir, samples[j]+"_%s_%d.pileups" %(the_other_locus_id, marker_pos_the_other_locus))
				genotype_the_other_locus, cov_the_other_locus, type_the_other_locus = analyze_individual_genotype(tmp_genoytpes[j], ref_allele_the_other_locus, alt_allele_the_other_locus)
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
					print "\t\t***",samples[j], genotype_the_other_locus, cov_the_other_locus, genotype_one_locus, cov_one_locus
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
			print "\tproposed genotype:", proposed_genotype_at_variant_site
			return proposed_genotype_at_variant_site
	else :
		if num_variant_pair_with_expected_gt_pattern_between_loci > len(variant_pairs)-num_variant_pair_with_expected_gt_pattern_between_loci :
			print "\tproposed genotype:", proposed_genotype_at_variant_site
			return proposed_genotype_at_variant_site

# given an individual genotype as defined in the INFO field in the VCF file, the function spit out the genotype and coverage at the site #
def analyze_individual_genotype(individual_genotype, ref_allele, alt_allele) :
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
def analyze_pileups(marker_pos, individual_out_pileups) :
	genotype = "./."
	cov_at_marker_pos = 0
	type = "NA"
	num_ref_allele, num_alt_allele = 0, 0
	fPILEUP = open(individual_out_pileups, 'r')
	stderr.write("%s\n" %(individual_out_pileups))
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

def read_assembly(assembled_genome_file) :
	""" read the assembled genome """
	stdout.write(timestamper() + " Reading the assembled genome %s\n" %(assembled_genome_file))
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

def main(alignment_file, assembled_genome_file, gatk_vcf_file, smt_genotypes_file, all_sample_bam_file, splitted_sample_bam_dir, out_dir, min_aln_len, min_idn) :
	assembled_genome_dict = read_assembly(assembled_genome_file)
	variants_dict, samples, smt_genotype_dict = get_genotypes(gatk_vcf_file, smt_genotypes_file)
	split_loci_candidates, cigar_dict = identify_splited_loci(alignment_file, assembled_genome_file, all_sample_bam_file, out_dir, splitted_sample_bam_dir, variants_dict, samples, smt_genotype_dict, assembled_genome_dict, min_aln_len, min_idn/100.0)

if __name__ == "__main__" :
	parser = ArgumentParser(description="Identify alleles/loci that are erroneously split by assemblers")
	parser.add_argument("-assembly", metavar="FILE", dest="assembled_genome_file", required=True, help="Specify assembled genome in fasta")
	parser.add_argument("-aln", metavar="FILE", dest="alignment_file", required=True, help="Specify the file with alignment results of aligning assembled genome against itself. Note: lastz alignment results in CIGAR format is supported only in current version")
	parser.add_argument("-gatk_vcf", metavar="FILE", dest="gatk_vcf_file", required=True, help="Specify the VCF file with genotypes of the mapping population obtained from GATK")
	parser.add_argument("-bam_file", metavar="FILE", dest="all_sample_bam_file", required=True, help="Specify the BAM file of individuals of the entire mapping population")
	parser.add_argument("-bam_dir", metavar="DIR", dest="splitted_sample_bam_dir", required=True, help="Specify the directory with individual BAM files (not merged as a single BAM file like what -bam_file requires)")
	parser.add_argument("-smt_vcf", metavar="FILE", dest="smt_genotypes_file", required=True, help="Specify the VCF file with genotypes of the mapping population obtained from SAMtools")
	parser.add_argument("-out_dir", metavar="DIR", dest="out_dir", required=True, help="Specifying the output directory under which multiple subfolders will be created")
	parser.add_argument("-min_len", metavar="INT", dest="min_aln_len", type=int, default=1000, help="Specify the minimum alignment length of a pair of loci/alleles to be considered. DEFAULT: 1000")
	parser.add_argument("-min_idn", metavar="INT", dest="min_idn", type=int, default=90, help="Specify the minimum alignment identity (e.g. 90 mean 90%% identity) of a pair of loci/alleles to be considered. DEFAULT: 90")

	args = parser.parse_args()

	check_files_existence(args.alignment_file, args.assembled_genome_file, args.gatk_vcf_file, args.smt_genotypes_file)
	check_dirs_existence(args.splitted_sample_bam_dir)
	make_dirs_if_needed(args.out_dir)

	main(args.alignment_file, args.assembled_genome_file, args.gatk_vcf_file, args.smt_genotypes_file, args.all_sample_bam_file, args.splitted_sample_bam_dir, args.out_dir, args.min_aln_len, args.min_idn)
