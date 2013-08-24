# the script is used to deal with alleles that were merged during the assembly #
# 1. identify a potential set of sites based on genotypes and coverage depth #
# 2. make target windows around each of the sites #
# 3. map reads from each individual on each of the targeted regions #
# 4. grap reads such that one end map to the targeted region and the other end map to some places else, hopefully to get two cluseters of reads #
# 5. each of the reads cluster indicate position of each of the alleles. In this step, we will figure out which allele goes where #

# update: apply a variance filter so that for a locus candidate, the abnormal coverage must be uniform across the locus #

import os, sys
import re
import argparse

min_cov_threshold, max_cov_threshold = 650, 1100				# will become options when the scripts is finished #
heterozygous_rate_at_a_site = 90
min_alignment_len = 1000

# read the assembly to get the assembled genome #
def Read_Assembly(assembly_file) :
	print "parsing sequences from reference file: %s" %(assembly_file)
	numseq = 0
	seq, header = "", ""
	assembled_genome_dict = {}
	fASSEMBLY = open(assembly_file, 'r')
	for line in fASSEMBLY :
		if line.startswith('>') :
			if seq != "" :
				if not assembled_genome_dict.has_key(header) :
					assembled_genome_dict[header] = seq
					seq = ""
					numseq += 1
			header = line.strip()[1:]
		else :
			seq += line.strip()
	if seq != "" :
		if not assembled_genome_dict.has_key(header) :
			assembled_genome_dict[header] = seq
			numseq += 1
			seq = ""
	fASSEMBLY.close()
	print "%d number of sequences were parsed" %(numseq)
	return assembled_genome_dict

# generate pileups for a given region #
def Pileup_Generator(bam_file, assembly_file, region, out_pileup) :
	print "pileup generator"
	pileup_cmd = "samtools mpileup -d 9999 -f %s -r %s %s > %s" %(assembly_file, region, bam_file, out_pileup)
	print pileup_cmd
	os.system(pileup_cmd)

def Merged_loci_Detection_Cov_Analsysis(pileup_file, scaffold, assembled_genome_dict, heterozygous_sites, seq_outdir, hist_outdir, plot_outdir) :
	##############################################
	# get the coverage of a scaffold from pileup #
	##############################################
	covs = Pileup_Reader(pileup_file)

	################################
	# scan for all possible region #
	################################
	tmp_covs_per_win = []
	merged_loci_candidates = []
	i, win_size, step = 0, 100, 100
	tmp_locus_start, tmp_locus_end = -1, -1

	while i <= len(covs) - win_size :
		tmp_covs_per_win = sorted(covs[i:i+win_size])
		median_cov_per_win = (tmp_covs_per_win[win_size/2]+tmp_covs_per_win[(win_size/2)-1])/2
		if median_cov_per_win >=  550 :
			if tmp_locus_start == -1 :
				tmp_locus_start = i
			else :
				tmp_locus_end = i
		else :
			if tmp_locus_start != -1 and tmp_locus_end != -1 :
				tmp_candidate_cov = covs[tmp_locus_start:tmp_locus_end+win_size]
				merged_loci_candidates.append("%d-%d" %(tmp_locus_start, tmp_locus_end+win_size))
				median_cov_locus = (sorted(tmp_candidate_cov)[len(tmp_candidate_cov)/2]+sorted(tmp_candidate_cov)[len(tmp_candidate_cov)/2-1])/2
			tmp_candidate_cov = []
			tmp_locus_start, tmp_locus_end = -1, -1
		i += step
	if tmp_locus_start != -1 :
		if tmp_locus_end == -1 :
			tmp_candidate_cov = covs[tmp_locus_start:len(covs)]
			merged_loci_candidates.append("%d-%d" %(tmp_locus_start, len(covs)))
		else :
			tmp_candidate_cov = covs[tmp_locus_start:tmp_locus_end+win_size]
			merged_loci_candidates.append("%d-%d" %(tmp_locus_start, tmp_locus_end+win_size))
		median_cov_locus = Compute_Median_Cov(tmp_candidate_cov)
	print "\tmerged loci candidates before further merging", merged_loci_candidates

	# try to further merge the adjacent two merged loci candidates, if necessary #
	updated_merged_loci_candidates = []
	if len(merged_loci_candidates) > 0 :
		pre_locus_start, pre_locus_end = 0, 0
		for i in range(len(merged_loci_candidates)) :
			cur_locus_start = int(re.split('-', merged_loci_candidates[i])[0])
			cur_locus_end = int(re.split('-', merged_loci_candidates[i])[1])
			if pre_locus_start == 0 and pre_locus_end == 0 :
				pre_locus_start = cur_locus_start
				pre_locus_end = cur_locus_end
			else :
				if cur_locus_start - pre_locus_end <= 100 :
					print "\t", scaffold, "further merged", cur_locus_start, pre_locus_end
					pre_locus_start = pre_locus_start
					pre_locus_end = cur_locus_end
				else :
					median_cov = Compute_Median_Cov(covs[pre_locus_start:pre_locus_end])
					if pre_locus_end - pre_locus_start >= min_alignment_len and median_cov >= min_cov_threshold and median_cov <= max_cov_threshold :
						print "\tcandidates:", scaffold, pre_locus_start, pre_locus_end, pre_locus_end-pre_locus_start, median_cov
						if Apply_Genotype_Filter(scaffold, pre_locus_start, pre_locus_end, heterozygous_sites) :
							updated_merged_loci_candidates.append("%d-%d" %(pre_locus_start, pre_locus_end))
							fOUT = open(os.path.join(seq_outdir, "%s_%d_%d_%dx.fasta" %(scaffold, pre_locus_start, pre_locus_end, median_cov)), 'w')
							fOUT.write(">%s_%d_%d_%dx\n%s\n" %(scaffold, pre_locus_start, pre_locus_end, median_cov, assembled_genome_dict[scaffold][pre_locus_start:pre_locus_end]))
							fOUT.close()
							Histogram_Generator(covs[pre_locus_start:pre_locus_end], scaffold, pre_locus_start, pre_locus_end, median_cov, hist_outdir, plot_outdir)
					pre_locus_start = cur_locus_start
					pre_locus_end = cur_locus_end
		median_cov = Compute_Median_Cov(covs[pre_locus_start:pre_locus_end])
		if pre_locus_end - pre_locus_start >= min_alignment_len and median_cov >= min_cov_threshold and median_cov <= max_cov_threshold :
			print "\tcandidates:", scaffold, pre_locus_start, pre_locus_end, pre_locus_end-pre_locus_start, median_cov
			if Apply_Genotype_Filter(scaffold, pre_locus_start, pre_locus_end, heterozygous_sites) :
				updated_merged_loci_candidates.append("%d-%d" %(pre_locus_start, pre_locus_end))
				fOUT = open(os.path.join(seq_outdir, "%s_%d_%d_%dx.fasta" %(scaffold, pre_locus_start, pre_locus_end, median_cov)), 'w')
				fOUT.write(">%s_%d_%d_%dx\n%s\n" %(scaffold, pre_locus_start, pre_locus_end, median_cov, assembled_genome_dict[scaffold][pre_locus_start:pre_locus_end]))
				fOUT.close()
				Histogram_Generator(covs[pre_locus_start:pre_locus_end], scaffold, pre_locus_start, pre_locus_end, median_cov, hist_outdir, plot_outdir)
		print "\t%s:" %(scaffold), updated_merged_loci_candidates

	if len(updated_merged_loci_candidates) == 0 :
		print "\t*not found"
	else :
		print "\tFinally, number of merged loci_candidates: %d" %(len(updated_merged_loci_candidates))
	return updated_merged_loci_candidates

# applying genotype filter to get merged loci candidates #
def Apply_Genotype_Filter(scaffold, locus_start, locus_end, heterozygous_sites) :
	num_het_markers_within_merged_loci_candidates = 0
	if heterozygous_sites.has_key(scaffold) :
		for i in range(len(heterozygous_sites[scaffold])) :
			if int(heterozygous_sites[scaffold][i]) >= locus_start and int(heterozygous_sites[scaffold][i]) <= locus_end :
				num_het_markers_within_merged_loci_candidates += 1
				print "\tmarkers:", scaffold, heterozygous_sites[scaffold][i]
		print "\tWithin the candidate merged loci, number of markers that at least 90%% of individuals are heterozygous: %d" %(num_het_markers_within_merged_loci_candidates)
		if num_het_markers_within_merged_loci_candidates != 0 :
			return 1
		else :
			return 0
	else :
		return 0

# parsing the vcf file to get markers and associated genotypes #
def Genotypes_Parser(genotypes_file) :
	heterozygous_sites = {}
	num_het_sites, total_sites = 0, 0
	fVCF = open(genotypes_file, 'r')
	for line in fVCF :
		if line.startswith('#') :
			if line.startswith('#CHROM') :
				samples = re.split('\t', line.strip())[9:]
		else :
			tmp_line = re.split('\t', line.strip())
			genotypes = tmp_line[9:]
			marker_chr = tmp_line[0]					# chromosome where the marker sits #
			marker_pos = tmp_line[1]					# position on the chromosome #
			total_sites += 1
			if total_sites % 10000 == 0 :
				print "%d sites being parsed ... " %(total_sites)
			num_het_at_a_site, num_missing_at_a_site = 0, 0
			for i in range(len(samples)) :
				ind_format_field = re.split(':', genotypes[i])
				ind_gt = ind_format_field[0]
				if ind_gt == "0/1" :
					num_het_at_a_site += 1
				elif ind_gt == "./." :
					num_missing_at_a_site += 1
			# we only consider sites where either (i) at least 90% of individuals are heterzygous #
			#							   or (ii) 25% of the individuals at the site do not have genotypes available #
			if num_missing_at_a_site != len(samples) :
				if float(num_het_at_a_site)/(len(samples) - num_missing_at_a_site) >= 0.9 or float(num_missing_at_a_site)/len(samples) >= 11 :
					if not heterozygous_sites.has_key(marker_chr) :
						heterozygous_sites[marker_chr] = [marker_pos]
					else :
						heterozygous_sites[marker_chr].append(marker_pos)
					num_het_sites += 1
	if total_sites % 10000 != 0 :
		print "%d sites being parsed ... " %(total_sites-(total_sites/10000)*10000)
	fVCF.close()
	print "number of markers where at least 90%% individual are heterozygous: %d" %(num_het_sites)
	return heterozygous_sites

# given a coverage array, this function calculate the median coverage #
def Compute_Median_Cov(covs) :
	sorted_covs = sorted(covs)
	if len(sorted_covs)%2 == 0 :
		median_cov = (sorted_covs[len(sorted_covs)/2-1]+sorted_covs[len(sorted_covs)/2])/2
	else :
		median_cov = sorted_covs[len(sorted_covs)/2]
	return median_cov

# generate coverage histogram for each of the merged loci candidates #
def Histogram_Generator(locus_candidate_cov, scaffold, locus_start, locus_end, median_cov_locus, hist_outdir, plot_outdir) :
	out_histfile = os.path.join(hist_outdir, "%s_%d_%d_%dx.hist" %(scaffold, locus_start, locus_end, median_cov_locus))
	fHIST = open(out_histfile, 'w')
	for i in range(1, len(locus_candidate_cov)+1) :
		fHIST.write("%d\t%d\n" %(i, locus_candidate_cov[i-1]))
	fHIST.close()
	out_covplot_file = os.path.join(plot_outdir, "%s_%d_%d_%dx.cov_plot.pdf" %(scaffold, locus_start, locus_end, median_cov_locus))
	Coverage_Plot_Generator(out_histfile, out_covplot_file)

# generate coverage plot for each of the merged loci candidates #
def Coverage_Plot_Generator(out_histfile, out_covplot_file) :
	Rscript = "/N/u/simozhan/Mason/project/rpgc/worm/src/cov_plot.r"
	if not os.path.exists(Rscript) :
		print "Error: cannot find your R script for generating the coverage plot: %s" %(Rscript)
		sys.exit(1)
	Rcmd = "R CMD BATCH --vanilla \"--args %s %s\" %s test" %(out_histfile, out_covplot_file, Rscript)
	os.system(Rcmd)

# reading the pileups and get the coverage #
def Pileup_Reader(pileup_file) :
	fPILEUP = open(pileup_file, 'r')
	covs = []
	last_pos, cur_pos = 0, 0
	num_gap = 0
	for line in fPILEUP :
		tmp_line = re.split('\t', line.strip())
		cur_pos = int(tmp_line[1])
		if last_pos == 0 :
			last_pos = cur_pos
			covs.append(int(tmp_line[3]))
		else :
			if cur_pos - last_pos > 1 :
				num_gap += 1
				for i in range(cur_pos-last_pos-1) :
					covs.append(0)
				print "\tGap: ", i, last_pos, cur_pos
				covs.append(int(tmp_line[3]))
			else :
				covs.append(int(tmp_line[3]))
			last_pos = cur_pos
	fPILEUP.close()
	return covs

def Dir_Creator(dir) :
	if not os.path.exists(dir) :
		os.makedirs(dir)

def Main(all_sample_bam_file, pileups_or_not, root_dir, assembly_file, genotypes_file) :
	#############################################
	# preparing subfolders under rood directory #
	#############################################
	pileups_dir = os.path.join(root_dir, "pileups")
	summary_dir = os.path.join(root_dir, "cov%d_%d_hs%d_len%d/summary" %(min_cov_threshold, max_cov_threshold, heterozygous_rate_at_a_site, min_alignment_len))
	seq_outdir = os.path.join(root_dir, "cov%d_%d_hs%d_len%d/seqs" %(min_cov_threshold, max_cov_threshold, heterozygous_rate_at_a_site, min_alignment_len))
	hist_outdir = os.path.join(root_dir, "cov%d_%d_hs%d_len%d/hists" %(min_cov_threshold, max_cov_threshold, heterozygous_rate_at_a_site, min_alignment_len))
	plot_outdir = os.path.join(root_dir, "cov%d_%d_hs%d_len%d/plots" %(min_cov_threshold, max_cov_threshold, heterozygous_rate_at_a_site, min_alignment_len))
	Dir_Creator(pileups_dir)				# directory where pileups of scaffolds will be put #
	Dir_Creator(seq_outdir)						# directory where sequences of candidates of merged loci will be put #
	Dir_Creator(hist_outdir)				# directory where coverage histograms of candidates of merged loci will be put #
	Dir_Creator(plot_outdir)				# directory where plots of candidates of merged loci will be put #
	Dir_Creator(summary_dir)				# directory where the summary file will be created #

	# read the assembly #
	assembled_genome_dict = Read_Assembly(assembly_file)
	heterozygous_sites = Genotypes_Parser(genotypes_file)
	# generate pileups for each scaffold if necessary #
	if pileups_or_not :
		for scaffold in assembled_genome_dict.iterkeys() :
			out_scaffold_pileup = os.path.join(pileups_dir, scaffold+".pileups")
			Pileup_Generator(all_sample_bam_file, assembly_file, scaffold, out_scaffold_pileup)
	# detect merged loci candidates for each of the scaffolds #
	merged_loci_candidates_dict = {}
	for scaffold in assembled_genome_dict.iterkeys() :
		out_scaffold_pileup = os.path.join(pileups_dir, scaffold+".pileups")
		if not merged_loci_candidates_dict.has_key(scaffold) :
			print ">", scaffold
			merged_loci_candidates_dict[scaffold] = Merged_loci_Detection_Cov_Analsysis(out_scaffold_pileup, scaffold, assembled_genome_dict, heterozygous_sites, seq_outdir, hist_outdir, plot_outdir)
			if len(merged_loci_candidates_dict[scaffold]) == 0:
				del merged_loci_candidates_dict[scaffold]

	summary_file = os.path.join(summary_dir, "merged_loci_candidates.summary")
	fSUMMARY = open(summary_file, 'w')
	num_scaffold_with_merged_loci_candidates = 0
	num_merged_loci_candidates = 0
	for scaffold in merged_loci_candidates_dict.iterkeys() :
		num_scaffold_with_merged_loci_candidates += 1
		num_merged_loci_candidates += len(merged_loci_candidates_dict[scaffold])
		fSUMMARY.write("%s\t%d\t" %(scaffold, len(assembled_genome_dict[scaffold])))
		for i in range(len(merged_loci_candidates_dict[scaffold])) :
			if i == 0 :
				fSUMMARY.write("%s" %(merged_loci_candidates_dict[scaffold][i]))
			else :
				fSUMMARY.write(";%s" %(merged_loci_candidates_dict[scaffold][i]))
		fSUMMARY.write("\n")
		print scaffold, len(assembled_genome_dict[scaffold]), merged_loci_candidates_dict[scaffold]
	fSUMMARY.close()
	print "%d number of scaffolds have merged loci candidates identified" %(num_scaffold_with_merged_loci_candidates)
	print "%d number of merged loci candidates being identified using coverage analysis" %(num_merged_loci_candidates)

	# debug purpose #
	#out_scaffold_pileup = os.path.join(pileups_dir, "scaffold_172.pileups")
	#out_scaffold_hist = os.path.join(pileups_dir, "scaffold_172.hist")
	#merged_loci_candidates_dict["scaffold_172"] = Merged_loci_Detection_Cov_Analsysis(out_scaffold_pileup, "scaffold_172", assembled_genome_dict)
	#print len(assembled_genome_dict["scaffold_172"])

if __name__ == "__main__" :
	all_sample_bam_file = sys.argv[1]		# the bam file from which you will generate the pileups #
	root_dir = sys.argv[2]					# the root directory under which many subfolderd will be created correspondingly #
	pileups_or_not = int(sys.argv[3])		# a flag tells the program whether or not to generate pileup files #
	assembly_file = sys.argv[4]				# the assembled genome in fasta file #
	genotypes_file =  sys.argv[5]			# the vcf file of genotypes of the mapping population #

	if not os.path.exists(all_sample_bam_file) :
		print "Error: cannot find the pileup file you provided: %s" %(all_sample_bam_file)
		sys.exit(1)
	if not os.path.exists(genotypes_file) :
		print "Error: cannot find the vcf file with genotypes information you provided: %s" %(genotypes_file)
		sys.exit(1)

	if pileups_or_not :
		print "generate pileups: \tTrue"
	else :
		print "generate pileups: \tFalse"

	if not os.path.exists(root_dir) :
		os.makedirs(root_dir)
	print "root direcoty: \t%s" %(root_dir)

	if not os.path.exists(assembly_file) :
		print "Error: cannot find the assembly file: %s" %(assembly_file)
		sys.exit(1)
	print "assembled genome: \t%s" %(assembly_file)

	Main(all_sample_bam_file, pileups_or_not, root_dir, assembly_file, genotypes_file)
