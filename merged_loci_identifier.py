# the script is used to deal with alleles that were merged during the assembly #
# 1. identify a potential set of sites based on genotypes and coverage depth #
# 2. make target windows around each of the sites #
# 3. map reads from each individual on each of the targeted regions #
# 4. grap reads such that one end map to the targeted region and the other end map to some places else, hopefully to get two cluseters of reads #
# 5. each of the reads cluster indicate position of each of the alleles. In this step, we will figure out which allele goes where #

# update: apply a variance filter so that for a locus candidate, the abnormal coverage must be uniform across the locus #

import re, shlex
from argparse import ArgumentParser
from os import makedirs
from os.path import exists, join, dirname, realpath
from sys import exit, stdout, stderr, argv
from datetime import datetime
from subprocess import Popen, PIPE

# read the assembly to get the assembled genome #
def read_assembly(assembly_file) :
	stdout.write("%s reading sequences from assemly file: %s\n" %(timestamper(), assembly_file))
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
	stdout.write("%s %d sequences parsed\n" %(timestamper(), numseq))
	return assembled_genome_dict

# generate pileups for a given region #
def generate_pileups(bam_file, assembly_file, region, out_pileup) :
	pileup_cmd = "samtools mpileup -d 9999 -f %s -r %s %s" %(assembly_file, region, bam_file)
	stdout.write(timestamper() + " " + pileup_cmd + "\n")
	proc = Popen(shlex.split(pileup_cmd), stdout=file(out_pileup, 'w'), stderr=PIPE, shell=False)
	proc_stderr = proc.communicate()[1]
	if proc.returncode != 0 :
		stderr.write(timestamper() + "[Pileups Generating Error] %s\n" %(proc_stderr))
		exit()

def identify_collapsed_loci(pileup_file, scaffold, assembled_genome_dict, heterozygous_sites, seq_outdir, hist_outdir, plot_outdir, min_cov_threshold, max_cov_threshold, min_alignment_len) :

	stdout.write(timestamper() + " %s:" %(scaffold))
	# get the coverage of a scaffold from pileup #
	covs = read_pileups(pileup_file)

	# scan for all possible region #
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
		median_cov_locus = calcualte_median_cov(tmp_candidate_cov)
	#print "\tmerged loci candidates before further merging", merged_loci_candidates

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
					#print "\t", scaffold, "further merged", cur_locus_start, pre_locus_end
					pre_locus_start = pre_locus_start
					pre_locus_end = cur_locus_end
				else :
					median_cov = calcualte_median_cov(covs[pre_locus_start:pre_locus_end])
					if pre_locus_end - pre_locus_start >= min_alignment_len and median_cov >= min_cov_threshold and median_cov <= max_cov_threshold :
						#print "\tcandidates:", scaffold, pre_locus_start, pre_locus_end, pre_locus_end-pre_locus_start, median_cov
						if apply_genotype_filter(scaffold, pre_locus_start, pre_locus_end, heterozygous_sites) :
							updated_merged_loci_candidates.append("%d-%d" %(pre_locus_start, pre_locus_end))
							fOUT = open(join(seq_outdir, "%s_%d_%d_%dx.fasta" %(scaffold, pre_locus_start, pre_locus_end, median_cov)), 'w')
							fOUT.write(">%s_%d_%d_%dx\n%s\n" %(scaffold, pre_locus_start, pre_locus_end, median_cov, assembled_genome_dict[scaffold][pre_locus_start:pre_locus_end]))
							fOUT.close()
							generate_histogram(covs[pre_locus_start:pre_locus_end], scaffold, pre_locus_start, pre_locus_end, median_cov, hist_outdir, plot_outdir)
					pre_locus_start = cur_locus_start
					pre_locus_end = cur_locus_end
		median_cov = calcualte_median_cov(covs[pre_locus_start:pre_locus_end])
		if pre_locus_end - pre_locus_start >= min_alignment_len and median_cov >= min_cov_threshold and median_cov <= max_cov_threshold :
			#print "\tcandidates:", scaffold, pre_locus_start, pre_locus_end, pre_locus_end-pre_locus_start, median_cov
			if apply_genotype_filter(scaffold, pre_locus_start, pre_locus_end, heterozygous_sites) :
				updated_merged_loci_candidates.append("%d-%d" %(pre_locus_start, pre_locus_end))
				fOUT = open(join(seq_outdir, "%s_%d_%d_%dx.fasta" %(scaffold, pre_locus_start, pre_locus_end, median_cov)), 'w')
				fOUT.write(">%s_%d_%d_%dx\n%s\n" %(scaffold, pre_locus_start, pre_locus_end, median_cov, assembled_genome_dict[scaffold][pre_locus_start:pre_locus_end]))
				fOUT.close()
				generate_histogram(covs[pre_locus_start:pre_locus_end], scaffold, pre_locus_start, pre_locus_end, median_cov, hist_outdir, plot_outdir)
		#print "\t%s:" %(scaffold), updated_merged_loci_candidates

	if len(updated_merged_loci_candidates) == 0 :
		stdout.write("\tno collapsed loci candidates found\n")
	else :
		stdout.write("\t%d collapsed loci candidates found\n" %(len(updated_merged_loci_candidates)))
	return updated_merged_loci_candidates

# applying genotype filter to get merged loci candidates #
def apply_genotype_filter(scaffold, locus_start, locus_end, heterozygous_sites) :
	num_het_markers_within_merged_loci_candidates = 0
	if heterozygous_sites.has_key(scaffold) :
		for i in range(len(heterozygous_sites[scaffold])) :
			if int(heterozygous_sites[scaffold][i]) >= locus_start and int(heterozygous_sites[scaffold][i]) <= locus_end :
				num_het_markers_within_merged_loci_candidates += 1
				#print "\tmarkers:", scaffold, heterozygous_sites[scaffold][i]
		#print "\tWithin the candidate merged loci, number of markers that at least 90%% of individuals are heterozygous: %d" %(num_het_markers_within_merged_loci_candidates)
		if num_het_markers_within_merged_loci_candidates != 0 :
			return 1
		else :
			return 0
	else :
		return 0

# parsing the vcf file to get markers and associated genotypes #
def get_genotypes(genotypes_file, het_prop_per_site) :
	stdout.write(timestamper() + " Getting genotype information from genotype file %s\n" %(genotypes_file))
	heterozygous_sites = {}
	samples = []
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
				stdout.write(timestamper() + " %d variant sites being parsed\n" %(total_sites))
			num_het_at_a_site, num_missing_at_a_site = 0, 0
			for i in range(len(samples)) :
				ind_format_field = re.split(':', genotypes[i])
				ind_gt = ind_format_field[0]
				if ind_gt == "0/1" :
					num_het_at_a_site += 1
				elif ind_gt == "./." :
					num_missing_at_a_site += 1
			"""
			we only consider sites where either (i) at least 90% of individuals are heterzygous #
												(ii) 25% of the individuals at the site do not have genotypes available #
			"""
			if num_missing_at_a_site != len(samples) :
				if float(num_het_at_a_site)/(len(samples) - num_missing_at_a_site) >= het_prop_per_site or float(num_missing_at_a_site)/len(samples) >= 0.25*len(samples) :
					if not heterozygous_sites.has_key(marker_chr) :
						heterozygous_sites[marker_chr] = [marker_pos]
					else :
						heterozygous_sites[marker_chr].append(marker_pos)
					num_het_sites += 1
	if total_sites % 10000 != 0 :
		stdout.write(timestamper() + " %d sites being parsed\n" %(total_sites-(total_sites/10000)*10000))
	fVCF.close()
	stdout.write(timestamper() + " %d variant sites where at least 90%% individual are heterozygous\n" %(num_het_sites))
	return heterozygous_sites

# given a coverage array, this function calculate the median coverage #
def calcualte_median_cov(covs) :
	sorted_covs = sorted(covs)
	if len(sorted_covs)%2 == 0 :
		median_cov = (sorted_covs[len(sorted_covs)/2-1]+sorted_covs[len(sorted_covs)/2])/2
	else :
		median_cov = sorted_covs[len(sorted_covs)/2]
	return median_cov

# generate coverage histogram for each of the merged loci candidates #
def generate_histogram(locus_candidate_cov, scaffold, locus_start, locus_end, median_cov_locus, hist_outdir, plot_outdir) :
	out_histfile = join(hist_outdir, "%s_%d_%d_%dx.hist" %(scaffold, locus_start, locus_end, median_cov_locus))
	fHIST = open(out_histfile, 'w')
	for i in range(1, len(locus_candidate_cov)+1) :
		fHIST.write("%d\t%d\n" %(i, locus_candidate_cov[i-1]))
	fHIST.close()
	out_covplot_file = join(plot_outdir, "%s_%d_%d_%dx.cov_plot.pdf" %(scaffold, locus_start, locus_end, median_cov_locus))
	#generate_cov_plot(out_histfile, out_covplot_file)

# generate coverage plot for each of the merged loci candidates #
def generate_cov_plot(out_histfile, out_covplot_file) :
	Rscript = "/N/u/simozhan/Mason/project/rpgc/worm/src/cov_plot.r"
	if not os.path.exists(Rscript) :
		print "Error: cannot find your R script for generating the coverage plot: %s" %(Rscript)
		sys.exit(1)
	Rcmd = "R CMD BATCH --vanilla \"--args %s %s\" %s test" %(out_histfile, out_covplot_file, Rscript)
	os.system(Rcmd)

# reading the pileups and get the coverage #
def read_pileups(pileup_file) :
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
				#print "\tGap: ", i, last_pos, cur_pos
				covs.append(int(tmp_line[3]))
			else :
				covs.append(int(tmp_line[3]))
			last_pos = cur_pos
	fPILEUP.close()
	return covs

def make_dirs_if_necessary(*dirs) :
	""" create directory if necessary """
	for dir in dirs :
		if not exists(dir) :
			makedirs(dir)

def timestamper() :
	""" generate a dash_separated time stamped string """
	return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

def check_file_existence(*files) :
	""" check the existence of file/files """
	for file in files :
		if not exists(file) :
			stderr.write(timestamper() + " [IO Error] Cannot find the file you provided: %s\n" %(file))
			exit()

def main(bam_file, pileup_off, outdir, assembly_file, genotypes_file, min_alignment_len, min_cov_threshold, max_cov_threshold, het_prop_per_site) :
	"""
	preparing subfolders under the output directory
	"""
	stdout.write("\n%s Setting up output directory: %s\n" %(timestamper(), outdir))
	"""
	pileups_dir: directory where pileups of scaffolds will be put
	seq_outdir: directory where sequences of candidates of merged loci will be put
	hist_outdir: directory where coverage histograms of candidates of merged loci will be put
	plot_outdir: directory where plots of candidates of merged loci will be put
	summary_dir: directory where the summary file will be created
	"""
	pileups_dir = join(outdir, "pileups")
	summary_dir = join(outdir, "cov%d_%d_hs%d_len%d/summary" %(min_cov_threshold, max_cov_threshold, het_prop_per_site, min_alignment_len))
	seq_outdir = join(outdir, "cov%d_%d_hs%d_len%d/seqs" %(min_cov_threshold, max_cov_threshold, het_prop_per_site, min_alignment_len))
	hist_outdir = join(outdir, "cov%d_%d_hs%d_len%d/hists" %(min_cov_threshold, max_cov_threshold, het_prop_per_site, min_alignment_len))
	""" will fix in the future version """
	#plot_outdir = join(outdir, "cov%d_%d_hs%d_len%d/plots" %(min_cov_threshold, max_cov_threshold, het_prop_per_site, min_alignment_len))
	make_dirs_if_necessary(pileups_dir, seq_outdir, hist_outdir, plot_outdir, summary_dir)				# #

	# read the assembly #
	assembled_genome_dict = read_assembly(assembly_file)

	# generate pileups for each scaffold if necessary #
	if pileup_off :
		for scaffold in assembled_genome_dict.iterkeys() :
			out_scaffold_pileup = join(pileups_dir, scaffold+".pileups")
			stdout.write(timestamper() + " generating pileups for %s\n" %(scaffold))
			generate_pileups(bam_file, assembly_file, scaffold, out_scaffold_pileup)

	# get sites with at least X% of individual in the mapping population are heterozygous #
	heterozygous_sites = get_genotypes(genotypes_file, het_prop_per_site)

	# detect merged loci candidates for each of the scaffolds #
	stdout.write("\n" + timestamper() + " Start idnetifying collapsed loci candidates by scaffold\n")
	merged_loci_candidates_dict = {}
	for scaffold in assembled_genome_dict.iterkeys() :
		out_scaffold_pileup = join(pileups_dir, scaffold+".pileups")
		if not merged_loci_candidates_dict.has_key(scaffold) :
			merged_loci_candidates_dict[scaffold] = identify_collapsed_loci(out_scaffold_pileup, scaffold, assembled_genome_dict, heterozygous_sites, seq_outdir, hist_outdir, plot_outdir, min_cov_threshold, max_cov_threshold, min_alignment_len)
			if len(merged_loci_candidates_dict[scaffold]) == 0:
				del merged_loci_candidates_dict[scaffold]

	summary_file = join(summary_dir, "merged_loci_candidates.summary")
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
	stdout.write("\n" + timestamper() + " %d number of scaffolds have merged loci candidates identified\n" %(num_scaffold_with_merged_loci_candidates))
	stdout.write("\n" + timestamper() + " %d number of merged loci candidates being identified using coverage analysis\n" %(num_merged_loci_candidates))

	""" debug purpose """
	#out_scaffold_pileup = os.path.join(pileups_dir, "scaffold_172.pileups")
	#out_scaffold_hist = os.path.join(pileups_dir, "scaffold_172.hist")
	#merged_loci_candidates_dict["scaffold_172"] = identify_collapsed_loci(out_scaffold_pileup, "scaffold_172", assembled_genome_dict)
	#print len(assembled_genome_dict["scaffold_172"])

if __name__ == "__main__" :

	parser = ArgumentParser(description="Identifying loci/alleles that are erroneously collapsed by assemblers")
	parser.add_argument("-bam", metavar="FILE", dest = "bam_file", help="BAM file of the mapping population (including all the individuals)")
	parser.add_argument("-outdir", metavar="DIR", dest="outdir", help="output directory under which subfolder will be created")
	parser.add_argument("-assembly", metavar="FILE", dest="assembly_file", help="the assembled genome in fasta file")
	parser.add_argument("-gt", metavar="FILE", dest="genotypes_file", help="specify the VCF file of genotypes of the mapping population")
	parser.add_argument("-pileupoff", action="store_true", dest="pileup_off", help="turn off generating pileups")
	parser.add_argument("-min_len", metavar="INT", dest="min_alignment_len", type=int, default=1000, help="the minimum length of loci to be considered. DEFAULT=1000")
	parser.add_argument("-het", metavar="INT", dest="het_prop_per_site", type=int, default=90, help="the minimum proportion of individuals in the mapping population that are heterozygous at a variant site. DEFAULT=90")
	parser.add_argument("-min_cov", metavar="INT", dest="min_cov_threshold", type=int, default=650, help="the minimum coverage of loci to be considered. DEFAULT=650")
	parser.add_argument("-max_cov", metavar="INT", dest="max_cov_threshold", type=int, default=1100, help="the maximum coverage of loci to be considered. DEFAULT=1100")
	args = parser.parse_args()

	check_file_existence(args.bam_file, args.genotypes_file, args.assembly_file)
	make_dirs_if_necessary(args.outdir)

	stdout.write("\n#python %s" %(argv[0]))
	for k, v in vars(args).iteritems() :
		stdout.write(" -%s %s" %(k, v))
	stdout.write("\n")
	stdout.write("#BAM file: \t\t%s\n" %(args.bam_file))
	stdout.write("#assembled genome: \t%s\n" %(args.assembly_file))
	stdout.write("#genotypes file: \t%s\n" %(args.genotypes_file))
	stdout.write("#output directory: \t%s\n" %(args.outdir))
	if not args.pileup_off :
		stdout.write("#generate pileups: \tTrue\n")
	else :
		stdout.write("#generate pileups: \tFalse\n")
	stdout.write("#min len of loci: \t%d\n" %(args.min_alignment_len))
	stdout.write("#min cov of loci: \t%d\n" %(args.min_cov_threshold))
	stdout.write("#max cov of loci: \t%d\n" %(args.max_cov_threshold))
	stdout.write("#min proportion of het indvs: \t%.2f\n" %(args.het_prop_per_site/100.0))

	main(args.bam_file, not args.pileup_off, args.outdir, args.assembly_file, args.genotypes_file, args.min_alignment_len, args.min_cov_threshold, args.max_cov_threshold, args.het_prop_per_site/100.0)
