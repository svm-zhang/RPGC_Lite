# The script is used to get the genotypes of each individual at each locus, and output to a file that is ready to be used by RQTL #
# At the same time, the program will identify candidate sites for incorrectly collapsed alleles (segments of DNA supposed to have two or more copies are incorrectly merged by assemblers)

import re
from argparse import ArgumentParser
from os.path import exists, dirname, realpath, join
from os import makedirs
from sys import exit, stdout, stderr
from datetime import datetime

#global_scaffolds = ["scaffold_72", "scaffold_32", "scaffold_60", "scaffold_19", "scaffold_45", "scaffold_17", "scaffold_74",  "scaffold_52", "scaffold_40", "scaffold_30"]

# get genotypes from the VCF file #
def Genotypes_Parser(vcf_file, module) :
	genotypes_dict = {}
	alleles_dict, variants_dict = {}, {}
	rilsID_to_Int_dict = {}
	num_variant_sites, num_variant_sites_per_chrom = 0, 0
	num_variants_kept = 0
	samples = []
	num_markers_per_chr_dict = {}
	chrom_id, last_chrom_id = "", ""
	fVCF = open(vcf_file, 'r')
	for line in fVCF :
		if line.startswith('#') :
			if line[1] != "#" :
				samples = re.split('\t', line.strip())[9:]
		else :
			tmp_line = re.split('\t', line.strip())
			chrom_id = tmp_line[0]
			if last_chrom_id == "" :
				last_chrom_id = chrom_id
			else :
				if last_chrom_id != chrom_id :
					num_markers_per_chr_dict[last_chrom_id] = num_variant_sites_per_chrom
					num_variant_sites_per_chrom = 0
					last_chrom_id = chrom_id
			num_variant_sites += 1
			variant_site_pos = tmp_line[1]
			ref_base = tmp_line[3]
			alt_base = tmp_line[4]
			num_het, num_missing = 0, 0
			genotypes = tmp_line[9:]
			for i in range(len(genotypes)) :
				individual_genotype = re.split('\:', genotypes[i])[0]
				if individual_genotype == "./." :
					num_missing += 1
				elif individual_genotype == "0/1" :
					num_het += 1
			if num_missing != len(samples) :
				if num_het/float(len(samples) - num_missing) <= 0.05 and num_missing <= 4 :
					if len(ref_base) == 1 and len(alt_base) == 1 :
						variants_dict[chrom_id+":"+variant_site_pos] = ref_base + "\t" + alt_base
						num_variants_kept += 1
						num_variant_sites_per_chrom += 1
						for i in range(len(genotypes)) :
							individual_genotype = re.split('\:', genotypes[i])[0]
							tmp_alleles = ""
							if individual_genotype == "1/1" :
								tmp_alleles = "%s_%s" %(alt_base, alt_base)
							elif individual_genotype == "0/0" :
								tmp_alleles = "%s_%s" %(ref_base, ref_base)
							elif individual_genotype == "./." :
								tmp_alleles = "?_?"
							else :
								tmp_alleles = "%s_%s" %(ref_base, alt_base)
							if not alleles_dict.has_key(samples[i]+"\t"+chrom_id) :
								alleles_dict[samples[i]+"\t"+chrom_id] = [tmp_alleles]
							else :
								alleles_dict[samples[i]+"\t"+chrom_id].append(tmp_alleles)
						###########################################
						# this is for faking the genotypes for f1 #
						###########################################
						if not alleles_dict.has_key("f1\t"+chrom_id) :
							alleles_dict["f1\t"+chrom_id] = ["%s_%s" %(ref_base, alt_base)]
						else :
							alleles_dict["f1\t"+chrom_id].append("%s_%s" %(ref_base, alt_base))
			if (num_variant_sites+1) % 10000 == 0 :
				stdout.write("\t%d variant sites being processed\n" %((((num_variant_sites+1)/10000)*10000)))
	num_markers_per_chr_dict[last_chrom_id] = num_variant_sites_per_chrom
	if (((num_variant_sites+1)/10000.0)-(num_variant_sites/10000))*10000 < 10000 :
		stdout.write("\t%d variant sites being processed\n" %((((num_variant_sites+1)/10000.0)-(num_variant_sites/10000))*10000))
	elif (num_variant_sites+1) % 10000 == 0 :
		stdout.write("\t%d variant sites being processed\n" %((((num_variant_sites+1)/10000)*10000)))
	fVCF.close()
	stdout.write("\t%d variants parsed\n" %(num_variant_sites))
	exit()
	#################
	# debug purpose #
	#################
	#num = 0
	#for key in alleles_dict.iterkeys() :
	#	if re.search("f1\t", key) :
	#		num += len(alleles_dict[key])
	#print num
	#exit()
	#if module == "prepare_mapping" :
	#	return genotypes_dict, samples, num_markers_per_chr_dict, num_variant_sites
	#elif module == "parse_for_Phasing" :
	return alleles_dict, variants_dict, ["f1"]+samples, num_markers_per_chr_dict, num_variant_sites

def outputter_gt_table(genotypes_dict, samples, num_variant_sites, out_prefix, out_format, thinflag) :
	if out_format == "mstmap" :
		if not thinflag :
			out_gt_table = out_prefix + "_%d_markers.mst" %(num_variant_sites)
		else :
			out_gt_table = out_prefix + "_%d_markers.mst" %(num_variant_sites)
	elif out_format == "rqtl" :
		if not thinflag :
			out_gt_table = out_prefix + "_%d_markers.rqtl" %(num_variant_sites)
		else :
			out_gt_table = out_prefix + "_%d_markers.rqtl" %(num_variant_sites)
	stdout.write("Output genotypes of %d individulas at %d sites into %s\n" %(len(samples), num_variant_sites, out_gt_table))
	outGTable = open(out_gt_table, 'w')
	if out_format == "mstmap" :
		outGTable.write("markers\t")
	elif out_format == "rqtl" :
		outGTable.write(",")
	for i in range(len(samples)) :
		if out_format == "mstmap" :
			outGTable.write("\t%s" %(samples[i]))
		else :
			outGTable.write(",%s" %(samples[i]))
	outGTable.write('\n')
	for i in range(num_variant_sites) :
		marker_pos = ""
		for individual, gt_info in sorted(genotypes_dict.iteritems()) :
			tmp = re.split(':', gt_info[i])
			if thinflag :
				gt = tmp[1]
			else :
				gt = tmp[2]
			if marker_pos == "" :
				if thinflag :
					marker_pos = tmp[0]
				else :
					marker_pos = re.sub('_', "at", tmp[0]) + "at" + tmp[1]
				if out_format == "rqtl" :
					outGTable.write("%s,aaxbb" %(marker_pos))
				else :
					outGTable.write("%s" %(marker_pos))
			if out_format == "rqtl" :
				if gt == "X" :
					outGTable.write(",H")
				else :
					outGTable.write(",%s" %(gt))
			else :
				outGTable.write("\t%s" %(gt))
		outGTable.write("\n")
		outGTable.flush()
	outGTable.close()

def thinning_genotyp_table(gt_dict, num_markers_per_chr_dict, num_variant_sites, out_prefix, out_format) :
	"""
	thinning out the genotype table and prepare the table ready for MSTMAP and rQTL
	"""

	stdout.write("Start thinning the genotypes\n")
	individual, last_individual = "", ""
	breaks_individual_dict = {}
	for individual, gt_info in sorted(gt_dict.iteritems()) :
		gt, last_gt = "", ""
		step = 0
		i = 0
		chrom_id, last_chrom_id = "", ""
		while i < len(gt_info) :
			tmp = re.split(":", gt_info[i])
			marker_index = re.sub("_", "at", tmp[0]) + "at" + tmp[1]
			chrom_id = tmp[0]
			gt = tmp[2]
			if last_chrom_id == "" :
				last_chrom_id = chrom_id
				step += num_markers_per_chr_dict[last_chrom_id]
			else :
				if last_chrom_id != chrom_id :
					step += num_markers_per_chr_dict[chrom_id]
					if num_markers_per_chr_dict[chrom_id] <= 10 :
					#if num_markers_per_chr_dict[chrom_id] <= 30 :
						tmp1 = re.split(":", gt_info[i-1])
						tmp_marker_index = re.sub("_", "at", tmp1[0]) + "at" + tmp1[1]
						breaks_individual_dict[individual].append("%s:%s" %(tmp_marker_index, last_gt))
						if i == step - 1 :
							breaks_individual_dict[individual].append("%s:%s" %(marker_index, gt))
						else :
							j = i
							while j < step - 1 :
								tmp2 = re.split(":", gt_info[j])
								tmp_marker_index2 = re.sub("_", "at", tmp2[0]) + "at" + tmp2[1]
								breaks_individual_dict[individual].append("%s:%s" %(tmp_marker_index2, tmp2[2]))
								j += 1
							i = j
					else :
						tmp1 = re.split(":", gt_info[i-1])
						tmp_marker_index = re.sub("_", "at", tmp1[0]) + "at" + tmp1[1]
						breaks_individual_dict[individual].append("%s:%s" %(tmp_marker_index, last_gt))
					last_gt = ""
				last_chrom_id = chrom_id
			if gt == "X" or gt == "U" :
				if last_gt == "" :			# in case the first genotype is X or U
					j = i + 1
					tmp_alternative_gt, tmp_alternative_gt_dict, tmp_pos_alternative_gt_dict = "", {}, {}
					while j < step and j < i+10 :
					#while j < step and j < i+5 :
						if not tmp_alternative_gt_dict.has_key(re.split(":", gt_info[j])[2]) :
							tmp_alternative_gt_dict[re.split(":", gt_info[j])[2]] = 1
							tmp_pos_alternative_gt_dict[re.split(":", gt_info[j])[2]] = j
						else :
							tmp_alternative_gt_dict[re.split(":", gt_info[j])[2]] += 1
							tmp_pos_alternative_gt_dict[re.split(":", gt_info[j])[2]] = j
						j += 1
					if not tmp_alternative_gt_dict.has_key(gt) :
						tmp_alternative_gt_dict[gt] = 0
						tmp_pos_alternative_gt_dict[gt] = i + 1
					tmp_max, tmp_max_gt = 0, ""
					tmp_max_gt = gt
					for tmp_gt in tmp_alternative_gt_dict.iterkeys() :
						if tmp_gt != last_gt :
							#if tmp_alternative_gt_dict[tmp_gt] > tmp_alternative_gt_dict[gt] and tmp_alternative_gt_dict[tmp_gt] >= 4 :
							if tmp_alternative_gt_dict[tmp_gt] > tmp_alternative_gt_dict[gt] and tmp_alternative_gt_dict[tmp_gt] >= 4 :
								tmp_max_gt = tmp_gt
								tmp_max = tmp_alternative_gt_dict[tmp_gt]
					last_gt = tmp_max_gt
					if j < step - 1 :
						i = tmp_pos_alternative_gt_dict[last_gt]
					else :
						i = j
				else :
					if last_gt == gt :
						i += 1
					else :
						j = i + 1
						tmp_num_diff = 0
						tmp_alternative_gt, tmp_alternative_gt_dict, tmp_pos_alternative_gt_dict = "", {}, {}
						while j < i+ 15 and j < step :
						#while j < i+ 5 and j < step :
							if not tmp_alternative_gt_dict.has_key(re.split(":", gt_info[j])[2]) :
								tmp_alternative_gt_dict[re.split(":", gt_info[j])[2]] = 1
								tmp_pos_alternative_gt_dict[re.split(":", gt_info[j])[2]] = j
							else :
								tmp_alternative_gt_dict[re.split(":", gt_info[j])[2]] += 1
								tmp_pos_alternative_gt_dict[re.split(":", gt_info[j])[2]] = j
							j += 1
						if not tmp_alternative_gt_dict.has_key(last_gt) :
							tmp_alternative_gt_dict[last_gt] = 0
							tmp_pos_alternative_gt_dict[last_gt] = i + 1
						tmp_max, tmp_max_gt = 0, ""
						tmp_max_gt = last_gt
						for tmp_gt in tmp_alternative_gt_dict.iterkeys() :
							if tmp_gt != last_gt :
								if tmp_alternative_gt_dict[tmp_gt] > tmp_alternative_gt_dict[last_gt] and tmp_alternative_gt_dict[tmp_gt] >= 4 :
									tmp_max_gt = tmp_gt
									tmp_max = tmp_alternative_gt_dict[tmp_gt]
						if tmp_max_gt != last_gt :
							tmp1 = re.split(":", gt_info[i-1])
							tmp_marker_index = re.sub("_", "at", tmp1[0]) + "at" + tmp1[1]
							if not breaks_individual_dict.has_key(individual) :
								breaks_individual_dict[individual] = ["%s:%s" %(tmp_marker_index, last_gt)]
							else :
								breaks_individual_dict[individual].append("%s:%s" %(tmp_marker_index, last_gt))
							last_gt = tmp_max_gt
							i = tmp_pos_alternative_gt_dict[tmp_max_gt]
						else :
							if j < step - 1 :
								i = tmp_pos_alternative_gt_dict[last_gt]
							else :
								i = j
			else :
				# first non-U/X genotype #
				if last_gt == "" :
					j = i + 1
					tmp_alternative_gt, tmp_alternative_gt_dict, tmp_pos_alternative_gt_dict = "", {}, {}
					#while j < step and j < i+10 :
					while j < step and j < i+3 :
						if not tmp_alternative_gt_dict.has_key(re.split(":", gt_info[j])[2]) :
							tmp_alternative_gt_dict[re.split(":", gt_info[j])[2]] = 1
							tmp_pos_alternative_gt_dict[re.split(":", gt_info[j])[2]] = j
						else :
							tmp_alternative_gt_dict[re.split(":", gt_info[j])[2]] += 1
							tmp_pos_alternative_gt_dict[re.split(":", gt_info[j])[2]] = j
						j += 1
					if not tmp_alternative_gt_dict.has_key(gt) :
						tmp_alternative_gt_dict[gt] = 0
						tmp_pos_alternative_gt_dict[gt] = i + 1
					tmp_max, tmp_max_gt = 0, ""
					tmp_max_gt = gt
					for tmp_gt in tmp_alternative_gt_dict.iterkeys() :
						if tmp_gt != last_gt :
							if tmp_alternative_gt_dict[tmp_gt] > tmp_alternative_gt_dict[gt] and tmp_alternative_gt_dict[tmp_gt] >= 2 :
								tmp_max_gt = tmp_gt
								tmp_max = tmp_alternative_gt_dict[tmp_gt]
					last_gt = tmp_max_gt
					if j < step - 1 :
						i = tmp_pos_alternative_gt_dict[last_gt]
					else :
						i = j
				else :
					if last_gt == gt :
						i += 1
					else :
						j = i + 1
						tmp_num_diff = 0
						tmp_alternative_gt, tmp_alternative_gt_dict, tmp_pos_alternative_gt_dict = "", {}, {}
						#while j < step and j < i + 12 :
						while j < step and j < i + 3 :
							if not tmp_alternative_gt_dict.has_key(re.split(":", gt_info[j])[2]) :
								tmp_alternative_gt_dict[re.split(":", gt_info[j])[2]] = 1
								tmp_pos_alternative_gt_dict[re.split(":", gt_info[j])[2]] = j
							else :
								tmp_alternative_gt_dict[re.split(":", gt_info[j])[2]] += 1
								tmp_pos_alternative_gt_dict[re.split(":", gt_info[j])[2]] = j
							j += 1
						if not tmp_alternative_gt_dict.has_key(last_gt) :
							tmp_alternative_gt_dict[last_gt] = 0
							tmp_pos_alternative_gt_dict[last_gt] = i + 1
						tmp_max, tmp_max_gt = 0, ""
						tmp_max_gt = last_gt
						for tmp_gt in tmp_alternative_gt_dict.iterkeys() :
							if tmp_gt != last_gt :
								#if tmp_alternative_gt_dict[tmp_gt] > tmp_alternative_gt_dict[last_gt] and tmp_alternative_gt_dict[tmp_gt] >= 6 :
								if tmp_alternative_gt_dict[tmp_gt] > tmp_alternative_gt_dict[last_gt] and tmp_alternative_gt_dict[tmp_gt] >= 2 :
									tmp_max_gt = tmp_gt
									tmp_max = tmp_alternative_gt_dict[tmp_gt]
						if tmp_max_gt != last_gt :
							tmp1 = re.split(":", gt_info[i-1])
							tmp_marker_index = re.sub("_", "at", tmp1[0]) + "at" + tmp1[1]
							if not breaks_individual_dict.has_key(individual) :
								breaks_individual_dict[individual] = ["%s:%s" %(tmp_marker_index, last_gt)]
							else :
								breaks_individual_dict[individual].append("%s:%s" %(tmp_marker_index, last_gt))
							last_gt = tmp_max_gt
							i = tmp_pos_alternative_gt_dict[tmp_max_gt]
						else :
							if j < step-1 :
								i = tmp_pos_alternative_gt_dict[last_gt]
							else :
								i = j

	overall_breaks_list = []
	for individual, breaks in sorted(breaks_individual_dict.iteritems()) :
		for i in range(len(breaks)) :
			if re.split(":", breaks[i])[0] not in overall_breaks_list :
				overall_breaks_list.append(re.split(":", breaks[i])[0])
	num_variant_sites = len(overall_breaks_list)

	breaks_per_chrom_dict = {}
	for i in range(len(sorted(overall_breaks_list))) :
		chrom_id = re.split("at", overall_breaks_list[i])[0] + "at" + re.split("at", overall_breaks_list[i])[1]
		marker_pos = int(re.split("at", overall_breaks_list[i])[2])
		if not breaks_per_chrom_dict.has_key(chrom_id) :
			breaks_per_chrom_dict[chrom_id] = [marker_pos]
		else :
			breaks_per_chrom_dict[chrom_id].append(marker_pos)
	for key in breaks_per_chrom_dict.iterkeys() :
		breaks_per_chrom_dict[key] = sorted(breaks_per_chrom_dict[key])

	thinned_gt_table_dict, thinned_num_markers_dict = {}, {}
	for individual, marker_info in sorted(breaks_individual_dict.iteritems()) :
		chrom_id, last_chrom_id = "", ""
		tmp_break_pos_list = []
		num_marker_per_chrom_after_thinning = 0
		j = 0
		for i in range(len(marker_info)) :
			tmp = re.split(":", marker_info[i])
			chrom_id = re.split("at", tmp[0])[0] + "at" + re.split("at", tmp[0])[1]
			if last_chrom_id == "" :
				last_chrom_id = chrom_id
			else :
				if last_chrom_id != chrom_id :
					num_marker_per_chrom_after_thinning = 0
					last_chrom_id = chrom_id
					j = 0
			break_pos = int(re.split("at", tmp[0])[2])
			break_gt = tmp[1]
			while j < len(breaks_per_chrom_dict[chrom_id]) :
				if breaks_per_chrom_dict[chrom_id][j] <= break_pos :
					if not thinned_gt_table_dict.has_key(individual) :
						thinned_gt_table_dict[individual] = ["%sat%d:%s" %(chrom_id, breaks_per_chrom_dict[chrom_id][j], break_gt)]
					else :
						thinned_gt_table_dict[individual].append("%sat%d:%s" %(chrom_id, breaks_per_chrom_dict[chrom_id][j], break_gt))
				else :
					break
				j += 1

	thinning_logfile = out_prefix + "_thinned.log"
	fLOG = open(thinning_logfile, 'w')
	for chrom in sorted(breaks_per_chrom_dict.iterkeys()) :
		fLOG.write("%s\t%s\t%s\n" %(chrom, num_markers_per_chr_dict[re.sub("at", "_", chrom)], len(breaks_per_chrom_dict[chrom])))
	fLOG.close()
	stdoutleft.write("\t%d variant sites kept\n" %(num_variant_sites))
	return thinned_gt_table_dict, num_variant_sites

def getKnownLinkageGroups(known_LGs_file) :
	known_LGs_dict, embedded_scaffolds_dict = {}, {}
	scaffolds_orientation_dict = {}
	scaffolds = []
	lg_id = ""
	fLG = open(known_LGs_file, 'r')
	for line in fLG :
		if line.startswith(">") :
			if len(scaffolds) != 0 :
				known_LGs_dict[lg_id] = scaffolds
				scaffolds = []
			lg_id = line.strip()[1:]
		else :
			tmp_line = re.split('\t', line.strip())
			if tmp_line[1] == "embedded" :
				if not embedded_scaffolds_dict.has_key(tmp_line[2]) :
					embedded_scaffolds_dict[tmp_line[2]] = [tmp_line[0]]
				else :
					embedded_scaffolds_dict[tmp_line[2]].append(tmp_line[0])
			elif tmp_line[1] != "not_embedded" :
				if tmp_line[1] != "na" :
					scaffolds.append(tmp_line[0])
					scaffolds_orientation_dict[tmp_line[0]] = tmp_line[1]
	if len(scaffolds) != 0 :
		known_LGs_dict[lg_id] = scaffolds
		scaffolds = []
	fLG.close()
	print embedded_scaffolds_dict
	print scaffolds_orientation_dict
	for key, value in known_LGs_dict.iteritems() :
		print key, known_LGs_dict[key]
	return known_LGs_dict, embedded_scaffolds_dict, scaffolds_orientation_dict

def getNumVariantsForLG(num_markers_per_chr_dict, scaffolds_list) :
	num_variant_sites_per_lg = 0
	for i in range(len(scaffolds_list)) :
		#if scaffolds_list[i] in global_scaffolds :
		for each_scaffold in num_markers_per_chr_dict.iterkeys() :
			if scaffolds_list[i] == each_scaffold :
				num_variant_sites_per_lg += num_markers_per_chr_dict[each_scaffold]
	return num_variant_sites_per_lg

# prepare genotypes file for phasing program #
def prepare_gtyps_for_Phasing(alleles_dict, samples, scaffolds_list, outfile_handle, embedded_scaffolds_dict, scaffolds_orientation_dict) :
	for i in range(len(samples)) :
		outfile_handle.write("#id_%s\n" %(samples[i]))
		first_alleles, second_alleles = [], []
		total_len_first, total_len_second = 0, 0
		for j in range(len(scaffolds_list)) :
		#	if scaffolds_list[j] in global_scaffolds :
			key = samples[i] + "\t" + scaffolds_list[j]
			tmp_first_alleles, tmp_second_alleles = [], []
			for k in range(len(alleles_dict[key])) :
				tmp = re.split("_", alleles_dict[key][k])
				tmp_first_alleles.append(tmp[0])
				tmp_second_alleles.append(tmp[1])
			if scaffolds_orientation_dict[scaffolds_list[j]] == '-' :
				tmp_first_alleles = tmp_first_alleles[::-1]
				tmp_second_alleles = tmp_second_alleles[::-1]
			first_alleles += tmp_first_alleles
			second_alleles += tmp_second_alleles
		for m in range(len(first_alleles)) :
			outfile_handle.write("%s" %(first_alleles[m]))
		outfile_handle.write("\n")
		for m in range(len(second_alleles)) :
			outfile_handle.write("%s" %(second_alleles[m]))
		outfile_handle.write("\n")
	outfile_handle.close()

def prepare_inputs_for_beagle(alleles_dict, samples, scaffolds_list, outfile, embedded_scaffolds_dict, scaffolds_orientation_dict, num_markers_per_chr_dict) :
	fOUT = open(outfile, 'w')
	fOUT.write('I id')
	for i in range(len(samples)) :
		for j in range(0,2) :
			fOUT.write(" %s" %(samples[i]))
	fOUT.write('\n')
	cumulative = 0
	for i in range(len(scaffolds_list)) :
		#if scaffolds_list[i] in global_scaffolds :
		tmp_alleles_dict = alleles_dict
		if scaffolds_orientation_dict[scaffolds_list[i]] == '-' :
			for j in range(len(samples)) :
				tmp_alleles_dict[samples[j]+'\t'+scaffolds_list[i]] = alleles_dict[samples[j]+'\t'+scaffolds_list[i]][::-1]
		for j in range(num_markers_per_chr_dict[scaffolds_list[i]]) :
			fOUT.write('M %s_%d' %(scaffolds_list[i], j))
			for k in range(len(samples)) :
				key = samples[k] + "\t" + scaffolds_list[i]
				first_allele = re.split('_', tmp_alleles_dict[key][j])[0]
				second_allele = re.split('_', tmp_alleles_dict[key][j])[1]
				#print key, tmp_alleles_dict[key][j]
				#print first_allele, second_allele
				fOUT.write(" %s %s" %(first_allele, second_allele))
			fOUT.write('\n')
			cumulative += 1
	fOUT.close()

# output indels information #
def outputter_variants(variants_dict, scaffolds_list, variants_outfile, scaffolds_orientation_dict) :
	fINDELOUT = open(variants_outfile, 'w')
	for i in range(len(scaffolds_list)) :
		#if scaffolds_list[i] in ["scaffold_19", "scaffold_60", "scaffold_52"] :
		#if scaffolds_list[i] in global_scaffolds :
		tmp_pos_list = []
		for key in variants_dict.iterkeys() :
			each_scaffold = re.split(":", key)[0]
			if each_scaffold == scaffolds_list[i] :
				pos = int(re.split(":", key)[1])
				tmp_pos_list.append(pos)
		if scaffolds_orientation_dict[scaffolds_list[i]] == '+' :
			tmp_pos_list = sorted(tmp_pos_list)
		else :
			tmp_pos_list = sorted(tmp_pos_list, reverse=True)
		for j in range(len(tmp_pos_list)) :
			if variants_dict.has_key(scaffolds_list[i]+":"+str(tmp_pos_list[j])) :
				fINDELOUT.write("%s\t%s\n" %(scaffolds_list[i]+":"+str(tmp_pos_list[j]), variants_dict[scaffolds_list[i]+":"+str(tmp_pos_list[j])]))
	fINDELOUT.close()

def timestamper() :
	"""
	generate a dash-separated time stamped string
	"""
	return datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

def make_dir_if_needed(dir) :
	"""make a directory if necessary """
	if not exists(dir) :
		makedirs(dir)

def command_parser() :
	""" parsing command line """
	parser = ArgumentParser(description="Parsing VCF file for multiple downstream analysis")
	parser.add_argument("-analysis", metavar="STR", dest="analysis", choices=["prepare_mapping", "prepare_phasing"], required=True, help="available analysis: prepare_mapping, prepare_phasing")
	parser.add_argument("-vcf", metavar="FILE", dest="vcf_file", required=True, help="VCF file generated from GATK")
	parser.add_argument("-p", metavar="STR", dest="out_prefix", required=True, help="prefix of output genotype table file in either MSTMAP or rQTL required input format")

	# add option group for analysis "prepare_mapping"
	prepare_mapping_args_group = parser.add_argument_group("Options for \"prepare_mapping\" analysis", "generate genotype tables that are ready for downstream MSTMAP and/or rQTL\n")
	prepare_mapping_args_group.add_argument("-thinoff", dest="thin_or_not", action='store_true', help="turn off thining the genotype table. Not suggested if you get a big table")
	#prepare_mapping_args_group.add_argument("-thin_mode", dest="thin_mode", choices=["sliding_window", "random"], default="sliding_window", help="specify how you want to thin your genotype table. Current version only supports sliding_window. Other modes are coming soon")
	prepare_mapping_args_group.add_argument("-gt_outfmt", metavar="STR", dest="gt_outfmt", choices=["mstmap", "rqtl", "both"], type=str, default="both", help="specify the format of the output genotype table. Current version supports: mstmap, rqtl, and both. By default, genotype tables in both formats will be generated")
	prepare_mapping_args_group.add_argument("-gt_proposal", metavar="FILE", dest="gt_proposal_file", help="file of genotypes proposed at variant sites")

	# add option group for analysis "prepare_phasing"
	prepare_phasing_args_group = parser.add_argument_group("Options for \"prepare_phasing\" analysis", "generate genotype file for downstream Phasing programs\n")
	prepare_phasing_args_group.add_argument("-known_LGs", metavar="FILE", dest="known_LGs_file", help="a file of a set of known linkage groups. If available, the output file generated for phasing program will be categorized by linkage group")
	prepare_phasing_args_group.add_argument("-phasing_outfmt", metavar="STR", dest="phasing_outfmt", choices=["fastphase", "beagle"], type=str, default="fastphase", help="specify the format of the output genotype file for the following Phasing programs. Current version supports: fastphase, beagle")

	args = parser.parse_args()

	# make sure vcf file provided exist #
	if not exists(args.vcf_file) :
		stderr.write(timestamper() + " [IO Error]: Cannot find the vcf file you provided %s\n" %(args.vcf_file))
		exit()

	# make sure the output path provided is reachable #
	make_dir_if_needed(dirname(realpath(args.out_prefix)))

	if args.analysis == "prepare_mapping" :
		run_prepare_mapping(args.analysis, args.vcf_file, args.gt_proposal_file, realpath(args.out_prefix), args.gt_outfmt, args.thin_or_not)
	elif args.analysis == "prepare_phasing" :
		run_prepare_phasing(args.analysis, args.vcf_file, args.known_LGs_file, realpath(args.out_prefix), args.phasing_outfmt)

def run_prepare_mapping(module, vcf_file, gt_proposal_file, out_prefix, gt_outfmt, thin_or_not) :
	stdout.write("###################################################################\n")
	stdout.write("#Analysis: %s\n" %(module))
	stdout.write("#VCF: %s\n" %(vcf_file))
	stdout.write("#Genotype Table output: %s\n" %(out_prefix))
	stdout.write("#Ouput Format: %s\n" %(gt_outfmt))
	if thin_or_not :
		stdout.write("#Thin the Table: No\n")
	else :
		stdout.write("#Thin the Table: Yes\n")
	stdout.write("###################################################################\n")
	gt_dict, samples, num_markers_per_chr_dict, num_variant_sites = prepare_genotype_for_mapping(vcf_file, gt_proposal_file)
	if not thin_or_not :
		thinned_gt_table_dict, num_variant_sites_after_thinning = thinning_genotyp_table(gt_dict, num_markers_per_chr_dict, num_variant_sites, out_prefix, gt_outfmt)
		if gt_outfmt == "both" :
			outputter_gt_table(thinned_gt_table_dict, samples, num_variant_sites_after_thinning, out_prefix, "mstmap", 1)
			outputter_gt_table(thinned_gt_table_dict, samples, num_variant_sites_after_thinning, out_prefix, "rqtl", 1)
		else :
			outputter_gt_table(thinned_gt_table_dict, samples, num_variant_sites_after_thinning, out_prefix, gt_outfmt, 1)
	else :
		if gt_outfmt == "both" :
			outputter_gt_table(gt_dict, samples, num_variant_sites, out_prefix, "mstmap", 0)
			outputter_gt_table(gt_dict, samples, num_variant_sites, out_prefix, "rqtl", 0)
		else :
			outputter_gt_table(gt_dict, samples, num_variant_sites, out_prefix, gt_outfmt, 0)

def prepare_genotype_for_mapping(vcf_file, gt_proposal_file = None) :
	proposed_genotypes_from_RPGC = {}
	if gt_proposal_file :
		stdout.write("Reading proposed genotypes from %s\n" %(gt_proposal_file))
		fRPGC = open(gt_proposal_file, 'r')
		for line in fRPGC :
			tmp_line = re.split('\t', line.strip())
			chrom_id = tmp_line[0]
			variant_pos = tmp_line[1]
			site = chrom_id + ':' + variant_pos
			proposed_genotypes_from_RPGC[site] = "\t".join(tmp_line[2:])
		fRPGC.close()

	stdout.write("Reading VCF file %s\n" %(vcf_file))
	samples = []
	num_markers_per_chr_dict = {}
	genotypes_dict = {}
	chrom_id, last_chrom_id = "", ""
	num_variant_sites_per_chrom, num_variant_sites = 0, 0
	num_informative_markers = 0
	fVCF = open(vcf_file, 'r')
	for line in fVCF :
		if line.startswith('#') :
			if line.startswith("#CHROM") :
				tmp_line = re.split('\t', line.strip())
				samples = tmp_line[9:]
		else :
			tmp_line = re.split('\t', line.strip())
			chrom_id = tmp_line[0]
			if last_chrom_id == "" :
				last_chrom_id = chrom_id
			else :
				if last_chrom_id != chrom_id :
					num_markers_per_chr_dict[last_chrom_id] = num_variant_sites_per_chrom
					num_variant_sites_per_chrom = 0
					last_chrom_id = chrom_id
			variant_pos = tmp_line[1]
			site = chrom_id + ':' + variant_pos
			ref_base = tmp_line[3]
			alt_base = tmp_line[4]
			if len(re.split('\,', ref_base)) > 1 or len(re.split('\,', alt_base)) > 1 :
				continue
			if proposed_genotypes_from_RPGC.has_key(site) :
				samples_genotypes = proposed_genotypes_from_RPGC[site]
				site_type = "proposed"
			else :
				samples_genotypes = '\t'.join(tmp_line[9:])
				site_type = "original"
			num_homo, num_samples_with_available_gt = 0, len(samples)
			tmp_genotypes_dict = {}
			for i in range(len(samples)) :
				tmp_genotypes = re.split('\t', samples_genotypes)
				individual_sample = re.split(':', tmp_genotypes[i])
				individual_sample_genotype = individual_sample[0]
				if individual_sample_genotype == "0/0" :
					one_letter_genotype = "A"
					num_homo += 1
				elif individual_sample_genotype == "1/1" :
					one_letter_genotype = "B"
					num_homo += 1
				elif individual_sample_genotype == "0/1" :
					one_letter_genotype = "X"
				else :
					one_letter_genotype = "U"

				tmp_genotypes_dict[samples[i]] = chrom_id+":"+variant_pos+":"+one_letter_genotype

			if site_type != "proposed" :
				if tmp_line[6] == "PASS" and num_samples_with_available_gt/float(len(samples)) >= 0.8 :
					if float(num_homo)/num_samples_with_available_gt >= 0.95 :
						num_informative_markers += 1
						num_variant_sites_per_chrom += 1
						num_variant_sites += 1
						for sample in tmp_genotypes_dict.iterkeys() :
							if not genotypes_dict.has_key(sample) :
								genotypes_dict[sample] = [tmp_genotypes_dict[sample]]
							else :
								genotypes_dict[sample].append(tmp_genotypes_dict[sample])
			else :
				num_informative_markers += 1
				num_variant_sites_per_chrom += 1
				num_variant_sites += 1
				for sample in tmp_genotypes_dict.iterkeys() :
					if not genotypes_dict.has_key(sample) :
						genotypes_dict[sample] = [tmp_genotypes_dict[sample]]
					else :
						genotypes_dict[sample].append(tmp_genotypes_dict[sample])
	num_markers_per_chr_dict[last_chrom_id] = num_variant_sites_per_chrom
	fVCF.close()
	stdout.write("\t%d variants parsed\n" %(num_informative_markers))

	return genotypes_dict, samples, num_markers_per_chr_dict, num_variant_sites

def run_prepare_phasing(module, vcf_file, known_LGs_file, out_prefix, phasing_outfmt) :
	stdout.write("###################################################################\n")
	stdout.write("#Analysis: %s\n" %(module))
	stdout.write("#VCF: %s\n" %(vcf_file))
	stdout.write("#output prefix: %s\n" %(out_prefix))
	stdout.write("#Ouput Format: %s\n" %(phasing_outfmt))
	stdout.write("###################################################################\n")
	alleles_dict, variants_dict, samples, num_markers_per_chr_dict, num_variant_sites = Genotypes_Parser(vcf_file, module)
	known_LGs_dict, embedded_scaffolds_dict = {}, {}
	outfile = ""
	if known_LGs_file != "" :
		known_LGs_dict, embedded_scaffolds_dict, scaffolds_orientation_dict = getKnownLinkageGroups(known_LGs_file)
		for each_lg, scaffolds_list in known_LGs_dict.iteritems() :
			#if each_lg == "lg_5" :
			num_variant_sites_per_lg = getNumVariantsForLG(num_markers_per_chr_dict, scaffolds_list)
			print each_lg, num_variant_sites_per_lg
			variants_outfile = out_prefix + "%s.variants" %(each_lg)
			outputter_variants(variants_dict, scaffolds_list, variants_outfile, scaffolds_orientation_dict)
			if phasing_outfmt == "fastphase" :
				outfile = out_prefix + "%s.inp" %(each_lg)
				fOUT = open(outfile, 'w')
				fOUT.write("%d\n" %(len(samples)))
				fOUT.write("%d\n" %(num_variant_sites_per_lg))
				prepare_gtyps_for_Phasing(alleles_dict, samples, scaffolds_list, fOUT, embedded_scaffolds_dict, scaffolds_orientation_dict)
			elif phasing_outfmt == "beagle" :
				outfile = out_prefix + "%s.bgl.unphased" %(each_lg)
				prepare_inputs_for_beagle(alleles_dict, samples, scaffolds_list, outfile, embedded_scaffolds_dict, scaffolds_orientation_dict, num_markers_per_chr_dict)
	else :
		if phasing_outfmt == "fastphase" :
			OUTFILE = Out_prefix + ".inp"
		scaffolds_list = []
		for each_scaffold in num_markers_per_chr_dict.iterkeys() :
			scaffolds_list.append(each_scaffold)
		fOUT = open(outfile, 'w')
		fOUT.write("%d\n" %(len(samples)))
		fOUT.write("%d\n" %(num_variant_sites))
		prepare_gtyps_for_Phasing(alleles_dict, samples, scaffolds_list, fOUT)

if __name__ == "__main__" :
	command_parser()
