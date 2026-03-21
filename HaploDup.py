#!/usr/bin/env python

import argparse
import concurrent.futures
import shlex
from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.GFF_lib import *
from lib_files.FASTA_lib import *
from lib_files.map_lib import *



def main() :

	###### Options and help ######
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta", default=False, dest="fasta",
						help="FASTA file(s) with genomic sequences", metavar="genome.fasta [Required]")
	parser.add_argument("-g", "--gff", default=False, dest="gff",
						help="Annotation file(s) in GFF3 format", metavar="annotation.gff3 [Required]")
	parser.add_argument("-a", "--annotation", default=False, dest="annotation",
						help="Text file with the functional annotation of each transcript", metavar="functional_annotation.txt")
	parser.add_argument("-c", "--correspondence", dest="corr", default=False ,
						help="Tab text file(s) of corresponding sequence names between haplotypes for each chromosome and, if present, the reference genome. Tab separated file with 3/4 columns: 1) Chromosome ID, 2) Haplotype 1 ID, 3) Haplotype 2 ID, 4) [optional] Reference ID [Required]", metavar="Hap1_to_Hap2.txt")

	# Add marker BED input
	parser.add_argument( "-b" , "--markers_bed" , dest="markers_hits" ,
						help="Marker position on input sequences", metavar="markers.bed")
	parser.add_argument("--markers_map", dest="marker_map",
					help="Table of markers order when sorted genomic position. Tab separated file with 3 columns: 1) Chromosome ID, 2) Marker sorting position, 3) Marker ID.", metavar="markers_list")

	# Add pseudomolecule structure AGP
	parser.add_argument("--agp" , dest="agp" ,
						help="AGP structure of input sequences", metavar="structure.agp")
	parser.add_argument("--input_groups" , dest="input_groups",
					help="Tab separated values file of group association of input sequences, it will be used for structure QC if '--legacy_groups' is not given. Two columns reporting the sequences ID and the association groups id. Sequence can be associated to multiple groups with multiple rows."  , metavar="input_groups.tsv")

	# Add pseudomolecule legacy structure AGP
	parser.add_argument("--legacy_agp" , dest="legacy_agp" ,
						help="AGP structure of input sequences based on older component sequences", metavar="legacy_structure.agp")
	parser.add_argument("--legacy_groups" , dest="legacy_groups",
					help="Tab separated values file of group association of components present the input sequences. Two columns reporting the component ID and the association group it belongs to. Sequence can be associated to multiple groups with multiple rows.", metavar="legacy_groups.tsv")


	# QC for rejected
	parser.add_argument("--rejected_list" , dest="rejected" ,
						help="List of rejected sequences for each chromosome, run QC of rejection", metavar="rejected.list")

	parser.add_argument("-o", "--out", dest="out", default="out",
						help="Output files prefix [default: out]", metavar="NAME")

	parser.add_argument("--hit_identity", dest="hit_identity", default="90",
						help="Mimimum genome mapping hit identity [default: 90]", metavar="N")
	parser.add_argument("--hit_coverage", dest="hit_length", default="3000",
						help="Mimimum genome mapping hit length [default: 3000]", metavar="N")

	parser.add_argument("--gene_identity", dest="gene_identity", default="95",
						help="Mimimum gene mapping identity [default: 95]", metavar="N")
	parser.add_argument("--gene_coverage", dest="gene_coverage", default="95",
						help="Mimimum gene mapping coverage [default: 95]", metavar="N")

	parser.add_argument("--unbalanced_ratio", dest="unbalanced_ratio" , default="0.33",
						help="Gene count ratio between haplotype to call the locus underrepresented [values range: 0-1]", metavar="N")

	parser.add_argument("-r", "--reference", default=False, dest="reference",
						help="FASTA file(s) of reference genomes sequences (haploid)", metavar="reference.fasta")

	parser.add_argument("--reuse_mappings" , dest="reuse_mappings", default=False, action="store_true",
						help="If set, alignments present in the output folder are reused and not overwritten by performing alignments again [Default: overwrite]")
	parser.add_argument("--reuse_dotplots" , dest="reuse_dotplots", default=False, action="store_true",
						help="If set, dotplots present in the output folder are reused and not overwritten [Default: overwrite]")
	parser.add_argument("--reuse_gmap" , dest="reuse_gmap", default=False, action="store_true",
						help="If set, CDS mapping with GMAP are reused and not overwritten by performing again the analysis [Default: overwrite]")
	parser.add_argument("--skip_chr_pair_reports" , dest="skip_chr_pair_reports", default=False, action="store_true",
					help="Skip generation of the Hap1 vs Hap2 overview report for each chromosome pair [Default: generate]")
	parser.add_argument("--skip_dotplots_by_chr" , dest="skip_dotplots", default=False, action="store_true",
						help="If set, skips the generation of all individual chromosome-vs-chromosome dotplots (including both paired and unpaired comparisons) and only produces the whole-genome dotplot (all chromosomes combined in one plot). [Default: generate all individual chr-vs-chr dotplots]")
	parser.add_argument("--only_paired_dotplots" , dest="only_paired_dotplots", default=False, action="store_true",
						help="If set, only generates per-chromosome dotplots for matched chromosome pairs (e.g. chr01.hap1 vs chr01.hap2). Off-diagonal comparisons (e.g. chr01 vs chr17) are skipped. The whole-genome dotplot is always generated. [Default: generate all chr-vs-chr combinations]")

	# TODO: allow to use a custom set of CDS sequences instead of annotations
	parser.add_argument("--cds", default=False, dest="cds",
						#help="CDS sequences to use to generate a temporary annotation", metavar="cds.fasta [Required]")
						help=argparse.SUPPRESS )

	# Mapper selection is hidden as only gmap alignment is supported
	parser.add_argument("-m" , "--mapper", dest="mapper", default="gmap",
						#help="Mapping tool to use [Default: gmap]" , metavar="[gmap|blat]" )
						help=argparse.SUPPRESS )
	parser.add_argument("-t", "--threads", dest="cores", default=4,
						help="Cores used in mapping process [default: 4]", metavar="N")
	parser.add_argument("--feature", dest="feature", default="CDS",
						help="If GFF is used, feature type to use for mapping. Choice of CDS or mRNA [default: CDS]", metavar="[CDS|mRNA]")
	parser.add_argument("-w", "--window", dest="window", default=10,
						help="Window size (number of genes) for search of blocks of genes with unbalanced count between alleles. 0 disables search [default: 10]", metavar="N")
	parser.add_argument("--allowed", dest="allowed", default=5,
						help="Allowed number of unbalance gene per window. [default: 5]", metavar="N")



	scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/support_scripts"
	print("Running HaploDup tool from HaploSync version " + get_version(), file=sys.stdout)
	print("To reproduce this run use the following command: " + " ".join( shlex.quote(x) for x in sys.argv), file=sys.stdout)
	print("----", file=sys.stdout)
	# Sanity Check

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()

	if not options.fasta :
		print("[ERROR] Genome FASTA file(s) missing", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	#if not ( options.gff or options.cds ) :
	#	print >> sys.stderr , "[ERROR] Annotation file missing"
	#	parser.print_help()
	#	sys.exit(1)

	if not options.corr :
		print("[ERROR] Correspondence file missing ", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	haplodup_dir = options.out + ".HaploDup_dir"

	if (options.reuse_mappings or options.reuse_dotplots or options.reuse_gmap) and (not os.path.exists(haplodup_dir)):
		print("[ERROR] Required to reuse output folder content but output folder do not exist (prefix used: " + options.out + ")", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	paths = set_paths(os.path.join(sys.path[0], 'HaploSync.conf.toml'))

	iden_threshold = float(options.gene_identity)
	cov_threshold = float(options.gene_coverage)

	hit_iden = options.hit_identity
	hit_len = options.hit_length
	
	unbalanced_ratio = options.unbalanced_ratio 

	print('[' + str(datetime.datetime.now()) + '] = Read inputs', file=sys.stdout)
	print('# Read inputs', file=sys.stderr)
	# Make directory

	coord_tables = {}
	plot_files = {}
	print('[' + str(datetime.datetime.now()) + "] = Performing HaploDup", file=sys.stdout)
	mkdir( "./" + haplodup_dir )

	print('[' + str(datetime.datetime.now()) + '] == Loading FASTA sequences', file=sys.stdout)
	print('## Loading FASTA sequences', file=sys.stderr)
	fasta_files = options.fasta
	fasta_files_list = fasta_files.split(",")
	fasta_dict = {}
	for file_name in fasta_files_list :
		print('### Loading FASTA file: ' + file_name, file=sys.stderr)
		if fasta_dict == {} :
			fasta_dict = read_fasta(file_name)
		else :
			fasta_dict.update(read_fasta(file_name))
	fasta_len_dict = get_length_from_fasta_db(fasta_dict)
	complete_fasta_file_name = write_fasta_from_db( fasta_dict , haplodup_dir + "/" + options.out + ".input.fasta" , False )

	if options.reference :
		reference = read_fasta(options.reference)
		reference_len = get_length_from_fasta_db(reference)
		reference_file_name = write_fasta_from_db( reference , haplodup_dir + "/" + options.out + ".ref.fasta" , False)

	print('[' + str(datetime.datetime.now()) + '] == Loading sequence correspondence information', file=sys.stdout)
	print('## Loading sequence correspondence information', file=sys.stderr)
	corr_files = options.corr
	corr_files_list = corr_files.split(",")
	pairs = []

	hap1_fasta = {}
	hap2_fasta = {}

	hap1_ids = []
	hap2_ids = []
	ref_ids = []
	ref_to_hap1 = {}
	hap1_to_ref = {}
	ref_to_hap2 = {}
	hap2_to_ref = {}
	hap2_to_hap1 = {}
	hap1_to_hap2 = {}

	chr_ids = []
	chr_to_hap1 = {}
	chr_to_hap2 = {}
	chr_to_ref = {}
	hap1_to_chr = {}
	hap2_to_chr = {}
	ref_to_chr = {}

	for file_name in corr_files_list :
		print('### Loading correspondence file: ' + file_name, file=sys.stderr)
		if not options.reference :
			for line in open(file_name) :
				try :
					chr, seq1 , seq2 = line.rstrip().split("\t")
				except :
					print('[ERROR] Quitting', file=sys.stdout)
					print('[ERROR] Correspondence file  ' + file_name + ' has an unexpected number of columns (expected 3, in order: Chr_id, Hap1_id and Hap2_id)', file=sys.stderr)
					sys.exit(1)
				else :
					if seq1 not in fasta_dict or seq2 not in fasta_dict :
						print('[ERROR] Quitting', file=sys.stdout)
						print('[ERROR] Correspondence file  ' + file_name + " contains sequences missing from the FASTA file:" + seq1, file=sys.stderr)
						sys.exit(1)
					if seq2 not in fasta_dict :
						print('[ERROR] Quitting', file=sys.stdout)
						print('[ERROR] Correspondence file  ' + file_name + " contains sequences missing from the FASTA file:" + seq2, file=sys.stderr)
						sys.exit(1)
					hap1_ids.append(seq1)
					hap1_fasta[seq1] = fasta_dict[seq1]
					hap2_ids.append(seq2)
					hap2_fasta[seq2] = fasta_dict[seq2]
					pairs.append([seq1 , seq2])
					hap1_to_hap2[seq1] = seq2
					hap2_to_hap1[seq2] = seq1
					chr_ids.append(chr)
					chr_to_hap1[chr] = seq1
					chr_to_hap2[chr] = seq2
					hap1_to_chr[seq1] = chr
					hap2_to_chr[seq2] = chr
		else :
			for line in open(file_name) :
				try :
					chr , seq1 , seq2 , ref = line.rstrip().split("\t")
				except :
					print('[ERROR] Quitting', file=sys.stdout)
					print('[ERROR] Correspondence file  ' + file_name + ' has an unexpected number of columns (expected 4, in order: Chr_id, Hap1_id, Hap2_id, Ref_id)', file=sys.stderr)
					sys.exit(1)
				else :
					if seq1 not in fasta_dict or seq2 not in fasta_dict :
						print('[ERROR] Quitting', file=sys.stdout)
						print('[ERROR] Correspondence file  ' + file_name + " contains sequences missing from the FASTA file:" + seq1, file=sys.stderr)
						sys.exit(1)
					if seq2 not in fasta_dict :
						print('[ERROR] Quitting', file=sys.stdout)
						print('[ERROR] Correspondence file  ' + file_name + " contains sequences missing from the FASTA file:" + seq2, file=sys.stderr)
						sys.exit(1)
					hap1_ids.append(seq1)
					hap1_fasta[seq1] = fasta_dict[seq1]
					hap2_ids.append(seq2)
					hap2_fasta[seq2] = fasta_dict[seq2]
					ref_ids.append(ref)
					pairs.append([seq1 , seq2])
					hap1_to_hap2[seq1] = seq2
					hap2_to_hap1[seq2] = seq1
					hap1_to_ref[seq1] = ref
					hap2_to_ref[seq2] = ref
					ref_to_hap2[ref] = seq2
					ref_to_hap1[ref] = seq1
					chr_ids.append(chr)
					chr_to_hap1[chr] = seq1
					chr_to_hap2[chr] = seq2
					chr_to_ref[chr] = ref
					hap1_to_chr[seq1] = chr
					hap2_to_chr[seq2] = chr
					ref_to_chr[ref] = chr

	# Generate separate databases for the 2 haplotypes
	fasta_db_1 = {}
	fasta_db_2 = {}
	fasta_1_len = {}
	fasta_2_len = {}
	for hap1_seq_id in sorted(hap1_to_hap2.keys()) :
		hap2_seq_id = hap1_to_hap2[hap1_seq_id]
		fasta_db_1[hap1_seq_id] = fasta_dict[hap1_seq_id]
		fasta_1_len[hap1_seq_id] = len(fasta_dict[hap1_seq_id])
		fasta_db_2[hap2_seq_id] = fasta_dict[hap2_seq_id]
		fasta_2_len[hap2_seq_id] = len(fasta_dict[hap2_seq_id])

	# Read GFF files
	# For now the the annotation in GFF is mandatory
	# TODO: Use CDS sequences to build a temporary reference annotation

	annotation_db = {}

	if options.gff :
		if options.feature.upper() == "CDS" :
			feat_to_extract = "CDS"
		elif options.feature.upper() == "MRNA" :
			feat_to_extract = "mRNA"
		else :
			print('[ERROR] Unknown mapping feature ('+options.feature+') requested', file=sys.stdout)
			print('[ERROR] Unknown mapping feature ('+options.feature+') requested', file=sys.stderr)
			sys.exit(1)

		print('[' + str(datetime.datetime.now()) + '] == Loading GFF3 annotation', file=sys.stdout)
		print('## Loading GFF3 annotation', file=sys.stderr)
		gff_files = options.gff
		gff_files_list = gff_files.split(",")
		gff_db = {}
		mRNA_db = {}
		for file_name in gff_files_list :
			print('### Loading GFF3 file: ' + file_name, file=sys.stderr)
			if gff_db == {} :
				gff_db , mRNA_db = read_gff3(file_name)
			else :
				gff_tmp , mrna_tmp = read_gff3(file_name)
				gff_db.update(gff_tmp)
				mRNA_db.update(mrna_tmp)
		#else :
		#	# options.cds has been set instead
		#	# TODO: make one loci annotation if CDS sequences are given
		mRNA_to_gene_db = get_gene2mRNA_from_db(gff_db)
		# mRNA_to_gene_db[mRNA_id] = gene_id

		if options.annotation :
			# Read the functional annotation associated to the GFF file
			print('[' + str(datetime.datetime.now()) + '] == Loading functional annotation', file=sys.stdout)
			print('## Loading functional annotation', file=sys.stderr)
			annotation_files = options.annotation
			annotation_files_list = annotation_files.split(",")
			for file_name in annotation_files_list :
				print('### Loading annotation file: ' + file_name, file=sys.stderr)
				for line in open(file_name) :
					if line == "" or line[0] == "#" :
						continue
					else :
						feat_name , description = line.rstrip().split("\t")[0:2]
						description = description.upper()
						if feat_name in mRNA_to_gene_db :
							# Annotation associated to mRNAs -> add to gene annotation (uniquely)
							gene_id =  mRNA_to_gene_db[feat_name]
							if not gene_id in annotation_db :
								annotation_db[gene_id] = {}
							if description not in annotation_db[gene_id] :
								annotation_db[gene_id][description] = []
						else :
							if not feat_name in annotation_db :
								annotation_db[feat_name] = {}
							if description not in annotation_db[feat_name] :
								annotation_db[feat_name][description] = []

			annotation_db = add_counts_to_dict(annotation_db)

	# Read marker BED input
	if options.markers_hits :
		print('[' + str(datetime.datetime.now()) + '] == Loading markers position', file=sys.stdout)
		print('## Loading markers position', file=sys.stderr)
		# Copy the file in haplodup_dir
		all_markers_file_name = "all_markers.bed"
		dup_markers_file_name = "duplicated_markers.bed"
		all_markers_file_name_fullpath = haplodup_dir + "/all_markers.bed"
		dup_markers_file_name_fullpath = haplodup_dir + "/duplicated_markers.bed"
		seq_all_markers_file = open(all_markers_file_name_fullpath , 'w')
		seq_duplicated_markers_file = open(dup_markers_file_name_fullpath , 'w')
		## Write markers files
		markers_db = {}
		for line in open(options.markers_hits) :
			chr_id , start , stop , marker_id = line.rstrip().split("\t")
			if chr_id not in markers_db :
				markers_db[chr_id] = {}
			if marker_id not in markers_db[chr_id] :
				markers_db[chr_id][marker_id] = []
			markers_db[chr_id][marker_id].append([chr_id , start , stop , marker_id])

		# Read translated marker coordinates
		for chr_id in sorted(markers_db.keys()) :
			for marker_id in sorted(markers_db[chr_id].keys()) :
				if len(markers_db[chr_id][marker_id]) == 1 :
					# Unique hit
					print("\t".join(markers_db[chr_id][marker_id][0]), file=seq_all_markers_file)
				else :
					# Multiple hits, report all in both files
					for hit in sorted(markers_db[chr_id][marker_id]) :
						print("\t".join(hit), file=seq_all_markers_file)
						print("\t".join(hit), file=seq_duplicated_markers_file)
		seq_all_markers_file.close()
		seq_duplicated_markers_file.close()
	else :
		all_markers_file_name = ""
		dup_markers_file_name = ""


	# Pseudomolecule structure AGP
	if options.agp :
		agp_db = read_agp(options.agp)
		# Convert agp structure in table
		structure_file = "structure.tsv"
		structure_file_fullpath = haplodup_dir + "/structure.tsv"
		agp_table = []
		for seq_id in agp_db :
			seq_agp = agp_db[seq_id]
			for start in sorted(seq_agp.keys()) :
				Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
				if Compnt_Type == "W" :
					agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation ])
		structure_file_fullpath = write_table(agp_table, structure_file_fullpath)
	else :
		structure_file = ""
		structure_file_fullpath = ""

	# Pseudomolecule legacy structure AGP
	if options.legacy_agp :
		legacy_agp = read_agp(options.legacy_agp)
		# convert the db in table and write it down
		legacy_structure_file = "legacy_structure.tsv"
		legacy_structure_file_fullpath = haplodup_dir + "/legacy_structure.tsv"
		legacy_agp_table = []
		for seq_id in legacy_agp :
			seq_agp = legacy_agp[seq_id]
			for start in sorted(seq_agp.keys()) :
				Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation = seq_agp[start]
				if Compnt_Type == "W" :
					legacy_agp_table.append([Obj_Name , Obj_start , Obj_End , CompntId , Orientation ])
		legacy_structure_file_fullpath = write_table(legacy_agp_table, legacy_structure_file_fullpath)
	else :
		legacy_structure_file = ""
		legacy_structure_file_fullpath = ""

	showcoords_path = paths["show-coords"]
	if showcoords_path == "" :
		minimap2_search=subprocess.Popen( "which show-coords" , shell=True, stdout=subprocess.PIPE, text=True)
		command_line , error = minimap2_search.communicate()
		command_line = command_line.rstrip()
	else :
		command_line = showcoords_path + "/show-coords"
	if not os.path.exists(command_line) :
		print("[ERROR] wrong or no path to show-coords (Mummer4)", file=sys.stderr)
		sys.exit(1)

	# Perform nucmer alignments
	hap1_ids = ",".join(sorted(fasta_db_1.keys()))
	if options.reference :
		reference_ids = ",".join(sorted(reference.keys()))
	hap2_ids = ",".join(sorted(fasta_db_2.keys()))

	print('[' + str(datetime.datetime.now()) + "] = Mapping sequences", file=sys.stdout)
	query_1_file = haplodup_dir + "/" + options.out + ".1.fasta"
	write_fasta_from_db( fasta_db_1 , query_1_file )
	query_2_file = haplodup_dir + "/" + options.out + ".2.fasta"
	write_fasta_from_db( fasta_db_2 , query_2_file )

	print('[' + str(datetime.datetime.now()) + "] == Hap1 vs Hap1", file=sys.stdout)
	print('[' + str(datetime.datetime.now()) + "] === Mapping", file=sys.stdout)
	if options.reuse_mappings :
		print('[' + str(datetime.datetime.now()) + "] === [SKIP] Mapping: --reuse_mappings set, reusing existing files", file=sys.stdout)
		outfile_prefix = "Hap1.on.Hap1"
	else :
		outfile_prefix = map_nucmer_dotplot("Hap1" , query_1_file , "Hap1" , query_1_file , haplodup_dir , options.cores , paths , False )
	print('[' + str(datetime.datetime.now()) + "] === Converting files", file=sys.stdout)
	coord_tables["Hap1_vs_Hap1"] = make_coords_table( outfile_prefix , haplodup_dir , command_line)

	print('[' + str(datetime.datetime.now()) + "] == Hap2 vs Hap2", file=sys.stdout)
	print('[' + str(datetime.datetime.now()) + "] === Mapping", file=sys.stdout)
	if options.reuse_mappings :
		print('[' + str(datetime.datetime.now()) + "] === [SKIP] Mapping: --reuse_mappings set, reusing existing files", file=sys.stdout)
		outfile_prefix = "Hap2.on.Hap2"
	else:
		outfile_prefix = map_nucmer_dotplot("Hap2" , query_2_file , "Hap2" , query_2_file , haplodup_dir , options.cores , paths , False  )
	print('[' + str(datetime.datetime.now()) + "] === Converting files", file=sys.stdout)
	coord_tables["Hap2_vs_Hap2"] = make_coords_table( outfile_prefix , haplodup_dir , command_line)

	print('[' + str(datetime.datetime.now()) + "] == Hap2 vs Hap1", file=sys.stdout)
	print('[' + str(datetime.datetime.now()) + "] === Mapping", file=sys.stdout)
	if options.reuse_mappings :
		print('[' + str(datetime.datetime.now()) + "] === [SKIP] Mapping: --reuse_mappings set, reusing existing files", file=sys.stdout)
		outfile_prefix = "Hap2.on.Hap1"
	else :
		outfile_prefix = map_nucmer_dotplot("Hap1" , query_1_file , "Hap2" , query_2_file , haplodup_dir , options.cores , paths , False  )
	print('[' + str(datetime.datetime.now()) + "] === Converting files", file=sys.stdout)
	coord_tables["Hap2_vs_Hap1"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap1_vs_Hap1"] ]

	print('[' + str(datetime.datetime.now()) + "] == Hap1 vs Hap2", file=sys.stdout)
	print('[' + str(datetime.datetime.now()) + "] === Mapping", file=sys.stdout)
	if options.reuse_mappings :
		print('[' + str(datetime.datetime.now()) + "] === [SKIP] Mapping: --reuse_mappings set, reusing existing files", file=sys.stdout)
		outfile_prefix = "Hap1.on.Hap2"
	else :
		outfile_prefix = map_nucmer_dotplot("Hap2" , query_2_file , "Hap1" , query_1_file , haplodup_dir , options.cores , paths , False  )
	print('[' + str(datetime.datetime.now()) + "] === Converting files", file=sys.stdout)
	coord_tables["Hap1_vs_Hap2"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap2_vs_Hap2"] ]


	if options.reference :
		print('[' + str(datetime.datetime.now()) + "] == Hap1 vs Reference", file=sys.stdout)
		print('[' + str(datetime.datetime.now()) + "] === Mapping", file=sys.stdout)
		if options.reuse_mappings :
			print('[' + str(datetime.datetime.now()) + "] === [SKIP] Mapping: --reuse_mappings set, reusing existing files", file=sys.stdout)
			outfile_prefix = "Hap1.on.Ref"
		else :
			outfile_prefix = map_nucmer_dotplot("Ref" , options.reference , "Hap1" , query_1_file , haplodup_dir , options.cores , paths , False )
			# outfile_prefix.delta and outfile_prefix.coords (show-cords -c) are generated
		print('[' + str(datetime.datetime.now()) + "] === Converting files", file=sys.stdout)
		coord_tables["Hap1_vs_Reference"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , "" ]

		print('[' + str(datetime.datetime.now()) + "] == Reference vs Hap1 ", file=sys.stdout)
		print('[' + str(datetime.datetime.now()) + "] === Mapping", file=sys.stdout)
		if options.reuse_mappings :
			print('[' + str(datetime.datetime.now()) + "] === [SKIP] Mapping: --reuse_mappings set, reusing existing files", file=sys.stdout)
			outfile_prefix = "Ref.on.Hap1"
		else:
			outfile_prefix = map_nucmer_dotplot( "Hap1" , query_1_file , "Ref" , options.reference , haplodup_dir , options.cores , paths , False )
		print('[' + str(datetime.datetime.now()) + "] === Converting files", file=sys.stdout)
		coord_tables["Reference_vs_Hap1"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap1_vs_Hap1"] ]

		print('[' + str(datetime.datetime.now()) + "] == Hap2 vs Reference", file=sys.stdout)
		print('[' + str(datetime.datetime.now()) + "] === Mapping", file=sys.stdout)
		if options.reuse_mappings :
			print('[' + str(datetime.datetime.now()) + "] === [SKIP] Mapping: --reuse_mappings set, reusing existing files", file=sys.stdout)
			outfile_prefix = "Hap2.on.Ref"
		else :
			outfile_prefix = map_nucmer_dotplot("Ref" , options.reference , "Hap2" , query_2_file , haplodup_dir , options.cores , paths , False )
		print('[' + str(datetime.datetime.now()) + "] === Converting files", file=sys.stdout)
		coord_tables["Hap2_vs_Reference"] = [ make_coords_table( outfile_prefix , haplodup_dir , command_line) , "" ]

		print('[' + str(datetime.datetime.now()) + "] == Reference vs Hap2", file=sys.stdout)
		print('[' + str(datetime.datetime.now()) + "] === Mapping", file=sys.stdout)
		if options.reuse_mappings :
			print('[' + str(datetime.datetime.now()) + "] === [SKIP] Mapping: --reuse_mappings set, reusing existing files", file=sys.stdout)
			outfile_prefix = "Ref.on.Hap2"
		else :
			outfile_prefix = map_nucmer_dotplot("Hap2" , query_2_file , "Ref" , options.reference , haplodup_dir , options.cores , paths , False )
		print('[' + str(datetime.datetime.now()) + "] === Converting files", file=sys.stdout)
		coord_tables["Reference_vs_Hap2"] = [make_coords_table( outfile_prefix , haplodup_dir , command_line) , coord_tables["Hap2_vs_Hap2"] ]


	if options.gff :
		gmap_results = haplodup_dir + "/CDS.on.genome.gmap.gff3"
		# Map genes, generate copy number counts, render by chromosome reports
		print('[' + str(datetime.datetime.now()) + "] = Generating files about gene map count for fusion dedup", file=sys.stdout)

		## Generate CDS sequences from first haplotype
		print('[' + str(datetime.datetime.now()) + "] == Generating CDS sequences", file=sys.stdout)
		if options.reuse_gmap :
			print('[' + str(datetime.datetime.now()) + "] == [SKIP] CDS extraction and GMAP mapping: --reuse_gmap set, reusing " + haplodup_dir + "/CDS.on.genome.gmap.gff3", file=sys.stdout)
			gmap_results = haplodup_dir + "/CDS.on.genome.gmap.gff3"

		else :
			print('[' + str(datetime.datetime.now()) + "] == Mapping CDSs on genome assembly", file=sys.stdout)
			CDS_file = get_sequence( gff_db , fasta_dict , haplodup_dir + "/new" , "CDS")
			index_dir = haplodup_dir + "/gmap_index"

			# The new gmap

			mkdir(index_dir)

			# Hap1
			## index
			print('[' + str(datetime.datetime.now()) + "] === Hap1 ", file=sys.stdout)
			print('[' + str(datetime.datetime.now()) + "] ==== Indexing genome ", file=sys.stdout)
			indexing_out_file = open( haplodup_dir + "/gmap_index.1.log" ,"w" )
			indexing_err_file = open( haplodup_dir + "/gmap_index.1.err" ,"w" )
			hap1_name = haplodup_dir + "/new.hap1.fasta"
			hap1_file = write_fasta_from_db( hap1_fasta , hap1_name)
			indexing_command = "gmap_build -D " + index_dir + " -d hap1.fasta " + hap1_name
			indexProcess = subprocess.Popen(indexing_command, shell=True, stdout=indexing_out_file , stderr=indexing_err_file)
			output, error = indexProcess.communicate()
			indexing_out_file.close()
			indexing_err_file.close()

			# Gmap CDS on results
			print('[' + str(datetime.datetime.now()) + "] ==== Mapping with Gmap ", file=sys.stdout)
			gmap_results_1 = haplodup_dir + "/CDS.on.hap1.gmap.gff3"
			gmap_1_gff3 = open( gmap_results_1 , "w" )
			gmap_err = open( gmap_results_1 + ".err" , "w" )
			gmapCommand = "gmap -D " + index_dir + " -d hap1.fasta -f 2 -n 500 -t " + str(options.cores) + " " + CDS_file

			gmapProcess = subprocess.Popen(gmapCommand, shell=True, stdout=gmap_1_gff3 , stderr=gmap_err)
			output, error = gmapProcess.communicate()
			gmap_1_gff3.close()
			gmap_err.close()


			# Hap2
			## Index
			print('[' + str(datetime.datetime.now()) + "] === Hap2 ", file=sys.stdout)
			print('[' + str(datetime.datetime.now()) + "] ==== Indexing genome ", file=sys.stdout)
			indexing_out_file = open( haplodup_dir + "/gmap_index.2.log" ,"w" )
			indexing_err_file = open( haplodup_dir + "/gmap_index.2.err" ,"w" )
			hap2_name = haplodup_dir + "/new.hap2.fasta"
			hap2_file = write_fasta_from_db( hap2_fasta , hap2_name)
			indexing_command = "gmap_build -D " + index_dir + " -d hap2.fasta " + hap2_name
			indexProcess = subprocess.Popen(indexing_command, shell=True, stdout=indexing_out_file , stderr=indexing_err_file)
			output, error = indexProcess.communicate()
			indexing_out_file.close()
			indexing_err_file.close()
			# Gmap CDS on results
			print('[' + str(datetime.datetime.now()) + "] ==== Mapping with Gmap ", file=sys.stdout)
			gmap_results_2 = haplodup_dir + "/CDS.on.hap2.gmap.gff3"
			gmap_2_gff3 = open( gmap_results_2 , "w" )
			gmap_err = open( gmap_results_2 + ".err" , "w" )
			gmapCommand = "gmap -D " + index_dir + " -d hap2.fasta -f 2 -n 500 -t " + str(options.cores) + " " + CDS_file

			gmapProcess = subprocess.Popen(gmapCommand, shell=True, stdout=gmap_2_gff3 , stderr=gmap_err)
			output, error = gmapProcess.communicate()
			gmap_2_gff3.close()
			gmap_err.close()



			# Concatenate Hap1 and Hap2 results
			gmap_results = haplodup_dir + "/CDS.on.genome.gmap.gff3"
			gmap_gff3 = open( gmap_results , "w" )
			for line in open( gmap_results_1 ) :
				print(line.rstrip(), file=gmap_gff3)
			for line in open( gmap_results_2 ) :
				print(line.rstrip(), file=gmap_gff3)
			gmap_gff3.close()

		# Extract valid alignments per locus
		print('[' + str(datetime.datetime.now()) + "] === Extracting valid alignments (identity > " + str(iden_threshold) + "% , coverage> " + str(cov_threshold)  + "%)", file=sys.stdout)
		gmap_hits_hap1 , gmap_hits_hap2 = read_gmap_results_Hap(gmap_results, cov_threshold , iden_threshold , "by_locus", mRNA_to_gene_db)
		test_1 = open( haplodup_dir + "/gmap_hits_hap1.txt" ,'w')
		for chr in sorted(gmap_hits_hap1.keys()) :
			for locus in sorted(gmap_hits_hap1[chr].keys()) :
				for hit in sorted(gmap_hits_hap1[chr][locus]) :
					print(chr + "\t" + "\t".join([ str(x) for x in hit ]) + "\t" + locus, file=test_1)
		test_1.close()
		## Make Hap1 gene table
		print('[' + str(datetime.datetime.now()) + "] === Extracting loci positions on Hap1", file=sys.stdout)
		# TODO: correct filtering procedure
		hap1_genes = gff3_filter2table_Hap(gff_db, "gene", "hap1")
		#print >> sys.stderr, "### hap1 gene"
		#print >> sys.stderr, hap1_genes

		test_2 = open( haplodup_dir + "/gmap_hits_hap2.txt" ,'w')
		for chr in sorted(gmap_hits_hap2.keys()) :
			for locus in sorted(gmap_hits_hap2[chr].keys()) :
				for hit in sorted(gmap_hits_hap2[chr][locus]) :
					print(chr + "\t" + "\t".join([ str(x) for x in hit ]) + "\t" + locus, file=test_2)
		test_2.close()
		print('[' + str(datetime.datetime.now()) + "] === Extracting loci positions on Hap2", file=sys.stdout)
		hap2_genes = gff3_filter2table_Hap(gff_db, "gene", "hap2")
		#print >> sys.stderr, "### hap2 gene"
		#print >> sys.stderr, hap2_genes

		# Decide hap tokens automatically from GFF seqnames (no argparse required)
		hap1_token, hap2_token = "hap1", "hap2"
		try:
			seqnames = set()
			for entry in list(gff_db.values()):
				try:
					seqnames.add(entry[1])
				except Exception:
					continue
			# prefer exact dot style
			if any('.hap1.' in s.lower() for s in seqnames):
				hap1_token, hap2_token = ".hap1.", ".hap2."
			# prefer underscore style (preserve case if present)
			elif any('_hap1_' in s.lower() for s in seqnames):
				if any('_Hap1_' in s for s in seqnames):
					hap1_token, hap2_token = "_Hap1_", "_Hap2_"
				else:
					hap1_token, hap2_token = "_hap1_", "_hap2_"
		except Exception:
			pass

		# Join results
		print('[' + str(datetime.datetime.now()) + "] === Counting intra-chromosome hits", file=sys.stdout)
		print('[' + str(datetime.datetime.now()) + "] ==== Hap1", file=sys.stdout)
		hit_counts_1 = do_count_hits_Hap(hap1_genes, gmap_hits_hap1, gmap_hits_hap2, hap1_token, {})
		hit_file_1 = print_hit_counts(hit_counts_1, haplodup_dir + "/diploid_gene_count_trace.hap1.txt")

		print('[' + str(datetime.datetime.now()) + "] ==== Hap2", file=sys.stdout)
		hit_counts_2 = do_count_hits_Hap(hap2_genes, gmap_hits_hap1, gmap_hits_hap2, hap2_token, {})
		hit_file_2 = print_hit_counts(hit_counts_2, haplodup_dir + "/diploid_gene_count_trace.hap2.txt")



	if options.rejected :
		# Read rejection relationships
		print('[' + str(datetime.datetime.now()) + '] = Running unused sequences QC analysis', file=sys.stdout)
		print('# Running unused sequences QC analysis', file=sys.stderr)

		structure_comparison_dir = options.out + ".structure_comparison"
		mkdir(structure_comparison_dir)

		all_unused_by_chr = {}
		all_unused_to_chr = {}
		#for seq_id in all_unused_by_seq_id.keys() :
		#	chr_id , orientation = all_unused_by_seq_id[seq_id]
		#	if chr_id not in all_unused_by_chr :
		#		all_unused_by_chr[chr_id] = []
		#	all_unused_by_chr[chr_id].append(seq_id+"|"+orientation)
		for line in open(options.rejected , 'r') :
			seq_id , hapx_id , orientation = line.rstrip().split("\t")
			if hapx_id in hap1_to_chr :
				chr_id = hap1_to_chr[hapx_id]
			elif hapx_id in hap2_to_chr :
				chr_id = hap2_to_chr[hapx_id]
			elif hapx_id in ref_to_chr :
				chr_id = ref_to_chr[hapx_id]
			elif hapx_id in chr_to_hap1 :
				chr_id = hapx_id
			else :
				print("[WARNING] Sequence " + seq_id + " is associated to " + hapx_id + " pseudomolecule that is not present in correspondence file. Analysis cannot be performed", file=sys.stderr)
				continue

			if chr_id not in all_unused_by_chr :
				all_unused_by_chr[chr_id] = []
			all_unused_by_chr[chr_id].append(seq_id+"|"+orientation)
			all_unused_to_chr[seq_id] = chr_id
			all_unused_to_chr[seq_id+"|"+orientation] = chr_id

		if options.marker_map and options.markers_hits :
			marker_map_by_seq = {}
			marker_map_by_id = {}
			# Read map
			for line in open( options.marker_map ) :
				if line == "" :
					continue
				if line[0] == "#" :
					continue
				try :
					seq_id , pos , marker_id = line.rstrip().split("\t")
				except :
					print("Error in markers map file " + options.marker_map + ", line do not match file format: " + line.rstrip(), file=sys.stdout)
					print("Error in markers map file " + options.marker_map + ", line do not match file format: " + line.rstrip(), file=sys.stdout)
					sys.exit(1)
				else :
					if seq_id not in marker_map_by_seq :
						marker_map_by_seq[seq_id] = []
					marker_map_by_seq[seq_id].append([ int(pos) , marker_id ] )
					marker_map_by_id[marker_id] = [ seq_id , int(pos) , marker_id ]
		else :
			marker_map_by_seq = ""
			print("[WARNING] Marker information missing", file=sys.stderr)

		if options.legacy_groups :
			associated_legacy_seqid_file = structure_comparison_dir + "/legacy_components.association.tsv"
			# all_agp_db = legacy_agp + agp_db >> all sequences related to the components >> grouping on components
			all_agp_db = dict(legacy_agp)
			all_agp_db.update(agp_db)
			associated_input_seqid_file , seq_group_db = make_seq_pair_from_groups( associated_legacy_seqid_file , options.legacy_groups , all_agp_db , hap1_ids.split(",") , hap2_ids.split(",") , list(fasta_dict.keys()) , {} , "legacy" )
			# associated_input_seqid_file >> [ ... , [ "hap1_legacy" , Tid , Tstart , Tstop , "Query_legacy" , Qid , Qstart , Qstop , group_id] , ... ]
			#
		else :
			# No legacy info available -> use what is known the input sequences
			associated_input_seqid_file = structure_comparison_dir + "/input.association.tsv"
			if options.input_groups :
				# override with user option options.input_groups
				if options.input_groups == 0 :
					# force mapping,
					associated_input_seqid_file = ""
				else :
					# file name is given -> read as [ id , group ] table
					# agp_db + feed all input sequences length >> whole sequence relationship for unplaced >> grouping on input sequences
					associated_input_seqid_file , seq_group_db = make_seq_pair_from_groups( associated_input_seqid_file ,options.input_groups , agp_db , hap1_ids.split(",") , hap2_ids.split(",") , list(fasta_dict.keys()) , fasta_len_dict , "input" )
					# associated_input_seqid_file >> [ ... , [ "hap1_HS" , Tid , Tstart , Tstop , "Query_HS" , Qid , Qstart , Qstop , group_id] , ... ]
			#else :
			#	# parse known_groups and unwanted_pairs
			#	associated_input_seqid_file = make_seq_pair_from_constrains(associated_input_seqid_file, known_groups, unwanted_pairs, alternative_sequences , agp_db, hap1_ids.split(","), hap2_ids.split(","), query_fasta_db.keys(), query_len, "input")
			#	# associated_input_seqid_file >> [ ... , [ "hap1_HS" , Tid , Tstart , Tstop , "Query_HS" , Qid , Qstart , Qstop , group_id] , ... ]

		reported_hits_on_seq = report_marker_usage( options.markers_hits , marker_map_by_seq , marker_map_by_id , agp_db , legacy_agp , hap1_to_chr , hap2_to_chr , all_unused_to_chr , structure_comparison_dir )
		# 	reported_hits_on_seq[seq_id]["id"] : seq_id
		# 	reported_hits_on_seq[seq_id]["chr"] : chr_id
		# 	reported_hits_on_seq[seq_id]["markers"] : [ ... , [chr_id , chr_pos , marker_id , seq_id, start , stop] , ... ]
		# 	reported_hits_on_seq[seq_id]["range"] : [marker_pos_min , marker_pos_max] ,
		# 	reported_hits_on_seq[seq_id]["orientation"] : ["+" or "-" or "."] }

		if options.markers_hits :
			# --markers_hits is typically run.markers.bed (translated positions), which has wrapper IDs
			# like "Hap1_unplaced_ctg123" but NOT the bare contig ID "ctg123" that rejected_QC needs.
			# Build a mapping: bare_contig_id -> (wrapper_id, agp_entry) from agp_db
			bare_to_wrapper = {}
			for wrapper_id in agp_db :
				for agp_start in sorted(agp_db[wrapper_id].keys()) :
					entry = agp_db[wrapper_id][agp_start]
					if entry[4] == "W" :
						bare_to_wrapper[entry[5]] = ( wrapper_id , entry )
			# For each unplaced sequence missing from markers_db, translate from its wrapper
			for unpl_id in sorted(all_unused_to_chr.keys()) :
				if "|" in unpl_id or unpl_id in markers_db :
					continue
				if unpl_id not in bare_to_wrapper :
					continue
				wrapper_id , entry = bare_to_wrapper[unpl_id]
				if wrapper_id not in markers_db :
					continue
				obj_start    = int(entry[1])
				compnt_start = int(entry[6])
				compnt_end   = int(entry[7])
				orient       = entry[8]
				markers_db[unpl_id] = {}
				for mk_id in markers_db[wrapper_id] :
					markers_db[unpl_id][mk_id] = []
					for hit in markers_db[wrapper_id][mk_id] :
						_ , h_start , h_stop , _ = hit
						hs , he = int(h_start) , int(h_stop)
						if orient == "+" :
							off = compnt_start - obj_start
							ns , ne = hs + off , he + off
						else :
							ns = compnt_end - (he - obj_start)
							ne = compnt_end - (hs - obj_start)
						markers_db[unpl_id][mk_id].append([ unpl_id , str(ns) , str(ne) , mk_id ])
			# Populate reported_hits_on_seq for unplaced sequences still missing (needed for usage flags).
			if marker_map_by_id :
				for unpl_id in sorted(all_unused_to_chr.keys()) :
					if "|" in unpl_id or unpl_id in reported_hits_on_seq :
						continue
					if unpl_id not in markers_db :
						continue
					unpl_chr = all_unused_to_chr[unpl_id]
					hits = []
					for mk_id in markers_db[unpl_id] :
						if mk_id not in marker_map_by_id :
							continue
						marker_chr , marker_pos , _ = marker_map_by_id[mk_id]
						for hit in markers_db[unpl_id][mk_id] :
							_ , h_start , h_stop , _ = hit
							hits.append([ marker_chr , marker_pos , mk_id , unpl_id , int(h_start) , int(h_stop) ])
					if hits :
						marker_pos_values = [ int(x[1]) for x in hits ]
						reported_hits_on_seq[unpl_id] = {
							"id"          : unpl_id ,
							"chr"         : unpl_chr ,
							"orientation" : "+" ,
							"markers"     : hits ,
							"range"       : [ min(marker_pos_values) , max(marker_pos_values) ]
						}
		print('[' + str(datetime.datetime.now()) + "] = Reading markers hits", file=sys.stdout)
		structure_plot_db = {"Rejected" : {}}

		# Pre-create output directories and collect all (chr_id, seq_id) work items
		work_items = []
		for chr_id in sorted(chr_ids) :
			print('[' + str(datetime.datetime.now()) + "] == Chr: " + chr_id, file=sys.stdout)
			print("## Chr: " + chr_id, file=sys.stderr)
			if not chr_id in all_unused_by_chr :
				continue
			qc_out_dir = structure_comparison_dir + "/" + chr_id
			mkdir(qc_out_dir)
			structure_plot_db["Rejected"][chr_id] = {}
			for seq_id in all_unused_by_chr[chr_id] :
				work_items.append((chr_id, seq_id))

		def _run_rejected_qc(chr_id, seq_id) :
			print("### seq_id: " + seq_id, file=sys.stderr)
			query_id_bare = seq_id[:-2]
			# Build a synthetic input_agp for this unplaced sequence.
			# In HaploSplit, input_agp = old_legacy_agp (original assembly AGP) which provides
			# entries keyed by the original sequence ID (e.g. "ctg123"), enabling the Rmd to
			# filter agp_legacy by unplacedID, with CompntId = legacy contig IDs (needed for
			# group color lookup via groups_by_sequence).
			# Here we reconstruct that by:
			#   1. Finding the wrapper sequence in agp_db that wraps query_id_bare
			#   2. Using the legacy_agp entry for that wrapper (which has legacy contig IDs as components)
			#   3. Remapping Obj_Name to query_id_bare so the Rmd filter matches unplacedID
			unpl_input_agp = {}
			wrapper_name = None
			for agp_key in agp_db :
				for agp_start in sorted(agp_db[agp_key].keys()) :
					entry = agp_db[agp_key][agp_start]
					if entry[4] == "W" and entry[5] == query_id_bare :
						wrapper_name = agp_key
						break
				if wrapper_name is not None :
					break
			if wrapper_name is not None and wrapper_name in legacy_agp :
				# Preferred: use legacy_agp structure (legacy contig IDs as components)
				for leg_start in sorted(legacy_agp[wrapper_name].keys()) :
					entry = legacy_agp[wrapper_name][leg_start]
					Obj_Name_e , Obj_start_e , Obj_End_e , PartNum_e , Compnt_Type_e , CompntId_e , CompntStart_e , CompntEnd_e , Orientation_e = entry
					if Compnt_Type_e == "W" :
						if query_id_bare not in unpl_input_agp :
							unpl_input_agp[query_id_bare] = {}
						unpl_input_agp[query_id_bare][int(Obj_start_e)] = [ query_id_bare , Obj_start_e , Obj_End_e , PartNum_e , "W" , CompntId_e , CompntStart_e , CompntEnd_e , Orientation_e ]
			elif wrapper_name is not None :
				# Fallback: wrapper not in legacy_agp, use self-referencing entry from agp_db
				for agp_start in sorted(agp_db[wrapper_name].keys()) :
					entry = agp_db[wrapper_name][agp_start]
					Obj_Name_e , Obj_start_e , Obj_End_e , PartNum_e , Compnt_Type_e , CompntId_e , CompntStart_e , CompntEnd_e , Orientation_e = entry
					if Compnt_Type_e == "W" and CompntId_e == query_id_bare :
						if query_id_bare not in unpl_input_agp :
							unpl_input_agp[query_id_bare] = {}
						unpl_input_agp[query_id_bare][int(CompntStart_e)] = [ query_id_bare , CompntStart_e , CompntEnd_e , PartNum_e , "W" , CompntId_e , CompntStart_e , CompntEnd_e , Orientation_e ]
			return chr_id , seq_id , rejected_QC( structure_comparison_dir , seq_id , fasta_dict , chr_id , fasta_db_1 , fasta_db_2 , chr_to_hap1 , chr_to_hap2 , coord_tables["Hap2_vs_Hap1"][0] , associated_input_seqid_file , associated_legacy_seqid_file , agp_db , legacy_agp , unpl_input_agp , seq_group_db , markers_db , reported_hits_on_seq , marker_map_by_seq , options.cores , paths)

		# Run QC reports in parallel — subprocess.communicate() releases the GIL so threads
		# genuinely overlap while waiting for Rscript processes to finish.
		with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.cores)) as executor :
			futures = { executor.submit(_run_rejected_qc, chr_id, seq_id) : (chr_id, seq_id) for chr_id, seq_id in work_items }
			for future in concurrent.futures.as_completed(futures) :
				chr_id , seq_id , result = future.result()
				structure_plot_db["Rejected"][chr_id][seq_id] = result

		rejected_index_file_name = "index.rejected_sequences.html"
		rejected_index_file_full_path = structure_comparison_dir + "/" + rejected_index_file_name
		rejected_index_file_full_path = make_index_from_report_db(rejected_index_file_name , "." , structure_comparison_dir ,  structure_plot_db  )


	if options.skip_chr_pair_reports :
		print('[' + str(datetime.datetime.now()) + '] = [SKIP] Chromosome pair overview reports: --skip_chr_pair_reports set', file=sys.stdout)
		print('# [SKIP] Chromosome pair overview reports: --skip_chr_pair_reports set', file=sys.stderr)
	else :
		print('[' + str(datetime.datetime.now()) + '] = Generating chromosome pair overview reports', file=sys.stdout)
		print('# Generating chromosome pair overview reports', file=sys.stderr)

		chr_pair_reports_dir = options.out + ".chr_pair_reports"
		mkdir(chr_pair_reports_dir)

		# Load marker_map if available and not already loaded (e.g. --rejected was not used)
		if options.marker_map and options.markers_hits :
			try :
				chr_pair_marker_map
			except NameError :
				chr_pair_marker_map = {}
				for line in open(options.marker_map) :
					if line == "" or line[0] == "#" :
						continue
					try :
						seq_id , pos , marker_id = line.rstrip().split("\t")
					except :
						pass
					else :
						if seq_id not in chr_pair_marker_map :
							chr_pair_marker_map[seq_id] = []
						chr_pair_marker_map[seq_id].append([int(pos) , marker_id])
		else :
			chr_pair_marker_map = {}

		# Use marker_map_by_seq if available from --rejected, otherwise use chr_pair_marker_map
		try :
			effective_marker_map = marker_map_by_seq
		except NameError :
			effective_marker_map = chr_pair_marker_map

		# Use markers_db if available, otherwise empty
		effective_marker_bed = markers_db if options.markers_hits else ""

		# Use seq_group_db if available from --rejected, otherwise rebuild from group file
		try :
			effective_groups = seq_group_db
		except NameError :
			group_file = options.legacy_groups if options.legacy_groups else (options.input_groups if options.input_groups and options.input_groups != 0 else None)
			if group_file :
				effective_groups = {}
				for seq_id , group_id in read_table(group_file) :
					if seq_id not in effective_groups :
						effective_groups[seq_id] = []
					effective_groups[seq_id].append(group_id)
			else :
				effective_groups = {}

		chr_pair_plot_db = {"Chr_Pair_Reports" : {"Reports" : {}}}

		def _run_chr_pair(chr_id) :
			chr_pair_out_dir = chr_pair_reports_dir + "/" + chr_id
			mkdir(chr_pair_out_dir)
			result = chr_pair_report(chr_pair_out_dir , chr_id , fasta_db_1 , fasta_db_2 , chr_to_hap1 , chr_to_hap2 , coord_tables["Hap2_vs_Hap1"][0] , agp_db , legacy_agp , effective_marker_bed , "" , effective_marker_map , effective_groups , options.cores , paths)
			result["html"] = chr_id + "/" + result["html"]
			return chr_id , result

		# Run chr pair reports in parallel — subprocess.communicate() releases the GIL so threads
		# genuinely overlap while waiting for Rscript processes to finish.
		with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.cores)) as executor :
			futures = { executor.submit(_run_chr_pair, chr_id) : chr_id for chr_id in sorted(chr_ids) }
			for future in concurrent.futures.as_completed(futures) :
				chr_id , result = future.result()
				print('[' + str(datetime.datetime.now()) + "] == Chr: " + chr_id, file=sys.stdout)
				print("## Chr: " + chr_id, file=sys.stderr)
				chr_pair_plot_db["Chr_Pair_Reports"]["Reports"][chr_id] = result

		chr_pair_index_file_name = "index.chr_pair_reports.html"
		make_index_from_report_db(chr_pair_index_file_name , "." , chr_pair_reports_dir , chr_pair_plot_db)


	print('[' + str(datetime.datetime.now()) + "] = Generating dotplots", file=sys.stdout)

	# Build paired-only maps for --only_paired_dotplots:
	# Each map is {tId: qId} covering only the matched pairs for that comparison.
	if options.only_paired_dotplots :
		_pm_h1h1 = {h: h for h in hap1_ids.split(",")}
		_pm_h2h2 = {h: h for h in hap2_ids.split(",")}
		_pm_h2h1 = dict(hap1_to_hap2)   # tIds=hap1, qIds=hap2
		_pm_h1h2 = dict(hap2_to_hap1)   # tIds=hap2, qIds=hap1
		if options.reference :
			_pm_refh1 = dict(ref_to_hap1)   # tIds=ref,  qIds=hap1
			_pm_h1ref = dict(hap1_to_ref)   # tIds=hap1, qIds=ref
			_pm_refh2 = dict(ref_to_hap2)   # tIds=ref,  qIds=hap2
			_pm_h2ref = dict(hap2_to_ref)   # tIds=hap2, qIds=ref
	else :
		_pm_h1h1 = _pm_h2h2 = _pm_h2h1 = _pm_h1h2 = None
		if options.reference :
			_pm_refh1 = _pm_h1ref = _pm_refh2 = _pm_h2ref = None

	print('[' + str(datetime.datetime.now()) + "] == Hap1 vs Hap1", file=sys.stdout)
	outfile_prefix = "Hap1.on.Hap1"
	plot_files["Hap1_vs_Hap1"] = {}
	plot_files["Hap1_vs_Hap1"]["Whole"] , plot_files["Hap1_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , hap1_ids , outfile_prefix, haplodup_dir, "Hap1_vs_Hap1", coord_tables["Hap1_vs_Hap1"] , options.reuse_dotplots , options.skip_dotplots , paired_only_map = _pm_h1h1)

	print('[' + str(datetime.datetime.now()) + "] == Hap2 vs Hap2", file=sys.stdout)
	outfile_prefix = "Hap2.on.Hap2"
	plot_files["Hap2_vs_Hap2"] = {}
	plot_files["Hap2_vs_Hap2"]["Whole"] , plot_files["Hap2_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids, hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Hap2", coord_tables["Hap2_vs_Hap2"] , options.reuse_dotplots , options.skip_dotplots , paired_only_map = _pm_h2h2)

	print('[' + str(datetime.datetime.now()) + "] == Hap2 vs Hap1", file=sys.stdout)
	outfile_prefix = "Hap2.on.Hap1"
	plot_files["Hap2_vs_Hap1"] = {}
	plot_files["Hap2_vs_Hap1"]["Whole"] , plot_files["Hap2_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Hap1", coord_tables["Hap2_vs_Hap1"][0] , options.reuse_dotplots , options.skip_dotplots , paired_only_map = _pm_h2h1)

	print('[' + str(datetime.datetime.now()) + "] == Hap1 vs Hap2", file=sys.stdout)
	outfile_prefix = "Hap1.on.Hap2"
	plot_files["Hap1_vs_Hap2"] = {}
	plot_files["Hap1_vs_Hap2"]["Whole"], plot_files["Hap1_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids , hap1_ids, outfile_prefix, haplodup_dir, "Hap1_vs_Hap2", coord_tables["Hap1_vs_Hap2"][0] , options.reuse_dotplots , options.skip_dotplots , paired_only_map = _pm_h1h2)

	if options.reference :
		print('[' + str(datetime.datetime.now()) + "] == Hap1 vs Reference", file=sys.stdout)
		outfile_prefix = "Hap1.on.Ref"
		plot_files["Hap1_vs_Reference"] = {}
		plot_files["Hap1_vs_Reference"]["Whole"] , plot_files["Hap1_vs_Reference"]["All_Dotplots"] = whole_genome_dotplot( reference_ids , hap1_ids , outfile_prefix, haplodup_dir, "Hap1_vs_Reference", coord_tables["Hap1_vs_Reference"][0] , options.reuse_dotplots , options.skip_dotplots , paired_only_map = _pm_refh1)

		print('[' + str(datetime.datetime.now()) + "] == Reference vs Hap1 ", file=sys.stdout)
		outfile_prefix = "Ref.on.Hap1"
		plot_files["Reference_vs_Hap1"] = {}
		plot_files["Reference_vs_Hap1"]["Whole"] , plot_files["Reference_vs_Hap1"]["All_Dotplots"] = whole_genome_dotplot( hap1_ids , reference_ids , outfile_prefix, haplodup_dir, "Reference_vs_Hap1", coord_tables["Reference_vs_Hap1"][0], options.reuse_dotplots , options.skip_dotplots , paired_only_map = _pm_h1ref)

		print('[' + str(datetime.datetime.now()) + "] == Hap2 vs Reference", file=sys.stdout)
		outfile_prefix = "Hap2.on.Ref"
		plot_files["Hap2_vs_Reference"] = {}
		plot_files["Hap2_vs_Reference"]["Whole"] , plot_files["Hap2_vs_Reference"]["All_Dotplots"] = whole_genome_dotplot( reference_ids , hap2_ids , outfile_prefix, haplodup_dir, "Hap2_vs_Reference", coord_tables["Hap2_vs_Reference"][0] , options.reuse_dotplots , options.skip_dotplots , paired_only_map = _pm_refh2)

		print('[' + str(datetime.datetime.now()) + "] == Reference vs Hap2", file=sys.stdout)
		outfile_prefix = "Ref.on.Hap2"
		plot_files["Reference_vs_Hap2"] = {}
		plot_files["Reference_vs_Hap2"]["Whole"] , plot_files["Reference_vs_Hap2"]["All_Dotplots"] = whole_genome_dotplot( hap2_ids , reference_ids , outfile_prefix, haplodup_dir, "Reference_vs_Hap2", coord_tables["Reference_vs_Hap2"][0] , options.reuse_dotplots , options.skip_dotplots , paired_only_map = _pm_h2ref)

	if options.gff :
		print('[' + str(datetime.datetime.now()) + "] = Generating reports", file=sys.stdout)

		report_tasks = []  # (comparison, queryID, fn) — fn() returns the result dict

		for comparison in list(coord_tables.keys()) :
			plot_files[comparison]["Reports"] = {}
			if isinstance(coord_tables[comparison], list):
				if len(coord_tables[comparison]) == 2 :
					coords_file , coords_file_self = coord_tables[comparison]
				else :
					print("[ERROR] QC comparison with unexpected data content", file=sys.stdout)
					print("[ERROR] QC comparison with unexpected data content: " + comparison + " >>> " + coord_tables[comparison], file=sys.stderr)
					sys.exit(1)
			else :
				if not coord_tables[comparison] == ""  :
					coords_file = coord_tables[comparison]
					coords_file_self = ""
				else :
					print("[ERROR] QC comparison with unexpected data content", file=sys.stdout)
					print("[ERROR] QC comparison " + comparison + " has no data content associated", file=sys.stderr)
					sys.exit(1)
			print('[' + str(datetime.datetime.now()) + "] == " + comparison, file=sys.stdout)
			outdir_name = haplodup_dir + "/" + comparison
			if comparison == "Hap1_vs_Reference" :
				for queryID in sorted(hap1_ids.split(",")) :
					refID = hap1_to_ref[queryID]
					hap1ID = queryID ; hap2ID = hap1_to_hap2[hap1ID]
					hap1Len = fasta_1_len[hap1ID] ; hap2Len = fasta_2_len[hap2ID]
					def _fn(qID=queryID, rID=refID, h1=hap1ID, h2=hap2ID, h1L=hap1Len, h2L=hap2Len, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_pair_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, "" , "" , "" , "" , h1 , h2 , h1L , h2L , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) ,
						         "pdf" : make_no_genes_pdf_report(cf , cfs, haplodup_dir, odn, qID, rID, hit_len , hit_iden) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Reference_vs_Hap1" :
				for queryID in sorted(ref_ids) :
					refID = ref_to_hap1[queryID]
					hap1ID = refID ; hap2ID = hap1_to_hap2[hap1ID]
					hap1Len = fasta_1_len[hap1ID] ; hap2Len = fasta_2_len[hap2ID]
					def _fn(qID=queryID, rID=refID, h1=hap1ID, h2=hap2ID, h1L=hap1Len, h2L=hap2Len, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_pair_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , h1L , h2L , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) ,
						         "pdf" : make_pair_pdf_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , "diploid_gene_count_trace.hap1.txt", hit_len , hit_iden, unbalanced_ratio) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Hap1_vs_Hap1" :
				for queryID in sorted(hap1_ids.split(",")) :
					refID = queryID
					hap1ID = queryID ; hap2ID = hap1_to_hap2[hap1ID]
					hap1Len = fasta_1_len[hap1ID] ; hap2Len = fasta_2_len[hap2ID]
					def _fn(qID=queryID, rID=refID, h1=hap1ID, h2=hap2ID, h1L=hap1Len, h2L=hap2Len, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_pair_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , h1L , h2L , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) ,
						         "pdf" : make_pair_pdf_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , "diploid_gene_count_trace.hap1.txt", hit_len , hit_iden, unbalanced_ratio) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Hap2_vs_Reference" :
				for queryID in sorted(hap2_ids.split(",")) :
					refID = hap2_to_ref[queryID]
					hap1ID = hap2_to_hap1[queryID] ; hap2ID = queryID
					hap1Len = fasta_1_len[hap1ID] ; hap2Len = fasta_2_len[hap2ID]
					def _fn(qID=queryID, rID=refID, h1=hap1ID, h2=hap2ID, h1L=hap1Len, h2L=hap2Len, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_pair_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, "" , "" , "" , "" , h1 , h2 , h1L , h2L , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) ,
						         "pdf" : make_no_genes_pdf_report(cf , cfs, haplodup_dir, odn, qID, rID, hit_len , hit_iden) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Reference_vs_Hap2" :
				for queryID in sorted(ref_ids) :
					refID = ref_to_hap2[queryID]
					hap1ID = hap2_to_hap1[refID] ; hap2ID = refID
					hap1Len = fasta_1_len[hap1ID] ; hap2Len = fasta_2_len[hap2ID]
					def _fn(qID=queryID, rID=refID, h1=hap1ID, h2=hap2ID, h1L=hap1Len, h2L=hap2Len, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_pair_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , h1L , h2L , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) ,
						         "pdf" : make_pair_pdf_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Hap2_vs_Hap1" :
				for hap2ID in sorted(hap2_ids.split(",")) :
					hap1ID = hap2_to_hap1[hap2ID] ; refID = hap1ID ; queryID = hap2ID
					hap1Len = fasta_1_len[hap1ID] ; hap2Len = fasta_2_len[hap2ID]
					def _fn(qID=queryID, rID=refID, h1=hap1ID, h2=hap2ID, h1L=hap1Len, h2L=hap2Len, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_pair_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , h1L , h2L , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) ,
						         "pdf" : make_pair_pdf_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , "diploid_gene_count_trace.hap1.txt", hit_len , hit_iden, unbalanced_ratio) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Hap1_vs_Hap2" :
				for hap1ID in sorted(hap1_ids.split(",")) :
					hap2ID = hap1_to_hap2[hap1ID] ; queryID = hap1ID ; refID = hap2ID
					hap1Len = fasta_1_len[hap1ID] ; hap2Len = fasta_2_len[hap2ID]
					def _fn(qID=queryID, rID=refID, h1=hap1ID, h2=hap2ID, h1L=hap1Len, h2L=hap2Len, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_pair_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , h1L , h2L , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) ,
						         "pdf" : make_pair_pdf_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Hap2_vs_Hap2" :
				for hap2ID in sorted(hap2_ids.split(",")) :
					hap1ID = hap2_to_hap1[hap2ID] ; refID = hap2ID ; queryID = hap2ID
					hap1Len = fasta_1_len[hap1ID] ; hap2Len = fasta_2_len[hap2ID]
					def _fn(qID=queryID, rID=refID, h1=hap1ID, h2=hap2ID, h1L=hap1Len, h2L=hap2Len, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_pair_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , h1L , h2L , "diploid_gene_count_trace.hap1.txt", "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) ,
						         "pdf" : make_pair_pdf_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , h1 , h2 , "diploid_gene_count_trace.hap2.txt", hit_len , hit_iden, unbalanced_ratio) }
					report_tasks.append((comparison, queryID, _fn))
			else :
				print("[ERROR] Report required for unknown comparison: " + comparison, file=sys.stderr)
				sys.exit(1)

		# Run all report tasks in parallel — subprocess.communicate() releases the GIL so threads
		# genuinely overlap while waiting for Rscript processes to finish.
		with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.cores)) as executor :
			futures = { executor.submit(fn) : (comp, qID) for comp, qID, fn in report_tasks }
			for future in concurrent.futures.as_completed(futures) :
				comp , qID = futures[future]
				plot_files[comp]["Reports"][qID] = future.result()

	else:
		print('[' + str(datetime.datetime.now()) + "] = Generating reports", file=sys.stdout)

		report_tasks = []

		for comparison in list(coord_tables.keys()) :
			plot_files[comparison]["Reports"] = {}
			if isinstance(coord_tables[comparison], list):
				if len(coord_tables[comparison]) == 2 :
					coords_file , coords_file_self = coord_tables[comparison]
				else :
					print("[ERROR] QC comparison with unexpected data content", file=sys.stdout)
					print("[ERROR] QC comparison with unexpected data content: " + comparison + " >>> " + coord_tables[comparison], file=sys.stderr)
					sys.exit(1)
			else :
				if not coord_tables[comparison] == ""  :
					coords_file = coord_tables[comparison]
					coords_file_self = ""
				else :
					print("[ERROR] QC comparison with unexpected data content", file=sys.stdout)
					print("[ERROR] QC comparison " + comparison + " has no data content associated", file=sys.stderr)
					sys.exit(1)
			print('[' + str(datetime.datetime.now()) + "] == " + comparison, file=sys.stdout)
			outdir_name = haplodup_dir + "/" + comparison
			mkdir(outdir_name)
			if comparison == "Hap1_vs_Reference" :
				for queryID in sorted(hap1_ids.split(",")) :
					refID = hap1_to_ref[queryID]
					def _fn(qID=queryID, rID=refID, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_no_genes_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, "" , "" , "" , "" , hit_len , hit_iden) ,
						         "pdf"  : make_no_genes_pdf_report( os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, "" , "" , "" , "" , hit_len , hit_iden) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Reference_vs_Hap1" :
				for queryID in sorted(ref_ids) :
					refID = ref_to_hap1[queryID]
					def _fn(qID=queryID, rID=refID, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_no_genes_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hit_len , hit_iden) ,
						         "pdf"  : make_no_genes_pdf_report( os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hit_len , hit_iden) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Hap1_vs_Hap1" :
				continue
			elif comparison == "Hap2_vs_Reference" :
				for queryID in sorted(hap2_ids.split(",")) :
					refID = hap2_to_ref[queryID]
					def _fn(qID=queryID, rID=refID, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_no_genes_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, "" , "" , "" , "" , hit_len , hit_iden) ,
						         "pdf"  : make_no_genes_pdf_report( os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, "" , "" , "" , "" , hit_len , hit_iden) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Reference_vs_Hap2" :
				for queryID in sorted(ref_ids) :
					refID = ref_to_hap2[queryID]
					def _fn(qID=queryID, rID=refID, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_no_genes_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hit_len , hit_iden) ,
						         "pdf"  : make_no_genes_pdf_report( os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hit_len , hit_iden) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Hap2_vs_Hap1" :
				for queryID in sorted(hap2_ids.split(",")) :
					refID = hap2_to_hap1[queryID]
					def _fn(qID=queryID, rID=refID, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_no_genes_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hit_len , hit_iden) ,
						         "pdf"  : make_no_genes_pdf_report( os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hit_len , hit_iden) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Hap1_vs_Hap2" :
				for queryID in sorted(hap1_ids.split(",")) :
					refID = hap1_to_hap2[queryID]
					def _fn(qID=queryID, rID=refID, cf=coords_file, cfs=coords_file_self, odn=outdir_name) :
						return { "html" : make_no_genes_html_report(os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hit_len , hit_iden) ,
						         "pdf"  : make_no_genes_pdf_report( os.path.basename(cf), os.path.basename(cfs), os.path.realpath(haplodup_dir), os.path.realpath(odn), qID, rID, structure_file , legacy_structure_file , all_markers_file_name , dup_markers_file_name , hit_len , hit_iden) }
					report_tasks.append((comparison, queryID, _fn))
			elif comparison == "Hap2_vs_Hap2" :
				continue
			else :
				print("[ERROR] Report required for unknown comparison: " + comparison, file=sys.stderr)
				sys.exit(1)

		# Run all report tasks in parallel — subprocess.communicate() releases the GIL so threads
		# genuinely overlap while waiting for Rscript processes to finish.
		with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.cores)) as executor :
			futures = { executor.submit(fn) : (comp, qID) for comp, qID, fn in report_tasks }
			for future in concurrent.futures.as_completed(futures) :
				comp , qID = futures[future]
				plot_files[comp]["Reports"][qID] = future.result()

	# Make Index
	html_index = make_index_from_report_db("index.html" , "." , haplodup_dir ,  plot_files  )


	# Windows with hotspots
	window_size = int(options.window)
	# On hit_counts_1 and hit_counts_2
	# hit_counts_1:
	# hit_counts_1[chr_hap1]= [
	#	...
	# 	[chr_hap1, start , end , id , h1_len , h2_len , ratio , descriptions , counts]
	#	...
	#	]
	if window_size > 0 :
		print('[' + str(datetime.datetime.now()) + '] === Searching for windows of genes with unbalanced count', file=sys.stdout)
		print('### Searching for windows of genes with unbalanced count', file=sys.stdout)
		max_unbalanced = int(options.allowed)
		ratio_threshold = float(unbalanced_ratio)
		inverted_ratio_threshold = 1/float(unbalanced_ratio)

		hotspot_windows = {}
		for chr in sorted(hit_counts_1.keys()) :
			hotspot_windows[chr] = []
			gene_counts = sorted(hit_counts_1[chr])
			# Sorted by chr id (always the same) and then by gene start then by gene stop
			gene_counts_len = len(gene_counts)
			# Make windows
			for i in range(0,gene_counts_len + 1 - window_size ) :
				start = gene_counts[i][1]
				stop = gene_counts[i+window_size-1][2]
				count_off = 0
				for j in range(i , i + window_size ) :
					if int(gene_counts[j][4]) == 0 or int(gene_counts[j][5]) == 0 :
						if (int(gene_counts[j][4]) == 0 and int(gene_counts[j][5]) > 1 ) or (int(gene_counts[j][4])>1 and int(gene_counts[j][5]) == 0) :
							# One of the two counts is null and the other is 2 or more
							count_off += 1
					else :
						# None of the count is null, check the ratio
						ratio = float(gene_counts[j][5]) / float(gene_counts[j][4])
						if ratio <= ratio_threshold or ratio >= inverted_ratio_threshold :
							count_off += 1
				if count_off > max_unbalanced :
					hotspot_windows[chr].append( [ chr , int(start) , int(stop )] )
		hotspot_hap1_file = haplodup_dir + "/" + options.out + ".hotspots.windows.hap1.txt"
		hotspot_hap1 = open(hotspot_hap1_file , "w+")
		for chr in sorted(hotspot_windows.keys()) :
			for element in sorted(hotspot_windows[chr]) :
				print("\t".join([ str(x) for x in element ]), file=hotspot_hap1)
		hotspot_hap1.close()

		hotspot_windows = {}
		for chr in sorted(hit_counts_2.keys()) :
			hotspot_windows[chr] = []
			gene_counts = sorted(hit_counts_2[chr])
			# Sorted by chr id (always the same) and then by gene start then by gene stop
			gene_counts_len = len(gene_counts)
			# Make windows
			for i in range(0,gene_counts_len + 1 - window_size ) :
				start = gene_counts[i][1]
				stop = gene_counts[i+window_size-1][2]
				count_off = 0
				for j in range(i , i + window_size ) :
					if int(gene_counts[j][4]) == 0 or int(gene_counts[j][5]) == 0 :
						if (int(gene_counts[j][4]) == 0 and int(gene_counts[j][5]) > 1 ) or (int(gene_counts[j][4])>1 and int(gene_counts[j][5]) == 0) :
							# One of the two counts is null and the other is 2 or more
							count_off += 1
					else :
						# None of the count is null, check the ratio
						ratio = float(gene_counts[j][5]) / float(gene_counts[j][4])
						if ratio <= ratio_threshold or ratio >= inverted_ratio_threshold :
							count_off += 1
				if count_off > max_unbalanced :
					hotspot_windows[chr].append( [ chr , int(start) , int(stop )] )
		hotspot_hap2_file = haplodup_dir + "/" + options.out + ".hotspots.windows.hap2.txt"
		hotspot_hap2 = open(hotspot_hap2_file , "w+")
		for chr in sorted(hotspot_windows.keys()) :
			for element in sorted(hotspot_windows[chr]) :
				print("\t".join([ str(x) for x in element ]), file=hotspot_hap2)
		hotspot_hap2.close()


	print("------------------------------", file=sys.stdout)
	print("- Done", file=sys.stdout)
	print("------------------------------", file=sys.stdout)
	print("##############################", file=sys.stderr)
	print("# Done", file=sys.stderr)
	print("##############################", file=sys.stderr)



if __name__ == '__main__':
	main()