#!/usr/bin/env python

import argparse
import shlex
from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *



def main() :

	###### Options and help ######
	parser = argparse.ArgumentParser()

	parser.add_argument("-1", "--hap1", dest="hap1",
					help="FASTA sequences of the first haplotype. [REQUIRED]", metavar="Hap1.fasta")
	parser.add_argument("-2", "--hap2", dest="hap2",
					help="FASTA sequences of the second haplotype. [REQUIRED]", metavar="Hap2.fasta")
	parser.add_argument("-U", "--unplaced", dest="unpl",
					help="FASTA sequences of the unplaced sequences.", metavar="Unplaced.fasta")
	parser.add_argument("-c", "--haplotype_correspondence", dest="corr",
					help="Tab separated values file of corresponding sequence names in the two haplotypes. [REQUIRED]", metavar="Hap1_to_Hap2.txt")


	parser.add_argument("--exclusion", dest="exclusion",
					help="Tab separated values file of unplaced sequences incompatible with a given input pseudomolecule", metavar="exclusion_list.tsv")
	parser.add_argument("--known", dest="known",
					help="Tab separated values file of unplaced sequences with known pseudomolecule association", metavar="known.tsv")

	parser.add_argument("-r", "--repeats",  dest="repeats",
					help="GFF/BED file of repetitive sequences [REQUIRED]", metavar="repeats.bed")
	parser.add_argument("--repeats_format",  dest="repeats_format", default="BED",
					help="Repeat annotation file format [default: BED]", metavar="[BED|GFF3|GFF|GTF]")
	parser.add_argument("-o", "--output", dest="output", default="out",
					help="Output files prefix [Default: out]", metavar="outprefix")
	parser.add_argument("--C12", "--Coverage12", dest="C12",
					help="File of per base sequencing coverage of both haplotypes (obtained with: 'bedtools genomecov -d -ibam') [REQUIRED]. To perform analysis, set '-C' and ['-s'|'--b1','--b2']", metavar="cov1.txt")
	parser.add_argument("--C1", "--c1" , "--Coverage1", dest="C1",
					help="File of the haplotype 1 per base sequencing coverage (obtained with: 'bedtools genomecov -d -ibam') [REQUIRED]. To perform analysis, set '-C' and ['-s'|'--b1','--b2']", metavar="cov1.txt")
	parser.add_argument("--C2", "--c2" , "--Coverage2", dest="C2",
					help="File of the haplotype 2 per base sequencing coverage (obtained with: 'bedtools genomecov -d -ibam') [REQUIRED]. To perform analysis, set '-C' and ['-s'|'--b1','--b2']", metavar="cov2.txt")
	parser.add_argument("-C", "--calculate_coverage" , dest="calculate_coverage", default=False, action="store_true",
					help="If set, run coverage calculation. [default: false] Requires -s")
	parser.add_argument("-s", "--sequenced_reads", dest="reads",
					help="FASTA/FASTQ (eventually gzipped) of sequenced reads. Required for '-C'", metavar="reads.fastq.gz")
	parser.add_argument("--map_threads", dest="map_threads", default="4",
					help="Number of threads to run mapping with. [default: 4]", metavar="N")
	parser.add_argument("--sequencing_technology", dest="tech", default="PacBio",
					help="Sequencing technology used for reads in '-s'. [default: PacBio]", metavar="[PacBio|ONT|Illumina_pe|Illumina_se]")
	parser.add_argument("--b1", dest="b1",
					help="BAM alignment of reads on the first haplotype. Required for '-C'", metavar="reads.on.hap1.bam")
	parser.add_argument("--b2", dest="b2",
					help="BAM alignment of reads on the second haplotype. Required for '-C'", metavar="reads.on.hap2.bam")
	parser.add_argument("-e" , "--expected_coverage", dest="expected_cov",
					help="Expected coverage threshold to be used for local ploidy evaluation. If not set, median value is calculated and used", metavar="100")
	parser.add_argument("-t" , "--temp" , dest="tmp", default="tmp_HaploFill" ,
					help="Path to temp folder. Required for resuming a process [default: ./tmp_HaploFill]", metavar="path/to/tmp")
	parser.add_argument("--resume" , dest="resume_step", default=0,
					help="Resume processing from given step", metavar="N")
	parser.add_argument("--stop" , dest="stop_step", default=11,
					help="Stop processing at given step", metavar="N")
	parser.add_argument("--overwrite" , dest="overwrite", default=False, action="store_true",
					help="If set, force to overwrite content in the temporary folder")
	parser.add_argument( "--flanking" , dest="flanking" , default=150000 ,
					help="Size of the flanking region around gaps [default: 150000]", metavar="N")
	parser.add_argument( "--coverage" , dest="coverage" , default=20 ,
					help="Minimum coverage percentage of the filler to match the supporting sequence(s) [default: 20]", metavar="0-100")
	parser.add_argument( "--nohomozygous" , dest="no_homozygous" , default=False, action="store_true",
					help="Do not search and output homozygous fillers")

	print("Running HaploFill tool from HaploSync version " + get_version(), file=sys.stdout)
	print("To reproduce this run use the following command: " + " ".join( shlex.quote(x) for x in sys.argv), file=sys.stdout)
	print("----", file=sys.stdout)
	scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/support_scripts"
	# Sanity Check

	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()

	if not ( options.hap1 and options.hap2 ) :
		print("[ERROR] FASTA file missing", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	if not options.corr :
		print("[ERROR] Sequence name correspondence file missing", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	if options.calculate_coverage and not ( options.reads or ( options.b1 and options.b2 ) ) :
		print("[ERROR] Coverage analysis requested but no reads file provided", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	paths = set_paths(os.path.join(sys.path[0], 'HaploSync.conf.toml'))
	minimap_path = paths["minimap2"]
	samtools_path = paths["samtools"]
	bedtools_path = paths["bedtools"]
	nucmer_path = paths["nucmer"]
	showcoords_path = paths["show-coords"]


	###### Main ######
	### Print intro

	temp_folder = options.tmp
	if not os.path.exists( options.tmp ) :
		mkdir(temp_folder)

	if options.overwrite :
		shutil.rmtree(temp_folder)
		mkdir(temp_folder)

	conf_file = temp_folder + "/conf.files.json"
	pairs_file = temp_folder + "/conf.pairs.json"
	status_file = temp_folder + "/status.json"

	if os.path.exists(conf_file) :
		files_and_folders = json.load( open( conf_file ) )
	else :
		files_and_folders = {}
	if os.path.exists(pairs_file) :
		pairs = json.load( open( pairs_file ) )
	else :
		pairs = {}

	if os.path.exists(status_file) :
		status = json.load( open(status_file) )
	else :
		status = {}
	if ("1-setup" not in status) or (status["1-setup"] == {}) :
		status["1-setup"] = {
							"1.1-split": "TODO" ,
							"1.2-pairs" : "TODO" ,
							"1.3-repeat": "TODO" ,
							"1.4-gap": "TODO"
							}
	if ("2-coverage" not in status) or (status["2-coverage"] == {}) :
		status["2-coverage"] = {
							"2.1-map1" : "TODO" ,
							"2.2-map2" : "TODO" ,
							"2.3-cov1" : "TODO" ,
							"2.4-cov2" : "TODO" ,
							"2.5-split" : "TODO"
							}
	if ("3-ploidy" not in status) or (status["3-ploidy"] == {}) :
		status["3-ploidy"] = {
							"3.1-median" : "TODO" ,
							"3.2-categorize" : "TODO"
							}
	if ("4-hap2hap" not in status) or (status["4-hap2hap"] == {}) :
		status["4-hap2hap"] = {
							"4.1-map" : "TODO" ,
							"4.2-uniquify" : "TODO",
							"4.3-pairing" : "TODO"
							}
	if ("5-upgradeable" not in status) or (status["5-upgradeable"] == {}) :
		status["5-upgradeable"] = {
							"5.1-Unmatched": "TODO" ,
							"5.2-preprocess_gaps" : "TODO" ,
							"5.3-get_sequences" : "TODO"
							}
	if ("6-filling" not in status) or (status["6-filling"] == {}) :
		status["6-filling"] = {
							"6.1-gather": "TODO" ,
							"6.2-map" : "TODO" ,
							"6.3-select" : "TODO"
							}


	unwanted_pairs = {}
	if options.exclusion :
		print('[' + str(datetime.datetime.now()) + "] = Reading mutually exclusive sequences pairs", file=sys.stdout)
		for line in open(options.exclusion) :
			if ( line.rstrip() == "" ) or ( line[0] == "#") :
				continue
			else:
				try :
					seq_id , chr_id = line.rstrip().split("\t")
				except :
					print("[ERROR] Line in exclude pair file (" + options.exclusion + ") does not contain two ids", file=sys.stdout)
					print("[ERROR] Line in exclude pair file (" + options.exclusion + ") does not contain two ids", file=sys.stderr)
					print("[ERROR] Line: ", file=sys.stderr)
					print("[ERROR] " + line.rstrip(), file=sys.stderr)
					sys.exit(1)
				if (seq_id + "|+") not in unwanted_pairs :
					unwanted_pairs[seq_id] = []
					unwanted_pairs[seq_id + "|+"] = []
					unwanted_pairs[seq_id + "|-"] = []
					unwanted_pairs[seq_id + "|."] = []
				unwanted_pairs[seq_id].append(chr_id)
				unwanted_pairs[seq_id + "|+"].append(chr_id)
				unwanted_pairs[seq_id + "|-"].append(chr_id)
				unwanted_pairs[seq_id + "|."].append(chr_id)

	# Known relationships
	known_chr_by_seqid = {}
	if options.known :
		print('[' + str(datetime.datetime.now()) + "] = Reading known groups of sequences in the same haplotype", file=sys.stdout)
		for line in open(options.known) :
			if ( line == "" ) or ( line[0] == "#") :
				continue
			else:
				try :
					seq_id , chr_id = line.rstrip().split("\t")
				except :
					print("[ERROR] Line in exclude pair file (" + options.known + ") does not contain two ids", file=sys.stdout)
					print("[ERROR] Line in exclude pair file (" + options.known + ") does not contain two ids", file=sys.stderr)
					print("[ERROR] Line: ", file=sys.stderr)
					print("[ERROR] " + line.rstrip(), file=sys.stderr)
					sys.exit(1)
				if seq_id not in known_chr_by_seqid :
					known_chr_by_seqid[seq_id] = []
					known_chr_by_seqid[seq_id + "|+" ] = []
					known_chr_by_seqid[seq_id + "|-" ] = []
					known_chr_by_seqid[seq_id + "|." ] = []
				known_chr_by_seqid[seq_id].append(chr_id)
				known_chr_by_seqid[seq_id + "|+" ].append(chr_id)
				known_chr_by_seqid[seq_id + "|-" ].append(chr_id)
				known_chr_by_seqid[seq_id + "|." ].append(chr_id)



	######  STEP 1: SETUP ######
	if int(options.resume_step) < 2 :
		###### Read files and split in folders ######
		print('[' + str(datetime.datetime.now()) + '] = STEP 1: Setup', file=sys.stdout)
		print('# STEP 1: Setup', file=sys.stderr)
		if "inputs" in files_and_folders :
			if not files_and_folders["inputs"]["1"] == options.hap1 :
				print('## [WARNING] Hap1 FASTA file differs from the original one. Original sequences will be used', file=sys.stderr)
			if not files_and_folders["inputs"]["2"] == options.hap2 :
				print('## [WARNING] Hap2 FASTA file differs from the original one. Original sequences will be used', file=sys.stderr)
			if options.unpl and not files_and_folders["inputs"]["U"] == options.unpl :
				print('## [WARNING] Unplaced sequences FASTA file differs from the original one. Original sequences will be used', file=sys.stderr)
		else :
			files_and_folders["inputs"] = {}
			files_and_folders["inputs"]["1"] = options.hap1
			files_and_folders["inputs"]["2"] = options.hap2
			files_and_folders["inputs"]["1_list"] = get_fasta_ids(options.hap1)
			files_and_folders["inputs"]["2_list"] = get_fasta_ids(options.hap2)
			if options.unpl :
				files_and_folders["inputs"]["U"] = options.unpl
				files_and_folders["inputs"]["U_list"] = get_fasta_ids(options.unpl)
		if "tool_paths" in files_and_folders :
			if not files_and_folders["tool_paths"] == paths :
				print('## [MEMO] Executable paths updated', file=sys.stderr)
		else :
			files_and_folders["tool_paths"] = paths

		##### Split fasta files
		if status["1-setup"]["1.1-split"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] == STEP 1.1 - Splitting sequences', file=sys.stdout)
			print('## STEP 1.1 - Splitting sequences', file=sys.stderr)
			files_and_folders["sequences"] = {}
			####  Hap1
			print('[' + str(datetime.datetime.now()) + '] === STEP 1.1.1 - Splitting Hap1 FASTA', file=sys.stdout)
			print('### STEP 1.1.1 - Splitting Hap1 FASTA', file=sys.stderr)
			files_and_folders = organize_fasta( options.hap1, files_and_folders , temp_folder , "1" )
			#print >> sys.stderr, files_and_folders
			####  Hap2
			print('[' + str(datetime.datetime.now()) + '] === STEP 1.1.2 - Splitting Hap2 FASTA', file=sys.stdout)
			print('### STEP 1.1.2 - Splitting Hap2 FASTA', file=sys.stderr)
			files_and_folders = organize_fasta( options.hap2, files_and_folders , temp_folder , "2" )
			#print >> sys.stderr, files_and_folders
			#### Unplaced
			if options.unpl :
				print('[' + str(datetime.datetime.now()) + '] === STEP 1.1.3 - Splitting Unplaced sequences FASTA', file=sys.stdout)
				print('### STEP 1.1.3 - Splitting Unplaced sequences FASTA', file=sys.stderr)
				files_and_folders = organize_fasta( options.unpl , files_and_folders , temp_folder , "U")
				#print >> sys.stderr, files_and_folders
			status["1-setup"]["1.1-split"] = "DONE"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
		else :
			print('[' + str(datetime.datetime.now()) + '] == STEP 1.1 - Splitting completed, skipping', file=sys.stdout)
			print('## STEP 1.1 - Splitting completed, skipping', file=sys.stderr)

		##### Set pairing information for pseudomolecules
		if status["1-setup"]["1.2-pairs"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] == STEP 1.2 - Sequence pairing information setup', file=sys.stdout)
			print('## STEP 1.2 - Sequence pairing information setup', file=sys.stderr)

			pair_id = -1

			for line in open(options.corr) :
				pair_id += 1
				pairs[pair_id] = line.rstrip().split("\t")
				hap1_id , hap2_id = pairs[pair_id]
				pair_folder = temp_folder + "/Pair_" + str(pair_id)
				mkdir(pair_folder)

				files_and_folders["sequences"][hap1_id]["mate_id"] = hap2_id
				files_and_folders["sequences"][hap2_id]["mate_id"] = hap1_id
				files_and_folders["sequences"][hap1_id]["pair_id"] = pair_id
				files_and_folders["sequences"][hap2_id]["pair_id"] = pair_id
				files_and_folders["sequences"][hap1_id]["pair_folder"] = pair_folder
				files_and_folders["sequences"][hap2_id]["pair_folder"] = pair_folder
			status["1-setup"]["1.2-pairs"] = "DONE"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
		else :
			print('[' + str(datetime.datetime.now()) + '] == STEP 1.2 - Sequence already paired, skipping', file=sys.stdout)
			print('## STEP 1.2 - Sequence already paired, skipping', file=sys.stderr)

		##### Split Repeats
		if status["1-setup"]["1.3-repeat"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] == STEP 1.3 - Setup repeats information', file=sys.stdout)
			print('## STEP 1.3 - Setup repeats information', file=sys.stderr)
			if not options.repeats_format == "BED" :
				print('### Converting repeats file', file=sys.stderr)
				print('[' + str(datetime.datetime.now()) + '] === Converting repeats file', file=sys.stdout)
				repeats_file = gfx2bed_print( options.repeats , temp_folder + "/repeats-tmp.bed" )
				print('### Splitting Repeats', file=sys.stderr)
				files_and_folders = organize_bed( repeats_file , files_and_folders , "repeat" , False)
			else :
				repeats_file = options.repeats
				print('### Splitting Repeats', file=sys.stderr)
				files_and_folders = organize_bed( repeats_file , files_and_folders , "repeat" , False )
			status["1-setup"]["1.3-repeat"]="DONE"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
		else :
			print('[' + str(datetime.datetime.now()) + '] == STEP 1.3 - Repeat information already setup, skipping', file=sys.stdout)
			print('## STEP 1.3 - Repeat information already setup, skipping', file=sys.stderr)

		##### Get gaps
		if status["1-setup"]["1.4-gap"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] == STEP 1.4 - Detecting gaps in sequences', file=sys.stdout)
			print('## STEP 1.4 - Detecting gaps in sequences', file=sys.stderr)
			for sequence_id in sorted(files_and_folders["sequences"].keys()) :
				files_and_folders["sequences"][sequence_id]["gap_file"] = get_gap_bed(files_and_folders["sequences"][sequence_id]["fasta_file"])
			status["1-setup"]["1.4-gap"]="DONE"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
		else :
			print('[' + str(datetime.datetime.now()) + '] == STEP 1.4 - Gap information already setup, skipping', file=sys.stdout)
			print('## STEP 1.4 - Gap information already setup, skipping', file=sys.stderr)

	else :
		print('# STEP 1: Setup - Skip', file=sys.stderr)
		## Sanity check
		if status["1-setup"]["1.1-split"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] [ERROR] STEP 1 Setup cannot be skipped: Step 1.1 not completed', file=sys.stdout)
			print('[ERROR] STEP 1 Setup cannot be skipped: Step 1.1 (sequence splitting) not completed', file=sys.stderr)
			sys.exit(110)
		if status["1-setup"]["1.2-pairs"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] [ERROR] STEP 1 Setup cannot be skipped: Step 1.2 not completed', file=sys.stdout)
			print('[ERROR] STEP 1 Setup cannot be skipped: Step 1.2 (haplotype sequence paring) not completed', file=sys.stderr)
			sys.exit(120)
		if status["1-setup"]["1.3-repeat"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] [ERROR] STEP 1 Setup cannot be skipped: Step 1.3 not completed', file=sys.stdout)
			print('[ERROR] STEP 1 Setup cannot be skipped: Step 1.3 (repeats setup) not completed', file=sys.stderr)
			sys.exit(130)
		if status["1-setup"]["1.4-gap"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] [ERROR] STEP 1 Setup cannot be skipped: Step 1.4 not completed', file=sys.stdout)
			print('[ERROR] STEP 1 Setup cannot be skipped: Step 1.4 (gap search) not completed', file=sys.stderr)
			sys.exit(140)


	if int(options.stop_step) == 1 :
		print("------------------------------", file=sys.stdout)
		print("- Done", file=sys.stdout)
		print("------------------------------", file=sys.stdout)
		sys.exit(0)



	###### STEP 2: COVERAGE ######
	if int(options.resume_step) < 3 :
		#####  Run coverage or read coverage files
		print('# STEP 2: Coverage extraction', file=sys.stderr)
		print('[' + str(datetime.datetime.now()) + '] = STEP 2: Coverage extraction', file=sys.stdout)
		hap1_cov = temp_folder + "/cov1.txt.gz"
		hap2_cov = temp_folder + "/cov2.txt.gz"

		if not options.calculate_coverage :
			status["2-coverage"]["2.1-map1"] = "DONE"
			status["2-coverage"]["2.2-map2"] = "DONE"
			print('[' + str(datetime.datetime.now()) + '] == STEP 2.1: Importing coverage files', file=sys.stdout)
			print('## Importing coverage files', file=sys.stderr)
			if not status["2-coverage"]["2.5-split"] == "DONE" :
				if not ( status["2-coverage"]["2.3-cov1"] == "DONE" and status["2-coverage"]["2.4-cov2"] == "DONE") and options.C12 :
					hap_cov_bed = temp_folder + "/cov.bed.gz"
					if not options.C1[-3:] == ".gz" :
						compress_file( options.C1 , hap_cov_bed )
						print('### ' + options.C1 + 'file uncompressed, copying and compressing', file=sys.stderr)
					else :
						copy_file(options.C1 , hap_cov_bed)
					status["2-coverage"]["2.3-cov1"] = "DONE"
					status["2-coverage"]["2.4-cov2"] = "DONE"
					save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				elif not status["2-coverage"]["2.3-cov1"] == "DONE" :
					hap1_cov_bed = temp_folder + "/cov1.bed.gz"
					if not options.C1[-3:] == ".gz" :
						compress_file( options.C1 , hap1_cov_bed )
						print('### ' + options.C1 + 'file uncompressed, copying and compressing', file=sys.stderr)
					else :
						copy_file(options.C1 , hap1_cov_bed)
					status["2-coverage"]["2.3-cov1"] = "DONE"
					save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				elif not status["2-coverage"]["2.4-cov2"] == "DONE" :
					hap2_cov_bed = temp_folder + "/cov2.bed.gz"
					if not options.C2[-3:] == ".gz" :
						compress_file( options.C2 , hap2_cov_bed )
						print('### ' + options.C2 + 'file uncompressed, copying and compressing', file=sys.stderr)
					else :
						copy_file( options.C2 , hap2_cov_bed )
					status["2-coverage"]["2.4-cov2"] = "DONE"
					save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				else :
					status["2-coverage"]["2.3-cov1"] = "FAILED"
					status["2-coverage"]["2.4-cov2"] = "FAILED"
					save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
					print("[ERROR] Input files and missing analysis incompatibility issue. Please run again the process", file=sys.stdout)
					print('[ERROR] Coverage BED files given as input do not allow to recover some of the missing data', file=sys.stderr)
					print('[ERROR] If no coverage data was given as input, please rerun feeding it with --C12 option or --C1 + --C2 oprions', file=sys.stderr)
					print('[ERROR] If --C12 option was present, please just rerun the same command line, status has been updated to allow recover', file=sys.stderr)
					sys.exit(20)

				print('[' + str(datetime.datetime.now()) + '] == STEP 2.2: Splitting by sequence', file=sys.stdout)
				print('## Splitting coverage traces', file=sys.stderr)
				#### Split given files
				if os.path.exists(temp_folder + "/cov.bed.gz"):
					### Hap1+2
					print('[' + str(datetime.datetime.now()) + '] === STEP 2.2.1: Both haplotypes ', file=sys.stdout)
					print('### Both haplotypes', file=sys.stderr)
					files_and_folders = split_coverage_bed( temp_folder + "/cov.bed.gz" , files_and_folders )
				elif os.path.exists(temp_folder + "/cov1.bed.gz") and os.path.exists(temp_folder + "/cov2.bed.gz") :
					### Hap1
					print('[' + str(datetime.datetime.now()) + '] === STEP 2.2.1: Hap1 ', file=sys.stdout)
					print('### Hap1', file=sys.stderr)
					files_and_folders = split_coverage_bed( temp_folder + "/cov1.bed.gz" , files_and_folders )
					### Hap2
					print('[' + str(datetime.datetime.now()) + '] === STEP 2.2.1: Hap2 ', file=sys.stdout)
					print('### Hap2', file=sys.stderr)
					files_and_folders = split_coverage_bed( temp_folder + "/cov2.bed.gz" , files_and_folders )
					status["2-coverage"]["2.5-split"] = "DONE"
				else :
					print('[' + str(datetime.datetime.now()) + '] === STEP 2.2.1: Checking each sequence ', file=sys.stdout)
					print('### Checking each sequence', file=sys.stderr)
					to_check = files_and_folders["inputs"]["1_list"] + files_and_folders["inputs"]["2_list"]
					OK = True
					for sequence_to_check in to_check :
						if "coverage_file" in files_and_folders["sequences"][sequence_to_check] :
							if os.path.exists(files_and_folders["sequences"][sequence_to_check]["coverage_file"]) :
								print('#### ' + sequence_to_check + " OK ", file=sys.stderr)
							else :
								print('#### [ERROR] ' + sequence_to_check + " missing coverage ", file=sys.stderr)
								OK = False
						else :
							OK = False
					if not OK :
						print('[ERROR] Some sequences are missing coverage information. See standard error for more insight.', file=sys.stdout)
						print('[ERROR] Some sequences failed the coverage information check.', file=sys.stderr)
						status["2-coverage"]["2.5-split"] = "FAILED"
						save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
						sys.exit(25)
					else :
						print('#### All sequences are OK', file=sys.stderr)
			else :
				print('[' + str(datetime.datetime.now()) + '] === STEP 2.1.1: Checking each sequence for coverage information', file=sys.stdout)
				print('### Checking each sequence', file=sys.stderr)
				to_check = files_and_folders["inputs"]["1_list"] + files_and_folders["inputs"]["2_list"]
				OK = True
				for sequence_to_check in to_check :
					if "coverage_file" in files_and_folders["sequences"][sequence_to_check] :
						if os.path.exists(files_and_folders["sequences"][sequence_to_check]["coverage_file"]) :
							print('#### ' + sequence_to_check + " OK ", file=sys.stderr)
						else :
							print('#### [ERROR] ' + sequence_to_check + " missing coverage ", file=sys.stderr)
							OK = False
					else :
						OK = False
				if not OK :
					print('[ERROR] Some sequences are missing coverage information. See standard error for more insight.', file=sys.stdout)
					print('[ERROR] Some sequences failed the coverage information check.', file=sys.stderr)
					status["2-coverage"]["2.5-split"] = "FAILED"
					save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
					sys.exit(25)
				else :
					print('#### All sequences are OK', file=sys.stderr)
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
		else :
			if not status["2-coverage"]["2.5-split"] == "DONE" :
				if not options.reads :
					if not status["2-coverage"]["2.1-map1"] == "DONE" :
						print('## Alignment file for Hap1 present', file=sys.stderr)
						print('### Importing Hap1 alignment', file=sys.stderr)
						original_bamfile1 = os.path.abspath(options.b1)
						bamfile1 = hap1_cov.rstrip(".txt.gz") + ".bam"
						os.symlink( original_bamfile1 , bamfile1 )
						print('#### Indexing Hap1 alignment', file=sys.stderr)
						index_bam( bamfile1 , samtools_path )
						status["2-coverage"]["2.1-map1"] = "DONE"
					if not status["2-coverage"]["2.2-map2"] == "DONE" :
						print('## Alignment file for Hap2 present', file=sys.stderr)
						print('### Importing Hap2 alignment', file=sys.stderr)
						original_bamfile2 = os.path.abspath(options.b2)
						bamfile2 = hap2_cov.rstrip(".txt.gz") + ".bam"
						os.symlink( original_bamfile2 , bamfile2 )
						print('#### Indexing Hap2 alignment', file=sys.stderr)
						index_bam( bamfile2 , samtools_path )
						status["2-coverage"]["2.2-map2"] = "DONE"
					save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				else :
					print('[' + str(datetime.datetime.now()) + '] == STEP 2.1: Mapping reads', file=sys.stdout)
					print('## Mapping reads', file=sys.stderr)
					parameters = " -a --MD -L "
					# Possible technologies [PacBio|ONT|Illumina_pe|Illumina_se]
					if options.tech == "PacBio" : parameters += " -x map-pb "
					elif options.tech == "ONT" : parameters += " -x map-ont "
					elif options.tech == "Illumina_se" : parameters += " -x sr "
					elif options.tech == "Illumina_pe" : parameters += " -x sr "
					else :
						print("[ERROR] Unrecognized ", file=sys.stderr)

					### Calculate coverage, save files
					## Hap 1
					print('[' + str(datetime.datetime.now()) + '] === STEP 2.1.1: Hap1', file=sys.stdout)
					print('### Hap1', file=sys.stderr)
					if status["2-coverage"]["2.1-map1"] == "DONE" :
						print('[' + str(datetime.datetime.now()) + '] ==== Alignment previously performed, skipping ', file=sys.stdout)
						bamfile1 = hap1_cov.rstrip(".txt.gz") + ".bam"
					else :
						samfile1 = hap1_cov.rstrip(".txt.gz") + ".sam"
						bamfile1 = hap1_cov.rstrip(".txt.gz") + ".bam"
						map_minimap( options.hap1 , options.reads.split(",") , int(options.map_threads) , parameters , samfile1 , minimap_path )
						sam2sorted_bam(samfile1 , bamfile1 , options.hap1 , int(options.map_threads) , samtools_path )
						status["2-coverage"]["2.1-map1"] = "DONE"
						save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
					## Hap 2
					print('[' + str(datetime.datetime.now()) + '] === STEP 2.1.2: Hap2', file=sys.stdout)
					print('### Hap2', file=sys.stderr)
					if status["2-coverage"]["2.2-map2"] == "DONE" :
						print('[' + str(datetime.datetime.now()) + '] ==== Alignment previously performed, skipping ', file=sys.stdout)
						bamfile2 = hap2_cov.rstrip(".txt.gz") + ".bam"
					else :
						samfile2 = hap2_cov.rstrip(".txt.gz") + ".sam"
						bamfile2 = hap2_cov.rstrip(".txt.gz") + ".bam"
						map_minimap( options.hap2 , options.reads.split(",") , int(options.map_threads) , parameters , samfile2 , minimap_path )
						sam2sorted_bam(samfile2 , bamfile2 , options.hap2 , int(options.map_threads) , samtools_path )
						status["2-coverage"]["2.2-map2"] = "DONE"
					save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

				# Extract and split
				print('[' + str(datetime.datetime.now()) + '] == STEP 2.2: Extracting per base coverage', file=sys.stdout)
				print('## Extracting per base coverage', file=sys.stderr)
				### Hap1
				if not status["2-coverage"]["2.3-cov1"] == "DONE" :
					print('[' + str(datetime.datetime.now()) + '] === STEP 2.2.1: Extracting Hap1 coverage', file=sys.stdout)
					print('### Hap1', file=sys.stderr)
					files_and_folders = write_coverage_bed(bamfile1 , files_and_folders , files_and_folders["inputs"]["1_list"] ,  bedtools_path , samtools_path)
					if options.C1 and not os.path.exists(options.C1) :
						copy_file( hap1_cov , os.path.dirname(options.C1) )
					status["2-coverage"]["2.3-cov1"] = "DONE"
				else :
					print('[' + str(datetime.datetime.now()) + '] === Extraction of Hap1 coverage previously perfored, skipping', file=sys.stdout)
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				### Hap2
				if not status["2-coverage"]["2.4-cov2"] == "DONE" :
					print('[' + str(datetime.datetime.now()) + '] === STEP 2.2.2: Extracting Hap2 coverage', file=sys.stdout)
					print('### Hap2', file=sys.stderr)
					files_and_folders = write_coverage_bed(bamfile2 , files_and_folders , files_and_folders["inputs"]["2_list"] , bedtools_path , samtools_path)
					if options.C2 and not os.path.exists(options.C2) :
						copy_file( hap2_cov , os.path.dirname(options.C2) )
					status["2-coverage"]["2.4-cov2"] = "DONE"
				else :
					print('[' + str(datetime.datetime.now()) + '] === Extraction of Hap2 coverage previously perfored, skipping', file=sys.stdout)
				status["2-coverage"]["2.5-split"] = "DONE"
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
			else :
				print('[' + str(datetime.datetime.now()) + '] == STEP 2.1: Checking each sequence for coverage information', file=sys.stdout)
				print('### Checking each sequence', file=sys.stderr)
				to_check = files_and_folders["inputs"]["1_list"] + files_and_folders["inputs"]["2_list"]
				OK = True
				for sequence_to_check in to_check :
					if "coverage_file" in files_and_folders["sequences"][sequence_to_check] :
						if os.path.exists(files_and_folders["sequences"][sequence_to_check]["coverage_file"]) :
							print('#### ' + sequence_to_check + " OK ", file=sys.stderr)
						else :
							print('#### [ERROR] ' + sequence_to_check + " missing coverage ", file=sys.stderr)
							OK = False
					else :
						OK = False
				if not OK :
					print('[ERROR] Some sequences are missing coverage information. See standard error for more insight.', file=sys.stdout)
					print('[ERROR] Some sequences failed the coverage information check.', file=sys.stderr)
					status["2-coverage"]["2.5-split"] = "FAILED"
					save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
					sys.exit(25)
				else :
					print('#### All sequences are OK', file=sys.stderr)
		save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

	else:
		print('[' + str(datetime.datetime.now()) + '] = STEP 2: Coverage extraction already performed, skipping', file=sys.stdout)
		print('# STEP 2: Coverage extraction - Skip', file=sys.stderr)
		print('[' + str(datetime.datetime.now()) + '] === STEP 2.1.1: Checking each sequence for coverage information', file=sys.stdout)
		print('### Checking each sequence', file=sys.stderr)
		to_check = files_and_folders["inputs"]["1_list"] + files_and_folders["inputs"]["2_list"]
		OK = True
		for sequence_to_check in to_check :
			if "coverage_file" in files_and_folders["sequences"][sequence_to_check] :
				if os.path.exists(files_and_folders["sequences"][sequence_to_check]["coverage_file"]) :
					print('#### ' + sequence_to_check + " OK ", file=sys.stderr)
				else :
					print('#### [ERROR] ' + sequence_to_check + " missing coverage ", file=sys.stderr)
					OK = False
			else :
				OK = False
		if not OK :
			print('[ERROR] Some sequences are missing coverage information. See standard error for more insight.', file=sys.stdout)
			print('[ERROR] Some sequences failed the coverage information check.', file=sys.stderr)
			status["2-coverage"]["2.5-split"] = "FAILED"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
			sys.exit(25)
		else :
			print('#### All sequences are OK', file=sys.stderr)
		save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)


	if int(options.stop_step) == 2 :
		print("------------------------------", file=sys.stdout)
		print("- Done", file=sys.stdout)
		print("------------------------------", file=sys.stdout)
		sys.exit(0)



	######  STEP 3: Local ploidy level classification ######
	if int(options.resume_step) < 4 :
		#status["3-ploidy"] = {
		#			"3.1-median" : "TODO" ,
		#			"3.2-categorize" : "TODO"
		#			}

		#####  Run coverage or read coverage files
		print('# STEP 3: Local ploidy level classification', file=sys.stderr)
		print('[' + str(datetime.datetime.now()) + '] = STEP 3: Local ploidy level classification', file=sys.stdout)

		smoothed_masked = {}
		masked_coverage = {}
		if ( status["3-ploidy"]["3.1-median"] == "TODO" ) or ( "coverage" not in files_and_folders ) :
			files_and_folders["coverage"] = {}
			status["3-ploidy"]["3.2-categorize"] = "TODO"
			print('## STEP 3.1: Calculating median coverage excluding repeats and gaps', file=sys.stderr)
			print('[' + str(datetime.datetime.now()) + '] == STEP 3.1: Calculating median coverage excluding repeats and gaps', file=sys.stdout)
			medianCoverage , masked_coverage , smoothed_masked = calculate_clean_median( files_and_folders )
			files_and_folders = write_masked_signal(files_and_folders, masked_coverage, smoothed_masked)
			files_and_folders = write_masked_coverage_bed(masked_coverage, smoothed_masked, files_and_folders)
			print('### Median coverage: ' + str(medianCoverage), file=sys.stderr)
			files_and_folders["coverage"]["calculated"] = medianCoverage
			print('[' + str(datetime.datetime.now()) + '] === Coverage calculated median coverage: ' + str(medianCoverage), file=sys.stdout)
			if options.expected_cov :
				print('### Expected coverage as given by user inputted: ' + str(options.expected_cov), file=sys.stderr)
				print('[' + str(datetime.datetime.now()) + '] === Expected coverage as given from user: ' + str(options.expected_cov), file=sys.stdout)
				if "given" in files_and_folders["coverage"] :
					# Check if the same stored
					if not options.expected_cov == files_and_folders["coverage"]["given"] :
						print('#### [WARNING] Expected coverage given differs from a previous run. It used to be ' + files_and_folders["coverage"]["given"], file=sys.stderr)
						print('#### [WARNING] Actual one will be used and stored', file=sys.stderr)
				files_and_folders["coverage"]["given"] = options.expected_cov
				try :
					plot_histogram( masked_coverage , "Coverage distribution" , "Coverage Depth (X-fold change)" , temp_folder + "/Coverage_distribution.hist.pdf" , medianCoverage , float(options.expected_cov) )
				except :
					plot_histogram( masked_coverage , "Coverage distribution" , "Coverage Depth (X-fold change)" , temp_folder + "/Coverage_distribution.hist.pdf" , medianCoverage )
			else :
				if "given" in files_and_folders["coverage"] :
					print('#### [MEMO] A user given expected coverage of ' + str( files_and_folders["coverage"]["given"] ) + " was stored in memory from a previous run", file=sys.stderr)
					print('#### [MEMO] Calculated one will be used and memory cleaned', file=sys.stderr)
					del(files_and_folders["coverage"]["given"])
				plot_histogram( masked_coverage , "Coverage distribution" , "Coverage Depth (X-fold change)" , temp_folder + "/Coverage_distribution.hist.pdf" , medianCoverage )
			status["3-ploidy"]["3.1-median"] = "DONE"
			print('## Saving status', file=sys.stderr)
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

		if status["3-ploidy"]["3.2-categorize"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] == STEP 3.2: Categorize status from coverage', file=sys.stdout)
			print('## STEP 3.2: Categorize status from coverage', file=sys.stderr)
			if "given" in files_and_folders["coverage"] :
				coverage_threshold = files_and_folders["coverage"]["given"]
			elif "calculated" in files_and_folders["coverage"] :
				coverage_threshold = files_and_folders["coverage"]["calculated"]
			else :
				print('[' + str(datetime.datetime.now()) + '] [ERROR] Something went wrong, no coverage threshold available', file=sys.stdout)
				print('[ERROR] Something went wrong, no coverage threshold available', file=sys.stderr)
				status["3-ploidy"]["3.1-median"] = "TODO"
				status["3-ploidy"]["3.2-categorize"] = "TODO"
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				exit(31)

			if not check_all_masked(files_and_folders) :
				print('[' + str(datetime.datetime.now()) + '] [ERROR] Something went wrong, coverage files are missing', file=sys.stdout)
				print('[ERROR] Something went wrong, coverage files are missing', file=sys.stderr)
				status["3-ploidy"]["3.1-median"] = "TODO"
				status["3-ploidy"]["3.2-categorize"] = "TODO"
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				exit(32)
			else :
				print('[' + str(datetime.datetime.now()) + '] === STEP 3.2.1: Loading existing coverage files', file=sys.stdout)
				print('### STEP 3.2.1: Loading existing coverage files', file=sys.stderr)
				if smoothed_masked == {} :
					masked_coverage , smoothed_masked = load_smoothed_masked(files_and_folders)
				files_and_folders = get_category( files_and_folders , smoothed_masked , coverage_threshold )
				status["3-ploidy"]["3.2-categorize"] = "DONE"
				##### Save configuration file
				print('## Saving status', file=sys.stderr)
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

	else :
		print('# STEP 3: Local ploidy level classification - Skip', file=sys.stderr)
		if not check_all_masked(files_and_folders) :
			print('[' + str(datetime.datetime.now()) + '] [ERROR] Something went wrong, coverage files are missing', file=sys.stdout)
			print('[ERROR] Something went wrong, coverage files are missing', file=sys.stderr)
			status["3-ploidy"]["3.1-median"] = "TODO"
			status["3-ploidy"]["3.2-categorize"] = "TODO"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
			exit(30)
		else :
			if not check_all_categories(files_and_folders) :
				print('[' + str(datetime.datetime.now()) + '] [ERROR] Something went wrong, category files are missing', file=sys.stdout)
				print('[ERROR] Something went wrong, category files are missing', file=sys.stderr)
				status["3-ploidy"]["3.2-categorize"] = "TODO"
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				exit(30)
		if "given" in files_and_folders["coverage"] :
			# Check if the same stored
			if not float(options.expected_cov) == files_and_folders["coverage"]["given"] :
				print('#### [WARNING] Expected coverage given differs from a previous run. It used to be ' + str( files_and_folders["coverage"]["given"] ) + '. Given threshold will be ignored', file=sys.stderr)
				print('#### [WARNING] If you want to use the given threshold instead of the stored one, please rerun from step 3', file=sys.stderr)


	if int(options.stop_step) == 3 :
		print("------------------------------", file=sys.stdout)
		print("- Done", file=sys.stdout)
		print("------------------------------", file=sys.stdout)
		sys.exit(0)



	###### STEP 4: Pairwise alignment of sequence pairs ######
	mapping_pairs_starts = {}
	mapping_pairs_stops = {}
	mapping_ranges = {}

	if int(options.resume_step) < 5 :
		print('# STEP 4: Pairwise alignment of sequence pairs', file=sys.stderr)
		print('[' + str(datetime.datetime.now()) + '] = STEP 4: Pairwise alignment of sequence pairs', file=sys.stdout)
		#status["4-hap2hap"] = {
		#				"4.1-map" : "TODO" ,
		#				"4.2-uniquify" : "TODO",
		#				"4.3-pairing" : "TODO"
		#				}

		if status["4-hap2hap"]["4.3-pairing"] == "TODO" :
			print('## STEP 4.1: Pairing sequences', file=sys.stderr)
			print('[' + str(datetime.datetime.now()) + '] == STEP 4.1: Pairing sequences', file=sys.stdout)
			pairwise_parameters = " -k19 -w10 -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N 5 --min-occ-floor=100 --for-only -t " + str(options.map_threads) + " "
			for pair_id in sorted(pairs.keys()) :
				hap1_id , hap2_id = pairs[pair_id]
				print('### Pair ' + str(pair_id) + ": " + hap2_id + " (Hap2) Vs. " + hap1_id + " (Hap 1)", file=sys.stderr)
				print('[' + str(datetime.datetime.now()) + '] == Pair ' + str(pair_id) + ": " + hap2_id + " (Hap2) Vs. " + hap1_id + " (Hap 1)", file=sys.stdout)
				print('#### Mapping ', file=sys.stderr)
				print('[' + str(datetime.datetime.now()) + '] === Mapping ', file=sys.stdout)
				alignment_file = files_and_folders["sequences"][hap1_id]["pair_folder"] + "/" + hap2_id + ".on." + hap1_id + ".paf"
				map_minimap( files_and_folders["sequences"][hap1_id]["fasta_file"] , [ files_and_folders["sequences"][hap2_id]["fasta_file"] ], int(options.map_threads) , pairwise_parameters , alignment_file , minimap_path , " | awk \'$5==\"+\"\' " )
				files_and_folders["sequences"][hap1_id]["pair_map_file"] = alignment_file
				files_and_folders["sequences"][hap2_id]["pair_map_file"] = alignment_file

				print('#### Uniquify alignments ', file=sys.stderr)
				print('[' + str(datetime.datetime.now()) + '] === Uniquify alignments ', file=sys.stdout)
				uniq_alignment_file = uniquify_paf(alignment_file , files_and_folders )
				files_and_folders["sequences"][hap1_id]["uniquified_map_file"] = uniq_alignment_file
				files_and_folders["sequences"][hap2_id]["uniquified_map_file"] = uniq_alignment_file

				print('#### Extract region pairing information', file=sys.stderr)
				print('[' + str(datetime.datetime.now()) + '] === Extract region pairing information', file=sys.stdout)
				pairing_file_hap1 , pairing_file_hap2 = paf2pair( uniq_alignment_file , hap1_id , hap2_id , files_and_folders)
				files_and_folders["sequences"][hap1_id]["pairing_file"] = pairing_file_hap1
				files_and_folders["sequences"][hap2_id]["pairing_file"] = pairing_file_hap2
				files_and_folders["sequences"][hap1_id]["mapping_pairs_starts"] , files_and_folders["sequences"][hap1_id]["mapping_pairs_stops"] , files_and_folders["sequences"][hap1_id]["mapping_ranges"] = write_pairs(hap1_id, files_and_folders, files_and_folders["sequences"][hap1_id]["pairing_file"])
				files_and_folders["sequences"][hap2_id]["mapping_pairs_starts"] , files_and_folders["sequences"][hap2_id]["mapping_pairs_stops"] , files_and_folders["sequences"][hap2_id]["mapping_ranges"] = write_pairs(hap2_id, files_and_folders, files_and_folders["sequences"][hap2_id]["pairing_file"])

			status["4-hap2hap"]["4.1-map"] = "DONE"
			status["4-hap2hap"]["4.2-uniquify"] = "DONE"
			status["4-hap2hap"]["4.3-pairing"] = "DONE"

		else :
			print('## STEP 4.1: Sequence pairing already performed, skipping', file=sys.stderr)
			print('[' + str(datetime.datetime.now()) + '] == STEP 4.1: Sequence pairing already performed, skipping', file=sys.stdout)
			if not check_all_paired(files_and_folders) :
				status["4-hap2hap"] = {
							"4.1-map" : "TODO" ,
							"4.2-uniquify" : "TODO",
							"4.3-pairing" : "TODO"
							}
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				print('[ERROR] - Paring information missing', file=sys.stderr)
				print('[ERROR] - Paring information missing. Please rerun the process', file=sys.stdout)
				exit(40)

		print('## STEP 4.2: Loading pairing information', file=sys.stderr)
		print('[' + str(datetime.datetime.now()) + '] == STEP 4.1: Loading pairing information', file=sys.stdout)
		mapping_pairs_starts , mapping_pairs_stops , mapping_ranges = read_pairs( files_and_folders )

		##### Save configuration file
		print('## Saving status', file=sys.stderr)
		save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

	else :
		print('# STEP 4: Pairwise alignment of sequence pairs - Loading pairing information', file=sys.stderr)
		print('[' + str(datetime.datetime.now()) + '] = STEP 4: Pairwise alignment of sequence pairs - Loading pairing information', file=sys.stdout)
		if status["4-hap2hap"]["4.3-pairing"] == "TODO" or (not check_all_paired(files_and_folders) ):
			status["4-hap2hap"] = {
						"4.1-map" : "TODO" ,
						"4.2-uniquify" : "TODO",
						"4.3-pairing" : "TODO"
						}
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
			print('[ERROR] - Paring information missing', file=sys.stderr)
			print('[ERROR] - Paring information missing. Please rerun the process', file=sys.stdout)
			exit(40)

		print('## STEP 4.1: Pairwise alignment of sequence pairs - Loading pairing information', file=sys.stderr)
		print('[' + str(datetime.datetime.now()) + '] == STEP 4.1: Pairwise alignment of sequence pairs - Loading pairing information', file=sys.stdout)
		mapping_pairs_starts , mapping_pairs_stops , mapping_ranges = read_pairs( files_and_folders )


	if int(options.stop_step) == 4 :
		print("------------------------------", file=sys.stdout)
		print("- Done", file=sys.stdout)
		print("------------------------------", file=sys.stdout)
		sys.exit(0)



	###### STEP 5: Identification of regions suitable for upgrade ######
	if int(options.resume_step) < 6 :
		##### Identify regions for each sequence
		print('[' + str(datetime.datetime.now()) + '] = STEP 5: Identification of regions suitable for upgrade', file=sys.stdout)
		print('# STEP 5: Identification of regions suitable for upgrade', file=sys.stderr)

		region_file_list = []

		if status["5-upgradeable"]["5.1-Unmatched"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] == STEP 5.1: Extracting unmatched regions', file=sys.stdout)
			print('## STEP 5.1: Extracting unmatched regions', file=sys.stderr)
			for chr in sorted(files_and_folders["sequences"].keys()) :
				if not files_and_folders["sequences"][chr]["hap"] == "U" :
					mate_id = files_and_folders["sequences"][chr]["mate_id"]
					print('[' + str(datetime.datetime.now()) + '] === ' + chr, file=sys.stdout)
					print('### ' + chr, file=sys.stderr)
					unmatched_regions_db , mate_open_edges_db = complement_regions(chr , files_and_folders , mapping_pairs_starts , mapping_ranges )

					complement_region_file_name = files_and_folders["sequences"][chr]["folder"] + "/" + chr + ".unpaired_regions.pkl.gz"
					complement_region_file = gzip.open(complement_region_file_name, 'wb')
					pickle.dump( unmatched_regions_db , complement_region_file , pickle.HIGHEST_PROTOCOL)
					complement_region_file.close()
					files_and_folders["sequences"][chr]["unpaired_regions"] = complement_region_file_name

					open_edges_file = files_and_folders["sequences"][mate_id]["folder"] + "/" + mate_id + ".open_edges.json.gz"
					with gzip.open(open_edges_file, 'wt') as f :
						json.dump( mate_open_edges_db , f , indent=4 )
					files_and_folders["sequences"][mate_id]["open_edges"] = open_edges_file
			status["5-upgradeable"]["5.1-Unmatched"]="DONE"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
		else :
			print('[' + str(datetime.datetime.now()) + '] == STEP 5.1: Validating unmatched regions', file=sys.stdout)
			print('## STEP 5.1: Validating unmatched regions', file=sys.stderr)
			if not check_unmatched_and_open_edges(files_and_folders ) :
				status["5-upgradeable"]["5.1-Unmatched"]="TODO"
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				print('[ERROR] Unpaired regions and/or open edges file missing. Status updated, rerun the script to recover', file=sys.stdout)
				print('[ERROR] Unpaired regions and/or open edges file missing. Status updated, rerun the script to recover', file=sys.stderr)
				exit(51)

		if status["5-upgradeable"]["5.2-preprocess_gaps"] == "TODO" :
			print('[' + str(datetime.datetime.now()) + '] == STEP 5.2: Paring gaps position on mate haplotype', file=sys.stdout)
			print('## STEP 5.2: Paring gaps position on mate haplotype', file=sys.stderr)
			for chr in sorted(files_and_folders["sequences"].keys()) :
				if not files_and_folders["sequences"][chr]["hap"] == "U" :
					mate_id = files_and_folders["sequences"][chr]["mate_id"]
					print('[' + str(datetime.datetime.now()) + '] === ' + chr, file=sys.stdout)
					print('### ' + chr, file=sys.stderr)
					if check_unmatched_and_open_edges(files_and_folders) :
						unmatched_regions_db = pickle.load( gzip.open( files_and_folders["sequences"][chr]["unpaired_regions"] , 'rb') )
						open_edges_db = json.load( gzip.open( files_and_folders["sequences"][chr]["open_edges"] , 'rt') )
					else :
						status["5-upgradeable"]["5.1-Unmatched"]="TODO"
						save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
						print('[ERROR] Unpaired regions and/or open edges file missing. Status updated, rerun the script to recover', file=sys.stdout)
						print('[ERROR] Unpaired regions and/or open edges file missing. Status updated, rerun the script to recover', file=sys.stderr)
						exit(52)

					raw_gap_list = bed_db_to_range_list( read_bed(files_and_folders["sequences"][chr]["gap_file"]) , chr )
					#print >> sys.stderr , gap_list
					processed_gap_file = files_and_folders["sequences"][chr]["folder"] + "/" + chr + ".processed_gaps.bed"
					gap_list , files_and_folders["sequences"][chr]["processed_gap_file"] , files_and_folders["sequences"][chr]["structure"] = preprocess_gap_list( chr , raw_gap_list , open_edges_db , files_and_folders["sequences"][chr]["length"] , processed_gap_file )
					processed_gap_db = files_and_folders["sequences"][chr]["folder"] + "/" + chr + ".processed_gaps.json.gz"
					gap_db , files_and_folders["sequences"][chr]["processed_gap_db"] = gap_mate_position(chr , gap_list , mapping_ranges , mapping_pairs_starts , unmatched_regions_db , processed_gap_db , files_and_folders )

			status["5-upgradeable"]["5.2-preprocess_gaps"] = "DONE"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
		else :
			print('[' + str(datetime.datetime.now()) + '] == STEP 5.2: Checking paring gaps position on mate haplotype results', file=sys.stdout)
			print('## STEP 5.2: Checking paring gaps position on mate haplotype results', file=sys.stderr)
			if not check_paring_gaps(files_and_folders) :
				status["5-upgradeable"]["5.2-preprocess_gaps"] = "TODO"
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				print('[ERROR] Gap paring information missing. Status updated, rerun the script to recover', file=sys.stdout)
				print('[ERROR] Gap paring information missing. Status updated, rerun the script to recover', file=sys.stderr)
				exit(52)

		if not status["5-upgradeable"]["5.3-get_sequences"] == "DONE" :
			# TODO parallelize on chromosomes
			#  Use file ".done" to manage missing analysis
			chr_list = []
			status["5-upgradeable"]["5.3-get_sequences"] = {}
			for chr in sorted(files_and_folders["sequences"].keys()) :
				if not files_and_folders["sequences"][chr]["hap"] == "U" :
					status["5-upgradeable"]["5.3-get_sequences"][chr] = "TODO"
					chr_list.append(chr)
			print('[' + str(datetime.datetime.now()) + '] == STEP 5.3: Gap flanking regions information extraction', file=sys.stdout)
			print('## STEP 5.3: Gap flanking regions information extraction', file=sys.stderr)

			for chr in sorted(chr_list) :
				region_of_interest_file = files_and_folders["sequences"][chr]["folder"] + "/" + chr + ".upgradable_regions.json.gz"
				done_file_name = region_of_interest_file + ".done"
				print('[' + str(datetime.datetime.now()) + '] === ' + chr, file=sys.stdout)
				print('### ' + chr, file=sys.stderr)
				if not os.path.exists(done_file_name) :
					mate_id = files_and_folders["sequences"][chr]["mate_id"]
					chr_sequence = read_fasta(files_and_folders["sequences"][chr]["fasta_file"])
					mate_sequence = read_fasta(files_and_folders["sequences"][mate_id]["fasta_file"])
					gap_db = json.load(open(files_and_folders["sequences"][chr]["processed_gap_db"]))
					mate_gap_db = json.load(open(files_and_folders["sequences"][mate_id]["processed_gap_db"]))

					region_file_list.append( extract_sequence_and_signals( chr , mate_id , files_and_folders , gap_db , mate_gap_db , region_of_interest_file , chr_sequence , mate_sequence , int(options.flanking) ) )

					files_and_folders["sequences"][chr]["region_of_interest_file"] = region_of_interest_file
					done_file = open(done_file_name , 'w+')
					done_file.close()
				else :
					print('#### Sequence already processed, skipping', file=sys.stderr)

			status["5-upgradeable"]["5.3-get_sequences"] = "DONE"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
		else :
			print('[' + str(datetime.datetime.now()) + '] == STEP 5.3: Gap flanking regions information validation', file=sys.stdout)
			print('## STEP 5.3: Gap flanking regions information validation', file=sys.stderr)
			missing = check_extracted_sequences( files_and_folders )
			if not missing == [] :
				status["5-upgradeable"]["5.3-get_sequences"] = "TODO"
				save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
				print('[ERROR] Gap paring information missing. Status updated, rerun the script to recover', file=sys.stdout)
				print('[ERROR] Gap paring information missing. Status updated, rerun the script to recover', file=sys.stderr)
				print("[ERROR] Sequences to be processed: " + " ; ".join([ str(x) for x in missing ]), file=sys.stderr)
				exit(53)

		##### Save configuration file
		print('## Saving status', file=sys.stderr)
		save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

	else :
		print('[' + str(datetime.datetime.now()) + '] = STEP 5: Identification of regions suitable for upgrade', file=sys.stdout)
		print('# STEP 5: Identification of regions suitable for upgrade, skipping', file=sys.stderr)
		print('[' + str(datetime.datetime.now()) + '] == STEP 5.1: Gap flanking regions information validation', file=sys.stdout)
		print('## STEP 5.1: Gap flanking regions information validation', file=sys.stderr)
		missing = check_extracted_sequences( files_and_folders )
		if not missing == [] :
			status["5-upgradeable"]["5.3-get_sequences"] = "TODO"
			print('[ERROR] Gap paring information missing. Status updated, rerun the script to recover', file=sys.stdout)
			print('[ERROR] Gap paring information missing. Status updated, rerun the script to recover', file=sys.stderr)
			print("[ERROR] Sequences to be processed: " + " ; ".join([ str(x) for x in missing ]), file=sys.stderr)

			print('[ERROR] Checking intermediate steps', file=sys.stdout)
			print('[ERROR] Checking intermediate steps', file=sys.stderr)
			if not check_paring_gaps(files_and_folders) :
				status["5-upgradeable"]["5.2-preprocess_gaps"] = "TODO"
				if not check_unmatched_and_open_edges(files_and_folders) :
					status["5-upgradeable"]["5.1-Unmatched"] = "TODO"
			else :
				status["5-upgradeable"]["5.1-Unmatched"] = "DONE"
			print('Status updated, rerun the script to recover', file=sys.stdout)
			print('Status updated, rerun the script to recover', file=sys.stderr)
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
			exit(5)


	if int(options.stop_step) == 5:
		print('## Saving status', file=sys.stderr)
		save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
		print("------------------------------", file=sys.stdout)
		print("- Done", file=sys.stdout)
		print("------------------------------", file=sys.stdout)
		sys.exit(0)



	###### STEP 6: Patch by mapping unplaced sequences or by filling with sequence on the alternative haplotype ######
	# status["6-filling"] = {
	#				"6.1-gather": "TODO" ,
	#				"6.2-map" : "TODO" ,
	#				"6.3-select" : "TODO"
	#				}

	##### Identify regions for each sequence
	print('[' + str(datetime.datetime.now()) + '] = STEP 6: Gap patching ', file=sys.stdout)
	print('# STEP 6: Gap patching', file=sys.stderr)
	print('[' + str(datetime.datetime.now()) + '] == STEP 6.1: Gathering gaps information ', file=sys.stdout)
	print('## STEP 6.1: Gathering gaps information', file=sys.stderr)
	# Load gap files and split the gap instances in different folders
	if status["6-filling"]["6.1-gather"] == "TODO" :

		region_file_list = collect_files_names( files_and_folders , "region_of_interest_file" )

		gap_db = {}
		gap_db_file = temp_folder + "/gap_db.json.gz"
		gap_flanking_sequences = {}
		gap_flanking_signal = {}
		gap_alt_sequences = {}
		gap_alt_signal = {}
		gap_folder = temp_folder + "/gaps"
		mkdir(gap_folder)
		gap_id = 0
		for gap_file in region_file_list :
			print("### Processing " + gap_file, file=sys.stderr)
			gap_descriptors_list = json.load( gzip.open( gap_file , 'rt' ) )
			#json.dump( gap_descriptors_list , sys.stderr , indent=4 )
			for el in sorted(gap_descriptors_list, key = lambda i: i['gap_region'])  :
				#	el["seq_id"] = seq_id
				#	el["mate_id"] = mate_id
				#	el["gap_region"] = [ gap_start , gap_stop ]
				#	el["flanking_upstream_content"] = flanking_upstream_content
				#	el["flanking_upstream_region"] = flanking_upstream_region
				#	el["flanking_downstream_content"] = flanking_downstream_content
				#	el["flanking_downstream_region"] = flanking_downstream_region
				#	el["gap_corr_region"] = [ gap_corr_start , gap_corr_stop ]
				#	el["gap_corr_region_content"] = corr_region_content
				#	el["upstream_corr_region_content"] = upstream_corr_region_content
				#	el["upstream_corr_region_region"] = upstream_corr_region_region
				#	el["downstream_corr_region_content"] = downstream_corr_region_content
				#	el["downstream_corr_region_region"] = downstream_corr_region_region
				#	el["patching_strategy"]
				#	#	el["patching_strategy"]["map_gap"] = "hybrid_left:alt:right" |...
				#	#	el["patching_strategy"]["map_alt"] = "rep:alt:right" | ...
				#	el["fasta_sequences"]
				#	#	el["fasta_sequences"]["flanking_upstream"] = str( seq_fasta[ int(flanking_upstream_region[0]) : int(flanking_upstream_region[1]) ] )
				#	#	el["fasta_sequences"]["flanking_downstream"] = str(seq_fasta[ int(flanking_downstream_region[0]) : int(flanking_downstream_region[1]) ])
				#	#	el["fasta_sequences"]["corr_seq"] = str(mate_fasta[ int(gap_corr_start) : int(gap_corr_stop) ])
				#	#	el["fasta_sequences"]["corr_seq_upstream"] = str(mate_fasta[ int(upstream_corr_region_region[0]) : int(upstream_corr_region_region[1]) ])
				#	#	el["fasta_sequences"]["corr_seq_downstream"] = str(mate_fasta[ int(downstream_corr_region_region[0]) : int(downstream_corr_region_region[1]) ])
				#	el["sequences_characteristics"]
				#	#	el["sequences_characteristics"]["flanking_upstream_content"] = analyze_signal( flanking_upstream_content )
				#	#	el["sequences_characteristics"]["flanking_downstream_content"] = analyze_signal( flanking_downstream_content )
				#	#	el["sequences_characteristics"]["gap_corr_region_content"] = analyze_signal( corr_region_content )
				#	#	el["sequences_characteristics"]["upstream_corr_region_content"] = analyze_signal( upstream_corr_region_content )
				#	#	el["sequences_characteristics"]["downstream_corr_region_content"] = analyze_signal( downstream_corr_region_content )

				# Generate folder for gap
				gap_id += 1
				print("#### Generating gap " + str(gap_id), file=sys.stderr)
				gap_instance_dir = gap_folder + "/" + str(gap_id)
				gap_instance_file = gap_instance_dir + "/" + "gap_instance.json.gz"
				mkdir(gap_instance_dir)

				# Extract sequence and signal information for every gap
				chr = el["seq_id"]
				start , stop = el["gap_region"]
				for part in list(files_and_folders["sequences"][chr]["structure"].keys()) :
					if (int(files_and_folders["sequences"][chr]["structure"][part][0]) == int(start)) and (int(files_and_folders["sequences"][chr]["structure"][part][1]) == int(stop) ):
						files_and_folders["sequences"][chr]["structure"][part][2] = gap_id

				el["patch"] = make_sequences_and_signals( el , gap_id , files_and_folders , gap_instance_dir )
				# Save gap fasta sequence and signal in db for mapping
				gap_flanking_sequences[ el["patch"]["target_1"]["id"] ] = el["patch"]["target_1"]["sequence"]
				gap_flanking_signal[ el["patch"]["target_1"]["id"] ] = el["patch"]["target_1"]["signal"]
				gap_alt_sequences[ el["patch"]["target_2"]["id"] ] = el["patch"]["target_2"]["sequence"]
				gap_alt_signal[ el["patch"]["target_2"]["id"] ] = el["patch"]["target_2"]["signal"]

				# Save tracking info in gap db
				gap_db[gap_id] = {
					"dir" : gap_instance_dir ,
					"file": gap_instance_file ,
					"sequence_id" : chr ,
					"coordinates" : el["gap_region"] ,
					"flanking_upstream_region" : el["flanking_upstream_region"] ,
					"flanking_downstream_region" : el["flanking_downstream_region"] ,
					"corresponding_sequence_id" : el["mate_id"]  ,
					"corresponding_coordinates" : el["gap_corr_region"] ,
					"upstream_corr_region_region" : el["upstream_corr_region_region"] ,
					"downstream_corr_region_region" : el["downstream_corr_region_region"]
				}
				gap_db[gap_id]["target_1_id"] = el["patch"]["target_1"]["id"]
				gap_db[gap_id]["target_1_strategy"] = el["patching_strategy"]["map_gap"] # "hybrid_left:alt:right"
				gap_db[gap_id]["target_2_id"] = el["patch"]["target_2"]["id"]
				gap_db[gap_id]["target_2_strategy"] = el["patching_strategy"]["map_alt"] # "left:alt:right"

				# Save info in gap folder
				with gzip.open( gap_instance_file , 'wt' ) as f :
					json.dump( el , f , indent=4 )

		print("### Saving results", file=sys.stderr)
		with gzip.open( gap_db_file , 'wt' ) as f :
			json.dump( gap_db , f , indent=4 , sort_keys=True)
		files_and_folders["gap_db_file"] = gap_db_file
		save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

		# make gap sequences fasta
		gap_flanking_sequences_fasta = temp_folder + "/gap_flanking_sequences.fasta"
		gap_alt_sequences_fasta = temp_folder + "/gap_alt_sequences.fasta"
		gap_flanking_signal_file = temp_folder + "/gap_flanking_sequences.signal.json"
		gap_alt_signal_file = temp_folder + "/gap_alt_sequences.signal.json"
		gap_flanking_sequences_fasta = write_fasta_from_db( gap_flanking_sequences , gap_flanking_sequences_fasta , compress = False )
		gap_alt_sequences_fasta = write_fasta_from_db( gap_alt_sequences , gap_alt_sequences_fasta , compress = False )
		json.dump( gap_flanking_signal , open( gap_flanking_signal_file , 'w+' ) , indent=4 , sort_keys=True)
		json.dump( gap_alt_signal , open( gap_alt_signal_file , 'w+' ) , indent=4 , sort_keys=True)
		touch(temp_folder + "/gap_sequences_and_singals.done")
		status["6-filling"]["6.1-gather"] = "DONE"
		save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
	else :
		# Check if the necessary files are present, otherwise reset status and exit
		print('[' + str(datetime.datetime.now()) + '] === Processing of gaps regions already performed, validating files integrity and loading results', file=sys.stdout)
		print('### Processing of gaps regions already performed, validating files integrity and loading results', file=sys.stderr)
		gap_flanking_sequences_fasta = temp_folder + "/gap_flanking_sequences.fasta"
		gap_alt_sequences_fasta = temp_folder + "/gap_alt_sequences.fasta"
		gap_flanking_signal_file = temp_folder + "/gap_flanking_sequences.signal.json"
		gap_alt_signal_file = temp_folder + "/gap_alt_sequences.signal.json"
		gap_db_file = temp_folder + "/gap_db.json.gz"

		if os.path.exists(temp_folder + "/gap_sequences_and_singals.done") and \
				os.path.exists(gap_flanking_sequences_fasta) and \
				os.path.exists(gap_alt_sequences_fasta) and \
				os.path.exists(gap_flanking_signal_file) and \
				os.path.exists(gap_alt_signal_file) and \
				os.path.exists(gap_db_file) and \
				os.path.exists(temp_folder + "/gaps") :

			gap_flanking_sequences = read_fasta(temp_folder + "/gap_flanking_sequences.fasta")
			gap_alt_sequences = read_fasta(temp_folder + "/gap_alt_sequences.fasta")
			gap_flanking_signal = json.load( open( gap_flanking_signal_file ) )
			gap_alt_signal = json.load( open( gap_alt_signal_file ) )
			gap_db = json.load( gzip.open( gap_db_file , 'rt' ) )

		else:
			status["6-filling"] = {
				"6.1-gather": "TODO" ,
				"6.2-map" : "TODO" ,
				"6.3-select" : "TODO"
				}
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

			print('[ERROR] Gap strucutre information missing. Status updated, rerun the script to recover', file=sys.stdout)
			print('[ERROR] Gap strucutre information missing. Status updated, rerun the script to recover', file=sys.stderr)
			exit(61)

	# map and uniquify
	if status["6-filling"]["6.2-map"] == "TODO" :
		print('[' + str(datetime.datetime.now()) + '] == STEP 6.2: Align unplaced on gaps', file=sys.stdout)
		print('## STEP 6.2: Align unplaced on gaps', file=sys.stderr)
		unplaced_on_gap_mappings_file = temp_folder + "/unplaced_on_gap_mappings.json.gz"
		if os.path.exists(unplaced_on_gap_mappings_file) :
			unplaced_on_gap_mappings = json.load( gzip.open( unplaced_on_gap_mappings_file , "rt") )
		else :
			unplaced_on_gap_mappings = {}

		unplaced_on_gap_mappings = map_on_gap( unplaced_on_gap_mappings , unplaced_on_gap_mappings_file , gap_flanking_sequences_fasta , gap_alt_sequences_fasta ,  gap_flanking_signal , gap_alt_signal, files_and_folders , int(options.map_threads) , temp_folder , gap_db , unwanted_pairs , known_chr_by_seqid )

		files_and_folders["unplaced_on_gap_mappings_file"] = unplaced_on_gap_mappings_file
		touch(unplaced_on_gap_mappings_file + ".done")
		status["6-filling"]["6.2-map"] = "DONE"
		save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

		# unplaced_on_gap_mappings
		#	#	unplaced_on_gap_mappings["target_1"] =
		#	#	unplaced_on_gap_mappings["target_X"]
		#	#	#	unplaced_on_gap_mappings["target_2"]["map_file"]
		#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"] = {}
		#	#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"][Tid] = {}
		#	#	#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"][Tid]["best_match"] = Qid
		#	#	#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"][Tid]["best_match_len"] = Q_len
		#	#	#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"][Tid]["longest_pass_path"] = longest_pass_path
		#	#	#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"][Tid]["longest_pass_graph"] = longest_pass_graph
		#	#	#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"][Tid]["longest_pass_path_matches"] = longest_pass_path_matches
		#	#	#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"][Tid]["longest_all_path"] = longest_all_path
		#	#	#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"][Tid]["longest_all_graph"] = longest_all_graph
		#	#	#	#	#	unplaced_on_gap_mappings["target_2"]["map_info"][Tid]["longest_all_path_matches"] = longest_all_path_matches
	else :
		# Check if the necessary files are present, otherwise reset status and exit
		print('[' + str(datetime.datetime.now()) + '] === Alignment already performed, checking files integrity', file=sys.stdout)
		print('### Alignment already performed, checking files integrity', file=sys.stderr)
		if "unplaced_on_gap_mappings_file" in files_and_folders and \
				os.path.exists(files_and_folders["unplaced_on_gap_mappings_file"]) and \
				os.path.exists(files_and_folders["unplaced_on_gap_mappings_file"] + ".done") :
			unplaced_on_gap_mappings_file = files_and_folders["unplaced_on_gap_mappings_file"]
			unplaced_on_gap_mappings = json.load( gzip.open( unplaced_on_gap_mappings_file , "rt") )
		else :
			status["6-filling"]["6.2-map"] = "TODO"
			status["6-filling"]["6.3-select"] = "TODO"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
			print('[ERROR] Alignments of unplaced on gaps is missing. Status updated, rerun the script to recover', file=sys.stdout)
			print('[ERROR] Alignments of unplaced on gaps is missing. Status updated, rerun the script to recover', file=sys.stderr)
			exit(62)

	### Select best patch for each gap
	### -> If none found: Trigger for diploid confirmation -> homozygous fill
	if status["6-filling"]["6.3-select"] == "TODO" :
		print('[' + str(datetime.datetime.now()) + '] == STEP 6.3: Selecting best filler', file=sys.stdout)
		print('## STEP 6.3: Selecting best filler', file=sys.stderr)
		findings_file = options.output + ".gap_filling_findings.txt"
		findings = open(findings_file , 'w+')
		for gap_id in sorted( gap_db.keys() ) :
			print('### Processing Gap:' + str(gap_id), file=sys.stderr)

			T1_sequence_id = gap_db[gap_id]["sequence_id"]
			T1_coordinates = gap_db[gap_id]["coordinates"]
			T1_upstream = gap_db[gap_id]["flanking_upstream_region"]
			T1_downstream = gap_db[gap_id]["flanking_downstream_region"]

			T2_sequence_id = gap_db[gap_id]["corresponding_sequence_id"]
			T2_coordinates = gap_db[gap_id]["corresponding_coordinates"]
			T2_upstream = gap_db[gap_id]["upstream_corr_region_region"]
			T2_downstream = gap_db[gap_id]["downstream_corr_region_region"]

			T1_id = gap_db[gap_id]["target_1_id"]
			T1_strategy = gap_db[gap_id]["target_1_strategy"]
			T2_id = gap_db[gap_id]["target_2_id"]
			T2_strategy = gap_db[gap_id]["target_2_strategy"]

			if T1_id in unplaced_on_gap_mappings["target_1"]["map_info"] :
				T1_best_match = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["best_match"][0:-2]
				T1_best_match_sequence_len = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["best_match_sequence_len"]
				T1_best_match_strand = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["best_match"][-1]
				T1_longest_pass_path = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_pass_path"]
				T1_longest_all_path = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_all_path"]
				T1_longest_pass_path_matches = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_pass_path_matches"]
				T1_longest_all_path_matches  = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_all_path_matches"]
			else :
				T1_best_match = ""
				T1_best_match_sequence_len = 0
				T1_best_match_strand = ""
				T1_longest_pass_path = ""
				T1_longest_all_path = ""
				T1_longest_pass_path_matches = 0
				T1_longest_all_path_matches = 0

			if T2_id in unplaced_on_gap_mappings["target_2"]["map_info"] :
				T2_best_match = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["best_match"][0:-2]
				T2_best_match_sequence_len = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["best_match_sequence_len"]
				T2_best_match_strand = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["best_match"][-1]
				T2_longest_pass_path = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_pass_path"]
				T2_longest_all_path = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_all_path"]
				T2_longest_pass_path_matches = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_pass_path_matches"]
				T2_longest_all_path_matches  = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_all_path_matches"]
			else:
				T2_best_match = ""
				T2_best_match_sequence_len = 0
				T2_best_match_strand = ""
				T2_longest_pass_path = ""
				T2_longest_all_path = ""
				T2_longest_pass_path_matches = 0
				T2_longest_all_path_matches = 0

			print("> Gap " + str(gap_id), file=findings)
			print("# Gap side (" + T1_id + ")", file=findings)
			print("## Sequence id                  :" + str(T1_sequence_id), file=findings)
			print("## Coordinates                  :" + str(T1_coordinates), file=findings)
			print("## Upstream flanking region     :" + str(T1_upstream), file=findings)
			print("## Downstream flanking region   :" + str(T1_downstream), file=findings)
			print("## T1_strategy                  :" + str(T1_strategy), file=findings)
			print("## T1_best_match                :" + str(T1_best_match ), file=findings)
			print("## T1_best_match_sequence_len   :" + str(T1_best_match_sequence_len ), file=findings)
			print("## T1_best_match_strand         :" + str(T1_best_match_strand ), file=findings)
			print("## T1_longest_pass_path         :" + str(T1_longest_pass_path ), file=findings)
			print("## T1_longest_all_path          :" + str(T1_longest_all_path ), file=findings)
			print("## T1_longest_pass_path_matches :" + str(T1_longest_pass_path_matches ), file=findings)
			print("## T1_longest_all_path_matches  :" + str(T1_longest_all_path_matches ), file=findings)

			print("# Alternative allele side (" + T2_id + ")", file=findings)
			print("## Sequence id                  :" + str(T2_sequence_id), file=findings)
			print("## Coordinates                  :" + str(T2_coordinates), file=findings)
			print("## Upstream flanking region     :" + str(T2_upstream), file=findings)
			print("## Downstream flanking region   :" + str(T2_downstream), file=findings)
			print("## T2_best_match                :" + str(T2_best_match ), file=findings)
			print("## T2_strategy                  :" + str(T2_strategy), file=findings)
			print("## T2_best_match_sequence_len   :" + str(T2_best_match_sequence_len ), file=findings)
			print("## T2_best_match_strand         :" + str(T2_best_match_strand ), file=findings)
			print("## T2_longest_pass_path         :" + str(T2_longest_pass_path ), file=findings)
			print("## T2_longest_all_path          :" + str(T2_longest_all_path ), file=findings)
			print("## T2_longest_pass_path_matches :" + str(T2_longest_pass_path_matches ), file=findings)
			print("## T2_longest_all_path_matches  :" + str(T2_longest_all_path_matches  ), file=findings)
			print("", file=findings)

			longest_pass_path_matches_delta = int(T2_longest_pass_path_matches) - int(T1_longest_pass_path_matches)
			longest_all_path_matches_delta  = int(T2_longest_all_path_matches) - int(T1_longest_all_path_matches)
			use_filler = "NONE"
			if float(T1_best_match_sequence_len) > 0 :
				T1_coverage = 100 * ( float(T1_longest_pass_path_matches) / float(T1_best_match_sequence_len) )
			else :
				T1_coverage = 0
			if float(T2_best_match_sequence_len) > 0 :
				T2_coverage = 100 * ( float(T2_longest_pass_path_matches) / float(T2_best_match_sequence_len) )
			else :
				T2_coverage = 0

			if T1_id == T2_id :
				use_filler = "T1"
			elif ( T1_strategy in [
									"hybrid_alt:rep" ,
									"hybrid_alt:right" ,
									"hybrid_left:alt" ,
									"hybrid_left:alt:rep" ,
									"hybrid_left:alt:right" ,
									"hybrid_rep:alt" ,
									"hybrid_rep:alt:rep" ,
									"hybrid_rep:alt:right"
								  ]
				) and (T1_id in unplaced_on_gap_mappings["target_1"]["map_info"]) :
				if longest_pass_path_matches_delta <= (0.5 * float(T1_longest_pass_path_matches)) :
					use_filler = "T1"
				elif longest_all_path_matches_delta <= (0.5 * float(T1_longest_all_path_matches)) :
					use_filler = "T1"
				else :
					use_filler = "T2"
			elif T2_id in unplaced_on_gap_mappings["target_2"]["map_info"] :
				if longest_pass_path_matches_delta <= (-0.5 * float(T2_longest_pass_path_matches)) :
					use_filler = "T1"
				elif longest_all_path_matches_delta <= (-0.5 * float(T2_longest_all_path_matches)) :
					use_filler = "T1"
				else :
					use_filler = "T2"
			elif T1_id in unplaced_on_gap_mappings["target_1"]["map_info"] :
				use_filler = "T1"
			else :
				use_filler = "NONE"

			# Minimum matching coverage percentage filter

			if not use_filler == "T1" and T1_coverage < float(options.coverage) :
				use_filler = "NONE"

			if not use_filler == "T2" and T2_coverage < float(options.coverage) :
				use_filler = "NONE"

			print('#### Selected supporting sequence: ' + use_filler, file=sys.stderr)

			if use_filler == "T2" :
				gap_db[gap_id]["best_match"] = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["best_match"][0:-2]
				gap_db[gap_id]["best_match_strand"] = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["best_match"][-1]
				gap_db[gap_id]["best_match_region"] = [ 0 , int(unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["best_match_sequence_len"]) ]
				gap_db[gap_id]["longest_pass_path"] = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_pass_path"]
				#gap_db[gap_id]["longest_pass_graph"] = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_pass_graph"]
				gap_db[gap_id]["longest_pass_path_matches"] = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_pass_path_matches"]
				gap_db[gap_id]["longest_all_path"] = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_all_path"]
				#gap_db[gap_id]["longest_all_graph"] = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_all_graph"]
				gap_db[gap_id]["longest_all_path_matches"]  = unplaced_on_gap_mappings["target_2"]["map_info"][T2_id]["longest_all_path_matches"]
			elif use_filler == "T1" :
				gap_db[gap_id]["best_match"] = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["best_match"][0:-2]
				gap_db[gap_id]["best_match_strand"] = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["best_match"][-1]
				gap_db[gap_id]["best_match_region"] = [ 0 , int(unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["best_match_sequence_len"]) ]
				gap_db[gap_id]["longest_pass_path"] = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_pass_path"]
				#gap_db[gap_id]["longest_pass_graph"] = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_pass_graph"]
				gap_db[gap_id]["longest_pass_path_matches"] = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_pass_path_matches"]
				gap_db[gap_id]["longest_all_path"] = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_all_path"]
				#gap_db[gap_id]["longest_all_graph"] = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_all_graph"]
				gap_db[gap_id]["longest_all_path_matches"]  = unplaced_on_gap_mappings["target_1"]["map_info"][T1_id]["longest_all_path_matches"]
			else :
				print('#### Checking alternative haplotype for diploid homozygous state: ', file=sys.stderr)
				# Assess diploid of alternative region
				# Disable adding the homozygous filler with options.no_homozygous
				if not options.no_homozygous :
					print('##### Option Enabled', file=sys.stderr)
					gap_info = json.load( gzip.open( gap_db[gap_id]["file"] , "rt") )
					corr_region_content = gap_info["sequences_characteristics"]["gap_corr_region_content"]
					if corr_region_content in ( "OK" , "DIP" ) :
						print('##### Homozygous region, filling with ' + str(gap_info["gap_corr_region"]), file=sys.stderr)
						# Alternative sequence has diploid coverage, set it as filler
						gap_db[gap_id]["best_match"] = gap_info["mate_id"]
						gap_db[gap_id]["best_match_region"] = gap_info["gap_corr_region"]
						gap_db[gap_id]["best_match_strand"] = "+"
					else :
						print('##### Alternative haplotype not diploid, homozygous, or reliable. No filling will be used.', file=sys.stderr)
				else :
					print('##### Option disabled', file=sys.stderr)

		findings.close()
		with gzip.open( gap_db_file , 'wt' ) as f :
			json.dump( gap_db , f , indent=4 , sort_keys=True)
		print('[' + str(datetime.datetime.now()) + '] === Saving resultd in BLOCK format for HaploMaker', file=sys.stdout)
		print('### Saving resultd in BLOCK format for HaploMaker', file=sys.stderr)
		# Save file to use as input of HaploMake (block format)
		haplomake_file_name = options.output + ".structure.block"
		# block file format:
		#	#	>OutSeqID_1
		#	#	SeqID_1	Start	Stop	Strand
		#	#	SeqID_2	Start	Stop	Strand
		#	#	...
		#	#	>OutSeqID_2
		#	#	SeqID_N	Start	Stop	Strand
		#	#	SeqID_M	Start	Stop	Strand
		#	#	...
		haplomake_file = open( haplomake_file_name , 'w' )
		used_seq_ids = []
		gap_counter = 0
		filled_counter = 0
		for chr in sorted(files_and_folders["sequences"].keys()) :
			if files_and_folders["sequences"][chr]["hap"] == "U" :
				continue
			used_seq_ids.append(chr)
			print(">" + str(chr), file=haplomake_file)
			for part_num in sorted([ int(x) for x in list(files_and_folders["sequences"][chr]["structure"].keys()) ]) :
				part_start , part_stop , part_id = files_and_folders["sequences"][chr]["structure"][str(part_num)]
				if part_id == "seq" :
					# component is the original sequence
					# Print it directly
					print("\t".join(str(x) for x in [chr , part_start , part_stop , "+"] ), file=haplomake_file)
				else :
					gap_counter += 1
					# component is a gap
					# Extract filling information and print block accordingly
					if str(part_id) in gap_db :
						gap_id = str(part_id)
					else :
						if int(part_id) in gap_db :
							gap_id = int(part_id)
						#else :
						#	# write error
					if "best_match" in gap_db[gap_id] :
						fill_name = gap_db[gap_id]["best_match"]
						fill_start , fill_stop = gap_db[gap_id]["best_match_region"]
						fill_strand = gap_db[gap_id]["best_match_strand"]
						used_seq_ids.append(fill_name)
						filled_counter += 1
						print("\t".join(str(x) for x in [fill_name , fill_start , fill_stop , fill_strand] ), file=haplomake_file)
		unplaced_counter = 0
		placed_unplaced_counter = 0
		# write sequences not used for patching purpose
		for chr in sorted(files_and_folders["sequences"].keys()) :
			if files_and_folders["sequences"][chr]["hap"] == "U" :
				unplaced_counter += 1
			if chr in used_seq_ids :
				placed_unplaced_counter += 1
				continue
			else:
				print(">" + chr, file=haplomake_file)
				print("\t".join(str(x) for x in [chr , 0 , files_and_folders["sequences"][chr]["length"] , "+"] ), file=haplomake_file)
		haplomake_file.close()

		print('#### Number of gaps to fill: ' + str(gap_counter), file=sys.stderr)
		print('#### Number of gaps filled: ' + str(filled_counter), file=sys.stderr)
		print('#### Number of original unplaced sequences: ' + str(unplaced_counter), file=sys.stderr)
		print('#### Number of unplaced sequences used as filler: ' + str(placed_unplaced_counter), file=sys.stderr)

		print('## Saving status', file=sys.stderr)
		###### Save configuration file
		print('## Saving status', file=sys.stderr)
		status["6-filling"]["6.3-select"] = "DONE"
		save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)

	else :
		print('[' + str(datetime.datetime.now()) + '] == STEP 6.3: Best filler selection already performed, checking files', file=sys.stdout)
		print('## STEP 6.3: Best filler selection already performed, checking files', file=sys.stderr)
		# Check if output for HaploMake is present, otherwise reset status and exit
		haplomake_file = options.output + ".structure.block"
		findings_file = options.output + ".gap_filling_findings.txt"
		if os.path.exists(haplomake_file + ".done" ) and os.path.exists(findings_file + ".done" ) :
			print('[' + str(datetime.datetime.now()) + '] === All result files already done, nothing to do', file=sys.stdout)
			print('### All result files already done, nothing to do', file=sys.stderr)
		else :
			status["6-filling"]["6.3-select"] = "TODO"
			save_status(files_and_folders, conf_file, pairs, pairs_file, status, status_file)
			print('[' + str(datetime.datetime.now()) + '] === Result files may be incomplete. Status updated, rerun the script to recover', file=sys.stdout)
			print('### Result files may be incomplete. Status updated, rerun the script to recover', file=sys.stderr)
			exit(63)

	###### END ######

	print("------------------------------", file=sys.stdout)
	print("- Done", file=sys.stdout)
	print("------------------------------", file=sys.stdout)
	print("##############################", file=sys.stderr)
	print("# Done", file=sys.stderr)
	print("##############################", file=sys.stderr)



if __name__ == '__main__':
	main()
