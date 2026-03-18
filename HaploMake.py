#!/usr/bin/env python

import argparse
from lib_files.HaploFunct import *
from lib_files.GFF_lib import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


def main() :

	###### Options and help ######
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta", dest="fasta",
					help="FASTA file(s) with all input sequences", metavar="genome.fasta [Required]")
	parser.add_argument("-s", "--structure", dest="structure",
					help="Text file describing the novel sequence structure to build. [REQUIRED]", metavar="structure.txt")
	parser.add_argument("--format", dest="format", default="BLOCK" ,
					help="Structure file format [Default: BLOCK]", metavar="[BLOCK|AGP|BED]")

	parser.add_argument("-o", "--out", dest="out", default="out" ,
					help="Output file names [Default: out]", metavar="outname")
	parser.add_argument("-p", "--prefix", dest="prefix", default="" ,
					help="New sequence names prefix to be added before the old name. Leave unset to keep the same names [Default: NONE]", metavar="new")

	parser.add_argument("-g", "--gff3", dest="gff3",
					help="Genome annotation in GFF3 format.", metavar="annotation.gff3")
	parser.add_argument("-b", "--bed", dest="bed",
					help="Region annotation in BED format. Accepted BED3 and BED6 formats, additional columns will be considered as annotation and reported unedited", metavar="annotation.bed")

	parser.add_argument("-a", "--agp", dest="agp",
					help="AGP file defining the coordinates of original sequences composing actual assembly", metavar="previous_to_actual.genome.agp")
	parser.add_argument("-m" , "--mapper", dest="mapper", default="blat",
					#help="Mapping tool to use [Default: blat]" , metavar="[blat|nucmer]")
					help=argparse.SUPPRESS )
	parser.add_argument("-c" , "--cores", dest="cores", default="1",
					#help="[For nucmer] Number of cores to use for mapping [Default: 1]" , metavar="N")
					help=argparse.SUPPRESS )
	parser.add_argument("--skipoverlap", dest="skipoverlap", default=False, action="store_true",
					help="Skip the search of overlap between flanking blocks")
	parser.add_argument("--overhang" , dest="overhang", default="10000" ,
					#help="Maximum overhang sequence to use for overlap search [default 10000]", metavar="N")
					help=argparse.SUPPRESS )
	parser.add_argument("--gap", dest="gap_size", default="1000",
					help="Size of the residual gap around inserted sequences if no overlap is detected between flanking regions of consecutive blocks [Default: 1,000bp]", metavar="N")
	parser.add_argument("--spacer", dest="spacer", default="10" ,
					help="Size of the spacer gap inserted between trimmed sequences when there is an overlap between flanking regions [Default: 10bp]", metavar="N")
	parser.add_argument("--ignoreids", dest="ignoreids", default=False, action="store_true",
					help="Ignore output sequence ids reported in structure file, use [-p|--prefix] with progressive numbers instead [Default for BED input]")

	parser.add_argument("--reverse", dest="mode", default=False, action="store_true",
					help="[For AGP structure files] Reverse direction of feature extraction, from new to old. FASTA input must contain the genomic sequences of the new version")

	parser.add_argument("-u" , "--unplaced" , dest="add_unplaced" , default=False, action="store_true",
					help="Add unplaced sequences into output files" )

	parser.add_argument("--noagp", dest="noagp", default=False, action="store_true",
					help="Avoid printing the AGP file output")
	parser.add_argument("--noprint", dest="noprint", default=False, action="store_true",
					help="Avoid printing the FASTA file output")

	print("Running HaploMake tool from HaploSync version " + get_version(), file=sys.stdout)
	print("To reproduce this run use the following command: " + " ".join( pipes.quote(x) for x in sys.argv), file=sys.stdout)
	print("----", file=sys.stdout)

	# Sanity Check

	paths = set_paths(os.path.join(sys.path[0], 'HaploSync.conf.toml'))
	minimap_path = paths["minimap2"]
	samtools_path = paths["samtools"]
	bedtools_path = paths["bedtools"]
	nucmer_path = paths["nucmer"]
	showcoords_path = paths["show-coords"]
	blat_path = paths["blat"]
	scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/support_scripts"

	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()
	gap_size = int(options.gap_size)

	if not options.fasta :
		print("[ERROR] Genome FASTA file missing", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	if not options.structure :
		print("[ERROR] Structure file missing", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	if (not options.skipoverlap) and options.mode :
		print("[ERROR] Reverse use of AGP structure [--reverse] incompatible with overlap search and correction.", file=sys.stderr)
		sys.exit(1)

	if options.noprint :
		no_fasta = True
	else :
		no_fasta = False

	if not options.gff3 :
		print("[MEMO] Annotation file missing: Coordinates conversion will not take place. You'll need to run it downstream of this process", file=sys.stdout)

	#if not options.agp :
	#	print >> sys.stdout, "[MEMO] No AGP file was given. Structure of original sequences will not be available"

	if options.skipoverlap :
		print("[MEMO] Adjacent sequences overlap analysis will be skipped. Sequence will be inserted entirely separated by gaps of " + str(gap_size) + "bp in length.\nThis procedure may allow the use of duplicated genomic content", file=sys.stderr)

	# Read inputs
	print('[' + str(datetime.datetime.now()) + '] = Read inputs', file=sys.stdout)
	print('# Read inputs', file=sys.stderr)
	print('[' + str(datetime.datetime.now()) + '] == Loading FASTA sequences', file=sys.stdout)
	print('## Loading FASTA sequences', file=sys.stderr)
	fasta_files = options.fasta
	fasta_files_list = fasta_files.split(",")
	fasta_dict = {}
	for file_name in fasta_files_list :
		print('### Loading ' + file_name, file=sys.stderr)
		if fasta_dict == {} :
			fasta_dict = read_fasta(file_name)
		else :
			fasta_dict.update(read_fasta(file_name))
	fasta_len_dict = get_length_from_fasta_db(fasta_dict)

	print('[' + str(datetime.datetime.now()) + '] == Listing structure files', file=sys.stdout)
	print('## Listing structure files', file=sys.stderr)
	structure_files = options.structure
	structure_files_list = structure_files.split(",")

	if options.gff3 :
		print('[' + str(datetime.datetime.now()) + '] == Loading GFF3 annotation', file=sys.stdout)
		print('## Loading GFF3 annotation', file=sys.stderr)
		annotation_gff3 , mRNA_db = read_gff3(options.gff3)
	else :
		annotation_gff3 = ""

	if options.bed :
		print('[' + str(datetime.datetime.now()) + '] == Loading BED file', file=sys.stdout)
		print('## Loading BED file', file=sys.stderr)
		bed_regions = read_bed_sorted_list( options.bed )
		# bed_regions[line] = [chrom , chromStart(0->) , chromEnd(1->) , name , score , strand , ... ]
	else :
		bed_regions = ""

	print('[' + str(datetime.datetime.now()) + '] == Loading structure(s)', file=sys.stdout)
	print('## Loading structure(s)', file=sys.stderr)
	structure_file_format = str(options.format).lower()
	print('### Structure file(s) format: ' + structure_file_format, file=sys.stderr)

	if options.skipoverlap :
		# Convert input structures to AGP and export results
		if structure_file_format == "bed" :
			print("[MEMO] Structure in BED format. A new sequence will be generated for each input file. Progressive numeric id will follow input order.", file=sys.stdout)
			id = 0
			agp_db = {}
			for file_name in structure_files_list :
				id += 1
				bed_db = read_bed_sorted_list(file_name)
				print('### Input file id: ' + str(id) + " | Input BED File name: " + file_name + " | Corresponfing sequence IDs: " + options.prefix + "_" + str(id), file=sys.stderr)
				if agp_db == {} :
					agp_db = bed_to_agp_onefile(bed_db , gap_size , options.prefix + "_" + str(id))
				else :
					agp_db.update(bed_to_agp_onefile(bed_db , gap_size , options.prefix + "_" + str(id)))

		elif structure_file_format == "agp" :
			#### Read AGP
			agp_db = {}
			for file_name in structure_files_list :
				if agp_db == {} :
					agp_db = read_agp(file_name)
				else :
					agp_db.update(read_agp(file_name))

			mode = "old_to_new"
			if options.mode :
				print("[WARNING] Reverse use of structure AGP file requested [--reverse]", file=sys.stdout)
				print("[WARNING] FASTA file will not be produced", file=sys.stdout)
				no_fasta = True
				mode = "new_to_old"
				# Invert AGP
				new_agp_db = invert_agp(agp_db)
				agp_db = clean_self_agp(new_agp_db)
				#agp_file_name = options.out + ".inverted.agp"
				#agp_file_name = write_agp( new_agp_db , agp_file_name )

		elif structure_file_format == "block" :
			agp_db = {}
			for file_name in structure_files_list :
				block_db = read_block(file_name)
				if agp_db == {} :
					agp_db = block_to_agp(block_db , options.gap_size )
				else :
					agp_db.update( block_to_agp(block_db) )

		else :
			print("[ERROR] Structure file format " + str(options.format) + " unknown.", file=sys.stderr)
			print("[ERROR] Structure file format " + str(options.format) + " unknown.", file=sys.stdout)
			parser.print_help()
			sys.exit(1)

		if options.add_unplaced :
			print('### Adding unplaced sequences to output', file=sys.stderr)
			used_sequences = []
			# agp_db[seq_id][int(start)] = [ Obj_Name , Obj_start , Obj_End , PartNum , Compnt_Type , CompntId , CompntStart , CompntEnd ,  Orientation ]
			for seq_id in list(agp_db.keys()) :
				for start in list(agp_db[seq_id].keys()) :
					if agp_db[seq_id][int(start)][4] == "W" :
						used_sequences.append(agp_db[seq_id][int(start)][5])

			for seq_id in sorted(fasta_len_dict.keys()) :
				if seq_id not in used_sequences :
					seq_length = fasta_len_dict(seq_id)
					agp_db[seq_id][1] = [ seq_id , 1 , seq_length , 1 , "W" , seq_id , 1 , seq_length ,  "+" ]

	else :
		# Convert structures to BLOCK,
		# identify overlaps between blocks
		# correct coordinates
		# convert to AGP
		# export results
		if structure_file_format == "block" :
			block_db = {}
			for file_name in structure_files_list :
				print('### Loading ' + file_name, file=sys.stderr)
				if block_db == {} :
					block_db = read_block(file_name)
				else :
					block_db.update( read_block(file_name) )

		elif structure_file_format == "agp" :
			if options.mode :
				print("[ERROR] Reverse use of AGP structure file incompatible with overlap search and correction.", file=sys.stdout)
				print("[ERROR] Reverse use of AGP structure file incompatible with overlap search and correction.", file=sys.stdout)
				sys.exit(1)
			agp_db = {}
			for file_name in structure_files_list :
				print('### Loading ' + file_name, file=sys.stderr)
				if agp_db == {} :
					agp_db = read_agp(file_name)
				else :
					agp_db.update(read_agp(file_name))
			block_db = agp_to_block( agp_db )

		elif structure_file_format == "bed" :
			id = 0
			agp_db = {}
			for file_name in structure_files_list :
				print('### Loading ' + file_name, file=sys.stderr)
				id += 1
				bed_db = read_bed_sorted_list(file_name)
				print('### Input file id: ' + str(id) + " | Input BED File name: " + file_name + " | Corresponfing sequence IDs: " + options.prefix + "_" + str(id), file=sys.stderr)
				if agp_db == {} :
					agp_db = bed_to_agp_onefile(bed_db , gap_size , options.prefix + "_" + str(id))
				else :
					agp_db.update(bed_to_agp_onefile(bed_db , gap_size , options.prefix + "_" + str(id)))
			block_db = agp_to_block( agp_db )

		else :
			print("[ERROR] Structure file format " + str(options.format) + " unknown.", file=sys.stderr)
			print("[ERROR] Structure file format " + str(options.format) + " unknown.", file=sys.stdout)
			parser.print_help()
			sys.exit(1)

		# Correct agp regions with smart overlap dodging
		temp_dir = options.out + ".temp_dir"
		mkdir(temp_dir)
		print('[' + str(datetime.datetime.now()) + '] = Refining structure(s) based on overlap', file=sys.stdout)
		print('# Refining structure(s) based on overlap', file=sys.stderr)
		agp_db , harmed_loci = dodge_overlaps( block_db , fasta_dict , fasta_len_dict , int(options.gap_size) , int(options.spacer) , int(options.cores) , options.mapper , annotation_gff3 , temp_dir , paths, options.add_unplaced)
		if not harmed_loci == [] :
			harmed_loci_file = open( options.out + ".loci_to_check.txt" , "w+")
			for name in harmed_loci :
				print(name, file=harmed_loci_file)
			harmed_loci_file.close()

		if options.ignoreids :
			agp_db = rename_agp_sequences(agp_db , options.prefix )


	# Convert legacy agp if given
	if options.agp :
		print('[' + str(datetime.datetime.now()) + '] = Translating legacy components position from legacy agp file', file=sys.stdout)
		print('# Translating legacy components position from legacy agp file', file=sys.stderr)
		old_agp = read_agp(options.agp)
		new_to_legacy_agp_db = agp_translate_agp(agp_db , old_agp)

	# Convert BED file if given
	if options.bed :
		print('[' + str(datetime.datetime.now()) + '] = Translating coordinates of features in bed the file', file=sys.stdout)
		print('# Translating coordinates of features in bed the file', file=sys.stderr)
		new_bed = translate_bed_sorted_list( bed_regions ,  agp_db )

	# Writing output files
	print('[' + str(datetime.datetime.now()) + '] = Writing output files', file=sys.stdout)
	print('# Writing output files', file=sys.stderr)

	if options.skipoverlap :
		if not options.noagp :
			agp_file_name = options.out + ".structure.agp"
			agp_file_name = write_agp( agp_db , agp_file_name )
		export_from_agp(options.out, no_fasta, agp_db, fasta_dict, "old_to_new" , fasta_len_dict, annotation_gff3)

	else :
		export_from_agp(options.out, no_fasta, agp_db, fasta_dict, "old_to_new", fasta_len_dict, annotation_gff3)
		agp_file_name = options.out + ".structure.agp"
		agp_file_name = write_agp( agp_db , agp_file_name )

	if options.agp :
		legacy_agp_file_name = options.out + ".legacy_structure.agp"
		legacy_agp_file_name = write_agp( new_to_legacy_agp_db , legacy_agp_file_name )

	if options.bed :
		out_bed_file_name = options.out + ".bed"
		out_bed_file = open(out_bed_file_name , 'w')
		for line in sorted(new_bed) :
			print("\t".join([str(x) for x in line]), file=out_bed_file)
		out_bed_file.close()

	###### END ######

	print("------------------------------", file=sys.stdout)
	print("- Done", file=sys.stdout)
	print("------------------------------", file=sys.stdout)
	print("##############################", file=sys.stderr)
	print("# Done", file=sys.stderr)
	print("##############################", file=sys.stderr)


if __name__ == '__main__':
	main()
