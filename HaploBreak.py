#!/usr/bin/env python


import argparse
from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *
import sys

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


def main() :

	###### Options and help ######

	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta", dest="fasta",
					help="Genome sequences in FASTA format. Use unmasked or soft masked sequences only. [REQUIRED]", metavar="genome.fasta")
	parser.add_argument("-b", "--breakpoints", dest="breaks",
					help="List of breakpoints to be applied [REQUIRED]. Tab separated format, 1) sequence id; 2)comma separated values of regions to break: start1-stop1[,start2-stop2, ... ,startN-stopN]", metavar="breakpoints.tsv")
	parser.add_argument("-o", "--out", dest="out", default="out" ,
					help="Output file names [Default: out]", metavar="outname")
	parser.add_argument("-p", "--prefix", dest="prefix", default="" ,
					help="New sequence names prefix to be added before the old name. Leave unset to keep the same names [Default: NONE]", metavar="new")
	parser.add_argument("-g", "--gff3", dest="gff3",
					help="Genome annotation in GFF3 format.", metavar="annotation.gff3")
	parser.add_argument("-a", "--agp", dest="agp",
					help="AGP file defining the coordinates of legacy sequences composing the input assembly", metavar="previous_to_actual.genome.agp")
	parser.add_argument("--bed", dest="bed",
					help="Region annotation in BED format. Accepted BED3 and BED6 formats, additional columns will be considered as annotation and reported unedited", metavar="annotation.bed")

	parser.add_argument("-d", "--distance", dest="distance", default="25000" ,
					help="Maximum distance (in bp) allowed to search precise breakpoint into junctions and trigger in FASTA gap search [Default: 25,000]", metavar="25000")
	parser.add_argument("-m", "--minimum", dest="min", default="200" ,
					help="Minimum gap size (in bp) on fasta for break point selection [Default: 200]", metavar="200")
	parser.add_argument("--allow_annotation_breaking", dest="break_annot", default=False, action="store_true",
					help="Allow breaking annotated genes if GFF3 is reported")

	scriptDirectory = os.path.dirname(os.path.realpath(__file__)) + "/support_scripts"
	print("Running HaploBreaker tool from HaploSync version " + get_version(), file=sys.stdout)
	print("To reproduce this run use the following command: " + " ".join( pipes.quote(x) for x in sys.argv), file=sys.stdout)
	print("----", file=sys.stdout)
	# Sanity Check

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()

	if not options.fasta :
		# noinspection PyCompatibility
		print("[ERROR] Genome FASTA file missing", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	if not options.breaks :
		print("[ERROR] Breakpoints file missing", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	if not options.gff3 :
		print("[NOTE] Annotation file missing: Coordinates conversion will not take place (you must run it downstream) and breakpoints selection will not trigger gene loci disruption protection", file=sys.stderr)

	if options.agp :
		print("[NOTE] AGP file given: breaking on path gap will be prioritized. Any breaking of original sequences will be notified. Output will be delivered also for original sequences", file=sys.stderr)
	else :
		print("[WARNING] No AGP file was given, \"N\" encoded gaps in the sequence will be used for searching safe breakpoints only", file=sys.stderr)
	max_distance = int(options.distance)

	##### Reads inputs
	#### Read sequence
	print("- Importing genome sequences", file=sys.stdout)
	fasta_in = read_fasta(options.fasta)
	sequences = {}
	for id in fasta_in :
		sequences[id] = remove_islets_fasta(fasta_in[id].upper())

	print("-- Calculating reference sequences lengths", file=sys.stdout)
	sequences_len = get_length_from_fasta_db(sequences)

	print("--- Number of reference sequences: " + str(len(list(sequences_len.keys()))), file=sys.stdout)

	#### Read breaking points
	print("- Importing breaking points", file=sys.stdout)
	breakpoint_db = {}
	for line in open( options.breaks , "r") :
		try :
			seq_id , breaks = line.rstrip().rstrip(",").split("\t")
		except :
			continue
		if seq_id not in sequences :
			print("[[WARNING]] Breakpoints reported for sequence missing from FASTA file: [" + seq_id + "]. Data discarded", file=sys.stderr)
			continue
		if breaks == "" :
			# No breakpoints
			continue
		else :
			# Sanity check
			breaks_list = breaks.split(",")
			for break_couple in breaks_list :
				try :
					start , stop = break_couple.rstrip().lstrip().split("-")
				except :
					print("[ERROR] Error for breakpoint  coordinates in " + seq_id + " : " + str(break_couple.rstrip()) + " found, not a valid pair", file=sys.stderr)
					sys.exit(1)
				if start.upper() == "BEGIN" or start.upper() == "START" :
					start = 0
				else :
					try :
						start = int(start)
					except :
						print("[ERROR] Non integer coordinates for breakpoint in " + seq_id + " : " + str(start) +" found", file=sys.stderr)
						sys.exit(1)
				if stop.upper() == "END" or stop.upper() == "STOP" :
					stop = sequences_len(seq_id)
				else :
					try:
						stop = int(stop)
					except :
						print("[ERROR] Non integer coordinates for breakpoint in " + seq_id + " : " + str(stop) +" found", file=sys.stderr)
						sys.exit(1)
				if seq_id not in breakpoint_db :
					breakpoint_db[seq_id] = []
				breakpoint_db[seq_id].append( [ start , stop ] )

	print("-- Breaking points imported:", file=sys.stdout)
	for seq_id in sorted(sequences.keys()) :
		break_line = ""
		if seq_id in breakpoint_db :
			for element in sorted(breakpoint_db[seq_id]) :
				break_line += "(" + str(element[0]) + ":" + str(element[1]) + "); "
		if not break_line == "" :
			print("--- Sequence " + str(seq_id) + ": " + break_line, file=sys.stdout)

	#### - IF GIVEN - Read annotation
	if options.gff3:
		print("- Importing GFF3", file=sys.stdout)
		annotation_gff3 , mRNA_db = read_gff3(options.gff3)
		gene_db = feature_ranges(annotation_gff3 , "gene")
	else :
		annotation_gff3 = ""
		gene_db = {}
		for seq_id in sorted(sequences.keys()) :
			gene_db[seq_id] = []

	#### List gaps in FASTA sequences
	print("- Extracting position of gaps in sequences", file=sys.stdout)

	if options.gff3 and not options.break_annot :
		#gap_list = get_gap_from_fasta_db(sequences , {})
		#print >> sys.stderr , "gap list generated without gene masking"
		#print >> sys.stderr , gap_list
		gap_list = get_gap_from_fasta_db( sequences , gene_db )
		#print >> sys.stderr , "gap list generated with gene masking"
		#print >> sys.stderr , gap_list
	else :
		gap_list = get_gap_from_fasta_db(sequences , {})
		print("gap list generated", file=sys.stderr)

	gaps_file_name = options.out + ".breakable_gaps.list.txt"
	gaps_file = open(gaps_file_name , 'w')
	for seq in gap_list :
		for element in sorted(gap_list[seq]) :
			print("\t".join([ str(x) for x in element]), file=gaps_file)
	gaps_file.close()

	gapped_agp_db = get_agp_from_gap_list(gap_list , sequences )
	ungapped_db = agp_ungapped_to_gap(gapped_agp_db , gene_db , int(options.min) )

	#### - IF GIVEN - Read AGP
	if options.agp :
		print("- Importing AGP and original sequence information", file=sys.stdout)
		agp_db , original_sequence_position , junctions , sequences_legacy , original_gap_list = get_original_from_agp( options.agp , sequences )
		sequences_legacy_length_db = get_length_from_fasta_db(sequences_legacy)
		#### Test uniqueness of old ids
		print("-- Testing uniqueness AGP ids", file=sys.stdout)
		overlapping_original_ranges = test_range_uniqueness(agp_db)
		if overlapping_original_ranges == {} :
			print("--- None of the original ids was used more than once", file=sys.stdout)
			translation_db = translate_from_AGP_whole_genome_reverse(agp_db)
			if options.gff3:
				print("--- Converting annotation", file=sys.stdout)
				legacy_annotation_gff3 = translate_gff3( annotation_gff3 , translation_db , options.out + ".broken_loci.given_to_legacy_transfer.txt" )
				write_gff3(legacy_annotation_gff3 , options.out + ".legacy.annotation.gff3" , get_length_from_fasta_db( sequences_legacy ) )
			else :
				legacy_annotation_gff3 = {}
				for seq_id in sequences_legacy :
					legacy_annotation_gff3[seq_id] = []
			legacy_gap_list = get_gap_from_fasta_db( sequences_legacy , legacy_annotation_gff3 )
			legacy_gaps_file_name = options.out + ".legacy_breakable_gaps.list.txt"
			legacy_gaps_file = open(legacy_gaps_file_name , 'w')
			for seq in sorted(legacy_gap_list.keys()) :
				for element in sorted(legacy_gap_list[seq]) :
					print("\t".join([ str(x) for x in element]), file=legacy_gaps_file)
			legacy_gaps_file.close()

		else :
			print("--- " + str(len(list(overlapping_original_ranges.keys()))) + " original ids were used multiple times", file=sys.stdout)
			# TODO: Uniquify and run on uniquified
			print("--- Please check and correct", file=sys.stdout)
			print("------------------------------", file=sys.stdout)
			print("- Quitting", file=sys.stdout)
			print("------------------------------", file=sys.stdout)
			print_overlapping_info(overlapping_original_ranges, ".error.overlapping_ranges.info")
			sys.exit(2)
			#print >> sys.stdout, "---- Generating intermediate unique files"
			# TODO: uniquify_agp in AGP_lib.py
			# TODO: Save uniquified files
			# TODO: updated variables in script

	##### Find gap to break
	print("- Searching for best breaking position", file=sys.stdout)

	regions_to_break = {}
	regions_to_break_old = {}
	blacklist = {}
	#print sys.stderr , breakpoint_db.keys()
	for seq_id in sorted(breakpoint_db.keys()) :
		blacklist[seq_id] = []
		print("-- Analysing sequence: " + str(seq_id), file=sys.stdout)
		print("## Analysing sequence: " + str(seq_id), file=sys.stderr)
		for break_couple in sorted(breakpoint_db[seq_id]) :
			start , stop = break_couple
			print("--- Breakpoints: " + str(start) + ":" + str(stop), file=sys.stdout)
			print("### Breakpoints: " + str(start) + ":" + str(stop), file=sys.stderr)
			if options.agp :
				#### AGP IS GIVEN - search in AGP junctions for breaking point
				# If the region is between two junction -> blacklist the legacy sequence in-between
				print("---- Searching in given AGP for original sequence junction", file=sys.stdout)
				print("#### Searching in given AGP for original sequence junction", file=sys.stderr)

				start_gap = find_nearest_gap(start , junctions[seq_id] , sequences_len[seq_id] , "both" , max_distance )
				stop_gap = find_nearest_gap(stop , junctions[seq_id] , sequences_len[seq_id] , "both" , max_distance )
				# Test the junctions
				#print >> sys.stderr, start_gap
				#print >> sys.stderr, stop_gap

				if ( not start_gap[0] == "sequence" ) and (not stop_gap[0] == "sequence" ) :
					# The breakpoints hit a junction or the end of the sequence
					# Check it is not the same
					print("----- Breakpoints lead to known junctions", file=sys.stdout)
					print("##### Breakpoints lead to known junctions on " + seq_id + ": [" + str(start_gap[1]) + ":" + str(start_gap[2]) + "] -> [" + str(stop_gap[1]) + ":" + str(stop_gap[2]) + "]", file=sys.stderr)

					if ( start_gap == stop_gap ) or ( start_gap[0] == "extremity" and stop_gap[0] == "extremity" ) :
						# Same junction/extremity or from one extremity to the other:

						if ( start_gap == stop_gap ) :
							print("------ Right and left breakpoints lead to the same junction.", file=sys.stdout)
						else :
							print("------ Breakpoints lead to sequence extremities", file=sys.stdout)

						# assign the junction to the nearest, break the sequence for the other
						start_sq_distance = min( (int(start)- start_gap[1])**2 , (int(start) - start_gap[2])**2 )
						stop_sq_distance = min( (int(stop)- start_gap[1])**2 , (int(stop) - start_gap[2])**2 )

						if start_sq_distance <= stop_sq_distance :
							print("------ Searching a gap in sequence compatible with right breakpoint.", file=sys.stdout)
							legacy_stop_gap_seq_id, legacy_stop = translate_coords( stop , translation_db[seq_id] )
							print("------ Breakpoint position in legacy sequence: " + seq_id + ":" + str(stop) + "->" + str(legacy_stop_gap_seq_id) + ":" + str(legacy_stop), file=sys.stdout)
							print("###### Breakpoint position in legacy sequence: " + seq_id + ":" + str(stop) + "->" + str(legacy_stop_gap_seq_id) + ":" + str(legacy_stop), file=sys.stderr)
							legacy_stop_gap = find_nearest_gap( legacy_stop , legacy_gap_list[legacy_stop_gap_seq_id] , sequences_legacy_length_db[legacy_stop_gap_seq_id] , "both" , sequences_legacy_length_db[legacy_stop_gap_seq_id] )

						else :
							print("------ Searching a gap in sequence compatible with left breakpoint.", file=sys.stdout)
							legacy_start_gap_seq_id , legacy_start = translate_coords( start , translation_db[seq_id] )
							print("------ Breakpoint position in legacy sequence: " + seq_id + ":" + str(start) + "->" + str(legacy_start_gap_seq_id) + ":" + str(legacy_start), file=sys.stdout)
							print("###### Breakpoint position in legacy sequence: " + seq_id + ":" + str(start) + "->" + str(legacy_start_gap_seq_id) + ":" + str(legacy_start), file=sys.stderr)
							legacy_start_gap = find_nearest_gap( legacy_start , legacy_gap_list[legacy_start_gap_seq_id] , sequences_legacy_length_db[legacy_start_gap_seq_id] , "both" , sequences_legacy_length_db[legacy_start_gap_seq_id] )

					else :
						# Different junctions/extremity
						# Blacklist the sequences in between
						print("------ Junctions to break found: [" + str(start_gap[1]) + ":" + str(start_gap[2]) + "] <-> [" + str(stop_gap[1]) + ":" + str(stop_gap[2]) + "]", file=sys.stdout)
						print("###### Junctions to break found: [" + str(start_gap[1]) + ":" + str(start_gap[2]) + "] <-> [" + str(stop_gap[1]) + ":" + str(stop_gap[2]) + "]", file=sys.stderr)
						for component_start in agp_db[seq_id] :
							if (start_gap[1] <= int(component_start) <= stop_gap[2]) and (agp_db[seq_id][component_start][4] == "W") :
								# if the component in the AGP is a sequence and falls withing the breakpoints -> blacklist
								blacklist[seq_id].append(agp_db[seq_id][component_start][5])
						print("###### Blacklisted for " + str(seq_id) + " : " + ",".join([str(x) for x in blacklist[seq_id] ]), file=sys.stderr)
						print("------ Blacklisted " + str(len(blacklist[seq_id])) + " sequences", file=sys.stdout)

				else :
					# At least one breakpoint falls inside a sequence.
					# Break the original in the nearest gap
					legacy_start_gap = False
					legacy_stop_gap = False

					if start_gap[0] == "sequence" :
						# No junction is compatible with the breakpoint
						# Search a breakpoint in sequence gaps
						print("----- No junction compatible with left breakpoint. Searching suitable gap in sequence", file=sys.stdout)
						print("##### No junction compatible with left breakpoint. Searching suitable gap in sequence", file=sys.stderr)
						legacy_start_gap_seq_id , legacy_start = translate_coords( start , translation_db[seq_id] )
						print("------ Breakpoint position in legacy sequence: " + seq_id + ":" + str(start) + "->" + str(legacy_start_gap_seq_id) + ":" + str(legacy_start), file=sys.stdout)
						print("###### Breakpoint position in legacy sequence: " + seq_id + ":" + str(start) + "->" + str(legacy_start_gap_seq_id) + ":" + str(legacy_start), file=sys.stderr)
						legacy_start_gap = find_nearest_gap( legacy_start , legacy_gap_list[legacy_start_gap_seq_id] , sequences_legacy_length_db[legacy_start_gap_seq_id] , "both" , sequences_legacy_length_db[legacy_start_gap_seq_id] )
					else :
						# junction is compatible with the breakpoint
						print("----- Junction compatible with left breakpoint found: " + seq_id + ":" + str(start_gap[1]) + "-" + str(start_gap[2]), file=sys.stdout)
						print("##### Junction compatible with left breakpoint found: " + seq_id + ":" + str(start_gap[1]) + "-" + str(start_gap[2]), file=sys.stderr)

					if stop_gap[0] == "sequence" :
						# No junction is compatible with the breakpoint
						# Search a breakpoint in sequence gaps
						print("----- No junction compatible with right breakpoint. Searching suitable gap in sequence", file=sys.stdout)
						print("##### No junction compatible with right breakpoint. Searching suitable gap in sequence", file=sys.stderr)
						legacy_stop_gap_seq_id, legacy_stop = translate_coords( stop , translation_db[seq_id] )
						print("------ Breakpoint position in legacy sequence: " + seq_id + ":" + str(stop) + "->" + str(legacy_stop_gap_seq_id) + ":" + str(legacy_stop), file=sys.stdout)
						print("###### Breakpoint position in legacy sequence: " + seq_id + ":" + str(stop) + "->" + str(legacy_stop_gap_seq_id) + ":" + str(legacy_stop), file=sys.stderr)
						legacy_stop_gap = find_nearest_gap( legacy_stop , legacy_gap_list[legacy_stop_gap_seq_id] , sequences_legacy_length_db[legacy_stop_gap_seq_id] , "both" , sequences_legacy_length_db[legacy_stop_gap_seq_id] )
					else :
						# junction is compatible with the breakpoint
						print("----- Junction compatible with right breakpoint found: " + seq_id + ":" + str(stop_gap[1]) + "-" + str(stop_gap[2]), file=sys.stdout)
						print("##### Junction compatible with right breakpoint found: " + seq_id + ":" + str(stop_gap[1]) + "-" + str(stop_gap[2]), file=sys.stderr)

					if legacy_start_gap and legacy_stop_gap and (legacy_start_gap == legacy_stop_gap) :
						print("----- [WARNING] : The same gap is the nearest to both breakpoints and the same gap. Manual inspection suggested", file=sys.stdout)
						print("[WARNING] : Breakpoints - " + str(start) + ":" + str(stop) + " leading to the same gap", file=sys.stderr)
						print("[WARNING] : Nearest gap - " + str(legacy_stop_gap[1]) + ":" + str(legacy_stop_gap[2]), file=sys.stderr)
						print("[WARNING] : Junctions in AGP:", file=sys.stderr)
						#for junction in junctions[seq_id] :
						#	print >> sys.stderr, junction
						#print >> sys.stderr, "[WARNING] : Usable gaps in sequence:"
						#for gap in gap_list[seq_id] :
						#	print >> sys.stderr, gap
						print("------ Forcing research to upstream of left and downstream of right", file=sys.stdout)
						print("###### Forcing research to upstream of left and downstream of right", file=sys.stderr)
						nearest = nearest_to_gap( legacy_start , legacy_stop , legacy_start_gap )
						if legacy_start < legacy_stop :
							if nearest == legacy_stop :
								legacy_start_gap = find_nearest_gap( legacy_start , legacy_gap_list[legacy_start_gap_seq_id] , sequences_legacy_length_db[legacy_start_gap_seq_id] , "upstream" , sequences_legacy_length_db[legacy_start_gap_seq_id] )
								print("------- Forced search: upstream gap of left breakpoint:" + str(legacy_start_gap), file=sys.stdout)
								print("------- Forced search: distance of left breakpoint from upstream gap:" + str(int(legacy_start) - int(legacy_start_gap[2])), file=sys.stdout)
								print("####### Forced search: upstream gap of left breakpoint:" + str(legacy_start_gap), file=sys.stderr)
								print("####### Forced search: distance of left breakpoint from upstream gap:" + str(int(legacy_start) - int(legacy_start_gap[2])), file=sys.stderr)
							else :
								legacy_stop_gap = find_nearest_gap( legacy_stop , legacy_gap_list[legacy_stop_gap_seq_id] , sequences_legacy_length_db[legacy_stop_gap_seq_id] , "downstream" , sequences_legacy_length_db[legacy_stop_gap_seq_id] )
								print("------- Forced search: downstream gap of left breakpoint:" + str(legacy_stop_gap), file=sys.stdout)
								print("------- Forced search: distance of right breakpoint from downstream gap" + str(int(legacy_stop_gap[1]) - int(legacy_stop)), file=sys.stdout)
								print("####### Forced search: downstream gap of left breakpoint:" + str(legacy_stop_gap), file=sys.stderr)
								print("####### Forced search: distance of right breakpoint from downstream gap" + str(int(legacy_stop_gap[1]) - int(legacy_stop)), file=sys.stderr)
						else:
							# Legacy sequence is in reverse direction
							if nearest == legacy_stop :
								legacy_start_gap = find_nearest_gap( legacy_start , legacy_gap_list[legacy_start_gap_seq_id] , sequences_legacy_length_db[legacy_start_gap_seq_id] , "downstream" , sequences_legacy_length_db[legacy_start_gap_seq_id] )
								print("------- Forced search: upstream gap of left breakpoint:" + str(legacy_start_gap), file=sys.stdout)
								print("------- Forced search: distance of left breakpoint from upstream gap:" + str(int(legacy_start) - int(legacy_start_gap[2])), file=sys.stdout)
								print("####### Forced search: upstream gap of left breakpoint:" + str(legacy_start_gap), file=sys.stderr)
								print("####### Forced search: distance of left breakpoint from upstream gap:" + str(int(legacy_start) - int(legacy_start_gap[2])), file=sys.stderr)
							else :
								legacy_stop_gap = find_nearest_gap( legacy_stop , legacy_gap_list[legacy_stop_gap_seq_id] , sequences_legacy_length_db[legacy_stop_gap_seq_id] , "upstream" , sequences_legacy_length_db[legacy_stop_gap_seq_id] )
								print("------- Forced search: downstream gap of left breakpoint:" + str(legacy_stop_gap), file=sys.stdout)
								print("------- Forced search: distance of right breakpoint from downstream gap" + str(int(legacy_stop_gap[1]) - int(legacy_stop)), file=sys.stdout)
								print("####### Forced search: downstream gap of left breakpoint:" + str(legacy_stop_gap), file=sys.stderr)
								print("####### Forced search: distance of right breakpoint from downstream gap" + str(int(legacy_stop_gap[1]) - int(legacy_stop)), file=sys.stderr)

					if legacy_start_gap :
						if legacy_start_gap_seq_id not in regions_to_break_old :
							regions_to_break_old[legacy_start_gap_seq_id] = []
						regions_to_break_old[legacy_start_gap_seq_id].append( legacy_start_gap )
					if legacy_stop_gap :
						if legacy_stop_gap_seq_id not in regions_to_break_old :
							regions_to_break_old[legacy_stop_gap_seq_id] = []
						regions_to_break_old[legacy_stop_gap_seq_id].append( legacy_stop_gap )

				if seq_id not in regions_to_break :
					regions_to_break[seq_id] = []
				regions_to_break[seq_id].append(start_gap)
				regions_to_break[seq_id].append(stop_gap)

			else :
				#### AGP IS NOT GIVEN
				## Only in-sequence-gap search
				print("---- No AGP given, searching in the sequence encoded gaps", file=sys.stdout)
				stop_gap_new = find_nearest_gap( stop ,  gap_list[seq_id] , sequences_len[seq_id] , "both" , sequences_len[seq_id] )
				start_gap_new = find_nearest_gap( start , gap_list[seq_id] , sequences_len[seq_id] , "both" , sequences_len[seq_id] )
				if seq_id not in regions_to_break :
					regions_to_break[seq_id] = []
				regions_to_break[seq_id].append( start_gap_new )
				regions_to_break[seq_id].append( stop_gap_new )

	broken_gaps_file_name = options.out + ".broken_gaps.given_sequences.list.txt"
	broken_gaps_file = open(broken_gaps_file_name , 'w')
	print("# Final list of gap break", file=sys.stderr)
	print("## Given sequences", file=sys.stderr)
	for seq_id in sorted(regions_to_break.keys()) :
		# uniquify the list of gap to break
		old_list = regions_to_break[seq_id]
		new_list = list(set(tuple(i) for i in old_list))
		regions_to_break[seq_id] = new_list
		print('>' + seq_id, file=broken_gaps_file)
		for element in regions_to_break[seq_id] :
			print("\t".join([ str(x) for x in element ]), file=broken_gaps_file)
	#print >> sys.stderr, regions_to_break
	broken_gaps_file.close()

	if not regions_to_break_old == {} :
		print("## Legacy sequences", file=sys.stderr)
		broken_gaps_legacy_file_name = options.out + ".broken_gaps.legacy_sequences.list.txt"
		broken_gaps_legacy_file = open(broken_gaps_legacy_file_name , 'w')
		for legacy_stop_gap_seq_id in sorted(regions_to_break_old.keys()) :
			# uniquify the list of gap to break
			old_list = regions_to_break_old[legacy_stop_gap_seq_id]
			new_list = list(set(tuple(i) for i in old_list))
			regions_to_break_old[legacy_stop_gap_seq_id] = new_list
			print(">" + legacy_stop_gap_seq_id, file=broken_gaps_legacy_file)
			for element in regions_to_break_old[legacy_stop_gap_seq_id] :
				print("\t".join([ str(x) for x in element ]), file=broken_gaps_legacy_file)
		#print >> sys.stdout, regions_to_break_old
		broken_gaps_legacy_file.close()


	##### Generate new AGP using breakpoints
	print("- Updating sequences with novel breaking points", file=sys.stdout)
	print("# Updating sequences with novel breaking points", file=sys.stderr)

	if options.agp :
		print("-- Splitting legacy sequences by gap breaking", file=sys.stdout)
		new_legacy_sequences , new_legacy_agp_db = get_broken_agp( sequences_legacy, regions_to_break_old , options.prefix)
		new_legacy_agp_file = write_agp( new_legacy_agp_db , options.out + ".new_legacy.agp" )
		#print >> sys.stderr, new_legacy_agp_db
		#legacy_to_new_legacy_agp_db = invert_agp( new_legacy_agp_db )
		#legacy_to_new_legacy_agp_file = write_agp( legacy_to_new_legacy_agp_db , options.out + ".legacy_to_new_legacy.agp" )
		#print >> sys.stdout, "-- Updating sequences to novel legacy sequences relationship"
		#legacy_to_new_agp_db = agp_translate_agp( agp_db , legacy_to_new_legacy_agp_db )
		#legacy_to_new_agp_file = write_agp( legacy_to_new_agp_db , options.out + ".legacy_to_new.agp" )
	else :
		print("-- Splitting sequences by gap breaking", file=sys.stdout)
		#print >> sys.stderr, regions_to_break
		new_sequences , new_agp_db = get_broken_agp(sequences, regions_to_break , options.prefix)
		new_agp_file = write_agp( new_agp_db , options.out + ".new.agp" )

	#### Print new sequences
	print("-- Printing new sequences", file=sys.stdout)
	if options.agp :
		write_fasta_from_db( new_legacy_sequences , options.out + ".new_legacy.fasta" )
	else :
		write_fasta_from_db( new_sequences , options.out + ".new.fasta" )

	#### Print blacklist
	print("-- Printing file of blacklisted sequences", file=sys.stdout)
	blacklist_file_name = options.out + ".blacklist.txt"
	blacklist_file = open(blacklist_file_name , 'w')
	for seq_id in sorted(blacklist.keys()) :
		info = ",".join(str(x) for x in blacklist[seq_id])
		if not info.rstrip().lstrip() == "" :
			blacklist_line = seq_id + "\t" + info
			print(blacklist_line, file=blacklist_file)
	blacklist_file.close()

	#### update and print new annotation files (if set)
	if options.gff3 :
		print("-- Updating annotation cooridnates", file=sys.stdout)
		print("--- Converting coordinates to updated sequences", file=sys.stdout)
		
		if options.agp :
			print("--- Converting coordinates to updated legacy sequences", file=sys.stdout)
			new_legacy_translate_db = translate_from_AGP_whole_genome(new_legacy_agp_db)
			#print >> sys.stderr, "##"
			#print >> sys.stderr, new_legacy_agp_db
			#print >> sys.stderr, "##"
			#print >> sys.stderr, legacy_annotation_gff3
			new_legacy_annotation_gff3 = translate_gff3( legacy_annotation_gff3 , new_legacy_translate_db , options.out + ".broken_loci.annotation.legacy_updated.txt" )
			print("---- Generating legacy GFF3 file for updated legacy sequences", file=sys.stdout)
			write_gff3(new_legacy_annotation_gff3 , options.out + ".new_legacy.annotation.gff3" , get_length_from_fasta_db( new_legacy_sequences ) )
		else:
			new_translate_db = translate_from_AGP_whole_genome(new_agp_db)
			new_annotation_gff3 = translate_gff3( annotation_gff3 , new_translate_db , options.out + ".broken_loci.annotation.txt" )
			print("---- Generating legacy GFF3 file for updated sequences", file=sys.stdout)
			write_gff3(new_annotation_gff3 , options.out + ".new.annotation.gff3" , get_length_from_fasta_db( new_sequences ) )

	# Convert BED file if given
	if options.bed :
		bed_regions = read_bed_sorted_list( options.bed )
		print('[' + str(datetime.datetime.now()) + '] = Translating coordinates of features in bed the file', file=sys.stdout)
		print('# Translating coordinates of features in bed the file', file=sys.stderr)
		if options.agp:
			new_legacy_bed = translate_bed_sorted_list( bed_regions ,  new_legacy_agp_db )
			out_bed_file_name = options.out + ".new_legacy.bed"
			out_bed_file = open(out_bed_file_name , 'w')
			for line in sorted(new_legacy_bed) :
				print("\t".join([str(x) for x in line]), file=out_bed_file)
			out_bed_file.close()
		else:
			new_bed = translate_bed_sorted_list( bed_regions ,  new_agp_db )
			out_bed_file_name = options.out + ".new.bed"
			out_bed_file = open(out_bed_file_name , 'w')
			for line in sorted(new_bed) :
				print("\t".join([str(x) for x in line]), file=out_bed_file)
			out_bed_file.close()

	##### Finished

	print("------------------------------", file=sys.stdout)
	print("- Done", file=sys.stdout)
	print("------------------------------", file=sys.stdout)
	print("##############################", file=sys.stderr)
	print("# Done", file=sys.stderr)
	print("##############################", file=sys.stderr)


if __name__ == '__main__':
	main()
