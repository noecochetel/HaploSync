#!/usr/bin/env python

import argparse
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *

gc.garbage.append(sys.stdout)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


def main():
	#### Main

	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--fasta", dest="fasta", nargs='+',
					help="FASTA file(s) of sequences to analyse. [Required]", metavar="seq.fasta")

	parser.add_argument("-g", "--genomesize" , dest="genome_size",
					help="Expected genome size. If given, it is used to calculate NG values", metavar="N")
	parser.add_argument("-s", "--nstep" , dest="nstep", default="10",
					help="Step between N values reported [Default: 10]", metavar="N")
	parser.add_argument("--ngstep" , dest="ngstep", default="1",
					help="Step between NG values reported [Default: 1]", metavar="N")
	parser.add_argument("-c", "--contigs", dest="contigs_stats", default=True, action="store_false",
					help="Avoid performing contigs identification and statistics [Default: perform]")

	parser.add_argument("--nolengths", dest="print_len", default=True, action="store_false",
					help="Avoid printing sequence length file [Default: print]")
	parser.add_argument("--printcontigs", dest="print_contigs", default=False, action="store_true",
					help="Print contigs sequences in [input_fasta].contigs.fasta [Default: do not print]")

	parser.add_argument("--gap", dest="gap_size", default="10",
					help="Minimum size of a Ns stretch to be considered a gap between contigs [Default: 10bp]", metavar="N")

	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit(1)

	options = parser.parse_args()
	gap_size = int(options.gap_size)

	if not options.fasta :
		print("[ERROR] Genome FASTA file missing", file=sys.stderr)
		parser.print_help()
		sys.exit(1)

	fasta_files_list = options.fasta

	Nsize_list = list(range( 0 , 100 , int(options.nstep)))
	NGsize_list = list(range( 0 , 100 , int(options.ngstep)))

	results_db = {}
	for sequence_file in fasta_files_list :
		print("# Processing " + sequence_file, file=sys.stderr)
		print("## Importing sequences", file=sys.stderr)
		fasta_in = read_fasta( sequence_file )
		reference = {}
		for id in fasta_in :
			reference[id] = remove_islets_fasta(fasta_in[id].upper())
		results_db[sequence_file] = {}

		print("## Analysing base content", file=sys.stderr)
		results_db[sequence_file]["base_content"] = gc_count( reference )

		print("## Analysing input sequences", file=sys.stderr)
		genome_seq_lengths = sorted([ len(reference[x]) for x in reference ] ,reverse=True)
		if options.print_len :
			print("## Printing sequences lengths in: " + sequence_file + ".len", file=sys.stderr)
			try :
				out_file = open( sequence_file + ".len" ,'w')
				for seq_id in sorted(reference.keys()) :
					print(seq_id + "\t" + str(len(reference[seq_id])), file=out_file)
				out_file.close()
			except IOError as e :
				print("[ERROR] Impossible to write on " + sequence_file + ".len" + "(Error " + str(e[0]) + ": " + e[1] + ")", file=sys.stderr)

		results_db[sequence_file]["len_stats"] = len_stats(genome_seq_lengths)

		results_db[sequence_file]["Nvalue_Length"] = {}
		results_db[sequence_file]["Nvalue_Index"] = {}
		for threshold in Nsize_list :
			N , L = Nvalue( int(threshold) , genome_seq_lengths)
			results_db[sequence_file]["Nvalue_Length"][threshold] = str('%.0f' % N)
			results_db[sequence_file]["Nvalue_Index"][threshold] = str('%.0f' % L)

		# If there the expected genome size is given, calculate also NG values
		if options.genome_size :
			results_db[sequence_file]["NGvalue_Length"] = {}
			results_db[sequence_file]["NGvalue_Index"] = {}
			for threshold in NGsize_list :
				N , L = NGvalue( int(threshold) , genome_seq_lengths , int(options.genome_size) )
				results_db[sequence_file]["NGvalue_Length"][threshold] = str('%.0f' % N)
				results_db[sequence_file]["NGvalue_Index"][threshold] = str('%.0f' % L)

		# Get gaps
		print("## Analysing gaps in input sequences", file=sys.stderr)
		gaps_in = get_gap_from_fasta_db(reference)
		# format:
		# gaps_in == gap_db[seqid] = [ [ id , start , stop ] , .. , [ ] ]
		results_db[sequence_file]["Gaps_stats"] = gap_stats( gaps_in , sum(genome_seq_lengths) , int(options.gap_size) )

		if options.contigs_stats :
			print("## Extracting contigs", file=sys.stderr)
			contigs_in = split_fasta_on_gap( reference , gaps_in , int(options.gap_size) )
			print("## Analysing contigs sequences", file=sys.stderr)
			contigs_seq_lengths = sorted([ len(contigs_in[x]) for x in contigs_in ] ,reverse=True)
			results_db[sequence_file]["contig_len_stats"] = len_stats(contigs_seq_lengths)

			results_db[sequence_file]["contig_Nvalue_Length"] = {}
			results_db[sequence_file]["contig_Nvalue_Index"] = {}
			for threshold in Nsize_list :
				N , L = Nvalue( int(threshold) , contigs_seq_lengths)
				results_db[sequence_file]["contig_Nvalue_Length"][threshold] = str('%.0f' % N)
				results_db[sequence_file]["contig_Nvalue_Index"][threshold] = str('%.0f' % L)

			# If there the expected genome size is given, calculate also NG values
			if options.genome_size :
				results_db[sequence_file]["contig_NGvalue_Length"] = {}
				results_db[sequence_file]["contig_NGvalue_Index"] = {}
				for threshold in NGsize_list :
					N , L = NGvalue( int(threshold) , contigs_seq_lengths , int(options.genome_size) )
					results_db[sequence_file]["contig_NGvalue_Length"][threshold] = str('%.0f' % N)
					results_db[sequence_file]["contig_NGvalue_Index"][threshold] = str('%.0f' % L)

			# Print contigs
			if options.print_contigs :
				print("## Printing contigs sequences in: " + sequence_file + ".contigs.fasta", file=sys.stderr)
				try :
					write_fasta_from_db(contigs_in , sequence_file + ".contigs.fasta" , False)
				except IOError as e :
					print("[ERROR] Impossible to write on " + sequence_file + ".contigs.fasta" + "(Error " + str(e[0]) + ": " + e[1] + ")", file=sys.stderr)

	print("# Printing statistics", file=sys.stderr)

	# print ids
	id_list = fasta_files_list
	print("# Input sequences statistics", file=sys.stdout)
	print("File name" + "\t" + "\t".join([ str(x) for x in id_list ] ), file=sys.stdout)

	print("Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Cumulative length" ) ), file=sys.stdout)
	print("Number of sequences" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Number of sequences" ) ), file=sys.stdout)
	print("Average sequence length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Average sequence length" ) ), file=sys.stdout)
	print("Median sequence length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Median sequence length" ) ), file=sys.stdout)
	print("Minimum sequence length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Minimum sequence length" ) ), file=sys.stdout)
	print("Maximum sequence length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Maximum sequence length" ) ), file=sys.stdout)

	print("# Sequence count by sequence size", file=sys.stdout)
	print("Sequences > 100b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100b - Count" ) ), file=sys.stdout)
	print("Sequences > 500b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500b - Count" ) ), file=sys.stdout)
	print("Sequences > 1Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Kb - Count" ) ), file=sys.stdout)
	print("Sequences > 5Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Kb - Count" ) ), file=sys.stdout)
	print("Sequences > 10Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 10Kb - Count" ) ), file=sys.stdout)
	print("Sequences > 50Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 50Kb - Count" ) ), file=sys.stdout)
	print("Sequences > 100Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100Kb - Count" ) ), file=sys.stdout)
	print("Sequences > 500Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500Kb - Count" ) ), file=sys.stdout)
	print("Sequences > 1Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Mb - Count" ) ), file=sys.stdout)
	print("Sequences > 5Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Mb - Count" ) ), file=sys.stdout)

	print("# Sequence cumulative length by sequence size", file=sys.stdout)
	print("Sequences > 100b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100b - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 500b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500b - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 1Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Kb - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 5Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Kb - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 10Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 10Kb - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 50Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 50Kb - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 100Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100Kb - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 500Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500Kb - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 1Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Mb - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 5Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Mb - Cumulative length" ) ), file=sys.stdout)

	print("# Sequence count fraction by sequence size", file=sys.stdout)
	print("Sequences > 100b ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100b ratio - Count" ) ), file=sys.stdout)
	print("Sequences > 500b ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500b ratio - Count" ) ), file=sys.stdout)
	print("Sequences > 1Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Kb ratio - Count" ) ), file=sys.stdout)
	print("Sequences > 5Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Kb ratio - Count" ) ), file=sys.stdout)
	print("Sequences > 10Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 10Kb ratio - Count" ) ), file=sys.stdout)
	print("Sequences > 50Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 50Kb ratio - Count" ) ), file=sys.stdout)
	print("Sequences > 100Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100Kb ratio - Count" ) ), file=sys.stdout)
	print("Sequences > 500Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500Kb ratio - Count" ) ), file=sys.stdout)
	print("Sequences > 1Mb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Mb ratio - Count" ) ), file=sys.stdout)
	print("Sequences > 5Mb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Mb ratio - Count" ) ), file=sys.stdout)

	print("# Sequence cumulative length fraction by sequence size", file=sys.stdout)
	print("Sequences > 100b ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100b ratio - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 500b ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500b ratio - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 1Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Kb ratio - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 5Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Kb ratio - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 10Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 10Kb ratio - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 50Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 50Kb ratio - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 100Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 100Kb ratio - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 500Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 500Kb ratio - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 1Mb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 1Mb ratio - Cumulative length" ) ), file=sys.stdout)
	print("Sequences > 5Mb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "len_stats" , "Sequences > 5Mb ratio - Cumulative length" ) ), file=sys.stdout)

	# N values, length and index
	print("# Sequence N-values length", file=sys.stdout)
	for value in Nsize_list :
		# Print N length
		print("N" + str(value) + " Length\t" + "\t".join( list_from_stats_db( results_db , id_list , "Nvalue_Length" , value ) ), file=sys.stdout)
	print("# Sequence N-values index", file=sys.stdout)
	for value in Nsize_list :
		# Print N index
		print("N" + str(value) + " Index\t" + "\t".join( list_from_stats_db( results_db , id_list , "Nvalue_Index" , value ) ), file=sys.stdout)
		
	# NG values, length and index
	if options.genome_size :
		print("# Sequence NG-values length (based on and expected genome size of " + options.genome_size + "bp)", file=sys.stdout)
		for value in NGsize_list :
			# Print NG length
			print("NG" + str(value) + " Length\t" + "\t".join( list_from_stats_db( results_db , id_list , "NGvalue_Length" , value ) ), file=sys.stdout)
		print("# Sequence NG-values index (based on and expected genome size of " + options.genome_size + "bp)", file=sys.stdout)
		for value in NGsize_list :
			# Print NG index
			print("NG" + str(value) + " Index\t" + "\t".join( list_from_stats_db( results_db , id_list , "NGvalue_Index" , value ) ), file=sys.stdout)

	# Composition
	print("# Sequences nucleotide content", file=sys.stdout)
	print("A count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "A_count" ) ), file=sys.stdout)
	print("C count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "C_count" ) ), file=sys.stdout)
	print("T count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "T_count" ) ), file=sys.stdout)
	print("G count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "G_count" ) ), file=sys.stdout)
	print("U count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "U_count" ) ), file=sys.stdout)
	print("N count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "N_count" ) ), file=sys.stdout)
	print("Other bases count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "Other_count" ) ), file=sys.stdout)
	print("GC percentage" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "GC_perc" ) ), file=sys.stdout)
	print("GC percentage without unknown bases" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "base_content" , "GC_perc_corr" ) ), file=sys.stdout)

	print("##############################", file=sys.stdout)
	print("# Gaps statistics", file=sys.stdout)
	print("File name" + "\t" + "\t".join([ str(x) for x in id_list ] ), file=sys.stdout)
	# Gaps
	print("Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Cumulative length" ) ), file=sys.stdout)
	print("Percentage of assembly length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Genome percentage" ) ), file=sys.stdout)
	print("Number of gaps" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Number of sequences" ) ), file=sys.stdout)
	print("Average gap length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Average sequence length" ) ), file=sys.stdout)
	print("Median gap length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Median sequence length" ) ), file=sys.stdout)
	print("Minimum gap length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Minimum sequence length" ) ), file=sys.stdout)
	print("Maximum gap length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Maximum sequence length" ) ), file=sys.stdout)

	print("# Gap count by sequence size", file=sys.stdout)
	print("Gaps > 100b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 100b - Count" ) ), file=sys.stdout)
	print("Gaps > 500b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 500b - Count" ) ), file=sys.stdout)
	print("Gaps > 1Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 1Kb - Count" ) ), file=sys.stdout)
	print("Gaps > 5Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 5Kb - Count" ) ), file=sys.stdout)
	print("Gaps > 10Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 10Kb - Count" ) ), file=sys.stdout)
	print("Gaps > 50Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 50Kb - Count" ) ), file=sys.stdout)
	print("Gaps > 100Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 100Kb - Count" ) ), file=sys.stdout)
	print("Gaps > 500Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 500Kb - Count" ) ), file=sys.stdout)
	print("Gaps > 1Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 1Mb - Count" ) ), file=sys.stdout)
	print("Gaps > 5Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 5Mb - Count" ) ), file=sys.stdout)

	print("# Gap cumulative length by sequence size", file=sys.stdout)
	print("Gaps > 100b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 100b - Cumulative length" ) ), file=sys.stdout)
	print("Gaps > 500b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 500b - Cumulative length" ) ), file=sys.stdout)
	print("Gaps > 1Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 1Kb - Cumulative length" ) ), file=sys.stdout)
	print("Gaps > 5Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 5Kb - Cumulative length" ) ), file=sys.stdout)
	print("Gaps > 10Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 10Kb - Cumulative length" ) ), file=sys.stdout)
	print("Gaps > 50Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 50Kb - Cumulative length" ) ), file=sys.stdout)
	print("Gaps > 100Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 100Kb - Cumulative length" ) ), file=sys.stdout)
	print("Gaps > 500Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 500Kb - Cumulative length" ) ), file=sys.stdout)
	print("Gaps > 1Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 1Mb - Cumulative length" ) ), file=sys.stdout)
	print("Gaps > 5Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "Gaps_stats" , "Sequences > 5Mb - Cumulative length" ) ), file=sys.stdout)

	# Contigs
	if options.contigs_stats :
		print("##############################", file=sys.stdout)
		print("# Contigs statistics", file=sys.stdout)
		print("File name" + "\t" + "\t".join([ str(x) for x in id_list ] ), file=sys.stdout)

		print("Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Cumulative length" ) ), file=sys.stdout)
		print("Number of contigs" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Number of sequences" ) ), file=sys.stdout)
		print("Average contig length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Average sequence length" ) ), file=sys.stdout)
		print("Median contig length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Median sequence length" ) ), file=sys.stdout)
		print("Minimum contig length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Minimum sequence length" ) ), file=sys.stdout)
		print("Maximum contig length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Maximum sequence length" ) ), file=sys.stdout)

		print("# Contig count by sequence size", file=sys.stdout)
		print("Contigs > 100b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100b - Count" ) ), file=sys.stdout)
		print("Contigs > 500b - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500b - Count" ) ), file=sys.stdout)
		print("Contigs > 1Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Kb - Count" ) ), file=sys.stdout)
		print("Contigs > 5Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Kb - Count" ) ), file=sys.stdout)
		print("Contigs > 10Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 10Kb - Count" ) ), file=sys.stdout)
		print("Contigs > 50Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 50Kb - Count" ) ), file=sys.stdout)
		print("Contigs > 100Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100Kb - Count" ) ), file=sys.stdout)
		print("Contigs > 500Kb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500Kb - Count" ) ), file=sys.stdout)
		print("Contigs > 1Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Mb - Count" ) ), file=sys.stdout)
		print("Contigs > 5Mb - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Mb - Count" ) ), file=sys.stdout)

		print("# Contig cumulative length by sequence size", file=sys.stdout)
		print("Contigs > 100b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100b - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 500b - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500b - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 1Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Kb - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 5Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Kb - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 10Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 10Kb - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 50Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 50Kb - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 100Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100Kb - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 500Kb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500Kb - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 1Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Mb - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 5Mb - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Mb - Cumulative length" ) ), file=sys.stdout)

		print("# Contig count ratio by sequence size", file=sys.stdout)
		print("Contigs > 100b ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100b ratio - Count" ) ), file=sys.stdout)
		print("Contigs > 500b ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500b ratio - Count" ) ), file=sys.stdout)
		print("Contigs > 1Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Kb ratio - Count" ) ), file=sys.stdout)
		print("Contigs > 5Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Kb ratio - Count" ) ), file=sys.stdout)
		print("Contigs > 10Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 10Kb ratio - Count" ) ), file=sys.stdout)
		print("Contigs > 50Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 50Kb ratio - Count" ) ), file=sys.stdout)
		print("Contigs > 100Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100Kb ratio - Count" ) ), file=sys.stdout)
		print("Contigs > 500Kb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500Kb ratio - Count" ) ), file=sys.stdout)
		print("Contigs > 1Mb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Mb ratio - Count" ) ), file=sys.stdout)
		print("Contigs > 5Mb ratio - Count" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Mb ratio - Count" ) ), file=sys.stdout)

		print("# Contig cumulative length ratio by sequence size", file=sys.stdout)
		print("Contigs > 100b ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100b ratio - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 500b ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500b ratio - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 1Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Kb ratio - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 5Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Kb ratio - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 10Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 10Kb ratio - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 50Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 50Kb ratio - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 100Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 100Kb ratio - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 500Kb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 500Kb ratio - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 1Mb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 1Mb ratio - Cumulative length" ) ), file=sys.stdout)
		print("Contigs > 5Mb ratio - Cumulative length" + "\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_len_stats" , "Sequences > 5Mb ratio - Cumulative length" ) ), file=sys.stdout)

		print("# Contig N-values length", file=sys.stdout)
		# N values, length and index
		for value in Nsize_list :
			# Print N length
			print("N" + str(value) + " Length\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_Nvalue_Length" , value ) ), file=sys.stdout)
		print("# Contig N-values index", file=sys.stdout)
		for value in Nsize_list :
			# Print N index
			print("N" + str(value) + " Index\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_Nvalue_Index" , value ) ), file=sys.stdout)
			
		# NG values, length and index
		if options.genome_size :
			print("# Contig NG-values length (based on and expected genome size of " + options.genome_size + "bp)", file=sys.stdout)
			for value in NGsize_list :
				# Print NG length
				print("NG" + str(value) + " Length\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_NGvalue_Length" , value ) ), file=sys.stdout)
			print("# Contig NG-values index (based on and expected genome size of " + options.genome_size + "bp)", file=sys.stdout)
			for value in NGsize_list :
				# Print NG index
				print("NG" + str(value) + " Index\t" + "\t".join( list_from_stats_db( results_db , id_list , "contig_NGvalue_Index" , value ) ), file=sys.stdout)


if __name__ == '__main__':
	main()
