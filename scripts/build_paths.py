#!/usr/bin/env python3
"""
build_paths.py — Marker QC, chimera detection, graph construction, and DAG
tiling path selection for both haplotypes.

Supports three modes:
  - Markers only  (-n + -m): genetic-map DAG with weight='marker_num'
  - Reference only (-g):     alignment DAG with weight='align' via minimap2/PAF
  - Both (-n + -m + -g):    markers build paths, nucmer resolves undirected orientations

Outputs:
  {out}.1.list               — Hap1 tiling paths (chr_id<TAB>seq|strand,...)
  {out}.2.list               — Hap2 tiling paths (optional, unless --N2)
  {out}.Un.list              — Comma-separated unplaced sequence IDs
  {out}.unused_sequences.list — Sequences assigned to a chr but not placed
  {out}.missing_orientation.hap1.txt
  {out}.missing_orientation.hap2.txt
  {out}.unknown_markers.txt  (if any)
"""

import argparse
import os
import subprocess
import sys
import datetime
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *
from lib_files.GFF_lib import *
from Bio.Seq import Seq


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__)

    # Required
    parser.add_argument('-i', '--input', dest='query', required=True,
        metavar='query.fasta',
        help='FASTA file of the input sequences')
    parser.add_argument('-n', '--map', dest='markers_hits',
        metavar='map.bed',
        help='Marker positions on input sequences (BED-like, 4 columns)')
    parser.add_argument('-m', '--markers', dest='marker_map',
        metavar='markers_list',
        help='Marker genetic map (chr_id, position, marker_id)')
    parser.add_argument('-o', '--out', dest='out', default='out',
        metavar='NAME',
        help='Output files prefix [default: out]')

    # Reference genome inputs
    parser.add_argument('-g', '--guide', dest='reference',
        metavar='reference.fasta',
        help='Guide/reference genome FASTA')
    parser.add_argument('-l', '--local_alignment', dest='hits', default='local.paf',
        metavar='local.paf',
        help='PAF alignment file (used when --align is not set) [default: local.paf]')
    parser.add_argument('--align', dest='map', action='store_true',
        help='Run minimap2 to generate PAF alignment against reference')
    parser.add_argument('-t', '--tool', dest='mapping', default=' --cs -x asm20 -r 1000',
        metavar='OPTS',
        help='Minimap2 options [default: --cs -x asm20 -r 1000]')
    parser.add_argument('--hitgap', dest='hitgap', default='100000',
        metavar='N',
        help='Max gap between alignment hits to merge [default: 100000]')
    parser.add_argument('--distance1', dest='distance1', default='2000000',
        metavar='N',
        help='Max Hap1 adjacency distance in reference coords [default: 2000000]')
    parser.add_argument('--distance2', dest='distance2', default='4000000',
        metavar='N',
        help='Max Hap2 adjacency distance in reference coords [default: 4000000]')
    parser.add_argument('--reuse_intermediate', dest='reuse_intermediate',
        action='store_true',
        help='Reuse existing nucmer intermediate files')

    # Optional inputs
    parser.add_argument('-a', '--input_agp', dest='input_agp',
        metavar='query.agp',
        help='AGP file of query sequence composition')
    parser.add_argument('--GFF3', dest='gff3',
        metavar='genes.gff3',
        help='Annotation in GFF3 format (used for chimera QC only)')
    parser.add_argument('-e', '--exclusion', dest='exclusion',
        metavar='exclusion.tsv',
        help='Mutually exclusive sequence pairs (tab-separated)')
    parser.add_argument('-k', '--known', dest='known',
        metavar='known.tsv',
        help='Sequences known to be in the same haplotype')
    parser.add_argument('--alternative_groups', dest='alternative',
        metavar='alternative_groups.tsv',
        help='Alternative haplotype sequence pairs')
    parser.add_argument('--R1', dest='Require1',
        metavar='1st.txt',
        help='Sequences required in haplotype 1')
    parser.add_argument('--R2', dest='Require2',
        metavar='2nd.txt',
        help='Sequences required in haplotype 2')
    parser.add_argument('--F1', dest='force_direction1',
        default=False, action='store_true',
        help='Force marker direction for --R1 sequences')
    parser.add_argument('--F2', dest='force_direction2',
        default=False, action='store_true',
        help='Force marker direction for --R2 sequences')
    parser.add_argument('--min1', dest='minR1', default='0', metavar='N',
        help='Minimum length for --R1 sequences [default: 0]')
    parser.add_argument('--min2', dest='minR2', default='0', metavar='N',
        help='Minimum length for --R2 sequences [default: 0]')
    parser.add_argument('--B1', '--blacklist1', dest='Blacklist1',
        metavar='blacklist1.txt',
        help='Blacklisted sequences for haplotype 1')
    parser.add_argument('--B2', '--blacklist2', dest='Blacklist2',
        metavar='blacklist2.txt',
        help='Blacklisted sequences for haplotype 2')
    parser.add_argument('--input_groups', dest='input_groups',
        metavar='input_groups.tsv',
        help='Sequence group associations (passed through; not used here)')

    # Behaviour flags
    parser.add_argument('-c', '--cores', dest='cores', default=4,
        metavar='N', type=int,
        help='Cores for parallel steps [default: 4]')
    parser.add_argument('-f', '--filter', dest='filter_hits',
        default=False, action='store_true',
        help='Remove sequences with intra-sequence marker duplication')
    parser.add_argument('--extended', dest='extended_region',
        default=False, action='store_true',
        help='Extend sequence association to all markers on the chromosome')
    parser.add_argument('--N2', dest='No2',
        default=False, action='store_true',
        help="Don't build the second haplotype path")
    parser.add_argument('--skip_chimeric_qc', dest='skip_qc',
        default=False, action='store_true',
        help='Skip chimeric sequence QC')
    parser.add_argument('--disable_marker_ploidy_check',
        dest='disable_marker_ploidy_check',
        default=False, action='store_true')
    parser.add_argument('--conflict', dest='conflict_resolution',
        default='exit',
        help='Conflict resolution policy: exit | ignore | release')
    parser.add_argument('-v', '--dry', dest='dry',
        default=False, action='store_true',
        help='Dry run: check marker uniqueness and quit')
    parser.add_argument('--allow_rearrangements', dest='rearrangements',
        default=False, action='store_true')
    parser.add_argument('-1', '--path1', dest='use1',
        metavar='path1.txt',
        help='Pre-defined Hap1 tiling path (skips DAG search; chr_id<TAB>seq|strand,...)')
    parser.add_argument('-2', '--path2', dest='use2',
        metavar='path2.txt',
        help='Pre-defined Hap2 tiling path (skips DAG search; chr_id<TAB>seq|strand,...)')
    parser.add_argument('--required_as_path', dest='forced_as_path',
        default=False, action='store_true',
        help='Treat --R1/--R2 lists as full tiling paths (used with -g to skip marker DAG)')
    parser.add_argument('--only_markers', dest='only_markers',
        default=False, action='store_true',
        help='Limit rejected-sequence QC to marker-bearing sequences only')
    parser.add_argument('--avoid_rejected_qc', dest='avoidrejectedqc',
        default=False, action='store_true',
        help='Skip QC of sequences assigned to a chromosome but not placed')
    parser.add_argument('--conf', dest='conf',
        default=None,
        help='Path to HaploSync.conf.toml [default: auto-detect]')

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    options = parser.parse_args()

    # Validate mode
    use_markers = bool(options.marker_map and options.markers_hits)
    use_reference = bool(options.reference)
    if not use_markers and not use_reference and not options.forced_as_path:
        print('[ERROR] Provide (-m AND -n) and/or -g/--guide', file=sys.stderr)
        sys.exit(1)
    if bool(options.marker_map) != bool(options.markers_hits):
        print('[ERROR] -m and -n must be given together', file=sys.stderr)
        sys.exit(1)
    if options.forced_as_path and not use_reference:
        print('[ERROR] --required_as_path requires -g/--guide', file=sys.stderr)
        sys.exit(1)

    # Locate conf file
    conf_path = options.conf or os.path.join(
        os.path.dirname(os.path.realpath(__file__)), '..', 'HaploSync.conf.toml')
    paths = set_paths(conf_path)

    print("Running build_paths.py from HaploSync version " + get_version(), file=sys.stdout)
    print("Command: " + " ".join(sys.argv), file=sys.stdout)
    print("----", file=sys.stdout)

    conflict_resolution = options.conflict_resolution.lower()

    # -------------------------------------------------------------------------
    # Load query FASTA
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Reading query sequences', file=sys.stdout)
    query_fasta_db = read_fasta(options.query)
    query_len = {}
    all_seq_length = {}
    for query_name in query_fasta_db:
        query_len[query_name] = int(len(query_fasta_db[query_name]))
        all_seq_length[query_name]        = query_len[query_name]
        all_seq_length[query_name + '|+'] = query_len[query_name]
        all_seq_length[query_name + '|-'] = query_len[query_name]
        all_seq_length[query_name + '|.'] = query_len[query_name]
    print('[' + str(datetime.datetime.now()) + '] === ' + str(len(query_len)) + ' query sequences', file=sys.stdout)

    # -------------------------------------------------------------------------
    # Load input AGP (optional)
    # -------------------------------------------------------------------------
    if options.input_agp:
        print('[' + str(datetime.datetime.now()) + '] = Reading query AGP', file=sys.stdout)
        query_agp_db = read_agp(options.input_agp)
    else:
        query_agp_db = {}

    # -------------------------------------------------------------------------
    # Load GFF3 (optional, for chimera QC)
    # -------------------------------------------------------------------------
    if options.gff3:
        print('[' + str(datetime.datetime.now()) + '] = Reading GFF3 annotation', file=sys.stdout)
        annotation_gff3, mRNA_db = read_gff3(options.gff3)
    else:
        annotation_gff3 = {}
        mRNA_db = {}

    # -------------------------------------------------------------------------
    # Load marker map and marker hits (markers mode only)
    # -------------------------------------------------------------------------
    reference_sequence_list = []
    marker_hits_by_seq = {}
    marker_hits_by_id = {}
    marker_map_by_seq = {}
    marker_map_by_id = {}
    unknown_markers = []

    if use_markers:
        print('[' + str(datetime.datetime.now()) + '] = Reading markers map', file=sys.stdout)
        for line in open(options.marker_map):
            if line == '' or line[0] == '#':
                continue
            try:
                seq_id, pos, marker_id = line.rstrip().split('\t')
            except:
                print('[ERROR] Marker map format error: ' + line.rstrip(), file=sys.stderr)
                sys.exit(1)
            if seq_id not in marker_map_by_seq:
                reference_sequence_list.append(seq_id)
                marker_map_by_seq[seq_id] = []
            marker_map_by_seq[seq_id].append([int(pos), marker_id])
            marker_map_by_id[marker_id] = [seq_id, int(pos), marker_id]

        print('[' + str(datetime.datetime.now()) + '] = Reading marker hits', file=sys.stdout)
        for line in open(options.markers_hits):
            if line == '' or line[0] == '#':
                continue
            try:
                seq_id, pos_1, pos_2, marker_id = line.rstrip().split('\t')
            except:
                print('[ERROR] Marker hits format error (expected: Seq_ID<TAB>Start<TAB>Stop<TAB>Marker_ID): ' + line.rstrip(), file=sys.stderr)
                sys.exit(1)
            if marker_id in marker_map_by_id:
                if seq_id not in marker_hits_by_seq:
                    marker_hits_by_seq[seq_id] = []
                marker_chr = marker_map_by_id[marker_id][0]
                marker_pos = marker_map_by_id[marker_id][1]
                start = min(int(pos_1), int(pos_2))
                stop  = max(int(pos_1), int(pos_2))
                marker_hits_by_seq[seq_id].append([int(start), int(stop), marker_id, marker_chr, int(marker_pos)])
                if marker_id not in marker_hits_by_id:
                    marker_hits_by_id[marker_id] = []
                marker_hits_by_id[marker_id].append([seq_id, int(start), int(stop)])
            else:
                unknown_markers.append(marker_id)

        if unknown_markers:
            fname = options.out + '.unknown_markers.txt'
            with open(fname, 'w') as f:
                for m in sorted(set(unknown_markers)):
                    print(m, file=f)
            print('[WARNING] ' + str(len(set(unknown_markers))) + ' unknown marker IDs written to ' + fname, file=sys.stderr)

    # -------------------------------------------------------------------------
    # Load reference FASTA (reference mode)
    # -------------------------------------------------------------------------
    reference = {}
    reference_len = {}
    if use_reference:
        print('[' + str(datetime.datetime.now()) + '] = Reading reference sequences', file=sys.stdout)
        reference = read_fasta(options.reference)
        for ref_name in reference:
            reference_len[ref_name] = int(len(reference[ref_name]))
            all_seq_length[ref_name]        = reference_len[ref_name]
            all_seq_length[ref_name + '|+'] = reference_len[ref_name]
            all_seq_length[ref_name + '|-'] = reference_len[ref_name]
            all_seq_length[ref_name + '|.'] = reference_len[ref_name]
        if not use_markers:
            reference_sequence_list = list(reference.keys())

    # -------------------------------------------------------------------------
    # Load optional constraint files
    # -------------------------------------------------------------------------
    unwanted_pairs = {}
    if options.exclusion:
        print('[' + str(datetime.datetime.now()) + '] = Reading exclusion pairs', file=sys.stdout)
        for line in open(options.exclusion):
            if line == '' or line[0] == '#':
                continue
            try:
                seq_1_id, seq_2_id = line.rstrip().split('\t')
            except:
                print('[ERROR] Exclusion pair file format error: ' + line.rstrip(), file=sys.stderr)
                sys.exit(1)
            for sid in (seq_1_id, seq_2_id):
                other = seq_2_id if sid == seq_1_id else seq_1_id
                for o in ('|+', '|-', '|.'):
                    if (sid + o) not in unwanted_pairs:
                        unwanted_pairs[sid + o] = []
                    unwanted_pairs[sid + o] += three_orientation_list(other, False)

    known_groups = {}
    known_groups_by_seqid = {}
    if options.known:
        print('[' + str(datetime.datetime.now()) + '] = Reading known haplotype groups', file=sys.stdout)
        for line in open(options.known):
            if line == '' or line[0] == '#':
                continue
            try:
                seq_id, group_id = line.rstrip().split('\t')
            except:
                print('[ERROR] Known groups file format error: ' + line.rstrip(), file=sys.stderr)
                sys.exit(1)
            if group_id not in known_groups:
                known_groups[group_id] = []
            known_groups[group_id].append(seq_id)
        for group_id in list(known_groups.keys()):
            for seq_id in known_groups[group_id]:
                if seq_id not in known_groups_by_seqid:
                    known_groups_by_seqid[seq_id]       = []
                    known_groups_by_seqid[seq_id + '|+'] = []
                    known_groups_by_seqid[seq_id + '|-'] = []
                    known_groups_by_seqid[seq_id + '|.'] = []
                for seq_2_id in known_groups[group_id]:
                    for key in (seq_id, seq_id+'|+', seq_id+'|-', seq_id+'|.'):
                        known_groups_by_seqid[key].append(seq_2_id)
                        known_groups_by_seqid[key] += three_orientation_list(seq_2_id, False)

    alternative_sequences = {}
    if options.alternative:
        for line in open(options.alternative):
            if line[0] == '#' or line.rstrip() == '':
                continue
            group_1, group_2 = line.rstrip().split('\t')
            group_1_list = group_1.split(',')
            group_2_list = group_2.split(',')
            for seq_id in group_1_list:
                if seq_id not in alternative_sequences:
                    alternative_sequences[seq_id] = []
                alternative_sequences[seq_id] += group_2_list
                alternative_sequences[seq_id] = list(set(alternative_sequences[seq_id]))
            for seq_id in group_2_list:
                if seq_id not in alternative_sequences:
                    alternative_sequences[seq_id] = []
                alternative_sequences[seq_id] += group_1_list
                alternative_sequences[seq_id] = list(set(alternative_sequences[seq_id]))

    # Required path lists
    forced_list_1 = defaultdict(list)
    forced_list_2 = defaultdict(list)
    forced_list_1_id = []
    forced_list_2_id = []
    discarded = []

    if options.Require1 or options.Require2:
        print('[' + str(datetime.datetime.now()) + '] = Reading required sequence lists', file=sys.stdout)

    if options.Require1:
        for line in open(options.Require1):
            try:
                Target_sequence, loaded_path = line.rstrip().split('\t')
                prompted_list = loaded_path.split(',')
                forced_list_1[Target_sequence] = []
                for id_ in prompted_list:
                    id_name = id_[:-2]
                    if query_len.get(id_name, 0) >= int(options.minR1):
                        forced_list_1[Target_sequence].append(id_)
                        forced_list_1_id.append(id_name)
                    else:
                        discarded.append(id_name)
            except:
                pass

    if options.Require2:
        for line in open(options.Require2):
            try:
                Target_sequence, loaded_path = line.rstrip().split('\t')
                prompted_list = loaded_path.split(',')
                forced_list_2[Target_sequence] = []
                for id_ in prompted_list:
                    id_name = id_[:-2]
                    if query_len.get(id_name, 0) >= int(options.minR2):
                        forced_list_2[Target_sequence].append(id_)
                        forced_list_2_id.append(id_name)
                    else:
                        discarded.append(id_name)
            except:
                pass
        doubled = [s for s in forced_list_2_id if s in forced_list_1_id]
        if doubled:
            print('[ERROR] Same sequence(s) required for both haplotypes: ' + ', '.join(doubled), file=sys.stderr)
            sys.exit(1)

    # Blacklists
    blacklist_1 = defaultdict(list)
    blacklist_2 = defaultdict(list)
    for ref_name in reference_sequence_list:
        blacklist_1[ref_name] = []
        blacklist_2[ref_name] = []

    if options.Blacklist1:
        for line in open(options.Blacklist1):
            try:
                Tid, Qids = line.rstrip().split('\t')
                if Tid in blacklist_1:
                    for name in Qids.split(','):
                        blacklist_1[Tid] += three_orientation_list(name)
            except:
                pass

    for key_in in sorted(blacklist_1.keys()):
        Required_somewhere_else = []
        for key_out in sorted(blacklist_1.keys()):
            if key_out != key_in:
                Required_somewhere_else += list(forced_list_1[key_out])
            Required_somewhere_else += list(forced_list_2[key_out])
        for name in list(set(Required_somewhere_else)):
            blacklist_1[key_in] += three_orientation_list(name)

    if options.Blacklist2:
        for line in open(options.Blacklist2):
            try:
                Tid, Qids = line.rstrip().split('\t')
                if Tid in blacklist_2:
                    for name in Qids.split(','):
                        blacklist_2[Tid] += three_orientation_list(name)
            except:
                pass

    for key_in in sorted(blacklist_2.keys()):
        Required_somewhere_else = []
        for key_out in sorted(blacklist_2.keys()):
            if key_out != key_in:
                Required_somewhere_else += list(forced_list_2[key_out])
            Required_somewhere_else += list(forced_list_1[key_out])
        for name in list(set(Required_somewhere_else)):
            blacklist_2[key_in] += three_orientation_list(name)

    # =========================================================================
    # Path building
    # =========================================================================
    best_1_paths = {}
    best_1_paths_edges = {}
    best_2_paths = {}
    best_2_paths_edges = {}
    all_used = []
    unused_by_chr = {}

    # Pre-load user-defined paths if provided (skips DAG search for those chromosomes)
    if options.use1:
        print('[' + str(datetime.datetime.now()) + '] = Loading pre-defined Hap1 path from ' + options.use1, file=sys.stdout)
        for line in open(options.use1):
            try:
                Target_sequence, loaded_path = line.rstrip().split('\t')
                prompted_list = loaded_path.split(',')
                best_1_paths_edges[Target_sequence] = prompted_list
                all_used += prompted_list
            except:
                pass

    if options.use2:
        print('[' + str(datetime.datetime.now()) + '] = Loading pre-defined Hap2 path from ' + options.use2, file=sys.stdout)
        for line in open(options.use2):
            try:
                Target_sequence, loaded_path = line.rstrip().split('\t')
                prompted_list = loaded_path.split(',')
                best_2_paths_edges[Target_sequence] = prompted_list
                all_used += prompted_list
            except:
                pass

    if use_markers:
        # -------------------------------------------------------------------------
        # Marker QC
        # -------------------------------------------------------------------------
        print('[' + str(datetime.datetime.now()) + '] = Testing markers uniqueness', file=sys.stdout)
        if not options.disable_marker_ploidy_check:
            ploidy = 1 if options.No2 else 2
        else:
            ploidy = 100

        multi_copy_markers, unique_marker_hits_by_id = check_marker_copies(marker_hits_by_id, ploidy)
        chimera_list, markers_itradup, unique_distinct_marker_hits_by_id = check_in_sequence_duplications(
            marker_hits_by_seq, unique_marker_hits_by_id, multi_copy_markers,
            query_fasta_db, annotation_gff3, query_agp_db,
            options.out, options.cores, paths, options.skip_qc, ploidy)

        if options.dry:
            print('[' + str(datetime.datetime.now()) + '] = Dry run completed', file=sys.stdout)
            sys.exit(0)

        print('[' + str(datetime.datetime.now()) + '] = Cleaning markers', file=sys.stdout)
        clean_marker_set_by_seq = clean_markers(
            unique_distinct_marker_hits_by_id, chimera_list,
            marker_map_by_seq, marker_map_by_id,
            options.out, query_fasta_db,
            {'show-coords': paths['show-coords'], 'nucmer': paths['nucmer']},
            options.cores, options.filter_hits, options.extended_region,
            forced_list_1, forced_list_2,
            options.force_direction1, options.force_direction2)

        # -------------------------------------------------------------------------
        # Build tiling paths (markers mode)
        # -------------------------------------------------------------------------
        print('[' + str(datetime.datetime.now()) + '] = Building tiling paths on markers', file=sys.stdout)
        clean_marker_set_by_chr = {}
        for seq_id in list(clean_marker_set_by_seq.keys()):
            chr_id = clean_marker_set_by_seq[seq_id]['chr'][0]
            if chr_id not in clean_marker_set_by_chr:
                clean_marker_set_by_chr[chr_id] = []
            clean_marker_set_by_chr[chr_id].append(clean_marker_set_by_seq[seq_id])

        # Pre-compute alternative-to-forced blacklists
        unwanted_sequences_alternative_to_forced = {}
        for Target_sequence in list(forced_list_1.keys()):
            for seq_id in forced_list_1[Target_sequence]:
                if seq_id in alternative_sequences:
                    for Target_sequence_2 in list(forced_list_1.keys()):
                        if Target_sequence_2 not in unwanted_sequences_alternative_to_forced:
                            unwanted_sequences_alternative_to_forced[Target_sequence_2] = []
                        unwanted_sequences_alternative_to_forced[Target_sequence_2] += alternative_sequences[seq_id]

        # --- Hap1 ---
        print('[' + str(datetime.datetime.now()) + '] = Hap1', file=sys.stdout)
        for chr_id in sorted(clean_marker_set_by_chr.keys()):
            if chr_id in best_1_paths_edges:
                print('[' + str(datetime.datetime.now()) + '] == ' + chr_id + ' pre-loaded via --path1, skipping DAG', file=sys.stdout)
                continue
            print('[' + str(datetime.datetime.now()) + '] == Processing ' + chr_id, file=sys.stdout)
            stop_pos = max([int(x[0]) for x in marker_map_by_seq[chr_id]]) + 1
            oriented_marker_set = add_orientation_to_ID(clean_marker_set_by_chr[chr_id])
            marker_set, validated_forced_list_1, validated_forced_list_2, validated_blacklist_1, validated_blacklist_2 = \
                validate_marker_set(oriented_marker_set, forced_list_1[chr_id], forced_list_2[chr_id],
                                    blacklist_1[chr_id], blacklist_2[chr_id], options.conflict_resolution)

            preferred_db = {}
            unwanted_to_remove_list = list(set(
                (unwanted_sequences_alternative_to_forced.get(chr_id, []))))
            unwanted_to_reuse_list = []
            new_forced = validated_forced_list_1[:]
            best_found = False

            while not best_found:
                marker_graph_1 = remove_sequence_from_graph(
                    unwanted_to_remove_list, chr_id, stop_pos, marker_set, new_forced, validated_blacklist_1)
                best_1_nodes = nx.dag_longest_path(marker_graph_1, weight='marker_num')
                best_1_edges, best_1_edges_names = get_marker_subgraph_from_path(marker_graph_1, best_1_nodes)
                unwanted_to_remove = search_unwanted(best_1_edges_names, unwanted_pairs, all_seq_length, validated_forced_list_1)

                if unwanted_to_remove == '':
                    all_preferred_used = True
                    to_reintegrate = []
                    preferred_to_remove = []
                    for seq_id in list(preferred_db.keys()):
                        if seq_id not in best_1_edges_names:
                            all_preferred_used = False
                            unwanted_to_reuse_list += three_orientation_list(seq_id, True)
                            preferred_to_remove.append(seq_id)
                            to_reintegrate += preferred_db[seq_id]
                    if all_preferred_used:
                        best_found = True
                    else:
                        new_forced = best_1_edges_names
                        for reintegrate_id in list(set(to_reintegrate)):
                            if reintegrate_id in unwanted_to_remove_list and reintegrate_id not in unwanted_to_reuse_list:
                                unwanted_to_remove_list.remove(reintegrate_id)
                        for seq_id in preferred_to_remove:
                            del preferred_db[seq_id]
                            unwanted_to_remove_list += three_orientation_list(seq_id, True)
                        unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                else:
                    updated_forced = [e for e in new_forced if e not in unwanted_to_remove['remove']]
                    new_forced = updated_forced
                    to_reintegrate = []
                    for key in unwanted_to_remove['remove']:
                        if key in preferred_db:
                            to_reintegrate += preferred_db[key]
                            del preferred_db[key]
                    if to_reintegrate:
                        unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                        for reintegrate_id in to_reintegrate:
                            if reintegrate_id in unwanted_to_remove_list and reintegrate_id not in unwanted_to_reuse_list:
                                unwanted_to_remove_list.remove(reintegrate_id)
                                new_forced.append(reintegrate_id)
                    unwanted_to_remove_list += unwanted_to_remove['remove']
                    unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                    if unwanted_to_remove['keep'] not in preferred_db:
                        preferred_db[unwanted_to_remove['keep']] = []
                    preferred_db[unwanted_to_remove['keep']] += unwanted_to_remove['remove']

            best_1_paths[chr_id] = best_1_edges
            used, best_1_paths_edges[chr_id] = make_list_from_marker_path(best_1_edges)
            all_used += used

        # --- Hap2 ---
        if not options.No2:
            paired_to_used = []
            if not options.rearrangements:
                for Target_sequence in list(best_1_paths_edges.keys()):
                    for seq_id in best_1_paths_edges[Target_sequence]:
                        if seq_id in known_groups_by_seqid:
                            paired_to_used += known_groups_by_seqid[seq_id]
                paired_to_used = list(set(paired_to_used))

            unwanted_sequences_alternative_to_used = {}
            for Target_sequence in list(best_1_paths_edges.keys()):
                unwanted_sequences_alternative_to_used[Target_sequence] = []
                for seq_id in best_1_paths_edges[Target_sequence]:
                    if seq_id in alternative_sequences:
                        for Target_sequence_2 in list(best_1_paths_edges.keys()):
                            if Target_sequence_2 != Target_sequence:
                                unwanted_sequences_alternative_to_used.setdefault(Target_sequence_2, [])
                                unwanted_sequences_alternative_to_used[Target_sequence_2] += alternative_sequences[seq_id]

            unwanted_sequences_alternative_to_forced_2 = {}
            for Target_sequence in list(forced_list_2.keys()):
                for seq_id in forced_list_2[Target_sequence]:
                    if seq_id in alternative_sequences:
                        for Target_sequence_2 in list(forced_list_2.keys()):
                            if Target_sequence_2 not in unwanted_sequences_alternative_to_forced_2:
                                unwanted_sequences_alternative_to_forced_2[Target_sequence_2] = []
                            unwanted_sequences_alternative_to_forced_2[Target_sequence_2] += alternative_sequences[seq_id]

            print('[' + str(datetime.datetime.now()) + '] = Hap2', file=sys.stdout)
            for chr_id in sorted(clean_marker_set_by_chr.keys()):
                if chr_id in best_2_paths_edges:
                    print('[' + str(datetime.datetime.now()) + '] == ' + chr_id + ' pre-loaded via --path2, skipping DAG', file=sys.stdout)
                    continue
                print('[' + str(datetime.datetime.now()) + '] == Processing ' + chr_id, file=sys.stdout)
                stop_pos = max([int(x[0]) for x in marker_map_by_seq[chr_id]]) + 1
                oriented_marker_set = add_orientation_to_ID(clean_marker_set_by_chr[chr_id])
                marker_set, validated_forced_list_1, validated_forced_list_2, validated_blacklist_1, validated_blacklist_2 = \
                    validate_marker_set(oriented_marker_set, forced_list_1[chr_id], forced_list_2[chr_id],
                                        blacklist_1[chr_id], blacklist_2[chr_id], options.conflict_resolution)

                preferred_db = {}
                unwanted_to_remove_list = list(paired_to_used)
                unwanted_to_remove_list += unwanted_sequences_alternative_to_used.get(chr_id, [])
                unwanted_to_remove_list += unwanted_sequences_alternative_to_forced_2.get(chr_id, [])
                for seq_id in best_1_paths_edges[chr_id]:
                    unwanted_to_remove_list += three_orientation_list(seq_id, True)
                    if seq_id in known_groups_by_seqid:
                        unwanted_to_remove_list += known_groups_by_seqid[seq_id]
                unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                unwanted_to_reuse_list = []
                new_forced = validated_forced_list_2[:]
                best_found = False

                while not best_found:
                    marker_graph_2 = remove_sequence_from_graph(
                        unwanted_to_remove_list, chr_id, stop_pos, marker_set, new_forced, validated_blacklist_2)
                    best_2_nodes = nx.dag_longest_path(marker_graph_2, weight='marker_num')
                    best_2_edges, best_2_edges_names = get_marker_subgraph_from_path(marker_graph_2, best_2_nodes)
                    unwanted_to_remove = search_unwanted(best_2_edges_names, unwanted_pairs, all_seq_length, validated_forced_list_2)

                    if unwanted_to_remove == '':
                        all_preferred_used = True
                        to_reintegrate = []
                        preferred_to_remove = []
                        for seq_id in list(preferred_db.keys()):
                            if seq_id not in best_2_edges_names:
                                all_preferred_used = False
                                unwanted_to_reuse_list += three_orientation_list(seq_id, True)
                                preferred_to_remove.append(seq_id)
                                to_reintegrate += preferred_db[seq_id]
                        if all_preferred_used:
                            best_found = True
                        else:
                            new_forced = best_2_edges_names
                            for reintegrate_id in list(set(to_reintegrate)):
                                if reintegrate_id in unwanted_to_remove_list and reintegrate_id not in unwanted_to_reuse_list:
                                    unwanted_to_remove_list.remove(reintegrate_id)
                            for seq_id in preferred_to_remove:
                                del preferred_db[seq_id]
                                unwanted_to_remove_list += three_orientation_list(seq_id, True)
                            unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                    else:
                        updated_forced = [e for e in new_forced if e not in unwanted_to_remove['remove']]
                        new_forced = updated_forced
                        to_reintegrate = []
                        for key in unwanted_to_remove['remove']:
                            if key in preferred_db:
                                to_reintegrate += preferred_db[key]
                                del preferred_db[key]
                        if to_reintegrate:
                            unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                            for reintegrate_id in to_reintegrate:
                                if reintegrate_id in unwanted_to_remove_list and reintegrate_id not in unwanted_to_reuse_list:
                                    unwanted_to_remove_list.remove(reintegrate_id)
                        unwanted_to_remove_list += unwanted_to_remove['remove']
                        unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                        if unwanted_to_remove['keep'] not in preferred_db:
                            preferred_db[unwanted_to_remove['keep']] = []
                        preferred_db[unwanted_to_remove['keep']] += unwanted_to_remove['remove']

                best_2_paths[chr_id] = best_2_edges
                used, best_2_paths_edges[chr_id] = make_list_from_marker_path(best_2_edges)
                all_used += used
                unused_by_chr[chr_id] = [e['id'] for e in marker_set if e['id'] not in all_used]

    else:
        # -------------------------------------------------------------------------
        # Reference-only path building
        # -------------------------------------------------------------------------
        print('[' + str(datetime.datetime.now()) + '] = Building tiling paths from reference alignment', file=sys.stdout)

        query_base = os.path.basename(options.query)
        hits_base = os.path.basename(options.hits)

        doubleQuery = {}
        doubleQuery_len = {}

        if options.map:
            # Run minimap2: forward pass
            print('[' + str(datetime.datetime.now()) + '] = Mapping query sequences on reference', file=sys.stdout)

            query_for = 'tmp.' + query_base + '.for'
            with open(query_for, 'w') as qf_file:
                for query_seq in query_fasta_db:
                    print('>' + query_seq + '|+', file=qf_file)
                    print(str(query_fasta_db[query_seq]).upper(), file=qf_file)
                    doubleQuery[query_seq + '|+'] = str(query_fasta_db[query_seq]).upper()
                    doubleQuery_len[query_seq + '|+'] = int(len(str(query_fasta_db[query_seq])))

            minimap_path = paths.get('minimap2', '')
            if minimap_path == '':
                minimap2_search = subprocess.Popen('which minimap2', shell=True, stdout=subprocess.PIPE, text=True)
                command_line, _ = minimap2_search.communicate()
                command_line = command_line.rstrip()
                if command_line == '':
                    print('[ERROR] Minimap expected to be in $PATH, not found', file=sys.stderr)
                    sys.exit(1)
            else:
                command_line = minimap_path + '/minimap2'

            if not os.path.exists(command_line):
                print('[ERROR] wrong or no path to minimap2', file=sys.stderr)
                sys.exit(1)

            minimap2_command = command_line + ' ' + options.mapping + ' '

            print('[' + str(datetime.datetime.now()) + '] == Mapping forward query sequences', file=sys.stdout)
            map_for_path = 'tmp.' + hits_base + '.for'
            mapping_cmd = (minimap2_command + ' --for-only -t ' + str(options.cores) +
                           ' ' + options.reference + ' ' + query_for + " | awk '$5==\"+\"'")
            print('[' + str(datetime.datetime.now()) + '] === Command line: ' + mapping_cmd + ' > ' + map_for_path, file=sys.stdout)
            with open(map_for_path, 'w') as map_for_file:
                map_proc = subprocess.Popen(mapping_cmd, shell=True, stdout=map_for_file)
                map_proc.communicate()

            # Reverse pass
            print('[' + str(datetime.datetime.now()) + '] == Reversing query sequences', file=sys.stdout)
            query_rev = 'tmp.' + query_base + '.rev'
            with open(query_rev, 'w') as qr_file:
                for query_seq in query_fasta_db:
                    print('>' + query_seq + '|-', file=qr_file)
                    print(str(Seq(query_fasta_db[query_seq]).reverse_complement()).upper(), file=qr_file)
                    doubleQuery[query_seq + '|-'] = str(Seq(query_fasta_db[query_seq]).reverse_complement()).upper()
                    doubleQuery_len[query_seq + '|-'] = int(len(str(query_fasta_db[query_seq])))

            print('[' + str(datetime.datetime.now()) + '] == Mapping reverse query sequences', file=sys.stdout)
            map_rev_path = 'tmp.' + hits_base + '.rev'
            mapping_cmd = (minimap2_command + ' --for-only -t ' + str(options.cores) +
                           ' ' + options.reference + ' ' + query_rev + " | awk '$5==\"+\"'")
            print('[' + str(datetime.datetime.now()) + '] === Command line: ' + mapping_cmd + ' > ' + map_rev_path, file=sys.stdout)
            with open(map_rev_path, 'w') as map_rev_file:
                map_proc = subprocess.Popen(mapping_cmd, shell=True, stdout=map_rev_file)
                map_proc.communicate()

            # Merge
            print('[' + str(datetime.datetime.now()) + '] == Merge forward and reverse alignments', file=sys.stdout)
            map_multi = 'tmp.' + hits_base + '.multimapping'
            with open(map_multi, 'w') as map_multi_file:
                merge_proc = subprocess.Popen('cat ' + map_for_path + ' ' + map_rev_path,
                                              shell=True, stdout=map_multi_file)
                merge_proc.communicate()

            cp_proc = subprocess.Popen('cp ' + map_multi + ' ' + options.hits, shell=True)
            cp_proc.communicate()

        else:
            # Use existing PAF
            map_multi = 'tmp.' + hits_base + '.multimapping'
            cp_proc = subprocess.Popen('cp ' + options.hits + ' ' + map_multi, shell=True)
            cp_proc.communicate()

            for query_seq in query_fasta_db:
                doubleQuery[query_seq + '|+'] = str(query_fasta_db[query_seq]).upper()
                doubleQuery_len[query_seq + '|+'] = int(len(doubleQuery[query_seq + '|+']))
                doubleQuery[query_seq + '|-'] = str(Seq(query_fasta_db[query_seq]).reverse_complement()).upper()
                doubleQuery_len[query_seq + '|-'] = int(len(doubleQuery[query_seq + '|-']))

        # Merge and unique hits
        unique_hits = hit_mu(map_multi, 'paf', int(options.hitgap), reference_len, doubleQuery_len)

        # Build hits dict per chromosome
        # unique_hits[id] = [Qid, Qlen, Qstart, Qstop, strand, Tid, Tlen, Tstart, Tstop, matches, hitLen, quality]
        hits = {}
        for uid in list(unique_hits.keys()):
            Qid, Qlen, Qstart, Qstop, strand, Tid, Tlen, Tstart, Tstop, matches, hitLen, quality = unique_hits[uid]
            if Tid not in hits:
                hits[Tid] = []
            hits[Tid].append([Qid, Tstart, Tstop, Qstart, Qstop, matches, hitLen])

        print('[' + str(datetime.datetime.now()) + '] = Generating graphs and best paths', file=sys.stdout)

        used_by_chr_hap1 = {}

        # --- Hap1 (reference-only) ---
        print('[' + str(datetime.datetime.now()) + '] = Hap1', file=sys.stdout)
        for Target_sequence in sorted(hits.keys()):
            if Target_sequence in best_1_paths_edges:
                print('[' + str(datetime.datetime.now()) + '] == ' + Target_sequence + ' pre-loaded via --path1, skipping DAG', file=sys.stdout)
                continue
            print('[' + str(datetime.datetime.now()) + '] == Processing ' + Target_sequence, file=sys.stdout)
            hits_1 = list(hits[Target_sequence])
            hits_1.append(['ChrStart', 0, 0, 0, 0, 0, 0])
            hits_1.append(['ChrStop', reference_len[Target_sequence], reference_len[Target_sequence], 0, 0, 0, 0])

            preferred_db = {}
            distance1 = int(options.distance1)
            unwanted_to_remove_list = list(blacklist_1[Target_sequence])
            unwanted_to_remove_list = list(set(unwanted_to_remove_list))
            unwanted_to_reuse_list = []

            if Target_sequence in forced_list_1 and forced_list_1[Target_sequence]:
                new_forced = forced_list_1[Target_sequence][:]
            else:
                new_forced = []

            best_found = False
            while not best_found:
                if new_forced:
                    hit_graph = make_forced_graph(hits_1, distance1, new_forced, unwanted_to_remove_list, 'required.err.txt')
                else:
                    hit_graph = make_graph(hits_1, distance1, unwanted_to_remove_list)

                while not nx.has_path(hit_graph, source=0, target=reference_len[Target_sequence]):
                    distance1 += 100000
                    print('[' + str(datetime.datetime.now()) + '] ==== Increased allowed gap to ' + str(distance1), file=sys.stdout)
                    if new_forced:
                        hit_graph = make_forced_graph(hits_1, distance1, new_forced, unwanted_to_remove_list, 'required.err.txt')
                    else:
                        hit_graph = make_graph(hits_1, distance1, unwanted_to_remove_list)

                print('[' + str(datetime.datetime.now()) + '] === Searching 1st path maximizing coverage: ' +
                      Target_sequence + ' (0:' + str(reference_len[Target_sequence]) + ')', file=sys.stdout)
                best_1_nodes = nx.dag_longest_path(hit_graph, weight='align')
                best_1_edges, best_1_edges_names = get_subgraph_from_path(hit_graph, best_1_nodes)

                unwanted_to_remove = search_unwanted(best_1_edges_names, unwanted_pairs, all_seq_length,
                                                     forced_list_1[Target_sequence])

                if unwanted_to_remove == '':
                    all_preferred_used = True
                    to_reintegrate = []
                    preferred_to_remove = []
                    for seq_id in list(preferred_db.keys()):
                        if seq_id not in best_1_edges_names:
                            all_preferred_used = False
                            unwanted_to_reuse_list += three_orientation_list(seq_id, True)
                            preferred_to_remove.append(seq_id)
                            to_reintegrate += preferred_db[seq_id]
                    if all_preferred_used:
                        best_found = True
                    else:
                        new_forced = best_1_edges_names
                        for reintegrate_id in list(set(to_reintegrate)):
                            if reintegrate_id in unwanted_to_remove_list and reintegrate_id not in unwanted_to_reuse_list:
                                unwanted_to_remove_list.remove(reintegrate_id)
                        for seq_id in preferred_to_remove:
                            del preferred_db[seq_id]
                            unwanted_to_remove_list += three_orientation_list(seq_id, True)
                        unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                else:
                    updated_forced = [e for e in new_forced if e not in unwanted_to_remove['remove']]
                    new_forced = updated_forced
                    to_reintegrate = []
                    for key in unwanted_to_remove['remove']:
                        if key in preferred_db:
                            to_reintegrate += preferred_db[key]
                            del preferred_db[key]
                    if to_reintegrate:
                        unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                        for reintegrate_id in to_reintegrate:
                            if reintegrate_id in unwanted_to_remove_list and reintegrate_id not in unwanted_to_reuse_list:
                                unwanted_to_remove_list.remove(reintegrate_id)
                    unwanted_to_remove_list += unwanted_to_remove['remove']
                    unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                    if unwanted_to_remove['keep'] not in preferred_db:
                        preferred_db[unwanted_to_remove['keep']] = []
                    preferred_db[unwanted_to_remove['keep']] += unwanted_to_remove['remove']

            best_1_paths[Target_sequence] = best_1_edges
            used, best_1_paths_edges[Target_sequence] = make_list_from_path(best_1_edges)
            all_used += used
            used_by_chr_hap1[Target_sequence] = used
            print('[' + str(datetime.datetime.now()) + '] === Used for ' + Target_sequence + ': ' +
                  ','.join(best_1_paths_edges[Target_sequence]), file=sys.stdout)

        # --- Hap2 (reference-only) ---
        if not options.No2:
            paired_to_used = []
            if not options.rearrangements:
                for Target_sequence in list(best_1_paths_edges.keys()):
                    for seq_id in best_1_paths_edges[Target_sequence]:
                        if seq_id in known_groups_by_seqid:
                            paired_to_used += known_groups_by_seqid[seq_id]
                paired_to_used = list(set(paired_to_used))

            unwanted_sequences_alternative_to_used = {}
            for Target_sequence in list(best_1_paths_edges.keys()):
                unwanted_sequences_alternative_to_used[Target_sequence] = []
                for seq_id in best_1_paths_edges[Target_sequence]:
                    if seq_id in alternative_sequences:
                        for Target_sequence_2 in list(best_1_paths_edges.keys()):
                            if Target_sequence_2 != Target_sequence:
                                unwanted_sequences_alternative_to_used.setdefault(Target_sequence_2, [])
                                unwanted_sequences_alternative_to_used[Target_sequence_2] += alternative_sequences[seq_id]

            print('[' + str(datetime.datetime.now()) + '] = Hap2', file=sys.stdout)
            for Target_sequence in sorted(hits.keys()):
                if Target_sequence in best_2_paths_edges:
                    print('[' + str(datetime.datetime.now()) + '] == ' + Target_sequence + ' pre-loaded via --path2, skipping DAG', file=sys.stdout)
                    continue
                print('[' + str(datetime.datetime.now()) + '] == Processing ' + Target_sequence, file=sys.stdout)
                # Exclude sequences already used in Hap1
                hits_2 = [x for x in hits[Target_sequence] if x[0] not in all_used]
                hits_2.append(['ChrStart', 0, 0, 0, 0, 0, 0])
                hits_2.append(['ChrStop', reference_len[Target_sequence], reference_len[Target_sequence], 0, 0, 0, 0])

                preferred_db = {}
                distance2 = int(options.distance2)
                unwanted_to_remove_list = list(blacklist_2[Target_sequence])
                unwanted_to_remove_list += paired_to_used
                unwanted_to_remove_list += unwanted_sequences_alternative_to_used.get(Target_sequence, [])
                for seq_id in best_1_paths_edges.get(Target_sequence, []):
                    unwanted_to_remove_list += three_orientation_list(seq_id, True)
                    if seq_id in known_groups_by_seqid:
                        unwanted_to_remove_list += known_groups_by_seqid[seq_id]
                unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                unwanted_to_reuse_list = []

                if Target_sequence in forced_list_2 and forced_list_2[Target_sequence]:
                    new_forced = forced_list_2[Target_sequence][:]
                else:
                    new_forced = []

                best_found = False
                while not best_found:
                    if new_forced:
                        hit_graph_2 = make_forced_graph(hits_2, distance2, new_forced, unwanted_to_remove_list,
                                                        'required2.err.txt')
                    else:
                        hit_graph_2 = make_graph(hits_2, distance2, unwanted_to_remove_list)

                    while not nx.has_path(hit_graph_2, source=0, target=reference_len[Target_sequence]):
                        distance2 += 1000000
                        print('[' + str(datetime.datetime.now()) + '] ==== Increased allowed gap to ' + str(distance2), file=sys.stdout)
                        if new_forced:
                            hit_graph_2 = make_forced_graph(hits_2, distance2, new_forced, unwanted_to_remove_list,
                                                            'required2.err.txt')
                        else:
                            hit_graph_2 = make_graph(hits_2, distance2, unwanted_to_remove_list)

                    print('[' + str(datetime.datetime.now()) + '] === Searching 2nd path maximizing coverage: ' +
                          Target_sequence + ' (0:' + str(reference_len[Target_sequence]) + ')', file=sys.stdout)
                    best_2_nodes = nx.dag_longest_path(hit_graph_2, weight='align')
                    best_2_edges, best_2_edges_names = get_subgraph_from_path(hit_graph_2, best_2_nodes)

                    unwanted_to_remove = search_unwanted(best_2_edges_names, unwanted_pairs, all_seq_length,
                                                         forced_list_2[Target_sequence])

                    if unwanted_to_remove == '':
                        all_preferred_used = True
                        to_reintegrate = []
                        preferred_to_remove = []
                        for seq_id in list(preferred_db.keys()):
                            if seq_id not in best_2_edges_names:
                                all_preferred_used = False
                                unwanted_to_reuse_list += three_orientation_list(seq_id, True)
                                preferred_to_remove.append(seq_id)
                                to_reintegrate += preferred_db[seq_id]
                        if all_preferred_used:
                            best_found = True
                        else:
                            new_forced = best_2_edges_names
                            for reintegrate_id in list(set(to_reintegrate)):
                                if reintegrate_id in unwanted_to_remove_list and reintegrate_id not in unwanted_to_reuse_list:
                                    unwanted_to_remove_list.remove(reintegrate_id)
                            for seq_id in preferred_to_remove:
                                del preferred_db[seq_id]
                                unwanted_to_remove_list += three_orientation_list(seq_id, True)
                            unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                    else:
                        updated_forced = [e for e in new_forced if e not in unwanted_to_remove['remove']]
                        new_forced = updated_forced
                        to_reintegrate = []
                        for key in unwanted_to_remove['remove']:
                            if key in preferred_db:
                                to_reintegrate += preferred_db[key]
                                del preferred_db[key]
                        if to_reintegrate:
                            unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                            for reintegrate_id in to_reintegrate:
                                if reintegrate_id in unwanted_to_remove_list and reintegrate_id not in unwanted_to_reuse_list:
                                    unwanted_to_remove_list.remove(reintegrate_id)
                        unwanted_to_remove_list += unwanted_to_remove['remove']
                        unwanted_to_remove_list = list(set(unwanted_to_remove_list))
                        if unwanted_to_remove['keep'] not in preferred_db:
                            preferred_db[unwanted_to_remove['keep']] = []
                        preferred_db[unwanted_to_remove['keep']] += unwanted_to_remove['remove']

                best_2_paths[Target_sequence] = best_2_edges
                used, best_2_paths_edges[Target_sequence] = make_list_from_path(best_2_edges)
                all_used += used

                # Track unused: sequences with hits on this chr not placed in either hap
                all_used_base = set(s[:-2] for s in all_used if len(s) > 2)
                unused_by_chr[Target_sequence] = [
                    x[0] for x in hits[Target_sequence]
                    if x[0][:-2] not in all_used_base
                ]

    # =========================================================================
    # Fill missing orientation
    # =========================================================================
    print('[' + str(datetime.datetime.now()) + '] = Filling missing orientation info', file=sys.stdout)

    if use_reference and use_markers:
        # Both mode: use nucmer to resolve |. sequences
        print('[' + str(datetime.datetime.now()) + '] = Updating sequence orientation using guide genome', file=sys.stdout)
        tmp_dir = options.out + '.tmp'
        mkdir(tmp_dir)
        orienation_db = {}
        fasta_files = {}
        fasta_len_dict = {}

        # Write reference chromosome FASTA files
        for ref_id in list(reference.keys()):
            ref_fasta = reference[ref_id]
            ref_file_name = tmp_dir + '/' + ref_id + '.fasta'
            fasta_files[ref_id] = ref_file_name
            with open(ref_file_name, 'w') as ref_file:
                print('>' + ref_id, file=ref_file)
                print(ref_fasta, file=ref_file)
            fasta_len_dict[ref_id] = len(ref_fasta)

        # Find undirected sequences (|.)
        undirected_sequences = get_no_orientation(best_1_paths_edges, best_2_paths_edges)

        # Write query FASTA files for undirected sequences (fwd + rev)
        for query_id in list(undirected_sequences.keys()):
            query_fasta_seq = query_fasta_db[query_id]
            query_file_name = tmp_dir + '/' + query_id + '.fasta'
            fasta_files[query_id] = query_file_name
            with open(query_file_name, 'w') as query_file:
                print('>' + query_id + '|+', file=query_file)
                print(query_fasta_seq.upper(), file=query_file)
                print('>' + query_id + '|-', file=query_file)
                print(str(Seq(query_fasta_seq).reverse_complement()).upper(), file=query_file)
            fasta_len_dict[query_id + '|+'] = len(query_fasta_seq)
            fasta_len_dict[query_id + '|-'] = len(query_fasta_seq)

        nucmer_path = paths.get('nucmer', '')
        showcoords_path = paths.get('show-coords', '')

        # Map each undirected sequence on its chromosome
        for query_id in list(undirected_sequences.keys()):
            target_id = undirected_sequences[query_id]
            map_prefix = target_id + '.on.' + query_id
            coords_file = tmp_dir + '/' + map_prefix + '.coords'
            if options.reuse_intermediate and os.path.exists(coords_file):
                print('[' + str(datetime.datetime.now()) + '] === Reusing ' + coords_file, file=sys.stdout)
            else:
                coords_file = map_nucmer(fasta_files[target_id], fasta_files[query_id],
                                         int(options.cores), coords_file,
                                         nucmer_path, showcoords_path,
                                         ' --forward ', ' -l -r -T -H ')
            map_results = read_nucmer_coords(coords_file)
            best_alignment = hits_best_tiling_path(map_results, fasta_len_dict)
            orienation_db[query_id] = best_orientation(best_alignment)

        best_1_paths, best_1_paths_edges = fill_orientation(
            best_1_paths, best_1_paths_edges,
            options.out + '.missing_orientation.hap1.txt', orienation_db)
        if not options.No2:
            best_2_paths, best_2_paths_edges = fill_orientation(
                best_2_paths, best_2_paths_edges,
                options.out + '.missing_orientation.hap2.txt', orienation_db)

    elif use_markers:
        # Markers-only: default undirected to '+'
        best_1_paths, best_1_paths_edges = fill_orientation(
            best_1_paths, best_1_paths_edges, options.out + '.missing_orientation.hap1.txt')
        if not options.No2:
            best_2_paths, best_2_paths_edges = fill_orientation(
                best_2_paths, best_2_paths_edges, options.out + '.missing_orientation.hap2.txt')
    # else: reference-only — orientation already determined by alignment strand

    # =========================================================================
    # Write output list files
    # =========================================================================
    print('[' + str(datetime.datetime.now()) + '] = Writing path list files', file=sys.stdout)

    with open(options.out + '.1.list', 'w') as f:
        for chr_id in sorted(best_1_paths_edges.keys()):
            print(chr_id + '\t' + ','.join(best_1_paths_edges[chr_id]), file=f)

    if not options.No2:
        with open(options.out + '.2.list', 'w') as f:
            for chr_id in sorted(best_2_paths_edges.keys()):
                print(chr_id + '\t' + ','.join(best_2_paths_edges[chr_id]), file=f)

    all_unused = [id_ for id_ in sorted(query_len.keys()) if (id_ + '|+') not in all_used]
    with open(options.out + '.Un.list', 'w') as f:
        print(','.join(all_unused), file=f)

    with open(options.out + '.unused_sequences.list', 'w') as f:
        for chr_id in sorted(unused_by_chr.keys()):
            print(chr_id + '\t' + ','.join(unused_by_chr[chr_id]), file=f)

    print('[' + str(datetime.datetime.now()) + '] = Done', file=sys.stdout)


if __name__ == '__main__':
    main()
