#!/usr/bin/env python3
"""
reconstruct.py — AGP + FASTA + correspondence reconstruction.

Reads the tiling-path list files produced by build_paths.py and writes:
  {out}.1.agp / .2.agp / .Un.agp
  {out}.1.fasta / .2.fasta / .Un.fasta
  {out}.correspondence.tsv
  {out}.conflicting_pseudomolecules.txt
  {out}.unplaced_to_pseudomolecule.txt
"""

import argparse
import os
import sys
import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *
from lib_files.GFF_lib import *


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__)

    # Required
    parser.add_argument('-i', '--input', dest='query', required=True,
        metavar='query.fasta',
        help='FASTA file of the input sequences')
    parser.add_argument('-o', '--out', dest='out', default='out',
        metavar='NAME',
        help='Output files prefix (must match build_paths -o) [default: out]')
    parser.add_argument('-p', '--prefix', dest='prefix', default='NEW',
        metavar='PREFIX',
        help='Prefix for output sequence IDs [default: NEW]')

    # Optional inputs
    parser.add_argument('-k', '--known', dest='known',
        metavar='known.tsv',
        help='Sequences known to be in the same haplotype (for relationship reporting)')
    parser.add_argument('--alternative_groups', dest='alternative',
        metavar='alternative_groups.tsv',
        help='Alternative haplotype sequence pairs (for relationship reporting)')

    # Assembly options
    parser.add_argument('--gapsize', dest='gap', default='1000',
        metavar='N',
        help='Gap size placeholder for AGP and FASTA (bp) [default: 1000]')
    parser.add_argument('--concatenate', dest='conc', default='',
        metavar='N',
        help='Concatenate unplaced sequences using this gap size (bp) [default: separate]')
    parser.add_argument('--N2', dest='No2',
        default=False, action='store_true',
        help="Don't reconstruct the second haplotype")

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    options = parser.parse_args()

    print("Running reconstruct.py from HaploSync version " + get_version(), file=sys.stdout)
    print("Command: " + " ".join(sys.argv), file=sys.stdout)
    print("----", file=sys.stdout)

    gap_length = int(options.gap)

    # -------------------------------------------------------------------------
    # Load query FASTA
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Reading query sequences', file=sys.stdout)
    query_fasta_db = read_fasta(options.query)
    query_len = {}
    for name in query_fasta_db:
        query_len[name] = int(len(query_fasta_db[name]))
    print('[' + str(datetime.datetime.now()) + '] === ' + str(len(query_len)) + ' query sequences', file=sys.stdout)

    # -------------------------------------------------------------------------
    # Load optional constraint files (for relationship reporting)
    # -------------------------------------------------------------------------
    known_groups = {}
    if options.known:
        print('[' + str(datetime.datetime.now()) + '] = Reading known haplotype groups', file=sys.stdout)
        for line in open(options.known):
            line = line.rstrip()
            if not line or line[0] == '#':
                continue
            try:
                seq_id, group_id = line.split('\t')
            except:
                print('[ERROR] Known groups format error: ' + line, file=sys.stderr)
                sys.exit(1)
            if group_id not in known_groups:
                known_groups[group_id] = []
            known_groups[group_id].append(seq_id)

    alternative_sequences = {}
    if options.alternative:
        for line in open(options.alternative):
            line = line.rstrip()
            if not line or line[0] == '#':
                continue
            group_1, group_2 = line.split('\t')
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

    # -------------------------------------------------------------------------
    # Read list files from build_paths
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Reading tiling path lists', file=sys.stdout)

    best_1_paths_edges = {}
    list_1_path = options.out + '.1.list'
    for line in open(list_1_path):
        line = line.rstrip('\n').rstrip('\r')
        if not line:
            continue
        parts = line.split('\t', 1)
        chr_id = parts[0]
        path_str = parts[1] if len(parts) > 1 else ''
        if path_str:
            best_1_paths_edges[chr_id] = path_str.split(',')

    best_2_paths_edges = {}
    if not options.No2:
        list_2_path = options.out + '.2.list'
        for line in open(list_2_path):
            line = line.rstrip('\n').rstrip('\r')
            if not line:
                continue
            parts = line.split('\t', 1)
            chr_id = parts[0]
            path_str = parts[1] if len(parts) > 1 else ''
            if path_str:
                best_2_paths_edges[chr_id] = path_str.split(',')

    all_unused = []
    un_list_path = options.out + '.Un.list'
    line = open(un_list_path).read().rstrip()
    if line:
        all_unused = line.split(',')

    unused_by_chr = {}
    rejected_path = options.out + '.unused_sequences.list'
    for line in open(rejected_path):
        line = line.rstrip('\n').rstrip('\r')
        if not line:
            continue
        parts = line.split('\t', 1)
        chr_id = parts[0]
        seq_str = parts[1] if len(parts) > 1 else ''
        unused_by_chr[chr_id] = seq_str.split(',') if seq_str else []

    print('[' + str(datetime.datetime.now()) + '] === Hap1: ' + str(len(best_1_paths_edges)) + ' chromosomes', file=sys.stdout)
    if not options.No2:
        print('[' + str(datetime.datetime.now()) + '] === Hap2: ' + str(len(best_2_paths_edges)) + ' chromosomes', file=sys.stdout)
    print('[' + str(datetime.datetime.now()) + '] === Unplaced: ' + str(len(all_unused)) + ' sequences', file=sys.stdout)

    # -------------------------------------------------------------------------
    # Write AGP files
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Writing AGP files (gap size ' + str(gap_length) + ' bp)', file=sys.stdout)

    agp_1_file_name  = options.out + '.1.agp'
    agp_2_file_name  = options.out + '.2.agp'
    agp_Un_file_name = options.out + '.Un.agp'

    # Mappings built during AGP writing — reused for FASTA, relationships, correspondence
    ref_to_hap1  = {}
    hap1_to_ref  = {}
    ref_to_hap2  = {}
    hap2_to_ref  = {}
    hap1_to_hap2 = {}
    hap2_to_hap1 = {}

    print('[' + str(datetime.datetime.now()) + '] == Hap1', file=sys.stdout)
    agp_1_file = open(agp_1_file_name, 'w')
    for Target_sequence in sorted(best_1_paths_edges.keys()):
        Obj_name = options.prefix + '_Hap1_' + Target_sequence
        ref_to_hap1[Target_sequence] = Obj_name
        hap1_to_ref[Obj_name] = Target_sequence
        print('[' + str(datetime.datetime.now()) + '] === ' + Obj_name, file=sys.stdout)
        make_agp_from_list(best_1_paths_edges[Target_sequence], query_len, gap_length, Obj_name, agp_1_file)
    agp_1_file.close()

    print('[' + str(datetime.datetime.now()) + '] == Hap2', file=sys.stdout)
    if not options.No2:
        agp_2_file = open(agp_2_file_name, 'w')
        for Target_sequence in sorted(best_2_paths_edges.keys()):
            Obj_name = options.prefix + '_Hap2_' + Target_sequence
            ref_to_hap2[Target_sequence] = Obj_name
            hap2_to_ref[Obj_name] = Target_sequence
            hap1_id = ref_to_hap1.get(Target_sequence)
            if hap1_id:
                hap2_to_hap1[Obj_name] = hap1_id
                hap1_to_hap2[hap1_id] = Obj_name
            print('[' + str(datetime.datetime.now()) + '] === ' + Obj_name, file=sys.stdout)
            make_agp_from_list(best_2_paths_edges[Target_sequence], query_len, gap_length, Obj_name, agp_2_file)
        agp_2_file.close()

    print('[' + str(datetime.datetime.now()) + '] == Unplaced', file=sys.stdout)
    agp_Un_file = open(agp_Un_file_name, 'w')
    if options.conc != '':
        conc_gap = int(options.conc)
        Obj_name = options.prefix + '_Un'
        make_agp_from_list(all_unused, query_len, conc_gap, Obj_name, agp_Un_file)
    else:
        for idx, comp in enumerate(all_unused, 1):
            Obj_name = options.prefix + '_Un_' + str(idx)
            make_agp_from_list([comp + '|+'], query_len, 0, Obj_name, agp_Un_file)
    agp_Un_file.close()

    # -------------------------------------------------------------------------
    # Build and write FASTA sequences
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Building FASTA sequences (gap size ' + str(gap_length) + ' bp)', file=sys.stdout)

    fasta_to_chr   = {}
    chr_to_fasta_1 = {}
    chr_to_fasta_2 = {}

    fasta_db_1 = {}
    print('[' + str(datetime.datetime.now()) + '] == Hap1', file=sys.stdout)
    for Target_sequence in sorted(best_1_paths_edges.keys()):
        Obj_name = options.prefix + '_Hap1_' + Target_sequence
        print('[' + str(datetime.datetime.now()) + '] === ' + Obj_name, file=sys.stdout)
        fasta_db_1[Obj_name] = make_fasta_from_list(best_1_paths_edges[Target_sequence], query_fasta_db, gap_length, Obj_name)
        fasta_to_chr[Obj_name] = Target_sequence
        chr_to_fasta_1[Target_sequence] = Obj_name

    fasta_db_2 = {}
    if not options.No2:
        print('[' + str(datetime.datetime.now()) + '] == Hap2', file=sys.stdout)
        for Target_sequence in sorted(best_2_paths_edges.keys()):
            Obj_name = options.prefix + '_Hap2_' + Target_sequence
            print('[' + str(datetime.datetime.now()) + '] === ' + Obj_name, file=sys.stdout)
            fasta_db_2[Obj_name] = make_fasta_from_list(best_2_paths_edges[Target_sequence], query_fasta_db, gap_length, Obj_name)
            fasta_to_chr[Obj_name] = Target_sequence
            chr_to_fasta_2[Target_sequence] = Obj_name

    fasta_db_un = {}
    if options.conc != '':
        conc_gap = int(options.conc)
        Obj_name = options.prefix + '_Un'
        fasta_db_un[Obj_name] = make_fasta_from_list(all_unused, query_fasta_db, conc_gap, Obj_name)
    else:
        for idx, comp in enumerate(all_unused, 1):
            Obj_name = options.prefix + '_Un_' + str(idx)
            fasta_db_un[Obj_name] = make_fasta_from_list([comp + '|+'], query_fasta_db, 0, Obj_name)

    print('[' + str(datetime.datetime.now()) + '] = Writing FASTA files', file=sys.stdout)
    write_fasta_from_db(fasta_db_1,  options.out + '.1.fasta')
    if not options.No2:
        write_fasta_from_db(fasta_db_2,  options.out + '.2.fasta')
    write_fasta_from_db(fasta_db_un, options.out + '.Un.fasta')

    # -------------------------------------------------------------------------
    # Relationship reporting
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Building relationship reports', file=sys.stdout)

    groups_by_seqid = {}
    for group_id in known_groups:
        for seq_id in known_groups[group_id]:
            if seq_id not in groups_by_seqid:
                groups_by_seqid[seq_id] = []
            groups_by_seqid[seq_id].append(group_id)

    groups_by_pseudomolecule  = {}
    pseudomolecule_by_group   = {}
    alternative_by_pseudomolecule = {}
    pseudomolecule_by_alternative = {}

    for Target_sequence in sorted(best_1_paths_edges.keys()):
        new_id = options.prefix + '_Hap1_' + Target_sequence
        groups_by_pseudomolecule[new_id] = []
        for seq_id in best_1_paths_edges[Target_sequence]:
            clean_name = seq_id[:-2]
            if clean_name in groups_by_seqid:
                seq_id_groups = groups_by_seqid[clean_name]
                groups_by_pseudomolecule[new_id] += seq_id_groups
                for group_id in seq_id_groups:
                    if group_id not in pseudomolecule_by_group:
                        pseudomolecule_by_group[group_id] = []
                    pseudomolecule_by_group[group_id].append(new_id)
            if clean_name in alternative_sequences:
                alt_ids = alternative_sequences[clean_name]
                alternative_by_pseudomolecule[new_id] = list(set(
                    alternative_by_pseudomolecule.get(new_id, []) + alt_ids))
                for alt_id in alt_ids:
                    if alt_id not in pseudomolecule_by_alternative:
                        pseudomolecule_by_alternative[alt_id] = []
                    pseudomolecule_by_alternative[alt_id].append(new_id)
                    pseudomolecule_by_alternative[alt_id] = list(set(pseudomolecule_by_alternative[alt_id]))

    for Target_sequence in sorted(best_2_paths_edges.keys()):
        new_id = options.prefix + '_Hap2_' + Target_sequence
        groups_by_pseudomolecule[new_id] = []
        for seq_id in best_2_paths_edges[Target_sequence]:
            clean_name = seq_id[:-2]
            if clean_name in groups_by_seqid:
                seq_id_groups = groups_by_seqid[clean_name]
                groups_by_pseudomolecule[new_id] += seq_id_groups
                for group_id in seq_id_groups:
                    if group_id not in pseudomolecule_by_group:
                        pseudomolecule_by_group[group_id] = []
                    pseudomolecule_by_group[group_id].append(new_id)
            if clean_name in alternative_sequences:
                alt_ids = alternative_sequences[clean_name]
                alternative_by_pseudomolecule[new_id] = list(set(
                    alternative_by_pseudomolecule.get(new_id, []) + alt_ids))
                for alt_id in alt_ids:
                    if alt_id not in pseudomolecule_by_alternative:
                        pseudomolecule_by_alternative[alt_id] = []
                    pseudomolecule_by_alternative[alt_id].append(new_id)
                    pseudomolecule_by_alternative[alt_id] = list(set(pseudomolecule_by_alternative[alt_id]))

    # Pseudomolecule vs pseudomolecule conflicts
    conflicting_pm_file = open(options.out + '.conflicting_pseudomolecules.txt', 'w')
    print('Chr_id\tConflicting_pseudomolecules', file=conflicting_pm_file)
    for pm_id in sorted(groups_by_pseudomolecule.keys()):
        pm_groups = groups_by_pseudomolecule[pm_id]
        conflicting = []
        for group_id in pm_groups:
            conflicting += pseudomolecule_by_group[group_id]
        conflicting = list(set(conflicting))
        if pm_id in conflicting:
            conflicting.remove(pm_id)
        print(pm_id + '\t' + ','.join(conflicting), file=conflicting_pm_file)
    conflicting_pm_file.close()

    # Unplaced sequence to pseudomolecule relationships
    unused_by_seq_id     = {}
    all_unused_by_seq_id = {}
    for chr_id in unused_by_chr:
        for seq_id_oriented in unused_by_chr[chr_id]:
            seq_id, orientation = seq_id_oriented.split('|')
            unused_by_seq_id[seq_id]     = [chr_id, orientation]
            all_unused_by_seq_id[seq_id] = [chr_id, orientation]

    conflicting_un_file = open(options.out + '.unplaced_to_pseudomolecule.txt', 'w')
    print('Seq_id\tConflicting_marker_Vs_associated\tMultiple_associated_pseudomolecules\t'
          'Conflicting_marker_Vs_alternative\tMultiple_alternative_pseudomolecules\t'
          'Expected_chr\tOrientation\tAssociated_pseudomolecules\tAlternative_pseudomolecules',
          file=conflicting_un_file)

    for seq_id in all_unused:
        if seq_id in unused_by_seq_id:
            chr_id, orientation = unused_by_seq_id[seq_id]
        else:
            chr_id, orientation = '.', '.'

        associated_pseudomolecules = []
        if seq_id in groups_by_seqid:
            for group_id in groups_by_seqid[seq_id]:
                if group_id in pseudomolecule_by_group:
                    associated_pseudomolecules += pseudomolecule_by_group[group_id]
            associated_pseudomolecules = list(set(associated_pseudomolecules))

        if len(associated_pseudomolecules) == 0:
            err_chr_assoc = '.'
            mult_assoc = '.'
            assoc_chromosomes = [['.']  ]
        elif len(associated_pseudomolecules) > 1:
            mult_assoc = 'TRUE'
            assoc_chromosomes = list(set([fasta_to_chr[x] for x in associated_pseudomolecules]))
            if len(assoc_chromosomes) > 1:
                err_chr_assoc = 'TRUE'
            else:
                err_chr_assoc = 'FALSE' if assoc_chromosomes[0] == chr_id else ('.' if chr_id == '.' else 'TRUE')
        else:
            mult_assoc = 'FALSE'
            assoc_chromosomes = [fasta_to_chr[associated_pseudomolecules[0]]]
            err_chr_assoc = 'FALSE' if assoc_chromosomes == [chr_id] else ('.' if chr_id == '.' else 'TRUE')

        alternative_pseudomolecules = list(set(pseudomolecule_by_alternative.get(seq_id, [])))

        if len(alternative_pseudomolecules) == 0:
            err_chr_alt = '.'
            mult_alt = '.'
            alt_chromosomes = [['.'] ]
        elif len(alternative_pseudomolecules) > 1:
            mult_alt = 'TRUE'
            alt_chromosomes = list(set([fasta_to_chr[x] for x in alternative_pseudomolecules]))
            if len(alt_chromosomes) > 1:
                err_chr_alt = 'TRUE'
            else:
                err_chr_alt = 'FALSE' if alt_chromosomes[0] == chr_id else ('.' if chr_id == '.' else 'TRUE')
        else:
            mult_alt = 'FALSE'
            alt_chromosomes = [fasta_to_chr[alternative_pseudomolecules[0]]]
            err_chr_alt = 'FALSE' if alt_chromosomes == [chr_id] else ('.' if chr_id == '.' else 'TRUE')

        print('\t'.join([
            seq_id, err_chr_assoc, mult_assoc, err_chr_alt, mult_alt, chr_id, orientation,
            ','.join(associated_pseudomolecules), ','.join(alternative_pseudomolecules)
        ]), file=conflicting_un_file)

        if seq_id not in all_unused_by_seq_id:
            if len(assoc_chromosomes) == 1 and len(alt_chromosomes) == 1:
                if assoc_chromosomes[0] == alt_chromosomes[0]:
                    if assoc_chromosomes[0] != ['.']:
                        all_unused_by_seq_id[seq_id] = [assoc_chromosomes[0], '+']
                else:
                    if assoc_chromosomes[0] != ['.'] and alt_chromosomes[0] == ['.']:
                        all_unused_by_seq_id[seq_id] = [assoc_chromosomes[0], '+']
                    elif assoc_chromosomes[0] == ['.'] and alt_chromosomes[0] != '.':
                        all_unused_by_seq_id[seq_id] = [alt_chromosomes[0], '+']

    conflicting_un_file.close()

    # -------------------------------------------------------------------------
    # Write correspondence file
    # -------------------------------------------------------------------------
    if not options.No2 and chr_to_fasta_1 and chr_to_fasta_2:
        corr_file_name = options.out + '.correspondence.tsv'
        corr_file = open(corr_file_name, 'w')
        for chr_id in sorted(chr_to_fasta_1.keys()):
            hap1_id = chr_to_fasta_1[chr_id]
            hap2_id = chr_to_fasta_2.get(chr_id)
            if hap2_id is None:
                continue
            print(chr_id + '\t' + hap1_id + '\t' + hap2_id, file=corr_file)
        corr_file.close()
        print('[' + str(datetime.datetime.now()) + '] = Written correspondence file: ' + corr_file_name, file=sys.stdout)


if __name__ == '__main__':
    main()
