#!/usr/bin/env python3
"""
rejected_qc.py — Generate per-unplaced-sequence QC reports.

Reads HaploSplit outputs (unused_sequences.list, correspondence, FASTA, AGP,
optional markers/groups) and generates one HTML report per unplaced sequence
via rejected_QC().

Outputs:
  {out}.structure_comparison/                              — per-sequence reports
  {out}.structure_comparison/index.rejected_sequences.html — index
"""

import argparse
import concurrent.futures
import datetime
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__)

    parser.add_argument('-u', '--unused', dest='unused', required=True,
        metavar='unused.list',
        help='Unused sequences list (chr_id TAB seq1|+,seq2|-,...)')
    parser.add_argument('-c', '--correspondence', dest='correspondence', required=True,
        metavar='correspondence.tsv',
        help='Correspondence file (chr_id TAB hap1_id TAB hap2_id)')
    parser.add_argument('-f', '--fasta', dest='fasta', required=True,
        metavar='hap1.fasta,hap2.fasta,un.fasta',
        help='Comma-separated FASTA files containing all sequences')
    parser.add_argument('-a', '--agp', dest='agp', required=True,
        metavar='combined.agp',
        help='Combined AGP structure of pseudomolecules')
    parser.add_argument('-b', '--markers_bed', dest='markers_bed', default=None,
        metavar='markers.bed',
        help='Translated marker positions (BED, 4 columns)')
    parser.add_argument('-m', '--markers_map', dest='marker_map', default=None,
        metavar='markers_map.tsv',
        help='Marker genetic map (chr_id, position, marker_id)')
    parser.add_argument('--legacy_agp', dest='legacy_agp', default=None,
        metavar='legacy.agp',
        help='Legacy AGP structure')
    parser.add_argument('--input_groups', dest='input_groups', default=None,
        metavar='input_groups.tsv',
        help='Input sequence group file (seq_id, group_id)')
    parser.add_argument('--legacy_groups', dest='legacy_groups', default=None,
        metavar='legacy_groups.tsv',
        help='Legacy component group file (seq_id, group_id)')
    parser.add_argument('-o', '--out', dest='out', default='out',
        metavar='PREFIX',
        help='Output prefix [default: out]')
    parser.add_argument('-t', '--threads', dest='cores', default=4, type=int,
        metavar='N',
        help='CPU cores [default: 4]')
    parser.add_argument('--conf', dest='conf', default=None,
        metavar='HaploSync.conf.toml',
        help='Path to HaploSync conf.toml [default: auto-detect]')

    options = parser.parse_args()

    conf_path = options.conf or os.path.join(
        os.path.dirname(os.path.realpath(__file__)), '..', 'HaploSync.conf.toml')
    paths = set_paths(conf_path)

    print('[' + str(datetime.datetime.now()) + '] = Running rejected_qc.py', file=sys.stdout)
    print('[' + str(datetime.datetime.now()) + '] = Command: ' + ' '.join(sys.argv), file=sys.stdout)

    # -------------------------------------------------------------------------
    # Pre-flight: verify Rscript and rmarkdown are available
    # -------------------------------------------------------------------------
    import shutil
    import subprocess as _sp
    _rscript_bin = shutil.which('Rscript') or 'Rscript'
    print('[rejected_qc] Rscript resolved to: ' + _rscript_bin, file=sys.stderr)
    try:
        _check = _sp.run([_rscript_bin, '-e', 'library(rmarkdown)'],
                         capture_output=True, text=True)
        if _check.returncode != 0:
            print('[rejected_qc] ERROR: rmarkdown check failed (exit ' + str(_check.returncode) + ')', file=sys.stderr)
            print('[rejected_qc]   stdout: ' + _check.stdout.strip(), file=sys.stderr)
            print('[rejected_qc]   stderr: ' + _check.stderr.strip(), file=sys.stderr)
            sys.exit(1)
        print('[rejected_qc] rmarkdown check passed.', file=sys.stderr)
    except FileNotFoundError:
        print('[rejected_qc] ERROR: Rscript not found at: ' + _rscript_bin, file=sys.stderr)
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Load FASTA
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] == Loading FASTA sequences', file=sys.stdout)
    fasta_dict = {}
    for fasta_file in options.fasta.split(','):
        fasta_file = fasta_file.strip()
        if not fasta_file:
            continue
        print('[' + str(datetime.datetime.now()) + '] === ' + fasta_file, file=sys.stdout)
        fasta_dict.update(read_fasta(fasta_file))
    fasta_len_dict = get_length_from_fasta_db(fasta_dict)

    # -------------------------------------------------------------------------
    # Load correspondence
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] == Loading correspondence', file=sys.stdout)
    chr_ids = []
    chr_to_hap1 = {}
    chr_to_hap2 = {}
    hap1_to_chr = {}
    hap2_to_chr = {}
    fasta_db_1 = {}
    fasta_db_2 = {}
    for line in open(options.correspondence):
        line = line.rstrip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        chr_id, hap1_id, hap2_id = parts[0], parts[1], parts[2]
        chr_ids.append(chr_id)
        chr_to_hap1[chr_id] = hap1_id
        chr_to_hap2[chr_id] = hap2_id
        hap1_to_chr[hap1_id] = chr_id
        hap2_to_chr[hap2_id] = chr_id
        if hap1_id in fasta_dict:
            fasta_db_1[hap1_id] = fasta_dict[hap1_id]
        if hap2_id in fasta_dict:
            fasta_db_2[hap2_id] = fasta_dict[hap2_id]

    # -------------------------------------------------------------------------
    # Load unused sequences list
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] == Loading unused sequences list', file=sys.stdout)
    all_unused_by_chr = {}
    all_unused_to_chr = {}
    for line in open(options.unused):
        line = line.rstrip('\n').rstrip('\r')
        if not line:
            continue
        parts = line.split('\t', 1)
        chr_id = parts[0]
        seq_str = parts[1] if len(parts) > 1 else ''
        for seq_entry in seq_str.split(','):
            seq_entry = seq_entry.strip()
            if not seq_entry:
                continue
            if chr_id not in all_unused_by_chr:
                all_unused_by_chr[chr_id] = []
            all_unused_by_chr[chr_id].append(seq_entry)
            seq_id_bare = seq_entry[:-2]
            all_unused_to_chr[seq_id_bare] = chr_id
            all_unused_to_chr[seq_entry] = chr_id

    # -------------------------------------------------------------------------
    # Load AGP
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] == Loading AGP', file=sys.stdout)
    agp_db = read_agp(options.agp)

    # -------------------------------------------------------------------------
    # Load legacy AGP
    # -------------------------------------------------------------------------
    if options.legacy_agp:
        print('[' + str(datetime.datetime.now()) + '] == Loading legacy AGP', file=sys.stdout)
        legacy_agp = read_agp(options.legacy_agp)
    else:
        legacy_agp = ""

    # -------------------------------------------------------------------------
    # Load markers BED
    # -------------------------------------------------------------------------
    if options.markers_bed:
        print('[' + str(datetime.datetime.now()) + '] == Loading markers BED', file=sys.stdout)
        markers_db = {}
        for line in open(options.markers_bed):
            line = line.rstrip()
            if not line:
                continue
            chr_id, start, stop, marker_id = line.split('\t')
            if chr_id not in markers_db:
                markers_db[chr_id] = {}
            if marker_id not in markers_db[chr_id]:
                markers_db[chr_id][marker_id] = []
            markers_db[chr_id][marker_id].append([chr_id, start, stop, marker_id])
    else:
        markers_db = ""

    # -------------------------------------------------------------------------
    # Load markers map
    # -------------------------------------------------------------------------
    if options.marker_map and options.markers_bed:
        print('[' + str(datetime.datetime.now()) + '] == Loading markers map', file=sys.stdout)
        marker_map_by_seq = {}
        marker_map_by_id = {}
        for line in open(options.marker_map):
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue
            try:
                seq_id, pos, marker_id = line.split('\t')
            except ValueError:
                continue
            if seq_id not in marker_map_by_seq:
                marker_map_by_seq[seq_id] = []
            marker_map_by_seq[seq_id].append([int(pos), marker_id])
            marker_map_by_id[marker_id] = [seq_id, int(pos), marker_id]
    else:
        marker_map_by_seq = ""
        marker_map_by_id = {}

    # -------------------------------------------------------------------------
    # Build output directory and association files
    # -------------------------------------------------------------------------
    structure_comparison_dir = options.out + '.structure_comparison'
    mkdir(structure_comparison_dir)

    hap1_ids_list = sorted(fasta_db_1.keys())
    hap2_ids_list = sorted(fasta_db_2.keys())
    all_seq_ids = list(fasta_dict.keys())

    if options.legacy_groups:
        print('[' + str(datetime.datetime.now()) + '] == Building legacy group associations', file=sys.stdout)
        associated_legacy_seqid_file = structure_comparison_dir + '/legacy_components.association.tsv'
        all_agp_db = dict(legacy_agp) if legacy_agp != "" else {}
        all_agp_db.update(agp_db)
        associated_input_seqid_file, seq_group_db = make_seq_pair_from_groups(
            associated_legacy_seqid_file, options.legacy_groups, all_agp_db,
            hap1_ids_list, hap2_ids_list, all_seq_ids, {}, 'legacy')
    else:
        associated_legacy_seqid_file = ""
        if options.input_groups:
            print('[' + str(datetime.datetime.now()) + '] == Building input group associations', file=sys.stdout)
            associated_input_seqid_file = structure_comparison_dir + '/input.association.tsv'
            associated_input_seqid_file, seq_group_db = make_seq_pair_from_groups(
                associated_input_seqid_file, options.input_groups, agp_db,
                hap1_ids_list, hap2_ids_list, all_seq_ids, fasta_len_dict, 'input')
        else:
            # Create an empty association file to avoid the NUCMER fallback in rejected_QC().
            # When both associated_input_seqid_file and associated_legacy_ids_file are "",
            # rejected_QC() falls into a nucmer alignment path that requires the nucmer binary.
            # Passing a non-empty path to an empty file makes it take path 2 instead,
            # producing plots without sequence-derived connection lines.
            associated_input_seqid_file = structure_comparison_dir + '/input.association.tsv'
            open(associated_input_seqid_file, 'w').close()
            seq_group_db = {}

    # -------------------------------------------------------------------------
    # Report marker usage
    # -------------------------------------------------------------------------
    if options.markers_bed:
        print('[' + str(datetime.datetime.now()) + '] == Reporting marker usage', file=sys.stdout)
        legacy_agp_for_marker = legacy_agp if legacy_agp != "" else {}
        reported_hits_on_seq = report_marker_usage(
            options.markers_bed, marker_map_by_seq, marker_map_by_id,
            agp_db, legacy_agp_for_marker,
            hap1_to_chr, hap2_to_chr, all_unused_to_chr, structure_comparison_dir)

        # Remap markers from wrapper IDs to bare contig IDs
        bare_to_wrapper = {}
        for wrapper_id in agp_db:
            for agp_start in sorted(agp_db[wrapper_id].keys()):
                entry = agp_db[wrapper_id][agp_start]
                if entry[4] == 'W':
                    bare_to_wrapper[entry[5]] = (wrapper_id, entry)
        for unpl_id in sorted(all_unused_to_chr.keys()):
            if '|' in unpl_id or unpl_id in markers_db:
                continue
            if unpl_id not in bare_to_wrapper:
                continue
            wrapper_id, entry = bare_to_wrapper[unpl_id]
            if wrapper_id in reported_hits_on_seq:
                reported_hits_on_seq[unpl_id] = reported_hits_on_seq[wrapper_id]
        # Also remap markers_db: translate.py keys the output BED by wrapper pseudomolecule
        # IDs (e.g. prefix_Un_1), but rejected_QC() looks up bare contig IDs (e.g. seq99).
        for unpl_id in sorted(all_unused_to_chr.keys()):
            if '|' in unpl_id or unpl_id in markers_db:
                continue
            if unpl_id not in bare_to_wrapper:
                continue
            wrapper_id, entry = bare_to_wrapper[unpl_id]
            if wrapper_id in markers_db:
                markers_db[unpl_id] = markers_db[wrapper_id]
    else:
        reported_hits_on_seq = ""

    # -------------------------------------------------------------------------
    # Run rejected QC
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Running unused sequences QC analysis', file=sys.stdout)
    structure_plot_db = {'Rejected': {}}

    work_items = []
    for chr_id in sorted(chr_ids):
        if chr_id not in all_unused_by_chr:
            continue
        qc_out_dir = structure_comparison_dir + '/' + chr_id
        mkdir(qc_out_dir)
        structure_plot_db['Rejected'][chr_id] = {}
        for seq_id in all_unused_by_chr[chr_id]:
            work_items.append((chr_id, seq_id))

    def _run_rejected_qc(chr_id, seq_id):
        print('### seq_id: ' + seq_id, file=sys.stderr)
        query_id_bare = seq_id[:-2]
        unpl_input_agp = {}
        wrapper_name = None
        for agp_key in agp_db:
            for agp_start in sorted(agp_db[agp_key].keys()):
                entry = agp_db[agp_key][agp_start]
                if entry[4] == 'W' and entry[5] == query_id_bare:
                    wrapper_name = agp_key
                    break
            if wrapper_name is not None:
                break
        if wrapper_name is not None and legacy_agp != "" and wrapper_name in legacy_agp:
            for leg_start in sorted(legacy_agp[wrapper_name].keys()):
                entry = legacy_agp[wrapper_name][leg_start]
                Obj_Name_e, Obj_start_e, Obj_End_e, PartNum_e, Compnt_Type_e, CompntId_e, CompntStart_e, CompntEnd_e, Orientation_e = entry
                if Compnt_Type_e == 'W':
                    if query_id_bare not in unpl_input_agp:
                        unpl_input_agp[query_id_bare] = {}
                    unpl_input_agp[query_id_bare][int(Obj_start_e)] = [
                        query_id_bare, Obj_start_e, Obj_End_e, PartNum_e, 'W',
                        CompntId_e, CompntStart_e, CompntEnd_e, Orientation_e]
        elif wrapper_name is not None:
            for agp_start in sorted(agp_db[wrapper_name].keys()):
                entry = agp_db[wrapper_name][agp_start]
                Obj_Name_e, Obj_start_e, Obj_End_e, PartNum_e, Compnt_Type_e, CompntId_e, CompntStart_e, CompntEnd_e, Orientation_e = entry
                if Compnt_Type_e == 'W' and CompntId_e == query_id_bare:
                    if query_id_bare not in unpl_input_agp:
                        unpl_input_agp[query_id_bare] = {}
                    unpl_input_agp[query_id_bare][int(CompntStart_e)] = [
                        query_id_bare, CompntStart_e, CompntEnd_e, PartNum_e, 'W',
                        CompntId_e, CompntStart_e, CompntEnd_e, Orientation_e]
        legacy_agp_arg = legacy_agp if legacy_agp != "" else ""
        return chr_id, seq_id, rejected_QC(
            structure_comparison_dir, seq_id, fasta_dict, chr_id,
            fasta_db_1, fasta_db_2, chr_to_hap1, chr_to_hap2,
            '', associated_input_seqid_file, associated_legacy_seqid_file,
            agp_db, legacy_agp_arg, unpl_input_agp,
            seq_group_db, markers_db, reported_hits_on_seq, marker_map_by_seq,
            options.cores, paths)

    with concurrent.futures.ThreadPoolExecutor(max_workers=options.cores) as executor:
        futures = {
            executor.submit(_run_rejected_qc, chr_id, seq_id): (chr_id, seq_id)
            for chr_id, seq_id in work_items
        }
        for future in concurrent.futures.as_completed(futures):
            chr_id, seq_id, result = future.result()
            structure_plot_db['Rejected'][chr_id][seq_id] = result
            # Forward Rscript errors to stderr so they appear in Nextflow logs
            query_id_bare = seq_id[:-2]
            err_file = structure_comparison_dir + '/.' + query_id_bare + '_qc.err'
            if os.path.exists(err_file):
                err_content = open(err_file).read().strip()
                if err_content:
                    print('[rejected_qc] Rscript stderr for ' + seq_id + ':\n' + err_content,
                          file=sys.stderr)

    rejected_index_file_full_path = structure_comparison_dir + '/index.rejected_sequences.html'
    make_index_from_report_db('index.rejected_sequences.html', '.', structure_comparison_dir, structure_plot_db)
    print('[' + str(datetime.datetime.now()) + '] = Done.', file=sys.stdout)


if __name__ == '__main__':
    main()
