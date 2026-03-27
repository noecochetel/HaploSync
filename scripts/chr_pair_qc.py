#!/usr/bin/env python3
"""
chr_pair_qc.py — Generate per-chromosome Hap1 vs Hap2 overview QC reports.

Reads HaploSplit outputs (correspondence, AGP, optional markers/groups) and
generates one HTML report per chromosome pair via chr_pair_report().

Outputs:
  {out}.chr_pair_reports/                            — per-chromosome reports
  {out}.chr_pair_reports/index.chr_pair_reports.html — index
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

    parser.add_argument('-c', '--correspondence', dest='correspondence', required=True,
        metavar='correspondence.tsv',
        help='Correspondence file (chr_id TAB hap1_id TAB hap2_id)')
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

    print('[' + str(datetime.datetime.now()) + '] = Running chr_pair_qc.py', file=sys.stdout)
    print('[' + str(datetime.datetime.now()) + '] = Command: ' + ' '.join(sys.argv), file=sys.stdout)

    # -------------------------------------------------------------------------
    # Pre-flight: verify Rscript and rmarkdown are available
    # -------------------------------------------------------------------------
    import shutil
    import subprocess as _sp
    _rscript_bin = shutil.which('Rscript') or 'Rscript'
    print('[chr_pair_qc] Rscript resolved to: ' + _rscript_bin, file=sys.stderr)
    try:
        _check = _sp.run([_rscript_bin, '-e', 'library(rmarkdown)'],
                         capture_output=True, text=True)
        if _check.returncode != 0:
            print('[chr_pair_qc] ERROR: rmarkdown check failed (exit ' + str(_check.returncode) + ')', file=sys.stderr)
            print('[chr_pair_qc]   stdout: ' + _check.stdout.strip(), file=sys.stderr)
            print('[chr_pair_qc]   stderr: ' + _check.stderr.strip(), file=sys.stderr)
            sys.exit(1)
        print('[chr_pair_qc] rmarkdown check passed.', file=sys.stderr)
    except FileNotFoundError:
        print('[chr_pair_qc] ERROR: Rscript not found at: ' + _rscript_bin, file=sys.stderr)
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Load correspondence
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] == Loading correspondence', file=sys.stdout)
    chr_ids = []
    chr_to_hap1 = {}
    chr_to_hap2 = {}
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
    else:
        marker_map_by_seq = {}

    # -------------------------------------------------------------------------
    # Load groups
    # -------------------------------------------------------------------------
    group_file = options.legacy_groups or options.input_groups
    if group_file:
        print('[' + str(datetime.datetime.now()) + '] == Loading groups', file=sys.stdout)
        groups_by_sequence = {}
        for seq_id, group_id in read_table(group_file):
            if seq_id not in groups_by_sequence:
                groups_by_sequence[seq_id] = []
            groups_by_sequence[seq_id].append(group_id)
    else:
        groups_by_sequence = {}

    # -------------------------------------------------------------------------
    # Generate chromosome pair reports
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Generating chromosome pair overview reports', file=sys.stdout)
    chr_pair_reports_dir = options.out + '.chr_pair_reports'
    mkdir(chr_pair_reports_dir)

    chr_pair_plot_db = {'Chr_Pair_Reports': {'Reports': {}}}

    def _run_chr_pair(chr_id):
        chr_pair_out_dir = chr_pair_reports_dir + '/' + chr_id
        mkdir(chr_pair_out_dir)
        result = chr_pair_report(
            chr_pair_out_dir, chr_id, {}, {}, chr_to_hap1, chr_to_hap2,
            '', agp_db, legacy_agp, markers_db, '', marker_map_by_seq,
            groups_by_sequence, options.cores, paths)
        result['html'] = chr_id + '/' + result['html']
        return chr_id, result

    with concurrent.futures.ThreadPoolExecutor(max_workers=options.cores) as executor:
        futures = {
            executor.submit(_run_chr_pair, chr_id): chr_id
            for chr_id in sorted(chr_ids)
        }
        for future in concurrent.futures.as_completed(futures):
            chr_id, result = future.result()
            print('[' + str(datetime.datetime.now()) + '] == Chr: ' + chr_id, file=sys.stdout)
            chr_pair_plot_db['Chr_Pair_Reports']['Reports'][chr_id] = result
            # Forward Rscript errors to stderr so they appear in Nextflow logs
            err_file = chr_pair_reports_dir + '/' + chr_id + '/.' + chr_id + '_chr_pair.err'
            if os.path.exists(err_file):
                err_content = open(err_file).read().strip()
                if err_content:
                    print('[chr_pair_qc] Rscript stderr for ' + chr_id + ':\n' + err_content,
                          file=sys.stderr)

    make_index_from_report_db('index.chr_pair_reports.html', '.', chr_pair_reports_dir, chr_pair_plot_db)
    print('[' + str(datetime.datetime.now()) + '] = Done.', file=sys.stdout)


if __name__ == '__main__':
    main()
