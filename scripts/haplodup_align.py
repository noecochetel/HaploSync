#!/usr/bin/env python3
"""
haplodup_align.py — Run pairwise nucmer alignments for HaploDup.

Reads FASTA and correspondence, runs all pairwise nucmer comparisons
(Hap1×Hap1, Hap2×Hap2, Hap1×Hap2, Hap2×Hap1, and optionally vs Reference),
and writes .delta files to {out}.HaploDup_dir/.

Outputs:
  {out}.HaploDup_dir/*.delta   — nucmer alignment files (consumed by HAPLODUP_REPORT)
"""

import argparse
import datetime
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

from lib_files.HaploFunct import *
from lib_files.FASTA_lib import *
from lib_files.map_lib import *


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__)

    parser.add_argument('-f', '--fasta', dest='fasta', required=True,
        metavar='hap1.fasta,hap2.fasta,...',
        help='Comma-separated FASTA file(s) with genomic sequences')
    parser.add_argument('-c', '--correspondence', dest='corr', required=True,
        metavar='correspondence.tsv',
        help='Correspondence file (chr TAB hap1 TAB hap2 [TAB ref])')
    parser.add_argument('-r', '--reference', dest='reference', default=None,
        metavar='reference.fasta',
        help='Reference genome FASTA (haploid, optional)')
    parser.add_argument('-o', '--out', dest='out', default='out',
        metavar='PREFIX',
        help='Output prefix [default: out]')
    parser.add_argument('-t', '--threads', dest='cores', default=4, type=int,
        metavar='N',
        help='CPU cores for nucmer [default: 4]')
    parser.add_argument('--conf', dest='conf', default=None,
        metavar='HaploSync.conf.toml',
        help='Path to HaploSync conf.toml [default: auto-detect]')

    options = parser.parse_args()

    conf_path = options.conf or os.path.join(
        os.path.dirname(os.path.realpath(__file__)), '..', 'HaploSync.conf.toml')
    paths = set_paths(conf_path)

    print('[' + str(datetime.datetime.now()) + '] = Running haplodup_align.py', file=sys.stdout)
    print('[' + str(datetime.datetime.now()) + '] = Command: ' + ' '.join(sys.argv), file=sys.stdout)

    haplodup_dir = options.out + '.HaploDup_dir'
    mkdir(haplodup_dir)

    # -------------------------------------------------------------------------
    # Load FASTA
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] == Loading FASTA sequences', file=sys.stdout)
    fasta_dict = {}
    for fasta_file in options.fasta.split(','):
        fasta_file = fasta_file.strip()
        if fasta_file:
            fasta_dict.update(read_fasta(fasta_file))

    if options.reference:
        reference = read_fasta(options.reference)
        reference_file = haplodup_dir + '/' + options.out + '.ref.fasta'
        write_fasta_from_db(reference, reference_file, False)

    # -------------------------------------------------------------------------
    # Load correspondence and split into hap1/hap2 FASTA dicts
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] == Loading correspondence', file=sys.stdout)
    hap1_fasta = {}
    hap2_fasta = {}

    for line in open(options.corr):
        line = line.rstrip()
        if not line:
            continue
        parts = line.split('\t')
        if options.reference:
            if len(parts) < 4:
                continue
            chr_id, seq1, seq2, ref = parts[:4]
        else:
            if len(parts) < 3:
                continue
            chr_id, seq1, seq2 = parts[:3]
        if seq1 in fasta_dict:
            hap1_fasta[seq1] = fasta_dict[seq1]
        if seq2 in fasta_dict:
            hap2_fasta[seq2] = fasta_dict[seq2]

    # Write split FASTA files (needed by HaploDup.py with --reuse_mappings)
    query_1_file = haplodup_dir + '/' + options.out + '.1.fasta'
    query_2_file = haplodup_dir + '/' + options.out + '.2.fasta'
    write_fasta_from_db(hap1_fasta, query_1_file)
    write_fasta_from_db(hap2_fasta, query_2_file)

    # -------------------------------------------------------------------------
    # Run nucmer alignments
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Running nucmer alignments', file=sys.stdout)

    print('[' + str(datetime.datetime.now()) + '] == Hap1 vs Hap1', file=sys.stdout)
    map_nucmer_dotplot('Hap1', query_1_file, 'Hap1', query_1_file, haplodup_dir, options.cores, paths, False)

    print('[' + str(datetime.datetime.now()) + '] == Hap2 vs Hap2', file=sys.stdout)
    map_nucmer_dotplot('Hap2', query_2_file, 'Hap2', query_2_file, haplodup_dir, options.cores, paths, False)

    print('[' + str(datetime.datetime.now()) + '] == Hap2 vs Hap1', file=sys.stdout)
    map_nucmer_dotplot('Hap1', query_1_file, 'Hap2', query_2_file, haplodup_dir, options.cores, paths, False)

    print('[' + str(datetime.datetime.now()) + '] == Hap1 vs Hap2', file=sys.stdout)
    map_nucmer_dotplot('Hap2', query_2_file, 'Hap1', query_1_file, haplodup_dir, options.cores, paths, False)

    if options.reference:
        ref_file = haplodup_dir + '/' + options.out + '.ref.fasta'
        print('[' + str(datetime.datetime.now()) + '] == Hap1 vs Reference', file=sys.stdout)
        map_nucmer_dotplot('Ref', ref_file, 'Hap1', query_1_file, haplodup_dir, options.cores, paths, False)

        print('[' + str(datetime.datetime.now()) + '] == Reference vs Hap1', file=sys.stdout)
        map_nucmer_dotplot('Hap1', query_1_file, 'Ref', ref_file, haplodup_dir, options.cores, paths, False)

        print('[' + str(datetime.datetime.now()) + '] == Hap2 vs Reference', file=sys.stdout)
        map_nucmer_dotplot('Ref', ref_file, 'Hap2', query_2_file, haplodup_dir, options.cores, paths, False)

        print('[' + str(datetime.datetime.now()) + '] == Reference vs Hap2', file=sys.stdout)
        map_nucmer_dotplot('Hap2', query_2_file, 'Ref', ref_file, haplodup_dir, options.cores, paths, False)

    print('[' + str(datetime.datetime.now()) + '] = Done.', file=sys.stdout)


if __name__ == '__main__':
    main()
