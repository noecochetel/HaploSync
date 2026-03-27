#!/usr/bin/env python3
"""
haplodup_gmap.py — Run GMAP gene mapping for HaploDup.

Reads FASTA, GFF3 annotation, and correspondence, extracts CDS sequences,
builds GMAP indices for Hap1 and Hap2, maps CDS to both haplotypes, and
concatenates results.

Outputs:
  {out}.HaploDup_dir/CDS.on.genome.gmap.gff3   — GMAP mapping results (consumed by HAPLODUP_REPORT)
"""

import argparse
import datetime
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

from lib_files.HaploFunct import *
from lib_files.GFF_lib import *
from lib_files.FASTA_lib import *


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__)

    parser.add_argument('-f', '--fasta', dest='fasta', required=True,
        metavar='hap1.fasta,hap2.fasta,...',
        help='Comma-separated FASTA file(s) with genomic sequences')
    parser.add_argument('-c', '--correspondence', dest='corr', required=True,
        metavar='correspondence.tsv',
        help='Correspondence file (chr TAB hap1 TAB hap2)')
    parser.add_argument('-g', '--gff', dest='gff', required=True,
        metavar='annotation.gff3',
        help='Gene annotation file(s) in GFF3 format (comma-separated)')
    parser.add_argument('-o', '--out', dest='out', default='out',
        metavar='PREFIX',
        help='Output prefix [default: out]')
    parser.add_argument('-t', '--threads', dest='cores', default=4, type=int,
        metavar='N',
        help='CPU cores for GMAP [default: 4]')
    parser.add_argument('--feature', dest='feature', default='CDS',
        metavar='[CDS|mRNA]',
        help='Feature type to use for mapping [default: CDS]')
    parser.add_argument('--conf', dest='conf', default=None,
        metavar='HaploSync.conf.toml',
        help='Path to HaploSync conf.toml [default: auto-detect]')

    options = parser.parse_args()

    conf_path = options.conf or os.path.join(
        os.path.dirname(os.path.realpath(__file__)), '..', 'HaploSync.conf.toml')
    paths = set_paths(conf_path)

    print('[' + str(datetime.datetime.now()) + '] = Running haplodup_gmap.py', file=sys.stdout)
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
        if len(parts) < 3:
            continue
        chr_id, seq1, seq2 = parts[:3]
        if seq1 in fasta_dict:
            hap1_fasta[seq1] = fasta_dict[seq1]
        if seq2 in fasta_dict:
            hap2_fasta[seq2] = fasta_dict[seq2]

    # -------------------------------------------------------------------------
    # Load GFF3 annotation
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] == Loading GFF3 annotation', file=sys.stdout)
    gff_db = {}
    mRNA_db = {}
    for gff_file in options.gff.split(','):
        gff_file = gff_file.strip()
        if not gff_file:
            continue
        gff_tmp, mrna_tmp = read_gff3(gff_file)
        gff_db.update(gff_tmp)
        mRNA_db.update(mrna_tmp)

    # -------------------------------------------------------------------------
    # Extract CDS and run GMAP
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Extracting CDS sequences', file=sys.stdout)
    CDS_file = get_sequence(gff_db, fasta_dict, haplodup_dir + '/new', 'CDS')

    index_dir = haplodup_dir + '/gmap_index'
    mkdir(index_dir)

    # Hap1 — index
    print('[' + str(datetime.datetime.now()) + '] == Hap1: indexing', file=sys.stdout)
    hap1_name = haplodup_dir + '/new.hap1.fasta'
    write_fasta_from_db(hap1_fasta, hap1_name)
    indexing_out = open(haplodup_dir + '/gmap_index.1.log', 'w')
    indexing_err = open(haplodup_dir + '/gmap_index.1.err', 'w')
    indexProcess = subprocess.Popen(
        'gmap_build -D ' + index_dir + ' -d hap1.fasta ' + hap1_name,
        shell=True, stdout=indexing_out, stderr=indexing_err)
    indexProcess.communicate()
    indexing_out.close()
    indexing_err.close()

    # Hap1 — GMAP
    print('[' + str(datetime.datetime.now()) + '] == Hap1: mapping with GMAP', file=sys.stdout)
    gmap_results_1 = haplodup_dir + '/CDS.on.hap1.gmap.gff3'
    gmap_1_gff3 = open(gmap_results_1, 'w')
    gmap_err = open(gmap_results_1 + '.err', 'w')
    gmapProcess = subprocess.Popen(
        'gmap -D ' + index_dir + ' -d hap1.fasta -f 2 -n 500 -t ' + str(options.cores) + ' ' + CDS_file,
        shell=True, stdout=gmap_1_gff3, stderr=gmap_err)
    gmapProcess.communicate()
    gmap_1_gff3.close()
    gmap_err.close()

    # Hap2 — index
    print('[' + str(datetime.datetime.now()) + '] == Hap2: indexing', file=sys.stdout)
    hap2_name = haplodup_dir + '/new.hap2.fasta'
    write_fasta_from_db(hap2_fasta, hap2_name)
    indexing_out = open(haplodup_dir + '/gmap_index.2.log', 'w')
    indexing_err = open(haplodup_dir + '/gmap_index.2.err', 'w')
    indexProcess = subprocess.Popen(
        'gmap_build -D ' + index_dir + ' -d hap2.fasta ' + hap2_name,
        shell=True, stdout=indexing_out, stderr=indexing_err)
    indexProcess.communicate()
    indexing_out.close()
    indexing_err.close()

    # Hap2 — GMAP
    print('[' + str(datetime.datetime.now()) + '] == Hap2: mapping with GMAP', file=sys.stdout)
    gmap_results_2 = haplodup_dir + '/CDS.on.hap2.gmap.gff3'
    gmap_2_gff3 = open(gmap_results_2, 'w')
    gmap_err = open(gmap_results_2 + '.err', 'w')
    gmapProcess = subprocess.Popen(
        'gmap -D ' + index_dir + ' -d hap2.fasta -f 2 -n 500 -t ' + str(options.cores) + ' ' + CDS_file,
        shell=True, stdout=gmap_2_gff3, stderr=gmap_err)
    gmapProcess.communicate()
    gmap_2_gff3.close()
    gmap_err.close()

    # Concatenate results
    print('[' + str(datetime.datetime.now()) + '] == Concatenating GMAP results', file=sys.stdout)
    gmap_results = haplodup_dir + '/CDS.on.genome.gmap.gff3'
    with open(gmap_results, 'w') as gmap_gff3:
        for line in open(gmap_results_1):
            print(line.rstrip(), file=gmap_gff3)
        for line in open(gmap_results_2):
            print(line.rstrip(), file=gmap_gff3)

    print('[' + str(datetime.datetime.now()) + '] = Done. Output: ' + gmap_results, file=sys.stdout)


if __name__ == '__main__':
    main()
