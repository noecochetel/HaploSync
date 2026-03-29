#!/usr/bin/env python3
"""
hapfill_coverage.py — Per-chromosome coverage extraction for HaploFill.

Processes a single chromosome from a sorted BAM file and writes:
  {out}.cov.txt.gz    — per-base coverage signal (tab-separated integers, gzipped)
  {out}.cov.bed.gz    — run-length encoded coverage ranges BED (gzipped)

These files are consumed by subsequent HaploFill steps (ploidy classification,
gap region preparation). Designed to be scattered one job per chromosome in
the HF_COVERAGE Nextflow module, enabling full parallelism across all
chromosomes of both haplotypes simultaneously.

Coverage tool options:
  bedtools (default) — bedtools genomecov -d, widely available, validated
  mosdepth           — 20-50x faster, lower RAM, opt-in for testing

Usage:
    hapfill_coverage.py -b reads.on.hap1.sorted.bam
                        -c chr01
                        -l 100000000
                        -f chr01.fasta
                        -o chr01
                        [-t bedtools|mosdepth]
"""

import argparse
import gzip
import os
import subprocess
import sys

import numpy as np


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        description='Per-chromosome coverage extraction for HaploFill',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    req = parser.add_argument_group('Required')
    req.add_argument('-b', '--bam',         required=True,
                     help='Sorted BAM file (full genome, indexed)',
                     metavar='reads.on.hap.sorted.bam')
    req.add_argument('-c', '--chr',         required=True,
                     help='Chromosome/sequence name to process',
                     metavar='chr01')
    req.add_argument('-l', '--chr_length',  required=True, type=int,
                     help='Length of the chromosome in bp',
                     metavar='N')
    req.add_argument('-f', '--chr_fasta',   required=True,
                     help='FASTA file for this chromosome (used by samtools view)',
                     metavar='chr01.fasta')
    req.add_argument('-o', '--out',         required=True,
                     help='Output prefix (outputs: {out}.cov.txt.gz, {out}.cov.bed.gz)',
                     metavar='PREFIX')

    opt = parser.add_argument_group('Optional')
    opt.add_argument('-t', '--coverage_tool', default='bedtools',
                     choices=['bedtools', 'mosdepth'],
                     help='Coverage extraction tool [default: bedtools]')

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Signal helpers (same format as HaploFunct.py)
# ---------------------------------------------------------------------------

def write_signal_file(signal, signal_file):
    """Write numpy coverage array as tab-separated gzipped text."""
    with gzip.open(signal_file, 'wt') as f:
        print('\t'.join(str(x) for x in signal), file=f)
    return signal_file


def signal2range_bed(signal, chr_name, bed_file):
    """Write run-length encoded coverage BED from numpy array."""
    with gzip.open(bed_file, 'wt') as out:
        value = None
        start = 0
        stop  = 1
        for pos in range(len(signal)):
            call = signal[pos]
            if value is None:
                value = call
                start = 0
                stop  = 1
            elif value != call:
                print('\t'.join(str(x) for x in [chr_name, start, stop, value]),
                      file=out)
                value = call
                start = pos
                stop  = pos + 1
            else:
                stop = pos + 1
        if value is not None:
            print('\t'.join(str(x) for x in [chr_name, start, stop, value]),
                  file=out)
    return bed_file


# ---------------------------------------------------------------------------
# Bedtools path
# ---------------------------------------------------------------------------

def run_bedtools_coverage(bam, chr_name, chr_length, chr_fasta, out_prefix):
    """
    Extract per-base coverage using bedtools genomecov -d.

    Steps:
      1. samtools view — subset BAM to this chromosome
      2. bedtools genomecov -d | gzip -1 — per-base coverage, fast compression
      3. Parse gzip BED into numpy int32 array (~400 MB per 100 Mb chr)
      4. Write signal file + range BED
    """
    bedtools = 'bedtools'
    samtools = 'samtools'

    # Ensure BAI index exists — create it if absent
    bai = bam + '.bai'
    if not os.path.exists(bai) and not os.path.exists(bam.replace('.bam', '.bai')):
        print('[hapfill_coverage] BAI index not found, indexing: samtools index {}'.format(bam),
              file=sys.stderr)
        subprocess.run([samtools, 'index', bam], check=True)

    chr_bam      = out_prefix + '.cov.bam'
    len_file     = out_prefix + '.len'
    cov_bed_gz   = out_prefix + '.cov.single_base.bed.gz'
    cov_err      = out_prefix + '.cov.err'
    signal_file  = out_prefix + '.cov.txt.gz'
    range_bed    = out_prefix + '.cov.bed.gz'

    # Write chromosome length file
    with open(len_file, 'w') as lf:
        print(chr_name + '\t' + str(chr_length), file=lf)

    # Subset BAM for this chromosome
    cmd = (samtools + ' view -h ' + bam + ' ' + chr_name
           + ' | ' + samtools + ' view -b -o ' + chr_bam
           + ' -T ' + chr_fasta)
    print('[hapfill_coverage] Subsetting BAM: ' + cmd, file=sys.stderr)
    subprocess.run(cmd, shell=True, check=True)

    # Index per-chromosome BAM
    subprocess.run(samtools + ' index ' + chr_bam, shell=True, check=True)

    # Extract per-base coverage (gzip -1: fast compression for intermediate file)
    cmd = (bedtools + ' genomecov -d -ibam ' + chr_bam
           + ' -g ' + len_file
           + ' 2> ' + cov_err
           + ' | gzip -1 > ' + cov_bed_gz)
    print('[hapfill_coverage] Running bedtools: ' + cmd, file=sys.stderr)
    subprocess.run(cmd, shell=True, check=True)

    # Parse into numpy array (12x less RAM than Python list of strings)
    # Use bounds-safe assignment: bedtools may use BAM-header length which can
    # differ by 1 from our computed chr_length (e.g. off-by-one at chromosome end)
    arr = np.zeros(chr_length, dtype=np.int32)
    with gzip.open(cov_bed_gz, 'rt') as fh:
        for line in fh:
            parts = line.rstrip().split('\t')
            idx = int(parts[1]) - 1
            if 0 <= idx < chr_length:
                arr[idx] = int(parts[2])

    write_signal_file(arr, signal_file)
    signal2range_bed(arr, chr_name, range_bed)

    print('[hapfill_coverage] Done (bedtools): ' + signal_file, file=sys.stderr)
    return signal_file, range_bed


# ---------------------------------------------------------------------------
# Mosdepth path
# ---------------------------------------------------------------------------

def run_mosdepth_coverage(bam, chr_name, chr_length, out_prefix):
    """
    Extract per-base coverage using mosdepth.

    Advantages over bedtools:
      - No BAM subsetting needed (-r region flag handles it)
      - Run-length encoded output: far fewer lines to parse
      - Typically 20-50x faster, lower peak RAM
      - Same output format as bedtools path — drop-in replacement

    mosdepth outputs intervals of equal depth (run-length encoded), so
    numpy slice assignment is used to fill each interval in one operation.
    """
    mosdepth = 'mosdepth'

    mosdepth_prefix = out_prefix + '.mosdepth'
    mosdepth_bed    = mosdepth_prefix + '.per-base.bed.gz'
    mosdepth_err    = out_prefix + '.mosdepth.err'
    signal_file     = out_prefix + '.cov.txt.gz'
    range_bed       = out_prefix + '.cov.bed.gz'

    # Run mosdepth — -r restricts to this chromosome, no BAM subsetting needed
    cmd = (mosdepth + ' --no-abbrev -r ' + chr_name
           + ' ' + mosdepth_prefix
           + ' ' + bam
           + ' 2> ' + mosdepth_err)
    print('[hapfill_coverage] Running mosdepth: ' + cmd, file=sys.stderr)
    subprocess.run(cmd, shell=True, check=True)

    # Parse run-length encoded BED into numpy array
    arr = np.zeros(chr_length, dtype=np.int32)
    with gzip.open(mosdepth_bed, 'rt') as fh:
        for line in fh:
            parts = line.rstrip().split('\t')
            if parts[0] != chr_name:
                continue
            start = int(parts[1])
            end   = int(parts[2])
            depth = int(parts[3])
            arr[start:end] = depth  # numpy slice assignment: one op per interval

    write_signal_file(arr, signal_file)
    signal2range_bed(arr, chr_name, range_bed)

    print('[hapfill_coverage] Done (mosdepth): ' + signal_file, file=sys.stderr)
    return signal_file, range_bed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = get_args()

    print('[hapfill_coverage] Processing chromosome: ' + args.chr, file=sys.stderr)
    print('[hapfill_coverage] Coverage tool: ' + args.coverage_tool, file=sys.stderr)
    print('[hapfill_coverage] BAM: ' + args.bam, file=sys.stderr)
    print('[hapfill_coverage] Length: ' + str(args.chr_length), file=sys.stderr)

    if args.coverage_tool == 'bedtools':
        run_bedtools_coverage(
            bam        = args.bam,
            chr_name   = args.chr,
            chr_length = args.chr_length,
            chr_fasta  = args.chr_fasta,
            out_prefix = args.out
        )
    else:
        run_mosdepth_coverage(
            bam        = args.bam,
            chr_name   = args.chr,
            chr_length = args.chr_length,
            out_prefix = args.out
        )


if __name__ == '__main__':
    main()
