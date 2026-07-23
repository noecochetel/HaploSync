#!/usr/bin/env bash
# Realign the gap_fill_fixtures BAMs after add_gap_fill_contig_scenario.py has
# appended TEST_Hap1_chr3/TEST_Hap2_chr3 to pm01.1.fasta/pm01.2.fasta and
# written gap_fill_fixtures/chr3_reads.pool.fasta.
#
# Run where minimap2 + samtools are available (the pipeline's own conda env).
# Mirrors the original fixture generation commands (see BAM header: minimap2
# -ax map-hifi ... reads.pool.fasta), just with chr3's reads folded in, so
# the whole genome (chr1+chr2+chr3) is re-aligned in one pass rather than
# merging separately-aligned BAMs.
set -euo pipefail
cd "$(dirname "$0")/gap_fill_fixtures"

cat ../reads.pool.fasta chr3_reads.pool.fasta > combined_reads.pool.fasta

minimap2 -ax map-hifi -t 2 pm01.1.fasta combined_reads.pool.fasta | samtools sort -o hap1.on.hap1.sorted.bam -
samtools index hap1.on.hap1.sorted.bam

minimap2 -ax map-hifi -t 2 pm01.2.fasta combined_reads.pool.fasta | samtools sort -o hap2.on.hap2.sorted.bam -
samtools index hap2.on.hap2.sorted.bam

rm combined_reads.pool.fasta
echo "Realigned hap1.on.hap1.sorted.bam / hap2.on.hap2.sorted.bam (chr1+chr2+chr3)"
