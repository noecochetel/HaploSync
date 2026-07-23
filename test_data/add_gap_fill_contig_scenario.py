#!/usr/bin/env python3
"""Add a real unplaced-contig gap-filling scenario to the gap_fill fixtures.

test_data/gap_fill_fixtures/pm01.1.fasta / pm01.Un.fasta already cover the
homozygous-fallback fill mechanism (TEST_Hap1_chr1's existing gap, resolved
via the mate haplotype since pm01.Un.fasta was empty). This script adds a
third, independent chromosome (TEST_Hap1_chr3 / TEST_Hap2_chr3) with a
deliberate gap in BOTH haplotypes at the same locus, plus a genuine unplaced
contig engineered to fill it via HaploFill's main contig-search mechanism
(STEP 6.2/6.3) - a dedicated new chromosome rather than reusing chr1/chr2, so
the existing reconstruct_pm test data (chr1/chr2, markers, genetic map) and
the chr1 homozygous-fallback test are both left untouched.

Why a dedicated, larger chromosome (20000bp) instead of just adding a gap to
the existing tiny chr1/chr2 (4000-6000bp): tile_reads() simulates reads with
a sliding window (window=800bp) that cannot start before position 0 or
extend past the chromosome end, so real simulated depth genuinely drops for
roughly the first/last ~500bp of any chromosome (confirmed empirically: on
the existing 4000bp chr2, both gap-flanking regions land at ~74.2-74.7%
diploid - just under HaploFill's 75% "reliable region" threshold, purely
because that fixed ~500bp edge effect is ~23% of a 4000bp chromosome). On a
20000bp chromosome with the gap placed centrally, the same fixed-size edge
effect becomes a small fraction of each ~9850bp flank, comfortably clearing
the 75% threshold.

Why the gap is in BOTH haplotypes, not just Hap1: with only Hap1 gapped,
HaploFill classifies the gap's own strategy as "hybrid" (Hap2's homologous
region is reliable, so it gets spliced into the search target alongside
Hap1's own flanks - see status_to_strategy()/make_sequences_and_signals() in
HaploFunct.py). That splice pads the search target with filler bases whose
own category signal doesn't cleanly classify as reliable, and classify_hit()
then rejects alignment hits that touch it - even a perfect-identity contig
alignment gets discarded before it can ever be scored. Gapping both
haplotypes at the same locus removes Hap2 as a homozygous-fallback source
*and* as a splice source, forcing the simpler, non-hybrid "gap-only" search
strategy (built purely from Hap1's own flanking regions, no splicing) - the
code path this fixture is actually meant to exercise.

Why FLANK=700 (not something smaller): analize_unplaced_hits() in
HaploFunct.py hardcodes all_threshold = 1000 - a candidate filler is only
ever recorded if its summed matched length (both flanks combined, since the
true gap content itself can't match anything under the "gap-only" strategy)
exceeds 1000bp, regardless of how high its identity/coverage percentage is.
FLANK=700 gives 2*700=1400bp of matched length, safely clearing that floor
(confirmed empirically - FLANK=300, i.e. 600bp matched, silently failed this
exact check despite a clean, coverage-threshold-passing 66.7% match).

This script only generates the FASTA/read-pool content (deterministic,
seeded, no external tools needed). Aligning the reads into the existing
BAMs requires minimap2 + samtools, which are only available where the
pipeline's own conda env is (see align_gap_fill_contig_scenario.sh, run on
that host after this script).

Run once from the test_data/ directory:
    python3 add_gap_fill_contig_scenario.py
"""
import random

random.seed(43)  # distinct from make_test_data.py's seed(42), independent sequence
BASES = "ACGT"

CHR3_LEN = 20000
GAP_START, GAP_END = 9850, 10150  # centered, far from both edges
# HaploFunct.py's analize_unplaced_hits() hardcodes all_threshold = 1000: a
# candidate filler is only ever recorded if its summed matched length (both
# flanks combined, since the true gap content itself can't match anything)
# exceeds 1000bp. FLANK=700 gives 2*700=1400bp of matched length, safely
# clearing that floor with margin (confirmed empirically - FLANK=300, i.e.
# 600bp matched, silently failed this exact check despite a clean 66.7%
# coverage-threshold-passing alignment).
FLANK = 700
UN_CONTIG_ID = "TEST_Un_chr3_filler"
READ_WINDOW = 800
READ_STEP = 150
ERROR_RATE = 0.005

FIXTURES = "gap_fill_fixtures"
HAP1_FASTA = f"{FIXTURES}/pm01.1.fasta"
HAP2_FASTA = f"{FIXTURES}/pm01.2.fasta"
UN_FASTA = f"{FIXTURES}/pm01.Un.fasta"
CORRESPONDENCE = f"{FIXTURES}/pm01.correspondence.tsv"
TRUE_HAP1_FASTA = "true_hap1.fasta"
TRUE_HAP2_FASTA = "true_hap2.fasta"
CHR3_READS_FASTA = f"{FIXTURES}/chr3_reads.pool.fasta"


def random_seq(length):
    return "".join(random.choice(BASES) for _ in range(length))


def mutate(seq, period=300, offset=150):
    """One fixed-position substitution per `period` bp -> hap2 from hap1 (mirrors make_test_data.py)."""
    seq = list(seq)
    for pos in range(offset, len(seq), period):
        alts = [b for b in BASES if b != seq[pos]]
        seq[pos] = random.choice(alts)
    return "".join(seq)


def tile_reads(seq, chrom, hap, window=READ_WINDOW, step=READ_STEP, error_rate=ERROR_RATE):
    reads = []
    last_start = max(len(seq) - window, 0)
    for start in range(0, last_start + 1, step):
        frag = list(seq[start:start + window])
        for i in range(len(frag)):
            if random.random() < error_rate:
                alts = [b for b in BASES if b != frag[i]]
                frag[i] = random.choice(alts)
        name = f"{hap}_{chrom}_{start}_{start + len(frag)}"
        reads.append((name, "".join(frag)))
    return reads


def append_fasta(path, records):
    with open(path, "a") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def main():
    hap1_chr3 = random_seq(CHR3_LEN)
    hap2_chr3 = mutate(hap1_chr3)

    true_gap_content = hap1_chr3[GAP_START:GAP_END]
    gap_mask = "N" * (GAP_END - GAP_START)
    gapped_hap1_chr3 = hap1_chr3[:GAP_START] + gap_mask + hap1_chr3[GAP_END:]
    gapped_hap2_chr3 = hap2_chr3[:GAP_START] + gap_mask + hap2_chr3[GAP_END:]
    un_contig_seq = hap1_chr3[GAP_START - FLANK:GAP_END + FLANK]

    append_fasta(TRUE_HAP1_FASTA, [("chr3", hap1_chr3)])
    append_fasta(TRUE_HAP2_FASTA, [("chr3", hap2_chr3)])
    append_fasta(HAP1_FASTA, [("TEST_Hap1_chr3", gapped_hap1_chr3)])
    append_fasta(HAP2_FASTA, [("TEST_Hap2_chr3", gapped_hap2_chr3)])
    append_fasta(UN_FASTA, [(UN_CONTIG_ID, un_contig_seq)])

    with open(CORRESPONDENCE, "a") as fh:
        fh.write("chr3\tTEST_Hap1_chr3\tTEST_Hap2_chr3\n")

    reads = tile_reads(hap1_chr3, "chr3", "hap1") + tile_reads(hap2_chr3, "chr3", "hap2")
    with open(CHR3_READS_FASTA, "w") as fh:
        for name, seq in reads:
            fh.write(f">{name}\n{seq}\n")

    print(f"TEST_Hap1_chr3: {CHR3_LEN}bp, gap [{GAP_START}:{GAP_END}) masked to N")
    print(f"{UN_CONTIG_ID}: {len(un_contig_seq)}bp ([{GAP_START - FLANK}:{GAP_END + FLANK}) of the true sequence)")
    print(f"True gap content ({len(true_gap_content)}bp): {true_gap_content}")
    print(f"Wrote {len(reads)} reads to {CHR3_READS_FASTA}")
    print("Next: run align_gap_fill_contig_scenario.sh where minimap2/samtools are available.")


if __name__ == "__main__":
    main()
