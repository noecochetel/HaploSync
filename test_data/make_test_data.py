#!/usr/bin/env python3
"""Generate a tiny synthetic diploid genome for a local HaploSync pipeline test.

Produces (all deterministic, seeded):
  - input_assembly.fasta   pooled hifiasm-like contigs (HaploSplit input)
  - genetic_map.tsv        genetic map (HaploSplit --markers_map)
  - markers_hits.bed       marker hits on contigs (HaploSplit --markers)
  - true_hap1.fasta        true, gap-free hap1 chromosomes (for read simulation only)
  - true_hap2.fasta        true, gap-free hap2 chromosomes (for read simulation only)
  - reads.pool.fasta       deterministic tiled "reads" for HaploFill coverage BAMs
  - repeats.empty.bed      empty repeats file (HaploFill --hapfill_repeats)
"""
import random

random.seed(42)
BASES = "ACGT"

CHR1_LEN = 6000
CHR2_LEN = 4000
GAP_LOCAL_START, GAP_LOCAL_END = 1700, 2000  # within contig2 (true 4700-5000)
READ_WINDOW = 800
READ_STEP = 150
ERROR_RATE = 0.005


def random_seq(length):
    return "".join(random.choice(BASES) for _ in range(length))


def mutate(seq, period=300, offset=150):
    """One fixed-position substitution per `period` bp -> hap2 from hap1."""
    seq = list(seq)
    for pos in range(offset, len(seq), period):
        alts = [b for b in BASES if b != seq[pos]]
        seq[pos] = random.choice(alts)
    return "".join(seq)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


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


def main():
    # --- True, gap-free haplotype chromosomes ---
    hap1_chr1 = random_seq(CHR1_LEN)
    hap1_chr2 = random_seq(CHR2_LEN)
    hap2_chr1 = mutate(hap1_chr1)
    hap2_chr2 = mutate(hap1_chr2)

    write_fasta("true_hap1.fasta", [("chr1", hap1_chr1), ("chr2", hap1_chr2)])
    write_fasta("true_hap2.fasta", [("chr1", hap2_chr1), ("chr2", hap2_chr2)])

    # --- Hifiasm-like pooled input contigs ---
    h1_chr1_p1 = hap1_chr1[0:3000]
    p2_seq = list(hap1_chr1[3000:6000])
    for i in range(GAP_LOCAL_START, GAP_LOCAL_END):
        p2_seq[i] = "N"
    h1_chr1_p2 = "".join(p2_seq)

    h2_chr1_full = hap2_chr1
    h1_chr2_full = hap1_chr2
    h2_chr2_full = hap2_chr2

    # Neutral, arbitrary contig names - deliberately not hinting at which
    # haplotype/chromosome each is expected to land in, so test failures can't
    # be masked by a name that happens to match the expected outcome.
    write_fasta("input_assembly.fasta", [
        ("contig1", h1_chr1_p1),
        ("contig2", h1_chr1_p2),
        ("contig3", h2_chr1_full),
        ("contig4", h1_chr2_full),
        ("contig5", h2_chr2_full),
    ])

    # --- Genetic map ---
    chr1_markers = [(300 + 600 * i, f"m{i + 1:02d}") for i in range(10)]   # 300..5700
    chr2_markers = [(300 + 600 * i, f"c2m{i + 1:02d}") for i in range(7)]  # 300..3900

    with open("genetic_map.tsv", "w") as fh:
        for pos, mk in chr1_markers:
            fh.write(f"chr1\t{pos}\t{mk}\n")
        for pos, mk in chr2_markers:
            fh.write(f"chr2\t{pos}\t{mk}\n")

    # --- Marker hits per contig ---
    hits = []
    for pos, mk in chr1_markers:
        if pos < 3000:
            hits.append(("contig1", pos, pos + 1, mk))
    for pos, mk in chr1_markers:
        if pos >= 3000:
            local = pos - 3000
            hits.append(("contig2", local, local + 1, mk))
    for pos, mk in chr1_markers:
        if mk == "m05":  # deliberately omitted -> contig1+contig2 out-score contig3 on chr1
            continue
        hits.append(("contig3", pos, pos + 1, mk))
    for pos, mk in chr2_markers:
        hits.append(("contig4", pos, pos + 1, mk))
    for pos, mk in chr2_markers:
        # Omit the last marker (not an interior one like chr1's m05): both chr2
        # contigs are single/unfragmented, so an interior omission leaves contig4
        # and contig5 with the identical [300:3900] marker range, which
        # markers_to_network() (bin/lib_files/HaploFunct.py) collapses onto one
        # graph edge - silently dropping one candidate regardless of marker
        # count. A shorter contig5 range avoids that collision so contig4
        # genuinely out-scores contig5 as intended.
        if mk == "c2m07":
            continue
        hits.append(("contig5", pos, pos + 1, mk))

    with open("markers_hits.bed", "w") as fh:
        for seq_id, start, stop, mk in hits:
            fh.write(f"{seq_id}\t{start}\t{stop}\t{mk}\n")

    # --- Empty repeats file ---
    open("repeats.empty.bed", "w").close()

    # --- Simulated reads (deterministic tiling, both true haplotypes) ---
    reads = []
    reads += tile_reads(hap1_chr1, "chr1", "hap1")
    reads += tile_reads(hap1_chr2, "chr2", "hap1")
    reads += tile_reads(hap2_chr1, "chr1", "hap2")
    reads += tile_reads(hap2_chr2, "chr2", "hap2")
    write_fasta("reads.pool.fasta", reads)

    print(f"Wrote {len(hits)} marker hits, {len(reads)} reads.")


if __name__ == "__main__":
    main()
