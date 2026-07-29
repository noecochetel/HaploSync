#!/usr/bin/env python3
"""Add synthetic gene (GFF3) and marker (BED) annotation to the bundled test
fixtures, so --gff3 (reconstruct_pm) and --hapmake_gff3/--hapmake_bed
(gap_fill) can be exercised by the test profiles. Without these, HaploDup's
chromosome-pair reports render with has_genes=FALSE (and, for gap_fill,
has_markers=FALSE too), skipping the gene-count and marker tracks entirely.

Produces (all deterministic, no randomness needed - positions are picked by
hand against the actual fixture sequence lengths and N-gap locations):
  test_data/genes.gff3                    - genes for --gff3, on
                                             input_assembly.fasta contigs
                                             (reconstruct_pm)
  test_data/gap_fill_fixtures/genes.gff3   - genes for --hapmake_gff3, on
                                             pm01.{1,2}.fasta pseudomolecules
                                             (gap_fill)
  test_data/gap_fill_fixtures/markers.bed  - markers for --hapmake_bed, on
                                             the same pm01.{1,2}.fasta
                                             pseudomolecules (gap_fill)

Each gene is a single-exon gene/mRNA/CDS triplet (sufficient for GMAP-based
CDS-vs-genome mapping). Positions were chosen to avoid every N-gap already
present in the underlying FASTA (verified with a scan of the actual files,
not assumed):
  input_assembly.fasta:              contig2         N-gap [1700,2000)
  gap_fill_fixtures/pm01.1.fasta:    TEST_Hap1_chr1  N-gaps [3000,4000), [5700,6000)
                                     TEST_Hap1_chr3  N-gap  [9850,10150)
  gap_fill_fixtures/pm01.2.fasta:    TEST_Hap2_chr3  N-gap  [9850,10150)
"""
import os

GENE_LEN = 300


def read_fasta_lengths(path):
    lengths = {}
    seq_id = None
    length = 0
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if seq_id:
                lengths[seq_id] = length
            seq_id = line[1:].split()[0]
            length = 0
        else:
            length += len(line)
    if seq_id:
        lengths[seq_id] = length
    return lengths


def genes_at(seq_id, positions, prefix):
    """positions: 0-based start coordinates. Returns 1-based inclusive
    (seq_id, start, stop, gene_id) tuples, each GENE_LEN bp long."""
    genes = []
    for i, pos in enumerate(positions, 1):
        genes.append((seq_id, pos + 1, pos + GENE_LEN, f"{prefix}_g{i:02d}"))
    return genes


def write_gff3(path, genes):
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for seq_id, start, stop, gene_id in genes:
            mrna_id = gene_id + ".1"
            cds_id = gene_id + ".1.CDS1"
            fh.write(f"{seq_id}\tsynthetic\tgene\t{start}\t{stop}\t.\t+\t.\tID={gene_id}\n")
            fh.write(f"{seq_id}\tsynthetic\tmRNA\t{start}\t{stop}\t.\t+\t.\tID={mrna_id};Parent={gene_id}\n")
            fh.write(f"{seq_id}\tsynthetic\tCDS\t{start}\t{stop}\t.\t+\t0\tID={cds_id};Parent={mrna_id}\n")


def write_bed(path, markers):
    with open(path, "w") as fh:
        for seq_id, start, stop, marker_id in markers:
            fh.write(f"{seq_id}\t{start}\t{stop}\t{marker_id}\n")


def main():
    script_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(script_dir)

    # --- reconstruct_pm: genes.gff3 on input_assembly.fasta contigs ---
    asm_lengths = read_fasta_lengths("input_assembly.fasta")
    recon_genes = []
    recon_genes += genes_at("contig1", [200, 1400, 2400], "c1")
    recon_genes += genes_at("contig2", [200, 2200, 2600], "c2")  # avoids N-gap [1700,2000)
    recon_genes += genes_at("contig3", [200, 1400, 2400, 3200, 5200, 5600], "c3")
    recon_genes += genes_at("contig4", [200, 2200], "c4")
    recon_genes += genes_at("contig5", [200, 2200], "c5")
    for seq_id, _, stop, _ in recon_genes:
        assert stop <= asm_lengths[seq_id], f"{seq_id} gene runs past sequence end"
    write_gff3("genes.gff3", recon_genes)

    # --- gap_fill: genes.gff3 + markers.bed on pm01.{1,2}.fasta pseudomolecules ---
    hap1_lengths = read_fasta_lengths("gap_fill_fixtures/pm01.1.fasta")
    hap2_lengths = read_fasta_lengths("gap_fill_fixtures/pm01.2.fasta")

    gf_genes = []
    # chr1: TEST_Hap1_chr1 has N-gaps at [3000,4000) and [5700,6000) - avoid
    gf_genes += genes_at("TEST_Hap1_chr1", [500, 4500, 6300], "GFH1c1")
    gf_genes += genes_at("TEST_Hap2_chr1", [500, 2500, 4500], "GFH2c1")
    # chr2: no gaps in either haplotype
    gf_genes += genes_at("TEST_Hap1_chr2", [500, 2500], "GFH1c2")
    gf_genes += genes_at("TEST_Hap2_chr2", [500, 2500], "GFH2c2")
    # chr3: both haplotypes gapped at the same [9850,10150) window - avoid
    gf_genes += genes_at("TEST_Hap1_chr3", [500, 5000, 15000, 19000], "GFH1c3")
    gf_genes += genes_at("TEST_Hap2_chr3", [500, 5000, 15000, 19000], "GFH2c3")
    for seq_id, _, stop, _ in gf_genes:
        lengths = hap1_lengths if seq_id.startswith("TEST_Hap1") else hap2_lengths
        assert stop <= lengths[seq_id], f"{seq_id} gene runs past sequence end"
    write_gff3("gap_fill_fixtures/genes.gff3", gf_genes)

    gf_markers = [
        ("TEST_Hap1_chr1", 1000, 1001, "gfm01"),
        ("TEST_Hap1_chr1", 6500, 6501, "gfm02"),
        ("TEST_Hap2_chr1", 1000, 1001, "gfm01"),
        ("TEST_Hap2_chr1", 5000, 5001, "gfm03"),
        ("TEST_Hap1_chr2", 1000, 1001, "gfm04"),
        ("TEST_Hap2_chr2", 1000, 1001, "gfm04"),
        ("TEST_Hap1_chr3", 1000, 1001, "gfm05"),
        ("TEST_Hap1_chr3", 18000, 18001, "gfm06"),
        ("TEST_Hap2_chr3", 1000, 1001, "gfm05"),
        ("TEST_Hap2_chr3", 18000, 18001, "gfm06"),
    ]
    write_bed("gap_fill_fixtures/markers.bed", gf_markers)

    print(f"Wrote {len(recon_genes)} reconstruct_pm genes to genes.gff3")
    print(f"Wrote {len(gf_genes)} gap_fill genes to gap_fill_fixtures/genes.gff3")
    print(f"Wrote {len(gf_markers)} gap_fill markers to gap_fill_fixtures/markers.bed")


if __name__ == "__main__":
    main()
