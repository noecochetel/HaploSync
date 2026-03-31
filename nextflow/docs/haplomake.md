# HaploMake

**Entry point:** `nextflow/haplomake.nf`

Also available as a post-pipeline convenience entry point:
- `nextflow/gap_fill.nf --run_haplomake` — runs automatically after HaploFill
- `nextflow/gap_fill.nf -entry HAPLOMAKE` — auto-reads structure block from `{outdir}/HaploFill/`

Constructs new pseudomolecule FASTA and AGP files from a structure description file. The structure file defines the ordered composition of each output sequence — which source sequences to include, in what orientation, and with what gap sizes between components. Three input formats are accepted: a HaploFill structure block (`BLOCK`), an AGP file (`AGP`), or a BED file (`BED`).

HaploMake is used both as the final step of the gap-filling pipeline and as a standalone sequence-editing tool for manual curation (e.g., splitting overassembled contigs before re-running HaploSplit).

---

## What it does

HaploMake reads a structure file and:

1. Extracts the specified slices from the source FASTA files
2. Concatenates components in defined order with gaps (`N` stretches) between them
3. Trims overlapping joins (unless `--skipoverlap`)
4. Produces a new FASTA with optionally renamed sequences (`--hapmake_prefix`)
5. Optionally translates coordinates: AGP, BED, and GFF3 into the new sequence space

---

## Structure file formats

### BLOCK (default)

The HaploFill structure block format, produced automatically by the gap-filling pipeline. One output sequence per `>` header, with component rows below it:

```
>NEW_Chr01
source_seq  start  end  strand
source_seq  start  end  strand
>NEW_Chr02
...
```

### AGP

Standard NCBI AGP format. The most convenient format for manual curation — take an existing AGP from HaploSplit or HaploMake, edit it to split or merge sequences, and pass it to HaploMake with `--format AGP`.

```
# object      obj_start  obj_end  part_num  comp_type  comp_name     comp_start  comp_end  orientation
NEW_Chr01     1          5000000  1         W          contig_001    1           5000000   +
NEW_Chr01     5000001    5001000  2         N          1000          scaffold    yes       na
NEW_Chr01     5001001    9500000  3         W          contig_002    1           4499000   +
NEW_Chr02     1          3000000  1         W          contig_003    1           3000000   -
```

Component types: `W` = sequence component, `N` = gap. For `W` rows: columns 6–9 are source sequence name, start, end, and orientation. For `N` rows: column 6 is gap length.

Use `--reverse` to go in the opposite direction: extract old coordinates from a new-to-old AGP (requires the new FASTA as input).

### BED

A BED file listing regions to extract. BED3 and BED6 formats are supported. Each input file produces one output sequence; the output sequence name is set by `--hapmake_prefix` + a progressive number.

```
# BED3: chrom  start(0-based)  end
contig_001    0       2500000
contig_002    500000  3000000

# BED6: chrom  start  end  name  score  strand
contig_001    0       2500000  region_A  0  +
contig_002    500000  3000000  region_B  0  -
```

---

## Entry points

### Standalone

```bash
nextflow run nextflow/haplomake.nf -profile mamba \
    --hap1_fasta hap1.fasta --hap2_fasta hap2.fasta \
    --structure_block myproject.structure.block \
    --out myproject --outdir results
```

### Within the gap-filling pipeline

Run automatically when `--run_haplomake` or `--run_haplodup` is set:

```bash
nextflow run nextflow/gap_fill.nf -profile mamba \
    --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \
    --hapfill_unplaced unplaced.fasta \
    --hapfill_correspondence correspondence.tsv \
    --hapfill_repeats repeats.bed \
    --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \
    --run_haplomake \
    --out myproject --outdir results
```

### After gap filling (convenience)

Reads the structure block automatically from `{outdir}/HaploFill/{out}.structure.block`. Use `--structure_block` to override:

```bash
nextflow run nextflow/gap_fill.nf -entry HAPLOMAKE -profile mamba \
    --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \
    --out myproject --outdir results
```

---

## Parameters

### Required (standalone)

| Parameter | Description |
|-----------|-------------|
| `--hap1_fasta` | Hap1 FASTA |
| `--hap2_fasta` | Hap2 FASTA |
| `--structure_block` | Structure block file |

### Optional inputs

| Parameter | Description |
|-----------|-------------|
| `--unplaced_fasta` | Unplaced sequences FASTA |

### HaploMake options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hapmake_prefix` | — | Sequence ID prefix for output sequences |
| `--hapmake_gap` | 1000 | Gap size in bp between components |
| `--hapmake_skipoverlap` | false | Skip overlap trimming at joins |
| `--hapmake_noagp` | false | Skip AGP output |
| `--hapmake_unplaced` | — | Override unplaced sequences FASTA |

### Coordinate translation (optional)

| Parameter | Description |
|-----------|-------------|
| `--hapmake_agp` | Input AGP to lift over into the new assembly space |
| `--hapmake_gff3` | Gene annotation GFF3 to translate |
| `--hapmake_bed` | BED file to translate |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--out` | `out` | Output files prefix |
| `--outdir` | `results` | Results directory |

---

## Output files

Written to `{outdir}/HaploMake/`:

| File | Description |
|------|-------------|
| `{out}.fasta` | New pseudomolecule FASTA |
| `{out}.structure.agp` | AGP structure of the new assembly |
| `{out}.legacy_structure.agp` | Lifted-over input AGP (if `--hapmake_agp`) |
| `{out}.loci_to_check.txt` | Regions needing manual review (if overlaps found) |

---

## Examples

```bash
# Standalone — from a HaploFill structure block
nextflow run nextflow/haplomake.nf -profile mamba \
    --hap1_fasta hap1.fasta --hap2_fasta hap2.fasta \
    --unplaced_fasta unplaced.fasta \
    --structure_block myproject.structure.block \
    --out myproject --outdir results

# Standalone — manual curation from an edited AGP (e.g., split an overassembled contig)
nextflow run nextflow/haplomake.nf -profile mamba \
    --hap1_fasta hap1.fasta --hap2_fasta hap2.fasta \
    --structure_block hap1_corrected.agp \
    --hapmake_prefix NEW \
    --out myproject_corrected --outdir results

# Standalone — with AGP and annotation translation
nextflow run nextflow/haplomake.nf -profile mamba \
    --hap1_fasta hap1.fasta --hap2_fasta hap2.fasta \
    --structure_block myproject.structure.block \
    --hapmake_agp previous.agp \
    --hapmake_gff3 annotation.gff3 \
    --out myproject --outdir results

# As part of the gap-filling pipeline
nextflow run nextflow/gap_fill.nf -profile mamba \
    --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \
    --hapfill_correspondence correspondence.tsv \
    --hapfill_repeats repeats.bed \
    --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \
    --run_haplomake \
    --hapmake_agp previous.agp \
    --hapmake_gff3 annotation.gff3 \
    --out myproject --outdir results

# HPC (SLURM)
nextflow run nextflow/haplomake.nf -profile hpc \
    -params-file nextflow/params_haplomake.yml
```
