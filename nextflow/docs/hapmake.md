# HaploMake

**Entry points:**
- `nextflow/gap_fill.nf --run_haplomake` — as part of the gap-filling pipeline
- `nextflow/gap_fill.nf -entry HAPLOMAKE` — standalone, reads structure block from `{outdir}/HaploFill/`
- `HaploMake.py` directly — for custom structure files (AGP, BED, or block)

Constructs new pseudomolecule FASTA and AGP files from a structure description file. The structure file defines the ordered composition of each output sequence — which source sequences to include, in what orientation, and with what gap sizes between components. Three input formats are accepted: a HaploFill structure block (`BLOCK`), an AGP file (`AGP`), or a BED file (`BED`).

HaploMake is used both as the final step of the gap-filling pipeline and as a standalone sequence-editing tool for manual curation.

---

## What it does

HaploMake reads a structure file and:

1. Extracts the specified slices from the source FASTA files
2. Concatenates components in defined order with gaps (`N` stretches) between them
3. Trims overlapping joins (unless `--skipoverlap`)
4. Produces a new FASTA with optionally renamed sequences (`-p`/`--prefix`)
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

A BED file listing regions to extract. BED3 and BED6 formats are supported. Each input file produces one output sequence; the output sequence name is set by `-p`/`--prefix` + a progressive number.

```
# BED3: chrom  start(0-based)  end
contig_001    0       2500000
contig_002    500000  3000000

# BED6: chrom  start  end  name  score  strand
contig_001    0       2500000  region_A  0  +
contig_002    500000  3000000  region_B  0  -
```

---

## Usage

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

### Standalone via Nextflow

Reads the structure block automatically from `{outdir}/HaploFill/{out}.structure.block`. Use `--structure_block` to provide a custom path (e.g., a manually edited AGP or block file):

```bash
# Using the structure block from a previous gap-fill run
nextflow run nextflow/gap_fill.nf -entry HAPLOMAKE -profile mamba \
    --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \
    --hapfill_unplaced unplaced.fasta \
    --out myproject --outdir results

# With a custom structure block (e.g., manually curated)
nextflow run nextflow/gap_fill.nf -entry HAPLOMAKE -profile mamba \
    --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \
    --structure_block myproject_curated.structure.block \
    --out myproject_curated --outdir results
```

### Standalone via HaploMake.py

For full control over input format (AGP, BED, or block), call `HaploMake.py` directly:

```bash
# From an AGP (typical manual curation workflow)
python3 HaploMake.py \
    -f hap1.fasta,hap2.fasta \
    -s edited_assembly.agp \
    --format AGP \
    -o myproject_corrected \
    -p NEW

# From a structure block
python3 HaploMake.py \
    -f hap1.fasta,hap2.fasta,unplaced.fasta \
    -s myproject.structure.block \
    -o myproject_corrected
```

---

## Parameters (Nextflow)

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

## HaploMake.py options

For direct use of `HaploMake.py`:

| Flag | Description |
|------|-------------|
| `-f FASTAS` | Comma-separated input FASTA files |
| `-s FILE` | Structure file (block, AGP, or BED) |
| `--format` | Format: `BLOCK` \| `AGP` \| `BED` (default: `BLOCK`) |
| `-o PREFIX` | Output prefix |
| `-p PREFIX` | Sequence ID prefix |
| `-a AGP` | Input AGP to lift over (legacy coordinate mapping) |
| `-g GFF3` | GFF3 annotation to translate |
| `-b BED` | BED file to translate |
| `--gap N` | Gap size in bp (default: 1000) |
| `--spacer N` | Spacer between trimmed sequences at overlaps (default: 10) |
| `--skipoverlap` | Skip overlap trimming |
| `--reverse` | Reverse AGP direction (new → old) |
| `-u` | Add unplaced sequences into output |
| `--noagp` | Skip AGP output |
| `--ignoreids` | Ignore output IDs in structure file, use prefix + numbering |

---

## Output files

Written to `{outdir}/HaploMake/` when run via Nextflow, or to the current directory when run directly:

| File | Description |
|------|-------------|
| `{out}.fasta` | New pseudomolecule FASTA |
| `{out}.structure.agp` | AGP structure of the new assembly |
| `{out}.legacy_structure.agp` | Lifted-over input AGP (if `-a`/`--hapmake_agp`) |
| `{out}.loci_to_check.txt` | Regions needing manual review (if overlaps found) |

---

## Examples

```bash
# Gap fill + build new assembly via Nextflow
nextflow run nextflow/gap_fill.nf -profile mamba \
    --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \
    --hapfill_correspondence correspondence.tsv \
    --hapfill_repeats repeats.bed \
    --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \
    --run_haplomake \
    --hapmake_agp previous.agp \
    --hapmake_gff3 annotation.gff3 \
    --out myproject --outdir results

# Split an overassembled contig from an edited AGP
python3 HaploMake.py \
    -f hap1.fasta,hap2.fasta \
    -s hap1_corrected.agp \
    --format AGP \
    -o myproject_corrected \
    -p NEW

# Rebuild from a HaploFill structure block, adding unplaced sequences
python3 HaploMake.py \
    -f hap1.fasta,hap2.fasta,unplaced.fasta \
    -s myproject.structure.block \
    -o myproject_v2 \
    -u
```
