# HaploSync: Output

## Introduction

This document describes the output produced by each `--step`. All paths are relative to `--outdir` (default `results`).

## `--step reconstruct_pm`

### `HaploSplit/`

- `{out}.1.fasta` / `{out}.2.fasta` — Hap1 / Hap2 pseudomolecule FASTA
- `{out}.Un.fasta` — unplaced sequences
- `{out}.1.agp` / `{out}.2.agp` / `{out}.Un.agp` — AGP structure of each pseudomolecule set
- `{out}.1.list` / `{out}.2.list` / `{out}.Un.list` — tiling paths (`ChrID<TAB>seq_id|strand,...`)
- `{out}.correspondence.tsv` — chromosome ID → Hap1 sequence ID → Hap2 sequence ID
- `{out}.markers.bed` — marker positions translated into pseudomolecule space (when `--markers`/`--markers_map` given)
- `{out}.legacy_structure.agp` — legacy coordinate mapping (when `--input_agp` given)
- `{out}.unused_sequences.list` — chromosome-assigned but unplaced sequences
- `{out}.conflicting_pseudomolecules.txt`, `{out}.unplaced_to_pseudomolecule.txt` — diagnostic reports

### `HaploDup/` (with `--run_haplodup`)

- `{out}.HaploDup_dir/` — nucmer `.delta` alignments, dotplots, per-chromosome-pair comparisons
- `{out}.HaploDup_dir/index.html` — summary HTML report (browse this for a visual overview)

## `--step gap_fill`

### `HaploFill/`

- `{out}.structure.block` — gap-fill instructions (input to HaploMake): `>SeqID` header + `SeqID\tstart\tstop\tstrand` rows per component
- `{out}.gap_filling_findings.txt` — per-gap diagnostic report (detected gap, mate coordinates, chosen filler strategy and source)

### `HaploMake/` (with `--run_haplomake`, implied by `--run_haplodup`)

- `{out}.fasta` — new gap-filled pseudomolecule FASTA
- `{out}.structure.agp` — AGP for the new assembly
- `{out}.legacy_structure.agp` — legacy coordinate mapping (when `--hapmake_agp` given)
- `{out}.contigs.legacy_structure.agp` — deeper legacy coordinate mapping (when `--hapmake_legacy_agp` given), e.g. tracing back to the original pre-HaploSplit contigs, produced by a second `--noprint`-only HaploMake pass
- `{out}.bed` — marker positions translated onto the new assembly (when `--hapmake_bed` given)
- `{out}.annotation.gff3` — annotation translated onto the new assembly (when `--hapmake_gff3` given)
- `{out}.dropped_loci.txt`, `{out}.multiple_copy_loci.txt` — diagnostic reports from the annotation translation (when `--hapmake_gff3` given)

### `HaploDup/` (with `--run_haplodup`)

Same layout as the reconstruct_pm stage, run against the new gap-filled assembly.

## Pipeline information

- `pipeline_info/execution_report_*.html` — Nextflow execution report (timing, resource usage per process)
- `pipeline_info/execution_timeline_*.html` — Gantt-chart-style timeline
- `pipeline_info/execution_trace_*.txt` — per-task trace (exit codes, duration, memory)
- `pipeline_info/pipeline_dag_*.html` — task dependency graph

## Every module's `versions.yml`

Each Nextflow process also emits a `versions.yml` capturing the exact tool versions it ran with (Python, and whichever of minimap2/samtools/bedtools/MUMmer4/GMAP that specific step invokes) — useful for reproducibility when comparing runs across machines or Nextflow versions.
