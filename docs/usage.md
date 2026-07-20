# HaploSync: Usage

## Introduction

HaploSync runs as two mandatory, separately-invoked stages selected with `--step`, because HaploSplit's tiling-path output almost always needs manual curation (fixing chimeras, adjusting paths, editing AGP coordinates) before it's safe to gap-fill, and Nextflow can't pause mid-run for that:

```bash
# 1. Reconstruct pseudomolecules
nextflow run . --step reconstruct_pm -profile <mamba/conda> \
    --input_fasta assembly.fasta \
    --markers markers.bed --markers_map genetic_map.tsv \
    --outdir results

# -> inspect/edit results/HaploSplit/* by hand (see the curation guide below)

# 2. Fill gaps using the (possibly hand-edited) reconstruct_pm output
nextflow run . --step gap_fill -profile <mamba/conda> \
    --hapfill_hap1 results/HaploSplit/out.1.fasta \
    --hapfill_hap2 results/HaploSplit/out.2.fasta \
    --hapfill_correspondence results/HaploSplit/out.correspondence.tsv \
    --hapfill_repeats repeats.bed \
    --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \
    --outdir results
```

See the [Assembly curation guide](../nextflow/docs/curation_guide.md) for the full iterative loop between the two stages, and [reconstruct_pm.md](../nextflow/docs/reconstruct_pm.md) / [gap_fill.md](../nextflow/docs/gap_fill.md) for the complete parameter reference of each stage (also available via `nextflow run . --step <step> --help`).

## `--step` values

| Value | Default | Description |
| --- | --- | --- |
| `reconstruct_pm` | yes | HaploSplit (tiling-path selection) → QC → optional HaploDup |
| `gap_fill` | no | HaploFill (gap patching) → optional HaploMake → optional HaploDup |

## Advanced: standalone reruns

Beyond the two `--step` values, several named entry points (`-entry <NAME>`) support rerunning a single piece of either stage without repeating everything before it — useful after further manual edits:

| `-entry` | Rereads from | Purpose |
| --- | --- | --- |
| `QC` | `{outdir}/HaploSplit/` | Rerun chromosome-pair / unplaced-sequence QC only |
| `RECONSTRUCT_PM_HAPLODUP` | `{outdir}/HaploSplit/` | Rerun HaploDup on the reconstruct_pm output only |
| `HAPLOMAKE` | `{outdir}/HaploFill/` | Rerun HaploMake on an existing HaploFill structure block |
| `GAPFILL_HAPLODUP` | `{outdir}/HaploMake/` | Rerun HaploDup on the gap_fill output only |
| `HAPLODUP_GENERIC` | explicit `--hap1_fasta`/`--hap2_fasta` | Run HaploDup on any pair of haplotype FASTAs |
| `HAPLOMAKE_GENERIC` | explicit `--fasta`/`--structure_block` | Run HaploMake from any BLOCK/AGP/BED structure file |

Each supports its own `--help`, e.g. `nextflow run . -entry HAPLODUP_GENERIC --help`.

## Running the bundled test genome

A tiny (~10 kb) synthetic diploid genome lives under `test_data/`, with a deliberately-injected gap in one haplotype so it exercises HaploFill's gap-filling logic. It's used both by `nf-core pipelines lint`/CI and for local smoke-testing:

```bash
nextflow run . -profile test,<mamba/conda> --outdir results_pm
nextflow run . -profile test_gapfill,<mamba/conda> --outdir results_gf
```

`-profile test_gapfill` runs `--step gap_fill` against a **static fixture** (the already-verified output of `-profile test`, committed under `test_data/gap_fill_fixtures/`) rather than chaining a live run — gap_fill always follows a manual-curation step that CI can't perform, so this mirrors how [nf-core/sarek](https://nf-co.re/sarek) tests its own staged `--step` reruns: each stage is tested independently against a fixed, known-good input.

The same two scenarios are also covered by the `nf-test` suite in `tests/`:

```bash
nf-test test --profile <conda/mamba>
```

## Updating tool citations

Software dependency citations live in [CITATIONS.md](../CITATIONS.md), used by the pipeline's completion summary.

## Core parameters

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Multiple profiles can be loaded, comma-separated, e.g. `-profile test,mamba`.

- `standard` — local execution, no environment management (tools must already be on `$PATH`)
- `conda` / `mamba` — local execution, Nextflow creates the conda environment automatically
- `hpc` — SLURM executor with mamba environment management
- `test` / `test_gapfill` / `test_full` — the bundled tiny test genome (see above)

### `-resume`

Restart a pipeline using cached results from a previous run. Nextflow only re-executes processes whose inputs changed, letting you resume after a manual-curation edit or a failure without rerunning everything.

```bash
nextflow run . --step reconstruct_pm -resume -params-file params.yml
```

### `-c`

Specify the path to a custom config file, e.g. to cap resource labels for a laptop:

```bash
nextflow run . -c my.config ...
```

## Running in the background

Add `-bg` to a `nextflow run` command to launch it in the background, or use `screen`/`tmux` for long HPC runs.

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory. If you set `NXF_OPTS` in your environment (e.g. `export NXF_OPTS='-Xms1g -Xmx4g'`) it can be used to override this.
