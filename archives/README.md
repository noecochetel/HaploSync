# Archives

This directory contains legacy scripts and documentation from the original HaploSync v1.0 repository that are not part of the current Nextflow-based workflows.

## Contents

### `scripts/`

Legacy Python tools not used by the current Nextflow pipelines:

| Script | Description |
|--------|-------------|
| `HaploBreak.py` | Correction of chimeric scaffolds |
| `HaploMap.py` | Pairwise haplotype comparison and mapping |
| `GFF_extract_features.py` | GFF feature extraction utility |
| `StructComp.py` | Structure comparison utility |
| `assembly_statistics.py` | Assembly statistics generator |
| `find_global_alignment.py` | Global alignment finder |

### `manual/`

Original step-by-step documentation for the v1.0 command-line workflow. Superseded by the Nextflow workflow documentation in `nextflow/docs/`.

### `HaploSync.conf.toml_empty`

Configuration template for the legacy TOML-based configuration system. Replaced by Nextflow params files (`nextflow/params_reconstruct_pm.yml`, `nextflow/params_gap_fill.yml`).

### `install.sh`

Legacy installation script. The current pipeline uses a conda environment defined in `nextflow/envs/haplosync.yml`, activated automatically via `-profile conda` or `-profile mamba`.
