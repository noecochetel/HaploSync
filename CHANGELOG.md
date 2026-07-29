# nf-core/haplosync: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Major rewrite of HaploSync as an nf-core-style Nextflow DSL2 pipeline, replacing the standalone Python 2 script collection from v1.0.

### `Added`

- Nextflow DSL2 pipeline scaffold restructured to nf-core conventions, with a unified `--step` entry point and standalone `haplosplit.nf`, `haplodup.nf`, `haplomake.nf` entry points for running individual tools.
- `HAPLOSYNC_GAP_FILL` pipeline (HaploFill → HaploMake → HaploDup) for closing gaps using an alternative haplotype's sequence, including coverage-based reliability signal, contig-based gap-fill test scenarios, and synthetic gene/marker annotation fixtures.
- `HAPLOSYNC_RECONSTRUCT_PM` pipeline (HAPLOSPLIT, QC, HAPLODUP sub-workflows) plus a standalone QC entry point for rejected-sequence structure comparison reports.
- nf-test based CI (`nf-test.yml`), nf-core linting workflow, and status badges.
- Assembly curation guide, Nextflow tips, and workflow diagrams added to the documentation.

### `Changed`

- Entire codebase converted from Python 2 to Python 3 (`2to3` conversion plus manual fixes for gzip/JSON binary-vs-text mode differences).
- HaploDup and HaploSplit split into discrete Nextflow modules (e.g. `HAPLODUP_ALIGN`, `HAPLODUP_GMAP`, `HAPLODUP_REPORT`; `BUILD_PATHS`, `RECONSTRUCT`, `TRANSLATE`) instead of monolithic scripts.
- Chromosome-pair and whole-genome dotplot reports reworked into parametrized Rmd templates rendered in parallel, replacing the large set of hand-duplicated per-variant report scripts.
- HaploMake's overlap-dodging logic retired in favor of a simpler, more reliable unplaced-sequence handling path.
- README rewritten and reorganized around the Nextflow pipeline; legacy v1.0 Python tools and docs moved to `archives/`.

### `Fixed`

- Gap-filling coverage bug and vectorized hot paths in HaploFill STEP 3.2/5.
- Silent nucmer failures and unwired thread count during HaploFill gap-filling; nucmer forced single-threaded to avoid a known mummer4 reliability issue.
- HaploDup input file name collisions and phantom/missing output declarations in the gap-fill pipeline.
- Parallel-render race condition in `make_pair_html_report` and `make_no_genes_html_report`.
- Marker chromosome filter for legacy components in `report_marker_usage`, and contig color coding in chromosome-pair reports when `--rejected` is not used.
- Various Nextflow pipeline correctness issues: correspondence file guards, AGP collisions, combined AGP handling for HaploDup, and config/script parsing under newer Nextflow releases.
- Removed a broken PDF-report code path (`make_pair_pdf_report`, `make_no_genes_pdf_report`) that referenced R scripts which no longer existed under those names; the HTML report remains the supported report format.

### `Dependencies`

- Added mosdepth as a declared conda dependency for coverage calculation.
- Bumped gmap and pinned mummer4/nucmer versions for Linux compatibility.
