#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters } from 'plugin/nf-schema'
include { HAPLOSYNC_RECONSTRUCT_PM } from './workflows/haplosync'
include { HAPLOSYNC_GAP_FILL       } from './workflows/haplosync'
include { ALIGN    as HD_ALIGN     } from './nextflow/modules/local/haplodup_align/main'
include { GMAP     as HD_GMAP      } from './nextflow/modules/local/haplodup_gmap/main'
include { REPORT   as HD_REPORT    } from './nextflow/modules/local/haplodup_report/main'
include { CHR_PAIR as QC_CHR_PAIR  } from './nextflow/modules/local/chr_pair_qc/main'
include { REJECTED as QC_REJECTED  } from './nextflow/modules/local/rejected_qc/main'
include { HM_MAKE                  } from './nextflow/modules/local/hapmake/main'

// --------------------------------------------------------------------------
// Typed parameters
//   Only the non-string parameters are declared here, because the strict syntax
//   parser converts command-line values to the declared type: without this a
//   bare `--run_haplodup` or `--cores 8` arrives as the string "true" or "8" and
//   fails schema validation. Defaults live in nextflow.config (checked against
//   nextflow_schema.json by nf-core lint), so no defaults are repeated here.
//   When adding a boolean/integer/number parameter, declare it in all three.
// --------------------------------------------------------------------------
params {
    // Reconstruct PM options
    run_alignment: Boolean
    hitgap: Integer
    distance1: Integer
    distance2: Integer
    reuse_intermediate: Boolean
    force_direction1: Boolean
    force_direction2: Boolean
    minR1: Integer?
    minR2: Integer?
    No2: Boolean
    gap: Integer
    conc: Integer?
    filter_hits: Boolean
    extended_region: Boolean
    allow_rearrangements: Boolean
    required_as_path: Boolean
    skip_chimeric_qc: Boolean
    disable_marker_ploidy_check: Boolean
    only_markers: Boolean
    avoid_rejected_qc: Boolean
    skip_chr_pair_reports: Boolean
    skip_unplaced_qc: Boolean

    // HaploDup options
    run_haplodup: Boolean
    hit_identity: Integer
    hit_coverage: Integer
    gene_identity: Integer
    gene_coverage: Integer
    unbalanced_ratio: Float
    haplodup_window: Integer
    haplodup_allowed: Integer
    reuse_mappings: Boolean
    reuse_dotplots: Boolean
    reuse_gmap: Boolean
    skip_dotplots_by_chr: Boolean
    only_paired_dotplots: Boolean

    // Gap Fill options
    hapfill_coverage: Integer?
    hapfill_flanking: Integer?
    hapfill_map_threads: Integer
    hapfill_nohomozygous: Boolean
    hapfill_overwrite: Boolean

    // HaploMake options
    run_haplomake: Boolean
    hapmake_gap: Integer
    hapmake_noagp: Boolean

    // Standalone HaploMake / HaploDup options
    hapmake_reverse: Boolean

    // Generic options
    help: Boolean
    cores: Integer
    monochrome_logs: Boolean
}

// --------------------------------------------------------------------------
// Help message
// --------------------------------------------------------------------------
def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                         HaploSync v2.0                              ║
    ║          Haplotype-resolved pseudomolecule assembly                 ║
    ╚══════════════════════════════════════════════════════════════════════╝

    HaploSync runs as two mandatory, separately-invoked stages, because
    HaploSplit's tiling-path output almost always needs manual curation
    (fixing chimeras, adjusting paths, editing AGP coordinates) before it's
    safe to gap-fill. Nextflow can't pause mid-run for that, so the two
    stages are selected with --step and chained by hand:

        1. nextflow run . --step reconstruct_pm ...
           -> inspect/edit results/HaploSplit/* by hand
        2. nextflow run . --step gap_fill --hapfill_hap1 <possibly-edited path> ...

    See docs/usage.md for the full manual-curation workflow.

    Usage:
        nextflow run . --step reconstruct_pm [options]
        nextflow run . --step gap_fill [options]
        nextflow run . -params-file params.yml

    ── Steps ───────────────────────────────────────────────────────────────
        --step reconstruct_pm   [default] HaploSplit -> QC -> optional HaploDup
        --step gap_fill                   HaploFill -> optional HaploMake -> optional HaploDup

        For step-specific options: nextflow run . --step <step> --help

    ── Advanced / manual reruns (also selected with --step) ────────────────
        --step qc                        Rerun QC only, reading an existing HaploSplit output
        --step reconstruct_pm_haplodup   Rerun HaploDup only, reading an existing HaploSplit output
        --step haplomake                 Rerun HaploMake only, reading an existing HaploFill output
        --step gapfill_haplodup          Rerun HaploDup only, reading an existing HaploMake output
        --step haplodup_generic          Run HaploDup directly on any pair of haplotype FASTAs
        --step haplomake_generic         Run HaploMake directly from any structure file (BLOCK/AGP/BED)

        Each supports --help for its own options, e.g.:
        nextflow run . --step haplodup_generic --help

    ── Output ──────────────────────────────────────────────────────────────
        --out               Output files prefix          [default: out]
        --outdir            Results directory            [default: results]

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── Profiles ────────────────────────────────────────────────────────────
        -profile standard   Local execution
        -profile conda      Local execution with conda
        -profile mamba      Local execution with mamba/micromamba
        -profile test          CI/smoke test: --step reconstruct_pm on the bundled tiny genome
        -profile test_gapfill  CI/smoke test: --step gap_fill on the bundled fixture output
    """.stripIndent()
}

def helpReconstructPm() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║              HaploSync — --step reconstruct_pm                      ║
    ║          Haplotype-resolved pseudomolecule assembly                 ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run . --step reconstruct_pm [options]
        nextflow run . --step reconstruct_pm -params-file params.yml

    ── Required ────────────────────────────────────────────────────────────
        --input_fasta       Draft assembly FASTA
        Provide (--markers AND --markers_map) and/or --guide_genome

    ── Genetic map ─────────────────────────────────────────────────────────
        --markers           Marker hits on input sequences (BED, 4 columns)
        --markers_map       Marker genetic map (chr, position, marker_id)

    ── Reference genome ────────────────────────────────────────────────────
        --guide_genome      Guide/reference genome FASTA
        --local_alignment   Existing PAF file, skips minimap2 (-l)
        --run_alignment     Run minimap2 to generate alignment (--align)
        --mapping_tool      Minimap2 options string      [default: '--cs -x asm20 -r 1000']
        --distance1         Max Hap1 sequence gap (bp)   [default: 2,000,000]
        --distance2         Max Hap2 sequence gap (bp)   [default: 4,000,000]
        --hitgap            Max gap to merge hits (bp)   [default: 100,000]
        --reuse_intermediate  Reuse existing nucmer files [default: false]

    ── Optional inputs ─────────────────────────────────────────────────────
        --input_agp         Input AGP structure file
        --gff3              Gene annotation GFF3
        --exclusion         Mutually exclusive sequence pairs (TSV)
        --known             Sequences known to be in same haplotype (TSV)
        --alternative_groups  Alternative haplotype sequence pairs (TSV)
        --Require1          Sequences required in Hap1 (TSV)
        --Require2          Sequences required in Hap2 (TSV)
        --Blacklist1        Blacklisted sequences for Hap1 (TSV)
        --Blacklist2        Blacklisted sequences for Hap2 (TSV)
        --input_groups      Sequence grouping file
        --legacy_groups     Legacy component group file
        --path1             Pre-defined Hap1 tiling path, skips DAG (-1)
        --path2             Pre-defined Hap2 tiling path, skips DAG (-2)

    ── Assembly behaviour ───────────────────────────────────────────────────
        --No2               Skip Hap2 reconstruction     [default: false]
        --gap               Gap size in bp               [default: 1000]
        --filter_hits       Remove seqs with intra-seq marker duplication
        --extended_region   Extend seq association to all chr markers
        --conflict_resolution  Conflict resolution: exit|ignore|release
        --allow_rearrangements  Allow rearrangements between haplotypes
        --required_as_path  Treat --R1/--R2 as full tiling paths

    ── QC ──────────────────────────────────────────────────────────────────
        --skip_chimeric_qc  Skip chimeric sequence QC    [default: false]
        --disable_marker_ploidy_check
        --only_markers      Limit rejected-seq QC to marker-bearing seqs
        --avoid_rejected_qc  Skip QC of chr-assigned but unplaced seqs
        --skip_chr_pair_reports  Skip chr pair overview reports [default: false]
        --skip_unplaced_qc  Skip unplaced sequence QC reports  [default: false]

    ── Output ──────────────────────────────────────────────────────────────
        --out               Output files prefix          [default: out]
        --prefix            Sequence ID prefix           [default: NEW]
        --outdir            Results directory            [default: results]

    ── HaploDup ────────────────────────────────────────────────────────────
        --run_haplodup      Run HaploDup after reconstruction [default: false]
        --haplodup_opts     Extra HaploDup flags as a quoted string

        For full HaploDup options: nextflow run . --step reconstruct_pm_haplodup --help

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── Examples ────────────────────────────────────────────────────────────
        # Genetic map mode
        nextflow run . -profile mamba --step reconstruct_pm \\
            --input_fasta assembly.fasta \\
            --markers markers.bed \\
            --markers_map genetic_map.tsv \\
            --out myproject --outdir results

        # Full pipeline with HaploDup
        nextflow run . -profile mamba --step reconstruct_pm \\
            --input_fasta assembly.fasta \\
            --markers markers.bed --markers_map genetic_map.tsv \\
            --run_haplodup \\
            --out myproject --outdir results

        # Resume after manual curation
        nextflow run . -profile mamba --step reconstruct_pm -resume -params-file params.yml
    """.stripIndent()
}

def helpGapFill() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                   HaploSync — --step gap_fill                       ║
    ║              Gap-filling and pseudomolecule patching                ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run . --step gap_fill [options]
        nextflow run . --step gap_fill -params-file params.yml

    Pipeline: HaploFill -> HaploMake -> HaploDup (optional)

    ── Required ────────────────────────────────────────────────────────────
        --hapfill_hap1          Hap1 FASTA (pseudomolecules to patch — normally
                                 the manually-curated output of --step reconstruct_pm)
        --hapfill_hap2          Hap2 FASTA (pseudomolecules to patch — normally
                                 the manually-curated output of --step reconstruct_pm)
        --hapfill_correspondence  Chromosome correspondence TSV
        --hapfill_repeats       Repeats BED file
        --hapfill_b1            BAM file aligned to Hap1
        --hapfill_b2            BAM file aligned to Hap2

    ── Optional inputs ─────────────────────────────────────────────────────
        --hapfill_unplaced      Unplaced sequences FASTA
        --hapfill_coverage      Per-base coverage threshold  [default: auto]
        --hapfill_flanking      Flanking region size (bp)    [default: auto]
        --hapfill_map_threads   Threads for minimap2 alignment [default: 4]
        --hapfill_nohomozygous  Skip homozygous gap filling  [default: false]

    ── Coverage tool ────────────────────────────────────────────────────────
        --coverage_tool         Coverage backend: bedtools|mosdepth [default: bedtools]

    ── HaploMake options ────────────────────────────────────────────────────
        --hapmake_prefix        Sequence ID prefix
        --hapmake_agp           AGP reference for legacy coordinate mapping
        --hapmake_gff3          Gene annotation GFF3 to translate
        --hapmake_bed           BED file to translate
        --hapmake_gap           Gap size in bp             [default: 1000]
        --hapmake_noagp         Skip AGP output            [default: false]
        --hapmake_legacy_agp    Deeper legacy AGP (e.g. pre-HaploSplit contigs),
                                 ported via a second, --noprint-only HaploMake pass

    ── HaploMake (optional) ─────────────────────────────────────────────────
        --run_haplomake         Build new FASTA/AGP from gap-fill result [default: false]

    ── HaploDup (optional) ──────────────────────────────────────────────────
        --run_haplodup          Run HaploDup on the gap-filled assembly  [default: false]
                                 Implies --run_haplomake

        For full HaploDup options: nextflow run . --step gapfill_haplodup --help

    ── Output ──────────────────────────────────────────────────────────────
        --out               Output files prefix          [default: out]
        --outdir            Results directory            [default: results]

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── RAM note ─────────────────────────────────────────────────────────────
        HF_COVERAGE runs one job per chromosome in parallel.
        Cap concurrency in nextflow.config to avoid RAM exhaustion:
          process { withName: 'HF_COVERAGE' { maxForks = 8 } }   // ~48 GB RAM

    ── Examples ────────────────────────────────────────────────────────────
        # Basic gap-fill run
        nextflow run . -profile mamba --step gap_fill \\
            --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \\
            --hapfill_correspondence correspondence.tsv \\
            --hapfill_repeats repeats.bed \\
            --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \\
            --out myproject --outdir results

        # Gap fill + build new assembly + HaploDup QC
        nextflow run . -profile mamba --step gap_fill \\
            --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \\
            --hapfill_correspondence correspondence.tsv \\
            --hapfill_repeats repeats.bed \\
            --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \\
            --run_haplodup \\
            --out myproject --outdir results
    """.stripIndent()
}

def helpQc() {
    log.info """
    Usage:
        nextflow run . --step qc [options]

    Reads HaploSplit outputs automatically from --outdir/HaploSplit/ using
    the --out prefix. Run --step reconstruct_pm first.

    ── Required ────────────────────────────────────────────────────────────
        --out               Output prefix (must match HaploSplit run) [default: out]
        --outdir            Results directory (must match HaploSplit run) [default: results]

    ── QC selection ────────────────────────────────────────────────────────
        --skip_chr_pair_reports  Skip chromosome pair overview reports
        --skip_unplaced_qc       Skip unplaced sequence QC reports

    ── Optional inputs ─────────────────────────────────────────────────────
        --markers_map       Marker genetic map (chr, position, marker_id)
        --input_groups      Sequence grouping file
        --legacy_groups     Legacy component group file
    """.stripIndent()
}

def helpReconstructPmHaplodup() {
    log.info """
    Usage:
        nextflow run . --step reconstruct_pm_haplodup [options]

    HaploDup reads HaploSplit outputs automatically from --outdir/HaploSplit/
    using the --out prefix. Run --step reconstruct_pm first.

    ── Required ────────────────────────────────────────────────────────────
        --out               Output prefix (must match HaploSplit run) [default: out]
        --outdir            Results directory (must match HaploSplit run) [default: results]

    ── Optional inputs ─────────────────────────────────────────────────────
        --reference             Reference genome for dotplots
        --markers_map           Marker genetic map for QC
        --input_groups          Sequence grouping file
        --legacy_groups         Legacy grouping file
        --functional_annotation Functional annotation per transcript

    ── Alignment thresholds ─────────────────────────────────────────────────
        --hit_identity          Min genome mapping hit identity  [default: 90]
        --hit_coverage          Min genome mapping hit length    [default: 3000]
        --gene_identity         Min gene mapping identity        [default: 95]
        --gene_coverage         Min gene mapping coverage        [default: 95]
        --unbalanced_ratio      Gene count ratio threshold       [default: 0.33]

    ── Gene mapping ─────────────────────────────────────────────────────────
        --haplodup_feature      GFF feature type [CDS|mRNA]     [default: CDS]
        --haplodup_window       Window size for unbalanced gene search [default: 10]
        --haplodup_allowed      Allowed unbalanced genes/window  [default: 5]

    ── Reuse / skip ─────────────────────────────────────────────────────────
        --reuse_mappings        Reuse existing genome alignments [default: false]
        --reuse_dotplots        Reuse existing dotplots          [default: false]
        --reuse_gmap            Reuse existing GMAP mappings     [default: false]
        --skip_dotplots_by_chr  Skip individual chr-vs-chr dotplots [default: false]
        --only_paired_dotplots  Only generate matched-pair dotplots [default: false]
    """.stripIndent()
}

def helpHaplomake() {
    log.info """
    Usage:
        nextflow run . --step haplomake [options]

    Reads the structure block from {outdir}/HaploFill/{out}.structure.block
    unless --structure_block is provided.

    ── Required ────────────────────────────────────────────────────────────
        --hapfill_hap1          Hap1 FASTA
        --hapfill_hap2          Hap2 FASTA
        --out               Output prefix   [default: out]
        --outdir            Results directory [default: results]

    ── Optional ────────────────────────────────────────────────────────────
        --hapfill_unplaced      Unplaced sequences FASTA
        --structure_block       Override path to .structure.block file
        --hapmake_prefix        Sequence ID prefix
        --hapmake_agp           AGP to lift over
        --hapmake_gff3          GFF3 annotation to translate
        --hapmake_bed           BED file to translate
        --hapmake_gap           Gap size in bp [default: 1000]
        --hapmake_noagp         Skip AGP output
    """.stripIndent()
}

def helpGapfillHaplodup() {
    log.info """
    Usage:
        nextflow run . --step gapfill_haplodup [options]

    HaploDup reads HaploMake outputs from --outdir/HaploMake/ using the
    --out prefix. Run --step gap_fill (with --run_haplomake) first.

    ── Required ────────────────────────────────────────────────────────────
        --hapfill_hap1          Original Hap1 FASTA (for correspondence)
        --hapfill_hap2          Original Hap2 FASTA (for correspondence)
        --hapfill_correspondence  Chromosome correspondence TSV
        --out               Output prefix (must match gap-fill run) [default: out]
        --outdir            Results directory (must match gap-fill run) [default: results]

    ── Optional inputs ─────────────────────────────────────────────────────
        --hapfill_unplaced      Unplaced sequences FASTA
        --gff3                  Gene annotation GFF3
    """.stripIndent()
}

def helpHaplodupGeneric() {
    log.info """
    Usage:
        nextflow run . --step haplodup_generic [options]

    Runs all-vs-all nucmer alignments between haplotypes, optionally maps
    gene models with GMAP, and generates HTML/PDF reports. Unlike the
    context-aware HaploDup reruns (--step reconstruct_pm_haplodup and
    --step gapfill_haplodup), this takes explicit file paths rather than
    discovering them from --outdir.

    ── Required ────────────────────────────────────────────────────────────
        --hap1_fasta        Hap1 pseudomolecule FASTA
        --hap2_fasta        Hap2 pseudomolecule FASTA
        --correspondence    Chromosome correspondence TSV

    ── Optional inputs ─────────────────────────────────────────────────────
        --unplaced_fasta    Unplaced sequences FASTA
        --agp               AGP file(s), comma-separated
        --gff3              Gene annotation GFF3 (enables GMAP mapping)
        --markers_bed       Marker positions BED
        --legacy_agp        Legacy AGP structure
        --reference         Reference genome FASTA (for additional dotplots)
        --markers_map       Marker genetic map
        --functional_annotation  Functional annotation per transcript
        --input_groups      Sequence grouping file
        --legacy_groups     Legacy component group file

    ── Alignment thresholds ─────────────────────────────────────────────────
        --hit_identity      Min genome alignment identity  [default: 90]
        --hit_coverage      Min genome alignment length    [default: 3000]
        --gene_identity     Min gene mapping identity      [default: 95]
        --gene_coverage     Min gene mapping coverage      [default: 95]
        --unbalanced_ratio  Gene count imbalance threshold [default: 0.33]
    """.stripIndent()
}

def helpHaplomakeGeneric() {
    log.info """
    Usage:
        nextflow run . --step haplomake_generic [options]

    Constructs new pseudomolecule FASTA and AGP from a structure file.
    The structure file can be a HaploFill block, an AGP, or a BED file.

    ── Required ────────────────────────────────────────────────────────────
        --fasta             Input FASTA file(s), comma-separated if multiple
        --structure_block   Structure file (.structure.block, .agp, or .bed)

    ── HaploMake options ────────────────────────────────────────────────────
        --hapmake_format        Structure file format [BLOCK|AGP|BED] [default: BLOCK]
        --hapmake_prefix        Sequence ID prefix (only used with BED input; with
                                 BLOCK/AGP the object names come from the structure file)
        --hapmake_agp           AGP to lift over into new assembly space
        --hapmake_gff3          GFF3 annotation to translate
        --hapmake_bed           BED file to translate
        --hapmake_gap           Gap size in bp             [default: 1000]
        --hapmake_noagp         Skip AGP output            [default: false]
        --hapmake_reverse       Reverse AGP direction (new -> old) [default: false]

    ── Output ──────────────────────────────────────────────────────────────
        Published to --outdir/HaploMake/ (files marked * only if requested):
        {out}.fasta, {out}.structure.agp, {out}.legacy_structure.agp *,
        {out}.annotation.gff3 *, {out}.bed *, {out}.dropped_loci.txt *,
        {out}.multiple_copy_loci.txt *
    """.stripIndent()
}

// Print the help text of the selected --step; unknown steps get the overview.
def helpForStep(step) {
    if      (step == 'reconstruct_pm')         { helpReconstructPm() }
    else if (step == 'gap_fill')               { helpGapFill() }
    else if (step == 'qc')                     { helpQc() }
    else if (step == 'reconstruct_pm_haplodup') { helpReconstructPmHaplodup() }
    else if (step == 'haplomake')              { helpHaplomake() }
    else if (step == 'gapfill_haplodup')       { helpGapfillHaplodup() }
    else if (step == 'haplodup_generic')       { helpHaplodupGeneric() }
    else if (step == 'haplomake_generic')      { helpHaplomakeGeneric() }
    else                                       { helpMessage() }
}

// --------------------------------------------------------------------------
// Default entry point — dispatches on --step
//   Log names: HAPLOSYNC_RECONSTRUCT_PM:HAPLOSPLIT:<PROCESS>, etc.
//              HAPLOSYNC_GAP_FILL:HAPLOFILL:<PROCESS>, etc.
//   The advanced reruns (QC, *_HAPLODUP, HAPLOMAKE*) are named workflows
//   defined below and also selected with --step: the strict syntax parser
//   (default since Nextflow 26.04) rejects -entry.
// --------------------------------------------------------------------------
workflow {

    // Help must run before validateParameters(): under the strict parser a bare
    // --help reaches the schema check as the string 'true' and fails validation.
    if (params.help) {
        helpForStep(params.step)
        exit 0
    }

    validateParameters()

    def wf = workflow
    def outdir = params.outdir
    workflow.onComplete {
        log.info (wf.success
            ? "\n[HaploSync] Pipeline completed successfully.\n  Results: ${outdir}"
            : "\n[HaploSync] Pipeline failed. Check logs for details.")
    }

    if (params.step == 'reconstruct_pm') {
        if (!params.input_fasta) {
            log.error "[ERROR] --input_fasta is required"
            helpReconstructPm()
            exit 1
        }
        def has_markers   = params.markers && params.markers_map
        def has_reference = params.guide_genome
        if (!has_markers && !has_reference) {
            log.error "[ERROR] Provide (--markers AND --markers_map) and/or --guide_genome"
            helpReconstructPm()
            exit 1
        }
        if ((params.markers && !params.markers_map) || (!params.markers && params.markers_map)) {
            log.error "[ERROR] --markers and --markers_map must be provided together"
            helpReconstructPm()
            exit 1
        }

        HAPLOSYNC_RECONSTRUCT_PM()

    } else if (params.step == 'gap_fill') {
        def required = [
            hapfill_hap1:            '--hapfill_hap1',
            hapfill_hap2:            '--hapfill_hap2',
            hapfill_correspondence:  '--hapfill_correspondence',
            hapfill_repeats:         '--hapfill_repeats',
            hapfill_b1:              '--hapfill_b1',
            hapfill_b2:              '--hapfill_b2'
        ]
        required.each { param, flag ->
            if (!params[param]) {
                log.error "[ERROR] ${flag} is required"
                helpGapFill()
                exit 1
            }
        }

        HAPLOSYNC_GAP_FILL()

    } else if (params.step == 'qc') {
        QC()

    } else if (params.step == 'reconstruct_pm_haplodup') {
        RECONSTRUCT_PM_HAPLODUP()

    } else if (params.step == 'haplomake') {
        HAPLOMAKE()

    } else if (params.step == 'gapfill_haplodup') {
        GAPFILL_HAPLODUP()

    } else if (params.step == 'haplodup_generic') {
        HAPLODUP_GENERIC()

    } else if (params.step == 'haplomake_generic') {
        HAPLOMAKE_GENERIC()

    } else {
        log.error "[ERROR] Unknown --step '${params.step}'. Valid values: reconstruct_pm, gap_fill, qc, reconstruct_pm_haplodup, haplomake, gapfill_haplodup, haplodup_generic, haplomake_generic"
        helpMessage()
        exit 1
    }
}

// --------------------------------------------------------------------------
// Named workflow: QC (--step qc, standalone rerun)
//   Reads HaploSplit outputs from --outdir/HaploSplit/.
// --------------------------------------------------------------------------
workflow QC {

    def hs_dir = "${params.outdir}/HaploSplit"
    def pfx    = "${hs_dir}/${params.out}"

    def required_files = [
        file("${pfx}.correspondence.tsv"),
        file("${pfx}.1.agp"),
        file("${pfx}.Un.agp")
    ]
    required_files.each { f ->
        if (!f.exists()) {
            log.error "[ERROR] Required HaploSplit output not found: ${f}\n         Run --step reconstruct_pm first or check --out / --outdir"
            exit 1
        }
    }

    def correspondence = Channel.fromPath("${pfx}.correspondence.tsv")

    def agp_ch = Channel.of(
            file("${pfx}.1.agp"),
            file("${pfx}.2.agp"),
            file("${pfx}.Un.agp")
        )
        .filter { it.exists() }
        .collect()

    def markers_bed_f  = file("${pfx}.markers.bed")
    def legacy_agp_f   = file("${pfx}.legacy_structure.agp")

    def markers_bed_ch = markers_bed_f.exists() ? Channel.value(markers_bed_f) : Channel.value([])
    def legacy_agp_ch  = legacy_agp_f.exists()  ? Channel.value(legacy_agp_f)  : Channel.value([])

    if (!params.skip_chr_pair_reports && !params.No2) {
        QC_CHR_PAIR(correspondence, agp_ch, markers_bed_ch, legacy_agp_ch)
    }

    if (!params.skip_unplaced_qc && !params.No2) {
        def unused_list_f = file("${pfx}.unused_sequences.list")
        if (!unused_list_f.exists()) {
            log.error "[ERROR] Required file not found for unplaced QC: ${unused_list_f}"
            exit 1
        }

        def hap1_fasta = Channel.fromPath("${pfx}.1.fasta")
        def hap2_fasta = file("${pfx}.2.fasta").exists()
                             ? Channel.fromPath("${pfx}.2.fasta")
                             : Channel.empty()
        def un_fasta   = Channel.fromPath("${pfx}.Un.fasta")

        def fasta_ch = hap1_fasta
            .mix(hap2_fasta)
            .mix(un_fasta)
            .collect()

        QC_REJECTED(
            Channel.value(unused_list_f),
            correspondence,
            fasta_ch,
            agp_ch,
            markers_bed_ch,
            legacy_agp_ch
        )
    }
}

// --------------------------------------------------------------------------
// Named workflow: RECONSTRUCT_PM_HAPLODUP (--step reconstruct_pm_haplodup, standalone rerun)
//   Reads HaploSplit outputs from --outdir/HaploSplit/.
// --------------------------------------------------------------------------
workflow RECONSTRUCT_PM_HAPLODUP {

    def hs_dir = "${params.outdir}/HaploSplit"
    def pfx    = "${hs_dir}/${params.out}"

    def required_files = [
        file("${pfx}.1.fasta"),
        file("${pfx}.Un.fasta"),
        file("${pfx}.correspondence.tsv"),
        file("${pfx}.1.agp"),
        file("${pfx}.Un.agp")
    ]
    required_files.each { f ->
        if (!f.exists()) {
            log.error "[ERROR] Required HaploSplit output not found: ${f}\n         Run --step reconstruct_pm first or check --out / --outdir"
            exit 1
        }
    }

    def hap1_fasta     = Channel.fromPath("${pfx}.1.fasta")
    def hap2_fasta     = file("${pfx}.2.fasta").exists()
                             ? Channel.fromPath("${pfx}.2.fasta")
                             : Channel.value([])
    def un_fasta       = Channel.fromPath("${pfx}.Un.fasta")
    def correspondence = Channel.fromPath("${pfx}.correspondence.tsv")

    def agp_ch = Channel.of(
            file("${pfx}.1.agp"),
            file("${pfx}.2.agp"),
            file("${pfx}.Un.agp")
        )
        .filter { it.exists() }
        .collect()

    def markers_bed_f  = file("${pfx}.markers.bed")
    def legacy_agp_f   = file("${pfx}.legacy_structure.agp")
    def annotation_f   = file("${pfx}.annotation.gff3")

    def markers_bed_ch = markers_bed_f.exists() ? Channel.value(markers_bed_f) : Channel.value([])
    def legacy_agp_ch  = legacy_agp_f.exists()  ? Channel.value(legacy_agp_f)  : Channel.value([])
    def annotation_ch  = annotation_f.exists()   ? Channel.value(annotation_f)  : Channel.value([])

    HD_ALIGN(hap1_fasta, hap2_fasta, un_fasta, correspondence)

    def run_gmap     = annotation_f.exists() && !params.No2
    def gmap_gff3_ch = Channel.value([])
    if (run_gmap) {
        HD_GMAP(hap1_fasta, hap2_fasta, un_fasta, correspondence, annotation_ch)
        gmap_gff3_ch = HD_GMAP.out.gmap_gff3
    }

    HD_REPORT(
        hap1_fasta,
        hap2_fasta,
        un_fasta,
        agp_ch,
        correspondence,
        markers_bed_ch,
        legacy_agp_ch,
        annotation_ch,
        HD_ALIGN.out.delta_files.collect(),
        gmap_gff3_ch
    )
}

// --------------------------------------------------------------------------
// Named workflow: HAPLOMAKE (--step haplomake, standalone rerun, gap-fill context)
//   Reads HaploFill structure block from --outdir/HaploFill/.
//   Use --structure_block to override with a custom path.
// --------------------------------------------------------------------------
workflow HAPLOMAKE {

    def block_path = params.structure_block
        ?: "${params.outdir}/HaploFill/${params.out}.structure.block"
    def block_file = file(block_path)
    if (!block_file.exists()) {
        log.error "[ERROR] Structure block not found: ${block_file}\n         Run --step gap_fill first or provide --structure_block"
        exit 1
    }

    def hap1_fasta = Channel.fromPath(params.hapfill_hap1)
    def hap2_fasta = Channel.fromPath(params.hapfill_hap2)
    def un_fasta   = params.hapfill_unplaced
                         ? Channel.fromPath(params.hapfill_unplaced)
                         : Channel.value([])

    HM_MAKE(hap1_fasta, hap2_fasta, un_fasta, Channel.value(block_file))
}

// --------------------------------------------------------------------------
// Named workflow: GAPFILL_HAPLODUP (--step gapfill_haplodup, standalone rerun, gap-fill context)
//   Reads HaploMake outputs from --outdir/HaploMake/.
// --------------------------------------------------------------------------
workflow GAPFILL_HAPLODUP {

    def hm_dir = "${params.outdir}/HaploMake"
    def pfx    = "${hm_dir}/${params.out}"

    def required_files = [ file("${pfx}.fasta") ]
    required_files.each { f ->
        if (!f.exists()) {
            log.error "[ERROR] Required HaploMake output not found: ${f}\n         Run --step gap_fill first or check --out / --outdir"
            exit 1
        }
    }

    def new_fasta      = Channel.fromPath("${pfx}.fasta")
    def un_fasta       = params.hapfill_unplaced
                             ? Channel.fromPath(params.hapfill_unplaced)
                             : Channel.value([])
    def correspondence = Channel.fromPath(params.hapfill_correspondence)

    def agp_ch = Channel.of(file("${pfx}.structure.agp"))
        .filter { it.exists() }
        .collect()

    def legacy_agp_f = file("${pfx}.legacy_structure.agp")
    def legacy_agp_ch = legacy_agp_f.exists() ? Channel.value(legacy_agp_f) : Channel.value([])
    def annotation_ch = (params.gff3 && file(params.gff3).exists())
                            ? Channel.fromPath(params.gff3)
                            : Channel.value([])

    HD_ALIGN(new_fasta, new_fasta, un_fasta, correspondence)

    def run_gmap     = params.gff3 && !params.No2
    def gmap_gff3_ch = Channel.value([])
    if (run_gmap) {
        HD_GMAP(new_fasta, new_fasta, un_fasta, correspondence, annotation_ch)
        gmap_gff3_ch = HD_GMAP.out.gmap_gff3
    }

    HD_REPORT(
        new_fasta,
        new_fasta,
        un_fasta,
        agp_ch,
        correspondence,
        Channel.value([]),
        legacy_agp_ch,
        annotation_ch,
        HD_ALIGN.out.delta_files.collect(),
        gmap_gff3_ch
    )
}

// --------------------------------------------------------------------------
// Named workflow: HAPLODUP_GENERIC (--step haplodup_generic, fully standalone, any haplotype FASTA pair)
// --------------------------------------------------------------------------
workflow HAPLODUP_GENERIC {

    def required = [
        hap1_fasta:      '--hap1_fasta',
        hap2_fasta:      '--hap2_fasta',
        correspondence:  '--correspondence'
    ]
    required.each { param, flag ->
        if (!params[param]) {
            log.error "[ERROR] ${flag} is required"
            exit 1
        }
    }

    def hap1_fasta     = Channel.fromPath(params.hap1_fasta)
    def hap2_fasta     = Channel.fromPath(params.hap2_fasta)
    def un_fasta       = params.unplaced_fasta
                             ? Channel.fromPath(params.unplaced_fasta)
                             : Channel.value([])
    def correspondence = Channel.fromPath(params.correspondence)

    def agp_ch = params.agp
        ? Channel.fromPath(params.agp.tokenize(',')).collect()
        : Channel.value([])

    def markers_bed_ch = params.markers_bed
        ? Channel.fromPath(params.markers_bed)
        : Channel.value([])
    def legacy_agp_ch  = params.legacy_agp
        ? Channel.fromPath(params.legacy_agp)
        : Channel.value([])
    def annotation_ch  = params.gff3
        ? Channel.fromPath(params.gff3)
        : Channel.value([])

    HD_ALIGN(hap1_fasta, hap2_fasta, un_fasta, correspondence)

    def run_gmap     = params.gff3 && !params.No2
    def gmap_gff3_ch = Channel.value([])
    if (run_gmap) {
        HD_GMAP(hap1_fasta, hap2_fasta, un_fasta, correspondence, annotation_ch)
        gmap_gff3_ch = HD_GMAP.out.gmap_gff3
    }

    HD_REPORT(
        hap1_fasta,
        hap2_fasta,
        un_fasta,
        agp_ch,
        correspondence,
        markers_bed_ch,
        legacy_agp_ch,
        annotation_ch,
        HD_ALIGN.out.delta_files.collect(),
        gmap_gff3_ch
    )
}

// --------------------------------------------------------------------------
// Process: HAPLOMAKE_GENERIC_PROC
//   Distinct from HM_MAKE (nextflow/modules/local/hapmake/main) — supports
//   --hapmake_format/--hapmake_reverse for lifting over arbitrary BLOCK/AGP/BED
//   structure files, not just a HaploFill .structure.block.
// --------------------------------------------------------------------------
process HAPLOMAKE_GENERIC_PROC {

    label 'process_medium'

    conda "${projectDir}/nextflow/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploMake", mode: 'copy'

    input:
    path fasta_files
    path structure_block

    output:
    path "${params.out}.fasta",                  emit: fasta
    path "${params.out}.structure.agp",          emit: agp,                optional: true
    path "${params.out}.legacy_structure.agp",   emit: legacy_agp,         optional: true
    path "${params.out}.annotation.gff3",        emit: gff3,               optional: true
    path "${params.out}.bed",                    emit: bed,                optional: true
    path "${params.out}.dropped_loci.txt",       emit: dropped_loci,       optional: true
    path "${params.out}.multiple_copy_loci.txt", emit: multiple_copy_loci, optional: true

    script:
    def fasta_list = fasta_files instanceof List ? fasta_files.join(',') : fasta_files
    def cmd = "hapmake.py"
    cmd    += " -f ${fasta_list}"
    cmd    += " -s ${structure_block}"
    cmd    += " -o ${params.out}"
    if (params.hapmake_format)       cmd += " --format ${params.hapmake_format}"
    if (params.hapmake_prefix)       cmd += " -p ${params.hapmake_prefix}"
    if (params.hapmake_agp)          cmd += " -a ${params.hapmake_agp}"
    if (params.hapmake_gff3)         cmd += " --gff3 ${params.hapmake_gff3}"
    if (params.hapmake_bed)          cmd += " -b ${params.hapmake_bed}"
    if (params.hapmake_gap)          cmd += " --gap ${params.hapmake_gap}"
    if (params.hapmake_noagp)        cmd += " --noagp"
    if (params.hapmake_reverse)      cmd += " --reverse"

    """
    ${cmd}
    """
}

// --------------------------------------------------------------------------
// Named workflow: HAPLOMAKE_GENERIC (--step haplomake_generic, fully standalone, any structure file)
// --------------------------------------------------------------------------
workflow HAPLOMAKE_GENERIC {

    if (!params.fasta) {
        log.error "[ERROR] --fasta is required"
        exit 1
    }
    if (!params.structure_block) {
        log.error "[ERROR] --structure_block is required"
        exit 1
    }

    def block_file = file(params.structure_block)
    if (!block_file.exists()) {
        log.error "[ERROR] Structure file not found: ${block_file}"
        exit 1
    }

    def fasta_ch = Channel.fromPath(params.fasta.tokenize(',')).collect()

    HAPLOMAKE_GENERIC_PROC(fasta_ch, Channel.value(block_file))
}
