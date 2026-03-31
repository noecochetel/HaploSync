#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPLOSYNC_GAP_FILL } from './workflows/gap_fill'
include { ALIGN  as HD_ALIGN } from './modules/local/haplodup_align/main'
include { GMAP   as HD_GMAP  } from './modules/local/haplodup_gmap/main'
include { REPORT as HD_REPORT} from './modules/local/haplodup_report/main'

// --------------------------------------------------------------------------
// Help messages
// --------------------------------------------------------------------------
def helpGapFill() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                    HaploSync — HaploFill v2.0                       ║
    ║              Gap-filling and pseudomolecule patching                ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run gap_fill.nf [options]
        nextflow run gap_fill.nf -params-file params.yml

    Pipeline: HaploFill → HaploMake → HaploDup (optional)

    ── Required ────────────────────────────────────────────────────────────
        --hapfill_hap1          Hap1 FASTA (pseudomolecules to patch)
        --hapfill_hap2          Hap2 FASTA (pseudomolecules to patch)
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
                                mosdepth is 20-50x faster, opt-in for validation

    ── HaploMake options ────────────────────────────────────────────────────
        --hapmake_prefix        Sequence ID prefix
        --hapmake_agp           AGP reference for legacy coordinate mapping
        --hapmake_gff3          Gene annotation GFF3 to translate
        --hapmake_bed           BED file to translate
        --hapmake_gap           Gap size in bp             [default: 1000]
        --hapmake_skipoverlap   Skip overlap trimming      [default: false]
        --hapmake_noagp         Skip AGP output            [default: false]
        --hapmake_unplaced      Override unplaced sequences FASTA

    ── HaploMake (optional) ─────────────────────────────────────────────────
        --run_haplomake         Build new FASTA/AGP from gap-fill result [default: false]
                                Produces: {out}.fasta, {out}.structure.agp

    ── HaploDup (optional) ──────────────────────────────────────────────────
        --run_haplodup          Run HaploDup on the gap-filled assembly  [default: false]
                                Implies --run_haplomake (HaploMake runs automatically)

        For full HaploDup options: nextflow run gap_fill.nf -entry HAPLODUP --help

    ── Output ──────────────────────────────────────────────────────────────
        --out               Output files prefix          [default: out]
        --outdir            Results directory            [default: results]

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── Profiles ────────────────────────────────────────────────────────────
        -profile standard   Local execution
        -profile conda      Local execution with conda
        -profile mamba      Local execution with mamba/micromamba
        -profile hpc        SLURM execution with mamba/micromamba

    ── RAM note ─────────────────────────────────────────────────────────────
        HF_COVERAGE runs one job per chromosome in parallel.
        Cap concurrency in nextflow.config to avoid RAM exhaustion:
          process { withName: 'HF_COVERAGE' { maxForks = 8 } }   // ~48 GB RAM
        With mosdepth (--coverage_tool mosdepth), maxForks can be much higher.

    ── Examples ────────────────────────────────────────────────────────────
        # Basic gap-fill run
        nextflow run gap_fill.nf -profile mamba \\
            --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \\
            --hapfill_correspondence correspondence.tsv \\
            --hapfill_repeats repeats.bed \\
            --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \\
            --out myproject --outdir results

        # With mosdepth coverage (faster)
        nextflow run gap_fill.nf -profile mamba \\
            --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \\
            --hapfill_correspondence correspondence.tsv \\
            --hapfill_repeats repeats.bed \\
            --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \\
            --coverage_tool mosdepth \\
            --out myproject --outdir results

        # Gap fill only — inspect .structure.block before building assembly
        nextflow run gap_fill.nf -profile mamba \\
            --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \\
            --hapfill_correspondence correspondence.tsv \\
            --hapfill_repeats repeats.bed \\
            --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \\
            --out myproject --outdir results

        # Gap fill + build new assembly
        nextflow run gap_fill.nf -profile mamba \\
            --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \\
            --hapfill_correspondence correspondence.tsv \\
            --hapfill_repeats repeats.bed \\
            --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \\
            --run_haplomake \\
            --out myproject --outdir results

        # Gap fill + build + HaploDup QC (--run_haplodup implies --run_haplomake)
        nextflow run gap_fill.nf -profile mamba \\
            --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \\
            --hapfill_unplaced unplaced.fasta \\
            --hapfill_correspondence correspondence.tsv \\
            --hapfill_repeats repeats.bed \\
            --hapfill_b1 hap1.bam --hapfill_b2 hap2.bam \\
            --run_haplodup \\
            --out myproject --outdir results

        # Using a params file
        nextflow run gap_fill.nf -profile mamba -params-file params.yml
    """.stripIndent()
}

def helpHaploDup() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║              HaploSync — HaploDup on gap-filled assembly            ║
    ║              Haplotype duplication QC and reporting                 ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run gap_fill.nf -entry HAPLODUP [options]
        nextflow run gap_fill.nf -entry HAPLODUP -params-file params.yml

    HaploDup reads HaploMake outputs from --outdir/HaploMake/ using the
    --out prefix. Run the default gap_fill.nf pipeline first, or point to
    an existing HaploMake output directory.

    ── Required ────────────────────────────────────────────────────────────
        --hapfill_hap1          Original Hap1 FASTA (for correspondence)
        --hapfill_hap2          Original Hap2 FASTA (for correspondence)
        --hapfill_correspondence  Chromosome correspondence TSV
        --out               Output prefix (must match gap-fill run)
                            [default: out]
        --outdir            Results directory (must match gap-fill run)
                            [default: results]

    ── Optional inputs ─────────────────────────────────────────────────────
        --hapfill_unplaced      Unplaced sequences FASTA
        --gff3                  Gene annotation GFF3

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── Examples ────────────────────────────────────────────────────────────
        nextflow run gap_fill.nf -entry HAPLODUP -profile mamba \\
            --hapfill_hap1 hap1.fasta --hapfill_hap2 hap2.fasta \\
            --hapfill_correspondence correspondence.tsv \\
            --out myproject --outdir results
    """.stripIndent()
}

// --------------------------------------------------------------------------
// Default entry point: HAPLOSYNC_GAP_FILL
//   HaploFill → HaploMake → optional HaploDup
//   Log names: HAPLOSYNC_GAP_FILL:HAPLOFILL:<PROCESS>
//              HAPLOSYNC_GAP_FILL:HAPLOMAKE:<PROCESS>
//              HAPLOSYNC_GAP_FILL:HAPLODUP:<PROCESS>
// --------------------------------------------------------------------------
workflow {

    if (params.help) {
        helpGapFill()
        exit 0
    }

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
}

// --------------------------------------------------------------------------
// Entry point: HAPLOMAKE (standalone)
//   Reads HaploFill structure block from --outdir/HaploFill/.
//   Use --structure_block to override with a custom path.
//   Processes called directly for clean log names: HAPLOMAKE:<PROCESS>
// --------------------------------------------------------------------------
workflow HAPLOMAKE {

    if (params.help) {
        log.info """
    Usage:
        nextflow run gap_fill.nf -entry HAPLOMAKE -profile mamba [options]

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
        --hapmake_skipoverlap   Skip overlap trimming
        --hapmake_noagp         Skip AGP output
        """.stripIndent()
        exit 0
    }

    def block_path = params.structure_block
        ?: "${params.outdir}/HaploFill/${params.out}.structure.block"
    def block_file = file(block_path)
    if (!block_file.exists()) {
        log.error "[ERROR] Structure block not found: ${block_file}\n         Run gap_fill.nf first or provide --structure_block"
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
// Entry point: HAPLODUP (standalone, gap-fill context)
//   Reads HaploMake outputs from --outdir/HaploMake/.
//   Processes called directly for clean log names: HAPLODUP:<PROCESS>
// --------------------------------------------------------------------------
workflow HAPLODUP {

    if (params.help) {
        helpHaploDup()
        exit 0
    }

    def hm_dir = "${params.outdir}/HaploMake"
    def pfx    = "${hm_dir}/${params.out}"

    def required_files = [ file("${pfx}.fasta") ]
    required_files.each { f ->
        if (!f.exists()) {
            log.error "[ERROR] Required HaploMake output not found: ${f}\n         Run gap_fill.nf first or check --out / --outdir"
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

workflow.onComplete {
    log.info (workflow.success
        ? "\n[HaploSync] Gap-fill pipeline completed successfully.\n  Results: ${params.outdir}"
        : "\n[HaploSync] Gap-fill pipeline failed. Check logs for details.")
}
