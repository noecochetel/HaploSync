#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HM_MAKE } from './modules/local/hapmake/main'

// --------------------------------------------------------------------------
// Help message
// --------------------------------------------------------------------------
def helpHaploMake() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                   HaploSync — HaploMake v2.0                        ║
    ║          Assembly construction from a structure description         ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run nextflow/haplomake.nf [options]
        nextflow run nextflow/haplomake.nf -params-file params.yml

    Constructs new pseudomolecule FASTA and AGP from a structure block file
    (produced by HaploFill or written manually).

    ── Required ────────────────────────────────────────────────────────────
        --hap1_fasta        Hap1 FASTA
        --hap2_fasta        Hap2 FASTA
        --structure_block   Structure block file (.structure.block)

    ── Optional inputs ─────────────────────────────────────────────────────
        --unplaced_fasta    Unplaced sequences FASTA

    ── HaploMake options ────────────────────────────────────────────────────
        --hapmake_prefix        Sequence ID prefix
        --hapmake_agp           AGP to lift over into new assembly space
        --hapmake_gff3          GFF3 annotation to translate
        --hapmake_bed           BED file to translate
        --hapmake_gap           Gap size in bp             [default: 1000]
        --hapmake_skipoverlap   Skip overlap trimming      [default: false]
        --hapmake_noagp         Skip AGP output            [default: false]
        --hapmake_unplaced      Override unplaced sequences FASTA

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

    ── Examples ────────────────────────────────────────────────────────────
        # From a HaploFill structure block
        nextflow run nextflow/haplomake.nf -profile mamba \\
            --hap1_fasta hap1.fasta --hap2_fasta hap2.fasta \\
            --structure_block myproject.structure.block \\
            --out myproject --outdir results

        # With unplaced sequences and annotation translation
        nextflow run nextflow/haplomake.nf -profile mamba \\
            --hap1_fasta hap1.fasta --hap2_fasta hap2.fasta \\
            --unplaced_fasta unplaced.fasta \\
            --structure_block myproject.structure.block \\
            --hapmake_agp previous.agp \\
            --hapmake_gff3 annotation.gff3 \\
            --out myproject --outdir results

        # Using a params file
        nextflow run nextflow/haplomake.nf -profile mamba \\
            -params-file nextflow/params_haplomake.yml
    """.stripIndent()
}

// --------------------------------------------------------------------------
// Entry point: HAPLOMAKE
// --------------------------------------------------------------------------
workflow {

    if (params.help) {
        helpHaploMake()
        exit 0
    }

    def required = [
        hap1_fasta:      '--hap1_fasta',
        hap2_fasta:      '--hap2_fasta',
        structure_block: '--structure_block'
    ]
    required.each { param, flag ->
        if (!params[param]) {
            log.error "[ERROR] ${flag} is required"
            helpHaploMake()
            exit 1
        }
    }

    def block_file = file(params.structure_block)
    if (!block_file.exists()) {
        log.error "[ERROR] Structure block not found: ${block_file}"
        exit 1
    }

    def hap1_fasta = Channel.fromPath(params.hap1_fasta)
    def hap2_fasta = Channel.fromPath(params.hap2_fasta)
    def un_fasta   = params.unplaced_fasta
                         ? Channel.fromPath(params.unplaced_fasta)
                         : Channel.value([])

    HM_MAKE(hap1_fasta, hap2_fasta, un_fasta, Channel.value(block_file))
}

workflow.onComplete {
    log.info (workflow.success
        ? "\n[HaploSync] HaploMake completed successfully.\n  Results: ${params.outdir}"
        : "\n[HaploSync] HaploMake failed. Check logs for details.")
}
