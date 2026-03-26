#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PSEUDOMOLECULE_GENERATION } from './workflows/pseudomolecule_generation'

// --------------------------------------------------------------------------
// Help
// --------------------------------------------------------------------------
def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                        HaploSync v2.0                               ║
    ║        Haplotype-resolved pseudomolecule assembly pipeline          ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf -profile <profile> [options]

    ── Required ────────────────────────────────────────────────────────────
        --input_fasta       Draft assembly FASTA
        --guide_genome      Guide/reference genome FASTA
                            (required if --markers is not set)
        --markers           Markers FASTA
                            (required if --guide_genome is not set)

    ── Optional inputs ─────────────────────────────────────────────────────
        --markers_map       Markers genetic map
        --input_agp         Input AGP structure file
        --input_groups      Sequence grouping file
        --legacy_groups     Legacy grouping file
        --gff3              Gene annotation GFF3
        --reference         Reference genome (used by HaploDup)

    ── Output ──────────────────────────────────────────────────────────────
        --out               Output files prefix          [default: out]
        --prefix            Sequence ID prefix           [default: NEW]
        --outdir            Results directory             [default: results]

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process         [default: 4]

    ── HaploSplit ──────────────────────────────────────────────────────────
        --haplosplit_opts   Extra HaploSplit flags as a quoted string
                            e.g. '--skip_chimeric_qc --reuse_intermediate'

    ── HaploDup ────────────────────────────────────────────────────────────
        --run_haplodup      Run HaploDup after HaploSplit [default: false]
        --haplodup_opts     Extra HaploDup flags as a quoted string
                            e.g. '--only_paired_dotplots --skip_chr_pair_reports'

    ── Profiles ────────────────────────────────────────────────────────────
        -profile standard   Local execution
        -profile conda      Local execution with conda
        -profile mamba      Local execution with mamba/micromamba
        -profile hpc        SLURM execution with mamba/micromamba

    ── Examples ────────────────────────────────────────────────────────────
        # Basic run with guide genome
        nextflow run main.nf -profile mamba \\
            --input_fasta assembly.fasta \\
            --guide_genome reference.fasta \\
            --out myproject --outdir results

        # With HaploDup and markers
        nextflow run main.nf -profile mamba \\
            --input_fasta assembly.fasta \\
            --guide_genome reference.fasta \\
            --markers markers.fasta \\
            --run_haplodup \\
            --out myproject --outdir results

        # Resume after manual curation
        nextflow run main.nf -profile mamba -resume \\
            --input_fasta assembly.fasta \\
            --guide_genome reference.fasta \\
            --out myproject --outdir results
    """.stripIndent()
}

// --------------------------------------------------------------------------
// Entry point
// --------------------------------------------------------------------------
workflow {

    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate required inputs
    if (!params.input_fasta) {
        log.error "[ERROR] --input_fasta is required"
        helpMessage()
        exit 1
    }
    if (!params.guide_genome && !params.markers) {
        log.error "[ERROR] At least one of --guide_genome or --markers must be provided"
        helpMessage()
        exit 1
    }

    PSEUDOMOLECULE_GENERATION()
}
