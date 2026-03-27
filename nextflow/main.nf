#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPLOSPLIT } from './workflows/pseudomolecule_generation'

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
        nextflow run main.nf -profile <profile> -params-file params.yml

    ── Required ────────────────────────────────────────────────────────────
        --input_fasta       Draft assembly FASTA
        --markers           Marker hits on input sequences (BED, 4 columns)
        --markers_map       Marker genetic map (chr, position, marker_id)

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
        --input_groups      Sequence grouping file (--input_groups)

    ── Assembly behaviour ───────────────────────────────────────────────────
        --gap               Gap size in bp               [default: 1000]
        --skip_chimeric_qc  Skip chimeric sequence QC    [default: false]

    ── Output ──────────────────────────────────────────────────────────────
        --out               Output files prefix          [default: out]
        --prefix            Sequence ID prefix           [default: NEW]
        --outdir            Results directory            [default: results]

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── Extra flags (per module) ─────────────────────────────────────────────
        --build_paths_opts  Extra build_paths.py flags as a quoted string
        --reconstruct_opts  Extra reconstruct.py flags as a quoted string
        --translate_opts    Extra translate.py flags as a quoted string
        --haplodup_opts     Extra HaploDup flags as a quoted string

    ── HaploDup ────────────────────────────────────────────────────────────
        --run_haplodup      Run HaploDup after reconstruction [default: false]
        --reference         Reference genome for HaploDup dotplots

    ── Profiles ────────────────────────────────────────────────────────────
        -profile standard   Local execution
        -profile conda      Local execution with conda
        -profile mamba      Local execution with mamba/micromamba
        -profile hpc        SLURM execution with mamba/micromamba

    ── Examples ────────────────────────────────────────────────────────────
        # Diploid assembly with genetic map (no reference)
        nextflow run main.nf -profile mamba \\
            --input_fasta assembly.fasta \\
            --markers markers.bed \\
            --markers_map genetic_map.tsv \\
            --input_agp assembly.agp \\
            --out myproject --outdir results

        # Using a params file
        nextflow run main.nf -profile mamba -params-file params.yml

        # Resume after manual curation
        nextflow run main.nf -profile mamba -resume -params-file params.yml
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
    if (!params.markers) {
        log.error "[ERROR] --markers (marker hits BED) is required"
        helpMessage()
        exit 1
    }
    if (!params.markers_map) {
        log.error "[ERROR] --markers_map (genetic map) is required"
        helpMessage()
        exit 1
    }

    HAPLOSPLIT()
}

workflow.onComplete {
    log.info (workflow.success
        ? "\n[HaploSync] Pipeline completed successfully.\n  Results: ${params.outdir}/HaploSplit"
        : "\n[HaploSync] Pipeline failed. Check logs for details.")
}
