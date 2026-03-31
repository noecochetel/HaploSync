#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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
        nextflow run nextflow/haplomake.nf -params-file params_haplomake.yml

    Constructs new pseudomolecule FASTA and AGP from a structure file.
    The structure file can be a HaploFill block, an AGP, or a BED file.

    ── Required ────────────────────────────────────────────────────────────
        --fasta             Input FASTA file(s), comma-separated if multiple
        --structure_block   Structure file (.structure.block, .agp, or .bed)

    ── HaploMake options ────────────────────────────────────────────────────
        --hapmake_format        Structure file format [BLOCK|AGP|BED] [default: BLOCK]
        --hapmake_prefix        Sequence ID prefix
        --hapmake_agp           AGP to lift over into new assembly space
        --hapmake_gff3          GFF3 annotation to translate
        --hapmake_bed           BED file to translate
        --hapmake_gap           Gap size in bp             [default: 1000]
        --hapmake_skipoverlap   Skip overlap trimming      [default: false]
        --hapmake_noagp         Skip AGP output            [default: false]
        --hapmake_reverse       Reverse AGP direction (new → old) [default: false]

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
            --fasta assembly.fasta \\
            --structure_block myproject.structure.block \\
            --out myproject_new --outdir results

        # Manual curation from an edited AGP
        nextflow run nextflow/haplomake.nf -profile mamba \\
            --fasta assembly.fasta \\
            --structure_block assembly_corrected.agp \\
            --hapmake_format AGP \\
            --hapmake_prefix NEW \\
            --out assembly_corrected --outdir results

        # Using a params file
        nextflow run nextflow/haplomake.nf -profile mamba \\
            -params-file nextflow/params_haplomake.yml
    """.stripIndent()
}

// --------------------------------------------------------------------------
// Process: HAPLOMAKE
// --------------------------------------------------------------------------
process HAPLOMAKE {

    label 'process_medium'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploMake", mode: 'copy'

    input:
    path fasta_files
    path structure_block

    output:
    path "${params.out}.fasta",                      emit: fasta
    path "${params.out}.structure.agp",              emit: agp,        optional: true
    path "${params.out}.legacy_structure.agp",       emit: legacy_agp, optional: true

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def fasta_list = fasta_files instanceof List ? fasta_files.join(',') : fasta_files
    def cmd = "python3 ${haplosync}/scripts/hapmake.py"
    cmd    += " -f ${fasta_list}"
    cmd    += " -s ${structure_block}"
    cmd    += " -o ${params.out}"
    if (params.hapmake_format)       cmd += " --format ${params.hapmake_format}"
    if (params.hapmake_prefix)       cmd += " -p ${params.hapmake_prefix}"
    if (params.hapmake_agp)          cmd += " -a ${params.hapmake_agp}"
    if (params.hapmake_gff3)         cmd += " --gff3 ${params.hapmake_gff3}"
    if (params.hapmake_bed)          cmd += " -b ${params.hapmake_bed}"
    if (params.hapmake_gap)          cmd += " --gap ${params.hapmake_gap}"
    if (params.hapmake_skipoverlap)  cmd += " --skipoverlap"
    if (params.hapmake_noagp)        cmd += " --noagp"
    if (params.hapmake_reverse)      cmd += " --reverse"

    """
    ${cmd}
    """
}

// --------------------------------------------------------------------------
// Entry point: HAPLOMAKE
// --------------------------------------------------------------------------
workflow {

    if (params.help) {
        helpHaploMake()
        exit 0
    }

    if (!params.fasta) {
        log.error "[ERROR] --fasta is required"
        helpHaploMake()
        exit 1
    }
    if (!params.structure_block) {
        log.error "[ERROR] --structure_block is required"
        helpHaploMake()
        exit 1
    }

    def block_file = file(params.structure_block)
    if (!block_file.exists()) {
        log.error "[ERROR] Structure file not found: ${block_file}"
        exit 1
    }

    def fasta_ch = Channel.fromPath(params.fasta.tokenize(',')).collect()

    HAPLOMAKE(fasta_ch, Channel.value(block_file))
}

workflow.onComplete {
    log.info (workflow.success
        ? "\n[HaploSync] HaploMake completed successfully.\n  Results: ${params.outdir}"
        : "\n[HaploSync] HaploMake failed. Check logs for details.")
}
