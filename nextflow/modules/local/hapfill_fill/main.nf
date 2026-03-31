/*
 * Process: HF_FILL
 *
 * HaploFill Step 6: gap patching.
 *
 * Step 6.1 — gather gap information from descriptor files
 * Step 6.2 — align unplaced sequences on gap flanking targets (minimap2)
 * Step 6.3 — select best filler per gap using coverage + alignment scoring
 *
 * Runs after HF_PAIR. Fast (~11 min in typical runs).
 *
 * Inputs:
 *   temp_dir — HaploFill temp folder with gap descriptors (from HF_PAIR)
 *
 * Outputs:
 *   {out}.structure.block          — gap fill instructions for HM_MAKE
 *   {out}.gap_filling_findings.txt — detailed per-gap report
 */

nextflow.enable.dsl = 2

process HF_FILL {

    label 'process_medium'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploFill", mode: 'copy'

    input:
    path temp_dir
    path hap1_fasta
    path hap2_fasta
    path un_fasta
    path correspondence
    path repeats

    output:
    path "${params.out}.structure.block",          emit: structure_block
    path "${params.out}.gap_filling_findings.txt", emit: findings

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def cmd = "python3 ${haplosync}/scripts/hapfill_fill.py"
    cmd    += " -1 ${hap1_fasta}"
    cmd    += " -2 ${hap2_fasta}"
    if (un_fasta)       cmd += " -U ${un_fasta}"
    cmd    += " -c ${correspondence}"
    cmd    += " -r ${repeats}"
    cmd    += " -o ${params.out}"
    cmd    += " -t ${params.out}_tmp"
    cmd    += " -C"
    if (params.hapfill_b1)          cmd += " --b1 ${params.hapfill_b1}"
    if (params.hapfill_b2)          cmd += " --b2 ${params.hapfill_b2}"
    if (params.hapfill_coverage)    cmd += " --coverage ${params.hapfill_coverage}"
    if (params.hapfill_nohomozygous) cmd += " --nohomozygous"

    """
    ${cmd}
    """
}
