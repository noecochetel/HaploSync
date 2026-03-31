/*
 * Process: HF_SETUP
 *
 * HaploFill Step 1: sequence splitting, pairing setup, repeat parsing,
 * gap detection. Produces the temp folder structure used by all subsequent
 * HaploFill modules.
 *
 * Outputs:
 *   {out}_tmp/   — HaploFill temp folder with per-sequence FASTA, gap BED,
 *                  repeat BED, and status/config JSON files
 */

nextflow.enable.dsl = 2

process HF_SETUP {

    label 'process_medium'

    conda "${projectDir}/envs/haplosync.yml"

    input:
    path hap1_fasta
    path hap2_fasta
    path un_fasta
    path correspondence
    path repeats

    output:
    path "${params.out}_tmp/", emit: temp_dir

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def cmd = "python3 ${haplosync}/scripts/hapfill_setup.py"
    cmd    += " -1 ${hap1_fasta}"
    cmd    += " -2 ${hap2_fasta}"
    if (un_fasta)      cmd += " -U ${un_fasta}"
    cmd    += " -c ${correspondence}"
    cmd    += " -r ${repeats}"
    cmd    += " -o ${params.out}"
    cmd    += " -t ${params.out}_tmp"
    if (params.hapfill_repeats_format) cmd += " --repeats_format ${params.hapfill_repeats_format}"
    if (params.hapfill_exclusion)      cmd += " --exclusion ${params.hapfill_exclusion}"
    if (params.hapfill_known)          cmd += " --known ${params.hapfill_known}"
    if (params.hapfill_overwrite)      cmd += " --overwrite"

    """
    ${cmd}
    """
}
