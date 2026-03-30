/*
 * Process: REJECTED_QC
 *
 * Generates per-unplaced-sequence QC reports from HaploSplit outputs.
 * Wraps scripts/rejected_qc.py.
 *
 * Outputs:
 *   {out}.structure_comparison/   — per-sequence HTML reports + index
 */

nextflow.enable.dsl = 2

process REJECTED {

    label 'process_medium'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploSplit", mode: 'copy'

    input:
    path unused_list
    path correspondence
    path fasta
    path agp
    path markers_bed
    path legacy_agp

    output:
    path "${params.out}.structure_comparison/", emit: structure_comparison

    script:
    def haplosync  = params.haplosync_dir ?: "${projectDir}/.."
    def agp_files  = agp instanceof List ? agp.join(' ') : agp
    def fasta_list = fasta instanceof List ? fasta.join(',') : fasta
    def cmd        = "cat ${agp_files} > combined_rejected.agp && python3 ${haplosync}/scripts/rejected_qc.py"
    cmd           += " -u ${unused_list}"
    cmd           += " -c ${correspondence}"
    cmd           += " -f ${fasta_list}"
    cmd           += " -a combined_rejected.agp"
    cmd           += " -o ${params.out}"
    cmd           += " -t ${task.cpus}"
    if (markers_bed)              cmd += " -b ${markers_bed}"
    if (params.markers_map)       cmd += " -m ${params.markers_map}"
    if (legacy_agp)               cmd += " --legacy_agp ${legacy_agp}"
    if (params.input_groups)      cmd += " --input_groups ${params.input_groups}"
    if (params.legacy_groups)     cmd += " --legacy_groups ${params.legacy_groups}"

    """
    ${cmd}
    """
}
