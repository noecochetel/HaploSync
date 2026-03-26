/*
 * Process: HAPLODUP
 *
 * Wraps HaploDup.py as a single Nextflow process.
 * Receives its inputs from HAPLOSPLIT via channel outputs.
 */

nextflow.enable.dsl = 2

process HAPLODUP {

    label 'process_high'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/haplodup", mode: 'copy'

    input:
    path hap1_fasta
    path hap2_fasta
    path un_fasta
    path agp
    path correspondence
    path markers_bed
    path legacy_agp
    path annotation

    output:
    path "${params.out}.HaploDup_dir/", emit: haplodup_dir
    path "${params.out}.html",          emit: report, optional: true

    script:
    def haplosync  = params.haplosync_dir ?: "${projectDir}/.."
    def fasta_arg  = "${hap1_fasta},${hap2_fasta},${un_fasta}"
    def agp_files  = agp instanceof List ? agp.join(' ') : agp
    def cmd        = "cat ${agp_files} > combined.agp && python3 ${haplosync}/HaploDup.py"
    cmd           += " -f ${fasta_arg}"
    cmd           += " -c ${correspondence}"
    cmd           += " --agp combined.agp"
    cmd           += " -o ${params.out}"
    cmd           += " -t ${task.cpus}"
    if (markers_bed)  cmd += " -b ${markers_bed}"
    if (params.markers_map) cmd += " --markers_map ${params.markers_map}"
    if (legacy_agp)   cmd += " --legacy_agp ${legacy_agp}"
    if (params.input_groups)  cmd += " --input_groups ${params.input_groups}"
    if (params.legacy_groups) cmd += " --legacy_groups ${params.legacy_groups}"
    if (annotation)   cmd += " -g ${annotation}"
    if (params.reference)     cmd += " -r ${params.reference}"
    if (params.haplodup_opts) cmd += " ${params.haplodup_opts}"

    """
    ${cmd}
    """
}
