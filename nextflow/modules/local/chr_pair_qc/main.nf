/*
 * Process: CHR_PAIR_QC
 *
 * Generates per-chromosome Hap1 vs Hap2 overview QC reports from
 * HaploSplit outputs. Wraps scripts/chr_pair_qc.py.
 *
 * Outputs:
 *   {out}.chr_pair_reports/   — per-chromosome HTML reports + index
 */

nextflow.enable.dsl = 2

process CHR_PAIR {

    label 'process_medium'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploSplit", mode: 'copy'

    input:
    path correspondence
    path agp
    path markers_bed
    path legacy_agp

    output:
    path "${params.out}.chr_pair_reports/", emit: chr_pair_reports

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def agp_files = agp instanceof List ? agp.join(' ') : agp
    def cmd       = "cat ${agp_files} > combined_chr_pair.agp && python3 ${haplosync}/scripts/chr_pair_qc.py"
    cmd          += " -c ${correspondence}"
    cmd          += " -a combined_chr_pair.agp"
    cmd          += " -o ${params.out}"
    cmd          += " -t ${task.cpus}"
    if (markers_bed)              cmd += " -b ${markers_bed}"
    if (params.markers_map)       cmd += " -m ${params.markers_map}"
    if (legacy_agp)               cmd += " --legacy_agp ${legacy_agp}"
    if (params.input_groups)      cmd += " --input_groups ${params.input_groups}"
    if (params.legacy_groups)     cmd += " --legacy_groups ${params.legacy_groups}"

    """
    ${cmd}
    """
}
