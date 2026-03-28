/*
 * Process: HAPLODUP_ALIGN
 *
 * Runs all pairwise nucmer alignments for HaploDup:
 *   Hap1×Hap1, Hap2×Hap2, Hap1×Hap2, Hap2×Hap1
 *   (and optionally vs Reference if params.reference is set).
 * Wraps scripts/haplodup_align.py.
 *
 * Outputs:
 *   {out}.HaploDup_dir/*.delta   — nucmer delta files for HAPLODUP_REPORT
 */

nextflow.enable.dsl = 2

process ALIGN {

    label 'process_high'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploDup", mode: 'copy'

    input:
    path hap1_fasta
    path hap2_fasta
    path un_fasta
    path correspondence

    output:
    path "${params.out}.HaploDup_dir/*.delta", emit: delta_files

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def fasta_list = "${hap1_fasta},${hap2_fasta},${un_fasta}"
    def cmd        = "python3 ${haplosync}/scripts/haplodup_align.py"
    cmd           += " -f ${fasta_list}"
    cmd           += " -c ${correspondence}"
    cmd           += " -o ${params.out}"
    cmd           += " -t ${task.cpus}"
    if (params.reference) cmd += " -r ${params.reference}"

    """
    ${cmd}
    """
}
