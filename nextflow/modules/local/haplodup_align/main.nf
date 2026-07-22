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

    conda "${projectDir}/nextflow/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploDup", mode: 'copy'

    input:
    path hap1_fasta, stageAs: 'hap1.fasta'
    path hap2_fasta, stageAs: 'hap2.fasta'
    path un_fasta
    path correspondence

    output:
    path "versions.yml", emit: versions
    path "${params.out}.HaploDup_dir/*.delta", emit: delta_files

    script:
    def fasta_list = "${hap1_fasta},${hap2_fasta},${un_fasta}"
    def cmd        = "haplodup_align.py"
    cmd           += " -f ${fasta_list}"
    cmd           += " -c ${correspondence}"
    cmd           += " -o ${params.out}"
    cmd           += " -t ${task.cpus}"
    if (params.reference) cmd += " -r ${params.reference}"

    """
    ${cmd}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //')
        mummer: \$(nucmer --version 2>&1 | tail -n1 | sed "s/.*version //")
    END_VERSIONS
    """
}
