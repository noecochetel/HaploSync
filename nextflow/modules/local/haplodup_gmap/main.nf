/*
 * Process: HAPLODUP_GMAP
 *
 * Extracts CDS sequences from GFF3 annotation, builds GMAP indices for
 * Hap1 and Hap2, maps CDS to both haplotypes, and concatenates results.
 * Runs in parallel with HAPLODUP_ALIGN.
 * Wraps scripts/haplodup_gmap.py.
 *
 * Outputs:
 *   {out}.HaploDup_dir/CDS.on.genome.gmap.gff3   — GMAP mapping results for HAPLODUP_REPORT
 */

nextflow.enable.dsl = 2

process GMAP {

    label 'process_high'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploDup", mode: 'copy'

    input:
    path hap1_fasta
    path hap2_fasta
    path un_fasta
    path correspondence
    path gff

    output:
    path "${params.out}.HaploDup_dir/CDS.on.genome.gmap.gff3", emit: gmap_gff3

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def fasta_list = "${hap1_fasta},${hap2_fasta},${un_fasta}"
    def cmd        = "python3 ${haplosync}/scripts/haplodup_gmap.py"
    cmd           += " -f ${fasta_list}"
    cmd           += " -c ${correspondence}"
    cmd           += " -g ${gff}"
    cmd           += " -o ${params.out}"
    cmd           += " -t ${task.cpus}"
    if (params.haplodup_feature) cmd += " --feature ${params.haplodup_feature}"

    """
    ${cmd}
    """
}
