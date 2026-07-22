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

    conda "${projectDir}/nextflow/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploDup", mode: 'copy'

    input:
    path hap1_fasta, stageAs: 'hap1.fasta'
    path hap2_fasta, stageAs: 'hap2.fasta'
    path un_fasta
    path correspondence
    path gff

    output:
    path "versions.yml", emit: versions
    path "${params.out}.HaploDup_dir/CDS.on.genome.gmap.gff3", emit: gmap_gff3

    script:
    def fasta_list = "${hap1_fasta},${hap2_fasta},${un_fasta}"
    def cmd        = "haplodup_gmap.py"
    cmd           += " -f ${fasta_list}"
    cmd           += " -c ${correspondence}"
    cmd           += " -g ${gff}"
    cmd           += " -o ${params.out}"
    cmd           += " -t ${task.cpus}"
    if (params.haplodup_feature) cmd += " --feature ${params.haplodup_feature}"

    """
    ${cmd}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //')
        gmap: \$(gmap --version 2>&1 | head -n1 | sed "s/.*version //;s/ .*//")
    END_VERSIONS
    """
}
