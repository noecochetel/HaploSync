/*
 * Process: TRANSLATE
 *
 * Translates coordinates from query sequence space to pseudomolecule space.
 * Each translation type (markers, legacy AGP, annotation) is independent and
 * can be enabled/disabled based on what inputs are available.
 * Wraps scripts/translate.py.
 */

nextflow.enable.dsl = 2

process TRANSLATE {

    label 'process_low'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploSplit", mode: 'copy'

    input:
    path hap1_agp
    path hap2_agp
    path un_agp
    // FASTA files are needed only for GFF3 translation (sequence length lookup)
    path hap1_fasta
    path hap2_fasta
    path un_fasta

    output:
    path "${params.out}.markers.bed",         emit: markers_bed,  optional: true
    path "${params.out}.legacy_structure.agp", emit: legacy_agp,  optional: true
    path "${params.out}.annotation.gff3",     emit: annotation,   optional: true
    path "${params.out}.broken_genes.txt",    emit: broken_genes, optional: true

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def cmd = "python3 ${haplosync}/scripts/translate.py"
    cmd    += " -o ${params.out}"
    if (params.markers)   cmd += " -n ${params.markers}"
    if (params.input_agp) cmd += " -a ${params.input_agp}"
    if (params.gff3)      cmd += " --GFF3 ${params.gff3}"
    if (params.No2)       cmd += " --N2"

    """
    ${cmd}
    """
}
