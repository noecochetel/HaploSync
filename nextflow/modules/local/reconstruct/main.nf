/*
 * Process: RECONSTRUCT
 *
 * Builds AGP files, FASTA sequences, correspondence TSV, and relationship
 * reports from the tiling path lists produced by BUILD_PATHS.
 * Wraps scripts/reconstruct.py.
 */

nextflow.enable.dsl = 2

process RECONSTRUCT {

    label 'process_medium'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploSplit", mode: 'copy'

    input:
    path hap1_list
    path hap2_list
    path un_list
    path unused_list

    output:
    path "${params.out}.1.fasta",                          emit: hap1_fasta
    path "${params.out}.2.fasta",                          emit: hap2_fasta,        optional: true
    path "${params.out}.Un.fasta",                         emit: un_fasta
    path "${params.out}.1.agp",                            emit: hap1_agp
    path "${params.out}.2.agp",                            emit: hap2_agp,          optional: true
    path "${params.out}.Un.agp",                           emit: un_agp
    path "${params.out}.correspondence.tsv",               emit: correspondence
    path "${params.out}.conflicting_pseudomolecules.txt",  emit: conflicting_pm
    path "${params.out}.unplaced_to_pseudomolecule.txt",   emit: unplaced_report

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def cmd = "python3 ${haplosync}/scripts/reconstruct.py"
    cmd    += " -i ${params.input_fasta}"
    cmd    += " -o ${params.out}"
    cmd    += " -p ${params.prefix}"
    if (params.gap)        cmd += " --gapsize ${params.gap}"
    if (params.conc)       cmd += " --concatenate ${params.conc}"
    if (params.known)      cmd += " -k ${params.known}"
    if (params.alternative_groups) cmd += " --alternative_groups ${params.alternative_groups}"
    if (params.No2)        cmd += " --N2"

    """
    ${cmd}
    """
}
