/*
 * Process: HAPLOSPLIT
 *
 * Wraps HaploSplit.py as a single Nextflow process.
 *
 * NOTE (Phase 1): HaploSplit.py is a stateful monolith — its internal
 * nucmer alignment phases are coupled to in-memory data structures and
 * cannot yet be split into parallel sub-processes without deeper Python
 * refactoring. Phase 2 will decompose ALIGN_CHR and ALIGN_UNPLACED as
 * parallel processes feeding a leaner RECONSTRUCT process.
 *
 * Outputs:
 *   hap1_fasta      — Haplotype 1 pseudomolecules
 *   hap2_fasta      — Haplotype 2 pseudomolecules
 *   un_fasta        — Unplaced sequences
 *   agp             — All AGP files (hap1, hap2, Un)
 *   correspondence  — Chr/Hap1/Hap2 correspondence (always written)
 *   markers_bed     — Marker hits BED (optional, if --markers provided)
 *   legacy_agp      — Legacy structure AGP (optional, if --input_agp provided)
 *   annotation      — Translated annotation GFF3 (optional, if --gff3 provided)
 */

nextflow.enable.dsl = 2

process HAPLOSPLIT {

    label 'process_high'

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/haplosplit", mode: 'copy'

    output:
    path "${params.out}.1.fasta",                  emit: hap1_fasta
    path "${params.out}.2.fasta",                  emit: hap2_fasta
    path "${params.out}.Un.fasta",                 emit: un_fasta
    path "${params.out}.*.agp",                    emit: agp
    path "${params.out}.correspondence.tsv",       emit: correspondence
    path "${params.out}.markers.bed",              emit: markers_bed,  optional: true
    path "${params.out}.legacy_structure.agp",     emit: legacy_agp,   optional: true
    path "${params.out}.annotation.gff3",          emit: annotation,   optional: true

    script:
    // Build the HaploSplit command from pipeline params
    def cmd = "python3 ${projectDir}/../HaploSplit.py"
    cmd    += " -i ${params.input_fasta}"
    if (params.guide_genome)  cmd += " -g ${params.guide_genome}"
    if (params.markers)       cmd += " -n ${params.markers}"
    if (params.markers_map)   cmd += " -m ${params.markers_map}"
    if (params.input_agp)     cmd += " -a ${params.input_agp}"
    if (params.input_groups)  cmd += " --input_groups ${params.input_groups}"
    if (params.legacy_groups) cmd += " --legacy_groups ${params.legacy_groups}"
    if (params.gff3)          cmd += " --GFF3 ${params.gff3}"
    if (params.reference)     cmd += " --reference ${params.reference}"
    cmd    += " -o ${params.out}"
    cmd    += " -p ${params.prefix}"
    cmd    += " -c ${task.cpus}"
    if (params.haplosplit_opts) cmd += " ${params.haplosplit_opts}"

    """
    ${cmd}
    """
}
