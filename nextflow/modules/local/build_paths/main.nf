/*
 * Process: BUILD_PATHS
 *
 * Marker QC, chimera detection, graph construction, and DAG tiling path
 * selection for both haplotypes. Wraps scripts/build_paths.py.
 *
 * Outputs four list files consumed by RECONSTRUCT:
 *   {out}.1.list                — Hap1 tiling paths
 *   {out}.2.list                — Hap2 tiling paths
 *   {out}.Un.list               — Unplaced sequence IDs
 *   {out}.unused_sequences.list — Chr-assigned but unplaced sequences
 */

nextflow.enable.dsl = 2

process BUILD_PATHS {

    label 'process_high'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploSplit", mode: 'copy'

    output:
    path "${params.out}.1.list",                   emit: hap1_list
    path "${params.out}.2.list",                   emit: hap2_list,    optional: true
    path "${params.out}.Un.list",                  emit: un_list
    path "${params.out}.unused_sequences.list",    emit: unused_list
    path "${params.out}.missing_orientation.hap1.txt", emit: missing_orient_hap1, optional: true
    path "${params.out}.missing_orientation.hap2.txt", emit: missing_orient_hap2, optional: true
    path "${params.out}.unknown_markers.txt",          emit: unknown_markers,     optional: true

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def cmd = "python3 ${haplosync}/scripts/build_paths.py"
    cmd    += " -i ${params.input_fasta}"
    if (params.markers)      cmd += " -n ${params.markers}"
    if (params.markers_map)  cmd += " -m ${params.markers_map}"
    cmd    += " -o ${params.out}"
    cmd    += " -c ${task.cpus}"
    if (params.guide_genome)       cmd += " -g ${params.guide_genome}"
    if (params.local_alignment)    cmd += " -l ${params.local_alignment}"
    if (params.run_alignment)      cmd += " --align"
    if (params.mapping_tool)       cmd += " -t '${params.mapping_tool}'"
    if (params.hitgap)             cmd += " --hitgap ${params.hitgap}"
    if (params.distance1)          cmd += " --distance1 ${params.distance1}"
    if (params.distance2)          cmd += " --distance2 ${params.distance2}"
    if (params.reuse_intermediate) cmd += " --reuse_intermediate"
    if (params.input_agp)     cmd += " -a ${params.input_agp}"
    if (params.gff3)          cmd += " --GFF3 ${params.gff3}"
    if (params.exclusion)     cmd += " -e ${params.exclusion}"
    if (params.known)         cmd += " -k ${params.known}"
    if (params.alternative_groups) cmd += " --alternative_groups ${params.alternative_groups}"
    if (params.Require1)      cmd += " --R1 ${params.Require1}"
    if (params.Require2)      cmd += " --R2 ${params.Require2}"
    if (params.force_direction1) cmd += " --F1"
    if (params.force_direction2) cmd += " --F2"
    if (params.minR1)         cmd += " --min1 ${params.minR1}"
    if (params.minR2)         cmd += " --min2 ${params.minR2}"
    if (params.Blacklist1)    cmd += " --B1 ${params.Blacklist1}"
    if (params.Blacklist2)    cmd += " --B2 ${params.Blacklist2}"
    if (params.input_groups)  cmd += " --input_groups ${params.input_groups}"
    if (params.No2)           cmd += " --N2"
    if (params.filter_hits)   cmd += " -f"
    if (params.extended_region) cmd += " --extended"
    if (params.skip_chimeric_qc) cmd += " --skip_chimeric_qc"
    if (params.disable_marker_ploidy_check) cmd += " --disable_marker_ploidy_check"
    if (params.conflict_resolution) cmd += " --conflict ${params.conflict_resolution}"
    if (params.allow_rearrangements) cmd += " --allow_rearrangements"
    if (params.path1)              cmd += " -1 ${params.path1}"
    if (params.path2)              cmd += " -2 ${params.path2}"
    if (params.required_as_path)   cmd += " --required_as_path"
    if (params.only_markers)       cmd += " --only_markers"
    if (params.avoid_rejected_qc)  cmd += " --avoid_rejected_qc"

    """
    ${cmd}
    """
}
