/*
 * Process: HF_PAIR
 *
 * HaploFill Steps 4-5: haplotype pairing + gap region preparation.
 *
 * Step 4 — pairwise hap-to-hap minimap2 alignment, PAF uniquification,
 *           pairing JSON files (one per chromosome pair)
 * Step 5 — unmatched region extraction, gap mate position mapping,
 *           flanking sequence extraction (upgradable_regions.json.gz)
 *
 * Runs after HF_PLOIDY. All chromosome pairs currently run as a single
 * job (HaploFill iterates internally). Future improvement: scatter per pair.
 *
 * Inputs:
 *   temp_dir — HaploFill temp folder with ploidy categories (from HF_PLOIDY)
 *
 * Outputs:
 *   {out}_tmp/ — updated temp folder with gap descriptors per chromosome
 */

nextflow.enable.dsl = 2

process HF_PAIR {

    label 'process_high'

    conda "${projectDir}/nextflow/envs/haplosync.yml"

    input:
    path temp_dir

    output:
    path "versions.yml", emit: versions
    path "${params.out}_tmp/", emit: temp_dir

    script:
    def cmd = "hapfill_pair.py"
    cmd    += " -1 ${params.hapfill_hap1}"
    cmd    += " -2 ${params.hapfill_hap2}"
    if (params.hapfill_unplaced)       cmd += " -U ${params.hapfill_unplaced}"
    cmd    += " -c ${params.hapfill_correspondence}"
    cmd    += " -r ${params.hapfill_repeats}"
    cmd    += " -o ${params.out}"
    cmd    += " -t ${params.out}_tmp"
    cmd    += " -C"
    if (params.hapfill_b1)             cmd += " --b1 ${params.hapfill_b1}"
    if (params.hapfill_b2)             cmd += " --b2 ${params.hapfill_b2}"
    if (params.hapfill_flanking)       cmd += " --flanking ${params.hapfill_flanking}"
    if (params.hapfill_map_threads)    cmd += " --map_threads ${params.hapfill_map_threads}"

    """
    ${cmd}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python3: \$(python3 --version | sed 's/Python //')
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}
