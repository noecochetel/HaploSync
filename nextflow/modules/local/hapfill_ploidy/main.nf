/*
 * Process: HF_PLOIDY
 *
 * HaploFill Step 3: local ploidy classification.
 *
 * Gather step — waits for all HF_COVERAGE jobs to complete, then:
 *   3.1 Calculates genome-wide median coverage (all chromosomes needed)
 *   3.2 Categorises each position as haploid / diploid / repetitive
 *
 * The coverage signal files produced by HF_COVERAGE are staged into the
 * per-chromosome subdirectories of temp_dir before running Step 3.
 *
 * Inputs:
 *   temp_dir      — HaploFill temp folder (from HF_SETUP)
 *   signal_files  — all .cov.txt.gz files (collected from HF_COVERAGE)
 *   range_beds    — all .cov.bed.gz files (collected from HF_COVERAGE)
 *
 * Outputs:
 *   {out}_tmp/   — updated temp folder with category BED files per chromosome
 */

nextflow.enable.dsl = 2

process HF_PLOIDY {

    label 'process_medium'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploFill", mode: 'copy'

    input:
    path temp_dir
    path signal_files
    path range_beds

    output:
    path "${params.out}_tmp/", emit: temp_dir

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."

    // Stage coverage files into the correct per-chromosome subdirectories
    // inside temp_dir before calling HaploFill Step 3
    def cmd = """
    # Place each signal/bed file into its chromosome subdirectory
    for f in *.cov.txt.gz *.cov.bed.gz ; do
        chr=\$(basename \$f | sed 's/\\.cov\\..*//')
        target="${params.out}_tmp/\${chr}"
        mkdir -p "\${target}"
        cp "\$f" "\${target}/"
    done

    # Update conf.files.json to add coverage_file paths.
    # HF_SETUP runs step 1 only, so write_coverage_bed() (step 2) is never
    # called and the coverage_file key is absent from conf.files.json.
    # HaploFill --resume 3 checks for this key and exits 25 if it is missing.
    python3 - <<'PYEOF'
import json, os
conf = "${params.out}_tmp/conf.files.json"
with open(conf) as fh:
    db = json.load(fh)
for seq, info in db["sequences"].items():
    folder = info.get("folder", "")
    cov = os.path.join(folder, seq + ".cov.txt.gz")
    if os.path.exists(cov):
        db["sequences"][seq]["coverage_file"] = cov
with open(conf, "w") as fh:
    json.dump(db, fh, indent=4)
PYEOF

    python3 ${haplosync}/scripts/hapfill_ploidy.py \\
        -1 ${params.hapfill_hap1} \\
        -2 ${params.hapfill_hap2} \\
        -c ${params.hapfill_correspondence} \\
        -r ${params.hapfill_repeats} \\
        -o ${params.out} \\
        -t ${params.out}_tmp
    """

    cmd
}
