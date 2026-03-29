/*
 * Process: HF_COVERAGE
 *
 * HaploFill Step 2: per-chromosome coverage extraction.
 *
 * Scattered one job per chromosome — all chromosomes of both haplotypes
 * run in parallel. Use 'maxForks' to cap concurrency and RAM:
 *   process.withName: 'HF_COVERAGE' { maxForks = 8 }
 *
 * Coverage tools:
 *   bedtools (default) — validated, widely available
 *   mosdepth           — 20-50x faster, set --coverage_tool mosdepth to use
 *
 * Inputs:
 *   bam          — sorted, indexed BAM (full genome)
 *   temp_dir     — HaploFill temp folder from HF_SETUP
 *   chr_name     — chromosome/sequence name (val, one per job)
 *   chr_length   — chromosome length in bp (val, one per job)
 *   chr_fasta    — per-chromosome FASTA (from temp_dir/{chr}/{chr}.fasta)
 *
 * Outputs:
 *   {chr}.cov.txt.gz  — per-base coverage signal
 *   {chr}.cov.bed.gz  — run-length encoded coverage BED
 */

nextflow.enable.dsl = 2

process HF_COVERAGE {

    label 'process_high'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploFill/${params.out}_tmp", mode: 'copy',
        saveAs: { filename -> "${chr_name}/${filename}" }

    input:
    tuple val(chr_name), val(chr_length), path(chr_fasta), path(bam), path(bai)

    output:
    tuple val(chr_name), path("${chr_name}.cov.txt.gz"), emit: signal
    tuple val(chr_name), path("${chr_name}.cov.bed.gz"),  emit: range_bed

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def tool      = params.coverage_tool ?: 'bedtools'
    def cmd = "python3 ${haplosync}/scripts/hapfill_coverage.py"
    cmd    += " -b ${bam}"
    cmd    += " -c ${chr_name}"
    cmd    += " -l ${chr_length}"
    cmd    += " -f ${chr_fasta}"
    cmd    += " -o ${chr_name}"
    cmd    += " -t ${tool}"

    """
    ${cmd}
    """
}
