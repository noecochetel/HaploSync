/*
 * Process: HM_MAKE
 *
 * HaploMake: construct new pseudomolecule sequences from a HaploFill
 * structure block file.
 *
 * Reads the .structure.block output of HF_FILL and builds:
 *   - New FASTA with gap-filled sequences
 *   - AGP metadata file
 *   - Optionally: legacy coordinate mapping, translated annotation
 *
 * Inputs:
 *   hap1_fasta      — Hap1 FASTA
 *   hap2_fasta      — Hap2 FASTA
 *   un_fasta        — Unplaced sequences FASTA (optional)
 *   structure_block — .structure.block from HF_FILL
 *
 * Outputs:
 *   {out}.fasta              — new gap-filled pseudomolecule sequences
 *   {out}.structure.agp      — AGP for new assembly
 *   {out}.legacy_structure.agp — legacy coordinate mapping (if --hapmake_agp)
 */

nextflow.enable.dsl = 2

process HM_MAKE {

    label 'process_medium'

    conda "${projectDir}/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploMake", mode: 'copy'

    input:
    path hap1_fasta
    path hap2_fasta
    path un_fasta
    path structure_block

    output:
    path "${params.out}.fasta",         emit: fasta
    path "${params.out}.structure.agp", emit: agp,           optional: true
    path "${params.out}.legacy_structure.agp", emit: legacy_agp, optional: true

    script:
    def haplosync = params.haplosync_dir ?: "${projectDir}/.."
    def fasta_list = un_fasta ? "${hap1_fasta},${hap2_fasta},${un_fasta}"
                               : "${hap1_fasta},${hap2_fasta}"
    def cmd = "python3 ${haplosync}/scripts/hapmake.py"
    cmd    += " -f ${fasta_list}"
    cmd    += " -s ${structure_block}"
    cmd    += " -o ${params.out}"
    if (params.hapmake_prefix)     cmd += " -p ${params.hapmake_prefix}"
    if (params.hapmake_agp)        cmd += " -a ${params.hapmake_agp}"
    if (params.hapmake_gff3)       cmd += " --gff3 ${params.hapmake_gff3}"
    if (params.hapmake_bed)        cmd += " -b ${params.hapmake_bed}"
    if (params.hapmake_gap)        cmd += " --gap ${params.hapmake_gap}"
    if (params.hapmake_skipoverlap) cmd += " --skipoverlap"
    if (params.hapmake_noagp)      cmd += " --noagp"
    if (params.hapmake_unplaced)   cmd += " -u ${params.hapmake_unplaced}"

    """
    ${cmd}
    """
}
