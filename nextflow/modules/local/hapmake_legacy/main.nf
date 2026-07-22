/*
 * Process: HM_MAKE_LEGACY
 *
 * HaploMake, second pass: ports a deeper "legacy" AGP (e.g. tracing back to
 * the original pre-HaploSplit contigs) through the same gap-fill structure
 * block, independently of HM_MAKE's own --hapmake_agp pass.
 *
 * HaploMake.py's -a/--agp flag only accepts a single AGP file per
 * invocation, so porting two distinct AGP layers (HaploSplit's own AGP,
 * handled by HM_MAKE, and a deeper legacy one, handled here) genuinely
 * requires two separate runs. This pass always uses --noprint since its
 * only purpose is producing the legacy AGP, not rebuilding the FASTA
 * (HM_MAKE already does that).
 *
 * Inputs:
 *   hap1_fasta      — Hap1 FASTA
 *   hap2_fasta      — Hap2 FASTA
 *   un_fasta        — Unplaced sequences FASTA (optional)
 *   structure_block — .structure.block from HF_FILL
 *
 * Outputs:
 *   {out}.contigs.legacy_structure.agp — legacy coordinate mapping
 *   {out}.contigs.structure.agp        — side artifact, same structure as HM_MAKE's own AGP
 */

nextflow.enable.dsl = 2

process HM_MAKE_LEGACY {

    label 'process_medium'

    conda "${projectDir}/nextflow/envs/haplosync.yml"

    publishDir "${params.outdir}/HaploMake", mode: 'copy'

    input:
    path hap1_fasta
    path hap2_fasta
    path un_fasta
    path structure_block

    output:
    path "versions.yml", emit: versions
    path "${params.out}.contigs.structure.agp",        emit: agp,        optional: true
    path "${params.out}.contigs.legacy_structure.agp", emit: legacy_agp, optional: true

    script:
    def fasta_list = un_fasta ? "${hap1_fasta},${hap2_fasta},${un_fasta}"
                               : "${hap1_fasta},${hap2_fasta}"
    def cmd = "hapmake.py"
    cmd    += " -f ${fasta_list}"
    cmd    += " -s ${structure_block}"
    cmd    += " -o ${params.out}.contigs"
    cmd    += " -a ${params.hapmake_legacy_agp}"
    cmd    += " --noprint"
    if (params.hapmake_gap)         cmd += " --gap ${params.hapmake_gap}"
    if (params.hapmake_skipoverlap) cmd += " --skipoverlap"

    """
    ${cmd}
    cat <<END_VERSIONS > versions.yml
"${task.process}":
    python3: \$(python3 --version | sed 's/Python //')
END_VERSIONS
    """
}
