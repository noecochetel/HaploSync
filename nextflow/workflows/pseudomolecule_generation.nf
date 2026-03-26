/*
 * Workflow: PSEUDOMOLECULE_GENERATION
 *
 * Steps:
 *   1. HAPLOSPLIT  — reconstruct haplotype pseudomolecules
 *   2. HAPLODUP    — deduplication QC (optional, --run_haplodup)
 *
 * After reviewing HaploDup output the user may manually curate
 * blacklists/exclusion files and re-run this workflow. Nextflow's
 * -resume mechanism will skip HaploSplit if its inputs are unchanged.
 */

nextflow.enable.dsl = 2

include { HAPLOSPLIT } from '../modules/local/haplosplit/main'
include { HAPLODUP   } from '../modules/local/haplodup/main'

workflow PSEUDOMOLECULE_GENERATION {

    // -----------------------------------------------------------------------
    // Step 1: HaploSplit
    // -----------------------------------------------------------------------
    HAPLOSPLIT()

    // -----------------------------------------------------------------------
    // Step 2: HaploDup (optional)
    // -----------------------------------------------------------------------
    if (params.run_haplodup) {
        HAPLODUP(
            HAPLOSPLIT.out.hap1_fasta,
            HAPLOSPLIT.out.hap2_fasta,
            HAPLOSPLIT.out.un_fasta,
            HAPLOSPLIT.out.agp,
            HAPLOSPLIT.out.correspondence,
            HAPLOSPLIT.out.markers_bed.ifEmpty([]),
            HAPLOSPLIT.out.legacy_agp.ifEmpty([]),
            HAPLOSPLIT.out.annotation.ifEmpty([])
        )
    }
}
