/*
 * Workflow: PSEUDOMOLECULE_GENERATION
 *
 * Steps:
 *   1. BUILD_PATHS  — marker QC + DAG tiling path selection
 *   2. RECONSTRUCT  — AGP + FASTA + correspondence from tiling paths
 *   3. TRANSLATE    — coordinate translation (markers BED, legacy AGP, GFF3)
 *                     optional; runs in parallel when inputs are provided
 *   4. HAPLODUP     — deduplication QC (optional, --run_haplodup)
 */

nextflow.enable.dsl = 2

include { BUILD_PATHS } from '../modules/local/build_paths/main'
include { RECONSTRUCT } from '../modules/local/reconstruct/main'
include { TRANSLATE   } from '../modules/local/translate/main'
include { HAPLODUP    } from '../modules/local/haplodup/main'

workflow HAPLOSPLIT {

    // -----------------------------------------------------------------------
    // Step 1: Build tiling paths (marker QC + DAG path selection)
    // -----------------------------------------------------------------------
    BUILD_PATHS()

    // -----------------------------------------------------------------------
    // Step 2: Reconstruct pseudomolecules (AGP + FASTA + correspondence)
    // -----------------------------------------------------------------------
    RECONSTRUCT(
        BUILD_PATHS.out.hap1_list,
        BUILD_PATHS.out.hap2_list.ifEmpty([]),
        BUILD_PATHS.out.un_list,
        BUILD_PATHS.out.unused_list
    )

    // -----------------------------------------------------------------------
    // Step 3: Coordinate translation (optional, only if inputs provided)
    // -----------------------------------------------------------------------
    def run_translate = params.markers || params.input_agp || params.gff3

    if (run_translate) {
        // Collect all three AGP files for the translation AGP lookup
        def all_agp = RECONSTRUCT.out.hap1_agp
            .mix(RECONSTRUCT.out.hap2_agp.ifEmpty(Channel.empty()))
            .mix(RECONSTRUCT.out.un_agp)
            .collect()

        TRANSLATE(
            RECONSTRUCT.out.hap1_agp,
            RECONSTRUCT.out.hap2_agp.ifEmpty([]),
            RECONSTRUCT.out.un_agp,
            RECONSTRUCT.out.hap1_fasta,
            RECONSTRUCT.out.hap2_fasta.ifEmpty([]),
            RECONSTRUCT.out.un_fasta
        )
    }

    // -----------------------------------------------------------------------
    // Step 4: HaploDup (optional)
    // -----------------------------------------------------------------------
    if (params.run_haplodup) {
        // Collect Hap1 + Hap2 + Un AGP into a single list for HAPLODUP
        def agp_ch = RECONSTRUCT.out.hap1_agp
            .mix(RECONSTRUCT.out.hap2_agp)
            .mix(RECONSTRUCT.out.un_agp)
            .collect()

        def markers_bed_ch = (run_translate && params.markers)
            ? TRANSLATE.out.markers_bed.ifEmpty([])
            : Channel.value([])
        def legacy_agp_ch  = (run_translate && params.input_agp)
            ? TRANSLATE.out.legacy_agp.ifEmpty([])
            : Channel.value([])
        def annotation_ch  = (run_translate && params.gff3)
            ? TRANSLATE.out.annotation.ifEmpty([])
            : Channel.value([])

        HAPLODUP(
            RECONSTRUCT.out.hap1_fasta,
            RECONSTRUCT.out.hap2_fasta,
            RECONSTRUCT.out.un_fasta,
            agp_ch,
            RECONSTRUCT.out.correspondence,
            markers_bed_ch,
            legacy_agp_ch,
            annotation_ch
        )
    }
}
