/*
 * Workflow: PSEUDOMOLECULE_GENERATION
 *
 * Steps:
 *   1. BUILD_PATHS    — marker QC + DAG tiling path selection
 *   2. RECONSTRUCT    — AGP + FASTA + correspondence from tiling paths
 *   3. TRANSLATE      — coordinate translation (markers BED, legacy AGP, GFF3)
 *                       optional; runs in parallel when inputs are provided
 *   4. CHR_PAIR_QC    — per-chromosome Hap1 vs Hap2 overview reports (optional)
 *   5. REJECTED_QC    — per-unplaced-sequence QC reports (optional)
 *   6. HAPLODUP       — deduplication QC (optional, --run_haplodup)
 */

nextflow.enable.dsl = 2

include { BUILD_PATHS     } from '../modules/local/build_paths/main'
include { RECONSTRUCT     } from '../modules/local/reconstruct/main'
include { TRANSLATE       } from '../modules/local/translate/main'
include { CHR_PAIR_QC     } from '../modules/local/chr_pair_qc/main'
include { REJECTED_QC     } from '../modules/local/rejected_qc/main'
include { HAPLODUP_ALIGN  } from '../modules/local/haplodup_align/main'
include { HAPLODUP_GMAP   } from '../modules/local/haplodup_gmap/main'
include { HAPLODUP_REPORT } from '../modules/local/haplodup_report/main'

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
    // Step 4: Chromosome pair overview QC (optional)
    // -----------------------------------------------------------------------
    def run_chr_pair_qc = !params.skip_chr_pair_reports && !params.No2

    if (run_chr_pair_qc) {
        def agp_for_qc = RECONSTRUCT.out.hap1_agp
            .mix(RECONSTRUCT.out.hap2_agp.ifEmpty(Channel.empty()))
            .mix(RECONSTRUCT.out.un_agp)
            .collect()

        def markers_bed_qc = (run_translate && params.markers)
            ? TRANSLATE.out.markers_bed.ifEmpty([])
            : Channel.value([])
        def legacy_agp_qc  = (run_translate && params.input_agp)
            ? TRANSLATE.out.legacy_agp.ifEmpty([])
            : Channel.value([])

        CHR_PAIR_QC(
            RECONSTRUCT.out.correspondence,
            agp_for_qc,
            markers_bed_qc,
            legacy_agp_qc
        )
    }

    // -----------------------------------------------------------------------
    // Step 5: Unplaced sequence QC (optional)
    // -----------------------------------------------------------------------
    def run_rejected_qc = !params.skip_unplaced_qc && !params.No2

    if (run_rejected_qc) {
        def agp_for_rqc = RECONSTRUCT.out.hap1_agp
            .mix(RECONSTRUCT.out.hap2_agp.ifEmpty(Channel.empty()))
            .mix(RECONSTRUCT.out.un_agp)
            .collect()

        def fasta_for_rqc = RECONSTRUCT.out.hap1_fasta
            .mix(RECONSTRUCT.out.hap2_fasta.ifEmpty(Channel.empty()))
            .mix(RECONSTRUCT.out.un_fasta)
            .collect()

        def markers_bed_rqc = (run_translate && params.markers)
            ? TRANSLATE.out.markers_bed.ifEmpty([])
            : Channel.value([])
        def legacy_agp_rqc  = (run_translate && params.input_agp)
            ? TRANSLATE.out.legacy_agp.ifEmpty([])
            : Channel.value([])

        REJECTED_QC(
            BUILD_PATHS.out.unused_list,
            RECONSTRUCT.out.correspondence,
            fasta_for_rqc,
            agp_for_rqc,
            markers_bed_rqc,
            legacy_agp_rqc
        )
    }

    // -----------------------------------------------------------------------
    // Step 6: HaploDup (optional) — three parallel/sequential modules
    //   HAPLODUP_ALIGN  — nucmer pairwise alignments (compute-heavy)
    //   HAPLODUP_GMAP   — GMAP gene mapping (compute-heavy, parallel with ALIGN)
    //   HAPLODUP_REPORT — dotplots + HTML/PDF reports + hotspots + index
    // -----------------------------------------------------------------------
    if (params.run_haplodup) {
        // Shared AGP channel for HAPLODUP_REPORT
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

        // Step 6a: Pairwise nucmer alignments
        HAPLODUP_ALIGN(
            RECONSTRUCT.out.hap1_fasta,
            RECONSTRUCT.out.hap2_fasta,
            RECONSTRUCT.out.un_fasta,
            RECONSTRUCT.out.correspondence
        )

        // Step 6b: GMAP gene mapping (only when GFF3 annotation is provided;
        //          runs in parallel with HAPLODUP_ALIGN)
        def run_gmap = params.gff3 && !params.No2
        def gmap_gff3_ch = Channel.value([])
        if (run_gmap) {
            HAPLODUP_GMAP(
                RECONSTRUCT.out.hap1_fasta,
                RECONSTRUCT.out.hap2_fasta,
                RECONSTRUCT.out.un_fasta,
                RECONSTRUCT.out.correspondence,
                annotation_ch
            )
            gmap_gff3_ch = HAPLODUP_GMAP.out.gmap_gff3
        }

        // Step 6c: Reports — waits for ALIGN (and GMAP when applicable)
        HAPLODUP_REPORT(
            RECONSTRUCT.out.hap1_fasta,
            RECONSTRUCT.out.hap2_fasta,
            RECONSTRUCT.out.un_fasta,
            agp_ch,
            RECONSTRUCT.out.correspondence,
            markers_bed_ch,
            legacy_agp_ch,
            annotation_ch,
            HAPLODUP_ALIGN.out.delta_files.collect(),
            gmap_gff3_ch
        )
    }
}
