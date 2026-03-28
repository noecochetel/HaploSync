/*
 * Workflows: HAPLOSPLIT, QC, HAPLODUP, HAPLOSYNC_RECONSTRUCT_PM
 *
 * HAPLOSPLIT            — Steps 1–3: tiling path selection, pseudomolecule
 *                         reconstruction, coordinate translation.
 *                         Emits FASTA/AGP/correspondence/annotation channels.
 *
 * QC                    — Steps 4–5: chromosome pair overview reports and
 *                         unplaced sequence QC reports.
 *                         Controlled by --skip_chr_pair_reports /
 *                         --skip_unplaced_qc flags.
 *
 * HAPLODUP              — Steps 6a–6c: pairwise nucmer alignments, GMAP gene
 *                         mapping, dotplot + HTML/PDF reports.
 *
 * HAPLOSYNC_RECONSTRUCT_PM — Pipeline wrapper: HAPLOSPLIT → QC → HAPLODUP
 *                            (HAPLODUP only when --run_haplodup).
 *
 * Nextflow log naming produced by this structure:
 *   HAPLOSYNC_RECONSTRUCT_PM:HAPLOSPLIT:BUILD_PATHS
 *   HAPLOSYNC_RECONSTRUCT_PM:QC:CHR_PAIR_QC
 *   HAPLOSYNC_RECONSTRUCT_PM:HAPLODUP:HAPLODUP_ALIGN
 */

nextflow.enable.dsl = 2

include { BUILD_PATHS as HS_BUILD_PATHS } from '../modules/local/build_paths/main'
include { RECONSTRUCT as HS_RECONSTRUCT } from '../modules/local/reconstruct/main'
include { TRANSLATE   as HS_TRANSLATE   } from '../modules/local/translate/main'
include { CHR_PAIR    as QC_CHR_PAIR    } from '../modules/local/chr_pair_qc/main'
include { REJECTED    as QC_REJECTED    } from '../modules/local/rejected_qc/main'
include { ALIGN       as HD_ALIGN       } from '../modules/local/haplodup_align/main'
include { GMAP        as HD_GMAP        } from '../modules/local/haplodup_gmap/main'
include { REPORT      as HD_REPORT      } from '../modules/local/haplodup_report/main'

// ---------------------------------------------------------------------------
// Sub-workflow: HAPLOSPLIT
//   Steps 1–3: tiling paths → reconstruction → coordinate translation.
//   Emits all channels needed by QC and HAPLODUP.
// ---------------------------------------------------------------------------
workflow HAPLOSPLIT {

    main:

    // Step 1: Marker QC + DAG tiling path selection
    HS_BUILD_PATHS()

    // Step 2: AGP + FASTA + correspondence from tiling paths
    HS_RECONSTRUCT(
        HS_BUILD_PATHS.out.hap1_list,
        HS_BUILD_PATHS.out.hap2_list.ifEmpty([]),
        HS_BUILD_PATHS.out.un_list,
        HS_BUILD_PATHS.out.unused_list
    )

    // Step 3: Coordinate translation (only when inputs provided)
    def run_translate = params.markers || params.input_agp || params.gff3

    def markers_bed_out = Channel.value([])
    def legacy_agp_out  = Channel.value([])
    def annotation_out  = Channel.value([])

    if (run_translate) {
        HS_TRANSLATE(
            HS_RECONSTRUCT.out.hap1_agp,
            HS_RECONSTRUCT.out.hap2_agp.ifEmpty([]),
            HS_RECONSTRUCT.out.un_agp,
            HS_RECONSTRUCT.out.hap1_fasta,
            HS_RECONSTRUCT.out.hap2_fasta.ifEmpty([]),
            HS_RECONSTRUCT.out.un_fasta
        )
        if (params.markers)   markers_bed_out = HS_TRANSLATE.out.markers_bed.ifEmpty([])
        if (params.input_agp) legacy_agp_out  = HS_TRANSLATE.out.legacy_agp.ifEmpty([])
        if (params.gff3)      annotation_out  = HS_TRANSLATE.out.annotation.ifEmpty([])
    }

    emit:
    hap1_fasta     = HS_RECONSTRUCT.out.hap1_fasta
    hap2_fasta     = HS_RECONSTRUCT.out.hap2_fasta
    un_fasta       = HS_RECONSTRUCT.out.un_fasta
    correspondence = HS_RECONSTRUCT.out.correspondence
    hap1_agp       = HS_RECONSTRUCT.out.hap1_agp
    hap2_agp       = HS_RECONSTRUCT.out.hap2_agp
    un_agp         = HS_RECONSTRUCT.out.un_agp
    unused_list    = HS_BUILD_PATHS.out.unused_list
    markers_bed    = markers_bed_out
    legacy_agp     = legacy_agp_out
    annotation     = annotation_out
}

// ---------------------------------------------------------------------------
// Sub-workflow: QC
//   Steps 4–5: chromosome pair overview QC + unplaced sequence QC.
//   --skip_chr_pair_reports and --skip_unplaced_qc control what runs.
//   --No2 disables both (no Hap2 to compare against).
// ---------------------------------------------------------------------------
workflow QC {

    take:
    correspondence
    agp_ch
    fasta_ch
    unused_list
    markers_bed_ch
    legacy_agp_ch

    main:

    // Step 4: Per-chromosome Hap1 vs Hap2 overview reports
    if (!params.skip_chr_pair_reports && !params.No2) {
        QC_CHR_PAIR(
            correspondence,
            agp_ch,
            markers_bed_ch,
            legacy_agp_ch
        )
    }

    // Step 5: Per-unplaced-sequence QC reports
    if (!params.skip_unplaced_qc && !params.No2) {
        QC_REJECTED(
            unused_list,
            correspondence,
            fasta_ch,
            agp_ch,
            markers_bed_ch,
            legacy_agp_ch
        )
    }
}

// ---------------------------------------------------------------------------
// Sub-workflow: HAPLODUP
//   Steps 6a–6c: nucmer alignments, GMAP gene mapping, reports.
//   GMAP only runs when GFF3 annotation is provided (params.gff3) and
//   Hap2 is present (!params.No2).
// ---------------------------------------------------------------------------
workflow HAPLODUP {

    take:
    hap1_fasta
    hap2_fasta
    un_fasta
    correspondence
    agp_ch
    markers_bed_ch
    legacy_agp_ch
    annotation_ch

    main:

    // Step 6a: Pairwise nucmer alignments (compute-heavy)
    HD_ALIGN(hap1_fasta, hap2_fasta, un_fasta, correspondence)

    // Step 6b: GMAP gene mapping (parallel with HD_ALIGN)
    def run_gmap     = params.gff3 && !params.No2
    def gmap_gff3_ch = Channel.value([])

    if (run_gmap) {
        HD_GMAP(hap1_fasta, hap2_fasta, un_fasta, correspondence, annotation_ch)
        gmap_gff3_ch = HD_GMAP.out.gmap_gff3
    }

    // Step 6c: Reports (waits for ALIGN and optionally GMAP)
    HD_REPORT(
        hap1_fasta,
        hap2_fasta,
        un_fasta,
        agp_ch,
        correspondence,
        markers_bed_ch,
        legacy_agp_ch,
        annotation_ch,
        HD_ALIGN.out.delta_files.collect(),
        gmap_gff3_ch
    )
}

// ---------------------------------------------------------------------------
// Pipeline wrapper: HAPLOSYNC_RECONSTRUCT_PM
//   Runs HAPLOSPLIT → QC → optional HAPLODUP (--run_haplodup).
//   Produces Nextflow log names:
//     HAPLOSYNC_RECONSTRUCT_PM:HAPLOSPLIT:<PROCESS>
//     HAPLOSYNC_RECONSTRUCT_PM:QC:<PROCESS>
//     HAPLOSYNC_RECONSTRUCT_PM:HAPLODUP:<PROCESS>
// ---------------------------------------------------------------------------
workflow HAPLOSYNC_RECONSTRUCT_PM {

    HAPLOSPLIT()

    def agp_ch = HAPLOSPLIT.out.hap1_agp
        .mix(HAPLOSPLIT.out.hap2_agp.ifEmpty(Channel.empty()))
        .mix(HAPLOSPLIT.out.un_agp.ifEmpty(Channel.empty()))
        .collect()

    def fasta_ch = HAPLOSPLIT.out.hap1_fasta
        .mix(HAPLOSPLIT.out.hap2_fasta.ifEmpty(Channel.empty()))
        .mix(HAPLOSPLIT.out.un_fasta.ifEmpty(Channel.empty()))
        .collect()

    QC(
        HAPLOSPLIT.out.correspondence,
        agp_ch,
        fasta_ch,
        HAPLOSPLIT.out.unused_list,
        HAPLOSPLIT.out.markers_bed,
        HAPLOSPLIT.out.legacy_agp
    )

    if (params.run_haplodup) {
        HAPLODUP(
            HAPLOSPLIT.out.hap1_fasta,
            HAPLOSPLIT.out.hap2_fasta,
            HAPLOSPLIT.out.un_fasta,
            HAPLOSPLIT.out.correspondence,
            agp_ch,
            HAPLOSPLIT.out.markers_bed,
            HAPLOSPLIT.out.legacy_agp,
            HAPLOSPLIT.out.annotation
        )
    }
}
