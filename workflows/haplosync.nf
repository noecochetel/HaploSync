/*
 * Composable sub-workflows shared by the two HaploSync pipeline stages,
 * dispatched from the repo-root main.nf via --step.
 *
 * Reconstruct PM side (HAPLOSPLIT -> QC -> optional PM_HAPLODUP):
 *   HAPLOSPLIT               — Steps 1-3: tiling path selection, pseudomolecule
 *                               reconstruction, coordinate translation.
 *   QC                       — Steps 4-5: chromosome pair overview reports and
 *                               unplaced sequence QC reports.
 *   PM_HAPLODUP               — Steps 6a-6c: pairwise nucmer alignments, GMAP gene
 *                               mapping, dotplot + HTML/PDF reports.
 *   HAPLOSYNC_RECONSTRUCT_PM  — Pipeline wrapper: HAPLOSPLIT -> QC -> PM_HAPLODUP
 *                               (PM_HAPLODUP only when --run_haplodup).
 *
 * Gap fill side (HAPLOFILL -> optional HAPLOMAKE -> optional GF_HAPLODUP):
 *   HAPLOFILL                — Steps 1-6: setup, per-chromosome coverage, ploidy,
 *                               haplotype pairing, gap patching.
 *   HAPLOMAKE                 — Constructs new FASTA/AGP from the HaploFill
 *                               structure block.
 *   GF_HAPLODUP                — Duplication QC on the new gap-filled assembly.
 *   HAPLOSYNC_GAP_FILL        — Pipeline wrapper: HAPLOFILL -> HAPLOMAKE -> GF_HAPLODUP.
 *
 * PM_HAPLODUP and GF_HAPLODUP are the same underlying steps (HD_ALIGN/HD_GMAP/HD_REPORT)
 * applied to different inputs (HaploSplit output vs. HaploMake output) — kept as
 * two names only because both used to be called plain "HAPLODUP" in separate files;
 * now that they live in one file they need distinct names.
 */

nextflow.enable.dsl = 2

include { BUILD_PATHS as HS_BUILD_PATHS } from '../nextflow/modules/local/build_paths/main'
include { RECONSTRUCT as HS_RECONSTRUCT } from '../nextflow/modules/local/reconstruct/main'
include { TRANSLATE   as HS_TRANSLATE   } from '../nextflow/modules/local/translate/main'
include { CHR_PAIR    as QC_CHR_PAIR    } from '../nextflow/modules/local/chr_pair_qc/main'
include { REJECTED    as QC_REJECTED    } from '../nextflow/modules/local/rejected_qc/main'
include { ALIGN       as HD_ALIGN       } from '../nextflow/modules/local/haplodup_align/main'
include { GMAP        as HD_GMAP        } from '../nextflow/modules/local/haplodup_gmap/main'
include { REPORT      as HD_REPORT      } from '../nextflow/modules/local/haplodup_report/main'
include { HF_SETUP    as HF_SETUP       } from '../nextflow/modules/local/hapfill_setup/main'
include { HF_COVERAGE as HF_COVERAGE    } from '../nextflow/modules/local/hapfill_coverage/main'
include { HF_PLOIDY   as HF_PLOIDY      } from '../nextflow/modules/local/hapfill_ploidy/main'
include { HF_PAIR     as HF_PAIR        } from '../nextflow/modules/local/hapfill_pair/main'
include { HF_FILL     as HF_FILL        } from '../nextflow/modules/local/hapfill_fill/main'
include { HM_MAKE     as HM_MAKE        } from '../nextflow/modules/local/hapmake/main'
include { HM_MAKE_LEGACY as HM_MAKE_LEGACY } from '../nextflow/modules/local/hapmake_legacy/main'

// ---------------------------------------------------------------------------
// Sub-workflow: HAPLOSPLIT
//   Steps 1-3: tiling paths -> reconstruction -> coordinate translation.
//   Emits all channels needed by QC and PM_HAPLODUP.
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
//   Steps 4-5: chromosome pair overview QC + unplaced sequence QC.
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
// Sub-workflow: PM_HAPLODUP
//   Steps 6a-6c: nucmer alignments, GMAP gene mapping, reports.
//   GMAP only runs when GFF3 annotation is provided (params.gff3) and
//   Hap2 is present (!params.No2).
// ---------------------------------------------------------------------------
workflow PM_HAPLODUP {

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
//   Runs HAPLOSPLIT -> QC -> optional PM_HAPLODUP (--run_haplodup).
//   Nextflow log names:
//     HAPLOSYNC_RECONSTRUCT_PM:HAPLOSPLIT:<PROCESS>
//     HAPLOSYNC_RECONSTRUCT_PM:QC:<PROCESS>
//     HAPLOSYNC_RECONSTRUCT_PM:PM_HAPLODUP:<PROCESS>
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
        PM_HAPLODUP(
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

// ---------------------------------------------------------------------------
// Sub-workflow: HAPLOFILL
//   Steps 1-6: setup -> per-chromosome coverage (scattered) -> ploidy ->
//              haplotype pairing -> gap filling.
// ---------------------------------------------------------------------------
workflow HAPLOFILL {

    take:
    hap1_fasta
    hap2_fasta
    un_fasta
    correspondence
    repeats
    bam_hap1
    bam_hap1_bai
    bam_hap2
    bam_hap2_bai

    main:

    // Step 1: Setup — sequence splitting, pairing, repeats, gap detection
    HF_SETUP(hap1_fasta, hap2_fasta, un_fasta, correspondence, repeats)

    // Build per-chromosome channel AFTER HF_SETUP completes, by scanning the
    // temp_dir that HF_SETUP produces. Each entry:
    //   tuple(chr_name, chr_length, chr_fasta, bam, bai)
    // BAM assignment is determined by which input FASTA (hap1 or hap2) the
    // sequence comes from. HaploFill requires coverage for ALL sequences in
    // both FASTAs (not just paired ones), so we must not filter by correspondence.
    def chr_channel = HF_SETUP.out.temp_dir
        .flatMap { dir ->
            def hap1_seqs = [] as Set
            def hap2_seqs = [] as Set
            file(params.hapfill_hap1).eachLine { line ->
                if (line.startsWith('>'))
                    hap1_seqs.add(line.substring(1).trim().split(/\s+/)[0])
            }
            file(params.hapfill_hap2).eachLine { line ->
                if (line.startsWith('>'))
                    hap2_seqs.add(line.substring(1).trim().split(/\s+/)[0])
            }
            def result = []
            dir.eachDir { chrDir ->
                if (chrDir.name == 'unplaced' || chrDir.name.startsWith('Pair_')) return
                def fasta = file("${chrDir}/${chrDir.name}.fasta")
                if (!fasta.exists()) return
                def bam
                def bai
                if (hap1_seqs.contains(chrDir.name)) {
                    bam = file(params.hapfill_b1)
                    bai = file(params.hapfill_b1 + '.bai')
                } else if (hap2_seqs.contains(chrDir.name)) {
                    bam = file(params.hapfill_b2)
                    bai = file(params.hapfill_b2 + '.bai')
                } else {
                    return
                }
                def length = 0L
                fasta.eachLine { ln -> if (!ln.startsWith('>')) length += ln.trim().size() }
                result << tuple(chrDir.name, length, fasta, bam, bai)
            }
            result
        }

    // Step 2: Coverage — scattered one job per chromosome
    // maxForks controls parallel jobs (RAM budget: ~4-6 GB/job with bedtools,
    // ~0.5 GB/job with mosdepth). Set in nextflow.config:
    //   process { withName: 'HF_COVERAGE' { maxForks = 8 } }
    HF_COVERAGE(chr_channel)

    // Step 3: Ploidy — gather all signal files, compute global median
    HF_PLOIDY(
        HF_SETUP.out.temp_dir,
        HF_COVERAGE.out.signal.map { it[1] }.collect(),
        HF_COVERAGE.out.range_bed.map { it[1] }.collect()
    )

    // Steps 4-5: Haplotype pairing + gap region preparation
    HF_PAIR(HF_PLOIDY.out.temp_dir)

    // Step 6: Gap filling — produces .structure.block
    // Stage hap1/hap2/unplaced/correspondence/repeats so HaploFill step 6.2
    // (map_nucmer_unplaced_on_target) can find the files by basename in the work dir.
    HF_FILL(HF_PAIR.out.temp_dir, hap1_fasta, hap2_fasta, un_fasta, correspondence, repeats)

    emit:
    structure_block = HF_FILL.out.structure_block
    findings        = HF_FILL.out.findings
}

// ---------------------------------------------------------------------------
// Sub-workflow: HAPLOMAKE
//   Constructs new FASTA/AGP from HaploFill structure block.
// ---------------------------------------------------------------------------
workflow HAPLOMAKE {

    take:
    hap1_fasta
    hap2_fasta
    un_fasta
    structure_block

    main:

    HM_MAKE(hap1_fasta, hap2_fasta, un_fasta, structure_block)

    emit:
    fasta      = HM_MAKE.out.fasta
    agp        = HM_MAKE.out.agp
    legacy_agp = HM_MAKE.out.legacy_agp
    bed        = HM_MAKE.out.bed
    gff3       = HM_MAKE.out.gff3
}

// ---------------------------------------------------------------------------
// Sub-workflow: GF_HAPLODUP (gap-fill context)
//   Duplication QC on the new gap-filled assembly (optional).
// ---------------------------------------------------------------------------
workflow GF_HAPLODUP {

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

    HD_ALIGN(hap1_fasta, hap2_fasta, un_fasta, correspondence)

    def run_gmap     = params.hapmake_gff3 && !params.No2
    def gmap_gff3_ch = Channel.value([])

    if (run_gmap) {
        HD_GMAP(hap1_fasta, hap2_fasta, un_fasta, correspondence, annotation_ch)
        gmap_gff3_ch = HD_GMAP.out.gmap_gff3
    }

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
// Pipeline wrapper: HAPLOSYNC_GAP_FILL
//   HAPLOFILL -> optional HAPLOMAKE -> optional GF_HAPLODUP
//   --run_haplomake  : run HaploMake after gap filling
//   --run_haplodup   : run HaploDup on the new assembly (implies --run_haplomake)
// ---------------------------------------------------------------------------
workflow HAPLOSYNC_GAP_FILL {

    def hap1_fasta     = Channel.fromPath(params.hapfill_hap1)
    def hap2_fasta     = Channel.fromPath(params.hapfill_hap2)
    def un_fasta       = params.hapfill_unplaced
                             ? Channel.fromPath(params.hapfill_unplaced)
                             : Channel.value([])
    def correspondence = Channel.fromPath(params.hapfill_correspondence)
    def repeats        = Channel.fromPath(params.hapfill_repeats)
    def bam_hap1       = Channel.fromPath(params.hapfill_b1)
    def bam_hap1_bai   = Channel.fromPath(params.hapfill_b1 + '.bai')
    def bam_hap2       = Channel.fromPath(params.hapfill_b2)
    def bam_hap2_bai   = Channel.fromPath(params.hapfill_b2 + '.bai')

    HAPLOFILL(
        hap1_fasta, hap2_fasta, un_fasta,
        correspondence, repeats,
        bam_hap1, bam_hap1_bai,
        bam_hap2, bam_hap2_bai
    )

    // HaploMake runs if explicitly requested or implicitly required by HaploDup
    def need_haplomake = params.run_haplomake || params.run_haplodup

    if (need_haplomake) {
        HAPLOMAKE(
            hap1_fasta,
            hap2_fasta,
            un_fasta,
            HAPLOFILL.out.structure_block
        )

        if (params.hapmake_legacy_agp) {
            HM_MAKE_LEGACY(
                hap1_fasta,
                hap2_fasta,
                un_fasta,
                HAPLOFILL.out.structure_block
            )
        }

        if (params.run_haplodup) {
            def agp_ch = HAPLOMAKE.out.agp.collect()

            def annotation_ch = params.hapmake_gff3
                ? HAPLOMAKE.out.gff3.ifEmpty([])
                : Channel.value([])

            GF_HAPLODUP(
                HAPLOMAKE.out.fasta,
                HAPLOMAKE.out.fasta,
                un_fasta,
                correspondence,
                agp_ch,
                HAPLOMAKE.out.bed.ifEmpty([]),
                HAPLOMAKE.out.legacy_agp.ifEmpty([]),
                annotation_ch
            )
        }
    }
}
