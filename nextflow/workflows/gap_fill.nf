/*
 * Workflow: HAPLOSYNC_GAP_FILL
 *
 * Gap-filling pipeline: HaploFill → HaploMake → HaploDup (optional).
 *
 * Sub-workflows:
 *   HAPFILL   — Steps 1–6: coverage extraction, ploidy classification,
 *               haplotype pairing, gap patching. Produces .structure.block.
 *   HAPMAKE   — Constructs new FASTA/AGP from .structure.block.
 *   HAPLODUP  — Duplication QC on the new assembly (optional, --run_haplodup).
 *
 * HF_COVERAGE is scattered one job per chromosome — all chromosomes of
 * both haplotypes run in parallel. Cap concurrency with:
 *   process { withName: 'HF_COVERAGE' { maxForks = 8 } }   // ~48 GB RAM
 *
 * Coverage tool (params.coverage_tool):
 *   bedtools  — default, validated, ~40 min/chromosome
 *   mosdepth  — opt-in, 20-50x faster, same output format
 *
 * Log names produced:
 *   HAPLOSYNC_GAP_FILL:HAPFILL:HF_SETUP
 *   HAPLOSYNC_GAP_FILL:HAPFILL:HF_COVERAGE
 *   HAPLOSYNC_GAP_FILL:HAPFILL:HF_PLOIDY
 *   HAPLOSYNC_GAP_FILL:HAPFILL:HF_PAIR
 *   HAPLOSYNC_GAP_FILL:HAPFILL:HF_FILL
 *   HAPLOSYNC_GAP_FILL:HAPMAKE:HM_MAKE
 *   HAPLOSYNC_GAP_FILL:HAPLODUP:HD_ALIGN
 *   HAPLOSYNC_GAP_FILL:HAPLODUP:HD_GMAP
 *   HAPLOSYNC_GAP_FILL:HAPLODUP:HD_REPORT
 */

nextflow.enable.dsl = 2

include { HF_SETUP    as HF_SETUP    } from '../modules/local/hapfill_setup/main'
include { HF_COVERAGE as HF_COVERAGE } from '../modules/local/hapfill_coverage/main'
include { HF_PLOIDY   as HF_PLOIDY   } from '../modules/local/hapfill_ploidy/main'
include { HF_PAIR     as HF_PAIR     } from '../modules/local/hapfill_pair/main'
include { HF_FILL     as HF_FILL     } from '../modules/local/hapfill_fill/main'
include { HM_MAKE     as HM_MAKE     } from '../modules/local/hapmake/main'
include { ALIGN       as HD_ALIGN    } from '../modules/local/haplodup_align/main'
include { GMAP        as HD_GMAP     } from '../modules/local/haplodup_gmap/main'
include { REPORT      as HD_REPORT   } from '../modules/local/haplodup_report/main'

// ---------------------------------------------------------------------------
// Sub-workflow: HAPFILL
//   Steps 1–6: setup → per-chromosome coverage (scattered) → ploidy →
//              haplotype pairing → gap filling.
// ---------------------------------------------------------------------------
workflow HAPFILL {

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

    // Prepare per-chromosome channel from temp_dir
    // Each entry: [chr_name, chr_length, chr_fasta, haplotype_bam, bai]
    // Built from params.hapfill_chr_list: list of [chr, length, hap(1|2)] tuples
    def chr_channel = Channel
        .fromPath("${params.outdir}/HaploFill/${params.out}_tmp/**/*.fasta",
                   followLinks: true)
        .map { fasta ->
            def chr_name   = fasta.baseName
            def chr_length = params.hapfill_chr_lengths[chr_name]
            def hap        = params.hapfill_chr_hap[chr_name]     // "1" or "2"
            def bam        = hap == "1" ? bam_hap1 : bam_hap2
            def bai        = hap == "1" ? bam_hap1_bai : bam_hap2_bai
            tuple(chr_name, chr_length, fasta, bam, bai)
        }

    // Step 2: Coverage — scattered one job per chromosome
    // maxForks controls parallel jobs (RAM budget: ~4-6 GB/job with bedtools,
    // ~0.5 GB/job with mosdepth). Set in nextflow.config:
    //   process { withName: 'HF_COVERAGE' { maxForks = 8 } }
    HF_COVERAGE(
        chr_channel.map { it[3] },   // bam
        chr_channel.map { it[4] },   // bai
        HF_SETUP.out.temp_dir,
        chr_channel.map { it[0] },   // chr_name
        chr_channel.map { it[1] },   // chr_length
        chr_channel.map { it[2] }    // chr_fasta
    )

    // Step 3: Ploidy — gather all signal files, compute global median
    HF_PLOIDY(
        HF_SETUP.out.temp_dir,
        HF_COVERAGE.out.signal.map { it[1] }.collect(),
        HF_COVERAGE.out.range_bed.map { it[1] }.collect()
    )

    // Steps 4-5: Haplotype pairing + gap region preparation
    HF_PAIR(HF_PLOIDY.out.temp_dir)

    // Step 6: Gap filling — produces .structure.block
    HF_FILL(HF_PAIR.out.temp_dir)

    emit:
    structure_block = HF_FILL.out.structure_block
    findings        = HF_FILL.out.findings
}

// ---------------------------------------------------------------------------
// Sub-workflow: HAPMAKE
//   Constructs new FASTA/AGP from HaploFill structure block.
// ---------------------------------------------------------------------------
workflow HAPMAKE {

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
}

// ---------------------------------------------------------------------------
// Sub-workflow: HAPLODUP (gap-fill context)
//   Duplication QC on the new gap-filled assembly (optional).
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

    HD_ALIGN(hap1_fasta, hap2_fasta, un_fasta, correspondence)

    def run_gmap     = params.gff3 && !params.No2
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
//   HAPFILL → HAPMAKE → optional HAPLODUP
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

    HAPFILL(
        hap1_fasta, hap2_fasta, un_fasta,
        correspondence, repeats,
        bam_hap1, bam_hap1_bai,
        bam_hap2, bam_hap2_bai
    )

    HAPMAKE(
        hap1_fasta,
        hap2_fasta,
        un_fasta,
        HAPFILL.out.structure_block
    )

    if (params.run_haplodup) {
        def agp_ch = HAPMAKE.out.agp.collect()

        HAPLODUP(
            HAPMAKE.out.fasta,
            HAPMAKE.out.fasta,
            un_fasta,
            correspondence,
            agp_ch,
            Channel.value([]),
            HAPMAKE.out.legacy_agp.ifEmpty([]),
            Channel.value([])
        )
    }
}
