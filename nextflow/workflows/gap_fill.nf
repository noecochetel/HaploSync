/*
 * Workflow: HAPLOSYNC_GAP_FILL
 *
 * Gap-filling pipeline: HaploFill → HaploMake (optional) → HaploDup (optional).
 *
 * Sub-workflows:
 *   HAPFILL   — Steps 1–6: coverage extraction, ploidy classification,
 *               haplotype pairing, gap patching. Produces .structure.block.
 *               Always runs.
 *   HAPMAKE   — Constructs new FASTA/AGP from .structure.block.
 *               Optional (--run_haplomake). Also runs when --run_haplodup
 *               is set (HaploDup requires the new assembly).
 *   HAPLODUP  — Duplication QC on the new assembly (--run_haplodup).
 *               Requires --run_haplomake (or implies it automatically).
 *
 * HF_COVERAGE is scattered one job per chromosome — all chromosomes of
 * both haplotypes run in parallel. The chromosome list is built dynamically
 * from HF_SETUP's output (flatMap on temp_dir) using the correspondence file
 * to assign each chromosome to the correct BAM. Cap concurrency with:
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
 *   HAPLOSYNC_GAP_FILL:HAPMAKE:HM_MAKE       (--run_haplomake or --run_haplodup)
 *   HAPLOSYNC_GAP_FILL:HAPLODUP:HD_ALIGN     (--run_haplodup)
 *   HAPLOSYNC_GAP_FILL:HAPLODUP:HD_GMAP      (--run_haplodup + gff3)
 *   HAPLOSYNC_GAP_FILL:HAPLODUP:HD_REPORT    (--run_haplodup)
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

    // Build per-chromosome channel AFTER HF_SETUP completes, by scanning the
    // temp_dir that HF_SETUP produces. Each entry:
    //   tuple(chr_name, chr_length, chr_fasta, bam, bai)
    // Hap assignment is read directly from the correspondence file.
    def chr_channel = HF_SETUP.out.temp_dir
        .flatMap { dir ->
            def hap_map = [:]
            file(params.hapfill_correspondence).eachLine { line ->
                def parts = line.trim().split('\t')
                if (parts.size() >= 2) {
                    hap_map[parts[0]] = '1'
                    hap_map[parts[1]] = '2'
                }
            }
            def result = []
            dir.eachDir { chrDir ->
                if (chrDir.name == 'unplaced' || chrDir.name.startsWith('Pair_')) return
                def fasta = file("${chrDir}/${chrDir.name}.fasta")
                if (!fasta.exists()) return
                def hap = hap_map[chrDir.name]
                if (!hap) return
                def length = 0L
                fasta.eachLine { ln -> if (!ln.startsWith('>')) length += ln.trim().size() }
                def bam = hap == '1' ? file(params.hapfill_b1) : file(params.hapfill_b2)
                def bai = hap == '1' ? file(params.hapfill_b1 + '.bai') : file(params.hapfill_b2 + '.bai')
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
//   HAPFILL → optional HAPMAKE → optional HAPLODUP
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

    HAPFILL(
        hap1_fasta, hap2_fasta, un_fasta,
        correspondence, repeats,
        bam_hap1, bam_hap1_bai,
        bam_hap2, bam_hap2_bai
    )

    // HaploMake runs if explicitly requested or implicitly required by HaploDup
    def need_haplomake = params.run_haplomake || params.run_haplodup

    if (need_haplomake) {
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
}
