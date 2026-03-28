#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HAPLOSYNC_RECONSTRUCT_PM              } from './workflows/pseudomolecule_generation'
include { ALIGN    as HD_ALIGN                  } from './modules/local/haplodup_align/main'
include { GMAP     as HD_GMAP                   } from './modules/local/haplodup_gmap/main'
include { REPORT   as HD_REPORT                 } from './modules/local/haplodup_report/main'
include { CHR_PAIR as QC_CHR_PAIR               } from './modules/local/chr_pair_qc/main'
include { REJECTED as QC_REJECTED               } from './modules/local/rejected_qc/main'

// --------------------------------------------------------------------------
// Help messages
// --------------------------------------------------------------------------
def helpHaploSplit() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                    HaploSync — HaploSplit v2.0                      ║
    ║          Haplotype-resolved pseudomolecule assembly                 ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf [options]
        nextflow run main.nf -params-file params.yml

    ── Required ────────────────────────────────────────────────────────────
        --input_fasta       Draft assembly FASTA
        Provide (--markers AND --markers_map) and/or --guide_genome

    ── Genetic map ─────────────────────────────────────────────────────────
        --markers           Marker hits on input sequences (BED, 4 columns)
        --markers_map       Marker genetic map (chr, position, marker_id)

    ── Reference genome ────────────────────────────────────────────────────
        --guide_genome      Guide/reference genome FASTA
        --local_alignment   Existing PAF file, skips minimap2 (-l)
        --run_alignment     Run minimap2 to generate alignment (--align)
        --mapping_tool      Minimap2 options string      [default: '--cs -x asm20 -r 1000']
        --distance1         Max Hap1 sequence gap (bp)   [default: 2,000,000]
        --distance2         Max Hap2 sequence gap (bp)   [default: 4,000,000]
        --hitgap            Max gap to merge hits (bp)   [default: 100,000]
        --reuse_intermediate  Reuse existing nucmer files [default: false]

    ── Optional inputs ─────────────────────────────────────────────────────
        --input_agp         Input AGP structure file
        --gff3              Gene annotation GFF3
        --exclusion         Mutually exclusive sequence pairs (TSV)
        --known             Sequences known to be in same haplotype (TSV)
        --alternative_groups  Alternative haplotype sequence pairs (TSV)
        --Require1          Sequences required in Hap1 (TSV)
        --Require2          Sequences required in Hap2 (TSV)
        --Blacklist1        Blacklisted sequences for Hap1 (TSV)
        --Blacklist2        Blacklisted sequences for Hap2 (TSV)
        --input_groups      Sequence grouping file
        --legacy_groups     Legacy component group file
        --path1             Pre-defined Hap1 tiling path, skips DAG (-1)
        --path2             Pre-defined Hap2 tiling path, skips DAG (-2)

    ── Assembly behaviour ───────────────────────────────────────────────────
        --No2               Skip Hap2 reconstruction     [default: false]
        --gap               Gap size in bp               [default: 1000]
        --filter_hits       Remove seqs with intra-seq marker duplication
        --extended_region   Extend seq association to all chr markers
        --conflict_resolution  Conflict resolution: exit|ignore|release
        --allow_rearrangements  Allow rearrangements between haplotypes
        --required_as_path  Treat --R1/--R2 as full tiling paths

    ── QC ──────────────────────────────────────────────────────────────────
        --skip_chimeric_qc  Skip chimeric sequence QC    [default: false]
        --disable_marker_ploidy_check
        --only_markers      Limit rejected-seq QC to marker-bearing seqs
        --avoid_rejected_qc  Skip QC of chr-assigned but unplaced seqs
        --skip_chr_pair_reports  Skip chr pair overview reports [default: false]
        --skip_unplaced_qc  Skip unplaced sequence QC reports  [default: false]

    ── Output ──────────────────────────────────────────────────────────────
        --out               Output files prefix          [default: out]
        --prefix            Sequence ID prefix           [default: NEW]
        --outdir            Results directory            [default: results]

    ── HaploDup ────────────────────────────────────────────────────────────
        --run_haplodup      Run HaploDup after reconstruction [default: false]
        --haplodup_opts     Extra HaploDup flags as a quoted string

        For full HaploDup options: nextflow run main.nf -entry HAPLODUP --help

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── Profiles ────────────────────────────────────────────────────────────
        -profile standard   Local execution
        -profile conda      Local execution with conda
        -profile mamba      Local execution with mamba/micromamba
        -profile hpc        SLURM execution with mamba/micromamba

    ── Examples ────────────────────────────────────────────────────────────
        # Genetic map mode
        nextflow run main.nf -profile mamba \\
            --input_fasta assembly.fasta \\
            --markers markers.bed \\
            --markers_map genetic_map.tsv \\
            --out myproject --outdir results

        # Reference genome mode
        nextflow run main.nf -profile mamba \\
            --input_fasta assembly.fasta \\
            --guide_genome reference.fasta --run_alignment \\
            --out myproject --outdir results

        # Full pipeline with HaploDup
        nextflow run main.nf -profile mamba \\
            --input_fasta assembly.fasta \\
            --markers markers.bed --markers_map genetic_map.tsv \\
            --run_haplodup \\
            --out myproject --outdir results

        # Using a params file
        nextflow run main.nf -profile mamba -params-file params.yml

        # Resume after manual curation
        nextflow run main.nf -profile mamba -resume -params-file params.yml

    """.stripIndent()
}

def helpHaploDup() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                    HaploSync — HaploDup v2.0                        ║
    ║              Haplotype duplication QC and reporting                 ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf -entry HAPLODUP [options]
        nextflow run main.nf -entry HAPLODUP -params-file params.yml

    HaploDup reads HaploSplit outputs automatically from --outdir/HaploSplit/
    using the --out prefix. Run HaploSplit first, or point to an existing
    HaploSplit output directory.

    The following files are discovered automatically when present:
        {out}.markers.bed           — marker positions in pseudomolecule space
        {out}.legacy_structure.agp  — legacy AGP structure
        {out}.annotation.gff3       — gene annotation

    ── Required ────────────────────────────────────────────────────────────
        --out               Output prefix (must match HaploSplit run)
                            [default: out]
        --outdir            Results directory (must match HaploSplit run)
                            [default: results]

    ── Optional inputs ─────────────────────────────────────────────────────
        --reference             Reference genome for dotplots (-r)
        --markers_map           Marker genetic map for QC (--markers_map)
        --input_groups          Sequence grouping file (--input_groups)
        --legacy_groups         Legacy grouping file (--legacy_groups)
        --functional_annotation Functional annotation per transcript (-a)

    ── Alignment thresholds ─────────────────────────────────────────────────
        --hit_identity          Min genome mapping hit identity  [default: 90]
        --hit_coverage          Min genome mapping hit length    [default: 3000]
        --gene_identity         Min gene mapping identity        [default: 95]
        --gene_coverage         Min gene mapping coverage        [default: 95]
        --unbalanced_ratio      Gene count ratio threshold       [default: 0.33]

    ── Gene mapping ─────────────────────────────────────────────────────────
        --haplodup_feature      GFF feature type [CDS|mRNA]     [default: CDS]
        --haplodup_window       Window size for unbalanced gene search [default: 10]
        --haplodup_allowed      Allowed unbalanced genes/window  [default: 5]

    ── Reuse / skip ─────────────────────────────────────────────────────────
        --reuse_mappings        Reuse existing genome alignments [default: false]
        --reuse_dotplots        Reuse existing dotplots          [default: false]
        --reuse_gmap            Reuse existing GMAP mappings     [default: false]
        --skip_dotplots_by_chr  Skip individual chr-vs-chr dotplots [default: false]
        --only_paired_dotplots  Only generate matched-pair dotplots [default: false]

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── Profiles ────────────────────────────────────────────────────────────
        -profile standard   Local execution
        -profile conda      Local execution with conda
        -profile mamba      Local execution with mamba/micromamba
        -profile hpc        SLURM execution with mamba/micromamba

    ── Examples ────────────────────────────────────────────────────────────
        # Run HaploDup on an existing HaploSplit output
        nextflow run main.nf -entry HAPLODUP -profile mamba \\
            --out myproject --outdir results

        # With reference genome for dotplots
        nextflow run main.nf -entry HAPLODUP -profile mamba \\
            --out myproject --outdir results \\
            --reference reference.fasta

        # Using a params file
        nextflow run main.nf -entry HAPLODUP -profile mamba \\
            -params-file params.yml
    """.stripIndent()
}

def helpQC() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                    HaploSync — Assembly QC v2.0                     ║
    ║        Chromosome pair reports and unplaced sequence QC             ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf -entry QC [options]
        nextflow run main.nf -entry QC -params-file params.yml

    QC reads HaploSplit outputs automatically from --outdir/HaploSplit/
    using the --out prefix. Run HaploSplit first, or point to an existing
    HaploSplit output directory.

    The following files are discovered automatically when present:
        {out}.markers.bed           — marker positions in pseudomolecule space
        {out}.legacy_structure.agp  — legacy AGP structure

    Both QC reports run by default. Use skip flags to select which to run.

    ── Required ────────────────────────────────────────────────────────────
        --out               Output prefix (must match HaploSplit run)
                            [default: out]
        --outdir            Results directory (must match HaploSplit run)
                            [default: results]

    ── QC selection ────────────────────────────────────────────────────────
        --skip_chr_pair_reports  Skip chromosome pair overview reports
        --skip_unplaced_qc       Skip unplaced sequence QC reports

    ── Optional inputs ─────────────────────────────────────────────────────
        --markers_map       Marker genetic map (chr, position, marker_id)
        --input_groups      Sequence grouping file
        --legacy_groups     Legacy component group file

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── Profiles ────────────────────────────────────────────────────────────
        -profile standard   Local execution
        -profile conda      Local execution with conda
        -profile mamba      Local execution with mamba/micromamba
        -profile hpc        SLURM execution with mamba/micromamba

    ── Examples ────────────────────────────────────────────────────────────
        # Run both QC reports
        nextflow run main.nf -entry QC -profile mamba \\
            --out myproject --outdir results

        # Run chromosome pair reports only
        nextflow run main.nf -entry QC -profile mamba \\
            --out myproject --outdir results \\
            --skip_unplaced_qc

        # Run unplaced sequence QC only
        nextflow run main.nf -entry QC -profile mamba \\
            --out myproject --outdir results \\
            --skip_chr_pair_reports

        # Using a params file
        nextflow run main.nf -entry QC -profile mamba \\
            -params-file params.yml
    """.stripIndent()
}

// --------------------------------------------------------------------------
// Default entry point: HAPLOSYNC_RECONSTRUCT_PM
//   Pseudomolecule reconstruction + QC + optional HaploDup.
//   Log names: HAPLOSYNC_RECONSTRUCT_PM:HAPLOSPLIT:<PROCESS>
//              HAPLOSYNC_RECONSTRUCT_PM:QC:<PROCESS>
//              HAPLOSYNC_RECONSTRUCT_PM:HAPLODUP:<PROCESS>
// --------------------------------------------------------------------------
workflow {

    if (params.help) {
        helpHaploSplit()
        exit 0
    }

    if (!params.input_fasta) {
        log.error "[ERROR] --input_fasta is required"
        helpHaploSplit()
        exit 1
    }
    def has_markers   = params.markers && params.markers_map
    def has_reference = params.guide_genome
    if (!has_markers && !has_reference) {
        log.error "[ERROR] Provide (--markers AND --markers_map) and/or --guide_genome"
        helpHaploSplit()
        exit 1
    }
    if ((params.markers && !params.markers_map) || (!params.markers && params.markers_map)) {
        log.error "[ERROR] --markers and --markers_map must be provided together"
        helpHaploSplit()
        exit 1
    }

    HAPLOSYNC_RECONSTRUCT_PM()
}

// --------------------------------------------------------------------------
// Entry point: HAPLODUP (standalone)
//   Reads HaploSplit outputs from --outdir/HaploSplit/.
//   Processes called directly for clean log names: HAPLODUP:<PROCESS>
// --------------------------------------------------------------------------
workflow HAPLODUP {

    if (params.help) {
        helpHaploDup()
        exit 0
    }

    def hs_dir = "${params.outdir}/HaploSplit"
    def pfx    = "${hs_dir}/${params.out}"

    def required_files = [
        file("${pfx}.1.fasta"),
        file("${pfx}.Un.fasta"),
        file("${pfx}.correspondence.tsv"),
        file("${pfx}.1.agp"),
        file("${pfx}.Un.agp")
    ]
    required_files.each { f ->
        if (!f.exists()) {
            log.error "[ERROR] Required HaploSplit output not found: ${f}\n         Run HaploSplit first or check --out / --outdir"
            exit 1
        }
    }

    def hap1_fasta     = Channel.fromPath("${pfx}.1.fasta")
    def hap2_fasta     = file("${pfx}.2.fasta").exists()
                             ? Channel.fromPath("${pfx}.2.fasta")
                             : Channel.value([])
    def un_fasta       = Channel.fromPath("${pfx}.Un.fasta")
    def correspondence = Channel.fromPath("${pfx}.correspondence.tsv")

    def agp_ch = Channel.of(
            file("${pfx}.1.agp"),
            file("${pfx}.2.agp"),
            file("${pfx}.Un.agp")
        )
        .filter { it.exists() }
        .collect()

    def markers_bed_f  = file("${pfx}.markers.bed")
    def legacy_agp_f   = file("${pfx}.legacy_structure.agp")
    def annotation_f   = file("${pfx}.annotation.gff3")

    def markers_bed_ch = markers_bed_f.exists() ? Channel.value(markers_bed_f) : Channel.value([])
    def legacy_agp_ch  = legacy_agp_f.exists()  ? Channel.value(legacy_agp_f)  : Channel.value([])
    def annotation_ch  = annotation_f.exists()   ? Channel.value(annotation_f)  : Channel.value([])

    HD_ALIGN(hap1_fasta, hap2_fasta, un_fasta, correspondence)

    def run_gmap     = annotation_f.exists() && !params.No2
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

// --------------------------------------------------------------------------
// Entry point: QC (standalone)
//   Reads HaploSplit outputs from --outdir/HaploSplit/.
//   Processes called directly for clean log names: QC:<PROCESS>
//   --skip_chr_pair_reports / --skip_unplaced_qc control what runs.
// --------------------------------------------------------------------------
workflow QC {

    if (params.help) {
        helpQC()
        exit 0
    }

    def hs_dir = "${params.outdir}/HaploSplit"
    def pfx    = "${hs_dir}/${params.out}"

    def required_files = [
        file("${pfx}.correspondence.tsv"),
        file("${pfx}.1.agp"),
        file("${pfx}.Un.agp")
    ]
    required_files.each { f ->
        if (!f.exists()) {
            log.error "[ERROR] Required HaploSplit output not found: ${f}\n         Run HaploSplit first or check --out / --outdir"
            exit 1
        }
    }

    def correspondence = Channel.fromPath("${pfx}.correspondence.tsv")

    def agp_ch = Channel.of(
            file("${pfx}.1.agp"),
            file("${pfx}.2.agp"),
            file("${pfx}.Un.agp")
        )
        .filter { it.exists() }
        .collect()

    def markers_bed_f  = file("${pfx}.markers.bed")
    def legacy_agp_f   = file("${pfx}.legacy_structure.agp")

    def markers_bed_ch = markers_bed_f.exists() ? Channel.value(markers_bed_f) : Channel.value([])
    def legacy_agp_ch  = legacy_agp_f.exists()  ? Channel.value(legacy_agp_f)  : Channel.value([])

    if (!params.skip_chr_pair_reports && !params.No2) {
        QC_CHR_PAIR(correspondence, agp_ch, markers_bed_ch, legacy_agp_ch)
    }

    if (!params.skip_unplaced_qc && !params.No2) {
        def unused_list_f = file("${pfx}.unused_sequences.list")
        if (!unused_list_f.exists()) {
            log.error "[ERROR] Required file not found for unplaced QC: ${unused_list_f}"
            exit 1
        }

        def hap1_fasta = Channel.fromPath("${pfx}.1.fasta")
        def hap2_fasta = file("${pfx}.2.fasta").exists()
                             ? Channel.fromPath("${pfx}.2.fasta")
                             : Channel.empty()
        def un_fasta   = Channel.fromPath("${pfx}.Un.fasta")

        def fasta_ch = hap1_fasta
            .mix(hap2_fasta)
            .mix(un_fasta)
            .collect()

        QC_REJECTED(
            Channel.value(unused_list_f),
            correspondence,
            fasta_ch,
            agp_ch,
            markers_bed_ch,
            legacy_agp_ch
        )
    }
}

workflow.onComplete {
    log.info (workflow.success
        ? "\n[HaploSync] Pipeline completed successfully.\n  Results: ${params.outdir}"
        : "\n[HaploSync] Pipeline failed. Check logs for details.")
}
