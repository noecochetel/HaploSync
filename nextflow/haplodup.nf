#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ALIGN  as HD_ALIGN  } from './modules/local/haplodup_align/main'
include { GMAP   as HD_GMAP   } from './modules/local/haplodup_gmap/main'
include { REPORT as HD_REPORT } from './modules/local/haplodup_report/main'

// --------------------------------------------------------------------------
// Help message
// --------------------------------------------------------------------------
def helpHaploDup() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                   HaploSync — HaploDup v2.0                         ║
    ║              Haplotype duplication QC and reporting                 ║
    ╚══════════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run nextflow/haplodup.nf [options]
        nextflow run nextflow/haplodup.nf -params-file params.yml

    Runs all-vs-all nucmer alignments between haplotypes, optionally maps
    gene models with GMAP, and generates HTML/PDF reports.

    ── Required ────────────────────────────────────────────────────────────
        --hap1_fasta        Hap1 pseudomolecule FASTA
        --hap2_fasta        Hap2 pseudomolecule FASTA
        --correspondence    Chromosome correspondence TSV

    ── Optional inputs ─────────────────────────────────────────────────────
        --unplaced_fasta    Unplaced sequences FASTA
        --agp               AGP file(s), comma-separated
        --gff3              Gene annotation GFF3 (enables GMAP mapping)
        --markers_bed       Marker positions BED
        --legacy_agp        Legacy AGP structure
        --reference         Reference genome FASTA (for additional dotplots)
        --markers_map       Marker genetic map
        --functional_annotation  Functional annotation per transcript
        --input_groups      Sequence grouping file
        --legacy_groups     Legacy component group file

    ── Alignment thresholds ─────────────────────────────────────────────────
        --hit_identity      Min genome alignment identity  [default: 90]
        --hit_coverage      Min genome alignment length    [default: 3000]
        --gene_identity     Min gene mapping identity      [default: 95]
        --gene_coverage     Min gene mapping coverage      [default: 95]
        --unbalanced_ratio  Gene count imbalance threshold [default: 0.33]

    ── Gene mapping ─────────────────────────────────────────────────────────
        --haplodup_feature  GFF feature type [CDS|mRNA]   [default: CDS]
        --haplodup_window   Window size for unbalanced gene search [default: 10]
        --haplodup_allowed  Allowed unbalanced genes/window [default: 5]

    ── Report options ───────────────────────────────────────────────────────
        --reuse_dotplots        Reuse existing dotplots    [default: false]
        --skip_dotplots_by_chr  Skip per-chromosome dotplots [default: false]
        --only_paired_dotplots  Only generate matched-pair dotplots [default: false]

    ── Output ──────────────────────────────────────────────────────────────
        --out               Output files prefix          [default: out]
        --outdir            Results directory            [default: results]

    ── Resources ───────────────────────────────────────────────────────────
        --cores             CPU cores per process        [default: 4]

    ── Profiles ────────────────────────────────────────────────────────────
        -profile standard   Local execution
        -profile conda      Local execution with conda
        -profile mamba      Local execution with mamba/micromamba
        -profile hpc        SLURM execution with mamba/micromamba

    ── Examples ────────────────────────────────────────────────────────────
        # Basic run — dotplots and chromosome board only
        nextflow run nextflow/haplodup.nf -profile mamba \\
            --hap1_fasta hap1.fasta --hap2_fasta hap2.fasta \\
            --correspondence correspondence.tsv \\
            --out myproject --outdir results

        # With gene annotation and reference genome
        nextflow run nextflow/haplodup.nf -profile mamba \\
            --hap1_fasta hap1.fasta --hap2_fasta hap2.fasta \\
            --correspondence correspondence.tsv \\
            --gff3 annotation.gff3 \\
            --reference reference.fasta \\
            --out myproject --outdir results

        # Using a params file
        nextflow run nextflow/haplodup.nf -profile mamba \\
            -params-file nextflow/params_haplodup.yml

        # HPC (SLURM)
        nextflow run nextflow/haplodup.nf -profile hpc \\
            -params-file nextflow/params_haplodup.yml
    """.stripIndent()
}

// --------------------------------------------------------------------------
// Entry point: HAPLODUP
// --------------------------------------------------------------------------
workflow {

    if (params.help) {
        helpHaploDup()
        exit 0
    }

    def required = [
        hap1_fasta:      '--hap1_fasta',
        hap2_fasta:      '--hap2_fasta',
        correspondence:  '--correspondence'
    ]
    required.each { param, flag ->
        if (!params[param]) {
            log.error "[ERROR] ${flag} is required"
            helpHaploDup()
            exit 1
        }
    }

    def hap1_fasta     = Channel.fromPath(params.hap1_fasta)
    def hap2_fasta     = Channel.fromPath(params.hap2_fasta)
    def un_fasta       = params.unplaced_fasta
                             ? Channel.fromPath(params.unplaced_fasta)
                             : Channel.value([])
    def correspondence = Channel.fromPath(params.correspondence)

    def agp_ch = params.agp
        ? Channel.fromPath(params.agp.tokenize(',')).collect()
        : Channel.value([])

    def markers_bed_ch = params.markers_bed
        ? Channel.fromPath(params.markers_bed)
        : Channel.value([])
    def legacy_agp_ch  = params.legacy_agp
        ? Channel.fromPath(params.legacy_agp)
        : Channel.value([])
    def annotation_ch  = params.gff3
        ? Channel.fromPath(params.gff3)
        : Channel.value([])

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

workflow.onComplete {
    log.info (workflow.success
        ? "\n[HaploSync] HaploDup completed successfully.\n  Results: ${params.outdir}"
        : "\n[HaploSync] HaploDup failed. Check logs for details.")
}
