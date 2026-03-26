#!/usr/bin/env nextflow

/*
 * HaploSync — Haplotype-resolved pseudomolecule assembly pipeline
 *
 * Usage:
 *   nextflow run main.nf -profile conda \
 *       --input_fasta assembly.fasta \
 *       --guide_genome reference.fasta \
 *       --out myproject \
 *       --outdir results
 *
 * For a full parameter list see nextflow.config or run:
 *   nextflow run main.nf --help
 */

nextflow.enable.dsl = 2

include { PSEUDOMOLECULE_GENERATION } from './workflows/pseudomolecule_generation'

// --------------------------------------------------------------------------
// Entry point
// --------------------------------------------------------------------------
workflow {

    // Validate required inputs
    if (!params.input_fasta) {
        error "[ERROR] --input_fasta is required"
    }
    if (!params.guide_genome && !params.markers) {
        error "[ERROR] At least one of --guide_genome or --markers must be provided"
    }

    PSEUDOMOLECULE_GENERATION()
}
