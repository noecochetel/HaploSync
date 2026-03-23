# HaploDup

HaploDup compares the gene content of two haplotypes, delivers interactive graphs, and identifies chimeric scaffolding where alternative haplotypes are included in the same pseudomolecule.

```
usage: HaploDup.py [-h] [-f genome.fasta [Required]] [-c Hap1_to_Hap2.txt [Required]]
                   [-g annotation.gff3] [-a functional_annotation.txt]
                   [-b markers.bed] [--markers_map markers_list]
                   [--agp structure.agp] [--input_groups input_groups.tsv]
                   [--legacy_agp legacy_structure.agp]
                   [--legacy_groups legacy_groups.tsv]
                   [--rejected_list rejected.list] [-o NAME]
                   [--hit_identity N] [--hit_coverage N] [--gene_identity N]
                   [--gene_coverage N] [--unbalanced_ratio N]
                   [-r reference.fasta] [--reuse_mappings] [--reuse_dotplots]
                   [--reuse_gmap] [--skip_chr_pair_reports]
                   [--skip_dotplots_by_chr] [--only_paired_dotplots] [-t N]
                   [--feature [CDS|mRNA]] [-w N] [--allowed N]
```

## Input and Arguments

#### Mandatory

* `-f | --fasta genome.fasta`: FASTA file of pseudomolecule sequences. A comma-separated list of files can be provided if the sequences are split across multiple files. The sequences produced by HaploSync, for example, are split by haplotype.

* `-c | --correspondence Hap1_to_Hap2.txt`: Tab-delimited file(s) relating sequence names across haplotypes and, if present, the reference genome.
  * Columns:
    1. Chromosome ID
    2. Haplotype 1 sequence ID
    3. Haplotype 2 sequence ID
    4. Reference sequence ID *(only if a reference is used)*
  * Ex. `Hap1_to_Hap2.txt`:
    ```text
    Chr01   Hap1_chr01   Hap2_chr01   Ref_chr01
    Chr02   Hap1_chr02   Hap2_chr02   Ref_chr02
    [...]
    ChrNN   Hap1_chrNN   Hap2_chrNN   Ref_chrNN
    ```

#### Optional

* `-g | --gff annotation.gff3`: Gene annotation for the pseudomolecules in GFF3 format. A comma-separated list of annotation files can be provided if the annotation is split across multiple files. If not provided, GMAP-based gene analysis is skipped and only dotplots are produced.

* `-o | --out NAME`: Output file prefix.
  * *default: `out`*

* `-t | --threads N`: Cores used in the mapping process.
  * *default: `4`*

* `-r | --reference reference.fasta`: FASTA file(s) of reference genome sequences (haploid).

* `-a | --annotation functional_annotation.txt`: Text file with the functional annotation of each transcript or gene.

* `-b | --markers_bed markers.bed`: Marker positions on input sequences.

* `--markers_map markers_list`: Table of genomic map markers with sorting order information. The table does not need to be sorted, but should contain information regarding the order of the markers in the map.
  * Tab-delimited file with 3 columns:
    1. Chromosome ID
    2. Marker position
    3. Marker ID

* `--agp structure.agp`: AGP structure of input sequences.

* `--legacy_agp legacy_structure.agp`: AGP structure of input sequences based on legacy sequences.

* `--rejected_list rejected.list`: A list of unplaced sequences associated with a chromosome but rejected from the tiling path.
  * Tab-delimited file with 3 columns:
    1. Unplaced sequence ID
    2. Associated pseudomolecule (either Haplotype 1 or Haplotype 2 sequence ID)
    3. Expected orientation
  * Setting this parameter will run QC of unplaced sequences versus the associated chromosome pseudomolecules.

* `--legacy_groups legacy_groups.tsv`: Tab-delimited file associating the input sequences composing legacy groups. One column reports the component ID, another names the group it belongs to. A sequence can be associated with multiple groups using multiple rows.

* `--input_groups input_groups.tsv`: Tab-delimited file describing groups of input sequences. Used to trace sequence relationships if `--legacy_groups` is not set. One column reports the sequence ID, another its group ID. A sequence can be associated with multiple groups using multiple rows.

* `--hit_identity N`: Minimum genome mapping hit identity.
  * *default: `90`*

* `--hit_coverage N`: Minimum genome mapping hit length.
  * *default: `3000`*

* `--gene_identity N`: Minimum gene mapping identity.
  * *default: `95`*

* `--gene_coverage N`: Minimum gene mapping coverage.
  * *default: `95`*

* `--feature [CDS|mRNA]`: Feature type to use for mapping when a GFF is provided.
  * *default: `CDS`*

* `-w | --window N`: Window size (number of genes) to search for blocks of genes with unbalanced allele counts. Set to 0 to disable.
  * *default: `10`*

* `--allowed N`: Allowed number of unbalanced genes per window.
  * *default: `5`*

* `--unbalanced_ratio N`: Gene count ratio between haplotypes to call the locus unbalanced.
  * Values range: [0-1]
  * *default: `0.33`*

* `--reuse_mappings`: Reuse alignments present in the output folder instead of re-running them.
  * *default: overwrite*

* `--reuse_dotplots`: Reuse dotplots present in the output folder instead of regenerating them.
  * *default: overwrite*

* `--reuse_gmap`: Reuse CDS mappings with GMAP instead of redoing the analysis.
  * *default: overwrite*

* `--skip_chr_pair_reports`: Skip generation of the Hap1 vs Hap2 overview report for each chromosome pair.
  * *default: generate reports*

* `--skip_dotplots_by_chr`: Skip the generation of individual chromosome pair dotplots. Only the whole-genome dotplot is produced.

* `--only_paired_dotplots`: When generating individual chromosome pair dotplots, only produce plots for matched pairs (e.g. chr01 vs chr01). Off-diagonal comparisons (e.g. chr01 vs chr17) are skipped. The whole-genome dotplot is always produced.

HaploDup can be run standalone or invoked automatically from HaploSplit via `--haplodup`. When called from HaploSplit, the correspondence file, combined AGP, and rejected sequence list are prepared automatically.

## Output

* `${NAME}.HaploDup_dir`: Results of HaploDup analysis for the newly produced pseudomolecules. It contains:
  * `index.html`: Open in a web browser to navigate all report files.
  * `index.chr_pair_reports.html`: Index of the chromosome pair overview reports (generated by default, suppressed with `--skip_chr_pair_reports`).
  * See [Pseudomolecules QC report](../Workflows/pseudomolecule_QC.md#2.3.1---Pseudomolecules-QC-report)
  * Other useful intermediate files and folders:
    * `${NAME}.ref.fasta`: Pseudomolecule sequences.
    * `gmap_index/`, `gmap_index.log`, `gmap_index.err`: GMAP index folder and logs.
    * `CDS.on.hap1.gmap.gff3`, `CDS.on.hap2.gmap.gff3`: Alignments of CDS or mRNA sequences on the genome.
    * `diploid_gene_count_trace.hap1.txt`, `diploid_gene_count_trace.hap2.txt`: Tables with hit counts of each gene annotated on the two haplotypes. Columns:
      1. Chromosome ID
      2. Locus start
      3. Locus end
      4. Locus ID
      5. Hap1 counts
      6. Hap2 counts
      7. Locus functional annotation
      8. Number of genes with the given functional description(s) across all annotated genes
    * `${NAME}.hotspots.windows.hap1.txt`, `${NAME}.hotspots.windows.hap2.txt`: Windows on the genome showing hotspots of unbalanced genes between haplotypes. Controlled by `-w | --window N`, `--unbalanced_ratio N`, and `--allowed N`.
    * `diploid_gene_count_trace.hap1.rda`, `diploid_gene_count_trace.hap2.rda`: RDA files used by R to generate the plots from trace counts.

* `${NAME}.structure_comparison`: If `--rejected_list` is given, this folder contains the QC reports for unplaced sequences.
  * `index.rejected_sequences.html`: Open in a web browser to navigate report files.
  * See [Unplaced sequence rejection QC report](../Workflows/pseudomolecule_QC.md#2.3.2---Unplaced-sequence-rejection-QC-report)

## How it works

HaploDup does the following procedure:

* Dotplots compare pseudomolecules and a guide/reference haploid genome.
* CDS sequences are extracted from the input annotation of both haplotypes:
  1. Map the CDS on both haplotypes.
  2. Count the alignments of each CDS on each pseudomolecule for both haplotypes.
  3. Compare the counts of the CDS alignments for each haplotype's pseudomolecules.
  4. Produce an interactive plot for each pseudomolecule to analyze differences in hit counts for each annotated gene.
* Count the markers found in each sequence.
* Relate marker positions in each sequence to the pseudomolecules or unplaced sequences.
* Relate the sequences and legacy structure to their positions in the pseudomolecules or unplaced sequences.

![HaploDup.diagram](HaploDup.diagram.png)

## Plots

### Pseudomolecules QC report

Plots can be viewed in a web browser using `${NAME}.HaploDup_dir/index.html`. Reports include interactive HTML plots. The whole-genome dotplot is also saved as a static PNG. The report will contain:

* An alignment of the query sequence on the target sequence, super-imposed on the alignment of the target pseudomolecule on itself.
* The target sequence structure, with two lines of sequences at most.
* Marker position and duplication status.
* Gene mapping count on both haplotypes.
* Unbalanced gene mapping count between haplotypes.

![Pseudomolecules_QC_report_1](Pseudomolecules_QC_report_1.png)

![Pseudomolecules_QC_report_2](Pseudomolecules_QC_report_2.png)

### Unplaced sequence rejection QC report

Reports can be viewed as interactive HTML plots from `${NAME}.structure_comparison/index.rejected_sequences.html`. Plot content depends on input. Generally, the report contains:

* A comparison of associated pseudomolecule structures in terms of markers and sequence content. The plot can show two levels of structural information (e.g. scaffolds and the contigs that compose them).
* Comparison of the unplaced sequence with the associated pseudomolecules in terms of markers and sequence content. The plot can show two levels of structural information (e.g. scaffolds input to HaploSplit and the legacy contigs that compose them).
* Marker content usage comparison.

![Unplaced_QC](Unplaced_QC.png)
