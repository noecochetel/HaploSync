#!/usr/bin/env python
import sys
import os
from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *
from lib_files.map_lib import *

def main():
    if len(sys.argv) == 1 or len(sys.argv) != 6:
        print("Usage: python StructComp.py hap1.agp hap2.agp hap1.fasta hap2.fasta outdir")
        print("Compares AGP structures and generates dotplots and HTML report for Hap1 vs Hap2.")
        sys.exit(1)

    hap1_agp = sys.argv[1]
    hap2_agp = sys.argv[2]
    hap1_fasta = sys.argv[3]
    hap2_fasta = sys.argv[4]
    out_dir = sys.argv[5]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Set paths using config file (same as HaploSplit)
    paths = set_paths(os.path.join(sys.path[0], 'HaploSync.conf.toml'))
    nucmer_path = paths["nucmer"]
    showcoords_path = paths["show-coords"]

    # Read AGP files
    agp1 = read_agp(hap1_agp)
    agp2 = read_agp(hap2_agp)

    # Write structure comparison table
    struct_comp_file = os.path.join(out_dir, "structure_comparison.tsv")
    with open(struct_comp_file, "w") as out:
        out.write("Chr\tComponent_Hap1\tOrientation_Hap1\tStart_Hap1\tEnd_Hap1\tComponent_Hap2\tOrientation_Hap2\tStart_Hap2\tEnd_Hap2\tStatus\n")
        for chr1 in agp1:
            if chr1 in agp2:
                comps1 = agp1[chr1]
                comps2 = agp2[chr1]
                for i in range(max(len(comps1), len(comps2))):
                    c1 = comps1[i] if i < len(comps1) else ("-", "-", "-", "-")
                    c2 = comps2[i] if i < len(comps2) else ("-", "-", "-", "-")
                    status = "MATCH" if c1 == c2 else "DIFF"
                    out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(
                        chr1, c1[0], c1[1], c1[2], c1[3], c2[0], c2[1], c2[2], c2[3], status))
            else:
                for c1 in agp1[chr1]:
                    out.write("{0}\t{1}\t{2}\t{3}\t{4}\t-\t-\t-\t-\tONLY_HAP1\n".format(
                        chr1, c1[0], c1[1], c1[2], c1[3]))
        for chr2 in agp2:
            if chr2 not in agp1:
                for c2 in agp2[chr2]:
                    out.write("{0}\t-\t-\t-\t-\t{1}\t{2}\t{3}\t{4}\tONLY_HAP2\n".format(
                        chr2, c2[0], c2[1], c2[2], c2[3]))

    # Prepare dotplot output prefix
    outfile_prefix = os.path.join(out_dir, "Hap1_vs_Hap2")
    coords_file = outfile_prefix + ".coords"

    # Run nucmer and show-coords (as in HaploSplit)
    coords_file = map_nucmer(hap1_fasta, hap2_fasta, 4, coords_file, nucmer_path, showcoords_path, " --forward ", " -l -r -T -H ")

    # Convert coords to table
    coords_table = make_coords_table(outfile_prefix, out_dir, showcoords_path)

    # Generate dotplot and HTML report
    hap1_ids = ",".join(read_fasta_ids(hap1_fasta))
    hap2_ids = ",".join(read_fasta_ids(hap2_fasta))
    whole_genome_dotplot(hap1_ids, hap2_ids, outfile_prefix, out_dir, "Hap1_vs_Hap2", coords_table)
    make_pair_html_report(os.path.basename(coords_file), "", out_dir, out_dir, "Hap1", "Hap2", "", "", "", "", "Hap1", "Hap2", 0, 0, "", "", "3000", "90", "0.33")

    print(("Structure comparison and dotplots written to {0}".format(out_dir)))

def read_fasta_ids(fasta_file):
    ids = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip())
    return ids

if __name__ == "__main__":
    main()