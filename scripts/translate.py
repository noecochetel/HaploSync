#!/usr/bin/env python3
"""
translate.py — Coordinate translation for markers, legacy AGP, and annotation.

Reads the AGP files produced by reconstruct.py and translates coordinates from
the original query sequence space to the new pseudomolecule coordinate space.

Outputs (each optional, produced only when the corresponding input is given):
  {out}.markers.bed          — Marker hit positions on pseudomolecules
  {out}.legacy_structure.agp — Input AGP structure mapped to pseudomolecule coords
  {out}.annotation.gff3      — Annotation in pseudomolecule coordinate space
"""

import argparse
import os
import sys
import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

from lib_files.HaploFunct import *
from lib_files.AGP_lib import *
from lib_files.FASTA_lib import *
from lib_files.GFF_lib import *


def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=__doc__)

    # Required
    parser.add_argument('-o', '--out', dest='out', default='out',
        metavar='NAME',
        help='Output files prefix (must match reconstruct -o) [default: out]')

    # Optional translation inputs
    parser.add_argument('-n', '--map', dest='markers_hits',
        metavar='map.bed',
        help='Marker positions on query sequences (BED-like, 4 columns)')
    parser.add_argument('-a', '--input_agp', dest='input_agp',
        metavar='query.agp',
        help='AGP file of query sequence composition (for legacy structure translation)')
    parser.add_argument('--GFF3', dest='gff3',
        metavar='genes.gff3',
        help='Annotation in GFF3 format (for coordinate translation)')
    parser.add_argument('--N2', dest='No2',
        default=False, action='store_true',
        help="Hap2 was not built (skip reading Hap2 AGP)")

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    options = parser.parse_args()

    if not options.markers_hits and not options.input_agp and not options.gff3:
        print('[ERROR] At least one of -n/--map, -a/--input_agp, or --GFF3 must be provided.',
              file=sys.stderr)
        sys.exit(1)

    print("Running translate.py from HaploSync version " + get_version(), file=sys.stdout)
    print("Command: " + " ".join(sys.argv), file=sys.stdout)
    print("----", file=sys.stdout)

    # -------------------------------------------------------------------------
    # Build combined AGP lookup from the three pseudomolecule AGP files
    # -------------------------------------------------------------------------
    print('[' + str(datetime.datetime.now()) + '] = Reading pseudomolecule AGP files', file=sys.stdout)

    agp_db = dict(read_agp(options.out + '.1.agp'))
    if not options.No2:
        agp_db.update(read_agp(options.out + '.2.agp'))
    agp_db.update(read_agp(options.out + '.Un.agp'))

    # -------------------------------------------------------------------------
    # Translate marker BED
    # -------------------------------------------------------------------------
    if options.markers_hits:
        print('[' + str(datetime.datetime.now()) + '] = Translating marker coordinates', file=sys.stdout)
        bed_regions = read_bed_sorted_list(options.markers_hits)
        new_bed = translate_bed_sorted_list(bed_regions, agp_db)
        out_bed = options.out + '.markers.bed'
        with open(out_bed, 'w') as f:
            for line in sorted(new_bed):
                print('\t'.join([str(x) for x in line]), file=f)
        print('[' + str(datetime.datetime.now()) + '] = Written: ' + out_bed, file=sys.stdout)

    # -------------------------------------------------------------------------
    # Translate legacy AGP
    # -------------------------------------------------------------------------
    if options.input_agp:
        print('[' + str(datetime.datetime.now()) + '] = Translating legacy AGP: ' + options.input_agp, file=sys.stdout)
        old_legacy_agp = read_agp(options.input_agp)
        legacy_agp = agp_translate_agp(agp_db, old_legacy_agp)
        legacy_agp_file = options.out + '.legacy_structure.agp'
        write_agp(legacy_agp, legacy_agp_file)
        print('[' + str(datetime.datetime.now()) + '] = Written: ' + legacy_agp_file, file=sys.stdout)

    # -------------------------------------------------------------------------
    # Translate GFF3 annotation
    # -------------------------------------------------------------------------
    if options.gff3:
        print('[' + str(datetime.datetime.now()) + '] = Translating annotation: ' + options.gff3, file=sys.stdout)
        annotation_gff3, mRNA_db = read_gff3(options.gff3)
        translation_db = translate_from_AGP_whole_genome(read_agp(options.out + '.1.agp'))
        if not options.No2:
            translation_db.update(translate_from_AGP_whole_genome(read_agp(options.out + '.2.agp')))
        translation_db.update(translate_from_AGP_whole_genome(read_agp(options.out + '.Un.agp')))

        new_gff3 = translate_gff3(annotation_gff3, translation_db, options.out + '.broken_genes.txt')

        output_seq_lengths = get_fasta_lengths_from_file(options.out + '.1.fasta')
        if not options.No2:
            output_seq_lengths.update(get_fasta_lengths_from_file(options.out + '.2.fasta'))
        output_seq_lengths.update(get_fasta_lengths_from_file(options.out + '.Un.fasta'))

        out_gff3 = options.out + '.annotation.gff3'
        write_gff3(new_gff3, out_gff3, output_seq_lengths)
        print('[' + str(datetime.datetime.now()) + '] = Written: ' + out_gff3, file=sys.stdout)


if __name__ == '__main__':
    main()
