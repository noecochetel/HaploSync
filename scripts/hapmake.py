#!/usr/bin/env python3
"""
hapmake.py — HaploMake wrapper for the HAPLOSYNC_GAP_FILL Nextflow pipeline.

Constructs new pseudomolecule sequences from a HaploFill structure block file.
Thin wrapper around HaploMake.py — all arguments are forwarded directly.

Outputs (controlled by HaploMake flags):
  {out}.fasta              — new assembled sequences
  {out}.structure.agp      — AGP metadata
  {out}.legacy_structure.agp — legacy coordinate mapping (if --agp provided)
  {out}.bed                — translated feature coordinates (if --bed provided)
  {out}.loci_to_check.txt  — regions needing manual review (if overlaps found)

Part of the HAPLOSYNC_GAP_FILL Nextflow pipeline.

Usage:
    hapmake.py [HaploMake arguments...]

All arguments are passed directly to HaploMake.py.
"""

import os
import subprocess
import sys


def main():
    haplosync_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    hapmake       = os.path.join(haplosync_dir, 'HaploMake.py')

    cmd = [sys.executable, hapmake] + sys.argv[1:]
    print('[hapmake] Running: ' + ' '.join(cmd), file=sys.stderr)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == '__main__':
    main()
