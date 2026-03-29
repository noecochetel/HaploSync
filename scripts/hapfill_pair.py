#!/usr/bin/env python3
"""
hapfill_pair.py — HaploFill Steps 4-5: haplotype pairing + gap region prep.

Runs HaploFill Steps 4 and 5 only:
  Step 4 — pairwise hap-to-hap minimap2 alignment, PAF uniquification,
            pairing JSON files
  Step 5 — unmatched region extraction, gap mate position mapping,
            flanking sequence extraction (upgradable_regions.json.gz)

Expects ploidy category files from Step 3 to be present in the temp folder.
Outputs gap descriptor files consumed by HF_FILL (Step 6).

Part of the HAPLOSYNC_GAP_FILL Nextflow pipeline.

Usage:
    hapfill_pair.py [HaploFill arguments...]

All arguments are passed directly to HaploFill.py with --resume 4 --stop 5.
"""

import os
import subprocess
import sys


def main():
    haplosync_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    hapfill       = os.path.join(haplosync_dir, 'HaploFill.py')

    args = sys.argv[1:]

    clean_args = []
    skip_next  = False
    for arg in args:
        if skip_next:
            skip_next = False
            continue
        if arg in ('--stop', '--resume'):
            skip_next = True
            continue
        if arg.startswith('--stop=') or arg.startswith('--resume='):
            continue
        clean_args.append(arg)

    cmd = [sys.executable, hapfill] + clean_args + ['--resume', '4', '--stop', '5']
    print('[hapfill_pair] Running: ' + ' '.join(cmd), file=sys.stderr)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == '__main__':
    main()
