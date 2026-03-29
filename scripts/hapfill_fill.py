#!/usr/bin/env python3
"""
hapfill_fill.py — HaploFill Step 6: gap patching.

Runs HaploFill Step 6 only:
  Step 6.1 — gather gap information from descriptor files
  Step 6.2 — align unplaced sequences on gap targets (minimap2 per gap)
  Step 6.3 — select best filler per gap, write output

Expects gap descriptor files from Steps 4-5 to be present in the temp folder.

Outputs:
  {prefix}.structure.block          — gap fill instructions for HaploMake
  {prefix}.gap_filling_findings.txt — detailed per-gap report

Part of the HAPLOSYNC_GAP_FILL Nextflow pipeline.

Usage:
    hapfill_fill.py [HaploFill arguments...]

All arguments are passed directly to HaploFill.py with --resume 6.
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

    cmd = [sys.executable, hapfill] + clean_args + ['--resume', '6']
    print('[hapfill_fill] Running: ' + ' '.join(cmd), file=sys.stderr)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == '__main__':
    main()
