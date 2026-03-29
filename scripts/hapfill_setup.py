#!/usr/bin/env python3
"""
hapfill_setup.py — HaploFill Step 1: setup.

Runs HaploFill Step 1 only (sequence splitting, pairing, repeat parsing,
gap detection). Outputs the temp folder structure used by all downstream
HaploFill modules (HF_COVERAGE, HF_PLOIDY, HF_PAIR, HF_FILL).

Part of the HAPLOSYNC_GAP_FILL Nextflow pipeline.

Usage:
    hapfill_setup.py [HaploFill arguments...]

All arguments are passed directly to HaploFill.py with --stop 1 appended.
"""

import os
import subprocess
import sys


def main():
    haplosync_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    hapfill       = os.path.join(haplosync_dir, 'HaploFill.py')

    args = sys.argv[1:]

    # Remove any pre-existing --stop or --resume flags (we control step range)
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

    cmd = [sys.executable, hapfill] + clean_args + ['--stop', '1']
    print('[hapfill_setup] Running: ' + ' '.join(cmd), file=sys.stderr)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == '__main__':
    main()
