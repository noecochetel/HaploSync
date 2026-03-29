#!/usr/bin/env python3
"""
hapfill_ploidy.py — HaploFill Step 3: local ploidy classification.

Runs HaploFill Step 3 only (median coverage calculation + per-position
ploidy categorisation). Expects coverage signal files (.cov.txt.gz) to
already be present in the temp folder — produced by HF_COVERAGE jobs.

Part of the HAPLOSYNC_GAP_FILL Nextflow pipeline.

Usage:
    hapfill_ploidy.py [HaploFill arguments...]

All arguments are passed directly to HaploFill.py with --resume 3 --stop 3.
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

    cmd = [sys.executable, hapfill] + clean_args + ['--resume', '3', '--stop', '3']
    print('[hapfill_ploidy] Running: ' + ' '.join(cmd), file=sys.stderr)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == '__main__':
    main()
