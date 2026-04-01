# Nextflow tips

This page covers the Nextflow commands and concepts that are most useful when running HaploSync. It is not a Nextflow tutorial — the [official documentation](https://www.nextflow.io/docs/latest/) covers everything in depth — but it collects the practical knowledge you will need day to day.

---

## Running a workflow

```bash
nextflow run nextflow/reconstruct_pm.nf -profile mamba -params-file params.yml
```

| Element | Purpose |
|---------|---------|
| `nextflow run <script>` | Entry point — path to the `.nf` file |
| `-profile mamba` | Activates the conda/mamba environment automatically |
| `-params-file params.yml` | Loads parameters from a YAML file instead of passing them on the command line |
| `-resume` | Resumes a previous run (see below) |
| `-bg` | Runs Nextflow in the background (see below) |

---

## Running in the background

```bash
nextflow run nextflow/reconstruct_pm.nf -profile mamba -bg -params-file params.yml
```

`-bg` detaches the Nextflow process from your terminal so the run continues even if you close your session. Output is redirected to `.nextflow.log` instead of the terminal.

This is useful on a server or HPC login node where you do not want to keep a terminal open. You can monitor progress at any time:

```bash
# Follow the log live
tail -f .nextflow.log

# Check whether Nextflow is still running
ps aux | grep nextflow
```

> **Tip:** Combine with `nohup` or a terminal multiplexer (`tmux`, `screen`) for extra safety on HPC systems where login sessions may be killed after inactivity.

---

## Resuming an interrupted run

```bash
nextflow run nextflow/reconstruct_pm.nf -profile mamba -resume -params-file params.yml
```

Nextflow caches every successfully completed task in the `work/` directory. With `-resume`, tasks whose inputs have not changed are skipped and their cached outputs reused. Only tasks that are new, failed, or whose inputs changed will re-run.

**Important caveats:**
- Resume only works if the `work/` directory from the previous run is still present.
- Changing a process script (even whitespace) busts its cache and forces re-execution of that step and all downstream steps.
- Changing a parameter that feeds into a cached step also busts that step's cache.

---

## The work directory

Nextflow executes each task in an isolated subdirectory under `work/`. The path is derived from a hash of the task inputs:

```
work/
└── 3f/
    └── a2b8c1d4e5.../
        ├── .command.sh      # the exact shell command that was run
        ├── .command.log     # combined stdout + stderr
        ├── .command.out     # stdout only
        ├── .command.err     # stderr only
        ├── .exitcode        # exit code
        └── <input symlinks and output files>
```

### Finding the work directory for a failed task

When a task fails, Nextflow prints the work directory path in the error message:

```
ERROR ~ Error executing process > 'HAPLOFILL:HF_FILL (chr1)'

Caused by:
  ...

Work dir:
  /path/to/work/3f/a2b8c1d4e5...
```

Navigate there and inspect `.command.log` or `.command.err` to see the full error output:

```bash
cat work/3f/a2b8c1d4e5.../.command.log
```

You can also re-run the exact failing command interactively:

```bash
cd work/3f/a2b8c1d4e5...
bash .command.sh
```

---

## Logs

### Main Nextflow log

`.nextflow.log` is written in the current directory every time you run Nextflow. It contains detailed information about task scheduling, caching decisions, and errors. Previous runs are rotated to `.nextflow.log.1`, `.nextflow.log.2`, etc.

```bash
# Check the last few lines after a failure
tail -50 .nextflow.log

# Search for a specific process name
grep "HF_FILL" .nextflow.log
```

### Execution reports

Nextflow can generate HTML reports summarising resource usage, task durations, and the execution timeline. These are useful for identifying bottlenecks or memory issues:

```bash
nextflow run nextflow/gap_fill.nf -profile mamba \
    -params-file params.yml \
    -with-report report.html \
    -with-timeline timeline.html \
    -with-trace trace.tsv
```

| Flag | Output |
|------|--------|
| `-with-report` | HTML report with CPU, memory, and duration per task |
| `-with-timeline` | Gantt chart of task execution |
| `-with-trace` | Tab-separated table of all task metrics |

---

## Profiles

Profiles are defined in `nextflow/nextflow.config` and control how tasks are executed.

| Profile | Use case |
|---------|---------|
| `mamba` | Local execution, conda environment managed by mamba |
| `conda` | Local execution, conda environment managed by conda |
| `hpc` | SLURM cluster execution |

```bash
# Local run
nextflow run nextflow/reconstruct_pm.nf -profile mamba -params-file params.yml

# HPC run
nextflow run nextflow/reconstruct_pm.nf -profile hpc -params-file params.yml
```

On HPC, resource requests (CPUs, memory, queue) are set per process label in `nextflow/nextflow.config`. Adjust them there if jobs are failing due to resource limits.

---

## Named workflow entry points

Some `.nf` files expose multiple named workflows via `-entry`. This allows running a subset of the pipeline:

```bash
# Run only HaploDup on gap-filled results (reads from {outdir}/HaploMake/)
nextflow run nextflow/gap_fill.nf -profile mamba -entry HAPLODUP -params-file params.yml
```

Standalone entry points for HaploMake and HaploDup are also available as dedicated scripts:

```bash
nextflow run nextflow/haplomake.nf -profile mamba -params-file params_haplomake.yml
nextflow run nextflow/haplodup.nf  -profile mamba -params-file params_haplodup.yml
```

---

## Cleaning up

The `work/` directory accumulates cached task directories and can grow large over time. Use `nextflow clean` to remove old runs while preserving the cache for the most recent one:

```bash
# List runs
nextflow log

# Remove all runs except the last
nextflow clean -but-last -f

# Remove a specific run by session ID or run name
nextflow clean <run-name> -f
```

> **Warning:** Cleaning removes the work directory entries for those runs. Any `-resume` from a cleaned run will re-execute all steps from scratch.
