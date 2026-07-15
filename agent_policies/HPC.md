# HPC and Slurm Workflows

This file is the source of truth for preparing, reviewing, and running long
simulations or Slurm workflows in FLOWPanel.jl.

## When and Where to Run

Login nodes are for editing, inspection, and light checks only. Jobs expected
to take more than roughly 20–30 laptop minutes belong on HPC. Request only the
resources and hardware features the case needs, and configure account, QOS,
constraints, and modules for the user's allocation instead of assuming
project-specific values.

Repository Slurm scripts are submitted from the top level of the
`FLOWPanel.jl/` checkout. Keep project, example, data, and output paths relative
to that working directory, for example `--project=.`, `examples/case.jl`, and
`data/case/`. Do not derive or export a `REPO_ROOT` from `BASH_SOURCE`, and do
not `cd` inside a batch script unless the user explicitly requests a launcher
that works from arbitrary directories.

## Standard Single-Node Julia Launcher

Use the following pattern for an ordinary single-process, threaded Julia job:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=flowpanel-case
#SBATCH --nodes=1
#SBATCH --ntasks=<N>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<HH:MM:SS>
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Submit from the top level of the FLOWPanel.jl checkout.
set -euo pipefail

THREADS=<N>

export JULIA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

echo "FLOWPanel case"
echo "  repo:    $(pwd)"
echo "  threads: $THREADS"

julia --project=. -t "$THREADS" examples/case.jl
```

Replace both occurrences of `<N>` with the CPU count appropriate for the case;
it may be lower or higher than the counts used by existing launchers. The
`--ntasks=<N>` convention is standard, but no particular value of `<N>` is a
repository-wide default. Slurm launchers must declare explicit memory, time,
output-log, and error-log directives. Use `set -euo pipefail` so a failed
workflow stage stops the job. Slurm opens log paths before the script runs, so
any requested log directory must already exist when the user submits the job.

Set the single `THREADS=<N>` variable explicitly to the same CPU count requested
by `#SBATCH --ntasks=<N>`, and export it consistently through
`JULIA_NUM_THREADS`, `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
`BLAS_NUM_THREADS`, and `MKL_NUM_THREADS`. Pass the same value to Julia with
`-t "$THREADS"` and print it near the start of the log for auditability.

The batch script already executes inside its allocation, so invoke Julia
directly. Avoid `srun` for ordinary single-process threaded Julia jobs. It is
appropriate only for a genuine multi-process launch step or when required by
the site, in which case its task and CPU request must match the allocation.

## Inputs and Checkpoints

Before submission, confirm that every input asset used by the script is present
on the cluster checkout. An untracked local CSV or mesh will not appear in a
fresh clone; copy it to the matching repository path or pass its cluster path
explicitly through the script's environment. Do not silently substitute a
different geometry or polar input merely to make a job start.

For multi-stage jobs, retain completed, validated stage artifacts in the run
directory and let the driver checkpoint-skip that stage on resubmission. The
rotor axial comparison treats both tagged CCBlade sectional CSVs, its polar
CSV, and its validation report as the complete XFOIL/CCBlade checkpoint. Set
`FORCE_CCBLADE=1` only when those polars must be deliberately rebuilt.

## Simulation Outputs

When launching a simulation from an example, especially a diagnostic or
reproduction run for an active investigation, leave VTK output enabled by
default. The I/O cost is usually small relative to the value of retaining
ParaView-ready state. Do not set `SAVE_VTK=false` unless the user asks for a
no-output run.

To avoid filling the disk across repeated iterations, write each new run over
the previous run's directory instead of creating a directory per attempt. Use
one of these patterns:

- Reuse the example's default `save_path`, typically `data/<run_name>/`, and
  let `simulate!` overwrite per-step files in place. Before launching, remove
  the previous run directory so stale steps past the current run length do not
  linger.
- When side-by-side comparison is needed, use one persistent sibling directory
  per scenario, such as `data/<run_name>_nocouple/`, and overwrite it on each
  rerun. Do not suffix directories with timestamps or attempt numbers.

When the user explicitly wants to preserve a previous run for comparison, ask
before overwriting and offer to move the old directory aside.

## Submission Boundary

Agents may prepare and syntax-check Slurm scripts, but must never run `sbatch`,
`srun`, `salloc`, or otherwise launch a supercomputer job. Submission is left
to the user. See BYU's [Slurm guidance](https://rc.byu.edu/wiki/?id=Slurm) and
[script generator](https://rc.byu.edu/documentation/slurm/script-generator).
