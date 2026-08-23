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

If a run should use Julia or Python, include the appropriate module with `module load julia python` etc. Note that Julia 1.12 is used on the HPC by default, but 1.11 can be requested during the `module load julia` command by specifying the julia version.

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

To avoid filling the disk across repeated iterations, care must be taken to
delete the paraview files (*.vtu, *.vts, etc.) from prior runs. Unless a simulation
is identified as obselete, keep the concise .csv logs. Don't let FLOWPanel take
up more than 200G at any time. This may require deleting paraview files while a
simulation is still live. Don't ever delete paraview files for the newest 36 timesteps (raised from 10 by Ryan, 2026-08-20)
unless a simulation is obselete.

When the user explicitly wants to preserve a previous run for comparison, ask
before overwriting and offer to move the old directory aside.

### Enforcing the 200 G cap

The policy above is enforced mechanically by `scripts/p018_vtk_sweeper.sh`
(campaign-generic despite the name; `PROTECT_FILE=` and `DATA_DIR=` are
environment-overridable). **Do not hand-roll VTK deletions.** A step is only
restartable if all four paths the warmstart loader reads exist, and `.vtm`
indices must be kept together with their pieces — a sweep that left index stubs
without pieces killed job 13036477 on 2026-08-04. The sweeper encodes that rule;
ad-hoc `rm` does not.

**Launch the `hpc-storage` subagent** rather than doing this inline whenever:

- a job that writes VTK is submitted, or one is running unattended;
- checkout usage crosses roughly 100 G — well below the cap, because with ~20
  concurrent writers the checkout grows ~24 GB/h and sweeping late is a race;
- a new sweep arm is about to start; or
- the cluster reports any disk-space or quota error.

Pass it the run names in flight if you know them. It performs one sweep cycle
and returns a before/after report plus a ledger line for the active campaign's
`ledger.md`; it has no write tools, so the main session appends that line.

For continuous coverage during a campaign, the *main session* installs a
`Monitor` that probes `du -sm /home/rander39/projects/FLOWPanel.jl` every ~10
minutes and wakes on ≥100 G, then launches `hpc-storage` for the cycle. Harden
the probe: require the specific line ending in `projects/FLOWPanel.jl` and treat
its absence as a probe failure rather than a small number, `tr -cd '0-9'` before
any arithmetic (the login banner emits ANSI codes on the first stdout line), and
allow two consecutive failures before alarming.

The VTK protect list
(`BRAINSTORM/018_dji9443_hover_convergence_campaign/vtk_protect_list.txt`, with
a matching cluster copy that is the one that takes effect) is Ryan's file:
agents read it and never write it. If a run needs protecting, ask him.

For BYU agent policy, read `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` from the login node.

See BYU's [Slurm guidance](https://rc.byu.edu/wiki/?id=Slurm) and
[script generator](https://rc.byu.edu/documentation/slurm/script-generator) if you are having trouble.
