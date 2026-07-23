#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-veldiag
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=02:30:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }
export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

# Diagnostic: identical DJI axial case as data/dji9443_panelwake/axial_rows072,
# but FORMULATION=velocity (VelocityThroughSources) with the same finite PanelWake.
# Discriminates wake-model deficit (also low) vs potential-pathway defect (recovers).
export RUN_NAME=dji9443_panelwake_diag/axial_velocity_rows072
export RHPC_MESH=40_40 RPM=5400 NT=36 WAKE_ROWS=72
export SETTLE_REVS=3 START_VINF=4.0 TERMINAL_VINF=4.0 FREESTREAM_DECREASE_REVS=0
export FORMULATION=velocity SAVE_VTK=true

echo "DJI 9443 velocity-formulation diagnostic (finite panel wake, axial rows=72)"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"
julia --project=. -t "$THREADS" examples/rotor_panel_wake_study.jl
