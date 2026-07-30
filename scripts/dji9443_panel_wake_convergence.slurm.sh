#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-panelwake
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }
export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
# Non-interactive submission shells don't put julia on PATH (Manifest pins 1.11.7);
# fall back to the site spack julia-1.11.7 binary (the shared /apps/juliaup launcher is broken).
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"
echo "DJI 9443 direct-potential panel-wake convergence"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"
julia --project=. -t "$THREADS" examples/run_dji9443_panel_wake_convergence.jl
