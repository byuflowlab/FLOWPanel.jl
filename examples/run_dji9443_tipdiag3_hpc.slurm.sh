#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-tipdiag3
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Phase 2d experiment E2 — cap-treatment control solves at n_airfoil=185 (flat tip cap;
# CapUMinTess 2 and 4 round-cap variants). If a cap-treatment change snaps the 185 tip
# circulation back onto the {145,249} family, the tip-cap tessellation is causally
# confirmed and the fix is a mesh-recipe prescription.
set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

echo "DJI 9443 Phase 2d tipdiag round 3 (E2 cap controls)"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"

run3() {
  local case="$1" label="$2"
  echo "=== tipdiag $label ($(date '+%F %T')) ==="
  if ! env PHASE2C_MODE=tipdiag "PHASE2C_CASE=$case" "PHASE2D_LABEL=$label" \
      julia --project=. -t "$THREADS" examples/dji9443_mesh_convergence.jl; then
    echo "!!! FAILED: $label"
  fi
}

run3 dji45_185c_flat dji45_185c_flat
run3 dji45_185c_ct2  dji45_185c_ct2
run3 dji45_185c_ct4  dji45_185c_ct4

echo "All Phase 2d round-3 runs finished ($(date '+%F %T'))."
