#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-tipdiag4
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Phase 2d step F — flat-cap Dirichlet ladder extension (n_airfoil 201, 249) to
# demonstrate the converged flat-cap configuration across the full chordwise ladder.
set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

echo "DJI 9443 Phase 2d step F flat-cap ladder (201, 249)"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"

for TAG in dji45_201c_flat dji45_249c_flat; do
  echo "=== tipdiag $TAG ($(date '+%F %T')) ==="
  if ! env PHASE2C_MODE=tipdiag "PHASE2C_CASE=$TAG" "PHASE2D_LABEL=$TAG" \
      julia --project=. -t "$THREADS" examples/dji9443_mesh_convergence.jl; then
    echo "!!! FAILED: $TAG"
  fi
done

echo "Phase 2d flat-ladder runs finished ($(date '+%F %T'))."
