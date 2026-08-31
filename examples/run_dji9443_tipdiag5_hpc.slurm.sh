#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-tipdiag5
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Phase 2d step F — round-CapUMinTess=4 production-candidate ladder (Ryan prefers round
# tip caps; ct4 is empirically equivalent to ct2/flat and finer than the defective ct3):
# Dirichlet 201/249 + the matched round-ct4-tip Neumann reference at 185.
set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

echo "DJI 9443 Phase 2d step F round-ct4 ladder (201, 249, 185u tip-round)"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"

for TAG in dji45_201c_ct4 dji45_249c_ct4 dji45_185u_tipround4; do
  echo "=== tipdiag $TAG ($(date '+%F %T')) ==="
  if ! env PHASE2C_MODE=tipdiag "PHASE2C_CASE=$TAG" "PHASE2D_LABEL=$TAG" \
      julia --project=. -t "$THREADS" examples/dji9443_mesh_convergence.jl; then
    echo "!!! FAILED: $TAG"
  fi
done

echo "Phase 2d round-ct4 ladder runs finished ($(date '+%F %T'))."
