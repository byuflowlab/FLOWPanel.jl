#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-mesh-249
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Phase 2c extended ladder — re-run only the finest n_airfoil=249 pair after the initial
# job (12890990) completed 185/201 but failed on 249c in grad_mu post-processing
# (basis=:quad degenerate agglomerate). The driver now uses the robust grad_mu path
# (basis=:tri, tri_robust=true), which is result-neutral for the TE circulation
# (verified bit-identical on dji81c). Each case runs in a fresh Julia process (~18 GB).
set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

echo "DJI 9443 Phase 2c extended ladder — n_airfoil=249 re-run"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"

for TAG in dji45_249c dji45_249u; do
  echo "=== solving $TAG ($(date '+%F %T')) ==="
  PHASE2C_MODE=solve PHASE2C_CASE="$TAG" \
    julia --project=. -t "$THREADS" examples/dji9443_mesh_convergence.jl
done

echo "n_airfoil=249 solves finished ($(date '+%F %T'))."
