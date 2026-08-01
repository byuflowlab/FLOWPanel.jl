#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-tipdiag2
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Phase 2d tip-cap discrepancy, round 2 — attached-wake causal perturbations on the
# anomalous 185c (and controls): Das strip-length (eta) sweep and semi-infinite-wake
# variants. If the 185 tip Γ is hypersensitive to the attached-strip geometry while
# 145/249 are not, the wake-attachment mechanism is confirmed.
set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

echo "DJI 9443 Phase 2d tipdiag round 2 (attached-wake perturbations)"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"

# run2 CASE LABEL EXTRA_ENV...
run2() {
  local case="$1" label="$2"; shift 2
  echo "=== tipdiag $label ($(date '+%F %T')) ==="
  if ! env PHASE2C_MODE=tipdiag "PHASE2C_CASE=$case" "PHASE2D_LABEL=$label" "$@" \
      julia --project=. -t "$THREADS" examples/dji9443_mesh_convergence.jl; then
    echo "!!! FAILED: $label"
  fi
}

run2 dji45_185c dji45_185c_eta0.1   PHASE2D_DAS_ETA=0.1
run2 dji45_185c dji45_185c_eta0.4   PHASE2D_DAS_ETA=0.4
run2 dji45_185c dji45_185c_eta1.0   PHASE2D_DAS_ETA=1.0
run2 dji45_185c dji45_185c_semiinf  PHASE2D_SEMIINF=1
run2 dji45_249c dji45_249c_eta0.4   PHASE2D_DAS_ETA=0.4
run2 dji45_249c dji45_249c_semiinf  PHASE2D_SEMIINF=1
run2 dji45_201c dji45_201c_semiinf  PHASE2D_SEMIINF=1
run2 dji45_185u dji45_185u_semiinf  PHASE2D_SEMIINF=1

echo "All Phase 2d round-2 runs finished ($(date '+%F %T'))."
