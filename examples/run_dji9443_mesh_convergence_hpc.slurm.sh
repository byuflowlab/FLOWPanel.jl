#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-mesh-ext
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=02:30:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Phase 2c extended chordwise ladder — the 45-span n_airfoil 185/201/249 solves deferred
# from the laptop (dense Backslash needs 10/12/18 GB). Submit from the top level of the
# cluster FLOWPanel.jl checkout. Each case runs in a FRESH Julia process so each large
# dense solve gets clean memory (peak ~18 GB for 249 capped).
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

echo "DJI 9443 Phase 2c extended chordwise ladder (185/201/249)"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"

for TAG in dji45_185c dji45_185u dji45_201c dji45_201u dji45_249c dji45_249u; do
  echo "=== solving $TAG ($(date '+%F %T')) ==="
  PHASE2C_MODE=solve PHASE2C_CASE="$TAG" \
    julia --project=. -t "$THREADS" examples/dji9443_mesh_convergence.jl
done

echo "All extended-ladder solves finished ($(date '+%F %T'))."
echo "Raw CSVs under data/dji_convergence_20260722/phase_02c_dji_mesh_convergence/raw/"
