#!/usr/bin/env bash
#SBATCH --job-name=fp-dji9443-tipdiag
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Phase 2d tip-cap discrepancy — instrumented solves (experiment A) and the
# kerneloffset/kernelcutoff regularization sweep (experiment C) on the anomalous
# 45-span n_airfoil=185 capped mesh, plus conditioning references on 185u/201c/249c.
# Submit from the top level of the cluster FLOWPanel.jl checkout. Each run is a FRESH
# Julia process so each large dense solve gets clean memory (185c base keeps an
# unfactored copy of G for an exact residual: ~20 GB matrices).
set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

echo "DJI 9443 Phase 2d tipdiag + regularization sweep"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"

# run_tipdiag CASE LABEL KOFF_PANEL KCUTOFF KEEP_A
# Empty KOFF_PANEL/KCUTOFF means "study default". A failed point is reported but does
# not abort the remaining sweep points.
run_tipdiag() {
  local case="$1" label="$2" kp="$3" kc="$4" keep="$5"
  echo "=== tipdiag $label ($(date '+%F %T')) ==="
  local env=(PHASE2C_MODE=tipdiag "PHASE2C_CASE=$case" "PHASE2D_LABEL=$label" "PHASE2D_KEEP_A=$keep")
  [[ -n "$kp" ]] && env+=("PHASE2D_KOFF_PANEL=$kp")
  [[ -n "$kc" ]] && env+=("PHASE2D_KCUTOFF=$kc")
  if ! env "${env[@]}" julia --project=. -t "$THREADS" examples/dji9443_mesh_convergence.jl; then
    echo "!!! FAILED: $label"
  fi
}

# --- Experiment A: instrumented base solves (study-default regularization) ---
run_tipdiag dji45_185c dji45_185c_base "" "" 1   # exact residual (keeps copy of G)
run_tipdiag dji45_185u dji45_185u_base "" "" 0   # Neumann reference instrumentation
run_tipdiag dji45_201c dji45_201c_base "" "" 0   # second anomalous mesh
run_tipdiag dji45_249c dji45_249c_base "" "" 0   # trustworthy-group fine mesh

# --- Experiment C: regularization sweep on 185c (R = 0.119; base kp = R*1e-10) ---
run_tipdiag dji45_185c dji45_185c_kp0    0.0      "" 0
run_tipdiag dji45_185c dji45_185c_kp1e-7 1.19e-7  "" 0   # R*1e-6
run_tipdiag dji45_185c dji45_185c_kp1e-5 1.19e-5  "" 0   # R*1e-4
run_tipdiag dji45_185c dji45_185c_kp1e-4 1.19e-4  "" 0   # R*1e-3 ~ tip panel scale
run_tipdiag dji45_185c dji45_185c_kc0    ""       0.0     0
run_tipdiag dji45_185c dji45_185c_kc1e-7 ""       1.19e-7 0

echo "All Phase 2d tipdiag runs finished ($(date '+%F %T'))."
echo "CSVs under data/dji_convergence_20260722/phase_02d_tip_cap_discrepancy/raw/"
