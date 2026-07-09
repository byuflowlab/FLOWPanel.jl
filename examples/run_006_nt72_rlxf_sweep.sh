#!/usr/bin/env bash
#
# Item 006 follow-up: NT=72 rotor-hover stable-wake runs with the R/2 downstream
# relaxation filter ON, sweeping the corrected-Pedrizzetti relaxation factor
# rlxf = 0.3 and 0.6. Intended to run on a separate machine.
#
# Two runs, sequential (each uses THREADS threads; avoid oversubscription):
#   A: RELAX_RLXF=0.3  -> data/rotor_hover_relax006_nt72_filt_rlxf03
#   B: RELAX_RLXF=0.6  -> data/rotor_hover_relax006_nt72_filt_rlxf06
#
# Usage (from the repo root):
#   bash examples/run_006_nt72_rlxf_sweep.sh
# Override thread count:
#   THREADS=4 bash examples/run_006_nt72_rlxf_sweep.sh
#
set -u

# Move to repo root (parent of this script's examples/ dir).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
cd "$SCRIPT_DIR/.." || exit 1
echo "Repo root: $(pwd)"

THREADS="${THREADS:-2}"

# Pin BLAS/OpenMP to the same thread budget as Julia to avoid oversubscription.
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"

# --- Shared schedule (validated 005 non-damping stable-wake baseline) -------
# SFS + viscous are ON by default in the example.
export RHPC_MESH=40_40
export NT=72                       # <-- doubled temporal resolution vs item 006
export RPM=6000
export BERNOULLI_ONLY=true
export SAVE_VTK=true
export SPINUP_REVS=1.5
export SPINUP_START_FRACTION=0.4
export MAGVINF_START=0.0
export MAGVINF_PEAK=5.0
export MAGVINF_END=0.0
export FREESTREAM_RAMP_REVS=1.0
export FREESTREAM_HOLD_REVS=1.5
export FREESTREAM_WITHDRAW_REVS=12
export SETTLE_REVS=12
export TRUNCATION_DEPTH_R=4
export NREVS=10                    # below schedule total; required_revs drives length (~28 rev)

# --- Item 006 R/2 downstream relaxation filter, ON for both runs ------------
export RELAX_FILTER_DOWNSTREAM_R=0.125

run_case () {
  local rlxf="$1"
  local tag="$2"
  local run_name="rotor_hover_relax006_nt72_filt${RELAX_FILTER_DOWNSTREAM_R}_${tag}"
  echo "==============================================================="
  echo "Run: $run_name  (RELAX_RLXF=$rlxf, RELAX_FILTER_DOWNSTREAM_R=$RELAX_FILTER_DOWNSTREAM_R, THREADS=$THREADS)"
  echo "Start: $(date)"
  echo "==============================================================="
  rm -rf "data/${run_name}"
  RELAX_RLXF="$rlxf" RUN_NAME="$run_name" \
    julia --project --threads "$THREADS" examples/rotor_hover_pressure_comparison.jl \
    > "data/006_nt72_${tag}.log" 2>&1
  local rc=$?
  echo "Done: $run_name  exit=$rc  $(date)"
  return $rc
}

# Run A then Run B, sequentially.
run_case 0.3 rlxf03
run_case 0.6 rlxf06

echo "BOTH_RUNS_DONE"
echo "CT history:"
echo "  data/rotor_hover_relax006_nt72_filt_rlxf03/rotor_hover_relax006_nt72_filt_rlxf03_CT_vs_rev.csv  (col CT_bernoulli)"
echo "  data/rotor_hover_relax006_nt72_filt_rlxf06/rotor_hover_relax006_nt72_filt_rlxf06_CT_vs_rev.csv  (col CT_bernoulli)"
echo "Plateau diagnostics: re-derive with"
echo "  julia --project examples/analyze_stable_wake_oscillation.jl data/<RUN_NAME> 12 6000"
