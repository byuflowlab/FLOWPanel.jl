#!/usr/bin/env bash
#SBATCH --job-name=fp-018-p15gate
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# BRAINSTORM/018 Phase-15 instrumentation pre-submission gate
# (ops_reference.md "Pre-submission gate" + phase_15_drift_source.md (b)):
#   stage 1: full test suite (src changed: WakeInventoryMonitor +
#            WakeHealthMonitor attribution columns)
#   stage 2: formulation_test (10 stages)
#   stage 3: NT=6 ~6-rev smoke, RHPC_MESH=40_40 --
#            arm von : velocity, diagnostics ON
#            arm voff: velocity, diagnostics OFF (bit-identity reference)
#            arm gon : green,    diagnostics ON
#   stage 4: bit-identity diff of the physics CSVs (von vs voff) --
#            wake_health/wake_inventory excluded (they differ by design)
#
# Submit from the top level of the FLOWPanel.jl checkout.
set -euo pipefail

# Manifest pins julia 1.11.7 -- same PATH fallback as the production
# launcher (module-default julia is 1.12.6 and fails precompile against
# this Manifest; the shared /apps/juliaup launcher is broken).
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

THREADS=8
export JULIA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export MPLBACKEND=Agg

echo "FLOWPanel phase-15 instrumentation gate"
echo "  repo:    $(pwd)"
echo "  threads: $THREADS"
echo "  julia:   $(julia --version)"

echo "=== stage 1: full runtests ==="
julia --project=. -t "$THREADS" -e 'include("test/runtests.jl")'

echo "=== stage 2: formulation_test ==="
julia --project=. -t "$THREADS" test/formulation_test.jl

echo "=== stage 3: NT=6 smoke arms ==="
export RHPC_MESH=40_40
export NT=6
export RPM=6000
export BERNOULLI_ONLY=true
export SAVE_VTK=true
export SPINUP_REVS=1.0
export SPINUP_START_FRACTION=0.4
export MAGVINF_START=0.0
export MAGVINF_PEAK=5.0
export MAGVINF_END=0.0
export FREESTREAM_RAMP_REVS=0.5
export FREESTREAM_HOLD_REVS=0.5
export FREESTREAM_WITHDRAW_REVS=1.0
export SETTLE_REVS=3
export NREVS=6
export TRUNCATION_DEPTH_R=4

run_arm () {
  local tag="$1" formulation="$2" inv="$3" attr="$4"
  local run_name="p15_smoke_${tag}"
  echo "--- arm ${tag}: formulation=${formulation} WAKE_INVENTORY=${inv} ATTRIBUTION=${attr}"
  rm -rf "data/${run_name}"
  RHPC_FORMULATION="$formulation" WAKE_INVENTORY="$inv" \
    WAKE_HEALTH_ATTRIBUTION="$attr" RUN_NAME="$run_name" \
    julia --project=. -t "$THREADS" examples/rotor_hover_pressure_comparison.jl \
    > "data/p15_smoke_${tag}.log" 2>&1
  echo "arm ${tag} done; final CT lines:"
  grep -E "CF = " "data/p15_smoke_${tag}.log" | tail -2
}

run_arm von  velocity true  true
run_arm voff velocity false false
run_arm gon  green    true  true

echo "=== stage 4: bit-identity diff (physics CSVs, von vs voff) ==="
status=0
for f in data/p15_smoke_voff/monitors/*.csv; do
  b=$(basename "$f")
  case "$b" in
    *wake_health*|*wake_inventory*) continue ;;
  esac
  g="data/p15_smoke_von/monitors/${b/voff/von}"
  if ! cmp -s "$f" "$g"; then
    echo "MISMATCH: $b"
    status=1
  fi
done
for pair in CT_vs_rev CT_per_rev; do
  a="data/p15_smoke_voff/p15_smoke_voff_${pair}.csv"
  b="data/p15_smoke_von/p15_smoke_von_${pair}.csv"
  if [ -f "$a" ] && ! cmp -s "$a" "$b"; then
    echo "MISMATCH: ${pair}"
    status=1
  fi
done
if [ "$status" -eq 0 ]; then
  echo "BIT-IDENTITY PASS: all physics CSVs identical (von vs voff)"
else
  echo "BIT-IDENTITY FAIL"
fi
echo "=== inventory CSV head (arm von) ==="
head -2 data/p15_smoke_von/monitors/*wake_inventory* | cut -c1-400
echo "=== wake_health attribution header (arm von) ==="
head -1 data/p15_smoke_von/monitors/*wake_health*
exit $status
