#!/usr/bin/env bash
#SBATCH --job-name=fp-pitch-pressure
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=36g
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Submit from the top level of the FLOWPanel.jl checkout with:
#   sbatch examples/run_pitching_wing_pressure_comparison.slurm.sh
#
# Override run settings without editing this file:
#   N_CYCLES=3 C_PER_DT=0.5 SAVE_VTK=true sbatch ...

set -euo pipefail

THREADS=48
RUN_DIR="${PRESSURE_COMPARISON_PATH:-data/pitching_wing_pressure_comparison}"

# Julia's thread count is fixed at process startup. Keep it and the native
# OpenMP pool aligned with the Slurm CPU allocation.
export JULIA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

# Preserve the example's useful diagnostic defaults while allowing lightweight
# batch runs when requested.
export PRESSURE_COMPARISON_PATH="$RUN_DIR"
export C_PER_DT="${C_PER_DT:-0.5}"
export N_CYCLES="${N_CYCLES:-3.0}"
export SAVE_VTK="${SAVE_VTK:-true}"
export PRESSURE_COMPARISON_PLOT="${PRESSURE_COMPARISON_PLOT:-true}"
export FMM_EXPANSION_ORDER="${FMM_EXPANSION_ORDER:-8}"
export FMM_ACCEPTANCE="${FMM_ACCEPTANCE:-0.4}"
export FMM_LEAF_SIZE="${FMM_LEAF_SIZE:-40}"
export PRESSURE_ITMAX_PER_PANEL="${PRESSURE_ITMAX_PER_PANEL:-2.0}"

echo "FLOWPanel pitching-wing pressure comparison"
echo "  repo:      $(pwd)"
echo "  output:    $PRESSURE_COMPARISON_PATH"
echo "  CPUs:      $THREADS"
echo "  Julia:     $JULIA_NUM_THREADS threads"
echo "  OpenMP:    $OMP_NUM_THREADS threads"
echo "  BLAS:      $BLAS_NUM_THREADS threads"
echo "  itmax/panel: $PRESSURE_ITMAX_PER_PANEL"
echo "  C_PER_DT:  $C_PER_DT"
echo "  N_CYCLES:  $N_CYCLES"
echo "  SAVE_VTK:  $SAVE_VTK"

julia --project=. -t "$THREADS" examples/pitching_wing_pressure_comparison.jl
