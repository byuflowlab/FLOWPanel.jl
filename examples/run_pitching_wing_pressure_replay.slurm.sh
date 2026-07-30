#!/usr/bin/env bash
#SBATCH --job-name=fp-pitch-replay
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=36g
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Submit from the top level of the synced FLOWPanel.jl checkout. The default
# replays all 989 saved samples. A contiguous smoke run can be requested with:
#   PRESSURE_REPLAY_STEPS=0:2 PRESSURE_REPLAY_PLOT=false \
#       sbatch examples/run_pitching_wing_pressure_replay.slurm.sh

set -euo pipefail

THREADS=48

export JULIA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

export PRESSURE_REPLAY_INPUT_PATH="${PRESSURE_REPLAY_INPUT_PATH:-data/pitching_wing_convergence/cases/wake_2__c_per_dt_0p25}"
export PRESSURE_REPLAY_OUTPUT_PATH="${PRESSURE_REPLAY_OUTPUT_PATH:-data/pitching_wing_pressure_replay}"
export PRESSURE_REPLAY_RUN_NAME="${PRESSURE_REPLAY_RUN_NAME:-pitching_wing}"
export PRESSURE_REPLAY_STEPS="${PRESSURE_REPLAY_STEPS:-0:988}"
export PRESSURE_REPLAY_PLOT="${PRESSURE_REPLAY_PLOT:-true}"
export PRESSURE_REPLAY_VERBOSE="${PRESSURE_REPLAY_VERBOSE:-true}"
export FMM_EXPANSION_ORDER="${FMM_EXPANSION_ORDER:-8}"
export FMM_ACCEPTANCE="${FMM_ACCEPTANCE:-0.4}"
export FMM_LEAF_SIZE="${FMM_LEAF_SIZE:-40}"
export PRESSURE_ITMAX_PER_PANEL="${PRESSURE_ITMAX_PER_PANEL:-2.0}"

echo "FLOWPanel pitching-wing pressure replay"
echo "  repo:        $(pwd)"
echo "  input:       $PRESSURE_REPLAY_INPUT_PATH"
echo "  output:      $PRESSURE_REPLAY_OUTPUT_PATH"
echo "  steps:       $PRESSURE_REPLAY_STEPS"
echo "  threads:     $THREADS"
echo "  plot:        $PRESSURE_REPLAY_PLOT"
echo "  FMM order:   $FMM_EXPANSION_ORDER"
echo "  FMM accept:  $FMM_ACCEPTANCE"
echo "  FMM leaf:    $FMM_LEAF_SIZE"
echo "  itmax/panel: $PRESSURE_ITMAX_PER_PANEL"

julia --project=. -t "$THREADS" examples/pitching_wing_pressure_replay.jl
