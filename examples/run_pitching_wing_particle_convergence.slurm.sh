#!/usr/bin/env bash
#SBATCH --job-name=fp-pitch-particle-conv
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=64g
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Submit from the top level of the FLOWPanel.jl checkout with:
#   sbatch examples/run_pitching_wing_particle_convergence.slurm.sh
#
# Cases run in the same adaptive order as the panel-wake convergence study:
# wake-length refinement first, followed by timestep refinement. If the driver
# stops before the allocation expires, submit this same command again;
# completed cases and cycle checkpoints will be reused automatically.
# Account, QOS, constraint, and Julia module directives are intentionally left
# to the allocation/site configuration.

set -euo pipefail

THREADS=48
RUN_ROOT="${PITCHCONV_ROOT:-data/pitching_wing_particle_convergence}"

export JULIA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

export PITCHCONV_ROOT="$RUN_ROOT"
export PITCHCONV_WAKE_MODEL=particle
export PITCHCONV_SECTION_ETA="${PITCHCONV_SECTION_ETA:-0.25,0.475,0.8,0.9}"
export PITCHCONV_BASELINE_N_SPAN="${PITCHCONV_BASELINE_N_SPAN:-13}"
export PITCHCONV_N_SPAN_VALUES="${PITCHCONV_N_SPAN_VALUES:-26,52,104}"
export PITCHCONV_WAKE_VALUES="${PITCHCONV_WAKE_VALUES:-1,2,4,8,16}"
export PITCHCONV_DT_VALUES="${PITCHCONV_DT_VALUES:-0.25,0.125,0.0625,0.03125}"
export PITCHCONV_BASELINE_C_PER_DT="${PITCHCONV_BASELINE_C_PER_DT:-0.5}"
export PITCHCONV_INITIAL_CYCLES="${PITCHCONV_INITIAL_CYCLES:-3}"
export PITCHCONV_MAX_CYCLES="${PITCHCONV_MAX_CYCLES:-10}"
export PITCHCONV_AMPLITUDE_TOL="${PITCHCONV_AMPLITUDE_TOL:-0.01}"
export PITCHCONV_FREQUENCY_TOL="${PITCHCONV_FREQUENCY_TOL:-0.005}"
export PITCHCONV_TIME_SAFETY="${PITCHCONV_TIME_SAFETY:-1.5}"
export PITCHCONV_TIME_RESERVE_SECONDS="${PITCHCONV_TIME_RESERVE_SECONDS:-900}"

echo "FLOWPanel adaptive pitching-wing particle-wake convergence"
echo "  repo:             $(pwd)"
echo "  output:           $PITCHCONV_ROOT"
echo "  wake model:       $PITCHCONV_WAKE_MODEL (PanelParticleWake)"
echo "  section etas:     $PITCHCONV_SECTION_ETA"
echo "  span panels:      $PITCHCONV_BASELINE_N_SPAN -> $PITCHCONV_N_SPAN_VALUES"
echo "  threads:          $THREADS"
echo "  wake spans:       $PITCHCONV_WAKE_VALUES"
echo "  c/dt baseline:    $PITCHCONV_BASELINE_C_PER_DT"
echo "  c/dt refinements: $PITCHCONV_DT_VALUES"
echo "  cycle range:      $PITCHCONV_INITIAL_CYCLES..$PITCHCONV_MAX_CYCLES"
echo "  amplitude tol:    $PITCHCONV_AMPLITUDE_TOL"
echo "  frequency tol:    $PITCHCONV_FREQUENCY_TOL"

julia --project=. -t "$THREADS" examples/pitching_wing_convergence.jl
