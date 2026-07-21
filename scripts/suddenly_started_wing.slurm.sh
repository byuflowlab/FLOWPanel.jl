#!/usr/bin/env bash
#SBATCH --job-name=fp-ssw-convergence
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Submit from the top level of the FLOWPanel.jl checkout with:
#   sbatch scripts/suddenly_started_wing.slurm.sh
#
# Account, QOS, constraint, and Julia module directives are intentionally left
# to the allocation/site configuration.

set -euo pipefail

THREADS=16

export JULIA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

export SSW_MODE="${SSW_MODE:-convergence}"
export SSW_OUTPUT="${SSW_OUTPUT:-data/suddenly_started_wing}"
export SSW_BACKEND="${SSW_BACKEND:-fmm}"
export SSW_VERBOSE="${SSW_VERBOSE:-true}"
export SSW_BACKSLASH_MAX_PANELS="${SSW_BACKSLASH_MAX_PANELS:-20000}"
export SSW_RESUME="${SSW_RESUME:-true}"
export SAVE_VTK="${SAVE_VTK:-true}"

# Keep plotting headless on compute nodes.
export GKSwstype="${GKSwstype:-100}"
export MPLBACKEND="${MPLBACKEND:-Agg}"

echo "FLOWPanel suddenly-started-wing convergence"
echo "  repo:      $(pwd)"
echo "  threads:   $THREADS"
echo "  mode:      $SSW_MODE"
echo "  backend:   $SSW_BACKEND"
echo "  output:    $SSW_OUTPUT"
echo "  save VTK:  $SAVE_VTK"
echo "  dense cutoff: $SSW_BACKSLASH_MAX_PANELS panels"

julia --project=. -t "$THREADS" examples/suddenly_started_wing.jl
