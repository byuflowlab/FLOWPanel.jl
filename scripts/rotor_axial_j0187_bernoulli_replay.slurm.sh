#!/usr/bin/env bash
# Submit manually from the top level of the FLOWPanel.jl checkout.
# Account/QOS/modules are site/project specific: pass them to sbatch or uncomment locally.
#
# Replay-only job: reruns the corrected steady PressureBernoulli monitor over
# the saved rotor_axial_j0187_ccblade VTK state. It reads the existing run
# directory and writes only into data/rotor_axial_j0187_ccblade/bernoulli_replay/.
# It must never invoke rotor_axial_j0187_panel.jl or modify the saved run.
#SBATCH --job-name=fp-axial-bern-replay
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

THREADS=16

export JULIA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export MPLBACKEND=Agg

echo "FLOWPanel bernoulli-only axial replay"
echo "  repo:    $(pwd)"
echo "  threads: $THREADS"

if [[ ! -d data/rotor_axial_j0187_ccblade ]]; then
    echo "ERROR: data/rotor_axial_j0187_ccblade not found; run this from the checkout holding the saved axial run." >&2
    exit 2
fi

# module load julia/<site-version>

julia --project=. -t "$THREADS" examples/rotor_axial_j0187_bernoulli_replay.jl
