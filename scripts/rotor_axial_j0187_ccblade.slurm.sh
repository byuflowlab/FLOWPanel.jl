#!/usr/bin/env bash
# Submit manually from the repository root after creating data/rotor_axial_j0187_ccblade/.
# Account/QOS/modules are site/project specific: pass them to sbatch or uncomment locally.
#SBATCH --time=08:00:00
#SBATCH --ntasks=48
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH -J "fp-axial-j0187"
#SBATCH --output=ccblade_comp_slurm-%j.out
#SBATCH --error=ccblade_comp_slurm-%j.err

set -euo pipefail
THREADS=48
RUN_DIR=data/rotor_axial_j0187_ccblade
mkdir -p "$RUN_DIR"

# module load julia/<site-version>
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export JULIA_NUM_THREADS="$THREADS"
export MPLBACKEND=Agg

echo "THREADS: $THREADS"
OMP_NUM_THREADS="$THREADS" julia --project=. -t $THREADS examples/rotor_axial_j0187_ccblade.jl
OMP_NUM_THREADS="$THREADS" julia --project=. -t $THREADS examples/rotor_axial_j0187_panel.jl
OMP_NUM_THREADS="$THREADS" julia --project=. -t $THREADS examples/rotor_axial_j0187_replay.jl
OMP_NUM_THREADS="$THREADS" julia --project=. -t $THREADS examples/rotor_axial_j0187_loading_comparison.jl
