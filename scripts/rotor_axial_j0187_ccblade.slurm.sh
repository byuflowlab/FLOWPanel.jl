#!/usr/bin/env bash
# Submit manually from the repository root after creating data/rotor_axial_j0187_ccblade/.
# Account/QOS/modules are site/project specific: pass them to sbatch or uncomment locally.
#SBATCH --job-name=fp-axial-j0187
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --constraint=bigmem
#SBATCH --output=data/rotor_axial_j0187_ccblade/slurm-%j.out
#SBATCH --error=data/rotor_axial_j0187_ccblade/slurm-%j.err
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --qos=YOUR_QOS

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
THREADS="${SLURM_NTASKS:?SLURM_NTASKS must be set by Slurm}"
RUN_DIR=data/rotor_axial_j0187_ccblade
mkdir -p "$RUN_DIR"

# module load julia/<site-version>
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export JULIA_NUM_THREADS="$THREADS"
export MPLBACKEND=Agg

srun --ntasks=1 --cpus-per-task="$THREADS" julia --project=. -t "$THREADS" examples/rotor_axial_j0187_ccblade.jl
srun --ntasks=1 --cpus-per-task="$THREADS" julia --project=. -t "$THREADS" examples/rotor_axial_j0187_panel.jl
srun --ntasks=1 --cpus-per-task="$THREADS" julia --project=. -t "$THREADS" examples/rotor_axial_j0187_replay.jl
srun --ntasks=1 --cpus-per-task="$THREADS" julia --project=. -t "$THREADS" examples/rotor_axial_j0187_loading_comparison.jl
