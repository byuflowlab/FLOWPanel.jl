#!/usr/bin/env bash
#SBATCH --job-name=ssw017relax
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=ssw017relax-%j.out
#SBATCH --error=ssw017relax-%j.err

set -euo pipefail
export JULIA_NUM_THREADS=6
export OMP_NUM_THREADS=6
export OPENBLAS_NUM_THREADS=6
export BLAS_NUM_THREADS=6
export MKL_NUM_THREADS=6
export MPLBACKEND=Agg
export SSW_NO_PLOT=true
export SSWSPS_OUTPUT=data/ssw_sheet_particle_split_relaxation
export SSWSPS_RESUME=false
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"
julia --project=. -t 6 examples/ssw_sheet_particle_split.jl
