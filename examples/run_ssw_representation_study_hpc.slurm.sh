#!/usr/bin/env bash
# Submit from the top level of the FLOWPanel.jl checkout.
#SBATCH --job-name=ssw017
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=ssw017-%j.out
#SBATCH --error=ssw017-%j.err

set -euo pipefail

THREADS=6
export JULIA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export MPLBACKEND=Agg
export MPLCONFIGDIR="${TMPDIR:-/tmp}/flowpanel-mpl"
export SSW_NO_PLOT=true
# Julia's Serialization format is not portable across the local 1.12 and
# cluster 1.11 runtimes. Always regenerate the six small Phase A snapshots on
# the execution host; the immutable CSV prediction remains the Phase B gate.
export SSWRP_REUSE_STATES=false
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

echo "BRAINSTORM 017 wing representation study"
echo "  repo:    $(pwd)"
echo "  threads: $THREADS"

# The immutable prediction is a hard gate for Phase B.
julia --project=. -t "$THREADS" examples/ssw_representation_probe.jl
test -s data/ssw_representation_probe/phase_b_prediction.csv
julia --project=. -t "$THREADS" examples/ssw_sheet_particle_split.jl
