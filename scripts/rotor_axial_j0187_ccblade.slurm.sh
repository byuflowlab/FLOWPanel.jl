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
REPO_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$REPO_ROOT"
RUN_DIR=data/rotor_axial_j0187_ccblade
mkdir -p "$RUN_DIR"

# This source geometry is intentionally not substituted: using a different CST
# file would invalidate the intended panel/CCBlade comparison.
export CST_CSV="${CST_CSV:-$REPO_ROOT/examples/rotor_hover_scan/processed/dji9443_brainstorm_item003.csv}"
if [[ ! -f "$CST_CSV" ]]; then
    echo "ERROR: required DJI 9443 CST input is missing: $CST_CSV" >&2
    echo "Copy examples/rotor_hover_scan/processed/dji9443_brainstorm_item003.csv into the checkout, or submit with CST_CSV=/path/to/the/file." >&2
    exit 2
fi

# module load julia/<site-version>
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export JULIA_NUM_THREADS="$THREADS"
export MPLBACKEND=Agg

echo "THREADS: $THREADS"

# CCBlade writes these only after both ncrit=4 and ncrit=9 have completed.
# Preserve them after a wall-time timeout so a panel rerun does not repeat the
# expensive XFOIL polar sweeps. Set FORCE_CCBLADE=1 to rebuild deliberately.
ccblade_outputs=(
    "$RUN_DIR/rotor_hover_ccblade_polars.csv"
    "$RUN_DIR/rotor_hover_ccblade_operating_point_validation_Vc4_J0p1867.csv"
    "$RUN_DIR/rotor_hover_ccblade_sectional_ncrit4_Vc4_J0p1867.csv"
    "$RUN_DIR/rotor_hover_ccblade_sectional_ncrit9_Vc4_J0p1867.csv"
)
ccblade_complete=true
for output in "${ccblade_outputs[@]}"; do
    [[ -s "$output" ]] || ccblade_complete=false
done

if [[ "${FORCE_CCBLADE:-0}" == "1" || "$ccblade_complete" != true ]]; then
    echo "Running CCBlade/XFOIL polar generation."
    OMP_NUM_THREADS="$THREADS" julia --project=. -t "$THREADS" examples/rotor_axial_j0187_ccblade.jl
else
    echo "Reusing completed CCBlade polar and sectional outputs; skipping XFOIL."
fi

OMP_NUM_THREADS="$THREADS" julia --project=. -t $THREADS examples/rotor_axial_j0187_panel.jl
OMP_NUM_THREADS="$THREADS" julia --project=. -t $THREADS examples/rotor_axial_j0187_replay.jl
OMP_NUM_THREADS="$THREADS" julia --project=. -t $THREADS examples/rotor_axial_j0187_loading_comparison.jl
