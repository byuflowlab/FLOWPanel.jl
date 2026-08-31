#!/usr/bin/env bash
# Local smoke tests for the 022 multi-rotor driver + launcher arms.
# Coarse mesh, NT=6, ~9-step schedule, 4 threads (local thread cap).
# Usage: bash scripts/p022m_smoke_local.sh <case>
#   case in: 1r_oge  2r_oge  2r_ige  4r_oge  4r_ige
set -euo pipefail
cd "$(dirname "$0")/.."

CASE="${1:?pass a smoke case: 1r_oge 2r_oge 2r_ige 4r_oge 4r_ige}"

export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 BLAS_NUM_THREADS=1
export RHPC_MESH=40_40
export NT=6
export BERNOULLI_ONLY=true
export SAVE_VTK=true
export RUN_MONITORS=true
export FLOWPANEL_FILAMENT_REG=gaussian
# Production wake carrier knobs, so the smoke exercises the same Das/sigma
# code paths the launcher will (uniform-D Das over N bodies).
export PARTICLE_SHEDDING=sigma_overlap
export OVERLAP=2.75
export P_PER_STEP=12
export NWAKEROWS=1
export DAS_UNIFORM_DSIGMA=3.4
export CORE_SPREADING_ACTIVE=true
export WAKE_CORE_BETA=1e9
export RELAX_RLXF=0.3
export MERGE_R_FACTOR=0.0055
export GS_LOG=true
# Short schedule: 0.5 + 0 + 0.5 + 0.5 = 1.5 revs x NT 6 = 9 steps, no spinup.
export SPINUP_REVS=0
export MAGVINF_PEAK=5.0
export FREESTREAM_RAMP_REVS=0.5
export FREESTREAM_HOLD_REVS=0
export FREESTREAM_WITHDRAW_REVS=0.5
export SETTLE_REVS=0.5
export NREVS=1.5
export ROTOR_SPACING_R=2.7

case "$CASE" in
  1r_oge) export NROTORS=1 GROUND_ENABLE=false TRUNC_RADIUS_R=1.5 TRUNCATION_DEPTH_R=4 ;;
  2r_oge) export NROTORS=2 GROUND_ENABLE=false TRUNC_RADIUS_R=1.5 TRUNCATION_DEPTH_R=4 ;;
  4r_oge) export NROTORS=4 GROUND_ENABLE=false TRUNC_RADIUS_R=1.5 TRUNCATION_DEPTH_R=4 ;;
  2r_ige|4r_ige)
    export NROTORS="${CASE:0:1}" GROUND_ENABLE=true
    export GROUND_H_R=1.5 TRUNCATION_DEPTH_R=4.5
    export GROUND_RADIUS_R=4.0 GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=none GROUND_DAMP_BAND_R=0
    export TRUNC_RADIUS_R=3.0 ;;
  *) echo "unknown smoke case $CASE" >&2; exit 2 ;;
esac

export RUN_NAME="smoke1t_$CASE"
rm -rf "data/$RUN_NAME"

echo "== smoke $CASE (NROTORS=$NROTORS ground=$GROUND_ENABLE) -> data/$RUN_NAME =="
julia --project=. -t 1 examples/rotor_hover_ground_effect.jl
echo "== smoke $CASE DONE =="
