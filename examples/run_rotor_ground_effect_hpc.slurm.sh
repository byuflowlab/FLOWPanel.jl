#!/usr/bin/env bash
#SBATCH --job-name=fp-022-ige
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# BRAINSTORM/022 — rotor hover in ground effect (driver:
# examples/rotor_hover_ground_effect.jl, fork of RHPC @ b251071).
#
# Submit ONE case per job, from the top level of the cluster checkout:
#
#   sbatch --job-name=fp-022-oge-fine   examples/run_rotor_ground_effect_hpc.slurm.sh p022_oge_fine
#   sbatch --job-name=fp-022-oge-coarse --time=24:00:00 examples/run_rotor_ground_effect_hpc.slurm.sh p022_oge_coarse
#   sbatch --job-name=fp-022-ige-fine   examples/run_rotor_ground_effect_hpc.slurm.sh p022_ige_fine
#   sbatch --job-name=fp-022-ige-coarse --time=24:00:00 examples/run_rotor_ground_effect_hpc.slurm.sh p022_ige_coarse
#
# Phase 3 below-ground policy arms (Ryan 2026-08-20, ruling 6), 48 h each:
#   sbatch --job-name=fp-022-ppdamp  --time=48:00:00 examples/run_rotor_ground_effect_hpc.slurm.sh p022_ige_coarse_damp
#   sbatch --job-name=fp-022-pptrunc --time=48:00:00 examples/run_rotor_ground_effect_hpc.slurm.sh p022_ige_coarse_trunc
#
# Phase 5 h/R sweep (NOT yet submitted; gated on the Phase 2 h/R=1 harvest).
# 48 h, not 24 h: the coarse rung measured ~87 s/step average (13207680), so the
# 1007-step schedule needs 25-31 h, and warm-start with GROUND_ENABLE=true errors
# out by design, so a walled run cannot be chained -- only re-run cold.
#   sbatch --job-name=fp-022-hr05 --time=48:00:00 examples/run_rotor_ground_effect_hpc.slurm.sh p022_hr05
#   sbatch --job-name=fp-022-hr15 --time=48:00:00 examples/run_rotor_ground_effect_hpc.slurm.sh p022_hr15
#   sbatch --job-name=fp-022-hr20 --time=48:00:00 examples/run_rotor_ground_effect_hpc.slurm.sh p022_hr20
#
# `logs/slurm/` must exist before submission. Memory: 45_185_ct4 dense
# Backslash G ~10 GB; the ground adds only an O(N) diagonal solve and a few
# thousand extra FMM sources — 64G has ample headroom.
#
# OPS RULE (018 lesson): every defining knob is exported UNCONDITIONALLY in
# its case arm — one flat block per tag, no ${VAR:-default} on defining knobs
# — so an inherited environment value can never silently win. Always verify
# the emitted *_case_metadata.toml and the driver banner after landing.

set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

CASE="${1:-}"
[[ -n "$CASE" ]] || { echo "ERROR: pass a case tag, e.g. p022_ige_fine" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg

# BRAINSTORM/025 (2026-08-20): the codebase default filament regularization
# changed to Gaussian. Pin this campaign to the legacy Vatistas family so
# results stay comparable; override at submission via --export if intended.
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-vatistas}"
# Manifest pins julia 1.11.7 (018 ops: a 1.12 re-resolve killed a whole trio).
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

# ---- 022 fixed operating point + carrier (identical for every case) ---------
# Operating point (Ryan 2026-08-18): RPM 6000, rho 1.16, R 0.1195. Carrier =
# 018 production carrier p018_ufront_n1_s040_visc knobs (sigma 0.04R, N=1
# uniform-D Das D=3.4sigma, stock rlxf 0.3, spreading-only viscosity).
export RPM=6000
export RHO=1.16
export ROTOR_R=0.1195
export RHPC_FORMULATION=velocity
export NT=36
export TRUNCATION_DEPTH_R=4
export RELAX_RLXF=0.3
export PARTICLE_SHEDDING=sigma_overlap
export OVERLAP=2.75
export P_PER_STEP=12
export MERGE_R_FACTOR=0.0055
export NWAKEROWS=1
export DAS_UNIFORM_DSIGMA=3.4
export CORE_SPREADING_ACTIVE=true
export WAKE_CORE_BETA=1e9
export BERNOULLI_ONLY=true
export RUN_KJ=false
export SAVE_VTK=true
# Validated staged-startup recipe (Item 005 stable-wake baseline, as in 018).
export SPINUP_REVS=1.5
export SPINUP_START_FRACTION=0.4
export MAGVINF_START=0.0
export MAGVINF_PEAK=5.0
export MAGVINF_END=0.0
export FREESTREAM_RAMP_REVS=1.0
export FREESTREAM_HOLD_REVS=1.5
export FREESTREAM_WITHDRAW_REVS=4
export SETTLE_REVS=20               # IGE may settle slower; judge per-rev blocks
export NREVS=10                     # schedule total wins
export CONVERGENCE_REVS=10
export CONVERGENCE_MEAN_TOL=0.005
export CONVERGENCE_PTP_TOL=0.02
# GS instrumentation always on for 022 (health gate: gs_nonconverged must be 0).
export GS_LOG=true
export GS_MAX_OUTER=50
export GS_TOL=1e-8

# ---- Per-case knobs (flat, unconditional) ------------------------------------
case "$CASE" in
  p022_oge_fine)
    export RHPC_MESH=45_185_ct4
    export GROUND_ENABLE=false
    export TRUNC_RADIUS_R=1.5 ;;
  p022_oge_coarse)
    export RHPC_MESH=56_57
    export GROUND_ENABLE=false
    export TRUNC_RADIUS_R=1.5 ;;
  p022_ige_fine)
    export RHPC_MESH=45_185_ct4
    export GROUND_ENABLE=true
    export GROUND_H_R=1.0
    export GROUND_RADIUS_R=4.0
    export GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=none
    export GROUND_DAMP_BAND_R=0
    export TRUNC_RADIUS_R=3.0 ;;
  p022_ige_coarse)
    export RHPC_MESH=56_57
    export GROUND_ENABLE=true
    export GROUND_H_R=1.0
    export GROUND_RADIUS_R=4.0
    export GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=none
    export GROUND_DAMP_BAND_R=0
    export TRUNC_RADIUS_R=3.0 ;;
  # --- Phase 3 (Ryan 2026-08-20, ruling 6): below-ground policy arms, 48 h.
  # _damp keeps policy=none ON PURPOSE so the below-ground census still counts
  # pass-throughs (the A/B observable vs cancelled 13207682's logged census).
  # _trunc = production candidate: damping + truncation floored at the ground.
  p022_ige_coarse_damp)
    export RHPC_MESH=56_57
    export GROUND_ENABLE=true
    export GROUND_H_R=1.0
    export GROUND_RADIUS_R=4.0
    export GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=none
    export GROUND_DAMP_BAND_R=0.1
    export TRUNC_RADIUS_R=3.0 ;;
  p022_ige_coarse_trunc)
    export RHPC_MESH=56_57
    export GROUND_ENABLE=true
    export GROUND_H_R=1.0
    export GROUND_RADIUS_R=4.0
    export GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=truncate
    export GROUND_DAMP_BAND_R=0.1
    export TRUNC_RADIUS_R=3.0 ;;
  # --- Phase 5 h/R sweep (coarse rung). TRUNCATION_DEPTH_R is measured from the
  # rotor plane, not the ground, so it must track h/R to hold the below-ground
  # allowance fixed at 3R: depth = GROUND_H_R + 3. These re-exports override the
  # shared block's TRUNCATION_DEPTH_R=4 (they run after it). h/R=1.0 is not
  # repeated here -- Phase 5 reuses p022_ige_coarse for that point.
  p022_hr05)
    export RHPC_MESH=56_57
    export GROUND_ENABLE=true
    export GROUND_H_R=0.5
    export TRUNCATION_DEPTH_R=3.5
    export GROUND_RADIUS_R=4.0
    export GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=none
    export GROUND_DAMP_BAND_R=0   # update to Phase-3 verdict before submitting
    export TRUNC_RADIUS_R=3.0 ;;
  p022_hr15)
    export RHPC_MESH=56_57
    export GROUND_ENABLE=true
    export GROUND_H_R=1.5
    export TRUNCATION_DEPTH_R=4.5
    export GROUND_RADIUS_R=4.0
    export GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=none
    export GROUND_DAMP_BAND_R=0   # update to Phase-3 verdict before submitting
    export TRUNC_RADIUS_R=3.0 ;;
  p022_hr20)
    export RHPC_MESH=56_57
    export GROUND_ENABLE=true
    export GROUND_H_R=2.0
    export TRUNCATION_DEPTH_R=5.0
    export GROUND_RADIUS_R=4.0
    export GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=none
    export GROUND_DAMP_BAND_R=0   # update to Phase-3 verdict before submitting
    export TRUNC_RADIUS_R=3.0 ;;
  *) echo "ERROR: unknown 022 case tag '$CASE' (Phase 3/4 rungs get their own explicit arms)" >&2; exit 2 ;;
esac

export RUN_NAME="$CASE"
[[ -n "${P022_RUN_NAME:-}" ]] && export RUN_NAME="$P022_RUN_NAME"

# Wipe-or-keep previous attempt (HPC disk policy; RHPC_KEEP_PREV=true to keep).
if [[ -d "data/$RUN_NAME" ]]; then
  if [[ "${RHPC_KEEP_PREV:-false}" == "true" ]]; then
    rm -rf "data/$RUN_NAME.prev"
    mv "data/$RUN_NAME" "data/$RUN_NAME.prev"
    echo "  previous run preserved at data/$RUN_NAME.prev"
  else
    rm -rf "data/$RUN_NAME"
  fi
fi

echo "BRAINSTORM/022 rotor hover ground effect — case $CASE"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"
echo "  mesh:$RHPC_MESH RPM:$RPM rho:$RHO R:$ROTOR_R NT:$NT depth:${TRUNCATION_DEPTH_R}R trunc_radius:${TRUNC_RADIUS_R}R rlxf:$RELAX_RLXF"
echo "  ground:${GROUND_ENABLE} h/R:${GROUND_H_R:-na} disc:${GROUND_RADIUS_R:-na}R panel:${GROUND_PANEL_LENGTH_R:-na}R policy:${GROUND_PARTICLE_POLICY:-na} damp_band:${GROUND_DAMP_BAND_R:-0}R"
echo "  overlap:$OVERLAP pps:$P_PER_STEP merge_r:$MERGE_R_FACTOR nrows:$NWAKEROWS das_uniform:$DAS_UNIFORM_DSIGMA visc:$CORE_SPREADING_ACTIVE settle:$SETTLE_REVS"
echo "  gs_log:$GS_LOG gs_cap:$GS_MAX_OUTER gs_tol:$GS_TOL"
echo "  started $(date '+%F %T')"

julia --project=. -t "$THREADS" examples/rotor_hover_ground_effect.jl

echo "Case $CASE finished ($(date '+%F %T'))."
echo "Artifacts: data/$RUN_NAME/${RUN_NAME}_CT_vs_rev.csv"
echo "           data/$RUN_NAME/${RUN_NAME}_CT_per_rev.csv"
echo "           data/$RUN_NAME/${RUN_NAME}_case_metadata.toml"
echo "           data/$RUN_NAME/${RUN_NAME}_gs_convergence.csv (IGE)"
echo "           data/$RUN_NAME/${RUN_NAME}_ground_diagnostics.csv (IGE)"
