#!/usr/bin/env bash
#SBATCH --job-name=fp-022m
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=96G
#SBATCH --time=72:00:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# BRAINSTORM/022 multi-rotor extension — 2 and 4 rotors, OGE and IGE (driver:
# examples/rotor_hover_ground_effect.jl with NROTORS support, 2026-08-25).
#
# Submit ONE case per job, from the top level of the cluster checkout, normally
# via the per-situation wrappers (they derive the h/R-coupled knobs):
#
#   scripts/p022m_submit_2r_oge.sh
#   scripts/p022m_submit_2r_ige.sh     # GROUND_H_R env-overridable, default 1.5
#   scripts/p022m_submit_4r_oge.sh
#   scripts/p022m_submit_4r_ige.sh
#
# Direct submission (OGE arms only; IGE arms REQUIRE the wrapper-derived
# IGE_H_R / IGE_TRUNC_DEPTH and fail fast without them):
#
#   sbatch --job-name=fp-022m-2r-oge --time=48:00:00 --mem=96G  examples/run_rotor_multi_ground_effect_hpc.slurm.sh p022m_2r_oge
#   sbatch --job-name=fp-022m-4r-oge --time=72:00:00 --mem=128G examples/run_rotor_multi_ground_effect_hpc.slurm.sh p022m_4r_oge
#
# `logs/slurm/` must exist before submission. Walltime is LOAD-BEARING: warm
# start is blocked for NROTORS>1 and for GROUND_ENABLE=true, so a walled run
# cannot be chained — only re-run cold. Verify s/step in the first hour.
#
# CARRIER (differs from run_rotor_ground_effect_hpc.slurm.sh ON PURPOSE,
# Ryan 2026-08-25): this campaign pins the 018/021 PRODUCTION carrier —
# FLOWPANEL_FILAMENT_REG=gaussian + the tuned per-pass FMM knobs (023/025,
# adopted 2026-08-21) — NOT the vatistas pin the single-rotor 022 arms carry
# for internal comparability. Do not mix the two families in one A/B.
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
[[ -n "$CASE" ]] || { echo "ERROR: pass a case tag, e.g. p022m_2r_ige" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg

# Production regularization family (defining knob for this campaign — pinned
# unconditionally, unlike the single-rotor 022 launcher's vatistas fallback).
export FLOWPANEL_FILAMENT_REG=gaussian
# Manifest pins julia 1.11.7 (018 ops: a 1.12 re-resolve killed a whole trio).
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

# ---- Fixed operating point + carrier (identical for every case) --------------
# Operating point: RPM 6000, rho 1.16, R 0.1195 (022). Wake carrier = 018
# production p018_ufront_n1_s040_visc knobs (sigma 0.04R, N=1 uniform-D Das
# D=3.4sigma, stock rlxf 0.3, spreading-only viscosity).
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
export BERNOULLI_ONLY=true          # mandatory for NROTORS>1 (and with ground)
export RUN_KJ=false
export SAVE_VTK=true
# Tuned per-pass FMM knobs (023/025 production since 2026-08-21).
export FMM_BODY_EXPANSION_ORDER=17
export FMM_BODY_ACCEPTANCE=0.7
export FMM_BODY_LEAF_SIZE=109
export FMM_WAKE_EXPANSION_ORDER=16
export FMM_WAKE_ACCEPTANCE=0.6
export FMM_WAKE_LEAF_SIZE=38
# Validated staged-startup recipe (Item 005 stable-wake baseline, as in 018).
export SPINUP_REVS=1.5
export SPINUP_START_FRACTION=0.4
export MAGVINF_START=0.0
export MAGVINF_PEAK=5.0
export MAGVINF_END=0.0
export FREESTREAM_RAMP_REVS=1.0
export FREESTREAM_HOLD_REVS=1.5
export FREESTREAM_WITHDRAW_REVS=4
export SETTLE_REVS=20
export NREVS=10                     # schedule total wins
export CONVERGENCE_REVS=10
export CONVERGENCE_MEAN_TOL=0.005
export CONVERGENCE_PTP_TOL=0.02
# GS instrumentation always on (health gate: gs_nonconverged must be 0). For
# multi-rotor OGE the block Gauss-Seidel couples the rotors even with no ground.
export GS_LOG=true
export GS_MAX_OUTER=50
export GS_TOL=1e-8
# Multi-rotor layout (Ryan 2026-08-25): 2.7R center-to-center, counter-rotating
# adjacent pairs (driver default directions).
export ROTOR_SPACING_R=2.7

# ---- Per-case knobs (flat, unconditional) ------------------------------------
# Mesh: RHPC_MESH is wrapper-overridable via P022M_MESH; the IGE default is the
# FINE mesh (45_185_ct4) because the coarse 56_57 blew up near the ground at
# rev 17.6 in single-rotor 022 (ledger). OGE defaults coarse (56_57) for cost.
# Confirm the mesh choice with Ryan at first production submission.
case "$CASE" in
  p022m_2r_oge)
    export NROTORS=2
    export RHPC_MESH="${P022M_MESH:-56_57}"
    export GROUND_ENABLE=false
    export TRUNC_RADIUS_R=1.5 ;;
  p022m_4r_oge)
    export NROTORS=4
    export RHPC_MESH="${P022M_MESH:-56_57}"
    export GROUND_ENABLE=false
    export TRUNC_RADIUS_R=1.5 ;;
  p022m_2r_ige)
    export NROTORS=2
    export RHPC_MESH="${P022M_MESH:-45_185_ct4}"
    export GROUND_ENABLE=true
    export GROUND_H_R="${IGE_H_R:?p022m IGE arms require IGE_H_R (use the submit wrapper)}"
    export TRUNCATION_DEPTH_R="${IGE_TRUNC_DEPTH:?p022m IGE arms require IGE_TRUNC_DEPTH = h/R + 3 (use the submit wrapper)}"
    export GROUND_RADIUS_R=4.0
    export GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=none
    export GROUND_DAMP_BAND_R=0.1 # Ryan 2026-08-26 (post Phase-3 autopsy): Phase-3 arms overflowed/blew up ~rev 15-18 with 32k below-ground particles; breach ruling applied — 0.1R damp band on IGE arms
    export TRUNC_RADIUS_R=3.0 ;;
  p022m_4r_ige)
    export NROTORS=4
    export RHPC_MESH="${P022M_MESH:-45_185_ct4}"
    export GROUND_ENABLE=true
    export GROUND_H_R="${IGE_H_R:?p022m IGE arms require IGE_H_R (use the submit wrapper)}"
    export TRUNCATION_DEPTH_R="${IGE_TRUNC_DEPTH:?p022m IGE arms require IGE_TRUNC_DEPTH = h/R + 3 (use the submit wrapper)}"
    export GROUND_RADIUS_R=4.0
    export GROUND_PANEL_LENGTH_R=0.15
    export GROUND_PARTICLE_POLICY=none
    export GROUND_DAMP_BAND_R=0.1 # Ryan 2026-08-26 (post Phase-3 autopsy): Phase-3 arms overflowed/blew up ~rev 15-18 with 32k below-ground particles; breach ruling applied — 0.1R damp band on IGE arms
    export TRUNC_RADIUS_R=3.0 ;;
  *) echo "ERROR: unknown 022m case tag '$CASE'" >&2; exit 2 ;;
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

echo "BRAINSTORM/022 multi-rotor ground effect — case $CASE"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"
echo "  nrotors:$NROTORS spacing:${ROTOR_SPACING_R}R (counter-rotating pairs)"
echo "  mesh:$RHPC_MESH RPM:$RPM rho:$RHO R:$ROTOR_R NT:$NT depth:${TRUNCATION_DEPTH_R}R trunc_radius:${TRUNC_RADIUS_R}R rlxf:$RELAX_RLXF"
echo "  ground:${GROUND_ENABLE} h/R:${GROUND_H_R:-na} disc:${GROUND_RADIUS_R:-na}R panel:${GROUND_PANEL_LENGTH_R:-na}R policy:${GROUND_PARTICLE_POLICY:-na} damp_band:${GROUND_DAMP_BAND_R:-0}R"
echo "  reg:$FLOWPANEL_FILAMENT_REG fmm_body:$FMM_BODY_EXPANSION_ORDER/$FMM_BODY_ACCEPTANCE/$FMM_BODY_LEAF_SIZE fmm_wake:$FMM_WAKE_EXPANSION_ORDER/$FMM_WAKE_ACCEPTANCE/$FMM_WAKE_LEAF_SIZE"
echo "  overlap:$OVERLAP pps:$P_PER_STEP merge_r:$MERGE_R_FACTOR nrows:$NWAKEROWS das_uniform:$DAS_UNIFORM_DSIGMA visc:$CORE_SPREADING_ACTIVE settle:$SETTLE_REVS"
echo "  gs_log:$GS_LOG gs_cap:$GS_MAX_OUTER gs_tol:$GS_TOL"
echo "  started $(date '+%F %T')"

julia --project=. -t "$THREADS" examples/rotor_hover_ground_effect.jl

echo "Case $CASE finished ($(date '+%F %T'))."
echo "Artifacts: data/$RUN_NAME/${RUN_NAME}_CT_vs_rev.csv (total + per-rotor columns)"
echo "           data/$RUN_NAME/${RUN_NAME}_CT_per_rev.csv"
echo "           data/$RUN_NAME/${RUN_NAME}_case_metadata.toml"
echo "           data/$RUN_NAME/${RUN_NAME}_gs_convergence.csv"
echo "           data/$RUN_NAME/${RUN_NAME}_ground_diagnostics.csv (IGE)"
