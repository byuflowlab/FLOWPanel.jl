#!/usr/bin/env bash
# Rotor-hover STABLE-WAKE driver (BRAINSTORM/005_rotor_hover_stable_converged_particle_wake.md).
#
# Goal of item 005: a STABLE (non-oscillating) particle-wake CT history -- a flat
# plateau at ANY value is the win; magnitude is 002/003/004's job. The lever is a
# smoothly staged startup: gradual RPM spin-up + a convecting freestream that is
# ramped up (to convect the shed wake clear and build a wake column) then slowly
# WITHDRAWN so wake self-induction sustains the column into hover.
#
# Phase 2 (plans/20260617_brainstorm_005.md): CHARACTERIZE the post-withdrawal
# hover CT oscillation (period/amplitude consistency, damping, truncation
# artifact, ramp sensitivity) and test damping schemes. Driver runs
# examples/rotor_hover_pressure_comparison.jl over an EXPERIMENT matrix. Each job
# writes data/<RUN_NAME>/<RUN_NAME>_CT_vs_rev.csv + per-step monitor force CSV +
# VTK; read the "Item 005 plateau diagnostics" block at the end of each log.
#
# EXPERIMENTS (set EXPERIMENT=...):
#   e2_depth     truncation-depth control (1/2/4) -- is it a boundary artifact?
#   e3_withdraw  withdraw-rate sweep (2.5/6/12) + pure-spinup control
#   e4_damping   all damping schemes on the chosen best baseline
#   e1_long      single long run (~30 settle rev) for period/damping
#
# THREAD BUDGET = 36 total. OMP/BLAS thread counts are tied to THREADS_PER_RUN
# (julia --threads N + OMP_NUM_THREADS=N). Keep CONCURRENCY*THREADS_PER_RUN <= 36.
#   3-config sweep:  EXPERIMENT=e2_depth CONCURRENCY=3 THREADS_PER_RUN=12 bash examples/run_rotor_hover_stable_wake.sh
#   long single run: EXPERIMENT=e1_long  CONCURRENCY=1 THREADS_PER_RUN=12 bash examples/run_rotor_hover_stable_wake.sh
#   e4 (13 jobs):    EXPERIMENT=e4_damping CONCURRENCY=6 THREADS_PER_RUN=6 bash examples/run_rotor_hover_stable_wake.sh
#
# Set the best baseline for e4 via env, e.g.
#   EXPERIMENT=e4_damping TRUNCATION_DEPTH_R=4 FREESTREAM_WITHDRAW_REVS=12 SETTLE_REVS=12 ...
#
# Written for stock macOS bash 3.2 (indexed arrays only, no `wait -n`).
set -euo pipefail

cd "$(dirname "$0")/.."   # repo root

# ---- Shared run config (env-overridable) ------------------------------------
RHPC_MESH="${RHPC_MESH:-40_40}"       # cheap mesh for tuning
NT="${NT:-18}"
RPM="${RPM:-6000}"
SAVE_VTK="${SAVE_VTK:-true}"          # CT CSV + particle .vtp histories
RUN_KJ="${RUN_KJ:-false}"             # off keeps cost down
BERNOULLI_ONLY="${BERNOULLI_ONLY:-true}"  # cheap monitor (no CG solve / no Hessian)

# Thread budget: 36 total. Default = 3 parallel jobs x 12 threads. OMP/BLAS thread
# counts are tied to THREADS_PER_RUN per the user directive (julia --threads N +
# OMP_NUM_THREADS=N). Keep CONCURRENCY * THREADS_PER_RUN <= 36.
CONCURRENCY="${CONCURRENCY:-3}"
THREADS_PER_RUN="${THREADS_PER_RUN:-12}"

# Which experiment to run (e1_long | e2_depth | e3_withdraw | e4_damping).
EXPERIMENT="${EXPERIMENT:-e2_depth}"
STABLE_PREFIX="${STABLE_PREFIX:-stable_wake_${EXPERIMENT}}"
LOG_DIR="data/${STABLE_PREFIX}_logs"

# ---- Best-known staged-startup baseline (applied to EVERY job) ---------------
# Validated in Phase 1: clean positive-CT startup, flat through ramp+hold.
# Per-experiment blocks below override SETTLE_REVS / TRUNCATION_DEPTH_R / NT and
# set the JOB_TABLE. Viscous core spreading + SFS stay ON (defaults in the .jl).
RECIPE_ENV=(
    "SPINUP_REVS=${SPINUP_REVS:-1.5}"
    "SPINUP_START_FRACTION=${SPINUP_START_FRACTION:-0.4}"
    "MAGVINF_PEAK=${MAGVINF_PEAK:-5.0}"
    "MAGVINF_END=${MAGVINF_END:-0.0}"
    "FREESTREAM_RAMP_REVS=${FREESTREAM_RAMP_REVS:-1.0}"
    "FREESTREAM_HOLD_REVS=${FREESTREAM_HOLD_REVS:-1.5}"
    "FREESTREAM_WITHDRAW_REVS=${FREESTREAM_WITHDRAW_REVS:-3.0}"
)

# ---- Per-experiment settle/truncation + job matrix --------------------------
# JOB_TABLE entries: "name|ENV1=val ENV2=val" (per-job single-variable overrides).
case "$EXPERIMENT" in
  e1_long)   # E1: many cycles to judge period consistency + damping (Q1/Q3)
    RECIPE_ENV+=("SETTLE_REVS=${SETTLE_REVS:-30.0}" "TRUNCATION_DEPTH_R=${TRUNCATION_DEPTH_R:-4}")
    NT="${NT_OVERRIDE:-12}"
    JOB_TABLE=( "e1_long|" )
    : "${CONCURRENCY:=1}"
    ;;
  e2_depth)  # E2: is the oscillation a truncation artifact? (Q4) -> 3 parallel
    RECIPE_ENV+=("SETTLE_REVS=${SETTLE_REVS:-12.0}")
    JOB_TABLE=(
        "depth_1|TRUNCATION_DEPTH_R=1"
        "depth_2|TRUNCATION_DEPTH_R=2"
        "depth_4|TRUNCATION_DEPTH_R=4"
    )
    ;;
  e3_withdraw)  # E3: does gentler withdrawal reduce the oscillation? (Q2)
    RECIPE_ENV+=("SETTLE_REVS=${SETTLE_REVS:-12.0}" "TRUNCATION_DEPTH_R=${TRUNCATION_DEPTH_R:-4}")
    JOB_TABLE=(
        "withdraw_2p5|FREESTREAM_WITHDRAW_REVS=2.5"
        "withdraw_6|FREESTREAM_WITHDRAW_REVS=6.0"
        "withdraw_12|FREESTREAM_WITHDRAW_REVS=12.0"
        "spinup_only|MAGVINF_PEAK=0.0001 FREESTREAM_WITHDRAW_REVS=0.5"
    )
    ;;
  e4_damping)  # E4: all damping schemes on the best non-damping baseline
    # Set the best E1-E3 baseline via env (SETTLE_REVS/TRUNCATION_DEPTH_R/
    # FREESTREAM_WITHDRAW_REVS) before running this experiment.
    RECIPE_ENV+=("SETTLE_REVS=${SETTLE_REVS:-12.0}" "TRUNCATION_DEPTH_R=${TRUNCATION_DEPTH_R:-4}")
    JOB_TABLE=(
        "baseline|"
        "nu_x3|WAKE_NU_FACTOR=3"
        "nu_x10|WAKE_NU_FACTOR=10"
        "sfs_strong|SFS_MAXC=2.0 SFS_RLXF=0.02"
        "sfs_off|SFS_OFF=true"
        "overlap_hi|OVERLAP=4.0"
        "corebeta_hi|WAKE_CORE_BETA=2.0"
        "nt_hi|NT=36"
        "pps_hi|P_PER_STEP=4"
        "relax_off|PARTICLE_RELAX=false"
        "merge_aggressive|MERGE_R_FACTOR=0.05"
        "kernoff_hi|KERNELOFFSET_TARGETS=3e-3"
        "bound_rlx|BOUND_STRENGTH_RLX=0.5"   # E4.8 body-strength low-pass (numerical, last resort)
    )
    ;;
  *)
    echo "Unknown EXPERIMENT=$EXPERIMENT (use e1_long|e2_depth|e3_withdraw|e4_damping)" >&2
    exit 1
    ;;
esac

# Optional subset: JOBS="depth_2 depth_4" selects matching job names.
JOBS="${JOBS:-}"

# ---- Job pool (bash 3.2 safe) -----------------------------------------------
mkdir -p "$LOG_DIR"
pids=()

wait_for_slot () {
    while [ "${#pids[@]}" -ge "$CONCURRENCY" ]; do
        local oldest="${pids[0]}"
        wait "$oldest" || echo "WARNING: job pid $oldest exited non-zero" >&2
        pids=("${pids[@]:1}")
    done
}

run_one () {
    local name=$1 overrides=$2
    local run_name="${STABLE_PREFIX}_${name}"
    local logf="${LOG_DIR}/${name}.log"
    rm -rf "data/${run_name}"     # disk hygiene: overwrite previous attempt
    echo "=== launching ${name}  (RUN_NAME=${run_name}, threads=${THREADS_PER_RUN})  overrides: ${overrides:-<none>} ==="
    # OMP/BLAS thread counts tied to this run's Julia thread count (user directive).
    # shellcheck disable=SC2086
    env JULIA_NUM_THREADS="$THREADS_PER_RUN" \
        OMP_NUM_THREADS="$THREADS_PER_RUN" \
        OPENBLAS_NUM_THREADS="$THREADS_PER_RUN" \
        MKL_NUM_THREADS="$THREADS_PER_RUN" \
        VECLIB_MAXIMUM_THREADS="$THREADS_PER_RUN" \
        RUN_NAME="$run_name" \
        RHPC_MESH="$RHPC_MESH" NT="$NT" RPM="$RPM" \
        SAVE_VTK="$SAVE_VTK" RUN_KJ="$RUN_KJ" BERNOULLI_ONLY="$BERNOULLI_ONLY" \
        "${RECIPE_ENV[@]}" \
        $overrides \
        julia --project --threads "$THREADS_PER_RUN" examples/rotor_hover_pressure_comparison.jl \
        > "$logf" 2>&1 &
    pids+=("$!")
}

selected () {
    [ -z "$JOBS" ] && { echo yes; return; }
    local j
    for j in $JOBS; do [ "$j" = "$1" ] && { echo yes; return; }; done
    echo no
}

echo "Stable-wake sweep: mesh=${RHPC_MESH} nt=${NT} rpm=${RPM}"
echo "Recipe: ${RECIPE_ENV[*]}"
echo "Concurrency=${CONCURRENCY} threads/run=${THREADS_PER_RUN}  logs -> ${LOG_DIR}"
echo

for entry in "${JOB_TABLE[@]}"; do
    name="${entry%%|*}"
    overrides="${entry#*|}"
    [ "$(selected "$name")" = yes ] || continue
    wait_for_slot
    run_one "$name" "$overrides"
done

if [ "${#pids[@]}" -gt 0 ]; then
    for p in "${pids[@]}"; do
        wait "$p" || echo "WARNING: job pid $p exited non-zero" >&2
    done
fi

echo
echo "All stable-wake runs complete. Per-job plateau diagnostics:"
for entry in "${JOB_TABLE[@]}"; do
    name="${entry%%|*}"
    [ "$(selected "$name")" = yes ] || continue
    logf="${LOG_DIR}/${name}.log"
    echo "--- ${name} ---"
    grep -A4 "Item 005 plateau diagnostics" "$logf" 2>/dev/null || echo "  (no diagnostics block found; check ${logf})"
done
