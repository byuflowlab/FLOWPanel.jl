#!/usr/bin/env bash
# Rotor-hover WAKE-DECAY AUDIT driver (BRAINSTORM/004_rotor_hover_wake_decay_audit.md).
#
# Runs examples/rotor_hover_pressure_comparison.jl once per audit variant: a fixed
# baseline with ONE wake-modeling mechanism toggled per run (diffusion, merging,
# SFS, core-size growth, body/wake coupling). Each run writes its own
# data/<RUN_NAME>/<RUN_NAME>_CT_vs_rev.csv plus per-step particle .vtp files; the
# aggregation step (examples/aggregate_wake_decay_audit.jl) then overlays the CT
# histories, extracts particle-count histories, and writes a summary table.
#
# EXECUTION MODES (workstation):
#   all-parallel:        CONCURRENCY=12 THREADS_PER_RUN=4  bash examples/run_rotor_hover_wake_decay_audit.sh
#   sequential, threaded: CONCURRENCY=1  THREADS_PER_RUN=48 bash examples/run_rotor_hover_wake_decay_audit.sh
#   any mix:             CONCURRENCY=4  THREADS_PER_RUN=12 bash examples/run_rotor_hover_wake_decay_audit.sh
#
# Quick smoke test (cheap, proves the pipeline end-to-end):
#   NREVS=2 JOBS="baseline merge_off" CONCURRENCY=2 THREADS_PER_RUN=2 \
#       bash examples/run_rotor_hover_wake_decay_audit.sh
#
# NOTE (sequencing with item 005): the audit is only meaningful once a STABLE,
# non-oscillating CT baseline exists. Item 005 determines that recipe (gradual RPM
# spin-up + ramped-then-withdrawn freestream + downstream particle truncation).
# Fill the STABLE-RECIPE CONFIG block below with 005's result BEFORE running the
# full matrix; until then this harness reproduces the plain (transient) baseline.
#
# Written for stock macOS bash 3.2 (indexed arrays only, no `wait -n`).
set -euo pipefail

cd "$(dirname "$0")/.."   # repo root

# ---- Shared run config (env-overridable) ------------------------------------
AUDIT_PREFIX="${AUDIT_PREFIX:-wake_decay_audit}"
RHPC_MESH="${RHPC_MESH:-40_40}"
NREVS="${NREVS:-15}"
NT="${NT:-36}"
RPM="${RPM:-6000}"
SAVE_VTK="${SAVE_VTK:-true}"          # needed: CT CSV + particle .vtp histories
RUN_KJ="${RUN_KJ:-false}"             # off keeps cost down (KJ cross-checked in 001)

CONCURRENCY="${CONCURRENCY:-1}"       # number of simultaneous julia processes
THREADS_PER_RUN="${THREADS_PER_RUN:-4}"   # JULIA_NUM_THREADS per process

# Aggregation output dir (separate from per-run data/<RUN_NAME> dirs).
AGG_DIR="data/${AUDIT_PREFIX}_aggregate"
LOG_DIR="${AGG_DIR}/logs"

# ---- STABLE-RECIPE CONFIG (fill in from item 005 before the full audit) ------
# These env assignments are appended to EVERY job so the whole matrix shares the
# same staged-startup baseline. Item 005 added the freestream-pulse and
# truncation ENV hooks used below (see BRAINSTORM/005 2026-06-17 section).
#
# IMPORTANT (005 result): the staged startup KILLS the impulsive startup ring,
# but withdrawing the freestream to hover does NOT produce a flat plateau -- a
# non-damping oscillation, intrinsic to the self-sustaining particle wake,
# re-appears once magVinf=0. So there is no fully-stable hover baseline to share;
# the recipe below is the best available (smooth through ramp+hold) and is the
# right starting point for these ablations, which are precisely the candidate
# dampers for that intrinsic oscillation. Judge each ablation by whether it
# SHRINKS the post-withdrawal (magVinf=0) settle-window oscillation amplitude.
STABLE_RECIPE_ENV=(
    "SPINUP_REVS=1.5"
    "SPINUP_START_FRACTION=0.4"
    "MAGVINF_PEAK=5.0"
    "MAGVINF_END=0.0"
    "FREESTREAM_RAMP_REVS=1.0"
    "FREESTREAM_HOLD_REVS=1.5"
    "FREESTREAM_WITHDRAW_REVS=3.0"
    "SETTLE_REVS=5.0"
    "TRUNCATION_DEPTH_R=4"
)

# ---- Audit job matrix: "name|ENV1=val ENV2=val" -----------------------------
# baseline = STABLE_RECIPE_ENV only; every other job = baseline + one toggle.
JOB_TABLE=(
    "baseline|"
    "merge_off|MERGE_PARTICLES=false"
    "sfs_off|SFS_OFF=true"
    "sfs_noclip|SFS_NO_BACKSCATTER_CLIP=true"
    "nu_low|WAKE_NU=1.0e-6"
    "nu_high|WAKE_NU=5.0e-5"
    "beta_nogrow|WAKE_CORE_BETA=1.0"
    "beta_high|WAKE_CORE_BETA=2.0"
    "body_on_wake_off|BODY_ON_WAKE=false"
    "pw_on_part_off|PANEL_WAKE_ON_PARTICLES=false"
    "no_hessian_to_part|PANEL_WAKE_HESSIAN_TO_PARTICLES=false"
    "part_hessian_off|PARTICLE_HESSIAN_SELF=false"
)

# Optional subset: JOBS="baseline merge_off" selects matching job names.
JOBS="${JOBS:-}"

# ---- Job pool (bash 3.2 safe) -----------------------------------------------
mkdir -p "$LOG_DIR"
pids=()

wait_for_slot () {
    # Block until fewer than CONCURRENCY background jobs remain by waiting on the
    # oldest tracked PID (FIFO drain; portable to bash 3.2 which lacks `wait -n`).
    while [ "${#pids[@]}" -ge "$CONCURRENCY" ]; do
        local oldest="${pids[0]}"
        wait "$oldest" || echo "WARNING: job pid $oldest exited non-zero" >&2
        pids=("${pids[@]:1}")
    done
}

run_one () {
    local name=$1 overrides=$2
    local run_name="${AUDIT_PREFIX}_${name}"
    local logf="${LOG_DIR}/${name}.log"
    # Disk hygiene: overwrite previous attempt of this variant in place.
    rm -rf "data/${run_name}"
    echo "=== launching ${name}  (RUN_NAME=${run_name}, threads=${THREADS_PER_RUN})  overrides: ${overrides:-<none>} ==="
    # shellcheck disable=SC2086
    env JULIA_NUM_THREADS="$THREADS_PER_RUN" \
        RUN_NAME="$run_name" \
        RHPC_MESH="$RHPC_MESH" NREVS="$NREVS" NT="$NT" RPM="$RPM" \
        SAVE_VTK="$SAVE_VTK" RUN_KJ="$RUN_KJ" \
        "${STABLE_RECIPE_ENV[@]}" \
        $overrides \
        julia --project examples/rotor_hover_pressure_comparison.jl \
        > "$logf" 2>&1 &
    pids+=("$!")
}

selected () {
    # Echo "yes" if $1 is in the JOBS allowlist (or JOBS is empty = all).
    [ -z "$JOBS" ] && { echo yes; return; }
    local j
    for j in $JOBS; do [ "$j" = "$1" ] && { echo yes; return; }; done
    echo no
}

echo "Wake-decay audit: mesh=${RHPC_MESH} nrevs=${NREVS} nt=${NT} rpm=${RPM}"
echo "Concurrency=${CONCURRENCY} threads/run=${THREADS_PER_RUN}  logs -> ${LOG_DIR}"
echo

for entry in "${JOB_TABLE[@]}"; do
    name="${entry%%|*}"
    overrides="${entry#*|}"
    [ "$(selected "$name")" = yes ] || continue
    wait_for_slot
    run_one "$name" "$overrides"
done

# Drain remaining jobs (guard empty array for bash 3.2 + set -u).
if [ "${#pids[@]}" -gt 0 ]; then
    for p in "${pids[@]}"; do
        wait "$p" || echo "WARNING: job pid $p exited non-zero" >&2
    done
fi

echo
echo "All audit runs complete. Aggregating..."
env AUDIT_PREFIX="$AUDIT_PREFIX" RPM="$RPM" NT="$NT" \
    julia --project examples/aggregate_wake_decay_audit.jl

echo
echo "Done. Outputs in ${AGG_DIR}/"
