#!/usr/bin/env bash
#SBATCH --job-name=p021-r12-champion-rerun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --array=0-13
#SBATCH --output=logs/slurm/r12-champion-%A_%a.out
#SBATCH --error=logs/slurm/r12-champion-%A_%a.err
# 021 R1-R2 champion-scheme re-run (Ryan 2026-09-19: "re-run all of them...
# re-tune under the new approach (saving different tuning parameters for
# different numbers of threads)... record the new accuracy... then the new
# benchmarks", approved for R1-R2 only). Array = (rung x thread ladder):
# tasks 0-6 = R1 at j 1/2/4/8/16/32/64, tasks 7-13 = R2 likewise — the full
# per-class thread ladder of the 2026-09-19 policy (budget caps prune inside
# the tuner; phase2.jl skips budgets with no tuned row gracefully).
#
# New scheme per (rung, j), each stage RE-TUNED at its own thread count:
#   1. fgstune    — FGS knob descent + tolerance STAIRCASE with
#                   sweep_order=dagteam (FGS_DAGTEAM_PRECISION, default f64 —
#                   the safe off-R4 rung, mathematically the lexicographic
#                   iterate; f32full only after a case's own staircase
#                   certifies, resubmit with FGS_DAGTEAM_PRECISION=f32full).
#                   bc_rel_l2 columns ARE the recorded new accuracy.
#   2. fgsprecond — stage-3 preconditioner sweep ladder (SWEEP_LADDER_1E6=1).
#   3. phase2 tuner — apply-knob descent over the FULL machine-class ladder
#                   0/16/32/64/128/500 GiB; budgets whose thread cap < j are
#                   skipped by the tuner (2 GiB/thread policy).
#   4. phase2      — full default CONFIGS table (backslash, krylov family,
#                   fgs, fgmres_fgs + nfcache variants) at the fresh knobs.
#
# Placement = champion interleave over socket 0 (nodes 0-3, all 64 cores the
# j=64 point needs); BLAS = julia threads (the historical phase-1/2 multi
# convention — kept so backslash_ldiv rows stay comparable; the FGS family is
# BLAS-insensitive and its tolerances are calibrated per environment anyway).
#
# Isolation: each task gets its own BENCH_CASE_ROOT and PHASE2_OUTDIR under
# the data root — the phase-1 knob CSVs carry no sweep_order column and
# stage3_winner selects on rung + knobs only, so per-run directories are the
# correctness boundary (see phase1_case.jl), and concurrent tasks must never
# share an NFS CSV (2026-08-18 append hazard).
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"

LADDER=(1 2 4 8 16 32 64)
RUNGS=(R1 R2)
J="${LADDER[$((SLURM_ARRAY_TASK_ID % 7))]}"
RUNG_SEL="${RUNGS[$((SLURM_ARRAY_TASK_ID / 7))]}"

NODE_CPUS=$(getconf _NPROCESSORS_ONLN)
ALLOC_CPUS=${SLURM_CPUS_ON_NODE:-0}
if [ "$NODE_CPUS" != "128" ] || [ "$ALLOC_CPUS" != "128" ]; then
  echo "ERROR: need an exclusive 128-core zen3 node (got node=$NODE_CPUS" >&2
  echo "       alloc=$ALLOC_CPUS); timings would not be comparable." >&2
  exit 1
fi

export JULIA_PKG_PRECOMPILE_AUTO=0
if command -v flock >/dev/null 2>&1; then
  ( if flock -E 99 -w 3600 9; then
      julia --project="$COLD_PROJECT" --startup-file=no \
        -e 'using Pkg; Pkg.precompile()' \
        || echo "WARNING: Pkg.precompile() failed; compiling in memory"
    else
      echo "WARNING: precompile lock unavailable; compiling in memory"
    fi
  ) 9>"$HOME/.julia/flowpanel-021-precompile.lock"
else
  echo "WARNING: flock(1) not found; skipping shared precompile"
fi

run="$COLD_DATA_ROOT/r12-champion-$RUNG_SEL-j$J-$SLURM_ARRAY_JOB_ID"
mkdir "$run"
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
numactl --hardware > "$run/numactl_hardware.txt" 2>&1 || true
echo "rung=$RUNG_SEL julia_threads=$J" > "$run/point.txt"

ILV="--interleave=0-3 --cpunodebind=0-3"
numactl $ILV numactl --show > "$run/numactl_show.txt" 2>&1 || true

export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3}"
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export JULIA_DEBUG="${JULIA_DEBUG:-loading}"

( while true; do
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] heartbeat: $RUNG_SEL j$J alive"
    sleep 300
  done ) &
HEARTBEAT_PID=$!
trap 'kill "$HEARTBEAT_PID" 2>/dev/null || true' EXIT

status() { printf '%s\n' "$2" > "$run/STATUS_$1"; }

# Shared env for every julia stage. BENCH_CASE_ROOT is per-task; KNOBS_MODE
# stays the default (= threading_mode), so writers and readers agree inside
# this task's tree. BLAS = j (see header).
common_env=(RUNG="$RUNG_SEL" EXPECT_JULIA_THREADS="$J" THREADING_MODE=multi
            OPENBLAS_NUM_THREADS="$J" OMP_NUM_THREADS="$J"
            CACHE_B=1 K_REPS=3 BENCH_CASE_ROOT="$run/case"
            FGS_SWEEP_ORDER=dagteam
            FGS_DAGTEAM_PRECISION="${FGS_DAGTEAM_PRECISION:-f64}")

jrun() { # jrun <logname> <extra env...> -- <script>
    local log="$1"; shift
    local extra=()
    while [ "$1" != "--" ]; do extra+=("$1"); shift; done
    shift
    env "${common_env[@]}" "${extra[@]}" numactl $ILV \
        julia --project="$COLD_PROJECT" --startup-file=no \
        --compiled-modules=existing -t "$J" "$1" > "$run/$log.log" 2>&1
}

ok=1
if jrun fgstune -- benchmark/rotor_hover_solver_phase1_fgstune.jl; then
    status fgstune ok
else
    status fgstune FAILED; ok=0
fi

if [ "$ok" = 1 ]; then
    if jrun fgsprecond SWEEP_LADDER_1E6=1 -- \
            benchmark/rotor_hover_solver_phase1_fgsprecond.jl; then
        status fgsprecond ok
    else
        status fgsprecond FAILED; ok=0
    fi
fi

if [ "$ok" = 1 ]; then
    if jrun p2tune MEM_BUDGETS=0:16:32:64:128:500 PHASE2_OUTDIR="$run/phase2" \
            TUNE_MAX_SECONDS=14400 -- \
            benchmark/rotor_hover_solver_phase2_tune.jl; then
        status p2tune ok
    else
        status p2tune FAILED; ok=0
    fi
fi

if [ "$ok" = 1 ]; then
    if jrun p2 MEM_BUDGETS=16:32:64:128:500 PHASE2_OUTDIR="$run/phase2" -- \
            benchmark/rotor_hover_solver_phase2.jl; then
        status p2 ok
    else
        status p2 FAILED; ok=0
    fi
fi

printf 'completed rung=%s j=%s ok=%s\n' "$RUNG_SEL" "$J" "$ok" > "$run/COMPLETED"
