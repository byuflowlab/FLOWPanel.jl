#!/usr/bin/env bash
#SBATCH --job-name=p021-r4-thread-scaling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=48:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/slurm/r4-thread-scaling-%A_%a.out
#SBATCH --error=logs/slurm/r4-thread-scaling-%A_%a.err
# 021 preliminary R4 thread-scaling study (Ryan 2026-09-19): before adopting
# the full memory ladder (0/16/32/64/128/500 GiB with per-class thread caps,
# ~2 GiB/thread) as general benchmark policy, measure thread efficiency of the
# two production solver families at R4 across j = 1/8/16/32/64 so plateaued
# tiers can be pruned. One array task per thread count; each point is RE-TUNED
# at its own j (tolerances and knobs never carry across thread counts — 021
# v23 trap list):
#
#   FGS arm    — retained champion knobs (P8/MAC0.4/leaf100/inner3), dagteam
#                f32full + colored twins staircase-calibrated AT THIS j
#                (AB_MODE=calibrate = the numerical gate: every accepted solve
#                is evaluator-certified BC rel-L2 <= 1e-6), then uninstrumented
#                A/B trials. The colored arm doubles as the fallback-ladder
#                scaling curve for free.
#   iLU-GMRES  — FMM apply knobs re-descended at this j by the phase-2 tuner
#                (necessary: FastMultipole 4c0f1b8f -> f4d6b671 changed the
#                apply path — 030 block assembly, one-division kernel,
#                parallel NF cache build — so the 2026-08-25 knobs are stale)
#                at budgets 0 (uncached endpoint) and 500 (node), then
#                measured as krylov_ilu + krylov_ilu_nfcache by
#                rotor_hover_solver_phase2.jl at the fresh knobs.
#
# Placement is the promoted champion's and is load-bearing (socket-membind
# collapsed dagteam to 1.09x): interleave + cpunodebind over NUMA nodes 0-3 =
# all of socket 0 on NPS4 zen3, which holds all 64 cores the j=64 point needs.
# BLAS is pinned to 1 in every stage (champion convention). HARDWARE_TAG
# carries the placement so these rows/traces can never replay into or be read
# as the historical both-socket phase-2 rows.
#
# Each array task writes to its OWN directory under COLD_DATA_ROOT
# (PHASE2_OUTDIR + per-task BENCH_CASE_ROOT): five tasks appending to one NFS
# tune_phase2.csv is the append hazard that destroyed R1's rows on 2026-08-18.
#
# The FGS arm at j=1 exercises the dagteam executor single-threaded, which no
# gate has covered; its failure is a FINDING, not a job failure — each arm
# records its own STATUS_* verdict and the task continues (judge by outputs,
# never sacct).
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"

LADDER=(1 8 16 32 64)
J="${LADDER[$SLURM_ARRAY_TASK_ID]}"

# ---- pinned-hardware assertions (021 ruling 2026-08-24, from p2_tune.sh) ----
NODE_CPUS=$(getconf _NPROCESSORS_ONLN)
ALLOC_CPUS=${SLURM_CPUS_ON_NODE:-0}
if [ "$NODE_CPUS" != "128" ] || [ "$ALLOC_CPUS" != "128" ]; then
  echo "ERROR: need an exclusive 128-core zen3 node (got node=$NODE_CPUS" >&2
  echo "       alloc=$ALLOC_CPUS); timings would not be comparable." >&2
  exit 1
fi

# ---- NFS precompile-lock guard (2026-08-25 stall, from p2_tune.sh) ----------
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

run="$COLD_DATA_ROOT/thread-scaling-j$J-$SLURM_ARRAY_JOB_ID"
mkdir "$run"
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
numactl --hardware > "$run/numactl_hardware.txt" 2>&1 || true
echo "julia_threads=$J" > "$run/point.txt"

ILV="--interleave=0-3 --cpunodebind=0-3"
numactl $ILV numactl --show > "$run/numactl_show.txt" 2>&1 || true

# Placement-qualified provenance tag: part of the tuning trace's HARD guard,
# so a socket-0-interleave timing can never replay into a both-socket descent.
export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3-blas1}"
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export JULIA_DEBUG="${JULIA_DEBUG:-loading}"

( while true; do
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] heartbeat: j$J alive"
    sleep 300
  done ) &
HEARTBEAT_PID=$!
trap 'kill "$HEARTBEAT_PID" 2>/dev/null || true' EXIT

status() { printf '%s\n' "$2" > "$run/STATUS_$1"; }

# ---- shared cold-harness parse + precompile (thread-count independent) ------
export RUNG=R4 CONFIGS=fgs STAGE=verify COLD_PREPARED_ONLY=1
export CONFIG_FILE="$PWD/benchmark/retained_r4_diagnostics.toml"
export DAGTEAM_PRECISION="${DAGTEAM_PRECISION:-f32full}"

export OUTDIR="$run/parse" BENCH_CASE_ROOT="$run/fixture-controls"
bash benchmark/run_cold_process.sh 1 1 \
    benchmark/cold_parse.jl > "$run/parse.log" 2>&1
bash benchmark/run_cold_process.sh 4 1 \
    benchmark/cold_precompile.jl > "$run/precompile.log" 2>&1

# ---- FGS arm: calibrate at this j, then A/B trials --------------------------
fgs_ok=1
export AB_MODE=calibrate
export OUTDIR="$run/fgs-calibrate/results" BENCH_CASE_ROOT="$run/fixture-calibrate"
mkdir -p "$run/fgs-calibrate"
if numactl $ILV bash benchmark/run_cold_process.sh "$J" 1 \
        benchmark/fgs_r4_dagteam_ab.jl > "$run/fgs-calibrate/process.log" 2>&1 \
        && [ -f "$run/fgs-calibrate/results/dagteam_selected.toml" ] \
        && [ -f "$run/fgs-calibrate/results/colored_selected.toml" ]; then
    status fgs_calibrate ok
else
    status fgs_calibrate FAILED
    fgs_ok=0
fi

if [ "$fgs_ok" = 1 ]; then
    export AB_MODE=trials
    export COLORED_CONFIG="$run/fgs-calibrate/results/colored_selected.toml"
    export DAGTEAM_CONFIG="$run/fgs-calibrate/results/dagteam_selected.toml"
    export OUTDIR="$run/fgs-trials/results" BENCH_CASE_ROOT="$run/fixture-trials"
    mkdir -p "$run/fgs-trials"
    if numactl $ILV bash benchmark/run_cold_process.sh "$J" 1 \
            benchmark/fgs_r4_dagteam_ab.jl > "$run/fgs-trials/process.log" 2>&1 \
            && [ -f "$run/fgs-trials/results/ab_summary.toml" ]; then
        status fgs_trials ok
    else
        status fgs_trials FAILED
        fgs_ok=0
    fi
fi
unset AB_MODE COLORED_CONFIG DAGTEAM_CONFIG COLD_PREPARED_ONLY CONFIG_FILE
unset CONFIGS STAGE OUTDIR BENCH_CASE_ROOT

# ---- iLU-GMRES arm: re-descend apply knobs at this j, then measure ----------
# PHASE2_OUTDIR isolates this task's tune_phase2.csv + traces (the measurement
# stage reads its knobs from the same directory); BENCH_CASE_ROOT isolates the
# frozen-b cache. Budget 0 first (phase2.jl requires the budget-0 row for the
# uncached knobs). 10 h descent backstop per budget: 2 budgets + FGS arm +
# measurement must fit the 48 h task.
ilu_ok=1
mkdir -p "$run/ilu"
# OPENBLAS/OMP at load time as well as BENCH_BLAS_THREADS: on some BLAS
# builds set_num_threads alone fails the common.jl pin assert (seen on the
# local pre-submit smoke).
ilu_env=(RUNG=R4 MEM_BUDGETS=0:500 EXPECT_JULIA_THREADS="$J"
         THREADING_MODE=multi BENCH_BLAS_THREADS=1
         OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 CACHE_B=1
         PHASE2_OUTDIR="$run/ilu" BENCH_CASE_ROOT="$run/ilu-case"
         TUNE_SEED=15:0.55:32 TUNE_SEED_B0=10:0.6:6
         TUNE_MAX_SECONDS=36000)
if env "${ilu_env[@]}" numactl $ILV \
        julia --project="$COLD_PROJECT" --startup-file=no \
        --compiled-modules=existing -t "$J" \
        benchmark/rotor_hover_solver_phase2_tune.jl \
        > "$run/ilu-tune.log" 2>&1; then
    status ilu_tune ok
else
    status ilu_tune FAILED
    ilu_ok=0
fi

if [ "$ilu_ok" = 1 ]; then
    if env RUNG=R4 CONFIGS=krylov_ilu,krylov_ilu_nfcache MEM_BUDGETS=500 \
            EXPECT_JULIA_THREADS="$J" THREADING_MODE=multi \
            BENCH_BLAS_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
            CACHE_B=1 PHASE2_OUTDIR="$run/ilu" \
            BENCH_CASE_ROOT="$run/ilu-case" numactl $ILV \
            julia --project="$COLD_PROJECT" --startup-file=no \
            --compiled-modules=existing -t "$J" \
            benchmark/rotor_hover_solver_phase2.jl \
            > "$run/ilu-measure.log" 2>&1; then
        status ilu_measure ok
    else
        status ilu_measure FAILED
        ilu_ok=0
    fi
fi

# COMPLETED = the task ran to the end; per-arm verdicts live in STATUS_* files.
# A failed arm at one j is a scaling finding, not a lost task.
printf 'completed j=%s fgs_ok=%s ilu_ok=%s\n' "$J" "$fgs_ok" "$ilu_ok" \
    > "$run/COMPLETED"
