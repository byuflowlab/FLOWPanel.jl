#!/usr/bin/env bash
#SBATCH --job-name=p021-fgs-dagedge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/r4-fgs-dagedge-%j.out
#SBATCH --error=logs/slurm/r4-fgs-dagedge-%j.err
# 021 :dagedge HPC benchmark (fgs_dagedge_benchmark_reset_prompt_20260924.md):
# edge-level partial pulls on a static schedule vs the dagteam+backoff Stage-2
# champion (3.24 s/solve at j=64). Gate-0 predicts a 4.25x sweep bound
# (289.5 -> 68-69 MB/sweep at theta=4KB), roughly flat in j; expected
# end-to-end ~1.7-1.8 s. ONE exclusive zen3 node, champion placement
# (socket 0, interleave 0-3), per-j calibrated configs from the verified
# 13777133 ladder, fixed 27x3 cold-START work (cold = zero-initial-guess
# solves; arms BATCHED per (j, block) process per Ryan 2026-09-23).
#
# RUN_MODE=perf (default) — submission 1, uninstrumented performance trials:
#   * primary-pN   : 3 paired A/B blocks @ j=64 backoff, dagteam vs
#                    dagedge(theta=4KB), arm order alternating across blocks
#   * ladder-jJ    : both executors @ j in {16,32} backoff (j=64 = primary)
#   * theta-j64    : dagedge theta in {0,4KB,16KB} @ j=64 backoff, one process
# RUN_MODE=profile — submission 2, diagnostics-instrumented attribution
#   (NEVER a performance trial; rankings come only from the perf run):
#   * prof-j64-b1/b2 : primary A/B pair instrumented, 2 blocks
#   * prof-jJ        : instrumented pair @ j in {16,32}
#
# Judge by outputs (STATUS_* files), never sacct; logs are output-buffered.
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"
RUN_MODE="${RUN_MODE:-perf}"
case "$RUN_MODE" in perf|profile) ;; *) echo "ERROR: RUN_MODE must be perf or profile" >&2; exit 1;; esac

# ---- pinned-hardware assertions (021 ruling 2026-08-24) ----------------------
NODE_CPUS=$(getconf _NPROCESSORS_ONLN)
ALLOC_CPUS=${SLURM_CPUS_ON_NODE:-0}
if [ "$NODE_CPUS" != "128" ] || [ "$ALLOC_CPUS" != "128" ]; then
  echo "ERROR: need an exclusive 128-core zen3 node (got node=$NODE_CPUS" >&2
  echo "       alloc=$ALLOC_CPUS); timings would not be comparable." >&2
  exit 1
fi
NUMA_NODES=$(numactl --hardware | awk '/^available:/{print $2}')
if [ "$NUMA_NODES" != "8" ]; then
  echo "ERROR: expected 8 NUMA nodes (NPS4 dual zen3), got $NUMA_NODES" >&2
  exit 1
fi

# ---- NFS precompile-lock guard (2026-08-25 stall) ----------------------------
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

# Resume (RESUME_FROM_JOB_ID=<old job id>): reuse that run dir; STATUS_*=ok
# stages skip, failed/missing stages re-run with their stale outputs moved
# aside by one_process.
run="$COLD_DATA_ROOT/fgs-dagedge-$RUN_MODE-${RESUME_FROM_JOB_ID:-$SLURM_JOB_ID}"
if [ -n "${RESUME_FROM_JOB_ID:-}" ] && [ -d "$run" ]; then
    prev="$run/logs.before.$SLURM_JOB_ID"
    mkdir -p "$prev"
    mv "$run/campaign_pins.toml" "$run/Manifest.toml" "$prev"/ 2>/dev/null || true
    rm -f "$run/COMPLETED"
else
    mkdir "$run"
fi
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
numactl --hardware > "$run/numactl_hardware.txt" 2>&1 || true

# Champion placement only (Stage 1 settled placement decisively)
ILV="--interleave=0-3 --cpunodebind=0-3"
numactl $ILV numactl --show > "$run/numactl_show_champion.txt" 2>&1 || true

export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3-blas1}"
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export JULIA_DEBUG="${JULIA_DEBUG:-loading}"
export RUNG=R4 CONFIGS=fgs STAGE=verify COLD_PREPARED_ONLY=1
export CONFIG_FILE="$PWD/benchmark/retained_r4_diagnostics.toml"
export DAGEDGE_SOLVES="${DAGEDGE_SOLVES:-5}"
CALIB_ROOT="${CALIB_ROOT:-$COLD_DATA_ROOT}"
CALIB_JOB="${CALIB_JOB:-13777133}"
dagteam_config() { # $1 = j
  echo "$CALIB_ROOT/thread-scaling-j$1-$CALIB_JOB/fgs-calibrate/results/dagteam_selected.toml"
}
for J in 16 32 64; do
  [ -f "$(dagteam_config "$J")" ] || { echo "ERROR: missing $(dagteam_config "$J")" >&2; exit 1; }
done

( while true; do
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] heartbeat: dagedge-$RUN_MODE alive"
    sleep 300
  done ) &
HEARTBEAT_PID=$!
( while true; do
    awk -v t="$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
      '/^cpu MHz/{n++; s+=$4; if($4>mx)mx=$4; if(mn==0||$4<mn)mn=$4}
       END{printf "%s,%d,%.0f,%.0f,%.0f\n", t, n, mn, s/n, mx}' \
      /proc/cpuinfo >> "$run/cpufreq.csv"
    sleep 60
  done ) &
CPUFREQ_PID=$!
trap 'kill "$HEARTBEAT_PID" "$CPUFREQ_PID" 2>/dev/null || true' EXIT

status() { printf '%s\n' "$2" > "$run/STATUS_$1"; }
stage_ok() { [ "$(cat "$run/STATUS_$1" 2>/dev/null || true)" = ok ]; }

# ---- shared parse + precompile ----------------------------------------------
export OUTDIR="$run/parse" BENCH_CASE_ROOT="$run/fixture-controls"
bash benchmark/run_cold_process.sh 1 1 \
    benchmark/cold_parse.jl > "$run/parse.log" 2>&1
bash benchmark/run_cold_process.sh 4 1 \
    benchmark/cold_precompile.jl > "$run/precompile.log" 2>&1

# one_process NAME J DIAG BLOCK LABEL ARMS
# (champion placement + backoff + per-j calibrated config are implicit;
#  DAGTEAM_WORKERS=0 throughout — Stage 2's champion is uncapped backoff)
# Failures are findings: record STATUS_<name>=FAILED and continue.
FAILED_COUNT=0
one_process() {
  local name=$1 j=$2 diag=$3 block=$4 label=$5 arms=$6
  if stage_ok "$name"; then
    echo "resume: $name already ok — skipping"
    return 0
  fi
  local dir="$run/$name"
  mkdir -p "$dir"
  for stale in results fixture; do
    if [ -e "$dir/$stale" ]; then
      mv "$dir/$stale" "$dir/$stale.before.$SLURM_JOB_ID"
    fi
  done
  if env DAGEDGE_DIAG="$diag" DAGEDGE_BLOCK="$block" \
        DAGEDGE_LABEL="$label" DAGEDGE_ARMS="$arms" \
        DAGTEAM_WORKERS=0 DAGTEAM_IDLE=backoff \
        DAGTEAM_CONFIG="$(dagteam_config "$j")" \
        OUTDIR="$dir/results" BENCH_CASE_ROOT="$dir/fixture" \
        numactl $ILV bash benchmark/run_cold_process.sh "$j" 1 \
        benchmark/fgs_r4_dagedge.jl > "$dir/process.log" 2>&1 \
      && [ -f "$dir/results/status.toml" ]; then
    status "$name" ok
  else
    status "$name" FAILED
    FAILED_COUNT=$((FAILED_COUNT + 1))
    echo "WARNING: $name FAILED (continuing; see $dir/process.log)"
  fi
}

AB="dagteam,dagedge:4096"
BA="dagedge:4096,dagteam"

if [ "$RUN_MODE" = perf ]; then
  # ---- primary A/B @ j=64 (3 paired blocks, alternating in-process order) ----
  for P in 1 2 3; do
    if [ $((P % 2)) = 1 ]; then arms="$AB"; else arms="$BA"; fi
    one_process "primary-p${P}" 64 0 "$P" primary "$arms"
  done
  # ---- flat-in-j ladder @ j in {16,32} ---------------------------------------
  for J in $(shuf -e 16 32); do
    one_process "ladder-j${J}" "$J" 0 1 ladder "$AB"
  done
  # ---- theta probe @ j=64 (one process; theta is a plan-build knob) ----------
  one_process "theta-j64" 64 0 1 theta "dagedge:0,dagedge:4096,dagedge:16384"
else
  # ---- instrumented attribution (never pooled with perf trials) --------------
  for B in 1 2; do
    if [ $((B % 2)) = 1 ]; then arms="$AB"; else arms="$BA"; fi
    one_process "prof-j64-b${B}" 64 1 "$B" profile "$arms"
  done
  for J in $(shuf -e 16 32); do
    one_process "prof-j${J}" "$J" 1 1 profile "$AB"
  done
fi

printf 'completed failed_count=%s\n' "$FAILED_COUNT" > "$run/COMPLETED"
