#!/usr/bin/env bash
#SBATCH --job-name=p021-fgs-stage2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/r4-fgs-stage2-%j.out
#SBATCH --error=logs/slurm/r4-fgs-stage2-%j.err
# 021 FGS scalability diagnostic — Stage 2 (mechanism): per-worker drain-loop
# aggregates + one bounded-backoff idle-policy A/B, on ONE exclusive zen3
# node. Stage 1 (job 13858983) reproduced both effects and localized them to
# the dagteam near-field sweep (nonself_product anti-scales 2.0→2.7→4.5 s
# over j=16/32/64); this stage separates the three candidate mechanisms:
#
#   * per-task busy time inflating with j at fixed task count  → memory
#     bandwidth / locality;
#   * flat busy but idle+lockmgmt shares growing with j        → scheduling
#     (narrow DAG width and/or queue-lock contention);
#   * a backoff recovery on top of that                        → the idle
#     lock hammering specifically (workers polling the SpinLock while empty).
#
# COST CUTS vs Stage 1 (per Ryan 2026-09-23): no j=1 rung, no placement A/B,
# no cap32, no accepted bridge, 1-2 blocks instead of 3, 12 h wall instead of
# 48. All arms run the champion placement (socket 0, interleave 0-3) and the
# per-thread calibrated configs from the verified 13777133 ladder.
#
# Stages (fixed-work arm, 27x3/tolerance=0, gates as Stage 1):
#   A aggdiag   — diag=1 (coarse phases + NEW per-worker aggregates):
#                 j16-b1, j32-b1, j64-b1, j64-b2 at w=0; cap-j64-b1/b2 at w=16
#   B anchors   — diag=0 spin: j16-b1, j32-b1 (w=0), cap-j64-b1 (w=16);
#                 the j64 w=0 anchors are section C's spin arms
#   C idle A/B  — 3 pairs spin-vs-backoff @ j=64 w=0, alternating arm order
#   D safety    — backoff @ j16 w=0 and @ j64 w=16 (must not regress)
#   E confirm   — 2 pairs spin-vs-backoff @ j=32 ONLY if backoff's paired
#                 median beats spin at j=64 (backoff_verdict.txt)
#
# Judge by outputs (STATUS_* files), never sacct; logs are output-buffered.
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"

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
run="$COLD_DATA_ROOT/fgs-stage2-${RESUME_FROM_JOB_ID:-$SLURM_JOB_ID}"
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

# Champion placement only (Stage 1 settled placement decisively: both-socket
# was +4.4 s at j=64)
ILV="--interleave=0-3 --cpunodebind=0-3"
numactl $ILV numactl --show > "$run/numactl_show_champion.txt" 2>&1 || true

export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3-blas1}"
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export JULIA_DEBUG="${JULIA_DEBUG:-loading}"
export RUNG=R4 CONFIGS=fgs STAGE=verify COLD_PREPARED_ONLY=1
export CONFIG_FILE="$PWD/benchmark/retained_r4_diagnostics.toml"
export STAGE1_SOLVES="${STAGE1_SOLVES:-5}"
CALIB_ROOT="${CALIB_ROOT:-$COLD_DATA_ROOT}"
CALIB_JOB="${CALIB_JOB:-13777133}"
dagteam_config() { # $1 = j
  echo "$CALIB_ROOT/thread-scaling-j$1-$CALIB_JOB/fgs-calibrate/results/dagteam_selected.toml"
}
for J in 16 32 64; do
  [ -f "$(dagteam_config "$J")" ] || { echo "ERROR: missing $(dagteam_config "$J")" >&2; exit 1; }
done

( while true; do
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] heartbeat: stage2 alive"
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

# one_process NAME J DIAG BLOCK LABEL WORKERS IDLE
# (champion placement + fixed arm + per-j calibrated config are implicit)
# Failures are findings: record STATUS_<name>=FAILED and continue.
FAILED_COUNT=0
one_process() {
  local name=$1 j=$2 diag=$3 block=$4 label=$5 workers=$6 idle=$7
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
  if env STAGE1_ARM=fixed STAGE1_DIAG="$diag" STAGE1_BLOCK="$block" \
        STAGE1_LABEL="$label" DAGTEAM_WORKERS="$workers" DAGTEAM_IDLE="$idle" \
        DAGTEAM_CONFIG="$(dagteam_config "$j")" \
        OUTDIR="$dir/results" BENCH_CASE_ROOT="$dir/fixture" \
        numactl $ILV bash benchmark/run_cold_process.sh "$j" 1 \
        benchmark/fgs_r4_dagteam_stage1.jl > "$dir/process.log" 2>&1 \
      && [ -f "$dir/results/status.toml" ]; then
    status "$name" ok
  else
    status "$name" FAILED
    FAILED_COUNT=$((FAILED_COUNT + 1))
    echo "WARNING: $name FAILED (continuing; see $dir/process.log)"
  fi
}

# ---- A: per-worker aggregates (diag=1), spin --------------------------------
for spec in $(shuf -e j16-b1 j32-b1 j64-b1 j64-b2); do
  J=${spec%%-b*}; J=${J#j}; B=${spec##*-b}
  one_process "aggdiag-${spec}" "$J" 1 "$B" aggdiag 0 spin
done
for B in 1 2; do
  one_process "aggdiag-cap-j64-b${B}" 64 1 "$B" aggdiag-cap 16 spin
done

# ---- B: uninstrumented spin anchors (overhead gate; j64 w=0 anchors are C's
# spin arms) -------------------------------------------------------------------
for J in $(shuf -e 16 32); do
  one_process "anchor-j${J}-b1" "$J" 0 1 anchor 0 spin
done
one_process "anchor-cap-j64-b1" 64 0 1 anchor-cap 16 spin

# ---- C: idle-policy A/B @ j=64 w=0 (paired, alternating order) ---------------
run_idle_pairs() { # $1 = j, $2 = npairs
  local j=$1 npairs=$2 pair
  for pair in $(seq 1 "$npairs"); do
    if [ $((pair % 2)) = 1 ]; then order="spin backoff"; else order="backoff spin"; fi
    for idle in $order; do
      one_process "idle-${idle}-j${j}-p${pair}" "$j" 0 "$pair" "idle-${idle}" 0 "$idle"
    done
  done
}
run_idle_pairs 64 3

# ---- D: backoff safety checks ------------------------------------------------
one_process "idle-backoff-j16-b1" 16 0 1 idle-backoff-safety 0 backoff
one_process "idle-backoff-cap-j64-b1" 64 0 1 idle-backoff-cap 16 backoff

# ---- E: confirm at 32 only if backoff helps at 64 ----------------------------
backoff_helps=$(julia --startup-file=no -e '
  using TOML, Statistics
  run = ARGS[1]
  med(idle) = begin
      t = Float64[]
      for pair in 1:3
          f = joinpath(run, "idle-$(idle)-j64-p$(pair)", "results", "stage1_summary.toml")
          isfile(f) && push!(t, TOML.parsefile(f)["median_solve_seconds"])
      end
      isempty(t) ? Inf : median(t)
  end
  println(med("backoff") < med("spin") ? "yes" : "no")' "$run" 2>/dev/null || echo no)
echo "backoff_helps=$backoff_helps" | tee "$run/backoff_verdict.txt"
if [ "$backoff_helps" = yes ]; then
  run_idle_pairs 32 2
else
  echo "j=32 idle pairs skipped (backoff did not beat spin at j=64)" >> "$run/backoff_verdict.txt"
fi

printf 'completed failed_count=%s\n' "$FAILED_COUNT" > "$run/COMPLETED"
