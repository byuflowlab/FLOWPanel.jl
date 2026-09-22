#!/usr/bin/env bash
#SBATCH --job-name=p021-fgs-stage1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=48:00:00
#SBATCH --output=logs/slurm/r4-fgs-stage1-%j.out
#SBATCH --error=logs/slurm/r4-fgs-stage1-%j.err
# 021 FGS scalability diagnostic — Stage 1 (plan C
# fgs_scalability_diagnostic_plan_20260921c.md): reproduce and decompose the
# dagteam 16->32 plateau and the 64-thread regression on ONE exclusive zen3
# node. Everything runs from benchmark/fgs_r4_dagteam_stage1.jl in fresh
# processes; this launcher sequences (per-thread calibrated configs come from
# the verified 13777133 ladder):
#
#   1. fixed-work ladder  — 27x3/tolerance=0 at j in {1,16,32,64}, 3 fresh-
#      process blocks x STAGE1_SOLVES warmed solves, thread order shuffled
#      within each block; uninstrumented (champion placement, BLAS=1).
#   2. diagnostics ladder — separate matched runs with the coarse-phase
#      diagnostics dict (incl. new dagteam spawn/join/wait/reduce splits);
#      analysis must verify <=5% overhead and unchanged scaling shape before
#      attributing anything.
#   3. accepted-accuracy bridge — calibrated production config at 16/32/64,
#      1 block (Stage-0's verified 40-row/j evidence carries the screen; this
#      ties the new pins to it).
#   4. placement A/B @ j=64 — champion (interleave+cpunodebind 0-3 = socket 0)
#      vs both-socket (0-7), fresh first-touch under each, 3 pairs,
#      alternating order.
#   5. worker-cap A/B @ j=64 — DAGTEAM_WORKERS=16 vs 0 on the champion
#      placement, 3 pairs, alternating; cap=32 pairs run only if cap=16's
#      median beats the paired baseline.
#
# Plateau and regression stay SEPARATE conclusions; decomposition happens in
# the analysis script, not here. Judge by outputs (STATUS_* files), never
# sacct; logs are output-buffered.
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
# both placements below assume NPS4 dual-socket zen3 = NUMA nodes 0-7
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
# aside by one_process. First-pass top-level logs go to logs.before.<new id>/.
run="$COLD_DATA_ROOT/fgs-stage1-${RESUME_FROM_JOB_ID:-$SLURM_JOB_ID}"
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

# Champion placement (verified 13777133 ladder): socket 0 only, memory
# interleaved over its 4 NUMA nodes. ALT spreads the same 64 threads over both
# sockets/8 nodes (H4/H6 discriminator). Fresh first-touch is automatic: every
# arm is a fresh process constructing its own state under its policy.
ILV="--interleave=0-3 --cpunodebind=0-3"
ALT="--interleave=0-7 --cpunodebind=0-7"
numactl $ILV numactl --show > "$run/numactl_show_champion.txt" 2>&1 || true
numactl $ALT numactl --show > "$run/numactl_show_alt.txt" 2>&1 || true

export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3-blas1}"
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export JULIA_DEBUG="${JULIA_DEBUG:-loading}"
export RUNG=R4 CONFIGS=fgs STAGE=verify COLD_PREPARED_ONLY=1
export CONFIG_FILE="$PWD/benchmark/retained_r4_diagnostics.toml"
export STAGE1_SOLVES="${STAGE1_SOLVES:-5}"
# per-thread calibrated dagteam configs from the verified thread-scaling ladder
CALIB_ROOT="${CALIB_ROOT:-$COLD_DATA_ROOT}"
CALIB_JOB="${CALIB_JOB:-13777133}"
dagteam_config() { # $1 = j
  echo "$CALIB_ROOT/thread-scaling-j$1-$CALIB_JOB/fgs-calibrate/results/dagteam_selected.toml"
}
for J in 1 16 32 64; do
  [ -f "$(dagteam_config "$J")" ] || { echo "ERROR: missing $(dagteam_config "$J")" >&2; exit 1; }
done

( while true; do
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] heartbeat: stage1 alive"
    sleep 300
  done ) &
HEARTBEAT_PID=$!
# under-load core-frequency observation (measurement contract; from outside
# the Julia processes so timing is never perturbed)
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

# one_process NAME PLACEMENT J ARM DIAG BLOCK LABEL WORKERS CONFIG
# Failures are findings: record STATUS_<name>=FAILED and continue.
FAILED_COUNT=0
one_process() {
  local name=$1 placement=$2 j=$3 arm=$4 diag=$5 block=$6 label=$7 workers=$8 config=$9
  if stage_ok "$name"; then
    echo "resume: $name already ok — skipping"
    return 0
  fi
  local dir="$run/$name"
  mkdir -p "$dir"
  # cold_preflight requires a NEW OUTDIR and a fresh BENCH_CASE_ROOT: move any
  # stale first-pass outputs of a non-ok stage aside before re-running
  for stale in results fixture; do
    if [ -e "$dir/$stale" ]; then
      mv "$dir/$stale" "$dir/$stale.before.$SLURM_JOB_ID"
    fi
  done
  if env STAGE1_ARM="$arm" STAGE1_DIAG="$diag" STAGE1_BLOCK="$block" \
        STAGE1_LABEL="$label" DAGTEAM_WORKERS="$workers" \
        DAGTEAM_CONFIG="$config" \
        OUTDIR="$dir/results" BENCH_CASE_ROOT="$dir/fixture" \
        numactl $placement bash benchmark/run_cold_process.sh "$j" 1 \
        benchmark/fgs_r4_dagteam_stage1.jl > "$dir/process.log" 2>&1 \
      && [ -f "$dir/results/status.toml" ]; then
    status "$name" ok
  else
    status "$name" FAILED
    FAILED_COUNT=$((FAILED_COUNT + 1))
    echo "WARNING: $name FAILED (continuing; see $dir/process.log)"
  fi
}

# ---- 1+2: fixed-work ladder, uninstrumented then diagnostics ----------------
for diag in 0 1; do
  tag=$([ "$diag" = 0 ] && echo ladder || echo ladderdiag)
  for block in 1 2 3; do
    for J in $(shuf -e 1 16 32 64); do
      one_process "${tag}-j${J}-b${block}" "$ILV" "$J" fixed "$diag" "$block" \
        "$tag" 0 "$(dagteam_config "$J")"
    done
  done
done

# ---- 3: accepted-accuracy bridge (new pins vs Stage-0 verified ladder) ------
for J in $(shuf -e 16 32 64); do
  one_process "accepted-j${J}-b1" "$ILV" "$J" accepted 0 1 accepted 0 \
    "$(dagteam_config "$J")"
done

# ---- 4: placement A/B @ j=64 (fixed work, alternating order) ----------------
CONF64="$(dagteam_config 64)"
for pair in 1 2 3; do
  if [ $((pair % 2)) = 1 ]; then order="champion alt"; else order="alt champion"; fi
  for armname in $order; do
    placement=$([ "$armname" = champion ] && echo "$ILV" || echo "$ALT")
    one_process "placement-${armname}-j64-p${pair}" "$placement" 64 fixed 0 \
      "$pair" "placement-$armname" 0 "$CONF64"
  done
done

# ---- 5: worker-cap A/B @ j=64 (champion placement, alternating order) -------
run_cap_pairs() { # $1 = cap
  local cap=$1 pair
  for pair in 1 2 3; do
    if [ $((pair % 2)) = 1 ]; then order="$cap 0"; else order="0 $cap"; fi
    for w in $order; do
      one_process "cap${cap}-w${w}-j64-p${pair}" "$ILV" 64 fixed 0 "$pair" \
        "cap${cap}-arm" "$w" "$CONF64"
    done
  done
}
run_cap_pairs 16

# cap=32 only if cap=16 helps: compare paired medians from the landed rows
cap16_helps=$(julia --startup-file=no -e '
  using TOML, Statistics
  run = ARGS[1]
  med(pat, w) = begin
      t = Float64[]
      for pair in 1:3
          f = joinpath(run, "cap16-w$(w)-j64-p$(pair)", "results", "stage1_summary.toml")
          isfile(f) && push!(t, TOML.parsefile(f)["median_solve_seconds"])
      end
      isempty(t) ? Inf : median(t)
  end
  println(med("cap16", 16) < med("cap16", 0) ? "yes" : "no")' "$run" 2>/dev/null || echo no)
echo "cap16_helps=$cap16_helps" | tee "$run/cap16_verdict.txt"
if [ "$cap16_helps" = yes ]; then
  run_cap_pairs 32
else
  echo "cap=32 pairs skipped (cap=16 did not beat baseline)" >> "$run/cap16_verdict.txt"
fi

printf 'completed failed_count=%s\n' "$FAILED_COUNT" > "$run/COMPLETED"
