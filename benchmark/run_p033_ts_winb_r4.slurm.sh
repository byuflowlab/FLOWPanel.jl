#!/usr/bin/env bash
#SBATCH --job-name=p033-ts-winb-r4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=48:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/slurm/p033-ts-winb-r4-%A_%a.out
#SBATCH --error=logs/slurm/p033-ts-winb-r4-%A_%a.err
# p033 R4 thread-scaling, Job B stage 2: warm-started rotor-hover legs across
# j = 1/8/16/32/64 under the 018-ported wake physics (WAKE_ENV_018=1), per the
# 021 winB protocol — restart from the FAMILY checkpoint at step 108 (rev 4,
# 36 steps). One array task per thread count; each task runs the proj2
# best-vs-best pair sequentially (fgs_proj2, then ilu_nfcache_proj2), one
# fresh Julia process per leg.
#
# Knobs are the FIXED certified R4 j64 set at every j (Ryan 2026-10-02 ruling:
# no re-tuning; recorded in every CSV row, labeled as transplants at harvest).
# Checkpoints come from run_p033_ts_ckpt_r4.slurm.sh (campaign-qualified names
# p033ts_R4_ckpt_{fgs,ilu}, read via the worktree's data/ symlink) — submit
# this job with --dependency=afterok:<ckpt job id>.
#
# Reporting caveat carried from the protocol: solver warm-start histories are
# not serialized in the checkpoint, so the first (order+1) restarted steps are
# effectively cold — flag in the harvest, never silently exclude.
#
# Required env on the submit line:
#   WSR4_PROJECT, CAMPAIGN_PINS, WSR4_DATA_ROOT, FGS_TOL_ABS,
#   KNOBS_P/KNOBS_MAC/KNOBS_LEAF
# Judge by STATUS_*/COMPLETED_* files + the CSV solved column, never sacct.
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${WSR4_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${WSR4_DATA_ROOT:?}" \
  "${FGS_TOL_ABS:?}" "${KNOBS_P:?}" "${KNOBS_MAC:?}" "${KNOBS_LEAF:?}"

LADDER=(1 8 16 32 64)
J="${LADDER[$SLURM_ARRAY_TASK_ID]}"

# ---- pinned-hardware assertions ----------------------------------------------
NODE_CPUS=$(getconf _NPROCESSORS_ONLN)
ALLOC_CPUS=${SLURM_CPUS_ON_NODE:-0}
if [ "$NODE_CPUS" != "128" ] || [ "$ALLOC_CPUS" != "128" ]; then
  echo "ERROR: need an exclusive 128-core zen3 node (got node=$NODE_CPUS" >&2
  echo "       alloc=$ALLOC_CPUS); timings would not be comparable." >&2
  exit 1
fi

# ---- NFS precompile-lock guard -----------------------------------------------
export JULIA_PKG_PRECOMPILE_AUTO=0
if command -v flock >/dev/null 2>&1; then
  ( if flock -E 99 -w 3600 9; then
      julia --project="$WSR4_PROJECT" --startup-file=no \
        -e 'using Pkg; Pkg.precompile()' \
        || echo "WARNING: Pkg.precompile() failed; compiling in memory"
    else
      echo "WARNING: precompile lock unavailable; compiling in memory"
    fi
  ) 9>"$HOME/.julia/flowpanel-021-precompile.lock"
else
  echo "WARNING: flock(1) not found; skipping shared precompile"
fi

run="$WSR4_DATA_ROOT/p033ts-winb-j$J-${RESUME_FROM_JOB_ID:-$SLURM_ARRAY_JOB_ID}"
if [ -n "${RESUME_FROM_JOB_ID:-}" ] && [ -d "$run" ]; then
    rm -f "$run/COMPLETED"
else
    mkdir "$run"
fi
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$WSR4_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
echo "julia_threads=$J knobs=fixed-transplant-j64" > "$run/point.txt"

ILV="--interleave=0-3 --cpunodebind=0-3"

export JULIA_NUM_THREADS="$J"
export BENCH_BLAS_THREADS=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
export EXPECT_JULIA_THREADS="$J" THREADING_MODE=multi
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3-blas1}"
export RUNG=R4 NT=36 PHASE=phase3wsr4 SNAPSHOT_STRENGTHS=1
export FGS_PRECISION="${FGS_PRECISION:-f64}"
export WAKE_ENV_018=1

# Das arc table + checkpoint reachability (relative data/ paths)
[ -e data ] || { echo "ERROR: no data/ symlink in $(pwd) — deploy step missing" >&2; exit 1; }
[ -f data/p018_cs_l3p4_rs1_te_downwash_te.csv ] || {
  echo "ERROR: Das arc table not reachable via data/ — jobs die ~1 min in" >&2
  exit 1; }
for fam in fgs ilu; do
  [ -d "data/p033ts_R4_ckpt_$fam" ] || {
    echo "ERROR: checkpoint data/p033ts_R4_ckpt_$fam missing — run" >&2
    echo "       run_p033_ts_ckpt_r4.slurm.sh first" >&2; exit 1; }
done

( while true; do
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] heartbeat: p033ts-winb j$J alive"
    sleep 300
  done ) &
HEARTBEAT_PID=$!
trap 'kill "$HEARTBEAT_PID" 2>/dev/null || true' EXIT

FAILED_COUNT=0
for arm in fgs_proj2 ilu_nfcache_proj2; do
  case "$arm" in fgs_*) fam=fgs ;; *) fam=ilu ;; esac
  armleg="${arm}_winB"
  if [ "$(cat "$run/STATUS_$armleg" 2>/dev/null || true)" = ok ]; then
    echo "resume: $armleg already ok — skipping"
    continue
  fi
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] === ARM=$arm LEG=winB (j=$J) ==="
  if env ARM="$arm" WSR4_LEG=winB OUTDIR_OVERRIDE="$run" \
        RESTART_NAME="p033ts_R4_ckpt_$fam" \
        RESTART_PATH="data/p033ts_R4_ckpt_$fam" \
        numactl $ILV \
        julia --project="$WSR4_PROJECT" --startup-file=no -t "$J" \
        benchmark/fgs_r4_warmstart_ab.jl > "$run/$armleg.log" 2>&1 \
      && [ -f "$run/COMPLETED_$armleg" ]; then
    echo "  $armleg ok"
  else
    FAILED_COUNT=$((FAILED_COUNT + 1))
    echo "WARNING: $armleg FAILED (continuing; see $run/$armleg.log)"
  fi
done

printf 'completed j=%s failed_count=%s\n' "$J" "$FAILED_COUNT" > "$run/COMPLETED"
