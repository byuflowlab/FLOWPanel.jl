#!/usr/bin/env bash
#SBATCH --job-name=p033-ts-ckpt-r4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/slurm/p033-ts-ckpt-r4-%A_%a.out
#SBATCH --error=logs/slurm/p033-ts-ckpt-r4-%A_%a.err
# p033 R4 thread-scaling, Job B stage 1: fresh family checkpoints under the
# 018-ported wake physics (WAKE_ENV_018=1). The ported physics changes the
# case, so the 2026-09-24 fgs_wsr4_R4_ckpt_* checkpoints are INVALID here —
# each family's cold arm re-marches revs 1-3 (108 steps, VTK on) at j=64.
# Checkpoints are physics-level and thread-independent: the winB job restarts
# from these at every j. RUN_NAME is campaign-qualified (p033ts_...) so the
# old checkpoints are never read or overwritten.
#
# Array task 0 = fgs family (ARM=fgs_cold), 1 = ilu family
# (ARM=ilu_nfcache_cold). One fresh Julia process per leg, champion placement
# (socket 0, interleave 0-3), BLAS pinned to 1 for both families.
#
# Required env on the submit line:
#   WSR4_PROJECT, CAMPAIGN_PINS, WSR4_DATA_ROOT, FGS_TOL_ABS,
#   KNOBS_P/KNOBS_MAC/KNOBS_LEAF (certified R4 budget-500 apply knobs)
# Judge by STATUS_*/COMPLETED_* files + the CSV solved column, never sacct.
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${WSR4_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${WSR4_DATA_ROOT:?}" \
  "${FGS_TOL_ABS:?}" "${KNOBS_P:?}" "${KNOBS_MAC:?}" "${KNOBS_LEAF:?}"

FAMILIES=(fgs ilu)
FAMILY="${FAMILIES[$SLURM_ARRAY_TASK_ID]}"
case "$FAMILY" in
  fgs) ARM=fgs_cold ;;
  ilu) ARM=ilu_nfcache_cold ;;
esac

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

run="$WSR4_DATA_ROOT/p033ts-ckpt-$FAMILY-${RESUME_FROM_JOB_ID:-$SLURM_ARRAY_JOB_ID}"
if [ -n "${RESUME_FROM_JOB_ID:-}" ] && [ -d "$run" ]; then
    rm -f "$run/COMPLETED"
else
    mkdir "$run"
fi
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$WSR4_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"

ILV="--interleave=0-3 --cpunodebind=0-3"

THREADS=64
export JULIA_NUM_THREADS="$THREADS"
export BENCH_BLAS_THREADS=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
export EXPECT_JULIA_THREADS="$THREADS" THREADING_MODE=multi
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3-blas1}"
export RUNG=R4 NT=36 PHASE=phase3wsr4 SNAPSHOT_STRENGTHS=1
export FGS_PRECISION="${FGS_PRECISION:-f64}"
export WAKE_ENV_018=1
# Campaign-qualified checkpoint name: never read/write the 2026-09-24
# checkpoints (different physics). The winB launcher must use the same name.
export RUN_NAME="p033ts_R4_ckpt_$FAMILY"

( while true; do
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] heartbeat: p033ts-ckpt-$FAMILY alive"
    sleep 300
  done ) &
HEARTBEAT_PID=$!
trap 'kill "$HEARTBEAT_PID" 2>/dev/null || true' EXIT

# RHPC resolves inputs/outputs via relative data/ paths
[ -e data ] || { echo "ERROR: no data/ symlink in $(pwd) — deploy step missing" >&2; exit 1; }
[ -f data/p018_cs_l3p4_rs1_te_downwash_te.csv ] || {
  echo "ERROR: Das arc table not reachable via data/ — jobs die ~1 min in" >&2
  exit 1; }

armleg="${ARM}_ckpt"
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] === ARM=$ARM LEG=ckpt (j=$THREADS) ==="
if env ARM="$ARM" WSR4_LEG=ckpt OUTDIR_OVERRIDE="$run" \
      numactl $ILV \
      julia --project="$WSR4_PROJECT" --startup-file=no -t "$THREADS" \
      benchmark/fgs_r4_warmstart_ab.jl > "$run/$armleg.log" 2>&1 \
    && [ -f "$run/COMPLETED_$armleg" ]; then
  echo "  $armleg ok"
  printf 'completed family=%s\n' "$FAMILY" > "$run/COMPLETED"
else
  echo "ERROR: $armleg FAILED (see $run/$armleg.log)" >&2
  exit 1
fi
