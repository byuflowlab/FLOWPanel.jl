#!/usr/bin/env bash
#SBATCH --job-name=p021-fgs-wsr4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=36:00:00
#SBATCH --output=logs/slurm/r4-fgs-wsr4-%j.out
#SBATCH --error=logs/slurm/r4-fgs-wsr4-%j.err
# 021 warm-start R4 head-to-head (fgs_warmstart_r4_reset_prompt_20260924.md,
# Job 1): FGS (dagteam+backoff default) vs krylov_ilu_nfcache
# (persistent_plan), warm-started, wake-on, 4 revolutions (144 steps @ NT=36),
# Windows A (steps 1-36) and B (steps 109-144) harvested WITH transients.
#
# One exclusive zen3 node; arms run SEQUENTIALLY, one fresh Julia process per
# arm (constant hardware across the comparison; cold = zero-initial-guess in
# the same process, Ryan 2026-09-23). Champion placement (socket 0,
# interleave 0-3) at j=64, BLAS pinned to 1 for BOTH solver families.
#
# Required env on the submit line:
#   WSR4_PROJECT      campaign Julia project (Manifest dev-pointed at the
#                     campaign worktrees)
#   CAMPAIGN_PINS     pins TOML (annotated tags + SHAs, three repos)
#   WSR4_DATA_ROOT    consolidated data root for outputs (never the worktree)
#   FGS_TOL_ABS       explicit FGS stopping tolerance for THIS fixture
#   KNOBS_P/KNOBS_MAC/KNOBS_LEAF
#                     certified Krylov APPLY knobs, copied from the
#                     2026-09-22 R4 budget-500 tune row (the tune_phase2.csv
#                     lives in the thread-scaling run dirs, not the deployed
#                     tree, so the values are passed explicitly and recorded
#                     in the CSV)
# Optional:
#   ARMS              subset/order override (default: all seven)
#   RESUME_FROM_JOB_ID  reuse a previous run dir; ok arms skip
#
# Judge by outputs (STATUS_*/COMPLETED_* files), never sacct; logs are
# output-buffered (may freeze hours mid-compute — check CPU/outputs).
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${WSR4_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${WSR4_DATA_ROOT:?}" \
  "${FGS_TOL_ABS:?}" "${KNOBS_P:?}" "${KNOBS_MAC:?}" "${KNOBS_LEAF:?}"

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

run="$WSR4_DATA_ROOT/fgs-wsr4-${RESUME_FROM_JOB_ID:-$SLURM_JOB_ID}"
if [ -n "${RESUME_FROM_JOB_ID:-}" ] && [ -d "$run" ]; then
    rm -f "$run/COMPLETED"
else
    mkdir "$run"
fi
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$WSR4_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
numactl --hardware > "$run/numactl_hardware.txt" 2>&1 || true

# Champion placement (021 v23 champion; also used by the cold j64 rows)
ILV="--interleave=0-3 --cpunodebind=0-3"
numactl $ILV numactl --show > "$run/numactl_show_champion.txt" 2>&1 || true

THREADS=64
export JULIA_NUM_THREADS="$THREADS"
# BLAS pinned to 1 for BOTH arms (ILU is never BLAS-swept; stated in provenance)
export BENCH_BLAS_THREADS=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
export EXPECT_JULIA_THREADS="$THREADS" THREADING_MODE=multi
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3-blas1}"
export RUNG=R4 NT=36 N_STEPS=144 PHASE=phase3wsr4 SNAPSHOT_STRENGTHS=1
export FGS_PRECISION="${FGS_PRECISION:-f64}"   # f32full needs re-certification

( while true; do
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] heartbeat: fgs-wsr4 alive"
    sleep 300
  done ) &
HEARTBEAT_PID=$!
trap 'kill "$HEARTBEAT_PID" 2>/dev/null || true' EXIT

# Under the git-archive ("rsync") deployment mode there is no .git, so the
# CSV commit columns read "unknown" and provenance rests on the pins file +
# content manifests. Verify the deployed content actually matches the pinned
# manifests before burning node hours (same check the cold harness does).
# CONTENT_MANIFESTS = space-separated "<dir>:<manifest>" pairs.
if [ -n "${CONTENT_MANIFESTS:-}" ]; then
  for pair in $CONTENT_MANIFESTS; do
    dir="${pair%%:*}"; man="${pair#*:}"
    ( cd "$dir" && sha256sum --quiet -c "$man" ) || {
      echo "ERROR: deployed content mismatch in $dir vs $man" >&2; exit 1; }
    echo "content verified: $dir"
  done
fi

# RHPC resolves inputs/outputs via relative data/ paths — the campaign
# worktree must carry the standard `data` symlink to the shared data root
# (agent_policies/HPC.md), and the Das arc table must be reachable through it
[ -e data ] || { echo "ERROR: no data/ symlink in $(pwd) — deploy step missing" >&2; exit 1; }
[ -f data/p018_cs_l3p4_rs1_te_downwash_te.csv ] || {
  echo "ERROR: Das arc table not reachable via data/ — jobs die ~1 min in" >&2
  exit 1; }

ARMS="${ARMS:-fgs_cold fgs_prev fgs_proj1 fgs_proj2 ilu_nfcache_cold ilu_nfcache_prev ilu_nfcache_proj1}"
FAILED_COUNT=0
for arm in $ARMS; do
  if [ "$(cat "$run/STATUS_$arm" 2>/dev/null || true)" = ok ]; then
    echo "resume: $arm already ok — skipping"
    continue
  fi
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] === ARM=$arm ==="
  if env ARM="$arm" OUTDIR_OVERRIDE="$run" RUN_NAME="fgs_wsr4_R4_$arm" \
        numactl $ILV \
        julia --project="$WSR4_PROJECT" --startup-file=no -t "$THREADS" \
        benchmark/fgs_r4_warmstart_ab.jl > "$run/$arm.log" 2>&1 \
      && [ -f "$run/COMPLETED_$arm" ]; then
    echo "  $arm ok"
  else
    FAILED_COUNT=$((FAILED_COUNT + 1))
    echo "WARNING: $arm FAILED (continuing; see $run/$arm.log)"
  fi
done

printf 'completed failed_count=%s\n' "$FAILED_COUNT" > "$run/COMPLETED"
