#!/usr/bin/env bash
#SBATCH --job-name=p033-ts-cold-r4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=12:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/slurm/p033-ts-cold-r4-%A_%a.out
#SBATCH --error=logs/slurm/p033-ts-cold-r4-%A_%a.err
# p033 R4 thread-scaling, Job A (cold): latest FGS (threaded setup, 033 B-I2)
# vs krylov_ilu_nfcache at FIXED champion knobs across j = 1/8/16/32/64 — NO
# per-j re-tuning (Ryan 2026-10-02 ruling; supersedes the re-tune layout of
# run_r4_thread_scaling.slurm.sh for THIS campaign). Knobs are the certified
# R4 j64 set and are recorded in every CSV row; the harvest labels them as
# fixed-across-j transplants.
#
#   FGS arm — benchmark/fgs_setup_ab.jl ARMS=old,new CERT_SOLVE=1: one-time
#             setup (serial vs threaded, full_ctor rows), bitwise certs, and
#             AB_SOLVE_K cold solves per arm (cold_solve rows; champion TOML
#             knobs + tolerance, dagteam f32full as retained).
#   ILU arm — benchmark/rotor_hover_solver_phase2.jl at pre-seeded knobs:
#             this launcher WRITES the tune_phase2.csv knob rows (budget 0 =
#             P15/MAC0.55/leaf21, budget 500 = P12/MAC0.55/leaf48 — the
#             2026-09-22 certified rows) so phase2_knobs() resolves without a
#             tuner pass. Componentized setup + t_solve_min land in phase2.csv.
#
# Placement, BLAS=1, per-task run dirs, STATUS_*/COMPLETED sentinels, and the
# judge-by-outputs rule all follow run_r4_thread_scaling.slurm.sh verbatim.
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"

LADDER=(1 8 16 32 64)
J="${LADDER[$SLURM_ARRAY_TASK_ID]}"

# ---- pinned-hardware assertions (021 ruling 2026-08-24) ----------------------
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

run="$COLD_DATA_ROOT/p033ts-cold-j$J-${RESUME_FROM_JOB_ID:-$SLURM_ARRAY_JOB_ID}"
if [ -n "${RESUME_FROM_JOB_ID:-}" ] && [ -d "$run" ]; then
    rm -f "$run/COMPLETED"
else
    mkdir "$run"
fi
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
numactl --hardware > "$run/numactl_hardware.txt" 2>&1 || true
echo "julia_threads=$J knobs=fixed-transplant-j64" > "$run/point.txt"

ILV="--interleave=0-3 --cpunodebind=0-3"
numactl $ILV numactl --show > "$run/numactl_show.txt" 2>&1 || true

export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3-blas1}"
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export JULIA_DEBUG="${JULIA_DEBUG:-loading}"

( while true; do
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] heartbeat: p033ts-cold j$J alive"
    sleep 300
  done ) &
HEARTBEAT_PID=$!
trap 'kill "$HEARTBEAT_PID" 2>/dev/null || true' EXIT

status() { printf '%s\n' "$2" > "$run/STATUS_$1"; }
stage_ok() { [ "$(cat "$run/STATUS_$1" 2>/dev/null || true)" = ok ]; }

# ---- FGS arm: setup A/B + certified cold solves at this j --------------------
fgs_ok=1
if stage_ok fgs_setup_ab; then
    echo "resume: fgs_setup_ab already ok — skipping"
else
    mkdir -p "$run/fgs-setup"
    if env RUNG=R4 ARMS=old,new AB_K=2 CERT_SOLVE=1 AB_SOLVE_K=3 \
            SKIP_B=0 CACHE_B=1 \
            EXPECT_JULIA_THREADS="$J" THREADING_MODE=multi \
            BENCH_BLAS_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
            BENCH_CASE_ROOT="$run/fgs-case" numactl $ILV \
            julia --project="$COLD_PROJECT" --startup-file=no \
            --compiled-modules=existing -t "$J" \
            benchmark/fgs_setup_ab.jl "$run/fgs-setup" \
            > "$run/fgs-setup.log" 2>&1; then
        status fgs_setup_ab ok
    else
        status fgs_setup_ab FAILED
        fgs_ok=0
    fi
fi

# ---- ILU arm: pre-seed the fixed knob rows, then measure ---------------------
# Row identity in tune_phase2.csv is (rung, mem_budget_gib, julia_threads);
# phase2_knobs() matches on (rung, budget) with latest-row-wins, so per-task
# PHASE2_OUTDIR isolation (one task = one j = one dir) keeps this unambiguous.
ilu_ok=1
mkdir -p "$run/ilu"
knobs_csv="$run/ilu/tune_phase2.csv"
if [ ! -s "$knobs_csv" ]; then
  cat > "$knobs_csv" <<EOF
rung,mesh_file,n_panels,mem_budget_gib,cached,expansion_order,multipole_acceptance,leaf_size,t_solve_warm,niter,bc_rel_l2_certified,bc_certified,t_cache_build,cache_bytes,mem_total_predicted,mem_total_measured,dense_bytes,cache_capped,tune_reps,tune_abandon_factor,tune_timed_out,t_tune,n_candidates,n_abandoned,cache_tree,threading_mode,julia_threads,blas_threads,commit,fm_commit,date,hardware_tag,notes
R4,dji9443_20260813_65_209_capped_captess4.msh,58192,0,false,15,0.55,21,0,0,0,false,0,0,0,0,0,false,0,0,false,0,0,0,true,multi,$J,1,preseeded,preseeded,$(date -u +%Y-%m-%d),$HARDWARE_TAG,pre-seeded fixed knobs (2026-09-22 certified budget-0 row; no tune) p033 threadscale
R4,dji9443_20260813_65_209_capped_captess4.msh,58192,500,true,12,0.55,48,0,0,0,false,0,0,0,0,0,false,0,0,false,0,0,0,true,multi,$J,1,preseeded,preseeded,$(date -u +%Y-%m-%d),$HARDWARE_TAG,pre-seeded fixed knobs (2026-09-22 certified budget-500 row; no tune) p033 threadscale
EOF
fi

if stage_ok ilu_measure; then
    echo "resume: ilu_measure already ok — skipping"
else
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

printf 'completed j=%s fgs_ok=%s ilu_ok=%s\n' "$J" "$fgs_ok" "$ilu_ok" \
    > "$run/COMPLETED"
