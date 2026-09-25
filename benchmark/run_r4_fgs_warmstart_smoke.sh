#!/usr/bin/env bash
# 021 warm-start campaign — LOCAL smoke (fgs_warmstart_r4_reset_prompt_20260924
# Job 1): every arm on R1, few steps, <=4 threads, functional (uncertified)
# knobs. Verifies arm plumbing (config construction, warm-start seeding,
# per-step CSV incl. t_project + setup split, snapshots, STATUS/COMPLETED
# sentinels) before anything is staged for HPC. NOT a performance or accuracy
# measurement.
set -euo pipefail
cd "$(dirname "$0")/.."

SMOKE_ROOT="${SMOKE_ROOT:-benchmark/results/wsr4_smoke_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$SMOKE_ROOT"
SMOKE_ABS="$(cd "$SMOKE_ROOT" && pwd)"
echo "smoke output root: $SMOKE_ABS"

# Leg-structured stage list mirroring the campaign launcher (checkpoint +
# restart layout, Ryan 2026-09-24). NT=4 makes a smoke "revolution" 4 steps:
# ckpt = 12 steps w/ VTK, winA = 4 steps, winB = restart at step 12 + 4 steps.
# Subset of warm arms — every leg type runs in both solver families.
STAGES="${STAGES:-fgs_cold:ckpt ilu_nfcache_cold:ckpt \
fgs_prev:winA ilu_nfcache_prev:winA \
fgs_cold:winB fgs_prev:winB ilu_nfcache_prev:winB}"

fail=0
for stage in $STAGES; do
  arm="${stage%%:*}"; lg="${stage#*:}"
  armleg="${arm}_${lg}"
  echo "=== ARM=$arm LEG=$lg ==="
  if env ARM="$arm" WSR4_LEG="$lg" RUNG=R1 NT=4 \
      EXPECT_JULIA_THREADS=4 THREADING_MODE=multi BENCH_BLAS_THREADS=1 \
      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLAS_NUM_THREADS=1 \
      MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 \
      FLOWPANEL_FILAMENT_REG=linegauss \
      KNOBS_P=17 KNOBS_MAC=0.5 KNOBS_LEAF=21 \
      FGS_P=6 FGS_MAC=0.3 FGS_LEAF=50 FGS_INNER=3 FGS_TOL_ABS=1e-6 \
      NFCACHE_MAX_GIB=4 \
      OUTDIR_OVERRIDE="$SMOKE_ABS" \
      julia --project=. -t 4 benchmark/fgs_r4_warmstart_ab.jl \
      > "$SMOKE_ABS/$armleg.log" 2>&1; then
    echo "  ok ($(grep -c "^R1" "$SMOKE_ABS/unsteady.csv" 2>/dev/null || echo '?') CSV rows so far)"
  else
    echo "  FAILED — tail of $SMOKE_ABS/$armleg.log:"
    tail -15 "$SMOKE_ABS/$armleg.log"
    fail=1
  fi
done

echo "=== sentinels ==="
ls -1 "$SMOKE_ABS" | grep -E "STATUS_|COMPLETED_" || true
[ "$fail" = 0 ] && echo "SMOKE PASS" || { echo "SMOKE FAIL"; exit 1; }
