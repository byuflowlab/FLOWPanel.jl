#!/usr/bin/env bash
# ORC jobs must follow BYU_ORC_AGENTS.md. Invoke from the campaign worktree.
set -euo pipefail
jthreads=${1:?Julia threads}; bthreads=${2:?BLAS threads}; driver=${3:?driver}
shift 3
[[ "$jthreads" =~ ^[1-9][0-9]*$ && "$bthreads" =~ ^[1-9][0-9]*$ ]] || exit 2
export JULIA_NUM_THREADS="$jthreads" EXPECT_JULIA_THREADS="$jthreads"
export THREADING_MODE=multi BENCH_BLAS_THREADS="$bthreads"
export OMP_NUM_THREADS="$bthreads" OPENBLAS_NUM_THREADS="$bthreads"
export BLAS_NUM_THREADS="$bthreads" MKL_NUM_THREADS="$bthreads"
export VECLIB_MAXIMUM_THREADS="$bthreads" BLIS_NUM_THREADS="$bthreads"
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
export FLOWPANEL_FILAMENT_REG=linegauss CACHE_B=0 SKIP_B=0 MEMORY_GIB=500
: "${COLD_PROJECT:?dedicated environment required}"
# The batch allocation owns the affinity mask; this launcher does not spawn srun.
exec julia --startup-file=no --project="$COLD_PROJECT" -t "$jthreads" "$driver" "$@"
