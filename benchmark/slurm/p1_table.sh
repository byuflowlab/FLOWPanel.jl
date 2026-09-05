#!/usr/bin/env bash
#SBATCH --job-name=p1-table
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=m12
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --mem=500G
#SBATCH --requeue
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# 021 Phase 1 solver x rung table for one rung + threading mode.
# Submit from the top level of the FLOWPanel.jl checkout, e.g.:
#   sbatch --job-name=p1-table-R4-multi \
#     --export=ALL,RUNG=R4,MODE=multi,CONFIGS=backslash_ldiv,krylov_gmres,krylov_jacobi,krylov_ilu,fgs,fgmres_fgs \
#     benchmark/slurm/p1_table.sh
# MODE=single runs the 1-thread/BLAS-1 mode (ruling 6). Dense configs
# (backslash_ldiv) only fit through R5 (ledger cost-ceiling rule).
set -euo pipefail

module load julia

# ---- NFS precompile-lock guard (added 2026-08-25) --------------------------
# See p2_tune.sh for the full diagnosis. Short version: the 2026-08-24 Phase 1+2
# launch fanned 13 jobs onto one NFS depot with auto-precompilation live; one
# node won the per-package pidfile lock and the other seven blocked ~12 h inside
# package loading, before any driver code ran, burning four jobs' full walltime.
# This flag only suppresses the implicit per-job attempt; the serialized
# precompile block immediately below builds the cache, so no manual pre-step
# is required before a fan-out.
export JULIA_PKG_PRECOMPILE_AUTO=0

# ---- Serialized precompile (added 2026-08-25) ------------------------------
# One writer at a time across the whole cluster. /home is NFSv4.1 with
# local_lock=none (verified 2026-08-25), so locks are delegated to the VAST
# server and flock(2) serializes ACROSS NODES; on an NFSv3 mount it would be
# node-local and this would silently allow one builder per node.
#
# No "is it warm?" probe: Pkg.precompile() is already a fast no-op on a warm
# depot, and probing one package (`using FLOWPanel`) would miss driver-only
# deps such as CSV/DataFrames. Warm => every job passes through in seconds.
# Cold => the first job builds and the rest find it warm.
#
# BOUNDED AND NON-FATAL BY DESIGN. The 2026-08-24 failure was an UNBOUNDED
# wait on a wedged lock holder. Here a timeout or a failed build only logs and
# falls through, because the run below uses --compiled-modules=existing and so
# can always compile in memory instead. A slow job beats a deadlocked one.
if ! command -v flock >/dev/null 2>&1; then
  # Degrade SAFELY and LOUDLY. With no lock available, running Pkg.precompile()
  # from every job at once would recreate the exact 2026-08-24 thundering herd,
  # so skip the shared build entirely and let --compiled-modules=existing
  # compile in memory: slower per job, but structurally unable to deadlock.
  echo "[$(date -u +%FT%TZ)] WARNING: flock(1) not found; SKIPPING the shared" \
       "precompile. This job compiles in memory (slower, but safe)."
else
  ( if flock -E 99 -w 3600 9; then
      echo "[$(date -u +%FT%TZ)] precompile lock acquired"
      julia --project=. --startup-file=no -e 'using Pkg; Pkg.precompile()' \
        || echo "[$(date -u +%FT%TZ)] WARNING: Pkg.precompile() failed; this" \
                "job will compile in memory instead"
    else
      rc=$?   # must be the FIRST command here: this is flock's exit status
      if [ "$rc" = 99 ]; then
        echo "[$(date -u +%FT%TZ)] WARNING: precompile lock timed out after" \
             "1 h; this job will compile in memory instead"
      else
        echo "[$(date -u +%FT%TZ)] WARNING: flock failed (exit $rc); this job" \
             "will compile in memory instead"
      fi
    fi
  ) 9>"$HOME/.julia/flowpanel-021-precompile.lock"
fi
echo "[$(date -u +%FT%TZ)] precompile stage done"
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

# ---- Pinned benchmark hardware (021 ruling 2026-08-24) --------------------
# Every 021 cost number must be comparable across rungs, phases and months, so
# the hardware is pinned in the SBATCH header above, NOT left to the submit
# line: m12 / zen3 / --exclusive / --mem=0. Measured node spec (scontrol,
# 2026-08-24): Sockets=2, CoresPerSocket=64, ThreadsPerCore=1 => CPUTot=128,
# RealMemory=524288 MB, MaxMemPerNode=524288.
#
# Why each pin:
#   --constraint=zen3  the ten 2026-08-22 jobs all ran zen3 only because the
#                      submitter typed it by hand; p1_table.sh's documented
#                      submit line omitted it entirely. One forgotten flag
#                      silently makes a rung incomparable.
#   --exclusive        WHOLE node allocated (128 CPUs) while Julia uses only
#                      THREADS=64. Ryan ruling 2026-08-24, after briefly trying
#                      a 64-CPU half-node allocation. Two reasons, both about
#                      measurement validity rather than speed:
#                        1. These partitions are OverSubscribe=YES:4. On a
#                           shared node a stranger's job contends for memory
#                           bandwidth, and this FMM/ILU workload is
#                           bandwidth-bound. That noise is NOT reproducible —
#                           it depends on who else happens to land there.
#                        2. NUMA placement. These are 2 sockets x 64 cores. A
#                           64-CPU cgroup can bind all 64 threads onto ONE
#                           socket, halving available memory bandwidth; on an
#                           exclusive node the 64 threads spread across both
#                           sockets. So the same thread count can time
#                           differently under the two allocations.
#                      It also matches the pre-2026-08-24 ladder, so
#                      cross-generation comparison isolates the FMM change
#                      instead of confounding it with allocation.
#                      Cost: ~2x fewer concurrent jobs. Accepted.
#   --mem=500G         effectively the whole 512 GiB node. Explicit rather than
#                      --mem=0 so the requirement is documented; covers R5's
#                      ~94 GiB dense G and the 32 GiB near-field cache cap.
#                      Replaces the 08-22 generation's inconsistent
#                      64G/128G/192G/256G spread.
#   --requeue          jobs run on the preemptible `standby` QOS (see below), so
#                      a preempted job must come back rather than vanish. Safe
#                      ONLY because the drivers now resume: each row is flushed
#                      as it completes and a restart skips every row already in
#                      the CSV, so a requeue never duplicates rows (2026-08-24).
#
# PARTITION/QOS ARE THE ONE PERMITTED OVERRIDE. These three partitions are
# spec-IDENTICAL (2 sockets x 64 zen3 cores, ThreadsPerCore=1, CPUTot=128,
# RealMemory=524288 MB, ActiveFeatures include m12,amd,zen3), so a job is
# equally comparable on any of them:
#     m12       public,  MaxTime  3 d, 136 nodes, --qos=normal  (contended)
#     physics2  private, MaxTime  7 d,  10 nodes, --qos=standby (preemptible)
#     m12pws    private, MaxTime 28 d,   1 node,  --qos=standby (preemptible)
# e.g.  sbatch --partition=physics2 --qos=standby --time=7-00:00:00 ...
# Everything else (--constraint, --exclusive, --mem) stays as pinned, and
# the CPU assert below is the backstop. Override --time and --job-name freely.
export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-2x64c-512g}"
ALLOC_CPUS=${SLURM_CPUS_ON_NODE:-0}       # CPUs allocated to this job on the node
NODE_CPUS=$(getconf _NPROCESSORS_ONLN)    # physical CPUs on the node
AFFINITY=$(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')
echo "hardware pin: partition=${SLURM_JOB_PARTITION:-?} node=${SLURMD_NODENAME:-?}" \
     "alloc_cpus=$ALLOC_CPUS node_cpus=$NODE_CPUS mem=${SLURM_MEM_PER_NODE:-?}MB" \
     "affinity=$AFFINITY tag=$HARDWARE_TAG"
# The node must be the pinned 128-core zen3 part AND this job must hold all of
# it, so the 64 Julia threads get uncontended memory bandwidth across both
# sockets. `affinity` is echoed (not asserted) so thread placement is on the
# record for every row rather than assumed.
# NB: do NOT use `nproc` here — GNU nproc honours OMP_NUM_THREADS, which Slurm
# presets to 1 under --cpus-per-task, so it under-reports (cost a full fan-out
# on 2026-08-24). SLURM_CPUS_ON_NODE is the reliable source.
if [ "$NODE_CPUS" != "128" ]; then
  echo "ERROR: expected a 128-core zen3 node, got $NODE_CPUS. This job is NOT on" >&2
  echo "       the campaign's pinned hardware; its timings would not be" >&2
  echo "       comparable. Refusing to run." >&2
  exit 1
fi
if [ "$ALLOC_CPUS" != "128" ]; then
  echo "ERROR: expected an exclusive 128-CPU allocation, got $ALLOC_CPUS. A shared" >&2
  echo "       node lets a co-tenant contend for memory bandwidth and can bind all" >&2
  echo "       64 threads to one socket. Timings would not be comparable." >&2
  exit 1
fi
# --------------------------------------------------------------------------

: "${RUNG:?set RUNG via --export}"
: "${MODE:?set MODE=multi|single via --export}"
: "${CONFIGS:?set CONFIGS via --export (colon-separated, e.g. fgs:fgmres_fgs)}"
CONFIGS="${CONFIGS//:/,}"   # colon-separated at sbatch --export (commas collide)
K_REPS="${K_REPS:-5}"

if [ "$MODE" = "multi" ]; then THREADS=64; else THREADS=1; fi
export JULIA_NUM_THREADS="$THREADS"

# BRAINSTORM/025: the codebase default filament regularization changed to
# Gaussian on 2026-08-20. This campaign was pinned to the legacy Vatistas
# family until 2026-08-22, when Ryan moved it to Gaussian; the pre-08-22 rungs
# are to be RERUN under Gaussian later (see log.md). Override at submission
# via --export if intended.
# NOTE inert for the Phase-1/2b drivers: phase1_case.jl is a frozen
# single-step solve that never sheds a wake, so no filament kernel is
# reached (measured 2026-08-22 — R1 solution bit-identical under both).
# It DOES bite the wake-carrying Phase-2 unsteady arms.
# LineGauss is the campaign default from 2026-09-05 (Ryan), following the
# codebase default flip of 2026-08-29 (task 052d): it is the exact
# segment-kernel/Gaussian-blob convolution, of which Gaussian is the
# infinite-line limit, so radius_inflation genuinely bounds the direct/
# expansion mismatch instead of leaving Gaussian's open along-line error
# channel (which cost 9e-4 relU on far-field just-shed wake particles).
# Rows dated before 2026-09-05 in this campaign ran GAUSSIAN; set
# FLOWPANEL_FILAMENT_REG=gaussian to reproduce them.
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

export EXPECT_JULIA_THREADS="$THREADS"
export THREADING_MODE="$MODE"
export CACHE_B=1
export PER_RUNG_DIR=1
export KNOBS_MODE=multi
export CONFIGS K_REPS

echo "021 Phase 1 table"
echo "  repo:    $(pwd)"
echo "  rung:    $RUNG  mode: $MODE  threads: $THREADS  k_reps: $K_REPS"
echo "  configs: $CONFIGS"

julia --project=. --startup-file=no --compiled-modules=existing -t "$THREADS" benchmark/rotor_hover_solver_phase1_table.jl

echo "table complete for $RUNG/$MODE"
