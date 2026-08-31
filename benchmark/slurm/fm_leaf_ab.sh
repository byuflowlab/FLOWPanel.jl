#!/usr/bin/env bash
#SBATCH --job-name=fm-leaf-ab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=m12
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --mem=500G
#SBATCH --requeue
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# 021 Phase 1 R4-R7 tuning pipeline for one rung (RUNG from --export):
# stage 1 tune_fmm (dense-free) -> stage 2 FGS tau-ladder -> stage 3
# preconditioner sweep ladder.
# Submit from the top level of the FLOWPanel.jl checkout, e.g.:
#   sbatch --job-name=p1-tune-R4 --mem=48G --export=ALL,RUNG=R4 benchmark/slurm/p1_tune.sh
# Override --time per rung at sbatch time (R6/R7 need more). Do NOT
# override --mem/--constraint/--exclusive: they are the pinned hardware.
#
# STAGES (2026-08-22): colon- or comma-separated subset of 1:2:3, default all
# three. The rung cost grows ~N^1.6 and R7 extrapolates to ~114 h against a
# 3-day partition ceiling, so R7 must be run one stage per job, chained
# (stage 2 consumes stage 1's tune.csv, stage 3 consumes stage 2's
# fgstune_selected.csv):
#   sbatch --job-name=p1-tune-R7-s1 --time=3-00:00:00 \
#     --export=ALL,RUNG=R7,STAGES=1 benchmark/slurm/p1_tune.sh
#   sbatch --job-name=p1-tune-R7-s2 --time=3-00:00:00 \
#     --dependency=afterany:<s1_id> \
#     --export=ALL,RUNG=R7,STAGES=2 benchmark/slurm/p1_tune.sh
#
# CHAIN WITH afterany, NOT afterok (Ryan ruling 2026-08-24). An afterok cascade
# turns one upstream failure into N silent CANCELLED-with-Reason=Dependency
# jobs, which has already been misdiagnosed twice on this campaign (08-22 and
# again on the 13306603 -> 13306604/13306605 chain). Under afterany the
# downstream job RUNS and hits the stage input assert below, which names the
# missing file and exits 1 — a diagnosable failure instead of a vanished job.
# The stage scripts APPEND rather than checkpoint-skip, so never re-run a
# stage that already wrote its CSVs into results/phase1/multi/<RUNG>/ — it
# duplicates rows. This supersedes p1_tune_s3.sh (STAGES=3).
set -euo pipefail

module load julia

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

THREADS=64
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
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-gaussian}"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

export EXPECT_JULIA_THREADS="$THREADS"
export THREADING_MODE=multi
export CACHE_B=1
export PER_RUNG_DIR=1
: "${RUNG:?set RUNG via --export=ALL,RUNG=Rn}"

# Diagnostic: sweep leaf_size at fixed p/mac and report solve cost, to test
# whether the PUBLISHED (pre-2026-08-24) ladder carries the same ~2x leaf error
# that R1 did. Writes nothing to the campaign results tree.
# NB benchmark/fm_leaf_ab.jl's bc_rel_l2 column is INVALID (wrong args to
# bc_error!); only its timings are sound.
#   sbatch --job-name=fm-leaf-ab-R2 --export=ALL,RUNG=R2 benchmark/slurm/fm_leaf_ab.sh
echo "leaf A/B sweep for $RUNG"
julia --project=. -t "$THREADS" benchmark/fm_leaf_ab.jl
echo "leaf A/B sweep complete for $RUNG"
