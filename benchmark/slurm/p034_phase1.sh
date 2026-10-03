#!/usr/bin/env bash
#SBATCH --job-name=p034-ph1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# 034 Phase 1 consistency march: ONE (rung, arm) per job (Ryan 2026-10-02:
# R1+R2, all four arms, 3 cycles, SINGLE-thread mode, modest NON-exclusive
# allocations — Phase 1 judges accuracy/identity, not timing; exclusive-node
# timing belongs to Phase 2). Submit from the top level of the CAMPAIGN
# WORKTREE (tag campaign/p034-phase1-YYYYMMDD; data/ symlinked to the
# consolidated data root by prep_campaign_worktree.sh):
#
#   sbatch --job-name=p034-ph1-R1-fgs \
#     --export=ALL,P034_RUNG=1,ARM=fgs \
#     benchmark/slurm/p034_phase1.sh
#
# Walltime guidance (override --time at submit): R1 backslash/ilu/fgs ~4 h,
# R1 gmres ~12 h; R2 backslash/ilu/fgs ~12 h, R2 gmres ~48 h (local -t 1
# sizing: R2 ~27 s/step x 495 steps; gmres ~8x the ilu per-step cost).
# Partition/QOS are chosen at submit time from a fresh slurm-availability
# probe — non-exclusive, so any zen3/zen4 CPU partition with a normal QOS
# works; timings are context only, HARDWARE_TAG still records the node class.
set -euo pipefail

module load julia

# ---- NFS precompile-lock guard (pattern from p2_unsteady.sh, 2026-08-25;
# see the orc precompile-race postmortem) --------------------------------
export JULIA_PKG_PRECOMPILE_AUTO=0
if ! command -v flock >/dev/null 2>&1; then
  echo "[$(date -u +%FT%TZ)] WARNING: flock(1) not found; SKIPPING the shared" \
       "precompile. This job compiles in memory (slower, but safe)."
else
  ( if flock -E 99 -w 3600 9; then
      echo "[$(date -u +%FT%TZ)] precompile lock acquired"
      julia --project=. --startup-file=no -e 'using Pkg; Pkg.precompile()' \
        || echo "[$(date -u +%FT%TZ)] WARNING: Pkg.precompile() failed; this" \
                "job will compile in memory instead"
    else
      rc=$?
      if [ "$rc" = 99 ]; then
        echo "[$(date -u +%FT%TZ)] WARNING: precompile lock timed out after" \
             "1 h; this job will compile in memory instead"
      else
        echo "[$(date -u +%FT%TZ)] WARNING: flock failed (exit $rc); this job" \
             "will compile in memory instead"
      fi
    fi
  ) 9>"$HOME/.julia/flowpanel-034-precompile.lock"
fi
echo "[$(date -u +%FT%TZ)] precompile stage done"
# ------------------------------------------------------------------------

ARM=${ARM:?set ARM (backslash|krylov_gmres|krylov_ilu_nfcache|fgs)}
P034_RUNG=${P034_RUNG:?set P034_RUNG (1-4)}

export THREADING_MODE=single
export EXPECT_JULIA_THREADS=1
export P034_RUNG
export P034_ARMS="$ARM"
export P034_NCYCLES=${P034_NCYCLES:-3.0}
export P034_OUTDIR=${P034_OUTDIR:-data/p034_phase1}
export HARDWARE_TAG="${HARDWARE_TAG:-orc-${SLURM_JOB_PARTITION:-?}-nonexclusive}"

echo "p034 Phase 1 consistency"
echo "  repo:    $(pwd)"
echo "  rung:    R$P034_RUNG   arm: $ARM   n_cycles: $P034_NCYCLES"
echo "  node:    ${SLURMD_NODENAME:-?} partition: ${SLURM_JOB_PARTITION:-?}"
echo "  outdir:  $P034_OUTDIR (judge from the CSVs there, never this log)"

julia --project=. --startup-file=no --compiled-modules=existing -t 1 \
    benchmark/p034_phase1_consistency.jl
