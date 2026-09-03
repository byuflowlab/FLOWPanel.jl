#!/usr/bin/env bash
#SBATCH --job-name=fp-051-price
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=00:45:00
# qos=test dropped after 13391706: test-QoS nodes gave ~24 GB/s effective
# gemv bandwidth and noisy timings; default QoS for representative hardware.
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# Task 051 Stage 3 (Job C): price the production p018 CPU solve.
#
# CPU-only (no GPU request -- the H200 queue wait does not apply): warm-starts
# the p018_cs_f1_l3p4 restart set exactly like fm051_parity.sh full mode, then
# times solve_formulation! cold + N warm repeats (benchmark/fm051_solve_pricing.jl)
# against the ~0.78 s U-only headroom (3.3 - [0.373 + 0.124 + 2.02]).
# Env block below is the same verbatim p018 recipe as fm051_parity.sh -- if the
# launcher recipe changes, re-copy BOTH scripts.
#
# Submit from the cluster checkout top level:
#   sbatch benchmark/slurm/fm051_solve_pricing.sh

# -u is deferred until after /etc/profile + module load: the cluster's
# /etc/profile.d/debuginfod.sh reads an unbound DEBUGINFOD_URLS and would
# kill the job under set -u (failed job 13306465; fm051_run.sh pattern)
set -eo pipefail

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
# ~/FLOWPanel-046 is the SYNCED tree (benchmark/slurm/fm051_parity_submit.sh);
# the standing ~/projects/FLOWPanel.jl checkout only supplies read-only data
# via symlinks (data/p018_cs_f1_l3p4, examples/data)
EXPECTED_REPO="$HOME/FLOWPanel-046"
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }
[[ -e data/p018_cs_f1_l3p4 && -e examples/data ]] || { echo "ERROR: data symlinks missing; run benchmark/slurm/fm051_parity_submit.sh locally first" >&2; exit 2; }

source /etc/profile
module load julia/1.11.7-6bmogfl
set -u

echo "=== node: $(hostname)"
echo "(CPU-only job: no GPU requested)"

# Julia environment.  FM051_ENV must resolve FLOWPanel + FLOWVPM +
# FastMultipole (devved) AND CUDA.jl; $HOME/fm051env (built by
# fm051_parity_submit.sh: fm048env recipe + FLOWPanel + VSPGeom; fm048env
# itself does NOT dev FLOWPanel).  Set FM051_ENV=. to run against the repo's
# own Manifest (works for FM051_SEAM=host; the :cuda arm needs CUDA.jl present,
# and the seam then reports "CUDA not functional" and degrades to host).
ENVDIR="${FM051_ENV:-$HOME/fm051env}"
export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" \
       OPENBLAS_NUM_THREADS="$THREADS" BLAS_NUM_THREADS="$THREADS" \
       MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg

# The harness arms/disarms the seam itself and refuses to start if this is set.
unset FLOWPANEL_GPU_INFLUENCE || true

# ---- production config: verbatim p018_cs_f1_l3p4 (launcher fixed block) -----
export RHPC_MESH=45_185_ct4
export RPM=5400
export BERNOULLI_ONLY=true
export RUN_KJ=false
export SAVE_VTK=true
export SPINUP_REVS=1.5
export SPINUP_START_FRACTION=0.4
export MAGVINF_START=0.0
export MAGVINF_PEAK=5.0
export MAGVINF_END=0.0
export FREESTREAM_RAMP_REVS=1.0
export FREESTREAM_HOLD_REVS=1.5
export FREESTREAM_WITHDRAW_REVS=4
export NREVS=10
export CONVERGENCE_MEAN_TOL=0.005
export CONVERGENCE_PTP_TOL=0.02
export DAS_KINEMATIC_ARC=false
# ---- p018_* group defaults ----
unset RELAX_FILTER_DOWNSTREAM_R || true
export RHPC_FORMULATION=velocity
export NT=36
export TRUNCATION_DEPTH_R=4
export RELAX_RLXF=0.3
export SETTLE_REVS=22
export CONVERGENCE_REVS=10
export PARTICLE_SHEDDING=sigma_overlap
export DAS_ETA_KINEMATIC=1.0
# ---- p018_cs_f1_l3p4 case arm ----
export OVERLAP=2.75
export P_PER_STEP=12
export MERGE_R_FACTOR=0.0055
export NWAKEROWS=1
export SIGMA_CHORD_FRACTION=0.313
export SIGMA_FLOOR_R=0
export DAS_SIGMA_LAMBDA=3.4
export DAS_CURVATURE_BETA=0.6
export CORE_SPREADING_ACTIVE=true
export WAKE_CORE_BETA=1e9

# ---- warm-start selection (p023 names; consumed by p018_mature_wake_timing.jl,
#      which fm051_pass_parity.jl loads into a private module) ----------------
export BENCH_RESTART_RUN="${BENCH_RESTART_RUN:-p018_cs_f1_l3p4}"
export RESTART_STEP="${RESTART_STEP:--1}"        # -1 = highest step in the body PVD
export BENCH_NSTEPS="${BENCH_NSTEPS:-1}"         # the parity step only
export BENCH_WRITE_OUTPUT=false                  # the parity harness writes no VTK

# ---- parity harness knobs --------------------------------------------------
export FM051_MODE="${FM051_MODE:-full}"
export FM051_SEAM="${FM051_SEAM:-host}"   # include-time requirement; the pricing script never arms the seam
export FM051_GATE="${FM051_GATE:-1e-3}"
export FM051_GATE_P1="${FM051_GATE_P1:-1e-11}"
export FM051_GATE_SFS="${FM051_GATE_SFS:-1e-10}"

mkdir -p logs/slurm benchmark/results

echo "FM051 Stage 3 solve pricing -- CPU"
echo "  repo:$PWD  env:$ENVDIR  threads:$THREADS  host:$(hostname)  job:${SLURM_JOB_ID:-none}"
echo "  restart run:$BENCH_RESTART_RUN  restart step:$RESTART_STEP"
echo "  started $(date '+%F %T')"

status=0
julia --project="$ENVDIR" -t "$THREADS" benchmark/fm051_solve_pricing.jl || status=$?

echo "Finished ($(date '+%F %T')) with status $status."
echo "Results: benchmark/results/fm051_solve_pricing.csv"
exit $status
