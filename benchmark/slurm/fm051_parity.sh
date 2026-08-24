#!/usr/bin/env bash
#SBATCH --job-name=fp-051-parity
#SBATCH --nodes=1
#SBATCH --gpus=h200:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# Task 051 stage 2/3 GATE: pass-by-pass CPU-vs-GPU parity of the FLOWPanel
# rectangular-influence seam, on the MATURE p018 wake.
#
# Runs benchmark/fm051_pass_parity.jl in FM051_MODE=full: warm-starts the
# retained p018_cs_f1_l3p4 restart set (mirror of simulate_warmstart! sections
# 2.5-5, kinematic replay included) and then replays each cross-influence pass
# of the next timestep twice -- once with the seam :off (FastMultipole.fmm!)
# and once with the seam armed (FASTMULTIPOLE direct_rectangular! on the H200)
# -- comparing per-target outputs.  Results land in benchmark/results/;
# data/ is only ever read.
#
# Job pattern: FLOWVPM scripts/fm051_run.sh (julia PINNED to 1.11.7 --
# 1.12.6 segfaults JIT-compiling the device step, job 13058191) crossed with
# benchmark/slurm/p023_mature_timing.sh (production p018 env block).
#
# Submit from the cluster checkout top level:
#   sbatch benchmark/slurm/fm051_parity.sh
# Host-kernel (no GPU) control run -- same pair math, no device:
#   sbatch --gpus=0 --export=ALL,FM051_SEAM=host benchmark/slurm/fm051_parity.sh
# Mini smoke (seconds, no restart data needed):
#   sbatch --gpus=0 --export=ALL,FM051_MODE=mini benchmark/slurm/fm051_parity.sh
#
# The env block below is a verbatim copy of the p018_cs_f1_l3p4 recipe from
# benchmark/slurm/p023_mature_timing.sh (itself copied from
# examples/run_dji9443_hover_ct_hpc.slurm.sh).  If the launcher recipe
# changes, re-copy -- drift here means the restored state is not the run that
# produced the restart files.

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
module load cuda julia/1.11.7-6bmogfl
set -u

echo "=== node: $(hostname)"
nvidia-smi -L || echo "(no GPU allocated -- host-kernel run)"
echo "CUDA_HOME=${CUDA_HOME:-unset}"

# Julia environment.  FM051_ENV must resolve FLOWPanel + FLOWVPM +
# FastMultipole (devved) AND CUDA.jl; $HOME/fm051env (built by
# fm051_parity_submit.sh: fm048env recipe + FLOWPanel + VSPGeom; fm048env
# itself does NOT dev FLOWPanel).  Set FM051_ENV=. to run against the repo's
# own Manifest (works for FM051_SEAM=host; the :cuda arm needs CUDA.jl present,
# and the seam then reports "CUDA not functional" and degrades to host).
ENVDIR="${FM051_ENV:-$HOME/fm051env}"
export FASTMULTIPOLE_FORCE_CUDA_LOAD=1
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
export FM051_SEAM="${FM051_SEAM:-cuda}"
export FM051_GATE="${FM051_GATE:-1e-3}"
export FM051_GATE_P1="${FM051_GATE_P1:-1e-11}"
export FM051_GATE_SFS="${FM051_GATE_SFS:-1e-10}"

mkdir -p logs/slurm benchmark/results

# ---- stage 0 (combined-job rule): rect-kernel profile ------------------------
# FM051_PROFILE=1 (default in full mode) reruns the FLOWVPM rect bench with
# its FM051_BENCH_PROFILE decomposition (pass-2 tag split, pass-1 occupancy
# probe) on the same H200 allocation before the parity harness, so the
# optimization data rides this job instead of costing a second queue wait.
FM051_PROFILE="${FM051_PROFILE:-$([ "$FM051_MODE" = full ] && echo 1 || echo 0)}"
if [[ "$FM051_PROFILE" == "1" && -d "$HOME/FLOWVPM-046" ]]; then
    echo "=== stage 0: rect-kernel profile (FLOWVPM-046 bench, FM051_BENCH_PROFILE=1) ==="
    # the bench requires the particle snapshot as its argument (job 13306588
    # omitted it and the stage errored at load; the cp below then shipped a
    # stale fm051_results.csv as "profile" output -- copy only on success now)
    if ( cd "$HOME/FLOWVPM-046" && \
      FM051_BENCH_PROFILE=1 julia --project="$ENVDIR" -t "$THREADS" \
          scripts/fm051_rect_bench.jl \
          "${VPM051_BIN:-$HOME/FLOWVPM-046/data/fm049/p018_710_particles.bin}" ); then
        cp -f "$HOME/FLOWVPM-046/fm051_results.csv" \
            benchmark/results/fm051_profile_results.csv
    else
        echo "WARNING: profile stage failed (status $?) -- continuing to the parity harness"
    fi
fi

# ---- stage 0c: FastMultipole rect device testsets (combined-job rule) -------
# Job 13309844 attribution: seam(:cuda) vs exact 1.235e-2 while cpu-fmm vs
# exact was 1.1e-10 -- the deviation is in the CUDA arm. The rect test file now
# includes an ON-SURFACE p018-scale device-vs-host parity case (the regime the
# old well-separated sweeps never covered). Running it here answers whether the
# defect lives in the FastMultipole device kernel or the seam orchestration.
FM051_DEVTESTS="${FM051_DEVTESTS:-1}"
if [[ "$FM051_DEVTESTS" == "1" && -f "$HOME/FastMultipole-046/test/direct_rectangular_test.jl" ]]; then
    echo "=== stage 0c: FastMultipole rect device testsets (incl. on-surface p018-scale) ==="
    ( cd "$HOME/FastMultipole-046" && \
      julia --project="$ENVDIR" -t "$THREADS" test/direct_rectangular_test.jl ) \
      || echo "WARNING: rect device testsets FAILED (status $?) -- continuing; treat as a FINDING"
fi

echo "FM051 stage 2/3 pass parity -- mode $FM051_MODE, seam $FM051_SEAM"
echo "  repo:$PWD  env:$ENVDIR  threads:$THREADS  host:$(hostname)  job:${SLURM_JOB_ID:-none}"
echo "  restart run:$BENCH_RESTART_RUN  restart step:$RESTART_STEP"
echo "  mesh:$RHPC_MESH formulation:$RHPC_FORMULATION RPM:$RPM NT:$NT nwakerows:$NWAKEROWS settle:$SETTLE_REVS"
echo "  started $(date '+%F %T')"

status=0
julia --project="$ENVDIR" -t "$THREADS" benchmark/fm051_pass_parity.jl || status=$?

echo "Finished ($(date '+%F %T')) with status $status."
echo "Results: benchmark/results/fm051_pass_parity_${FM051_MODE}.csv"
exit $status
