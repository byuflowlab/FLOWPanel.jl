#!/usr/bin/env bash
#SBATCH --job-name=fp-023-fmmtune
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# BRAINSTORM/023: per-phase timing of production-config timesteps warm-started
# from the retained mature restart set of p018_cs_f1_l3p4 (step ~1034, rev
# ~28.5, ~181k particles). Runs benchmark/p018_mature_wake_timing.jl; results
# land in benchmark/results/, scratch VTK in node tempdir — data/ is untouched
# except read-only restart loads.
#
# Submit from the cluster checkout top level:
#   sbatch benchmark/slurm/p023_mature_timing.sh
# Thread-scaling spot check (optional):
#   sbatch --ntasks=32 --export=ALL,THREADS=32,BENCH_NSTEPS=3,BENCH_PROFILE=false \
#          benchmark/slurm/p023_mature_timing.sh
#
# The env block below is a verbatim copy of the p018_cs_f1_l3p4 recipe from
# examples/run_dji9443_hover_ct_hpc.slurm.sh (fixed block + p018_* group
# defaults + case arm), with SETTLE_REVS=22 (the F1 trio was submitted with
# P018_SETTLE_REVS=22 => 30-rev schedule, t_range of 1081 steps). If the
# launcher recipe changes, re-copy — drift here invalidates the timing.

set -euo pipefail
THREADS="${THREADS:-64}"
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

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
unset RELAX_FILTER_DOWNSTREAM_R
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

# ---- benchmark knobs ----
export BENCH_RESTART_RUN="${BENCH_RESTART_RUN:-p018_cs_f1_l3p4}"
export RESTART_STEP="${RESTART_STEP:--1}"
export BENCH_NSTEPS="${BENCH_NSTEPS:-5}"
export BENCH_WRITE_OUTPUT="${BENCH_WRITE_OUTPUT:-true}"
export BENCH_PROFILE="${BENCH_PROFILE:-true}"

echo "BRAINSTORM/023 mature-wake timing — restart run $BENCH_RESTART_RUN"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"
echo "  mesh:$RHPC_MESH formulation:$RHPC_FORMULATION RPM:$RPM NT:$NT depth:${TRUNCATION_DEPTH_R}R rlxf:$RELAX_RLXF"
echo "  overlap:$OVERLAP pps:$P_PER_STEP merge_r:$MERGE_R_FACTOR nwakerows:$NWAKEROWS settle:$SETTLE_REVS"
echo "  sigma_chord:$SIGMA_CHORD_FRACTION das_lambda:$DAS_SIGMA_LAMBDA das_beta:$DAS_CURVATURE_BETA visc:$CORE_SPREADING_ACTIVE"
echo "  nsteps:$BENCH_NSTEPS write_output:$BENCH_WRITE_OUTPUT profile:$BENCH_PROFILE restart_step:$RESTART_STEP"
echo "  started $(date '+%F %T')"

julia --project=. -t "$THREADS" benchmark/p023_fmm_tune.jl

echo "Finished ($(date '+%F %T')). Results: benchmark/results/p023_fmm_tune.csv + _history.csv"
