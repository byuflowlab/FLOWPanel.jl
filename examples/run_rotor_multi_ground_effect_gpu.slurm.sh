#!/usr/bin/env bash
#SBATCH --job-name=fp-022g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=02:15:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# BRAINSTORM/022 Phase-6 GPU carrier owned by 052b; extended 2026-08-31 for the
# Phase-5 h/R sweep (p022g_hr* arms, production mode, unified repo per HPC.md
# ruling). GPU placement is supplied by the submit wrapper.
set -euo pipefail

THREADS=16
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
PROJECT="${P022G_PROJECT_OVERRIDE:-/home/rander39/projects/envs/$(uname -m)}"
CASE="${1:-}"
MODE="${P022G_MODE:-smoke}"
[[ "$PWD" == "$EXPECTED_REPO" || "${P022G_SETUP_ONLY:-0}" == 1 ]] || {
  echo "ERROR: run from $EXPECTED_REPO (current: $PWD)" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# linegauss per Ryan's 2026-08-31 sweep ruling (CPU anchors ran vatistas —
# expect a small family-systematic; see plans/p022_hr_sweep_gpu_20260831/)
# P022G_VPM_ARRAYTYPE=array falls back to a host-resident particle field
# (GPU still runs panel influence): diagnostic escape hatch for the missing
# device-resident source_to_buffer! overload (2026-08-31 smoke finding).
export FLOWPANEL_FILAMENT_REG=linegauss VPM_ARRAYTYPE="${P022G_VPM_ARRAYTYPE:-cuarray}"
export FASTMULTIPOLE_FORCE_CUDA_LOAD=1
export FLOWPANEL_GPU_INFLUENCE=cuda GPU_ALLOW_FALLBACK=false
export FLOWPANEL_STEP_TIMERS=1 FLOWPANEL_GPU_TIMERS=1
# Sharing is intentionally off until a separate direct-matrix audit records
# parity for the +1 and -1 handedness classes. A validated production mapping
# is ROTOR_OPERATOR_OWNERS=1,2,1,2 with ROTOR_OPERATOR_CLASSES_VALIDATED=true.
unset SHARE_ROTOR_OPERATOR ROTOR_OPERATOR_OWNERS ROTOR_OPERATOR_CLASSES_VALIDATED
export RPM=6000 RHO=1.16 ROTOR_R=0.1195 NT=36 RHPC_FORMULATION=velocity
export RHPC_MESH=45_185_ct4 RELAX_RLXF=0.3 PARTICLE_SHEDDING=sigma_overlap
export OVERLAP=2.75 P_PER_STEP=12 MERGE_R_FACTOR=0.0055 NWAKEROWS=1
export DAS_UNIFORM_DSIGMA=3.4 CORE_SPREADING_ACTIVE=true WAKE_CORE_BETA=1e9
export BERNOULLI_ONLY=true RUN_KJ=false SAVE_VTK=true
export FMM_BODY_EXPANSION_ORDER=${FMM_BODY_EXPANSION_ORDER:-17} FMM_BODY_ACCEPTANCE=${FMM_BODY_ACCEPTANCE:-0.7} FMM_BODY_LEAF_SIZE=${FMM_BODY_LEAF_SIZE:-109}
export FMM_WAKE_EXPANSION_ORDER=${FMM_WAKE_EXPANSION_ORDER:-16} FMM_WAKE_ACCEPTANCE=${FMM_WAKE_ACCEPTANCE:-0.6} FMM_WAKE_LEAF_SIZE=${FMM_WAKE_LEAF_SIZE:-38}
export SPINUP_REVS=1.5 SPINUP_START_FRACTION=0.4
export MAGVINF_START=0.0 MAGVINF_PEAK=5.0 MAGVINF_END=0.0
case "$MODE" in
  smoke)
    export SPINUP_REVS=0.1 FREESTREAM_RAMP_REVS=0.05 FREESTREAM_HOLD_REVS=0.05
    export FREESTREAM_WITHDRAW_REVS=0.05 SETTLE_REVS=0.10 NREVS=0.35
    ;;
  probe)
    export SPINUP_REVS=0.25 FREESTREAM_RAMP_REVS=0.25 FREESTREAM_HOLD_REVS=0.25
    export FREESTREAM_WITHDRAW_REVS=0.25 SETTLE_REVS=1.0 NREVS=2.0
    ;;
  accept)
    [[ "${P022G_CONFIRM_ACCEPTANCE:-}" == "YES" ]] || {
      echo "ERROR: full acceptance requires P022G_CONFIRM_ACCEPTANCE=YES" >&2; exit 2; }
    export SPINUP_REVS=1.5 FREESTREAM_RAMP_REVS=1.0 FREESTREAM_HOLD_REVS=1.5
    export FREESTREAM_WITHDRAW_REVS=4.0 SETTLE_REVS=3.5 NREVS=10
    ;;
  production)
    # Full p022 fine schedule (1007 steps ~ 28 revs), matching CPU anchor
    # job 13207681 (examples/run_rotor_ground_effect_hpc.slurm.sh).
    [[ "${P022G_CONFIRM_PRODUCTION:-}" == "YES" ]] || {
      echo "ERROR: production requires P022G_CONFIRM_PRODUCTION=YES" >&2; exit 2; }
    export SPINUP_REVS=1.5 FREESTREAM_RAMP_REVS=1.0 FREESTREAM_HOLD_REVS=1.5
    export FREESTREAM_WITHDRAW_REVS=4.0 SETTLE_REVS=20 NREVS=10
    ;;
  *) echo "ERROR: P022G_MODE must be smoke, probe, accept, or production" >&2; exit 2 ;;
esac
export CONVERGENCE_REVS=10 CONVERGENCE_MEAN_TOL=0.005 CONVERGENCE_PTP_TOL=0.02
# GS knobs are default-guarded so sbatch --export=ALL,GS_...=... can override
# them (they were previously hard-exported here, silently clobbering sweeps).
export GS_LOG="${GS_LOG:-true}" GS_MAX_OUTER="${GS_MAX_OUTER:-50}" GS_TOL="${GS_TOL:-1e-8}"
export GS_VERBOSE="${GS_VERBOSE:-false}" ROTOR_SPACING_R=2.7
export TRUNCATION_DEPTH_R=4.5

case "$CASE" in
  p022g_1r_oge) export NROTORS=1 GROUND_ENABLE=false TRUNC_RADIUS_R=1.5 ;;
  p022g_2r_oge) export NROTORS=2 GROUND_ENABLE=false TRUNC_RADIUS_R=1.5 ;;
  p022g_4r_oge) export NROTORS=4 GROUND_ENABLE=false TRUNC_RADIUS_R=1.5 ;;
  p022g_1r_ige) export NROTORS=1 GROUND_ENABLE=true GROUND_H_R=1.5 TRUNC_RADIUS_R=3.0
                  export GROUND_RADIUS_R=4.0 GROUND_PANEL_LENGTH_R=0.15
                  export GROUND_PARTICLE_POLICY=none GROUND_DAMP_BAND_R=0.1 ;;
  p022g_2r_ige) export NROTORS=2 GROUND_ENABLE=true GROUND_H_R=1.5 TRUNC_RADIUS_R=3.0
                  export GROUND_RADIUS_R=4.0 GROUND_PANEL_LENGTH_R=0.15
                  export GROUND_PARTICLE_POLICY=none GROUND_DAMP_BAND_R=0.1 ;;
  p022g_4r_ige) export NROTORS=4 GROUND_ENABLE=true GROUND_H_R=1.5 TRUNC_RADIUS_R=3.0
                  export GROUND_RADIUS_R=4.0 GROUND_PANEL_LENGTH_R=0.15
                  export GROUND_PARTICLE_POLICY=none GROUND_DAMP_BAND_R=0.1 ;;
  # ---- Phase-5 h/R sweep arms (2026-08-31): single rotor, p022_ige_fine
  # config (job 13207681). Every defining knob exported UNCONDITIONALLY per
  # the 018 ops rule. TRUNCATION_DEPTH_R = h/R + 3 (022 convention).
  p022g_hr05) export NROTORS=1 GROUND_ENABLE=true GROUND_H_R=0.5 TRUNC_RADIUS_R=3.0
              export GROUND_RADIUS_R=4.0 GROUND_PANEL_LENGTH_R=0.15
              export GROUND_PARTICLE_POLICY=none GROUND_DAMP_BAND_R=0
              export TRUNCATION_DEPTH_R=3.5 ;;
  p022g_hr10) export NROTORS=1 GROUND_ENABLE=true GROUND_H_R=1.0 TRUNC_RADIUS_R=3.0
              export GROUND_RADIUS_R=4.0 GROUND_PANEL_LENGTH_R=0.15
              export GROUND_PARTICLE_POLICY=none GROUND_DAMP_BAND_R=0
              export TRUNCATION_DEPTH_R=4 ;;
  p022g_hr15) export NROTORS=1 GROUND_ENABLE=true GROUND_H_R=1.5 TRUNC_RADIUS_R=3.0
              export GROUND_RADIUS_R=4.0 GROUND_PANEL_LENGTH_R=0.15
              export GROUND_PARTICLE_POLICY=none GROUND_DAMP_BAND_R=0
              export TRUNCATION_DEPTH_R=4.5 ;;
  p022g_hr20) export NROTORS=1 GROUND_ENABLE=true GROUND_H_R=2.0 TRUNC_RADIUS_R=3.0
              export GROUND_RADIUS_R=4.0 GROUND_PANEL_LENGTH_R=0.15
              export GROUND_PARTICLE_POLICY=none GROUND_DAMP_BAND_R=0
              export TRUNCATION_DEPTH_R=5 ;;
  *) echo "ERROR: unknown 052b case '$CASE'" >&2; exit 2 ;;
esac

export RUN_NAME="${P022_RUN_NAME:-$CASE}"
echo "052b GPU mode=$MODE case=$CASE run=$RUN_NAME nrotors=$NROTORS ground=$GROUND_ENABLE"
echo "filament_reg=$FLOWPANEL_FILAMENT_REG ground_h_r=${GROUND_H_R:-n/a} damp_band=${GROUND_DAMP_BAND_R:-n/a} particle_policy=${GROUND_PARTICLE_POLICY:-n/a} trunc_depth=$TRUNCATION_DEPTH_R"
if [[ "$MODE" == accept ]]; then
  echo "mesh=$RHPC_MESH schedule=1.0+1.5+4.0+3.5 acceptance_steps=360 spinup_steps=54 total_steps=414"
elif [[ "$MODE" == production ]]; then
  echo "mesh=$RHPC_MESH schedule=1.0+1.5+4.0+20+10 spinup=1.5 total_steps=1007"
fi
echo "array=$VPM_ARRAYTYPE influence=$FLOWPANEL_GPU_INFLUENCE fallback=$GPU_ALLOW_FALLBACK operator_owners=fresh no_gpu_S=true"
echo "body_fmm=$FMM_BODY_EXPANSION_ORDER/$FMM_BODY_ACCEPTANCE/$FMM_BODY_LEAF_SIZE wake_fmm=$FMM_WAKE_EXPANSION_ORDER/$FMM_WAKE_ACCEPTANCE/$FMM_WAKE_LEAF_SIZE"
[[ "${P022G_SETUP_ONLY:-0}" == 1 ]] && exit 0

# Unified env, arch-dispatched (2026-09-01, P018 gpu052 launcher pattern):
# x86 uses the julia 1.11.7 module; ARM (gh200 / mgh nodes) uses the aarch64
# tarball julia + the ARM depot, with CUDA supplied by CUDA.jl artifacts and
# the node driver (no x86 module tree on ARM).
case "$(uname -m)" in
  aarch64)
    JULIA_BIN="${P022G_JULIA_OVERRIDE:-$HOME/julia/julia-1.11.7/bin/julia}"
    export JULIA_DEPOT_PATH="$HOME/fm052depot-gh200"
    ;;
  *)
    module load cuda julia/1.11.7-6bmogfl
    JULIA_BIN="${P022G_JULIA_OVERRIDE:-$(command -v julia)}"
    ;;
esac
[[ -x "$JULIA_BIN" ]] || { echo "ERROR: julia not found at $JULIA_BIN" >&2; exit 2; }
JULIA_VERSION="$("$JULIA_BIN" --version)"
[[ "$JULIA_VERSION" == "julia version 1.11.7" ]] || {
  echo "ERROR: required Julia 1.11.7, got '$JULIA_VERSION'" >&2; exit 2; }
[[ -d "$PROJECT" ]] || { echo "ERROR: project env $PROJECT missing" >&2; exit 2; }
command -v nvidia-smi >/dev/null || { echo "ERROR: nvidia-smi unavailable" >&2; exit 2; }
export P022G_GPU_MODEL="$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
[[ -z "${P022G_REQUIRED_GPU_MODEL:-}" || "$P022G_GPU_MODEL" == *"$P022G_REQUIRED_GPU_MODEL"* ]] || {
  echo "ERROR: GPU '$P022G_GPU_MODEL' does not match required '$P022G_REQUIRED_GPU_MODEL'" >&2; exit 2; }
export P022G_FLOWPANEL_REV="$(git rev-parse HEAD)"
export P022G_FASTMULTIPOLE_REV="$(git -C ../FastMultipole rev-parse HEAD 2>/dev/null || echo unknown)"
export P022G_FLOWVPM_REV="$(git -C ../FLOWVPM.jl rev-parse HEAD 2>/dev/null || echo unknown)"
echo "revisions FLOWPanel=$P022G_FLOWPANEL_REV FastMultipole=$P022G_FASTMULTIPOLE_REV FLOWVPM=$P022G_FLOWVPM_REV"
echo "gpu=$P022G_GPU_MODEL"
"$JULIA_BIN" --project="$PROJECT" -e 'import FastMultipole; FastMultipole.load_cuda_radix_lifecycle!() || error(FastMultipole.cuda_radix_status()); C=getglobal(FastMultipole,:CUDA); Base.invokelatest(C.functional) || error("CUDA.functional() is false")'

if [[ -d "data/$RUN_NAME" ]]; then
  [[ "${P022G_EXISTING_RESULT:-error}" == preserve ]] || {
    echo "ERROR: data/$RUN_NAME exists; set P022G_EXISTING_RESULT=preserve to move it" >&2; exit 2; }
  previous="data/${RUN_NAME}.prev.$(date +%Y%m%dT%H%M%S)"
  mv "data/$RUN_NAME" "$previous"
  echo "preserved previous result at $previous"
fi

gpu_memory_log="${TMPDIR:-/tmp}/p022g_gpu_memory_${SLURM_JOB_ID:-$$}.csv"
host_memory_log="${TMPDIR:-/tmp}/p022g_host_memory_${SLURM_JOB_ID:-$$}.txt"
nvidia-smi --query-gpu=timestamp,memory.used --format=csv,noheader,nounits -l 1 >"$gpu_memory_log" &
gpu_monitor_pid=$!
cleanup_gpu_monitor() {
  kill "$gpu_monitor_pid" 2>/dev/null || true
  wait "$gpu_monitor_pid" 2>/dev/null || true
  if [[ -d "data/$RUN_NAME" ]]; then
    mkdir -p "data/$RUN_NAME/run_artifacts"
    cp -f "$gpu_memory_log" "data/$RUN_NAME/run_artifacts/" 2>/dev/null || true
    cp -f "$host_memory_log" "data/$RUN_NAME/run_artifacts/" 2>/dev/null || true
  fi
}
trap cleanup_gpu_monitor EXIT

case_start_epoch="$(date +%s)"
/usr/bin/time -v -o "$host_memory_log" "$JULIA_BIN" --project="$PROJECT" -t "$THREADS" examples/rotor_hover_ground_effect.jl
case_elapsed_s="$(( $(date +%s) - case_start_epoch ))"
cleanup_gpu_monitor
trap - EXIT
peak_device_mib="$(awk -F, 'BEGIN{m=0}{gsub(/ /,"",$2); if ($2+0>m)m=$2+0}END{print m}' "$gpu_memory_log")"
peak_host_kib="$(awk -F: '/Maximum resident set size/{gsub(/ /,"",$2); print $2}' "$host_memory_log")"
device_total_mib="$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1 | tr -d ' ')"
device_reserve_fraction="$(awk -v p="$peak_device_mib" -v t="$device_total_mib" 'BEGIN{print (t-p)/t}')"
metadata="data/$RUN_NAME/${RUN_NAME}_case_metadata.toml"
if [[ -f "$metadata" ]]; then
  {
    echo "flowpanel_revision = \"$P022G_FLOWPANEL_REV\""
    echo "fastmultipole_revision = \"$P022G_FASTMULTIPOLE_REV\""
    echo "flowvpm_revision = \"$P022G_FLOWVPM_REV\""
    echo "gpu_model = \"$P022G_GPU_MODEL\""
    echo "peak_device_memory_mib = $peak_device_mib"
    echo "peak_host_memory_kib = ${peak_host_kib:-0}"
    echo "device_memory_reserve_fraction = $device_reserve_fraction"
    echo "authoritative_case_elapsed_s = $case_elapsed_s"
    echo "p022g_mode = \"$MODE\""
  } >>"$metadata"
fi
if [[ "$MODE" == accept ]]; then
  (( case_elapsed_s <= 7200 )) || { echo "ERROR: case elapsed $case_elapsed_s s exceeds 7200 s" >&2; exit 1; }
  awk -v r="$device_reserve_fraction" 'BEGIN{exit !(r >= 0.20)}' || {
    echo "ERROR: device memory reserve $device_reserve_fraction is below 20%" >&2; exit 1; }
fi
if [[ "$MODE" == accept || "$MODE" == production ]]; then
  "$JULIA_BIN" -e 'using TOML; m=TOML.parsefile(ARGS[1]); get(m,"all_finite",false) || error("nonfinite case"); get(m,"gs_nonconverged",1)==0 || error("nonconverged block solve"); get(m,"gs_iters_max",typemax(Int))<=50 || error("GS cap exceeded"); get(m,"block_gs_normalized_residual",Inf)<=1e-8 || error("final normalized residual failed"); get(m,"gpu_route_fallback_total",1)==0 || error("GPU fallback recorded"); get(m,"gpu_route_total_hits",0)>0 || error("no GPU route hits")' "$metadata"
fi
echo "Case $CASE finished; peak device ${peak_device_mib} MiB; peak host ${peak_host_kib:-unknown} KiB"
