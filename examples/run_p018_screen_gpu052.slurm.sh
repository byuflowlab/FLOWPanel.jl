#!/usr/bin/env bash
#SBATCH --job-name=fp-052-scr-gpu
#SBATCH --nodes=1
#SBATCH --partition=m13h
#SBATCH --gres=gpu:h200:1
#SBATCH --cpus-per-task=64
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# GPU launcher for SCREEN cases on the item-052 NEW stack (cross-pass stages
# A-F + LineGauss + leak fixes; 2026-08-29). Same wrapper pattern as the 018
# silo's run_p018_screen_gpu.slurm.sh, but runs from the 052 silo trees
# (~/F*-052-<arch>) against the 052 depot project ~/fm052env-<arch>, whose
# Manifest dev-paths point at those silos. Case table stays single-sourced in
# this silo's run_p018_screen_hpc.slurm.sh.
#
#   sbatch [resource overrides] examples/run_p018_screen_gpu052.slurm.sh <arch> <case_tag>
#
# arch: h200 (default; in-file #SBATCH header matches m13h) | h100 | gh200.
# H200 via eng: override --partition=eng --qos=eng on the sbatch line
# (052 eng+m13h parallel-submit pattern; scancel the loser once one starts).
# NOTE: the silo dispatcher pins FLOWPANEL_FILAMENT_REG=vatistas by default —
# pass FLOWPANEL_FILAMENT_REG=linegauss via --export for LineGauss runs.
#
# P018_JULIA_TIMEOUT_S (default: wall minus 10 min) bounds the julia run so
# the GPU-path gate below always executes before Slurm kills the job; exit
# 124 from timeout is treated as a valid partial run (probe semantics).

set -euo pipefail

ARCH="${1:-h200}"
CASE="${2:-}"
[[ -n "$CASE" ]] || { echo "ERROR: usage: ... run_p018_screen_gpu052.slurm.sh <arch> <case_tag>" >&2; exit 2; }

case "$ARCH" in
  gh200)
    # ARM node: no x86 module tree; CUDA from CUDA.jl artifacts + node driver.
    export P018_JULIA="${P018_JULIA_OVERRIDE:-$HOME/julia/julia-1.11.7/bin/julia}"
    export JULIA_DEPOT_PATH="$HOME/fm052depot-gh200"
    ;;
  h200|h100)
    module load cuda julia/1.11.7-6bmogfl
    export P018_JULIA="${P018_JULIA_OVERRIDE:-$(command -v julia)}"
    ;;
  *) echo "ERROR: unknown arch '$ARCH' (gh200|h200|h100)" >&2; exit 2 ;;
esac

export P018_REPO="${P018_REPO_OVERRIDE:-$HOME/projects_unified/FLOWPanel.jl}"
export P018_PROJECT="${P018_PROJECT_OVERRIDE:-$HOME/projects_unified/envs/$(uname -m)}"
export P018_THREADS="${SLURM_CPUS_PER_TASK:-64}"

[[ -x "$P018_JULIA" ]] || { echo "ERROR: julia not found at $P018_JULIA" >&2; exit 2; }
[[ -d "$P018_REPO" ]]  || { echo "ERROR: silo tree $P018_REPO missing" >&2; exit 2; }
[[ -d "$P018_PROJECT" ]] || { echo "ERROR: depot project $P018_PROJECT missing" >&2; exit 2; }

# GPU env: source the maintained 052 tuning bundle from the unified tree.
FM052_COMMON="$HOME/projects_unified/FLOWVPM.jl/scripts/fm052_common.sh"
if [[ -f "$FM052_COMMON" ]]; then
  # shellcheck source=/dev/null
  source "$FM052_COMMON"
  if declare -p FM052_GPU_ENV >/dev/null 2>&1; then
    for kv in "${FM052_GPU_ENV[@]}"; do export "$kv"; done
  fi
fi
# screen-campaign hook: let sbatch env override the fm052 bundle reserve
[[ -n "${SCR_GPU_RESERVE_GIB:-}" ]] && export RHPC_SOLVER_S_GPU_RESERVE_GIB="$SCR_GPU_RESERVE_GIB"
export VPM_ARRAYTYPE="${VPM_ARRAYTYPE:-cuarray}"
export FLOWPANEL_GPU_INFLUENCE="${FLOWPANEL_GPU_INFLUENCE:-cuda}"
export RHPC_SOLVER_S="${RHPC_SOLVER_S:-true}"   # driver refuses S_GPU without S
export RHPC_SOLVER_S_GPU="${RHPC_SOLVER_S_GPU:-true}"
export FASTMULTIPOLE_FORCE_CUDA_LOAD=1

echo "=== 052 GPU launcher: arch=$ARCH case=$CASE filament_reg=${FLOWPANEL_FILAMENT_REG:-<src default>} ==="
echo "  julia:$P018_JULIA project:$P018_PROJECT repo:$P018_REPO threads:$P018_THREADS"
nvidia-smi -L

cd "$P018_REPO"
mkdir -p logs/slurm

# Bound the julia run so the gate always runs before the wall (10 min margin).
WALL_S=$(( $(squeue -h -j "$SLURM_JOB_ID" -O TimeLimit | awk -F'[-:]' \
  'NF==4{print (($1*24+$2)*60+$3)*60+$4} NF==3{print ($1*60+$2)*60+$3} NF==2{print $1*60+$2}') ))
TIMEOUT_S="${P018_JULIA_TIMEOUT_S:-$(( WALL_S > 1200 ? WALL_S - 600 : WALL_S/2 ))}"
echo "  julia timeout: ${TIMEOUT_S}s (wall ${WALL_S}s)"

# Run the CPU dispatcher in a subshell wrapped by `timeout`; it inherits the
# P018_*/GPU env and the case table stays single-sourced. Dispatcher exit 124
# = internal timeout = acceptable partial run.
set +e
timeout --signal=TERM --kill-after=60 "$TIMEOUT_S" \
  bash -c 'source examples/run_p018_screen_hpc.slurm.sh "$1"' _ "$CASE"
RC=$?
set -e
if [[ $RC -eq 124 ]]; then
  echo "NOTE: julia run hit internal timeout ${TIMEOUT_S}s -- gating on partial log (probe semantics)."
elif [[ $RC -ne 0 ]]; then
  echo "ERROR: dispatcher exited rc=$RC" >&2
fi

# ---- GPU-path gate (052c pattern) -------------------------------------------
OUTFILE=$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *StdOut=//p')
ERRFILE=$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *StdErr=//p')
[[ "$ERRFILE" == "$OUTFILE" ]] && ERRFILE=/dev/null
GPU_N=$(cat "$OUTFILE" "$ERRFILE" | grep -c source_influence_s_gpu_gemv || true)
CPU_N=$(cat "$OUTFILE" "$ERRFILE" | grep -c 'source_influence_s_gemv'   || true)
# Word-match capital NaN only; exclude the per-step CT table and the
# "Bernoulli vs KJ:" summary (NaN-by-design KJ columns when RUN_KJ=false).
# Also exclude "plateau mean=" diagnostics: an empty settle window (screen
# probes) printed NaN there and failed clean run 13542825 on 2026-08-31; the
# drivers now print "n/a" but this keeps old driver output from tripping.
NAN_N=$(cat "$OUTFILE" "$ERRFILE" | grep -w 'NaN' | grep -v ' | ' | grep -v 'vs KJ:' | grep -vc 'plateau mean=' || true)
echo "GATE: gpu_gemv=$GPU_N cpu_gemv=$CPU_N nan_lines=$NAN_N dispatcher_rc=$RC"
GATE_RC=0
[[ $GPU_N -gt 0 ]] || { echo "GATE FAIL: GPU source-influence path never ran" >&2; GATE_RC=1; }
[[ $CPU_N -eq 0 ]] || { echo "GATE FAIL: CPU source-influence path ran $CPU_N times" >&2; GATE_RC=1; }
[[ $NAN_N -eq 0 ]] || { echo "GATE FAIL: NaN in log" >&2; GATE_RC=1; }
[[ $RC -eq 0 || $RC -eq 124 ]] || GATE_RC=1
echo "=== 052 GPU launcher done: gate_rc=$GATE_RC ==="
exit $GATE_RC
