#!/usr/bin/env bash
#SBATCH --job-name=fp-018-gpu
#SBATCH --nodes=1
#SBATCH --partition=mgh
#SBATCH --gres=gpu:gh200:1
#SBATCH --cpus-per-task=72
#SBATCH --mem=192G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# GPU launcher for the 018 campaign (BRAINSTORM/018): wraps the CPU case
# dispatcher run_dji9443_hover_ct_hpc.slurm.sh with the item-052 GPU stack so
# any p018_* case runs on a GPU node. The case table stays single-sourced in
# the CPU script; this payload only selects arch, env, and tree.
#
#   sbatch [resource overrides] examples/run_dji9443_hover_ct_gpu.slurm.sh <arch> <case_tag>
#
# arch: gh200 (default; in-file #SBATCH header matches) | h200 | h100.
# No constraint in the header: mgh is all-ARM so it's redundant there, and a
# hardcoded arm constraint rejects x86 partitions (eng refused exactly this on
# 2026-09-05). For non-gh200 arches override partition/gres/cpus on the sbatch
# command line (052 submit-wrapper style), e.g. H200 via eng:
#   sbatch --partition=eng --qos=eng --gres=gpu:h200:1 \
#          --cpus-per-task=64 ... examples/run_dji9443_hover_ct_gpu.slurm.sh h200 p018_final
#
# Runs from the per-arch silo tree ~/FLOWPanel-018-gpu-<arch> against the 052
# shared depot project ~/fm052env-<arch> (GPU seam + precompiled CUDA stack).
# Julia is pinned to 1.11.7 (1.12.6 segfaults in host LLVM JIT with the GPU
# stack); on gh200 that is the manual aarch64 install, the module is x86-only.
#
# P018_JULIA_TIMEOUT_S (default: wall minus 10 min) bounds the julia run so the
# GPU-path gate below always executes before Slurm kills the job; exit 124 from
# timeout is treated as a valid partial run (probe semantics, 052 pattern).

set -euo pipefail

ARCH="${1:-gh200}"
CASE="${2:-}"
[[ -n "$CASE" ]] || { echo "ERROR: usage: ... run_dji9443_hover_ct_gpu.slurm.sh <arch> <case_tag>" >&2; exit 2; }

case "$ARCH" in
  gh200)
    # No module load: the ARM node has no x86 module tree; CUDA comes from
    # CUDA.jl artifacts in the ARM depot + the node driver (052a recipe).
    export P018_JULIA="${P018_JULIA_OVERRIDE:-$HOME/julia/julia-1.11.7/bin/julia}"   # aarch64 manual install (052a)
    export JULIA_DEPOT_PATH="$HOME/fm052depot-gh200"          # aarch64 depot (052a)
    ;;
  h200|h100)
    module load cuda julia/1.11.7-6bmogfl
    export P018_JULIA="${P018_JULIA_OVERRIDE:-$(command -v julia)}"
    ;;
  *) echo "ERROR: unknown arch '$ARCH' (gh200|h200|h100)" >&2; exit 2 ;;
esac

# P018_REPO_OVERRIDE / P018_PROJECT_OVERRIDE (Ryan order 2026-08-29): point a
# submission at an alternate silo tree + env (e.g. the -lg LineGauss snapshot)
# without disturbing the default silos mid-ladder.
export P018_REPO="${P018_REPO_OVERRIDE:-$HOME/FLOWPanel-018-gpu-$ARCH}"
# p018env-<arch> = copy of fm052env-<arch> with FLOWPanel dev-path repointed
# at the 018 silo tree (the 052 envs/trees stay untouched for 052c).
export P018_PROJECT="${P018_PROJECT_OVERRIDE:-$HOME/p018env-$ARCH}"
export P018_THREADS="${SLURM_CPUS_PER_TASK:-64}"

[[ -x "$P018_JULIA" ]] || { echo "ERROR: julia not found at $P018_JULIA" >&2; exit 2; }
[[ -d "$P018_REPO" ]]  || { echo "ERROR: silo tree $P018_REPO missing" >&2; exit 2; }
[[ -d "$P018_PROJECT" ]] || { echo "ERROR: depot project $P018_PROJECT missing" >&2; exit 2; }

# GPU env: source the maintained 052 tuning bundle from the unified tree
# (legacy silos are slated for deletion); fall back to inline exports.
FM052_COMMON="$HOME/projects/FLOWVPM.jl/scripts/fm052_common.sh"
if [[ -f "$FM052_COMMON" ]]; then
  # shellcheck source=/dev/null
  source "$FM052_COMMON"
  if declare -p FM052_GPU_ENV >/dev/null 2>&1; then
    for kv in "${FM052_GPU_ENV[@]}"; do export "$kv"; done
  fi
fi
export VPM_ARRAYTYPE="${VPM_ARRAYTYPE:-cuarray}"
export FLOWPANEL_GPU_INFLUENCE="${FLOWPANEL_GPU_INFLUENCE:-cuda}"
export RHPC_SOLVER_S="${RHPC_SOLVER_S:-true}"   # driver refuses S_GPU without S
export RHPC_SOLVER_S_GPU="${RHPC_SOLVER_S_GPU:-true}"
export FASTMULTIPOLE_FORCE_CUDA_LOAD=1

echo "=== 018 GPU launcher: arch=$ARCH case=$CASE ==="
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
  bash -c 'source examples/run_dji9443_hover_ct_hpc.slurm.sh "$1"' _ "$CASE"
RC=$?
set -e
if [[ $RC -eq 124 ]]; then
  echo "NOTE: julia run hit internal timeout ${TIMEOUT_S}s -- gating on partial log (probe semantics)."
elif [[ $RC -ne 0 ]]; then
  echo "ERROR: dispatcher exited rc=$RC" >&2
fi

# ---- GPU-path gate (052c pattern) -------------------------------------------
# The gemv markers are julia @info lines and land on STDERR; scan both streams
# (they are separate files when submitted with --output and --error).
OUTFILE=$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *StdOut=//p')
ERRFILE=$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *StdErr=//p')
[[ "$ERRFILE" == "$OUTFILE" ]] && ERRFILE=/dev/null
GPU_N=$(cat "$OUTFILE" "$ERRFILE" | grep -c source_influence_s_gpu_gemv || true)
CPU_N=$(cat "$OUTFILE" "$ERRFILE" | grep -c 'source_influence_s_gemv'   || true)
# Word-match capital NaN only: the dispatcher banner legitimately prints
# lowercase ":nan" for unset das_* knobs. Exclude the per-step CT comparison
# table (" | "-separated) AND the driver's "Bernoulli vs KJ:" summary line —
# both carry KJ columns that are NaN by design when RUN_KJ=false. (2026-08-27:
# the summary line is not " | "-separated, so it slipped the table filter and
# failed the gate on three otherwise-clean production runs.) Also exclude the
# "plateau mean=" diagnostics: an empty settle window (screen probes) printed
# NaN there and failed clean run 13542825 on 2026-08-31; the drivers now print
# "n/a" but this keeps old driver output from tripping the gate.
NAN_N=$(cat "$OUTFILE" "$ERRFILE" | grep -w 'NaN' | grep -v ' | ' | grep -v 'vs KJ:' | grep -vc 'plateau mean=' || true)
echo "GATE: gpu_gemv=$GPU_N cpu_gemv=$CPU_N nan_lines=$NAN_N dispatcher_rc=$RC"
GATE_RC=0
[[ $GPU_N -gt 0 ]] || { echo "GATE FAIL: GPU source-influence path never ran" >&2; GATE_RC=1; }
[[ $CPU_N -eq 0 ]] || { echo "GATE FAIL: CPU source-influence path ran $CPU_N times" >&2; GATE_RC=1; }
[[ $NAN_N -eq 0 ]] || { echo "GATE FAIL: NaN in log" >&2; GATE_RC=1; }
[[ $RC -eq 0 || $RC -eq 124 ]] || GATE_RC=1
echo "=== 018 GPU launcher done: gate_rc=$GATE_RC ==="
exit $GATE_RC
