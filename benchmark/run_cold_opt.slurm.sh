#!/usr/bin/env bash
#SBATCH --job-name=p021-cold-opt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=06:00:00
#SBATCH --output=logs/slurm/cold-opt-%j.out
#SBATCH --error=logs/slurm/cold-opt-%j.err
# Initialized-FGS optimization campaign driver (see fgs_cold_README.md and
# BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_initialized_cpu_optimization_plan_20260911.md).
# Read BYU_ORC_AGENTS.md before submission. Submit from the pinned worktree.
set -euo pipefail
set +u # Site profile reads optional interactive-shell variables.
source /etc/profile
set -u
# The pinned environment includes CUDA extensions; precompile needs its toolkit
# even though this campaign executes CPU solvers only.
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"
export RUNG=${COLD_OPT_RUNG:-R2} CONFIGS=fgs
screen_set=${COLD_OPT_SCREEN_SET:-inner:1,2,3,5}
stage=${COLD_OPT_STAGE:-all}
[[ "$stage" == all || "$stage" == controls_smoke || "$stage" == screen_profile ]] || exit 2
if [[ "$stage" == screen_profile ]]; then
 : "${CONFIG_FILE:?selected smoke settings required}"
 [[ -f "$CONFIG_FILE" ]] || exit 2
fi
run="$COLD_DATA_ROOT/opt-$SLURM_JOB_ID"
mkdir "$run"
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
command -v ptxas > "$run/ptxas_path.txt"
lscpu > "$run/lscpu.txt"
# Use distinct physical cores from the allocation, keeping every arm on this node.
cpulist=$(python3 - <<'PY'
import os
from pathlib import Path
seen=set(); cores=[]
for c in sorted(os.sched_getaffinity(0)):
 p=Path('/sys/devices/system/cpu')/f'cpu{c}'/'topology'
 key=((p/'physical_package_id').read_text(),(p/'core_id').read_text())
 if key not in seen:
  cores.append(c); seen.add(key)
assert len(cores)>=64, 'fewer than 64 physical cores in allocation'
print(','.join(map(str,cores[:64])))
PY
)
printf '%s\n' "$cpulist" > "$run/cpu_affinity.txt"
run_process() {
 local generation=$1 jt=$2 bt=$3 driver=$4
 export OUTDIR="$run/$generation" BENCH_CASE_ROOT="$run/fixture-$generation"
 taskset -c "$cpulist" bash benchmark/run_cold_process.sh "$jt" "$bt" "$driver" > "$run/$generation.log" 2>&1
}
if [[ "$stage" != screen_profile ]]; then
# Parse without loading packages before spending time on precompilation.
taskset -c "$cpulist" bash benchmark/run_cold_process.sh 1 1 benchmark/cold_parse.jl > "$run/parse.log" 2>&1
# Environment/precompilation process is isolated from tests and timed processes.
taskset -c "$cpulist" bash benchmark/run_cold_process.sh 4 1 benchmark/cold_precompile.jl > "$run/precompile.log" 2>&1
run_process controls-j1 1 1 test/runtests_benchmark_cold.jl
run_process controls-j4 4 1 test/runtests_benchmark_cold.jl
export STAGE=baseline
unset CONFIG_FILE
run_process smoke-j4-b1 4 1 benchmark/rotor_hover_solver_cold_smoke.jl
export CONFIG_FILE="$run/smoke-j4-b1/selected.toml"
if [[ "$stage" == controls_smoke ]]; then
 printf 'controls and smoke completed\n' > "$run/COMPLETED"
 exit 0
fi
fi
sha256sum "$CONFIG_FILE" > "$run/selected.sha256"
# Prepared-only seed baselines at BLAS=1: the plan excludes construction.
export COLD_PREPARED_ONLY=1
export STAGE=verify
run_process baseline-j4-b1 4 1 benchmark/rotor_hover_solver_cold.jl
run_process baseline-j64-b1 64 1 benchmark/rotor_hover_solver_cold.jl
# Step-1 screen: explicit roster, per-candidate tolerance recalibration; a
# failed roster point is recorded and skipped inside the process.
(
 export STAGE=screen SCREEN_SET="$screen_set"
 unset CONFIG_FILE
 run_process screen-j64-b1 64 1 benchmark/rotor_hover_solver_cold.jl
)
# Accumulated CPU profile of the calibrated seed for denser attribution.
export STAGE=verify COLD_INVESTIGATION=1 COLD_PROFILE_REPS=10
run_process profile-fgs-j64-b1 64 1 benchmark/rotor_hover_solver_phase2_profile.jl
sha256sum -c "$run/selected.sha256"
printf 'completed\n' > "$run/COMPLETED"
