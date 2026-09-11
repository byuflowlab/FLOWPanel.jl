#!/usr/bin/env bash
#SBATCH --job-name=p021-cold-pilot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/cold-pilot-%j.out
#SBATCH --error=logs/slurm/cold-pilot-%j.err
# Read BYU_ORC_AGENTS.md before submission. Submit from this pinned worktree.
set -euo pipefail
set +u # Site profile reads optional interactive-shell variables.
source /etc/profile
set -u
module load julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"
export RUNG=R1 CONFIGS=fgs:krylov_ilu
stage=${COLD_PILOT_STAGE:-all}
[[ "$stage" == all || "$stage" == controls_smoke || "$stage" == timing_profiles ]] || exit 2
if [[ "$stage" == timing_profiles ]]; then
 : "${CONFIG_FILE:?selected smoke settings required}"
 [[ -f "$CONFIG_FILE" ]] || exit 2
fi
pilot="$COLD_DATA_ROOT/pilot-$SLURM_JOB_ID"
mkdir "$pilot"
cp "$CAMPAIGN_PINS" "$pilot/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$pilot/Manifest.toml"
lscpu > "$pilot/lscpu.txt"
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
printf '%s\n' "$cpulist" > "$pilot/cpu_affinity.txt"
run_process() {
 local generation=$1 jt=$2 bt=$3 driver=$4
 export OUTDIR="$pilot/$generation" BENCH_CASE_ROOT="$pilot/fixture-$generation"
 taskset -c "$cpulist" bash benchmark/run_cold_process.sh "$jt" "$bt" "$driver" > "$pilot/$generation.log" 2>&1
}
if [[ "$stage" != timing_profiles ]]; then
# Environment/precompilation process is isolated from tests and timed processes.
taskset -c "$cpulist" bash benchmark/run_cold_process.sh 4 1 benchmark/cold_precompile.jl > "$pilot/precompile.log" 2>&1
run_process controls-j1 1 1 test/runtests_benchmark_cold.jl
run_process controls-j4 4 1 test/runtests_benchmark_cold.jl
export STAGE=baseline
unset CONFIG_FILE
run_process smoke-j4-b1 4 1 benchmark/rotor_hover_solver_cold_smoke.jl
export CONFIG_FILE="$pilot/smoke-j4-b1/selected.toml"
if [[ "$stage" == controls_smoke ]]; then
 printf 'controls and smoke completed\n' > "$pilot/COMPLETED"
 exit 0
fi
fi
export STAGE=verify
sha256sum "$CONFIG_FILE" > "$pilot/selected.sha256"
run_process timing-j4-b1 4 1 benchmark/rotor_hover_solver_cold.jl
run_process timing-j64-b1 64 1 benchmark/rotor_hover_solver_cold.jl
run_process timing-j64-b64 64 64 benchmark/rotor_hover_solver_cold.jl
export COLD_INVESTIGATION=1
for kind in fgs krylov_ilu; do
 export CONFIGS="$kind"
 run_process "profile-$kind-j64-b1" 64 1 benchmark/rotor_hover_solver_phase2_profile.jl
done
sha256sum -c "$pilot/selected.sha256"
printf 'completed\n' > "$pilot/COMPLETED"
