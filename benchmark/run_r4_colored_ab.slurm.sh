#!/usr/bin/env bash
#SBATCH --job-name=p021-r4-colored-v21
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/r4-colored-v21-%j.out
#SBATCH --error=logs/slurm/r4-colored-v21-%j.err
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"

run="$COLD_DATA_ROOT/colored-v21-$SLURM_JOB_ID"
mkdir "$run"
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
getconf CLK_TCK > "$run/clock_ticks_per_second.txt"
lscpu --extended=CPU,NODE,SOCKET,CORE,ONLINE > "$run/lscpu_extended.txt"
numactl --hardware > "$run/numactl_hardware.txt" 2>&1 || true

python3 - "$run" <<'PY'
import os, pathlib, sys
out = pathlib.Path(sys.argv[1])
allowed = os.sched_getaffinity(0)
by_socket = {}
seen = set()
for cpu in sorted(allowed):
    p = pathlib.Path('/sys/devices/system/cpu') / f'cpu{cpu}' / 'topology'
    socket = int((p/'physical_package_id').read_text())
    core = int((p/'core_id').read_text())
    if (socket, core) not in seen:
        by_socket.setdefault(socket, []).append(cpu)
        seen.add((socket, core))
ordered = [cpu for socket in sorted(by_socket) for cpu in by_socket[socket]]
assert len(ordered) >= 64, 'fewer than 64 physical cores in allocation'
for n in (1,4,16,64):
    (out/f'cpu_affinity_j{n}.txt').write_text(','.join(map(str,ordered[:n]))+'\n')
PY

export RUNG=R4 CONFIGS=fgs STAGE=verify COLD_PREPARED_ONLY=1
export CONFIG_FILE="$PWD/benchmark/retained_r4_diagnostics.toml"
export COLD_AB_REPS=${COLD_AB_REPS:-10}
export COLD_AB_BATCHES=${COLD_AB_BATCHES:-4}
allcpus=$(<"$run/cpu_affinity_j64.txt")
export OUTDIR="$run/parse" BENCH_CASE_ROOT="$run/fixture-controls"
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 1 1 \
    test/runtests_r4_colored_ab_driver.jl > "$run/controls-colored-ab-driver.log" 2>&1
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 1 1 \
    benchmark/cold_parse.jl > "$run/parse.log" 2>&1
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 4 1 \
    benchmark/cold_precompile.jl > "$run/precompile.log" 2>&1
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 1 1 \
    test/runtests_benchmark_cold.jl > "$run/controls-benchmark-j1.log" 2>&1
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 4 1 \
    test/runtests_benchmark_cold.jl > "$run/controls-benchmark-j4.log" 2>&1
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 4 1 \
    test/runtests_unit_solver.jl > "$run/controls-flowpanel-solver.log" 2>&1
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 4 1 \
    test/runtests_unit_fgs_history.jl > "$run/controls-flowpanel-history.log" 2>&1
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 4 1 -e \
    'using FastMultipole, Test; fm=pkgdir(FastMultipole); include(joinpath(fm,"test","gravitational.jl")); include(joinpath(fm,"test","solve_test.jl")); include(joinpath(fm,"test","fgs_coloring_test.jl"))' \
    > "$run/controls-fastmultipole.log" 2>&1

# Separate calibration for the colored order (accumulation ordering differs;
# the retained lexicographic tolerance does not carry). Runs at j64.
export AB_MODE=calibrate
export OUTDIR="$run/calibrate/results" BENCH_CASE_ROOT="$run/fixture-calibrate"
mkdir -p "$run/calibrate"
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 64 1 \
    benchmark/fgs_r4_colored_ab.jl > "$run/calibrate/process.log" 2>&1
export COLORED_CONFIG="$run/calibrate/results/colored_selected.toml"
[ -f "$COLORED_CONFIG" ]

# Uninstrumented A/B performance trials, alternating batches, per thread arm.
export AB_MODE=trials
for jt in 1 4 16 64; do
    cpulist=$(<"$run/cpu_affinity_j$jt.txt")
    out="$run/j$jt-trials"
    mkdir "$out"
    taskset -c "$cpulist" numactl --show > "$out/numactl_show.txt" 2>&1 || true
    export OUTDIR="$out/results" BENCH_CASE_ROOT="$run/fixture-j$jt"
    taskset -c "$cpulist" bash benchmark/run_cold_process.sh "$jt" 1 \
        benchmark/fgs_r4_colored_ab.jl > "$out/process.log" 2>&1
done

# Budget attribution only (excluded from rankings): does the nearfield chain
# actually run parallel under :colored at j64?
export AB_MODE=activity
export OUTDIR="$run/j64-activity/results" BENCH_CASE_ROOT="$run/fixture-activity"
mkdir -p "$run/j64-activity"
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 64 1 \
    benchmark/fgs_r4_colored_ab.jl > "$run/j64-activity/process.log" 2>&1

printf 'completed\n' > "$run/COMPLETED"
