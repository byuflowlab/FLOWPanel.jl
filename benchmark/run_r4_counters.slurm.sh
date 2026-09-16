#!/usr/bin/env bash
# Follow BYU_ORC_AGENTS.md; this job measures only its own child processes.
#SBATCH --job-name=p021-r4-counters-v17
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=06:00:00
#SBATCH --output=logs/slurm/r4-counters-v17-%j.out
#SBATCH --error=logs/slurm/r4-counters-v17-%j.err
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"
run="$COLD_DATA_ROOT/counters-v17-$SLURM_JOB_ID"
mkdir "$run"
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
lscpu --extended=CPU,NODE,SOCKET,CORE,ONLINE > "$run/lscpu_extended.txt"
getconf CLK_TCK > "$run/clock_ticks_per_second.txt"
numactl --hardware > "$run/numactl_hardware.txt" 2>&1 || true
perf list > "$run/perf_events.txt" 2>&1
{
    date -u
    cat /proc/sys/kernel/perf_event_paranoid
    for tool in perf likwid-perfctr AMDuProfPcm; do command -v "$tool" || true; done
    # Inventory only: uncore events cannot be attributed to this process.
    ls /sys/bus/event_source/devices
} > "$run/counter_capabilities.txt" 2>&1
mkfifo "$run/smoke-control.fifo" "$run/smoke-ack.fifo"
perf stat -D -1 --control="fifo:$run/smoke-control.fifo,$run/smoke-ack.fifo" \
    -x, -o "$run/perf-smoke.csv" \
    -e '{cycles:u,instructions:u,cache-references:u,cache-misses:u}' -e task-clock \
    -- env JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
        julia --startup-file=no --project="$COLD_PROJECT" test/r4_perf_control_smoke.jl "$run/smoke-control.fifo" "$run/smoke-ack.fifo" \
    > "$run/perf-smoke.log" 2>&1
rm "$run/smoke-control.fifo" "$run/smoke-ack.fifo"
python3 - "$run" <<'PY'
import os, pathlib, sys
out = pathlib.Path(sys.argv[1]); cores = {}
for cpu in sorted(os.sched_getaffinity(0)):
    p = pathlib.Path('/sys/devices/system/cpu') / f'cpu{cpu}' / 'topology'
    key = (int((p/'physical_package_id').read_text()), int((p/'core_id').read_text()))
    cores.setdefault(key, cpu)
ordered = [cores[k] for k in sorted(cores)]
assert len(ordered) >= 64
for n in (4, 64):
    (out/f'cpu_affinity_j{n}.txt').write_text(','.join(map(str, ordered[:n]))+'\n')
PY

export RUNG=R4 CONFIGS=fgs STAGE=verify COLD_PREPARED_ONLY=1
export CONFIG_FILE="$PWD/benchmark/retained_r4_diagnostics.toml"
export OUTDIR="$run/parse" BENCH_CASE_ROOT="$run/fixture-controls"
bash benchmark/run_cold_process.sh 1 1 test/runtests_r4_counters_driver.jl \
    > "$run/controls-counter-driver.log" 2>&1
bash benchmark/run_cold_process.sh 4 1 -e \
    'using FastMultipole, Test; fm=pkgdir(FastMultipole); include(joinpath(fm,"test","gravitational.jl")); include(joinpath(fm,"test","solve_test.jl")); include(joinpath(fm,"test","fgs_coloring_test.jl"))' \
    > "$run/controls-fastmultipole.log" 2>&1
bash benchmark/run_cold_process.sh 4 1 test/runtests_unit_solver.jl \
    > "$run/controls-flowpanel-solver.log" 2>&1
bash benchmark/run_cold_process.sh 4 1 test/runtests_unit_fgs_history.jl \
    > "$run/controls-flowpanel-history.log" 2>&1

for jt in 4 64; do
    out="$run/j$jt-b1"
    mkdir "$out"
    export COLD_PERF_CONTROL="$out/perf-control.fifo" COLD_PERF_ACK="$out/perf-ack.fifo"
    mkfifo "$COLD_PERF_CONTROL" "$COLD_PERF_ACK"
    export OUTDIR="$out/results" BENCH_CASE_ROOT="$run/fixture-j$jt"
    cpulist=$(<"$run/cpu_affinity_j$jt.txt")
    taskset -c "$cpulist" numactl --show > "$out/numactl_show.txt" 2>&1 || true
    # Events begin disabled. Julia enables only around the warmed prepared solve.
    taskset -c "$cpulist" perf stat -D -1 \
        --control="fifo:$COLD_PERF_CONTROL,$COLD_PERF_ACK" \
        -x, -o "$out/perf-stat.csv" \
        -e '{cycles:u,instructions:u,cache-references:u,cache-misses:u}' -e task-clock \
        -- bash benchmark/run_cold_process.sh "$jt" 1 benchmark/fgs_r4_counters.jl \
        > "$out/process.log" 2>&1
    rm "$COLD_PERF_CONTROL" "$COLD_PERF_ACK"
done
printf 'completed\n' > "$run/COMPLETED"
