#!/usr/bin/env bash
#SBATCH --job-name=p021-r4-dagteam-v23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=10:00:00
#SBATCH --output=logs/slurm/r4-dagteam-v23-%j.out
#SBATCH --error=logs/slurm/r4-dagteam-v23-%j.err
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"

run="$COLD_DATA_ROOT/dagteam-v23-$SLURM_JOB_ID"
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

# Champion placement from gate 2d (job 13773581): threads on socket 0
# (NUMA nodes 0-3), memory interleaved across its four controllers.
ILV="--interleave=0-3 --cpunodebind=0-3"

export RUNG=R4 CONFIGS=fgs STAGE=verify COLD_PREPARED_ONLY=1
export CONFIG_FILE="$PWD/benchmark/retained_r4_diagnostics.toml"
export DAGTEAM_PRECISION="${DAGTEAM_PRECISION:-f32full}"
export COLD_AB_REPS=${COLD_AB_REPS:-10}
export COLD_AB_BATCHES=${COLD_AB_BATCHES:-4}
allcpus=$(<"$run/cpu_affinity_j64.txt")
export OUTDIR="$run/parse" BENCH_CASE_ROOT="$run/fixture-controls"
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 1 1 \
    test/runtests_r4_dagteam_ab_driver.jl > "$run/controls-dagteam-ab-driver.log" 2>&1
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
# Gate-1 dagteam correctness harness runs standalone so its self-include
# preamble (compat overloads) activates.
fmdir=$(awk -F'"' '/\[packages.FastMultipole\]/{f=1} f&&/^path/{print $2; exit}' "$CAMPAIGN_PINS")
[ -d "$fmdir/test" ]
taskset -c "$allcpus" bash benchmark/run_cold_process.sh 4 1 \
    "$fmdir/test/fgs_dagteam_gate1_test.jl" > "$run/controls-dagteam-gate1.log" 2>&1

# Separate calibration for BOTH twins at the champion thread count and
# placement (coloring changes accumulation ordering; dagteam changes it
# further and, at reduced precision, perturbs the operator — neither
# tolerance carries). The staircase + independent evaluator (BC rel-L2 <=
# 1e-6) is the numerical gate for DAGTEAM_PRECISION. Dagteam iterates are
# thread-count invariant (bitwise per run), so one calibration carries.
export AB_MODE=calibrate
export OUTDIR="$run/calibrate/results" BENCH_CASE_ROOT="$run/fixture-calibrate"
mkdir -p "$run/calibrate"
numactl $ILV bash benchmark/run_cold_process.sh 16 1 \
    benchmark/fgs_r4_dagteam_ab.jl > "$run/calibrate/process.log" 2>&1
export COLORED_CONFIG="$run/calibrate/results/colored_selected.toml"
export DAGTEAM_CONFIG="$run/calibrate/results/dagteam_selected.toml"
[ -f "$COLORED_CONFIG" ]
[ -f "$DAGTEAM_CONFIG" ]

# Uninstrumented A/B performance trials, alternating batches.
#   j16-interleave: the decision arm (champion count + champion placement).
#   j16-native:     control anchor to the v21 accepted environment
#                   (taskset physical-core list, default placement).
#   j64-interleave: robustness arm at the full socket.
export AB_MODE=trials
for arm in j16-interleave j16-native j64-interleave; do
    out="$run/$arm-trials"
    mkdir "$out"
    export OUTDIR="$out/results" BENCH_CASE_ROOT="$run/fixture-$arm"
    case "$arm" in
    j16-interleave)
        numactl $ILV --show > "$out/numactl_show.txt" 2>&1 || true
        numactl $ILV bash benchmark/run_cold_process.sh 16 1 \
            benchmark/fgs_r4_dagteam_ab.jl > "$out/process.log" 2>&1 ;;
    j16-native)
        cpulist=$(<"$run/cpu_affinity_j16.txt")
        taskset -c "$cpulist" numactl --show > "$out/numactl_show.txt" 2>&1 || true
        taskset -c "$cpulist" bash benchmark/run_cold_process.sh 16 1 \
            benchmark/fgs_r4_dagteam_ab.jl > "$out/process.log" 2>&1 ;;
    j64-interleave)
        numactl $ILV bash benchmark/run_cold_process.sh 64 1 \
            benchmark/fgs_r4_dagteam_ab.jl > "$out/process.log" 2>&1 ;;
    esac
done

# Budget attribution only (excluded from rankings): where does the dagteam
# solve spend its stages at the champion configuration?
export AB_MODE=activity
export OUTDIR="$run/j16-activity/results" BENCH_CASE_ROOT="$run/fixture-activity"
mkdir -p "$run/j16-activity"
numactl $ILV bash benchmark/run_cold_process.sh 16 1 \
    benchmark/fgs_r4_dagteam_ab.jl > "$run/j16-activity/process.log" 2>&1

printf 'completed\n' > "$run/COMPLETED"
