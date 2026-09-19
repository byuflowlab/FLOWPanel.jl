#!/usr/bin/env bash
#SBATCH --job-name=p021-numa-dgemv
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=00:30:00
#SBATCH --output=logs/slurm/p021-numa-dgemv-%j.out
#SBATCH --error=logs/slurm/p021-numa-dgemv-%j.err
# NUMA-placement dgemv microbenchmark (021, 2026-09-18). Diagnostic, not a
# campaign: pure Julia + BLAS, touches no repo code. Pinning copied from
# benchmark/run_r4_chunked_ab.slurm.sh (first 64 physical cores socket-major
# -> cpubind NUMA nodes 0-3; default first-touch mempolicy except arm c).
set -euo pipefail
set +u
source /etc/profile
set -u
module load julia/1.11.7-6bmogfl
export JULIA_PKG_PRECOMPILE_AUTO=0

run="$HOME/projects/FLOWPanel.jl/data/p021-cold-20260910/numa-bench-$SLURM_JOB_ID"
mkdir -p "$run"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
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
for n in (1, 64):
    (out/f'cpu_affinity_j{n}.txt').write_text(','.join(map(str,ordered[:n]))+'\n')
PY

cpu1=$(<"$run/cpu_affinity_j1.txt")
cpus64=$(<"$run/cpu_affinity_j64.txt")
bench="$HOME/projects/FLOWPanel.jl/benchmark/numa_dgemv_bench.jl"
export OUTDIR="$run"

# arm a: 1 thread, single-thread first-touch (lex analogue)
ARM=a taskset -c "$cpu1" julia --startup-file=no -O2 -t 1 "$bench" \
    > "$run/arm_a.log" 2>&1

# arm b: 64 threads static, single-thread first-touch (chunked v22 analogue)
ARM=b taskset -c "$cpus64" julia --startup-file=no -O2 -t 64 "$bench" \
    > "$run/arm_b.log" 2>&1

# arm c: as b, pages interleaved across NUMA nodes 0-3
taskset -c "$cpus64" numactl --show > "$run/numactl_show_arm_c.txt" 2>&1 || true
ARM=c taskset -c "$cpus64" numactl --interleave=0-3 \
    julia --startup-file=no -O2 -t 64 "$bench" > "$run/arm_c.log" 2>&1

# arm d: 64 threads, parallel chunk-affine first-touch
ARM=d taskset -c "$cpus64" julia --startup-file=no -O2 -t 64 "$bench" \
    > "$run/arm_d.log" 2>&1

# summary
{ head -n1 "$run/result_arma.csv"; for arm in a b c d; do tail -n1 "$run/result_arm$arm.csv"; done; } \
    > "$run/summary.csv"
printf 'completed\n' > "$run/COMPLETED"
