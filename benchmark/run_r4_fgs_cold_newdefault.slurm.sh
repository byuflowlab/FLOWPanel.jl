#!/usr/bin/env bash
#SBATCH --job-name=p021-fgs-cold-newdefault
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=normal
#SBATCH --time=08:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/slurm/r4-fgs-cold-newdefault-%A_%a.out
#SBATCH --error=logs/slurm/r4-fgs-cold-newdefault-%A_%a.err
# 021 Job 2 (fgs_warmstart_r4_reset_prompt_20260924.md): cold R4 FGS under the
# NEW production default (sweep_order=:dagteam, dagteam_idle=:backoff) across
# j in {1,8,16,32,64}, for direct comparison against the certified iLU-GMRES
# columns already harvested (2026-09-22 R4 table). The krylov_ilu /
# krylov_ilu_nfcache columns are NOT rerun.
#
# Same contract as the existing rows: per-j staircase-calibrated tolerance,
# certified BC evaluator on every accepted solve, t_solve_min on cold
# (zero-initial-guess) isolated solves, champion knobs P8/MAC0.4/leaf100,
# f32full (its cold-R4/zen3 certification is exactly this environment),
# champion placement, BLAS=1. The per-j calibrated configs are REUSED from the
# 13777133 thread-scaling ladder with dagteam_idle=backoff injected — valid
# because the backoff idle policy changes scheduling only, never arithmetic
# (FastMultipole solve_dagteam.jl contract; the dagedge campaign reused the
# same calibrations the same way and its cross-executor delta was 2.24e-7
# against a 1e-5 tripwire). Every trial is still independently gated by the
# certified evaluator, so a transfer failure would be caught, not absorbed.
#
# New run dirs alongside thread-scaling-j<J>-13777133 (never write into the
# old ones). Judge by outputs (STATUS_*/status.toml), never sacct.
#
# Required env: COLD_PROJECT, CAMPAIGN_PINS, COLD_DATA_ROOT
# Optional: CALIB_ROOT (default COLD_DATA_ROOT), CALIB_JOB (default 13777133)
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"

LADDER=(1 8 16 32 64)
J="${LADDER[$SLURM_ARRAY_TASK_ID]}"

NODE_CPUS=$(getconf _NPROCESSORS_ONLN)
ALLOC_CPUS=${SLURM_CPUS_ON_NODE:-0}
if [ "$NODE_CPUS" != "128" ] || [ "$ALLOC_CPUS" != "128" ]; then
  echo "ERROR: need an exclusive 128-core zen3 node (got node=$NODE_CPUS" >&2
  echo "       alloc=$ALLOC_CPUS); timings would not be comparable." >&2
  exit 1
fi

export JULIA_PKG_PRECOMPILE_AUTO=0
if command -v flock >/dev/null 2>&1; then
  ( if flock -E 99 -w 3600 9; then
      julia --project="$COLD_PROJECT" --startup-file=no \
        -e 'using Pkg; Pkg.precompile()' \
        || echo "WARNING: Pkg.precompile() failed; compiling in memory"
    else
      echo "WARNING: precompile lock unavailable; compiling in memory"
    fi
  ) 9>"$HOME/.julia/flowpanel-021-precompile.lock"
else
  echo "WARNING: flock(1) not found; skipping shared precompile"
fi

CALIB_ROOT="${CALIB_ROOT:-$COLD_DATA_ROOT}"
CALIB_JOB="${CALIB_JOB:-13777133}"
CALIB="$CALIB_ROOT/thread-scaling-j$J-$CALIB_JOB/fgs-calibrate/results/dagteam_selected.toml"
[ -f "$CALIB" ] || { echo "ERROR: missing calibrated config $CALIB" >&2; exit 1; }

run="$COLD_DATA_ROOT/fgs-cold-newdefault-j$J-${RESUME_FROM_JOB_ID:-$SLURM_ARRAY_JOB_ID}"
if [ -n "${RESUME_FROM_JOB_ID:-}" ] && [ -d "$run" ]; then
    rm -f "$run/COMPLETED"
else
    mkdir "$run"
fi
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
numactl --hardware > "$run/numactl_hardware.txt" 2>&1 || true
echo "julia_threads=$J calib=$CALIB" > "$run/point.txt"

ILV="--interleave=0-3 --cpunodebind=0-3"
numactl $ILV numactl --show > "$run/numactl_show.txt" 2>&1 || true

export HARDWARE_TAG="${HARDWARE_TAG:-orc-m12-zen3-socket0-ilv0-3-blas1}"
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-linegauss}"
export JULIA_DEBUG="${JULIA_DEBUG:-loading}"

# Inject the new default idle policy into the calibrated config (a pure
# scheduling knob; the calibrated tolerance and knobs are untouched). Assert
# the champion knob set while at it so a stale/foreign config cannot slip in.
CONFIG_INJECTED="$run/dagteam_backoff_config.toml"
julia --startup-file=no -e '
    using TOML
    doc = TOML.parsefile(ARGS[1])
    cs = haskey(doc, "configs") ? doc["configs"] : [doc]
    length(cs) == 1 || error("expected exactly one calibrated config")
    c = cs[1]
    (c["kind"] == "fgs" && c["sweep_order"] == "dagteam" && c["rung"] == "R4" &&
     c["P"] == 8 && c["MAC"] == 0.4 && c["leaf"] == 100) ||
        error("calibrated config is not the R4 dagteam champion: $(c)")
    c["tolerance"] > 0 || error("calibrated config carries no tolerance")
    c["dagteam_idle"] = "backoff"
    open(io -> TOML.print(io, Dict("configs" => [c]); sorted=true), ARGS[2], "w")
' "$CALIB" "$CONFIG_INJECTED"

status() { printf '%s\n' "$2" > "$run/STATUS_$1"; }

export RUNG=R4 CONFIGS=fgs STAGE=verify COLD_PREPARED_ONLY=1
export CONFIG_FILE="$CONFIG_INJECTED"

export OUTDIR="$run/parse" BENCH_CASE_ROOT="$run/fixture-controls"
bash benchmark/run_cold_process.sh 1 1 \
    benchmark/cold_parse.jl > "$run/parse.log" 2>&1
bash benchmark/run_cold_process.sh 4 1 \
    benchmark/cold_precompile.jl > "$run/precompile.log" 2>&1

export OUTDIR="$run/results" BENCH_CASE_ROOT="$run/fixture"
if numactl $ILV bash benchmark/run_cold_process.sh "$J" 1 \
        benchmark/rotor_hover_solver_cold.jl > "$run/process.log" 2>&1 \
    && [ -f "$run/results/selected.toml" ]; then
  status verify ok
else
  status verify FAILED
  echo "WARNING: j$J verify FAILED (see $run/process.log)"
fi

printf 'completed j=%s\n' "$J" > "$run/COMPLETED"
