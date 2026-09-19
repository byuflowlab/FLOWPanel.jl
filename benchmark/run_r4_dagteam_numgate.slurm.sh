#!/usr/bin/env bash
#SBATCH --job-name=p021-r4-dagteam-numgate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --constraint=zen3
#SBATCH --exclusive
#SBATCH --qos=test
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm/r4-dagteam-numgate-%j.out
#SBATCH --error=logs/slurm/r4-dagteam-numgate-%j.err
# 021 v23 numerical gate (spec gate 3): staircase-calibrate the colored and
# dagteam twins of the retained R4 lex config at the champion thread count and
# placement. Every accepted solve passes the independent evaluator at BC
# rel-L2 <= 1e-6; a DAGTEAM_PRECISION rung that cannot certify fails this job,
# which is the early accuracy verdict ahead of the full A/B campaign. Lean by
# design (parse + precompile + calibrate only) to fit qos=test; the A/B job
# repeats calibration self-contained, so this job's failure or timeout never
# blocks it.
set -euo pipefail
set +u
source /etc/profile
set -u
module load cuda/12.8.1-zkkfiog julia/1.11.7-6bmogfl
: "${COLD_PROJECT:?}" "${CAMPAIGN_PINS:?}" "${COLD_DATA_ROOT:?}"

run="$COLD_DATA_ROOT/dagteam-numgate-$SLURM_JOB_ID"
mkdir "$run"
cp "$CAMPAIGN_PINS" "$run/campaign_pins.toml"
cp "$COLD_PROJECT/Manifest.toml" "$run/Manifest.toml"
module list > "$run/modules.txt" 2>&1
lscpu > "$run/lscpu.txt"
numactl --hardware > "$run/numactl_hardware.txt" 2>&1 || true

ILV="--interleave=0-3 --cpunodebind=0-3"

export RUNG=R4 CONFIGS=fgs STAGE=verify COLD_PREPARED_ONLY=1
export CONFIG_FILE="$PWD/benchmark/retained_r4_diagnostics.toml"
export DAGTEAM_PRECISION="${DAGTEAM_PRECISION:-f32full}"

export OUTDIR="$run/parse" BENCH_CASE_ROOT="$run/fixture-controls"
bash benchmark/run_cold_process.sh 1 1 \
    benchmark/cold_parse.jl > "$run/parse.log" 2>&1
bash benchmark/run_cold_process.sh 4 1 \
    benchmark/cold_precompile.jl > "$run/precompile.log" 2>&1

export AB_MODE=calibrate
export OUTDIR="$run/calibrate/results" BENCH_CASE_ROOT="$run/fixture-calibrate"
mkdir -p "$run/calibrate"
numactl $ILV bash benchmark/run_cold_process.sh 16 1 \
    benchmark/fgs_r4_dagteam_ab.jl > "$run/calibrate/process.log" 2>&1
[ -f "$run/calibrate/results/dagteam_selected.toml" ]
[ -f "$run/calibrate/results/colored_selected.toml" ]

printf 'completed\n' > "$run/COMPLETED"
