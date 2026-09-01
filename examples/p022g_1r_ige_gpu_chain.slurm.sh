#!/usr/bin/env bash
#SBATCH --job-name=fp-022g-1r-ige-chain
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# 052b single-allocation gate chain for the 1-rotor IGE GPU case
# (2026-09-01, Ryan-approved plan): smoke -> probe -> projection gate ->
# 414-step accept, all inside one H200 allocation so the queue wait is paid
# once. Each stage is a fresh carrier invocation
# (run_rotor_multi_ground_effect_gpu.slurm.sh) with its own RUN_NAME; any
# stage failure aborts the chain. The projection gate refuses the accept
# stage unless median late-probe step time projects the 414-step case under
# the 7200 s budget with 10% reserve and 1.5x growth margin, and the probe
# left >=20% device memory free.
set -uo pipefail

CARRIER=examples/run_rotor_multi_ground_effect_gpu.slurm.sh
# Generalized 2026-09-01: case tag as $1 (default keeps the original 1r IGE
# behavior). Submit other cases with e.g.
#   sbatch -J fp-022g-2r-ige-chain -p ... p022g_1r_ige_gpu_chain.slurm.sh p022g_2r_ige
CASE="${1:-p022g_1r_ige}"
export P022G_REQUIRED_GPU_MODEL="${P022G_REQUIRED_GPU_MODEL:-H200}"
export P022G_EXISTING_RESULT=preserve
export P022G_CONFIRM_ACCEPTANCE=YES

mkdir -p logs/chain
JOB=${SLURM_JOB_ID:-manual}

run_stage() { # mode run_name
  local mode=$1 name=$2
  # NOTE: must be a separate `local` line — bash expands all args of one
  # `local` before running it, so ${mode} on the same line is unset (fatal
  # under set -u; killed job 13549590 in 2 s).
  local log=logs/chain/${CASE}_${mode}_${JOB}.log
  echo "=== chain stage $mode ($name) start $(date -Is)"
  P022G_MODE=$mode P022_RUN_NAME=$name bash "$CARRIER" "$CASE" 2>&1 | tee "$log"
  local rc=${PIPESTATUS[0]}
  echo "=== chain stage $mode exit=$rc $(date -Is)"
  return $rc
}

run_stage smoke ${CASE}_smoke || { echo "CHAIN ABORT: smoke failed"; exit 1; }
run_stage probe ${CASE}_probe || { echo "CHAIN ABORT: probe failed"; exit 1; }

plog=logs/chain/${CASE}_probe_${JOB}.log
med=$(grep -o 'step_timer total_step step=[0-9]* [0-9.eE+-]* s' "$plog" \
  | awk '{print $4}' | tail -30 | sort -g \
  | awk '{a[NR]=$1} END{if (NR==0) print "NA"; else print a[int((NR+1)/2)]}')
[[ "$med" != NA ]] || { echo "CHAIN ABORT: no step_timer total_step lines in probe log"; exit 1; }
proj=$(awk -v m="$med" 'BEGIN{printf "%.0f", 414*m*1.5}')
echo "projection: median_late_probe_step=${med}s projected_accept=${proj}s (limit 6480s)"
(( proj <= 6480 )) || { echo "CHAIN ABORT: projected accept ${proj}s exceeds 6480s"; exit 1; }

meta=data/${CASE}_probe/${CASE}_probe_case_metadata.toml
reserve=$(awk -F' = ' '/device_memory_reserve_fraction/{print $2}' "$meta" 2>/dev/null | tail -1)
echo "probe device memory reserve fraction: ${reserve:-unknown}"
awk -v r="${reserve:-0}" 'BEGIN{exit !(r >= 0.20)}' \
  || { echo "CHAIN ABORT: probe device reserve ${reserve:-unknown} below 0.20"; exit 1; }

if [[ "${P022G_CHAIN_STOP_AFTER:-}" == probe ]]; then
  echo "CHAIN STOP: P022G_CHAIN_STOP_AFTER=probe (gates evaluated, accept skipped)"
  exit 0
fi

run_stage accept "$CASE" || { echo "CHAIN FAIL: accept failed"; exit 1; }
echo "CHAIN COMPLETE $(date -Is)"
