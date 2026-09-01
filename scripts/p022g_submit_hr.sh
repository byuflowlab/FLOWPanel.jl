#!/usr/bin/env bash
# Submit a 022 Phase-5 h/R sweep arm on GPU (carrier:
# examples/run_rotor_multi_ground_effect_gpu.slurm.sh, case p022g_hr<NN>).
#
#   P022G_GRES=gpu:h200:1 [P022G_PARTITION=m13h] [P022G_QOS=gpu] \
#   [P022G_MODE=smoke|probe|production] [P022G_TIME=HH:MM:SS] \
#   scripts/p022g_submit_hr.sh <05|10|15|20>
#
# Submit from /home/rander39/projects_unified/FLOWPanel.jl.
set -euo pipefail
HR="${1:-}"
case "$HR" in 05|10|15|20) ;; *) echo "ERROR: usage: p022g_submit_hr.sh <05|10|15|20>" >&2; exit 2 ;; esac
: "${P022G_GRES:?set P022G_GRES, e.g. gpu:h200:1}"
mkdir -p logs/slurm
sbatch ${P022G_PARTITION:+--partition="$P022G_PARTITION"} \
       ${P022G_QOS:+--qos="$P022G_QOS"} \
       ${P022G_TIME:+--time="$P022G_TIME"} \
       --gres="$P022G_GRES" --job-name="fp-022g-hr$HR" \
       examples/run_rotor_multi_ground_effect_gpu.slurm.sh "p022g_hr$HR"
