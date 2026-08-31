#!/usr/bin/env bash
set -euo pipefail
: "${P022G_GRES:?set P022G_GRES, e.g. gpu:gh200:1}"
sbatch ${P022G_PARTITION:+--partition="$P022G_PARTITION"} ${P022G_QOS:+--qos="$P022G_QOS"} --gres="$P022G_GRES" --job-name=fp-022g-1r-ige examples/run_rotor_multi_ground_effect_gpu.slurm.sh p022g_1r_ige
