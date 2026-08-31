#!/usr/bin/env bash
# BRAINSTORM/022 multi-rotor: TWO rotors IN ground effect. Run ON the cluster.
# h/R is env-overridable: GROUND_H_R=1.0 scripts/p022m_submit_2r_ige.sh
# (default 1.5). TRUNCATION_DEPTH_R tracks h/R + 3 (022 convention: fixed 3R
# below-ground wake allowance). Optional: P022M_MESH (default 45_185_ct4 fine —
# coarse 56_57 blew up near ground in single-rotor 022), RHPC_KEEP_PREV.
set -euo pipefail
cd /home/rander39/projects/FLOWPanel.jl
export PATH=/apps/slurm/latest/bin:$PATH   # sbatch absent in non-interactive ssh
mkdir -p logs/slurm

HR="${GROUND_H_R:-1.5}"
DEPTH=$(awk -v h="$HR" 'BEGIN{printf "%g", h + 3}')
TAG="hr$(printf '%s' "$HR" | tr -d '.')"

sbatch --job-name="fp-022m-2r-ige-$TAG" --time=72:00:00 --mem=96G \
  --export=ALL,IGE_H_R="$HR",IGE_TRUNC_DEPTH="$DEPTH",P022_RUN_NAME="p022m_2r_ige_$TAG" \
  examples/run_rotor_multi_ground_effect_hpc.slurm.sh p022m_2r_ige

squeue -u "$USER" -o "%.10i %.24j %.8T %.10M %.9l %R" | tail -5
echo "Verify the banner (nrotors:2, ground:true h/R:$HR depth:${DEPTH}R, reg:gaussian) once the job starts."
