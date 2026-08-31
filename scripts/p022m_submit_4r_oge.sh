#!/usr/bin/env bash
# BRAINSTORM/022 multi-rotor: FOUR rotors, OUT of ground effect. Run ON the
# cluster. Optional env: P022M_MESH (default 56_57 for OGE), RHPC_KEEP_PREV.
set -euo pipefail
cd /home/rander39/projects/FLOWPanel.jl
export PATH=/apps/slurm/latest/bin:$PATH   # sbatch absent in non-interactive ssh
mkdir -p logs/slurm

sbatch --job-name=fp-022m-4r-oge --time=72:00:00 --mem=128G \
  --export=ALL,P022_RUN_NAME=p022m_4r_oge \
  examples/run_rotor_multi_ground_effect_hpc.slurm.sh p022m_4r_oge

squeue -u "$USER" -o "%.10i %.24j %.8T %.10M %.9l %R" | tail -5
echo "Verify the banner (nrotors:4, ground:false, reg:gaussian) once the job starts."
