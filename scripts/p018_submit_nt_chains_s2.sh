#!/usr/bin/env bash
# 018 NT-ladder restart chains (_s2), staged 2026-08-24.
# Chains all 8 NT72/NT144 arms from their last complete restart sets to 30 revs
# (P018_SETTLE_REVS=22). Pins pre-2026-08-21 physics per the chain warning
# (ledger 2026-08-21 entry): vatistas + production FMM knobs (driver defaults
# body 8/0.4/20, wake 4/0.4/50 — verified: the original arms ran the pre-08-21
# launcher/driver with no FMM env, and their metadata records no knobs).
# Memory 96G (was 64G): jobs 13185010/13245449 died SIGBUS at ~67G RSS on the
# 64G request (ledger 2026-08-24 forensics). NT72 need ~35h -> 48h; NT144
# cannot finish 30 revs in 72h (~163 s/step) — s2 reaches ~rev 22, an s3
# chain from s2's last snapshot will be needed.
# RESTART_STEP values = min over {body .vtu, wake .1/.2 .vts, particles .vtp}
# max indices, verified aligned on 2026-08-24.
# Run ON the cluster from the repo top level:
#   bash scripts/p018_submit_nt_chains_s2.sh
set -euo pipefail
cd /home/rander39/projects/FLOWPanel.jl
export PATH=/apps/slurm/latest/bin:$PATH
mkdir -p logs/slurm

PIN="FLOWPANEL_FILAMENT_REG=vatistas,FMM_BODY_EXPANSION_ORDER=8,FMM_BODY_ACCEPTANCE=0.4,FMM_BODY_LEAF_SIZE=20,FMM_WAKE_EXPANSION_ORDER=4,FMM_WAKE_ACCEPTANCE=0.4,FMM_WAKE_LEAF_SIZE=50"

submit () { # $1=run $2=restart_step $3=walltime
  sbatch --job-name=fp-018-${1#p018_}_s2 --time=$3 --mem=96G \
    --export=ALL,P018_RUN_NAME=${1}_s2,P018_SETTLE_REVS=22,RESTART_STEP=$2,RESTART_NAME=$1,RESTART_PATH=data/$1,$PIN \
    examples/run_dji9443_hover_ct_hpc.slurm.sh "$1"
}

submit p018_csarc_nt72_l2p4     1453 48:00:00
submit p018_csarc_nt72_l3p4     1459 48:00:00
submit p018_csarc_nt72_l4p8     1518 48:00:00
submit p018_csarc_n0_nt72_l4p8  1449 48:00:00
submit p018_csarc_nt144_l2p4    1648 72:00:00
submit p018_csarc_nt144_l3p4    1682 72:00:00
submit p018_csarc_nt144_l4p8    1663 72:00:00
submit p018_csarc_n0_nt144_l4p8 1648 72:00:00

squeue -u rander39 -o "%.10i %.28j %.2t %.10M"
echo "REMINDER: verify each banner in logs/slurm/slurm-fp-018-*_s2-*.out within"
echo "minutes of start: expect filament_reg:vatistas fmm_body:8/0.4/20"
echo "fmm_wake:4/0.4/50, correct NT/rlxf/pps/lambda, settle:22, and a"
echo "'restart' line naming the parent run and step."
