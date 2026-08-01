#!/usr/bin/env bash
#
# Phase 2e deployment: copy the three assets the DJI 9443 unsteady hover CT
# study needs onto the cluster checkout, verify every one by md5, and verify
# that the cluster's src/ matches this checkout (the Phase 2c stale-cluster-src
# lesson: a stale cluster source silently produced a monitor bug).
#
# Run from the top level of the LOCAL FLOWPanel.jl checkout:
#
#   bash scripts/deploy_phase02e.sh
#
# Use the `orc` SSH alias, which authenticates non-interactively. The raw host
# rander39@ssh.rc.byu.edu demands keyboard-interactive 2FA and fails under
# BatchMode. `squeue`/`sbatch` are not on the non-interactive PATH, hence the
# `bash -lc` wrapper on the remote command below.
#
# Nothing is submitted; this only stages files and reports verification.

set -euo pipefail

REMOTE="${REMOTE:-orc}"
REMOTE_REPO="${REMOTE_REPO:-/home/rander39/projects/FLOWPanel.jl}"

# Files that must be copied: the untracked production mesh, the modified driver,
# and the new launcher.
COPY_FILES=(
  "examples/data/dji9443_20260725_45_185_capped_captess4.msh"
  "examples/rotor_hover_pressure_comparison.jl"
  "examples/run_dji9443_hover_ct_hpc.slurm.sh"
)

# Files that must already MATCH on the cluster (not copied; verified only).
VERIFY_FILES=(
  "examples/dji9443_trailing_edge.jl"
  "src/FLOWPanel_formulation.jl"
  "src/FLOWPanel_simulate.jl"
  "src/FLOWPanel_liftingbody.jl"
)

md5_local() { md5 -q "$1" 2>/dev/null || md5sum "$1" | awk '{print $1}'; }

echo "== local md5 =="
for f in "${COPY_FILES[@]}" "${VERIFY_FILES[@]}"; do
  printf '  %s  %s\n' "$(md5_local "$f")" "$f"
done

echo
echo "== copying to $REMOTE:$REMOTE_REPO =="
for f in "${COPY_FILES[@]}"; do
  echo "  -> $f"
  scp -q "$f" "$REMOTE:$REMOTE_REPO/$f"
done

echo
echo "== remote md5 + prerequisites =="
ssh "$REMOTE" bash -lc "cd '$REMOTE_REPO' && \
  mkdir -p logs/slurm && \
  chmod +x examples/run_dji9443_hover_ct_hpc.slurm.sh && \
  echo '--- md5 ---' && \
  md5sum ${COPY_FILES[*]} ${VERIFY_FILES[*]} && \
  echo '--- branch ---' && git rev-parse --abbrev-ref HEAD && git log --oneline -1 && \
  echo '--- logs/slurm ---' && ls -d logs/slurm && \
  echo '--- queue (study cap is 3 active jobs) ---' && squeue -u rander39"

echo
echo "Compare the two md5 blocks above line by line. Every entry must match."
echo "If a VERIFY_FILES entry differs, the cluster src is stale — reconcile it"
echo "(git pull / commit + push) BEFORE submitting anything."
