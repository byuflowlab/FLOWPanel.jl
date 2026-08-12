#!/usr/bin/env bash
#
# BRAINSTORM 019 deployment: stage the max_dtZ tripwire column (monitor src +
# driver env knob) and the Campaign E screen case arms onto the cluster
# checkout, md5-verifying every copied file and spot-verifying that the rest
# of the cluster src/ matches this checkout (Phase 2c stale-cluster-src
# lesson). Nothing is submitted; this only stages files and reports.
#
# Run from the top level of the LOCAL FLOWPanel.jl checkout:
#
#   bash scripts/deploy_p019.sh
#
# Uses the `orc` SSH alias (non-interactive; raw host demands 2FA).

set -euo pipefail

REMOTE="${REMOTE:-orc}"
REMOTE_REPO="${REMOTE_REPO:-/home/rander39/projects/FLOWPanel.jl}"

COPY_FILES=(
  "src/FLOWPanel_simulate_monitors.jl"
  "examples/rotor_hover_pressure_comparison.jl"
  "examples/run_p018_screen_hpc.slurm.sh"
)

# Verified only (must already match; the screen jobs exercise these paths).
VERIFY_FILES=(
  "src/FLOWPanel_wake.jl"
  "src/FLOWPanel_simulate.jl"
  "src/FLOWPanel_liftingbody.jl"
  "examples/dji9443_trailing_edge.jl"
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
  chmod +x examples/run_p018_screen_hpc.slurm.sh && \
  echo '--- md5 ---' && \
  md5sum ${COPY_FILES[*]} ${VERIFY_FILES[*]} && \
  echo '--- branch ---' && git rev-parse --abbrev-ref HEAD && git log --oneline -1 && \
  echo '--- queue (study cap is 10 active jobs, Ryan 2026-08-04) ---' && squeue -u rander39"

echo
echo "Compare the two md5 blocks line by line; every entry must match."
echo "If a VERIFY_FILES entry differs, the cluster src is stale — reconcile"
echo "BEFORE submitting. Submission (staged under the 10-job cap):"
echo "  sbatch --job-name=fp-019-scr-<tag> examples/run_p018_screen_hpc.slurm.sh scr_p019_<tag>"
