#!/usr/bin/env bash
#
# BRAINSTORM/022 deployment: copy the ground-effect driver + launcher to the
# cluster checkout, md5-verify them, and verify that every src/ file the
# driver depends on MATCHES this checkout (018 lesson: stale/drifted cluster
# src silently produces wrong results). Nothing is submitted.
#
# Run from the top level of the LOCAL FLOWPanel.jl checkout:
#
#   bash scripts/deploy_022.sh
#
# Uses the `orc` SSH alias (needs a live ControlMaster socket; the raw host
# demands 2FA under BatchMode).

set -euo pipefail

REMOTE="${REMOTE:-orc}"
REMOTE_REPO="${REMOTE_REPO:-/home/rander39/projects/FLOWPanel.jl}"

COPY_FILES=(
  "examples/rotor_hover_ground_effect.jl"
  "examples/run_rotor_ground_effect_hpc.slurm.sh"
)

# Verified only (never copied): src the driver exercises, TE detector, and the
# two meshes the four 022 cases use (45_185_ct4 was deployed by 018; 56_57 is
# a legacy tracked mesh).
VERIFY_FILES=(
  "examples/dji9443_trailing_edge.jl"
  "examples/data/dji9443_20260725_45_185_capped_captess4.msh"
  "examples/data/dji9443_56_57.msh"
  "src/FLOWPanel_solver.jl"
  "src/FLOWPanel_nonliftingbody.jl"
  "src/FLOWPanel_formulation.jl"
  "src/FLOWPanel_simulate.jl"
  "src/FLOWPanel_liftingbody.jl"
  "src/FLOWPanel_wake.jl"
  "src/FLOWPanel_instrumentation.jl"
  "Manifest.toml"
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
  chmod +x examples/run_rotor_ground_effect_hpc.slurm.sh && \
  echo '--- md5 ---' && \
  md5sum ${COPY_FILES[*]} ${VERIFY_FILES[*]} && \
  echo '--- branch ---' && git rev-parse --abbrev-ref HEAD && git log --oneline -1 && \
  echo '--- queue ---' && squeue -u rander39"

echo
echo "Compare the two md5 blocks line by line. Every entry must match."
echo "A VERIFY_FILES mismatch means cluster src has drifted — RECONCILE WITH"
echo "THE LIVE 018 STATE (do not blindly overwrite) before submitting 022."
