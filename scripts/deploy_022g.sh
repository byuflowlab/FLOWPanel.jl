#!/usr/bin/env bash
set -euo pipefail

# Target: the unified checkout (HPC.md ruling 2026-08-31; the 022g silo was
# never created and is superseded).
REMOTE="${P022G_REMOTE:-orc}"
REMOTE_DIR="${P022G_REMOTE_DIR:-/home/rander39/projects_unified/FLOWPanel.jl}"
[[ "$REMOTE_DIR" == /home/rander39/projects_unified/FLOWPanel.jl ]] || {
  echo "ERROR: refusing non-unified deployment target '$REMOTE_DIR'" >&2; exit 2; }

files=(
  examples/rotor_hover_ground_effect.jl
  examples/run_rotor_multi_ground_effect_gpu.slurm.sh
  src/FLOWPanel_gpu_influence.jl src/FLOWPanel_instrumentation.jl
  src/FLOWPanel_simulate.jl src/FLOWPanel_solver.jl src/FLOWPanel_wake.jl
  src/FLOWPanel_warmstart.jl
  scripts/p022g_submit_1r_oge.sh scripts/p022g_submit_2r_oge.sh
  scripts/p022g_submit_4r_oge.sh scripts/p022g_submit_1r_ige.sh
  scripts/p022g_submit_2r_ige.sh scripts/p022g_submit_4r_ige.sh
  scripts/p022g_submit_hr.sh scripts/p022_harvest.py
)
ssh "$REMOTE" "test -d '$REMOTE_DIR/.git'"
rsync -azR --checksum "${files[@]}" "$REMOTE:$REMOTE_DIR/"
for file in "${files[@]}"; do
  local_sum="$(shasum -a 256 "$file" | awk '{print $1}')"
  remote_sum="$(ssh "$REMOTE" "sha256sum '$REMOTE_DIR/$file'" | awk '{print $1}')"
  [[ "$local_sum" == "$remote_sum" ]] || { echo "checksum mismatch: $file" >&2; exit 3; }
done
echo "052b deployment verified at $REMOTE:$REMOTE_DIR"
