#!/usr/bin/env bash
#
# BRAINSTORM 020 Phase-2R deployment: stage the corrected frozen-gradient
# geometric integrator (FLOWVPM euler_exp + FLOWPanel plumbing + WAKE_EXPINT
# knob) and the decisive Phase-2R screen arm onto the cluster checkouts,
# md5-verifying every copied file and verifying the already-deployed dirty
# FLOWVPM src set still matches local (Phase 2c stale-cluster-src lesson).
# Nothing is submitted; this only stages files and reports.
#
# Run from the top level of the LOCAL FLOWPanel.jl checkout:
#   bash scripts/deploy_p020.sh

set -euo pipefail

REMOTE="${REMOTE:-orc}"
REMOTE_FP="${REMOTE_FP:-/home/rander39/projects/FLOWPanel.jl}"
REMOTE_VPM="${REMOTE_VPM:-/home/rander39/projects/FLOWVPM.jl}"
LOCAL_VPM="../FLOWVPM.jl"

FP_COPY=(
  "src/FLOWPanel_wake.jl"
  "examples/rotor_hover_pressure_comparison.jl"
  "examples/run_p018_screen_hpc.slurm.sh"
  "test/runtests_unit_simulate.jl"
)
VPM_COPY=(
  "src/FLOWVPM_timeintegration.jl"
  "src/FLOWVPM_viscous.jl"
  # Runaway-safe sparse binning is part of the screen's failure containment;
  # the remote copies must not retain the older dense-grid implementation.
  "src/FLOWVPM_merging.jl"
  "src/FLOWVPM_splitting.jl"
  "test/runtests_expint.jl"
  "test/runtests.jl"
)
# Dirty-in-local FLOWVPM src the cluster already runs: must MATCH local.
VPM_VERIFY=(
  "src/FLOWVPM_subfilterscale_models.jl"
)

md5_local() { md5 -q "$1" 2>/dev/null || md5sum "$1" | awk '{print $1}'; }

echo "== local md5 (FLOWPanel) =="
for f in "${FP_COPY[@]}"; do printf '  %s  %s\n' "$(md5_local "$f")" "$f"; done
echo "== local md5 (FLOWVPM) =="
for f in "${VPM_COPY[@]}" "${VPM_VERIFY[@]}"; do
  printf '  %s  %s\n' "$(md5_local "$LOCAL_VPM/$f")" "$f"
done

echo
echo "== copying =="
for f in "${FP_COPY[@]}"; do
  echo "  FP  -> $f"; scp -q "$f" "$REMOTE:$REMOTE_FP/$f"
done
for f in "${VPM_COPY[@]}"; do
  echo "  VPM -> $f"; scp -q "$LOCAL_VPM/$f" "$REMOTE:$REMOTE_VPM/$f"
done

echo
echo "== remote md5 + prerequisites =="
ssh "$REMOTE" bash -l -s <<REMOTE_EOF
set -euo pipefail
cd '$REMOTE_FP'
mkdir -p logs/slurm
chmod +x examples/run_p018_screen_hpc.slurm.sh
echo '--- FLOWPanel md5 ---'
md5sum ${FP_COPY[*]}
cd '$REMOTE_VPM'
echo '--- FLOWVPM md5 ---'
md5sum ${VPM_COPY[*]} ${VPM_VERIFY[*]}
echo '--- active jobs ---'
squeue -u rander39 -h | wc -l
REMOTE_EOF

echo
echo "Compare md5 blocks line by line; every entry must match local."
echo "Then run the cluster test gate BEFORE submitting:"
echo "  ssh $REMOTE bash -lc \"cd $REMOTE_VPM && julia --project=. -e 'using Test; import FLOWVPM; include(\\\"test/runtests_expint.jl\\\")'\""
echo "Submission (from $REMOTE_FP):"
echo "  sbatch --job-name=fp-020r-geom examples/run_p018_screen_hpc.slurm.sh scr_p020r_geom_s020v"
