#!/usr/bin/env bash
# Rotor-hover mesh-convergence driver: 40_40 vs 56_57 DJI9443 meshes.
#
#   bash examples/run_rotor_hover_convergence_sweep.sh          # baseline only
#   SWEEP=1 bash examples/run_rotor_hover_convergence_sweep.sh   # + grad_mu / kerneloffset sweep
#
# Phase A writes data/rhc_sweep_baseline/convergence_history.csv (both meshes,
# default settings) so you can judge convergence. Phase B (SWEEP=1) re-runs both
# meshes under several grad_mu reconstructions and kerneloffset_targets values,
# each into its own data/<RUN_NAME>/convergence_history.csv.
#
# Written for stock macOS bash 3.2 (no associative arrays).
set -euo pipefail

cd "$(dirname "$0")/.."   # repo root

RPM=5400
MESHES="40_40 56_57 80_81"

# Mesh file + trailing-edge seed indices (0-based ParaView IDs) per mesh tag,
# copied from examples/rotor_hover_pressure_comparison.jl. Echoes "MSH TE1 TE2".
mesh_spec () {
    case "$1" in
        40_40) echo "dji9443_new_40_40.msh 1614,1574,45 3324,3284,1755" ;;
        56_57) echo "dji9443_56_57.msh 6370,6314,3255 3117,3061,0" ;;
        80_81) echo "dji9443_80_81.msh 12898,12818,6549 6351,6271,3" ;;
        *) echo "unknown mesh tag: $1" >&2; return 1 ;;
    esac
}

# run_one <run_name> <mesh-tag> [KEY=VAL ...extra env...]
run_one () {
    local run_name=$1 mesh=$2; shift 2
    local spec msh te1 te2
    spec=$(mesh_spec "$mesh")
    msh=$(echo "$spec" | cut -d' ' -f1)
    te1=$(echo "$spec" | cut -d' ' -f2)
    te2=$(echo "$spec" | cut -d' ' -f3)
    echo "=== ${run_name} :: ${mesh} :: $* ==="
    env RUN_NAME="$run_name" \
        MSH_FILE="$msh" \
        TE_INDICES_1="$te1" \
        TE_INDICES_2="$te2" \
        RPM="$RPM" \
        "$@" \
        julia --project -t 3 examples/rotor_hover_convergence.jl
}

# ---- Phase A: baseline (default grad_mu, default kerneloffset_targets) --------
for m in $MESHES; do
    run_one rhc_sweep_baseline "$m"
done

echo
echo "Baseline convergence history:"
cat data/rhc_sweep_baseline/convergence_history.csv

# ---- Phase B: sweep (opt-in) -------------------------------------------------
if [ "${SWEEP:-0}" = "1" ]; then
    # grad_mu reconstructions (default kerneloffset_targets).
    # Each entry: "<tag>|<extra env, space separated>"
    GM_CONFIGS="
quad_default|GRAD_MU_BASIS=quad
quad_nogrow|GRAD_MU_BASIS=quad GRAD_MU_QUAD_GROW=false
tri|GRAD_MU_BASIS=tri
tri_robust|GRAD_MU_BASIS=tri GRAD_MU_TRI_ROBUST=true
"
    echo "$GM_CONFIGS" | while IFS='|' read -r tag extra; do
        [ -z "$tag" ] && continue
        for m in $MESHES; do
            # shellcheck disable=SC2086
            run_one "rhc_sweep_gm_${tag}" "$m" $extra
        done
    done

    # kerneloffset_targets sweep (default grad_mu).
    for ko in 1e-2 1e-3 1e-4; do
        for m in $MESHES; do
            run_one "rhc_sweep_ko_${ko}" "$m" "KERNELOFFSET_TARGETS=${ko}"
        done
    done

    echo
    echo "Sweep convergence histories:"
    for d in data/rhc_sweep_gm_* data/rhc_sweep_ko_*; do
        [ -f "$d/convergence_history.csv" ] || continue
        echo "--- $d"
        cat "$d/convergence_history.csv"
    done
fi
