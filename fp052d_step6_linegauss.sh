#!/usr/bin/env bash
#SBATCH --job-name=fp052d6lg
#SBATCH --gpus=h200:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=00:40:00
#SBATCH --output=%x-%j.out
# 052d Step-6: LineGauss production switch validation. Re-run the 5e-style
# xverify short sim (NREVS=0.05) with FLOWPANEL_FILAMENT_REG=linegauss so both
# host and device use reg 4 (LineGauss). Expect the panel_cross_xverify relU
# trajectory to drop from ~9e-4 to ~5e-5 (cross-pass truncation) at np~3544.
# (Fix decided by Ryan 2026-08-29; env var scoped to this job — the shared
# fm052_common.sh production array is untouched.)
set -eo pipefail
source /etc/profile
module load cuda julia/1.11.7-6bmogfl
echo "=== node: $(hostname)"; nvidia-smi -L

FMDIR="${FP052D_FMDIR:-$HOME/FastMultipole-052-h200}"
FPDIR="${FP052D_FPDIR:-$HOME/FLOWPanel-052-h200}"
ENVDIR="${FP052D_ENV:-$HOME/fm052env-h200}"
export FASTMULTIPOLE_FORCE_CUDA_LOAD=1
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

cd "$FPDIR"
DUMPDIR="$FMDIR/relU_dumps_6lg_${SLURM_JOB_ID}"
XLOG="fp052d6lg_xverify_${SLURM_JOB_ID}.log"
source "$HOME/FLOWVPM-052-h200/scripts/fm052_common.sh" || true
set +e
timeout 25m env "${FM052_PRODUCTION_ENV[@]}" "${FM052_GPU_ENV[@]}" \
  FLOWPANEL_FILAMENT_REG=linegauss \
  PANEL_INFLUENCE_FMM=1 PANEL_INFLUENCE_FMM_DEVICE=1 \
  PANEL_INFLUENCE_FMM_XVERIFY=1 \
  PANEL_FMM_DUMP_DIR="$DUMPDIR" PANEL_FMM_DUMP_NP=3500 \
  NREVS=0.05 RUN_NAME="fp052d6lg_xverify_${SLURM_JOB_ID}" \
  julia --project="$ENVDIR" --threads="$JULIA_NUM_THREADS" \
  examples/rotor_hover_pressure_comparison.jl > "$XLOG" 2>&1
xrc=$?
set -e
echo "xverify run exit code: $xrc (124 = timeout, OK if xverify lines landed)"
grep -E "filament|linegauss|LineGauss" "$XLOG" | head -5
grep -E "panel_cross_xverify|panel_fmm_dump|panel_cross_eval_dump|panel_cross config" "$XLOG" | tail -20
ls -la "$DUMPDIR" 2>/dev/null || true
NX=$(grep -c "panel_cross_xverify" "$XLOG" || true)
echo "xverify lines: $NX"
[ "$NX" -ge 1 ] && echo "fp052d step-6 complete: XVERIFY LINES PRESENT" \
  || { echo "fp052d step-6: NO XVERIFY OUTPUT"; tail -30 "$XLOG"; exit 1; }
