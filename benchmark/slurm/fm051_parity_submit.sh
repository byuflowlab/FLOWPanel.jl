#!/bin/bash
# Task 051 stage 2/3 LOCAL DRIVER (Job B) — RUN BY THE USER, never by an
# agent (user submits all cluster jobs; agents only author this script).
#
# Syncs the FLOWPanel working tree (branch fastmultipole, incl. the untracked
# gpu-influence seam + parity harness) plus the unified FLOWVPM/FastMultipole
# trees to the BYU cluster, sets up ~/fm051env (fm048env pattern + FLOWPanel
# + VSPGeom), then submits the full-mode pass-parity job
# (benchmark/slurm/fm051_parity.sh). Pattern: FLOWVPM scripts/fm051_submit.sh.
#
# Cluster layout this creates/expects (resolved 2026-08-22):
#   ~/FLOWPanel-046           synced FLOWPanel tree (this script)
#   ~/FLOWVPM-046             synced FLOWVPM tree   (shared with fm051_submit.sh)
#   ~/FastMultipole-046       synced FastMultipole tree (shared)
#   ~/fm051env                NEW env: dev's the three trees + CUDA + VSPGeom
#                             (fm048env is left untouched for 048/049 replay)
#   ~/projects/FLOWPanel.jl   standing checkout — used READ-ONLY via symlinks:
#     data/p018_cs_f1_l3p4    (383M restart set) -> symlinked into the tree
#     examples/data           (meshes, 45_185_ct4 present) -> symlinked
#
#   bash benchmark/slurm/fm051_parity_submit.sh    # from the FLOWPanel.jl root
set -euo pipefail
REMOTE=orc
FPDIR=FLOWPanel-046
VPMDIR=FLOWVPM-046
FMDIR=FastMultipole-046
ENVDIR='$HOME/fm051env'
VPMLOCAL=../FLOWVPM.jl
FMLOCAL=../FastMultipole
STANDING='$HOME/projects/FLOWPanel.jl'

[ -f benchmark/fm051_pass_parity.jl ] || { echo "run from the FLOWPanel.jl repo root"; exit 1; }

ssh "$REMOTE" "mkdir -p $FPDIR/examples $FPDIR/data $FPDIR/logs/slurm $FPDIR/benchmark/results $VPMDIR $FMDIR"

# FLOWPanel: src/test/benchmark (minus bulky results) + top-level example
# scripts only; meshes and restart data stay in the standing checkout and are
# reached through symlinks (fm051_parity.sh reads, never writes, data/).
rsync -az --delete --exclude .git \
    src test Project.toml \
    "$REMOTE:$FPDIR/"
rsync -az --delete --exclude results \
    benchmark \
    "$REMOTE:$FPDIR/"
rsync -az --include='*.jl' --exclude='*' \
    examples/ \
    "$REMOTE:$FPDIR/examples/"

# unified FLOWVPM + FastMultipole trees (same targets as fm051_submit.sh, so
# Job A and Job B always see identical kernel/seam code)
rsync -az --delete --exclude .git \
    "$VPMLOCAL/src" "$VPMLOCAL/ext" "$VPMLOCAL/test" "$VPMLOCAL/scripts" "$VPMLOCAL/Project.toml" \
    "$REMOTE:$VPMDIR/"
rsync -az --delete --exclude .git \
    "$FMLOCAL/src" "$FMLOCAL/test" "$FMLOCAL/Project.toml" \
    "$REMOTE:$FMDIR/"

# read-only symlinks into the standing checkout (idempotent)
ssh "$REMOTE" "bash -lc 'ln -sfn $STANDING/data/p018_cs_f1_l3p4 \$HOME/$FPDIR/data/p018_cs_f1_l3p4 \
  && ln -sfn $STANDING/examples/data \$HOME/$FPDIR/examples/data'"

# fm051env: fm023env local-toolkit CUDA recipe (compute nodes have no
# internet; the three CUDA JLLs get local=true, precompile deferred to the
# GPU node)
ssh "$REMOTE" 'mkdir -p fm051env && cat > fm051env/LocalPreferences.toml <<EOF
[CUDA_Compiler_jll]
local = "true"

[CUDA_Driver_jll]
local = "true"

[CUDA_Runtime_jll]
local = "true"
EOF'

# login-node env setup; develop/add calls are idempotent
ssh "$REMOTE" "bash -lc 'module load julia/1.11.7-6bmogfl \
  && export JULIA_PKG_PRECOMPILE_AUTO=0 \
  && julia --project=$ENVDIR -e \"using Pkg; Pkg.develop(path=\\\"\$HOME/$FMDIR\\\"); Pkg.develop(path=\\\"\$HOME/$VPMDIR\\\"); Pkg.develop(path=\\\"\$HOME/$FPDIR\\\"); Pkg.add([\\\"CUDA\\\", \\\"VSPGeom\\\", \\\"Test\\\", \\\"Random\\\", \\\"SHA\\\", \\\"Statistics\\\", \\\"StaticArrays\\\", \\\"LinearAlgebra\\\", \\\"Printf\\\", \\\"Profile\\\"]); Pkg.instantiate()\" \
  && cd $FPDIR \
  && sbatch benchmark/slurm/fm051_parity.sh'"

echo "Submitted. Poll with:  ssh orc 'bash -lc \"squeue -u \\\$USER\"'"
