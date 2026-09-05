#!/usr/bin/env bash
# Prepare a campaign worktree from an annotated campaign tag (Campaign
# Reproducibility policy, ~/.claude/CLAUDE.md).
#
#   scripts/prep_campaign_worktree.sh <tag> <worktree-dir> [data-root]
#
# Creates <worktree-dir> on branch <tag>-wt from <tag>, then commits one
# campaign-branch change: the tracked data/ fixture directory is replaced by a
# symlink to <data-root> (default ~/projects/FLOWPanel.jl/data) so relative
# "data/..." reads and writes resolve to the consolidated data root. Also
# creates logs/slurm. The resulting worktree has no uncommitted state; run
# outputs land outside it. Julia env setup (Manifest dev-paths pointed at the
# campaign worktrees) is site-specific and stays manual — see the campaign's
# provenance file.
set -euo pipefail

TAG=${1:?usage: prep_campaign_worktree.sh <tag> <worktree-dir> [data-root]}
DIR=${2:?usage: prep_campaign_worktree.sh <tag> <worktree-dir> [data-root]}
DATA_ROOT=${3:-$HOME/projects/FLOWPanel.jl/data}

git rev-parse --verify --quiet "refs/tags/$TAG" >/dev/null \
  || { echo "ERROR: tag '$TAG' not found in $(git rev-parse --show-toplevel)" >&2; exit 1; }

git worktree add -b "${TAG//\//-}-wt" "$DIR" "$TAG"
cd "$DIR"
if [ -e data ] && [ ! -L data ]; then
  git rm -rq data 2>/dev/null || rm -rf data
fi
ln -sfn "$DATA_ROOT" data
git add data
git commit -qm "campaign worktree ($TAG): data/ -> $DATA_ROOT symlink"
mkdir -p logs/slurm

echo "worktree : $DIR"
echo "branch   : ${TAG//\//-}-wt"
echo "tag      : $TAG ($(git rev-parse --short "$TAG^{commit}"))"
echo "HEAD     : $(git rev-parse --short HEAD) (tag + data-symlink commit)"
echo "data     : $(readlink data)"
