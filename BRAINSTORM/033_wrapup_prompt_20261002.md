# 033 wrap-up prompt (2026-10-02) — verify branches, commits, prep pushes

Context: BRAINSTORM 033 session 7 finished 2026-10-01 — B-I2 and B-I3 are
complete, clear-context reviewed, and committed. Entry point for the full
state: `BRAINSTORM/033_reset_prompt_20261001.md`. This session is a short
verification/wrap-up only — do NOT start new measurements.

## 1. Verify branches and HEADs (expected state)

| repo | path | expected branch | expected HEAD |
|---|---|---|---|
| FLOWPanel.jl | ~/Dropbox/research/projects/FLOWPanel.jl | `fastmultipole` | this prompt's commit, directly atop `d9d05ea` (033 session 7) |
| FastMultipole | ~/Dropbox/research/projects/FastMultipole | `flowpanel-20260817` | `6456c221` (dynamic per-leaf repack) |
| FLOWVPM.jl | ~/Dropbox/research/projects/FLOWVPM.jl | `flowpanel` | `eebb984` (unchanged this session; tags only) |

Run `git branch --show-current` + `git log --oneline -3` in each. Caveat:
a session-start snapshot once reported FLOWPanel on branch `flowpanel`,
but all session-7 commits landed on `fastmultipole` — confirm which branch
is checked out and that `d9d05ea` is on it before pushing anything.

## 2. Verify everything is committed

`git status -s` in each repo. Expected leftovers in FLOWPanel (unrelated,
leave alone): untracked `BRAINSTORM/018_dji9443_hover_convergence_campaign/
ntladder_r4_plan_20261001.md` and `ntladder_reset_prompt_20261001.md`
(item 018, not ours to commit). FastMultipole and FLOWVPM should be clean.
Session-7 commits to confirm present:
- FLOWPanel `fastmultipole`: `1c03f68`, `6e15e6a`, `143c8d5`, `d9d05ea`
- FastMultipole `flowpanel-20260817`: `33c3e20b`, `29fcdcb2`, `6456c221`

## 3. Tags (exist on local + orc; GitHub origin pushes were Ryan-gated)

`campaign/p033-bi2prof-20261001` (+`...-20261001b` in FLOWPanel only),
`campaign/p033-bi2-20261001`, `campaign/p033-bi2b-20261001` — in all three
repos (verify with `git tag -l 'campaign/p033-*'`).

## 4. Push commands (hand to Ryan — https creds were unavailable to agents)

```bash
# FLOWPanel.jl
cd ~/Dropbox/research/projects/FLOWPanel.jl
git push origin fastmultipole && git push origin --tags

# FastMultipole
cd ~/Dropbox/research/projects/FastMultipole
git push origin flowpanel-20260817 && git push origin --tags

# FLOWVPM.jl (tags only; branch unchanged)
cd ~/Dropbox/research/projects/FLOWVPM.jl
git push origin --tags
```

Prefer pushing only the campaign tags instead of `--tags` if the tag
namespace is noisy: `git push origin 'refs/tags/campaign/p033-*'`.

Known divergence (do NOT force-push): the `orc` remotes' branches have
commits not in the local checkouts (orc-side merges, e.g. unified-052) —
branch pushes to orc fast-forward-fail and are unnecessary; the campaign
tags are already on orc, which is all the worktrees need.
