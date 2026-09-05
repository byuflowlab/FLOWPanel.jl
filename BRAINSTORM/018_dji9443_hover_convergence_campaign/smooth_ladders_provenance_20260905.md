# Provenance — smooth-conversion / SFS NT ladders (relaunch 2026-09-05)

First 018 campaign under the worktree policy (Ryan 2026-09-05: silos deprecated,
campaigns run from git worktrees of the unified repos; see Campaign
Reproducibility in ~/.claude/CLAUDE.md).

## Why a relaunch

The first submission (13592695–13592707, 2026-09-05) ran from the rsync'd
silos after a whole-file driver mirror; the mirror pulled the local BRAINSTORM
026 `rk3` kwarg into silos whose src predates it, and every started job died
in ~1 min at `PanelWake` construction (`got unsupported keyword argument
"rk3"`). All 12 were cancelled. This failure mode is exactly what the worktree
policy eliminates.

## Pinned state (exact silo snapshots, byte-verified vs silo src)

Worktrees live under `/home/rander39/wt018/` on ORC; pin commits snapshot the
018-gpu silo trees that ran every completed ladder arm (reference, `_nm`,
fixed-rlxf, `tp`, `tp`+rlxf), so the new arms are code-identical to their
comparators. FLOWPanel src was identical across arches (one pin); FLOWVPM and
FastMultipole carry real per-arch kernel patches (FLOWVPM_fmm_radix.jl;
containers.jl, translate_batched_cuda.jl, translate_batched_resident.jl) and
are pinned per arch.

| Package | Arch | Branch | SHA | Worktree |
|---|---|---|---|---|
| FLOWPanel.jl | h200 + gh200 | `p018pin-fp-20260905` | `8f3ca07ca501aac2a3fe333b94afcff4b0a1bb6a` | `wt018/FLOWPanel-pin-h200`, `wt018/FLOWPanel-pin-gh200` (detached, same SHA) |
| FLOWVPM.jl | h200 | `p018pin-vpm-h200-20260905` | `113dd6660aafcc283397307b81285e922bffaf6c` | `wt018/FLOWVPM-h200` |
| FLOWVPM.jl | gh200 | `p018pin-vpm-gh200-20260905` | `3df427e29bce7dda8ec85ba2b56f1c1cd8b1c2e0` | `wt018/FLOWVPM-gh200` |
| FastMultipole | h200 | `p018pin-fm-h200-20260905` | `1d2e91382fd4630ade1c8ac2ee3bb76bd14d90da` | `wt018/FastMultipole-h200` |
| FastMultipole | gh200 | `p018pin-fm-gh200-20260905` | `50a943db07f7a5ccaf1cf198035333dd0037412a` | `wt018/FastMultipole-gh200` |

The FLOWPanel pin's driver `examples/rotor_hover_pressure_comparison.jl` is the
local commit `593ffc2` version (contains the SIGMA_CHORD_FRACTION smooth-mode
edit, no `rk3` kwarg). FastMultipole silos' `LocalPreferences.toml` copied
into the worktrees uncommitted (env config, recorded here).

## Environments

`~/p018wtenv-h200`, `~/p018wtenv-gh200`: copies of `p018env-<arch>` with the
three dev-paths rewritten to the worktrees (verified in setup output). Depots
unchanged (`~/fm052depot-gh200` for ARM; default for x86). Jobs pass
`P018_REPO_OVERRIDE=$HOME/wt018/FLOWPanel-pin-<arch>` and
`P018_PROJECT_OVERRIDE=$HOME/p018wtenv-<arch>`.

Run outputs land in the consolidated data root
`~/projects/FLOWPanel.jl/data` via each worktree's `data` symlink; Slurm logs
in `~/projects/FLOWPanel.jl/logs/slurm/` (absolute paths).

## Jobs

Recipe per `smooth_ladders_handoff_20260902.md` §1–2 (4 ladders × NT36/72/144,
NT36→mgh, NT72→eng `--constraint=intel` — the in-script `arm` constraint must
be overridden on eng — NT144→m13h). Submitted 2026-09-05 from the worktrees:

| Ladder | NT36 (mgh) | NT72 (eng) | NT144 (m13h) |
|---|---|---|---|
| 1 `_3r_sv` | ~~13592722~~ → 13592738 | 13592723 | 13592724 |
| 2 `_3r_sv_h2p0` | ~~13592725~~ → 13592759 | 13592726 | 13592727 |
| 3 `_3r_sfs3nb` | ~~13592728~~ → 13592760 | 13592729 | 13592730 |
| 4 `_3r_sv_s1p5` | 13592731 | 13592732 | 13592733 |

13592738 is the GPU-path canary (SurfaceVorticityConversion vs device-resident
CuArray pfield, never exercised before).

Post-submission fix (recorded): the initial pin (`8f3ca07c`) inherited the
cluster repo's tracked `data/` fixture CSVs, so the worktree `data` was a real
directory and the consolidated-root symlink landed *inside* it; the three NT36
jobs that started before the fix (13592722/725/728) died at the
DAS_ARC_TABLE existence check. Pin advanced to
`f46c3fe311dc1be407bc1a0a9cfe9025fde4c26e` (both arches), which removes the
tracked `data/` and commits the symlink; the struck-through jobs were
resubmitted from it. NT72/NT144 jobs had not started and read the fixed
worktree at run time.

Smoke gate (precondition 2) PASSED 2026-09-05: `p018_smoke_sv2_local`
completed 6/6 steps; conversion ledger 5 conversions, 8385 particles,
residual_norm 1.19e-17; wake particle VTPs grow 0 → 6751 → 8342 points.
The "fmm! sources or targets empty" warnings were pre-first-shed steps.
