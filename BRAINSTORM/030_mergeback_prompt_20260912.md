# 030 — merge-back prompt (two-stage), 2026-09-12

Close out BRAINSTORM 030: two-stage merge-back of the shipped Phases 1+2.
Read FIRST: the 030 item `BRAINSTORM/030_generic_influence_block_assembly.md`
in the live FLOWPanel checkout (`~/Dropbox/research/projects/FLOWPanel.jl`) —
especially the Log entries of 2026-09-10 (sessions 1-3) and 2026-09-11
(retirement) — then `agent_policies/WORKFLOW.md` and `agent_policies/TESTING.md`.
Ryan has authorized merge-back (2026-09-12); Phases 3+3b are retired (already
reverted on the dev branches — merge as-is, do NOT resurrect them).

Branches: dev = `030-block-assembly` in `~/wt030/FastMultipole` (tip
`00bd5f5d`; payload `dd70fb19` hook+probe+calloc Matrices+grav opt-in,
`fff72d29` speedup) and in `~/wt030/FLOWPanel.jl` (tip `15bbbc6`; payload
`6798197` AbstractBody opt-in + 24 tests; 031 staging files ride along).
Production = the live FastMultipole checkout's working branch
(`~/Dropbox/research/projects/FastMultipole`) and the live FLOWPanel
`fastmultipole` branch. Identify the exact production branch names/SHAs from
the live checkouts before starting and record them.

## Stage 1 — merge production INTO the dev branches (wt030 worktrees only)

Fetch the live checkouts from inside each worktree (the worktrees were cut
from the live repos, so `git fetch <path-to-live-checkout> <branch>` works
without adding remotes; adding a path remote is fine too), then in each wt030
worktree `git merge <production-tip>`. Production has moved since the
worktrees were cut (e.g. 052e work in FastMultipole, post-8dce66c commits in
FLOWPanel), so expect real conflicts; resolve minimally and cite every
conflicted file and resolution in your report. Do NOT touch the live
checkouts in this stage.

## Stage 1 validation — all in wt030, ALL green before stage 2

1. FastMultipole cache/FGS driver (runtests prelude
   gravitational/vortex/vortex_filament/panels + `nearfield_cache_test.jl` +
   `solve_test.jl` + `transform_solver_test.jl` + `fgs_coloring_test.jl`,
   run `julia --project=$HOME/wt030/FLOWPanel.jl -t 4 <driver>`; was
   597299/597299 pre-merge — re-derive the merged total. Driver gotcha:
   `_rodrigues` lives in `transform_tree_test.jl`; define it in the driver
   if that file isn't included).
2. FLOWPanel `julia --project=. -t 4 test/runtests_unit_solver.jl` (was
   **Solvers | 489 489** pre-merge).
3. ONE full broad regression:
   `julia --project -e 'include("test/runtests.jl")'` — this has never run
   against the 030 code; budget for it and pipe output to a FILE (no
   `timeout` binary on this Mac — background watchdog pattern).

Manifest note: make sure the wt030 FLOWPanel environment resolves
FastMultipole to the wt030 FastMultipole worktree during validation, so you
test the merged PAIR together. Machine: 8 cores/16 GiB, max 4 threads, check
`ps aux | grep julia` for competing agent jobs before running anything. If
validation fails, fix forward on the dev branches (or stop and report if the
failure traces to production-side changes you don't understand) — production
stays untouched either way.

## Stage 2 — merge validated dev branches into production

Only after stage 1 is fully green. Preconditions — verify immediately before
merging, STOP AND REPORT instead of merging if any fails:
1. No queued/running job uses either live checkout (ask Ryan about HPC
   chains if unsure).
2. Uncommitted changes in live FastMultipole belong to the 052e agent — they
   must be committed/stashed by their owner, not you.
3. The live FLOWPanel worktree is dirty with other agents' in-flight edits —
   if the merge would conflict with uncommitted local modifications, stop
   and report.

Because production was already merged into dev in stage 1, the stage-2
merges should be conflict-free (near-fast-forward). Any conflict here means
production moved again during your session: re-run stage 1 on the new tip
rather than resolving conflicts directly in production. After merging, run
the FastMultipole driver and `runtests_unit_solver.jl` once more in the live
checkouts as a smoke check (the full `runtests.jl` already validated the
identical tree in stage 1).

Merge with `git merge` (no rebase, keep history). Commit but do NOT push.

## Pending items — ASK RYAN before acting

- (a) Optional upstream `tree.jl` guard: depth/degeneracy stop for
  target-tree subdivision so coincident targets fail loudly instead of
  hanging (gotcha documented in the 030 Log).
- (b) Notebook entry for 030 (propose content + verbosity per the notebook
  policy — none has ever been written for this item).
- (c) Deleting the `~/wt030` worktrees after a verified merge.
- Item approval checkboxes are Ryan's alone. Phase 4 (GPU seam) stays
  parked — do not start. Item 031 (quadrupole exploration) is staged but NOT
  active — do not start it.

## Wrap-up

Append ONE dated Log entry to the live 030 item at session end recording the
production branch names, merge SHAs (both stages), and verbatim test totals.
Final report: merge commits, conflicts and resolutions, test totals,
deviations.

## Context snapshot (do not re-derive)

- 030 history: Phases 1+2 SHIPPED (near-field cache hook, 2.90x serial on the
  cheap-kernel build; attribution and gotchas in the 030 Log and
  `030_reset_prompt_20260910b.md`). Phase 3 (route-A FGS migration) and
  Phase 3b (far-pair dipole prototype) were built, measured, and RETIRED by
  Ryan 2026-09-11: FGS builds are projection-bookkeeping-bound (hook
  0.99-1.02x, 1.9x SLOWER for NK=2 bodies) and the dipole far field is both
  too inaccurate (error 0.14/eta^2) and rarely applicable on shedding
  fixtures. Reverts: FLOWPanel `0c444a4`+`3fdf5af`, FastMultipole `00bd5f5d`;
  prototypes recoverable from `bc139c3`/`36428aa`/`417489d5`.
- Successor exploration staged as item 031
  (`BRAINSTORM/031_quadrupole_panel_farfield.md`, live checkout + copy on the
  dev branch) with entry prompt `031_reset_prompt_20260911.md` — parked until
  Ryan activates it.
