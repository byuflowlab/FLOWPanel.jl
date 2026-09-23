# Archiver support for /home/rander39/campaigns (2026-09-22)

Ryan-approved 2026-09-22 (ruling 1 in
`fgs_scalability_stage1_reset_prompt_20260922c.md`): teach the storage
archiver about `/home/rander39/campaigns`, which held ~212 G (57% of home
usage) invisible to `run_archiver.sh --all-checkouts`.

## Survey (orc, 2026-09-22)

- `/home/rander39`: 365 G of the 400 G cap; campaigns/ = 212 G.
- Dominant: `p032-rootomit-20260918` = **193 G** — four p018 rotor runs
  with full VTK under its deployed `FLOWPanel.jl/data/` (60/55/39/38 G).
- All 18 deployed FLOWPanel.jl trees in campaigns pass the archiver's
  per-checkout gate (Project.toml with FLOWPanel UUID + `data/`).
- p032-rootomit, p032-reopen, and p026-derisk carry **per-run symlinks**
  inside `data/` pointing at shared reference runs in the projects data
  root; p021/052b058 campaign trees have real dirs only.

## Changes (commits `5c4f7d1` + `6388a92` on `fastmultipole`)

1. `CHECKOUT_GLOBS` default extended with `$HOME/campaigns/*/*` — deployed
   trees sit three levels down and the old two-level globs never reached
   them.
2. New `ARCHIVER_SKIP` marker file: present in a checkout or its parent
   campaign dir → `--all-checkouts` discovery skips it (`SKIP-MARKED`
   line). Explicit `--root` overrides (a human decision). Marker placed on
   `/home/rander39/campaigns/p021-fgs-stage1-20260922` (code+env for
   queued job 13858983 — must never be archived).
3. New `ALIAS-RUN` guard in the run loop: a run entry that is itself a
   symlink is skipped — archiving it would file the same run under a
   second slug and **delete VTK through the symlink** in a tree the pass
   never classified. Per-run analogue of the checkout-level `ALIAS-SKIP`
   dedup. This hazard only became reachable once campaigns/*/* was
   discoverable, so it ships in the same change.

## Verification

- `scripts/tests/run_archiver_test.sh`: new T13b (11 checks) — discovery
  of a campaigns/<c>/FLOWPanel.jl tree, symlinked-`data/` alias dedup,
  both marker locations, `--root` override, and the ALIAS-RUN case
  (symlinked run skipped on `--apply`, not tarred, no deletion through
  the link). All pass.
- **Pre-existing suite failure noted for Ryan**: T5 ("resume-delete
  refuses when bytes disagree") expects exit 8 but gets 9 — fails
  identically on the committed pre-edit archiver (verified via stash), so
  it is not from this change. Left unfixed (archiver semantics change =
  Ryan-gated).
- Both scripts mirrored to
  `/home/rander39/projects/FLOWPanel.jl/scripts/` with md5sum verified on
  both sides.

## Follow-up

- `hpc-storage` apply cycle launched 2026-09-22 (this session) — report
  appended when it returns.
- Still NOT approved (left, re-reported): RECENT
  `scr_p020r_geom_s020v_om15` (7 GB).
