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

## Apply cycle report (hpc-storage, 2026-09-22 ~20:19 UTC)

- Dry run over 10 checkouts discovered via the new `--all-checkouts`
  (main clone, flowpanel-021, wt052, 7 campaigns/*/* trees);
  `ARCHIVER_SKIP` on p021-fgs-stage1-20260922 respected. STALE_COUNT=0,
  VERIFY_FAIL_COUNT=0, CHECKOUT_LOCKED_COUNT=0 system-wide;
  `PROTECTED p022lg_hr10` skipped; no LIVE runs in the candidate set
  (13829231_0 and 13858983 untouched).
- **Applied (detached, verification pending at session reset)**:
  `campaigns_p032-rootomit-20260918_FLOWPanel.jl/p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15`
  (23,764 files, 39,659 MB src, keep steps [2155–2159], expected free
  ~39,612 MB). Worker pid 500292 on login03, log
  `/home/rander39/archiver_worker_p032rootomit_20260922_201936.log`;
  confirmed alive and progressing (`.partial` growing). Completion check:
  `grep -E 'ARCHIVE_MODE=|TOTAL_FREED_MB=|STALE_COUNT=|VERIFY_FAIL_COUNT='`
  on that log + `ARCHIVED.txt` breadcrumb in the run dir.
- Archive quota at cycle start: 4.549 T / 93,402 files (of 20 TiB / 1 M).
- Full dry-run plan: `/home/rander39/archiver_dryrun_all_20260922_201506.log`.

### Awaiting Ryan (nothing else this cycle can safely do)

| run | class | VTK | quiet |
|---|---|---|---|
| p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15_g25 (p032-rootomit) | RECENT | 38,157 MB | 17 h |
| p018_csarc_n2_nt72_l3p0_3r_srlx_g25_omi1 (p032-rootomit) | RECENT | 61,283 MB | 16 h |
| p018_csarc_n2_nt72_l3p0_3r_srlx_g25_omi1_mo4 (p032-rootomit) | RECENT | 55,788 MB | 16 h |
| scr_p020r_geom_s020v_om15 (projects) | now ARCHIVE-eligible (~36 h) | 7,010 MB | held per standing instruction |

Approval command (three p032 runs, 155.2 GB):
`./scripts/run_archiver.sh --include-recent --only p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15_g25,p018_csarc_n2_nt72_l3p0_3r_srlx_g25_omi1,p018_csarc_n2_nt72_l3p0_3r_srlx_g25_omi1_mo4 --root /home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl --apply`

After the in-flight archive completes, home lands ~326 GiB — above the
300 G escalation threshold, but the residual is entirely the RECENT
approval queue: it needs Ryan's decision, not the sweeper ladder.

### Ledger line

```
2026-09-22 20:19 UTC — hpc-storage cycle: home 365.0G→(pending, worker still tarring) of 400G cap.
Archived (apply, detached, pid 500292): campaigns_p032-rootomit-20260918_FLOWPanel.jl/p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15
  (23764 files, 39659 MB src, keep steps [2155-2159], expected free ~39612 MB) — tarball verification pending.
Left untouched per explicit instruction: projects_FLOWPanel.jl/scr_p020r_geom_s020v_om15 (7010 MB, now ARCHIVE-eligible, still unapproved).
Approval queue (RECENT, not archived): p032-rootomit {sfs3nb_om15_g25 38157MB/17h, srlx_g25_omi1 61283MB/16h, srlx_g25_omi1_mo4 55788MB/16h} = 155.2GB, awaiting Ryan.
STALE_COUNT=0, VERIFY_FAIL_COUNT=0 (10 checkouts, --all-checkouts). Archive quota 4.549T / 93,402 files at start.
```

(These runs are 032 arms — 032's thread may want this line mirrored into
its ledger; not written there without Ryan's say-so.)

## Update (2026-09-22, later session): first archive VERIFIED, approvals granted

- Worker pid 500292 finished cleanly:
  `p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15` archived and stripped,
  `TOTAL_FREED_MB=39612`, `VERIFY_FAIL_COUNT=0`, `STALE_COUNT=0`,
  `DF_AFTER=342G`.
- **Ryan APPROVED (2026-09-22): (1)** archiving the three RECENT
  p032-rootomit runs (`..._sfs3nb_om15_g25`, `..._srlx_g25_omi1`,
  `..._srlx_g25_omi1_mo4`; 155.2 GB) via
  `--include-recent --only ...`, and **(2)** archiving
  `scr_p020r_geom_s020v_om15` (7 GB, projects checkout). Both launched
  as detached workers via `hpc-storage`:

| worker | pid | root | targets | log |
|---|---|---|---|---|
| 1 | 537577 | campaigns/p032-rootomit-20260918/FLOWPanel.jl | the three RECENT-approved runs (155.2 GB) | `/home/rander39/archiver_worker_p032rootomit_20260922_202851.log` |
| 2 | 539720 | projects/FLOWPanel.jl | scr_p020r_geom_s020v_om15 (~7 GB) | `/home/rander39/archiver_worker_scrp020r_20260922_202919.log` |

- Both confirmed alive at launch: correct banner (`keep=5 ... apply=true`),
  growing `.partial` tarballs on the archive side; log 1 shows
  `RECENT-APPROVED ..._sfs3nb_om15_g25 quiet=17h -- archiving on Ryan's
  explicit --only`. Pre-launch squeue check: only 13829231_0 and 13858983
  running, neither matching any target; the p021-fgs-stage1-20260922
  checkout untouched.
- Archive quota at this launch: 4.591 T / 93,407 files (no quota limit
  set on /nobackup/archive per lfs quota).
- **Completion is NOT yet verified** (61 GB-class tars run long). A clean
  finish per log = `ARCHIVE_MODE=APPLY`, `VERIFY_FAIL_COUNT=0`,
  `STALE_COUNT=0`, plus per-run `ARCHIVE ... freed=...MB` lines. Expected
  end state: home ≈ 342 G − (~155 G + ~7 G − kept-step residue) ≈ well
  under 200 G.

### Ledger line (launch, second wave)

```
2026-09-22 20:29 UTC — Ryan-approved archiving, two detached workers on orc (run_archiver.sh, keep=5, apply=true):
  (1) pid 537577, root campaigns/p032-rootomit-20260918/FLOWPanel.jl: RECENT-approved p018_csarc_n2_nt72_l3p0_3r_{sfs3nb_om15_g25,srlx_g25_omi1,srlx_g25_omi1_mo4} (155.2 GB source, quiet 16-17h, approved per Ryan's explicit --only/--include-recent);
  (2) pid 539720, root projects/FLOWPanel.jl: ARCHIVE-class scr_p020r_geom_s020v_om15 (~7 GB, quiet >=24h, Ryan-approved).
  Both detached (nohup, nice 19, ZSTD_THREADS=2), logs on /home, alive and writing .partial tarballs at check; neither target matched live squeue (13829231_0, 13858983); p021-fgs-stage1-20260922 checkout untouched. Completion/verify PENDING — check logs before claiming freed MB.
```
