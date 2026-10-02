# Reset prompt — merge-law production merge + omission re-open pipeline (2026-09-19d)

You are picking up work in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (branch
**fastmultipole**, HEAD `af92740` = tag `campaign/p032-rootomit-20260918`)
+ sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl` (branch
`flowpanel`, HEAD `8d4a3b4`). Read `CLAUDE.md` + the policies it names.
Supersedes `032_reset_prompt_20260919c.md` (BOTH its tasks DONE). Read
these records FIRST — do NOT re-derive state:

- `BRAINSTORM/032_followup_provenance_20260919.md` — C-ctrl/C-oldlaw
  Outcomes + verdict (both COMPLETED 467/467; C-ctrl CT 0.0714872 →
  omission is a general rescue incl. ctrl channel; C-oldlaw CT 0.0716405
  inside B15's ±0.038% band → §22 merge redesign unnecessary under
  omission). Run dirs moved to `orc:~/projects/FLOWPanel.jl/data/` +
  symlinks. Notebook 2-line follow-up APPROVED and WRITTEN
  (`journals/20260901.md` end of `# 20260919` 032 section).
- `BRAINSTORM/032_omission_reopen_prep_20260919.md` — Ryan-approved
  (2026-09-19 AskUserQuestion) prep for THREE omission arms: P1 022 hr10
  IGE, P2 018 NT72 `_3r` (sfs3nb recipe), P3 020 Phase-2R discriminator;
  NT144 CONDITIONAL on P2 success. Includes recovered verbatim sbatch
  lines (sacct SubmitLine) + pre-flight results: P2 ready as-is; P3 needs
  ~2-line CPU-screen-launcher override patch; P1 needs omission knob
  ported into `rotor_hover_ground_effect.jl` + launcher overrides + new
  tag. 026 NT144 un-park and feature-A move-vs-delete remain DECLINED.

## TASK 1 — finish the merge-law production merge (IN FLIGHT at reset)

Ryan directive 2026-09-19: **the §22 merge redesign (FLOWVPM `8d4a3b4`)
is production default for all future runs** unless explicitly specified
(memory `feedback_new_merge_law_production.md`). Ryan chose **full merge
of flowpanel into orc's `unified-052`** (resolves the parked
branch-divergence ruling). State at reset:

- Merge DONE and COMMITTED on a temp worktree branch:
  `orc:/tmp/rander39_vpm_mergetest`, branch `mergetest-unified052`,
  merge commit **`d896145`** (flowpanel `8d4a3b4` → `unified-052`
  `3315b22`). Resolution audited hunk-by-hunk: local flowpanel lineage
  won ALL 44 conflict hunks in 11 files (unified-052's side was ports of
  local `8b00dbd`/`d07b3c1` + silo snapshot `4f6e805`); 052-only content
  (launcher consolidation `186bff4`, gh200 dispatcher fix) retained via
  auto-merge; legacy `FLOWVPM_splitting.jl` stays REMOVED (local
  `8b0b70d` ruling); include list = flowpanel's (resolution_split +
  filament_edges, no splitting); zero conflict markers; scripts kept 755.
- FLOWVPM test suite was RUNNING at reset (temp env
  `orc:/tmp/rander39_vpm_testenv` = copy of p032 env with FLOWVPM
  dev-path → the mergetest tree) — the ssh died with the session, so
  ASSUME IT DIED: rerun
  `julia --project=/tmp/rander39_vpm_testenv -t4 -e 'using Pkg; Pkg.test("FLOWVPM")'`
  (login node, ≤4 threads; GPU subtests may skip off-GPU — fine).
- On PASS, and only with the queue quiet
  (`squeue -u rander39` empty — it was at reset): in
  `orc:~/projects/FLOWVPM.jl` (currently on `unified-052`):
  `git merge --ff-only mergetest-unified052`, then remove the temp
  worktree (`git worktree remove /tmp/rander39_vpm_mergetest --force`)
  and `/tmp/rander39_vpm_testenv`. Also refresh orc's STALE `flowpanel`
  ref (at `eaf257c` wip snapshot): tag it first
  (`archive/flowpanel-wip-20260902`), then push local flowpanel
  `--force-with-lease` to `orc` remote (local FLOWVPM has remote `orc` =
  `orc:projects/FLOWVPM.jl`). On FAIL: report to Ryan, leave live
  checkout untouched (it still works, just old-law).
- Untracked cruft in the live checkout (`data/`,
  `src/FLOWVPM_fmm_radix.jl*.pre052_untracked`) is pre-existing — not
  yours.
- ANSWER ALREADY GIVEN to Ryan: the three prepped arms all use the new
  law (p032 env → FLOWVPM 8d4a3b4 worktree); only C-oldlaw's dedicated
  env is old-law, by design.

## TASK 2 — omission re-open pipeline (all submissions Ryan-gated)

Work `032_omission_reopen_prep_20260919.md` in order P3 → P2 → P1:
finish P3/P1 pre-flight (launcher/driver patches are code edits — do
them, but commits + submissions need Ryan), verify on orc what commit
`wt018/FLOWPanel-pin-h200` is at (scout said `f46c3fe`, tag not in local
repo), check the 020 Phase-2R integrator survives in the p032 pins, then
bring Ryan one AskUserQuestion with the ready-to-submit sbatch lines.

## Archiving state (COMPLETED + reported before reset)

- hpc-storage archived 8/10 approved runs (5×026 slate + 3×032),
  verified, freed 80.6 G → /home at **396.4 of 400 G (3.6 G headroom)**.
  Ledger lines appended to `026_.../ledger.md` and the 032 item Log.
  Logs: `orc:~/archiver_worker_p026p032_20260919.log`.
- OWED: the two new arms (scr_p032om15_ctrllg_fs +
  scr_p032om15_explg_fs_oldlaw, ~31 G) were RECENT-HOT (<2 h quiet) and
  correctly refused — re-dispatch hpc-storage for them once quiet ≥2 h
  (last writes ~16:16Z). Real dirs in `~/projects/FLOWPanel.jl/data/`,
  symlinked from the p032 worktree — don't break the symlinks.
- RYAN DECISION PENDING: 212 GB of ARCHIVE-class
  `p018_csarc_*_3r_*` runs (52.6+72.9+87.5 GB, quiet ≥24 h, unprotected)
  is reclaimable but was NOT authorized this cycle; with 3.6 G headroom
  this is the obvious next lever — ask Ryan.
- Local staging dirs `~/scr_p026s9r2_explg_fs_last50steps/` and
  `~/scr_p032om15_explg_fs_steps278_327/` are DELETED (Ryan-approved).

## Owed / parked (carried)

- Ryan-gated commit bundle: 026 ledger line, provenance edits (incl.
  `032_followup_provenance_20260919.md`, the reopen-prep file, this
  file), 032 item-file Log updates. Also INDEX.md rows for 032 outcome
  refresh + missing `031_quadrupole_panel_farfield.md` row.
- github pushes of branches+tags (orc pushes are now partly done via
  Task 1; github still parked).
- k=3 cap retune + merged-σ/clamp telemetry gap (A4); 021 silo cleanup;
  `scr_p026gpuv_split` retry once quiet ≥24 h;
  `data/scratch_p032_smoke{A,B,C,D}` disposable.
- Sweep verdicts for the record (scouts, 2026-09-19): 020 death is
  tail-seated field-coupled runaway (omission = discriminator, not
  rescue); 005/006 ripple is truncation-seated (NO re-open); 017 is
  500k-cap-blocked (unranked, not approved); 018 Phase-16 λ-ladder died
  TIMEOUT+M2 (the adequacy-gate deaths are the `_3r` arms @1550/1876).

## Ground rules

Local ≤4 threads. `ssh orc` needs a live ControlMaster socket (`! ssh
orc echo ok` if 2FA blocks); `bash -lc` for Slurm; strip MOTD/ANSI from
ssh output. All submissions, commits, and notebook writes Ryan-gated.
Pre-existing dirty files (018/026 docs not from this arc, rotor_multi
slurm script, data TOML) are NOT yours. Judge runs by outputs, never
sacct.
