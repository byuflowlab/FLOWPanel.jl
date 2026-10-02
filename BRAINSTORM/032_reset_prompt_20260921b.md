# Reset prompt — P2-rerun matrix A1–A5 live; P3 PASS harvested; P2 diagnosis revised (2026-09-21 pm)

You are picking up work in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl`
(branch **fastmultipole**) + sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl`
(branch `flowpanel`, `8d4a3b4`). Read `CLAUDE.md` + the policies it names.
Supersedes `032_reset_prompt_20260921.md` (its NEXT ACTIONS all executed).
Context files, read as needed, do NOT re-derive:
`BRAINSTORM/032_p2rerun_provenance_20260921.md` (THE key file: revised P2
diagnosis, arm matrix, pins, acceptance criteria),
`BRAINSTORM/032_omission_reopen_prep_20260919.md`,
`BRAINSTORM/032_followup_provenance_20260919.md`.

## What happened this session (2026-09-21, established facts)

1. **P3 discriminator PASS.** Warm-start job 13829223 COMPLETED all 324
   steps, finite, normal exit — outlived original death (243) and both
   modern controls (211/~300). Verdict per prep file: root shed CO-DRIVES
   the 020 testbed death; 020 Phase 3 stays STOPPED. GOTCHA found:
   `simulate_warmstart!` REINITIALIZED the monitor CSVs — pre-splice
   history (steps 0–293) is LOST (no backups; original job 13774448 timed
   out before CT diagnostics). The in-run "Phase 2e NOT CONVERGED
   (CT 0.030±141%)" printout is an ARTIFACT (rev-block over empty data);
   real resumed-segment CT ≈ 0.0601, p2p ≈ 0.072, 30 steps only. A citable
   CT needs the Ryan-pre-authorized clean 12 h rerun (Ryan-gated, not
   requested yet). The monitor-clobber is a code finding worth a fix/log.
2. **P2 ignition diagnosis REVISED** (tiebreaker harvest, correct x-axis
   cylinder convention sqrt(y²+z²), VTP totals cross-validated vs
   monitor04): single Γ-ignition, last stable step 1729 (rev 24.0,
   CT 0.133 already elevated), CT −1084.7 at 1750. max|Γ| 0→13.8 (1729,
   rank-1 ~2000× the 99.9pct) →1563 (1740) →6487 (1750). Wake blasted
   upstream+radially; truncation culling (158k→7.3k particles) is a
   SYMPTOM, `3r` exculpated as cause. Old "frozen offender at 3.95R
   outside domain" = coordinate artifact; offender is INSIDE; no culling
   bug. "Merge frenzy trigger @2069" dead (320 steps post-ignition).
   Seeds: deep wake x≈2.8–3.5R incl. root-radius column r_perp/R≈0.17;
   seed σ down to σ/R≈0.0006 (~40× below shed σ) with SIGMA_FLOOR_FRAC=0,
   SIGMA_CEIL=Inf → small-σ Γ-ignition channel, lever = σ guard.
   Secondary note: ignition rev 24.3 is 2 revs after P018_SETTLE_REVS=22
   ended — withdrawal-transient trigger plausible, noted not tested.
3. **P2-rerun matrix SUBMITTED** (Ryan "Go", 2026-09-21): jobs
   **13842791–13842795** (`fp-018gpu-p2rr-a1..a5`), h200 via
   eng/qos=eng/intel, 16 h (A1–A3) / 24 h (A4–A5), from pinned wt
   `orc:/home/rander39/campaigns/p032-rootomit-20260918/FLOWPanel.jl`
   (tag `campaign/p032-rootomit-20260918` = `af92740`; FLOWVPM `8d4a3b4`;
   env `.../p032-rootomit-20260918/env`) — SAME pins as failed P2, so A1
   is exact-P2-plus-guard. All arms: case `p018_csarc_n2_nt72_l3p0`,
   g25 guard (SIGMA_FLOOR_FRAC=0.25, SIGMA_CEIL=0.030),
   TRUNCATION_RADIUS_R=3.0, MAX_PARTICLES=1500000, P018_SETTLE_REVS=22.
   | Arm | Job | Run name | Delta |
   |---|---|---|---|
   | A1 ctrl | 13842791 | `p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15_g25` | exact P2 + guard (SFS_THREELEVEL=true, OMIT=0.15→3/41) |
   | A2 base | 13842792 | `..._3r_srlx_g25_omi1` | SFS_RLXF=0.0025031, OMIT=0.12→1/41 (innermost only, local-probe verified) |
   | A3 | 13842793 | `..._omi1_mo4` | A2 + MERGE_OVERLAP=4 |
   | A4 | 13842794 | `..._omi1_split` | A2 + WAKE_SPLIT_VISCOUS=true, FRACs 0.587/0.73/0.3 |
   | A5 | 13842795 | `..._omi1_split_mo4` | A4 + MERGE_OVERLAP=4 |
   Exact submit lines mirror sacct SubmitLine of 13774449 (in provenance
   file). k=3 σ-cap deliberately NOT reused (retune still owed, 026).

## NEXT ACTIONS (in order)

1. **Banner verification** (an hpc-monitor agent was dispatched at submit
   time; its report may not have landed — redo cheaply if needed): for
   each started arm confirm (a) omission "3/41" (A1) / "1/41" (A2–A5),
   (b) "guard=on" with floor frac 0.25 / ceil 0.030 in the Particle
   diagnostics line, (c) SFS label (threelevel A1; DynamicSFS
   rlxf=0.0025031 A2–A5), (d) MERGE_OVERLAP / WAKE_SPLIT lines (A3–A5),
   (e) mesh 45_185_ct4, NT72, RPM 5400, (f) no ERRORs, GPU gate healthy.
   Logs: `<wt>/logs/slurm/slurm-fp-018gpu-p2rr-a*-<jobid>.out|.err`.
   Any wrong-config arm: report to Ryan immediately (scancel is
   Ryan-gated).
2. **Babysit + harvest A1–A5** against the provenance acceptance:
   survive step 1750, complete 2160 steps, finite CT, bounded monitor04
   max|Γ|/σ². Judge by outputs, never sacct. Delegate to
   hpc-monitor/harvester. Decision tree: A1 survives → σ guard is THE
   rescue lever for the NT72 class → re-raise the NT144 offer with Ryan
   (still gated). A1 dies but A2 survives → confounded (SFS × omission
   width) → escalate to Ryan before more arms. Splitting arms: also watch
   particle count vs 1.5M cap and merged-σ behavior. Remember guard CT
   offset +0.39% — compare slopes vs unguarded history, not levels.
3. **P3 provenance note** still owed when docs unlock: record the splice
   (13774448→13829223, restart step 293, warm-start not bit-exact,
   monitor history lost) in the reopen provenance file.

## BLOCKED on Ryan (do NOT launder via peers)

1. **hpc-storage re-dispatch** — STILL never run: (a) archive
   `scr_p032om15_ctrllg_fs` + `scr_p032om15_explg_fs_oldlaw` (~31 G);
   (b) 212 GB reclaim of OLD `p018_csarc_*_3r_*` runs — **exclude
   anything with `om15`/`omi1`/`srlx_g25` in the name** (the five NEW runs
   match the old glob!). Also reconcile df ~287 G/2.0 T vs earlier
   396.4/400 quota figure.
2. **Task-1 orc cleanup** one-liner (hung PID 1628006, /tmp vpm worktree):
   `! ssh orc 'kill 1628006; cd ~/projects/FLOWVPM.jl && git worktree remove /tmp/rander39_vpm_mergetest --force; rm -rf /tmp/rander39_vpm_testenv /tmp/rander39_vpm_test_20260919d.log'`

## Owed / parked (carried, Ryan-gated)

- Docs commit bundle: 026 ledger line; provenance edits; 032 item-file Log
  update (five 09-19 reopen verdicts + P3 splice + THIS matrix); INDEX.md
  rows (032 refresh + missing `031_quadrupole_panel_farfield.md`); this
  session's new files (`032_p2rerun_provenance_20260921.md`, this prompt).
- ALL pushes (github + orc, incl. deferred FLOWVPM force-with-lease).
- Notebook entries (×4 from 021 + the 032 reopen/rerun arc) — offer only.
- P3 clean 12 h rerun (tie-breaker for a citable CT) — pre-authorized in
  principle, still ask before submitting.
- `simulate_warmstart!` monitor-clobber fix (code change, Ryan-gated).
- k=3 cap retune + merged-σ/clamp telemetry (026 A4); 021 silo cleanup;
  `scr_p026gpuv_split` retry once quiet ≥24 h; 026 feature-A
  move-vs-delete DECLINED; NT144 offer OFF unless A1 rescues.
- 021 entry point: `fgs_acceleration_reset_prompt_20260919e.md` (gated).

## Ground rules

Local ≤4 threads; macOS has NO `timeout` (use poll-and-kill wrappers).
ssh orc has a live ControlMaster socket this session — after a reset it
may need Ryan to run `! ssh orc echo ok` (2FA). Strip MOTD/ANSI. Judge
runs by outputs, never sacct. Delegate monitoring/harvest/scouting to the
repo subagents (`hpc-monitor`, `harvester`, `brainstorm-scout`,
`code-scout`); keep conclusions inline. New submissions/commits/remote
git state/notebook writes Ryan-gated (the A1–A5 submission was explicitly
approved and is done). orc login julia is 1.12 — never let it touch a
campaign Manifest (spack `julia/1.11.7-6bmogfl` for env ops). VTP
particle coordinates: rotor/cylinder axis is X — radial is sqrt(y²+z²),
never sqrt(x²+y²) (this exact mistake produced the withdrawn P2 verdict).
Pre-existing dirty files (018/026 docs, rotor_multi slurm script,
pressure-comparison TOML, BRAINSTORM docs incl. this arc) are expected —
commit only with Ryan's approval.
