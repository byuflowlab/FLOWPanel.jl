# RESET PROMPT — 018 mechanism tests wave 1, MONITOR + SCORE (2026-09-15b, post-launch)

You are a clean-context agent on BRAINSTORM item 018 (DJI-9443 hover
convergence). Wave 1 of the Γ-mechanism isolation tests is **already
launched** — your job is to monitor, banner-verify, score, and report. Do NOT
resubmit or redesign anything that is queued/running. Read, in order:
`mechanism_tests_provenance_20260915.md` (this directory — the authoritative
record of pins, arms, env, and the warmstart incident), then skim
`mechanism_tests_reset_prompt_20260915.md` (the launch charter: task
definitions T1/T2/T4/T5, scoring recipes, verdict criteria, ORC gotchas) and
`gamma_distribution_status_20260915.md` (the signature being tested). Standing
rules: never tune to CT_exp; notebook/ledger writes need Ryan's approval
(both stale since 2026-09-08; backlog listed at the top of the launch
charter, now plus this wave's launch+incident).

## In-flight jobs (submitted 2026-09-15 evening; poll `sacct`, monitors died with the old session)

`sacct -j 13712093,13712094,13712095,13712096,13712097,13712154,13712155 -X --format=JobID,JobName%28,State,Elapsed`

| job | test | run name | notes |
|---|---|---|---|
| 13712093 | T2 merge-cadence | `p018_csarc_n2_nt72_l3p0_3r_srlx_mrg2_g25` | RUNNING at reset; banner ALREADY VERIFIED (NT:72, rlxf:0.16334, merge_every:2, guard=on, DynamicSFS rlxf=0.0025031) |
| 13712154 | T1 chain: gate → fidelity check → T1a | `…_rgate` then `p018_csarc_l3p0_3r_g25_s2` | wrapper `orc:~/p018_mech_tests_20260915/t1chain.slurm.sh`; 9 h wall |
| 13712155 | T1b (afterok:13712154) | `p018_csarc_n2_nt72_l3p0_3r_srlx_g25_s2` | NT72 rev 30→60 restart; only runs if the whole chain (incl. scripted gate check) exits 0 |
| 13712094 | T4a | `p018_csarc_l3p0_3r_csfs_g25` | ConstantSFS Cs=0.14 |
| 13712095 | T4b | `p018_csarc_n2_nt72_l3p0_3r_csfs_g25` | ConstantSFS Cs=0.14 |
| 13712096 | T5a | `p018_csarc_l3p0_3r_nv_g25` | Inviscid (`_nv` case arm) |
| 13712097 | T5b | `p018_csarc_n2_nt72_l3p0_3r_srlx_nv_g25` | Inviscid |

History you will see in sacct (expected, do not re-diagnose): 13712091 first
T1 chain FAILED (warmstart `dsigma2_*` version skew — root-caused and FIXED in
FLOWPanel commit `5cdf058`, see provenance §Pins), 13712092 first T1b
auto-CANCELLED (dependency never satisfied). Also 13711487/88 (022g, another
session) and 13711596 (021) are not yours.

## Your tasks

1. **Banner-verify every job as it starts** (mandatory, ops_reference rule).
   Logs: `orc:~/wt018/FLOWPanel-mechtests/logs/slurm/slurm-fp-018mech-*-<jobid>.out`.
   Expect: T4a/b `SFS=ConstantSFS(Cs=0.14…)` (grep `SFS=`); T5a/b Inviscid /
   `visc:false` and NO CoreSpreading line; T1 chain stage banners
   (`T1CHAIN STAGE n`) + `simulate_warmstart!: resuming from step 1044`
   (gate), `1079` (T1a); T1b `resuming from step 2159`. All: `guard=on`,
   correct NT/rlxf, repo `~/wt018/FLOWPanel-mechtests`. A wrong banner ⇒
   scancel that job immediately, then diagnose.
2. **T1 gate verdict**: when 13712154 finishes (or after its stage 2), grep
   the log for `GATECHECK` and `T1CHAIN`. The scripted check
   (`orc:~/p018_mech_tests_20260915/gate_check.py`) compares per-step force
   (auto-picked thrust component) rgate-vs-original over overlapping steps,
   tol 2e-4 rel (skip first 3), + wake n/Σ|Γ| at step 1079 (informational).
   PASS ⇒ chain proceeds to T1a and T1b releases automatically. FAIL ⇒ T1b
   stays blocked — inspect `max rel dCF`; a marginal miss (few×1e-4, decaying
   seam transient, VTPs are Float32 so restart is NOT replay-exact) is a
   judgment call to bring to Ryan with the trace; a gross miss ⇒ restart
   unreliable ⇒ STOP T1, report to Ryan (fresh 60-rev runs are the fallback,
   needs his sign-off). Note: launcher gate also fails jobs on
   `gpu_gemv=0`/NaN — read the last `GATE:` line to tell WHICH gate tripped.
3. **Score each completed run** (recipes + verdict criteria in the launch
   charter §Scoring; baselines = g25 pair 13704962/13704963):
   per-rev near/far Σ|Γ| time series (`orc:~/p018_gamma_dist_20260915/p018_gamma_ts.py`),
   phase-avg per-bin ratios rev 29→30 (`p018_gamma_dist_avg.py` — for T1
   extensions also do late windows, e.g. rev 44→45, 59→60), CT̄ via
   `python3 scripts/p018_analyze.py m1 --revs A B <run>` from
   `orc:~/projects/FLOWPanel.jl` (one run per call; windows 21 25 / 26 30 /
   21 30; T1 also 31 40 / 41 50 / 51 60). Primary discriminator: NT72's late
   far-wake re-acceleration + NT36/72 far-wake divergence; mechanism
   implicated if the signature shrinks >~50%. T1 verdict: does NT72's far
   wake saturate by rev ~45–60 and the climb shrink in late windows?
4. **Deliverable**: dated status file in this directory
   (`mechanism_tests_status_<date>.md`): three metrics per run side-by-side
   with g25 baselines; verdict per mechanism (merging cadence / SFS dynamics
   / viscous / transient length); wave-2 recommendation (T3 merging-off is
   pre-approved follow-up if T2 is positive). Offer (don't write) notebook
   entries incl. the backlog.

## Environment / provenance facts (already done — do not redo)

- Tags `campaign/p018-mech-tests-20260915`: FLOWPanel `5cdf058`
  (worktree `~/wt018/FLOWPanel-mechtests`; = expguard pin `d5dd772` +
  examples knobs `5d8711c` [MERGE_EVERY, `_nv` arms, banner merge_every] +
  warmstart hasproperty fix), FLOWVPM `7468712` and FastMultipole `3da58a1a`
  (expguard worktrees reused, clean). Env `~/p018wtenv-expguard-gh200`,
  FLOWPanel dev-path re-pointed at the mechtests worktree.
- T4's Cs=0.14 provenance: extracted from g25 VTP per-particle `C` (column 1
  of 3; cols 2–3 are Lagrangian-avg storage — do NOT flatten), phase-avg rev
  29→30: NT36 mean 0.1446 / NT72 0.1308 (|Γ|-wtd ~0.24, ~54% clipped to 0).
  Script `orc:~/p018_gamma_dist_20260915/p018_sfs_C_avg.py`, log `sfs_C_avg.log`.
- Restart is env-only (`RESTART_STEP/NAME/PATH`); `restart_reconstruct_required`
  in metadata is an inert TOML-serialization marker, not a health flag.

## Ops gotchas (compressed; details in launch charter + ops_reference.md)

- `ssh orc 'bash -lc "…"'` for slurm/python; heredoc `ssh orc 'bash -ls' <<'EOF'`
  for anything with quotes/$; MOTD banner glues to line 1 — lead with `echo`.
- mgh = 2 nodes × 1 GH200; ~35 GPU-h queued ⇒ wave completes ~a day out.
  Walls: T1 chain 9h, T1b 16h, NT72 30-rev runs 14h, NT36 5h.
- Never cat VTP/CSV bytes; print summaries. Axial axis is x, R=0.11995 m.
- New runs' full VTP sets are live while running; the sweeper culls landed
  runs to newest-36 — if a completed test run's full set matters (e.g. for
  the Γ time series, which needs many steps), run the Γ scripts PROMPTLY
  after landing and/or harvest to
  `/nobackup/archive/usr/rander39/FLOWPanel_runs/` (g25 originals also live
  in-run-dir + tarred there). `vtk_protect_list.txt` is Ryan's — never write.
- Job logs land in the mechtests worktree `logs/slurm/`, run data in
  `~/projects/FLOWPanel.jl/data/<run>/` (worktree `data` is a symlink).
- p018_analyze reconstructs CT from the force monitor if a wall-killed run
  has no CT CSVs. Judge health from monitors CSVs, exit 0 ≠ health.
- ≤4 local threads; local FLOWPanel/FLOWVPM checkouts are DIRTY with other
  sessions' 021/026 work — never commit or run from them; cluster campaign
  work happens only in the `~/wt018/*` worktrees.
