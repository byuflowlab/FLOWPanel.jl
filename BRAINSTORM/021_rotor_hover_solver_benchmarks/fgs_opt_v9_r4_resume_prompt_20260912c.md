# BRAINSTORM 021: v9 R4 diagnostics resume prompt C (2026-09-12, third reset)

Entry point for the next agent. Supersedes
`fgs_opt_v9_r4_resume_prompt_20260912b.md` (its steps 1–2 are DONE; its
mission, Ryan's binding decisions, cost notes, and the "Session gotchas" in
`fgs_opt_v9_r4_resume_prompt_20260912.md` all remain valid — read both).
Governing handoff: `fgs_opt_r4_diagnostics_handoff_20260912.md` (binding:
stopping boundary, gates, diagnostics sequence). Pins and all five job
launches: `fgs_opt_v9_provenance_20260912.md`. Read `~/.claude/CLAUDE.md`,
repo `CLAUDE.md`, `agent_policies/HPC.md`.

## Mission (unchanged)

Finish stage B, run finalist confirmation + reprofile, deliver validated
tables, text profiles, and ranked >=5% implementation opportunities
(measured vs hypothesis). **Do not implement solver optimizations. No
notebook entry without Ryan's approval.** Config-only reruns reuse clean v9
(exec `f03ab18`, worktree
`/home/rander39/campaigns/p021-cold-opt-20260912-v9/FLOWPanel.jl`); any
harness code change needs a new pinned generation (v10).

## State at reset

**Job 13657404 (R4 all-stage) DONE + harvested** (see prompt B). Winner
inner=3, runner-up inner=5; seed baselines j4/b1 17.016 s, j64/b1 11.099 s;
R4 profile: dense nonself GEMV ~77% of snapshots, `influence!` temporaries
(FastMultipole solve.jl:1336) top allocation site, scatter ~6%, leaf solves
~5%, fmm! ~1.3%.

**Job 13660643 (stage A: leaf axis) DONE, validated, harvested, AUDITED.**
Evidence: `fgs_opt_evidence_20260912/opt-13660643/` (`harvest_summary.md`,
`screen_rank.csv` — audited against raw CSVs, mappings verified). All 8
candidates completed/eligible, gates PASS (worst BC rel-L2 6.68e-7, all
certified_fmm, repeat deltas 0). Median prepared j64/b1: i3/l100 11.160 <
i5/l100 11.283 < i3/l50 11.576 < i3/l25 11.776 < i5/l50 11.975 < i5/l25
12.390 < i5/l200 12.761 < i3/l200 12.852. **leaf=100 optimal on the leaf
axis; retained = the two bases themselves.** Recalibrated tolerances came
out bit-identical to job 1 (deterministic calibration). In-job seed
baselines: j4/b1 17.129, j64/b1 10.968. Benign quirks: root
`selected.sha256` cites job-1's selected.toml (it IS the CONFIG_FILE);
`convergence.csv` header is `iteration,residual` (outer iters = data rows);
i3/l100 ran 28 outers here vs 27 in job 1 (same tolerance).

**Stage B (13661797 = i3/l100 base on m12-1-17, 13661798 = i5/l100 base on
m12-1-25) DONE 2026-09-12 ~20:00, validated, rsync-harvested to
`fgs_opt_evidence_20260912/opt-1366179{7,8}/`.** Per job the roster was base
+ P:{6,10} + MAC:{0.3,0.5} one-factor neighbors, self-anchored. Results:
base wins in BOTH jobs (13661797: base 11.014 < P10 11.969 < MAC0.3 14.971,
anchors j4 17.261/j64 11.212; 13661798: base 11.126 < P10 11.792 < MAC0.3
15.376, anchors 16.844/11.081). **P6 and MAC0.5 FAILED calibration in both
jobs** ("FGS staircase has no certified crossing with a decreasing
successor", fgs_cold_common.jl:442 — looser far field can't certify 1e-6 at
R4; findings, outputs retained). All completed candidates pass gates.
Recalibrated tolerances bit-identical across nodes/jobs (deterministic).
**Finalists = i3/l100 (= the seed) and i5/l100. Retained config = seed
P8/MAC0.4/leaf100/inner3 → step-4 reprofile NOT needed** (seed already
profiled at 20 reps in job 13657404). A harvester subagent was building
`screen_rank.csv`+`harvest_summary.md` for both stage-B evidence dirs at
reset — if absent or unaudited, redo/audit against raw CSVs (failed
candidates: config in `requested_config.toml`, error in `status.toml`).

**Confirmation RUNNING at reset (jobs 13663309 on m12-1-17, 13663310 on
m12-1-25; submitted 2026-09-12 20:07, 6 h wall, ~3 h expected).** Env in
provenance file: `screen_profile`, `CONFIG_FILE=confirm_r4_20260912.toml`
(md5 `800aab6a535e91f5aebc0d2b8efd0e88`; configs = [i3/l100 tol
3.479128881193055e-7, i5/l100 tol 2.5712316195808637e-7], timed
sequentially >=10 trials each in one process per baseline stage), trivial
screen roster (`COLD_OPT_SCREEN_BASE_FILE=bases_r4_leaf_20260912.toml`,
`COLD_OPT_SCREEN_SET=leaf:100` → dedupes to the 2 finalists), end profile
stage profiles BOTH confirm configs at stock 10 reps. Outputs:
`.../opt-1366330{9,10}/` (opt-13663309, opt-13663310). DO NOT resubmit; the
old poller died with this reset — start your own (plain Bash ssh loop,
`COMPLETED` marker + sacct terminal fallback, >=60 s between scheduler
queries; judge by outputs, not sacct).

## Verified harness facts (save re-derivation)

- `cold_run` (fgs_cold_common.jl:654) iterates ALL configs of a CONFIG_FILE
  sequentially in one process, for verify/baseline AND profile paths. So a
  multi-config `confirm.toml` times seed→C1→C2 as consecutive >=10-trial
  batches in one matched process, and the launcher's end profile stage
  profiles EVERY config in the file (~20–35 min each — budget walltime).
  Duplicate cold_ids collide on dir names → distinct configs only. Verify
  stage stops at first failed candidate.
- `SCREEN_SET` multi-axis syntax `"P:6,10;MAC:0.3,0.5"` confirmed
  (`cold_screen_axes`); saved-base tolerance resets to 0 for base AND
  neighbors (recalibration everywhere).
- sbatch needs a login shell: `ssh orc 'bash -lc "... sbatch ..."'`;
  `--test-only` first. zsh gotcha: don't use unquoted `===` echo separators.

## Next actions

1. Poll 13663309+13663310 to completion. Validate gates per handoff (BC
   rel-L2 <=1e-6, certified evaluator or explicit direct fallback, repeat
   agreement <=1e-8, convergence; map configs via each dir's `config.toml`,
   never dir order; two configs per baseline stage — both must be checked).
   rsync each output root to `fgs_opt_evidence_20260912/opt-<jobid>/`;
   delegate summary tables to `harvester` (sonnet; opt-13660643 files as
   templates) and AUDIT against raw CSVs. The confirmation comparison is
   i3/l100 vs i5/l100 medians within each job (two independent
   replications); repeat ambiguous comparisons; never loosen gates.
2. (Reprofile is NOT needed — retained config = the seed, already profiled
   at 20 reps in opt-13657404.)
3. Harvest confirmation the same way; then deliver the final
   diagnostics package: validated config/timing tables, text profiles,
   ranked >=5% implementation opportunities with supporting R4
   stacks/allocations (measured vs hypothesis), extra measurements needed.
   Measured picture the ranking must reflect: dense nonself GEMV ~77%,
   influence! allocation temporaries (solve.jl:1336), scatter ~6%, leaf
   solves ~5%, FMM ~1.3%. Then STOP. Offer a notebook entry (ask detail
   level); do not write it without approval.

## Cost/scale notes

- Stage-B pace: baselines ~14 min each, ~20 min/candidate, ~35 min profile.
- Storage: campaign outputs ~150 MB/job, no VTK; home ~192 GiB used.
- All ssh via `orc` alias (needs live ControlMaster; if cold ask Ryan to run
  `ssh orc -fN`; never retry into MFA).
