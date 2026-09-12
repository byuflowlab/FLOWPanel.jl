# BRAINSTORM 021: v9 R4 diagnostics resume prompt (2026-09-12)

Context reset during the first R4 job. This file is the entry point for the
next agent. It supersedes nothing — the governing handoff is still
`fgs_opt_r4_diagnostics_handoff_20260912.md` (read it in full; its
stopping boundary, gates, and diagnostics sequence remain binding).

## Mission (unchanged)

Validate R4 controls, harvest the R4 baselines + inner screen + profiles,
run the bounded near/far tuning pass and finalist confirmation, reprofile the
retained configuration, and hand Ryan validated profiles, tuning results, and
a ranked list of >=5% implementation opportunities. **Do not implement solver
optimizations.** Skip R3 and further R2 runs. Inherited settings are
provisional. No notebook entry without Ryan's approval.

## State at reset

Deployment is DONE and recorded in `fgs_opt_v9_provenance_20260912.md`
(same directory — read it; it has all pins, md5s, and launch settings):

- v9 source tag `campaign/p021-cold-source-20260912-v9` (`721235e`) and exec
  tag `campaign/p021-cold-exec-20260912-v9` (`f03ab18`) pushed to origin.
- ORC worktree `/home/rander39/campaigns/p021-cold-opt-20260912-v9/FLOWPanel.jl`
  (clean, HEAD f03ab18); env + corrected pins.toml in the same root.
- Dependency pins unchanged from v8 (FastMultipole `ef10643`, FLOWVPM
  `05c658f` at `/home/rander39/campaigns/p021-cold-20260910-v1/`).

Jobs:

- **13657038 FAILED** (2m49s): first pins.toml lacked the top-level
  `[packages]` table required by `cold_packages()` in
  `benchmark/fgs_cold_common.jl:266`. Config-only error, fixed; outputs
  retained. Not a harness defect.
- **13657404 RUNNING** — all-stage R4 (`COLD_OPT_RUNG=R4 COLD_OPT_STAGE=all`,
  `COLD_OPT_SCREEN_SET=inner:1,2,5,10` + always-run inner=3 seed,
  `COLD_MIN_REPS=10`, `COLD_OPT_PROFILE_REPS=20`), 12 h walltime, m12 zen3
  exclusive 64 CPU / 500G. Submitted ~08:55, output root
  `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/opt-13657404/`.
  Progress at reset (~3h15m elapsed): parse, precompile, controls j1+j4
  (passed), smoke, baseline-j4-b1, baseline-j64-b1 all complete; screen 4/5
  candidate `status.toml` written, 5th (inner=1, slowest) in progress; 20-rep
  profile stage still ahead. Rough stage pace: fixture+baseline ~45 min each
  stage-group, screen ~15–25 min/candidate.

## Next actions

1. Poll job 13657404 (ssh orc; sacct + stage logs + `COMPLETED` marker, ≥60 s
   between scheduler queries). Judge by outputs, not sacct alone
   ([[feedback_slurm_failed_status_unreliable]]).
2. On completion: validate per handoff step 4–6 (gates, leaf statuses,
   configuration mapping via requested/config TOMLs — not directory order),
   harvest medians (rank by median prepared total time, not min), and pull the
   evidence to a durable local dir
   `BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_opt_evidence_20260912/opt-13657404/`
   (follow the opt-13653852 evidence layout: harvest_summary.md,
   screen_rank.csv, profile text). Delegate scraping to `harvester`.
3. Read the R4 CPU/allocation profiles independently (do NOT assume R2's
   dense-nonself-GEMV attribution carries over; overlapping stack counts are
   not additive shares). Deserialize `.jls` only on a compute allocation.
4. Bounded near/far tuning (handoff step 7): saved top-two bases + leaf
   `{25,50,100,200}`, then P/MAC neighbors, via `COLD_OPT_SCREEN_BASE_FILE`
   (+`COLD_OPT_SCREEN_SET`, screen stage). Retain two by median each stage.
   No full Cartesian sweep.
5. Finalist confirmation: >=10 unprofiled trials, alternating
   baseline/candidate batches (stock launcher baselines alone do not provide
   this). Reprofile the retained configuration (profile driver + calibrated
   CONFIG_FILE, separate clean allocation).
6. Deliver: validated tables, text profiles, ranked >=5% implementation
   opportunities (measured vs hypothesis distinguished), then STOP.

Config-only reruns may reuse clean v9; any harness code change requires a new
pinned generation (v10) — never move v9 tags or edit the running worktree.

## Session gotchas (cost real time — avoid repeats)

- **pins.toml schema**: `[packages.<Name>]` sections with `path`/`tag`/`sha`;
  `cold_packages()` reads `TOML.parsefile(CAMPAIGN_PINS)["packages"]` and
  verifies annotated tag, SHA, clean status, worktree-ness.
- **Don't filter ssh output with `grep -v '^['`** — it silently strips TOML
  section headers (that's what produced the bad pins.toml). Login banner noise
  is `[1m...` tips; filter on `tip:|fortune` instead.
- **Monitor-tool ssh fails in its sandbox** (no ControlMaster access) and
  fabricates "left queue"; poll with plain Bash ssh + ScheduleWakeup instead.
- ORC cannot push to GitHub (https, no askpass): create tags on ORC, then
  `git fetch orc:/home/rander39/projects/FLOWPanel.jl tag <tag>` locally and
  push from the local machine.
- `ssh orc` needs the live ControlMaster socket (`ssh -O check orc`); if cold,
  ask Ryan to run `ssh orc -fN` (2FA). Don't retry into MFA.
- Local v9 worktree: `/private/tmp/flowpanel-cold-opt-20260912-v9` (clean; has
  the launcher and `fgs_cold_common.jl` for reference reads without ssh).
- Availability at submission: m12 had 28 idle fitting zen3 nodes, immediate
  start at 12 h walltime, maxtime 3 d, `--qos=normal`. Reprobe before new
  submissions (`slurm-availability` skill, `--cpus 64 --mem-gb 500`).
- Storage: home ~192 GiB used (df proxy), R2 evidence 84 MiB, no VTK in this
  campaign. Refresh before big new outputs.
