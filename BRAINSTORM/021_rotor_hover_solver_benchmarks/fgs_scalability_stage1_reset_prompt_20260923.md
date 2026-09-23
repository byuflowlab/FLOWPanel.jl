# Reset prompt: 021 FGS scalability Stage 1 — babysit + harvest (2026-09-23, supersedes fgs_scalability_stage1_reset_prompt_20260922d.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 is CLOSED AT PROMOTION; active follow-on = FGS scalability
diagnostic, plan C (`fgs_scalability_diagnostic_plan_20260921c.md`). Read
`CLAUDE.md` first and `agent_policies/HPC.md` before HPC work.

**Completed on 2026-09-23 (this session):**

1. **TASK A CLOSED — both Ryan-approved archive workers VERIFIED COMPLETE**
   (recorded in `archiver_campaigns_support_20260922.md`, updates section
   dated 2026-09-23): pid 537577 freed 154.8 GB (the three RECENT
   p032-rootomit runs, breadcrumbs in place), pid 539720 freed 6.8 GB
   (`scr_p020r_geom_s020v_om15`). VERIFY_FAIL_COUNT=0, STALE_COUNT=0, no
   `.partial` remnants. Home: **312 G used of 400 G cap** (88 G headroom).
   Still above the 300 G escalation threshold, but the RECENT approval
   queue is empty — residual is live/protected runs; nothing sweepable
   without a new Ryan ruling. 032 ledger mirror OFFERED to Ryan, not yet
   written.
2. **Job 13829231_0 (`p021-r4-thread-scaling`) TIMED OUT at its 48 h wall
   ~2026-09-23 06:50 MDT — this is FINE and needs NO resubmit.** Its FINAL
   table was already harvested (per
   `fgs_scalability_stage1_reset_prompt_20260922c.md` §owed: "nothing owed
   from it"); the Stage 1 plan reuses its records as evidence, it is not
   superseded work. At kill it was heartbeating on tail work (j1) beyond
   the captured table. Its run dir under
   `campaigns/p021-thread-scaling-20260919/FLOWPanel.jl` becomes
   archive-eligible for the storage flow once quiet.

## TASK — babysit job 13858983, then harvest (only remaining task)

**Job 13858983** (`p021-fgs-stage1`, m12-3-29, zen3 exclusive 128c/500G,
48 h wall, started 2026-09-22 19:44 MDT). Run dir:
`/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fgs-stage1-13858983`.
Deploy tree `/home/rander39/campaigns/p021-fgs-stage1-20260922/`
(ARCHIVER_SKIP-marked — never touch); pins in
`fgs_scalability_stage1_provenance_20260922.md`.

Progress at last check (2026-09-23 ~06:55 MDT, elapsed 11:11, ~36:49
remaining): **34 STATUS_* files, all ok, zero failures** — ladder 12/12,
ladderdiag 12/12, accepted 3/3, placement 5/5; now inside the **cap
family** (final phase): `cap16-w16-j64-p1` ok, `cap16-w0-j64-p1`
executing. Cap = worker-cap A/B at j=64 (w=16 caps dagteam near-field
sweep workers, w=0 uncapped baseline; 3 pairs with alternating arm order;
launcher writes `cap16_verdict.txt` and runs 6 more cap32 pairs ONLY if
cap16's paired median beats baseline — see
`benchmark/run_r4_fgs_stage1.slurm.sh:209-239`). Expect completion well
inside the wall.

1. Monitor via `hpc-monitor` only (never inline). Judge by outputs
   (STATUS_* files), never sacct; task logs output-buffered — judge
   liveness by output growth, not log mtime. Launcher writes
   `COMPLETED` with `failed_count=N` at the end; individual stage FAILEDs
   are findings, not aborts.
2. On COMPLETED: harvest via `harvester`; then run
   `benchmark/fgs_stage1_analysis.jl <run_dir>` — FIRST real execution,
   expect to debug it against the actual CSVs (only parse-checked).
   Also read `cap16_verdict.txt` for the worker-cap outcome.
3. **Gate**: if the 16→32 plateau or the 32→64 regression fails to
   reproduce, report deltas and STOP — no interventions on a
   non-reproduced effect. Plateau and regression are SEPARATE conclusions.
4. Resume path if the job dies: resubmit with
   `RESUME_FROM_JOB_ID=13858983` in the `--export` list (same
   COLD_PROJECT=`.../p021-fgs-stage1-20260922/env`,
   CAMPAIGN_PINS=`.../pins.toml`, COLD_DATA_ROOT as above, from the deploy
   tree top level); STATUS_*=ok stages skip.

## Owed (carried)

- **032 ledger mirror** of the p032-rootomit archive lines — offered,
  awaiting Ryan's say-so (lines in `archiver_campaigns_support_20260922.md`).
- **Origin push** of branches+tags in all three repos once Ryan re-auths
  GitHub on the laptop (`gh auth login -h github.com`). Local-only commits
  include archiver work (`5c4f7d1`, `6388a92`) and 021 status commits
  through `7a8097b`, plus this session's status-note edits (uncommitted).
- Archiver **T5 pre-existing failure** (exit 9 vs 8) awaits Ryan's ruling.
- fp-018gpu-p2rr-* (032 arms) left the queue ~2026-09-22 08:00 with no run
  dirs found — 032's thread owns this; flag only.
- Standing Ryan-gated ledger unchanged (see
  `fgs_scalability_reset_prompt_20260922.md` §"Standing Ryan-gated
  ledger").

## Traps (prior traps bind; key ones)

- Fixed-work rows: `solved=false`/`eligible=false` BY CONSTRUCTION — gate
  on certified-accepted + iterations==27 + 1e-8 repeat, never filter on
  eligible/solved. Accepted solves run ONE extra fmm+influence+residual
  vs fixed work.
- `diag_*` = −1 marks uninstrumented rows; never pool instrumented with
  uninstrumented. Worker cap does not pin threads to cores (participation
  test, not topology test — caveat, not bug).
- Cap verdict: compare PAIRED medians only (arm order alternates by pair
  to control node drift); a cap16 win localizes the 64-thread loss to
  sweep participation but does NOT distinguish DAG width vs queue
  contention vs memory traffic.
- `ssh orc` needs a live ControlMaster socket (`ssh orc -fN` if cold) AND
  `bash -lc "..."` for slurm/module commands. Local runs NEVER >4 threads.
  Never edit source while a job uses its deployment — rsync manifests
  make deployed trees tamper-evident: any write fails the job at
  `cold_packages`.
- Home-disk figures: quota usage (du, 312 G) ≠ filesystem free space (df,
  325 G) — earlier logs mixed these; cite usage-vs-cap.
- Run-dir names case-sensitive; judge runs by outputs, never sacct.

## House rules (binding)

Monitoring `hpc-monitor`; harvesting `harvester`; storage `hpc-storage`;
notebook writes Ryan-gated (offer, don't write); dated status/provenance
files in BRAINSTORM/021. HPC submission remains approved for STAGE 1 ONLY
(the 13858983 resume path included); Stage 2/3 Ryan-gated. No resubmit of
13829231_0 (moot — data already harvested).
