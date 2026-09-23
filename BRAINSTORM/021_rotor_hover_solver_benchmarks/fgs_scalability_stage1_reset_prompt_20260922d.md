# Reset prompt: 021 FGS scalability Stage 1 — babysit + harvest (2026-09-22, supersedes fgs_scalability_stage1_reset_prompt_20260922c.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 is CLOSED AT PROMOTION; active follow-on = FGS scalability
diagnostic, plan C (`fgs_scalability_diagnostic_plan_20260921c.md`). Read
`CLAUDE.md` first and `agent_policies/HPC.md` before HPC work.

**All three tasks of the 20260922c prompt were EXECUTED on 2026-09-22:**

1. **Rsync deployment DONE.** Tagged triple
   (`campaign/p021-fgs-stage1-20260922`: FLOWPanel `c8621af`, FastMultipole
   `649405cb`, FLOWVPM `8d4a3b4`) shipped via `git archive` into fresh dirs
   under `/home/rander39/campaigns/p021-fgs-stage1-20260922/`. Per-repo
   sha256 manifests verified on orc (`sha256sum --quiet -c` clean in all
   three trees; manifest self-hashes match both sides). `pins.toml` in
   `deployment = "rsync"` mode; env built with julia/1.11.7-6bmogfl,
   Manifest dev-paths confirmed at the deploy dirs. R4 meshes present.
   Full record: `fgs_scalability_stage1_provenance_20260922.md`
   (deployment checklist all ticked except origin push).
2. **Stage 1 SUBMITTED: job 13858983** (m12 zen3 exclusive 128c/500G,
   `--qos=normal`, 48 h wall; submitted 2026-09-22 ~19:50 MDT, estimated
   start 2026-09-23T03:27Z). Run dir:
   `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/fgs-stage1-13858983`.
   Submitted from the deploy tree top level with
   COLD_PROJECT=`.../p021-fgs-stage1-20260922/env`,
   CAMPAIGN_PINS=`.../pins.toml`, COLD_DATA_ROOT as above.
3. **Archiver campaigns support DONE** (Ryan-approved): commits `5c4f7d1`
   (campaigns/*/* discovery + `ARCHIVER_SKIP` marker) and `6388a92`
   (`ALIAS-RUN`: per-run symlinks inside `data/` skipped — they exist in
   p032-rootomit/p032-reopen/p026-derisk and would otherwise be
   double-archived and deleted-through). Test suite extended (T13b, 11
   checks, green); both scripts mirrored to
   `/home/rander39/projects/FLOWPanel.jl/scripts/` md5-verified.
   `ARCHIVER_SKIP` marker placed on the p021-fgs-stage1-20260922 deploy
   dir. Note: `archiver_campaigns_support_20260922.md`.
   **Pre-existing test failure flagged to Ryan, NOT fixed**: T5
   (resume-delete byte-mismatch) exits 9 where the test expects 8 —
   identical before the edits; archiver semantics are Ryan-gated.

## TASK A — collect the hpc-storage apply cycle (launched, outcome UNKNOWN)

An `hpc-storage` apply cycle was launched at session end 2026-09-22 (home
was at 365 G / 400 G; campaigns/ = 212 G, dominated by
`p032-rootomit-20260918` = 193 G of VTK in four p018 rotor runs). The
session reset before its report landed.

- Via `hpc-monitor`: look for detached worker logs
  `/home/rander39/archiver_worker_*_2026092*.log` (canonical launch in
  HPC.md) and live `run_archiver.sh` processes; collect every worker's
  final `STALE_COUNT=`, `VERIFY_FAIL_COUNT=`, and exit status. Do NOT
  report success without all of them. If no logs/processes exist, the
  cycle may not have launched — run a fresh `hpc-storage` cycle (context:
  constraints below).
- ARCHIVED-STALE ⇒ stop and report to Ryan; never `--resume-delete`
  without his approval. `RECENT scr_p020r_geom_s020v_om15` (7 GB) is NOT
  approved — leave it, re-report it.
- Never touch `/home/rander39/campaigns/p021-fgs-stage1-20260922`
  (ARCHIVER_SKIP-marked; job 13858983 loads from it), nor pins/env of any
  campaign dir a queued/running job uses. 032 arms (`*_m13h`, A1 eng) may
  be RUNNING — LIVE/RECENT classification is the guard.
- Append the resulting ledger line where the storage agent indicates, and
  record before/after in a dated status file in BRAINSTORM/021.

## TASK B — babysit job 13858983, then harvest

1. Monitor via `hpc-monitor` only (never inline). Judge by outputs
   (STATUS_* files in the run dir), never sacct; task logs are
   output-buffered — judge liveness by CPU/outputs, not log mtime.
   Runtime ~ ladder 24 procs + bridge 3 + placement 6 + cap 6–12; 48 h
   wall. Individual stage FAILEDs are findings, not aborts (launcher
   continues and writes `COMPLETED` with a failed_count).
2. On completion: harvest via `harvester`; then run
   `benchmark/fgs_stage1_analysis.jl <run_dir>` — FIRST real execution,
   expect to debug it against the actual CSVs (it is only parse-checked).
3. **Gate**: if the 16→32 plateau or the 32→64 regression fails to
   reproduce, report deltas and STOP — no interventions on a
   non-reproduced effect. Plateau and regression are SEPARATE conclusions.
4. Resume path if the job dies: resubmit with
   `RESUME_FROM_JOB_ID=13858983` in the `--export` list (same
   COLD_PROJECT/CAMPAIGN_PINS/COLD_DATA_ROOT, from the deploy tree top
   level); STATUS_*=ok stages skip.

## Owed (carried)

- **Origin push** of branches+tags in all three repos once Ryan re-auths
  GitHub on the laptop (`gh auth login -h github.com`). Local-only
  commits now include the archiver work (`5c4f7d1`, `6388a92`) and the
  021 provenance/status commits through `7a8097b`.
- Archiver **T5 pre-existing failure** (exit 9 vs 8) awaits Ryan's ruling.
- fp-018gpu-p2rr-* (032 arms) left the queue ~2026-09-22 08:00 with no
  run dirs found — 032's thread owns this; flag only.
- Standing Ryan-gated ledger unchanged (see
  `fgs_scalability_reset_prompt_20260922.md` §"Standing Ryan-gated
  ledger").

## Traps (prior traps bind; key ones)

- Fixed-work rows: `solved=false`/`eligible=false` BY CONSTRUCTION — gate
  on certified-accepted + iterations==27 + 1e-8 repeat, never filter on
  eligible/solved. Accepted solves run ONE extra
  fmm+influence+residual vs fixed work.
- `diag_*` = −1 marks uninstrumented rows; never pool instrumented with
  uninstrumented. Worker cap does not pin threads to cores (caveat, not
  bug).
- `ssh orc` needs a live ControlMaster socket (`ssh orc -fN` if cold) AND
  `bash -lc "..."` for slurm/module commands. Local runs NEVER >4
  threads. Never edit source while a job uses its deployment — the
  rsync manifests make the deployed trees tamper-evident: any write under
  them fails the job at `cold_packages`.
- Run-dir names case-sensitive; judge runs by outputs, never sacct.

## House rules (binding)

Monitoring `hpc-monitor`; harvesting `harvester`; storage `hpc-storage`;
notebook writes Ryan-gated (offer, don't write); dated status/provenance
files in BRAINSTORM/021. HPC submission remains approved for STAGE 1 ONLY
(the 13858983 resume path included); Stage 2/3 Ryan-gated.
