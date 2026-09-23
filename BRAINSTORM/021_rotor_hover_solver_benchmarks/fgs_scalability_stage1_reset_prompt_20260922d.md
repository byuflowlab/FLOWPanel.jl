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

## TASK A — verify the two in-flight Ryan-approved archive workers

Storage history so far (full detail:
`archiver_campaigns_support_20260922.md`): the first worker (pid 500292)
COMPLETED CLEAN — `p018_csarc_n2_nt72_l3p0_3r_sfs3nb_om15` archived,
`TOTAL_FREED_MB=39612`, `VERIFY_FAIL_COUNT=0`, `STALE_COUNT=0`, home
342 G after. Ryan then APPROVED (2026-09-22) archiving the three RECENT
p032-rootomit runs AND `scr_p020r_geom_s020v_om15`; both were launched
detached and confirmed alive/progressing, **completion NOT yet verified**:

| pid | log | targets |
|---|---|---|
| 537577 | `/home/rander39/archiver_worker_p032rootomit_20260922_202851.log` | `p018_csarc_n2_nt72_l3p0_3r_{sfs3nb_om15_g25,srlx_g25_omi1,srlx_g25_omi1_mo4}` (155.2 GB, root campaigns/p032-rootomit-20260918/FLOWPanel.jl) |
| 539720 | `/home/rander39/archiver_worker_scrp020r_20260922_202919.log` | `scr_p020r_geom_s020v_om15` (~7 GB, root projects/FLOWPanel.jl) |

- Verify via `hpc-monitor`: clean finish per log = `ARCHIVE_MODE=APPLY`,
  `VERIFY_FAIL_COUNT=0`, `STALE_COUNT=0`, per-run `ARCHIVE ... freed=`
  lines (worker 1 must show THREE), `ARCHIVED.txt` breadcrumbs. Expected
  end state ≈ well under 200 G on home. Record the final freed figures
  and home headline in the status note; offer Ryan a ledger mirror for
  032.
- A lingering `.partial` tarball with a dead pid = interrupted transfer:
  nothing was deleted; relaunch the same `--only` apply (the approval
  stands). `VERIFY-FAIL` or `ARCHIVED-STALE` ⇒ stop, report to Ryan,
  never `--resume-delete` without his ruling.
- Never touch `/home/rander39/campaigns/p021-fgs-stage1-20260922`
  (ARCHIVER_SKIP-marked; job 13858983 loads from it), nor pins/env of any
  campaign dir a queued/running job uses.

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
