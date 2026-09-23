# Reset prompt: 021 FGS scalability Stage 1 — rsync deploy + submit (2026-09-22, supersedes fgs_scalability_stage1_reset_prompt_20260922b.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Where this stands

BRAINSTORM 021 is CLOSED AT PROMOTION; active follow-on = FGS scalability
diagnostic, plan C (`fgs_scalability_diagnostic_plan_20260921c.md`). Stage 0
complete (`fgs_scalability_stage0_audit_20260922.md`). Read `CLAUDE.md`
first and `agent_policies/HPC.md` before HPC work. Provenance so far:
`fgs_scalability_stage1_provenance_20260922.md` (pins recorded; deployment
checklist open).

**Stage-1 code is SMOKED, COMMITTED, and TAGGED (2026-09-22 session):**

- Annotated tag `campaign/p021-fgs-stage1-20260922` on the triple:
  - FLOWPanel.jl branch `fastmultipole`, tag at `c8621af` (branch tip is
    `5f0e53b` after provenance + R2-j1 harvest commits — deploy the TAG)
  - FastMultipole branch `flowpanel-20260817`, tag at `649405cb`
  - FLOWVPM.jl branch `flowpanel`, tag at `8d4a3b4`
- Smokes all green: `runtests_unit_solver.jl` full pass; 4-thread dagteam
  worker-cap equivalence 37/37 (caps 0/1/2/4/99 → team size
  `clamp(cap,1,nthreads)`, bit-identical solutions across caps, diagnostics
  keys subsets of `:nonself_product_ns`, zero for other sweep orders);
  `cold_check_config` whitelist smoke pass. Analysis script
  `benchmark/fgs_stage1_analysis.jl` committed (parse-checked only).
- Owed R2-j1 CLOSED (`fgs_scalability_r2j1_harvest_20260922.md`): all
  13829232 arms harvested. 13829231_0 (R4 thread-scaling) may still be
  RUNNING — its FINAL table is already harvested, nothing owed from it;
  just don't disturb it.
- **Origin push is BLOCKED**: laptop GitHub creds dead (gh token invalid,
  no ssh key, keychain empty). Ryan has ruled (below); push to origin
  remains owed later, after he re-auths (`gh auth login -h github.com`).

## Ryan's rulings (2026-09-22, verbatim intent)

1. **APPROVED: teach the storage agent (archiver) about
   `/home/rander39/campaigns`** (209.2 GiB, 57% of home usage, currently
   invisible to `run_archiver.sh --all-checkouts` because the top dir has
   no Project.toml).
2. **APPROVED: use rsync to update the HPC version of the code** (the
   committed Stage-1 code) in lieu of the blocked origin push.
3. **Then launch the jobs** — Stage-1 HPC submission (already pre-approved
   for Stage 1 ONLY; Stage 2/3 remain Ryan-gated).

## TASK 1 — rsync-deploy the pinned triple to orc

The harness has a sanctioned rsync mode: `cold_packages` in
`benchmark/fgs_cold_common.jl` (~lines 277–296) accepts pins with
`deployment = "rsync"` + `content_manifest` + `content_manifest_sha256`,
verifies the manifest hash and runs `sha256sum --quiet -c` over the
deployed tree at job start, and skips the git-worktree checks. Use it.

- Deploy dir (fresh, following the established layout):
  `/home/rander39/campaigns/p021-fgs-stage1-20260922/{FLOWPanel.jl,FastMultipole,FLOWVPM.jl}`
  plus `env/` and `pins.toml`. Schema example:
  `/home/rander39/campaigns/p021-r12-champion-20260919/pins.toml`
  (that one is git_worktree mode; add the three rsync keys per package).
- **NEVER rsync onto an existing git checkout or old deployment** — the
  2026-08-24 provenance collision (HPC.md) came from exactly that. Fresh
  directories only.
- Ship exactly the tagged tracked content, not the dirty working trees
  (all three live checkouts carry other campaigns' uncommitted edits).
  Cleanest: per repo
  `git archive campaign/p021-fgs-stage1-20260922 | ssh orc 'mkdir -p <dir> && tar -x -C <dir>'`
  (this IS the rsync-style content push Ryan approved; plain rsync from a
  clean `git worktree add --detach <tmp> <tag>` export is equally fine).
- Content manifest per repo: sha256sum-format file over ALL deployed files
  (relative paths, generated from the same export), placed on orc at an
  absolute path (e.g. `<deploy>/MANIFEST.<name>.sha256`); record each
  manifest's own sha256 in pins.toml. Runtime enforces both.
- pins.toml per package: `path` (the orc deploy dir — must realpath-match
  what Julia loads), `tag = "campaign/p021-fgs-stage1-20260922"`, `sha`
  (from the table above), `deployment = "rsync"`, `content_manifest`,
  `content_manifest_sha256`.
- `env/`: copy the Project.toml pattern from
  `/home/rander39/campaigns/p021-r12-champion-20260919/env`, then
  `Pkg.develop(path=...)` the three deployed trees +
  `Pkg.instantiate()`; confirm `pkgdir(FLOWPanel)` etc. resolve to the
  deploy dirs. This env is COLD_PROJECT.
- The R4 mesh (`examples/data/dji9443_20260813_*_capped_captess4.msh`) is
  git-tracked, so `git archive` ships it — verify it exists in the deployed
  tree anyway (HPC.md inputs rule).
- Update the deployment checklist in
  `fgs_scalability_stage1_provenance_20260922.md` (note Ryan's rsync
  approval + the manifest sha256s + deploy paths) and commit.

## TASK 2 — submit Stage 1

1. `ssh orc -fN` if the ControlMaster socket is cold. NOTE: non-login ssh
   lacks slurm bins — wrap remote commands in `bash -lc "..."`.
2. slurm-availability skill: `--cpus 128 --mem-gb 500 --eta` (match the
   job's REAL ask: 128 cpus, exclusive zen3). Confirm an m12/zen3 node is
   reachable under `--qos=normal`.
3. From the FLOWPanel **deploy dir** top level: `mkdir -p logs/slurm`, then
   ```
   sbatch --export=ALL,COLD_PROJECT=/home/rander39/campaigns/p021-fgs-stage1-20260922/env,CAMPAIGN_PINS=/home/rander39/campaigns/p021-fgs-stage1-20260922/pins.toml,COLD_DATA_ROOT=/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910 \
     benchmark/run_r4_fgs_stage1.slurm.sh
   ```
   (CALIB_JOB defaults to 13777133; the launcher pre-checks the four
   calibrated `dagteam_selected.toml` files — verified present 2026-09-22.)
4. Record the job id in the provenance file, commit, then monitor via
   `hpc-monitor` only. Runtime ~ ladder 24 procs + bridge 3 + placement 6 +
   cap 6–12, 48 h wall. Failures are findings (STATUS_* files); judge by
   outputs, never sacct.
5. When it completes: harvest via `harvester`, run
   `benchmark/fgs_stage1_analysis.jl <run_dir>` (first real execution —
   expect to debug it against the actual CSVs). Gate: if plateau or
   regression fails to reproduce, report deltas and STOP (no
   interventions on a non-reproduced effect).

## TASK 3 — archiver support for /home/rander39/campaigns (approved)

- First SURVEY: `du -sh /home/rander39/campaigns/*/` and inside the big
  ones — the 209 G is presumably VTK under the deployed trees' own `data/`
  dirs (old runs that wrote relative paths). Know the shape before coding.
- The deployed `FLOWPanel.jl` trees inside campaigns DO have Project.toml,
  so `--root <campaign>/FLOWPanel.jl` may already work per-tree; the gap is
  `--all-checkouts` discovery. Teach it to scan
  `/home/rander39/campaigns/*/` (keep the Project.toml+data/ gate per
  tree; dedupe by realpath — some may symlink to the projects data root).
- **Mandatory after any archiver edit**: run
  `scripts/tests/run_archiver_test.sh` locally, mirror the script to
  `/home/rander39/projects/FLOWPanel.jl/scripts/` with md5sum verify on
  both sides (HPC.md rule).
- EXCLUDE the new `p021-fgs-stage1-20260922` deploy dir from any
  archiving (it's code + env, not runs) and never touch pins/env of any
  campaign dir a queued/running job uses.
- Then launch `hpc-storage` for an apply cycle. Also still awaiting Ryan
  (NOT approved yet): the RECENT item
  `scr_p020r_geom_s020v_om15` (7 GB, quiet 23 h) — leave unless he rules.

## Owed (carried)

- Origin push of branches+tags in all three repos once Ryan re-auths
  GitHub on the laptop.
- fp-018gpu-p2rr-* (032 arms) left the queue ~2026-09-22 08:00 with no
  run dirs found by the storage scan — 032's thread owns this; flag only.
- Standing Ryan-gated ledger unchanged (see
  `fgs_scalability_reset_prompt_20260922.md` §"Standing Ryan-gated ledger").

## Traps (prior traps bind; key ones + new)

- Fixed-work rows: `solved=false`/`eligible=false` BY CONSTRUCTION — gate
  on accepted+iterations==27+repeat, never filter on eligible/solved.
- Accepted solves run ONE extra fmm+influence+residual vs fixed work.
- `diag_*` = −1 marks uninstrumented rows; never pool instrumented with
  uninstrumented. Worker cap does not pin threads to cores (caveat, not
  bug). Plateau and regression = SEPARATE conclusions.
- Run-dir names case-sensitive (`r12-champion-R2-j1-...`); task logs
  output-buffered (judge liveness by CPU/outputs, never log mtime); judge
  runs by outputs, never sacct; local runs NEVER >4 threads; never edit
  source while a job uses its deployment; `ssh orc` needs a live
  ControlMaster socket AND `bash -lc` for slurm commands.
- rsync mode: pin `path` must realpath-match the loaded pkgdir; manifest
  must cover exactly the deployed files (regenerate if you re-deploy).

## House rules (binding)

Monitoring `hpc-monitor`; harvesting `harvester`; storage `hpc-storage`;
notebook writes Ryan-gated (offer, don't write); dated status/provenance
files in BRAINSTORM/021. HPC submission approved FOR STAGE 1 ONLY.
