# R4 follow-up: context-reset handoff

Prepared at Ryan's request; last updated **2026-09-16 02:19 UTC** (2026-09-15 local MDT).
Resume execution and evidence analysis, not just planning. This supersedes the
older clear-context prompt's live-state/deployment details.

## Keep the parent context small

Ryan explicitly requests cheaper subagents for mechanical tasks and **Astra
subagents for tasks requiring substantial reasoning**. Delegate bounded,
independent work; require compact summaries with an indexed artifact list.
Do not pour raw logs, full campaign histories, or source manifests into the
parent context. Read summaries first and inspect only relevant evidence.

- `gpt-5.6-luna`: status/log checks, file inventory, simple harvesting and
  tabulation. Parent must check critical conclusions; this model needed
  corrections when implementing a schema-sensitive audit script.
- `gpt-5.6-terra`: Julia verification or nontrivial mechanical code/test tasks.
- `gpt-6-astra`: bounded reasoning-heavy reviews: counter measurement scope,
  interpretation of stage/scaling evidence, numerical acceptance, and choosing
  a justified next experiment. Give each a narrow question and a file index.
- Prefer `fork_turns="none"` with explicit context and ≤40–60-line reports.
  Model overrides require a fresh/limited fork. Do useful parent work alongside
  agents; do not spawn agents merely to wait.
- Observe repo roles `.claude/agents/{hpc-monitor,harvester,test-runner,
  brainstorm-scout}.md`. Never assume old agents survive the reset.
- The last Luna monitor turn hit an account/model usage limit. Do not repeatedly
  retry an exhausted model. If still unavailable, report this and do a small
  bounded status check inline or use an available appropriate model.

## Immediate priority: harvest the completed ladder

**Fresh verification at 2026-09-16 02:19:05 UTC:** Ryan re-established SSH;
the existing connection now works. Job **13694724 COMPLETED, exit 0:0**, elapsed
**08:52:49**, ending **2026-09-15 23:10:21 UTC**. The root `COMPLETED` file and
all six `j1/4/8/16/32/64-b1/results/status.toml` files exist remotely. The j64
log reaches the 58,192-panel case (direct source assembly 34.4 s), without a
reported error in its tail. Full j4+ numerical/provenance auditing is still
pending: scheduler completion and marker existence alone are insufficient.

Start by harvesting all remaining arm results, root completion and scheduler
logs, verify hashes, then run the prepared local audit. Do not resubmit v15.
The earlier 19:59 UTC RUNNING/j4 observation below is historical only.

Before Ryan reopened SSH, a deployment transfer returned
`Permission denied (keyboard-interactive)`; retries were stopped. That access
block is now resolved, but the v16 deployment remains partial exactly as
recorded below. No new job was submitted this turn.

Use the existing ControlMaster with BatchMode and bounded connection timeout.
Sandbox socket denial previously required approved `require_escalated` SSH;
do not confuse it with authentication failure. If the socket later expires,
stop rather than retry into MFA. Slurm timestamps default to MDT; the fresh
accounting command explicitly set `TZ=UTC`.

## Required policy reads and compact document index

Read `/Users/ryan/.claude/CLAUDE.md`, repository `AGENTS.md` and `CLAUDE.md`, and
`agent_policies/{WORKFLOW,TESTING,HPC}.md`. Before cluster work, read current
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` and applicable `ai-docs/`
references (slurm, storage, performance; development/software for deployment).
At least 60 seconds between periodic scheduler queries; use longer backoff for
long runs. Local computations use at most four threads.

Have a scout read the governing sections, rather than loading the whole history:

1. `fgs_r4_followup_validation_20260915.md`: v15 fix/provenance plus this turn's
   appended j1 analysis; its 19:59 snapshot is historical.
2. `fgs_opt_r4_diagnostics_handoff_20260912.md`: **2026-09-14 follow-up request**
   defines required controls, diagnostics, deliverables, and conditional next
   experiment. Older assignment sections are historical.
3. `fgs_opt_r4_diagnostics_package_20260912.md`: scientific interpretation and
   retained configuration; needs final measured update.
4. `fgs_r4_followup_validation_resume_20260914.md`: older deployment and cleanup
   background, only consult targeted sections when needed.

## Existing v15 campaign

- Silo: `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo`
- Output: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/diag-v15-13694724/`
- Scheduler logs: silo `FLOWPanel.jl/logs/slurm/r4-diag-v15-13694724.out/.err`
- Environment/pins: silo `env/` and `pins.toml`
- FLOWPanel: `campaign/p021-cold-source-20260915-v15`
  (`39ec4e3630bc6c04f0865a1d3feecce7130904fe`), worktree
  `/private/tmp/flowpanel-p021-r4-diag-v15`
- FastMultipole: `campaign/p021-r4-diag-source-20260914-v10`
  (`87cbc8460b51f24ddf34cc5f41a1d1b6682bf04a`), worktree
  `/private/tmp/fastmultipole-p021-r4-diag-v10`
- FLOWVPM: `campaign/p021-cold-exec-20260910-v1`
  (`05c658f7804ec5f9b68d4cb9826a9f97cfecb373`), remote worktree
  `/home/rander39/campaigns/p021-cold-20260910-v1/FLOWVPM.jl`

The original v15 source, environment, pins and data symlink were NOT modified.
Its requested resources are exclusive Zen3, 64 CPUs, 500 GiB, 12 hours; BLAS=1.
It runs j1/4/8/16/32/64 serially in one allocation, 40 trials per arm in batches
uninstrumented/instrumented/uninstrumented/instrumented (10 each).

## Verified local evidence and partial conclusions

All relative paths here are under `BRAINSTORM/021_rotor_hover_solver_benchmarks/`.

`fgs_r4_followup_evidence_20260914/diag-v15-13694724/` contains 24 root files
and 19 files from completed `j1-b1/`. **All 43 match remote SHA256 digests**;
see `remote-sha256.txt`, `sha256-verification.txt`, and `harvest_summary.md`.
No j4+ results were harvested this turn. Root `COMPLETED` is now confirmed
remotely but has not been harvested; the saved local audit is a partial snapshot.

| j1 result | Value |
|---|---:|
| Uninstrumented median (20 trials) | 37.5594 s |
| Uninstrumented range | 37.3886–37.6857 s |
| Instrumented median (20 trials) | 37.4689 s |
| Instrumented range | 37.0734–38.2301 s |
| Median difference | -0.0905 s (-0.241%) |
| Certified FMM BC rel-L2 | 4.77952513e-7 |
| Direct BC rel-L2, equivalence controls | 4.78150591e-7 |
| Direct/FMM disagreement | 5.4676859e-9 |
| Repeat/instrumentation solution difference | 0 |
| Iterations / inner sweeps | 27 / 81 |

All j1 timed gates and profile validation pass; instrumentation histories are
identical. Negative median overhead is not a speedup claim; ranges overlap and
independent process replication is absent. Instrumented medians: FMM 26.73485 s,
nonself products 8.16198 s, initialization 1.00337 s, scatter 0.79205 s, leaf
solve 0.59771 s. Internal total 37.46101 s, exclusive sum 37.44830 s,
unaccounted 0.01268 s; outer solver adds median 0.00788 s. Per-trial additive
checks pass within CSV decimal rounding (maximum discrepancy <1e-7 s).

Census: 1,068 leaves, 2,862,850,032 Float64 matrix bytes, 95,390 directed
edges, 48,627 undirected conflicts. Independent edge/degree checks pass.
Hypothetical coloring: 79 groups, sizes 1–27 (median 16), no same-color
conflicts. Execution remains **serial lexicographic** in the leaf loop.

Reusable local analysis is saved under that evidence directory's `analysis/`:
- `audit_r4_diag.py`: numerical/completeness audit; accepts run-directory arg.
  Exit 1 is expected while arms/root completion are missing. This is not
  provenance verification. It was corrected against real j1 schema and tested
  with injected numerical/stage failures. It accounts for CSV rounding.
- `audit-partial.txt`: j1 PASS, campaign incomplete.
- `analyze_j1.py`, `j1_summary.md`: census, stage and profile details.

Profile fractions are not wall fractions. Existing thread CPU ticks span reset
and validation too, not just solves. NaN direct metrics mean unevaluated.
Generic perf PMU probes succeeded on the compute node; no workload-scoped
hardware counters have yet been measured, so bandwidth saturation is unresolved.

## Prepared v16 counter/activity generation — NOT submitted

Two clean local worktrees were committed and annotated:

| Package | Worktree | Tag | SHA |
|---|---|---|---|
| FLOWPanel | `/private/tmp/flowpanel-p021-r4-counters-v16` | `campaign/p021-r4-counters-source-20260915-v16` | `546037672ac7e41eba5c64cedf6a2e4d6b7394e6` |
| FastMultipole | `/private/tmp/fastmultipole-p021-r4-activity-v11` | `campaign/p021-r4-activity-source-20260915-v11` | `adb9967d5b696cf9ab05c557aa2a649e100004dc` |

FLOWVPM pin is unchanged. Files to review (use a bounded Astra reviewer for
measurement reasoning):
- FLOWPanel `benchmark/fgs_r4_counters.jl`: separate baseline, coarse CPU-stage
  observation, then perf-gated prepared solve; zero resets, certified FMM,
  direct baseline check, finite solutions, exact history/repeat agreement,
  memory gates. No performance comparisons from these diagnostic runs.
- `benchmark/run_r4_counters.slurm.sh`: exclusive Zen3/64 CPU/500G/6h request,
  j4 and j64, perf FIFO smoke before fixtures, solver/history/FM controls.
  Perf starts disabled and is enabled only around one warmed solve.
- `test/runtests_r4_counters_driver.jl`: actual AST/world-age/helper tests.
- `test/r4_perf_control_smoke.py`: installed perf FIFO smoke, 15s ack timeout.
- FLOWPanel `src/FLOWPanel_solver.jl`: forwards optional stage observer.
- FastMultipole `src/solve.jl`, `test/solve_test.jl`: coarse paired observations
  outside leaf hot loops, unchanged ordering, baseline/observer equivalence.

Local driver AST/protocol/bash checks and FastMultipole solver/coloring tests
passed at j1 and j4, BLAS=1. Logs are copied to
`fgs_r4_followup_evidence_20260914/v16-deployment/`. The isolated test environment
is `/private/tmp/p021-r4-counters-verify-v1`. Linux perf smoke, actual Julia FIFO
control on Linux, and full FLOWPanel cluster controls have **not** run. The
launcher gates expensive fixture work on its smoke/controls, but review the
measurement boundaries and protocol robustness before submission. Generic
cache counters cannot establish DRAM bandwidth saturation. The coarse
`nearfield_update` activity span combines leaf/product/scatter and updates;
/proc tick resolution and sequential snapshots limit short-stage conclusions.

### Partial deployment: exact state

New generation directory, separate from every v15 loaded source path:
`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/counters-v16/`.

- Created `FLOWPanel.jl/`, `FastMultipole/`, `env/` under this directory.
- FLOWPanel rsync completed: 1,490 tracked regular files, 511,597,474 bytes.
- FastMultipole rsync then failed at SSH authentication; **not verified or
  established as transferred**. Do not treat this generation as deployable.
- Environment, pins, source manifests, data symlink, and logs directory have
  **not** been installed remotely. No remote full-source verification yet.
- No running job refers to these new paths; v15 paths remain untouched.

Local deployment material is ready in `.../v16-deployment/`:
`flowpanel.files0`, `fastmultipole.files0` (NUL-separated rsync file lists),
`flowpanel.sha256`, `fastmultipole.sha256`, `Project.toml`, `Manifest.toml`,
`pins.toml`, and local test logs. Source manifests exclude canonical `data/`
and symlinks, matching prior campaign scope.

| Manifest | Files | Bytes | SHA256 of manifest |
|---|---:|---:|---|
| flowpanel.sha256 | 1490 | 511597474 | `1bd61e1fc87281ba2c4e9d30429d22eedaa32fb3158a61adee189f1b33032ee7` |
| fastmultipole.sha256 | 6441 | 123522579 | `93ac52068a66db5aad376e4e03556b50f6626e4fab901d7f1202c3570413341b` |

The prepared Manifest points FLOWPanel and FastMultipole at these new nested
paths, and FLOWVPM at its unchanged campaign worktree. Finish rsync using
`--from0 --files-from=<package>.files0`, without `--delete`; copy env files to
`counters-v16/env/`, pins/manifests to `counters-v16/`. Create only the new
FLOWPanel generation's `data` symlink to the canonical root and its
`logs/slurm/` directory. Verify every remote source hash, manifest digest,
loaded paths, VPM clean tag/SHA, assets, storage headroom and allocation before
submission. Record evidence first. No need to redeploy or modify v15.

If v16 needs any executable correction, make a new recorded generation/tag;
do not silently change annotated pins or claim new contents under old pins.

## Next actions

1. SSH is restored and 13694724 is confirmed COMPLETED. Inspect and harvest
   every arm; do not resubmit or rerun the ladder. Diagnose any artifact or
   numerical inconsistency before deciding on further work.
2. Harvest all completed results and logs; verify hashes and provenance.
   Require root COMPLETED plus all six per-arm statuses and all numerical gates.
3. Review and finish the separate counter generation (or correct with fresh
   pins), refresh resource availability/storage, run its smoke and controls,
   and collect workload counters plus stage activity. Do not modify live source.
4. Delegate bulk tabulation; use bounded Astra reasoning to reconcile wall-time
   budgets, overhead, affinity/NUMA, block/dependency census and scaling. Follow
   the governing request for the conditional colored-sweep experiment or a
   concrete justified next experiment. Do not infer memory saturation or parallel
   GEMVs from configured thread count/profile percentages alone.
5. Update the diagnostics package/provenance with measured conclusions and
   limitations. Preserve convergence, finite solutions, certified authoritative
   FMM, BC rel-L2 <=1e-6, repeat <=1e-8, direct/FMM <=1e-7 when evaluated,
   BLAS=1, zero-reset and constructor-free timing. Keep profiling separate.
6. Complete the already-authorized bounded cleanup below. Stop after evidence
   handoff; no notebook entry without Ryan's separate approval.

## Authorization and mandatory cleanup

Ryan authorized these tests and explicitly allowed rsync into a new silo
**so long as it is deleted when done**. Do not ask for redundant deployment or
cleanup authorization. This overrides ordinary Git-only/no-new-silo rules for
this campaign; local pins remain clean annotated-tag worktrees.

After **all jobs using any generation in the silo are terminal**, and **all
results, scheduler logs, source/environment provenance, and failed-run evidence
are hash-verified outside it**, delete only:

`/home/rander39/FLOWPanel-p021-r4-diag-v10-silo`

This includes the nested unfinished/new counter generation. Verify removal.
Never follow its data symlinks or delete the canonical data root. Preserve
shared dependencies and other silos. Prior failed-run evidence and v14 passed
control evidence are already local and must be retained. Cleanup conditions
are currently unmet. No notebook entry was written.

The live repository has unrelated modified/untracked files, including other
campaigns that changed during this session. Do not revert, clean, or commit
those. Do not rely on earlier agents/monitors being alive.
