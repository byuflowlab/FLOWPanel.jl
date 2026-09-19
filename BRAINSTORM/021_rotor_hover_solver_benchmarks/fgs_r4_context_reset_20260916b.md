# R4 follow-up: context-reset handoff (2026-09-16b)

> **SUPERSEDED 2026-09-17** by `fgs_r4_context_reset_20260917.md` — the
> diagnostics pass completed (job 13733332, after v18/v19 failures), §7 is
> updated, the silo is deleted, and the campaign lineage is merged back.
> The 20260917 doc directs the colored-sweep experiment.

Prepared 2026-09-16 ~04:40 UTC. Supersedes `fgs_r4_context_reset_20260916.md`
(its campaign-state and next-actions sections are stale: its v17 job 13711596
FAILED and was replaced by v18 job 13712587). Resume execution and evidence
analysis, not just planning.

## Keep the parent context small

Delegate bounded mechanical work to the repo subagents in `.claude/agents/`
(`hpc-monitor` read-only status/log checks, `harvester` tabulation,
`test-runner`, `brainstorm-scout`, `code-scout`), choosing the cheapest model
that is reliable; require ≤40-line summaries with an indexed artifact list.
Keep inline: physics/numerics reasoning, conclusions, code edits, job
submission, anything needing Ryan. Do not pour raw logs or full campaign
histories into the parent context. No prior agents, monitors, or cron watches
survive the reset — re-arm your own job watch.

## Required policy reads

`~/.claude/CLAUDE.md`, repo `CLAUDE.md` (+ `agent_policies/{WORKFLOW,TESTING,
HPC}.md` as routed). Before cluster work: current
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` (+ its slurm/storage
subdocs). ≥60 s between scheduler queries, longer backoff while waiting on a
long run. Local jobs ≤4 threads. `ssh orc` uses a live ControlMaster socket
with BatchMode — if auth fails, STOP (never retry into MFA); a sandbox socket
denial needs approved escalated SSH and is not an auth failure. Slurm
timestamps default to MDT; set `TZ=UTC` explicitly. Note: `sbatch/sacct` need
`source /etc/profile` in non-interactive ssh.

## Campaign state (verified 2026-09-16 04:30 UTC)

### DONE — v15 thread ladder (job 13694724): complete, harvested, audited

Unchanged from the prior handoff. All six arms j∈{1,4,8,16,32,64} PASS every
gate; evidence SHA256-verified locally at
`fgs_r4_followup_evidence_20260914/diag-v15-13694724/`. **Measured conclusions
are §7 of `fgs_opt_r4_diagnostics_package_20260912.md`** — read §7 first.
Headline: serial lexicographic leaf-sweep chain owns ~85% of j64 wall; FMM far
field dominates j1 (71%) but scales 24×; whole-solve ceiling 3.43×; Test-5
condition SATISFIED — the colored sweep is the justified next experiment. Do
NOT rerun the v15 ladder.

### DEAD — v17 counters job 13711596: FAILED at smoke gate, postmortem closed

40-s failure (exit 1): `test/r4_perf_control_smoke.jl` had a triple-quoted
string directly above `using Test`; Julia parses it as a docstring for the
`using` statement and errors. No arms ran. Root cause of the testing gap: the
local v17 gate ran only `runtests_r4_counters_driver.jl`; the smoke script's
sole caller is the Slurm launcher, so it was never executed pre-submit. Failed
evidence harvested + SHA256-verified:
`fgs_r4_followup_evidence_20260914/counters-v17-13711596-FAILED/`. Do not
rerun v17 or touch its deployed sources.

### RUNNING/PENDING — v18 counters/activity job 13712587

Submitted 2026-09-16 ~04:20 UTC; walltime cap 6 h; expected ~1.5–3 h once
started (precompile + smoke + controls + j4/j64 perf-gated arms). Judge the
run by outputs, never sacct exit status.

- v18 = v17 + two changes (commit `6033dbdc06df3519e6c623245bcd747dad064555`,
  tag `campaign/p021-r4-counters-source-20260915-v18`): smoke docstring →
  comment; launcher v17→v18 naming. This time the smoke script WAS executed
  locally against a mock perf-FIFO responder (exit 0, correct
  disable/enable/disable ack sequence) and the driver suite re-passed; sibling
  scripts scanned, no repeats of the pattern.
- Other pins unchanged: FastMultipole
  `campaign/p021-r4-activity-source-20260915-v11` (`adb9967d`); FLOWVPM
  `campaign/p021-cold-exec-20260910-v1` (`05c658f7`, remote worktree verified
  clean at pin pre-submit). Local worktrees:
  `/private/tmp/flowpanel-p021-r4-counters-v17` (directory keeps its v17 name,
  now checked out clean at the v18 tag) and
  `/private/tmp/fastmultipole-p021-r4-activity-v11`.
- Deployment: fresh generation
  `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/counters-v18/` (cp -al seed
  from verified v17 + delta rsync of git-tracked lists); remote `sha256sum -c`
  PASS (FLOWPanel 1,491 files, FastMultipole 6,441); Manifest dev-paths
  verified at v18 paths; data symlink → canonical data root; storage preflight
  365 G/400 G. Full provenance:
  `fgs_r4_followup_evidence_20260914/v18-deployment/submission-provenance-13712587.md`.
- Output: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/counters-v18-13712587/`
  — expect `COMPLETED`, `campaign_pins.toml`, `Manifest.toml`, hardware/env
  captures, `perf-smoke.{csv,log}`, `controls-*.log`, and `j4-b1/`, `j64-b1/`
  each with `perf-stat.csv`, `process.log`, `numactl_show.txt`, `results/`.
  Scheduler logs:
  `counters-v18/FLOWPanel.jl/logs/slurm/r4-counters-v18-13712587.{out,err}`.
- Perf events are disabled except around one warmed prepared solve per arm.
  **No performance comparisons from these diagnostic runs.**

## Next actions, in order

1. Watch 13712587 to a terminal state (background sacct poll ≥300 s spacing).
   Then harvest the whole run dir + both scheduler logs into
   `fgs_r4_followup_evidence_20260914/counters-v18-13712587/`, generate remote
   SHA256s, verify locally (`remote-sha256-full.txt` + `sha256sum -c`; pattern:
   the v15 and v17-FAILED harvests).
2. Verify gates from outputs: smoke PASS, all four controls PASS, both arms'
   status/results present, solver converged, certified FMM authoritative, BC
   rel-L2 ≤1e-6, repeat ≤1e-8, finite solutions, BLAS=1, zero resets.
   Diagnose any inconsistency before further work; retain failed evidence.
3. Analyze counters + stage activity (delegate tabulation; keep reasoning
   inline): reconcile with §7's wall budget. Hard limits: generic cache
   counters CANNOT establish DRAM bandwidth saturation; the coarse
   `nearfield_update` activity span mixes leaf/product/scatter/updates; /proc
   tick resolution limits short stages. If saturation stays unresolved, record
   the limitation — the serial-execution evidence already justifies the
   colored-sweep experiment.
4. Update `fgs_opt_r4_diagnostics_package_20260912.md` §7 (replace its
   "pending 13711596" paragraph with measured findings + limitations, noting
   the v17 failure and that data came from v18/13712587) and append the
   live-state update to `fgs_r4_followup_validation_20260915.md` (its tail =
   2026-09-16 04:25 entry on the failure/resubmission).
5. **Mandatory authorized cleanup** (Ryan pre-authorized; do not re-ask): once
   ALL jobs using any generation in
   `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo` are terminal AND all
   results/scheduler logs/provenance are hash-verified locally, delete that
   silo directory entirely (counters-v16/, -v17/, -v18/), verify removal.
   NEVER follow its data symlinks or touch the canonical data root
   `/home/rander39/projects/FLOWPanel.jl/data`; preserve other silos, campaign
   worktrees, and shared deps. v15/v17-FAILED evidence is already local —
   retain it.
6. Stop after the evidence handoff. The colored-sweep experiment is the
   selected next experiment but is a NEW implementation/campaign phase: fresh
   v19+ pinned generation, separate calibration, all gates — do not start it
   inside this diagnostics pass without Ryan. **No notebook entry without
   Ryan's separate approval** (a ladder-findings entry is owed and drafted on
   request).

## Failure handling

- If 13712587 dies or a gate fails: harvest evidence first, diagnose root
  cause before any rerun (pattern: the 13711596 postmortem in the v18
  provenance file). Any executable fix = new committed+annotated tag and a
  fresh silo generation (counters-v19/); never move deployed tags, never edit
  deployed sources, never claim new contents under old pins. Config-only
  reruns may reuse the verified counters-v18 deployment.
- Tags/commits are durable in the main repos (`git tag -l 'campaign/p021-*'`).
- The live repo has unrelated modified/untracked files from other campaigns —
  do not revert, clean, or commit them. The v18 fix commit lives only on the
  worktree branch `campaign/p021-r4-counters-v17-wt`; merge-back to the
  campaign line is future work, not this pass.

## Compact document index

All under `BRAINSTORM/021_rotor_hover_solver_benchmarks/` unless noted.

1. `fgs_opt_r4_diagnostics_package_20260912.md` — the deliverable; §7 = this
   phase's measured update (read §7 + §6 gates; §1–§5 as needed).
2. `fgs_opt_r4_diagnostics_handoff_20260912.md` — governing 2026-09-14
   follow-up request (top section only; tests 1–4 complete, test 5 condition
   satisfied, counters = last piece of test 4).
3. `fgs_r4_followup_validation_20260915.md` — chronological validation log;
   tail = 2026-09-16 04:25 failure/resubmission entry.
4. `fgs_r4_followup_evidence_20260914/` — all local evidence:
   `diag-v15-13694724/` (ladder + analysis), `counters-v17-13711596-FAILED/`
   (failed-run harvest), `v17-deployment/` and `v18-deployment/` (provenance,
   manifests, prepare_deployment.py), older diag-v10/fm-v14 dirs.
5. Memory: `project_021_solver_benchmarks.md`.
