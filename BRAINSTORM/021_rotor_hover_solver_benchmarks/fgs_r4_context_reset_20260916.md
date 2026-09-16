# R4 follow-up: context-reset handoff (2026-09-16)

Prepared 2026-09-16 ~03:25 UTC. Supersedes `fgs_r4_context_reset_20260915b.md`
(its live state and its v16 deployment details are historical). Resume
execution and evidence analysis, not just planning.

## Keep the parent context small

Delegate bounded mechanical work to the repo subagents in `.claude/agents/`
(`hpc-monitor` read-only status/log checks, `harvester` tabulation,
`test-runner`, `brainstorm-scout`, `code-scout`), choosing the cheapest model
that is reliable; require ≤40-line summaries with an indexed artifact list.
Keep inline: physics/numerics reasoning, conclusions, code edits, job
submission, anything needing Ryan. Do not pour raw logs or full campaign
histories into the parent context. No prior agents or monitors survive the
reset — re-arm your own job watch.

## Required policy reads

`~/.claude/CLAUDE.md`, repo `CLAUDE.md` (+ `agent_policies/{WORKFLOW,TESTING,
HPC}.md` as routed). Before cluster work: current
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` (+ its slurm/storage
subdocs). ≥60 s between scheduler queries, longer backoff while waiting on a
long run. Local jobs ≤4 threads. `ssh orc` uses a live ControlMaster socket
with BatchMode — if auth fails, STOP (never retry into MFA); a sandbox socket
denial needs approved escalated SSH and is not an auth failure. Slurm
timestamps default to MDT; set `TZ=UTC` explicitly.

## Campaign state (verified 2026-09-16 03:20 UTC)

### DONE — v15 thread ladder (job 13694724): complete, harvested, audited

- All six arms j∈{1,4,8,16,32,64} PASS every gate (27 iterations invariant,
  identical instrumented/uninstrumented histories, BC rel-L2 4.78e-7, repeat
  delta 0, instrumentation overhead −2.1%…+0.9% inside batch spread).
- Local evidence (all 131 files SHA256-verified against remote):
  `BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_r4_followup_evidence_20260914/diag-v15-13694724/`
  with `analysis/{audit-full.txt, ladder_summary.md, analyze_ladder.py,
  analyze_j1.py, j1_summary.md, audit_r4_diag.py}`.
- **Measured conclusions are written as §7 of
  `fgs_opt_r4_diagnostics_package_20260912.md`** — read §7 before anything
  else. Headline: the old "nonself 77–80%" was a main-task-sample artifact;
  FMM far field dominates j1 wall (26.7 s/71%) but scales 24×; the serial
  lexicographic leaf-sweep chain (nonself products 7.89 s/72% at 1.04×,
  scatter 7.9% flat, leaf solve 5% flat) owns ~85% of j64 wall; whole-solve
  ceiling 3.43×. Census: 1,068 leaves, 2.86 GB matrices, 48,627 conflicts,
  79-color schedule (median size 16). **Test-5 condition SATISFIED: the
  colored sweep is the justified next experiment** (separately calibrated
  config, all gates, ranked by total time to accepted accuracy — iteration
  count may move off 27). Do NOT resubmit or rerun the v15 ladder.

### RUNNING/PENDING — v17 counters/activity job 13711596

At 03:20 UTC: PENDING (Resources), Slurm start estimate 05:30 UTC, walltime
cap 6 h; expected runtime ~1.5–3 h once started (precompile + smoke +
controls + j4/j64 perf-gated arms). Judge the run by its outputs, never by
sacct exit status (see memory: Slurm FAILED unreliable).

- Pins (submission provenance =
  `fgs_r4_followup_evidence_20260914/v17-deployment/submission-provenance-13711596.md`):
  FLOWPanel `campaign/p021-r4-counters-source-20260915-v17`
  (`b0b6eec150183c174f67fee80574c3ca95e77751`, worktree
  `/private/tmp/flowpanel-p021-r4-counters-v17`; = v16 `5460376` + bounded
  poll(2)+read(2) perf-FIFO acks with 15 s deadline, Julia perf smoke
  `test/r4_perf_control_smoke.jl`, real-FIFO driver tests — all passed
  locally, log `v17-deployment/counter-driver-local-v17.log`);
  FastMultipole `campaign/p021-r4-activity-source-20260915-v11`
  (`adb9967d`, unchanged from v16 prep); FLOWVPM
  `campaign/p021-cold-exec-20260910-v1` (`05c658f7`, remote worktree verified
  clean at pin pre-submit).
- Deployment (all hash-verified pre-submit, `sha256sum -c` PASS on 1,491 +
  6,441 files): silo generation
  `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo/counters-v17/`
  {FLOWPanel.jl, FastMultipole, env/{Project,Manifest}.toml, pins.toml,
  flowpanel.sha256, fastmultipole.sha256}; data symlink →
  `/home/rander39/projects/FLOWPanel.jl/data`; the partial `counters-v16/`
  was never used as a generation (only seeded the rsync).
- Output: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/counters-v17-13711596/`
  — expect `COMPLETED`, `campaign_pins.toml`, `Manifest.toml`, hardware/env
  captures, `perf-smoke.{csv,log}`, `controls-*.log` (counter driver,
  FastMultipole solver/coloring, FLOWPanel solver/history), and `j4-b1/`,
  `j64-b1/` each with `perf-stat.csv`, `process.log`, `numactl_show.txt`,
  `results/` (stage activity CSVs + status). Scheduler logs:
  `counters-v17/FLOWPanel.jl/logs/slurm/r4-counters-v17-13711596.{out,err}`.
- The launcher gates fixtures on the perf smoke and controls; perf events are
  disabled except around one warmed prepared solve per arm. **No performance
  comparisons from these diagnostic runs.**

## Next actions, in order

1. Watch 13711596 to a terminal state (background sacct poll ≥300 s spacing).
   Then harvest the whole run dir + both scheduler logs into
   `fgs_r4_followup_evidence_20260914/counters-v17-13711596/`, generate remote
   SHA256s, verify locally (pattern: v15 harvest in `diag-v15-13694724/`,
   `remote-sha256-full.txt` + `sha256sum -c`).
2. Verify gates from outputs: smoke PASS, all four controls PASS, both arms'
   status/results present, solver converged, certified FMM authoritative, BC
   rel-L2 ≤1e-6, repeat ≤1e-8, finite solutions, BLAS=1, zero resets.
   Diagnose any inconsistency before further work; retain failed evidence.
3. Analyze counters + stage activity (delegate tabulation; keep reasoning
   inline): reconcile with §7's wall budget. Hard limits to respect: generic
   cache counters CANNOT establish DRAM bandwidth saturation; the coarse
   `nearfield_update` activity span mixes leaf/product/scatter/updates; /proc
   tick resolution limits short stages. If saturation stays unresolved, record
   the limitation — the serial-execution evidence already justifies the
   colored-sweep experiment.
4. Update `fgs_opt_r4_diagnostics_package_20260912.md` §7 (replace its
   "pending 13711596" paragraph with measured findings + limitations) and
   append the live-state update to `fgs_r4_followup_validation_20260915.md`.
5. **Mandatory authorized cleanup** (Ryan pre-authorized; do not re-ask): once
   ALL jobs using any generation in
   `/home/rander39/FLOWPanel-p021-r4-diag-v10-silo` are terminal AND all
   results/scheduler logs/provenance are hash-verified locally, delete that
   silo directory entirely (includes counters-v16/ and counters-v17/), verify
   removal. NEVER follow its data symlinks or touch the canonical data root
   `/home/rander39/projects/FLOWPanel.jl/data`; preserve other silos,
   campaign worktrees, and shared deps. v14/v15 evidence is already local —
   retain it.
6. Stop after the evidence handoff. The colored-sweep experiment is the
   selected next experiment but is a NEW implementation/campaign phase: fresh
   v18+ pinned generation, separate calibration, all gates — do not start it
   inside this diagnostics pass without Ryan. **No notebook entry without
   Ryan's separate approval** (a ladder-findings entry is owed and drafted on
   request).

## Failure handling

- If 13711596 dies or a gate fails: harvest evidence first, diagnose root
  cause before any rerun. Any executable fix = new committed+annotated tag
  and a fresh silo generation (counters-v18/); never move existing tags,
  never edit v15/v17 deployed sources, never claim new contents under old
  pins. Config-only reruns may reuse the verified counters-v17 deployment.
- Local worktrees `/private/tmp/flowpanel-p021-r4-counters-v17` and
  `/private/tmp/fastmultipole-p021-r4-activity-v11` are clean at their tags;
  tags/commits are durable in the main repos (`git tag -l 'campaign/p021-*'`).
- The live repo has unrelated modified/untracked files from other campaigns —
  do not revert, clean, or commit them.

## Compact document index

All under `BRAINSTORM/021_rotor_hover_solver_benchmarks/` unless noted.

1. `fgs_opt_r4_diagnostics_package_20260912.md` — the deliverable; §7 = this
   phase's measured update (read §7 + §6 gates; §1–§5 as needed).
2. `fgs_opt_r4_diagnostics_handoff_20260912.md` — governing **2026-09-14
   follow-up request** (top section only; tests 1–4 complete, test 5
   condition satisfied, counters = last piece of test 4).
3. `fgs_r4_followup_validation_20260915.md` — chronological validation log;
   tail = 2026-09-16 03:15 update.
4. `fgs_r4_followup_evidence_20260914/` — all local evidence:
   `diag-v15-13694724/` (ladder + analysis), `v17-deployment/`
   (submission provenance, manifests, prepare_deployment.py, local test
   logs), older diag-v10/fm-v14 dirs.
5. Memory: `project_021_solver_benchmarks.md` (updated 2026-09-16).
