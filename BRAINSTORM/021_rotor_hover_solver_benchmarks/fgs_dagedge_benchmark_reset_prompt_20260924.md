# Reset prompt: 021 :dagedge HPC benchmark — performance run, then profile run (2026-09-24, supersedes fgs_edgepull_reset_prompt_20260924.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Mission (Ryan 2026-09-24)

Launch the HPC benchmark campaign for `FGSSolver(sweep_order=:dagedge)` —
the edge-level partial-pull executor, implemented and smoked green
2026-09-24. **Submission of exactly two runs is PRE-APPROVED by Ryan
(2026-09-24):**

1. a **performance run** — uninstrumented cold-solve trials, dagedge vs the
   dagteam+backoff champion; then, after it is submitted (don't wait for
   completion to prepare it),
2. a **profile run** — diagnostics-instrumented arms to verify where the
   time went (attribute costs/improvements: busy/idle/wait/reduce split,
   per-worker imbalance, task counts).

Anything beyond these two submissions (wider ladders, follow-up stages)
needs a fresh Ryan go. Read `CLAUDE.md`, `agent_policies/HPC.md`, and the
top-level `BYU_ORC_AGENTS.md` before any HPC work. This is an **official
campaign**: worktrees + annotated-tag pins + provenance file BEFORE
submitting (rules below). Local runs never >4 threads.

## What :dagedge is (context)

Edge-level split of the :dagteam lower pull: big lower edges (block bytes
= sizeof(TM)·n_i·n_j ≥ `dagedge_theta`, default 4096) become independent
partial tasks started the moment their source leaf publishes; small edges
stay in one per-leaf aggregate; the worker delivering the last partial slot
finalizes inline (fixed-order reduce + cached LU + publish). NO shared
ready queue: a deterministic list-scheduling simulation at plan build
assigns static per-worker task lists; workers wait on published-flag
predicates with the existing :spin/:backoff policies. `DagEdgePlan` wraps
an unmodified `DagTeamPlan`; iterate-preserving vs :dagteam/:lexicographic
(same certification, NO tolerance recalibration), bitwise deterministic
across workers/threads/idle at fixed θ. Design + validation record:
`BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_dagedge_design_20260924.md`.

## Why (evidence chain, all verified)

- Stage 1+2 (results `fgs_scalability_stage1_results_20260923.md` /
  `fgs_scalability_stage2_results_20260924.md`): near-field GS sweep is
  critical-path bound. Champion = **dagteam `dagteam_idle=:backoff` @ j=64
  = 3.24 s/solve**; backoff ladder 4.155/3.478/3.262 s at j16/32/64 —
  plateau does not close ⇒ pure DAG-width starvation.
- Gate-0 (`fgs_lshortening_gate0_20260924.md`): edge split shortens the
  byte-weighted critical path 289.5 → 68–69 MB/sweep = **4.25× sweep bound,
  flat in j** (θ=4KB: 4.19× at ~31.6k tasks/sweep); LU-chain floor 7.3×.
  Expected end-to-end at R4 j64: 3.24 → **~1.7–1.8 s**.
- Prototype smokes ALL GREEN (2026-09-24, local -t 4):
  `FastMultipole/test/fgs_dagedge_test.jl` 200/200 (≤1e-14 vs lex AND
  dagteam at θ∈{0,4KB,∞}; bitwise across reruns / workers {1,2,4} /
  `-t {1,2,4}` / backoff); dagteam gate-1 regression 28/28; FLOWPanel
  plumbing 20/20. Local R2 j4 overhead probe
  (`benchmark/fgs_dagedge_sweep_probe.jl`): 5.5× task count costs 7%/sweep
  (~0.3 µs/task) — no scheduler pathology; the win only appears at high j.

## Pins (create the campaign from these)

- FLOWPanel `83b8482` (branch `fastmultipole`) — :dagedge plumbing +
  gate-0 artifacts.
- FastMultipole `90a60cc3` (branch `flowpanel-20260817`) — :dagedge
  implementation.
- FLOWVPM per the new-merge-law default (≥ `8d4a3b4`; Stage 2 used
  `8d4a3b4`).
- **NOT yet pushed to origin** (gh re-auth owed): get the commits onto ORC
  the way Stage 2 did — see the deployment section of
  `fgs_scalability_stage2_results_20260924.md` /
  `fgs_acceleration_provenance_20260919*.md` for the working transfer path
  (ORC unified repo; `ssh orc` needs a live ControlMaster socket +
  `bash -lc`, 2FA otherwise). If Ryan has re-authed gh meanwhile, push
  branches + old tag `campaign/p021-fgs-stage2-20260923` + the new tag too.
- Tag convention: `campaign/p021-fgs-dagedge-20260924` (or the actual
  date), annotated, one per repo, cited with SHA in the provenance file.
  Worktrees under `/home/rander39/campaigns/*` from the tags; Manifest
  dev-paths → campaign worktrees; no uncommitted state in worktrees.

## Performance run (submission 1)

Reuse the Stage-2 harness family (launcher/analysis committed at FLOWPanel
`90f7452`; results doc names files): R4 champion knobs
(`benchmark/retained_r4_champion.toml`: P8/MAC0.4/leaf100/f32full), fixed
27×3 cold-START solve work (**cold = zero-initial-guess solves**, NOT fresh
process per arm — batch arms per (j, placement) process,
[[feedback-cold-means-cold-start]]), champion placement from Stage 1,
uninstrumented trials (diagnostics=nothing).

Arms (keep it one m12-class job like Stage 2's 13875511, 12 h):

- **Primary A/B:** dagteam(:backoff) vs dagedge(:backoff, θ=4KB) at j=64.
- **Flat-in-j check:** both executors at j ∈ {16, 32, 64} — gate-0 predicts
  dagedge's sweep time roughly flat in j while dagteam's plateaus.
- **θ probe (cheap, same processes):** dagedge θ ∈ {0, 4KB, 16KB} at j=64.
  θ is a plan-build knob — separate solver constructions in the same
  process are fine (j/placement stay process-level).
- Log at construction: `plan.ntasks`, big/small/back task counts,
  `sim_makespan`, `sim_edge_L` (bytes) per arm — free schedule metadata.

Accuracy gate: every accepted solve passes the independent evaluator as in
Stage 2 (f32full is certified by the evaluator, never bit-compared; dagedge
vs dagteam solutions are mathematically equal, not bitwise).

Decision thresholds (from gate-0 + design doc): sweep ~2.0 → ~0.5 s/solve
expected, end-to-end 3.24 → ~1.7–1.8 s. **Treat < ~3× sweep gain as
scheduler overhead eating the prize** → that's exactly what the profile run
must then attribute.

## Profile run (submission 2 — prepare while #1 queues, submit after it)

Same worktrees/pins/fixture, diagnostics-instrumented (pass the diagnostics
dict): harvest the `:dagteam_*` keys (dagedge reuses them —
busy_lower/busy_back/idle/wait/reduce/spawn/join, busy_max/min imbalance,
n_lower/n_back, empty_pops; lockmgmt is identically 0 for dagedge).
Arms: the primary A/B pair + the j ladder, instrumented. **Profile runs are
NEVER performance trials** — Stage 2 measured instrumentation itself
shifting timing (j64-w0 spin got 8.5% FASTER instrumented); rankings come
only from run #1, the profile run only explains them. Key questions:

- Where did the sweep time go: busy_lower vs idle (dependency waits) vs
  boundary wait/reduce? Does idle collapse vs dagteam at j=64?
- Per-task overhead at 31.6k tasks/sweep (busy_lower per n_lower vs the
  R2-extrapolated ~0.3 µs/task) and per-worker imbalance (busy_max/min) —
  static-schedule timing skew is the known risk; the recorded fallback is
  work-stealing deques (design doc), and NUMA first-touch repack of split
  Lmat blocks is the recorded next optimization.

## Harvest & reporting

Delegate monitoring to `hpc-monitor`, harvesting to `harvester` (compact
tables), storage/archiving to `hpc-storage` (judge runs by outputs, never
sacct; VTK/archive policy per HPC.md). Results as dated
`fgs_dagedge_benchmark_results_YYYYMMDD.md` in BRAINSTORM/021 + provenance
file filled before/at submission. Run outputs to the consolidated data
root, never into worktrees.

## Traps

- Never edit FastMultipole/FLOWPanel live checkouts while any queued or
  running job uses them — campaigns run from worktrees; `ARCHIVER_SKIP`
  worktrees under `/home/rander39/campaigns/*` are untouchable.
- `dagedge_theta` and `dagteam_workers` are construction-time (plan/list
  build); `dagteam_idle` semantics: for dagedge it's the dependency-wait
  pause policy (no queue, no qhint).
- dagteam_workers is clamped to the job's `-t`; the static schedule is
  built for that team size.
- Laptop `runtests_benchmark_cold.jl` failure is pre-existing/environmental
  (BLAS pin) — don't chase it.
- Don't resurrect prio/size-selective splitting (gate-0 killed it); the
  θ-cutoff full split is the design. Upper/backward filler + serial
  boundary reduction are out of scope (slate #4, ~0.2 s ceiling).

## Standing Ryan gates (unchanged, surface don't act)

- backoff as :dagteam production default (recommended); Stage 3
  cancellation (recommended).
- Notebook entry for Stages 1+2 + gate-0 + dagedge prototype (+ this
  campaign once harvested): OFFER, don't write.
- Origin push of branches + tags after `gh auth login -h github.com`
  (all three repos).

## House rules (binding)

Dated status/provenance files in BRAINSTORM/021; notebook Ryan-gated;
delegation per CLAUDE.md; local ≤4 threads; the two submissions above are
the ONLY pre-approved HPC submissions.
