# R4 colored-sweep experiment: context-reset handoff (2026-09-17)

Prepared 2026-09-17. Supersedes `fgs_r4_context_reset_20260916b.md` (its
diagnostics pass is COMPLETE: job 13733332 succeeded, §7 updated, silo
deleted, merges done). Ryan has directed the colored-sweep experiment
(2026-09-17): this doc is the execution handoff for that new campaign phase.

## Keep the parent context small

Delegate bounded mechanical work to the repo subagents in `.claude/agents/`
(`hpc-monitor` read-only status/logs, `harvester` tabulation, `test-runner`,
`brainstorm-scout`, `code-scout`), cheapest reliable model, ≤40-line summaries
with indexed artifacts. Keep inline: physics/numerics reasoning, conclusions,
code edits, job submission, anything needing Ryan. Re-arm your own job watch
(lesson: parse sacct output defensively — cluster banners pollute stdout; use
`-n` AND filter for known state words, never `head -1`).

## Required policy reads

`~/.claude/CLAUDE.md`, repo `CLAUDE.md` (+ `agent_policies/{WORKFLOW,TESTING,
HPC}.md` as routed). Before cluster work: current
`/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` (+ slurm/storage
subdocs). ≥60 s between scheduler queries. Local jobs ≤4 threads. `ssh orc`
needs a live ControlMaster socket with BatchMode — auth failure = STOP (never
MFA); a sandbox denial needs escalated SSH and is not an auth failure. `TZ=UTC`
for Slurm timestamps; `sbatch/sacct` need `source /etc/profile` over
non-interactive ssh.

## State at handoff (2026-09-17)

- **Diagnostics phase closed.** Measured conclusions = §7 of
  `fgs_opt_r4_diagnostics_package_20260912.md` (READ §7 FIRST, then §6
  gates). Headline: `nearfield_update` measured at avg ≈1.00 active thread
  at BOTH j4 and j64 (non-scaling ~9.5 s span, direct /proc evidence, job
  13733332); FMM far field scales 24×; whole-solve ceiling 3.43×; DRAM
  saturation formally unresolved (generic-counter limitation, recorded).
- Evidence: `fgs_r4_followup_evidence_20260914/` — `diag-v15-13694724/`
  (ladder), `counters-v20-13733332/` (counters/activity + `analysis/`),
  `counters-v{17,18,19}-*-FAILED/` (failure harvests), `v*-deployment/`
  (provenance). All SHA256-verified local. Cluster silo
  `FLOWPanel-p021-r4-diag-v10-silo` DELETED (pre-authorized).
- **Merges landed (2026-09-17, NOT pushed):** FLOWPanel `dad2ceb` merges
  `campaign/p021-r4-counters-v17-wt` (v10–v20 lineage incl. perf-FIFO fixes)
  into `fastmultipole`; FastMultipole merges `adb9967d` (campaign v11 tag:
  FGS stage instrumentation + activity observers + solve! diagnostics
  kwargs) into `flowpanel-20260817`. Validated locally: FLOWPanel
  `runtests_unit_solver.jl` (419 pass), `runtests_unit_fgs_history.jl`,
  `runtests_r4_counters_driver.jl`, and the FMM chain
  gravitational+solve_test+fgs_coloring_test (157 + 2216 colored-sweep
  cases) all PASS.
- Local campaign worktrees still exist:
  `/private/tmp/flowpanel-p021-r4-counters-v17` (at v20 tag) and
  `/private/tmp/fastmultipole-p021-r4-activity-v11` (at v11 tag). Now that
  both lines are merged into the dev branches, cut NEW worktrees from NEW
  tags for v21 — do not reuse these.

## The experiment

**Question:** does a conflict-free colored leaf-sweep order parallelize the
~9.5 s serial `nearfield_update` chain (85% of j64 wall) and beat serial
lexicographic on total time to accepted accuracy?

**Implementation status — mostly DONE, do not re-implement:**
`FastMultipole.solve!` already accepts `sweep_order ∈ {:lexicographic,
:colored}`; `color_leaves` (src/solve.jl ~985) builds the conflict graph
(leaf L conflicts with every leaf whose own RHS rows overlap rows L writes
through its direct near-field entries, symmetrized) and greedy first-fit
colors it in ascending leaf order; `leaves_by_color` drives the sweep. The
2216-case `test/fgs_coloring_test.jl` suite passes on the merged line and ran
as a v20 cluster control. On the R4 production tree the census measured
**79 colors** (1,068 leaves, 48,627 undirected conflicts, color sizes 1–27,
median 16, verified conflict-free). 79 is an OUTPUT of greedy first-fit on
that tree's conflict graph, not a tunable: it is bounded below by the graph's
clique number (dense row-overlap clusters near the root/large leaves) and
above by max-degree+1; a different mesh/tree gives a different count. Expect
sync cost ≈ 79 color barriers × 81 inner sweeps ≈ 6.4k barriers/solve.

**What remains is the campaign phase:**

1. Harness plumbing: the cold harness (merged `benchmark/fgs_cold_common.jl`)
   carries `sweep_order` as a config axis; verify `"colored"` value plumbs
   end-to-end (config → `cold_make` → solve!). The superseded 004ce84 WIP
   snapshot had a `"sweep_order" => ["colored"]` draft axis — re-derive, do
   not resurrect the snapshot.
2. **Separate calibration.** Coloring changes accumulation ordering, so the
   iterate path differs: iteration count may move off 27 and tolerances may
   need retuning. Rank ONLY by total time to accepted accuracy under the §6
   gates (BC rel-L2 ≤1e-6, certified-FMM authoritative, repeat ≤1e-8, finite,
   BLAS=1) — never by per-sweep speed.
3. Design: A/B lexicographic vs colored on the retained R4 config, thread
   arms at least j∈{1,4,16,64}, uninstrumented performance trials in
   alternating batches (v15 pattern: 40 trials/arm, batches of 10), plus j1
   equivalence controls (direct-vs-FMM, both orders). Include one
   activity-instrumented arm pair (j64) to verify the chain actually runs
   parallel under :colored (avg active threads over `nearfield_update`) —
   budget attribution only, excluded from rankings.
4. Success yardsticks: j64 uninstrumented wall vs the 10.96 s lexicographic
   median; chain Amdahl bound if coloring were free ≈ (10.96 − 9.3) + 9.3/16
   ≈ 2.2 s (median 16-way) — real gain will be less (barriers, small colors,
   cache effects). Also collect j4: coloring may pay earlier at low j.
5. Campaign hygiene (all standing rules apply): commit + annotated tags
   (`campaign/p021-r4-colored-source-20260917-v21` in FLOWPanel; new
   FastMultipole tag only if its source changes — else retain `...-v11`,
   but note both repos now carry the merged lineage on dev branches, so cut
   v21 tags from the MERGED branches). Fresh git worktrees from the tags —
   per Ryan 2026-09-05 (no-more-silos): deploy on the cluster as **pinned
   git worktrees, not rsync silos**; Manifest dev-paths at the campaign
   worktrees; pins in the provenance file pre-submit; outputs to the
   canonical data root `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/`.
6. **Pre-submit gate protocol (v17–v19 lesson chain, mandatory):** execute
   every script the launcher will invoke, locally, under the cluster's Julia
   minor line (julia-1.11.x — cluster module is julia/1.11.7), and if any
   external binary protocol is involved (perf etc.), validate against the
   real binary on the login node, not a mock. Uninstrumented perf trials
   don't need perf at all — prefer dropping perf from the v21 launcher
   entirely.
7. Storage preflight (<400 G target), ≥300 s sacct spacing with defensive
   parsing, judge runs by outputs. Storage state 2026-09-17: an hpc-storage
   archive sweep + the silo deletion brought /home from 725 G to ~590–625 G;
   the detached archiver worker (log
   `/home/rander39/archiver_apply_main_20260916_193658.log`) finished clean
   (DF_AFTER=580G, zero locked checkouts) but usage is STILL over the 400 G
   cap — run another `hpc-storage` sweep before submitting long jobs.

## On failure

Harvest evidence FIRST (pattern: `counters-v*-FAILED/` harvests +
`sha256sum -c`), root-cause before any rerun (reproduce locally when
possible — v18/v19 postmortems are the template), fix = new commit +
annotated tag + fresh deployment; never edit deployed sources or move tags.

## Stop conditions

Report A/B results to Ryan before deriving follow-on experiments. No
notebook entry without Ryan's separate approval (the diagnostics
ladder-findings entry is still owed and drafted on request). Merged-but-
unpushed branches: do not push without Ryan.

## Document index

All under `BRAINSTORM/021_rotor_hover_solver_benchmarks/`:
1. `fgs_opt_r4_diagnostics_package_20260912.md` — §7 = measured profile +
   counters/activity findings; §6 = gates.
2. `fgs_r4_followup_validation_20260915.md` — chronological log; tail =
   2026-09-17 completion entry.
3. `fgs_r4_followup_evidence_20260914/` — all evidence + provenance.
4. Memory: `project_021_solver_benchmarks.md`.
