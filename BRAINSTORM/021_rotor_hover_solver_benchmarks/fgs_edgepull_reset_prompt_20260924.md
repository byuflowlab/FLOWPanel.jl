# Reset prompt: 021 edge-level partial pulls — design + prototype (2026-09-24, supersedes fgs_lshortening_reset_prompt_20260924.md)

Copy everything below the line into a fresh agent session started in
`~/Dropbox/research/projects/FLOWPanel.jl`.

---

## Mission (Ryan 2026-09-24)

Implement **edge-level partial pulls** for
`FGSSolver(sweep_order=:dagteam)` — the L-shortening slate's item #1, now
sized and confirmed as the next step by gate-0. Design doc + local
prototype + correctness smokes first. **Any HPC benchmark run is a NEW
campaign and Ryan-gated** (worktree/annotated-tag pins, provenance file,
explicit go). Read `CLAUDE.md` + `agent_policies/HPC.md` before HPC work;
local runs never >4 threads.

## Why (evidence chain, all verified)

- Stage 1+2 (jobs 13858983/13875511/13878514, results in
  `fgs_scalability_stage1_results_20260923.md` /
  `fgs_scalability_stage2_results_20260924.md`): the near-field GS sweep is
  **critical-path bound**, not bandwidth bound. `dagteam_idle=:backoff` @
  j=64 = 3.24 s is the current best; backoff ladder 4.155/3.478/3.262 s at
  j16/32/64 — adding threads saturates ⇒ pure DAG-width starvation.
- **Gate-0** (`fgs_lshortening_gate0_20260924.md`, script
  `benchmark/fgs_dag_L_profile.jl`, evidence CSVs in
  `fgs_lshortening_gate0_20260924/`): byte-weighted critical path of the R4
  champion plan (1,068 leaves, 48,167 lower edges, depth 279) is
  L_node = 289.5 MB/sweep, 92% pull-GEMV bytes. Edge-level split ⇒
  L = 68.1 MB, **4.25× sweep speedup bound, flat in j** (LU-chain floor
  7.3×). Expected end-to-end at R4 j64: 3.24 → ~1.7–1.8 s.
- **The catch that shapes the design:** the path is broad. Greedy
  path-driven splitting converges only at 815/1,068 leaves; prio-top-20%
  gives just 1.44×. Best practical point: **θ=4KB min-edge-bytes cutoff =
  4.19× at ~31.6k tasks/sweep** (~30× today; mean edge task ≈16 KB ≈1–2 µs
  of streaming, ×81 sweeps/solve ≈ 2.6M tasks/solve). Today's single
  SpinLock ready-queue (2–6% lockmgmt at 1,068 tasks) will not survive
  this — **the scheduler IS the design problem.**

## What to build

Semantics (iterate-preserving — mathematically identical to
:dagteam/:lexicographic, same certification, NO tolerance recalibration):
today leaf i waits for ALL lower preds then does one aggregated GEMV
(`dagteam_do_lower!`); instead compute partials y_ij = L_ij·x_j as
independent tasks the moment pred j publishes; when all of leaf i's
partials exist, reduce them **in fixed ascending-j order** (same trick as
the boundary reduction), then the leaf LU solve + publish.

Design decisions to settle in the doc (gate-0's recommendations):
1. **Scheduler**: promote the parked static-list-schedule idea to primary
   candidate — per-worker static task lists (exact byte costs known; DAG
   repeats 81×/solve so schedule once at plan build) + per-leaf atomic
   arrival counters, NO shared ready queue; enables NUMA owner-local
   first-touch of the split Lmat blocks (locality worth +4.4 s at Stage 1).
   Fallback: sharded queues / work-stealing deques with backoff idling.
2. **Aggregation cutoff**: edges < θ≈4KB stay in a per-leaf aggregated
   residual GEMV (gate-0: costs only 4.25→4.19). Make θ a kwarg.
3. **Buffers**: per-edge y_ij output space (Σ over split edges of nof(i)
   floats; allocate like `q` does per upper block). Repack already stores L
   in per-(i,j) column blocks — split points exist (`colofs`).
4. **Determinism contract**: fixed per-row arithmetic order + fixed
   ascending-source reduction order; scheduling must never affect values.
5. Keep :spin/:backoff/`dagteam_workers` semantics coherent (or make the
   static schedule its own `sweep_order=:dagedge` / plan variant — decide
   and justify).

## Where the code lives

- FastMultipole checkout `~/Dropbox/research/projects/FastMultipole`,
  branch `flowpanel-20260817` @ `053c8de7` (dev-pathed from FLOWPanel's
  Manifest). `src/solve_dagteam.jl` (606 lines): `build_dagteam_plan`
  (repack, `colofs`, prio = byte-weighted downstream path),
  `DagTeamPlan` struct in `src/containers.jl:1156` (preds/nsucc/uppers/
  prio/ptot/Lmat/Umat/red/offset/x/b/u/lsum/q/lus/readyQ/backQ/qlock/
  qhint/stats), executor + `dagteam_do_lower!`/`dagteam_pop_ready!` in
  solve_dagteam.jl, construction trigger in `src/solve.jl:760–777`.
- FLOWPanel plumbing: `src/FLOWPanel_solver.jl:1481–1557` (FGSSolver kwargs
  → FastGaussSeidel; forward new kwargs only when requested, pattern at
  1542).
- Gate-0 analysis: `benchmark/fgs_dag_L_profile.jl` (runs local:
  `RUNG=R4 SKIP_B=1 THREADING_MODE=multi EXPECT_JULIA_THREADS=4
  BENCH_BLAS_THREADS=8 julia --project=. -t 4 ...`; fixture ~2.5 min;
  BENCH_BLAS_THREADS=8 is the laptop OpenBLAS quirk). Graph dump CSVs in
  `fgs_lshortening_gate0_20260924/` allow offline schedule prototyping
  without the fixture (see `dagL_offline.jl` there).

## Validation contract (before any benchmark ask)

- 4-thread local smoke: solutions **bit-identical** to `:lexicographic`
  and to current `:dagteam` (f64; f32full certified by evaluator, not bit
  compare). Pattern: prior session's `stage2_smoke.jl` (46 tests) — R1/R2
  rungs are laptop-sized via `phase1_case.jl`.
- Determinism repeat: same answer across repeated solves and across
  nthreads ∈ {1,2,4}.
- Static-schedule sanity: every edge task's source publish precedes it in
  its worker list order dependencies (assert at plan build).
- `test/runtests_benchmark_cold.jl` fails on the laptop (BLAS pin,
  pre-existing environmental) — don't chase it.

## Benchmark (LATER, Ryan-gated new campaign)

Baseline = dagteam+backoff @ j64 (3.24 s), champion placement, R4 27×3
fixed work, cold-START solves (zero initial guess — Ryan's meaning of
"cold"). Iterate-preserving ⇒ fixed-work contract reusable directly. Build
the batched-arm harness with it ([[feedback-cold-means-cold-start]]: arms
share a process per (j, placement); runtime-toggleable knobs). Predicted
observable: sweep ~2.0 → ~0.5 s/solve; treat < ~3× sweep gain as scheduler
overhead eating the prize and profile before widening the campaign.

## Live loose ends

- Resume job 13878514: DONE, harvested, j32 addendum appended to the Stage 2
  results note; run dirs `fgs-stage2-13875511` (+ stage-1 13858983)
  archive-eligible once quiet — route `hpc-storage`.
- **Ryan-gated standing:** backoff as production default (recommended);
  Stage 3 cancellation (recommended); notebook entry for Stages 1+2 +
  gate-0 (offer, don't write).
- Origin push owed after `gh auth login -h github.com`: branches + tag
  `campaign/p021-fgs-stage2-20260923`, all three repos.
- Uncommitted: `benchmark/fgs_dag_L_profile.jl`, gate-0 note + evidence
  dir, Stage-2 results addendum (commit with the next 021 commit).

## Traps

- Never edit FastMultipole source while any queued/running job uses its
  live checkout — check the queue first; campaign worktrees under
  `/home/rander39/campaigns/*` are ARCHIVER_SKIP, never touch.
- `ssh orc` needs a live ControlMaster socket + `bash -lc`; judge runs by
  outputs never sacct.
- prio/size-based *selective* splitting is a dead end (gate-0 measured it)
  — don't resurrect it as a "cheaper first version"; the θ-cutoff full
  split is the design.
- Upper/backward filler products and the serial boundary reduction are NOT
  in scope (slate #4, ~0.2 s ceiling; only after edge pulls land).
- Gate-0's L is an infinite-processor bound in bytes: flop-efficiency loss
  of small GEMVs and per-task overhead are unmodeled — that's what the
  prototype must measure.

## House rules (binding)

Monitoring `hpc-monitor`; harvesting `harvester`; storage `hpc-storage`;
BRAINSTORM catch-up `brainstorm-scout`; notebook Ryan-gated; dated
status/provenance files in BRAINSTORM/021. NO HPC submission approval is
currently outstanding — every submission needs Ryan's explicit go.
