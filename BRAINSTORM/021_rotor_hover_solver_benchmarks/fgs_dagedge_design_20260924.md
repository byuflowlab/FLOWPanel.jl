# 021 edge-level partial pulls: `sweep_order=:dagedge` design (2026-09-24)

Design record for the L-shortening slate item #1 (sized by
`fgs_lshortening_gate0_20260924.md`: edge-level split of the :dagteam lower
pull ⇒ byte-weighted critical path 289.5 → 68–69 MB/sweep, **4.25× sweep
bound flat in j**; θ=4KB cutoff keeps 4.19× at ~31.6k tasks/sweep). This doc
settles the five design questions from `fgs_edgepull_reset_prompt_20260924.md`
and specifies the prototype implemented in FastMultipole
(branch `flowpanel-20260817`) + FLOWPanel plumbing. Local prototype +
correctness smokes only — **any benchmark is a new, Ryan-gated campaign**.

## Semantics (iterate-preserving)

Identical iterate to `:dagteam`/`:lexicographic` (mathematically, not
bitwise — same contract :dagteam already has vs lex; see Determinism below).
Today leaf $i$ waits for ALL lower predecessors, then runs one aggregated
GEMV over `Lmat[i]` (n_i × ptot_i). Instead:

- each **big** lower edge $(i,j)$ (block bytes $= s_{TM} n_i n_j \ge \theta$)
  becomes an independent task computing the partial
  $y_{ij} = L_{ij} x_j$ the moment source $j$ publishes;
- all **small** edges of leaf $i$ ($< \theta$) stay aggregated in ONE
  per-leaf task that waits for all of $i$'s small predecessors and
  accumulates their per-edge GEMVs in fixed ascending-$j$ order;
- when all of leaf $i$'s partial slots have arrived, the **finalize** step
  reduces the slots in fixed ascending slot order (big edges ascending $j$,
  then the small aggregate last), forms $x_i = b_i - (Lx)_i - u_i$, solves
  the cached diagonal LU in place, and publishes.

The upper/filler side (backward products `q_j = U_j x_j`, serial boundary
reduction into the frozen `u`) is UNCHANGED (slate #4, out of scope). The
outer-iteration RHS invariant, `lsum` bookkeeping, warm-start priming
(`dagteam_initialize!`), and residual/relaxation paths are reused verbatim
from :dagteam.

## Decision 1 — Scheduler: static per-worker task lists (primary)

Gate-0's recommendation is adopted: **no shared ready queue.** At plan build
we run a deterministic event-driven greedy list-scheduling simulation over
the edge-level DAG with exact byte costs (the DAG repeats 81×/solve, costs
are known exactly), producing one static task list per worker, ordered by
simulated start time. At runtime each worker walks its own list in order,
waiting on each task's dependency predicate before executing it:

- big-edge / back-product task → `published[j] ≥ s` (one atomic flag read;
  `published[j]` stores the sweep number in which leaf $j$ last published,
  monotone within an inner-sweep block — **no per-sweep flag reset**);
- small-aggregate task → all small predecessors' flags `≥ s`, polled in
  ascending order with a resumable index (reader-side polling; publishers
  do NO successor scatter);
- root leaves (no predecessors) get their finalize scheduled as an explicit
  static task with no dependency.

**Finalize runs inline on the worker that delivers the last arriving slot**
(per-leaf atomic arrival counter, `fetch-add`; old+1 == nslots triggers it).
This keeps the reduce+LU immediately after the last partial — the assumption
under gate-0's L bound — instead of waiting for a statically assigned worker
to reach it. The simulation models this (finalize cost charged to the
sim-last-arriving worker).

Deadlock-freedom: each worker list is ordered by simulated start time and
the simulation only starts a task after its dependency's simulated
completion, so a feasible execution exists. At runtime, suppose all workers
were blocked; take the blocked head task with minimal simulated start time —
every task in its dependency's ancestry has strictly smaller simulated start
and is either done or at/before some worker's head, giving a blocked head
with smaller simulated start: contradiction. (Asserted mechanically at plan
build: the sim verifies dependency-completion ≤ start for every scheduled
task and per-list start-time monotonicity.)

Priorities for the greedy sim: byte-weighted downstream critical path over
the edge DAG (task cost + remaining path through its target's finalize),
computed leaf-descending as in `build_dagteam_plan`. Back products get
priority 0 (filler), exactly their :dagteam role.

Fallback (NOT implemented now, kept as the recorded alternative if the
prototype shows static-schedule timing skew eating the prize): sharded
per-worker deques with work stealing + backoff idling.

NUMA: per-edge tasks pin naturally to their list's worker every sweep; the
owner-local first-touch repack of split `Lmat` blocks (Stage-1 locality,
worth ~+4.4 s) is an HPC-phase optimization — the prototype uses contiguous
column views into the existing `Lmat[i]`, which already gives stable
worker→block affinity. Recorded as the first knob to add before the campaign.

## Decision 2 — Aggregation cutoff θ

`dagedge_theta::Int` kwarg, **bytes**, default 4096 (gate-0: 4.25→4.19,
task count 49k→31.6k). Edge bytes $= \mathrm{sizeof}(TM)\, n_i n_j$.
Degenerate ends behave sensibly: θ=0 → full split; θ=typemax → every leaf is
one small-aggregate task ≈ node-level dagteam on a static schedule.
Small edges are NOT column-repacked: the aggregate task runs per-edge
GEMV-accumulates over the existing ascending-$j$ column blocks (each block
is contiguous). If profiling shows tiny-GEMV overhead matters, a small-block
column repack is the recorded optimization.

## Decision 3 — Buffers

Per-leaf partial matrix `Y[i] :: Matrix{TS}(n_i, nslots_i)`, slot = one big
edge (ascending $j$) or the small aggregate (last). Memory
$\Sigma_i n_i\,\mathrm{nslots}_i \le \Sigma_{\text{edges}} n_i$ ≈ 10–20 MB
at R4 f32 — trivial next to the split coefficient copy. Split points come
from the repack's per-target column offsets (recomputed from `preds` +
`offset`; `colofs` is derivable, not stored). Everything else (x, b, u,
lsum, q, LU cache, worker scratch) is reused from the wrapped `DagTeamPlan`.

## Decision 4 — Determinism contract

Each slot is written by exactly one task; the reduce runs only after all
arrivals, in fixed ascending slot order; each GEMV is single-threaded with
fixed per-row arithmetic order; `x_j` reads are gated by the sequentially
consistent `published` flags. Hence values are **independent of scheduling**:
bitwise-identical across repeated solves, worker counts, thread counts, and
idle policies, at every precision mode. As with :dagteam vs :lexicographic,
:dagedge is mathematically (not bitwise) equal to both — the partial-sum
grouping differs from one aggregated GEMV — and **changing θ changes the
grouping**, so runs are bitwise-repeatable only at fixed θ. Same
certification story as :dagteam (f32 modes certified by the independent
evaluator); NO tolerance recalibration.

## Decision 5 — Own sweep_order, shared kwarg family

New `sweep_order=:dagedge` (not a :dagteam flag): the executor is a
different machine (static lists, no queues) and A/Bs must name it
explicitly. It shares `dagteam_precision` / `dagteam_workers` /
`dagteam_idle` (`:spin`/`:backoff` now select the dependency-wait pause
policy — same semantics, no `qhint` since there is no queue) and adds
`dagedge_theta`. `DagEdgePlan` wraps an unmodified `DagTeamPlan` (built by
the existing constructor), so repack, validation, initialize, boundary
reduction, LU caches, and diagnostics plumbing are inherited, and the
:dagteam path is untouched. Diagnostics reuse the `:dagteam_*` keys
(harvest tooling keeps working); `lockmgmt` is identically 0, `idle_ns` is
dependency-wait time.

## Task/count bookkeeping

ntasks per sweep = #big-edge + #small-agg + #back + n_leaves (one finalize
per leaf; roots' finalize is the scheduled root task). Per-sweep reset =
slot-arrival counters + `ndone` only (published flags are monotone in the
sweep number; reset once per inner block). Team lifecycle, epoch barrier,
and coordinator-side wait/reduce timers mirror :dagteam.

## Validation (this session, local, ≤4 threads)

Test `FastMultipole/test/fgs_dagedge_test.jl`, patterned on
`fgs_dagteam_gate1_test.jl`:

1. sweep-level equivalence vs :lexicographic AND vs :dagteam, zero+nonzero
   starts, 1 and 3 sweeps, rel-dev ≤ 1e-12 (f64);
2. bitwise determinism: repeated identical sweep blocks; and
   `dagteam_workers ∈ {1,2,4}` bitwise-identical solutions at fixed θ;
   cross-`-t {1,2,4}` bitwise check via a small driver;
3. θ ∈ {0, 4096, typemax} each ≤1e-12 vs lex, each self-bitwise;
4. schedule feasibility asserts at plan build (dep-before-start, list
   monotonicity, every task scheduled exactly once);
5. end-to-end `solve!` vs :lexicographic; transformed-solver fixture;
   :f32conv/:f32full sanity at single-precision tolerance;
6. FLOWPanel plumbing smoke: `FGSSolver(sweep_order=:dagedge)` solve on a
   small body matches :dagteam.

## Prototype status (implemented + smoked this session, local -t 4)

- FastMultipole (`flowpanel-20260817`): `src/solve_dagedge.jl` (plan build +
  schedule sim + executor), `DagEdgeTask`/`DagEdgePlan` in `containers.jl`,
  `dagteam_gemv!` widened to views + accumulate variant, solve.jl dispatch,
  `dagedge_theta` kwarg. FLOWPanel: `FGSSolver(sweep_order=:dagedge,
  dagedge_theta=…)` plumbing + metadata.
- `FastMultipole/test/fgs_dagedge_test.jl`: **200/200 pass** — ≤1e-14 vs
  :lexicographic AND :dagteam at θ ∈ {0, 4KB, ∞}, zero+nonzero starts;
  bitwise across reruns, workers {1,2,4}, :spin/:backoff; end-to-end solve!;
  transformed fixture; f32 modes. Cross-`-t {1,2,4}` sha256 of strengths
  identical. `fgs_dagteam_gate1_test.jl` regression: 28/28 (no :dagteam
  change). FLOWPanel `runtests_unit_solver.jl` sweep-order plumbing
  (extended with :dagedge vs :dagteam solve match): 20/20; full file green.
- Overhead sanity (`benchmark/fgs_dagedge_sweep_probe.jl`, R2 f32full
  θ=4KB, j=4, 50 raw sweeps): 328 leaves / 10,623 edges → 5,829 tasks/sweep
  (5.5× dagteam). dagteam 4.398 ms/sweep vs dagedge 4.731 ms/sweep (−7%) —
  at j=4 both are work-bound (sim makespan 55.2 MB vs edge-L 9.7 MB), so the
  edge split can't win here; the point is per-task overhead ≈ 0.3 µs, which
  extrapolates to ~2–3 ms/sweep at R4's 31.6k tasks against a ~25→6 ms/sweep
  prize. Rel dev vs dagteam after 50 f32full sweeps: 4.8e-07 (single-
  precision-consistent).

## Benchmark (LATER, Ryan-gated new campaign)

Baseline dagteam+backoff @ j64 = 3.24 s, champion placement, R4 27×3 fixed
work, cold-START (zero-initial-guess) solves batched per (j, placement)
process. Predicted sweep ~2.0 → ~0.5 s/solve; treat < ~3× sweep gain as
scheduler overhead and profile before widening. Unmodeled risks the
prototype/benchmark must measure: per-task wait/flag overhead at ~31.6k
tasks/sweep, flop-efficiency loss of small GEMVs, static-schedule timing
skew (fallback: work-stealing deques), NUMA repack upside.
