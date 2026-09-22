# R4 FGS scalability diagnosis — plan 20260921b (rev 3, standalone)

My best version of the plan. It adopts most of the 20260921 revision's
measurement contract and causal gates (which I think are correct), and
differs where stated. Differences from the 20260921 plan are marked **[Δ]**
with the reasoning inline; everything else can be assumed shared.

## Objective, and how much this diagnosis is worth

Explain both the 16–32-thread plateau and the 64-thread regression of
`dagteam + f32full` at R4, and identify changes that reduce certified total
solve time. A useful workaround (e.g. a thread cap) and an established
mechanism are distinct outcomes; report either without overstating.

**[Δ] Effort bound tied to the alternative solver.** The 2026-09-21 harvest
shows `krylov_ilu_nfcache` at R4 delivering 4.25/3.56/2.73/2.41 s at
j=8/16/32/64 per-step (cache build excluded, ~8.5 GB cache), beating dagteam
at every measured j ≥ 8. The practical ceiling any FGS fix must aim at is
therefore ~2.4–3.6 s, not open-ended. Consequence: Stages 0–2 are worth
doing regardless (cheap, and thread-cap/placement knowledge transfers to
production); Stage 3 and the Jacobi pilot are worth building **only if**
Stages 0–2 indicate FGS headroom plausibly reaching the nfcache envelope,
or if FGS's lighter memory footprint is itself the requirement. State this
trade explicitly in the final report.

**Scope:** R4 only (58,192 panels, P8/MAC0.4/leaf100/inner3, cached leaf LU,
champion precision, BLAS=1). Existing R2/coloring records may be cited as
evidence; no new R2, coloring, or Krylov runs.

## Hypothesis slate (Stage 0 maintains this as an evidence table)

- H1 dependency width (DAG/coloring limits ready work)
- H2 queue/scheduling/waiting cost growing with j
- H3 serial sections (boundary reduction, copies)
- H4 memory bandwidth / NUMA traffic
- H5 GC/allocation pressure
- H6 core-occupancy headroom (GC/OS/runtime threads with no spare core)
- H7 iteration growth in accepted solves
- H8 placement policy wrong at high j
- H9 **[Δ]** DVFS: all-core boost roll-off from 16 to 64 active cores —
  measurable, uninterventionable, and a confound for crediting every other
  fix; must be quantified, not assumed away.
- H10 **[Δ]** SMT siblings: if the champion CPU set at j=64 forces sibling
  co-scheduling on the actual topology, the regression has a mundane
  candidate that precedes all other interpretation.

Existing evidence enters the table at zero cost: the m12 ladder
(33.07/5.97/4.50/4.42/6.42 s), the 4.475→6.223 s 16→64 accepted pair,
placement sensitivity (5.951→4.475 s at j=16 with interleave), BLAS
inertness, colored-variant regressing harder, **and the recorded R2-j64
blow-up (~9×, worse for small rungs = less work between syncs) as
prior-weighting for H2** — corroboration only, never cause. Krylov scaling
on the same nodes prioritizes probes but rules nothing out (different
kernels and working sets).

## Stage 0 — desk audit (small, local, ≤4 threads)

As in the 20260921 plan: verify the reported numbers against raw rows
(13777133/13829231 + handoffs); extract niter, allocations, GC time,
placement, calibration per row; flat niter demotes H7 (without proving
identical work); flat alloc bytes do not clear H5. Inspect the current
executor/harness only enough to answer: exact timed-solve boundary; whether
27×3 fixed work is enforceable and assertable; whether a worker-cap knob
exists and whether parked workers can avoid polling and barriers; which
dependency edges and barriers exist; what coarse timers/GC stats are
reusable.

**[Δ] Keep a cheap static DAG summary in Stage 0** (task counts, level
sizes, block sizes, max width vs j) — not as proof (level width ≠ dynamic
ready width; block size ≠ validated cost) but as an asymmetric screen: if
max width < 64 over most of the sweep, H1 is strongly promoted and Stage 2
ordering changes; if width is ample, H1 loses priority. Weighted bounds
($T_j \geq \max(W/j, L)$, validated against measured one-worker sweep time)
stay in Stage 3.

**[Δ] Resolve H10 on paper now:** champion CPU-set semantics vs the
diagnostic node's actual topology (sockets, NPS mode, SMT). Record the
answer; if j=64 implies siblings, Stage 1 gains one mandatory paired rung.

Exit: baseline manifest, reusable knobs/instruments list, updated
hypothesis table, and the smallest Stage 1 that discriminates the leaders.

## Measurement contract

Adopted from the 20260921 plan essentially whole: exclusive certified zen3
node; pinned physical cores, explicit nested CPU sets, topology-resolved
(never assume NUMA 0–3 = one socket or 64 threads = 64 physical cores);
verify page placement after warmup, not just the numactl request; fixed
memory policy across the ladder; fresh-process blocks with randomized rung
order; process medians and paired ratios (inner solves are repeated
observations, not replicates); three blocks as a screen, more for decisive
claims, CI method stated; predeclared 5% total-time screening threshold;
identical timing boundary as the baseline metric with construction/warmup/
reset/validation reported separately and never silently moved outside the
boundary; GC enabled by default, GC-disabled solves diagnostic-only
(headroom checked, restored in `finally`); fixed work = exactly 27×3 with
actual FMM/sweep counts asserted and recorded; accepted-accuracy arm kept
separate with per-thread calibration and BC rel-L2 ≤ 1e-6 certification;
first-touch discipline for placement comparisons (serial first-touch under
"local alloc" concentrates pages and is not a locality test).

**[Δ] Add per-rung effective-frequency sampling** (aperf/mperf or perf,
sampled during the timed region) as a standing column. Report raw AND
frequency-normalized speedups in every ladder. Interventions are judged
against the frequency-adjusted residual; the DVFS-owned fraction of the
plateau is reported as an operating-point fact (H9), not a defect.

## Stage 1 — reproduce + all cheap rungs and riders in ONE allocation

**[Δ] Front-load every cheap uninstrumented measurement into the single
exclusive-node job.** Rationale: once the node is allocated, an extra rung
costs minutes; deferring it to a conditional later stage costs a full
job/human round-trip. The 20260921 plan defers headroom/placement/co-run
probes to Stage 2/3 conditionals — I would run all of the following in the
Stage 1 job, since none requires code changes:

1. Fixed-work ladder at **1 / 16 / 32 / 48 / 62 / 64** (8 optional if the
   knee needs resolving). 48 localizes regression onset; 62-vs-64 reads H6
   directly. Per-solve GC stats and frequency sampling on every rung.
2. If Stage 0 flagged siblings (H10): the paired sibling-vs-spanning (or
   SMT-off) comparison at j=64, before anything else is interpreted.
3. Placement A/B at j=64 (champion interleave vs one explicit alternative,
   state rebuilt/first-touched under each policy) — H8.
4. Co-run probe: two independently prepared 32-thread processes on
   disjoint verified memory-controller domains, solo vs simultaneous —
   properly caveated (slowdown ⇒ shared-resource interference of some kind;
   no slowdown says little about one 64-thread process; never a substitute
   for bandwidth counters).
5. Accepted-accuracy arm at 16/32/64 (existing matched evidence may replace
   screening rows, not final certification).
6. In matched separate runs: coarse exclusive phase timers (init, FMM,
   influence map, residual, near-field sweeps, final update; sweep-internal
   only where naturally coarse: team lifecycle, reduction, copies). Check
   the timer sum against elapsed; ≤5% overhead with preserved scaling
   shape, else reduce granularity. No per-poll timers, no traces.

Analysis exactly as in the 20260921 plan: phase differences
$T_p(32)-T_p(16)$ and $T_p(64)-T_p(32)$ including negative contributions,
plus shortfall-from-ideal $T_p(32)-T_p(16)/2$ (a flat phase explains the
plateau at zero absolute difference); matched sample means for additive
decomposition; optimistic whole-solve saving quoted for any serial phase.
All of it also frequency-normalized (H9 column subtracted).

**Gate:** if either effect fails to reproduce, report which, compare
hardware/placement/pins/calibration/timing contract against the m12
records, and stop with a bounded reproduction report — non-reproduction on
a different exclusive node does not itself establish co-tenancy. Continue
only on a reproducible effect.

## Stage 2 — code-change interventions, gated on implementation effort

**[Δ] The gating resource here is engineering effort and correctness risk,
not node time** — so conditionality applies to code changes, not rungs.
Select from the 20260921 plan's probe table based on Stage 1; do not run a
Cartesian product. At j=64 with a j=16 control; confirm winners at 32.
Order by cost:

1. Near-field worker cap (16, then 32 only if useful) with full thread
   count for FMM; parked workers must not poll or join active barriers;
   parking/wakeup cost accounted. Recovery localizes loss to sweep
   participation without alone separating H1/H2/H4.
2. Bounded backoff replacing empty-queue busy polling (dependency
   visibility, progress, and arithmetic order preserved; starvation/
   deadlock and low-thread regression checked).
3. GC intervention if H5 is alive (one supported gcthreads setting or
   targeted allocation removal, compute threads fixed).
4. Target-parallel boundary reduction — only if Stage 1 timers show the
   reduction owning enough time to justify it.

Every semantics-preserving variant passes the 1e-8 repeat-solution gate and
reproduces baseline residual history; algorithm-changing variants certify
independently. A thread cap may cure the regression while leaving the
plateau unexplained — keep the conclusions separate.

**Decision gate (with the nfcache bound in view):** if measured loss +
successful interventions explain the effects within uncertainty, stop and
report. If material ambiguity remains, proceed to Stage 3 only with a named
unresolved question, a probe that could change the recommendation, **and a
credible path for FGS to approach the nfcache envelope (or a stated
memory-footprint rationale)** — otherwise recommend the operating point +
solver-selection trade and stop.

## Stage 3 — resolve remaining ambiguity only

As in the 20260921 plan, unchanged in substance:

- Per-worker preallocated padded aggregates first (useful work by kernel,
  lock/queue wait, empty wait, boundary wait); sampled task timeline only
  if aggregates cannot separate "no ready work" from "failure to dispatch";
  occupancy/imbalance/availability derived without summing overlapping
  worker time; observer overhead and scaling shape re-validated.
- Weighted DAG bound $T_j \geq \max(W/j, L)$ with measured task costs, $W$
  checked against one-worker useful sweep time, $L$ compared against the
  multithreaded sweep (never the whole solve), sensitivity bounds on costs.
- Frozen-input replay (actual blocks, static balanced schedule, gather/copy
  included, consumption checked, matched precision/layout/placement/warmup)
  as an empirical throughput reference — built only if it decides
  dependency/executor vs kernel-throughput. Bandwidth claims only with
  phase-scoped memory-controller counters and a replay shown to sit at a
  stable throughput limit; otherwise the saturation claim stays unresolved.

## Algorithm follow-up — only with demonstrated headroom

Block-Jacobi pilot exactly as specified in the 20260921 plan (same outer
loop/partition/matrices/cached LU/precision; separate write buffer;
relaxation 1.0/0.8/0.6/0.4; 300-iteration cap; screen 16/64, extend
promising settings to 32), judged on certified total solve time. **[Δ] Its
go/no-go additionally cites the nfcache bound:** a Jacobi variant projected
to land above ~2.4–3.6 s at high j is not worth certifying unless the
memory-footprint argument applies. No ILU–Krylov replacement in scope.

## Validation, logistics, deliverable

- Narrow FGS/solver tests (per `agent_policies/TESTING.md`) + ≤4-thread
  local smoke incl. multi-thread progress checks for cap/backoff before any
  HPC submission; high-thread correctness/performance confirmed only on HPC.
- Campaign rules: clean annotated-tagged worktrees for FLOWPanel and all
  dev dependencies; Manifests pointed at worktrees; loaded SHAs recorded;
  outputs to the consolidated data root; source frozen while jobs run.
- Budget Stage 1 from a measured pilot (startup + solves × reps), not
  assumptions; checkpoint completed rows; spend reruns on ambiguous or
  winning comparisons.
- Deliverable: one compact report — provenance and raw-row locations; raw
  and frequency-normalized ladders with CIs; phase decomposition of both
  effects; rider results (GC, headroom, SMT, placement, co-run); paired
  intervention effects with uncertainty; each hypothesis classified
  supported / demoted / unresolved with basis; and a recommendation stated
  as production advice: either an FGS operating point (threads, placement,
  any adopted change) or the explicit finding that FGS cannot approach the
  nfcache envelope at R4 and solver selection should proceed on the
  memory-vs-speed trade.

**Done** when both effects have a quantitative explanation and a controlled
recovery consistent with it, and any proposed improvement wins on certified
total time — or when the report honestly states the workaround, the
remaining ambiguity, and the smallest discriminating experiment. Never
manufacture a root cause to satisfy the stopping rule.
