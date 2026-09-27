# FGS scaling review and recommended next experiments

Date: 2026-09-26. Status: recommendation and experiment plan; no implementation or new benchmark campaign performed by this review.

User constraints: target up to **64 shared-memory threads**, with larger problems of interest; distributed memory is out of scope. Include both iterate-preserving improvements and changes to ordering or stationary iteration, provided they meet the same accuracy target. **No GMRES/Krylov outer iteration.**

**TL;DR: multiple workers per leaf product are worth testing.** Start with small persistent teams of 2–4 workers, retaining the current Gauss–Seidel dependencies. For the target of up to 64 shared-memory threads, this is a more focused next experiment than another fine-grained edge scheduler.

## What the existing work establishes

The champion is **`:dagteam` plus `dagteam_idle=:backoff`**. Backoff prevents idle workers from repeatedly acquiring the shared queue lock. FLOWPanel now defaults to these settings; the tuned R4 benchmark additionally uses `:f32full`, whereas the precision default remains `:f64`.

The September 24 results supersede the older status in `BRAINSTORM/INDEX.md`:

| Threads | Dagteam + spin | Dagteam + backoff |
|---:|---:|---:|
| 16 | 4.399 s | 4.155 s |
| 32 | 4.328 s | 3.478 s |
| 64 | 6.149 s | 3.262 s |

Thus, backoff fixed the high-thread-count regression, but doubling from 32 to 64 still yields only **1.07×**. See [Stage 2 results](fgs_scalability_stage2_results_20260924.md), especially the September 24 addendum. (Table note: the 32-thread spin entry of 4.328 s is the paired-comparison median; Stage 2's main table shows a 4.387 s anchor for the same point — both are correct, from different pairings.)

The dependency-chain explanation is substantially right. The current executor waits for **all lower predecessors**, then one worker performs the leaf's aggregated lower product and cached LU solve before releasing successors. The serial chain therefore contains substantial parallelizable matrix-vector work—not just inherently sequential decisions. See `FastMultipole/src/solve_dagteam.jl:317–356` in the sibling FastMultipole repository.

At R4:

- There are 1,068 leaves, with mean/max DAG level widths of 3.8/6.
- The lower-work byte-weighted work/span ratio is only 2.8.
- Lower matrix-vector products account for **92% of the modeled critical-path bytes**.
- Leaf sizes vary substantially: median 54 unknowns, maximum 1,450.

These support attacking work *inside* ready leaves. They are structural estimates, however, not measured execution-time bounds or universal limits for larger meshes. See [DAG analysis](fgs_lshortening_gate0_20260924.md).

## First recommendation: cooperative leaf products

Use two levels of parallelism: several ready leaves execute concurrently, and a small persistent team cooperates on sufficiently expensive leaf products.

**Modeled ceiling (byte model — conditional estimate, not an established bound).** Cooperating workers divide only the GEMV portion of the critical path: the node critical path is 289.5 MB = 265.6 MB GEMV + 23.9 MB LU (see [DAG analysis](fgs_lshortening_gate0_20260924.md)), so w workers give L_w ≥ 265.6/w + 23.9 MB. w=2 → 156.7 MB → ≤1.85× sweep; w=4 → 90.3 MB → ≤3.21× sweep. At R4 j64 (3.24 s solve, ~2.0 s sweep) that projects to roughly ~2.3 s (w=2) and ~1.9 s (w=4) — illustrative only, carrying all of gate-0's own caveats (infinite processors, bytes-as-cost, no per-task overhead). Splitting also changes which path is critical. **Pre-check before building any executor:** extend `benchmark/fgs_dag_L_profile.jl` (~50 lines) to recompute the critical path over the whole split graph (node GEMV cost divided among w subtasks) with a per-subtask overhead parameter swept over a plausible range; calibrate that parameter later against measured product, LU, synchronization, and reduction times from the prototype.

**Decision threshold** (mirroring dagedge's "<~3× sweep gain ⇒ overhead ate the prize"): a team-of-4 candidate must beat dagteam+backoff in paired uninstrumented j64 solves by more than ~15% wall clock, or this thread stops.

The initial experiment should:

- Compare **1, 2, and 4 workers per leaf**, within a fixed total worker budget. Keep BLAS single-threaded.
- Reuse workers throughout the inner-sweep block; avoid spawning Julia tasks for every leaf or interaction.
- Compare two product layouts: output-row partitions with disjoint writes, and contiguous column partitions with private partial vectors and a fixed-order reduction. Before integrating either into the executor, microbenchmark both kernels standalone at w=2/4 on representative leaf shapes (short-wide median, large tail), and let those measurements pick or gate the executor variants.
- Complete the product before one worker performs the existing LU solve and publishes successors.
- Select cooperative execution using measured product cost, including synchronization and reduction costs; retain ordinary dagteam execution for small products. **Split-size threshold rule (corrected 2026-09-26):** under ideal product scaling, split only when product cost exceeds $h_w/(1-1/w)$, where $h_w$ is the cooperative implementation's own measured **elapsed overhead for the whole team** (dispatch + sync + reduction on the elapsed path). The earlier $w c_{\mathrm{ovh}}/(1-1/w)$ expression assumes serialized coordination with $h_w=w c_{\mathrm{ovh}}$; summed worker cost is not generally elapsed overhead. Confirm that the complete cooperative product is faster for the measured shape/layout; ideal $1/w$ scaling is only a screening assumption. Do not seed $h_w$ from dagedge's 8.2 µs — that figure is average whole-task busy time including useful computation, reductions, atomics, and some waits (see the caveat below), not a fixed overhead. Which leaves split is an empirical outcome, not predictable in advance. See item 033 A-T3 for the elapsed/work distinction and the whole-solve acceptance requirement.
- Expose team size as a runtime toggle (like `dagteam_workers` alongside `dagteam_idle`) so paired same-process arms are possible; batch arms per (threads, placement) process. Per the Stage-2 "cold means cold-start" correction, "cold" means zero-initial-guess solves, not a fresh process per arm.
- Defer NUMA team placement to a follow-up: Stage 2 showed the j64 regression was lock contention, not bandwidth, so NUMA placement is premature complexity for the first A/B.

**Do not assume row partitioning wins.** It avoids partial-vector reductions, but median leaves have only 54 rows — though 2–4 workers can still partition 54 rows usefully, so row splitting is not disqualified for short products either. Column partitioning can divide a wide product into substantial contiguous pieces even when its output is short. Its reduction changes floating-point grouping, so the contract should be mathematical iterate equivalence plus accuracy certification. Keep both layouts experimental until the standalone kernel microbenchmarks above rule.

Also, do not restrict this to a handful of hot leaves. The earlier graph analysis found that accelerating one path exposes other paths; splitting only the top 1% of priority-ranked leaves barely changed the bound.

The distinction from the failed `:dagedge` experiment is granularity. That prototype created tens of thousands of small tasks per sweep and lost: **3.766 versus 3.540 seconds** in the paired 64-thread comparison. Its measured lower-task busy time increased 4.15×. The reported 8.2 µs/task includes computation, reductions, atomics, and some waits—it does not isolate scheduler overhead. Small cooperative teams should be tested precisely because they can divide a large aggregated product without that task explosion. See [Dagedge results](fgs_dagedge_benchmark_results_20260924.md) and the task-body timers in `FastMultipole/src/solve_dagedge.jl:328–435`.

## Other ways to expose useful parallelism

**For larger problems, investigate spatial subdomains with explicit interface treatment.** Independent subdomain work can expose more concurrency than the current global leaf ordering. Two distinct approaches deserve consideration:

1. **Preserve the GS recurrence through a different triangular solver.** Partition the lower system, solve subdomain interiors concurrently, and resolve their coupling through a reduced interface system. This trades extra setup, storage, and arithmetic for a shorter sequential path. Parallel GS using a SPIKE-based triangular solver provides a concrete precedent; its applicability here depends on interface size and fill, which should be estimated before implementation. See [Torun et al., Partitioning and Reordering for SPIKE-Based Distributed-Memory Parallel Gauss-Seidel](https://user.ceng.metu.edu.tr/~manguoglu/PDFs/dmpGS_SISC.pdf). The relevant idea is the algebraic decomposition; distributed-memory implementation is outside this plan.

   Detailed handoff: [Exact triangular-solver feasibility, build order, and test plan](fgs_exact_triangular_solver_plan_20260926.md). This starts with interface-growth and fill analysis, then a tiny algebra reference, before any production integration or campaign. Its Stage A (graph-only, laptop-scale) can run concurrently with the cooperative prototype rather than only after its failure.

2. **Change the stationary iteration to local GS plus coarse correction.** Perform parallel subdomain updates and use a coarse residual correction to communicate global error. This addresses the weakness of plain chunked GS/Jacobi: the prior experiment shortened sweeps but increased iterations from 27 to 44. A coarse correction could recover convergence while retaining parallel local work. This is a longer-term research direction, with no guaranteed convergence for this operator. Multigrid combined with multipole evaluation has precedent for boundary integral equations, but needs formulation-specific validation here. See [Caspar, Fast Solution of Boundary Integral Equations by Using Multigrid Methods and Multipole Evaluation Techniques](https://www.witpress.com/elibrary/wit-transactions-on-modelling-and-simulation/19/7526).

   Detailed phase guidance: [Overlapping local GS with coarse correction — principles, experiment order, and stopping rules](fgs_overlapping_gs_coarse_correction_strategy_20260926.md). This incorporates interface halos, weighted reconciliation, full-residual coarse correction, interface-to-interior sensitivity modes, and the role of cached local inverses.

   Implementation handoff: [Step-by-step overlapping-GS/coarse-correction prototype plan](fgs_overlapping_gs_coarse_implementation_handoff_20260926.md), including operator contracts, reference algebra, conditional enrichment, tests, and stop gates.

   Comparative prior (to be decided by measurement): coarse correction is the preferred research bet for this operator — with ~45 predecessors per leaf, spatial partitions cut edges everywhere, suggesting the exact plan's interface will be large with heavy fill, while a coarse space sidesteps that density by representing only selected global patterns — but predecessor count does not determine separator size or fill, so the exact plan's inexpensive graph-only Stage A should still run and settle whether the interface is actually prohibitive.

Neither requires a Krylov outer iteration.

Give lower priority to:

- **Another generic coloring implementation:** color-major GS was already tested; many color boundaries limited its scaling. The existing `:colored` method genuinely changes ordering; it is not merely lexicographic DAG-level scheduling. See `FastMultipole/src/solve.jl:1148–1277` and [colored benchmark](fgs_r4_followup_evidence_20260914/colored-v21-13738665/analysis/ab_summary.md).
- **Cross-sweep pipelining:** possible in principle with versioned state, but upper dependencies still constrain progress and only three inner sweeps are currently used. This would require generation-specific storage/lifetime tracking and cannot cross the outer FMM/residual boundary unchanged.
- **Parallelizing the upper boundary reduction:** worthwhile eventually, but its measured roughly 0.2 s/solve is a smaller opportunity than the lower products.

## Validation and decision criteria

The next bounded study should compare cooperative products against **dagteam + backoff**, first on R4 and then on larger meshes.

**Competitive context.** Cold R4 j64, krylov_ilu_nfcache = 2.41 s already beats dagteam+backoff 3.24 s ([dagedge results](fgs_dagedge_benchmark_results_20260924.md), "Thread-scalability thread PARKED"). In the warm-started, order-matched head-to-head ([warm-start R4 results](fgs_warmstart_r4_results_20260925.md), "Cumulative cost and crossover"), fgs_proj2 is a per-step statistical tie with ilu_nfcache_proj2 (8.246 vs 8.216 s/step median, gap ~0.4% within noise) while carrying ~210 s more setup — **no FGS break-even at any horizon**. The no-Krylov constraint is Ryan's, but the record should be plain: FGS currently trails ILU at every measured operating point, cold and warm. That is the bar any FGS speedup ultimately competes against.

- Measure product speed, synchronization, LU time, and the resulting **time-weighted critical path** before predicting a solve-level gain.
- Test zero and nonzero starts, one and multiple inner sweeps, repeatability, and mathematical agreement with the existing recurrence.
- Rank by **time to independently certified accuracy**, including iteration count. Report setup cost and memory separately.
- Benchmark 1/8/16/32/64 threads on HPC; keep local checks at four threads or fewer.
- Measure both fixed-problem scaling and how available parallelism changes with mesh size. More leaves do not automatically mean a wider dependency graph.
- Keep production defaults unchanged until a candidate wins paired, uninstrumented solves. No public API change is needed for the initial experimental prototype.

The immediate priority is **coarse cooperative leaf GEMVs**. If measured synchronization costs consume their benefit, move to subdomain/interface or coarse-correction methods rather than spending another round solely on queue scheduling.

Any future official campaign must follow the repository's HPC and reproducibility policies (see `agent_policies/HPC.md`); saving this plan does not launch or authorize a new campaign.
