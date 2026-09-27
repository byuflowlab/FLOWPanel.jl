# 021: implementation handoff — overlapping local GS plus coarse correction

Date: 2026-09-26. Status: implementation instructions for a future agent; this document does not report implemented code or executed tests.

Read the [strategy and stopping rules](fgs_overlapping_gs_coarse_correction_strategy_20260926.md) first. The [conceptual explanation](fgs_local_gs_coarse_correction_explanation_20260926.md) provides additional motivation. This handoff specifies a bounded first implementation, not every possible extension of domain decomposition.

## 1. Deliverable and nonnegotiable boundaries

Build a **private, benchmarkable stationary-solver prototype** that reuses FGS geometry, near-field matrices, and leaf LU factors, but replaces the global GS chain with overlapping local sweeps and a full-residual coarse correction. Establish whether the numerical idea works before optimizing its threading.

The first deliverable is an experimental driver usable on tiny matrix fixtures and existing FLOWPanel benchmark fixtures. Do not add a production public solver API or change defaults in this handoff. Promotion into a supported solver is separate work after the evidence exists.

- Target at most 64 shared-memory threads; local execution at most four. BLAS uses one thread in the first implementation.
- No Krylov outer solver, asynchronous updates, adaptive basis changes during a solve, full-system inverse, or exact global interface solve.
- Float64 first. Reduced precision, explicit small-block inverses, and more than two levels are later experiments, not prerequisites.
- This is a changed stationary iteration. Do **not** require its intermediate iterates to equal dagteam. Require correct fixed points, consistent residuals, and the same independent physical acceptance target.
- Full residuals and the coarse operator must use the same fixed linear operator. Do not implement the full-system correction as merely another `sweep_order` inside a frozen-far-field inner block.
- Follow current repository instructions before implementing; inspect dirty state and do not overwrite user work. Future campaigns follow `agent_policies/HPC.md`; this handoff authorizes no submission.

## 2. Read these implementation anchors before editing

In the sibling FastMultipole repository:

- `src/solve_dagteam.jl`: `build_dagteam_plan`, `dagteam_initialize!`, `dagteam_inner_sweeps!`, lower/upper layouts, leaf ranges, and existing LU references.
- `src/solve.jl`: full operator/FMM evaluation, mapping between solve vectors and system strengths, transformations, residuals, and cache lifecycle.
- `test/fgs_dagteam_gate1_test.jl`: fixtures and zero/nonzero-start conventions. Its iterate-equivalence assertions are appropriate for dagteam, **not** for this new iteration.

In FLOWPanel: `src/FLOWPanel_solver.jl`, `test/runtests_unit_solver.jl`, and the existing 021 benchmark fixtures. Read `agent_policies/TESTING.md` before running checks and the HPC policies before any remote work.

Use new private implementation/test files rather than modifying the dagteam executor. Suggested files in FastMultipole are `src/solve_coarsegs.jl` and `test/fgs_coarsegs_test.jl`; suggested FLOWPanel experiment driver is `benchmark/fgs_coarsegs_probe.jl`. These names are proposed new files, not claims that they exist.

## 2b. Gate −1: operator-application budget (before any smoother code)

The cycle in §6 uses 2–3 full operator applications per iteration versus ~1 per outer iteration today. Before writing smoother code, on the target fixture: (a) **measure** the cost of one complete operator application as this method defines it — the ~37 ms/outer figure inferred from FMM time in the [gate-0 analysis](fgs_lshortening_gate0_20260924.md) is FMM work only, not necessarily a full residual evaluation; (b) **project** total accepted-solve cost as (applications/cycle) × measured application cost × cycles + local work + coarse work + setup amortization; (c) **compare** against the measured baselines — dagteam+backoff 3.24 s cold at R4 j64 and the [warm-start ILU/FGS numbers](fgs_warmstart_r4_results_20260925.md) — not against hypothetical cooperative-team projections. Fewer than 27 cycles is not required; realistic projected total cost is the criterion (illustrative floor: 2 × 37 ms × 27 ≈ 2.0 s of application cost alone). Any go/no-go margin (e.g. projected ≥15% under baseline) is a chosen threshold, not a model consequence.

## 3. Define the operator adapter first

The prototype works with an effective square system $Ax=b$ in a single documented ordering. Inputs must provide:

| Input | Required meaning |
|---|---|
| `apply_A!(y, x)` | Overwrite `y` with the complete linear operator applied to `x`, including near/far and applicable formulation/bound-wake terms. No accumulated stale buffers. |
| `b` | Complete RHS in the same coordinates; prescribed exterior contributions included exactly once. |
| Leaf ranges and centers | Map each leaf to its unknown/equation rows, plus coordinates for partitioning. |
| `near_block(i, j)` | Read-only access to the cached near-field block mapping source leaf `j` strengths to target leaf `i` equations; absent block means zero in the local approximation only. |
| Leaf LU factors | Factors of the corresponding self blocks, with matching ordering, precision, and sign. |
| Validity token | Explicit identity/version of operator, geometry, ordering, and constraints used to build the plan. Initially rebuild for each new fixture; do not invent automatic cross-timestep reuse. |

The first real fixture must satisfy the existing dagteam assumption that leaf strength rows match leaf RHS rows. If transformations or constraints prevent the adapter from representing this consistently, stop and document the missing mapping. Do not ignore constraints or assume a constant mode is a physical nullspace.

For tiny fixtures, `apply_A!` is ordinary dense multiplication. Before using FMM, check its adapter against an assembled reference on zero, basis, and seeded random vectors; check linearity, repeatability, input preservation, and zero/nonzero-start behavior. Use accuracy consistent with the configured FMM approximation and record the deviation. Do not demand machine-precision equality from an intentionally approximate FMM operator.

Matrix-free probes can overwrite strengths and field buffers. Save/restore affected state, or isolate scratch operator state. Verify that building the coarse plan does not change the physical initial guess or RHS. Any dynamic evaluation choices that make `apply_A!` depend nonlinearly on strengths must be frozen or replaced for this prototype.

**Gate 0:** the adapter is understood and tested. Do not debug a new iteration and uncertain operator bookkeeping simultaneously.

## 3b. First diagnostic: coarse correction on the existing unmodified dagteam sweep

Once Gate 0 passes and a minimal coarse basis exists (§6 construction, without the decomposition of §4 or the local smoother of §5), run the cheapest probe of the central unknown first: keep the current dagteam sweep as the smoother and bolt the coarse step onto the **outer** iteration — after the FMM residual, solve $A_c c=P^T r$, update $x\leftarrow x+Pc$, continue. No halos, partition-of-unity weights, or threaded local runtime are needed.

Scope it honestly:

- Gate 0 applies here too: after the coarse update, strengths, source buffers, and cached near/far contributions must be made consistent before the next sweep — real integration work, not "zero smoother risk".
- Success is total cost, not iteration count alone: 27→~15 outer iterations wins only after the extra residual evaluations and coarse-solve/prolongation work per iteration are counted (Gate −1 accounting).
- Failure rules out only this configuration: global GS leaves different surviving error than partitioned local GS. On failure, record a diagnosis — basis coverage of the surviving error, projection stability, cost — and let it steer (or stop) the partitioned-smoother stages below.

Attribution rationale: the near-field-only local smoother in §5 changes two things at once (partitioning AND dropped cross-subdomain near-field coupling; the chunked-GS 27→44 result is direct evidence the latter hurts), so failures of the full design are hard to diagnose without this decoupled step. Relevant precedent, as hypothesis only: Hackbusch's multigrid of the second kind for Fredholm integral equations (see the [strategy doc](fgs_overlapping_gs_coarse_correction_strategy_20260926.md) references).

## 4. Build a deterministic decomposition

Represent a plan with core leaf lists, extended leaf lists, local/global row maps, core ownership, overlap multiplicities, worker scratch, near-block references, coarse basis/factorization, and diagnostic counters.

Initial partition algorithm:

1. Assign each leaf the work proxy $n_i^2+n_i\sum_{j\ne i,\,N_{ij}\ne0} n_j$, where $N$ is the cached near-field operator including self blocks.
2. Recursively split leaf centers on the coordinate with greatest extent, using weighted-median cuts. Resolve ties by coordinate index and then original leaf index. Restrict the initial partition-count study to powers of two and leave no empty partition.
3. Keep original leaf order within each resulting core. Never split one leaf's unknowns.
4. Form the undirected halo graph from the union of directed nonself near-field edges. Initial `halo_layers=1` includes every graph neighbor of every core leaf. Also support 0 for the control; 2 is a conditional experiment only.
5. Give each duplicated unknown weight equal to the reciprocal of the number of extended subdomains containing it. Also support core-owner-only publication as an explicit comparison mode.

Keep partition count independent of worker count. Fix a decomposition when comparing different worker counts. Report extended/core sizes, duplicated unknown ratio, predicted work imbalance, and allocation estimates before construction.

If one graph layer already expands most subdomains toward the whole problem, stop this partition/halo candidate. Do not allocate huge local matrices or silently discard couplings to make it fit. A selective-halo policy would be a separately documented follow-up.

## 5. Implement the serial local-correction reference

Use $N$ only as a local approximation; keep full $A$ in residual evaluation. For subdomain $s$, extract the near-field submatrix $N_s$ over its core plus halo through block references, not a dense full-subdomain allocation.

Given a **common full residual** $r$, initialize every local correction $d_s=0$. For each requested local sweep, visit its leaves in original order and compute

$$
d_{s,i}\leftarrow N_{ii}^{-1}\left(r_i-\sum_{j\ne i,\,j\in s}N_{ij}d_{s,j}\right).
$$

Use newly computed local values where available and the previous local values otherwise. Diagonal application uses the existing leaf LU. Start with one sweep. Two sweeps are a bounded comparison, not an automatic default.

After **all** local corrections are computed, combine them in ascending subdomain order with the chosen weights and apply one global update $x\leftarrow x+\omega d$. Do not let an earlier subdomain update the global state seen by a later subdomain. That would accidentally turn the intended parallel iteration into a sequential Schwarz variant.

First compare `omega=1`. If stability requires it, test 1/2 and 1/4 and freeze the selected value for subsequent comparisons. Any post-smoothing call uses a freshly computed full residual and newly zeroed local correction vectors.

Tests: no-overlap and overlap cases, unequal leaf sizes, an isolated leaf, directed/asymmetric couplings, a strong edge cut by a partition, correct weight sums, and zero residual giving zero correction. Compare with an independently assembled dense local-block reference. Use well-conditioned synthetic cases for tight Float64 assertions; physical convergence is a later gate.

**Gate 1:** local corrections match their defined algebra. Local-only convergence is measured, not guaranteed. Preserve this implementation as the serial reference for threading.

## 6. Add a minimal full-system coarse space

Start with one constant strength pattern per connected core component, separating bodies/unknown types where the formulation requires it. Normalize columns and remove linear dependence with a rank-revealing factorization. The first scalar fixtures use identity residual weighting, $W=I$. Physically motivated weighting is a later isolated comparison.

For each basis column $p_j$, compute $Ap_j$ with the operator adapter and assemble

$$
A_c=P^T AP.
$$

Factor $A_c$ with pivoted LU. Report dimension, estimated conditioning, and construction time. A singular or unusably ill-conditioned coarse matrix is a failed candidate; do not add an unexplained diagonal shift or silently drop physical constraints.

Initially retain sparse/local-support $P$ where possible and **do not retain dense $AP$** after construction. Recompute full residuals explicitly. This avoids silently committing the prototype to $O(Nm)$ additional storage. Record the operator-application count; caching $AP$ is a later cost comparison.

Private proposed entry points, with ordinary Julia keyword arguments rather than a new public configuration framework:

- `build_coarsegs_plan(...)`: decomposition, basis, coarse factorization, scratch, and diagnostics.
- `coarsegs_local_correction!(d, plan, r)`: defined local GS/overlap operation; does not update `x`.
- `coarsegs_cycle!(x, plan, b, apply_A!)`: one complete cycle and phase diagnostics.
- `solve_coarsegs!(x, plan, b, apply_A!; rtol, atol, maxcycles, callback)`: stationary driver with convergence/failure status and history.

One cycle is explicitly:

1. Form $r=b-Ax$.
2. Apply and combine local corrections; update $x$.
3. Recompute $r=b-Ax$.
4. Solve $A_c c=P^T r$ and update $x\leftarrow x+Pc$.
5. If enabled, recompute the full residual and apply one local post-sweep.
6. Recompute the final full residual for stopping and reporting. A valid final residual may be reused as the next cycle's input residual.

Initial stopping is $\|r\|_2\le\mathrm{atol}+\mathrm{rtol}\|b\|_2$, with explicit maximum-cycle and nonfinite-result failure statuses. Choose tolerances explicitly in fixtures. Existing dagteam calibrated internal tolerances are not automatically transferable to this different residual norm; real solves must pass the same independent BC target.

Required tests: a coarse-range error $e=Pc$ is eliminated by the coarse step alone; random nonzero starting guesses converge on controlled well-conditioned matrices; the known exact solution remains fixed; scalar scaling of RHS behaves consistently; whole-cycle results agree with an independent dense implementation. Fixed-point and final-solution agreement replace dagteam sweep-by-sweep equivalence.

**Gate 2:** a correct two-level serial reference and a reproducible comparison with local-only and global GS. If even the basic coarse exactness test fails, do not tune the method.

## 7. Conduct the bounded numerical screen before new machinery

On small representative panel fixtures, study partition counts 1/2/4/8, clipped by available leaves. Compare:

- Current global GS/dagteam reference.
- Disjoint local GS.
- One-layer halo local GS.
- The identical halo smoother plus the basic coarse correction.
- The same two-level method with one post-sweep, only if needed.

Use at least three seeded homogeneous-error starts and two distinct representative physical RHS vectors, including a nonzero initial guess. Record complete residual/error histories, transient growth, cycles, operator applications, overlap ratio, coarse dimension/conditioning, and phase times. Keep method changes separate so a gain can be attributed.

A useful result is a two-level cycle whose contraction deteriorates much less with increasing partition count than local-only GS. Do not demand that it reproduce the old 27-iteration count; compare accepted total cost.

**Gate 3:** if the simple basis works, carry it forward as the leading candidate. If it fails, diagnose whether surviving errors are poorly represented by $P$ or whether the coarse projection amplifies errors that $P$ represents well. The next two steps are conditional on that diagnosis.

## 8. Add interface-response enrichment only if needed

Define interface patches by connected groups of cut-edge incident leaves shared by each pair of core subdomains. Bundle these leaves; do not create a mode for each cut edge or panel. Use actual surface connectivity where available; document when only the near-field graph is used.

Start with a constant pattern on each patch. If justified, add up to two independent linear coordinate variations along it. Normalize local coordinates and remove degenerate modes. Keep the initial bulk/patch coarse modes as a control against errors that do not create interface disagreement.

For each local patch pattern $q$, prescribe it on that patch and zero on other local interface traces. Define $I$ as the extended subdomain's remaining unknowns and compute an extension using

$$
N_{II}h=-N_{I\Gamma}q.
$$

On tiny fixtures, solve this exactly for the reference. On real cases, apply a fixed small number of the same local GS sweeps to this response equation, initially two. Do not allocate/factor a dense entire-subdomain matrix. This extension is approximate and near-field-based; the resulting coarse operator still uses full $A$.

Assemble each extended pattern into a global column with the established overlap weights. Rank-reveal the combined basis, recompute the entire coarse operator, and compare against an equally sized basis formed by subdividing the geometric patches. This asks whether sensitivity information earns its construction cost.

If low-order patterns fail to represent surviving error, allow **one bounded enrichment round**: collect homogeneous-error snapshots, restrict them to interface patches, use local QR/SVD to select at most two additional independent patterns per patch, extend them, and validate on different seeds/RHS vectors. Persistent duplicate disagreement can supplement these snapshots but must not be the only source.

Report basis support, numerical rank, construction work, and coarse-dimension growth. Stop enrichment if it mostly adds dependent/already-damped directions, pushes cost toward an exact interface solve, or produces no useful improvement on held-out cases. Do not automatically increase the number of modes until a selected training solve converges.

## 9. Use a projection-stability diagnostic when justified

If basis-approximation tests are good but the actual coarse update is unstable, compare a QR-based residual-minimizing correction over the same fixed $P$:

$$
c=\operatorname*{arg\,min}_z\|r-APz\|_2.
$$

This diagnostic may store dense $AP$ on small fixtures. Do not form normal equations. It is not a Krylov outer iteration. It cannot increase the chosen residual norm in exact arithmetic, but the whole cycle may still fail.

If it works while $P^T AP$ fails, record a restriction/projection problem; a more suitable Petrov–Galerkin restriction is future work unless the QR cost itself is demonstrably affordable. If it also fails, additional restriction complexity is not justified without revisiting the smoother and basis.

## 10. Thread the retained method and connect real fixtures

Only now parallelize independent subdomain correction work and, where safe, independent response construction. Keep a fixed-size worker pool with task-owned scratch. Avoid nested BLAS threading and per-leaf task creation. Never key scratch solely by `threadid()`; existing dagteam code documents interactive-pool indexing concerns.

Publish local corrections only after each is complete. Combine overlapping outputs in fixed subdomain order so scheduling does not alter the arithmetic. Initially keep the small coarse solve and overlap reduction serial and measure their cost. End all worker activity before a threaded FMM operator application begins; idle spinning workers must not compete with FMM or recreate the old lock-oversaturation problem.

Verify repeated-cycle equality against the serial reference and repeatability at 1/2/4 workers for a fixed plan, basis, and Float64 arithmetic. Build the basis once when comparing worker counts so setup ordering does not confound the comparison.

The FLOWPanel experimental driver should build the ordinary fixture/cache, construct the adapter and private plan, run the new solver, restore/write strengths through existing mappings, and call the independent evaluator. Keep the current FGSSolver behavior untouched. Include cold and warm initial guesses, transformed geometry where supported, and operator-state restoration after basis probes.

Run the targeted new tests and existing dagteam regression. Run FLOWPanel solver tests if its adapter/wrapper code changes. Broaden only when the actual changes touch further contracts, following `agent_policies/TESTING.md`.

**Gate 4:** numerical and physical acceptance plus a correct threaded implementation. No speed claim before this gate.

## 11. Benchmark, decide, and stop

Produce a local pilot report first. A campaign at 1/8/16/32/64 threads is a later approved action using pinned worktrees. Hold geometry, precision, operator settings, and independent acceptance target fixed; compare same-thread-count dagteam+backoff and candidate in paired alternating-arm runs. Keep instrumented diagnostic runs separate from ranking runs.

Begin with R4, then larger meshes in the same geometry/refinement family. Freeze a partition/basis per geometry before final timing. Report:

- Total accepted solve time, cycles, full operator applications, and residual history.
- Local solve, overlap reduction, full operator, coarse solve, and prolongation costs.
- Plan construction, response construction, coarse factorization, peak/retained memory, and amortization across reusable RHS vectors.
- Behavior versus both partition count and panel count; worker occupancy alone is not evidence of scalability.

Stop a candidate when its convergence gain does not offset duplicated local work, extra residual evaluations, coarse work, and setup. Compare total time directly; a faster cycle can still lose through more cycles. Reject accuracy failures outright. If only some geometries/sizes benefit, report that limited operating region.

Explicit inverses of small leaf/augmented blocks, cached $AP$, optimized overlap reduction, mixed precision, and multilevel recursion are **follow-ups only after a numerically useful candidate exists**. A small inverse may accelerate local application but does not supply missing coarse modes. Cached response matrices $H_Q$ are preferred to full inverses for larger subdomains.

If the bounded simple-basis, interface-enrichment, and indicated projection-diagnostic comparisons fail, stop and write a negative result. Do not rescue the work by switching to Krylov, dense global inversion, unlimited halos, or unbounded basis growth.

## 12. What the next agent must leave behind

For each gate, leave a short report containing the tested hypothesis, exact configuration, correctness/accuracy status, contraction and cost evidence, and a continue/stop reason. Keep proposed bounds and measured speedups distinct.

The complete successful handoff consists of the operator adapter, serial reference, deterministic decomposition/local correction, retained coarse-space construction, private stationary driver, focused tests, and a benchmarkable FLOWPanel fixture adapter. Record unsupported formulation cases explicitly. A production solver API and defaults are not part of this first implementation.

If this document conflicts with an older exploratory suggestion, use this document for the bounded first prototype, the strategy document for guiding principles, and current repository policies for execution and campaign rules. Ask only about genuinely missing user intent; discover existing APIs and data mappings from the source.
