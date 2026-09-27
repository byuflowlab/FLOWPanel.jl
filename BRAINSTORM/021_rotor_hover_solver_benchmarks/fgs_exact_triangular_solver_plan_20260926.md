# 021: feasibility and prototype plan for an exact alternative GS triangular solve

Date: 2026-09-26. Status: planning only; no solver implementation or tests performed for this document.

Parent: [FGS scaling review and experiment plan](fgs_cooperative_leaf_scaling_plan_20260926.md).
Motivation: [explanation of what a different triangular solver can change](fgs_triangular_solver_scaling_explanation_20260926.md).

## 1. Objective and boundaries

Determine whether the existing GS lower solve can be replaced by **parallel subdomain solves plus a smaller exact interface solve**, improving time to solution on up to 64 shared-memory threads as panel count grows.

The central question is not whether there are more leaves or more scheduled tasks. It is whether we can shorten the global dependency chain without paying more in interface work, memory traffic, or setup than we save.

This is a feasibility investigation, not a promised optimization. A well-supported rejection is a successful outcome. The implementer should follow the stages below and stop at a failed gate rather than build the entire solver first.

Keep the current leaf order, leaf blocks, near/far split, outer FMM iteration, relaxation, stopping criteria, and number of inner sweeps fixed during the first comparison. Preserve the mathematical GS update; bitwise agreement across different algorithms is not required. Do not add Krylov iterations, approximate interface solves, dropped couplings, new GS ordering, or coarse corrections to rescue this experiment. Those would answer a different question.

Use `:dagteam` with `dagteam_idle=:backoff` as the baseline. Start in Float64. The tuned R4 performance baseline uses `:f32full`, whereas the production precision default is `:f64`; compare matching precision modes and identify them explicitly.

## 2. The algorithm to investigate

### Existing recurrence

One inner sweep is the block lower-triangular solve

$$
T x^{s+1}=r^s,\qquad T=D+L,\qquad r^s=b-Ux^s.
$$

Here, $D$ contains the existing dense leaf self blocks, $L$ contains the existing strictly lower near-field blocks, and $b$ includes the contribution held fixed during the inner-sweep block. The current dagteam executor evaluates this through leaf forward substitution. Its source-major upper products supply the next sweep's frozen upper contribution.

### Partition without changing the recurrence

Group **contiguous ranges of existing leaf indices** into $K$ subdomains. Do not split a leaf. This makes the first implementation easy to audit and preserves the existing GS order. The tree ordering may provide spatial locality, but measure interface size rather than assume it.

Let $B$ retain every block of $T$ whose source and target are in the same subdomain. Thus $B$ is block diagonal over subdomains, but each diagonal subdomain block is itself a leaf-block lower-triangular system. Let $E=T-B$ contain all cross-subdomain lower couplings.

Define the interface set $C$ as **all unknowns of every source leaf with at least one outgoing cross-subdomain lower edge**. Include the whole leaf, even if individual coefficient entries happen to vanish. This is a structural definition; never drop small numerical entries. Store $C$ in original global order.

Let $R_C$ select those entries from a full vector, and let $E_C$ contain the corresponding columns of $E$. Then

$$
E=E_C R_C.
$$

Define

$$
H=B^{-1}E_C,\qquad S=I+R_C H.
$$

For each sweep:

1. Solve $B y=r^s$ independently in every subdomain.
2. Solve the reduced interface system $S z=R_C y$.
3. Recover the full solution as $x^{s+1}=y-Hz$.

The identity follows by substituting $z=R_Cx$ into $Bx+E_Cz=r$. This is an exact reformulation in exact arithmetic; it does not freeze cross-subdomain values at the old iterate.

Because the original system is lower triangular in **leaf blocks**, and the partition preserves that order, $S$ is lower triangular in selected leaf blocks with identity diagonal blocks. The reduced solve still has dependencies. Its size and dependency depth are the main feasibility questions.

### Storage choice for the first practical prototype

Use a tiny dense reference implementation of $H$ to verify the derivation, but **do not store full dense $H$ for a real panel case**.

The practical prototype stores only the necessary interface rows $R_CH$, represented as structurally nonzero leaf blocks of $S$. Build these by subdomain forward substitution on batches of incoming interface columns, retaining the selected output rows and discarding other temporary response rows.

Recover the full solution using a second set of parallel subdomain solves:

$$
x^{s+1}=B^{-1}(r^s-E_Cz).
$$

This replaces storage and application of full $H$ with another local solve. Include that second solve in every cost estimate. Use existing leaf LU factors for all subdomain solves; do not form an explicit inverse or dense LU of an entire subdomain.

After recovery, use the existing upper-product machinery to prepare the next sweep. The first prototype need not overlap this work, but must account for losing the baseline's within-sweep upper-product overlap.

## 3. Stage A — graph-only feasibility before numerical construction

Build a standalone analysis script, proposed name `benchmark/fgs_triangular_interface_profile.jl` in FLOWPanel. Reuse the fixture/plan-building approach of `benchmark/fgs_dag_L_profile.jl`. Do not assemble a global dense matrix.

Inputs: an existing FGS plan, leaf ordering and sizes, lower predecessor lists, and requested partition counts. Initial partition counts are $K\in\{1,2,4,8,16,32,64\}$, excluding counts larger than the number of leaves.

Choose contiguous cut positions using cumulative per-leaf weight

$$
w_i=n_i\,\mathrm{ptot}_i+n_i^2,
$$

with deterministic ties and nonempty ranges. This is a cheap initial work proxy, not a calibrated timing model. Keep the partition count separate from the number of runtime threads; the same partition must be comparable at several thread counts.

For each partition:

- Mark cross-subdomain edges and interface source leaves.
- Propagate each incoming interface-column block through the target subdomain's local lower DAG to predict the structural nonzero blocks of $S$. Count fill caused by this propagation, not just original cross edges.
- Record interface unknowns and leaves, interface fraction, subdomain imbalance, predicted $S$ storage in Float64/Float32, construction scratch, and extra resident bytes relative to the existing dagteam plan.
- Record work and weighted path estimates for both local solves, the interface solve, cross-product application, and upper work. Include the sequential ordering of these phases; do not add parallel worker times as if they were elapsed time.
- Record original lower-DAG work/span and depth alongside the transformed estimates. A small interface fraction alone is not a pass.

Run the structural analysis first on existing R2 and R4 fixtures. Then extend the same geometry/refinement family toward approximately 2× and 4× R4's panel count, recording actual counts and all meshing parameters. If a larger fixture is unaffordable, report that limit and leave the large-problem conclusion unresolved. Do not silently substitute unrelated geometry or independent duplicated bodies to manufacture parallelism.

Outputs: one CSV row per geometry/partition, a short report of growth with panel count, and a selected candidate partition for the next stage. Keep graph-only byte/work estimates clearly separate from measured seconds.

**Gate A:** reject partitions whose interface is essentially the original problem with extra work, whose predicted storage is unaffordable for the intended run, or whose reduced path removes no meaningful bottleneck. If no partition shows a plausible advantage across increasing problem sizes, stop and report why. A marginal model is inconclusive, not evidence of speedup. Do not proceed directly to a campaign from this stage.

## 4. Stage B — prove the algebra on tiny systems

Build a standalone dense reference test in FastMultipole, proposed name `test/fgs_triangular_interface_test.jl`, before touching the production executor. Use seeded, well-conditioned block-lower matrices with nonsymmetric invertible diagonal blocks. For each case compare:

- Direct solution of $Tx=r$.
- Ordinary leaf forward substitution.
- Interface solution with explicit $H$ and recovery $y-Hz$.
- Interface solution with recovery $B^{-1}(r-E_Cz)$.

Use normalized forward error and normalized residual. For the well-conditioned Float64 fixtures, require both at or below $10^{-12}$. Report conditioning when a deliberately difficult fixture is added; do not loosen the core fixture tolerance to hide an algebraic error.

Required cases:

- One subdomain and an empty interface.
- Independent subdomains with no cross edges.
- A chain crossing every partition boundary.
- A source feeding several subdomains and several sources feeding one target.
- Unequal leaf sizes, single-leaf subdomains, and a nearly all-interface graph.
- A case where local propagation creates interface fill absent from the original cross-edge graph.
- Multiple RHS vectors, including zero and nonzero values.

Check structural identities $T=B+E_C R_C$, the selected interface ordering, and the identity diagonal of $S$. Check that the symbolic pattern from Stage A contains every numerical nonzero generated by the reference construction.

**Gate B:** all identities and solution comparisons pass before a threaded implementation begins.

## 5. Stage C — practical standalone executor and measured cost

Build an internal experimental plan in FastMultipole, proposed implementation file `src/solve_triangular_interface.jl`. It should contain partition ranges, local lower-block access, interface index maps, sparse leaf-block storage for $S$, temporary vectors, and references to existing leaf LU factors. No public solver keyword is needed yet.

Implement in this order:

1. Sequential construction of the reduced blocks, with bounded batches of interface columns and no global dense matrix or full $H$ allocation.
2. Sequential execution of the two-local-solves algorithm, checked against Stage B and dagteam lower solves.
3. Parallel execution across subdomains for construction, the first solve, cross-product application, and recovery. Use at most the available Julia thread budget and single-thread BLAS. Each local subdomain solve initially remains serial in its existing leaf order.
4. A deterministic serial block forward solve for $S$ as the reference reduced solver. Measure it explicitly; do not assume that making local work parallel makes this phase negligible.

Reuse buffers and task assignments across sweeps. Give each output a single owner and use fixed source accumulation order. Require repeatability for a fixed plan across repeated runs and 1/2/4 workers. Different partitions or precision modes need numerical equivalence, not bitwise identity.

Measure construction, first local solve, interface solve, cross products, recovery, and allocations independently. Record peak memory as well as retained storage. Replay the same frozen RHS and coefficients for candidate and baseline. Include whole-sweep measurements; isolated fast local solves are insufficient.

**Gate C:** proceed only if measured costs support a useful whole-sweep improvement at an affordable memory cost. If the interface solve dominates, report that the original bottleneck moved rather than disappeared. Recursive interface elimination is a separate follow-up design, not an automatic expansion of this prototype. Do not substitute an approximate solve or alter GS ordering.

## 6. Stage D — connect to the FGS recurrence and certify accuracy

Only after Gates A–C, add an explicitly experimental dispatch, proposed name `sweep_order=:triinterface`, with `triinterface_parts` specifying the partition count. Validate positive counts not exceeding the number of leaves. Leave production defaults and the existing dagteam path unchanged. Keep construction failure visible; do not silently fall back and label the timing as the new method.

Integrate at the inner-sweep layer. Inspect these existing contracts before editing:

- `FastMultipole/src/solve_dagteam.jl`: plan construction, initialization, lower pulls, upper products/reduction, and `dagteam_inner_sweeps!`.
- `FastMultipole/src/solve.jl`: outer FMM/residual/relaxation boundaries and LU/plan lifetime.
- `FLOWPanel/src/FLOWPanel_solver.jl`: wrapper configuration and solver tests.

The new executor must initialize the upper contribution from the actual starting strengths, solve with the frozen upper contribution, rebuild the next upper contribution from the new strengths, and return the same RHS/residual bookkeeping expected by the outer iteration.

In particular, preserve the saved lower-product information used at the outer boundary. For the first correctness implementation, explicitly recompute the required final lower products using the existing matrices and completed strengths. Include this cost in end-to-end timing. Optimize that reconstruction only as a separately checked change; do not assume triangular-system residuals are numerically identical to the baseline's saved products.

Invalidate/rebuild the interface plan whenever its coefficients, leaf membership/order, partition, or precision changes. Reuse it for new RHS vectors only when the underlying operator remains valid under the existing cache contract.

Tests, in order:

1. One and three inner sweeps, zero and nonzero starts, compared with dagteam and lexicographic references in Float64.
2. Multiple outer iterations, relaxation, callbacks/history, and a transformed-solver fixture following the existing dagteam tests.
3. The existing `fgs_dagteam_gate1_test.jl` regression and FLOWPanel `test/runtests_unit_solver.jl`; run FGS history tests if that path is touched.
4. Float32 state/storage using the same converted operator as the matched dagteam mode. Use the independent BC evaluator for acceptance; internal residual agreement is insufficient.
5. Cold-start and representative warm-start solves at the existing acceptance target. Record iterations rather than forcing them to match when floating-point grouping changes.

**Gate D:** mathematical sweep equivalence on controlled fixtures, valid outer-solver bookkeeping, and independently accepted physical solves. Accuracy failures stop performance promotion.

## 7. Stage E — scaling and amortization study

Future official measurements follow the repository HPC and reproducibility policies (see `agent_policies/HPC.md`); this document authorizes no job submission, and local work uses no more than four threads.

Use matched dagteam+backoff and candidate runs at 1/8/16/32/64 threads, with BLAS=1, documented placement, identical geometry, precision, solver parameters, and acceptance target. Retain same-process alternating-arm pairs within a thread-count/placement configuration where the harness permits. Keep diagnostic runs separate from ranking runs and check instrumentation effects.

First compare on R4, then on the larger same-family meshes that survived Stage A. Freeze the selected partition per geometry before the final comparison; report any partition search cost separately. Report:

- Accepted solve time, sweep time, iteration count, and strong-scaling speedup relative to each method's own one-thread run.
- Absolute speedup against dagteam at the same thread count, not just improved scaling from a slower baseline.
- Setup time, extra memory, interface size/fill, and measured phase breakdown versus panel count.
- Cold and warm behavior, and the number of reused solves needed to amortize extra setup: extra setup time divided by positive per-solve saving. If there is no saving, there is no break-even.

The final verdict must answer Ryan's concern explicitly: **as leaf count grows, did useful concurrency grow, did the limiting interface grow, or did we merely accelerate each step of another long chain?** State the measured range; do not extrapolate beyond it as established scaling.

Promotion requires reproducible accepted-solve improvement exceeding observed run variation, affordable memory, and acceptable setup amortization. If only one geometry or size wins, report a conditional result instead of changing the general default.

## 8. Reading order and deliverables

Read this document, then the current dagteam implementation, then the [Stage 2 results](fgs_scalability_stage2_results_20260924.md), [lower-DAG sizing](fgs_lshortening_gate0_20260924.md), and [dagedge benchmark](fgs_dagedge_benchmark_results_20260924.md). The dagedge result is the warning against trusting graph-only speedup estimates.

Algorithm reference: [Torun et al., Partitioning and Reordering for SPIKE-Based Distributed-Memory Parallel Gauss-Seidel](https://user.ceng.metu.edu.tr/~manguoglu/PDFs/dmpGS_SISC.pdf), especially the triangular SPIKE derivation. Borrow the exact interface elimination idea; MPI and the paper's matrix reordering are outside this first experiment.

Deliverables are staged: (A) graph/fill feasibility report; (B) tiny algebra reference and tests; (C) standalone exact executor with measured costs; (D) integrated experimental sweep with accuracy evidence; (E) scaling, memory, and amortization verdict. Stop at the first failed gate and leave a concise explanation for the next agent.
