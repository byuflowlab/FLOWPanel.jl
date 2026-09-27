# 021: overlapping local GS with coarse correction — principles and experiment order

Date: 2026-09-26. Status: proposed research phase; documentation only. No implementation, tests, or campaign have been performed for this document.

This consolidates the [local GS/coarse-correction explanation](fgs_local_gs_coarse_correction_explanation_20260926.md) and the subsequent discussion of duplicated interface leaves, halos, reconciliation, sensitivity propagation, and local explicit inverses. It expands the stationary-iteration alternative in the [overall FGS scaling plan](fgs_cooperative_leaf_scaling_plan_20260926.md). The [exact triangular-solver plan](fgs_exact_triangular_solver_plan_20260926.md) is a separate approach with a different numerical contract.

**Target:** faster time to independently accepted accuracy on up to 64 shared-memory threads, including larger panel problems. No distributed memory and no GMRES/Krylov outer iteration.

Implementation handoff: [Step-by-step prototype construction, tests, and decision gates](fgs_overlapping_gs_coarse_implementation_handoff_20260926.md). Use that document to carry out the bounded first implementation; use this document for the governing principles and broader alternatives.

## 1. Objective and governing principles

The proposed method exchanges the long global GS chain for many concurrent local GS solves, then repairs the resulting convergence penalty through overlap and a small globally coupled correction.

The intended division of work is:

- **Local GS** removes error that can be treated cheaply inside each subdomain.
- **Interface halos** preserve strong interactions cut by the partition and reduce artificial boundary effects.
- **Coarse correction** coordinates the error that remains difficult for the parallel local iteration, including global loading patterns and persistent interface-related error.

These are hypotheses to test, not guaranteed properties of this panel operator. The method changes the GS iterates; it must retain the original linear problem and acceptance target.

The following principles govern the entire phase:

1. **Optimize accepted time to solution, not worker occupancy or sweep time alone.** Earlier chunked GS reduced sweep work but increased iterations from 27 to 44. More parallel activity is not sufficient evidence of progress.
2. **Choose interfaces from coupling, not consecutive leaf indices.** A first/last leaf pair is only a toy model. A real leaf can have many near-field neighbors in other subdomains. Current dagteam scheduling is dynamic, so persistent subdomains must be introduced explicitly; a thread itself has no fixed physical interface.
3. **Overlap and coarse correction have complementary roles.** Overlap can repair local cut interactions, but need not remove a large-scale convergence bottleneck. A coarse basis must represent error that actually survives the smoother.
4. **Agreement between duplicated strengths is not accuracy.** Averaging can force copies to agree even when their common value is wrong. Use disagreements to diagnose and construct modes; use the full equation residual to drive and certify correction.
5. **Keep global physics in the coarse equations.** For the preferred design, the coarse operator and residual use the full panel operator, including near/far contributions and the applicable bound-wake/formulation terms.
6. **Treat nonsymmetry explicitly.** Geometric smoothness, SPD convergence bounds, and stability of a transpose-based restriction cannot be assumed. Measure contraction, transient growth, coarse conditioning, and physical accuracy.
7. **Compress the useful responses, not necessarily the inverse.** A few interface-to-interior response modes may supply most of the benefit at much lower cost than a full subdomain inverse or exact interface solve.
8. **Bound each experiment.** Establish which component fails before adding overlap, modes, damping, or levels. A negative feasibility result is a useful outcome.

The recorded dagteam+backoff baseline and the failed fine-edge decomposition remain essential controls: [Stage 2](fgs_scalability_stage2_results_20260924.md), [dagedge results](fgs_dagedge_benchmark_results_20260924.md). Their results do not establish performance for this new iteration.

## 2. What should be computed

### Parallel local correction with halos

For one fixed linear solve, write $Ax=b$ and $r=b-Ax$. Here $A$ is the actual operator in the solver's coordinates, with required constraints/gauges respected. Prescribed contributions belong in $b$.

Partition the core leaves into balanced persistent subdomains. Initially preserve the existing order within each subdomain. Add a halo around leaves incident to cut near-field edges. Include relevant couplings in both directions; do not treat the lower GS predecessor chain as the entire physical interaction graph.

All subdomains start from a common global iterate/residual snapshot. Each computes a correction over its core plus halo, using a fixed small number of local GS sweeps and existing leaf factors. Combine the corrections with weights that sum to one on each global unknown, optionally followed by damping. An owner-only publication rule is a comparison arm: halo values assist the solve, but only the core owner publishes each unknown.

Weighted averaging of absolute solutions is equivalent to averaging corrections only when they share the same starting global iterate and consistent weights. Do not reconcile stale, differently based updates as if they were equivalent.

For a strong pair or small cluster straddling a partition, an augmented local matrix must include **both the unknowns and their equations**, including all internal cross blocks. Remove those internal interactions from any frozen exterior RHS to avoid double counting. A joint solve is a possible refinement; do not begin by densely factoring entire large subdomains.

### Full-residual coarse correction

After the local update, compute a consistent full residual. Let columns of $P$ be fine-grid correction patterns and let $R$ restrict residuals to coarse equations. Form and factor

$$
A_c=RAP.
$$

Apply

$$
A_c c=Rr,\qquad x\leftarrow x+Pc.
$$

The coarse unknowns may describe interface patterns, but **the coarse RHS is derived from the entire residual**, not just the disagreement or residual on duplicated leaves. For example, $R=P^T W$ with a suitable positive equation weighting includes interior residuals wherever the extended patterns have support. This is an initial restriction to test, not a universally stable choice.

If the current error lies in the range of $P$, an undamped consistent coarse correction removes it exactly in exact arithmetic, provided $A_c$ is nonsingular. The practical objective is to approximate the errors left by local GS well enough that the cycle count does not rise sharply with partition count.

A genuinely condensed interface solve would have RHS $r_\Gamma-A_{\Gamma I}A_{II}^{-1}r_I$. This also contains interior residual information. Discarding $r_I$ is not equivalent unless the interior equations are already satisfied. Exact condensation of the full dense operator is outside this phase.

### Interface-response basis

For a local interior/interface split with exterior values held fixed, an interface perturbation induces

$$
\delta x_I=-A_{II}^{-1}A_{I\Gamma}\delta z.
$$

Choose a few interface patterns as columns of $Q$. Compute only their interior responses:

$$
A_{II}H_Q=-A_{I\Gamma}Q.
$$

Each combined interface pattern and interior response becomes a candidate column of $P$. Assemble overlapping columns into global vectors consistently, then remove dependencies. Use a few modes per **connected interface patch**, not one mode per interface panel or per individual cut edge. Patches may bundle many leaves along the boundary between subdomains.

An exact local solve on a tiny fixture supplies the reference. For practical subdomains, fixed local GS response sweeps using the cached near-field operator can supply approximate extensions. Label them approximate; do not claim that they exactly balance the full dense interior equations. Their usefulness is determined by the resulting full-system coarse correction.

Retain a small set of bulk loading modes as well. The dense integral operator can have difficult global error with little duplicate disagreement, and the sparse-PDE assumption that all residual is confined to interfaces does not apply automatically.

### Whole stationary cycle

Start with one local sweep, a full-residual coarse correction, and zero or one local post-sweep. Recompute or consistently update the residual after every change needed by the next stage. Do not use a stale frozen far field while claiming a full-system coarse correction.

With fixed smoother and transfers, the cycle is stationary and requires no Krylov outer solver. Keep basis adaptation in setup between solves during this phase. The error propagator is

$$
E_{\mathrm{cycle}}=E_{\mathrm{post}}\left(I-P(RAP)^{-1}RA\right)E_{\mathrm{pre}}.
$$

Asymptotic convergence requires $\rho(E_{\mathrm{cycle}})<1$, but practical short solves also depend on transient growth and initial error. Always measure complete trajectories.

## 3. Attempt the paths in this order

### Gate −1. Operator-application budget (before any smoother code)

The proposed cycle needs 2–3 full operator applications per iteration versus ~1 per outer iteration today, so application cost bounds the method from below. On the target fixture:

1. **Measure** the cost of one complete operator application as this method defines it. The ~37 ms/outer figure inferred from FMM time (1.0 s FMM / 27 outer iterations, [gate-0 analysis](fgs_lshortening_gate0_20260924.md)) is FMM work only, not necessarily a full residual evaluation.
2. **Project** total accepted-solve cost as (applications/cycle) × measured application cost × cycles + local work + coarse work + setup amortization.
3. **Compare** against the measured baselines — dagteam+backoff 3.24 s cold at R4 j64, and the warm-start ILU/FGS numbers ([warm-start R4 results](fgs_warmstart_r4_results_20260925.md)) for the production regime — not against hypothetical cooperative-team projections.

Fewer than 27 cycles is not required to beat 3.24 s; what matters is projected total cost under realistic cycle counts. Illustrative floor for intuition only: 2 applications × 37 ms × 27 cycles ≈ 2.0 s of application cost alone. A go/no-go margin (e.g. "projected ≥15% under baseline") is a chosen experiment threshold, not a model consequence — state it as such when picking one. This gate makes §4's "count all operator applications" numeric.

### A0. Cheap first diagnostic: coarse correction on the existing unmodified dagteam sweep

Before any halos, partition-of-unity weights, or threaded local runtime, bolt the coarse step onto the **outer** iteration with the current dagteam sweep kept as the smoother: after the FMM residual, solve $A_c c=P^T r$, update $x\leftarrow x+Pc$, continue. This probes the central unknown — whether a small basis captures this operator's slow error — at minimal cost. Scope it honestly:

- The handoff's operator-adapter correctness gate (Gate 0) applies here too: after the coarse update, strengths, source buffers, and cached near/far contributions must be made consistent before the next sweep. This integration/state handling is real work, not "zero smoother risk".
- Success is measured in total cost, not iteration count alone: cutting 27→~15 outer iterations is a win only after the extra residual evaluations and coarse-solve/prolongation work per iteration are counted (Gate −1 accounting applies).
- Failure does not kill the partitioned branch: global GS leaves different surviving error than partitioned local GS, so a null result here rules out only this configuration. On failure, record a diagnosis — basis coverage of the surviving error, projection stability, and cost — and let that diagnosis steer (or stop) the partitioned-smoother stages rather than declaring the branch dead.

Attribution rationale: the near-field-only local smoother in the handoff changes two things at once — partitioning AND dropped cross-subdomain near-field coupling (the chunked-GS 27→44 result is direct evidence the latter hurts) — so failures of the full design are hard to diagnose without this decoupled step.

### A. Establish a small numerical reference and failure diagnosis

Use a small representative assembled operator first, preserving the actual formulation, constraints, and nonsymmetry. Compare global GS, disjoint local GS, and the proposed corrections before creating a new threaded runtime.

Use multiple independent homogeneous-error starts and representative physical RHS vectors. A homogeneous test applies the iteration to $Ae=0$, so the evolving vector is known error rather than a residual being mistaken for error. Preserve the constrained solution space. Include cold and nonzero starts.

Confirm local correction signs, indexing, residual bookkeeping, and zero-residual invariance. Construct a trial $e=Pc$ and verify that one consistent undamped coarse step removes it to numerical accuracy. If this identity fails, fix the implementation before tuning modes.

**Decision:** identify whether local partitioning leaves coherent error that a small basis might represent. Failure of local-only convergence does not by itself reject the two-level method, but rapid amplification must be understood. Do not proceed based only on a favorable physical RHS.

### B. Add a modest halo and establish a reconciliation baseline

Compare disjoint subdomains with one near-field graph layer of halo. Try two layers only if the first still leaves error concentrated at the cut interactions. Record overlap inflation, $\sum_s n_s^{\mathrm{extended}}/N$, and work imbalance before building expensive data.

Start with partition-of-unity averaging; compare owner-only publication if duplicate reconciliation is costly or convergence is poor. Keep local sweep count fixed initially. If required, use a bounded damping comparison of 1, 1/2, and 1/4; freeze the selected value during each solve. Do not endlessly tune damping to conceal failure.

**Decision:** retain overlap only when its convergence benefit compensates for duplicated work or enables a better two-level cycle. If larger halos approach whole-domain duplication, destroy useful concurrency, or offer negligible improvement, stop increasing halo depth. Weak global contraction after local boundary error is removed is a reason to add coarse correction, not more halo.

### C. Add the simplest full-system coarse correction as a control

Use connected geometric patch constants, plus separate modes where bodies, unknown types, or constraints require them. Do not assume constants are physical null modes. Coarse patches need not coincide one-to-one with compute subdomains.

Build $A_c=RAP$ with the full operator and a suitable scaled restriction. Start in Float64 with a direct coarse factorization. Compare local+halo alone against the same smoother plus this coarse step. Test one post-sweep separately rather than changing several components together.

**Decision:** if this simple basis stabilizes cycle count across partition counts, preserve it as a serious candidate; greater sophistication must beat it in total cost. If it fails, distinguish poor basis coverage from an unstable coarse projection before moving on.

### D. Add interface-response modes — preferred refinement

Start with constant and independent linear variations on each connected interface patch, extending them into adjacent subdomains through the response solve above. Add them to, or compare them against, the patch/bulk control basis under a matched coarse-dimension budget. Rank-reveal and remove duplicates.

For each proposed mode, ask whether it represents error that survived the halo smoother. If simple interface shapes miss it, retain a few persistent disagreement/error snapshots, extract independent patch modes, and extend them. Test on independent snapshots and RHS vectors so the basis is not merely fitted to one solve. Raw disagreement alone is insufficient training data.

**Decision:** continue enriching only when each increase materially improves complete-cycle contraction or accepted time. Stop this path if the required mode count grows toward the full interface size, the coarse problem becomes dominant, or larger bases mostly fit already-damped error. That outcome says a compact interface representation has not been found; do not silently turn the method into an exact interface solver.

### E. Diagnose projection stability only when indicated

If the basis approximates surviving errors well but the actual coarse step amplifies them or has a poorly conditioned coarse matrix, investigate restriction/scaling. As a diagnostic, compare a fixed-space residual-minimizing correction:

$$
c=\operatorname*{arg\,min}_z\|W^{1/2}(r-APz)\|_2.
$$

Use QR of $W^{1/2}AP$, not normal equations. This cannot increase that residual norm in exact arithmetic, but it does not guarantee convergence of the full cycle. It may require costly $N\times m$ storage. It is a fixed coarse-space correction, not a Krylov outer method.

**Decision:** if QR correction succeeds where the original restriction fails, the basis merits further restriction/Petrov–Galerkin work. If neither works, return to the smoother/basis diagnosis. Do not keep expanding a poorly targeted basis to compensate for projection instability.

### F. Optimize local application only after convergence is credible

Consider explicit inverses for small leaf or augmented-interface matrices. An inverse replaces triangular substitutions with independent output-row products, potentially helping local throughput. It does not remove inter-subdomain coupling, improve the mathematical local solve beyond LU, or by itself create a coarse correction.

For already-dense small blocks, inverse storage is comparable to LU storage, but construction costs and numerical behavior differ. Compare repeated application time and residual quality. For larger subdomains prefer cached $H_Q$ responses, computed by multiple-RHS LU solves or fixed local sweeps, to full inverse storage. Never construct the full panel-system inverse for this experiment.

**Decision:** discard inverse-based application if setup amortization, precision, or measured throughput is worse. This is a separable kernel optimization; its failure does not reject halo/coarse correction.

### G. Integrate and assess actual scaling

Only after the numerical gates, integrate an opt-in experimental stationary cycle. Preserve the existing champion as a control, independent BC evaluation, warm-start semantics, and operator/cache invalidation. Compare matched Float64 first, then the tuned Float32 mode with renewed certification. Production defaults remain unchanged.

Evaluate the retained candidates at increasing partition counts independently of thread count, then at 1/8/16/32/64 threads on HPC. Compare R4 and larger meshes from the same geometry/refinement family. Include setup, cold/warm starts, multiple RHS vectors, and memory; do not change the accuracy target to obtain a speedup.

Future official campaigns follow the repository HPC and reproducibility policies (see `agent_policies/HPC.md`); local runs use at most four threads, and this document does not authorize implementation or job submission.

**Decision:** consider another level only if a useful two-level method demonstrably loses scalability because coarse dimension/cost grows. A recursive hierarchy is a separate design step. Do not add levels to a method whose coarse space has not yet demonstrated that it removes the relevant error.

## 4. Measurement and stopping rules for the entire phase

Use three distinct judgments; none substitutes for the others:

| Question | Evidence | Stop or redirect when |
|---|---|---|
| Is the method correct? | Full residual bookkeeping, coarse exactness check, constraints, independent BC acceptance | Any check fails; numerical tuning cannot compensate for an inconsistent operator or RHS. |
| Is the decomposition numerically effective? | Error/residual histories, surviving-error coverage, transient growth, cycle count versus partition count and mesh size | Halo, a simple coarse basis, and one targeted enrichment/projection diagnosis fail to produce useful contraction. Record the failing modes before stopping. |
| Is it faster at useful scale? | Paired uninstrumented accepted solves, setup amortization, phase timings, memory | Extra cycles, FMM applications, overlap, or coarse work erase the local parallel savings. Better thread utilization alone does not justify continuing. |

For a rough ranking of convergent candidates, compare measured cycle cost $t_c$ and residual contraction $q$ through $t_c/(-\log q)$. Use a stable portion of the history and report transients. This estimates time per logarithmic residual reduction; it is not a substitute for accepted solve timing. When $q\geq1$, it provides no convergent rate.

An even simpler break-even condition is

$$
\frac{n_{\mathrm{candidate}}}{n_{\mathrm{baseline}}}
<\frac{t_{\mathrm{baseline\ cycle}}}{t_{\mathrm{candidate\ cycle}}},
$$

provided each cycle time includes all required work. With unequal cycle definitions, compare total solves directly. Count all full operator applications: moving work into residual evaluation does not make it disappear. Gate −1 in §3 turns this into a numeric budget gate against the measured baselines.

Constructing $A_c$ through $m$ separate operator applications can be expensive. Caching $AP$ allows $r\leftarrow r-APc$ after correction but costs $O(Nm)$ storage. Include coarse factorization, transfer storage, local response construction, and additional factors in the memory/setup accounting. Reuse is valid only while the operator, constraints, partition, and basis remain valid. Moving geometry or changing coupling may erase the expected amortization.

Dense coarse storage, factorization, and application cost $O(m^2)$, $O(m^3)$, and $O(m^2)$ respectively. Do not assume a fixed coarse dimension will remain sufficient as mesh size increases. Conversely, do not assume mode count must track every interface unknown; discovering a compact useful representation is the point of this phase.

Declare a winner only when paired accepted-solve savings exceed observed variation and persist on the intended sizes/starts. If only one case wins, report a conditional operating region. If no candidate wins after the bounded numerical comparisons, stop and document whether the obstacle was local instability, inadequate coarse coverage, projection instability, or cost. Do not automatically switch to Krylov, full inverses, or exact interface elimination.

## 5. Recommended initial bet and supporting references

**Initial bet:** one-layer coupling-graph halos, weighted local GS corrections, a simple full-operator patch coarse space as the control, followed by a small interface-response enrichment if needed. Use the full residual throughout; retain bulk modes for errors invisible to duplicate disagreement. Optimize inverses or response application only after the cycle has demonstrated useful convergence.

The expected benefit is a shorter local dependency chain without the earlier large iteration penalty. It remains an untested hypothesis for this solver. Neither overlap alone nor a small coarse solve alone guarantees scalability.

**Comparative prior (decided by measurement, not asserted):** coarse correction is the preferred research bet for this operator. With ~45 predecessors per leaf, spatial partitions cut edges everywhere, suggesting the exact triangular-solver plan's interface will be large with heavy fill; the coarse space sidesteps that density by representing only selected global patterns. But predecessor count does not determine separator size or fill, so the exact plan's inexpensive graph-only Stage A should still run and settle whether the interface is actually prohibitive. Amortization also fits this branch: setup ($m$ operator applications for $A_c$, enrichment) amortizes over many warm-started solves in a time-march when the operator is reusable across steps.

- [FreeFEM domain-decomposition documentation](https://doc.freefem.org/documentation/ffddm/introduction-to-the-domain-decomposition-method.html): overlapping corrections and partition-of-unity/ownership constructions.
- [Ciaramella and Vanzan, Spectral Coarse Spaces for the Substructured Parallel Schwarz Method](https://link.springer.com/article/10.1007/s10915-022-01840-9) (J. Sci. Comput.): interface-based coarse corrections and stationary Schwarz methods. Its PDE conclusions must not be asserted as results for this dense panel operator. Their companion paper on [substructured two-grid/multigrid domain decomposition methods](https://link.springer.com/article/10.1007/s11075-022-01268-0) (Numer. Algorithms) is a distinct, also-relevant reference.
- Hackbusch's **multigrid of the second kind** (multigrid for second-kind Fredholm integral equations; *Multi-Grid Methods and Applications*, ch. 16, and *Integral Equations: Theory and Numerical Treatment*) is relevant precedent alongside Caspar — framed here as a **hypothesis**, not established structure: dominant self-influence blocks alone do not establish second-kind character or explain convergence for this particular (transformed, constrained) operator; whether the surviving error is coarse-representable is exactly what these experiments must measure.
- [Southworth and Manteuffel, On Compatible Transfer Operators in Nonsymmetric Algebraic Multigrid](https://epubs.siam.org/doi/10.1137/23M1586069): why restriction/projection stability matters in addition to coarse-space coverage.
- [Druinsky and Toledo, How Accurate is inv(A)*b?](https://arxiv.org/abs/1201.6035): explicit inverses should be assessed by actual accuracy and cost rather than rejected categorically.

Required handoff from each experiment: configuration, the isolated hypothesis tested, correctness/accuracy status, contraction and cost evidence, and a clear continue/stop decision with its reason. No production implementation should be started merely because the conceptual decomposition looks parallel.
