# Local GS plus coarse correction: principles and promising designs

> **Superseded by the [strategy doc](fgs_overlapping_gs_coarse_correction_strategy_20260926.md)** — kept for background; the QR diagnostic and error-propagator material now live there and in the implementation handoff.

Date: 2026-09-26. Conceptual explanation and recommendations only; no solver implementation or tests performed.

Related: [overall scaling plan](fgs_cooperative_leaf_scaling_plan_20260926.md), [exact triangular-solver plan](fgs_exact_triangular_solver_plan_20260926.md). Target: up to 64 shared-memory threads. No Krylov outer iteration.

## 1. What the coarse solve is supposed to accomplish

The aim is to let many subdomains perform local GS concurrently while using a small, globally coupled solve to correct the error those local iterations remove slowly.

For intuition, suppose each region has mostly the right detailed strength distribution, but the loading amplitude varies incorrectly from region to region. Local solves can improve each region conditional on the neighboring values, yet repeatedly exchange an unresolved large-scale imbalance. A coarse solve represents a handful of changes to each region and solves for all their amplitudes together. A global loading error can then be corrected in one coarse step instead of many local-update cycles.

For a boundary integral operator, interactions are already global: it would be inaccurate to claim that information travels only one neighboring subdomain per iteration. The issue is that the inexpensive local approximate inverse may handle certain coupled error patterns poorly. The coarse solve treats those patterns globally.

This differs from the exact interface method. The exact method retains enough information to reproduce a GS sweep. A coarse space deliberately retains only selected error patterns, changes the stationary iteration, and relies on local smoothing to handle the rest. That freedom can make its global solve much smaller.

## 2. The correction in equations

Write the actual linear system for one fixed solve as

$$
Ax=b.
$$

Here $A$ includes the full near- and far-field panel operator, applicable bound-wake coupling, and any formulation constraints/transformation. Contributions that are prescribed for this solve belong in $b$. Let $e=x_\star-x$ be the error and $r=b-Ax=Ae$ the residual.

Choose a matrix $P$ with a small number $m$ of columns. Each column is a possible correction pattern in the fine-grid unknowns. A coarse vector $c$ produces the full correction $Pc$. Choose a restriction $R$ that maps the fine residual to $m$ coarse equations, and construct

$$
A_c=RAP.
$$

The coarse step is

$$
A_c c=Rr,\qquad x\leftarrow x+Pc.
$$

Assuming $A_c$ is nonsingular and solved accurately, this step eliminates any error lying exactly in the range of $P$. Indeed, if $e=Pa$, then $Rr=RAPa=A_ca$, so $c=a$.

Thus the main question is **whether $P$ represents the errors that survive local GS**, rather than whether the coarse mesh resembles a low-resolution picture of the body. This is the approximation property needed from the coarse space.

The coarse operator must also be stable: representing the right patterns is insufficient if the chosen residual restriction makes the coarse solve nearly singular or amplifies complementary error.

## 3. A complete stationary cycle

An initial design would use a multiplicative two-level cycle:

1. Apply a fixed small number of parallel local GS correction sweeps, starting with one.
2. Compute the full residual at the updated strengths.
3. Restrict that residual, solve the small coarse system, and prolong/add its correction.
4. Optionally apply one parallel local GS post-sweep to remove local error left or introduced by the correction.
5. Repeat, checking the full problem's residual and independent BC acceptance.

Subdomains execute concurrently; within each subdomain the existing leaf ordering and cached diagonal LU solves can be reused. This need not introduce dense exact subdomain factorizations. A local sweep can use the existing near-field blocks as an approximate local inverse acting on a full residual.

Use a common residual snapshot for parallel subdomain corrections. For overlapping subdomains, combine corrections with a partition of unity or restricted ownership so duplicated unknowns are not updated multiple times at full weight. A small overlap halo is worth comparing with disjoint regions because the partition cuts strong near-field interactions. Define the halo from the near-field coupling graph, not solely surface adjacency.

Let $E_s$ and $E_p$ denote the error propagators of pre- and post-smoothing. The full cycle has error propagator

$$
E_{\mathrm{cycle}}=E_p\left(I-P(RAP)^{-1}RA\right)E_s.
$$

For a fixed linear stationary cycle, asymptotic convergence requires $\rho(E_{\mathrm{cycle}})<1$. Good scalability asks for a contraction factor that does not deteriorate substantially as the number of subdomains or mesh resolution increases. For a nonsymmetric, nonnormal system, also measure transient growth; the asymptotic spectral radius alone does not characterize practical short solves.

Keep the basis and smoother fixed during each solve. Building a better basis in setup does not require a Krylov outer solver. The repeated two-level cycle itself is the solver.

The established two-level Schwarz/multigrid framework combines these complementary components. Its convergence theorems for particular elliptic problems are motivation, not a theorem for FLOWPanel's panel formulation. See [Ciaramella and Vanzan, Substructured Two-Grid and Multigrid Domain Decomposition Methods (Numer. Algorithms)](https://link.springer.com/article/10.1007/s11075-022-01268-0) — a different paper from their spectral-coarse-spaces one (J. Sci. Comput., cited in the strategy doc); both are valid — and the [PETSc multigrid interface](https://petsc.org/release/manualpages/PC/PCMG/).

## 4. The most promising coarse spaces

### First: geometric aggregation, followed by measured enrichment

Start with surface patches made from groups of existing leaves. Distinguish the **compute partition** from the **coarse approximation partition**. Sixty-four worker subdomains do not require exactly 64 coarse unknowns. Several coarse patches can sit inside each compute subdomain, and overlap may span subdomain boundaries.

A simple initial basis has one strength-offset mode per connected patch, with appropriate separation of unknown types and bodies. This gives the coarse solve independent amplitudes with which to redistribute loading. Constants are candidates, not an assertion that this operator has a constant nullspace. Existing gauge or compatibility constraints must still be enforced.

Then add variation within patches: local surface-coordinate linear modes, or smaller patches along directions where the error varies. On a blade, spanwise/chordwise coordinates may be more appropriate than arbitrary Cartesian axes. Remove dependent modes on nearly degenerate patches. Preserve distinctions between separate blades, components, sharp features, and distinct surface sheets instead of merging them merely because their points are nearby.

For scale intuition, 64 compute subdomains with four constant coarse patches each gives 256 coarse unknowns. That is an illustration, not a chosen optimum. Adding two independent linear modes per patch would raise it to 768, so basis enrichment must be judged against construction, storage, and solve costs.

Piecewise constants are useful as the simplest diagnostic, but their discontinuities at aggregate boundaries can be a poor approximation to the surviving error. Overlapping/blended basis functions or carefully filtered aggregate modes are candidates for improvement. Do not blindly apply an SPD smoothed-aggregation recipe to this nonsymmetric operator; measure whether the filtering is stable and retains rank and useful modes.

### Second: learn the error that local GS fails to remove

This is the most promising route to making the method robust rather than relying only on geometric intuition.

For setup diagnostics, choose several independent initial vectors and apply the proposed local stationary iteration to the homogeneous problem $Ae=0$. The iterates themselves are known error vectors. Use a modest fixed number of steps and retain snapshots; the components that remain or grow reveal weaknesses of the smoother. Respect any nullspace constraints throughout.

Use those snapshots to enrich the patch basis. Restrict them to patches, remove redundant directions with rank-revealing QR/SVD, and retain a few important modes on each patch. Check the enriched method on independent initial errors and physical RHS vectors. A single long relaxed vector may reveal only one dominant mode and is not sufficient coverage.

The homogeneous experiment avoids confusing residual patterns with error patterns: a small residual component does not necessarily mean small solution error. For difficult nonnormal behavior, inspect transient snapshots as well as the eventual surviving vector.

This follows the principle of adaptive algebraic multigrid: identify components that relaxation does not eliminate and teach the coarse space to represent them. See [Adaptive Algebraic Multigrid Methods](https://www.osti.gov/servlets/purl/875948). The cited analysis focuses on SPD systems; applying the principle here requires the nonsymmetric stability checks below.

My strongest candidate is therefore **local GS with modest overlap, a geometric aggregate basis, and a small amount of enrichment from observed slow error**. The first prototype should begin without enrichment so we can measure what each addition fixes.

### Third: an actual coarse panel mesh

A geometrically coarsened panel model could directly represent large-scale loading with familiar panel unknowns. It may be attractive if compatible coarse meshes and transfers already exist, and can reduce setup cost compared with many fine-operator applications.

The difficulty is consistency: coarsening must preserve trailing edges, shedding relationships, narrow gaps, body topology, and strength-transfer meaning. A coarse model assembled independently generally differs from $RAP$. If it misses a critical coupling, the exact coarse-space elimination identity no longer holds.

For that reason I would first establish convergence using aggregation on the existing fine discretization and $A_c=RAP$. A rediscretized coarse panel model is a later cost optimization or a separately validated alternative.

### Why not use the existing FMM multipoles directly?

The FMM tree is useful for finding spatial groups, but its expansion coefficients summarize how sources generate fields. A solver coarse space must describe strength errors and how to correct them while satisfying boundary conditions. Those are different roles. FMM moments can motivate candidate modes, but the FMM hierarchy is not automatically a multigrid solver hierarchy.

## 5. The coarse operator must contain global physics

For the preferred full-system design, use $A_c=RAP$ with the **full** operator. This makes every coarse loading pattern interact with all other patches, including distant bodies or blades. That global coupling is the reason the coarse correction might recover the convergence lost by splitting local GS.

A matrix-free construction is straightforward conceptually: apply $A$ to each coarse basis vector and restrict the result. Dense storage and a reusable direct factorization of $A_c$ can be reasonable when $m$ is small. Use a fixed linear operator application, or an accurate enough approximation with its error checked; strength-dependent evaluation choices must not silently invalidate the linear identities.

Do not insert this coarse correction inside an inner block with a stale far-field contribution and then claim it corrects the current full-system residual. Start at a point where a consistent full residual can be formed after the local sweep. A coarse correction built only for the frozen near-field subproblem is a different experiment and should be labeled as such.

A useful bookkeeping identity is

$$
r_{\mathrm{new}}=r-(AP)c.
$$

Caching $AP$ can avoid an extra full operator application immediately after coarse correction, but costs $O(Nm)$ storage. Alternatively, discard $AP$ after constructing $A_c$ and recompute the full residual when needed. Count both options honestly; the FMM cost cannot be treated as free.

## 6. Nonsymmetry: the main reason to test rather than assume

FLOWPanel's effective collocation/transformed systems cannot generally be treated as SPD. Start with an appropriately scaled restriction such as $R=P^T W$, where $W$ is a positive residual weighting compatible with the formulation. Panel areas can be a candidate for suitable scalar collocation equations, but are not a universal choice for every transformed equation set.

With nonsymmetric $A$, $P^T WAP$ need not be positive definite or even well conditioned. A coarse correction can eliminate represented errors while amplifying other components. Monitor coarse conditioning and complete-cycle contraction. Selecting damping for the smoother or coarse step can help in some cases, but damping is not a general cure for an inadequate coarse space or unstable projection. See [Southworth and Manteuffel, On Compatible Transfer Operators in Nonsymmetric Algebraic Multigrid](https://epubs.siam.org/doi/10.1137/23M1586069).

A useful diagnostic alternative for a fixed basis is residual-minimizing coarse correction:

$$
c=\operatorname*{arg\,min}_{z}\|W^{1/2}(r-APz)\|_2.
$$

Use a QR factorization of $W^{1/2}AP$, not normal equations. This coarse step cannot increase that residual norm in exact arithmetic and still removes error in the coarse range when the relevant columns are independent. It does not guarantee convergence of the whole smoother-plus-coarse cycle, and its $N\times m$ storage/setup cost can be substantial. It is a fixed coarse-space projection, not a Krylov outer iteration.

I would use it to distinguish a poor basis from an unstable restriction, rather than immediately commit to its production cost. If the basis works only with more careful restriction, a Petrov–Galerkin design informed by left as well as right slow modes is a subsequent option.

## 7. What scalable behavior would mean

For fixed mesh size, increasing the number of local subdomains should shorten local work without causing a large rise in cycle count. For increasing mesh size, the basis must keep representing the error left by the smoother. A fixed 64-dimensional coarse space is not guaranteed to do that as local subdomains become larger or geometry becomes more complicated.

The convergence objective is an approximately stable error-reduction factor per cycle; the performance objective is low time to independently accepted accuracy. Neither follows merely from high worker occupancy.

The coarse solve is a new global operation. Dense coarse storage scales as $O(m^2)$, its factorization as $O(m^3)$, and its repeated solve as $O(m^2)$. If the required coarse dimension grows too large, the next step is a multilevel hierarchy so the intermediate coarse problem is itself treated by local work plus a smaller coarse correction. That is a later extension, not necessary for the first 64-thread experiment.

Operator-dependent setup can be expensive: constructing a basis with hundreds of columns through separate full-operator applications can cost hundreds of matvecs. Reuse helps across RHS vectors only while the operator and transfers remain valid. Moving geometry, changing relative body placement, or changing wake coupling can require rebuilding. Report setup amortization separately from per-solve speed.

## 8. How I would discriminate among the designs

Start with a small representative assembled operator so the numerical idea can be tested without a new threaded runtime. Compare local-only iteration, local plus patch constants, local plus linear/enriched patch modes, and the existing global GS reference.

Perform three especially informative tests:

1. **Coarse-space exactness:** initialize an error $e=Pc$. One undamped, consistently constructed coarse correction should remove it to numerical accuracy. Failure indicates an operator/transfer/constraint implementation problem.
2. **Complementarity:** apply local smoothing to independent random and physically structured errors, then measure how well the coarse basis represents the survivors and how much the actual coarse step removes. Separate a basis-approximation failure from a projection-stability failure.
3. **Partition robustness:** increase subdomain count while holding the problem fixed, then increase mesh size. Measure complete-cycle contraction, transient growth, accepted iteration count, full operator applications, time, and memory.

Ablate overlap and coarse enrichment separately. Include cold and warm starts and more than one RHS. Compare all methods at the same independent BC accuracy target; an internally small residual in a rounded or approximate operator is insufficient.

The most encouraging result would be that local-only iteration degrades as subdomain count rises, while a modest enriched coarse space keeps the two-level cycle count approximately stable. If cycle count still grows, examine the surviving error before adding more local sweeps: those sweeps may be spending time on components already well controlled.

My recommended sequence is: **geometric aggregate correction with the full operator; modest overlap; enrichment from measured surviving error; only then more sophisticated restriction or a multilevel extension.** This is more likely to reveal a useful coarse space than assuming that one global mode per worker is enough, and it directly tests whether coarse correction can repair the convergence penalty observed in the earlier chunked-GS experiment.
