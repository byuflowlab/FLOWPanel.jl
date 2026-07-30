# `src/FLOWPanel_solver.jl` Summary

Solver backend definitions: everything that turns an assembled boundary-value
problem (a body's normals, control points, and boundary conditions) into panel
strengths. Covers direct (matrix-full) solvers, matrix-free Krylov solvers,
the FastGaussSeidel (FMM-native) solver, a diagonal flat-ground solver, and the
multi-body coupling schemes that sit on top of all of them.

> **Flagged for confirmation:** two items below (§Build/Solve Time
> Measurement, §Accuracy Improvements #1–2) describe behavior I verified by
> reading the code directly, not by running it — the `BackslashCoupled`
> kernel-offset omission and the `FGSSolver` return-value mismatch. If either
> is already known/intentional, ignore the flag; otherwise they look like real
> bugs worth a second pair of eyes before relying on those paths.

## The Underlying Physics

All solvers here discretize the same continuous problem: incompressible,
irrotational flow admits a scalar potential $\phi$ with $\nabla^2\phi = 0$
(Laplace's equation). Green's second identity turns this into a boundary
integral equation — the potential anywhere is fully determined by a
source-strength layer $\sigma$ and a doublet-strength layer $\mu$ smeared
over the body surface (`docs/src/potentialflow.md` derives this from first
principles; `docs/dirichlet_potential_theory.md` states FLOWPanel's exact
discrete form).

Panel discretization turns this continuous boundary integral into a dense
linear system

$$G \, q = \text{RHS}$$

where `_G!` (`src/FLOWPanel_solver.jl:137-196`) builds $G$ panel-pair by
panel-pair — $G_{ij}$ is the influence of source panel $j$ on target control
point $i$, either as induced potential $\phi$ (Dirichlet bodies, `DBC=true`)
or induced normal velocity $\mathbf u \cdot \hat{\mathbf n}$ (Neumann bodies,
`DBC=false`). $q$ is the unknown panel strength vector (source $\sigma$ for
non-lifting bodies, doublet $\mu$ for lifting bodies via the vortex-ring/
doublet-sheet formulation). Every solver in this file is a different strategy
for solving (or approximating a solve of) this same $Gq = \text{RHS}$ system
— they differ in whether $G$ is ever assembled, how exactly vs.
approximately the solve is done, and whether one body or several coupled
bodies are solved at once.

### `Backslash` — direct, per-body

**Physics:** assembles the *entire* dense $G$ once and solves exactly (up to
floating-point/LU conditioning) via Gaussian elimination. This is the
"ground truth" discretization of the boundary integral equation — no
convergence tolerance, no truncated far-field approximation.

**Workflow** (construct once, reuse across solves):

```julia
# src/FLOWPanel_solver.jl:221-231
function Backslash(body::AbstractBody{<:Any,<:Any,TF}) where TF
    G = zeros(TF, body.ncells, body.ncells)
    ...
    update_geometry = true
    _G!(G, body, body; kerneloffset=body.kerneloffset_panel, update_geometry=update_geometry)
    Glu = lu!(G)

    return Backslash{TF,typeof(Glu)}(G, Glu, rhs, Uext, phi_ext)
end
```

then per solve, `_solve!` (`:730-752` Neumann, `:754-772` Dirichlet) only
rebuilds/refactors $G$ when `update_G=true`; otherwise it reuses the cached
LU and just does the triangular solve:

```julia
# src/FLOWPanel_solver.jl:745-749 (Neumann _solve!)
    rhs = solver.rhs
    rhs .= zero(eltype(rhs))
    calc_bc_noflowthrough!(rhs, body.velocity, body.normals)

    ts = @elapsed ldiv!(view(body.strength, :, strength_index), solver.Glu, rhs)
```

Cost: $O(N^2)$ to assemble $G$, $O(N^3)$ to factor — the classic panel-method
scaling wall that motivates every other solver in this file.

### `BackslashCoupled` — direct, multi-body

**Physics:** identical idea, but for a `Tuple` of bodies solved as one block
system. $G$ becomes $(\sum N_i) \times (\sum N_i)$, assembled block-by-block
so that off-diagonal blocks capture body-to-body interference (e.g. wing/
ground-plane or wing/aileron coupling) directly inside the matrix, rather
than through outer iteration:

```julia
# src/FLOWPanel_solver.jl:1223-1235
        t_build = @elapsed begin
            for (bi, source) in enumerate(bodies)
                c = offsets[bi]+1 : offsets[bi+1] # columns of sources
                for (ti, target) in enumerate(bodies)
                    r = offsets[ti]+1 : offsets[ti+1] # rows of targets
                    _G!(view(solver.G, r, c), target,
                        source;
                        update_geometry=false)
                end
            end
        # Factorize G matrix and cache it in solver
            solver.Glu = lu!(solver.G)
        end
```

**Concern:** unlike `Backslash`'s constructor and `_solve!` (both above),
this `_G!` call passes **no explicit `kerneloffset`**. It silently falls back
to `_G!`'s bare default, `kerneloffset=source_system.kerneloffset` — the
generic field — rather than `kerneloffset_panel`, which is what every other
panel-panel $G$ build in this file uses. `Backslash` and `BackslashCoupled`
are supposed to be the same physics at different scales (single body vs.
coupled tuple); if `kerneloffset` and `kerneloffset_panel` are ever tuned to
different values, these two solvers will silently regularize near-singular
self/near-panel terms differently, producing systematically different
answers for what should be an equivalent discretization.

RHS assembly dispatches per-body by boundary-condition type:

```julia
# src/FLOWPanel_solver.jl:1073-1093
function boundary_condition!(body::AbstractBody{<:Any,<:Any,<:Any,true}, RHS, backend; optargs...)
    calc_bc_dirichlet(RHS, body, backend; optargs...)
end
function boundary_condition!(body::AbstractBody{<:Any,<:Any,<:Any,false}, RHS, backend; optargs...)
    calc_bc_noflowthrough!(RHS, body.velocity, body.normals)
end
```

`calc_bc_dirichlet` **accumulates** onto RHS (`RHS .-= body.potential`), so
`solver.rhs` must be zeroed every call — the code does this
(`fill!(solver.rhs, 0)` at `:1214`) and calls it out in a comment as not
"guaranteed on subsequent calls if the same solver object is reused."

### `KrylovSolver` — matrix-free GMRES, per-body

**Physics:** same $Gq=\text{RHS}$ system, but $G$ is never assembled or
stored. GMRES only needs the *action* of $G$ on a trial vector, and each such
matrix-vector product is computed by re-evaluating induced velocity/potential
via `influence!` — which can use the Fast Multipole Method for
near-linear-time evaluation instead of $O(N^2)$ direct summation. This is the
entire point of going matrix-free: avoid ever forming the $O(N^2)$ (storage)
/ $O(N^3)$ (factor) dense operator, at the cost of an iterative convergence
tolerance instead of an exact factorization.

The solver struct **is** the linear operator — it's callable:

```julia
# src/FLOWPanel_solver.jl:373-392 (Neumann operator application)
function (solver::KrylovSolver{<:AbstractBody{<:Any, <:Any, <:Any, false}})(C, B, α, β)
    ...
    solver.body.velocity .= 0
    influence!(solver.body, solver.body, solver.backend; velocity=true)
    # dot product with normals
    solver.body.velocity .*= solver.normals
    ...
    C .*= β
    C .+= α .* view(solver.body.velocity, 1, :)
end
```

wrapped as a `LinearOperators.LinearOperator` and handed to `Krylov.jl`:

```julia
# src/FLOWPanel_solver.jl:421-446
    prod! = (y, x, α, β) -> solver(y, x, α, β)
    tb = @elapsed begin
        A = LinearOperators.LinearOperator(TF2, nrows, ncols, symmetric, hermitian, prod!)
    end
    ...
    workspace = Krylov.krylov_workspace(Val(solver.method), A, RHS)
    x0_norm = LA.norm(solver.unabbreviated_strengths)  # TEMPORARY diagnostic
    Krylov.warm_start!(workspace, solver.unabbreviated_strengths)
```

**Physics note — warm-starting is unconditional, not optional.** Every call
seeds GMRES from `solver.unabbreviated_strengths` (the previous converged
solve). This is not a flag to toggle; it falls out for free from reusing the
same `KrylovSolver` object across repeated calls (e.g. successive outer
Gauss-Seidel sweeps, or successive AOAs in a sweep that reuses solvers). The
first-ever call is unaffected (starts from the constructor's all-zero
default).

### `KrylovCoupled` — matrix-free GMRES, multi-body

Same matrix-free idea as `KrylovSolver`, but the unknown vector is the
stacked strengths of *all* bodies at once (`_coupled_offsets` /
`_write_coupled_unknowns!` / `_collect_coupled_operator!` pack/unpack the
flat vector into/out of each body's `strength`/`potential`/`velocity`). Its
`solve!` (`:595-656`) follows the same zero-before-`influence!`,
zero-before-`boundary_condition!` accumulation pattern as `BackslashCoupled`
above, with the same "must re-zero on reuse" caveat called out in its own
comment (`:624-628`).

### `FGSSolver` — FastGaussSeidel, FMM-native

**Physics:** a block Gauss-Seidel fixed-point iteration built directly on top
of the FMM's own octree, rather than routed through the generic
`influence!`/backend abstraction the other solvers use. Panels are grouped
into tree leaves; each leaf's near-field ("self") interactions are
pre-factored as small dense blocks and solved exactly every sweep; far-field
("non-self") interactions between distant leaves are evaluated through
multipole-to-local (M2L) translations and treated as a fixed source term,
updated each outer sweep as other leaves' strengths change. This converges by
residual/MSE reduction, not a single exact solve — and unlike GMRES it isn't
guaranteed to converge (`FastMultipole.solve!` emits `@warn "FastGaussSeidel
did not converge..."` if it exhausts `max_iterations`).

```julia
# src/FLOWPanel_solver.jl:868-880
    FastMultipole.solve!(body, solver.fgs;
        max_iterations=solver.max_iterations,
        inner_iterations=solver.inner_iterations,
        tolerance=solver.tolerance,
        rlx=solver.rlx,
        scalar_potential=dirichlet_bc,
        gradient=!dirichlet_bc,
        hessian=false,
        reverse_pass=solver.reverse_pass,
        verbose=solver.verbose,
        final_update=false
    )
```

It additionally supports polynomial-extrapolation warm-starting across
repeated solves (e.g. across timesteps), via a rolling history buffer:

```julia
# src/FLOWPanel_solver.jl:804-817
@inline function project_solution!(body::AbstractBody, solver::FGSSolver)
    solver.project_solution || return false
    n = solver.solution_history_nsaved
    n < 2 && return false
    order = min(solver.project_solution_order, n - 1)
    H = solver.solution_history
    c0 = order + 1
    @views @. body.strength = c0 * H[:, :, 1]
    @inbounds for j in 1:order
        c = ifelse(isodd(j), -binomial(order + 1, j + 1), binomial(order + 1, j + 1))
        @views @. body.strength += c * H[:, :, j + 1]
    end
    return true
end
```

**Concern — return-value mismatch.** `_solve!(body, solver::FGSSolver; ...)`
(`:841-895`) ends with `save_solution!(body, solver)`, which itself
`return nothing`s (`:821-831`) — so `_solve!` returns `nothing`, not a
`(tb, ts)` tuple. There is no `FGSSolver`-specific outer `solve!(body,
solver)` method, so a caller going through the public API hits the generic
wrapper:

```julia
# src/FLOWPanel_solver.jl:234-256 (Dirichlet) / :258-276 (Neumann)
        tb, ts = _solve!(body, solver; backend, update_G, optargs...)
```

Destructuring `tb, ts = nothing` errors. As far as I can tell from this file
alone, `FGSSolver` is only ever driven through `_solve!` directly or through
`FastMultipole.solve!`, never through the public `solve!(body, solver)` path
— if that's correct, this is a latent bug rather than an active one, but it
means `FGSSolver` currently cannot participate in the outer coupled
Gauss-Seidel loop (`solve!(bodies::Tuple, solvers::Tuple)` below), which
calls exactly this generic per-body `solve!`.

### `FlatGroundSolver` — diagonal special case

**Physics:** for flat, *exactly coplanar* constant-source ground panels, a
source panel's induced velocity normal to its own panel plane, evaluated at
another point in that same plane, is exactly zero by the up/down
antisymmetry of a source panel's velocity field about its own plane. So every
off-diagonal entry of $G$ vanishes identically — this is an exact geometric
fact for coplanar panels, not an approximation — leaving only the self-term
(the panel's own one-sided limit, $-0.5$):

```julia
# src/FLOWPanel_solver.jl:1013-1023
function _solve!(body::NonLiftingBody{ConstantSource, 1, TF}, solver::FlatGroundSolver{TF}; optargs...) where TF
    calc_bc_noflowthrough!(solver.rhs, body.velocity, body.normals)
    for i in 1:body.ncells
        # For a flat ground problem with constant source panels, the influence of each panel on itself is -0.5
        # ...and there is no influence from other panels.
        body.strength[i, 1] = -solver.rhs[i] / (-0.5)
    end
    return nothing
end
```

No matrix is ever formed — each panel solves independently, $O(N)$.

**Concern:** this exactness depends entirely on the ground panels being
*exactly* coplanar. If a `FlatGround` body is ever translated/rotated
slightly out of plane (geometry generation error, or future code that
composes ground panels with a body's own kinematics), the off-diagonal terms
are no longer zero and this solver would silently return a wrong answer —
there's no check that the body is actually planar before dispatching here.

### Outer coupled Gauss-Seidel loop — `solve!(bodies::Tuple, solvers::Tuple)`

**Physics:** drives any *mix* of independently-typed per-body solvers
(`Backslash`, `KrylovSolver`, `FGSSolver`, ...) — one solver per body — as an
outer fixed-point iteration. This is a different coupling strategy than
`BackslashCoupled`/`KrylovCoupled`: instead of one big block system, each
body solves its *own* system every sweep, picking up the other bodies'
induced velocity as an RHS update:

```julia
# src/FLOWPanel_solver.jl:934-953
    for iter in 1:max_outer_iterations
        iter_update_G = iter == 1 && update_G
        for (i, (body, solver)) in enumerate(zip(bodies, solvers))
            body.velocity .= prev_velocity[i]
            sources = tuple((bodies[j] for j in eachindex(bodies) if j != i)...)
            t_influence = 0.0  # TEMPORARY diagnostic
            if !isempty(sources)
                t_influence = @elapsed influence!((body,), sources, backends[i];
                    scalar_potential=false, velocity=true, optargs...)
            end
            tb, ts = solve!(body, solver; backend=backends[i], update_G=iter_update_G)
            println("  [outer diag] iter=$iter body=$i t_influence=$t_influence tb=$tb ts=$ts")
            t_build += tb
            t_solve += ts
        end
        ...
        if max_delta < outer_tolerance
            converged = true
            break
        end
    end
```

Converges when the largest per-panel strength change across all bodies drops
below `outer_tolerance` (default `1e-6`), capped at `max_outer_iterations`
(default `50`). `update_G` only forces a rebuild once, on the first outer
iteration (`iter_update_G = iter == 1 && update_G`) — matrix-full solvers
whose cached $G$ has gone stale (geometry/wake changed since last build) need
this; matrix-free solvers ignore the flag.

## Build/Solve Time Measurement — How It's Actually Computed

Every `solve!`/`_solve!` returns `(tb, ts)` — "build" and "solve" wall time,
via `@elapsed`. `examples/benchmarks.jl` sums these into a
`solver,nps,AOA,CL,CD,t_build,t_solve,total_time` CSV row per AOA per solver,
intending to compare solver families head-to-head. **The comparison is not
apples-to-apples**, for two separate reasons found by reading each solver's
timed region directly:

**1. What counts as "build" differs by solver family in a way that doesn't
reflect actual cost.**

| Solver | `tb` measures | `ts` measures |
|---|---|---|
| `Backslash` | $G$ assembly ($O(N^2)$) + LU factor ($O(N^3)$) — real, often-dominant cost — *when `update_G=true`* | triangular solve only (cheap, $O(N^2)$) |
| `BackslashCoupled` | same, for the block system — *when `update_G=true`* | triangular solve only |
| `KrylovSolver`/`KrylovCoupled` | just wrapping a `LinearOperators.LinearOperator` (`:422-430`) — essentially free, no computation happens | **all** GMRES iterations, each performing a full `influence!` evaluation — this is where the real per-solve cost lives |
| `FGSSolver` | n/a — see return-value concern above | n/a |

So a `t_build` column showing near-zero for `KrylovSolver` doesn't mean
building is cheap for that solver — matrix-free solvers don't have a
separable "build" phase in the same sense; comparing `t_build` across rows
answers "how expensive was matrix *assembly*," which only exists for
matrix-full solvers, while comparing `total_time` (or `t_solve` alone) is the
only fair cross-family comparison.

**2. A more consequential gap: cross-body coupling evaluation is untimed
everywhere it appears**, so `t_build + t_solve` undercounts true wall time —
by an amount that varies by solver and body count, which further undermines
comparing `total_time` across rows:

- In the outer Gauss-Seidel loop (excerpt above), `t_influence` — the cost of
  each body picking up every *other* body's induced velocity — is computed,
  printed, but **never added to `t_build` or `t_solve`** (`:941-951`).
- In `BackslashCoupled`'s `solve!`, the coupled `influence!` call and the
  `boundary_condition!` RHS assembly both happen **before** `t_solve = 0.0`
  is even initialized (`:1200`, `:1215`, vs. `:1216`) — only the final
  `ldiv!` (and, conditionally, the $G$ rebuild) are timed.
- The single-body Dirichlet `solve!` wrapper has the same gap: its
  `influence!` self-assembly call (`:249`) happens before `_solve!` is
  invoked, so it's outside `(tb, ts)` too.

Cross-body/self `influence!` evaluation is frequently the most expensive part
of a coupled solve (it's exactly the operation FMM acceleration exists for),
so excluding it from the returned timing means **the reported `total_time`
can systematically understate wall-clock cost, and by different amounts for
different solver/body-count combinations** — which is precisely the kind of
number `benchmarks.jl` uses to argue one solver is faster than another.

## Possible Improvements

### Accuracy, first

1. ~~**Fix `BackslashCoupled`'s missing `kerneloffset_panel`**~~ **FIXED :
   2026-07-30 :** `_G!` now passes `kerneloffset=source.kerneloffset_panel`
   explicitly, matching `Backslash`. Regression test:
   `test/runtests_unit_solver.jl`, `"BackslashCoupled uses kerneloffset_panel,
   not stale kerneloffset"`.
2. ~~**Fix `FGSSolver`'s `_solve!` return value**~~ **FIXED : 2026-07-30 :**
   `_solve!` now times the `FastMultipole.solve!` call and returns `(tb, ts)`
   like every other `_solve!`. Verified via the existing `"FGSSolver
   boundary-condition preparation"` test, which exercises the public
   `solve!(body, solver)` path that previously crashed on `tb, ts = nothing`.
3. **Audit compounded iterative tolerances.** A `KrylovSolver`
   (`atol=rtol=1e-6` default) driven by the outer Gauss-Seidel loop
   (`outer_tolerance=1e-6` default) has two nested convergence tolerances
   whose combined effective accuracy isn't obviously $10^{-6}$ — worth
   characterizing empirically (e.g. residual vs. a `Backslash` reference)
   rather than assuming the tightest quoted tolerance applies end-to-end.
   Same concern applies to `FGSSolver`'s own `tolerance`/`max_iterations`
   when it's driven through the outer loop.
4. ~~**Reconcile `update_G` defaults across call paths**~~ **AUDITED + FIXED :
   2026-07-30 :** the four defaults are not equally arbitrary — audited
   every `solve!(body, ...)`/`solve!(bodies, ...)` call site in `src/` and
   `examples/`. `BackslashCoupled`'s `update_G=true` default and the outer
   Gauss-Seidel loop's `update_G=false` default are both *required*, not
   inconsistent: `BackslashCoupled`'s constructor only allocates a dummy
   identity-matrix `G` (`Glu = lu!(G)  # dummy init; will be overwritten on
   first update_G=true`), so its first solve genuinely needs `true`; the
   outer loop's per-body solvers are constructed the normal (eager-build)
   way, so trusting the pre-built `G` at `false` is correct. The one real,
   accidental inconsistency was single-body Dirichlet (`false`) vs.
   single-body Neumann (`true`) `solve!` — both wrap `Backslash`, whose
   constructor eagerly builds+factors `G` the same way regardless of
   Dirichlet/Neumann, so there was no construction-time reason for these to
   differ. No audited call site relies on the Neumann `true` default for
   correctness (every one either reconstructs its solver fresh right before
   the single solve, or passes `update_G` explicitly) — changed the Neumann
   default to `false` to match, eliminating a redundant `G` rebuild on every
   such solve. See the `# IMPROVEMENT` comment at the Neumann `solve!`
   wrapper and the new `"Backslash Neumann solve! defaults update_G=false"`
   test.
5. **Guard `FlatGroundSolver`'s coplanarity assumption** — even a cheap
   assertion that all panels share a plane/normal before dispatching to the
   diagonal solve would turn a silent wrong-answer failure mode into a loud
   one.
6. ~~**Verify `FGSSolver` warm-start extrapolation isn't compounding
   under-converged iterations.**~~ **CONFIRMED : 2026-07-30 (investigation,
   fix deferred) :** the risk is real, not hypothetical, confirmed by
   reading both sides of the call rather than by running it (Julia
   execution was blocked by an unrelated, concurrently-in-progress
   FastMultipole.jl branch issue at investigation time). `_solve!`
   (`FGSSolver` branch, this file) calls `FastMultipole.solve!(body,
   solver.fgs; ...)` and unconditionally calls `save_solution!(body,
   solver)` right after — there is no gate on convergence. On the
   FastMultipole side (`FastMultipole.jl/src/solve.jl:968-971`), the
   convergence check exists (`if mse > tolerance; @warn "FastGaussSeidel did
   not converge..."; end`) but `mse` is a local variable never returned —
   `solve!` ends on `buffer_to_system_strength!(...)` with no explicit
   `return`, so FLOWPanel has no programmatic signal to gate on, only the
   `@warn` text. **A real fix needs a coordinated, cross-package change**
   (`FastMultipole.solve!` would need to return `mse`/a converged flag, and
   `FGSSolver`'s `_solve!` would need to skip or flag `save_solution!` when
   not converged) — deferred rather than attempted unilaterally from the
   FLOWPanel.jl side alone, especially while FastMultipole.jl is under
   active, separate development.

### Speed, second

1. **Make cross-body `influence!` timed and FMM-accelerated by default.**
   `BackslashCoupled`, `KrylovCoupled`, and the outer Gauss-Seidel loop all
   default `backend=DirectBackend()` for the RHS/coupling evaluation
   (`O(N^2)`) while the rest of the pipeline defaults to FMM — for larger
   multi-body cases this untimed, unaccelerated step is a plausible hidden
   bottleneck. Defaulting it to `FastMultipoleBackend()` (or at least
   including its cost in `(t_build, t_solve)`) would both speed up large
   cases and make the benchmark numbers trustworthy.
   **DONE : 2026-07-30 :** the timing half is fixed — cross-body
   `influence!`/`boundary_condition!` cost is now folded into `(tb, ts)` for
   `BackslashCoupled`, `KrylovCoupled`, and the outer loop (previously
   discarded/only printed). For the default-backend half: turns out
   `KrylovCoupled` already defaulted its coupling evaluation to
   `FastMultipoleBackend()` (via its constructor's own `backend` field,
   threaded through `solve!(bodies::Tuple, solver::KrylovCoupled;
   backend=solver.backend, ...)`) — only `BackslashCoupled`'s `solve!` and
   the outer Gauss-Seidel loop's `solve!` actually defaulted to
   `DirectBackend()`. Benchmarked Direct vs. FMM coupling cost first (two
   spheres 5 units apart, 256/1024/2304 total panels, default FMM settings):
   FMM relative error vs. direct was ~1e-8 to 1e-16 (negligible), while
   Direct became ~2x then ~20x slower than FMM as panel count grew past
   ~1000; FMM was only slower (~1.5x) at the smallest count (~250 total,
   tree-build overhead not yet worth it). Audited every `solve!` call site in
   `src/`, `examples/`, and `test/` reaching these two methods — all pass
   `backend` explicitly, so none relied on the old Direct default. Changed
   both defaults to `FastMultipoleBackend()`; see the `# IMPROVEMENT`
   comments at each site and the two new
   `"... solve! defaults to FastMultipoleBackend"` regression tests.
2. **Reserve `Backslash`/`BackslashCoupled` for small $N$ or as a validation
   baseline.** Their $O(N^2)$ assembly / $O(N^3)$ factorization is the
   scaling wall the other solvers exist to avoid; for large panel counts
   `FGSSolver` or `KrylovSolver`+FMM should be the default recommendation,
   with the direct solvers kept around explicitly for correctness checks
   (which is how `benchmarks.jl`'s commented-out sections already seem to
   be using them).
3. **Benchmark the opt-in Jacobi preconditioner.** `KrylovSolver`'s
   block-Jacobi preconditioner (`preconditioner_cell_size`, `:321,:330-335`)
   is disabled by default — enabling it trades a one-time build cost (itself
   currently untimed, same gap as above) for potentially far fewer GMRES
   iterations. Worth quantifying rather than leaving off by default.
   **UPDATE : 2026-07-30 :** `KrylovCoupled` previously had no preconditioner
   option at all; it now takes the same `preconditioner_cell_size` kwarg as
   `KrylovSolver` (mirrored field/constructor/`M=...` wiring), so both Krylov
   solver families can be compared with/without preconditioning on equal
   footing — this is exactly the comparison the in-progress Krylov
   accuracy/timing study (`examples/benchmarks.jl`) is set up to quantify.
   **CAVEAT found while testing this :** `JacobiPreconditioner`'s default
   `derivatives_switches=(scalar_potential=true, gradient=true)`
   (`FastMultipole.jl/src/solve.jl:1088`) builds each cell's self-influence
   matrix consistent with a Dirichlet (potential) formulation; applying it to
   a **Neumann** (velocity-dot-normal) coupled system produces badly wrong,
   non-converged strengths (confirmed empirically: two source-only
   `NonLiftingBody` targets, preconditioned vs. unpreconditioned answers
   differed by ~100%, tangency residual ~0.3-0.7 instead of <1e-6). This
   looks like a FastMultipole-side gap in `JacobiPreconditioner` for Neumann
   targets, not a `KrylovCoupled` wiring bug — out of scope to fix from
   `FLOWPanel.jl` alone. Not a blocker for the wing+aileron Krylov study,
   since that case uses the Dirichlet (`Union{ConstantSource,ConstantDoublet}`)
   formulation, which is the configuration already confirmed to work. The
   `"KrylovCoupled + JacobiPreconditioner"` test uses Dirichlet bodies
   accordingly (mirroring `"KrylovSolver + JacobiPreconditioner"`'s existing,
   proven-good configuration).
4. **Consider accelerating the outer Gauss-Seidel fixed point.** The loop
   iterates on `max_delta` with no extrapolation (Aitken/Anderson-style)
   toward the fixed point — for slowly-converging multi-body configurations,
   a cheap acceleration step could cut `max_outer_iterations` needed, which
   directly reduces the (currently untimed but real) repeated
   `influence!` cost per body per iteration.
5. ~~**Check whether `FGSSolver`'s default backend is rebuilt per call.**~~
   **RESOLVED : 2026-07-30 (investigation, no code change) :** confirmed not
   a bug. `FastMultipoleBackend` (`src/FLOWPanel_fmm.jl:36-40`) is a plain
   3-field config struct (`expansion_order`, `multipole_acceptance`,
   `leaf_size`) with no tree or other precomputed state — constructing a
   fresh one every `_solve!` call is a cheap 3-field allocation, and it does
   *not* force any extra tree reconstruction beyond what already happens:
   the actual octree is built downstream inside `influence!`/
   `FastMultipole.solve!` from the body's current positions on every call
   regardless of whether the `FastMultipoleBackend` object itself is new or
   reused. No fix needed.
6. ~~**Remove/guard unconditional debug output in hot paths**~~
   **FIXED : 2026-07-30 :** removed `println("Tuple of bodies")` (was
   `:922`), `println("BackslashCoupled")` (was `:1176`), and the two
   `[gmres diag/...]`/`[outer diag]` per-call prints (§Known Rough Edges);
   see the `# IMPROVEMENT` comments at each site in `src/FLOWPanel_solver.jl`.

## Known Rough Edges

- Several `_solve!`/outer-loop diagnostic `println`s are explicitly marked
  "TEMPORARY" in comments (warm-start investigation) and not yet removed —
  see §Speed Improvements #6.
- Two large commented-out `solve!` methods near the end of the file
  (`RigidWakeBody{TK,1,TF}` direct solve, `RigidWakeBody{Union{VortexRing,
  UniformVortexSheet},3,TF}` prescribed-circulation solve) are dead code from
  an earlier API, kept for reference — not part of the current dispatch
  surface.
- `update_G` defaults differently by call path — see §Accuracy Improvements
  #4 for the consolidated list and the risk it poses.
