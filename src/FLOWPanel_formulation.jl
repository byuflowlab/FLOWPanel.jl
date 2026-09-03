#=##############################################################################
# DESCRIPTION
    Swappable wake→body solve formulations for the Dirichlet finite-wake solve.

    The 052b production formulation (`VelocityThroughSources`) transfers free-wake
    influence to the body through velocity only, converted to source strengths.
    `docs/wake_solve_schemes.md` derives two alternatives implemented here:

    - `GreenReconstruction`: reconstruct the wake body-trace q by solving the
      body-only Green system (I−B)q = Sσ from sampled wake velocities, then
      solve the explicit-potential system G·μE = −S·σ0 − q.
    - `HybridWakePotential`: evaluate retained panel-wake potential directly,
      Green-reconstruct only the particle contribution, and preserve
      target-specific potential/velocity coupling in a multibody solve.
    - `TraceCorrected`: solve G·μ̃ = −S(σ0+σ) + W·c and apply the affine Kutta
      relation γ = C·μ̃ − c downstream through the body's wake-strength
      correction channel (see `set_wake_correction!`).

    Contrary to a remark in an earlier revision of the theory note, (I−B) is
    not a low-rank update of the factored G (only B = G − WC is), so the Green
    solve carries its own one-time bordered LU.
=###############################################################################

################################################################################
# FORMULATION TYPES (immutable option holders; runtime state lives separately)
################################################################################

abstract type AbstractSolveFormulation end

"""
    VelocityThroughSources(; max_outer_iterations=50, outer_tolerance=1e-8,
                           dirichlet_residual_scale=1,
                           neumann_residual_scale=1,
                           require_outer_convergence=false)

Production 052b wake→body solve formulation: the free wake acts on the body through
velocity only, converted to source strengths inside `solve!`. The outer-solve
options mirror [`HybridWakePotential`](@ref), allowing production baseline
runs to hard-fail on a nonconverged multibody solve. It has no runtime state.
"""
struct VelocityThroughSources <: AbstractSolveFormulation
    max_outer_iterations::Int
    outer_tolerance::Float64
    dirichlet_residual_scale::Float64
    neumann_residual_scale::Float64
    require_outer_convergence::Bool

    function VelocityThroughSources(; max_outer_iterations::Int=50,
            outer_tolerance::Real=1e-8,
            dirichlet_residual_scale::Real=1,
            neumann_residual_scale::Real=1,
            require_outer_convergence::Bool=false)
        max_outer_iterations >= 1 || error("max_outer_iterations must be >= 1")
        isfinite(outer_tolerance) && outer_tolerance > 0 || error(
            "outer_tolerance must be finite and positive")
        isfinite(dirichlet_residual_scale) && dirichlet_residual_scale > 0 || error(
            "dirichlet_residual_scale must be finite and positive")
        isfinite(neumann_residual_scale) && neumann_residual_scale > 0 || error(
            "neumann_residual_scale must be finite and positive")
        new(max_outer_iterations, Float64(outer_tolerance),
            Float64(dirichlet_residual_scale), Float64(neumann_residual_scale),
            require_outer_convergence)
    end
end

"""
    GreenReconstruction(; gauge=:area_mean, recompute_interval=1,
                        green_solver=nothing)

Explicit-potential formulation with a Green-reconstructed wake trace: solve the
body-only Green system `(I−B)q = Sσ` for the wake potential trace `q` (σ from
sampled wake velocities), then solve `G·μE = −S·σ0 − q`. The near-null constant
mode of `(I−B)` is removed by `gauge`:

- `:area_mean`: area-weighted zero-mean trace (bordered Lagrange row).
- `:lsq`: min-norm least squares with the same area constraint row (QR;
  requires `green_solver=nothing`).

`green_solver` selects how the Green system is solved:

- `nothing` (default): one-time dense assembly + bordered LU — the
  `Backslash`-analogue.
- a [`KrylovSolver`](@ref): matrix-free bordered solve honoring that solver's
  `method`, `itmax`, `atol`, `rtol`, and `backend` (Direct or FMM); a
  preconditioner, if set, is ignored (it targets `G`, not `(I−B)`).
- an [`FGSSolver`](@ref): matrix-free relaxed Picard iteration
  `q ← (1−rlx)·q + rlx·(Sσ + B·q)` with per-sweep area-mean projection,
  honoring `rlx`, `tolerance`, and `max_iterations`; the solver's internal
  FastGaussSeidel object and FMM tuning fields are unused, so it may be
  constructed with `build_fgs=false` when used only here. Warm-started from
  the previous step's `q`; products use the step's solve backend.

`recompute_interval` re-solves for `q` every that many steps, reusing the
cached trace in between (lagging on the wake-convection timescale).
"""
struct GreenReconstruction{TS<:Union{Nothing,KrylovSolver,FGSSolver}} <: AbstractSolveFormulation
    gauge::Symbol
    recompute_interval::Int
    green_solver::TS

    function GreenReconstruction(; gauge::Symbol=:area_mean,
            recompute_interval::Int=1,
            green_solver::Union{Nothing,KrylovSolver,FGSSolver}=nothing)
        gauge in (:area_mean, :lsq) ||
            error("Invalid GreenReconstruction gauge :$gauge; "*
                  "expected :area_mean or :lsq.")
        gauge === :lsq && !isnothing(green_solver) &&
            error("gauge=:lsq requires green_solver=nothing (dense QR route).")
        recompute_interval >= 1 ||
            error("recompute_interval must be >= 1; got $recompute_interval.")
        return new{typeof(green_solver)}(gauge, recompute_interval, green_solver)
    end
end

"""
    HybridWakePotential(; gauge=:area_mean, recompute_interval=1,
                        max_outer_iterations=50, outer_tolerance=1e-8,
                        dirichlet_residual_scale=1,
                        neumann_residual_scale=1,
                        require_outer_convergence=false)

Experimental 052e multibody explicit-potential formulation for mixed
panel/particle wakes. It is a comparison/regression API, not a 052b production
acceptance dependency.
Retained `PanelWake` potential is evaluated directly. The particle contribution
is isolated from the sampled wake velocity and projected independently onto
each Dirichlet body's harmonic boundary trace through the gauge-fixed Green
system `(I-B)q=Sσ`. During the coupled solve, Dirichlet sources are formed from
freestream/kinematics only; other bodies and both wake traces enter through
scalar potential. Neumann bodies (for example a paneled ground) retain the
complete wake velocity.

The dense body-local Green factorization is intentionally independent for each
body. Shared rigid-body factors are not enabled until callers have separately
validated exact operator parity.
"""
struct HybridWakePotential <: AbstractSolveFormulation
    gauge::Symbol
    recompute_interval::Int
    max_outer_iterations::Int
    outer_tolerance::Float64
    dirichlet_residual_scale::Float64
    neumann_residual_scale::Float64
    require_outer_convergence::Bool

    function HybridWakePotential(; gauge::Symbol=:area_mean,
            recompute_interval::Int=1, max_outer_iterations::Int=50,
            outer_tolerance::Real=1e-8,
            dirichlet_residual_scale::Real=1,
            neumann_residual_scale::Real=1,
            require_outer_convergence::Bool=false)
        gauge in (:area_mean, :lsq) || error(
            "Invalid HybridWakePotential gauge :$gauge; expected :area_mean or :lsq.")
        recompute_interval >= 1 || error("recompute_interval must be >= 1")
        max_outer_iterations >= 1 || error("max_outer_iterations must be >= 1")
        isfinite(outer_tolerance) && outer_tolerance > 0 || error(
            "outer_tolerance must be finite and positive")
        isfinite(dirichlet_residual_scale) && dirichlet_residual_scale > 0 || error(
            "dirichlet_residual_scale must be finite and positive")
        isfinite(neumann_residual_scale) && neumann_residual_scale > 0 || error(
            "neumann_residual_scale must be finite and positive")
        new(gauge, recompute_interval, max_outer_iterations,
            Float64(outer_tolerance), Float64(dirichlet_residual_scale),
            Float64(neumann_residual_scale), require_outer_convergence)
    end
end

"""
    TraceCorrected(; estimator=:green, gauge=:area_mean, quadrature=:trapezoid,
                   interior_path=:straight, path_depth=1.0,
                   recompute_interval=1, relaxation=1.0)

!!! warning "Deprecated (2026-07-29)"
    Use [`GreenReconstruction`](@ref) instead. Deprecated for three reasons:
    (1) Gate A0 of the DJI convergence campaign
    (`logs/dji_convergence_20260722/phase_02b_formulation_attribution.md`)
    found the whole affine-`c` channel moves integrated circulation by only
    ~0.74% and rejected it as the Dirichlet–Neumann lever; (2) its `:green`
    estimator equals `GreenReconstruction` to displayed precision (the Stage 6
    consistency identity), making it redundant; (3) its downstream velocity
    channel is vacuous — `influence!` never applies `wake_strength_shift`
    under any backend (including `DirectBackend`, which also routes through
    the FastMultipole source buffers), so post-solve surface velocities and
    Bernoulli pressures/forces silently omit the affine −c attached-strip
    term. The shed circulation (`shed_wake!` via `_get_wakestrength_mu`) and
    the solve itself are unaffected. Kept functional for archived scripts
    (`examples/dji9443_formulation_attribution.jl`, `test/formulation_test.jl`)
    and old-run replay; construction warns once. If a corrected-velocity
    variant is ever needed, `_add_affine_attached_velocity!`
    (`src/FLOWPanel_kutta.jl`) is the ready-made exact fix.

Velocity-only trace-corrected formulation: solve the production system with one
added right-hand-side vector, `G·μ̃ = −S(σ0+σ) + W·c`, and apply the affine
Kutta relation `γ = C·μ̃ − c` downstream (shedding, influence evaluation, VTK,
circulation output) through the body's wake-strength correction channel.

`c` estimators:

- `:green`: Kutta trace of the Green-reconstructed body trace, `c = C·q` with
  `q` from the same `(I−B)q = Sσ` solve as [`GreenReconstruction`](@ref)
  (gauge selected by `gauge`; solver route selected by `green_solver` — see
  [`GreenReconstruction`](@ref) for the `nothing`/`KrylovSolver`/`FGSSolver`
  options). Robust for particle and rolled-up wakes.
- `:line_integral`: interior line integral of wake-only velocity from the
  lower to the upper trailing-edge control point. `quadrature` is `:trapezoid`
  (endpoint reuse, no extra influence evaluations) or `:simpson` (one batched
  wake evaluation at edge midpoints). `interior_path=:deformed` routes the
  path through an interior point `path_depth`·(TE chord length) upstream of
  the trailing edge to avoid the near-singular straight chord.

`recompute_interval` recomputes `c` every that many steps; `relaxation` blends
`c ← (1−α)c_prev + α c_new` when recomputing (α = 1 disables).
"""
struct TraceCorrected{TS<:Union{Nothing,KrylovSolver,FGSSolver}} <: AbstractSolveFormulation
    estimator::Symbol
    gauge::Symbol
    quadrature::Symbol
    interior_path::Symbol
    path_depth::Float64
    recompute_interval::Int
    relaxation::Float64
    green_solver::TS

    function TraceCorrected(; estimator::Symbol=:green,
            gauge::Symbol=:area_mean,
            quadrature::Symbol=:trapezoid,
            interior_path::Symbol=:straight,
            path_depth::Real=1.0,
            recompute_interval::Int=1,
            relaxation::Real=1.0,
            green_solver::Union{Nothing,KrylovSolver,FGSSolver}=nothing)
        estimator in (:green, :line_integral) ||
            error("Invalid TraceCorrected estimator :$estimator; "*
                  "expected :green or :line_integral.")
        gauge in (:area_mean, :lsq) ||
            error("Invalid TraceCorrected gauge :$gauge; "*
                  "expected :area_mean or :lsq.")
        quadrature in (:trapezoid, :simpson) ||
            error("Invalid TraceCorrected quadrature :$quadrature; "*
                  "expected :trapezoid or :simpson.")
        interior_path in (:straight, :deformed) ||
            error("Invalid TraceCorrected interior_path :$interior_path; "*
                  "expected :straight or :deformed.")
        path_depth > 0 || error("path_depth must be > 0; got $path_depth.")
        recompute_interval >= 1 ||
            error("recompute_interval must be >= 1; got $recompute_interval.")
        0 < relaxation <= 1 ||
            error("relaxation must be in (0, 1]; got $relaxation.")
        gauge === :lsq && !isnothing(green_solver) &&
            error("gauge=:lsq requires green_solver=nothing (dense QR route).")
        return new{typeof(green_solver)}(estimator, gauge, quadrature,
                   interior_path, Float64(path_depth), recompute_interval,
                   Float64(relaxation), green_solver)
    end
end

"""
    DirectWakePotential(; recompute_interval=1)

Direct fixed-wake scalar-potential formulation — the production version of the
Task-3 diagnostic (`debug/dirichlet_solve/task3.md`). Each step solves the
explicit-potential system

    G_Δ · μE = −S·σ0 − q_f ,      σ0 = −u_nonwake · n

on the marched finite wake, where `q_f` is the scalar potential evaluated
*directly* from the constant-doublet `PanelWake` panels at exact body centroids
(via the FMM `direct!`/scalar-potential path) — NOT a reconstructed trace as in
[`GreenReconstruction`](@ref). The wake velocity is not converted into body
sources; no mean is removed from `q_f`. The finite-body `G_Δ` factorization is
the [`Backslash`](@ref) solver's own LU, applied manually (single-body `solve!`
clears preassembled external potential).

Requires a single Dirichlet source+doublet `RigidWakeBody` with
`semiinfinite_wake=false`, a `Backslash` body solver, and a finite `PanelWake`
with `include_final_filament=false` — every active wake source must expose a
scalar potential, so particle and mixed wakes are rejected. Runs on either the
direct or FMM backend (no `DirectBackend` guard: only a scalar potential is
needed, which the FMM `PanelWake` path already supplies).

`recompute_interval` re-evaluates `q_f` every that many steps, reusing the
cached potential in between; `recompute_interval=1` (recompute every step) is
required for production convergence runs.
"""
struct DirectWakePotential <: AbstractSolveFormulation
    recompute_interval::Int

    function DirectWakePotential(; recompute_interval::Int=1)
        recompute_interval >= 1 ||
            error("recompute_interval must be >= 1; got $recompute_interval.")
        return new(recompute_interval)
    end
end

################################################################################
# RUNTIME STATE (concrete parametric types, built once by
# `initialize_formulation` where body/solver types are known)
################################################################################

"Flattened shedding-edge index map: edge e has upper/lower panel indices
(`lower[e] <= 0` marks an unpaired edge), plus the owning shedding surface and
in-surface edge index for `Das` lookups."
struct SheddingEdgeMap
    upper::Vector{Int}
    lower::Vector{Int}
    isurf::Vector{Int}
    iedge::Vector{Int}
end

function _shedding_edge_map(body::RigidWakeBody)
    upper = Int[]; lower = Int[]; isurf = Int[]; iedge = Int[]
    for (s, shed) in enumerate(body.shedding)
        for i in axes(shed, 2)
            push!(upper, shed[1, i])
            push!(lower, shed[4, i])
            push!(isurf, s)
            push!(iedge, i)
        end
    end
    return SheddingEdgeMap(upper, lower, isurf, iedge)
end

"Kutta map application: `out[e] = x[upper_e] − x[lower_e]` (0 contribution from
an unpaired lower panel)."
function _apply_kutta_map!(out::AbstractVector, x::AbstractVector,
        edges::SheddingEdgeMap)
    for e in eachindex(edges.upper)
        lo = edges.lower[e]
        out[e] = x[edges.upper[e]] - (lo > 0 ? x[lo] : zero(eltype(x)))
    end
    return out
end

"Common supertype of the Green-system `(I−B)q = Sσ` solve routes (dense
factorization, matrix-free Krylov, relaxed Picard). Every subtype carries the
reconstructed trace in a `q::Vector` field."
abstract type AbstractGreenState end

"Dense gauge-fixed body-only Green solve (the `Backslash`-analogue): one-time
assembly of B plus a bordered LU (`gauge=:area_mean`) or QR (`gauge=:lsq`)."
struct GreenSolveState{TF, TFACT} <: AbstractGreenState
    fact::TFACT           # :area_mean → LU of bordered [(I−B) a; aᵀ 0]
                          # :lsq       → QR of [(I−B); aᵀ]
    gauge::Symbol
    q::Vector{TF}         # length N
    rhs_b::Vector{TF}     # length N+1
    sol_b::Vector{TF}     # length N+1 (used by :area_mean)
end

"Matrix-free bordered Green solve reusing a `KrylovSolver`'s options and
backend. The (N+1)-sized bordered operator makes the solution identical to the
dense `:area_mean` route up to the Krylov tolerance. The operator closure and
the scratch window both act on `body` (captured at initialization), so the
state is self-consistent regardless of which body instance reaches the solve."
struct GreenKrylovState{TF, TB<:RigidWakeBody, TKS<:KrylovSolver, TOP, TW} <: AbstractGreenState
    body::TB              # body the operator closure evaluates on
    ks::TKS               # honored: method, itmax, atol, rtol, backend
    a_hat::Vector{TF}     # unit-norm panel-area vector (gauge row/column)
    q::Vector{TF}         # length N
    rhs_b::Vector{TF}     # length N+1
    sol_b::Vector{TF}     # length N+1: previous [q; λ] (warm start)
    Bx::Vector{TF}        # product work buffer
    op::TOP               # bordered LinearOperator
    workspace::TW         # Krylov workspace, allocated once
    have_warm::Base.RefValue{Bool}
end

"Matrix-free relaxed Picard iteration `q ← (1−rlx)·q + rlx·(Sσ + B·q)` with
per-sweep area-mean projection, reusing an `FGSSolver`'s `rlx`, `tolerance`,
and `max_iterations` (its FastGaussSeidel object and FMM tuning fields are
unused). Warm-started from the previous step's `q`."
struct GreenPicardState{TF, TFGS<:FGSSolver} <: AbstractGreenState
    fgs::TFGS
    a::Vector{TF}         # raw panel areas
    suma::TF
    q::Vector{TF}
    q_new::Vector{TF}
    Bx::Vector{TF}
end

struct GreenReconstructionState{TF, TG<:AbstractGreenState}
    green::TG
    u_prewake::Matrix{TF} # 3×N snapshot of body.velocity before the wake pass
    sigma::Vector{TF}     # wake-only source coefficients σ = −u_f·n
    sigma0::Vector{TF}    # non-wake source coefficients σ0 = −u_0·n
    Ssigma::Vector{TF}    # source-potential work buffer
    rhs::Vector{TF}
    last_recompute::Base.RefValue{Int}
end

struct HybridBodyState{TF,TG<:AbstractGreenState}
    green::TG
    u_prewake::Matrix{TF}
    panel_velocity::Matrix{TF}
    particle_velocity::Matrix{TF}
    sigma_particle::Vector{TF}
    Ssigma::Vector{TF}
    q_panel::Vector{TF}
    q_total::Vector{TF}
    q_hodge::Vector{TF}
    green_residual::Base.RefValue{TF}
    gauge_defect::Base.RefValue{TF}
    green_hodge_mismatch::Base.RefValue{TF}
    tangential_projection_defect::Base.RefValue{TF}
end

struct HybridWakePotentialState
    bodies::Vector{Any}
    body_indices::Vector{Int}
    last_recompute::Base.RefValue{Int}
end

struct TraceCorrectedState{TF, TG<:Union{AbstractGreenState, Nothing},
                           TP}
    edges::SheddingEdgeMap
    W::Matrix{TF}         # N×M unit-strength attached-panel potential columns
    u_prewake::Matrix{TF}
    sigma::Vector{TF}     # wake-only σ (used by the :green estimator)
    green::TG
    probes::TP            # FastMultipole.ProbeSystem or nothing
    n_probes_per_edge::Int
    c::Vector{TF}
    c_new::Vector{TF}
    work::Vector{TF}      # length N (holds W·c)
    last_recompute::Base.RefValue{Int}
end

"Runtime state for [`DirectWakePotential`](@ref): reusable per-step buffers for
the direct fixed-wake potential solve. The finite-body `G_Δ` LU is the
`Backslash` solver's own factorization (reused, not held here)."
struct DirectWakePotentialState{TF}
    u_prewake::Matrix{TF} # 3×N snapshot of body.velocity before the wake pass
    sigma0::Vector{TF}    # non-wake source coefficients σ0 = −u_nonwake·n
    sigma::Vector{TF}     # wake-induced σ (computed for completeness, unused in RHS)
    q_wake::Vector{TF}    # direct wake scalar potential q_f at body centroids
    Ssigma::Vector{TF}    # source-potential work buffer (S·σ0)
    rhs::Vector{TF}       # linear-system RHS  −S·σ0 − q_f
    last_recompute::Base.RefValue{Int}
end

################################################################################
# VALIDATION AND INITIALIZATION
################################################################################

_single_body_solver(body_solvers) =
    body_solvers isa Tuple ? body_solvers[1] : body_solvers

function _validate_formulation_common(f::AbstractSolveFormulation,
        systems_tuple::Tuple, wakes_tuple::Tuple, body_solvers)
    fname = nameof(typeof(f))
    length(systems_tuple) == 1 ||
        error("$fname currently supports a single body; "*
              "got $(length(systems_tuple)) systems.")
    body = systems_tuple[1]
    body isa RigidWakeBody{<:Any, 2, <:Any, true} ||
        error("$fname requires a Dirichlet (DBC=true) source+doublet (NK=2) "*
              "RigidWakeBody; got $(typeof(body)).")
    body.semiinfinite_wake &&
        error("$fname requires semiinfinite_wake=false "*
              "(finite attached wake handing off to a marched free wake).")
    solver = _single_body_solver(body_solvers)
    solver isa Backslash ||
        error("$fname currently requires a Backslash solver; "*
              "got $(typeof(solver)).")
    any(!isnothing, wakes_tuple) ||
        error("$fname requires a wake system; got none.")
    if any(any(view(shed, 4, :) .<= 0) for shed in body.shedding)
        @warn "$fname: body has unpaired shedding edges; their Kutta trace "*
              "correction is forced to zero (an explicit gauge convention)." maxlog=1
    end
    return body, solver
end

"""
    initialize_formulation(formulation, systems_tuple, wakes_tuple,
                           body_solvers, backend_solve, backend_system)

Validate the formulation against the simulation setup and build its runtime
state (`nothing` for [`VelocityThroughSources`](@ref)). Called once by
`simulate!` before the time loop; the state carries all one-time assembled
operators, which assumes rigidly invariant body+`Das` geometry (the same
assumption `Backslash` reuse already makes).
"""
initialize_formulation(::VelocityThroughSources, systems_tuple, wakes_tuple,
    body_solvers, backend_solve, backend_system) = nothing

function initialize_formulation(f::HybridWakePotential, systems_tuple,
        wakes_tuple, body_solvers, backend_solve, backend_system)
    any(!isnothing, wakes_tuple) || error(
        "HybridWakePotential requires at least one wake system")
    solvers_tuple = body_solvers isa Tuple ? body_solvers : (body_solvers,)
    length(solvers_tuple) == length(systems_tuple) || error(
        "HybridWakePotential requires one solver per body")
    states = Any[]
    indices = Int[]
    for (i, (body, solver)) in enumerate(zip(systems_tuple, solvers_tuple))
        has_dirichlet_bc(body) || continue
        body isa RigidWakeBody{<:Any,2,<:Any,true} || error(
            "HybridWakePotential Dirichlet target $i must be a source+doublet RigidWakeBody")
        solver isa Backslash || error(
            "HybridWakePotential Dirichlet target $i requires a Backslash solver")
        TF = eltype(body.strength)
        N = body.ncells
        green = _build_green_state(body, f.gauge, nothing)
        push!(states, HybridBodyState{TF,typeof(green)}(green,
            zeros(TF, 3, N), zeros(TF, 3, N), zeros(TF, 3, N),
            zeros(TF, N), zeros(TF, N), zeros(TF, N), zeros(TF, N),
            zeros(TF, N), Ref(TF(Inf)), Ref(TF(Inf)), Ref(TF(Inf)),
            Ref(TF(Inf))))
        push!(indices, i)
    end
    isempty(indices) && error(
        "HybridWakePotential requires at least one Dirichlet body")
    return HybridWakePotentialState(states, indices, Ref(-1))
end

function initialize_formulation(f::GreenReconstruction, systems_tuple,
        wakes_tuple, body_solvers, backend_solve, backend_system)
    body, _ = _validate_formulation_common(f, systems_tuple, wakes_tuple,
        body_solvers)
    TF = eltype(body.strength)
    N = body.ncells
    green = _build_green_state(body, f.gauge, f.green_solver)
    return GreenReconstructionState{TF, typeof(green)}(green,
        zeros(TF, 3, N), zeros(TF, N), zeros(TF, N), zeros(TF, N),
        zeros(TF, N), Ref(-1))
end

function initialize_formulation(f::TraceCorrected, systems_tuple, wakes_tuple,
        body_solvers, backend_solve, backend_system)
    @warn "TraceCorrected is deprecated (2026-07-29): use GreenReconstruction "*
        "instead. Its :green estimator is identical to GreenReconstruction, "*
        "the affine-c channel was rejected as a lever by Gate A0, and its "*
        "post-solve velocity channel is vacuous (influence! never applies "*
        "wake_strength_shift under any backend, so surface velocities omit "*
        "the -c attached-strip term). See the TraceCorrected docstring." maxlog=1
    body, _ = _validate_formulation_common(f, systems_tuple, wakes_tuple,
        body_solvers)
    backend_system isa DirectBackend ||
        error("TraceCorrected requires backend_system isa DirectBackend: the "*
              "attached-wake strength correction does not yet flow through "*
              "FMM source buffers. Got $(typeof(backend_system)).")
    TF = eltype(body.strength)
    N = body.ncells
    edges = _shedding_edge_map(body)
    M = length(edges.upper)
    W = zeros(TF, N, M)
    _assemble_W!(W, body, edges)
    green = f.estimator === :green ?
        _build_green_state(body, f.gauge, f.green_solver) : nothing
    n_ppe = 0
    if f.estimator === :line_integral
        n_segment_probes = f.quadrature === :simpson ? 1 : 0
        if f.interior_path === :deformed
            # interior waypoint + optional per-segment midpoints (2 segments)
            n_ppe = 1 + 2*n_segment_probes
        else
            n_ppe = n_segment_probes
        end
    end
    probes = n_ppe > 0 ? FastMultipole.ProbeSystem(n_ppe*M, TF) : nothing
    return TraceCorrectedState{TF, typeof(green), typeof(probes)}(edges, W,
        zeros(TF, 3, N), zeros(TF, N), green, probes, n_ppe,
        zeros(TF, M), zeros(TF, M), zeros(TF, N), Ref(-1))
end

function initialize_formulation(f::DirectWakePotential, systems_tuple,
        wakes_tuple, body_solvers, backend_solve, backend_system)
    body, _ = _validate_formulation_common(f, systems_tuple, wakes_tuple,
        body_solvers)
    # every active wake must be a finite PanelWake whose sources all expose a
    # scalar potential (no trailing vector-potential-only filament, no
    # particles)
    for w in wakes_tuple
        isnothing(w) && continue
        w isa PanelWake ||
            error("DirectWakePotential requires a finite PanelWake with "*
                  "include_final_filament=false; got $(typeof(w)). Particle "*
                  "and mixed wakes cannot supply a complete scalar potential.")
        w.include_final_filament &&
            error("DirectWakePotential requires PanelWake "*
                  "include_final_filament=false: the trailing semi-infinite "*
                  "filament is vector-potential-only.")
        for src in get_sources(w)
            FastMultipole.has_vector_potential(src) &&
                error("DirectWakePotential wake source $(typeof(src)) exposes "*
                      "only a vector potential; a scalar potential is required.")
        end
    end
    TF = eltype(body.strength)
    N = body.ncells
    return DirectWakePotentialState{TF}(zeros(TF, 3, N), zeros(TF, N),
        zeros(TF, N), zeros(TF, N), zeros(TF, N), zeros(TF, N), Ref(-1))
end

################################################################################
# OPERATOR ASSEMBLY
################################################################################

_panel_areas(body::AbstractBody) = calc_areas(body.nodes, body.cells)

"Assemble the body-only interior-trace doublet operator B (including the +μ/2
self limit) by suppressing the attached-wake contribution of shedding panels."
function _assemble_B!(B, body::RigidWakeBody)
    was_active = _operator_mode_begin!(body)
    body.suppress_attached_wake[] = true
    try
        _G!(B, body, body; core_size=body.core_size_panel,
            update_geometry=false)
    finally
        body.suppress_attached_wake[] = false
        _operator_mode_end!(body, was_active)
    end
    return B
end

"Assemble the N×M attached-wake potential columns W: `W[i, e]` is the potential
at body centroid i of the attached transition panel of shedding edge e carrying
unit strength. Uses the exact in-matrix attached-panel kernel (difference of
`induced` with the attached wake on and off) under `core_size_panel`."
function _assemble_W!(W, body::RigidWakeBody, edges::SheddingEdgeMap)
    TF = eltype(W)
    old_strength = copy(body.strength)
    old_offset = body.core_size
    was_active = _operator_mode_begin!(body)
    body.core_size = body.core_size_panel
    body.strength .= zero(TF)
    switch = FastMultipole.DerivativesSwitch(true, false, false)
    CPs = body.controlpoints
    try
        for (e, up) in enumerate(edges.upper)
            body.strength[up, 2] = one(TF)
            for i in 1:body.ncells
                target = SVector{3,TF}(CPs[1,i], CPs[2,i], CPs[3,i])
                phi_full, _, _ = induced(target, body, up, switch;
                    core_size=body.core_size_panel)
                body.suppress_attached_wake[] = true
                phi_body, _, _ = induced(target, body, up, switch;
                    core_size=body.core_size_panel)
                body.suppress_attached_wake[] = false
                W[i, e] = phi_full - phi_body
            end
            body.strength[up, 2] = zero(TF)
        end
    finally
        body.strength .= old_strength
        body.core_size = old_offset
        body.suppress_attached_wake[] = false
        _operator_mode_end!(body, was_active)
    end
    return W
end

function _build_green_solve_state(body::RigidWakeBody, gauge::Symbol)
    TF = eltype(body.strength)
    N = body.ncells
    B = Matrix{TF}(undef, N, N)
    _assemble_B!(B, body)
    a = TF.(_panel_areas(body))

    if gauge === :area_mean
        K = Matrix{TF}(undef, N+1, N+1)
        @views begin
            K[1:N, 1:N] .= .-B
            for i in 1:N
                K[i, i] += one(TF)
            end
            K[1:N, N+1] .= a
            K[N+1, 1:N] .= a
        end
        K[N+1, N+1] = zero(TF)
        fact = lu!(K)
    else # :lsq
        A = Matrix{TF}(undef, N+1, N)
        @views begin
            A[1:N, :] .= .-B
            for i in 1:N
                A[i, i] += one(TF)
            end
            A[N+1, :] .= a
        end
        fact = LA.qr(A)
    end
    return GreenSolveState{TF, typeof(fact)}(fact, gauge, zeros(TF, N),
        zeros(TF, N+1), zeros(TF, N+1))
end

"Solve the gauge-fixed body-only Green system for the trace `q` given the
source potential `Sσ` at body centroids."
function _green_solve_q!(gs::GreenSolveState, Ssigma::AbstractVector)
    N = length(gs.q)
    gs.rhs_b[1:N] .= Ssigma
    gs.rhs_b[N+1] = zero(eltype(gs.rhs_b))
    if gs.gauge === :area_mean
        ldiv!(gs.sol_b, gs.fact, gs.rhs_b)
        gs.q .= view(gs.sol_b, 1:N)
    else # :lsq (min-residual with the area row; overdetermined (N+1)×N)
        gs.q .= gs.fact \ gs.rhs_b
    end
    return gs.q
end
_green_solve_q!(gs::GreenSolveState, Ssigma::AbstractVector, body, backend) =
    _green_solve_q!(gs, Ssigma)

"Route dispatch for the Green-system solve state: `nothing` → dense bordered
factorization (`Backslash`-analogue); `KrylovSolver` → matrix-free bordered
Krylov; `FGSSolver` → matrix-free relaxed Picard."
_build_green_state(body::RigidWakeBody, gauge::Symbol, ::Nothing) =
    _build_green_solve_state(body, gauge)
_build_green_state(body::RigidWakeBody, gauge::Symbol, ks::KrylovSolver) =
    _build_green_krylov_state(body, ks)
_build_green_state(body::RigidWakeBody, gauge::Symbol, fgs::FGSSolver) =
    _build_green_picard_state(body, fgs)

"Run `f()` with the body prepared for body-only linear operator products:
strengths/potential saved and restored, affine wake correction off (operator
mode), attached-wake influence suppressed."
function _with_green_scratch(f, body::RigidWakeBody)
    old_strength = copy(body.strength)
    old_potential = copy(body.potential)
    was_active = _operator_mode_begin!(body)
    old_suppress = body.suppress_attached_wake[]
    body.suppress_attached_wake[] = true
    try
        return f()
    finally
        body.strength .= old_strength
        body.potential .= old_potential
        body.suppress_attached_wake[] = old_suppress
        _operator_mode_end!(body, was_active)
    end
end

"Matrix-free product `out = B·x` (body-only interior-trace doublet operator)
through `influence!` with the selected backend. Must run inside
`_with_green_scratch` (which owns the flag toggles and state restoration)."
function _green_B_product!(out::AbstractVector, body::RigidWakeBody,
        x::AbstractVector, backend)
    body.strength[:, 1] .= zero(eltype(body.strength))
    body.strength[:, 2] .= x
    body.potential .= zero(eltype(body.potential))
    influence!(body, body, backend; scalar_potential=true, velocity=false)
    out .= body.potential
    return out
end

function _build_green_krylov_state(body::RigidWakeBody, ks::KrylovSolver)
    TF = eltype(body.strength)
    N = body.ncells
    ks.body.ncells == N ||
        error("green_solver KrylovSolver was built for a body with "*
              "$(ks.body.ncells) cells; the formulation body has $N.")
    isnothing(ks.preconditioner) ||
        @warn "green_solver KrylovSolver preconditioner is ignored for the "*
              "Green system: it was built for the wake-coupled operator G, "*
              "not (I−B)." maxlog=1

    a = TF.(_panel_areas(body))
    a_hat = a ./ LA.norm(a)
    q = zeros(TF, N)
    rhs_b = zeros(TF, N+1)
    sol_b = zeros(TF, N+1)
    Bx = zeros(TF, N)
    backend = ks.backend

    # bordered operator: [ (I−B)  â ; âᵀ  0 ] — same solution as the dense
    # :area_mean route (up to the â vs a column scaling of λ, which is not
    # part of q)
    prod! = (y, x, α, β) -> begin
        xq = view(x, 1:N)
        _green_B_product!(Bx, body, xq, backend)
        if iszero(β)
            @views y[1:N] .= α .* (xq .- Bx .+ x[N+1] .* a_hat)
            y[N+1] = α * LA.dot(a_hat, xq)
        else
            @views y[1:N] .= β .* y[1:N] .+ α .* (xq .- Bx .+ x[N+1] .* a_hat)
            y[N+1] = β * y[N+1] + α * LA.dot(a_hat, xq)
        end
        y
    end
    op = LinearOperators.LinearOperator(TF, N+1, N+1, false, false, prod!)
    workspace = Krylov.krylov_workspace(Val(ks.method), op, rhs_b)
    return GreenKrylovState{TF, typeof(body), typeof(ks), typeof(op), typeof(workspace)}(
        body, ks, a_hat, q, rhs_b, sol_b, Bx, op, workspace, Ref(false))
end

function _green_solve_q!(gs::GreenKrylovState, Ssigma::AbstractVector,
        body::RigidWakeBody, backend)
    N = length(gs.q)
    gs.rhs_b[1:N] .= Ssigma
    gs.rhs_b[N+1] = zero(eltype(gs.rhs_b))
    ks = gs.ks
    # the operator closure mutates gs.body (captured at initialization); the
    # scratch/flag window must target that same instance
    _with_green_scratch(gs.body) do
        if gs.have_warm[]
            Krylov.krylov_solve!(gs.workspace, gs.op, gs.rhs_b, gs.sol_b;
                atol=ks.atol, rtol=ks.rtol, itmax=ks.itmax)
        else
            Krylov.krylov_solve!(gs.workspace, gs.op, gs.rhs_b;
                atol=ks.atol, rtol=ks.rtol, itmax=ks.itmax)
        end
    end
    stats = gs.workspace.stats
    stats.solved ||
        @warn "Green Krylov solve did not converge: $(stats.status) after "*
              "$(stats.niter) iterations." maxlog=5
    gs.sol_b .= gs.workspace.x
    gs.q .= view(gs.sol_b, 1:N)
    gs.have_warm[] = true
    return gs.q
end

function _build_green_picard_state(body::RigidWakeBody, fgs::FGSSolver)
    TF = eltype(body.strength)
    N = body.ncells
    a = TF.(_panel_areas(body))
    return GreenPicardState{TF, typeof(fgs)}(fgs, a, sum(a),
        zeros(TF, N), zeros(TF, N), zeros(TF, N))
end

function _green_solve_q!(gs::GreenPicardState, Ssigma::AbstractVector,
        body::RigidWakeBody, backend)
    fgs = gs.fgs
    rlx = fgs.rlx
    converged = false
    defect = convert(eltype(gs.q), Inf)
    niter = 0
    _with_green_scratch(body) do
        for _ in 1:fgs.max_iterations
            niter += 1
            _green_B_product!(gs.Bx, body, gs.q, backend)
            gs.q_new .= Ssigma .+ gs.Bx
            # area-mean projection kills the constant mode (B's unit
            # eigenvector), on which the Picard map does not contract
            gs.q_new .-= LA.dot(gs.a, gs.q_new)/gs.suma
            defect = LA.norm(gs.q_new .- gs.q) /
                max(LA.norm(gs.q_new), eps(eltype(gs.q)))
            @. gs.q = (1 - rlx)*gs.q + rlx*gs.q_new
            if defect <= fgs.tolerance
                converged = true
                break
            end
        end
    end
    converged ||
        @warn "Green Picard iteration did not converge: relative defect "*
              "$defect after $niter sweeps (tolerance $(fgs.tolerance))." maxlog=5
    return gs.q
end

"Evaluate the body source-panel potential `Sσ` at body centroids for an
arbitrary source coefficient vector, preserving body strengths and potential."
function _source_potential!(out::AbstractVector, body::AbstractBody,
        sigma::AbstractVector, backend)
    old_strength = copy(body.strength)
    old_potential = copy(body.potential)
    was_active = _operator_mode_begin!(body)
    try
        body.strength[:, 1] .= sigma
        body.strength[:, 2] .= zero(eltype(body.strength))
        body.potential .= zero(eltype(body.potential))
        influence!(body, body, backend; scalar_potential=true, velocity=false)
        out .= body.potential
    finally
        body.strength .= old_strength
        body.potential .= old_potential
        _operator_mode_end!(body, was_active)
    end
    return out
end

"Evaluate the free wake's scalar potential `q_f` directly at body centroids,
preserving body potential. Only wake sources that expose a scalar potential are
collected (particle fields are excluded upstream)."
function _wake_potential!(out::AbstractVector, body::AbstractBody,
        wakes_tuple::Tuple, backend)
    sources = _collect_wake_scalar_sources(wakes_tuple)
    old_potential = copy(body.potential)
    try
        body.potential .= zero(eltype(body.potential))
        if length(sources) > 0
            influence!((body,), sources, backend; scalar_potential=true,
                velocity=false, precalc=true)
        end
        out .= body.potential
    finally
        body.potential .= old_potential
    end
    return out
end

"Split the accumulated control-point velocity into wake-only σ and non-wake σ0
source coefficients using the pre-wake snapshot."
function _split_sigma!(sigma::AbstractVector, sigma0::AbstractVector,
        body::AbstractBody, u_prewake::AbstractMatrix)
    normals = body.normals
    velocity = body.velocity
    for i in eachindex(sigma)
        s0 = zero(eltype(sigma))
        st = zero(eltype(sigma))
        for d in 1:3
            s0 -= u_prewake[d, i]*normals[d, i]
            st -= velocity[d, i]*normals[d, i]
        end
        sigma0[i] = s0
        sigma[i] = st - s0
    end
    return sigma, sigma0
end

################################################################################
# PER-STEP HOOKS
################################################################################

"""
    formulation_prewake!(formulation, state, systems_tuple)

Called by `_steady_aerodynamics!` immediately before the wake-influence pass.
Snapshotting formulations record the pre-wake control-point velocity so the
wake-only contribution can be isolated afterwards.
"""
formulation_prewake!(::VelocityThroughSources, state, systems_tuple) = nothing
function formulation_prewake!(::HybridWakePotential,
        state::HybridWakePotentialState, systems_tuple)
    for (bs, i) in zip(state.bodies, state.body_indices)
        bs.u_prewake .= systems_tuple[i].velocity
    end
    return nothing
end
formulation_prewake!(::AbstractSolveFormulation, state, systems_tuple) =
    (state.u_prewake .= systems_tuple[1].velocity; nothing)

"""
    solve_formulation!(formulation, state, systems, systems_tuple, wakes_tuple,
                       body_solvers; backend_solve, backend_wake, i_step)

Formulation-dispatched replacement for the body solve inside
`_steady_aerodynamics!`. The default method is exactly the production solve.
"""
function solve_formulation!(f::VelocityThroughSources, state, systems,
        systems_tuple, wakes_tuple, body_solvers;
        backend_solve, backend_wake, i_step::Int=0)
    t0 = time()   # task 052: env-gated solve timer (FLOWPANEL_GPU_TIMERS)
    if systems isa Tuple
        # GS_VERBOSE=true prints the per-iteration normalized block residual
        # (solve!'s own verbose path). Env-gated, default off: each verbose
        # residual measurement costs ~one influence sweep per outer iteration.
        gs_verbose = lowercase(get(ENV, "GS_VERBOSE", "false")) in ("true", "1")
        solve!(systems, body_solvers; backend=backend_solve,
            max_outer_iterations=f.max_outer_iterations,
            outer_tolerance=f.outer_tolerance,
            dirichlet_residual_scale=f.dirichlet_residual_scale,
            neumann_residual_scale=f.neumann_residual_scale,
            require_outer_convergence=f.require_outer_convergence,
            verbose=gs_verbose)
    else
        # Preserve the exact historical single-body solve call. Outer-loop
        # controls apply only to tuple/block-GS orchestration.
        solve!(systems, body_solvers; backend=backend_solve)
    end
    _gpu_timer_log("solve step=$(i_step)", time() - t0)
    return nothing
end

function _wake_panel_velocity!(out::AbstractMatrix, body::AbstractBody,
        wakes_tuple::Tuple, backend)
    sources = ()
    for w in wakes_tuple
        isnothing(w) && continue
        pw = w isa PanelParticleWake ? w.panel_wake : w
        pw isa PanelWake || continue
        sources = (sources..., get_sources(pw)...)
    end
    old_velocity = copy(body.velocity)
    try
        body.velocity .= zero(eltype(body.velocity))
        isempty(sources) || influence!((body,), sources, backend;
            scalar_potential=false, velocity=true, precalc=true)
        out .= body.velocity
    finally
        body.velocity .= old_velocity
    end
    return out
end

function _green_diagnostics!(state::HybridBodyState, body::RigidWakeBody)
    Bq = similar(state.green.q)
    _with_green_scratch(body) do
        _green_B_product!(Bq, body, state.green.q, DirectBackend())
    end
    areas = eltype(Bq).(_panel_areas(body))
    # Dense :area_mean uses the bordered equation
    # (I-B)q + a*lambda = Ssigma. Include the compatibility multiplier in
    # the reported linear-system residual instead of mislabelling it as error.
    lambda = state.green.gauge === :area_mean ? state.green.sol_b[end] :
        zero(eltype(Bq))
    defect = state.green.q .- Bq .+ areas .* lambda .- state.Ssigma
    state.green_residual[] = LA.norm(defect) /
        max(LA.norm(state.Ssigma), eps(eltype(defect)))
    state.gauge_defect[] = abs(LA.dot(areas, state.green.q)) /
        max(LA.norm(areas) * LA.norm(state.green.q), eps(eltype(defect)))
    return nothing
end

_has_active_particles(wakes_tuple::Tuple) = any(w ->
    w isa PanelParticleWake && w.pfield.np > 0, wakes_tuple)

"""
    surface_hodge_trace!(q, body, velocity)

Independent edge-based least-squares surface-Hodge reconstruction. Adjacent
panel traces are fitted so `q[j]-q[i]` matches the line integral of the
sampled tangential velocity between their control points, with an
area-weighted zero-mean gauge. Returns the relative tangential projection
defect. This diagnostic does not participate in the solve.
"""
function surface_hodge_trace!(q::AbstractVector, body::AbstractBody,
        velocity::AbstractMatrix)
    owners = Dict{Tuple{Int,Int},Int}()
    pairs = Tuple{Int,Int}[]
    for panel in axes(body.cells, 2)
        c = view(body.cells, :, panel)
        for (a, b) in ((c[1], c[2]), (c[2], c[3]), (c[3], c[1]))
            edge = a < b ? (a, b) : (b, a)
            if haskey(owners, edge)
                push!(pairs, (owners[edge], panel))
            else
                owners[edge] = panel
            end
        end
    end
    TF = eltype(q)
    m = length(pairs)
    A = zeros(TF, m + 1, body.ncells)
    rhs = zeros(TF, m + 1)
    for (row, (i, j)) in enumerate(pairs)
        A[row, i] = -one(TF)
        A[row, j] = one(TF)
        dx = view(body.controlpoints, :, j) .- view(body.controlpoints, :, i)
        u = (view(velocity, :, i) .+ view(velocity, :, j)) ./ 2
        n = view(body.normals, :, i) .+ view(body.normals, :, j)
        nn = LA.norm(n)
        ut = nn > eps(TF) ? u .- n .* (LA.dot(u, n) / (nn * nn)) : u
        rhs[row] = LA.dot(ut, dx)
    end
    areas = TF.(_panel_areas(body))
    A[end, :] .= areas ./ max(LA.norm(areas), eps(TF))
    q .= A \ rhs
    defect = LA.norm(view(A, 1:m, :) * q - view(rhs, 1:m)) /
        max(LA.norm(view(rhs, 1:m)), eps(TF))
    return defect
end

function solve_formulation!(f::HybridWakePotential,
        state::HybridWakePotentialState, systems, systems_tuple, wakes_tuple,
        body_solvers; backend_solve, backend_wake, i_step::Int=0)
    saved_velocity = [copy(body.velocity) for body in systems_tuple]
    saved_potential = [copy(body.potential) for body in systems_tuple]
    recompute = state.last_recompute[] < 0 ||
        i_step - state.last_recompute[] >= f.recompute_interval
    has_particles = _has_active_particles(wakes_tuple)
    try
        for (bs, i) in zip(state.bodies, state.body_indices)
            body = systems_tuple[i]
            solve_backend = backend_solve isa Tuple ||
                backend_solve isa AbstractVector ? backend_solve[i] : backend_solve
            clear_wake_correction!(body)
            if recompute
                if has_particles
                    _wake_panel_velocity!(bs.panel_velocity, body, wakes_tuple,
                        backend_wake)
                    # total sampled wake velocity minus retained panel-wake
                    # velocity leaves the particle-only contribution
                    for k in eachindex(bs.sigma_particle)
                        sigma = zero(eltype(bs.sigma_particle))
                        for d in 1:3
                            up = saved_velocity[i][d, k] - bs.u_prewake[d, k] -
                                bs.panel_velocity[d, k]
                            bs.particle_velocity[d, k] = up
                            sigma -= up * body.normals[d, k]
                        end
                        bs.sigma_particle[k] = sigma
                    end
                    _source_potential!(bs.Ssigma, body, bs.sigma_particle,
                        solve_backend)
                    _green_solve_q!(bs.green, bs.Ssigma, body, solve_backend)
                    _green_diagnostics!(bs, body)
                    bs.tangential_projection_defect[] = surface_hodge_trace!(
                        bs.q_hodge, body, bs.particle_velocity)
                    bs.green_hodge_mismatch[] = LA.norm(
                        bs.green.q .- bs.q_hodge) / max(LA.norm(bs.green.q),
                        LA.norm(bs.q_hodge), eps(eltype(bs.q_hodge)))
                else
                    fill!(bs.panel_velocity, 0)
                    fill!(bs.particle_velocity, 0)
                    fill!(bs.sigma_particle, 0)
                    fill!(bs.Ssigma, 0)
                    fill!(bs.green.q, 0)
                    bs.green_residual[] = 0
                    bs.gauge_defect[] = 0
                    fill!(bs.q_hodge, 0)
                    bs.green_hodge_mismatch[] = 0
                    bs.tangential_projection_defect[] = 0
                end
                _wake_potential!(bs.q_panel, body, wakes_tuple, backend_wake)
                bs.q_total .= bs.q_panel .+ bs.green.q
            end
            # Frozen Dirichlet sources see only freestream/kinematics. The
            # prescribed wake trace is carried in potential and preserved by
            # the corrected tuple solve.
            body.velocity .= bs.u_prewake
            # `body.potential` on entry is accumulated simulation workspace
            # from prior passes/steps, not an independently prescribed
            # incident trace. The hybrid RHS owns its wake trace explicitly.
            body.potential .= bs.q_total
        end
        recompute && (state.last_recompute[] = i_step)
        solvers_tuple = body_solvers isa Tuple ? body_solvers : (body_solvers,)
        solve!(systems_tuple, solvers_tuple; backend=backend_solve,
            max_outer_iterations=f.max_outer_iterations,
            outer_tolerance=f.outer_tolerance,
            dirichlet_residual_scale=f.dirichlet_residual_scale,
            neumann_residual_scale=f.neumann_residual_scale,
            require_outer_convergence=f.require_outer_convergence)
    finally
        for (i, body) in enumerate(systems_tuple)
            body.velocity .= saved_velocity[i]
            body.potential .= saved_potential[i]
        end
    end
    return nothing
end

function solve_formulation!(f::GreenReconstruction,
        state::GreenReconstructionState, systems, systems_tuple, wakes_tuple,
        body_solvers; backend_solve, backend_wake, i_step::Int=0)
    body = systems_tuple[1]
    solver = _single_body_solver(body_solvers)

    # this formulation stores the explicit-potential μE; the Kutta map applies
    # to it directly with no affine correction
    clear_wake_correction!(body)

    _split_sigma!(state.sigma, state.sigma0, body, state.u_prewake)

    # reconstruct the wake trace q from the wake-only source coefficients
    if state.last_recompute[] < 0 ||
            i_step - state.last_recompute[] >= f.recompute_interval
        _source_potential!(state.Ssigma, body, state.sigma, backend_solve)
        _green_solve_q!(state.green, state.Ssigma, body, backend_solve)
        state.last_recompute[] = i_step
    end

    # explicit-potential solve G·μE = −S·σ0 − q, reusing the existing LU
    # (the ordinary single-body solve! is bypassed because it discards
    # preassembled external potential)
    _source_potential!(state.Ssigma, body, state.sigma0, backend_solve)
    state.rhs .= .-state.Ssigma .- state.green.q
    solver.rhs .= state.rhs
    ldiv!(view(body.strength, :, 2), solver.Glu, solver.rhs)
    # sources carry σ0 only: the free wake's field is represented by the
    # actual wake system in this formulation
    body.strength[:, 1] .= state.sigma0
    return nothing
end

function solve_formulation!(f::DirectWakePotential,
        state::DirectWakePotentialState, systems, systems_tuple, wakes_tuple,
        body_solvers; backend_solve, backend_wake, i_step::Int=0)
    body = systems_tuple[1]
    solver = _single_body_solver(body_solvers)

    # this formulation stores the explicit-potential μE directly; no affine
    # Kutta correction is applied
    clear_wake_correction!(body)

    # σ0 = non-wake normal velocity. The wake-induced σ is computed for
    # completeness but NOT stored on the body: the free wake's field is carried
    # by the actual wake system, entering the RHS through the direct q_f below.
    _split_sigma!(state.sigma, state.sigma0, body, state.u_prewake)

    # direct scalar potential q_f from the finite panel wake at body centroids
    # (recompute_interval lags this on the wake-convection timescale)
    if state.last_recompute[] < 0 ||
            i_step - state.last_recompute[] >= f.recompute_interval
        _wake_potential!(state.q_wake, body, wakes_tuple, backend_wake)
        state.last_recompute[] = i_step
    end

    # explicit-potential solve G_Δ·μE = −S·σ0 − q_f, reusing the Backslash LU
    # (single-body solve! is bypassed because it discards preassembled external
    # potential; see task3.md)
    _source_potential!(state.Ssigma, body, state.sigma0, backend_solve)
    state.rhs .= .-state.Ssigma .- state.q_wake
    solver.rhs .= state.rhs
    ldiv!(view(body.strength, :, 2), solver.Glu, solver.rhs)
    # sources carry σ0 only
    body.strength[:, 1] .= state.sigma0
    return nothing
end

function solve_formulation!(f::TraceCorrected, state::TraceCorrectedState,
        systems, systems_tuple, wakes_tuple, body_solvers;
        backend_solve, backend_wake, i_step::Int=0)
    body = systems_tuple[1]
    solver = _single_body_solver(body_solvers)

    # operator mode for the entire solve: the previous step's correction must
    # not contaminate influence evaluations
    body.wake_correction_active[] = false

    _update_c!(f, state, body, wakes_tuple; backend_solve, backend_wake, i_step)

    # production-equivalent solve with the single added RHS term +W·c
    # (replicates the single-body Dirichlet solve!, which cannot be reused
    # because it clears preassembled external potential)
    # This formulation drives the raw `_solve!` kernel directly (the public
    # single-body solve! cannot be reused — it clears preassembled external
    # potential), so it owns the top-level step boundary itself (021 Phase 3).
    begin_step_solution!(solver)

    potential_old = copy(body.potential)
    try
        set_strengths!(body)
        body.potential .= zero(eltype(body.potential))
        influence!(body, body, backend_solve; scalar_potential=true,
            velocity=false)
        LA.mul!(state.work, state.W, state.c)
        body.potential .-= state.work     # rhs = −potential = −S(σ0+σ) + W·c
        _solve!(body, solver; backend=backend_solve)
        note_step_solve!(solver)
    finally
        body.potential .= potential_old
    end

    # physical mode: γ = C·μ̃ − c for every downstream consumer
    set_wake_correction!(body, state.c)

    # Only now is the formulation complete: finalizing inside the `try` above
    # would commit the step before the physical wake correction is installed.
    finalize_step_solution!(body, solver)
    return nothing
end

################################################################################
# KUTTA-TRACE ESTIMATORS
################################################################################

function _update_c!(f::TraceCorrected, state::TraceCorrectedState,
        body::RigidWakeBody, wakes_tuple; backend_solve, backend_wake,
        i_step::Int)
    first_compute = state.last_recompute[] < 0
    if !first_compute && i_step - state.last_recompute[] < f.recompute_interval
        return state.c
    end

    if f.estimator === :green
        _split_sigma!(state.sigma, state.work, body, state.u_prewake)
        _source_potential!(state.work, body, state.sigma, backend_solve)
        _green_solve_q!(state.green, state.work, body, backend_solve)
        _apply_kutta_map!(state.c_new, state.green.q, state.edges)
    else
        _line_integral_c!(f, state, body, wakes_tuple; backend_wake)
    end

    # unpaired shedding edges carry no potential difference; force the gauge
    # convention c = 0 there
    for e in eachindex(state.c_new)
        state.edges.lower[e] > 0 || (state.c_new[e] = zero(eltype(state.c_new)))
    end

    if first_compute
        state.c .= state.c_new
    else
        α = f.relaxation
        @. state.c = (1 - α)*state.c + α*state.c_new
    end
    state.last_recompute[] = i_step
    return state.c
end

_wake_velocity_at(state::TraceCorrectedState, body, i::Int) = SVector{3}(
    body.velocity[1, i] - state.u_prewake[1, i],
    body.velocity[2, i] - state.u_prewake[2, i],
    body.velocity[3, i] - state.u_prewake[3, i])

function _zero_probes!(probes)
    TF = eltype(probes.scalar_potential)
    zero_v = zero(eltype(probes.gradient))
    zero_h = zero(eltype(probes.hessian))
    for k in eachindex(probes.position)
        probes.gradient[k] = zero_v
        probes.hessian[k] = zero_h
        probes.scalar_potential[k] = zero(TF)
    end
    return probes
end

"Mean upstream wake direction of shedding edge e (unit vector), from the two
`Das` columns bracketing the edge."
function _edge_wake_direction(body::RigidWakeBody, edges::SheddingEdgeMap,
        e::Int)
    Das = body.Das[edges.isurf[e]]
    i = edges.iedge[e]
    d = SVector{3}(Das[1, i] + Das[1, i+1],
                   Das[2, i] + Das[2, i+1],
                   Das[3, i] + Das[3, i+1])
    nrm = LA.norm(d)
    return nrm > 0 ? d/nrm : d
end

function _line_integral_c!(f::TraceCorrected, state::TraceCorrectedState,
        body::RigidWakeBody, wakes_tuple; backend_wake)
    edges = state.edges
    CPs = body.controlpoints
    TF = eltype(state.c_new)

    # batched wake-only velocity at any extra quadrature nodes
    if state.n_probes_per_edge > 0
        _set_line_integral_probes!(f, state, body)
        _zero_probes!(state.probes)
        wake_sources = _collect_wake_sources(wakes_tuple)
        if length(wake_sources) > 0
            influence!((state.probes,), wake_sources, backend_wake;
                precalc=false, scalar_potential=false, gradient=true,
                hessian=false)
        end
    end

    for e in eachindex(edges.upper)
        lo = edges.lower[e]
        if lo <= 0
            state.c_new[e] = zero(TF)
            continue
        end
        up = edges.upper[e]
        xu = SVector{3,TF}(CPs[1, up], CPs[2, up], CPs[3, up])
        xl = SVector{3,TF}(CPs[1, lo], CPs[2, lo], CPs[3, lo])
        uu = _wake_velocity_at(state, body, up)
        ul = _wake_velocity_at(state, body, lo)

        if f.interior_path === :straight
            d = xu - xl
            if f.quadrature === :trapezoid
                state.c_new[e] = LA.dot(0.5*(uu + ul), d)
            else # :simpson, midpoint probe at slot e
                um = state.probes.gradient[e]
                state.c_new[e] = LA.dot(ul + 4*um + uu, d)/6
            end
        else # :deformed two-segment path xl → xi → xu
            k0 = state.n_probes_per_edge*(e - 1)
            xi = state.probes.position[k0 + 1]
            ui = state.probes.gradient[k0 + 1]
            if f.quadrature === :trapezoid
                state.c_new[e] = LA.dot(0.5*(ul + ui), xi - xl) +
                                 LA.dot(0.5*(ui + uu), xu - xi)
            else # :simpson, per-segment midpoints in slots k0+2, k0+3
                um1 = state.probes.gradient[k0 + 2]
                um2 = state.probes.gradient[k0 + 3]
                state.c_new[e] = LA.dot(ul + 4*um1 + ui, xi - xl)/6 +
                                 LA.dot(ui + 4*um2 + uu, xu - xi)/6
            end
        end
    end
    return state.c_new
end

function _set_line_integral_probes!(f::TraceCorrected,
        state::TraceCorrectedState, body::RigidWakeBody)
    edges = state.edges
    CPs = body.controlpoints
    TF = eltype(state.c_new)
    for e in eachindex(edges.upper)
        lo = edges.lower[e]
        up = edges.upper[e]
        # unpaired edges get placeholder probes at the upper control point
        # (results unused; c is forced to zero)
        xu = SVector{3,TF}(CPs[1, up], CPs[2, up], CPs[3, up])
        xl = lo > 0 ? SVector{3,TF}(CPs[1, lo], CPs[2, lo], CPs[3, lo]) : xu
        if f.interior_path === :straight
            # single midpoint probe per edge (Simpson)
            state.probes.position[e] = 0.5*(xl + xu)
        else
            xm = 0.5*(xl + xu)
            h = LA.norm(xu - xl)
            xi = xm - TF(f.path_depth)*h*_edge_wake_direction(body, edges, e)
            k0 = state.n_probes_per_edge*(e - 1)
            state.probes.position[k0 + 1] = xi
            if f.quadrature === :simpson
                state.probes.position[k0 + 2] = 0.5*(xl + xi)
                state.probes.position[k0 + 3] = 0.5*(xi + xu)
            end
        end
    end
    return state.probes
end
