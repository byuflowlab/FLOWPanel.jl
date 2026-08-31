#=##############################################################################
# DESCRIPTION
    Independently configurable wake attachment and Kutta closure for a single
    finite-wake, fully paired, Dirichlet RigidWakeBody (BRAINSTORM item 015,
    Phase 3). Implements the architecture approved in
    BRAINSTORM/015_pressure_continuity_kutta_condition/phase_02_architecture.md:

    * wake attachment axis: `RigidTransitionAttachment` (Route A, the existing
      body-owned TE → TE+Das transition) or `TEAnchoredAttachment` (Route B, a
      live wake panel from the current TE to the first convected wake row, with
      no independent rigid `Das`);
    * Kutta closure axis: `JumpKutta` (γ = Cμ inside the coupled solve) or
      `PressureContinuityKutta` (nonlinear paired-centroid pressure matching
      through the shared correction coordinate c = Cμ − γ).

    The exact legacy pair (RigidTransitionAttachment + JumpKutta) never reaches
    this file's runtime: `simulate!`/`steady!` branch into the pre-existing
    call sequence before any Kutta state is allocated.

    First implementation (user amendment 2026-07-29): only the `Backslash`
    body solver is supported; `KrylovSolver`/`FGSSolver` support is deferred.

# SIGN/CHANNEL CONVENTIONS
    The correction channel reuses the production affine machinery proven by
    `TraceCorrected` (src/FLOWPanel_formulation.jl):

    * during a coupled solve the correction is NOT active
      (`wake_correction_active[] == false`), so the assembled operator is the
      strictly linear G = B + W·C; the affine part enters through the
      right-hand side as `body.potential .-= W*c`;
    * after a solve, `set_wake_correction!(body, c)` installs γ = Cμ − c for
      every downstream consumer (velocity reconstruction through the
      DirectBackend index path, `shed_wake!` deposits, force recovery).

    Velocity reconstruction is backend-complete: the main influence passes run
    with the correction INACTIVE through any backend (the proportional γ = Cμ
    part flows through both the direct and FMM paths), and the affine part —
    sourced only by the M trailing-edge transition strips, one wake row of
    panels — is completed by the exact direct add-on
    `_add_affine_attached_velocity!` (O(M×N_targets), evaluated with the same
    triangle construction as the index-path `_induced_wake`). No FastMultipole
    buffer or multipole-moment changes are needed, and no backend restriction
    applies. (`TraceCorrected`'s own DirectBackend restriction is unchanged;
    it could adopt the same helper in a separate item.)
=###############################################################################

################################################################################
# PUBLIC STRATEGY TYPES
################################################################################

"Supertype of wake-attachment strategies (`RigidTransitionAttachment`,
`TEAnchoredAttachment`)."
abstract type AbstractWakeAttachment end

"Route A: retain the existing body-owned rigid `TE → TE + Das` attached
transition panel. This is the default and is never `Das`-free."
struct RigidTransitionAttachment <: AbstractWakeAttachment end

"Route B: TE-anchored live wake panel. The wake owns a node row at the
geometric trailing edge; the live panel connects it to the first convected
downstream wake row. There is no independent rigid `Das`, and every
`set_Das_*` option is rejected."
struct TEAnchoredAttachment <: AbstractWakeAttachment end

"Supertype of Kutta closures (`JumpKutta`, `PressureContinuityKutta`)."
abstract type AbstractKuttaClosure end

"Deterministic jump closure γ = Cμ enforced inside the same coupled body
solve (the existing production closure; equivalently c = 0)."
struct JumpKutta <: AbstractKuttaClosure end

"Supertype of trial-safe pressure providers used by
[`PressureContinuityKutta`](@ref)."
abstract type AbstractKuttaPressureProvider end

"""
    BroydenKutta(; max_iterations=30, max_backtracks=12, max_restarts=2,
                 max_scaled_step=2.0, min_linesearch=2.0^-20)

Configuration of the primary safeguarded dense inverse-Broyden strategy for
[`PressureContinuityKutta`](@ref). Not exported; construct as
`FLOWPanel.BroydenKutta(...)`.
"""
struct BroydenKutta
    max_iterations::Int
    max_backtracks::Int
    max_restarts::Int
    max_scaled_step::Float64
    min_linesearch::Float64

    function BroydenKutta(; max_iterations::Int=30, max_backtracks::Int=12,
            max_restarts::Int=2, max_scaled_step::Real=2.0,
            min_linesearch::Real=2.0^-20)
        max_iterations >= 1 || throw(ArgumentError("BroydenKutta max_iterations must be >= 1"))
        max_backtracks >= 0 || throw(ArgumentError("BroydenKutta max_backtracks must be >= 0"))
        max_restarts >= 0 || throw(ArgumentError("BroydenKutta max_restarts must be >= 0"))
        (isfinite(max_scaled_step) && max_scaled_step > 0) ||
            throw(ArgumentError("BroydenKutta max_scaled_step must be finite and positive"))
        (isfinite(min_linesearch) && 0 < min_linesearch <= 1) ||
            throw(ArgumentError("BroydenKutta min_linesearch must be in (0, 1]"))
        return new(max_iterations, max_backtracks, max_restarts,
            Float64(max_scaled_step), Float64(min_linesearch))
    end
end

"""
    NewtonKrylovKutta(; max_newton=12, max_krylov=30,
                      fd_relative_step=sqrt(eps(Float64)),
                      max_scaled_step=2.0, max_backtracks=12,
                      min_linesearch=2.0^-20)

Configuration of the finite-difference Newton–Krylov fallback strategy for
[`PressureContinuityKutta`](@ref). Not exported; construct as
`FLOWPanel.NewtonKrylovKutta(...)`.
"""
struct NewtonKrylovKutta
    max_newton::Int
    max_krylov::Int
    fd_relative_step::Float64
    max_scaled_step::Float64
    max_backtracks::Int
    min_linesearch::Float64

    function NewtonKrylovKutta(; max_newton::Int=12, max_krylov::Int=30,
            fd_relative_step::Real=sqrt(eps(Float64)),
            max_scaled_step::Real=2.0, max_backtracks::Int=12,
            min_linesearch::Real=2.0^-20)
        max_newton >= 1 || throw(ArgumentError("NewtonKrylovKutta max_newton must be >= 1"))
        max_krylov >= 1 || throw(ArgumentError("NewtonKrylovKutta max_krylov must be >= 1"))
        (isfinite(fd_relative_step) && fd_relative_step > 0) ||
            throw(ArgumentError("NewtonKrylovKutta fd_relative_step must be finite and positive"))
        (isfinite(max_scaled_step) && max_scaled_step > 0) ||
            throw(ArgumentError("NewtonKrylovKutta max_scaled_step must be finite and positive"))
        max_backtracks >= 0 || throw(ArgumentError("NewtonKrylovKutta max_backtracks must be >= 0"))
        (isfinite(min_linesearch) && 0 < min_linesearch <= 1) ||
            throw(ArgumentError("NewtonKrylovKutta min_linesearch must be in (0, 1]"))
        return new(max_newton, max_krylov, Float64(fd_relative_step),
            Float64(max_scaled_step), max_backtracks, Float64(min_linesearch))
    end
end

_is_kutta_strategy(::BroydenKutta) = true
_is_kutta_strategy(::NewtonKrylovKutta) = true
_is_kutta_strategy(::Any) = false

################################################################################
# DIAGNOSTICS AND ERRORS
################################################################################

"""
    KuttaDiagnostics

Per-physical-step record of an accepted (or failed) Kutta-closure solve.
`status` is one of `:converged` (ordinary pressure convergence),
`:startup_jump` (deterministic Route B cold-start jump solve),
`:jump_fallback` (pressure failure followed by a committed fresh jump solve),
or `:failed` (thrown; carried inside [`KuttaConvergenceError`](@ref)).
"""
struct KuttaDiagnostics{TF<:AbstractFloat}
    status::Symbol
    method::Symbol            # :none | :broyden | :newton_krylov
    route::Symbol             # :A | :B
    outer_iterations::Int
    body_solves::Int
    backtracks::Int
    restarts::Int
    r_inf::TF                 # dimensional pressure residual ∞-norm
    r_inf_scaled::TF
    r_l2_edge::TF             # edge-length-weighted pressure L2 norm
    r_l2_edge_scaled::TF
    dc_inf::TF                # dimensional accepted correction step ∞-norm
    dc_inf_scaled::TF
    pressure_scale::TF
    correction_scale::TF
    U_s::TF
    L_s::TF
    inner_status::Symbol      # :converged | :nonfinite | :not_run
    inner_iterations::Int     # 1 for a direct solve
    inner_residual::TF        # NaN for a direct solve
    fallback_entered::Bool
    startup::Symbol           # :none | :startup_jump
    disposition::Symbol       # :converged | :startup_jump | :jump_fallback | :error
    i_step::Int
    c::Vector{TF}             # accepted correction (restart state)
end

"""
    KuttaConvergenceError

Thrown when a `PressureContinuityKutta` step cannot be completed under
`on_failure=:error` (or when the explicit `:jump` fallback itself fails). The
pre-step state has been fully restored before the throw. Carries the final
uncommitted [`KuttaDiagnostics`](@ref) record and a `cause` symbol
(`:nonconvergence`, `:inner_solve`, `:nonfinite`, `:commit`,
`:jump_fallback`).
"""
struct KuttaConvergenceError{TF<:AbstractFloat} <: Exception
    diagnostics::KuttaDiagnostics{TF}
    cause::Symbol
    msg::String
end

function Base.showerror(io::IO, err::KuttaConvergenceError)
    print(io, "KuttaConvergenceError (cause=$(err.cause)): ", err.msg)
    d = err.diagnostics
    print(io, "\n  step $(d.i_step), route $(d.route), method $(d.method): ",
        "$(d.outer_iterations) outer iterations, $(d.body_solves) body solves, ",
        "$(d.backtracks) backtracks, $(d.restarts) restarts")
    print(io, "\n  ‖r‖∞/P_s = $(d.r_inf_scaled) (dimensional $(d.r_inf)), ",
        "‖Δc‖∞/C_s = $(d.dc_inf_scaled)")
end

################################################################################
# PRESSURE-CONTINUITY CLOSURE
################################################################################

"""
    PressureContinuityKutta(provider;
        primary=BroydenKutta(), fallback=NewtonKrylovKutta(),
        pressure_scale=:auto, correction_scale=:auto,
        pressure_tolerance=1e-6, correction_tolerance=1e-6,
        on_failure=:error, store_diagnostics=true)

Nonlinear Kutta closure selecting the circulation for which the two panels of
each paired shedding edge have equal centroid pressure, evaluated by
`provider` (see [`SteadyBernoulliProvider`](@ref)). Solved with safeguarded
dense inverse Broyden, falling back to finite-difference Newton–Krylov, on
the shared correction coordinate c = Cμ − γ.

`pressure_scale`/`correction_scale` may be `:auto` (frozen per physical step
from the initial clean trial: P_s = ½ρU_s², C_s = U_s·L_s) or explicit
finite positive values. `on_failure` is `:error` (restore and throw
[`KuttaConvergenceError`](@ref)) or `:jump` (restore, commit a fresh c = 0
jump solve, warn). Accepted per-step diagnostics are retrievable with
[`kutta_diagnostics`](@ref) when `store_diagnostics=true`.
"""
struct PressureContinuityKutta{P<:AbstractKuttaPressureProvider, S, F} <: AbstractKuttaClosure
    provider::P
    primary::S
    fallback::F
    pressure_scale::Float64      # NaN encodes :auto
    correction_scale::Float64    # NaN encodes :auto
    pressure_tolerance::Float64
    correction_tolerance::Float64
    on_failure::Symbol
    store_diagnostics::Bool
    diagnostics::Vector{KuttaDiagnostics{Float64}}
end

function PressureContinuityKutta(provider::AbstractKuttaPressureProvider;
        primary=BroydenKutta(),
        fallback=NewtonKrylovKutta(),
        pressure_scale::Union{Symbol, Real}=:auto,
        correction_scale::Union{Symbol, Real}=:auto,
        pressure_tolerance::Real=1e-6,
        correction_tolerance::Real=1e-6,
        on_failure::Symbol=:error,
        store_diagnostics::Bool=true)

    _validate_kutta_provider_interface(provider)
    _is_kutta_strategy(primary) || throw(ArgumentError(
        "PressureContinuityKutta primary strategy must satisfy the nonlinear-"*
        "strategy protocol (BroydenKutta or NewtonKrylovKutta); got $(typeof(primary))."))
    _is_kutta_strategy(fallback) || throw(ArgumentError(
        "PressureContinuityKutta fallback strategy must satisfy the nonlinear-"*
        "strategy protocol (BroydenKutta or NewtonKrylovKutta); got $(typeof(fallback))."))

    ps = _kutta_scale_value(:pressure_scale, pressure_scale)
    cs = _kutta_scale_value(:correction_scale, correction_scale)
    (isfinite(pressure_tolerance) && pressure_tolerance > 0) || throw(ArgumentError(
        "PressureContinuityKutta pressure_tolerance must be finite and positive; got $(pressure_tolerance)."))
    (isfinite(correction_tolerance) && correction_tolerance > 0) || throw(ArgumentError(
        "PressureContinuityKutta correction_tolerance must be finite and positive; got $(correction_tolerance)."))
    on_failure in (:error, :jump) || throw(ArgumentError(
        "PressureContinuityKutta on_failure must be :error or :jump; got $(repr(on_failure))."))

    return PressureContinuityKutta{typeof(provider), typeof(primary), typeof(fallback)}(
        provider, primary, fallback, ps, cs,
        Float64(pressure_tolerance), Float64(correction_tolerance),
        on_failure, store_diagnostics, KuttaDiagnostics{Float64}[])
end

"Convert an `:auto`/explicit scale argument to its internal Float64 encoding
(NaN = auto); explicit values must be finite and strictly positive."
function _kutta_scale_value(name::Symbol, scale::Union{Symbol, Real})
    if scale isa Symbol
        scale === :auto || throw(ArgumentError(
            "PressureContinuityKutta $(name) must be :auto or a finite positive value; got $(repr(scale))."))
        return NaN
    end
    (isfinite(scale) && scale > 0) || throw(ArgumentError(
        "PressureContinuityKutta $(name) must be :auto or a finite positive value; got $(scale)."))
    return Float64(scale)
end

"""
    kutta_diagnostics(kutta_closure)

Return the accepted per-step [`KuttaDiagnostics`](@ref) records in
physical-step order. Returns an empty collection for [`JumpKutta`](@ref) or
when the closure was constructed with `store_diagnostics=false` (failure
reporting through [`KuttaConvergenceError`](@ref) is unaffected).
"""
kutta_diagnostics(closure::PressureContinuityKutta) = closure.diagnostics
kutta_diagnostics(::JumpKutta) = KuttaDiagnostics{Float64}[]

################################################################################
# PRESSURE-PROVIDER INTERFACE
################################################################################

"""
    pressure_requirements(provider)

Tuple of `Symbol` requirements the runtime must be able to supply before trial
evaluation (e.g. `(:relative_surface_velocity,)`).
"""
function pressure_requirements end

"""
    evaluate_pressure!(pressure, provider, trial_view)

Evaluate one pressure per paired-edge panel into `pressure` (a 2×M matrix;
row 1 = upper panels, row 2 = lower panels) from the read-only completed
trial state `trial_view` (a [`KuttaTrialView`](@ref)). Must be deterministic
in the trial state and must mutate only `pressure` and provider scratch
explicitly designated as trial-local. A pressure gauge constant is permitted
only if it cancels in paired differences.
"""
function evaluate_pressure! end

"""
    commit_pressure!(provider, trial_view)

Called exactly once with the accepted trial. No-op for stateless providers.
"""
function commit_pressure! end

"""
    rollback_pressure!(provider)

Return provider trial state to the frozen pre-step state; called before every
trial evaluation and after every rejected or failed trial. No-op for
stateless providers.
"""
function rollback_pressure! end

function _validate_kutta_provider_interface(provider)
    for (f, sig) in ((pressure_requirements, Tuple{typeof(provider)}),
                     (evaluate_pressure!, Tuple{AbstractMatrix, typeof(provider), KuttaTrialView}),
                     (commit_pressure!, Tuple{typeof(provider), KuttaTrialView}),
                     (rollback_pressure!, Tuple{typeof(provider)}))
        hasmethod(f, sig) || throw(ArgumentError(
            "pressure provider $(typeof(provider)) does not implement $(nameof(f)) "*
            "for the required signature; see the pressure-provider interface docs."))
    end
    reqs = pressure_requirements(provider)
    unsupported = Tuple(r for r in reqs if r !== :relative_surface_velocity)
    isempty(unsupported) || throw(ArgumentError(
        "pressure provider $(typeof(provider)) declares requirements "*
        "$(unsupported) that this runtime cannot supply; supported: "*
        "(:relative_surface_velocity,)."))
    return nothing
end

"Read-only view of one completed trial handed to pressure providers: the body,
the paired-edge map, and the completed exterior relative surface velocity
(3×ncells, kinematics already subtracted)."
struct KuttaTrialView{TF<:AbstractFloat, TB<:RigidWakeBody}
    body::TB
    edges::SheddingEdgeMap
    velocity::Matrix{TF}
end

"""
    SteadyBernoulliProvider(rho)

Steady-Bernoulli pressure provider for [`PressureContinuityKutta`](@ref):
gauge pressure `p = −½ρ‖u_t‖²` at panel centroids from the tangential
projection of the completed relative exterior surface velocity (the common
reference `p_ref + ½ρU_ref²` cancels exactly in paired differences). Requires
only `:relative_surface_velocity`; stateless (commit and rollback are
no-ops). Does not wrap or reuse the `PressureBernoulli` monitor.
"""
struct SteadyBernoulliProvider{T<:Real} <: AbstractKuttaPressureProvider
    rho::T

    function SteadyBernoulliProvider(rho::T) where {T<:Real}
        (isfinite(rho) && rho > 0) || throw(ArgumentError(
            "SteadyBernoulliProvider density must be finite and positive; got $(rho)."))
        return new{T}(rho)
    end
end

pressure_requirements(::SteadyBernoulliProvider) = (:relative_surface_velocity,)
commit_pressure!(::SteadyBernoulliProvider, trial_view) = nothing
rollback_pressure!(::SteadyBernoulliProvider) = nothing

function evaluate_pressure!(pressure::AbstractMatrix,
        provider::SteadyBernoulliProvider, trial::KuttaTrialView)
    body = trial.body
    edges = trial.edges
    U = trial.velocity
    normals = body.normals
    rho = provider.rho
    for e in eachindex(edges.upper)
        pressure[1, e] = _steady_bernoulli_gauge(U, normals, edges.upper[e], rho)
        pressure[2, e] = _steady_bernoulli_gauge(U, normals, edges.lower[e], rho)
    end
    return pressure
end

"Gauge pressure −½ρ‖u_t‖² at panel `i` from the tangential projection of the
relative surface velocity (matches the steady `PressureBernoulli` velocity
convention without reusing the monitor)."
@inline function _steady_bernoulli_gauge(U::AbstractMatrix,
        normals::AbstractMatrix, i::Int, rho)
    un = U[1, i]*normals[1, i] + U[2, i]*normals[2, i] + U[3, i]*normals[3, i]
    ut2 = (U[1, i] - un*normals[1, i])^2 +
          (U[2, i] - un*normals[2, i])^2 +
          (U[3, i] - un*normals[3, i])^2
    return -0.5*rho*ut2
end

################################################################################
# CONFIGURATION VALIDATION
################################################################################

const _KUTTA_DOMAIN_MSG = "the Phase 3 pressure-continuity/TE-anchored runtime "*
    "supports exactly one fully paired finite-wake Dirichlet RigidWakeBody with "*
    "one PanelWake, VelocityThroughSources, a Backslash solver, and "*
    "bound_strength_rlx == 1; all other configurations run only through the "*
    "default RigidTransitionAttachment + JumpKutta path"

"Is this the exact legacy default pair that must branch into the pre-Phase-3
call sequence before any new state is allocated?"
_is_legacy_kutta(wake_attachment::AbstractWakeAttachment,
    kutta_closure::AbstractKuttaClosure) =
    wake_attachment isa RigidTransitionAttachment && kutta_closure isa JumpKutta

"""
    _validate_kutta_configuration(entry, systems_tuple, wakes_tuple,
                                  body_solvers, formulation, backend_system,
                                  wake_attachment, kutta_closure; options...)

Validate a non-default attachment/closure configuration against the Phase 3
support boundary before any body, wake, provider, or solver state is mutated.
`entry` is `:simulate` or `:steady`. Returns `(body, solver, wake)` where
`wake` is `nothing` for `steady!`.
"""
function _validate_kutta_configuration(entry::Symbol, systems_tuple::Tuple,
        wakes_tuple::Tuple, body_solvers, formulation,
        backend_system,
        wake_attachment::AbstractWakeAttachment,
        kutta_closure::AbstractKuttaClosure;
        bound_strength_rlx::Real=1,
        set_Das_eta_kinematic::Real=NaN,
        set_Das_eta_freestream::Real=NaN,
        set_Das_min_kinematic_displacement::Real=0.0,
        set_Das_kinematic_arc::Bool=true,
        set_Das_refresh::Bool=false)

    # entry-point boundary
    if entry === :steady && wake_attachment isa TEAnchoredAttachment
        throw(ArgumentError("steady! rejects TEAnchoredAttachment: a steady "*
            "solve cannot infer the convected downstream row that defines the "*
            "Route B live panel. Use RigidTransitionAttachment, or simulate!."))
    end

    # single Dirichlet paired RigidWakeBody
    length(systems_tuple) == 1 || throw(ArgumentError(
        "unsupported configuration: $(length(systems_tuple)) bodies; $(_KUTTA_DOMAIN_MSG)."))
    body = systems_tuple[1]
    body isa RigidWakeBody{<:Any, 2, <:Any, true} || throw(ArgumentError(
        "unsupported body $(typeof(body)): $(_KUTTA_DOMAIN_MSG)."))
    body.semiinfinite_wake && throw(ArgumentError(
        "unsupported semiinfinite_wake=true body: $(_KUTTA_DOMAIN_MSG)."))
    body.nsheddings > 0 || throw(ArgumentError(
        "unsupported body without shedding edges: $(_KUTTA_DOMAIN_MSG)."))
    for (s, shed) in enumerate(body.shedding)
        for i in axes(shed, 2)
            shed[4, i] > 0 || throw(ArgumentError(
                "unsupported unpaired shedding edge $(i) on surface $(s) "*
                "(lower panel index $(shed[4, i])): every shedding edge must "*
                "have one valid upper/lower panel pair."))
        end
    end

    # formulation
    formulation isa VelocityThroughSources || throw(ArgumentError(
        "unsupported formulation $(typeof(formulation)): non-default "*
        "wake_attachment/kutta_closure requires VelocityThroughSources (the "*
        "correction channel would collide with other formulations)."))

    # solver (first implementation: Backslash only, per user amendment)
    solver = _single_body_solver(body_solvers)
    if solver isa FGSSolver || solver isa KrylovSolver
        throw(ArgumentError(
            "unsupported solver $(typeof(solver)): KrylovSolver/FGSSolver "*
            "support for non-default wake_attachment/kutta_closure is "*
            "deferred until the Backslash implementation is validated; use a "*
            "Backslash solver."))
    end
    solver isa Backslash || throw(ArgumentError(
        "unsupported solver $(typeof(solver)): $(_KUTTA_DOMAIN_MSG)."))

    # relaxation invalidates both closures' post-solve consistency
    bound_strength_rlx == 1 || throw(ArgumentError(
        "unsupported bound_strength_rlx=$(bound_strength_rlx): post-solve "*
        "strength relaxation changes Cμ after the coupled closure; non-default "*
        "wake_attachment/kutta_closure requires bound_strength_rlx == 1."))

    # wake boundary
    wake = nothing
    if entry === :simulate
        active = Tuple(w for w in wakes_tuple if !isnothing(w))
        length(active) == 1 || throw(ArgumentError(
            "unsupported wake count $(length(active)): $(_KUTTA_DOMAIN_MSG)."))
        wake = active[1]
        wake isa PanelWake || throw(ArgumentError(
            "unsupported wake $(typeof(wake)): $(_KUTTA_DOMAIN_MSG) "*
            "(PanelParticleWake remains legacy-only)."))
    end

    # Route A assembles the attachment operator W and the Backslash
    # factorization exactly once per run (the §8.1 per-geometry-version
    # refresh is not implemented); a per-step Das refresh would silently
    # solve against stale operators.
    if wake_attachment isa RigidTransitionAttachment && set_Das_refresh
        throw(ArgumentError(
            "unsupported set_Das_refresh=true with a non-default "*
            "wake_attachment/kutta_closure: the runtime assembles the "*
            "attachment operator and factorization once per run, so a "*
            "per-step Das refresh would silently solve with stale operators; "*
            "deferred until per-geometry-version reassembly is implemented."))
    end

    # Route B rejects every set_Das_* option (it never reads Das as user input)
    if wake_attachment isa TEAnchoredAttachment
        isnan(set_Das_eta_kinematic) || throw(ArgumentError(
            "TEAnchoredAttachment rejects set_Das_eta_kinematic=$(set_Das_eta_kinematic): "*
            "Route B owns its live-panel geometry and reads no user Das."))
        isnan(set_Das_eta_freestream) || throw(ArgumentError(
            "TEAnchoredAttachment rejects set_Das_eta_freestream=$(set_Das_eta_freestream): "*
            "Route B owns its live-panel geometry and reads no user Das."))
        set_Das_min_kinematic_displacement == 0 || throw(ArgumentError(
            "TEAnchoredAttachment rejects set_Das_min_kinematic_displacement="*
            "$(set_Das_min_kinematic_displacement): Route B owns its live-panel "*
            "geometry and reads no user Das."))
        set_Das_kinematic_arc || throw(ArgumentError(
            "TEAnchoredAttachment rejects non-default set_Das_kinematic_arc=false: "*
            "Route B owns its live-panel geometry and reads no user Das."))
        !set_Das_refresh || throw(ArgumentError(
            "TEAnchoredAttachment rejects set_Das_refresh=true: Route B owns "*
            "its live-panel geometry and reads no user Das."))
    end

    return body, solver, wake
end

################################################################################
# RUNTIME STATE
################################################################################

"Restorable snapshot of every value a rejected or repeated trial could mutate
(architecture §7.1), captured once per physical step after the frozen wake-
influence pass."
mutable struct KuttaSnapshot{TF<:AbstractFloat}
    strength::Matrix{TF}
    potential::Vector{TF}
    velocity_frozen::Matrix{TF}       # freestream + kinematic + wake pass
    Das::Vector{Matrix{TF}}
    wake_correction::Vector{Vector{TF}}
    wake_shift::Vector{TF}
    wake_correction_active::Bool
    wake_nodes::Vector{Array{TF, 3}}
    wake_strength::Vector{Array{TF, 3}}
    wake_velocity::Vector{Array{TF, 3}}
    nwakes::Int
    overflowed::Bool
    live_rows::Int
    live_step_id::Int
end

function KuttaSnapshot{TF}(body::RigidWakeBody, wake) where {TF<:AbstractFloat}
    nwake_surfs = isnothing(wake) ? 0 : length(wake.nodes)
    return KuttaSnapshot{TF}(
        similar(body.strength), similar(body.potential),
        similar(body.velocity),
        [similar(D) for D in body.Das],
        [similar(c) for c in body.wake_strength_correction],
        similar(body.wake_strength_shift), false,
        [isnothing(wake) ? Array{TF, 3}(undef, 0, 0, 0) : similar(wake.nodes[i]) for i in 1:nwake_surfs],
        [isnothing(wake) ? Array{TF, 3}(undef, 0, 0, 0) : similar(wake.strength[i]) for i in 1:nwake_surfs],
        [isnothing(wake) ? Array{TF, 3}(undef, 0, 0, 0) : similar(wake.velocity[i]) for i in 1:nwake_surfs],
        0, false, 0, -1)
end

"Restorable record of one completed nonlinear trial (architecture §7.2). The
accepted record is committed without re-running the trial."
mutable struct KuttaTrialRecord{TF<:AbstractFloat}
    c::Vector{TF}
    mu::Vector{TF}            # solved doublet strengths
    residual::Vector{TF}      # dimensional paired pressure residual
    r_inf::TF
    r_l2_edge::TF
    finite::Bool
    inner_status::Symbol
    valid::Bool               # has this record been filled this step?
end

KuttaTrialRecord{TF}(N::Int, M::Int) where {TF<:AbstractFloat} =
    KuttaTrialRecord{TF}(zeros(TF, M), zeros(TF, N), zeros(TF, M),
        TF(NaN), TF(NaN), false, :not_run, false)

function _copy_record!(dst::KuttaTrialRecord, src::KuttaTrialRecord)
    dst.c .= src.c
    dst.mu .= src.mu
    dst.residual .= src.residual
    dst.r_inf = src.r_inf
    dst.r_l2_edge = src.r_l2_edge
    dst.finite = src.finite
    dst.inner_status = src.inner_status
    dst.valid = src.valid
    return dst
end

"""
Internal runtime for one non-default wake-attachment/Kutta-closure
configuration (single body, single optional `PanelWake`, `Backslash` solver).
Initialized only when the configuration differs from the legacy default.
"""
mutable struct KuttaRuntime{TF<:AbstractFloat, TB<:RigidWakeBody,
        TS<:Backslash, TW, C<:AbstractKuttaClosure, TGLU}
    closure::C
    route::Symbol             # :A | :B
    entry::Symbol             # :simulate | :steady
    body::TB
    solver::TS
    wake::TW                  # PanelWake or Nothing
    edges::SheddingEdgeMap
    # attachment-aware operator
    W::Matrix{TF}             # N×M attached/live-panel RHS potential columns
    Glu::TGLU                 # active factorization (route A: solver's; route B: runtime-owned)
    Gb::Matrix{TF}            # route B assembly buffer (0×0 for route A)
    geometry_version::Int
    live_active::Bool         # route B: live panel participates in the solve
    # per-step frozen data
    edge_lengths::Vector{TF}
    U_s::TF
    L_s::TF
    P_s::TF
    C_s::TF
    # nonlinear state
    c::Vector{TF}             # accepted correction (warm start across steps)
    snapshot::KuttaSnapshot{TF}
    trial::KuttaTrialRecord{TF}
    best::KuttaTrialRecord{TF}
    accepted::KuttaTrialRecord{TF}
    current::KuttaTrialRecord{TF}    # last accepted nonlinear iterate
    Jinv::Matrix{TF}          # dense M×M inverse-Broyden approximation
    pressure::Matrix{TF}      # 2×M paired centroid pressure workspace
    work_N::Vector{TF}        # length-N workspace (holds W·c)
    x::Vector{TF}             # scaled correction iterate
    dx::Vector{TF}
    r_scaled::Vector{TF}
    r_prev::Vector{TF}
    work_M::Vector{TF}
    work_M2::Vector{TF}
    # counters (reset per step)
    body_solves::Int
    backtracks::Int
    restarts::Int
    outer_iterations::Int
    fallback_entered::Bool
    method::Symbol
    accepted_dc_inf_scaled::TF   # last accepted scaled correction step
end

"Build the Kutta runtime for a validated non-default configuration. Route A
assembles its attachment operator once (rigid-geometry invariance, matching
`Backslash` reuse); Route B assembles per step."
function _initialize_kutta(entry::Symbol, body::RigidWakeBody, solver::Backslash,
        wake, wake_attachment::AbstractWakeAttachment,
        kutta_closure::AbstractKuttaClosure)
    TF = eltype(body.strength)
    N = body.ncells
    edges = _shedding_edge_map(body)
    M = length(edges.upper)
    route = wake_attachment isa TEAnchoredAttachment ? :B : :A

    W = zeros(TF, N, M)
    Gb = route === :B ? Matrix{TF}(undef, N, N) : Matrix{TF}(undef, 0, 0)
    if route === :A
        _assemble_W!(W, body, edges)
        Glu = solver.Glu
    else
        # placeholder factorization; replaced at every Route B geometry update
        Gb .= zero(TF)
        _G!(Gb, body, body; core_size=body.core_size_panel,
            update_geometry=false)
        Glu = lu!(Gb)
    end

    rt = KuttaRuntime{TF, typeof(body), typeof(solver), typeof(wake),
            typeof(kutta_closure), typeof(Glu)}(
        kutta_closure, route, entry, body, solver, wake, edges,
        W, Glu, Gb, 0, false,
        zeros(TF, M), TF(NaN), TF(NaN), TF(NaN), TF(NaN),
        zeros(TF, M),
        KuttaSnapshot{TF}(body, wake),
        KuttaTrialRecord{TF}(N, M), KuttaTrialRecord{TF}(N, M),
        KuttaTrialRecord{TF}(N, M), KuttaTrialRecord{TF}(N, M),
        Matrix{TF}(LA.I, M, M),
        zeros(TF, 2, M), zeros(TF, N),
        zeros(TF, M), zeros(TF, M), zeros(TF, M), zeros(TF, M),
        zeros(TF, M), zeros(TF, M),
        0, 0, 0, 0, false, :none, TF(NaN))

    _kutta_update_edge_lengths!(rt)

    # a warm start restores the committed correction onto the body before the
    # first resumed solve; seed the nonlinear warm-start coordinate from it
    if body.wake_correction_active[]
        k = 0
        for corr in body.wake_strength_correction
            for v in corr
                k += 1
                rt.c[k] = v
            end
        end
    end
    return rt
end

"TE edge length ℓ_e of each paired shedding edge (frozen per trial sequence;
used for the edge-weighted L2 norm and the automatic correction scale)."
function _kutta_update_edge_lengths!(rt::KuttaRuntime)
    body = rt.body
    for e in eachindex(rt.edges.upper)
        s = rt.edges.isurf[e]
        i = rt.edges.iedge[e]
        shed = body.shedding[s]
        panel = shed[1, i]
        na = shed[2, i]
        nb = shed[3, i]
        va = view(body.nodes, :, body.cells[na, panel])
        vb = view(body.nodes, :, body.cells[nb, panel])
        rt.edge_lengths[e] = sqrt((va[1]-vb[1])^2 + (va[2]-vb[2])^2 + (va[3]-vb[3])^2)
    end
    return rt
end

################################################################################
# SNAPSHOT CAPTURE/RESTORE
################################################################################

function _kutta_capture_snapshot!(rt::KuttaRuntime)
    body = rt.body
    snap = rt.snapshot
    snap.strength .= body.strength
    snap.potential .= body.potential
    snap.velocity_frozen .= body.velocity
    for (D, Ds) in zip(body.Das, snap.Das)
        Ds .= D
    end
    for (c, cs) in zip(body.wake_strength_correction, snap.wake_correction)
        cs .= c
    end
    snap.wake_shift .= body.wake_strength_shift
    snap.wake_correction_active = body.wake_correction_active[]
    wake = rt.wake
    if !isnothing(wake)
        for i in eachindex(wake.nodes)
            snap.wake_nodes[i] .= wake.nodes[i]
            snap.wake_strength[i] .= wake.strength[i]
            snap.wake_velocity[i] .= wake.velocity[i]
        end
        snap.nwakes = wake.nwakes[]
        snap.overflowed = wake.overflowed[]
        snap.live_rows = wake.live_rows[]
        snap.live_step_id = wake.live_step_id[]
    end
    return rt
end

function _kutta_restore_snapshot!(rt::KuttaRuntime)
    body = rt.body
    snap = rt.snapshot
    body.strength .= snap.strength
    body.potential .= snap.potential
    body.velocity .= snap.velocity_frozen
    for (D, Ds) in zip(body.Das, snap.Das)
        D .= Ds
    end
    for (c, cs) in zip(body.wake_strength_correction, snap.wake_correction)
        c .= cs
    end
    body.wake_strength_shift .= snap.wake_shift
    body.wake_correction_active[] = snap.wake_correction_active
    wake = rt.wake
    if !isnothing(wake)
        for i in eachindex(wake.nodes)
            wake.nodes[i] .= snap.wake_nodes[i]
            wake.strength[i] .= snap.wake_strength[i]
            wake.velocity[i] .= snap.wake_velocity[i]
        end
        wake.nwakes[] = snap.nwakes
        wake.overflowed[] = snap.overflowed
        wake.live_rows[] = snap.live_rows
        wake.live_step_id[] = snap.live_step_id
    end
    rollback_pressure!(_kutta_provider(rt))
    return rt
end

_kutta_provider(rt::KuttaRuntime) =
    rt.closure isa PressureContinuityKutta ? rt.closure.provider : _NoProvider()

struct _NoProvider <: AbstractKuttaPressureProvider end
pressure_requirements(::_NoProvider) = ()
evaluate_pressure!(pressure::AbstractMatrix, ::_NoProvider, trial_view) = pressure
commit_pressure!(::_NoProvider, trial_view) = nothing
rollback_pressure!(::_NoProvider) = nothing

################################################################################
# TRIAL TRANSACTION
################################################################################

"""
One complete residual transaction (architecture §7.2): restore the frozen
snapshot, apply the trial correction `c` through the RHS channel, solve the
coupled linear body problem on the cached factorization, reconstruct the
complete exterior relative surface velocity once, evaluate pressure, and form
the paired residual into `rt.trial`. Never touches wake topology, provider
history, or accepted state.

`suppress_attached` runs the whole transaction with the attached/live wake
suppressed (Route B cold startup). `reconstruct=true, evaluate=false` completes
the exterior velocity without a pressure evaluation (used to freeze the
automatic scales from the initial clean trial).
"""
function _kutta_trial!(rt::KuttaRuntime{TF}, c::AbstractVector;
        backend_solve, backend_system, grad_mu_options,
        reconstruct::Bool=true, evaluate::Bool=true,
        suppress_attached::Bool=false) where {TF}

    body = rt.body
    trial = rt.trial
    trial.valid = false
    trial.c .= c

    _kutta_restore_snapshot!(rt)

    # ----- coupled linear solve: operator linear, correction on the RHS -----
    body.wake_correction_active[] = false
    body.suppress_attached_wake[] = suppress_attached
    inner_ok = false
    try
        _set_core_sizes!((body,), :core_size_panel)
        set_strengths!(body)                 # σ = −u·n from frozen velocity; μ = 0
        body.potential .= zero(TF)
        influence!(body, body, backend_solve; scalar_potential=true,
            velocity=false)
        if !suppress_attached
            LA.mul!(rt.work_N, rt.W, trial.c)
            body.potential .-= rt.work_N     # rhs = −S·σ + W·c
        end
        _kutta_ldiv!(rt, body)
        rt.body_solves += 1
        inner_ok = all(isfinite, view(body.strength, :, 2))
    finally
        body.suppress_attached_wake[] = false
    end
    trial.inner_status = inner_ok ? :converged : :nonfinite
    trial.mu .= view(body.strength, :, 2)
    # the trial's potential is scratch; physical potential returns with the
    # next restore (mirrors the production solve!'s save/restore semantics)
    body.potential .= rt.snapshot.potential

    if !inner_ok
        trial.finite = false
        trial.r_inf = TF(NaN)
        trial.r_l2_edge = TF(NaN)
        trial.valid = true
        return trial
    end

    if !reconstruct && !evaluate     # jump solves need neither
        trial.finite = true
        trial.residual .= zero(TF)
        trial.r_inf = zero(TF)
        trial.r_l2_edge = zero(TF)
        trial.valid = true
        return trial
    end

    # ----- one complete exterior velocity reconstruction -----
    _kutta_reconstruct_body_velocity!(rt, trial.c; backend_system,
        grad_mu_options, suppress_attached)

    if !evaluate                     # scale-freezing pre-pass
        trial.finite = true
        trial.residual .= zero(TF)
        trial.r_inf = zero(TF)
        trial.r_l2_edge = zero(TF)
        trial.valid = true
        return trial
    end

    # ----- pressure evaluation and paired residual -----
    provider = _kutta_provider(rt)
    rollback_pressure!(provider)
    trial_view = KuttaTrialView{TF, typeof(body)}(body, rt.edges, body.velocity)
    evaluate_pressure!(rt.pressure, provider, trial_view)

    finite = true
    for e in eachindex(rt.edges.upper)
        r = rt.pressure[1, e] - rt.pressure[2, e]
        trial.residual[e] = r
        finite &= isfinite(r)
    end
    trial.finite = finite
    trial.r_inf = finite ? maximum(abs, trial.residual) : TF(NaN)
    trial.r_l2_edge = finite ? _edge_weighted_l2(trial.residual, rt.edge_lengths) : TF(NaN)
    trial.valid = true
    return trial
end

"Solve the assembled Dirichlet system on the runtime's active factorization
(`rhs = −body.potential`, matching the production `_solve!` convention)."
function _kutta_ldiv!(rt::KuttaRuntime, body::RigidWakeBody)
    rhs = rt.solver.rhs
    rhs .= .-body.potential
    ldiv!(view(body.strength, :, 2), rt.Glu, rhs)
    return nothing
end

"Reconstruct the complete exterior relative surface velocity on the body only
(wake probes accumulate exactly once per step, at commit): body-on-body
influence through any backend with the correction INACTIVE (the proportional
γ = Cμ part flows through both the direct and FMM paths), then the exact
direct affine add-on for the trial correction, plus the +½∇μ tangential
half-jump. The add-on must use `core_size_panel`: the self-pair
conditioning rule evaluates body-on-body sources — including the attached
transition strips whose proportional part this add-on completes — at
`core_size_panel`, so both halves of one strip strength must share that
regularization (wake probes are cross pairs and use `core_size_targets`)."
function _kutta_reconstruct_body_velocity!(rt::KuttaRuntime, c::AbstractVector;
        backend_system, grad_mu_options, suppress_attached::Bool=false)
    body = rt.body
    body.velocity .= rt.snapshot.velocity_frozen
    body.suppress_attached_wake[] = suppress_attached
    was_active = _operator_mode_begin!(body)
    try
        _set_core_sizes!((body,), :core_size_targets)
        influence!((body,), (body,), backend_system; precalc=false,
            scalar_potential=false,
            velocity=true,
            velocity_gradient=(false,),
            direct_conditioning=_self_panel_core_size_conditioning())
    finally
        body.suppress_attached_wake[] = false
        _operator_mode_end!(body, was_active)
    end
    suppress_attached || _add_affine_attached_velocity!(body.velocity,
        body.controlpoints, body, c; core_size=body.core_size_panel)
    if has_grad_mu(body)
        compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
            body.cells, body.neighbor,
            view(body.strength, :, get_Gammai(body)),
            _bound_surface_vorticity_te_info(body);
            scale=0.5,
            nodes=body.nodes,
            grad_mu_options=grad_mu_options)
    end
    return nothing
end

################################################################################
# DIRECT AFFINE ATTACHED-WAKE ADD-ON
################################################################################

"""
Velocity induced at `target` by panel `i_panel`'s attached transition strip at
UNIT doublet strength: the exact finite-wake construction of the index-path
`_induced_wake` (two triangles from the panel's TE nodes and `Das` columns via
`shedding_full`), at an explicit `core_size`. Returns a zero vector for a
panel with no shedding edge.
"""
function _unit_attached_strip_velocity(target, body::RigidWakeBody{<:Any,<:Any,TF},
        i_panel::Int, TK, switch, core_size) where {TF}
    idx_1 = body.shedding_full[1, i_panel]
    idx_1 > 0 || return zero(SVector{3,TF})

    v1n = body.cells[idx_1, i_panel]
    v1 = SVector{3,TF}(body.nodes[1, v1n], body.nodes[2, v1n], body.nodes[3, v1n])
    idx_2 = body.shedding_full[2, i_panel]
    v2n = body.cells[idx_2, i_panel]
    v2 = SVector{3,TF}(body.nodes[1, v2n], body.nodes[2, v2n], body.nodes[3, v2n])

    i_surf = body.shedding_full[3, i_panel]
    das_col_1 = body.shedding_full[5, i_panel]
    das_col_2 = body.shedding_full[6, i_panel]
    Das = body.Das[i_surf]
    vw1 = SVector{3,TF}(v1[1] + Das[1, das_col_1], v1[2] + Das[2, das_col_1],
        v1[3] + Das[3, das_col_1])
    vw2 = SVector{3,TF}(v2[1] + Das[1, das_col_2], v2[2] + Das[2, das_col_2],
        v2[3] + Das[3, das_col_2])

    strength_vec = SVector{1,TF}(one(TF))

    control_point = (v1 + v2 + vw1) * TF(0.333333333333333)
    R, _ = rotate_to_panel(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3],
        vw1[1], vw1[2], vw1[3])
    _, u, _ = _induced(target, (v1, v2, vw1), control_point, strength_vec, TK,
        core_size, R, switch)

    control_point = (vw1 + v2 + vw2) * TF(0.333333333333333)
    R, _ = rotate_to_panel(vw1[1], vw1[2], vw1[3], v2[1], v2[2], v2[3],
        vw2[1], vw2[2], vw2[3])
    _, du, _ = _induced(target, (vw1, v2, vw2), control_point, strength_vec, TK,
        core_size, R, switch)

    return u + du
end

"Per-panel affine strip weights for a correction vector `c` — the same
decomposition `set_wake_correction!` installs (upper −c/2, lower +c/2 per
paired edge), as a `(panel, weight)` list over nonzero entries."
function _kutta_affine_strips(body::RigidWakeBody{<:Any,<:Any,TF},
        c::AbstractVector) where {TF}
    strips = Tuple{Int, TF}[]
    k = 0
    for shed in body.shedding
        for i in axes(shed, 2)
            k += 1
            ce = TF(c[k])
            iszero(ce) && continue
            lo = shed[4, i]
            @assert lo > 0 "affine attached-wake add-on requires paired edges"
            push!(strips, (shed[1, i], -ce/2))
            push!(strips, (lo, ce/2))
        end
    end
    return strips
end

"""
    _add_affine_attached_velocity!(velocity, targets, body, c; core_size)

Accumulate the affine attached-wake velocity of the correction `c` — the
γ = Cμ − c shift carried by the M trailing-edge transition strips (one wake
row of panels) — onto `velocity[:, j]` at the points `targets[:, j]`. This is
the direct, backend-independent evaluation of the term the index-path
`wake_strength_shift` supplies (and the FMM buffer path does not): the main
reconstruction influence runs with the correction inactive through any
backend, and this O(M×N_targets) add-on completes it exactly.
"""
function _add_affine_attached_velocity!(velocity::AbstractMatrix,
        targets::AbstractMatrix, body::RigidWakeBody{<:Any,<:Any,TF},
        c::AbstractVector; core_size=body.core_size_targets) where {TF}
    iszero(c) && return velocity
    strips = _kutta_affine_strips(body, c)
    isempty(strips) && return velocity
    TK = get_wake_kernel(body)
    switch = FastMultipole.DerivativesSwitch(false, true, false)
    Threads.@threads for j in axes(targets, 2)
        target = SVector{3,TF}(targets[1, j], targets[2, j], targets[3, j])
        u = zero(SVector{3,TF})
        for (i_panel, w) in strips
            u += w*_unit_attached_strip_velocity(target, body, i_panel, TK,
                switch, core_size)
        end
        velocity[1, j] += u[1]
        velocity[2, j] += u[2]
        velocity[3, j] += u[3]
    end
    return velocity
end

"Affine attached-wake add-on at the wake probe nodes (rows 1:nwakes+1),
accumulated into `wake.velocity` — the same targets the body→wake-probes
influence pass feeds."
function _add_affine_attached_velocity!(wake::PanelWake,
        body::RigidWakeBody{<:Any,<:Any,TF}, c::AbstractVector;
        core_size=body.core_size_targets) where {TF}
    iszero(c) && return wake
    nrows = wake.nwakes[] + 1
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        vel = wake.velocity[i_surf]
        _add_affine_attached_velocity!(
            reshape(view(vel, :, 1:nrows, :), 3, :),
            reshape(view(nodes, :, 1:nrows, :), 3, :),
            body, c; core_size)
    end
    return wake
end

"Edge-length-weighted pressure L2 norm √(Σℓ_e r_e² / Σℓ_e)."
function _edge_weighted_l2(r::AbstractVector, ℓ::AbstractVector)
    num = zero(eltype(r))
    den = zero(eltype(ℓ))
    for e in eachindex(r)
        num += ℓ[e]*r[e]^2
        den += ℓ[e]
    end
    return sqrt(num/den)
end

################################################################################
# FROZEN SCALES
################################################################################

"""
Freeze the pressure and correction scales for this physical step from the
initial clean finite trial (architecture §9.1): `U_s` is the RMS of the
completed relative surface speed over all body centroids (a deterministic
serial reduction), `L_s` the median characteristic length of the paired TE
panels, `P_s = ½ρU_s²`, `C_s = U_s·L_s`. Explicit user scales are honored
and equally frozen. Degenerate automatic scales throw an `ArgumentError`
requesting explicit positive scales.
"""
function _kutta_freeze_scales!(rt::KuttaRuntime{TF}) where {TF}
    closure = rt.closure::PressureContinuityKutta
    body = rt.body

    # deterministic serial RMS of the relative surface speed
    N = body.ncells
    acc = zero(TF)
    for i in 1:N
        acc += body.velocity[1, i]^2 + body.velocity[2, i]^2 + body.velocity[3, i]^2
    end
    rt.U_s = sqrt(acc/N)

    areas = _panel_areas(body)
    lengths = TF[]
    for e in eachindex(rt.edges.upper)
        push!(lengths, TF(sqrt(areas[rt.edges.upper[e]])))
        push!(lengths, TF(sqrt(areas[rt.edges.lower[e]])))
    end
    sort!(lengths)
    n = length(lengths)
    rt.L_s = isodd(n) ? lengths[(n+1)÷2] : (lengths[n÷2] + lengths[n÷2+1])/2

    rho = closure.provider isa SteadyBernoulliProvider ? closure.provider.rho : one(TF)

    if isnan(closure.pressure_scale)
        P_s = TF(0.5)*TF(rho)*rt.U_s^2
        (isfinite(P_s) && P_s > 0) || throw(ArgumentError(
            "degenerate automatic pressure scale P_s=$(P_s) (U_s=$(rt.U_s)): "*
            "supply an explicit positive pressure_scale to PressureContinuityKutta."))
        rt.P_s = P_s
    else
        rt.P_s = TF(closure.pressure_scale)
    end

    if isnan(closure.correction_scale)
        C_s = rt.U_s*rt.L_s
        (isfinite(C_s) && C_s > 0) || throw(ArgumentError(
            "degenerate automatic correction scale C_s=$(C_s) (U_s=$(rt.U_s), "*
            "L_s=$(rt.L_s)): supply an explicit positive correction_scale to "*
            "PressureContinuityKutta."))
        rt.C_s = C_s
    else
        rt.C_s = TF(closure.correction_scale)
    end
    return rt
end

################################################################################
# NONLINEAR DRIVERS
################################################################################

"Do all three §9.4 acceptance gates pass for the current trial with scaled
step `dx_inf_scaled`?"
function _kutta_gates_pass(rt::KuttaRuntime, trial::KuttaTrialRecord,
        dx_inf_scaled::Real)
    closure = rt.closure::PressureContinuityKutta
    return trial.finite &&
        trial.inner_status === :converged &&
        trial.r_inf/rt.P_s <= closure.pressure_tolerance &&
        dx_inf_scaled <= closure.correction_tolerance
end

"Track the best finite restorable trial."
function _kutta_track_best!(rt::KuttaRuntime, trial::KuttaTrialRecord)
    if trial.finite && trial.inner_status === :converged &&
            (!rt.best.valid || !rt.best.finite || trial.r_inf < rt.best.r_inf)
        _copy_record!(rt.best, trial)
    end
    return rt
end

"""
Safeguarded dense inverse-Broyden primary driver (architecture §9.2) on the
scaled coordinates x = c/C_s, r̂ = r/P_s. Returns `:converged`,
`:fallback` (recoverable stagnation → try Newton–Krylov from the best
trial), or `:failed`.
"""
function _kutta_broyden!(rt::KuttaRuntime{TF}, trial_fun!) where {TF}
    strategy = (rt.closure::PressureContinuityKutta).primary::BroydenKutta
    M = length(rt.x)
    rt.method = :broyden

    # initial clean trial at the warm-started correction
    rt.x .= rt.c ./ rt.C_s
    trial = trial_fun!(rt.x)
    _kutta_track_best!(rt, trial)
    if !(trial.finite && trial.inner_status === :converged)
        return :failed        # the clean trial itself is unusable
    end
    rt.r_scaled .= trial.residual ./ rt.P_s
    _copy_record!(rt.current, trial)
    # a defined zero step satisfies the correction gate only if the other
    # gates pass on the clean trial
    if _kutta_gates_pass(rt, trial, zero(TF))
        rt.accepted_dc_inf_scaled = zero(TF)
        _copy_record!(rt.accepted, trial)
        return :converged
    end

    # conservative initial inverse-Jacobian scale from one secant probe along
    # the uniform direction (Phase 1 §6 recommendation): captures the sign and
    # magnitude of the dominant diagonal so the first proposals descend
    jinv_scale = one(TF)
    let delta = TF(1e-2)
        rt.work_M .= rt.x .+ delta/sqrt(TF(M))
        probe = trial_fun!(rt.work_M)
        if probe.finite && probe.inner_status === :converged
            _kutta_track_best!(rt, probe)
            s = zero(TF)
            for e in 1:M
                s += (probe.residual[e]/rt.P_s - rt.r_scaled[e])/sqrt(TF(M))
            end
            s /= delta
            if isfinite(s) && abs(s) > sqrt(eps(TF))
                jinv_scale = one(TF)/s
            end
        end
    end
    rt.Jinv .= jinv_scale .* Matrix{TF}(LA.I, M, M)
    restarts_left = strategy.max_restarts

    for iter in 1:strategy.max_iterations
        rt.outer_iterations += 1

        # proposed scaled step dx = −Jinv·r̂, capped in the ∞-norm
        LA.mul!(rt.dx, rt.Jinv, rt.r_scaled)
        rt.dx .*= -one(TF)
        dx_inf = maximum(abs, rt.dx)
        if dx_inf > strategy.max_scaled_step
            rt.dx .*= strategy.max_scaled_step/dx_inf
            dx_inf = strategy.max_scaled_step
        end

        # the current iterate already satisfies the pressure and inner gates
        # and the proposed step is below the correction tolerance: accept the
        # current iterate (no further step can improve within tolerance)
        if _kutta_gates_pass(rt, rt.current, dx_inf)
            rt.accepted_dc_inf_scaled = dx_inf
            _copy_record!(rt.accepted, rt.current)
            return :converged
        end

        # Armijo backtracking on ‖r̂‖₂
        r0 = LA.norm(rt.r_scaled)
        lambda = one(TF)
        accepted = false
        n_back = 0
        while true
            rt.work_M .= rt.x .+ lambda .* rt.dx
            trial = trial_fun!(rt.work_M)
            usable = trial.finite && trial.inner_status === :converged
            if usable
                _kutta_track_best!(rt, trial)
                rnew = LA.norm(trial.residual ./ rt.P_s)
                if rnew <= (1 - TF(1e-4)*lambda)*r0
                    accepted = true
                    break
                end
            end
            n_back += 1
            rt.backtracks += 1
            (n_back > strategy.max_backtracks || lambda/2 < strategy.min_linesearch) && break
            lambda /= 2
        end

        if !accepted
            # repeated rejected steps → restart the inverse approximation
            if restarts_left > 0
                restarts_left -= 1
                rt.restarts += 1
                rt.Jinv .= jinv_scale .* Matrix{TF}(LA.I, M, M)
                # re-evaluate the base point so r_scaled matches a clean state
                trial = trial_fun!(rt.x)
                _kutta_track_best!(rt, trial)
                (trial.finite && trial.inner_status === :converged) || return :fallback
                rt.r_scaled .= trial.residual ./ rt.P_s
                _copy_record!(rt.current, trial)
                continue
            end
            return :fallback
        end

        # accepted step: update iterate and inverse approximation
        rt.work_M2 .= trial.residual ./ rt.P_s .- rt.r_scaled    # Δr̂
        rt.dx .*= lambda                                          # actual step
        rt.x .+= rt.dx
        rt.r_scaled .= trial.residual ./ rt.P_s
        _copy_record!(rt.current, trial)
        dx_inf_scaled = maximum(abs, rt.dx)

        if _kutta_gates_pass(rt, trial, dx_inf_scaled)
            rt.accepted_dc_inf_scaled = dx_inf_scaled
            _copy_record!(rt.accepted, trial)
            return :converged
        end

        # good-Broyden inverse update:
        # Jinv += (dx − Jinv·Δr̂)(dxᵀJinv)/(dxᵀ Jinv Δr̂), guarded
        Jdr = rt.Jinv*rt.work_M2
        denom = LA.dot(rt.dx, Jdr)
        if abs(denom) > sqrt(eps(TF))*LA.norm(rt.dx)*LA.norm(Jdr) && isfinite(denom)
            xtJ = (rt.Jinv')*rt.dx
            num = rt.dx .- Jdr
            for j in 1:M, i in 1:M
                rt.Jinv[i, j] += num[i]*xtJ[j]/denom
            end
        elseif restarts_left > 0
            restarts_left -= 1
            rt.restarts += 1
            rt.Jinv .= jinv_scale .* Matrix{TF}(LA.I, M, M)
        else
            return :fallback
        end
    end
    return :fallback
end

"""
Finite-difference Newton–Krylov fallback (architecture §9.3), started from
the best finite trial. Directional Jacobian products are finite differences
of complete residual transactions.
"""
function _kutta_newton_krylov!(rt::KuttaRuntime{TF}, trial_fun!) where {TF}
    strategy = (rt.closure::PressureContinuityKutta).fallback::NewtonKrylovKutta
    M = length(rt.x)
    rt.method = :newton_krylov
    rt.fallback_entered = true

    rt.best.valid && rt.best.finite || return :failed
    rt.x .= rt.best.c ./ rt.C_s

    # base residual at the fallback start point; a base point that already
    # satisfies the pressure and inner gates is the fallback's initial trial
    # and is accepted with a defined zero step
    trial = trial_fun!(rt.x)
    _kutta_track_best!(rt, trial)
    (trial.finite && trial.inner_status === :converged) || return :failed
    rt.r_scaled .= trial.residual ./ rt.P_s
    _copy_record!(rt.current, trial)
    if _kutta_gates_pass(rt, trial, zero(TF))
        rt.accepted_dc_inf_scaled = zero(TF)
        _copy_record!(rt.accepted, trial)
        return :converged
    end

    for newton in 1:strategy.max_newton
        rt.outer_iterations += 1
        r_base = copy(rt.r_scaled)
        x_base = copy(rt.x)

        # matrix-free FD Jacobian-vector product of full transactions
        function jv!(out, v, alpha, beta)
            vn = LA.norm(v)
            if iszero(vn)
                out .= beta == 0 ? zero(TF) : beta .* out
                return out
            end
            eps_fd = strategy.fd_relative_step*(one(TF) + LA.norm(x_base))/vn
            rt.work_M .= x_base .+ eps_fd .* v
            t = trial_fun!(rt.work_M)
            if !(t.finite && t.inner_status === :converged)
                out .= TF(NaN)
                return out
            end
            jv = (t.residual ./ rt.P_s .- r_base) ./ eps_fd
            if beta == 0
                out .= alpha .* jv
            else
                out .= alpha .* jv .+ beta .* out
            end
            return out
        end
        op = LinearOperators.LinearOperator(TF, M, M, false, false, jv!)
        rhs = .-r_base
        dx, stats = Krylov.gmres(op, rhs; itmax=strategy.max_krylov)
        # an unconverged Krylov solve (including one poisoned by a NaN FD
        # product from an invalid trial) must not be accepted as a direction
        (stats.solved && all(isfinite, dx)) || return :failed

        dx_inf = maximum(abs, dx)
        if dx_inf > strategy.max_scaled_step
            dx .*= strategy.max_scaled_step/dx_inf
            dx_inf = strategy.max_scaled_step
        end

        # accept the current iterate when it already satisfies the pressure
        # and inner gates and the proposed Newton step is within tolerance
        if _kutta_gates_pass(rt, rt.current, dx_inf)
            rt.accepted_dc_inf_scaled = dx_inf
            _copy_record!(rt.accepted, rt.current)
            return :converged
        end

        # same backtracking/acceptance rules as the primary
        r0 = LA.norm(r_base)
        lambda = one(TF)
        accepted = false
        n_back = 0
        local last_trial = trial
        while true
            rt.work_M .= x_base .+ lambda .* dx
            t = trial_fun!(rt.work_M)
            if t.finite && t.inner_status === :converged
                _kutta_track_best!(rt, t)
                rnew = LA.norm(t.residual ./ rt.P_s)
                if rnew <= (1 - TF(1e-4)*lambda)*r0
                    accepted = true
                    last_trial = t
                    break
                end
            end
            n_back += 1
            rt.backtracks += 1
            (n_back > strategy.max_backtracks || lambda/2 < strategy.min_linesearch) && break
            lambda /= 2
        end
        accepted || return :failed

        rt.x .= x_base .+ lambda .* dx
        rt.r_scaled .= last_trial.residual ./ rt.P_s
        _copy_record!(rt.current, last_trial)
        dx_inf_scaled = lambda*maximum(abs, dx)
        if _kutta_gates_pass(rt, last_trial, dx_inf_scaled)
            rt.accepted_dc_inf_scaled = dx_inf_scaled
            _copy_record!(rt.accepted, last_trial)
            return :converged
        end
    end
    return :failed
end

################################################################################
# ROUTE B GEOMETRY
################################################################################

"""
Route B trailing-edge update: anchor wake node row 1 at the bare body TE and
refresh the ephemeral attachment-geometry cache `body.Das .= nodes[:,2,·] − TE`
so the body's attached transition panel coincides exactly with the live wake
panel (TE → first convected downstream row). Overwrites whatever
`rotate_Das!` did between steps; the cache is never persisted as authoritative
wake state. During cold startup (no convected row yet) the cache is zeroed and
the attached wake is suppressed instead of evaluated.
"""
function _update_TE_route_b!(rt::KuttaRuntime)
    wake = rt.wake
    body = rt.body
    isnothing(wake) && return rt
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        shedding = body.shedding[i_surf]
        Das = body.Das[i_surf]
        for i_shed in axes(shedding, 2)
            i_panel = shedding[1, i_shed]
            n = body.cells[shedding[3, i_shed], i_panel]  # nib
            for d in 1:3
                te = body.nodes[d, n]
                nodes[d, 1, i_shed] = te
                Das[d, i_shed] = rt.live_active ? nodes[d, 2, i_shed] - te : zero(te)
            end
        end
        # final node column (nia of the last edge = last Das column)
        i_panel = shedding[1, end]
        n = body.cells[shedding[2, end], i_panel]
        jn = size(nodes, 3)
        jd = size(Das, 2)
        for d in 1:3
            te = body.nodes[d, n]
            nodes[d, 1, jn] = te
            Das[d, jd] = rt.live_active ? nodes[d, 2, jn] - te : zero(te)
        end
    end
    return rt
end

"Reassemble the Route B attachment-aware operator and RHS columns for the
current live-panel geometry (once per physical step / geometry version). At
cold startup the operator is the body-only system (attached wake suppressed)."
function _kutta_update_route_b_operators!(rt::KuttaRuntime{TF}) where {TF}
    body = rt.body
    body.suppress_attached_wake[] = !rt.live_active
    try
        rt.Gb .= zero(TF)
        _G!(rt.Gb, body, body; core_size=body.core_size_panel,
            update_geometry=false)
    finally
        body.suppress_attached_wake[] = false
    end
    rt.Glu = lu!(rt.Gb)
    if rt.live_active
        _assemble_W!(rt.W, body, rt.edges)
    else
        rt.W .= zero(TF)
    end
    rt.geometry_version += 1
    return rt
end

"Route B end-of-step topology advancement bookkeeping, called by `simulate!`
after the unchanged `propagate!`/`shed_wake!` pair: the accepted live block has
been shifted into old-wake storage with its strength retained, and the fresh
row-1 deposit is the reserved next live slot (a warm-start prior, excluded
from old-wake source views and overwritten at the next commit). Not a second
commit; the closure is not called."
function _kutta_advance_topology!(rt::KuttaRuntime, i_step::Int)
    rt.route === :B || return rt
    wake = rt.wake
    isnothing(wake) && return rt
    wake.live_rows[] = 1
    wake.live_step_id[] = i_step + 1
    rt.live_active = true
    return rt
end

################################################################################
# STEP ORCHESTRATOR
################################################################################

function _reset_kutta_counters!(rt::KuttaRuntime{TF}) where {TF}
    rt.body_solves = 0
    rt.backtracks = 0
    rt.restarts = 0
    rt.outer_iterations = 0
    rt.fallback_entered = false
    rt.method = :none
    rt.accepted_dc_inf_scaled = TF(NaN)
    rt.trial.valid = false
    rt.best.valid = false
    rt.accepted.valid = false
    rt.current.valid = false
    return rt
end

"""
    _kutta_step!(rt, systems_tuple, wakes_tuple, frames, uinf; ...)

Non-default replacement for `_steady_aerodynamics!`: runs the identical
pre-solve stages (reset, freestream, kinematics, trailing-edge update, frozen
wake influence), freezes the step snapshot, executes the configured closure
(jump directly; pressure through the trial transaction and nonlinear
drivers), and commits the accepted state exactly once — including the single
body→wake-probes influence pass and the exterior half-jump, in the legacy
relative order.
"""
function _kutta_step!(rt::KuttaRuntime{TF}, systems_tuple::Tuple,
        wakes_tuple::Tuple, frames, uinf;
        backend_wake=nothing, backend_solve, backend_system,
        needs_induced_vorticity::Bool=false,
        grad_mu_options=(;),
        i_step::Int=0,
        wakerow_no_hessian_to_particles::Bool=false,
        body_hessian_to_particles::Bool=false,
        body_gradient_core_size::Float64=NaN,
        body_on_wake::Bool=true,
        panel_wake_on_particles::Bool=true,
        particle_hessian_self::Bool=true,
        particle_body_overlap_policy::ParticleBodyOverlapPolicy=ParticleBodyOverlapPolicy()) where {TF}

    body = rt.body
    closure = rt.closure
    normalized_gm = _normalize_grad_mu_options(grad_mu_options;
        default_basis=:quad)

    # ----- legacy pre-solve stages (identical statements/order) -----
    wake_probes, targets, wake_sources = _sa_collect(systems_tuple, wakes_tuple)
    _sa_reset_freestream_kinematic!(systems_tuple, wakes_tuple, frames, uinf)

    if rt.route === :A
        for (sys, w) in zip(systems_tuple, wakes_tuple)
            !isnothing(w) && update_TE!(w, sys)
        end
    else
        rt.live_active = !isnothing(rt.wake) && rt.wake.live_rows[] > 0
        _update_TE_route_b!(rt)
    end

    check_particle_body_overlap!(particle_body_overlap_policy,
        systems_tuple, wakes_tuple, i_step)

    _sa_wake_influence!(targets, wake_sources, backend_wake;
        needs_induced_vorticity, wakerow_no_hessian_to_particles,
        panel_wake_on_particles, particle_hessian_self)

    # ----- per-step attachment operators and frozen data -----
    rt.route === :B && _kutta_update_route_b_operators!(rt)
    _kutta_update_edge_lengths!(rt)
    _kutta_capture_snapshot!(rt)
    _reset_kutta_counters!(rt)

    commit_kwargs = (; i_step, targets, systems_tuple, backend_system,
        normalized_gm, needs_induced_vorticity, body_on_wake,
        body_hessian_to_particles, body_gradient_core_size)
    trial_kwargs = (; backend_solve, backend_system,
        grad_mu_options=normalized_gm)

    # ----- Route B cold startup: deterministic jump solve, no live panel -----
    if rt.route === :B && !rt.live_active
        rt.work_M .= zero(TF)
        trial = _kutta_trial!(rt, rt.work_M; trial_kwargs...,
            reconstruct=false, evaluate=false, suppress_attached=true)
        trial.inner_status === :converged || throw(KuttaConvergenceError(
            _kutta_diagnostics_record(rt, i_step, :failed, :error),
            :inner_solve, "Route B cold-startup jump solve produced a "*
            "non-finite body solution."))
        _kutta_commit!(rt, trial; startup=true, commit_kwargs...)
        if closure isa PressureContinuityKutta && closure.store_diagnostics
            push!(closure.diagnostics,
                _kutta_diagnostics_record(rt, i_step, :startup_jump, :startup_jump))
        end
        return nothing
    end

    # ----- jump closure: single coupled solve at c = 0 -----
    if closure isa JumpKutta
        rt.work_M .= zero(TF)
        trial = _kutta_trial!(rt, rt.work_M; trial_kwargs...,
            reconstruct=false, evaluate=false)
        trial.inner_status === :converged || throw(KuttaConvergenceError(
            _kutta_diagnostics_record(rt, i_step, :failed, :error),
            :inner_solve, "jump-closure body solve produced a non-finite "*
            "solution."))
        _kutta_commit!(rt, trial; commit_kwargs...)
        return nothing
    end

    # ----- pressure-continuity closure -----
    pclosure = closure::PressureContinuityKutta

    # initial clean trial (warm-started correction) completes the exterior
    # velocity from which the automatic scales are frozen for this step
    trial = _kutta_trial!(rt, rt.c; trial_kwargs..., evaluate=false)
    if !(trial.finite && trial.inner_status === :converged)
        _kutta_fail!(rt, i_step, :inner_solve,
            "the initial clean trial produced a non-finite body solution";
            trial_kwargs, commit_kwargs)
        return nothing
    end
    try
        _kutta_freeze_scales!(rt)
    catch
        # a degenerate automatic scale aborts the step; the initial clean
        # trial has already mutated the body, so restore the pre-step state
        # before the ArgumentError escapes (§7.4 determinism)
        _kutta_restore_snapshot!(rt)
        rethrow()
    end

    trial_fun! = x -> begin
        rt.trial.c .= x .* rt.C_s
        _kutta_trial!(rt, rt.trial.c; trial_kwargs...)
    end

    status = _kutta_broyden!(rt, trial_fun!)
    if status === :fallback
        status = _kutta_newton_krylov!(rt, trial_fun!)
        status === :converged || (status = :failed)
    end

    if status === :converged
        _kutta_commit!(rt, rt.accepted; commit_kwargs...)
        if pclosure.store_diagnostics
            push!(pclosure.diagnostics,
                _kutta_diagnostics_record(rt, i_step, :converged, :converged))
        end
        return nothing
    end

    _kutta_fail!(rt, i_step, :nonconvergence,
        "the pressure-continuity solve did not satisfy the pressure, "*
        "correction-step, and inner-solve gates";
        trial_kwargs, commit_kwargs)
    return nothing
end

################################################################################
# COMMIT AND FAILURE
################################################################################

"""
One-time commit of an accepted trial record (architecture §7.3): restore the
frozen snapshot, install the accepted strengths and correction, write the
Route B live strength once, run the single accepted-field reconstruction
(body → bodies and wake probes) with the correction in physical mode, and
call `commit_pressure!` once. Never re-runs the trial.
"""
function _kutta_commit!(rt::KuttaRuntime{TF}, record::KuttaTrialRecord;
        startup::Bool=false, i_step::Int=0, targets, systems_tuple,
        backend_system, normalized_gm, needs_induced_vorticity::Bool,
        body_on_wake::Bool, body_hessian_to_particles::Bool,
        body_gradient_core_size::Float64) where {TF}
    try
        _kutta_commit_inner!(rt, record; startup, targets, systems_tuple,
            backend_system, normalized_gm, needs_induced_vorticity,
            body_on_wake, body_hessian_to_particles,
            body_gradient_core_size)
    catch err
        # §7.4: commit cannot complete — restore the complete pre-step
        # snapshot (which also rolls back the provider) and fail
        # deterministically. The half-installed accepted state is never left
        # on the body or wake.
        _kutta_restore_snapshot!(rt)
        throw(KuttaConvergenceError(
            _kutta_diagnostics_record(rt, i_step, :failed, :error), :commit,
            "commit could not complete: "*sprint(showerror, err)))
    end
    return nothing
end

function _kutta_commit_inner!(rt::KuttaRuntime{TF}, record::KuttaTrialRecord;
        startup::Bool, targets, systems_tuple, backend_system,
        normalized_gm, needs_induced_vorticity::Bool,
        body_on_wake::Bool, body_hessian_to_particles::Bool,
        body_gradient_core_size::Float64) where {TF}

    body = rt.body
    _kutta_restore_snapshot!(rt)

    # accepted strengths: σ from the frozen velocity, μ from the record
    set_strengths!(body)
    view(body.strength, :, 2) .= record.mu

    # physical-mode correction γ = Cμ − c for every downstream consumer
    if any(!iszero, record.c)
        set_wake_correction!(body, record.c)
    else
        clear_wake_correction!(body)
    end
    rt.c .= record.c

    # Route B: the single live-strength write of the accepted circulation
    if rt.route === :B && !startup && !isnothing(rt.wake)
        wake = rt.wake
        k = 0
        for (i_surf, shed) in enumerate(body.shedding)
            for i in axes(shed, 2)
                k += 1
                wake.strength[i_surf][1, 1, i] = _get_wakestrength_Gamma(body, i, i_surf)
            end
        end
    end

    # single accepted-field exposure: body → (bodies, wake probes) with the
    # correction inactive (operator mode, any backend), completed by the
    # direct affine add-on on the body and the wake probes, then the exterior
    # half-jump — the legacy post-solve order
    body.velocity .= rt.snapshot.velocity_frozen
    needs_induced_vorticity && _add_bound_surface_vorticity!(systems_tuple;
        grad_mu_options=normalized_gm)
    _set_core_sizes!(systems_tuple, :core_size_targets)
    body.suppress_attached_wake[] = startup
    was_active = _operator_mode_begin!(body)
    try
        _sa_body_influence!(targets, systems_tuple, backend_system;
            needs_induced_vorticity, body_on_wake, body_hessian_to_particles,
            body_gradient_core_size)
    finally
        body.suppress_attached_wake[] = false
        _operator_mode_end!(body, was_active)
    end
    if !startup && any(!iszero, record.c)
        _add_affine_attached_velocity!(body.velocity, body.controlpoints,
            body, record.c; core_size=body.core_size_panel)
        if !isnothing(rt.wake) && body_on_wake
            _add_affine_attached_velocity!(rt.wake, body, record.c;
                core_size=body.core_size_targets)
        end
    end
    _sa_half_jump!(systems_tuple, normalized_gm)

    provider = _kutta_provider(rt)
    trial_view = KuttaTrialView{TF, typeof(body)}(body, rt.edges, body.velocity)
    commit_pressure!(provider, trial_view)

    _copy_record!(rt.accepted, record)
    return nothing
end

"""
Failure path (architecture §7.4). `on_failure=:error`: restore the complete
pre-step snapshot and throw [`KuttaConvergenceError`](@ref). `:jump`: restore,
run and commit a fresh c = 0 jump solve (which must itself succeed), warn,
and record the fallback disposition. The last nonlinear trial is never left
installed and never relabeled as the jump solution.
"""
function _kutta_fail!(rt::KuttaRuntime{TF}, i_step::Int, cause::Symbol,
        msg::String; trial_kwargs, commit_kwargs) where {TF}
    closure = rt.closure::PressureContinuityKutta

    if closure.on_failure === :error
        _kutta_restore_snapshot!(rt)
        throw(KuttaConvergenceError(
            _kutta_diagnostics_record(rt, i_step, :failed, :error), cause, msg))
    end

    # explicit physical fallback: fresh coupled jump solve at c = 0
    rt.work_M .= zero(TF)
    trial = _kutta_trial!(rt, rt.work_M; trial_kwargs...,
        reconstruct=false, evaluate=false)
    if trial.inner_status !== :converged
        _kutta_restore_snapshot!(rt)
        throw(KuttaConvergenceError(
            _kutta_diagnostics_record(rt, i_step, :failed, :error),
            :jump_fallback, "the explicit jump fallback solve failed after: "*msg))
    end
    _kutta_commit!(rt, trial; commit_kwargs...)
    @warn "PressureContinuityKutta did not converge at step $(i_step); "*
        "committed a fresh jump (c = 0) solution instead (on_failure=:jump). "*
        "Cause: $(cause)." maxlog=10
    if closure.store_diagnostics
        push!(closure.diagnostics,
            _kutta_diagnostics_record(rt, i_step, :jump_fallback, :jump_fallback))
    end
    return nothing
end

################################################################################
# METADATA AND WARM START
################################################################################

_kutta_scale_metadata(scale::Float64) = isnan(scale) ? "auto" : scale

"Manifest metadata for a non-default Kutta configuration (architecture
§10.2): concrete names, provider parameters, strategy safeguards, tolerances,
and policies. Describes committed configuration only."
function _kutta_manifest_dict(rt::KuttaRuntime)
    d = Dict{String, Any}(
        "wake_attachment" => rt.route === :B ? "TEAnchoredAttachment" :
            "RigidTransitionAttachment",
        "kutta_closure" => String(nameof(typeof(rt.closure))),
    )
    closure = rt.closure
    if closure isa PressureContinuityKutta
        provider = closure.provider
        d["provider"] = String(nameof(typeof(provider)))
        provider isa SteadyBernoulliProvider && (d["provider_rho"] = Float64(provider.rho))
        d["pressure_scale"] = _kutta_scale_metadata(closure.pressure_scale)
        d["correction_scale"] = _kutta_scale_metadata(closure.correction_scale)
        d["pressure_tolerance"] = closure.pressure_tolerance
        d["correction_tolerance"] = closure.correction_tolerance
        d["on_failure"] = String(closure.on_failure)
        d["store_diagnostics"] = closure.store_diagnostics
        p = closure.primary
        if p isa BroydenKutta
            d["primary"] = Dict{String, Any}("name" => "BroydenKutta",
                "max_iterations" => p.max_iterations,
                "max_backtracks" => p.max_backtracks,
                "max_restarts" => p.max_restarts,
                "max_scaled_step" => p.max_scaled_step,
                "min_linesearch" => p.min_linesearch)
        end
        f = closure.fallback
        if f isa NewtonKrylovKutta
            d["fallback"] = Dict{String, Any}("name" => "NewtonKrylovKutta",
                "max_newton" => f.max_newton,
                "max_krylov" => f.max_krylov,
                "fd_relative_step" => f.fd_relative_step,
                "max_scaled_step" => f.max_scaled_step,
                "max_backtracks" => f.max_backtracks,
                "min_linesearch" => f.min_linesearch)
        end
    end
    return d
end

"Per-step accepted Kutta state for the metadata step record (architecture
§10.2-10.3): the accepted correction, accepted circulation, live-block
identifiers, and disposition — everything a Route B or A/pressure warm start
must restore. Called after commit, so the body carries the accepted
correction in physical mode."
function _kutta_step_dict(rt::KuttaRuntime)
    body = rt.body
    gamma = Float64[]
    for (i_surf, shed) in enumerate(body.shedding)
        for i in axes(shed, 2)
            push!(gamma, Float64(_get_wakestrength_Gamma(body, i, i_surf)))
        end
    end
    d = Dict{String, Any}(
        "c" => Float64.(rt.c),
        "gamma_accepted" => gamma,
        "route" => String(rt.route),
    )
    if !isnothing(rt.wake)
        d["live_rows"] = rt.wake.live_rows[]
        d["live_step_id"] = rt.wake.live_step_id[]
    end
    if rt.closure isa PressureContinuityKutta
        closure = rt.closure::PressureContinuityKutta
        if !isempty(closure.diagnostics)
            d["status"] = String(closure.diagnostics[end].status)
        end
    end
    return d
end

"""
Warm-start compatibility validation and accepted-state restoration
(architecture §10.3), called by `simulate_warmstart!` after body/wake state
is loaded and before the end-of-step replay. Validates that the saved
configuration matches the requested one, reinstalls the committed correction
(so the replayed `shed_wake!` deposits γ = Cμ − c), and restores Route B
live-block metadata. Returns `nothing` on the legacy path.
"""
function _kutta_warmstart_restore!(systems_tuple::Tuple, wakes_tuple::Tuple,
        metadata, restart_step::Int,
        wake_attachment::AbstractWakeAttachment,
        kutta_closure::AbstractKuttaClosure)
    saved = isnothing(metadata) ? nothing : get(metadata, "kutta", nothing)
    if _is_legacy_kutta(wake_attachment, kutta_closure)
        isnothing(saved) || throw(ArgumentError(
            "the saved run used a non-default wake_attachment/kutta_closure "*
            "($(saved["wake_attachment"])/$(saved["kutta_closure"])); pass the "*
            "same configuration to simulate_warmstart!."))
        return nothing
    end
    isnothing(saved) && throw(ArgumentError(
        "the saved run has no [kutta] metadata but a non-default "*
        "wake_attachment/kutta_closure was requested; a legacy run cannot be "*
        "resumed under a different Kutta configuration."))
    want_attachment = wake_attachment isa TEAnchoredAttachment ?
        "TEAnchoredAttachment" : "RigidTransitionAttachment"
    want_closure = String(nameof(typeof(kutta_closure)))
    (saved["wake_attachment"] == want_attachment &&
        saved["kutta_closure"] == want_closure) || throw(ArgumentError(
        "saved Kutta configuration ($(saved["wake_attachment"])/"*
        "$(saved["kutta_closure"])) does not match the requested "*
        "($(want_attachment)/$(want_closure))."))

    # locate the step record for the restart step
    steps = get(metadata, "step", Any[])
    idx = findfirst(s -> Int(s["i_step"]) == restart_step && haskey(s, "kutta"),
        steps)
    isnothing(idx) && throw(ArgumentError(
        "no [step.kutta] record found for restart_step=$(restart_step); the "*
        "saved run is missing the accepted Kutta state needed to resume "*
        "(incomplete new-format state is rejected rather than reseeded)."))
    krec = steps[idx]["kutta"]

    body = systems_tuple[1]
    c = Float64.(krec["c"])
    if any(!iszero, c)
        set_wake_correction!(body, c)
    else
        clear_wake_correction!(body)
    end

    if wake_attachment isa TEAnchoredAttachment
        (haskey(krec, "live_rows") && haskey(krec, "live_step_id")) ||
            throw(ArgumentError("saved Route B step record is missing "*
                "live_rows/live_step_id; incomplete state is rejected."))
        wake = wakes_tuple[1]
        wake.live_rows[] = Int(krec["live_rows"])
        wake.live_step_id[] = Int(krec["live_step_id"])
    end
    return nothing
end

"Assemble the diagnostics record for the current runtime counters and the most
relevant trial (accepted if valid, else best, else last)."
function _kutta_diagnostics_record(rt::KuttaRuntime{TF}, i_step::Int,
        status::Symbol, disposition::Symbol) where {TF}
    rec = rt.accepted.valid ? rt.accepted : (rt.best.valid ? rt.best : rt.trial)
    startup = status === :startup_jump ? :startup_jump : :none
    return KuttaDiagnostics{Float64}(status, rt.method, rt.route,
        rt.outer_iterations, rt.body_solves, rt.backtracks, rt.restarts,
        Float64(rec.valid ? rec.r_inf : NaN),
        Float64(rec.valid && !isnan(rt.P_s) ? rec.r_inf/rt.P_s : NaN),
        Float64(rec.valid ? rec.r_l2_edge : NaN),
        Float64(rec.valid && !isnan(rt.P_s) ? rec.r_l2_edge/rt.P_s : NaN),
        Float64(!isnan(rt.accepted_dc_inf_scaled) && !isnan(rt.C_s) ?
            rt.accepted_dc_inf_scaled*rt.C_s : NaN),
        Float64(rt.accepted_dc_inf_scaled),
        Float64(rt.P_s), Float64(rt.C_s), Float64(rt.U_s), Float64(rt.L_s),
        rec.valid ? rec.inner_status : :not_run, 1, NaN,
        rt.fallback_entered, startup, disposition, i_step,
        Float64.(rec.valid ? rec.c : rt.c))
end
