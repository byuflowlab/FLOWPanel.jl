#--- monitor dependency traits ---#

abstract type AbstractMonitor end

# Each monitor declares which body fields it populates (`monitor_provides`)
# and which it consumes (`monitor_requires`). `audit_monitors` checks at the
# start of `simulate!` that every requirement is met by an earlier monitor
# in the tuple, surfacing ordering bugs before the time loop starts.
#
# Field symbols correspond to per-body arrays touched by monitors:
#   :P  → body.P (pressure)
#   :F  → body.F (distributed surface force)
#
# Add a new monitor's contract by defining methods on these two traits.

monitor_provides(::Any) = ()
monitor_requires(::Any) = ()

# Monitors that need ∇u at body control points return true here so
# simulate! can flip body.needs_velocity_gradient[] before the time loop.
monitor_requires_body_hessian(::Any) = false

"""
    audit_monitors(monitors)

Walk `monitors` left-to-right and verify that every monitor's
`monitor_requires` symbols are produced by some earlier monitor's
`monitor_provides`. Throws `ArgumentError` on the first violation; returns
`monitors` unchanged on success. Called automatically by `simulate!`.
"""
function audit_monitors(monitors)
    provided = Set{Symbol}()
    for (i, m) in enumerate(monitors)
        for r in monitor_requires(m)
            if !(r in provided)
                throw(ArgumentError(
                    "Monitor at position $(i) ($(nameof(typeof(m)))) requires " *
                    "body.$(r) to be populated, but no earlier monitor in the " *
                    "tuple provides it. Insert a monitor whose " *
                    "`monitor_provides` includes :$(r) (e.g. PressureBernoulli " *
                    "for :P, ForceMonitor for :F) before position $(i)."))
            end
        end
        for p in monitor_provides(m)
            push!(provided, p)
        end
    end
    return monitors
end

"""
    PressureBernoulli(rho; unsteady=false, correct_kuttacondition=true, clip=nothing,
                      backend=FastMultipoleBackend(expansion_order=10,
                                                   multipole_acceptance=0.4,
                                                   leaf_size=100))

Monitor that populates `body.P` on every body in the simulation by evaluating
the Bernoulli equation each step. With `unsteady=false` (default) the steady
form ``P = \\tfrac{1}{2} \\rho (U_\\infty^2 - U^2)`` is used; with
`unsteady=true` the term ``-\\rho \\, \\partial \\phi / \\partial t`` is added.
The monitor evaluates scalar potential at the body control points, finite
differences it between calls, and subtracts the control-point convective term.
Only source systems with a scalar potential are included in that potential
evaluation.

`correct_kuttacondition` and `clip` are forwarded to `calcfield_P!`.

Pressure-dependent monitors (e.g. `ForceMonitor`) must appear *after* this
monitor in the `monitors` tuple passed to `simulate!`.
"""
mutable struct PressureBernoulli{TF, TC, TB} <: AbstractMonitor
    rho::TF
    unsteady::Bool
    correct_kuttacondition::Bool
    clip::TC
    backend::TB
    phi_dot::Vector{Vector{Float64}}
    potential_history::Vector{Vector{Float64}}
    probes::Vector{FastMultipole.ProbeSystem{Float64}}
end

monitor_provides(::PressureBernoulli) = (:P,)

function PressureBernoulli(rho::Real; unsteady::Bool=false,
                          correct_kuttacondition::Bool=true,
                          clip=nothing,
                          backend=FastMultipoleBackend(expansion_order=10,
                                                       multipole_acceptance=0.4,
                                                       leaf_size=100))
    return PressureBernoulli{typeof(rho), typeof(clip), typeof(backend)}(
        rho, unsteady, correct_kuttacondition, clip, backend,
        Vector{Float64}[], Vector{Float64}[], FastMultipole.ProbeSystem{Float64}[])
end

function (m::PressureBernoulli)(systems, wakes,
                               frames::AbstractVector{<:ReferenceFrame},
                               uinf, i_step::Int, dt::Real)
    systems_tuple = _systems_tuple(systems)
    Uinf_mag = norm(uinf)
    scalar_sources = m.unsteady ? _bernoulli_scalar_sources(systems_tuple, wakes) : ()
    m.unsteady && _pressure_bernoulli_ensure_storage!(m, systems_tuple)
    for (i_body, body) in enumerate(systems_tuple)
        fill!(body.P, zero(eltype(body.P)))
        phi_dot = m.unsteady ?
            _pressure_bernoulli_phi_dot!(m, body, i_body, scalar_sources, uinf, dt) :
            nothing
        calcfield_P!(body.P, body, body.velocity, Uinf_mag, m.rho, phi_dot;
                     correct_kuttacondition=m.correct_kuttacondition,
                     clip=m.clip)
    end
end

function _pressure_bernoulli_ensure_storage!(m::PressureBernoulli, systems_tuple::Tuple)
    nbodies = length(systems_tuple)
    if length(m.phi_dot) != nbodies ||
       length(m.potential_history) != nbodies ||
       length(m.probes) != nbodies
        m.phi_dot = Vector{Float64}[]
        m.potential_history = Vector{Float64}[]
        m.probes = FastMultipole.ProbeSystem{Float64}[]
        sizehint!(m.phi_dot, nbodies)
        sizehint!(m.potential_history, nbodies)
        sizehint!(m.probes, nbodies)
        for body in systems_tuple
            push!(m.phi_dot, zeros(Float64, body.ncells))
            push!(m.potential_history, zeros(Float64, body.ncells))
            push!(m.probes, FastMultipole.ProbeSystem(body.ncells, Float64))
        end
        return nothing
    end

    for (i, body) in enumerate(systems_tuple)
        if length(m.phi_dot[i]) != body.ncells ||
           length(m.potential_history[i]) != body.ncells ||
           length(m.probes[i].scalar_potential) != body.ncells
            m.phi_dot[i] = zeros(Float64, body.ncells)
            m.potential_history[i] = zeros(Float64, body.ncells)
            m.probes[i] = FastMultipole.ProbeSystem(body.ncells, Float64)
        end
    end

    return nothing
end

function _bernoulli_scalar_sources(systems_tuple::Tuple, wakes)
    wake_sources = _collect_wake_scalar_sources(wakes)
    return _filter_scalar_potential_sources((systems_tuple..., wake_sources...))
end

function _filter_scalar_potential_sources(sources::Tuple)
    result = ()
    for source in sources
        if !FastMultipole.has_vector_potential(source)
            result = (result..., source)
        end
    end
    return result
end

function _pressure_bernoulli_phi_dot!(m::PressureBernoulli, body::AbstractBody,
                                      i_body::Int, scalar_sources::Tuple,
                                      uinf, dt::Real)
    dt > 0 || throw(ArgumentError("PressureBernoulli(unsteady=true) requires a positive runtime dt; got $(dt)."))

    probes = m.probes[i_body]
    history = m.potential_history[i_body]
    phi_dot = m.phi_dot[i_body]
    inv_dt = inv(Float64(dt))

    @inbounds for p in 1:body.ncells
        probes.position[p] = SVector{3, Float64}(
            body.controlpoints[1, p],
            body.controlpoints[2, p],
            body.controlpoints[3, p],
        )
        probes.scalar_potential[p] = 0.0
    end

    if length(scalar_sources) > 0
        influence!((probes,), scalar_sources, m.backend;
            precalc=false, scalar_potential=true, gradient=false,
            hessian=(false,))
    end

    @inbounds for p in 1:body.ncells
        phi = probes.scalar_potential[p] +
              uinf[1] * body.controlpoints[1, p] +
              uinf[2] * body.controlpoints[2, p] +
              uinf[3] * body.controlpoints[3, p]
        Dphi_Dt = (history[p] + phi) * inv_dt
        vk1 = body.velocity_kinematic[1, p]
        vk2 = body.velocity_kinematic[2, p]
        vk3 = body.velocity_kinematic[3, p]
        grad_phi_1 = body.velocity[1, p] + vk1
        grad_phi_2 = body.velocity[2, p] + vk2
        grad_phi_3 = body.velocity[3, p] + vk3
        phi_dot[p] = Dphi_Dt - (vk1 * grad_phi_1 + vk2 * grad_phi_2 + vk3 * grad_phi_3)
        history[p] = -phi
    end

    return phi_dot
end

abstract type AbstractPressurePreconditioner end

struct JacobiPressurePreconditioner{TS} <: AbstractPressurePreconditioner
    states::TS
end

JacobiPressurePreconditioner() =
    JacobiPressurePreconditioner(PressureJacobiState{Vector{Float64}}[])

struct NoPressurePreconditioner <: AbstractPressurePreconditioner end
struct IncompleteCholeskyPressurePreconditioner <: AbstractPressurePreconditioner end
struct AMGPressurePreconditioner <: AbstractPressurePreconditioner end

struct PressureJacobiState{TV}
    invdiag::TV
end

function LA.ldiv!(y::AbstractVector, M::PressureJacobiState, x::AbstractVector)
    @inbounds for i in eachindex(x)
        y[i] = M.invdiag[i] * x[i]
    end
    return y
end

"""
    PressureLaplace(bodies, rho; atol=1e-8, rtol=1e-8, itmax=200,
                    preconditioner=JacobiPressurePreconditioner(),
                    reference_panel=1, reference_pressure=0.0,
                    cache=true, verbose=false)

Monitor that populates `body.P` by solving a sparse panel-centered surface
pressure Poisson equation. The monitor uses `velocity_dot` as a rolling
time-difference buffer of the tangent-projected **inertial** fluid velocity
`u_inertial = tangent(body.velocity) + body.velocity_kinematic`: between calls it
stores `-u_inertial_old`, then during the call it becomes
`(u_inertial_new - u_inertial_old) / dt`, which is the panel-following rate
`d/dt[u(x_p(t), t)]`. The material derivative used in the RHS is then
`Du/Dt = d/dt[u(x_p,t)] + (u_rel · ∇) u_inertial`, where
`u_rel = body.velocity` is the body-relative slip velocity (on an impermeable
surface `u_rel · n = 0` in the continuous limit).
`velocity_dot` is initialized to zero, so on the first call the unsteady term
is `u_inertial / dt` rather than a true finite difference; if this warm-up
transient matters, pre-populate
`monitor.velocity_dot[i]` with `-(tangent(body.velocity) .+ body.velocity_kinematic)`
before the first call.

The convective term `(u_rel · ∇) u_inertial` uses the analytic spatial
Jacobian `∇u_inertial = body.velocity_gradient + [Ω]_×`, where
`body.velocity_gradient` is the multipole Hessian populated by FastMultipole
during the per-step `influence!` calls (the monitor signals this need via
`monitor_requires_body_hessian(::PressureLaplace) = true`, which simulate!
turns into `body.needs_velocity_gradient[] = true` before the time loop), and
`[Ω]_×` is the skew-symmetric tensor of the body's net angular velocity
accumulated in `body.angular_velocity` by `kinematic_velocity!`. The previous
least-squares surface-gradient estimator has been removed; the new ∇u is
mesh-independent and accurate even on high-AR triangulations.

`bodies` must match the bodies later passed to `simulate!`, in the same order,
so the monitor can preallocate one velocity-difference matrix, sparse operator,
preconditioner, Krylov workspace, and scratch arrays per body. Shared panel
edges are stored as one `4 × nedges` integer matrix per body in `edges`; rows
are `(node_a, node_b, panel_i, panel_j)`, where `(node_a, node_b)` is the
undirected shared edge and `panel_i`, `panel_j` are the adjacent panels.
The timestep `dt` is supplied by `simulate!` at monitor-call time and is not
stored in the monitor.

# Fields
- `rho::Float64` — fluid density used to scale the pressure RHS and output
- `atol::Float64` — absolute residual tolerance for the CG solve
- `rtol::Float64` — relative residual tolerance for the CG solve
- `itmax::Int` — maximum CG iterations per call
- `preconditioner::TP` — preconditioner applied during the CG solve (e.g. `JacobiPressurePreconditioner`)
- `reference_panel::Int` — index of the panel whose pressure is pinned to `reference_pressure` (gauge condition)
- `reference_pressure::Float64` — pressure value assigned to the reference panel
- `cache::Bool` — when `true`, the sparse operator and preconditioner are reused if the geometry signature is unchanged
- `verbose::Bool` — when `true`, prints per-step rebuild status and panel count to stdout
- `velocity_dot::Vector{Matrix{Float64}}` — rolling buffer; initialized to zero, then between calls holds `-u_inertial_old` with `u_inertial = tangent(body.velocity) + body.velocity_kinematic` (3 × ncells per body); during a call becomes `(u_inertial_new - u_inertial_old) / dt`
- `u_inertial::Vector{Matrix{Float64}}` — scratch buffer for the inertial fluid velocity `tangent(body.velocity) + body.velocity_kinematic` (3 × ncells per body)
- `L::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}}` — sparse FV surface Laplacian per body, gauge-fixed at `reference_panel`
- `b::Vector{Vector{Float64}}` — RHS vector per body (length ncells); rebuilt each call from the tangential material acceleration
- `p::Vector{Vector{Float64}}` — solution vector per body (length ncells); written to `body.P` after each solve
- `acceleration::Vector{Matrix{Float64}}` — material acceleration `Du/Dt` scratch buffer (3 × ncells per body)
- `tangential::Vector{Matrix{Float64}}` — tangential projection of `acceleration` (3 × ncells per body)
- `edges::Vector{Matrix{Int}}` — shared interior edges per body (4 × nedges); rows are `(node_a, node_b, panel_i, panel_j)`
- `workspace::Vector{Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}}` — preallocated Krylov CG workspace per body
- `geometry_signature::Vector{UInt64}` — hash of panel count, connectivity, edge lengths, and control-point distances per body; triggers a rebuild when changed
"""
mutable struct PressureLaplace{TP} <: AbstractMonitor
    rho::Float64
    atol::Float64
    rtol::Float64
    itmax::Int
    preconditioner::TP
    reference_panel::Int
    reference_pressure::Float64
    cache::Bool
    verbose::Bool
    velocity_dot::Vector{Matrix{Float64}}
    L::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}}
    b::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}
    acceleration::Vector{Matrix{Float64}}
    tangential::Vector{Matrix{Float64}}
    edges::Vector{Matrix{Int}}
    workspace::Vector{Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}}
    geometry_signature::Vector{UInt64}
    u_inertial::Vector{Matrix{Float64}}
end

monitor_provides(::PressureLaplace) = (:P,)
monitor_requires_body_hessian(::PressureLaplace) = true

function PressureLaplace(bodies, rho::Real;
                         atol::Real=1e-8,
                         rtol::Real=1e-8,
                         itmax::Integer=200,
                         preconditioner::AbstractPressurePreconditioner=JacobiPressurePreconditioner(),
                         reference_panel::Integer=1,
                         reference_pressure::Real=0.0,
                         cache::Bool=true,
                         verbose::Bool=false)
    reference_panel >= 1 || throw(ArgumentError("reference_panel must be at least 1; got $(reference_panel)."))
    systems_tuple = _systems_tuple(bodies)
    isempty(systems_tuple) && throw(ArgumentError("PressureLaplace requires at least one body."))
    nbodies = length(systems_tuple)
    velocity_dot = Matrix{Float64}[]
    Ls = SparseArrays.SparseMatrixCSC{Float64, Int}[]
    bs = Vector{Float64}[]
    ps = Vector{Float64}[]
    acceleration = Matrix{Float64}[]
    tangential = Matrix{Float64}[]
    edges = Matrix{Int}[]
    workspace = Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}[]
    signatures = UInt64[]
    u_inertial = Matrix{Float64}[]
    sizehint!(velocity_dot, nbodies)
    sizehint!(Ls, nbodies)
    sizehint!(bs, nbodies)
    sizehint!(ps, nbodies)
    sizehint!(acceleration, nbodies)
    sizehint!(tangential, nbodies)
    sizehint!(edges, nbodies)
    sizehint!(workspace, nbodies)
    sizehint!(signatures, nbodies)
    sizehint!(u_inertial, nbodies)

    for body in systems_tuple
        body.ncells > 0 || throw(ArgumentError("PressureLaplace requires bodies with at least one panel."))
        reference_panel <= body.ncells || throw(ArgumentError(
            "reference_panel=$(reference_panel) exceeds body.ncells=$(body.ncells)."))
        calc_normals!(body)
        calc_controlpoints!(body; off=abs(body.CPoffset))
        body_edges = _pressure_panel_edges(body)
        L = _assemble_pressure_laplacian(body, Int(reference_panel), body_edges)
        b = zeros(Float64, body.ncells)
        push!(velocity_dot, zeros(Float64, size(body.velocity)))
        push!(Ls, L)
        push!(bs, b)
        push!(ps, zeros(Float64, body.ncells))
        push!(acceleration, zeros(Float64, 3, body.ncells))
        push!(tangential, zeros(Float64, 3, body.ncells))
        push!(edges, body_edges)
        push!(workspace, Krylov.krylov_workspace(Val(:cg), L, b))
        push!(signatures, _pressure_geometry_signature(body, Int(reference_panel), body_edges))
        push!(u_inertial, zeros(Float64, size(body.velocity)))
    end
    preconditioner = build_pressure_preconditioner(preconditioner, Ls)

    return PressureLaplace{typeof(preconditioner)}(
        Float64(rho), Float64(atol), Float64(rtol), Int(itmax),
        preconditioner, Int(reference_panel), Float64(reference_pressure),
        cache, verbose, velocity_dot, Ls, bs, ps,
        acceleration, tangential,
        edges, workspace, signatures, u_inertial)
end

function PressureLaplace(rho::Real, dt::Real; optargs...)
    throw(ArgumentError(
        "PressureLaplace requires the simulation body or tuple of bodies for preallocation, " *
        "and receives dt from simulate! at runtime. Use PressureLaplace(bodies, rho; ...)."))
end

function PressureLaplace(rho::Real; optargs...)
    throw(ArgumentError(
        "PressureLaplace requires the simulation body or tuple of bodies for preallocation. " *
        "Use PressureLaplace(bodies, rho; ...)."))
end

function (m::PressureLaplace)(systems, wakes,
                              frames::AbstractVector{<:ReferenceFrame},
                              uinf, i_step::Int, dt::Real)
    dt > 0 || throw(ArgumentError("PressureLaplace requires a positive runtime dt; got $(dt)."))
    systems_tuple = _systems_tuple(systems)
    length(systems_tuple) == length(m.b) || throw(ArgumentError(
        "PressureLaplace was constructed for $(length(m.b)) bodies, got $(length(systems_tuple)) bodies."))
    for (i, body) in enumerate(systems_tuple)
        _pressure_laplace_body!(m, body, i, i_step, dt)
    end
end

function build_pressure_preconditioner(pc::NoPressurePreconditioner, Ls::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}})
    return pc
end

function build_pressure_preconditioner(::NoPressurePreconditioner, L::SparseArrays.SparseMatrixCSC{Float64, Int})
    return nothing
end

function build_pressure_preconditioner(::JacobiPressurePreconditioner, Ls::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}})
    states = PressureJacobiState{Vector{Float64}}[]
    sizehint!(states, length(Ls))
    for L in Ls
        push!(states, build_pressure_preconditioner(JacobiPressurePreconditioner(), L))
    end
    return JacobiPressurePreconditioner(states)
end

function build_pressure_preconditioner(pc::IncompleteCholeskyPressurePreconditioner, Ls::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}})
    throw(ArgumentError("IncompleteCholeskyPressurePreconditioner is reserved but not implemented yet."))
end

function build_pressure_preconditioner(pc::AMGPressurePreconditioner, Ls::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}})
    throw(ArgumentError("AMGPressurePreconditioner is reserved but not implemented yet."))
end

function build_pressure_preconditioner(::JacobiPressurePreconditioner, L::SparseArrays.SparseMatrixCSC{Float64, Int})
    d = Vector(LA.diag(L))
    invdiag = similar(d, Float64)
    @inbounds for i in eachindex(d)
        invdiag[i] = iszero(d[i]) ? 1.0 : inv(d[i])
    end
    return PressureJacobiState(invdiag)
end

build_pressure_preconditioner(pc::AbstractPressurePreconditioner, L, cache) =
    build_pressure_preconditioner(pc, L)

function rebuild_pressure_preconditioner!(pc::NoPressurePreconditioner,
                                          L::SparseArrays.SparseMatrixCSC{Float64, Int}, i_body::Int)
    return pc
end

function rebuild_pressure_preconditioner!(pc::JacobiPressurePreconditioner,
                                          L::SparseArrays.SparseMatrixCSC{Float64, Int}, i_body::Int)
    state = pc.states[i_body]
    d = Vector(LA.diag(L))
    length(state.invdiag) == length(d) || throw(ArgumentError(
        "Jacobi pressure preconditioner size changed from $(length(state.invdiag)) to $(length(d)). Reconstruct PressureLaplace for changed panel counts."))
    @inbounds for i in eachindex(d)
        state.invdiag[i] = iszero(d[i]) ? 1.0 : inv(d[i])
    end
    return state
end

pressure_preconditioner_argument(pc::JacobiPressurePreconditioner, i_body::Int) = pc.states[i_body]
pressure_preconditioner_argument(::NoPressurePreconditioner, i_body::Int) = nothing

function _pressure_laplace_body!(m::PressureLaplace, body::AbstractBody,
                                 i_body::Int, i_step::Int, dt::Real)
    m.reference_panel <= body.ncells || throw(ArgumentError(
        "reference_panel=$(m.reference_panel) exceeds body.ncells=$(body.ncells)."))

    # Keep normals and control points current; bodies may have moved this step.
    calc_normals!(body)
    calc_controlpoints!(body; off=abs(body.CPoffset))

    # Panel count changes would invalidate all preallocated arrays.
    topology_changed = body.ncells != length(m.b[i_body]) ||
                       size(m.velocity_dot[i_body]) != size(body.velocity) ||
                       size(m.u_inertial[i_body]) != size(body.velocity)
    topology_changed && throw(ArgumentError(
        "PressureLaplace does not support panel-count changes after construction. Reconstruct the monitor for the new body sizes."))

    # Recompute L only when geometry has changed (rigid motion leaves L unchanged).
    sig = _pressure_geometry_signature(body, m.reference_panel, m.edges[i_body])
    rebuild = !m.cache || m.geometry_signature[i_body] != sig

    if rebuild
        m.L[i_body] = _assemble_pressure_laplacian(body, m.reference_panel, m.edges[i_body])
        rebuild_pressure_preconditioner!(m.preconditioner, m.L[i_body], i_body)
        m.workspace[i_body] = Krylov.krylov_workspace(Val(:cg), m.L[i_body], m.b[i_body])
        m.geometry_signature[i_body] = sig
    end

    # velocity_dot currently holds -u_old; this call turns it into (u_new - u_old)/dt.
    _pressure_velocity_dot!(m, body, i_body, dt)
    _pressure_rhs!(m, body, i_body)   # build b from material acceleration
    _pressure_solve!(m, i_body)       # solve L p = b with CG
    copyto!(body.P, m.p[i_body])
    # Store -u_new so the next call can form the finite difference.
    _pressure_store_negative_velocity!(m, body, i_body)

    if m.verbose
        println("\t\tPressureLaplace[step=$(i_step+1)]: rebuild=$(rebuild), panels=$(body.ncells)")
    end

    return body.P
end

function _pressure_velocity_dot!(m::PressureLaplace, body::AbstractBody, i_body::Int, dt::Real)
    # On entry velocity_dot holds -u_inertial_old (written by the previous call to
    # _pressure_store_negative_velocity!, or zero on the very first call).
    # Adding u_inertial_new = tangent(body.velocity) + body.velocity_kinematic
    # and dividing by dt gives the panel-following rate
    # d/dt[u(x_p,t)] ≈ (u_new - u_old)/dt.
    vdot = m.velocity_dot[i_body]
    u    = body.velocity
    uk   = body.velocity_kinematic
    nrm  = body.normals
    inv_dt = inv(Float64(dt))
    @inbounds for p in 1:body.ncells
        nx, ny, nz = nrm[1, p], nrm[2, p], nrm[3, p]
        ux, uy, uz = u[1, p], u[2, p], u[3, p]
        un = ux*nx + uy*ny + uz*nz
        vdot[1, p] = (vdot[1, p] + (ux - un*nx) + uk[1, p]) * inv_dt
        vdot[2, p] = (vdot[2, p] + (uy - un*ny) + uk[2, p]) * inv_dt
        vdot[3, p] = (vdot[3, p] + (uz - un*nz) + uk[3, p]) * inv_dt
    end

    return m.velocity_dot[i_body]
end

function _pressure_store_negative_velocity!(m::PressureLaplace, body::AbstractBody, i_body::Int)
    # Stash -u_inertial_new so the next call to _pressure_velocity_dot! can form
    # u_inertial_new - u_inertial_old.
    vdot = m.velocity_dot[i_body]
    u    = body.velocity
    uk   = body.velocity_kinematic
    nrm  = body.normals
    @inbounds for p in 1:body.ncells
        nx, ny, nz = nrm[1, p], nrm[2, p], nrm[3, p]
        ux, uy, uz = u[1, p], u[2, p], u[3, p]
        un = ux*nx + uy*ny + uz*nz
        vdot[1, p] = -((ux - un*nx) + uk[1, p])
        vdot[2, p] = -((uy - un*ny) + uk[2, p])
        vdot[3, p] = -((uz - un*nz) + uk[3, p])
    end

    return m.velocity_dot[i_body]
end

function _pressure_solve!(m::PressureLaplace, i_body::Int)
    M = pressure_preconditioner_argument(m.preconditioner, i_body)
    if M === nothing
        Krylov.krylov_solve!(m.workspace[i_body], m.L[i_body], m.b[i_body]; atol=m.atol, rtol=m.rtol, itmax=m.itmax)
    else
        Krylov.krylov_solve!(m.workspace[i_body], m.L[i_body], m.b[i_body]; M, ldiv=true, atol=m.atol, rtol=m.rtol, itmax=m.itmax)
    end
    m.p[i_body] .= m.workspace[i_body].x
    return m.p[i_body]
end

function _pressure_panel_edges(body::AbstractBody)
    # Build the interior-edge table: one column per shared edge, skipping boundary edges
    # (those adjacent to only one panel).  Row layout: (node_a, node_b, panel_i, panel_j).
    edge_to_cells = _calc_edge_to_cells(body.cells)
    nedges = count(refs -> length(refs) == 2, values(edge_to_cells))
    edges = Matrix{Int}(undef, 4, nedges)
    k = 0
    for ((edge_a, edge_b), refs) in edge_to_cells
        length(refs) == 2 || continue
        k += 1
        edges[1, k] = edge_a
        edges[2, k] = edge_b
        edges[3, k] = refs[1][1]
        edges[4, k] = refs[2][1]
    end
    return edges
end

function _pressure_geometry_signature(body::AbstractBody, reference_panel::Int,
                                      edges::Matrix{Int})
    # Hash everything that affects the Laplacian weights: topology (ncells, cells, neighbor),
    # the gauge choice (reference_panel), and the per-edge co-normal FV weights.
    # Rigid translation/rotation changes control-point positions and node positions equally,
    # so the co-normal weights are unchanged and the signature stays the same — no rebuild needed.
    sig = hash((body.ncells, size(body.cells), body.cells, size(body.neighbor), body.neighbor, reference_panel))
    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        w, ell, nu1, nu2, nu3, n1, n2, n3 = _pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)
        sig = hash((w, ell, nu1, nu2, nu3, n1, n2, n3), sig)
    end
    return sig
end

function _assemble_pressure_laplacian(body::AbstractBody, reference_panel::Int)
    return _assemble_pressure_laplacian(body, reference_panel, _pressure_panel_edges(body))
end

function _assemble_pressure_laplacian(body::AbstractBody, reference_panel::Int,
                                      edges::Matrix{Int})
    n = body.ncells
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    sizehint!(rows, 4n)
    sizehint!(cols, 4n)
    sizehint!(vals, 4n)

    # For each shared interior edge between panels i and j, the panel-centered
    # FV weight is w = ℓ (ν_ij · r_ij) / |r_ij|², where ν_ij is the averaged-normal
    # surface co-normal oriented from panel i to panel j. Each edge contributes
    # four symmetric entries: +w on both diagonals, -w on both off-diagonals.
    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        i == j && continue
        w = _pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]
        append!(rows, (i, j, i, j))
        append!(cols, (i, j, j, i))
        append!(vals, (w, w, -w, -w))
    end

    L = SparseArrays.sparse(rows, cols, vals, n, n)
    # Gauge-fix to make L SPD so CG applies; skipped if reference_panel <= 0.
    return reference_panel > 0 ? _pressure_apply_gauge(L, reference_panel) : L
end

function _pressure_apply_gauge(L, reference_panel::Int)
    # Pin p[reference_panel] = reference_pressure by replacing the entire row and
    # column with the identity constraint L_rr = 1, L_rj = L_ir = 0 (j ≠ r).
    # Both row and column are zeroed to preserve symmetry, which is required for CG.
    Iidx, Jidx, V = SparseArrays.findnz(L)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    sizehint!(rows, length(V) + 1)
    sizehint!(cols, length(V) + 1)
    sizehint!(vals, length(V) + 1)
    @inbounds for k in eachindex(V)
        i = Iidx[k]
        j = Jidx[k]
        if i != reference_panel && j != reference_panel
            push!(rows, i); push!(cols, j); push!(vals, V[k])
        end
    end
    push!(rows, reference_panel)
    push!(cols, reference_panel)
    push!(vals, 1.0)
    return SparseArrays.sparse(rows, cols, vals, size(L, 1), size(L, 2))
end

function _pressure_rhs!(m::PressureLaplace, body::AbstractBody, i_body::Int)
    b = m.b[i_body]
    acceleration = m.acceleration[i_body]
    edges = m.edges[i_body]
    rho = m.rho
    reference_panel = m.reference_panel
    reference_pressure = m.reference_pressure
    fill!(b, 0.0)

    # Step 1: compute Du/Dt = ∂u/∂t + (u_t·∇_s)u at every panel control point.
    _pressure_material_acceleration!(m, body, i_body)

    # Step 2: accumulate the edge-integrated source ρ ∮ a_t·ν dℓ into b, matching
    # the panel-centered co-normal FV Laplacian assembled in L. Each side of the
    # midpoint acceleration is projected with the same averaged edge normal used
    # to form ν_ij, so the operator and RHS use the same surface metric.
    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        w, ell, nu1, nu2, nu3, n1, n2, n3 = _pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)

        ai_n = acceleration[1, i] * n1 + acceleration[2, i] * n2 + acceleration[3, i] * n3
        aj_n = acceleration[1, j] * n1 + acceleration[2, j] * n2 + acceleration[3, j] * n3
        ai = (acceleration[1, i] - ai_n * n1) * nu1 +
             (acceleration[2, i] - ai_n * n2) * nu2 +
             (acceleration[3, i] - ai_n * n3) * nu3
        aj = (acceleration[1, j] - aj_n * n1) * nu1 +
             (acceleration[2, j] - aj_n * n2) * nu2 +
             (acceleration[3, j] - aj_n * n3) * nu3
        flux = ell * 0.5 * (ai + aj)
        b[i] += rho * flux
        b[j] -= rho * flux
        # When reference_pressure ≠ 0 the gauge column removal shifts the RHS.
        if !iszero(reference_pressure)
            i == reference_panel && (b[j] += w * reference_pressure)
            j == reference_panel && (b[i] += w * reference_pressure)
        end
    end

    # Step 3: enforce the gauge constraint b[reference_panel] = reference_pressure.
    b[reference_panel] = reference_pressure
    return b
end

function _pressure_edge_conormal_weight(body::AbstractBody, edge_a::Int, edge_b::Int,
                                        i::Int, j::Int)
    r1 = body.controlpoints[1, j] - body.controlpoints[1, i]
    r2 = body.controlpoints[2, j] - body.controlpoints[2, i]
    r3 = body.controlpoints[3, j] - body.controlpoints[3, i]
    r2norm = r1*r1 + r2*r2 + r3*r3
    r2norm > eps(Float64) || throw(ArgumentError(
        "PressureLaplace found coincident control points for panels $(i) and $(j)."))

    e1 = body.nodes[1, edge_b] - body.nodes[1, edge_a]
    e2 = body.nodes[2, edge_b] - body.nodes[2, edge_a]
    e3 = body.nodes[3, edge_b] - body.nodes[3, edge_a]
    ell = sqrt(e1*e1 + e2*e2 + e3*e3)
    ell > eps(Float64) || throw(ArgumentError(
        "PressureLaplace found a degenerate shared edge ($(edge_a), $(edge_b))."))
    t1, t2, t3 = e1 / ell, e2 / ell, e3 / ell

    n1 = body.normals[1, i] + body.normals[1, j]
    n2 = body.normals[2, i] + body.normals[2, j]
    n3 = body.normals[3, i] + body.normals[3, j]
    nmag = sqrt(n1*n1 + n2*n2 + n3*n3)
    nmag > eps(Float64) || throw(ArgumentError(
        "PressureLaplace found opposing panel normals across edge ($(edge_a), $(edge_b))."))
    n1, n2, n3 = n1 / nmag, n2 / nmag, n3 / nmag

    # ν = t_edge × n_avg is tangent to the averaged surface and perpendicular
    # to the shared edge. Orient it from panel i toward panel j.
    nu1 = t2 * n3 - t3 * n2
    nu2 = t3 * n1 - t1 * n3
    nu3 = t1 * n2 - t2 * n1
    numag = sqrt(nu1*nu1 + nu2*nu2 + nu3*nu3)
    numag > eps(Float64) || throw(ArgumentError(
        "PressureLaplace found an edge tangent parallel to the averaged normal across edge ($(edge_a), $(edge_b))."))
    nu1, nu2, nu3 = nu1 / numag, nu2 / numag, nu3 / numag

    proj = nu1*r1 + nu2*r2 + nu3*r3
    if proj < 0
        nu1, nu2, nu3 = -nu1, -nu2, -nu3
        proj = -proj
    end
    proj > sqrt(eps(Float64)) * sqrt(r2norm) || throw(ArgumentError(
        "PressureLaplace found a non-orthogonal panel pair with nonpositive co-normal distance between panels $(i) and $(j)."))

    w = ell * proj / r2norm
    return w, ell, nu1, nu2, nu3, n1, n2, n3
end

function _pressure_material_acceleration!(m::PressureLaplace, body::AbstractBody, i_body::Int)
    n = body.ncells
    acceleration = m.acceleration[i_body]
    velocity_dot = m.velocity_dot[i_body]
    u_inertial = m.u_inertial[i_body]

    # Form u_inertial = tangent(body.velocity) + body.velocity_kinematic. The
    # tangent projection strips the small numerical normal residual from the
    # solver out of velocity_dot and (now obsolete) downstream uses; u_inertial
    # is still maintained for documentation / debugging visibility.
    nrm = body.normals
    uk  = body.velocity_kinematic
    @inbounds for p in 1:n
        nx, ny, nz = nrm[1, p], nrm[2, p], nrm[3, p]
        ux, uy, uz = body.velocity[1, p], body.velocity[2, p], body.velocity[3, p]
        un = ux*nx + uy*ny + uz*nz
        u_inertial[1, p] = (ux - un*nx) + uk[1, p]
        u_inertial[2, p] = (uy - un*ny) + uk[2, p]
        u_inertial[3, p] = (uz - un*nz) + uk[3, p]
    end

    # Start with the panel-following rate d/dt[u_inertial(x_p,t)] stored in velocity_dot.
    copyto!(acceleration, velocity_dot)

    # ∇u_inertial = ∇u_induced + [Ω]_×. The induced part is the analytic
    # multipole Hessian accumulated into body.velocity_gradient by
    # FastMultipole. The kinematic part has spatial Jacobian [Ω]_× because
    # u_kinematic(x) = U_O + Ω × (x - x_O) is rigid-body — its derivative w.r.t.
    # x is the skew-symmetric tensor of the angular velocity Ω. Ω is the net
    # angular velocity accumulated during kinematic_velocity!.
    ω = body.angular_velocity
    ωx, ωy, ωz = ω[1], ω[2], ω[3]
    G_ind = body.velocity_gradient

    # (u_rel · ∇) u_inertial, with u_rel = body.velocity tangent-projected
    # (impermeability is enforced exactly here so the convective term has no
    # spurious normal component from solver residuals).
    @inbounds for i in 1:n
        ur1, ur2, ur3 = body.velocity[1, i], body.velocity[2, i], body.velocity[3, i]
        nx, ny, nz = body.normals[1, i], body.normals[2, i], body.normals[3, i]
        un = ur1 * nx + ur2 * ny + ur3 * nz
        ut1 = ur1 - un * nx
        ut2 = ur2 - un * ny
        ut3 = ur3 - un * nz

        # ∇u_inertial = G_ind[:,:,i] + skew(Ω); add ([Ω]_× · u_t)_k = (Ω × u_t)_k
        # directly to avoid materializing the skew matrix.
        cross_x = ωy * ut3 - ωz * ut2
        cross_y = ωz * ut1 - ωx * ut3
        cross_z = ωx * ut2 - ωy * ut1

        acceleration[1, i] += ut1 * G_ind[1, 1, i] + ut2 * G_ind[1, 2, i] + ut3 * G_ind[1, 3, i] + cross_x
        acceleration[2, i] += ut1 * G_ind[2, 1, i] + ut2 * G_ind[2, 2, i] + ut3 * G_ind[2, 3, i] + cross_y
        acceleration[3, i] += ut1 * G_ind[3, 1, i] + ut2 * G_ind[3, 2, i] + ut3 * G_ind[3, 3, i] + cross_z
    end

    return acceleration
end

function _pressure_edge_length(nodes, a::Int, b::Int)
    dx = nodes[1, b] - nodes[1, a]
    dy = nodes[2, b] - nodes[2, a]
    dz = nodes[3, b] - nodes[3, a]
    return sqrt(dx*dx + dy*dy + dz*dz)
end

function _pressure_controlpoint_distance(controlpoints, i::Int, j::Int)
    dx = controlpoints[1, j] - controlpoints[1, i]
    dy = controlpoints[2, j] - controlpoints[2, i]
    dz = controlpoints[3, j] - controlpoints[3, i]
    return sqrt(dx*dx + dy*dy + dz*dz)
end


"""
    WingNormalization(rho, Sref, Lref)

Default normalization for `ForceMonitor`: divides force by `0.5 ρ |U∞|² Sref`
and moment by `0.5 ρ |U∞|² Sref Lref`.
"""
struct WingNormalization{TF} <: AbstractMonitor
    rho::TF
    Sref::TF
    Lref::TF
end

function (n::WingNormalization)(CF, CM, systems, frames, uinf)
    qinf = n.rho / 2 * (uinf[1]^2 + uinf[2]^2 + uinf[3]^2)
    return CF / (qinf * n.Sref), CM / (qinf * n.Sref * n.Lref)
end

"""
    NoNormalization()

Pass-through normalization for `ForceMonitor`: returns dimensional force and
moment unchanged.
"""
struct NoNormalization end

function (::NoNormalization)(CF, CM, systems, frames, uinf)
    return CF, CM
end

"""
    RotorNormalization(rho, D, i_frame)

Rotor thrust/torque coefficient normalization.  At each call, shaft speed `n`
(rev/s) is read from `frames[i_frame].ω / 2π`.  Force is divided by `ρ n² D⁴`
(giving CT) and moment by `ρ n² D⁵` (giving CQ), where `D` is rotor diameter.
"""
struct RotorNormalization{TF}
    rho::TF
    D::TF       # rotor diameter
    i_frame::Int
end

function (c::RotorNormalization)(CF, CM, systems, frames, uinf)
    n = frames[c.i_frame].ω / (2 * pi)   # rev/s
    ref = c.rho * n^2 * c.D^4
    return CF / ref, CM / (ref * c.D)
end

struct RotorNormalization2{TF}
    rho::TF
    D::TF       # rotor diameter
    i_frame::Int
end

function (c::RotorNormalization2)(CF, CM, systems, frames, uinf)
    n = frames[c.i_frame].ω # rad/s
    ref = c.rho * pi * c.D^2 * 0.25 * (n * c.D / 2)^2
    return CF / ref, CM / (ref * c.D*0.5)
end

"""
    ForceMonitor(nt, i_system; i_frame=-1, rho=1.0, Sref=1.0, Lref=1.0, TF=Float64, normalization=WingNormalization(...), correct_kuttacondition=true, verbose=false)

Storage monitor for integrated force and moment coefficient histories over `nt`
simulation steps.  At each step the monitor calls `calcfield_F!` to populate
`body.F` (reading the pressure field `body.P` set by an earlier
`PressureBernoulli`), then integrates force and moment on `systems[i_system]`
in the coordinate system of `frames[i_frame]`. A `PressureBernoulli` must
therefore appear **before** this monitor in the `monitors` tuple passed to
`simulate!`.

`normalization` is called as `normalization(CF, CM, systems, frames, uinf)` and
must return a `(CF_norm, CM_norm)` tuple.  The default is `WingNormalization`
which divides by `0.5 ρ |U∞|² Sref` (and `… Lref` for moments).

If `verbose=true`, the normalized CF and CM for each step are printed to stdout
with a single `\\t` indent.
"""
struct ForceMonitor{TF, TN} <: AbstractMonitor
    force::Matrix{TF}
    moment::Matrix{TF}
    i_system::Int
    i_frame::Int
    normalization::TN
    correct_kuttacondition::Bool
    verbose::Bool
end

monitor_requires(::ForceMonitor) = (:P,)
monitor_provides(::ForceMonitor) = (:F,)

"""
    KuttaJoukowskiForce(body, nt, i_system; rho=1.225, backend=DirectBackend(),
                        i_frame=-1, TF=Float64,
                        normalization=NoNormalization(), verbose=false)

Diagnostic force monitor that integrates the Kutta–Joukowski force on each panel
edge and reports the **force on the body**:
``F = ρ Σᵢ Σⱼ γᵢ \\, (Vᵢⱼ × Δsᵢⱼ) = -ρ Σᵢ Σⱼ γᵢ \\, (Δsᵢⱼ × Vᵢⱼ)``
summed over panels `i` and their three edges `j`, where `γᵢ` is the panel
circulation (column of `body.strength` named `"gamma"` or `"mu"`), `Δsᵢⱼ` is the
directed edge vector along the filament, and `Vᵢⱼ` is the inertial-frame fluid
velocity (induced + freestream) at the edge midpoint, evaluated via a
`FastMultipole.ProbeSystem`. The sign matches the pressure-integral
`ForceMonitor` convention (force on the body, not on the fluid).

Provides an independent cross-check against the pressure-integral force
returned by `ForceMonitor`/`calcfield_F!`. The body must use a `ConstantDoublet`
or `VortexRing` kernel.

If `i_frame < 0`, force is reported in global coordinates. Otherwise the summed
force is rotated into the coordinate system of `frames[i_frame]`, matching
`ForceMonitor`.

If `verbose=true`, the normalized CF for each step is printed to stdout with a
single `\\t` indent.
"""
struct KuttaJoukowskiForce{TF, TB, TN} <: AbstractMonitor
    force::Matrix{TF}
    i_system::Int
    i_frame::Int
    i_strength::Int
    rho::TF
    backend::TB
    probes::FastMultipole.ProbeSystem{TF}
    edge_node_a::Vector{Int}
    edge_node_b::Vector{Int}
    panel_of_probe::Vector{Int}
    normalization::TN
    verbose::Bool
end

function KuttaJoukowskiForce(body::AbstractBody, nt::Int, i_system::Int;
                              rho::Real=1.225,
                              backend::AbstractBackend=DirectBackend(),
                              i_frame::Int=-1,
                              TF=Float64,
                              normalization=NoNormalization(),
                              verbose::Bool=false)

    names = strength_names(body)
    i_strength = something(findfirst(==("gamma"), names),
                            findfirst(==("mu"), names),
                            0)
    i_strength == 0 && throw(ArgumentError(
        "KuttaJoukowskiForce requires a body with a ConstantDoublet or "*
        "VortexRing kernel; got strengths $(names)."))

    cells = body.cells
    @assert size(cells, 1) == 3 "KuttaJoukowskiForce assumes triangular panels (3 nodes per cell), got $(size(cells, 1))."
    ncells = body.ncells
    n_probes = 3 * ncells

    edge_node_a = Vector{Int}(undef, n_probes)
    edge_node_b = Vector{Int}(undef, n_probes)
    panel_of_probe = Vector{Int}(undef, n_probes)
    @inbounds for i in 1:ncells
        n1 = cells[1, i]; n2 = cells[2, i]; n3 = cells[3, i]
        k = 3*(i-1)
        edge_node_a[k+1] = n1; edge_node_b[k+1] = n2
        edge_node_a[k+2] = n2; edge_node_b[k+2] = n3
        edge_node_a[k+3] = n3; edge_node_b[k+3] = n1
        panel_of_probe[k+1] = i
        panel_of_probe[k+2] = i
        panel_of_probe[k+3] = i
    end

    probes = FastMultipole.ProbeSystem(n_probes, TF)
    force = zeros(TF, 3, nt)

    return KuttaJoukowskiForce{TF, typeof(backend), typeof(normalization)}(
        force, i_system, i_frame, i_strength, TF(rho), backend, probes,
        edge_node_a, edge_node_b, panel_of_probe, normalization, verbose)
end

function (m::KuttaJoukowskiForce{TF})(systems, wakes,
                                       frames::AbstractVector{<:ReferenceFrame},
                                       uinf, i_step::Int, dt::Real) where {TF}
    body = systems[m.i_system]

    # 1. Update probe positions to edge midpoints; reset accumulators.
    zero_v = zero(SVector{3, TF})
    zero_h = zero(SMatrix{3, 3, TF, 9})
    @inbounds for k in eachindex(m.edge_node_a)
        a = m.edge_node_a[k]; b = m.edge_node_b[k]
        m.probes.position[k] = SVector{3, TF}(
            0.5 * (body.nodes[1, a] + body.nodes[1, b]),
            0.5 * (body.nodes[2, a] + body.nodes[2, b]),
            0.5 * (body.nodes[3, a] + body.nodes[3, b]),
        )
        m.probes.gradient[k] = zero_v
        m.probes.scalar_potential[k] = zero(TF)
        m.probes.hessian[k] = zero_h
    end

    # 2. Collect all sources: every body, plus every wake's source list.
    wake_sources = _collect_wake_sources(wakes)
    all_sources  = (systems..., wake_sources...)

    # 3. Compute induced velocity at probes.
    influence!((m.probes,), all_sources, m.backend;
                precalc=false,
                scalar_potential=false,
                gradient=true,
                hessian=(false,))

    # 4. Add freestream manually (no apply_freestream! method for ProbeSystem).
    uinf_sv = SVector{3, TF}(uinf[1], uinf[2], uinf[3])
    @inbounds for k in eachindex(m.probes.gradient)
        m.probes.gradient[k] += uinf_sv
    end

    # 5. Sum Kutta–Joukowski contributions. Per the K–J theorem the force on the
    # body from a vortex filament segment of length Δs (along the filament) and
    # circulation γ in fluid velocity V is F_body = ρ γ (V × Δs) = -ρ γ (Δs × V).
    # The opposite sign (ρ γ Δs × V) is the force on the fluid; we report
    # body-side force here to match the pressure-integral ForceMonitor convention.
    Fsum = zero(SVector{3, TF})
    @inbounds for k in eachindex(m.edge_node_a)
        a = m.edge_node_a[k]; b = m.edge_node_b[k]
        Δs = SVector{3, TF}(
            body.nodes[1, b] - body.nodes[1, a],
            body.nodes[2, b] - body.nodes[2, a],
            body.nodes[3, b] - body.nodes[3, a],
        )
        γ = body.strength[m.panel_of_probe[k], m.i_strength]
        V = m.probes.gradient[k]
        Fsum -= m.rho * γ * cross(Δs, V)
    end

    if m.i_frame >= 0
        _, R_f2g = frame_global_transform(frames, m.i_frame)
        Fsum = transpose(R_f2g) * Fsum
    end

    # 6. Optional normalization (moment is not computed; pass zero).
    Mzero = zero(SVector{3, TF})
    CF, _ = m.normalization(Fsum, Mzero, systems, frames, uinf)
    m.force[:, i_step + 1] .= CF

    if m.verbose
        println("\t\tKuttaJoukowskiForce[i_system=$(m.i_system), i_frame=$(m.i_frame), step=$(i_step+1)]:")
        println("\t\t\tCF = ($(round(CF[1], sigdigits=4)), $(round(CF[2], sigdigits=4)), $(round(CF[3], sigdigits=4)))")
    end
end


function ForceMonitor(nt::Int, i_system::Int;
                       i_frame=-1, TF=Float64,
                       normalization=WingNormalization(TF(rho), TF(Sref), TF(Lref)),
                       correct_kuttacondition::Bool=true,
                       verbose::Bool=false)
    force = zeros(TF, 3, nt)
    moment = zeros(TF, 3, nt)
    return ForceMonitor{TF, typeof(normalization)}(force, moment, i_system, i_frame, normalization, correct_kuttacondition, verbose)
end

function (monitor::ForceMonitor)(systems, wakes,
                                  frames::AbstractVector{<:ReferenceFrame},
                                  uinf, i_step::Int, dt::Real)
    systems_tuple = _systems_tuple(systems)
    body = systems_tuple[monitor.i_system]

    # Populate body.F from body.P. Reset first so multiple ForceMonitors in the same
    # chain each see only their own pressure, not accumulated values from earlier monitors.
    fill!(body.F, zero(eltype(body.F)))
    calcfield_F!(body; correct_kuttacondition=monitor.correct_kuttacondition)

    # total force in global frame
    Ftot = calcfield_Ftot(body)
    Fvec = FastMultipole.SVector{3}(Ftot[1], Ftot[2], Ftot[3])

    if monitor.i_frame < 0
        # global frame: moment about the origin
        Mtot = calcfield_Mtot(body, zeros(3))
        Mvec = FastMultipole.SVector{3}(Mtot[1], Mtot[2], Mtot[3])
    else
        # frame-local: moment about frame origin, rotated into frame axes
        origin_global, R_f2g = frame_global_transform(frames, monitor.i_frame)
        Mtot = calcfield_Mtot(body, collect(origin_global))
        R_g2f = transpose(R_f2g)
        Fvec = R_g2f * Fvec
        Mvec = R_g2f * FastMultipole.SVector{3}(Mtot[1], Mtot[2], Mtot[3])
    end

    # normalise to coefficients
    CF_norm, CM_norm = monitor.normalization(Fvec, Mvec, systems_tuple, frames, uinf)
    monitor.force[:, i_step + 1] .= CF_norm
    monitor.moment[:, i_step + 1] .= CM_norm

    if monitor.verbose
        println("\t\tForceMonitor[i_system=$(monitor.i_system), step=$(i_step+1)]:")
        println("\t\t\tCF = ($(round(CF_norm[1], sigdigits=4)), $(round(CF_norm[2], sigdigits=4)), $(round(CF_norm[3], sigdigits=4)))")
        println("\t\t\tCM = ($(round(CM_norm[1], sigdigits=4)), $(round(CM_norm[2], sigdigits=4)), $(round(CM_norm[3], sigdigits=4)))")
    end
end
