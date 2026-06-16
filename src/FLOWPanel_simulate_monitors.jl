#--- monitor dependency traits ---#

abstract type AbstractMonitor end

# Each monitor declares which per-step fields it populates (`monitor_provides`)
# and which it consumes (`monitor_requires`). `audit_monitors` checks at the
# start of `simulate!` that every requirement is met by an earlier monitor
# in the tuple, surfacing ordering bugs before the time loop starts.
#
# Field symbols correspond to monitor-owned per-body arrays:
#   :P  → pressure
#   :F  → distributed surface force
#
# Add a new monitor's contract by defining methods on these two traits.

monitor_provides(::Any) = ()
monitor_requires(::Any) = ()

# Monitors that need ∇u at body control points return true here so
# simulate! can flip body.needs_velocity_gradient[] before the time loop.
monitor_requires_body_hessian(::Any) = false

# Monitors that need volumetric induced vorticity at body control points return
# true here so simulate! can request FastMultipole extra_outputs=3.
monitor_requires_induced_vorticity(::Any) = false

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
                    ":$(r) to be produced, but no earlier monitor in the " *
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

mutable struct MonitorContext
    fields::Dict{Tuple{Symbol, Int}, Any}
    time::Float64
    has_time::Bool
end

MonitorContext() = MonitorContext(Dict{Tuple{Symbol, Int}, Any}(), 0.0, false)

function monitor_set_time!(ctx::MonitorContext, t::Real)
    ctx.time = Float64(t)
    ctx.has_time = true
    return ctx
end

monitor_time(ctx::MonitorContext, i_step::Int, dt::Real) =
    ctx.has_time ? ctx.time : Float64(i_step * dt)

function monitor_register!(ctx::MonitorContext, field::Symbol, i_body::Integer, value)
    ctx.fields[(field, Int(i_body))] = value
    return value
end

function monitor_field(ctx::MonitorContext, field::Symbol, i_body::Integer)
    key = (field, Int(i_body))
    haskey(ctx.fields, key) || throw(ArgumentError(
        "Monitor field :$(field) for body $(i_body) is unavailable. " *
        "Check monitor ordering and `monitor_provides`/`monitor_requires`."))
    return ctx.fields[key]
end

function _run_monitor!(m, ctx::MonitorContext, systems, wakes,
                       frames::AbstractVector{<:ReferenceFrame},
                       uinf, i_step::Int, dt::Real, t=nothing)
    !isnothing(t) && monitor_set_time!(ctx, t)
    m(systems, wakes, frames, uinf, i_step, dt)
    _register_monitor_outputs!(ctx, m, _systems_tuple(systems))
    return nothing
end

_register_monitor_outputs!(ctx::MonitorContext, m, systems_tuple::Tuple) = nothing

write_vtk_fields!(vtk, monitor, body, i_system::Int, i_step::Int) = nothing
write_vtk_fields!(vtk, monitor, body, i_system::Int, i_step::Int,
                  field_names::VTKFieldNameAllocator, i_monitor::Int) = nothing

write_monitor_csv!(monitor, dir::AbstractString, name::AbstractString,
                   i_monitor::Int, ctx::MonitorContext, systems_tuple::Tuple,
                   i_step::Int, dt::Real; overwrite::Bool=false) = nothing

function _monitor_csv_path(dir::AbstractString, name::AbstractString,
                           i_monitor::Int, stem::AbstractString, i_system::Int)
    filename = "$(name)_monitor$(lpad(i_monitor, 2, '0'))_$(stem)_system$(i_system).csv"
    return joinpath(dir, filename)
end

function _monitor_csv_open(writer::Function, path::AbstractString, i_step::Int,
                           overwrite::Bool, header::AbstractString)
    mode = (overwrite || i_step == 0) ? "w" : "a"
    dir = dirname(path)
    !isdir(dir) && mkpath(dir)
    open(path, mode) do io
        mode == "w" && println(io, header)
        writer(io)
    end
    return nothing
end

_monitor_csv_bool(x::Bool) = x ? "true" : "false"

"""
    PressureBernoulli(rho; unsteady=false, correct_kuttacondition=true, clip=nothing,
                      backend=FastMultipoleBackend(expansion_order=10,
                                                   multipole_acceptance=0.4,
                                                   leaf_size=100))

Monitor that owns pressure arrays for every body in the simulation by
evaluating the Bernoulli equation each step. With `unsteady=false` (default) the steady
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
    pressure::Vector{Vector{Float64}}
    phi_dot::Vector{Vector{Float64}}
    potential_history::Vector{Vector{Float64}}
    probes::Vector{FastMultipole.ProbeSystem{Float64}}
    vtk_fields::Tuple{Vararg{Symbol}}
    file::Bool
end

monitor_provides(::PressureBernoulli) = (:P,)

function PressureBernoulli(rho::Real; unsteady::Bool=false,
                          correct_kuttacondition::Bool=true,
                          clip=nothing,
                          backend=FastMultipoleBackend(expansion_order=10,
                                                       multipole_acceptance=0.4,
                                                       leaf_size=100),
                          vtk_fields::Tuple{Vararg{Symbol}}=(:pressure,),
                          file::Bool=true)
    return PressureBernoulli{typeof(rho), typeof(clip), typeof(backend)}(
        rho, unsteady, correct_kuttacondition, clip, backend,
        Vector{Float64}[], Vector{Float64}[], Vector{Float64}[],
        FastMultipole.ProbeSystem{Float64}[], vtk_fields, file)
end

function (m::PressureBernoulli)(systems, wakes,
                               frames::AbstractVector{<:ReferenceFrame},
                               uinf, i_step::Int, dt::Real)
    systems_tuple = _systems_tuple(systems)
    Uinf_mag = norm(uinf)
    scalar_sources = m.unsteady ? _bernoulli_scalar_sources(systems_tuple, wakes) : ()
    _pressure_bernoulli_ensure_storage!(m, systems_tuple)
    for (i_body, body) in enumerate(systems_tuple)
        pressure = m.pressure[i_body]
        fill!(pressure, 0.0)
        phi_dot = m.unsteady ?
            _pressure_bernoulli_phi_dot!(m, body, i_body, scalar_sources, uinf, dt) :
            nothing
        calcfield_P!(pressure, body, body.velocity, Uinf_mag, m.rho, phi_dot;
                     correct_kuttacondition=m.correct_kuttacondition,
                     clip=m.clip)
    end
end

function _register_monitor_outputs!(ctx::MonitorContext, m::PressureBernoulli, systems_tuple::Tuple)
    for i in eachindex(systems_tuple)
        monitor_register!(ctx, :P, i, m.pressure[i])
    end
    return nothing
end

function write_vtk_fields!(vtk, m::PressureBernoulli, body, i_system::Int, i_step::Int)
    (:pressure in m.vtk_fields) || return nothing
    i_system <= length(m.pressure) || return nothing
    length(m.pressure[i_system]) == body.ncells || return nothing
    vtk["gauge pressure", VTKCellData()] = m.pressure[i_system]
    return nothing
end

function write_vtk_fields!(vtk, m::PressureBernoulli, body, i_system::Int, i_step::Int,
                           field_names::VTKFieldNameAllocator, i_monitor::Int)
    (:pressure in m.vtk_fields) || return nothing
    i_system <= length(m.pressure) || return nothing
    pressure = m.pressure[i_system]
    length(pressure) == body.ncells || return nothing
    vtk[_vtk_monitor_field_name!(field_names, "gauge pressure", m, i_monitor), VTKCellData()] = pressure
    return nothing
end

function _pressure_bernoulli_ensure_storage!(m::PressureBernoulli, systems_tuple::Tuple)
    nbodies = length(systems_tuple)
    if length(m.pressure) != nbodies ||
       length(m.phi_dot) != nbodies ||
       length(m.potential_history) != nbodies ||
       length(m.probes) != nbodies
        m.pressure = Vector{Float64}[]
        m.phi_dot = Vector{Float64}[]
        m.potential_history = Vector{Float64}[]
        m.probes = FastMultipole.ProbeSystem{Float64}[]
        sizehint!(m.pressure, nbodies)
        sizehint!(m.phi_dot, nbodies)
        sizehint!(m.potential_history, nbodies)
        sizehint!(m.probes, nbodies)
        for body in systems_tuple
            push!(m.pressure, zeros(Float64, body.ncells))
            push!(m.phi_dot, zeros(Float64, body.ncells))
            push!(m.potential_history, zeros(Float64, body.ncells))
            push!(m.probes, FastMultipole.ProbeSystem(body.ncells, Float64))
        end
        return nothing
    end

    for (i, body) in enumerate(systems_tuple)
        if length(m.pressure[i]) != body.ncells ||
           length(m.phi_dot[i]) != body.ncells ||
           length(m.potential_history[i]) != body.ncells ||
           length(m.probes[i].scalar_potential) != body.ncells
            m.pressure[i] = zeros(Float64, body.ncells)
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
    PressureLaplace(bodies, rho; unsteady=false, atol=1e-8, rtol=1e-8, itmax=1000,
                    preconditioner=JacobiPressurePreconditioner(),
                    reference_panel=1, reference_pressure=0.0,
                    rebuild_every_step=false, verbose=false, gradient_mode=:raw_hessian,
                    acceleration_form=:material_derivative)

Monitor that owns pressure by solving a sparse panel-centered surface
pressure Poisson equation. The monitor uses `velocity_dot` as a rolling
time-difference buffer of the tangent-projected **inertial** fluid velocity
`u_inertial = tangent(body.velocity) + body.velocity_kinematic`: between calls it
stores `-u_inertial_old`, then during the call it becomes
`(u_inertial_new - u_inertial_old) / dt`, which is the panel-following rate
`d/dt[u(x_p(t), t)]`. With `unsteady=true`, this finite-difference term is
included in the RHS. With `unsteady=false` (default), the rolling buffer is
still updated but the RHS uses only the convective part of
`Du/Dt = ∂u/∂t + (u_rel · ∇)u`.
`velocity_dot` is initialized to zero, so on the first call the unsteady term
is `u_inertial / dt` rather than a true finite difference; if this warm-up
transient matters, pre-populate
`monitor.velocity_dot[i]` with `-(tangent(body.velocity) .+ body.velocity_kinematic)`
before the first call.

The diagnostic material-acceleration reconstruction has two gradient modes. The
default `gradient_mode=:raw_hessian` uses the analytic spatial Jacobian
`∇u_inertial = body.velocity_gradient + [Ω]_×`, where `body.velocity_gradient`
is the multipole Hessian populated by FastMultipole during per-step
`influence!` calls. The opt-in `gradient_mode=:surface_velocity` reconstructs
`∇ₛu_inertial` directly from the final surface velocity field; this is useful
for Dirichlet lifting bodies whose `body.velocity` includes a postprocessed
`∇ₛμ` contribution that is not represented in the raw panel Hessian.

The pressure RHS has two velocity-only acceleration forms. The default
`acceleration_form=:material_derivative` uses the direct
`∂t u + (u · ∇)u` edge form. `acceleration_form=:lamb_vector` uses the Lamb
decomposition `(u · ∇)u = ∇(|u|²/2) + ω × u` with the same tangent relative
velocity used by the material-derivative edge form, where `ω` is the volumetric
induced vorticity accumulated in `body.induced_vorticity` by FastMultipole
`extra_outputs=3`. Both forms use only velocity and its derivatives; neither
requires a scalar potential.

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
- `unsteady::Bool` — include the rolling finite-difference velocity term in the RHS
- `atol::Float64` — absolute residual tolerance for the CG solve
- `rtol::Float64` — relative residual tolerance for the CG solve
- `itmax::Int` — maximum CG iterations per call
- `preconditioner::TP` — preconditioner applied during the CG solve (e.g. `JacobiPressurePreconditioner`)
- `reference_panel::Int` — index of the panel whose pressure is pinned to `reference_pressure` (gauge condition)
- `reference_pressure::Float64` — pressure value assigned to the reference panel
- `rebuild_every_step::Bool` — when `true`, rebuild the sparse operator, preconditioner, and CG workspace on every call; by default they are reused
- `verbose::Bool` — when `true`, prints per-step rebuild status and panel count to stdout
- `gradient_mode::Symbol` — `:raw_hessian` for the current Hessian path or `:surface_velocity` for a surface reconstruction of the final velocity field
- `acceleration_form::Symbol` — `:material_derivative` for the direct acceleration edge form or `:lamb_vector` for the Lamb-vector decomposition
- `velocity_dot::Vector{Matrix{Float64}}` — rolling buffer; initialized to zero, then between calls holds `-u_inertial_old` with `u_inertial = tangent(body.velocity) + body.velocity_kinematic` (3 × ncells per body); during a call becomes `(u_inertial_new - u_inertial_old) / dt`
- `u_inertial::Vector{Matrix{Float64}}` — scratch buffer for the inertial fluid velocity `tangent(body.velocity) + body.velocity_kinematic` (3 × ncells per body)
- `surface_velocity_gradient::Vector{Array{Float64,3}}` — scratch buffer for `∇ₛu_inertial` in `:surface_velocity` mode (3 × 3 × ncells per body)
- `L::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}}` — sparse FV surface Laplacian per body, gauge-fixed at `reference_panel`
- `pressure_operator::Vector{PressureLaplacianOperator}` — matrix-free FV surface Laplacian per body for CG matvecs
- `b::Vector{Vector{Float64}}` — RHS vector per body (length ncells); rebuilt each call from the selected acceleration form
- `p::Vector{Vector{Float64}}` — owned pressure solution vector per body (length ncells)
- `acceleration::Vector{Matrix{Float64}}` — material acceleration `Du/Dt` scratch buffer (3 × ncells per body)
- `tangential::Vector{Matrix{Float64}}` — tangential projection of `acceleration` (3 × ncells per body)
- `edges::Vector{Matrix{Int}}` — shared interior edges per body (4 × nedges); rows are `(node_a, node_b, panel_i, panel_j)`
- `workspace::Vector{Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}}` — preallocated Krylov CG workspace per body
"""
struct PressureLaplacianOperator
    n::Int
    reference_panel::Int
    row_diagonal::Vector{Float64}
    row_neighbors::Vector{Vector{Int}}
    row_weights::Vector{Vector{Float64}}
end

Base.size(A::PressureLaplacianOperator) = (A.n, A.n)
Base.size(A::PressureLaplacianOperator, d::Int) = d <= 2 ? A.n : 1
Base.eltype(::Type{PressureLaplacianOperator}) = Float64
Base.eltype(::PressureLaplacianOperator) = Float64
LA.issymmetric(::PressureLaplacianOperator) = true
LA.ishermitian(::PressureLaplacianOperator) = true
Base.transpose(A::PressureLaplacianOperator) = A
Base.adjoint(A::PressureLaplacianOperator) = A

(A::PressureLaplacianOperator)(y, x, α, β) = _pressure_laplacian_mul!(y, x, α, β, A)
LA.mul!(y, A::PressureLaplacianOperator, x) =
    _pressure_laplacian_mul!(y, x, one(Float64), zero(Float64), A)
LA.mul!(y, A::PressureLaplacianOperator, x, α, β) =
    _pressure_laplacian_mul!(y, x, α, β, A)

mutable struct PressureLaplace{TP} <: AbstractMonitor
    rho::Float64
    unsteady::Bool
    atol::Float64
    rtol::Float64
    itmax::Int
    preconditioner::TP
    reference_panel::Int
    reference_pressure::Float64
    rebuild_every_step::Bool
    verbose::Bool
    gradient_mode::Symbol
    acceleration_form::Symbol
    velocity_dot::Vector{Matrix{Float64}}
    L::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}}
    pressure_operator::Vector{PressureLaplacianOperator}
    b::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}
    acceleration::Vector{Matrix{Float64}}
    tangential::Vector{Matrix{Float64}}
    edges::Vector{Matrix{Int}}
    workspace::Vector{Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}}
    u_inertial::Vector{Matrix{Float64}}
    surface_velocity_gradient::Vector{Array{Float64,3}}
    vtk_fields::Tuple{Vararg{Symbol}}
    last_rebuild::Vector{Bool}
    file::Bool
end

monitor_provides(::PressureLaplace) = (:P,)
monitor_requires_body_hessian(m::PressureLaplace) =
    m.gradient_mode == :raw_hessian && m.acceleration_form == :material_derivative
monitor_requires_induced_vorticity(m::PressureLaplace) = m.acceleration_form == :lamb_vector

function PressureLaplace(bodies, rho::Real;
                         unsteady::Bool=false,
                         atol::Real=1e-8,
                         rtol::Real=1e-8,
                         itmax::Integer=1000,
                         preconditioner::AbstractPressurePreconditioner=JacobiPressurePreconditioner(),
                         reference_panel::Integer=1,
                         reference_pressure::Real=0.0,
                         rebuild_every_step::Bool=false,
                         verbose::Bool=false,
                         gradient_mode::Symbol=:raw_hessian,
                         acceleration_form::Symbol=:material_derivative,
                         vtk_fields::Tuple{Vararg{Symbol}}=(:pressure,),
                         file::Bool=true)
    reference_panel >= 1 || throw(ArgumentError("reference_panel must be at least 1; got $(reference_panel)."))
    gradient_mode in (:raw_hessian, :surface_velocity) || throw(ArgumentError(
        "gradient_mode must be :raw_hessian or :surface_velocity; got $(gradient_mode)."))
    acceleration_form in (:material_derivative, :lamb_vector) || throw(ArgumentError(
        "acceleration_form must be :material_derivative or :lamb_vector; got $(acceleration_form)."))
    systems_tuple = _systems_tuple(bodies)
    isempty(systems_tuple) && throw(ArgumentError("PressureLaplace requires at least one body."))
    nbodies = length(systems_tuple)
    velocity_dot = Matrix{Float64}[]
    Ls = SparseArrays.SparseMatrixCSC{Float64, Int}[]
    pressure_operators = PressureLaplacianOperator[]
    bs = Vector{Float64}[]
    ps = Vector{Float64}[]
    acceleration = Matrix{Float64}[]
    tangential = Matrix{Float64}[]
    edges = Matrix{Int}[]
    workspace = Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}[]
    u_inertial = Matrix{Float64}[]
    surface_velocity_gradient = Array{Float64,3}[]
    last_rebuild = Bool[]
    sizehint!(velocity_dot, nbodies)
    sizehint!(Ls, nbodies)
    sizehint!(pressure_operators, nbodies)
    sizehint!(bs, nbodies)
    sizehint!(ps, nbodies)
    sizehint!(acceleration, nbodies)
    sizehint!(tangential, nbodies)
    sizehint!(edges, nbodies)
    sizehint!(workspace, nbodies)
    sizehint!(u_inertial, nbodies)
    sizehint!(surface_velocity_gradient, nbodies)
    sizehint!(last_rebuild, nbodies)

    for body in systems_tuple
        body.ncells > 0 || throw(ArgumentError("PressureLaplace requires bodies with at least one panel."))
        reference_panel <= body.ncells || throw(ArgumentError(
            "reference_panel=$(reference_panel) exceeds body.ncells=$(body.ncells)."))
        calc_normals!(body)
        calc_controlpoints!(body)
        body_edges = _pressure_panel_edges(body)
        L = _assemble_pressure_laplacian(body, Int(reference_panel), body_edges)
        A = _pressure_laplacian_operator(body, Int(reference_panel), body_edges)
        b = zeros(Float64, body.ncells)
        push!(velocity_dot, zeros(Float64, size(body.velocity)))
        push!(Ls, L)
        push!(pressure_operators, A)
        push!(bs, b)
        push!(ps, zeros(Float64, body.ncells))
        push!(acceleration, zeros(Float64, 3, body.ncells))
        push!(tangential, zeros(Float64, 3, body.ncells))
        push!(edges, body_edges)
        push!(workspace, Krylov.krylov_workspace(Val(:cg), A, b))
        push!(u_inertial, zeros(Float64, size(body.velocity)))
        push!(surface_velocity_gradient, zeros(Float64, 3, 3, body.ncells))
        push!(last_rebuild, false)
    end
    preconditioner = build_pressure_preconditioner(preconditioner, Ls)

    return PressureLaplace{typeof(preconditioner)}(
        Float64(rho), unsteady, Float64(atol), Float64(rtol), Int(itmax),
        preconditioner, Int(reference_panel), Float64(reference_pressure),
        rebuild_every_step, verbose, gradient_mode, acceleration_form, velocity_dot, Ls, pressure_operators, bs, ps,
        acceleration, tangential,
        edges, workspace, u_inertial, surface_velocity_gradient, vtk_fields, last_rebuild, file)
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
    calc_controlpoints!(body)

    # Panel count changes would invalidate all preallocated arrays.
    topology_changed = body.ncells != length(m.b[i_body]) ||
                       size(m.velocity_dot[i_body]) != size(body.velocity) ||
                       size(m.u_inertial[i_body]) != size(body.velocity) ||
                       size(m.surface_velocity_gradient[i_body], 3) != body.ncells
    topology_changed && throw(ArgumentError(
        "PressureLaplace does not support panel-count changes after construction. Reconstruct the monitor for the new body sizes."))

    rebuild = m.rebuild_every_step
    m.last_rebuild[i_body] = rebuild

    if rebuild
        m.L[i_body] = _assemble_pressure_laplacian(body, m.reference_panel, m.edges[i_body])
        m.pressure_operator[i_body] =
            _pressure_laplacian_operator(body, m.reference_panel, m.edges[i_body])
        rebuild_pressure_preconditioner!(m.preconditioner, m.L[i_body], i_body)
        m.workspace[i_body] =
            Krylov.krylov_workspace(Val(:cg), m.pressure_operator[i_body], m.b[i_body])
    end

    # velocity_dot currently holds -u_old; this call turns it into (u_new - u_old)/dt.
    _pressure_velocity_dot!(m, body, i_body, dt)
    _pressure_rhs!(m, body, i_body)   # build b from material acceleration
    _pressure_solve!(m, i_body)       # solve L p = b with CG
    # Store -u_new so the next call can form the finite difference.
    _pressure_store_negative_velocity!(m, body, i_body)

    if m.verbose
        ws = m.workspace[i_body]
        println("\t\tPressureLaplace($(m.acceleration_form))[step=$(i_step+1)]: rebuild=$(rebuild), panels=$(body.ncells), CG iters=$(ws.stats.niter), solved=$(ws.stats.solved)")
    end

    return m.p[i_body]
end

function write_monitor_csv!(m::PressureLaplace, dir::AbstractString, name::AbstractString,
                            i_monitor::Int, ctx::MonitorContext, systems_tuple::Tuple,
                            i_step::Int, dt::Real; overwrite::Bool=false)
    m.file || return nothing
    t = monitor_time(ctx, i_step, dt)
    for (i_system, body) in enumerate(systems_tuple)
        i_system <= length(m.p) || continue
        path = _monitor_csv_path(dir, name, i_monitor, "pressure_laplace", i_system)
        header = "step,time,system,panels,rebuild,cg_iters,cg_solved"
        _monitor_csv_open(path, i_step, overwrite, header) do io
            ws = m.workspace[i_system]
            println(io, join((i_step, t, i_system, body.ncells,
                              _monitor_csv_bool(m.last_rebuild[i_system]),
                              ws.stats.niter,
                              _monitor_csv_bool(ws.stats.solved)), ","))
        end
    end
    return nothing
end

function _register_monitor_outputs!(ctx::MonitorContext, m::PressureLaplace, systems_tuple::Tuple)
    for i in eachindex(systems_tuple)
        monitor_register!(ctx, :P, i, m.p[i])
    end
    return nothing
end

function write_vtk_fields!(vtk, m::PressureLaplace, body, i_system::Int, i_step::Int)
    (:pressure in m.vtk_fields) || return nothing
    i_system <= length(m.p) || return nothing
    length(m.p[i_system]) == body.ncells || return nothing
    vtk["gauge pressure", VTKCellData()] = m.p[i_system]
    return nothing
end

function write_vtk_fields!(vtk, m::PressureLaplace, body, i_system::Int, i_step::Int,
                           field_names::VTKFieldNameAllocator, i_monitor::Int)
    (:pressure in m.vtk_fields) || return nothing
    i_system <= length(m.p) || return nothing
    pressure = m.p[i_system]
    length(pressure) == body.ncells || return nothing
    vtk[_vtk_monitor_field_name!(field_names, "gauge pressure", m, i_monitor), VTKCellData()] = pressure
    return nothing
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
    A = m.pressure_operator[i_body]
    if M === nothing
        Krylov.krylov_solve!(m.workspace[i_body], A, m.b[i_body]; atol=m.atol, rtol=m.rtol, itmax=m.itmax)
    else
        Krylov.krylov_solve!(m.workspace[i_body], A, m.b[i_body]; M, ldiv=true, atol=m.atol, rtol=m.rtol, itmax=m.itmax)
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

function _pressure_laplacian_operator(body::AbstractBody, reference_panel::Int)
    return _pressure_laplacian_operator(body, reference_panel, _pressure_panel_edges(body))
end

function _pressure_laplacian_operator(body::AbstractBody, reference_panel::Int,
                                      edges::Matrix{Int})
    n = body.ncells
    row_diagonal = zeros(Float64, n)
    row_neighbors = [Int[] for _ in 1:n]
    row_weights = [Float64[] for _ in 1:n]

    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        i == j && continue
        w = _pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]

        if i != reference_panel
            row_diagonal[i] += w
            if j != reference_panel
                push!(row_neighbors[i], j)
                push!(row_weights[i], w)
            end
        end
        if j != reference_panel
            row_diagonal[j] += w
            if i != reference_panel
                push!(row_neighbors[j], i)
                push!(row_weights[j], w)
            end
        end
    end

    if reference_panel > 0
        row_diagonal[reference_panel] = 1.0
        empty!(row_neighbors[reference_panel])
        empty!(row_weights[reference_panel])
    end

    return PressureLaplacianOperator(
        n, reference_panel, row_diagonal, row_neighbors, row_weights)
end

function _pressure_laplacian_mul!(y, x, α, β, A::PressureLaplacianOperator)
    length(y) == A.n || throw(DimensionMismatch(
        "PressureLaplace operator output length $(length(y)) does not match $(A.n)."))
    length(x) == A.n || throw(DimensionMismatch(
        "PressureLaplace operator input length $(length(x)) does not match $(A.n)."))

    if iszero(β)
        fill!(y, zero(eltype(y)))
    else
        y .*= β
    end

    @inbounds Threads.@threads for i in 1:A.n
        acc = A.row_diagonal[i] * x[i]
        neighbors = A.row_neighbors[i]
        weights = A.row_weights[i]
        for k in eachindex(neighbors)
            acc -= weights[k] * x[neighbors[k]]
        end
        y[i] += α * acc
    end
    return y
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
    # Refresh velocity_dot/u_inertial and diagnostic material-acceleration caches.
    _pressure_material_acceleration!(m, body, i_body)
    velocity_dot = m.unsteady ? m.velocity_dot[i_body] : nothing
    if m.acceleration_form == :material_derivative
        return _pressure_rhs_from_edge_material_derivative!(m.b[i_body], m, body, i_body,
            velocity_dot)
    elseif m.acceleration_form == :lamb_vector
        return _pressure_rhs_from_lamb_vector!(m.b[i_body], m, body, i_body,
            velocity_dot)
    else
        throw(ArgumentError("Unknown PressureLaplace acceleration_form $(m.acceleration_form)."))
    end
end

function _pressure_apply_reference_pressure_rhs!(b::AbstractVector,
                                                 reference_panel::Int,
                                                 reference_pressure::Real,
                                                 w::Real, i::Int, j::Int)
    if !iszero(reference_pressure)
        i == reference_panel && (b[j] += w * reference_pressure)
        j == reference_panel && (b[i] += w * reference_pressure)
    end
    return b
end

function _pressure_rhs_from_edge_material_derivative!(b::AbstractVector, m::PressureLaplace,
                                                      body::AbstractBody, i_body::Int,
                                                      velocity_dot::Union{Nothing,AbstractMatrix})
    edges = m.edges[i_body]
    rho = m.rho
    reference_panel = m.reference_panel
    reference_pressure = m.reference_pressure
    fill!(b, 0.0)

    # Accumulate a two-point material-derivative flux compatible with the FV
    # Laplacian. When enabled, the unsteady part uses midpoint velocity_dot; the
    # convective part approximates [(u_rel · ∇)u]·r by the edge directional
    # difference u_rel,edge · (body.velocity_j - body.velocity_i). This keeps
    # the RHS in the ∂t u + (u·∇)u form without coupling panel-center gradients
    # to a two-point pressure operator.
    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        w = _pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]
        r1 = body.controlpoints[1, j] - body.controlpoints[1, i]
        r2 = body.controlpoints[2, j] - body.controlpoints[2, i]
        r3 = body.controlpoints[3, j] - body.controlpoints[3, i]
        udot_edge_dot_r = 0.0
        if velocity_dot !== nothing
            udot_edge_dot_r = 0.5 * (
                (velocity_dot[1, i] + velocity_dot[1, j]) * r1 +
                (velocity_dot[2, i] + velocity_dot[2, j]) * r2 +
                (velocity_dot[3, i] + velocity_dot[3, j]) * r3)
        end

        nx_i, ny_i, nz_i = body.normals[1, i], body.normals[2, i], body.normals[3, i]
        nx_j, ny_j, nz_j = body.normals[1, j], body.normals[2, j], body.normals[3, j]
        ui_n = body.velocity[1, i] * nx_i + body.velocity[2, i] * ny_i + body.velocity[3, i] * nz_i
        uj_n = body.velocity[1, j] * nx_j + body.velocity[2, j] * ny_j + body.velocity[3, j] * nz_j
        urel1 = 0.5 * (body.velocity[1, i] - ui_n * nx_i + body.velocity[1, j] - uj_n * nx_j)
        urel2 = 0.5 * (body.velocity[2, i] - ui_n * ny_i + body.velocity[2, j] - uj_n * ny_j)
        urel3 = 0.5 * (body.velocity[3, i] - ui_n * nz_i + body.velocity[3, j] - uj_n * nz_j)
        du1 = body.velocity[1, j] - body.velocity[1, i]
        du2 = body.velocity[2, j] - body.velocity[2, i]
        du3 = body.velocity[3, j] - body.velocity[3, i]
        convective_edge_dot_r = urel1 * du1 + urel2 * du2 + urel3 * du3

        flux = rho * w * (udot_edge_dot_r + convective_edge_dot_r)
        b[i] += flux
        b[j] -= flux
        _pressure_apply_reference_pressure_rhs!(b, reference_panel, reference_pressure, w, i, j)
    end

    # Step 3: enforce the gauge constraint b[reference_panel] = reference_pressure.
    b[reference_panel] = reference_pressure
    return b
end

function _pressure_rhs_from_lamb_vector!(b::AbstractVector, m::PressureLaplace,
                                         body::AbstractBody, i_body::Int,
                                         velocity_dot::Union{Nothing,AbstractMatrix})
    edges = m.edges[i_body]
    rho = m.rho
    reference_panel = m.reference_panel
    reference_pressure = m.reference_pressure
    fill!(b, 0.0)

    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        w = _pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]
        r1 = body.controlpoints[1, j] - body.controlpoints[1, i]
        r2 = body.controlpoints[2, j] - body.controlpoints[2, i]
        r3 = body.controlpoints[3, j] - body.controlpoints[3, i]

        udot_edge_dot_r = 0.0
        if velocity_dot !== nothing
            udot_edge_dot_r = 0.5 * (
                (velocity_dot[1, i] + velocity_dot[1, j]) * r1 +
                (velocity_dot[2, i] + velocity_dot[2, j]) * r2 +
                (velocity_dot[3, i] + velocity_dot[3, j]) * r3)
        end

        wx_i, wy_i, wz_i = body.induced_vorticity[1, i], body.induced_vorticity[2, i], body.induced_vorticity[3, i]
        wx_j, wy_j, wz_j = body.induced_vorticity[1, j], body.induced_vorticity[2, j], body.induced_vorticity[3, j]
        wx = 0.5 * (wx_i + wx_j)
        wy = 0.5 * (wy_i + wy_j)
        wz = 0.5 * (wz_i + wz_j)

        nx_i, ny_i, nz_i = body.normals[1, i], body.normals[2, i], body.normals[3, i]
        nx_j, ny_j, nz_j = body.normals[1, j], body.normals[2, j], body.normals[3, j]
        ui_n = body.velocity[1, i] * nx_i + body.velocity[2, i] * ny_i + body.velocity[3, i] * nz_i
        uj_n = body.velocity[1, j] * nx_j + body.velocity[2, j] * ny_j + body.velocity[3, j] * nz_j
        uit1 = body.velocity[1, i] - ui_n * nx_i
        uit2 = body.velocity[2, i] - ui_n * ny_i
        uit3 = body.velocity[3, i] - ui_n * nz_i
        ujt1 = body.velocity[1, j] - uj_n * nx_j
        ujt2 = body.velocity[2, j] - uj_n * ny_j
        ujt3 = body.velocity[3, j] - uj_n * nz_j
        urel1 = 0.5 * (uit1 + ujt1)
        urel2 = 0.5 * (uit2 + ujt2)
        urel3 = 0.5 * (uit3 + ujt3)

        qi = 0.5 * (uit1^2 + uit2^2 + uit3^2)
        qj = 0.5 * (ujt1^2 + ujt2^2 + ujt3^2)
        kinetic_jump = qj - qi

        lamb1 = wy * urel3 - wz * urel2
        lamb2 = wz * urel1 - wx * urel3
        lamb3 = wx * urel2 - wy * urel1
        lamb_edge_dot_r = lamb1 * r1 + lamb2 * r2 + lamb3 * r3

        flux = rho * w * (udot_edge_dot_r + kinetic_jump + lamb_edge_dot_r)
        b[i] += flux
        b[j] -= flux
        _pressure_apply_reference_pressure_rhs!(b, reference_panel, reference_pressure, w, i, j)
    end

    b[reference_panel] = reference_pressure
    return b
end

function _pressure_rhs_from_acceleration!(b::AbstractVector, m::PressureLaplace,
                                          body::AbstractBody, i_body::Int,
                                          acceleration::AbstractMatrix)
    edges = m.edges[i_body]
    rho = m.rho
    reference_panel = m.reference_panel
    reference_pressure = m.reference_pressure
    fill!(b, 0.0)

    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        w = _pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]
        r1 = body.controlpoints[1, j] - body.controlpoints[1, i]
        r2 = body.controlpoints[2, j] - body.controlpoints[2, i]
        r3 = body.controlpoints[3, j] - body.controlpoints[3, i]
        aedge_dot_r = 0.5 * (
            (acceleration[1, i] + acceleration[1, j]) * r1 +
            (acceleration[2, i] + acceleration[2, j]) * r2 +
            (acceleration[3, i] + acceleration[3, j]) * r3)
        flux = rho * w * aedge_dot_r
        b[i] += flux
        b[j] -= flux
        _pressure_apply_reference_pressure_rhs!(b, reference_panel, reference_pressure, w, i, j)
    end

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

    G_surf = nothing
    if m.gradient_mode == :surface_velocity
        G_surf = m.surface_velocity_gradient[i_body]
        bad_panel_mask = panel_aspect_ratio_mask(body.nodes, body.cells; threshold=10.0)
        te_info = hasproperty(body, :shedding_full) ?
                  view(body.shedding_full, 1:2, :) :
                  zeros(Int, 2, body.ncells)
        compute_surface_velocity_gradient!(G_surf, u_inertial, body.controlpoints,
            body.normals, body.cells, body.neighbor, te_info;
            bad_panel_mask=any(bad_panel_mask) ? bad_panel_mask : nothing,
            nodes=body.nodes)
    end

    # In raw_hessian mode, ∇u_inertial = ∇u_induced + [Ω]_×. The induced part is
    # the analytic multipole Hessian accumulated into body.velocity_gradient by
    # FastMultipole. The kinematic part has spatial Jacobian [Ω]_× because
    # u_kinematic(x) = U_O + Ω × (x - x_O) is rigid-body. In surface_velocity
    # mode, ∇ₛu_inertial is reconstructed from the final surface velocity field,
    # including postprocessed terms such as ∇ₛμ and the kinematic velocity, so no
    # separate [Ω]_× term is added there.
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

        # In raw_hessian mode, add ([Ω]_× · u_t)_k = (Ω × u_t)_k directly to
        # avoid materializing the skew matrix.
        cross_x = ωy * ut3 - ωz * ut2
        cross_y = ωz * ut1 - ωx * ut3
        cross_z = ωx * ut2 - ωy * ut1

        if G_surf === nothing
            acceleration[1, i] += ut1 * G_ind[1, 1, i] + ut2 * G_ind[1, 2, i] + ut3 * G_ind[1, 3, i] + cross_x
            acceleration[2, i] += ut1 * G_ind[2, 1, i] + ut2 * G_ind[2, 2, i] + ut3 * G_ind[2, 3, i] + cross_y
            acceleration[3, i] += ut1 * G_ind[3, 1, i] + ut2 * G_ind[3, 2, i] + ut3 * G_ind[3, 3, i] + cross_z
        else
            acceleration[1, i] += ut1 * G_surf[1, 1, i] + ut2 * G_surf[1, 2, i] + ut3 * G_surf[1, 3, i]
            acceleration[2, i] += ut1 * G_surf[2, 1, i] + ut2 * G_surf[2, 2, i] + ut3 * G_surf[2, 3, i]
            acceleration[3, i] += ut1 * G_surf[3, 1, i] + ut2 * G_surf[3, 2, i] + ut3 * G_surf[3, 3, i]
        end
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
    NoSectionalNormalization()

Pass-through normalization for `SpanwiseLoadingMonitor`: returns dimensional
sectional loads and bin forces unchanged.
"""
struct NoSectionalNormalization end

function (::NoSectionalNormalization)(load_components, force_components, systems, frames, uinf)
    return load_components, force_components
end

"""
    FreestreamSectionalNormalization(rho, Lref)

Normalize spanwise per-length sectional loads by
`0.5*rho*|uinf|^2*Lref`. This is intended for fixed-wing steady or
quasi-steady freestream cases.
"""
struct FreestreamSectionalNormalization{TF}
    rho::TF
    Lref::TF
end

function (n::FreestreamSectionalNormalization)(load_components, force_components,
                                               systems, frames, uinf)
    qinf = n.rho / 2 * (uinf[1]^2 + uinf[2]^2 + uinf[3]^2)
    ref = qinf * n.Lref
    return load_components / ref, force_components / ref
end

"""
    RotorSectionalNormalization(rho, R, i_frame; omega_scale=:tip)

Normalize spanwise sectional loads using a rotor dynamic pressure scale from
`frames[i_frame].ω` and radius `R`. The default `omega_scale=:tip` uses
`q = 0.5*rho*(abs(ω)*R)^2`.
"""
struct RotorSectionalNormalization{TF}
    rho::TF
    R::TF
    i_frame::Int
    omega_scale::Symbol
end

function RotorSectionalNormalization(rho::Real, R::Real, i_frame::Int; omega_scale::Symbol=:tip)
    omega_scale == :tip || throw(ArgumentError(
        "RotorSectionalNormalization only supports omega_scale=:tip; got $(omega_scale)."))
    R > 0 || throw(ArgumentError("RotorSectionalNormalization requires positive R; got $(R)."))
    TF = promote_type(typeof(rho), typeof(R))
    return RotorSectionalNormalization{TF}(TF(rho), TF(R), i_frame, omega_scale)
end

function (n::RotorSectionalNormalization)(load_components, force_components,
                                          systems, frames, uinf)
    ω = frames[n.i_frame].ω
    q = n.rho / 2 * (abs(ω) * n.R)^2
    return load_components / q, force_components / q
end

function _spanwise_unit_vector(v, name::AbstractString; TF=Float64)
    length(v) == 3 || throw(ArgumentError("$(name) must have length 3."))
    sv = SVector{3, TF}(v[1], v[2], v[3])
    isapprox(norm(sv), one(TF); atol=sqrt(eps(TF)), rtol=sqrt(eps(TF))) ||
        throw(ArgumentError("$(name) must be unit length; got norm $(norm(sv))."))
    return sv
end

function _spanwise_default_axis(component_vectors::AbstractVector)
    length(component_vectors) >= 2 || throw(ArgumentError(
        "SpanwiseLoadingMonitor requires span_axis when fewer than two components are provided."))
    axis = cross(component_vectors[1], component_vectors[2])
    mag = norm(axis)
    mag > sqrt(eps(eltype(axis))) || throw(ArgumentError(
        "Cannot infer span_axis from collinear first two components."))
    return axis / mag
end

function _spanwise_validate_axis(axis, component_vectors)
    for c in component_vectors
        isapprox(dot(axis, c), zero(eltype(axis)); atol=sqrt(eps(eltype(axis))),
                 rtol=sqrt(eps(eltype(axis)))) ||
            throw(ArgumentError("span_axis must be orthogonal to every component direction."))
    end
    return axis
end

"""
    SpanwiseLoadingMonitor(nbins, i_system; i_frame=-1, components,
                           span_axis=nothing,
                           normalization=NoSectionalNormalization(),
                           per_length=false, file=true, TF=Float64,
                           verbose=false, vtk_fields=(:bin_id,))

Consume distributed panel forces from an earlier force monitor, bin them along a
span axis in the requested reference frame, and store latest-step sectional
loading. `components` is a named tuple of unit vectors expressed in the selected
frame, e.g. `components=(lift=Lhat, drag=Dhat)`.

`select` optionally restricts which panels are binned: pass a predicate
`cp_frame -> Bool` evaluated on each control point expressed in the selected
frame (`i_frame`). This is useful to isolate one blade of a multi-blade rotor,
e.g. `select = cp -> cp[2] > 0`. The default (`nothing`) bins every panel.
"""
mutable struct SpanwiseLoadingMonitor{TF, TN, TS} <: AbstractMonitor
    nbins::Int
    i_system::Int
    i_frame::Int
    component_names::Tuple{Vararg{Symbol}}
    components::Matrix{TF}
    span_axis::SVector{3, TF}
    normalization::TN
    per_length::Bool
    file::Bool
    verbose::Bool
    vtk_fields::Tuple{Vararg{Symbol}}
    bin_center::Vector{TF}
    bin_width::Vector{TF}
    load_components::Matrix{TF}
    force_components::Matrix{TF}
    counts::Vector{Int}
    panel_bin_id::Vector{Int}
    bin_force::Matrix{TF}
    select::TS                      # optional predicate cp_frame -> Bool, or nothing
end

monitor_requires(::SpanwiseLoadingMonitor) = (:F,)

function SpanwiseLoadingMonitor(nbins::Int, i_system::Int;
                                i_frame::Int=-1,
                                components,
                                span_axis=nothing,
                                normalization=NoSectionalNormalization(),
                                per_length::Bool=false,
                                file::Bool=true,
                                TF=Float64,
                                select=nothing,
                                verbose::Bool=false,
                                vtk_fields::Tuple{Vararg{Symbol}}=(:bin_id,))
    nbins > 0 || throw(ArgumentError("SpanwiseLoadingMonitor requires nbins > 0; got $(nbins)."))
    components isa NamedTuple || throw(ArgumentError("components must be a NamedTuple of unit vectors."))
    !isempty(keys(components)) || throw(ArgumentError("components must contain at least one direction."))

    names = Tuple(Symbol.(keys(components)))
    vectors = [_spanwise_unit_vector(v, "component $(name)"; TF) for (name, v) in pairs(components)]
    axis = isnothing(span_axis) ?
        _spanwise_default_axis(vectors) :
        _spanwise_unit_vector(span_axis, "span_axis"; TF)
    axis = _spanwise_validate_axis(axis, vectors)

    component_matrix = zeros(TF, 3, length(vectors))
    for (i, v) in enumerate(vectors)
        component_matrix[:, i] .= v
    end

    return SpanwiseLoadingMonitor{TF, typeof(normalization), typeof(select)}(
        nbins, i_system, i_frame, names, component_matrix, axis, normalization,
        per_length, file, verbose, vtk_fields,
        zeros(TF, nbins), zeros(TF, nbins),
        zeros(TF, length(vectors), nbins), zeros(TF, length(vectors), nbins),
        zeros(Int, nbins), Int[], zeros(TF, 3, nbins), select)
end

function _spanwise_ensure_storage!(m::SpanwiseLoadingMonitor{TF}, ncells::Int) where TF
    length(m.panel_bin_id) == ncells || (m.panel_bin_id = zeros(Int, ncells))
    size(m.bin_force) == (3, m.nbins) || (m.bin_force = zeros(TF, 3, m.nbins))
    return nothing
end

function _spanwise_frame_transform(frames, i_frame::Int)
    if i_frame < 0
        origin = SVector{3, Float64}(0.0, 0.0, 0.0)
        R_g2f = SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0,
                                          0.0, 1.0, 0.0,
                                          0.0, 0.0, 1.0)
        return origin, R_g2f
    end
    origin_global, R_f2g = frame_global_transform(frames, i_frame)
    return origin_global, transpose(R_f2g)
end

function _spanwise_fill_empty_bins!(bin_force::AbstractMatrix, counts::AbstractVector)
    populated = findall(!=(0), counts)
    isempty(populated) && throw(ArgumentError(
        "SpanwiseLoadingMonitor cannot compute loading because all bins are empty."))
    length(populated) == length(counts) && return bin_force

    first_pop = first(populated)
    for b in 1:first_pop-1
        bin_force[:, b] .= bin_force[:, first_pop]
    end
    last_pop = last(populated)
    for b in last_pop+1:length(counts)
        bin_force[:, b] .= bin_force[:, last_pop]
    end
    for k in 1:length(populated)-1
        left = populated[k]
        right = populated[k + 1]
        gap = right - left
        gap <= 1 && continue
        for b in left+1:right-1
            α = (b - left) / gap
            bin_force[:, b] .= (1 - α) .* bin_force[:, left] .+ α .* bin_force[:, right]
        end
    end
    return bin_force
end

function write_monitor_csv!(m::SpanwiseLoadingMonitor, dir::AbstractString, name::AbstractString,
                            i_monitor::Int, ctx::MonitorContext, systems_tuple::Tuple,
                            i_step::Int, dt::Real; overwrite::Bool=false)
    m.file || return nothing
    path = _monitor_csv_path(dir, name, i_monitor, "spanwise", m.i_system)
    header = join(("step", "time", "bin", "bin_center", "bin_width", "count",
                   string.(m.component_names)...), ",")
    t = monitor_time(ctx, i_step, dt)
    _monitor_csv_open(path, i_step, overwrite, header) do io
        for b in 1:m.nbins
            vals = Any[i_step, t, b, m.bin_center[b], m.bin_width[b], m.counts[b]]
            append!(vals, m.load_components[:, b])
            println(io, join(vals, ","))
        end
    end
    return nothing
end

function _run_monitor!(m::SpanwiseLoadingMonitor{TF}, ctx::MonitorContext, systems, wakes,
                       frames::AbstractVector{<:ReferenceFrame},
                       uinf, i_step::Int, dt::Real, t=nothing) where TF
    !isnothing(t) && monitor_set_time!(ctx, t)
    systems_tuple = _systems_tuple(systems)
    body = systems_tuple[m.i_system]
    F = monitor_field(ctx, :F, m.i_system)
    size(F) == (3, body.ncells) || throw(ArgumentError(
        "SpanwiseLoadingMonitor requires distributed force with size (3, $(body.ncells)); got $(size(F))."))
    body.ncells > 0 || throw(ArgumentError(
        "SpanwiseLoadingMonitor requires at least one panel."))
    _spanwise_ensure_storage!(m, body.ncells)

    origin, R_g2f = _spanwise_frame_transform(frames, m.i_frame)
    span_coord = Vector{TF}(undef, body.ncells)
    selected = Vector{Bool}(undef, body.ncells)
    @inbounds for p in 1:body.ncells
        cp_global = SVector{3, TF}(body.controlpoints[1, p], body.controlpoints[2, p], body.controlpoints[3, p])
        cp_frame = R_g2f * (cp_global - origin)
        span_coord[p] = dot(m.span_axis, cp_frame)
        selected[p] = m.select === nothing ? true : Bool(m.select(cp_frame))
    end

    any(selected) || throw(ArgumentError(
        "SpanwiseLoadingMonitor `select` excluded every panel; nothing to bin."))
    smin = minimum(span_coord[p] for p in 1:body.ncells if selected[p])
    smax = maximum(span_coord[p] for p in 1:body.ncells if selected[p])
    span = smax - smin
    width = span > sqrt(eps(TF)) ? span / m.nbins : one(TF)
    fill!(m.counts, 0)
    fill!(m.bin_force, zero(TF))
    fill!(m.panel_bin_id, 0)
    @inbounds for p in 1:body.ncells
        selected[p] || continue
        b = span > sqrt(eps(TF)) ? clamp(floor(Int, (span_coord[p] - smin) / width) + 1, 1, m.nbins) : 1
        m.panel_bin_id[p] = b
        m.counts[b] += 1
        f_global = SVector{3, TF}(F[1, p], F[2, p], F[3, p])
        f_frame = R_g2f * f_global
        m.bin_force[1, b] += f_frame[1]
        m.bin_force[2, b] += f_frame[2]
        m.bin_force[3, b] += f_frame[3]
    end

    for b in 1:m.nbins
        m.bin_width[b] = width
        m.bin_center[b] = span > sqrt(eps(TF)) ? smin + (b - 0.5) * width : smin
    end

    _spanwise_fill_empty_bins!(m.bin_force, m.counts)

    raw_force_components = zeros(TF, size(m.force_components))
    raw_load_components = zeros(TF, size(m.load_components))
    @inbounds for b in 1:m.nbins
        f = SVector{3, TF}(m.bin_force[1, b], m.bin_force[2, b], m.bin_force[3, b])
        l = m.per_length ? f / m.bin_width[b] : f
        for c in axes(m.components, 2)
            comp = SVector{3, TF}(m.components[1, c], m.components[2, c], m.components[3, c])
            raw_force_components[c, b] = dot(f, comp)
            raw_load_components[c, b] = dot(l, comp)
        end
    end

    loads, forces = m.normalization(raw_load_components, raw_force_components,
                                    systems_tuple, frames, uinf)
    m.load_components .= loads
    m.force_components .= forces

    if m.verbose
        println("\t\tSpanwiseLoadingMonitor[i_system=$(m.i_system), i_frame=$(m.i_frame), step=$(i_step+1)]:")
        println("\t\t\tspan range = ($(round(smin, sigdigits=4)), $(round(smax, sigdigits=4)))")
    end
    return nothing
end

function write_vtk_fields!(vtk, m::SpanwiseLoadingMonitor, body, i_system::Int, i_step::Int)
    (:bin_id in m.vtk_fields) || return nothing
    i_system == m.i_system || return nothing
    length(m.panel_bin_id) == body.ncells || return nothing
    vtk["spanwise bin id", VTKCellData()] = m.panel_bin_id
    return nothing
end

function write_vtk_fields!(vtk, m::SpanwiseLoadingMonitor, body, i_system::Int, i_step::Int,
                           field_names::VTKFieldNameAllocator, i_monitor::Int)
    (:bin_id in m.vtk_fields) || return nothing
    i_system == m.i_system || return nothing
    length(m.panel_bin_id) == body.ncells || return nothing
    vtk[_vtk_monitor_field_name!(field_names, "spanwise bin id", m, i_monitor), VTKCellData()] = m.panel_bin_id
    return nothing
end

"""
    ForceMonitor(nt, i_system; i_frame=-1, rho=1.0, Sref=1.0, Lref=1.0, TF=Float64, normalization=WingNormalization(...), correct_kuttacondition=true, verbose=false)

Storage monitor for integrated force and moment coefficient histories over `nt`
simulation steps. At each step the monitor reads the latest monitor-provided
pressure field, writes its own distributed-force array, then integrates force
and moment on `systems[i_system]`
in the coordinate system of `frames[i_frame]`. A `PressureBernoulli` must
therefore appear **before** this monitor in the `monitors` tuple passed to
`simulate!`.

`normalization` is called as `normalization(CF, CM, systems, frames, uinf)` and
must return a `(CF_norm, CM_norm)` tuple.  The default is `WingNormalization`
which divides by `0.5 ρ |U∞|² Sref` (and `… Lref` for moments).

If `verbose=true`, the normalized CF and CM for each step are printed to stdout
with a single `\\t` indent.
"""
mutable struct ForceMonitor{TF, TN} <: AbstractMonitor
    force::Matrix{TF}
    moment::Matrix{TF}
    distributed_force::Matrix{TF}
    i_system::Int
    i_frame::Int
    normalization::TN
    correct_kuttacondition::Bool
    verbose::Bool
    vtk_fields::Tuple{Vararg{Symbol}}
    file::Bool
end

monitor_requires(::ForceMonitor) = (:P,)
monitor_provides(::ForceMonitor) = (:F,)

"""
    SurfaceVorticityForce(body, nt, i_system; rho=1.225, i_frame=-1,
                          TF=Float64, normalization=NoNormalization(),
                          correct_kuttacondition=true,
                          grad_mu_options=(;),
                          verbose=false)

Diagnostic force monitor that reconstructs the surface vortex sheet from the
panel strength gradient and integrates the body-side Kutta-Joukowski force at
control points:
``kappa = -n x grad_s(mu)`` and ``dF = rho (V_cp x kappa) dS``.
The leading minus sign follows FLOWPanel's stored doublet/vortex-ring strength
convention, which is opposite the common exterior-minus-interior potential-jump
convention used in much of the panel-method literature.

The strength column is selected from `strength_names(body)`, preferring
`"gamma"` and falling back to `"mu"`. `body.velocity` is used exactly as stored
at the time the monitor runs. The monitor writes its own distributed-force
array, then integrates force and moment with the same frame rotation and
normalization convention as `ForceMonitor`.
"""
struct SurfaceVorticityForce{TF, TN} <: AbstractMonitor
    force::Matrix{TF}
    moment::Matrix{TF}
    distributed_force::Matrix{TF}
    grad_mu::Matrix{Float64}
    areas::Vector{Float64}
    i_system::Int
    i_frame::Int
    i_strength::Int
    rho::TF
    normalization::TN
    correct_kuttacondition::Bool
    grad_mu_options::NamedTuple
    verbose::Bool
    vtk_fields::Tuple{Vararg{Symbol}}
    file::Bool
end

monitor_provides(::SurfaceVorticityForce) = (:F,)

function SurfaceVorticityForce(body::AbstractBody, nt::Int, i_system::Int;
                               rho::Real=1.225,
                               i_frame::Int=-1,
                               TF=Float64,
                               normalization=NoNormalization(),
                               correct_kuttacondition::Bool=true,
                               grad_mu_options=(;),
                               verbose::Bool=false,
                               vtk_fields::Tuple{Vararg{Symbol}}=(:distributed_force,),
                               file::Bool=true)
    names = strength_names(body)
    i_strength = something(findfirst(==("gamma"), names),
                            findfirst(==("mu"), names),
                            0)
    i_strength == 0 && throw(ArgumentError(
        "SurfaceVorticityForce requires a body with a ConstantDoublet or "*
        "VortexRing kernel; got strengths $(names)."))

    force = zeros(TF, 3, nt)
    moment = zeros(TF, 3, nt)
    distributed_force = zeros(TF, 3, body.ncells)
    grad_mu = zeros(Float64, 3, body.ncells)
    areas = calc_areas(body)
    normalized_grad_mu_options = _normalize_grad_mu_options(grad_mu_options;
        default_basis=:quad)
    return SurfaceVorticityForce{TF, typeof(normalization)}(
        force, moment, distributed_force, grad_mu, areas, i_system, i_frame, i_strength, TF(rho),
        normalization, correct_kuttacondition, normalized_grad_mu_options,
        verbose, vtk_fields, file)
end

_surface_vorticity_te_info(body::AbstractLiftingBody) = view(body.shedding_full, 1:2, :)
_surface_vorticity_te_info(body::AbstractBody) = zeros(Int, 2, body.ncells)

function _surface_vorticity_correct_kutta!(F::AbstractMatrix, body::AbstractLiftingBody)
    @inbounds for shedding in body.shedding
        for (pi, nia, nib, pj, nja, njb) in eachcol(shedding)
            if pj != -1
                f1 = 0.5 * (F[1, pi] + F[1, pj])
                f2 = 0.5 * (F[2, pi] + F[2, pj])
                f3 = 0.5 * (F[3, pi] + F[3, pj])
                F[1, pi] = f1; F[2, pi] = f2; F[3, pi] = f3
                F[1, pj] = f1; F[2, pj] = f2; F[3, pj] = f3
            end
        end
    end
    return F
end

_surface_vorticity_correct_kutta!(F::AbstractMatrix, body::AbstractBody) = F

"""
    BoundCirculationMonitor(body, nt, i_system; i_frame, radial_dimension, R,
                            section_tol=nothing, TF=Float64, verbose=false)

Storage monitor for rotor-section bound circulation histories. Each shedding
chain in `body.shedding` is treated as one blade, and each current trailing-edge
edge midpoint is treated as one section station. The monitor stores dimensional
circulation in memory only:

- `r_over_R[section, blade]`
- `circulation_te[section, blade, step]`
- `circulation_slice[section, blade, step]`
- `valid_section[section, blade]`

The trailing-edge estimate uses the wake-strength convention
`upper_strength - lower_strength`. The slice estimate transforms panel geometry
to `frames[i_frame]`, selects panels in a radial tolerance band around the
station and on the same signed blade side, and sums signed vortex-ring/doublet
edge crossings through the section plane.
"""
struct BoundCirculationMonitor{TF} <: AbstractMonitor
    r_over_R::Matrix{TF}
    circulation_te::Array{TF, 3}
    circulation_slice::Array{TF, 3}
    valid_section::BitMatrix
    i_system::Int
    i_frame::Int
    radial_dimension::Int
    radius::TF
    section_tol::Union{Nothing, TF}
    i_strength::Int
    verbose::Bool
    file::Bool
end

function BoundCirculationMonitor(body::AbstractBody, nt::Int, i_system::Int; kwargs...)
    names = strength_names(body)
    i_strength = something(findfirst(==("gamma"), names),
                            findfirst(==("mu"), names),
                            0)
    i_strength == 0 && throw(ArgumentError(
        "BoundCirculationMonitor requires a body with a ConstantDoublet or " *
        "VortexRing kernel; got strengths $(names)."))
    throw(ArgumentError(
        "BoundCirculationMonitor requires an AbstractLiftingBody with shedding chains."))
end

function BoundCirculationMonitor(body::AbstractLiftingBody, nt::Int, i_system::Int;
                                 i_frame::Int,
                                 radial_dimension::Int,
                                 R::Real,
                                 section_tol=nothing,
                                 TF=Float64,
                                 verbose::Bool=false,
                                 file::Bool=true)
    1 <= radial_dimension <= 3 ||
        throw(ArgumentError("BoundCirculationMonitor radial_dimension must be in 1:3; got $(radial_dimension)."))
    R > 0 || throw(ArgumentError("BoundCirculationMonitor requires positive R; got $(R)."))

    names = strength_names(body)
    i_strength = something(findfirst(==("gamma"), names),
                            findfirst(==("mu"), names),
                            0)
    i_strength == 0 && throw(ArgumentError(
        "BoundCirculationMonitor requires a body with a ConstantDoublet or " *
        "VortexRing kernel; got strengths $(names)."))
    isempty(body.shedding) && throw(ArgumentError(
        "BoundCirculationMonitor requires at least one shedding chain."))

    n_blades = length(body.shedding)
    max_sections = maximum(size(s, 2) for s in body.shedding)
    max_sections > 0 || throw(ArgumentError(
        "BoundCirculationMonitor requires shedding chains with at least one section."))

    r_over_R = fill(TF(NaN), max_sections, n_blades)
    circulation_te = fill(TF(NaN), max_sections, n_blades, nt)
    circulation_slice = fill(TF(NaN), max_sections, n_blades, nt)
    valid_section = falses(max_sections, n_blades)
    for (i_blade, shedding) in pairs(body.shedding)
        valid_section[1:size(shedding, 2), i_blade] .= true
    end

    tol = isnothing(section_tol) ? nothing : TF(section_tol)
    if !isnothing(tol) && !(tol > 0)
        throw(ArgumentError("BoundCirculationMonitor section_tol must be positive; got $(section_tol)."))
    end

    return BoundCirculationMonitor{TF}(
        r_over_R, circulation_te, circulation_slice, valid_section, i_system,
        i_frame, radial_dimension, TF(R), tol, i_strength, verbose, file)
end

function _bound_circulation_frame_point(body::AbstractBody, node::Integer,
                                        origin_global, R_g2f)
    p = SVector{3}(body.nodes[1, node], body.nodes[2, node], body.nodes[3, node])
    return R_g2f * (p - origin_global)
end

function write_monitor_csv!(m::BoundCirculationMonitor, dir::AbstractString, name::AbstractString,
                            i_monitor::Int, ctx::MonitorContext, systems_tuple::Tuple,
                            i_step::Int, dt::Real; overwrite::Bool=false)
    m.file || return nothing
    path = _monitor_csv_path(dir, name, i_monitor, "bound_circulation", m.i_system)
    header = "step,time,blade,section,r_over_R,circulation_te,circulation_slice"
    t = monitor_time(ctx, i_step, dt)
    step_index = i_step + 1
    _monitor_csv_open(path, i_step, overwrite, header) do io
        for blade in axes(m.valid_section, 2), section in axes(m.valid_section, 1)
            m.valid_section[section, blade] || continue
            println(io, join((i_step, t, blade, section,
                              m.r_over_R[section, blade],
                              m.circulation_te[section, blade, step_index],
                              m.circulation_slice[section, blade, step_index]), ","))
        end
    end
    return nothing
end

function _bound_circulation_frame_controlpoint(body::AbstractBody, panel::Integer,
                                               origin_global, R_g2f)
    p = SVector{3}(body.controlpoints[1, panel],
                   body.controlpoints[2, panel],
                   body.controlpoints[3, panel])
    return R_g2f * (p - origin_global)
end

function _bound_circulation_te_midpoint(body::AbstractLiftingBody, shedding_col,
                                        origin_global, R_g2f)
    _, nia, nib, _, _, _ = shedding_col
    a = _bound_circulation_frame_point(body, nia, origin_global, R_g2f)
    b = _bound_circulation_frame_point(body, nib, origin_global, R_g2f)
    return 0.5 * (a + b)
end

function _bound_circulation_infer_tol(stations::AbstractVector{TF}) where {TF}
    vals = sort!(collect(stations))
    diffs = TF[]
    for i in 2:length(vals)
        d = abs(vals[i] - vals[i - 1])
        d > eps(TF) && push!(diffs, d)
    end
    isempty(diffs) && return one(TF)
    sort!(diffs)
    n = length(diffs)
    med = isodd(n) ? diffs[(n + 1) ÷ 2] : 0.5 * (diffs[n ÷ 2] + diffs[n ÷ 2 + 1])
    return 0.5 * med
end

function _bound_circulation_same_side(point, station, radial_dimension::Int)
    side_dim = 0
    side_mag = zero(eltype(point))
    for d in 1:3
        d == radial_dimension && continue
        mag = abs(station[d])
        if mag > side_mag
            side_mag = mag
            side_dim = d
        end
    end
    side_dim == 0 && return true
    abs(station[side_dim]) <= eps(eltype(point)) && return true
    return point[side_dim] * station[side_dim] >= zero(eltype(point))
end

function _bound_circulation_edge_crossing_sign(ra, rb, station_r)
    if (ra < station_r && station_r < rb)
        return 1
    elseif (rb < station_r && station_r < ra)
        return -1
    else
        return 0
    end
end

function _bound_circulation_slice(body::AbstractBody, station, tol, radial_dimension::Int,
                                  i_strength::Int, origin_global, R_g2f)
    station_r = station[radial_dimension]
    Γ = zero(eltype(body.strength))
    @inbounds for panel in 1:body.ncells
        cp = _bound_circulation_frame_controlpoint(body, panel, origin_global, R_g2f)
        abs(cp[radial_dimension] - station_r) <= tol || continue
        _bound_circulation_same_side(cp, station, radial_dimension) || continue

        γ = body.strength[panel, i_strength]
        nodes = body.cells[:, panel]
        for i_edge in 1:3
            na = nodes[i_edge]
            nb = nodes[i_edge == 3 ? 1 : i_edge + 1]
            pa = _bound_circulation_frame_point(body, na, origin_global, R_g2f)
            pb = _bound_circulation_frame_point(body, nb, origin_global, R_g2f)
            s = _bound_circulation_edge_crossing_sign(
                pa[radial_dimension], pb[radial_dimension], station_r)
            Γ += s * γ
        end
    end
    return Γ
end

function (m::BoundCirculationMonitor{TF})(systems, wakes,
                                           frames::AbstractVector{<:ReferenceFrame},
                                           uinf, i_step::Int, dt::Real) where {TF}
    body = systems[m.i_system]
    body isa AbstractLiftingBody || throw(ArgumentError(
        "BoundCirculationMonitor requires systems[$(m.i_system)] to be an AbstractLiftingBody."))

    origin_global, R_f2g = frame_global_transform(frames, m.i_frame)
    R_g2f = transpose(R_f2g)

    stations_by_blade = Vector{Vector{SVector{3, TF}}}(undef, length(body.shedding))
    all_station_r = TF[]
    for (i_blade, shedding) in pairs(body.shedding)
        stations = Vector{SVector{3, TF}}(undef, size(shedding, 2))
        for (i_section, col) in enumerate(eachcol(shedding))
            station = SVector{3, TF}(_bound_circulation_te_midpoint(
                body, col, origin_global, R_g2f)...)
            stations[i_section] = station
            push!(all_station_r, station[m.radial_dimension])
        end
        stations_by_blade[i_blade] = stations
    end

    tol = isnothing(m.section_tol) ?
          _bound_circulation_infer_tol(all_station_r) :
          m.section_tol

    @inbounds for (i_blade, shedding) in pairs(body.shedding)
        for (i_section, col) in enumerate(eachcol(shedding))
            station = stations_by_blade[i_blade][i_section]
            pi, _, _, pj, _, _ = col
            upper = body.strength[pi, m.i_strength]
            lower = pj != -1 ? body.strength[pj, m.i_strength] : zero(upper)

            m.r_over_R[i_section, i_blade] = station[m.radial_dimension] / m.radius
            m.circulation_te[i_section, i_blade, i_step + 1] = upper - lower
            m.circulation_slice[i_section, i_blade, i_step + 1] =
                _bound_circulation_slice(body, station, tol, m.radial_dimension,
                                         m.i_strength, origin_global, R_g2f)
        end
    end

    if m.verbose
        Γte = view(m.circulation_te, :, :, i_step + 1)
        Γslice = view(m.circulation_slice, :, :, i_step + 1)
        println("\t\tBoundCirculationMonitor[i_system=$(m.i_system), step=$(i_step+1)]:")
        println("\t\t\tΓ_TE range = ($(minimum(skipmissing(vec(Γte)))), $(maximum(skipmissing(vec(Γte)))))")
        println("\t\t\tΓ_slice range = ($(minimum(skipmissing(vec(Γslice)))), $(maximum(skipmissing(vec(Γslice)))))")
    end
end

"""
    KuttaJoukowskiForce(body, nt, i_system; rho=1.225, backend=DirectBackend(),
                        i_frame=-1, TF=Float64,
                        normalization=NoNormalization(), verbose=false)

Diagnostic force monitor that integrates the Kutta–Joukowski force on each panel
edge and reports the **force on the body**:
``F = ρ Σᵢ Σⱼ γᵢ \\, (Vᵢⱼ × Δsᵢⱼ) = -ρ Σᵢ Σⱼ γᵢ \\, (Δsᵢⱼ × Vᵢⱼ)``
summed over panels `i` and their three edges `j`, where `γᵢ` is the panel
circulation (column of `body.strength` named `"gamma"` or `"mu"`), `Δsᵢⱼ` is the
directed edge vector along the filament, and `Vᵢⱼ` is the body-relative/apparent
fluid velocity at the edge midpoint: induced velocity plus freestream, minus
the local rigid-body kinematic velocity. The induced velocity is evaluated via
a `FastMultipole.ProbeSystem`. The sign matches the pressure-integral
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
    velocity_kinematic::Matrix{TF}
    normalization::TN
    verbose::Bool
    file::Bool
end

function KuttaJoukowskiForce(body::AbstractBody, nt::Int, i_system::Int;
                              rho::Real=1.225,
                              backend::AbstractBackend=DirectBackend(),
                              i_frame::Int=-1,
                              TF=Float64,
                              normalization=NoNormalization(),
                              verbose::Bool=false,
                              file::Bool=true)

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
    velocity_kinematic = zeros(TF, 3, n_probes)

    return KuttaJoukowskiForce{TF, typeof(backend), typeof(normalization)}(
        force, i_system, i_frame, i_strength, TF(rho), backend, probes,
        edge_node_a, edge_node_b, panel_of_probe, velocity_kinematic,
        normalization, verbose, file)
end

function _kutta_joukowski_edge_kinematic_velocity!(
    velocity_kinematic::AbstractMatrix{TF},
    probes::FastMultipole.ProbeSystem{TF},
    i_system::Int,
    frames::AbstractVector{<:ReferenceFrame},
) where {TF}
    fill!(velocity_kinematic, zero(TF))
    identity = SMatrix{3, 3, TF, 9}(
        one(TF), zero(TF), zero(TF),
        zero(TF), one(TF), zero(TF),
        zero(TF), zero(TF), one(TF),
    )
    _kutta_joukowski_edge_kinematic_velocity!(
        velocity_kinematic, probes, i_system, frames, 1,
        zero(SVector{3, TF}), identity)
    return velocity_kinematic
end

function _write_force_moment_csv!(path::AbstractString, force::AbstractMatrix,
                                  moment::AbstractMatrix, i_step::Int, t::Real,
                                  overwrite::Bool)
    header = "step,time,CFx,CFy,CFz,CMx,CMy,CMz"
    col = i_step + 1
    _monitor_csv_open(path, i_step, overwrite, header) do io
        println(io, join((i_step, t,
                          force[1, col], force[2, col], force[3, col],
                          moment[1, col], moment[2, col], moment[3, col]), ","))
    end
    return nothing
end

function write_monitor_csv!(m::SurfaceVorticityForce, dir::AbstractString, name::AbstractString,
                            i_monitor::Int, ctx::MonitorContext, systems_tuple::Tuple,
                            i_step::Int, dt::Real; overwrite::Bool=false)
    m.file || return nothing
    path = _monitor_csv_path(dir, name, i_monitor, "surface_vorticity_force", m.i_system)
    return _write_force_moment_csv!(path, m.force, m.moment, i_step,
                                    monitor_time(ctx, i_step, dt), overwrite)
end

function _kutta_joukowski_edge_kinematic_velocity!(
    velocity_kinematic::AbstractMatrix{TF},
    probes::FastMultipole.ProbeSystem{TF},
    i_system::Int,
    frames::AbstractVector{<:ReferenceFrame},
    i_frame::Int,
    dx_parent_to_global,
    R_parent_to_global,
) where {TF}
    frame = frames[i_frame]

    origin_global = R_parent_to_global * frame.x + dx_parent_to_global
    v_global = R_parent_to_global * frame.v
    ω_global = R_parent_to_global * frame.ω_axis * frame.ω

    if i_system in frame.dependent_index
        @inbounds for k in eachindex(probes.position)
            dv = v_global + cross(ω_global, probes.position[k] - origin_global)
            velocity_kinematic[1, k] += dv[1]
            velocity_kinematic[2, k] += dv[2]
            velocity_kinematic[3, k] += dv[3]
        end
    end

    dx_parent_to_global = origin_global
    R_parent_to_global = R_parent_to_global * frame.R
    for i in frame.child_index
        _kutta_joukowski_edge_kinematic_velocity!(
            velocity_kinematic, probes, i_system, frames, i,
            dx_parent_to_global, R_parent_to_global)
    end
    return velocity_kinematic
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
    _kutta_joukowski_edge_kinematic_velocity!(
        m.velocity_kinematic, m.probes, m.i_system, frames)

    # 2. Collect all sources: every body, plus every wake's source list.
    wake_sources = _collect_wake_sources(wakes)
    all_sources  = (systems..., wake_sources...)

    # 3. Compute induced velocity at probes.
    old_offsets = [sys.kerneloffset for sys in systems]
    try
        for (i, sys) in pairs(systems)
            sys.kerneloffset = i == m.i_system ? sys.kerneloffset_panel : sys.kerneloffset_targets
        end
        influence!((m.probes,), all_sources, m.backend;
                    precalc=false,
                    scalar_potential=false,
                    gradient=true,
                    hessian=(false,))
    finally
        for (sys, offset) in zip(systems, old_offsets)
            sys.kerneloffset = offset
        end
    end

    # 4. Add freestream manually (no apply_freestream! method for ProbeSystem).
    uinf_sv = SVector{3, TF}(uinf[1], uinf[2], uinf[3])
    @inbounds for k in eachindex(m.probes.gradient)
        m.probes.gradient[k] += uinf_sv
        m.probes.gradient[k] -= SVector{3, TF}(
            m.velocity_kinematic[1, k],
            m.velocity_kinematic[2, k],
            m.velocity_kinematic[3, k],
        )
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

function (m::SurfaceVorticityForce{TF})(systems, wakes,
                                         frames::AbstractVector{<:ReferenceFrame},
                                         uinf, i_step::Int, dt::Real) where {TF}
    systems_tuple = _systems_tuple(systems)
    body = systems_tuple[m.i_system]
    size(m.distributed_force) == (3, body.ncells) || throw(ArgumentError(
        "SurfaceVorticityForce was constructed for $(size(m.distributed_force, 2)) panels, got $(body.ncells)."))

    F = m.distributed_force
    fill!(F, zero(eltype(F)))
    fill!(m.grad_mu, 0.0)
    calc_areas!(body, m.areas)

    compute_mu_gradient!(m.grad_mu, body.controlpoints, body.normals,
        body.cells, body.neighbor, view(body.strength, :, m.i_strength),
        _surface_vorticity_te_info(body);
        scale=1.0,
        nodes=body.nodes,
        grad_mu_options=m.grad_mu_options)

    @inbounds for i in 1:body.ncells
        n = SVector{3, Float64}(body.normals[1, i], body.normals[2, i], body.normals[3, i])
        grad = SVector{3, Float64}(m.grad_mu[1, i], m.grad_mu[2, i], m.grad_mu[3, i])
        kappa = cross(n, grad)
        V = SVector{3, Float64}(body.velocity[1, i], body.velocity[2, i], body.velocity[3, i])
        Fi = m.rho * cross(V, kappa) * m.areas[i]
        F[1, i] = Fi[1]
        F[2, i] = Fi[2]
        F[3, i] = Fi[3]
    end

    m.correct_kuttacondition && _surface_vorticity_correct_kutta!(F, body)

    Ftot = calcfield_Ftot(body, F)
    Fvec = SVector{3, TF}(Ftot[1], Ftot[2], Ftot[3])

    if m.i_frame < 0
        Mtot = calcfield_Mtot(body, zeros(3), body.controlpoints, F)
        Mvec = SVector{3, TF}(Mtot[1], Mtot[2], Mtot[3])
    else
        origin_global, R_f2g = frame_global_transform(frames, m.i_frame)
        Mtot = calcfield_Mtot(body, collect(origin_global), body.controlpoints, F)
        R_g2f = transpose(R_f2g)
        Fvec = R_g2f * Fvec
        Mvec = R_g2f * SVector{3, TF}(Mtot[1], Mtot[2], Mtot[3])
    end

    CF_norm, CM_norm = m.normalization(Fvec, Mvec, systems_tuple, frames, uinf)
    m.force[:, i_step + 1] .= CF_norm
    m.moment[:, i_step + 1] .= CM_norm

    if m.verbose
        println("\t\tSurfaceVorticityForce[i_system=$(m.i_system), i_frame=$(m.i_frame), step=$(i_step+1)]:")
        println("\t\t\tCF = ($(round(CF_norm[1], sigdigits=4)), $(round(CF_norm[2], sigdigits=4)), $(round(CF_norm[3], sigdigits=4)))")
        println("\t\t\tCM = ($(round(CM_norm[1], sigdigits=4)), $(round(CM_norm[2], sigdigits=4)), $(round(CM_norm[3], sigdigits=4)))")
    end
end

function _register_monitor_outputs!(ctx::MonitorContext, m::SurfaceVorticityForce, systems_tuple::Tuple)
    monitor_register!(ctx, :F, m.i_system, m.distributed_force)
    return nothing
end

function write_vtk_fields!(vtk, m::SurfaceVorticityForce, body, i_system::Int, i_step::Int)
    (:distributed_force in m.vtk_fields) || return nothing
    i_system == m.i_system || return nothing
    size(m.distributed_force) == (3, body.ncells) || return nothing
    vtk["F", VTKCellData()] = m.distributed_force
    return nothing
end

function write_vtk_fields!(vtk, m::SurfaceVorticityForce, body, i_system::Int, i_step::Int,
                           field_names::VTKFieldNameAllocator, i_monitor::Int)
    (:distributed_force in m.vtk_fields) || return nothing
    i_system == m.i_system || return nothing
    size(m.distributed_force) == (3, body.ncells) || return nothing
    vtk[_vtk_monitor_field_name!(field_names, "F", m, i_monitor), VTKCellData()] = m.distributed_force
    return nothing
end


function ForceMonitor(nt::Int, i_system::Int;
                       i_frame=-1, rho=1.0, Sref=1.0, Lref=1.0, TF=Float64,
                       normalization=WingNormalization(TF(rho), TF(Sref), TF(Lref)),
                       correct_kuttacondition::Bool=true,
                       verbose::Bool=false,
                       vtk_fields::Tuple{Vararg{Symbol}}=(:distributed_force,),
                       file::Bool=true)
    force = zeros(TF, 3, nt)
    moment = zeros(TF, 3, nt)
    distributed_force = zeros(TF, 3, 0)
    return ForceMonitor{TF, typeof(normalization)}(
        force, moment, distributed_force, i_system, i_frame, normalization,
        correct_kuttacondition, verbose, vtk_fields, file)
end

function _run_monitor!(monitor::ForceMonitor, ctx::MonitorContext, systems, wakes,
                       frames::AbstractVector{<:ReferenceFrame},
                       uinf, i_step::Int, dt::Real, t=nothing)
    !isnothing(t) && monitor_set_time!(ctx, t)
    systems_tuple = _systems_tuple(systems)
    body = systems_tuple[monitor.i_system]
    pressure = monitor_field(ctx, :P, monitor.i_system)

    if size(monitor.distributed_force) != (3, body.ncells)
        monitor.distributed_force = zeros(eltype(monitor.force), 3, body.ncells)
    end
    fill!(monitor.distributed_force, zero(eltype(monitor.distributed_force)))
    calcfield_F!(monitor.distributed_force, body, calc_areas(body), body.normals, pressure;
                 correct_kuttacondition=monitor.correct_kuttacondition)

    # total force in global frame
    Ftot = calcfield_Ftot(body, monitor.distributed_force)
    Fvec = FastMultipole.SVector{3}(Ftot[1], Ftot[2], Ftot[3])

    if monitor.i_frame < 0
        # global frame: moment about the origin
        Mtot = calcfield_Mtot(body, zeros(3), body.controlpoints, monitor.distributed_force)
        Mvec = FastMultipole.SVector{3}(Mtot[1], Mtot[2], Mtot[3])
    else
        # frame-local: moment about frame origin, rotated into frame axes
        origin_global, R_f2g = frame_global_transform(frames, monitor.i_frame)
        Mtot = calcfield_Mtot(body, collect(origin_global), body.controlpoints, monitor.distributed_force)
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
    monitor_register!(ctx, :F, monitor.i_system, monitor.distributed_force)
    return nothing
end

function write_monitor_csv!(m::ForceMonitor, dir::AbstractString, name::AbstractString,
                            i_monitor::Int, ctx::MonitorContext, systems_tuple::Tuple,
                            i_step::Int, dt::Real; overwrite::Bool=false)
    m.file || return nothing
    path = _monitor_csv_path(dir, name, i_monitor, "force", m.i_system)
    return _write_force_moment_csv!(path, m.force, m.moment, i_step,
                                    monitor_time(ctx, i_step, dt), overwrite)
end

function write_monitor_csv!(m::KuttaJoukowskiForce, dir::AbstractString, name::AbstractString,
                            i_monitor::Int, ctx::MonitorContext, systems_tuple::Tuple,
                            i_step::Int, dt::Real; overwrite::Bool=false)
    m.file || return nothing
    path = _monitor_csv_path(dir, name, i_monitor, "kutta_joukowski_force", m.i_system)
    header = "step,time,CFx,CFy,CFz"
    col = i_step + 1
    _monitor_csv_open(path, i_step, overwrite, header) do io
        println(io, join((i_step, monitor_time(ctx, i_step, dt),
                          m.force[1, col], m.force[2, col], m.force[3, col]), ","))
    end
    return nothing
end

function (monitor::ForceMonitor)(systems, wakes,
                                  frames::AbstractVector{<:ReferenceFrame},
                                  uinf, i_step::Int, dt::Real)
    throw(ArgumentError("ForceMonitor requires monitor-owned pressure context; run it through simulate! or replay after a pressure monitor."))
end

function write_vtk_fields!(vtk, m::ForceMonitor, body, i_system::Int, i_step::Int)
    (:distributed_force in m.vtk_fields) || return nothing
    i_system == m.i_system || return nothing
    size(m.distributed_force) == (3, body.ncells) || return nothing
    vtk["F", VTKCellData()] = m.distributed_force
    return nothing
end

function write_vtk_fields!(vtk, m::ForceMonitor, body, i_system::Int, i_step::Int,
                           field_names::VTKFieldNameAllocator, i_monitor::Int)
    (:distributed_force in m.vtk_fields) || return nothing
    i_system == m.i_system || return nothing
    size(m.distributed_force) == (3, body.ncells) || return nothing
    vtk[_vtk_monitor_field_name!(field_names, "F", m, i_monitor), VTKCellData()] = m.distributed_force
    return nothing
end
