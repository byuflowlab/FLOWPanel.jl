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
    PressureBernoulli(rho; unsteady=false, correct_kuttacondition=false, clip=nothing,
                      allow_partial=false,
                      backend=FastMultipoleBackend(expansion_order=10,
                                                   multipole_acceptance=0.4,
                                                   leaf_size=100))

Monitor that owns pressure arrays for every body in the simulation by
evaluating the Bernoulli equation each step.

With `unsteady=false` (default) the monitor evaluates the *body-relative steady
loading formulation* ``P = \\tfrac{1}{2} \\rho (U_\\infty^2 - |u_{rel,t}|^2)``,
where ``u_{rel,t}`` is the tangential projection of the body-relative surface
velocity (`body.velocity`, i.e. with the rigid-body kinematic velocity already
subtracted). This is valid for flows that are steady in the body frame (e.g. a
rotor at constant rotation rate). In a rotating frame the complete steady
pressure relation also carries a centrifugal/reference-potential term
``\\tfrac{1}{2}\\rho|\\mathbf{w}|^2`` which is omitted here, so the steady
pressure is defined only up to that rotating-frame reference contribution; the
omitted term is symmetric across a blade section and loading-neutral, but the
steady field is not the complete absolute pressure for every rotating body.

With `unsteady=true` the term ``-\\rho \\, \\partial \\phi / \\partial t`` is
added and the kinetic term uses the reconstructed total *inertial* surface
velocity instead. The monitor evaluates the exterior scalar-potential trace at
the body control points, finite differences it with zero/BE/variable-step-BDF2
startup, and uses the Arbitrary Lagrangian--Eulerian (ALE) identity
``\\partial_t\\phi=D_g\\phi-\\mathbf{w}\\cdot\\nabla\\phi``.

Vector-potential-capable body sources are unsupported and throw in unsteady
mode. Vector-potential-only wake sources also throw by default because their
scalar-potential contribution is unavailable. `allow_partial=true` explicitly
enables a mixed partial diagnostic: their velocity remains in kinetic energy
but is excluded from ``\\mathbf{w}\\cdot\\nabla\\phi``. This mode warns once per
monitor.

`correct_kuttacondition=true` opts into heuristic averaging of pressures on
trailing-edge panel pairs. It and `clip` are forwarded to `calcfield_P!`.

Pressure-dependent monitors (e.g. `ForceMonitor`) must appear *after* this
monitor in the `monitors` tuple passed to `simulate!`.
"""
mutable struct PressureBernoulli{UN, TF, TC, TB} <: AbstractMonitor
    rho::TF
    unsteady::Bool
    correct_kuttacondition::Bool
    allow_partial::Bool
    clip::TC
    backend::TB
    pressure::Vector{Vector{Float64}}
    phi_dot::Vector{Vector{Float64}}
    potential_history::Vector{Vector{Float64}}
    potential_history_older::Vector{Vector{Float64}}
    probes::Vector{FastMultipole.ProbeSystem{Float64}}
    inertial_velocity::Vector{Matrix{Float64}}
    history_count::Vector{Int}
    previous_dt::Vector{Float64}
    older_dt::Vector{Float64}
    last_step::Vector{Int}
    warned_excluded_sources::Bool
    vtk_fields::Tuple{Vararg{Symbol}}
    file::Bool
end

monitor_provides(::PressureBernoulli) = (:P,)

function PressureBernoulli(rho::Real; unsteady::Bool=false,
                          correct_kuttacondition::Bool=false,
                          allow_partial::Bool=false,
                          clip=nothing,
                          backend=FastMultipoleBackend(expansion_order=10,
                                                       multipole_acceptance=0.4,
                                                       leaf_size=100),
                          vtk_fields::Tuple{Vararg{Symbol}}=(:pressure,),
                          file::Bool=true)
    return PressureBernoulli{unsteady, typeof(rho), typeof(clip), typeof(backend)}(
        rho, unsteady, correct_kuttacondition, allow_partial, clip, backend,
        Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[],
        FastMultipole.ProbeSystem{Float64}[], Matrix{Float64}[], Int[], Float64[],
        Float64[], Int[], false, vtk_fields, file)
end

function (m::PressureBernoulli)(systems, wakes,
                               frames::AbstractVector{<:ReferenceFrame},
                               uinf, i_step::Int, dt::Real)
    systems_tuple = _systems_tuple(systems)
    Uinf_mag = norm(uinf)
    scalar_sources, excluded_sources = m.unsteady ?
        _bernoulli_source_partition(m, systems_tuple, wakes) : ((), ())
    if m.unsteady && !isempty(excluded_sources) && !m.warned_excluded_sources
        @warn "PressureBernoulli(unsteady=true, allow_partial=true) excludes vector-potential-only wake sources from scalar-potential history and the scalar ALE contraction. Continuing with a partial diagnostic: their velocity is retained only in inertial kinetic energy."
        m.warned_excluded_sources = true
    end
    _pressure_bernoulli_ensure_storage!(m, systems_tuple)
    for (i_body, body) in enumerate(systems_tuple)
        pressure = m.pressure[i_body]
        fill!(pressure, 0.0)
        phi_dot = m.unsteady ?
            _pressure_bernoulli_phi_dot!(m, body, i_body, scalar_sources,
                excluded_sources, uinf, i_step, dt) :
            nothing
        # Steady mode evaluates the body-relative steady loading formulation,
        # so the kinetic energy uses the relative surface trace; unsteady mode
        # uses the inertial trace, compensated by the ALE phi_dot term. Both
        # project out the normal component so that finite Dirichlet normal
        # leakage is not counted as physical kinetic energy.
        pressure_velocity = m.unsteady ?
            _pressure_fill_inertial_surface_velocity!(m.inertial_velocity[i_body], body) :
            _pressure_fill_relative_surface_velocity!(m.inertial_velocity[i_body], body)
        calcfield_P!(pressure, body, pressure_velocity, Uinf_mag, m.rho, phi_dot;
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
    if length(m.pressure) != nbodies || length(m.inertial_velocity) != nbodies ||
       (m.unsteady && (length(m.phi_dot) != nbodies ||
                       length(m.potential_history) != nbodies ||
                       length(m.potential_history_older) != nbodies ||
                       length(m.probes) != nbodies))
        m.pressure = Vector{Float64}[]
        m.phi_dot = Vector{Float64}[]
        m.potential_history = Vector{Float64}[]
        m.potential_history_older = Vector{Float64}[]
        m.probes = FastMultipole.ProbeSystem{Float64}[]
        m.inertial_velocity = Matrix{Float64}[]
        empty!(m.history_count); empty!(m.previous_dt); empty!(m.older_dt); empty!(m.last_step)
        sizehint!(m.pressure, nbodies)
        sizehint!(m.phi_dot, nbodies)
        sizehint!(m.potential_history, nbodies); sizehint!(m.potential_history_older, nbodies)
        sizehint!(m.probes, nbodies); sizehint!(m.inertial_velocity, nbodies)
        for body in systems_tuple
            push!(m.pressure, zeros(Float64, body.ncells))
            push!(m.inertial_velocity, zeros(Float64, size(body.velocity)))
            if m.unsteady
                push!(m.phi_dot, zeros(Float64, body.ncells))
                push!(m.potential_history, zeros(Float64, body.ncells))
                push!(m.potential_history_older, zeros(Float64, body.ncells))
                push!(m.probes, FastMultipole.ProbeSystem(body.ncells, Float64))
                push!(m.history_count, 0); push!(m.previous_dt, NaN)
                push!(m.older_dt, NaN); push!(m.last_step, typemin(Int))
            end
        end
        return nothing
    end

    for (i, body) in enumerate(systems_tuple)
        if length(m.pressure[i]) != body.ncells ||
           size(m.inertial_velocity[i]) != size(body.velocity) ||
           (m.unsteady && (length(m.phi_dot[i]) != body.ncells ||
                           length(m.potential_history[i]) != body.ncells ||
                           length(m.potential_history_older[i]) != body.ncells ||
                           length(m.probes[i].scalar_potential) != body.ncells))
            m.pressure[i] = zeros(Float64, body.ncells)
            m.inertial_velocity[i] = zeros(Float64, size(body.velocity))
            if m.unsteady
                m.phi_dot[i] = zeros(Float64, body.ncells)
                m.potential_history[i] = zeros(Float64, body.ncells)
                m.potential_history_older[i] = zeros(Float64, body.ncells)
                m.probes[i] = FastMultipole.ProbeSystem(body.ncells, Float64)
                m.history_count[i] = 0; m.previous_dt[i] = NaN
                m.older_dt[i] = NaN; m.last_step[i] = typemin(Int)
            end
        end
    end

    return nothing
end

function _bernoulli_source_partition(m::PressureBernoulli, systems_tuple::Tuple, wakes)
    vector_bodies = Tuple(body for body in systems_tuple
        if FastMultipole.has_vector_potential(body))
    isempty(vector_bodies) || throw(ArgumentError(
        "PressureBernoulli(unsteady=true) does not support vector-potential-capable " *
        "body sources because their retained exterior scalar trace cannot be separated."))

    wake_sources = _collect_wake_sources(wakes)
    excluded_sources = Tuple(source for source in wake_sources
        if FastMultipole.has_vector_potential(source))
    if !isempty(excluded_sources) && !m.allow_partial
        throw(ArgumentError(
            "PressureBernoulli(unsteady=true) encountered vector-potential-only wake " *
            "sources. Pass allow_partial=true to request the warned partial diagnostic."))
    end
    all_sources = (systems_tuple..., wake_sources...)
    return _filter_scalar_potential_sources(all_sources), excluded_sources
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
                                      excluded_sources::Tuple,
                                      uinf, i_step::Int, dt::Real)
    dt > 0 || throw(ArgumentError("PressureBernoulli(unsteady=true) requires a positive runtime dt; got $(dt)."))

    probes = m.probes[i_body]
    history = m.potential_history[i_body]
    phi_dot = m.phi_dot[i_body]
    older = m.potential_history_older[i_body]

    @inbounds for p in 1:body.ncells
        probes.position[p] = SVector{3, Float64}(
            body.controlpoints[1, p],
            body.controlpoints[2, p],
            body.controlpoints[3, p],
        )
        probes.scalar_potential[p] = 0.0
        probes.gradient[p] = zero(eltype(probes.gradient))
    end

    if length(scalar_sources) > 0
        influence!((probes,), scalar_sources, m.backend;
            precalc=false, scalar_potential=true, gradient=false,
            hessian=(false,))
    end

    # Vector-potential-only wake velocity belongs to total kinetic energy but
    # not to the gradient of the retained scalar potential in the ALE term.
    if length(excluded_sources) > 0
        influence!((probes,), excluded_sources, m.backend;
            precalc=false, scalar_potential=false, gradient=true,
            hessian=(false,))
    end

    @inbounds for p in 1:body.ncells
        phi = probes.scalar_potential[p] +
              uinf[1] * body.controlpoints[1, p] +
              uinf[2] * body.controlpoints[2, p] +
              uinf[3] * body.controlpoints[3, p]
        # Exact-control-point body kernels return FLOWPanel's canonical
        # interior limit. Bernoulli history requires the exterior trace.
        has_grad_mu(body) && (phi -= body.strength[p, get_Gammai(body)])
        count = m.history_count[i_body]
        consecutive = i_step == m.last_step[i_body] + 1
        valid_history = consecutive && count > 0
        Dphi_Dt = if !valid_history
            0.0
        elseif count == 1
            (phi - history[p]) / m.previous_dt[i_body]
        else
            _variable_bdf2(phi, history[p], older[p], m.previous_dt[i_body], m.older_dt[i_body])
        end
        vk1 = body.velocity_kinematic[1, p]
        vk2 = body.velocity_kinematic[2, p]
        vk3 = body.velocity_kinematic[3, p]
        excluded = probes.gradient[p]
        grad_phi_1 = body.velocity[1, p] + vk1 - excluded[1]
        grad_phi_2 = body.velocity[2, p] + vk2 - excluded[2]
        grad_phi_3 = body.velocity[3, p] + vk3 - excluded[3]
        phi_dot[p] = valid_history ?
            Dphi_Dt - (vk1 * grad_phi_1 + vk2 * grad_phi_2 + vk3 * grad_phi_3) : 0.0
        older[p] = history[p]
        history[p] = phi
    end

    if i_step != m.last_step[i_body] + 1
        m.history_count[i_body] = 1
    else
        m.history_count[i_body] = min(m.history_count[i_body] + 1, 2)
    end
    m.older_dt[i_body] = m.previous_dt[i_body]
    m.previous_dt[i_body] = Float64(dt)
    m.last_step[i_body] = i_step

    return phi_dot
end

@inline function _variable_bdf2(current, previous, older, h::Real, hprev::Real)
    (isfinite(h) && h > 0 && isfinite(hprev) && hprev > 0) ||
        throw(ArgumentError("Variable-step BDF2 requires positive finite intervals; got h=$(h), hprev=$(hprev)."))
    return ((2h + hprev) / (h * (h + hprev))) * current -
           ((h + hprev) / (h * hprev)) * previous +
           (h / (hprev * (h + hprev))) * older
end

function _pressure_fill_relative_surface_velocity!(out::AbstractMatrix, body::AbstractBody)
    @inbounds for p in 1:body.ncells
        nx, ny, nz = body.normals[1,p], body.normals[2,p], body.normals[3,p]
        qx, qy, qz = body.velocity[1,p], body.velocity[2,p], body.velocity[3,p]
        qn = qx*nx + qy*ny + qz*nz
        out[1,p] = qx - qn*nx
        out[2,p] = qy - qn*ny
        out[3,p] = qz - qn*nz
    end
    return out
end

function _pressure_fill_inertial_surface_velocity!(out::AbstractMatrix, body::AbstractBody)
    @inbounds for p in 1:body.ncells
        nx, ny, nz = body.normals[1,p], body.normals[2,p], body.normals[3,p]
        qx, qy, qz = body.velocity[1,p], body.velocity[2,p], body.velocity[3,p]
        qn = qx*nx + qy*ny + qz*nz
        out[1,p] = qx - qn*nx + body.velocity_kinematic[1,p]
        out[2,p] = qy - qn*ny + body.velocity_kinematic[2,p]
        out[3,p] = qz - qn*nz + body.velocity_kinematic[3,p]
    end
    return out
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
                    rebuild_every_step=false, verbose=false, gradient_mode=:unvalidated,
                    acceleration_form=:material_derivative)

Monitor that owns pressure by solving a sparse panel-centered surface
pressure Poisson equation. With `unsteady=true`, it differences the inertial
surface trace `u_inertial = tangent(body.velocity) + body.velocity_kinematic`.
The first sample seeds history with zero derivative, the next uses backward
Euler, and later samples use variable-step BDF2. Steady configurations allocate
no temporal-history arrays.

No formulation has yet passed the complete unsteady pitching-wing acceptance
gate. Omitting `gradient_mode` therefore emits a warning and uses corrected
`:edge_difference` only as a compatibility fallback; production validation
should select a mode explicitly.

The diagnostic `gradient_mode=:corrected_hessian` uses the inertial kernel
Jacobian in `body.velocity_gradient` and adds the surface derivative of the
postprocessed exterior `-½∇ₛμ_code` jump. It is unsupported for production
because the exterior hypersingular surface limit is not implemented. No
`[Ω]×` term is added: grid kinematics modify velocity values, not the stored
kernel gradient. The complete ALE acceleration is projected tangentially once
before conservative edge divergence. `:surface_velocity` reconstructs the
full trace gradient with node-aware paired quads when available. The default
`:edge_difference` is Hessian-free and uses exterior volumetric
vorticity (with bound-sheet κ removed) to correct nonsymmetric gradients.
`:raw_hessian` is accepted as a compatibility alias
for `:corrected_hessian`.

The pressure RHS has two acceleration forms. The default
`acceleration_form=:material_derivative` uses the selected corrected gradient
mode. The deprecated `acceleration_form=:lamb_vector` uses the Lamb
decomposition `(u · ∇)u = ∇(|u|²/2) + ω × u` with the same tangent relative
velocity used by the material-derivative edge form, where `ω` is the volumetric
induced vorticity accumulated in `body.induced_vorticity`. Lamb mode remains a
diagnostic because its complete ALE surface identity is not established.

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
- `gradient_mode::Symbol` — explicit `:edge_difference`, `:corrected_hessian`, or `:surface_velocity`; omission uses a warned, unvalidated edge fallback
- `acceleration_form::Symbol` — `:material_derivative` for the direct acceleration edge form or `:lamb_vector` for the Lamb-vector decomposition
- `lamb_vorticity::Symbol` — vorticity fed to the `:lamb_vector` RHS: `:induced` (default; `body.induced_vorticity` as accumulated by the solver = wake-induced curl + bound-sheet κ), `:no_bound` (induced − κ), `:bound_only` (κ alone), or `:hessian_curl` (curl of the FastMultipole velocity Hessian + κ)
- `kappa_basis::Symbol` — μ-gradient basis used when the monitor itself computes κ for the `lamb_vorticity` variants (`:quad` default, matching what `simulate!`/`replay` accumulate into `body.induced_vorticity`; `:tri` available for sensitivity studies — the rotor-hover lamb CT was observed to move by ~2.5e-3 between the two)
- `project_edge_du::Bool` — deprecated compatibility field; corrected edge mode never blanket-projects `Δu`
- `velocity_dot::Vector{Matrix{Float64}}` — current zero/BE/BDF2 panel-following derivative; empty for steady monitors
- `u_inertial::Vector{Matrix{Float64}}` — scratch buffer for the inertial fluid velocity `tangent(body.velocity) + body.velocity_kinematic` (3 × ncells per body)
- `surface_velocity_gradient::Vector{Array{Float64,3}}` — scratch buffer for `∇ₛu_inertial` in `:surface_velocity` mode (3 × 3 × ncells per body)
- `L::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}}` — sparse FV surface Laplacian per body, gauge-fixed at `reference_panel`
- `pressure_operator::Vector{PressureLaplacianOperator}` — matrix-free FV surface Laplacian per body for CG matvecs
- `b::Vector{Vector{Float64}}` — RHS vector per body (length ncells); rebuilt each call from the selected acceleration form
- `p::Vector{Vector{Float64}}` — owned pressure solution vector per body (length ncells)
- `residual::Vector{Vector{Float64}}` — reusable `L*p-b` residual workspace per body
- `acceleration::Vector{Matrix{Float64}}` — material acceleration `Du/Dt` scratch buffer (3 × ncells per body)
- `tangential::Vector{Matrix{Float64}}` — tangential projection of `acceleration` (3 × ncells per body)
- `edges::Vector{Matrix{Int}}` — shared interior edges per body (4 × nedges); rows are `(node_a, node_b, panel_i, panel_j)`
- `workspace::Vector{Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}}` — preallocated Krylov CG workspace per body
- `absolute_residual`, `relative_residual` — independently evaluated `L*p-b` solve diagnostics per body
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

mutable struct PressureLaplace{UN, GM, AF, LV, TP} <: AbstractMonitor
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
    lamb_vorticity::Symbol
    kappa_basis::Symbol
    project_edge_du::Bool
    velocity_dot::Vector{Matrix{Float64}}
    velocity_history::Vector{Matrix{Float64}}
    velocity_history_older::Vector{Matrix{Float64}}
    history_count::Vector{Int}
    previous_dt::Vector{Float64}
    older_dt::Vector{Float64}
    last_step::Vector{Int}
    L::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}}
    pressure_operator::Vector{PressureLaplacianOperator}
    b::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}
    residual::Vector{Vector{Float64}}
    acceleration::Vector{Matrix{Float64}}
    tangential::Vector{Matrix{Float64}}
    edges::Vector{Matrix{Int}}
    workspace::Vector{Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}}
    u_inertial::Vector{Matrix{Float64}}
    omega_used::Vector{Matrix{Float64}}
    surface_velocity_gradient::Vector{Array{Float64,3}}
    jump_velocity::Vector{Matrix{Float64}}
    jump_velocity_gradient::Vector{Array{Float64,3}}
    vtk_fields::Tuple{Vararg{Symbol}}
    last_rebuild::Vector{Bool}
    absolute_residual::Vector{Float64}
    relative_residual::Vector{Float64}
    convergence_warned::Bool
    file::Bool
end

monitor_provides(::PressureLaplace) = (:P,)
monitor_requires_body_hessian(::PressureLaplace{UN,GM,AF,LV}) where {UN,GM,AF,LV} =
    (GM == :corrected_hessian && AF == :material_derivative) ||
    (AF == :lamb_vector && LV == :hessian_curl)
monitor_requires_induced_vorticity(::PressureLaplace{UN,GM,AF}) where {UN,GM,AF} =
    AF == :lamb_vector || GM == :edge_difference

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
                         gradient_mode::Symbol=:unvalidated,
                         acceleration_form::Symbol=:material_derivative,
                         lamb_vorticity::Symbol=:induced,
                         kappa_basis::Symbol=:quad,
                         project_edge_du::Bool=false,
                         vtk_fields::Tuple{Vararg{Symbol}}=(:pressure,),
                         file::Bool=true)
    reference_panel >= 1 || throw(ArgumentError("reference_panel must be at least 1; got $(reference_panel)."))
    if gradient_mode == :unvalidated
        @warn "PressureLaplace has no formulation that has passed all manufactured, steady, direct/FMM, and unsteady pitching-wing gates. Using :edge_difference as an unvalidated compatibility fallback; select gradient_mode explicitly for diagnostic or validation work." maxlog=1
        gradient_mode = :edge_difference
    end
    gradient_mode in (:corrected_hessian, :raw_hessian, :surface_velocity, :edge_difference) || throw(ArgumentError(
        "gradient_mode must be :corrected_hessian, :surface_velocity, or :edge_difference; got $(gradient_mode)."))
    gradient_mode == :raw_hessian && (gradient_mode = :corrected_hessian)
    if gradient_mode == :corrected_hessian
        @warn "PressureLaplace gradient_mode=:corrected_hessian is an unsupported diagnostic: the exterior hypersingular surface limit of the panel Hessian is not implemented. On-sheet self Hessians and edge-singular neighboring derivatives do not define a convergent exterior surface gradient; use :surface_velocity for gradient-based pressure recovery." maxlog=1
    end
    acceleration_form in (:material_derivative, :lamb_vector) || throw(ArgumentError(
        "acceleration_form must be :material_derivative or :lamb_vector; got $(acceleration_form)."))
    if acceleration_form == :lamb_vector
        @warn "PressureLaplace(acceleration_form=:lamb_vector) is deprecated and retained only as a diagnostic for one compatibility cycle. It lacks a complete ALE surface derivation; use acceleration_form=:material_derivative with :corrected_hessian, :edge_difference, or :surface_velocity." maxlog=1
    end
    lamb_vorticity in (:induced, :no_bound, :bound_only, :hessian_curl) || throw(ArgumentError(
        "lamb_vorticity must be :induced, :no_bound, :bound_only, or :hessian_curl; got $(lamb_vorticity)."))
    kappa_basis in (:quad, :tri) || throw(ArgumentError(
        "kappa_basis must be :quad or :tri; got $(kappa_basis)."))
    systems_tuple = _systems_tuple(bodies)
    isempty(systems_tuple) && throw(ArgumentError("PressureLaplace requires at least one body."))
    nbodies = length(systems_tuple)
    velocity_dot = Matrix{Float64}[]
    velocity_history = Matrix{Float64}[]
    velocity_history_older = Matrix{Float64}[]
    history_count = Int[]; previous_dt = Float64[]; older_dt = Float64[]; last_step = Int[]
    Ls = SparseArrays.SparseMatrixCSC{Float64, Int}[]
    pressure_operators = PressureLaplacianOperator[]
    bs = Vector{Float64}[]
    ps = Vector{Float64}[]
    residuals = Vector{Float64}[]
    acceleration = Matrix{Float64}[]
    tangential = Matrix{Float64}[]
    edges = Matrix{Int}[]
    workspace = Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}[]
    u_inertial = Matrix{Float64}[]
    omega_used = Matrix{Float64}[]
    surface_velocity_gradient = Array{Float64,3}[]
    jump_velocity = Matrix{Float64}[]
    jump_velocity_gradient = Array{Float64,3}[]
    last_rebuild = Bool[]
    absolute_residual = Float64[]
    relative_residual = Float64[]
    sizehint!(velocity_dot, nbodies)
    sizehint!(velocity_history, nbodies); sizehint!(velocity_history_older, nbodies)
    sizehint!(Ls, nbodies)
    sizehint!(pressure_operators, nbodies)
    sizehint!(bs, nbodies)
    sizehint!(ps, nbodies)
    sizehint!(residuals, nbodies)
    sizehint!(acceleration, nbodies)
    sizehint!(tangential, nbodies)
    sizehint!(edges, nbodies)
    sizehint!(workspace, nbodies)
    sizehint!(u_inertial, nbodies)
    sizehint!(omega_used, nbodies)
    sizehint!(surface_velocity_gradient, nbodies)
    sizehint!(jump_velocity, nbodies); sizehint!(jump_velocity_gradient, nbodies)
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
        if unsteady
            push!(velocity_dot, zeros(Float64, size(body.velocity)))
            push!(velocity_history, zeros(Float64, size(body.velocity)))
            push!(velocity_history_older, zeros(Float64, size(body.velocity)))
            push!(history_count, 0); push!(previous_dt, NaN); push!(older_dt, NaN)
            push!(last_step, typemin(Int))
        end
        push!(Ls, L)
        push!(pressure_operators, A)
        push!(bs, b)
        push!(ps, zeros(Float64, body.ncells))
        push!(residuals, zeros(Float64, body.ncells))
        push!(acceleration, zeros(Float64, 3, body.ncells))
        push!(tangential, zeros(Float64, 3, body.ncells))
        push!(edges, body_edges)
        push!(workspace, Krylov.krylov_workspace(Val(:cg), A, b))
        push!(u_inertial, zeros(Float64, size(body.velocity)))
        push!(omega_used, zeros(Float64, 3, body.ncells))
        push!(surface_velocity_gradient, zeros(Float64, 3, 3, body.ncells))
        push!(jump_velocity, zeros(Float64, 3, body.ncells))
        push!(jump_velocity_gradient, zeros(Float64, 3, 3, body.ncells))
        push!(last_rebuild, false)
        push!(absolute_residual, Inf)
        push!(relative_residual, Inf)
    end
    preconditioner = build_pressure_preconditioner(preconditioner, Ls)

    return PressureLaplace{unsteady, gradient_mode, acceleration_form, lamb_vorticity, typeof(preconditioner)}(
        Float64(rho), unsteady, Float64(atol), Float64(rtol), Int(itmax),
        preconditioner, Int(reference_panel), Float64(reference_pressure),
        rebuild_every_step, verbose, gradient_mode, acceleration_form,
        lamb_vorticity, kappa_basis, project_edge_du, velocity_dot, velocity_history,
        velocity_history_older, history_count, previous_dt, older_dt, last_step,
        Ls, pressure_operators, bs, ps, residuals,
        acceleration, tangential,
        edges, workspace, u_inertial, omega_used, surface_velocity_gradient,
        jump_velocity, jump_velocity_gradient, vtk_fields, last_rebuild,
        absolute_residual, relative_residual, false, file)
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
                       size(m.u_inertial[i_body]) != size(body.velocity) ||
                       size(m.surface_velocity_gradient[i_body], 3) != body.ncells ||
                       (m.unsteady && size(m.velocity_dot[i_body]) != size(body.velocity))
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

    _pressure_fill_inertial_surface_velocity!(m.u_inertial[i_body], body)
    _pressure_velocity_dot!(m, body, i_body, i_step, dt)
    _pressure_rhs!(m, body, i_body)   # build b from material acceleration
    _pressure_solve!(m, i_body)       # solve L p = b with CG

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
        header = "step,time,system,panels,rebuild,cg_iters,cg_solved,absolute_residual,relative_residual"
        _monitor_csv_open(path, i_step, overwrite, header) do io
            ws = m.workspace[i_system]
            println(io, join((i_step, t, i_system, body.ncells,
                              _monitor_csv_bool(m.last_rebuild[i_system]),
                              ws.stats.niter,
                              _monitor_csv_bool(ws.stats.solved),
                              m.absolute_residual[i_system],
                              m.relative_residual[i_system]), ","))
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

_pressure_velocity_dot!(m::PressureLaplace{false}, body::AbstractBody,
                        i_body::Int, i_step::Int, dt::Real) = nothing

function _pressure_velocity_dot!(m::PressureLaplace{true}, body::AbstractBody,
                                 i_body::Int, i_step::Int, dt::Real)
    dt > 0 || throw(ArgumentError("PressureLaplace history requires positive dt; got $(dt)."))
    vdot = m.velocity_dot[i_body]
    current = m.u_inertial[i_body]
    previous = m.velocity_history[i_body]
    older = m.velocity_history_older[i_body]
    count = m.history_count[i_body]
    consecutive = i_step == m.last_step[i_body] + 1
    if !consecutive || count == 0
        fill!(vdot, 0.0)
    elseif count == 1
        @. vdot = (current - previous) / m.previous_dt[i_body]
    else
        h, hp = m.previous_dt[i_body], m.older_dt[i_body]
        @inbounds for idx in eachindex(vdot)
            vdot[idx] = _variable_bdf2(current[idx], previous[idx], older[idx], h, hp)
        end
    end
    copyto!(older, previous)
    copyto!(previous, current)
    m.history_count[i_body] = consecutive ? min(count + 1, 2) : 1
    m.older_dt[i_body] = m.previous_dt[i_body]
    m.previous_dt[i_body] = Float64(dt)
    m.last_step[i_body] = i_step
    return vdot
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
    residual = m.residual[i_body]
    LA.mul!(residual, A, m.p[i_body])
    residual .-= m.b[i_body]
    absolute_residual = LA.norm(residual)
    relative_residual = absolute_residual / max(LA.norm(m.b[i_body]), eps(Float64))
    m.absolute_residual[i_body] = absolute_residual
    m.relative_residual[i_body] = relative_residual
    if !m.workspace[i_body].stats.solved && !m.convergence_warned
        @warn "PressureLaplace CG did not converge; pressure from capped iterations is under-converged and must not be treated as converged." body=i_body iterations=m.workspace[i_body].stats.niter itmax=m.itmax absolute_residual relative_residual
        m.convergence_warned = true
    end
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
    velocity_dot = m.unsteady ? m.velocity_dot[i_body] : nothing
    if m.acceleration_form == :material_derivative
        _pressure_material_acceleration!(m, body, i_body)
        if m.gradient_mode == :edge_difference
            return _pressure_rhs_from_edge_material_derivative!(m.b[i_body], m, body, i_body,
                velocity_dot)
        end
        return _pressure_rhs_from_acceleration!(m.b[i_body], m, body, i_body,
            m.tangential[i_body])
    elseif m.acceleration_form == :lamb_vector
        return _pressure_rhs_from_lamb_vector!(m.b[i_body], m, body, i_body,
            velocity_dot)
    else
        throw(ArgumentError("Unknown PressureLaplace acceleration_form $(m.acceleration_form)."))
    end
end

# Select which surface vorticity the lamb-vector RHS ingests. The baseline
# (:induced) copies body.induced_vorticity (wake-induced curl + bound-sheet
# kappa, as accumulated by simulate!/replay) and is bit-identical to the
# pre-kwarg behavior; the other variants isolate individual contributions for
# reliability studies.
function _pressure_fill_omega_used!(m::PressureLaplace, body::AbstractBody, i_body::Int)
    omega = m.omega_used[i_body]
    # kappa_basis must match the basis used when κ was accumulated into
    # body.induced_vorticity (simulate!/replay default: :quad) or the
    # :no_bound decomposition would subtract a different κ than was added.
    kappa_options = (; basis=m.kappa_basis)
    if m.lamb_vorticity == :induced
        copyto!(omega, body.induced_vorticity)
    elseif m.lamb_vorticity == :no_bound
        copyto!(omega, body.induced_vorticity)
        _subtract_bound_surface_vorticity!(omega, body; grad_mu_options=kappa_options)
    elseif m.lamb_vorticity == :bound_only
        _bound_surface_vorticity!(omega, body; grad_mu_options=kappa_options)
    elseif m.lamb_vorticity == :hessian_curl
        # curl(u) from the FastMultipole Hessian (layout G[k, l, i]), plus the
        # bound-sheet kappa that the Hessian of off-sheet sources cannot see.
        G = body.velocity_gradient
        @inbounds for i in axes(omega, 2)
            omega[1, i] = G[3, 2, i] - G[2, 3, i]
            omega[2, i] = G[1, 3, i] - G[3, 1, i]
            omega[3, i] = G[2, 1, i] - G[1, 2, i]
        end
        _add_bound_surface_vorticity_into!(omega, body; grad_mu_options=kappa_options)
    else
        throw(ArgumentError("Unknown PressureLaplace lamb_vorticity $(m.lamb_vorticity)."))
    end
    return omega
end

# Exterior volumetric vorticity used only by the corrected edge identity.
# simulate!/replay accumulate wake/particle curl and then add the bound sheet
# κ to body.induced_vorticity for historical monitor compatibility.  κ is a
# surface jump, not ordinary exterior volumetric curl, so remove it here while
# retaining genuine wake and particle contributions.
function _pressure_fill_edge_exterior_vorticity!(m::PressureLaplace,
                                                  body::AbstractBody,
                                                  i_body::Int)
    omega = m.omega_used[i_body]
    copyto!(omega, body.induced_vorticity)
    _subtract_bound_surface_vorticity!(omega, body;
        grad_mu_options=(; basis=m.kappa_basis))
    return omega
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
    omega = _pressure_fill_edge_exterior_vorticity!(m, body, i_body)
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
        uI = m.u_inertial[i_body]
        du1 = uI[1, j] - uI[1, i]
        du2 = uI[2, j] - uI[2, i]
        du3 = uI[3, j] - uI[3, i]
        wx = 0.5 * (omega[1,i] + omega[1,j])
        wy = 0.5 * (omega[2,i] + omega[2,j])
        wz = 0.5 * (omega[3,i] + omega[3,j])
        # rᵀGq = qᵀGr + r·(ω×q), required when G is nonsymmetric.
        vort1 = wy*urel3 - wz*urel2
        vort2 = wz*urel1 - wx*urel3
        vort3 = wx*urel2 - wy*urel1
        convective_edge_dot_r = urel1*du1 + urel2*du2 + urel3*du3 +
                                vort1*r1 + vort2*r2 + vort3*r3

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
    omega_used = _pressure_fill_omega_used!(m, body, i_body)
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

        wx_i, wy_i, wz_i = omega_used[1, i], omega_used[2, i], omega_used[3, i]
        wx_j, wy_j, wz_j = omega_used[1, j], omega_used[2, j], omega_used[3, j]
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
        w, ell, nu1, nu2, nu3 = _pressure_edge_conormal_weight(
            body, edge_a, edge_b, i, j)[1:5]
        aedge_dot_nu = 0.5 * (
            (acceleration[1, i] + acceleration[1, j]) * nu1 +
            (acceleration[2, i] + acceleration[2, j]) * nu2 +
            (acceleration[3, i] + acceleration[3, j]) * nu3)
        flux = rho * ell * aedge_dot_nu
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
    u_inertial = m.u_inertial[i_body]
    if m.unsteady
        copyto!(acceleration, m.velocity_dot[i_body])
    else
        fill!(acceleration, 0.0)
    end

    G_surf = nothing
    if m.gradient_mode == :surface_velocity
        G_surf = m.surface_velocity_gradient[i_body]
        te_info = hasproperty(body, :shedding_full) ?
                  view(body.shedding_full, 1:2, :) :
                  zeros(Int, 2, body.ncells)
        # Structured split-quad meshes use the node-aware agglomerated
        # reconstruction. Genuinely unstructured all-triangle meshes have no
        # logical pairs and retain the triangle reconstruction.
        partners = _quad_pairing_propagate(body.nodes, body.cells, body.neighbor,
            body.normals, te_info)
        use_quads = any(!iszero, partners)
        grad_options = use_quads ? (; basis=:quad) : (; basis=:tri, tri_robust=true)
        compute_surface_velocity_gradient!(G_surf, u_inertial, body.controlpoints,
            body.normals, body.cells, body.neighbor, te_info;
            nodes=body.nodes, grad_mu_options=grad_options)
    end

    G_jump = m.jump_velocity_gradient[i_body]
    fill!(G_jump, 0.0)
    if m.gradient_mode == :corrected_hessian
        jump = m.jump_velocity[i_body]
        fill!(jump, 0.0)
        if has_grad_mu(body) && body.ncells >= 3 && size(m.edges[i_body], 2) > 0
            te_info = hasproperty(body, :shedding_full) ?
                      view(body.shedding_full, 1:2, :) : zeros(Int, 2, body.ncells)
            grad_options = (; basis=m.kappa_basis)
            compute_mu_gradient!(jump, body.controlpoints, body.normals, body.cells,
                body.neighbor, view(body.strength, :, get_Gammai(body)), te_info;
                scale=0.5, nodes=body.nodes, grad_mu_options=grad_options)
            compute_surface_velocity_gradient!(G_jump, jump, body.controlpoints,
                body.normals, body.cells, body.neighbor, te_info;
                nodes=body.nodes, grad_mu_options=grad_options)
        end
    end

    # body.velocity_gradient is already the inertial induced-velocity gradient:
    # kinematic_velocity! subtracts grid velocity only from values, not from this
    # buffer. Add only the derivative of the postprocessed exterior half jump.
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

        if m.gradient_mode == :corrected_hessian
            acceleration[1, i] += ut1 * (G_ind[1,1,i] + G_jump[1,1,i]) + ut2 * (G_ind[1,2,i] + G_jump[1,2,i]) + ut3 * (G_ind[1,3,i] + G_jump[1,3,i])
            acceleration[2, i] += ut1 * (G_ind[2,1,i] + G_jump[2,1,i]) + ut2 * (G_ind[2,2,i] + G_jump[2,2,i]) + ut3 * (G_ind[2,3,i] + G_jump[2,3,i])
            acceleration[3, i] += ut1 * (G_ind[3,1,i] + G_jump[3,1,i]) + ut2 * (G_ind[3,2,i] + G_jump[3,2,i]) + ut3 * (G_ind[3,3,i] + G_jump[3,3,i])
        elseif G_surf !== nothing
            acceleration[1, i] += ut1 * G_surf[1, 1, i] + ut2 * G_surf[1, 2, i] + ut3 * G_surf[1, 3, i]
            acceleration[2, i] += ut1 * G_surf[2, 1, i] + ut2 * G_surf[2, 2, i] + ut3 * G_surf[2, 3, i]
            acceleration[3, i] += ut1 * G_surf[3, 1, i] + ut2 * G_surf[3, 2, i] + ut3 * G_surf[3, 3, i]
        end
    end

    # Apply the surface projection once, after temporal and convective terms
    # have been summed. A blanket P*G*P would omit curvature/product-rule terms.
    tangent = m.tangential[i_body]
    @inbounds for i in 1:n
        nx, ny, nz = body.normals[1,i], body.normals[2,i], body.normals[3,i]
        an = acceleration[1,i]*nx + acceleration[2,i]*ny + acceleration[3,i]*nz
        tangent[1,i] = acceleration[1,i] - an*nx
        tangent[2,i] = acceleration[2,i] - an*ny
        tangent[3,i] = acceleration[3,i] - an*nz
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
                           binning=:control_point,
                           verbose=false, vtk_fields=(:bin_id,))

Consume distributed panel forces from an earlier force monitor, bin them along a
span axis in the requested reference frame, and store latest-step sectional
loading. `components` is a named tuple of unit vectors expressed in the selected
frame, e.g. `components=(lift=Lhat, drag=Dhat)`.

`select` optionally restricts which panels are binned: pass a predicate
`cp_frame -> Bool` evaluated on each control point expressed in the selected
frame (`i_frame`). This is useful to isolate one blade of a multi-blade rotor,
e.g. `select = cp -> cp[2] > 0`. The default (`nothing`) bins every panel.

`binning=:control_point` assigns each selected panel to one bin using its
control point. `binning=:span_overlap` distributes each selected panel's force
across every bin overlapped by the panel's vertex-projected span interval.
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
    binning::Symbol
    select::TS                      # optional predicate cp_frame -> Bool, or nothing
end

monitor_requires(::SpanwiseLoadingMonitor) = (:F,)
monitor_provides(::SpanwiseLoadingMonitor) = (:sectional_F,)

function SpanwiseLoadingMonitor(nbins::Int, i_system::Int;
                                i_frame::Int=-1,
                                components,
                                span_axis=nothing,
                                normalization=NoSectionalNormalization(),
                                per_length::Bool=false,
                                file::Bool=true,
                                TF=Float64,
                                select=nothing,
                                binning::Symbol=:control_point,
                                verbose::Bool=false,
                                vtk_fields::Tuple{Vararg{Symbol}}=(:bin_id,))
    nbins > 0 || throw(ArgumentError("SpanwiseLoadingMonitor requires nbins > 0; got $(nbins)."))
    components isa NamedTuple || throw(ArgumentError("components must be a NamedTuple of unit vectors."))
    !isempty(keys(components)) || throw(ArgumentError("components must contain at least one direction."))
    binning in (:control_point, :span_overlap) ||
        throw(ArgumentError("SpanwiseLoadingMonitor binning must be :control_point or :span_overlap; got $(binning)."))

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
        zeros(Int, nbins), Int[], zeros(TF, 3, nbins), binning, select)
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

_spanwise_bin_index(s::Real, smin::Real, width::Real, nbins::Int) =
    clamp(floor(Int, (s - smin) / width) + 1, 1, nbins)

function _spanwise_panel_vertex_range(body::AbstractBody, panel::Integer,
                                      origin, R_g2f, span_axis)
    smin = Inf
    smax = -Inf
    @inbounds for i_node in axes(body.cells, 1)
        node = body.cells[i_node, panel]
        p_global = SVector{3}(body.nodes[1, node], body.nodes[2, node], body.nodes[3, node])
        p_frame = R_g2f * (p_global - origin)
        s = dot(span_axis, p_frame)
        smin = min(smin, s)
        smax = max(smax, s)
    end
    return smin, smax
end

function _spanwise_accumulate_control_point!(m::SpanwiseLoadingMonitor{TF},
                                             f_frame, panel::Int, span_coord::TF,
                                             smin::TF, span::TF, width::TF) where TF
    b = span > sqrt(eps(TF)) ? _spanwise_bin_index(span_coord, smin, width, m.nbins) : 1
    m.panel_bin_id[panel] = b
    m.counts[b] += 1
    m.bin_force[1, b] += f_frame[1]
    m.bin_force[2, b] += f_frame[2]
    m.bin_force[3, b] += f_frame[3]
    return nothing
end

function _spanwise_accumulate_overlap!(m::SpanwiseLoadingMonitor{TF}, body::AbstractBody,
                                       f_frame, panel::Int, span_coord::TF,
                                       origin, R_g2f, smin::TF, span::TF,
                                       width::TF) where TF
    cp_bin = span > sqrt(eps(TF)) ? _spanwise_bin_index(span_coord, smin, width, m.nbins) : 1
    m.panel_bin_id[panel] = cp_bin
    pmin, pmax = _spanwise_panel_vertex_range(body, panel, origin, R_g2f, m.span_axis)
    pspan = pmax - pmin
    if !(pspan > sqrt(eps(TF))) || !(span > sqrt(eps(TF)))
        return _spanwise_accumulate_control_point!(
            m, f_frame, panel, span_coord, smin, span, width)
    end

    b_first = _spanwise_bin_index(pmin, smin, width, m.nbins)
    b_last = _spanwise_bin_index(pmax, smin, width, m.nbins)
    @inbounds for b in b_first:b_last
        bmin = smin + (b - 1) * width
        bmax = smin + b * width
        overlap = min(pmax, bmax) - max(pmin, bmin)
        overlap > zero(TF) || continue
        w = overlap / pspan
        m.counts[b] += 1
        m.bin_force[1, b] += w * f_frame[1]
        m.bin_force[2, b] += w * f_frame[2]
        m.bin_force[3, b] += w * f_frame[3]
    end
    return nothing
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
    if m.binning == :span_overlap
        ranges = (_spanwise_panel_vertex_range(body, p, origin, R_g2f, m.span_axis)
                  for p in 1:body.ncells if selected[p])
        first_range = first(ranges)
        smin = first_range[1]
        smax = first_range[2]
        for r in ranges
            smin = min(smin, r[1])
            smax = max(smax, r[2])
        end
    else
        smin = minimum(span_coord[p] for p in 1:body.ncells if selected[p])
        smax = maximum(span_coord[p] for p in 1:body.ncells if selected[p])
    end
    span = smax - smin
    width = span > sqrt(eps(TF)) ? span / m.nbins : one(TF)
    fill!(m.counts, 0)
    fill!(m.bin_force, zero(TF))
    fill!(m.panel_bin_id, 0)
    @inbounds for p in 1:body.ncells
        selected[p] || continue
        f_global = SVector{3, TF}(F[1, p], F[2, p], F[3, p])
        f_frame = R_g2f * f_global
        if m.binning == :span_overlap
            _spanwise_accumulate_overlap!(
                m, body, f_frame, p, span_coord[p], origin, R_g2f, smin, span, width)
        else
            _spanwise_accumulate_control_point!(
                m, f_frame, p, span_coord[p], smin, span, width)
        end
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

    # Publish the raw (dimensional, frame-coords) sectional binning so downstream
    # monitors (e.g. DragPolarMonitor) can work in physical units, bypassing the
    # normalization applied to load_components/force_components above.
    monitor_register!(ctx, :sectional_F, m.i_system,
        (bin_center = m.bin_center, bin_width = m.bin_width,
         bin_force = m.bin_force, panel_bin_id = m.panel_bin_id,
         counts = m.counts, span_axis = m.span_axis, i_frame = m.i_frame))
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

`select` optionally restricts which panels contribute to the integrated force
and moment: pass a predicate `cp_frame -> Bool` evaluated on each control point
expressed in `frames[i_frame]` (global coordinates if `i_frame < 0`), e.g.
`select = cp -> sqrt(cp[2]^2 + cp[3]^2) > 0.1R` to exclude a rotor hub.
Excluded panels have their columns of `distributed_force` zeroed before
integration (so VTK output reflects the exclusion too). The mask is computed on
first use and cached while `ncells` is unchanged — the body is assumed rigid in
the selected frame. Default `nothing` includes every panel.

If `verbose=true`, the normalized CF and CM for each step are printed to stdout
with a single `\\t` indent.
"""
mutable struct ForceMonitor{TF, TN, TS} <: AbstractMonitor
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
    select::TS                      # optional predicate cp_frame -> Bool, or nothing
    select_mask::Vector{Bool}       # cached per-panel mask (empty until first use)
    source::Symbol                  # :pressure (compute F from :P) or :context_force (read :F)
end

# The requires-trait is instance-dependent: a :pressure ForceMonitor computes its
# own force from an upstream pressure field (:P), while a :context_force
# ForceMonitor integrates a per-panel force field an upstream monitor registered
# as :F (e.g. DragPolarMonitor's inviscid+drag container).
monitor_requires(m::ForceMonitor) = m.source === :pressure ? (:P,) : (:F,)
monitor_provides(::ForceMonitor) = (:F,)

"""
    _monitor_select_mask!(mask, select, body, frames, i_frame)

Fill `mask` with `select(cp_frame)` evaluated on each control point expressed
in `frames[i_frame]` (global coordinates if `i_frame < 0`). The body is assumed
rigid in the selected frame, so callers may cache the result while `ncells` is
unchanged.
"""
function _monitor_select_mask!(mask::Vector{Bool}, select, body::AbstractBody,
                               frames, i_frame::Int)
    resize!(mask, body.ncells)
    origin, R_g2f = _spanwise_frame_transform(frames, i_frame)
    @inbounds for p in 1:body.ncells
        cp_global = SVector{3, Float64}(body.controlpoints[1, p],
                                        body.controlpoints[2, p],
                                        body.controlpoints[3, p])
        mask[p] = Bool(select(R_g2f * (cp_global - origin)))
    end
    return mask
end

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
                        normalization=NoNormalization(), select=nothing,
                        verbose=false)

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

`select` optionally restricts which panels contribute: pass a predicate
`cp_frame -> Bool` evaluated on each control point expressed in
`frames[i_frame]` (global coordinates if `i_frame < 0`). Edges belonging to
excluded panels are skipped in the summation. The mask is computed on first use
and cached while `ncells` is unchanged — the body is assumed rigid in the
selected frame. Default `nothing` includes every panel.

If `verbose=true`, the normalized CF for each step is printed to stdout with a
single `\\t` indent.
"""
struct KuttaJoukowskiForce{TF, TB, TN, TS} <: AbstractMonitor
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
    select::TS                      # optional predicate cp_frame -> Bool, or nothing
    select_mask::Vector{Bool}       # cached per-panel mask (empty until first use)
end

function KuttaJoukowskiForce(body::AbstractBody, nt::Int, i_system::Int;
                              rho::Real=1.225,
                              backend::AbstractBackend=DirectBackend(),
                              i_frame::Int=-1,
                              TF=Float64,
                              normalization=NoNormalization(),
                              verbose::Bool=false,
                              select=nothing,
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

    return KuttaJoukowskiForce{TF, typeof(backend), typeof(normalization), typeof(select)}(
        force, i_system, i_frame, i_strength, TF(rho), backend, probes,
        edge_node_a, edge_node_b, panel_of_probe, velocity_kinematic,
        normalization, verbose, file, select, Bool[])
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
    if m.select !== nothing && length(m.select_mask) != body.ncells
        _monitor_select_mask!(m.select_mask, m.select, body, frames, m.i_frame)
    end
    use_mask = m.select !== nothing

    Fsum = zero(SVector{3, TF})
    @inbounds for k in eachindex(m.edge_node_a)
        use_mask && !m.select_mask[m.panel_of_probe[k]] && continue
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
                       select=nothing,
                       file::Bool=true,
                       source::Symbol=:pressure)
    source in (:pressure, :context_force) || throw(ArgumentError(
        "ForceMonitor source must be :pressure or :context_force; got $(source)."))
    force = zeros(TF, 3, nt)
    moment = zeros(TF, 3, nt)
    distributed_force = zeros(TF, 3, 0)
    return ForceMonitor{TF, typeof(normalization), typeof(select)}(
        force, moment, distributed_force, i_system, i_frame, normalization,
        correct_kuttacondition, verbose, vtk_fields, file, select, Bool[], source)
end

function _run_monitor!(monitor::ForceMonitor, ctx::MonitorContext, systems, wakes,
                       frames::AbstractVector{<:ReferenceFrame},
                       uinf, i_step::Int, dt::Real, t=nothing)
    !isnothing(t) && monitor_set_time!(ctx, t)
    systems_tuple = _systems_tuple(systems)
    body = systems_tuple[monitor.i_system]

    if size(monitor.distributed_force) != (3, body.ncells)
        monitor.distributed_force = zeros(eltype(monitor.force), 3, body.ncells)
    end
    if monitor.source === :context_force
        # Integrate a per-panel force field an upstream monitor registered as :F.
        F = monitor_field(ctx, :F, monitor.i_system)
        size(F) == (3, body.ncells) || throw(ArgumentError(
            "ForceMonitor(source=:context_force) requires distributed force with " *
            "size (3, $(body.ncells)); got $(size(F))."))
        copyto!(monitor.distributed_force, F)
    else
        pressure = monitor_field(ctx, :P, monitor.i_system)
        fill!(monitor.distributed_force, zero(eltype(monitor.distributed_force)))
        calcfield_F!(monitor.distributed_force, body, calc_areas(body), body.normals, pressure;
                     correct_kuttacondition=monitor.correct_kuttacondition)
    end

    if monitor.select !== nothing
        if length(monitor.select_mask) != body.ncells
            _monitor_select_mask!(monitor.select_mask, monitor.select, body,
                                  frames, monitor.i_frame)
        end
        @inbounds for p in 1:body.ncells
            if !monitor.select_mask[p]
                monitor.distributed_force[1, p] = zero(eltype(monitor.distributed_force))
                monitor.distributed_force[2, p] = zero(eltype(monitor.distributed_force))
                monitor.distributed_force[3, p] = zero(eltype(monitor.distributed_force))
            end
        end
    end

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

"""
    DragPolarMonitor(nt, i_system, polar, chord; i_frame=-1, rho=1.0,
                     inflow_method=:both, backend=DirectBackend(),
                     n_probe=16, probe_radius_chords=1.0, TF=Float64,
                     file=true, verbose=false)

Add sectional viscous (profile-drag) forces from an airfoil drag polar on top of
the inviscid panel loading. Consumes the per-panel force `:F` and the raw
sectional binning `:sectional_F` produced by an upstream `ForceMonitor` and
`SpanwiseLoadingMonitor`, and **re-registers `:F`** as a new per-panel force
container equal to the inviscid force plus the panel-distributed drag. Place a
`ForceMonitor(...; source=:context_force)` after this monitor to integrate total
thrust/torque including drag.

`polar` is a callable `cd = polar(cl, span_coord)` and `chord` a callable
`c = chord(span_coord)`, both taking the span coordinate in frame units (the
same units as the span axis, e.g. meters). `cl` is passed as the magnitude of
the in-plane sectional lift coefficient.

For each span bin the monitor determines the effective in-plane inflow direction
two ways and compares them:
- `:surface` — area-weighted average of `body.velocity` (the body-relative
  apparent fluid velocity `V_induced + U∞ − V_kinematic`) over the bin's panels;
- `:probe` — average of the induced+freestream−kinematic velocity over a ring of
  `n_probe` off-body probe points (radius `probe_radius_chords × chord`) placed
  in the section plane centred on the quarter-chord, evaluated with `backend`.

`inflow_method` selects which drives the drag magnitude: `:surface`, `:probe`,
or `:both` (default; uses the probe inflow and records the surface inflow purely
for the comparison stored in `inflow_angle_diff`).

The lift direction is the in-plane part of the sectional force; the drag
direction is perpendicular to lift, in the section plane, signed to point mostly
along the inflow.
"""
mutable struct DragPolarMonitor{TF, TB, TP, TC} <: AbstractMonitor
    i_system::Int
    i_frame::Int
    rho::TF
    nt::Int
    polar::TP
    chord::TC
    inflow_method::Symbol
    backend::TB
    n_probe::Int
    probe_radius_chords::TF
    probes::FastMultipole.ProbeSystem{TF}
    velocity_kinematic::Matrix{TF}
    initialized::Bool
    nbins::Int
    # --- per-bin state (frame coords unless noted), sized on first call ---
    distributed_force::Matrix{TF}   # 3 × ncells GLOBAL: upstream :F + panel-distributed drag
    drag_force::Matrix{TF}          # 3 × nbins sectional viscous force (frame coords, per bin)
    lift_dir::Matrix{TF}            # 3 × nbins
    drag_dir_surface::Matrix{TF}    # 3 × nbins
    drag_dir_probe::Matrix{TF}      # 3 × nbins
    inflow_surface::Matrix{TF}      # 3 × nbins in-plane inflow (frame coords, with magnitude)
    inflow_probe::Matrix{TF}        # 3 × nbins
    qc_frame::Matrix{TF}            # 3 × nbins quarter-chord (frame coords)
    n_probe_skipped::Vector{Int}    # nbins probes dropped for being too near the surface
    # --- histories ---
    inflow_angle_diff::Matrix{TF}   # nbins × nt angle (deg) between surface & probe inflow
    cl_history::Matrix{TF}          # nbins × nt
    cd_history::Matrix{TF}          # nbins × nt
    bin_center_history::Vector{TF}  # nbins span coordinate of the last call
    force::Matrix{TF}               # 3 × nt viscous-only total force (frame coords, dimensional)
    moment::Matrix{TF}              # 3 × nt viscous-only total moment about frame origin
    file::Bool
    verbose::Bool
end

monitor_requires(::DragPolarMonitor) = (:F, :sectional_F)
monitor_provides(::DragPolarMonitor) = (:F,)

function DragPolarMonitor(nt::Int, i_system::Int, polar, chord;
                          i_frame::Int=-1, rho::Real=1.0,
                          inflow_method::Symbol=:both,
                          backend::AbstractBackend=DirectBackend(),
                          n_probe::Int=16,
                          probe_radius_chords::Real=1.0,
                          TF=Float64, file::Bool=true, verbose::Bool=false)
    inflow_method in (:surface, :probe, :both) || throw(ArgumentError(
        "DragPolarMonitor inflow_method must be :surface, :probe, or :both; got $(inflow_method)."))
    n_probe > 2 || throw(ArgumentError("DragPolarMonitor requires n_probe > 2; got $(n_probe)."))
    empty_probes = FastMultipole.ProbeSystem(0, TF)
    z3 = zeros(TF, 3, 0)
    return DragPolarMonitor{TF, typeof(backend), typeof(polar), typeof(chord)}(
        i_system, i_frame, TF(rho), nt, polar, chord, inflow_method, backend,
        n_probe, TF(probe_radius_chords), empty_probes, zeros(TF, 3, 0),
        false, 0,
        z3, copy(z3), copy(z3), copy(z3), copy(z3), copy(z3), copy(z3), copy(z3),
        Int[],
        zeros(TF, 0, 0), zeros(TF, 0, 0), zeros(TF, 0, 0), TF[],
        zeros(TF, 3, nt), zeros(TF, 3, nt), file, verbose)
end

@inline function _dp_normalize(v::SVector{3, TF}) where TF
    n = norm(v)
    return n > sqrt(eps(TF)) ? v / n : zero(SVector{3, TF})
end

function _dragpolar_frame_transform(frames, i_frame::Int, ::Type{TF}) where TF
    if i_frame < 0
        origin = zero(SVector{3, TF})
        R_f2g = SMatrix{3, 3, TF, 9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
        return origin, R_f2g
    end
    origin_global, R_f2g = frame_global_transform(frames, i_frame)
    return SVector{3, TF}(origin_global), SMatrix{3, 3, TF, 9}(R_f2g)
end

function _dragpolar_ensure_storage!(m::DragPolarMonitor{TF}, nbins::Int, ncells::Int) where TF
    if !m.initialized || m.nbins != nbins
        m.nbins = nbins
        m.drag_force        = zeros(TF, 3, nbins)
        m.lift_dir          = zeros(TF, 3, nbins)
        m.drag_dir_surface  = zeros(TF, 3, nbins)
        m.drag_dir_probe    = zeros(TF, 3, nbins)
        m.inflow_surface    = zeros(TF, 3, nbins)
        m.inflow_probe      = zeros(TF, 3, nbins)
        m.qc_frame          = zeros(TF, 3, nbins)
        m.n_probe_skipped   = zeros(Int, nbins)
        m.bin_center_history = zeros(TF, nbins)
        m.inflow_angle_diff = fill(TF(NaN), nbins, m.nt)
        m.cl_history        = fill(TF(NaN), nbins, m.nt)
        m.cd_history        = fill(TF(NaN), nbins, m.nt)
        m.probes            = FastMultipole.ProbeSystem(nbins * m.n_probe, TF)
        m.velocity_kinematic = zeros(TF, 3, nbins * m.n_probe)
        m.initialized = true
    end
    if size(m.distributed_force) != (3, ncells)
        m.distributed_force = zeros(TF, 3, ncells)
    end
    return nothing
end

# In-plane projector: drop the component along the (unit) span axis.
@inline _dp_inplane(v::SVector{3, TF}, span_axis::SVector{3, TF}) where TF =
    v - dot(v, span_axis) * span_axis

# Most-separated pair of a bin's panel nodes, projected to the bin-center span
# station. Returns the two endpoints in frame coordinates and their distance.
function _dragpolar_chord_endpoints(body, panels::AbstractVector{<:Integer},
                                    origin, R_g2f, span_axis::SVector{3, TF},
                                    bin_center::TF) where TF
    pts = SVector{3, TF}[]
    seen = Set{Int}()
    @inbounds for p in panels
        for i_node in axes(body.cells, 1)
            node = body.cells[i_node, p]
            node in seen && continue
            push!(seen, node)
            x = R_g2f * (SVector{3, TF}(body.nodes[1, node], body.nodes[2, node],
                                        body.nodes[3, node]) - origin)
            # project onto the bin-center plane (normal = span_axis)
            x = x - (dot(x, span_axis) - bin_center) * span_axis
            push!(pts, x)
        end
    end
    length(pts) >= 2 || return (zero(SVector{3, TF}), zero(SVector{3, TF}), zero(TF))
    best = zero(TF); ia = 1; ib = 2
    @inbounds for i in 1:length(pts)-1, j in i+1:length(pts)
        d = norm(pts[i] - pts[j])
        if d > best
            best = d; ia = i; ib = j
        end
    end
    return pts[ia], pts[ib], best
end

# Drag direction: perpendicular to lift, in the section plane, pointing mostly
# along the in-plane inflow. Falls back to span × lift when the inflow is nearly
# parallel to lift (degenerate).
function _dragpolar_drag_dir(v_inflow_ip::SVector{3, TF}, lhat::SVector{3, TF},
                             span_axis::SVector{3, TF}) where TF
    vhat = _dp_normalize(v_inflow_ip)
    resid = vhat - dot(vhat, lhat) * lhat
    if norm(resid) > sqrt(eps(TF))
        return _dp_normalize(resid)
    end
    fallback = _dp_normalize(cross(span_axis, lhat))
    return dot(fallback, vhat) >= 0 ? fallback : -fallback
end

function _run_monitor!(m::DragPolarMonitor{TF}, ctx::MonitorContext, systems, wakes,
                       frames::AbstractVector{<:ReferenceFrame},
                       uinf, i_step::Int, dt::Real, t=nothing) where TF
    !isnothing(t) && monitor_set_time!(ctx, t)
    systems_tuple = _systems_tuple(systems)
    body = systems_tuple[m.i_system]
    sec = monitor_field(ctx, :sectional_F, m.i_system)
    F_up = monitor_field(ctx, :F, m.i_system)
    size(F_up) == (3, body.ncells) || throw(ArgumentError(
        "DragPolarMonitor requires distributed force with size (3, $(body.ncells)); got $(size(F_up))."))

    nbins = length(sec.bin_center)
    _dragpolar_ensure_storage!(m, nbins, body.ncells)
    col = i_step + 1

    origin, R_f2g = _dragpolar_frame_transform(frames, m.i_frame, TF)
    R_g2f = transpose(R_f2g)
    span_axis = _dp_normalize(SVector{3, TF}(sec.span_axis[1], sec.span_axis[2], sec.span_axis[3]))
    areas = calc_areas(body)

    # panels per bin
    bin_panels = [Int[] for _ in 1:nbins]
    @inbounds for p in 1:body.ncells
        b = sec.panel_bin_id[p]
        (1 <= b <= nbins) && push!(bin_panels[b], p)
    end

    do_probe = m.inflow_method in (:probe, :both)

    # --- Phase 1: surface inflow + section geometry per populated bin ---
    fill!(m.inflow_surface, zero(TF))
    fill!(m.qc_frame, zero(TF))
    fill!(m.n_probe_skipped, 0)
    populated = falses(nbins)
    chat = [zero(SVector{3, TF}) for _ in 1:nbins]   # in-plane chord unit dir
    nhat = [zero(SVector{3, TF}) for _ in 1:nbins]   # in-plane normal-to-chord unit dir
    @inbounds for b in 1:nbins
        panels = bin_panels[b]
        isempty(panels) && continue
        populated[b] = true
        bc = TF(sec.bin_center[b])
        m.bin_center_history[b] = bc

        # area-weighted surface inflow (frame coords), in-plane
        vacc = zero(SVector{3, TF}); asum = zero(TF)
        for p in panels
            a = TF(areas[p])
            v = R_g2f * SVector{3, TF}(body.velocity[1, p], body.velocity[2, p], body.velocity[3, p])
            vacc += a * v; asum += a
        end
        v_surf = asum > 0 ? vacc / asum : zero(SVector{3, TF})
        v_surf_ip = _dp_inplane(v_surf, span_axis)
        m.inflow_surface[:, b] .= v_surf_ip

        # chord geometry
        e1, e2, cgeom = _dragpolar_chord_endpoints(body, panels, origin, R_g2f, span_axis, bc)
        if cgeom > sqrt(eps(TF))
            infl_dir = _dp_normalize(v_surf_ip)
            # leading edge = endpoint further upstream (smaller projection on inflow)
            le, te = dot(e1, infl_dir) <= dot(e2, infl_dir) ? (e1, e2) : (e2, e1)
            qc = le + TF(0.25) * (te - le)
            m.qc_frame[:, b] .= qc
            ch = _dp_normalize(_dp_inplane(te - le, span_axis))
            chat[b] = ch
            nhat[b] = _dp_normalize(cross(span_axis, ch))
        end
    end

    # fill empty bins' geometry/inflow from nearest populated bin
    @inbounds for b in 1:nbins
        populated[b] && continue
        nearest = 0; bestdist = typemax(Int)
        for bb in 1:nbins
            populated[bb] || continue
            d = abs(bb - b)
            if d < bestdist; bestdist = d; nearest = bb; end
        end
        nearest == 0 && continue
        m.inflow_surface[:, b] .= m.inflow_surface[:, nearest]
        m.qc_frame[:, b] .= m.qc_frame[:, nearest]
        m.bin_center_history[b] = TF(sec.bin_center[b])
        chat[b] = chat[nearest]; nhat[b] = nhat[nearest]
    end

    # --- Phase 2: probe inflow ---
    if do_probe && any(populated)
        zero_v = zero(SVector{3, TF}); zero_h = zero(SMatrix{3, 3, TF, 9})
        @inbounds for b in 1:nbins
            bc = m.bin_center_history[b]
            radius = m.probe_radius_chords * TF(m.chord(bc))
            qc = SVector{3, TF}(m.qc_frame[1, b], m.qc_frame[2, b], m.qc_frame[3, b])
            cb = chat[b]; nb = nhat[b]
            for kk in 1:m.n_probe
                θ = TF(2π) * (kk - 1) / m.n_probe
                x_frame = qc + radius * (cos(θ) * cb + sin(θ) * nb)
                idx = (b - 1) * m.n_probe + kk
                m.probes.position[idx] = origin + R_f2g * x_frame
                m.probes.gradient[idx] = zero_v
                m.probes.scalar_potential[idx] = zero(TF)
                m.probes.hessian[idx] = zero_h
            end
        end

        _kutta_joukowski_edge_kinematic_velocity!(m.velocity_kinematic, m.probes, m.i_system, frames)

        wake_sources = _collect_wake_sources(wakes)
        all_sources = (systems_tuple..., wake_sources...)
        old_offsets = [sys.kerneloffset for sys in systems_tuple]
        try
            for (i, sys) in pairs(systems_tuple)
                sys.kerneloffset = i == m.i_system ? sys.kerneloffset_panel : sys.kerneloffset_targets
            end
            influence!((m.probes,), all_sources, m.backend;
                       precalc=false, scalar_potential=false, gradient=true, hessian=(false,))
        finally
            for (sys, offset) in zip(systems_tuple, old_offsets)
                sys.kerneloffset = offset
            end
        end

        uinf_sv = SVector{3, TF}(uinf[1], uinf[2], uinf[3])
        fill!(m.inflow_probe, zero(TF))
        @inbounds for b in 1:nbins
            populated[b] || continue
            panels = bin_panels[b]
            radius = m.probe_radius_chords * TF(m.chord(m.bin_center_history[b]))
            skip_r = TF(0.1) * radius / m.probe_radius_chords   # 0.1 * chord
            vacc = zero(SVector{3, TF}); nused = 0; nskip = 0
            for kk in 1:m.n_probe
                idx = (b - 1) * m.n_probe + kk
                pos = m.probes.position[idx]
                # drop probes that fall too close to the bin's surface
                dmin = typemax(TF)
                for p in panels
                    cp = SVector{3, TF}(body.controlpoints[1, p], body.controlpoints[2, p], body.controlpoints[3, p])
                    dmin = min(dmin, norm(pos - cp))
                end
                if dmin < skip_r
                    nskip += 1
                    continue
                end
                v = m.probes.gradient[idx] + uinf_sv -
                    SVector{3, TF}(m.velocity_kinematic[1, idx], m.velocity_kinematic[2, idx], m.velocity_kinematic[3, idx])
                vacc += v; nused += 1
            end
            m.n_probe_skipped[b] = nskip
            v_probe = nused > 0 ? R_g2f * (vacc / nused) : zero(SVector{3, TF})
            m.inflow_probe[:, b] .= _dp_inplane(v_probe, span_axis)
        end
        # fill empty bins from nearest populated
        @inbounds for b in 1:nbins
            populated[b] && continue
            nearest = 0; bestdist = typemax(Int)
            for bb in 1:nbins
                populated[bb] || continue
                d = abs(bb - b)
                if d < bestdist; bestdist = d; nearest = bb; end
            end
            nearest != 0 && (m.inflow_probe[:, b] .= m.inflow_probe[:, nearest])
        end
    end

    # --- Phase 3: drag per bin, panel distribution, histories ---
    copyto!(m.distributed_force, F_up)
    Ftot = zero(SVector{3, TF}); Mtot = zero(SVector{3, TF})
    @inbounds for b in 1:nbins
        Fb = SVector{3, TF}(sec.bin_force[1, b], sec.bin_force[2, b], sec.bin_force[3, b])
        Fb_ip = _dp_inplane(Fb, span_axis)
        lift_mag = norm(Fb_ip)
        lhat = _dp_normalize(Fb_ip)
        m.lift_dir[:, b] .= lhat

        v_surf_ip = SVector{3, TF}(m.inflow_surface[1, b], m.inflow_surface[2, b], m.inflow_surface[3, b])
        v_probe_ip = SVector{3, TF}(m.inflow_probe[1, b], m.inflow_probe[2, b], m.inflow_probe[3, b])
        dhat_surf = _dragpolar_drag_dir(v_surf_ip, lhat, span_axis)
        m.drag_dir_surface[:, b] .= dhat_surf
        dhat_probe = do_probe ? _dragpolar_drag_dir(v_probe_ip, lhat, span_axis) : zero(SVector{3, TF})
        m.drag_dir_probe[:, b] .= dhat_probe

        # comparison angle (deg) between the two inflow directions
        if do_probe && norm(v_surf_ip) > sqrt(eps(TF)) && norm(v_probe_ip) > sqrt(eps(TF))
            cang = clamp(dot(_dp_normalize(v_surf_ip), _dp_normalize(v_probe_ip)), -one(TF), one(TF))
            m.inflow_angle_diff[b, col] = acosd(cang)
        end

        # choose inflow driving the drag magnitude
        v_use = m.inflow_method === :surface ? v_surf_ip : v_probe_ip
        dhat_use = m.inflow_method === :surface ? dhat_surf : dhat_probe
        mag = norm(v_use)
        q = TF(0.5) * m.rho * mag^2
        c = TF(m.chord(m.bin_center_history[b]))
        w = TF(sec.bin_width[b])
        denom = q * c * w
        cl = denom > sqrt(eps(TF)) ? lift_mag / denom : zero(TF)
        cd = TF(m.polar(cl, m.bin_center_history[b]))
        m.cl_history[b, col] = cl
        m.cd_history[b, col] = cd

        Db = denom * cd * dhat_use   # frame coords
        m.drag_force[:, b] .= Db
        Ftot += Db
        qc = SVector{3, TF}(m.qc_frame[1, b], m.qc_frame[2, b], m.qc_frame[3, b])
        Mtot += cross(qc, Db)

        # distribute drag over the bin's panels, area-weighted, in global coords
        panels = bin_panels[b]
        if !isempty(panels)
            asum = zero(TF)
            for p in panels; asum += TF(areas[p]); end
            if asum > 0
                Db_g = R_f2g * Db
                for p in panels
                    frac = TF(areas[p]) / asum
                    m.distributed_force[1, p] += Db_g[1] * frac
                    m.distributed_force[2, p] += Db_g[2] * frac
                    m.distributed_force[3, p] += Db_g[3] * frac
                end
            end
        end
    end

    m.force[:, col] .= Ftot
    m.moment[:, col] .= Mtot

    if m.verbose
        totdrag = sqrt(Ftot[1]^2 + Ftot[2]^2 + Ftot[3]^2)
        angs = filter(!isnan, m.inflow_angle_diff[:, col])
        meanang = isempty(angs) ? NaN : sum(angs) / length(angs)
        println("\t\tDragPolarMonitor[i_system=$(m.i_system), step=$(col)]:")
        println("\t\t\t|F_drag| = $(round(totdrag, sigdigits=4)), mean inflow angle diff = $(round(meanang, sigdigits=4)) deg")
    end

    monitor_register!(ctx, :F, m.i_system, m.distributed_force)
    return nothing
end

function write_monitor_csv!(m::DragPolarMonitor, dir::AbstractString, name::AbstractString,
                            i_monitor::Int, ctx::MonitorContext, systems_tuple::Tuple,
                            i_step::Int, dt::Real; overwrite::Bool=false)
    m.file || return nothing
    path = _monitor_csv_path(dir, name, i_monitor, "dragpolar", m.i_system)
    header = "step,time,bin,bin_center,cl,cd,Dx,Dy,Dz,angle_diff_deg,n_probe_skipped"
    t = monitor_time(ctx, i_step, dt)
    col = i_step + 1
    _monitor_csv_open(path, i_step, overwrite, header) do io
        for b in 1:m.nbins
            println(io, join((i_step, t, b, m.bin_center_history[b],
                              m.cl_history[b, col], m.cd_history[b, col],
                              m.drag_force[1, b], m.drag_force[2, b], m.drag_force[3, b],
                              m.inflow_angle_diff[b, col], m.n_probe_skipped[b]), ","))
        end
    end
    return nothing
end
