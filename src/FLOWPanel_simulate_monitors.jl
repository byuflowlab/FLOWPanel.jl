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
    PressureBernoulli(rho; unsteady=false, correct_kuttacondition=true, clip=nothing)

Monitor that populates `body.P` on every body in the simulation by evaluating
the Bernoulli equation each step. With `unsteady=false` (default) the steady
form ``P = \\tfrac{1}{2} \\rho (U_\\infty^2 - U^2)`` is used; with
`unsteady=true` the term ``-\\rho \\, \\partial \\phi / \\partial t`` is added,
using `body.dphidt` populated by `simulate!`.

`correct_kuttacondition` and `clip` are forwarded to `calcfield_P!`.

Pressure-dependent monitors (e.g. `ForceMonitor`) must appear *after* this
monitor in the `monitors` tuple passed to `simulate!`.
"""
struct PressureBernoulli{TF, TC} <: AbstractMonitor
    rho::TF
    unsteady::Bool
    correct_kuttacondition::Bool
    clip::TC
end

monitor_provides(::PressureBernoulli) = (:P,)

function PressureBernoulli(rho::Real; unsteady::Bool=false,
                          correct_kuttacondition::Bool=true,
                          clip=nothing)
    return PressureBernoulli{typeof(rho), typeof(clip)}(rho, unsteady,
                                                       correct_kuttacondition,
                                                       clip)
end

function (m::PressureBernoulli)(systems, wakes,
                               frames::AbstractVector{<:ReferenceFrame},
                               uinf, i_step::Int, dt::Real)
    systems_tuple = _systems_tuple(systems)
    Uinf_mag = norm(uinf)
    for body in systems_tuple
        fill!(body.P, zero(eltype(body.P)))
        dphidt = m.unsteady ? body.dphidt : zero(body.dphidt)
        calcfield_P!(body.P, body, body.velocity, Uinf_mag, m.rho, dphidt;
                     correct_kuttacondition=m.correct_kuttacondition,
                     clip=m.clip)
    end
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
time-difference buffer of the **inertial** fluid velocity
`u_inertial = body.velocity + body.velocity_kinematic`: between calls it
stores `-u_inertial_old`, then during the call it becomes
`(u_inertial_new - u_inertial_old) / dt`, which is the panel-following rate
`d/dt[u(x_p(t), t)]`. The material derivative used in the RHS is then
`Du/Dt = d/dt[u(x_p,t)] + (u_rel · ∇) u_inertial`, where
`u_rel = body.velocity` is the body-relative slip velocity (on an impermeable
surface `u_rel · n = 0` in the continuous limit).
`velocity_dot` is initialized to zero, so on the first call the unsteady term
is `u_inertial / dt` rather than a true finite difference; if this warm-up
transient matters, pre-populate
`monitor.velocity_dot[i] .= -(body.velocity .+ body.velocity_kinematic)`
before the first call.

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
- `velocity_dot::Vector{Matrix{Float64}}` — rolling buffer; initialized to zero, then between calls holds `-u_inertial_old` (3 × ncells per body); during a call becomes `(u_inertial_new - u_inertial_old) / dt`
- `u_inertial::Vector{Matrix{Float64}}` — scratch buffer for the inertial fluid velocity `body.velocity + body.velocity_kinematic` (3 × ncells per body); used as the scalar field for the surface gradient computation
- `L::Vector{SparseArrays.SparseMatrixCSC{Float64, Int}}` — sparse FV surface Laplacian per body, gauge-fixed at `reference_panel`
- `b::Vector{Vector{Float64}}` — RHS vector per body (length ncells); rebuilt each call from the tangential material acceleration
- `p::Vector{Vector{Float64}}` — solution vector per body (length ncells); written to `body.P` after each solve
- `acceleration::Vector{Matrix{Float64}}` — material acceleration `Du/Dt` scratch buffer (3 × ncells per body)
- `tangential::Vector{Matrix{Float64}}` — tangential projection of `acceleration` (3 × ncells per body)
- `grad::Vector{Array{Float64, 3}}` — surface gradient of each velocity component (3 × 3 × ncells per body; first index is velocity component, second is gradient direction)
- `grad_comp::Vector{Matrix{Float64}}` — scratch buffer for one velocity-component gradient (3 × ncells per body)
- `ATA::Vector{Matrix{Float64}}` — 3 × 3 normal-equation accumulator for the per-panel least-squares gradient solve
- `ATb::Vector{Vector{Float64}}` — length-3 RHS accumulator for the per-panel least-squares gradient solve
- `stencil::Vector{Vector{Int}}` — reusable neighbor-index scratch list for the gradient stencil of the current panel
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
    grad::Vector{Array{Float64, 3}}
    grad_comp::Vector{Matrix{Float64}}
    ATA::Vector{Matrix{Float64}}
    ATb::Vector{Vector{Float64}}
    stencil::Vector{Vector{Int}}
    edges::Vector{Matrix{Int}}
    workspace::Vector{Krylov.CgWorkspace{Float64, Float64, Vector{Float64}}}
    geometry_signature::Vector{UInt64}
    u_inertial::Vector{Matrix{Float64}}
end

monitor_provides(::PressureLaplace) = (:P,)

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
    grad = Array{Float64, 3}[]
    grad_comp = Matrix{Float64}[]
    ATA = Matrix{Float64}[]
    ATb = Vector{Float64}[]
    stencil = Vector{Int}[]
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
    sizehint!(grad, nbodies)
    sizehint!(grad_comp, nbodies)
    sizehint!(ATA, nbodies)
    sizehint!(ATb, nbodies)
    sizehint!(stencil, nbodies)
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
        push!(grad, zeros(Float64, 3, 3, body.ncells))
        push!(grad_comp, zeros(Float64, 3, body.ncells))
        push!(ATA, zeros(Float64, 3, 3))
        push!(ATb, zeros(Float64, 3))
        st = Int[]
        sizehint!(st, 10)
        push!(stencil, st)
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
        acceleration, tangential, grad, grad_comp, ATA, ATb, stencil,
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
    # Adding u_inertial_new = body.velocity + body.velocity_kinematic and dividing
    # by dt gives the panel-following rate d/dt[u(x_p,t)] ≈ (u_new - u_old)/dt.
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
    # the gauge choice (reference_panel), and the per-edge metric quantities ℓ and d.
    # Rigid translation/rotation changes control-point positions and node positions equally,
    # so ℓ and d are unchanged and the signature stays the same — no rebuild needed.
    sig = hash((body.ncells, size(body.cells), body.cells, size(body.neighbor), body.neighbor, reference_panel))
    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        ell = _pressure_edge_length(body.nodes, edge_a, edge_b)
        d = _pressure_controlpoint_distance(body.controlpoints, i, j)
        sig = hash((ell, d), sig)
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

    # For each shared interior edge between panels i and j, the FV weight is
    # w = ℓ/d (edge length / control-point distance).  Each edge contributes
    # four symmetric entries: +w on both diagonals, -w on both off-diagonals.
    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        i == j && continue
        d = _pressure_controlpoint_distance(body.controlpoints, i, j)
        d > 0 || continue
        w = _pressure_edge_length(body.nodes, edge_a, edge_b) / d
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
    tangential = m.tangential[i_body]
    acceleration = m.acceleration[i_body]
    edges = m.edges[i_body]
    rho = m.rho
    reference_panel = m.reference_panel
    reference_pressure = m.reference_pressure
    fill!(b, 0.0)

    # Step 1: compute Du/Dt = ∂u/∂t + (u_t·∇_s)u at every panel control point.
    _pressure_material_acceleration!(m, body, i_body)

    # Step 2: project acceleration onto the local tangent plane (remove normal component).
    # The surface Poisson equation only involves the tangential divergence.
    @inbounds for i in 1:body.ncells
        nx, ny, nz = body.normals[1, i], body.normals[2, i], body.normals[3, i]
        an = acceleration[1, i] * nx + acceleration[2, i] * ny + acceleration[3, i] * nz
        tangential[1, i] = acceleration[1, i] - an * nx
        tangential[2, i] = acceleration[2, i] - an * ny
        tangential[3, i] = acceleration[3, i] - an * nz
    end

    # Step 3: accumulate the edge-integrated source ρ ∮ a_t·n̂ dℓ into b, matching
    # the SPD -∇² FV Laplacian assembled in L (diag +w, off-diag -w, w = ℓ/d).
    # The continuous identity ∫_i ∇²p = ∮_{∂i} ∇p·n̂ dℓ discretizes to a midpoint-
    # average edge flux f = ℓ · ((a_t,i + a_t,j)/2)·ê_ij with outward-from-i = +ê_ij
    # and outward-from-j = -ê_ij, so b[i] += ρ f and b[j] -= ρ f.
    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        r1 = body.controlpoints[1, j] - body.controlpoints[1, i]
        r2 = body.controlpoints[2, j] - body.controlpoints[2, i]
        r3 = body.controlpoints[3, j] - body.controlpoints[3, i]
        d = sqrt(r1*r1 + r2*r2 + r3*r3)
        d > 0 || continue
        ell = _pressure_edge_length(body.nodes, edge_a, edge_b)
        e1, e2, e3 = r1 / d, r2 / d, r3 / d          # unit vector i→j
        ai = tangential[1, i] * e1 + tangential[2, i] * e2 + tangential[3, i] * e3
        aj = tangential[1, j] * e1 + tangential[2, j] * e2 + tangential[3, j] * e3
        flux = ell * 0.5 * (ai + aj)
        b[i] += rho * flux
        b[j] -= rho * flux
        # When reference_pressure ≠ 0 the gauge column removal shifts the RHS.
        if !iszero(reference_pressure)
            w = ell / d
            i == reference_panel && (b[j] += w * reference_pressure)
            j == reference_panel && (b[i] += w * reference_pressure)
        end
    end

    # Step 4: enforce the gauge constraint b[reference_panel] = reference_pressure.
    b[reference_panel] = reference_pressure
    return b
end

function _pressure_material_acceleration!(m::PressureLaplace, body::AbstractBody, i_body::Int)
    n = body.ncells
    acceleration = m.acceleration[i_body]
    velocity_dot = m.velocity_dot[i_body]
    grad = m.grad[i_body]
    grad_comp = m.grad_comp[i_body]
    u_inertial = m.u_inertial[i_body]

    # Form u_inertial = tangent(body.velocity) + body.velocity_kinematic.
    # body.velocity is the body-relative slip velocity; in the continuous limit it is purely
    # tangent, but the solver only enforces u·n = 0 to its tolerance. Stripping the normal
    # residual here keeps it out of the panel-following finite difference and the surface
    # gradient that feeds the convective term.
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
    te_info = hasproperty(body, :shedding_full) ? view(body.shedding_full, 1:2, :) : nothing

    # Compute the surface gradient of each inertial velocity component (u, v, w separately),
    # then pack into grad[comp, direction, panel].
    for comp in 1:3
        scalar = view(u_inertial, comp, :)
        fill!(grad_comp, 0.0)
        _pressure_surface_gradient!(grad_comp, body.controlpoints, body.normals,
            body.cells, body.neighbor, scalar, te_info,
            m.ATA[i_body], m.ATb[i_body], m.stencil[i_body])
        @inbounds for i in 1:n, k in 1:3
            grad[comp, k, i] = grad_comp[k, i]
        end
    end

    # Add the convective term (u_rel · ∇) u_inertial, using the body-relative slip
    # velocity u_rel = body.velocity. On an impermeable surface u_rel · n = 0
    # continuously; the explicit tangent projection here removes any small numerical
    # normal component from the solver.
    @inbounds for i in 1:n
        ur1, ur2, ur3 = body.velocity[1, i], body.velocity[2, i], body.velocity[3, i]
        un = ur1 * body.normals[1, i] + ur2 * body.normals[2, i] + ur3 * body.normals[3, i]
        ut1 = ur1 - un * body.normals[1, i]
        ut2 = ur2 - un * body.normals[2, i]
        ut3 = ur3 - un * body.normals[3, i]
        for comp in 1:3
            acceleration[comp, i] += ut1 * grad[comp, 1, i] +
                                     ut2 * grad[comp, 2, i] +
                                     ut3 * grad[comp, 3, i]
        end
    end

    return acceleration
end

function _pressure_surface_gradient!(grad_out, controlpoints, normals, cells,
                                     neighbors, scalar, te_info,
                                     ATA::Matrix{Float64},
                                     ATb::Vector{Float64},
                                     stencil::Vector{Int})
    # Per-panel weighted least-squares gradient estimate in 3D.
    # For each panel i, solve min ‖A g - b‖ where A rows are displacement vectors
    # to neighbors and b entries are scalar differences.  The normal-direction
    # penalty (large weight on n·g = 0) constrains the gradient to lie in the
    # tangent plane, and a small Tikhonov term regularizes isolated panels.
    n = size(cells, 2)
    @inbounds for i in 1:n
        empty!(stencil)

        # At a trailing edge, exclude the neighbor that shares the TE edge to
        # avoid mixing pressure across the discontinuity.
        is_te = te_info !== nothing && te_info[1, i] > 0 && te_info[2, i] > 0
        te_v1 = is_te ? cells[te_info[1, i], i] : -1
        te_v2 = is_te ? cells[te_info[2, i], i] : -1
        for k in 1:3
            j = neighbors[k, i]
            j <= 0 && continue
            if is_te
                has_v1 = cells[1, j] == te_v1 || cells[2, j] == te_v1 || cells[3, j] == te_v1
                has_v2 = cells[1, j] == te_v2 || cells[2, j] == te_v2 || cells[3, j] == te_v2
                has_v1 && has_v2 && continue
            end
            push!(stencil, j)
        end

        # Expand to a 2-ring stencil when the 1-ring is under-populated or at a TE,
        # mirroring compute_mu_gradient! in postprocess.jl. With only 1-2 face
        # neighbors, ATA is rank-deficient in-plane and the gradient is controlled
        # by Tikhonov noise; pulling in neighbors-of-neighbors restores rank.
        # Preserve TE-edge exclusion when growing across the second ring so we
        # never cross a trailing edge.
        if is_te || length(stencil) < 3
            n_current = length(stencil)
            for s in 1:n_current
                s_idx = stencil[s]
                for k in 1:3
                    nn_idx = neighbors[k, s_idx]
                    nn_idx <= 0 && continue
                    nn_idx == i && continue
                    nn_idx in stencil && continue
                    if is_te
                        has_v1 = cells[1, nn_idx] == te_v1 || cells[2, nn_idx] == te_v1 || cells[3, nn_idx] == te_v1
                        has_v2 = cells[1, nn_idx] == te_v2 || cells[2, nn_idx] == te_v2 || cells[3, nn_idx] == te_v2
                        has_v1 && has_v2 && continue
                    end
                    push!(stencil, nn_idx)
                end
            end
        end

        fill!(ATA, 0.0)
        fill!(ATb, 0.0)
        mean_sq_dist = 0.0
        for j in stencil
            dx = controlpoints[1, j] - controlpoints[1, i]
            dy = controlpoints[2, j] - controlpoints[2, i]
            dz = controlpoints[3, j] - controlpoints[3, i]
            ds = scalar[j] - scalar[i]
            ATA[1, 1] += dx * dx; ATA[1, 2] += dx * dy; ATA[1, 3] += dx * dz
            ATA[2, 1] += dy * dx; ATA[2, 2] += dy * dy; ATA[2, 3] += dy * dz
            ATA[3, 1] += dz * dx; ATA[3, 2] += dz * dy; ATA[3, 3] += dz * dz
            ATb[1] += dx * ds
            ATb[2] += dy * ds
            ATb[3] += dz * ds
            mean_sq_dist += dx^2 + dy^2 + dz^2
        end
        mean_sq_dist = isempty(stencil) ? 1.0 : mean_sq_dist / length(stencil)

        # Penalty term: add penalty * (n ⊗ n) to ATA so that any gradient
        # component along the normal is heavily penalized.  Scaled by mean_sq_dist
        # so the penalty is dimensionally consistent with the geometry.
        nx, ny, nz = normals[1, i], normals[2, i], normals[3, i]
        penalty = 1e4 * mean_sq_dist
        ATA[1, 1] += penalty * nx * nx; ATA[1, 2] += penalty * nx * ny; ATA[1, 3] += penalty * nx * nz
        ATA[2, 1] += penalty * ny * nx; ATA[2, 2] += penalty * ny * ny; ATA[2, 3] += penalty * ny * nz
        ATA[3, 1] += penalty * nz * nx; ATA[3, 2] += penalty * nz * ny; ATA[3, 3] += penalty * nz * nz

        # Small Tikhonov regularization to handle panels with no neighbors.
        reg = 1e-10 * mean_sq_dist
        ATA[1, 1] += reg; ATA[2, 2] += reg; ATA[3, 3] += reg

        A = SMatrix{3, 3, Float64, 9}(
            ATA[1, 1], ATA[2, 1], ATA[3, 1],
            ATA[1, 2], ATA[2, 2], ATA[3, 2],
            ATA[1, 3], ATA[2, 3], ATA[3, 3],
        )
        rhs = SVector{3, Float64}(ATb[1], ATb[2], ATb[3])
        g = A \ rhs
        grad_out[1, i] = g[1]
        grad_out[2, i] = g[2]
        grad_out[3, i] = g[3]
    end
    return grad_out
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
                        TF=Float64, normalization=NoNormalization(), verbose=false)

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

If `verbose=true`, the normalized CF for each step is printed to stdout with a
single `\\t` indent.
"""
struct KuttaJoukowskiForce{TF, TB, TN} <: AbstractMonitor
    force::Matrix{TF}
    i_system::Int
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
        force, i_system, i_strength, TF(rho), backend, probes,
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

    # 6. Optional normalization (moment is not computed; pass zero).
    Mzero = zero(SVector{3, TF})
    CF, _ = m.normalization(Fsum, Mzero, systems, frames, uinf)
    m.force[:, i_step + 1] .= CF

    if m.verbose
        println("\t\tKuttaJoukowskiForce[i_system=$(m.i_system), step=$(i_step+1)]:")
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

    # Populate body.F from body.P (set by a PressureBernoulli earlier in the chain).
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
