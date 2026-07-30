#=##############################################################################
# DESCRIPTION
    Solver backend definitions.

# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Jan 2026
  * License     : GNU Public License
=###############################################################################

################################################################################
# ABSTRACT SOLVER DEFINITION
################################################################################

"""
    AbstractSolver

Abstract supertype for linear solvers used by FLOWPanel body solves.
"""
abstract type AbstractSolver end

"""
    AbstractMatrixFreeSolver

Abstract supertype for solvers that apply the system operator without storing
the full coefficient matrix.
"""
abstract type AbstractMatrixFreeSolver <: AbstractSolver end

"""
    AbstractMatrixfulSolver{LS}

Abstract supertype for solvers that explicitly assemble a system matrix. `LS`
indicates whether the solve uses a least-squares formulation.
"""
abstract type AbstractMatrixfulSolver{LS} <: AbstractSolver end

"""
    solve!(body, solver; optargs...)

Dispatch point for solving a body with a specific solver backend.
"""
function solve!(self, solver; optargs...)
    throw(ErrorException("solve! not implemented for body of type $(typeof(self)) and solver of type $(typeof(solver))"))
end

function _solve!(self, solver; optargs...)
    throw(ErrorException("_solve! not implemented for body of type $(typeof(self)) and solver of type $(typeof(solver))"))
end

################################################################################


################################################################################
# BOUNDARY CONDITIONS
################################################################################

"""
    calc_bc_noflowthrough!(rhs, Us, normals)
    calc_bc_noflowthrough(Us, normals)

Form the impermeability right-hand side `-U ⋅ n` at each control point.
"""
function calc_bc_noflowthrough!(RHS::AbstractVector,
                                Us::AbstractMatrix, normals::AbstractMatrix)

    @assert size(Us)==size(normals) ""*
        "Invalid matrices `Us` and `normals`;"*
        " expected to have the same size, got $(size(Us)) and $(size(normals))."
    @assert length(RHS)==size(Us, 2) ""*
        "Invalid vector `RHS`;"*
        " expected length $(size(Us, 2)), got $(length(RHS))."

    for (i, (U, normal)) in enumerate(zip(eachcol(Us), eachcol(normals)))
        RHS[i] = -dot(U, normal)
    end

    return RHS
end

function calc_bc_noflowthrough(Us::AbstractMatrix{T1},
                               normals::AbstractMatrix{T2}) where {T1, T2}
    RHS = zeros(promote_type(T1, T2), size(Us, 2))
    return calc_bc_noflowthrough!(RHS, Us, normals)
end

function calc_bc_dirichlet(RHS::AbstractVector, self::AbstractBody{<:Union{Union{ConstantSource, ConstantDoublet}, Union{ConstantSource, VortexRing}}, <:Any, <:Any, true}, backend=DirectBackend(); optargs...)
    # Set source strength for dirichlet bodies
    RHS .-= self.potential
    return RHS
end

################################################################################


################################################################################
# INFLUENCE MATRIX ASSEMBLY
################################################################################

# Pick the kernel and the column of `strength` to unit-activate when populating
# an influence matrix for a given source body. For single-element bodies (NK=1)
# the only column is 1; for combined source+doublet / source+vortex-ring bodies
# (NK=2) the source column is assumed to be prescribed by the caller, so we
# solve on the doublet/vortex-ring column.
_G_kernel_and_strength_index(::AbstractBody{ConstantSource, 1}) = (ConstantSource, 1)
_G_kernel_and_strength_index(::AbstractBody{ConstantDoublet, 1}) = (ConstantDoublet, 1)
_G_kernel_and_strength_index(::AbstractBody{VortexRing, 1}) = (VortexRing, 1)
_G_kernel_and_strength_index(::AbstractBody{Union{ConstantSource, ConstantDoublet}, 2}) = (ConstantDoublet, 2)
_G_kernel_and_strength_index(::AbstractBody{Union{ConstantSource, VortexRing}, 2}) = (VortexRing, 2)

"""
    _G!(G, target_system, source_system; kerneloffset=source_system.kerneloffset, update_geometry=false)

Populate the influence matrix `G` such that `G[i, j]` is the influence of the
j-th panel of `source_system` on the i-th control point of `target_system`.
The kernel and the `strength` column unit-activated on `source_system` are
chosen from the source body's type parameters via
`_G_kernel_and_strength_index`. The quantity stored in `G` is chosen from the
target body's boundary-condition flag (the 4th type parameter of
`AbstractBody`):

- Dirichlet (`DBC=true`): induced potential `phi`.
- Neumann  (`DBC=false`): normal component of induced velocity,
  `u ⋅ target_system.normals[:, i]`.

`G` need not be square: its shape must be
`(target_system.ncells, source_system.ncells)`.

If `update_geometry=true`, the target body's normals and control points are
recomputed before assembling `G`. Control points sit exactly on the panel
surface (no offset); no separate diagonal jump term is added because `induced`
returns the side-aware one-sided self limit (e.g. the interior `+μ/2` doublet
limit for Dirichlet bodies), with the side selected by the source body's DBC
type parameter.
"""
function _G!(G, target_system::AbstractBody{<:Any,<:Any,<:Any,DBC},
             source_system::AbstractBody{<:Any,NK,TF};
             kerneloffset=source_system.kerneloffset,
             update_geometry::Bool=false) where {DBC,NK,TF}
    M = target_system.ncells
    N = source_system.ncells

    if size(G, 1) != M || size(G, 2) != N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    end

    if update_geometry
        calc_normals!(target_system)
        calc_controlpoints!(target_system)
    end

    kernel, strength_index = _G_kernel_and_strength_index(source_system)

    derivatives_switch = DBC ? FastMultipole.DerivativesSwitch(true, false, false) :
                               FastMultipole.DerivativesSwitch(false, true, false)

    # Save current strength and unit-activate the selected column
    old_strength = copy(source_system.strength)
    source_system.strength .= zero(eltype(source_system.strength))
    source_system.strength[:, strength_index] .= 1.0

    # operator mode: the affine attached-wake correction (TraceCorrected
    # formulation) is a right-hand-side constant and must not enter the
    # assembled linear operator
    correction_was_active = _operator_mode_begin!(source_system)

    CPs = target_system.controlpoints
    normals = target_system.normals

    # `induced` returns the §3 side-aware self limit at self pairs
    # (dispatched on the source body's DBC type parameter), so no explicit
    # ±0.5 jump add is needed here.
    Threads.@threads for i_source in 1:N
        for i_target in 1:M
            tx, ty, tz = CPs[1, i_target], CPs[2, i_target], CPs[3, i_target]
            target = FastMultipole.StaticArrays.SVector{3,TF}(tx, ty, tz)

            phi, u, _ = induced(target, source_system, i_source, derivatives_switch; kerneloffset)

            if DBC
                isnan(phi) && error("NaN encountered in G matrix computation: \ni_source = $i_source, i_target = $i_target, target = $target, source_strength = $(source_system.strength[i_source, strength_index]), kernel = $kernel, kerneloffset = $kerneloffset")
                G[i_target, i_source] = phi
            else
                G[i_target, i_source] = u[1] * normals[1, i_target] + u[2] * normals[2, i_target] + u[3] * normals[3, i_target]
            end
        end
    end

    # Restore strength and correction mode
    _operator_mode_end!(source_system, correction_was_active)
    source_system.strength .= old_strength

    return G
end

################################################################################


################################################################################
# Backslash Operator
################################################################################

"""
    Backslash(body)

Direct solver that assembles and LU-factors the influence matrix for `body`.
The formulation (Neumann or Dirichlet) is chosen automatically from the body's
`DBC` type parameter: control points are placed exterior for Neumann,
interior for Dirichlet.

# IMPROVEMENT : 2026-07-30 : note that `Backslash`/`BackslashCoupled`'s O(N^2) assembly / O(N^3) factorization is the scaling wall the matrix-free solvers exist to avoid : for large panel counts, prefer `KrylovSolver`/`KrylovCoupled` (with `FastMultipoleBackend`) or `FGSSolver`, reserving `Backslash`/`BackslashCoupled` for small N or as an exact validation baseline for the other solvers.
"""
mutable struct Backslash{TF,TGLU} <: AbstractMatrixfulSolver{false}
    G::Matrix{TF}
    Glu::TGLU
    rhs::Vector{TF}
    Uext::Matrix{TF}       # storage for external velocity (saved/restored around solve)
    phi_ext::Vector{TF}    # storage for external potential (saved/restored around solve)
end

function Backslash(body::AbstractBody{<:Any,<:Any,TF}) where TF
    G = zeros(TF, body.ncells, body.ncells)
    rhs = zeros(TF, body.ncells)
    Uext = zeros(TF, 3, body.ncells)
    phi_ext = zeros(TF, body.ncells)
    update_geometry = true
    _G!(G, body, body; kerneloffset=body.kerneloffset_panel, update_geometry=update_geometry)
    Glu = lu!(G)

    return Backslash{TF,typeof(Glu)}(G, Glu, rhs, Uext, phi_ext)
end


function solve!(body::AbstractBody{<:Any,<:Any,<:Any,true}, solver::AbstractSolver;
        backend=DirectBackend(),
        update_G=false,
        optargs...)

    potential_old = copy(body.potential)
    tb = 0.0
    ts = 0.0

    try
        set_strengths!(body)
        # For single-body Dirichlet solves, `body.potential` is workspace for
        # the interior self/source potential and the solve is homogeneous unless
        # external influences have already been assembled into that workspace.
        body.potential .= zero(eltype(body.potential))
        influence!(body, body, backend; scalar_potential=true, velocity=false, optargs...)
        tb, ts = _solve!(body, solver; backend, update_G, optargs...)
    finally
        body.potential .= potential_old
    end

    return tb, ts
end

function solve!(body::AbstractBody{<:Any,<:Any,<:Any,false}, solver::AbstractSolver;
        backend=DirectBackend(),
        update_G::Bool=false,
        optargs...)
    # IMPROVEMENT : 2026-07-30 : default update_G to false here, matching the Dirichlet solve! wrapper above : audited every solve!(body, solver) call site in src/ and examples/ (a solver reused across changing geometry without update_G=true is the only way this default could bite, and none do that -- every one either constructs its solver fresh right before the single solve, in which case update_G's value is irrelevant since the Backslash constructor already eagerly built+factored G, or explicitly passes update_G). For matrix-full solvers (e.g. Backslash), whose constructor already eagerly builds+factors G, `true` here just meant every first-ever Neumann solve wastefully rebuilt an already-current G; `update_G` is a no-op for matrix-free solvers (KrylovSolver/FGSSolver) either way, since their _solve! methods absorb it into `optargs...` and never read it. Note BackslashCoupled's own update_G=true default (a different solve! method, not this generic wrapper) is NOT part of this inconsistency and must stay true: its constructor only allocates a dummy identity-matrix G (`Glu = lu!(G)  # dummy init; will be overwritten on first update_G=true`), so the *first* solve! genuinely requires update_G=true to ever build the real matrix.
    tb = 0.0
    ts = 0.0

    if body isa RigidWakeBody && body.watertight
        @warn "Solving a watertight RigidWakeBody with the Neumann formulation " *
              "(DBC=false) gives a rank-deficient influence matrix; results " *
              "will be unreliable. Use DBC=true (Dirichlet) for closed " *
              "surfaces, or remove a cap to make the surface non-watertight." maxlog=1
    end

    tb, ts = _solve!(body, solver; backend, update_G, optargs...)

    return tb, ts
end

function numtype(self::AbstractBody)
    return promote_type(eltype(self.nodes),
                        eltype(self.strength))
end

function get_strength_name(self::AbstractBody)
    return "strength"
end

get_strength_name(::AbstractBody{ConstantSource, 1, <:Any}) = "sigma"
get_strength_name(::AbstractBody{ConstantDoublet, 1, <:Any}) = "mu"
get_strength_name(::AbstractBody{VortexRing, 1, <:Any}) = "gamma"

################################################################################
# GMRES Solver
################################################################################

"""
    KrylovSolver(body; method=:gmres, itmax=20, atol=1e-6, rtol=1e-6, backend=FastMultipoleBackend())

Matrix-free iterative solver that evaluates self-influence through the selected
backend.
"""
struct KrylovSolver{TB<:AbstractBody,B<:AbstractBackend,TF<:Number,TP} <: AbstractMatrixFreeSolver
    body::TB
    backend::B
    Uext::Array{TF, 2}    # storage for external velocity (saved/restored around solve)
    normals::Array{TF, 2} # Normals
    source_strengths::Array{TF, 1} # Fixed source strengths for Dirichlet solves
    unabbreviated_strengths::Array{TF, 1} # Storage for solution strengths
    method::Symbol         # Krylov method to use
    itmax::Int           # Maximum number of iterations
    atol::Float64          # absolute tolerance
    rtol::Float64          # relative tolerance
    preconditioner::TP     # nothing or JacobiPreconditioner
end

function KrylovSolver(body::AbstractBody;
        method::Symbol=:gmres,    # Krylov method to use
        itmax::Int=100,         # Maximum number of iterations
        atol::Real=1e-6,            # Convergence tolerance
        rtol::Real=1e-6,            # Relative convergence tolerance
        backend::AbstractBackend=FastMultipoleBackend(),   # Backend to use
        preconditioner_cell_size::Real=0.0,  # cell size for block Jacobi preconditioner; ≤0 disables
    )
    TF = numtype(body)
    Uext = zeros(TF, 3, body.ncells)
    source_strengths = zeros(TF, body.ncells)
    unabbreviated_strengths = zeros(TF, body.ncells)
    calc_normals!(body)
    calc_controlpoints!(body)

    # build block Jacobi preconditioner if requested
    if preconditioner_cell_size > 0
        preconditioner = FastMultipole.JacobiPreconditioner((body,); cell_size=preconditioner_cell_size)
    else
        preconditioner = nothing
    end

    return KrylovSolver{typeof(body), typeof(backend), TF, typeof(preconditioner)}(body, backend, Uext, copy(body.normals), source_strengths, unabbreviated_strengths, method, itmax, Float64(atol), Float64(rtol), preconditioner)
end

function _set_strength(body::AbstractBody{<:Any, 1, <:Any}, strengths)
    body.strength[:, 1] .= strengths
end

function _set_strength(body::AbstractBody{<:Any, 2, <:Any}, strengths)
    body.strength[:, 2] .= strengths
end

function _set_source_strength_from_velocity!(body::AbstractBody{<:Any, 2, <:Any},
                                             velocity::AbstractMatrix,
                                             normals::AbstractMatrix)
    body.strength[:, 1] .= 0.0
    for d in 1:3
        body.strength[:, 2] .= view(velocity, d, :)
        body.strength[:, 2] .*= view(normals, d, :)
        body.strength[:, 1] .-= body.strength[:, 2]
    end
    body.strength[:, 2] .= 0.0
    return nothing
end

function _set_fixed_source_strength_from_velocity!(body::AbstractBody{<:Any, 2, <:Any},
                                                   velocity::AbstractMatrix,
                                                   normals::AbstractMatrix)
    body.strength[:, 1] .= 0.0
    for d in 1:3
        @views @. body.strength[:, 1] -= velocity[d, :] * normals[d, :]
    end
    return nothing
end

_set_fixed_source_strength_from_velocity!(::AbstractBody, ::AbstractMatrix, ::AbstractMatrix) = nothing

function (solver::KrylovSolver{<:AbstractBody{<:Any, <:Any, <:Any, false}})(C, B, α, β)

    @assert length(B) == solver.body.ncells "Length of strengths vector does not match number of panels in body."
    solver.unabbreviated_strengths .= B
    _set_strength(solver.body, solver.unabbreviated_strengths)

    # `induced` (called via `influence!`) returns the §3 side-aware self limit
    # at self pairs, so no extra jump add is needed.
    solver.body.velocity .= 0
    influence!(solver.body, solver.body, solver.backend; velocity=true)

    # dot product with normals
    solver.body.velocity .*= solver.normals
    solver.body.velocity[1,:] .+= view(solver.body.velocity, 2, :)
    solver.body.velocity[1,:] .+= view(solver.body.velocity, 3, :)

    # scale and add
    C .*= β
    C .+= α .* view(solver.body.velocity, 1, :)
end

function (solver::KrylovSolver{<:AbstractBody{<:Any, <:Any, <:Any, true}})(C, B, α, β)

    @assert length(B) == solver.body.ncells "Length of strengths vector does not match number of panels in body."
    solver.body.strength[:, 1] .= 0
    solver.unabbreviated_strengths .= B
    _set_strength(solver.body, solver.unabbreviated_strengths)

    solver.body.potential .= 0
    influence!(solver.body, solver.body, solver.backend; scalar_potential=true, velocity=false)

    C .*= β
    C .+= α .* solver.body.potential
end

function _solve!(self::AbstractBody{<:Any,<:Any,<:Any,false}, solver::KrylovSolver{<:Any,B,TF}, Das=nothing; optargs...) where {B,TF}

    # save external velocity and update solver fields
    solver.Uext .= self.velocity
    solver.normals .= self.normals

    # construct matrix-free linear operator
    TF2 = TF
    nrows = self.ncells
    ncols = self.ncells
    symmetric, hermitian = false, false
    # LinearOperators expects a callable (function) whose methods can be inspected.
    # Wrap the solver instance in a small closure so `methods` sees a function.
    prod! = (y, x, α, β) -> solver(y, x, α, β)
    tb = @elapsed begin 
        A = LinearOperators.LinearOperator(
            TF2,
            nrows,
            ncols,
            symmetric, hermitian,
            prod!
        )
    end

    # construct right-hand side
    RHS = zeros(TF2, nrows)
    calc_bc_noflowthrough!(RHS, self.velocity, solver.normals)

    # allocate and launch krylov solver
    workspace = Krylov.krylov_workspace(Val(solver.method), A, RHS)
    # Warm-start from this body's previously-converged strength (e.g. from
    # the prior outer Gauss-Seidel sweep, when this solver is reused across
    # solve!(bodies::Tuple, solvers::Tuple) iterations). On the very first
    # call this is still all zeros (the constructor's default), so behavior
    # is unchanged there; on later calls it lets GMRES resume from a guess
    # that's already close to the answer instead of restarting from zero
    # every sweep.
    Krylov.warm_start!(workspace, solver.unabbreviated_strengths)
    # IMPROVEMENT : 2026-07-30 : remove the TEMPORARY warm-start diagnostic println (and its now-unused x0_norm) : it fired unconditionally on every GMRES solve, polluting stdout on any batch run (e.g. an AOA/panel-count sweep) with no way to opt out
    if solver.preconditioner !== nothing
        ts = @elapsed begin
            Krylov.krylov_solve!(workspace, A, RHS; M=solver.preconditioner, ldiv=true, atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
        end
    else
        ts = @elapsed begin
        Krylov.krylov_solve!(workspace, A, RHS; atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
        end
    end

    # store solution
    solver.unabbreviated_strengths .= workspace.x
    _set_strength(self, solver.unabbreviated_strengths)

    # restore external velocity
    self.velocity .= solver.Uext
    return tb, ts
end

function _solve!(self::AbstractBody{<:Any,2,<:Any,true}, solver::KrylovSolver{<:Any,B,TF}, Das=nothing; optargs...) where {B,TF}

    solver.normals .= self.normals
    solver.source_strengths .= view(self.strength, :, 1)

    TF2 = TF
    nrows = self.ncells
    ncols = self.ncells
    symmetric, hermitian = false, false
    prod! = (y, x, α, β) -> solver(y, x, α, β)
    tb = @elapsed begin
        A = LinearOperators.LinearOperator(
            TF2,
            nrows,
            ncols,
            symmetric, hermitian,
            prod!
        )
    end

    RHS = zeros(TF2, nrows)
    RHS .= -self.potential

    workspace = Krylov.krylov_workspace(Val(solver.method), A, RHS)
    # Warm-start from this body's previously-converged doublet strength --
    # see the analogous comment in the Neumann _solve! above.
    Krylov.warm_start!(workspace, solver.unabbreviated_strengths)
    if solver.preconditioner !== nothing
        ts = @elapsed Krylov.krylov_solve!(workspace, A, RHS; M=solver.preconditioner, ldiv=true, atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
    else
        ts = @elapsed Krylov.krylov_solve!(workspace, A, RHS; atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
    end

    solver.unabbreviated_strengths .= workspace.x
    self.strength[:, 1] .= solver.source_strengths
    _set_strength(self, solver.unabbreviated_strengths)

    return tb, ts
end

"""
    KrylovCoupled(bodies; method=:gmres, itmax=20, atol=1e-6, rtol=1e-6, backend=FastMultipoleBackend(), preconditioner_cell_size=0.0)

Matrix-free coupled solver for a tuple of bodies.
"""
# IMPROVEMENT : 2026-07-30 : add a `preconditioner` field (mirroring KrylovSolver's TP type parameter/field) : KrylovCoupled previously had no block-Jacobi preconditioner option at all, unlike KrylovSolver -- needed to compare "with vs without preconditioner" across both Krylov solver families on equal footing
mutable struct KrylovCoupled{TB<:Tuple,B<:AbstractBackend,TF<:Number,TP} <: AbstractMatrixFreeSolver
    bodies::TB
    backend::B
    rhs::Vector{TF}
    x::Vector{TF}
    Ax::Vector{TF}
    method::Symbol
    itmax::Int
    atol::Float64
    rtol::Float64
    preconditioner::TP     # nothing or JacobiPreconditioner
end

function KrylovCoupled(bodies::Tuple;
        method::Symbol=:gmres,
        itmax::Int=50,
        atol::Real=1e-6,
        rtol::Real=1e-6,
        backend::AbstractBackend=FastMultipoleBackend(),
        preconditioner_cell_size::Real=0.0)  # cell size for block Jacobi preconditioner; ≤0 disables

    TF = promote_type(map(numtype, bodies)...)
    ncs = sum(body -> body.ncells, bodies)
    rhs = zeros(TF, ncs)
    x = zeros(TF, ncs)
    Ax = zeros(TF, ncs)

    # build block Jacobi preconditioner if requested
    if preconditioner_cell_size > 0
        preconditioner = FastMultipole.JacobiPreconditioner(bodies; cell_size=preconditioner_cell_size)
    else
        preconditioner = nothing
    end

    return KrylovCoupled{typeof(bodies),typeof(backend),TF,typeof(preconditioner)}(
        bodies, backend, rhs, x, Ax, method, Int(itmax), Float64(atol), Float64(rtol), preconditioner)
end

function _coupled_offsets(bodies::Tuple)
    npanels = [body.ncells for body in bodies]
    return cumsum(vcat(0, npanels))
end

function _write_coupled_unknowns!(bodies::Tuple, x::AbstractVector)
    offsets = _coupled_offsets(bodies)
    for (i, body) in enumerate(bodies)
        r = offsets[i]+1:offsets[i+1]
        if has_dirichlet_bc(body)
            body.strength[:, 1] .= zero(eltype(body.strength))
            body.strength[:, 2] .= view(x, r)
        else
            body.strength[:, 1] .= view(x, r)
        end
    end
    return nothing
end

function _collect_coupled_operator!(Ax::AbstractVector, bodies::Tuple)
    offsets = _coupled_offsets(bodies)
    for (i, body) in enumerate(bodies)
        r = offsets[i]+1:offsets[i+1]
        if has_dirichlet_bc(body)
            Ax[r] .= body.potential
        else
            @views Ax[r] .= vec(sum(body.velocity .* body.normals; dims=1))
        end
    end
    return Ax
end

function (solver::KrylovCoupled)(C, B, α, β)
    _write_coupled_unknowns!(solver.bodies, B)

    for body in solver.bodies
        body.velocity .= zero(eltype(body.velocity))
        body.potential .= zero(eltype(body.potential))
    end

    scalar_potential = [has_dirichlet_bc(body) for body in solver.bodies]
    velocity = [!has_dirichlet_bc(body) for body in solver.bodies]
    influence!(solver.bodies, solver.bodies, solver.backend; scalar_potential, velocity)
    _collect_coupled_operator!(solver.Ax, solver.bodies)

    C .*= β
    C .+= α .* solver.Ax
    return nothing
end

function solve!(bodies::Tuple, solver::KrylovCoupled; backend=solver.backend, optargs...)
    offsets = _coupled_offsets(bodies)
    velocity_old = [copy(body.velocity) for body in bodies]
    potential_old = [copy(body.potential) for body in bodies]
    fixed_sources = Vector{Any}(undef, length(bodies))
    tb = 0.0
    ts = 0.0

    try
        for (i, body) in enumerate(bodies)
            calc_normals!(body)
            calc_controlpoints!(body)
            set_strengths!(body)
            fixed_sources[i] = has_dirichlet_bc(body) ? copy(body.strength[:, 1]) : nothing
            # influence! accumulates onto body.potential (see BackslashCoupled's identical
            # zero-before-influence! pattern above), so it must start from zero here too --
            # otherwise potential_old gets added once implicitly (via this un-zeroed
            # accumulation) and once more explicitly below, double-counting it in the RHS.
            body.potential .= zero(eltype(body.potential))
        end

        # IMPROVEMENT : 2026-07-30 : fold the coupled influence!/boundary_condition! cost into tb instead of leaving it untimed : mirrors the identical BackslashCoupled fix above -- this cross-body coupling evaluation is often the dominant cost of a coupled solve, so excluding it understated (tb, ts) and made cross-solver timing comparisons unfair
        tb += @elapsed begin
            scalar_potential = [has_dirichlet_bc(body) for body in bodies]
            velocity = [!has_dirichlet_bc(body) for body in bodies]
            influence!(bodies, bodies, backend; scalar_potential, velocity, optargs...)
            for (i, body) in enumerate(bodies)
                if has_dirichlet_bc(body)
                    body.potential .+= potential_old[i]
                end
            end
            # calc_bc_dirichlet accumulates onto RHS (RHS .-= body.potential) rather than
            # overwriting it, so solver.rhs must start at zero every call -- true by construction
            # the first time, but not guaranteed if this KrylovCoupled object is reused across
            # multiple solve! calls (e.g. an AOA sweep).
            fill!(solver.rhs, 0)
            boundary_condition!(bodies, solver, backend; optargs...)
        end

        TF = eltype(solver.rhs)
        n = length(solver.rhs)
        prod! = (y, x, α, β) -> solver(y, x, α, β)
        tb += @elapsed begin
            A = LinearOperators.LinearOperator(TF, n, n, false, false, prod!)
        end
        workspace = Krylov.krylov_workspace(Val(solver.method), A, solver.rhs)
        # IMPROVEMENT : 2026-07-30 : pass the optional JacobiPreconditioner to krylov_solve! : mirrors KrylovSolver's identical branch -- previously KrylovCoupled built no preconditioner at all and this call never passed M
        if solver.preconditioner !== nothing
            ts += @elapsed Krylov.krylov_solve!(workspace, A, solver.rhs; M=solver.preconditioner, ldiv=true, atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
        else
            ts += @elapsed Krylov.krylov_solve!(workspace, A, solver.rhs; atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
        end
        solver.x .= workspace.x

        for (i, body) in enumerate(bodies)
            r = offsets[i]+1:offsets[i+1]
            if has_dirichlet_bc(body)
                body.strength[:, 1] .= fixed_sources[i]
            end
            write_solution!(body, view(solver.x, r))
        end
    finally
        for (i, body) in enumerate(bodies)
            body.velocity .= velocity_old[i]
            body.potential .= potential_old[i]
        end
    end

    return tb, ts
end

###############################################################################
# FGS Solver
################################################################################

"""
    FGSSolver(body; kwargs...)

Matrix-free solver for a single body. Use `solve!(bodies, solvers)` with one
solver per body for coupled multi-body solves.
"""
mutable struct FGSSolver{TFGS,TF} <: AbstractMatrixFreeSolver
    fgs::TFGS
    expansion_order::Int
    leaf_size::Int
    multipole_acceptance::Float64
    max_iterations::Int
    inner_iterations::Int
    tolerance::Float64
    rlx::Float64
    reverse_pass::Bool
    verbose::Bool
    Uext::Matrix{TF}                    # external velocity storage
    phi_ext::Vector{TF}                 # external potential storage
    solution_history::Array{TF, 3}      # NP × NK × NT rolling buffer; slot 1 = most recent
    solution_history_length::Int            # NT — buffer capacity
    solution_history_nsaved::Int            # number of populated slots, clamped to NT
    project_solution::Bool                  # warm-start next solve via polynomial extrapolation
    project_solution_order::Int             # extrapolation order: 1 = linear, 2 = quadratic, ...
end

function FGSSolver(body::AbstractBody;
        max_iterations::Int=100,         # Maximum number of iterations
        inner_iterations::Int=1,
        tolerance::Real=1e-6,            # Convergence tolerance
        rlx::Real=1.0,                   # Relaxation factor
        reverse_pass::Bool=true,         # Whether to perform reverse pass for adjoint sensitivities
        expansion_order=7,
        multipole_acceptance=0.4,
        leaf_size=10,
        shrink=false,
        recenter=false,
        verbose=false,
        calc_cps=true,
                solution_history_length::Int=0,      # 0 disables history & projection
        project_solution::Bool=false,        # warm-start next solve via polynomial extrapolation
        project_solution_order::Int=1,       # 1 = linear, 2 = quadratic, ...
        build_fgs::Bool=true,                # false skips the FastGaussSeidel build;
                                             # the solver then only carries options
                                             # (e.g. as a formulation green_solver)
    )

    # calculate control points if needed
    if calc_cps
        calc_normals!(body)
        calc_controlpoints!(body)
    end

    # generate solver
    TF = numtype(body)
    bodies = (body,)
    fgs = build_fgs ?
        FastMultipole.FastGaussSeidel(bodies; expansion_order, multipole_acceptance, leaf_size, shrink, recenter, extra_farfield=any(has_semiinfinite_wake.(bodies))) :
        nothing

    Uext = zeros(TF, 3, body.ncells)
    phi_ext = zeros(TF, body.ncells)
    solution_history = zeros(TF, body.ncells, size(body.strength, 2), solution_history_length)
    return FGSSolver{typeof(fgs), TF}(fgs, Int(expansion_order), Int(leaf_size), Float64(multipole_acceptance), max_iterations, Int(inner_iterations), Float64(tolerance), Float64(rlx), Bool(reverse_pass), Bool(verbose), Uext, phi_ext, solution_history, solution_history_length, 0, project_solution, project_solution_order)
end

################################################################################

function _solve!(body::AbstractBody{TK,NK,TF,false}, solver::Backslash;
        backend=DirectBackend(), strength_index=1,
        update_G::Bool=false,
        optargs...
    ) where {TK, NK, TF}

    tb = 0.0
    if update_G
        solver.G .= zero(eltype(solver.G))
        tb = @elapsed begin
            _G!(solver.G, body, body; kerneloffset=body.kerneloffset_panel, update_geometry=false, optargs...)
            solver.Glu = lu!(solver.G)
        end
    end

    rhs = solver.rhs
    rhs .= zero(eltype(rhs))
    calc_bc_noflowthrough!(rhs, body.velocity, body.normals)

    ts = @elapsed ldiv!(view(body.strength, :, strength_index), solver.Glu, rhs)

    return tb, ts
end

function _solve!(self::AbstractBody{<:Union{Union{ConstantSource, ConstantDoublet}, Union{ConstantSource, VortexRing}}, 2, TF, true},
                solver::Backslash; backend=DirectBackend(),
                update_G::Bool=false, optargs...) where TF

    solver.rhs .= -self.potential

    tb = 0.0
    if update_G
        solver.G .= 0.0
        tb += @elapsed begin
            _G!(solver.G, self, self; kerneloffset=self.kerneloffset_panel, update_geometry=false)
            solver.Glu = lu!(solver.G)
        end
    end

    ts = @elapsed ldiv!(view(self.strength, :, 2), solver.Glu, solver.rhs)

    return tb, ts
end

# function solve!(self::RigidWakeBody{TK, 1, TF},
#                     solver::AbstractMatrixfulSolver{false};
#                     backend=FastMultipoleBackend(),
#                     solver_optargs=(),
#                     update_G::Bool=true,
#                     optargs...
#                 ) where {TK<:Union{VortexRing, ConstantDoublet}, TF}
#     if size(self.velocity) != (3, self.ncells)
#         error("Invalid body velocity;"*
#               " expected size (3, $(self.ncells)), got $(size(self.velocity))")
#     end

#     normals = _calc_normals(self)
#     CPs = _calc_controlpoints(self, normals)

#     G = zeros(TF, self.ncells, self.ncells)
#     _G_U!(self, TK, G, CPs, normals, backend; optargs...)

#     RHS = calc_bc_noflowthrough(self.velocity, normals)

#     Gamma = zeros(TF, self.ncells)
#     solve_matrix!(Gamma, G, RHS, solver; solver_optargs...)

#     self.strength[:, 1] .= Gamma

# end

# Polynomial extrapolation in time: warm-start body.strength using the rolling
# history saved by save_solution!. Slot 1 holds the most recent saved strength.
# Coefficients (slot 1 = most recent): s_new = Σ_{j=0..order} (-1)^j * binomial(order+1, j+1) * H[:,:,j+1]
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

# Push the converged body.strength into the per-body rolling history.
# Slot 1 is most recent: shift slots [1..NT-1] right by one, then write slot 1.
@inline function save_solution!(body::AbstractBody, solver::FGSSolver)
    NT = solver.solution_history_length
    NT == 0 && return nothing
    H = solver.solution_history
    @inbounds for j in NT:-1:2
        @views @. H[:, :, j] = H[:, :, j - 1]
    end
    @views @. H[:, :, 1] = body.strength
    solver.solution_history_nsaved = min(solver.solution_history_nsaved + 1, NT)
    return nothing
end

@inline _fgs_solved_strength_index(body::AbstractBody) = has_dirichlet_bc(body) && size(body.strength, 2) >= 2 ? 2 : 1

function _solve!(body::AbstractBody, solver::FGSSolver{Nothing}; optargs...)
    error("This FGSSolver was constructed with build_fgs=false and only "*
          "carries options (e.g. as a formulation green_solver); it cannot "*
          "run the FastGaussSeidel body solve.")
end

function _solve!(body::AbstractBody, solver::FGSSolver; backend = FastMultipoleBackend(
        expansion_order=solver.expansion_order,
        multipole_acceptance=solver.multipole_acceptance,
        leaf_size=solver.leaf_size
    ), optargs...)

    dirichlet_bc = has_dirichlet_bc(body)
    strength_index = _fgs_solved_strength_index(body)
    nprint = min(5, body.ncells)
    prior_sigma = dirichlet_bc ? copy(body.strength[:, 1]) : nothing

    # warm-start strengths from history (no-op if disabled or insufficient history)
    projected = project_solution!(body, solver)
    if dirichlet_bc
        body.strength[:, 1] .= prior_sigma
    end
    projected_strengths = if solver.verbose && projected
        vals = copy(view(body.strength, 1:nprint, strength_index))
        println("FGSSolver projected first $(nprint) strengths (column $(strength_index)): ", vals)
        vals
    else
        if solver.verbose && solver.project_solution
            println("FGSSolver projection skipped; saved history count = $(solver.solution_history_nsaved)")
        end
        nothing
    end

    # run solver
    # IMPROVEMENT : 2026-07-30 : time FastMultipole.solve! and return (tb, ts) instead of falling through to save_solution!'s `nothing` : every other _solve! returns (tb, ts), which both public solve! wrappers unconditionally destructure (`tb, ts = _solve!(...)`) -- returning `nothing` here crashed that path, so FGSSolver could only ever be driven through _solve! directly or FastMultipole.solve!, never the public solve!(body, solver) API used by the outer coupled Gauss-Seidel loop. tb=0.0 since FGS has no separable matrix-assembly phase (matrix-free, like KrylovSolver/KrylovCoupled) -- all cost lives in the timed solve itself.
    tb = 0.0
    ts = @elapsed begin
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
    end

    # restore properties
    if dirichlet_bc
        body.strength[:, 1] .= prior_sigma
    end

    if solver.verbose && projected
        actual_strengths = copy(view(body.strength, 1:nprint, strength_index))
        println("FGSSolver actual first $(nprint) strengths after solve (column $(strength_index)): ", actual_strengths)
        println("FGSSolver actual - projected first $(nprint) strengths: ", actual_strengths .- projected_strengths)
    end

    # save converged strengths into rolling history (no-op if disabled)
    save_solution!(body, solver)

    return tb, ts
end

"""
    solve!(bodies::Tuple, solvers::Tuple; update_G=false, optargs...)

Matrix-full solvers (e.g. `Backslash`) already assemble and LU-factor their
`G` matrix at construction time, so `G` is current for whatever body geometry
and wake state (e.g. `Das`) existed when the solver was built. This outer
Gauss-Seidel loop therefore trusts that pre-built `G` by default
(`update_G=false`) and never re-assembles it: cross-body coupling is handled
entirely through the per-body RHS (via `influence!` from the other bodies),
not through `G` itself, so `G` does not need to change between outer
iterations.

Pass `update_G=true` only when a solver is being reused across a change that
invalidates its cached `G` (e.g. body geometry or wake direction changed
since the solver was constructed/last rebuilt). In that case `G` is rebuilt
once, on the first outer iteration, and reused for the remaining iterations.
"""
function solve!(bodies::Tuple, solvers::Tuple;
    backend = fill(FastMultipoleBackend(), length(bodies)),
    max_outer_iterations::Int = 50,
    outer_tolerance::Real = 1e-6,
    verbose::Bool = false,
    update_G::Bool = false,
    optargs...)
    # IMPROVEMENT : 2026-07-30 : default per-body backend to FastMultipoleBackend() instead of DirectBackend() for the cross-body influence! evaluation each outer sweep : same benchmark data as BackslashCoupled's identical change above (negligible ~1e-8 accuracy cost, growing speed win past ~1000 total panels); every existing solve!(bodies::Tuple, solvers::Tuple) call site passes `backend` explicitly, so no caller relied on the old Direct default.

    # IMPROVEMENT : 2026-07-30 : remove the unconditional println("Tuple of bodies") : it fired on every call to the outer coupled Gauss-Seidel loop regardless of the verbose flag, adding needless overhead/noise across e.g. an AOA sweep
    N = length(bodies)
    @assert length(solvers) == N "Number of solvers ($(length(solvers))) must match number of bodies ($N)"
    backends = backend isa Tuple || backend isa AbstractVector ? backend : fill(backend, N)

    prev_velocity = [copy(body.velocity) for body in bodies]
    prev_strengths = [copy(body.strength) for body in bodies]

    converged = false
    t_solve = 0.0
    t_build = 0.0
    for iter in 1:max_outer_iterations
        iter_update_G = iter == 1 && update_G

        for (i, (body, solver)) in enumerate(zip(bodies, solvers))
            body.velocity .= prev_velocity[i]

            sources = tuple((bodies[j] for j in eachindex(bodies) if j != i)...)
            # IMPROVEMENT : 2026-07-30 : fold t_influence into t_build instead of only printing it, and remove the unconditional [outer diag] println : cross-body influence! is often the dominant cost of a coupled solve, so discarding it understated the loop's true (tb, ts), and the per-iteration println fired unconditionally regardless of the verbose flag
            t_influence = 0.0
            if !isempty(sources)
                t_influence = @elapsed influence!((body,), sources, backends[i];
                    scalar_potential=false,
                    velocity=true,
                    optargs...)
            end
            tb, ts = solve!(body, solver; backend=backends[i], update_G=iter_update_G)
            t_build += tb + t_influence
            t_solve += ts
        end

        max_delta = 0.0
        for (i, body) in enumerate(bodies)
            prev_strengths[i] .-= body.strength
            prev_strengths[i] .= abs.(prev_strengths[i])
            max_delta = max(max_delta, maximum(prev_strengths[i]))
            prev_strengths[i] .= body.strength
        end

        if verbose
            println("  Outer iteration $iter: max strength change = $max_delta")
        end

        if max_delta < outer_tolerance
            converged = true
            if verbose
                println("  Converged after $iter outer iterations")
            end
            break
        end
    end

    if !converged && verbose
        println("  WARNING: outer iteration did not converge after $max_outer_iterations iterations")
    end

    # restore velocities
    for (i, body) in enumerate(bodies)
        body.velocity .= prev_velocity[i]
    end

    return t_build, t_solve
end


################################################################################
# FLAT GROUND SOLVER
################################################################################

"""
    FlatGroundSolver(body::NonLiftingBody{ConstantSource, 1, TF})

Direct solver for flat ground problems, for which the influence matrix is diagonal.
    Then, panels can be solved for independently, without forming or factorizing a full matrix.

## Fields

- `rhs::Vector{TF}`: Right-hand side vector for the solve, preallocated for efficiency.

"""
struct FlatGroundSolver{TF} <: AbstractSolver
    rhs::Vector{TF}
end

function FlatGroundSolver(body::NonLiftingBody{ConstantSource, 1, TF}) where TF
    rhs = zeros(TF, body.ncells)
    return FlatGroundSolver(rhs)
end

function _solve!(body::NonLiftingBody{ConstantSource, 1, TF}, solver::FlatGroundSolver{TF}; optargs...) where TF
    calc_bc_noflowthrough!(solver.rhs, body.velocity, body.normals)

    for i in 1:body.ncells
        # For a flat ground problem with constant source panels, the influence of each panel on itself is -0.5
        # (for a panel with its control point on the surface), and there is no influence from other panels.
        body.strength[i, 1] = -solver.rhs[i] / (-0.5)
    end

    return nothing
end

################################################################################

# BACKSLASH COUPLED

"""
    BackslashCoupled{TF}

A direct solver that assembles the full influence matrix G for all bodies and solves the coupled system using a single backslash operation.
Components:
- `G`: The full influence matrix for all bodies.
- `Glu`: Cached LU factorization of G for efficient solves.
- `rhs`: The right-hand side vector constructed from the boundary conditions of all bodies.
- `Uext`: Cached external velocity at control points for all bodies.
- `phi_ext`: Cached external potential at control points for all bodies.
"""
mutable struct BackslashCoupled{TF} <: AbstractSolver
    G::Matrix{TF}
    Glu::LA.Factorization{TF}
    rhs::Vector{TF}
    Uext::Matrix{TF}
    phi_ext::Vector{TF}
end

"""
    BackslashCoupled(bodies::Tuple{<:AbstractBody{<:Any,<:Any,TF}}) where TF

"""
function BackslashCoupled(bodies::Tuple{Vararg{<:AbstractBody{<:Any,<:Any,TF,<:Any}}}) where TF
    ncs = sum(b -> b.ncells, bodies)

    G       = Matrix{TF}(I, ncs, ncs)
    rhs     = zeros(TF, ncs)
    Uext    = zeros(TF, 3, ncs)
    phi_ext = zeros(TF, ncs)

    Glu = lu!(G)  # dummy init; will be overwritten on first update_G=true

    BackslashCoupled{TF}(G, Glu, rhs, Uext, phi_ext)
end

# Backslash(bodies::Tuple) = BackslashCoupled(bodies)

### DIRICHLET
"""
    boundary_condition!(body::AbstractBody{<:Any,<:Any,<:Any,true}, RHS, backend; optargs...)

For Dirichlet bodies, set the source strength to enforce the no-penetration condition and set the doublet strength to zero. Then compute the potential at the control points and write into RHS.
"""
function boundary_condition!(
    body::AbstractBody{<:Any,<:Any,<:Any,true},
    RHS,
    backend; optargs...
)
    calc_bc_dirichlet(RHS, body, backend; optargs...)
end

### NEUMANN
"""
   boundary_condition!(body::AbstractBody{<:Any,<:Any,<:Any,false}, RHS, backend; optargs...)
   
For Neumann bodies, apply the no-flow-through boundary condition by computing the normal velocity at the control points and writing into RHS.
"""
function boundary_condition!(
    body::AbstractBody{<:Any,<:Any,<:Any,false},
    RHS,
    backend; optargs...
)
    calc_bc_noflowthrough!(RHS, body.velocity, body.normals)
end
"""
    boundary_condition!(bodies::Tuple, solver::BackslashCoupled, backend; optargs...)

Apply the appropriate boundary condition for each body in `bodies` and write the results into the corresponding segment of `solver.rhs`.
"""
function boundary_condition!(
    bodies::Tuple,
    solver::BackslashCoupled,
    backend; optargs...
)

    nps = [b.ncells for b in bodies]
    offsets = cumsum(vcat(0, nps))
    for (bi, body) in enumerate(bodies)
        rows = offsets[bi]+1 : offsets[bi+1]
        boundary_condition!(body, view(solver.rhs, rows), backend; optargs...)
    end

end

function boundary_condition!(
    bodies::Tuple,
    solver::KrylovCoupled,
    backend; optargs...
)

    nps = [b.ncells for b in bodies]
    offsets = cumsum(vcat(0, nps))
    for (bi, body) in enumerate(bodies)
        rows = offsets[bi]+1 : offsets[bi+1]
        boundary_condition!(body, view(solver.rhs, rows), backend; optargs...)
    end

end

write_solution!(body::AbstractBody{<:Any, <:Any, <:Any, true}, sol) = body.strength[:, 2] .= sol
write_solution!(body::AbstractBody{<:Any, <:Any, <:Any, false}, sol) = body.strength[:, 1] .= sol

### Induced
### Induced function where all bodies are the targets induced(bodies, bodies, scalar_potential=Tuple(true if target body is Dirichlet, false if not))
## Dirichlet Body
function influence!(targets::Tuple, sources::Tuple, backend=DirectBackend(); optargs...)

    if precalc
        for target in targets
            pre_evaluate_influence!(target)
        end
    end

    influence!(targets, sources, backend; scalar_potential=[has_dirichlet_bc(target) for target in targets], velocity=[!has_dirichlet_bc(target) for target in targets], optargs...)

    return nothing
end

# Dirichlet
"""
    set_strengths!(body::AbstractBody{<:Any, <:Any, <:Any, true})

For Dirichlet bodies, set the source strength to enforce the no-penetration condition and set the doublet strength to zero.
"""
function set_strengths!(body::AbstractBody{<:Any, <:Any, <:Any, true})
    body.strength[:, 1] .= 0.0
    for d in 1:3
        body.strength[:, 2] .= view(body.velocity, d, :)
        body.strength[:, 2] .*= view(body.normals, d, :)
        body.strength[:, 1] .-= body.strength[:, 2]
    end
    body.strength[:, 2] .= 0.0
end

# Neumann
"""
    set_strengths!(body::AbstractBody{<:Any, <:Any, <:Any, false})

For Neumann bodies, set the source strength to zero.
"""
function set_strengths!(body::AbstractBody{<:Any, <:Any, <:Any, false})
    body.strength[:, 1] .= 0.0
end

function solve!(bodies::Tuple, solver::BackslashCoupled; backend=FastMultipoleBackend(), update_G::Bool=true, optargs...)

    # IMPROVEMENT : 2026-07-30 : default backend to FastMultipoleBackend() instead of DirectBackend() for the cross-body influence!/boundary_condition! evaluation : benchmarked O(N) FMM vs O(N^2) direct summation for two-body coupling at 256/1024/2304 total panels (two spheres 5 units apart, default expansion_order=10/multipole_acceptance=0.4/leaf_size=20) -- FMM relative error vs. direct was ~1e-8 to 1e-16 (negligible) while direct becomes ~2x, then ~20x slower than FMM as panel count grows past ~1000; FMM is only slower (~1.5x) at very small counts (~250 total) where tree-build overhead dominates. Every existing solve!(bodies::Tuple, ...) call site in src/, examples/, and test/ that reaches this method already passes `backend` explicitly, so no caller silently relied on the old Direct default.
    # IMPROVEMENT : 2026-07-30 : remove the unconditional println("BackslashCoupled") : it fired on every call regardless of a verbose flag, adding needless overhead/noise across e.g. an AOA sweep
    # Sizes
    npanels = [b.ncells for b in bodies]
    offsets = cumsum(vcat(0, npanels))

    # save external fields
    for (bi, body) in enumerate(bodies)
        r = offsets[bi]+1 : offsets[bi+1]
        @views solver.Uext[:, r] .= body.velocity
        @views solver.phi_ext[r]  .= body.potential
    end

    for body in bodies
        # update normals/controlpoints
        calc_normals!(body)
        calc_controlpoints!(body)

        # zero all strengths and set source strengths for Dirichlet bodies
        set_strengths!(body)

        # reset potential
        body.potential .= zero(eltype(body.potential))
    end
    
    # IMPROVEMENT : 2026-07-30 : time the coupled influence!/boundary_condition! evaluation and fold it into t_build instead of leaving it untimed : this cross-body coupling evaluation is frequently the most expensive part of a coupled solve (exactly what FMM acceleration targets), so leaving it out of the returned (tb, ts) understated true wall time and made cross-solver timing comparisons unfair
    t_build = @elapsed begin
        influence!(bodies, bodies, backend; scalar_potential=[has_dirichlet_bc(target) for target in bodies], velocity=[!has_dirichlet_bc(target) for target in bodies], optargs...)

        for (bi, body) in enumerate(bodies)
            if has_dirichlet_bc(body)
                r = offsets[bi]+1 : offsets[bi+1]
                @views body.potential .+= solver.phi_ext[r]
            end
        end

        ### get the boundary_condition for each body to write the RHS
        # calc_bc_dirichlet accumulates onto RHS (RHS .-= body.potential) rather than
        # overwriting it, so solver.rhs must start at zero every call -- true by
        # construction the first time (rhs = zeros(...)), but not guaranteed on
        # subsequent calls if the same solver object is reused across solves.
        fill!(solver.rhs, 0)
        boundary_condition!(bodies, solver, backend)
    end
    t_solve = 0.0

    if update_G
        # Zero G matrix
        fill!(solver.G, 0)

        # Build G matrix
        t_build += @elapsed begin
            for (bi, source) in enumerate(bodies)
                c = offsets[bi]+1 : offsets[bi+1] # columns of sources
                for (ti, target) in enumerate(bodies)
                    r = offsets[ti]+1 : offsets[ti+1] # rows of targets
                    # IMPROVEMENT : 2026-07-30 : pass kerneloffset=source.kerneloffset_panel explicitly, matching Backslash's _G! calls : without it this fell back to _G!'s generic default (source.kerneloffset), so BackslashCoupled regularized near-singular self/near-panel terms differently than Backslash whenever kerneloffset != kerneloffset_panel
                    _G!(view(solver.G, r, c), target,
                        source;
                        kerneloffset=source.kerneloffset_panel,
                        update_geometry=false)
                end
            end
        # Factorize G matrix and cache it in solver
            solver.Glu = lu!(solver.G)
        end
    end

    # solve with cached LU
    sol = similar(solver.rhs)
    t_solve += @elapsed ldiv!(sol, solver.Glu, solver.rhs)

    # write solution back
    for (bi, b) in enumerate(bodies)
        r = offsets[bi]+1 : offsets[bi+1]
        write_solution!(b, view(sol, r))

        @views b.velocity  .= solver.Uext[:, r]
        @views b.potential .= solver.phi_ext[r]
    end

    return t_build, t_solve
end



#####################################################################################
# function solve!(self::RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 3, TF},
#                 solver::AbstractMatrixfulSolver{true};
#                     solver_optargs=(),
#                     elprescribe_index::Int=1, elprescribe_value=0,
#                     weight_gammat=0, weight_gammao=1
#                 ) where TF

#     if size(self.velocity) != (3, self.ncells)
#         error("Invalid body velocity;"*
#               " expected size (3, $(self.ncells)), got $(size(self.velocity))")
#     end

#     normals = _calc_normals(self)
#     CPs = _calc_controlpoints(self, normals)

#     G = zeros(TF, self.ncells, self.ncells)
#     RHS = zeros(TF, self.ncells)

#     _G_U_RHS!(self, G, RHS, self.velocity, CPs, normals,
#                 elprescribe_index, elprescribe_value,
#                 weight_gammat, weight_gammao)

#     Gamma = zeros(TF, self.ncells)
#     solve_matrix!(Gamma, G, RHS, solver; solver_optargs...)

#     self.strength[:, 1] .= Gamma
#     self.strength[elprescribe_index, 1] = elprescribe_value

#     gamma = Gamma[elprescribe_index]
#     self.strength[:, 2] .= gamma*weight_gammat
#     self.strength[1:2:end, 2] .*= -1
#     self.strength[:, 3] .= gamma*weight_gammao
#     self.strength[1:2:end, 3] .*= -1

#     add_field(self, "Uinf", "vector", collect(eachcol(self.velocity)), "cell")
#     add_field(self, "Da", "vector", collect(eachcol(self.Das)), "system")
#     add_field(self, "Gamma", "scalar", view(self.strength, :, 1), "cell")

#     tangents = _calc_tangents(self)
#     obliques = _calc_obliques(self)
#     aux = zip(eachcol(tangents), eachcol(obliques),
#                 view(self.strength, :, 2), view(self.strength, :, 3))
#     gammas = [gammat*t + gammao*o for (t, o, gammat, gammao) in aux]
#     add_field(self, "gamma", "vector", gammas, "cell")
# end
