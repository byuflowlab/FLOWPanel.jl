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
recomputed before assembling `G`, with the sign of `target_system.CPoffset`
forced to match the formulation (negative/interior for Dirichlet, positive/
exterior for Neumann). The original `CPoffset` is restored on return.
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
        CPoffset_old = target_system.CPoffset
        target_system.CPoffset = abs(CPoffset_old) * (DBC ? -1 : 1)
        calc_normals!(target_system)
        calc_controlpoints!(target_system)
        target_system.CPoffset = CPoffset_old
    end

    kernel, strength_index = _G_kernel_and_strength_index(source_system)

    derivatives_switch = DBC ? FastMultipole.DerivativesSwitch(true, false, false) :
                               FastMultipole.DerivativesSwitch(false, true, false)

    # Save current strength and unit-activate the selected column
    old_strength = copy(source_system.strength)
    source_system.strength .= zero(eltype(source_system.strength))
    source_system.strength[:, strength_index] .= 1.0

    CPs = target_system.controlpoints
    normals = target_system.normals

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

    # Restore strength
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
"""
struct Backslash{TF,TGLU} <: AbstractMatrixfulSolver{false}
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

    _G!(G, body, body; kerneloffset=body.kerneloffset, update_geometry=true)
    Glu = lu!(G)

    return Backslash{TF,typeof(Glu)}(G, Glu, rhs, Uext, phi_ext)
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
        itmax::Int=20,         # Maximum number of iterations
        atol::Real=1e-6,            # Convergence tolerance
        rtol::Real=1e-6,            # Relative convergence tolerance
        backend::AbstractBackend=FastMultipoleBackend(),   # Backend to use
        preconditioner_cell_size::Real=0.0,  # cell size for block Jacobi preconditioner; ≤0 disables
    )
    TF = numtype(body)
    Uext = zeros(TF, 3, body.ncells)
    source_strengths = zeros(TF, body.ncells)
    unabbreviated_strengths = zeros(TF, body.ncells)
    normals = _calc_normals(body)

    # build block Jacobi preconditioner if requested
    if preconditioner_cell_size > 0
        preconditioner = FastMultipole.JacobiPreconditioner((body,); cell_size=preconditioner_cell_size)
    else
        preconditioner = nothing
    end

    return KrylovSolver{typeof(body), typeof(backend), TF, typeof(preconditioner)}(body, backend, Uext, normals, source_strengths, unabbreviated_strengths, method, itmax, Float64(atol), Float64(rtol), preconditioner)
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

function (solver::KrylovSolver{<:AbstractBody{<:Any, <:Any, <:Any, false}})(C, B, α, β)

    @assert length(B) == solver.body.ncells "Length of strengths vector does not match number of panels in body."
    solver.unabbreviated_strengths .= B
    _set_strength(solver.body, solver.unabbreviated_strengths)

    # get induced velocity at control points
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

function solve!(self::AbstractBody{<:Any,<:Any,<:Any,false}, solver::KrylovSolver{<:Any,B,TF}, Das=nothing; optargs...) where {B,TF}

    # save external velocity and update solver fields
    solver.Uext .= self.velocity
    solver.normals .= _calc_normals(self)
    calc_controlpoints!(self, solver.normals)

    # construct matrix-free linear operator
    TF2 = TF
    nrows = self.ncells
    ncols = self.ncells
    symmetric, hermitian = false, false
    # LinearOperators expects a callable (function) whose methods can be inspected.
    # Wrap the solver instance in a small closure so `methods` sees a function.
    prod! = (y, x, α, β) -> solver(y, x, α, β)
    A = LinearOperators.LinearOperator(
            TF2,
            nrows,
            ncols,
            symmetric, hermitian,
            prod!
        )

    # construct right-hand side
    RHS = zeros(TF2, nrows)
    calc_bc_noflowthrough!(RHS, self.velocity, solver.normals)

    # allocate and launch krylov solver
    workspace = Krylov.krylov_workspace(Val(solver.method), A, RHS)
    if solver.preconditioner !== nothing
        Krylov.krylov_solve!(workspace, A, RHS; M=solver.preconditioner, ldiv=true, atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
    else
        Krylov.krylov_solve!(workspace, A, RHS; atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
    end
    @show workspace.stats

    # store solution
    solver.unabbreviated_strengths .= workspace.x
    _set_strength(self, solver.unabbreviated_strengths)

    # restore external velocity
    self.velocity .= solver.Uext
end

function solve!(self::AbstractBody{<:Any,2,<:Any,true}, solver::KrylovSolver{<:Any,B,TF}, Das=nothing; update_G=false, optargs...) where {B,TF}

    solver.Uext .= self.velocity

    CPoffset_old = self.CPoffset
    self.CPoffset = -abs(CPoffset_old)

    solver.normals .= calc_normals!(self)
    calc_controlpoints!(self)

    _set_source_strength_from_velocity!(self, solver.Uext, solver.normals)
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
    self.potential .= 0
    influence!(self, self, solver.backend; scalar_potential=true, velocity=false)
    RHS .= -self.potential

    workspace = Krylov.krylov_workspace(Val(solver.method), A, RHS)
    if solver.preconditioner !== nothing
        ts = Krylov.krylov_solve!(workspace, A, RHS; M=solver.preconditioner, ldiv=true, atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
    else
        ts = Krylov.krylov_solve!(workspace, A, RHS; atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
    end
    @show workspace.stats

    solver.unabbreviated_strengths .= workspace.x
    self.strength[:, 1] .= solver.source_strengths
    _set_strength(self, solver.unabbreviated_strengths)

    self.CPoffset = CPoffset_old
    self.velocity .= solver.Uext
    self.potential .= 0
    return tb, ts
end


###############################################################################
# FGS Solver
################################################################################

"""
    FGSSolver(body; kwargs...)
    FGSSolver(bodies; kwargs...)

Block fixed-point / Gauss-Seidel-style solver for one or more coupled bodies.
"""
struct FGSSolver{TFGS,TF} <: AbstractMatrixFreeSolver
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
    Uext::Vector{Matrix{TF}}            # vector of per-body external velocity storage
    phi_ext::Vector{Vector{TF}}         # vector of per-body external potential storage
end

# Single-body convenience constructor
FGSSolver(body::AbstractBody; kwargs...) = FGSSolver((body,); kwargs...)

function FGSSolver(bodies::Tuple;
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
    )

    # save and set CPoffset to negative for Dirichlet bodies (interior solve)
    CPoffset_olds = zeros(Float64, length(bodies))
    for (i, body) in enumerate(bodies)
        CPoffset_olds[i] = body.CPoffset
        body.CPoffset = abs(body.CPoffset) * (-1)^has_dirichlet_bc(body)
    end

    # calculate control points if needed
    if calc_cps
        for body in bodies
            calc_normals!(body)
            calc_controlpoints!(body)
        end
    end

    # generate solver
    TF = promote_type(numtype.(bodies)...)
    fgs = FastMultipole.FastGaussSeidel(bodies; expansion_order, multipole_acceptance, leaf_size, shrink, recenter, extra_farfield=any(has_semiinfinite_wake.(bodies)))

    # restore CPoffsets
    for (body, CPoffset_old) in zip(bodies, CPoffset_olds)
        body.CPoffset = CPoffset_old
    end

    Uext = [zeros(TF, 3, body.ncells) for body in bodies]
    phi_ext = [zeros(TF, body.ncells) for body in bodies]
    return FGSSolver{typeof(fgs), TF}(fgs, Int(expansion_order), Int(leaf_size), Float64(multipole_acceptance), max_iterations, Int(inner_iterations), Float64(tolerance), Float64(rlx), Bool(reverse_pass), Bool(verbose), Uext, phi_ext)
end

#--- test solve! ---#

# Single-body convenience: wrap and delegate to tuple method
function solve!(self::AbstractBody, solver::FGSSolver; kwargs...)
    solve!((self,), solver; kwargs...)
end

################################################################################

function solve!(body::NonLiftingBody{TK,NK,TF,false}, solver::Backslash;
        backend=DirectBackend(), strength_index=1,
        update_G::Bool=false,
        optargs...
    ) where {TK, NK, TF}

    normals = calc_normals!(body)
    calc_controlpoints!(body)

    Glu = solver.Glu
    if update_G
        solver.G .= zero(eltype(solver.G))
        _G!(solver.G, body, body; optargs...)
        Glu = lu!(solver.G)
    end

    rhs = solver.rhs
    rhs .= zero(eltype(rhs))
    calc_bc_noflowthrough!(rhs, body.velocity, normals)

    ldiv!(view(body.strength, :, strength_index), Glu, rhs)

    return nothing
end

function solve!(self::NonLiftingBody{<:Union{ConstantSource, ConstantDoublet}, 2, TF, true},
                solver::Backslash; backend=DirectBackend(),
                update_G::Bool=false, optargs...) where TF

    solver.Uext .= self.velocity
    solver.phi_ext .= self.potential

    CPoffset_old = self.CPoffset
    self.CPoffset = -abs(CPoffset_old)

    calc_normals!(self)
    calc_controlpoints!(self)

    self.strength[:, 1] .= 0.0
    for d in 1:3
        self.strength[:, 2] .= view(self.velocity, d, :)
        self.strength[:, 2] .*= view(self.normals, d, :)
        self.strength[:, 1] .-= self.strength[:, 2]
    end
    self.strength[:, 2] .= 0.0

    self.potential .= 0
    influence!(self, self, backend; scalar_potential=true, velocity=false, optargs...)
    solver.rhs .= -self.potential

    Glu = solver.Glu
    if update_G
        solver.G .= 0.0
        _G!(solver.G, self, self; kerneloffset=self.kerneloffset)
        Glu = lu!(solver.G)
    end

    ldiv!(view(self.strength, :, 2), Glu, solver.rhs)

    self.CPoffset = CPoffset_old
    self.velocity .= solver.Uext
    self.potential .= solver.phi_ext

    return nothing
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

function solve!(self::RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing}, 2, TF}, solver::Backslash; backend=DirectBackend(), update_G=false, optargs...) where TF

    println("Backslash")
    solver.Uext .= self.velocity
    solver.phi_ext .= self.potential

    CPoffset_old = self.CPoffset
    self.CPoffset = -abs(CPoffset_old)

    calc_normals!(self)
    calc_controlpoints!(self)

    self.strength[:, 1] .= 0.0
    for d in (1,2,3)
        self.strength[:, 2] .= view(self.velocity, d, :)
        self.strength[:, 2] .*= view(self.normals, d, :)
        self.strength[:, 1] .-= self.strength[:, 2]
    end
    self.strength[:, 2] .= 0.0

    self.potential .= 0
    influence!(self, self, backend; scalar_potential=true, velocity=false, optargs...)
    solver.rhs .= self.potential
    solver.rhs .*= -1.0

    if update_G
        G = solver.G
        G .= 0.0
        t_build = @elapsed begin
            _G!(G, self, self; kerneloffset=self.kerneloffset)
            Glu = lu!(G)
        end
    else
        Glu = solver.Glu
        t_build = 0.0
    end

    t_solve = @elapsed begin
        ldiv!(view(self.strength, :, 2), Glu, solver.rhs)
    end

    self.CPoffset = CPoffset_old
    self.velocity .= solver.Uext
    self.potential .= solver.phi_ext
    return t_build, t_solve
end

function solve!(bodies::Tuple, solver::FGSSolver; backend = FastMultipoleBackend(
        expansion_order=solver.expansion_order,
        multipole_acceptance=solver.multipole_acceptance,
        leaf_size=solver.leaf_size
    ), optargs...)

    println("FGS")

    for (i, body) in enumerate(bodies)
        solver.Uext[i] .= body.velocity
        solver.phi_ext[i] .= body.potential
    end

    # Set CPoffset: negative (interior) for Dirichlet, positive (exterior) for Neumann
    CPoffset_olds = Tuple(body.CPoffset for body in bodies)
    for body in bodies
        body.CPoffset = abs(body.CPoffset) * (-1)^has_dirichlet_bc(body)
    end

    # update control points
    for body in bodies
        normals = calc_normals!(body)
        calc_controlpoints!(body, normals)
    end

    # store source strengths per body
    prior_sigmas = [copy(body.strength[:, 1]) for body in bodies]

    # induced potential due to source strengths for dirichlet bodies (interior solve)
    influence!(bodies, bodies, backend;
        scalar_potential=[has_dirichlet_bc(b) for b in bodies],
        velocity=[false for b in bodies],
        optargs...)

    # run solver
    FastMultipole.solve!(bodies, solver.fgs;
        max_iterations=solver.max_iterations,
        inner_iterations=solver.inner_iterations,
        tolerance=solver.tolerance,
        rlx=solver.rlx,
        scalar_potential=[has_dirichlet_bc(b) for b in bodies],
        gradient=[!has_dirichlet_bc(b) for b in bodies],
        hessian=false,
        reverse_pass=solver.reverse_pass,
        verbose=solver.verbose,
        final_update=false
    )

    # restore properties
    for (i, body) in enumerate(bodies)
        if has_dirichlet_bc(body)
            body.strength[:, 1] .= prior_sigmas[i]
        end
        body.CPoffset = CPoffset_olds[i]
        body.velocity .= solver.Uext[i]
        body.potential .= solver.phi_ext[i]
    end
end

function solve!(bodies::Tuple, solvers::Tuple;
    backend = fill(DirectBackend(), length(bodies)),
    max_outer_iterations::Int = 50,
    outer_tolerance::Real = 1e-8,
    verbose::Bool = false,
    optargs...)

    println("Tuple of bodies")

    N = length(bodies)
    @assert length(solvers) == N "Number of solvers ($(length(solvers))) must match number of bodies ($N)"

    prev_velocity = [copy(body.velocity) for body in bodies]
    prev_strengths = [copy(body.strength) for body in bodies]
    scalar_potential_flags = [has_dirichlet_bc(body) for body in bodies]
    velocity_flags = [!has_dirichlet_bc(body) for body in bodies]

    converged = false
    t_solve = 0.0
    t_build = 0.0
    for iter in 1:max_outer_iterations
        for (i, (body, solver)) in enumerate(zip(bodies, solvers))
            body.velocity .= prev_velocity[i]

            for (j, source) in enumerate(bodies)
                j == i && continue
                influence!(body, source, backend[j];
                    scalar_potential=scalar_potential_flags[i],
                    velocity=velocity_flags[i],
                    optargs...)
            end

            update_G = iter == 1
            tb, ts = solve!(body, solver; backend=backend[i], update_G, optargs...)
            t_build += tb
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

function solve!(body::NonLiftingBody{ConstantSource, 1, TF}, solver::FlatGroundSolver{TF}; optargs...) where TF
    normals = calc_normals!(body)
    calc_controlpoints!(body)

    calc_bc_noflowthrough!(solver.rhs, body.velocity, normals)

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
mutable struct BackslashCoupled{TF}
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
function set_strengths(body::AbstractBody{<:Any, <:Any, <:Any, true})
    body.strength[:, 1] .= 0.0
    for d in 1:3
        body.strength[:, 2] .= view(body.velocity, d, :)
        body.strength[:, 2] .*= view(body.normals, d, :)
        body.strength[:, 1] .-= body.strength[:, 2]
    end
    body.strength[:, 2] .= 0.0
    body.potential .= 0
end

# Neumann
"""
    set_strengths!(body::AbstractBody{<:Any, <:Any, <:Any, false})

For Neumann bodies, set the source strength to zero.
"""
function set_strengths(body::AbstractBody{<:Any, <:Any, <:Any, false})
    body.strength[:, 1] .= 0.0
end

function solve!(bodies::Tuple, solver::BackslashCoupled; backend=DirectBackend(), update_G::Bool=false, optargs...)

    println("BackslashCoupled")
    # Sizes
    npanels = [b.ncells for b in bodies]
    offsets = cumsum(vcat(0, npanels))

    # save external fields
    for (bi, body) in enumerate(bodies)
        r = offsets[bi]+1 : offsets[bi+1]
        @views solver.Uext[:, r] .= body.velocity
        @views solver.phi_ext[r]  .= body.potential
    end

    # flip CP offset
    CPoffset_old = map(b -> b.CPoffset, bodies)
    for b in bodies
        b.CPoffset = abs(b.CPoffset) * (-1)^(has_dirichlet_bc(b))
    end

    for body in bodies
        calc_normals!(body)
        calc_controlpoints!(body)

        ### Zero all strengths AND Set source strengths for Dirichlet bodies (function) set_strengths(body)
        set_strengths(body)
    end

    influence!(bodies, bodies, backend; scalar_potential=[has_dirichlet_bc(target) for target in bodies], velocity=[!has_dirichlet_bc(target) for target in bodies], optargs...)

    ### get the boundary_condition for each body to write the RHS
    boundary_condition!(bodies, solver, backend)

    if update_G
        # Zero G matrix
        fill!(solver.G, 0)
        
        # Build G matrix
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
        end

        # Factorize G matrix and cache it in solver
        solver.Glu = lu!(solver.G)
    end

    # solve with cached LU
    sol = similar(solver.rhs)
    t_solve = @elapsed ldiv!(sol, solver.Glu, solver.rhs)

    # write solution back
    for (bi, b) in enumerate(bodies)
        r = offsets[bi]+1 : offsets[bi+1]
        write_solution!(b, view(sol, r))

        b.CPoffset = CPoffset_old[bi]
        @views b.velocity  .= solver.Uext[:, r]
        @views b.potential .= solver.phi_ext[r]
    end

    return t_build, t_solve
end



#####################################################################################
# function solve!(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 1, TF},
#                 solver::AbstractMatrixfulSolver{true};
#                     solver_optargs=(),
#                     elprescribe::AbstractArray{Tuple{Int, Float64}}=[(1, 0.0)],
#                     GPUArray=Array{TF},
#                     update_G::Bool=true,
#                     optargs...
#                 ) where TF<:Real
#     if size(self.velocity) != (3, self.ncells)
#         error("Invalid body velocity;"*
#               " expected size (3, $(self.ncells)), got $(size(self.velocity))")
#     end

#     normals = _calc_normals(self)
#     CPs = _calc_controlpoints(self, normals)

#     G, RHS = _G_U_RHS(self, self.velocity, CPs, normals, elprescribe;
#                                                 GPUArray=GPUArray,
#                                                 optargs...)

#     Gamma = GPUArray(undef, self.ncells-length(elprescribe))
#     solve_matrix!(Gamma, G, RHS, solver; solver_optargs...)

#     if !(GPUArray <: Array)
#         Gamma = Array{TF}(Gamma)
#     end

#     set_solution(self, nothing, Gamma, elprescribe, self.velocity)
# end

# function _G_U_RHS(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 1}, args...; optargs...)
#     return _G_U_RHS_leastsquares(self, args...; optargs...)
# end

# function _G_U_RHS(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 2}, args...; optargs...)
#     @warn "_G_U_RHS called for RigidWakeBody{VortexRing, 2} as though `2` indicates the least-squares solver;
#     this is deprecated and may be removed in the future."
#     return _G_U_RHS_leastsquares(self, args...; optargs...)
# end

# function _G_U_RHS!(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 1}, args...; optargs...)
#     return _G_U_RHS_leastsquares!(self, args...; optargs...)
# end

# function _G_U_RHS!(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 2}, args...; optargs...)
#     @warn "_G_U_RHS! called for RigidWakeBody{VortexRing, 2} as though `2` indicates the least-squares solver;
#     this is deprecated and may be removed in the future."
#     return _G_U_RHS_leastsquares!(self, args...; optargs...)
# end

# function _G_U_RHS_leastsquares(self::AbstractBody,
#                                 Uinfs::AbstractMatrix{T1}, CPs, normals,
#                                 elprescribe::AbstractArray{Tuple{Int, T2}},
#                                 args...;
#                                 GPUArray=Array{promote_type(T1, T2)},
#                                 optargs...
#                                 ) where {T1, T2}

#     T = promote_type(T1, T2)

#     n = self.ncells
#     npres = length(elprescribe)

#     G = zeros(T, n, n)
#     Gred = zeros(T, n, n-npres)
#     tGred = zeros(T, n-npres, n)
#     gpuGred = GPUArray(undef, size(Gred))
#     Gls = GPUArray(undef, n-npres, n-npres)
#     RHS = zeros(T, n)
#     RHSls = GPUArray(undef, n-npres)

#     _G_U_RHS_leastsquares!(self, G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
#                 Uinfs, CPs, normals,
#                 elprescribe,
#                 args...; optargs...)

#     return Gls, RHSls
# end

# function _G_U_RHS_leastsquares!(self::AbstractBody,
#                                 G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
#                                 Uinfs, CPs, normals,
#                                 elprescribe::AbstractArray{Tuple{Int, T}};
#                                 onlycomputeG=false,
#                                 optargs...
#                                 ) where {T<:Number}

#     n = self.ncells
#     npres = length(elprescribe)

#     @assert size(G, 1)==n && size(G, 2)==n ""*
#         "Invalid $(size(G, 1))x$(size(G, 2)) matrix G; expected $(n)x$(n)"
#     @assert size(Gred, 1)==n && size(Gred, 2)==n-npres ""*
#         "Invalid $(size(Gred, 1))x$(size(Gred, 2)) matrix Gred; expected $(n)x$(n-npres)"
#     @assert size(tGred, 1)==n-npres && size(tGred, 2)==n ""*
#         "Invalid $(size(tGred, 1))x$(size(tGred, 2)) matrix tGred; expected $(n-npres)x$(n)"
#     @assert size(Gls, 1)==n-npres && size(Gls, 2)==n-npres ""*
#         "Invalid $(size(Gls, 1))x$(size(Gls, 2)) matrix Gls; expected $(n)x$(n-npres)"

#     @assert length(RHS)==n "Invalid RHS length $(length(RHS)); expected $(n)"
#     @assert length(RHSls)==n-npres "Invalid RHSls length $(length(RHSls)); expected $(n-pres)"

#     sort!(elprescribe, by = x -> x[1])

#     calc_bc_noflowthrough!(RHS, Uinfs, normals)

#     _G_U!(self, G, CPs, normals; optargs...)

#     if onlycomputeG
#         return Gls, RHSls
#     end

#     for (eli, elval) in elprescribe
#         for i in 1:length(RHS)
#             RHS[i] -= elval*G[i, eli]
#         end
#     end

#     prev_eli = 0
#     for (i, (eli, elval)) in enumerate(elprescribe)

#         Gred[:, (prev_eli+2-i):(eli-i)] .= view(G, :, (prev_eli+1):(eli-1))

#         if i==length(elprescribe) && eli!=size(G, 2)
#             Gred[:, (eli-i+1):end] .= view(G, :, eli+1:size(G, 2))
#         end

#         prev_eli = eli
#     end

#     if typeof(gpuGred) <: Array
#         permutedims!(tGred, Gred, [2, 1])
#         LA.mul!(RHSls, tGred, RHS)
#         LA.mul!(Gls, tGred, Gred)

#     else
#         copyto!(gpuGred, Gred)
#         tGred = transpose(gpuGred)
#         LA.mul!(RHSls, tGred, typeof(RHSls)(RHS))
#         LA.mul!(Gls, tGred, gpuGred)

#     end

#     return Gls, RHSls
# end

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
