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
Abstract type `AbstractSolver{MF<:Bool}` where `MF` indicates whether the solver 
explicitly forms a matrix (false) or uses a matrix-free approach (true).

Implementations of `<:AbstractSolver` are expected to provide methods for `solve()`
where the solver is passed as the third argument:
```julia
    function solve(self::AbstractBody, Uinfs::Array{<:Real, 2}, solver::AbstractSolver, args...; optargs...)
        .
        .
        .
    end
```
"""
abstract type AbstractSolver end

abstract type AbstractMatrixFreeSolver <: AbstractSolver end

"LS indicates whether the solver uses least-squares solution."
abstract type AbstractMatrixfulSolver{LS} <: AbstractSolver end

function solve2!(self, solver; optargs...)
    throw(ErrorException("solve2! not implemented for body of type $(typeof(self)) and solver of type $(typeof(solver))"))
end

################################################################################


################################################################################
# Backslash Operator
################################################################################

struct Backslash{TF,LS} <: AbstractMatrixfulSolver{LS}
    G::Matrix{TF}    # Coefficient matrix
end

function Backslash(self::AbstractBody{TK, 1, <:Any};
        backend::AbstractBackend=FastMultipoleBackend(),  # Backend to use for matrix construction
        kernel=TK,                      # Kernel type for matrix construction
        TFG=numtype(self), # Type for G matrix
        least_squares::Bool=true,     # Whether to use least squares solution
        optargs...                    # Additional optional arguments
    ) where TK

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    G = zeros(TFG, self.ncells, self.ncells)
    _G_U!(self, kernel, G, CPs, normals, backend; optargs...)

    return Backslash{TFG, least_squares}(G)
end

function numtype(self::AbstractBody)
    return promote_type(eltype(self.nodes),
                        eltype(self.strength),
                        Float64)
end

function get_strength_name(self::AbstractBody)
    return "strength"
end

get_strength_name(::AbstractBody{ConstantSource, 1, <:Any}) = "sigma"
get_strength_name(::AbstractBody{ConstantDoublet, 1, <:Any}) = "mu"
get_strength_name(::AbstractBody{VortexRing, 1, <:Any}) = "gamma"

# interface with existing methods
solve_matrix!(y, A, b, ::Backslash) = solve_backslash!(y, A, b)

#--- Dirichlet formulation ---#

struct BackslashDirichlet{TF} <: AbstractMatrixfulSolver{false}
    G::Matrix{TF}
    rhs::Vector{TF}
    Uext::Matrix{TF}       # storage for external velocity (saved/restored around solve)
    phi_ext::Vector{TF}    # storage for external potential (saved/restored around solve)
end

function BackslashDirichlet(body::AbstractBody{<:Any,<:Any,TF}) where TF
    G = zeros(TF, body.ncells, body.ncells)
    rhs = zeros(TF, body.ncells)
    Uext = zeros(TF, 3, body.ncells)
    phi_ext = zeros(TF, body.ncells)
    return BackslashDirichlet{TF}(G, rhs, Uext, phi_ext)
end

################################################################################
# GMRES Solver
################################################################################

struct KrylovSolver{TB<:AbstractBody,B<:AbstractBackend,TF<:Number} <: AbstractMatrixFreeSolver
    body::TB
    backend::B
    Uext::Array{TF, 2}    # storage for external velocity (saved/restored around solve)
    normals::Array{TF, 2} # Normals
    unabbreviated_strengths::Array{TF, 1} # Storage for unabbreviated strengths
    elprescribe::Vector{Tuple{Int,Float64}} # Prescribed element indices and values
    method::Symbol         # Krylov method to use
    itmax::Int           # Maximum number of iterations
    atol::Float64          # absolute tolerance
    rtol::Float64          # relative tolerance
end

function KrylovSolver(body::AbstractBody;
        method::Symbol=:gmres,    # Krylov method to use
        itmax::Int=20,         # Maximum number of iterations
        atol::Real=1e-6,            # Convergence tolerance
        rtol::Real=1e-6,            # Relative convergence tolerance
        backend::AbstractBackend=FastMultipoleBackend(),   # Backend to use
        elprescribe="automatic"      # Prescribed element indices and values
    )
    TF = numtype(body)
    Uext = zeros(TF, 3, body.ncells)
    unabbreviated_strengths = zeros(TF, body.ncells)
    normals = _calc_normals(body)
    elprescribe = elprescribe == "automatic" ? calc_elprescribe(body) : elprescribe
    return KrylovSolver{typeof(body), typeof(backend), TF}(body, backend, Uext, normals, unabbreviated_strengths, elprescribe, method, itmax, Float64(atol), Float64(rtol))
end

function _set_strength(body::AbstractBody, strengths, C, elprescribe=Tuple{Int,Float64}[])
    # check vector lengths
    @assert length(strengths) == body.ncells "Length of strengths vector does not match number of panels in body."
    @assert length(C) + length(elprescribe) == body.ncells "Length of abbreviated strengths vector plus number of prescribed strengths does not match number of panels in body."

    # populate unabbreviated strengths
    Ui = 1
    Ci = 1
    for (i, val) in elprescribe
        strengths[i] = val
        rng = Ui:i-1
        Crng = Ci:Ci+length(rng)-1
        if length(rng) > 0
            strengths[rng] .= view(C, Crng)
        end
        Ui = i + 1
        Ci += length(rng)
    end

    # fill in remaining strengths after last prescribed element
    if Ui <= body.ncells
        rng = Ui:body.ncells
        Crng = Ci:Ci+length(rng)-1
        strengths[rng] .= view(C, Crng)
        @assert Ci + length(rng) - 1 == length(C) "Length of abbreviated strengths vector does not match number of non-prescribed panels in body."
    end

    # set strengths in body
    _set_strength(body, strengths)
end

function _set_strength(body::AbstractBody{<:Any, 1, <:Any}, strengths)
    body.strength[:, 1] .= strengths
end

function (solver::KrylovSolver)(C, B, α, β)

    # set strengths in body
    # NOTE: C is an abbreviated vector if prescribed strengths are used,
    #       effectively skipping the prescribed strengths;
    #       in that case, we need to set the strengths accordingly
    _set_strength(solver.body, solver.unabbreviated_strengths, B, solver.elprescribe)

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

function solve2!(self::AbstractBody, solver::KrylovSolver{<:Any,B,TF}, Das=nothing; optargs...) where {B,TF}

    # save external velocity and update solver fields
    solver.Uext .= self.velocity
    solver.normals .= _calc_normals(self)
    calc_controlpoints!(self, solver.normals)

    # construct matrix-free linear operator
    TF2 = TF
    nrows = self.ncells
    ncols = self.ncells - length(solver.elprescribe)
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

    # verify solver compatibility
    if solver.method == :gmres
        @assert nrows == ncols "GMRES solver requires a square matrix; got $(nrows)x$(ncols)."
    end

    # construct right-hand side
    RHS = zeros(TF2, nrows)
    calc_bc_noflowthrough!(RHS, Uinfs, solver.normals)

    # allocate and launch krylov solver
    workspace = Krylov.krylov_workspace(Val(solver.method), A, RHS)
    Krylov.krylov_solve!(workspace, A, RHS; atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
    @show workspace.stats
    
    # store solution
    set_solution(self, solver.unabbreviated_strengths, workspace.x, solver.elprescribe, Uinfs)

    # restore external velocity
    self.velocity .= solver.Uext
end

###############################################################################
# FGS Solver
################################################################################

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

    # save and set CPoffset to negative for each body (interior solve)
    CPoffset_olds = Tuple(body.CPoffset for body in bodies)
    for body in bodies
        body.CPoffset = -abs(body.CPoffset)
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
    println("building FastMultipole.fgs:")
    @time fgs = FastMultipole.FastGaussSeidel(bodies; expansion_order, multipole_acceptance, leaf_size, shrink, recenter, extra_farfield=any(has_semiinfinite_wake.(bodies)))
    println("done.")

    # restore CPoffsets
    for (body, CPoffset_old) in zip(bodies, CPoffset_olds)
        body.CPoffset = CPoffset_old
    end

    Uext = [zeros(TF, 3, body.ncells) for body in bodies]
    phi_ext = [zeros(TF, body.ncells) for body in bodies]
    return FGSSolver{typeof(fgs), TF}(fgs, Int(expansion_order), Int(leaf_size), Float64(multipole_acceptance), max_iterations, Int(inner_iterations), Float64(tolerance), Float64(rlx), Bool(reverse_pass), Bool(verbose), Uext, phi_ext)
end

#--- test solve! ---#

function solve2!(self::AbstractBody{<:Any,1,<:Any}, solver::FGSSolver; optargs...)

    # solve system
    FastMultipole.solve!((self,), solver.fgs; max_iterations=solver.max_iterations, inner_iterations=solver.inner_iterations, tolerance=solver.tolerance, rlx=solver.rlx, verbose=solver.verbose)

    # store solution
    set_solution(self, self.strength, self.strength, Tuple{Int,Float64}[], Uinfs)
end

################################################################################