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

function solve2!(self::AbstractBody, Uinfs::Array{<:Real, 2}, solver::AbstractSolver, args...; optargs...)
    throw(ErrorException("solve2! not implemented for body of type $(typeof(self)) and solver of type $(typeof(solver))"))
end

################################################################################


################################################################################
# Backslash Operator
################################################################################

struct Backslash{TF,LS} <: AbstractMatrixfulSolver{LS}
    G::Matrix{TF}    # Coefficient matrix
end

function Backslash(self::AbstractBody;
        TFG=numtype(self), # Type for G matrix
        least_squares::Bool=true,     # Whether to use least squares solution
        optargs...                    # Additional optional arguments
    )

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    G = zeros(TFG, self.ncells, self.ncells)
    _G_U!(self, G, CPs, normals; optargs...)

    return Backslash{TFG, least_squares}(G)
end

function numtype(self::AbstractBody)
    return promote_type(eltype(self.grid._nodes),
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


################################################################################
# LU Decomposition Solver
################################################################################

# struct LUDecomposition{TF,TP,LS<:Bool} <: AbstractSolver{false}
#     G::Matrix{TF}    # Coefficient matrix
#     store_LU::Bool   # Whether to store the LU decomposition
#     ALU::Matrix{TP}  # LU decomposition of G (using primal type TP)
# end

# function LUDecomposition(self::AbstractBody; 
#         store_LU::Bool=true,         # Whether to store the LU decomposition
#         TFG=eltype(self.grid._nodes) # Type for G matrix
#     )

#     # Compute normals and control points
#     normals = _calc_normals(self)
#     CPs = _calc_controlpoints(self, normals)

#     # Compute geometric matrix (left-hand-side influence matrix)
#     G = zeros(TFG, self.ncells, self.ncells)
#     _G_U!(self, G, CPs, normals, Das, Dbs; optargs...)
# end

# solve_matrix!(y, A, b, ::LU; Alu=nothing) = solve_ludiv!(y, A, b; Alu)

################################################################################
# GMRES Solver
################################################################################

struct KrylovSolver{TB<:AbstractBody,B<:AbstractBackend,TF<:Number} <: AbstractMatrixFreeSolver
    body::TB
    backend::B
    Uind::Array{TF, 2}    # induced velocity storage
    CPs::Array{TF, 2}     # Control points
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
        # restart::Int=50,           # Number of iterations between restarts
        backend::AbstractBackend=FastMultipoleBackend(),   # Backend to use
        elprescribe="automatic"      # Prescribed element indices and values
    )
    TF = numtype(body)
    Uind = zeros(TF, 3, body.ncells)
    unabbreviated_strengths = zeros(TF, body.ncells)
    normals = _calc_normals(body)
    CPs = _calc_controlpoints(body, normals)
    elprescribe = elprescribe == "automatic" ? calc_elprescribe(body) : elprescribe
    return KrylovSolver{typeof(body), typeof(backend), TF}(body, backend, Uind, CPs, normals, unabbreviated_strengths, elprescribe, method, itmax, Float64(atol), Float64(rtol))#, restart)
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
    solver.Uind .= zero(solver.Uind) # reset induced velocity
    _Uind!(solver.body, solver.CPs, solver.Uind, solver.backend)

    # dot product with normals
    solver.Uind .*= solver.normals
    solver.Uind[1,:] .+= view(solver.Uind, 2, :)
    solver.Uind[1,:] .+= view(solver.Uind, 3, :)
    
    # scale and add
    C .*= β
    C .+= α .* view(solver.Uind, 1, :)
end

function solve2!(self::AbstractBody, Uinfs::Array{<:Real, 2}, solver::KrylovSolver{<:Any,B,TF}, Das=nothing, Dbs=nothing; optargs...) where {B,TF}
    
    # update solver fields
    solver.normals .= _calc_normals(self)
    solver.CPs .= _calc_controlpoints(self, solver.normals) # TODO: avoid allocating here?

    # construct matrix-free linear operator
    TF2 = promote_type(eltype(Uinfs), TF)
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
end

# solve_matrix!(y, A, b, ::GMRES; Avalue=nothing, optargs...) =
#     solve_gmres!(y, A, b; Avalue=Avalue, optargs...)

###############################################################################
# FGS Solver
################################################################################

struct FGSSolver{TFGS,TF<:Number} <: AbstractMatrixFreeSolver
    fgs::TFGS
    max_iterations::Int
    tolerance::Float64
end

function FGSSolver(body::AbstractBody; 
        max_iterations::Int=100,         # Maximum number of iterations
        tolerance::Real=1e-6,            # Convergence tolerance
        expansion_order=7,
        multipole_acceptance=0.4,
        leaf_size=10,
        shrink=false,
        recenter=false
    )

    # generate solver
    TF = numtype(body)
    fgs = FastMultipole.FastGaussSeidel((body,), (body,); expansion_order, multipole_acceptance, leaf_size, shrink, recenter)

    return FGSSolver{typeof(fgs), TF}(fgs, max_iterations, Float64(tolerance))
end

#--- test solve! ---#

function solve2!(self::AbstractBody, Uinfs::Array{<:Real, 2}, solver::FGSSolver{<:Any,TF}; optargs...) where {TF}
    
    # construct right-hand side
    # TF2 = promote_type(eltype(Uinfs), TF)
    # RHS = zeros(TF2, self.ncells)
    # normals = _calc_normals(self)
    # calc_bc_noflowthrough!(RHS, Uinfs, normals)

    # apply freestream
    self.velocity .= Uinfs

    # solve system
    FastMultipole.solve!(self, solver.fgs; max_iterations=10, tolerance=1e-3)

    # store solution
    set_solution(self, self.strength, self.strength, Tuple{Int,Float64}[], Uinfs)
end

################################################################################