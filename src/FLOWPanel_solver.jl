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

Implementations of `<:AbstractSolver` are expected to provide methods for `solve!()`
where the solver is passed as the third argument:
```julia
    function solve!(self::AbstractBody, solver::AbstractSolver; optargs...)
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

function solve!(self, solver; optargs...)
    throw(ErrorException("solve! not implemented for body of type $(typeof(self)) and solver of type $(typeof(solver))"))
end

################################################################################


################################################################################
# Backslash Operator
################################################################################

struct BackslashNeumann{TF,TGLU,LS} <: AbstractMatrixfulSolver{LS}
    G::Matrix{TF}    # Coefficient matrix
    Glu::TGLU
    rhs::Vector{TF}
end

function BackslashNeumann(body::AbstractBody{TK,1,TF}) where {TK,TF}
    
    # update control points (exterior)
    calc_normals!(body)
    calc_controlpoints!(body; off=abs(body.CPoffset))
    
    # populate G
    G = zeros(TF, body.ncells, body.ncells)
    _G_U!(body, TK, G, body.controlpoints, body.normals; kerneloffset=body.kerneloffset)

    # factorization
    Glu = lu!(G)

    # construct solver
    return BackslashNeumann{TF,typeof(Glu),false}(G, Glu, zeros(TF, body.ncells))
end

function Backslash(body::NonLiftingBody)
    return BackslashNeumann(body)
end

function Backslash(body::AbstractLiftingBody)
    return BackslashDirichlet(body)
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

#--- Dirichlet formulation ---#

struct BackslashDirichlet{TF,TGLU} <: AbstractMatrixfulSolver{false}
    G::Matrix{TF}
    Glu::TGLU
    rhs::Vector{TF}
    Uext::Matrix{TF}       # storage for external velocity (saved/restored around solve)
    phi_ext::Vector{TF}    # storage for external potential (saved/restored around solve)
end

function BackslashDirichlet(body::AbstractBody{<:Any,<:Any,TF}) where TF
    G = zeros(TF, body.ncells, body.ncells)
    rhs = zeros(TF, body.ncells)
    Uext = zeros(TF, 3, body.ncells)
    phi_ext = zeros(TF, body.ncells)

    # update control points (interior)
    calc_normals!(body)
    calc_controlpoints!(body; off=-abs(body.CPoffset))

    # populate G
    _G_phi!(body, ConstantDoublet, G, body.controlpoints; kerneloffset=body.kerneloffset)
    
    # factorization
    Glu = lu!(G)

    return BackslashDirichlet{TF,typeof(Glu)}(G, Glu, rhs, Uext, phi_ext)
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

function solve!(self::AbstractBody, solver::KrylovSolver{<:Any,B,TF}, Das=nothing; optargs...) where {B,TF}

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
    calc_bc_noflowthrough!(RHS, self.velocity, solver.normals)

    # allocate and launch krylov solver
    workspace = Krylov.krylov_workspace(Val(solver.method), A, RHS)
    Krylov.krylov_solve!(workspace, A, RHS; atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
    @show workspace.stats
    
    # store solution
    set_solution(self, solver.unabbreviated_strengths, workspace.x, solver.elprescribe, self.velocity)

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

calc_elprescribe(::NonLiftingBody{ConstantSource, 1}) = Tuple{Int,Float64}[]
calc_elprescribe(body::NonLiftingBody{VortexRing, 1}) = body.watertight ? [(1, 0.0)] : Tuple{Int,Float64}[]
calc_elprescribe(body::NonLiftingBody{ConstantDoublet, 1}) = body.watertight ? [(1, 0.0)] : Tuple{Int,Float64}[]

calc_elprescribe(::RigidWakeBody{ConstantSource, 1}) = Tuple{Int,Float64}[]
calc_elprescribe(body::RigidWakeBody{VortexRing, 1}) = body.watertight ? [(1, 0.0)] : Tuple{Int,Float64}[]
calc_elprescribe(body::RigidWakeBody{ConstantDoublet, 1}) = body.watertight ? [(1, 0.0)] : Tuple{Int,Float64}[]

function solve!(body::NonLiftingBody{TK,NK,TF}, solver::BackslashNeumann{<:Any, <:Any, false};
        backend=DirectBackend(), strength_index=1,
        update_G::Bool=false,
        optargs...
    ) where {TK, NK, TF}

    normals = calc_normals!(body)
    calc_controlpoints!(body)

    Glu = solver.Glu
    if update_G
        solver.G .= zero(eltype(solver.G))
        _G_U!(body, TK, solver.G, body.controlpoints, normals; optargs...)
        Glu = lu!(solver.G)
    end

    rhs = solver.rhs
    rhs .= zero(eltype(rhs))
    calc_bc_noflowthrough!(rhs, body.velocity, normals)

    ldiv!(view(body.strength, :, strength_index), Glu, rhs)

    return nothing
end

function solve!(self::RigidWakeBody{TK, 1, TF},
                    solver::AbstractMatrixfulSolver{false};
                    backend=FastMultipoleBackend(),
                    solver_optargs=(),
                    update_G::Bool=true,
                    optargs...
                ) where {TK<:Union{VortexRing, ConstantDoublet}, TF}
    if size(self.velocity) != (3, self.ncells)
        error("Invalid body velocity;"*
              " expected size (3, $(self.ncells)), got $(size(self.velocity))")
    end

    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    G = zeros(TF, self.ncells, self.ncells)
    _G_U!(self, TK, G, CPs, normals, backend; optargs...)

    RHS = calc_bc_noflowthrough(self.velocity, normals)

    Gamma = zeros(TF, self.ncells)
    solve_matrix!(Gamma, G, RHS, solver; solver_optargs...)

    self.strength[:, 1] .= Gamma

end

function solve!(self::RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing}, 2, TF}, solver::BackslashDirichlet; backend=DirectBackend(), update_G=false, optargs...) where TF

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
    influence!(self, self, backend; scalar_potential=true, optargs...)
    solver.rhs .= self.potential
    solver.rhs .*= -1.0

    if update_G
        G = solver.G
        G .= 0.0

        _G_phi!(self, ConstantDoublet, G, self.controlpoints; kerneloffset=self.kerneloffset)
        Glu = lu!(G)
    else
        Glu = solver.Glu
    end

    ldiv!(view(self.strength, :, 2), Glu, solver.rhs)

    self.CPoffset = CPoffset_old
    self.velocity .= solver.Uext
    self.potential .= solver.phi_ext
end

function solve!(bodies::Tuple, solver::FGSSolver; backend = FastMultipoleBackend(
        expansion_order=solver.expansion_order,
        multipole_acceptance=solver.multipole_acceptance,
        leaf_size=solver.leaf_size
    ), optargs...)

    for (i, body) in enumerate(bodies)
        solver.Uext[i] .= body.velocity
        solver.phi_ext[i] .= body.potential
    end

    # Set CPoffset: negative (interior) for Dirichlet, positive (exterior) for Neumann
    CPoffset_olds = Tuple(body.CPoffset for body in bodies)
    for body in bodies
        body.CPoffset = abs(body.CPoffset) * (-1)^has_dirichlet_bc(body)
    end

    for body in bodies
        normals = calc_normals!(body)
        calc_controlpoints!(body, normals)
    end

    # Prepare initial strengths per BC type
    sigmas = Tuple(begin
        if has_dirichlet_bc(body)
            # Dirichlet: compute sigma = -U·n, zero doublet
            body.strength[:, 1] .= 0.0
            for d in (1,2,3)
                body.strength[:, 2] .= view(body.velocity, d, :)
                body.strength[:, 2] .*= view(body.normals, d, :)
                body.strength[:, 1] .-= body.strength[:, 2]
            end
            sigma = copy(body.strength[:,1])
            body.strength[:, 2] .= 0.0
            sigma
        else
            # Neumann: zero all strengths
            body.strength .= 0.0
            nothing
        end
    end for body in bodies)

    # Evaluate initial influence from the initial strength distribution
    for body in bodies
        body.potential .= 0
        body.velocity .= 0
    end
    influence!(bodies, bodies, backend; scalar_potential=true, velocity=true, optargs...)

    # For Neumann bodies, the RHS is the freestream velocity at exterior CPs
    for (i, body) in enumerate(bodies)
        if !has_dirichlet_bc(body)
            body.velocity .+= solver.Uext[i]
        end
    end

    any_neumann = any(!has_dirichlet_bc(b) for b in bodies)
    any_dirichlet = any(has_dirichlet_bc(b) for b in bodies)

    FastMultipole.solve!(bodies, solver.fgs;
        max_iterations=solver.max_iterations,
        inner_iterations=solver.inner_iterations,
        tolerance=solver.tolerance,
        rlx=solver.rlx,
        scalar_potential=any_dirichlet,
        gradient=any_neumann,
        hessian=false,
        reverse_pass=solver.reverse_pass,
        verbose=solver.verbose,
        final_update=false
    )

    for (i, body) in enumerate(bodies)
        if has_dirichlet_bc(body)
            body.strength[:, 1] .= sigmas[i]
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

    N = length(bodies)
    @assert length(solvers) == N "Number of solvers ($(length(solvers))) must match number of bodies ($N)"

    Uinit = [copy(body.velocity) for body in bodies]
    prev_strengths = [copy(body.strength) for body in bodies]

    converged = false
    for iter in 1:max_outer_iterations

        for (i, (body, solver)) in enumerate(zip(bodies, solvers))
            body.velocity .= Uinit[i]

            for (j, source) in enumerate(bodies)
                j == i && continue
                influence!(body, source, backend[j];
                    scalar_potential=false,
                    velocity=true,
                    optargs...)
            end

            solve!(body, solver; backend=backend[i], optargs...)
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
        body.velocity .= Uinit[i]
    end

    return nothing
end

function solve!(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 1, TF},
                solver::AbstractMatrixfulSolver{true};
                    solver_optargs=(),
                    elprescribe::AbstractArray{Tuple{Int, Float64}}=[(1, 0.0)],
                    GPUArray=Array{TF},
                    update_G::Bool=true,
                    optargs...
                ) where TF<:Real
    if size(self.velocity) != (3, self.ncells)
        error("Invalid body velocity;"*
              " expected size (3, $(self.ncells)), got $(size(self.velocity))")
    end

    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    G, RHS = _G_U_RHS(self, self.velocity, CPs, normals, elprescribe;
                                                GPUArray=GPUArray,
                                                optargs...)

    Gamma = GPUArray(undef, self.ncells-length(elprescribe))
    solve_matrix!(Gamma, G, RHS, solver; solver_optargs...)

    if !(GPUArray <: Array)
        Gamma = Array{TF}(Gamma)
    end

    set_solution(self, nothing, Gamma, elprescribe, self.velocity)
end

function _G_U_RHS(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 1}, args...; optargs...)
    return _G_U_RHS_leastsquares(self, args...; optargs...)
end

function _G_U_RHS(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 2}, args...; optargs...)
    @warn "_G_U_RHS called for RigidWakeBody{VortexRing, 2} as though `2` indicates the least-squares solver;
    this is deprecated and may be removed in the future."
    return _G_U_RHS_leastsquares(self, args...; optargs...)
end

function _G_U_RHS!(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 1}, args...; optargs...)
    return _G_U_RHS_leastsquares!(self, args...; optargs...)
end

function _G_U_RHS!(self::RigidWakeBody{<:Union{VortexRing, ConstantDoublet}, 2}, args...; optargs...)
    @warn "_G_U_RHS! called for RigidWakeBody{VortexRing, 2} as though `2` indicates the least-squares solver;
    this is deprecated and may be removed in the future."
    return _G_U_RHS_leastsquares!(self, args...; optargs...)
end

function _G_U_RHS_leastsquares(self::AbstractBody,
                                Uinfs::AbstractMatrix{T1}, CPs, normals,
                                elprescribe::AbstractArray{Tuple{Int, T2}},
                                args...;
                                GPUArray=Array{promote_type(T1, T2)},
                                optargs...
                                ) where {T1, T2}

    T = promote_type(T1, T2)

    n = self.ncells
    npres = length(elprescribe)

    G = zeros(T, n, n)
    Gred = zeros(T, n, n-npres)
    tGred = zeros(T, n-npres, n)
    gpuGred = GPUArray(undef, size(Gred))
    Gls = GPUArray(undef, n-npres, n-npres)
    RHS = zeros(T, n)
    RHSls = GPUArray(undef, n-npres)

    _G_U_RHS_leastsquares!(self, G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
                Uinfs, CPs, normals,
                elprescribe,
                args...; optargs...)

    return Gls, RHSls
end

function _G_U_RHS_leastsquares!(self::AbstractBody,
                                G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
                                Uinfs, CPs, normals,
                                elprescribe::AbstractArray{Tuple{Int, T}};
                                onlycomputeG=false,
                                optargs...
                                ) where {T<:Number}

    n = self.ncells
    npres = length(elprescribe)

    @assert size(G, 1)==n && size(G, 2)==n ""*
        "Invalid $(size(G, 1))x$(size(G, 2)) matrix G; expected $(n)x$(n)"
    @assert size(Gred, 1)==n && size(Gred, 2)==n-npres ""*
        "Invalid $(size(Gred, 1))x$(size(Gred, 2)) matrix Gred; expected $(n)x$(n-npres)"
    @assert size(tGred, 1)==n-npres && size(tGred, 2)==n ""*
        "Invalid $(size(tGred, 1))x$(size(tGred, 2)) matrix tGred; expected $(n-npres)x$(n)"
    @assert size(Gls, 1)==n-npres && size(Gls, 2)==n-npres ""*
        "Invalid $(size(Gls, 1))x$(size(Gls, 2)) matrix Gls; expected $(n)x$(n-npres)"

    @assert length(RHS)==n "Invalid RHS length $(length(RHS)); expected $(n)"
    @assert length(RHSls)==n-npres "Invalid RHSls length $(length(RHSls)); expected $(n-pres)"

    sort!(elprescribe, by = x -> x[1])

    calc_bc_noflowthrough!(RHS, Uinfs, normals)

    _G_U!(self, G, CPs, normals; optargs...)

    if onlycomputeG
        return Gls, RHSls
    end

    for (eli, elval) in elprescribe
        for i in 1:length(RHS)
            RHS[i] -= elval*G[i, eli]
        end
    end

    prev_eli = 0
    for (i, (eli, elval)) in enumerate(elprescribe)

        Gred[:, (prev_eli+2-i):(eli-i)] .= view(G, :, (prev_eli+1):(eli-1))

        if i==length(elprescribe) && eli!=size(G, 2)
            Gred[:, (eli-i+1):end] .= view(G, :, eli+1:size(G, 2))
        end

        prev_eli = eli
    end

    if typeof(gpuGred) <: Array
        permutedims!(tGred, Gred, [2, 1])
        LA.mul!(RHSls, tGred, RHS)
        LA.mul!(Gls, tGred, Gred)

    else
        copyto!(gpuGred, Gred)
        tGred = transpose(gpuGred)
        LA.mul!(RHSls, tGred, typeof(RHSls)(RHS))
        LA.mul!(Gls, tGred, gpuGred)

    end

    return Gls, RHSls
end

function solve!(self::RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 3, TF},
                solver::AbstractMatrixfulSolver{true};
                    solver_optargs=(),
                    elprescribe_index::Int=1, elprescribe_value=0,
                    weight_gammat=0, weight_gammao=1
                ) where TF

    if size(self.velocity) != (3, self.ncells)
        error("Invalid body velocity;"*
              " expected size (3, $(self.ncells)), got $(size(self.velocity))")
    end

    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    G = zeros(TF, self.ncells, self.ncells)
    RHS = zeros(TF, self.ncells)

    _G_U_RHS!(self, G, RHS, self.velocity, CPs, normals,
                elprescribe_index, elprescribe_value,
                weight_gammat, weight_gammao)

    Gamma = zeros(TF, self.ncells)
    solve_matrix!(Gamma, G, RHS, solver; solver_optargs...)

    self.strength[:, 1] .= Gamma
    self.strength[elprescribe_index, 1] = elprescribe_value

    gamma = Gamma[elprescribe_index]
    self.strength[:, 2] .= gamma*weight_gammat
    self.strength[1:2:end, 2] .*= -1
    self.strength[:, 3] .= gamma*weight_gammao
    self.strength[1:2:end, 3] .*= -1

    add_field(self, "Uinf", "vector", collect(eachcol(self.velocity)), "cell")
    add_field(self, "Da", "vector", collect(eachcol(self.Das)), "system")
    add_field(self, "Gamma", "scalar", view(self.strength, :, 1), "cell")

    tangents = _calc_tangents(self)
    obliques = _calc_obliques(self)
    aux = zip(eachcol(tangents), eachcol(obliques),
                view(self.strength, :, 2), view(self.strength, :, 3))
    gammas = [gammat*t + gammao*o for (t, o, gammat, gammao) in aux]
    add_field(self, "gamma", "vector", gammas, "cell")
end
