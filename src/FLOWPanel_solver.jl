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
    # assignment, not accumulation: the coupled solvers reuse a persistent rhs
    # across solves, and `.-=` doubled the Dirichlet rows on every solve after
    # the first (021 Phase 1 finding, 2026-08-14)
    RHS .= .-self.potential
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
    _G!(G, target_system, source_system; core_size=source_system.core_size, update_geometry=false)

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
             core_size=source_system.core_size,
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

            phi, u, _ = induced(target, source_system, i_source, derivatives_switch; core_size)

            if DBC
                isnan(phi) && error("NaN encountered in G matrix computation: \ni_source = $i_source, i_target = $i_target, target = $target, source_strength = $(source_system.strength[i_source, strength_index]), kernel = $kernel, core_size = $core_size")
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
"""
mutable struct Backslash{TF,TGLU} <: AbstractMatrixfulSolver{false}
    G::Matrix{TF}
    Glu::TGLU              # aliases G's memory (lu!); refreshed in place on update_G=true
    rhs::Vector{TF}
    Uext::Matrix{TF}       # storage for external velocity (saved/restored around solve)
    phi_ext::Vector{TF}    # storage for external potential (saved/restored around solve)
end

function Backslash(body::AbstractBody{<:Any,<:Any,TF}) where TF
    G = zeros(TF, body.ncells, body.ncells)
    rhs = zeros(TF, body.ncells)
    Uext = zeros(TF, 3, body.ncells)
    phi_ext = zeros(TF, body.ncells)

    calc_normals!(body)
    calc_controlpoints!(body)
    _G!(G, body, body; core_size=body.core_size_panel, update_geometry=false)
    Glu = lu!(G)

    return Backslash{TF,typeof(Glu)}(G, Glu, rhs, Uext, phi_ext)
end


function solve!(body::AbstractBody{<:Any,<:Any,<:Any,true}, solver::AbstractSolver;
        backend=DirectBackend(),
        optargs...)

    potential_old = copy(body.potential)

    try
        set_strengths!(body)
        # For single-body Dirichlet solves, `body.potential` is workspace for
        # the interior self/source potential and the solve is homogeneous unless
        # external influences have already been assembled into that workspace.
        body.potential .= zero(eltype(body.potential))
        influence!(body, body, backend; scalar_potential=true, velocity=false, optargs...)
        _solve!(body, solver; backend, optargs...)
    finally
        body.potential .= potential_old
    end

    return nothing
end

function solve!(body::AbstractBody{<:Any,<:Any,<:Any,false}, solver::AbstractSolver;
        backend=DirectBackend(),
        optargs...)

    if body isa RigidWakeBody && body.watertight
        @warn "Solving a watertight RigidWakeBody with the Neumann formulation " *
              "(DBC=false) gives a rank-deficient influence matrix; results " *
              "will be unreliable. Use DBC=true (Dirichlet) for closed " *
              "surfaces, or remove a cap to make the surface non-watertight." maxlog=1
    end

    _solve!(body, solver; backend, optargs...)

    return nothing
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
    KrylovOperator(body, backend)

Matrix-free application of the body's self-influence operator (`apply_G!`'s
code path) as the 5-argument mul! callable LinearOperators.jl expects.
"""
struct KrylovOperator{TB<:AbstractBody,B<:AbstractBackend,TF}
    body::TB
    backend::B
    normals::Matrix{TF}            # snapshot refreshed at the start of each solve
    strengths_scratch::Vector{TF}
    cache_tree::Bool               # reuse the FMM plan across applies within a solve
    cache_nearfield::Bool          # dense near-field cache on the plan (implies cache_tree)
    nearfield_cache_max_bytes::Int # size cap for that cache (FastMultipole default is 4 GiB)
    nearfield_cache_max_build_time::Float64  # wall-clock cap; the cache build is serial
    plan_slot::Base.RefValue{Any}  # (FmmPlan, key) or nothing; cleared around each solve
end

function KrylovOperator(body::AbstractBody, backend::AbstractBackend;
                        cache_tree::Bool=false, cache_nearfield::Bool=false,
                        nearfield_cache_max_bytes::Integer=FastMultipole.NEARFIELD_CACHE_DEFAULT_MAX_BYTES,
                        nearfield_cache_max_build_time::Real=Inf)
    TF = numtype(body)
    return KrylovOperator{typeof(body),typeof(backend),TF}(
        body, backend, copy(body.normals), zeros(TF, body.ncells),
        cache_tree, cache_nearfield, Int(nearfield_cache_max_bytes),
        Float64(nearfield_cache_max_build_time), Ref{Any}(nothing))
end

function (op::KrylovOperator{<:AbstractBody{<:Any, <:Any, <:Any, false}})(C, B, α, β)
    @assert length(B) == op.body.ncells "Length of strengths vector does not match number of panels in body."
    Gx = _apply_neumann_G!(op.body, B, op.backend, op.normals, op.strengths_scratch;
                           plan_slot=op.cache_tree ? op.plan_slot : nothing,
                           cache_nearfield=op.cache_nearfield,
                           nearfield_cache_max_bytes=op.nearfield_cache_max_bytes,
                           nearfield_cache_max_build_time=op.nearfield_cache_max_build_time)
    C .*= β
    C .+= α .* Gx
end

function (op::KrylovOperator{<:AbstractBody{<:Any, <:Any, <:Any, true}})(C, B, α, β)
    @assert length(B) == op.body.ncells "Length of strengths vector does not match number of panels in body."
    Gx = _apply_dirichlet_G!(op.body, B, op.backend, op.strengths_scratch;
                             plan_slot=op.cache_tree ? op.plan_slot : nothing,
                             cache_nearfield=op.cache_nearfield,
                             nearfield_cache_max_bytes=op.nearfield_cache_max_bytes,
                             nearfield_cache_max_build_time=op.nearfield_cache_max_build_time)
    C .*= β
    C .+= α .* Gx
end

# Krylov methods whose workspace takes a `memory` argument
const _KRYLOV_MEMORY_METHODS = (:gmres, :fgmres, :fom, :dqgmres, :diom)

"""
    KrylovSolver(body; method=:gmres, itmax=20, atol=1e-6, rtol=1e-6,
                 backend=FastMultipoleBackend(), memory=20,
                 preconditioner=nothing, preconditioner_cell_size=0.0,
                 warmstart=false, record_history=false)

Matrix-free iterative solver that evaluates self-influence through the selected
backend. The linear operator, right-hand-side vector, and Krylov workspace are
built once at construction and reused across solves.

Options:
- `method`: any Krylov.jl method symbol (`:gmres`, `:fgmres`, ...). Use
  `:fgmres` with an [`FGSPreconditioner`](@ref) (flexible right
  preconditioning).
- `memory`: restart/memory parameter for methods that take one.
- `preconditioner`: `nothing`, a `FastMultipole.JacobiPreconditioner`, or an
  [`FGSPreconditioner`](@ref). Both are applied as *right* preconditioners
  (`N=`, 021 ruling 11a) so the monitored/stopping residual is the true one.
  May be swapped between solves.
- `preconditioner_cell_size > 0` builds a block-Jacobi preconditioner
  (back-compatible shorthand, ignored when `preconditioner` is given).
- `warmstart`: seed each solve with the previous solution (`x0`); off by
  default. See BRAINSTORM 021 Phase 0 W3.
- `record_history`: capture per-iteration residuals + wall-clock timestamps
  into `solver.history` (021 ruling 4).
- `cache_tree`: build the FMM trees/interaction lists once per SOLVE and
  reuse them across the solve's operator applies (021 Phase 2b); off by
  default. Scope is strictly per-solve: the plan is built on the first apply
  and dropped when the solve returns, so geometry and the active
  `core_size` (which feeds the FMM buffer radii via `radius_inflation`)
  are frozen for the plan's whole lifetime by construction. Only meaningful
  with a `FastMultipoleBackend`.
- `cache_nearfield`: additionally freeze the plan's near field as dense
  influence blocks (`FastMultipole.NearfieldInfluenceCache`) built lazily on
  the solve's first operator apply alongside the plan, so every subsequent
  apply evaluates the near field as packed BLAS matvecs instead of analytic
  kernel calls (021 Phase 2b). Implies (and requires) `cache_tree=true`; the
  cache shares the plan's lifetime and validity contract. Solution shifts at
  the rtol level are expected (BLAS summation order differs from the
  kernel's).
- `persistent_plan`: keep the plan (and any near-field cache on it) alive
  ACROSS solves instead of dropping it when each solve returns (the
  originally deferred nearfield-cache commit 4, extended with rigid-motion
  support). Implies (and requires) `cache_tree=true`. Contract: between
  solves the body geometry either stays frozen or moves RIGIDLY with every
  rigid step mirrored by [`transform_solver_geometry!`](@ref) (`simulate!`
  does this automatically from the frame kinematics), and the active
  core_size that governed the plan's buffer radii is restored before
  every solve (the production Dirichlet path does this — the plan is built
  lazily inside `_krylov_launch!` under the panel-offset state, which
  `_set_core_sizes!` re-establishes each solve). The Dirichlet
  scalar-potential operator (and its near-field cache blocks) is EXACTLY
  invariant under rigid motion, so build cost amortizes over the run.

After each solve, `solver.niter` and `solver.solved` report Krylov's stats.
"""
mutable struct KrylovSolver{TB<:AbstractBody,B<:AbstractBackend,TF<:Number,TP,TK,TA,TW} <: AbstractMatrixFreeSolver
    body::TB
    backend::B
    Uext::Array{TF, 2}    # storage for external velocity (saved/restored around solve)
    source_strengths::Array{TF, 1} # Fixed source strengths for Dirichlet solves
    unabbreviated_strengths::Array{TF, 1} # Storage for solution strengths
    kop::TK                # KrylovOperator (owns the normals snapshot + scratch)
    A::TA                  # LinearOperators wrapper around kop (built once)
    rhs::Vector{TF}        # preallocated right-hand side
    workspace::TW          # Krylov workspace (built once, reused)
    method::Symbol         # Krylov method to use
    itmax::Int             # Maximum number of iterations
    atol::Float64          # absolute tolerance
    rtol::Float64          # relative tolerance
    memory::Int            # restart/memory parameter (methods that take one)
    preconditioner::TP     # nothing, JacobiPreconditioner, or FGSPreconditioner
    warmstart::Bool        # seed solves with x_prev (off by default)
    x_prev::Vector{TF}     # previous solution (warmstart seed)
    have_x_prev::Bool      # whether x_prev holds a valid solution
    record_history::Bool   # capture per-iteration convergence history
    history::ConvergenceHistory
    niter::Int             # iterations of the last solve
    solved::Bool           # convergence flag of the last solve
    cache_tree::Bool       # per-solve FMM plan reuse (see docstring)
    cache_nearfield::Bool  # dense near-field cache on the per-solve plan (see docstring)
    persistent_plan::Bool  # cross-solve plan (+cache) persistence (see docstring)
end

function KrylovSolver(body::AbstractBody;
        method::Symbol=:gmres,    # Krylov method to use
        itmax::Int=20,         # Maximum number of iterations
        atol::Real=1e-6,            # Convergence tolerance
        rtol::Real=1e-6,            # Relative convergence tolerance
        backend::AbstractBackend=FastMultipoleBackend(),   # Backend to use
        memory::Int=20,             # restart/memory parameter
        preconditioner=nothing,     # preconditioner object (see docstring)
        preconditioner_cell_size::Real=0.0,  # cell size for block Jacobi preconditioner; ≤0 disables
        warmstart::Bool=false,      # seed solves with the previous solution
        record_history::Bool=false, # capture per-iteration convergence history
        cache_tree::Bool=false,     # per-solve FMM plan reuse (021 Phase 2b)
        cache_nearfield::Bool=false, # dense near-field cache on the plan (021 Phase 2b)
        # Caps for that cache. FastMultipole's 4 GiB default is below what the
        # larger 021 rungs need (R4 ≈ 4.5 GiB), and the build is serial, so the
        # time cap is wall-clock rather than a per-thread budget.
        nearfield_cache_max_bytes::Integer=FastMultipole.NEARFIELD_CACHE_DEFAULT_MAX_BYTES,
        nearfield_cache_max_build_time::Real=Inf,
        persistent_plan::Bool=false, # cross-solve plan persistence (see docstring)
    )
    TF = numtype(body)
    cache_tree |= cache_nearfield | persistent_plan  # both live on the plan, so they imply cache_tree
    cache_tree && !(backend isa FastMultipoleBackend) && throw(ArgumentError(
        "cache_tree=true (or cache_nearfield/persistent_plan=true) requires a FastMultipoleBackend (got $(typeof(backend)))"))
    Uext = zeros(TF, 3, body.ncells)
    source_strengths = zeros(TF, body.ncells)
    unabbreviated_strengths = zeros(TF, body.ncells)
    calc_normals!(body)
    calc_controlpoints!(body)

    # build block Jacobi preconditioner if requested and no object was given
    if preconditioner === nothing && preconditioner_cell_size > 0
        preconditioner = FastMultipole.JacobiPreconditioner((body,); cell_size=preconditioner_cell_size)
    end

    # matrix-free operator + persistent workspace
    kop = KrylovOperator(body, backend; cache_tree, cache_nearfield,
                         nearfield_cache_max_bytes, nearfield_cache_max_build_time)
    n = body.ncells
    prod! = (y, x, α, β) -> kop(y, x, α, β)
    A = LinearOperators.LinearOperator(TF, n, n, false, false, prod!)
    rhs = zeros(TF, n)
    workspace = method in _KRYLOV_MEMORY_METHODS ?
        Krylov.krylov_workspace(Val(method), A, rhs; memory) :
        Krylov.krylov_workspace(Val(method), A, rhs)

    x_prev = zeros(TF, n)
    history = ConvergenceHistory(:krylov_precnorm)

    return KrylovSolver{typeof(body), typeof(backend), TF, typeof(preconditioner),
                        typeof(kop), typeof(A), typeof(workspace)}(
        body, backend, Uext, source_strengths, unabbreviated_strengths,
        kop, A, rhs, workspace, method, itmax, Float64(atol), Float64(rtol),
        Int(memory), preconditioner, warmstart, x_prev, false,
        record_history, history, 0, false, cache_tree, cache_nearfield,
        persistent_plan)
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

# Launch the (persistent) Krylov workspace on solver.rhs: preconditioner
# routing, optional warmstart, optional history capture. Assumes solver.rhs
# holds the right-hand side; leaves the solution in solver.workspace.x and
# updates solver.niter / solver.solved.
function _krylov_launch!(solver::KrylovSolver)
    ws = solver.workspace
    A = solver.A
    rhs = solver.rhs

    # cache_tree scope = exactly this launch: clear any stale plan so the
    # first apply rebuilds under the CURRENT geometry + active core_size
    # (the buffer radii depend on it via radius_inflation). persistent_plan
    # opts out: the caller guarantees frozen-or-rigidly-mirrored geometry and
    # restored core_size state across solves (see the KrylovSolver
    # docstring), so the surviving plan stays valid.
    solver.cache_tree && !solver.persistent_plan && (solver.kop.plan_slot[] = nothing)

    # convergence-history capture (021 ruling 4): iteration + timestamp +
    # Krylov's internal residual (preconditioned norm when a left
    # preconditioner is active)
    if solver.record_history
        reset!(solver.history; metric=:krylov_precnorm)
        callback = ws -> begin
            # ws.stats.niter is only finalized after the solve; count locally
            res = isempty(ws.stats.residuals) ? NaN : ws.stats.residuals[end]
            record!(solver.history, length(solver.history) + 1, res)
            false
        end
    else
        callback = ws -> false
    end

    common = (; atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax,
              history=solver.record_history, callback)
    use_x0 = solver.warmstart && solver.have_x_prev

    P = solver.preconditioner
    if P === nothing
        use_x0 ? Krylov.krylov_solve!(ws, A, rhs, solver.x_prev; common...) :
                 Krylov.krylov_solve!(ws, A, rhs; common...)
    elseif P isa FGSPreconditioner
        # flexible right preconditioning: the residual Krylov monitors stays
        # the true one, and FGMRES tolerates a per-iteration-varying apply
        use_x0 ? Krylov.krylov_solve!(ws, A, rhs, solver.x_prev; N=P, ldiv=true, common...) :
                 Krylov.krylov_solve!(ws, A, rhs; N=P, ldiv=true, common...)
    else
        # right preconditioning for the Jacobi path too (021 ruling 11a): the
        # monitored/stopping residual stays the TRUE residual — left routing
        # stops on a preconditioned norm, which reported "converged" at a much
        # worse true residual in the Phase 0 smoke
        use_x0 ? Krylov.krylov_solve!(ws, A, rhs, solver.x_prev; N=P, ldiv=true, common...) :
                 Krylov.krylov_solve!(ws, A, rhs; N=P, ldiv=true, common...)
    end

    solver.niter = ws.stats.niter
    solver.solved = ws.stats.solved

    # persist the solution as the next warmstart seed
    solver.x_prev .= ws.x
    solver.have_x_prev = true

    # drop the plan at solve end: releases the Cache buffers (the dominant
    # allocation) and guarantees the next solve cannot reuse a plan built
    # under a different core_size/geometry (persistent_plan keeps it — see
    # the entry clear above for the contract)
    solver.cache_tree && !solver.persistent_plan && (solver.kop.plan_slot[] = nothing)

    return nothing
end

"""
    transform_solver_geometry!(solver, body, R, t)

Mirror one RIGID body motion `x -> R*x + t` into any FMM state the solver
persists across solves, so that state remains valid without a rebuild:

- `KrylovSolver` with `persistent_plan=true`: transforms the live `FmmPlan`
  (trees, interaction lists, and any near-field cache) via
  `FastMultipole.transform_plan!`. No-op while no plan is alive yet (the
  first solve builds it lazily from current geometry).
- `FGSSolver`: transforms the `FastGaussSeidel` trees via
  `FastMultipole.transform_solver!` — this fixes the unsteady staleness bug
  where far-field expansions were silently formed about construction-time
  branch centers as rotation accumulated (BRAINSTORM 021
  rigid_motion_tree_reuse item). The dense self/nonself matrices are
  untouched (scalar rows are exactly rotation-invariant).
- every other solver: no-op (`Backslash`'s dense Dirichlet/Neumann operators
  are rotation-invariant; plain `KrylovSolver` rebuilds trees per apply).

`simulate!` calls this automatically for each body/solver pair after
`propagate_kinematics!` using the frame tree's per-step affine deltas; call
it manually only when driving rigid motion outside `simulate!`.
"""
transform_solver_geometry!(solver, body, R, t) = nothing

function transform_solver_geometry!(solver::KrylovSolver, body, R, t)
    (solver.persistent_plan && solver.kop.plan_slot[] !== nothing) || return nothing
    FastMultipole.transform_plan!(solver.kop.plan_slot[][1], (body,), R, t)
    return nothing
end

"""
    transform_body_solvers!(body_solvers, systems_tuple, transforms)

Forward each body's per-step rigid transform (from `propagate_kinematics!`)
to its solver via [`transform_solver_geometry!`](@ref). `transforms` is the
per-system vector of `(R, t)` affine deltas; identity entries are skipped.
Solver tuples that do not align one-to-one with the bodies (coupled solvers)
are left untouched — their per-apply rebuilds remain correct.
"""
function transform_body_solvers!(body_solvers, systems_tuple::Tuple, transforms)
    transforms === nothing && return nothing
    solvers = body_solvers isa Tuple ? body_solvers : (body_solvers,)
    length(solvers) == length(systems_tuple) || return nothing
    for (solver, body, (R, t)) in zip(solvers, systems_tuple, transforms)
        # skip identity deltas (static bodies) — cheap exact check: the
        # kinematics emit exact I / zero for frames with no motion
        R == I && iszero(t) && continue
        transform_solver_geometry!(solver, body, R, t)
    end
    return nothing
end

function _solve!(self::AbstractBody{<:Any,<:Any,<:Any,false}, solver::KrylovSolver{<:Any,B,TF}, Das=nothing; optargs...) where {B,TF}

    # save external velocity and refresh the operator's normals snapshot
    solver.Uext .= self.velocity
    solver.kop.normals .= self.normals

    # construct right-hand side
    calc_bc_noflowthrough!(solver.rhs, self.velocity, solver.kop.normals)

    # launch krylov solver on the persistent workspace
    _krylov_launch!(solver)

    # store solution
    solver.unabbreviated_strengths .= solver.workspace.x
    _set_strength(self, solver.unabbreviated_strengths)

    # restore external velocity
    self.velocity .= solver.Uext
    return nothing
end

function _solve!(self::AbstractBody{<:Any,2,<:Any,true}, solver::KrylovSolver{<:Any,B,TF}, Das=nothing; optargs...) where {B,TF}

    solver.kop.normals .= self.normals
    solver.source_strengths .= view(self.strength, :, 1)

    # construct right-hand side
    solver.rhs .= -self.potential

    # launch krylov solver on the persistent workspace
    _krylov_launch!(solver)

    solver.unabbreviated_strengths .= solver.workspace.x
    self.strength[:, 1] .= solver.source_strengths
    _set_strength(self, solver.unabbreviated_strengths)

    return nothing
end

"""
    KrylovCoupledOperator(bodies, backend, Ax)

Matrix-free application of the coupled multi-body influence operator as the
5-argument mul! callable LinearOperators.jl expects.
"""
struct KrylovCoupledOperator{TB<:Tuple,B<:AbstractBackend,TF}
    bodies::TB
    backend::B
    Ax::Vector{TF}
end

function (op::KrylovCoupledOperator)(C, B, α, β)
    _write_coupled_unknowns!(op.bodies, B)

    for body in op.bodies
        body.velocity .= zero(eltype(body.velocity))
        body.potential .= zero(eltype(body.potential))
    end

    scalar_potential = [has_dirichlet_bc(body) for body in op.bodies]
    velocity = [!has_dirichlet_bc(body) for body in op.bodies]
    influence!(op.bodies, op.bodies, op.backend; scalar_potential, velocity)
    _collect_coupled_operator!(op.Ax, op.bodies)

    C .*= β
    C .+= α .* op.Ax
    return nothing
end

"""
    KrylovCoupled(bodies; method=:gmres, itmax=20, atol=1e-6, rtol=1e-6,
                  backend=FastMultipoleBackend(), memory=20, warmstart=false)

Matrix-free coupled solver for a tuple of bodies. The linear operator and
Krylov workspace are built once at construction and reused across solves;
`solver.niter`/`solver.solved` report the last solve's stats, and
`warmstart=true` seeds each solve with the previous coupled solution.
"""
mutable struct KrylovCoupled{TB<:Tuple,B<:AbstractBackend,TF<:Number,TK,TA,TW} <: AbstractMatrixFreeSolver
    bodies::TB
    backend::B
    rhs::Vector{TF}
    x::Vector{TF}          # last coupled solution (warmstart seed)
    Ax::Vector{TF}
    kop::TK                # KrylovCoupledOperator
    A::TA                  # LinearOperators wrapper (built once)
    workspace::TW          # Krylov workspace (built once, reused)
    method::Symbol
    itmax::Int
    atol::Float64
    rtol::Float64
    memory::Int
    warmstart::Bool        # seed solves with the previous solution
    have_x_prev::Bool      # whether x holds a valid previous solution
    niter::Int             # iterations of the last solve
    solved::Bool           # convergence flag of the last solve
end

function KrylovCoupled(bodies::Tuple;
        method::Symbol=:gmres,
        itmax::Int=20,
        atol::Real=1e-6,
        rtol::Real=1e-6,
        backend::AbstractBackend=FastMultipoleBackend(),
        memory::Int=20,
        warmstart::Bool=false)

    TF = promote_type(map(numtype, bodies)...)
    ncs = sum(body -> body.ncells, bodies)
    rhs = zeros(TF, ncs)
    x = zeros(TF, ncs)
    Ax = zeros(TF, ncs)

    kop = KrylovCoupledOperator{typeof(bodies),typeof(backend),TF}(bodies, backend, Ax)
    prod! = (y, xv, α, β) -> kop(y, xv, α, β)
    A = LinearOperators.LinearOperator(TF, ncs, ncs, false, false, prod!)
    workspace = method in _KRYLOV_MEMORY_METHODS ?
        Krylov.krylov_workspace(Val(method), A, rhs; memory) :
        Krylov.krylov_workspace(Val(method), A, rhs)

    return KrylovCoupled{typeof(bodies),typeof(backend),TF,typeof(kop),typeof(A),typeof(workspace)}(
        bodies, backend, rhs, x, Ax, kop, A, workspace,
        method, Int(itmax), Float64(atol), Float64(rtol), Int(memory),
        warmstart, false, 0, false)
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

function solve!(bodies::Tuple, solver::KrylovCoupled; backend=solver.backend, optargs...)
    offsets = _coupled_offsets(bodies)
    velocity_old = [copy(body.velocity) for body in bodies]
    potential_old = [copy(body.potential) for body in bodies]
    fixed_sources = Vector{Any}(undef, length(bodies))

    try
        for (i, body) in enumerate(bodies)
            calc_normals!(body)
            calc_controlpoints!(body)
            set_strengths!(body)
            fixed_sources[i] = has_dirichlet_bc(body) ? copy(body.strength[:, 1]) : nothing
        end

        scalar_potential = [has_dirichlet_bc(body) for body in bodies]
        velocity = [!has_dirichlet_bc(body) for body in bodies]
        influence!(bodies, bodies, backend; scalar_potential, velocity, optargs...)
        for (i, body) in enumerate(bodies)
            if has_dirichlet_bc(body)
                body.potential .+= potential_old[i]
            end
        end
        boundary_condition!(bodies, solver, backend; optargs...)

        workspace = solver.workspace
        if solver.warmstart && solver.have_x_prev
            Krylov.krylov_solve!(workspace, solver.A, solver.rhs, solver.x;
                atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
        else
            Krylov.krylov_solve!(workspace, solver.A, solver.rhs;
                atol=solver.atol, rtol=solver.rtol, itmax=solver.itmax)
        end
        solver.niter = workspace.stats.niter
        solver.solved = workspace.stats.solved
        solver.x .= workspace.x
        solver.have_x_prev = true

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

    return nothing
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
    cache_leaf_lu::Bool
    sweep_order::Symbol
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
        cache_leaf_lu::Bool=true,
        sweep_order::Symbol=:lexicographic,  # :colored = parallel per-color sweeps (021 Phase 2b; changes the GS iteration)
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
        FastMultipole.FastGaussSeidel(bodies; expansion_order, multipole_acceptance, leaf_size, cache_leaf_lu, sweep_order, shrink, recenter, extra_farfield=any(has_semiinfinite_wake.(bodies))) :
        nothing

    Uext = zeros(TF, 3, body.ncells)
    phi_ext = zeros(TF, body.ncells)
    solution_history = zeros(TF, body.ncells, size(body.strength, 2), solution_history_length)
    return FGSSolver{typeof(fgs), TF}(fgs, Int(expansion_order), Int(leaf_size), Float64(multipole_acceptance), Bool(cache_leaf_lu), Symbol(sweep_order), max_iterations, Int(inner_iterations), Float64(tolerance), Float64(rlx), Bool(reverse_pass), Bool(verbose), Uext, phi_ext, solution_history, solution_history_length, 0, project_solution, project_solution_order)
end

################################################################################

function _solve!(body::AbstractBody{TK,NK,TF,false}, solver::Backslash;
        backend=DirectBackend(), strength_index=1,
        update_G::Bool=false,
        optargs...
    ) where {TK, NK, TF}

    if update_G
        solver.G .= zero(eltype(solver.G))
        _G!(solver.G, body, body; core_size=body.core_size_panel, update_geometry=false, optargs...)
        # write the refreshed factorization back so later direct consumers of
        # solver.Glu (Kutta Route A, Green-family formulations) never see a
        # stale-pivot factorization (Glu.factors aliases solver.G)
        solver.Glu = lu!(solver.G)
    end

    rhs = solver.rhs
    rhs .= zero(eltype(rhs))
    calc_bc_noflowthrough!(rhs, body.velocity, body.normals)

    ldiv!(view(body.strength, :, strength_index), solver.Glu, rhs)

    return nothing
end

function _solve!(self::AbstractBody{<:Union{Union{ConstantSource, ConstantDoublet}, Union{ConstantSource, VortexRing}}, 2, TF, true},
                solver::Backslash; backend=DirectBackend(),
                update_G::Bool=false, optargs...) where TF

    solver.rhs .= -self.potential

    if update_G
        solver.G .= 0.0
        _G!(solver.G, self, self; core_size=self.core_size_panel, update_geometry=false)
        # write back (see Neumann _solve! comment)
        solver.Glu = lu!(solver.G)
    end

    ldiv!(view(self.strength, :, 2), solver.Glu, solver.rhs)

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

# rigid-motion hook (see transform_solver_geometry! docstring above): fixes
# the FGS unsteady staleness bug — far-field expansions about
# construction-time branch centers under accumulated rotation
function transform_solver_geometry!(solver::FGSSolver, body, R, t)
    FastMultipole.transform_solver!(solver.fgs, (body,), R, t)
    return nothing
end

transform_solver_geometry!(solver::FGSSolver{Nothing}, body, R, t) = nothing

function _solve!(body::AbstractBody, solver::FGSSolver; backend = FastMultipoleBackend(
        expansion_order=solver.expansion_order,
        multipole_acceptance=solver.multipole_acceptance,
        leaf_size=solver.leaf_size
    ), callback=nothing, optargs...)

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
        final_update=false,
        callback
    )

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
end

"""
    _solve_history!(body, solver::FGSSolver, history::ConvergenceHistory; optargs...)

Run the production FGS `_solve!` while recording the per-outer-iteration
residual (the exact value the loop's tolerance check uses; max-abs, labeled
`:fgs_maxabs`) with wall-clock timestamps into `history` (021 ruling 4). The
solve itself is bit-identical to `_solve!` without recording.
"""
function _solve_history!(body::AbstractBody, solver::FGSSolver,
                         history::ConvergenceHistory; optargs...)
    reset!(history; metric=:fgs_maxabs)
    _solve!(body, solver;
            callback=(iteration, residual) -> record!(history, iteration, residual),
            optargs...)
    return history
end

################################################################################
# FGS PRECONDITIONER (021 roster config 1f)
################################################################################

"""
    FGSPreconditioner(body; sweeps=1, inner_iterations=2, rlx=1.0,
                      expansion_order=7, multipole_acceptance=0.4,
                      leaf_size=10, cache_leaf_lu=true,
                      shrink=false, recenter=false)

Fast-Gauss–Seidel preconditioner for [`KrylovSolver`](@ref): approximately
applies `G⁻¹` through a fixed number (`sweeps`) of FGS outer iterations, each
one FMM far-field pass plus `inner_iterations` per-leaf dense sweeps. The
FastGaussSeidel tree and leaf matrices are built once at construction and
reused for every apply; geometry must not change afterwards.

Pass it as `KrylovSolver(body; method=:fgmres, preconditioner=...)`. It is
routed as a *flexible right* preconditioner (`N=`, `ldiv=true`), so plain
GMRES's assumptions are never violated and the monitored residual stays the
true one. Each apply starts from zero strengths, making the map linear in its
input for a fixed sweep count.

By default, each dense leaf block is LU-factorized once at construction and
reused on every sweep. Pass `cache_leaf_lu=false` to retain the lower-memory
per-sweep factorization path.
"""
struct FGSPreconditioner{TS<:FGSSolver,TB<:AbstractBody,TF}
    solver::TS                 # fixed-sweep FGSSolver (max_iterations=sweeps, tolerance=0)
    body::TB
    strength_save::Matrix{TF}  # body-state save/restore buffers (side-effect-free apply)
    velocity_save::Matrix{TF}
    potential_save::Vector{TF}
end

function FGSPreconditioner(body::AbstractBody;
        sweeps::Int=1,             # FGS outer iterations per apply (each includes one fmm! pass)
        inner_iterations::Int=2,   # per-leaf Gauss-Seidel sweeps per outer iteration
        rlx::Real=1.0,             # under-relaxation inside the sweeps
        expansion_order=7,
        multipole_acceptance=0.4,
        leaf_size=10,
        cache_leaf_lu::Bool=true,
        sweep_order::Symbol=:lexicographic,
        shrink=false,
        recenter=false,
    )

    fgssolver = FGSSolver(body;
        max_iterations=sweeps, inner_iterations, tolerance=0.0, rlx,
        reverse_pass=false, verbose=false, expansion_order,
        multipole_acceptance, leaf_size, cache_leaf_lu, sweep_order, shrink, recenter)

    TF = numtype(body)
    return FGSPreconditioner{typeof(fgssolver),typeof(body),TF}(
        fgssolver, body,
        zeros(TF, size(body.strength)),
        zeros(TF, 3, body.ncells),
        zeros(TF, body.ncells))
end

"""
    ldiv!(y, P::FGSPreconditioner, x)

Approximate `y = G⁻¹ x` with `P.solver.max_iterations` FGS sweeps: the input
is packed into the body's boundary-condition channel (Neumann: velocity such
that `-u·n = x`; Dirichlet: `potential = -x` with sources zeroed), the
fixed-sweep FGS solve is run from zero strengths, and the solved strength
column is copied out. The body's `strength`/`velocity`/`potential` are saved
and restored, so the apply is side-effect-free.
"""
function LA.ldiv!(y::AbstractVector, P::FGSPreconditioner, x::AbstractVector)
    body = P.body
    dirichlet = has_dirichlet_bc(body)
    strength_index = _fgs_solved_strength_index(body)

    P.strength_save .= body.strength
    P.velocity_save .= body.velocity
    P.potential_save .= body.potential

    try
        # zero start: makes the fixed-sweep map linear in x
        body.strength .= zero(eltype(body.strength))

        if dirichlet
            # FGS right-hand side is -potential, so potential = -x gives rhs = x
            body.potential .= .-x
        else
            # FGS right-hand side is -u·n, so u = -x n gives rhs = x (unit normals)
            for d in 1:3
                @views @. body.velocity[d, :] = -x * body.normals[d, :]
            end
        end

        _solve!(body, P.solver)

        y .= view(body.strength, :, strength_index)
    finally
        body.strength .= P.strength_save
        body.velocity .= P.velocity_save
        body.potential .= P.potential_save
    end

    return y
end

################################################################################
# FMM-DIRECT-LIST ILU(0) PRECONDITIONER
################################################################################

"""
    ILUPreconditioner(body; leaf_size=10, multipole_acceptance=0.4,
        max_pattern_entries=512 * body.ncells, equilibrate=false,
        diagonal_shift=0.0, keep_matrix=false)

Build an ILU(0) right preconditioner from the directed panel pairs in a
dedicated `FastMultipole.Barba()` direct-interaction list.  The sparse
operator is assembled and factored in tree order; `ldiv!` transparently
permutes vectors to and from the body's panel order.

Construction is proportional to the direct-list pattern.  It never forms or
scans an `N×N` object, and rejects a pattern larger than
`max_pattern_entries` before allocating sparse triplets or evaluating panel
kernels.  Geometry must remain fixed after construction.
"""
mutable struct ILUPreconditioner{TF,TFACT,TMAT,TSTATS}
    fact::TFACT
    permutation::Vector{Int}          # tree index -> original body index
    inverse_permutation::Vector{Int}  # original body index -> tree index
    rhs_tree::Vector{TF}
    solution_tree::Vector{TF}
    row_scale::Vector{TF}
    matrix::TMAT                      # `nothing` unless keep_matrix=true
    leaf_size::Int
    multipole_acceptance::Float64
    interaction_list_method::Symbol
    equilibrate::Bool
    diagonal_shift::TF
    max_pattern_entries::Int
    stats::TSTATS
end

"""
Return tree/list data, the guarded number of expanded panel pairs, and the
tree/list build times in seconds.  Restores the body's normals/control points
on exit; callers that assemble afterwards must refresh geometry themselves.
"""
function _ilu_direct_pattern(body::AbstractBody, leaf_size::Int,
                             multipole_acceptance::Real,
                             max_pattern_entries::Int)
    leaf_size > 0 || throw(ArgumentError("leaf_size must be positive"))
    0 <= multipole_acceptance <= 1 || throw(ArgumentError(
        "multipole_acceptance must lie in [0, 1] (1 accepts all eligible " *
        "expansions; 0 accepts none)"))
    max_pattern_entries >= body.ncells || throw(ArgumentError(
        "max_pattern_entries=$max_pattern_entries is smaller than N=$(body.ncells); " *
        "the required diagonal alone has N entries"))

    normals_save = copy(body.normals)
    controlpoints_save = copy(body.controlpoints)
    local target_tree, source_tree, direct_list, diagonal_seen, requested
    t_tree = 0.0
    t_lists = 0.0
    try
        calc_normals!(body)
        calc_controlpoints!(body)
        method = FastMultipole.Barba()
        TF = numtype(body)
        derivatives = has_dirichlet_bc(body) ?
            FastMultipole.DerivativesSwitch(true, false, false, (body,)) :
            FastMultipole.DerivativesSwitch(false, true, false, (body,))
        leaf_sizes = FastMultipole.StaticArrays.SVector{1,Int}(leaf_size)

        t0 = time_ns()
        target_tree = FastMultipole.Tree((body,), true, derivatives, TF;
            leaf_size=leaf_sizes, shrink=true, interaction_list_method=method)
        source_tree = FastMultipole.Tree((body,), false, derivatives, TF;
            leaf_size=leaf_sizes, shrink=true, interaction_list_method=method)
        t_tree = (time_ns() - t0) / 1e9

        t0 = time_ns()
        _, direct_list = FastMultipole.build_interaction_lists(
            target_tree.branches, source_tree.branches, leaf_sizes,
            multipole_acceptance, true, true, true, method)
        t_lists = (time_ns() - t0) / 1e9

    # A repeated branch pair would create repeated expanded panel pairs.  This
    # check is proportional to the branch list, not to N².
        length(Set(Tuple(pair) for pair in direct_list)) == length(direct_list) ||
            error("FastMultipole Barba direct_list contains duplicate branch pairs")

        target_sort = target_tree.sort_index_list[1]
        source_sort = source_tree.sort_index_list[1]

        # Size the pattern BEFORE walking it. This pass is arithmetic over
        # branch pairs (cheap), whereas the diagonal pass below is O(entries);
        # doing it first means an over-cap request reports the TOTAL it needs
        # instead of the running subtotal at which it happened to trip, so the
        # limit can be set in one shot rather than by bisection.
        requested = 0
        for pair in direct_list
            requested += length(target_tree.branches[pair[1]].bodies_index[1]) *
                         length(source_tree.branches[pair[2]].bodies_index[1])
        end
        requested <= max_pattern_entries || throw(ArgumentError(
            "Barba direct-list pattern requests $requested entries " *
            "($(round(requested / body.ncells; digits=1)) per row) for N=$(body.ncells), " *
            "exceeding max_pattern_entries=$max_pattern_entries " *
            "($(round(max_pattern_entries / body.ncells; digits=1)) per row). Increase the " *
            "explicit limit or adjust leaf_size/MAC; construction stopped before sparse " *
            "allocation and kernel evaluation."))

        diagonal_seen = falses(body.ncells)
        for pair in direct_list
            target_range = target_tree.branches[pair[1]].bodies_index[1]
            source_range = source_tree.branches[pair[2]].bodies_index[1]
            for it in target_range, js in source_range
                target_sort[it] == source_sort[js] &&
                    (diagonal_seen[target_sort[it]] = true)
            end
        end
        requested += count(!, diagonal_seen)
        requested <= max_pattern_entries || throw(ArgumentError(
            "Barba direct-list pattern plus required diagonal requests $requested entries, " *
            "exceeding max_pattern_entries=$max_pattern_entries for N=$(body.ncells); " *
            "construction stopped before sparse allocation and kernel evaluation."))
    finally
        body.normals .= normals_save
        body.controlpoints .= controlpoints_save
    end

    return target_tree, source_tree, direct_list, diagonal_seen, requested,
           t_tree, t_lists
end

"Assemble the direct-list operator in one common (source-tree) ordering."
function _assemble_ilu_operator(body::AbstractBody{<:Any,NK,TF,DBC},
        target_tree, source_tree, direct_list, diagonal_seen, requested;
        core_size=body.core_size_panel) where {NK,TF,DBC}
    n = body.ncells
    permutation = copy(source_tree.sort_index_list[1])
    inverse_permutation = invperm(permutation)
    target_sort = target_tree.sort_index_list[1]
    source_sort = source_tree.sort_index_list[1]

    rows = Vector{Int}(undef, requested)
    cols = Vector{Int}(undef, requested)
    vals = Vector{TF}(undef, requested)

    derivatives = DBC ? FastMultipole.DerivativesSwitch(true, false, false) :
                        FastMultipole.DerivativesSwitch(false, true, false)
    _, strength_index = _G_kernel_and_strength_index(body)
    old_strength = copy(body.strength)
    correction_was_active = _operator_mode_begin!(body)
    body.strength .= zero(eltype(body.strength))
    body.strength[:, strength_index] .= one(eltype(body.strength))

    # Per-pair triplet offsets so the kernel evaluations can run threaded
    # (mirroring `_G!`; `induced` is pure, so concurrent evaluations of the
    # same source panel from different pairs are safe).
    offsets = Vector{Int}(undef, length(direct_list))
    k = 0
    for (p, pair) in enumerate(direct_list)
        offsets[p] = k
        k += length(target_tree.branches[pair[1]].bodies_index[1]) *
             length(source_tree.branches[pair[2]].bodies_index[1])
    end
    try
        CPs = body.controlpoints
        normals = body.normals
        Threads.@threads for p in eachindex(direct_list)
            pair = direct_list[p]
            target_range = target_tree.branches[pair[1]].bodies_index[1]
            source_range = source_tree.branches[pair[2]].bodies_index[1]
            kp = offsets[p]
            for js in source_range
                j_original = source_sort[js]
                j_tree = inverse_permutation[j_original]
                for it in target_range
                    i_original = target_sort[it]
                    i_tree = inverse_permutation[i_original]
                    target = FastMultipole.StaticArrays.SVector{3,TF}(
                        CPs[1, i_original], CPs[2, i_original], CPs[3, i_original])
                    phi, u, _ = induced(target, body, j_original, derivatives;
                                        core_size)
                    value = DBC ? phi :
                        u[1] * normals[1, i_original] +
                        u[2] * normals[2, i_original] +
                        u[3] * normals[3, i_original]
                    isfinite(value) || error("nonfinite ILU operator entry at " *
                        "body pair ($i_original, $j_original)")
                    kp += 1
                    rows[kp], cols[kp], vals[kp] = i_tree, j_tree, value
                end
            end
        end

        # Barba+self_induced normally supplies every diagonal.  Add any that
        # is absent rather than assuming a particular FastMultipole version.
        for original in eachindex(diagonal_seen)
            diagonal_seen[original] && continue
            tree_index = inverse_permutation[original]
            target = FastMultipole.StaticArrays.SVector{3,TF}(
                CPs[1, original], CPs[2, original], CPs[3, original])
            phi, u, _ = induced(target, body, original, derivatives; core_size)
            value = DBC ? phi : dot(u, view(normals, :, original))
            isfinite(value) || error("nonfinite ILU diagonal at body panel $original")
            k += 1
            rows[k], cols[k], vals[k] = tree_index, tree_index, value
        end
    finally
        _operator_mode_end!(body, correction_was_active)
        body.strength .= old_strength
    end
    k == requested || error("internal ILU pattern count mismatch: assembled $k of $requested")
    S = SparseArrays.sparse(rows, cols, vals, n, n)
    SparseArrays.nnz(S) == requested || error(
        "expanded Barba direct-list contains duplicate panel pairs: requested " *
        "$requested entries but sparse pattern has $(SparseArrays.nnz(S))")
    return S, permutation, inverse_permutation
end

function _ilu_row_scale(S::SparseArrays.SparseMatrixCSC{TF,Int}) where TF
    rowmax = zeros(TF, size(S, 1))
    for col in axes(S, 2), p in SparseArrays.nzrange(S, col)
        row = S.rowval[p]
        rowmax[row] = max(rowmax[row], abs(S.nzval[p]))
    end
    all(isfinite, rowmax) || error("nonfinite row norm in ILU operator")
    any(iszero, rowmax) && error(
        "ILU operator has an empty/zero row; try a larger leaf_size or different MAC")
    return inv.(rowmax)
end

function _check_ilu_pivots(fact, scale)
    pivots = [fact.u_nzval[fact.u_colptr[j + 1] - 1] for j in 1:fact.n]
    all(isfinite, pivots) || error(
        "ILU(0) produced a nonfinite pivot; try equilibrate=true or diagonal_shift")
    tolerance = 100 * eps(real(float(one(eltype(pivots))))) * max(one(scale), scale)
    bad = findfirst(p -> abs(p) <= tolerance, pivots)
    bad === nothing || error("ILU(0) has a missing/small pivot at tree row $bad " *
        "(|pivot|=$(abs(pivots[bad])) ≤ $tolerance); try equilibrate=true or " *
        "a nonzero diagonal_shift")
    return nothing
end

function ILUPreconditioner(body::AbstractBody;
        leaf_size::Int=10,
        multipole_acceptance::Real=0.4,
        max_pattern_entries::Integer=512 * body.ncells,
        equilibrate::Bool=false,
        diagonal_shift::Real=0.0,
        keep_matrix::Bool=false)
    total_start = time_ns()
    n = body.ncells
    TF = numtype(body)
    max_entries = Int(max_pattern_entries)

    # Constructor-side geometry refreshes are restored along with all operator
    # workspaces, even if the pattern guard or factorization throws.
    strength_save = copy(body.strength)
    velocity_save = copy(body.velocity)
    potential_save = copy(body.potential)
    normals_save = copy(body.normals)
    controlpoints_save = copy(body.controlpoints)

    local S, permutation, inverse_permutation, fact, row_scale
    t_assembly = 0.0
    t_factorization = 0.0
    try
        # Argument validation, tree/list construction, duplicate-branch check,
        # and the linear pattern-size guard all live in the shared helper.
        target_tree, source_tree, direct_list, diagonal_seen, requested,
            t_tree, t_lists = _ilu_direct_pattern(
                body, leaf_size, multipole_acceptance, max_entries)

        # The helper restores geometry on exit; refresh again so assembly
        # evaluates kernels on the same geometry the trees were built from.
        calc_normals!(body)
        calc_controlpoints!(body)

        t0 = time_ns()
        S, permutation, inverse_permutation = _assemble_ilu_operator(
            body, target_tree, source_tree, direct_list, diagonal_seen, requested)
        diagonal_shift_tf = convert(TF, diagonal_shift)
        if !iszero(diagonal_shift_tf)
            for i in 1:n
                S[i, i] += diagonal_shift_tf
            end
        end
        t_assembly = (time_ns() - t0) / 1e9

        row_scale = equilibrate ? _ilu_row_scale(S) : ones(TF, n)
        Sfactor = equilibrate ? SparseArrays.spdiagm(0 => row_scale) * S : S
        # ILUZero assumes every U column ends in a diagonal entry.
        all(i -> Sfactor[i, i] != 0, 1:n) || error(
            "ILU operator has a missing/zero diagonal; try equilibrate=true or " *
            "a nonzero diagonal_shift")
        t0 = time_ns()
        fact = ILUZero.ilu0(Sfactor)
        _check_ilu_pivots(fact, maximum(abs, Sfactor.nzval))
        t_factorization = (time_ns() - t0) / 1e9

        row_nnz = zeros(Int, n)
        for row in S.rowval
            row_nnz[row] += 1
        end
        stats = Dict{String,Any}(
            "n" => n,
            "nnz" => SparseArrays.nnz(S),
            "nnz_per_panel" => SparseArrays.nnz(S) / n,
            "nnz_fraction" => SparseArrays.nnz(S) / (float(n) * float(n)),
            "max_row_nnz" => maximum(row_nnz),
            "factor_nnz" => SparseArrays.nnz(fact),
            "tree_time" => t_tree,
            "interaction_list_time" => t_lists,
            "assembly_time" => t_assembly,
            "factorization_time" => t_factorization,
            "total_time" => (time_ns() - total_start) / 1e9,
            "retained_bytes" => 0,
        )
        kept = keep_matrix ? S : nothing
        P = ILUPreconditioner{TF,typeof(fact),typeof(kept),typeof(stats)}(
            fact, permutation, inverse_permutation, zeros(TF, n), zeros(TF, n),
            row_scale, kept, leaf_size, Float64(multipole_acceptance), :Barba,
            equilibrate, convert(TF, diagonal_shift), max_entries, stats)
        stats["retained_bytes"] = Base.summarysize(P)
        stats["total_time"] = (time_ns() - total_start) / 1e9
        return P
    finally
        body.strength .= strength_save
        body.velocity .= velocity_save
        body.potential .= potential_save
        body.normals .= normals_save
        body.controlpoints .= controlpoints_save
    end
end

"Apply the ILU factors with body↔tree permutation and optional row scaling."
function LA.ldiv!(y::AbstractVector, P::ILUPreconditioner, x::AbstractVector)
    n = length(P.permutation)
    length(x) == n && length(y) == n || throw(DimensionMismatch(
        "ILUPreconditioner expects vectors of length $n"))
    @inbounds for tree_index in 1:n
        P.rhs_tree[tree_index] =
            P.row_scale[tree_index] * x[P.permutation[tree_index]]
    end
    LA.ldiv!(P.solution_tree, P.fact, P.rhs_tree)
    @inbounds for tree_index in 1:n
        y[P.permutation[tree_index]] = P.solution_tree[tree_index]
    end
    return y
end

function solve!(bodies::Tuple, solvers::Tuple;
    backend = fill(DirectBackend(), length(bodies)),
    max_outer_iterations::Int = 50,
    outer_tolerance::Real = 1e-8,
    verbose::Bool = false,
    history::Union{Nothing,ConvergenceHistory} = nothing,
    optargs...)

    # println("Tuple of bodies")

    history === nothing || reset!(history; metric=:blockgs_maxdelta)

    N = length(bodies)
    @assert length(solvers) == N "Number of solvers ($(length(solvers))) must match number of bodies ($N)"
    backends = backend isa Tuple || backend isa AbstractVector ? backend : fill(backend, N)

    prev_velocity = [copy(body.velocity) for body in bodies]
    prev_strengths = [copy(body.strength) for body in bodies]

    converged = false
    for iter in 1:max_outer_iterations

        for (i, (body, solver)) in enumerate(zip(bodies, solvers))
            body.velocity .= prev_velocity[i]

            sources = tuple((bodies[j] for j in eachindex(bodies) if j != i)...)
            if !isempty(sources)
                influence!((body,), sources, backends[i];
                    scalar_potential=false,
                    velocity=true,
                    optargs...)
            end

            solve!(body, solver; backend=backends[i])
        end

        max_delta = 0.0
        for (i, body) in enumerate(bodies)
            prev_strengths[i] .-= body.strength
            prev_strengths[i] .= abs.(prev_strengths[i])
            max_delta = max(max_delta, maximum(prev_strengths[i]))
            prev_strengths[i] .= body.strength
        end

        history === nothing || record!(history, iter, max_delta)

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

    return nothing
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
    built::Bool          # whether G/Glu hold a real (non-dummy) factorization
    t_build::Float64     # seconds spent assembling+factorizing G in the last build
    t_solve::Float64     # seconds spent in ldiv! in the last solve
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

    Glu = lu!(G)  # dummy init; replaced on the first solve (built=false) or update_G=true

    BackslashCoupled{TF}(G, Glu, rhs, Uext, phi_ext, false, 0.0, 0.0)
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

function solve!(bodies::Tuple, solver::BackslashCoupled; backend=DirectBackend(), update_G::Bool=false, optargs...)

    # println("BackslashCoupled")
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

    influence!(bodies, bodies, backend; scalar_potential=[has_dirichlet_bc(target) for target in bodies], velocity=[!has_dirichlet_bc(target) for target in bodies], optargs...)

    for (bi, body) in enumerate(bodies)
        if has_dirichlet_bc(body)
            r = offsets[bi]+1 : offsets[bi+1]
            @views body.potential .+= solver.phi_ext[r]
        end
    end

    ### get the boundary_condition for each body to write the RHS
    boundary_condition!(bodies, solver, backend)

    # The constructor's Glu is a dummy identity factorization: build on the
    # first solve even without update_G (a first solve against the dummy would
    # silently return the RHS as the solution).
    if update_G || !solver.built
        # Zero G matrix
        fill!(solver.G, 0)

        # Build and factorize G matrix, caching the factorization in solver
        solver.t_build = @elapsed begin
            for (bi, source) in enumerate(bodies)
                c = offsets[bi]+1 : offsets[bi+1] # columns of sources
                for (ti, target) in enumerate(bodies)
                    r = offsets[ti]+1 : offsets[ti+1] # rows of targets
                    _G!(view(solver.G, r, c), target,
                        source;
                        update_geometry=false)
                end
            end
            solver.Glu = lu!(solver.G)
        end
        solver.built = true
    end

    # solve with cached LU
    sol = similar(solver.rhs)
    solver.t_solve = @elapsed ldiv!(sol, solver.Glu, solver.rhs)

    # write solution back
    for (bi, b) in enumerate(bodies)
        r = offsets[bi]+1 : offsets[bi+1]
        write_solution!(b, view(sol, r))

        @views b.velocity  .= solver.Uext[:, r]
        @views b.potential .= solver.phi_ext[r]
    end

    return nothing
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
