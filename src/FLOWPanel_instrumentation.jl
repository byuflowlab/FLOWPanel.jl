#=##############################################################################
# DESCRIPTION
    Solver instrumentation: shared operator application, true-residual
    evaluation, and convergence-history capture (BRAINSTORM 021 Phase 0).

    Every solver in the benchmark roster is judged by the same physical
    residual r = G x - b computed here, never by a solver's internal metric
    (FGS reports a max-abs residual labeled MSE; Krylov a preconditioned norm).

# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Aug 2026
  * License     : GNU Public License
=###############################################################################

################################################################################
# SHARED OPERATOR APPLICATION
################################################################################

# Compute G*x for a Neumann body (normal-velocity operator) using the same
# influence path the matrix-free solvers use. Overwrites `body.velocity` and
# the solved strength column; the result is returned as a view into row 1 of
# `body.velocity`. `strengths_scratch` receives a copy of `x` (handles abstract
# vector inputs from Krylov workspaces).
function _apply_neumann_G!(body::AbstractBody{<:Any,<:Any,<:Any,false},
                           x::AbstractVector, backend::AbstractBackend,
                           normals::AbstractMatrix,
                           strengths_scratch::AbstractVector;
                           plan_slot=nothing, cache_nearfield::Bool=false)

    strengths_scratch .= x
    _set_strength(body, strengths_scratch)

    # `induced` (called via `influence!`) returns the side-aware self limit
    # at self pairs, so no extra jump add is needed.
    body.velocity .= 0
    # plan_slot: FMM-plan reuse across applies (KrylovSolver cache_tree);
    # only forwarded when active so non-FMM backends never see the kwarg
    if plan_slot === nothing
        influence!(body, body, backend; velocity=true)
    else
        influence!(body, body, backend; velocity=true, plan_slot, cache_nearfield)
    end

    # dot product with normals, accumulated into row 1
    body.velocity .*= normals
    body.velocity[1, :] .+= view(body.velocity, 2, :)
    body.velocity[1, :] .+= view(body.velocity, 3, :)

    return view(body.velocity, 1, :)
end

# Compute G*x for a Dirichlet body (doublet-potential operator): sources are
# zeroed so only the solved column contributes. Overwrites `body.potential`,
# `body.strength[:, 1]`, and the solved strength column; the result is returned
# as `body.potential`.
function _apply_dirichlet_G!(body::AbstractBody{<:Any,<:Any,<:Any,true},
                             x::AbstractVector, backend::AbstractBackend,
                             strengths_scratch::AbstractVector;
                             plan_slot=nothing, cache_nearfield::Bool=false)

    body.strength[:, 1] .= 0
    strengths_scratch .= x
    _set_strength(body, strengths_scratch)

    body.potential .= 0
    if plan_slot === nothing
        influence!(body, body, backend; scalar_potential=true, velocity=false)
    else
        influence!(body, body, backend; scalar_potential=true, velocity=false,
                   plan_slot, cache_nearfield)
    end

    return body.potential
end

"""
    apply_G!(y, body, x, backend=DirectBackend(); scratch=similar(y))

Compute `y = G x`, where `G` is the body's self-influence operator in the
body's boundary-condition formulation (normal velocity for Neumann bodies,
doublet potential for Dirichlet bodies with sources zeroed). This is the same
code path the matrix-free solvers apply, so it is valid for every solver.

Overwrites the body's solved strength column and `velocity`/`potential`
workspace; use [`true_residual!`](@ref) for a side-effect-free evaluation.
"""
function apply_G!(y::AbstractVector, body::AbstractBody{<:Any,<:Any,<:Any,false},
                  x::AbstractVector, backend::AbstractBackend=DirectBackend();
                  scratch::AbstractVector=similar(y, body.ncells))
    Gx = _apply_neumann_G!(body, x, backend, body.normals, scratch)
    y .= Gx
    return y
end

function apply_G!(y::AbstractVector, body::AbstractBody{<:Any,<:Any,<:Any,true},
                  x::AbstractVector, backend::AbstractBackend=DirectBackend();
                  scratch::AbstractVector=similar(y, body.ncells))
    Gx = _apply_dirichlet_G!(body, x, backend, scratch)
    y .= Gx
    return y
end

################################################################################
# TRUE-RESIDUAL EVALUATOR (021 ruling 3)
################################################################################

"""
    assemble_rhs!(b, body)

Assemble the right-hand side of the body's linear system into `b`, in the
solve's native physical units.

- Neumann bodies: the impermeability RHS `-U ⋅ n` from the *current*
  `body.velocity` and `body.normals` — valid whenever `body.velocity` holds
  the external velocity at the control points.
- Dirichlet bodies: the RHS is `-body.potential` *at the instant the solver
  was invoked*, when `body.potential` holds the assembled source-influence
  workspace (see the Dirichlet `solve!`). Outside a solve that state is not
  reproducible from the body alone, so benchmark harnesses must capture `b`
  during their frozen-RHS setup; this method simply reads the current
  `body.potential`.
"""
function assemble_rhs!(b::AbstractVector, body::AbstractBody{<:Any,<:Any,<:Any,false})
    calc_bc_noflowthrough!(b, body.velocity, body.normals)
    return b
end

function assemble_rhs!(b::AbstractVector, body::AbstractBody{<:Any,<:Any,<:Any,true})
    b .= .-body.potential
    return b
end

"""
    true_residual!(r, body, x, b; backend=DirectBackend()) -> (rms, rmax)

Evaluate the unpreconditioned true residual `r = G x - b` of the body's linear
system for candidate solution `x` and right-hand side `b`, in physical units,
through one direct influence evaluation (021 ruling 3). Identical code path
for every solver; matrix-free solvers get their residual through the same one
influence pass.

Side-effect-free: the body's `strength`, `velocity`, and `potential` are
saved and restored around the evaluation.

Returns `(rms, rmax)` where `rms = ‖r‖₂ / √n` and `rmax = ‖r‖∞`.
"""
function true_residual!(r::AbstractVector, body::AbstractBody,
                        x::AbstractVector, b::AbstractVector;
                        backend::AbstractBackend=DirectBackend())

    strength_old = copy(body.strength)
    velocity_old = copy(body.velocity)
    potential_old = copy(body.potential)

    try
        apply_G!(r, body, x, backend)
        r .-= b
    finally
        body.strength .= strength_old
        body.velocity .= velocity_old
        body.potential .= potential_old
    end

    rms = LA.norm(r) / sqrt(length(r))
    rmax = maximum(abs, r)
    return rms, rmax
end

################################################################################
# CONVERGENCE-HISTORY CAPTURE (021 ruling 4)
################################################################################

"""
    ConvergenceHistory(metric::Symbol=:unset)

Per-iteration convergence record shared by all iterative solvers: iteration
number, wall-clock timestamp (`time_ns()`), and the solver's internal residual
metric, labeled by `metric` so downstream plots never conflate incompatible
norms:

- `:krylov_precnorm` — Krylov.jl residual history (preconditioned norm when a
  left preconditioner is active; true 2-norm otherwise)
- `:fgs_maxabs` — FastGaussSeidel residual (max-abs, despite its MSE print label)
- `:blockgs_maxdelta` — tuple block-Gauss–Seidel max strength change

Call [`reset!`](@ref) before a solve to stamp `t0_ns`, and [`record!`](@ref)
once per iteration.
"""
mutable struct ConvergenceHistory
    iter::Vector{Int}
    t_ns::Vector{UInt64}
    residual_internal::Vector{Float64}
    metric::Symbol
    t0_ns::UInt64
end

ConvergenceHistory(metric::Symbol=:unset) =
    ConvergenceHistory(Int[], UInt64[], Float64[], metric, UInt64(0))

"""
    reset!(history::ConvergenceHistory; metric=history.metric)

Empty the record, relabel its metric, and stamp `t0_ns = time_ns()`.
"""
function reset!(history::ConvergenceHistory; metric::Symbol=history.metric)
    empty!(history.iter)
    empty!(history.t_ns)
    empty!(history.residual_internal)
    history.metric = metric
    history.t0_ns = time_ns()
    return history
end

"""
    record!(history::ConvergenceHistory, iter, residual)

Append one iteration with the current wall-clock timestamp.
"""
function record!(history::ConvergenceHistory, iter::Integer, residual::Real)
    push!(history.iter, Int(iter))
    push!(history.t_ns, time_ns())
    push!(history.residual_internal, Float64(residual))
    return history
end

Base.length(history::ConvergenceHistory) = length(history.iter)

"""
    elapsed_seconds(history::ConvergenceHistory)

Per-iteration wall-clock times in seconds relative to the last `reset!`.
"""
elapsed_seconds(history::ConvergenceHistory) = (history.t_ns .- history.t0_ns) ./ 1e9
