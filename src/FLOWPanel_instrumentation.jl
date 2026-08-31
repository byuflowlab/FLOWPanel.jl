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
# ENV-GATED STEP TIMERS (task 052)
################################################################################

const _STEP_TIMERS = Ref{Union{Nothing,Bool}}(nothing)
const _STEP_TIMER_TLS_KEY = :flowpanel_step_timer_accumulator

"Whether FLOWPANEL_STEP_TIMERS is armed (read lazily, cached)."
function step_timers_enabled()
    if _STEP_TIMERS[] === nothing
        _STEP_TIMERS[] = lowercase(get(ENV, "FLOWPANEL_STEP_TIMERS", "false")) in
            ("1", "true", "on")
    end
    return _STEP_TIMERS[]::Bool
end

"Synchronize CUDA at an instrumented timing boundary when the CUDA seam is live."
function _step_timer_sync()
    step_timers_enabled() || return nothing
    GPU_INFLUENCE[] === :cuda || return nothing
    isdefined(FastMultipole, :CUDA) || return nothing
    getglobal(FastMultipole, :CUDA).synchronize()
    return nothing
end

function _step_timer_start()
    _step_timer_sync()
    return time_ns()
end

function _step_timer_elapsed(t0::UInt64)
    _step_timer_sync()
    return (time_ns() - t0) / 1e9
end

"Add time to the current step's exclusive or nested accumulator."
function _step_timer_add!(label::Symbol, seconds::Real; nested::Bool=false)
    step_timers_enabled() || return nothing
    state = get(task_local_storage(), _STEP_TIMER_TLS_KEY, nothing)
    state === nothing && return nothing
    bucket = nested ? state.nested : state.exclusive
    bucket[label] = get(bucket, label, 0.0) + Float64(seconds)
    return nothing
end

function _step_timer_measure(f::F, label::Symbol; nested::Bool=false) where F
    step_timers_enabled() || return f()
    t0 = _step_timer_start()
    try
        return f()
    finally
        _step_timer_add!(label, _step_timer_elapsed(t0); nested)
    end
end

function _step_timer_begin_step!()
    step_timers_enabled() || return nothing
    state = (exclusive=Dict{Symbol,Float64}(), nested=Dict{Symbol,Float64}())
    task_local_storage(_STEP_TIMER_TLS_KEY, state)
    _step_timer_sync()
    return (state=state, start_ns=time_ns())
end

function _step_timer_finish_step!(token, i_step::Integer)
    token === nothing && return nothing
    _step_timer_sync()
    total = (time_ns() - token.start_ns) / 1e9
    labels = (
        :controls_setup, :wake_influence, :solve, :body_influence,
        :remaining_aerodynamics, :monitors, :io,
        :wake_propagation_maintenance, :rigid_kinematics, :shedding,
    )
    classified = sum(get(token.state.exclusive, label, 0.0) for label in labels)
    residual = total - classified
    for label in labels
        @info "step_timer $(label) step=$(i_step) $(get(token.state.exclusive, label, 0.0)) s"
    end
    for (label, seconds) in sort!(collect(token.state.nested); by=first)
        @info "step_timer_nested $(label) step=$(i_step) $(seconds) s"
    end
    @info "step_timer total_step step=$(i_step) $(total) s"
    @info "step_timer unclassified_residual step=$(i_step) $(residual) s"
    task_local_storage(_STEP_TIMER_TLS_KEY, nothing)
    return nothing
end

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
                           plan_slot=nothing, cache_nearfield::Bool=false,
                           nearfield_cache_max_bytes::Integer=FastMultipole.NEARFIELD_CACHE_DEFAULT_MAX_BYTES,
                           nearfield_cache_max_build_time::Real=Inf)

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
        influence!(body, body, backend; velocity=true, plan_slot, cache_nearfield,
                   nearfield_cache_max_bytes, nearfield_cache_max_build_time)
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
                             plan_slot=nothing, cache_nearfield::Bool=false,
                             nearfield_cache_max_bytes::Integer=FastMultipole.NEARFIELD_CACHE_DEFAULT_MAX_BYTES,
                             nearfield_cache_max_build_time::Real=Inf)

    body.strength[:, 1] .= 0
    strengths_scratch .= x
    _set_strength(body, strengths_scratch)

    body.potential .= 0
    if plan_slot === nothing
        influence!(body, body, backend; scalar_potential=true, velocity=false)
    else
        influence!(body, body, backend; scalar_potential=true, velocity=false,
                   plan_slot, cache_nearfield,
                   nearfield_cache_max_bytes, nearfield_cache_max_build_time)
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
- `:blockgs_normalized_residual` — tuple block-Gauss–Seidel maximum normalized
  Dirichlet-potential/Neumann-normal-velocity residual

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
