#=##############################################################################
# DESCRIPTION
    FMM abstraction details.
# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Jan 2026
  * License     : GNU Public License
=###############################################################################


################################################################################
# ABSTRACT N-BODY BACKEND
################################################################################

"""
    AbstractBackend

Abstract supertype for influence-evaluation backends used by [`influence!`](@ref).
"""
abstract type AbstractBackend end

"""
    DirectBackend()

Backend that evaluates source-target interactions directly.
"""
struct DirectBackend <: AbstractBackend end

"""
    FastMultipoleBackend
    FastMultipoleBackend(; expansion_order=10, multipole_acceptance=0.4, leaf_size=20)

Backend configuration for accelerated influence evaluation via `FastMultipole`.
"""
struct FastMultipoleBackend <: AbstractBackend
    expansion_order::Int
    multipole_acceptance::Float64
    leaf_size::Int
end

FastMultipoleBackend(; expansion_order=10, multipole_acceptance=0.4, leaf_size=20) =
    FastMultipoleBackend(expansion_order, multipole_acceptance, leaf_size)

"""
    influence!(targets, sources, backend; scalar_potential=false, velocity=false, velocity_gradient=false, precalc=false, postcalc=false, optargs...)

Accumulate influence from `sources` onto `targets` using the requested backend.
Targets and sources may be individual systems or tuples of systems.
"""
function influence!(targets::Tuple, sources::Tuple, backend::AbstractBackend; optargs...)
    error("influence! not implemented for targets $(typeof(targets)), sources $(typeof(sources)), and backend $(typeof(backend))")
end


#-------- high-level interface -------#

has_semiinfinite_wake(self) = false

function influence!(target_bodies::Tuple, source_bodies::Tuple, backend::FastMultipoleBackend;
                     scalar_potential=false, velocity=false,
                     velocity_gradient=false, precalc=false, postcalc=false,
                     plan_slot=nothing, cache_nearfield::Bool=false,
                     nearfield_cache_max_bytes::Integer=FastMultipole.NEARFIELD_CACHE_DEFAULT_MAX_BYTES,
                     nearfield_cache_max_build_time::Real=Inf, optargs...)

    # apply pre-calculations per system
    if precalc
        for target in target_bodies
            pre_evaluate_influence!(target)
        end
    end

    # 051 stage 2 seam: brute-force rectangular influence (host or CUDA) when
    # armed via FLOWPANEL_GPU_INFLUENCE; a no-op returning false otherwise.
    # Handles only the wake-influence and post-solve body-influence passes;
    # the panel solve and anything unrecognized fall through to fmm! below.
    # See src/FLOWPanel_gpu_influence.jl.
    if _gpu_rect_influence!(target_bodies, source_bodies;
            scalar_potential, velocity, velocity_gradient, optargs...)
        postcalc && post_evaluate_influence!(target_bodies, source_bodies, backend, nothing)
        return nothing
    end

    # determine if extra_farfield is needed based
    extra_farfield = false
    for body in source_bodies
        if has_semiinfinite_wake(body)
            extra_farfield = true
            break
        end
    end

    # plan_slot (a Ref): reuse trees/lists/buffers across calls with frozen
    # geometry — only source STRENGTHS may change between calls (the FmmPlan
    # contract; the buffer radii fold in the ACTIVE core_size, so any
    # offset flip requires the owner to clear the slot). The slot stores
    # (plan, key); a changed derivative request rebuilds rather than trusting
    # a mismatched plan.
    outputs = if plan_slot === nothing
        FastMultipole.fmm!(target_bodies, source_bodies;
            expansion_order=backend.expansion_order,
            multipole_acceptance=backend.multipole_acceptance,
            leaf_size_source=backend.leaf_size,
            scalar_potential, gradient=velocity,
            hessian=velocity_gradient,
            extra_farfield,
            shrink=true,
            optargs...)
    else
        key = (scalar_potential, velocity, velocity_gradient, extra_farfield)
        if plan_slot[] === nothing || plan_slot[][2] != key
            plan = FastMultipole.FmmPlan(target_bodies, source_bodies;
                expansion_order=backend.expansion_order,
                multipole_acceptance=backend.multipole_acceptance,
                leaf_size_source=backend.leaf_size,
                scalar_potential, gradient=velocity,
                hessian=velocity_gradient,
                shrink=true)
            # cache_nearfield: freeze the near field as dense blocks owned by
            # the fresh plan (same lifetime/validity contract); fmm! picks the
            # cache up from the plan automatically. The caps are forwarded so
            # callers can size them to the case: FastMultipole's 4 GiB default
            # is below what the larger 021 ladder rungs need (R4 ≈ 4.5 GiB),
            # and the build is serial, so max_build_time is wall-clock.
            cache_nearfield && FastMultipole.build_nearfield_cache!(plan,
                target_bodies, source_bodies;
                max_bytes=nearfield_cache_max_bytes,
                max_build_time=nearfield_cache_max_build_time)
            plan_slot[] = (plan, key)
        end
        FastMultipole.fmm!(target_bodies, source_bodies, plan_slot[][1];
            extra_farfield, optargs...)
    end

    if postcalc
        post_evaluate_influence!(target_bodies, source_bodies, backend, outputs)
    end

    return nothing
end

function influence!(target_bodies::Tuple, source_bodies::Tuple, backend::DirectBackend;
                     scalar_potential=false, velocity=false,
                     velocity_gradient=false, precalc=false, postcalc=false, optargs...)

    # apply pre-calculations per system
    if precalc
        for target in target_bodies
            pre_evaluate_influence!(target)
        end
    end

    FastMultipole.direct!(target_bodies, source_bodies;
        scalar_potential=scalar_potential,
        gradient=velocity,
        hessian=velocity_gradient, 
        optargs...)

    if postcalc
        post_evaluate_influence!(target_bodies, source_bodies, backend, nothing)
    end

    return nothing
end

to_tuple(val::Tuple) = val
to_tuple(val) = (val,)
function influence!(target, source, backend; optargs...)
    influence!(to_tuple(target), to_tuple(source), backend; optargs...)
end

"""
    pre_evaluate_influence!(system)

Hook for systems that need preprocessing before influence evaluation. The
default implementation does nothing.
"""
function pre_evaluate_influence!(system)
    # default behavior
    return nothing
end

"""
    post_evaluate_influence!(targets, sources, backend, outputs)

Hook for systems that need postprocessing after influence evaluation. Tuple
dispatch walks all target-source combinations; the default scalar
implementation does nothing.
"""
function post_evaluate_influence!(targets::Tuple, sources::Tuple, backend::AbstractBackend, outputs)
    for (i_t, target) in enumerate(targets)
        for (i_s, source) in enumerate(sources)
            post_evaluate_influence!(target, source, backend, outputs; i_target=i_t, i_source=i_s)
        end
    end
    return nothing
end

function post_evaluate_influence!(target, source, backend::AbstractBackend, outputs; i_target::Int=1, i_source::Int=1)
    # default behavior
    return nothing
end
