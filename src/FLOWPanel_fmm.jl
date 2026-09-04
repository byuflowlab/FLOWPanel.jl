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
                     production_route=nothing,
                     nearfield_cache_max_bytes::Integer=FastMultipole.NEARFIELD_CACHE_DEFAULT_MAX_BYTES,
                     nearfield_cache_max_build_time::Real=Inf,
                     sfs_stage::Tuple=(1, 1), optargs...)

    # apply pre-calculations per system
    if precalc
        for target in target_bodies
            pre_evaluate_influence!(target; sfs_stage)
        end
    end

    # 051 stage 2 seam: brute-force rectangular influence (host or CUDA) when
    # armed via FLOWPANEL_GPU_INFLUENCE; a no-op returning false otherwise.
    # Handles only the wake-influence and post-solve body-influence passes;
    # the panel solve and anything unrecognized fall through to fmm! below.
    # See src/FLOWPanel_gpu_influence.jl.
    if _gpu_rect_influence!(target_bodies, source_bodies;
            scalar_potential, velocity, velocity_gradient, production_route,
            optargs...)
        postcalc && post_evaluate_influence!(target_bodies, source_bodies, backend, nothing; sfs_stage)
        return nothing
    end

    # Passes the seam does not handle run the host fmm! path below, which
    # scalar-indexes particle storage — disallowed on a device-backed
    # particle field. Swap such fields for their refreshed host mirrors
    # (the mirror is cached per field, so FmmPlan reuse and pfield===source
    # identity checks still hold across calls) and push the accumulated
    # outputs back to any device field among the targets afterwards.
    original_targets = target_bodies
    target_bodies = map(_influence_host_system, target_bodies)
    source_bodies = map(_influence_host_system, source_bodies)

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
        post_evaluate_influence!(target_bodies, source_bodies, backend, outputs; sfs_stage)
    end

    _influence_push_back!(original_targets, target_bodies)

    return nothing
end

# Device-wake seam for the non-rectangular influence paths: host view of a
# device-backed particle field (identity for every other system), plus the
# corresponding target write-back.
_influence_host_system(x) = x
function _influence_host_system(pf::FLOWVPM.ParticleField)
    _gpu_pfield_on_device(pf) || return pf
    return _gpu_sync_mirror_from_device!(_gpu_pfield_mirror(pf), pf)
end
function _influence_push_back!(originals::Tuple, hosted::Tuple)
    for (orig, host) in zip(originals, hosted)
        host === orig || _gpu_sync_device_from_mirror!(orig, host)
    end
    return nothing
end

function influence!(target_bodies::Tuple, source_bodies::Tuple, backend::DirectBackend;
                     scalar_potential=false, velocity=false,
                     velocity_gradient=false, precalc=false, postcalc=false,
                     production_route=nothing, sfs_stage::Tuple=(1, 1), optargs...)

    # apply pre-calculations per system
    if precalc
        for target in target_bodies
            pre_evaluate_influence!(target; sfs_stage)
        end
    end

    # same device-wake seam as the FastMultipoleBackend path: direct! also
    # scalar-indexes particle storage on the host
    original_targets = target_bodies
    target_bodies = map(_influence_host_system, target_bodies)
    source_bodies = map(_influence_host_system, source_bodies)

    FastMultipole.direct!(target_bodies, source_bodies;
        scalar_potential=scalar_potential,
        gradient=velocity,
        hessian=velocity_gradient,
        optargs...)

    if postcalc
        post_evaluate_influence!(target_bodies, source_bodies, backend, nothing; sfs_stage)
    end

    _influence_push_back!(original_targets, target_bodies)

    return nothing
end

to_tuple(val::Tuple) = val
to_tuple(val) = (val,)
function influence!(target, source, backend; optargs...)
    influence!(to_tuple(target), to_tuple(source), backend; optargs...)
end

"""
    pre_evaluate_influence!(system; sfs_stage=(1, 1))

Hook for systems that need preprocessing before influence evaluation. The
default implementation does nothing. `sfs_stage` carries the RK stage weights
`(a, b)` to SFS hooks that gate on the stage (FLOWVPM's Dynamic/Constant SFS
run their procedures only at `a == 1 || a == 0`); the default `(1, 1)`
preserves the single-evaluation-per-step behavior.
"""
function pre_evaluate_influence!(system; sfs_stage::Tuple=(1, 1))
    # default behavior
    return nothing
end

"""
    post_evaluate_influence!(targets, sources, backend, outputs)

Hook for systems that need postprocessing after influence evaluation. Tuple
dispatch walks all target-source combinations; the default scalar
implementation does nothing.
"""
function post_evaluate_influence!(targets::Tuple, sources::Tuple, backend::AbstractBackend, outputs;
        sfs_stage::Tuple=(1, 1))
    for (i_t, target) in enumerate(targets)
        for (i_s, source) in enumerate(sources)
            post_evaluate_influence!(target, source, backend, outputs; i_target=i_t, i_source=i_s, sfs_stage)
        end
    end
    return nothing
end

function post_evaluate_influence!(target, source, backend::AbstractBackend, outputs;
        i_target::Int=1, i_source::Int=1, sfs_stage::Tuple=(1, 1))
    # default behavior
    return nothing
end
