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

abstract type AbstractBackend end

struct DirectBackend <: AbstractBackend end

struct FastMultipoleBackend <: AbstractBackend
    expansion_order::Int
    multipole_acceptance::Float64
    leaf_size::Int
end

FastMultipoleBackend(; expansion_order=10, multipole_acceptance=0.4, leaf_size=20) = 
    FastMultipoleBackend(expansion_order, multipole_acceptance, leaf_size)

function influence!(targets::Tuple, sources::Tuple, backend::AbstractBackend; optargs...)
    error("influence! not implemented for targets $(typeof(targets)), sources $(typeof(sources)), and backend $(typeof(backend))")
end


#-------- high-level interface -------#

has_semiinfinite_wake(self) = false

function influence!(target_bodies::Tuple, source_bodies::Tuple, backend::FastMultipoleBackend;
                     scalar_potential=false, velocity=false,
                     velocity_gradient=false, precalc=false, optargs...)

    # apply pre-calculations per system
    if precalc
        for target in targets
            pre_evaluate_influence!(target)
        end
    end

    # determine if extra_farfield is needed based
    extra_farfield = false
    for body in source_bodies
        if has_semiinfinite_wake(body)
            extra_farfield = true
            break
        end
    end

    FastMultipole.fmm!(target_bodies, source_bodies;
        expansion_order=backend.expansion_order,
        multipole_acceptance=backend.multipole_acceptance,
        leaf_size_source=backend.leaf_size,
        scalar_potential, gradient=velocity,
        hessian=velocity_gradient,
        extra_farfield,
        shrink=true,
        optargs...)

    return nothing
end

function influence!(target_bodies::Tuple, source_bodies::Tuple, backend::DirectBackend;
                     scalar_potential=false, velocity=false,
                     velocity_gradient=false, precalc=false, optargs...)

    # apply pre-calculations per system
    if precalc
        for target in targets
            pre_evaluate_influence!(target)
        end
    end

    FastMultipole.direct!(target_bodies, source_bodies;
        scalar_potential=scalar_potential,
        gradient=velocity,
        hessian=velocity_gradient, 
        optargs...)

    return nothing
end

to_tuple(val::Tuple) = val
to_tuple(val) = (val,)
function influence!(target, source, backend; optargs...)
    influence!(to_tuple(target), to_tuple(source), backend; optargs...)
end

function pre_evaluate_influence!(system)
    # default behavior
    return nothing
end
