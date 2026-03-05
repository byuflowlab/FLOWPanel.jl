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

function evaluate_influence!(targets, sources, backend::AbstractBackend)
    error("evaluate_influence! not implemented for backend $(typeof(backend))")
end

################################################################################
# DIRECT BACKEND
################################################################################

struct DirectBackend <: AbstractBackend
end

function evaluate_influence!(targets::Tuple, sources::Tuple, ::DirectBackend; 
        scalar_potential=true, gradient=true, hessian=false
    )
    FastMultipole.direct!(targets, sources; scalar_potential, gradient, hessian)
end

################################################################################
# FASTMULTIPOLE BACKEND
################################################################################

struct FastMultipoleBackend <: AbstractBackend
    expansion_order::Int
    multipole_acceptance::Float64
    leaf_size::Int
end

function FastMultipoleBackend(; expansion_order::Int=5,
                                    multipole_acceptance::Float64=0.5,
                                    leaf_size::Int=10)
    return FastMultipoleBackend(expansion_order,
                                multipole_acceptance,
                                leaf_size)
end

function evaluate_influence!(targets::Tuple, sources::Tuple, backend::FastMultipoleBackend; 
        scalar_potential=true, gradient=true, hessian=false
    )
    # unpack backend
    expansion_order = backend.expansion_order
    multipole_acceptance = backend.multipole_acceptance
    leaf_size = backend.leaf_size

    # call FMM
    FastMultipole.fmm!(targets, sources; 
        expansion_order, multipole_acceptance, leaf_size_source=leaf_size,
        scalar_potential, gradient, hessian
    )
end

#--- overload FastMultipole functions ---#

################################################################################