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

################################################################################
# DIRECT BACKEND
################################################################################

struct DirectBackend <: AbstractBackend
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

#--- overload FastMultipole functions ---#

################################################################################