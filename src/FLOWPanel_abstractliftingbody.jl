#=##############################################################################
# DESCRIPTION
    Lifting paneled body types definition.
# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Sep 2018
  * License     : MIT License
=###############################################################################

################################################################################
# ABSTRACT LIFTING BODY BODY TYPE
################################################################################
"""
    AbstractLiftingBody

Abstract supertype for lifting-surface body models with shedding-edge and
wake-direction data in addition to the fields required by [`AbstractBody`](@ref).
"""
abstract type AbstractLiftingBody{E, N, TF, DBC} <: AbstractBody{E, N, TF, DBC} end

"""
    solve(body::AbstractLiftingBody, Uinfs, Das)

Solve a lifting-body boundary-value problem using control-point velocities
`Uinfs` and wake-direction data `Das`.
"""
function solve(self::AbstractLiftingBody, Uinfs::AbstractMatrix,
               Das::AbstractMatrix)
    error("solve(...) for body type $(typeof(self)) has not been implemented yet!")
end

##### COMMON FUNCTIONS  ########################################################

##### COMMON INTERNAL FUNCTIONS  ###############################################

function extra_reset!(body::AbstractLiftingBody)
    for vte in body.velocity_te
        vte .= 0.0
    end
end

function extra_apply_freestream!(body::AbstractLiftingBody, uinf)
    for i in eachindex(body.velocity_te)
        eachcol(body.velocity_te[i]) .+= Ref(uinf)
    end
end

##### END OF ABSTRACT LIFTING BODY #############################################
