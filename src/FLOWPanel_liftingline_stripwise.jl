#=##############################################################################
# DESCRIPTION
    Stripwise elements for non-linear lifting line method

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################

import CSV
import DataFrames: DataFrame
import AirfoilPrep

abstract type StripwiseElement end

################################################################################
# JETFOIL ELEMENT STRUCT
################################################################################
struct Jetfoil{TA<:AbstractVector, TDu<:AbstractVector, TDl<:AbstractVector, 
                TR<:AbstractVector, TM<:AbstractVector, TT<:AbstractVector, 
                Tl<:AbstractArray{<:Number, 6}, 
                Td<:AbstractArray{<:Number, 6}, 
                Tm<:AbstractArray{<:Number, 6}} <: StripwiseElement

    alpha::TA                           # (deg) angle of attack
    deltau::TDu                         # (deg) upper flap deflection
    deltal::TDl                         # (deg) lower flap deflection

    Re::TR                              # Reynolds number
    Mach::TM                            # Mach number
    CT::TT                              # Thrust coefficient

    cl::Tl                              # Lift coefficient
    cd::Td                              # Drag coefficient
    cm::Tm                              # Pitching moment coefficient

end


################################################################################
# AIRFOIL ELEMENT STRUCT
################################################################################
struct Airfoil{TA<:AbstractVector, TR<:AbstractVector, TM<:AbstractVector,
                    Tl<:AbstractArray{<:Number, 3}, 
                    Td<:AbstractArray{<:Number, 3}, 
                    Tm<:AbstractArray{<:Number, 3}} <: StripwiseElement

    alpha::TA                           # (deg) angle of attack

    Re::TR                              # Reynolds number
    Mach::TM                            # Mach number

    cl::Tl                              # Lift coefficient
    cd::Td                              # Drag coefficient
    cm::Tm                              # Pitching moment coefficient

end

################################################################################
# SIMPLE AIRFOIL ELEMENT STRUCT
################################################################################
struct SimpleAirfoil{TA<:AbstractVector,
                        Tl<:AbstractArray{<:Number, 1}, 
                        Td<:AbstractArray{<:Number, 1}, 
                        Tm<:AbstractArray{<:Number, 1},
                        Sl<:math.Akima, Sd<:math.Akima, Sm<:math.Akima
                        } <: StripwiseElement

    alpha::TA                           # (deg) angle of attack

    cl::Tl                              # Lift coefficient
    cd::Td                              # Drag coefficient
    cm::Tm                              # Pitching moment coefficient

    # Pre-computed Akima spline
    spl_cl::Sl
    spl_cd::Sd
    spl_cm::Sm

    function SimpleAirfoil(alpha::TA, cl::Tl, cd::Td, cm::Tm) where {TA, Tl, Td, Tm}

        spl_cl = math.Akima(alpha, cl)
        spl_cd = math.Akima(alpha, cd)
        spl_cm = math.Akima(alpha, cm)

        new{TA, Tl, Td, Tm, typeof(spl_cl), typeof(spl_cd), typeof(spl_cm)}(alpha, cl, cd, cm, spl_cl, spl_cd, spl_cm)
    end

end

function SimpleAirfoil(file_name::String; path::String="")
    data = Matrix(CSV.read(joinpath(path, file_name), DataFrame))
    
    alpha = data[:, 1]
    cl = data[:, 2]
    cd = data[:, 3]
    cm = data[:, 4]

    return SimpleAirfoil(alpha, cl, cd, cm)
end

(self::SimpleAirfoil)(alpha) = (self.spl_cl(alpha), self.spl_cd(alpha), self.spl_cm(alpha))

calc_cl(self::SimpleAirfoil, alpha) = self.spl_cl(alpha)
calc_cd(self::SimpleAirfoil, alpha) = self.spl_cd(alpha)
calc_cm(self::SimpleAirfoil, alpha) = self.spl_cm(alpha)

function extrapolate(self::SimpleAirfoil, args...; optargs...)
    
    alpha, cl, cd, cm = extrapolate(self.alpha*pi/180, self.cl, self.cd, self.cm, 
                                                        args...; optargs...)
    alpha *= 180/pi

    return SimpleAirfoil(alpha, cl, cd, cm)
end

"""
Blend two stripwise elements using a given weight, where `weight=0` simply
returns `airfoil0` and `weight=1` returns `airfoil1`
"""
function blend(airfoil0::SimpleAirfoil, airfoil1::SimpleAirfoil, weight::Number)

    # Convert StripwiseElements to Polar objects
    polar0 = AirfoilPrep.prepy.Polar(-1, airfoil0.alpha, airfoil0.cl, airfoil0.cd, airfoil0.cm)
    polar1 = AirfoilPrep.prepy.Polar(-1, airfoil1.alpha, airfoil1.cl, airfoil1.cd, airfoil1.cm)

    # Blend the two polars
    polar2 = polar0.blend(polar1, weight)

    return SimpleAirfoil(polar2.alpha, polar2.cl, polar2.cd, polar2.cm)
end



##### INTERNAL/COMMON FUNCTIONS ################################################

"""
    viterna(alpha, cl, cd, cr75, nalpha=50)

**COPY/PASTE FROM DuctAPE.jl, modified to also return cm**

Viterna extrapolation.  Follows Viterna paper and somewhat follows NREL version 
of AirfoilPrep, but with some modifications for better robustness and smoothness.

# Arguments
- `alpha::Vector{Float64}`: angles of attack
- `cl::Vector{Float64}`: correspnding lift coefficients
- `cd::Vector{Float64}`: correspnding drag coefficients
- `nalpha::Int64`: number of discrete points (angles of attack) to include in extrapolation

# Returns
- `alpha::Vector{Float64}`: angle of attack from -180 to 180
- `cl::Vector{Float64}`: correspnding extrapolated lift coefficients
- `cd::Vector{Float64}`: correspnding extrapolated drag coefficients
"""
function extrapolate(alpha, cl, cd, cm, AR=5.0, nalpha=50, mincd=0.0001)

    # estimate cdmax
    cdmaxAR = 1.11 + 0.018 * AR
    cdmax = max(maximum(cd), cdmaxAR)

    # find clmax
    i_ps = argmax(cl)  # positive stall
    cl_ps = cl[i_ps]
    cd_ps = cd[i_ps]
    a_ps = alpha[i_ps]

    # and clmin
    i_bs = alpha .< a_ps  # before stall
    i_ns = argmin(cl[i_bs])  # negative stall
    cl_ns = cl[i_bs][i_ns]
    cd_ns = cd[i_bs][i_ns]
    a_ns = alpha[i_bs][i_ns]

    # coefficients in method
    B1pos = cdmax
    A1pos = B1pos / 2.0 * ones(nalpha)
    sa = sin(a_ps)
    ca = cos(a_ps)
    A2pos = (cl_ps - cdmax * sa * ca) * sa / ca^2
    B2pos = (cd_ps - cdmax * sa^2) / ca * ones(nalpha)

    B1neg = cdmax
    A1neg = B1neg / 2.0
    sa = sin(a_ns)
    ca = cos(a_ns)
    A2neg = (cl_ns - cdmax * sa * ca) * sa / ca^2 * ones(nalpha)
    B2neg = (cd_ns - cdmax * sa^2) / ca * ones(nalpha)

    # angles of attack to extrapolate to
    apos = range(alpha[end], pi; length=nalpha + 1)
    apos = apos[2:end]  # don't duplicate point
    aneg = range(-pi, alpha[1]; length=nalpha + 1)
    aneg = aneg[1:(end - 1)]  # don't duplicate point

    # high aoa adjustments
    adjpos = ones(nalpha)
    idx = findall(apos .>= pi / 2)
    adjpos[idx] .= -0.7
    A1pos[idx] .*= -1
    B2pos[idx] .*= -1

    # idx = findall(aneg .<= -alpha[end])

    adjneg = ones(nalpha)
    idx = findall(aneg .<= -pi / 2)
    adjneg[idx] .= 0.7
    A2neg[idx] .*= -1
    B2neg[idx] .*= -1

    # extrapolate
    clpos = @. adjpos * (A1pos * sin(2 * apos) + A2pos * cos(apos)^2 / sin(apos))
    cdpos = @. B1pos * sin(apos)^2 + B2pos * cos(apos)
    clneg = @. adjneg * (A1neg * sin(2 * aneg) + A2neg * cos(aneg)^2 / sin(aneg))
    cdneg = @. B1neg * sin(aneg)^2 + B2neg * cos(aneg)

    # # override region between -alpha_high and alpha_low (if it exists)
    # idx = findall(-alpha[end] .<= aneg .<= alpha[1])
    # @. clneg[idx] = -cl[end]*0.7 + (aneg[idx]+alpha[end])/(alpha[1]+alpha[end])*(cl[1]+cl[end]*0.7)
    # @. cdneg[idx] = cd[1] + (aneg[idx]-alpha[1])/(-alpha[end]-alpha[1])*(cd[end]-cd[1])

    # override with linear variation at ends
    idx = findall(apos .>= pi - a_ps)
    @. clpos[idx] = (apos[idx] - pi) / a_ps * cl_ps * 0.7
    idx = findall(aneg .<= -pi - a_ns)
    @. clneg[idx] = (aneg[idx] + pi) / a_ns * cl_ns * 0.7

    # concatenate
    alphafull = [aneg; alpha; apos]
    clfull = [clneg; cl; clpos]
    cdfull = [cdneg; cd; cdpos]

    # don't allow negative drag
    cdfull = max.(cdfull, mincd)

    # Obtain extrapolated cm from AirfoilPrep
    ap = AirfoilPrep.prepy.Polar(-1, alpha*180/pi, cl, cd, cm)
    ap_extrap = ap.extrapolate(cdmax, nalpha=nalpha)
    cmfull = math.akima(ap_extrap.alpha, ap_extrap.cm, alphafull*180/pi)

    return alphafull, clfull, cdfull, cmfull
end
##### END OF STRIPWISE ELEMENTS ################################################