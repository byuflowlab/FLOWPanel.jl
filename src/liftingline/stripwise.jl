#=##############################################################################
# DESCRIPTION
    Stripwise elements for non-linear lifting line method

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################

abstract type StripwiseElement <: AbstractElement end


################################################################################
# AIRFOIL ELEMENT STRUCT
################################################################################
struct GeneralAirfoil{N,
                        Tdim<:NTuple{N, Int},
                        Tp<:NTuple{N, <:AbstractVector},
                        Tn<:NTuple{N, <:AbstractString},
                        Tl<:AbstractArray{<:Number, N}, 
                        Td<:AbstractArray{<:Number, N}, 
                        Tm<:AbstractArray{<:Number, N},
                        } <: StripwiseElement

    dims::Tdim                          # Number of grid nodes in each dimension
    parameters::Tp                      # Axes of uniform grid where `parameters[i]` contains the values of the i-th dimension
    names::Tn                           # Name of each parameter dimension

    cl::Tl                              # Lift coefficient grid data
    cd::Td                              # Drag coefficient grid data
    cm::Tm                              # Pitching moment coefficient grid data

    interp1d::Function                  # Base 1D interpolator function

    function GeneralAirfoil(parameters::NTuple{N, <:AbstractVector},
                            names::NTuple{N, <:AbstractString},
                            cl::AbstractArray{<:Number, N}, 
                            cd::AbstractArray{<:Number, N}, 
                            cm::AbstractArray{<:Number, N};
                            interp1d=math.akima) where {N}

        dims = Tuple(length(values) for values in parameters)
                            
        return new{N, 
                    typeof(dims), typeof(parameters), typeof(names), 
                    typeof(cl), typeof(cd), typeof(cm)}(
                        dims, parameters, names,
                        cl, cd, cm, 
                        interp1d
                    )

    end

end

"""

## Arguments
* `sweep_name`          :   Sweep name (folder to read CSV files from)
* `path`                :   Where to read sweep from (rest of path to sweep)
* `cl_name`             :   Name of cl column
* `cd_name`             :   Name of cd column
* `cm_name`             :   Name of cm column

NOTE: Here we assume that the first column in each CSV file corresponds to the 
angle of attack (deg).
"""
function GeneralAirfoil(sweep_name::String; path::String="",
                            cl_name = "cl",
                            cd_name = "cd",
                            cm_name = "cm",
                            file_extension = ".csv",
                            verbose=true, tab_lvl=0, optargs...)

    sweep_path = joinpath(path, sweep_name)                 # Where to read sweep from
    data_names = (cl_name, cd_name, cm_name)                # Names of data columns
    ndata = length(data_names)                              # Number of data columns

    # Read sweep
    files = [f for f in readdir(sweep_path) if splitext(f)[end]==file_extension]

    df = DataFrame()
    for f in files
        data = CSV.read(joinpath(sweep_path, f), DataFrame)

        append!(df, data)
    end

    # Determine parameters                                  # Names of parameter columns
    parameter_names = [n for n in names(df) if !(n in data_names)]
    nparameters = length(parameter_names)                   # Number of parameters

    # Fetch table of parameter values
    df_parameters = df[!, parameter_names]

    # Reduce parameters to ranges of full coverage across the dataset
    # NOTE: this is needed to ensure the data is on a uniform grid
    parameters = reduce_parameters(df_parameters)

    # Display parameter ranges for user to verify
    if verbose

        println("\t"^tab_lvl * "Read sweep under $(sweep_path) and found "*
                "the following parameter ranges of full coverage")

        for (name, vals) in zip(parameter_names, parameters)
            println("\t"^(tab_lvl+1) * "$(name): $(vals)")
        end

    end

    # Discard any single-valued parameters
    to_discard = [i for (i, (name, vals)) in enumerate(zip(parameter_names, parameters)) if length(vals) <= 1]

    if verbose

        println("\t"^tab_lvl * "Discarding the following parameters due to"*
                " being single-valued: " * join((parameter_names[i] for i in to_discard), ", "))

    end

    reverse!(to_discard)

    for i in to_discard
        deleteat!(parameter_names, i)
        parameters = tuple(deleteat!(collect(parameters), i)...)
    end

    # Format list of parameters into a matrix
    dims = Tuple(length(values) for values in parameters)
    parameters_matrix = Array{NTuple{length(dims), Float64}, length(dims)}(undef, dims...)

    for I in CartesianIndices(dims)
        parameters_matrix[I] = Tuple(parameters[di][i] for (di, i) in enumerate(Tuple(I)))
    end

    # Fetch data as a uniform grid
    cl = [df[broadcast(&, [getproperty(df, Symbol(name)) .== val for (name, val) in zip(parameter_names, entry)]...), cl_name][1] for entry in parameters_matrix]
    cd = [df[broadcast(&, [getproperty(df, Symbol(name)) .== val for (name, val) in zip(parameter_names, entry)]...), cd_name][1] for entry in parameters_matrix]
    cm = [df[broadcast(&, [getproperty(df, Symbol(name)) .== val for (name, val) in zip(parameter_names, entry)]...), cm_name][1] for entry in parameters_matrix]

    # Format data as a uniform grid
    
    # cl = [df[I, cl_name] for I in LinearIndices(dims)]
    # cl = reshape(cl, dims)

    # cd = [df[I, cd_name] for I in LinearIndices(dims)]
    # cd = reshape(cd, dims)

    # cm = [df[I, cm_name] for I in LinearIndices(dims)]
    # cm = reshape(cm, dims)


    return GeneralAirfoil(parameters, tuple(parameter_names...), cl, cd, cm; optargs...)
end



function Base.show(io::IO, self::GeneralAirfoil{N}) where {N}

    println(io, "GeneralAirfoil with $(N) dimensional parameters")

    for (i, (name, vals, dim)) in enumerate(zip(self.names, self.parameters, self.dims))
        if i != N
            print(io, "├─")
        else
            print(io, "└─")
        end
        println(io, rpad(name, 20, " ") * lpad(dim, 3, " ") * " values" * " [$(vals[1]), ..., $(vals[end])]")
    end

end


(self::GeneralAirfoil)(args...; optargs...) = (
    calc_cl(self, args...; optargs...), 
    calc_cd(self, args...; optargs...), 
    calc_cm(self, args...; optargs...)
)

calc_claero(self::GeneralAirfoil, args...; optargs...) = calc_cl(self, args...; optargs...)

function calc_cl(self::GeneralAirfoil{1}, args...; optargs...)
    self.interp1d(self.parameters..., self.cl, args...; optargs...)
end
function calc_cd(self::GeneralAirfoil{1}, args...; optargs...)
    self.interp1d(self.parameters..., self.cd, args...; optargs...)
end
function calc_cm(self::GeneralAirfoil{1}, args...; optargs...)
    self.interp1d(self.parameters..., self.cm, args...; optargs...)
end

function calc_cl(self::GeneralAirfoil{2}, args...; optargs...)
    math.interp2d(self.interp1d, self.parameters..., self.cl, args...; optargs...)[1]
end
function calc_cd(self::GeneralAirfoil{2}, args...; optargs...)
    math.interp2d(self.interp1d, self.parameters..., self.cd, args...; optargs...)[1]
end
function calc_cm(self::GeneralAirfoil{2}, args...; optargs...)
    math.interp2d(self.interp1d, self.parameters..., self.cm, args...; optargs...)[1]
end

function calc_cl(self::GeneralAirfoil{3}, args...; optargs...)
    math.interp3d(self.interp1d, self.parameters..., self.cl, args...; optargs...)[1]
end
function calc_cd(self::GeneralAirfoil{3}, args...; optargs...)
    math.interp3d(self.interp1d, self.parameters..., self.cd, args...; optargs...)[1]
end
function calc_cm(self::GeneralAirfoil{3}, args...; optargs...)
    math.interp3d(self.interp1d, self.parameters..., self.cm, args...; optargs...)[1]
end

function calc_cl(self::GeneralAirfoil{4}, args...; optargs...)
    math.interp4d(self.interp1d, self.parameters..., self.cl, args...; optargs...)[1]
end
function calc_cd(self::GeneralAirfoil{4}, args...; optargs...)
    math.interp4d(self.interp1d, self.parameters..., self.cd, args...; optargs...)[1]
end
function calc_cm(self::GeneralAirfoil{4}, args...; optargs...)
    math.interp4d(self.interp1d, self.parameters..., self.cm, args...; optargs...)[1]
end

function calc_cl(self::GeneralAirfoil{5}, args...; optargs...)
    interp5d(self.interp1d, self.parameters..., self.cl, args...; optargs...)[1]
end
function calc_cd(self::GeneralAirfoil{5}, args...; optargs...)
    interp5d(self.interp1d, self.parameters..., self.cd, args...; optargs...)[1]
end
function calc_cm(self::GeneralAirfoil{5}, args...; optargs...)
    interp4d(self.interp1d, self.parameters..., self.cm, args...; optargs...)[1]
end

function extrapolate(self::GeneralAirfoil, args...; optargs...)

    new_parameters, new_cl, new_cd, new_cm = extrapolate_ndim(self.parameters,
                                                                self.cl,
                                                                self.cd,
                                                                self.cm,
                                                                args...; 
                                                                optargs...)
    return GeneralAirfoil(new_parameters, self.names,
                            new_cl, new_cd, new_cm;
                            interp1d=self.interp1d)
end

"""
Apply Viterna extrapolation to N-dimensional lookup tables
"""
function extrapolate_ndim(parameters::NTuple{N, <:AbstractVector}, 
                            cl::AbstractArray{R1, N},
                            cd::AbstractArray{R2, N},
                            cm::AbstractArray{R3, N}, 
                            args...; 
                            optargs...) where {N, R1, R2, R3}

    R = promote_type(R1, R2, R3)

    # Base case
    if N == 1

        alpha, cl, cd, cm = extrapolate(parameters[1]*pi/180, cl, cd, cm, args...; optargs...)

        alpha *= 180/pi

        return (alpha, ), cl, cd, cm

    # Recursively extrapolate cl, cd, and cm curves along inner dimensions
    else

        dims = Tuple(length(values) for values in parameters)

        # Slice the parameters along the outer inner dimension (d==2)
        sliced_parameters = parameters[[i for i in eachindex(parameters) if i != 2]]

        # Iterate over the values of the outer inner dimension extrapolating
        new_alphas = nothing
        new_cl, new_cd, new_cm = nothing, nothing, nothing

        for (i, val) in enumerate(parameters[2])

            # Extrapolate this slice
            (this_parameters,
                    this_cl, this_cd, this_cm) = extrapolate_ndim(
                                                            sliced_parameters, 
                                                            selectdim(cl, 2, i), 
                                                            selectdim(cd, 2, i), 
                                                            selectdim(cm, 2, i), 
                                                            args...; optargs...)

            # Initialize storage memory
            if i==1
                new_alphas = this_parameters[1]
                new_cl = zeros(R, length(new_alphas), dims[2:end]...)
                new_cd = zeros(R, length(new_alphas), dims[2:end]...)
                new_cm = zeros(R, length(new_alphas), dims[2:end]...)
            end

            # Store extrapolation at this slice
            selectdim(new_cl, 2, i) .= this_cl
            selectdim(new_cd, 2, i) .= this_cd
            selectdim(new_cm, 2, i) .= this_cm

        end

        new_parameters = (new_alphas, parameters[2:end]...)

        return new_parameters, new_cl, new_cd, new_cm

    end

end

"""
Blend two stripwise elements using a given weight, where `weight=0` simply
returns `airfoil0` and `weight=1` returns `airfoil1`
"""
function blend(airfoil0::GeneralAirfoil, airfoil1::GeneralAirfoil, weight::Number)

    alphas, cls, cds, cms = blend(airfoil0.parameters[1], 
                                        airfoil0.cl, airfoil0.cd, airfoil0.cm, 
                                        airfoil1.parameters[1], 
                                        airfoil1.cl, airfoil1.cd, airfoil1.cm, 
                                        weight)

    return SimpleAirfoil(alphas, cls, cds, cms)
end

function blend_ndim(parameters1::NTuple{N, <:AbstractVector}, 
                            cl1::AbstractArray{R1, N},
                            cd1::AbstractArray{R2, N},
                            cm1::AbstractArray{R3, N},
                    parameters2::NTuple{N, <:AbstractVector}, 
                            cl2::AbstractArray{R4, N},
                            cd2::AbstractArray{R5, N},
                            cm2::AbstractArray{R6, N},
                            args...; optargs...) where {N, R1, R2, R3, R4, R5, R6}

    
    R = promote_type(R1, R2, R3, R4, R5, R6)

    # Base case
    if N == 1

        alpha, cl, cd, cm = blend(parameters1[1], cl1, cd1, cm1, 
                                    parameters2[1], cl2, cd2, cm2, 
                                    args...; optargs...)

        return (alpha, ), cl, cd, cm

    # Recursively blend cl, cd, and cm curves along inner dimensions
    else

        dims = Tuple(length(values) for values in parameters)

        # Slice the parameters along the outer inner dimension (d==2)
        sliced_parameters = parameters[[i for i in eachindex(parameters) if i != 2]]

        # Iterate over the values of the outer inner dimension blending it
        new_alphas = nothing
        new_cl, new_cd, new_cm = nothing, nothing, nothing

        for (i, val) in enumerate(parameters[2])

            # blend this slice
            (this_parameters,
                    this_cl, this_cd, this_cm) = extrapolate_ndim(
                                                            sliced_parameters, 
                                                            selectdim(cl, 2, i), 
                                                            selectdim(cd, 2, i), 
                                                            selectdim(cm, 2, i), 
                                                            args...; optargs...)

            # Initialize storage memory
            if i==1
                new_alphas = this_parameters[1]
                new_cl = zeros(R, length(new_alphas), dims[2:end]...)
                new_cd = zeros(R, length(new_alphas), dims[2:end]...)
                new_cm = zeros(R, length(new_alphas), dims[2:end]...)
            end

            # Store extrapolation at this slice
            selectdim(new_cl, 2, i) .= this_cl
            selectdim(new_cd, 2, i) .= this_cd
            selectdim(new_cm, 2, i) .= this_cm

        end

        new_parameters = (new_alphas, parameters[2:end]...)

        return new_parameters, new_cl, new_cd, new_cm

    end

end


function plot_slice(self::GeneralAirfoil{N}, slice; 
                        alphas = range(-180, 180, step=1)
                        ) where N

    @assert length(slice) == N-1 ""*
        "Invalid slice dimensions; expected $(N-1) dimension, got $(length(slice))"

    for (si, i) in enumerate(slice)
        @assert 0 < i <= self.dims[si+1] ""*
            "Slice $(slice) out of bounds $(self.dims[2:end])"
    end

    fig = plt.figure(figsize = [7, 0.75*5*3]*7/9)
    axs = fig.subplots(3, 1)

    fig.suptitle("Data slice at " * 
                join(("$(name) = $(self.parameters[ni][slice[ni-1]])" for (ni, name) in enumerate(self.names) if ni != 1), ", "))

    slice_parameters = (self.parameters[si+1][i] for (si, i) in enumerate(slice))

    ax = axs[1]
    ax.plot(self.parameters[1], self.cl[:, slice...], "ok", label="Data")
    ax.plot(alphas, [calc_cl(self, a, slice_parameters...) for a in alphas], "-", label="N-dimensional spline")

    ax.set_ylabel(L"c_\ell")

    ax = axs[2]
    ax.plot(self.parameters[1], self.cd[:, slice...], "ok", label="Data")
    ax.plot(alphas, [calc_cd(self, a, slice_parameters...) for a in alphas], "-", label="N-dimensional spline")

    ax.set_ylabel(L"c_d")

    ax = axs[3]
    ax.plot(self.parameters[1], self.cm[:, slice...], "ok", label="Data")
    ax.plot(alphas, [calc_cm(self, a, slice_parameters...) for a in alphas], "-", label="N-dimensional spline")

    ax.set_ylabel(L"c_m")


    for ax in axs
        ax.set_xlabel(self.names[1])
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        ax.legend(loc="best", frameon=false, fontsize=8)
    end
        
    fig.tight_layout()

    return fig, axs
end

################################################################################
# SIMPLE AIRFOIL ELEMENT STRUCT
################################################################################
struct SimpleAirfoil{Ta<:AbstractVector,
                        Tl<:AbstractArray{<:Number, 1}, 
                        Td<:AbstractArray{<:Number, 1}, 
                        Tm<:AbstractArray{<:Number, 1},
                        R<:Number,
                        Sl<:math.Akima, Sd<:math.Akima, Sm<:math.Akima
                        } <: StripwiseElement

    alpha::Ta                           # (deg) angle of attack

    cl::Tl                              # Lift coefficient
    cd::Td                              # Drag coefficient
    cm::Tm                              # Pitching moment coefficient

    alpha0::R                           # (deg) AOA at zero lift

    # Pre-computed Akima spline
    spl_cl::Sl
    spl_cd::Sd
    spl_cm::Sm

    function SimpleAirfoil(alpha::Ta, cl::Tl, cd::Td, cm::Tm) where {Ta, Tl, Td, Tm}

        # Spline data
        spl_cl = math.Akima(alpha, cl)
        spl_cd = math.Akima(alpha, cd)
        spl_cm = math.Akima(alpha, cm)

        # Find AOA at zero lift
        f(u, p) = [spl_cl(u[1])]
        u0 = [0.0]
        prob = SimpleNonlinearSolve.NonlinearProblem{false}(f, u0)
        result = SimpleNonlinearSolve.solve(prob, SimpleNonlinearSolve.SimpleNewtonRaphson(), abstol = 1e-9)
        alpha0 = result.u[1]

        new{Ta, Tl, Td, Tm, 
            typeof(alpha0),
            typeof(spl_cl), typeof(spl_cd), typeof(spl_cm)
            }(alpha, cl, cd, cm, alpha0, spl_cl, spl_cd, spl_cm)
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

calc_claero(self::SimpleAirfoil, alpha) = self.spl_cl(alpha)
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

    alphas, cls, cds, cms = blend(airfoil0.alpha, 
                                        airfoil0.cl, airfoil0.cd, airfoil0.cm, 
                                        airfoil1.alpha, 
                                        airfoil1.cl, airfoil1.cd, airfoil1.cm, 
                                        weight)

    return SimpleAirfoil(alphas, cls, cds, cms)
end



##### INTERNAL/COMMON FUNCTIONS ################################################

"""
Calculate swept sectional lift coefficient as in Goates 2022, Eq. (28).
"""
function calc_sweptcl(airfoil::StripwiseElement, sweep::Number, alpha_Λ::Number, 
                                                            args...; optargs...)

    alpha = alpha_Λ + airfoil.alpha0*( 1 - 1/cosd(sweep) )

    return calc_cl(airfoil, alpha, args...; optargs...)

end

"""
Calculate swept sectional pitching moment coefficient as in Goates 2022, Eq. 
(30).
"""
function calc_sweptcm(airfoil::StripwiseElement, sweep::Number, alpha_Λ::Number, 
                                                            args...; optargs...)

    return calc_cm(airfoil, alpha_Λ, args...; optargs...) / cosd(sweep)

end

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

    # --- Begin Julia replacement for AirfoilPrep cm extrapolation ---

    # Build alpha_cm and cm_ext (degrees)
    cm1_alpha = floor(alpha[1]*180.0/pi / 10.0) * 10.0
    cm2_alpha = ceil(alpha[end]*180.0/pi  / 10.0) * 10.0

    alpha_num = abs(Int((-180.0 - cm1_alpha) / 10.0 - 1))
    alpha_cm1 = collect(range(-180.0, stop=cm1_alpha, length=max(alpha_num, 0)))
    alpha_cm2_len = Int((180.0 - cm2_alpha)/10.0 + 1)
    alpha_cm2 = collect(range(cm2_alpha, stop=180.0, length=max(alpha_cm2_len, 0)))

    # alpha_cm in degrees, include original alpha converted to degrees
    alpha_cm = vcat(alpha_cm1, (alpha .* 180.0/pi), alpha_cm2)

    cm1 = zeros(length(alpha_cm1))
    cm2 = zeros(length(alpha_cm2))
    cm_ext = vcat(cm1, cm, cm2)

    cl_high = cl[end]
    cd_high = cd[end]
    cm_high = cm[end]

    # Only attempt extrapolation if there are any non-zero cm entries
    if count(x -> x != 0, cm) > 0
        # compute cmCoef and cm0
        cmCoef, cm0 = _cm_coeff(alpha .* 180.0/pi, cl, cm, cl_high, cd_high, cm_high)

        # Interpolate cl and cd onto the alpha_cm grid (degrees)
        cl_cm = math.linear(alpha .* 180.0/pi, cl, alpha_cm)
        cd_cm = math.linear(alpha .* 180.0/pi, cd, alpha_cm)

        alpha_low_deg = alpha[1] * 180.0/pi
        alpha_high_deg = alpha[end] * 180.0/pi

        for i in 1:length(alpha_cm)
            cm_new = _get_cm(i, cmCoef, alpha_cm, cl_cm, cd_cm, alpha_low_deg, alpha_high_deg, cm0)
            if cm_new === nothing
                # skip (mirrors Python 'None' handling)
            else
                cm_ext[i] = cm_new
            end
        end
    end

    # Interpolate cm_ext onto alphafull (alphafull is radians); convert alphafull to degrees
    cmfull = nothing
    try
        cmfull = math.akima(alpha_cm, cm_ext, alphafull .* 180.0/pi)
    catch
        cmfull = zeros(length(clfull))
    end
    # --- End Julia replacement ---

    return alphafull, clfull, cdfull, cmfull
end

function blend(alphas1::AbstractArray, 
                cls1::AbstractArray, cds1::AbstractArray, cms1::AbstractArray, 
                alphas2::AbstractArray, 
                cls2::AbstractArray, cds2::AbstractArray, cms2::AbstractArray, 
                weight::Number)

    # Determine bounds of minimum set
    alpha_min = max(minimum(alphas1), minimum(alphas2))
    alpha_max = min(maximum(alphas1), maximum(alphas2))

    # Reduce alphas to minimum unique set
    new_alphas = unique!(sort!(vcat([[a for a in alphas if alpha_min <= a && a<=alpha_max] for alphas in (alphas1, alphas2)]...)))

    # Interpolate data to new alphas
    new_cls1 = math.linear(alphas1, cls1, new_alphas)
    new_cds1 = math.linear(alphas1, cds1, new_alphas)
    new_cms1 = math.linear(alphas1, cms1, new_alphas)

    new_cls2 = math.linear(alphas2, cls2, new_alphas)
    new_cds2 = math.linear(alphas2, cds2, new_alphas)
    new_cms2 = math.linear(alphas2, cms2, new_alphas)

    # Linearly interpolate data
    new_cls = new_cls1 + weight*(new_cls2 - new_cls1)
    new_cds = new_cds1 + weight*(new_cds2 - new_cds1)
    new_cms = new_cms1 + weight*(new_cms2 - new_cms1)

    return new_alphas, new_cls, new_cds, new_cms
end



# Helper: compute cm value for a single alpha index (port of Python __getCM)
function _get_cm(i::Integer, cmCoef::Number, alpha_cm::AbstractVector,
                cl_ext::AbstractVector, cd_ext::AbstractVector,
                alpha_low_deg::Number, alpha_high_deg::Number,
                cm0::Number)

    a = alpha_cm[i]

    # If inside the original cm-provided range -> nothing to do
    if a >= alpha_low_deg && a <= alpha_high_deg
        return nothing
    end

    # Typical region
    if a > -165.0 && a < 165.0
        if abs(a) < 0.01
            return cm0
        else
            if a > 0.0
                x = cmCoef * tan((a * pi/180.0) - pi/2) + 0.25
                return cm0 - x * (cl_ext[i] * cos(a * pi/180.0) + cd_ext[i] * sin(a * pi/180.0))
            else
                x = cmCoef * tan((-a * pi/180.0) - pi/2) + 0.25
                return -(cm0 - x * (-cl_ext[i] * cos(-a * pi/180.0) + cd_ext[i] * sin(-a * pi/180.0)))
            end
        end
    else
        if isapprox(a, 165.0; atol=1e-8)
            return -0.4
        elseif isapprox(a, 170.0; atol=1e-8)
            return -0.5
        elseif isapprox(a, 175.0; atol=1e-8)
            return -0.25
        elseif isapprox(a, 180.0; atol=1e-8)
            return 0.0
        elseif isapprox(a, -165.0; atol=1e-8)
            return 0.35
        elseif isapprox(a, -170.0; atol=1e-8)
            return 0.4
        elseif isapprox(a, -175.0; atol=1e-8)
            return 0.2
        elseif isapprox(a, -180.0; atol=1e-8)
            return 0.0
        else
            @warn "Angle encountered for which there is no CM table value (near +/-180 deg): $a"
            return nothing
        end
    end
end

# Helper: compute cmCoef and cm0 (port of Python __CMCoeff)
function _cm_coeff(alpha_deg::AbstractVector, cl_vec::AbstractVector, cm_vec::AbstractVector,
                    cl_high::Number, cd_high::Number, cm_high::Number)

    found_zero_lift = false
    cm0 = 0.0

    for i in 1:(length(cm_vec)-1)
        try
            if abs(alpha_deg[i]) < 20.0 && cl_vec[i] <= 0 && cl_vec[i+1] >= 0
                p = -cl_vec[i] / (cl_vec[i+1] - cl_vec[i])
                cm0 = cm_vec[i] + p * (cm_vec[i+1] - cm_vec[i])
                found_zero_lift = true
                break
            end
        catch
            # ignore indexing issues (mirrors Python 'try/except')
        end
    end

    if !found_zero_lift
        p = -cl_vec[1] / (cl_vec[2] - cl_vec[1])
        cm0 = cm_vec[1] + p * (cm_vec[2] - cm_vec[1])
    end

    alpha_high_rad = alpha_deg[end] * pi / 180.0
    XM = (-cm_high + cm0) / (cl_high * cos(alpha_high_rad) + cd_high * sin(alpha_high_rad))
    cmCoef = (XM - 0.25) / tan(alpha_high_rad - pi/2)

    return cmCoef, cm0
end


"""
Determine the minimum sets of overlapping parameters (meaning, it gets
rid of any parameter entries that are not available across the datapoints).
This reduces the parameter values to have full coverage over the dataset.
"""
function reduce_parameters(df)

    # Base case
    if size(df, 2) == 1

        return (unique(df[!, 1]), )

    # Recursively intersect the range of parameter values of each parameter
    # to find full coverage
    else

        outer_layer = unique(df[!, 1])
        reduced_values = (unique(df[!, 2]), )
        
        for val in outer_layer
            
            df_at_this_value = df[df[!, 1] .== val, 2:end]

            aux = reduce_parameters(df_at_this_value)
            
            reduced_values = (intersect(reduced_values[1], aux[1]), aux[2:end]...)
            
        end

        return (outer_layer, reduced_values...)
    end
    
end
##### END OF STRIPWISE ELEMENTS ################################################