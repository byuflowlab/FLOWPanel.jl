#=##############################################################################
# DESCRIPTION
    Postprocessing functions of swept wing example and experimental data
=###############################################################################

import PyPlot as plt
import CSV
import DataFrames: DataFrame
import OrderedCollections: OrderedDict
import PyPlot: @L_str
include(joinpath(pnl.examples_path, "plotformat.jl"))

expdata_path = joinpath(pnl.examples_path, "data")

dot(A, B) = sum(a*b for (a,b) in zip(A, B))
norm(A) = sqrt(mapreduce(x->x^2, +, A))
function cross!(out, A, B)
    out[1] = A[2]*B[3] - A[3]*B[2]
    out[2] = A[3]*B[1] - A[1]*B[3]
    out[3] = A[1]*B[2] - A[2]*B[1]
end
cross(A,B) = (out = zeros(3); cross!(out, A, B); return out)

function plot_Cps(body::Union{pnl.NonLiftingBody, pnl.AbstractLiftingBody}, controlpoints, spanposs, b;
                        spandirection=[0, 1, 0], AOA=nothing,
                        _fig=nothing, _axs=nothing,
                        ttl=nothing, xscaling=1,
                        stl=".-", xlims=[-0.1, 1.1], ylims=[1.0, -0.7],
                        plot_optargs=(label="FLOWPanel",),
                        plot_exp=true, lbl_exp="Experimental",
                        plot_vsp=true, lbl_vsp="VSPAERO",
                        out=[])

    npos = length(spanposs)

    # Plot slices of the wing along the span
    fig = _fig==nothing ? plt.figure(figsize=[7, 5*0.75]*2/3 .* [2, ceil(npos/2)]) : _fig
    axs = _axs==nothing ? fig.subplots(ceil(Int, npos/2), 2) : _axs
    axs = _axs==nothing ? [axs[i, j] for j in 1:size(axs, 2), i in 1:size(axs, 1)] : axs

    for (axi, (ax, spanpos)) in enumerate(zip(axs, spanposs))

        ax.set_title(ttl==nothing ? L"$2y/b = $"*"$(spanpos)" : ttl)

        # Plot experimental
        if plot_exp && AOA != nothing

            data_exp = weber_Cps["$AOA"]
            keyi = findfirst(str -> parse(Float64, str)==abs(spanpos), collect(keys(data_exp)))

            if keyi==nothing
                println("Experimental data at spanpos=$(abs(spanpos)) not found.")
            else
                key = collect(keys(data_exp))[keyi]
                ax.plot(weber_xoc_up, data_exp[key][1], "ok", markersize=5, label=lbl_exp)
                ax.plot(weber_xoc_lo, data_exp[key][2], "ok", markersize=5)
            end
        end


            println(spanpos)

        # Plot VSPAERO
        if plot_vsp && AOA == 4.2

            rowstart =  spanpos==0.041 ? 1139 :
                        spanpos==0.163 ? 1169 :
                        spanpos==0.245 ? 1184 :
                        spanpos==0.510 ? 1199 :
                        nothing

            if !isnothing(rowstart)

                data_vsp_cp = CSV.read(vsp_file, DataFrame; skipto=rowstart+7, limit=1)
                data_vsp_loc = CSV.read(vsp_file, DataFrame; skipto=rowstart+12, limit=3)

                cp_vsp = [val for val in data_vsp_cp[1, 2:end]]
                x_vsp = [val for val in data_vsp_loc[1, 2:end]]
                x_vsp .-= minimum(x_vsp)
                x_vsp ./= maximum(x_vsp) .- minimum(x_vsp)

                ax.plot(x_vsp, cp_vsp, "-",
                        color=color_vsp, alpha=alpha_vsp, label=lbl_vsp, linewidth=1, )

            end
        end

        # Position of grid columns that are the closest to target
        points, Cps = pnl.slicefield(body, controlpoints, "Cp", spanpos*b/2, spandirection, false)

        chordposs = points[1, :]
        chordposs .-= minimum(chordposs)

        # Plot FLOWPanel
        ax.plot(chordposs*xscaling, Cps, stl; clip_on=false, plot_optargs...)

        if xlims!=nothing; ax.set_xlim(xlims); end;
        if ylims!=nothing; ax.set_ylim(ylims); end;

        if xlims!=nothing; ax.set_xticks(0:0.2:1); end;
        if ylims!=nothing; ax.set_yticks(ylims[1]:-0.5:ylims[2]); end

        if axi >= length(axs)-1
            ax.set_xlabel(L"x/c")
        end
        if axi%2==1
            ax.set_ylabel(L"Pressure coefficient $C_p$")
        end

        if axi==1
            ax.legend(loc="best", fontsize=10, frameon=false, reverse=true)
        end

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        push!(out, [spanpos, chordposs*xscaling, Cps])
    end

    return fig, axs

end

function plot_Cps(body, spanposs, b; optargs...)
    normals_b = pnl._calc_normals(body)
    controlpoints_b = pnl._calc_controlpoints(body, normals_b)

    return plot_Cps(body, controlpoints_b, spanposs, b; optargs...)
end

function plot_deltaCps(body::Union{pnl.NonLiftingBody, pnl.AbstractLiftingBody}, controlpoints, spanposs, b;
                        spandirection=[0, 1, 0], AOA=nothing,
                        _fig=nothing, _axs=nothing,
                        ttl=nothing, xscaling=1,
                        stl=".-", xlims=[-0.1, 1.1], ylims=[1.0, -0.7],
                        plot_optargs=(label="FLOWPanel",),
                        plot_exp=true, lbl_exp="Experimental",
                        out=[])

    npos = length(spanposs)

    # Plot slices of the wing along the span
    fig = _fig==nothing ? plt.figure(figsize=[7, 5*0.75]*2/3 .* [2, ceil(npos/2)]) : _fig
    axs = _axs==nothing ? fig.subplots(ceil(Int, npos/2), 2) : _axs
    axs = _axs==nothing ? [axs[i, j] for j in 1:size(axs, 2), i in 1:size(axs, 1)] : axs

    for (axi, (ax, spanpos)) in enumerate(zip(axs, spanposs))

        ax.set_title(ttl==nothing ? L"$2y/b = $"*"$(spanpos)" : ttl)

        # Plot experimental
        if plot_exp && AOA==4.2

            spanpos_str = calc_id(abs(spanpos), 3, 2)

            try
                data = CSV.read(joinpath(expdata_path, "weber1958-fig2-y2b-$(spanpos_str).csv"), DataFrame)
                ax.plot(data[!, 1], data[!, 2], "ok", markersize=5, label=lbl_exp)
                ax.plot(data[!, 1], data[!, 2], "ok", markersize=5)
            catch e
                println(e)
                println("Experimental data at spanpos=$(abs(spanpos)) not found.")
            end
        end

        # Position of grid columns that are the closest to target
        points, Cps = pnl.slicefield(body, controlpoints, "Cp", spanpos*b/2, spandirection, false)

        chordposs = points[1, :]
        chordposs .-= minimum(chordposs)

        # Calculate difference between upper and lower surface
        np = Int(length(Cps)/2)
        diffchordposs = [(chordposs[i] + chordposs[end-(i-1)])/2 for i in 1:np]
        diffCps = [Cps[end-(i-1)] - Cps[i] for i in 1:np]

        # Plot FLOWPanel
        ax.plot(diffchordposs*xscaling, diffCps, stl; plot_optargs...)

        if xlims!=nothing; ax.set_xlim(xlims); end;
        if ylims!=nothing; ax.set_ylim(ylims); end;

        if xlims!=nothing; ax.set_xticks(0:0.2:1); end;
        if ylims!=nothing; ax.set_yticks(ylims[1]:-0.5:ylims[2]); end

        if axi >= length(axs)-1
            ax.set_xlabel(L"x/c")
        end
        if axi%2==1
            ax.set_ylabel(L"$\Delta C_p$")
        end

        if axi==1
            ax.legend(loc="best", fontsize=10, frameon=false, reverse=true)
        end

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        push!(out, [spanpos, chordposs*xscaling, Cps])
    end

    return fig, axs

end

function plot_deltaCps(body, spanposs, b; optargs...)
    normals_b = pnl._calc_normals(body)
    controlpoints_b = pnl._calc_controlpoints(body, normals_b)

    return plot_deltaCps(body, controlpoints_b, spanposs, b; optargs...)
end



function plot_loading(body::Union{pnl.NonLiftingBody, pnl.AbstractLiftingBody}, Lhat, Dhat, b;
                        spandirection=[0, 1, 0], dimspan=2, dimchord=1,
                        to_plot=[1, 2, 3], yscalings=[1.0, 1.0, 1.0],
                        _fig=nothing, _axs=nothing,
                        stl=".-", xlims=[-1, 1], ylims=nothing,
                        plot_optargs=(label="FLOWPanel",),
                        plot_exp=true, lbl_exp="Experimental", AOA=nothing,
                        plot_vsp=true, lbl_vsp="VSPAERO",
                        out=[])

    # Calculate sectional force along the span
    fs, spanpos = pnl.calcfield_sectionalforce(body; spandirection=spandirection,
                                                    dimspan=dimspan, dimchord=dimchord)

    # Decompose the force into lift, drag, and sideslip
    lds = pnl.decompose(fs, Lhat, Dhat)

    # Plot loading
    fig = _fig==nothing ? plt.figure(figsize=[7, 5*0.75]*2/3 .* [length(to_plot), 1]) : _fig
    axs = _axs==nothing ? fig.subplots(1, length(to_plot)) : _axs
    axs = _axs==nothing ? length(to_plot)==1 ? [axs] : [axs[i, j] for j in 1:size(axs, 2), i in 1:size(axs, 1)] : axs

    for (axi, (ax, pi)) in enumerate(zip(axs, to_plot))

        # Plot experimental
        if plot_exp && pi != 3

            vals_exp =  [[cls_web, cds_web], [cns_web, cts_web]][1][pi]

            rowi = findfirst(a -> a==AOA, alphas_web)
            if rowi != nothing

                for f in [-1, 1]
                    ax.plot(f*y2b_web, vals_exp[rowi, :], "o--k", label=lbl_exp^(f==1),
                                                linewidth=0.5, alpha=1.0, clip_on=true)
                end
            else
                println("Experimental data at AOA=$(AOA) not found; valid AOAs are $(alphas_web).")
            end

        end

        # Plot VSPAERO
        if plot_vsp

            rowstart =  AOA==2.1 ? 482 :
                        AOA==4.2 ? 544 :
                        AOA==6.3 ? 606 :
                        AOA==8.4 ? 668 :
                        nothing

            if !isnothing(rowstart)

                data_vsp = CSV.read(vsp_file, DataFrame; skipto=rowstart+20, limit=24)

                ypos_vsp = [val for val in data_vsp[1, 2:end]]
                cl_vsp = [val for val in data_vsp[9, 2:end]]
                cd_vsp = [val for val in data_vsp[7, 2:end]]

                ypos_vsp .-= 0.5
                ypos_vsp .*= 2.0

                ax.plot(ypos_vsp, pi==1 ? cl_vsp : cd_vsp, "-",
                        color=color_vsp, alpha=alpha_vsp*0.75, label=lbl_vsp, linewidth=1)

            end
        end

        # Plot FLOWPanel
        ax.plot(spanpos*2/b, lds[pi, :]*yscalings[pi], stl; plot_optargs...)

        if xlims!=nothing; ax.set_xlim(xlims); end;
        if ylims!=nothing; ax.set_ylim(ylims[axi][1:2]); end;

        if xlims!=nothing; ax.set_xticks(xlims[1]:0.25:xlims[2]); end;
        if ylims!=nothing; ax.set_yticks(ylims[axi][1]:ylims[axi][3]:ylims[axi][2]); end

        ax.set_xlabel(L"Span position $2y/b$")
        ax.set_ylabel("Sectional "*[L"lift $\ell$", L"drag $d$", L"sideslip $s$"][pi]*" (N/m)")

        if axi==1
            ax.legend(loc="lower center", fontsize=10, frameon=false, reverse=true)
        end

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        push!(out, [spanpos*2/b, lds[pi, :]*yscalings[pi]])
    end

    return fig, axs

end

function plot_loading(multibody::pnl.MultiBody, args...; _fig=nothing, _axs=nothing,
                        plots_optargs=[(label="FLOWPanel",) for i in 1:10],
                        plot_exp=true, plot_vsp=true,
                        optargs...)

    fig, axs = _fig, _axs

    for (bi, body) in enumerate(multibody.bodies)

        fig, axs = plot_loading(body, args...; optargs..., _fig=fig, _axs=axs,
                                    plot_optargs=plots_optargs[bi],
                                    plot_exp=plot_exp*(bi==1), plot_vsp=plot_vsp*(bi==1))
    end

    return fig, axs

end


function calc_id(val, nd, f)
    num = ceil(Int, val*10^f)
    num = "0"^(num==0 ? nd-1 : nd - 1 - floor(Int, log10(num)))*"$num"

    return num
end


# Pressure distributions in Table 2
weber_xoc_up = [0, 0.01, 0.03, 0.08, 0.15, 0.225, 0.35, 0.50 ,0.65, 0.75, 0.85, 0.95]
weber_xoc_lo = [0.01, 0.03, 0.08, 0.15, 0.225, 0.35, 0.50 ,0.65, 0.75, 0.85, 0.95]
weber_Cps = OrderedDict(  # weber_Cps[AOA][2y/b][upper/lower][xoc]
        "4.2" => OrderedDict(
            "0" => [
                        # Upper surface
                        [0.980, 0.190, -0.035, -0.185, -0.275, nothing, -0.380,
                        -0.325, -0.240, -0.170, nothing, -0.035],
                        # Lower surface
                        [0.745, 0.460, 0.210, 0.090, 0.005, -0.095, -0.095,
                        -0.060, -0.035, -0.005, 0.025]
                       ],
            "0.041" => [
                        # Upper surface
                        [0.425, -0.235, -0.370, -0.375, -0.400, -0.415, -0.385,
                        -0.280, -0.180, -0.120, -0.040, 0.050],
                        # Lower surface
                        [0.480, 0.315, 0.130, -0.005, -0.055, -0.115,
                        -0.095, -0.050, nothing, 0.030, 0.060]
                       ],
            "0.082" => [
                        # Upper surface
                        [0.270, -0.355, -0.485, -0.465, -0.465, -0.455,
                        -0.400, -0.275, -0.155, -0.090, -0.020, 0.060],
                        # Lower surface
                        [0.480, 0.310, 0.085, -0.025, -0.080, -0.120,
                        -0.085, -0.020, 0, 0.040, 0.075]
                       ],
            "0.163" => [
                        # Upper surface
                        [0.165, -0.485, -0.550, -0.540, -0.520, -0.465,
                        -0.400, -0.235, nothing, -0.060, -0.005, 0.060],
                        # Lower surface
                        [0.455, 0.280, 0.085, -0.030, -0.085, -0.125,
                        -0.070, -0.010, 0.020, 0.040, 0.065]
                       ],
            "0.245" => [
                        # Upper surface
                        [0.220, -0.505, -0.615, -0.535, -0.505, -0.465,
                        -0.370, -0.205, -0.100, -0.050, 0.015, 0.070],
                        # Lower surface
                        [0.450, 0.300, 0.085, -0.025, -0.075, -0.110,
                        -0.065, 0, 0.035, 0.055, 0.085]
                       ],
            "0.367" => [
                        # Upper surface
                        [0.270, -0.565, -0.595, -0.545, -0.515, -0.460,
                        -0.345, -0.185, -0.100, -0.030, 0.015, 0.080],
                        # Lower surface
                        [0.465, 0.295, 0.115, -0.015, -0.065, -0.105,
                        -0.035, 0.005, 0.025, 0.060, 0.080]
                       ],
            "0.510" => [
                        # Upper surface
                        [0.380, -0.515, -0.620, -0.575, -0.505, -0.445, -0.345,
                        -0.190, -0.095, -0.030, 0.025, 0.085],
                        # Lower surface
                        [0.480, 0.355, 0.130, 0, -0.060, -0.095, -0.035,
                        0, 0.030, 0.055, 0.075]
                       ],
            "0.653" => [
                        # Upper surface
                        [0.170, -0.575, -0.640, -0.595, -0.510, -0.445,
                        -0.355, -0.205, -0.090, -0.015, 0.025, 0.085],
                        # Lower surface
                        [0.475, 0.335, 0.125, 0.010, -0.065, -0.110,
                        -0.040, 0.010, 0.025, 0.055, 0.085]
                       ]
        )
)

# Weber's experimental loading distribution from Table 3
alphas_web = [2.1, 4.2, 6.3, 8.4, 10.5]
y2b_web = [0.0, 0.041, 0.082, 0.163, 0.245, 0.367, 0.510, 0.653, 0.898, 0.949]
cls_web = [
            0.118 0.121 0.126 0.129 0.129 0.129 0.131 0.125 NaN 0.087;      # AOA = 2.1deg
            0.235 0.241 0.248 0.253 0.251 0.251 0.251 0.246 0.192 0.171;    # AOA = 4.2deg
            0.351 0.358 0.367 0.374 0.375 0.373 0.377 0.365 NaN 0.256;      # AOA = 6.3deg
            0.466 0.476 0.483 0.494 0.494 0.493 0.493 0.48 NaN 0.34;        # AOA = 8.4deg
            0.577 0.589 0.597 0.607 0.611 0.605 0.599 0.587 0.415 0.401     # AOA = 10.5deg
         ]
cds_web = [
            0.044 0.014 0.007 0.002 0.0 -0.001 0.0 -0.001 -0.009 -0.01;     # AOA = 0deg
            0.047 0.016 0.01 0.004 0.002 0.001 0.002 0.001 NaN -0.009;
            0.059 0.025 0.016 0.009 0.007 0.006 0.006 0.004 -0.002 -0.007;
            0.078 0.039 0.028 0.017 0.015 0.012 0.012 0.009 NaN -0.005;
            0.104 0.058 0.045 0.029 0.026 0.023 0.022 0.016 NaN -0.001;
            0.138 0.084 0.065 0.044 0.041 0.037 0.035 0.026 0.009 0.004
         ]

cns_web = [
            0.12 0.121 0.127 0.129 0.129 0.129 0.131 0.125 NaN 0.087;
            0.239 0.243 0.249 0.254 0.251 0.25 0.251 0.246 0.192 0.17;
            0.357 0.36 0.368 0.374 0.374 0.372 0.376 0.364 NaN 0.254;
            0.476 0.479 0.484 0.493 0.492 0.491 0.491 0.477 NaN 0.336;
            0.592 0.595 0.599 0.605 0.608 0.601 0.595 0.582 0.41 0.395
         ]

cts_web = [
            0.044 0.014 0.007 0.002 0.0 -0.001 0.0 -0.001 -0.009 -0.01;     # AOA = 0deg
            0.043 0.012 0.005 -0.001 -0.003 -0.004 -0.003 -0.004 NaN -0.012;
            0.042 0.007 -0.002 -0.01 -0.012 -0.012 -0.013 -0.014 -0.016 -0.02;
            0.039 -0.001 -0.013 -0.024 -0.027 -0.029 -0.03 -0.031 NaN -0.033;
            0.035 -0.012 -0.027 -0.044 -0.047 -0.05 -0.05 -0.053 NaN -0.05;
            0.03 -0.025 -0.045 -0.067 -0.071 -0.074 -0.075 -0.082 -0.067 -0.069
         ]

# Add the zero-AOA drag that they substracted in the experiment
cds_web = mapslices(x-> x .+ cds_web[1, :], cds_web[2:end, :]; dims=2)
cts_web = mapslices(x-> x .+ cts_web[1, :], cts_web[2:end, :]; dims=2)

# Integrated coefficients from Table 4B
CLs_web = [0.121, 0.238, 0.350, 0.456, 0.559]
CDs_web = [nothing, 0.005, 0.012, 0.022, 0.035]

# File with VSP solution
vsp_file = joinpath(pnl.examples_path, "data", "weber_vspaero.csv")
color_vsp = "tab:orange"
alpha_vsp = 0.3
