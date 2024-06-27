using Pkg
Pkg.activate("/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/")
using FLOWPanel
using FLOWPanel.StaticArrays
import FLOWPanel.FastMultipole as fmm
using LinearAlgebra
using BSON
using PythonPlot
include("auxiliary_functions.jl")
include("benchmark_wing_1b.jl")

function direct_results(wing)
    # reset wing
    reset_potential_velocity!(wing)

    # direct calculation
    t_direct = @elapsed fmm.direct!(wing; velocity_gradient=false)
    ϕ_direct = deepcopy(wing.potential)
    v_direct = deepcopy(wing.velocity)

    return t_direct, ϕ_direct, v_direct
end

function fmm_cost_error(wing, expansion_order, multipole_threshold, leaf_size, ϕ_direct, v_direct)

    # reset wing
    reset_potential_velocity!(wing)

    # fmm calculation
    t_fmm = @elapsed fmm.fmm!(wing; velocity_gradient=false, expansion_order, multipole_threshold, leaf_size)
    ϕ_fmm = wing.potential[:,:,1]
    v_fmm = wing.velocity[:,:,1]

    # error calculation
    ε_ϕ = maximum(abs.(ϕ_direct - ϕ_fmm))
    ε_v_rel = maximum(norm.(v_direct - v_fmm) ./ norm.(v_direct))
    ε_v_abs = maximum(norm.(v_direct - v_fmm))

    return ε_ϕ, ε_v_rel, ε_v_abs, t_fmm, ϕ_fmm, v_fmm
end

function explore_fmm_parameters(wing, expansion_orders, leaf_sizes, multipole_thresholds)
    εs_ϕ = zeros(length(expansion_orders), length(leaf_sizes), length(multipole_thresholds))
    εs_v_abs = deepcopy(εs_ϕ)
    εs_v_rel = deepcopy(εs_ϕ)
    ts_fmm = deepcopy(εs_ϕ)
    ϕs_fmm = Array{Matrix{Float64},3}(undef,length(expansion_orders), length(leaf_sizes), length(multipole_thresholds))
    vs_fmm = Array{Matrix{SVector{3,Float64}},3}(undef,length(expansion_orders), length(leaf_sizes), length(multipole_thresholds))

    t_direct, ϕ_direct, v_direct = direct_results(wing)

    for (i_p,p) in enumerate(expansion_orders)
        for (i_l,l) in enumerate(leaf_sizes)
            for (i_m,m) in enumerate(multipole_thresholds)
                @show (length(wing.panels),p,l,m)
                ε_ϕ, ε_v_rel, ε_v_abs, t_fmm, ϕ_fmm, v_fmm = fmm_cost_error(wing, p, m, l, ϕ_direct, v_direct)
                εs_ϕ[i_p,i_l,i_m] = ε_ϕ
                εs_v_rel[i_p,i_l,i_m] = ε_v_rel
                εs_v_abs[i_p,i_l,i_m] = ε_v_abs
                ts_fmm[i_p,i_l,i_m] = t_fmm
                @show typeof(ϕ_fmm) typeof(v_fmm)
                ϕs_fmm[i_p,i_l,i_m] = ϕ_fmm
                vs_fmm[i_p,i_l,i_m] = v_fmm
            end
        end
    end

    return εs_ϕ, εs_v_abs, εs_v_rel, t_direct, ts_fmm, ϕ_direct, v_direct, ϕs_fmm, vs_fmm
end


function get_results(;
        expansion_orders = [1,2,3,4,5,6,7,8,9,10,12,14,16,18,20],
        leaf_sizes = [10,15,20,25,50,100,200,300,400],
        multipole_thresholds = [0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75],
        nc0 = 5, ns0 = 50,
        ms = 1:8,
        file="fmm_tuning_2"
    )

    for m in ms
        println("\nm = $m")
        println("-----")
        local wing = prepare_wing(;nc=nc0*m, ns=ns0*m)
        result = explore_fmm_parameters(wing, expansion_orders, leaf_sizes, multipole_thresholds)
        n = 2*nc0*ns0*m^2
        BSON.@save file*"_m$m.bson" result n expansion_orders leaf_sizes multipole_thresholds
    end
end

function clabel_fmt(x)
    str = "\$\\varepsilon = 10^{" * string(Int(round(x))) * "}\$"
    return LaTeXString(str)
end

function colorbar_format(x,y)
    str = string(round(10^x, sigdigits=1))
    return str
end

function plot_results(results_bson_prefix::String, ms; contour_type="abs", nc0=5, ns0=50)

    # cbar_labels = (L"\log_{10} \varepsilon_\phi", L"\log_{10} \varepsilon_v", "", L"\log_{10} t_{\text{fmm}}")

    # cbar_label = contour_type == "abs" ? L"\left| \log_{10} \left| \vec{v}_{fmm} - \vec{v}_{direct} \right| \right|_\infty" : L"\left| \log_{10} \left| \frac{\vec{v}_{fmm} - \vec{v}_{direct}}{|\vec{v}|} \right| \right|_\infty"
    cbar_label = "fmm cost, seconds"

    for m in ms

        # load results
        BSON.@load results_bson_prefix*"_m$m.bson" result n expansion_orders leaf_sizes multipole_thresholds
        n = 2*nc0*ns0*m^2
        εs_ϕ, εs_v_abs, εs_v_rel, t_direct, ts_fmm, ϕ_direct, v_direct, ϕs_fmm, vs_fmm = result

        # log scale
        log10_leaf_sizes = log10.(leaf_sizes)
        # log10_εs_v_abs = log10.(εs_v_abs)
        # log10_εs_v_rel = log10.(εs_v_rel)
        log10_εs = contour_type == "abs" ? log10.(εs_v_abs) : log10.(εs_v_rel)
        log10_ts_fmm = log10.(ts_fmm)

        # for labels
        vmin = minimum(log10_ts_fmm)
        vmax = maximum(log10_ts_fmm)
        vmin_contour = minimum(log10_εs)
        vmax_contour = maximum(log10_εs)

        # cs_levels
        if m==1
            cs_levels = [-14, -6, -2]
        elseif m==2
            cs_levels = [-14, -6, -2]
        elseif m==3
            cs_levels = [-14, -10, -6, -2]
        else
            cs_levels = collect(Int(floor(vmin_contour)):Int(ceil(vmax_contour)))
        end
        cs_levels = [-9, -6, -3]

        # new plot for each expansion order
        for (i_p,expansion_order) in enumerate(expansion_orders)

            fig = figure("cost")
            fig.clear()
            fig.add_subplot(111, xlabel="leaf size", ylabel="multipole threshold")
            ax = fig.get_axes()[0]
            # ax.set_xscale("log")

            data = transpose(view(log10_ts_fmm,i_p,:,:))
            cs = ax.contourf(log10_leaf_sizes, multipole_thresholds, data, levels=range(vmin,vmax,100))

            #TODO: xticks
            ax.set_xticks(1.0:0.5:2.5)
            ax.set_xticklabels([10,32,100,316])
            fig.colorbar(cs, label=cbar_label, ticks=range(round(vmin), stop=round(vmax), step=0.5), format=colorbar_format)
            fig.tight_layout()
            fig.savefig("/Users/ryan/Dropbox/research/notebooks/img/20240621_fmm_tuning_wing/fmm_cost_n$(n)_p$(expansion_order).png")

            # add error contours
            data_contour = transpose(view(log10_εs,i_p,:,:))
            cs_contour = ax.contour(log10_leaf_sizes, multipole_thresholds, data_contour, colors="black", linewidths=1.0, linestyles="dashed", levels=cs_levels)
            ax.clabel(cs_contour, levels=cs_levels, fmt=clabel_fmt)
            fig.savefig("/Users/ryan/Dropbox/research/notebooks/img/20240621_fmm_tuning_wing/fmm_cost_error_contour_" * contour_type * "_n$(n)_p$(expansion_order).png")
        end # p
    end # m
end

function optimize_results(results_bson_prefix::String, ms; tolerances=[1e-3, 1e-6, 1e-9], ns0=50, nc0=5)
    optimal_p = []
    optimal_m = []
    optimal_l = []
    for m in ms
        BSON.@load results_bson_prefix*"_m$m.bson" result n expansion_orders leaf_sizes multipole_thresholds
        p_stars_v, theta_stars_v, n_stars_v = optimize_results(result, tolerances, expansion_orders, leaf_sizes, multipole_thresholds)
        push!(optimal_p, p_stars_v)
        push!(optimal_m, theta_stars_v)
        push!(optimal_l, n_stars_v)
    end
    ns = ns0*nc0*2*ms .^2
    return optimal_p, optimal_m, optimal_l, ns
end

function optimize_results(result, tolerances, expansion_orders, leaf_sizes, multipole_thresholds)
    p_stars_v = zeros(Int64, length(tolerances))
    theta_stars_v = zeros(length(tolerances))
    n_stars_v = zeros(Int64, length(tolerances))

    for (i_tol, tol) in enumerate(tolerances)
        t_fmm_v = Inf
        for (i_p,p) in enumerate(expansion_orders)
            for (i_l,l) in enumerate(leaf_sizes)
                for (i_m,m) in enumerate(multipole_thresholds)
                    # unpack results
                    εs_ϕ, εs_v_abs, εs_v_rel, t_direct, ts_fmm, ϕ_direct, v_direct, ϕs_fmm, vs_fmm = result

                    err_v = εs_v_abs[i_p,i_l,i_m]
                    t_fmm = ts_fmm[i_p,i_l,i_m]
                    if err_v <= tol && t_fmm < t_fmm_v # meets error tolerance in less time
                        p_stars_v[i_tol] = p
                        n_stars_v[i_tol] = l
                        theta_stars_v[i_tol] = m
                        t_fmm_v = t_fmm
                    end
                end
            end
        end
    end


    return p_stars_v, theta_stars_v, n_stars_v
end

results_bson = "fmm_tuning_2"
ms = 1:8
# results = get_results(; ms=1:8, file=results_bson)

# plot_results(results_bson, ms; contour_type="abs")

optimal_p, optimal_m, optimal_l, ns = optimize_results(results_bson, ms; tolerances=[1e-3, 1e-6, 1e-9], ns0=50, nc0=5)

