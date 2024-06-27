using Pkg
Pkg.activate("/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/")
using FLOWPanel
using FLOWPanel.StaticArrays
import FLOWPanel.FastMultipole as fmm
using LinearAlgebra
using BSON
using PythonPlot
include("auxiliary_functions.jl")
include("benchmark_wing.jl")

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
    ϕ_fmm = deepcopy(wing.potential)
    v_fmm = deepcopy(wing.velocity)

    # error calculation
    ε_ϕ = maximum(abs.(ϕ_direct - ϕ_fmm))
    ε_v = maximum(norm.(v_direct - v_fmm))

    return ε_ϕ, ε_v, t_fmm
end

function explore_fmm_parameters(wing, expansion_orders, leaf_sizes, multipole_thresholds)
    εs_ϕ = zeros(length(expansion_orders), length(leaf_sizes), length(multipole_thresholds))
    εs_v = deepcopy(εs_ϕ)
    ts_fmm = deepcopy(εs_ϕ)

    t_direct, ϕ_direct, v_direct = direct_results(wing)

    for (i_p,p) in enumerate(expansion_orders)
        for (i_l,l) in enumerate(leaf_sizes)
            for (i_m,m) in enumerate(multipole_thresholds)
                @show (p,l,m)
                ε_ϕ, ε_v, t_fmm = fmm_cost_error(wing, p, m, l, ϕ_direct, v_direct)
                εs_ϕ[i_p,i_l,i_m] = ε_ϕ
                εs_v[i_p,i_l,i_m] = ε_v
                ts_fmm[i_p,i_l,i_m] = t_fmm
            end
        end
    end

    return εs_ϕ, εs_v, t_direct, ts_fmm
end


function get_results(;
        expansion_orders = [1,2,3,4,5,6,7,8,9,10,12,14,16,18,20],
        leaf_sizes = [5,10,15,20,25,50,100,200,300,400],
        multipole_thresholds = [0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7],
        nc0 = 5, ns0 = 50,
        ms = 1:8,
        file="fmm_tuning.bson"
    )
    n_panels = [2*nc0*ns0*m^2 for m in ms]
    results = []

    for m in ms
        println("\nm = $m")
        println("-----")
        local wing = prepare_wing(;nc=nc0*m, ns=ns0*m)
        result = explore_fmm_parameters(wing, expansion_orders, leaf_sizes, multipole_thresholds)
        push!(results, result)
        BSON.@save file results n_panels expansion_orders leaf_sizes multipole_thresholds
    end

    return results
end

function clabel_fmt(x)
    str = "\$\\varepsilon = 10^{" * string(Int(round(x))) * "}\$"
    return LaTeXString(str)
end

function plot_results(results_bson::String)
    BSON.@load results_bson results n_panels expansion_orders leaf_sizes multipole_thresholds

    log10_leaf_sizes = log10.(leaf_sizes)
    cbar_labels = (L"\log_{10} \varepsilon_\phi", L"\log_{10} \varepsilon_v", "", L"\log_{10} t_{\text{fmm}}")

    # plot types
    for (i_type, type_name) in enumerate(("potential_error", "velocity_error", "t_direct", "t_fmm"))
        if i_type != 3
            for (i_n, n) in enumerate(n_panels)
                vmin = minimum(log10.(results[i_n][i_type]))
                vmax = maximum(log10.(results[i_n][i_type]))
                vmin_contour = minimum(log10.(results[i_n][2]))
                vmax_contour = maximum(log10.(results[i_n][2]))
                for (i_p,expansion_order) in enumerate(expansion_orders)
                    fig = figure("cost")
                    fig.clear()
                    # n == n_panels[1] && (fig.suptitle("expansion order: $expansion_order"))
                    fig.add_subplot(111, xlabel="leaf size", ylabel="multipole threshold")
                    ax = fig.get_axes()[0]
                    # expansion_order == expansion_orders[1] && ax.annotate(L"N"*" panels: $n", (-0.65,0.5), xycoords = "axes fraction", rotation=90, va="center", fontweight="bold", fontsize=10)
                    # ax.set_xscale("log")

                    data = transpose(log10.(results[i_n][i_type][i_p,:,:]))
                    cs = ax.contourf(log10_leaf_sizes, multipole_thresholds, data, levels=range(vmin,vmax,100))
                    ax.set_xticks(1.0:0.5:2.5)
                    ax.set_xticklabels([10,32,100,316])
                    fig.colorbar(cs, label=cbar_labels[i_type], ticks=range(round(vmin),round(vmax)))
                    fig.tight_layout()
                    fig.savefig("/Users/ryan/Dropbox/research/notebooks/img/20240618_fmm_tuning_wing/$(type_name)_n$(n)_p$(expansion_order).png")

                    if i_type == 4 # add error contours
                        data_contour = transpose(log10.(results[i_n][2][i_p,:,:]))
                        cs_levels = [-16,-14,-12,-10,-8,-6,-4,-3,-2,-1,0,1]
                        cs_contour = ax.contour(log10_leaf_sizes, multipole_thresholds, data_contour, colors="black", linewidths=1.0, linestyles="dashed", levels=cs_levels)
                        ax.clabel(cs_contour, levels=cs_levels, fmt=clabel_fmt)
                        fig.savefig("/Users/ryan/Dropbox/research/notebooks/img/20240618_fmm_tuning_wing/$(type_name)_error_contour__n$(n)_p$(expansion_order).png")
                    end

                end
            end
        end
    end

end

function optimize_results(results_bson::String; tolerances=[1e-3, 1e-6, 1e-9])
    BSON.@load results_bson results n_panels expansion_orders leaf_sizes multipole_thresholds
    return optimize_results(results, tolerances, expansion_orders, leaf_sizes, multipole_thresholds)
end

function optimize_results(results, tolerances, expansion_orders, leaf_sizes, multipole_thresholds)
    p_stars_ϕ = zeros(Int64, length(results), length(tolerances))
    theta_stars_ϕ = zeros(length(results), length(tolerances))
    n_stars_ϕ = zeros(Int64, length(results), length(tolerances))
    p_stars_v = zeros(Int64, length(results), length(tolerances))
    theta_stars_v = zeros(length(results), length(tolerances))
    n_stars_v = zeros(Int64, length(results), length(tolerances))

    for (i,result) in enumerate(results)
        for (i_tol, tol) in enumerate(tolerances)
            t_fmm_ϕ = Inf
            t_fmm_v = Inf
            for (i_p,p) in enumerate(expansion_orders)
                for (i_l,l) in enumerate(leaf_sizes)
                    for (i_m,m) in enumerate(multipole_thresholds)
                        # unpack results
                        err_ϕ = result[1][i_p,i_l,i_m]
                        err_v = result[2][i_p,i_l,i_m]
                        t_fmm = result[4][i_p,i_l,i_m]
                        if err_ϕ <= tol && t_fmm < t_fmm_ϕ # meets error tolerance in less time
                            p_stars_ϕ[i,i_tol] = p
                            n_stars_ϕ[i,i_tol] = l
                            theta_stars_ϕ[i,i_tol] = m
                            t_fmm_ϕ = t_fmm
                        end

                        if err_v <= tol && t_fmm < t_fmm_v # meets error tolerance in less time
                            p_stars_v[i,i_tol] = p
                            n_stars_v[i,i_tol] = l
                            theta_stars_v[i,i_tol] = m
                            t_fmm_v = t_fmm
                        end
                    end
                end
            end
        end
    end

    return p_stars_ϕ, theta_stars_ϕ, n_stars_ϕ, p_stars_v, theta_stars_v, n_stars_v
end

results_bson = "fmm_tuning.bson"
# results = get_results(; file=results_bson)

# plot_results(results_bson)

tolerances = [1e-3, 1e-6, 1e-9]
p_stars_ϕ, theta_stars_ϕ, n_stars_ϕ, p_stars_v, theta_stars_v, n_stars_v = optimize_results(results_bson; tolerances)

