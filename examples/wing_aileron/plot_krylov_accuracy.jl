using CSV, DataFrames
using CairoMakie
using Statistics: mean

# ------------------------------------------------------------
# Load CSV (KrylovCoupled / KrylovCoupled-FMM vs experimental CL)
# ------------------------------------------------------------
df = CSV.read("examples/wing_aileron/krylov_solvers.csv", DataFrame)
df = filter(r -> r.CL != "ERROR", df)
df.CL = Float64.(df.CL)
df.AOA = Float64.(df.AOA)

# --- Experimental data ---
CL_exp = [
 -9.58790170132325  -0.34195250659630627
 -7.841209829867704 -0.20211081794194818
 -5.988657844990541 -0.05435356200527597
 -4.136105860113403  0.09340369393140424
 -2.3894139886577754 0.2358839050132029
 -0.536862003780719  0.38100263852242966
  1.3156899810964084 0.5261213720316649
  3.168241965973529  0.6686015831134591
  5.02079395085066   0.8084432717678132
  6.926275992438619  0.9535620052770526
  8.672967863894144  1.0907651715039615
 10.472589792060482  1.2279683377308745
 12.325141776937613  1.3546174142480258
 14.230623818525515  1.4759894459102951
 15.077504725897917  1.528759894459108
 16.400756143667294  1.3308707124010595
 18.465028355387517  1.2807387862796877
]

aoa_exp = CL_exp[:,1]
cl_exp  = CL_exp[:,2]

function experimental_CL(aoa; tol=1e-3)
    idx = findfirst(a -> isapprox(a, aoa; atol=tol), aoa_exp)
    idx === nothing && error("No experimental CL found for AOA=$aoa")
    return cl_exp[idx]
end

df.CL_exp = experimental_CL.(df.AOA)
df.sq_err = (df.CL .- df.CL_exp) .^ 2

palette = [:steelblue, :orange, :seagreen, :firebrick, :purple, :goldenrod, :teal]

# ------------------------------------------------------------
# CL vs AOA for a single solver, one color per panel count,
# each point labeled with its panel count (nps)
# ------------------------------------------------------------
function plot_solver(df, solver_name)
    sub = filter(r -> r.solver == solver_name, df)
    nps_list = sort(filter(n -> n >= 25000, unique(sub.nps)))

    fig = Figure(size = (800, 550))
    ax = Axis(fig[1, 1],
        xlabel = "Angle of Attack (deg)",
        ylabel = "CL",
        # title  = solver_name,
    )

    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.topspinevisible = false
    ax.rightspinevisible = false

    for (i, n) in enumerate(nps_list)
        sn = sort(filter(r -> r.nps == n, sub), :AOA)
        color = palette[mod1(i, length(palette))]
        lines!(ax, sn.AOA, sn.CL; color = color, label = "nps=$n")
        scatter!(ax, sn.AOA, sn.CL; color = color)
        # stagger label offsets per panel-count group so labels for nearby
        # points (different nps, similar AOA/CL) don't stack on top of
        # each other
        # voffset = 6 + 12 * (i - 1)
        # for row in eachrow(sn)
        #     text!(ax, row.AOA, row.CL; text = string(row.nps),
        #         fontsize = 10, color = color,
        #         offset = (6, voffset))
        # end
    end

    lines!(ax, aoa_exp, cl_exp;
        color = :black, linestyle = :dash, label = "Experimental")
    scatter!(ax, aoa_exp, cl_exp; color = :black, marker = :xcross)
    
    axislegend(ax, position = :lt)
    return fig
end

fig_coupled = plot_solver(df, "KrylovCoupled")
fig_fmm     = plot_solver(df, "KrylovCoupled-FMM")

save("examples/wing_aileron/krylov_accuracy_coupled_high_nps.png", fig_coupled)
save("examples/wing_aileron/krylov_accuracy_fmm_high_nps.png", fig_fmm)

# ------------------------------------------------------------
# MSE (vs experimental CL) per solver, per panel count
# ------------------------------------------------------------
mse_data = combine(groupby(df, [:solver, :nps]), :sq_err => mean => :mse)
sort!(mse_data, [:solver, :nps])

println(mse_data)

fig_coupled
