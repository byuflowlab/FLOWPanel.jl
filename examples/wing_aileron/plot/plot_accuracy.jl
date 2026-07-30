using CSV, DataFrames
using CairoMakie
using Statistics: mean

# ------------------------------------------------------------
# Load CSV
# ------------------------------------------------------------
df = CSV.read("examples/wing_aileron/results/krylov_solvers.csv", DataFrame)

# ------------------------------------------------------------
# Experimental (wind-tunnel) reference CL, indexed by AOA
# Full point set from examples/benchmarks.jl's `experimental` array
# ------------------------------------------------------------
experimental = [
    -9.58790170132325   -0.34195250659630627
    -7.841209829867704  -0.20211081794194818
    -5.988657844990541  -0.05435356200527597
    -4.136105860113403   0.09340369393140424
    -2.3894139886577754  0.2358839050132029
    -0.536862003780719   0.38100263852242966
     1.3156899810964084  0.5261213720316649
     3.168241965973529   0.6686015831134591
     5.02079395085066    0.8084432717678132
     6.926275992438619   0.9535620052770526
     8.672967863894144   1.0907651715039615
    10.472589792060482   1.2279683377308745
    12.325141776937613   1.3546174142480258
    14.230623818525515   1.4759894459102951
    15.077504725897917   1.528759894459108
    16.400756143667294   1.3308707124010595
]

# ------------------------------------------------------------
# Collapse repeated runs (same solver/nps/AOA) to a single mean CL,
# then attach the matching experimental CL and squared error
# ------------------------------------------------------------
data = combine(groupby(df, [:solver, :nps, :AOA]), :CL => mean => :CL)
data.CL_exp = experimental_CL.(data.AOA)
data.sq_err = (data.CL .- data.CL_exp) .^ 2

palette = [:steelblue, :orange, :seagreen, :firebrick, :purple, :goldenrod, :teal]

# ------------------------------------------------------------
# CL vs AOA, split into one plot per solver, one line per panel count
# ------------------------------------------------------------
function plot_CL_vs_AOA(data, solver_name, title)
    sub = filter(r -> r.solver == solver_name, data)
    nps_list = sort(unique(sub.nps))

    fig = Figure(size = (800, 500))
    ax = Axis(fig[1, 1],
        xlabel = "Angle of Attack (deg)",
        ylabel = "CL",
        title  = title,
    )

    for (i, n) in enumerate(nps_list)
        sn = sort(filter(r -> r.nps == n, sub), :AOA)
        color = palette[mod1(i, length(palette))]
        lines!(ax, sn.AOA, sn.CL; color = color, label = "nps=$n")
        scatter!(ax, sn.AOA, sn.CL; color = color)
    end

    exp_sorted = sortslices(experimental, dims = 1)
    lines!(ax, exp_sorted[:, 1], exp_sorted[:, 2];
        color = :black, linestyle = :dash, label = "Experimental")
    scatter!(ax, exp_sorted[:, 1], exp_sorted[:, 2]; color = :black, marker = :xcross)

    axislegend(ax, position = :lt)
    return fig
end

fig_coupled   = plot_CL_vs_AOA(data, "BackslashCoupled",   "BackslashCoupled: CL vs AOA")
fig_iterative = plot_CL_vs_AOA(data, "BackslashIterative", "BackslashIterative: CL vs AOA")

save("examples/wing_aileron/figures/accuracy_coupled.png", fig_coupled)
save("examples/wing_aileron/figures/accuracy_iterative.png", fig_iterative)

# ------------------------------------------------------------
# Average MSE (vs experimental CL) over the AOAs, per panel count
# ------------------------------------------------------------
mse_data = combine(groupby(data, [:solver, :nps]), :sq_err => mean => :mse)

fig_mse = Figure(size = (800, 500))
ax_mse = Axis(fig_mse[1, 1],
    xlabel = "Panel Count",
    ylabel = "Mean Squared Error (CL vs experimental)",
    title  = "Average MSE over AOAs, per Panel Count",
    xscale = log10,
)

for (i, s) in enumerate(sort(unique(mse_data.solver)))
    sm = sort(filter(r -> r.solver == s, mse_data), :nps)
    color = palette[mod1(i, length(palette))]
    lines!(ax_mse, sm.nps, sm.mse; color = color, label = s)
    scatter!(ax_mse, sm.nps, sm.mse; color = color)
end

axislegend(ax_mse, position = :lt)

save("examples/wing_aileron/figures/accuracy_mse.png", fig_mse)

fig_mse
