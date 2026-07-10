using CSV, DataFrames
using CairoMakie

# ------------------------------------------------------------

# Load CSV

# ------------------------------------------------------------

df = CSV.read("examples/wing_aileron/second.csv", DataFrame)

# ------------------------------------------------------------

# Filter 15152-panel rows for both solvers

# ------------------------------------------------------------

df_coupled = filter(r -> r.nps == 15152 && r.solver == "BackslashCoupled", df)
df_iter    = filter(r -> r.nps == 15152 && r.solver == "BackslashIterative", df)

# Ensure same ordering

sort!(df_coupled, :AOA)
sort!(df_iter, :AOA)

AOA = df_coupled.AOA

build_c = df_coupled.t_build
solve_c = df_coupled.t_solve

build_i = df_iter.t_build
solve_i = df_iter.t_solve

# ------------------------------------------------------------

# Plot

# ------------------------------------------------------------

fig = Figure()
ax = Axis(fig[1,1],
xlabel = "Angle of Attack (deg)",
ylabel = "Time (s)",
title  = "15152 Panels: Build + Solve Time"
)

x = collect(1:length(AOA))   # integer positions for Makie
w = 0.35

# --- Coupled (left, stacked) ---

b1 = barplot!(ax, x .- w/2, build_c;
width = w,
color = (:steelblue, 0.85)
)

b2 = barplot!(ax, x .- w/2, solve_c;
width = w,
offset = build_c,
color = (:lightskyblue, 0.85)
)

# --- Iterative (right, stacked) ---
b3 = barplot!(ax, x .+ w/2, build_i;
    width = w,
    color = (:orange, 0.85)
)

b4 = barplot!(ax, x .+ w/2, solve_i;
    width = w,
    offset = build_i,
    color = (:gold, 0.85)
)

# X ticks as AOA values
ax.xticks = (x, string.(AOA))

axislegend(ax,
    [b1, b2, b3, b4],
    ["Coupled build", "Coupled solve",
    "Iterative build", "Iterative solve"]
)

fig
