## Plot sweptwing refinement convergence: C_L vs refinement rows
# Two plots are produced:
#   1. chordwise refinement (C_L vs n_rfl) at fixed n_span = 24
#   2. spanwise  refinement (C_L vs n_span) at fixed n_rfl = 40
# Each has two source_half settings (positive / negative); both are plotted.

import PythonPlot as plt
using LaTeXStrings
using DelimitedFiles

include("plotformat.jl")

datadir = joinpath(@__DIR__, "..", "data", "sweptwing000")

function make_plot(csv_path, xcol, xlabel, title, out_path)
    data, header = readdlm(csv_path, ','; header=true)
    header = vec(header)

    col(name) = findfirst(==(name), header)
    x    = Float64.(data[:, col(xcol)])
    CL   = Float64.(data[:, col("CL")])
    half = strip.(string.(data[:, col("source_half")]))

    fig, ax = plt.subplots(figsize=(7, 4.5))

    for (key, label, color) in (("positive", "mesh: "*L"+y", "steelblue"),
                                ("negative", "mesh: "*L"-y", "darkorange"))
        mask = half .== key
        ax.plot(x[mask], CL[mask]; label=label,
                color=color, lw=1.5, marker="o", ms=5)
    end

    ax.set_xlabel(xlabel)
    ax.set_ylabel(L"$C_L$")
    ax.set_title(title)
    ax.legend(framealpha=0.8)
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.grid(true, lw=0.4, alpha=0.5)
    fig.tight_layout()

    fig.savefig(out_path, dpi=150)
    println("Saved $(out_path)")
end

# 1. chordwise-refinement convergence (n_span = 24)
make_plot(joinpath(datadir, "mirrored_convergence_n_span24.csv"),
          "n_rfl", "chordwise rows",
          "Swept wing chordwise-refinement convergence (n_span = 24)",
          joinpath(datadir, "mirrored_convergence_n_span24.png"))

# 2. spanwise-refinement convergence (n_rfl = 40)
make_plot(joinpath(datadir, "mirrored_convergence_n_rfl40.csv"),
          "n_span", "spanwise rows",
          "Swept wing spanwise-refinement convergence (n_rfl = 40)",
          joinpath(datadir, "mirrored_convergence_n_rfl40.png"))
