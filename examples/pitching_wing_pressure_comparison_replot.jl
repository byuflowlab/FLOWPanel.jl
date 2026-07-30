# Replot an existing pitching-wing pressure comparison without rerunning Julia's
# panel-wake simulation.
#
# Usage:
#   julia --project examples/pitching_wing_pressure_comparison_replot.jl \
#       data/pitching_wing_pressure_comparison

import CSV
using DataFrames
import PythonPlot as plt

const METHODS = (
    (:bernoulli, "Unsteady Bernoulli"),
    (:edge_difference, "Laplace corrected edge"),
    (:corrected_hessian, "Laplace corrected Hessian"),
    (:surface_velocity, "Laplace surface velocity"),
)

function replot_pitching_wing_pressure_comparison(path)
    csv_path = joinpath(path, "pitching_wing_pressure_comparison.csv")
    isfile(csv_path) || throw(ArgumentError("Missing comparison CSV: $(csv_path)"))
    data = CSV.read(csv_path, DataFrame)
    cycles = data.t_over_T

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for (method, label) in METHODS
        ax.plot(cycles, data[!, Symbol("CL_$(method)")]; linewidth=1.25, label=label)
    end
    ax.set_xlabel("t/T")
    ax.set_ylabel("C_L")
    ax.grid(true, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    load_path = joinpath(path, "pitching_wing_pressure_loads.png")
    fig.savefig(load_path, dpi=170)
    plt.pyplot.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for (method, label) in METHODS[2:end]
        ax.plot(cycles, data[!, Symbol("pressure_l2_relative_$(method)")];
            linewidth=1.25, label=label)
    end
    ax.set_xlabel("t/T")
    ax.set_ylabel("relative pressure L2 error")
    ax.set_yscale("log")
    ax.grid(true, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    error_path = joinpath(path, "pitching_wing_pressure_errors.png")
    fig.savefig(error_path, dpi=170)
    plt.pyplot.close(fig)

    println("Wrote $(load_path)")
    println("Wrote $(error_path)")
    return (; loads=load_path, errors=error_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    path = length(ARGS) >= 1 ? ARGS[1] :
        get(ENV, "PRESSURE_COMPARISON_PATH", joinpath("data", "pitching_wing_pressure_comparison"))
    replot_pitching_wing_pressure_comparison(path)
end
