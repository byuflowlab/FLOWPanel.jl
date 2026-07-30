#=##############################################################################
# DESCRIPTION
#   Regenerate the pitching-wing plots from a completed run's CSV output,
#   without rerunning the simulation. Reads `pitching_wing_unsteady.csv` (and,
#   if present, `pitching_wing_static_polar.csv`) and re-emits the convergence,
#   α-hysteresis, static-polar, and experiment-comparison figures with the
#   current plotting style. The comparison figures overlay, per spanwise η, the
#   unsteady loop on the dynamic experiment and the static polar on the static
#   experiment (data/pitching_wing_exp).
#
# USAGE
#   julia --project examples/pitching_wing_replot.jl [run_dir] [n_loops]
#
#   run_dir  directory holding the CSVs (default: data/pitching_wing)
#   n_loops  cycles to overlay in the hysteresis loop (default: 2)
#
#   Set MPLBACKEND=Agg to run headless (no window).
=###############################################################################

import FLOWPanel as pnl
include(joinpath(@__DIR__, "pitching_wing.jl"))
include(joinpath(@__DIR__, "..", "data", "pitching_wing_exp", "load.jl"))
import DelimitedFiles: readdlm

# Parse the section η values from the `cl_eta_0p2500`-style CSV header columns.
function _replot_section_eta(header)
    etas = Float64[]
    for name in header
        startswith(name, "cl_eta_") || continue
        push!(etas, parse(Float64, replace(name[length("cl_eta_")+1:end], "p" => ".")))
    end
    return etas
end

function replot_pitching_wing(dir=joinpath("data", "pitching_wing"); n_loops::Integer=2)
    unsteady_csv = joinpath(dir, "pitching_wing_unsteady.csv")
    isfile(unsteady_csv) || error("Missing unsteady CSV: $(unsteady_csv). " *
        "Run the pitching-wing example to completion first.")

    header = split(strip(readline(unsteady_csv)), ',')
    section_eta = _replot_section_eta(header)
    isempty(section_eta) && error("No cl_eta_* columns found in $(unsteady_csv).")

    raw = readdlm(unsteady_csv, ',', Float64; skipstart=1)
    time = raw[:, 1]
    t_over_T = raw[:, 2]
    alpha = raw[:, 3]
    CL = raw[:, 4]
    section = permutedims(raw[:, 5:4+length(section_eta)])  # nη × ntime
    period = time[2] / t_over_T[2]

    println("Replotting from $(unsteady_csv)")
    println("  $(size(raw, 1)) samples, $(round(t_over_T[end], sigdigits=3)) cycles, " *
            "η = $(section_eta)")

    conv = plot_pitching_wing_convergence(time, period, CL, section_eta, section; path=dir)
    hyst = plot_pitching_wing_hysteresis(time, period, alpha, CL, section_eta, section;
        path=dir, n_loops)
    println("  Wrote convergence plots: $(conv.CL), $(conv.section_cl)")
    println("  Wrote hysteresis plot: $(hyst)")

    # Dynamic comparison: the unsteady c_ℓ-vs-α loop (same cycle window as the
    # hysteresis plot) overlaid on the dynamic experiment, one subplot per η.
    compare_dyn = nothing
    if isdir(PITCHING_WING_EXP_DATA_DIR)
        w = _pitching_wing_last_cycles_window(time, period, n_loops)
        compare_dyn = plot_pitching_wing_exp_comparison(section_eta, alpha[w],
            section[:, w], :dynamic;
            path=dir, filename="pitching_wing_cl_vs_alpha_dynamic.png")
        println("  Wrote dynamic comparison plot: $(compare_dyn)")
    end

    # Static polar CSV drives both the static polar plots and the static
    # comparison (static-sim vs static-experiment, over the shared α range).
    static_csv = joinpath(dir, "pitching_wing_static_polar.csv")
    static = nothing
    compare_stat = nothing
    if isfile(static_csv)
        section_eta_static = _replot_section_eta(split(strip(readline(static_csv)), ','))
        nη = length(section_eta_static)
        s = readdlm(static_csv, ',', Float64; skipstart=1)
        rows = [(; alpha_deg=s[i, 1], CL=s[i, 2], section_cl=s[i, 3:2+nη])
                for i in axes(s, 1)]
        static = plot_pitching_wing_static_polar(rows, section_eta_static; path=dir)
        println("  Wrote static polar plots: $(static.section_cl), $(static.CL)")
        if isdir(PITCHING_WING_EXP_DATA_DIR)
            compare_stat = plot_pitching_wing_exp_comparison(section_eta_static,
                s[:, 1], permutedims(s[:, 3:2+nη]), :static;
                path=dir, filename="pitching_wing_cl_vs_alpha_static.png")
            println("  Wrote static comparison plot: $(compare_stat)")
        end
    end

    compare = (; dynamic=compare_dyn, static=compare_stat)
    return (; convergence=conv, hysteresis=hyst, comparison=compare, static)
end

if abspath(PROGRAM_FILE) == @__FILE__
    dir = length(ARGS) >= 1 ? ARGS[1] : joinpath("data", "pitching_wing")
    n_loops = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2
    replot_pitching_wing(dir; n_loops)
end
