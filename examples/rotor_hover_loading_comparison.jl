## Rotor hover radial loading comparison: CCBlade vs FLOWPanel replay stats.
##
## Defaults compare inviscid CCBlade sectional loading against the replayed
## FLOWPanel pressure-loading stats from rotor_hover_pressure_comparison.
##
## Loading is plotted NONdimensionally as dCT/d(r/R) = (dT/dr) R / (rho n^2 D^4)
## with n taken from each dataset's own RPM (BEM_RPM / PANEL_RPM), so datasets
## at different operating points cannot be compared apples-to-oranges silently.
## Each legend entry carries its curve's integrated CT as a built-in sanity
## check against the known integrated values (panel ~0.0506 converged at
## RPM=5400 on the convergence mesh, inviscid BEM ~0.113 including the spurious
## r/R=0.785 spike). PANEL_RPM defaults to the omega recorded in the replay
## metadata TOML next to the panel run (6000 RPM for
## rotor_hover_pressure_comparison).

import CSV
import DataFrames
using DataFrames: DataFrame
import PythonPlot as plt
import TOML
using Printf: @sprintf

const DEFAULT_BEM_CSV = joinpath(@__DIR__, "..", "data", "rotor_hover_ccblade",
    "rotor_hover_ccblade_sectional_inviscid.csv")
const DEFAULT_PANEL_CSV = joinpath(@__DIR__, "..", "data", "rotor_hover_pressure_comparison",
    "spanwise_loading_replay", "rotor_hover_pressure_comparison_spanwise_loading_stats.csv")
const DEFAULT_OUT_DIR = joinpath(@__DIR__, "..", "data", "rotor_hover_ccblade")

const RHO = parse(Float64, get(ENV, "RHO", "1.179"))
const R = parse(Float64, get(ENV, "R", "0.119"))
const D = 2R
const BEM_RPM = parse(Float64, get(ENV, "BEM_RPM", "5400"))

# Fewer samples than this in the panel stats window means the medians come from
# a short (likely transient) window rather than converged hover.
const MIN_WINDOW_SAMPLES = parse(Int, get(ENV, "MIN_WINDOW_SAMPLES", "20"))

function env_list(name)
    s = strip(get(ENV, name, ""))
    isempty(s) && return String[]
    return [strip(x) for x in split(s, ",") if !isempty(strip(x))]
end

function polar_label(path)
    base = basename(path)
    stem = replace(base, ".csv" => "")
    label = replace(stem, "rotor_hover_ccblade_sectional_" => "")
    startswith(label, "ncrit") && return replace(label, "ncrit" => "ncrit=")
    return label
end

function require_columns(df, cols, path)
    missing = [c for c in cols if !(c in propertynames(df))]
    isempty(missing) || error("$(path) is missing required columns: $(missing)")
    return nothing
end

"Infer the panel run's RPM from the replay metadata TOML two levels above the stats CSV."
function infer_panel_rpm(panel_csv)
    haskey(ENV, "PANEL_RPM") && return parse(Float64, ENV["PANEL_RPM"])
    run_dir = dirname(dirname(abspath(panel_csv)))
    candidates = filter(f -> endswith(f, ".metadata.toml"), readdir(run_dir; join=true))
    for path in candidates
        metadata = TOML.parsefile(path)
        haskey(metadata, "step") || continue
        for step in metadata["step"], frame in get(step, "frame", [])
            omega = Float64(get(frame, "omega", 0.0))
            omega > 0 && return omega * 60 / (2π)
        end
    end
    error("Could not infer panel RPM from metadata in $(run_dir); set PANEL_RPM explicitly.")
end

ct_norm(rpm) = RHO * (rpm / 60)^2 * D^4

trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i in 1:length(x)-1)

function load_bem(path)
    df = CSV.read(path, DataFrame)
    require_columns(df, [:r_over_R, :dTdr_total, :Tp], path)
    return df
end

function plot_bem_lines!(ax, bem_csv, extra_bem_csvs, component)
    norm_bem = ct_norm(BEM_RPM)
    for (i, path) in enumerate((bem_csv, extra_bem_csvs...))
        df = load_bem(path)
        color = i == 1 ? "k" : nothing
        linewidth = i == 1 ? 1.7 : 1.1
        alpha = i == 1 ? 1.0 : 0.75
        y = (component == :normal ? df.dTdr_total : 2.0 .* df.Tp) .* R ./ norm_bem
        ct = trapz(df.r_over_R, y)
        label = @sprintf("CCBlade %s, RPM=%.0f (int %.4f)", polar_label(path), BEM_RPM, ct)
        ax.plot(df.r_over_R, y, "-"; color, linewidth, alpha, label)
    end
end

function panel_component_available(stats, component)
    prefix = component == :normal ? "dTdr" : "dFtdr"
    return Symbol("median_$(prefix)_total_equiv") in propertynames(stats)
end

function plot_panel_source!(ax, stats, source, component, panel_rpm)
    prefix = component == :normal ? "dTdr" : "dFtdr"
    med_col = Symbol("median_$(prefix)_total_equiv")
    min_col = Symbol("min_$(prefix)_total_equiv")
    max_col = Symbol("max_$(prefix)_total_equiv")
    q25_col = Symbol("q25_$(prefix)_total_equiv")
    q75_col = Symbol("q75_$(prefix)_total_equiv")

    sub = stats[stats.source .== String(source), :]
    DataFrames.nrow(sub) == 0 && (@warn "No panel rows for source $(source)"; return nothing)
    scale = R / ct_norm(panel_rpm)

    blade_colors = Dict(1 => "tab:blue", 2 => "tab:orange")
    blade_styles = Dict(1 => "-", 2 => "--")
    blade_bar_styles = Dict(1 => "solid", 2 => "dashed")
    for blade in sort(unique(sub.blade))
        bdf = sub[sub.blade .== blade, :]
        color = get(blade_colors, blade, nothing)
        linestyle = get(blade_styles, blade, "-")
        barstyle = get(blade_bar_styles, blade, "solid")
        finite = isfinite.(bdf[!, med_col])
        ct = any(finite) ? trapz(bdf.r_over_R[finite], scale .* bdf[!, med_col][finite]) : NaN
        ax.plot(bdf.r_over_R, scale .* bdf[!, med_col]; color, linestyle, marker="o",
            linewidth=1.8, markersize=3.5,
            label=@sprintf("%s blade %d median, RPM=%.0f (int %.4f)",
                           source, blade, panel_rpm, ct))
        ax.vlines(bdf.r_over_R, scale .* bdf[!, min_col], scale .* bdf[!, max_col];
            color, linestyles=barstyle, linewidth=0.6, alpha=0.45)
        ax.vlines(bdf.r_over_R, scale .* bdf[!, q25_col], scale .* bdf[!, q75_col];
            color, linestyles=barstyle, linewidth=2.0, alpha=0.8)
    end
    return nothing
end

function panel_window_note(stats)
    ns = maximum(skipmissing(stats.n_samples))
    ns >= MIN_WINDOW_SAMPLES && return @sprintf("panel window: %d samples", ns)
    note = @sprintf("panel window: %d samples — TRANSIENT, not converged hover", ns)
    @warn "Panel stats CSV medians come from only $(ns) samples; treat the panel curve " *
          "as a transient snapshot, not converged hover loading."
    return note
end

function plot_component(stats, bem_csv, extra_bem_csvs, out_dir, panel_sources, component,
                        panel_rpm)
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    plot_bem_lines!(ax, bem_csv, extra_bem_csvs, component)

    has_panel = panel_component_available(stats, component)
    if has_panel
        for source in panel_sources
            plot_panel_source!(ax, stats, source, component, panel_rpm)
        end
    elseif component == :tangential
        @warn "Panel CSV lacks dFtdr columns; writing tangential plot with BEM only."
    else
        @warn "Panel CSV lacks dTdr columns; writing normal plot with BEM only."
    end

    if component == :normal
        ylabel = "dCT/d(r/R) = (dT/dr) R / (rho n^2 D^4)"
        title = "Rotor hover normal loading"
        filename = "rotor_hover_loading_comparison_normal.png"
    else
        ylabel = "dCFt/d(r/R) = (dFt/dr) R / (rho n^2 D^4)"
        title = "Rotor hover tangential loading"
        filename = "rotor_hover_loading_comparison_tangential.png"
    end
    ax.set_xlabel("r/R")
    ax.set_ylabel(ylabel)
    ax.set_title("$(title)\n$(panel_window_note(stats))", fontsize=10)
    ax.grid(true, alpha=0.35)
    ax.legend(fontsize=7)
    fig.tight_layout()
    path = joinpath(out_dir, filename)
    fig.savefig(path, dpi=170)
    plt.close()
    return path
end

function plot_loading_comparison(; bem_csv=get(ENV, "CCBLADE_CSV", DEFAULT_BEM_CSV),
                                  panel_csv=get(ENV, "PANEL_CSV", DEFAULT_PANEL_CSV),
                                  out_dir=get(ENV, "OUT_DIR", DEFAULT_OUT_DIR),
                                  panel_sources=("laplace_lamb", "bernoulli"),
                                  extra_bem_csvs=env_list("EXTRA_BEM_CSVS"))
    isfile(bem_csv) || error("BEM CSV not found: $(bem_csv)")
    isfile(panel_csv) || error("Panel CSV not found: $(panel_csv)")
    mkpath(out_dir)

    panel_rpm = infer_panel_rpm(panel_csv)
    println("BEM RPM = $(BEM_RPM), panel RPM = $(panel_rpm)")

    stats = CSV.read(panel_csv, DataFrame)
    require_columns(stats, [:source, :blade, :r_over_R, :n_samples], panel_csv)
    normal_path = plot_component(stats, bem_csv, extra_bem_csvs, out_dir, panel_sources,
        :normal, panel_rpm)
    tangential_path = plot_component(stats, bem_csv, extra_bem_csvs, out_dir, panel_sources,
        :tangential, panel_rpm)
    println("Wrote $(normal_path)")
    println("Wrote $(tangential_path)")
    return (; normal_path, tangential_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_loading_comparison()
end
