## BRAINSTORM item 003 extension: compare the hover BEM prediction at baseline
## (physical, Re ~ 10k-53k) vs 10x Reynolds section polars, and report XFOIL
## solves that failed to converge in either run.
##
## CSV-only: reads the outputs of examples/rotor_hover_ccblade.jl for the
## baseline and for the RE_SCALE=10 SUFFIX=_re10x rerun. Failures are derived
## uniformly for both runs as the attempted alpha grid (-12:0.5:20, 65 angles)
## minus the converged angles recorded in each polars CSV; for the re10x run
## this reproduces the directly logged rotor_hover_ccblade_xfoil_failures_re10x.csv.

import CSV
using DataFrames: DataFrame
import PythonPlot as plt
using Printf: @printf, @sprintf

const DATA_DIR = get(ENV, "DATA_DIR",
    joinpath(@__DIR__, "..", "data", "rotor_hover_ccblade"))
const SUFFIX = get(ENV, "SUFFIX", "_re10x")
const RUNS = ("baseline" => "", "re10x" => SUFFIX)

const RPM = parse(Float64, get(ENV, "RPM", "5400"))
const R = 0.119
const B = 2
const RHO = 1.179
const D = 2R
const CT_NORM = RHO * (RPM / 60)^2 * D^4
const CT_PANEL = 0.0506
const CT_EXPERIMENT = 0.072

# attempted XFOIL alpha grid in xfoil_polar (both sweep branches combined)
const ATTEMPTED = sort!(unique!(vcat(collect(0.0:0.5:20.0), collect(0.0:-0.5:-12.0))))

polar_legend(ncrit::Integer) = ncrit == 0 ? "inviscid" : "ncrit=$(ncrit)"
polar_label(ncrit::Integer) = ncrit == 0 ? "inviscid" : "ncrit$(ncrit)"
polar_legend(label::AbstractString) = startswith(label, "ncrit") ?
    replace(label, "ncrit" => "ncrit=") : label

function ensure_polarset!(df)
    :polarset in propertynames(df) && return df
    df.polarset = [polar_label(ncrit) for ncrit in df.ncrit]
    return df
end

trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i in 1:length(x)-1)

function load_run(suffix)
    radial = CSV.read(joinpath(DATA_DIR, "rotor_hover_ccblade_radial$(suffix).csv"), DataFrame)
    polars = CSV.read(joinpath(DATA_DIR, "rotor_hover_ccblade_polars$(suffix).csv"), DataFrame)
    ensure_polarset!(radial)
    ensure_polarset!(polars)
    return radial, polars
end

hover_rows(radial, polarset) = begin
    Vc_hover = radial.Vc[argmin(abs.(radial.Vc))]
    sort(radial[(radial.polarset .== polarset) .& (radial.Vc .== Vc_hover), :], :r_over_R)
end

data = Dict(name => load_run(suffix) for (name, suffix) in RUNS)

# ------------------------------------------------------------- CT table -------
println("Integrated hover CT (trapezoid of B*Np), RPM=$(RPM):")
println("  polar set   baseline    re10x")
re10x_polarsets = sort(unique(data["re10x"][1].polarset))
for polarset in sort(unique(data["baseline"][1].polarset))
    cts = map(("baseline", "re10x")) do name
        radial = data[name][1]
        polarset in radial.polarset || return NaN
        sub = hover_rows(radial, polarset)
        trapz(sub.r, B .* sub.Np) / CT_NORM
    end
    @printf("  %-24s   %8.4f   %8.4f\n", polar_legend(polarset), cts...)
end
@printf("  reference: panel %.4f, experiment %.4f\n", CT_PANEL, CT_EXPERIMENT)

# ------------------------------------------------- spanwise loading overlay ---
fig, ax = plt.subplots(figsize=(7, 4.6))
colors = Dict(zip(re10x_polarsets, ("tab:blue", "tab:orange", "tab:green", "tab:red",
                                    "tab:purple", "tab:brown", "tab:pink", "tab:olive")))
for polarset in re10x_polarsets
    for (name, style, lw) in (("baseline", "-", 1.2), ("re10x", "--", 1.6))
        radial = data[name][1]
        polarset in radial.polarset || continue
        sub = hover_rows(radial, polarset)
        ax.plot(sub.r_over_R, B .* sub.Np, style; color=colors[polarset], linewidth=lw,
                label="$(polar_legend(polarset)) $(name)")
    end
end
ax.set_xlabel("r/R"); ax.set_ylabel("dT/dr (both blades) [N/m]")
ax.set_title(@sprintf("BEM hover loading: baseline vs 10x Reynolds (RPM=%.0f)", RPM))
ax.grid(true, alpha=0.4); ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
overlay_png = joinpath(DATA_DIR, "rotor_hover_ccblade_loading_re10x_vs_baseline.png")
fig.savefig(overlay_png, dpi=150)
plt.close()
println("\nWrote $(overlay_png)")

# ------------------------------------------- XFOIL failure report (both) ------
report = DataFrame(run=String[], polarset=String[], ncrit=Int[], section=Int[],
                   r_over_R=Float64[], Re=Float64[], alpha_deg=Float64[])
for (name, _) in RUNS
    polars = data[name][2]
    for polarset in sort(unique(polars.polarset)), sec in sort(unique(polars.section))
        sub = polars[(polars.polarset .== polarset) .& (polars.section .== sec), :]
        isempty(sub.alpha_deg) && continue
        recorded = Set(round.(sub.alpha_deg, digits=3))
        for a in ATTEMPTED
            round(a, digits=3) in recorded && continue
            push!(report, (name, polarset, sub.ncrit[1], sec, sub.r_over_R[1], sub.Re[1], a))
        end
    end
end
report_csv = joinpath(DATA_DIR, "rotor_hover_ccblade_xfoil_failures_report.csv")
CSV.write(report_csv, report)
println("Wrote $(report_csv)")

println("\nUnconverged XFOIL solves (out of $(length(ATTEMPTED)) attempted alphas/section, 23 sections):")
for (name, _) in RUNS
    println("  $(name):")
    sub_run = report[report.run .== name, :]
    for polarset in sort(unique(data[name][2].polarset))
        sub = sub_run[sub_run.polarset .== polarset, :]
        if isempty(sub.section)
            @printf("    %-24s: all converged\n", polar_legend(polarset))
        else
            secs = join(["$(s)@r/R=$(round(sub.r_over_R[findfirst(==(s), sub.section)], digits=2)) (x$(count(==(s), sub.section)))"
                         for s in sort(unique(sub.section))], ", ")
            @printf("    %-24s: %3d failed; sections: %s\n",
                    polar_legend(polarset), size(sub, 1), secs)
        end
    end
end
