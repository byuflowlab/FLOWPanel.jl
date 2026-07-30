## BRAINSTORM item 003 follow-up: diagnostics for the inviscid BEM thrust
## overprediction and the loading spike at r/R~0.79.
##
## Reads the CSVs already written by examples/rotor_hover_ccblade.jl (no XFOIL
## or CCBlade rerun) and produces radial diagnostic plots plus a summary table.
## Established diagnosis this script visualizes:
##   (a) the spike at r/R=0.785 is a spurious high-induction BEM root
##       (u=97.7 m/s vs ~5 m/s at neighbors) worth ~0.022 of CT;
##   (b) the remaining inviscid excess over the viscous polar sets is the
##       inviscid cl level itself (no viscous decambering at Re 10k-53k).

import CSV
import DataFrames
using DataFrames: DataFrame
import PythonPlot as plt
import FLOWMath
using Printf: @printf, @sprintf

const DATA_DIR = get(ENV, "DATA_DIR",
    joinpath(@__DIR__, "..", "data", "rotor_hover_ccblade"))
const RADIAL_CSV = joinpath(DATA_DIR, "rotor_hover_ccblade_radial.csv")
const POLARS_CSV = joinpath(DATA_DIR, "rotor_hover_ccblade_polars.csv")

# Operating point (must match rotor_hover_ccblade.jl)
const RPM  = parse(Float64, get(ENV, "RPM", "5400"))
const R    = 0.119
const B    = 2
const RHO  = 1.179
const NREV = RPM / 60
const D    = 2R
const CT_PANEL = 0.0506
const CT_EXPERIMENT = 0.072

polar_legend(ncrit::Integer) = ncrit == 0 ? "inviscid" : "ncrit=$(ncrit)"
polar_label(ncrit::Integer) = ncrit == 0 ? "inviscid" : "ncrit$(ncrit)"
polar_legend(label::AbstractString) = startswith(label, "ncrit") ?
    replace(label, "ncrit" => "ncrit=") : label

function ensure_polarset!(df)
    :polarset in propertynames(df) && return df
    df.polarset = [polar_label(ncrit) for ncrit in df.ncrit]
    return df
end

cumtrapz(x, y) = [i == 1 ? 0.0 :
                  sum(0.5 * (y[j] + y[j+1]) * (x[j+1] - x[j]) for j in 1:i-1)
                  for i in eachindex(x)]
trapz(x, y) = cumtrapz(x, y)[end]

# ------------------------------------------------------------------ load ------
radial = CSV.read(RADIAL_CSV, DataFrame)
polars = CSV.read(POLARS_CSV, DataFrame)
ensure_polarset!(radial)
ensure_polarset!(polars)
polarsets = sort(unique(radial.polarset))
Vc_hover = radial.Vc[argmin(abs.(radial.Vc))]
hover(polarset) = sort(radial[(radial.polarset .== polarset) .& (radial.Vc .== Vc_hover), :],
                       :r_over_R)
println("polar sets: $(polarsets); hover Vc = $(Vc_hover) m/s")

# reference viscous set for the cl-gap comparison (prefer ncrit9)
ref_candidates = [p for p in reverse(polarsets) if p != "inviscid"]
REF_POLARSET = "ncrit9" in polarsets ? "ncrit9" :
    (isempty(ref_candidates) ? polarsets[end] : ref_candidates[1])

# ------------------------------------------------- simple radial plots --------
function plot_radial(col, ylabel, fname; yscale=nothing)
    fig, ax = plt.subplots(figsize=(6, 4))
    for polarset in polarsets
        sub = hover(polarset)
        ax.plot(sub.r_over_R, sub[!, col], "-o"; markersize=3,
                label=polar_legend(polarset))
    end
    isnothing(yscale) || ax.set_yscale(yscale)
    ax.set_xlabel("r/R"); ax.set_ylabel(ylabel)
    ax.set_title(@sprintf("hover (Vc=%.4f m/s)", Vc_hover))
    ax.grid(true); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(joinpath(DATA_DIR, fname), dpi=150)
    plt.close()
    println("Wrote $(fname)")
end

plot_radial(:alpha_deg, "Angle of attack alpha [deg]", "radial_alpha.png")
plot_radial(:W,         "Resultant velocity W [m/s]",  "radial_W.png")
plot_radial(:cl,        "Section cl",                  "radial_cl.png")
plot_radial(:cd,        "Section cd",                  "radial_cd.png")

# -------------------------------------------- cl at fixed (reference) alpha ---
# For each section, evaluate every polar set's cl at the *same* alpha (the
# REF_POLARSET hover operating alpha). Differences are pure polar lift level,
# independent of the BEM inflow solution.
ref = hover(REF_POLARSET)
let fig_ax = plt.subplots(figsize=(6, 4))
    fig, ax = fig_ax
    for polarset in polarsets
        cl_at = Float64[]
        for (i, roR) in enumerate(ref.r_over_R)
            p = sort(polars[(polars.polarset .== polarset) .&
                            (isapprox.(polars.r_over_R, roR; atol=1e-9)), :], :alpha_deg)
            push!(cl_at, FLOWMath.linear(p.alpha_deg, p.cl, ref.alpha_deg[i]))
        end
        ax.plot(ref.r_over_R, cl_at, "-o"; markersize=3, label=polar_legend(polarset))
    end
    ax.set_xlabel("r/R")
    ax.set_ylabel("cl at $(polar_legend(REF_POLARSET)) operating alpha")
    ax.set_title("Polar lift level at identical alpha")
    ax.grid(true); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(joinpath(DATA_DIR, "radial_cl_gap.png"), dpi=150)
    plt.close()
    println("Wrote radial_cl_gap.png")
end

# ----------------------------------- operating point on polar (spike zone) ----
# Sections 17-19 bracket the inviscid spike (r/R = 0.745, 0.785, 0.826). The
# spike section's operating point sits on a smooth part of its polar -> the
# anomaly is in the inflow solution, not the section data.
spike_sections = [s for s in (17, 18, 19) if s in polars.section]
let fig_axs = plt.subplots(1, length(spike_sections); figsize=(4.2 * length(spike_sections), 4))
    fig, axs = fig_axs
    for (k, sec) in enumerate(spike_sections)
        ax = length(spike_sections) == 1 ? axs : axs[k - 1]
        roR = first(polars[polars.section .== sec, :r_over_R])
        for polarset in polarsets
            p = sort(polars[(polars.polarset .== polarset) .& (polars.section .== sec), :],
                     :alpha_deg)
            line = ax.plot(p.alpha_deg, p.cl, "-"; linewidth=1.0,
                           label=polar_legend(polarset))
            color = line[0].get_color()
            op = hover(polarset)[sec, :]
            ax.plot([op.alpha_deg], [op.cl], "*"; color, markersize=13,
                    markeredgecolor="k", markeredgewidth=0.4)
        end
        ax.set_title(@sprintf("section %d, r/R=%.3f", sec, roR))
        ax.set_xlabel("alpha [deg]"); ax.set_ylabel("cl")
        ax.grid(true); k == 1 && ax.legend(fontsize=7, title="star = op point")
    end
    fig.tight_layout()
    fig.savefig(joinpath(DATA_DIR, "polars_oppoint_spike.png"), dpi=150)
    plt.close()
    println("Wrote polars_oppoint_spike.png")
end

# ---------------------------------------------------- cumulative thrust -------
despike(y) = begin
    y2 = copy(y)
    i = argmax(y2)
    if 1 < i < length(y2)
        y2[i] = 0.5 * (y2[i-1] + y2[i+1])
    end
    y2
end

ct_norm = RHO * NREV^2 * D^4
let fig_ax = plt.subplots(figsize=(6, 4.5))
    fig, ax = fig_ax
    for polarset in polarsets
        sub = hover(polarset)
        cum = cumtrapz(sub.r, B .* sub.Np) ./ ct_norm
        ax.plot(sub.r_over_R, cum, "-o"; markersize=3, label=polar_legend(polarset))
        if polarset == "inviscid"
            cum2 = cumtrapz(sub.r, despike(B .* sub.Np)) ./ ct_norm
            ax.plot(sub.r_over_R, cum2, "--"; color="k", linewidth=1.2,
                    label="inviscid despiked")
        end
    end
    ax.axhline(CT_PANEL; color="k", linestyle=":", linewidth=1,
               label=@sprintf("panel %.4f", CT_PANEL))
    ax.axhline(CT_EXPERIMENT; color="gray", linestyle=":", linewidth=1,
               label=@sprintf("experiment %.3f", CT_EXPERIMENT))
    ax.set_xlabel("r/R"); ax.set_ylabel("cumulative CT(r)")
    ax.set_title("Cumulative thrust buildup (hover)")
    ax.grid(true); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(joinpath(DATA_DIR, "radial_cumCT.png"), dpi=150)
    plt.close()
    println("Wrote radial_cumCT.png")
end

# ------------------------------------------------ momentum consistency --------
# Blade-element loading B*Np vs annulus momentum 4*pi*rho*r*u^2 (Vx->0 limit).
# A converged physical BEM solution satisfies these approximately (up to
# tip-loss F); the spurious root's u=97.7 m/s throws the momentum curve off by
# ~2 orders of magnitude at one station.
let fig_ax = plt.subplots(figsize=(6.5, 4.5))
    fig, ax = fig_ax
    for polarset in ("inviscid", REF_POLARSET)
        polarset in polarsets || continue
        sub = hover(polarset)
        line = ax.plot(sub.r_over_R, B .* sub.Np, "-o"; markersize=3,
                       label="blade element B*Np, $(polar_legend(polarset))")
        color = line[0].get_color()
        ax.plot(sub.r_over_R, 4pi .* RHO .* sub.r .* sub.u .^ 2, "--x";
                color, markersize=4,
                label="momentum 4*pi*rho*r*u^2, $(polar_legend(polarset))")
    end
    ax.set_yscale("log")
    ax.set_xlabel("r/R"); ax.set_ylabel("dT/dr [N/m]")
    ax.set_title("Blade-element vs annulus-momentum loading (hover)")
    ax.grid(true, which="both", alpha=0.4); ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(joinpath(DATA_DIR, "radial_momentum_check.png"), dpi=150)
    plt.close()
    println("Wrote radial_momentum_check.png")
end

# ------------------------------------------------------------- summary --------
println("\n" * "="^70)
println("Integrated hover CT (trapezoid of B*Np):")
for polarset in polarsets
    sub = hover(polarset)
    ct = trapz(sub.r, B .* sub.Np) / ct_norm
    if polarset == "inviscid"
        ct2 = trapz(sub.r, despike(B .* sub.Np)) / ct_norm
        i = argmax(sub.Np)
        @printf("  %-9s: CT = %.4f  (despiked %.4f; spike at r/R=%.3f contributes %.4f)\n",
                polar_legend(polarset), ct, ct2, sub.r_over_R[i], ct - ct2)
    else
        @printf("  %-9s: CT = %.4f\n", polar_legend(polarset), ct)
    end
end
@printf("  reference: panel %.4f, experiment %.4f\n", CT_PANEL, CT_EXPERIMENT)

if "inviscid" in polarsets
    println("\nInviscid / $(polar_legend(REF_POLARSET)) cl ratio at the $(polar_legend(REF_POLARSET)) operating alpha:")
    for (i, roR) in enumerate(ref.r_over_R)
        p0 = sort(polars[(polars.polarset .== "inviscid") .&
                         (isapprox.(polars.r_over_R, roR; atol=1e-9)), :], :alpha_deg)
        cl0 = FLOWMath.linear(p0.alpha_deg, p0.cl, ref.alpha_deg[i])
        @printf("  r/R=%.3f  alpha=%6.2f deg  cl_inviscid=%.3f  cl_%s=%.3f  ratio=%.2f\n",
                roR, ref.alpha_deg[i], cl0, REF_POLARSET, ref.cl[i], cl0 / ref.cl[i])
    end
end
println("="^70)
