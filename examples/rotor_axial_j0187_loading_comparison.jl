## Compare final-cycle panel dCT/d(r/R) with CCBlade ncrit=4 and 9 at J=0.1867.

import CSV
using DataFrames: DataFrame, nrow
import PythonPlot as plt
using Printf: @sprintf

const RUN_DIR = get(ENV, "RUN_DIR", joinpath(@__DIR__, "..", "data", "rotor_axial_j0187_ccblade"))
const STATS_CSV = get(ENV, "PANEL_STATS_CSV", joinpath(RUN_DIR, "spanwise_loading_replay", "rotor_axial_j0187_ccblade_spanwise_loading_stats.csv"))
const BEM_DIR = get(ENV, "BEM_DIR", RUN_DIR)
const OUT_DIR = get(ENV, "OUT_DIR", joinpath(RUN_DIR, "comparison"))
const RHO = parse(Float64, get(ENV, "RHO", "1.179")); const R = 0.119; const RPM = 5400.0
const CT_SCALE = RHO * (RPM / 60)^2 * (2R)^4
mkpath(OUT_DIR)
trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i in 1:length(x)-1)

stats = CSV.read(STATS_CSV, DataFrame)
panel = stats[stats.source .== get(ENV, "PANEL_SOURCE", "laplace_lamb"), :]
nrow(panel) > 0 || error("No selected panel source in $(STATS_CSV)")
# Each blade is stored as a total-equivalent curve; average blade realizations.
rs = sort(unique(panel.r_over_R))
med = [sum(panel[(panel.r_over_R .== r), :].median_dTdr_total_equiv) / count(==(r), panel.r_over_R) for r in rs]
q25 = [sum(panel[(panel.r_over_R .== r), :].q25_dTdr_total_equiv) / count(==(r), panel.r_over_R) for r in rs]
q75 = [sum(panel[(panel.r_over_R .== r), :].q75_dTdr_total_equiv) / count(==(r), panel.r_over_R) for r in rs]
panel_curve = med .* R ./ CT_SCALE
panel_ct = trapz(rs, panel_curve)

rows = [(case="panel_cycle_mean", CT=panel_ct, relative_to_panel=0.0)]
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.plot(rs, panel_curve, "o-"; color="tab:blue", label=@sprintf("panel cycle mean (CT %.5f)", panel_ct))
ax.vlines(rs, q25 .* R ./ CT_SCALE, q75 .* R ./ CT_SCALE; color="tab:blue", linewidth=2.2, alpha=0.7, label="panel q25–q75")
for (label, color) in (("ncrit4", "black"), ("ncrit9", "tab:orange"))
    path = joinpath(BEM_DIR, "rotor_hover_ccblade_sectional_$(label)_Vc4_J0p1867.csv")
    bem = CSV.read(path, DataFrame)
    curve = bem.dTdr_total .* R ./ CT_SCALE
    ct = trapz(bem.r_over_R, curve)
    push!(rows, (case="ccblade_$(label)", CT=ct, relative_to_panel=(ct - panel_ct) / panel_ct))
    ax.plot(bem.r_over_R, curve, "-"; color, linewidth=1.7, label=@sprintf("CCBlade %s (CT %.5f)", replace(label, "ncrit" => "ncrit="), ct))
end
ax.set_xlabel("r/R"); ax.set_ylabel("dCT/d(r/R)"); ax.set_title("DJI 9443 axial flow: Vc=4 m/s, J=0.1867")
ax.grid(true, alpha=0.35); ax.legend(fontsize=8); fig.tight_layout()
fig.savefig(joinpath(OUT_DIR, "spanwise_total_thrust_dCTdr.png"), dpi=170); plt.close()
CSV.write(joinpath(OUT_DIR, "CT_comparison.csv"), DataFrame(rows))
fig, ax = plt.subplots(figsize=(5.5, 4.0)); ax.bar([r.case for r in rows], [r.CT for r in rows])
ax.set_ylabel("CT"); ax.set_title("CT comparison, Vc=4 m/s / J=0.1867"); fig.tight_layout()
fig.savefig(joinpath(OUT_DIR, "CT_comparison.png"), dpi=170); plt.close()
