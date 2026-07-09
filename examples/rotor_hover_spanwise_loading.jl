## Spanwise loading: FLOWPanel particle-wake rotor vs CCBlade BEM (BRAINSTORM 003).
##
## Reuses the validated rotor_hover_pressure_comparison.jl setup via the
## RHPC_SETUP_ONLY hook (builds rotor + particle wake + frames + base monitors,
## skips its own simulate!), then runs simulate! with three SpanwiseLoadingMonitors
## placed after the Bernoulli, Laplace-lamb, and SurfaceVorticity force monitors.
##
## A single blade is isolated with select=cp->cp[2]>0 in the rotor frame, so the
## monitor's linear binning (span along +y, thrust along the rotation axis -x,
## torque-consuming tangential load along +z) is correct. The last revolution is
## averaged and overlaid on the CCBlade BEM
## spanwise loading (data/rotor_hover_ccblade/, produced by rotor_hover_ccblade.jl).
##
## Defaults: RPM=5400, full 10-rev / 36-step run (360 steps). For a quick local
## smoke test set e.g. NREVS=0.3 NT=12. Intended to be committed and run on a
## workstation; download data/rotor_hover_spanwise_loading/ afterwards.

# ---- configure the included setup via ENV (set before include) ----
get!(ENV, "RPM", "5400")
get!(ENV, "RHPC_SETUP_ONLY", "true")
get!(ENV, "RUN_NAME", "rotor_hover_spanwise_loading")
get!(ENV, "SAVE_VTK", "true")          # needed: monitor CSVs require a save path

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm
import CSV
using DataFrames: DataFrame
import DataFrames
import PythonPlot as plt
using Printf: @printf, @sprintf

# Build rotor + particle wake + frames + base monitors (skips its own simulate!).
include(joinpath(@__DIR__, "rotor_hover_pressure_comparison.jl"))

const NBLADES = 2
nbins = parse(Int, get(ENV, "SPANWISE_NBINS", "25"))
bem_dir = joinpath(@__DIR__, "..", "data", "rotor_hover_ccblade")
norm_rotor = pnl.RotorNormalization(rho, 2 * R, 1)

# ---- third force method: surface vorticity (Bernoulli + Laplace-lamb come from include) ----
force_monitor_surface_vorticity = pnl.SurfaceVorticityForce(rotor, length(t_range), 1;
    rho, i_frame=1, normalization=norm_rotor,
    correct_kuttacondition=p_correct_kuttacondition_flag, verbose=true)

# ---- spanwise loading monitors (one per force method) ----
# Rotor frame (i_frame=1): blade spans +y, rotation axis (=thrust dir) is -x.
# select isolates the +y blade so signed-y binning maps cleanly to radius.
one_blade = cp -> cp[2] > 0
make_span() = pnl.SpanwiseLoadingMonitor(nbins, 1;
    i_frame=1,
    span_axis=[0.0, 1.0, 0.0],
    components=(thrust=[-1.0, 0.0, 0.0], tangential=[0.0, 0.0, 1.0]),
    per_length=true,
    select=one_blade,
    verbose=false)
span_bernoulli = make_span()
span_laplace_lamb = make_span()
span_surface_vorticity = make_span()

# Each spanwise monitor must immediately follow the force monitor whose :F it reads.
monitors = (
    pressure_bernoulli,            force_monitor_bernoulli,            span_bernoulli,
    pressure_laplace_lamb,         force_monitor_laplace_lamb,         span_laplace_lamb,
    force_monitor_surface_vorticity,                                    span_surface_vorticity,
    bound_circulation,
)

println("\nSpanwise loading run: RPM=$(RPM), $(length(t_range)) steps, nbins=$(nbins)")
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    set_Das_eta_kinematic=NaN,
    set_Das_min_kinematic_displacement,
    monitors,
    body_solvers, backend, backend_wake,
    wakerow_no_hessian_to_particles,
    body_hessian_to_particles,
    body_on_wake,
    panel_wake_on_particles,
    particle_hessian_self,
    particle_relax,
    diagnose_particle_gamma,
    diagnose_particle_influence,
    diagnostic_vertical=particle_diagnostic_vertical,
    verbose=true,
    path=save_path, name,
)

# ------------------------------------------------------------------ post ------
# The spanwise monitors wrote per-step CSVs to save_path/monitors/. Average the
# last revolution (last `nt` steps) per bin to get a clean radial profile.
monitors_dir = joinpath(save_path, "monitors")
span_csvs = sort(filter(f -> occursin("_spanwise_system1.csv", f), readdir(monitors_dir)))
length(span_csvs) == 3 || error("Expected 3 spanwise monitor CSVs in $(monitors_dir); found $(length(span_csvs)): $(span_csvs)")
# Ascending monitorNN order matches the tuple order: bernoulli, laplace_lamb, surface_vorticity.
method_labels = ["bernoulli", "laplace_lamb", "surface_vorticity"]

"Average the last `nrev_steps` steps of a spanwise monitor CSV; return per-blade loads."
function average_last_rev(csv_path, nrev_steps)
    df = CSV.read(csv_path, DataFrame)
    (:thrust in propertynames(df) && :tangential in propertynames(df)) ||
        error("Spanwise CSV $(csv_path) must contain thrust and tangential columns")
    laststep = maximum(df.step)
    firststep = max(minimum(df.step), laststep - nrev_steps + 1)
    win = df[(df.step .>= firststep) .& (df.step .<= laststep), :]
    bins = sort(unique(win.bin))
    rR = Float64[]; dTdr = Float64[]; dFtdr = Float64[]
    for b in bins
        sub = win[win.bin .== b, :]
        push!(rR, sum(sub.bin_center) / DataFrames.nrow(sub) / R)
        push!(dTdr, sum(sub.thrust) / DataFrames.nrow(sub))   # per-blade dT/dr [N/m]
        push!(dFtdr, sum(sub.tangential) / DataFrames.nrow(sub))   # per-blade dFt/dr [N/m]
    end
    return (; r_over_R=rR, dTdr, dFtdr)
end

function signed_tangential(rR, dFtdr)
    r = rR .* R
    torque = 0.0
    for i in 1:length(r)-1
        torque += 0.5 * (r[i] * dFtdr[i] + r[i+1] * dFtdr[i+1]) * (r[i+1] - r[i])
    end
    return torque < 0 ? -dFtdr : dFtdr
end

panel = Dict{String,NamedTuple}()
for (label, csv) in zip(method_labels, span_csvs)
    avg = average_last_rev(joinpath(monitors_dir, csv), nt)
    # Sign convention: report dT/dr > 0 for thrust. Flip if the net integral is negative.
    dTdr = sum(avg.dTdr) < 0 ? -avg.dTdr : avg.dTdr
    dFtdr = signed_tangential(avg.r_over_R, avg.dFtdr)
    panel[label] = (; r_over_R=avg.r_over_R, dTdr, dFtdr)
end

# consolidated panel CSV (total = both-blade = per-blade x NBLADES, to match BEM dTdr_total)
ref_rR = panel["bernoulli"].r_over_R
panel_csv = DataFrame(r_over_R = ref_rR,
    dTdr_total_bernoulli         = NBLADES .* panel["bernoulli"].dTdr,
    dTdr_total_laplace_lamb      = NBLADES .* panel["laplace_lamb"].dTdr,
    dTdr_total_surface_vorticity = NBLADES .* panel["surface_vorticity"].dTdr,
    dFtdr_total_bernoulli         = NBLADES .* panel["bernoulli"].dFtdr,
    dFtdr_total_laplace_lamb      = NBLADES .* panel["laplace_lamb"].dFtdr,
    dFtdr_total_surface_vorticity = NBLADES .* panel["surface_vorticity"].dFtdr)
panel_csv_path = joinpath(save_path, "rotor_hover_panel_spanwise_loading.csv")
CSV.write(panel_csv_path, panel_csv)
println("\nWrote $(panel_csv_path)")

# integrated-thrust sanity check: ∫ dTdr_total dr ≈ T  (vs CT·ρ·n²·D⁴)
n_rev = RPM / 60
ct_ref = rho * n_rev^2 * (2 * R)^4
for label in method_labels
    rR = panel[label].r_over_R
    dTdr = panel[label].dTdr
    r = rR .* R
    T = 0.0
    for i in 1:length(r)-1
        T += 0.5 * (NBLADES*dTdr[i] + NBLADES*dTdr[i+1]) * (r[i+1] - r[i])
    end
    @printf("  panel %-18s  integrated T = %7.4f N  -> CT = %.5f\n", label, T, T/ct_ref)
end

# ---- overlay plot: BEM (per ncrit) + panel (3 methods) ----
fig, ax = plt.subplots(figsize=(7, 5))
bem_files = isdir(bem_dir) ? sort(filter(f -> occursin("rotor_hover_ccblade_sectional_", f) && endswith(f, ".csv"), readdir(bem_dir))) : String[]
for f in bem_files
    df = CSV.read(joinpath(bem_dir, f), DataFrame)
    polar = replace(replace(f, "rotor_hover_ccblade_sectional_" => ""), ".csv" => "")
    ax.plot(df.r_over_R, df.dTdr_total, "-"; linewidth=1.0, alpha=0.7,
            label="BEM $(polar)")
end
isempty(bem_files) && @warn "No CCBlade sectional CSVs found in $(bem_dir); run rotor_hover_ccblade.jl first for the overlay."
for (label, ls) in zip(method_labels, ("-o", "-s", "-^"))
    rR = panel[label].r_over_R
    dTdr = panel[label].dTdr
    ax.plot(rR, NBLADES .* dTdr, ls; linewidth=2.0, markersize=3, label="panel $(label)")
end
ax.set_xlabel("r/R"); ax.set_ylabel("dT/dr (both blades) [N/m]")
ax.set_title(@sprintf("Spanwise loading: panel vs BEM, RPM=%.0f, hover", RPM))
ax.grid(true); ax.legend(fontsize=7)
fig.tight_layout()
overlay_path = joinpath(save_path, "spanwise_loading_panel_vs_bem.png")
fig.savefig(overlay_path, dpi=150)
plt.close()
println("Wrote $(overlay_path)")
