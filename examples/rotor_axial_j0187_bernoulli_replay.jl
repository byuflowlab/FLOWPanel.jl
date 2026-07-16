## Bernoulli-only replay of the axial J=0.1867 comparison run.
##
## The original HPC run computed steady PressureBernoulli with the ef1fe1e
## inertial kinetic-energy regression (first-order blade loading cancelled).
## This script replays the saved VTK state through the corrected steady
## Bernoulli monitor only (no PressureLaplace, no CG solves) and produces:
##   1. CT vs revolution with the run's Laplace(Du/Dt) history for context and
##      CCBlade ncrit=4/9 total CT as horizontal lines;
##   2. spanwise dCT/d(r/R) averaged over the final NAVG_SAMPLES replay samples
##      (default 72 = final 2 revolutions at 36 steps/rev), overlaying both
##      CCBlade sectional curves;
##   3. a saved-vs-replayed forensic CSV identifying which Bernoulli
##      formulation produced the original run output.
##
## Read-only with respect to the saved run: it never re-solves, re-sheds, or
## overwrites anything under SAVE_PATH outside OUT_DIR.

# Load the shared replay helpers (SpanwiseReplayBinning/Snapshot, metadata and
# blade-geometry inference, flatness stats, stats CSV writer) without running
# the two-source replay driver itself.
ENV["SPANWISE_REPLAY_RUN"] = "false"
include(joinpath(@__DIR__, "rotor_hover_spanwise_loading_replay.jl"))

# ------------------------------------------------------------------ config --
run_name = get(ENV, "RUN_NAME", "rotor_axial_j0187_ccblade")
save_path = get(ENV, "SAVE_PATH", joinpath("data", "rotor_axial_j0187_ccblade"))
out_dir = get(ENV, "OUT_DIR", joinpath(save_path, "bernoulli_replay"))
metadata = TOML.parsefile(metadata_path(save_path, run_name))

R = env_float("R", something(infer_metadata_value(metadata, "R"), 0.119))
rho = env_float("RHO", something(infer_metadata_value(metadata, "rho"), 1.179))
rpm = env_float("RPM", something(infer_rpm(metadata), 5400.0))
axial_dimension = env_int("AXIAL_DIMENSION", 1)
radial_dimension = env_int("RADIAL_DIMENSION", 2)
navg_samples = env_int("NAVG_SAMPLES", 72)
p_correct_kuttacondition = env_bool("P_CORRECT_KUTTA", false)
spanwise_binning = parse_spanwise_binning(get(ENV, "SPANWISE_BINNING", "span_overlap"))
vc_target = env_float("VC_TARGET", 4.0)

# Metadata per-step uinf takes precedence inside replay; this only guards
# against runs whose metadata is missing [[step]] uinf entries.
uinf_fallback = let s = get(ENV, "UINF_FALLBACK", "4.0,0.0,0.0")
    vals = parse.(Float64, strip.(split(s, ",")))
    length(vals) == 3 || error("UINF_FALLBACK must be three comma-separated numbers")
    SVector{3, Float64}(vals...)
end

ct_vs_j_csv = get(ENV, "CT_VS_J_CSV", joinpath(save_path, "rotor_hover_ccblade_CT_vs_J.csv"))
sectional_csvs = (
    ncrit4 = get(ENV, "CCBLADE_SECTIONAL_NCRIT4",
        joinpath(save_path, "rotor_hover_ccblade_sectional_ncrit4_Vc4_J0p1867.csv")),
    ncrit9 = get(ENV, "CCBLADE_SECTIONAL_NCRIT9",
        joinpath(save_path, "rotor_hover_ccblade_sectional_ncrit9_Vc4_J0p1867.csv")),
)
run_ct_csv = get(ENV, "RUN_CT_CSV", joinpath(save_path, run_name * "_CT_vs_rev.csv"))

mkpath(out_dir)
ct_scale = rho * (rpm / 60)^2 * (2R)^4

backend = pnl.FastMultipoleBackend(;
    expansion_order=env_int("FMM_EXPANSION_ORDER", 8),
    multipole_acceptance=env_float("FMM_ACCEPTANCE", 0.4),
    leaf_size=env_int("FMM_LEAF_SIZE", 20),
)

trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i in 1:length(x)-1)

# ----------------------------------------------------------------- monitors --
monitor_store = Dict{Symbol, Any}()

monitor_factory = (systems, wakes, frames, t_range) -> begin
    body = systems[1]
    nbins = env_int("NBINS", size(first(body.shedding), 2))
    e_axial = unit_axis(axial_dimension)
    thrust_dir = -e_axial
    blade_dirs, r_min_shed, r_max_shed = infer_blade_geometry(body, frames, 1, axial_dimension)
    if length(body.shedding) == 2 ||
       any((r_max_shed .- r_min_shed) .<= sqrt(eps(Float64))) ||
       has_duplicate_blade_dirs(blade_dirs)
        blade_dirs = fallback_blade_dirs(length(body.shedding), axial_dimension, radial_dimension)
    end
    tan_dirs = [cross(e_axial, d) / norm(cross(e_axial, d)) for d in blade_dirs]
    binning = SpanwiseReplayBinning(blade_dirs, tan_dirs, e_axial, thrust_dir, R,
        fill(NaN, length(blade_dirs)), fill(NaN, length(blade_dirs)),
        nbins, 1, spanwise_binning)
    nt = length(t_range)

    pressure = pnl.PressureBernoulli(rho; unsteady=false,
        correct_kuttacondition=p_correct_kuttacondition, backend)
    force = pnl.ForceMonitor(nt, 1; i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2 * R, 1),
        correct_kuttacondition=p_correct_kuttacondition, verbose=false,
        file=false, vtk_fields=())
    snapshot = SpanwiseReplaySnapshot(:bernoulli, binning, nt)
    monitor_store[:force] = force
    monitor_store[:snapshot] = snapshot
    return (pressure, force, snapshot)
end

println("Bernoulli-only axial replay")
println("  run_name:   $(run_name)")
println("  save_path:  $(save_path)")
println("  rho/R/RPM:  $(rho) / $(R) / $(rpm)")
println("  avg window: last $(navg_samples) samples")
println("  output:     $(out_dir)")

# Single pass over every saved step: CT history and spanwise samples together.
result = pnl.replay(save_path, run_name;
    monitor_factory,
    recompute=(:auto,),
    steps=:all,
    backend,
    backend_wake=backend,
    backend_system=backend,
    Uinf=t -> uinf_fallback,
    verbose=true)

force = monitor_store[:force]
snapshot = monitor_store[:snapshot]
ct = -force.force[axial_dimension, :]
revs = result.t_range .* rpm ./ 60

# ------------------------------------------------------------ window + stats --
nsamples = length(ct)
navg = clamp(navg_samples, 1, nsamples)
window_idxs = collect((nsamples - navg + 1):nsamples)
window_stats = flatness_stats(ct[window_idxs], result.t_range[window_idxs])
@printf("\nAveraging window: samples %d:%d (rev %.3f:%.3f)\n",
    first(window_idxs), last(window_idxs),
    revs[first(window_idxs)], revs[last(window_idxs)])
@printf("  bernoulli CT mean %.8g, ptp %.4g (%.2f%%), drift %.4g (%.2f%%)\n",
    window_stats.mean, window_stats.ptp, 100 * window_stats.rel_ptp,
    window_stats.drift, 100 * window_stats.rel_drift)
passes_flatness(window_stats) ||
    @warn "Bernoulli CT is not flat over the averaging window (ptp tol 5%, drift tol 2.5%). Averaged loading may not be converged."

# --------------------------------------------------------------- CCBlade CT --
function bem_ct_from_ct_vs_j(path, ncrit, vc)
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    all(c -> c in propertynames(df), (:ncrit, :Vc, :CT)) || return nothing
    sub = df[(df.ncrit .== ncrit) .& (abs.(df.Vc .- vc) .< 1e-9), :]
    size(sub, 1) == 0 && return nothing
    return Float64(sub.CT[1])
end

function bem_ct_from_sectional(path)
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    return trapz(df.r_over_R, df.dTdr_total .* R ./ ct_scale)
end

bem_cts = Dict{String, Float64}()
for (label, ncrit) in (("ncrit4", 4), ("ncrit9", 9))
    ct_bem = bem_ct_from_ct_vs_j(ct_vs_j_csv, ncrit, vc_target)
    if isnothing(ct_bem)
        ct_bem = bem_ct_from_sectional(getproperty(sectional_csvs, Symbol(label)))
    end
    if isnothing(ct_bem)
        @warn "No CCBlade CT found for $(label) (checked $(ct_vs_j_csv) and the tagged sectional CSV)."
    else
        bem_cts[label] = ct_bem
    end
end

# ------------------------------------------------------- run CSV for context --
run_ct = isfile(run_ct_csv) ? CSV.read(run_ct_csv, DataFrame) : nothing
isnothing(run_ct) &&
    @warn "Run CT CSV not found at $(run_ct_csv); CT plot will omit the run's Laplace history and the forensic CSV will be skipped."

# ------------------------------------------------------ plot 1: CT history --
ct_plot = joinpath(out_dir, "CT_vs_revolution_bernoulli.png")
let
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(revs, ct, "-"; color="tab:blue", linewidth=1.6,
        label=@sprintf("steady Bernoulli replay (fixed), window mean CT %.5f", window_stats.mean))
    if !isnothing(run_ct)
        ax.plot(run_ct.revolution, run_ct.CT_laplace_matderiv, "-"; color="tab:green",
            linewidth=1.0, alpha=0.8, label="run PressureLaplace Du/Dt")
        ax.plot(run_ct.revolution, run_ct.CT_bernoulli, "-"; color="tab:red",
            linewidth=1.0, alpha=0.5, label="run Bernoulli (broken steady form)")
    end
    bem_styles = Dict("ncrit4" => ("black", "--"), "ncrit9" => ("tab:orange", "--"))
    for (label, ct_bem) in sort(collect(bem_cts))
        color, linestyle = bem_styles[label]
        ax.axhline(ct_bem; color, linestyle, linewidth=1.3,
            label=@sprintf("CCBlade %s CT %.5f", replace(label, "ncrit" => "ncrit="), ct_bem))
    end
    ax.set_xlabel("revolution")
    ax.set_ylabel("CT")
    ax.set_title(@sprintf("DJI 9443 axial flow, Vc=%.3g m/s (J=0.1867): CT history", vc_target))
    ax.grid(true, alpha=0.35)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ct_plot, dpi=170)
    plt.close()
end
println("Wrote $(ct_plot)")

ct_history_csv = joinpath(out_dir, "CT_vs_revolution_bernoulli.csv")
open(ct_history_csv, "w") do io
    println(io, "step,revolution,CT_bernoulli_replay_fixed")
    for (i, idx) in enumerate(result.steps)
        println(io, "$(idx),$(revs[i]),$(ct[i])")
    end
end
println("Wrote $(ct_history_csv)")

# ------------------------------------------------ plot 2: spanwise loading --
blade_count = length(snapshot.binning.blade_dirs)
stats_path = joinpath(out_dir, run_name * "_spanwise_loading_stats.csv")
stats = write_stats_csv(stats_path, Dict(:bernoulli => snapshot), (:bernoulli,),
    window_idxs, R, blade_count)
println("Wrote $(stats_path)")

span_plot = joinpath(out_dir, "spanwise_dCTdr_bernoulli.png")
let
    sub = stats[stats.source .== "bernoulli", :]
    # Average blade realizations at each radial station (each blade row already
    # carries the total-equivalent curve).
    rs = sort(unique(sub.r_over_R))
    station_mean(col) = [sum(sub[sub.r_over_R .== r, col]) / count(==(r), sub.r_over_R) for r in rs]
    med = station_mean(:median_dTdr_total_equiv) .* R ./ ct_scale
    q25 = station_mean(:q25_dTdr_total_equiv) .* R ./ ct_scale
    q75 = station_mean(:q75_dTdr_total_equiv) .* R ./ ct_scale
    panel_ct = trapz(rs, med)

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(rs, med, "o-"; color="tab:blue", markersize=3.5, linewidth=1.8,
        label=@sprintf("panel steady Bernoulli, final %d samples (CT %.5f)", navg, panel_ct))
    ax.vlines(rs, q25, q75; color="tab:blue", linewidth=2.0, alpha=0.7,
        label="panel q25–q75")
    for (label, color) in (("ncrit4", "black"), ("ncrit9", "tab:orange"))
        path = getproperty(sectional_csvs, Symbol(label))
        isfile(path) || (@warn "CCBlade sectional CSV not found at $(path)"; continue)
        bem = CSV.read(path, DataFrame)
        curve = bem.dTdr_total .* R ./ ct_scale
        ct_bem = trapz(bem.r_over_R, curve)
        ax.plot(bem.r_over_R, curve, "-"; color, linewidth=1.7,
            label=@sprintf("CCBlade %s (CT %.5f)", replace(label, "ncrit" => "ncrit="), ct_bem))
    end
    ax.set_xlabel("r/R")
    ax.set_ylabel("dCT/d(r/R)")
    ax.set_title(@sprintf("DJI 9443 axial flow, Vc=%.3g m/s (J=0.1867): spanwise loading", vc_target))
    ax.grid(true, alpha=0.35)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(span_plot, dpi=170)
    plt.close()
end
println("Wrote $(span_plot)")

# --------------------------------------------------------------- forensics --
# The run's CSV row k corresponds to saved VTU index k-1 (the run wrote
# rev = (k-1)*dt*RPM/60). Comparing the fixed replay against the saved history
# identifies which steady formulation produced the original output.
if !isnothing(run_ct)
    forensic_csv = joinpath(out_dir, "bernoulli_forensic.csv")
    run_by_step = Dict(Int(row.step) - 1 => row for row in eachrow(run_ct))
    open(forensic_csv, "w") do io
        println(io, "step,rev,CT_bernoulli_replay_fixed,CT_bernoulli_run_csv,CT_laplace_matderiv_run_csv,ratio_fixed_over_run")
        for (i, idx) in enumerate(result.steps)
            row = get(run_by_step, idx, nothing)
            isnothing(row) && continue
            ratio = abs(row.CT_bernoulli) > eps() ? ct[i] / row.CT_bernoulli : NaN
            println(io, "$(idx),$(revs[i]),$(ct[i]),$(row.CT_bernoulli),$(row.CT_laplace_matderiv),$(ratio)")
        end
    end
    println("Wrote $(forensic_csv)")

    tail = window_idxs
    run_tail = [get(run_by_step, result.steps[i], nothing) for i in tail]
    run_tail = [r for r in run_tail if !isnothing(r)]
    if !isempty(run_tail)
        run_bern = mean(r.CT_bernoulli for r in run_tail)
        run_lap = mean(r.CT_laplace_matderiv for r in run_tail)
        fixed = window_stats.mean
        println("\nForensic summary (window means):")
        @printf("  fixed steady Bernoulli replay: %.6g\n", fixed)
        @printf("  run Bernoulli (saved CSV):     %.6g\n", run_bern)
        @printf("  run Laplace Du/Dt (saved CSV): %.6g\n", run_lap)
        if abs(run_bern) < 0.5 * abs(fixed)
            println("  -> run Bernoulli is collapsed relative to the fixed replay: the run used the post-ef1fe1e inertial steady form.")
        else
            println("  -> run Bernoulli is NOT collapsed relative to the fixed replay; re-examine the assumed regression window.")
        end
    end
end

println("\nBernoulli-only axial replay complete.")
