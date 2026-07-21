## Bernoulli-only replay of the axial J=0.1867 comparison run.
##
## The original HPC run computed steady PressureBernoulli with the ef1fe1e
## inertial kinetic-energy regression (first-order blade loading cancelled).
## This script replays the saved VTK state through the corrected steady
## Bernoulli monitor only (no PressureLaplace, no CG solves) and produces:
##   1. CT vs revolution with the run's Laplace(Du/Dt) history for context and
##      CCBlade total CT per BEM_NCRITS polar set as horizontal lines;
##   2. spanwise dCT/d(r/R) averaged over the final NAVG_SAMPLES replay samples
##      (default 72 = final 2 revolutions at 36 steps/rev), overlaying the
##      CCBlade sectional curves for the same polar sets;
##   3. a saved-vs-replayed forensic CSV identifying which Bernoulli
##      formulation produced the original run output.
##
## PLOTS_ONLY=true regenerates the two plots from the CSVs a previous replay
## already wrote into OUT_DIR (no VTK state, metadata TOML, or FMM evaluation
## needed) — the mode to use on a machine that only holds the copied run CSVs.
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
plots_only = env_bool("PLOTS_ONLY", false)

# In plots-only mode the run metadata TOML is typically not available (it stays
# with the VTK state on the machine that ran the simulation).
metadata = plots_only ? nothing : TOML.parsefile(metadata_path(save_path, run_name))
meta_value(key, default) = isnothing(metadata) ?
    default : something(infer_metadata_value(metadata, key), default)

R = env_float("R", meta_value("R", 0.119))
rho = env_float("RHO", meta_value("rho", 1.179))
rpm = env_float("RPM", something(isnothing(metadata) ? nothing : infer_rpm(metadata), 5400.0))
axial_dimension = env_int("AXIAL_DIMENSION", 1)
radial_dimension = env_int("RADIAL_DIMENSION", 2)
navg_samples = env_int("NAVG_SAMPLES", 72)
p_correct_kuttacondition = env_bool("P_CORRECT_KUTTA", false)
spanwise_binning = parse_spanwise_binning(get(ENV, "SPANWISE_BINNING", "span_overlap"))
vc_target = env_float("VC_TARGET", 4.0)
op_tag = get(ENV, "OP_TAG", "Vc4_J0p1867")

# CCBlade polar sets to overlay, by `polarset` label (polar_label in
# rotor_hover_ccblade.jl: "ncrit$(n)", "inviscid", "mixed_cl-..._cd-...").
# Plotted in the listed order.
bem_sets = String.(strip.(split(get(ENV, "BEM_SETS",
    "ncrit1,ncrit2,ncrit4,ncrit9,inviscid,mixed_cl-inviscid_cd-ncrit1"), ",")))

# One fixed color per polar set, shared by both plots (identity follows the
# entity). Panel = tab:blue, run Laplace = tab:green, broken run Bernoulli =
# tab:red are reserved.
const BEM_COLORS = Dict(
    "ncrit1" => "tab:purple", "ncrit2" => "tab:cyan",
    "ncrit4" => "black", "ncrit9" => "tab:orange",
    "inviscid" => "tab:pink", "mixed_cl-inviscid_cd-ncrit1" => "tab:olive")
const BEM_FALLBACK_COLORS = ("tab:gray", "tab:brown", "gold", "navy")
bem_color(label) = get(BEM_COLORS, label,
    BEM_FALLBACK_COLORS[mod1(hash(label) % 4 + 1, length(BEM_FALLBACK_COLORS))])

# Readable legend name for a polarset label.
function bem_pretty(label)
    label == "inviscid" && return "inviscid"
    m = match(r"^mixed_cl-(.+)_cd-(.+)$", label)
    isnothing(m) || return "mixed: cl " * bem_pretty(m.captures[1]) *
                           ", cd " * bem_pretty(m.captures[2])
    return replace(label, "ncrit" => "ncrit=")
end

# Metadata per-step uinf takes precedence inside replay; this only guards
# against runs whose metadata is missing [[step]] uinf entries.
uinf_fallback = let s = get(ENV, "UINF_FALLBACK", "4.0,0.0,0.0")
    vals = parse.(Float64, strip.(split(s, ",")))
    length(vals) == 3 || error("UINF_FALLBACK must be three comma-separated numbers")
    SVector{3, Float64}(vals...)
end

mkpath(out_dir)
ct_scale = rho * (rpm / 60)^2 * (2R)^4

trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i in 1:length(x)-1)

# --------------------------------------------------- CCBlade file resolution --
# Sectional CSVs may live at the top level (possibly suffixed, e.g. `_all` from
# a local rerun) or in an `analysis/` copy of the HPC outputs.
bem_search_dirs = (save_path, joinpath(save_path, "analysis"))

function find_sectional_csv(label)
    pattern = Regex("^rotor_hover_ccblade_sectional_$(label)_$(op_tag)\\w*\\.csv\$")
    for dir in bem_search_dirs
        isdir(dir) || continue
        for f in sort(readdir(dir))
            occursin(pattern, f) && return joinpath(dir, f)
        end
    end
    return nothing
end

# Combine every CT_vs_J CSV found (base + suffixed reruns); later files win so
# a fresh local rerun overrides the older HPC rows for the same polar set.
function bem_ct_table()
    cts = Dict{String, Float64}()
    files = String[]
    haskey(ENV, "CT_VS_J_CSV") && isfile(ENV["CT_VS_J_CSV"]) && push!(files, ENV["CT_VS_J_CSV"])
    for dir in bem_search_dirs
        isdir(dir) || continue
        for f in sort(readdir(dir))
            occursin(r"^rotor_hover_ccblade_CT_vs_J\w*\.csv$", f) && push!(files, joinpath(dir, f))
        end
    end
    for path in files
        df = CSV.read(path, DataFrame)
        all(c -> c in propertynames(df), (:polarset, :Vc, :CT)) || continue
        for row in eachrow(df)
            abs(row.Vc - vc_target) < 1e-9 || continue
            cts[String(row.polarset)] = Float64(row.CT)
        end
    end
    return cts
end

function bem_ct_from_sectional(path)
    isnothing(path) && return nothing
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    return trapz(df.r_over_R, df.dTdr_total .* R ./ ct_scale)
end

sectional_paths = Dict(label => find_sectional_csv(label) for label in bem_sets)
bem_cts = let table = bem_ct_table()
    cts = Dict{String, Float64}()
    for label in bem_sets
        ct_bem = get(table, label, nothing)
        isnothing(ct_bem) && (ct_bem = bem_ct_from_sectional(sectional_paths[label]))
        if isnothing(ct_bem)
            @warn "No CCBlade CT found for polar set $(label) (no CT_vs_J row at Vc=$(vc_target) and no sectional CSV)."
        else
            cts[label] = ct_bem
        end
    end
    cts
end

run_ct_csv = let candidates = [get(ENV, "RUN_CT_CSV", ""),
        joinpath(save_path, run_name * "_CT_vs_rev.csv"),
        joinpath(save_path, "analysis", run_name * "_CT_vs_rev.csv")]
    i = findfirst(p -> !isempty(p) && isfile(p), candidates)
    isnothing(i) ? nothing : candidates[i]
end

# ------------------------------------------------- replay or CSV reload -----
ct_history_csv = joinpath(out_dir, "CT_vs_revolution_bernoulli.csv")
stats_path = joinpath(out_dir, run_name * "_spanwise_loading_stats.csv")

if plots_only

println("Bernoulli-only axial replay (plots-only mode)")
println("  run_name:   $(run_name)")
println("  save_path:  $(save_path)")
println("  rho/R/RPM:  $(rho) / $(R) / $(rpm)")
println("  BEM sets:   $(join(bem_sets, ", "))")
println("  output:     $(out_dir)")

isfile(ct_history_csv) || error("PLOTS_ONLY=true requires an existing $(ct_history_csv)")
isfile(stats_path) || error("PLOTS_ONLY=true requires an existing $(stats_path)")
hist = CSV.read(ct_history_csv, DataFrame)
replay_steps = Int.(hist.step)
ct = Float64.(hist.CT_bernoulli_replay_fixed)
revs = Float64.(hist.revolution)
t_range = revs .* 60 ./ rpm
stats = CSV.read(stats_path, DataFrame)

else

backend = pnl.FastMultipoleBackend(;
    expansion_order=env_int("FMM_EXPANSION_ORDER", 8),
    multipole_acceptance=env_float("FMM_ACCEPTANCE", 0.4),
    leaf_size=env_int("FMM_LEAF_SIZE", 20),
)

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
println("  BEM sets:   $(join(bem_sets, ", "))")
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
t_range = result.t_range
replay_steps = result.steps
revs = t_range .* rpm ./ 60

end # plots_only

# ------------------------------------------------------------ window + stats --
nsamples = length(ct)
navg = clamp(navg_samples, 1, nsamples)
window_idxs = collect((nsamples - navg + 1):nsamples)
window_stats = flatness_stats(ct[window_idxs], t_range[window_idxs])
@printf("\nAveraging window: samples %d:%d (rev %.3f:%.3f)\n",
    first(window_idxs), last(window_idxs),
    revs[first(window_idxs)], revs[last(window_idxs)])
@printf("  bernoulli CT mean %.8g, ptp %.4g (%.2f%%), drift %.4g (%.2f%%)\n",
    window_stats.mean, window_stats.ptp, 100 * window_stats.rel_ptp,
    window_stats.drift, 100 * window_stats.rel_drift)
passes_flatness(window_stats) ||
    @warn "Bernoulli CT is not flat over the averaging window (ptp tol 5%, drift tol 2.5%). Averaged loading may not be converged."

# ------------------------------------------------------- run CSV for context --
run_ct = isnothing(run_ct_csv) ? nothing : CSV.read(run_ct_csv, DataFrame)
isnothing(run_ct) &&
    @warn "Run CT CSV not found (searched SAVE_PATH and SAVE_PATH/analysis); CT plot will omit the run's Laplace history and the forensic CSV will be skipped."

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
    for label in bem_sets
        haskey(bem_cts, label) || continue
        ax.axhline(bem_cts[label]; color=bem_color(label), linestyle="--", linewidth=1.3,
            label=@sprintf("CCBlade %s CT %.5f", bem_pretty(label), bem_cts[label]))
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

if !plots_only
    open(ct_history_csv, "w") do io
        println(io, "step,revolution,CT_bernoulli_replay_fixed")
        for (i, idx) in enumerate(replay_steps)
            println(io, "$(idx),$(revs[i]),$(ct[i])")
        end
    end
    println("Wrote $(ct_history_csv)")
end

# ------------------------------------------------ plot 2: spanwise loading --
if !plots_only
    blade_count = length(snapshot.binning.blade_dirs)
    stats = write_stats_csv(stats_path, Dict(:bernoulli => snapshot), (:bernoulli,),
        window_idxs, R, blade_count)
    println("Wrote $(stats_path)")
end

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
    for label in bem_sets
        path = sectional_paths[label]
        if isnothing(path)
            @warn "CCBlade sectional CSV not found for polar set $(label); skipping its spanwise curve."
            continue
        end
        bem = CSV.read(path, DataFrame)
        curve = bem.dTdr_total .* R ./ ct_scale
        ct_bem = trapz(bem.r_over_R, curve)
        ax.plot(bem.r_over_R, curve, "-"; color=bem_color(label), linewidth=1.7,
            label=@sprintf("CCBlade %s (CT %.5f)", bem_pretty(label), ct_bem))
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
# identifies which steady formulation produced the original output. Skipped in
# plots-only mode: nothing about the replay changed.
if !plots_only && !isnothing(run_ct)
    forensic_csv = joinpath(out_dir, "bernoulli_forensic.csv")
    run_by_step = Dict(Int(row.step) - 1 => row for row in eachrow(run_ct))
    open(forensic_csv, "w") do io
        println(io, "step,rev,CT_bernoulli_replay_fixed,CT_bernoulli_run_csv,CT_laplace_matderiv_run_csv,ratio_fixed_over_run")
        for (i, idx) in enumerate(replay_steps)
            row = get(run_by_step, idx, nothing)
            isnothing(row) && continue
            ratio = abs(row.CT_bernoulli) > eps() ? ct[i] / row.CT_bernoulli : NaN
            println(io, "$(idx),$(revs[i]),$(ct[i]),$(row.CT_bernoulli),$(row.CT_laplace_matderiv),$(ratio)")
        end
    end
    println("Wrote $(forensic_csv)")

    run_tail = [get(run_by_step, replay_steps[i], nothing) for i in window_idxs]
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
