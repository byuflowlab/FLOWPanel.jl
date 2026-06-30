## Rotor hover spanwise loading replay.
##
## Replays saved VTK output from rotor_hover_pressure_comparison.jl, reruns the
## Bernoulli and PressureLaplace(lamb_vector) pressure/force/spanwise monitor
## chains, selects a converged averaging window, and writes per-blade spanwise
## loading statistics against CCBlade ncrit=4.

import FLOWPanel as pnl
import CSV
import DataFrames
using DataFrames: DataFrame
using LinearAlgebra: dot, norm
import PythonPlot as plt
using Printf: @printf, @sprintf
using Statistics: mean, median, quantile
import TOML
using FLOWPanel.FastMultipole.StaticArrays: SMatrix, SVector

import FLOWPanel: monitor_requires, monitor_provides, _run_monitor!

function env_bool(name, default)
    return parse(Bool, get(ENV, name, string(default)))
end

function env_float(name, default)
    return parse(Float64, get(ENV, name, string(default)))
end

function env_int(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function parse_steps(value::AbstractString)
    v = lowercase(strip(value))
    v == "all" && return :all
    selected = Int[]
    for part in split(v, ",")
        token = strip(part)
        isempty(token) && continue
        if occursin(":", token)
            bounds = split(token, ":")
            length(bounds) in (2, 3) || error("Invalid STEPS range: $(token)")
            if length(bounds) == 2
                append!(selected, parse(Int, bounds[1]):parse(Int, bounds[2]))
            else
                stride = parse(Int, bounds[2])
                stride != 0 || error("Invalid zero stride in STEPS range: $(token)")
                append!(selected, parse(Int, bounds[1]):stride:parse(Int, bounds[3]))
            end
        else
            push!(selected, parse(Int, token))
        end
    end
    isempty(selected) && error("STEPS did not select any saved step")
    return selected
end

function selected_replay_steps(save_path, run_name)
    haskey(ENV, "STEPS") && return parse_steps(ENV["STEPS"])
    haskey(ENV, "NSTEPS_REPLAY") || return :all
    _, idxs = pnl._read_pvd_steps(joinpath(save_path, run_name * "_body1.pvd"))
    n = min(env_int("NSTEPS_REPLAY", length(idxs)), length(idxs))
    n > 0 || error("No saved steps found for $(run_name)")
    return idxs[(end - n + 1):end]
end

function concrete_replay_steps(save_path, run_name, steps)
    steps !== :all && return steps
    _, idxs = pnl._read_pvd_steps(joinpath(save_path, run_name * "_body1.pvd"))
    isempty(idxs) && error("No saved steps found for $(run_name)")
    return idxs
end

function step_description(steps)
    steps === :all && return "all"
    isempty(steps) && return "none"
    contiguous = length(steps) == 1 || all(diff(steps) .== 1)
    contiguous && return "$(first(steps)):$(last(steps)) ($(length(steps)) samples)"
    return "$(steps) ($(length(steps)) samples)"
end

function metadata_path(save_path, run_name)
    current = joinpath(save_path, run_name * ".metadata.toml")
    isfile(current) && return current
    legacy = joinpath(save_path, run_name * ".replay.toml")
    isfile(legacy) && return legacy
    error("No replay metadata found at $(current) or $(legacy)")
end

function infer_metadata_value(metadata, key::AbstractString)
    for section in ("monitor", "body", "simulation")
        haskey(metadata, section) || continue
        items = metadata[section] isa AbstractVector ? metadata[section] : Any[metadata[section]]
        for item in items
            item isa AbstractDict || continue
            haskey(item, key) && return item[key]
            if haskey(item, "normalization") && item["normalization"] isa AbstractDict
                normdict = item["normalization"]
                haskey(normdict, key) && return normdict[key]
            end
        end
    end
    return nothing
end

function infer_rpm(metadata)
    haskey(metadata, "step") || return nothing
    steps = metadata["step"]
    isempty(steps) && return nothing
    for step in steps
        haskey(step, "frame") || continue
        for frame in step["frame"]
            haskey(frame, "omega") || continue
            omega = Float64(frame["omega"])
            omega > 0 && return omega * 60 / (2π)
        end
    end
    return nothing
end

function unit_axis(dim::Int)
    1 <= dim <= 3 || error("Axis dimension must be in 1:3; got $(dim)")
    return SVector{3, Float64}(ntuple(i -> i == dim ? 1.0 : 0.0, 3))
end

function frame_transform(frames, i_frame::Int)
    if i_frame < 0
        return SVector(0.0, 0.0, 0.0),
               SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0,
                                          0.0, 1.0, 0.0,
                                          0.0, 0.0, 1.0)
    end
    origin, R_f2g = pnl.frame_global_transform(frames, i_frame)
    return origin, transpose(R_f2g)
end

function point_frame(body, node::Integer, origin, R_g2f)
    p = SVector{3, Float64}(body.nodes[1, node], body.nodes[2, node], body.nodes[3, node])
    return R_g2f * (p - origin)
end

function infer_blade_geometry(body, frames, i_frame::Int, axial_dimension::Int)
    body isa pnl.AbstractLiftingBody || error("Replay body must be a lifting body with shedding chains.")
    isempty(body.shedding) && error("Replay body has no shedding chains; cannot infer blades.")
    origin, R_g2f = frame_transform(frames, i_frame)
    e_axial = unit_axis(axial_dimension)
    dirs = Vector{SVector{3, Float64}}(undef, length(body.shedding))
    r_min = zeros(Float64, length(body.shedding))
    r_max = zeros(Float64, length(body.shedding))
    for (i_blade, shedding) in pairs(body.shedding)
        accum = SVector(0.0, 0.0, 0.0)
        mids = SVector{3, Float64}[]
        npts = 0
        for col in eachcol(shedding)
            _, nia, nib, _, _, _ = col
            mid = 0.5 * (point_frame(body, nia, origin, R_g2f) +
                         point_frame(body, nib, origin, R_g2f))
            accum += mid - dot(mid, e_axial) * e_axial
            push!(mids, mid)
            npts += 1
        end
        npts > 0 || error("Blade $(i_blade) shedding chain has no stations.")
        d = accum / npts
        nd = norm(d)
        nd > eps(Float64) || error("Could not infer a nonzero radial direction for blade $(i_blade).")
        dirs[i_blade] = d / nd
        radial = [dot(mid, dirs[i_blade]) for mid in mids]
        r_min[i_blade] = minimum(radial)
        r_max[i_blade] = maximum(radial)
    end
    return dirs, r_min, r_max
end

function fallback_blade_dirs(n_blades::Int, axial_dimension::Int, radial_dimension::Int)
    e_axial = unit_axis(axial_dimension)
    e1 = unit_axis(radial_dimension)
    e1 = e1 - dot(e1, e_axial) * e_axial
    norm(e1) > eps(Float64) || error("RADIAL_DIMENSION must differ from AXIAL_DIMENSION.")
    e1 = e1 / norm(e1)
    e2 = SVector(e_axial[2] * e1[3] - e_axial[3] * e1[2],
                 e_axial[3] * e1[1] - e_axial[1] * e1[3],
                 e_axial[1] * e1[2] - e_axial[2] * e1[1])
    dirs = SVector{3, Float64}[]
    if n_blades == 2
        push!(dirs, e1, -e1)
    else
        for i in 0:(n_blades - 1)
            θ = 2π * i / n_blades
            push!(dirs, cos(θ) * e1 + sin(θ) * e2)
        end
    end
    return dirs
end

function panel_radial_ranges(body, frames, i_frame::Int, blade_dirs, axial_dimension::Int)
    origin, R_g2f = frame_transform(frames, i_frame)
    e_axial = unit_axis(axial_dimension)
    vals = [Float64[] for _ in blade_dirs]
    @inbounds for p in 1:body.ncells
        cp = SVector{3, Float64}(body.controlpoints[1, p], body.controlpoints[2, p], body.controlpoints[3, p])
        cp_frame = R_g2f * (cp - origin)
        blade = nearest_blade(cp_frame, blade_dirs, e_axial)
        blade == 0 && continue
        for i_node in axes(body.cells, 1)
            node = body.cells[i_node, p]
            pg = SVector{3, Float64}(body.nodes[1, node], body.nodes[2, node], body.nodes[3, node])
            pf = R_g2f * (pg - origin)
            r = dot(pf, blade_dirs[blade])
            r >= 0 && push!(vals[blade], r)
        end
    end
    all(v -> !isempty(v), vals) || error("Could not assign panels to every blade sector.")
    return [minimum(v) for v in vals], [maximum(v) for v in vals]
end

function parse_spanwise_binning(value::AbstractString)
    symbol = Symbol(lowercase(strip(value)))
    symbol in (:control_point, :span_overlap) ||
        error("SPANWISE_BINNING must be control_point or span_overlap; got $(value)")
    return symbol
end

function has_duplicate_blade_dirs(blade_dirs)
    for i in 1:length(blade_dirs)-1, j in i+1:length(blade_dirs)
        dot(blade_dirs[i], blade_dirs[j]) > 0.98 && return true
    end
    return false
end

function nearest_blade(point, blade_dirs, e_axial)
    radial = point - dot(point, e_axial) * e_axial
    nr = norm(radial)
    nr <= eps(Float64) && return 0
    rhat = radial / nr
    best_i = 0
    best_dot = -Inf
    for (i, d) in pairs(blade_dirs)
        s = dot(rhat, d)
        if s > best_dot
            best_dot = s
            best_i = i
        end
    end
    return best_i
end

mutable struct SpanwiseReplayBinning
    blade_dirs::Vector{SVector{3, Float64}}
    e_axial::SVector{3, Float64}
    thrust_dir::SVector{3, Float64}
    radius::Float64
    r_min::Vector{Float64}
    r_max::Vector{Float64}
    nbins::Int
    i_frame::Int
    binning::Symbol
end

mutable struct SpanwiseReplaySnapshot <: pnl.AbstractMonitor
    source::Symbol
    binning::SpanwiseReplayBinning
    values::Array{Float64, 3}
    counts::Array{Int, 3}
end

monitor_requires(::SpanwiseReplaySnapshot) = (:F,)
monitor_provides(::SpanwiseReplaySnapshot) = ()

function SpanwiseReplaySnapshot(source::Symbol, binning::SpanwiseReplayBinning, nt::Int)
    nblades = length(binning.blade_dirs)
    return SpanwiseReplaySnapshot(source, binning,
        fill(NaN, nt, nblades, binning.nbins),
        zeros(Int, nt, nblades, binning.nbins))
end

function _run_monitor!(m::SpanwiseReplaySnapshot, ctx::pnl.MonitorContext,
                       systems, wakes,
                       frames::AbstractVector{<:pnl.ReferenceFrame},
                       uinf, i_step::Int, dt::Real, t=nothing)
    body = systems[1]
    F = pnl.monitor_field(ctx, :F, 1)
    if any(!isfinite, m.binning.r_min) || any(!isfinite, m.binning.r_max)
        r_min, r_max = panel_radial_ranges(body, frames, m.binning.i_frame,
            m.binning.blade_dirs, findfirst(!iszero, Tuple(m.binning.e_axial)))
        m.binning.r_min .= r_min
        m.binning.r_max .= r_max
    end
    origin, R_g2f = frame_transform(frames, m.binning.i_frame)
    sums = zeros(Float64, length(m.binning.blade_dirs), m.binning.nbins)
    counts = zeros(Int, length(m.binning.blade_dirs), m.binning.nbins)

    @inbounds for p in 1:body.ncells
        cp = SVector{3, Float64}(body.controlpoints[1, p], body.controlpoints[2, p], body.controlpoints[3, p])
        cp_frame = R_g2f * (cp - origin)
        blade = nearest_blade(cp_frame, m.binning.blade_dirs, m.binning.e_axial)
        blade == 0 && continue
        r = dot(cp_frame, m.binning.blade_dirs[blade])
        r0 = m.binning.r_min[blade]
        r1 = m.binning.r_max[blade]
        width_blade = (r1 - r0) / m.binning.nbins
        width_blade > eps(Float64) || continue
        r0 <= r <= r1 || continue
        f = SVector{3, Float64}(F[1, p], F[2, p], F[3, p])
        f_frame = R_g2f * f
        thrust = dot(f_frame, m.binning.thrust_dir)
        if m.binning.binning == :span_overlap
            pmin = Inf
            pmax = -Inf
            for i_node in axes(body.cells, 1)
                node = body.cells[i_node, p]
                pg = SVector{3, Float64}(body.nodes[1, node], body.nodes[2, node], body.nodes[3, node])
                pf = R_g2f * (pg - origin)
                rp = dot(pf, m.binning.blade_dirs[blade])
                pmin = min(pmin, rp)
                pmax = max(pmax, rp)
            end
            pspan = pmax - pmin
            if pspan > sqrt(eps(Float64))
                b_first = clamp(floor(Int, (pmin - r0) / width_blade) + 1, 1, m.binning.nbins)
                b_last = clamp(floor(Int, (pmax - r0) / width_blade) + 1, 1, m.binning.nbins)
                for bin in b_first:b_last
                    bmin = r0 + (bin - 1) * width_blade
                    bmax = r0 + bin * width_blade
                    overlap = min(pmax, bmax) - max(pmin, bmin)
                    overlap > 0 || continue
                    sums[blade, bin] += (overlap / pspan) * thrust
                    counts[blade, bin] += 1
                end
                continue
            end
        end
        bin = clamp(floor(Int, (r - r0) / width_blade) + 1, 1, m.binning.nbins)
        sums[blade, bin] += thrust
        counts[blade, bin] += 1
    end

    out_i = i_step + 1
    @inbounds for blade in axes(sums, 1), bin in axes(sums, 2)
        counts[blade, bin] == 0 && continue
        width_blade = (m.binning.r_max[blade] - m.binning.r_min[blade]) / m.binning.nbins
        m.values[out_i, blade, bin] = sums[blade, bin] / width_blade
        m.counts[out_i, blade, bin] = counts[blade, bin]
    end
    return nothing
end

function make_blade_select(blade_dirs, blade, e_axial)
    return cp -> nearest_blade(SVector{3, Float64}(cp[1], cp[2], cp[3]), blade_dirs, e_axial) == blade
end

function flatness_stats(ct, times)
    μ = mean(ct)
    ptp = maximum(ct) - minimum(ct)
    rel_ptp = abs(μ) > eps(Float64) ? ptp / abs(μ) : Inf
    x = times .- mean(times)
    denom = sum(abs2, x)
    slope = denom > 0 ? sum(x .* (ct .- μ)) / denom : 0.0
    drift = slope * (maximum(times) - minimum(times))
    rel_drift = abs(μ) > eps(Float64) ? abs(drift) / abs(μ) : Inf
    return (; mean=μ, ptp, rel_ptp, drift, rel_drift)
end

passes_flatness(stats; ptp_tol=0.05, drift_tol=0.025) =
    isfinite(stats.rel_ptp) && isfinite(stats.rel_drift) &&
    stats.rel_ptp <= ptp_tol && stats.rel_drift <= drift_tol

function window_length_samples(times, rpm, nrevs)
    length(times) >= 2 || return 1
    dt = median(diff(times))
    samples = round(Int, nrevs * 60 / rpm / dt)
    return clamp(samples, 1, length(times))
end

function select_window(ct, times, rpm, requested_nrevs)
    attempts = requested_nrevs >= 2 ? (requested_nrevs, 1) : (requested_nrevs,)
    for nrevs in attempts
        n = window_length_samples(times, rpm, nrevs)
        n <= length(ct) || continue
        for stop in length(ct):-1:n
            idxs = (stop - n + 1):stop
            stats = flatness_stats(ct[idxs], times[idxs])
            passes_flatness(stats) && return collect(idxs), nrevs, stats, false
        end
    end
    n = window_length_samples(times, rpm, 1)
    idxs = collect((length(ct) - n + 1):length(ct))
    return idxs, 1, flatness_stats(ct[idxs], times[idxs]), true
end

function finite_quantiles(v)
    vals = collect(skipmissing(filter(isfinite, v)))
    isempty(vals) && return (NaN, NaN, NaN, NaN, NaN, NaN, 0)
    return (mean(vals), minimum(vals), quantile(vals, 0.25), median(vals),
            quantile(vals, 0.75), maximum(vals), length(vals))
end

function write_stats_csv(path, snapshots, sources, window_idxs, R, blade_count)
    rows = NamedTuple[]
    for source in sources
        snap = snapshots[source]
        nbins = snap.binning.nbins
        for blade in axes(snap.values, 2), bin in 1:nbins
            width = (snap.binning.r_max[blade] - snap.binning.r_min[blade]) / nbins
            r_m = snap.binning.r_min[blade] + (bin - 0.5) * width
            vals = snap.values[window_idxs, blade, bin]
            μ, lo, q25, med, q75, hi, ns = finite_quantiles(vals)
            push!(rows, (;
                source=String(source), blade, bin,
                r_over_R=r_m / R,
                r_m,
                bin_width_m=width,
                n_samples=ns,
                mean_dTdr_blade=μ,
                min_dTdr_blade=lo,
                q25_dTdr_blade=q25,
                median_dTdr_blade=med,
                q75_dTdr_blade=q75,
                max_dTdr_blade=hi,
                mean_dTdr_total_equiv=blade_count * μ,
                min_dTdr_total_equiv=blade_count * lo,
                q25_dTdr_total_equiv=blade_count * q25,
                median_dTdr_total_equiv=blade_count * med,
                q75_dTdr_total_equiv=blade_count * q75,
                max_dTdr_total_equiv=blade_count * hi,
            ))
        end
    end
    df = DataFrame(rows)
    CSV.write(path, df)
    return df
end

function load_ccblade(path)
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    (:r_over_R in propertynames(df) && :dTdr_total in propertynames(df)) ||
        error("CCBlade CSV $(path) must contain r_over_R and dTdr_total columns")
    return df
end

function plot_source(path, stats, source::Symbol, ccblade)
    sub = stats[stats.source .== String(source), :]
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    blade_colors = Dict(1 => "tab:blue", 2 => "tab:orange")
    blade_styles = Dict(1 => "-", 2 => "--")
    blade_bar_styles = Dict(1 => "solid", 2 => "dashed")
    for blade in sort(unique(sub.blade))
        bdf = sub[sub.blade .== blade, :]
        r = bdf.r_over_R
        med = bdf.median_dTdr_total_equiv
        color = get(blade_colors, blade, nothing)
        linestyle = get(blade_styles, blade, "-")
        barstyle = get(blade_bar_styles, blade, "solid")
        ax.plot(r, med; color, linestyle, marker="o", linewidth=1.8,
                markersize=3.5, label="blade $(blade) panel median")
        ax.vlines(r, bdf.min_dTdr_total_equiv, bdf.max_dTdr_total_equiv;
                  color, linestyles=barstyle, linewidth=0.6, alpha=0.45)
        ax.vlines(r, bdf.q25_dTdr_total_equiv, bdf.q75_dTdr_total_equiv;
                  color, linestyles=barstyle, linewidth=2.0, alpha=0.8)
    end
    if !isnothing(ccblade)
        ax.plot(ccblade.r_over_R, ccblade.dTdr_total, "k-"; linewidth=1.4,
                label="CCBlade ncrit=4")
    end
    ax.set_xlabel("r/R")
    ax.set_ylabel("dT/dr total-equivalent [N/m]")
    ax.set_title(source == :laplace_lamb ? "PressureLaplace Lamb spanwise loading" :
                                      "PressureBernoulli spanwise loading")
    ax.grid(true, alpha=0.35)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close()
    return path
end

function plot_ct_history(path, revs, histories, sources)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for source in sources
        ax.plot(revs, histories[source], "-o"; linewidth=1.5, markersize=3,
                label=String(source))
    end
    ax.set_xlabel("revolution")
    ax.set_ylabel("CT")
    ax.set_title("Rotor hover CT history")
    ax.grid(true, alpha=0.35)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close()
    return path
end

function parse_rev_window(value::AbstractString)
    v = lowercase(strip(value))
    v in ("", "auto") && return nothing
    parts = occursin(":", v) ? split(v, ":") : split(v, ",")
    length(parts) == 2 || error("REV_WINDOW must be 'auto' or 'rev_start:rev_stop', got $(value)")
    r0 = parse(Float64, strip(parts[1]))
    r1 = parse(Float64, strip(parts[2]))
    r1 >= r0 || error("REV_WINDOW stop must be >= start; got $(value)")
    return (r0, r1)
end

function prompt_rev_window(ct_plot_path, revs)
    if haskey(ENV, "REV_WINDOW")
        return parse_rev_window(ENV["REV_WINDOW"])
    elseif haskey(ENV, "REV_START") || haskey(ENV, "REV_STOP")
        haskey(ENV, "REV_START") && haskey(ENV, "REV_STOP") ||
            error("Set both REV_START and REV_STOP, or use REV_WINDOW=start:stop.")
        return (parse(Float64, ENV["REV_START"]), parse(Float64, ENV["REV_STOP"]))
    end

    println("\nCT history plot written to $(ct_plot_path)")
    @printf("Available saved revolution range: %.6g:%.6g\n", minimum(revs), maximum(revs))
    print("Enter revolution window as start:stop, or press Enter for automatic window selection: ")
    flush(stdout)
    try
        return parse_rev_window(readline(stdin))
    catch err
        err isa EOFError || rethrow()
        println("\nNo interactive input available; using automatic window selection.")
        return nothing
    end
end

function steps_for_rev_window(steps, revs, window)
    isnothing(window) && return steps
    r0, r1 = window
    selected = [step for (step, rev) in zip(steps, revs) if r0 <= rev <= r1]
    isempty(selected) && error("No saved steps fall inside revolution window $(r0):$(r1)")
    return selected
end

function selected_window_indices_from_revs(revs, window)
    isnothing(window) && return nothing
    r0, r1 = window
    idxs = findall(rev -> r0 <= rev <= r1, revs)
    isempty(idxs) && error("No replay samples fall inside revolution window $(r0):$(r1)")
    return idxs
end

function make_force_monitor_factory(force_store, rho, radius,
                                    p_correct_kuttacondition, backend)
    normalization = pnl.RotorNormalization(rho, 2 * radius, 1)
    return (systems, wakes, frames, t_range) -> begin
        nt = length(t_range)
        monitors = Any[]
        empty!(force_store)
        for source in (:laplace_lamb, :bernoulli)
            pressure = source == :laplace_lamb ?
                pnl.PressureLaplace(systems[1], rho; acceleration_form=:lamb_vector, verbose=false) :
                pnl.PressureBernoulli(rho; unsteady=false,
                    correct_kuttacondition=p_correct_kuttacondition, backend)
            force = pnl.ForceMonitor(nt, 1; i_frame=1, normalization,
                correct_kuttacondition=p_correct_kuttacondition, verbose=false,
                file=false, vtk_fields=())
            push!(monitors, pressure, force)
            force_store[source] = force
        end
        return Tuple(monitors)
    end
end

run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
save_path = get(ENV, "SAVE_PATH", joinpath("data", run_name))
metadata = TOML.parsefile(metadata_path(save_path, run_name))

R = env_float("R", something(infer_metadata_value(metadata, "R"), 0.119))
rho = env_float("RHO", something(infer_metadata_value(metadata, "rho"), 1.179))
rpm = env_float("RPM", something(infer_rpm(metadata), 6000.0))
axial_dimension = env_int("AXIAL_DIMENSION", 1)
radial_dimension = env_int("RADIAL_DIMENSION", 2)
nrevs_avg = env_int("NREVS_AVG", 2)
p_correct_kuttacondition = env_bool("P_CORRECT_KUTTA", false)
spanwise_binning = parse_spanwise_binning(get(ENV, "SPANWISE_BINNING", "span_overlap"))
ccblade_csv = get(ENV, "CCBLADE_CSV",
    joinpath("data", "rotor_hover_ccblade", "rotor_hover_ccblade_sectional_ncrit4.csv"))
out_dir = get(ENV, "OUT_DIR", joinpath(save_path, "spanwise_loading_replay"))
mkpath(out_dir)

plot_steps = selected_replay_steps(save_path, run_name)
normalization = pnl.RotorNormalization(rho, 2 * R, 1)
backend = pnl.FastMultipoleBackend(;
    expansion_order=env_int("FMM_EXPANSION_ORDER", 8),
    multipole_acceptance=env_float("FMM_ACCEPTANCE", 0.4),
    leaf_size=env_int("FMM_LEAF_SIZE", 20),
)

snapshots = Dict{Symbol, SpanwiseReplaySnapshot}()
forces = Dict{Symbol, pnl.ForceMonitor}()
ct_forces = Dict{Symbol, pnl.ForceMonitor}()

monitor_factory = (systems, wakes, frames, t_range) -> begin
    body = systems[1]
    nbins = env_int("NBINS", size(first(body.shedding), 2))
    e_axial = unit_axis(axial_dimension)
    thrust_dir = -e_axial
    blade_dirs, r_min_shedding, r_max_shedding = infer_blade_geometry(body, frames, 1, axial_dimension)
    if length(body.shedding) == 2 ||
       any((r_max_shedding .- r_min_shedding) .<= sqrt(eps(Float64))) ||
       has_duplicate_blade_dirs(blade_dirs)
        blade_dirs = fallback_blade_dirs(length(body.shedding), axial_dimension, radial_dimension)
    end
    r_min = fill(NaN, length(blade_dirs))
    r_max = fill(NaN, length(blade_dirs))
    binning = SpanwiseReplayBinning(blade_dirs, e_axial, thrust_dir, R,
        r_min, r_max, nbins, 1, spanwise_binning)
    nt = length(t_range)

    function make_span()
        return pnl.SpanwiseLoadingMonitor(nbins, 1;
            i_frame=1,
            span_axis=unit_axis(radial_dimension),
            components=(thrust=thrust_dir,),
            per_length=true,
            binning=spanwise_binning,
            normalization=pnl.NoSectionalNormalization(),
            file=false,
            vtk_fields=(),
            verbose=false)
    end

    monitors = Any[]
    for source in (:laplace_lamb, :bernoulli)
        pressure = source == :laplace_lamb ?
            pnl.PressureLaplace(systems[1], rho; acceleration_form=:lamb_vector, verbose=false) :
            pnl.PressureBernoulli(rho; unsteady=false,
                correct_kuttacondition=p_correct_kuttacondition, backend)
        force = pnl.ForceMonitor(nt, 1; i_frame=1, normalization,
            correct_kuttacondition=p_correct_kuttacondition, verbose=false,
            file=false, vtk_fields=())
        snapshot = SpanwiseReplaySnapshot(source, binning, nt)
        push!(monitors, pressure, force)
        push!(monitors, make_span())
        push!(monitors, snapshot)
        forces[source] = force
        snapshots[source] = snapshot
    end
    return Tuple(monitors)
end

println("Rotor hover spanwise loading replay")
println("  run_name:  $(run_name)")
println("  save_path: $(save_path)")
println("  CT steps:  $(step_description(plot_steps))")
println("  rho/R/RPM: $(rho) / $(R) / $(rpm)")
println("  output:    $(out_dir)")

sources = (:laplace_lamb, :bernoulli)
ct_result = pnl.replay(save_path, run_name;
    monitor_factory=make_force_monitor_factory(ct_forces, rho, R,
        p_correct_kuttacondition, backend),
    recompute=(:auto,),
    steps=plot_steps,
    backend,
    backend_wake=backend,
    backend_system=backend,
    verbose=true)

ct_steps = concrete_replay_steps(save_path, run_name, plot_steps)
ct_revs = ct_result.t_range .* rpm ./ 60
ct_histories = Dict(source => -ct_forces[source].force[axial_dimension, :]
                    for source in sources)
ct_plot = joinpath(out_dir, run_name * "_CT_vs_revolution.png")
plot_ct_history(ct_plot, ct_revs, ct_histories, sources)
rev_window = prompt_rev_window(ct_plot, ct_revs)
steps = steps_for_rev_window(ct_steps, ct_revs, rev_window)
println("  replay steps: $(step_description(steps))")

result = pnl.replay(save_path, run_name;
    monitor_factory,
    recompute=(:auto,),
    steps,
    backend,
    backend_wake=backend,
    backend_system=backend,
    verbose=true)

ct_laplace = -forces[:laplace_lamb].force[axial_dimension, :]
result_revs = result.t_range .* rpm ./ 60
manual_window_idxs = selected_window_indices_from_revs(result_revs, rev_window)
if isnothing(manual_window_idxs)
    window_idxs, selected_nrevs, primary_stats, forced_latest = select_window(
        ct_laplace, result.t_range, rpm, nrevs_avg)
else
    window_idxs = manual_window_idxs
    selected_nrevs = max(1, round(Int, ceil(maximum(result_revs[window_idxs]) -
                                             minimum(result_revs[window_idxs]))))
    primary_stats = flatness_stats(ct_laplace[window_idxs], result.t_range[window_idxs])
    forced_latest = false
end

first_i = first(window_idxs)
last_i = last(window_idxs)
println("\nSelected averaging window:")
@printf("  samples: %d:%d of %d, saved steps %d:%d, rev %.3f:%.3f (%d rev target)\n",
    first_i, last_i, length(result.t_range),
    result.steps[first_i], result.steps[last_i],
    result.t_range[first_i] * rpm / 60,
    result.t_range[last_i] * rpm / 60,
    selected_nrevs)
@printf("  laplace_lamb CT mean %.8g, ptp %.4g (%.2f%%), drift %.4g (%.2f%%)\n",
    primary_stats.mean, primary_stats.ptp, 100 * primary_stats.rel_ptp,
    primary_stats.drift, 100 * primary_stats.rel_drift)
forced_latest && @warn "No converged 2-rev or 1-rev CT window found; using latest 1-rev window."

for source in sources
    ct = -forces[source].force[axial_dimension, window_idxs]
    stats = flatness_stats(ct, result.t_range[window_idxs])
    @printf("  %-13s CT mean %.8g, ptp %.4g (%.2f%%), drift %.4g (%.2f%%)\n",
        String(source), stats.mean, stats.ptp, 100 * stats.rel_ptp,
        stats.drift, 100 * stats.rel_drift)
    passes_flatness(stats) || @warn "$(source) CT does not meet flatness checks over the selected window."
end

blade_count = length(snapshots[:laplace_lamb].binning.blade_dirs)
stats_path = joinpath(out_dir, run_name * "_spanwise_loading_stats.csv")
stats = write_stats_csv(stats_path, snapshots, sources, window_idxs, R, blade_count)
println("\nWrote $(stats_path)")

ccblade = load_ccblade(ccblade_csv)
isnothing(ccblade) && @warn "CCBlade ncrit=4 CSV not found at $(ccblade_csv); plots will omit BEM overlay."

laplace_plot = joinpath(out_dir, run_name * "_spanwise_loading_laplace_lamb.png")
bernoulli_plot = joinpath(out_dir, run_name * "_spanwise_loading_bernoulli.png")
plot_source(laplace_plot, stats, :laplace_lamb, ccblade)
plot_source(bernoulli_plot, stats, :bernoulli, ccblade)
println("Wrote $(laplace_plot)")
println("Wrote $(bernoulli_plot)")

finite_cols = [:mean_dTdr_blade, :q25_dTdr_blade, :median_dTdr_blade, :q75_dTdr_blade]
all(all(isfinite, skipmissing(stats[!, c])) for c in finite_cols) ||
    @warn "Some spanwise statistic columns contain non-finite values; inspect bins with n_samples=0."
