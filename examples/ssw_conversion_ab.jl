#=##############################################################################
# DESCRIPTION
#   BRAINSTORM 016 Phase 4 §3 case 1: suddenly-started AR=6 wing A/B between
#   the legacy edge-jump panel-to-particle conversion and the opt-in smooth
#   surface-vorticity conversion. Within each pair ONLY the conversion
#   strategy differs: identical geometry, wake, particle sigma target, solver,
#   timestep, maintenance, and monitors. The smooth sigma default resolves to
#   the legacy unsteady-row sigma (pps_overlap*dt_star/pps_n), pinning the
#   effective particle resolution across the pair.
#
#   Pre-registered interpretation (committed before the first run):
#   - Both arms of every pair must settle (ssw_settled_stats gates). A smooth
#     arm that destabilizes or fails to settle where legacy settles is a
#     NEGATIVE result for the smooth strategy at these settings, recorded as
#     such.
#   - The primary observable is tail delta-CL (smooth - legacy) at matched
#     settings, with its trend over the dt* ladder; secondary observables are
#     the settled Gamma_TE and spanwise-loading profiles, particle count, and
#     the conversion ledger residual (must stay at round-off).
#   - The sigma/overlap arm exposes resolution masquerading as smoothing: if
#     the smooth-vs-legacy delta moves materially with conversion sigma or
#     overlap at fixed dt*, the comparison is resolution-confounded and is
#     reported per-sigma rather than as one number.
#   - The attribution arm (:upstream vs :split) is a wing-side datapoint for
#     the rotor T4 null only; interpretation is deferred to the rotor session.
#   - A measurable delta, non-convergence, or sensitivity merely transferred
#     to another parameter are all valid conclusions (Phase 4 §3).
#
# RUN
#   julia --project -t 6 examples/ssw_conversion_ab.jl
#
# ENVIRONMENT
#   SSWCAB_DTS       dt_star ladder (default: 0.25,0.125,0.0625)
#   SSWCAB_TEND      t_end_star (default: 20)
#   SSWCAB_NSPAN     spanwise strips (default: 24)
#   SSWCAB_ROWS      panel_rows sheet buffer (default: 8)
#   SSWCAB_ETA       Das eta (default: 1.0 — on the log-law plateau)
#   SSWCAB_ARMS      comma list of primary,sigma,attribution (default: all)
#   SSWCAB_OUTPUT    output root (default: data/ssw_conversion_ab)
#   SSWCAB_RESUME    reuse completed histories (default: true)
#   SSWCAB_SMOKE     single coarse pair only (default: false)
=###############################################################################

if !isdefined(@__MODULE__, :SSWConfig)
    include(joinpath(@__DIR__, "suddenly_started_wing.jl"))
end

import Printf: @printf

_cab_list(name, default) = parse.(Float64, strip.(split(get(ENV, name, default), ',')))

function _cab_means(result, name)
    direct = name == :gamma ? (:gamma_te in keys(result) ? vec(sum(
        result.gamma_te[:, result.tail_indices]; dims=2) ./ length(result.tail_indices)) :
        nothing) : (:spanwise_loading in keys(result) ? vec(sum(
        result.spanwise_loading[:, result.tail_indices]; dims=2) ./
        length(result.tail_indices)) : nothing)
    !isnothing(direct) && return direct
    artifact = name == :gamma ? result.gamma_artifact : result.loading_artifact
    isnothing(artifact) && throw(ArgumentError("cached case $(result.tag) lacks $(name) artifact"))
    return artifact.mean
end

"Conversion-ledger totals: from the in-memory result for fresh runs, from the
sibling conversion_diagnostics.csv for resumed ones, NaN-filled for legacy."
function _cab_conversion_stats(result)
    :conversion_stats in keys(result) && return result.conversion_stats
    path = joinpath(result.path, "conversion_diagnostics.csv")
    isfile(path) || return (; conversions=0, particles=0, residual_norm=NaN)
    fields = split(readlines(path)[2], ',')
    return (; conversions=parse(Int, fields[1]),
        particles=parse(Int, fields[2]),
        residual_norm=parse(Float64, fields[9]))
end

function _cab_row(label, smooth, legacy)
    gamma = _cab_means(smooth, :gamma)
    gamma_ref = _cab_means(legacy, :gamma)
    loading = _cab_means(smooth, :loading)
    loading_ref = _cab_means(legacy, :loading)
    stats = _cab_conversion_stats(smooth)
    return (; label,
        dt_star=smooth.config.dt_star,
        conversion_sigma_over_c=_ssw_conversion_sigma_over_c(smooth.config),
        conversion_overlap=smooth.config.conversion_overlap,
        attribution=smooth.config.conversion_attribution,
        legacy_tail_CL=legacy.tail_CL.mean,
        smooth_tail_CL=smooth.tail_CL.mean,
        delta_CL_percent=100 * (smooth.tail_CL.mean - legacy.tail_CL.mean) /
            abs(legacy.tail_CL.mean),
        gamma_error_percent=100 * maximum(abs.(gamma .- gamma_ref)) /
            max(maximum(abs, gamma_ref), eps(Float64)),
        loading_error_percent=100 * maximum(abs.(loading .- loading_ref)) /
            max(maximum(abs, loading_ref), eps(Float64)),
        legacy_settled=legacy.tail_CL.settled && legacy.tail_gamma.settled,
        smooth_settled=smooth.tail_CL.settled && smooth.tail_gamma.settled,
        legacy_CL_ripple_percent=100legacy.tail_CL.ripple,
        smooth_CL_ripple_percent=100smooth.tail_CL.ripple,
        legacy_particles=legacy.n_particles,
        smooth_particles=smooth.n_particles,
        conversions=stats.conversions,
        conversion_residual=stats.residual_norm,
        legacy_tag=legacy.tag, smooth_tag=smooth.tag)
end

function _write_cab_summary(rows, path)
    open(path, "w") do io
        names = keys(first(rows))
        println(io, join(names, ","))
        for row in rows
            println(io, join((getproperty(row, n) for n in names), ","))
        end
    end
    return path
end

function run_ssw_conversion_ab()
    base = SSWConfig(
        AR=6.0,
        n_span=parse(Int, get(ENV, "SSWCAB_NSPAN", "24")),
        n_airfoil=21,
        eta=parse(Float64, get(ENV, "SSWCAB_ETA", "1.0")),
        dt_star=0.125,
        t_end_star=parse(Float64, get(ENV, "SSWCAB_TEND", "20.0")),
        wake_model=:particle,
        panel_rows=parse(Int, get(ENV, "SSWCAB_ROWS", "8")),
        # FMM by default (the SSWConfig default, per Ryan 2026-08-03); the
        # backend is shared by both arms of every pair, so the A/B stays
        # backend-consistent either way. SSWCAB_BACKEND=direct for cross-checks.
        backend_kind=Symbol(get(ENV, "SSWCAB_BACKEND", "fmm")),
        output_root=get(ENV, "SSWCAB_OUTPUT",
            joinpath("data", "ssw_conversion_ab")),
        save_vtk=false,
        verbose=_envbool("SSW_VERBOSE", false),
    )
    dts = _cab_list("SSWCAB_DTS", "0.25,0.125,0.0625")
    arms = Symbol.(strip.(split(get(ENV, "SSWCAB_ARMS", "primary,sigma,attribution"), ',')))
    resume = _envbool("SSWCAB_RESUME", true)
    smoke = _envbool("SSWCAB_SMOKE", false)
    mkpath(base.output_root)

    run_case(c) = begin
        while true
            case_dir = joinpath(c.output_root, _ssw_case_tag(c))
            complete = all(isfile, (joinpath(case_dir, "history.csv"),
                joinpath(case_dir, "gamma_te.csv"),
                joinpath(case_dir, "spanwise_loading.csv"),
                joinpath(case_dir, "settledness.csv")))
            result = if resume && complete
                load_ssw_result(c)
            else
                relaxation, recorder = ssw_recording_relaxation()
                run_suddenly_started_wing(c; relaxation,
                    relaxation_recorder=recorder)
            end
            result.tail_CL.settled && result.tail_gamma.settled && return result
            c.t_end_star >= 40 && return result
            c = ssw_with(c; t_end_star=min(40.0, c.t_end_star + 5.0))
            @printf("  extending unsettled case to t*=%.1f\n", c.t_end_star)
        end
    end

    legacy_at = Dict{Float64,Any}()
    function legacy_case(dt)
        haskey(legacy_at, dt) && return legacy_at[dt]
        @printf("legacy  dt*=%g\n", dt)
        legacy_at[dt] = run_case(ssw_with(base; dt_star=dt))
        return legacy_at[dt]
    end
    smooth_case(dt; kwargs...) = begin
        @printf("smooth  dt*=%g %s\n", dt, isempty(kwargs) ? "" : string((; kwargs...)))
        run_case(ssw_with(base; dt_star=dt, conversion=:smooth, kwargs...))
    end

    rows = NamedTuple[]
    if smoke
        dt = maximum(dts)
        push!(rows, _cab_row("smoke", smooth_case(dt), legacy_case(dt)))
    else
        if :primary in arms
            for dt in dts
                push!(rows, _cab_row("primary", smooth_case(dt), legacy_case(dt)))
            end
        end
        if :sigma in arms
            dt = 0.125
            base_sigma = _ssw_conversion_sigma_over_c(ssw_with(base; dt_star=dt))
            for factor in (0.5, 2.0)
                push!(rows, _cab_row("sigma_x$(factor)",
                    smooth_case(dt; conversion_sigma_over_c=factor * base_sigma),
                    legacy_case(dt)))
            end
            for overlap in (1.0, 1.6)
                push!(rows, _cab_row("overlap_$(overlap)",
                    smooth_case(dt; conversion_overlap=overlap),
                    legacy_case(dt)))
            end
        end
        if :attribution in arms
            dt = 0.125
            push!(rows, _cab_row("attribution_split",
                smooth_case(dt; conversion_attribution=:split),
                legacy_case(dt)))
        end
    end

    csv_path = joinpath(base.output_root, "conversion_ab_summary.csv")
    _write_cab_summary(rows, csv_path)

    plot_path = nothing
    if SSW_PLOTTING
        primary = [r for r in rows if r.label == "primary"]
        if !isempty(primary)
            fig = plt.figure("ssw_conversion_ab", figsize=(6.5, 4.2))
            fig.clear()
            ax = fig.add_subplot(1, 1, 1)
            pts = sort([(r.dt_star, r.delta_CL_percent) for r in primary])
            ax.semilogx(first.(pts), last.(pts), "o-"; label="smooth - legacy")
            ax.axhline(0.0; color="k", lw=0.8)
            ax.set(xlabel="dt*", ylabel="tail ΔCL (%)",
                title="AR=$(Int(base.AR)) wing: smooth vs legacy conversion, " *
                    "rows=$(base.panel_rows), eta=$(base.eta)")
            ax.grid(true; which="both", alpha=0.3)
            ax.legend(fontsize=8)
            fig.tight_layout()
            plot_path = joinpath(base.output_root, "conversion_ab.png")
            fig.savefig(plot_path; dpi=180, bbox_inches="tight")
        end
    end

    println("\nconversion A/B summary:")
    @printf("%18s %8s %10s %10s %10s %8s %8s %12s\n",
        "label", "dt*", "dCL(%)", "dGamma(%)", "dLoad(%)", "L.set", "S.set",
        "conv.resid")
    for r in rows
        @printf("%18s %8.4g %+10.4f %10.4f %10.4f %8s %8s %12.3e\n",
            r.label, r.dt_star, r.delta_CL_percent, r.gamma_error_percent,
            r.loading_error_percent, string(r.legacy_settled),
            string(r.smooth_settled), r.conversion_residual)
    end
    return (; rows, csv_path, plot_path)
end

"""
    run_ssw_conversion_probes()

Phase 4 §3 "induced velocity at fixed probes": run one matched legacy/smooth
pair for a short horizon IN MEMORY (no VTK) and compare the induced velocity
of the full hybrid wake (panel rows + retained filament + particles) at fixed
probe points. Production `simulate!` sets `overflowed[]`, so the retained
filament genuinely radiates here (unlike static fixtures — review finding R4).
Writes `probe_comparison.csv` next to the A/B summary.
"""
function run_ssw_conversion_probes(; dt_star=0.25, t_end_star=5.0)
    base = SSWConfig(
        AR=6.0,
        n_span=parse(Int, get(ENV, "SSWCAB_NSPAN", "24")),
        n_airfoil=21,
        eta=parse(Float64, get(ENV, "SSWCAB_ETA", "1.0")),
        dt_star=dt_star,
        t_end_star=t_end_star,
        wake_model=:particle,
        panel_rows=parse(Int, get(ENV, "SSWCAB_ROWS", "8")),
        backend_kind=Symbol(get(ENV, "SSWCAB_BACKEND", "fmm")),
        output_root=get(ENV, "SSWCAB_OUTPUT",
            joinpath("data", "ssw_conversion_ab")),
        save_vtk=false,
        verbose=false,
    )
    mkpath(base.output_root)

    function final_wake(config)
        sim = prepare_suddenly_started_wing(config)
        pnl.simulate!((sim.wing,), (sim.wake,), sim.frames, sim.maneuver!,
            sim.Uinf, sim.t_range;
            body_solvers=(sim.solver,), backend=sim.backend,
            monitors=sim.monitors, path=nothing, name="probe",
            set_Das_eta_freestream=NaN, grad_mu_options=SSW_GRAD_MU_OPTIONS,
            verbose=false)
        return sim.wake
    end

    # Fixed probes: downstream of the trailing edge on/off the wake centerline
    # and outboard of the tip, in chords (span is y, freestream ~x).
    semib = base.AR * base.c / 2
    probes_x = SVector{3,Float64}[]
    for x in (1.5, 3.0, 5.0), z in (0.25, 1.0), y in (0.0, 0.6 * semib, 1.1 * semib)
        push!(probes_x, SVector(x * base.c, y, z * base.c))
    end

    function probe_velocity(wake)
        probes = pnl.FastMultipole.ProbeSystem(length(probes_x), Float64)
        for i in eachindex(probes_x)
            probes.position[i] = probes_x[i]
            probes.scalar_potential[i] = 0.0
            probes.gradient[i] = zero(SVector{3,Float64})
            probes.hessian[i] = zero(pnl.FastMultipole.SMatrix{3,3,Float64,9})
        end
        pnl.influence!((probes,),
            (pnl.get_sources(wake.panel_wake)..., wake.pfield),
            pnl.DirectBackend(); scalar_potential=false, gradient=true,
            hessian=(false,))
        return copy(probes.gradient)
    end

    println("probe pair: legacy")
    u_legacy = probe_velocity(final_wake(base))
    println("probe pair: smooth")
    u_smooth = probe_velocity(final_wake(ssw_with(base; conversion=:smooth)))

    ref = sqrt(sum(LA.norm.(u_legacy) .^ 2) / length(u_legacy))
    csv_path = joinpath(base.output_root, "probe_comparison.csv")
    open(csv_path, "w") do io
        println(io, "x_over_c,y_over_c,z_over_c,u_legacy,u_smooth," *
            "delta_norm,delta_over_rms_legacy")
        for (i, p) in enumerate(probes_x)
            d = LA.norm(u_smooth[i] - u_legacy[i])
            @printf(io, "%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e\n",
                p[1] / base.c, p[2] / base.c, p[3] / base.c,
                LA.norm(u_legacy[i]), LA.norm(u_smooth[i]), d, d / ref)
        end
    end
    deltas = [LA.norm(u_smooth[i] - u_legacy[i]) / ref
              for i in eachindex(probes_x)]
    @printf("probe deltas (relative to legacy RMS): max=%.4e mean=%.4e\n",
        maximum(deltas), sum(deltas) / length(deltas))
    return (; csv_path, max_delta=maximum(deltas),
        mean_delta=sum(deltas) / length(deltas))
end

if abspath(PROGRAM_FILE) == @__FILE__
    if _envbool("SSWCAB_PROBES", false)
        run_ssw_conversion_probes()
    else
        run_ssw_conversion_ab()
    end
end
