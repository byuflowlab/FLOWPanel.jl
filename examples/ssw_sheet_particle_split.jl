#=##############################################################################
# DESCRIPTION
#   BRAINSTORM 014 follow-up: does circulation/lift on a SIMPLE RECTANGULAR
#   WING depend on how much of the wake is panel sheet vs vortex particles?
#
#   The rotor shows CT rising monotonically with sheet extent (0.034 -> 0.110
#   across discretizations); existing wing evidence (two single points,
#   code_audit Tasks 1/2b) says particles ~ sheets to <0.4%. This script runs
#   the missing systematic sweep: a `panel_rows` ladder on the suddenly-started
#   AR=6 wing with a PanelParticleWake, against the pure-PanelWake control.
#
#   Pre-registered interpretation:
#   - FLAT (±0.5%) at wing sigma  => particle representation fine on wings;
#     rotor effect is rotor-specific. The overlap wave then tests whether
#     sigma/chord is the bridge variable.
#   - RISING toward the panel control => generic particle under-induction;
#     slope vs buffer length (and sigma) routes the mainline to 008/011/012.
#
# RUN
#   julia --project -t 6 examples/ssw_sheet_particle_split.jl
#
# ENVIRONMENT
#   SSWSPS_METHOD    sigma_pps, sigma_overlap, or all (default: all)
#   SSWSPS_SIGMAS    sigma/c list (default: 0.08,0.15,0.3,0.6,1.2)
#   SSWSPS_PPS       primary SigmaPPS value (default: 2)
#   SSWSPS_BUFFERS   sheet-buffer lengths in chords (default: 0.25,1,4)
#   SSWSPS_OVERLAP   SigmaOverlap target overlap (default: 1.3)
#   SSWSPS_ETA       Das eta (default: 1.0 — on the log-law plateau)
#   SSWSPS_DT        dt_star (default: 0.125)
#   SSWSPS_TEND     t_end_star (default: 20)
#   SSWSPS_NSPAN     spanwise strips (default: 24)
#   SSWSPS_OUTPUT    output root (default: data/ssw_sheet_particle_split)
#   SSWSPS_RESUME    reuse completed histories (default: true)
=###############################################################################

if !isdefined(@__MODULE__, :SSWConfig)
    include(joinpath(@__DIR__, "suddenly_started_wing.jl"))
end

import Printf: @printf

_sps_list(name, default) = parse.(Float64, strip.(split(get(ENV, name, default), ',')))

mutable struct SSWMatchedGeometry <: pnl.AbstractParticleFunctionalPolicy
    uinf::SVector{3,Float64}
    targets::Vector{SVector{3,Float64}}
    sum_correction::Float64
    max_correction::Float64
    samples::Int
end

SSWMatchedGeometry(uinf) = SSWMatchedGeometry(SVector{3,Float64}(uinf),
    SVector{3,Float64}[], 0.0, 0.0, 0)

"""
After the normal relaxed Euler update, replace each pre-existing particle's
position by its freestream-convected position. Newly appearing particles are
registered on their first update by exactly undoing that Euler displacement;
their strengths, smoothing radii, and all subsequent strength dynamics remain
untouched.
"""
function pnl.apply_particle_policy!(policy::SSWMatchedGeometry, pfield,
        ctx::pnl.ParticleMaintenanceContext)
    np = pfield.np
    for i in eachindex(policy.targets)
        policy.targets[i] += ctx.dt * policy.uinf
    end
    for i in length(policy.targets) + 1:np
        x_after = SVector{3,Float64}(pnl.FLOWVPM.get_X(pfield, i))
        velocity = SVector{3,Float64}(
            view(pfield.particles, pnl.FLOWVPM.U_INDEX, i))
        push!(policy.targets,
            x_after - ctx.dt * velocity + ctx.dt * policy.uinf)
    end
    for i in 1:np
        x = pnl.FLOWVPM.get_X(pfield, i)
        correction = LA.norm(SVector{3}(x) - policy.targets[i])
        policy.sum_correction += correction
        policy.max_correction = max(policy.max_correction, correction)
        policy.samples += 1
        x .= policy.targets[i]
    end
    return nothing
end

ssw_matched_geometry_stats(policy::SSWMatchedGeometry) =
    (; mean_position_correction=policy.samples == 0 ? 0.0 :
            policy.sum_correction / policy.samples,
        max_position_correction=policy.max_correction,
        samples=policy.samples)

function _sps_artifact_means(result, name)
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

function _sps_wake_extent(result)
    return :wake_extent_over_c in keys(result) ?
        result.wake_extent_over_c : NaN
end

function _sps_case_row(result, control, method, sigma, pps, overlap, buffer)
    gamma = _sps_artifact_means(result, :gamma)
    gamma_control = _sps_artifact_means(control, :gamma)
    loading = _sps_artifact_means(result, :loading)
    loading_control = _sps_artifact_means(control, :loading)
    dcl_percent = 100 * (result.tail_CL.mean - control.tail_CL.mean) /
        abs(control.tail_CL.mean)
    gamma_error_percent = 100 * maximum(abs.(gamma .- gamma_control)) /
        max(maximum(abs, gamma_control), eps(Float64))
    loading_error_percent = 100 * maximum(abs.(loading .- loading_control)) /
        max(maximum(abs, loading_control), eps(Float64))
    trailing_length = result.config.dt_star
    unsteady_length = result.config.AR / result.config.n_span
    if method == :sigma_pps
        h_trailing_over_sigma = trailing_length / (pps * sigma)
        h_unsteady_over_sigma = unsteady_length / (pps * sigma)
    else
        h_trailing_over_sigma = trailing_length /
            (max(1, ceil(Int, overlap * trailing_length / sigma)) * sigma)
        h_unsteady_over_sigma = unsteady_length /
            (max(1, ceil(Int, overlap * unsteady_length / sigma)) * sigma)
    end
    eligible = method in (:sigma_overlap, :sigma_overlap_growing)
    settled = result.tail_CL.settled && result.tail_gamma.settled
    admissible = eligible && settled &&
        abs(dcl_percent) <= 0.5 && gamma_error_percent <= 1.0
    relaxation = :relaxation_stats in keys(result) ? result.relaxation_stats :
        (; mean_relative_change=NaN, max_relative_change=NaN, samples=0)
    return (; method, sigma_over_c=sigma, pps, overlap, buffer_over_c=buffer,
        panel_rows=result.config.panel_rows,
        tail_CL=result.tail_CL.mean, delta_CL_percent=dcl_percent,
        CL_drift_percent=100result.tail_CL.drift,
        CL_ripple_percent=100result.tail_CL.ripple,
        gamma_drift_percent=100result.tail_gamma.drift,
        gamma_ripple_percent=100result.tail_gamma.ripple,
        gamma_error_percent, loading_error_percent,
        settled, eligible, admissible,
        n_particles=result.n_particles,
        wake_extent_over_c=_sps_wake_extent(result),
        wake_extent_estimated=:wake_extent_estimated in keys(result) ?
            result.wake_extent_estimated : false,
        sigma_over_d=sigma / buffer, d_over_sigma=buffer / sigma,
        h_trailing_over_sigma, h_unsteady_over_sigma,
        relaxation_mean_percent=100relaxation.mean_relative_change,
        relaxation_max_percent=100relaxation.max_relative_change,
        tag=result.tag)
end

function _write_sps_summary(rows, path)
    open(path, "w") do io
        names = keys(first(rows))
        println(io, join(names, ","))
        for row in rows
            println(io, join((getproperty(row, n) for n in names), ","))
        end
    end
    return path
end

function run_ssw_sheet_particle_split()
    base = SSWConfig(
        AR=6.0,
        n_span=parse(Int, get(ENV, "SSWSPS_NSPAN", "24")),
        n_airfoil=21,
        eta=parse(Float64, get(ENV, "SSWSPS_ETA", "1.0")),
        wake_core_over_c=parse(Float64, get(ENV, "SSWSPS_CORE", "0.001")),
        dt_star=parse(Float64, get(ENV, "SSWSPS_DT", "0.125")),
        t_end_star=parse(Float64, get(ENV, "SSWSPS_TEND", "20.0")),
        backend_kind=:direct,
        output_root=get(ENV, "SSWSPS_OUTPUT",
            joinpath("data", "ssw_sheet_particle_split")),
        save_vtk=false,
        verbose=_envbool("SSW_VERBOSE", false),
    )
    sigmas = _sps_list("SSWSPS_SIGMAS", "0.08,0.15,0.3,0.6,1.2")
    buffers = _sps_list("SSWSPS_BUFFERS", "0.25,1,4")
    primary_pps = parse(Int, get(ENV, "SSWSPS_PPS", "2"))
    overlap = parse(Float64, get(ENV, "SSWSPS_OVERLAP", "1.3"))
    method_mode = Symbol(lowercase(get(ENV, "SSWSPS_METHOD", "all")))
    method_mode in (:all, :sigma_pps, :sigma_overlap) ||
        throw(ArgumentError("SSWSPS_METHOD must be all, sigma_pps, or sigma_overlap"))
    prediction_path = get(ENV, "SSWSPS_PREDICTION",
        joinpath("data", "ssw_representation_probe", "phase_b_prediction.csv"))
    !_envbool("SSWSPS_SMOKE", false) && !isfile(prediction_path) &&
        throw(ArgumentError("Phase B requires immutable Phase A prediction at " *
            prediction_path))
    resume = _envbool("SSWSPS_RESUME", true)
    mkpath(base.output_root)

    run_case(c; particle_maintenance=pnl.ParticleMaintenance()) = begin
        while true
            case_dir = joinpath(c.output_root, _ssw_case_tag(c))
            complete = all(isfile, (joinpath(case_dir, "history.csv"),
                joinpath(case_dir, "gamma_te.csv"),
                joinpath(case_dir, "spanwise_loading.csv"),
                joinpath(case_dir, "settledness.csv")))
            if resume && complete
                result = load_ssw_result(c)
            elseif c.wake_model == :particle
                relaxation, recorder = ssw_recording_relaxation()
                result = run_suddenly_started_wing(c; relaxation,
                    relaxation_recorder=recorder, particle_maintenance)
            else
                result = run_suddenly_started_wing(c)
            end
            result.tail_CL.settled && result.tail_gamma.settled && return result
            c.t_end_star >= 40 && return result
            c = ssw_with(c; t_end_star=min(40.0, c.t_end_star + 5.0))
            @printf("  extending unsettled case to t*=%.1f\n", c.t_end_star)
        end
    end

    # pure-panel control (sheet spans every step)
    println("control: pure PanelWake")
    control = run_case(ssw_with(base; wake_model=:panel))

    cases = NamedTuple[]
    smoke = _envbool("SSWSPS_SMOKE", false)
    if smoke
        push!(cases, (; method=:sigma_pps, sigma=first(sigmas),
            pps=primary_pps, overlap, buffer=first(buffers)))
    end
    if !smoke && method_mode in (:all, :sigma_pps)
        for sigma in sigmas, buffer in buffers
            push!(cases, (; method=:sigma_pps, sigma, pps=primary_pps, overlap,
                buffer))
        end
        for sigma in (0.08, 1.2), pps in (1, 4), buffer in (0.25, 4.0)
            any(c -> c.method == :sigma_pps && c.sigma == sigma &&
                c.pps == pps && c.buffer == buffer, cases) ||
                push!(cases, (; method=:sigma_pps, sigma, pps, overlap, buffer))
        end
    end
    if !smoke && method_mode in (:all, :sigma_overlap)
        for sigma in sigmas, buffer in buffers
            push!(cases, (; method=:sigma_overlap, sigma, pps=primary_pps,
                overlap, buffer))
        end
    end

    results = NamedTuple[]
    function execute_case(case)
        rows = max(1, round(Int, case.buffer / base.dt_star))
        @printf("case method=%s sigma/c=%g pps=%d overlap=%g buffer/c=%g rows=%d\n",
            string(case.method), case.sigma, case.pps, case.overlap, case.buffer, rows)
        config = ssw_with(base; wake_model=:particle, panel_rows=rows,
            shed_method=case.method == :sigma_overlap_growing ?
                :sigma_overlap : case.method, sigma_over_c=case.sigma,
            pps_n=case.pps, pps_overlap=case.overlap)
        result = run_case(config)
        push!(results, merge(case, (; result)))
        return nothing
    end
    for case in cases
        execute_case(case)
    end

    # Fixed overlap is only a marginal theory path. If its two finest
    # refinements fail to reduce both observables, follow h ∝ σ^1.5, anchored
    # at overlap=1.3 for σ/c=0.3. SigmaOverlap derives and capacity-limits the
    # particle count through the configured ParticleField allocation.
    fixed_rows = [_sps_case_row(case.result, control, case.method, case.sigma,
        case.pps, case.overlap, case.buffer) for case in results
        if case.method == :sigma_overlap]
    needs_growing = false
    if !smoke && method_mode in (:all, :sigma_overlap)
        for buffer in buffers
            row_at(sigma) = only([r for r in fixed_rows
                if r.buffer_over_c == buffer && r.sigma_over_c == sigma])
            r008, r015, r030 = row_at(0.08), row_at(0.15), row_at(0.3)
            improves = abs(r015.delta_CL_percent) < abs(r030.delta_CL_percent) &&
                r015.gamma_error_percent < r030.gamma_error_percent &&
                abs(r008.delta_CL_percent) < abs(r015.delta_CL_percent) &&
                r008.gamma_error_percent < r015.gamma_error_percent
            needs_growing |= !improves
        end
    end
    if needs_growing
        growing_prediction = get(ENV, "SSWSPS_GROWING_PREDICTION",
            joinpath("data", "ssw_representation_probe",
                "phase_b_prediction_sigma_overlap_growing.csv"))
        isfile(growing_prediction) || throw(ArgumentError(
            "growing-overlap path triggered; create its immutable frozen " *
            "prediction first with examples/ssw_growing_overlap_prediction.jl"))
        println("fixed-overlap fine points failed joint CL/Gamma reduction; " *
            "adding h proportional to sigma^1.5")
        for sigma in sigmas, buffer in buffers
            growing_overlap = overlap * sqrt(0.3 / sigma)
            case = (; method=:sigma_overlap_growing, sigma, pps=primary_pps,
                overlap=growing_overlap, buffer)
            push!(cases, case)
            execute_case(case)
        end
    end

    summary_rows = [_sps_case_row(case.result, control, case.method, case.sigma,
        case.pps, case.overlap, case.buffer) for case in results]
    csv_path = joinpath(base.output_root, "split_summary.csv")
    _write_sps_summary(summary_rows, csv_path)

    # plot: tail ratio vs sheet-buffer length, one series per overlap
    plot_path = nothing
    if SSW_PLOTTING
        fig = plt.figure("ssw_split", figsize=(7, 4.5))
        fig.clear()
        ax = fig.add_subplot(1, 1, 1)
        for method in (:sigma_pps, :sigma_overlap), sigma in sigmas
            pts = sort([(r.buffer_over_c, r.delta_CL_percent) for r in summary_rows
                if r.method == method && r.sigma_over_c == sigma &&
                   (method == :sigma_overlap || r.pps == primary_pps)])
            isempty(pts) && continue
            ax.semilogx(first.(pts), last.(pts), "o-";
                label="$(method), σ/c=$(sigma)")
        end
        ax.axhspan(-0.5, 0.5; color="k", alpha=0.08, label="CL tolerance")
        ax.set(xlabel="sheet-buffer length before particles (chords)",
            ylabel="tail ΔCL vs panel (%)",
            title="AR=6 wing: sheet-vs-particle split, eta=$(base.eta), dt*=$(base.dt_star)")
        ax.grid(true; which="both", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        plot_path = joinpath(base.output_root, "sheet_particle_split.png")
        fig.savefig(plot_path; dpi=180, bbox_inches="tight")
    end

    println("\nadmissibility summary:")
    @printf("%14s %8s %6s %8s %10s %10s %9s\n",
        "method", "sigma/c", "pps", "buffer", "dCL(%)", "dGamma(%)", "admiss.")
    for row in summary_rows
        @printf("%14s %8.3g %6d %8.3g %+10.4f %10.4f %9s\n",
            string(row.method), row.sigma_over_c, row.pps, row.buffer_over_c,
            row.delta_CL_percent, row.gamma_error_percent, string(row.admissible))
    end
    return (; control, results, summary_rows, csv_path, plot_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_ssw_sheet_particle_split()
end
