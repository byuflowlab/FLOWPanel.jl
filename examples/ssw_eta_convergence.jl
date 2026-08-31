#=##############################################################################
# DESCRIPTION
#   BRAINSTORM 014: eta x dt convergence matrix for the first-wake-row offset
#   on the suddenly-started wing (Das = eta*dt*U), at finite aspect ratio so
#   R.T. Jones' finite-AR indicial function (NACA Rept. 681) provides an
#   external sanity band.
#
#   Question this answers: is eta ordinary truncation error (all eta collapse
#   to one curve as dt -> 0) or a model parameter (curves converge in dt to
#   eta-dependent limits)? The PRIMARY criterion is self-convergence across
#   the matrix; the Jones curve is a sanity band only (elliptic wing vs this
#   rectangular mesh, and lifting-line-approximate).
#
# RUN
#   julia --project -t 6 examples/ssw_eta_convergence.jl
#
# ENVIRONMENT
#   SSWETA_AR        aspect ratio (default: 6)
#   SSWETA_ETAS      comma list (default: 0.125,0.25,0.5,1.0,2.0,4.0)
#   SSWETA_DTS       comma list of U*dt/c, halving ladder recommended
#                    (default: 0.25,0.125,0.0625,0.03125)
#   SSWETA_NSPAN     spanwise strips (default: 24)
#   SSWETA_NAIRFOIL  airfoil contour points (default: 21)
#   SSWETA_BACKEND   direct | fmm (default: direct)
#   SSWETA_TEND      final tU/c (default: 7)
#   SSWETA_OUTPUT    output root (default: data/ssw_eta_convergence)
#   SSWETA_MESHCHECK run span/airfoil refinement at one matrix corner
#                    (default: true)
#   SSWETA_TMIN      exclude t* < TMIN from all comparison metrics (default: 1;
#                    the first ~2 samples carry an impulsive-start ring in the
#                    unsteady Bernoulli that would otherwise dominate them)
#   SSWETA_RESUME    reuse completed histories (default: true). NOTE the case
#                    tag carries SSW_MESH_REVISION, so a mesh change cannot
#                    silently resume histories from the previous geometry.
#   SSWETA_LOGLAW_EXTRA
#                    comma list of "eta:dt_star" completed cases outside the
#                    matrix to fold into the log-law fit only, resolved against
#                    SSWETA_OUTPUT (e.g. "8:0.25,16:0.25" for large-Das probes)
=###############################################################################

if !isdefined(@__MODULE__, :SSWConfig)
    include(joinpath(@__DIR__, "suddenly_started_wing.jl"))
end

import Printf: @printf, @sprintf

_parse_list(name, default) = parse.(Float64, split(get(ENV, name, default), ','))

function _sswe_base_config()
    backend = Symbol(lowercase(get(ENV, "SSWETA_BACKEND", "direct")))
    return SSWConfig(
        AR=parse(Float64, get(ENV, "SSWETA_AR", "6.0")),
        n_span=parse(Int, get(ENV, "SSWETA_NSPAN", "24")),
        n_airfoil=parse(Int, get(ENV, "SSWETA_NAIRFOIL", "21")),
        t_end_star=parse(Float64, get(ENV, "SSWETA_TEND", "7.0")),
        backend_kind=backend,
        output_root=get(ENV, "SSWETA_OUTPUT",
            joinpath("data", "ssw_eta_convergence")),
        save_vtk=false,
        verbose=_envbool("SSW_VERBOSE", false),
    )
end

"Run (or resume) one case; returns the standard SSW result named tuple."
function _sswe_case(base, eta, dt_star; resume=true)
    c = ssw_with(base; eta, dt_star)
    history = joinpath(c.output_root, _ssw_case_tag(c), "history.csv")
    if resume && isfile(history)
        return load_ssw_result(c)
    end
    return run_suddenly_started_wing(c)
end

_sswe_tmin() = parse(Float64, get(ENV, "SSWETA_TMIN", "1.0"))

"""
Curve-difference metrics of the CL/CL_steady ratio between two cases, on the
coarser case's time grid restricted to t* >= tmin. The dt of `a` must be an
integer multiple of the dt of `b` (equal dt allowed).
"""
function _sswe_delta(a, b; tmin=_sswe_tmin())
    coarse, fine = a.config.dt_star >= b.config.dt_star ? (a, b) : (b, a)
    stride = coarse.config.dt_star / fine.config.dt_star
    isinteger(stride) || throw(ArgumentError(
        "dt ratio $(stride) is not an integer; use a halving ladder"))
    idx = Int(stride):Int(stride):length(fine.ratio)
    length(idx) == length(coarse.ratio) || throw(ArgumentError(
        "subsampled fine history does not align with coarse history"))
    coarse.time_star ≈ fine.time_star[idx] || throw(ArgumentError(
        "subsampled times do not match"))
    keep = coarse.time_star .>= tmin
    any(keep) || throw(ArgumentError("tmin=$(tmin) leaves no samples"))
    delta = fine.ratio[idx][keep] .- coarse.ratio[keep]
    ref = fine.ratio[idx][keep]
    return (; max_absolute=maximum(abs, delta),
        relative_l2=LA.norm(delta) / max(LA.norm(ref), eps(Float64)))
end

"RMS/max deviation from the case's reference curve over t* >= tmin."
function _sswe_ref_error(r; tmin=_sswe_tmin())
    keep = r.time_star .>= tmin
    err = r.ratio[keep] .- r.wagner[keep]
    return (; rms=LA.norm(err) / sqrt(length(err)), max=maximum(abs, err))
end

_sswe_tail_mean(r; window=1.0) = begin
    t_end = maximum(r.time_star)
    sel = r.time_star .>= t_end - window
    sum(r.ratio[sel]) / count(sel)
end

"""
    _sswe_fit_slope(x, y) -> (slope, rms_residual)

Least-squares slope of `y` on `x`, plus the RMS residual about the fit. Returns
`(NaN, NaN)` for fewer than two points or degenerate `x`.
"""
function _sswe_fit_slope(x, y)
    length(x) >= 2 || return (NaN, NaN)
    xbar, ybar = sum(x) / length(x), sum(y) / length(y)
    sxx = sum((x .- xbar).^2)
    sxx > 0 || return (NaN, NaN)
    slope = sum((x .- xbar) .* (y .- ybar)) / sxx
    resid = y .- (ybar .+ slope .* (x .- xbar))
    return (slope, LA.norm(resid) / sqrt(length(resid)))
end

"""
Write `das_analysis.csv`: the three claims BRAINSTORM/014 draws from this matrix,
computed rather than hand-derived from the console table.

* `das_diagonal` — cases sharing a physical offset `Das = eta*dt_star` (constant
  along the matrix's anti-diagonals), ordered coarse dt first. Reports successive
  tail increments, their ratios, and the Richardson limit. Tests "fixed-`Das`
  converges in dt".
* `fixed_eta` — successive tail increments down the dt ladder at fixed eta.
  Tests "at fixed eta there is NO dt -> 0 limit" (the increments should not
  shrink, because halving dt also halves `Das`).
* `log_law` — least-squares fit of tail ratio against `log2(Das)` over every case
  supplied, reported as %/doubling. Tests the logarithmic `Das` dependence.

`extra` supplies cases outside the matrix (e.g. large-`Das` probes) as
`(eta, dt_star) => result`; they participate in the `log_law` fit only.
"""
function _sswe_write_das_analysis(path, etas, dts, results; extra=Dict())
    tail(r) = _sswe_tail_mean(r)
    open(path, "w") do io
        println(io, "kind,das,eta,dt_star,tail,delta,ratio,richardson,slope_pct_per_doubling,rms_resid,n")

        # --- fixed-Das diagonals -------------------------------------------
        groups = Dict{Float64, Vector{Tuple{Float64, Float64}}}()  # das => [(dt, eta)]
        for eta in etas, dt in dts
            push!(get!(groups, eta * dt, Tuple{Float64, Float64}[]), (dt, eta))
        end
        for das in sort(collect(keys(groups)))
            g = sort(groups[das]; rev=true)          # coarse dt first
            length(g) >= 2 || continue
            ts = [tail(results[(eta, dt)]) for (dt, eta) in g]
            deltas = diff(ts)
            for (i, (dt, eta)) in enumerate(g)
                delta = i == 1 ? NaN : deltas[i-1]
                ratio = i <= 2 ? NaN : deltas[i-1] / deltas[i-2]
                # Richardson limit from the last increment and its ratio, valid
                # only once a ratio exists and the series is contracting.
                rich = NaN
                if i == length(g) && length(deltas) >= 2
                    rr = deltas[end] / deltas[end-1]
                    0 < rr < 1 && (rich = ts[end] + deltas[end] * rr / (1 - rr))
                end
                @printf(io, "das_diagonal,%.8g,%.8g,%.8g,%.16e,%.16e,%.16e,%.16e,NaN,NaN,%d\n",
                    das, eta, dt, ts[i], delta, ratio, rich, length(g))
            end
        end

        # --- fixed-eta dt ladders ------------------------------------------
        for eta in etas
            ts = [tail(results[(eta, dt)]) for dt in dts]   # dts sorted coarse first
            deltas = diff(ts)
            for (i, dt) in enumerate(dts)
                delta = i == 1 ? NaN : deltas[i-1]
                ratio = i <= 2 ? NaN : deltas[i-1] / deltas[i-2]
                @printf(io, "fixed_eta,%.8g,%.8g,%.8g,%.16e,%.16e,%.16e,NaN,NaN,NaN,%d\n",
                    eta * dt, eta, dt, ts[i], delta, ratio, length(dts))
            end
        end

        # --- log law -------------------------------------------------------
        # Three views, because the tail is a function of Das only in the
        # dt -> 0 limit: at finite dt, cases sharing a Das still differ (that is
        # exactly what the das_diagonal rows measure). A naive fit over all
        # (eta, dt) therefore mixes the two dependencies and is reported last,
        # as a diagnostic only.

        # (i) dt-converged: slope of the Richardson limits vs log2(Das).
        rich_x = Float64[]; rich_y = Float64[]
        for das in sort(collect(keys(groups)))
            g = sort(groups[das]; rev=true)
            length(g) >= 3 || continue
            ts = [tail(results[(eta, dt)]) for (dt, eta) in g]
            d = diff(ts)
            rr = d[end] / d[end-1]
            0 < rr < 1 || continue
            push!(rich_x, log2(das))
            push!(rich_y, ts[end] + d[end] * rr / (1 - rr))
        end
        slope, rms = _sswe_fit_slope(rich_x, rich_y)
        @printf(io, "log_law_richardson,NaN,NaN,NaN,NaN,NaN,NaN,NaN,%.16e,%.16e,%d\n",
            100 * slope, 100 * rms, length(rich_x))

        # (ii) widest Das range at a single fixed dt (the coarsest), including
        #      any large-Das probes. Per-doubling increments are emitted as
        #      `delta` so non-constancy of the "log law" is visible.
        dt0 = maximum(dts)
        row = [(eta * dt0, tail(r)) for ((eta, dt), r) in
               merge(Dict((eta, dt) => results[(eta, dt)]
                          for eta in etas for dt in dts), extra) if dt == dt0]
        sort!(row; by=first)
        for (i, (das, t)) in enumerate(row)
            delta = i == 1 ? NaN : t - row[i-1][2]
            @printf(io, "log_law_fixed_dt,%.8g,%.8g,%.8g,%.16e,%.16e,NaN,NaN,NaN,NaN,%d\n",
                das, das / dt0, dt0, t, delta, length(row))
        end
        slope, rms = _sswe_fit_slope(log2.(first.(row)), last.(row))
        @printf(io, "log_law_fixed_dt_fit,NaN,NaN,%.8g,NaN,NaN,NaN,NaN,%.16e,%.16e,%d\n",
            dt0, 100 * slope, 100 * rms, length(row))

        # (iii) diagnostic only — mixes Das and dt dependence.
        all_cases = merge(Dict((eta, dt) => results[(eta, dt)]
                               for eta in etas for dt in dts), extra)
        xs = Float64[]; ys = Float64[]
        for ((eta, dt), r) in all_cases
            das = eta * dt
            das > 0 || continue
            push!(xs, log2(das)); push!(ys, tail(r))
        end
        slope, rms = _sswe_fit_slope(xs, ys)
        @printf(io, "log_law_all_mixed,NaN,NaN,NaN,NaN,NaN,NaN,NaN,%.16e,%.16e,%d\n",
            100 * slope, 100 * rms, length(xs))
    end
    return path
end

function run_ssw_eta_convergence()
    base = _sswe_base_config()
    etas = _parse_list("SSWETA_ETAS", "0.125,0.25,0.5,1.0,2.0,4.0")
    dts = sort(_parse_list("SSWETA_DTS", "0.25,0.125,0.0625,0.03125"); rev=true)
    resume = _envbool("SSWETA_RESUME", true)
    mkpath(base.output_root)
    _, ref_name = ssw_reference(base)
    println("eta x dt matrix at AR=$(base.AR) " *
        "(ns=$(base.n_span), na=$(base.n_airfoil), $(base.backend_kind)); " *
        "reference: $ref_name")

    # cheapest-first so early failures cost little
    results = Dict{Tuple{Float64, Float64}, Any}()   # (eta, dt) => result
    for dt in dts, eta in etas
        println("case eta=$(eta), dt*=$(dt)")
        results[(eta, dt)] = _sswe_case(base, eta, dt; resume)
    end

    # ---- per-case summary (reference metrics masked to t* >= tmin) ----
    summary_path = joinpath(base.output_root, "summary.csv")
    open(summary_path, "w") do io
        println(io, "eta,dt_star,tag,CLsteady,ref_rms,ref_max,tail_ratio,elapsed")
        for dt in dts, eta in etas
            r = results[(eta, dt)]
            e = _sswe_ref_error(r)
            @printf(io, "%.8g,%.8g,%s,%.16e,%.16e,%.16e,%.16e,%.8g\n",
                eta, dt, r.tag, r.steady_CL, e.rms, e.max,
                _sswe_tail_mean(r), r.elapsed)
        end
    end

    # ---- dt-convergence per eta (successive halvings + distance to finest) ----
    analysis_path = joinpath(base.output_root, "analysis.csv")
    open(analysis_path, "w") do io
        println(io, "kind,eta,dt_coarse,dt_fine,max_absolute,relative_l2")
        for eta in etas
            for i in 1:length(dts)-1
                ch = _sswe_delta(results[(eta, dts[i])], results[(eta, dts[i+1])])
                @printf(io, "dt_pair,%.8g,%.8g,%.8g,%.16e,%.16e\n",
                    eta, dts[i], dts[i+1], ch.max_absolute, ch.relative_l2)
            end
            for dt in dts[1:end-1]
                ch = _sswe_delta(results[(eta, dt)], results[(eta, dts[end])])
                @printf(io, "vs_finest,%.8g,%.8g,%.8g,%.16e,%.16e\n",
                    eta, dt, dts[end], ch.max_absolute, ch.relative_l2)
            end
        end
        # eta-spread at each dt: max pairwise curve difference (same grid)
        for dt in dts
            worst = (; max_absolute=0.0, relative_l2=0.0)
            for i in 1:length(etas), j in i+1:length(etas)
                ch = _sswe_delta(results[(etas[i], dt)], results[(etas[j], dt)])
                worst = (; max_absolute=max(worst.max_absolute, ch.max_absolute),
                    relative_l2=max(worst.relative_l2, ch.relative_l2))
            end
            @printf(io, "eta_spread,NaN,%.8g,%.8g,%.16e,%.16e\n",
                dt, dt, worst.max_absolute, worst.relative_l2)
        end
    end

    # ---- Das-diagonal / fixed-eta / log-law analysis ----
    # SSWETA_LOGLAW_EXTRA pulls already-completed cases outside the matrix into
    # the log-law fit only (e.g. large-Das probes run in a prior invocation),
    # given as a comma list of "eta:dt_star" resolved against this output root.
    extra = Dict{Tuple{Float64, Float64}, Any}()
    for spec in split(get(ENV, "SSWETA_LOGLAW_EXTRA", ""), ',')
        isempty(strip(spec)) && continue
        e_str, dt_str = split(strip(spec), ':')
        eta_x, dt_x = parse(Float64, e_str), parse(Float64, dt_str)
        c = ssw_with(base; eta=eta_x, dt_star=dt_x)
        history = joinpath(c.output_root, _ssw_case_tag(c), "history.csv")
        isfile(history) || error("SSWETA_LOGLAW_EXTRA case eta=$(eta_x), " *
            "dt*=$(dt_x) has no history at $history; run it first")
        extra[(eta_x, dt_x)] = load_ssw_result(c)
    end
    das_analysis_path = _sswe_write_das_analysis(
        joinpath(base.output_root, "das_analysis.csv"), etas, dts, results; extra)

    # ---- optional mesh check at one matrix corner ----
    mesh_check = nothing
    if _envbool("SSWETA_MESHCHECK", true)
        eta0, dt0 = etas[argmin(abs.(etas .- 0.25))], dts[min(2, length(dts))]
        println("mesh check at eta=$(eta0), dt*=$(dt0)")
        r0 = results[(eta0, dt0)]
        r_span = _sswe_case(ssw_with(base; n_span=2base.n_span), eta0, dt0; resume)
        r_foil = _sswe_case(ssw_with(base; n_airfoil=2base.n_airfoil - 1), eta0,
            dt0; resume)
        mesh_check = (; eta=eta0, dt_star=dt0,
            span=_sswe_delta(r0, r_span),
            airfoil=_sswe_delta(r0, r_foil))
        open(analysis_path, "a") do io
            @printf(io, "mesh_span_x2,%.8g,%.8g,%.8g,%.16e,%.16e\n", eta0, dt0,
                dt0, mesh_check.span.max_absolute, mesh_check.span.relative_l2)
            @printf(io, "mesh_airfoil_x2,%.8g,%.8g,%.8g,%.16e,%.16e\n", eta0,
                dt0, dt0, mesh_check.airfoil.max_absolute,
                mesh_check.airfoil.relative_l2)
        end
    end

    plots = _plot_ssw_eta_matrix(base, etas, dts, results)

    # ---- console verdict table ----
    println("\ntail CL/CLsteady (rows: eta; cols: dt*):")
    @printf("%8s", "eta\\dt*")
    foreach(dt -> @printf("%12.5g", dt), dts)
    println()
    for eta in etas
        @printf("%8.4g", eta)
        foreach(dt -> @printf("%12.6f", _sswe_tail_mean(results[(eta, dt)])), dts)
        println()
    end
    println("\n$(ref_name) tail value: ",
        @sprintf("%.6f", ssw_reference(base)[1](base.t_end_star)))
    return (; base, etas, dts, results, summary_path, analysis_path,
        das_analysis_path, mesh_check, plots)
end

function _plot_ssw_eta_matrix(base, etas, dts, results)
    mkpath(base.output_root)
    ref_fun, ref_name = ssw_reference(base)

    # (1) ratio curves at finest dt, all eta, vs reference
    fig = plt.figure("sswe_curves", figsize=(11, 4.5))
    fig.clear()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    dtf = dts[end]
    tfine = range(0.0, base.t_end_star; length=400)
    ax1.plot(tfine, ref_fun.(tfine), "k:"; linewidth=2, label=ref_name)
    for eta in etas
        r = results[(eta, dtf)]
        ax1.plot(r.time_star, r.ratio; linewidth=1.2, label="eta=$(eta)")
        ax2.plot(r.time_star, r.ratio .- r.wagner; linewidth=1.2,
            label="eta=$(eta)")
    end
    ax1.set(xlabel="t* = t U/c", ylabel="CL/CLsteady",
        title="AR=$(base.AR), dt*=$(dtf)")
    ax2.set(xlabel="t* = t U/c", ylabel="ratio - reference",
        title="Deviation from $ref_name")
    for ax in (ax1, ax2)
        ax.grid(true; alpha=0.3)
        ax.legend(fontsize=7, loc="best")
    end
    fig.tight_layout()
    curves_path = joinpath(base.output_root, "eta_curves.png")
    fig.savefig(curves_path; dpi=180, bbox_inches="tight")

    # (2) dt-convergence per eta + eta-spread vs dt
    fig = plt.figure("sswe_convergence", figsize=(11, 4.5))
    fig.clear()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    for eta in etas
        pair_dts = dts[1:end-1]
        deltas = [_sswe_delta(results[(eta, dts[i])],
            results[(eta, dts[i+1])]).max_absolute for i in 1:length(dts)-1]
        ax1.loglog(pair_dts, deltas, "o-"; label="eta=$(eta)")
    end
    spread = [maximum(_sswe_delta(results[(etas[i], dt)],
            results[(etas[j], dt)]).max_absolute
        for i in 1:length(etas) for j in i+1:length(etas)) for dt in dts]
    ax2.loglog(dts, spread, "s-"; color="k")
    ax1.set(xlabel="coarser dt* of pair", ylabel="max |Δ ratio|",
        title="dt self-convergence per eta")
    ax2.set(xlabel="dt*", ylabel="max pairwise |Δ ratio| across eta",
        title="eta-spread vs dt (must -> 0 if eta is truncation error)")
    for ax in (ax1, ax2)
        ax.grid(true; which="both", alpha=0.3)
    end
    ax1.legend(fontsize=7, loc="best")
    fig.tight_layout()
    conv_path = joinpath(base.output_root, "eta_dt_convergence.png")
    fig.savefig(conv_path; dpi=180, bbox_inches="tight")
    return (; curves_path, conv_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_ssw_eta_convergence()
end
