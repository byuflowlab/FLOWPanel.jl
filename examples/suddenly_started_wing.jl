#=##############################################################################
# DESCRIPTION
#   Wagner-function validation using a suddenly started, open-tip NACA 0012
#   wing and FLOWPanel's unsteady Neumann vortex-ring formulation.
#
#   The reference case is VortexLattice.jl/test/suddenly_started_wing_uvlm.jl.
#   Its executable timestep is U*dt/c = 1/8 (the 1/16 header comment is stale).
#
# RUN
#   SSW_MODE=single julia --project examples/suddenly_started_wing.jl
#   SSW_MODE=convergence julia --project -t auto examples/suddenly_started_wing.jl
#
# ENVIRONMENT
#   SSW_MODE       single | coarse | convergence (default: single)
#   SSW_OUTPUT     output root (default: data/suddenly_started_wing)
#   SSW_BACKEND    direct | fmm (default: fmm)
#   SSW_FMM_ORDER  FMM expansion order (default: 10)
#   SSW_FMM_THETA  FMM multipole acceptance threshold (default: 0.4)
#   SSW_FMM_LEAF   FMM source leaf size (default: 100)
#   SSW_NAIRFOIL   contour request for single mode (default: 21)
#   SSW_NSPAN      spanwise strips for single mode (default: 12)
#   SSW_DT_STAR    U*dt/c for single mode (default: 0.125)
#   SSW_T_END_STAR final tU/c (default: 7)
#   SSW_BACKSLASH_MAX_PANELS dense-solver cutoff (default: 10000)
#   SSW_RESUME     reuse completed convergence histories (default: true)
#   SSW_VERBOSE    true | false (default: true)
#   SAVE_VTK       true | false (default: true)
=###############################################################################

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra as LA
import Printf: @printf
import PythonPlot as plt

const DEFAULT_SSW_OUTPUT = joinpath("data", "suddenly_started_wing")
const SSW_GRAD_MU_OPTIONS = (; basis=:tri)

Base.@kwdef struct SSWConfig
    c::Float64 = 1.0
    AR::Float64 = 100.0
    U::Float64 = 1.0
    alpha_deg::Float64 = 5.0
    rho::Float64 = 1.0
    t_end_star::Float64 = 7.0
    dt_star::Float64 = 1 / 8
    n_airfoil::Int = 21
    n_span::Int = 12
    eta::Float64 = 0.25
    wake_core_over_c::Float64 = 1e-3
    kerneloffset_over_c::Float64 = 1e-6
    kernelcutoff_over_c::Float64 = 1e-12
    output_root::String = DEFAULT_SSW_OUTPUT
    save_vtk::Bool = true
    backend_kind::Symbol = :fmm
    fmm_expansion_order::Int = 10
    fmm_multipole_acceptance::Float64 = 0.4
    fmm_leaf_size::Int = 100
    verbose::Bool = true
    # Dense factorization is substantially faster and more robust than the
    # current unpreconditioned matrix-free solve for this extreme-AR geometry.
    backslash_max_panels::Int = 10_000
end

"Jones' two-exponential approximation to Wagner's function, in semichord time."
ssw_wagner(tstar::Real) = 1 - 0.165exp(-0.09tstar) - 0.335exp(-0.6tstar)

function ssw_time_range(config::SSWConfig)
    config.dt_star > 0 || throw(ArgumentError("dt_star must be positive"))
    nsteps = round(Int, config.t_end_star / config.dt_star)
    isapprox(nsteps * config.dt_star, config.t_end_star; atol=100eps(Float64), rtol=0) ||
        throw(ArgumentError("t_end_star must be an integer multiple of dt_star"))
    dt = config.dt_star * config.c / config.U
    return range(0.0; step=dt, length=nsteps + 1)
end

function _ssw_directions(config::SSWConfig)
    alpha = deg2rad(config.alpha_deg)
    drag = SVector{3}(cos(alpha), 0.0, sin(alpha))
    span = SVector{3}(0.0, 1.0, 0.0)
    lift = LA.cross(drag, span)
    return drag, lift
end

"Closed-trailing-edge NACA 0012 contour, ordered TE-lower-LE-upper-TE."
function naca0012_contour(n::Integer=21; thickness::Real=0.12)
    n >= 21 || throw(ArgumentError("n_airfoil must be at least 21"))
    isodd(n) || throw(ArgumentError("n_airfoil must be odd"))
    n_half = cld(n, 2)
    beta = range(0.0, pi; length=n_half)
    x = 0.5 .* (1 .- cos.(beta))
    yt = 5thickness .* (0.2969 .* sqrt.(x) .- 0.1260 .* x .-
        0.3516 .* x.^2 .+ 0.2843 .* x.^3 .- 0.1036 .* x.^4)
    lower = hcat(reverse(x), -reverse(yt))
    upper = hcat(x[2:end-1], yt[2:end-1])
    return vcat(lower, upper)
end

"Triangular open-tip extrusion of a NACA 0012 contour."
function suddenly_started_wing_mesh(c::Real, b::Real;
        n_span::Integer=12, n_airfoil::Integer=21)
    n_span >= 1 || throw(ArgumentError("n_span must be positive"))
    contour = naca0012_contour(n_airfoil)
    n_section = size(contour, 1)
    ys = range(-b / 2, b / 2; length=n_span + 1)
    node_index(i, j) = i + (j - 1) * n_section

    nodes = zeros(Float64, 3, n_section * length(ys))
    for (j, y) in enumerate(ys), i in 1:n_section
        nodes[:, node_index(i, j)] .= (c * contour[i, 1], y, c * contour[i, 2])
    end

    cells = zeros(Int, 3, 2 * n_section * n_span)
    k = 0
    for j in 1:n_span, i in 1:n_section
        ip = i == n_section ? 1 : i + 1
        n11, n21 = node_index(i, j), node_index(ip, j)
        n12, n22 = node_index(i, j + 1), node_index(ip, j + 1)
        cells[:, k += 1] .= (n11, n21, n22)
        cells[:, k += 1] .= (n11, n22, n12)
    end
    return nodes, cells
end

function _ssw_shedding(nodes, cells, c)
    tol = max(1e-8 * abs(c), 100eps(Float64) * max(abs(c), 1.0))
    te_nodes = findall(i -> isapprox(nodes[1, i], c; atol=tol, rtol=0), axes(nodes, 2))
    sort!(te_nodes; by=i -> nodes[2, i])
    length(te_nodes) >= 2 || error("could not identify the trailing-edge chain")
    lower = [c - tol, minimum(nodes[2, te_nodes]) - tol, minimum(nodes[3, te_nodes]) - tol]
    upper = [c + tol, maximum(nodes[2, te_nodes]) + tol, maximum(nodes[3, te_nodes]) + tol]
    return pnl.calc_shedding_from_seed(nodes, cells, te_nodes[1], te_nodes[2];
        end_node=te_nodes[end], bbox=(lower, upper), normal_jump_tol=1.0,
        max_turn_angle=pi / 2)
end

function build_suddenly_started_wing(config::SSWConfig; semiinfinite_wake::Bool=false)
    b = config.AR * config.c
    nodes, cells = suddenly_started_wing_mesh(config.c, b;
        n_span=config.n_span, n_airfoil=config.n_airfoil)
    bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
    options = (;
        kerneloffset=config.kerneloffset_over_c * config.c,
        kernelcutoff=config.kernelcutoff_over_c * config.c,
        semiinfinite_wake,
        watertight=false,
    )

    # The constructor may re-wind cells. Derive shedding from the constructed
    # base body, never from the raw mesh.
    base = bodytype(nodes, cells, zeros(Int, 6, 0); options...)
    shedding = _ssw_shedding(base.nodes, base.cells, config.c)
    return bodytype(copy(base.nodes), copy(base.cells), [shedding]; options...)
end

function _set_ssw_Das!(body, displacement)
    for Das in body.Das, j in axes(Das, 2)
        Das[:, j] .= displacement
    end
    return body
end

function _ssw_backend(config::SSWConfig)
    config.backend_kind == :direct && return pnl.DirectBackend()
    if config.backend_kind == :fmm
        config.fmm_expansion_order > 0 ||
            throw(ArgumentError("fmm_expansion_order must be positive"))
        config.fmm_multipole_acceptance >= 0 ||
            throw(ArgumentError("fmm_multipole_acceptance must be nonnegative"))
        config.fmm_leaf_size > 0 || throw(ArgumentError("fmm_leaf_size must be positive"))
        return pnl.FastMultipoleBackend(
            expansion_order=config.fmm_expansion_order,
            multipole_acceptance=config.fmm_multipole_acceptance,
            leaf_size=config.fmm_leaf_size,
        )
    end
    throw(ArgumentError(
        "backend_kind must be :direct or :fmm; got $(config.backend_kind)"))
end

function _ssw_solver(body, backend, max_backslash::Int)
    if body.ncells <= max_backslash
        return pnl.Backslash(body), :backslash
    end
    return pnl.KrylovSolver(body; backend, method=:gmres, atol=1e-9,
        rtol=1e-9, itmax=1000), :krylov
end

function prepare_suddenly_started_wing(config::SSWConfig)
    t_range = ssw_time_range(config)
    drag_hat, lift_hat = _ssw_directions(config)
    full_uinf = config.U * drag_hat
    Uinf(t) = t <= first(t_range) ? zero(full_uinf) : full_uinf
    maneuver!(frames, systems, wakes, t) = nothing

    wing = build_suddenly_started_wing(config; semiinfinite_wake=false)
    dt = step(t_range)
    _set_ssw_Das!(wing, config.eta * dt * full_uinf)
    wake = pnl.PanelWake(wing;
        nwakerows=length(t_range) - 1,
        core_size=config.wake_core_over_c * config.c,
        include_final_filament=false,
        shed_with_induced_velocity=true)
    frames = pnl.ReferenceFrame(wing)
    backend = _ssw_backend(config)
    solver, solver_name = _ssw_solver(wing, backend, config.backslash_max_panels)
    pressure = pnl.PressureBernoulli(config.rho; unsteady=true,
        backend, correct_kuttacondition=false)
    force = pnl.ForceMonitor(length(t_range), 1;
        normalization=pnl.NoNormalization(), i_frame=-1,
        correct_kuttacondition=false, verbose=false)

    # FLOWPanel solves at every sample, including t=0. Prime Bernoulli history
    # with a genuinely zero-flow body solve, then give the unattached wake row
    # the first interval's full convection velocity. After propagation/shedding,
    # update_TE! resets the new upstream edge while the downstream edge retains
    # this displacement, avoiding a degenerate first wake panel.
    function initialize_wake_convection!(systems, wakes, frames, uinf, i_step, dt)
        i_step == 0 || return nothing
        panel_wake = wakes[1]
        for velocity in panel_wake.velocity
            velocity[:, 1, :] .= full_uinf
        end
        return nothing
    end
    monitors = (pressure, force, initialize_wake_convection!)
    return (; wing, wake, frames, maneuver!, Uinf, full_uinf, t_range,
        drag_hat, lift_hat, backend, solver, solver_name, pressure, force,
        initialize_wake_convection!, monitors)
end

function _ssw_steady_cl(config::SSWConfig, backend; path=nothing, name="steady")
    wing = build_suddenly_started_wing(config; semiinfinite_wake=true)
    drag_hat, lift_hat = _ssw_directions(config)
    full_uinf = config.U * drag_hat
    _set_ssw_Das!(wing, full_uinf / LA.norm(full_uinf))
    solver, solver_name = _ssw_solver(wing, backend, config.backslash_max_panels)
    pressure = pnl.PressureBernoulli(config.rho; correct_kuttacondition=false,
        backend)
    force = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization(),
        i_frame=-1, correct_kuttacondition=false, verbose=false)
    elapsed = @elapsed pnl.steady!(wing, pnl.ReferenceFrame(wing), full_uinf;
        body_solvers=solver, backend, monitors=(pressure, force), path, name,
        grad_mu_options=SSW_GRAD_MU_OPTIONS, verbose=false)
    qS = 0.5 * config.rho * config.U^2 * config.AR * config.c^2
    cl = LA.dot(force.force[:, 1], lift_hat) / qS
    cd = LA.dot(force.force[:, 1], drag_hat) / qS
    return (; cl, cd, elapsed, solver_name)
end

function _ssw_case_tag(config::SSWConfig)
    dt_tag = replace(string(config.dt_star), "." => "p")
    return "na$(config.n_airfoil)_ns$(config.n_span)_dt$(dt_tag)_$(config.backend_kind)"
end

function ssw_with(config::SSWConfig; kwargs...)
    names = fieldnames(SSWConfig)
    values = ntuple(i -> getfield(config, names[i]), length(names))
    base = NamedTuple{names}(values)
    return SSWConfig(; merge(base, (; kwargs...))...)
end

function _write_ssw_case_csv(result)
    mkpath(result.path)
    filename = joinpath(result.path, "history.csv")
    open(filename, "w") do io
        println(io, "time_star,CL,CD,CL_over_CLsteady,wagner,error")
        for i in eachindex(result.time_star)
            @printf(io, "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                result.time_star[i], result.CL[i], result.CD[i], result.ratio[i],
                result.wagner[i], result.ratio[i] - result.wagner[i])
        end
    end
    return filename
end

function load_ssw_result(config::SSWConfig)
    tag = _ssw_case_tag(config)
    path = joinpath(config.output_root, tag)
    csv_path = joinpath(path, "history.csv")
    isfile(csv_path) || throw(ArgumentError("no completed history at $csv_path"))
    rows = [parse.(Float64, split(line, ',')) for line in readlines(csv_path)[2:end]]
    isempty(rows) && error("empty history at $csv_path")
    time_star = getindex.(rows, 1)
    CL = getindex.(rows, 2)
    CD = getindex.(rows, 3)
    ratio = getindex.(rows, 4)
    wagner = getindex.(rows, 5)
    error = ratio .- wagner
    steady_candidates = [cl / r for (cl, r) in zip(CL, ratio)
        if isfinite(cl) && isfinite(r) && abs(r) > eps(Float64)]
    steady_CL = isempty(steady_candidates) ? NaN : first(steady_candidates)
    panels = 2 * (config.n_airfoil - 1) * config.n_span
    solver_name = panels <= config.backslash_max_panels ? :backslash : :krylov
    return (; config, tag, path, panels, solver_name, time_star, CL, CD, ratio,
        wagner, steady_CL, steady_CD=NaN,
        rms_error=LA.norm(error) / sqrt(length(error)),
        max_error=maximum(abs, error), elapsed=NaN, steady_elapsed=NaN,
        wake_rows=length(time_star), csv_path)
end

function run_suddenly_started_wing(config::SSWConfig=SSWConfig())
    sim = prepare_suddenly_started_wing(config)
    tag = _ssw_case_tag(config)
    case_path = joinpath(config.output_root, tag)
    vtk_path = config.save_vtk ? case_path : nothing

    config.verbose && println("Suddenly-started wing: $tag, panels=$(sim.wing.ncells), solver=$(sim.solver_name)")
    steady = _ssw_steady_cl(config, sim.backend;
        path=vtk_path, name=tag * "_steady")
    elapsed = @elapsed pnl.simulate!((sim.wing,), (sim.wake,), sim.frames,
        sim.maneuver!, sim.Uinf, sim.t_range;
        body_solvers=(sim.solver,), backend=sim.backend,
        monitors=sim.monitors, path=vtk_path, name=tag,
        set_Das_eta_freestream=NaN, grad_mu_options=SSW_GRAD_MU_OPTIONS,
        verbose=config.verbose)

    qS = 0.5 * config.rho * config.U^2 * config.AR * config.c^2
    CL_all = vec(sim.lift_hat' * sim.force.force) ./ qS
    CD_all = vec(sim.drag_hat' * sim.force.force) ./ qS
    # VortexLattice reports one result per completed interval, at dt:dt:t_end.
    indices = 2:length(sim.t_range)
    time_star = collect(sim.t_range[indices]) .* config.U ./ config.c
    CL = CL_all[indices]
    CD = CD_all[indices]
    ratio = CL ./ steady.cl
    wagner = ssw_wagner.(time_star)
    error = ratio .- wagner
    rms_error = LA.norm(error) / sqrt(length(error))
    max_error = maximum(abs, error)
    result = (; config, tag, path=case_path, panels=sim.wing.ncells,
        solver_name=sim.solver_name, time_star, CL, CD, ratio, wagner,
        steady_CL=steady.cl, steady_CD=steady.cd, rms_error, max_error,
        elapsed, steady_elapsed=steady.elapsed, wake_rows=sim.wake.nwakes[])
    csv_path = _write_ssw_case_csv(result)

    @printf("  CLsteady=%+.8f, Wagner RMS=%8.4e, max=%8.4e, elapsed=%.2fs\n",
        result.steady_CL, result.rms_error, result.max_error, result.elapsed)
    return merge(result, (; csv_path))
end

function ssw_curve_change(coarse, fine)
    coarse.time_star == fine.time_star || throw(ArgumentError("curve times differ"))
    delta = fine.ratio .- coarse.ratio
    scale = max.(abs.(fine.ratio), sqrt(eps(Float64)))
    return (; max_absolute=maximum(abs, delta),
        max_relative=maximum(abs.(delta) ./ scale),
        relative_l2=LA.norm(delta) / max(LA.norm(fine.ratio), eps(Float64)))
end

function _write_ssw_summary(results, output_root)
    mkpath(output_root)
    filename = joinpath(output_root, "convergence.csv")
    open(filename, "w") do io
        println(io, "tag,backend,n_airfoil,n_span,dt_star,panels,solver,CLsteady,wagner_rms,wagner_max,successive_max_absolute,successive_max_relative,successive_relative_l2,elapsed")
        previous = Dict{Tuple{Symbol, Float64}, Any}()
        for r in sort(collect(results); by=x -> (string(x.config.backend_kind), x.config.dt_star, x.panels))
            key = (r.config.backend_kind, r.config.dt_star)
            change = haskey(previous, key) && previous[key].time_star == r.time_star ?
                ssw_curve_change(previous[key], r) :
                (; max_absolute=NaN, max_relative=NaN, relative_l2=NaN)
            @printf(io, "%s,%s,%d,%d,%.8g,%d,%s,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.8g\n",
                r.tag, string(r.config.backend_kind), r.config.n_airfoil,
                r.config.n_span, r.config.dt_star, r.panels, string(r.solver_name),
                r.steady_CL, r.rms_error, r.max_error, change.max_absolute,
                change.max_relative, change.relative_l2, r.elapsed)
            previous[key] = r
        end
    end
    return filename
end

function plot_ssw_results(results; output_root=first(results).config.output_root,
        emphasize_tag=nothing)
    isempty(results) && return nothing
    mkpath(output_root)
    fig = plt.figure("ssw_wagner", figsize=(11, 4.5))
    fig.clear()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    tfine = range(0.0, maximum(first(results).time_star); length=400)
    ax1.plot(tfine, ssw_wagner.(tfine), "k:"; linewidth=2, label="Wagner")
    for (i, r) in enumerate(results)
        label = "$(r.config.n_airfoil)/$(r.config.n_span), $(r.config.backend_kind), dt*=$(r.config.dt_star)"
        emphasized = isnothing(emphasize_tag) ? i == length(results) : r.tag == emphasize_tag
        ax1.plot(r.time_star, r.ratio; linewidth=emphasized ? 2.2 : 1.0,
            label)
        ax2.plot(r.time_star, r.ratio .- r.wagner;
            linewidth=emphasized ? 2.2 : 1.0, label)
    end
    ax1.set(xlabel="t*=t Uinf/c", ylabel="CL/CLsteady",
        title="Suddenly-started NACA 0012 wing")
    ax2.set(xlabel="t*=t Uinf/c", ylabel="FLOWPanel - Wagner",
        title="Validation error")
    for ax in (ax1, ax2)
        ax.grid(true; alpha=0.3)
        ax.legend(fontsize=7, loc="best")
    end
    fig.tight_layout()
    filename = joinpath(output_root, "wagner_comparison.png")
    fig.savefig(filename; dpi=180, bbox_inches="tight")
    return filename
end

function plot_ssw_convergence(results; output_root=first(results).config.output_root,
        backend_kind=first(results).config.backend_kind,
        dt_star=first(results).config.dt_star)
    spatial = [r for r in results if r.config.backend_kind == backend_kind &&
        r.config.dt_star == dt_star]
    sort!(spatial; by=r -> r.panels)
    # Keep only the finest-dimensional case at duplicate panel counts.
    spatial = [spatial[i] for i in eachindex(spatial)
        if i == length(spatial) || spatial[i + 1].panels != spatial[i].panels]
    isempty(spatial) && return nothing
    panels = getproperty.(spatial, :panels)
    rms = getproperty.(spatial, :rms_error)
    maxerr = getproperty.(spatial, :max_error)
    successive_absolute = fill(NaN, length(spatial))
    successive_l2 = fill(NaN, length(spatial))
    for i in 2:length(spatial)
        spatial[i - 1].time_star == spatial[i].time_star || continue
        change = ssw_curve_change(spatial[i - 1], spatial[i])
        successive_absolute[i] = change.max_absolute
        successive_l2[i] = change.relative_l2
    end

    fig = plt.figure("ssw_convergence", figsize=(10, 4.5))
    fig.clear()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    ax1.loglog(panels, rms, "o-"; label="Wagner RMS")
    ax1.loglog(panels, maxerr, "s-"; label="Wagner max")
    ax2.semilogx(panels[2:end], 100 .* successive_absolute[2:end], "o-";
        label="successive max absolute")
    ax2.semilogx(panels[2:end], 100 .* successive_l2[2:end], "s-";
        label="successive relative L2")
    ax2.axhline(2.0; color="k", linestyle=":", label="2% target")
    ax1.set(xlabel="surface panels", ylabel="error", title="Wagner error")
    ax2.set(xlabel="surface panels", ylabel="percent", title="Mesh convergence")
    for ax in (ax1, ax2)
        ax.grid(true; alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    end
    fig.tight_layout()
    filename = joinpath(output_root, "convergence.png")
    fig.savefig(filename; dpi=180, bbox_inches="tight")
    return filename
end

function run_ssw_convergence(config::SSWConfig=SSWConfig(); tolerance=0.02,
        span_values=(12, 24, 48, 96, 192), airfoil_values=(21, 41, 81),
        resume::Bool=true)
    cache = Dict{Tuple{Int, Int, Float64, Symbol}, Any}()
    results = Any[]
    run_case(na, ns, dtstar=config.dt_star, kind=config.backend_kind) =
            get!(cache, (na, ns, Float64(dtstar), kind)) do
        c = ssw_with(config; n_airfoil=na, n_span=ns, dt_star=dtstar,
            backend_kind=kind)
        history_path = joinpath(c.output_root, _ssw_case_tag(c), "history.csv")
        r = resume && isfile(history_path) ? load_ssw_result(c) :
            run_suddenly_started_wing(c)
        push!(results, r)
        _write_ssw_summary(results, config.output_root)
        plot_ssw_results(results; output_root=config.output_root)
        r
    end

    direct_crosscheck = nothing
    backend_change = nothing
    if config.backend_kind == :fmm
        direct_crosscheck = run_case(first(airfoil_values), first(span_values),
            config.dt_star, :direct)
        fmm_coarse = run_case(first(airfoil_values), first(span_values))
        backend_change = ssw_curve_change(direct_crosscheck, fmm_coarse)
    end

    previous = nothing
    accepted_span = last(span_values)
    for ns in span_values
        current = run_case(first(airfoil_values), ns)
        if previous !== nothing
            change = ssw_curve_change(previous, current)
            @printf("  span %d -> %d: max relative=%7.3f%%, L2=%7.3f%%\n",
                previous.config.n_span, ns, 100 * change.max_relative,
                100 * change.relative_l2)
            if change.max_absolute <= tolerance && change.relative_l2 <= tolerance
                accepted_span = ns
                break
            end
        end
        previous = current
    end

    previous = nothing
    accepted_airfoil = last(airfoil_values)
    for na in airfoil_values
        current = run_case(na, accepted_span)
        if previous !== nothing
            change = ssw_curve_change(previous, current)
            @printf("  airfoil %d -> %d: max relative=%7.3f%%, L2=%7.3f%%\n",
                previous.config.n_airfoil, na, 100 * change.max_relative,
                100 * change.relative_l2)
            if change.max_absolute <= tolerance && change.relative_l2 <= tolerance
                accepted_airfoil = na
                break
            end
        end
        previous = current
    end

    accepted = run_case(accepted_airfoil, accepted_span)
    confirm = run_case(accepted_airfoil, 2 * accepted_span)
    joint_change = ssw_curve_change(accepted, confirm)
    joint_converged = joint_change.max_absolute <= tolerance &&
        joint_change.relative_l2 <= tolerance
    accepted = confirm
    if !joint_converged
        @warn "joint span confirmation exceeded tolerance" joint_change
    end

    temporal = run_case(accepted.config.n_airfoil, accepted.config.n_span,
        accepted.config.dt_star / 2)
    temporal_ratio = temporal.ratio[2:2:end]
    length(temporal_ratio) == length(accepted.ratio) ||
        error("half-timestep history does not align with baseline")
    temporal_delta = temporal_ratio .- accepted.ratio
    temporal_change = (; max_absolute=maximum(abs, temporal_delta),
        max_relative=maximum(abs.(temporal_delta) ./
            max.(abs.(temporal_ratio), sqrt(eps(Float64)))),
        relative_l2=LA.norm(temporal_delta) / LA.norm(temporal_ratio))

    summary_path = _write_ssw_summary(results, config.output_root)
    plot_path = plot_ssw_results(results; output_root=config.output_root,
        emphasize_tag=accepted.tag)
    convergence_plot_path = plot_ssw_convergence(results;
        output_root=config.output_root, backend_kind=config.backend_kind,
        dt_star=config.dt_star)
    return (; results, accepted, confirm, temporal, joint_change,
        joint_converged, temporal_change, direct_crosscheck, backend_change,
        summary_path, plot_path, convergence_plot_path)
end

_envbool(name, default) = lowercase(get(ENV, name, string(default))) in
    ("1", "true", "yes", "on")

function _ssw_config_from_env()
    backend = Symbol(lowercase(get(ENV, "SSW_BACKEND", "fmm")))
    return SSWConfig(
        n_airfoil=parse(Int, get(ENV, "SSW_NAIRFOIL", "21")),
        n_span=parse(Int, get(ENV, "SSW_NSPAN", "12")),
        dt_star=parse(Float64, get(ENV, "SSW_DT_STAR", "0.125")),
        t_end_star=parse(Float64, get(ENV, "SSW_T_END_STAR", "7.0")),
        output_root=get(ENV, "SSW_OUTPUT", DEFAULT_SSW_OUTPUT),
        save_vtk=_envbool("SAVE_VTK", true),
        backend_kind=backend,
        fmm_expansion_order=parse(Int, get(ENV, "SSW_FMM_ORDER", "10")),
        fmm_multipole_acceptance=parse(Float64,
            get(ENV, "SSW_FMM_THETA", "0.4")),
        fmm_leaf_size=parse(Int, get(ENV, "SSW_FMM_LEAF", "100")),
        verbose=_envbool("SSW_VERBOSE", true),
        backslash_max_panels=parse(Int,
            get(ENV, "SSW_BACKSLASH_MAX_PANELS", "10000")),
    )
end

function main()
    mode = Symbol(lowercase(get(ENV, "SSW_MODE", "single")))
    config = _ssw_config_from_env()
    if mode == :single
        result = run_suddenly_started_wing(config)
        plot_ssw_results([result]; output_root=config.output_root)
        return result
    elseif mode == :coarse
        coarse = ssw_with(config; n_airfoil=21, n_span=12)
        result = run_suddenly_started_wing(coarse)
        plot_ssw_results([result]; output_root=config.output_root)
        return result
    elseif mode == :convergence
        return run_ssw_convergence(config; resume=_envbool("SSW_RESUME", true))
    end
    error("SSW_MODE must be single, coarse, or convergence; got $mode")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
