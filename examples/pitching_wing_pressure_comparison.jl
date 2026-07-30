#=##############################################################################
# DESCRIPTION
#   Compare unsteady pressure monitors on the panel-wake pitching-wing case.
#
#   The default comparison is unsteady Bernoulli against PressureLaplace using
#   corrected edge differences, corrected Hessians, and surface-velocity
#   reconstruction.  All methods run in one simulation and therefore see the
#   same body and panel wake state at every step.
#
#   Julia's thread count must be selected before launching Julia, for example:
#       JULIA_NUM_THREADS=16 OMP_NUM_THREADS=16 \
#           julia --project examples/pitching_wing_pressure_comparison.jl
#
#   OMP_NUM_THREADS is also accepted as FLOWPANEL_OMP_NUM_THREADS when it is
#   useful to set the value before FLOWPanel and its native dependencies load.
=###############################################################################

# Set this before importing FLOWPanel.  Julia itself reads JULIA_NUM_THREADS at
# process startup and cannot be reconfigured reliably from this script.
if !haskey(ENV, "OMP_NUM_THREADS") && haskey(ENV, "FLOWPANEL_OMP_NUM_THREADS")
    ENV["OMP_NUM_THREADS"] = ENV["FLOWPANEL_OMP_NUM_THREADS"]
end

import FLOWPanel as pnl
import FLOWPanel: _run_monitor!, monitor_field
import LinearAlgebra: norm
import Printf: @printf
import Statistics: mean

include(joinpath(@__DIR__, "pitching_wing.jl"))

const PRESSURE_COMPARISON_METHODS = (
    (:bernoulli, "Unsteady Bernoulli"),
    (:edge_difference, "Laplace corrected edge"),
    (:corrected_hessian, "Laplace Hessian (unsupported diagnostic)"),
    (:surface_velocity, "Laplace surface velocity"),
)

_env_bool(name, default) = parse(Bool, get(ENV, name, string(default)))
_env_int(name, default) = parse(Int, get(ENV, name, string(default)))
_env_float(name, default) = parse(Float64, get(ENV, name, string(default)))

mutable struct PressureComparisonRecorder <: pnl.AbstractMonitor
    method::Symbol
    pressure_l2_relative::Vector{Float64}
    pressure_linf::Vector{Float64}
    converged::Vector{Bool}
    iterations::Vector{Int}
    absolute_residual::Vector{Float64}
    relative_residual::Vector{Float64}
    rhs_l2::Vector{Float64}
    gradient_l2::Vector{Float64}
    reference_history::Vector{Vector{Float64}}
    laplace::Any
end

monitor_requires(::PressureComparisonRecorder) = (:P,)
monitor_provides(::PressureComparisonRecorder) = ()

function PressureComparisonRecorder(method::Symbol, nt::Integer;
                                    reference_history=Vector{Vector{Float64}}(),
                                    laplace=nothing)
    return PressureComparisonRecorder(method, fill(NaN, nt), fill(NaN, nt),
        trues(nt), zeros(Int, nt), fill(NaN, nt), fill(NaN, nt), fill(NaN, nt), fill(NaN, nt),
        reference_history, laplace)
end

function _run_monitor!(m::PressureComparisonRecorder, ctx::pnl.MonitorContext,
                       systems, wakes, frames::AbstractVector{<:pnl.ReferenceFrame}, uinf,
                       i_step::Int, dt::Real,
                       t=nothing)
    pressure = collect(monitor_field(ctx, :P, 1))
    if m.method == :bernoulli
        push!(m.reference_history, pressure)
    else
        i_step + 1 <= length(m.reference_history) || throw(ArgumentError(
            "Bernoulli pressure must be recorded before $(m.method)."))
        reference = m.reference_history[i_step + 1]
        delta = pressure .- reference
        m.pressure_l2_relative[i_step + 1] = norm(delta) / max(norm(reference), eps())
        m.pressure_linf[i_step + 1] = maximum(abs, delta)
        laplace = m.laplace
        m.converged[i_step + 1] = laplace.workspace[1].stats.solved
        m.iterations[i_step + 1] = laplace.workspace[1].stats.niter
        m.absolute_residual[i_step + 1] = laplace.absolute_residual[1]
        m.relative_residual[i_step + 1] = laplace.relative_residual[1]
        m.rhs_l2[i_step + 1] = norm(laplace.b[1])
        if laplace.gradient_mode == :surface_velocity
            m.gradient_l2[i_step + 1] = norm(laplace.surface_velocity_gradient[1])
        elseif laplace.gradient_mode == :corrected_hessian
            m.gradient_l2[i_step + 1] = sqrt(norm(systems[1].velocity_gradient)^2 +
                                                 norm(laplace.jump_velocity_gradient[1])^2)
        end
    end
    return nothing
end

function _comparison_csv(path, t_range, period, alpha_history, forces,
                         pressure_metrics, diagnostics;
                         methods=PRESSURE_COMPARISON_METHODS)
    mkpath(dirname(path))
    names = String["time", "t_over_T", "alpha_deg"]
    for (method, _) in methods
        push!(names, "CL_$(method)")
        push!(names, "CD_$(method)")
    end
    for (method, _) in methods[2:end]
        push!(names, "pressure_l2_relative_$(method)")
        push!(names, "pressure_linf_$(method)")
        append!(names, ("cg_solved_$(method)", "cg_iterations_$(method)", "absolute_residual_$(method)",
                        "relative_residual_$(method)", "rhs_l2_$(method)",
                        "gradient_l2_$(method)"))
    end
    open(path, "w") do io
        println(io, join(names, ","))
        for i in eachindex(t_range)
            row = Any[t_range[i], t_range[i] / period, alpha_history[i]]
            for (method, _) in methods
                push!(row, forces[method][3, i])
                push!(row, forces[method][1, i])
            end
            for (method, _) in methods[2:end]
                push!(row, pressure_metrics[method].l2[i])
                push!(row, pressure_metrics[method].linf[i])
                d = diagnostics[method]
                append!(row, (d.converged[i], d.iterations[i], d.absolute_residual[i],
                              d.relative_residual[i], d.rhs_l2[i], d.gradient_l2[i]))
            end
            println(io, join(row, ","))
        end
    end
    return path
end

function _comparison_summary(t_range, period, forces, pressure_metrics, diagnostics;
                             skip_first_cycle::Bool=true,
                             methods=PRESSURE_COMPARISON_METHODS)
    indices = skip_first_cycle ? findall(t -> t >= period, t_range) : collect(eachindex(t_range))
    isempty(indices) && (indices = collect(eachindex(t_range)))
    rows = NamedTuple[]
    reference = forces[:bernoulli][3, indices]
    for (method, label) in methods
        valid = method == :bernoulli ? indices : [i for i in indices if diagnostics[method].converged[i]]
        cl = forces[method][3, valid]
        ref = forces[:bernoulli][3, valid]
        relative_cl = method == :bernoulli ? 0.0 :
            (isempty(valid) ? NaN : norm(cl - ref) / max(norm(ref), eps()))
        l2 = method == :bernoulli ? 0.0 :
            (isempty(valid) ? NaN : mean(pressure_metrics[method].l2[valid]))
        linf = method == :bernoulli ? 0.0 :
            (isempty(valid) ? NaN : maximum(pressure_metrics[method].linf[valid]))
        cl_mean = isempty(valid) ? NaN : mean(cl)
        cl_pp = isempty(valid) ? NaN : maximum(cl) - minimum(cl)
        push!(rows, (; method, label, converged_samples=length(valid), total_samples=length(indices),
                     CL_mean=cl_mean, CL_peak_to_peak=cl_pp,
                     relative_CL_L2=relative_cl, pressure_L2_mean=l2,
                     pressure_Linf_max=linf))
    end
    return rows
end

function _write_comparison_summary(path, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "method,label,converged_samples,total_samples,CL_mean,CL_peak_to_peak,relative_CL_L2,pressure_L2_mean,pressure_Linf_max")
        for row in rows
            println(io, join((row.method, row.label, row.converged_samples, row.total_samples,
                              row.CL_mean, row.CL_peak_to_peak,
                              row.relative_CL_L2, row.pressure_L2_mean,
                              row.pressure_Linf_max), ","))
        end
    end
    return path
end

function _plot_comparison(path, t_range, period, forces, pressure_metrics;
                          methods=PRESSURE_COMPARISON_METHODS,
                          stem="pitching_wing_pressure")
    plt = Core.eval(@__MODULE__, :(import PythonPlot as pressure_comparison_plt; pressure_comparison_plt))
    return Base.invokelatest(_plot_comparison_impl, plt, path, collect(t_range), period,
        forces, pressure_metrics, methods, stem)
end

function _plot_comparison_impl(plt, path, t_range, period, forces, pressure_metrics,
                               methods, stem)
    mkpath(path)
    cycles = collect(t_range) ./ period
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for (method, label) in methods
        ax.plot(cycles, forces[method][3, :]; linewidth=1.25, label=label)
    end
    ax.set_xlabel("t/T"); ax.set_ylabel("C_L"); ax.grid(true, alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); load_path = joinpath(path, "$(stem)_loads.png")
    fig.savefig(load_path, dpi=170); plt.pyplot.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for (method, label) in methods[2:end]
        ax.plot(cycles, pressure_metrics[method].l2; linewidth=1.25, label=label)
    end
    ax.set_xlabel("t/T"); ax.set_ylabel("relative pressure L2 error")
    ax.set_yscale("log"); ax.grid(true, alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); error_path = joinpath(path, "$(stem)_errors.png")
    fig.savefig(error_path, dpi=170); plt.pyplot.close(fig)
    return (; loads=load_path, errors=error_path)
end

function run_pitching_wing_pressure_comparison(; path=get(ENV, "PRESSURE_COMPARISON_PATH",
                                                   joinpath("data", "pitching_wing_pressure_comparison")),
        plot=_env_bool("PRESSURE_COMPARISON_PLOT", true), save_vtk=_env_bool("SAVE_VTK", false),
        n_cycles=_env_float("N_CYCLES", 3.0), c_per_dt=_env_float("C_PER_DT", 0.5),
        n_span=_env_int("N_SPAN", 13), n_airfoil=_env_int("N_AIRFOIL", 161),
        n_endcap=_env_int("N_ENDCAP", 9), panel_wake_rows=nothing,
        pressure_itmax=nothing,
        backend=pnl.FastMultipoleBackend(expansion_order=_env_int("FMM_EXPANSION_ORDER", 8),
            multipole_acceptance=_env_float("FMM_ACCEPTANCE", 0.4),
            leaf_size=_env_int("FMM_LEAF_SIZE", 40)))
    setup0 = prepare_pitching_wing(; include_static_polar=false, save_vtk,
        path, n_cycles, c_per_dt, n_span, n_airfoil, n_endcap, panel_wake_rows, backend)
    nt = length(setup0.t_range)
    normalization = pnl.WingNormalization(setup0.setup.rho, setup0.setup.Sref, setup0.setup.c)
    pressure_itmax === nothing && (pressure_itmax = max(1000,
        ceil(Int, _env_float("PRESSURE_ITMAX_PER_PANEL", 2.0) * setup0.wing.ncells)))
    bernoulli = pnl.PressureBernoulli(setup0.setup.rho;
        unsteady=true, allow_partial=true, backend)
    laplace = Dict(method => pnl.PressureLaplace((setup0.wing,), setup0.setup.rho;
        unsteady=true, gradient_mode=method, reference_panel=1, verbose=false,
        itmax=pressure_itmax)
        for (method, _) in PRESSURE_COMPARISON_METHODS[2:end])
    pressures = merge((bernoulli=bernoulli,), laplace)
    forces = Dict{Symbol,Any}()
    recorders = Dict{Symbol,PressureComparisonRecorder}()
    reference_history = Vector{Vector{Float64}}()
    monitors = Any[]
    for (method, _) in PRESSURE_COMPARISON_METHODS
        pressure = pressures[method]
        force = pnl.ForceMonitor(nt, 1; normalization, i_frame=-1,
            correct_kuttacondition=true, verbose=false, file=false)
        recorder = PressureComparisonRecorder(method, nt; reference_history,
            laplace=method == :bernoulli ? nothing : pressure)
        forces[method] = force.force
        recorders[method] = recorder
        push!(monitors, pressure); push!(monitors, force); push!(monitors, recorder)
    end

    println("Pitching-wing pressure comparison")
    @printf("  Julia threads = %d, OMP_NUM_THREADS = %s\n", Threads.nthreads(), get(ENV, "OMP_NUM_THREADS", "unset"))
    @printf("  panels = %d, steps = %d, dt = %.6g s\n", setup0.wing.ncells, nt, setup0.setup.dt)
    @printf("  PressureLaplace itmax = %d\n", pressure_itmax)
    pnl.simulate!((setup0.wing,), (setup0.wake,), setup0.frames, setup0.maneuver!, setup0.Uinf,
        setup0.t_range; body_solvers=(setup0.solver,), backend, monitors=Tuple(monitors),
        path=save_vtk ? path : nothing, name="pitching_wing_pressure_comparison",
        set_Das_eta_freestream=NaN, verbose=false)

    pressure_metrics = Dict(method => (; l2=recorders[method].pressure_l2_relative,
        linf=recorders[method].pressure_linf) for (method, _) in PRESSURE_COMPARISON_METHODS[2:end])
    diagnostics = Dict(method => (; converged=recorders[method].converged,
        iterations=recorders[method].iterations,
        absolute_residual=recorders[method].absolute_residual,
        relative_residual=recorders[method].relative_residual,
        rhs_l2=recorders[method].rhs_l2,
        gradient_l2=recorders[method].gradient_l2)
        for (method, _) in PRESSURE_COMPARISON_METHODS[2:end])
    alpha_history = [setup0.setup.alpha_mean_deg + setup0.setup.alpha_amp_deg *
        sin(setup0.setup.omega * t) for t in setup0.t_range]
    csv_path = _comparison_csv(joinpath(path, "pitching_wing_pressure_comparison.csv"),
        setup0.t_range, setup0.setup.period, alpha_history, forces, pressure_metrics, diagnostics)
    rows = _comparison_summary(setup0.t_range, setup0.setup.period, forces, pressure_metrics, diagnostics)
    summary_path = _write_comparison_summary(joinpath(path, "pitching_wing_pressure_comparison_summary.csv"), rows)
    plot_paths = plot ? _plot_comparison(path, setup0.t_range, setup0.setup.period, forces, pressure_metrics) : nothing

    println("\nmethod                         solved      CL mean       CL p-p       rel CL L2       p L2 mean       p Linf max")
    for row in rows
        @printf("%-30s %4d/%-4d % .6e  % .6e  % .6e  % .6e  % .6e\n", row.label,
            row.converged_samples, row.total_samples, row.CL_mean, row.CL_peak_to_peak, row.relative_CL_L2,
            row.pressure_L2_mean, row.pressure_Linf_max)
    end
    println("\nWrote comparison CSV: $(csv_path)")
    println("Wrote summary CSV: $(summary_path)")
    plot !== false && println("Wrote plots: $(plot_paths.loads), $(plot_paths.errors)")
    return (; setup=setup0, pressures, forces, pressure_metrics, diagnostics, rows, csv_path, summary_path, plot_paths)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_pitching_wing_pressure_comparison()
end
