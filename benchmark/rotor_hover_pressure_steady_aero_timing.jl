#!/usr/bin/env julia

# Benchmark the substeps of `_steady_aerodynamics!` over two full post-restart
# rotor-hover timesteps. The second step is treated as the steady-state truth;
# the first step pays one-time FMM allocation / branch warm-up costs.
#
# This is a near-copy of rotor_hover_pressure_restart_timing.jl with the single
# `steady_aerodynamics_total` timer replaced by a line-by-line inlined version
# of `_steady_aerodynamics!` (src/FLOWPanel_simulate.jl:281) where each top-
# level block has its own timer.

import FLOWPanel as pnl
using Printf
using Dates

const ORIGINAL_RUN_NAME = "rotor_hover_pressure_comparison"
const BENCH_NAME = "rotor_hover_pressure_steady_aero_timing"
const RESULTS_DIR = joinpath("benchmark", "results")
const PLOT_DIR = RESULTS_DIR
const VTK_COMPRESS = parse(Bool, get(ENV, "VTK_COMPRESS", "false"))

function default_restart_step()
    pvd = joinpath("data", ORIGINAL_RUN_NAME, ORIGINAL_RUN_NAME * "_body1.pvd")
    _, idxs = pnl._read_pvd_steps(pvd)
    isempty(idxs) && error("No saved restart steps found in $(pvd)")
    return idxs[end]
end

const RESTART_STEP_DEFAULT = parse(Int, get(ENV, "RESTART_STEP", string(default_restart_step())))
const MEASURED_STEPS = (RESTART_STEP_DEFAULT + 1):(RESTART_STEP_DEFAULT + 2)
const FIRST_STEP = first(MEASURED_STEPS)
const LAST_STEP = last(MEASURED_STEPS)

const SA_SUBSTEPS = (
    "sa_reset",
    "sa_collect",
    "sa_freestream",
    "sa_kinematic_velocity",
    "sa_update_te",
    "sa_wake_influence",
    "sa_set_core_size_panel",
    "sa_solve",
    "sa_add_bound_vorticity",
    "sa_set_core_size_targets",
    "sa_body_influence",
    "sa_diagnose_particle_influence",
    "sa_mu_half_jump",
)

function set_setup_env!()
    ENV["RHPC_SETUP_ONLY"] = "1"
    ENV["RUN_NAME"] = ORIGINAL_RUN_NAME
    ENV["RUN_MONITORS"] = "true"
    ENV["RUN_KJ"] = "false"
    ENV["SAVE_VTK"] = "true"
    ENV["RHPC_MESH"] = get(ENV, "RHPC_MESH", "40_40")
    # Need a time vector long enough for the final measured step to include
    # propagation/shed.
    ENV["NREVS"] = get(ENV, "NREVS", "10.11111111111111")
    return nothing
end

function canonical_setup!()
    set_setup_env!()
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "examples", "rotor_hover_pressure_comparison.jl"))
    m = @__MODULE__
    g(sym) = Core.eval(m, sym)
    return (
        systems=g(:systems),
        wakes=g(:wakes),
        frames=g(:frames),
        maneuver=g(:maneuver!),
        Uinf=g(:Uinf),
        t_range=g(:t_range),
        monitors=g(:monitors),
        body_solvers=g(:body_solvers),
        backend=g(:backend),
        backend_wake=g(:backend_wake),
        wakerow_no_hessian_to_particles=g(:wakerow_no_hessian_to_particles),
        body_hessian_to_particles=g(:body_hessian_to_particles),
        body_on_wake=g(:body_on_wake),
        panel_wake_on_particles=g(:panel_wake_on_particles),
        particle_hessian_self=g(:particle_hessian_self),
        particle_relax=g(:particle_relax),
        diagnose_particle_gamma=g(:diagnose_particle_gamma),
        diagnose_particle_influence=g(:diagnose_particle_influence),
        diagnostic_vertical=g(:particle_diagnostic_vertical),
        set_Das_min_kinematic_displacement=g(:set_Das_min_kinematic_displacement),
    )
end

elapsed_s(f) = begin
    t0 = time_ns()
    f()
    (time_ns() - t0) / 1e9
end

function add_time!(times, phase, i_step, dt)
    get!(times, phase, Dict{Int, Float64}())
    times[phase][i_step] = get(times[phase], i_step, 0.0) + dt
    return nothing
end

function particle_count(wakes)
    total = 0
    for w in wakes
        if w isa pnl.PanelParticleWake
            total += w.pfield.np
        end
    end
    return total
end

function monitor_label(monitor, i)
    return @sprintf("monitor_%02d_%s", i, nameof(typeof(monitor)))
end

function load_restart!(cfg; restart_step=RESTART_STEP_DEFAULT)
    systems_tuple = pnl._systems_tuple(cfg.systems)
    wakes_tuple = pnl._wakes_tuple(cfg.systems, cfg.wakes)
    restart_path = joinpath("data", ORIGINAL_RUN_NAME)

    for sys in systems_tuple
        pnl.calc_normals!(sys)
        pnl.calc_controlpoints!(sys)
    end

    for (i, sys) in enumerate(systems_tuple)
        pnl._load_body_vtk!(sys, restart_path, ORIGINAL_RUN_NAME * "_body$(i)", restart_step)
    end
    for (i, w) in enumerate(wakes_tuple)
        if w isa pnl.PanelParticleWake
            pnl._load_panel_particle_wake_vtk!(w, restart_path, ORIGINAL_RUN_NAME * "_wake$(i)", restart_step)
        elseif w isa pnl.PanelWake
            pnl._load_panel_wake_vtk!(w, restart_path, ORIGINAL_RUN_NAME * "_wake$(i)", restart_step)
        elseif !isnothing(w)
            error("Unsupported wake type $(typeof(w))")
        end
    end
    pnl._load_frame_state_toml!(cfg.frames, restart_path, ORIGINAL_RUN_NAME, restart_step)

    for sys in systems_tuple
        pnl.calc_normals!(sys)
        pnl.calc_controlpoints!(sys)
    end

    return systems_tuple, wakes_tuple
end

function replay_restart_end!(cfg, systems_tuple, wakes_tuple; restart_step=RESTART_STEP_DEFAULT)
    times = Dict{String, Dict{Int, Float64}}()
    dt_end = cfg.t_range[restart_step + 2] - cfg.t_range[restart_step + 1]

    add_time!(times, "setup_replay_propagate", restart_step,
        elapsed_s() do
            for w in wakes_tuple
                if w isa pnl.PanelParticleWake
                    pnl.propagate!(w, dt_end; relax=cfg.particle_relax, step=restart_step,
                        frames=cfg.frames, diagnose_particle_gamma=cfg.diagnose_particle_gamma,
                        diagnostic_vertical=cfg.diagnostic_vertical)
                elseif !isnothing(w)
                    pnl.propagate!(w, dt_end; step=restart_step, frames=cfg.frames)
                end
            end
        end)

    add_time!(times, "setup_replay_kinematics", restart_step,
        elapsed_s() do
            pnl.propagate_kinematics!(systems_tuple, cfg.frames, dt_end)
            for sys in systems_tuple
                pnl.calc_normals!(sys)
                pnl.calc_controlpoints!(sys)
            end
        end)

    add_time!(times, "setup_replay_shed", restart_step,
        elapsed_s() do
            for (sys, w) in zip(systems_tuple, wakes_tuple)
                !isnothing(w) && pnl.shed_wake!(w, sys)
            end
        end)

    return times
end

function write_outputs!(cfg, systems_tuple, wakes_tuple, uinf, i_step, t, path, name, first_output_step)
    times = Dict{String, Float64}()
    overwrite = i_step == first_output_step

    if overwrite || !isfile(pnl._metadata_toml_path(path, name))
        times["output_metadata_manifest"] = elapsed_s() do
            pnl._write_metadata_toml(path, name, systems_tuple, wakes_tuple, cfg.frames, cfg.t_range,
                cfg.body_solvers, cfg.backend_wake, cfg.backend, cfg.backend, cfg.monitors;
                start_step=first_output_step,
                set_Das_eta_kinematic=NaN,
                set_Das_eta_freestream=NaN,
                set_Das_min_kinematic_displacement=cfg.set_Das_min_kinematic_displacement,
                clean_files=true)
        end
    end

    times["output_body_vtk"] = elapsed_s() do
        for (i, sys) in enumerate(systems_tuple)
            body_name = name * "_body$(i)"
            pnl.write_vtk(joinpath(path, body_name), sys, i_step, t;
                monitors=cfg.monitors, i_system=i, overwrite, compress=VTK_COMPRESS)
        end
    end

    times["output_wake_vtk"] = elapsed_s() do
        for (i, w) in enumerate(wakes_tuple)
            if !isnothing(w)
                wake_name = name * "_wake$(i)"
                pnl.write_vtk(joinpath(path, wake_name), w, i_step, t; overwrite, compress=VTK_COMPRESS)
            end
        end
    end

    times["output_metadata_append"] = elapsed_s() do
        pnl._append_metadata_step_toml(path, name, cfg.frames, i_step, t; uinf)
    end

    times["output_total"] = sum(values(times))
    return times
end

# Inlined, instrumented copy of `_steady_aerodynamics!`
# (src/FLOWPanel_simulate.jl:281). The control flow is preserved verbatim;
# each top-level block is wrapped in its own timer. Returns the sum of
# substep elapsed times so the caller can record a synthetic
# `steady_aerodynamics_total`.
function timed_steady_aerodynamics!(times, i_step, cfg, systems_tuple, wakes_tuple, uinf;
        needs_induced_vorticity::Bool, update_trailing_edges::Bool=true)

    backend_wake = cfg.backend_wake
    backend_solve = cfg.backend
    backend_system = cfg.backend
    body_hessian_to_particles      = cfg.body_hessian_to_particles
    body_on_wake                   = cfg.body_on_wake
    panel_wake_on_particles        = cfg.panel_wake_on_particles
    particle_hessian_self          = cfg.particle_hessian_self
    wakerow_no_hessian_to_particles = cfg.wakerow_no_hessian_to_particles
    diagnose_particle_influence    = cfg.diagnose_particle_influence
    diagnostic_vertical            = cfg.diagnostic_vertical
    FV_PF = pnl.FLOWVPM.ParticleField

    total = 0.0

    dt_phase = elapsed_s() do
        for w in wakes_tuple
            !isnothing(w) && pnl.reset!(w)
        end
        for sys in systems_tuple
            pnl.reset!(sys)
        end
    end
    add_time!(times, "sa_reset", i_step, dt_phase); total += dt_phase

    local wake_probes, targets, wake_sources
    dt_phase = elapsed_s() do
        wake_probes = pnl._collect_wake_probes(wakes_tuple)
        targets = (systems_tuple..., wake_probes...)
        wake_sources = pnl._collect_wake_sources(wakes_tuple)
    end
    add_time!(times, "sa_collect", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        pnl.apply_freestream!(systems_tuple, uinf)
        for w in wakes_tuple
            !isnothing(w) && pnl.apply_freestream!(w, uinf)
        end
    end
    add_time!(times, "sa_freestream", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        pnl.kinematic_velocity!(systems_tuple, cfg.frames)
    end
    add_time!(times, "sa_kinematic_velocity", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        if update_trailing_edges
            for (sys, w) in zip(systems_tuple, wakes_tuple)
                !isnothing(w) && pnl.update_TE!(w, sys)
            end
        end
    end
    add_time!(times, "sa_update_te", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        if length(wake_sources) > 0
            if wakerow_no_hessian_to_particles || !panel_wake_on_particles
                panel_sources  = Tuple(s for s in wake_sources if !(s isa FV_PF))
                pfield_sources = Tuple(s for s in wake_sources if   s isa FV_PF)
                panel_targets = !panel_wake_on_particles ?
                    Tuple(t for t in targets if !(t isa FV_PF)) : targets
                if length(panel_sources) > 0 && length(panel_targets) > 0
                    pnl.influence!(panel_targets, panel_sources, backend_wake; precalc=true,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(sys isa FV_PF ? false : pnl.requires_hessian(sys) for sys in panel_targets),
                        extra_outputs=pnl._induced_vorticity_extra_outputs(panel_targets, needs_induced_vorticity))
                end
                if length(pfield_sources) > 0
                    particle_targets = Tuple(t for t in targets if t isa FV_PF)
                    nonparticle_targets = Tuple(t for t in targets if !(t isa FV_PF))
                    if particle_hessian_self || length(particle_targets) == 0
                        pnl.influence!(targets, pfield_sources, backend_wake; precalc=true,
                            postcalc=true,
                            scalar_potential=false,
                            velocity=true,
                            velocity_gradient=Tuple(pnl.requires_hessian(sys) for sys in targets),
                            extra_outputs=pnl._induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
                    else
                        if length(nonparticle_targets) > 0
                            pnl.influence!(nonparticle_targets, pfield_sources, backend_wake; precalc=true,
                                postcalc=false,
                                scalar_potential=false,
                                velocity=true,
                                velocity_gradient=Tuple(pnl.requires_hessian(sys) for sys in nonparticle_targets),
                                extra_outputs=pnl._induced_vorticity_extra_outputs(nonparticle_targets, needs_induced_vorticity))
                        end
                        pnl.influence!(particle_targets, pfield_sources, backend_wake; precalc=true,
                            postcalc=true,
                            scalar_potential=false,
                            velocity=true,
                            velocity_gradient=Tuple(false for _ in particle_targets),
                            extra_outputs=pnl._induced_vorticity_extra_outputs(particle_targets, needs_induced_vorticity))
                    end
                end
            elseif !particle_hessian_self && any(source isa FV_PF for source in wake_sources)
                panel_sources = Tuple(s for s in wake_sources if !(s isa FV_PF))
                pfield_sources = Tuple(s for s in wake_sources if s isa FV_PF)
                if length(panel_sources) > 0
                    pnl.influence!(targets, panel_sources, backend_wake; precalc=true,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(pnl.requires_hessian(sys) for sys in targets),
                        extra_outputs=pnl._induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
                end
                particle_targets = Tuple(t for t in targets if t isa FV_PF)
                nonparticle_targets = Tuple(t for t in targets if !(t isa FV_PF))
                if length(pfield_sources) > 0
                    if length(nonparticle_targets) > 0
                        pnl.influence!(nonparticle_targets, pfield_sources, backend_wake; precalc=true,
                            postcalc=false,
                            scalar_potential=false,
                            velocity=true,
                            velocity_gradient=Tuple(pnl.requires_hessian(sys) for sys in nonparticle_targets),
                            extra_outputs=pnl._induced_vorticity_extra_outputs(nonparticle_targets, needs_induced_vorticity))
                    end
                    if length(particle_targets) > 0
                        pnl.influence!(particle_targets, pfield_sources, backend_wake; precalc=true,
                            postcalc=true,
                            scalar_potential=false,
                            velocity=true,
                            velocity_gradient=Tuple(false for _ in particle_targets),
                            extra_outputs=pnl._induced_vorticity_extra_outputs(particle_targets, needs_induced_vorticity))
                    end
                end
            else
                wake_postcalc = any(source isa FV_PF for source in wake_sources)
                pnl.influence!(targets, wake_sources, backend_wake; precalc=true,
                    postcalc=wake_postcalc,
                    scalar_potential=false,
                    velocity=true,
                    velocity_gradient=Tuple(pnl.requires_hessian(sys) for sys in targets),
                    extra_outputs=pnl._induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
            end
        end
    end
    add_time!(times, "sa_wake_influence", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        pnl._set_core_sizes!(systems_tuple, :core_size_panel)
    end
    add_time!(times, "sa_set_core_size_panel", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        pnl.solve!(cfg.systems, cfg.body_solvers; backend=backend_solve)
    end
    add_time!(times, "sa_solve", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        needs_induced_vorticity && pnl._add_bound_surface_vorticity!(systems_tuple)
    end
    add_time!(times, "sa_add_bound_vorticity", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        pnl._set_core_sizes!(systems_tuple, :core_size_targets)
    end
    add_time!(times, "sa_set_core_size_targets", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        if !body_on_wake
            pnl.influence!(systems_tuple, systems_tuple, backend_system; precalc=false,
                scalar_potential=false,
                velocity=true,
                velocity_gradient=Tuple(pnl.requires_hessian(sys) for sys in systems_tuple),
                extra_outputs=pnl._induced_vorticity_extra_outputs(systems_tuple, needs_induced_vorticity),
                direct_conditioning=pnl._self_panel_core_size_conditioning())
        else
            pnl.influence!(targets, systems_tuple, backend_system; precalc=false,
                scalar_potential=false,
                velocity=true,
                velocity_gradient=Tuple((sys isa FV_PF && !body_hessian_to_particles) ? false : pnl.requires_hessian(sys) for sys in targets),
                extra_outputs=pnl._induced_vorticity_extra_outputs(targets, needs_induced_vorticity),
                direct_conditioning=pnl._self_panel_core_size_conditioning())
        end
    end
    add_time!(times, "sa_body_influence", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        if diagnose_particle_influence
            pnl._diagnose_particle_influence!(wakes_tuple, systems_tuple, backend_wake, backend_system;
                needs_induced_vorticity, particle_hessian_self, diagnostic_vertical)
        end
    end
    add_time!(times, "sa_diagnose_particle_influence", i_step, dt_phase); total += dt_phase

    dt_phase = elapsed_s() do
        for body in systems_tuple
            if pnl.has_grad_mu(body)
                pnl.compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
                    body.cells, body.neighbor,
                    view(body.strength, :, pnl.get_Gammai(body)),
                    pnl._bound_surface_vorticity_te_info(body);
                    scale=0.5,
                    nodes=body.nodes)
            end
        end
    end
    add_time!(times, "sa_mu_half_jump", i_step, dt_phase); total += dt_phase

    add_time!(times, "steady_aerodynamics_total", i_step, total)
    return total
end

function run_two_steps!(cfg; path=nothing, name=BENCH_NAME)
    systems_tuple, wakes_tuple = load_restart!(cfg)
    setup_times = replay_restart_end!(cfg, systems_tuple, wakes_tuple)

    pnl.audit_monitors(cfg.monitors)
    needs_grad = any(pnl.monitor_requires_body_hessian, cfg.monitors)
    needs_induced_vorticity = any(pnl.monitor_requires_induced_vorticity, cfg.monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = needs_grad
    end

    if !isnothing(path)
        mkpath(path)
        pnl._clean_monitor_csv_files!(path, name)
    end

    times = Dict{String, Dict{Int, Float64}}()
    totals = Dict{Int, Float64}()
    np_before = particle_count(wakes_tuple)
    np_by_step_before = Dict{Int, Int}()
    np_by_step_after = Dict{Int, Int}()

    for i_step in MEASURED_STEPS
        t = cfg.t_range[i_step + 1]
        dt = cfg.t_range[i_step + 2] - cfg.t_range[i_step + 1]
        np_by_step_before[i_step] = particle_count(wakes_tuple)
        step_total_start = time_ns()

        add_time!(times, "controls", i_step,
            elapsed_s() do
                Base.invokelatest(cfg.maneuver, cfg.frames, systems_tuple, wakes_tuple, t)
            end)

        uinf = Base.invokelatest(cfg.Uinf, t)
        timed_steady_aerodynamics!(times, i_step, cfg, systems_tuple, wakes_tuple, uinf;
            needs_induced_vorticity, update_trailing_edges=true)

        monitor_context = pnl.MonitorContext()
        pnl.monitor_set_time!(monitor_context, t)
        monitor_csv_dir = isnothing(path) ? nothing : joinpath(path, "monitors")
        monitor_compute_total = 0.0
        monitor_csv_total = 0.0
        for (i_monitor, monitor) in enumerate(cfg.monitors)
            label = monitor_label(monitor, i_monitor)
            dt_monitor = elapsed_s() do
                pnl._run_monitor!(monitor, monitor_context, systems_tuple, wakes_tuple,
                    cfg.frames, uinf, i_step, dt, t)
            end
            monitor_compute_total += dt_monitor
            add_time!(times, label, i_step, dt_monitor)
            if !isnothing(monitor_csv_dir)
                dt_csv = elapsed_s() do
                    pnl.write_monitor_csv!(monitor, monitor_csv_dir, name, i_monitor,
                        monitor_context, systems_tuple, i_step, dt;
                        overwrite=i_step == FIRST_STEP)
                end
                monitor_csv_total += dt_csv
                add_time!(times, "csv_" * label, i_step, dt_csv)
            end
        end
        add_time!(times, "monitors_compute_total", i_step, monitor_compute_total)
        !isnothing(path) && add_time!(times, "monitors_csv_total", i_step, monitor_csv_total)

        if !isnothing(path)
            output_times = write_outputs!(cfg, systems_tuple, wakes_tuple, uinf, i_step, t,
                path, name, FIRST_STEP)
            for (phase, dt_phase) in output_times
                add_time!(times, phase, i_step, dt_phase)
            end
        end

        add_time!(times, "propagate_total", i_step,
            elapsed_s() do
                for w in wakes_tuple
                    if w isa pnl.PanelParticleWake
                        pnl.propagate!(w, dt; relax=cfg.particle_relax, step=i_step,
                            frames=cfg.frames, diagnose_particle_gamma=cfg.diagnose_particle_gamma,
                            diagnostic_vertical=cfg.diagnostic_vertical)
                    elseif !isnothing(w)
                        pnl.propagate!(w, dt; step=i_step, frames=cfg.frames)
                    end
                end
                pnl.propagate_kinematics!(systems_tuple, cfg.frames, dt)
                for sys in systems_tuple
                    pnl.calc_normals!(sys)
                    pnl.calc_controlpoints!(sys)
                end
            end)

        add_time!(times, "shed_total", i_step,
            elapsed_s() do
                for (sys, w) in zip(systems_tuple, wakes_tuple)
                    !isnothing(w) && pnl.shed_wake!(w, sys)
                end
            end)

        np_by_step_after[i_step] = particle_count(wakes_tuple)
        totals[i_step] = (time_ns() - step_total_start) / 1e9
        @printf("  %s step %d: %.3f s, particles %d -> %d\n",
            isnothing(path) ? "compute" : "output ", i_step, totals[i_step],
            np_by_step_before[i_step], np_by_step_after[i_step])
        flush(stdout)
    end

    return (
        times=times,
        totals=totals,
        setup_times=setup_times,
        np_before=np_before,
        np_after=particle_count(wakes_tuple),
    )
end

function phase_value(times, phase, step)
    return get(get(times, phase, Dict{Int, Float64}()), step, 0.0)
end

function rows_for_thread(threads, compute_result, output_result)
    rows = Vector{NamedTuple}()
    major = vcat(
        ["controls"],
        collect(SA_SUBSTEPS),
        ["steady_aerodynamics_total",
         "monitors_compute_total",
         "monitors_csv_total",
         "output_total",
         "propagate_total",
         "shed_total"],
    )
    detail_phases = sort!(collect(union(keys(compute_result.times), keys(output_result.times))))
    phases = unique(vcat(major, detail_phases))

    # `steady_aerodynamics_total` is already the sum of sa_* substeps, so it
    # must not be re-added when computing the grand total per step.
    total_phases = filter(p -> p != "steady_aerodynamics_total", major)
    total_steps = Dict(step => 0.0 for step in MEASURED_STEPS)
    for phase in total_phases
        source = startswith(phase, "output") || startswith(phase, "monitors_csv") ? output_result : compute_result
        for step in MEASURED_STEPS
            total_steps[step] += phase_value(source.times, phase, step)
        end
    end
    grand_total = sum(values(total_steps))

    for phase in phases
        source = startswith(phase, "output") || startswith(phase, "csv_") || phase == "monitors_csv_total" ? output_result : compute_result
        vals = [phase_value(source.times, phase, step) for step in MEASURED_STEPS]
        total = sum(vals)
        push!(rows, (
            threads=threads, phase=phase,
            step_first_s=vals[1], step_last_s=vals[2],
            total_s=total, mean_s=total / length(MEASURED_STEPS),
            percent_of_total=grand_total > 0 ? 100 * total / grand_total : 0.0,
            np_before=compute_result.np_before, np_after=compute_result.np_after,
        ))
    end

    diff_vals = [output_result.totals[step] - compute_result.totals[step] for step in MEASURED_STEPS]
    push!(rows, (
        threads=threads, phase="output_by_total_difference",
        step_first_s=diff_vals[1], step_last_s=diff_vals[2],
        total_s=sum(diff_vals), mean_s=sum(diff_vals) / length(diff_vals),
        percent_of_total=grand_total > 0 ? 100 * sum(diff_vals) / grand_total : 0.0,
        np_before=compute_result.np_before, np_after=compute_result.np_after,
    ))

    sort!(rows, by=r -> r.total_s, rev=true)
    push!(rows, (
        threads=threads, phase="TOTAL",
        step_first_s=total_steps[FIRST_STEP], step_last_s=total_steps[LAST_STEP],
        total_s=grand_total, mean_s=grand_total / length(MEASURED_STEPS),
        percent_of_total=100.0,
        np_before=compute_result.np_before, np_after=compute_result.np_after,
    ))
    return rows
end

function write_csv(path, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "threads,phase,step_$(FIRST_STEP)_s,step_$(LAST_STEP)_s,total_s,mean_s,percent_of_total,np_before,np_after")
        for r in rows
            println(io, join((
                r.threads, r.phase,
                r.step_first_s, r.step_last_s,
                r.total_s, r.mean_s, r.percent_of_total,
                r.np_before, r.np_after), ","))
        end
    end
    return path
end

function read_csv_rows(path)
    lines = readlines(path)
    rows = Vector{NamedTuple}()
    for line in lines[2:end]
        isempty(strip(line)) && continue
        cols = split(line, ",")
        push!(rows, (
            threads=parse(Int, cols[1]), phase=String(cols[2]),
            step_first_s=parse(Float64, cols[3]),
            step_last_s=parse(Float64, cols[4]),
            total_s=parse(Float64, cols[5]),
            mean_s=parse(Float64, cols[6]),
            percent_of_total=parse(Float64, cols[7]),
            np_before=parse(Int, cols[8]),
            np_after=parse(Int, cols[9]),
        ))
    end
    return rows
end

function markdown_table(rows)
    buf = IOBuffer()
    println(buf, "| threads | phase | step_$(FIRST_STEP)_s | step_$(LAST_STEP)_s | total_s | mean_s | percent_of_total | np_before | np_after |")
    println(buf, "|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows
        @printf(buf, "| %d | %s | %.3f | %.3f | %.3f | %.3f | %.1f | %d | %d |\n",
            r.threads, r.phase, r.step_first_s, r.step_last_s,
            r.total_s, r.mean_s, r.percent_of_total, r.np_before, r.np_after)
    end
    return String(take!(buf))
end

function maybe_write_combined_artifacts()
    paths = [joinpath(RESULTS_DIR, "$(BENCH_NAME)_threads$(n).csv") for n in (1, 3)]
    all(isfile, paths) || return nothing
    rows = vcat(read_csv_rows.(paths)...)
    sort!(rows, by=r -> (r.threads, r.phase == "TOTAL" ? Inf : -r.total_s))

    md_path = joinpath(RESULTS_DIR, "$(BENCH_NAME)_summary.md")
    open(md_path, "w") do io
        println(io, "# Rotor Hover Pressure Steady-Aero Substep Timing")
        println(io)
        println(io, "Generated: $(Dates.now())")
        println(io)
        println(io, "Step $(LAST_STEP) (`step_$(LAST_STEP)_s`) is the steady-state measurement; ")
        println(io, "step $(FIRST_STEP) (`step_$(FIRST_STEP)_s`) includes one-time FMM warm-up costs.")
        println(io)
        print(io, markdown_table(rows))
    end
    println("Wrote combined Markdown summary: $(md_path)")

    try
        @eval using PythonPlot
        pp = Core.eval(@__MODULE__, :PythonPlot)

        phases_of_interest = vcat(collect(SA_SUBSTEPS),
            ["propagate_total", "shed_total",
             "monitors_compute_total", "monitors_csv_total", "output_total",
             "output_by_total_difference"])
        thread_counts = [1, 3]

        means_by_thread = Dict(n => Dict{String,Float64}() for n in thread_counts)
        for r in rows
            r.phase in phases_of_interest || continue
            means_by_thread[r.threads][r.phase] = r.mean_s
        end

        out = joinpath(PLOT_DIR, "$(BENCH_NAME)_stacked.png")
        Base.invokelatest(pp, phases_of_interest, thread_counts, means_by_thread, out) do pp, phases, tcs, means, outpath
            cmap = pp.get_cmap("tab20")
            colors = [cmap(i / length(phases)) for i in 0:length(phases)-1]
            fig, ax = pp.subplots(; figsize=(8, 6))
            x = 1:length(tcs)
            bottoms = zeros(length(tcs))
            for (pi, phase) in enumerate(phases)
                heights = [get(means[n], phase, 0.0) for n in tcs]
                ax.bar(x, heights; bottom=bottoms, color=colors[pi], label=phase)
                bottoms .+= heights
            end
            ax.set_xticks(collect(x))
            ax.set_xticklabels(["threads=$(n)" for n in tcs])
            ax.set_ylabel("Mean time per step (s)")
            ax.set_title("Rotor hover steady-aero substep timing")
            ax.legend(loc="upper right", fontsize=6)
            fig.tight_layout()
            fig.savefig(outpath, dpi=180)
            println("Wrote plot: $(outpath)")
        end
    catch err
        @warn "Could not generate PythonPlot bar plots" exception=(err, catch_backtrace())
    end
    return md_path
end

function main()
    threads = Threads.nthreads()
    output_path = joinpath("/private/tmp", "flowpanel_rotor_hover_steady_aero_bench_t$(threads)")
    output_name = BENCH_NAME * "_threads$(threads)"

    println("Rotor hover steady-aero substep timing benchmark")
    println("  threads: $(threads)")
    println("  measured steps: $(FIRST_STEP):$(LAST_STEP) (step $(LAST_STEP) is truth)")
    println("  output path: $(output_path)")

    if parse(Bool, get(ENV, "BENCH_PREFLIGHT", "false"))
        println("\nPreflight only: setup + restart load, no measured timesteps.")
        cfg = canonical_setup!()
        systems_tuple, wakes_tuple = load_restart!(cfg)
        println("  loaded bodies: $(length(systems_tuple))")
        println("  body nodes/cells: " *
            join(["$(size(sys.nodes, 2))/$(sys.ncells)" for sys in systems_tuple], ", "))
        println("  particle count: $(particle_count(wakes_tuple))")
        return nothing
    end

    if parse(Bool, get(ENV, "BENCH_PLOT_ONLY", "false"))
        println("\nPlot-only mode: reading existing CSVs and regenerating artifacts.")
        maybe_write_combined_artifacts()
        return nothing
    end

    println("\nSetting up compute-only restart...")
    compute_cfg = canonical_setup!()
    compute_result = run_two_steps!(compute_cfg; path=nothing, name=output_name)

    println("\nSetting up output restart...")
    output_cfg = canonical_setup!()
    output_result = run_two_steps!(output_cfg; path=output_path, name=output_name)

    rows = rows_for_thread(threads, compute_result, output_result)
    csv_path = joinpath(RESULTS_DIR, "$(BENCH_NAME)_threads$(threads).csv")
    write_csv(csv_path, rows)
    println("\nWrote per-run CSV: $(csv_path)")
    println(markdown_table(rows))
    maybe_write_combined_artifacts()
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isempty(PROGRAM_FILE)
    main()
end
