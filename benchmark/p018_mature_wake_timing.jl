#!/usr/bin/env julia

# BRAINSTORM/023: per-phase timing of production-config 018 timesteps restarted
# from a MATURE wake (~181k particles at rev ~28.5 of p018_cs_f1_l3p4).
#
# Successor to rotor_hover_pressure_steady_aero_timing.jl (2026-06, 40_40 mesh,
# 39k particles) with two structural fixes:
#   1. The warm-start restore mirrors simulate_warmstart! sections 2.5-5
#      (src/FLOWPanel_warmstart.jl), including the MANDATORY kinematic replay
#      that re-rotates the un-persisted `Das` shed offsets and the
#      `_restore_wake_continuation!` / `_kutta_warmstart_restore!` calls. The
#      old harness loaded frames from TOML only, which misplaces the first
#      wake row (see the section 2.5 comment in FLOWPanel_warmstart.jl).
#   2. _steady_aerodynamics! is no longer inlined: the stage helpers
#      (_sa_collect, _sa_reset_freestream_kinematic!, _sa_wake_influence!,
#      solve_formulation!, _sa_body_influence!, _sa_half_jump!) now exist as
#      standalone functions in src/FLOWPanel_simulate.jl, so each is timed
#      directly and the control flow cannot drift from the source.
#
# Environment (beyond the full production RHPC env exported by the wrapper
# benchmark/slurm/p023_mature_timing.sh):
#   BENCH_RESTART_RUN   restart run name under data/ (default p018_cs_f1_l3p4)
#   RESTART_STEP        step to restart from (-1 = highest in the body PVD)
#   BENCH_NSTEPS        measured steps (default 5; step 1 reported separately
#                       as JIT/FMM warm-up — judge steady phases on steps 2+)
#   BENCH_WRITE_OUTPUT  also time VTK/metadata/monitor-CSV writes (default true;
#                       writes go to BENCH_OUT_PATH, never to data/)
#   BENCH_OUT_PATH      output-arm scratch dir (default under tempdir())
#   BENCH_PROFILE       run one extra step under Profile.@profile and write a
#                       flat report (default false)
#
# Legacy-Kutta production configs only (RHPC passes no wake_attachment /
# kutta_closure). Asserted below.

import FLOWPanel as pnl
using Printf
using Profile

const RESTART_RUN = get(ENV, "BENCH_RESTART_RUN", "p018_cs_f1_l3p4")
const BENCH_NAME = get(ENV, "BENCH_NAME", "p018_mature_wake_timing")
const RESULTS_DIR = joinpath("benchmark", "results")
const NSTEPS = parse(Int, get(ENV, "BENCH_NSTEPS", "5"))
const WRITE_OUTPUT = parse(Bool, get(ENV, "BENCH_WRITE_OUTPUT", "true"))
const OUT_PATH = get(ENV, "BENCH_OUT_PATH",
    joinpath(tempdir(), "flowpanel_p023_bench_t$(Threads.nthreads())"))
const DO_PROFILE = parse(Bool, get(ENV, "BENCH_PROFILE", "false"))
const VTK_COMPRESS = parse(Bool, get(ENV, "VTK_COMPRESS", "true"))

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

function particle_count(wakes_tuple)
    total = 0
    for w in wakes_tuple
        if w isa pnl.PanelParticleWake
            total += w.pfield.np
        end
    end
    return total
end

function canonical_setup!()
    ENV["RHPC_SETUP_ONLY"] = "1"
    ENV["RUN_NAME"] = RESTART_RUN
    Base.include(@__MODULE__, joinpath(@__DIR__, "..", "examples",
        "rotor_hover_pressure_comparison.jl"))
    m = @__MODULE__
    g(sym) = Core.eval(m, sym)
    # BENCH_FMM_BODY / BENCH_FMM_WAKE ("p,mac,leaf") override the RHPC
    # backends independently — the driver's FMM_* env vars set both backends
    # at once, which cannot express the per-pass tuned knobs from
    # benchmark/p023_fmm_tune.jl (RHPC itself is frozen).
    parse_backend(s) = begin
        p, mac, leaf = split(s, r"[:,]")
        pnl.FastMultipoleBackend(parse(Int, p), parse(Float64, mac), parse(Int, leaf))
    end
    backend_body = haskey(ENV, "BENCH_FMM_BODY") ?
        parse_backend(ENV["BENCH_FMM_BODY"]) : g(:backend)
    backend_wake_ = haskey(ENV, "BENCH_FMM_WAKE") ?
        parse_backend(ENV["BENCH_FMM_WAKE"]) : g(:backend_wake)
    haskey(ENV, "BENCH_FMM_BODY") && println("  backend override (body): $(backend_body)")
    haskey(ENV, "BENCH_FMM_WAKE") && println("  backend override (wake): $(backend_wake_)")
    return (
        systems=g(:systems),
        wakes=g(:wakes),
        frames=g(:frames),
        maneuver=g(Symbol("maneuver!")),
        Uinf=g(:Uinf),
        t_range=g(:t_range),
        monitors=g(:monitors),
        body_solvers=g(:body_solvers),
        backend=backend_body,
        backend_wake=backend_wake_,
        formulation=g(:formulation),
        wakerow_no_hessian_to_particles=g(:wakerow_no_hessian_to_particles),
        body_hessian_to_particles=g(:body_hessian_to_particles),
        body_gradient_core_size=g(:body_gradient_core_size),
        body_on_wake=g(:body_on_wake),
        panel_wake_on_particles=g(:panel_wake_on_particles),
        particle_hessian_self=g(:particle_hessian_self),
        particle_relax=g(:particle_relax),
        bound_strength_rlx=g(:bound_strength_rlx),
        diagnose_particle_gamma=g(:diagnose_particle_gamma),
        diagnose_particle_influence=g(:diagnose_particle_influence),
        diagnostic_vertical=g(:particle_diagnostic_vertical),
        set_Das_refresh=g(:set_Das_refresh),
        set_Das_min_kinematic_displacement=g(:set_Das_min_kinematic_displacement),
    )
end

function default_restart_step()
    pvd = joinpath("data", RESTART_RUN, RESTART_RUN * "_body1.pvd")
    _, idxs = pnl._read_pvd_steps(pvd)
    isempty(idxs) && error("No saved restart steps found in $(pvd)")
    return idxs[end]
end

# Mirror of simulate_warmstart! sections 2-5 (src/FLOWPanel_warmstart.jl:363-497)
# for the legacy-Kutta, set_Das_refresh=false production configuration, with
# setup timers. Any behavioral deviation from simulate_warmstart! is a bug.
function warmstart_restore!(cfg, restart_step)
    @assert !cfg.set_Das_refresh "production 018 configs run set_Das_refresh=false"
    systems_tuple = pnl._systems_tuple(cfg.systems)
    wakes_tuple = pnl._wakes_tuple(cfg.systems, cfg.wakes)
    frames = cfg.frames
    t_range = cfg.t_range
    rpath = joinpath("data", RESTART_RUN)
    rname = RESTART_RUN
    setup_times = Dict{String, Float64}()

    @assert restart_step + NSTEPS + 1 < length(t_range) "t_range too short: " *
        "restart_step=$(restart_step) + $(NSTEPS) measured steps needs " *
        ">= $(restart_step + NSTEPS + 2) entries, have $(length(t_range)) " *
        "(check SETTLE_REVS matches the parent run)"

    # section 2: pre-loop geometry init (Das eta re-accumulation skipped —
    # RHPC passes set_Das_eta_kinematic=NaN when set_Das_refresh=false)
    for sys in systems_tuple
        pnl.calc_normals!(sys)
        pnl.calc_controlpoints!(sys)
    end

    # section 2.5: kinematic replay 0..restart_step-1 (rotates Das), then this
    # step's maneuver!; frame manifest is authoritative for the frames.
    setup_times["setup_kinematic_replay"] = elapsed_s() do
        for i in 0:(restart_step - 1)
            Base.invokelatest(cfg.maneuver, frames, systems_tuple, wakes_tuple, t_range[i+1])
            pnl.propagate_kinematics!(systems_tuple, frames, t_range[i+2] - t_range[i+1])
        end
        Base.invokelatest(cfg.maneuver, frames, systems_tuple, wakes_tuple, t_range[restart_step+1])
    end
    try
        frames_manifest = deepcopy(frames)
        pnl._load_frame_state_toml!(frames_manifest, rpath, rname, restart_step)
        for (fr, fm) in zip(frames, frames_manifest)
            dev = max(maximum(abs, fr.R .- fm.R), maximum(abs, fr.x .- fm.x))
            dev > 1e-8 && @warn "replayed frame state deviates from manifest (max dev $(dev))"
        end
        for i in eachindex(frames)
            frames[i] = frames_manifest[i]
        end
    catch err
        println("frame manifest unavailable ($(sprint(showerror, err))); using replayed frames")
    end

    # section 3: on-disk state at restart_step + wake continuation state
    setup_times["setup_load_state"] = elapsed_s() do
        for (i, sys) in enumerate(systems_tuple)
            pnl._load_body_vtk!(sys, rpath, rname * "_body$(i)", restart_step)
        end
        for (i, w) in enumerate(wakes_tuple)
            isnothing(w) && continue
            wake_name = rname * "_wake$(i)"
            if w isa pnl.PanelParticleWake
                pnl._load_panel_particle_wake_vtk!(w, rpath, wake_name, restart_step)
            elseif w isa pnl.PanelWake
                pnl._load_panel_wake_vtk!(w, rpath, wake_name, restart_step)
            else
                error("unsupported wake type $(typeof(w))")
            end
        end
    end
    metadata = pnl._read_metadata_toml(rpath, rname)
    pnl._restore_wake_continuation!(wakes_tuple, metadata, restart_step)

    # section 4/4.5: derived geometry + legacy Kutta restore
    for sys in systems_tuple
        pnl.calc_normals!(sys)
        pnl.calc_controlpoints!(sys)
    end
    pnl._kutta_warmstart_restore!(systems_tuple, wakes_tuple, metadata,
        restart_step, pnl.RigidTransitionAttachment(), pnl.JumpKutta())

    # section 5: replay the saved step's skipped end-of-step actions
    # (propagate! deliberately WITHOUT the relax kwarg — mirror of
    # FLOWPanel_warmstart.jl:478-488)
    dt_end = t_range[restart_step + 2] - t_range[restart_step + 1]
    setup_times["setup_replay_end_of_step"] = elapsed_s() do
        for w in wakes_tuple
            !isnothing(w) && pnl.propagate!(w, dt_end; step=restart_step, frames)
        end
        pnl.propagate_kinematics!(systems_tuple, frames, dt_end)
        for sys in systems_tuple
            pnl.calc_normals!(sys)
            pnl.calc_controlpoints!(sys)
        end
        for (sys, w) in zip(systems_tuple, wakes_tuple)
            !isnothing(w) && pnl.shed_wake!(w, sys)
        end
    end

    return systems_tuple, wakes_tuple, setup_times
end

# One full production timestep, phase-timed. Mirror of the simulate! loop body
# (src/FLOWPanel_simulate.jl:1181-1366) on the legacy-Kutta path, with
# _steady_aerodynamics! (:654) expanded into its stage helpers.
function timed_step!(times, cfg, systems_tuple, wakes_tuple, formulation_state,
        i_step; needs_induced_vorticity, path, name, first_output_step)
    t = cfg.t_range[i_step + 1]
    dt = i_step < length(cfg.t_range) - 1 ?
        cfg.t_range[i_step+2] - cfg.t_range[i_step+1] :
        cfg.t_range[i_step+1] - cfg.t_range[i_step]
    step_start = time_ns()

    add_time!(times, "controls", i_step, elapsed_s() do
        Base.invokelatest(cfg.maneuver, cfg.frames, systems_tuple, wakes_tuple, t)
    end)
    uinf = Base.invokelatest(cfg.Uinf, t)

    # --- _steady_aerodynamics!, stage by stage ---
    normalized_grad_mu_options = pnl._normalize_grad_mu_options((;); default_basis=:quad)
    local targets, wake_sources
    add_time!(times, "sa_collect", i_step, elapsed_s() do
        _, targets, wake_sources = pnl._sa_collect(systems_tuple, wakes_tuple)
    end)
    add_time!(times, "sa_reset_freestream_kinematic", i_step, elapsed_s() do
        pnl._sa_reset_freestream_kinematic!(systems_tuple, wakes_tuple, cfg.frames, uinf)
    end)
    add_time!(times, "sa_update_te", i_step, elapsed_s() do
        for (sys, w) in zip(systems_tuple, wakes_tuple)
            !isnothing(w) && pnl.update_TE!(w, sys)
        end
    end)
    add_time!(times, "sa_formulation_prewake", i_step, elapsed_s() do
        pnl.formulation_prewake!(cfg.formulation, formulation_state, systems_tuple)
    end)
    add_time!(times, "sa_wake_influence", i_step, elapsed_s() do
        pnl._sa_wake_influence!(targets, wake_sources, cfg.backend_wake;
            needs_induced_vorticity,
            wakerow_no_hessian_to_particles=cfg.wakerow_no_hessian_to_particles,
            panel_wake_on_particles=cfg.panel_wake_on_particles,
            particle_hessian_self=cfg.particle_hessian_self)
    end)
    add_time!(times, "sa_set_core_size_panel", i_step, elapsed_s() do
        pnl._set_core_sizes!(systems_tuple, :core_size_panel)
    end)
    add_time!(times, "sa_solve", i_step, elapsed_s() do
        pnl.solve_formulation!(cfg.formulation, formulation_state, cfg.systems,
            systems_tuple, wakes_tuple, cfg.body_solvers;
            backend_solve=cfg.backend, backend_wake=cfg.backend_wake, i_step)
    end)
    add_time!(times, "sa_add_bound_vorticity", i_step, elapsed_s() do
        needs_induced_vorticity && pnl._add_bound_surface_vorticity!(systems_tuple;
            grad_mu_options=normalized_grad_mu_options)
    end)
    add_time!(times, "sa_set_core_size_targets", i_step, elapsed_s() do
        pnl._set_core_sizes!(systems_tuple, :core_size_targets)
    end)
    add_time!(times, "sa_body_influence", i_step, elapsed_s() do
        pnl._sa_body_influence!(targets, systems_tuple, cfg.backend;
            needs_induced_vorticity,
            body_on_wake=cfg.body_on_wake,
            body_hessian_to_particles=cfg.body_hessian_to_particles,
            body_gradient_core_size=cfg.body_gradient_core_size)
    end)
    add_time!(times, "sa_diagnose_particle_influence", i_step, elapsed_s() do
        cfg.diagnose_particle_influence &&
            pnl._diagnose_particle_influence!(wakes_tuple, systems_tuple,
                cfg.backend_wake, cfg.backend;
                needs_induced_vorticity,
                particle_hessian_self=cfg.particle_hessian_self,
                diagnostic_vertical=cfg.diagnostic_vertical)
    end)
    add_time!(times, "sa_half_jump", i_step, elapsed_s() do
        pnl._sa_half_jump!(systems_tuple, normalized_grad_mu_options)
    end)
    @assert cfg.bound_strength_rlx == 1 "bound_strength_rlx != 1 not carried in this harness"

    # --- monitors ---
    monitor_context = pnl.MonitorContext()
    pnl.monitor_set_time!(monitor_context, t)
    monitor_csv_dir = isnothing(path) ? nothing : joinpath(path, "monitors")
    monitor_compute = 0.0
    monitor_csv = 0.0
    for (i_monitor, monitor) in enumerate(cfg.monitors)
        label = @sprintf("monitor_%02d_%s", i_monitor, nameof(typeof(monitor)))
        dt_m = elapsed_s() do
            pnl._run_monitor!(monitor, monitor_context, systems_tuple, wakes_tuple,
                cfg.frames, uinf, i_step, dt, t)
        end
        monitor_compute += dt_m
        add_time!(times, label, i_step, dt_m)
        if !isnothing(monitor_csv_dir)
            monitor_csv += elapsed_s() do
                pnl.write_monitor_csv!(monitor, monitor_csv_dir, name, i_monitor,
                    monitor_context, systems_tuple, i_step, dt;
                    overwrite=i_step == first_output_step)
            end
        end
    end
    add_time!(times, "monitors_compute_total", i_step, monitor_compute)
    add_time!(times, "monitors_csv_total", i_step, monitor_csv)

    # --- save state ---
    if !isnothing(path)
        add_time!(times, "output_metadata", i_step, elapsed_s() do
            if i_step == first_output_step || !isfile(pnl._metadata_toml_path(path, name))
                pnl._write_metadata_toml(path, name, systems_tuple, wakes_tuple,
                    cfg.frames, cfg.t_range, cfg.body_solvers, cfg.backend_wake,
                    cfg.backend, cfg.backend, cfg.monitors;
                    start_step=first_output_step,
                    set_Das_eta_kinematic=NaN, set_Das_eta_freestream=NaN,
                    set_Das_min_kinematic_displacement=cfg.set_Das_min_kinematic_displacement,
                    clean_files=true)
            end
            pnl._append_metadata_step_toml(path, name, cfg.frames, i_step, t;
                uinf, wakes=wakes_tuple)
        end)
        add_time!(times, "output_body_vtk", i_step, elapsed_s() do
            for (i, sys) in enumerate(systems_tuple)
                pnl.write_vtk(joinpath(path, name * "_body$(i)"), sys, i_step, t;
                    monitors=cfg.monitors, i_system=i,
                    overwrite=i_step == first_output_step, compress=VTK_COMPRESS)
            end
        end)
        add_time!(times, "output_wake_vtk", i_step, elapsed_s() do
            for (i, w) in enumerate(wakes_tuple)
                !isnothing(w) && pnl.write_vtk(joinpath(path, name * "_wake$(i)"),
                    w, i_step, t; overwrite=i_step == first_output_step,
                    compress=VTK_COMPRESS)
            end
        end)
    end

    # --- propagate + shed ---
    add_time!(times, "propagate_total", i_step, elapsed_s() do
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
    add_time!(times, "shed_total", i_step, elapsed_s() do
        for (sys, w) in zip(systems_tuple, wakes_tuple)
            !isnothing(w) && pnl.shed_wake!(w, sys)
        end
    end)

    return (time_ns() - step_start) / 1e9
end

function write_results(times, totals, setup_times, np_by_step, measured_steps)
    mkpath(RESULTS_DIR)
    threads = Threads.nthreads()
    csv_path = joinpath(RESULTS_DIR, "$(BENCH_NAME)_t$(threads).csv")
    phases = sort!(collect(keys(times)))
    steady_steps = length(measured_steps) > 1 ? measured_steps[2:end] : measured_steps
    open(csv_path, "w") do io
        println(io, "phase," * join(("step_$(s)_s" for s in measured_steps), ",") *
            ",steady_mean_s,steady_pct")
        steady_total = sum(totals[s] for s in steady_steps)
        for phase in phases
            vals = [get(times[phase], s, 0.0) for s in measured_steps]
            steady = sum(get(times[phase], s, 0.0) for s in steady_steps) / length(steady_steps)
            pct = 100 * steady * length(steady_steps) / steady_total
            println(io, phase * "," * join((@sprintf("%.4f", v) for v in vals), ",") *
                @sprintf(",%.4f,%.2f", steady, pct))
        end
        println(io, "STEP_TOTAL," * join((@sprintf("%.4f", totals[s]) for s in measured_steps), ",") *
            @sprintf(",%.4f,100.00", steady_total / length(steady_steps)))
        println(io, "N_PARTICLES," * join((string(np_by_step[s]) for s in measured_steps), ",") * ",,")
        for (k, v) in sort!(collect(setup_times))
            println(io, k * "," * @sprintf("%.2f", v) * repeat(",", length(measured_steps)) * ",")
        end
    end
    println("\nWrote $(csv_path)")

    # ranked steady-state table (steps 2+; step 1 pays JIT/FMM warm-up)
    println("\nSteady-state phase ranking (mean over steps $(steady_steps)):")
    ranked = sort!([(phase, sum(get(times[phase], s, 0.0) for s in steady_steps) / length(steady_steps))
                    for phase in phases], by=x -> -x[2])
    steady_mean_total = sum(totals[s] for s in steady_steps) / length(steady_steps)
    println("| phase | mean_s | % of step |")
    println("|---|---:|---:|")
    for (phase, v) in ranked
        v < 0.005 && continue
        @printf("| %s | %.3f | %.1f |\n", phase, v, 100 * v / steady_mean_total)
    end
    @printf("| STEP_TOTAL | %.3f | 100.0 |\n", steady_mean_total)
    return csv_path
end

function main()
    threads = Threads.nthreads()
    println("BRAINSTORM/023 mature-wake per-phase timing")
    println("  restart run: $(RESTART_RUN)  threads: $(threads)  nsteps: $(NSTEPS)")
    println("  write output: $(WRITE_OUTPUT) -> $(OUT_PATH)")
    flush(stdout)

    setup_s = @elapsed cfg = canonical_setup!()
    println("  RHPC setup (incl. solver factorization): $(round(setup_s, digits=1)) s")
    flush(stdout)

    restart_step = parse(Int, get(ENV, "RESTART_STEP", "-1"))
    restart_step < 0 && (restart_step = Base.invokelatest(default_restart_step))
    println("  restarting from step $(restart_step)")

    systems_tuple, wakes_tuple, setup_times =
        Base.invokelatest(warmstart_restore!, cfg, restart_step)
    println("  restored particle count: $(particle_count(wakes_tuple))")
    for (k, v) in sort!(collect(setup_times))
        println("    $(k): $(round(v, digits=1)) s")
    end
    flush(stdout)

    pnl.audit_monitors(cfg.monitors)
    needs_grad = any(pnl.monitor_requires_body_hessian, cfg.monitors)
    needs_induced_vorticity = any(pnl.monitor_requires_induced_vorticity, cfg.monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = needs_grad
    end
    println("  needs_body_hessian: $(needs_grad)  needs_induced_vorticity: $(needs_induced_vorticity)")

    formulation_state = Base.invokelatest(pnl.initialize_formulation, cfg.formulation,
        systems_tuple, wakes_tuple, cfg.body_solvers, cfg.backend, cfg.backend)

    path = WRITE_OUTPUT ? OUT_PATH : nothing
    name = BENCH_NAME
    if !isnothing(path)
        mkpath(path)
        pnl._clean_monitor_csv_files!(path, name)
    end

    measured_steps = (restart_step + 1):(restart_step + NSTEPS)
    times = Dict{String, Dict{Int, Float64}}()
    totals = Dict{Int, Float64}()
    np_by_step = Dict{Int, Int}()
    for i_step in measured_steps
        np_by_step[i_step] = particle_count(wakes_tuple)
        totals[i_step] = Base.invokelatest(timed_step!, times, cfg, systems_tuple,
            wakes_tuple, formulation_state, i_step;
            needs_induced_vorticity, path, name, first_output_step=first(measured_steps))
        phase_sum = sum(get(d, i_step, 0.0) for d in values(times)) -
            get(get(times, "monitors_compute_total", Dict{Int,Float64}()), i_step, 0.0) -
            get(get(times, "monitors_csv_total", Dict{Int,Float64}()), i_step, 0.0)
        @printf("  step %d: %.1f s total, %.1f s in timed phases, %d particles\n",
            i_step, totals[i_step], phase_sum, np_by_step[i_step])
        flush(stdout)
    end

    write_results(times, totals, setup_times, np_by_step, collect(measured_steps))

    if DO_PROFILE
        i_step = last(measured_steps) + 1
        println("\nProfiling step $(i_step)...")
        Profile.init(n=10^8, delay=0.005)
        profile_times = Dict{String, Dict{Int, Float64}}()
        Profile.@profile Base.invokelatest(timed_step!, profile_times, cfg,
            systems_tuple, wakes_tuple, formulation_state, i_step;
            needs_induced_vorticity, path, name, first_output_step=first(measured_steps))
        prof_path = joinpath(RESULTS_DIR, "$(BENCH_NAME)_t$(threads)_profile_flat.txt")
        open(prof_path, "w") do io
            Profile.print(IOContext(io, :displaysize => (200, 240));
                format=:flat, sortedby=:count, mincount=500, combine=true)
        end
        println("Wrote $(prof_path)")
        # tree report: preserves the call hierarchy so nearfield direct work
        # vs farfield expansion passes inside FastMultipole.fmm! can be split
        # (the flat report loses the parents of inlined leaf arithmetic)
        tree_path = joinpath(RESULTS_DIR, "$(BENCH_NAME)_t$(threads)_profile_tree.txt")
        open(tree_path, "w") do io
            Profile.print(IOContext(io, :displaysize => (200, 240));
                format=:tree, mincount=2000, maxdepth=40, combine=true)
        end
        println("Wrote $(tree_path)")
    end
    println("\nDone.")
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isempty(PROGRAM_FILE)
    main()
end
