#--- Das initialization helpers ---#

function _accumulate_Das!(sys::AbstractLiftingBody, eta; min_displacement=0.0)
    for ishedding in eachindex(sys.Das)
        Das = sys.Das[ishedding]
        Vte = sys.velocity_te[ishedding]
        for j in axes(Das, 2)
            speed = sqrt(Vte[1, j]^2 + Vte[2, j]^2 + Vte[3, j]^2)
            if speed > zero(speed)
                displacement_length = max(abs(eta) * speed, min_displacement)
                signed_length = eta < zero(eta) ? -displacement_length : displacement_length
                scale = signed_length / speed
                Das[1, j] += Vte[1, j] * scale
                Das[2, j] += Vte[2, j] * scale
                Das[3, j] += Vte[3, j] * scale
            end
        end
    end
end

_accumulate_Das!(::AbstractBody, eta; min_displacement=0.0) = nothing

#--- wake tuple helpers ---#

function _systems_tuple(systems::Tuple)
    for system in systems
        system isa AbstractBody || throw(ArgumentError("Each system must be an AbstractBody, got $(typeof(system))"))
    end
    return systems
end
_systems_tuple(system::AbstractBody) = (system,)
_systems_tuple(systems) = throw(ArgumentError("systems must be an AbstractBody or a Tuple of AbstractBody objects, got $(typeof(systems))"))

function _wakes_tuple(systems::Tuple, wakes)
    wakes isa Tuple || throw(ArgumentError("wakes must be a Tuple when systems is a Tuple, got $(typeof(wakes))"))
    length(wakes) == length(systems) || throw(ArgumentError("Number of wakes ($(length(wakes))) must match number of systems ($(length(systems)))"))
    for wake in wakes
        (isnothing(wake) || wake isa AbstractFreeWake) || throw(ArgumentError("Each wake must be an AbstractFreeWake or nothing, got $(typeof(wake))"))
    end
    return wakes
end

function _wakes_tuple(system::AbstractBody, wake)
    (isnothing(wake) || wake isa AbstractFreeWake) || throw(ArgumentError("wake must be an AbstractFreeWake or nothing when systems is an AbstractBody, got $(typeof(wake))"))
    return (wake,)
end

_induced_vorticity_extra_outputs(targets::Tuple, enabled::Bool) =
    Tuple(enabled && target isa AbstractBody ? 3 : 0 for target in targets)

function _bound_surface_vorticity_te_info(body::AbstractBody)
    return hasproperty(body, :shedding_full) ?
        view(body.shedding_full, 1:2, :) :
        zeros(Int, 2, body.ncells)
end

function _add_bound_surface_vorticity!(systems::Tuple)
    for body in systems
        _add_bound_surface_vorticity!(body)
    end
    return nothing
end

function _add_bound_surface_vorticity!(body::AbstractBody)
    has_grad_mu(body) || return nothing

    # Accumulate the bound surface vorticity κ = n × ∇sμ into
    # body.induced_vorticity on top of any wake-induced contribution already
    # there. Build κ in a scratch buffer so compute_mu_gradient! (which writes
    # ∇μ, not accumulates) does not clobber existing values.
    TF = eltype(body.induced_vorticity)
    kappa = zeros(TF, 3, body.ncells)
    compute_mu_gradient!(kappa, body.controlpoints, body.normals,
        body.cells, body.neighbor,
        view(body.strength, :, get_Gammai(body)),
        _bound_surface_vorticity_te_info(body);
        scale=-1.0)

    @inbounds for i in axes(kappa, 2)
        nx, ny, nz = body.normals[1, i], body.normals[2, i], body.normals[3, i]
        gx = kappa[1, i]
        gy = kappa[2, i]
        gz = kappa[3, i]
        body.induced_vorticity[1, i] += ny * gz - nz * gy
        body.induced_vorticity[2, i] += nz * gx - nx * gz
        body.induced_vorticity[3, i] += nx * gy - ny * gx
    end

    return nothing
end

function _validate_body_solvers(systems::Tuple, body_solvers)
    if body_solvers isa Tuple
        length(body_solvers) == length(systems) || throw(ArgumentError("Number of body_solvers ($(length(body_solvers))) must match number of systems ($(length(systems)))"))
    elseif !(body_solvers isa AbstractSolver)
        throw(ArgumentError("body_solvers must be a Tuple or an AbstractSolver, got $(typeof(body_solvers))"))
    end
    return body_solvers
end

function _validate_body_solvers(system::AbstractBody, body_solvers)
    if body_solvers isa Tuple
        throw(ArgumentError("body_solvers must be an AbstractSolver when systems is an AbstractBody, got $(typeof(body_solvers))"))
    elseif !(body_solvers isa AbstractSolver)
        throw(ArgumentError("body_solvers must be a Tuple or an AbstractSolver, got $(typeof(body_solvers))"))
    end
    return body_solvers
end

function _validate_influence_backend(name::Symbol, backend)
    backend isa AbstractBackend || throw(ArgumentError("$(name) must be an AbstractBackend for influence! calls, got $(typeof(backend))"))
    return backend
end

function _validate_solve_backend(systems, body_solvers, backend)
    if systems isa AbstractBody && (backend isa Tuple || backend isa AbstractVector)
        throw(ArgumentError("backend_solve must be a single backend when systems is an AbstractBody, got $(typeof(backend))"))
    elseif body_solvers isa Tuple && (backend isa Tuple || backend isa AbstractVector)
        length(backend) == length(_systems_tuple(systems)) || throw(ArgumentError("backend_solve length ($(length(backend))) must match number of systems ($(length(_systems_tuple(systems))))"))
    elseif body_solvers isa AbstractSolver && (backend isa Tuple || backend isa AbstractVector)
        throw(ArgumentError("backend_solve must be a single backend when body_solvers is an AbstractSolver, got $(typeof(backend))"))
    end
    return backend
end

function _collect_wake_probes(wakes::Tuple)
    result = ()
    for w in wakes
        if !isnothing(w)
            result = (result..., get_probes(w)...)
        end
    end
    return result
end

_collect_wake_probes(wake::AbstractFreeWake) = _collect_wake_probes((wake,))
_collect_wake_probes(::Nothing) = ()

function _set_kerneloffsets!(systems::Tuple, field::Symbol)
    for system in systems
        system.kerneloffset = getfield(system, field)
    end
    return nothing
end

function _self_panel_kerneloffset_conditioning()
    before! = function (source_buffer, source_system, i_source_system, target_buffer, i_target_system)
        source_system.kerneloffset = source_system.kerneloffset_panel
        return nothing
    end
    after! = function (source_buffer, source_system, i_source_system, target_buffer, i_target_system)
        source_system.kerneloffset = source_system.kerneloffset_targets
        return nothing
    end
    return FastMultipole.DirectConditioningRule(FastMultipole.SelfPairs(), before!, after!)
end

function _collect_wake_sources(wakes::Tuple)
    result = ()
    for w in wakes
        if !isnothing(w)
            result = (result..., get_sources(w)...)
        end
    end
    return result
end

_collect_wake_sources(wake::AbstractFreeWake) = _collect_wake_sources((wake,))
_collect_wake_sources(::Nothing) = ()

function _particle_j_rms(pfield::FLOWVPM.ParticleField)
    np = pfield.np
    np == 0 && return 0.0
    J = view(pfield.particles, FLOWVPM.J_INDEX, 1:np)
    return sqrt(sum(abs2, J) / np)
end

function _set_particle_j_zero!(pfield::FLOWVPM.ParticleField)
    pfield.np == 0 && return nothing
    view(pfield.particles, FLOWVPM.J_INDEX, 1:pfield.np) .= 0
    return nothing
end

function _diagnostic_particle_j_from_sources!(pfield::FLOWVPM.ParticleField, sources::Tuple,
        backend; velocity_gradient::Bool=true, optargs...)
    length(sources) == 0 && return 0.0
    pfield.np == 0 && return 0.0

    particles = copy(view(pfield.particles, :, 1:pfield.np))
    try
        _set_particle_j_zero!(pfield)
        influence!((pfield,), sources, backend;
            scalar_potential=false,
            velocity=false,
            velocity_gradient=(velocity_gradient,),
            precalc=false,
            postcalc=false,
            optargs...)
        return _particle_j_rms(pfield)
    finally
        view(pfield.particles, :, 1:size(particles, 2)) .= particles
    end
end

function _diagnose_particle_influence!(wakes_tuple::Tuple, systems_tuple::Tuple,
        backend_wake, backend_system; needs_induced_vorticity::Bool=false,
        particle_hessian_self::Bool=true, diagnostic_vertical=(0.0, 0.0, 1.0))
    for (iw, w) in enumerate(wakes_tuple)
        w isa PanelParticleWake || continue
        pfield = w.pfield
        pfield.np == 0 && continue

        wake_sources = _collect_wake_sources((w,))
        panel_sources = Tuple(s for s in wake_sources if !(s isa FLOWVPM.ParticleField))
        particle_sources = Tuple(s for s in wake_sources if s isa FLOWVPM.ParticleField)

        j_total = _particle_j_rms(pfield)
        j_body = _diagnostic_particle_j_from_sources!(pfield, systems_tuple, backend_system;
            direct_conditioning=_self_panel_kerneloffset_conditioning())
        j_panelwake = _diagnostic_particle_j_from_sources!(pfield, panel_sources, backend_wake)
        j_particles = _diagnostic_particle_j_from_sources!(pfield, particle_sources, backend_wake;
            velocity_gradient=particle_hessian_self)

        println("particle influence wake=$(iw) np=$(pfield.np) |J_total|=$(j_total) |J_body|=$(j_body) |J_panelwake|=$(j_panelwake) |J_particles|=$(j_particles)")
    end
    return nothing
end

# Wake sources that have a well-defined scalar potential. Excludes vortex
# particle fields (which carry only a vector potential).
_scalar_potential_sources(w::PanelWake) = get_sources(w)
_scalar_potential_sources(w::PanelParticleWake) = get_sources(w.panel_wake)
_scalar_potential_sources(::Nothing) = ()

function _collect_wake_scalar_sources(wakes::Tuple)
    result = ()
    for w in wakes
        result = (result..., _scalar_potential_sources(w)...)
    end
    return result
end

function initialize_Das!(systems, frames, Uinf::Function, t0, dt0;
        set_Das_eta_kinematic=NaN,
        set_Das_eta_freestream=NaN,
        set_Das_min_kinematic_displacement=0.0)
    if isnan(set_Das_eta_freestream) && isnan(set_Das_eta_kinematic)
        return systems
    end

    systems_tuple = _systems_tuple(systems)

    if !isnan(set_Das_eta_freestream)
        uinf0 = Uinf(t0)
        for sys in systems_tuple
            extra_reset!(sys)
            extra_apply_freestream!(sys, uinf0)
            _accumulate_Das!(sys, dt0 * set_Das_eta_freestream)
        end
    end

    if !isnan(set_Das_eta_kinematic)
        for sys in systems_tuple
            extra_reset!(sys)
        end
        kinematic_velocity!(systems_tuple, frames)
        for sys in systems_tuple
            _accumulate_Das!(sys, dt0 * set_Das_eta_kinematic;
                min_displacement=set_Das_min_kinematic_displacement)
        end
    end

    # reset velocity fields modified during Das computation
    for sys in systems_tuple
        reset!(sys)
    end

    return systems
end

function _steady_aerodynamics!(systems, systems_tuple::Tuple, wakes_tuple::Tuple,
        frames, uinf, body_solvers; backend_wake=nothing, backend_solve,
        backend_system, needs_induced_vorticity::Bool=false,
        update_trailing_edges::Bool=false,
        wakerow_no_hessian_to_particles::Bool=false,
        body_hessian_to_particles::Bool=false,
        body_on_wake::Bool=true,
        panel_wake_on_particles::Bool=true,
        particle_hessian_self::Bool=true,
        diagnose_particle_influence::Bool=false,
        diagnostic_vertical=(0.0, 0.0, 1.0))
    for w in wakes_tuple
        !isnothing(w) && reset!(w)
    end
    for sys in systems_tuple
        reset!(sys)
    end

    wake_probes = _collect_wake_probes(wakes_tuple)
    targets = (systems_tuple..., wake_probes...)
    wake_sources = _collect_wake_sources(wakes_tuple)

    apply_freestream!(systems_tuple, uinf)
    for w in wakes_tuple
        !isnothing(w) && apply_freestream!(w, uinf)
    end

    kinematic_velocity!(systems_tuple, frames)

    if update_trailing_edges
        for (sys, w) in zip(systems_tuple, wakes_tuple)
            !isnothing(w) && update_TE!(w, sys)
        end
    end

    if length(wake_sources) > 0
        # Diagnostic gates:
        #   wakerow_no_hessian_to_particles: ablate the panel-wake-row ->
        #     particle Hessian (vortex rings + closing filaments produce
        #     |∇U| ~ 1/r^2 at edges that linear panels would smooth into a sheet).
        #   panel_wake_on_particles=false: drop panel-wake-row -> particle
        #     velocity entirely; particles still receive particle-on-particle
        #     velocity, and panel wake still acts on bodies and panel-wake nodes.
        if wakerow_no_hessian_to_particles || !panel_wake_on_particles
            panel_sources  = Tuple(s for s in wake_sources if !(s isa FLOWVPM.ParticleField))
            pfield_sources = Tuple(s for s in wake_sources if   s isa FLOWVPM.ParticleField)
            panel_targets = !panel_wake_on_particles ?
                Tuple(t for t in targets if !(t isa FLOWVPM.ParticleField)) : targets
            if length(panel_sources) > 0 && length(panel_targets) > 0
                influence!(panel_targets, panel_sources, backend_wake; precalc=true,
                    scalar_potential=false,
                    velocity=true,
                    velocity_gradient=Tuple(sys isa FLOWVPM.ParticleField ? false : requires_hessian(sys) for sys in panel_targets),
                    extra_outputs=_induced_vorticity_extra_outputs(panel_targets, needs_induced_vorticity))
            end
            if length(pfield_sources) > 0
                particle_targets = Tuple(t for t in targets if t isa FLOWVPM.ParticleField)
                nonparticle_targets = Tuple(t for t in targets if !(t isa FLOWVPM.ParticleField))
                if particle_hessian_self || length(particle_targets) == 0
                    influence!(targets, pfield_sources, backend_wake; precalc=true,
                        postcalc=true,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(requires_hessian(sys) for sys in targets),
                        extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
                else
                    if length(nonparticle_targets) > 0
                        influence!(nonparticle_targets, pfield_sources, backend_wake; precalc=true,
                            postcalc=false,
                            scalar_potential=false,
                            velocity=true,
                            velocity_gradient=Tuple(requires_hessian(sys) for sys in nonparticle_targets),
                            extra_outputs=_induced_vorticity_extra_outputs(nonparticle_targets, needs_induced_vorticity))
                    end
                    influence!(particle_targets, pfield_sources, backend_wake; precalc=true,
                        postcalc=true,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(false for _ in particle_targets),
                        extra_outputs=_induced_vorticity_extra_outputs(particle_targets, needs_induced_vorticity))
                end
            end
        elseif !particle_hessian_self && any(source isa FLOWVPM.ParticleField for source in wake_sources)
            panel_sources = Tuple(s for s in wake_sources if !(s isa FLOWVPM.ParticleField))
            pfield_sources = Tuple(s for s in wake_sources if s isa FLOWVPM.ParticleField)
            if length(panel_sources) > 0
                influence!(targets, panel_sources, backend_wake; precalc=true,
                    scalar_potential=false,
                    velocity=true,
                    velocity_gradient=Tuple(requires_hessian(sys) for sys in targets),
                    extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
            end
            particle_targets = Tuple(t for t in targets if t isa FLOWVPM.ParticleField)
            nonparticle_targets = Tuple(t for t in targets if !(t isa FLOWVPM.ParticleField))
            if length(pfield_sources) > 0
                if length(nonparticle_targets) > 0
                    influence!(nonparticle_targets, pfield_sources, backend_wake; precalc=true,
                        postcalc=false,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(requires_hessian(sys) for sys in nonparticle_targets),
                        extra_outputs=_induced_vorticity_extra_outputs(nonparticle_targets, needs_induced_vorticity))
                end
                if length(particle_targets) > 0
                    influence!(particle_targets, pfield_sources, backend_wake; precalc=true,
                        postcalc=true,
                        scalar_potential=false,
                        velocity=true,
                        velocity_gradient=Tuple(false for _ in particle_targets),
                        extra_outputs=_induced_vorticity_extra_outputs(particle_targets, needs_induced_vorticity))
                end
            end
        else
            wake_postcalc = any(source isa FLOWVPM.ParticleField for source in wake_sources)
            influence!(targets, wake_sources, backend_wake; precalc=true,
                postcalc=wake_postcalc,
                scalar_potential=false,
                velocity=true,
                velocity_gradient=Tuple(requires_hessian(sys) for sys in targets),
                extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
        end
    end

    _set_kerneloffsets!(systems_tuple, :kerneloffset_panel)
    solve!(systems, body_solvers; backend=backend_solve)

    needs_induced_vorticity && _add_bound_surface_vorticity!(systems_tuple)

    _set_kerneloffsets!(systems_tuple, :kerneloffset_targets)
    if !body_on_wake
        # body-on-body only; skip the body-on-wake-probes pass so the wake
        # never receives body-induced velocity this step.
        influence!(systems_tuple, systems_tuple, backend_system; precalc=false,
            scalar_potential=false,
            velocity=true,
            velocity_gradient=Tuple(requires_hessian(sys) for sys in systems_tuple),
            extra_outputs=_induced_vorticity_extra_outputs(systems_tuple, needs_induced_vorticity),
            direct_conditioning=_self_panel_kerneloffset_conditioning())
    else
        influence!(targets, systems_tuple, backend_system; precalc=false,
            scalar_potential=false,
            velocity=true,
            velocity_gradient=Tuple((sys isa FLOWVPM.ParticleField && !body_hessian_to_particles) ? false : requires_hessian(sys) for sys in targets),
            extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity),
            direct_conditioning=_self_panel_kerneloffset_conditioning())
    end

    if diagnose_particle_influence
        _diagnose_particle_influence!(wakes_tuple, systems_tuple, backend_wake, backend_system;
            needs_induced_vorticity, particle_hessian_self, diagnostic_vertical)
    end

    # Add the +½∇μ tangential half-jump on each surface so body.velocity is
    # the EXTERIOR surface limit (matching OLD calcfield_U!). The kernel-
    # induced velocity at the on-surface centroid is the PV (continuous
    # through the doublet sheet); the exterior limit requires this extra
    # tangential half-jump that depends on neighbor strengths, which
    # _self_limit cannot supply locally.
    for body in systems_tuple
        if has_grad_mu(body)
            compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
                body.cells, body.neighbor,
                view(body.strength, :, get_Gammai(body)),
                _bound_surface_vorticity_te_info(body);
                scale=0.5,
                nodes=body.nodes)
        end
    end

    return nothing
end

"""
    steady!(systems, frames, uinf; body_solvers, backend=FastMultipoleBackend(...),
            backend_solve=backend, backend_system=backend, monitors=(), i_run=1,
            dt=1.0, path=nothing, name="default_steady", verbose=false)

Run a one-shot steady panel solve for one body or a tuple of bodies. The solve
uses wake geometry already stored on the bodies, if any, and does not create,
propagate, or shed external free wakes.
"""
function steady!(systems, frames, uinf;
        name="default_steady",
        path=nothing,
        body_solvers,
        backend=FastMultipoleBackend(;
                expansion_order=10,
                multipole_acceptance=0.4,
                leaf_size=100,
            ),
        backend_solve=backend,
        backend_system=backend,
        monitors=(),
        i_run::Int=1,
        dt::Real=1.0,
        wakerow_no_hessian_to_particles::Bool=false,
        body_hessian_to_particles::Bool=false,
        body_on_wake::Bool=true,
        panel_wake_on_particles::Bool=true,
        particle_hessian_self::Bool=true,
        diagnose_particle_influence::Bool=false,
        diagnostic_vertical=(0.0, 0.0, 1.0),
        verbose=false
    )
    i_run >= 1 || throw(ArgumentError("i_run must be >= 1, got $(i_run)."))
    dt > 0 || throw(ArgumentError("dt must be positive, got $(dt)."))

    systems_tuple = _systems_tuple(systems)
    wakes_tuple = Tuple(nothing for _ in systems_tuple)
    _validate_body_solvers(systems, body_solvers)
    _validate_influence_backend(:backend_system, backend_system)
    _validate_solve_backend(systems, body_solvers, backend_solve)
    audit_monitors(monitors)

    i_step = i_run - 1
    verbose && println("\tsteady run $(i_run)")

    needs_grad = any(monitor_requires_body_hessian, monitors)
    needs_induced_vorticity = any(monitor_requires_induced_vorticity, monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = needs_grad
    end

    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end

    for sys in systems_tuple
        calc_normals!(sys)
        calc_controlpoints!(sys)
    end

    _steady_aerodynamics!(systems, systems_tuple, wakes_tuple, frames, uinf,
        body_solvers; backend_solve, backend_system, needs_induced_vorticity,
        wakerow_no_hessian_to_particles,
        body_hessian_to_particles,
        body_on_wake,
        panel_wake_on_particles,
        particle_hessian_self,
        diagnose_particle_influence,
        diagnostic_vertical)

    monitor_context = MonitorContext()
    monitor_set_time!(monitor_context, i_step * dt)
    for monitor in monitors
        _run_monitor!(monitor, monitor_context, systems_tuple, wakes_tuple, frames, uinf, i_step, dt)
    end

    if !isnothing(path)
        t = i_step * dt
        for (i, sys) in enumerate(systems_tuple)
            body_name = name * "_body$(i)"
            write_vtk(joinpath(path, body_name), sys, i_step, t;
                      monitors=monitors, i_system=i, overwrite=i_run==1)
        end
    end

    return systems
end

"""
    simulate!(systems, wakes, frames, maneuver!, Uinf, t_range; body_solvers, backend=FastMultipoleBackend(...), backend_wake=backend, backend_solve=backend, backend_system=backend, monitors=(), ...)

Advance one or more coupled body-wake systems through `t_range`, solving the
aerodynamics, updating wakes, optionally writing VTK output, and calling any
registered monitors.
"""
function simulate!(systems, wakes, frames, maneuver!::Function, Uinf::Function, t_range;
        name="default_sim", path=joinpath("data", "default_simulation"),
        body_solvers,
        backend=FastMultipoleBackend(;
                expansion_order=10,
                multipole_acceptance=0.4,
                leaf_size=100,
            ),
        backend_wake=backend,
        backend_solve=backend,
        backend_system=backend,
        monitors=(),
        set_Das_eta_kinematic=NaN,
        set_Das_eta_freestream=NaN,
        set_Das_min_kinematic_displacement=0.0,
        start_step::Int=0,
        wakerow_no_hessian_to_particles::Bool=false,
        body_hessian_to_particles::Bool=false,
        body_on_wake::Bool=true,
        panel_wake_on_particles::Bool=true,
        particle_hessian_self::Bool=true,
        particle_relax::Bool=true,
        diagnose_particle_gamma::Bool=false,
        diagnose_particle_influence::Bool=false,
        diagnostic_vertical=(0.0, 0.0, 1.0),
        verbose=false
    )
    @assert 0 <= start_step < length(t_range) "start_step ($(start_step)) must be in [0, $(length(t_range))-1)"
    systems_tuple = _systems_tuple(systems)
    wakes_tuple = _wakes_tuple(systems, wakes)
    _validate_body_solvers(systems, body_solvers)
    _validate_influence_backend(:backend_wake, backend_wake)
    _validate_influence_backend(:backend_system, backend_system)
    _validate_solve_backend(systems, body_solvers, backend_solve)
    audit_monitors(monitors)

    # Flip body.needs_velocity_gradient based on monitor contracts so that
    # requires_hessian(body) propagates the right HS flag into the per-step
    # influence! calls. Done once before the time loop; bodies are mutable
    # so the Ref persists across steps.
    needs_grad = any(monitor_requires_body_hessian, monitors)
    needs_induced_vorticity = any(monitor_requires_induced_vorticity, monitors)
    for sys in systems_tuple
        sys.needs_velocity_gradient[] = needs_grad
    end

    # create save path if it does not exist
    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end

    # update control points and normals according to Neumann/Dirichlet BCs
    for sys in systems_tuple
        calc_normals!(sys)
        calc_controlpoints!(sys)
    end

    if !isnan(set_Das_eta_freestream) || !isnan(set_Das_eta_kinematic)
        initialize_Das!(systems_tuple, frames, Uinf, t_range[1], t_range[2] - t_range[1];
            set_Das_eta_kinematic, set_Das_eta_freestream,
            set_Das_min_kinematic_displacement)
    end

    # begin simulation
    i_step = start_step
    for t in @view t_range[start_step+1:end]
        if verbose
            println("\tstep $(i_step)/$(length(t_range)-1) at time $(t)")
            # flush(stdout)
        end

        #------- controls -------#

        # update frames based on maneuver
        # (RPMs, tilting systems, prescribed trajectory, etc.)
        dynamics_toggle = maneuver!(frames, systems_tuple, wakes_tuple, t)

        #------- aerodynamics -------#

        uinf = Uinf(t)

        # dt for this step
        dt = i_step < length(t_range) - 1 ? t_range[i_step+2] - t_range[i_step+1] : t_range[i_step+1] - t_range[i_step]

        _steady_aerodynamics!(systems, systems_tuple, wakes_tuple, frames, uinf,
            body_solvers; backend_wake, backend_solve, backend_system,
            needs_induced_vorticity, update_trailing_edges=true,
            wakerow_no_hessian_to_particles,
            body_hessian_to_particles,
            body_on_wake,
            panel_wake_on_particles,
            particle_hessian_self,
            diagnose_particle_influence,
            diagnostic_vertical)

        #------- other solvers -------#

        # e.g. structures, acoustics, dynamics, etc.

        #------- update state -------#


        #------- monitors -------#

        # Run monitors before VTK write so monitor-owned fields can be passed to
        # downstream monitors and output files.
        monitor_context = MonitorContext()
        monitor_set_time!(monitor_context, t)
        for monitor in monitors
            _run_monitor!(monitor, monitor_context, systems_tuple, wakes_tuple, frames, uinf, i_step, dt, t)
        end

        #------- save state -------#

        if !isnothing(path)
            metadata_path = _metadata_toml_path(path, name)
            if i_step == start_step || !isfile(metadata_path)
                _write_metadata_toml(path, name, systems_tuple, wakes_tuple, frames, t_range,
                    body_solvers, backend_wake, backend_solve, backend_system, monitors;
                    start_step=start_step,
                    set_Das_eta_kinematic=set_Das_eta_kinematic,
                    set_Das_eta_freestream=set_Das_eta_freestream,
                    set_Das_min_kinematic_displacement=set_Das_min_kinematic_displacement)
            end

            for (i, sys) in enumerate(systems_tuple)
                body_name = name * "_body$(i)"
                write_vtk(joinpath(path, body_name), sys, i_step, t;
                          monitors=monitors, i_system=i, overwrite=i_step==0)
            end

            for (i, w) in enumerate(wakes_tuple)
                if !isnothing(w)
                    wake_name = name * "_wake$(i)"
                    write_vtk(joinpath(path, wake_name), w, i_step, t; overwrite=i_step==0)
                end
            end

            _append_metadata_step_toml(path, name, frames, i_step, t; uinf)
        end

        #------- propagate system -------#

        if i_step < length(t_range) - 1

            #--- state evolution ---#

            # propagate wake
            for w in wakes_tuple
                if w isa PanelParticleWake
                    propagate!(w, dt; relax=particle_relax,
                        step=i_step, frames, diagnose_particle_gamma, diagnostic_vertical)
                elseif !isnothing(w)
                    propagate!(w, dt; step=i_step, frames)
                end
            end

            # propagate rigid-body kinematics
            propagate_kinematics!(systems_tuple, frames, dt)

            # update control points and normals according to Neumann/Dirichlet BCs
            for sys in systems_tuple
                calc_normals!(sys)
        calc_controlpoints!(sys)
            end

            #--- shed new wake ---#

            for (sys, w) in zip(systems_tuple, wakes_tuple)
                !isnothing(w) && shed_wake!(w, sys)
            end

        end

        # increment step
        i_step += 1
    end

    verbose && println()
end

get_Gammai(::AbstractBody{TK,NK,TF}) where {TK, NK, TF} = NK==2 ? 2 : 1
has_grad_mu(::AbstractBody{TK,NK,TF}) where {TK, NK, TF} = TK == ConstantDoublet || TK == VortexRing || TK == Union{ConstantSource, ConstantDoublet} || TK == Union{ConstantSource, VortexRing}

#------- wake shedding -------#

function update_TE!(wake::PanelWake, system::AbstractBody)

    # update first row based on system
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        shedding = system.shedding[i_surf]
        Das = system.Das[i_surf]

        # loop over shedding panels
        for i_shed in axes(shedding, 2)
            i_panel = shedding[1, i_shed]
            idx_1 = shedding[3, i_shed] # nib (second node of the shedding edge)
            v1 = FastMultipole.SVector{3}(
                system.nodes[1, system.cells[idx_1, i_panel]],
                system.nodes[2, system.cells[idx_1, i_panel]],
                system.nodes[3, system.cells[idx_1, i_panel]],
            )
            nodes[:, 1, i_shed] .= v1 .+ view(Das, :, i_shed) # nib = vertex i_shed in Das
        end

        # final node of this edge (nia of last edge = last vertex in Das)
        i_panel = shedding[1, end]
        idx_2 = shedding[2, end] # nia (first node of the shedding edge)
        v2 = FastMultipole.SVector{3}(
            system.nodes[1, system.cells[idx_2, i_panel]],
            system.nodes[2, system.cells[idx_2, i_panel]],
            system.nodes[3, system.cells[idx_2, i_panel]],
        )
        nodes[:, 1, end] .= v2 .+ view(Das, :, size(Das, 2)) # nia of last edge = last vertex
    end
end

update_TE!(wake::PanelParticleWake, system::AbstractBody) = update_TE!(wake.panel_wake, system)

# update_TE!(w::ParticleWake, sys) = update_TE!(w.panel_wake, sys)

function shed_wake!(wake::PanelWake, system::AbstractBody)

    # check storage dimensions
    @assert length(system.Das) == length(system.shedding) == length(wake.nodes) "Length of system.Das ($(length(system.Das))) must match number of surfaces in wake ($(length(wake.nodes)))"

    # shift panels back a row
    n_rows = size(wake.nodes[1], 2)
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        strength = wake.strength[i_surf]

        for j_row in min(wake.nwakes[]+1, n_rows-1):-1:1
            nodes[:, j_row+1, :] .= nodes[:, j_row, :]
        end
        for j_row in min(wake.nwakes[], n_rows-1):-1:1
            strength[:, j_row+1, :] .= strength[:, j_row, :]
        end
    end

    # update first row of strengths
    for i_surf in eachindex(wake.strength)
        for i_shed in axes(wake.strength[i_surf], 3)
            strengthi, strengthj = _get_wakestrength_mu(system, i_shed, i_surf)
            # Negate the strength because we are swapping the order of the nodes
            # to form a contiguous sequence (v1=nib, v4=nia instead of v1=nia, v4=nib)
            wake.strength[i_surf][1, 1, i_shed] = strengthi - strengthj
        end
    end

    # update nwakes
    if wake.nwakes[] == n_rows - 1 # about to overflow
        wake.overflowed[] = true
    end
    wake.nwakes[] = min(wake.nwakes[] + 1, n_rows - 1) # ensure we don't exceed storage

end

#------- PanelParticleWake shedding -------#

function shed_wake!(wake::PanelParticleWake, system::AbstractBody)
    pw = wake.panel_wake
    n_rows = size(pw.nodes[1], 2)
    buffer_full = pw.nwakes[] >= n_rows - 1

    if buffer_full
        # new particles
        _convert_to_particles!(wake)
    end

    # Shift panel rows (existing PanelWake method)
    shed_wake!(pw, system)
end
