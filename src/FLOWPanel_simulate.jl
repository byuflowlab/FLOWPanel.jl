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
        _set_formulation_geometry!(sys, true)
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

        # reset potential/velocity
        for w in wakes_tuple
            !isnothing(w) && reset!(w)
        end
        for sys in systems_tuple
            reset!(sys)
        end

        # get probes
        wake_probes = _collect_wake_probes(wakes_tuple)
        targets = (systems_tuple..., wake_probes...)
        wake_sources = _collect_wake_sources(wakes_tuple)

        # freestream
        uinf = Uinf(t)
        apply_freestream!(systems_tuple, uinf)
        for w in wakes_tuple
            !isnothing(w) && apply_freestream!(w, uinf)
        end

        # kinematics
        kinematic_velocity!(systems_tuple, frames)

        # dt for this step
        dt = i_step < length(t_range) - 1 ? t_range[i_step+2] - t_range[i_step+1] : t_range[i_step+1] - t_range[i_step]

        # snap first row of wake nodes to the trailing edge
        for (sys, w) in zip(systems_tuple, wakes_tuple)
            !isnothing(w) && update_TE!(w, sys)
        end

        # apply wake velocity to body surface
        if length(wake_sources) > 0
            influence!(targets, wake_sources, backend_wake; precalc=true,
                scalar_potential=false,
                velocity=true,
                velocity_gradient=Tuple(requires_hessian(sys) for sys in targets),
                extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity))
        end

        # solve systems
        solve!(systems, body_solvers; backend=backend_solve, update_cps_normals=false)

        # update control points (normals should not have changed)
        for sys in systems_tuple
            calc_controlpoints!(sys; off=abs(sys.CPoffset))
        end

        # system-on-all influence
        influence!(targets, systems_tuple, backend_system; precalc=false,
            scalar_potential=false,
            velocity=true,
            velocity_gradient=Tuple(requires_hessian(sys) for sys in targets),
            extra_outputs=_induced_vorticity_extra_outputs(targets, needs_induced_vorticity))

        #------- other solvers -------#

        # e.g. structures, acoustics, dynamics, etc.

        #------- update state -------#


        #------- monitors -------#

        # Run monitors before VTK write so monitor-populated fields
        # (e.g. PressureBernoulli → body.P, ForceMonitor → body.F) land in
        # the output files.
        for monitor in monitors
            monitor(systems_tuple, wakes_tuple, frames, uinf, i_step, dt)
        end

        #------- save state -------#

        if !isnothing(path)
            for (i, sys) in enumerate(systems_tuple)
                body_name = name * "_body$(i)"
                write_vtk(joinpath(path, body_name), sys, i_step, t; overwrite=i_step==0)
            end

            for (i, w) in enumerate(wakes_tuple)
                if !isnothing(w)
                    wake_name = name * "_wake$(i)"
                    write_vtk(joinpath(path, wake_name), w, i_step, t; overwrite=i_step==0)
                end
            end

            # frame-state companion file for simulate_warmstart!
            _write_frame_state_toml(path, name, frames, i_step, t; truncate=(i_step==0))
        end

        #------- propagate system -------#

        if i_step < length(t_range) - 1

            #--- state evolution ---#

            # propagate wake
            for w in wakes_tuple
                !isnothing(w) && propagate!(w, dt; step=i_step, frames)
            end

            # propagate rigid-body kinematics
            propagate_kinematics!(systems_tuple, frames, dt)

            # update control points and normals according to Neumann/Dirichlet BCs
            for sys in systems_tuple
                _set_formulation_geometry!(sys, true)
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
