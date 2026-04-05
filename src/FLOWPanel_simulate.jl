"""
    ForcesMonitor(nt, TF=Float64; frame=Body())

Simple storage monitor for integrated force and moment histories over `nt`
simulation steps.
"""
struct ForcesMonitor{TF,F}
    CF::Vector{FastMultipole.SVector{3,TF}}
    CM::Vector{FastMultipole.SVector{3,TF}}
    frame::F
end

"""
    ForcesMonitor(nt, TF=Float64; frame=Body())

Construct a force/moment history monitor with `nt` output slots.
"""
function ForcesMonitor(nt::Int, TF=Float64; frame=Body())
    CF = zeros(FastMultipole.SVector{3,TF}, nt)
    CM = zeros(FastMultipole.SVector{3,TF}, nt)

    return ForcesMonitor{TF,typeof(frame)}(CF, CM, frame)
end

function (monitor::ForcesMonitor)(systems::Tuple, wakes::Tuple, i_step::Int)
    CF, CM = body_forces(systems[1].surfaces, systems[1].properties,
                            systems[1].reference[], systems[1].freestream[],
                            systems[1].symmetric, monitor.frame)
    monitor.CF[i_step + 1] = CF
    monitor.CM[i_step + 1] = CM
end

"""
    FrameForcesMonitor

Storage type for frame-resolved force and moment histories.
"""
struct FrameForcesMonitor{TF,F}
    CF::Vector{FastMultipole.SVector{3,TF}}
    CM::Vector{FastMultipole.SVector{3,TF}}
    frame::F
end

#--- dphidt dispatch helpers ---#

_store_neg_potential!(sys::AbstractLiftingBody) = (sys.dphidt .= .-sys.potential)
_store_neg_potential!(::AbstractBody) = nothing

_compute_dphidt!(sys::AbstractLiftingBody, dt) = (sys.dphidt .= (sys.dphidt .+ sys.potential) ./ dt)
_compute_dphidt!(::AbstractBody, dt) = nothing

_get_dphidt(sys::AbstractLiftingBody) = sys.dphidt
_get_dphidt(::AbstractBody) = nothing

#--- wake tuple helpers ---#

function _collect_wake_probes(wakes::Tuple)
    result = ()
    for w in wakes
        if !isnothing(w)
            result = (result..., get_probes(w)...)
        end
    end
    return result
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

#--- single-body backward-compat wrapper ---#

"""
    simulate!(system, wake, frames, maneuver!, Uinf, t_range; body_solver=BackslashDirichlet(system), optargs...)
    simulate!(systems, wakes, frames, maneuver!, Uinf, t_range; body_solvers, backend=FastMultipoleBackend(...), rho=1.225, monitors=(), ...)

Advance one or more coupled body-wake systems through `t_range`, solving the
aerodynamics, updating wakes, optionally writing VTK output, and calling any
registered monitors.
"""
function simulate!(system::AbstractBody{TK,NK,TF}, wake::AbstractFreeWake, frames, maneuver!::Function, Uinf::Function, t_range;
        body_solver=BackslashDirichlet(system), optargs...
    ) where {TK, NK, TF}
    # wrap maneuver to match tuple signature
    _maneuver!(frames, systems, wakes, t) = maneuver!(frames, systems[1], wakes[1], t)
    simulate!((system,), (wake,), frames, _maneuver!, Uinf, t_range;
        body_solvers=(body_solver,), optargs...)
end

#--- primary tuple-based simulate! ---#

function simulate!(systems::Tuple, wakes::Tuple, frames, maneuver!::Function, Uinf::Function, t_range;
        name="default_sim", path="./default_simulation",
        body_solvers::Tuple,
        backend=FastMultipoleBackend(;
                expansion_order=10,
                multipole_acceptance=0.4,
                leaf_size=100,
            ),
        rho=1.225,
        monitors=(),
        cp_correct_kuttacondition=true,
        cp_clip=nothing,
        verbose=false
    )
    # create save path if it does not exist
    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end

    # update control points and normals
    for sys in systems
        calc_normals!(sys)
        calc_controlpoints!(sys; off=abs(sys.CPoffset))
    end

    # begin simulation
    i_step = 0
    for t in t_range
        verbose && println("\tstep $(i_step)/$(length(t_range)-1) at time $(t)")

        #------- controls -------#

        # update frames based on maneuver
        # (RPMs, tilting systems, prescribed trajectory, etc.)
        dynamics_toggle = maneuver!(frames, systems, wakes, t)

        #------- aerodynamics -------#

        # store -φ_old for dφ/dt computation (before reset wipes potential)
        for sys in systems
            _store_neg_potential!(sys)
        end

        # reset potential/velocity
        for w in wakes
            !isnothing(w) && reset!(w)
        end
        for sys in systems
            reset!(sys)
        end

        # get probes
        wake_probes = _collect_wake_probes(wakes)
        targets = (systems..., wake_probes...)
        wake_sources = _collect_wake_sources(wakes)

        # freestream
        uinf = Uinf(t)
        apply_freestream!(systems, uinf)
        for w in wakes
            !isnothing(w) && apply_freestream!(w, uinf)
        end

        # kinematics
        kinematic_velocity!(systems, frames)

        # dt for this step
        dt = i_step < length(t_range) - 1 ? t_range[i_step+2] - t_range[i_step+1] : t_range[i_step+1] - t_range[i_step]

        # snap first row of wake nodes to the trailing edge
        for (sys, w) in zip(systems, wakes)
            !isnothing(w) && update_TE!(w, sys)
        end

        # apply wake velocity to body surface
        if length(wake_sources) > 0
            influence!(targets, wake_sources, backend; scalar_potential=false, gradient=true, hessian=Tuple(requires_hessian(sys) for sys in targets))
        end

        # solve systems with cross-body coupling
        solve!(systems, body_solvers; backend=fill(backend, length(systems)))

        # update control points (normals should not have changed)
        for sys in systems
            calc_controlpoints!(sys; off=abs(sys.CPoffset))
        end

        # system-on-all influence
        influence!(targets, systems, backend; scalar_potential=false, gradient=true, hessian=Tuple(requires_hessian(sys) for sys in targets))

        #--- forces and moments ---#

        # compute dφ/dt: dphidt holds -φ_old, add φ_new and divide by dt
        for sys in systems
            _compute_dphidt!(sys, dt)
        end

        calcfield_Cp!(systems, norm(uinf); dphidt=Tuple(_get_dphidt(s) for s in systems), correct_kuttacondition=fill(cp_correct_kuttacondition, length(systems)), clip=cp_clip)
        calcfield_F!(systems, norm(uinf), rho)

        #------- other solvers -------#

        # e.g. structures, acoustics, dynamics, etc.

        #------- update state -------#


        #------- save state -------#

        if !isnothing(path)
            for (i, sys) in enumerate(systems)
                body_name = length(systems) == 1 ? name : name * "_body$(i)"
                write_vtk(joinpath(path, body_name), sys, i_step, t; overwrite=i_step==0)
            end

            for (i, w) in enumerate(wakes)
                if !isnothing(w)
                    wake_name = length(systems) == 1 ? name * "_wake" : name * "_wake$(i)"
                    write_vtk(joinpath(path, wake_name), w, i_step, t; overwrite=i_step==0)
                end
            end
        end

        for monitor in monitors
            monitor(systems, wakes, i_step)
        end

        #------- propagate system -------#

        if i_step < length(t_range) - 1

            #--- state evolution ---#

            # propagate wake
            for w in wakes
                !isnothing(w) && propagate!(w, dt)
            end

            # propagate rigid-body kinematics
            propagate_kinematics!(systems, frames, dt)

            # update control points and normals
            for sys in systems
                calc_normals!(sys)
                calc_controlpoints!(sys; off=abs(sys.CPoffset))
            end

            #--- shed new wake ---#

            for (sys, w) in zip(systems, wakes)
                !isnothing(w) && shed_wake!(w, sys)
            end

        end

        # increment step
        i_step += 1
    end
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
