struct ForcesMonitor{TF,F}
    CF::Vector{FastMultipole.SVector{3,TF}}
    CM::Vector{FastMultipole.SVector{3,TF}}
    frame::F
end

function ForcesMonitor(nt::Int, TF=Float64; frame=Body())
    CF = zeros(FastMultipole.SVector{3,TF}, nt)
    CM = zeros(FastMultipole.SVector{3,TF}, nt)

    return ForcesMonitor{TF,typeof(frame)}(CF, CM, frame)
end

function (monitor::ForcesMonitor)(system::AbstractBody, wake, i_step::Int)
    CF, CM = body_forces(system.surfaces, system.properties,
                            system.reference[], system.freestream[], 
                            system.symmetric, monitor.frame)
    monitor.CF[i_step + 1] = CF
    monitor.CM[i_step + 1] = CM
end

struct FrameForcesMonitor{TF,F}
    CF::Vector{FastMultipole.SVector{3,TF}}
    CM::Vector{FastMultipole.SVector{3,TF}}
    frame::F
end

function simulate!(system::AbstractBody{TK,NK,TF}, wake::AbstractFreeWake, frames#=::AbstractVector{<:ReferenceFrame}=#, maneuver!::Function, Uinf::Function, t_range; # Ωinf=(t)->SVector{3}(0.0, 0.0, 0.0);
        name="default_sim", path="./default_simulation",
        # vtk_args=(trailing_vortices=false, write_wakes=false), vtk_postshed=false,
        # fmm_wake_args=(), fmm_vehicle_args=(),
        # derivatives=false, nonlinear_analysis=false, nonlinear_args=(),
        # particle_trailing_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        # particle_unsteady_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        body_solver=BackslashDirichlet(system), 
        backend=FastMultipoleBackend(;
                expansion_order=10,
                multipole_acceptance=0.4,
                leaf_size=100,
            ),
        # trailing_vortices=fill(false, length(system.surfaces)),
        # shedding_surfaces=fill(true, length(system.surfaces)),
        rho=1.225,
        monitors=(),
        cp_correct_kuttacondition=true,
        cp_clip=nothing,
        verbose=false
    ) where {TK, NK, TF}
    # create save path if it does not exist
    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end

    # update control points and normals
    calc_normals!(system)
    calc_controlpoints!(system; off=abs(system.CPoffset))

    # begin simulation
    i_step = 0
    for t in t_range
        verbose && println("\tstep $(i_step)/$(length(t_range)-1) at time $(t)")

        # particle field
        # FLOWVPM._reset_particles(wake)
        # FLOWVPM._reset_particles_sfs(wake)

        #------- controls -------#

        # update frames based on maneuver
        # (RPMs, tilting systems, prescribed trajectory, etc.)
        dynamics_toggle = maneuver!(frames, system, wake, t)
        
        #------- aerodynamics -------#
        
        # store -φ_old for dφ/dt computation (before reset wipes potential)
        system.dphidt .= .-system.potential

        # reset potential/velocity
        reset!(wake)
        reset!(system)

        # get probes
        targets = (system, get_probes(wake)...)
        wake_sources = Tuple(w for w in get_sources(wake))

        # freestream
        uinf = Uinf(t)
        apply_freestream!(system, uinf)
        apply_freestream!(wake, uinf)

        # kinematics
        kinematic_velocity!((system,), frames)

        # dt for this step
        dt = i_step < length(t_range) - 1 ? t_range[i_step+2] - t_range[i_step+1] : t_range[i_step+1] - t_range[i_step]

        # snap first row of wake nodes to the trailing edge
        update_TE!(wake, system)

        # apply wake potential to body surface
        evaluate_influence!(targets, wake_sources, backend; scalar_potential=true, gradient=true, hessian=Tuple(requires_hessian(sys) for sys in targets))

        # solve system (shouldn't modify system velocity, but will update system strength)
        solve!(system, body_solver; backend)

        # update control points (normals should not have changed)
        calc_controlpoints!(system; off=abs(system.CPoffset))

        # system-on-all influence
        evaluate_influence!(targets, (system,), backend; scalar_potential=true, gradient=true, hessian=Tuple(requires_hessian(sys) for sys in targets))

        #--- forces and moments ---#

        # compute dφ/dt: dphidt holds -φ_old, add φ_new and divide by dt
        system.dphidt .= (system.dphidt .+ system.potential) ./ dt

        calcfield_Cp!(system, norm(uinf); dphidt=system.dphidt, correct_kuttacondition=cp_correct_kuttacondition, clip=cp_clip)
        calcfield_F!(system, norm(uinf), rho)
        
        #------- other solvers -------#
        
        # e.g. structures, acoustics, dynamics, etc.
        
        #------- update state -------#
        

        #------- save state -------#

        if !isnothing(path)
            # panel body
            write_vtk(joinpath(path, name), system, i_step, t; overwrite=i_step==0)

            # wake
            write_vtk(joinpath(path, name*"_wake"), wake, i_step, t; overwrite=i_step==0)
        end

        for monitor in monitors
            monitor(system, wake, i_step)
        end

        #------- propagate system -------#

        if i_step < length(t_range) - 1

            #--- state evolution ---#
            
            # propagate wake (and save )
            propagate!(wake, dt)

            # dynamics function
            # if dynamics_toggle
            #     apply_dynamics!(system, frames)
            # end

            # # calculate next step's wake trailing edge
            # this_V = nothing # ignore wake- and vehicle-induced velocity for wake shedding location update
            # update_vpm_shedding_TE!(wakes, ref, fs, dt, additional_velocity, this_V) # uses current step's freestream

            # # store trailing edge location for next step's wsl
            # store_trailing_edge!(wake_shedding_locations, current_surfaces)
            
            # propagate rigid-body kinematics
            propagate_kinematics!(system, frames, dt)

            # update control points and normals
            calc_normals!(system)
            calc_controlpoints!(system; off=abs(system.CPoffset))

            # # next step's freestream
            # idx = i_step == length(t_range) - 1 ? i_step + 1 : i_step + 2
            # uinf = Uinf(t_range[idx])
            # # Ω = Ωinf(t_range[idx])
            # # fs = Freestream(frames[1], ref, vinf)
            # # fs = velocity_to_freestream(vinf, Ω)
            # # system.freestream[] = fs

            # # next step's freestream
            # uinf = Uinf(t_range[i_step + 2])

            # # next step's dt
            # dt = t_range[i_step + 2] - t_range[i_step + 1]

            #--- shed new wake ---#

            shed_wake!(wake, system)

            # shed_wake!(wake, system,  dt, Γ_wake, dΓdt,
            #     particle_trailing_methods, particle_unsteady_methods)

            # update wake shedding locations based on wake and vehicle
            # accounts for vehicle-induced, wake-induced, freestream,
            # and kinematic velocities
            # update_vpm_shedding_LE!(wakes, ref, fs, dt, additional_velocity, V)

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

