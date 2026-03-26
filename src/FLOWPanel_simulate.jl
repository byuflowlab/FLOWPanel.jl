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
        eta=0.3, 
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

        # set wake shedding locations
        dt = i_step < length(t_range) - 1 ? t_range[i_step+2] - t_range[i_step+1] : t_range[i_step+1] - t_range[i_step]
        update_wake_shedding_locations!(system, dt, eta)

        # snap first row of wake nodes to the trailing edge
        update_TE!(wake, system)

        # apply wake potential to body surface
        evaluate_influence!(targets, wake_sources, backend; scalar_potential=false, gradient=true, hessian=Tuple(requires_hessian(sys) for sys in targets))

        # solve system (shouldn't modify system velocity, but will update system strength)
        solve2!(system, system.velocity, body_solver; backend)

        # update control points (normals should not have changed)
        calc_controlpoints!(system; off=abs(system.CPoffset))

        # system-on-all influence
        evaluate_influence!(targets, (system,), backend; scalar_potential=false, gradient=true, hessian=Tuple(requires_hessian(sys) for sys in targets))

        #--- forces and moments ---#

        calcfield_Cp!(system, norm(uinf); correct_kuttacondition=cp_correct_kuttacondition, clip=cp_clip)
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

            # update wake shedding locations / wake leading edge
            update_wake_shedding_locations!(system, dt, eta)

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

function update_wake_shedding_locations!(system, dt, eta)
    for i in eachindex(system.Das)

        # update Das (vertex-based, size (3, nshed+1))
        system.Das[i] .= system.velocity_te[i]

        # uinf and kinematic velocities are already included in Vte
        system.Das[i] .*= eta * dt
    end
end

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

#------- ParticleWake shedding -------#

# function shed_wake!(wake::ParticleWake, system::AbstractBody)
#     pw = wake.panel_wake

#     # At this point in simulate!, propagate! has already convected pw.nodes,
#     # and update_wake_shedding_locations! has updated system.Das for the next step.
#     #
#     # pw.nodes[:, 1, :] = OLD TE+Das convected by velocity*dt (trailing edge of shed strip)
#     # body TE + new Das = NEXT step's filament position (leading edge of shed strip)
#     # Shedding circulation is read directly from system.strength via _get_wakestrength_mu

#     for i_surf in eachindex(pw.nodes)
#         convected_nodes = pw.nodes[i_surf]  # (3, 1, nshed+1) — convected old TE+Das
#         strength = pw.strength[i_surf]
#         shedding = system.shedding[i_surf]
#         Das = system.Das[i_surf]
#         n_cols = size(strength, 3)

#         # Check if this surface wraps on itself
#         r1 = SVector{3}(view(convected_nodes, :, 1, 1))
#         rend = SVector{3}(view(convected_nodes, :, 1, n_cols + 1))
#         wraps = norm(r1 - rend) < 5 * eps()

#         # Read shedding circulation directly from solved body strengths
#         if wraps
#             si, sj = _get_wakestrength_mu(system, n_cols, i_surf)
#             Γ_last = si - sj
#         else
#             Γ_last = zero(eltype(strength))
#         end

#         for icol in 1:n_cols
#             si, sj = _get_wakestrength_mu(system, icol, i_surf)
#             Γ = si - sj

#             # New TE+Das (leading edge of shed strip = next step's filament position)
#             i_panel = shedding[1, icol]
#             idx_1 = shedding[3, icol]  # nib
#             r_le = SVector{3}(
#                 system.nodes[1, system.cells[idx_1, i_panel]] + Das[1, icol],
#                 system.nodes[2, system.cells[idx_1, i_panel]] + Das[2, icol],
#                 system.nodes[3, system.cells[idx_1, i_panel]] + Das[3, icol],
#             )

#             # Convected old TE+Das (trailing edge of shed strip)
#             r_te = SVector{3}(view(convected_nodes, :, 2, icol))

#             # Trailing (streamwise) particles: net spanwise gradient
#             _shed_particles!(wake.pfield, r_le, r_te, Γ - Γ_last, wake.method)

#             # Unsteady (spanwise) particles: time variation along trailing edge
#             r_te_next = SVector{3}(view(convected_nodes, :, 1, icol + 1))
#             Γ_prev = wake.prev_strength[i_surf][icol]
#             _shed_particles!(wake.pfield, r_te, r_te_next, Γ - Γ_prev, wake.method)

#             # Update tracking
#             wake.prev_strength[i_surf][icol] = Γ
#             Γ_last = Γ
#         end

#         if !wraps
#             # Right trailing (streamwise) particles
#             si, sj = _get_wakestrength_mu(system, n_cols, i_surf)
#             Γ = si - sj

#             # Last new TE+Das node (nia of last shedding edge + Das)
#             i_panel = shedding[1, end]
#             idx_2 = shedding[2, end]  # nia
#             das_col = size(Das, 2)
#             r_le = SVector{3}(
#                 system.nodes[1, system.cells[idx_2, i_panel]] + Das[1, das_col],
#                 system.nodes[2, system.cells[idx_2, i_panel]] + Das[2, das_col],
#                 system.nodes[3, system.cells[idx_2, i_panel]] + Das[3, das_col],
#             )
#             r_te = SVector{3}(view(convected_nodes, :, 1, n_cols + 1))
#             _shed_particles!(wake.pfield, r_le, r_te, -Γ, wake.method)
#         end
#     end

#     # No-op shift (nwakerows=0, nwakes stays 0)
#     shed_wake!(pw, system)
# end

#------- PanelParticleWake simulate! -------#

# function simulate!(system::AbstractBody{TK,NK,TF}, wake::PanelParticleWake, frames, maneuver!::Function, Uinf::Function, t_range;
#         name="default_sim", path="./default_simulation",
#         eta=0.3,
#         body_solver=BackslashDirichlet(system),
#         backend=FastMultipoleBackend(;
#                 expansion_order=10,
#                 multipole_acceptance=0.4,
#                 leaf_size=100,
#             ),
#         rho=1.225,
#         monitors=(),
#         cp_correct_kuttacondition=true,
#         cp_clip=nothing,
#         verbose=false
#     ) where {TK, NK, TF}

#     # create save path if it does not exist
#     if !isnothing(path) && !isdir(path)
#         mkpath(path)
#     end

#     pw = wake.panel_wake
#     pf = wake.pfield

#     # no wake panels to begin with
#     pw.nwakes[] = 0

#     # freestream for initial step
#     uinf = Uinf(t_range[1])
#     apply_freestream!(system, uinf)

#     # velocity due to kinematics
#     kinematic_velocity!((system,), frames)

#     # set wake shedding locations
#     dt = t_range[2] - t_range[1]
#     update_wake_shedding_locations!(system, uinf, dt, eta)

#     # initialize first row of wake nodes at shedding locations
#     update_TE!(wake, system)

#     # update control points and normals
#     calc_normals!(system)
#     calc_controlpoints!(system; off=abs(system.CPoffset))

#     # begin simulation
#     i_step = 0
#     for t in t_range
#         verbose && println("\tstep $(i_step)/$(length(t_range)-1) at time $(t), particles=$(pf.np)")

#         #------- controls -------#

#         dynamics_toggle = maneuver!(frames, system, wake, t)

#         #------- aerodynamics -------#

#         # reset potential/velocity
#         reset!(wake)
#         reset!(system)

#         # get probes (delegates to panel_wake)
#         wake_targets = get_probes(wake)
#         wake_sources = get_sources(wake)

#         # freestream
#         uinf = Uinf(t)
#         apply_freestream!(system, uinf)
#         apply_freestream!(wake, uinf)

#         # kinematics
#         kinematic_velocity!((system,), frames)

#         # snap newly shed wake to the trailing edge
#         update_TE!(wake, system)

#         # wake-on-all velocity (particles + panel wake as sources)
#         has_particles = pf.np > 0
#         if has_particles
#             targets = (system, wake_targets...)
#             evaluate_influence!(
#                 targets,
#                 (wake_sources...),
#                 backend;
#                 scalar_potential=false,
#                 gradient=true,
#                 hessian=Tuple(requires_hessian(sys) for sys in targets)
#             )
#         else
#             evaluate_influence!(
#                 (system, wake_probes...),
#                 (pw,),
#                 backend;
#                 scalar_potential=(false, false),
#                 gradient=(true, true),
#                 hessian=(requires_hessian(system), false)
#             )
#         end

#         # solve system
#         solve2!(system, system.velocity, body_solver; backend)

#         # update control points (normals should not have changed)
#         calc_controlpoints!(system; off=abs(system.CPoffset))

#         # body-on-all influence
#         if has_particles
#             evaluate_influence!(
#                 (system, wake_probes, pf),
#                 (system,),
#                 backend;
#                 scalar_potential=(false, false, false),
#                 gradient=(true, true, true),
#                 hessian=(requires_hessian(system), false, false)
#             )
#         else
#             evaluate_influence!(
#                 (system, wake_probes),
#                 (system,),
#                 backend;
#                 scalar_potential=(false, false),
#                 gradient=true,
#                 hessian=(requires_hessian(system), false)
#             )
#         end

#         #--- forces and moments ---#

#         calcfield_Cp!(system, norm(uinf); correct_kuttacondition=cp_correct_kuttacondition, clip=cp_clip)
#         calcfield_F!(system, norm(uinf), rho)

#         #------- save state -------#

#         if !isnothing(path)
#             write_vtk(joinpath(path, name), system, i_step, t; overwrite=i_step==0)
#             write_vtk(joinpath(path, name*"_wake"), wake, i_step, t; overwrite=i_step==0)
#         end

#         for monitor in monitors
#             monitor(system, wake, i_step)
#         end

#         #------- propagate system -------#

#         if i_step < length(t_range) - 1

#             # propagate wake (panels + particles)
#             propagate!(wake, dt)

#             # propagate rigid-body kinematics
#             propagate_kinematics!(system, frames, dt)

#             # update control points and normals
#             calc_normals!(system)
#             calc_controlpoints!(system; off=abs(system.CPoffset))

#             # next step's freestream
#             uinf = Uinf(t_range[i_step + 2])

#             # next step's dt
#             dt = t_range[i_step + 2] - t_range[i_step + 1]

#             # update wake shedding locations
#             update_wake_shedding_locations!(system, uinf, dt, eta)

#             # shed new wake (saves last row → shifts → converts to particles)
#             shed_wake!(wake, system)

#         end

#         # increment step
#         i_step += 1
#     end
# end
