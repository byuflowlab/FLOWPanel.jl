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

function simulate!(system::AbstractBody, frames#=::AbstractVector{<:ReferenceFrame}=#, maneuver!::Function, Vinf::Function, t_range;
        # particle_trailing_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        # particle_unsteady_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        wake_args=(), kwargs...)

    # construct particle field
    # n_particles_per_step = get_max_particles(system, particle_trailing_methods, particle_unsteady_methods)
    # wake = ParticleField(n_particles_per_step * length(t_range), TF; Uinf=Vinf, wake_args...)
    wake = PanelWake(system)
    
    # begin simulation
    simulate!(system, wake, frames, maneuver!, Vinf, t_range) 
        # particle_trailing_methods, particle_unsteady_methods, kwargs...)

    return wake
end

function simulate!(system::AbstractBody{TK,NK,TF}, wake::PanelWake, frames#=::AbstractVector{<:ReferenceFrame}=#, maneuver!::Function, Uinf::Function, t_range; # Ωinf=(t)->SVector{3}(0.0, 0.0, 0.0);
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
                shrink=true,
                recenter=false,
            ),
        # trailing_vortices=fill(false, length(system.surfaces)),
        # shedding_surfaces=fill(true, length(system.surfaces)),
        rho=1.225,
        monitors=(),
    ) where {TK, NK, TF}
    # create save path if it does not exist
    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end

    # empty wake shedding locations
    # empty_wake_shedding_locations = fill(nothing, length(system.surfaces))

    # no wake panels to begin with
    wake.nwakes[] = 0

    # freestream for initial step
    uinf = Uinf(t_range[1])

    # set wake shedding locations
    dt = t_range[2] - t_range[1]
    update_wake_shedding_locations!(system, uinf, frames, dt, eta)

    # begin simulation
    i_step = 0
    for t in t_range
        println("\tstep $(i_step)/$(length(t_range)-1) at time $(t)")
        
        #------- reset system -------#

        for i_surf in eachindex(wake.velocity)
            wake.velocity[i_surf] .= zero(eltype(wake.velocity[i_surf]))
        end
        system.potential .= zero(eltype(system.potential)) # reset potential at control points to zero before solving for new circulation
        system.velocity .= zero(eltype(system.velocity)) # reset velocity at control points to zero before solving n-body problem

        # particle field
        # FLOWVPM._reset_particles(wake)
        # FLOWVPM._reset_particles_sfs(wake)

        #------- controls -------#

        # update frames based on maneuver
        # (RPMs, tilting systems, prescribed trajectory, etc.)
        dynamics_toggle = maneuver!(frames, system, wake, t)

        # update kinematic velocity due to rigid body motion
        # (structural deflections should be remembered from the previous step)
        # NOTE: this skips the top level frame, which is captured in system.fs
        # CORRECTION: just changed this to not skip the top level frame
        # current_surfaces = system.surfaces
        # kinematic_velocity!(Vcp, Vh, Vv, Vte, current_surfaces, frames; skip_top_level=false)
        
        #------- aerodynamics -------#

        uinf = Uinf(t)
        solve!(system, wake, uinf, t;
            body_solver, backend)

        #--- forces and moments ---#

        add_field(system, "U", "vector", collect(eachcol(system.velocity)), "cell")
        if has_grad_mu(system)
            calcfield_Ugradmu(system; Gammai=get_Gammai(system))
            addfields(system, "Ugradmu", "U")
        end
        calcfield_Cp(system, norm(uinf))
        calcfield_F(system, norm(uinf), rho)
        
        #------- other solvers -------#
        
        # e.g. structures, acoustics, dynamics, etc.
        
        #------- update state -------#
        

        #------- save state -------#

        if !isnothing(path)
            # panel body
            save(system, joinpath(path, name * "_vehicle_step_$i_step"))

            # particle field

            # visualize wake
            write_vtk(joinpath(path, name), wake, i_step, t)
            
        end

        for monitor in monitors
            monitor(system, wake, i_step)
        end

        #------- propagate system -------#

        if i_step < length(t_range) - 1

            #--- state evolution ---#
            
            # propagate wake
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
            
            # # propagate rigid-body kinematics
            # propagate_kinematics!(system, frames, dt)

            # next step's freestream
            idx = i_step == length(t_range) - 1 ? i_step + 1 : i_step + 2
            uinf = Uinf(t_range[idx])
            # Ω = Ωinf(t_range[idx])
            # fs = Freestream(frames[1], ref, vinf)
            # fs = velocity_to_freestream(vinf, Ω)
            # system.freestream[] = fs

            # next step's freestream
            uinf = Uinf(t_range[i_step + 2])

            # next step's dt
            dt = t_range[i_step + 2] - t_range[i_step + 1]

            # update wake shedding locations / wake leading edge
            update_wake_shedding_locations!(system, uinf, frames, dt, eta)

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
has_grad_mu(::AbstractBody{TK,NK,TF}) where {TK, NK, TF} = TK == ConstantDoublet || TK == VortexRing || TK == Union{ConstantSource, ConstantDoublet}

function update_wake_shedding_locations!(system, uinf, frames, dt, eta)
    for i in eachindex(system.Das)
        for j in axes(system.Das[i], 2)
            system.Das[i][:, j] .= uinf .* (eta * dt)
            system.Dbs[i][:, j] .= uinf .* (eta * dt)
        end
    end
end

#------- wake shedding -------#

function shed_wake!(wake::PanelWake, system::AbstractBody)

    # check storage dimensions
    @assert length(system.Das) == length(system.shedding) == length(wake.nodes) "Length of system.Das ($(length(system.Das))) must match number of surfaces in wake ($(length(wake.nodes)))"
    
    # shift panels back a row
    n_rows = size(wake.nodes[1], 1)
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        for j_row in min(wake.nwakes[]+1, n_rows-1):-1:1
            nodes[j_row+1, :] .= nodes[j_row, :]
        end
    end

    # update nwakes
    wake.nwakes[] = min(wake.nwakes[] + 1, n_rows - 1) # ensure we don't exceed storage

    # update first row based on system
    for i_surf in eachindex(wake.nodes)
        nodes = wake.nodes[i_surf]
        shedding = system.shedding[i_surf]
        Das = system.Das[i_surf]

        # loop over shedding panels
        for i_shed in axes(shedding, 2)
            i_panel = shedding[1, i_shed]
            idx_1 = shedding[2, i_shed]
            v1 = FastMultipole.SVector{3}(
                system.grid._nodes[1, system.cells[idx_1, i_panel]],
                system.grid._nodes[2, system.cells[idx_1, i_panel]],
                system.grid._nodes[3, system.cells[idx_1, i_panel]],
            )
            nodes[:, 1, i_shed] .= v1 .+ view(Das, :, i_shed) # shift by Das to get shedding location
        end

        # final node of this edge
        i_panel = shedding[1, end]
        idx_2 = shedding[3, end]
        v2 = FastMultipole.SVector{3}(
            system.grid._nodes[1, system.cells[idx_2, i_panel]],
            system.grid._nodes[2, system.cells[idx_2, i_panel]],
            system.grid._nodes[3, system.cells[idx_2, i_panel]],
        )
        nodes[:, 1, end] .= v2 .+ view(system.Das[i_surf], :, size(system.Das[i_surf], 2)) # shift trailing edge by Das
    end
end

# abstract type WakeSheddingMethod end

# struct NoShed <: WakeSheddingMethod end

# struct SigmaPPS{TF} <: WakeSheddingMethod
#     sigma::TF
#     p_per_step::Int
# end

# struct SigmaOverlap{TF} <: WakeSheddingMethod
#     sigma::TF
#     overlap::TF
# end

# struct OverlapPPS{TF} <: WakeSheddingMethod
#     overlap::TF
#     p_per_step::Int
# end

# function shed_wake!(pfield::FLOWVPM.ParticleField, system, dt, Γ, dΓdt,
#         shedding_trailing::AbstractVector{<:WakeSheddingMethod}, shedding_unsteady::AbstractVector{<:WakeSheddingMethod})
#     # shed trailing edge particles
#     shed_trailing_edge!(pfield, system.surfaces, system.wakes, Γ, shedding_trailing)

#     # shed unsteady particles
#     shed_unsteady!(pfield, system.surfaces, system.wakes, dΓdt, dt, shedding_unsteady)
# end

# function shed_trailing_edge!(pfield::FLOWVPM.ParticleField, surfaces, wakes, Γ, shedding_methods)
#     # loop over surfaces
#     iΓ = 0
#     for isurf = eachindex(surfaces)
#         surface = surfaces[isurf]
#         wake = wakes[isurf]
#         method = shedding_methods[isurf]
#         nc, ns = size(surface)
#         Γlast = zero(eltype(Γ))
#         for j in 1:ns
#             # strength
#             iΓ += nc

#             # get vertices
#             panel = wake[1, j]
#             r2 = top_left(panel)
#             r1 = bottom_left(panel)

#             # shed left particles
#             Γthis = Γ[iΓ]
#             shed_particles!(pfield, r1, r2, Γthis - Γlast, method)

#             # recurse
#             Γlast = Γthis
#         end

#         # get vertices
#         panel = wake[1, end]
#         r1 = top_right(panel)
#         r2 = bottom_right(panel)

#         # shed right particles
#         shed_particles!(pfield, r1, r2, Γlast, method)
#     end
# end

# function shed_unsteady!(pfield::FLOWVPM.ParticleField, surfaces, wakes, dΓdt, dt, shedding_methods)
#     # loop over surfaces
#     iΓ = 0
#     for isurf = eachindex(surfaces)
#         surface = surfaces[isurf]
#         wake = wakes[isurf]
#         method = shedding_methods[isurf]
#         nc, ns = size(surface)
#         for j in 1:ns
#             # strength
#             iΓ += nc

#             # get vertices
#             panel = wake[1, j]
#             r2 = bottom_left(panel)
#             r1 = bottom_right(panel)
#             Γ = dΓdt[iΓ] * dt

#             # shed unsteady particles
#             shed_particles!(pfield, r1, r2, Γ, method)
#         end
#     end
# end

# function shed_particles!(pfield, r1, r2, Γ, method::OverlapPPS)
#     # shed particles with overlap and p_per_step
#     overlap = method.overlap
#     p_per_step = method.p_per_step
#     sigma = norm(r2 - r1) * overlap / p_per_step
#     return shed_particles!(pfield, r1, r2, Γ, SigmaPPS(sigma, p_per_step))
# end

# function shed_particles!(pfield, r1, r2, Γ, method::SigmaOverlap)
#     # shed particles with sigma and overlap
#     sigma = method.sigma
#     overlap = method.overlap
#     p_per_step = ceil(Int, overlap * norm(r2 - r1) / sigma)
#     return shed_particles!(pfield, r1, r2, Γ, SigmaPPS(sigma, p_per_step))
# end

# function shed_particles!(pfield, r1, r2, Γ, method::SigmaPPS)
#     # shed particles with sigma and p_per_step
#     sigma = method.sigma
#     p_per_step = method.p_per_step

#     # add particles
#     distance_vector = (r2 - r1) / p_per_step
#     Xp = r1 + distance_vector * 0.5
#     Γp = Γ * distance_vector
#     for i in 1:p_per_step
#         FLOWVPM.add_particle(pfield, Xp, Γp, sigma; circulation=Γ)
#         Xp += distance_vector
#     end
# end

# function shed_particles!(pfield, r1, r2, Γ, method::NoShed)
#     # do not shed particles
#     return nothing
# end

# function get_max_particles(surface::AbstractMatrix{<:SurfacePanel}, method::Union{<:SigmaPPS, <:OverlapPPS})
#     pps = method.p_per_step
#     _, ns = size(surface)
#     np = (ns + 1) * pps # pps particles at each trailing edge vertex
#     return np
# end

# function get_max_particles(surface::AbstractMatrix{<:SurfacePanel}, method::SigmaOverlap)
#     # estimate p_per_step
#     pps = 8

#     return get_max_particles(surface, SigmaPPS(method.sigma, pps))
# end
    
# function get_max_particles(surface::AbstractMatrix{<:SurfacePanel}, method::NoShed)
#     # no particles shed
#     return 0
# end

# function get_max_particles(system::System, particle_trailing_methods)
#     np = 0
#     for (isurf, surface) in enumerate(system.surfaces)
#         np += get_max_particles(surface, particle_trailing_methods[isurf])
#     end
#     return np
# end

# function get_max_particles(system, particle_trailing_methods, particle_unsteady_methods)
#     np_trailing = get_max_particles(system, particle_trailing_methods)
#     np_unsteady = get_max_particles(system, particle_unsteady_methods)
#     return np_trailing + np_unsteady
# end
