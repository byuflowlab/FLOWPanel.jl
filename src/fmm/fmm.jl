import DelimitedFiles
#------- FastMultipole compatibility functions -------#

Base.eltype(::AbstractPanels{<:Any,TF,<:Any,<:Any}) where TF = TF

function FastMultipole.strength_dims(system::AbstractPanels{<:Any,<:Any,NK,<:Any}) where NK
    return NK
end

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::AbstractPanels{<:Any,<:Any,NK,NS}, i_body) where {NK,NS}
    # position
    panel = system.panels[i_body]

    # centroid and radius
    buffer[1:3, i_buffer] .= panel.control_point
    buffer[4, i_buffer] = panel.radius

    # strength
    strength = panel.strength
    for i in 1:NK
        buffer[4+i, i_buffer] = strength[i]
    end

    # vertices
    i_offset = 0
    # for i_vertex in NS:-1:1
    for i_vertex in 1:NS
        vertex = panel.vertices[i_vertex]
        for i in 1:3
            buffer[4+NK+i+i_offset, i_buffer] = vertex[i]
        end
        i_offset += 3
    end

    # regularization radius
    σ = 1e-2
    buffer[5+NK+3*NS,i_buffer] = σ
end

function FastMultipole.data_per_body(system::AbstractPanels{<:Any,<:Any,NK,NS}) where {NK,NS}
    return 5+NK+3*NS
end

function FastMultipole.reset!(system::AbstractPanels{<:Any,TF,<:Any,<:Any}) where TF
    system.potential .= zero(TF)
    system.velocity .= zero(eltype(system.velocity))
end

function FastMultipole.get_position(system::AbstractPanels, i)
    return system.panels[i].control_point
end

FastMultipole.get_n_bodies(system::AbstractPanels) = length(system.panels)

FastMultipole.body_to_multipole!(system::AbstractPanels{ConstantSource,<:Any,<:Any,<:Any}, args...) = FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Source}, system, args...)

FastMultipole.body_to_multipole!(system::AbstractPanels{ConstantNormalDoublet,<:Any,<:Any,<:Any}, args...) = FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Dipole}, system, args...)

FastMultipole.body_to_multipole!(system::AbstractPanels{ConstantSourceNormalDoublet,<:Any,<:Any,<:Any}, args...) = FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.SourceDipole}, system, args...)

function FastMultipole.buffer_to_target_system!(target_system::AbstractPanels, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    # get values
    TF = eltype(target_buffer)
    scalar_potential = PS ? FastMultipole.get_scalar_potential(target_buffer, i_buffer) : zero(TF)
    velocity = VS ? FastMultipole.get_velocity(target_buffer, i_buffer) : zero(SVector{3,TF})
    velocity_gradient = GS ? FastMultipole.get_velocity_gradient(target_buffer, i_buffer) : zero(SMatrix{3,3,TF,9})

    # update system
    target_system.potential[i_target] += scalar_potential
    target_system.velocity[i_target] += velocity
end

# function convolve_kernel!(target_system, target_index, panel, kernel::AbstractUnrotatedKernel, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {PS,VS,GS}
#     # iterate over targets
#     for i_target in target_index
#         target = FastMultipole.get_position(target_system, i_target)
#         potential, velocity, gradient = _induced(target, panel, kernel, derivatives_switch)
#         if PS
#             FastMultipole.set_scalar_potential!(target_system, potential, i_target)
#         end
#         if VS
#             FastMultipole.set_velocity!(target_system, velocity, i_target)
#         end
#         if GS
#             FastMultipole.set_velocity_gradient!(target_system, velocity_gradient, i_target)
#         end
#     end
# end

function convolve_kernel!(target_system, target_index, source_system::AbstractPanels{<:Any,<:Any,<:Any,NS}, source_buffer, i_panel, kernel::AbstractRotatedKernel, derivatives_switch::DerivativesSwitch{PS,VS,GS}) where {NS,PS,VS,GS}
    # get vertices
    # # unpack vertices container
    # vertices = source_system.vertices

    # # populate
    # for i in 1:NS
    #     vertices[i] = FastMultipole.get_vertex(source_buffer, source_system, i_panel, i)
    # end
    vertices = SVector{NS}(FastMultipole.get_vertex(source_buffer, source_system, i_panel, i) for i in 1:NS)

    # other properties
    control_point = FastMultipole.get_position(source_buffer, i_panel)
    strength = FastMultipole.get_strength(source_buffer, source_system, i_panel)
    sigma = source_system.sigma
    normal = FastMultipole.get_normal(source_buffer, source_system, i_panel)

    # rotate into source panel frame
    R = rotate_to_panel(normal, vertices)

    # iterate over targets
    for i_target in target_index

        # get target
        target = FastMultipole.get_position(target_system, i_target)

        # compute induced velocity
        potential, velocity, velocity_gradient = _induced(target, vertices, control_point, strength, kernel, R, derivatives_switch)

        # regularize
        potential, velocity, velocity_gradient = regularize(potential, velocity, velocity_gradient, target, vertices, kernel, sigma)

        if PS
            FastMultipole.set_scalar_potential!(target_system, i_target, potential)
        end
        if VS
            FastMultipole.set_velocity!(target_system, i_target, velocity)
        end
        if GS
            FastMultipole.set_velocity_gradient!(target_system, i_target, velocity_gradient)
        end
    end

end

function FastMultipole.direct!(target_system, target_index, derivatives_switch, source_system::AbstractPanels{K,<:Any,<:Any,<:Any}, source_buffer, source_index) where K
    kernel = K()
    for i_panel in source_index
        convolve_kernel!(target_system, target_index, source_system, source_buffer, i_panel, kernel, derivatives_switch)
    end

    return nothing
end
