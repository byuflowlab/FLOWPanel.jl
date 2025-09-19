#------- FastMultipole compatibility functions -------#

Base.eltype(::AbstractPanels{<:Any,TF,<:Any,<:Any}) where TF = TF

function FastMultipole.strength_dims(system::AbstractPanels{<:Any,<:Any,NK,<:Any}) where NK
    return NK
end

# function get_r̃(panel::Panel, core_size, ε, ::Type{UniformSourceNormalDoublet})
    
#     # panel properties
#     ρ = panel.radius
#     A = get_area(panel.vertices)
#     μ = abs(panel.strength[2])
#     f1 = 1 / (core_size * core_size * core_size)
#     f2 = 2 * π * ε / (μ * A)

#     # solve for d
#     d_upper = cbrt(μ * A / (2 * π * ε))
#     d = Roots.find_zero(d -> exp(-d*d*d * f1) - d * d * d * f2, (zero(d_upper), d_upper), Roots.Brent())
    
#     # compute r̃
#     r̃ = ρ + d
    
#     return r̃
# end

# function get_r̃(panel::Panel, core_size, ε, ::Type{ConstantNormalDoublet})
    
#     # panel properties
#     ρ = panel.radius
#     A = get_area(panel.vertices)
#     μ = abs(panel.strength[1])
#     f1 = 1 / (core_size * core_size * core_size)
#     f2 = 2 * π * ε / (μ * A)

#     # solve for d
#     d_upper = cbrt(μ * A / (2 * π * ε))
#     d = Roots.find_zero(d -> exp(-d*d*d * f1) - d * d * d * f2, (zero(d_upper), d_upper), Roots.Brent())
#     # @show d, ρ
#     @show core_size, ε
    
#     # compute r̃
#     r̃ = ρ + d
    
#     return r̃
# end

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::AbstractPanels{TK,<:Any,NK,NS}, i_body) where {TK,NK,NS}
    # position
    panel = system.panels[i_body]

    # centroid
    buffer[1:3, i_buffer] .= panel.control_point

    # regularization radius
    # r̃ = get_r̃(panel, system.sigma[i_body], system.ε[i_body], TK)
    # buffer[4, i_buffer] = r̃
    buffer[4, i_buffer] = panel.radius + system.sigma
    # @show panel.radius / r̃

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

    # normal
    for i in 1:3
        buffer[4+NK+3*NS+i, i_buffer] = panel.normal[i]
    end

    # # regularization radius
    # σ = system.sigma[i_body]
    # buffer[5+NK+3*NS,i_buffer] = σ
end

function FastMultipole.data_per_body(system::AbstractPanels{<:Any,<:Any,NK,NS}) where {NK,NS}
    return 5+NK+3*NS+3
end

function FastMultipole.reset!(system::AbstractPanels{<:Any,TF,<:Any,<:Any}) where TF
    system.potential .= zero(TF)
    system.velocity .= zero(eltype(system.velocity))
end

function FastMultipole.get_position(system::AbstractPanels, i)
    return system.panels[i].control_point
end

function FastMultipole.has_vector_potential(::AbstractPanels) 
    return false
end

FastMultipole.get_n_bodies(system::AbstractPanels) = length(system.panels)

FastMultipole.body_to_multipole!(system::AbstractPanels{UniformSource,<:Any,<:Any,<:Any}, args...) = FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Source}, system, args...)

FastMultipole.body_to_multipole!(system::AbstractPanels{ConstantNormalDoublet,<:Any,<:Any,<:Any}, args...) = FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Dipole}, system, args...)

FastMultipole.body_to_multipole!(system::AbstractPanels{UniformSourceNormalDoublet,<:Any,<:Any,<:Any}, args...) = FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.SourceDipole}, system, args...)

function FastMultipole.buffer_to_target_system!(target_system::AbstractPanels, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    # get values
    TF = eltype(target_buffer)
    scalar_potential = PS ? FastMultipole.get_scalar_potential(target_buffer, i_buffer) : zero(TF)
    velocity = VS ? FastMultipole.get_gradient(target_buffer, i_buffer) : zero(SVector{3,TF})
    velocity_gradient = GS ? FastMultipole.get_hessian(target_buffer, i_buffer) : zero(SMatrix{3,3,TF,9})

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
        potential, velocity, velocity_gradient = _induced(target, vertices, control_point, strength, kernel, sigma, R, derivatives_switch)

        # regularize
        # potential, velocity, velocity_gradient = regularize(potential, velocity, velocity_gradient, target, vertices, kernel, sigma)

        if PS
            FastMultipole.set_scalar_potential!(target_system, i_target, potential)
        end
        if VS
            FastMultipole.set_gradient!(target_system, i_target, velocity)
        end
        if GS
            FastMultipole.set_hessian!(target_system, i_target, velocity_gradient)
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

#--- solver specific functions ---#

function FastMultipole.extra_data_per_target_body(system::AbstractPanels)
    return 3
end

function FastMultipole.extra_target_to_buffer!(buffer, system::AbstractPanels, i_body, i_sorted)
    # get normal
    normal = system.panels[i_sorted].normal

    # populate buffer
    for i in 1:3
        buffer[18+i, i_body] = normal[i]
    end
end

function FastMultipole.influence!(influence, target_buffer, panels::AbstractPanels, source_buffer)
    
    # loop over targets
    for i_target in eachindex(influence)

        # unpack data
        normal = SVector{3}(target_buffer[19, i_target], target_buffer[20, i_target], target_buffer[21, i_target])
        velocity = FastMultipole.get_gradient(target_buffer, i_target)

        # flow tangency boundary condition
        influence[i_target] = dot(velocity, normal)
    end
end

function FastMultipole.value_to_strength!(source_buffer, ::AbstractPanels{UniformSource,<:Any,<:Any,<:Any}, i_body, value)
    source_buffer[5, i_body] = value
end

function FastMultipole.value_to_strength!(source_buffer, ::AbstractPanels{ConstantNormalDoublet,<:Any,<:Any,<:Any}, i_body, value)
    source_buffer[5, i_body] = value
end

function FastMultipole.value_to_strength!(source_buffer, ::AbstractPanels{UniformSourceNormalDoublet,<:Any,<:Any,<:Any}, i_body, value)
    source_buffer[6, i_body] = value
end

function FastMultipole.strength_to_value(strength, ::AbstractPanels{UniformSource,<:Any,<:Any,<:Any})
    return strength[1]
end

function FastMultipole.strength_to_value(strength, ::AbstractPanels{ConstantNormalDoublet,<:Any,<:Any,<:Any})
    return strength[1]
end

function FastMultipole.strength_to_value(strength, ::AbstractPanels{UniformSourceNormalDoublet,<:Any,<:Any,<:Any})
    return strength[2]
end

function FastMultipole.buffer_to_system_strength!(system::AbstractPanels{<:Any,<:Any,1,<:Any}, i_body, source_buffer, i_buffer)

    # unpack panel
    panel = system.panels[i_body]
    vertices = panel.vertices
    control_point = panel.control_point
    normal = panel.normal
    # strength = panel.strength
    radius = panel.radius

    # new strength
    new_strength = SVector{1}(source_buffer[5, i_buffer])

    # assemble updated panel
    new_panel = eltype(system.panels)(vertices, control_point, normal, new_strength, radius)
    system.panels[i_body] = new_panel
end

function FastMultipole.buffer_to_system_strength!(system::AbstractPanels{<:Any,<:Any,2,<:Any}, i_body, source_buffer, i_buffer)

    # unpack panel
    panel = system.panels[i_body]
    vertices = panel.vertices
    control_point = panel.control_point
    normal = panel.normal
    source_strength = panel.strength[1]
    radius = panel.radius

    # new strength
    new_strength = SVector{2}(source_strength, source_buffer[6, i_buffer])

    # assemble updated panel
    new_panel = eltype(system.panels)(vertices, control_point, normal, new_strength, radius)
    system.panels[i_body] = new_panel
end

function FastMultipole.target_influence_to_buffer!(target_buffer, i_buffer, derivatives_switch, target_system::AbstractPanels, i_target)
    # unpack velocity
    velocity = target_system.velocity[i_target]

    # populate buffer
    for i in 1:3
        target_buffer[4+i, i_buffer] = velocity[i]
    end
end