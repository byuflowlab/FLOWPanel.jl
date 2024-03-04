Base.getindex(sys::AbstractPanels, i) = sys.panels[i], sys.potential[i], sys.velocity[i]

Base.getindex(sys::AbstractPanels, i, ::FLOWFMM.Position) = sys.panels[i].control_point

Base.getindex(sys::AbstractPanels, i, ::FLOWFMM.Radius) = sys.panels[i].radius

Base.getindex(sys::AbstractPanels, i, ::FLOWFMM.ScalarPotential) = sys.potential[i]

Base.getindex(sys::AbstractPanels, i, ::FLOWFMM.Velocity) = sys.velocity[i]

Base.getindex(sys::AbstractPanels, i, ::FLOWFMM.VelocityGradient) = zero(SMatrix{3,3,Float64,3})

Base.getindex(sys::AbstractPanels{ConstantNormalDoublet(),<:Any}, i, ::FLOWFMM.ScalarStrength) = sys.panels[i].strength[1]
Base.getindex(sys::AbstractPanels{ConstantSource(),<:Any}, i, ::FLOWFMM.ScalarStrength) = sys.panels[i].strength[1]
function Base.getindex(sys::AbstractPanels, i, ::FLOWFMM.ScalarStrength)
    @show typeof(sys) typeof(sys) <: AbstractPanels{ConstantSource(),<:Any}
    sys.panels[i].strength
end

# Base.getindex(sys::AbstractPanels{ConstantNormalDoublet(),<:Any}, i, ::FLOWFMM.ScalarStrength) = sys.panels[i].strength[1]

# Base.getindex(sys::AbstractPanels{ConstantSourceNormalDoublet(),<:Any}, i, ::FLOWFMM.ScalarStrength) = sys.panels[i].strength[1]

# Base.getindex(sys::AbstractPanels, i, ::FLOWFMM.VectorStrength) = zero(SVector{3,Float64})

Base.getindex(sys::AbstractPanels, i, ::FLOWFMM.Vertex, i_vertex) = sys.panels[i].vertices[i_vertex]

Base.getindex(sys::AbstractPanels, i, ::FLOWFMM.Normal) = sys.panels[i].normal

function Base.setindex!(sys::AbstractPanels, val, i)
    panel, potential, velocity = val
    sys.panels[i] = panel
    sys.potential[i] = potential
    sys.velocity[i] = velocity
end

Base.setindex!(sys::AbstractPanels, val, i, ::FLOWFMM.ScalarPotential) = sys.potential[i] = val

Base.setindex!(::AbstractPanels, val, i, ::FLOWFMM.VectorPotential) = nothing

Base.setindex!(sys::AbstractPanels, val, i, ::FLOWFMM.Velocity) = sys.velocity[i] = val

Base.setindex!(::AbstractPanels, val, i, ::FLOWFMM.VelocityGradient) = nothing

FLOWFMM.get_n_bodies(sys::AbstractPanels) = length(sys.panels)

FLOWFMM.buffer_element(sys::AbstractPanels) = deepcopy(sys.panels[1]), zero(eltype(sys.potential)), zero(eltype(sys.velocity))

Base.eltype(::AbstractPanels{<:Any,TF,<:Any,<:Any}) where TF = TF

FLOWFMM.B2M!(system::AbstractPanels{ConstantSource(),<:Any,1,4}, branch, bodies_index, harmonics, expansion_order) = FLOWFMM.B2M!_quadpanel(system, branch, bodies_index, harmonics, expansion_order, FLOWFMM.UniformSourcePanel())

FLOWFMM.B2M!(system::AbstractPanels{ConstantSource(),<:Any,1,3}, branch, bodies_index, harmonics, expansion_order) = FLOWFMM.B2M!_tripanel(system, branch, bodies_index, harmonics, expansion_order, FLOWFMM.UniformSourcePanel())

FLOWFMM.B2M!(system::AbstractPanels{ConstantNormalDoublet(),<:Any,1,4}, branch, bodies_index, harmonics, expansion_order) = FLOWFMM.B2M!_quadpanel(system, branch, bodies_index, harmonics, expansion_order, FLOWFMM.UniformNormalDipolePanel())

FLOWFMM.B2M!(system::AbstractPanels{ConstantNormalDoublet(),<:Any,1,3}, branch, bodies_index, harmonics, expansion_order) = FLOWFMM.B2M!_tripanel(system, branch, bodies_index, harmonics, expansion_order, FLOWFMM.UniformNormalDipolePanel())

function FLOWFMM.direct!(target_system::AbstractPanels, target_index, source_system::AbstractPanels{kernel,<:Any,<:Any,<:Any}, source_index) where kernel
    for i_source in source_index
        panel = source_system.panels[i_source]

        # rotate into source panel frame
        Rprime, Rxprime, Ryprime, Rzprime = rotate_to_panel(panel)

        # iterate over targets
        for i_target in target_index
            potential, velocity, _ = _induced(target_system[i_target, FLOWFMM.POSITION], panel, kernel, Rprime, Rxprime, Ryprime, Rzprime; toggle_potential=true, toggle_velocity=true, toggle_hessian=false)
            target_system[i_target, FLOWFMM.SCALAR_POTENTIAL] = target_system[i_target, FLOWFMM.SCALAR_POTENTIAL] + potential
            target_system[i_target, FLOWFMM.VELOCITY] = target_system[i_target, FLOWFMM.VELOCITY] + velocity
        end
    end
    return nothing
end

function FLOWFMM.direct!(target_system, target_index, source_system::AbstractPanels{kernel,<:Any,<:Any,<:Any}, source_index) where kernel
    for panel in source_system.panels
        # rotate into source panel frame
        Rprime, Rxprime, Ryprime, Rzprime = rotate_to_panel(panel)

        # iterate over targets
        for i_target in target_index
            potential, velocity, gradient = _induced(target_system[i_target, FLOWFMM.POSITION], panel, kernel, Rprime, Rxprime, Ryprime, Rzprime; toggle_potential=true, toggle_velocity=true, toggle_hessian=true)
            target_system[i_target, FLOWFMM.SCALAR_POTENTIAL] = target_system[i_target, FLOWFMM.SCALAR_POTENTIAL] + potential
            target_system[i_target, FLOWFMM.VELOCITY] = target_system[i_target, FLOWFMM.VELOCITY] + velocity
            target_system[i_target, FLOWFMM.VELOCITY_GRADIENT] = target_system[i_target, FLOWFMM.VELOCITY_GRADIENT] + gradient
        end
    end
    return nothing
end
