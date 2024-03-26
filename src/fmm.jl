Base.getindex(sys::AbstractPanels, i, ::FastMultipole.Body) = sys.panels[i], sys.potential[i], sys.velocity[i]

Base.getindex(sys::AbstractPanels, i, ::FastMultipole.Position) = sys.panels[i].control_point

Base.getindex(sys::AbstractPanels, i, ::FastMultipole.Radius) = sys.panels[i].radius

Base.getindex(sys::AbstractPanels, i, ::FastMultipole.ScalarPotential) = sys.potential[i]

Base.getindex(sys::AbstractPanels, i, ::FastMultipole.Velocity) = sys.velocity[i]

Base.getindex(sys::AbstractPanels, i, ::FastMultipole.VelocityGradient) = zero(SMatrix{3,3,Float64,3})

Base.getindex(sys::AbstractPanels{ConstantNormalDoublet(),<:Any}, i, ::FastMultipole.ScalarStrength) = sys.panels[i].strength[1]
Base.getindex(sys::AbstractPanels{ConstantSource(),<:Any}, i, ::FastMultipole.ScalarStrength) = sys.panels[i].strength[1]
function Base.getindex(sys::AbstractPanels, i, ::FastMultipole.ScalarStrength)
    @show typeof(sys) typeof(sys) <: AbstractPanels{ConstantSource(),<:Any}
    sys.panels[i].strength
end

# Base.getindex(sys::AbstractPanels{ConstantNormalDoublet(),<:Any}, i, ::FastMultipole.ScalarStrength) = sys.panels[i].strength[1]

# Base.getindex(sys::AbstractPanels{ConstantSourceNormalDoublet(),<:Any}, i, ::FastMultipole.ScalarStrength) = sys.panels[i].strength[1]

# Base.getindex(sys::AbstractPanels, i, ::FastMultipole.VectorStrength) = zero(SVector{3,Float64})

Base.getindex(sys::AbstractPanels, i, ::FastMultipole.Vertex, i_vertex) = sys.panels[i].vertices[i_vertex]

Base.getindex(sys::AbstractPanels, i, ::FastMultipole.Normal) = sys.panels[i].normal

function Base.setindex!(sys::AbstractPanels, val, i, ::FastMultipole.Body)
    panel, potential, velocity = val
    sys.panels[i] = panel
    sys.potential[i] = potential
    sys.velocity[i] = velocity
end

Base.setindex!(sys::AbstractPanels, val, i, ::FastMultipole.ScalarPotential) = sys.potential[i] = val

Base.setindex!(::AbstractPanels, val, i, ::FastMultipole.VectorPotential) = nothing

Base.setindex!(sys::AbstractPanels, val, i, ::FastMultipole.Velocity) = sys.velocity[i] = val

Base.setindex!(::AbstractPanels, val, i, ::FastMultipole.VelocityGradient) = nothing

FastMultipole.get_n_bodies(sys::AbstractPanels) = length(sys.panels)

FastMultipole.buffer_element(sys::AbstractPanels) = deepcopy(sys.panels[1]), zero(eltype(sys.potential)), zero(eltype(sys.velocity))

Base.eltype(::AbstractPanels{<:Any,TF,<:Any,<:Any}) where TF = TF

FastMultipole.B2M!(system::AbstractPanels{ConstantSource(),<:Any,1,4}, branch, bodies_index, harmonics, expansion_order) = FastMultipole.B2M!_quadpanel(system, branch, bodies_index, harmonics, expansion_order, FastMultipole.UniformSourcePanel())

FastMultipole.B2M!(system::AbstractPanels{ConstantSource(),<:Any,1,3}, branch, bodies_index, harmonics, expansion_order) = FastMultipole.B2M!_tripanel(system, branch, bodies_index, harmonics, expansion_order, FastMultipole.UniformSourcePanel())

FastMultipole.B2M!(system::AbstractPanels{ConstantNormalDoublet(),<:Any,1,4}, branch, bodies_index, harmonics, expansion_order) = FastMultipole.B2M!_quadpanel(system, branch, bodies_index, harmonics, expansion_order, FastMultipole.UniformNormalDipolePanel())

FastMultipole.B2M!(system::AbstractPanels{ConstantNormalDoublet(),<:Any,1,3}, branch, bodies_index, harmonics, expansion_order) = FastMultipole.B2M!_tripanel(system, branch, bodies_index, harmonics, expansion_order, FastMultipole.UniformNormalDipolePanel())

@inline function convolve_kernel!(target_system, target_index, panel, kernel::AbstractRotatedKernel, derivatives_switch::DerivativesSwitch{PS,<:Any,VS,GS}) where {PS,VS,GS}
    # rotate into source panel frame
    Rprime, Rxprime, Ryprime, Rzprime = rotate_to_panel(panel)

    # iterate over targets
    for i_target in target_index
        potential, velocity, gradient = _induced(target_system[i_target, FastMultipole.POSITION], panel, kernel, Rprime, Rxprime, Ryprime, Rzprime, derivatives_switch)
        if PS
            target_system[i_target, FastMultipole.SCALAR_POTENTIAL] += potential
        end
        # target_system[i_target, FastMultipole.SCALAR_POTENTIAL] = target_system[i_target, FastMultipole.SCALAR_POTENTIAL] + potential
        if VS
            target_system[i_target, FastMultipole.VELOCITY] += velocity
        end
        # target_system[i_target, FastMultipole.VELOCITY] = target_system[i_target, FastMultipole.VELOCITY] + velocity
        if GS
            target_system[i_target, FastMultipole.VELOCITY_GRADIENT] += gradient
        end
        # target_system[i_target, FastMultipole.VELOCITY_GRADIENT] = target_system[i_target, FastMultipole.VELOCITY_GRADIENT] + gradient
    end
end

@inline function convolve_kernel!(target_system, target_index, panel, kernel::AbstractUnrotatedKernel, derivatives_switch::DerivativesSwitch{PS,<:Any,VS,GS}) where {PS,VS,GS}
    # iterate over targets
    for i_target in target_index
        potential, velocity, gradient = _induced(target_system[i_target, FastMultipole.POSITION], panel, kernel, derivatives_switch::DerivativesSwitch{PS,<:Any,VS,GS}) where {PS,VS,GS}
        if PS
            target_system[i_target, FastMultipole.SCALAR_POTENTIAL] = target_system[i_target, FastMultipole.SCALAR_POTENTIAL] + potential
        end
        if VS
            target_system[i_target, FastMultipole.VELOCITY] = target_system[i_target, FastMultipole.VELOCITY] + velocity
        end
        if GS
            target_system[i_target, FastMultipole.VELOCITY_GRADIENT] = target_system[i_target, FastMultipole.VELOCITY_GRADIENT] + gradient
        end
    end
end

@inline function FastMultipole.direct!(target_system, target_index, derivatives_switch, source_system::AbstractPanels{kernel,<:Any,<:Any,<:Any}, source_index) where kernel
    for panel in source_system.panels
        convolve_kernel!(target_system, target_index, panel, kernel, derivatives_switch)
    end
    return nothing
end