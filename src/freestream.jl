function apply_freestream!(panel_array::AbstractPanels, freestream::AbstractVector{<:Number})
    # unpack panels
    (; panels, potential, velocity) = panel_array

    # potential
    for (i,panel) in enumerate(panels)
        x = panel.control_point
        potential[i] = -dot(x, freestream)
    end

    # velocity
    for i in eachindex(velocity)
        velocity[i] += freestream
    end

    return nothing
end