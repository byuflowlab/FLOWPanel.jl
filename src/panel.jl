struct Panel{TF,NK,NS}
    vertices::SVector{NS,SVector{3,TF}}
    control_point::SVector{3,TF}
    normal::SVector{3,TF}
    strength::SVector{NK,TF}
    radius::TF
end

function get_radius(centroid, vertices...)
    radius = 0.0
    for v in vertices
        radius = max(radius, norm(v-centroid))
    end
    return radius
end

function get_area(vertices::SVector{NS,<:Any}) where NS
    area = 0.0
    for i in 1:NS-2
        area += norm(cross(vertices[i] - vertices[NS], vertices[i+1] - vertices[NS]))
    end
    return area / 2
end

#####
##### unstructured grid
#####

struct UnstructuredGrid{TK,TF,NK,NS} <: AbstractPanels{TK,TF,NK,NS}
    # WriteVTK unstructured grid
    points::Vector{SVector{3,TF}}
    meshcells::Vector{MeshCell{VTKCellType,SVector{NS,Int}}}

    # additional panel data
    control_points::Vector{SVector{3,TF}}
    normals::Vector{SVector{3,TF}}
    strengths::Vector{SVector{NK,TF}}
    potential::Vector{TF}
    velocity::Vector{SVector{3,TF}}

    # efficient panel access assembly
    panels::Vector{Panel{TF,NK,NS}}

    # wake shed locations
    wake_points::Vector{SVector{3,TF}}
end

function PanelArray(points::Vector{<:AbstractVector}, meshcells::Vector{<:MeshCell}, kernel::AbstractKernel)

    n_cells = length(meshcells)
    TF = eltype(eltype(points))
    NK = kernel_multiplicity(kernel)
    NS = meshcells[1].ctype.nodes

    for meshcell in meshcells
        @assert meshcell.ctype.nodes == NS "Meshcells inconsistent: make sure all meshcells have the same number of sides."
    end

    panels = Vector{Panel{TF,NK,NS}}(undef,n_cells)
    for (i,meshcell) in enumerate(meshcells)
        vertices = points[meshcell.connectivity]
        control_point, normal, area = compute_centroid_normal_area(vertices...)
        radius = get_radius(control_point, vertices...)
        panels[i] = Panel(vertices, control_point, normal, SVector{NK,TF}(one(TF) for _ in 1:NK), radius)
    end

    control_points = [panel.control_point for panel in panels]
    normals = [panel.normal for panel in panels]
    strengths = [panel.strength for panel in panels]
    potential = zeros(TF,n_cells)
    velocity = zeros(SVector{3,TF},n_cells)

    wake_points = SVector{3,TF}[]

    return UnstructuredGrid{kernel,TF,NK,NS}(points, meshcells, control_points, normals, strengths, potential, velocity, panels, wake_points)
end

#####
##### structured grid
#####

struct StructuredGrid{TK,TF,NK,NS} <: AbstractPanels{TK,TF,NK,NS}
    # WriteVTK structured grid
    corner_grid::Array{SVector{3,TF},3}

    # additional panel data
    control_points::Array{SVector{3,TF},3}
    normals::Array{SVector{3,TF},3}
    strengths::Array{SVector{NK,TF},3}
    potential::Array{TF,3}
    velocity::Array{SVector{3,TF},3}

    # efficient panel access assembly
    panels::Vector{Panel{TF,NK,NS}}

    # wake shed locations
    wake_points::Vector{SVector{3,TF}}
end

function PanelArray(corner_grid::Array{<:Number,4}, kernel::AbstractKernel; invert_normals=false)
    _, nx, ny, nz = size(corner_grid)
    @assert nz == 1 "corner grid must be a single surface; $nz detected"
    new_corner_grid = Array{SVector{3,eltype(corner_grid)},3}(undef, nx, ny, 1)
    for iy in 1:ny
        for ix in 1:nx
            new_corner_grid[ix,iy,1] = SVector{3,eltype(corner_grid)}(corner_grid[1,ix,iy,1], corner_grid[2,ix,iy,1], corner_grid[3,ix,iy,1])
        end
    end
    return PanelArray(new_corner_grid, kernel; invert_normals)
end

function PanelArray(corner_grid::AbstractArray{SVector{3,TF}}, kernel::AbstractKernel; invert_normals=false) where TF
    # meta parameters
    NK = kernel_multiplicity(kernel)
    NS = 4 # quad panels
    nxp1, nyp1, _ = size(corner_grid)
    nx = nxp1-1
    ny = nyp1-1

    # normal vectors defined implicitly by the grid chirality
    if invert_normals # tranpose corner grid
        new_corner_grid = Array{SVector{3,TF}}(undef,nyp1,nxp1,1)
        corner_grid_t = transpose(view(corner_grid,:,:,1))
        for i in eachindex(corner_grid_t)
            new_corner_grid[i] = corner_grid_t[i]
        end
        corner_grid = new_corner_grid
        new_nx = ny
        ny = nx
        nx = new_nx
    end

    # initialize containers
    strengths = Array{SVector{NK,TF}}(undef,nx,ny,1)
    for i in eachindex(strengths) # unit strengths
        strengths[i] = SVector{NK,TF}(1.0 for _ in 1:NK)
    end
    normals = zeros(SVector{3,TF},nx,ny,1)
    control_points = zeros(SVector{3,TF},nx,ny,1)
    potential = zeros(TF,nx,ny,1)
    velocity = zeros(SVector{3,TF},nx,ny,1)
    # velocity_gradient = zeros(SMatrix{3,3,TF,9},nx,ny,1)

    # update centroids and normals and create fast-access panels
    i = 1
    panels = Vector{Panel{TF,NK,NS}}(undef,nx*ny)
    for iy in 1:ny
        for ix in 1:nx
            # counter-clockwise array access
            v_corner_1 = corner_grid[ix, iy, 1]
            v_corner_2 = corner_grid[ix+1, iy, 1]
            v_corner_3 = corner_grid[ix+1, iy+1, 1]
            v_corner_4 = corner_grid[ix, iy+1, 1]

            # get centroid and normals
            control_point, normal, area = compute_centroid_normal_area(v_corner_1, v_corner_2, v_corner_3, v_corner_4)
            control_points[ix,iy,1] = control_point
            normals[ix,iy,1] = normal

            radius = get_radius(control_point, v_corner_1, v_corner_2, v_corner_3, v_corner_4)

            # create fast-access panel
            panels[i] = Panel{TF,NK,NS}(
                SVector{NS}(v_corner_1, v_corner_2, v_corner_3, v_corner_4),
                control_point, normal, strengths[ix,iy,1], radius
            )
            i += 1
        end
    end

    wake_points = SVector{3,TF}[]

    # return panels
    return StructuredGrid{kernel,TF,NK,NS}(corner_grid, control_points, normals, strengths, potential, velocity, panels, wake_points)
end

function panels_2_vector_strengths!(strengths, panels::AbstractVector{Panel{TF,NK,NS}}) where {TF,NK,NS}
    n_panels = length(panels)
    @assert length(strengths) == length(panels) * NK "size of strengths vector inconsistent with panels"
    strengths_reshaped = reshape(strengths, NK, n_panels)
    for (i,panel) in enumerate(panels)
        strengths_reshaped[:,i] .= panel.strength
    end
    return nothing
end

function vector_2_panels_strengths!(panels::AbstractVector{Panel{TF,NK,NS}}, strengths; panel_indices=1:length(panels), relaxation=0) where {TF,NK,NS}
    n_panels = length(panels)
    @assert length(strengths) == length(panels) * NK "size of strengths vector inconsistent with panels"
    strengths_reshaped = reshape(strengths, NK, n_panels)
    for i in panel_indices
        # (; vertices, control_point, normal, radius) = panels[i]
        vertices = panels[i].vertices
        control_point = panels[i].control_point
        normal = panels[i].normal
        radius = panels[i].radius
        old_strength = panels[i].strength
        new_strength = SVector{NK,TF}(strengths_reshaped[j,i] for j in 1:NK)
        panels[i] = Panel(vertices, control_point, normal, new_strength*(1-relaxation) + old_strength*relaxation, radius)
    end
    return nothing
end

function set_unit_strength!(panels::AbstractPanels{<:Any,TF,<:Any,<:Any}; panel_indices=1:length(panels.panels)) where TF
    for i in panel_indices
        panels.strengths[i] = SVector{1,TF}(1.0)
    end
    set_unit_strength!(panels.panels)
end

function set_unit_strength!(panels::AbstractVector{Panel{TF,1,NS}}; panel_indices=1:length(panels)) where {TF,NS}
    for i in panel_indices
        # (; vertices, control_point, normal, radius) = panels[i]
        vertices = panels[i].vertices
        control_point = panels[i].control_point
        normal = panels[i].normal
        radius = panels[i].radius

        new_strength = SVector{1,TF}(1.0)
        panels[i] = Panel(vertices, control_point, normal, new_strength, radius)
    end
end

function reset_potential_velocity!(panels)
    for i in eachindex(panels.velocity)
        panels.velocity[i] = zero(eltype(panels.velocity))
        panels.potential[i] = zero(eltype(panels.potential))
    end
end

function grid_2_panels_strength!(panel_array::AbstractPanels{<:Any,TF,NK,NS}; panel_indices=1:length(panel_array.panels)) where {TF,NK,NS}
    # (; strengths, panels) = panel_array
    strengths = panel_array.strengths
    panels = panel_array.panels
    for i_panel in panel_indices
        # current panel values
        # (; vertices, control_point, normal, radius) = panels[i_panel]
        vertices = panels[i_panel].vertices
        control_point = panels[i_panel].control_point
        normal = panels[i_panel].normal
        radius = panels[i_panel].radius

        # create fast-access panel
        panels[i_panel] = Panel{TF,NK,NS}(
            vertices, control_point, normal, strengths[i_panel], radius
        )
    end
end

function panels_2_grid_strength!(panel_array::AbstractPanels{<:Any,TF,NK,NS}; panel_indices=1:length(panel_array.panels)) where {TF,NK,NS}

    strengths = panel_array.strengths
    panels = panel_array.panels
    for i_panel in panel_indices
        # current panel values
        strength = panels[i_panel].strength

        # update grid
        strengths[i_panel] = strength
    end

end

#=
function structured_grid_corner_indices(centroid_index::Int, nx::Int)
    icolumn, irow = divrem(centroid_index, nx)
    if irow == 0
        icolumn -= 1
        irow = nx
    end

    # counterclockwise
    corneri1 = icolumn*(nx+1) + irow
    corneri2 = corneri1 + 1
    corneri3 = corneri2 + nx + 1
    corneri4 = corneri3 - 1

    #     # clockwise
    #     corneri1 = icolumn*(nx+1) + irow
    #     corneri2 = corneri1 + nx + 1
    #     corneri3 = corneri2 + 1
    #     corneri4 = corneri1 + 1

    return corneri1, corneri2, corneri3, corneri4
end

function get_vertices(corner_grid::Array{<:SVector,3}, i)
    vertex_1_index, vertex_2_index, vertex_3_index, vertex_4_index = structured_grid_corner_indices(i, size(corner_grid,1)-1)
    return panels.corner_grid[indices[1]], panels.corner_grid[indices[2]], panels.corner_grid[indices[3]], panels.corner_grid[indices[4]]
end

function structured_grid_2_panels!(panels, corner_grid::Array{<:SVector,3}, control_points::Array{<:SVector,3}, normals::Array{<:SVector,3}, strengths::Array{<:Number,3})
    for i in eachindex(panels)
        vertices = get_vertices()
    end
end
=#
