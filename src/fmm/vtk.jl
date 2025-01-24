function vtk(name::AbstractString, unstructured_grid::UnstructuredGrid, time=0.0; save_path="")
    # (; points, meshcells, control_points, normals, strengths, potential, velocity) = unstructured_grid
    points = unstructured_grid.points
    meshcells = unstructured_grid.meshcells
    control_points = unstructured_grid.control_points
    normals = unstructured_grid.normals
    strengths = unstructured_grid.strengths
    potential = unstructured_grid.potential
    velocity = unstructured_grid.velocity
    vtk_grid(joinpath(save_path,name), points, meshcells) do vtk
        vtk["control points"] = control_points
        vtk["normals"] = normals
        vtk["strengths"] = strengths
        vtk["potential"] = potential
        vtk["velocity"] = velocity
        # vtk["time"] = time
    end
end

function vtk(name::AbstractString, structured_grid::StructuredGrid, time=0.0; save_path="")
    # (; corner_grid, control_points, normals, strengths, potential, velocity) = structured_grid
    corner_grid = structured_grid.corner_grid
    control_points = structured_grid.control_points
    normals = structured_grid.normals
    strengths = structured_grid.strengths
    potential = structured_grid.potential
    velocity = structured_grid.velocity
    vtk_grid(joinpath(save_path,name), corner_grid) do vtk
        vtk["control points"] = control_points
        vtk["normals"] = normals
        vtk["strengths"] = strengths
        vtk["potential"] = potential
        vtk["velocity"] = velocity
        # vtk["time"] = time
    end
end