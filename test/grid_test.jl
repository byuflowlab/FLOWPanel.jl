function create_sphere_structured(kernel;
        radius=1.0, n_phi=20, n_theta=10, d_theta=5*pi/180,
        invert_normals=true,
    )
    corner_grid = zeros(SVector{3,Float64},n_phi+1,n_theta+1,1)
    for (itheta,theta) in enumerate(range(d_theta,stop=pi-d_theta,length=n_theta+1))
        for (iphi,phi) in enumerate(range(0,stop=2*pi,length=n_phi+1))
            corner_grid[iphi,itheta] = SVector{3}(
                radius * sin(theta) * cos(phi),
                radius * sin(theta) * sin(phi),
                radius * cos(theta)
            )
        end
    end

    return PanelArray(corner_grid, kernel; invert_normals)
end

function create_sphere_unstructured(kernel;
        radius=1.0, n_phi=4, n_theta=3,
        invert_normals=true,
    )
    n_points = 2 + (n_theta-1) * n_phi
    n_mesh_cells = n_phi*2 + n_phi*(n_theta-2)*2
    points = Vector{SVector{3,Float64}}(undef,n_points)
    meshcells = Vector{WriteVTK.MeshCell{WriteVTK.VTKCellType, SVector{3,Int64}}}(undef,n_mesh_cells)

    d_phi = 2*pi/n_phi
    ipoint_theta = 0
    for (itheta,theta) in enumerate(range(0,stop=pi,length=n_theta+1))
        phi_range = itheta == 1 || itheta == n_theta+1 ? (0.0:0.0) : range(0,stop=(n_phi-1)*d_phi,length=n_phi)
        stheta, ctheta = sincos(theta)
        for (iphi,phi) in enumerate(phi_range)
            phi -= d_phi / 2 * 0^iseven(itheta)
            sphi, cphi = sincos(phi)
            points[ipoint_theta + iphi] = SVector{3}(stheta*cphi, stheta*sphi, ctheta)
        end
        ipoint_theta += length(phi_range)
    end
    
    # add top
    i_meshcell = 1
    for ipoint in 2:1+n_phi
        ipoint_p1 = ipoint < n_phi+1 ? ipoint + 1 : 2
        meshcells[i_meshcell] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(1,ipoint,ipoint_p1))
        i_meshcell += 1
    end

    # add middle sections
    for itheta in 2:n_theta-1
        for iphi in 1:n_phi
            ipoint1 = n_phi * (itheta-2) + iphi + 1
            if iseven(itheta)
                ipoint2 = ipoint1+n_phi
                ipoint3 = iphi < n_phi ? ipoint2 + 1 : n_phi * (itheta-1) + 2
                meshcells[i_meshcell] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(ipoint1,ipoint2,ipoint3))
                i_meshcell += 1
                ipoint2 = ipoint3
                ipoint3 = iphi < n_phi ? ipoint1+1 : n_phi * (itheta-2) + 2
                meshcells[i_meshcell] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(ipoint1,ipoint2,ipoint3))
                i_meshcell += 1
            else
                ipoint2 = ipoint1 + n_phi
                ipoint3 = iphi < n_phi ? ipoint1 + 1 : n_phi * (itheta-2) + 2
                meshcells[i_meshcell] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(ipoint1,ipoint2,ipoint3))
                i_meshcell += 1
                ipoint1 = ipoint3
                ipoint3 = iphi < n_phi ? ipoint2+1 : n_phi * (itheta-1) + 2
                meshcells[i_meshcell] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(ipoint1,ipoint2,ipoint3))
                i_meshcell += 1
            end
        end
    end

    # add bottom
    point1 = length(points)
    for ipoint in length(points)-1:-1:length(points)-n_phi
        point2 = ipoint
        point3 = ipoint == point1 - n_phi ? point1-1 : ipoint-1
        meshcells[i_meshcell] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, SVector{3}(point1, point2, point3))
        i_meshcell += 1
    end

    # create panels
    sphere = PanelArray(points, meshcells, kernel)
    return sphere
end

@testset "create grids" begin

sphere_structured = create_sphere_structured(ConstantSource(); n_phi=20, n_theta=10, d_theta=5*pi/180)
sphere_unstructured = create_sphere_unstructured(ConstantSource(); n_phi=20, n_theta=18)

#####
##### apply freestream
#####
freestream = SVector{3}(1.0,0,0)
apply_freestream!(sphere_structured, freestream)
apply_freestream!(sphere_unstructured, freestream)

vtk("unstructured_sphere", sphere_unstructured)
vtk("structured_sphere", sphere_structured)

end
