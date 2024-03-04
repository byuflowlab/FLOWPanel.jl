struct Quaternion{TF}
    real::TF
    pure::SVector{3,TF}
end

function rotate(v, q::Quaternion)
    t = 2*cross(q.pure, v)
    return v + q.real*t + cross(q.pure,t)
end

"perform the reverse rotation"
function antirotate(v, q::Quaternion)
    q = Quaternion(q.real, -q.pure)
    return rotate(v, q)
end

function get_quaternion(axis, angle)
    @assert isapprox(norm(axis), 1.0; atol=1e-12) "axis must be a unit vector"
    stheta_over_2, ctheta_over_2 = sincos(angle/2)
    real_part = ctheta_over_2
    pure_part = axis * stheta_over_2
    quat = Quaternion(real_part, pure_part)
    return quat
end

function rotate(v, axis, angle)
    return rotate(v, get_quaternion(axis, angle))
end

function antirotate(v, axis, angle)
    return rotate(v, get_quaternion(axis, -angle))
end

function create_panel(;
        vertices = SVector{4}(
            SVector{3}(-0.5,-0.5,0.0),
            SVector{3}(0.5,-0.5,0.0),
            SVector{3}(0.5,0.5,0.0),
            SVector{3}(-0.5,0.5,0.0),
        ),
        control_point = SVector{3}(0.0,0.0,0.0),
        normal = SVector{3}(0,0,1.0),
        strength = SVector{1}(1.0),
        invert_normals = false,
    )

    radius = FLOWPanel.get_radius(control_point, vertices...)
    
    if invert_normals
        vertices = reverse(vertices)
        normal = -1 * normal
    end

    return Panel(vertices, control_point, normal, strength, radius)
end

function create_random_panel(seed=123;
        vertices = SVector{3}(
            rand(SVector{3}),
            rand(SVector{3}),
            rand(SVector{3}),
        ),
        strength = SVector{1}(1.0),
        invert_normals = false,
    )
    control_point = (vertices[1]+vertices[2]+vertices[3])/3
    normal = cross(vertices[2]-vertices[1], vertices[3]-vertices[1])
    normal /= norm(normal)

    radius = FLOWPanel.get_radius(control_point, vertices...)
    
    if invert_normals
        vertices = reverse(vertices)
        normal = -1 * normal
    end

    return Panel(vertices, control_point, normal, strength, radius)
end

function create_tripanel(;
        vertices = SVector{4}(
            rand(SVector{3,Float64}),
            rand(SVector{3,Float64}),
            rand(SVector{3,Float64}),
        ),
        strength = SVector{1}(1.0),
        invert_normals = false,
        seed = 123,
    )
    Random.seed!(seed)
    control_point = (vertices[1] + vertices[2] + vertices[3])/3
    normal = cross(vertices[2]-vertices[1], vertices[3]-vertices[1])
    normal /= norm(normal)

    radius = FLOWPanel.get_radius(control_point, vertices...)
    
    if invert_normals
        vertices = reverse(vertices)
        normal = -1 * normal
    end

    return Panel(vertices, control_point, normal, strength, radius)
end

function test_velocity(target, panel, kernel=ConstantSource())
    this_potential(target) = induced(target, panel, kernel)[1]
    return -ForwardDiff.gradient(this_potential, target)
end

function test_gradient(target, panel, kernel=ConstantSource())
    this_velocity(target) = induced(target, panel, kernel)[2]
    return ForwardDiff.jacobian(this_velocity, target)
end

function create_sphere_structured(kernel, center=SVector{3}(0.0,0.0,0.0);
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
            ) + center
        end
    end

    return PanelArray(corner_grid, kernel; invert_normals)
end

function create_sphere_unstructured(kernel, center=SVector{3}(0.0,0.0,0.0);
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
            points[ipoint_theta + iphi] = SVector{3}(stheta*cphi, stheta*sphi, ctheta) + center
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
