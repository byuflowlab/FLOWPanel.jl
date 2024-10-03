#####
##### general geometric helper functions
#####
@inline function compute_centroid(p1,p2,p3)
    return (p1+p2+p3) ./ 3
end

function compute_centroid_normal_area(points...)
    n_points = length(points)
    p1 = points[1]
    centroid = @SVector zeros(eltype(p1),3)
    normal = @SVector zeros(eltype(p1),3)
    area = zero(eltype(p1))
    for i_triangle in 2:n_points-1
        v1 = points[i_triangle] - p1
        v2 = points[i_triangle+1] - p1
        this_cross = cross(v1,v2)
        this_cross_norm = norm(this_cross)
        this_area = this_cross_norm/2
        centroid = centroid + this_area * compute_centroid(points[1], points[i_triangle], points[i_triangle+1])
        normal = normal + this_cross * this_area
        area += this_area
    end
    return centroid ./ area, normal ./ norm(normal), area
end

#####
##### generators for common or test geometries
#####
function generate_cube(side_length, centroid, kernel)
    p1x = centroid[1] - side_length/2
    p1y = centroid[2] - side_length/2
    p1z = centroid[3] - side_length/2
    points = [
        SVector{3}(p1x,p1y,p1z),
        SVector{3}(p1x+side_length,p1y,p1z),
        SVector{3}(p1x+side_length,p1y+side_length,p1z),
        SVector{3}(p1x,p1y+side_length,p1z),
        SVector{3}(p1x,p1y,p1z+side_length),
        SVector{3}(p1x+side_length,p1y,p1z+side_length),
        SVector{3}(p1x+side_length,p1y+side_length,p1z+side_length),
        SVector{3}(p1x,p1y+side_length,p1z+side_length),
    ]
    cells = [
        MeshCell(VTKCellTypes.VTK_QUAD,[1,4,3,2]),
        MeshCell(VTKCellTypes.VTK_QUAD,[1,2,6,5]),
        MeshCell(VTKCellTypes.VTK_QUAD,[5,6,7,8]),
        MeshCell(VTKCellTypes.VTK_QUAD,[5,8,4,1]),
        MeshCell(VTKCellTypes.VTK_QUAD,[4,8,7,3]),
        MeshCell(VTKCellTypes.VTK_QUAD,[3,7,6,2]),
    ]
    return PanelArray(points, cells, kernel)
end

function generate_plane()

end