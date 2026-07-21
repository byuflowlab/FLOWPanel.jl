#!/usr/bin/env julia

"""
Convert an ASCII STL-like file whose `outer loop` blocks each contain four
vertices into an ordinary, consistently wound triangle STL.

Usage:

    julia --project=. debug/dirichlet_solve/triangulate_quad_stl.jl INPUT OUTPUT
"""

import GeoIO
import Meshes
import Printf: @printf

const Point3 = NTuple{3,Float64}
const Triangle = NTuple{3,Int}

subtract(a::Point3, b::Point3) = ntuple(i -> a[i] - b[i], 3)
cross3(a::Point3, b::Point3) = (
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1],
)
dot3(a::Point3, b::Point3) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

function parse_quad_loops(path::AbstractString)
    loops = Vector{NTuple{4,Point3}}()
    current = Point3[]
    in_loop = false

    for (line_number, raw_line) in enumerate(eachline(path))
        words = split(strip(raw_line))
        isempty(words) && continue
        keyword = lowercase(words[1])

        if keyword == "outer" && length(words) == 2 && lowercase(words[2]) == "loop"
            in_loop && error("nested outer loop at $path:$line_number")
            empty!(current)
            in_loop = true
        elseif keyword == "vertex"
            in_loop || error("vertex outside an outer loop at $path:$line_number")
            length(words) == 4 || error("invalid vertex at $path:$line_number")
            xyz = try
                tuple(parse.(Float64, words[2:4])...)
            catch
                error("invalid vertex coordinates at $path:$line_number")
            end
            all(isfinite, xyz) || error("non-finite vertex at $path:$line_number")
            push!(current, xyz)
        elseif keyword == "endloop"
            in_loop || error("endloop without outer loop at $path:$line_number")
            length(words) == 1 || error("invalid endloop at $path:$line_number")
            length(current) == 4 || error(
                "expected four vertices in loop ending at $path:$line_number; " *
                "found $(length(current))")
            push!(loops, tuple(current...))
            in_loop = false
        end
    end

    in_loop && error("unterminated outer loop in $path")
    isempty(loops) && error("no outer loops found in $path")
    return loops
end

function triangulate(loops)
    point_ids = Dict{Point3,Int}()
    points = Point3[]
    triangles = Triangle[]
    removed = 0

    point_id(point) = get!(point_ids, point) do
        push!(points, point)
        length(points)
    end

    for quad in loops
        ids = point_id.(quad)
        for local_ids in ((1, 2, 3), (1, 3, 4))
            tri = (ids[local_ids[1]], ids[local_ids[2]], ids[local_ids[3]])
            a, b, c = points[tri[1]], points[tri[2]], points[tri[3]]
            normal = cross3(subtract(b, a), subtract(c, a))
            if dot3(normal, normal) == 0.0
                removed += 1
            else
                push!(triangles, tri)
            end
        end
    end

    # Discard vertices that appeared only in removed triangles and compact IDs.
    used = sort!(unique(Iterators.flatten(triangles)))
    old_to_new = zeros(Int, length(points))
    for (new, old) in enumerate(used)
        old_to_new[old] = new
    end
    points = points[used]
    triangles = [map(i -> old_to_new[i], tri) for tri in triangles]
    return points, triangles, removed
end

edge_key(a::Int, b::Int) = minmax(a, b)
edge_direction(a::Int, b::Int) = a < b ? 1 : -1

function orient_outward!(points, triangles)
    edge_uses = Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
    for (triangle_id, (a, b, c)) in enumerate(triangles)
        for (u, v) in ((a, b), (b, c), (c, a))
            push!(get!(edge_uses, edge_key(u, v), Tuple{Int,Int}[]),
                (triangle_id, edge_direction(u, v)))
        end
    end

    bad_edges = count(uses -> length(uses) != 2, values(edge_uses))
    bad_edges == 0 || error(
        "mesh is not watertight: $bad_edges edges do not have exactly two triangles")

    neighbors = [Tuple{Int,Int,Int}[] for _ in triangles]
    for uses in values(edge_uses)
        (first_id, first_dir), (second_id, second_dir) = uses
        push!(neighbors[first_id], (second_id, first_dir, second_dir))
        push!(neighbors[second_id], (first_id, second_dir, first_dir))
    end

    orientation = zeros(Int, length(triangles))
    components = Vector{Vector{Int}}()
    for seed in eachindex(triangles)
        orientation[seed] != 0 && continue
        orientation[seed] = 1
        component = Int[]
        queue = [seed]
        next_in_queue = 1
        while next_in_queue <= length(queue)
            triangle_id = queue[next_in_queue]
            next_in_queue += 1
            push!(component, triangle_id)
            for (neighbor, this_dir, neighbor_dir) in neighbors[triangle_id]
                required = -orientation[triangle_id] * this_dir * neighbor_dir
                if orientation[neighbor] == 0
                    orientation[neighbor] = required
                    push!(queue, neighbor)
                elseif orientation[neighbor] != required
                    error("mesh is non-orientable or has inconsistent connectivity")
                end
            end
        end
        push!(components, component)
    end

    for triangle_id in eachindex(triangles)
        if orientation[triangle_id] == -1
            a, b, c = triangles[triangle_id]
            triangles[triangle_id] = (a, c, b)
        end
    end

    for component in components
        signed_volume6 = 0.0
        reference = points[triangles[first(component)][1]]
        for triangle_id in component
            a, b, c = (points[i] for i in triangles[triangle_id])
            signed_volume6 += dot3(subtract(a, reference),
                cross3(subtract(b, reference), subtract(c, reference)))
        end
        signed_volume6 == 0.0 && error("closed mesh component has zero signed volume")
        if signed_volume6 < 0.0
            for triangle_id in component
                a, b, c = triangles[triangle_id]
                triangles[triangle_id] = (a, c, b)
            end
        end
    end
    return triangles
end

function coordinate_bounds(points)
    lower = ntuple(axis -> minimum(point[axis] for point in points), 3)
    upper = ntuple(axis -> maximum(point[axis] for point in points), 3)
    return lower, upper
end

function write_ascii_stl(path, points, triangles)
    solid_name = replace(splitext(basename(path))[1], r"\s+" => "_")
    open(path, "w") do io
        println(io, "solid ", solid_name)
        for (a_id, b_id, c_id) in triangles
            a, b, c = points[a_id], points[b_id], points[c_id]
            normal = cross3(subtract(b, a), subtract(c, a))
            normal_length = sqrt(dot3(normal, normal))
            normal = ntuple(i -> normal[i] / normal_length, 3)
            @printf(io, "  facet normal %.17g %.17g %.17g\n", normal...)
            println(io, "    outer loop")
            for point in (a, b, c)
                @printf(io, "      vertex %.17g %.17g %.17g\n", point...)
            end
            println(io, "    endloop")
            println(io, "  endfacet")
        end
        println(io, "endsolid ", solid_name)
    end
end

function load_with_geoio(path)
    geometry = GeoIO.load(path).geometry
    mesh_vertices = collect(Meshes.vertices(geometry))
    points = Point3[tuple(Float64.(Meshes.coordinates(vertex))...) for vertex in mesh_vertices]
    connectivity = getfield(Meshes.topology(geometry), :connec)
    triangles = Triangle[]
    for cell in connectivity
        ids = Tuple(Int.(collect(getfield(cell, :indices))))
        length(ids) == 3 || error("GeoIO reloaded a non-triangle cell")
        push!(triangles, ids)
    end
    return points, triangles
end

function validate_mesh(points, triangles; expected_bounds=nothing)
    isempty(points) && error("mesh has no vertices")
    isempty(triangles) && error("mesh has no triangles")
    length(unique(points)) == length(points) || error("mesh contains duplicate vertices")
    all(tri -> length(unique(tri)) == 3, triangles) || error("mesh has repeated triangle vertices")
    for tri in triangles
        a, b, c = (points[i] for i in tri)
        normal = cross3(subtract(b, a), subtract(c, a))
        dot3(normal, normal) > 0.0 ||
            error("mesh contains a zero-area triangle")
    end
    expected_bounds === nothing || coordinate_bounds(points) == expected_bounds ||
        error("coordinate bounds changed during conversion")
    oriented = copy(triangles)
    orient_outward!(points, oriented)
    oriented == triangles || error("mesh winding is not consistently outward")
    return true
end

function convert_quad_stl(input_path::AbstractString, output_path::AbstractString)
    loops = parse_quad_loops(input_path)
    input_points = unique(Iterators.flatten(loops))
    bounds = coordinate_bounds(input_points)
    points, triangles, removed = triangulate(loops)
    orient_outward!(points, triangles)
    validate_mesh(points, triangles; expected_bounds=bounds)
    write_ascii_stl(output_path, points, triangles)

    loaded_points, loaded_triangles = load_with_geoio(output_path)
    validate_mesh(loaded_points, loaded_triangles; expected_bounds=bounds)
    length(loaded_triangles) == length(triangles) ||
        error("GeoIO triangle count changed on reload")
    length(loaded_points) == length(points) ||
        error("GeoIO vertex count changed on reload")

    @printf("input loops: %d\n", length(loops))
    @printf("output triangles: %d\n", length(triangles))
    @printf("unique vertices: %d\n", length(points))
    @printf("removed degeneracies: %d\n", removed)
    @printf("bounds: (%.17g, %.17g, %.17g) to (%.17g, %.17g, %.17g)\n",
        bounds[1]..., bounds[2]...)
    println("validation: passed (GeoIO reload, watertight, outward, nondegenerate)")
    return (; loops=length(loops), triangles=length(triangles), vertices=length(points),
        removed, bounds)
end

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: julia --project=. " *
        "debug/dirichlet_solve/triangulate_quad_stl.jl INPUT OUTPUT")
    convert_quad_stl(args[1], args[2])
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
