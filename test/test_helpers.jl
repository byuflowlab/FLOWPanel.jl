using Test
using LinearAlgebra: norm, dot, I
import FLOWPanel as pnl
import GeometricTools as gt
import FastMultipole
using StaticArrays: SVector

const NODES_1TRI = Float64[
    0 1 0;
    0 0 1;
    0 0 0;
]
const CELLS_1TRI = reshape(Int[1, 2, 3], 3, 1)

const NODES_2TRI = Float64[
    0 1 1 0;
    0 0 0 0;
    0 0 1 1;
]
const CELLS_2TRI = Int[
    1 1;
    2 3;
    3 4;
]

const NODES_OCT = Float64[
     1 -1  0  0  0  0;
     0  0  1 -1  0  0;
     0  0  0  0  1 -1;
]
const CELLS_OCT = Int[
    1 3 2 4 3 2 4 1;
    3 2 4 1 1 3 2 4;
    5 5 5 5 6 6 6 6;
]

function exact_linear_system(n::Int; seed_shift::Int=0)
    A = Matrix{Float64}(I, n, n)
    for i in 1:n, j in 1:n
        A[i, j] += 0.05 * sin(i + 3j + seed_shift)
    end
    y_exact = collect(range(0.25, 1.25; length=n))
    b = A * y_exact
    return A, y_exact, b
end

function random_spd_system(n::Int)
    M = zeros(Float64, n, n)
    for i in 1:n, j in 1:n
        M[i, j] = sin(0.3 * i + 0.7 * j) + cos(0.5 * i - 0.2 * j)
    end
    A = M' * M + n * I
    y_exact = [sin(0.1 * i) + 0.01 * i for i in 1:n]
    b = A * y_exact
    return Matrix(A), y_exact, b
end

make_nonlifting(::Type{E}, nodes=NODES_2TRI, cells=CELLS_2TRI; kwargs...) where {E} =
    pnl.NonLiftingBody{E}(copy(nodes), copy(cells); kwargs...)

function make_octa_source_body()
    body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_OCT), copy(CELLS_OCT); 
                CPoffset=1e-14,
                kerneloffset=1e-12)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.strength[:, 1] .= 1.0
    return body
end

function make_plate_vortex_body()
    nodes = Float64[
        0 1 1 0;
        0 0 1 1;
        0 0 0 0;
    ]
    cells = Int[
        1 1;
        2 3;
        3 4;
    ]
    shedding = [reshape(Int[1, 2, 3, 2, 3, 2], 6, 1)]
    body = pnl.RigidWakeBody{pnl.VortexRing}(nodes, cells, shedding; check_mesh=false, watertight=false)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.strength[:, 1] .= 1.0
    return body
end

const EXTERNAL_TARGETS = Float64[
    2.0  -1.5   0.0   0.5   1.25;
    0.2   0.4  -1.2   1.0  -0.75;
    1.1   1.3   1.5   2.0   1.7;
]

function make_sphere_source_body(; radius=1.0, ntheta=24, nphi=48, theta_pad=0.15)
    p_min = [theta_pad, 0.0, 0.0]
    p_max = [pi - theta_pad, 2pi, 0.0]
    ndivs = [ntheta, nphi, 0]
    grid = gt.Grid(p_min, p_max, ndivs; loop_dim=2)
    gt.transform!(grid, X -> gt.spherical3D(vcat(radius, X[1:2])))
    tri = gt.GridTriangleSurface(grid, 1)
    return pnl.NonLiftingBody{pnl.ConstantSource}(tri)
end

function solve_source_body!(body; uinf=[1.0, 0.0, 0.0], rho=1.0, backend=pnl.DirectBackend())
    body.velocity .= repeat(uinf, 1, body.ncells)
    solver = pnl.Backslash(body)
    pnl.solve!(body, solver)

    tgt = deepcopy(body)
    tgt.velocity .= 0
    tgt.potential .= 0
    pnl.influence!(tgt, body, backend; velocity=true)
    for d in 1:3
        tgt.velocity[d, :] .+= uinf[d]
    end

    body.velocity .= tgt.velocity
    pnl.calcfield_Cp!(body, norm(uinf))
    pnl.calcfield_F!(body, norm(uinf), rho)
    return body
end

function translated_nonlifting_target(shift)
    body = make_nonlifting(pnl.ConstantSource, NODES_OCT .+ shift, CELLS_OCT)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    return body
end

function translated_rigid_target(shift)
    body = make_plate_vortex_body()
    body.nodes .+= shift
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.velocity .= 0
    body.potential .= 0
    return body
end

function evaluate_velocity_and_potential(target_body, source_body, backend)
    target_body.velocity .= 0
    target_body.potential .= 0
    scalar_potential = !FastMultipole.has_vector_potential(source_body)
    pnl.influence!(target_body, source_body, backend; scalar_potential=scalar_potential, velocity=true)
    return copy(target_body.velocity), copy(target_body.potential)
end

function make_basic_triangle_surface()
    vertices = [
        gt.Meshes.Point(0.0, 0.0, 0.0),
        gt.Meshes.Point(1.0, 0.0, 0.0),
        gt.Meshes.Point(1.0, 0.0, 1.0),
        gt.Meshes.Point(0.0, 0.0, 1.0),
    ]
    triangles = [gt.Meshes.connect((1, 2, 3)), gt.Meshes.connect((1, 3, 4))]
    mesh = gt.Meshes.SimpleMesh(vertices, triangles)
    return gt.GridTriangleSurface(mesh)
end

function make_octa_triangle_surface()
    vertices = [gt.Meshes.Point(NODES_OCT[1, i], NODES_OCT[2, i], NODES_OCT[3, i]) for i in axes(NODES_OCT, 2)]
    triangles = [gt.Meshes.connect(Tuple(CELLS_OCT[:, i])) for i in axes(CELLS_OCT, 2)]
    mesh = gt.Meshes.SimpleMesh(vertices, triangles)
    return gt.GridTriangleSurface(mesh)
end

function translate_nodes!(nodes, vector=SVector(0.0, 0.0, 0.0))
    nodes .+= vector
    return nodes
end

function generate_body(
    meshfile::String,
    chord::Float64,
    span::Float64,
    bodytype::Type{<:pnl.RigidWakeBody},
    scaling::Float64 = 1.0,
    flip::Int64 = 1,
    Vinf::AbstractVector{<:Real} = zeros(3),
    firstnode=-1,
    secondnode=-1
)
    magVinf = norm(Vinf)

    # Read Gmsh mesh
    msh = GeoIO.load(meshfile).geometry

    # Transform the mesh: scale
    msh = msh |> Meshes.Scale(scaling)

    # Wrap into Grid object
    grid = pnl.gt.GridTriangleSurface(msh)

    # Generate TE shedding matrix
    shedding = zeros(Int, 6, 0)

    # Generate the paneled body
    CPoffset = 1e-6
    body = bodytype(grid, [shedding]; CPoffset, flip_normals=false)
    # pnl.write_vtk("spaced_nasa", body)

    # Recompute shedding from the finalized cell winding used by `body`.
    if firstnode == -1 || secondnode == -1
        @warn "firstnode and secondnode not provided; TE shedding will be disabled. This may cause inaccurate results for lifting bodies."
    else
        shedding = pnl.calc_shedding_from_seed(
            body.nodes,
            body.cells,
            firstnode, secondnode
        )
        body = bodytype(
            body.nodes,
            body.cells,
            [shedding];
            CPoffset,
            flip_normals=false,
            ensure_winding=false
        )
    end

    # initialize wake doublets
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(Vinf/magVinf, 1, size(body.Das[i],2))
    end

    pnl.apply_freestream!(body, Vinf)

    return body
end

function make_seeded_te_mesh()
    nodes = Float64[
        1 1 1 0 0 0 0 0;
        0 1 2 0 1 2 0 2;
        0 0 0 1 1 1 -1 -1;
    ]
    cells = Int[
        4 5 7 8;
        2 3 1 2;
        1 2 2 3;
    ]
    return nodes, cells
end

function make_relaxed_seed_mesh(second_edge_z)
    nodes = Float64[
        0 1 2 0 1 0 1;
        0 0 0 1 1 -1 -1;
        0 0 0 0 0 0 second_edge_z;
    ]
    cells = Int[
        4 6 5 7;
        1 2 2 3;
        2 1 3 2;
    ]
    return nodes, cells
end

function postprocess!(bodies, Vinf, rho, backend, chords, span)
    Dhat = Vinf / norm(Vinf)
    Shat = [0, 1, 0]
    Lhat = cross(Dhat, Shat)
    Sref = 0.0
    for chord in chords
        Sref += chord * span
    end

    pnl.calcfield_U!(bodies, Vinf; backend)
    pnl.apply_freestream!(bodies, Vinf)
    pnl.calcfield_Cp!(bodies, magVinf; correct_kuttacondition=fill(true, length(bodies)))
    pnl.calcfield_F!(bodies, magVinf, rho)
    LDS = pnl.calcfield_LDS!(zeros(3,3), bodies, Lhat, Dhat, cross(Lhat, Dhat))

    # Force coefficients
    nondim = 0.5 * rho * norm(Vinf)^2 * Sref
    CL = sign(dot(LDS[:,1], Lhat)) * norm(LDS[:,1]) / nondim
    CD = sign(dot(LDS[:,2], Dhat)) * norm(LDS[:,2]) / nondim

    return CL, CD
end

function get_chord_span(nodes::AbstractMatrix{<:Real})
    @assert size(nodes, 1) == 3 "nodes must be 3×N (x,y,z)"

    x = @view nodes[1, :]
    y = @view nodes[2, :]
    z = @view nodes[3, :]

    chord = maximum(x) - minimum(x)
    span  = maximum(y) - minimum(y)

    return chord, span
end

"""
Assumes freestream has already been applied.
"""
function flow_tangency_residuals(bodies::Tuple)
    for body in bodies
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body; off=1e-10)
    end

    pnl.influence!(bodies, bodies; scalar_potential=false, velocity=true)

    res = zeros(length(bodies))
    for (i, body) in enumerate(bodies)
        r = 0.0
        for (vel, normal) in zip(eachcol(body.velocity), eachcol(body.normals))
            vx, vy, vz = vel
            nx, ny, nz = normal
            r1 = vx * nx + vy * ny + vz * nz
            r += r1 * r1
        end
        res[i] = r / size(body.normals, 2)
    end

    return res
end

function nonlifting_flow_tangency_max_residual(body::pnl.NonLiftingBody{<:Any, <:Any, <:Any, false})
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body; off=0.0)
    pnl.influence!(body, body, pnl.DirectBackend(); velocity=true)

    r = 0.0
    for (vel, normal) in zip(eachcol(body.velocity), eachcol(body.normals))
        vx, vy, vz = vel
        nx, ny, nz = normal
        r = max(r, abs(vx * nx + vy * ny + vz * nz))
    end

    return r
end

function nonlifting_flow_tangency_max_residual(body::pnl.NonLiftingBody{<:Any, <:Any, <:Any, true})
    pnl.calc_normals!(body)
    return maximum(abs.(vec(body.strength[:, 1]) .+ vec(sum(body.velocity .* body.normals; dims=1))))
end

function flow_potential_residuals(bodies::Tuple; cp_off=-1e-10)
    for body in bodies
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body; off=cp_off)
        body.potential .= 0.0
    end

    pnl.influence!(bodies, bodies; scalar_potential=true, velocity=false)
    res = zeros(length(bodies))
    for (i, body) in enumerate(bodies)
        r = 0.0
        for potential in body.potential
            r += potential * potential
        end
        res[i] = r / length(body.potential)
    end

    return res
end
