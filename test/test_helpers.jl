using Test
using LinearAlgebra: norm, dot, I
import FLOWPanel as pnl
import GeometricTools as gt
import FastMultipole

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
    1 1 2 2 1 1 2 2;
    3 5 3 5 4 6 4 6;
    5 3 6 6 3 5 3 4;
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
    body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_OCT), copy(CELLS_OCT))
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
    pnl.influence!(target_body, source_body, backend; scalar_potential=true, velocity=true)
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

function translate_nodes!(nodes, vector=SVector(0.0, 0.0, 0.0))
    nodes .+= vector
    return nodes
end