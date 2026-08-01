using Test
using LinearAlgebra: norm, dot, I
import FLOWPanel as pnl
import FastMultipole
import GeoIO
import Meshes
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
    shedding = [pnl.calc_shedding_from_seed(nodes, cells, 1, 3)]
    body = pnl.RigidWakeBody{pnl.VortexRing}(nodes, cells, shedding; DBC=false, check_mesh=false, watertight=false)
    body.Das[1] .= repeat([1.0, 0.0, 0.0], 1, size(body.Das[1], 2))
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.strength[:, 1] .= 1.0
    return body
end

"""
    make_dirichlet_diamond_body(; nspan=1, thick=0.06, das=0.3)

Small diamond-airfoil wedge (2 chordwise panels per side, `nspan` span cells,
sharp TE with shared nodes) as a fully paired Dirichlet source+doublet
`RigidWakeBody` with a finite attached wake. The chordwise resolution gives
the trailing-edge circulation genuine leverage on the paired centroid
pressures, which the BRAINSTORM-015 pressure-continuity Kutta tests require
(a flat plate with a coplanar wake is degenerate: zero leverage).
"""
function make_dirichlet_diamond_body(; nspan::Int=1, thick=0.06, das=0.3)
    ys = range(0, 1; length=nspan+1)
    nodes = Float64[]
    # node layout per span station: LE, mid-upper, TE, mid-lower
    for y in ys
        append!(nodes, [0.0, y, 0.0])
        append!(nodes, [0.5, y, thick])
        append!(nodes, [1.0, y, 0.0])
        append!(nodes, [0.5, y, -thick])
    end
    nodes = reshape(nodes, 3, :)
    idx(j, k) = (j-1)*4 + k
    cells = Int[]
    for j in 1:nspan
        le1, up1, te1, lo1 = idx(j, 1), idx(j, 2), idx(j, 3), idx(j, 4)
        le2, up2, te2, lo2 = idx(j+1, 1), idx(j+1, 2), idx(j+1, 3), idx(j+1, 4)
        append!(cells, [le1, up1, up2]); append!(cells, [le1, up2, le2])
        append!(cells, [up1, te1, te2]); append!(cells, [up1, te2, up2])
        append!(cells, [le1, le2, lo2]); append!(cells, [le1, lo2, lo1])
        append!(cells, [lo1, lo2, te2]); append!(cells, [lo1, te2, te1])
    end
    cells = reshape(cells, 3, :)
    shedding = [pnl.calc_shedding_from_seed(nodes, cells, idx(1, 3), idx(2, 3))]
    body = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.VortexRing}}(
        nodes, cells, shedding; check_mesh=false, watertight=false,
        semiinfinite_wake=false)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat([das, 0.0, 0.0], 1, size(body.Das[i], 2))
    end
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    return body
end

const EXTERNAL_TARGETS = Float64[
    2.0  -1.5   0.0   0.5   1.25;
    0.2   0.4  -1.2   1.0  -0.75;
    1.1   1.3   1.5   2.0   1.7;
]

function make_sphere_source_body(; radius=1.0, ntheta=24, nphi=48, theta_pad=0.15)
    theta = collect(range(theta_pad, pi - theta_pad; length=ntheta + 1))
    phi = collect(range(0, 2pi; length=nphi + 1))[1:end-1]
    nodes = zeros(Float64, 3, length(theta) * length(phi))
    lin = LinearIndices((length(theta), length(phi)))
    for ti in eachindex(theta), pi_i in eachindex(phi)
        th = theta[ti]
        ph = phi[pi_i]
        nodes[:, lin[ti, pi_i]] .= radius .* [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)]
    end

    cells = zeros(Int, 3, 2 * ntheta * nphi)
    ci = 0
    for ti in 1:ntheta, pi_i in 1:nphi
        pn = pi_i == nphi ? 1 : pi_i + 1
        n11 = lin[ti, pi_i]
        n12 = lin[ti, pn]
        n21 = lin[ti + 1, pi_i]
        n22 = lin[ti + 1, pn]
        cells[:, ci += 1] .= (n11, n21, n22)
        cells[:, ci += 1] .= (n11, n22, n12)
    end

    return pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells; watertight=true)
end

function solve_source_body!(body; uinf=[1.0, 0.0, 0.0], rho=1.0, backend=pnl.DirectBackend())
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
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
    Ps = zeros(body.ncells)
    Fs = zeros(3, body.ncells)
    areas = pnl.calc_areas(body)
    pnl.calcfield_P!(Ps, body, body.velocity, norm(uinf), rho, nothing;
        correct_kuttacondition=false)
    pnl.calcfield_F!(Fs, body, areas, body.normals, Ps; correct_kuttacondition=false)
    return body, Ps, Fs
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

make_basic_triangle_surface() = (copy(NODES_2TRI), copy(CELLS_2TRI))
make_octa_triangle_surface() = (copy(NODES_OCT), copy(CELLS_OCT))

function meshes_to_nodes_cells(mesh)
    vertices = collect(Meshes.vertices(mesh))
    nodes = reduce(hcat, (collect(Meshes.coordinates(v)) for v in vertices))
    connec = getfield(Meshes.topology(mesh), :connec)
    cells = reduce(hcat, (collect(getfield(c, :indices)) for c in connec))
    return Float64.(nodes), Int.(cells)
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

    nodes, cells = meshes_to_nodes_cells(msh)

    # Generate TE shedding matrix
    shedding = zeros(Int, 6, 0)

    # Generate the paneled body
    body = bodytype(nodes, cells, [shedding]; flip_normals=false, watertight=true)
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
    pnl.calcfield_P!(bodies, magVinf, rho; correct_kuttacondition=fill(true, length(bodies)))
    for body in bodies
        pnl.calcfield_Cp(body, magVinf, rho)
    end
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
        pnl.calc_controlpoints!(body)
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
    pnl.calc_controlpoints!(body)
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

has_dirichlet_body(body::pnl.AbstractBody{<:Any, <:Any, <:Any, DBC}) where DBC = DBC

"""
Evaluate the solved Neumann boundary condition by adding the induced velocity
field to the current external velocity field and projecting onto panel normals.
"""
function flow_tangency_max_residuals(bodies::Tuple; backend=pnl.DirectBackend(), cp_off=nothing)
    Uext = [copy(body.velocity) for body in bodies]

    for body in bodies
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
    end

    pnl.influence!(bodies, bodies, backend; scalar_potential=false, velocity=true)

    res = zeros(length(bodies))
    for (i, body) in enumerate(bodies)
        r = 0.0
        for (vel, normal) in zip(eachcol(body.velocity), eachcol(body.normals))
            r = max(r, abs(dot(vel, normal)))
        end
        res[i] = r
    end

    for (body, velocity) in zip(bodies, Uext)
        body.velocity .= velocity
    end

    return res
end

"""
Evaluate the solved Dirichlet boundary condition by recomputing the interior
perturbation potential induced by the solved source/doublet strengths.
"""
function interior_potential_max_residuals(bodies::Tuple; backend=pnl.DirectBackend(), cp_off=nothing)
    phi_ext = [copy(body.potential) for body in bodies]

    for body in bodies
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.potential .= 0.0
    end

    pnl.influence!(bodies, bodies, backend; scalar_potential=true, velocity=false)

    res = [maximum(abs.(body.potential .+ potential)) for (body, potential) in zip(bodies, phi_ext)]

    for (body, potential) in zip(bodies, phi_ext)
        body.potential .= potential
    end

    return res
end

function assert_boundary_residuals(
    bodies::Tuple;
    backend=pnl.DirectBackend(),
    tangency_atol=1e-10,
    potential_atol=1e-10,
)
    tangency_residuals = flow_tangency_max_residuals(bodies; backend)
    potential_residuals = interior_potential_max_residuals(bodies; backend)

    for (i, body) in enumerate(bodies)
        if has_dirichlet_body(body)
            @test potential_residuals[i] < potential_atol
        else
            @test tangency_residuals[i] < tangency_atol
        end
    end

    return (; tangency_residuals, potential_residuals)
end

function flow_potential_residuals(bodies::Tuple; cp_off=nothing)
    for body in bodies
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
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

"""
    make_conversion_fixture(; nwakerows, wrap, max_particles=2000)

Deterministic `PanelParticleWake` whose panel-wake rows are written directly,
with **no time stepping**, so a single panel-to-particle conversion can be
exercised in isolation (BRAINSTORM 016).

The wake carries three shedding columns. `wrap=false` lays the sheet flat in
the x-y plane with rows at constant `x`; `wrap=true` closes the chain into a
triangular ring by making the last node column an exact copy of the first, so
the conversion's wrap detection fires and no root/tip closure is deposited.
Strengths vary in both the streamwise and spanwise directions, which keeps
every trailing, unsteady, and root/tip contribution distinct and nonzero.

`nwakes[]` is set to `nwakerows`, i.e. the buffer is full and the outgoing row
is the final active row -- the state in which `shed_wake!` converts.
"""
function make_conversion_fixture(; nwakerows::Int, wrap::Bool, max_particles=2000)
    body = make_dirichlet_diamond_body(; nspan=3)
    wake = pnl.PanelParticleWake(body; nwakerows=nwakerows, max_particles=max_particles)
    pw = wake.panel_wake
    nodes = pw.nodes[1]
    strength = pw.strength[1]
    n_node_rows = size(nodes, 2)
    n_node_cols = size(nodes, 3)

    for irow in 1:n_node_rows, icol in 1:n_node_cols
        if wrap
            theta = 2pi * (icol - 1) / (n_node_cols - 1)
            nodes[1, irow, icol] = cos(theta)
            nodes[2, irow, icol] = sin(theta)
            nodes[3, irow, icol] = 0.25 * (irow - 1)
        else
            nodes[1, irow, icol] = 0.5 * (irow - 1)
            nodes[2, irow, icol] = icol - 1
            nodes[3, irow, icol] = 0.0
        end
    end
    # Close the ring exactly; the 5*eps() wrap test must not depend on how
    # accurately cos/sin round-trips at 2pi.
    wrap && (nodes[:, :, n_node_cols] .= nodes[:, :, 1])

    for irow in 1:size(strength, 2), icol in 1:size(strength, 3)
        strength[1, irow, icol] = 0.1 * irow + 0.3 * icol
    end

    pw.nwakes[] = nwakerows
    return wake
end
