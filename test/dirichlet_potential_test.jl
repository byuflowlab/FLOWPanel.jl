using Test
import FLOWPanel as pnl
using LinearAlgebra: I, rank

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

# @testset "PanelWake satisfies exterior Green identity" begin
    backend = pnl.DirectBackend()

    AOA = 5.0
    c = 2.0
    b = 8.0
    Vinf = [cosd(AOA), 0.0, sind(AOA)]

    meshfile = joinpath(pnl.examples_path, "data", "wing_ar4_naca0016_5.msh")
    msh = GeoIO.load(meshfile).geometry |> Meshes.Scale(1.0)
    nodes, cells = pnl.meshes2nodes_cells(msh)

    trailingedge = zeros(3, 10000)
    trailingedge[1, :] .= c
    trailingedge[2, :] .= range(-b / 2, stop=b / 2, length=size(trailingedge, 2))
    shedding = pnl.calc_shedding(nodes, cells, trailingedge; tolerance=0.001 * b)

    body = pnl.RigidWakeBody{Union{pnl.ConstantSource,pnl.ConstantDoublet}}(
        nodes,
        cells,
        # pnl.noshedding;
        shedding;
        cp_outer=true,
        ensure_winding=false,
        semiinfinite_wake=false
    )
    body.Das[1] .= repeat((Vinf .+ [0.0, 0.0, 1.0]) / pnl.norm(Vinf), 1, size(body.Das[1], 2))

    wake = pnl.PanelWake(body; nwakerows=7)
    # wake = nothing
    frames = pnl.ReferenceFrame(body)
    solver = pnl.Backslash(body)
    maneuver = (frames, systems, wakes, t) -> nothing

    mktempdir() do simpath
        pnl.simulate!(body, wake, frames, maneuver, t -> Vinf, collect(0.0:0.05:0.25);
            body_solvers=solver,
            backend=backend,
            path=simpath,
        )
    end

    pnl.calc_normals!(body)
    body.cp_outer = true
    pnl.calc_controlpoints!(body)
    body.potential .= 0.0
    pnl.influence!(body, pnl.get_sources(wake), backend; scalar_potential=true, velocity=false)
    wake_phi_exterior = copy(body.potential)

    body.velocity .= 0.0
    pnl.influence!(body, pnl.get_sources(wake), backend; scalar_potential=false, velocity=true)
    wake_velocity_exterior = copy(body.velocity)

    body.velocity .= wake_velocity_exterior
    pnl.set_strengths!(body)
    body.strength[:, 2] .= 0.0
    body.potential .= 0.0
    pnl.influence!(body, body, backend; scalar_potential=true, velocity=false)
    phi_source = copy(body.potential)

    G = zeros(body.ncells, body.ncells)
    pnl._G!(G, body, body; kerneloffset=body.kerneloffset, update_geometry=false)
    green_matrix = -G
    wake_phi_green = green_matrix \ phi_source

    centered(v) = v .- sum(v) / length(v)

    body.cp_outer = false
    pnl.calc_controlpoints!(body)
    body.potential .= 0.0
    pnl.influence!(body, pnl.get_sources(wake), backend; scalar_potential=true, velocity=false)
    wake_phi_interior = copy(body.potential)

    body.velocity .= 0.0
    pnl.influence!(body, pnl.get_sources(wake), backend; scalar_potential=false, velocity=true)
    wake_velocity_interior = copy(body.velocity)

    G_interior = zeros(body.ncells, body.ncells)
    pnl._G!(G_interior, body, body; kerneloffset=body.kerneloffset, update_geometry=false)

    body.velocity .= wake_velocity_interior
    pnl.set_strengths!(body)
    body.strength[:, 2] .= 0.0
    body.potential .= 0.0
    pnl.influence!(body, body, backend; scalar_potential=true, velocity=false)
    phi_source_interior = copy(body.potential)

    wake_phi_green_interior = G_interior \ phi_source_interior
    # wake_phi_green_interior = (I - G_interior) \ phi_source_interior

    @test wake.nwakes[] == 5
    @test maximum(abs.(wake_phi_exterior)) > 0
    @test maximum(abs.(phi_source)) > 0
    # @test isapprox(centered(wake_phi_green), centered(wake_phi_exterior); atol=1e-4, rtol=7e-2)
    @test maximum(abs.(wake_phi_interior)) > 0
    @test maximum(abs.(phi_source_interior)) > 0
    # @test isapprox(centered(wake_phi_green_interior), centered(wake_phi_interior); atol=1e-4, rtol=7e-2)
# end

# NOTE: an influence matrix constructed with exterior control points is not full rank
@show rank(green_matrix), size(green_matrix)[1]
@show rank(I-G_interior), size(I-G_interior)[1]

# NOTE: G_interior == G + I
@show maximum(abs.(G_interior - (G + I)))
@show maximum(abs.(green_matrix - (I - G_interior)))

# new test
body.cp_outer = false # interior
pnl.calc_normals!(body)
pnl.calc_controlpoints!(body)
G_interior2 = zeros(body.ncells, body.ncells)
pnl._G!(G_interior2, body, body)
@show rank(G_interior2), size(G_interior2)[1]

# new test
body.velocity .= wake_velocity_interior
solver = pnl.Backslash(body)
pnl.solve!(body, solver; backend)
phi_sum = phi_source_interior - body.strength[:, 2]

#--- new new test ---#

# remove kutta wake panels for simplicity
body2 = pnl.RigidWakeBody{Union{pnl.ConstantSource,pnl.ConstantDoublet}}(
        nodes,
        cells,
        pnl.noshedding;
        # shedding;
        cp_outer=false,
        ensure_winding=false,
        semiinfinite_wake=false
    )
pnl.calc_normals!(body2)
pnl.calc_controlpoints!(body2)

# wake-induced potential and velocity at body control points
body2.velocity .= 0.0
body2.potential .= 0.0
pnl.influence!(body2, pnl.get_sources(wake), backend; scalar_potential=true, velocity=true)

# choose source and doublet strengths
body2.strength[:, 1] .= vec(sum(body2.normals .* .-body2.velocity, dims=1))  # source strength = normal velocity
body2.strength[:, 2] .= body2.potential  # doublet strength = potential
potential_ext = copy(body2.potential)
velocity_ext = copy(body2.velocity)
body2.velocity .= 0.0
body2.potential .= 0.0
pnl.influence!(body2, body2, backend; scalar_potential=true, velocity=true)
potential_green_int = copy(body2.potential) .* 2

# Indirect exterior-potential reconstruction.
#
# The interior and exterior double-layer limits differ by the potential jump:
# G_interior ≈ G_exterior + I. Therefore the exterior Green operator can be
# formed from interior-collocated influences as A = I - G_interior. This
# operator still has the physical constant-potential nullspace, so add a
# zero-mean gauge row/column correction before solving.
n = body.ncells
mean_gauge = fill(1 / n, n, n)
interior_jump_operator = Matrix(I, n, n) - G_interior
interior_jump_operator_gauged = interior_jump_operator + mean_gauge
phi_exterior_indirect = interior_jump_operator_gauged \ phi_source_interior

@test rank(interior_jump_operator) == n - 1
@test rank(interior_jump_operator_gauged) == n
@test isapprox(centered(phi_exterior_indirect), centered(wake_phi_exterior); atol=1e-4, rtol=7e-2)
