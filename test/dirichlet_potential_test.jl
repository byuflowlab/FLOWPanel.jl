using Test
import FLOWPanel as pnl
using LinearAlgebra: I

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

@testset "PanelWake satisfies exterior Green identity" begin
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
        shedding;
        CPoffset=1e-6,
        ensure_winding=false,
        semiinfinite_wake=false
    )
    body.Das[1] .= repeat(Vinf / pnl.norm(Vinf), 1, size(body.Das[1], 2))

    wake = pnl.PanelWake(body; nwakerows=7)
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
    pnl.calc_controlpoints!(body; off=abs(body.CPoffset))
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

    pnl.calc_controlpoints!(body; off=-abs(body.CPoffset))
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

    wake_phi_green_interior = (I - G_interior) \ phi_source_interior

    @test wake.nwakes[] == 5
    @test maximum(abs.(wake_phi_exterior)) > 0
    @test maximum(abs.(phi_source)) > 0
    @test isapprox(centered(wake_phi_green), centered(wake_phi_exterior); atol=1e-4, rtol=7e-2)
    @test maximum(abs.(wake_phi_interior)) > 0
    @test maximum(abs.(phi_source_interior)) > 0
    @test isapprox(centered(wake_phi_green_interior), centered(wake_phi_interior); atol=1e-4, rtol=7e-2)
end
