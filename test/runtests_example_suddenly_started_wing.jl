using Test
import FLOWPanel as pnl
using LinearAlgebra: norm

if !isdefined(@__MODULE__, :SSWConfig)
    include(joinpath(pnl.examples_path, "suddenly_started_wing.jl"))
end

function _ssw_boundary_edges(cells)
    counts = Dict{Tuple{Int,Int},Int}()
    for cell in eachcol(cells)
        for (a, b) in ((cell[1], cell[2]), (cell[2], cell[3]), (cell[3], cell[1]))
            edge = a < b ? (a, b) : (b, a)
            counts[edge] = get(counts, edge, 0) + 1
        end
    end
    return [edge for (edge, count) in counts if count == 1]
end

@testset "Suddenly-started wing geometry and formulation" begin
    contour = naca0012_contour(21)
    @test size(contour, 2) == 2
    @test minimum(contour[:, 1]) ≈ 0.0 atol=100eps()
    @test maximum(contour[:, 1]) ≈ 1.0 atol=100eps()
    @test maximum(contour[:, 2]) ≈ -minimum(contour[:, 2]) rtol=1e-13
    @test contour[1, 1] ≈ 1.0 atol=100eps()
    @test contour[1, 2] ≈ 0.0 atol=1e-12
    @test count(i -> isapprox(contour[i, 1], 0.0; atol=100eps()), axes(contour, 1)) == 1

    config = SSWConfig(n_span=2, n_airfoil=21, save_vtk=false,
        backend_kind=:direct, verbose=false)
    nodes, cells = suddenly_started_wing_mesh(config.c, config.AR * config.c;
        n_span=config.n_span, n_airfoil=config.n_airfoil)
    n_section = size(contour, 1)
    @test size(nodes) == (3, n_section * (config.n_span + 1))
    @test size(cells) == (3, 2 * n_section * config.n_span)
    @test collect(extrema(nodes[2, :])) ≈
        [-config.AR * config.c / 2, config.AR * config.c / 2]

    watertight, open_cells = pnl.iswatertight(nodes, cells; return_open_cells=true)
    @test !watertight
    @test !isempty(open_cells)
    boundary_edges = _ssw_boundary_edges(cells)
    @test length(boundary_edges) == 2 * n_section
    @test all(edge -> begin
        y = nodes[2, collect(edge)]
        all(isapprox.(abs.(y), config.AR * config.c / 2; atol=100eps()))
    end, boundary_edges)

    body = build_suddenly_started_wing(config)
    @test body isa pnl.RigidWakeBody{pnl.VortexRing,1,Float64,false}
    @test !pnl.has_dirichlet_bc(body)
    @test !body.watertight
    @test !body.semiinfinite_wake
    @test length(body.shedding) == 1
    @test size(only(body.shedding), 2) == config.n_span
    @test body.nsheddings == config.n_span

    fmm_config = ssw_with(config; backend_kind=:fmm,
        fmm_expansion_order=8, fmm_multipole_acceptance=0.25,
        fmm_leaf_size=64)
    backend = _ssw_backend(fmm_config)
    @test backend.expansion_order == 8
    @test backend.multipole_acceptance == 0.25
    @test backend.leaf_size == 64

    shedding = only(body.shedding)
    te_midpoints = map(eachcol(shedding)) do edge
        panel, nia, nib = edge[1:3]
        node_ids = body.cells[[nia, nib], panel]
        (body.nodes[:, node_ids[1]] + body.nodes[:, node_ids[2]]) / 2
    end
    @test all(p -> isapprox(p[1], config.c; atol=1e-8 * config.c), te_midpoints)
    @test issorted(getindex.(te_midpoints, 2)) || issorted(getindex.(te_midpoints, 2); rev=true)
end

@testset "Suddenly-started wing time, inflow, Wagner function, and wake capacity" begin
    config = SSWConfig(n_span=1, n_airfoil=21, t_end_star=0.25,
        dt_star=0.125, save_vtk=false, backend_kind=:direct, verbose=false)
    time = collect(ssw_time_range(config))
    @test time ≈ [0.0, 0.125, 0.25] .* config.c ./ config.U
    @test diff(time) ≈ fill(config.dt_star * config.c / config.U, 2)

    sim = prepare_suddenly_started_wing(config)
    @test sim.Uinf(0.0) == zero(sim.Uinf(0.0))
    post_start = sim.Uinf(time[2])
    @test norm(post_start) ≈ config.U
    @test post_start[1] ≈ config.U * cosd(config.alpha_deg)
    @test post_start[2] == 0
    @test abs(post_start[3]) ≈ config.U * sind(config.alpha_deg)

    @test ssw_wagner(0.0) ≈ 0.5
    @test ssw_wagner(1.0) ≈ 1 - 0.165exp(-0.09) - 0.335exp(-0.6)
    @test ssw_wagner(Inf) == 1.0
    @test issorted(ssw_wagner.(0.0:0.25:7.0))

    @test length(sim.t_range) == length(time)
    @test sim.t_range ≈ time
    @test sim.wake isa pnl.PanelWake
    @test size(only(sim.wake.nodes), 2) - 1 == length(time) - 1
    @test !sim.wake.include_final_filament
    @test sim.wake.core_size ≈ config.c * config.wake_core_over_c
    @test !sim.wake.overflowed[]
    @test pnl.monitor_provides(sim.pressure) == (:P,)
    @test pnl.monitor_requires(sim.force) == (:P,)
    @test pnl.audit_monitors((sim.pressure, sim.force)) !== nothing
end

@testset "Suddenly-started wing tiny unsteady smoke" begin
    config = SSWConfig(n_span=2, n_airfoil=21, t_end_star=0.125,
        dt_star=0.125, save_vtk=false, backend_kind=:direct, verbose=false,
        backslash_max_panels=10_000)
    sim = prepare_suddenly_started_wing(config)

    pnl.simulate!((sim.wing,), (sim.wake,), sim.frames, sim.maneuver!, sim.Uinf, sim.t_range;
        body_solvers=(sim.solver,),
        backend=sim.backend,
        monitors=sim.monitors,
        path=nothing,
        set_Das_eta_freestream=NaN,
        grad_mu_options=(; basis=:tri),
        verbose=false,
    )

    @test all(isfinite, sim.force.force[:, 2])
    @test norm(sim.force.force[:, 2]) > 0
    @test sim.wake.nwakes[] == 1
    @test !sim.wake.overflowed[]
end
