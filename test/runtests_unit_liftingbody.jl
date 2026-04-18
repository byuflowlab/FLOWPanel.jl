using Test
import FLOWPanel as pnl
import GeometricTools as gt
using StaticArrays: SVector

include("test_helpers.jl")

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

@testset verbose=true "RigidWakeBody" begin
    @testset "construction" begin
        body = make_plate_vortex_body()
        @test body.nnodes == 4
        @test body.ncells == 2
        @test length(body.shedding) == 1
        @test size(body.Das[1]) == (3, size(body.shedding[1], 2) + 1)
        @test size(body.strength) == (body.ncells, 1)
    end

    @testset "shedding edge consistency" begin
        body = make_plate_vortex_body()
        @test size(body.shedding[1], 1) == 6
        @test size(body.shedding_full) == (6, body.ncells)
        te_panels = findall(!=(-1), vec(body.shedding_full[3, :]))
        @test !isempty(te_panels)
        @test all(body.shedding_full[4, te_panels] .> 0)
    end

    @testset "extra_reset!" begin
        body = make_plate_vortex_body()
        for vte in body.velocity_te
            vte .= 3.0
        end
        pnl.reset!(body)
        @test all(all(vte .== 0) for vte in body.velocity_te)
    end

    @testset "seeded shedding trace" begin
        nodes, cells = make_seeded_te_mesh()
        bbox = ([0.8, -0.1, -0.1], [1.1, 2.1, 0.1])

        trace = pnl.trace_trailing_edge(nodes, cells, 1, 2; bbox=bbox, end_node=3)
        @test trace.nodes == [1, 2, 3]
        @test length(trace.edges) == 2

        shedding = pnl.calc_shedding_from_seed(nodes, cells, 1, 2; bbox=bbox, end_node=3)
        expected = Int[
            3 4;
            3 3;
            2 2;
            1 2;
            3 3;
            2 2;
        ]
        @test shedding == expected

        bbox_svec = [SVector(0.8, -0.1, -0.1), SVector(1.1, 2.1, 0.1)]
        @test pnl.calc_shedding_from_seed(nodes, cells, 1, 2; bbox=bbox_svec, end_node=3) == expected

        body = pnl.RigidWakeBody{pnl.VortexRing}(nodes, cells, shedding; check_mesh=false, watertight=false)
        @test size(body.shedding[1]) == (6, 2)
        @test size(body.Das[1]) == (3, 3)
    end

    @testset "seeded shedding validation" begin
        nodes, cells = make_seeded_te_mesh()
        @test_throws ErrorException pnl.calc_shedding_from_seed(nodes, cells, 4, 2)
        @test_throws ErrorException pnl.calc_shedding_from_seed(nodes, cells, 1, 2; end_node=8)
    end

    @testset "seeded shedding relaxed bootstrap" begin
        bbox = ([0.0, -0.1, -0.1], [2.0, 0.1, 0.1])

        nodes, cells = make_relaxed_seed_mesh(-1.0)
        trace = pnl.trace_trailing_edge(nodes, cells, 1, 2; bbox=bbox, end_node=3)
        @test trace.nodes == [1, 2, 3]
        @test length(trace.edges) == 2
        @test trace.edges[1].normal_jump < 0.2
        @test trace.edges[2].normal_jump >= 0.2

        nodes_strict, cells_strict = make_relaxed_seed_mesh(0.0)
        @test_throws ErrorException pnl.trace_trailing_edge(nodes_strict, cells_strict, 1, 2; bbox=bbox, end_node=3)
    end

    @testset "seeded shedding backslashcoupled" begin
        nodes, cells = make_seeded_te_mesh()
        bbox = ([0.8, -0.1, -0.1], [1.1, 2.1, 0.1])

        trace = pnl.trace_trailing_edge(nodes, cells, 1, 2; bbox=bbox, end_node=3)
        @test trace.nodes == [1, 2, 3]
        @test length(trace.edges) == 2

        shedding = pnl.calc_shedding_from_seed(nodes, cells, 1, 2; bbox=bbox, end_node=3)
        expected = Int[
            3 4;
            3 3;
            2 2;
            1 2;
            3 3;
            2 2;
        ]
        @test shedding == expected

        bbox_svec = [SVector(0.8, -0.1, -0.1), SVector(1.1, 2.1, 0.1)]
        @test pnl.calc_shedding_from_seed(nodes, cells, 1, 2; bbox=bbox_svec, end_node=3) == expected

        body = pnl.RigidWakeBody{Union{<:pnl.ConstantSource, <:pnl.ConstantDoublet}}(nodes, cells, shedding; check_mesh=false, watertight=false)
        
        magVinf = 10.0
        AOA = 0.0
        Vinf = magVinf * [cosd(AOA), sind(AOA), 0.0]
        pnl.apply_freestream!(body, Vinf)
        
        for i in eachindex(body.Das)
            body.Das[i] .= repeat(Vinf/magVinf, 1, size(body.Das[i],2))
        end
    
        backend = pnl.DirectBackend()
        solver = pnl.BackslashCoupled((body,))
        @test size(solver.G) == (body.ncells, body.ncells)
        println("Solving body...")

        bodies = (body,)
        # benchmarks(out_file, bodies, solver; backend)
        pnl.solve!((body,), solver; backend, update_G=true)
        for (i, body) in enumerate(bodies)
            println("Strength column 1:")
            println("  max = ", maximum(bodies[i].strength[:, 1]))
            println("  min = ", minimum(bodies[i].strength[:, 1]))

            println("Strength column 2:")
            println("  max = ", maximum(bodies[i].strength[:, 2]))
            println("  min = ", minimum(bodies[i].strength[:, 2]))
        end
    end
end
