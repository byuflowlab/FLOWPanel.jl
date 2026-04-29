using Test
import FLOWPanel as pnl
import GeometricTools as gt
using LinearAlgebra: diag
import Meshes
using StaticArrays: SVector

@testset verbose=true "Solvers" begin
    @testset "Backslash construction" begin
        nlb = make_octa_source_body()
        rwb = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(
            Float64[
                0 1 1 0;
                0 0 1 1;
                0 0 0 0;
            ],
            Int[
                1 1;
                2 3;
                3 4;
            ];
            check_mesh=false,
            watertight=false,
        )
        @test pnl.Backslash(nlb) isa pnl.Backslash
        @test pnl.Backslash(rwb) isa pnl.Backslash
    end

    @testset "Backslash Neumann construction and solve" begin
        body = make_octa_source_body()
        solver = pnl.Backslash(body)
        @test size(solver.G) == (body.ncells, body.ncells)
        @test all(diag(solver.G) .!= 0)
        @test solver.Glu !== nothing
        @test length(solver.rhs) == body.ncells

        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        pnl.solve!(body, solver)
        @test any(abs.(body.strength[:, 1]) .> 0)
        max_residual = nonlifting_flow_tangency_max_residual(body)
        @test max_residual < 1e-12
    end

    @testset "Backslash Dirichlet construction and solve" begin
        body = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(
            Float64[
                0 1 1 0;
                0 0 1 1;
                0 0 0 0;
            ],
            Int[
                1 1;
                2 3;
                3 4;
            ],
            [reshape(Int[1, 2, 3, 2, 3, 2], 6, 1)];
            check_mesh=false,
            watertight=false,
        )
        body.Das[1] .= repeat([1.0, 0.0, 0.2], 1, size(body.Das[1], 2)) ./ sqrt(1.0^2 + 0.2^2)
        pnl.calc_normals!(body)
        normals = copy(body.normals)
        solver = pnl.Backslash(body)
        i1, i2, i3 = body.cells[:, 1]
        centroid1 = (body.nodes[:, i1] + body.nodes[:, i2] + body.nodes[:, i3]) ./ 3
        inward = body.controlpoints[:, 1] - centroid1
        @test dot(inward, normals[:, 1]) < 0
        @test size(solver.G) == (body.ncells, body.ncells)
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        body.velocity[3, :] .= 0.2
        pnl.solve!(body, solver)
        @test any(abs.(body.strength[:, 2]) .> 0)
        assert_boundary_residuals((body,); potential_atol=1e-7)
    end

    @testset "Backslash Dirichlet construction and solve" begin
        body = make_nonlifting(Union{pnl.ConstantSource, pnl.ConstantDoublet}; DBC=true)
        pnl.calc_normals!(body)
        normals = copy(body.normals)
        solver = pnl.Backslash(body)
        i1, i2, i3 = body.cells[:, 1]
        centroid1 = (body.nodes[:, i1] + body.nodes[:, i2] + body.nodes[:, i3]) ./ 3
        inward = body.controlpoints[:, 1] - centroid1
        @test dot(inward, normals[:, 1]) < 0
        @test size(solver.G) == (body.ncells, body.ncells)

        body = pnl.NonLiftingBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(copy(NODES_OCT), copy(CELLS_OCT); DBC=true)
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        body.velocity[3, :] .= 0.2
        normals = pnl.calc_normals!(body)
        expected_sigma = -vec(sum(body.velocity .* normals; dims=1))
        solver = pnl.Backslash(body)
        pnl.solve!(body, solver)
        @test any(abs.(body.strength[:, 2]) .> 0)
        @test isapprox(vec(body.strength[:, 1]), expected_sigma; atol=1e-12)
        assert_boundary_residuals((body,); potential_atol=1e-10)
    end

    @testset "KrylovSolver solve" begin
        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0

        solver = pnl.KrylovSolver(body; backend=pnl.DirectBackend(), atol=1e-8, rtol=1e-8, itmax=50)
        pnl.solve!(body, solver)

        @test any(abs.(body.strength[:, 1]) .> 0)
        assert_boundary_residuals((body,); tangency_atol=1e-7)
    end

    @testset "KrylovSolver Dirichlet solve" begin
        body_direct = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(
            Float64[
                0 1 1 0;
                0 0 1 1;
                0 0 0 0;
            ],
            Int[
                1 1;
                2 3;
                3 4;
            ];
            check_mesh=false,
            watertight=false,
        )
        body_krylov = deepcopy(body_direct)

        for body in (body_direct, body_krylov)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
            body.velocity[3, :] .= 0.2
        end

        direct = pnl.Backslash(body_direct)
        pnl.solve!(body_direct, direct; backend=pnl.DirectBackend())

        krylov = pnl.KrylovSolver(body_krylov; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=20)
        pnl.solve!(body_krylov, krylov; backend=pnl.DirectBackend())

        @test isapprox(vec(body_krylov.strength[:, 1]), vec(body_direct.strength[:, 1]); atol=1e-12)
        @test isapprox(vec(body_krylov.strength[:, 2]), vec(body_direct.strength[:, 2]); atol=1e-12)
        assert_boundary_residuals((body_direct,); potential_atol=1e-7)
        assert_boundary_residuals((body_krylov,); potential_atol=1e-7)
    end

    @testset "KrylovSolver + JacobiPreconditioner" begin
        make_dirichlet_body() = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(
            Float64[
                0 1 1 0;
                0 0 1 1;
                0 0 0 0;
            ],
            Int[
                1 1;
                2 3;
                3 4;
            ];
            check_mesh=false,
            watertight=false,
        )

        body_no_pc = make_dirichlet_body()
        body_pc = make_dirichlet_body()

        for body in (body_no_pc, body_pc)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
            body.velocity[3, :] .= 0.2
        end

        solver_no_pc = pnl.KrylovSolver(body_no_pc; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50)
        solver_pc = pnl.KrylovSolver(body_pc; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50, preconditioner_cell_size=0.5)

        @test solver_pc.preconditioner !== nothing
        @test length(solver_pc.preconditioner.cell_body_indices) > 0

        pnl.solve!(body_no_pc, solver_no_pc)
        pnl.solve!(body_pc, solver_pc)

        @test isapprox(vec(body_pc.strength[:, 1]), vec(body_no_pc.strength[:, 1]); atol=1e-6)
        @test isapprox(vec(body_pc.strength[:, 2]), vec(body_no_pc.strength[:, 2]); atol=1e-6)
        assert_boundary_residuals((body_no_pc,); potential_atol=1e-7)
        assert_boundary_residuals((body_pc,); potential_atol=1e-6)
    end

    @testset "FGSSolver construction" begin
        body1 = make_octa_source_body()
        body2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        solver1 = pnl.FGSSolver(body1)
        solver2 = pnl.FGSSolver((body1, body2))
        @test solver1.max_iterations == 100
        @test solver1.tolerance == 1e-6
        @test solver2.max_iterations == 100
        @test solver2.tolerance == 1e-6
    end

    @testset "Multi-body solve" begin
        body1 = make_octa_source_body()
        body2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        body1.velocity .= 0
        body2.velocity .= 0
        body1.velocity[1, :] .= 1.0
        body2.velocity[1, :] .= 1.0
        solver1 = pnl.Backslash(body1)
        solver2 = pnl.Backslash(body2)
        pnl.solve!((body1, body2), (solver1, solver2); backend=(pnl.DirectBackend(), pnl.DirectBackend()))
        @test any(abs.(body1.strength[:, 1]) .> 0)
        @test any(abs.(body2.strength[:, 1]) .> 0)
        assert_boundary_residuals((body1, body2); tangency_atol=1e-7)
    end

    @testset "numtype" begin
        body_source = make_nonlifting(pnl.ConstantSource)
        @test pnl.numtype(body_source) == Float64
        body32 = pnl.NonLiftingBody{pnl.ConstantSource}(Float32.(NODES_2TRI), copy(CELLS_2TRI))
        @test pnl.numtype(body32) == Float32
    end

    @testset "BackslashCoupled nonlifting" begin
        body1 = make_octa_source_body()
        body2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        magVinf = 10.0
        AOA = 0.0
        Vinf = magVinf * [cosd(AOA), sind(AOA), 0.0]
        pnl.apply_freestream!(body1, Vinf)
        pnl.apply_freestream!(body2, Vinf)

        bodies = (body1, body2)
        backend = pnl.DirectBackend()
        solver = pnl.BackslashCoupled(bodies)
        nps = sum(body.ncells for body in bodies)
        @test size(solver.G) == (nps, nps)

        pnl.solve!(bodies, solver; backend, update_G=true)
        @test any(abs.(body1.strength[:, 1]) .> 0)
        @test any(abs.(body2.strength[:, 1]) .> 0)
        assert_boundary_residuals(bodies; backend, tangency_atol=1e-8)
    end

    @testset "backslashcoupled" begin
        nodes, cells = make_seeded_te_mesh()
        bbox = ([0.8, -0.1, -0.1], [1.1, 2.1, 0.1])
        trace = pnl.trace_trailing_edge(nodes, cells, 1, 2; bbox=bbox, end_node=3)

        @test trace.nodes == [1, 2, 3]
        @test length(trace.edges) == 2
        
        nodes2 = translate_nodes!(copy(nodes), SVector(2.0, 0.0, 0.0))
        bbox2 = ([0.8, -0.1, -0.1], [5.1, 6.1, 4.1])
        trace2 = pnl.trace_trailing_edge(nodes2, cells, 1, 2; bbox=bbox2, end_node=3, debug=true)
        
        @test length(trace2.edges) == 2
        @test trace2.nodes == [1, 2, 3]

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
        # @show typeof(shedding)
        # @show shedding

        bbox_svec = [SVector(0.8, -0.1, -0.1), SVector(1.1, 2.1, 0.1)]
        @test pnl.calc_shedding_from_seed(nodes, cells, 1, 2; bbox=bbox_svec, end_node=3) == expected

        shedding2 = pnl.calc_shedding_from_seed(nodes2, cells, 1, 2; bbox=bbox2, end_node=3)
        bbox_svec2 = [SVector(0.8, -0.1, -0.1), SVector(5.1, 6.1, 4.1)]

        body = pnl.RigidWakeBody{Union{<:pnl.ConstantSource, <:pnl.ConstantDoublet}}(nodes, cells, shedding; check_mesh=false, watertight=true)
        body2 = pnl.RigidWakeBody{Union{<:pnl.ConstantSource, <:pnl.ConstantDoublet}}(nodes2, cells, shedding2; check_mesh=false, watertight=false)
        # body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells)
        # body2 = pnl.NonLiftingBody{pnl.ConstantSource}(nodes2, cells)

        magVinf = 10.0
        AOA = 0.0
        Vinf = magVinf * [cosd(AOA), sind(AOA), 0.0]
        pnl.apply_freestream!(body, Vinf)
        pnl.apply_freestream!(body2, Vinf)

        for i in eachindex(body.Das)
            body.Das[i] .= repeat(Vinf/magVinf, 1, size(body.Das[i],2))
        end
    
        for i in eachindex(body2.Das)
            body2.Das[i] .= repeat(Vinf/magVinf, 1, size(body2.Das[i],2))
        end

        bodies = (body, body2)
        backend = pnl.DirectBackend()
        solver = pnl.BackslashCoupled(bodies)
        nps = sum(b.ncells for b in bodies)
        @test size(solver.G) == (nps, nps)
        println("Solving body...")


        # benchmarks(out_file, bodies, solver; backend)
        pnl.solve!(bodies, solver; backend, update_G=true)
        for body in bodies
            body.velocity .= 0.0
        end
        pnl.apply_freestream!(body, Vinf)
        pnl.apply_freestream!(body2, Vinf)
        assert_boundary_residuals(bodies; backend, potential_atol=1e-4)
    end

    @testset "backslash" begin
        nodes, cells = make_seeded_te_mesh()
        bbox = ([0.8, -0.1, -0.1], [1.1, 2.1, 0.1])
        trace = pnl.trace_trailing_edge(nodes, cells, 1, 2; bbox=bbox, end_node=3)
        chord, span = get_chord_span(nodes)

        @test trace.nodes == [1, 2, 3]
        @test length(trace.edges) == 2
        
        nodes2 = translate_nodes!(copy(nodes), SVector(2.0, 0.0, 0.0))
        chord2, span2 = get_chord_span(nodes2)
        bbox2 = ([0.8, -0.1, -0.1], [5.1, 6.1, 4.1])
        trace2 = pnl.trace_trailing_edge(nodes2, cells, 1, 2; bbox=bbox2, end_node=3, debug=true)
        
        @test length(trace2.edges) == 2
        @test trace2.nodes == [1, 2, 3]

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

        shedding2 = pnl.calc_shedding_from_seed(nodes2, cells, 1, 2; bbox=bbox2, end_node=3)
        bbox_svec2 = [SVector(0.8, -0.1, -0.1), SVector(5.1, 6.1, 4.1)]

        body = pnl.RigidWakeBody{Union{<:pnl.ConstantSource, <:pnl.ConstantDoublet}}(nodes, cells, shedding; check_mesh=false, watertight=false)
        # body2 = pnl.RigidWakeBody{Union{<:pnl.ConstantSource, <:pnl.ConstantDoublet}}(nodes2, cells, shedding2; check_mesh=false, watertight=false)
        # body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells)
        body2 = pnl.NonLiftingBody{pnl.ConstantSource}(nodes2, cells)

        magVinf = 10.0
        AOA = 0.0
        Vinf = magVinf * [cosd(AOA), sind(AOA), 0.0]
        chords = [chord, chord2]
        pnl.apply_freestream!(body, Vinf)
        pnl.apply_freestream!(body2, Vinf)

        for i in eachindex(body.Das)
            body.Das[i] .= repeat(Vinf/magVinf, 1, size(body.Das[i],2))
        end
    
        # for i in eachindex(body2.Das)
        #     body2.Das[i] .= repeat(Vinf/magVinf, 1, size(body2.Das[i],2))
        # end

        # bodies = (body, body2)
        bodies = (body, body2)
        backend = pnl.DirectBackend()
        solver = (pnl.Backslash(body), pnl.Backslash(body2))
        nps = sum(b.ncells for b in bodies)
        println("Solving body...")

        # benchmarks(out_file, bodies, solver; backend)
        pnl.solve!(bodies, solver)
        # for (i, body) in enumerate(bodies) 
        #     for j in 1:size(body.strength, 2)
        #         println("Strength column $(j):")
        #         println("  max = ", maximum(bodies[i].strength[:, j]))
        #         println("  min = ", minimum(bodies[i].strength[:, j]))
        #     end
        # end
        for body in bodies
            body.velocity .= 0.0
        end
        pnl.apply_freestream!(body, Vinf)
        pnl.apply_freestream!(body2, Vinf)
        assert_boundary_residuals(bodies; backend, tangency_atol=1e-7, potential_atol=1e-5)

        # @show CL, CD = postprocess!(bodies, Vinf, rho, backend, chords, span)
    end

    @testset "backslash_meshes" begin
        import GeoIO

        run_names = ["nasa_wing.msh", "nasa_surface_spaced_repaired.msh"]
        # run_names = ["nasa_surface_spaced.msh"]
        file_path       = "examples"
        paraview        = true                      # Whether to visualize with Paraview
        out_file = joinpath(pnl.examples_path, "wing_aileron", "coupled_timing_results.csv")

        files = [joinpath(pnl.examples_path, "wing_aileron", name) for name in run_names]
        nodes1 = [42, 43]
        nodes2 = [34, 35]
        nodes1 = [42, 19]
        nodes2 = [34, 3]


        # ----------------- SIMULATION PARAMETERS --------------------------------------
        m              = 0.0254
        AOA             = 0.0                      # (deg) freestream angle of attack
        magVinf         = 117.3 * m * 12                      # (m/s) freestream velocity
        rho             = 1.225                     # (kg/m^3) air density

        # ----------------- GEOMETRY DESCRIPTION ---------------------------------------
        c_body1 = 10 * m
        b = 60 * m                            # (m) span length
        c_body2 = 2
        AR_body1 = b / c_body1                             # (m) span length
        AR_body2 = b / c_body2                             # (m) span length

        chords = [c_body1, c_body2]
        ARs = [AR_body1, AR_body2] 
        Sref = b * (c_body1 + c_body2)

        scaling = 1.0
        # ----------------- SOLVER SETTINGS -------------------------------------------

        # Body and wake model
        kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}               # Kernel type to use
        # body type
        bodytype = pnl.RigidWakeBody{kernel}

        # Processing
        clip_Cp         = 1 - 342.0/magVinf         # Clip pressure coefficients that are lower than this threshold
    
        Vinf = magVinf * [cosd(AOA), sind(AOA), 0.0]

        bodies = tuple([generate_body(file, chord, b, bodytype, scaling, 1, Vinf, firstnode, secondnode)
                        for (file, chord, firstnode, secondnode) in zip(files, chords, nodes1, nodes2)]...)

        # bodies = (bodies[2],)
        #------------------- SOLVE BODY ----------------------------------------------
        backend = pnl.DirectBackend()
        solver = pnl.BackslashCoupled(bodies)
        println("Solving body...")

        pnl.solve!(bodies, solver; backend, update_G=true)

        for body in bodies
            body.velocity .= 0.0
            pnl.apply_freestream!(body, Vinf)
        end

        assert_boundary_residuals(bodies; backend, potential_atol=1e-4)

        # write vtk files
        for i in eachindex(bodies)
            pnl.write_vtk("check_mesh_body_$(i)", bodies[i])
        end
    end

end

println("done.")
