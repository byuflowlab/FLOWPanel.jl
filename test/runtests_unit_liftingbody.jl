using Test
import FLOWPanel as pnl
using StaticArrays: SVector, SMatrix

@testset verbose=true "RigidWakeBody" begin
    @testset "construction" begin
        body = make_plate_vortex_body()
        @test body.nnodes == 4
        @test body.ncells == 2
        @test length(body.shedding) == 1
        @test size(body.Das[1]) == (3, size(body.shedding[1], 2) + 1)
        @test size(body.strength) == (body.ncells, 1)

        nodes = Float64[
            0 1 1 0;
            0 0 1 1;
            0 0 0 0;
        ]
        cells = Int[
            1 1;
            2 4;
            3 3;
        ]
        shedding = [reshape(Int[1, 2, 3, 2, 3, 2], 6, 1)]

        body_fixed = pnl.RigidWakeBody{pnl.VortexRing}(nodes, cells, shedding; check_mesh=false, watertight=false)
        pnl.calc_normals!(body_fixed)
        @test body_fixed.normals[3, 1] > 0
        @test body_fixed.normals[3, 2] > 0

        body_flipped = pnl.RigidWakeBody{pnl.VortexRing}(nodes, cells, shedding; check_mesh=false, watertight=false, flip_normals=true)
        pnl.calc_normals!(body_flipped)
        @test body_flipped.normals[3, 1] < 0
        @test body_flipped.normals[3, 2] < 0

        open_nodes, open_cells = make_basic_triangle_surface()
        open_grid_body = pnl.RigidWakeBody{pnl.VortexRing}(open_nodes, open_cells, zeros(Int, 6, 0); watertight=false)
        @test open_grid_body.watertight == false

        closed_nodes, closed_cells = make_octa_triangle_surface()
        closed_grid_body = pnl.RigidWakeBody{pnl.VortexRing}(closed_nodes, closed_cells, zeros(Int, 6, 0); watertight=true)
        @test closed_grid_body.watertight == true
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

    @testset "initialize_Das! matches simulate inline logic" begin
        body_helper = make_plate_vortex_body()
        body_inline = make_plate_vortex_body()
        systems_helper = (body_helper,)
        systems_inline = (body_inline,)
        frames_helper = pnl.ReferenceFrame(body_helper;
            origin=SVector{3}(0.0, 0.0, 0.0),
            v=SVector{3}(0.0, 0.0, 0.0),
            ω_axis=SVector{3}(0.0, 0.0, 1.0),
            ω=2.0,
            R=SMatrix{3,3}(1.0, 0.0, 0.0,
                           0.0, 1.0, 0.0,
                           0.0, 0.0, 1.0),
            name="test",
            child_index=Int[],
            dependent_index=[1])
        frames_inline = pnl.ReferenceFrame(body_inline;
            origin=SVector{3}(0.0, 0.0, 0.0),
            v=SVector{3}(0.0, 0.0, 0.0),
            ω_axis=SVector{3}(0.0, 0.0, 1.0),
            ω=2.0,
            R=SMatrix{3,3}(1.0, 0.0, 0.0,
                           0.0, 1.0, 0.0,
                           0.0, 0.0, 1.0),
            name="test",
            child_index=Int[],
            dependent_index=[1])
        Uinf(t) = SVector{3}(1.0, 0.25, -0.1)
        t0 = 0.0
        dt0 = 0.05
        eta_freestream = 0.3
        eta_kinematic = 0.2

        # `set_Das_kinematic_arc=false` selects the legacy tangent construction,
        # which is what the inline logic below reproduces. The default is now the
        # arc construction; it is compared separately at the end of this testset.
        pnl.initialize_Das!(systems_helper, frames_helper, Uinf, t0, dt0;
            set_Das_eta_freestream=eta_freestream,
            set_Das_eta_kinematic=eta_kinematic,
            set_Das_kinematic_arc=false)

        uinf0 = Uinf(t0)
        for sys in systems_inline
            pnl.extra_reset!(sys)
            pnl.extra_apply_freestream!(sys, uinf0)
            pnl._accumulate_Das!(sys, dt0 * eta_freestream)
        end
        for sys in systems_inline
            pnl.extra_reset!(sys)
        end
        pnl.kinematic_velocity!(systems_inline, frames_inline)
        for sys in systems_inline
            pnl._accumulate_Das!(sys, dt0 * eta_kinematic)
        end
        for sys in systems_inline
            pnl.reset!(sys)
        end

        @test length(body_helper.Das) == length(body_inline.Das)
        @test all(body_helper.Das[i] == body_inline.Das[i] for i in eachindex(body_helper.Das))

        # The default (arc) construction follows the trailing edge's swept path
        # instead of its tangent. Here θ = eta_kinematic*ω*dt0 = 0.02 rad, so the
        # two agree closely but are not identical.
        body_arc = make_plate_vortex_body()
        frames_arc = pnl.ReferenceFrame(body_arc;
            origin=SVector{3}(0.0, 0.0, 0.0),
            v=SVector{3}(0.0, 0.0, 0.0),
            ω_axis=SVector{3}(0.0, 0.0, 1.0),
            ω=2.0,
            R=SMatrix{3,3}(1.0, 0.0, 0.0,
                           0.0, 1.0, 0.0,
                           0.0, 0.0, 1.0),
            name="test",
            child_index=Int[],
            dependent_index=[1])
        pnl.initialize_Das!((body_arc,), frames_arc, Uinf, t0, dt0;
            set_Das_eta_freestream=eta_freestream,
            set_Das_eta_kinematic=eta_kinematic)
        @test all(isapprox(body_arc.Das[i], body_helper.Das[i]; rtol=1e-2)
                  for i in eachindex(body_arc.Das))
        @test !all(body_arc.Das[i] == body_helper.Das[i] for i in eachindex(body_arc.Das))
    end

    @testset "kinematic Das minimum displacement" begin
        eta = 0.1
        min_displacement = 0.2

        body_default = make_plate_vortex_body()
        body_default.velocity_te[1] .= [3.0 0.01;
                                        4.0 0.0;
                                        0.0 0.0]
        initial_default = copy(body_default.Das[1])
        expected_default = initial_default .+ copy(body_default.velocity_te[1]) .* eta
        pnl._accumulate_Das!(body_default, eta)
        @test body_default.Das[1] ≈ expected_default

        body_floored = make_plate_vortex_body()
        body_floored.velocity_te[1] .= body_default.velocity_te[1]
        initial_floored = copy(body_floored.Das[1])
        pnl._accumulate_Das!(body_floored, eta; min_displacement)

        @test body_floored.Das[1][:, 1] ≈ initial_floored[:, 1] .+ [0.3, 0.4, 0.0]
        @test body_floored.Das[1][:, 2] ≈ initial_floored[:, 2] .+ [min_displacement, 0.0, 0.0]

        body_zero = make_plate_vortex_body()
        initial_zero = copy(body_zero.Das[1])
        pnl._accumulate_Das!(body_zero, eta; min_displacement)
        @test body_zero.Das[1] ≈ initial_zero
    end

    @testset "seeded shedding trace" begin
        nodes, cells = make_seeded_te_mesh()
        final_cells = pnl.ensure_consistent_winding(nodes, cells; watertight=false)
        bbox = ([0.8, -0.1, -0.1], [1.1, 2.1, 0.1])

        trace = pnl.trace_trailing_edge(nodes, final_cells, 1, 2; bbox=bbox, end_node=3)
        @test trace.nodes == [1, 2, 3]
        @test length(trace.edges) == 2

        shedding = pnl.calc_shedding_from_seed(nodes, final_cells, 1, 2; bbox=bbox, end_node=3)
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
        @test pnl.calc_shedding_from_seed(nodes, final_cells, 1, 2; bbox=bbox_svec, end_node=3) == expected

        body = pnl.RigidWakeBody{pnl.VortexRing}(nodes, final_cells, shedding;
                                                 check_mesh=false, watertight=false,
                                                 ensure_winding=false)
        @test size(body.shedding[1]) == (6, 2)
        @test size(body.Das[1]) == (3, 3)
        @test body.cells == final_cells
        for i in axes(body.shedding[1], 2)
            pi, nia, nib, pj, nja, njb = body.shedding[1][:, i]
            @test Set(body.cells[[nia, nib], pi]) == Set([trace.nodes[i], trace.nodes[i+1]])
            @test Set(body.cells[[nja, njb], pj]) == Set([trace.nodes[i], trace.nodes[i+1]])
        end
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

end
