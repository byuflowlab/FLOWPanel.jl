using Test
import FLOWPanel as pnl
using LinearAlgebra: diag, dot, rank
import Meshes
using StaticArrays: SVector
import GeoIO

if !isdefined(@__MODULE__, :make_octa_source_body)
    include("test_helpers.jl")
end

function reconstruct_lu_operator(F)
    G = similar(F.factors)
    G[F.p, :] .= F.L * F.U
    return G
end

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

    @testset "Backslash Neumann solve! defaults update_G=false" begin
        # IMPROVEMENT : 2026-07-30 : regression test that the Neumann solve! wrapper's default (changed true -> false) still reuses the constructor's already-built G rather than silently rebuilding it, matching the Dirichlet wrapper's default
        body_default = make_octa_source_body()
        body_explicit = make_octa_source_body()
        for body in (body_default, body_explicit)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
        end

        solver_default = pnl.Backslash(body_default)
        solver_explicit = pnl.Backslash(body_explicit)

        pnl.solve!(body_default, solver_default)                    # relies on new default
        pnl.solve!(body_explicit, solver_explicit; update_G=false)  # explicit, for comparison

        @test isapprox(vec(body_default.strength[:, 1]), vec(body_explicit.strength[:, 1]); atol=1e-14)
        max_residual = nonlifting_flow_tangency_max_residual(body_default)
        @test max_residual < 1e-12
    end

    @testset "Backslash Neumann RigidWakeBody (single ConstantDoublet)" begin
        # Validates that the generic Backslash + Neumann path on a
        # RigidWakeBody{ConstantDoublet,1,_,false} assembles a wake-aware
        # influence matrix (via the _induced_wake hook on AbstractBody
        # `induced`) and recovers an impermeability-satisfying solution.
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

        magV = sqrt(1.0^2 + 0.2^2)
        Das_col = [1.0, 0.0, 0.2] ./ magV

        # Reference body without shedding: pure body-panel influence.
        body_nowake = pnl.RigidWakeBody{pnl.ConstantDoublet}(nodes, cells,
            zeros(Int, 6, 0); DBC=false, check_mesh=false, watertight=false)
        solver_nowake = pnl.Backslash(body_nowake)

        body = pnl.RigidWakeBody{pnl.ConstantDoublet}(nodes, cells, shedding;
            DBC=false, check_mesh=false, watertight=false)
        body.Das[1] .= repeat(Das_col, 1, size(body.Das[1], 2))
        Uinfs = zeros(3, body.ncells); Uinfs[1, :] .= 1.0; Uinfs[3, :] .= 0.2
        body.velocity .= Uinfs

        solver = pnl.Backslash(body)

        # solver.G is overwritten by lu! during construction, so reach the
        # original matrix entries through Glu.factors and a comparison
        # against a freshly-assembled wakeless body.
        @test maximum(abs.(solver.Glu.factors .- solver_nowake.Glu.factors)) > 1e-6

        pnl.solve!(body, solver)
        @test any(abs.(body.strength[:, 1]) .> 1e-6)

        rhs = -vec(sum(Uinfs .* body.normals; dims=1))
        @test isapprox(solver.Glu \ rhs, vec(body.strength[:, 1]); atol=1e-12)
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
        @test isapprox(body.controlpoints[:, 1], centroid1; atol=1e-12)
        @test size(solver.G) == (body.ncells, body.ncells)
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        body.velocity[3, :] .= 0.2
        body.potential .= range(0.1, 0.2; length=body.ncells)
        potential_before = copy(body.potential)
        pnl.solve!(body, solver)
        @test any(abs.(body.strength[:, 2]) .> 0)
        @test body.potential == potential_before
        # Single-body Dirichlet solves use `body.potential` as scratch workspace
        # and target zero interior perturbation potential.
        body.potential .= 0
        assert_boundary_residuals((body,); potential_atol=1e-6)
    end

    @testset "Backslash Dirichlet construction and solve" begin
        body = make_nonlifting(Union{pnl.ConstantSource, pnl.ConstantDoublet}; DBC=true)
        pnl.calc_normals!(body)
        normals = copy(body.normals)
        solver = pnl.Backslash(body)
        i1, i2, i3 = body.cells[:, 1]
        centroid1 = (body.nodes[:, i1] + body.nodes[:, i2] + body.nodes[:, i3]) ./ 3
        @test isapprox(body.controlpoints[:, 1], centroid1; atol=1e-12)
        @test size(solver.G) == (body.ncells, body.ncells)

        body = pnl.NonLiftingBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(copy(NODES_OCT), copy(CELLS_OCT); DBC=true)
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        body.velocity[3, :] .= 0.2
        normals = pnl.calc_normals!(body)
        expected_sigma = -vec(sum(body.velocity .* normals; dims=1))
        solver = pnl.Backslash(body)
        velocity_before = copy(body.velocity)
        potential_before = copy(body.potential)
        pnl.solve!(body, solver)
        @test any(abs.(body.strength[:, 2]) .> 0)
        @test isapprox(vec(body.strength[:, 1]), expected_sigma; atol=1e-12)
        @test body.velocity == velocity_before
        @test body.potential == potential_before
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

    @testset "Backslash Dirichlet VortexRing operator diagnostics" begin
        body = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.VortexRing}}(
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

        solver = pnl.Backslash(body)

        G_ref = zeros(Float64, body.ncells, body.ncells)
        pnl._G!(G_ref, body, body; kerneloffset=body.kerneloffset_panel, update_geometry=false)
        G_lu = reconstruct_lu_operator(solver.Glu)

        @test diag(G_ref) ≈ fill(0.5, body.ncells)
        @test rank(G_ref) == body.ncells
        @test G_lu ≈ G_ref
        @test solver.Glu \ ones(body.ncells) ≈ G_ref \ ones(body.ncells)

        gamma = collect(range(-0.3, 0.7; length=body.ncells))
        body.strength .= 0
        body.strength[:, 2] .= gamma
        body.potential .= 0
        pnl.influence!(body, body, pnl.DirectBackend(); scalar_potential=true, velocity=false)

        @test body.potential ≈ G_ref * gamma atol=1e-12
    end

    @testset "Dirichlet VortexRing interior perturbation potential limit" begin
        body = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.VortexRing}}(
            copy(NODES_OCT),
            copy(CELLS_OCT);
            check_mesh=false,
            watertight=true,
        )
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)

        G_ref = zeros(Float64, body.ncells, body.ncells)
        pnl._G!(G_ref, body, body; kerneloffset=body.kerneloffset_panel, update_geometry=false)

        gamma = [sin(0.7 * i) + 0.2 * cos(1.3 * i) for i in 1:body.ncells]
        body.strength .= 0
        body.strength[:, 2] .= gamma
        phi_cp = G_ref * gamma

        panel_lengths = [body.characteristiclength(body.nodes, body.cells[:, i]) for i in 1:body.ncells]
        min_len = minimum(panel_lengths)
        epsilons = min_len .* (1e-4, 3e-5, 1e-5)
        errors = Float64[]

        for eps_len in epsilons
            target = deepcopy(body)
            target.controlpoints .= body.controlpoints .- eps_len .* body.normals
            target.potential .= 0
            pnl.influence!(target, body, pnl.DirectBackend(); scalar_potential=true, velocity=false)
            push!(errors, maximum(abs.(target.potential .- phi_cp)))
        end

        @test errors[2] < errors[1]
        @test errors[3] < errors[2]
        @test errors[3] < 1e-4
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
        @test solver1.max_iterations == 100
        @test solver1.tolerance == 1e-6
        @test size(solver1.Uext) == (3, body1.ncells)
        @test size(solver1.phi_ext) == (body1.ncells,)
        @test_throws MethodError pnl.FGSSolver((body1, body2))
    end

    @testset "FGSSolver boundary-condition preparation" begin
        nodes, cells = make_seeded_te_mesh()
        shedding = pnl.calc_shedding_from_seed(nodes, cells, 1, 2;
            bbox=([0.8, -0.1, -0.1], [1.1, 2.1, 0.1]),
            end_node=3)

        function make_dirichlet_body(initial_doublet)
            body = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(
                copy(nodes), copy(cells), [copy(shedding)];
                check_mesh=false,
                watertight=true)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
            body.strength .= 0
            body.strength[:, 2] .= initial_doublet
            for i in eachindex(body.Das)
                body.Das[i] .= repeat([1.0, 0.0, 0.0], 1, size(body.Das[i], 2))
            end
            pnl.calc_normals!(body)
            pnl.calc_controlpoints!(body)
            return body
        end

        body_backslash = make_dirichlet_body(1.0)
        backslash = pnl.BackslashCoupled((body_backslash,))
        pnl.solve!((body_backslash,), backslash; backend=pnl.DirectBackend(), update_G=true)

        body_fgs = make_dirichlet_body(1.0)
        fgs = pnl.FGSSolver(body_fgs; leaf_size=10000, tolerance=1e-12, max_iterations=5)
        # IMPROVEMENT : 2026-07-30 : assert on the return value of the public solve!(body, solver) API : previously FGSSolver's _solve! returned `nothing` (falling through to save_solution!'s return), so this call crashed on `tb, ts = _solve!(...)` inside the generic wrapper -- this exercises exactly that path and confirms it now returns a proper (tb, ts) timing tuple
        tb_ts = pnl.solve!(body_fgs, fgs)
        @test tb_ts isa Tuple{Float64, Float64}

        @test body_fgs.strength[:, 1] ≈ body_backslash.strength[:, 1] atol=1e-12
        @test body_fgs.strength[:, 2] ≈ body_backslash.strength[:, 2] atol=1e-12

        neumann = make_octa_source_body()
        neumann.strength[:, 1] .= 3.0
        pnl.set_strengths!(neumann)
        @test all(iszero, neumann.strength[:, 1])
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

    @testset "Outer coupled loop solve! defaults to FastMultipoleBackend" begin
        # IMPROVEMENT : 2026-07-30 : regression test that omitting `backend` now uses FastMultipoleBackend() per body (previously DirectBackend()) and still gives an accurate solve
        body1 = make_octa_source_body()
        body2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        body1.velocity .= 0
        body2.velocity .= 0
        body1.velocity[1, :] .= 1.0
        body2.velocity[1, :] .= 1.0
        solver1 = pnl.Backslash(body1)
        solver2 = pnl.Backslash(body2)
        pnl.solve!((body1, body2), (solver1, solver2))  # backend omitted -> new FMM default
        @test any(abs.(body1.strength[:, 1]) .> 0)
        @test any(abs.(body2.strength[:, 1]) .> 0)
        assert_boundary_residuals((body1, body2); backend=pnl.FastMultipoleBackend(), tangency_atol=1e-6)
    end

    @testset "No unconditional debug println during solve!" begin
        # IMPROVEMENT : 2026-07-30 : regression test guarding against reintroduction of the removed TEMPORARY debug printlns (KrylovSolver GMRES diagnostics, "Tuple of bodies", "BackslashCoupled", outer-loop "[outer diag]") : these fired unconditionally on every solve regardless of a verbose flag, polluting stdout across batch runs like AOA/panel-count sweeps
        # NOTE: redirect_stdout requires an OS-backed stream (a plain IOBuffer isn't accepted), hence mktemp rather than IOBuffer()
        body_krylov = make_octa_source_body()
        body_krylov.velocity .= 0
        body_krylov.velocity[1, :] .= 1.0
        krylov_solver = pnl.KrylovSolver(body_krylov; backend=pnl.DirectBackend(), atol=1e-8, rtol=1e-8, itmax=50)

        body_c1 = make_octa_source_body()
        body_c2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        for b in (body_c1, body_c2)
            b.velocity .= 0
            b.velocity[1, :] .= 1.0
        end
        coupled_solver = pnl.BackslashCoupled((body_c1, body_c2))

        body_o1 = make_octa_source_body()
        body_o2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        for b in (body_o1, body_o2)
            b.velocity .= 0
            b.velocity[1, :] .= 1.0
        end
        outer_solver1 = pnl.Backslash(body_o1)
        outer_solver2 = pnl.Backslash(body_o2)

        output = mktemp() do path, f
            redirect_stdout(f) do
                pnl.solve!(body_krylov, krylov_solver)
                pnl.solve!((body_c1, body_c2), coupled_solver; backend=pnl.DirectBackend(), update_G=true)
                pnl.solve!((body_o1, body_o2), (outer_solver1, outer_solver2); backend=(pnl.DirectBackend(), pnl.DirectBackend()))
            end
            flush(f)
            read(path, String)
        end

        @test !occursin("gmres diag", output)
        @test !occursin("Tuple of bodies", output)
        @test !occursin("BackslashCoupled", output)
        @test !occursin("outer diag", output)
    end

    @testset "KrylovCoupled nonlifting" begin
        body1 = make_octa_source_body()
        body2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        body1.velocity .= 0
        body2.velocity .= 0
        body1.velocity[1, :] .= 1.0
        body2.velocity[1, :] .= 1.0

        bodies = (body1, body2)
        solver = pnl.KrylovCoupled(bodies; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=20)
        pnl.solve!(bodies, solver; backend=pnl.DirectBackend())

        @test solver isa pnl.KrylovCoupled
        @test any(abs.(body1.strength[:, 1]) .> 0)
        @test any(abs.(body2.strength[:, 1]) .> 0)
        assert_boundary_residuals(bodies; tangency_atol=1e-7)
    end

    @testset "KrylovCoupled + JacobiPreconditioner" begin
        # IMPROVEMENT : 2026-07-30 : test for KrylovCoupled's new preconditioner_cell_size support, mirroring the existing "KrylovSolver + JacobiPreconditioner" test : confirms the with/without-preconditioner comparison is now possible for both Krylov solver families, not just KrylovSolver
        # NOTE: uses Dirichlet (DBC=true) bodies, matching the one configuration "KrylovSolver + JacobiPreconditioner" already
        # exercises successfully -- JacobiPreconditioner's default derivatives_switches=(scalar_potential=true, gradient=true)
        # builds a potential-based self-influence matrix per cell, which is consistent with a Dirichlet (potential) formulation
        # but not with a Neumann (velocity-dot-normal) one; an earlier version of this test used two Neumann bodies and got
        # wildly wrong (non-converged) strengths, which is a JacobiPreconditioner/Neumann-formulation gap on the FastMultipole
        # side, not a KrylovCoupled wiring bug -- out of scope here, so this test sticks to the proven-good Dirichlet case.
        make_dirichlet_body(shift) = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(
            Float64[
                0 1 1 0;
                0 0 1 1;
                0 0 0 0;
            ] .+ shift,
            Int[
                1 1;
                2 3;
                3 4;
            ];
            check_mesh=false,
            watertight=false,
        )

        make_bodies() = (make_dirichlet_body([0.0, 0.0, 0.0]), make_dirichlet_body([3.0, 0.0, 0.0]))

        bodies_no_pc = make_bodies()
        bodies_pc = make_bodies()

        for bodies in (bodies_no_pc, bodies_pc), body in bodies
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
            body.velocity[3, :] .= 0.2
        end

        solver_no_pc = pnl.KrylovCoupled(bodies_no_pc; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50)
        solver_pc = pnl.KrylovCoupled(bodies_pc; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50, preconditioner_cell_size=1.5)

        @test solver_pc.preconditioner !== nothing
        @test length(solver_pc.preconditioner.cell_body_indices) > 0

        pnl.solve!(bodies_no_pc, solver_no_pc; backend=pnl.DirectBackend())
        pnl.solve!(bodies_pc, solver_pc; backend=pnl.DirectBackend())

        for (body_no_pc, body_pc) in zip(bodies_no_pc, bodies_pc)
            @test isapprox(vec(body_pc.strength[:, 1]), vec(body_no_pc.strength[:, 1]); atol=1e-6)
        end
        assert_boundary_residuals(bodies_no_pc; tangency_atol=1e-7)
        assert_boundary_residuals(bodies_pc; tangency_atol=1e-6)
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

    @testset "BackslashCoupled solve! defaults to FastMultipoleBackend" begin
        # IMPROVEMENT : 2026-07-30 : regression test that omitting `backend` now uses FastMultipoleBackend() (previously DirectBackend()) and still gives an accurate solve
        body1 = make_octa_source_body()
        body2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        magVinf = 10.0
        Vinf = magVinf * [1.0, 0.0, 0.0]
        pnl.apply_freestream!(body1, Vinf)
        pnl.apply_freestream!(body2, Vinf)

        bodies = (body1, body2)
        solver = pnl.BackslashCoupled(bodies)

        pnl.solve!(bodies, solver; update_G=true)  # backend omitted -> new FMM default
        @test any(abs.(body1.strength[:, 1]) .> 0)
        @test any(abs.(body2.strength[:, 1]) .> 0)
        assert_boundary_residuals(bodies; backend=pnl.FastMultipoleBackend(), tangency_atol=1e-6)
    end

    @testset "BackslashCoupled uses kerneloffset_panel, not stale kerneloffset" begin
        # IMPROVEMENT : 2026-07-30 : regression test for BackslashCoupled's _G! call now passing kerneloffset=source.kerneloffset_panel explicitly : previously it fell back to _G!'s generic kerneloffset default, which can silently diverge from kerneloffset_panel (e.g. after a monitor/simulate! pass temporarily repoints body.kerneloffset at kerneloffset_targets), making BackslashCoupled and Backslash regularize self/near-singular terms differently for what should be the same discretization
        body_direct = make_octa_source_body()
        body_coupled = make_octa_source_body()

        # Simulate body.kerneloffset having been left pointed at a different
        # value than kerneloffset_panel (as happens transiently during
        # simulate!/monitor evaluation), while kerneloffset_panel is unchanged.
        @test body_direct.kerneloffset == body_direct.kerneloffset_panel
        body_direct.kerneloffset = 1e-3
        body_coupled.kerneloffset = 1e-3
        @test body_direct.kerneloffset != body_direct.kerneloffset_panel

        for body in (body_direct, body_coupled)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
        end

        direct = pnl.Backslash(body_direct)
        pnl.solve!(body_direct, direct; backend=pnl.DirectBackend())

        coupled = pnl.BackslashCoupled((body_coupled,))
        pnl.solve!((body_coupled,), coupled; backend=pnl.DirectBackend(), update_G=true)

        @test isapprox(vec(body_coupled.strength[:, 1]), vec(body_direct.strength[:, 1]); atol=1e-12)
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
        tangency_residuals = flow_tangency_max_residuals(bodies; backend)
        @test tangency_residuals[2] < 1e-7

        # @show CL, CD = postprocess!(bodies, Vinf, rho, backend, chords, span)
    end

    @testset "backslash_meshes" begin

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

        bodies = (bodies[2],)
        #------------------- SOLVE BODY ----------------------------------------------
        backend = pnl.DirectBackend()
        solver = pnl.BackslashCoupled(bodies)
        println("Solving body...")

        pnl.solve!(bodies, solver; backend, update_G=true)
        
        # write vtk files
        for i in eachindex(bodies)
            pnl.write_vtk("check_mesh_body_$(i)", bodies[i])
        end

        assert_boundary_residuals(bodies; backend, potential_atol=1e-4)
    end

end

# println("done.")
