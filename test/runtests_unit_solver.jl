using Test
import FLOWPanel as pnl
using LinearAlgebra: diag, dot, rank, ldiv!
import Meshes
import SparseArrays
using StaticArrays: SVector

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

    @testset "Backslash Glu write-back on update_G (regression)" begin
        # Regression for the stale-Glu aliasing bug: lu!(G) aliases solver.G as
        # Glu.factors, so a local re-factorization under update_G=true used to
        # leave solver.Glu with fresh factors but stale pivots — corrupting any
        # later direct consumer (Kutta Route A, Green-family formulations).
        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        solver = pnl.Backslash(body)
        pnl.solve!(body, solver)
        strengths0 = copy(body.strength[:, 1])

        # geometry unchanged: an update_G re-solve must reproduce the solution
        # and leave solver.Glu a consistent (factors, ipiv) pair
        pnl._solve!(body, solver; update_G=true)
        @test isapprox(body.strength[:, 1], strengths0; atol=1e-12)

        G_ref = zeros(body.ncells, body.ncells)
        pnl._G!(G_ref, body, body; kerneloffset=body.kerneloffset_panel)
        @test reconstruct_lu_operator(solver.Glu) ≈ G_ref atol=1e-10

        rhs = [sin(1.1 * i) for i in 1:body.ncells]
        @test isapprox(solver.Glu \ rhs, G_ref \ rhs; atol=1e-10)
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

    @testset "KrylovSolver fgmres matches gmres" begin
        body_g = make_octa_source_body()
        body_f = make_octa_source_body()
        for body in (body_g, body_f)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
        end

        solver_g = pnl.KrylovSolver(body_g; method=:gmres, backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50)
        solver_f = pnl.KrylovSolver(body_f; method=:fgmres, backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50)
        pnl.solve!(body_g, solver_g)
        pnl.solve!(body_f, solver_f)

        @test isapprox(vec(body_f.strength[:, 1]), vec(body_g.strength[:, 1]); atol=1e-8)
        assert_boundary_residuals((body_f,); tangency_atol=1e-7)
    end

    @testset "KrylovSolver persistent workspace" begin
        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        solver = pnl.KrylovSolver(body; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50)

        ws = solver.workspace
        pnl.solve!(body, solver)
        strengths1 = copy(body.strength[:, 1])
        pnl.solve!(body, solver)

        @test solver.workspace === ws            # no per-solve reallocation
        @test isapprox(body.strength[:, 1], strengths1; atol=1e-12)
        @test solver.niter > 0
        @test solver.solved
    end

    @testset "FGSPreconditioner linearity and side effects" begin
        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0

        P = pnl.FGSPreconditioner(body; sweeps=1, inner_iterations=2, leaf_size=10000)
        P_uncached = pnl.FGSPreconditioner(body; sweeps=1, inner_iterations=2,
            leaf_size=10000, cache_leaf_lu=false)

        @test P.solver.cache_leaf_lu
        @test P.solver.fgs.cache_leaf_lu
        @test P.solver.fgs.leaf_lu_cache !== nothing
        @test !P_uncached.solver.cache_leaf_lu
        @test !P_uncached.solver.fgs.cache_leaf_lu
        @test P_uncached.solver.fgs.leaf_lu_cache === nothing
        @test pnl._preconditioner_metadata_dict(P)["cache_leaf_lu"] == true
        @test pnl._preconditioner_metadata_dict(P_uncached)["cache_leaf_lu"] == false

        strength0 = copy(body.strength)
        velocity0 = copy(body.velocity)
        potential0 = copy(body.potential)

        n = body.ncells
        x1 = [sin(0.9 * i) for i in 1:n]
        x2 = [cos(1.7 * i) - 0.3 for i in 1:n]
        α, β = 0.7, -1.9

        y1, y2, y12 = zeros(n), zeros(n), zeros(n)
        y1_uncached = zeros(n)
        ldiv!(y1, P, x1)
        ldiv!(y1_uncached, P_uncached, x1)
        ldiv!(y2, P, x2)
        ldiv!(y12, P, α .* x1 .+ β .* x2)

        @test isapprox(y12, α .* y1 .+ β .* y2; atol=1e-10)
        @test isapprox(y1, y1_uncached; atol=1e-12)

        # side-effect-free apply
        @test body.strength == strength0
        @test body.velocity == velocity0
        @test body.potential == potential0

        # single-leaf FGS apply is an exact G⁻¹: check against the dense operator
        G = zeros(n, n)
        pnl._G!(G, body, body; kerneloffset=body.kerneloffset_panel)
        @test isapprox(G * y1, x1; atol=1e-8)
    end

    @testset "KrylovSolver fgmres + FGSPreconditioner (config 1f)" begin
        body_ref = make_octa_source_body()
        body_pc = make_octa_source_body()
        body_no_pc = make_octa_source_body()
        for body in (body_ref, body_pc, body_no_pc)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
            body.velocity[3, :] .= 0.2
        end

        ref = pnl.Backslash(body_ref)
        pnl.solve!(body_ref, ref)

        solver_no_pc = pnl.KrylovSolver(body_no_pc; method=:gmres, backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=100)
        pnl.solve!(body_no_pc, solver_no_pc)

        P = pnl.FGSPreconditioner(body_pc; sweeps=1, inner_iterations=2, leaf_size=10000)
        solver_pc = pnl.KrylovSolver(body_pc; method=:fgmres, preconditioner=P, backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=100)
        pnl.solve!(body_pc, solver_pc)

        @test isapprox(vec(body_pc.strength[:, 1]), vec(body_ref.strength[:, 1]); atol=1e-6)
        @test solver_pc.solved
        @test solver_pc.niter <= solver_no_pc.niter
        assert_boundary_residuals((body_pc,); tangency_atol=1e-6)
    end

    @testset "FMM-direct-list ILU(0) construction and apply" begin
        for body in (make_octa_source_body(), make_dirichlet_diamond_body(),
                     make_plate_vortex_body())
            n = body.ncells
            body.strength .= reshape(sin.(1:length(body.strength)), size(body.strength))
            body.velocity .= reshape(cos.(1:length(body.velocity)), size(body.velocity))
            body.potential .= range(-0.2, 0.3; length=n)
            state0 = (strength=copy(body.strength), velocity=copy(body.velocity),
                      potential=copy(body.potential), normals=copy(body.normals),
                      controlpoints=copy(body.controlpoints))

            # One leaf makes the directed Barba pattern dense. ILU(0) then
            # drops no fill and must agree with a direct solve.
            P = pnl.ILUPreconditioner(body; leaf_size=10_000, keep_matrix=true)
            @test P.interaction_list_method == :Barba
            @test P.stats["nnz"] == n^2
            @test P.stats["factor_nnz"] == n^2
            @test P.stats["max_row_nnz"] == n
            @test all(i -> P.matrix[i, i] != 0, 1:n)
            @test sort(P.permutation) == collect(1:n)
            @test P.inverse_permutation[P.permutation] == collect(1:n)

            G = zeros(n, n)
            pnl._G!(G, body, body; kerneloffset=body.kerneloffset_panel)
            @test P.matrix ≈ G[P.permutation, P.permutation] atol=1e-12

            x1 = [sin(0.9i) for i in 1:n]
            x2 = [cos(1.7i) - 0.3 for i in 1:n]
            y1, y2, y12 = zeros(n), zeros(n), zeros(n)
            ldiv!(y1, P, x1)
            ldiv!(y2, P, x2)
            ldiv!(y12, P, 0.7 .* x1 .- 1.9 .* x2)
            @test y12 ≈ 0.7 .* y1 .- 1.9 .* y2 atol=1e-10
            @test y1 ≈ G \ x1 atol=1e-10

            # Constructor and apply are side-effect-free.
            @test body.strength == state0.strength
            @test body.velocity == state0.velocity
            @test body.potential == state0.potential
            @test body.normals == state0.normals
            @test body.controlpoints == state0.controlpoints
        end
    end

    @testset "ILU(0) sparse pattern, guard, scaling, and Krylov routing" begin
        # Nontrivial tree ordering: stored entries are exactly the corresponding
        # dense _G! entries, with no off-pattern storage.
        body = make_sphere_source_body(ntheta=4, nphi=8)
        P = pnl.ILUPreconditioner(body; leaf_size=4,
            multipole_acceptance=1.0, keep_matrix=true)
        target_tree, source_tree, direct_list, _, requested =
            pnl._ilu_direct_pattern(body, 4, 1.0, 512 * body.ncells)
        expanded = Set{Tuple{Int,Int}}()
        for pair in direct_list
            tr = target_tree.branches[pair[1]].bodies_index[1]
            sr = source_tree.branches[pair[2]].bodies_index[1]
            for it in tr, js in sr
                push!(expanded, (target_tree.sort_index_list[1][it],
                                 source_tree.sort_index_list[1][js]))
            end
        end
        for i in 1:body.ncells
            push!(expanded, (i, i))
        end
        @test length(expanded) == requested
        I, J, V = SparseArrays.findnz(P.matrix)
        @test length(Set(zip(I, J))) == length(V)
        stored_original = Set((P.permutation[I[k]], P.permutation[J[k]])
                              for k in eachindex(V))
        @test stored_original == expanded
        @test all(i -> P.matrix[i, i] != 0, 1:body.ncells)
        @test P.stats["nnz"] < body.ncells^2

        # The linear guard fires before triplet allocation/kernel evaluation
        # and restores constructor-visible body state.
        guard_body = make_sphere_source_body(ntheta=4, nphi=8)
        state0 = copy(guard_body.strength)
        err = try
            pnl.ILUPreconditioner(guard_body; leaf_size=10_000,
                max_pattern_entries=guard_body.ncells)
            nothing
        catch ex
            ex
        end
        @test err isa ArgumentError
        @test occursin("before sparse allocation and kernel evaluation", sprint(showerror, err))
        @test guard_body.strength == state0

        # Equilibrated DS applies D to the right-hand side and still represents
        # an approximation to S^{-1}; with a dense pattern it is exact.
        eq_body = make_octa_source_body()
        Peq = pnl.ILUPreconditioner(eq_body; leaf_size=10_000,
            equilibrate=true, diagonal_shift=1e-8, keep_matrix=true)
        rhs = [sin(i) for i in 1:eq_body.ncells]
        sol = similar(rhs)
        ldiv!(sol, Peq, rhs)
        shifted = copy(Peq.matrix)
        @test shifted * sol[Peq.permutation] ≈ rhs[Peq.permutation] atol=1e-9

        # Metadata retains every construction knob and accounting field.
        md = pnl._preconditioner_metadata_dict(P)
        for key in ("leaf_size", "multipole_acceptance", "interaction_list_method",
                    "equilibrate", "diagonal_shift", "max_pattern_entries",
                    "nnz", "nnz_per_panel", "max_row_nnz", "factor_nnz")
            @test haskey(md, key)
        end

        body_ref = make_octa_source_body()
        body_pc = make_octa_source_body()
        for b in (body_ref, body_pc)
            b.velocity .= 0
            b.velocity[1, :] .= 1.0
            b.velocity[3, :] .= 0.2
        end
        ref = pnl.Backslash(body_ref)
        pnl.solve!(body_ref, ref)
        solver = pnl.KrylovSolver(body_pc; method=:gmres,
            preconditioner=pnl.ILUPreconditioner(body_pc; leaf_size=10_000),
            backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50)
        pnl.solve!(body_pc, solver)
        @test solver.solved
        @test body_pc.strength[:, 1] ≈ body_ref.strength[:, 1] atol=1e-8
    end

    @testset "ILU(0) construction ladder remains linear" begin
        samples = Tuple{Int,Int}[]
        for (ntheta, nphi) in ((8, 16), (12, 22), (16, 32))
            body = make_sphere_source_body(; ntheta, nphi)
            P = pnl.ILUPreconditioner(body; leaf_size=10,
                multipole_acceptance=1.0)
            push!(samples, (P.stats["nnz"], P.stats["retained_bytes"]))
        end
        for i in 2:length(samples)
            @test samples[i][1] <= 2.5 * samples[i-1][1]
            @test samples[i][2] <= 2.5 * samples[i-1][2]
        end
    end

    @testset "ILU(0) direct list agrees with fmm!" begin
        # Keep this last: fmm! exercises the complete evaluator whereas ILU
        # construction intentionally uses only its tree/list portion.
        body = make_octa_source_body()
        target_tree, source_tree, direct_list, _, _ =
            pnl._ilu_direct_pattern(body, 4, 1.0, 512 * body.ncells)
        _, _, fmm_target_tree, fmm_source_tree, _, fmm_direct_list, _, _ =
            pnl.FastMultipole.fmm!(deepcopy(body);
                scalar_potential=false, gradient=true, hessian=false,
                leaf_size=4, multipole_acceptance=1.0,
                interaction_list_method=pnl.FastMultipole.Barba(),
                nearfield=true, farfield=true, self_induced=true)
        @test Set(Tuple.(direct_list)) == Set(Tuple.(fmm_direct_list))
        @test target_tree.sort_index_list == fmm_target_tree.sort_index_list
        @test source_tree.sort_index_list == fmm_source_tree.sort_index_list
    end

    @testset "KrylovSolver history capture" begin
        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        solver = pnl.KrylovSolver(body; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50, record_history=true)
        pnl.solve!(body, solver)

        h = solver.history
        @test length(h) > 0
        @test length(h) >= solver.niter
        @test issorted(h.t_ns)
        @test h.metric == :krylov_precnorm
        @test all(h.t_ns .>= h.t0_ns)

        # history resets on the next solve rather than accumulating
        len1 = length(h)
        pnl.solve!(body, solver)
        @test length(solver.history) <= len1 + 1
    end

    @testset "KrylovSolver warmstart (Neumann)" begin
        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        solver = pnl.KrylovSolver(body; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50)

        @test !solver.warmstart          # off by default
        pnl.solve!(body, solver)
        niter_cold = solver.niter
        strengths_cold = copy(body.strength[:, 1])
        @test solver.have_x_prev

        solver.warmstart = true
        pnl.solve!(body, solver)         # identical RHS, seeded with the exact solution
        @test solver.niter <= niter_cold
        @test isapprox(body.strength[:, 1], strengths_cold; atol=1e-8)

        # changed RHS: warmstarted solve still converges to the right answer
        body.velocity .= 0
        body.velocity[2, :] .= 1.0
        body_ref = make_octa_source_body()
        body_ref.velocity .= 0
        body_ref.velocity[2, :] .= 1.0
        ref = pnl.Backslash(body_ref)
        pnl.solve!(body_ref, ref)

        pnl.solve!(body, solver)
        @test isapprox(body.strength[:, 1], body_ref.strength[:, 1]; atol=1e-7)
    end

    @testset "KrylovSolver warmstart (Dirichlet)" begin
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
        body = make_dirichlet_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        body.velocity[3, :] .= 0.2

        solver = pnl.KrylovSolver(body; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50)
        pnl.solve!(body, solver)
        niter_cold = solver.niter
        mu_cold = copy(body.strength[:, 2])

        solver.warmstart = true
        pnl.solve!(body, solver)
        @test solver.niter <= niter_cold
        @test isapprox(body.strength[:, 2], mu_cold; atol=1e-8)
        assert_boundary_residuals((body,); potential_atol=1e-7)
    end

    @testset "KrylovSolver cache_tree (021 Phase 2b)" begin
        # FMM-plan reuse across a solve's operator applies must be bitwise
        # equivalent to the per-apply tree-rebuild path, and per-solve scoping
        # must survive kerneloffset changes between solves.
        backend = pnl.FastMultipoleBackend(; expansion_order=8,
            multipole_acceptance=0.4, leaf_size=20)
        kw = (; backend, method=:gmres, atol=1e-12, rtol=1e-10, itmax=200)

        function fresh_sphere()
            body = make_sphere_source_body()
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
            return body
        end

        # --- bitwise equivalence, cached vs uncached -----------------------
        body_ref = fresh_sphere()
        solver_ref = pnl.KrylovSolver(body_ref; kw...)
        pnl.solve!(body_ref, solver_ref)

        body_c = fresh_sphere()
        solver_c = pnl.KrylovSolver(body_c; kw..., cache_tree=true)
        pnl.solve!(body_c, solver_c)

        @test solver_c.niter == solver_ref.niter
        @test solver_c.niter > 1              # premise: the plan was reused
        @test body_c.strength == body_ref.strength   # bitwise
        @test solver_c.kop.plan_slot[] === nothing   # dropped at solve end

        # --- allocation reduction ------------------------------------------
        alloc_ref = @allocated pnl.solve!(body_ref, solver_ref)
        alloc_c = @allocated pnl.solve!(body_c, solver_c)
        @test alloc_ref > 10_000_000          # premise: rebuild path allocates
        @test alloc_c < 0.6 * alloc_ref

        # --- per-solve invalidation under geometry change ------------------
        # the plan is rebuilt every solve by construction; verify that a
        # geometry change between solves is honored bitwise (a stale-plan bug
        # would reproduce the OLD geometry's solution instead). The offset
        # knob is inert on this fixture (measured: kerneloffset changes leave
        # the solve bit-identical), so geometry is the discriminating input.
        function perturb!(body)
            body.nodes[:, 1] .+= 0.05
            pnl.calc_normals!(body)
            pnl.calc_controlpoints!(body)
            return body
        end
        perturb!(body_c)
        pnl.solve!(body_c, solver_c)
        body_f = perturb!(fresh_sphere())
        solver_f = pnl.KrylovSolver(body_f; kw...)
        pnl.solve!(body_f, solver_f)
        @test body_c.strength == body_f.strength     # bitwise
        # premise: the perturbation actually changed the solution
        @test body_c.strength != body_ref.strength

        # --- guards + metadata ---------------------------------------------
        @test_throws ArgumentError pnl.KrylovSolver(fresh_sphere();
            backend=pnl.DirectBackend(), cache_tree=true)
        md = pnl._solver_metadata_dict(solver_c)
        @test md["cache_tree"] === true
        @test pnl._solver_metadata_dict(solver_ref)["cache_tree"] === false
    end

    @testset "KrylovSolver cache_nearfield (021 Phase 2b)" begin
        # Dense near-field cache: the cached operator apply must reproduce the
        # kernel apply to rtol 1e-12 (NOT bitwise — BLAS sums in a different
        # order), on a body WITH shedding panels so the attached-wake term
        # `_induced_wake` is exercised (this is the plan's V0 linearity check:
        # the cache is built by single-panel unit-strength probes, so any
        # strength-like input outside buffer rows 5:4+strength_dims, or any
        # nonlinearity, breaks the comparison).
        backend = pnl.FastMultipoleBackend(; expansion_order=8,
            multipole_acceptance=0.4, leaf_size=16)

        # --- V0: Dirichlet apply_G exactness incl. attached wake -----------
        body = make_dirichlet_diamond_body(nspan=40)
        @test body.nsheddings > 0             # premise: shedding panels present
        @test any(!iszero, vcat(vec.(body.Das)...))  # premise: wake direction set
        n = body.ncells
        x = sin.(0.7 .* (1:n)) .+ 0.1
        scratch = zeros(n)

        y_ref = copy(pnl._apply_dirichlet_G!(body, x, backend, scratch))
        slot = Ref{Any}(nothing)
        y_cached = copy(pnl._apply_dirichlet_G!(body, x, backend, scratch;
            plan_slot=slot, cache_nearfield=true))

        # premise guards: a plan with a real near field was built and carries
        # the cache (bitwise inequality proves the BLAS path actually ran)
        plan = slot[][1]
        @test length(plan.direct_list) > 0
        @test plan.nearfield_cache[] isa FastMultipole.NearfieldInfluenceCache
        @test any(!iszero, y_ref)
        @test y_cached != y_ref
        @test isapprox(y_cached, y_ref; rtol=1e-12)

        # strength-change reuse through the SAME plan/cache
        x2 = cos.(0.3 .* (1:n)) .- 0.05
        y_cached2 = copy(pnl._apply_dirichlet_G!(body, x2, backend, scratch;
            plan_slot=slot, cache_nearfield=true))
        @test slot[][1] === plan              # premise: plan (and cache) reused
        y_ref2 = copy(pnl._apply_dirichlet_G!(body, x2, backend, scratch))
        @test !isapprox(y_ref2, y_ref; rtol=1e-6)  # premise: answer changed
        @test isapprox(y_cached2, y_ref2; rtol=1e-12)

        # --- Neumann solve equality ----------------------------------------
        kw = (; backend, method=:gmres, atol=1e-12, rtol=1e-10, itmax=200)
        function fresh_sphere()
            body = make_sphere_source_body()
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
            return body
        end

        body_ref = fresh_sphere()
        solver_ref = pnl.KrylovSolver(body_ref; kw...)
        pnl.solve!(body_ref, solver_ref)

        body_c = fresh_sphere()
        solver_c = pnl.KrylovSolver(body_c; kw..., cache_nearfield=true)
        @test solver_c.cache_tree === true    # implied by cache_nearfield
        pnl.solve!(body_c, solver_c)

        @test solver_c.solved
        @test solver_c.niter > 1              # premise: the cache was reused
        @test body_c.strength != body_ref.strength   # premise: BLAS path ran
        @test isapprox(body_c.strength, body_ref.strength; rtol=1e-10)
        @test solver_c.kop.plan_slot[] === nothing   # dropped at solve end
        # tangency at FMM truncation level (p=8 ≈ 1e-8 field error; the solve
        # itself converged to rtol 1e-10 — measured uncached residual is the
        # same 3e-8, so this is backend truncation, not a cache artifact)
        assert_boundary_residuals((body_c,); tangency_atol=1e-6)

        # --- guards + metadata ---------------------------------------------
        @test_throws ArgumentError pnl.KrylovSolver(fresh_sphere();
            backend=pnl.DirectBackend(), cache_nearfield=true)
        md = pnl._solver_metadata_dict(solver_c)
        @test md["cache_nearfield"] === true
        @test md["cache_tree"] === true
        @test pnl._solver_metadata_dict(solver_ref)["cache_nearfield"] === false
    end

    @testset "KrylovSolver persistent_plan + rigid motion (021 tree reuse)" begin
        # Cross-solve plan (+ near-field cache) persistence under rigid
        # motion: the Dirichlet scalar operator is exactly invariant, so a
        # plan/cache built before the motion, transformed via
        # transform_plan!/transform_solver_geometry!, must keep reproducing
        # the operator to rtol 1e-12 — shedding panels included so the
        # attached-wake term is exercised.
        backend = pnl.FastMultipoleBackend(; expansion_order=8,
            multipole_acceptance=0.4, leaf_size=16)

        # rigid motion: rotation by 30 deg about an off-axis direction, about
        # an off-body origin, plus a translation
        axis = [0.3, -1.0, 0.7] ./ norm([0.3, -1.0, 0.7])
        Rrot = pnl.Rodrigues(axis, 30.0 * pi / 180)
        origin = [0.1, 0.2, -0.05]
        dx = [0.02, -0.03, 0.05]
        t_affine = origin + dx - Rrot * origin
        function move!(body)
            pnl.rotate_translate!(body, origin, Rrot, dx)
            pnl.rotate_Das!(body, Rrot)
            pnl.calc_normals!(body)
            pnl.calc_controlpoints!(body)
            return body
        end

        # --- Dirichlet apply exactness through a transformed cache ---------
        body = make_dirichlet_diamond_body(nspan=40)
        @test body.nsheddings > 0             # premise: shedding panels present
        n = body.ncells
        x = sin.(0.7 .* (1:n)) .+ 0.1
        scratch = zeros(n)

        slot = Ref{Any}(nothing)
        y0 = copy(pnl._apply_dirichlet_G!(body, x, backend, scratch;
            plan_slot=slot, cache_nearfield=true))
        plan = slot[][1]
        @test length(plan.direct_list) > 0    # premise: near field exercised
        # (no m2l premise: kerneloffset radius inflation makes this fixture
        # all-direct at any MAC — measured m2l=0 up to mac=2.0. The panel
        # far-field-under-transform path is covered by FastMultipole's
        # source-panel transform_tree! test.)
        @test plan.nearfield_cache[] isa FastMultipole.NearfieldInfluenceCache
        @test any(!iszero, y0)

        move!(body)
        FastMultipole.transform_plan!(plan, (body,), Rrot, t_affine)

        # reused cache on the moved geometry
        y1 = copy(pnl._apply_dirichlet_G!(body, x, backend, scratch;
            plan_slot=slot, cache_nearfield=true))
        @test slot[][1] === plan              # premise: plan+cache reused

        # scalar G is exactly invariant under rigid motion (same tree, same
        # arithmetic path — differences only from rotated coordinates)
        @test isapprox(y1, y0; rtol=1e-12)

        # reused cache == cache REBUILT fresh on the moved geometry
        plan.nearfield_cache[] = nothing
        FastMultipole.build_nearfield_cache!(plan, (body,), (body,))
        y1_rebuilt = copy(pnl._apply_dirichlet_G!(body, x, backend, scratch;
            plan_slot=slot, cache_nearfield=true))
        @test isapprox(y1, y1_rebuilt; rtol=1e-12)

        # cross-check against a fresh-tree uncached apply at backend
        # truncation level (a fresh tree on moved geometry has different
        # topology, so 1e-12 is not expected here)
        y1_fresh = copy(pnl._apply_dirichlet_G!(body, x, backend, scratch))
        @test isapprox(y1, y1_fresh; rtol=1e-6)

        # --- persistent_plan solve-level plumbing --------------------------
        function fresh_diamond()
            b = make_dirichlet_diamond_body(nspan=40)
            b.velocity .= 0
            b.velocity[1, :] .= 1.0
            return b
        end
        kw = (; backend, method=:gmres, atol=1e-12, rtol=1e-10, itmax=200)

        body_p = fresh_diamond()
        solver_p = pnl.KrylovSolver(body_p; kw..., cache_nearfield=true,
            persistent_plan=true)
        @test solver_p.cache_tree === true    # implied by persistent_plan
        pnl.solve!(body_p, solver_p)
        @test solver_p.kop.plan_slot[] !== nothing  # plan SURVIVED the solve
        plan_p = solver_p.kop.plan_slot[][1]
        mu0 = copy(body_p.strength)

        # second solve on frozen geometry reuses the same plan object
        pnl.solve!(body_p, solver_p)
        @test solver_p.kop.plan_slot[][1] === plan_p
        @test body_p.strength == mu0          # deterministic replay

        # rigid motion mirrored by the hook: the co-rotated problem (BC
        # velocity rotates with the body) must reproduce the same mu through
        # the persistent transformed plan
        move!(body_p)
        body_p.velocity .= Rrot * body_p.velocity
        pnl.transform_solver_geometry!(solver_p, body_p, Rrot, t_affine)
        @test solver_p.kop.plan_slot[][1] === plan_p  # transformed, not dropped
        pnl.solve!(body_p, solver_p)
        @test solver_p.kop.plan_slot[][1] === plan_p
        @test any(!iszero, body_p.strength[:, 2])
        @test isapprox(body_p.strength, mu0; rtol=1e-9)

        # --- FGS hook plumbing ---------------------------------------------
        body_f = fresh_diamond()
        solver_f = pnl.FGSSolver(body_f; expansion_order=4,
            multipole_acceptance=0.5, leaf_size=50, max_iterations=2,
            tolerance=1e-3, verbose=false)
        @test solver_f.fgs.transformed[] === false
        pnl.transform_solver_geometry!(solver_f, body_f, Rrot, t_affine)
        @test solver_f.fgs.transformed[] === true

        # transform_body_solvers! skips identity deltas and forwards real ones
        body_g = fresh_diamond()
        solver_g = pnl.FGSSolver(body_g; expansion_order=4,
            multipole_acceptance=0.5, leaf_size=50, max_iterations=2,
            tolerance=1e-3, verbose=false)
        identity_delta = pnl._identity_transforms(1)
        pnl.transform_body_solvers!((solver_g,), (body_g,), identity_delta)
        @test solver_g.fgs.transformed[] === false
        pnl.transform_body_solvers!((solver_g,), (body_g,),
            [(FastMultipole.SMatrix{3,3,Float64,9}(Rrot),
              FastMultipole.SVector{3,Float64}(t_affine))])
        @test solver_g.fgs.transformed[] === true

        # --- kinematics deltas match the applied motion --------------------
        body_k = fresh_diamond()
        nodes_before = copy(body_k.nodes)
        frames = [pnl.ReferenceFrame(
            FastMultipole.SVector{3}(0.1, -0.2, 0.3),
            FastMultipole.SVector{3}(0.4, 0.1, -0.2),
            FastMultipole.SVector{3}(axis...),
            2.5,
            FastMultipole.SMatrix{3,3,Float64,9}(1.0,0,0,0,1.0,0,0,0,1.0),
            FastMultipole.SMatrix{3,3,Float64,9}(1.0,0,0,0,1.0,0,0,0,1.0),
            "rotor", 0, Int[], [1])]
        dt = 0.07
        deltas = pnl.propagate_kinematics!((body_k,), frames, dt)
        Rk, tk = deltas[1]
        @test !(Rk ≈ FastMultipole.SMatrix{3,3,Float64,9}(1.0,0,0,0,1.0,0,0,0,1.0))
        @test maximum(abs.(Rk * nodes_before .+ tk .- body_k.nodes)) < 1e-12

        # --- guards + metadata ---------------------------------------------
        @test_throws ArgumentError pnl.KrylovSolver(fresh_diamond();
            backend=pnl.DirectBackend(), persistent_plan=true)
        md = pnl._solver_metadata_dict(solver_p)
        @test md["persistent_plan"] === true
        @test pnl._solver_metadata_dict(pnl.KrylovSolver(fresh_diamond();
            kw...))["persistent_plan"] === false
    end

    @testset "FGSSolver sweep_order plumbing (021 Phase 2b)" begin
        body = make_sphere_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0

        # default stays lexicographic; :colored constructs, solves, and reports
        solver_lex = pnl.FGSSolver(body; expansion_order=6, leaf_size=50,
            multipole_acceptance=0.4, max_iterations=100, inner_iterations=2,
            tolerance=1e-8, verbose=false)
        @test solver_lex.sweep_order === :lexicographic
        @test pnl._solver_metadata_dict(solver_lex)["sweep_order"] == "lexicographic"

        solver_col = pnl.FGSSolver(body; expansion_order=6, leaf_size=50,
            multipole_acceptance=0.4, max_iterations=100, inner_iterations=2,
            tolerance=1e-8, verbose=false, sweep_order=:colored)
        @test solver_col.sweep_order === :colored
        @test length(solver_col.fgs.leaves_by_color) > 1   # premise: coloring built
        pnl.solve!(body, solver_col)
        @test any(abs.(body.strength[:, 1]) .> 0)
        @test pnl._solver_metadata_dict(solver_col)["sweep_order"] == "colored"

        P = pnl.FGSPreconditioner(body; sweeps=2, inner_iterations=2,
            expansion_order=6, multipole_acceptance=0.4, leaf_size=50,
            sweep_order=:colored)
        @test pnl._preconditioner_metadata_dict(P)["sweep_order"] == "colored"
    end

    @testset "KrylovCoupled warmstart" begin
        body1 = make_octa_source_body()
        body2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        for body in (body1, body2)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
        end
        bodies = (body1, body2)
        solver = pnl.KrylovCoupled(bodies; backend=pnl.DirectBackend(), atol=1e-10, rtol=1e-10, itmax=50, warmstart=true)

        ws = solver.workspace
        pnl.solve!(bodies, solver)
        niter_cold = solver.niter
        strengths1 = copy(body1.strength[:, 1])

        pnl.solve!(bodies, solver)
        @test solver.workspace === ws
        @test solver.niter <= niter_cold
        @test isapprox(body1.strength[:, 1], strengths1; atol=1e-8)
        assert_boundary_residuals(bodies; tangency_atol=1e-7)
    end

    @testset "FGSSolver construction" begin
        body1 = make_octa_source_body()
        body2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        solver1 = pnl.FGSSolver(body1)
        solver_uncached = pnl.FGSSolver(body2; cache_leaf_lu=false)
        @test solver1.max_iterations == 100
        @test solver1.tolerance == 1e-6
        @test solver1.cache_leaf_lu
        @test solver1.fgs.cache_leaf_lu
        @test solver1.fgs.leaf_lu_cache !== nothing
        @test !solver_uncached.cache_leaf_lu
        @test !solver_uncached.fgs.cache_leaf_lu
        @test solver_uncached.fgs.leaf_lu_cache === nothing
        @test pnl._solver_metadata_dict(solver1)["cache_leaf_lu"] == true
        @test pnl._solver_metadata_dict(solver_uncached)["cache_leaf_lu"] == false
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
        pnl.solve!(body_fgs, fgs)

        @test body_fgs.strength[:, 1] ≈ body_backslash.strength[:, 1] atol=1e-12
        @test body_fgs.strength[:, 2] ≈ body_backslash.strength[:, 2] atol=1e-12

        neumann = make_octa_source_body()
        neumann.strength[:, 1] .= 3.0
        pnl.set_strengths!(neumann)
        @test all(iszero, neumann.strength[:, 1])
    end

    @testset "FGSSolver Neumann RigidWakeBody (VortexRing, 1 column)" begin
        # Regression for the missing FastMultipole.strength_to_value overload:
        # FGSSolver (and FGSPreconditioner) crashed on the uncapped Neumann
        # rotor body type RigidWakeBody{VortexRing,1,_,false} (021 Phase 0).
        body_ref = make_plate_vortex_body()
        body_fgs = make_plate_vortex_body()
        for body in (body_ref, body_fgs)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
            body.velocity[3, :] .= 0.2
            body.strength .= 0
        end

        solver_ref = pnl.Backslash(body_ref)
        pnl.solve!(body_ref, solver_ref)

        solver_fgs = pnl.FGSSolver(body_fgs; leaf_size=10000, tolerance=1e-12, max_iterations=5)
        pnl.solve!(body_fgs, solver_fgs)

        @test isapprox(body_fgs.strength[:, 1], body_ref.strength[:, 1]; atol=1e-10)
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

    @testset "Coupled Dirichlet rhs idempotence (regression)" begin
        # calc_bc_dirichlet used to ACCUMULATE (`RHS .-= potential`) into the
        # coupled solvers' persistent rhs, which is never zeroed between
        # solves: every solve after the first doubled the Dirichlet rows
        # (x -> 2 x_ref). Found via 021 Phase 1 (2026-08-14).
        make_dbc_body() = begin
            body = pnl.NonLiftingBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(
                copy(NODES_OCT), copy(CELLS_OCT); DBC=true)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
            body.velocity[3, :] .= 0.2
            pnl.calc_normals!(body)
            body
        end

        ref_body = make_dbc_body()
        pnl.solve!(ref_body, pnl.Backslash(ref_body))
        x_ref = copy(ref_body.strength[:, 2])
        @test any(abs.(x_ref) .> 0)

        body = make_dbc_body()
        solver = pnl.BackslashCoupled((body,))
        pnl.solve!((body,), solver)
        x1 = copy(body.strength[:, 2])
        @test isapprox(x1, x_ref; atol=1e-10)
        body.strength[:, 2] .= 0
        pnl.solve!((body,), solver)
        @test isapprox(copy(body.strength[:, 2]), x1; atol=1e-12)

        # KrylovCoupled routes its rhs through the same calc_bc_dirichlet path
        kbody = make_dbc_body()
        ksolver = pnl.KrylovCoupled((kbody,); backend=pnl.DirectBackend(),
                                    atol=1e-12, rtol=1e-12, itmax=50)
        pnl.solve!((kbody,), ksolver)
        k1 = copy(kbody.strength[:, 2])
        @test isapprox(k1, x_ref; atol=1e-8)
        kbody.strength[:, 2] .= 0
        pnl.solve!((kbody,), ksolver)
        @test isapprox(copy(kbody.strength[:, 2]), k1; atol=1e-8)
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
