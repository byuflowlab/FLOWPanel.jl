using Test
import FLOWPanel as pnl
import GeometricTools as gt
using LinearAlgebra: diag

@testset verbose=true "Solvers" begin
    @testset "Backslash dispatch" begin
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
        @test pnl.Backslash(nlb) isa pnl.BackslashNeumann
        @test pnl.Backslash(rwb) isa pnl.BackslashDirichlet
    end

    @testset "BackslashNeumann construction and solve" begin
        body = make_octa_source_body()
        solver = pnl.BackslashNeumann(body)
        @test size(solver.G) == (body.ncells, body.ncells)
        @test all(diag(solver.G) .!= 0)
        @test solver.Glu !== nothing
        @test length(solver.rhs) == body.ncells

        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        solver = pnl.Backslash(body)
        pnl.solve!(body, solver)
        @test any(abs.(body.strength[:, 1]) .> 0)

        target = deepcopy(body)
        target.velocity .= 0
        target.potential .= 0
        pnl.influence!(target, body, pnl.DirectBackend(); velocity=true)
        target.velocity[1, :] .+= 1.0
        residual = abs.(sum(target.velocity .* body.normals; dims=1))
        @test maximum(residual) < 0.2
    end

    @testset "BackslashDirichlet construction and solve" begin
        body = make_nonlifting(Union{pnl.ConstantSource, pnl.ConstantDoublet}; DBC=true)
        pnl.calc_normals!(body)
        normals = copy(body.normals)
        solver = pnl.BackslashDirichlet(body)
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
        solver = pnl.BackslashDirichlet(body)
        pnl.solve!(body, solver)
        @test any(abs.(body.strength[:, 2]) .> 0)
        @test isapprox(vec(body.strength[:, 1]), expected_sigma; atol=1e-12)
    end

    @testset "KrylovSolver construction" begin
        body = make_octa_source_body()
        solver = pnl.KrylovSolver(body)
        @test solver.method == :gmres
        @test solver.itmax == 20

        backend = pnl.DirectBackend()
        custom = pnl.KrylovSolver(body; atol=1e-8, rtol=1e-8, backend=backend)
        @test custom.atol == 1e-8
        @test custom.rtol == 1e-8
        @test custom.backend === backend
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
    end

    @testset "calc_elprescribe and numtype" begin
        body_source = make_nonlifting(pnl.ConstantSource)
        body_vortex_closed = pnl.NonLiftingBody{pnl.VortexRing}(copy(NODES_2TRI), copy(CELLS_2TRI); watertight=true)
        body_vortex_open = pnl.NonLiftingBody{pnl.VortexRing}(copy(NODES_2TRI), copy(CELLS_2TRI); watertight=false)
        @test pnl.calc_elprescribe(body_source) == Tuple{Int, Float64}[]
        @test pnl.calc_elprescribe(body_vortex_closed) == [(1, 0.0)]
        @test pnl.calc_elprescribe(body_vortex_open) == Tuple{Int, Float64}[]

        @test pnl.numtype(body_source) == Float64
        body32 = pnl.NonLiftingBody{pnl.ConstantSource}(Float32.(NODES_2TRI), copy(CELLS_2TRI))
        @test pnl.numtype(body32) == Float32
    end
end
