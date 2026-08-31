#=##############################################################################
Unit tests for the solver instrumentation layer (BRAINSTORM 021 Phase 0):
shared operator application `apply_G!`, the true-residual evaluator
`true_residual!` (021 ruling 3: identical code path for every solver, verified
against a dense hand computation), `ConvergenceHistory`, and the
`BackslashCoupled` first-solve build + timer fields.
=###############################################################################

using Test
import FLOWPanel as pnl
using LinearAlgebra: norm
import Random

if !isdefined(@__MODULE__, :make_octa_source_body)
    include("test_helpers.jl")
end

"Assemble the dense influence matrix the way `Backslash` does."
function dense_G(body)
    n = body.ncells
    G = zeros(n, n)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    pnl._G!(G, body, body; core_size=body.core_size_panel)
    return G
end

@testset verbose=true "Instrumentation" begin

    rng = Random.MersenneTwister(20260811)

    @testset "apply_G! matches dense G*x (Neumann)" begin
        body = make_octa_source_body()
        G = dense_G(body)
        x = randn(rng, body.ncells)
        y = zeros(body.ncells)

        pnl.apply_G!(y, body, x)
        @test isapprox(y, G * x; atol=1e-12)

        # FMM backend agrees with the dense operator to FMM accuracy
        y_fmm = zeros(body.ncells)
        pnl.apply_G!(y_fmm, body, x, pnl.FastMultipoleBackend())
        @test isapprox(y_fmm, G * x; rtol=1e-6, atol=1e-8)
    end

    @testset "apply_G! matches dense G*x (Dirichlet)" begin
        body = make_dirichlet_diamond_body()
        G = dense_G(body)
        x = randn(rng, body.ncells)
        y = zeros(body.ncells)

        pnl.apply_G!(y, body, x)
        @test isapprox(y, G * x; atol=1e-12)

        y_fmm = zeros(body.ncells)
        pnl.apply_G!(y_fmm, body, x, pnl.FastMultipoleBackend())
        @test isapprox(y_fmm, G * x; rtol=1e-6, atol=1e-8)
    end

    @testset "true_residual! matches dense hand computation and is side-effect-free" begin
        for make_body in (make_octa_source_body, make_dirichlet_diamond_body)
            body = make_body()
            G = dense_G(body)
            n = body.ncells
            x = randn(rng, n)
            b = randn(rng, n)

            # snapshot state, evaluate, verify restoration
            body.velocity .= randn(rng, 3, n)
            body.potential .= randn(rng, n)
            strength0 = copy(body.strength)
            velocity0 = copy(body.velocity)
            potential0 = copy(body.potential)

            r = zeros(n)
            rms, rmax = pnl.true_residual!(r, body, x, b)

            r_dense = G * x - b
            @test isapprox(r, r_dense; atol=1e-12)
            @test isapprox(rms, norm(r_dense) / sqrt(n); atol=1e-14)
            @test isapprox(rmax, maximum(abs, r_dense); atol=1e-14)

            @test body.strength == strength0
            @test body.velocity == velocity0
            @test body.potential == potential0
        end
    end

    @testset "Backslash solution has round-off true residual (Neumann)" begin
        body = make_octa_source_body()
        solver = pnl.Backslash(body)
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        pnl.solve!(body, solver)

        n = body.ncells
        b = zeros(n)
        pnl.assemble_rhs!(b, body)
        r = zeros(n)
        rms, rmax = pnl.true_residual!(r, body, vec(body.strength[:, 1]), b)
        @test rms <= 1e-10 * (norm(b) / sqrt(n))
    end

    @testset "Backslash solution has round-off true residual (Dirichlet)" begin
        # Replicate the Dirichlet solve! internals so the RHS can be captured
        # at the instant the solver is invoked (b = -potential workspace).
        body = make_dirichlet_diamond_body()
        solver = pnl.Backslash(body)
        n = body.ncells
        body.velocity .= 0
        body.velocity[1, :] .= 1.0

        pnl.set_strengths!(body)
        body.potential .= 0
        pnl.influence!(body, body, pnl.DirectBackend(); scalar_potential=true, velocity=false)
        b = zeros(n)
        pnl.assemble_rhs!(b, body)
        pnl._solve!(body, solver)

        r = zeros(n)
        rms, rmax = pnl.true_residual!(r, body, vec(body.strength[:, 2]), b)
        @test rms <= 1e-10 * max(norm(b) / sqrt(n), 1e-16)
    end

    @testset "ConvergenceHistory record/reset" begin
        h = pnl.ConvergenceHistory(:krylov_precnorm)
        pnl.reset!(h)
        @test length(h) == 0
        @test h.t0_ns > 0

        for (i, res) in enumerate([1.0, 0.1, 0.01])
            pnl.record!(h, i, res)
        end
        @test length(h) == 3
        @test h.iter == [1, 2, 3]
        @test h.residual_internal == [1.0, 0.1, 0.01]
        @test issorted(h.t_ns)
        @test all(h.t_ns .>= h.t0_ns)
        @test all(pnl.elapsed_seconds(h) .>= 0)

        pnl.reset!(h; metric=:fgs_maxabs)
        @test length(h) == 0
        @test h.metric == :fgs_maxabs
    end

    @testset "BackslashCoupled builds on first solve (no update_G)" begin
        # Regression for the identity-LU bug: inside simulate!, solve! is
        # called without update_G; the solver must still assemble G on its
        # first solve rather than "solving" against the ctor's dummy identity.
        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0

        body_ref = make_octa_source_body()
        body_ref.velocity .= 0
        body_ref.velocity[1, :] .= 1.0
        solver_ref = pnl.Backslash(body_ref)
        pnl.solve!(body_ref, solver_ref)

        solver = pnl.BackslashCoupled((body,))
        @test !solver.built
        pnl.solve!((body,), solver)   # no update_G
        @test solver.built
        @test isapprox(body.strength[:, 1], body_ref.strength[:, 1]; atol=1e-10)
        @test solver.t_build > 0
        @test solver.t_solve > 0

        # second solve reuses the cached factorization and reproduces the answer
        t_build_first = solver.t_build
        pnl.solve!((body,), solver)
        @test isapprox(body.strength[:, 1], body_ref.strength[:, 1]; atol=1e-10)
        @test solver.t_build == t_build_first   # untouched: no rebuild happened
    end

end
