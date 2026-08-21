#=##############################################################################
Unit tests for iterative-solver convergence-history capture (BRAINSTORM 021
Phase 0 W4 / ruling 4): the FastGaussSeidel per-outer-iteration callback
(non-breaking kwarg added to FastMultipole `solve!`, decision D2), FLOWPanel's
`_solve_history!` wrapper, and the tuple block-Gauss–Seidel `history` kwarg.
=###############################################################################

using Test
import FLOWPanel as pnl

if !isdefined(@__MODULE__, :make_octa_source_body)
    include("test_helpers.jl")
end

@testset verbose=true "Solver convergence history" begin

    @testset "FGS callback leaves the solve bit-identical" begin
        body_plain = make_octa_source_body()
        body_hist = make_octa_source_body()
        for body in (body_plain, body_hist)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
        end

        make_solver(body) = pnl.FGSSolver(body; leaf_size=10000, tolerance=1e-12,
                                          max_iterations=5)

        solver_plain = make_solver(body_plain)
        pnl.solve!(body_plain, solver_plain)

        solver_hist = make_solver(body_hist)
        h = pnl.ConvergenceHistory()
        pnl._solve_history!(body_hist, solver_hist, h)

        @test body_hist.strength[:, 1] == body_plain.strength[:, 1]
    end

    @testset "FGS history contents" begin
        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        solver = pnl.FGSSolver(body; leaf_size=10000, tolerance=1e-12, max_iterations=5)

        h = pnl.ConvergenceHistory()
        pnl._solve_history!(body, solver, h)

        @test h.metric == :fgs_maxabs
        @test 1 <= length(h) <= solver.max_iterations
        @test h.iter == collect(1:length(h))
        @test issorted(h.t_ns)
        @test all(h.t_ns .>= h.t0_ns)
        # the recorded residual is the loop's own tolerance metric: the last
        # entry must be at (or below) the first for a converging solve
        @test h.residual_internal[end] <= h.residual_internal[1]
        # single-leaf FGS solves exactly in one sweep: the loop records the
        # initial residual, then the converged one, and breaks
        @test h.residual_internal[end] <= solver.tolerance
    end

    @testset "block-GS history kwarg" begin
        body1 = make_octa_source_body()
        body2 = translated_nonlifting_target([3.0, 0.0, 0.0])
        for body in (body1, body2)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
        end
        solvers = (pnl.Backslash(body1), pnl.Backslash(body2))

        h = pnl.ConvergenceHistory()
        pnl.solve!((body1, body2), solvers;
                   backend=(pnl.DirectBackend(), pnl.DirectBackend()),
                   outer_tolerance=1e-10, history=h)

        @test h.metric == :blockgs_maxdelta
        @test length(h) >= 2                       # at least one coupling update + convergence check
        @test h.iter == collect(1:length(h))
        @test issorted(h.t_ns)
        @test h.residual_internal[end] < 1e-10     # converged on max |Δstrength|

        # history capture must not perturb the solution
        body1_ref = make_octa_source_body()
        body2_ref = translated_nonlifting_target([3.0, 0.0, 0.0])
        for body in (body1_ref, body2_ref)
            body.velocity .= 0
            body.velocity[1, :] .= 1.0
        end
        solvers_ref = (pnl.Backslash(body1_ref), pnl.Backslash(body2_ref))
        pnl.solve!((body1_ref, body2_ref), solvers_ref;
                   backend=(pnl.DirectBackend(), pnl.DirectBackend()),
                   outer_tolerance=1e-10)
        @test body1.strength[:, 1] == body1_ref.strength[:, 1]
        @test body2.strength[:, 1] == body2_ref.strength[:, 1]
    end

end
