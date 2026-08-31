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

    # --- BRAINSTORM 021 Phase 3 -------------------------------------------
    # One tuple `solve!` is ONE top-level step, however many block-GS outer
    # iterations it takes. Before the fix, warm-start history advanced once per
    # OUTER ITERATION, so it held intra-step repeats and the extrapolated guess
    # could alternate around the answer instead of leading it.

    @testset "one-body tuple solves exactly once" begin
        fresh() = (b = make_octa_source_body(); b.velocity .= 0;
                   b.velocity[1, :] .= 1.0; b)

        body = fresh()
        solver = pnl.FGSSolver(body; tolerance=1e-10, max_iterations=100,
                               solution_history_length=4, project_solution=true,
                               project_solution_order=1)
        h = pnl.ConvergenceHistory()
        pnl.solve!((body,), (solver,); backend=(pnl.DirectBackend(),), history=h)

        @test length(h) == 1                       # one outer entry, by contract
        @test pnl.step_nsolves(solver) == 1
        @test pnl.step_niter_first(solver) == solver.niter
        @test pnl.step_solved(solver)
        @test solver.solution_history_nsaved == 1  # advanced exactly once

        pnl.solve!((body,), (solver,); backend=(pnl.DirectBackend(),))
        @test pnl.step_nsolves(solver) == 1
        @test solver.solution_history_nsaved == 2

        # same for the Krylov side, with extrapolation configured
        bodyK = fresh()
        solverK = pnl.KrylovSolver(bodyK; backend=pnl.DirectBackend(), atol=1e-12,
                                   rtol=1e-12, itmax=50, warmstart=true,
                                   warmstart_order=2)
        hK = pnl.ConvergenceHistory()
        pnl.solve!((bodyK,), (solverK,); backend=(pnl.DirectBackend(),), history=hK)
        @test length(hK) == 1
        @test pnl.step_nsolves(solverK) == 1
        @test solverK.x_history_nsaved == 1
    end

    @testset "two-body tuple: one history advance per step" begin
        make_pair() = begin
            b1 = make_octa_source_body()
            b2 = translated_nonlifting_target([3.0, 0.0, 0.0])
            for b in (b1, b2)
                b.velocity .= 0
                b.velocity[1, :] .= 1.0
            end
            (b1, b2)
        end

        body1, body2 = make_pair()
        solver1 = pnl.FGSSolver(body1; tolerance=1e-10, max_iterations=100,
                                solution_history_length=8)
        solver2 = pnl.FGSSolver(body2; tolerance=1e-10, max_iterations=100,
                                solution_history_length=8)
        h = pnl.ConvergenceHistory()
        pnl.solve!((body1, body2), (solver1, solver2);
                   backend=(pnl.DirectBackend(), pnl.DirectBackend()),
                   outer_tolerance=1e-10, history=h)

        @test length(h) >= 2                        # genuinely coupled: multiple
        for s in (solver1, solver2)                 # outer iterations
            @test s.solution_history_nsaved == 1    # ...but ONE history advance
            @test pnl.step_nsolves(s) == length(h)  # one per-body solve per outer
            @test pnl.step_niter_first(s) >= 0
            @test pnl.step_solved(s)
        end

        # a second step advances each history by exactly one again
        pnl.solve!((body1, body2), (solver1, solver2);
                   backend=(pnl.DirectBackend(), pnl.DirectBackend()),
                   outer_tolerance=1e-10)
        @test solver1.solution_history_nsaved == 2
        @test solver2.solution_history_nsaved == 2

        # Krylov: x_prev/x_history commit only at tuple completion, so a step
        # taking several outer iterations still advances the history by one.
        b1, b2 = make_pair()
        kw = (; backend=pnl.DirectBackend(), atol=1e-12, rtol=1e-12, itmax=50,
              warmstart=true, warmstart_order=2)
        k1 = pnl.KrylovSolver(b1; kw...)
        k2 = pnl.KrylovSolver(b2; kw...)
        hK = pnl.ConvergenceHistory()
        pnl.solve!((b1, b2), (k1, k2);
                   backend=(pnl.DirectBackend(), pnl.DirectBackend()),
                   outer_tolerance=1e-10, history=hK)
        @test length(hK) >= 2
        for s in (k1, k2)
            @test s.x_history_nsaved == 1
            @test pnl.step_nsolves(s) == length(hK)
        end
    end

    @testset "solver_optargs may not hijack the step boundary" begin
        body = make_octa_source_body()
        body.velocity .= 0
        body.velocity[1, :] .= 1.0
        solver = pnl.FGSSolver(body; tolerance=1e-8)
        @test_throws ArgumentError pnl.solve!((body,), (solver,);
            backend=(pnl.DirectBackend(),),
            solver_optargs=(; finalize_step=true))
    end

end
