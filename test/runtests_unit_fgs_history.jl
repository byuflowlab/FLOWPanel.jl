using Test
import FLOWPanel as pnl

include("test_helpers.jl")

@testset verbose=true "FGSSolver history & projection" begin

    @testset "construction with history kwargs" begin
        body = make_octa_source_body()
        solver = pnl.FGSSolver(body; solution_history_length=4,
                               project_solution=true, project_solution_order=2)

        @test solver.solution_history_length == 4
        @test solver.solution_history_nsaved == 0
        @test solver.project_solution == true
        @test solver.project_solution_order == 2
        @test size(solver.solution_history) == (body.ncells, size(body.strength, 2), 4)

        # default kwargs disable the feature with no allocation cost
        solver_off = pnl.FGSSolver(body)
        @test solver_off.solution_history_length == 0
        @test solver_off.project_solution == false
        @test size(solver_off.solution_history, 3) == 0
    end

    @testset "save_solution! shifts and writes slot 1" begin
        body = make_octa_source_body()
        solver = pnl.FGSSolver(body; solution_history_length=3)
        H = solver.solution_history

        # first save
        body.strength .= 1.0
        pnl.save_solution!(body, solver)
        @test solver.solution_history_nsaved == 1
        @test all(H[:, :, 1] .== 1.0)

        # second save: prior contents move to slot 2
        body.strength .= 2.0
        pnl.save_solution!(body, solver)
        @test solver.solution_history_nsaved == 2
        @test all(H[:, :, 1] .== 2.0)
        @test all(H[:, :, 2] .== 1.0)

        # third save: full
        body.strength .= 3.0
        pnl.save_solution!(body, solver)
        @test solver.solution_history_nsaved == 3
        @test all(H[:, :, 1] .== 3.0)
        @test all(H[:, :, 2] .== 2.0)
        @test all(H[:, :, 3] .== 1.0)

        # fourth save: oldest (1.0) drops off, nsaved clamped to NT
        body.strength .= 4.0
        pnl.save_solution!(body, solver)
        @test solver.solution_history_nsaved == 3
        @test all(H[:, :, 1] .== 4.0)
        @test all(H[:, :, 2] .== 3.0)
        @test all(H[:, :, 3] .== 2.0)
    end

    @testset "project_solution! recovers polynomial" begin
        body = make_octa_source_body()

        # quadratic in time: s(t) = a + b*t + c*t^2
        a, b, c = 0.7, -0.3, 0.15
        s(t) = a + b*t + c*t^2

        solver = pnl.FGSSolver(body; solution_history_length=4,
                               project_solution=true, project_solution_order=2)
        H = solver.solution_history
        # slot 1 = most recent → t=3, slot 2 → t=2, slot 3 → t=1
        H[:, :, 1] .= s(3)
        H[:, :, 2] .= s(2)
        H[:, :, 3] .= s(1)
        solver.solution_history_nsaved = 3

        body.strength .= 0
        pnl.project_solution!(body, solver)
        @test all(isapprox.(body.strength, s(4); atol=1e-12))

        # linear case: with only 2 saved points, falls back to order=1 even if requested 2
        solver.solution_history_nsaved = 2
        body.strength .= 0
        pnl.project_solution!(body, solver)
        # 2*s(3) - s(2) is the linear extrapolation prediction at t=4
        expected = 2 * s(3) - s(2)
        @test all(isapprox.(body.strength, expected; atol=1e-12))

        # nsaved < 2 → no-op, body.strength untouched
        solver.solution_history_nsaved = 1
        body.strength .= 99.0
        pnl.project_solution!(body, solver)
        @test all(body.strength .== 99.0)

        solver.solution_history_nsaved = 0
        body.strength .= 77.0
        pnl.project_solution!(body, solver)
        @test all(body.strength .== 77.0)
    end

    @testset "project_solution! disabled when project_solution=false" begin
        body = make_octa_source_body()
        solver = pnl.FGSSolver(body; solution_history_length=3,
                               project_solution=false, project_solution_order=2)
        H = solver.solution_history
        H[:, :, 1] .= 5.0
        H[:, :, 2] .= 3.0
        solver.solution_history_nsaved = 2

        body.strength .= 42.0
        pnl.project_solution!(body, solver)
        @test all(body.strength .== 42.0)
    end

    @testset "allocation-free hot path" begin
        body = make_octa_source_body()
        solver = pnl.FGSSolver(body; solution_history_length=4,
                               project_solution=true, project_solution_order=2)

        # populate history so projection has work to do
        H = solver.solution_history
        H[:, :, 1] .= 1.0
        H[:, :, 2] .= 2.0
        H[:, :, 3] .= 3.0
        solver.solution_history_nsaved = 3
        body.strength .= 0

        # warm up to compile both paths
        pnl.project_solution!(body, solver)
        pnl.save_solution!(body, solver)

        # measure
        @test (@allocated pnl.project_solution!(body, solver)) == 0
        @test (@allocated pnl.save_solution!(body, solver)) == 0
    end

end
