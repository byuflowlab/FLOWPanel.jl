using Test
import FLOWPanel as pnl

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

mutable struct SimNoopSolver <: pnl.AbstractSolver
    backend_seen::Any
end
SimNoopSolver() = SimNoopSolver(nothing)

function _record_sim_noop_backend!(solver::SimNoopSolver, backend)
    solver.backend_seen = backend
    return nothing
end

pnl.solve!(::pnl.AbstractBody{<:Any, <:Any, <:Any, true}, solver::SimNoopSolver; backend=nothing, kwargs...) =
    _record_sim_noop_backend!(solver, backend)
pnl.solve!(::pnl.AbstractBody{<:Any, <:Any, <:Any, false}, solver::SimNoopSolver; backend=nothing, kwargs...) =
    _record_sim_noop_backend!(solver, backend)

mutable struct SimCoupledNoopSolver <: pnl.AbstractSolver
    called::Bool
    backend_seen::Any
end
SimCoupledNoopSolver() = SimCoupledNoopSolver(false, nothing)

function pnl.solve!(::Tuple, solver::SimCoupledNoopSolver; backend=nothing, kwargs...)
    solver.called = true
    solver.backend_seen = backend
    return nothing
end

mutable struct SimMarkerBackend <: pnl.AbstractBackend
    calls::Int
end
SimMarkerBackend() = SimMarkerBackend(0)

function pnl.influence!(::Tuple, ::Tuple, backend::SimMarkerBackend; kwargs...)
    backend.calls += 1
    return nothing
end

@testset "simulate! input normalization" begin
    Uinf = t -> [1.0, 0.0, 0.0]
    maneuver = (frames, systems, wakes, t) -> nothing
    t_range = [0.0, 0.05]

    body = make_plate_vortex_body()
    wake = pnl.PanelWake(body; nwakerows=2)
    frames = pnl.ReferenceFrame(body)
    solver = SimNoopSolver()
    pnl.simulate!(body, wake, frames, maneuver, Uinf, t_range;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        path=nothing,
    )
    @test solver.backend_seen isa pnl.DirectBackend

    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    solver = SimNoopSolver()
    pnl.simulate!(body, nothing, frames, maneuver, Uinf, t_range;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        path=nothing,
    )
    @test solver.backend_seen isa pnl.DirectBackend

    body = make_plate_vortex_body()
    wake = pnl.PanelWake(body; nwakerows=2)
    frames = pnl.ReferenceFrame(body)
    @test_throws ArgumentError pnl.simulate!(body, wake, frames, maneuver, Uinf, t_range;
        body_solvers=(SimNoopSolver(),),
        backend=pnl.DirectBackend(),
        path=nothing,
    )

    body1 = make_plate_vortex_body()
    body2 = make_plate_vortex_body()
    body2.nodes[2, :] .+= 2.0
    pnl.calc_normals!(body2)
    pnl.calc_controlpoints!(body2)
    wake1 = pnl.PanelWake(body1; nwakerows=2)
    wake2 = pnl.PanelWake(body2; nwakerows=2)
    frames = pnl.ReferenceFrame(body1)
    solver = SimCoupledNoopSolver()
    solve_backend = pnl.DirectBackend()
    pnl.simulate!((body1, body2), (wake1, wake2), frames, maneuver, Uinf, t_range;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        backend_solve=solve_backend,
        path=nothing,
    )
    @test solver.called
    @test solver.backend_seen === solve_backend
end

@testset "simulate! backend split and validation" begin
    Uinf = t -> [1.0, 0.0, 0.0]
    maneuver = (frames, systems, wakes, t) -> nothing
    t_range = [0.0, 0.05]

    body = make_plate_vortex_body()
    wake = pnl.PanelWake(body; nwakerows=2)
    frames = pnl.ReferenceFrame(body)
    solver = SimNoopSolver()
    backend_wake = SimMarkerBackend()
    backend_system = SimMarkerBackend()
    backend_solve = :solve_backend

    pnl.simulate!(body, wake, frames, maneuver, Uinf, t_range;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        backend_wake=backend_wake,
        backend_solve=backend_solve,
        backend_system=backend_system,
        path=nothing,
    )
    @test backend_wake.calls == 2
    @test backend_system.calls == 2
    @test solver.backend_seen === backend_solve

    body1 = make_plate_vortex_body()
    body2 = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body1)
    @test_throws ArgumentError pnl.simulate!((body1, body2), (nothing,), frames, maneuver, Uinf, t_range;
        body_solvers=(SimNoopSolver(), SimNoopSolver()),
        backend=pnl.DirectBackend(),
        path=nothing,
    )

    @test_throws ArgumentError pnl.simulate!((body1, body2), (nothing, nothing), frames, maneuver, Uinf, t_range;
        body_solvers=(SimNoopSolver(),),
        backend=pnl.DirectBackend(),
        path=nothing,
    )

    @test_throws ArgumentError pnl.simulate!(body1, nothing, frames, maneuver, Uinf, t_range;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        backend_solve=(pnl.DirectBackend(),),
        path=nothing,
    )
end
