using Test
import FLOWPanel as pnl
import FastMultipole
import FLOWVPM

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
    kwargs_seen::Vector{Any}
end
SimMarkerBackend() = SimMarkerBackend(0, Any[])

function pnl.influence!(::Tuple, ::Tuple, backend::SimMarkerBackend; kwargs...)
    backend.calls += 1
    push!(backend.kwargs_seen, kwargs)
    return nothing
end

mutable struct SimWakeRecorder <: pnl.AbstractMonitor
    wakes_seen::Any
    i_step_seen::Int
    dt_seen::Float64
end
SimWakeRecorder() = SimWakeRecorder(nothing, -1, NaN)

function (m::SimWakeRecorder)(systems, wakes, frames, uinf, i_step::Int, dt::Real)
    m.wakes_seen = wakes
    m.i_step_seen = i_step
    m.dt_seen = Float64(dt)
    return nothing
end

mutable struct SimDistributedForce <: pnl.AbstractMonitor
    force::Matrix{Float64}
    i_system::Int
end

pnl.monitor_provides(::SimDistributedForce) = (:F,)

function (m::SimDistributedForce)(systems, wakes, frames, uinf, i_step::Int, dt::Real)
    return nothing
end

function pnl._register_monitor_outputs!(ctx::pnl.MonitorContext, m::SimDistributedForce, systems_tuple::Tuple)
    pnl.monitor_register!(ctx, :F, m.i_system, m.force)
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
    @test all(!haskey(kwargs, :postcalc) || kwargs[:postcalc] == false for kwargs in backend_wake.kwargs_seen)
    @test backend_system.calls == 2
    @test all(haskey(kwargs, :direct_conditioning) for kwargs in backend_system.kwargs_seen)
    @test all(kwargs[:direct_conditioning] isa FastMultipole.DirectConditioningRule for kwargs in backend_system.kwargs_seen)
    @test solver.backend_seen === backend_solve

    body = make_plate_vortex_body()
    wake = pnl.PanelParticleWake(body; nwakerows=2, particle_kerneloffset=0.1,
        SFS=FLOWVPM.ConstantSFS(FLOWVPM.Estr_fmm))
    frames = pnl.ReferenceFrame(body)
    solver = SimNoopSolver()
    backend_wake = SimMarkerBackend()
    backend_system = SimMarkerBackend()
    pnl.simulate!(body, wake, frames, maneuver, Uinf, t_range;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        backend_wake=backend_wake,
        backend_system=backend_system,
        path=nothing,
    )
    @test backend_wake.calls == 2
    @test all(kwargs[:postcalc] == true for kwargs in backend_wake.kwargs_seen)
    @test backend_system.calls == 2
    @test body.kerneloffset == body.kerneloffset_targets
    @test body.kerneloffset_targets == 0.1

    body = make_plate_vortex_body()
    wake = pnl.PanelParticleWake(body; nwakerows=2, particle_kerneloffset=0.1,
        SFS=FLOWVPM.ConstantSFS(FLOWVPM.Estr_fmm))
    frames = pnl.ReferenceFrame(body)
    solver = SimNoopSolver()
    backend_system = SimMarkerBackend()
    pnl.simulate!(body, wake, frames, maneuver, Uinf, t_range;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        backend_system=backend_system,
        body_hessian_to_particles=true,
        path=nothing,
    )
    @test all(kwargs[:velocity_gradient][end] == true for kwargs in backend_system.kwargs_seen)

    body = make_plate_vortex_body()
    wake = pnl.PanelParticleWake(body; nwakerows=2, particle_kerneloffset=0.1,
        SFS=FLOWVPM.ConstantSFS(FLOWVPM.Estr_fmm))
    frames = pnl.ReferenceFrame(body)
    solver = SimNoopSolver()
    backend_wake = SimMarkerBackend()
    pnl.simulate!(body, wake, frames, maneuver, Uinf, t_range;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        backend_wake=backend_wake,
        particle_hessian_self=false,
        path=nothing,
    )
    @test backend_wake.calls > 2
    @test any(length(kwargs[:velocity_gradient]) == 1 && kwargs[:velocity_gradient][1] == false
        for kwargs in backend_wake.kwargs_seen)

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

@testset "steady! body-only solve and validation" begin
    uinf = [1.0, 0.0, 0.0]

    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    solver = SimNoopSolver()
    backend_system = SimMarkerBackend()
    backend_solve = :steady_solve_backend
    recorder = SimWakeRecorder()

    pnl.steady!(body, frames, uinf;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        backend_solve=backend_solve,
        backend_system=backend_system,
        monitors=(recorder,),
        i_run=2,
        dt=0.25,
    )
    @test solver.backend_seen === backend_solve
    @test backend_system.calls == 1
    @test all(haskey(kwargs, :direct_conditioning) for kwargs in backend_system.kwargs_seen)
    @test all(kwargs[:direct_conditioning] isa FastMultipole.DirectConditioningRule for kwargs in backend_system.kwargs_seen)
    @test body.kerneloffset == body.kerneloffset_targets
    @test recorder.wakes_seen === (nothing,)
    @test recorder.i_step_seen == 1
    @test recorder.dt_seen == 0.25

    body1 = make_plate_vortex_body()
    body2 = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body1)
    @test_throws ArgumentError pnl.steady!((body1, body2), frames, uinf;
        body_solvers=(SimNoopSolver(),),
        backend=pnl.DirectBackend(),
    )
    @test_throws ArgumentError pnl.steady!(body1, frames, uinf;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        backend_solve=(pnl.DirectBackend(),),
    )
    @test_throws ArgumentError pnl.steady!(body1, frames, uinf;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        i_run=0,
    )
    @test_throws ArgumentError pnl.steady!(body1, frames, uinf;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        dt=0.0,
    )
end

@testset "steady! monitors" begin
    uinf = [1.0, 0.0, 0.0]

    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    solver = SimNoopSolver()
    pressure = pnl.PressureBernoulli(1.0)
    force = pnl.ForceMonitor(3, 1; normalization=pnl.NoNormalization())
    force.force .= NaN
    force.moment .= NaN

    pnl.steady!(body, frames, uinf;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        backend_system=SimMarkerBackend(),
        monitors=(pressure, force),
        i_run=2,
        dt=1.0,
    )
    @test all(isnan, force.force[:, 1])
    @test all(isfinite, force.force[:, 2])
    @test all(isnan, force.force[:, 3])

    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    @test_throws ArgumentError pnl.steady!(body, frames, uinf;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        monitors=(pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization()),),
    )

    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    laplace = pnl.PressureLaplace((body,), 1.0; reference_panel=1)
    recorder = SimWakeRecorder()
    @test !body.needs_velocity_gradient[]
    pnl.steady!(body, frames, uinf;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        backend_system=SimMarkerBackend(),
        monitors=(laplace, recorder),
        i_run=1,
        dt=0.1,
    )
    @test body.needs_velocity_gradient[]
    @test recorder.wakes_seen === (nothing,)

    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    steady_outdir = mktempdir()
    steady_stale = joinpath(steady_outdir, "monitors", "steadyclean_monitor09_force_system1.csv")
    mkpath(dirname(steady_stale))
    write(steady_stale, "stale\n")
    pnl.steady!(body, frames, uinf;
        name="steadyclean",
        path=steady_outdir,
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        backend_system=SimMarkerBackend(),
        monitors=(),
        i_run=1,
        dt=0.1,
    )
    @test !isfile(steady_stale)

    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    outdir = mktempdir()
    csv = joinpath(outdir, "monitors", "moncsv_monitor02_spanwise_system1.csv")
    stale = joinpath(outdir, "monitors", "moncsv_monitor99_force_system1.csv")
    mkpath(dirname(stale))
    write(stale, "stale\n")
    provider = SimDistributedForce([0.0 0.0; 0.0 0.0; 1.0 1.0], 1)
    spanwise = pnl.SpanwiseLoadingMonitor(1, 1;
        components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]))
    pnl.simulate!(body, nothing, frames, (frames, systems, wakes, t) -> nothing,
        t -> uinf, [2.0, 2.25];
        name="moncsv",
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        backend_system=SimMarkerBackend(),
        monitors=(provider, spanwise),
        path=outdir,
    )
    rows = readlines(csv)
    @test startswith(rows[2], "0,2.0,1,")
    @test startswith(rows[3], "1,2.25,1,")
    @test !isfile(stale)

    stale_keep = joinpath(outdir, "monitors", "keep_monitor99_force_system1.csv")
    write(stale_keep, "stale\n")
    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    pnl.simulate!(body, nothing, frames, (frames, systems, wakes, t) -> nothing,
        t -> uinf, [2.0, 2.25];
        name="keep",
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        backend_system=SimMarkerBackend(),
        monitors=(provider, spanwise),
        path=outdir,
        clean_files=false,
    )
    @test isfile(stale_keep)

    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    spanwise = pnl.SpanwiseLoadingMonitor(1, 1;
        components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]))
    pnl.simulate!(body, nothing, frames, (frames, systems, wakes, t) -> nothing,
        t -> uinf, [2.0, 2.25];
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        backend_system=SimMarkerBackend(),
        monitors=(provider, spanwise),
        path=nothing,
    )
end
