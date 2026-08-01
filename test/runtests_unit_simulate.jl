using Test
import FLOWPanel as pnl
import FastMultipole
import FLOWVPM
using LinearAlgebra: dot, norm, cross

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
        grad_mu_options=(; basis=:tri),
    )
    @test solver.backend_seen isa pnl.DirectBackend

    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body)
    solver = SimNoopSolver()
    pnl.simulate!(body, nothing, frames, maneuver, Uinf, t_range;
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        path=nothing,
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
    laplace = pnl.PressureLaplace((body,), 1.0; reference_panel=1,
        gradient_mode=:corrected_hessian)
    recorder = SimWakeRecorder()
    @test !body.needs_velocity_gradient[]
    pnl.steady!(body, frames, uinf;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        backend_system=SimMarkerBackend(),
        monitors=(laplace, recorder),
        i_run=1,
        dt=0.1,
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
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
        grad_mu_options=(; basis=:tri),
    )
end

@testset "steady! moving-frame steady Bernoulli uses relative trace" begin
    # Regression: steady PressureBernoulli on a body in a rotating frame must
    # exclude the kinematic velocity from the kinetic-energy term (the
    # ef1fe1e inertial form cancelled first-order blade loading).
    rho = 1.3
    uinf = [1.0, 0.0, 0.0]
    body = make_plate_vortex_body()
    frames = pnl.ReferenceFrame(body;
        ω_axis=FastMultipole.SVector{3}(0.0, 0.0, 1.0), ω=5.0)
    pressure = pnl.PressureBernoulli(rho; backend=pnl.DirectBackend(),
        correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())

    pnl.steady!(body, frames, uinf;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        backend_system=SimMarkerBackend(),
        monitors=(pressure, force),
        grad_mu_options=(; basis=:tri),
    )

    # kinematic_velocity! ran: body.velocity is body-relative, and the rotation
    # produced a nonzero kinematic velocity.
    @test any(!iszero, body.velocity_kinematic)

    # Steady pressure matches the relative-trace formula panelwise.
    Uinf2 = norm(uinf)^2
    expected = map(1:body.ncells) do p
        n = body.normals[:, p]
        q = body.velocity[:, p]
        qt = q - dot(q, n) * n
        0.5 * rho * (Uinf2 - norm(qt)^2)
    end
    @test pressure.pressure[1] ≈ expected
    @test all(isfinite, force.force[:, 1])
    @test any(!iszero, force.force[:, 1])

    # A rerun with the kinematic velocity zeroed but identical relative
    # velocity gives the same steady pressure.
    saved_velocity = copy(body.velocity)
    body.velocity_kinematic .= 0.0
    body.velocity .= saved_velocity
    static_pressure = pnl.PressureBernoulli(rho; backend=pnl.DirectBackend(),
        correct_kuttacondition=false)
    static_pressure((body,), (nothing,), pnl.ReferenceFrame(body), uinf, 0, 1.0)
    @test static_pressure.pressure[1] ≈ pressure.pressure[1]
end

@testset "kinematic Das arc construction" begin
    SV3 = FastMultipole.SVector{3, Float64}

    # --- unit level: _rigid_back_displacement -------------------------------

    origin = SV3(0.0, 0.0, 0.0)
    te = SV3(0.75, 0.0, 0.0)

    # Pure translation: the arc form must reduce to the tangent form exactly,
    # so translating bodies (the wing cases) are bit-for-bit unaffected.
    v = SV3(1.3, -0.4, 0.2)
    ω0 = SV3(0.0, 0.0, 0.0)
    τ = 0.05
    Δ_trans = pnl._rigid_back_displacement(te, origin, v, ω0, τ)
    @test Δ_trans ≈ -v * τ atol=1e-14

    # Pure rotation about z through a finite angle: the displaced point must lie
    # on the circle the trailing edge actually sweeps, i.e. the radius is
    # preserved exactly and the azimuth is rotated by exactly -θ.
    Ω = 2.0
    ω = SV3(0.0, 0.0, Ω)
    for τ_test in (0.01, 0.1, 0.35)          # θ = 0.02, 0.2, 0.7 rad
        θ = Ω * τ_test
        Δ = pnl._rigid_back_displacement(te, origin, SV3(0.0, 0.0, 0.0), ω, τ_test)
        p = te + Δ
        @test hypot(p[1], p[2]) ≈ hypot(te[1], te[2]) rtol=1e-13   # radius preserved
        @test atan(p[2], p[1]) ≈ -θ rtol=1e-12                     # exact arc angle
        @test p[3] ≈ te[3] atol=1e-14
    end

    # The legacy tangent construction instead lands at radius r*sqrt(1+θ²). Check
    # that the arc form removes exactly that error, and that the two agree to
    # O(θ²) so nothing changes in the small-angle limit.
    for τ_test in (0.001, 0.01, 0.35)
        θ = Ω * τ_test
        Δ_arc = pnl._rigid_back_displacement(te, origin, SV3(0.0, 0.0, 0.0), ω, τ_test)
        Δ_tan = -cross(ω, te - origin) * τ_test
        r = hypot(te[1], te[2])
        @test hypot((te + Δ_tan)[1], (te + Δ_tan)[2]) ≈ r * sqrt(1 + θ^2) rtol=1e-12
        @test norm(Δ_arc - Δ_tan) ≤ r * θ^2                  # second-order difference
    end

    # --- integration level: initialize_Das! ---------------------------------

    Uinf_zero = t -> SV3(0.0, 0.0, 0.0)
    axis = SV3(0.0, 0.0, 1.0)
    Ω_body = 3.0
    dt = 0.2                       # θ = eta*Ω*dt; large on purpose
    eta = 1.0

    function _das_endpoints(; arc::Bool)
        body = make_plate_vortex_body()
        body.Das[1] .= 0.0
        frames = pnl.ReferenceFrame(body;
            origin=SV3(0.0, 0.0, 0.0), v=SV3(0.0, 0.0, 0.0),
            ω_axis=axis, ω=Ω_body,
            R=FastMultipole.SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
            name="vehicle", child_index=Int[], dependent_index=[1])
        pnl.initialize_Das!((body,), frames, Uinf_zero, 0.0, dt;
            set_Das_eta_kinematic=eta, set_Das_kinematic_arc=arc)
        return body, frames
    end

    body_arc, _ = _das_endpoints(arc=true)
    body_tan, _ = _das_endpoints(arc=false)

    shedding = body_arc.shedding[1]
    θ_body = eta * Ω_body * dt
    radii_ok_arc = true
    radii_ok_tan = true
    for j in axes(shedding, 2)
        node_idx = body_arc.cells[shedding[3, j], shedding[1, j]]
        n = SV3(body_arc.nodes[1, node_idx], body_arc.nodes[2, node_idx], body_arc.nodes[3, node_idx])
        r = hypot(n[1], n[2])
        r == 0 && continue
        pa = n + SV3(body_arc.Das[1][1, j], body_arc.Das[1][2, j], body_arc.Das[1][3, j])
        pt = n + SV3(body_tan.Das[1][1, j], body_tan.Das[1][2, j], body_tan.Das[1][3, j])
        radii_ok_arc &= isapprox(hypot(pa[1], pa[2]), r; rtol=1e-12)
        radii_ok_tan &= isapprox(hypot(pt[1], pt[2]), r * sqrt(1 + θ_body^2); rtol=1e-12)
    end
    # The arc construction keeps every shed node on the swept circle; the legacy
    # tangent construction inflates its radius by sqrt(1+θ²).
    @test radii_ok_arc
    @test radii_ok_tan
    @test !isapprox(body_arc.Das[1], body_tan.Das[1]; rtol=1e-6)   # they really differ here

    # `arc=false` must reproduce the legacy path exactly.
    body_legacy = make_plate_vortex_body()
    body_legacy.Das[1] .= 0.0
    frames_legacy = pnl.ReferenceFrame(body_legacy;
        origin=SV3(0.0, 0.0, 0.0), v=SV3(0.0, 0.0, 0.0),
        ω_axis=axis, ω=Ω_body,
        R=FastMultipole.SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=[1])
    pnl.kinematic_velocity!((body_legacy,), frames_legacy)
    pnl._accumulate_Das!(body_legacy, eta * dt)
    @test body_legacy.Das[1] ≈ body_tan.Das[1] atol=1e-14

    # Small angle: arc and tangent must converge on each other.
    dt_small = 1e-3
    function _das_small(; arc::Bool)
        body = make_plate_vortex_body()
        body.Das[1] .= 0.0
        frames = pnl.ReferenceFrame(body;
            origin=SV3(0.0, 0.0, 0.0), v=SV3(0.0, 0.0, 0.0),
            ω_axis=axis, ω=Ω_body,
            R=FastMultipole.SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
            name="vehicle", child_index=Int[], dependent_index=[1])
        pnl.initialize_Das!((body,), frames, Uinf_zero, 0.0, dt_small;
            set_Das_eta_kinematic=eta, set_Das_kinematic_arc=arc)
        return body
    end
    d_arc = _das_small(arc=true).Das[1]
    d_tan = _das_small(arc=false).Das[1]
    # The two differ by the second-order arc-vs-chord term: |Δ_arc - Δ_tan| ≈ rθ²/2
    # against |Δ_tan| ≈ rθ, so the relative difference is θ/2.
    θ_small = eta * Ω_body * dt_small
    @test norm(d_arc - d_tan) / norm(d_tan) ≈ θ_small / 2 rtol=1e-3

    # Minimum-displacement floor still applies to the kinematic increment, with
    # the same semantics as the legacy path: nodes with zero kinematic velocity
    # (i.e. sitting on the rotation axis) have no defined direction and are left
    # untouched by both. Compare per-column magnitudes against legacy.
    function _das_floored(; arc::Bool)
        body = make_plate_vortex_body()
        body.Das[1] .= 0.0
        frames = pnl.ReferenceFrame(body;
            origin=SV3(0.0, 0.0, 0.0), v=SV3(0.0, 0.0, 0.0),
            ω_axis=axis, ω=Ω_body,
            R=FastMultipole.SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
            name="vehicle", child_index=Int[], dependent_index=[1])
        pnl.initialize_Das!((body,), frames, Uinf_zero, 0.0, dt_small;
            set_Das_eta_kinematic=eta, set_Das_min_kinematic_displacement=min_disp,
            set_Das_kinematic_arc=arc)
        return body.Das[1]
    end
    min_disp = 1.0e3      # far above any physical displacement here
    d_floor_arc = _das_floored(arc=true)
    d_floor_tan = _das_floored(arc=false)
    n_floored = 0
    for j in axes(d_floor_arc, 2)
        @test norm(d_floor_arc[:, j]) ≈ norm(d_floor_tan[:, j]) rtol=1e-12
        if norm(d_floor_tan[:, j]) > 0
            @test norm(d_floor_arc[:, j]) ≈ min_disp rtol=1e-12
            n_floored += 1
        end
    end
    @test n_floored > 0   # the floor actually fired somewhere
end

@testset "simulate! set_Das_refresh re-derives Das from the current state" begin
    # Time-growing freestream: with the default frozen Das the offset keeps its
    # t=0 magnitude for the whole run; with set_Das_refresh=true it is re-derived
    # each step from the current |Uinf|. BRAINSTORM 014.
    Uinf = t -> [1.0 + t, 0.0, 0.0]
    maneuver = (frames, systems, wakes, t) -> nothing
    t_range = [0.0, 0.1, 0.2]
    eta = 0.5
    dt = 0.1

    # finite attached wake: the semi-infinite kernel requires |Das| = 1 and
    # would reject an eta-scaled offset
    function _das_refresh_body()
        nodes = Float64[0 1 1 0; 0 0 1 1; 0 0 0 0]
        cells = Int[1 1; 2 3; 3 4]
        shedding = [pnl.calc_shedding_from_seed(nodes, cells, 1, 3)]
        body = pnl.RigidWakeBody{pnl.VortexRing}(nodes, cells, shedding;
            DBC=false, check_mesh=false, watertight=false,
            semiinfinite_wake=false)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        return body
    end

    function _run_das_refresh(; refresh::Bool)
        body = _das_refresh_body()
        body.Das[1] .= 0.0
        frames = pnl.ReferenceFrame(body)
        pnl.simulate!(body, nothing, frames, maneuver, Uinf, t_range;
            body_solvers=SimNoopSolver(),
            backend=pnl.DirectBackend(),
            path=nothing,
            set_Das_eta_freestream=eta,
            set_Das_refresh=refresh,
            grad_mu_options=(; basis=:tri),
        )
        return body.Das[1]
    end

    das_frozen = _run_das_refresh(refresh=false)
    das_refreshed = _run_das_refresh(refresh=true)

    # frozen: |Das| = eta*dt*|Uinf(0)| everywhere, along +x
    # refreshed: last refresh happens at the final step t = 0.2
    for j in axes(das_frozen, 2)
        @test das_frozen[:, j] ≈ [eta * dt * 1.0, 0.0, 0.0] atol=1e-13
        @test das_refreshed[:, j] ≈ [eta * dt * 1.2, 0.0, 0.0] atol=1e-13
    end

    # refresh with both etas NaN is a documented no-op (guard clause)
    body = _das_refresh_body()
    body.Das[1] .= 0.123
    frames = pnl.ReferenceFrame(body)
    pnl.simulate!(body, nothing, frames, maneuver, Uinf, t_range;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        path=nothing,
        set_Das_refresh=true,
        grad_mu_options=(; basis=:tri),
    )
    @test all(body.Das[1] .== 0.123)
end
