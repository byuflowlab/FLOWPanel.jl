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

mutable struct SimGeometryOrderSolver <: pnl.AbstractSolver
    transform_calls::Int
    geometry_current::Bool
end
SimGeometryOrderSolver() = SimGeometryOrderSolver(0, true)
pnl.solve!(::pnl.AbstractBody{<:Any, <:Any, <:Any, true},
    ::SimGeometryOrderSolver; kwargs...) = nothing
pnl.solve!(::pnl.AbstractBody{<:Any, <:Any, <:Any, false},
    ::SimGeometryOrderSolver; kwargs...) = nothing
function pnl.transform_solver_geometry!(solver::SimGeometryOrderSolver, body, R, t)
    reference = deepcopy(body)
    pnl.calc_normals!(reference)
    pnl.calc_controlpoints!(reference)
    solver.transform_calls += 1
    solver.geometry_current &= isapprox(body.normals, reference.normals;
        rtol=0, atol=1e-14)
    solver.geometry_current &= isapprox(body.controlpoints, reference.controlpoints;
        rtol=0, atol=1e-14)
    return nothing
end

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

@testset "simulate! refreshes geometry before persistent solvers" begin
    body = make_plate_vortex_body()
    axis = FastMultipole.SVector{3}(0.37, -0.61, 0.70)
    axis /= norm(axis)
    frames = pnl.ReferenceFrame(body;
        origin=FastMultipole.SVector{3}(0.0, 0.0, 0.0),
        v=FastMultipole.SVector{3}(0.1, -0.2, 0.05),
        ω_axis=axis, ω=4.0, dependent_index=[1])
    solver = SimGeometryOrderSolver()

    pnl.simulate!(body, nothing, frames, (args...) -> nothing,
        t -> FastMultipole.SVector{3}(0.0, 0.0, 0.0), [0.0, 0.05];
        body_solvers=solver, backend=pnl.DirectBackend(), path=nothing,
        grad_mu_options=(; basis=:tri))

    @test solver.transform_calls == 1
    @test solver.geometry_current
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

@testset "simulate! smooth surface-vorticity conversion smoke (BRAINSTORM 016)" begin
    # Phase 4 §2: the smooth strategy must survive the real time-marching
    # loop, not just static fixtures -- conversions fire when the panel buffer
    # overflows, and the external ledger totals must close.
    Uinf = t -> [1.0, 0.0, 0.0]
    maneuver = (frames, systems, wakes, t) -> nothing
    t_range = [0.0:0.05:0.3;]      # 6 steps: overflows the 2-row buffer

    body = make_plate_vortex_body()
    body.strength[1, 1] = 1.0
    body.strength[2, 1] = 2.0      # spanwise variation -> nonzero area deposit
    # Freestream-only convection keeps the sheet geometry regular: the noop
    # solver leaves O(1) body strengths whose induced velocity would fold the
    # coarse synthetic wake, which is not what this smoke test is probing.
    wake = pnl.PanelParticleWake(body; nwakerows=2, max_particles=10_000,
        freestream_convection=true,
        conversion=pnl.SurfaceVorticityConversion(0.05))
    frames = pnl.ReferenceFrame(body)
    pnl.simulate!(body, wake, frames, maneuver, Uinf, t_range;
        body_solvers=SimNoopSolver(),
        backend=pnl.DirectBackend(),
        path=nothing,
        grad_mu_options=(; basis=:tri),
    )

    @test wake.conversion_count[] >= 1
    diag = wake.conversion_diagnostics
    @test length(diag.records) == wake.conversion_count[]
    @test wake.pfield.np > 0
    @test all(isfinite, wake.pfield.particles[1:3, 1:wake.pfield.np])
    @test all(isfinite,
        wake.pfield.particles[FLOWVPM.GAMMA_INDEX, 1:wake.pfield.np])
    @test norm(diag.total_residual) < 1e-10
    @test diag.total_deposited ≈ diag.total_expected atol=1e-10
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

@testset "prescribed per-station Das lengths (chord-proportional, BRAINSTORM 018)" begin
    SV3 = FastMultipole.SVector{3, Float64}
    Uinf_zero = t -> SV3(0.0, 0.0, 0.0)
    axis = SV3(0.0, 0.0, 1.0)
    Ω_body = 3.0
    dt = 0.2
    eta = 1.0

    _rot_frames(body) = pnl.ReferenceFrame(body;
        origin=SV3(0.0, 0.0, 0.0), v=SV3(0.0, 0.0, 0.0),
        ω_axis=axis, ω=Ω_body,
        R=FastMultipole.SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=[1])

    # Reference: the tangent eta path, to borrow its per-station directions.
    body_eta = make_plate_vortex_body()
    body_eta.Das[1] .= 0.0
    pnl.initialize_Das!((body_eta,), _rot_frames(body_eta), Uinf_zero, 0.0, dt;
        set_Das_eta_kinematic=eta, set_Das_kinematic_arc=false)

    nstations = size(body_eta.Das[1], 2)
    lens = collect(range(0.05, 0.30; length=nstations))

    body = make_plate_vortex_body()
    body.Das[1] .= 0.0
    pnl.initialize_Das!((body,), _rot_frames(body), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),))

    # stations with zero TE kinematic velocity (on the rotation axis) carry no
    # direction and are skipped by both paths; assert on the moving ones
    moving = [j for j in 1:nstations if norm(body_eta.Das[1][:, j]) > 0]
    @test !isempty(moving)
    for j in moving
        dj = SV3(body.Das[1][1, j], body.Das[1][2, j], body.Das[1][3, j])
        de = SV3(body_eta.Das[1][1, j], body_eta.Das[1][2, j], body_eta.Das[1][3, j])
        # magnitude is exactly the prescribed length
        @test norm(dj) ≈ lens[j] rtol=1e-12
        # direction matches the tangent eta path exactly
        @test dj / norm(dj) ≈ de / norm(de) atol=1e-12
    end

    # min_displacement floor applies to prescribed lengths too
    body_floor = make_plate_vortex_body()
    body_floor.Das[1] .= 0.0
    floor_len = 0.2
    pnl.initialize_Das!((body_floor,), _rot_frames(body_floor), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),),
        set_Das_min_kinematic_displacement=floor_len)
    for j in moving
        dj = body_floor.Das[1][:, j]
        @test norm(dj) ≈ max(lens[j], floor_len) rtol=1e-12
    end

    # mutual exclusion with the eta path
    body_err = make_plate_vortex_body()
    @test_throws ErrorException pnl.initialize_Das!((body_err,), _rot_frames(body_err),
        Uinf_zero, 0.0, dt;
        set_Das_eta_kinematic=eta, set_Das_station_lengths=((lens,),))

    # station-count mismatch is loud
    body_err2 = make_plate_vortex_body()
    @test_throws ErrorException pnl.initialize_Das!((body_err2,), _rot_frames(body_err2),
        Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens[1:end-1],),))

    # frozen semantics: prescribed Das only rotates afterwards (magnitude kept)
    Rz = FastMultipole.SMatrix{3,3,Float64,9}(cos(0.3), sin(0.3), 0.0,
        -sin(0.3), cos(0.3), 0.0, 0.0, 0.0, 1.0)
    pnl.rotate_Das!(body, Rz)
    for j in moving
        @test norm(body.Das[1][:, j]) ≈ lens[j] rtol=1e-12
    end
end
@testset "uniform-d_front Das law identity (BRAINSTORM 018 phase_13 s4b)" begin
    SV3 = FastMultipole.SVector{3, Float64}
    Uinf_zero = t -> SV3(0.0, 0.0, 0.0)
    axis = SV3(0.0, 0.0, 1.0)
    Ω_body = 3.0
    dt = 0.2

    _rot_frames(body) = pnl.ReferenceFrame(body;
        origin=SV3(0.0, 0.0, 0.0), v=SV3(0.0, 0.0, 0.0),
        ω_axis=axis, ω=Ω_body,
        R=FastMultipole.SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=[1])

    # eta=1 reference: |Das_eta_j| = 1*dt*|Vte_j| = one step of TE travel,
    # exactly the travel term of the s4b construction.
    body_eta = make_plate_vortex_body()
    body_eta.Das[1] .= 0.0
    pnl.initialize_Das!((body_eta,), _rot_frames(body_eta), Uinf_zero, 0.0, dt;
        set_Das_eta_kinematic=1.0, set_Das_kinematic_arc=false)
    nstations = size(body_eta.Das[1], 2)
    travel = [norm(body_eta.Das[1][:, j]) for j in 1:nstations]
    moving = [j for j in 1:nstations if travel[j] > 0]
    @test !isempty(moving)

    # d_front = |Das| + (N-1)*travel must equal D*sigma at every moving station
    N = 3
    sigma = maximum(travel)          # keeps every lens strictly positive
    D = 3.4
    lens = D * sigma .- (N - 1) .* travel
    @test minimum(lens[moving]) > 0
    body = make_plate_vortex_body()
    body.Das[1] .= 0.0
    pnl.initialize_Das!((body,), _rot_frames(body), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),))
    for j in moving
        d_front = norm(body.Das[1][:, j]) + (N - 1) * travel[j]
        @test d_front ≈ D * sigma rtol=1e-12
    end

    # N=1 degenerate case: no travel term, |Das| = D*sigma span-uniform
    body1 = make_plate_vortex_body()
    body1.Das[1] .= 0.0
    pnl.initialize_Das!((body1,), _rot_frames(body1), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((fill(D * sigma, nstations),),))
    for j in moving
        @test norm(body1.Das[1][:, j]) ≈ D * sigma rtol=1e-12
    end

    # Chord–sigma co-scaling law (018 Phase 16): |Das|_j = lambda*sigma_j with
    # sigma_j prescribed per station — Das/sigma must come out span-uniform at
    # exactly lambda on the moving stations.
    lambda = 3.4
    sigmas = collect(range(0.02, 0.08; length=nstations))
    body_cs = make_plate_vortex_body()
    body_cs.Das[1] .= 0.0
    pnl.initialize_Das!((body_cs,), _rot_frames(body_cs), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lambda .* sigmas,),))
    for j in moving
        @test norm(body_cs.Das[1][:, j]) / sigmas[j] ≈ lambda rtol=1e-12
    end
end

@testset "endpoint-on-arc Das (018 Phase 16 F1b Route B)" begin
    SV3 = FastMultipole.SVector{3, Float64}
    Uinf_zero = t -> SV3(0.0, 0.0, 0.0)
    axis = SV3(0.0, 0.0, 1.0)
    Ω_body = 3.0
    dt = 0.2

    _rot_frames(body) = pnl.ReferenceFrame(body;
        origin=SV3(0.0, 0.0, 0.0), v=SV3(0.0, 0.0, 0.0),
        ω_axis=axis, ω=Ω_body,
        R=FastMultipole.SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=[1])

    _te_positions(body) = begin
        shed = body.shedding[1]
        n = size(shed, 2) + 1
        tes = Vector{SV3}(undef, n)
        for j in 1:size(shed, 2)
            idx = body.cells[shed[3, j], shed[1, j]]
            tes[j] = SV3(body.nodes[1, idx], body.nodes[2, idx], body.nodes[3, idx])
        end
        idx = body.cells[shed[2, end], shed[1, end]]
        tes[end] = SV3(body.nodes[1, idx], body.nodes[2, idx], body.nodes[3, idx])
        tes
    end

    body_probe = make_plate_vortex_body()
    nstations = size(body_probe.Das[1], 2)
    tes = _te_positions(body_probe)
    radii = [sqrt(te[1]^2 + te[2]^2) for te in tes]
    moving = [j for j in 1:nstations if radii[j] > 0]
    @test !isempty(moving)
    lens = collect(range(0.05, 0.30; length=nstations))
    zero_drift = zeros(3, nstations)

    # (1) kinematic-only reduction (u = 0): endpoint lies ON the swept circle
    # (radius preserved to machine precision — with u = 0 the path speed is
    # constant so the arc-length integration is exact), chord = 2r·sin(L/2r),
    # and the result differs from the tangent path at finite theta.
    body = make_plate_vortex_body()
    body.Das[1] .= 0.0
    pnl.initialize_Das!((body,), _rot_frames(body), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),),
        set_Das_station_drifts=((zero_drift,),))
    body_tan = make_plate_vortex_body()
    body_tan.Das[1] .= 0.0
    pnl.initialize_Das!((body_tan,), _rot_frames(body_tan), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),))
    for j in moving
        r = radii[j]
        endpoint = tes[j] + SV3(body.Das[1][1, j], body.Das[1][2, j], body.Das[1][3, j])
        @test sqrt(endpoint[1]^2 + endpoint[2]^2) ≈ r rtol=1e-10
        @test endpoint[3] ≈ tes[j][3] atol=1e-12
        chord = norm(body.Das[1][:, j])
        @test chord ≈ 2r * sin(lens[j] / (2r)) rtol=1e-10
        # tangent endpoint sits OFF the circle at r*sqrt(1+theta^2)
        theta = lens[j] / r
        end_tan = tes[j] + SV3(body_tan.Das[1][1, j], body_tan.Das[1][2, j],
            body_tan.Das[1][3, j])
        @test sqrt(end_tan[1]^2 + end_tan[2]^2) ≈ r * sqrt(1 + theta^2) rtol=1e-10
    end

    # (2) axial drift: exact helix — u ⟂ the rotation plane keeps path speed
    # constant, so tau = L/sqrt((Ωr)² + w²) exactly; z advances by w·tau and
    # the in-plane radius is preserved.
    w = 0.8
    drift_ax = vcat(zeros(2, nstations), fill(w, 1, nstations))
    body_hx = make_plate_vortex_body()
    body_hx.Das[1] .= 0.0
    pnl.initialize_Das!((body_hx,), _rot_frames(body_hx), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),),
        set_Das_station_drifts=((drift_ax,),))
    for j in moving
        r = radii[j]
        τ = lens[j] / sqrt((Ω_body * r)^2 + w^2)
        endpoint = tes[j] + SV3(body_hx.Das[1][1, j], body_hx.Das[1][2, j],
            body_hx.Das[1][3, j])
        @test sqrt(endpoint[1]^2 + endpoint[2]^2) ≈ r rtol=1e-10
        @test endpoint[3] - tes[j][3] ≈ w * τ rtol=1e-10
    end

    # (3) chord never exceeds the prescribed arc length (generic drift)
    drift_gen = 0.3 .* randn(3, nstations)
    body_gen = make_plate_vortex_body()
    body_gen.Das[1] .= 0.0
    pnl.initialize_Das!((body_gen,), _rot_frames(body_gen), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),),
        set_Das_station_drifts=((drift_gen,),))
    for j in moving
        # chord ≤ traversed path length; the explicit arc-length stepper has
        # O(δs) first-order error when the speed varies along the path, so
        # allow a few percent
        @test norm(body_gen.Das[1][:, j]) <= lens[j] * 1.05
        @test all(isfinite, body_gen.Das[1][:, j])
    end

    # (4) min_displacement floor applies to the ARC length
    body_floor = make_plate_vortex_body()
    body_floor.Das[1] .= 0.0
    floor_len = 0.2
    pnl.initialize_Das!((body_floor,), _rot_frames(body_floor), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),),
        set_Das_station_drifts=((zero_drift,),),
        set_Das_min_kinematic_displacement=floor_len)
    for j in moving
        r = radii[j]
        L = max(lens[j], floor_len)
        @test norm(body_floor.Das[1][:, j]) ≈ 2r * sin(L / (2r)) rtol=1e-10
    end

    # (5) drift exactly cancelling the kinematic TE velocity at zero lag:
    # speed0 = 0 => station skipped (left untouched), matching legacy handling
    j0 = moving[end]
    drift_cancel = copy(zero_drift)
    vkin = cross(SV3(0.0, 0.0, Ω_body), tes[j0])   # blade-point velocity ω×r
    drift_cancel[:, j0] .= vkin                    # u = -Vte = +ω×r
    body_stag = make_plate_vortex_body()
    body_stag.Das[1] .= 0.0
    pnl.initialize_Das!((body_stag,), _rot_frames(body_stag), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),),
        set_Das_station_drifts=((drift_cancel,),))
    @test norm(body_stag.Das[1][:, j0]) == 0.0
    for j in moving
        j == j0 && continue
        @test norm(body_stag.Das[1][:, j]) > 0
    end

    # (6) loud errors: drifts without lengths; drift shape mismatch
    body_err = make_plate_vortex_body()
    @test_throws ErrorException pnl.initialize_Das!((body_err,),
        _rot_frames(body_err), Uinf_zero, 0.0, dt;
        set_Das_station_drifts=((zero_drift,),))
    body_err2 = make_plate_vortex_body()
    @test_throws ErrorException pnl.initialize_Das!((body_err2,),
        _rot_frames(body_err2), Uinf_zero, 0.0, dt;
        set_Das_station_lengths=((lens,),),
        set_Das_station_drifts=((zero_drift[:, 1:end-1],),))

    # (7) frozen semantics: rotates afterwards with magnitude kept
    Rz = FastMultipole.SMatrix{3,3,Float64,9}(cos(0.3), sin(0.3), 0.0,
        -sin(0.3), cos(0.3), 0.0, 0.0, 0.0, 1.0)
    mags = [norm(body.Das[1][:, j]) for j in moving]
    pnl.rotate_Das!(body, Rz)
    for (i, j) in enumerate(moving)
        @test norm(body.Das[1][:, j]) ≈ mags[i] rtol=1e-12
    end
end

@testset "WakeHealthMonitor max_dtZ column (dtz kwarg)" begin
    using StaticArrays: SVector

    # Tiny benign particle field: default rVPM formulation f=0, g=1/5, so
    # Z = (1/5) * Gamma'*(J-projection)*Gamma / |Gamma|^2. With Gamma along x
    # and only J[1] = a nonzero, Z = a/5 under both transposed and classic
    # schemes, so max_dtZ = dt*a/5 analytically.
    function make_dtz_pfield()
        pf = FLOWVPM.ParticleField(10, Float64;
            fmm=FLOWVPM.FMM(autotune_reg_error=false))
        FLOWVPM.add_particle(pf, SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0), 0.05)
        FLOWVPM.add_particle(pf, SVector(1.0, 0.0, 0.0),
            SVector(0.0, 0.5, 0.0), 0.04)
        FLOWVPM.add_particle(pf, SVector(2.0, 0.0, 0.0),
            SVector(0.0, 0.0, 0.0), 0.03)   # zero-Gamma: Z must be 0
        FLOWVPM.get_J(pf, 1)[1] = 5.0        # strongest contraction
        FLOWVPM.get_J(pf, 2) .= 0.0
        FLOWVPM.get_J(pf, 2)[5] = -2.0       # expansive for particle 2
        FLOWVPM.get_U(pf, 1) .= (1.0, 2.0, 2.0)
        return pf
    end

    dt = 1e-3
    legacy_header = "step,time,n_particles,max_u,min_sigma,min_sigma_ratio," *
                    "max_gamma_over_sigma2,wall_s"

    pf = make_dtz_pfield()
    stats = pnl._wake_health_stats(pf, NaN, dt)
    @test length(stats) == 6
    @test stats[1] == 3.0
    @test stats[2] ≈ 3.0                      # max|u|
    @test stats[6] ≈ dt * 5.0 / 5 rtol=1e-12  # analytic dt*Z of particle 1
    @test isfinite(stats[6]) && abs(stats[6]) < 1

    # (a) dtz=false: CSV bit-identical structure to the legacy monitor
    mktempdir() do dir
        m = pnl.WakeHealthMonitor()           # default dtz=false
        @test m.dtz == false
        push!(m.stats, stats)
        m.wall_s = 0.5
        pnl.write_monitor_csv!(m, dir, "t", 1, pnl.MonitorContext(), (), 0, dt)
        lines = readlines(joinpath(dir, "t_monitor01_wake_health_system1.csv"))
        @test lines[1] == legacy_header
        @test count(==(','), lines[2]) == 7   # 8 columns, unchanged
    end

    # (b) dtz=true: column present, finite, benign magnitude
    mktempdir() do dir
        m = pnl.WakeHealthMonitor(dtz=true)
        push!(m.stats, stats)
        m.wall_s = 0.5
        pnl.write_monitor_csv!(m, dir, "t", 1, pnl.MonitorContext(), (), 0, dt)
        lines = readlines(joinpath(dir, "t_monitor01_wake_health_system1.csv"))
        @test lines[1] == legacy_header * ",max_dtZ"
        cols = split(lines[2], ',')
        @test length(cols) == 9
        val = parse(Float64, cols[end])
        @test isfinite(val) && abs(val) < 1
        @test val ≈ stats[6]
    end

    # (c) cross-check against FLOWVPM's own sigma update: the Euler step
    # applies sigma -= dt*sigma*Z with the same J, so the largest fractional
    # contraction -Delta sigma/sigma must equal the reported max_dtZ.
    pf2 = make_dtz_pfield()
    stats2 = pnl._wake_health_stats(pf2, NaN, dt)
    sig_pre = [FLOWVPM.get_sigma(pf2, i)[1] for i in 1:pf2.np]
    FLOWVPM._euler(pf2, dt)
    contraction = maximum((sig_pre[i] - FLOWVPM.get_sigma(pf2, i)[1]) / sig_pre[i]
                          for i in 1:pf2.np)
    @test contraction ≈ stats2[6] rtol=1e-10
end

@testset "WakeHealthMonitor min_sr attribution columns (attribution kwarg)" begin
    using StaticArrays: SVector

    pf = FLOWVPM.ParticleField(10, Float64;
        fmm=FLOWVPM.FMM(autotune_reg_error=false))
    FLOWVPM.add_particle(pf, SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), 0.05)
    FLOWVPM.add_particle(pf, SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), 0.04)
    FLOWVPM.add_particle(pf, SVector(2.0, 0.5, -0.5), SVector(0.0, 0.0, 1.0), 0.03)

    m = pnl.WakeHealthMonitor(attribution=true, sigma_ref=0.1)
    @test m.attribution == true

    # type-7 p1 of sorted [0.03, 0.04, 0.05]: h = 1.02 -> 0.03 + 0.02*0.01 =
    # 0.0302; ratio = 0.0302/0.1. Argmin is particle 3's position.
    attr = pnl._wake_health_attribution!(m.sigma_buf, pf, 0.1)
    @test attr[1] ≈ 0.302 rtol=1e-12
    @test attr[2] == 2.0 && attr[3] == 0.5 && attr[4] == -0.5

    # No sigma reference -> ratio NaN, position still reported
    attr_nan = pnl._wake_health_attribution!(m.sigma_buf, pf, NaN)
    @test isnan(attr_nan[1]) && attr_nan[2] == 2.0

    # attribution=true: four columns appended after the legacy set (dtz off)
    legacy_header = "step,time,n_particles,max_u,min_sigma,min_sigma_ratio," *
                    "max_gamma_over_sigma2,wall_s"
    mktempdir() do dir
        push!(m.stats, pnl._wake_health_stats(pf, 0.1, 1e-3))
        push!(m.attr, attr)
        m.wall_s = 0.5
        pnl.write_monitor_csv!(m, dir, "t", 1, pnl.MonitorContext(), (), 0, 1e-3)
        lines = readlines(joinpath(dir, "t_monitor01_wake_health_system1.csv"))
        @test lines[1] == legacy_header *
                          ",p1_sigma_ratio,argmin_x,argmin_y,argmin_z"
        cols = split(lines[2], ',')
        @test length(cols) == 12
        @test parse(Float64, cols[9]) ≈ 0.302 rtol=1e-12
        @test parse(Float64, cols[10]) == 2.0
    end
    # (default attribution=false bit-identity is asserted by the dtz testset
    # above: legacy header, 8 columns.)
end

@testset "WakeInventoryMonitor banded inventory (BRAINSTORM 018 phase 15)" begin
    using StaticArrays: SVector

    R = 2.0
    pf = FLOWVPM.ParticleField(10, Float64;
        fmm=FLOWVPM.FMM(autotune_reg_error=false))
    # band 1 (z/R 0-0.5) inboard: xi = 0.25, r = 0
    FLOWVPM.add_particle(pf, SVector(0.5, 0.0, 0.0), SVector(1.0, 0.0, 0.0), 0.05)
    # internal boundary xi = 0.5 exactly -> assigned to the LOWER band (band 1), r = 0
    FLOWVPM.add_particle(pf, SVector(1.0, 0.0, 0.0), SVector(0.5, 0.0, 0.0), 0.07)
    # band 1 outboard: xi = 0.25, r = 1.2 >= r_split
    FLOWVPM.add_particle(pf, SVector(0.5, 2.4, 0.0), SVector(-1.0, 0.0, 0.0), 0.04)
    # band 4 (z/R 1-2): xi = 1.5; antiparallel pair -> sum|G| = 4, vector sum = 0
    FLOWVPM.add_particle(pf, SVector(3.0, 0.0, 0.0), SVector(0.0, 2.0, 0.0), 0.03)
    FLOWVPM.add_particle(pf, SVector(3.0, 0.5, 0.0), SVector(0.0, -2.0, 0.0), 0.06)
    # outside: xi = 4 > outer edge 3.5
    FLOWVPM.add_particle(pf, SVector(8.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0), 0.02)
    # static particle in band 1: must be SKIPPED
    FLOWVPM.add_particle(pf, SVector(0.5, 0.0, 0.1), SVector(9.0, 0.0, 0.0), 0.01;
                         static=true)

    m = pnl.WakeInventoryMonitor(R)
    # cells: rin, rout, z0p5_1, z1_2, z2_3, z3_3p5, outside
    @test length(m.cell_names) == 7
    @test m.cell_names[1] == "z0p0_0p5_rin"
    @test m.cell_names[2] == "z0p0_0p5_rout"
    @test m.cell_names[end] == "outside"

    vals = fill(NaN, 9*7)
    bufs = [Float64[] for _ in 1:7]
    pnl._wake_inventory_stats!(vals, bufs, pf, m)
    n_of(c) = vals[(c - 1)*9 + 1]

    @test n_of(1) == 2                    # rin: xi=0.25 + boundary xi=0.5; static skipped
    @test n_of(2) == 1                    # rout
    @test n_of(3) == 0 && n_of(5) == 0 && n_of(6) == 0
    @test n_of(4) == 2
    @test n_of(7) == 1                    # outside
    @test sum(n_of(c) for c in 1:7) == 6  # partition of the non-static field

    # cancellation-aware vs magnitude inventory (band 4)
    b4 = (4 - 1)*9
    @test vals[b4 + 2] ≈ 4.0              # sum|Gamma|
    @test abs(vals[b4 + 4]) < 1e-14       # vector Gamma_y cancels
    # sigma stats over [0.03, 0.06]: mean = p50 = 0.045; type-7 p5/p95
    @test vals[b4 + 6] ≈ 0.045 && vals[b4 + 8] ≈ 0.045
    @test vals[b4 + 7] ≈ 0.0315 rtol=1e-12
    @test vals[b4 + 9] ≈ 0.0585 rtol=1e-12
    # empty cell -> NaN sigma stats, zero sums
    @test isnan(vals[(3 - 1)*9 + 6]) && vals[(3 - 1)*9 + 2] == 0.0

    # CSV: 2 + 9*7 columns, counts written as integers
    mktempdir() do dir
        push!(m.stats, vals)
        pnl.write_monitor_csv!(m, dir, "t", 1, pnl.MonitorContext(), (), 0, 1e-3)
        lines = readlines(joinpath(dir, "t_monitor01_wake_inventory_system1.csv"))
        hdr = split(lines[1], ',')
        @test length(hdr) == 2 + 9*7
        @test hdr[3] == "z0p0_0p5_rin_n"
        @test hdr[end] == "outside_sigma_p95"
        cols = split(lines[2], ',')
        @test length(cols) == 2 + 9*7
        @test cols[3] == "2"              # Int-formatted count
    end
end

@testset "PanelParticleWake expint kwarg (BRAINSTORM 020 Stage B)" begin
    # Off (default): stock forward Euler — the step! branch selects _euler and
    # the recorded integration is unchanged from prior behavior.
    body = make_plate_vortex_body()
    wake_off = pnl.PanelParticleWake(body; nwakerows=2, particle_kerneloffset=0.1)
    @test wake_off.pfield.integration === FLOWVPM.euler
    @test wake_off.pfield_optargs.integration === FLOWVPM.euler

    # On: integration is the exponential integrator and is echoed into the
    # reproduction metadata (pfield.integration is the single source of truth
    # for the step! branch).
    body2 = make_plate_vortex_body()
    wake_on = pnl.PanelParticleWake(body2; nwakerows=2, particle_kerneloffset=0.1,
        expint=true)
    @test wake_on.pfield.integration === FLOWVPM.euler_exp
    @test wake_on.pfield_optargs.integration === FLOWVPM.euler_exp
end

@testset "simulate! convert-at-shed N=0 smoke (BRAINSTORM 024)" begin
    # nwakerows=0 must survive the real time-marching loop with both
    # conversion strategies: every step sheds AND converts (no free sheet
    # enters any solve), particles stay finite, and the smooth ledger closes.
    Uinf = t -> [1.0, 0.0, 0.0]
    maneuver = (frames, systems, wakes, t) -> nothing
    t_range = [0.0:0.05:0.3;]

    for conversion_kwargs in (
            (;),                                       # LegacyEdgeJumpConversion
            (; conversion=pnl.SurfaceVorticityConversion(0.05;
                attribution=:downstream)),
        )
        body = make_plate_vortex_body()
        body.strength[1, 1] = 1.0
        body.strength[2, 1] = 2.0
        wake = pnl.PanelParticleWake(body; nwakerows=0, max_particles=10_000,
            freestream_convection=true, conversion_kwargs...)
        frames = pnl.ReferenceFrame(body)
        pnl.simulate!(body, wake, frames, maneuver, Uinf, t_range;
            body_solvers=SimNoopSolver(),
            backend=pnl.DirectBackend(),
            path=nothing,
            grad_mu_options=(; basis=:tri),
        )

        pw = wake.panel_wake
        @test pw.convert_at_shed
        @test pw.nwakes[] == 0            # no free sheet survives any step
        @test pw.overflowed[]
        np = wake.pfield.np
        @test np > 0
        @test all(isfinite, wake.pfield.particles[1:3, 1:np])
        @test all(isfinite, wake.pfield.particles[FLOWVPM.GAMMA_INDEX, 1:np])
        @test all(>(0), wake.pfield.particles[FLOWVPM.SIGMA_INDEX, 1:np])
        if wake.conversion isa pnl.SurfaceVorticityConversion
            # one conversion per shed, from the very first step
            @test wake.conversion_count[] == length(t_range) - 1
            diag = wake.conversion_diagnostics
            @test norm(diag.total_residual) < 1e-10
            @test diag.total_deposited ≈ diag.total_expected atol=1e-10
        end
    end
end
