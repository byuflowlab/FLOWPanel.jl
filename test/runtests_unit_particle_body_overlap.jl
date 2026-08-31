using Test
using Logging
import FLOWPanel as pnl
import FastMultipole
import FLOWVPM
using StaticArrays: SVector

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

const _OVERLAP_GAMMA = SVector(1.0, 0.0, 0.0)

function _overlap_wake(body; pfield=nothing, arraytype=Matrix)
    return pnl.PanelParticleWake(body; nwakerows=2, max_particles=32,
        particle_core_size=0.1, pfield, arraytype)
end

function _add_overlap_particles!(pfield)
    FLOWVPM.add_particle(pfield, SVector(0.25, 0.25, 0.20),
        _OVERLAP_GAMMA, 0.25)
    FLOWVPM.add_particle(pfield, SVector(4.0, 4.0, 4.0),
        _OVERLAP_GAMMA, 0.10)
    return pfield
end

@testset "particle/body exact point-triangle distance" begin
    a = SVector(0.0, 0.0, 0.0)
    b = SVector(1.0, 0.0, 0.0)
    c = SVector(0.0, 1.0, 0.0)

    @test sqrt(pnl._point_triangle_distance2(
        SVector(0.25, 0.25, 0.20), a, b, c)) ≈ 0.20
    @test sqrt(pnl._point_triangle_distance2(
        SVector(0.50, -0.30, 0.0), a, b, c)) ≈ 0.30
    @test sqrt(pnl._point_triangle_distance2(
        SVector(-0.30, -0.40, 0.0), a, b, c)) ≈ 0.50
end

@testset "particle/body overlap report and shared-field deduplication" begin
    body1 = make_plate_vortex_body()
    body2 = make_plate_vortex_body()
    body2.nodes[1, :] .+= 10.0
    pnl.calc_normals!(body2)
    pnl.calc_controlpoints!(body2)

    wake1 = _overlap_wake(body1)
    _add_overlap_particles!(wake1.pfield)
    wake2 = _overlap_wake(body2; pfield=wake1.pfield)

    report = pnl.particle_body_overlap((body1, body2), (wake1, wake2);
        core_ratio=1.0)
    @test report.checked_fields == 1
    @test report.checked_particles == 2
    @test report.overlap_count == 1
    @test report.min_distance ≈ 0.20
    @test report.min_distance_over_sigma ≈ 0.80
    @test report.body_index == 1
    @test report.wake_index == 1
    @test report.particle_index == 1
    @test report.particle_position == SVector(0.25, 0.25, 0.20)
    @test report.particle_sigma == 0.25
    @test report.threshold_ratio == 1.0
    @test report.triangle_evaluations <=
        report.checked_particles * (body1.ncells + body2.ncells)

    empty_report = pnl.particle_body_overlap(body1, nothing)
    @test empty_report.checked_fields == 0
    @test empty_report.checked_particles == 0
    @test empty_report.overlap_count == 0
    @test isinf(empty_report.min_distance)

    invalid_wake = _overlap_wake(body1)
    _add_overlap_particles!(invalid_wake.pfield)
    invalid_wake.pfield.particles[1, 1] = NaN
    @test_throws ArgumentError pnl.particle_body_overlap(body1, invalid_wake)
    invalid_wake.pfield.particles[1, 1] = 0.25
    invalid_wake.pfield.particles[FLOWVPM.SIGMA_INDEX, 1] = NaN
    @test_throws ArgumentError pnl.particle_body_overlap(body1, invalid_wake)
end

@testset "particle/body overlap policy" begin
    @test_throws ArgumentError pnl.ParticleBodyOverlapPolicy(action=:ignore)
    @test_throws ArgumentError pnl.ParticleBodyOverlapPolicy(core_ratio=0)
    @test_throws ArgumentError pnl.ParticleBodyOverlapPolicy(every=0)

    body = make_plate_vortex_body()
    wake = _overlap_wake(body)
    _add_overlap_particles!(wake.pfield)

    @test pnl.check_particle_body_overlap!(
        pnl.ParticleBodyOverlapPolicy(), body, wake, 0) === nothing
    warn_policy = pnl.ParticleBodyOverlapPolicy(action=:warn, every=2)
    @test pnl.check_particle_body_overlap!(warn_policy, body, wake, 1) === nothing
    @test_logs (:warn, r"particle/core overlaps a panel body") begin
        report = pnl.check_particle_body_overlap!(warn_policy, body, wake, 2)
        @test report.overlap_count == 1
    end

    error_policy = pnl.ParticleBodyOverlapPolicy(action=:error)
    err = try
        pnl.check_particle_body_overlap!(error_policy, body, wake, 0)
        nothing
    catch caught
        caught
    end
    @test err isa pnl.ParticleBodyOverlapError
    @test err.report.overlap_count == 1
    @test occursin("distance/sigma=0.8", sprint(showerror, err))
end

@testset "simulate! invokes overlap gate before wake influence" begin
    body = make_plate_vortex_body()
    wake = _overlap_wake(body)
    _add_overlap_particles!(wake.pfield)
    frames = pnl.ReferenceFrame(body)
    solver = pnl.Backslash(body)
    maneuver = (frames, systems, wakes, t) -> nothing
    Uinf = t -> SVector(0.0, 0.0, 0.0)

    @test_throws pnl.ParticleBodyOverlapError pnl.simulate!(
        body, wake, frames, maneuver, Uinf, [0.0, 0.01];
        body_solvers=solver,
        backend=pnl.DirectBackend(),
        path=nothing,
        grad_mu_options=(; basis=:tri),
        particle_body_overlap_policy=pnl.ParticleBodyOverlapPolicy(
            action=:error))

    # The production telemetry seam observes the already-computed report and
    # formulation state after the step; it must not trigger a second scan.
    body2 = make_plate_vortex_body()
    wake2 = _overlap_wake(body2)
    frames2 = pnl.ReferenceFrame(body2)
    solver2 = pnl.Backslash(body2)
    events = NamedTuple[]
    pnl.simulate!(body2, wake2, frames2, maneuver, Uinf, [0.0, 0.01];
        body_solvers=solver2,
        backend=pnl.DirectBackend(),
        path=nothing,
        grad_mu_options=(; basis=:tri),
        particle_body_overlap_policy=pnl.ParticleBodyOverlapPolicy(
            action=:warn),
        step_telemetry_callback=event -> push!(events, event))
    @test length(events) == 2
    @test events[1].i_step == 0
    @test events[1].overlap_report isa pnl.ParticleBodyOverlapReport
    @test events[1].overlap_report.checked_particles == 0
    @test events[1].formulation isa pnl.VelocityThroughSources
end

# Opt-in CUDA gate. The host report and the device field's refreshed host
# mirror must agree exactly on attribution and to roundoff on distances.
if parse(Bool, get(ENV, "FLOWPANEL_TEST_PARTICLE_BODY_OVERLAP_CUDA", "false"))
    @testset "particle/body CUDA host-mirror equivalence" begin
        @test FastMultipole.load_cuda_radix_lifecycle!()
        CUDAmod = getglobal(FastMultipole, :CUDA)
        @test Base.invokelatest(CUDAmod.functional)

        body = make_plate_vortex_body()
        host_wake = _overlap_wake(body)
        _add_overlap_particles!(host_wake.pfield)
        host_report = pnl.particle_body_overlap(body, host_wake)

        device_wake = _overlap_wake(body; arraytype=CUDAmod.CuArray)
        mirror = pnl._gpu_pfield_mirror(device_wake.pfield)
        _add_overlap_particles!(mirror)
        pnl._gpu_sync_device_from_mirror!(device_wake.pfield, mirror)
        device_report = pnl.particle_body_overlap(body, device_wake)

        @test device_report.checked_fields == host_report.checked_fields
        @test device_report.checked_particles == host_report.checked_particles
        @test device_report.overlap_count == host_report.overlap_count
        @test device_report.min_distance ≈ host_report.min_distance rtol=0 atol=1e-12
        @test device_report.min_distance_over_sigma ≈
            host_report.min_distance_over_sigma rtol=0 atol=1e-12
        @test device_report.body_index == host_report.body_index
        @test device_report.wake_index == host_report.wake_index
        @test device_report.particle_index == host_report.particle_index
        @test device_report.particle_position ≈ host_report.particle_position
        @test device_report.particle_sigma ≈ host_report.particle_sigma
    end
end
