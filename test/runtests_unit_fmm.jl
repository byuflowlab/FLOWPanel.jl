using Test
import FLOWPanel as pnl
import FastMultipole

if !isdefined(@__MODULE__, :make_octa_source_body)
    include("test_helpers.jl")
end

@testset verbose=true "FMM Backends" begin
    @testset "backend construction" begin
        @test pnl.DirectBackend() isa pnl.DirectBackend
        fmm = pnl.FastMultipoleBackend()
        @test fmm.expansion_order == 10
        @test fmm.multipole_acceptance == 0.4
        @test fmm.leaf_size == 20

        custom = pnl.FastMultipoleBackend(expansion_order=10)
        @test custom.expansion_order == 10
    end

    @testset "Direct vs FMM NonLiftingBody" begin
        source_body = make_octa_source_body()
        target_body = translated_nonlifting_target([3.0, -1.0, 2.5])
        direct = pnl.DirectBackend()
        fmm = pnl.FastMultipoleBackend(expansion_order=10)

        vel_direct, phi_direct = evaluate_velocity_and_potential(deepcopy(target_body), source_body, direct)
        vel_fmm, phi_fmm = evaluate_velocity_and_potential(deepcopy(target_body), source_body, fmm)

        @test all(isfinite, vel_direct)
        @test all(isfinite, vel_fmm)
        @test all(isfinite, phi_direct)
        @test all(isfinite, phi_fmm)
        @test maximum(abs.(vel_direct .- vel_fmm)) <= 1e-3
        @test maximum(abs.(phi_direct .- phi_fmm)) <= 1e-3
    end

    @testset "relative-error metadata" begin
        source_body = make_octa_source_body()
        target_body = translated_nonlifting_target([3.0, -1.0, 2.5])
        backend = pnl.FastMultipoleBackend(expansion_order=5)

        @test_nowarn pnl.influence!((target_body,), (source_body,), backend;
            scalar_potential=true,
            velocity=true,
            error_tolerance=FastMultipole.PowerRelativeGradient(1e-2))
    end

    @testset "Direct vs FMM RigidWakeBody" begin
        source_body = make_plate_vortex_body()
        target_body = translated_rigid_target([2.0, -1.5, 1.0])
        direct = pnl.DirectBackend()
        fmm = pnl.FastMultipoleBackend(expansion_order=10)

        vel_direct, _ = evaluate_velocity_and_potential(deepcopy(target_body), source_body, direct)
        vel_fmm, _ = evaluate_velocity_and_potential(deepcopy(target_body), source_body, fmm)

        @test all(isfinite, vel_direct)
        @test all(isfinite, vel_fmm)
        @test maximum(abs.(vel_direct .- vel_fmm)) <= 1e-3
    end

    @testset "self panel conditioning splits body offsets" begin
        body1 = make_plate_vortex_body()
        body2 = translated_rigid_target([1.7, -0.4, 0.6])
        for body in (body1, body2)
            body.kerneloffset_panel = 1e-8
            body.kerneloffset_targets = 0.15
            body.kerneloffset = body.kerneloffset_targets
        end

        direct = pnl.DirectBackend()
        combined1 = deepcopy(body1)
        combined2 = deepcopy(body2)
        for body in (combined1, combined2)
            body.velocity .= 0
            body.kerneloffset = body.kerneloffset_targets
        end
        pnl.influence!((combined1, combined2), (combined1, combined2), direct;
            velocity=true,
            direct_conditioning=pnl._self_panel_kerneloffset_conditioning())

        explicit1 = deepcopy(body1)
        explicit2 = deepcopy(body2)
        explicit1.velocity .= 0
        explicit2.velocity .= 0
        source1 = deepcopy(body1)
        source2 = deepcopy(body2)

        source1.kerneloffset = source1.kerneloffset_panel
        pnl.influence!((explicit1,), (source1,), direct; velocity=true)
        source2.kerneloffset = source2.kerneloffset_targets
        pnl.influence!((explicit1,), (source2,), direct; velocity=true)

        source1.kerneloffset = source1.kerneloffset_targets
        pnl.influence!((explicit2,), (source1,), direct; velocity=true)
        source2.kerneloffset = source2.kerneloffset_panel
        pnl.influence!((explicit2,), (source2,), direct; velocity=true)

        @test combined1.kerneloffset == combined1.kerneloffset_targets
        @test combined2.kerneloffset == combined2.kerneloffset_targets
        @test combined1.velocity ≈ explicit1.velocity atol=1e-12 rtol=1e-12
        @test combined2.velocity ≈ explicit2.velocity atol=1e-12 rtol=1e-12
    end
end
