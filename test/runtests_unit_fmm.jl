using Test
import FLOWPanel as pnl
import GeometricTools as gt

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
end
