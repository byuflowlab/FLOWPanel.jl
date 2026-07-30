using Test
import FLOWPanel as pnl
import FastMultipole
import FLOWVPM

if !isdefined(@__MODULE__, :make_octa_source_body)
    include("test_helpers.jl")
end

function sfs_backend_test_field(; SFS=FLOWVPM.ConstantSFS(FLOWVPM.Estr_fmm))
    pfield = FLOWVPM.ParticleField(4, Float64;
        transposed=true,
        fmm=FLOWVPM.FMM(;
            p=4,
            ncrit=2,
            theta=0.0,
            shrink_recenter=true,
            relative_tolerance=1e-10,
            absolute_tolerance=1e-10,
            autotune_p=false,
            autotune_ncrit=false,
            autotune_reg_error=false,
            min_ncrit=1),
        SFS)

    FLOWVPM.add_particle(pfield, (0.0, 0.0, 0.0), (1.0, 0.2, -0.1), 0.50)
    FLOWVPM.add_particle(pfield, (0.2, 0.1, 0.0), (-0.2, 1.1, 0.3), 0.45)
    FLOWVPM.add_particle(pfield, (1.5, 0.2, 0.1), (0.5, -0.4, 0.8), 0.55)
    FLOWVPM.add_particle(pfield, (-0.3, 1.0, 0.4), (-0.7, 0.3, 1.2), 0.60)
    return pfield
end

function sfs_values(pfield)
    return [copy(FLOWVPM.get_SFS(pfield, i)) for i in 1:FLOWVPM.get_np(pfield)]
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

    @testset "particle SFS postcalc hook" begin
        precalc = sfs_backend_test_field()
        precalc.particles[FLOWVPM.U_INDEX, 1:precalc.np] .= 1.25
        precalc.particles[FLOWVPM.J_INDEX, 1:precalc.np] .= -2.5
        pnl.pre_evaluate_influence!(precalc)
        @test all(==(1.25), precalc.particles[FLOWVPM.U_INDEX, 1:precalc.np])
        @test all(iszero, precalc.particles[FLOWVPM.J_INDEX, 1:precalc.np])

        reference = sfs_backend_test_field()
        FLOWVPM.UJ_direct(reference; sfs=false, reset=true, reset_sfs=true)
        FLOWVPM.Estr_direct!(reference)
        reference.SFS(reference, FLOWVPM.AfterUJ())
        reference_sfs = sfs_values(reference)

        direct = sfs_backend_test_field()
        pnl.influence!(direct, direct, pnl.DirectBackend();
            velocity=true,
            velocity_gradient=true,
            postcalc=true)
        for i in 1:FLOWVPM.get_np(direct)
            @test FLOWVPM.get_SFS(direct, i) ≈ reference_sfs[i] rtol=1e-12 atol=1e-14
        end

        fmm = sfs_backend_test_field()
        pnl.influence!(fmm, fmm, pnl.FastMultipoleBackend(
                expansion_order=3,
                multipole_acceptance=0.0,
                leaf_size=2);
            velocity=true,
            velocity_gradient=true,
            postcalc=true,
            error_tolerance=FastMultipole.PowerRelativeGradient{1e-10, 1e-10, true}(),
            tune=true)
        for i in 1:FLOWVPM.get_np(fmm)
            @test FLOWVPM.get_SFS(fmm, i) ≈ reference_sfs[i] rtol=1e-12 atol=1e-14
        end

        no_postcalc = sfs_backend_test_field()
        pnl.influence!(no_postcalc, no_postcalc, pnl.DirectBackend();
            velocity=true,
            velocity_gradient=true,
            postcalc=false)
        @test all(sfs -> all(iszero, sfs), sfs_values(no_postcalc))

        no_sfs = sfs_backend_test_field(; SFS=FLOWVPM.noSFS)
        pnl.influence!(no_sfs, no_sfs, pnl.DirectBackend();
            velocity=true,
            velocity_gradient=true,
            postcalc=true)
        @test all(sfs -> all(iszero, sfs), sfs_values(no_sfs))
    end
end
