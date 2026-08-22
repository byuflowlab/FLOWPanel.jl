using Test
import FLOWPanel as pnl
import FastMultipole
import FLOWVPM
import LinearAlgebra

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

@testset verbose=true "triaxial panel FmmPlan rigid equivariance" begin
    body = make_sphere_source_body()
    body.nodes[1, :] .*= 1.0
    body.nodes[2, :] .*= 0.7
    body.nodes[3, :] .*= 0.45
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.kerneloffset = body.kerneloffset_panel
    for i in 1:body.ncells
        body.strength[i, 1] = sin(0.017 * i) + 0.3 * cos(0.031 * i)
    end

    plan = FastMultipole.FmmPlan((body,), (body,);
        expansion_order=8, multipole_acceptance=0.4, leaf_size_source=20,
        scalar_potential=true, gradient=false, hessian=false, shrink=true)
    @test body.ncells == 2304
    @test !isempty(plan.m2l_list)
    @test !isempty(plan.direct_list)

    function direct_potential!(body)
        body.potential .= 0
        pnl.influence!((body,), (body,), pnl.DirectBackend();
            scalar_potential=true, velocity=false)
        return copy(body.potential)
    end
    function planned_potential!(body, plan)
        body.potential .= 0
        FastMultipole.fmm!((body,), (body,), plan)
        return copy(body.potential)
    end

    phi_direct_0 = direct_potential!(body)
    phi_fmm_0 = planned_potential!(body, plan)
    refnorm_0 = LinearAlgebra.norm(phi_direct_0)
    err_0 = LinearAlgebra.norm(phi_fmm_0 - phi_direct_0) / refnorm_0
    @test all(isfinite, phi_direct_0)
    @test all(isfinite, phi_fmm_0)
    @test refnorm_0 > 0
    @test err_0 < 1e-7

    axis = [0.37, -0.61, 0.70]
    axis ./= LinearAlgebra.norm(axis)
    theta = 0.73
    K = [0.0 -axis[3] axis[2];
         axis[3] 0.0 -axis[1];
         -axis[2] axis[1] 0.0]
    R = cos(theta) .* Matrix{Float64}(LinearAlgebra.I, 3, 3) .+
        (1 - cos(theta)) .* (axis * axis') .+ sin(theta) .* K
    t = [0.31, -0.27, 0.19]
    body.nodes .= R * body.nodes .+ t
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    FastMultipole.transform_plan!(plan, (body,), R, t)

    phi_direct_1 = direct_potential!(body)
    phi_fmm_1 = planned_potential!(body, plan)
    refnorm_1 = LinearAlgebra.norm(phi_direct_1)
    err_1 = LinearAlgebra.norm(phi_fmm_1 - phi_direct_1) / refnorm_1
    @test all(isfinite, phi_direct_1)
    @test all(isfinite, phi_fmm_1)
    @test refnorm_1 > 0
    @test err_1 <= 1.05 * err_0 + 1e-12
    @test LinearAlgebra.norm(phi_direct_1 - phi_direct_0) / refnorm_0 <= 1e-12
    @test LinearAlgebra.norm(phi_fmm_1 - phi_fmm_0) /
          LinearAlgebra.norm(phi_fmm_0) <= 1e-12
end

@testset verbose=true "kerneloffset-aware FMM radius" begin
    @testset "radius_inflation formulas" begin
        tol = pnl.FMM_RADIUS_TOL[]
        # source/doublet regularization is compactly supported: exact beyond koff
        @test pnl.radius_inflation(pnl.ConstantSource, 0.05, tol) == 0.05
        @test pnl.radius_inflation(pnl.ConstantDoublet, 0.05, tol) == 0.05
        @test pnl.radius_inflation(Union{pnl.ConstantSource, pnl.ConstantDoublet}, 0.05, tol) == 0.05
        # VortexRing inflation follows the active FilamentRegularization
        # family (BRAINSTORM 025); assert each family's rule explicitly
        old_family = pnl.FILAMENT_REGULARIZATION[]
        try
            # compact-support (default): exact beyond koff, tol-independent
            pnl.set_filament_regularization!(pnl.CompactRegularization)
            @test pnl.radius_inflation(pnl.VortexRing, 0.05, 1e-6) == 0.05
            @test pnl.radius_inflation(Union{pnl.ConstantSource, pnl.VortexRing}, 0.05, 1e-6) == 0.05
            # Gaussian: GRADIENT-AWARE radius (025 review finding 1) solving
            # e^{-z}(1+2z) = tol with z = (Δr/koff)²/2; exceeds the
            # velocity-only √(2 ln 1/tol) radius (≈5.90 koff at 1e-6)
            pnl.set_filament_regularization!(pnl.GaussianRegularization)
            Δr_g = pnl.radius_inflation(pnl.VortexRing, 0.05, 1e-6)
            z_g = (Δr_g / 0.05)^2 / 2
            @test isapprox(exp(-z_g) * (1 + 2z_g), 1e-6; rtol=1e-4)
            @test Δr_g > 0.05 * sqrt(2 * log(1e6))
            # Vatistas n=2 (legacy): rel error ½(koff/h)^4 ≤ tol at h ≥ koff·(2/tol)^(1/4)
            pnl.set_filament_regularization!(pnl.VatistasRegularization)
            @test pnl.radius_inflation(pnl.VortexRing, 0.05, 1e-6) ≈ 0.05 * (2e6)^0.25
            @test pnl.radius_inflation(Union{pnl.ConstantSource, pnl.VortexRing}, 0.05, 1e-6) ≈
                  0.05 * (2e6)^0.25
            # Inf disables (pre-fix behavior)
            @test pnl.radius_inflation(pnl.VortexRing, 0.05, Inf) == 0.0
        finally
            pnl.FILAMENT_REGULARIZATION[] = old_family
        end
        @test pnl.radius_inflation(pnl.ConstantSource, 0.05, Inf) == 0.0
    end

    @testset "panel buffer radius includes inflation" begin
        for (body, TK) in ((make_octa_source_body(), pnl.ConstantSource),
                           (make_plate_vortex_body(), pnl.VortexRing))
            body.kerneloffset_panel = 1e-3
            body.kerneloffset_targets = 2e-3
            buffer = zeros(FastMultipole.data_per_body(body), 1)

            # geometric radius reference: inflation disabled
            body.kerneloffset = 0.0
            FastMultipole.source_system_to_buffer!(buffer, 1, body, 1)
            r_geom = buffer[4, 1]

            # inflation follows the ACTIVE kerneloffset (the pass's offset,
            # set by _set_kerneloffsets! right before each influence call)
            for koff in (body.kerneloffset_panel, body.kerneloffset_targets)
                body.kerneloffset = koff
                FastMultipole.source_system_to_buffer!(buffer, 1, body, 1)
                @test buffer[4, 1] ≈ r_geom +
                      pnl.radius_inflation(TK, koff, pnl.fmm_radius_tolerance(body))
            end

            # FMM_RADIUS_TOL[] = Inf recovers the geometric radius
            old_tol = pnl.FMM_RADIUS_TOL[]
            try
                pnl.FMM_RADIUS_TOL[] = Inf
                body.kerneloffset = body.kerneloffset_targets
                FastMultipole.source_system_to_buffer!(buffer, 1, body, 1)
                @test buffer[4, 1] == r_geom
            finally
                pnl.FMM_RADIUS_TOL[] = old_tol
            end
        end
    end

    @testset "large-koff FMM matches regularized direct operator" begin
        # Regression for the 021 Phase 1 floor: with a kerneloffset comparable
        # to the panel size, geometric radii let the MAC admit expansions where
        # the singular and regularized kernels still disagree; the inflated
        # radius must push those pairs into the direct near field.
        function make_grid_vortex_body(shift; n=6, koff=0.4)
            xs = range(0.0, 1.0; length=n+1)
            nodes = zeros(3, (n+1)^2)
            for (j, y) in enumerate(xs), (i, x) in enumerate(xs)
                nodes[:, (j-1)*(n+1) + i] .= (x + shift[1], y + shift[2], shift[3])
            end
            cells = zeros(Int, 3, 2*n^2)
            c = 0
            for j in 1:n, i in 1:n
                p1 = (j-1)*(n+1) + i; p2 = p1 + 1; p3 = p1 + n + 2; p4 = p1 + n + 1
                cells[:, c += 1] .= (p1, p2, p3)
                cells[:, c += 1] .= (p1, p3, p4)
            end
            body = pnl.RigidWakeBody{pnl.VortexRing}(nodes, cells, pnl.noshedding;
                DBC=false, check_mesh=false, watertight=false, kerneloffset=koff)
            pnl.calc_normals!(body)
            pnl.calc_controlpoints!(body)
            for i in 1:body.ncells   # varied strengths so ring cancellations don't hide errors
                body.strength[i, 1] = 1.0 + 0.05 * i
            end
            body.velocity .= 0
            body.potential .= 0
            return body
        end
        koff = 0.4
        make_source = () -> make_grid_vortex_body((0.0, 0.0, 0.0); koff)
        target = make_grid_vortex_body((1.5, 0.2, 0.3); koff)
        direct = pnl.DirectBackend()
        fmm = pnl.FastMultipoleBackend(expansion_order=12, multipole_acceptance=0.5, leaf_size=4)

        vel_direct, _ = evaluate_velocity_and_potential(deepcopy(target), make_source(), direct)
        vel_fmm, _ = evaluate_velocity_and_potential(deepcopy(target), make_source(), fmm)
        err = maximum(abs.(vel_fmm .- vel_direct)) / maximum(abs.(vel_direct))

        # premise guard: with the inflation disabled the same comparison must
        # show the O((koff/d)^4) expansion mismatch, else this test is vacuous
        old_tol = pnl.FMM_RADIUS_TOL[]
        err_uninflated = try
            pnl.FMM_RADIUS_TOL[] = Inf
            vel_fmm_inf, _ = evaluate_velocity_and_potential(deepcopy(target), make_source(), fmm)
            maximum(abs.(vel_fmm_inf .- vel_direct)) / maximum(abs.(vel_direct))
        finally
            pnl.FMM_RADIUS_TOL[] = old_tol
        end

        @test err <= 10 * pnl.FMM_RADIUS_TOL[]
        @test err_uninflated > 10 * pnl.FMM_RADIUS_TOL[]
        @test err < 0.1 * err_uninflated
    end

    @testset "sliver panel near-plane far-field potential (edge-extension guard)" begin
        # Regression for the 2026-08-14 direct-kernel defect: the tan_term guard
        # used an ABSOLUTE 1e-12 tolerance (units length²) to detect targets on a
        # panel-side extension, zeroing one edge's solid-angle term for far
        # targets within ~√(1e-12/(ri·ds)) rad of the extension line of a sliver
        # panel — spurious potential ~10³× the truth (measured 3,600× on the DJI
        # TE strip; surfaced as a fake p-saturated "FMM residual floor").
        # DJI-TE-strip-like sliver + target near the long-edge extension, near-plane:
        v1 = pnl.SVector{3}(0.0, 0.0, 0.0)
        v2 = pnl.SVector{3}(5.36e-3, 0.0, 0.0)
        v3 = pnl.SVector{3}(1.0e-4, 2.0e-4, 0.0)
        xt = pnl.SVector{3}(0.034, 3.0e-6, 2.4e-6)

        phi_quad = let acc = 0.0, nref = 400
            nrm = LinearAlgebra.cross(v2 - v1, v3 - v1)
            nhat = nrm / LinearAlgebra.norm(nrm)
            dA = 0.5 * LinearAlgebra.norm(nrm) / nref^2
            for a in 0:nref-1, b in 0:nref-1-a, (da, db) in ((1/3, 1/3), (2/3, 2/3))
                da == 2/3 && a + b >= nref - 1 && continue
                u = (a + da) / nref; v = (b + db) / nref
                ξ = v1 + u * (v2 - v1) + v * (v3 - v1)
                r = xt - ξ
                acc += sum(nhat .* r) / LinearAlgebra.norm(r)^3 * dA
            end
            -acc / (4π)
        end

        R, _ = pnl.rotate_to_panel(v1[1], v1[2], v1[3], v2[1], v2[2], v2[3], v3[1], v3[2], v3[3])
        centroid = (v1 + v2 + v3) / 3
        phi, _, _ = pnl._induced(xt, (v1, v2, v3), centroid, pnl.SVector{1}(1.0),
            pnl.ConstantDoublet, 0.0, R, FastMultipole.DerivativesSwitch(true, false, false))
        @test isapprox(phi, phi_quad; rtol=0.05)
    end

    @testset "PanelWake buffer radius includes core_size reach" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelWake(body; nwakerows=2, core_size=1e-2)
        # give the (empty) wake one row of plausible geometry so the buffer
        # writer has a panel to read; the writer touches only nodes/strengths
        wake.nwakes[] = 1
        for isurf in eachindex(wake.nodes)
            ncols = size(wake.nodes[isurf], 3)
            for icol in 1:ncols, irow in 1:2
                wake.nodes[isurf][:, irow, icol] .= (icol - 1.0, irow - 1.0, 0.0)
            end
        end
        buffer = zeros(FastMultipole.data_per_body(wake), 1)
        FastMultipole.source_system_to_buffer!(buffer, 1, wake, 1)
        r_expected_geom = 0.5 * sqrt(2.0)   # half-diagonal of the unit quad above
        TK = typeof(wake).parameters[1]
        @test buffer[4, 1] ≈ r_expected_geom +
              pnl.radius_inflation(TK, wake.core_size, pnl.fmm_radius_tolerance(wake))
    end
end
