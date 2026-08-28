using Test
import FLOWPanel as pnl
using StaticArrays, LinearAlgebra

# BRAINSTORM 025: selectable filament regularization families.
# Verifies each family's velocity against the independent infinite-filament
# closed forms, gradients against central finite differences, far-field
# agreement with the singular kernel at the family's radius_inflation
# distance, compact-support C^1 continuity at the support boundary, and the
# selection plumbing. Derivations: BRAINSTORM/025_kernel_regularization_update/
# phase_01_theory.md.

# infinite-filament transverse profiles in units of Γ/(2π) (phase_00/01 docs)
u_inf_vatistas(h, rc) = h / sqrt(h^4 + rc^4)
u_inf_compact(h, rc) = h < rc ? h / (h^2 + (h - rc)^2) : 1 / h
u_inf_gaussian(h, rc) = (1 - exp(-h^2 / (2rc^2))) / h

const REG_FAMILIES = (
    (pnl.CompactRegularization, u_inf_compact),
    (pnl.GaussianRegularization, u_inf_gaussian),
    (pnl.VatistasRegularization, u_inf_vatistas),
    # the Gaussian transverse profile is exactly LineGauss's infinite-line
    # limit (052d DERIVATION.md) — same closed form applies
    (pnl.LineGaussRegularization, u_inf_gaussian),
)

@testset "filament regularization families" begin

    rc = 0.05
    old_family = pnl.FILAMENT_REGULARIZATION[]
    try

        @testset "selection plumbing" begin
            pnl.set_filament_regularization!(:vatistas)
            @test pnl.FILAMENT_REGULARIZATION[] == pnl.VatistasRegularization
            pnl.set_filament_regularization!(:gaussian)
            @test pnl.FILAMENT_REGULARIZATION[] == pnl.GaussianRegularization
            pnl.set_filament_regularization!(:compact)
            @test pnl.FILAMENT_REGULARIZATION[] == pnl.CompactRegularization
            pnl.set_filament_regularization!(:linegauss)
            @test pnl.FILAMENT_REGULARIZATION[] == pnl.LineGaussRegularization
            @test_throws ArgumentError pnl.set_filament_regularization!(:bogus)
        end

        # long straight filament along z through the origin; unit circulation.
        # For segment half-length L >> h the mid-plane transverse speed
        # approaches the infinite-filament profile with O((h/L)^2) error.
        L = 200.0
        v1 = SVector(0.0, 0.0, -L)
        v2 = SVector(0.0, 0.0, L)
        segment_velocity(target) = pnl._bound_vortex_velocity(v1 - target, v2 - target, true, rc)

        @testset "velocity vs infinite-filament closed form ($(family))" for (family, u_inf) in REG_FAMILIES
            pnl.set_filament_regularization!(family)
            for h in (0.2rc, 0.7rc, 1.0rc, 1.5rc, 4rc, 40rc)
                target = SVector(h, 0.0, 0.0)
                u = segment_velocity(target)
                u_sing = pnl._bound_vortex_velocity(v1 - target, v2 - target, false, rc)
                # transverse (y) component only; sign must match the singular
                # kernel's (family selection must never flip circulation)
                u_expect = u_inf(h, rc) / (2π)
                @test isapprox(abs(u[2]), u_expect; rtol=1e-3)
                @test sign(u[2]) == sign(u_sing[2])
                @test abs(u[1]) < 1e-12 && abs(u[3]) < 1e-12
            end
        end

        # (Ryan directive 2026-08-20, test (a)): analytic gradient == central
        # FD of the velocity, for ALL families, at radii spanning inside-core
        # (0.3rc, 0.5rc), near the core boundary (0.9rc, 1.1rc), and far
        # field (3rc, 8rc), including skew (non-perpendicular) targets
        @testset "velocity↔gradient FD consistency, inside/near-rc/far ($(family))" for (family, _) in REG_FAMILIES
            pnl.set_filament_regularization!(family)
            for target in (SVector(0.3rc, 0.0, 0.1), SVector(0.9rc, 0.02, -0.2),
                           SVector(1.1rc, -0.03, 0.15), SVector(3rc, 0.1, 0.3),
                           SVector(0.5rc, 0.5rc, 1.0), SVector(8rc, 0.2, -0.5))
                g = pnl._bound_vortex_gradient(v1 - target, v2 - target, true, rc)
                Δ = 1e-7
                for j in 1:3
                    e = SVector{3}(1:3 .== j) * Δ
                    dudx = (segment_velocity(target + e) - segment_velocity(target - e)) / (2Δ)
                    for i in 1:3
                        # code convention: gradient[i, j] = ∂u_i/∂x_j (verified
                        # against the legacy Vatistas implementation by FD)
                        @test isapprox(g[i, j], dudx[i]; rtol=2e-5, atol=1e-8)
                    end
                end
            end
        end

        # (Ryan directive 2026-08-20, test (b)): at and beyond each family's
        # radius_inflation distance, BOTH velocity and gradient match the
        # singular kernel within the family's contract: Gaussian uses the
        # gradient-aware radius (e^{-z}(1+2z)=tol, review finding 1) so both
        # fields are within tol there; Vatistas' shipped velocity-derived
        # radius leaves velocity ≤ tol AND gradient ≤ 1.25·tol (coefficient
        # 2.5(rc/h)^4 vs velocity ½(rc/h)^4, phase_01); compact is exact.
        @testset "far-field contract at radius_inflation ($(family))" for (family, _) in REG_FAMILIES
            pnl.set_filament_regularization!(family)
            tol = 1e-5
            Δr = pnl.radius_inflation(pnl.VortexRing, rc, tol)
            grad_bound = family == pnl.VatistasRegularization ? 1.25 * tol : tol
            for h in (Δr, 2Δr)
                target = SVector(h, 0.0, 0.0)
                u_reg = segment_velocity(target)
                u_sing = pnl._bound_vortex_velocity(v1 - target, v2 - target, false, rc)
                @test norm(u_reg - u_sing) / norm(u_sing) <= 1.02 * tol
                g_reg = pnl._bound_vortex_gradient(v1 - target, v2 - target, true, rc)
                g_sing = pnl._bound_vortex_gradient(v1 - target, v2 - target, false, rc)
                @test norm(g_reg - g_sing) / norm(g_sing) <= 1.05 * grad_bound
            end
        end

        # (Ryan directive 2026-08-20): monotone convergence u_reg → u_sing as
        # h grows (relative velocity error nonincreasing beyond the core)
        @testset "monotone far-field convergence ($(family))" for (family, _) in REG_FAMILIES
            pnl.set_filament_regularization!(family)
            hs = rc .* (1.05, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0)
            errs = map(hs) do h
                target = SVector(h, 0.0, 0.0)
                u_reg = segment_velocity(target)
                u_sing = pnl._bound_vortex_velocity(v1 - target, v2 - target, false, rc)
                norm(u_reg - u_sing) / norm(u_sing)
            end
            @test all(errs[i+1] <= errs[i] + 1e-14 for i in 1:length(errs)-1)
            @test errs[end] <= max(errs[1] * 1e-3, 5e-16)
        end

        @testset "compact support is exact beyond rc" begin
            pnl.set_filament_regularization!(pnl.CompactRegularization)
            for h in (1.0001rc, 2rc, 10rc)
                target = SVector(h, 0.05, 0.2)
                u_reg = segment_velocity(target)
                u_sing = pnl._bound_vortex_velocity(v1 - target, v2 - target, false, rc)
                @test u_reg ≈ u_sing rtol=1e-14
                g_reg = pnl._bound_vortex_gradient(v1 - target, v2 - target, true, rc)
                g_sing = pnl._bound_vortex_gradient(v1 - target, v2 - target, false, rc)
                @test g_reg ≈ g_sing rtol=1e-14
            end
        end

        @testset "compact C^1 continuity at the support boundary" begin
            pnl.set_filament_regularization!(pnl.CompactRegularization)
            ε = 1e-8
            t_in = SVector(rc*(1 - ε), 0.0, 0.0)
            t_out = SVector(rc*(1 + ε), 0.0, 0.0)
            @test segment_velocity(t_in) ≈ segment_velocity(t_out) rtol=1e-6
            g_in = pnl._bound_vortex_gradient(v1 - t_in, v2 - t_in, true, rc)
            g_out = pnl._bound_vortex_gradient(v1 - t_out, v2 - t_out, true, rc)
            @test g_in ≈ g_out rtol=1e-5
        end

        @testset "on-filament and endpoint guards ($(family))" for (family, _) in REG_FAMILIES
            pnl.set_filament_regularization!(family)
            @test segment_velocity(v1) == zero(SVector{3,Float64})
            u0 = segment_velocity(SVector(0.0, 0.0, 0.3))  # on the line, interior
            @test all(isfinite, u0)
            g0 = pnl._bound_vortex_gradient(v1 - SVector(1e-14, 0.0, 0.3),
                                            v2 - SVector(1e-14, 0.0, 0.3), true, rc)
            @test all(isfinite, g0)
        end

        @testset "radius_inflation per family" begin
            tol = 1e-6
            pnl.set_filament_regularization!(pnl.CompactRegularization)
            @test pnl.radius_inflation(pnl.VortexRing, rc, tol) == rc
            pnl.set_filament_regularization!(pnl.GaussianRegularization)
            # gradient-aware radius: defining property e^{-z}(1+2z) = tol
            # with z = (Δr/rc)²/2 (review finding 1); larger than the
            # velocity-only √(2 ln 1/tol) radius
            Δr_g = pnl.radius_inflation(pnl.VortexRing, rc, tol)
            z_g = (Δr_g / rc)^2 / 2
            @test isapprox(exp(-z_g) * (1 + 2z_g), tol; rtol=1e-4)
            @test Δr_g > rc * sqrt(2 * log(1/tol))
            pnl.set_filament_regularization!(pnl.VatistasRegularization)
            @test pnl.radius_inflation(pnl.VortexRing, rc, tol) ≈ rc * (2/tol)^0.25
            @test pnl.radius_inflation(pnl.VortexRing, rc, Inf) == 0.0
            # LineGauss: Gaussian gradient-aware fixed point + 0.35rc pad
            # (segment-distance calibration, 052d k01 T7/T7b dense scans)
            pnl.set_filament_regularization!(pnl.LineGaussRegularization)
            @test pnl.radius_inflation(pnl.VortexRing, rc, tol) ≈ Δr_g + 0.35rc
            @test pnl.radius_inflation(pnl.VortexRing, rc, Inf) == 0.0
            # source/doublet rule unchanged
            @test pnl.radius_inflation(pnl.ConstantSource, rc, tol) == rc
        end

    finally
        pnl.FILAMENT_REGULARIZATION[] = old_family
    end

end
