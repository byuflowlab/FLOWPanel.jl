using Test
import FLOWPanel as pnl
import FLOWVPM
const TOML = pnl.TOML
import FastMultipole
using LinearAlgebra: norm, cross
using StaticArrays: SVector

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

include(joinpath(@__DIR__, "data", "legacy_wake_conversion_reference.jl"))

@testset verbose=true "Free Wakes" begin
    @testset "SigmaOverlap chooses particle count from fixed sigma" begin
        pfield = FLOWVPM.ParticleField(10, Float64;
            fmm=FLOWVPM.FMM(autotune_reg_error=false))
        r1 = SVector(0.0, 0.0, 0.0)
        r2 = SVector(1.0, 0.0, 0.0)
        sigma = 0.2
        target_overlap = 0.7

        pnl._shed_particles!(pfield, r1, r2, 1.0, pnl.SigmaOverlap(sigma, target_overlap))

        @test pfield.np == 4
        @test all(pfield.particles[FLOWVPM.SIGMA_INDEX, 1:pfield.np] .== sigma)
        @test sigma * pfield.np / norm(r2 - r1) >= target_overlap
        @test pfield.particles[1, 1:pfield.np] ≈ [0.125, 0.375, 0.625, 0.875]
    end

    @testset "StationSigmaOverlap sheds per-station sigma (018 Phase 16)" begin
        # A constant per-station vector must be indistinguishable from the
        # span-uniform SigmaOverlap it resolves to — particle-by-particle.
        ov = 0.7
        for wrap in (false, true), nwakerows in (1, 2)
            ref = make_conversion_fixture(; nwakerows, wrap,
                method_trailing=pnl.SigmaOverlap(0.2, ov),
                method_unsteady=pnl.SigmaOverlap(0.2, ov))
            pnl._convert_to_particles!(ref)

            ncols1 = size(ref.panel_wake.nodes[1], 3)
            station = make_conversion_fixture(; nwakerows, wrap,
                method_trailing=pnl.StationSigmaOverlap([fill(0.2, ncols1)], ov),
                method_unsteady=pnl.StationSigmaOverlap([fill(0.2, ncols1)], ov))
            pnl._convert_to_particles!(station)

            np = ref.pfield.np
            @test station.pfield.np == np
            @test station.pfield.particles[:, 1:np] == ref.pfield.particles[:, 1:np]
        end

        # Varying sigma: trailing filaments at node column j (fixture geometry
        # puts column j at y = j-1) must carry sigmas[j]; unsteady filaments
        # (spanning y in (j-1, j)) the mean of columns j and j+1.
        wake = make_conversion_fixture(; nwakerows=1, wrap=false)
        ncols1 = size(wake.panel_wake.nodes[1], 3)
        sig = collect(range(0.1, 0.3; length=ncols1))
        wake = make_conversion_fixture(; nwakerows=1, wrap=false,
            method_trailing=pnl.StationSigmaOverlap([sig], 0.7),
            method_unsteady=pnl.StationSigmaOverlap([sig], 0.7))
        pnl._convert_to_particles!(wake)
        pf = wake.pfield
        @test pf.np > 0
        for i in 1:pf.np
            y = pf.particles[2, i]
            s = pf.particles[FLOWVPM.SIGMA_INDEX, i]
            if abs(y - round(y)) < 1e-9
                @test s ≈ sig[Int(round(y)) + 1]              # trailing
            else
                j = Int(floor(y)) + 1
                @test s ≈ 0.5 * (sig[j] + sig[j + 1])         # unsteady
            end
        end
        # The spanwise sigma range must actually be exercised.
        sigmas_shed = pf.particles[FLOWVPM.SIGMA_INDEX, 1:pf.np]
        @test minimum(sigmas_shed) < 0.15 && maximum(sigmas_shed) > 0.25

        # Wrong-length sigma vector fails loudly at conversion.
        bad = make_conversion_fixture(; nwakerows=1, wrap=false,
            method_trailing=pnl.StationSigmaOverlap([sig[1:end-1]], 0.7),
            method_unsteady=pnl.NoShed())
        @test_throws ErrorException pnl._convert_to_particles!(bad)

        # Constructor validation.
        @test_throws ArgumentError pnl.StationSigmaOverlap(Vector{Float64}[], 0.7)
        @test_throws ArgumentError pnl.StationSigmaOverlap([[0.1, -0.2]], 0.7)
        @test_throws ArgumentError pnl.StationSigmaOverlap([[0.1, NaN]], 0.7)
        @test_throws ArgumentError pnl.StationSigmaOverlap([[0.1, 0.2]], 0.0)

        # Unresolved direct shedding is an error; identity dispatch returns
        # the very same object for span-uniform methods.
        @test_throws ArgumentError pnl._shed_particles!(nothing,
            SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), 1.0,
            pnl.StationSigmaOverlap([[0.1, 0.2]], 0.7))
        m = pnl.OverlapPPS(1.3, 2)
        @test pnl._station_method(m, 1, 3, 4) === m
    end

    @testset "Viscous scheme actually runs under the euler stepper (018 erratum)" begin
        # FLOWPanel steps its particle wake with FLOWVPM._euler (see step!), so
        # the pfield must DECLARE integration=euler: viscousdiffusion branches
        # on pfield.integration, and under the FLOWVPM default (rungekutta3) a
        # CoreSpreading scheme silently no-ops — zeroed RK stage weights mean
        # no core spreading and no beta resets. Every 018 campaign run before
        # 2026-08-03 was functionally inviscid for this reason.
        body = make_dirichlet_diamond_body(; nspan=3)
        nu, sgm0 = 1.2e-5, 0.05
        wake = pnl.PanelParticleWake(body; nwakerows=2, max_particles=10,
            viscous=FLOWVPM.CoreSpreading(nu, sgm0, FLOWVPM.zeta_fmm; beta=1e6))
        pf = wake.pfield
        @test pf.integration === FLOWVPM.euler

        sigma0 = 0.04
        FLOWVPM.add_particle(pf, SVector(0.0, 0.0, 0.0),
            SVector(1e-3, 0.0, 0.0), sigma0; circulation=1e-3)
        dt = 1e-3
        FLOWVPM.viscousdiffusion(pf, dt)
        @test pf.particles[FLOWVPM.SIGMA_INDEX, 1] ≈ sqrt(sigma0^2 + 2 * nu * dt)
    end

    @testset "RK3 wake integrator wiring (BRAINSTORM 026 Phase 1b)" begin
        body = make_dirichlet_diamond_body(; nspan=3)

        # rk3=true declares the RK3 scheme on the pfield; combining with
        # expint is rejected loudly
        wake = pnl.PanelParticleWake(body; nwakerows=2, max_particles=10,
            rk3=true)
        @test wake.pfield.integration === FLOWVPM.rungekutta3
        @test_throws ArgumentError pnl.PanelParticleWake(body; nwakerows=2,
            max_particles=10, rk3=true, expint=true)

        # propagate! on an rk3 wake demands the stage-UJ closure
        @test_throws ArgumentError pnl.propagate!(wake, 1e-3)

        # _rk3_convect! consistency: a constant velocity field with zero
        # gradient must advance X by exactly dt*U (the low-storage weights
        # sum to 1), leave Gamma/sigma untouched, and with CoreSpreading
        # under the RK3 aux weights reproduce sigma^2 += 2*nu*dt exactly.
        nu, sigma0 = 1.2e-5, 0.04
        wake_cs = pnl.PanelParticleWake(body; nwakerows=2, max_particles=10,
            rk3=true,
            viscous=FLOWVPM.CoreSpreading(nu, sigma0, FLOWVPM.zeta_fmm; beta=1e6))
        pf = wake_cs.pfield
        FLOWVPM.add_particle(pf, SVector(0.0, 0.0, 0.0),
            SVector(1e-3, 0.0, 0.0), sigma0; circulation=1e-3)
        U0 = [0.3, -0.2, 0.1]
        set_stage_velocity! = (pfield) -> begin
            pfield.particles[FLOWVPM.U_INDEX, 1:pfield.np] .= U0
            pfield.particles[FLOWVPM.J_INDEX, 1:pfield.np] .= 0.0
            nothing
        end
        set_stage_velocity!(pf)
        stage_calls = Ref(0)
        stage_UJ = (pfield, a, b) -> begin
            stage_calls[] += 1
            set_stage_velocity!(pfield)
            nothing
        end
        dt = 1e-3
        gamma0 = copy(pf.particles[FLOWVPM.GAMMA_INDEX, 1])
        pnl._rk3_convect!(pf, dt; relax=false, stage_UJ)
        @test stage_calls[] == 2   # stage 1 reuses the step's evaluation
        @test pf.particles[1:3, 1] ≈ U0 .* dt rtol = 1e-14
        @test pf.particles[FLOWVPM.GAMMA_INDEX, 1] ≈ gamma0 rtol = 1e-14
        @test pf.particles[FLOWVPM.SIGMA_INDEX, 1] ≈
              sqrt(sigma0^2 + 2 * nu * dt) rtol = 1e-12
    end

    @testset "Legacy edge-jump conversion golden reference (BRAINSTORM 016)" begin
        # Characterization baseline for the Item 016 refactor, which moves the
        # legacy conversion verbatim behind an opt-in strategy dispatch. The
        # fixture writes panel-wake rows directly and converts once, so this
        # pins `_convert_to_particles!` itself with no time stepping in the way.
        # Equality is exact by intent: the legacy path must not move at all.
        for wrap in (false, true), nwakerows in (1, 2, 3)
            ref = LEGACY_CONVERSION_REFERENCE[(wrap, nwakerows)]
            wake = make_conversion_fixture(; nwakerows=nwakerows, wrap=wrap)
            pnl._convert_to_particles!(wake)
            pf = wake.pfield
            np = pf.np

            @test np == ref.np
            @test pf.particles[1:3, 1:np] == ref.X
            @test pf.particles[FLOWVPM.GAMMA_INDEX, 1:np] == ref.G
            @test pf.particles[FLOWVPM.SIGMA_INDEX, 1:np] == ref.S
            @test pf.particles[FLOWVPM.CIRCULATION_INDEX, 1:np] == ref.C
        end

        # A wrapping chain has no physical root or tip, so the streamwise
        # closure filament deposited on the open chain must be absent.
        for nwakerows in (1, 2, 3)
            @test LEGACY_CONVERSION_REFERENCE[(true, nwakerows)].np <
                  LEGACY_CONVERSION_REFERENCE[(false, nwakerows)].np
        end
    end

    @testset "Panel-particle conversion strategy axis (BRAINSTORM 016)" begin
        body = make_dirichlet_diamond_body(; nspan=3)

        # Omitting `conversion` and asking for the legacy strategy explicitly
        # must be indistinguishable, down to the resolved line policies.
        implicit = pnl.PanelParticleWake(body; nwakerows=2, max_particles=100)
        explicit = pnl.PanelParticleWake(body; nwakerows=2, max_particles=100,
            conversion=pnl.LegacyEdgeJumpConversion())

        @test implicit.conversion === pnl.LegacyEdgeJumpConversion()
        @test typeof(implicit) === typeof(explicit)
        @test implicit.method_trailing == explicit.method_trailing == pnl.OverlapPPS(1.3, 2)
        @test implicit.method_unsteady == explicit.method_unsteady == pnl.OverlapPPS(1.3, 2)
        @test implicit.conversion_workspace === nothing
        @test implicit.conversion_diagnostics === nothing
        @test implicit.conversion_count[] == 0

        # An explicitly supplied legacy line policy is still honored.
        custom = pnl.PanelParticleWake(body; nwakerows=2, max_particles=100,
            method_trailing=pnl.SigmaPPS(0.3, 4))
        @test custom.method_trailing == pnl.SigmaPPS(0.3, 4)
        @test custom.method_unsteady == pnl.OverlapPPS(1.3, 2)

        # Smooth strategy: parameter validation.
        c = pnl.SurfaceVorticityConversion(0.25)
        @test c.sigma == 0.25
        @test c.overlap == 1.3
        @test c.diagnose_nearfield === false
        @test c isa pnl.SurfaceVorticityConversion{Float64}
        # Mixed-type keywords promote to a single float type.
        @test pnl.SurfaceVorticityConversion(1; overlap=2) isa
              pnl.SurfaceVorticityConversion{Float64}

        @test_throws ArgumentError pnl.SurfaceVorticityConversion(0.0)
        @test_throws ArgumentError pnl.SurfaceVorticityConversion(-1.0)
        @test_throws ArgumentError pnl.SurfaceVorticityConversion(NaN)
        @test_throws ArgumentError pnl.SurfaceVorticityConversion(0.25; overlap=0.0)
        @test_throws ArgumentError pnl.SurfaceVorticityConversion(0.25; rank_rtol=-1.0)
        @test_throws ArgumentError pnl.SurfaceVorticityConversion(0.25; geometry_rtol=0.0)

        # The smooth strategy owns its own line sampling, so supplying a legacy
        # line policy is a configuration error rather than a silent no-op --
        # even when its value equals the legacy default.
        @test_throws ArgumentError pnl.PanelParticleWake(body; nwakerows=2,
            max_particles=100, conversion=c, method_trailing=pnl.OverlapPPS(1.3, 2))
        @test_throws ArgumentError pnl.PanelParticleWake(body; nwakerows=2,
            max_particles=100, conversion=c, method_unsteady=pnl.OverlapPPS(1.3, 2))

        # The unresolved sentinel must never reach the shedding path.
        @test_throws ArgumentError pnl._shed_particles!(nothing,
            SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), 1.0,
            pnl.DefaultWakeSheddingMethod())
    end

    @testset "Surface-vorticity reconstruction core (BRAINSTORM 016)" begin
        rtol = sqrt(eps(Float64))

        # Helper: reconstruct on a panel with unit normal `n` whose stored
        # strength is the exact affine field muhat(x) = dot(g, x), sampled at
        # two centroid offsets `d1`, `d2` lying in the tangent plane. An exact
        # reconstruction must return `g` itself.
        function reconstruct_affine(n, g, d1, d2)
            return pnl._reconstruct_surface_gradient(
                SVector{3,Float64}(n), SVector{3,Float64}(d1), dot(g, d1),
                SVector{3,Float64}(d2), dot(g, d2), rtol, rtol, 1.0)
        end

        @testset "constant and affine planar fields" begin
            n = SVector(0.0, 0.0, 1.0)
            d1 = SVector(0.5, 0.0, 0.0)   # streamwise
            d2 = SVector(0.0, 0.8, 0.0)   # spanwise

            # Constant field: zero gradient, full rank geometry.
            r = reconstruct_affine(n, SVector(0.0, 0.0, 0.0), d1, d2)
            @test r.rank == 2
            @test norm(r.gradient) < 1e-14

            for g in (SVector(2.0, 0.0, 0.0),      # streamwise only
                      SVector(0.0, -3.0, 0.0),     # spanwise only
                      SVector(2.0, -3.0, 0.0))     # combined
                r = reconstruct_affine(n, g, d1, d2)
                @test r.rank == 2
                @test r.gradient ≈ g
                @test isfinite(r.condition)
                @test all(r.observable)
                # Only the tangential part is observable; a normal component
                # must never appear in the reconstruction.
                @test abs(dot(r.gradient, n)) < 1e-14
            end
        end

        @testset "package sign convention" begin
            # muhat = -mu, so kappa = n x grad(mu) = -n x grad(muhat).
            n = SVector(0.0, 0.0, 1.0)
            grad_muhat = SVector(1.0, 0.0, 0.0)
            @test pnl._surface_vorticity(n, grad_muhat) ≈ SVector(0.0, -1.0, 0.0)
            # Antisymmetry in the strength sign.
            @test pnl._surface_vorticity(n, -grad_muhat) ≈
                  -pnl._surface_vorticity(n, grad_muhat)
        end

        @testset "rigid-rotation covariance" begin
            # A rigid rotation of the panel, its stencil, and its strength field
            # must rotate the reconstructed physical vector -- even though the
            # deterministic tangent basis may pick a different Cartesian axis.
            theta = 0.7
            axis = SVector(1.0, 2.0, -0.5) / norm(SVector(1.0, 2.0, -0.5))
            K = [0.0 -axis[3] axis[2]; axis[3] 0.0 -axis[1]; -axis[2] axis[1] 0.0]
            R = I + sin(theta) * K + (1 - cos(theta)) * K^2

            n = SVector(0.0, 0.0, 1.0)
            d1 = SVector(0.5, 0.1, 0.0)
            d2 = SVector(-0.2, 0.8, 0.0)
            g = SVector(2.0, -3.0, 0.0)

            base = reconstruct_affine(n, g, d1, d2)
            rot = reconstruct_affine(R * n, R * g, R * d1, R * d2)

            @test rot.rank == base.rank
            @test rot.gradient ≈ R * base.gradient
            @test rot.singular_values ≈ base.singular_values
            # And the derived surface vorticity is covariant too.
            @test pnl._surface_vorticity(SVector{3}(R * n), rot.gradient) ≈
                  R * pnl._surface_vorticity(n, base.gradient)
        end

        @testset "rank deficiency is diagnostic, not fatal" begin
            n = SVector(0.0, 0.0, 1.0)
            g = SVector(2.0, -3.0, 0.0)

            # Rank 1: both stencil legs parallel. The observable component is
            # recovered; the orthogonal one is flagged unobservable and the
            # minimum-norm solution carries no content there.
            d1 = SVector(0.5, 0.0, 0.0)
            d2 = SVector(1.0, 0.0, 0.0)
            r = reconstruct_affine(n, g, d1, d2)
            @test r.rank == 1
            @test r.condition == Inf
            @test r.observable == SVector(true, false)
            @test dot(r.gradient, SVector(1.0, 0.0, 0.0)) ≈ 2.0   # observable part
            @test abs(dot(r.gradient, SVector(0.0, 1.0, 0.0))) < 1e-12  # min-norm
            # Minimum norm means no larger solution is returned.
            @test norm(r.gradient) <= norm(g) + 1e-12

            # Rank 0: both legs vanish in the tangent plane (purely normal
            # offsets), so nothing about the surface gradient is observable.
            r0 = pnl._reconstruct_surface_gradient(n,
                SVector(0.0, 0.0, 0.5), 0.0, SVector(0.0, 0.0, -0.3), 0.0,
                rtol, rtol, 1.0)
            @test r0.rank == 0
            @test norm(r0.gradient) == 0
            @test r0.observable == SVector(false, false)
            @test r0.condition == Inf
        end

        @testset "vanishing metric scale is rejected" begin
            n = SVector(0.0, 0.0, 1.0)
            @test_throws pnl.WakeGeometryError pnl._reconstruct_surface_gradient(
                n, SVector(0.0, 0.0, 0.0), 0.0, SVector(0.0, 0.0, 0.0), 0.0,
                rtol, rtol, 1.0)
            @test_throws pnl.WakeGeometryError pnl._reconstruct_surface_gradient(
                n, SVector(NaN, 0.0, 0.0), 0.0, SVector(0.0, 1.0, 0.0), 0.0,
                rtol, rtol, 1.0)
        end

        @testset "bilinear geometry, area, and subdivision" begin
            # Unit square in the package vertex order.
            v1 = SVector(0.0, 0.0, 0.0)
            v2 = SVector(1.0, 0.0, 0.0)
            v3 = SVector(1.0, 1.0, 0.0)
            v4 = SVector(0.0, 1.0, 0.0)

            @test pnl._bilinear_position(v1, v2, v3, v4, 0.0, 0.0) ≈ v1
            @test pnl._bilinear_position(v1, v2, v3, v4, 1.0, 0.0) ≈ v2
            @test pnl._bilinear_position(v1, v2, v3, v4, 1.0, 1.0) ≈ v3
            @test pnl._bilinear_position(v1, v2, v3, v4, 0.0, 1.0) ≈ v4

            ref = pnl._validate_wake_panel(v1, v2, v3, v4, rtol)
            @test ref ≈ 1.0

            # Whole-panel area recovered exactly by the Gauss rule.
            @test pnl._subcell_area(v1, v2, v3, v4, 0.0, 1.0, 0.0, 1.0, rtol, ref) ≈ 1.0

            # Subcell areas partition the panel exactly.
            nxi, neta = pnl._subdivision_counts(v1, v2, v3, v4, 0.25)
            @test (nxi, neta) == (4, 4)
            total = sum(pnl._subcell_area(v1, v2, v3, v4,
                            (i - 1) / nxi, i / nxi, (j - 1) / neta, j / neta, rtol, ref)
                        for i in 1:nxi, j in 1:neta)
            @test total ≈ 1.0

            # Resolution floor: at least one subdivision per direction.
            @test pnl._subdivision_counts(v1, v2, v3, v4, 1e6) == (1, 1)

            # Nonuniform / stretched panel: counts follow the longer edge pair.
            w2 = SVector(3.0, 0.0, 0.0)
            w3 = SVector(3.0, 1.0, 0.0)
            @test pnl._subdivision_counts(v1, w2, w3, v4, 1.0) == (3, 1)

            # Warped (nonplanar) but consistently oriented quad is supported,
            # and its area exceeds the flat projection.
            p3 = SVector(1.0, 1.0, 0.5)
            refw = pnl._validate_wake_panel(v1, v2, p3, v4, rtol)
            @test pnl._subcell_area(v1, v2, p3, v4, 0.0, 1.0, 0.0, 1.0, rtol, refw) > 1.0

            # Refining the quadrature converges on the warped panel.
            coarse = pnl._subcell_area(v1, v2, p3, v4, 0.0, 1.0, 0.0, 1.0, rtol, refw)
            fine = sum(pnl._subcell_area(v1, v2, p3, v4,
                           (i - 1) / 16, i / 16, (j - 1) / 16, j / 16, rtol, refw)
                       for i in 1:16, j in 1:16)
            finer = sum(pnl._subcell_area(v1, v2, p3, v4,
                            (i - 1) / 32, i / 32, (j - 1) / 32, j / 32, rtol, refw)
                        for i in 1:32, j in 1:32)
            @test abs(finer - fine) < abs(fine - coarse)
        end

        @testset "degenerate and folded panels are rejected" begin
            v1 = SVector(0.0, 0.0, 0.0)
            v2 = SVector(1.0, 0.0, 0.0)
            v3 = SVector(1.0, 1.0, 0.0)
            v4 = SVector(0.0, 1.0, 0.0)

            # All vertices coincide.
            @test_throws pnl.WakeGeometryError pnl._validate_wake_panel(v1, v1, v1, v1, rtol)
            # Collapsed to a line: nonzero extent but vanishing metric.
            @test_throws pnl.WakeGeometryError pnl._validate_wake_panel(
                v1, v2, v2, v1, rtol)
            # Bow-tie: swapping two vertices folds the panel and reverses the
            # normal partway across it.
            @test_throws pnl.WakeGeometryError pnl._validate_wake_panel(
                v1, v2, v4, v3, rtol)
            # Non-finite vertex.
            @test_throws pnl.WakeGeometryError pnl._validate_wake_panel(
                v1, v2, v3, SVector(0.0, NaN, 0.0), rtol)
        end
    end


    @testset "Surface-vorticity conversion transaction (BRAINSTORM 016)" begin
        sigma = 0.2
        smooth(s=sigma; attribution=:upstream) =
            pnl.SurfaceVorticityConversion(s; attribution=attribution)
        # The fixture builds this body internally; rebuild an identical one so
        # the single-row upstream *strength* has a source.
        shed_body() = make_dirichlet_diamond_body(; nspan=3)

        fixture(; nwakerows, wrap, s=sigma, attribution=:upstream, kwargs...) =
            make_conversion_fixture(; nwakerows=nwakerows, wrap=wrap,
                max_particles=40000, conversion=smooth(s; attribution), kwargs...)

        # A stretched, warped, non-uniform sheet with a non-affine strength
        # field -- the geometry on which a centroid-stencil reconstruction
        # would *not* conserve circulation.
        stretched_nodes(irow, icol) = [0.4 * (1.35^(irow - 1) - 1) / 0.35,
                                       0.8 * (icol - 1) + 0.15 * (icol - 1)^2,
                                       0.05 * (irow - 1) * (icol - 1)]
        stretched_mu(irow, icol) =
            0.1irow + 0.3icol + 0.07 * irow * icol - 0.04 * icol^2
        stretched(; nwakerows, wrap, s=sigma, attribution=:upstream) =
            fixture(; nwakerows=nwakerows, wrap=wrap, s=s, attribution=attribution,
                node_fun=stretched_nodes, strength_fun=stretched_mu)

        # The exact filament content a conversion must transfer, computed from
        # stored strengths and nodes alone: the attributed share of each
        # streamwise-neighbour face plus every spanwise face jump times its edge
        # vector. Independent of the deposition code.
        # `startup` is the handoff state *before* the conversion (the flag is
        # already flipped by the time this is called): before the first handoff
        # the aft face is the sheet's physical trailing boundary and is
        # deposited whole.
        function expected_transfer(wake; alpha=1.0, startup=true)
            pw = wake.panel_wake
            N = pw.nwakes[]
            beta = startup ? 1.0 : 1 - alpha
            total = SVector(0.0, 0.0, 0.0)
            for i_surf in eachindex(pw.nodes)
                nodes = pw.nodes[i_surf]; st = pw.strength[i_surf]
                nc = size(nodes, 3) - 1
                node(r, c) = SVector{3}(nodes[1,r,c], nodes[2,r,c], nodes[3,r,c])
                wraps = norm(node(N,1) - node(N,nc+1)) < 5eps()
                muA(j) = N >= 2 ? st[1,N-1,j] :
                    (-)(pnl._get_wakestrength_mu(shed_body(), j, i_surf)...)
                for j in 1:nc
                    total += alpha * (muA(j) - st[1,N,j]) * (node(N,j+1) - node(N,j))
                    total += beta * (st[1,N,j] - st[1,N+1,j]) *
                             (node(N+1,j+1) - node(N+1,j))
                end
                # every streamwise face, each counted exactly once
                for c in 1:(wraps ? nc : nc + 1)
                    left = c == 1 ? (wraps ? st[1,N,nc] : 0.0) : st[1,N,c-1]
                    right = c == nc + 1 ? 0.0 : st[1,N,c]
                    total += (right - left) * (node(N+1,c) - node(N,c))
                end
            end
            return total
        end

        netG(w) = SVector{3}(
            sum(w.pfield.particles[FLOWVPM.GAMMA_INDEX[k], 1:w.pfield.np])
            for k in 1:3)

        @testset "smooth conversion requires its final-filament sources" begin
            body = shed_body()
            conversion = smooth()
            @test_throws ArgumentError pnl.PanelParticleWake(body;
                conversion, unsteady_filament=false)
            @test_throws ArgumentError pnl.PanelParticleWake(body;
                conversion, include_final_filament=false)
            @test_throws ArgumentError pnl.PanelParticleWake(body;
                conversion, unsteady_filament=false,
                include_final_filament=false)

            # These remain valid legacy PanelWake options. In particular, the
            # new validation must not narrow the historical conversion API.
            legacy = pnl.PanelParticleWake(body;
                unsteady_filament=false, include_final_filament=false)
            @test !legacy.panel_wake.unsteady_filament
            @test !legacy.panel_wake.include_final_filament

            # Establish the source-level reason for both guards independently
            # of the smooth deposition code. A panel ring contributes +mu_G on
            # its aft edge; only an emitted FilamentWrapper can cancel it.
            function aft_face_content(w)
                pw = w.panel_wake
                emits_filament = any(s -> s isa pnl.FilamentWrapper,
                                     pnl.get_sources(pw))
                total = zero(SVector{3,Float64})
                for is in eachindex(pw.nodes)
                    nd = pw.nodes[is]; st = pw.strength[is]
                    N = pw.nwakes[]; nc = size(nd, 3) - 1
                    for j in 1:nc
                        edge = SVector{3}(nd[:,N+1,j+1]) -
                               SVector{3}(nd[:,N+1,j])
                        filament = emits_filament ?
                            pnl._final_filament_strength(pw, is, N, j) : 0.0
                        total += (st[1,N,j] + filament) * edge
                    end
                end
                return total
            end

            function source_fixture(; unsteady_filament=true,
                                    include_final_filament=true)
                prescribed = make_conversion_fixture(; nwakerows=3, wrap=false)
                w = pnl.PanelParticleWake(body; nwakerows=3,
                    unsteady_filament, include_final_filament)
                w.panel_wake.nwakes[] = prescribed.panel_wake.nwakes[]
                for is in eachindex(w.panel_wake.nodes)
                    w.panel_wake.nodes[is] .= prescribed.panel_wake.nodes[is]
                    w.panel_wake.strength[is] .= prescribed.panel_wake.strength[is]
                end
                return w
            end
            startup = source_fixture(; unsteady_filament=true)
            steady_panel = source_fixture(; unsteady_filament=false)
            no_filament = source_fixture(; include_final_filament=false)
            no_filament.panel_wake.particle_handoff_active[] = true
            no_filament.panel_wake.particle_handoff_weight[] = 1.0

            @test norm(aft_face_content(startup)) > 0.1
            @test aft_face_content(steady_panel) ≈ zero(SVector{3,Float64}) atol=1e-12
            @test norm(aft_face_content(no_filament)) > 0.1
        end

        @testset "circulation is conserved exactly on arbitrary geometry" begin
            # The headline invariant: what leaves the panel sheet arrives in the
            # particle field, to round-off, on ANY mesh -- stretched, warped,
            # non-affine -- and at every `nwakerows`. A centroid-stencil
            # reconstruction satisfies this only on a uniform mesh.
            for wrap in (false, true), nwakerows in (1, 2, 3),
                    build in (fixture, stretched)
                wake = build(; nwakerows=nwakerows, wrap=wrap)
                pnl._convert_to_particles!(wake, shed_body())
                rec = wake.conversion_diagnostics.records[end]

                @test rec.residual_abs < 1e-12
                @test rec.residual_rel < 1e-12
                @test rec.deposited_total ≈ rec.expected_total atol=1e-12
                @test netG(wake) ≈ rec.deposited_total atol=1e-12
                # ...and against an independent recomputation from stored state
                @test rec.expected_total ≈ expected_transfer(wake) atol=1e-12
            end
        end

        @testset "kappa is the face content per unit area" begin
            wake = fixture(; nwakerows=3, wrap=false)
            # Steady-state conversion: a startup conversion additionally
            # deposits the whole aft face (covered by its own testset), which
            # would double the streamwise part here.
            wake.panel_wake.particle_handoff_active[] = true
            pnl._convert_to_particles!(wake, shed_body())
            rec = wake.conversion_diagnostics.records[end]

            # Uniform affine fixture: muhat = 0.2x + 0.3y + c on z = 0, n = +z.
            kappa_exact = -cross(SVector(0.0, 0.0, 1.0), SVector(0.2, 0.3, 0.0))
            interior = rec.panels[2]
            @test interior.kappa_conservative ≈ kappa_exact atol=1e-12
            # On a uniform mesh the divergence form IS the centered difference.
            @test interior.kappa_difference < 1e-12

            # Root/tip panels carry half the spanwise part: the other half is
            # the boundary line. That is exactly what makes the total close.
            for p in (rec.panels[1], rec.panels[end])
                @test p.kappa_conservative[1] ≈ kappa_exact[1] / 2 atol=1e-12
                @test p.kappa_conservative[2] ≈ kappa_exact[2] atol=1e-12
                @test p.kappa_difference ≈ abs(kappa_exact[1]) / 2 atol=1e-12
            end

            # Every area particle is kappa_j * dA_p with the panel's own kappa,
            # and the scalar circulation is |Gamma| / sqrt(dA).
            pf = wake.pfield
            @test all(pf.particles[FLOWVPM.SIGMA_INDEX, 1:pf.np] .== sigma)
            for i in 1:pf.np
                cls = wake.conversion_workspace.classes[i]
                cls == pnl.PerimeterLineParticle && continue
                G = SVector{3}(pf.particles[FLOWVPM.GAMMA_INDEX, i])
                k = cls == pnl.InteriorSurfaceParticle ?
                    interior.kappa_conservative : rec.panels[1].kappa_conservative
                dA = norm(G) / norm(k)
                @test pf.particles[FLOWVPM.CIRCULATION_INDEX, i] ≈ norm(G) / sqrt(dA)
            end
            @test sum(p.kappa_conservative * p.area for p in rec.panels) ≈
                  rec.total_area_strength atol=1e-12
        end

        @testset "the two kappa rules converge under refinement" begin
            # Ryan's question: does the divergence form converge to the
            # reconstruction? On a uniform mesh they are identical; on a
            # smoothly graded one they differ by O(row-to-row stretch), which
            # shrinks at first order as the sheet is refined.
            #
            # Refine by shrinking the grading exponent toward 1 at fixed row
            # count -- the same limit as dt -> 0 in a real convecting wake.
            function graded_gap(ratio)
                nf(irow, icol) = [0.4 * (ratio^(irow - 1) - 1) / max(ratio - 1, 1e-12),
                                  1.0 * (icol - 1), 0.0]
                w = fixture(; nwakerows=3, wrap=false, node_fun=nf)
                w.panel_wake.particle_handoff_active[] = true   # steady state
                pnl._convert_to_particles!(w, shed_body())
                r = w.conversion_diagnostics.records[end]
                # interior panel only: root/tip differ by design, not by grid
                return r.panels[2].kappa_difference, r.residual_abs
            end

            d1, res1 = graded_gap(1.40)
            d2, res2 = graded_gap(1.20)
            d3, res3 = graded_gap(1.10)
            d4, res4 = graded_gap(1.05)

            @test d1 > d2 > d3 > d4          # monotone in the stretch ratio
            @test d4 < 0.4 * d1              # and vanishing as the mesh evens out
            # halving the excess stretch roughly halves the gap (first order)
            @test 1.7 < d3 / d4 < 2.3
            # conservation, meanwhile, is exact at every level
            @test all(<(1e-12), (res1, res2, res3, res4))
        end

        @testset "affine field deposits the exact integrated strength end to end" begin
            # Phase 4 §1: the unit-level affine checks pin the reconstruction;
            # this pins the whole `_convert_to_particles!` chain. The default
            # fixture is the exact affine field muhat = 0.2x + 0.3y + c on the
            # uniform planar sheet, so each panel's area particles must sum to
            # kappa_exact * area, grouped from the particle field itself.
            wake = fixture(; nwakerows=3, wrap=false)
            wake.panel_wake.particle_handoff_active[] = true   # steady state
            pnl._convert_to_particles!(wake, shed_body())
            rec = wake.conversion_diagnostics.records[end]
            ws = wake.conversion_workspace
            pf = wake.pfield

            # Group area particles by spanwise panel: the uniform fixture's
            # panel j spans y in [j-1, j] and subcell centroids are interior.
            sums = [zero(SVector{3,Float64}) for _ in rec.panels]
            for i in 1:pf.np
                ws.classes[i] == pnl.PerimeterLineParticle && continue
                j = clamp(floor(Int, ws.positions[i][2]) + 1, 1, length(sums))
                sums[j] += SVector{3}(pf.particles[FLOWVPM.GAMMA_INDEX, i])
            end
            for (j, p) in enumerate(rec.panels)
                @test sums[j] ≈ p.deposited_strength atol=1e-12
            end

            # Interior panel: exact quadrature-integrated strength.
            kappa_exact = -cross(SVector(0.0, 0.0, 1.0), SVector(0.2, 0.3, 0.0))
            @test rec.panels[2].area ≈ 0.5 atol=1e-12          # 0.5 x 1 panel
            @test sums[2] ≈ kappa_exact * rec.panels[2].area atol=1e-12
            @test rec.residual_abs < 1e-12
        end

        @testset "kappa converges on nonuniform and warped sheets" begin
            # Phase 4 §1: measured convergence of the deposited surface
            # vorticity under panel (row) refinement, against the analytic
            # pointwise kappa of a smooth non-affine field, on (a) a
            # nonuniformly graded planar sheet and (b) a warped sheet. The
            # grading ratio approaches 1 with the scale (a smooth grid map),
            # the refinement limit of a real convecting wake.
            x0 = 0.8
            mu(x) = sin(1.7x)
            dmu(x) = 1.7 * cos(1.7x)

            function kappa_at(scale; f=x -> 0.0)
                r = 1 + 0.75 * scale
                g(i) = (r^(i - 1) - 1) / (r - 1)
                xnode(irow) = x0 + scale * (g(irow) - g(3))
                xcent(irow) = 0.5 * (xnode(irow) + xnode(irow + 1))
                w = fixture(; nwakerows=3, wrap=false, s=10.0,
                    node_fun=(irow, icol) -> [xnode(irow), 1.0 * (icol - 1),
                                              f(xnode(irow))],
                    strength_fun=(irow, icol) -> mu(xcent(irow)))
                w.panel_wake.particle_handoff_active[] = true   # steady state
                pnl._convert_to_particles!(w, shed_body())
                return norm(w.conversion_diagnostics.records[end].panels[2].kappa_conservative)
            end

            warp(x) = 0.4 * sin(2x)
            dwarp(x) = 0.8 * cos(2x)
            cases = (
                ("nonuniform planar", scale -> kappa_at(scale), abs(dmu(x0))),
                ("warped", scale -> kappa_at(scale; f=warp),
                 abs(dmu(x0)) / sqrt(1 + dwarp(x0)^2)),   # |dmu/ds| on z=f(x)
            )
            for (label, k, exact) in cases
                errs = [abs(k(s) - exact) / exact for s in (0.4, 0.2, 0.1)]
                orders = [log2(errs[i] / errs[i+1]) for i in 1:2]
                @test errs[1] > errs[2] > errs[3]     # monotone under refinement
                # first-order face rule (curvature/averaging terms admit a band)
                @test all(o -> 0.6 < o < 2.2, orders)
                @info "016 Phase 4 kappa convergence" label errs orders
            end
        end

        @testset "transfer matches the wake's own filament ledger" begin
            # The conservation test above compares the transaction against its
            # own `expected_total`, so it cannot see a face the deposition code
            # forgets entirely -- which is exactly how the startup edge was lost.
            # This check is independent of the deposition code: every panel
            # ring's perimeter vector sums to zero, so the whole field-relevant
            # content of the panel wake reduces to any retained final filament
            # that the wake actually returns from `get_sources`,
            #
            #     S = sum_j filament_j * (node(N+1,j+1) - node(N+1,j)),
            #
            # and what the particles gain must equal what S loses.
            function S(w)
                pw = w.panel_wake
                tot = SVector(0.0, 0.0, 0.0)
                emits_filament = any(s -> s isa pnl.FilamentWrapper,
                                     pnl.get_sources(pw))
                emits_filament || return tot
                for is in eachindex(pw.nodes)
                    nd = pw.nodes[is]; N = pw.nwakes[]; nc = size(nd, 3) - 1
                    for j in 1:nc
                        f = pnl._final_filament_strength(pw, is, N, j)
                        e = SVector{3}(nd[1,N+1,j+1], nd[2,N+1,j+1], nd[3,N+1,j+1]) -
                            SVector{3}(nd[1,N+1,j],   nd[2,N+1,j],   nd[3,N+1,j])
                        tot += f * e
                    end
                end
                return tot
            end
            gained(w, n0) = SVector{3}(
                sum(w.pfield.particles[FLOWVPM.GAMMA_INDEX[k], n0+1:w.pfield.np])
                for k in 1:3)
            # `shed_wake!` copies node row 1 -> 2 and leaves row 1 stale; the
            # real loop re-places it via update_TE!/propagate! at the start of
            # the next step. Emulate that so the new row-1 panel is not
            # degenerate. Only node row 1 moves, and a complete ring contributes
            # nothing to S, so this cannot perturb the invariant.
            function refresh_te!(w, dx)
                for nd in w.panel_wake.nodes
                    nd[:, 1, :] .= nd[:, 2, :]
                    nd[1, 1, :] .-= dx
                end
            end

            for attr in (:upstream, :downstream, :split), nwakerows in (1, 2, 3),
                    wrap in (false, true)
                wake = stretched(; nwakerows=nwakerows, wrap=wrap, attribution=attr)
                body = shed_body()
                for step in 1:2      # 1 = startup branch, 2 = steady state
                    step > 1 && refresh_te!(wake, 0.35)
                    s_before = S(wake)
                    n0 = wake.pfield.np
                    pnl.shed_wake!(wake, body)
                    @test (s_before - S(wake)) ≈ gained(wake, n0) atol=1e-12
                end
            end
        end

        @testset "the startup edge is deposited, not deleted" begin
            # Before the first handoff the aft-most face is the sheet's physical
            # trailing boundary -- the starting vortex. There is no earlier
            # conversion to have taken the upstream share, so it must be
            # deposited whole regardless of attribution.
            for attr in (:upstream, :downstream, :split)
                wake = fixture(; nwakerows=3, wrap=false, attribution=attr)
                pw = wake.panel_wake
                N = pw.nwakes[]; nodes = pw.nodes[1]; st = pw.strength[1]
                nc = size(nodes, 3) - 1
                node(r, c) = SVector{3}(nodes[1,r,c], nodes[2,r,c], nodes[3,r,c])

                # Net circulation sitting uncancelled on the aft face right now.
                startup = SVector(0.0, 0.0, 0.0)
                for j in 1:nc
                    fil = pnl._final_filament_strength(pw, 1, N, j)
                    startup += (st[1,N,j] + fil) * (node(N+1,j+1) - node(N+1,j))
                end
                @test norm(startup) > 0.1          # the fixture has a real one

                pnl._convert_to_particles!(wake, shed_body())
                rec = wake.conversion_diagnostics.records[end]
                @test rec.startup_edge_deposited
                @test rec.attribution === attr
                @test rec.deposited_total ≈ expected_transfer(wake;
                    alpha=Float64(pnl._ATTRIBUTION_ALPHA[attr]),
                    startup=true) atol=1e-12
                # the aft face's full content is inside what was deposited
                @test sum(p.downstream_face for p in rec.panels) ≈ startup atol=1e-12

                # A second conversion is no longer a startup, so it takes only
                # its attributed share.
                pnl._convert_to_particles!(wake, shed_body())
                rec2 = wake.conversion_diagnostics.records[end]
                @test !rec2.startup_edge_deposited
            end
        end

        @testset "the retained filament matches the attribution" begin
            # Whatever fraction of the interface the conversion did not deposit
            # must stay on the panel side, or the next step double-counts or
            # leaks. `:downstream` must reproduce the legacy filament exactly.
            for attr in (:upstream, :downstream, :split)
                wake = stretched(; nwakerows=3, wrap=false, attribution=attr)
                pnl.shed_wake!(wake, shed_body())
                pw = wake.panel_wake
                N = pw.nwakes[]; st = pw.strength[1]
                alpha = pnl._ATTRIBUTION_ALPHA[attr]
                @test pw.particle_handoff_active[]
                @test pw.particle_handoff_weight[] == alpha
                for j in 1:size(st, 3)
                    expected = -(alpha * st[1, N, j] + (1 - alpha) * st[1, N+1, j])
                    @test pnl._final_filament_strength(pw, 1, N, j) ≈ expected
                end
            end

            # `:downstream` is bit-identical to the legacy unsteady filament.
            dw = stretched(; nwakerows=3, wrap=false, attribution=:downstream)
            pnl.shed_wake!(dw, shed_body())
            pw = dw.panel_wake; N = pw.nwakes[]
            for j in 1:size(pw.strength[1], 3)
                @test pnl._final_filament_strength(pw, 1, N, j) ==
                      -pw.strength[1][1, N+1, j]
            end
        end

        @testset "split attribution is second-order in row spacing" begin
            # Depositing half of each streamwise face makes the difference
            # centered, so the leading grid-stretching error cancels. This is
            # the justification for keeping `:split` available at all.
            # A strictly LINEAR field, so the only error left is the one under
            # test: the streamwise face rule's sensitivity to row-to-row
            # stretching. A curved field adds an opposite-signed truncation term
            # that cancels against it at some resolution and destroys the
            # measured order.
            mu_f(x) = 0.8 + 1.3x
            dmu_f(x) = 1.3

            # A genuine refinement of a graded mesh: the row-to-row ratio must
            # approach 1 as the spacing shrinks (a smooth grid-generating map).
            # Holding the ratio fixed while shrinking the scale is NOT a
            # refinement -- the leading stretching error is then O(1).
            function kappa_error(attr, scale)
                ratio = 1 + 0.5 * scale
                xnode(irow) = scale * (ratio^(irow - 1) - 1) / (ratio - 1)
                xcent(irow) = 0.5 * (xnode(irow) + xnode(irow + 1))
                # Spanwise-constant strength isolates the streamwise face rule.
                w = fixture(; nwakerows=3, wrap=false, attribution=attr, s=10.0,
                    node_fun=(irow, icol) -> [xnode(irow), 1.0 * (icol - 1), 0.0],
                    strength_fun=(irow, icol) -> mu_f(xcent(irow)))
                # Emulate a steady-state conversion (not the startup one), so the
                # attribution's own weights are what act.
                w.panel_wake.particle_handoff_active[] = true
                pnl._convert_to_particles!(w, shed_body())
                rec = w.conversion_diagnostics.records[end]
                N = w.panel_wake.nwakes[]
                exact = -dmu_f(xcent(N))
                # kappa is uniform per panel; take an interior column.
                return abs(rec.panels[2].kappa_conservative[2] - exact) / abs(exact)
            end

            for attr in (:upstream, :split)
                e1 = kappa_error(attr, 0.4)
                e2 = kappa_error(attr, 0.2)
                order = log2(e1 / e2)
                if attr === :upstream
                    @test 0.7 < order < 1.4        # first order
                else
                    @test 1.6 < order < 2.5        # second order
                    @test e2 < 0.2 * kappa_error(:upstream, 0.2)
                end
            end
        end

        @testset "attribution is validated" begin
            @test pnl.SurfaceVorticityConversion(0.2).attribution === :upstream
            for a in (:upstream, :downstream, :split)
                @test pnl.SurfaceVorticityConversion(0.2; attribution=a).attribution === a
            end
            @test_throws ArgumentError pnl.SurfaceVorticityConversion(0.2; attribution=:middle)
            @test_throws ArgumentError pnl.SurfaceVorticityConversion(0.2; attribution=:none)
        end

        @testset "only true root and tip lines survive" begin
            open_wake = fixture(; nwakerows=3, wrap=false)
            pnl._convert_to_particles!(open_wake, shed_body())
            orec = open_wake.conversion_diagnostics.records[end]
            st = open_wake.panel_wake.strength[1]
            n_cols = size(open_wake.panel_wake.nodes[1], 3) - 1

            # The closure lines are the LEGACY filaments, unchanged: in
            # divergence form these are the two streamwise faces with no
            # neighbour to share with, so the stored cell value is correct and
            # no edge extrapolation is needed.
            @test orec.root_strength == st[1, 3, 1]
            @test orec.tip_strength == -st[1, 3, n_cols]
            @test orec.root_owned && orec.tip_owned
            @test orec.n_perimeter > 0
            @test orec.n_roottip > 0 && orec.n_interior > 0

            wrap_wake = fixture(; nwakerows=3, wrap=true)
            pnl._convert_to_particles!(wrap_wake, shed_body())
            wrec = wrap_wake.conversion_diagnostics.records[end]
            @test wrec.n_perimeter == 0
            @test !wrec.root_owned && !wrec.tip_owned
            @test wrec.root_strength == 0 && wrec.tip_strength == 0
            @test wrec.n_roottip == 0
            @test count(==(pnl.PerimeterLineParticle),
                        wrap_wake.conversion_workspace.classes) == 0
        end

        @testset "no upstream geometry is read at nwakerows == 1" begin
            # The whole reason the divergence form was adopted: at N == 1 the
            # upstream row's nodes do not exist yet (shed_wake! writes only
            # strength; update_TE! writes the nodes next step). Only the
            # upstream *strength* is consulted, so moving the body's trailing
            # edge or Das must not change a single deposited particle.
            base = fixture(; nwakerows=1, wrap=false)
            moved = fixture(; nwakerows=1, wrap=false)
            b1 = shed_body()
            b2 = shed_body()
            b2.nodes .+= 7.5                 # translate the whole body
            b2.Das[1] .*= 3.0                # and stretch the wake offset
            pnl._convert_to_particles!(base, b1)
            pnl._convert_to_particles!(moved, b2)

            @test base.pfield.np == moved.pfield.np
            @test base.pfield.particles[1:3, 1:base.pfield.np] ==
                  moved.pfield.particles[1:3, 1:moved.pfield.np]
            @test base.pfield.particles[FLOWVPM.GAMMA_INDEX, 1:base.pfield.np] ==
                  moved.pfield.particles[FLOWVPM.GAMMA_INDEX, 1:moved.pfield.np]

            # The upstream *strength*, by contrast, must matter.
            b3 = shed_body()
            b3.strength[:, end] .= collect(1.0:size(b3.strength, 1))
            strong = fixture(; nwakerows=1, wrap=false)
            pnl._convert_to_particles!(strong, b3)
            @test strong.conversion_diagnostics.records[end].expected_handoff !=
                  base.conversion_diagnostics.records[end].expected_handoff

            # ...and it must refuse to guess when the body is absent.
            orphan = fixture(; nwakerows=1, wrap=false)
            @test_throws pnl.WakeConversionStateError pnl._convert_to_particles!(orphan)
            @test orphan.pfield.np == 0

            # With an upstream row present, the body is irrelevant entirely.
            for nwakerows in (2, 3)
                a = fixture(; nwakerows=nwakerows, wrap=false)
                b = fixture(; nwakerows=nwakerows, wrap=false)
                bb = shed_body(); bb.strength .= 3.7; bb.Das[1] .*= 2.5
                pnl._convert_to_particles!(a, shed_body())
                pnl._convert_to_particles!(b, bb)
                @test a.pfield.particles[1:3, 1:a.pfield.np] ==
                      b.pfield.particles[1:3, 1:b.pfield.np]
                @test a.pfield.particles[FLOWVPM.GAMMA_INDEX, 1:a.pfield.np] ==
                      b.pfield.particles[FLOWVPM.GAMMA_INDEX, 1:b.pfield.np]
            end
        end

        @testset "active source count stays N; the stale row has no influence" begin
            # Phase 4 §2. Under the approved Option B there is no ghost
            # storage: the "ghost" is the pre-shift final active row, which IS
            # a legitimate source. What must be invisible is the stale strength
            # row N+1 -- it backs only the legacy unsteady filament, so once
            # the smooth handoff retargets the filament onto row N
            # (`:upstream`), nothing may read it.
            targets = (SVector(0.7, 1.3, 0.9), SVector(-0.5, 0.4, -1.1),
                       SVector(2.0, 2.5, 0.3))
            function probe(w)
                probes = FastMultipole.ProbeSystem(length(targets), Float64)
                for i in eachindex(targets)
                    probes.position[i] = targets[i]
                    probes.scalar_potential[i] = 0.0
                    probes.gradient[i] = zero(SVector{3,Float64})
                    probes.hessian[i] = zero(FastMultipole.SMatrix{3,3,Float64,9})
                end
                pnl.influence!((probes,), pnl.get_sources(w.panel_wake),
                    pnl.DirectBackend(); scalar_potential=false, gradient=true,
                    hessian=(false,))
                return copy(probes.gradient)
            end

            wake = stretched(; nwakerows=3, wrap=false)
            pw = wake.panel_wake
            N = pw.nwakes[]
            ncols = size(pw.strength[1], 3)

            # N active panel sources despite N+1 stored strength rows.
            @test size(pw.strength[1], 2) == N + 1
            @test FastMultipole.get_n_bodies(pw) == N * ncols

            # R4, pinned: a static fixture emits NO retained filament until
            # `overflowed[]` is set; production sets it in simulate!/warm start.
            fil = pnl.FilamentWrapper(pw)
            @test FastMultipole.get_n_bodies(fil) == 0
            pw.overflowed[] = true
            @test FastMultipole.get_n_bodies(fil) == ncols

            pnl._convert_to_particles!(wake, shed_body())
            @test FastMultipole.get_n_bodies(pw) == N * ncols  # unchanged

            stale = stretched(; nwakerows=3, wrap=false)
            active = stretched(; nwakerows=3, wrap=false)
            pnl._convert_to_particles!(stale, shed_body())
            pnl._convert_to_particles!(active, shed_body())
            for w in (stale, active)
                w.panel_wake.overflowed[] = true
            end
            stale.panel_wake.strength[1][1, N + 1, :] .+= 100.0
            active.panel_wake.strength[1][1, N, :] .+= 100.0

            base_field = probe(wake)
            @test probe(stale) == base_field          # row N+1 is invisible
            @test !(probe(active) ≈ base_field)       # the metric has teeth
        end

        @testset "legacy and smooth differ only by the handoff share" begin
            # Both strategies redistribute the same filament set; they differ in
            # WHICH streamwise-neighbour edge each conversion owns. The legacy
            # unsteady filament takes the downstream jump (against the stale row
            # past the buffer); the smooth path takes the upstream handoff,
            # because `particle_handoff_active` moves the retained filament.
            # Over a run these telescope to the same thing; per conversion they
            # differ by exactly this, on any geometry.
            for nwakerows in (2, 3), wrap in (false, true)
                legacy = make_conversion_fixture(; nwakerows=nwakerows, wrap=wrap,
                                                 max_particles=5000)
                pnl._convert_to_particles!(legacy)
                sm = fixture(; nwakerows=nwakerows, wrap=wrap)
                pnl._convert_to_particles!(sm, shed_body())

                pw = sm.panel_wake
                nodes = pw.nodes[1]; st = pw.strength[1]
                N = nwakerows; nc = size(nodes, 3) - 1
                node(r, c) = SVector{3}(nodes[1,r,c], nodes[2,r,c], nodes[3,r,c])
                # Legacy deposits the aft face; on a startup conversion the
                # smooth path deposits it too, so the two differ by exactly the
                # attributed share of the upstream face.
                delta = SVector(0.0, 0.0, 0.0)
                for j in 1:nc
                    delta += (st[1,N-1,j] - st[1,N,j]) * (node(N,j+1) - node(N,j))
                end
                @test netG(sm) - netG(legacy) ≈ delta atol=1e-12
            end
        end

        @testset "constant strength deposits nothing" begin
            # Constant muhat: every face jump and every handoff vanishes, so
            # kappa is exactly zero and the whole area deposit is elided.
            wake = fixture(; nwakerows=3, wrap=false,
                           strength_fun=(irow, icol) -> 0.0)
            pnl._convert_to_particles!(wake, shed_body())
            rec = wake.conversion_diagnostics.records[end]
            @test wake.pfield.np == 0
            @test rec.n_elided > 0
            @test iszero(rec.total_area_strength)
            @test rec.residual_abs < 1e-14

            # A nonzero constant still deposits no area vorticity, but its
            # root/tip closure is a genuine closed constant-strength loop.
            cw = fixture(; nwakerows=3, wrap=false,
                         strength_fun=(irow, icol) -> 0.7)
            pnl._convert_to_particles!(cw, shed_body())
            crec = cw.conversion_diagnostics.records[end]
            @test iszero(crec.total_area_strength)
            @test crec.root_strength ≈ 0.7
            @test crec.tip_strength ≈ -0.7
            @test crec.residual_abs < 1e-14
        end

        @testset "handoff flag flips and retargets the final filament" begin
            wake = fixture(; nwakerows=3, wrap=false)
            pw = wake.panel_wake
            @test pw.particle_handoff_active[] == false
            @test pnl._final_filament_strength(pw, 1, 3, 1) == -pw.strength[1][1, 4, 1]

            pnl._convert_to_particles!(wake, shed_body())
            @test pw.particle_handoff_active[] == true
            @test pnl._final_filament_strength(pw, 1, 3, 1) == -pw.strength[1][1, 3, 1]

            legacy = make_conversion_fixture(; nwakerows=3, wrap=false)
            pnl._convert_to_particles!(legacy)
            @test legacy.panel_wake.particle_handoff_active[] == false
            @test legacy.conversion_count[] == 0
        end

        @testset "repeated conversions keep the ledger" begin
            wake = stretched(; nwakerows=3, wrap=false)
            body = shed_body()
            pnl._convert_to_particles!(wake, body)
            np1 = wake.pfield.np
            first_record = wake.conversion_diagnostics.records[1]
            pnl._convert_to_particles!(wake, body)

            diag = wake.conversion_diagnostics
            @test length(diag.records) == 2
            @test diag.records[1] === first_record
            @test diag.records[2].ordinal == 2
            @test wake.conversion_count[] == 2
            @test wake.pfield.np == 2 * np1
            @test diag.total_particles == wake.pfield.np
            @test norm(diag.total_residual) < 1e-12
            @test diag.total_deposited ≈ diag.total_expected atol=1e-12
            @test diag.total_area ≈ 2 * first_record.total_area
        end

        @testset "failures are transactional" begin
            state_snapshot(w) = (copy(w.panel_wake.nodes[1]),
                                 copy(w.panel_wake.strength[1]),
                                 w.panel_wake.nwakes[],
                                 w.panel_wake.particle_handoff_active[],
                                 w.conversion_count[],
                                 w.pfield.np,
                                 copy(w.pfield.particles[:, 1:w.pfield.np]))

            small = make_conversion_fixture(; nwakerows=3, wrap=false,
                max_particles=4, conversion=pnl.SurfaceVorticityConversion(0.05))
            before = state_snapshot(small)
            err = try
                pnl._convert_to_particles!(small, shed_body()); nothing
            catch e; e end
            @test err isa pnl.PanelParticleCapacityError
            @test err.available == 4
            @test err.requested > err.available
            @test state_snapshot(small) == before
            @test isempty(small.conversion_diagnostics.records)

            bad = fixture(; nwakerows=3, wrap=false)
            bad.panel_wake.nodes[1][2, 3, 2] = NaN
            before = state_snapshot(bad)
            @test_throws pnl.WakeGeometryError pnl._convert_to_particles!(bad, shed_body())
            @test isequal(state_snapshot(bad), before)   # the poisoned NaN is preserved

            early = fixture(; nwakerows=3, wrap=false)
            early.panel_wake.nwakes[] = 2
            before = state_snapshot(early)
            @test_throws pnl.WakeConversionStateError pnl._convert_to_particles!(early, shed_body())
            @test state_snapshot(early) == before

            diagnostic = pnl.PanelParticleWake(shed_body();
                nwakerows=2, max_particles=100,
                conversion=pnl.SurfaceVorticityConversion(0.2; diagnose_nearfield=true))
            @test diagnostic.conversion.diagnose_nearfield
        end

        @testset "diagnostics record the mandated fields" begin
            wake = stretched(; nwakerows=3, wrap=false)
            pnl._convert_to_particles!(wake, shed_body())
            rec = wake.conversion_diagnostics.records[end]

            @test rec.capacity_requested == wake.pfield.np
            @test rec.capacity_available >= rec.capacity_requested
            @test rec.handoff_active_before == false
            @test rec.handoff_active_after == true
            @test rec.step_id == wake.panel_wake.live_step_id[]
            @test sum(p.area for p in rec.panels) ≈ rec.total_area
            @test sum(p.handoff for p in rec.panels) ≈ rec.expected_handoff

            for p in rec.panels
                @test p.i_surf == 1
                @test p.rank == 2                       # both diagnostic legs observable
                @test isfinite(p.condition)
                @test all(p.observable)
                @test p.singular_values[1] >= p.singular_values[2] > p.rank_threshold
                @test p.n_xi >= 1 && p.n_eta >= 1
                @test p.area > 0
                @test p.geometry_scale > 0
                @test isfinite(p.min_jacobian) && p.min_jacobian > 0
                @test p.n_requested == p.n_xi * p.n_eta
                @test p.n_elided == 0
                @test p.deposited_strength ≈ p.kappa_conservative * p.area atol=1e-12
                @test p.kappa_difference ≈
                      norm(p.kappa_conservative - p.kappa_reconstruction)
            end
            @test (rec.panels[1].class, rec.panels[end].class) ==
                  (pnl.RootTipSurfaceParticle, pnl.RootTipSurfaceParticle)
            @test rec.panels[2].class == pnl.InteriorSurfaceParticle
            @test rec.nearfield === nothing
        end

        @testset "optional near-field gradient summaries are panel-only" begin
            conversion = pnl.SurfaceVorticityConversion(0.2;
                diagnose_nearfield=true)
            diagnosed = make_conversion_fixture(; nwakerows=3, wrap=false,
                max_particles=40000, conversion)
            reference = make_conversion_fixture(; nwakerows=3, wrap=false,
                max_particles=40000, conversion)
            contaminated = make_conversion_fixture(; nwakerows=3, wrap=false,
                max_particles=40000, conversion)
            FLOWVPM.add_particle(contaminated.pfield, [100.0, -50.0, 20.0],
                [1e8, -2e8, 3e8], 0.01)

            pnl._convert_to_particles!(diagnosed, shed_body())
            pnl._convert_to_particles!(contaminated, shed_body())
            diag = diagnosed.conversion_diagnostics.records[end].nearfield
            cdiag = contaminated.conversion_diagnostics.records[end].nearfield
            @test diag isa pnl.SurfaceVorticityNearFieldDiagnostics
            @test diag == cdiag

            # Independent direct evaluation against the still-unconverted
            # reference panel wake. Existing particles are deliberately absent
            # from the source tuple.
            ws = diagnosed.conversion_workspace
            probes = FastMultipole.ProbeSystem(length(ws.positions), Float64)
            for i in eachindex(ws.positions)
                probes.position[i] = ws.positions[i]
                probes.scalar_potential[i] = 0.0
                probes.gradient[i] = zero(SVector{3,Float64})
                probes.hessian[i] = zero(FastMultipole.SMatrix{3,3,Float64,9})
            end
            pnl.influence!((probes,), pnl.get_sources(reference.panel_wake),
                pnl.DirectBackend(); scalar_potential=false, gradient=false,
                hessian=(true,))
            vals(class) = sort([norm(probes.hessian[i]) for i in eachindex(ws.classes)
                                if ws.classes[i] == class])
            function expected_summary(v)
                n = length(v); avg = sum(v)/n
                med = isodd(n) ? v[(n+1)÷2] : (v[n÷2] + v[n÷2+1])/2
                return (n, first(v), last(v), avg, sqrt(sum(abs2,v)/n), med,
                        v[ceil(Int, 0.95n)])
            end
            for (summary, class) in (
                    (diag.interior, pnl.InteriorSurfaceParticle),
                    (diag.roottip, pnl.RootTipSurfaceParticle),
                    (diag.perimeter, pnl.PerimeterLineParticle))
                expected = expected_summary(vals(class))
                @test summary.count == expected[1]
                @test collect((summary.minimum, summary.maximum, summary.mean,
                               summary.rms, summary.median, summary.p95)) ≈
                      collect(expected[2:end]) rtol=1e-13 atol=1e-13
            end

            # A failure in direct diagnostic evaluation is still preflight:
            # no particle or handoff mutation can escape.
            snapshot(w) = (copy(w.panel_wake.nodes[1]),
                copy(w.panel_wake.strength[1]), w.panel_wake.nwakes[],
                w.panel_wake.particle_handoff_active[], w.conversion_count[],
                w.pfield.np, copy(w.pfield.particles[:, 1:w.pfield.np]))
            bad = make_conversion_fixture(; nwakerows=3, wrap=false,
                max_particles=40000, conversion)
            bad.panel_wake.strength[1][1, end, 2] = NaN
            before = snapshot(bad)
            @test_throws pnl.WakeGeometryError pnl._convert_to_particles!(bad, shed_body())
            @test isequal(snapshot(bad), before)
            @test isempty(bad.conversion_diagnostics.records)
        end

        @testset "subdivision follows sigma / overlap" begin
            # Refining the target spacing changes only how finely the same
            # vorticity is sampled: area, integrated strength, and the
            # conservation invariant are all untouched.
            coarse = stretched(; nwakerows=3, wrap=false, s=0.5)
            fine = stretched(; nwakerows=3, wrap=false, s=0.125)
            pnl._convert_to_particles!(coarse, shed_body())
            pnl._convert_to_particles!(fine, shed_body())

            @test fine.pfield.np > coarse.pfield.np
            crec = coarse.conversion_diagnostics.records[end]
            frec = fine.conversion_diagnostics.records[end]
            # The 2x2 Gauss area of a *warped* panel is only second-order
            # accurate, so refining the subdivision nudges the area (and hence
            # kappa = V / A) by a converging amount -- rtol, not atol.
            @test frec.total_area ≈ crec.total_area rtol=1e-7
            @test frec.total_area_strength ≈ crec.total_area_strength rtol=1e-7
            # the face content itself is subdivision-independent and exact
            @test frec.expected_total ≈ crec.expected_total atol=1e-12
            @test frec.residual_abs < 1e-12
            for (fp, cp) in zip(frec.panels, crec.panels)
                @test fp.kappa_conservative ≈ cp.kappa_conservative rtol=1e-7
            end
        end

        @testset "resolution floor: one particle at the bilinear centroid" begin
            # Phase 4 §1: for a nonzero full-rank field, when both subdivision
            # counts hit the floor (1, 1) the single area particle per panel
            # sits at the bilinear panel centroid -- not just "one particle".
            wake = fixture(; nwakerows=3, wrap=false, s=1e6)
            pnl._convert_to_particles!(wake, shed_body())
            rec = wake.conversion_diagnostics.records[end]
            ws = wake.conversion_workspace
            pw = wake.panel_wake
            N = pw.nwakes[]
            nodes = pw.nodes[1]
            ncols = size(pw.strength[1], 3)

            @test all(p.n_xi == 1 && p.n_eta == 1 for p in rec.panels)
            area_ids = [i for i in 1:wake.pfield.np
                        if ws.classes[i] != pnl.PerimeterLineParticle]
            @test length(area_ids) == ncols            # one per panel

            # Bilinear map at (1/2, 1/2) is the mean of the four corners.
            centroids = [0.25 * (SVector{3}(nodes[:, N, j]) +
                                 SVector{3}(nodes[:, N + 1, j]) +
                                 SVector{3}(nodes[:, N + 1, j + 1]) +
                                 SVector{3}(nodes[:, N, j + 1]))
                         for j in 1:ncols]
            positions = sort([SVector{3}(wake.pfield.particles[1:3, i])
                              for i in area_ids]; by=p -> p[2])
            for (x, c) in zip(positions, centroids)
                @test x ≈ c atol=1e-13
            end
        end


        @testset "static sheet-to-hybrid near-field acceptance gate" begin
            # No propagation or shedding: compare the prescribed full sheet to
            # retained rows + final filament + staged particles after removing
            # the outgoing row analytically.
            function probe_velocity(wake, positions; hybrid=false)
                probes = FastMultipole.ProbeSystem(length(positions), Float64)
                for i in eachindex(positions)
                    probes.position[i] = positions[i]
                    probes.scalar_potential[i] = 0.0
                    probes.gradient[i] = zero(SVector{3,Float64})
                    probes.hessian[i] = zero(FastMultipole.SMatrix{3,3,Float64,9})
                end
                sources = hybrid ?
                    (pnl.get_sources(wake.panel_wake)..., wake.pfield) :
                    pnl.get_sources(wake.panel_wake)
                pnl.influence!((probes,), sources, pnl.DirectBackend();
                    scalar_potential=false, gradient=true, hessian=(false,))
                return copy(probes.gradient)
            end

            function static_errors(attribution)
                w = stretched(; nwakerows=3, wrap=false, s=0.2,
                    attribution=attribution)
                pw = w.panel_wake; N = pw.nwakes[]; nodes = pw.nodes[1]
                positions = SVector{3,Float64}[]
                family = Symbol[]; ratio_of = Float64[]
                # handoff, interior, and true root/tip families on both sides
                # of the sheet at the required d/sigma ladder.
                for ratio in (0.25,0.5,1.0,2.0,4.0,8.0)
                    for (label,j,xi) in ((:root,1,0.5),(:interior,2,0.5),
                                         (:tip,3,0.5),(:handoff,2,0.0))
                        v1 = SVector{3}(nodes[:,N,j]); v2 = SVector{3}(nodes[:,N+1,j])
                        v3 = SVector{3}(nodes[:,N+1,j+1]); v4 = SVector{3}(nodes[:,N,j+1])
                        x = xi == 0 ? 0.5*(v1+v4) : 0.25*(v1+v2+v3+v4)
                        n = cross(v2-v1,v4-v1); n /= norm(n)
                        for side in (-1,1)
                            push!(positions,x + side*ratio*0.2*n)
                            push!(family,label); push!(ratio_of,ratio)
                        end
                    end
                end
                sheet = probe_velocity(w,positions)
                pnl._convert_to_particles!(w,shed_body())
                pw.nwakes[] = N-1
                hybrid = probe_velocity(w,positions; hybrid=true)
                metrics = Dict{Tuple{Symbol,Float64},Tuple{Float64,Float64}}()
                for label in (:handoff,:interior,:root,:tip),
                    ratio in (0.25,0.5,1.0,2.0,4.0,8.0)
                    ids = findall(i -> family[i] == label && ratio_of[i] == ratio,
                                  eachindex(positions))
                    delta = [norm(hybrid[i]-sheet[i]) for i in ids]
                    ref_rms = sqrt(sum(norm(sheet[i])^2 for i in ids)/length(ids))
                    metrics[(label,ratio)] =
                        (sqrt(sum(abs2,delta)/length(delta))/ref_rms,
                         maximum(delta)/maximum(norm(sheet[i]) for i in ids))
                end
                G = sum(SVector{3}(w.pfield.particles[FLOWVPM.GAMMA_INDEX,i])
                        for i in 1:w.pfield.np; init=zero(SVector{3,Float64}))
                I = sum(cross(SVector{3}(w.pfield.particles[FLOWVPM.X_INDEX,i]),
                              SVector{3}(w.pfield.particles[FLOWVPM.GAMMA_INDEX,i]))
                        for i in 1:w.pfield.np; init=zero(SVector{3,Float64}))
                return metrics,G,I,w.pfield.np
            end

            upstream,Gup,Iup,npup = static_errors(:upstream)
            split,Gsp,Isp,npsp = static_errors(:split)
            @test all(isfinite, Iterators.flatten(values(upstream)))
            @test all(isfinite, Iterators.flatten(values(split)))
            @test npup == npsp > 0
            @test all(isfinite,Gup) && all(isfinite,Iup)
            @test all(isfinite,Gsp) && all(isfinite,Isp)
            # The binding decision rule is conjunctive. Split fails its first
            # clause on this representative coarse smooth stretched sheet: it
            # does not lower handoff RMS at every standoff. Keep :upstream.
            @test !all(split[(:handoff,d)][1] < upstream[(:handoff,d)][1]
                       for d in (0.25,0.5,1.0,2.0,4.0,8.0))
            @test pnl.SurfaceVorticityConversion(0.2).attribution == :upstream

            function adjacent_field_error(radius)
                s = 0.2
                body = shed_body()
                hybrid = stretched(; nwakerows=3, wrap=false, s=s,
                    attribution=:upstream)
                hpw = hybrid.panel_wake
                # Prescribed sheet has a physical downstream boundary: no
                # older-row cancellation exists beyond the initial final row.
                for st in hpw.strength
                    st[:,end,:] .= 0
                end
                pure = pnl.PanelWake(body; nwakerows=5,
                    core_size=hpw.core_size)
                for isurf in eachindex(pure.nodes)
                    pure.nodes[isurf][:,1:4,:] .= hpw.nodes[isurf]
                    pure.strength[isurf][:,1:4,:] .= hpw.strength[isurf]
                end
                pure.nwakes[] = 3

                # Two adjacent topology transactions, still with no propagation
                # or solve. The larger pure buffer retains both rows that the
                # hybrid converts to particles.
                pnl.shed_wake!(hybrid,body); pnl.shed_wake!(pure,body)
                for wake in (hpw,pure)
                    for nd in wake.nodes
                        nd[:,1,:] .= nd[:,2,:]
                        nd[1,1,:] .-= 0.35
                    end
                end
                pnl.shed_wake!(hybrid,body); pnl.shed_wake!(pure,body)
                # `shed_wake!` leaves row 1 duplicated until the following
                # update_TE!; place that row for the static field snapshot.
                for wake in (hpw,pure)
                    for nd in wake.nodes
                        nd[:,1,:] .= nd[:,2,:]
                        nd[1,1,:] .-= 0.35
                    end
                end
                @test length(hybrid.conversion_diagnostics.records) == 2
                @test all(c != pnl.PerimeterLineParticle ||
                          abs(hybrid.conversion_workspace.positions[i][2] -
                              extrema(hpw.nodes[1][2,:,:])[1]) < 1e-12 ||
                          abs(hybrid.conversion_workspace.positions[i][2] -
                              extrema(hpw.nodes[1][2,:,:])[2]) < 1e-12
                          for (i,c) in enumerate(hybrid.conversion_workspace.classes))

                center = SVector(1.5,1.5,0.0)
                directions = (SVector(1.0,0.2,0.4), SVector(-0.3,1.0,0.5),
                              SVector(0.2,-0.4,1.0), SVector(-0.8,-0.2,0.6))
                positions = [center + radius*d/norm(d) for d in directions]
                probes = FastMultipole.ProbeSystem(length(positions),Float64)
                function evaluate(sources)
                    for i in eachindex(positions)
                        probes.position[i]=positions[i]
                        probes.scalar_potential[i]=0.0
                        probes.gradient[i]=zero(SVector{3,Float64})
                        probes.hessian[i]=zero(FastMultipole.SMatrix{3,3,Float64,9})
                    end
                    pnl.influence!((probes,),sources,pnl.DirectBackend();
                        scalar_potential=false,gradient=true,hessian=(false,))
                    return copy(probes.gradient)
                end
                sheet=evaluate(pnl.get_sources(pure))
                field=evaluate((pnl.get_sources(hpw)...,hybrid.pfield))
                return sqrt(sum(norm(field[i]-sheet[i])^2 for i in eachindex(field))/
                            length(field))
            end
            adjacent_errors = [adjacent_field_error(r) for r in (10.0,20.0,40.0)]
            @test adjacent_errors[3] < adjacent_errors[2] < adjacent_errors[1]
        end
    end

    @testset "PanelWake first row shedding velocity option" begin
        body = make_plate_vortex_body()
        dt = 0.25
        uinf = [1.0, 2.0, 3.0]
        induced_first = [10.0, 20.0, 30.0]
        induced_downstream = [100.0, 200.0, 300.0]

        wake_total = pnl.PanelWake(body; nwakerows=2)
        wake_total.nwakes[] = 1
        pnl.apply_freestream!(wake_total, uinf)
        wake_total.velocity[1][:, 1, :] .+= induced_first
        wake_total.velocity[1][:, 2, :] .+= induced_downstream
        pnl.propagate!(wake_total, dt)

        @test wake_total.nodes[1][:, 1, 1] ≈ dt .* (uinf .+ induced_first)
        @test wake_total.nodes[1][:, 2, 1] ≈ dt .* (uinf .+ induced_downstream)

        wake_freestream = pnl.PanelWake(body; nwakerows=2,
            shed_with_induced_velocity=false)
        wake_freestream.nwakes[] = 1
        pnl.apply_freestream!(wake_freestream, uinf)
        wake_freestream.velocity[1][:, 1, :] .+= induced_first
        wake_freestream.velocity[1][:, 2, :] .+= induced_downstream
        pnl.propagate!(wake_freestream, dt)

        @test wake_freestream.nodes[1][:, 1, 1] ≈ dt .* uinf
        @test wake_freestream.nodes[1][:, 2, 1] ≈ dt .* (uinf .+ induced_downstream)
    end

    @testset "PanelWake freestream convection mode" begin
        body = make_plate_vortex_body()
        dt = 0.25
        uinf = [1.0, 2.0, 3.0]
        induced_first = [10.0, 20.0, 30.0]
        induced_downstream = [100.0, 200.0, 300.0]
        nsteps = 4

        # default preserves legacy behavior
        @test !pnl.PanelWake(body; nwakerows=2).freestream_convection

        wake = pnl.PanelWake(body; nwakerows=nsteps+1,
            freestream_convection=true)
        wake.nwakes[] = 1
        # seed a non-degenerate sheet so the translation test has content
        for i_surf in eachindex(wake.nodes)
            for icol in axes(wake.nodes[i_surf], 3),
                    irow in axes(wake.nodes[i_surf], 2)
                wake.nodes[i_surf][:, irow, icol] .=
                    (0.3 * irow, 0.7 * icol, 0.11 * irow * icol)
            end
        end
        initial = [copy(n) for n in wake.nodes]

        for _ in 1:nsteps
            # induced velocity is present but must be ignored entirely
            pnl.reset!(wake)
            pnl.apply_freestream!(wake, uinf)
            wake.velocity[1][:, 1, :] .+= induced_first
            wake.velocity[1][:, 2, :] .+= induced_downstream
            pnl.propagate!(wake, dt)
        end

        for i_surf in eachindex(wake.nodes)
            nrows = wake.nwakes[] + 1
            moved = view(wake.nodes[i_surf], :, 1:nrows, :)
            expected = view(initial[i_surf], :, 1:nrows, :) .+
                nsteps .* dt .* reshape(uinf, 3, 1, 1)
            @test moved ≈ expected atol=1e-14
            # untouched rows beyond nwakes must not move
            @test wake.nodes[i_surf][:, nrows+1:end, :] ==
                initial[i_surf][:, nrows+1:end, :]
        end

        # every node translated by the same vector => the sheet is rigidly
        # displaced (shape, and hence coplanarity, preserved)
        offsets = [wake.nodes[1][:, i, j] .- initial[1][:, i, j]
                   for i in 1:wake.nwakes[]+1, j in axes(wake.nodes[1], 3)]
        @test all(o -> isapprox(o, offsets[1]; atol=1e-14), offsets)

        # legacy modes are unchanged by the new branch
        for swiv in (true, false)
            legacy = pnl.PanelWake(body; nwakerows=2,
                shed_with_induced_velocity=swiv)
            legacy.nwakes[] = 1
            pnl.apply_freestream!(legacy, uinf)
            legacy.velocity[1][:, 1, :] .+= induced_first
            legacy.velocity[1][:, 2, :] .+= induced_downstream
            pnl.propagate!(legacy, dt)
            expected_row1 = swiv ? dt .* (uinf .+ induced_first) : dt .* uinf
            @test legacy.nodes[1][:, 1, 1] ≈ expected_row1
            @test legacy.nodes[1][:, 2, 1] ≈ dt .* (uinf .+ induced_downstream)
        end
    end

    @testset "PanelWake optional final filament" begin
        body = make_plate_vortex_body()
        default_wake = pnl.PanelWake(body; nwakerows=2)
        panel_only = pnl.PanelWake(body; nwakerows=2,
            include_final_filament=false)
        @test length(pnl.get_sources(default_wake)) == 2
        @test pnl.get_sources(panel_only) == (panel_only,)
        @test default_wake.include_final_filament
        @test !panel_only.include_final_filament
    end

    @testset "PanelParticleWake forwards shedding velocity option" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; shed_with_induced_velocity=false)

        @test !wake.panel_wake.shed_with_induced_velocity
    end

    @testset "PanelParticleWake stores particle kernel offset" begin
        body = make_plate_vortex_body()
        default_wake = pnl.PanelParticleWake(body)
        offset_wake = pnl.PanelParticleWake(body; particle_core_size=0.25)

        @test isnan(default_wake.particle_core_size)
        @test offset_wake.particle_core_size == 0.25
        @test body.core_size_targets == 0.25
    end

    @testset "PanelParticleWake forwards SFS model" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; SFS=FLOWVPM.SFS_Cd_twolevel_nobackscatter)

        @test wake.pfield.SFS === FLOWVPM.SFS_Cd_twolevel_nobackscatter
        @test FLOWVPM.isLES(wake.pfield)
    end

    @testset "PanelParticleWake forwards viscous model" begin
        body = make_plate_vortex_body()
        viscous = FLOWVPM.CoreSpreading(1.5e-5, 0.01, FLOWVPM.zeta_fmm)
        wake = pnl.PanelParticleWake(body; viscous)

        @test wake.pfield.viscous === viscous
        @test FLOWVPM.iscorespreading(wake.pfield.viscous)
    end

    @testset "PanelWake final filament strength option" begin
        body = make_plate_vortex_body()

        function buffered_filament_strength(; unsteady_filament)
            wake = pnl.PanelWake(body; nwakerows=2, unsteady_filament)
            wake.nwakes[] = 2
            wake.overflowed[] = true
            wake.strength[1][1, 2, 1] = 11.0
            wake.strength[1][1, 3, 1] = 22.0
            wake.nodes[1][:, 3, 1] .= [0.0, 0.0, 0.0]
            wake.nodes[1][:, 3, 2] .= [1.0, 0.0, 0.0]

            buffer = zeros(12, 1)
            FastMultipole.source_system_to_buffer!(buffer, 1, pnl.FilamentWrapper(wake), 1)
            return buffer[5, 1]
        end

        @test buffered_filament_strength(unsteady_filament=true) == -22.0
        @test buffered_filament_strength(unsteady_filament=false) == -11.0
    end

    @testset "PanelWake steady final filament cancels last wake row influence" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelWake(body; nwakerows=1, core_size=1e-3, unsteady_filament=false)
        wake.nwakes[] = 1
        wake.overflowed[] = true

        Γ = 2.5
        wake.strength[1][1, 1, 1] = Γ
        wake.strength[1][1, 2, 1] = -3.0 # ignored when unsteady_filament=false

        v1 = SVector(0.0, 0.0, 0.0)
        v2 = SVector(1.0, 0.0, 0.0)
        v3 = SVector(1.0, 1.0, 0.0)
        v4 = SVector(0.0, 1.0, 0.0)
        wake.nodes[1][:, 1, 1] .= v1
        wake.nodes[1][:, 2, 1] .= v2
        wake.nodes[1][:, 2, 2] .= v3
        wake.nodes[1][:, 1, 2] .= v4

        switch = FastMultipole.DerivativesSwitch(false, true, true)
        panel_buffer = zeros(FastMultipole.data_per_body(wake), 1)
        FastMultipole.source_system_to_buffer!(panel_buffer, 1, wake, 1)

        filament = pnl.FilamentWrapper(wake)
        filament_buffer = zeros(FastMultipole.data_per_body(filament), 1)
        FastMultipole.source_system_to_buffer!(filament_buffer, 1, filament, 1)

        function filament_velocity(target, a, b, strength)
            return pnl._bound_vortex_velocity(target - a, target - b, true, wake.core_size) * strength
        end

        function filament_gradient(target, a, b, strength)
            return pnl._bound_vortex_gradient(a - target, b - target, true, wake.core_size) * strength
        end

        for target in (
            SVector(0.25, 0.35, 0.7),
            SVector(-0.4, 0.8, 0.5),
            SVector(1.0 + 1e-6, 0.45, 1e-6),
        )
            _, panel_velocity, panel_gradient = pnl.induced(target, wake, panel_buffer, 1, switch)
            final_velocity = filament_velocity(target, v2, v3, filament_buffer[5, 1])
            final_gradient = filament_gradient(target, v2, v3, filament_buffer[5, 1])
            combined_velocity = panel_velocity + final_velocity
            combined_gradient = panel_gradient + final_gradient

            expected_velocity =
                filament_velocity(target, v1, v2, Γ) +
                filament_velocity(target, v3, v4, Γ) +
                filament_velocity(target, v4, v1, Γ)
            expected_gradient =
                filament_gradient(target, v1, v2, Γ) +
                filament_gradient(target, v3, v4, Γ) +
                filament_gradient(target, v4, v1, Γ)

            @test combined_velocity ≈ expected_velocity atol=1e-12 rtol=1e-12
            @test combined_gradient ≈ expected_gradient atol=1e-10 rtol=1e-10
        end
    end

    @testset "PanelParticleWake forwards unsteady filament option" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; unsteady_filament=false)

        @test !wake.panel_wake.unsteady_filament
    end

    @testset "Kernel offset override regularizes body-to-particle influence" begin
        body = make_plate_vortex_body()
        body.core_size = 1e-8
        wake = pnl.PanelParticleWake(body; particle_core_size=0.1)
        FLOWVPM.add_particle(wake.pfield, [0.5, 1e-4, 1e-4], [0.0, 0.0, 0.0], 0.01)

        function particle_speed(core_size)
            FLOWVPM._reset_particles(wake.pfield)
            old = body.core_size
            body.core_size = core_size
            try
                pnl.influence!((wake.pfield,), (body,), pnl.DirectBackend();
                    velocity=true, velocity_gradient=(false,))
                return norm(view(wake.pfield.particles, FLOWVPM.U_INDEX, 1))
            finally
                body.core_size = old
            end
        end

        baseline = particle_speed(body.core_size_panel)
        regularized = particle_speed(body.core_size_targets)

        @test regularized < baseline
    end

    @testset "Self panel conditioning leaves particle targets on target offset" begin
        body = make_plate_vortex_body()
        body.core_size_panel = 1e-8
        body.core_size_targets = 0.1
        body.core_size = body.core_size_targets
        wake = pnl.PanelParticleWake(body; particle_core_size=body.core_size_targets)
        FLOWVPM.add_particle(wake.pfield, [0.5, 1e-4, 1e-4], [0.0, 0.0, 0.0], 0.01)

        FLOWVPM._reset_particles(wake.pfield)
        body.velocity .= 0
        pnl.influence!((body, wake.pfield), (body,), pnl.DirectBackend();
            velocity=true,
            velocity_gradient=(false, false),
            direct_conditioning=pnl._self_panel_core_size_conditioning())
        conditioned_velocity = copy(view(wake.pfield.particles, FLOWVPM.U_INDEX, 1))

        FLOWVPM._reset_particles(wake.pfield)
        body.core_size = body.core_size_targets
        pnl.influence!((wake.pfield,), (body,), pnl.DirectBackend();
            velocity=true, velocity_gradient=(false,))
        target_offset_velocity = copy(view(wake.pfield.particles, FLOWVPM.U_INDEX, 1))

        @test body.core_size == body.core_size_targets
        @test conditioned_velocity ≈ target_offset_velocity atol=1e-12 rtol=1e-12
    end

    @testset "ParticleMaintenance separates mixed policies" begin
        maintenance = pnl.ParticleMaintenance((
            pnl.MinGamma(1e-3),
            pnl.MergeParticles(every=2, r_hash=0.25),
            pnl.GlobalBox([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]),
        ))

        @test maintenance.trim_policies isa Tuple{<:pnl.MinGamma,<:pnl.GlobalBox}
        @test maintenance.functional_policies isa Tuple{<:pnl.MergeParticles}
        @test_throws ArgumentError pnl.ParticleMaintenance((pnl.MinGamma(1e-3), :not_a_policy))
    end

    @testset "PanelParticleWake particle maintenance trims by strength" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.MinGamma(1e-3),)))

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [1e-4, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [3.0, 0.0, 0.0], [5e-4, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 2
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 1)) == [2e-3, 0.0, 0.0]
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 2)) == [1e-3, 0.0, 0.0]
    end

    @testset "PanelParticleWake particle maintenance trims by strength range" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.MinGamma(1e-3), pnl.MaxGamma(3e-3))))

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [5e-4, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [4e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 1
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 1)) == [2e-3, 0.0, 0.0]
    end

    @testset "PanelParticleWake particle maintenance trims by global box" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.GlobalBox([0.0, -1.0, -1.0], [2.0, 1.0, 1.0]),)))

        FLOWVPM.add_particle(wake.pfield, [-1.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [3e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [3.0, 0.0, 0.0], [4e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 2
        @test sort([FLOWVPM.get_X(wake.pfield, i)[1] for i in 1:wake.pfield.np]) == [0.0, 2.0]
    end

    @testset "PanelParticleWake particle maintenance trims by global cylinder" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.GlobalCylinder([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0),)))

        FLOWVPM.add_particle(wake.pfield, [-0.1, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.5, 0.0], [3e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 1.1, 0.0], [4e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 1.0, 0.0], [5e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.1, 0.0, 0.0], [6e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        retained = [SVector{3}(FLOWVPM.get_X(wake.pfield, i)) for i in 1:wake.pfield.np]
        @test wake.pfield.np == 3
        @test SVector(0.0, 0.0, 0.0) in retained
        @test SVector(1.0, 0.5, 0.0) in retained
        @test SVector(2.0, 1.0, 0.0) in retained
        @test_throws ArgumentError pnl.GlobalCylinder([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0)
        @test_throws ArgumentError pnl.GlobalCylinder([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], -1.0)
    end

    @testset "PanelParticleWake particle maintenance trims by frame box" begin
        body = make_plate_vortex_body()
        frames = pnl.ReferenceFrame(body; origin=SVector(10.0, 0.0, 0.0))
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.FrameBox(1, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]),)))

        FLOWVPM.add_particle(wake.pfield, [10.5, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [12.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)

        @test_throws ArgumentError pnl.propagate!(wake, 0.0; relax=false, step=1)

        pnl.propagate!(wake, 0.0; relax=false, step=1, frames)

        @test wake.pfield.np == 1
        @test collect(FLOWVPM.get_X(wake.pfield, 1)) == [10.5, 0.0, 0.0]
    end

    @testset "PanelParticleWake particle maintenance merges particles" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.MergeParticles(every=1, r=0.5, r_hash=0.5),)))

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [0.1, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 1
    end

    @testset "RelaxationPlaneFilter validation and frame refresh" begin
        body = make_plate_vortex_body()

        @test_throws ArgumentError pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
        @test_throws ArgumentError pnl.RelaxationPlaneFilter([0.0, 0.0], [1.0, 0.0, 0.0])
        @test_throws ArgumentError pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [Inf, 0.0, 0.0])

        global_filter = pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        global_wake = pnl.PanelParticleWake(body;
            relaxation=FLOWVPM.Relaxation(FLOWVPM.relax_correctedpedrizzetti, 1, 0.3, global_filter))
        pnl.propagate!(global_wake, 0.0; relax=false, step=1)

        frame_filter = pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]; i_frame=1)
        frame_wake = pnl.PanelParticleWake(body;
            relaxation=FLOWVPM.Relaxation(FLOWVPM.relax_correctedpedrizzetti, 1, 0.3, frame_filter))
        @test_throws ArgumentError pnl.propagate!(frame_wake, 0.0; relax=false, step=1)

        bad_frame_filter = pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]; i_frame=99)
        bad_frame_wake = pnl.PanelParticleWake(body;
            relaxation=FLOWVPM.Relaxation(FLOWVPM.relax_correctedpedrizzetti, 1, 0.3, bad_frame_filter))
        frames = pnl.ReferenceFrame(body)
        @test_throws ArgumentError pnl.propagate!(bad_frame_wake, 0.0; relax=false, step=1, frames)
    end

    @testset "Unified wake metadata round trip" begin
        body = make_plate_vortex_body()
        relaxation = pnl.plane_filtered_relaxation(
            FLOWVPM.relaxation_correctedpedrizzetti,
            SVector(0.0, 0.0, -2.0),
            SVector(0.0, 0.0, -2.0);
            i_frame=1)
        wake = pnl.PanelParticleWake(body;
            nwakerows=2,
            max_particles=128,
            method_trailing=pnl.NoShed(),
            method_unsteady=pnl.SigmaOverlap(0.25, 4.0),
            particle_maintenance=pnl.ParticleMaintenance((
                pnl.MinGamma(1e-3),
                pnl.GlobalCylinder([0.0, -1.0, -1.0], [2.0, 0.0, 0.0], 1.5),
                pnl.MergeParticles(every=2, r=0.4, r_hash=0.3, sigma_relative=false, max_sigma_ratio=1.7, skip_static=false),
            )),
            viscous=FLOWVPM.CoreSpreading(1.5e-5, 0.01, FLOWVPM.zeta_fmm; beta=1.5),
            SFS=FLOWVPM.SFS_Cd_twolevel_nobackscatter,
            relaxation,
        )
        frames = pnl.ReferenceFrame(body)
        path = mktempdir()
        pnl._write_metadata_toml(path, "run", (body,), (wake,), frames, [0.0, 0.1],
            (pnl.Backslash(body),), pnl.DirectBackend(), pnl.DirectBackend(), pnl.DirectBackend(), ())

        metadata = TOML.parsefile(joinpath(path, "run.metadata.toml"))
        @test metadata["wake"][1]["method_trailing"]["type"] == "NoShed"
        @test metadata["wake"][1]["pfield_optargs"]["viscous"]["type"] == "FLOWVPM.CoreSpreading"
        @test metadata["wake"][1]["pfield_optargs"]["SFS"]["type"] == "FLOWVPM.SFS_Cd_twolevel_nobackscatter"
        @test metadata["wake"][1]["particle_maintenance"]["trim_policies"][2]["type"] == "GlobalCylinder"
        @test metadata["wake"][1]["particle_maintenance"]["trim_policies"][2]["origin"] == [0.0, -1.0, -1.0]
        @test metadata["wake"][1]["particle_maintenance"]["trim_policies"][2]["extrude"] == [2.0, 0.0, 0.0]
        @test metadata["wake"][1]["particle_maintenance"]["trim_policies"][2]["radius"] == 1.5
        @test metadata["wake"][1]["pfield_optargs"]["relaxation"]["filter"]["type"] == "RelaxationPlaneFilter"
        @test metadata["wake"][1]["pfield_optargs"]["relaxation"]["filter"]["point"] == [0.0, 0.0, -2.0]
        @test metadata["wake"][1]["pfield_optargs"]["relaxation"]["filter"]["normal"] == [0.0, 0.0, -1.0]
        @test metadata["wake"][1]["pfield_optargs"]["relaxation"]["filter"]["i_frame"] == 1

        reconstructed = pnl._construct_wakes_from_manifest((body,), metadata)
        @test reconstructed[1] isa pnl.PanelParticleWake
        @test reconstructed[1].method_trailing isa pnl.NoShed
        @test reconstructed[1].method_unsteady isa pnl.SigmaOverlap
        @test reconstructed[1].particle_maintenance.trim_policies[2] isa pnl.GlobalCylinder
        @test reconstructed[1].particle_maintenance.functional_policies[1] isa pnl.MergeParticles
        @test reconstructed[1].pfield.SFS === FLOWVPM.SFS_Cd_twolevel_nobackscatter
        reconstructed_filter = reconstructed[1].pfield.relaxation.filter
        @test reconstructed_filter isa pnl.RelaxationPlaneFilter
        @test reconstructed_filter.point_local == SVector(0.0, 0.0, -2.0)
        @test reconstructed_filter.normal_local == SVector(0.0, 0.0, -1.0)
        @test reconstructed_filter.i_frame == 1

        # Station-indexed shedding (018 Phase 16) survives the round trip with
        # its per-surface sigma vectors intact.
        wake_cs = pnl.PanelParticleWake(body;
            nwakerows=2, max_particles=128,
            method_trailing=pnl.StationSigmaOverlap([[0.1, 0.2, 0.3]], 2.75),
            method_unsteady=pnl.NoShed())
        path_cs = mktempdir()
        pnl._write_metadata_toml(path_cs, "run", (body,), (wake_cs,), frames, [0.0, 0.1],
            (pnl.Backslash(body),), pnl.DirectBackend(), pnl.DirectBackend(), pnl.DirectBackend(), ())
        metadata_cs = TOML.parsefile(joinpath(path_cs, "run.metadata.toml"))
        @test metadata_cs["wake"][1]["method_trailing"]["type"] == "StationSigmaOverlap"
        reconstructed_cs = pnl._construct_wakes_from_manifest((body,), metadata_cs)
        @test reconstructed_cs[1].method_trailing isa pnl.StationSigmaOverlap
        @test reconstructed_cs[1].method_trailing.sigmas == [[0.1, 0.2, 0.3]]
        @test reconstructed_cs[1].method_trailing.overlap == 2.75

        step1_frames = pnl.ReferenceFrame(body; ω_axis=SVector(0.0, 1.0, 0.0), ω=1.0)
        step2_frames = pnl.ReferenceFrame(body; ω_axis=SVector(0.0, 1.0, 0.0), ω=2.0)
        replacement_frames = pnl.ReferenceFrame(body; ω_axis=SVector(0.0, 1.0, 0.0), ω=3.0)

        pnl._append_metadata_step_toml(path, "run", step1_frames, 1, 0.1)
        pnl._append_metadata_step_toml(path, "run", step2_frames, 2, 0.2)
        pnl._append_metadata_step_toml(path, "run", replacement_frames, 2, 0.25)

        metadata = TOML.parsefile(joinpath(path, "run.metadata.toml"))
        @test length(metadata["step"]) == 2
        @test pnl._metadata_step_frames(metadata, 1)[1]["omega"] == 1.0
        @test pnl._metadata_step_frames(metadata, 2)[1]["omega"] == 3.0
        @test pnl._metadata_step_frames(metadata, 3) === nothing
    end

    @testset "Convert-at-shed nwakerows=0 (BRAINSTORM 024)" begin
        # N=0 mode: no free wake-panel row survives a solve; shed_wake! shifts,
        # writes the fresh row, pins it to TE+Das, converts it to particles in
        # the same call, and resets nwakes[] to zero. Storage is the N=1 layout
        # with convert_at_shed=true.

        body = make_dirichlet_diamond_body(; nspan=3)
        body.strength[:, 2] .= range(0.4, 1.6; length=size(body.strength, 1))

        #--- construction and validation ---#
        w = pnl.PanelParticleWake(body; nwakerows=0, max_particles=100)
        @test w.panel_wake.convert_at_shed
        @test size(w.panel_wake.nodes[1], 2) == 2      # N=1 storage underneath
        @test pnl._logical_nwakerows(w.panel_wake) == 0
        @test !pnl.PanelParticleWake(body; nwakerows=2,
            max_particles=100).panel_wake.convert_at_shed
        @test pnl._logical_nwakerows(pnl.PanelParticleWake(body; nwakerows=2,
            max_particles=100).panel_wake) == 2

        @test_throws ArgumentError pnl.PanelParticleWake(body; nwakerows=-1,
            max_particles=100)
        @test_throws ArgumentError pnl.PanelWake(body; nwakerows=0)
        # The retained final filament at the Das line is the sheet/particle
        # interface carrier: N=0 cannot opt out of it.
        @test_throws ArgumentError pnl.PanelParticleWake(body; nwakerows=0,
            max_particles=100, unsteady_filament=false)
        @test_throws ArgumentError pnl.PanelParticleWake(body; nwakerows=0,
            max_particles=100, include_final_filament=false)
        # Smooth conversion: only :downstream conserves circulation at N=0
        # (upstream/split would strand the per-step unsteady jump).
        @test_throws ArgumentError pnl.PanelParticleWake(body; nwakerows=0,
            max_particles=100, conversion=pnl.SurfaceVorticityConversion(0.1))
        @test_throws ArgumentError pnl.PanelParticleWake(body; nwakerows=0,
            max_particles=100,
            conversion=pnl.SurfaceVorticityConversion(0.1; attribution=:split))
        @test pnl.PanelParticleWake(body; nwakerows=0, max_particles=100,
            conversion=pnl.SurfaceVorticityConversion(0.1; attribution=:downstream)
            ) isa pnl.PanelParticleWake

        #--- shed-time conversion equals a matched N=1 conversion ---#
        # Drive the N=0 wake through one real shed_wake! call, then build an
        # N=1 wake holding the identical conversion-time state (same fresh-row
        # geometry TE+Das -> convected line, same fresh and history strengths)
        # and call _convert_to_particles! directly. Particle output and the
        # post-shed filament bookkeeping must agree exactly.
        function n0_vs_n1(; conversion_kwargs...)
            ncols = size(body.shedding[1], 2) + 1
            S_prev = collect(range(0.2, 0.5; length=ncols - 1))

            w0 = pnl.PanelParticleWake(body; nwakerows=0, max_particles=2000,
                conversion_kwargs...)
            pw0 = w0.panel_wake
            for c in 1:ncols     # pre-shed row 1 = last step's convected line
                pw0.nodes[1][:, 1, c] .= (0.7, c - 1.0, 0.05)
            end
            pw0.strength[1][1, 1, :] .= S_prev
            pnl.shed_wake!(w0, body)

            @test pw0.nwakes[] == 0        # no free sheet survives
            @test pw0.overflowed[]         # filament/VTK guards active

            w1 = pnl.PanelParticleWake(body; nwakerows=1, max_particles=2000,
                conversion_kwargs...)
            pw1 = w1.panel_wake
            pnl.update_TE!(pw1, body)                  # row 1 = TE + Das line
            for c in 1:ncols
                pw1.nodes[1][:, 2, c] .= (0.7, c - 1.0, 0.05)
            end
            for j in 1:(ncols - 1)         # fresh strengths, shed_wake!'s expression
                si, sj = pnl._get_wakestrength_mu(body, j, 1)
                pw1.strength[1][1, 1, j] = si - sj
            end
            pw1.strength[1][1, 2, :] .= S_prev
            pw1.nwakes[] = 1
            pnl._convert_to_particles!(w1, body)

            np = w0.pfield.np
            @test np > 0
            @test w1.pfield.np == np
            @test view(w0.pfield.particles, :, 1:np) ==
                  view(w1.pfield.particles, :, 1:np)

            # post-shed wake state: fresh row + carried history + pinned line
            @test pw0.strength[1][1, 1, :] == pw1.strength[1][1, 1, :]
            @test pw0.strength[1][1, 2, :] == S_prev
            @test pw0.nodes[1][:, 1, :] == pw1.nodes[1][:, 1, :]
            @test pw0.nodes[1][:, 2, :] == pw1.nodes[1][:, 2, :]

            # final filament: sits on the Das line (node row 1), cancelling the
            # full just-shed strength — the sheet/particle interface carrier.
            @test FastMultipole.get_n_bodies(pnl.FilamentWrapper(pw0)) == ncols - 1
            @test FastMultipole.get_n_bodies(pw0) == 0     # no panel sources
            for j in 1:(ncols - 1)
                @test pnl._final_filament_strength(pw0, 1, pw0.nwakes[], j) ==
                      -pw0.strength[1][1, 1, j]
            end
            return w0
        end

        w0_legacy = n0_vs_n1()   # LegacyEdgeJumpConversion defaults
        w0_smooth = n0_vs_n1(;
            conversion=pnl.SurfaceVorticityConversion(0.08; attribution=:downstream))
        @test w0_smooth.conversion_count[] == 1
        @test w0_smooth.panel_wake.particle_handoff_active[]
        @test w0_smooth.panel_wake.particle_handoff_weight[] == 0.0
        # ledger closes exactly (the :downstream attribution deposits the whole
        # unsteady downstream face; the zero upstream face stays at the rigid row)
        diag = w0_smooth.conversion_diagnostics
        @test norm(diag.total_residual) < 1e-12
        @test diag.total_deposited ≈ diag.total_expected atol=1e-12

        #--- second shed: strength history and repeated conversion ---#
        pw = w0_legacy.panel_wake
        S1 = copy(pw.strength[1][1, 1, :])
        np1 = w0_legacy.pfield.np
        body.strength[:, 2] .*= 1.1                    # new "solved" state
        pw.nodes[1][:, 1, :] .+= [0.05, 0.0, 0.0]      # emulate convection
        pnl.shed_wake!(w0_legacy, body)
        @test pw.nwakes[] == 0
        @test pw.strength[1][1, 2, :] == S1            # history carried
        @test w0_legacy.pfield.np > np1                # converted again
        @test all(isfinite,
            w0_legacy.pfield.particles[:, 1:w0_legacy.pfield.np])
        body.strength[:, 2] ./= 1.1                    # restore

        #--- VTK write / warmstart-load round trip at nwakes[] == 0 ---#
        dir = mktempdir()
        pnl.write_vtk(joinpath(dir, "w"), w0_legacy, 3, 0.3; overwrite=true)
        w2 = pnl.PanelParticleWake(body; nwakerows=0, max_particles=2000)
        pnl._load_panel_particle_wake_vtk!(w2, dir, "w", 3)
        @test w2.panel_wake.nwakes[] == 0
        @test w2.panel_wake.overflowed[]               # idx >= 1 => post-first-shed
        @test w2.panel_wake.nodes[1][:, 1, :] == pw.nodes[1][:, 1, :]
        @test w2.pfield.np == w0_legacy.pfield.np
        @test view(w2.pfield.particles, 1:3, 1:w2.pfield.np) ==
              view(w0_legacy.pfield.particles, 1:3, 1:w0_legacy.pfield.np)
        @test view(w2.pfield.particles, FLOWVPM.GAMMA_INDEX, 1:w2.pfield.np) ==
              view(w0_legacy.pfield.particles, FLOWVPM.GAMMA_INDEX,
                   1:w0_legacy.pfield.np)
    end
end
