using Test
import FLOWPanel as pnl
import FLOWVPM
const TOML = pnl.TOML
import FastMultipole
using LinearAlgebra: norm
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
        offset_wake = pnl.PanelParticleWake(body; particle_kerneloffset=0.25)

        @test isnan(default_wake.particle_kerneloffset)
        @test offset_wake.particle_kerneloffset == 0.25
        @test body.kerneloffset_targets == 0.25
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
        body.kerneloffset = 1e-8
        wake = pnl.PanelParticleWake(body; particle_kerneloffset=0.1)
        FLOWVPM.add_particle(wake.pfield, [0.5, 1e-4, 1e-4], [0.0, 0.0, 0.0], 0.01)

        function particle_speed(kerneloffset)
            FLOWVPM._reset_particles(wake.pfield)
            old = body.kerneloffset
            body.kerneloffset = kerneloffset
            try
                pnl.influence!((wake.pfield,), (body,), pnl.DirectBackend();
                    velocity=true, velocity_gradient=(false,))
                return norm(view(wake.pfield.particles, FLOWVPM.U_INDEX, 1))
            finally
                body.kerneloffset = old
            end
        end

        baseline = particle_speed(body.kerneloffset_panel)
        regularized = particle_speed(body.kerneloffset_targets)

        @test regularized < baseline
    end

    @testset "Self panel conditioning leaves particle targets on target offset" begin
        body = make_plate_vortex_body()
        body.kerneloffset_panel = 1e-8
        body.kerneloffset_targets = 0.1
        body.kerneloffset = body.kerneloffset_targets
        wake = pnl.PanelParticleWake(body; particle_kerneloffset=body.kerneloffset_targets)
        FLOWVPM.add_particle(wake.pfield, [0.5, 1e-4, 1e-4], [0.0, 0.0, 0.0], 0.01)

        FLOWVPM._reset_particles(wake.pfield)
        body.velocity .= 0
        pnl.influence!((body, wake.pfield), (body,), pnl.DirectBackend();
            velocity=true,
            velocity_gradient=(false, false),
            direct_conditioning=pnl._self_panel_kerneloffset_conditioning())
        conditioned_velocity = copy(view(wake.pfield.particles, FLOWVPM.U_INDEX, 1))

        FLOWVPM._reset_particles(wake.pfield)
        body.kerneloffset = body.kerneloffset_targets
        pnl.influence!((wake.pfield,), (body,), pnl.DirectBackend();
            velocity=true, velocity_gradient=(false,))
        target_offset_velocity = copy(view(wake.pfield.particles, FLOWVPM.U_INDEX, 1))

        @test body.kerneloffset == body.kerneloffset_targets
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
end
