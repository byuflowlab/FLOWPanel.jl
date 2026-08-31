#=##############################################################################
Mechanical unit tests for the BRAINSTORM-015 Phase 3 Route B
(TEAnchoredAttachment) runtime: cold startup, live-block ownership, source-
view exclusion / no double counting, topology advancement, metadata and
warm-start restoration, and the B2 live-block seam (architecture phase_02
§12 categories 7, 8, 9, 20, 21).
=###############################################################################

using Test
import FLOWPanel as pnl
import FastMultipole
import LinearAlgebra as _KLB

if !isdefined(@__MODULE__, :make_dirichlet_diamond_body)
    include("test_helpers.jl")
end

const KUTTAB_UINF = t -> [1.0, 0.0, 0.3]
const KUTTAB_MANEUVER = (frames, systems, wakes, t) -> nothing
const KUTTAB_TRANGE = collect(range(0.0; step=0.05, length=5))
const KUTTAB_TOML = pnl.TOML

function kuttab_case(; nwakerows=8)
    body = make_dirichlet_diamond_body()
    wake = pnl.PanelWake(body; nwakerows)
    frames = pnl.ReferenceFrame(body; name="vehicle")
    solver = pnl.Backslash(body)
    return body, wake, frames, solver
end

function run_kuttab!(body, wake, frames, solver; t_range=KUTTAB_TRANGE,
        path=nothing, kwargs...)
    return pnl.simulate!((body,), (wake,), frames, KUTTAB_MANEUVER,
        KUTTAB_UINF, t_range; body_solvers=(solver,),
        backend=pnl.DirectBackend(), path, name="run",
        grad_mu_options=(; basis=:tri), kwargs...)
end

"Sum the influence of every panel a PanelWake exposes as an old-wake source,
via the same buffer path the backends use."
function wake_source_potential(wake::pnl.PanelWake, target)
    n = FastMultipole.get_n_bodies(wake)
    nd = FastMultipole.data_per_body(wake)
    buffer = zeros(nd, max(n, 1))
    phi = 0.0
    for i in 1:n
        FastMultipole.source_system_to_buffer!(buffer, i, wake, i)
        p, _, _ = pnl.induced(target, wake, buffer, i,
            FastMultipole.DerivativesSwitch(true, true, false))
        phi += p
    end
    return phi
end

@testset verbose=true "Route B / live block (BRAINSTORM 015)" begin

    # --- category 8 (mechanical): source-view exclusion arithmetic ---
    @testset "old-wake source view excludes the live block" begin
        body = make_dirichlet_diamond_body(; nspan=2)
        wake = pnl.PanelWake(body; nwakerows=6)
        ncols = size(wake.strength[1], 3)
        wake.nwakes[] = 3

        # legacy: no live block, all rows are sources
        @test wake.live_rows[] == 0
        @test FastMultipole.get_n_bodies(wake) == 3*ncols
        @test pnl.global_to_matrix_index(wake, 1) == (1, 1, 1)

        # one reserved live row: row 1 disappears from the source view
        wake.live_rows[] = 1
        @test FastMultipole.get_n_bodies(wake) == 2*ncols
        @test pnl.global_to_matrix_index(wake, 1) == (1, 2, 1)
        # round trip through the inverse map for every exposed source
        for i in 1:FastMultipole.get_n_bodies(wake)
            isurf, irow, icol = pnl.global_to_matrix_index(wake, i)
            @test irow >= 2      # never the live row
            @test pnl.matrix_to_global_index(wake, isurf, irow, icol) == i
        end

        # B2 seam: a multi-row live block uses the same arithmetic (category 21)
        wake.live_rows[] = 2
        @test FastMultipole.get_n_bodies(wake) == 1*ncols
        for i in 1:FastMultipole.get_n_bodies(wake)
            isurf, irow, icol = pnl.global_to_matrix_index(wake, i)
            @test irow >= 3
            @test pnl.matrix_to_global_index(wake, isurf, irow, icol) == i
        end
        wake.live_rows[] = 0
    end

    # --- category 8 (physical): no double counting ---
    @testset "live-block exclusion removes exactly the live row" begin
        body, wake, frames, solver = kuttab_case()
        run_kuttab!(body, wake, frames, solver;
            wake_attachment=pnl.TEAnchoredAttachment(),
            kutta_closure=pnl.JumpKutta())
        @test wake.live_rows[] == 1
        target = FastMultipole.SVector{3}(0.4, 0.5, 0.8)

        # brute force over ALL rows vs the excluded view + live row alone
        phi_view = wake_source_potential(wake, target)
        wake.live_rows[] = 0
        phi_all = wake_source_potential(wake, target)
        wake.live_rows[] = 1
        # the difference must be exactly the live (row 1) panels' influence,
        # i.e. the block the body-side attachment operator represents instead
        nd = FastMultipole.data_per_body(wake)
        wake0 = pnl.PanelWake(body; nwakerows=2)
        wake0.nwakes[] = 1
        for i_surf in eachindex(wake0.nodes)
            wake0.nodes[i_surf][:, 1:2, :] .= wake.nodes[i_surf][:, 1:2, :]
            wake0.strength[i_surf][:, 1, :] .= wake.strength[i_surf][:, 1, :]
        end
        phi_live = wake_source_potential(wake0, target)
        @test isapprox(phi_all - phi_view, phi_live; rtol=1e-12, atol=1e-14)

        # the body-side ephemeral cache places the attached panel exactly on
        # the live panel: Das = first convected row − TE
        shed = body.shedding[1]
        for e in axes(shed, 2)
            i_panel = shed[1, e]
            n = body.cells[shed[3, e], i_panel]
            te = body.nodes[:, n]
            @test isapprox(body.Das[1][:, e],
                wake.nodes[1][:, 2, e] .- te; atol=1e-13)
        end
    end

    # --- categories 7, 9: cold startup, activation, one write, advancement ---
    @testset "cold startup and per-step live-strength writes" begin
        body, wake, frames, solver = kuttab_case()
        clo = pnl.PressureContinuityKutta(pnl.SteadyBernoulliProvider(1.2))
        path = mktempdir()
        run_kuttab!(body, wake, frames, solver; path,
            wake_attachment=pnl.TEAnchoredAttachment(), kutta_closure=clo)
        diags = pnl.kutta_diagnostics(clo)
        @test length(diags) == length(KUTTAB_TRANGE)
        # t0 is the deterministic jump startup; the first pressure solve is t1
        @test diags[1].status === :startup_jump
        @test diags[1].startup === :startup_jump
        @test all(d -> d.status === :converged, diags[2:end])
        @test wake.live_rows[] == 1
        @test wake.live_step_id[] == length(KUTTAB_TRANGE) - 1

        # category 9: the committed circulation of each step, recorded in the
        # metadata step records, appears exactly once in the wake history:
        # after the final solve (no trailing shed), row r holds step
        # (last − r + 1)'s accepted circulation
        metadata = KUTTAB_TOML.parsefile(joinpath(path, "run.metadata.toml"))
        steps = metadata["step"]
        nsteps = length(steps)
        @test nsteps == length(KUTTAB_TRANGE)
        for (r, s) in enumerate(reverse(steps))
            r > wake.nwakes[] && break
            gamma = Float64.(s["kutta"]["gamma_accepted"])
            @test isapprox(vec(wake.strength[1][1, r, :]), gamma; atol=1e-12)
        end
        # manifest carries the configuration
        @test metadata["kutta"]["wake_attachment"] == "TEAnchoredAttachment"
        @test metadata["kutta"]["kutta_closure"] == "PressureContinuityKutta"
        @test metadata["kutta"]["provider"] == "SteadyBernoulliProvider"
    end

    # --- category 20: metadata round trip and warm-start restoration ---
    @testset "warm start restores accepted state and skips startup" begin
        for (attachment, closurefun, label) in (
                (pnl.TEAnchoredAttachment(), () -> pnl.JumpKutta(), "B/jump"),
                (pnl.TEAnchoredAttachment(),
                    () -> pnl.PressureContinuityKutta(pnl.SteadyBernoulliProvider(1.2)),
                    "B/pressure"),
                (pnl.RigidTransitionAttachment(),
                    () -> pnl.PressureContinuityKutta(pnl.SteadyBernoulliProvider(1.2)),
                    "A/pressure"))
            # full reference run
            bA, wA, fA, sA = kuttab_case()
            cA = closurefun()
            run_kuttab!(bA, wA, fA, sA;
                wake_attachment=attachment, kutta_closure=cA)
            # half run + warm-started continuation
            bB, wB, fB, sB = kuttab_case()
            pathB = mktempdir()
            run_kuttab!(bB, wB, fB, sB; t_range=KUTTAB_TRANGE[1:3], path=pathB,
                wake_attachment=attachment, kutta_closure=closurefun())
            bC, wC, fC, sC = kuttab_case()
            cC = closurefun()
            pnl.simulate_warmstart!((bC,), (wC,), fC, KUTTAB_MANEUVER,
                KUTTAB_UINF, KUTTAB_TRANGE;
                body_solvers=(sC,), backend=pnl.DirectBackend(), path=pathB,
                name="run", grad_mu_options=(; basis=:tri),
                wake_attachment=attachment, kutta_closure=cC)
            @test isapprox(bA.strength, bC.strength; rtol=1e-10, atol=1e-12)
            @test wA.nwakes[] == wC.nwakes[]
            @test isapprox(wA.strength[1], wC.strength[1]; rtol=1e-10, atol=1e-12)
            @test isapprox(wA.nodes[1], wC.nodes[1]; rtol=1e-10, atol=1e-12)
            @test wA.live_rows[] == wC.live_rows[]
            @test wA.live_step_id[] == wC.live_step_id[]
            # the resumed run never repeats :startup_jump (category 7)
            if cC isa pnl.PressureContinuityKutta
                @test all(d -> d.status !== :startup_jump,
                    pnl.kutta_diagnostics(cC))
            end
        end

        # configuration mismatch is rejected rather than reseeded
        bB, wB, fB, sB = kuttab_case()
        pathB = mktempdir()
        run_kuttab!(bB, wB, fB, sB; t_range=KUTTAB_TRANGE[1:3], path=pathB,
            wake_attachment=pnl.TEAnchoredAttachment(),
            kutta_closure=pnl.JumpKutta())
        bC, wC, fC, sC = kuttab_case()
        @test_throws ArgumentError pnl.simulate_warmstart!((bC,), (wC,), fC,
            KUTTAB_MANEUVER, KUTTAB_UINF, KUTTAB_TRANGE;
            body_solvers=(sC,), backend=pnl.DirectBackend(), path=pathB,
            name="run", grad_mu_options=(; basis=:tri))
        # ... and a legacy run cannot be resumed as Route B
        bD, wD, fD, sD = kuttab_case()
        pathD = mktempdir()
        run_kuttab!(bD, wD, fD, sD; t_range=KUTTAB_TRANGE[1:3], path=pathD)
        bE, wE, fE, sE = kuttab_case()
        @test_throws ArgumentError pnl.simulate_warmstart!((bE,), (wE,), fE,
            KUTTAB_MANEUVER, KUTTAB_UINF, KUTTAB_TRANGE;
            body_solvers=(sE,), backend=pnl.DirectBackend(), path=pathD,
            name="run", grad_mu_options=(; basis=:tri),
            wake_attachment=pnl.TEAnchoredAttachment(),
            kutta_closure=pnl.JumpKutta())

        # an unsupported configuration is rejected up front with an
        # ArgumentError (not a raw field error after mutation): resuming with
        # a PanelParticleWake under a non-default configuration
        bF = make_dirichlet_diamond_body()
        ppw = pnl.PanelParticleWake(bF; nwakerows=2, max_particles=100)
        state_before = copy(bF.strength)
        @test_throws "PanelParticleWake" pnl.simulate_warmstart!((bF,),
            (ppw,), pnl.ReferenceFrame(bF; name="vehicle"),
            KUTTAB_MANEUVER, KUTTAB_UINF, KUTTAB_TRANGE;
            body_solvers=(pnl.Backslash(bF),), backend=pnl.DirectBackend(),
            path=pathD, name="run", grad_mu_options=(; basis=:tri),
            wake_attachment=pnl.TEAnchoredAttachment(),
            kutta_closure=pnl.JumpKutta())
        @test bF.strength == state_before          # rejected before mutation
    end

    # --- replay refuses non-default-configuration manifests ---
    @testset "replay rejects a kutta manifest" begin
        body, wake, frames, solver = kuttab_case()
        path = mktempdir()
        run_kuttab!(body, wake, frames, solver; path,
            wake_attachment=pnl.TEAnchoredAttachment(),
            kutta_closure=pnl.JumpKutta())
        # replay has no attachment/closure plumbing: silently replaying this
        # run as legacy would reconstruct wrong fields
        @test_throws "replay does not yet support" pnl.replay(path, "run")
    end

    # --- category 21: equal-strength internal-edge cancellation ---
    @testset "B2 seam: equal-strength subdivision cancels internally" begin
        body = make_dirichlet_diamond_body()
        # one straight two-row wake sheet with equal strengths...
        wake2 = pnl.PanelWake(body; nwakerows=3)
        wake2.nwakes[] = 2
        # ...and the geometrically identical single-row sheet
        wake1 = pnl.PanelWake(body; nwakerows=3)
        wake1.nwakes[] = 1
        gamma = 0.37
        shed = body.shedding[1]
        for e in axes(shed, 2)
            i_panel = shed[1, e]
            te = body.nodes[:, body.cells[shed[3, e], i_panel]]
            for (w, rows) in ((wake2, 3), (wake1, 2))
                for r in 1:rows
                    w.nodes[1][:, r, e] .= te .+ (r - 1)/(rows - 1) .* [0.8, 0.0, 0.1]
                end
            end
        end
        # last node column
        te_end = body.nodes[:, body.cells[shed[2, end], shed[1, end]]]
        for (w, rows) in ((wake2, 3), (wake1, 2))
            for r in 1:rows
                w.nodes[1][:, r, end] .= te_end .+ (r - 1)/(rows - 1) .* [0.8, 0.0, 0.1]
            end
        end
        wake2.strength[1][1, 1:2, :] .= gamma
        wake1.strength[1][1, 1, :] .= gamma

        for target in (FastMultipole.SVector{3}(1.7, 0.4, 0.6),
                       FastMultipole.SVector{3}(0.9, 1.2, -0.4))
            phi2 = wake_source_potential(wake2, target)
            phi1 = wake_source_potential(wake1, target)
            @test isapprox(phi2, phi1; rtol=1e-10, atol=1e-13)
        end
    end
end
