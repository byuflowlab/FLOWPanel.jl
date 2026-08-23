#=##############################################################################
Mechanical unit tests for the BRAINSTORM-015 Phase 3 wake-attachment /
Kutta-closure runtime (architecture phase_02 §12, categories 1-6, 10-17, 19,
22 as applicable to the Backslash-only first implementation). These are
architecture tests; they claim no physical validation.
=###############################################################################

using Test
import FLOWPanel as pnl
import LinearAlgebra as _KLA

if !isdefined(@__MODULE__, :make_dirichlet_diamond_body)
    include("test_helpers.jl")
end

# ---- typed mock providers (concretely typed per §2.2) ----

"Provider returning NaN pressures: every trial is rejected as non-finite."
struct KuttaNaNProvider <: pnl.AbstractKuttaPressureProvider end
pnl.pressure_requirements(::KuttaNaNProvider) = (:relative_surface_velocity,)
pnl.evaluate_pressure!(P::AbstractMatrix, ::KuttaNaNProvider,
    ::pnl.KuttaTrialView) = (P .= NaN; P)
pnl.commit_pressure!(::KuttaNaNProvider, ::pnl.KuttaTrialView) = nothing
pnl.rollback_pressure!(::KuttaNaNProvider) = nothing

"Steady-Bernoulli-equivalent provider counting interface calls."
mutable struct KuttaCountingProvider <: pnl.AbstractKuttaPressureProvider
    rho::Float64
    evals::Int
    commits::Int
    rollbacks::Int
end
KuttaCountingProvider(rho) = KuttaCountingProvider(Float64(rho), 0, 0, 0)
pnl.pressure_requirements(::KuttaCountingProvider) = (:relative_surface_velocity,)
function pnl.evaluate_pressure!(P::AbstractMatrix, p::KuttaCountingProvider,
        trial::pnl.KuttaTrialView)
    p.evals += 1
    pnl.evaluate_pressure!(P, pnl.SteadyBernoulliProvider(p.rho), trial)
    return P
end
pnl.commit_pressure!(p::KuttaCountingProvider, ::pnl.KuttaTrialView) =
    (p.commits += 1; nothing)
pnl.rollback_pressure!(p::KuttaCountingProvider) = (p.rollbacks += 1; nothing)

"Steady-Bernoulli provider that turns non-finite after `fail_after` commits:
the run settles onto nonzero committed state before the failing step, so a
rollback test can distinguish restored from never-written."
mutable struct KuttaLateFailProvider <: pnl.AbstractKuttaPressureProvider
    rho::Float64
    fail_after::Int
    commits::Int
end
KuttaLateFailProvider(rho, fail_after) =
    KuttaLateFailProvider(Float64(rho), fail_after, 0)
pnl.pressure_requirements(::KuttaLateFailProvider) = (:relative_surface_velocity,)
function pnl.evaluate_pressure!(P::AbstractMatrix, p::KuttaLateFailProvider,
        trial::pnl.KuttaTrialView)
    if p.commits >= p.fail_after
        P .= NaN
    else
        pnl.evaluate_pressure!(P, pnl.SteadyBernoulliProvider(p.rho), trial)
    end
    return P
end
pnl.commit_pressure!(p::KuttaLateFailProvider, ::pnl.KuttaTrialView) =
    (p.commits += 1; nothing)
pnl.rollback_pressure!(::KuttaLateFailProvider) = nothing

"Steady-Bernoulli provider whose `commit_pressure!` throws: exercises the
§7.4 commit-cannot-complete guard."
struct KuttaCommitBombProvider <: pnl.AbstractKuttaPressureProvider
    rho::Float64
end
pnl.pressure_requirements(::KuttaCommitBombProvider) = (:relative_surface_velocity,)
pnl.evaluate_pressure!(P::AbstractMatrix, p::KuttaCommitBombProvider,
    trial::pnl.KuttaTrialView) =
    pnl.evaluate_pressure!(P, pnl.SteadyBernoulliProvider(p.rho), trial)
pnl.commit_pressure!(::KuttaCommitBombProvider, ::pnl.KuttaTrialView) =
    error("commit bomb")
pnl.rollback_pressure!(::KuttaCommitBombProvider) = nothing

# ---- shared drivers ----

const KUTTA_UINF = t -> [1.0, 0.0, 0.3]
const KUTTA_MANEUVER = (frames, systems, wakes, t) -> nothing
const KUTTA_TRANGE = collect(range(0.0; step=0.05, length=4))

function kutta_case(; nwakerows=8)
    body = make_dirichlet_diamond_body()
    wake = pnl.PanelWake(body; nwakerows)
    frames = pnl.ReferenceFrame(body; name="vehicle")
    solver = pnl.Backslash(body)
    return body, wake, frames, solver
end

function run_kutta!(body, wake, frames, solver; t_range=KUTTA_TRANGE, kwargs...)
    return pnl.simulate!((body,), (wake,), frames, KUTTA_MANEUVER, KUTTA_UINF,
        t_range; body_solvers=(solver,), backend=pnl.DirectBackend(),
        path=nothing, grad_mu_options=(; basis=:tri), kwargs...)
end

pressure_closure(provider=pnl.SteadyBernoulliProvider(1.2); kwargs...) =
    pnl.PressureContinuityKutta(provider; kwargs...)

@testset verbose=true "Kutta closure (BRAINSTORM 015)" begin

    # --- category 2: exact default dispatch into the legacy path ---
    @testset "default pair is the legacy path" begin
        b1, w1, f1, s1 = kutta_case()
        run_kutta!(b1, w1, f1, s1)
        b2, w2, f2, s2 = kutta_case()
        run_kutta!(b2, w2, f2, s2;
            wake_attachment=pnl.RigidTransitionAttachment(),
            kutta_closure=pnl.JumpKutta())
        # explicit legacy pair must produce bitwise-identical results to the
        # omitted-kwargs default
        @test b1.strength == b2.strength
        @test b1.velocity == b2.velocity
        @test w1.strength[1] == w2.strength[1]
        @test w1.nodes[1] == w2.nodes[1]
        # the legacy path allocates no live block
        @test w1.live_rows[] == 0
        @test w1.live_step_id[] == -1
        @test pnl._is_legacy_kutta(pnl.RigidTransitionAttachment(), pnl.JumpKutta())
        @test !pnl._is_legacy_kutta(pnl.TEAnchoredAttachment(), pnl.JumpKutta())
        @test !pnl._is_legacy_kutta(pnl.RigidTransitionAttachment(),
            pressure_closure())
    end

    # --- category 1: all four combinations run and stay finite ---
    @testset "four-combination matrix runs" begin
        for (att, clo, label) in (
                (pnl.RigidTransitionAttachment(), pressure_closure(), "A/pressure"),
                (pnl.TEAnchoredAttachment(), pnl.JumpKutta(), "B/jump"),
                (pnl.TEAnchoredAttachment(), pressure_closure(), "B/pressure"))
            body, wake, frames, solver = kutta_case()
            run_kutta!(body, wake, frames, solver;
                wake_attachment=att, kutta_closure=clo)
            @test all(isfinite, body.strength)
            @test all(isfinite, wake.strength[1])
            if att isa pnl.TEAnchoredAttachment
                @test wake.live_rows[] == 1
                @test wake.live_step_id[] == length(KUTTA_TRANGE) - 1
            end
            if clo isa pnl.PressureContinuityKutta
                diags = pnl.kutta_diagnostics(clo)
                @test length(diags) == length(KUTTA_TRANGE)
                @test all(d -> d.status in (:converged, :startup_jump), diags)
            end
        end
        @test isempty(pnl.kutta_diagnostics(pnl.JumpKutta()))
    end

    # --- the closure is not a silent no-op: leverage + independent residual ---
    @testset "pressure closure changes the solution" begin
        b1, w1, f1, s1 = kutta_case()
        run_kutta!(b1, w1, f1, s1)                       # A/jump baseline
        b2, w2, f2, s2 = kutta_case()
        clo = pressure_closure()
        run_kutta!(b2, w2, f2, s2;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo)
        diags = pnl.kutta_diagnostics(clo)
        # the accepted correction is genuinely nonzero on every step
        # (calibrated: max|c| ≈ 8.6e-3·C_s on this fixture)
        @test all(d -> maximum(abs, d.c) > 1e-3*d.correction_scale, diags)
        # and it moves the body solution well above solver noise
        # (calibrated: ≈2.2% of max|strength|)
        @test maximum(abs, b1.strength .- b2.strength) >
            1e-4*maximum(abs, b1.strength)
        # independent pressure-continuity verification: recompute the paired
        # centroid Δp directly from the post-run committed velocity with the
        # provider's own law — NOT from the driver's gate fields
        rho = 1.2
        P_s = diags[end].pressure_scale
        shed = b2.shedding[1]
        for i in axes(shed, 2)
            uu = _KLA.norm(view(b2.velocity, :, shed[1, i]))
            ul = _KLA.norm(view(b2.velocity, :, shed[4, i]))
            r_indep = rho/2*(ul^2 - uu^2)
            @test abs(r_indep)/P_s <= 2*clo.pressure_tolerance
        end
        # the jump baseline does NOT satisfy pressure continuity here (the
        # residual the closure removes is real): recompute on the jump run
        r_jump = maximum(axes(shed, 2)) do i
            uu = _KLA.norm(view(b1.velocity, :, shed[1, i]))
            ul = _KLA.norm(view(b1.velocity, :, shed[4, i]))
            abs(rho/2*(ul^2 - uu^2))
        end
        @test r_jump/P_s > 10*clo.pressure_tolerance
    end

    # --- the trial/commit affine add-on follows the self-pair panel offset ---
    @testset "reconstruction add-on uses core_size_panel" begin
        # distinct offsets: the proportional γ = Cμ strip part runs through the
        # self-pair conditioning at core_size_panel, so the affine −c part
        # must use the same offset or the two halves of one strip strength are
        # regularized differently.
        # Pinned to the Vatistas family (BRAINSTORM 025): the discriminator
        # needs a kernel whose value at control-point distances depends on the
        # offset; under the compact-support default both offsets are exactly
        # singular there and the fixture loses its teeth (the correctness
        # assertion dv == add_panel is family-independent and stays).
        old_family = pnl.FILAMENT_REGULARIZATION[]
        pnl.set_filament_regularization!(pnl.VatistasRegularization)
        body = make_dirichlet_diamond_body()
        body.core_size_panel = 5e-2
        body.core_size_targets = 1e-8
        solver = pnl.Backslash(body)
        rt = pnl._initialize_kutta(:steady, body, solver, nothing,
            pnl.RigidTransitionAttachment(), pressure_closure())
        pnl.reset!(body)
        pnl.apply_freestream!((body,), KUTTA_UINF(0.0))
        pnl._kutta_update_edge_lengths!(rt)
        pnl._kutta_capture_snapshot!(rt)
        pnl._reset_kutta_counters!(rt)
        gm = pnl._normalize_grad_mu_options((; basis=:tri); default_basis=:quad)
        c = [0.05]
        pnl._kutta_trial!(rt, c; backend_solve=pnl.DirectBackend(),
            backend_system=pnl.DirectBackend(), grad_mu_options=gm)
        # fixed strengths: reconstruct at c and at 0; the difference is exactly
        # the affine add-on at whatever offset the reconstruction chose
        pnl._kutta_reconstruct_body_velocity!(rt, c;
            backend_system=pnl.DirectBackend(), grad_mu_options=gm)
        v_c = copy(body.velocity)
        pnl._kutta_reconstruct_body_velocity!(rt, [0.0];
            backend_system=pnl.DirectBackend(), grad_mu_options=gm)
        dv = v_c .- body.velocity
        add_panel = zeros(3, body.ncells)
        pnl._add_affine_attached_velocity!(add_panel, body.controlpoints,
            body, c; core_size=body.core_size_panel)
        add_targets = zeros(3, body.ncells)
        pnl._add_affine_attached_velocity!(add_targets, body.controlpoints,
            body, c; core_size=body.core_size_targets)
        # the two offsets genuinely disagree here — the test has teeth
        @test maximum(abs, add_panel .- add_targets) > 1e-6
        @test isapprox(dv, add_panel; rtol=1e-12, atol=1e-13)
        @test !isapprox(dv, add_targets; rtol=1e-12, atol=1e-13)
        pnl.FILAMENT_REGULARIZATION[] = old_family
    end

    # --- §9.4/§11.4: an already-converged base point is accepted ---
    @testset "early-accept gate edge case" begin
        # with wide-open tolerances the warm-start iterate passes the pressure
        # and inner gates immediately and the proposed step is within
        # tolerance: the driver must accept the base point instead of failing
        # the Armijo/step gates on an unimprovable iterate
        body, wake, frames, solver = kutta_case()
        clo = pressure_closure(pnl.SteadyBernoulliProvider(1.2);
            pressure_tolerance=1e3, correction_tolerance=1e3)
        run_kutta!(body, wake, frames, solver;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo)
        diags = pnl.kutta_diagnostics(clo)
        @test all(d -> d.status === :converged, diags)
        @test all(d -> maximum(abs, d.c) == 0, diags)   # base point (c=0) accepted
        @test all(d -> !d.fallback_entered, diags)
    end

    # --- constructor validation (§2.1) ---
    @testset "constructor validation" begin
        provider = pnl.SteadyBernoulliProvider(1.2)
        @test_throws ArgumentError pnl.SteadyBernoulliProvider(0.0)
        @test_throws ArgumentError pnl.SteadyBernoulliProvider(NaN)
        @test_throws ArgumentError pressure_closure(provider; pressure_scale=-1.0)
        @test_throws ArgumentError pressure_closure(provider; correction_scale=0.0)
        @test_throws ArgumentError pressure_closure(provider; pressure_scale=:bogus)
        @test_throws ArgumentError pressure_closure(provider; pressure_tolerance=0.0)
        @test_throws ArgumentError pressure_closure(provider; correction_tolerance=-1e-3)
        @test_throws ArgumentError pressure_closure(provider; on_failure=:silently_continue)
        @test_throws ArgumentError pressure_closure(provider; primary=:not_a_strategy)
        @test_throws ArgumentError pressure_closure(provider; fallback="nope")
        @test_throws ArgumentError pnl.BroydenKutta(; max_iterations=0)
        @test_throws ArgumentError pnl.NewtonKrylovKutta(; fd_relative_step=0.0)
        # explicit scales are accepted
        c = pressure_closure(provider; pressure_scale=2.0, correction_scale=0.5)
        @test c.pressure_scale == 2.0 && c.correction_scale == 0.5
    end

    # --- categories 3, 4, 5, 6, 16, 17: support-boundary validation ---
    @testset "support-boundary validation" begin
        body, wake, frames, solver = kutta_case()
        pc = pressure_closure()
        A = pnl.RigidTransitionAttachment()
        B = pnl.TEAnchoredAttachment()
        J = pnl.JumpKutta()

        # unpaired shedding edge rejected (category 3)
        unpaired = make_dirichlet_diamond_body()
        unpaired.shedding[1][4, 1] = -1
        wu = pnl.PanelWake(unpaired; nwakerows=4)
        @test_throws ArgumentError run_kutta!(unpaired, wu, frames,
            pnl.Backslash(unpaired); wake_attachment=A, kutta_closure=pc)

        # Neumann body rejected (category 4)
        nb = make_plate_vortex_body()
        wn = pnl.PanelWake(nb; nwakerows=4)
        @test_throws ArgumentError run_kutta!(nb, wn,
            pnl.ReferenceFrame(nb; name="vehicle"), pnl.Backslash(nb);
            wake_attachment=A, kutta_closure=pc)

        # non-default formulation rejected (category 4)
        @test_throws "VelocityThroughSources" run_kutta!(body, wake, frames,
            solver; wake_attachment=A, kutta_closure=pc,
            formulation=pnl.TraceCorrected())

        # deferred solvers rejected with the deferral message (16/17)
        ks = pnl.KrylovSolver(body; backend=pnl.DirectBackend())
        @test_throws "deferred" run_kutta!(body, wake, frames, ks;
            wake_attachment=A, kutta_closure=pc)
        fgs_nothing = pnl.FGSSolver(body; build_fgs=false)
        @test_throws "deferred" run_kutta!(body, wake, frames, fgs_nothing;
            wake_attachment=A, kutta_closure=pc)

        # bound_strength_rlx != 1 rejected (category 6)
        @test_throws "bound_strength_rlx" run_kutta!(body, wake, frames,
            solver; wake_attachment=A, kutta_closure=pc, bound_strength_rlx=0.9)

        # multiple bodies rejected
        mb1 = make_dirichlet_diamond_body()
        mb2 = make_dirichlet_diamond_body()
        @test_throws "bodies" pnl.simulate!((mb1, mb2),
            (pnl.PanelWake(mb1; nwakerows=4), pnl.PanelWake(mb2; nwakerows=4)),
            frames, KUTTA_MANEUVER, KUTTA_UINF, KUTTA_TRANGE;
            body_solvers=(pnl.Backslash(mb1), pnl.Backslash(mb2)),
            backend=pnl.DirectBackend(), path=nothing,
            grad_mu_options=(; basis=:tri),
            wake_attachment=A, kutta_closure=pc)

        # a wakeless simulate! run is rejected (no PanelWake to own the block)
        @test_throws "wake count" run_kutta!(body, nothing, frames, solver;
            wake_attachment=A, kutta_closure=pc)

        # Route A rejects a per-step Das refresh: the attachment operator and
        # factorization are assembled once per run
        @test_throws "set_Das_refresh" run_kutta!(body, wake, frames, solver;
            wake_attachment=A, kutta_closure=pc, set_Das_refresh=true)

        # Route B rejects every set_Das_* option (category 5)
        for kw in ((; set_Das_eta_kinematic=0.2),
                   (; set_Das_eta_freestream=0.3),
                   (; set_Das_min_kinematic_displacement=0.01),
                   (; set_Das_kinematic_arc=false),
                   (; set_Das_refresh=true))
            @test_throws "Das" run_kutta!(body, wake, frames, solver;
                wake_attachment=B, kutta_closure=J, kw...)
        end

        # PanelParticleWake remains legacy-only
        ppw = pnl.PanelParticleWake(body; nwakerows=2, max_particles=100)
        @test_throws "PanelParticleWake" run_kutta!(body, ppw, frames, solver;
            wake_attachment=B, kutta_closure=J)

        # rejection happens BEFORE any state mutation: a body/wake pair with
        # settled nonzero state comes through an invalid configuration attempt
        # untouched
        bv, wv, fv, sv = kutta_case()
        run_kutta!(bv, wv, fv, sv)                       # legacy run: nonzero state
        strength_v = copy(bv.strength)
        velocity_v = copy(bv.velocity)
        Das_v = deepcopy(bv.Das)
        wnodes_v = copy(wv.nodes[1])
        wstrength_v = copy(wv.strength[1])
        nwakes_v = wv.nwakes[]
        @test_throws "bound_strength_rlx" run_kutta!(bv, wv, fv, sv;
            wake_attachment=A, kutta_closure=pc, bound_strength_rlx=0.9)
        @test bv.strength == strength_v
        @test bv.velocity == velocity_v
        @test all(D == Ds for (D, Ds) in zip(bv.Das, Das_v))
        @test wv.nodes[1] == wnodes_v
        @test wv.strength[1] == wstrength_v
        @test wv.nwakes[] == nwakes_v

        # semi-infinite wake rejected
        semi = make_dirichlet_diamond_body()
        semibody = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.VortexRing}}(
            copy(semi.nodes), copy(semi.cells), deepcopy(semi.shedding);
            check_mesh=false, watertight=false, semiinfinite_wake=true)
        for i in eachindex(semibody.Das)
            semibody.Das[i] .= repeat([1.0, 0.0, 0.0], 1, size(semibody.Das[i], 2))
        end
        pnl.calc_normals!(semibody); pnl.calc_controlpoints!(semibody)
        ws = pnl.PanelWake(semibody; nwakerows=4)
        @test_throws ArgumentError run_kutta!(semibody, ws,
            pnl.ReferenceFrame(semibody; name="vehicle"),
            pnl.Backslash(semibody); wake_attachment=A, kutta_closure=pc)
    end

    # --- category 12: repeated-trial determinism ---
    @testset "repeated-trial determinism" begin
        body = make_dirichlet_diamond_body()
        solver = pnl.Backslash(body)
        rt = pnl._initialize_kutta(:steady, body, solver, nothing,
            pnl.RigidTransitionAttachment(), pressure_closure())
        pnl.reset!(body)
        pnl.apply_freestream!((body,), KUTTA_UINF(0.0))
        pnl._kutta_update_edge_lengths!(rt)
        pnl._kutta_capture_snapshot!(rt)
        pnl._reset_kutta_counters!(rt)
        gm = pnl._normalize_grad_mu_options((; basis=:tri); default_basis=:quad)
        tk = (; backend_solve=pnl.DirectBackend(),
            backend_system=pnl.DirectBackend(), grad_mu_options=gm)
        c = [0.03]
        t1 = pnl._kutta_trial!(rt, c; tk...)
        r1 = copy(t1.residual); mu1 = copy(t1.mu)
        t2 = pnl._kutta_trial!(rt, c; tk...)
        @test t2.residual == r1
        @test t2.mu == mu1
        # trials leave the snapshot state restorable: a c=0 trial equals the
        # jump solution regardless of trial history
        t0a = pnl._kutta_trial!(rt, [0.0]; tk...)
        r0 = copy(t0a.residual)
        pnl._kutta_trial!(rt, [0.5]; tk...)
        t0b = pnl._kutta_trial!(rt, [0.0]; tk...)
        @test t0b.residual == r0
    end

    # --- categories 10, 11: failure rollback and explicit jump fallback ---
    @testset "failure rollback (:error)" begin
        body, wake, frames, solver = kutta_case()
        strength0 = copy(body.strength)
        wstrength0 = copy(wake.strength[1])
        clo = pressure_closure(KuttaNaNProvider())
        @test_throws pnl.KuttaConvergenceError run_kutta!(body, wake, frames,
            solver; t_range=KUTTA_TRANGE[1:2],
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo)
        # pre-step state fully restored before the throw (wake NODES are
        # legitimately written by the pre-solve update_TE! stage before the
        # snapshot is captured, so they are not asserted here; the round-trip
        # testset below covers node restoration)
        @test body.strength == strength0
        @test wake.strength[1] == wstrength0
        @test !body.wake_correction_active[]
    end

    @testset "failure rollback at a settled step (:error)" begin
        # the provider fails only from the third step on, so the failing
        # step's restored snapshot carries genuinely nonzero committed state —
        # this distinguishes "restored" from "never written"
        body, wake, frames, solver = kutta_case()
        provider = KuttaLateFailProvider(1.2, 2)
        clo = pressure_closure(provider)
        @test_throws pnl.KuttaConvergenceError run_kutta!(body, wake, frames,
            solver; wake_attachment=pnl.RigidTransitionAttachment(),
            kutta_closure=clo)
        @test provider.commits == 2                # steps 1-2 committed, step 3 threw
        @test all(isfinite, body.strength)
        @test any(!iszero, body.strength)          # step-2 accepted state survives
        @test body.wake_correction_active[]        # step-2 committed correction intact
        @test any(!iszero, body.wake_strength_shift)
        @test any(!iszero, wake.strength[1])       # committed wake rows intact
        @test all(isfinite, wake.strength[1])
    end

    @testset "commit failure is deterministic (:commit)" begin
        body, wake, frames, solver = kutta_case()
        strength0 = copy(body.strength)
        clo = pressure_closure(KuttaCommitBombProvider(1.2))
        err = try
            run_kutta!(body, wake, frames, solver;
                t_range=KUTTA_TRANGE[1:2],
                wake_attachment=pnl.RigidTransitionAttachment(),
                kutta_closure=clo)
            nothing
        catch e
            e
        end
        @test err isa pnl.KuttaConvergenceError
        @test err.cause === :commit
        # the half-installed accepted state was rolled back
        @test body.strength == strength0
        @test !body.wake_correction_active[]
    end

    @testset "snapshot restore round trip (nonzero state)" begin
        # every KuttaSnapshot field, asserted against a settled nonzero state
        body, wake, frames, solver = kutta_case()
        clo = pressure_closure()
        run_kutta!(body, wake, frames, solver;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo)
        rt = pnl._initialize_kutta(:simulate, body, solver, wake,
            pnl.RigidTransitionAttachment(), clo)
        pnl._kutta_capture_snapshot!(rt)
        saved_strength = copy(body.strength)
        saved_potential = copy(body.potential)
        saved_velocity = copy(body.velocity)
        saved_Das = deepcopy(body.Das)
        saved_corr = deepcopy(body.wake_strength_correction)
        saved_shift = copy(body.wake_strength_shift)
        saved_active = body.wake_correction_active[]
        saved_wnodes = deepcopy(wake.nodes)
        saved_wstrength = deepcopy(wake.strength)
        saved_wvelocity = deepcopy(wake.velocity)
        saved_nwakes = wake.nwakes[]
        saved_overflowed = wake.overflowed[]
        saved_live_rows = wake.live_rows[]
        saved_live_step = wake.live_step_id[]
        @test any(!iszero, saved_strength)         # the state is genuinely nonzero
        # scribble on every restorable field
        body.strength .= 123.0; body.potential .= -7.0; body.velocity .= 3.0
        for D in body.Das; D .= 0.5; end
        for c in body.wake_strength_correction; c .= 9.0; end
        body.wake_strength_shift .= 4.0
        body.wake_correction_active[] = !saved_active
        for i in eachindex(wake.nodes)
            wake.nodes[i] .= 1.0; wake.strength[i] .= 2.0; wake.velocity[i] .= 3.0
        end
        wake.nwakes[] = saved_nwakes + 1
        wake.overflowed[] = !saved_overflowed
        wake.live_rows[] = 5; wake.live_step_id[] = 99
        pnl._kutta_restore_snapshot!(rt)
        @test body.strength == saved_strength
        @test body.potential == saved_potential
        @test body.velocity == saved_velocity
        @test all(D == Ds for (D, Ds) in zip(body.Das, saved_Das))
        @test all(c == cs for (c, cs) in
            zip(body.wake_strength_correction, saved_corr))
        @test body.wake_strength_shift == saved_shift
        @test body.wake_correction_active[] == saved_active
        @test all(wake.nodes[i] == saved_wnodes[i] for i in eachindex(wake.nodes))
        @test all(wake.strength[i] == saved_wstrength[i] for i in eachindex(wake.strength))
        @test all(wake.velocity[i] == saved_wvelocity[i] for i in eachindex(wake.velocity))
        @test wake.nwakes[] == saved_nwakes
        @test wake.overflowed[] == saved_overflowed
        @test wake.live_rows[] == saved_live_rows
        @test wake.live_step_id[] == saved_live_step
    end

    @testset "explicit jump fallback (:jump)" begin
        # A/pressure with an always-failing provider and on_failure=:jump must
        # commit fresh c = 0 solves — bitwise the legacy A/jump trajectory
        b1, w1, f1, s1 = kutta_case()
        run_kutta!(b1, w1, f1, s1)
        b2, w2, f2, s2 = kutta_case()
        clo = pressure_closure(KuttaNaNProvider(); on_failure=:jump)
        @test_logs (:warn, r"committed a fresh jump") match_mode=:any begin
            run_kutta!(b2, w2, f2, s2;
                wake_attachment=pnl.RigidTransitionAttachment(),
                kutta_closure=clo)
        end
        @test b1.strength == b2.strength           # bitwise the legacy trajectory
        @test w1.strength[1] == w2.strength[1]
        diags = pnl.kutta_diagnostics(clo)
        @test length(diags) == length(KUTTA_TRANGE)
        @test all(d -> d.status === :jump_fallback, diags)
        @test all(d -> d.disposition === :jump_fallback, diags)
    end

    # --- one-time provider commit per accepted step ---
    @testset "one commit per accepted step" begin
        body, wake, frames, solver = kutta_case()
        provider = KuttaCountingProvider(1.2)
        clo = pressure_closure(provider)
        run_kutta!(body, wake, frames, solver;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo)
        @test provider.commits == length(KUTTA_TRANGE)
        @test provider.evals > provider.commits      # trials outnumber commits
        @test provider.rollbacks >= provider.evals   # rollback before every trial
    end

    # --- categories 13, 14: fallback path and simultaneous gates ---
    @testset "Broyden-to-Newton-Krylov fallback and gates" begin
        body, wake, frames, solver = kutta_case()
        # a primary that cannot converge (one iteration, no restarts) forces
        # the fallback from the best finite trial
        clo = pressure_closure(pnl.SteadyBernoulliProvider(1.2);
            primary=pnl.BroydenKutta(; max_iterations=1, max_restarts=0))
        run_kutta!(body, wake, frames, solver;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo)
        diags = pnl.kutta_diagnostics(clo)
        @test all(d -> d.status === :converged, diags)
        @test all(d -> d.fallback_entered, diags)
        @test all(d -> d.method === :newton_krylov, diags)
        # simultaneous gates on every accepted record (category 14)
        for d in diags
            @test d.r_inf_scaled <= clo.pressure_tolerance
            @test d.dc_inf_scaled <= clo.correction_tolerance
            @test d.inner_status === :converged
        end
        # the fallback and the healthy primary agree on the accepted correction
        b2, w2, f2, s2 = kutta_case()
        clo2 = pressure_closure()
        run_kutta!(b2, w2, f2, s2;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo2)
        d1 = pnl.kutta_diagnostics(clo)
        d2 = pnl.kutta_diagnostics(clo2)
        # the corrections are nonzero (leverage testset), so the atol below is
        # a genuine relative check, not vacuously satisfied by c ≈ 0
        @test all(d -> maximum(abs, d.c) > 1e-3*d.correction_scale, d2)
        @test all(isapprox(a.c, b.c; atol=1e-6) for (a, b) in zip(d1, d2))
        @test all(d -> d.method === :broyden && !d.fallback_entered, d2)
    end

    # --- category 15: frozen scales and degenerate-scale rejection ---
    @testset "scale freezing and degenerate scales" begin
        body, wake, frames, solver = kutta_case()
        clo = pressure_closure(pnl.SteadyBernoulliProvider(1.2);
            pressure_scale=3.0, correction_scale=0.7)
        run_kutta!(body, wake, frames, solver;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo)
        for d in pnl.kutta_diagnostics(clo)
            @test d.pressure_scale == 3.0
            @test d.correction_scale == 0.7
        end
        # automatic scales are recorded and positive
        b2, w2, f2, s2 = kutta_case()
        clo2 = pressure_closure()
        run_kutta!(b2, w2, f2, s2;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo2)
        for d in pnl.kutta_diagnostics(clo2)
            @test d.pressure_scale > 0 && d.correction_scale > 0
            @test d.U_s > 0 && d.L_s > 0
            @test d.pressure_scale ≈ 0.5*1.2*d.U_s^2
            @test d.correction_scale ≈ d.U_s*d.L_s
        end
        # zero-velocity state has a degenerate automatic pressure scale
        b3 = make_dirichlet_diamond_body()
        f3 = pnl.ReferenceFrame(b3; name="vehicle")
        @test_throws ArgumentError pnl.steady!((b3,), f3, [0.0, 0.0, 0.0];
            body_solvers=(pnl.Backslash(b3),), backend=pnl.DirectBackend(),
            grad_mu_options=(; basis=:tri),
            wake_attachment=pnl.RigidTransitionAttachment(),
            kutta_closure=pressure_closure())
    end

    # --- affine add-on ≡ production attached-strip influence (pins signs) ---
    # NOTE: influence! never applies wake_strength_shift, under ANY backend —
    # even DirectBackend routes through the FastMultipole source buffers,
    # whose _induced_wake variant has no shift. The reference here is the
    # production full-minus-suppressed unit-strip influence (the velocity
    # analogue of _assemble_W!'s construction), weighted by the same ∓c/2
    # per-panel decomposition set_wake_correction! installs.
    @testset "affine add-on matches production strip influence" begin
        body = make_dirichlet_diamond_body()
        pnl._set_core_sizes!((body,), :core_size_targets)
        influence_kwargs = (; precalc=false, scalar_potential=false,
            velocity=true, velocity_gradient=(false,),
            direct_conditioning=pnl._self_panel_core_size_conditioning())
        c = [0.07]

        v_ref = zeros(3, body.ncells)
        old_strength = copy(body.strength)
        for (i_panel, w) in pnl._kutta_affine_strips(body, c)
            body.strength .= 0.0
            body.strength[i_panel, 2] = 1.0
            body.velocity .= 0.0
            pnl.influence!((body,), (body,), pnl.DirectBackend(); influence_kwargs...)
            v_full = copy(body.velocity)
            body.velocity .= 0.0
            body.suppress_attached_wake[] = true
            pnl.influence!((body,), (body,), pnl.DirectBackend(); influence_kwargs...)
            body.suppress_attached_wake[] = false
            v_ref .+= w .* (v_full .- body.velocity)
        end
        body.strength .= old_strength

        v_addon = zeros(3, body.ncells)
        pnl._add_affine_attached_velocity!(v_addon, body.controlpoints,
            body, c; core_size=body.core_size_panel)
        @test isapprox(v_addon, v_ref; rtol=1e-11, atol=1e-13)
        @test _KLA.norm(v_addon) > 1e-4    # the term is genuinely nonzero

        # probe-node method agrees with the matrix method on the same points
        wake = pnl.PanelWake(body; nwakerows=3)
        wake.nwakes[] = 2
        for r in 1:3, e in axes(wake.nodes[1], 3)
            wake.nodes[1][:, r, e] .= [1.0 + 0.4r, 0.4e, 0.2r]
        end
        pnl.reset!(wake)
        pnl._add_affine_attached_velocity!(wake, body, c;
            core_size=body.core_size_targets)
        pts = reshape(view(wake.nodes[1], :, 1:3, :), 3, :)
        v_pts = zeros(3, size(pts, 2))
        pnl._add_affine_attached_velocity!(v_pts, pts, body, c;
            core_size=body.core_size_targets)
        @test isapprox(reshape(view(wake.velocity[1], :, 1:3, :), 3, :), v_pts;
            rtol=1e-12, atol=1e-14)
    end

    # --- category 18: FMM backends agree with direct (solve and full) ---
    @testset "direct/FMM backend agreement" begin
        b1, w1, f1, s1 = kutta_case()
        clo1 = pressure_closure()
        run_kutta!(b1, w1, f1, s1;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo1)
        # FMM solve backend only
        b2, w2, f2, s2 = kutta_case()
        clo2 = pressure_closure()
        run_kutta!(b2, w2, f2, s2;
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo2,
            backend=pnl.DirectBackend(),
            backend_solve=pnl.FastMultipoleBackend())
        @test isapprox(b1.strength, b2.strength; rtol=1e-6, atol=1e-9)
        # full FMM (all backend slots) — enabled by the direct affine add-on
        for (att, label) in ((pnl.RigidTransitionAttachment(), "A"),
                             (pnl.TEAnchoredAttachment(), "B"))
            bd, wd, fd, sd = kutta_case()
            clod = pressure_closure()
            run_kutta!(bd, wd, fd, sd; wake_attachment=att, kutta_closure=clod)
            bf, wf, ff, sf = kutta_case()
            clof = pressure_closure()
            run_kutta!(bf, wf, ff, sf; wake_attachment=att, kutta_closure=clof,
                backend=pnl.FastMultipoleBackend())
            @test isapprox(bd.strength, bf.strength; rtol=1e-6, atol=1e-9)
            @test isapprox(wd.strength[1], wf.strength[1]; rtol=1e-6, atol=1e-9)
            df = pnl.kutta_diagnostics(clof)
            @test all(d -> d.status in (:converged, :startup_jump), df)
        end
    end

    # --- category 19: steady! support boundary ---
    @testset "steady! Route A pressure and TE-anchored rejection" begin
        body = make_dirichlet_diamond_body()
        frames = pnl.ReferenceFrame(body; name="vehicle")
        clo = pressure_closure()
        pnl.steady!((body,), frames, [1.0, 0.0, 0.3];
            body_solvers=(pnl.Backslash(body),), backend=pnl.DirectBackend(),
            grad_mu_options=(; basis=:tri),
            wake_attachment=pnl.RigidTransitionAttachment(), kutta_closure=clo)
        @test all(isfinite, body.strength)
        diags = pnl.kutta_diagnostics(clo)
        @test length(diags) == 1 && diags[1].status === :converged
        # steady! rejects the TE-anchored route
        b2 = make_dirichlet_diamond_body()
        @test_throws "rejects TEAnchoredAttachment" pnl.steady!((b2,),
            pnl.ReferenceFrame(b2; name="vehicle"), [1.0, 0.0, 0.3];
            body_solvers=(pnl.Backslash(b2),), backend=pnl.DirectBackend(),
            grad_mu_options=(; basis=:tri),
            wake_attachment=pnl.TEAnchoredAttachment(),
            kutta_closure=pnl.JumpKutta())
    end

    # --- category 22: concrete typing of every new struct ---
    @testset "concrete typing" begin
        body = make_dirichlet_diamond_body()
        wake = pnl.PanelWake(body; nwakerows=4)
        solver = pnl.Backslash(body)
        clo = pressure_closure()
        rt = pnl._initialize_kutta(:simulate, body, solver, wake,
            pnl.RigidTransitionAttachment(), clo)
        diag = pnl._kutta_diagnostics_record(rt, 0, :converged, :converged)
        err = pnl.KuttaConvergenceError(diag, :nonconvergence, "msg")
        view = pnl.KuttaTrialView{Float64, typeof(body)}(body, rt.edges,
            body.velocity)
        instances = Any[
            pnl.RigidTransitionAttachment(), pnl.TEAnchoredAttachment(),
            pnl.JumpKutta(), pnl.BroydenKutta(), pnl.NewtonKrylovKutta(),
            pnl.SteadyBernoulliProvider(1.2), clo, diag, err, view,
            rt, rt.snapshot, rt.trial, rt.edges,
        ]
        for x in instances
            T = typeof(x)
            @test isconcretetype(T)
            for (fname, ftype) in zip(fieldnames(T), fieldtypes(T))
                @test isconcretetype(ftype) || ftype === Union{}
                @test ftype !== Any
                @test !(ftype isa UnionAll)
                @test ftype !== Function
            end
        end
        # the abstract API roots remain abstract by design
        @test isabstracttype(pnl.AbstractWakeAttachment)
        @test isabstracttype(pnl.AbstractKuttaClosure)
        @test isabstracttype(pnl.AbstractKuttaPressureProvider)
    end
end
