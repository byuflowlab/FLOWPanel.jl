#=##############################################################################
Lightweight verification of the swappable wake→body solve formulations
(`VelocityThroughSources`, `GreenReconstruction`, `TraceCorrected`).

Frozen small pitching-wing state, DirectBackend throughout. The substantive
physical validation lives in debug/dirichlet_solve (Tasks 4/5); this script
checks the mechanical contracts:

1. default passthrough is bit-identical and validation errors fire
2. the wake-strength correction channel: γ = C·μ̃ − c, exterior potential
   shift, operator-mode guard in G assembly
3. W columns match the attached-wake block of G; TraceCorrected with c = 0 is
   bit-identical to production; direct residual with prescribed c
4. TraceCorrected(:line_integral) end-to-end: c vs the trailing-edge wake
   scalar-potential oracle; first shed row equals C·μ̃ − c
5. Green machinery: (I−B)q = Sσ residual, gauge, constant eigenmode, gauge
   invariance of C·q
6. GreenReconstruction vs TraceCorrected(:green) exact consistency identity

Run: julia --project test/formulation_test.jl
=###############################################################################

using Test
import LinearAlgebra as LA
import LinearAlgebra: norm, dot, ldiv!
using FastMultipole.StaticArrays: SVector

import FLOWPanel as pnl

include(joinpath(pnl.examples_path, "pitching_wing.jl"))

# ------------------------------------------------------------------ setup ----
const C_CHORD = 1.0
const B_SPAN = 4.0
const UMAG = 100.0
const ALPHA = 5.0 * pi / 180
const UINF = SVector{3}(UMAG*cos(ALPHA), 0.0, UMAG*sin(ALPHA))
const DIRECT = pnl.DirectBackend()

function small_body(; semiinfinite_wake=false, das_length_c=0.05)
    body = build_pitching_wing_body(C_CHORD, B_SPAN;
        n_span=7, n_airfoil=41, n_endcap=5, endcap=:round, semiinfinite_wake)
    set_wake_Das!(body, UINF; magnitude=das_length_c*C_CHORD)
    return body
end

"Fabricated flat prescribed wake with per-strip strength `gamma`."
function flat_wake(body, gamma; nfree=8, free_row_length_c=0.5)
    wake = pnl.PanelWake(body; nwakerows=nfree, include_final_filament=false)
    pnl.update_TE!(wake, body)
    Das = view(body.Das[1], :, 1)
    direction = Das ./ norm(Das) .* (free_row_length_c*C_CHORD)
    for nodes in wake.nodes
        first_row = copy(view(nodes, :, 1, :))
        for row in 1:nfree+1
            view(nodes, :, row, :) .= first_row .+ (row - 1) .* direction
        end
    end
    for strength in wake.strength
        for row in 1:nfree
            view(strength, 1, row, :) .= gamma
        end
    end
    wake.nwakes[] = nfree
    return wake
end

"Sparse-free Kutta map application."
kutta_map(body, x) = [x[shed[1, i]] - (shed[4, i] > 0 ? x[shed[4, i]] : 0.0)
                      for shed in body.shedding for i in axes(shed, 2)]

nedges(body) = sum(size(shed, 2) for shed in body.shedding; init=0)

function run_aero!(body, wake, solver; formulation=pnl.VelocityThroughSources(),
        formulation_state=nothing, i_step=0)
    frames = pnl.ReferenceFrame(body)
    pnl._steady_aerodynamics!(body, (body,), (wake,), frames, UINF, solver;
        backend_wake=DIRECT, backend_solve=DIRECT, backend_system=DIRECT,
        update_trailing_edges=false,
        formulation, formulation_state, i_step)
    return body
end

function lu_product(factorization, x)
    permuted = factorization.L * (factorization.U * x)
    product = similar(permuted)
    product[factorization.p] .= permuted
    return product
end

body0 = small_body()
solver0 = pnl.Backslash(body0)
M = nedges(body0)
N = body0.ncells
println("test body: $N cells, $M shedding edges")

# a mild fabricated strip circulation for prescribed wakes
gamma0 = [0.1*UMAG*C_CHORD*sin(pi*(i - 0.5)/M) for i in 1:M]

# ------------------------------------------------------------------------------
@testset "Stage 1: default passthrough + validation" begin
    wake = flat_wake(body0, gamma0)

    bodyA = deepcopy(body0)
    run_aero!(bodyA, deepcopy(wake), solver0)
    bodyB = deepcopy(body0)
    run_aero!(bodyB, deepcopy(wake), solver0;
        formulation=pnl.VelocityThroughSources(), formulation_state=nothing)
    @test bodyA.strength == bodyB.strength
    @test bodyA.velocity == bodyB.velocity

    # validation errors
    semiinf = small_body(; semiinfinite_wake=true, das_length_c=1.0)
    @test_throws ErrorException pnl.initialize_formulation(
        pnl.TraceCorrected(), (semiinf,), (wake,), solver0,
        DIRECT, DIRECT)
    fmm = pnl.FastMultipoleBackend()
    @test_throws ErrorException pnl.initialize_formulation(
        pnl.TraceCorrected(), (body0,), (wake,), solver0, DIRECT, fmm)
    @test_throws ErrorException pnl.initialize_formulation(
        pnl.GreenReconstruction(), (body0,), (nothing,), solver0, DIRECT, DIRECT)
    @test_throws ErrorException pnl.TraceCorrected(estimator=:bogus)
    @test_throws ErrorException pnl.GreenReconstruction(gauge=:bogus)
end

# ------------------------------------------------------------------------------
@testset "Stage 2: correction channel + operator guard" begin
    body = deepcopy(body0)
    mu = randn(N)
    body.strength[:, 2] .= mu
    gamma_plain = [pnl._get_wakestrength_Gamma(body, i, s)
                   for (s, shed) in enumerate(body.shedding)
                   for i in axes(shed, 2)]
    @test gamma_plain ≈ kutta_map(body, mu) atol=1e-14

    c_syn = 0.01 .* randn(M)
    pnl.set_wake_correction!(body, c_syn)
    gamma_corr = [pnl._get_wakestrength_Gamma(body, i, s)
                  for (s, shed) in enumerate(body.shedding)
                  for i in axes(shed, 2)]
    @test gamma_corr ≈ gamma_plain .- c_syn atol=1e-13

    # exterior potential shifts by −Σ c_e φ_e with φ_e the unit attached-panel
    # potential (difference of induced with the attached wake on/off)
    probe = SVector{3}(2.5*C_CHORD, 0.3*B_SPAN, 0.8*C_CHORD)
    switch = pnl.FastMultipole.DerivativesSwitch(true, false, false)
    function total_potential(bdy)
        phi = 0.0
        for j in 1:bdy.ncells
            p, _, _ = pnl.induced(probe, bdy, j, switch;
                kerneloffset=bdy.kerneloffset_panel)
            phi += p
        end
        return phi
    end
    phi_corr = total_potential(body)
    pnl.clear_wake_correction!(body)
    phi_plain = total_potential(body)

    phi_unit = zeros(M)
    saved = copy(body.strength)
    body.strength .= 0.0
    e = 0
    for (s, shed) in enumerate(body.shedding), i in axes(shed, 2)
        e += 1
        up = shed[1, i]
        body.strength[up, 2] = 1.0
        p_full, _, _ = pnl.induced(probe, body, up, switch;
            kerneloffset=body.kerneloffset_panel)
        body.suppress_attached_wake[] = true
        p_body, _, _ = pnl.induced(probe, body, up, switch;
            kerneloffset=body.kerneloffset_panel)
        body.suppress_attached_wake[] = false
        phi_unit[e] = p_full - p_body
        body.strength[up, 2] = 0.0
    end
    body.strength .= saved
    @test phi_corr - phi_plain ≈ -dot(c_syn, phi_unit) rtol=1e-10

    # operator guard: G assembly identical with correction active
    pnl.set_wake_correction!(body, c_syn)
    G_corr = zeros(N, N)
    pnl._G!(G_corr, body, body; kerneloffset=body.kerneloffset_panel)
    @test body.wake_correction_active[]      # restored after assembly
    pnl.clear_wake_correction!(body)
    G_plain = zeros(N, N)
    pnl._G!(G_plain, body, body; kerneloffset=body.kerneloffset_panel)
    @test G_corr == G_plain
end

# ------------------------------------------------------------------------------
@testset "Stage 3: W assembly + Wc RHS injection" begin
    body = deepcopy(body0)
    edges = pnl._shedding_edge_map(body)
    W = zeros(N, M)
    pnl._assemble_W!(W, body, edges)
    @test all(isfinite, W)

    # W columns == attached-wake block of G (full minus wake-suppressed)
    G_full = zeros(N, N); B_only = zeros(N, N)
    pnl._G!(G_full, body, body; kerneloffset=body.kerneloffset_panel)
    pnl._assemble_B!(B_only, body)
    dG = G_full - B_only
    for e in 1:M
        @test W[:, e] ≈ dG[:, edges.upper[e]] atol=1e-13
    end
    # non-shedding columns carry no wake block
    shed_cols = Set(vcat(edges.upper, filter(>(0), edges.lower)))
    for j in 1:N
        j in shed_cols && continue
        @test maximum(abs, view(dG, :, j)) < 1e-13
    end

    # c = 0 → bit-identical to production on the same frozen state
    wake = flat_wake(body0, zeros(M))
    f = pnl.TraceCorrected(estimator=:line_integral, quadrature=:trapezoid)
    state = pnl.initialize_formulation(f, (deepcopy(body0),), (wake,), solver0,
        DIRECT, DIRECT)
    bodyP = deepcopy(body0)
    run_aero!(bodyP, deepcopy(wake), solver0)
    bodyT = deepcopy(body0)
    run_aero!(bodyT, deepcopy(wake), solver0; formulation=f,
        formulation_state=state)
    @test maximum(abs, state.c) == 0.0
    @test bodyT.strength == bodyP.strength

    # direct residual with a prescribed nonzero c:
    # G·μ̃ = −S(σ0+σ) + W·c
    wake1 = flat_wake(body0, gamma0)
    f2 = pnl.TraceCorrected(estimator=:line_integral, quadrature=:trapezoid,
        recompute_interval=typemax(Int)÷2)
    state2 = pnl.initialize_formulation(f2, (deepcopy(body0),), (wake1,),
        solver0, DIRECT, DIRECT)
    c_rand = 0.05 .* randn(M)
    state2.c .= c_rand
    state2.last_recompute[] = 0                # suppress recompute (lagged c)
    bodyR = deepcopy(body0)
    run_aero!(bodyR, deepcopy(wake1), solver0; formulation=f2,
        formulation_state=state2, i_step=1)
    @test state2.c == c_rand                   # untouched by the lag
    mu = copy(view(bodyR.strength, :, 2))
    sig = copy(view(bodyR.strength, :, 1))     # σ0+σ set by set_strengths!
    Ssig = zeros(N)
    pnl._source_potential!(Ssig, bodyR, sig, DIRECT)
    rhs = -Ssig .+ W*c_rand
    residual = norm(G_full*mu - rhs)/max(norm(rhs), eps())
    @test residual <= 1e-10
    # physical mode active with the affine Kutta strength
    @test bodyR.wake_correction_active[]
    gam = [pnl._get_wakestrength_Gamma(bodyR, i, s)
           for (s, shed) in enumerate(bodyR.shedding) for i in axes(shed, 2)]
    @test gam ≈ kutta_map(bodyR, mu) .- c_rand atol=1e-12
end

# ------------------------------------------------------------------------------
@testset "Stage 4: line-integral estimator end-to-end" begin
    wake = flat_wake(body0, gamma0)

    # oracle: wake scalar-potential difference at the paired TE control points
    edges = pnl._shedding_edge_map(body0)
    probes = pnl.FastMultipole.ProbeSystem(2M, Float64)
    for e in 1:M
        up, lo = edges.upper[e], edges.lower[e]
        probes.position[e]     = SVector{3}(body0.controlpoints[:, up]...)
        probes.position[M + e] = SVector{3}(body0.controlpoints[:, lo]...)
        probes.scalar_potential[e] = 0.0
        probes.scalar_potential[M + e] = 0.0
    end
    pnl.influence!((probes,), pnl.get_sources(wake), DIRECT;
        precalc=true, scalar_potential=true, gradient=false, hessian=false)
    c_oracle = [probes.scalar_potential[e] - probes.scalar_potential[M + e]
                for e in 1:M]

    results = Dict{String, Vector{Float64}}()
    for (label, f) in (
            ("trapezoid", pnl.TraceCorrected(estimator=:line_integral,
                quadrature=:trapezoid)),
            ("simpson", pnl.TraceCorrected(estimator=:line_integral,
                quadrature=:simpson)),
            ("deformed", pnl.TraceCorrected(estimator=:line_integral,
                quadrature=:simpson, interior_path=:deformed)),
            )
        state = pnl.initialize_formulation(f, (deepcopy(body0),), (wake,),
            solver0, DIRECT, DIRECT)
        body = deepcopy(body0)
        run_aero!(body, deepcopy(wake), solver0; formulation=f,
            formulation_state=state)
        @test all(isfinite, state.c)
        results[label] = copy(state.c)
        rel = norm(state.c - c_oracle)/max(norm(c_oracle), eps())
        println("  c($label) vs TE-potential oracle: rel L2 = $(round(rel; sigdigits=3))")

        if label == "trapezoid"
            # first shed row carries the corrected γ = C·μ̃ − c
            mu = copy(view(body.strength, :, 2))
            pnl.shed_wake!(deepcopy(wake), body)  # non-mutating check target
            wake2 = deepcopy(wake)
            pnl.shed_wake!(wake2, body)
            first_row = vec(copy(view(wake2.strength[1], 1, 1, :)))
            @test first_row ≈ kutta_map(body, mu) .- state.c atol=1e-12
        end
    end
    # quadrature refinement should not blow up the estimate; report only
    spread = norm(results["trapezoid"] - results["simpson"])
    println("  trapezoid-vs-simpson spread (quadrature error bar): $(round(spread; sigdigits=3))")
    println("  straight-vs-deformed spread (leakage diagnostic):  $(round(norm(results["simpson"] - results["deformed"]); sigdigits=3))")
end

# ------------------------------------------------------------------------------
@testset "Stage 5: Green machinery" begin
    body = deepcopy(body0)
    wake = flat_wake(body0, gamma0)

    # wake-only σ at centroids from a direct wake-velocity evaluation
    body.velocity .= 0.0
    pnl.influence!((body,), pnl.get_sources(wake), DIRECT;
        precalc=true, scalar_potential=false, velocity=true)
    sigma = [-dot(view(body.velocity, :, i), view(body.normals, :, i))
             for i in 1:N]
    body.velocity .= 0.0

    Ssigma = zeros(N)
    pnl._source_potential!(Ssigma, body, sigma, DIRECT)

    B = zeros(N, N)
    pnl._assemble_B!(B, body)
    a = pnl._panel_areas(body)
    ones_response = B*ones(N)
    println("  ‖B·1 − 1‖∞ (constant eigenmode defect): $(round(maximum(abs, ones_response .- 1); sigdigits=3))")

    gs_am = pnl._build_green_solve_state(body, :area_mean)
    q_am = copy(pnl._green_solve_q!(gs_am, Ssigma))
    gs_ls = pnl._build_green_solve_state(body, :lsq)
    q_ls = copy(pnl._green_solve_q!(gs_ls, Ssigma))

    # bordered system: (I−B)q + λa = Sσ with aᵀq = 0
    lambda = gs_am.sol_b[end]
    res_am = norm((q_am - B*q_am) .+ lambda .* a - Ssigma)/max(norm(Ssigma), eps())
    @test res_am <= 1e-10
    @test abs(dot(a, q_am)) <= 1e-8*norm(q_am)*norm(a)
    println("  gauge multiplier λ (compatibility defect): $(round(lambda; sigdigits=3))")

    # C·q is gauge-invariant on paired edges
    Cq_am = kutta_map(body, q_am)
    Cq_ls = kutta_map(body, q_ls)
    @test norm(Cq_am - Cq_ls) <= 1e-8*max(norm(Cq_am), eps())

    # report (not gate) the discrete Green-identity error vs the directly
    # evaluated wake potential trace
    q_direct = zeros(N)
    old_pot = copy(body.potential)
    body.potential .= 0.0
    pnl.influence!((body,), pnl.get_sources(wake), DIRECT;
        precalc=true, scalar_potential=true, velocity=false)
    q_direct .= body.potential
    body.potential .= old_pot
    q_direct .-= dot(a, q_direct)/sum(a)          # match the area-mean gauge
    println("  ‖q_green − q_direct‖/‖q_direct‖ (P1/P2 representation error): "*
        "$(round(norm(q_am - q_direct)/max(norm(q_direct), eps()); sigdigits=3))")
end

# ------------------------------------------------------------------------------
@testset "Stage 6: GreenReconstruction vs TraceCorrected(:green)" begin
    wake = flat_wake(body0, gamma0)

    fG = pnl.GreenReconstruction()
    stG = pnl.initialize_formulation(fG, (deepcopy(body0),), (wake,), solver0,
        DIRECT, DIRECT)
    bodyG = deepcopy(body0)
    run_aero!(bodyG, deepcopy(wake), solver0; formulation=fG,
        formulation_state=stG)
    muE = copy(view(bodyG.strength, :, 2))
    @test !bodyG.wake_correction_active[]
    @test bodyG.strength[:, 1] ≈ stG.sigma0 atol=0.0   # σ0-only sources

    fT = pnl.TraceCorrected(estimator=:green)
    stT = pnl.initialize_formulation(fT, (deepcopy(body0),), (wake,), solver0,
        DIRECT, DIRECT)
    bodyT = deepcopy(body0)
    run_aero!(bodyT, deepcopy(wake), solver0; formulation=fT,
        formulation_state=stT)
    mut = copy(view(bodyT.strength, :, 2))

    # both formulations reconstruct the same q on the same frozen state
    q = stG.green.q
    @test q ≈ stT.green.q rtol=1e-12
    @test stT.c ≈ kutta_map(body0, q) atol=1e-12

    # exact consistency identity: μ̃ − μE − q = −G⁻¹ r with
    # r = Sσ − (I−B)q (the gauge-row remainder λa)
    body = deepcopy(body0)
    B = zeros(N, N)
    pnl._assemble_B!(B, body)
    Ssigma = zeros(N)
    pnl._source_potential!(Ssigma, body, stG.sigma, DIRECT)
    r = Ssigma - (q - B*q)
    correction = ldiv!(zeros(N), solver0.Glu, copy(r))
    lhs = mut - muE - q
    @test norm(lhs + correction) <= 1e-9*max(norm(mut), norm(muE))

    # physical attached-wake strengths agree up to the same remainder
    gamG = kutta_map(bodyG, muE)
    gamT = [pnl._get_wakestrength_Gamma(bodyT, i, s)
            for (s, shed) in enumerate(bodyT.shedding) for i in axes(shed, 2)]
    gam_gap = norm(gamT - gamG)/max(norm(gamG), eps())
    println("  ‖γ_trace − γ_green‖/‖γ‖ (gauge-remainder gap): $(round(gam_gap; sigdigits=3))")
    @test gam_gap <= 1e-6

    # both linear systems satisfied
    G_full = zeros(N, N)
    pnl._G!(G_full, body, body; kerneloffset=body.kerneloffset_panel)
    Ssigma0 = zeros(N)
    pnl._source_potential!(Ssigma0, body, stG.sigma0, DIRECT)
    resE = norm(G_full*muE + Ssigma0 + q)/max(norm(Ssigma0), eps())
    W = stT.W
    Ssig_tot = zeros(N)
    pnl._source_potential!(Ssig_tot, body, stG.sigma0 .+ stG.sigma, DIRECT)
    resT = norm(G_full*mut + Ssig_tot - W*stT.c)/max(norm(Ssig_tot), eps())
    println("  explicit residual: $(round(resE; sigdigits=3)); trace-corrected residual: $(round(resT; sigdigits=3))")
    @test resE <= 1e-9
    @test resT <= 1e-9
end

# ------------------------------------------------------------------------------
@testset "Stage 7: Green solver routes (Krylov/FGS/Backslash)" begin
    body = deepcopy(body0)
    wake = flat_wake(body0, gamma0)

    # wake-only σ and its source potential, as in Stage 5
    body.velocity .= 0.0
    pnl.influence!((body,), pnl.get_sources(wake), DIRECT;
        precalc=true, scalar_potential=false, velocity=true)
    sigma = [-dot(view(body.velocity, :, i), view(body.normals, :, i))
             for i in 1:N]
    body.velocity .= 0.0
    Ssigma = zeros(N)
    pnl._source_potential!(Ssigma, body, sigma, DIRECT)

    fmm = pnl.FastMultipoleBackend(expansion_order=14,
        multipole_acceptance=0.4, leaf_size=50)

    # (1) B-product: FMM gate agrees with Direct; gate actually removes wake
    x = randn(N)
    Bx_direct = zeros(N); Bx_fmm = zeros(N); Gx_direct = zeros(N)
    pnl._with_green_scratch(body) do
        pnl._green_B_product!(Bx_direct, body, x, DIRECT)
        pnl._green_B_product!(Bx_fmm, body, x, fmm)
    end
    # same product without suppression = G·x (wake included)
    old_strength = copy(body.strength); old_potential = copy(body.potential)
    pnl._green_B_product!(Gx_direct, body, x, DIRECT)
    body.strength .= old_strength; body.potential .= old_potential
    @test norm(Bx_fmm - Bx_direct)/max(norm(Bx_direct), eps()) <= 1e-6
    @test norm(Gx_direct - Bx_direct) > 1e-8 * norm(Bx_direct)  # gate flips wake

    # (2) route agreement on C·q (gauge-invariant)
    gs_lu = pnl._build_green_solve_state(body, :area_mean)
    q_lu = copy(pnl._green_solve_q!(gs_lu, Ssigma))
    Cq_lu = kutta_map(body, q_lu)
    a = pnl._panel_areas(body)

    ks_d = pnl.KrylovSolver(body; method=:gmres, itmax=500, atol=1e-12,
        rtol=1e-12, backend=DIRECT)
    gs_kd = pnl._build_green_state(body, :area_mean, ks_d)
    q_kd = copy(pnl._green_solve_q!(gs_kd, Ssigma, body, DIRECT))
    @test norm(kutta_map(body, q_kd) - Cq_lu) <= 1e-8*max(norm(Cq_lu), eps())
    @test abs(dot(a, q_kd)) <= 1e-8*norm(q_kd)*norm(a)
    # warm-started second solve reproduces the same q
    q_kd2 = copy(pnl._green_solve_q!(gs_kd, Ssigma, body, DIRECT))
    @test norm(q_kd2 - q_kd) <= 1e-8*max(norm(q_kd), eps())

    ks_f = pnl.KrylovSolver(body; method=:gmres, itmax=500, atol=1e-12,
        rtol=1e-12, backend=fmm)
    gs_kf = pnl._build_green_state(body, :area_mean, ks_f)
    q_kf = copy(pnl._green_solve_q!(gs_kf, Ssigma, body, fmm))
    rel_kf = norm(kutta_map(body, q_kf) - Cq_lu)/max(norm(Cq_lu), eps())
    println("  Krylov(FMM) vs LU C·q rel diff: $(round(rel_kf; sigdigits=3))")
    @test rel_kf <= 1e-4

    fgs = pnl.FGSSolver(body; rlx=1.0, tolerance=1e-10, max_iterations=2000,
        calc_cps=false, build_fgs=false)
    gs_p = pnl._build_green_state(body, :area_mean, fgs)
    q_p = copy(pnl._green_solve_q!(gs_p, Ssigma, body, DIRECT))
    rel_p = norm(kutta_map(body, q_p) - Cq_lu)/max(norm(Cq_lu), eps())
    println("  Picard(FGS) vs LU C·q rel diff: $(round(rel_p; sigdigits=3))")
    @test rel_p <= 1e-4
    @test abs(dot(a, q_p)) <= 1e-8*max(norm(q_p)*norm(a), eps())

    # scratch discipline: body state untouched by the iterative solves
    @test body.strength == body0.strength
    @test !body.suppress_attached_wake[]

    # (3) end-to-end: GreenReconstruction through each route
    fL = pnl.GreenReconstruction()
    stL = pnl.initialize_formulation(fL, (deepcopy(body0),), (wake,), solver0,
        DIRECT, DIRECT)
    bodyL = deepcopy(body0)
    run_aero!(bodyL, deepcopy(wake), solver0; formulation=fL,
        formulation_state=stL)
    gamL = kutta_map(bodyL, view(bodyL.strength, :, 2))

    for (label, gsolver) in (("krylov", ks_d), ("fgs", fgs))
        f = pnl.GreenReconstruction(green_solver=gsolver)
        st = pnl.initialize_formulation(f, (deepcopy(body0),), (wake,),
            solver0, DIRECT, DIRECT)
        bodyX = deepcopy(body0)
        run_aero!(bodyX, deepcopy(wake), solver0; formulation=f,
            formulation_state=st)
        gamX = kutta_map(bodyX, view(bodyX.strength, :, 2))
        rel = norm(gamX - gamL)/max(norm(gamL), eps())
        println("  E2E GreenReconstruction($label) vs LU γ rel diff: $(round(rel; sigdigits=3))")
        @test rel <= 1e-4
    end

    # invalid combos
    @test_throws ErrorException pnl.GreenReconstruction(gauge=:lsq,
        green_solver=ks_d)
    @test_throws ErrorException pnl.TraceCorrected(gauge=:lsq,
        green_solver=ks_d)
end

# ------------------------------------------------------------------------------
@testset "Stage 8: DirectWakePotential — Task-3 manual-solve equivalence" begin
    wake = flat_wake(body0, gamma0)

    # independent Task-3 reference solve (task3.md): σ0 = −U∞·n, q_f evaluated
    # directly from the prescribed wake, μE = G_Δ \ (−S·σ0 − q_f)
    body_ref = deepcopy(body0)
    pnl.calc_normals!(body_ref); pnl.calc_controlpoints!(body_ref)
    sigma0 = [-dot(UINF, view(body_ref.normals, :, i)) for i in 1:N]
    Ssig = zeros(N); pnl._source_potential!(Ssig, body_ref, sigma0, DIRECT)
    qf = zeros(N); pnl._wake_potential!(qf, body_ref, (wake,), DIRECT)
    rhs = -Ssig .- qf
    mu_ref = zeros(N); ldiv!(mu_ref, solver0.Glu, rhs)

    # production DirectWakePotential through the full aero step
    f = pnl.DirectWakePotential()
    st = pnl.initialize_formulation(f, (deepcopy(body0),), (wake,), solver0,
        DIRECT, DIRECT)
    bodyD = deepcopy(body0)
    run_aero!(bodyD, deepcopy(wake), solver0; formulation=f,
        formulation_state=st)

    muD = view(bodyD.strength, :, 2)
    @test norm(muD - mu_ref) <= 1e-10*max(norm(mu_ref), eps())
    # sources carry σ0 only
    @test norm(view(bodyD.strength, :, 1) - sigma0) <= 1e-12*max(norm(sigma0), eps())
    # single-shot linear residual gate ≤1e-10 (task3 documented tolerance)
    resid = lu_product(solver0.Glu, collect(muD)) - rhs
    rel_resid = norm(resid)/max(norm(rhs), eps())
    println("  DirectWakePotential single-shot rel residual: $(round(rel_resid; sigdigits=3))")
    @test rel_resid <= 1e-10
    @test all(isfinite, bodyD.strength)

    # state does not corrupt wake geometry/strengths
    wref = flat_wake(body0, gamma0)
    @test all(isequal(a, b) for (a, b) in zip(wake.nodes, wref.nodes))
    @test all(isequal(a, b) for (a, b) in zip(wake.strength, wref.strength))

    # recompute_interval determinism: a fresh state reproduces the same μE
    st2 = pnl.initialize_formulation(f, (deepcopy(body0),), (wake,), solver0,
        DIRECT, DIRECT)
    bodyD2 = deepcopy(body0)
    run_aero!(bodyD2, deepcopy(wake), solver0; formulation=f,
        formulation_state=st2)
    @test view(bodyD2.strength, :, 2) == muD

    # no affine wake correction is left active
    @test !bodyD.wake_correction_active[]
end

# ------------------------------------------------------------------------------
@testset "Stage 9: DirectWakePotential invalid-config rejection" begin
    wake = flat_wake(body0, gamma0)
    f = pnl.DirectWakePotential()

    # include_final_filament=true (trailing filament is vector-potential-only)
    wake_ff = pnl.PanelWake(body0; nwakerows=8)   # default include_final_filament=true
    wake_ff.nwakes[] = 4
    @test_throws ErrorException pnl.initialize_formulation(
        f, (deepcopy(body0),), (wake_ff,), solver0, DIRECT, DIRECT)

    # particle wake (mixed / vector-potential-only sources)
    ppw = pnl.PanelParticleWake(body0; nwakerows=3, max_particles=64,
        include_final_filament=false)
    @test_throws ErrorException pnl.initialize_formulation(
        f, (deepcopy(body0),), (ppw,), solver0, DIRECT, DIRECT)

    # semi-infinite body
    semiinf = small_body(; semiinfinite_wake=true, das_length_c=1.0)
    @test_throws ErrorException pnl.initialize_formulation(
        f, (semiinf,), (wake,), solver0, DIRECT, DIRECT)

    # no wake
    @test_throws ErrorException pnl.initialize_formulation(
        f, (deepcopy(body0),), (nothing,), solver0, DIRECT, DIRECT)

    # non-Backslash solver
    fgs = pnl.FGSSolver(body0; rlx=1.0, tolerance=1e-8, max_iterations=10,
        calc_cps=false, build_fgs=false)
    @test_throws ErrorException pnl.initialize_formulation(
        f, (deepcopy(body0),), (wake,), fgs, DIRECT, DIRECT)

    # multi-body
    @test_throws ErrorException pnl.initialize_formulation(
        f, (deepcopy(body0), deepcopy(body0)), (wake, wake), solver0,
        DIRECT, DIRECT)

    # constructor rejects recompute_interval < 1
    @test_throws ErrorException pnl.DirectWakePotential(recompute_interval=0)

    # FMM backend is accepted (no DirectBackend guard)
    fmm = pnl.FastMultipoleBackend()
    st_fmm = pnl.initialize_formulation(f, (deepcopy(body0),), (wake,), solver0,
        DIRECT, fmm)
    @test st_fmm isa pnl.DirectWakePotentialState
end

# ------------------------------------------------------------------------------
@testset "Stage 10: ConstantDoublet vs VortexRing frozen wake" begin
    nfree = 8
    wake_cd = flat_wake(body0, gamma0; nfree=nfree)          # ConstantDoublet
    @test wake_cd isa pnl.PanelWake{pnl.ConstantDoublet}

    # identical geometry/strength VortexRing wake
    wake_vr = pnl.PanelWake(body0, pnl.VortexRing; nwakerows=nfree,
        include_final_filament=false)
    pnl.update_TE!(wake_vr, body0)
    for (nvr, ncd) in zip(wake_vr.nodes, wake_cd.nodes)
        nvr .= ncd
    end
    for (svr, scd) in zip(wake_vr.strength, wake_cd.strength)
        svr .= scd
    end
    wake_vr.nwakes[] = nfree
    @test wake_vr isa pnl.PanelWake{pnl.VortexRing}

    # the quantity DirectWakePotential actually consumes is the wake scalar
    # potential q_f at body centroids; a constant-doublet panel and a vortex
    # ring share the same dipole potential, so q_f must agree to machine
    # precision
    qf_cd = zeros(N); pnl._wake_potential!(qf_cd, deepcopy(body0), (wake_cd,), DIRECT)
    qf_vr = zeros(N); pnl._wake_potential!(qf_vr, deepcopy(body0), (wake_vr,), DIRECT)
    rel_q = norm(qf_cd - qf_vr)/max(norm(qf_cd), eps())
    println("  ConstantDoublet vs VortexRing q_f rel diff: $(round(rel_q; sigdigits=3))")
    @test rel_q <= 1e-10

    # induced velocity agrees only up to the discrete-kernel difference between
    # the two representations (dipole-panel vs vortex-filament evaluation) — a
    # far looser bound than the potential, and not what the formulation uses
    bcd = deepcopy(body0); bcd.velocity .= 0.0
    pnl.influence!((bcd,), pnl.get_sources(wake_cd), DIRECT;
        precalc=true, scalar_potential=false, velocity=true)
    bvr = deepcopy(body0); bvr.velocity .= 0.0
    pnl.influence!((bvr,), pnl.get_sources(wake_vr), DIRECT;
        precalc=true, scalar_potential=false, velocity=true)
    rel_u = norm(bcd.velocity - bvr.velocity)/max(norm(bcd.velocity), eps())
    println("  ConstantDoublet vs VortexRing induced-velocity rel diff: $(round(rel_u; sigdigits=3))")
    @test rel_u <= 1e-2

    # resulting bound circulation γ = C·μE agrees through DirectWakePotential
    f = pnl.DirectWakePotential()
    stc = pnl.initialize_formulation(f, (deepcopy(body0),), (wake_cd,), solver0,
        DIRECT, DIRECT)
    bc = deepcopy(body0)
    run_aero!(bc, deepcopy(wake_cd), solver0; formulation=f, formulation_state=stc)
    stv = pnl.initialize_formulation(f, (deepcopy(body0),), (wake_vr,), solver0,
        DIRECT, DIRECT)
    bv = deepcopy(body0)
    run_aero!(bv, deepcopy(wake_vr), solver0; formulation=f, formulation_state=stv)
    gam_c = kutta_map(bc, view(bc.strength, :, 2))
    gam_v = kutta_map(bv, view(bv.strength, :, 2))
    rel_g = norm(gam_c - gam_v)/max(norm(gam_c), eps())
    println("  ConstantDoublet vs VortexRing γ rel diff: $(round(rel_g; sigdigits=3))")
    @test rel_g <= 1e-8
end

println("\nformulation_test.jl: all stages passed")
