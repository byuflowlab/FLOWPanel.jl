#=##############################################################################
Unit tests for the implicit-Householder reduced Green solve
(`GreenHouseholderState`), the production dense `:area_mean` route adopted by
052e.2b Tier 0B-R (ADOPT ruling 2026-09-07; validated reference harness:
FastMultipole/MATRIX_OPERATOR_REFACTOR/scripts/tier0br_052e2b_householder_parity.jl).

Checks parity of trace, multiplier, gauge defect, and full-coordinate residual
against the retained bordered debug/reference route (`:area_mean_bordered`),
including on a deliberately incompatible RHS, plus state-reuse discipline and
route selection/validation. Tolerances follow the Tier 0B-R preregistration:
τ(N) = 1e3·√N·eps for parity, τ_g(N) = 1e2·√N·eps for the gauge defect.
=###############################################################################

using Test
import FLOWPanel as pnl
import LinearAlgebra as _GLA

if !isdefined(@__MODULE__, :make_dirichlet_diamond_body)
    include("test_helpers.jl")
end

@testset "Green Householder reduction (052e.2b)" begin
    body = make_dirichlet_diamond_body(nspan=12)
    N = body.ncells
    a = pnl._panel_areas(body)

    # deterministic smooth RHS (stands in for Sσ; parity needs no physics)
    CPs = body.controlpoints
    b = [sin(2.1*CPs[1, i]) + 0.3*cos(3.7*CPs[2, i]) + 0.2*CPs[3, i]
         for i in 1:N]

    tau = 1e3*sqrt(N)*eps()
    tau_g = 1e2*sqrt(N)*eps()

    # route selection: :area_mean → Householder, :area_mean_bordered → bordered
    gs_h = pnl._build_green_solve_state(body, :area_mean)
    gs_b = pnl._build_green_solve_state(body, :area_mean_bordered)
    @test gs_h isa pnl.GreenHouseholderState
    @test gs_b isa pnl.GreenSolveState

    q_h = copy(pnl._green_solve_q!(gs_h, b))
    lam_h = pnl._green_lambda(gs_h)
    q_b = copy(pnl._green_solve_q!(gs_b, b))
    lam_b = pnl._green_lambda(gs_b)

    # trace and multiplier parity vs the bordered reference
    scale = sqrt(sum(a .* q_b.^2)/sum(a))
    @test maximum(abs, q_h - q_b) <= tau*scale
    @test abs(lam_h - lam_b) <=
        tau*max(abs(lam_b), _GLA.norm(b)/_GLA.norm(a))

    # gauge defect: aᵀq = 0 up to roundoff
    @test abs(_GLA.dot(a, q_h)) <=
        tau_g*max(_GLA.norm(a)*_GLA.norm(q_h), eps())

    # full-coordinate residual (I−B)q + λa = b against an independent dense B
    B = zeros(N, N)
    pnl._assemble_B!(B, body)
    res = _GLA.norm(q_h - B*q_h .+ lam_h .* a - b) /
        max(_GLA.norm(b), eps())
    @test res <= 1e-10

    # incompatible RHS (constant shift): λ absorbs it; parity must hold
    bp = b .+ 1e-2*sum(a .* abs.(b))/sum(a)
    q_hp = copy(pnl._green_solve_q!(gs_h, bp))
    lam_hp = pnl._green_lambda(gs_h)
    q_bp = copy(pnl._green_solve_q!(gs_b, bp))
    lam_bp = pnl._green_lambda(gs_b)
    scale_p = sqrt(sum(a .* q_bp.^2)/sum(a))
    @test maximum(abs, q_hp - q_bp) <= tau*scale_p
    @test abs(lam_hp - lam_bp) <=
        tau*max(abs(lam_bp), _GLA.norm(bp)/_GLA.norm(a))

    # state reuse: re-solving the first RHS is bit-identical (buffer hygiene)
    @test pnl._green_solve_q!(gs_h, b) == q_h
    @test pnl._green_lambda(gs_h) == lam_h

    # constructor validation of the new gauge routes
    @test pnl.GreenReconstruction(gauge=:area_mean_bordered) isa
        pnl.GreenReconstruction
    @test pnl.HybridWakePotential(gauge=:area_mean_bordered) isa
        pnl.HybridWakePotential
    ks = pnl.KrylovSolver(body; method=:gmres, backend=pnl.DirectBackend())
    @test_throws ErrorException pnl.GreenReconstruction(
        gauge=:area_mean_bordered, green_solver=ks)
end
