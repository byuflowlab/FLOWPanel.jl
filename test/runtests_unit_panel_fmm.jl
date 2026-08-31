# runtests_unit_panel_fmm.jl — 052d Phase 1 parity testset.
#
# Panels -> particles through the host tree-FMM route (PANEL_INFLUENCE_FMM /
# set_panel_influence_fmm!) vs the dense rectangular seam, through the REAL
# call site (`_sa_body_influence!`, pass 3). Asserts:
#   - U and J relative RMS error vs the dense reference converges
#     monotonically (with slack) as p sweeps {4, 6, 8, 10};
#   - the p=10 error meets a production-shaped threshold;
#   - a p=1 canary error is large (proves the farfield expansion is actually
#     exercised — an all-nearfield configuration would be vacuously exact);
#   - body (non-particle) targets are bit-identical dense vs routed (the I3
#     split: only the particle leg changes);
#   - the U-only production shape (BODY_HESSIAN_TO_PARTICLES=false) leaves J
#     untouched.
#
# Run: julia --project=. -e 'cd(mktempdir()); include(".../runtests_unit_panel_fmm.jl")'
#
# KNOWN FLAKE (local Julia 1.12.5, 2026-08-27): ~50% of runs segfault at the
# first JIT of `_sa_body_influence!` during the DENSE reference stage (the
# 052d route is inert there; thread-count independent; no BoundsError under
# --check-bounds=yes, which also runs 3/3 clean). This matches the Julia 1.12
# JIT instability that pinned the cluster to julia 1.11.7
# (MATRIX_OPERATOR_REFACTOR/052-plan-2026-08-24.md:20). Rerun, or run with
# --check-bounds=yes; surviving runs are deterministic (identical errors).

using Test
import Random
import FLOWPanel as pnl
import FLOWPanel.FLOWVPM as vpm
using LinearAlgebra: norm

if !isdefined(@__MODULE__, :make_dirichlet_diamond_body)
    include(joinpath(@__DIR__, "test_helpers.jl"))
end

relrms(a, b) = norm(vec(a) .- vec(b)) / max(norm(vec(b)), eps())

function panel_fmm_fixture(; nspan=8, n_particles=1500)
    Random.seed!(52)
    body = make_dirichlet_diamond_body(; nspan=nspan)      # 8*nspan panels
    body.strength .= 0.1 .* randn(body.ncells, 2)
    body.needs_velocity_gradient[] = true
    wake = pnl.PanelParticleWake(body; nwakerows=2, max_particles=4*n_particles)
    pfield = wake.pfield
    # cloud spanning nearfield to true farfield so the MAC admits M2L cells
    for _ in 1:n_particles
        X = [0.5 + 6.0*rand() - 1.0, 1.5*rand() - 0.25, 1.5*(rand() - 0.5)]
        vpm.add_particle(pfield, X, 0.05 .* randn(3), 0.05 + 0.1*rand())
    end
    # pass-3 target tuple: body control points (dense, I3 split invariant)
    # + the particle field (the routed leg)
    targets = (body, pfield)
    return body, wake, pfield, targets
end

function run_pass3!(body, pfield, targets; hessian_to_particles=true)
    body.velocity .= 0
    body.velocity_gradient .= 0
    body.potential .= 0
    pfield.particles[vpm.U_INDEX, 1:pfield.np] .= 0
    pfield.particles[vpm.J_INDEX, 1:pfield.np] .= 0
    pnl._set_core_sizes!((body,), :core_size_targets)
    pnl._sa_body_influence!(targets, (body,),
        pnl.FastMultipoleBackend(; expansion_order=10, multipole_acceptance=0.4,
            leaf_size=20);
        body_hessian_to_particles=hessian_to_particles)
    return (bodyU = copy(body.velocity),
            partU = copy(pfield.particles[vpm.U_INDEX, 1:pfield.np]),
            partJ = copy(pfield.particles[vpm.J_INDEX, 1:pfield.np]))
end

@testset verbose=true "052d panels->particles host-FMM route" begin
    body, wake, pfield, targets = panel_fmm_fixture()
    old_env = Dict(k => get(ENV, k, nothing) for k in
        ("PANEL_INFLUENCE_FMM_P", "PANEL_INFLUENCE_FMM_THETA",
         "PANEL_INFLUENCE_FMM_LEAF"))
    try
        ENV["PANEL_INFLUENCE_FMM_THETA"] = "0.4"
        ENV["PANEL_INFLUENCE_FMM_LEAF"] = "10"   # small case: force farfield

        pnl.set_gpu_influence!(:host)
        pnl.set_panel_influence_fmm!(:off)
        println("052d stage: dense reference"); flush(stdout)
        ref = run_pass3!(body, pfield, targets)
        @test norm(ref.partU) > 0
        @test norm(ref.partJ) > 0

        pnl.set_panel_influence_fmm!(:on)

        # p=1 canary: the route must actually change the particle leg
        ENV["PANEL_INFLUENCE_FMM_P"] = "1"
        println("052d stage: canary p=1"); flush(stdout)
        canary = run_pass3!(body, pfield, targets)
        err_canary = relrms(canary.partU, ref.partU)
        @test err_canary > 1e-7

        ps = (4, 6, 8, 10)
        errU = Float64[]
        errJ = Float64[]
        for p in ps
            ENV["PANEL_INFLUENCE_FMM_P"] = string(p)
            println("052d stage: sweep p=", p); flush(stdout)
            got = run_pass3!(body, pfield, targets)
            push!(errU, relrms(got.partU, ref.partU))
            push!(errJ, relrms(got.partJ, ref.partJ))
            # I3 split invariant: body targets take the dense path bit-for-bit
            @test got.bodyU == ref.bodyU
        end
        println("052d parity errors: ps=", ps,
            " errU=", join(errU, ","), " errJ=", join(errJ, ","),
            " err_canary=", err_canary)

        # monotone convergence (with slack for nearfield-dominated floors)
        for i in 2:length(ps)
            @test errU[i] <= errU[i-1] * 1.05 + 1e-12
            @test errJ[i] <= errJ[i-1] * 1.05 + 1e-12
        end
        @test errU[1] < 1e-2            # p=4 already usable
        @test errU[end] < 1e-5          # p=10 production-shaped threshold
        @test errJ[end] < 1e-4          # J converges slower (one more derivative)
        @test errU[end] < err_canary    # expansion order actually matters

        # U-only production shape: J rows must stay untouched
        ENV["PANEL_INFLUENCE_FMM_P"] = "8"
        println("052d stage: u-only"); flush(stdout)
        got_uonly = run_pass3!(body, pfield, targets; hessian_to_particles=false)
        @test norm(got_uonly.partJ) == 0
        @test relrms(got_uonly.partU, ref.partU) < 1e-4
    finally
        pnl.set_panel_influence_fmm!(:off)
        pnl.set_gpu_influence!(:off)
        for (k, v) in old_env
            v === nothing ? delete!(ENV, k) : (ENV[k] = v)
        end
    end
end
