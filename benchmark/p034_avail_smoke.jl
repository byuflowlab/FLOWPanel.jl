#=##############################################################################
034 Phase 0 deliverable 5: availability smoke (021 W1–W6 analog) on the
pitching wing, provisional ladder rung 1 (n_span=7, n_airfoil=89, n_endcap=5 →
1920 cells), panel wake, formulation VelocityThroughSources (pinned).

All four campaign arms run a short unsteady march (~10 steps) to completion in
ONE process (cold = zero-initial-guess, NOT fresh-process — Ryan 2026-09-23;
each arm still gets a FRESH body/wake/monitor set so nothing leaks across arms):

  backslash          — pnl.Backslash (dense LU, the incumbent)
  krylov_gmres       — pnl.KrylovSolver, unpreconditioned GMRES
  krylov_ilu_nfcache — KrylovSolver + near-field ILU + persistent plan +
                       dense near-field cache (021 production iterative winner)
  fgs                — pnl.FGSSolver, dagteam+backoff defaults (021/033 champion
                       family; knobs here are SEEDS, recorded in notes)

Per step each arm gets one certified bc_error! pass (benchmark/common.jl)
against the fixed t=0 RHS scale rms_b_t0 (normalization constant only — the
true b moves every step; pass/fail thresholds are Phase 1's job). Availability
criterion, judged FROM THE CSVs: run completes, every bc pass certifies
(error_success), all values finite.

Run (local smoke; BENCH_BLAS_THREADS=8 per the macOS OpenBLAS note in
p034_phase0_lhs_rebuild.jl — this build's BLAS pin is neither observable nor
effective, so the strict single-mode assert is satisfied by recording 8):
    THREADING_MODE=single EXPECT_JULIA_THREADS=1 BENCH_BLAS_THREADS=8 \
        julia --project -t 1 benchmark/p034_avail_smoke.jl
Smoke only — nothing here is publishable.
=###############################################################################

import FLOWPanel as pnl
import LinearAlgebra
using Printf

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "..", "examples", "pitching_wing.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "BRAINSTORM",
    "034_pitching_wing_solver_benchmarks", "data", "phase0_avail_smoke")
mkpath(OUTDIR)

banner = assert_and_banner()
open(joinpath(OUTDIR, "banner.txt"), "w") do io
    println(io, banner.text)
end

const RUNG1 = (; n_span=7, n_airfoil=89, n_endcap=5)   # 1920 cells (provisional)
const N_CYCLES_SMOKE = 0.06                            # ~10 steps at c_per_dt=0.5
const BC_TARGET_REL = 1e-6

################################################################################
# per-step instrumentation: wrap the pinned formulation (TimedFormulation
# pattern from benchmark/rotor_hover_solver_unsteady.jl:539, leaned down)
################################################################################

mutable struct SmokeFormulation{F} <: pnl.AbstractSolveFormulation
    inner::F
    t_solve::Vector{Float64}
    niter::Vector{Int}
    bcerr_rel::Vector{Float64}
    bcerr_relmax::Vector{Float64}
    bcerr_eps::Vector{Float64}
    bcerr_cert::Vector{Bool}
    t_bcerr::Vector{Float64}
    rms_b::Float64              # fixed t=0 normalization scale
    vel_pre::Matrix{Float64}
end

pnl.initialize_formulation(f::SmokeFormulation, args...) =
    pnl.initialize_formulation(f.inner, args...)
pnl.formulation_prewake!(f::SmokeFormulation, state, systems_tuple) =
    pnl.formulation_prewake!(f.inner, state, systems_tuple)

_niter_of(s) = s isa pnl.KrylovSolver || s isa pnl.FGSSolver ? s.niter : -1

function pnl.solve_formulation!(f::SmokeFormulation, state, systems,
        systems_tuple, wakes_tuple, body_solvers; kwargs...)
    # bc_error!'s entry contract: body.velocity holds the apparent velocity at
    # the control points — snapshot BEFORE the solve, restore for the pass
    f.vel_pre .= systems_tuple[1].velocity
    t0 = time_ns()
    out = pnl.solve_formulation!(f.inner, state, systems, systems_tuple,
                                 wakes_tuple, body_solvers; kwargs...)
    push!(f.t_solve, (time_ns() - t0) / 1e9)
    push!(f.niter, _niter_of(body_solvers[1]))

    body = systems_tuple[1]
    x = Vector(view(body.strength, :, 2))     # Dirichlet solution column
    body.velocity .= f.vel_pre
    e = bc_error!(body, x; rms_b=f.rms_b, target_rel=BC_TARGET_REL,
                  backend=:fmm)
    push!(f.bcerr_rel, e.rel_l2)
    push!(f.bcerr_relmax, e.rel_max)
    push!(f.bcerr_eps, e.epsilon_requested)
    push!(f.bcerr_cert, e.error_success)
    push!(f.t_bcerr, e.t_eval)
    return out
end

################################################################################
# arms (021 constructions seeded as hypotheses, never carried as measurements)
################################################################################

krylov_kw = (; method=:gmres, itmax=500, atol=1e-14, rtol=1e-8, memory=50,
             warmstart=false)

# FGS absolute tolerance is derived per arm from the measured t=0 RHS scale
# (see rms_b_t0 below): tol_abs = BC_TARGET_REL * rms_b_t0.
function make_arm_solver(config::String, body, fgs_tol_abs::Float64)
    config == "backslash" && return pnl.Backslash(body),
        "dense LU; niter=-1 by convention"
    config == "krylov_gmres" && return pnl.KrylovSolver(body; krylov_kw...),
        "rtol=1e-8;atol=1e-14;memory=50 (021 seed)"
    config == "krylov_ilu_nfcache" && return pnl.KrylovSolver(body;
            krylov_kw...,
            preconditioner=pnl.ILUPreconditioner(body; leaf_size=10,
                multipole_acceptance=1.0,
                max_pattern_entries=8192 * body.ncells),
            cache_tree=true, cache_nearfield=true, persistent_plan=true,
            nearfield_cache_max_bytes=2 * 1024^3),
        "ilu leaf10/mac1.0;cache_nearfield;persistent_plan (021 production seed)"
    config == "fgs" && return pnl.FGSSolver(body;
            expansion_order=4, multipole_acceptance=0.5, leaf_size=50,
            inner_iterations=2, max_iterations=300, tolerance=fgs_tol_abs,
            rlx=1.0, shrink=true, recenter=false, reverse_pass=false,
            solution_history_length=0, project_solution=false),
        "p=4;mac=0.5;leaf=50;inner=2;tol_abs=$fgs_tol_abs;" *
        "dagteam+backoff defaults;f64 (knobs are SEEDS, not 021 measurements)"
    error("unknown arm $config")
end

################################################################################
# run all four arms
################################################################################

arms = ["backslash", "krylov_gmres", "krylov_ilu_nfcache", "fgs"]

steps_csv = joinpath(OUTDIR, "steps.csv")
summary_csv = joinpath(OUTDIR, "summary.csv")
steps_io = open(steps_csv, "w")
println(steps_io, "arm,step,t_solve,niter,bcerr_rel,bcerr_relmax," *
                  "bcerr_eps,bcerr_cert,t_bcerr")
summary_io = open(summary_csv, "w")
println(summary_io, "arm,n_panels,n_steps,completed,t_setup_s,t_sim_s," *
                    "mem_state_bytes,rms_b_t0,max_bcerr_rel,all_cert," *
                    "final_CL,notes")

for arm in arms
    println("\n===== arm: $arm =====")
    sim = prepare_pitching_wing(; RUNG1...,
        n_cycles=N_CYCLES_SMOKE,
        include_static_polar=false,
        save_vtk=false,
        # default Backslash is constructed inside prepare and replaced below;
        # at 1920 panels the extra ctor is negligible and keeps prepare stock
    )
    wing = sim.wing
    nsteps = length(sim.t_range)

    # t=0 RHS scale: RMS(φ_σ) with the freestream BC, via one direct body-only
    # pass (bc_error! with x=0 and rms_b=1 returns rel_l2 = RMS(φ_σ))
    pnl.apply_freestream!(wing, Vector(sim.Uinf(0.0)))
    rms_b_t0 = bc_error!(wing, zeros(wing.ncells); rms_b=1.0,
                         backend=:direct).rel_l2
    fgs_tol_abs = BC_TARGET_REL * rms_b_t0

    t_setup = @elapsed begin
        solver, notes = make_arm_solver(arm, wing, fgs_tol_abs)
    end
    mem = solver_state_bytes(solver)

    smoke = SmokeFormulation(pnl.VelocityThroughSources(),
        Float64[], Int[], Float64[], Float64[], Float64[], Bool[], Float64[],
        rms_b_t0, zeros(3, wing.ncells))

    completed = false
    t_sim = NaN
    err_msg = ""
    try
        t_sim = @elapsed pnl.simulate!((wing,), (sim.wake,), sim.frames,
            sim.maneuver!, sim.Uinf, sim.t_range;
            body_solvers=(solver,),
            backend=sim.backend,
            monitors=sim.monitors,
            path=nothing,
            name="p034_avail_$arm",
            set_Das_eta_freestream=NaN,
            formulation=smoke,
            verbose=false,
        )
        completed = true
    catch err
        err_msg = sprint(showerror, err)
        @warn "arm $arm FAILED" err_msg
    end

    for i in eachindex(smoke.t_solve)
        println(steps_io, join([arm, i - 1, smoke.t_solve[i], smoke.niter[i],
            smoke.bcerr_rel[i], smoke.bcerr_relmax[i], smoke.bcerr_eps[i],
            smoke.bcerr_cert[i], smoke.t_bcerr[i]], ","))
    end
    final_CL = sim.force_monitor.force[3, nsteps]
    max_bcerr = isempty(smoke.bcerr_rel) ? NaN : maximum(smoke.bcerr_rel)
    all_cert = !isempty(smoke.bcerr_cert) && all(smoke.bcerr_cert)
    note_field = replace(completed ? notes : notes * ";ERROR=" * err_msg,
                         "," => ";")
    println(summary_io, join([arm, wing.ncells, nsteps, completed, t_setup,
        t_sim, mem, rms_b_t0, max_bcerr, all_cert, final_CL,
        "\"$note_field\""], ","))
    flush(steps_io); flush(summary_io)
end

close(steps_io); close(summary_io)
println("\nWrote $steps_csv and $summary_csv")
