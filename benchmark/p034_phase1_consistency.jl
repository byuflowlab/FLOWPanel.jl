#=##############################################################################
034 Phase 1: consistency/calibration march (gate: all four arms agree on
certified BC <= 1e-6 and on CL hysteresis on 2 rungs; then the ladder and the
per-rung solver/FMM settings FREEZE).

One rung per invocation; each requested arm runs the SAME unsteady pitching
march (panel wake, VelocityThroughSources pinned) in one process with a fresh
body/wake/monitor set per arm (cold = zero-initial-guess, Ryan 2026-09-23).
Per step, each arm gets one certified bc_error! pass recording BOTH

  - bcerr_rel : relative L2 against the FIXED t=0 RHS scale rms_b_t0
    (normalization constant — the true b moves every step). The Phase 1 gate
    judges max-over-steps bcerr_rel <= 1e-6 per arm per rung, FROM THE CSV.
  - the 021 Phase 3 absolute order statistics (max/min/quartiles/rms) against
    the tolerance THE ARM PROMISED (bcerr_tol: FGS tolerance_abs; Krylov
    atol + rtol*||b||_step from the solver's own rhs; Backslash NaN) — the
    per-step self-consistency guard, reported alongside.

Identity signal: CL (and CM) hysteresis from the ForceMonitor — REPORTED,
never thresholded (021 ruling; cross-arm deltas are the Phase 1 deliverable).

Outputs under <outdir>/R<rung>/: banner.txt, steps_<arm>.csv (per-step),
history_<arm>.csv (time, t_over_T, alpha_deg, CL, CM, n_wake_rows),
summary.csv (one row per arm). Judge from the CSVs, never stdout.

Launch (local smoke — macOS BLAS caveat, ledger 2026-10-02):
    THREADING_MODE=single EXPECT_JULIA_THREADS=1 BENCH_BLAS_THREADS=8 \
        julia --project -t 1 benchmark/p034_phase1_consistency.jl
HPC (campaign worktree, strict pin):
    THREADING_MODE=single EXPECT_JULIA_THREADS=1 julia --project -t 1 ...
    THREADING_MODE=multi  EXPECT_JULIA_THREADS=<n> julia --project -t <n> ...

Env knobs:
  P034_RUNG      1-4 (default 1)
  P034_ARMS      comma list (default "backslash,krylov_gmres,krylov_ilu_nfcache,fgs")
  P034_NCYCLES   march length in pitch cycles (default 3.0; smoke uses ~0.1)
  P034_OUTDIR    output root (default BRAINSTORM .../data/phase1_consistency)
  P034_FGS_KNOBS "p,mac,leaf,inner,tolf" (default = Phase 1 R1 tune winner;
                 tol_abs = tolf * 1e-6 * rms_b_t0; RETUNE PER RUNG)
  P034_BCERR_EVERY  measure BC every k-th step (default 1; NaN rows otherwise)
=###############################################################################

import FLOWPanel as pnl
import LinearAlgebra
using Printf

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "..", "examples", "pitching_wing.jl"))

const LADDER = ((; n_span=7,  n_airfoil=89,  n_endcap=5),    # R1 1920
                (; n_span=13, n_airfoil=161, n_endcap=9),    # R2 6688
                (; n_span=19, n_airfoil=233, n_endcap=13),   # R3 14336
                (; n_span=27, n_airfoil=337, n_endcap=19))   # R4 30168
const RUNG_I = parse(Int, get(ENV, "P034_RUNG", "1"))
const RUNG = LADDER[RUNG_I]
const N_CYCLES = parse(Float64, get(ENV, "P034_NCYCLES", "3.0"))
const BC_TARGET_REL = 1e-6
const BCERR_EVERY = parse(Int, get(ENV, "P034_BCERR_EVERY", "1"))
const ARMS = split(get(ENV, "P034_ARMS",
    "backslash,krylov_gmres,krylov_ilu_nfcache,fgs"), ",")
#   Default = Phase 1 R1 tune winner (data/phase1_fgs_tune/summary.csv,
#   2026-10-02): p=6/mac=0.3/leaf=150/inner=5/tolf=1.0 — max bcerr_rel 9.3e-8
#   (10.7x margin), med t_solve tied with the p8 candidate, finer sweep
#   granularity (inner=5) for later warmstart work-matching.
const FGS_KNOBS = let v = split(get(ENV, "P034_FGS_KNOBS", "6,0.3,150,5,1.0"), ",")
    (; p=parse(Int, v[1]), mac=parse(Float64, v[2]), leaf=parse(Int, v[3]),
       inner=parse(Int, v[4]), tolf=parse(Float64, v[5]))
end

const OUTDIR = joinpath(get(ENV, "P034_OUTDIR",
    joinpath(@__DIR__, "..", "BRAINSTORM",
             "034_pitching_wing_solver_benchmarks", "data",
             "phase1_consistency")), "R$RUNG_I")
mkpath(OUTDIR)

banner = assert_and_banner()
open(joinpath(OUTDIR, "banner.txt"), "w") do io
    println(io, banner.text)
end

################################################################################
# arms (021 constructions as seeds; FGS knobs = Phase 1 tune winner via env)
################################################################################

krylov_kw = (; method=:gmres, itmax=500, atol=1e-14, rtol=1e-8, memory=50,
             warmstart=false)

function make_arm_solver(config, body, rms_b_t0)
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
    if config == "fgs"
        k = FGS_KNOBS
        tol_abs = k.tolf * BC_TARGET_REL * rms_b_t0
        return pnl.FGSSolver(body;
                expansion_order=k.p, multipole_acceptance=k.mac,
                leaf_size=k.leaf, inner_iterations=k.inner,
                max_iterations=300, tolerance=tol_abs,
                rlx=1.0, shrink=true, recenter=false, reverse_pass=false,
                solution_history_length=0, project_solution=false),
            "p=$(k.p);mac=$(k.mac);leaf=$(k.leaf);inner=$(k.inner);" *
            "tol_abs=$tol_abs(tolf=$(k.tolf));dagteam+backoff defaults;f64 " *
            "(Phase 1 R1-tuned; per-rung confirmation in this phase)"
    end
    error("unknown arm $config")
end

# What each arm PROMISED, in absolute units (021 Phase 3 _effective_tol):
# FGS = tolerance_abs; Krylov = atol + rtol*||b|| from the solver's OWN rhs
# (body.potential is overwritten by the solve); Backslash = NaN (no promise).
function effective_tol(solver)
    solver isa pnl.FGSSolver && return solver.tolerance
    solver isa pnl.KrylovSolver &&
        return krylov_kw.atol + krylov_kw.rtol * LinearAlgebra.norm(solver.rhs)
    return NaN
end

_niter_of(s) = s isa pnl.KrylovSolver || s isa pnl.FGSSolver ? s.niter : -1

################################################################################
# per-step instrumentation (021 Phase 3 TimedFormulation, leaned down)
################################################################################

mutable struct ConsistencyFormulation{F} <: pnl.AbstractSolveFormulation
    inner::F
    t_solve::Vector{Float64}
    niter::Vector{Int}
    bcerr_rel::Vector{Float64}      # rel L2 vs fixed rms_b_t0 (gate metric)
    bcerr_relmax::Vector{Float64}
    bcerr_max::Vector{Float64}      # absolute order stats (arm-promise guard)
    bcerr_min::Vector{Float64}
    bcerr_q1::Vector{Float64}
    bcerr_med::Vector{Float64}
    bcerr_q3::Vector{Float64}
    bcerr_tol::Vector{Float64}
    bcerr_eps::Vector{Float64}
    bcerr_cert::Vector{Bool}
    t_bcerr::Vector{Float64}
    n_wake_rows::Vector{Int}
    rms_b::Float64
    phi::Vector{Float64}
    vel_pre::Matrix{Float64}
    wake::pnl.PanelWake
    istep::Base.RefValue{Int}
end

pnl.initialize_formulation(f::ConsistencyFormulation, args...) =
    pnl.initialize_formulation(f.inner, args...)
pnl.formulation_prewake!(f::ConsistencyFormulation, state, systems_tuple) =
    pnl.formulation_prewake!(f.inner, state, systems_tuple)

function pnl.solve_formulation!(f::ConsistencyFormulation, state, systems,
        systems_tuple, wakes_tuple, body_solvers; kwargs...)
    # bc_error!'s entry contract: body.velocity holds the apparent velocity at
    # the control points — snapshot BEFORE the solve, restore for the pass
    f.vel_pre .= systems_tuple[1].velocity
    t0 = time_ns()
    out = pnl.solve_formulation!(f.inner, state, systems, systems_tuple,
                                 wakes_tuple, body_solvers; kwargs...)
    push!(f.t_solve, (time_ns() - t0) / 1e9)
    push!(f.niter, _niter_of(body_solvers[1]))
    push!(f.n_wake_rows, f.wake.nwakes[])

    f.istep[] += 1
    if (f.istep[] - 1) % BCERR_EVERY != 0
        for v in (f.bcerr_rel, f.bcerr_relmax, f.bcerr_max, f.bcerr_min,
                  f.bcerr_q1, f.bcerr_med, f.bcerr_q3, f.bcerr_tol,
                  f.bcerr_eps, f.t_bcerr)
            push!(v, NaN)
        end
        push!(f.bcerr_cert, false)
        return out
    end

    body = systems_tuple[1]
    x = Vector(view(body.strength, :, 2))
    tol = effective_tol(body_solvers[1])
    body.velocity .= f.vel_pre
    e = bc_error!(body, x; rms_b=f.rms_b, target_rel=BC_TARGET_REL,
                  backend=:fmm, phi_out=f.phi)
    ap = sort!(abs.(f.phi))
    np_ = length(ap)
    qi(q) = ap[clamp(round(Int, q * (np_ - 1)) + 1, 1, np_)]
    push!(f.bcerr_rel, e.rel_l2)
    push!(f.bcerr_relmax, e.rel_max)
    push!(f.bcerr_max, ap[end]);  push!(f.bcerr_min, ap[1])
    push!(f.bcerr_q1, qi(0.25));  push!(f.bcerr_med, qi(0.50))
    push!(f.bcerr_q3, qi(0.75))
    push!(f.bcerr_tol, tol)
    push!(f.bcerr_eps, e.epsilon_requested)
    push!(f.bcerr_cert, e.error_success)
    push!(f.t_bcerr, e.t_eval)
    return out
end

################################################################################
# run the requested arms
################################################################################

summary_csv = joinpath(OUTDIR, "summary.csv")
summary_io = open(summary_csv, "w")
println(summary_io, "arm,rung,n_panels,n_steps,n_cycles,completed,t_setup_s," *
    "t_sim_s,mem_state_bytes,rms_b_t0,max_bcerr_rel,all_cert,meets_gate," *
    "n_promise_viol,final_CL,threading_mode,julia_threads,blas_threads," *
    "commit,notes")

for arm in ARMS
    println("\n===== Phase 1 consistency: arm=$arm rung=R$RUNG_I " *
            "n_cycles=$N_CYCLES =====")
    sim = prepare_pitching_wing(; RUNG...,
        n_cycles=N_CYCLES, include_static_polar=false, save_vtk=false)
    wing = sim.wing
    nsteps = length(sim.t_range)

    # fixed t=0 RHS scale (normalization constant for bcerr_rel)
    pnl.apply_freestream!(wing, Vector(sim.Uinf(0.0)))
    rms_b_t0 = bc_error!(wing, zeros(wing.ncells); rms_b=1.0,
                         backend=:direct).rel_l2

    t_setup = @elapsed begin
        solver, notes = make_arm_solver(arm, wing, rms_b_t0)
    end
    mem = solver_state_bytes(solver)

    form = ConsistencyFormulation(pnl.VelocityThroughSources(),
        Float64[], Int[],
        Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
        Float64[], Float64[], Float64[], Bool[], Float64[], Int[],
        rms_b_t0, zeros(wing.ncells), zeros(3, wing.ncells),
        sim.wake, Ref(0))

    completed = false
    t_sim = NaN
    err_msg = ""
    try
        t_sim = @elapsed pnl.simulate!((wing,), (sim.wake,), sim.frames,
            sim.maneuver!, sim.Uinf, sim.t_range;
            body_solvers=(solver,), backend=sim.backend,
            monitors=sim.monitors, path=nothing,
            name="p034_ph1_$arm", set_Das_eta_freestream=NaN,
            formulation=form, verbose=false)
        completed = true
    catch err
        err_msg = sprint(showerror, err)
        @warn "arm $arm FAILED" err_msg
    end

    # per-step CSV (solver columns; monitor identity merged by index — the
    # maneuver callback fires once per t_range point including t=0)
    nrec = length(form.t_solve)
    steps_csv = joinpath(OUTDIR, "steps_$arm.csv")
    open(steps_csv, "w") do io
        println(io, "arm,step,time,t_over_T,alpha_deg,t_solve,niter," *
            "bcerr_rel,bcerr_relmax,bcerr_max,bcerr_min,bcerr_q1,bcerr_med," *
            "bcerr_q3,bcerr_tol,bcerr_eps,bcerr_cert,t_bcerr,n_wake_rows," *
            "CL,CM")
        for i in 1:nrec
            t = sim.t_range[i]
            alpha = sim.setup.alpha_mean_deg +
                    sim.setup.alpha_amp_deg * sin(sim.setup.omega * t)
            println(io, join([arm, i - 1, t, t / sim.setup.period,
                alpha, form.t_solve[i], form.niter[i],
                form.bcerr_rel[i], form.bcerr_relmax[i], form.bcerr_max[i],
                form.bcerr_min[i], form.bcerr_q1[i], form.bcerr_med[i],
                form.bcerr_q3[i], form.bcerr_tol[i], form.bcerr_eps[i],
                form.bcerr_cert[i], form.t_bcerr[i], form.n_wake_rows[i],
                sim.force_monitor.force[3, i],
                sim.force_monitor.moment[2, i]], ","))
        end
    end

    measured = filter(isfinite, form.bcerr_rel)
    max_bcerr = isempty(measured) ? NaN : maximum(measured)
    cert_meas = [form.bcerr_cert[i] for i in eachindex(form.bcerr_cert)
                 if isfinite(form.bcerr_rel[i])]
    all_cert = !isempty(cert_meas) && all(cert_meas)
    meets = completed && nrec == nsteps && all_cert && max_bcerr <= BC_TARGET_REL
    # arm-promise guard violations (band 1 + eps/tol, 021 ruling 2026-08-25)
    n_viol = count(i -> isfinite(form.bcerr_tol[i]) && isfinite(form.bcerr_max[i]) &&
                        form.bcerr_max[i] > form.bcerr_tol[i] *
                            (1 + form.bcerr_eps[i] / form.bcerr_tol[i]),
                   1:nrec)
    final_CL = nrec > 0 ? sim.force_monitor.force[3, nrec] : NaN
    note_field = replace(completed ? notes : notes * ";ERROR=" * err_msg,
                         "," => ";")
    println(summary_io, join([arm, "R$RUNG_I", wing.ncells, nsteps, N_CYCLES,
        completed, t_setup, t_sim, mem, rms_b_t0, max_bcerr, all_cert, meets,
        n_viol, final_CL, banner.threading_mode, banner.julia_threads,
        banner.blas_threads, banner.commit, "\"$note_field\""], ","))
    flush(summary_io)
    println("  wrote $steps_csv")
end

close(summary_io)
println("\nWrote $summary_csv")
