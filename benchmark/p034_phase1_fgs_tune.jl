#=##############################################################################
034 Phase 1: FGS knob retune on rung R1 (local smoke).

Phase 0 carried flag: the 021 coarse seed (p=4/mac=0.5/leaf=50/inner=2) carries
an apply-accuracy BC floor ~3-6e-6 rel that GROWS with wake rows — it cannot
certify the Phase 1 target BC <= 1e-6. Per the 021 ruling (2026-08-17 family),
the solver role must use the tau=1e-6-tuned config. This script marches a short
unsteady run (default ~16 steps, enough to see the wake-row growth trend) for
each candidate knob set and records per-step certified bc_error! so the winner
is chosen FROM THE CSV: smallest t_solve median subject to
max-over-steps bcerr_rel <= 1e-6 with every pass certified.

Run (local smoke; macOS BLAS caveat per ledger 2026-10-02):
    THREADING_MODE=single EXPECT_JULIA_THREADS=1 BENCH_BLAS_THREADS=8 \
        julia --project -t 1 benchmark/p034_phase1_fgs_tune.jl
Optional env: P034_RUNG (1-4, default 1), P034_NCYCLES (default 0.1),
P034_GRID (semicolon list "p,mac,leaf,inner", overrides the built-in grid).
Smoke only — nothing here is publishable; per-rung confirmation happens in the
Phase 1 HPC campaign.
=###############################################################################

import FLOWPanel as pnl
import LinearAlgebra
using Printf

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "..", "examples", "pitching_wing.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "BRAINSTORM",
    "034_pitching_wing_solver_benchmarks", "data", "phase1_fgs_tune")
mkpath(OUTDIR)

banner = assert_and_banner()
open(joinpath(OUTDIR, "banner.txt"), "w") do io
    println(io, banner.text)
end

const LADDER = ((; n_span=7,  n_airfoil=89,  n_endcap=5),    # R1 1920
                (; n_span=13, n_airfoil=161, n_endcap=9),    # R2 6688
                (; n_span=19, n_airfoil=233, n_endcap=13),   # R3 14336
                (; n_span=27, n_airfoil=337, n_endcap=19))   # R4 30168
const RUNG_I = parse(Int, get(ENV, "P034_RUNG", "1"))
const RUNG = LADDER[RUNG_I]
const N_CYCLES = parse(Float64, get(ENV, "P034_NCYCLES", "0.1"))
const BC_TARGET_REL = 1e-6

# Candidate grid: (expansion_order, multipole_acceptance, leaf_size,
# inner_iterations, tol_factor) with tol_abs = tol_factor * 1e-6 * rms_b_t0.
# First entry = Phase 0 seed (known apply-accuracy floor, kept as the control
# row). 021 tau=1e-6 rotor-rung winners (fgstune_verify.csv: R1 6/0.3/150/5,
# R2 8/0.4/100/10) and the FGSSolver constructor defaults (7/0.4/10) are SEEDS
# — hypotheses for this case, never measurements. 021 also tuned tol_abs below
# the raw 1e-6*rms_b (0.7-1.9e-7), hence the tol_factor dimension: once the
# apply error drops, the internal-residual tolerance can become the binding
# term in the certified BC.
default_grid = [
    (4, 0.5, 50, 2, 1.0),    # Phase 0 seed — CONTROL (expected to fail)
    (7, 0.4, 10, 2, 1.0),    # constructor defaults
    (6, 0.3, 150, 5, 1.0),   # 021 R1 tau=1e-6 winner (rotor) as seed
    (8, 0.4, 100, 10, 1.0),  # 021 R2 tau=1e-6 winner (rotor) as seed
    (7, 0.4, 10, 2, 0.3),
    (6, 0.3, 150, 5, 0.3),
    (7, 0.4, 10, 2, 0.1),
]
grid = if haskey(ENV, "P034_GRID")
    map(split(ENV["P034_GRID"], ";")) do s
        v = split(s, ",")
        (parse(Int, v[1]), parse(Float64, v[2]), parse(Int, v[3]),
         parse(Int, v[4]), parse(Float64, v[5]))
    end
else
    default_grid
end

################################################################################
# per-step instrumentation (SmokeFormulation pattern from p034_avail_smoke.jl)
################################################################################

mutable struct TuneFormulation{F} <: pnl.AbstractSolveFormulation
    inner::F
    t_solve::Vector{Float64}
    niter::Vector{Int}
    bcerr_rel::Vector{Float64}
    bcerr_relmax::Vector{Float64}
    bcerr_cert::Vector{Bool}
    t_bcerr::Vector{Float64}
    rms_b::Float64
    vel_pre::Matrix{Float64}
end

pnl.initialize_formulation(f::TuneFormulation, args...) =
    pnl.initialize_formulation(f.inner, args...)
pnl.formulation_prewake!(f::TuneFormulation, state, systems_tuple) =
    pnl.formulation_prewake!(f.inner, state, systems_tuple)

function pnl.solve_formulation!(f::TuneFormulation, state, systems,
        systems_tuple, wakes_tuple, body_solvers; kwargs...)
    f.vel_pre .= systems_tuple[1].velocity
    t0 = time_ns()
    out = pnl.solve_formulation!(f.inner, state, systems, systems_tuple,
                                 wakes_tuple, body_solvers; kwargs...)
    push!(f.t_solve, (time_ns() - t0) / 1e9)
    push!(f.niter, body_solvers[1].niter)

    body = systems_tuple[1]
    x = Vector(view(body.strength, :, 2))
    body.velocity .= f.vel_pre
    e = bc_error!(body, x; rms_b=f.rms_b, target_rel=BC_TARGET_REL,
                  backend=:fmm)
    push!(f.bcerr_rel, e.rel_l2)
    push!(f.bcerr_relmax, e.rel_max)
    push!(f.bcerr_cert, e.error_success)
    push!(f.t_bcerr, e.t_eval)
    return out
end

################################################################################
# sweep the grid
################################################################################

steps_csv = joinpath(OUTDIR, "steps.csv")
summary_csv = joinpath(OUTDIR, "summary.csv")
steps_io = open(steps_csv, "w")
println(steps_io, "config,step,t_solve,niter,bcerr_rel,bcerr_relmax," *
                  "bcerr_cert,t_bcerr")
summary_io = open(summary_csv, "w")
println(summary_io, "config,p,mac,leaf,inner,tolf,rung,n_panels,n_steps," *
                    "completed,t_setup_s,rms_b_t0,max_bcerr_rel,med_t_solve," *
                    "med_niter,all_cert,meets_target,final_CL,notes")

for (p, mac, leaf, inner, tolf) in grid
    cfg = "p$(p)_mac$(mac)_leaf$(leaf)_in$(inner)_tf$(tolf)"
    println("\n===== fgs config: $cfg (rung R$RUNG_I) =====")
    sim = prepare_pitching_wing(; RUNG...,
        n_cycles=N_CYCLES, include_static_polar=false, save_vtk=false)
    wing = sim.wing
    nsteps = length(sim.t_range)

    pnl.apply_freestream!(wing, Vector(sim.Uinf(0.0)))
    rms_b_t0 = bc_error!(wing, zeros(wing.ncells); rms_b=1.0,
                         backend=:direct).rel_l2
    tol_abs = tolf * BC_TARGET_REL * rms_b_t0

    t_setup = @elapsed solver = pnl.FGSSolver(wing;
        expansion_order=p, multipole_acceptance=mac, leaf_size=leaf,
        inner_iterations=inner, max_iterations=300, tolerance=tol_abs,
        rlx=1.0, shrink=true, recenter=false, reverse_pass=false,
        solution_history_length=0, project_solution=false)

    tune = TuneFormulation(pnl.VelocityThroughSources(),
        Float64[], Int[], Float64[], Float64[], Bool[], Float64[],
        rms_b_t0, zeros(3, wing.ncells))

    completed = false
    err_msg = ""
    try
        pnl.simulate!((wing,), (sim.wake,), sim.frames,
            sim.maneuver!, sim.Uinf, sim.t_range;
            body_solvers=(solver,), backend=sim.backend,
            monitors=sim.monitors, path=nothing,
            name="p034_fgstune_$cfg", set_Das_eta_freestream=NaN,
            formulation=tune, verbose=false)
        completed = true
    catch err
        err_msg = sprint(showerror, err)
        @warn "config $cfg FAILED" err_msg
    end

    for i in eachindex(tune.t_solve)
        println(steps_io, join([cfg, i - 1, tune.t_solve[i], tune.niter[i],
            tune.bcerr_rel[i], tune.bcerr_relmax[i], tune.bcerr_cert[i],
            tune.t_bcerr[i]], ","))
    end
    med(v) = isempty(v) ? NaN : sort(v)[cld(length(v), 2)]
    max_bcerr = isempty(tune.bcerr_rel) ? NaN : maximum(tune.bcerr_rel)
    all_cert = !isempty(tune.bcerr_cert) && all(tune.bcerr_cert)
    meets = completed && all_cert && max_bcerr <= BC_TARGET_REL
    final_CL = sim.force_monitor.force[3, nsteps]
    note = replace(completed ? "tol_abs=$tol_abs" :
                   "tol_abs=$tol_abs;ERROR=" * err_msg, "," => ";")
    println(summary_io, join([cfg, p, mac, leaf, inner, tolf, "R$RUNG_I",
        wing.ncells, nsteps, completed, t_setup, rms_b_t0, max_bcerr,
        med(tune.t_solve), med(tune.niter), all_cert, meets, final_CL,
        "\"$note\""], ","))
    flush(steps_io); flush(summary_io)
end

close(steps_io); close(summary_io)
println("\nWrote $steps_csv and $summary_csv")
