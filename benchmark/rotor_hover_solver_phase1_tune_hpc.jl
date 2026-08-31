#=##############################################################################
BRAINSTORM 021 Phase 1 Stage 6 — dense-free tune_fmm stage for R4–R7
(phase_01_hpc_procedure.md step 2).

The local TUNE stage (rotor_hover_solver_phase1.jl TUNE=1) needs a Backslash
x_ref; dense G does not fit at R6–R7. Here the solved strength state for
tune_fmm comes from an FGS solve at conservative knobs, verified by the
certified BC evaluator (reference-free) before tuning. Output: a
"tuned"-variant row in tune.csv (same schema as the local stage), which
phase1_case.jl reads as TUNED[rung] for R4–R7.

Run (HPC): RUNG=R4 CACHE_B=1 EXPECT_JULIA_THREADS=<N> THREADING_MODE=multi \
    julia --project -t <N> benchmark/rotor_hover_solver_phase1_tune_hpc.jl
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

target_rel = 1e-6
safety = 0.1

# ---- solved strength state via FGS (conservative knobs; no dense G) ----
# p=10/MAC 0.4/leaf 150/inner 10: every rung R1–R3 met 1e-6 comfortably at
# p=10-class applies; certified check below is the guard, not trust.
println("--- FGS solve for tune_fmm strength state ---")
t_solve = @elapsed begin
    solver = pnl.FGSSolver(rotor; expansion_order=10, multipole_acceptance=0.4,
        leaf_size=150, inner_iterations=10, max_iterations=100,
        tolerance=0.1 * target_rel * rms_b, rlx=1.0, shrink=true,
        recenter=false, reverse_pass=false, verbose=false)
    reset_cold!()
    pnl._solve!(rotor, solver)
end
x = copy(vec(rotor.strength[:, solution_column]))
solver = nothing; GC.gc()

# certified BC check with DEFAULT evaluator knobs (tuned knobs don't exist yet;
# with error_tolerance active, evaluator knobs affect cost, not certification)
e = bc_error!(rotor, x; rms_b, target_rel, safety, backend=:fmm)
println("FGS state: $(round(t_solve; digits=1)) s, certified BC rel = $(e.rel_l2) " *
        "(certified=$(e.error_success))")
(e.error_success && e.rel_l2 <= target_rel) ||
    error("tune-state solve failed certification (rel=$(e.rel_l2), " *
          "success=$(e.error_success)) — fix before tuning")

# ---- tune_fmm (same call as the local TUNE stage) ----
error_tolerance = FastMultipole.PowerAbsolutePotential(safety * target_rel * rms_b)
rotor.strength[:, solution_column] .= x     # tuning depends on strengths
println("\n--- tune_fmm (PowerAbsolutePotential($(safety * target_rel * rms_b))) ---")
t_tune = @elapsed tuned, _cache = FastMultipole.tune_fmm((rotor,), (rotor,);
    error_tolerance, multipole_acceptances=range(0.3, 0.7, step=0.1),
    max_expansion_order=20,
    scalar_potential=true, gradient=false, hessian=false, verbose=true)
p_seed = tuned.expansion_order
mac_seed = tuned.multipole_acceptance
leaf_seed = tuned.leaf_size_source[1]
println("seed (accuracy-only): p=$p_seed MAC=$mac_seed leaf=$leaf_seed " *
        "($(round(t_tune; digits=1)) s)")

# ---- Stage 1b: descend on MEASURED cost (2026-08-24) ------------------------
# Stock tune_fmm optimizes a SINGLE evaluation against an error tolerance. It
# does not measure the cost of an iterative solve, so it accepts any knob set
# that is accurate enough — including a large leaf_size that is ~2x more
# expensive per apply, which is exactly what Phase 1 exists to measure.
# Measured at R1 (job 13407693, fixed p=17/mac=0.5, current FM): leaf 9 = 15.5 s
# vs leaf 60 = 32.6 s at identical iteration count — a 2.1x cost difference that
# stock tune_fmm is blind to. The old campaign's picks show the same blindness
# (leaf 9 at R1 but 59/58 at R2/R3), so this is not a regression in the tuner,
# it is the wrong objective.
# tune_fmm_perturb (FastMultipole src/autotune.jl:337) coordinate-descends over
# (p, mac, leaf) on BENCHMARKED wall-clock subject to the same error_tolerance.
# It requires an accuracy-valid seed and errors if the seed misses tolerance —
# hence stock tune_fmm above feeds it. Set PERTURB_TUNE=0 to fall back to the
# accuracy-only seed (recorded as variant "tuned_seed_only").
#
# Two defects were fixed on 2026-08-24 after the descent stalled at R1 (leaf 45
# when the solve is fastest near leaf 9):
#
#   (a) NOISE. reps defaulted to 2. Raised per rung below. Ryan's cap is 5
#       trials anywhere; 1 is not used because min-of-2 is what discards the
#       first-touch/page-fault cost of the first call at a new leaf.
#
#   (b) BIAS (the bigger one). tune_fmm_perturb used to time each candidate
#       through FastMultipole's `Cache` path, which REBUILDS both octrees and
#       the interaction lists inside every timed call. A steady iterative solve
#       does not pay that per apply: FLOWPanel builds one FmmPlan and reuses
#       its trees/lists across every Krylov iteration
#       (src/FLOWPanel_fmm.jl, plan_slot branch). Build cost RISES as leaf
#       shrinks, so charging a full build to every apply penalizes small leaves
#       — exactly the direction of the stall. `tree_amortization` now says how
#       many applies share one build; the candidate cost is
#       t_build/tree_amortization + t_apply, and tree_amortization=1 reproduces
#       the old objective (correct for an unsteady run whose geometry moves
#       every step).
#
# TREE_AMORTIZATION default Inf: charge NOTHING for the tree build. Ryan's
# ruling 2026-08-24 — the panels-on-panels operator has frozen geometry, so its
# tree is built once a priori and reused across every solve iteration AND every
# timestep. Its build is a one-off that must not influence the choice of apply
# knobs. (A finite n prices a build reused over exactly n applies; n=1 prices a
# rebuild per apply, which is what tuning a PARTICLE field will need, since the
# particle tree must be rebuilt every timestep.)
_rung_num = something(tryparse(Int, replace(rung, "R" => "")), 1)
tune_reps = parse(Int, get(ENV, "TUNE_REPS", _rung_num <= 4 ? "5" : "2"))
1 <= tune_reps <= 5 || error("TUNE_REPS must be in 1:5 (Ryan's cap); got $tune_reps")
tree_amortization = parse(Float64, get(ENV, "TUNE_TREE_AMORTIZATION", "Inf"))
# wall-clock guard so a pathological rung cannot eat the job. Free at R1
# (<1 s per evaluation), real at R7 (419k panels).
_default_max_seconds = Dict("R1" => 900, "R2" => 1800, "R3" => 3600,
    "R4" => 7200, "R5" => 14400, "R6" => 21600, "R7" => 43200)
tune_max_seconds = parse(Float64,
    get(ENV, "TUNE_MAX_SECONDS", string(get(_default_max_seconds, rung, 3600))))
# Early abandonment (Ryan 2026-08-24): stop a trial the moment its running min
# passes 1.3x the fastest complete, error-satisfying candidate so far. Such a
# point can no longer be accepted, so the remaining reps are pure waste — this
# is what makes reps=5 affordable. TUNE_ABANDON_FACTOR=Inf disables it.
tune_abandon_factor = parse(Float64, get(ENV, "TUNE_ABANDON_FACTOR", "1.3"))

use_perturb = get(ENV, "PERTURB_TUNE", "1") == "1"
t_perturb = 0.0
tune_timed_out = false
if use_perturb
    println("\n--- tune_fmm_perturb (cost objective, seeded from above) ---")
    println("    reps=$tune_reps tree_amortization=$tree_amortization " *
            "max_seconds=$tune_max_seconds abandon_factor=$tune_abandon_factor")
    t_perturb = @elapsed tuned_c, _hist, tinfo = FastMultipole.tune_fmm_perturb(
        (rotor,), (rotor,);
        expansion_order=p_seed, multipole_acceptance=mac_seed,
        leaf_size_source=leaf_seed, error_tolerance,
        max_expansion_order=20,
        reps=tune_reps, tree_amortization=tree_amortization,
        max_seconds=tune_max_seconds, abandon_factor=tune_abandon_factor,
        scalar_potential=true, gradient=false, hessian=false, verbose=true)
    p0 = tuned_c.expansion_order
    mac0 = tuned_c.multipole_acceptance
    leaf0 = tuned_c.leaf_size_source isa Integer ? tuned_c.leaf_size_source :
            tuned_c.leaf_size_source[1]
    tune_timed_out = tinfo.timed_out
    println("cost-tuned: p=$p0 MAC=$mac0 leaf=$leaf0 " *
            "($(round(t_perturb; digits=1)) s, $(tinfo.n_candidates) candidates, " *
            "$(tinfo.n_abandoned) abandoned early)  " *
            "[seed was $p_seed/$mac_seed/$leaf_seed]")
    tune_timed_out && @warn "tune_fmm_perturb TIMED OUT for $rung — the knobs " *
        "below are best-so-far, NOT a converged minimum. Raise TUNE_MAX_SECONDS " *
        "and re-run (purge tune.csv first) before trusting this rung's cost rows."
else
    p0, mac0, leaf0 = p_seed, mac_seed, leaf_seed
    tune_reps = 0; tree_amortization = 0.0
    println("PERTURB_TUNE=0 — keeping the accuracy-only seed")
end
t_tune += t_perturb
tune_variant = use_perturb ? "tuned" : "tuned_seed_only"

# measure the tuned point: certified BC of the SAME x, evaluated at the tuned
# knobs (records achieved evaluation error + cost; direct floor infeasible here)
et = bc_error!(rotor, x; rms_b, target_rel, safety,
               max_expansion_order=20, multipole_acceptance=mac0,
               leaf_size=leaf0, backend=:fmm)
meets = et.error_success && et.rel_l2 <= target_rel

tune_header = "rung,mesh_file,n_panels,variant,expansion_order," *
    "multipole_acceptance,leaf_size,t_tune,t_eval,rms_residual,rel_rms," *
    "dt_rel_vs_tuned,meets_target,radius_tol,threading_mode,julia_threads," *
    "tune_reps,tree_amortization,tune_timed_out," *
    "commit,fm_commit,date"
csv_path = joinpath(results_dir, "tune.csv")
fresh = !isfile(csv_path) || filesize(csv_path) == 0
# header guard (2026-08-24): this driver APPENDS, and phase1_case.jl reads the
# knobs back by column NAME, so a pre-fix tune.csv would silently mix rows
# tuned under two different objectives. Purge the rung's CSVs and re-run.
if !fresh
    existing = strip(readline(csv_path))
    existing == tune_header || error(
        "$csv_path has a stale header (pre-2026-08-24 cost-tuner schema).\n" *
        "  found:  $existing\n  expect: $tune_header\n" *
        "Purge this rung's CSVs (keep bcache_R*.bin) and re-run the tune stage.")
end
open(csv_path, "a") do io
    fresh && println(io, tune_header)
    println(io, join(_csv_cell.([rung, msh_name, rotor.ncells, tune_variant,
        p0, mac0, leaf0, t_tune, et.t_eval, et.rel_l2 * rms_b, et.rel_l2,
        0.0, meets, pnl.FMM_RADIUS_TOL[], banner.threading_mode,
        banner.julia_threads, tune_reps, tree_amortization, tune_timed_out,
        banner.commit, banner.fm_commit,
        time_string()]), ","))
end
println("\n$rung tune stage done (meets_target=$meets). Row appended to $csv_path")
meets || exit(1)
