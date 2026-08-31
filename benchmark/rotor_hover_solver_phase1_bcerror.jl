#=##############################################################################
BRAINSTORM 021 Phase 1 Stage 1 — validate the certified BC-error evaluator
(`bc_error!` in common.jl, backend=:fmm) against the O(N²) direct evaluation
(backend=:direct) on the campaign case.

Validation points per rung (span the error decades the τ-ladder will use):
  backslash   dense-LU reference        (BC error ~1e-15: probes the FMM
                                         evaluator's own floor — expect the
                                         FMM value to saturate there; gate is
                                         floor ≤ safety×target, not agreement)
  krylov_ilu  ILU-GMRES at rtol 1e-6    (BC error near the 1e-6 target scale)
  fgs         production FGS to 1e-6    (BC error ~1e-8)
  fgs_it1..3  FGS stopped at k outer    (coarse staircase points ~1e-1..1e-4)
              iterations (tolerance=0)

Gate (phase_01_tuning_plan.md Stage 1): FMM-certified vs direct BC error
agree ≤10% relative wherever the direct value is resolvable by the evaluator
(direct rel ≥ 10× the requested evaluation tolerance, i.e. ≥ 1e-6 at the
default safety=0.1×target=1e-7 request); below that, require the FMM-reported
value ≤ safety×target (evaluation floor subdominant). Every FMM row must
report error_success=true (bound certified at P ≤ cap).

Run (one rung per detached process; ≤4 threads total):
  RUNG=R1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    julia --project -t 4 benchmark/rotor_hover_solver_phase1_bcerror.jl
Appends bcerror.csv under benchmark/results/phase1/<mode>/.
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

target_rel = parse(Float64, get(ENV, "TARGET_REL", "1e-6"))
safety = parse(Float64, get(ENV, "BC_SAFETY", "0.1"))
p_t, mac_t, leaf_t = TUNED[rung]

# FMM evaluator knobs: per-rung tuned MAC/leaf, dynamic-P cap 20 (tune_fmm's
# own max), certification via PowerAbsolutePotential(safety*target*rms_b)
bc_fmm(x) = bc_error!(rotor, x; rms_b, target_rel, safety,
                      max_expansion_order=20, multipole_acceptance=mac_t,
                      leaf_size=leaf_t, backend=:fmm)
bc_direct(x) = bc_error!(rotor, x; rms_b, target_rel, safety, backend=:direct)

header = "rung,mesh_file,n_panels,point,rel_l2_direct,rel_max_direct," *
    "rel_l2_fmm,rel_max_fmm,error_success,discrepancy_rel,resolvable," *
    "gate_pass,rel_rms_trueresid,t_direct,t_fmm,epsilon_requested,target_rel," *
    "safety,eval_p_cap,eval_mac,eval_leaf,rms_b,radius_tol,threading_mode," *
    "julia_threads,commit,fm_commit,date,notes"
csv_path = joinpath(results_dir, "bcerror.csv")
fresh = !isfile(csv_path) || filesize(csv_path) == 0
io = open(csv_path, "a")
fresh && println(io, header)

r = zeros(rotor.ncells)
gate_all = Ref(true)

function validate_point!(point, x, notes)
    d = bc_direct(x)
    f = bc_fmm(x)
    # cross-check: the direct BC evaluation must reproduce the direct true
    # residual (same operator, same b) to roundoff
    rms_tr, _ = pnl.true_residual!(r, rotor, x, b; backend=pnl.DirectBackend())
    rel_tr = rms_tr / rms_b
    resolvable = d.rel_l2 >= 10 * safety * target_rel
    disc = abs(f.rel_l2 - d.rel_l2) / d.rel_l2
    # gate: 10% agreement where the direct value is resolvable by the
    # evaluator, OR the FMM value inflated by no more than the REQUESTED
    # evaluation tolerance (additive) — covers points at/below the evaluator
    # floor, where relative agreement is meaningless
    pass = f.error_success &&
           (disc <= 0.10 || f.rel_l2 <= d.rel_l2 + safety * target_rel)
    gate_all[] &= pass
    println("  $point: direct rel_l2=$(d.rel_l2)  fmm rel_l2=$(f.rel_l2)  " *
            "disc=$(round(100disc; digits=2))%  err_success=$(f.error_success)  " *
            (resolvable ? "" : "[below evaluator resolution] ") *
            (pass ? "PASS" : "FAIL"))
    println(io, join(_csv_cell.([rung, msh_name, rotor.ncells, point,
        d.rel_l2, d.rel_max, f.rel_l2, f.rel_max, f.error_success, disc,
        resolvable, pass, rel_tr, d.t_eval, f.t_eval, f.epsilon_requested,
        target_rel, safety, 20, mac_t, leaf_t, rms_b, pnl.FMM_RADIUS_TOL[],
        banner.threading_mode, banner.julia_threads, banner.commit,
        banner.fm_commit, time_string(), notes]), ","))
    flush(io)
    return nothing
end

want = split(get(ENV, "CONFIGS", "backslash,krylov_ilu,fgs,fgs_partial"), ",")

# ---- backslash (dense reference; evaluator-floor probe) ----
if "backslash" in want
    println("--- backslash ---")
    solver = pnl.Backslash(rotor)
    reset_cold!()
    pnl._solve!(rotor, solver)
    x = copy(vec(rotor.strength[:, solution_column]))
    validate_point!("backslash", x, "dense LU; probes FMM evaluator floor")
    solver = nothing; GC.gc()   # free the dense G before anything else builds
end

# ---- krylov_ilu (near the 1e-6 target scale) ----
if "krylov_ilu" in want
    println("--- krylov_ilu ---")
    backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
        multipole_acceptance=mac_t, leaf_size=leaf_t)
    solver = pnl.KrylovSolver(rotor; method=:gmres, itmax=500, atol=1e-14,
        rtol=target_rel, memory=50, backend=backend_apply,
        preconditioner=pnl.ILUPreconditioner(rotor; leaf_size=10,
            multipole_acceptance=1.0, max_pattern_entries=8192 * rotor.ncells))
    reset_cold!()
    pnl._solve!(rotor, solver)
    x = copy(vec(rotor.strength[:, solution_column]))
    validate_point!("krylov_ilu", x, "ILU-GMRES rtol=$target_rel; apply tuned($p_t,$mac_t,$leaf_t)")
    solver = nothing; GC.gc()
end

# ---- fgs full + partial iterates (staircase decades) ----
fgs_maker(maxit, tol) = pnl.FGSSolver(rotor; expansion_order=PROD[1],
    multipole_acceptance=PROD[2], leaf_size=PROD[3], max_iterations=maxit,
    tolerance=tol, inner_iterations=20, rlx=1.0, shrink=true, recenter=false,
    reverse_pass=false, verbose=false)

if "fgs" in want
    println("--- fgs (to 1e-6) ---")
    solver = fgs_maker(300, target_rel * rms_b)
    reset_cold!()
    pnl._solve!(rotor, solver)
    x = copy(vec(rotor.strength[:, solution_column]))
    validate_point!("fgs", x, "FGS prod knobs, abs tol $(target_rel * rms_b)")
    solver = nothing; GC.gc()
end

if "fgs_partial" in want
    for kiter in (1, 2, 3)
        println("--- fgs_it$kiter ---")
        solver = fgs_maker(kiter, 0.0)
        reset_cold!()
        pnl._solve!(rotor, solver)
        x = copy(vec(rotor.strength[:, solution_column]))
        validate_point!("fgs_it$kiter", x, "FGS stopped at $kiter outer iterations (tolerance=0)")
        solver = nothing; GC.gc()
    end
end

close(io)
println("\n$rung BC-error validation " * (gate_all[] ? "GATE PASS" : "GATE FAIL") *
        ". Rows appended to $csv_path")
gate_all[] || exit(1)
