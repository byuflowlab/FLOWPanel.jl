#=##############################################################################
BRAINSTORM 021 Phase 1 Stage 5 — the solver × rung table at tuned knobs,
with the certified BC error (the primary metric, Ryan 2026-08-15 ruling 1)
as a standard post-solve column.

Configs (CONFIGS env filters; default all):
  backslash_ldiv  dense-LU per-solve cost (ctor = assembly+LU = setup;
                  solve = rhs copy + triangular ldiv, update_G=false)
  krylov_gmres / krylov_jacobi / krylov_ilu
                  shared per-rung tuned apply knobs (ruling 3)
  fgs             the SHARED FGS knob set (Stage 3 winner τ from
                  fgsprecond.csv; ruling 3), run to BC 1e-6 with the
                  staircase-derived stopping tolerance; DUAL-ROW reporting
                  (ruling 4) when it overshoots: row_kind=target_1e-6 and
                  row_kind=last_above_1e-6 (error snapped to half decade).
                  If the τ=1e-6-tuned config differs from the shared set, an
                  extra fgs_tuned1e6 row is emitted for Ryan's comparison.
  fgmres_fgs      FGMRES + FGSPreconditioner at the shared set's fixed
                  sweeps/inner (Stage 3 winner row), tuned apply knobs

Requires fgstune_{staircase,selected}.csv and fgsprecond.csv for the rung
(Stages 2–3) unless CONFIGS excludes the FGS-family configs.

Run (chunk per rung/config-group; ≤4 threads):
  RUNG=R1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    julia --project -t 4 benchmark/rotor_hover_solver_phase1_table.jl
Appends solvetable.csv under benchmark/results/phase1/<mode>/.
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

using LinearAlgebra: lu!

target_rel = 1e-6
# knob-selection + adaptive-min-of-k helpers extracted 2026-08-18 (pure code
# motion) so Phase 2 scripts share them; requires the globals above.
include(joinpath(@__DIR__, "phase1_knobs.jl"))
safety = 0.1
p_t, mac_t, leaf_t = TUNED[rung]

bc_fmm(x) = bc_error!(rotor, x; rms_b, target_rel, safety,
                      max_expansion_order=20, multipole_acceptance=mac_t,
                      leaf_size=leaf_t, backend=:fmm)

want = split(get(ENV, "CONFIGS",
    "backslash_ldiv,krylov_gmres,krylov_jacobi,krylov_ilu,fgs,fgmres_fgs"), ",")

# ---- output ----
header = "rung,mesh_file,n_panels,config,row_kind,t_setup,t_solve_min,k_reps," *
    "niter,bc_rel_l2_certified,bc_certified,bc_rel_l2_snapped,rel_rms_direct," *
    "solver_knobs,apply_p,apply_mac,apply_leaf,rms_b,radius_tol," *
    "threading_mode,julia_threads,commit,fm_commit,date,notes"
csv_path = joinpath(results_dir, "solvetable.csv")
fresh = !isfile(csv_path) || filesize(csv_path) == 0
io = open(csv_path, "a")
fresh && println(io, header)

# ---- Resume support (2026-08-24) --------------------------------------------
# Benchmarks run on the preemptible `standby` QOS, where a preempted job is
# requeued and restarts from scratch. This driver APPENDS, so without a skip a
# restart would duplicate every row it had already written. Row identity is
# (rung, config, row_kind); `emit!` flushes each row the moment its measurement
# completes, so a row present in the CSV means that measurement is DONE and must
# not be repeated. Columns 1/4/5 never contain commas (only `notes`, the last
# column, is ever quoted), so the naive split is safe — same assumption
# validate_runs_csv makes.
landed = Set{Tuple{String,String}}()
if !fresh
    for line in Iterators.drop(eachline(csv_path), 1)
        isempty(strip(line)) && continue
        c = split(line, ",")
        length(c) >= 5 && String(c[1]) == rung &&
            push!(landed, (String(c[4]), String(c[5])))
    end
end
has_landed(config, row_kind) = (config, row_kind) in landed
if !isempty(landed)
    println("RESUME: $(length(landed)) row(s) already landed for $rung — " *
            "skipping " * join(sort(["$a/$b" for (a, b) in landed]), ", "))
end

# Ryan ruling 2026-08-18: the per-row direct O(N^2) residual is RETIRED
# (redundant with the certified BC column, validated to <=1% at R1-R3; at
# R5-R7/1T it would dominate job cost). DIRECT_RESID=1 re-enables it.
compute_direct_resid = get(ENV, "DIRECT_RESID", "0") == "1"
r = zeros(rotor.ncells)

function emit!(config, row_kind, t_setup, t_solve, niter, x, knobs, notes;
               snapped=nothing, k_used=k_reps)
    e = bc_fmm(x)
    rms_d = compute_direct_resid ?
        pnl.true_residual!(r, rotor, x, b; backend=pnl.DirectBackend())[1] : NaN
    println("  $config/$row_kind: setup=$(round(t_setup; digits=2))s " *
            "solve=$(round(t_solve; digits=4))s niter=$niter " *
            "bc=$(e.rel_l2) (certified=$(e.error_success))")
    println(io, join(_csv_cell.([rung, msh_name, rotor.ncells, config,
        row_kind, t_setup, t_solve, k_used, niter, e.rel_l2, e.error_success,
        snapped === nothing ? "" : snapped,
        compute_direct_resid ? rms_d / rms_b : "", knobs,
        p_t, mac_t, leaf_t, rms_b, pnl.FMM_RADIUS_TOL[],
        banner.threading_mode, banner.julia_threads, banner.commit,
        banner.fm_commit, time_string(), notes]), ","))
    flush(io)
end
common_note = "t excludes N^2 source assembly (restored); tree rebuild per " *
              "Krylov apply included"

# ---- backslash_ldiv ----
if "backslash_ldiv" in want && has_landed("backslash_ldiv", "standard")
    println("--- backslash_ldiv --- SKIPPED (already landed)")
elseif "backslash_ldiv" in want
    println("--- backslash_ldiv ---")
    t_setup = @elapsed solver = pnl.Backslash(rotor)
    t_solve, kk = adaptive_min_of_k(() -> pnl._solve!(rotor, solver);
                                    setup! = reset_cold!)
    x = copy(vec(rotor.strength[:, solution_column]))
    emit!("backslash_ldiv", "standard", t_setup, t_solve, -1, x,
          "dense LU (ctor=assembly+factorization)",
          "solve = rhs copy + triangular ldiv (update_G=false); " * common_note;
          k_used=kk)
    solver = nothing; GC.gc()
end

# ---- Krylov trio at shared tuned apply knobs ----
backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)
krylov_kw = (; itmax=500, atol=1e-14, rtol=target_rel, memory=50,
             backend=backend_apply)
krylov_configs = [
    ("krylov_gmres",  () -> pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...)),
    ("krylov_jacobi", () -> pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...,
        preconditioner=pnl.FastMultipole.JacobiPreconditioner((rotor,); cell_size=R/4))),
    ("krylov_ilu",    () -> pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...,
        preconditioner=pnl.ILUPreconditioner(rotor; leaf_size=10,
            multipole_acceptance=1.0, max_pattern_entries=8192 * rotor.ncells))),
]
for (name, make) in krylov_configs
    name in want || continue
    if has_landed(name, "standard")
        println("--- $name --- SKIPPED (already landed)")
        continue
    end
    println("--- $name ---")
    t_setup = @elapsed solver = make()
    t_solve, kk = adaptive_min_of_k(() -> pnl._solve!(rotor, solver);
                                    setup! = reset_cold!)
    x = copy(vec(rotor.strength[:, solution_column]))
    emit!(name, "standard", t_setup, t_solve, solver.niter, x,
          "rtol=1e-6;atol=1e-14;memory=50", common_note; k_used=kk)
    solver = nothing; GC.gc()
end

# ---- FGS family at the shared knob set (Stage 3 winner) ----
if "fgs" in want || "fgmres_fgs" in want
    winner = stage3_winner()
    winner === nothing && error("fgsprecond.csv has no winner for $rung — run Stages 2-3 first")
    println("shared FGS knob set (Stage 3 winner τ=$(winner.tau)): " *
            "p=$(winner.p) MAC=$(winner.mac) leaf=$(winner.leaf) " *
            "inner=$(winner.inner) sweeps=$(winner.sweeps)")
    sc = staircase_for(winner.p, winner.mac, winner.leaf, winner.inner)
    sc === nothing && error("No staircase rows for the winner knobs — rerun fgstune")
    i_cross = findfirst(t -> t[4] <= target_rel, sc)
    i_cross === nothing && error("Winner knobs never reach BC 1e-6 in the staircase")

    # Dual-row emission (ruling 4) is decided from the staircase alone, so the
    # expected fgs row set is known BEFORE any measurement — that is what makes
    # the resume check exact rather than a guess.
    bc_target_row = sc[i_cross][4]
    want_dual = bc_target_row <= target_rel / sqrt(10) && i_cross > 1
    fgs_expected = want_dual ? [("fgs", "target_1e-6"), ("fgs", "last_above_1e-6")] :
                               [("fgs", "target_1e-6")]
    if "fgs" in want && all(rk -> has_landed(rk...), fgs_expected)
        println("--- fgs (shared set) --- SKIPPED (all expected rows landed)")
    elseif "fgs" in want
        println("--- fgs (shared set) ---")
        tol_abs = margin_tol(sc, i_cross)
        make_fgs(maxit, tol) = pnl.FGSSolver(rotor; expansion_order=winner.p,
            multipole_acceptance=winner.mac, leaf_size=winner.leaf,
            inner_iterations=winner.inner, max_iterations=maxit, tolerance=tol,
            rlx=1.0, shrink=true, recenter=false, reverse_pass=false, verbose=false)
        t_setup = @elapsed solver = make_fgs(300, tol_abs)
        hist = pnl.ConvergenceHistory(:fgs_maxabs)
        reset_cold!(); pnl._solve_history!(rotor, solver, hist)
        niter = length(hist)
        t_solve, kk = adaptive_min_of_k(() -> pnl._solve!(rotor, solver);
                                        setup! = reset_cold!)
        x = copy(vec(rotor.strength[:, solution_column]))
        knobs = "p=$(winner.p);mac=$(winner.mac);leaf=$(winner.leaf);" *
                "inner=$(winner.inner);tol_abs=$tol_abs"
        emit!("fgs", "target_1e-6", t_setup, t_solve, niter, x, knobs,
              "shared set = Stage 3 winner tau=$(winner.tau); " * common_note;
              k_used=kk)

        # dual-row (ruling 4): last iterate above 1e-6, when the target row
        # overshoots the target by more than a half decade
        if want_dual && !has_landed("fgs", "last_above_1e-6")
            klast = sc[i_cross-1][1]
            solver2 = make_fgs(klast, 0.0)
            t2, kk2 = adaptive_min_of_k(() -> pnl._solve!(rotor, solver2);
                                        setup! = reset_cold!)
            x2 = copy(vec(rotor.strength[:, solution_column]))
            e2 = bc_fmm(x2)
            emit!("fgs", "last_above_1e-6", t_setup, t2, klast, x2, knobs,
                  "ruling-4 dual row; fixed $(klast) outer iterations";
                  snapped=snap_half_decade(e2.rel_l2), k_used=kk2)
            solver2 = nothing
        end
        solver = nothing; GC.gc()

        # comparison row when the τ=1e-6-tuned config differs from the shared set
        cfg6 = stage2_selected(1e-6)
        if cfg6 !== nothing && (cfg6.p, cfg6.mac, cfg6.leaf, cfg6.inner) !=
                               (winner.p, winner.mac, winner.leaf, winner.inner) &&
           !has_landed("fgs_tuned1e6", "comparison")
            println("--- fgs_tuned1e6 (differs from shared set; comparison row) ---")
            sc6 = staircase_for(cfg6.p, cfg6.mac, cfg6.leaf, cfg6.inner)
            i6 = sc6 === nothing ? nothing : findfirst(t -> t[4] <= target_rel, sc6)
            i6 === nothing && error("tuned-1e-6 staircase never crosses target")
            tol6 = margin_tol(sc6, i6)
            t_setup6 = @elapsed solver = pnl.FGSSolver(rotor;
                expansion_order=cfg6.p, multipole_acceptance=cfg6.mac,
                leaf_size=cfg6.leaf, inner_iterations=cfg6.inner,
                max_iterations=300, tolerance=tol6, rlx=1.0, shrink=true,
                recenter=false, reverse_pass=false, verbose=false)
            hist6 = pnl.ConvergenceHistory(:fgs_maxabs)
            reset_cold!(); pnl._solve_history!(rotor, solver, hist6)
            t6, kk6 = adaptive_min_of_k(() -> pnl._solve!(rotor, solver);
                                        setup! = reset_cold!)
            x6 = copy(vec(rotor.strength[:, solution_column]))
            emit!("fgs_tuned1e6", "comparison", t_setup6, t6, length(hist6), x6,
                  "p=$(cfg6.p);mac=$(cfg6.mac);leaf=$(cfg6.leaf);inner=$(cfg6.inner);tol_abs=$tol6",
                  "tau=1e-6-tuned config (differs from ruling-3 shared set); " * common_note;
                  k_used=kk6)
            solver = nothing; GC.gc()
        end
    end

    if "fgmres_fgs" in want && has_landed("fgmres_fgs", "standard")
        println("--- fgmres_fgs --- SKIPPED (already landed)")
    elseif "fgmres_fgs" in want
        println("--- fgmres_fgs (shared set, fixed sweeps) ---")
        t_setup = @elapsed begin
            P = pnl.FGSPreconditioner(rotor; sweeps=winner.sweeps,
                inner_iterations=winner.inner, rlx=1.0,
                expansion_order=winner.p, multipole_acceptance=winner.mac,
                leaf_size=winner.leaf, shrink=true, recenter=false)
            solver = pnl.KrylovSolver(rotor; method=:fgmres, krylov_kw...,
                preconditioner=P)
        end
        t_solve, kk = adaptive_min_of_k(() -> pnl._solve!(rotor, solver);
                                        setup! = reset_cold!)
        x = copy(vec(rotor.strength[:, solution_column]))
        emit!("fgmres_fgs", "standard", t_setup, t_solve, solver.niter, x,
              "sweeps=$(winner.sweeps);inner=$(winner.inner);p=$(winner.p);" *
              "mac=$(winner.mac);leaf=$(winner.leaf);rtol=1e-6",
              "shared set = Stage 3 winner tau=$(winner.tau); " * common_note;
              k_used=kk)
        P = nothing; solver = nothing; GC.gc()
    end
end

close(io)
println("\n$rung solvetable done. Rows appended to $csv_path")
