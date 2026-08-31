#=##############################################################################
BRAINSTORM 021 Phase 2b — FastGaussSeidel leaf-LU cache A/B.

Uses the frozen RHPC/Dirichlet case and Phase-1 τ=1e-6 settings. Each run
emits distinct cached/uncached rows; it never overwrites the historical
uncached baseline. Local smoke stays at R1 and ≤4 threads. Published runs use
K_REPS≥5 on the dedicated HPC node in separate 1T/BLAS-1 and recorded 64T
modes.

Run locally:
  RUNG=R1 K_REPS=5 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    julia --project -t 4 benchmark/rotor_hover_solver_phase2_leaf_lu.jl
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

using LinearAlgebra: ldiv!

const TARGET_REL = 1e-6
const SAFETY = 0.1

# Phase-1 verified τ=1e-6 rows (fgstune_margin_verify_R*.csv and the latest
# matching fgsprecond.csv rows). Extend only after Phase 1 freezes later rungs.
const FGS_CACHE_AB = Dict(
    "R1" => (p=6, mac=0.3, leaf=150, inner=5, tol=9.92983722e-8,
             precond_sweeps=9),
    "R2" => (p=8, mac=0.4, leaf=100, inner=10, tol=5.10785171e-8,
             precond_sweeps=6),
    "R3" => (p=6, mac=0.3, leaf=100, inner=5, tol=1.38387014e-7,
             precond_sweeps=13),
)
haskey(FGS_CACHE_AB, rung) || error(
    "No frozen Phase-1 FGS settings for $rung; finish Phase 1 before extending the cache A/B")
fgs_knobs = FGS_CACHE_AB[rung]
p_t, mac_t, leaf_t = TUNED[rung]

outdir = joinpath(@__DIR__, "results", "phase2", "leaf_lu_cache",
                  banner.threading_mode)
mkpath(outdir)
write(joinpath(outdir, "banner.txt"), banner.text * "\n")
csv_path = joinpath(outdir, "cache_ab.csv")
header = "rung,mesh_file,n_panels,config,cache_leaf_lu,p,mac,leaf,inner," *
    "precond_sweeps,setup_total_s,factor_build_s,factor_cache_bytes," *
    "leaf_pass_min_s,leaf_pass_alloc_bytes,solve_min_s,solve_alloc_bytes," *
    "solver_state_bytes,iterations,bc_rel_l2,bc_certified,k_reps," *
    "threading_mode,julia_threads,blas_threads,commit,fm_commit,date,notes"
fresh = !isfile(csv_path) || filesize(csv_path) == 0
io = open(csv_path, "a")
fresh && println(io, header)

bc_fmm(x) = bc_error!(rotor, x; rms_b, target_rel=TARGET_REL, safety=SAFETY,
    max_expansion_order=20, multipole_acceptance=mac_t,
    leaf_size=leaf_t, backend=:fmm)

"One forward solve of every dense leaf block, isolated from FMM/residual work."
function leaf_pass!(fgs)
    ms = fgs.self_matrices
    cache = fgs.leaf_lu_cache
    for i_leaf in eachindex(ms.sizes)
        _, rhs = FastMultipole.get_matrix_vector(ms, i_leaf)
        leaf_strengths = view(fgs.strengths, fgs.strengths_by_leaf[i_leaf])
        FastMultipole.solve_leaf!(leaf_strengths, ms, cache, i_leaf)
    end
    return nothing
end

function leaf_metrics(fgs)
    # Deterministic nonsingular RHS; a production solve overwrites it later.
    for i_leaf in eachindex(fgs.self_matrices.sizes)
        _, rhs = FastMultipole.get_matrix_vector(fgs.self_matrices, i_leaf)
        rhs .= sin.(eachindex(rhs))
    end
    t, _ = min_of_k(() -> leaf_pass!(fgs); k=k_reps, warmup=1)
    alloc = @allocated leaf_pass!(fgs)
    return t, alloc
end

function cache_metrics(fgs)
    cache = fgs.leaf_lu_cache
    return cache === nothing ? (0.0, 0) : (cache.build_time, cache.bytes)
end

function emit!(config, cache_leaf_lu, solver, fgs, setup_total, solve_min,
               solve_alloc, iterations, x, notes)
    factor_build, factor_bytes = cache_metrics(fgs)
    leaf_time, leaf_alloc = leaf_metrics(fgs)
    e = bc_fmm(x)
    println(io, join(_csv_cell.([rung, msh_name, rotor.ncells, config,
        cache_leaf_lu, fgs_knobs.p, fgs_knobs.mac, fgs_knobs.leaf,
        fgs_knobs.inner, config == "fgmres_fgs" ? fgs_knobs.precond_sweeps : 0,
        setup_total, factor_build, factor_bytes, leaf_time, leaf_alloc,
        solve_min, solve_alloc, solver_state_bytes(solver), iterations,
        e.rel_l2, e.error_success, k_reps, banner.threading_mode,
        banner.julia_threads, banner.blas_threads, banner.commit,
        banner.fm_commit, time_string(), notes]), ","))
    flush(io)
    println("  $config cache=$cache_leaf_lu: setup=$(round(setup_total; digits=3))s " *
            "factor=$(round(factor_build; digits=4))s solve=$(round(solve_min; digits=4))s " *
            "leaf=$(round(leaf_time; digits=6))s iter=$iterations BC=$(e.rel_l2)")
end

for cache_leaf_lu in (false, true)
    println("--- FGS cache_leaf_lu=$cache_leaf_lu ---")
    make_fgs = () -> pnl.FGSSolver(rotor;
        expansion_order=fgs_knobs.p, multipole_acceptance=fgs_knobs.mac,
        leaf_size=fgs_knobs.leaf, inner_iterations=fgs_knobs.inner,
        max_iterations=300, tolerance=fgs_knobs.tol, rlx=1.0,
        shrink=true, recenter=false, reverse_pass=false, verbose=false,
        cache_leaf_lu)
    # Compile the constructor path before measuring setup. The timed object is
    # still a fresh tree/matrix/factor build.
    warm_fgs = make_fgs()
    warm_fgs = nothing
    GC.gc()
    setup_total = @elapsed fgs_solver = make_fgs()
    hist = pnl.ConvergenceHistory(:fgs_maxabs)
    reset_cold!()
    pnl._solve_history!(rotor, fgs_solver, hist)
    solve_min, _ = min_of_k(() -> pnl._solve!(rotor, fgs_solver);
                            k=k_reps, warmup=1, setup! = reset_cold!)
    reset_cold!()
    solve_alloc = @allocated pnl._solve!(rotor, fgs_solver)
    x = copy(vec(rotor.strength[:, solution_column]))
    write_history_csv(joinpath(outdir,
        "history_fgs_cache$(cache_leaf_lu)_$(rung).csv"), hist)
    emit!("fgs", cache_leaf_lu, fgs_solver, fgs_solver.fgs, setup_total,
          solve_min, solve_alloc, length(hist), x,
          "Phase-1 frozen τ=1e-6 settings; constructor compiled before setup timing; cold solve; source assembly restored")
    fgs_solver = nothing
    GC.gc()

    println("--- FGMRES+FGS cache_leaf_lu=$cache_leaf_lu ---")
    backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
        multipole_acceptance=mac_t, leaf_size=leaf_t)
    make_fgmres_fgs = () -> begin
        P = pnl.FGSPreconditioner(rotor; sweeps=fgs_knobs.precond_sweeps,
            inner_iterations=fgs_knobs.inner, rlx=1.0,
            expansion_order=fgs_knobs.p, multipole_acceptance=fgs_knobs.mac,
            leaf_size=fgs_knobs.leaf, shrink=true, recenter=false,
            cache_leaf_lu)
        krylov = pnl.KrylovSolver(rotor; method=:fgmres, itmax=500,
            atol=1e-14, rtol=TARGET_REL, memory=50, backend=backend_apply,
            preconditioner=P, record_history=true)
        return P, krylov
    end
    warm_P, warm_krylov = make_fgmres_fgs()
    warm_P = nothing
    warm_krylov = nothing
    GC.gc()
    setup_total = @elapsed P, krylov = make_fgmres_fgs()
    solve_min, _ = min_of_k(() -> pnl._solve!(rotor, krylov);
                            k=k_reps, warmup=1, setup! = reset_cold!)
    reset_cold!()
    solve_alloc = @allocated pnl._solve!(rotor, krylov)
    x = copy(vec(rotor.strength[:, solution_column]))
    write_history_csv(joinpath(outdir,
        "history_fgmres_fgs_cache$(cache_leaf_lu)_$(rung).csv"), krylov.history)
    emit!("fgmres_fgs", cache_leaf_lu, krylov, P.solver.fgs, setup_total,
          solve_min, solve_alloc, krylov.niter, x,
          "Phase-1 frozen τ=1e-6 preconditioner/apply settings; constructor compiled before setup timing; cold solve; source assembly restored")
    P = nothing
    krylov = nothing
    GC.gc()
end

close(io)
println("\nCache A/B rows appended to $csv_path")
