#=##############################################################################
Re-select the FGS tau winners from the immutable staircase CSV and verify a
post-crossing stopping margin for R1-R3.  The selected tolerance is the
geometric mean of the BC-crossing residual and the first later finite residual
that is strictly lower.  There is deliberately no fallback: a short staircase
is an error and must be extended.

The driver standardizes this campaign on Julia 4T / BLAS 1T and launches one
sequential process per rung.  Output is the new sidecar
benchmark/results/phase1/multi/fgstune_margin_verify.csv.
=###############################################################################

if get(ENV, "MARGIN_WORKER", "0") != "1"
    out = joinpath(@__DIR__, "results", "phase1", "multi",
                   "fgstune_margin_verify.csv")
    isfile(out) && error("Refusing to overwrite existing evidence sidecar $out")
    project = dirname(@__DIR__)
    for rung in ("R1", "R2", "R3")
        cmd = `$(Base.julia_cmd()) --project=$project --threads=4 $(@__FILE__)`
        env = copy(ENV)
        env["MARGIN_WORKER"] = "1"
        env["RUNG"] = rung
        env["EXPECT_JULIA_THREADS"] = "4"
        env["THREADING_MODE"] = "multi"
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["GOTO_NUM_THREADS"] = "1"
        env["OMP_NUM_THREADS"] = "1"
        println("\n=== FGS margin verification: $rung (Julia 4T / BLAS 1T) ===")
        run(setenv(cmd, env))
    end
    open(out, "w") do io
        wrote_header = false
        for rung in ("R1", "R2", "R3")
            rung_path = joinpath(dirname(out), "fgstune_margin_verify_$(rung).csv")
            for (i, line) in enumerate(eachline(rung_path))
                if i == 1
                    wrote_header || println(io, line)
                    wrote_header = true
                else
                    println(io, line)
                end
            end
        end
    end
    exit()
end

include(joinpath(@__DIR__, "common.jl"))
import LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
LinearAlgebra.BLAS.get_num_threads() == 1 || error("BLAS must be pinned to one thread")
include(joinpath(@__DIR__, "phase1_case.jl"))
# Case setup loads and exercises the direct backend, after which libblastrampoline
# may restore its default worker count.  Campaign policy is harness-owned: pin
# again after setup and assert before every production solve.
LinearAlgebra.BLAS.set_num_threads(1)
LinearAlgebra.BLAS.get_num_threads() == 1 || error("post-setup BLAS pin failed")
using CSV

const TAUS_MARGIN = (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
const STAIR_PATH = joinpath(@__DIR__, "results", "phase1", "multi",
                            "fgstune_staircase.csv")
const SIDE_PATH = joinpath(@__DIR__, "results", "phase1", "multi",
                           "fgstune_margin_verify_$(rung).csv")
isfile(SIDE_PATH) && error("Refusing to overwrite existing rung evidence $SIDE_PATH")

rows = [r for r in CSV.File(STAIR_PATH) if r.rung == rung]
isempty(rows) && error("no staircase rows for $rung")
groups = Dict{Tuple{Int,Float64,Int,Int},Vector{eltype(rows)}}()
for row in rows
    key = (row.p, row.mac, row.leaf, row.inner)
    push!(get!(groups, key, eltype(rows)[]), row)
end
foreach(v -> sort!(v; by=r -> r.iter), values(groups))

function winner_and_margin(tau)
    winner = nothing
    best_time = Inf
    i_cross = 0
    for (key, trajectory) in groups
        i = findfirst(r -> r.bc_certified && isfinite(r.bc_rel_l2) &&
                           r.bc_rel_l2 <= tau, trajectory)
        if i !== nothing && trajectory[i].t_wall < best_time
            winner = (key, trajectory)
            best_time = trajectory[i].t_wall
            i_cross = i
        end
    end
    winner === nothing && error("no $rung candidate reaches tau=$tau")
    key, trajectory = winner
    r_cross = trajectory[i_cross].mse_internal
    isfinite(r_cross) && r_cross > 0 || error(
        "$rung tau=$tau crossing residual is not positive and finite; extend staircase")
    i_successor = nothing
    for j in i_cross+1:length(trajectory)
        r = trajectory[j].mse_internal
        if isfinite(r) && 0 < r < r_cross
            i_successor = j
            break
        end
    end
    i_successor === nothing && error(
        "$rung tau=$tau has no lower finite post-crossing residual; extend staircase")
    r_successor = trajectory[i_successor].mse_internal
    return key, trajectory[i_cross], trajectory[i_successor],
           sqrt(r_cross * r_successor), best_time
end

target_rel = 1e-6
safety = 0.1
p_eval, mac_eval, leaf_eval = TUNED[rung]
bc_fmm_margin(x) = bc_error!(rotor, x; rms_b, target_rel, safety,
    max_expansion_order=20, multipole_acceptance=mac_eval,
    leaf_size=leaf_eval, backend=:fmm)

io = open(SIDE_PATH, "w")
println(io, "rung,tau,p,mac,leaf,inner,cross_iter,cross_residual," *
    "successor_iter,successor_residual,stopping_tol_abs,iterations,t_solve," *
    "bc_rel_l2,bc_certified,meets_tau,staircase_t_cross,julia_threads," *
    "blas_threads,commit,fm_commit,julia_version,date")

gate = true
for tau in TAUS_MARGIN
    # `bc_error!`'s FMM evaluator restores libblastrampoline's default after
    # returning on this Julia build.  Re-pin at the harness boundary before
    # every measured solve; FastMultipole itself remains free of global BLAS
    # mutations.
    LinearAlgebra.BLAS.set_num_threads(1)
    LinearAlgebra.BLAS.get_num_threads() == 1 || error("pre-solve BLAS pin failed")
    (p, mac, leaf, inner), crossrow, successor, tol_abs, stair_time =
        winner_and_margin(tau)
    solver = pnl.FGSSolver(rotor; expansion_order=p,
        multipole_acceptance=mac, leaf_size=leaf, inner_iterations=inner,
        max_iterations=max(80, successor.iter + 5), tolerance=tol_abs,
        rlx=1.0, shrink=true, recenter=false, reverse_pass=false, verbose=false)
    hist = pnl.ConvergenceHistory(:fgs_maxabs)
    # Constructor kernels can also initialize/reset the BLAS backend.
    LinearAlgebra.BLAS.set_num_threads(1)
    LinearAlgebra.BLAS.get_num_threads() == 1 || error("post-constructor BLAS pin failed")
    reset_cold!()
    t_solve = @elapsed pnl._solve_history!(rotor, solver, hist)
    LinearAlgebra.BLAS.get_num_threads() == 1 ||
        error("FGS solve changed the BLAS thread count")
    x = copy(vec(rotor.strength[:, solution_column]))
    LinearAlgebra.BLAS.set_num_threads(1)
    err = bc_fmm_margin(x)
    meets = err.error_success && err.rel_l2 <= tau
    global gate &= meets
    println("$rung tau=$tau tol=$tol_abs iterations=$(length(hist)) " *
            "BC=$(err.rel_l2) $(meets ? "PASS" : "FAIL")")
    println(io, join(_csv_cell.((rung, tau, p, mac, leaf, inner,
        crossrow.iter, crossrow.mse_internal, successor.iter,
        successor.mse_internal, tol_abs, length(hist), t_solve,
        err.rel_l2, err.error_success, meets, stair_time, Threads.nthreads(),
        banner.blas_threads, banner.commit, banner.fm_commit,
        banner.julia_version, time_string())), ","))
    flush(io)
    solver = nothing
    GC.gc()
end
close(io)
gate || exit(1)
