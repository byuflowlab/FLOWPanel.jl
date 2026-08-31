#=##############################################################################
BRAINSTORM 021 Phase 1 Stage 2 — purpose-built FGS τ-ladder tuner
(Ryan 2026-08-15 rulings 2/3; phase_01_tuning_plan.md Stage 2).

tune_fmm cannot tune FGS (no model of leaf-as-GS-block-size or
inner_iterations; measured picking anti-optimal leaf, fgs_check.csv). This
tuner measures instead of modeling:

- One instrumented cold solve per knob candidate {p, MAC, leaf, inner}
  (tolerance≈0 so it never stops early; maxit-capped): the FGS outer-iteration
  callback records wall time + internal residual, and a mid-solve
  `buffer_to_system_strength!` pulls each iterate's strengths (O(N) copy,
  timed AFTER the timestamp; with rlx=1.0 the buffers hold exactly the
  iterate whose residual the callback reports).
- Post-hoc (untimed), every iterate gets a certified BC error (`bc_error!`,
  FMM + PowerAbsolutePotential, cap P=20) — one solve yields the whole
  (cost, certified BC error) staircase.
- Cost-to-τ for τ ∈ {1e-1 … 1e-6} = wall time at the first iterate with
  BC ≤ τ. Objective (carried judgment, ruling 5): per-solve wall time,
  setup (ctor) timed separately.
- Knob search: coordinate descent seeded at the Phase 0 production knobs
  (10, 0.4, 150, inner 20), objective = cost to τ=1e-6 (the solver role);
  per-τ winners are then selected over ALL evaluated candidates' staircases,
  plus a few low-p seeds added to the pool so coarse-τ optima aren't hostage
  to the 1e-6 descent path.

Outputs (benchmark/results/phase1/<mode>/):
  fgstune_staircase.csv  one row per candidate × iterate
  fgstune_selected.csv   one row per τ: winning knobs, crossing iteration,
                         preconditioner sweep count (= iteration − 1),
                         cost-to-τ and internal residual at the BC crossing
  fgstune_margin_verify.csv  per-τ winner re-verified on a fresh cold solve
                         with a correctly directed post-crossing margin

Run (one rung per detached process; ≤4 threads total):
  RUNG=R1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    julia --project -t 4 benchmark/rotor_hover_solver_phase1_fgstune.jl
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

target_rel = 1e-6
safety = 0.1
maxit = parse(Int, get(ENV, "MAXIT", "40"))
descent_rounds = parse(Int, get(ENV, "DESCENT_ROUNDS", "2"))
p_t, mac_t, leaf_t = TUNED[rung]

const TAUS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]

# certified BC evaluator (shared knobs for every iterate; ε at the finest τ so
# one evaluation certifies every ladder level)
bc_fmm(x) = bc_error!(rotor, x; rms_b, target_rel, safety,
                      max_expansion_order=20, multipole_acceptance=mac_t,
                      leaf_size=leaf_t, backend=:fmm)

# ---- instrumented candidate evaluation ----

struct Staircase
    knobs::NTuple{4,Float64}        # p, mac, leaf, inner
    t_setup::Float64
    iters::Vector{Int}              # callback iteration index (1 = cold iterate)
    t_wall::Vector{Float64}         # seconds from solve start to callback
    mse::Vector{Float64}            # internal residual (max-abs) at iterate
    bc::Vector{Float64}             # certified BC rel L2 at iterate
    bc_success::Vector{Bool}
end

"First index with bc ≤ τ, or nothing."
crossing(sc::Staircase, tau) = findfirst(i -> sc.bc_success[i] && sc.bc[i] <= tau,
                                         eachindex(sc.bc))
cost_to(sc::Staircase, tau) = (i = crossing(sc, tau); i === nothing ? Inf : sc.t_wall[i])

"""
Return `(successor_index, tolerance)` for a safe post-crossing stopping margin.

The successor is the first finite internal residual *after* the certified BC
crossing that is strictly lower than the crossing residual.  Its geometric
mean with the crossing residual is below the crossing, so threshold stopping
must advance past the marginal crossing iterate.  A missing successor means
the measured staircase is too short; extending `MAXIT` is mandatory.
"""
function stopping_margin(sc::Staircase, i_cross::Int)
    r_cross = sc.mse[i_cross]
    isfinite(r_cross) && r_cross > 0 || error(
        "BC crossing at iterate $(sc.iters[i_cross]) has no positive finite " *
        "internal residual; extend the staircase (increase MAXIT)")
    i_successor = nothing
    for j in i_cross+1:length(sc.mse)
        if isfinite(sc.mse[j]) && 0 < sc.mse[j] < r_cross
            i_successor = j
            break
        end
    end
    i_successor === nothing && error(
        "No finite internal residual below crossing residual $r_cross after " *
        "iterate $(sc.iters[i_cross]); extend the staircase (increase MAXIT)")
    return i_successor, sqrt(r_cross * sc.mse[i_successor])
end

function run_candidate(p::Int, mac::Float64, leaf::Int, inner::Int)
    println("  candidate p=$p MAC=$mac leaf=$leaf inner=$inner ...")
    t_setup = @elapsed solver = pnl.FGSSolver(rotor; expansion_order=p,
        multipole_acceptance=mac, leaf_size=leaf, inner_iterations=inner,
        max_iterations=maxit, tolerance=0.1 * target_rel * rms_b, rlx=1.0,
        shrink=true, recenter=false, reverse_pass=false, verbose=false)

    snaps = Vector{Float64}[]
    ts = Float64[]; mses = Float64[]; its = Int[]
    t0 = Ref(UInt64(0))
    callback = (iteration, mse) -> begin
        t = time_ns()                       # timestamp BEFORE snapshot overhead
        push!(ts, (t - t0[]) / 1e9)
        push!(mses, Float64(mse)); push!(its, iteration)
        FastMultipole.buffer_to_system_strength!((rotor,), solver.fgs.source_tree)
        push!(snaps, copy(vec(rotor.strength[:, solution_column])))
        nothing
    end
    reset_cold!()
    t0[] = time_ns()
    pnl._solve!(rotor, solver; callback)
    t_end = (time_ns() - t0[]) / 1e9

    # if maxit exhausted, the last sweep phase never got a callback — append
    # the post-solve final iterate as an extra point
    xfinal = copy(vec(rotor.strength[:, solution_column]))
    if isempty(snaps) || norm(xfinal .- snaps[end]) > 0
        push!(ts, t_end); push!(mses, NaN); push!(its, its[end] + 1)
        push!(snaps, xfinal)
    end
    solver = nothing; GC.gc()

    bcs = Float64[]; succ = Bool[]
    for x in snaps
        e = bc_fmm(x)
        push!(bcs, e.rel_l2); push!(succ, e.error_success)
    end
    return Staircase((p, mac, leaf, inner), t_setup, its, ts, mses, bcs, succ)
end

# ---- CSV plumbing ----

stair_header = "rung,mesh_file,n_panels,p,mac,leaf,inner,t_setup,iter,t_wall," *
    "mse_internal,bc_rel_l2,bc_certified,target_rel,safety,eval_mac,eval_leaf," *
    "rms_b,radius_tol,threading_mode,julia_threads,commit,fm_commit,date"
io_stair = let path = joinpath(results_dir, "fgstune_staircase.csv")
    fresh = !isfile(path) || filesize(path) == 0
    io = open(path, "a"); fresh && println(io, stair_header); io
end
function write_staircase!(sc::Staircase)
    p, mac, leaf, inner = sc.knobs
    for i in eachindex(sc.iters)
        println(io_stair, join(_csv_cell.([rung, msh_name, rotor.ncells,
            Int(p), mac, Int(leaf), Int(inner), sc.t_setup, sc.iters[i],
            sc.t_wall[i], sc.mse[i], sc.bc[i], sc.bc_success[i], target_rel,
            safety, mac_t, leaf_t, rms_b, pnl.FMM_RADIUS_TOL[],
            banner.threading_mode, banner.julia_threads, banner.commit,
            banner.fm_commit, time_string()]), ","))
    end
    flush(io_stair)
end

# ---- candidate pool + coordinate descent on cost-to-1e-6 ----

cache = Dict{NTuple{4,Float64},Staircase}()
function evaluate!(p, mac, leaf, inner)
    key = (Float64(p), Float64(mac), Float64(leaf), Float64(inner))
    haskey(cache, key) && return cache[key]
    sc = run_candidate(Int(p), Float64(mac), Int(leaf), Int(inner))
    cache[key] = sc
    write_staircase!(sc)
    c6 = cost_to(sc, 1e-6)
    println("    -> $(length(sc.iters)) iterates, cost-to-1e-6 = " *
            (isfinite(c6) ? "$(round(c6; digits=2)) s" : "unreached") *
            ", setup $(round(sc.t_setup; digits=2)) s")
    return sc
end

P_SET = [6, 8, 10, 12, 14, 17]
MAC_SET = [0.3, 0.4, 0.5, 0.6]
LEAF_SET = [50, 100, 150, 200, 300]
INNER_SET = [5, 10, 20, 40]

current = (10.0, 0.4, 150.0, 20.0)     # Phase 0 production seed
println("--- coordinate descent (objective: cost to τ=1e-6) ---")
for round in 1:descent_rounds
    changed = false
    for (ci, values) in enumerate((P_SET, MAC_SET, LEAF_SET, INNER_SET))
        best_v = current[ci]; best_c = Inf
        for v in values
            cand = ntuple(j -> j == ci ? Float64(v) : current[j], 4)
            sc = evaluate!(cand...)
            c = cost_to(sc, 1e-6)
            c < best_c && (best_c = c; best_v = Float64(v))
        end
        best_v != current[ci] && (changed = true)
        global current = ntuple(j -> j == ci ? best_v : current[j], 4)
    end
    println("  round $round -> p=$(Int(current[1])) MAC=$(current[2]) " *
            "leaf=$(Int(current[3])) inner=$(Int(current[4]))" *
            (changed ? "" : " (converged)"))
    changed || break
end

# coarse-τ seed candidates DELETED (Ryan 2026-08-18): vestigial under the
# 2026-08-17 ruling — the preconditioner sweep count comes from the 1e-6
# config's OWN staircase, so coarse-τ-optimal other configs are unused.

# ---- per-τ selection over the full pool ----

sel_header = "rung,mesh_file,n_panels,tau,p,mac,leaf,inner,iter_cross," *
    "precond_sweeps,t_to_tau,t_setup,mse_at_cross,bc_at_cross,n_candidates," *
    "rms_b,threading_mode,julia_threads,commit,fm_commit,date,notes"
io_sel = let path = joinpath(results_dir, "fgstune_selected.csv")
    fresh = !isfile(path) || filesize(path) == 0
    io = open(path, "a"); fresh && println(io, sel_header); io
end

println("\n--- per-τ winners (over $(length(cache)) candidates) ---")
winners = Dict{Float64,Tuple{Staircase,Int}}()
for tau in TAUS
    best = nothing; best_c = Inf
    for sc in values(cache)
        c = cost_to(sc, tau)
        c < best_c && (best_c = c; best = sc)
    end
    if best === nothing
        println("  τ=$tau: UNREACHED by every candidate")
        continue
    end
    i = crossing(best, tau)
    winners[tau] = (best, i)
    p, mac, leaf, inner = best.knobs
    println("  τ=$tau: p=$(Int(p)) MAC=$mac leaf=$(Int(leaf)) inner=$(Int(inner)) " *
            "iter=$(best.iters[i]) t=$(round(best_c; digits=2))s bc=$(best.bc[i])")
    println(io_sel, join(_csv_cell.([rung, msh_name, rotor.ncells, tau,
        Int(p), mac, Int(leaf), Int(inner), best.iters[i],
        best.iters[i] - 1, best_c, best.t_setup, best.mse[i], best.bc[i],
        length(cache), rms_b, banner.threading_mode, banner.julia_threads,
        banner.commit, banner.fm_commit, time_string(),
        "mse_at_cross is diagnostic; production tolerance is recorded in " *
        "fgstune_margin_verify.csv; " *
        "precond_sweeps = fixed FGS sweep count reaching ~tau from zero seed"]), ","))
end
close(io_sel); close(io_stair)

# ---- validation gate: post-crossing margin, fresh cold solve ----

ver_header = "rung,tau,p,mac,leaf,inner,cross_iter,cross_residual," *
    "successor_iter,successor_residual,stopping_tol_abs,iterations," *
    "t_solve,bc_rel_l2,bc_certified,meets_tau,commit,fm_commit,date"
io_ver = let path = joinpath(results_dir, "fgstune_margin_verify.csv")
    fresh = !isfile(path) || filesize(path) == 0
    io = open(path, "a"); fresh && println(io, ver_header); io
end
println("\n--- verification: fresh cold solves at production stopping ---")
gate_all = true
for tau in TAUS
    haskey(winners, tau) || continue
    sc, i = winners[tau]
    p, mac, leaf, inner = sc.knobs
    j, tol_abs = stopping_margin(sc, i)
    solver = pnl.FGSSolver(rotor; expansion_order=Int(p),
        multipole_acceptance=mac, leaf_size=Int(leaf),
        inner_iterations=Int(inner), max_iterations=maxit,
        tolerance=tol_abs, rlx=1.0, shrink=true, recenter=false,
        reverse_pass=false, verbose=false)
    hist = pnl.ConvergenceHistory(:fgs_maxabs)
    reset_cold!()
    t_solve = @elapsed pnl._solve_history!(rotor, solver, hist)
    x = copy(vec(rotor.strength[:, solution_column]))
    e = bc_fmm(x)
    meets = e.error_success && e.rel_l2 <= tau
    global gate_all &= meets
    println("  τ=$tau: iters=$(length(hist)) t=$(round(t_solve; digits=2))s " *
            "bc=$(e.rel_l2) $(meets ? "PASS" : "FAIL")")
    println(io_ver, join(_csv_cell.([rung, tau, Int(p), mac, Int(leaf),
        Int(inner), sc.iters[i], sc.mse[i], sc.iters[j], sc.mse[j],
        tol_abs, length(hist), t_solve, e.rel_l2,
        e.error_success, meets, banner.commit, banner.fm_commit,
        time_string()]), ","))
    solver = nothing; GC.gc()
end
close(io_ver)

println("\n$rung FGS τ-ladder tuning done — verification " *
        (gate_all ? "GATE PASS" : "GATE FAIL") * ". CSVs under $results_dir")
gate_all || exit(1)
