#=##############################################################################
BRAINSTORM 021 Phase 1 — margin-verification replay (companion to the
2026-08-15 fgstune.jl edit, which switched the validation gate to a
post-crossing stopping margin: tolerance = geometric mean of the BC-crossing
internal residual and its first LOWER successor, so nondeterministic
threshold stopping can neither stop early above τ (the R3 zero-margin
τ=1e-5 FAIL) nor sit exactly on the crossing value).

This script re-runs ONLY the per-τ verification solves, reading the winners
from fgstune_selected.csv and the residual staircases from
fgstune_staircase.csv (both of record) — the 35-candidate tuning pool per
rung is NOT re-run. Writes the same fgstune_margin_verify.csv the edited
tuner would.

Run per rung:
  RUNG=R1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    julia --project -t 4 benchmark/rotor_hover_solver_phase1_margin_verify.jl
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

target_rel = 1e-6
safety = 0.1
maxit = parse(Int, get(ENV, "MAXIT", "60"))
p_t, mac_t, leaf_t = TUNED[rung]

const TAUS = [parse(Float64, t) for t in
              split(get(ENV, "TAUS", "1e-1,1e-2,1e-3,1e-4,1e-5,1e-6"), ",")]

bc_fmm(x) = bc_error!(rotor, x; rms_b, target_rel, safety,
                      max_expansion_order=20, multipole_acceptance=mac_t,
                      leaf_size=leaf_t, backend=:fmm)

function read_rows(path)
    isfile(path) || error("missing $path")
    lines = readlines(path)
    cols = Dict(String(c) => i for (i, c) in enumerate(split(lines[1], ",")))
    rows = [split(l, ",") for l in lines[2:end] if !isempty(strip(l))]
    return cols, rows
end

# winners (latest row per τ for this rung)
sel_cols, sel_rows = read_rows(joinpath(results_dir, "fgstune_selected.csv"))
winners = Dict{Float64,NamedTuple}()
for c in sel_rows
    c[sel_cols["rung"]] == rung || continue
    tau = parse(Float64, c[sel_cols["tau"]])
    winners[tau] = (; p=parse(Int, c[sel_cols["p"]]),
        mac=parse(Float64, c[sel_cols["mac"]]),
        leaf=parse(Int, c[sel_cols["leaf"]]),
        inner=parse(Int, c[sel_cols["inner"]]))
end

# staircases (latest run per knob set: later rows overwrite duplicate iters)
st_cols, st_rows = read_rows(joinpath(results_dir, "fgstune_staircase.csv"))
function staircase_for(p, mac, leaf, inner)
    pts = Dict{Int,NTuple{3,Float64}}()   # iter => (t_wall, mse, bc)
    for c in st_rows
        c[st_cols["rung"]] == rung || continue
        parse(Int, c[st_cols["p"]]) == p || continue
        parse(Float64, c[st_cols["mac"]]) == mac || continue
        parse(Int, c[st_cols["leaf"]]) == leaf || continue
        parse(Int, c[st_cols["inner"]]) == inner || continue
        pts[parse(Int, c[st_cols["iter"]])] = (
            parse(Float64, c[st_cols["t_wall"]]),
            something(tryparse(Float64, c[st_cols["mse_internal"]]), NaN),
            parse(Float64, c[st_cols["bc_rel_l2"]]))
    end
    its = sort(collect(keys(pts)))
    return its, [pts[i] for i in its]
end

"Post-crossing margin: (successor_iter_index, geomean tolerance), or nothing
when the staircase lacks a finite residual below the crossing (then the
staircase must be extended). Mirrors the edited tuner's stopping_margin."
function stopping_margin(mses, i_cross)
    r_cross = mses[i_cross]
    isfinite(r_cross) && r_cross > 0 || error("no positive finite residual at crossing")
    for j in i_cross+1:length(mses)
        if isfinite(mses[j]) && 0 < mses[j] < r_cross
            return j, sqrt(r_cross * mses[j])
        end
    end
    return nothing
end

# staircase-extension writer (same header/order as fgstune_staircase.csv)
io_stair = open(joinpath(results_dir, "fgstune_staircase.csv"), "a")

"""
Extend a winner's staircase past its BC crossing: one instrumented cold solve
at 10× tighter internal tolerance (0.01·target·rms_b) so post-crossing
residuals exist for the margin rule. Appends the new rows to
fgstune_staircase.csv (provenance: same schema, later rows supersede
duplicate iters) and returns (its, pts) like `staircase_for`.
"""
function extend_staircase!(w)
    println("    staircase too short for margin — extending (tol 0.01·target·rms_b)")
    t_setup = @elapsed solver = pnl.FGSSolver(rotor; expansion_order=w.p,
        multipole_acceptance=w.mac, leaf_size=w.leaf, inner_iterations=w.inner,
        max_iterations=maxit, tolerance=0.01 * target_rel * rms_b, rlx=1.0,
        shrink=true, recenter=false, reverse_pass=false, verbose=false)
    snaps = Vector{Float64}[]; ts = Float64[]; mses = Float64[]; its = Int[]
    t0 = Ref(UInt64(0))
    callback = (iteration, mse) -> begin
        t = time_ns()
        push!(ts, (t - t0[]) / 1e9)
        push!(mses, Float64(mse)); push!(its, iteration)
        FastMultipole.buffer_to_system_strength!((rotor,), solver.fgs.source_tree)
        push!(snaps, copy(vec(rotor.strength[:, solution_column])))
        nothing
    end
    reset_cold!()
    t0[] = time_ns()
    pnl._solve!(rotor, solver; callback)
    solver = nothing; GC.gc()
    pts = NTuple{3,Float64}[]
    for (k, x) in enumerate(snaps)
        e = bc_fmm(x)
        push!(pts, (ts[k], mses[k], e.rel_l2))
        println(io_stair, join(_csv_cell.([rung, msh_name, rotor.ncells,
            w.p, w.mac, w.leaf, w.inner, t_setup, its[k], ts[k], mses[k],
            e.rel_l2, e.error_success, target_rel, safety, mac_t, leaf_t,
            rms_b, pnl.FMM_RADIUS_TOL[], banner.threading_mode,
            banner.julia_threads, banner.commit, banner.fm_commit,
            time_string()]), ","))
    end
    flush(io_stair)
    return its, pts
end

ver_header = "rung,tau,p,mac,leaf,inner,cross_iter,cross_residual," *
    "successor_iter,successor_residual,stopping_tol_abs,iterations," *
    "t_solve,bc_rel_l2,bc_certified,meets_tau,commit,fm_commit,date"
io_ver = let path = joinpath(results_dir, "fgstune_margin_verify.csv")
    fresh = !isfile(path) || filesize(path) == 0
    io = open(path, "a"); fresh && println(io, ver_header); io
end

println("--- margin verification: fresh cold solves ---")
gate_all = true
for tau in TAUS
    haskey(winners, tau) || (println("  τ=$tau: no winner row, skipped"); continue)
    w = winners[tau]
    its, pts = staircase_for(w.p, w.mac, w.leaf, w.inner)
    isempty(its) && error("no staircase rows for τ=$tau winner knobs")
    bcs = [p[3] for p in pts]; mses = [p[2] for p in pts]
    i_cross = findfirst(<=(tau), bcs)
    i_cross === nothing && error("staircase for τ=$tau winner never crosses τ")
    margin = stopping_margin(mses, i_cross)
    if margin === nothing
        its, pts = extend_staircase!(w)
        bcs = [p[3] for p in pts]; mses = [p[2] for p in pts]
        i_cross = findfirst(<=(tau), bcs)
        i_cross === nothing && error("extended staircase never crosses τ=$tau")
        margin = stopping_margin(mses, i_cross)
        margin === nothing && error("extended staircase still lacks a post-crossing residual — raise MAXIT")
    end
    j, tol_abs = margin
    solver = pnl.FGSSolver(rotor; expansion_order=w.p,
        multipole_acceptance=w.mac, leaf_size=w.leaf,
        inner_iterations=w.inner, max_iterations=maxit,
        tolerance=tol_abs, rlx=1.0, shrink=true, recenter=false,
        reverse_pass=false, verbose=false)
    hist = pnl.ConvergenceHistory(:fgs_maxabs)
    reset_cold!()
    t_solve = @elapsed pnl._solve_history!(rotor, solver, hist)
    x = copy(vec(rotor.strength[:, solution_column]))
    e = bc_fmm(x)
    meets = e.error_success && e.rel_l2 <= tau
    global gate_all &= meets
    println("  τ=$tau: knobs p=$(w.p)/$(w.mac)/$(w.leaf)/$(w.inner) " *
            "tol=$tol_abs iters=$(length(hist)) t=$(round(t_solve; digits=2))s " *
            "bc=$(e.rel_l2) $(meets ? "PASS" : "FAIL")")
    println(io_ver, join(_csv_cell.([rung, tau, w.p, w.mac, w.leaf, w.inner,
        its[i_cross], mses[i_cross], its[j], mses[j], tol_abs, length(hist),
        t_solve, e.rel_l2, e.error_success, meets, banner.commit,
        banner.fm_commit, time_string()]), ","))
    solver = nothing; GC.gc()
end
close(io_ver); close(io_stair)

println("\n$rung margin verification " * (gate_all ? "GATE PASS" : "GATE FAIL") *
        ". Rows in fgstune_margin_verify.csv")
gate_all || exit(1)
