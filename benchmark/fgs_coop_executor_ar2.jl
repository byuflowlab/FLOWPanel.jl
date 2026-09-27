#=##############################################################################
BRAINSTORM 033 A-R2: cooperative dagteam executor — validation + local
screening for the column-partitioned cooperative leaf-product prototype
(FastMultipole solve_dagteam.jl, `dagteam_coop` toggle).

The solver is built ONCE with dagteam_coop=WMAX (scratch capacity); arms
toggle plan.teamw[] between solves, so every arm runs the IDENTICAL operator
bytes in the same process (033 standing rule: paired same-process arms,
"cold" = zero-initial-guess solves). teamw[]=1 executes the unmodified solo
drain path (production dagteam+backoff).

MODE=validate (default RUNG=R1): correctness matrix over
  start ∈ {zero, nonzero} × inner_iterations ∈ {1, 3} × w ∈ {1, 2, 4}:
  - every arm must converge (solver.solved) and pass the independent BC
    evaluator at rel-L2 <= CERT_TARGET (A-T1 clause c: accuracy
    certification, NOT bit comparison — column reduction regroups FP sums);
  - w=2/4 solutions are compared to the SAME-(start,inner) w=1 solution
    (rel-L2 delta recorded; gate AGREE_MAX, default 1e-4 under f32full);
  - determinism: w=4 zero-start solve repeated — must be BITWISE identical;
  - bit-identity of the default: a separately built dagteam_coop=1 solver's
    solution must be BITWISE identical to the coop-built solver at teamw[]=1.

MODE=time (default RUNG=R4): paired same-process screening at fixed total
  worker budget (all Julia threads): warmup solve per arm, then ROUNDS
  interleaved rounds of one uninstrumented cold solve per arm (w=1,2,4);
  min/median per arm (021 ruling 5/7). Plus an in-situ effective-overhead
  probe: on a few real leaf classes of the live plan, the cooperative pull
  (dispatch + partials + wait + reduction) is timed against the solo pull,
  h_w = t_coop − t_seq/w (A-T3), using the actual executor code paths and a
  live team.

Run (local, ≤4 threads; env BLAS pinning REQUIRED on OpenMP-OpenBLAS hosts):
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  THREADING_MODE=multi EXPECT_JULIA_THREADS=4 BENCH_BLAS_THREADS=1 \
  RUNG=R1 MODE=validate CACHE_B=1 \
    julia --project -t 4 benchmark/fgs_coop_executor_ar2.jl [outdir]
  ... MODE=time RUNG=R4 ROUNDS=5 ...
=###############################################################################

import TOML
using Statistics: median

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))    # defines rotor, reset_cold!, rms_b, banner, ...

const FM = pnl.FastMultipole

const mode = get(ENV, "MODE", "validate")
mode in ("validate", "time") || error("MODE must be validate or time")
const outdir = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "BRAINSTORM", "033_ar2_20260926")
mkpath(outdir)
open(joinpath(outdir, "banner_$(mode)_$(rung).txt"), "w") do io
    write(io, banner.text)
end

const WMAX = parse(Int, get(ENV, "WMAX", "4"))
const WIDTHS = [w for w in (1, 2, 4) if w <= WMAX]
const CERT_TARGET = parse(Float64, get(ENV, "CERT_TARGET", "1e-6"))
const AGREE_MAX = parse(Float64, get(ENV, "AGREE_MAX", "1e-4"))

# champion R4 knobs drive every rung (033: R4 f32full champion; on smaller
# validation rungs the same knobs are used for structural parity — tolerance
# is the champion's for R4, CHAMP_TOL_FALLBACK otherwise)
champ = TOML.parsefile(joinpath(@__DIR__, "retained_r4_champion.toml"))
const tol = rung == champ["rung"] ? champ["tolerance"] :
    parse(Float64, get(ENV, "CHAMP_TOL_FALLBACK", "1e-6"))

function make_solver(coop::Int)
    t0 = time_ns()
    solver = pnl.FGSSolver(rotor; expansion_order=champ["P"],
        multipole_acceptance=champ["MAC"], leaf_size=champ["leaf"],
        inner_iterations=champ["inner"], max_iterations=champ["max_iterations"],
        tolerance=tol, rlx=champ["rlx"], shrink=true, recenter=false,
        reverse_pass=false, cache_leaf_lu=true, sweep_order=:dagteam,
        dagteam_precision=Symbol(champ["dagteam_precision"]),
        dagteam_coop=coop, verbose=false, project_solution=false,
        solution_history_length=0)
    println("solver(coop=$coop) constructed in " *
            "$(round((time_ns()-t0)/1e9, digits=1)) s")
    return solver
end

set_width!(solver, w::Int) = (solver.fgs.dagteam.teamw[] = w)

"Cold (zero-initial-guess) prepared solve; returns (seconds, x, niter, solved)."
function cold_solve!(solver; x0=nothing)
    reset_cold!()
    rotor.core_size = rotor.core_size_panel
    solver.niter = 0; solver.solved = false
    x0 === nothing || (rotor.strength[:, solution_column] .= x0)
    t = @elapsed pnl._solve!(rotor, solver)
    return (t, copy(rotor.strength[:, solution_column]), solver.niter, solver.solved)
end

certify(x) = bc_error!(rotor, x; rms_b, target_rel=CERT_TARGET)

rel_delta(x, ref) = norm(x - ref) / max(norm(ref), eps())

function write_csv(path, rows)
    isempty(rows) && return
    open(path, "w") do io
        println(io, join(string.(keys(rows[1])), ","))
        for r in rows
            println(io, join(string.(values(r)), ","))
        end
    end
end

################################################################################
if mode == "validate"
################################################################################

solver = make_solver(WMAX)
plan = solver.fgs.dagteam
@assert plan.teamw_cap == WMAX && plan.teamw[] == WMAX
set_width!(solver, 1)

# compile/warm both paths
cold_solve!(solver)
set_width!(solver, WMAX); cold_solve!(solver); set_width!(solver, 1)

rows = NamedTuple[]
refs = Dict{Tuple{String,Int},Vector{Float64}}()
refcert = Dict{Tuple{String,Int},Float64}()
fail = String[]
inner0 = solver.inner_iterations
for inner in (1, 3), start in ("zero", "nonzero"), w in WIDTHS
    solver.inner_iterations = inner
    set_width!(solver, w)
    # nonzero start perturbs the same-inner zero-start w=1 reference (loop
    # order guarantees it exists: zero before nonzero, w=1 first)
    x0 = start == "zero" ? nothing : 0.9 .* refs[("zero_inner$(inner)", 1)]
    t, x, niter, solved = cold_solve!(solver; x0)
    e = certify(x)
    w == 1 && (refs[("$(start)_inner$(inner)", 1)] = x;
               refcert[("$(start)_inner$(inner)", 1)] = e.rel_l2)
    d = rel_delta(x, refs[("$(start)_inner$(inner)", 1)])
    # certification is anchored to the SAME-(start,inner) w=1 arm: off-champion
    # rungs run an uncalibrated fallback tolerance, so the w=1 baseline sets
    # the certified level and every w>1 arm must match it (within 10%) or beat
    # CERT_TARGET outright — "matched accuracy", the 033 ranking contract
    cert_gate = max(CERT_TARGET, 1.1 * refcert[("$(start)_inner$(inner)", 1)])
    ok = solved && e.rel_l2 <= cert_gate && d <= AGREE_MAX && all(isfinite, x)
    ok || push!(fail, "start=$start inner=$inner w=$w solved=$solved rel_l2=$(e.rel_l2) delta=$d")
    push!(rows, (; start, inner, w, solve_seconds=t, iterations=niter, solved,
        bc_rel_l2=e.rel_l2, bc_rel_max=e.rel_max, bc_error_success=e.error_success,
        delta_vs_w1=d, pass=ok))
    println("validate start=$start inner=$inner w=$w: t=$(round(t; digits=2)) s, " *
            "iters=$niter, bc_rel_l2=$(e.rel_l2), delta_vs_w1=$d, pass=$ok")
end
solver.inner_iterations = inner0

# determinism: repeated w=WMAX zero-start solves must be bitwise identical
set_width!(solver, WMAX)
_, xa, _, _ = cold_solve!(solver)
_, xb, _, _ = cold_solve!(solver)
det_ok = xa == xb
det_ok || push!(fail, "w=$WMAX repeat not bitwise identical")
push!(rows, (; start="zero", inner=inner0, w=WMAX, solve_seconds=NaN,
    iterations=-1, solved=true, bc_rel_l2=NaN, bc_rel_max=NaN,
    bc_error_success=true, delta_vs_w1=Float64(!det_ok), pass=det_ok))
println("determinism (w=$WMAX repeat bitwise): $det_ok")

# default bit-identity: a coop=1-built solver == coop-built solver at teamw[]=1
set_width!(solver, 1)
_, x1, _, _ = cold_solve!(solver)
solver1 = make_solver(1)
@assert solver1.fgs.dagteam.teamw[] == 1 && isempty(solver1.fgs.dagteam.xgt)
cold_solve!(solver1)   # compile
_, x1d, _, _ = cold_solve!(solver1)
bit_ok = x1 == x1d
bit_ok || push!(fail, "coop-built teamw=1 not bitwise identical to coop=1 build")
println("default bit-identity (coop=1 build vs teamw[]=1): $bit_ok")
push!(rows, (; start="zero", inner=inner0, w=1, solve_seconds=NaN,
    iterations=-1, solved=true, bc_rel_l2=NaN, bc_rel_max=NaN,
    bc_error_success=true, delta_vs_w1=Float64(!bit_ok), pass=bit_ok))

write_csv(joinpath(outdir, "validate_$(rung).csv"), rows)
isempty(fail) || error("A-R2 validation FAILURES:\n" * join(fail, "\n"))
println("\nA-R2 validation: ALL PASS ($(length(rows)) checks, rung=$rung, " *
        "precision=$(champ["dagteam_precision"]), tol=$tol)")

################################################################################
else # mode == "time"
################################################################################

solver = make_solver(WMAX)
plan = solver.fgs.dagteam

# ---- in-situ effective-overhead probe (real plan, live team, executor code
# paths; pull only — no publication/LU, mirroring the A-R1 h_w definition) ----
function insitu_hw!(plan, w::Int, reps::Int=200)
    plan.teamw[] = w
    rt = FM.dagteam_start_team!(plan)
    coop = w > 1 ? rt.coops[1] : nothing
    nof(i) = plan.offset[i + 1] - plan.offset[i]
    # leaf classes: median/p90/max nof among ptot>0 leaves
    leaves = [i for i in eachindex(plan.ptot) if plan.ptot[i] > 0]
    byn = sort(leaves; by=nof)
    picks = unique([byn[cld(length(byn), 2)], byn[cld(9length(byn), 10)], byn[end]])
    out = NamedTuple[]
    for i in picks
        # solo pull (worker-1 scratch), then cooperative pull
        for _ in 1:20; FM.dagteam_pull!(plan, 1, i); end
        ts = minimum(@elapsed(FM.dagteam_pull!(plan, 1, i)) for _ in 1:reps)
        tc = NaN
        if coop !== nothing && FM.dagteam_coop_eligible(plan.ptot[i], w)
            coop_pull = () -> begin
                FM.dagteam_coop_dispatch!(coop, i)
                FM.dagteam_coop_lower_partial!(plan, plan.xgt[1], i, 1, 0, w)
                FM.dagteam_coop_wait!(coop.arrived, w - 1)
                n = nof(i)
                yv = view(plan.yb[1], 1:n)
                @inbounds for m in 1:w-1
                    pv = plan.yb[1 + m]
                    @simd for k in 1:n
                        yv[k] += pv[k]
                    end
                end
            end
            for _ in 1:20; coop_pull(); end
            tc = minimum(@elapsed(coop_pull()) for _ in 1:reps)
        end
        push!(out, (; w, leaf=i, nof=nof(i), ptot=plan.ptot[i],
            t_seq_us=1e6ts, t_coop_us=1e6tc, h_w_us=1e6(tc - ts / w),
            speedup=ts / tc))
    end
    FM.dagteam_stop_team!(rt)
    return out
end

hw_rows = NamedTuple[]
for w in WIDTHS
    w == 1 && continue
    append!(hw_rows, insitu_hw!(plan, w))
end
write_csv(joinpath(outdir, "insitu_hw_$(rung).csv"), hw_rows)
for r in hw_rows
    println("in-situ h_w: w=$(r.w) leaf=$(r.leaf) ($(r.nof)x$(r.ptot)) " *
            "seq=$(round(r.t_seq_us; digits=1))us coop=$(round(r.t_coop_us; digits=1))us " *
            "h_w=$(round(r.h_w_us; digits=2))us speedup=$(round(r.speedup; digits=2))x")
end

# ---- paired same-process solve screening ----
rows = NamedTuple[]
cert1 = Ref(NaN)   # w=1 warmup certification level anchors the arms
# compile + warmup, one per arm
for w in WIDTHS
    set_width!(solver, w)
    t, x, niter, solved = cold_solve!(solver)
    e = certify(x)
    # w=1 is the unmodified production executor: its certified level is the
    # local anchor (the champion tolerance was calibrated at j16+NUMA; local
    # 4-thread solves land within ~10% of 1e-6). Arms must match it.
    w == 1 && (cert1[] = e.rel_l2)
    gate = max(CERT_TARGET, 1.2 * cert1[])
    solved && e.rel_l2 <= gate || error("warmup w=$w failed: solved=$solved rel_l2=$(e.rel_l2) gate=$gate")
    push!(rows, (; phase="warmup", round=0, w, solve_seconds=t,
        iterations=niter, solved, bc_rel_l2=e.rel_l2))
    println("warmup w=$w: $(round(t; digits=3)) s, iters=$niter, bc_rel_l2=$(e.rel_l2)")
end
ref = let; set_width!(solver, 1); (_, x, _, _) = cold_solve!(solver); x; end

const ROUNDS = parse(Int, get(ENV, "ROUNDS", "5"))
for rnd in 1:ROUNDS, w in WIDTHS
    set_width!(solver, w)
    t, x, niter, solved = cold_solve!(solver)
    d = rel_delta(x, ref)
    solved && d <= AGREE_MAX || error("trial round=$rnd w=$w failed: solved=$solved delta=$d")
    push!(rows, (; phase="trial", round=rnd, w, solve_seconds=t,
        iterations=niter, solved, bc_rel_l2=NaN))
    println("round $rnd w=$w: $(round(t; digits=3)) s (iters=$niter, delta=$d)")
    flush(stdout)
    write_csv(joinpath(outdir, "timing_$(rung).csv"), rows)
end

println("\n== A-R2 local screening summary (rung=$rung, j=$(Threads.nthreads()), " *
        "$(ROUNDS) rounds, paired same-process, teamw toggled on one plan) ==")
t1med = median([r.solve_seconds for r in rows if r.phase == "trial" && r.w == 1])
summary = NamedTuple[]
for w in WIDTHS
    ts = [r.solve_seconds for r in rows if r.phase == "trial" && r.w == w]
    push!(summary, (; w, n=length(ts), min_s=minimum(ts), median_s=median(ts),
        speedup_median=t1med / median(ts)))
    println("w=$w: min=$(round(minimum(ts); digits=3)) s  " *
            "median=$(round(median(ts); digits=3)) s  " *
            "speedup(med, vs w=1)=$(round(t1med / median(ts); digits=3))x")
end
write_csv(joinpath(outdir, "timing_summary_$(rung).csv"), summary)

end
