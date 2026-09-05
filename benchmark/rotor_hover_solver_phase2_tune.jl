#=##############################################################################
BRAINSTORM 021 Phase 2 — apply-knob tuning against a REAL SOLVE, swept over a
near-field-cache MEMORY BUDGET (W1+W2, Ryan 2026-08-24).

Replaces two things:
  * the TUNE_CACHED=1 stage of rotor_hover_solver_phase2b_nearfield_cache.jl,
    which tuned on stock `tune_fmm` — the accuracy-only cost-model objective,
    i.e. the very defect fixed for the uncached path on 2026-08-24;
  * Phase 1's ownership of the shared apply knobs. Phase 1 is now a
    verification pass (solver-to-solver agreement); tuning lives here, because
    Phase 2 is what publishes benchmarks.

## Why a real solve

`tune_fmm_perturb` minimizes ONE `fmm!` apply, but the campaign publishes SOLVE
wall clock, and the two were MEASURED to disagree below leaf ~40 (R1, p17,
MAC0.5: the real solve falls monotonically to leaf 9 at 15.5 s, while the apply
proxy read flat between leaf 29 and 43). This driver hands
`tune_fmm_perturb` a `cost` closure that builds the actual `KrylovSolver` at
the candidate knobs and times `_solve!`, so the objective IS the published
quantity. It also removes every modeling decision the proxy needed:
`tree_amortization` is moot, because whatever the benchmark rebuilds per apply,
the closure rebuilds per apply too.

## Why a memory budget, why TOTAL memory, and why it is an AXIS

An uncapped cached tune walks to a DENSE operator. Measured at R1 with
cache_capped=false (i.e. the cost optimum genuinely wanted it): the cached
optimum sat at leaf 275-342, 687-862 MB, against a dense 8N^2 of 514 MB. If
assembly is free — and Ryan's ruling is that it IS free for tuning, the
panels-on-self operator being built once and reused across every iteration and
every timestep — then matrix-ful always beats matrix-free wherever memory
allows. So a memory cap does not guard the answer, it DETERMINES it, and a
hand-picked constant would be doing the science.

Ryan's ruling: sweep it, and sweep the TOTAL memory a node must hold, not just
the cache — "memory rungs should correspond to the total memory available, not
just that used by the cache." A budget therefore covers

    body + solver state (workspace, preconditioner, scratch) + FMM plan + cache

which puts every solver on ONE axis, `backslash_ldiv`'s dense 8N^2 included,
and turns "which solver is fastest" into "which is fastest on a node with X
GB" — the decision a user actually faces. It also matters mechanically: the
plan grows as the leaf SHRINKS, partly offsetting the cache's growth, so
pricing the cache alone would misrank small-leaf candidates.

Whether a candidate is CACHED is then an outcome, not an input: the cache is
used iff the configuration still fits with it. A rung or budget that admits no
cache is tuned uncached automatically — which R6 and R7 need, their cache
floors being 24 and 78 GiB (`ledger.md`, W0 map).

Use `benchmark/nfcache_feasibility.jl` (W0) to see what a rung needs before
choosing a ladder.
## The build cost

Excluded from the objective and MEASURED anyway, per Ryan's ruling. The closure
builds the plan+cache in an UNTIMED warm-up solve and times only warm solves
afterwards, so exclusion is by construction rather than by subtraction. The
build time and bytes are recorded per winning point.

Run:
  RUNG=R1 CACHE_B=1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    MEM_BUDGETS=0:1:4:16 julia --project -t 4 \
    benchmark/rotor_hover_solver_phase2_tune.jl

Env:
  RUNG                required
  MEM_BUDGETS         colon/comma list of TOTAL memory budgets in GiB
                      (body + solver state + FMM plan + cache), default
                      "4:8:16:32:64:128:256:500" (top = the pinned --mem=500G
                      on 512 GB zen3 nodes). NOT cache size — a budget is
                      what a node must hold. Whether a candidate is cached is
                      decided per candidate by what fits.
  TUNE_REPS           min-of-reps per candidate (default 5 at R1-R4, 2 above;
                      asserted <= 5, Ryan's cap)
  TUNE_ABANDON_FACTOR default 1.3
  TUNE_MAX_SECONDS    per-budget wall-clock backstop (default: per-rung table)
  NFCACHE_MAX_BUILD_TIME  hard cap on a single cache build (default Inf).
                      The build is SERIAL, so this can bind before memory does.
  FLOOR_LEAF          leaf used to price the cheapest configuration a rung
                      admits, for the below-floor check (default 9).
  CACHE_TREE          per-solve FmmPlan reuse for the budget-0 rows. MUST match
                      what rotor_hover_solver_phase2.jl does for its uncached
                      configs (currently false = tree rebuilt per apply).

## Interruption transparency (Ryan's ruling, 2026-08-24)

Every candidate the tuner measures is appended and flushed to

    tune_trace_<rung>_b<budget>.csv

beside tune_phase2.csv. On restart the trace is replayed into
`tune_fmm_perturb`'s memo, so a resumed descent walks the identical path at no
measurement cost and continues from where the interrupted one stopped: a
standby preemption costs ONE candidate, not the whole descent (4-8 h at R5).
`--requeue` alone therefore suffices, and TUNE_MAX_SECONDS is back to being a
backstop rather than a preemption mitigation.

A memoized `t` is a TIMING, so the trace records its provenance (rung, budget,
reps, abandon factor, hardware tag, fm commit) and REFUSES to load under a
mismatch rather than silently corrupting the descent. To restart a descent
fresh, move the trace file aside. Note the replayed timings come from the
pre-preemption node; jobs are pinned to spec-identical zen3 hardware and the
winner is re-timed here and again by rotor_hover_solver_phase2.jl, so this can
only affect which path the descent walks, never a published number.

Appends tune_phase2.csv under benchmark/results/phase2/<mode>[/<rung>].
Deliberately a NEW filename: the old tune_cached.csv holds accuracy-only rows
from the superseded objective and must not be silently mixed with these.
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

target_rel = 1e-6
safety = 0.1
include(joinpath(@__DIR__, "phase1_knobs.jl"))

import FLOWPanel.FastMultipole as FM

isnan(rms_b) && error("SKIP_B=1 is set: this driver SOLVES and needs the " *
                      "frozen b. Unset SKIP_B.")

# Descent starting point. It only has to be FEASIBLE and error-satisfying —
# the descent finds the minimum from there — so an untuned rung falls back to
# the reference triple rather than failing. A small leaf is the safe default:
# it fits every memory budget, so budget-0 and budget-capped sweeps both start
# legally.
const REF_START = (17, 0.5, 21)
p_0, mac_0, leaf_0 = haskey(TUNED, rung) ? TUNED[rung] : REF_START
haskey(TUNED, rung) || println("no TUNED[$rung]; starting the descent from " *
    "the reference triple $REF_START")

# ---- explicit descent SEED (2026-09-05, Ryan: LineGauss relaunch) -----------
# TUNED/tune.csv carries the PHASE 1 accuracy-stage knobs, which are a poor
# start for the cost objective (the R1-R3 hardcoded fallbacks are flagged
# ~2x too large in leaf, above). The 2026-08-26..29 Gaussian Phase 2 descents
# measured the real cost optima, and they split HARD by cache regime:
#   uncached (budget 0): leaf -> 6      (R4 6, R5 6, R3 3; MAC 0.60-0.65)
#   cached  (budget >0): leaf -> 32     (R4 32/32, R5 32/32/32; MAC 0.50-0.60)
# A single seed cannot serve both (5x apart), so budget 0 takes its own.
# These are SEEDS ONLY: the descent still searches, and the BC <= target
# certification gate still decides admissibility, so a stale or wrong seed
# costs optimality, never correctness. Recorded in the notes column so a row
# never hides where its descent started.
#
# Format: "P:MAC:LEAF" (e.g. TUNE_SEED=15:0.55:32). Unset => the TUNED/
# REF_START behaviour above, unchanged.
function _parse_seed(name)
    raw = get(ENV, name, "")
    isempty(raw) && return nothing
    parts = split(raw, r"[:,]")
    length(parts) == 3 || error("$name must be \"P:MAC:LEAF\", got \"$raw\"")
    return (parse(Int, parts[1]), parse(Float64, parts[2]), parse(Int, parts[3]))
end
const TUNE_SEED    = _parse_seed("TUNE_SEED")
const TUNE_SEED_B0 = _parse_seed("TUNE_SEED_B0")
seed_for(bgib) = (bgib == 0 && TUNE_SEED_B0 !== nothing) ? TUNE_SEED_B0 :
                 (TUNE_SEED !== nothing ? TUNE_SEED : (p_0, mac_0, leaf_0))
TUNE_SEED === nothing && TUNE_SEED_B0 === nothing ||
    println("descent seed: budget>0 -> $(seed_for(16)), budget 0 -> $(seed_for(0))")

# ---- knobs ------------------------------------------------------------------
const _REPS_DEFAULT = rung in ("R1", "R2", "R3", "R4") ? 5 : 2
tune_reps = parse(Int, get(ENV, "TUNE_REPS", string(_REPS_DEFAULT)))
tune_reps <= 5 || error("TUNE_REPS must be <= 5 (Ryan's cap); got $tune_reps")
abandon_factor = parse(Float64, get(ENV, "TUNE_ABANDON_FACTOR", "1.3"))
# Raised 5x on 2026-08-25 (Ryan's call). The previous ladder was inherited
# from the Phase 1 APPLY-PROXY tuner, where it was sized in Phase 16 to sit
# "below the typical preemption window" — a salvage device from the era when
# the tuner had no memo and a preemption cost a whole budget's descent. Phase
# 16's candidate-level trace removed that premise (a preemption now costs ONE
# candidate), but the numbers were never re-sized for a REAL-SOLVE candidate,
# which is a full solve x TUNE_REPS rather than one cheap apply.
#
# MEASURED consequence, 2026-08-25: every descent that had landed by then hit
# the timer instead of converging — R1 b0 at 6 candidates, R1 b16 at 16, R2 b0
# at 6, R3 b0 at 5 — against a descent that admits up to max_iters=20 ACCEPTED
# moves plus every rejected candidate it explores. None was a tuned optimum.
#
# Safe to change mid-campaign: max_seconds is NOT part of the trace provenance
# guard (that is rung/budget/reps/abandon_factor/hardware_tag), so existing
# traces still replay and a relaunch resumes rather than restarts.
const _MAXSEC = Dict("R1" => 18000.0, "R2" => 36000.0, "R3" => 72000.0,
                     "R4" => 108000.0, "R5" => 216000.0, "R6" => 432000.0,
                     "R7" => 432000.0)
max_seconds = parse(Float64, get(ENV, "TUNE_MAX_SECONDS", string(_MAXSEC[rung])))
max_build_time = parse(Float64, get(ENV, "NFCACHE_MAX_BUILD_TIME", "Inf"))
cache_tree_uncached = get(ENV, "CACHE_TREE", "0") == "1"
# Leaves scanned to price a rung's memory FLOOR — the cheapest configuration it
# admits — so an impossible budget is reported up front rather than discovered
# as a descent failure. A scan rather than a single leaf because cache bytes
# fall with leaf while the FMM plan grows with it, so the total has an interior
# minimum (see the floor check for the measured numbers).
const FLOOR_LEAVES = [parse(Int, x) for x in
                      split(get(ENV, "FLOOR_LEAVES", "2:4:9:21"), r"[:,]")]
# TOTAL memory budgets in GiB (body + solver state + plan + cache), NOT cache
# size. Ryan 2026-08-24: express the ladder as the memory of machines people
# actually have, so the axis answers "what can you solve on a laptop / a
# workstation / a supercomputer node" rather than sweeping arbitrary powers of
# two:
#
#     0    cache FORBIDDEN — the matrix-free method endpoint, not a machine.
#          Memory is unconstrained here; this is what tunes the uncached
#          configs (krylov_*, fgs, fgmres_fgs), which otherwise would have no
#          tuned knobs at all, since at R1-R5 every real machine can afford a
#          cache and no budget would ever return an uncached winner.
#     16   laptop
#     128  workstation
#     500  supercomputer node (the pinned --mem=500G on 512 GB zen3)
#
# What the classes MEAN varies by rung, which is the point (W0 map, ledger.md):
#   R1-R3  all three saturate — the whole near field fits anywhere
#   R4-R5  laptop constrained, workstation and node saturate
#   R6-R7  laptop cannot cache AT ALL (floors 24 and 79 GiB), workstation
#          constrained, node saturates
const MACHINE_LABEL = Dict(0.0 => "matrix-free endpoint (cache forbidden)",
                           16.0 => "laptop", 128.0 => "workstation",
                           500.0 => "supercomputer node")
budgets_gib = [parse(Float64, x) for x in
               split(get(ENV, "MEM_BUDGETS", "0:16:128:500"), r"[:,]")]

n = rotor.ncells
dense_bytes = 8 * n^2
println("\n$rung: $n panels, dense 8N^2 = $(round(dense_bytes/1024^3; digits=3)) GiB")
println("budgets (GiB): $budgets_gib   reps=$tune_reps  " *
        "abandon=$abandon_factor  max_seconds=$max_seconds")

outdir = joinpath(@__DIR__, "results", "phase2", banner.threading_mode)
get(ENV, "PER_RUNG_DIR", "0") == "1" && (outdir = joinpath(outdir, rung))
mkpath(outdir)
csv_path = joinpath(outdir, "tune_phase2.csv")
header = "rung,mesh_file,n_panels,mem_budget_gib,cached,expansion_order," *
    "multipole_acceptance,leaf_size,t_solve_warm,niter,bc_rel_l2_certified," *
    "bc_certified,t_cache_build,cache_bytes,mem_total_predicted," *
    "mem_total_measured,dense_bytes,cache_capped," *
    "tune_reps,tune_abandon_factor,tune_timed_out,t_tune,n_candidates," *
    "n_abandoned,cache_tree,threading_mode,julia_threads,blas_threads," *
    "commit,fm_commit,date,hardware_tag,notes"
fresh = !isfile(csv_path) || filesize(csv_path) == 0
if !fresh
    first(eachline(csv_path)) == header || error(
        "$csv_path has a pre-W2 schema; move it aside before appending")
end

# ---- resume: (rung, mem_budget_gib) is the row identity ---------------------
# standby QOS + --requeue restarts a preempted job from scratch, which would
# otherwise append duplicates of whatever had already landed.
landed = Set{String}()
if !fresh
    for line in Iterators.drop(eachline(csv_path), 1)
        isempty(strip(line)) && continue
        c = split(line, ",")
        length(c) >= 4 && String(c[1]) == rung && push!(landed, String(c[4]))
    end
end

io = open(csv_path, "a")
fresh && println(io, header)
hardware_tag = get(ENV, "HARDWARE_TAG", banner.hardware_tag)

# ---- candidate-level resume: the tuning trace ------------------------------
# Ryan 2026-08-24: "checkpoint every benchmark during tuning so a re-launched
# job picks up where the interrupted one stopped, then just use --requeue."
#
# The row-level resume above is per (rung, budget), so a preemption inside a
# descent throws away that whole descent — hours at R4/R5, repeatedly, under
# standby QOS. Every candidate `tune_fmm_perturb` measures is therefore
# appended here and flushed immediately. On restart the trace is replayed into
# the tuner's memo: because the neighbour stencil and the acceptance test are
# deterministic functions of memo lookups, the descent walks the identical path
# at no measurement cost and continues from the point of interruption.
#
# HONEST LIMITATION (ledger): a replayed `t` was measured on the pre-preemption
# node. All jobs are pinned --constraint=zen3 --exclusive --mem=500G on
# spec-identical hardware, and the winner is independently re-timed below (and
# again by rotor_hover_solver_phase2.jl), so this can only affect WHICH PATH the
# descent walks, never a published number.
const TRACE_HEADER = "expansion_order,multipole_acceptance,leaf_size,t," *
    "success,abandoned,rung,mem_budget_gib,tune_reps,tune_abandon_factor," *
    "hardware_tag,fm_commit,filament_reg"
trace_path(bgib) = joinpath(outdir, "tune_trace_$(rung)_b$(bgib).csv")

_tnum(x) = @sprintf("%.17g", x)   # full precision: a truncated t would make the
                                  # replayed descent differ from the original

"""
Rebuild a memo from `path`, refusing any trace whose provenance disagrees with
this run. A memoized `t` is a TIMING: replaying one measured under different
reps, a different abandonment rule or on different hardware would silently
corrupt the descent, so a mismatch is an error naming the file, never a
silently-ignored line.
"""
function load_trace(path, bgib)
    memo = Dict{Any, @NamedTuple{t::Float64, success::Bool, abandoned::Bool}}()
    isfile(path) && filesize(path) > 0 || return memo
    lines = collect(eachline(path))
    strip(lines[1]) == TRACE_HEADER || error(
        "$path has an unrecognised trace schema; move it aside before resuming")
    expect = (rung, string(bgib), string(tune_reps), string(abandon_factor),
              hardware_tag, string(pnl.FILAMENT_REGULARIZATION[]))
    for (ln, line) in enumerate(Iterators.drop(lines, 1))
        isempty(strip(line)) && continue
        c = String.(split(line, ","))
        length(c) == 13 || error("$path line $(ln+1): expected 13 fields, " *
                                 "got $(length(c))")
        got = (c[7], c[8], c[9], c[10], c[11], c[13])
        got == expect || error(
            "$path line $(ln+1): trace provenance disagrees with this run — " *
            "(rung, mem_budget_gib, tune_reps, tune_abandon_factor, " *
            "hardware_tag, filament_reg) is $got in the trace but $expect now. " *
            "A memoized "*
            "t is a timing and must not be replayed across a change of " *
            "reps, abandonment rule or hardware. Move the file aside to " *
            "start this descent fresh.")
        c[12] == banner.fm_commit || @warn "$path line $(ln+1): trace was " *
            "written under FastMultipole $(c[12]) but this run is " *
            "$(banner.fm_commit); the replayed timings are from a different " *
            "build of the objective"
        k = (parse(Int, c[1]), round(parse(Float64, c[2]); digits=3),
             parse(Int, c[3]))
        memo[k] = (; t=parse(Float64, c[4]), success=parse(Bool, c[5]),
                     abandoned=parse(Bool, c[6]))
    end
    return memo
end

bc_check(x) = bc_error!(rotor, x; rms_b, target_rel, safety,
                        max_expansion_order=20, multipole_acceptance=mac_0,
                        leaf_size=leaf_0, backend=:fmm)

# ---- memory accounting -----------------------------------------------------
# Ryan 2026-08-24: "memory rungs should correspond to the total memory
# available - not just that used by the cache." So a budget is what a NODE must
# hold, and the cache is only one term in it. Counted:
#
#   body            mesh/strength/field arrays (Base.summarysize)
#   solver state    Krylov workspace, preconditioner, scratch, histories
#                   (`solver_state_bytes`, the ruling-8 comparable metric -
#                   note it EXCLUDES the body by construction, hence the
#                   separate term)
#   plan            FMM trees, buffers, interaction lists - grows as leaf
#                   SHRINKS, so it partly offsets the cache's growth
#   cache           the near-field blocks, when one is used
#
# The budget is enforced against the PREDICTED total, because the plan is live
# during a solve even when `cache_tree=false` drops it afterwards and it would
# not show up in a post-hoc measurement. Measured total is reported alongside.
#
# This also puts `backslash_ldiv` (dense G, 8N^2) on the SAME axis: the
# question becomes "given a node with X GB, which solver is fastest", which is
# the decision a user actually faces.
const BODY_BYTES = Base.summarysize(rotor)
println("body = $(round(BODY_BYTES/1024^3; digits=3)) GiB")

"""
Plan size and prospective cache size at these knobs, WITHOUT building a cache.
One plan build per candidate, reused for both numbers, so an over-budget
candidate costs a tree build rather than a solve (or an OOM).
"""
function plan_and_cache_bytes(P, mac, leaf)
    plan = FM.FmmPlan((rotor,), (rotor,); scalar_potential=true, gradient=false,
        hessian=false, expansion_order=P, multipole_acceptance=mac,
        leaf_size_source=leaf, shrink=true)
    plan_bytes = Base.summarysize(plan)
    est = FM.estimate_nearfield_cache(plan.target_tree, plan.source_tree,
        plan.direct_list, plan.derivatives_switches, (rotor,); sample=false)
    plan = nothing; GC.gc()
    return plan_bytes, est.bytes
end

"""
Build the fastest configuration that FITS `budget_bytes` at these knobs and
return `(t_warm, success, info)`.

Whether the near-field cache is used is an OUTCOME, not an input: it is used
iff the whole configuration still fits the total-memory budget with it. So a
rung or a budget that admits no cache is tuned UNCACHED automatically, which is
exactly what R6/R7 need (their cache FLOOR is 24 and 78 GiB respectively) —
Ryan: "if memory bandwidths prohibit a cache, then that particular run will
need to run without a cache... we'll need to tune those rungs that don't admit
a cache just so - without a cache."

`t_warm` is min-of-`reps` over solves that start from a COLD solution (zeroed
strengths, no warm start, so the iteration count is honest) but, when cached, a
WARM cache: the plan and its blocks were built by an untimed warm-up solve
beforehand. That is Ryan's ruling implemented structurally — the build is
excluded because it never enters a timed region, not because it was subtracted.

`success` is the certified BC verdict on the converged solution — the same
reference-free criterion Phase 1 uses, and strictly stronger than the per-apply
error tolerance the old `fmm!` proxy used.
"""
function solve_cost(P, mac, leaf, budget_bytes; allow_cache=true)
    plan_bytes, cache_bytes = plan_and_cache_bytes(P, mac, leaf)
    backend = pnl.FastMultipoleBackend(; expansion_order=P,
        multipole_acceptance=mac, leaf_size=leaf)

    # base solver state (workspace + preconditioner); the plan is built lazily,
    # so this does not include it yet
    probe = pnl.KrylovSolver(rotor; method=:gmres, itmax=500, atol=1e-14,
        rtol=target_rel, memory=50, backend)
    base_bytes = solver_state_bytes(probe)
    probe = nothing; GC.gc()

    fixed = BODY_BYTES + base_bytes + plan_bytes
    cached = allow_cache && fixed + cache_bytes <= budget_bytes
    predicted = fixed + (cached ? cache_bytes : 0)
    if predicted > budget_bytes
        return (Inf, false, (; cached=false, cache_bytes=0, plan_bytes,
            base_bytes, predicted, measured=-1, t_build=NaN, niter=-1, bc=NaN,
            reason="uncached configuration alone needs " *
                   "$(round(predicted/1024^3; digits=2)) GiB"))
    end

    solver = pnl.KrylovSolver(rotor; method=:gmres, itmax=500, atol=1e-14,
        rtol=target_rel, memory=50, backend,
        cache_tree = cached ? true : cache_tree_uncached,
        cache_nearfield = cached, persistent_plan = cached,
        nearfield_cache_max_bytes = cached ?
            max(cache_bytes, budget_bytes - fixed) : typemax(Int),
        nearfield_cache_max_build_time = max_build_time)

    # untimed warm-up: when cached, this is the solve that builds plan + cache
    reset_cold!()
    t_build = 0.0
    c_bytes = 0
    try
        pnl._solve!(rotor, solver)
    catch e
        solver = nothing; GC.gc()
        return (Inf, false, (; cached, cache_bytes, plan_bytes, base_bytes,
            predicted, measured=-1, t_build=NaN, niter=-1, bc=NaN,
            reason=sprint(showerror, e)))
    end
    if cached
        slot = solver.kop.plan_slot[]
        cache = slot === nothing ? nothing : slot[1].nearfield_cache[]
        cache === nothing && error("cache_nearfield=true but no cache on the " *
            "plan after the warm-up solve — the build path did not run")
        t_build = cache.build_time
        c_bytes = cache.bytes
    end
    niter = solver.niter
    converged = solver.solved
    e = bc_check(vec(rotor.strength[:, solution_column]))
    certified = e.error_success && e.rel_l2 <= target_rel && converged
    # measured total; when the plan is not retained (cache_tree=false) it is
    # still live during a solve, so add it back rather than under-reporting
    measured = BODY_BYTES + solver_state_bytes(solver) +
               ((cached || cache_tree_uncached) ? 0 : plan_bytes)

    t_min = Inf
    for _ in 1:tune_reps
        reset_cold!()
        t_min = min(t_min, @elapsed pnl._solve!(rotor, solver))
    end
    solver = nothing; GC.gc()
    return (t_min, certified,
            (; cached, cache_bytes=c_bytes, plan_bytes, base_bytes, predicted,
               measured, t_build, niter, bc=e.rel_l2,
               reason = certified ? "" :
                        (!converged ? "krylov did not converge" :
                         !e.error_success ? "fmm error not certified" :
                         "bc_rel_l2 $(e.rel_l2) > $target_rel")))
end

# ---- sweep ------------------------------------------------------------------
prev_signature = Ref{Any}(nothing)
for bgib in budgets_gib
    tag = string(bgib)
    if tag in landed
        println("\n=== budget $(bgib) GiB: already landed, skipping ===")
        continue
    end
    # budget 0 is the FORCED-UNCACHED sentinel, not a 0-byte machine: the
    # cache is forbidden and memory is otherwise unconstrained.
    allow_cache = bgib > 0
    budget_bytes = allow_cache ? round(Int, bgib * 1024^3) : typemax(Int)
    label = get(MACHINE_LABEL, bgib, "")
    println("\n=== " * (allow_cache ? "TOTAL memory budget $(bgib) GiB" :
            "matrix-free endpoint (cache forbidden)") *
            (isempty(label) || !allow_cache ? "" : " — $label") * " ===")

    # Floor check: can the UNCACHED configuration fit at all? If not, no
    # configuration exists at this budget and the descent would otherwise fail
    # from deep inside as "the cost function rejected the starting parameters".
    if allow_cache
        # Floor check: is ANY configuration affordable at this budget?
        #
        # The floor is NOT at a fixed leaf. Cache bytes fall as the leaf
        # shrinks (MEASURED R6: 23.8 GiB at leaf 9 -> 19.3 at leaf 2, still
        # falling but clearly asymptoting; R7: 78.2 -> 67.0), while the FMM
        # plan GROWS as the leaf shrinks. Total memory therefore has an
        # interior minimum, and pricing it at one arbitrary leaf — as this
        # check originally did at leaf 9 — can both overstate the floor (by
        # missing a cheaper small-leaf point) and understate it (by ignoring
        # plan growth). So scan and take the minimum of the actual total.
        best_unc, best_cac, best_leaf = Inf, Inf, 0
        for lf in FLOOR_LEAVES
            pb, cb = plan_and_cache_bytes(p_0, mac_0, lf)
            probe = pnl.KrylovSolver(rotor; method=:gmres, itmax=500,
                atol=1e-14, rtol=target_rel, memory=50,
                backend=pnl.FastMultipoleBackend(; expansion_order=p_0,
                    multipole_acceptance=mac_0, leaf_size=lf))
            unc = BODY_BYTES + solver_state_bytes(probe) + pb
            probe = nothing; GC.gc()
            unc < best_unc && (best_unc = unc)
            if unc + cb < best_cac
                best_cac = unc + cb
                best_leaf = lf
            end
        end
        if best_unc > budget_bytes
            msg = "budget $(bgib) GiB is BELOW this rung's floor: even the " *
                  "cheapest UNCACHED configuration needs " *
                  "$(round(best_unc/1024^3; digits=2)) GiB " *
                  "(body + solver state + plan), minimised over leaves " *
                  "$(FLOOR_LEAVES)"
            println("  $msg")
            println(io, join(_csv_cell.([rung, msh_name, n, bgib, false,
                "", "", "", "", "", "", false, "", "", best_unc, "",
                dense_bytes, true, tune_reps, abandon_factor, false, 0.0, 0, 0,
                cache_tree_uncached, banner.threading_mode,
                banner.julia_threads, banner.blas_threads, banner.commit,
                banner.fm_commit, time_string(), hardware_tag, msg]), ","))
            flush(io)
            continue
        end
        println("  floor: uncached $(round(best_unc/1024^3; digits=2)) GiB, " *
                "cheapest cached $(round(best_cac/1024^3; digits=2)) GiB " *
                "(at leaf $best_leaf) — " *
                (best_cac > budget_bytes ?
                 "NO cache fits this budget, tuning will run matrix-free" :
                 "a cache fits"))
    end

    function cost(; expansion_order, multipole_acceptance, leaf_size_source)
        t, ok, info = solve_cost(expansion_order, multipole_acceptance,
                                 leaf_size_source, budget_bytes; allow_cache)
        println("    P=$expansion_order MAC=$multipole_acceptance " *
                "leaf=$leaf_size_source -> " *
                (isfinite(t) ?
                 "t=$(round(t; digits=3)) s niter=$(info.niter) " *
                 "$(info.cached ? "CACHED" : "uncached") " *
                 "mem=$(round(info.predicted/1024^3; digits=2))/$(bgib) GiB " *
                 "bc=$(info.bc)" :
                 "INFEASIBLE ($(info.reason))"))
        return (t, ok)
    end

    tpath = trace_path(bgib)
    memo = load_trace(tpath, bgib)
    n_replayed = length(memo)
    n_replayed > 0 && println("  resuming from $(basename(tpath)): " *
        "$n_replayed candidate(s) will be replayed from the trace")
    tfresh = !isfile(tpath) || filesize(tpath) == 0
    trace_io = open(tpath, "a")
    tfresh && println(trace_io, TRACE_HEADER)
    # the checkpoint: one flushed line per freshly measured candidate, so a
    # preemption now costs at most the candidate in flight
    # string(), NOT _csv_cell: load_trace compares these against string(bgib)
    # etc., and _csv_cell's %.9g renders 0.0 as "0", which would never match
    # filament_reg (2026-09-05): the family changes the per-edge kernel cost
    # (LineGauss 4 erf + 1 exp vs Gaussian 1 expm1) AND radius_inflation, so it
    # changes the very objective these timings measure. It belongs in the HARD
    # guard alongside reps/abandon/hardware — without it a Gaussian trace would
    # replay silently into a LineGauss descent.
    fam_str = string(pnl.FILAMENT_REGULARIZATION[])
    provenance = "," * join(string.([rung, bgib, tune_reps, abandon_factor,
                                     hardware_tag, banner.fm_commit, fam_str]), ",")
    any(occursin(",", x) for x in (rung, hardware_tag, banner.fm_commit, fam_str)) &&
        error("a comma in rung/hardware_tag/fm_commit/filament_reg would " *
              "corrupt the trace")
    on_measure = function (P, mac, leaf, t, success, abandoned)
        # a per-system leaf tuple would not round-trip through the trace's one
        # integer column, and the memo key would then never match on replay
        leaf isa Integer || error("tuning trace expects a scalar leaf size; " *
                                  "got $(typeof(leaf)) = $leaf")
        println(trace_io, "$P,$(_tnum(mac)),$leaf," *
                "$(_tnum(t)),$success,$abandoned" * provenance)
        flush(trace_io)
    end

    local tuned, hist, tinfo, t_tune
    try
        seed_p, seed_mac, seed_leaf = seed_for(bgib)
        t_tune = @elapsed ((tuned, hist, tinfo) = FM.tune_fmm_perturb(
            (rotor,), (rotor,); expansion_order=seed_p,
            multipole_acceptance=seed_mac, leaf_size_source=seed_leaf,
            reps=tune_reps, abandon_factor, max_seconds, cost, memo, on_measure,
            verbose=true))
    catch e
        @warn "budget $(bgib) GiB: tuning FAILED" exception=(e, catch_backtrace())
        close(trace_io)
        continue
    end
    close(trace_io)

    leaf_w = maximum(tuned.leaf_size_source)
    # re-measure the winner so every reported number belongs to the point being
    # published, not to whatever candidate happened to be evaluated last
    t_w, ok_w, info_w = solve_cost(tuned.expansion_order,
        tuned.multipole_acceptance, leaf_w, budget_bytes; allow_cache)

    # capped: was the descent stopped by the BUDGET rather than by cost? True
    # iff the next leaf up the ladder would not fit as a cached configuration
    # while this one does. Reported so a capped winner is never read as an
    # unconstrained optimum.
    cache_capped = false
    if info_w.cached
        leaf_next = max(1, round(Int, leaf_w * 1.5))
        pb_n, cb_n = plan_and_cache_bytes(tuned.expansion_order,
            tuned.multipole_acceptance, leaf_next)
        cache_capped = BODY_BYTES + info_w.base_bytes + pb_n + cb_n > budget_bytes
    end

    println("  WINNER $(bgib) GiB$(isempty(label) ? "" : " ($label)"): " *
            "p=$(tuned.expansion_order) " *
            "mac=$(tuned.multipole_acceptance) leaf=$leaf_w " *
            "$(info_w.cached ? "CACHED" : "UNCACHED")  " *
            "t_warm=$(round(t_w; digits=3)) s niter=$(info_w.niter) " *
            "mem pred=$(round(info_w.predicted/1024^3; digits=2)) " *
            "meas=$(round(info_w.measured/1024^3; digits=2)) GiB " *
            "(cache $(round(info_w.cache_bytes/1024^3; digits=3)) GiB, " *
            "build $(round(info_w.t_build; digits=2)) s) " *
            "capped=$cache_capped certified=$ok_w timed_out=$(tinfo.timed_out)")
    tinfo.timed_out && @warn "budget $(bgib) GiB: descent TIMED OUT — the " *
        "returned knobs are best-so-far, NOT a converged minimum"
    ok_w || @warn "budget $(bgib) GiB: winner is NOT certified ($(info_w.reason))"

    # Budgets that produce an identical uncached winner are the SAME point;
    # flag it so the memory curve is not read as having independent evidence
    # at each budget.
    sig = (info_w.cached, tuned.expansion_order, tuned.multipole_acceptance, leaf_w)
    dup = sig == prev_signature[]
    prev_signature[] = sig
    dup && println("  (identical configuration to the previous budget — this " *
                   "budget adds no new operating point)")

    println(io, join(_csv_cell.([rung, msh_name, n, bgib, info_w.cached,
        tuned.expansion_order, tuned.multipole_acceptance, leaf_w,
        t_w, info_w.niter, info_w.bc, ok_w, info_w.t_build, info_w.cache_bytes,
        info_w.predicted, info_w.measured, dense_bytes, cache_capped,
        tune_reps, abandon_factor, tinfo.timed_out, t_tune, tinfo.n_candidates,
        tinfo.n_abandoned, info_w.cached ? true : cache_tree_uncached,
        banner.threading_mode, banner.julia_threads, banner.blas_threads,
        banner.commit, banner.fm_commit, time_string(), hardware_tag,
        (isempty(label) ? "" : "$label; ") *
        (n_replayed > 0 ? "resumed_from_trace n_replayed=$n_replayed; " : "") *
        "seed=$(seed_for(bgib)); " *
        "real-solve objective; TOTAL-memory budget (body+solver state+plan" *
        (info_w.cached ? "+cache" : "") * "); cache used iff it fits; " *
        "warm cache, cold solution; build excluded by construction" *
        (dup ? "; DUPLICATE of previous budget" : "")]), ","))
    flush(io)
end
close(io)
println("\nwrote $csv_path")
