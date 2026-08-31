#=##############################################################################
BRAINSTORM 021 Phase 2a — setup vs per-step cost, componentized, at the
Phase-1 frozen per-rung settings (knobs read via phase1_knobs.jl exactly as
the Phase 1 table does; ruling-3 shared sets).

Adds over solvetable.csv (ruling 7 / Stage 2a):
  - componentized one-time setup: dense assembly vs LU (Backslash; the LU
    component is DERIVED as ctor_total - scratch_assembly to avoid paying the
    O(N^3) factorization twice at R5), preconditioner build (Jacobi ctor; ILU
    tree/lists/assembly/factorization from its stats dict; FGSPreconditioner
    ctor), FGS leaf-LU cache build time/bytes;
  - memory: solver_state_bytes + @allocated on a cold solve;
  - convergence-history sidecars with wall-clock timestamps for every
    iterative config (equal-time-per-iteration validation happens in the
    analysis stage, ruling 4).
Per-step cost = isolated cold solve, adaptive min-of-k (ruling 12: strictly
cold — warmstart measurements are Phase 3's). The full-timestep share leg
lives in rotor_hover_solver_unsteady.jl. No ruling-4 dual rows here (cost
table, not the Phase-1 accuracy report).

Run (local smoke; published rows are HPC-only per ruling 5):
  RUNG=R1 CONFIGS=backslash_ldiv:fgs EXPECT_JULIA_THREADS=4 \
    THREADING_MODE=multi CACHE_B=1 \
    julia --project -t 4 benchmark/rotor_hover_solver_phase2.jl
CONFIGS is colon- or comma-separated. Appends phase2.csv under
benchmark/results/phase2/<mode>[/<rung> with PER_RUNG_DIR=1].
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

using LinearAlgebra: lu!

target_rel = 1e-6
safety = 0.1
include(joinpath(@__DIR__, "phase1_knobs.jl"))

# ---- apply knobs -----------------------------------------------------------
# Phase 2 now OWNS tuning (Ryan 2026-08-24): Phase 1 was descoped to
# solver-to-solver verification, so the shared apply knobs come from
# tune_phase2.csv (rotor_hover_solver_phase2_tune.jl), keyed by memory budget.
# Budget 0 = the uncached configs; every other budget is one point on the
# memory axis for the *_nfcache configs. Latest row per (rung, budget) wins.
function phase2_knobs(budget_gib)
    knobs_csv = joinpath(outdir, "tune_phase2.csv")
    isfile(knobs_csv) || return nothing
    lines = readlines(knobs_csv)
    cols = Dict(String(c) => i for (i, c) in enumerate(split(lines[1], ",")))
    out = nothing
    for l in lines[2:end]
        c = split(l, ",")
        length(c) >= length(cols) || continue
        c[cols["rung"]] == rung || continue
        parse(Float64, c[cols["mem_budget_gib"]]) == budget_gib || continue
        isempty(strip(c[cols["leaf_size"]])) && continue   # infeasible-budget row
        out = (p=parse(Int, c[cols["expansion_order"]]),
               mac=parse(Float64, c[cols["multipole_acceptance"]]),
               leaf=parse(Int, c[cols["leaf_size"]]),
               capped=lowercase(strip(c[cols["cache_capped"]])) == "true",
               # A budget's winner may legitimately be UNCACHED (the cache did
               # not fit the TOTAL memory budget). Follow the tuned row rather
               # than assuming a *_nfcache config is cached.
               cached=lowercase(strip(c[cols["cached"]])) == "true")
    end
    return out
end


bc_fmm(x) = bc_error!(rotor, x; rms_b, target_rel, safety,
                      max_expansion_order=20, multipole_acceptance=mac_t,
                      leaf_size=leaf_t, backend=:fmm)

want = split(get(ENV, "CONFIGS",
    "backslash_ldiv,krylov_gmres,krylov_jacobi,krylov_ilu,fgs,fgmres_fgs," *
    "krylov_gmres_nfcache,krylov_ilu_nfcache,fgmres_fgs_nfcache"),
    r"[:,]")

# results go under phase2/ (phase1_case's results_dir/knobs_dir stay pointed
# at phase1/ for the frozen knob CSVs)
outdir = joinpath(@__DIR__, "results", "phase2", banner.threading_mode)
get(ENV, "PER_RUNG_DIR", "0") == "1" && (outdir = joinpath(outdir, rung))
mkpath(outdir)
write(joinpath(outdir, "banner.txt"), banner.text * "\n")

# Memory-budget ladder (TOTAL memory: body + solver state + plan + cache),
# expressed as the machines people actually have (Ryan 2026-08-24):
# 16 = laptop, 128 = workstation, 500 = supercomputer node. Must match what
# rotor_hover_solver_phase2_tune.jl swept. Budget 0 is NOT in this list: it is
# the forced-uncached endpoint, which the non-nfcache configs already occupy.
const MACHINE_LABEL = Dict(16.0 => "laptop", 128.0 => "workstation",
                           500.0 => "supercomputer node")
nfcache_budgets = [parse(Float64, x) for x in
                   split(get(ENV, "MEM_BUDGETS", "16:128:500"), r"[:,]")]

# The UNCACHED apply knobs come from the tuner's budget-0 row — the explicit
# cache-forbidden point. It must exist: at R1-R5 every real machine can afford
# a cache, so no positive budget ever returns an uncached winner.
_k0 = phase2_knobs(0.0)
if _k0 === nothing
    haskey(TUNED, rung) || error("no budget-0 (cache-forbidden) row in " *
        "tune_phase2.csv for $rung and no TUNED fallback — run " *
        "rotor_hover_solver_phase2_tune.jl with 0 in MEM_BUDGETS first")
    @warn "no budget-0 tune_phase2.csv row for $rung; falling back to " *
          "TUNED[$rung] = $(TUNED[rung]). These are PRE-FIX accuracy-only " *
          "knobs, not valid campaign knobs."
    p_t, mac_t, leaf_t = TUNED[rung]
else
    p_t, mac_t, leaf_t = _k0.p, _k0.mac, _k0.leaf
end
# Total memory reported per row = body + solver state (ruling-8 metric, which
# excludes the body by construction). Same accounting as the tuner's budgets.
const BODY_BYTES = Base.summarysize(rotor)
println("uncached apply knobs: p=$p_t mac=$mac_t leaf=$leaf_t")
println("memory-budget ladder (TOTAL GiB): " *
        join(("$b" * (haskey(MACHINE_LABEL, b) ? " ($(MACHINE_LABEL[b]))" : "")
              for b in nfcache_budgets), ", "))

header = "rung,mesh_file,n_panels,config,row_kind," *
    "t_setup_total,t_setup_assembly,t_setup_factorization," *
    "t_setup_tree,t_setup_lists,t_setup_precond,t_setup_leaf_lu," *
    "leaf_lu_bytes,t_solve_min,k_reps,niter,solve_alloc_bytes," *
    "solver_state_bytes,bc_rel_l2_certified,bc_certified,solver_knobs," *
    "apply_p,apply_mac,apply_leaf,mem_budget_gib,nfcache_build_time," *
    "nfcache_bytes,nfcache_capped,mem_total_measured,rms_b,radius_tol," *
    "threading_mode," *
    "julia_threads,blas_threads,commit,fm_commit,date,notes"
csv_path = joinpath(outdir, "phase2.csv")
fresh = !isfile(csv_path) || filesize(csv_path) == 0
if !fresh
    first(eachline(csv_path)) == header || error(
        "$csv_path has a pre-W3 schema (no mem_budget/nfcache_build_time " *
        "columns); move it aside before appending")
end
io = open(csv_path, "a")
fresh && println(io, header)

# ---- Resume support (2026-08-24) --------------------------------------------
# Jobs run on the preemptible `standby` QOS with --requeue, so a preempted job
# restarts from scratch. This driver APPENDS; without a skip, a restart would
# duplicate every row it had already written. Row identity is
# (rung, config, row_kind) and each row is flushed as its measurement finishes,
# so a row present in the CSV means that measurement is DONE. Columns 1/4/5
# never contain commas (only the trailing notes cell is ever quoted).
# A *_nfcache config now emits ONE ROW PER MEMORY BUDGET, so the budget is
# part of the row identity; without it a requeue would skip every budget after
# the first one landed.
landed = Set{Tuple{String,String,String}}()
_bcol = findfirst(==("mem_budget_gib"), String.(split(header, ",")))
if !fresh
    for line in Iterators.drop(eachline(csv_path), 1)
        isempty(strip(line)) && continue
        c = split(line, ",")
        length(c) >= _bcol && String(c[1]) == rung &&
            push!(landed, (String(c[4]), String(c[5]), String(c[_bcol])))
    end
end
has_landed(config, row_kind="standard", budget="") =
    (config, row_kind, string(budget)) in landed
if !isempty(landed)
    println("RESUME: $(length(landed)) row(s) already landed for $rung — " *
            "skipping " * join(sort(["$a/$b" for (a, b) in landed]), ", "))
end

common_note = "cold isolated solve (ruling 12); t excludes N^2 source " *
              "assembly (restored); tree rebuild per Krylov apply included"

"""
Shared measurement block: cold instrumented history solve (when a recorder is
given), adaptive min-of-k cold solve, cold-solve allocations, certified BC,
CSV row + history sidecar. `setup` is a NamedTuple with any of
(total, assembly, factorization, tree, lists, precond, leaf_lu, leaf_lu_bytes)
— missing components print blank.
"""
function measure_and_emit!(config, solver, setup::NamedTuple, knobs, notes;
                           history_solve! = nothing, niter_from = nothing,
                           apply_knobs = (p_t, mac_t, leaf_t),
                           mem_budget = "", nfcache_build = "",
                           nfcache_bytes = "", nfcache_capped = "",
                           mem_total = "", row_kind = "standard")
    niter = -1
    if history_solve! !== nothing
        reset_cold!()
        hist = history_solve!()
        write_history_csv(joinpath(outdir, "history_$(config)_$(rung).csv"), hist)
        niter = length(hist)
    end
    t_solve, kk = adaptive_min_of_k(() -> pnl._solve!(rotor, solver);
                                    setup! = reset_cold!)
    reset_cold!()
    solve_alloc = @allocated pnl._solve!(rotor, solver)
    niter_from !== nothing && (niter = niter_from())
    x = copy(vec(rotor.strength[:, solution_column]))
    e = bc_fmm(x)
    g(k) = get(setup, k, "")
    println("  $config: setup=$(round(get(setup, :total, NaN); digits=2))s " *
            "solve=$(round(t_solve; digits=4))s (k=$kk) niter=$niter " *
            "alloc=$(solve_alloc) bc=$(e.rel_l2) (certified=$(e.error_success))")
    println(io, join(_csv_cell.([rung, msh_name, rotor.ncells, config,
        row_kind, g(:total), g(:assembly), g(:factorization), g(:tree),
        g(:lists), g(:precond), g(:leaf_lu), g(:leaf_lu_bytes),
        t_solve, kk, niter, solve_alloc, solver_state_bytes(solver),
        e.rel_l2, e.error_success, knobs, apply_knobs[1], apply_knobs[2],
        apply_knobs[3], mem_budget, nfcache_build, nfcache_bytes,
        nfcache_capped, mem_total, rms_b,
        pnl.FMM_RADIUS_TOL[], banner.threading_mode, banner.julia_threads,
        banner.blas_threads, banner.commit, banner.fm_commit, time_string(),
        notes]), ","))
    flush(io)
    return t_solve
end

# ---- backslash_ldiv: assembly measured on a scratch G, LU derived ----------
function run_backslash!()
    println("--- backslash_ldiv ---")
    G = zeros(rotor.ncells, rotor.ncells)
    t_assembly = @elapsed pnl._G!(G, rotor, rotor;
                                  core_size=rotor.core_size_panel)
    G = nothing; GC.gc()   # scratch freed before the ctor builds its own G
    t_total = @elapsed solver = pnl.Backslash(rotor)
    setup = (; total=t_total, assembly=t_assembly,
             factorization=max(t_total - t_assembly, 0.0))
    measure_and_emit!("backslash_ldiv", solver, setup,
        "dense LU (ctor=assembly+factorization)",
        "assembly timed on a scratch G; factorization = ctor_total - " *
        "assembly (derived, avoids double LU); " * common_note)
    solver = nothing; GC.gc()
end

# ---- Krylov trio at the shared tuned apply knobs ---------------------------
backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)
krylov_kw = (; itmax=500, atol=1e-14, rtol=target_rel, memory=50,
             backend=backend_apply, record_history=true)

# ---- cached-economics apply knobs for the *_nfcache variants ---------------
# Per Ryan's Phase 2 caching ruling (2026-08-24): the cached rows are measured
# WARM and the cache build is reported in its OWN column, never folded into
# t_solve.
#
# The previous comment here claimed the cache was per-SOLVE state because "a
# rotating rotor invalidates it every step". That is STALE. `persistent_plan`
# landed afterwards with rigid-motion support, and the Dirichlet
# scalar-potential operator and its cache blocks are EXACTLY invariant under
# rigid motion (FLOWPanel_solver.jl `persistent_plan` docstring;
# FastMultipole `transform_plan!`). A rotating rotor therefore reuses one cache
# for the whole run: the build is a genuine one-off, and charging it to every
# step would misprice it by the run length.
#
# So each *_nfcache config is measured as:
#   setup   : construct solver + ONE untimed warm-up solve, which is what
#             lazily builds the plan and its near-field blocks
#   t_solve : min-of-k solves from a COLD solution (reset_cold! zeroes the
#             strengths, warmstart is off, so niter is the honest cold count)
#             against a WARM cache
#   columns : nfcache_build_time / nfcache_bytes, separate
# The build is excluded by construction, not by subtraction.
#
# One row per memory budget: caching is a continuum between matrix-free and
# matrix-ful, and where a config lands on it is set by the memory it is allowed
# (see rotor_hover_solver_phase2_tune.jl for why the budget is an axis and not
# a constant).
nfcache_state = Dict{Float64,Any}()
function nfcache_setup(budget_gib)
    haskey(nfcache_state, budget_gib) && return nfcache_state[budget_gib]
    ck = phase2_knobs(budget_gib)
    ck === nothing && return nothing
    backend = pnl.FastMultipoleBackend(; expansion_order=ck.p,
        multipole_acceptance=ck.mac, leaf_size=ck.leaf)
    nfcache_state[budget_gib] = (backend, ck)
    return nfcache_state[budget_gib]
end

"""
Build the plan + near-field cache on `solver` with an untimed warm-up solve and
return `(t_warmup, cache)`. Must be called before any timed solve so that
t_solve is warm.
"""
function prime_nfcache!(solver)
    reset_cold!()
    t = @elapsed pnl._solve!(rotor, solver)
    slot = solver.kop.plan_slot[]
    slot === nothing && error("cache_nearfield/persistent_plan set but the " *
        "plan slot is empty after the warm-up solve")
    cache = slot[1].nearfield_cache[]
    cache === nothing && error("plan carries no near-field cache after the " *
        "warm-up solve — the build path did not run")
    return t, cache
end

function run_krylov!(name; budget_gib=0.0)
    nfcache = endswith(name, "_nfcache")
    ck = nothing; backend_nf = nothing
    if nfcache
        st = nfcache_setup(budget_gib)
        if st === nothing
            println("--- $name @ $(budget_gib) GiB --- SKIPPED (no tuned row " *
                    "for this budget; it may be below the rung's cache floor)")
            return
        end
        backend_nf, ck = st
        nfcache = ck.cached   # the tuner may have found no cache fits here
        println("--- $name @ $(budget_gib) GiB total (p=$(ck.p) mac=$(ck.mac) " *
                "leaf=$(ck.leaf), $(ck.cached ? "CACHED" : "UNCACHED — no " *
                "cache fits this budget")$(ck.capped ? ", CAPPED" : "")) ---")
    else
        println("--- $name ---")
    end
    t_precond = 0.0
    local P = nothing
    if name == "krylov_jacobi"
        t_precond = @elapsed P = pnl.FastMultipole.JacobiPreconditioner(
            (rotor,); cell_size=R/4)
    elseif name in ("krylov_ilu", "krylov_ilu_nfcache")
        t_precond = @elapsed P = pnl.ILUPreconditioner(rotor; leaf_size=10,
            multipole_acceptance=1.0, max_pattern_entries=8192 * rotor.ncells)
    end
    # persistent_plan keeps the plan AND its near-field blocks alive across
    # solves, so the build happens once in the warm-up below instead of inside
    # every timed solve (Ryan's ruling: build separate from per-step cost).
    # max_bytes is a loose guard here: the tuner already established that the
    # whole configuration fits `budget_gib`, so the cache alone cannot exceed it
    extra = nfcache ?
        (; backend=backend_nf, cache_nearfield=true, persistent_plan=true,
           nearfield_cache_max_bytes=round(Int, budget_gib * 1024^3)) :
        (ck === nothing ? (;) : (; backend=backend_nf))
    t_ctor = @elapsed solver = P === nothing ?
        pnl.KrylovSolver(rotor; method=:gmres, krylov_kw..., extra...) :
        pnl.KrylovSolver(rotor; method=:gmres, krylov_kw..., preconditioner=P,
                         extra...)
    t_warmup = 0.0; nfc = nothing
    nfcache && ((t_warmup, nfc) = prime_nfcache!(solver))
    setup = if P isa pnl.ILUPreconditioner
        st = P.stats
        (; total=t_precond + t_ctor + t_warmup, tree=st["tree_time"],
           lists=st["interaction_list_time"], assembly=st["assembly_time"],
           factorization=st["factorization_time"], precond=t_precond)
    else
        (; total=t_precond + t_ctor + t_warmup, precond=t_precond)
    end
    nfnote = nfcache ?
        "WARM cache: plan + near-field blocks built once in setup " *
        "(persistent_plan), so t_solve is the per-step cost with the build " *
        "EXCLUDED; the build is nfcache_build_time. " : ""
    note = nfcache ?
        nfnote * "cold isolated solve on a warm cache (ruling 12 applies to " *
        "the SOLUTION: reset_cold! zeroes strengths, warmstart off); " *
        "t excludes N^2 source assembly (restored)" :
        common_note
    measure_and_emit!(name, solver, setup,
        "rtol=1e-6;atol=1e-14;memory=50" *
        (nfcache ? ";cache_nearfield;persistent_plan" : ""),
        note; niter_from=() -> solver.niter,
        history_solve! = () -> (pnl._solve!(rotor, solver); solver.history),
        apply_knobs = nfcache ? (ck.p, ck.mac, ck.leaf) : (p_t, mac_t, leaf_t),
        mem_budget = ck === nothing ? "" : budget_gib,
        nfcache_build = nfcache ? nfc.build_time : "",
        nfcache_bytes = nfcache ? nfc.bytes : "",
        nfcache_capped = ck === nothing ? "" : ck.capped,
        mem_total = ck === nothing ? "" :
            BODY_BYTES + solver_state_bytes(solver))
    solver = nothing; P = nothing; GC.gc()
end

# ---- FGS family at the ruling-3 shared knob set ----------------------------
function fgs_shared_set()
    winner = stage3_winner()
    winner === nothing && error(
        "fgsprecond.csv has no winner for $rung — run Phase 1 Stages 2-3 first")
    sc = staircase_for(winner.p, winner.mac, winner.leaf, winner.inner)
    sc === nothing && error("No staircase rows for the winner knobs")
    i_cross = findfirst(t -> t[4] <= target_rel, sc)
    i_cross === nothing && error("Winner knobs never reach BC 1e-6 in the staircase")
    return winner, margin_tol(sc, i_cross)
end

function run_fgs!(winner, tol_abs)
    println("--- fgs (shared set) ---")
    t_total = @elapsed solver = pnl.FGSSolver(rotor;
        expansion_order=winner.p, multipole_acceptance=winner.mac,
        leaf_size=winner.leaf, inner_iterations=winner.inner,
        max_iterations=300, tolerance=tol_abs, rlx=1.0, shrink=true,
        recenter=false, reverse_pass=false, verbose=false)
    cache = solver.fgs.leaf_lu_cache
    setup = cache === nothing ? (; total=t_total) :
        (; total=t_total, leaf_lu=cache.build_time, leaf_lu_bytes=cache.bytes)
    measure_and_emit!("fgs", solver, setup,
        "p=$(winner.p);mac=$(winner.mac);leaf=$(winner.leaf);" *
        "inner=$(winner.inner);tol_abs=$tol_abs",
        "shared set = Stage 3 winner tau=$(winner.tau); tree/leaf-matrix " *
        "split not exposed by the FGS ctor; " * common_note;
        history_solve! = () -> (h = pnl.ConvergenceHistory(:fgs_maxabs);
                                pnl._solve_history!(rotor, solver, h); h))
    solver = nothing; GC.gc()
end

function run_fgmres_fgs!(winner)
    println("--- fgmres_fgs (shared set, fixed sweeps) ---")
    t_precond = @elapsed P = pnl.FGSPreconditioner(rotor;
        sweeps=winner.sweeps, inner_iterations=winner.inner, rlx=1.0,
        expansion_order=winner.p, multipole_acceptance=winner.mac,
        leaf_size=winner.leaf, shrink=true, recenter=false)
    t_ctor = @elapsed solver = pnl.KrylovSolver(rotor; method=:fgmres,
        krylov_kw..., preconditioner=P)
    cache = P.solver.fgs.leaf_lu_cache
    setup = (; total=t_precond + t_ctor, precond=t_precond,
             (cache === nothing ? (;) :
              (; leaf_lu=cache.build_time, leaf_lu_bytes=cache.bytes))...)
    measure_and_emit!("fgmres_fgs", solver, setup,
        "sweeps=$(winner.sweeps);inner=$(winner.inner);p=$(winner.p);" *
        "mac=$(winner.mac);leaf=$(winner.leaf);rtol=1e-6",
        "shared set = Stage 3 winner tau=$(winner.tau); " * common_note;
        niter_from=() -> solver.niter,
        history_solve! = () -> (pnl._solve!(rotor, solver); solver.history))
    P = nothing; solver = nothing; GC.gc()
end

function run_fgmres_fgs_nfcache!(winner; budget_gib=0.0)
    st = nfcache_setup(budget_gib)
    if st === nothing
        println("--- fgmres_fgs_nfcache @ $(budget_gib) GiB --- SKIPPED " *
                "(no tuned row for this budget)")
        return
    end
    backend_nf, ck = st
    if !ck.cached
        println("--- fgmres_fgs_nfcache @ $(budget_gib) GiB --- SKIPPED (no " *
                "cache fits this budget; the uncached point is fgmres_fgs)")
        return
    end
    println("--- fgmres_fgs_nfcache @ $(budget_gib) GiB total (shared set " *
            "precond, cached-knob apply) ---")
    t_precond = @elapsed P = pnl.FGSPreconditioner(rotor;
        sweeps=winner.sweeps, inner_iterations=winner.inner, rlx=1.0,
        expansion_order=winner.p, multipole_acceptance=winner.mac,
        leaf_size=winner.leaf, shrink=true, recenter=false)
    t_ctor = @elapsed solver = pnl.KrylovSolver(rotor; method=:fgmres,
        krylov_kw..., preconditioner=P, backend=backend_nf,
        cache_nearfield=true, persistent_plan=true,
        nearfield_cache_max_bytes=round(Int, budget_gib * 1024^3))
    t_warmup, nfc = prime_nfcache!(solver)
    cache = P.solver.fgs.leaf_lu_cache
    setup = (; total=t_precond + t_ctor + t_warmup, precond=t_precond,
             (cache === nothing ? (;) :
              (; leaf_lu=cache.build_time, leaf_lu_bytes=cache.bytes))...)
    measure_and_emit!("fgmres_fgs_nfcache", solver, setup,
        "sweeps=$(winner.sweeps);inner=$(winner.inner);p=$(winner.p);" *
        "mac=$(winner.mac);leaf=$(winner.leaf);rtol=1e-6;cache_nearfield;" *
        "persistent_plan",
        "shared set = Stage 3 winner tau=$(winner.tau) (precond knobs " *
        "unchanged; only the operator apply is cached); WARM cache built once " *
        "in setup, build reported separately as nfcache_build_time; " *
        "cold solution (reset_cold!), t excludes N^2 source assembly";
        niter_from=() -> solver.niter,
        history_solve! = () -> (pnl._solve!(rotor, solver); solver.history),
        apply_knobs=(ck.p, ck.mac, ck.leaf), mem_budget=budget_gib,
        nfcache_build=nfc.build_time, nfcache_bytes=nfc.bytes,
        nfcache_capped=ck.capped,
        mem_total=BODY_BYTES + solver_state_bytes(solver))
    P = nothing; solver = nothing; GC.gc()
end

"""
D4 additivity check (R1 only, ADDITIVITY_CHECK=1). Ryan 2026-08-24: "couldn't
we get both by merely performing the warm test? A cold solve would merely add
the cache build time, right?" This verifies that once, so every other rung can
rely on t_cold = t_build + t_warm instead of paying for a second timed pass.
Emits row_kind="additivity"; a mismatch means the cached rows cannot be read
as cold-comparable and the convention must be revisited.
"""
function run_additivity_check!(budget_gib)
    st = nfcache_setup(budget_gib)
    (st === nothing || !st[2].cached) && return
    backend_nf, ck = st
    println("--- additivity check @ $(budget_gib) GiB ---")
    # cold: persistent_plan OFF, so the plan + cache are rebuilt inside the
    # timed solve, exactly as the pre-W3 driver measured it
    solver_cold = pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...,
        backend=backend_nf, cache_nearfield=true,
        nearfield_cache_max_bytes=round(Int, budget_gib * 1024^3))
    reset_cold!(); pnl._solve!(rotor, solver_cold)          # JIT warm-up only
    t_cold = minimum(begin reset_cold!()
                           @elapsed pnl._solve!(rotor, solver_cold) end
                     for _ in 1:3)
    solver_cold = nothing; GC.gc()

    solver_warm = pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...,
        backend=backend_nf, cache_nearfield=true, persistent_plan=true,
        nearfield_cache_max_bytes=round(Int, budget_gib * 1024^3))
    _, nfc = prime_nfcache!(solver_warm)
    t_warm = minimum(begin reset_cold!()
                           @elapsed pnl._solve!(rotor, solver_warm) end
                     for _ in 1:3)
    predicted = nfc.build_time + t_warm
    rel = (t_cold - predicted) / t_cold
    println("  t_cold=$(round(t_cold; digits=3)) s  " *
            "t_build+t_warm=$(round(predicted; digits=3)) s  " *
            "(build $(round(nfc.build_time; digits=3)) + warm " *
            "$(round(t_warm; digits=3)))  rel discrepancy $(round(rel; digits=4))")
    abs(rel) < 0.10 || @warn "ADDITIVITY FAILED: t_cold and t_build+t_warm " *
        "differ by $(round(100*rel; digits=1))% — the warm rows cannot be " *
        "read as cold-comparable; revisit the convention before publishing"
    measure_and_emit!("krylov_gmres_nfcache", solver_warm,
        (; total=nfc.build_time), "additivity check",
        "D4: t_cold=$t_cold vs t_build+t_warm=$predicted, rel=$rel";
        niter_from=() -> solver_warm.niter,
        apply_knobs=(ck.p, ck.mac, ck.leaf), mem_budget=budget_gib,
        nfcache_build=nfc.build_time, nfcache_bytes=nfc.bytes,
        nfcache_capped=ck.capped, row_kind="additivity",
        mem_total=BODY_BYTES + solver_state_bytes(solver_warm))
    solver_warm = nothing; GC.gc()
end

# ---- dispatch --------------------------------------------------------------
# Every config below emits exactly one row_kind ("standard"), so a landed row
# means that config is finished and is skipped outright on a requeue.
_todo(name) = name in want && !has_landed(name)
# *_nfcache configs are keyed by (config, row_kind, budget): one row per point
# on the memory axis.
_todo_nf(name, b) = name in want && !has_landed(name, "standard", b)

const _NF = ("krylov_gmres_nfcache", "krylov_ilu_nfcache", "fgmres_fgs_nfcache")
for name in ("backslash_ldiv", "krylov_gmres", "krylov_jacobi", "krylov_ilu",
             "fgs", "fgmres_fgs")
    name in want && has_landed(name) && println("--- $name --- SKIPPED (already landed)")
end
_todo("backslash_ldiv") && run_backslash!()
for name in ("krylov_gmres", "krylov_jacobi", "krylov_ilu")
    _todo(name) && run_krylov!(name)
end
for b in nfcache_budgets, name in ("krylov_gmres_nfcache", "krylov_ilu_nfcache")
    if name in want && has_landed(name, "standard", b)
        println("--- $name @ $(b) GiB --- SKIPPED (already landed)")
    elseif _todo_nf(name, b)
        run_krylov!(name; budget_gib=b)
    end
end
if _todo("fgs") || _todo("fgmres_fgs") ||
        any(b -> _todo_nf("fgmres_fgs_nfcache", b), nfcache_budgets)
    winner, tol_abs = fgs_shared_set()
    println("shared FGS knob set (Stage 3 winner τ=$(winner.tau)): " *
            "p=$(winner.p) MAC=$(winner.mac) leaf=$(winner.leaf) " *
            "inner=$(winner.inner) sweeps=$(winner.sweeps)")
    _todo("fgs") && run_fgs!(winner, tol_abs)
    _todo("fgmres_fgs") && run_fgmres_fgs!(winner)
    for b in nfcache_budgets
        _todo_nf("fgmres_fgs_nfcache", b) &&
            run_fgmres_fgs_nfcache!(winner; budget_gib=b)
    end
end

# D4: verify t_cold = t_build + t_warm ONCE (R1 by default), then rely on it.
if get(ENV, "ADDITIVITY_CHECK", rung == "R1" ? "1" : "0") == "1"
    # smallest budget whose tuned winner actually uses a cache — there is
    # nothing to check on an uncached point
    b_add = findfirst(b -> (k = phase2_knobs(b); k !== nothing && k.cached),
                      sort(nfcache_budgets))
    if b_add !== nothing
        b = sort(nfcache_budgets)[b_add]
        has_landed("krylov_gmres_nfcache", "additivity", string(b)) ||
            run_additivity_check!(b)
    end
end

close(io)
println("\n$rung phase2 done. Rows appended to $csv_path")
