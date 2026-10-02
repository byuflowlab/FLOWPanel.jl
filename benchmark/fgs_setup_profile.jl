#=##############################################################################
BRAINSTORM 033 B-T1: FGS setup-cost attribution (profiling only — NO solves,
NO solver-behavior changes).

Decomposes FGSSolver construction (R4 champion knobs) into the actual phases
of FastMultipole.FastGaussSeidel, timed non-invasively by REPLAYING the
constructor's own call sequence in this harness (each phase is the same
internal call the constructor makes, in the same order, with the same
arguments — see FastMultipole/src/solve.jl `FastGaussSeidel`). The real
end-to-end `pnl.FGSSolver` constructor is also timed as a cross-check that
the phase sum accounts for the total.

Phases (constructor order):
  trees      : source Tree + target-tree replay + assert_shared_topology
  lists      : build_interaction_lists + the two stable double-sorts
  nonself    : nonself_influence_matrices (near-field non-self block assembly,
               per-source-body unit-strength probing) + old_influence_storage
  bookkeep   : add_self_interactions, index_by_source, map_by_leaf,
               map_by_branch, small vectors
  self       : self_influence_matrices (leaf self blocks, same probing)
  lu         : build_leaf_lu_cache (per-leaf Float64 lu!)
  dagplan    : build_dagteam_plan (DAG build, Lmat/Umat split repack with
               TM conversion, Float32 leaf LU refactorization for :f32full,
               scratch allocation)
  dag_lu_f32 : standalone re-timing of build_leaf_lu_cache_as(Float32, ·) —
               an inside-view of dagplan's LU component (NOT added to the sum)

Usage (B-T1 protocol; run once per (rung, threads) process):
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 THREADING_MODE=multi \
  EXPECT_JULIA_THREADS=4 BENCH_BLAS_THREADS=1 RUNG=R4 SKIP_B=1 \
  DECOMP_K=2 FULL_K=1 SETUP_THREADS=4 \
    julia --project=benchmark -t 4 benchmark/fgs_setup_profile.jl <outdir>

SETUP_THREADS (B-I2, 2026-10-01): 0 (default) = original serial probe replay;
n>0 = B-I1 threaded influence-matrix population (nonself/self) with n threads,
and the full-ctor cross-check passes threaded_setup=true. The dagplan phase is
additionally decomposed into dagplan_{edges,alloc,repack,prio,scratch,lu} rows
(informational; excluded from phase_sum).

Every pass of every phase is a CSV row (judge takes min per phase across
passes); pass 0 is the compile warmup and is written but flagged warmup=1.
=###############################################################################

import TOML

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

const FM = pnl.FastMultipole

outdir = isempty(ARGS) ? pwd() : ARGS[1]
mkpath(outdir)
open(joinpath(outdir, "banner_$(rung)_j$(banner.julia_threads).txt"), "w") do io
    println(io, banner.text)
end

# Champion R4 knobs. For rungs other than R4 the same knobs are transplanted
# (P8/MAC0.4/leaf100/f32full) — labeled as such; they are NOT that rung's
# tuned champion, but keep the phase decomposition knob-identical across rungs.
champ = TOML.parsefile(joinpath(@__DIR__, "retained_r4_champion.toml"))
champ["rung"] == rung || @warn "champion TOML is $(champ["rung"]); transplanting its knobs onto $rung"

P            = champ["P"]
MAC          = Float64(champ["MAC"])
leaf         = champ["leaf"]
precision    = Symbol(champ["dagteam_precision"])
shrink_      = true
recenter_    = false

decomp_k = parse(Int, get(ENV, "DECOMP_K", "2"))
full_k   = parse(Int, get(ENV, "FULL_K", "1"))

# B-I2 (2026-10-01): SETUP_THREADS=<n> replays the ctor with B-I1's threaded
# influence-matrix population (nonself/self probes); 0 = original serial
# replay. The full-ctor cross-check passes threaded_setup accordingly.
setup_threads_env = parse(Int, get(ENV, "SETUP_THREADS", "0"))

################################################################################
# Phase-decomposed replay of FastMultipole.FastGaussSeidel((rotor,); ...)
################################################################################

# one decomposed pass; returns Vector{(phase, t_s, bytes)}
function decomposed_pass()
    rows = Tuple{String,Float64,Int}[]
    target_systems = source_systems = (rotor,)
    leaf_size = FM.to_vector(leaf, length(source_systems))
    interaction_list_method = FM.Barba()
    derivatives_switches = FM.DerivativesSwitch(true, true, false, target_systems)
    TF = pnl.numtype(rotor)

    # --- trees ---
    t0 = time_ns(); a0 = Base.gc_bytes()
    switches = FM.DerivativesSwitch(true, true, true, source_systems)
    source_tree = FM.Tree(source_systems, FM.SourceTree(), switches;
        expansion_order=P, leaf_size, shrink=shrink_, recenter=recenter_,
        interaction_list_method)
    switches = FM.DerivativesSwitch(true, true, true, target_systems)
    target_tree = FM.Tree(source_tree, target_systems, switches;
        shrink=shrink_, recenter=recenter_)
    FM.assert_shared_topology(target_tree, source_tree)
    push!(rows, ("trees", (time_ns()-t0)/1e9, Int(Base.gc_bytes()-a0)))

    # --- interaction lists + canonical sorts ---
    t0 = time_ns(); a0 = Base.gc_bytes()
    farfield, nearfield, self_induced = true, true, false
    m2l_list, direct_list = FM.build_interaction_lists(target_tree.branches,
        source_tree.branches, leaf_size, MAC, farfield, nearfield,
        self_induced, interaction_list_method)
    m2l_list = FM.sort_by_target(FM.sort_by_source(m2l_list, source_tree.branches),
                                 target_tree.branches)
    direct_list = FM.sort_by_target(FM.sort_by_source(direct_list, source_tree.branches),
                                    target_tree.branches)
    push!(rows, ("lists", (time_ns()-t0)/1e9, Int(Base.gc_bytes()-a0)))

    # --- non-self influence matrices ---
    t0 = time_ns(); a0 = Base.gc_bytes()
    nonself_matrices, sorted_list = FM.nonself_influence_matrices(
        target_tree.buffers, source_tree.buffers, source_systems,
        target_tree, source_tree, direct_list, derivatives_switches;
        setup_threads=setup_threads_env)
    old_influence_storage = similar(nonself_matrices.rhs)
    push!(rows, ("nonself", (time_ns()-t0)/1e9, Int(Base.gc_bytes()-a0)))

    # --- index/bookkeeping ---
    t0 = time_ns(); a0 = Base.gc_bytes()
    full_direct_list = FM.add_self_interactions(direct_list, source_tree)
    index_map = FM.index_by_source(sorted_list, source_tree.leaf_index)
    strengths = zeros(TF, FM.get_n_bodies(source_systems))
    strengths_by_leaf = FM.map_by_leaf(source_tree)
    targets_by_branch = FM.map_by_branch(target_tree)
    push!(rows, ("bookkeep", (time_ns()-t0)/1e9, Int(Base.gc_bytes()-a0)))

    # --- self influence matrices ---
    t0 = time_ns(); a0 = Base.gc_bytes()
    self_matrices = FM.self_influence_matrices(target_tree.buffers,
        source_tree.buffers, source_systems, target_tree, source_tree,
        derivatives_switches; setup_threads=setup_threads_env)
    push!(rows, ("self", (time_ns()-t0)/1e9, Int(Base.gc_bytes()-a0)))

    # --- leaf LU cache (Float64) ---
    t0 = time_ns(); a0 = Base.gc_bytes()
    leaf_lu_cache = FM.build_leaf_lu_cache(self_matrices;
        setup_threads=setup_threads_env)
    push!(rows, ("lu", (time_ns()-t0)/1e9, Int(Base.gc_bytes()-a0)))

    # --- dagteam plan (split repack + f32 LU + scratch) ---
    # B-I2: setup_diagnostics decomposes the plan build into its internal
    # stages (edges/alloc/repack/prio/scratch/lu); emitted as dagplan_* rows,
    # informational only (NOT added to phase_sum — dagplan already covers them)
    dag_diag = Dict{Symbol,UInt64}()
    t0 = time_ns(); a0 = Base.gc_bytes()
    dagteam = FM.build_dagteam_plan(precision, nonself_matrices,
        sorted_list, index_map, source_tree, target_tree,
        strengths_by_leaf, targets_by_branch, self_matrices, leaf_lu_cache;
        nworkers=Threads.nthreads(), idle_policy=:backoff, coop=1,
        setup_threads=setup_threads_env, setup_diagnostics=dag_diag)
    # NOTE: coop=1 (solo/production, bit-identical) mirrors what the A-R2-
    # modified constructor passes; the kwarg requires the A-R2 signature of
    # build_dagteam_plan (uncommitted on flowpanel-20260817) — drop it to run
    # against pre-A-R2 FastMultipole.
    push!(rows, ("dagplan", (time_ns()-t0)/1e9, Int(Base.gc_bytes()-a0)))
    for key in (:dagplan_edges_ns, :dagplan_alloc_ns, :dagplan_repack_ns,
                :dagplan_prio_ns, :dagplan_scratch_ns, :dagplan_lu_ns)
        haskey(dag_diag, key) || continue
        push!(rows, (string(key)[1:end-3], dag_diag[key]/1e9, 0))
    end

    # --- inside view: the f32 LU component of dagplan (informational only) ---
    t0 = time_ns(); a0 = Base.gc_bytes()
    lus32 = FM.build_leaf_lu_cache_as(Float32, self_matrices;
        setup_threads=setup_threads_env)
    push!(rows, ("dag_lu_f32", (time_ns()-t0)/1e9, Int(Base.gc_bytes()-a0)))

    # keep references alive through timing
    return rows, (m2l_list, full_direct_list, old_influence_storage,
                  strengths, dagteam, lus32)
end

################################################################################
# Run
################################################################################

csv_path = joinpath(outdir,
    "fgs_setup_profile_$(rung)_j$(banner.julia_threads)_s$(setup_threads_env).csv")
open(csv_path, "w") do io
    println(io, "rung,n_panels,julia_threads,blas_threads,setup_threads,kind,pass,warmup,phase,t_s,gc_bytes")

    # decomposed passes (pass 0 = compile warmup, still recorded)
    for pass in 0:decomp_k
        GC.gc(); GC.gc()
        rows, _keep = decomposed_pass()
        tot = sum(r[2] for r in rows
                  if r[1] != "dag_lu_f32" && !startswith(r[1], "dagplan_"))
        for (phase, t, bytes) in rows
            println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
                "$(banner.blas_threads),$setup_threads_env,decomp,$pass,$(pass == 0 ? 1 : 0)," *
                "$phase,$t,$bytes")
        end
        println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
            "$(banner.blas_threads),$setup_threads_env,decomp,$pass,$(pass == 0 ? 1 : 0)," *
            "phase_sum,$tot,0")
        println("decomp pass $pass: phase_sum = $(round(tot, digits=2)) s")
        flush(io)
    end

    # full real constructor, end-to-end (warmup already paid above via the
    # same code paths; still run FULL_K+0 timed reps, first flagged)
    for pass in 1:full_k
        GC.gc(); GC.gc()
        t0 = time_ns(); a0 = Base.gc_bytes()
        solver = pnl.FGSSolver(rotor; expansion_order=P,
            multipole_acceptance=MAC, leaf_size=leaf,
            inner_iterations=champ["inner"], max_iterations=champ["max_iterations"],
            tolerance=champ["tolerance"], rlx=champ["rlx"], shrink=shrink_,
            recenter=recenter_, reverse_pass=false, cache_leaf_lu=true,
            sweep_order=:dagteam, dagteam_precision=precision,
            threaded_setup=setup_threads_env > 0,
            verbose=false, project_solution=false, solution_history_length=0)
        t = (time_ns()-t0)/1e9
        println(io, "$rung,$(rotor.ncells),$(banner.julia_threads)," *
            "$(banner.blas_threads),$setup_threads_env,full_ctor,$pass,0,full_ctor,$t," *
            "$(Int(Base.gc_bytes()-a0))")
        println("full ctor pass $pass: $(round(t, digits=2)) s")
        flush(io)
        solver = nothing
    end
end
println("CSV written to $csv_path")
