#=##############################################################################
BRAINSTORM 021: R1 FastGaussSeidel determinism attribution probe.

With no PROBE_WORKER setting this file launches, sequentially, the matrix
  Julia threads in {1,4} x BLAS threads in {1,4}.
No process uses more than four Julia or BLAS threads.  Each worker constructs
the frozen R1 case from phase1_case.jl and records compact, exact comparisons
under BRAINSTORM/021_rotor_hover_solver_benchmarks/results/fgs_determinism/.

Run:
  julia --project benchmark/fgs_determinism_probe.jl

Useful worker/debug overrides:
  GOTO_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    PROBE_WORKER=1 PROBE_JULIA_THREADS=4 PROBE_BLAS_THREADS=1 \
    EXPECT_JULIA_THREADS=4 THREADING_MODE=multi julia --project -t 4 \
    benchmark/fgs_determinism_probe.jl
=###############################################################################

if get(ENV, "PROBE_WORKER", "0") != "1"
    project = dirname(@__DIR__)
    script = @__FILE__
    for jt in (1, 4), bt in (1, 4)
        mode = jt == 1 ? "single" : "multi"
        cmd = `$(Base.julia_cmd()) --project=$project --threads=$jt $script`
        env = copy(ENV)
        env["PROBE_WORKER"] = "1"
        env["PROBE_JULIA_THREADS"] = string(jt)
        env["PROBE_BLAS_THREADS"] = string(bt)
        env["EXPECT_JULIA_THREADS"] = string(jt)
        env["THREADING_MODE"] = mode
        env["OPENBLAS_NUM_THREADS"] = string(bt)
        env["GOTO_NUM_THREADS"] = string(bt)
        env["OMP_NUM_THREADS"] = string(bt)
        println("\n=== FGS determinism worker: Julia $jt, BLAS $bt ===")
        run(setenv(cmd, env))
    end
    exit()
end

include(joinpath(@__DIR__, "common.jl"))
import LinearAlgebra

const PROBE_JT = parse(Int, get(ENV, "PROBE_JULIA_THREADS", "0"))
const PROBE_BT = parse(Int, get(ENV, "PROBE_BLAS_THREADS", "0"))
Threads.nthreads() == PROBE_JT || error("worker Julia-thread mismatch")
PROBE_JT in (1, 4) || error("PROBE_JULIA_THREADS must be 1 or 4")
PROBE_BT in (1, 4) || error("PROBE_BLAS_THREADS must be 1 or 4")
LinearAlgebra.BLAS.set_num_threads(PROBE_BT)
LinearAlgebra.BLAS.get_num_threads() == PROBE_BT || error("BLAS-thread mismatch")

# phase1_case owns the reproducible R1 geometry and frozen right-hand side.
ENV["RUNG"] = "R1"
include(joinpath(@__DIR__, "phase1_case.jl"))
# `phase1_case` calls the campaign banner, whose legacy single-mode contract
# pins BLAS to one.  Restore the requested matrix coordinate afterwards.
LinearAlgebra.BLAS.set_num_threads(PROBE_BT)
LinearAlgebra.BLAS.get_num_threads() == PROBE_BT || error("BLAS-thread mismatch after case setup")

const OUTDIR = joinpath(dirname(@__DIR__), "BRAINSTORM",
    "021_rotor_hover_solver_benchmarks", "results", "fgs_determinism")
mkpath(OUTDIR)
const CONFIG = get(ENV, "PROBE_LABEL", "j$(PROBE_JT)_b$(PROBE_BT)")
const N_SOLVES = parse(Int, get(ENV, "PROBE_SOLVE_REPS", "5"))
const N_CTORS = parse(Int, get(ENV, "PROBE_CTOR_REPS", "3"))
const N_OUTER = parse(Int, get(ENV, "PROBE_OUTER_ITERATIONS", "10"))
const N_TIMED = parse(Int, get(ENV, "PROBE_TIMED_REPS", "5"))

bits(x::AbstractArray{Float64}) = collect(reinterpret(UInt64, vec(copy(x))))
bits(x::Float64) = reinterpret(UInt64, x)

function firstdiff(a, b)
    length(a) == length(b) || return min(length(a), length(b)) + 1
    return findfirst(i -> !isequal(a[i], b[i]), eachindex(a))
end

function topology_signature(tree)
    out = Int64[]
    append!(out, Int64.(tree.leaf_index))
    for level in tree.levels_index
        append!(out, (Int64(first(level)), Int64(last(level))))
    end
    for v in tree.sort_index_list
        append!(out, Int64.(v))
    end
    for v in tree.inverse_sort_index_list
        append!(out, Int64.(v))
    end
    for b in tree.branches
        append!(out, Int64.(b.n_bodies))
        for r in b.bodies_index
            append!(out, (Int64(first(r)), Int64(last(r))))
        end
        append!(out, (b.n_branches, first(b.branch_index), last(b.branch_index),
                      b.i_parent, b.i_leaf))
    end
    return out
end

function geometry_signature(tree)
    out = Float64[]
    for b in tree.branches
        append!(out, b.center)
        push!(out, b.radius)
        append!(out, b.box)
        append!(out, (b.min_potential, b.min_gradient))
    end
    return bits(out)
end

function ctor_signature(solver)
    f = solver.fgs
    return (
        source_topology=topology_signature(f.source_tree),
        target_topology=topology_signature(f.target_tree),
        source_geometry=geometry_signature(f.source_tree),
        target_geometry=geometry_signature(f.target_tree),
        m2l=copy(f.m2l_list),
        direct=copy(f.direct_list),
        full_direct=copy(f.full_direct_list),
        index_map=copy(f.index_map),
        strengths_by_leaf=copy(f.strengths_by_leaf),
        targets_by_branch=copy(f.targets_by_branch),
        self_sizes=copy(f.self_matrices.sizes),
        self_matrix=bits(f.self_matrices.data),
        nonself_sizes=copy(f.nonself_matrices.sizes),
        nonself_matrix=bits(f.nonself_matrices.data),
    )
end

make_solver() = pnl.FGSSolver(rotor; expansion_order=6,
    multipole_acceptance=0.3, leaf_size=150, inner_iterations=5,
    max_iterations=N_OUTER, tolerance=-1.0, rlx=1.0, shrink=true,
    recenter=false, reverse_pass=false, verbose=false)

struct Trajectory
    residual::Vector{UInt64}
    strengths::Vector{Vector{UInt64}}
    source_buffers::Vector{Vector{UInt64}}
end

function cold_trajectory!(solver)
    residual = UInt64[]
    strengths = Vector{Vector{UInt64}}()
    source_buffers = Vector{Vector{UInt64}}()
    callback = function (_, r)
        push!(residual, bits(Float64(r)))
        push!(strengths, bits(solver.fgs.strengths))
        push!(source_buffers, bits(solver.fgs.source_tree.buffers[1]))
        return nothing
    end
    reset_cold!()
    pnl._solve!(rotor, solver; callback)
    push!(strengths, bits(solver.fgs.strengths))
    push!(source_buffers, bits(solver.fgs.source_tree.buffers[1]))
    return Trajectory(residual, strengths, source_buffers)
end

rows = NamedTuple[]
function record!(scope, replicate, component, pass, first_index, iteration=0; note="")
    push!(rows, (; config=CONFIG, julia_threads=PROBE_JT, blas_threads=PROBE_BT,
        scope, replicate, component, pass, first_index=something(first_index, 0),
        iteration, note))
end

println("--- constructor fingerprints ($N_CTORS fresh constructions) ---")
solvers = [make_solver() for _ in 1:N_CTORS]
signatures = ctor_signature.(solvers)

# Compare the generic O(n log n) lexicographic sort used by the fix with
# FastMultipole's stable O(n + branches) counting sorts.  Source-then-target
# produces the same (target, source) lexicographic order because both passes
# are stable.
f_order = solvers[1].fgs
raw_m2l, _ = FastMultipole.build_interaction_lists(
    f_order.target_tree.branches, f_order.source_tree.branches,
    f_order.source_tree.leaf_size, f_order.multipole_acceptance,
    true, true, false, f_order.interaction_list_method)
generic_order(list) = sort!(copy(list); by=ij -> (ij[1], ij[2]))
fm_order(list) = FastMultipole.sort_by_target(
    FastMultipole.sort_by_source(copy(list), f_order.source_tree.branches),
    f_order.target_tree.branches)
generic_order(raw_m2l) # compilation warmups
fm_order(raw_m2l)
generic_times = [@elapsed generic_order(raw_m2l) for _ in 1:5]
fm_times = [@elapsed fm_order(raw_m2l) for _ in 1:5]
generic_order(raw_m2l) == fm_order(raw_m2l) ||
    error("generic and FastMultipole canonical M2L orders differ")
open(joinpath(OUTDIR, "ordering_$(CONFIG).csv"), "w") do io
    println(io, "config,n_interactions,method,rep,t_seconds,minimum")
    for (method, samples) in (("generic_lexicographic", generic_times),
                              ("fm_stable_source_then_target", fm_times))
        for (rep, t) in enumerate(samples)
            println(io, join(_csv_cell.((CONFIG, length(raw_m2l), method,
                rep, t, minimum(samples))), ","))
        end
    end
end
println("ordering benchmark: generic=$(minimum(generic_times)) s, " *
        "FastMultipole two-pass=$(minimum(fm_times)) s, n=$(length(raw_m2l))")

for rep in 2:N_CTORS
    for name in propertynames(signatures[1])
        a = getproperty(signatures[1], name)
        candidate = getproperty(signatures[rep], name)
        d = firstdiff(a, candidate)
        record!("constructor", rep, String(name), d === nothing, d)
    end
end

println("--- repeated cold trajectories ($N_SOLVES solves, $N_OUTER fixed outer iterations) ---")
solver = solvers[1]
trajectories = [cold_trajectory!(solver) for _ in 1:N_SOLVES]
for rep in 2:N_SOLVES
    dres = firstdiff(trajectories[1].residual, trajectories[rep].residual)
    record!("solve", rep, "residual_sequence", dres === nothing, dres,
            something(dres, 0))
    n = min(length(trajectories[1].strengths), length(trajectories[rep].strengths))
    differing_iter = 0
    differing_index = nothing
    for i in 1:n
        d = firstdiff(trajectories[1].strengths[i], trajectories[rep].strengths[i])
        if d !== nothing
            differing_iter = i - 1 # callback 1 is the cold state; final is N_OUTER
            differing_index = d
            break
        end
    end
    if differing_index === nothing && length(trajectories[1].strengths) != length(trajectories[rep].strengths)
        differing_iter = n
        differing_index = 1
    end
    record!("solve", rep, "iterate_strengths", differing_index === nothing,
            differing_index, differing_iter)
    differing_iter = 0
    differing_index = nothing
    nbuf = min(length(trajectories[1].source_buffers),
               length(trajectories[rep].source_buffers))
    for i in 1:nbuf
        d = firstdiff(trajectories[1].source_buffers[i],
                      trajectories[rep].source_buffers[i])
        if d !== nothing
            differing_iter = i - 1
            differing_index = d
            break
        end
    end
    record!("solve", rep, "source_buffers", differing_index === nothing,
            differing_index, differing_iter)
end

# Independently replay the four far-field stages from fixed inputs.  Each
# downstream stage is restored from replicate 1's input, so an upstream bit
# difference cannot contaminate attribution.
if any(!r.pass for r in rows if r.scope == "solve")
    println("--- fixed-input far-field stage fingerprints ---")
    f = solver.fgs
    source_tree = f.source_tree
    target_tree = f.target_tree
    p = source_tree.expansion_order
    lamb = Val(FastMultipole.has_vector_potential((rotor,)))
    switches = FastMultipole.DerivativesSwitch(true, false, false, (rotor,))
    nthreads = Threads.nthreads()

    reset_cold!()
    FastMultipole.system_to_buffer!(source_tree.buffers, (rotor,),
                                    source_tree.sort_index_list)

    upward = Vector{Vector{UInt64}}()
    for _ in 1:5
        source_tree.expansions .= 0
        FastMultipole.upward_pass_multithread!(source_tree, (rotor,), p,
                                               lamb, nthreads)
        push!(upward, bits(source_tree.expansions))
    end
    fixed_upward = copy(source_tree.expansions)

    m2l = Vector{Vector{UInt64}}()
    for _ in 1:5
        source_tree.expansions .= fixed_upward
        target_tree.expansions .= 0
        FastMultipole.horizontal_pass_multithread!(target_tree, source_tree,
            f.m2l_list, lamb, p, nothing, f.interaction_list_method, nthreads)
        push!(m2l, bits(target_tree.expansions))
    end
    fixed_m2l = copy(target_tree.expansions)

    l2l = Vector{Vector{UInt64}}()
    for _ in 1:5
        target_tree.expansions .= fixed_m2l
        FastMultipole.downward_pass_multithread_1!(target_tree, p, lamb, nthreads)
        push!(l2l, bits(target_tree.expansions))
    end
    fixed_l2l = copy(target_tree.expansions)

    l2b = Vector{Vector{UInt64}}()
    for _ in 1:5
        target_tree.expansions .= fixed_l2l
        FastMultipole.reset!(target_tree.buffers)
        FastMultipole.downward_pass_multithread_2!(target_tree,
            target_tree.buffers, switches, p, lamb, nthreads)
        push!(l2b, bits(target_tree.buffers[1]))
    end

    for (component, samples) in (("upward", upward), ("m2l", m2l),
                                 ("l2l", l2l), ("l2b", l2b))
        for rep in 2:length(samples)
            d = firstdiff(samples[1], samples[rep])
            record!("farfield_stage", rep, component, d === nothing, d)
        end
    end
end

# A captured representative leaf isolates getrf/getrs from all FMM stages.
println("--- representative leaf factorization/solve ---")
leaf = findfirst(s -> s[1] >= 2 && s[1] == s[2], solver.fgs.self_matrices.sizes)
leaf === nothing && error("no representative square leaf matrix")
mat, _ = FastMultipole.get_matrix_vector(solver.fgs.self_matrices, leaf)
mat0 = copy(mat)
rhs0 = [sin(Float64(i)) for i in axes(mat0, 1)]
leaf_solutions = [bits(copy(mat0) \ rhs0) for _ in 1:10]
LinearAlgebra.BLAS.get_num_threads() == PROBE_BT ||
    error("BLAS thread count changed during representative leaf solves")
for rep in 2:length(leaf_solutions)
    d = firstdiff(leaf_solutions[1], leaf_solutions[rep])
    record!("leaf_lu", rep, "factor_and_solve", d === nothing, d; note="leaf=$leaf,size=$(size(mat0,1))")
end

# Performance contract: tolerance disabled, ten fixed outer iterations,
# one warmup and five timed genuine-cold solves.
println("--- fixed-iteration timing (1 warmup + $N_TIMED timed cold solves) ---")
reset_cold!(); pnl._solve!(rotor, solver)
times = Float64[]
for _ in 1:N_TIMED
    reset_cold!()
    t0 = time_ns()
    pnl._solve!(rotor, solver)
    push!(times, (time_ns() - t0) / 1e9)
end

csv = joinpath(OUTDIR, "probe_$(CONFIG).csv")
open(csv, "w") do io
    println(io, "config,julia_threads,blas_threads,scope,replicate,component,pass,first_index,iteration,note")
    for r in rows
        println(io, join(_csv_cell.((r.config, r.julia_threads, r.blas_threads,
            r.scope, r.replicate, r.component, r.pass, r.first_index,
            r.iteration, r.note)), ","))
    end
end

timing = joinpath(OUTDIR, "timing_$(CONFIG).csv")
open(timing, "w") do io
    println(io, "config,julia_threads,blas_threads,outer_iterations,rep,t_seconds,minimum,commit,fm_commit,julia_version,date")
    for (rep, t) in enumerate(times)
        println(io, join(_csv_cell.((CONFIG, PROBE_JT, PROBE_BT, N_OUTER,
            rep, t, minimum(times), banner.commit, banner.fm_commit,
            banner.julia_version, time_string())), ","))
    end
end

overall = all(r.pass for r in rows if r.scope in ("constructor", "solve", "leaf_lu"))
println("$CONFIG exact determinism: $(overall ? "PASS" : "FAIL"); " *
        "minimum fixed-$N_OUTER time $(minimum(times)) s")
println("wrote $csv and $timing")
