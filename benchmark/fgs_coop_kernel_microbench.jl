#=##############################################################################
BRAINSTORM 033 A-R1: standalone cooperative-leaf-product kernel microbenchmarks.

Times the leaf lower-product y = L*x (the dagteam_pull! GEMV) executed
sequentially vs cooperatively by a persistent spin-barrier team of w=2/4
workers, in two layouts:

  - row:    disjoint output row ranges, one BLAS GEMV per worker on a
            row-block view (bitwise vs sequential is RECORDED, not assumed:
            BLAS may pick different blocking for a row-block view);
  - col:    contiguous column blocks into private partial vectors, then a
            fixed ascending-order reduction by the timing worker (A-T1:
            exact-arithmetic equivalent, deterministic, NOT bitwise).

Leaf shapes are real (nof, ptot) pairs read from the R4 DAG graph CSV
(BRAINSTORM/033_atheory_20260926/fgs_dag_L_graph_R4.csv), selected at nof
quantiles plus the mean-54 leaf and the 1450 tail leaf. A zero-work "null"
shape isolates the pure dispatch+barrier round-trip. Per-shape LU ldiv!
(the leaf self-solve) is timed for the critical-path model.

Matrices are synthesized dense Float32 (szTM=4 szTS=4, matching the R4
f32full champion) at the real shapes; matrix CONTENT does not affect BLAS
GEMV cost, only shape does, so synthesis is admissible for kernel timing
(structure evidence stays with the 033_atheory CSVs).

Workers are persistent tasks reusing a sense-reversing spin barrier — no
per-product task spawning (the dagedge anti-pattern). The timing thread
participates as worker 1, matching the intended executor integration.

Protocol (021 rulings): hand timing via time_ns(), min + median over reps
after warmup; thread/BLAS assertion + banner via benchmark/common.jl.
Judge from the CSV, never stdout.

Run (local, ≤4 threads, BLAS pinned single). NOTE: on hosts with OpenMP-backed
OpenBLAS a runtime set_num_threads(1) is reset to the core count on first real
work (common.jl's probe then hard-errors), so the env pinning below is REQUIRED:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  THREADING_MODE=multi EXPECT_JULIA_THREADS=4 BENCH_BLAS_THREADS=1 \
    julia --project=benchmark -t 4 benchmark/fgs_coop_kernel_microbench.jl [outdir]
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))

using LinearAlgebra: mul!, ldiv!, lu!
import Random

const banner = assert_and_banner()

const outdir = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "BRAINSTORM", "033_ar1_20260926")
mkpath(outdir)
open(joinpath(outdir, "banner.txt"), "w") do io; write(io, banner.text); end

const TM = Float32              # matrix storage (f32full champion)
const TS = Float32              # sweep-state type

################################################################################
# Shape selection from the real R4 leaf graph
################################################################################

const graph_csv = get(ENV, "GRAPH_CSV",
    joinpath(@__DIR__, "..", "BRAINSTORM", "033_atheory_20260926",
             "fgs_dag_L_graph_R4.csv"))

function load_shapes(path)
    rows = Tuple{Int,Int}[]     # (nof, ptot)
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        f = split(line, ',')
        push!(rows, (parse(Int, f[2]), parse(Int, f[3])))
    end
    rows
end

"""Leaf with nof nearest the target (ties: larger ptot), ptot > 0."""
function nearest_leaf(rows, target_nof)
    best = (typemax(Int), -1, 0, 0)
    for (nof, ptot) in rows
        ptot == 0 && continue
        d = abs(nof - target_nof)
        if d < best[1] || (d == best[1] && ptot > best[3])
            best = (d, nof, ptot, 0)
        end
    end
    (best[2], best[3])
end

const leafrows = load_shapes(graph_csv)
const nofs = sort([n for (n, p) in leafrows])
quant(q) = nofs[clamp(round(Int, q * length(nofs)), 1, length(nofs))]

shape_list = Tuple{String,Int,Int}[]                    # (label, nof, ptot)
push!(shape_list, ("null", 0, 0))                       # pure barrier cost
for (label, target) in [("p25", quant(0.25)), ("p50", quant(0.50)),
                        ("mean54", 54), ("p75", quant(0.75)),
                        ("p90", quant(0.90)), ("p99", quant(0.99)),
                        ("max", maximum(nofs))]
    nof, ptot = nearest_leaf(leafrows, target)
    push!(shape_list, (label, nof, ptot))
end

println("shapes: ", shape_list); flush(stdout)

################################################################################
# Persistent spin-barrier team
################################################################################

"""
Persistent team: w-1 spawned worker tasks plus the timing thread (worker 1).
Sense-reversing epoch protocol: main bumps `go`, workers spin on it, run the
current kernel closure for their worker id, bump `done`; main runs its own
slice then spins until all done. Workers never yield (pure spin), matching a
spin-barrier executor's elapsed path.
"""
mutable struct Team
    w::Int
    go::Threads.Atomic{Int}
    done::Threads.Atomic{Int}
    quit::Threads.Atomic{Int}
    kernel::Base.RefValue{Any}          # (worker_id) -> nothing
    tasks::Vector{Task}
end

function Team(w::Int)
    t = Team(w, Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
             Threads.Atomic{Int}(0), Ref{Any}(nothing), Task[])
    started = Threads.Atomic{Int}(0)
    for wid in 2:w
        task = Threads.@spawn begin
            epoch = 0
            Threads.atomic_add!(started, 1)
            spins = 0
            while true
                # Spin, but stay GC-cooperative and yield rarely so a task
                # co-scheduled with a spinning peer can still migrate; once
                # each worker owns a thread the yield branch never fires on
                # the hot path.
                while t.go[] == epoch
                    t.quit[] == 1 && return
                    ccall(:jl_cpu_pause, Cvoid, ())
                    GC.safepoint()
                    spins += 1
                    if spins >= 10_000
                        spins = 0
                        yield()
                    end
                end
                epoch = t.go[]
                (t.kernel[])(wid)
                Threads.atomic_add!(t.done, 1)
                spins = 0
            end
        end
        push!(t.tasks, task)
    end
    while started[] < w - 1
        yield()
    end
    t
end

"""One cooperative product: release team, run worker 1's slice, spin-join.
Returns elapsed ns for the full team round-trip (the A-T3 elapsed path)."""
@inline function team_round!(t::Team, post::F) where F
    t0 = time_ns()
    Threads.atomic_add!(t.go, 1)
    (t.kernel[])(1)
    if t.w > 1
        target = (t.go[] - 0) * (t.w - 1)   # cumulative dones expected
        spins = 0
        while t.done[] < target
            ccall(:jl_cpu_pause, Cvoid, ())
            GC.safepoint()
            spins += 1
            if spins >= 10_000
                spins = 0
                yield()
            end
        end
    end
    post()                                   # reduction (col layout) or no-op
    time_ns() - t0
end

function shutdown!(t::Team)
    t.quit[] = 1
    foreach(wait, t.tasks)
end

################################################################################
# Kernels
################################################################################

"""Contiguous partition of 1:n into w ranges (empty ranges allowed)."""
function partition(n, w)
    ranges = UnitRange{Int}[]
    base, rem = divrem(n, w)
    lo = 1
    for k in 1:w
        len = base + (k <= rem ? 1 : 0)
        push!(ranges, lo:lo+len-1)
        lo += len
    end
    ranges
end

struct LeafCase
    A::Matrix{TM}
    x::Vector{TS}
    y::Vector{TS}                        # shared output (row layout / seq)
    ybufs::Vector{Vector{TS}}            # private partials (col layout)
    rowparts::Dict{Int,Vector{UnitRange{Int}}}
    colparts::Dict{Int,Vector{UnitRange{Int}}}
end

function LeafCase(nof, ptot, wmax)
    Random.seed!(20260926 + nof)
    A = rand(TM, nof, ptot) .- TM(0.5)
    x = rand(TS, ptot) .- TS(0.5)
    y = zeros(TS, nof)
    ybufs = [zeros(TS, nof) for _ in 1:wmax]
    rowparts = Dict(w => partition(nof, w) for w in (1, 2, 4))
    colparts = Dict(w => partition(ptot, w) for w in (1, 2, 4))
    LeafCase(A, x, y, ybufs, rowparts, colparts)
end

seq_product!(c::LeafCase) = mul!(c.y, c.A, c.x)

"""Row-partitioned slice for worker wid: disjoint output rows, full x."""
@inline function row_slice!(c::LeafCase, w, wid)
    r = c.rowparts[w][wid]
    isempty(r) && return nothing
    mul!(view(c.y, r), view(c.A, r, :), c.x)
    nothing
end

"""Column-partitioned slice: private partial over a contiguous column block."""
@inline function col_slice!(c::LeafCase, w, wid)
    cr = c.colparts[w][wid]
    yb = c.ybufs[wid]
    if isempty(cr)
        fill!(yb, zero(TS))
    else
        mul!(yb, view(c.A, :, cr), view(c.x, cr))
    end
    nothing
end

"""Fixed ascending-order reduction (A-T1 determinism requirement)."""
@inline function col_reduce!(c::LeafCase, w)
    y = c.y
    copyto!(y, c.ybufs[1])
    @inbounds for k in 2:w
        yb = c.ybufs[k]
        @simd for i in eachindex(y)
            y[i] += yb[i]
        end
    end
    nothing
end

################################################################################
# Timing
################################################################################

"""min/median/mean over per-rep elapsed ns (after warmup), as µs."""
function time_reps(f::F, nreps) where F
    ts = Vector{Float64}(undef, nreps)
    for i in 1:min(50, nreps)                 # warmup
        f()
    end
    for i in 1:nreps
        ts[i] = Float64(f())
    end
    sort!(ts)
    (min=ts[1]/1e3, med=ts[cld(nreps,2)]/1e3, mean=sum(ts)/nreps/1e3)
end

nreps_for(nof, ptot) = clamp(round(Int, 2e8 / max(nof*ptot, 1000)), 200, 20_000)

################################################################################
# Run
################################################################################

results = NamedTuple[]
noop() = nothing

for (label, nof, ptot) in shape_list
    isnull = label == "null"
    c = isnull ? nothing : LeafCase(nof, ptot, 4)
    nreps = isnull ? 20_000 : nreps_for(nof, ptot)

    # Sequential baseline (no team machinery at all)
    if !isnull
        stats = time_reps(nreps) do
            t0 = time_ns(); seq_product!(c); time_ns() - t0
        end
        push!(results, (; label, nof, ptot, layout="seq", w=1, nreps,
                        stats.min, stats.med, stats.mean, bitwise=1))
        yref = copy(c.y)

        # LU self-solve at this leaf size (critical-path term)
        S = rand(TM, nof, nof) .+ TM(nof) .* Matrix{TM}(LinearAlgebra.I, nof, nof)
        F = lu!(S)
        rhs = rand(TS, nof)
        lustats = time_reps(nreps) do
            t0 = time_ns(); ldiv!(F, rhs); time_ns() - t0
        end
        push!(results, (; label, nof, ptot, layout="lu_ldiv", w=1, nreps,
                        lustats.min, lustats.med, lustats.mean, bitwise=1))

        for w in (1, 2, 4), layout in ("row", "col")
            team = Team(w)
            if layout == "row"
                team.kernel[] = wid -> row_slice!(c, w, wid)
                post = noop
            else
                team.kernel[] = wid -> col_slice!(c, w, wid)
                post = () -> col_reduce!(c, w)
            end
            # Correctness (A-T1 contract checks) before timing. NOTE: BLAS
            # gemv picks different internal blocking for a row-block view
            # than for the full matrix, so even ROW partitioning is not
            # bitwise vs sequential in general — record bitwise as DATA and
            # enforce only exact-arithmetic tolerance + determinism.
            fill!(c.y, TS(NaN))
            team_round!(team, post)
            bitwise = (c.y == yref)
            maximum(abs.(c.y .- yref)) <= 1e-4 * max(1, maximum(abs.(yref))) ||
                error("$layout layout accuracy failure ($label w=$w)")
            y1 = copy(c.y)
            team_round!(team, post)
            c.y == y1 || error("A-T1 VIOLATION: $layout layout not " *
                               "deterministic ($label w=$w)")
            stats = time_reps(nreps) do
                team_round!(team, post)
            end
            shutdown!(team)
            push!(results, (; label, nof, ptot, layout, w, nreps,
                            stats.min, stats.med, stats.mean,
                            bitwise=Int(bitwise)))
        end
    else
        # Null shape: pure dispatch + barrier round-trip
        for w in (1, 2, 4)
            team = Team(w)
            team.kernel[] = wid -> nothing
            stats = time_reps(nreps) do
                team_round!(team, noop)
            end
            shutdown!(team)
            push!(results, (; label, nof, ptot, layout="barrier", w, nreps,
                            stats.min, stats.med, stats.mean, bitwise=1))
        end
    end
    println("done: $label ($nof x $ptot)"); flush(stdout)
end

################################################################################
# CSV
################################################################################

csvpath = joinpath(outdir, "fgs_coop_kernel_microbench.csv")
open(csvpath, "w") do io
    println(io, "shape,nof,ptot,layout,w,nreps,t_min_us,t_med_us,t_mean_us,bitwise," *
                "julia_threads,blas_threads,commit,fm_commit")
    for r in results
        @printf(io, "%s,%d,%d,%s,%d,%d,%.3f,%.3f,%.3f,%d,%d,%d,%s,%s\n",
                r.label, r.nof, r.ptot, r.layout, r.w, r.nreps,
                r.min, r.med, r.mean, r.bitwise,
                banner.julia_threads, banner.blas_threads,
                banner.commit, banner.fm_commit)
    end
end
println("wrote $csvpath  (", length(results), " rows)")
