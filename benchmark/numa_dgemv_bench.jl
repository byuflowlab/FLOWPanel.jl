# NUMA-placement dgemv microbenchmark (021, 2026-09-18).
# Mimics the FGS nonself influence-cache sweep: many small Float64 blocks,
# one dgemv each per pass, streamed every pass. Pure Base + LinearAlgebra —
# no FLOWPanel/FastMultipole involvement.
#
# ARM (env):
#   a = 1 thread sweep, single-thread first-touch      (lex analogue)
#   b = 64-thread static sweep, single-thread first-touch (chunked v22 analogue)
#   c = same code as b; launcher wraps it in `numactl --interleave=0-3`
#   d = 64-thread static sweep, parallel chunk-affine first-touch
#       (each thread allocates+fills its own static partition)
#
# Page placement is captured (numastat -p / numa_maps) between fill and sweep.

using LinearAlgebra, Statistics

BLAS.set_num_threads(1)

const ARM   = ENV["ARM"]
const OUT   = ENV["OUTDIR"]
const NMAT  = parse(Int, get(ENV, "NMAT",  "1068"))
const MROWS = parse(Int, get(ENV, "MROWS", "650"))
const NCOLS = parse(Int, get(ENV, "NCOLS", "550"))
const NWARM = parse(Int, get(ENV, "NWARM", "5"))
const NPASS = parse(Int, get(ENV, "NPASS", "20"))

mkpath(OUT)
const threaded = ARM != "a"
const nt = threaded ? Threads.nthreads() : 1

blocks = Vector{Matrix{Float64}}(undef, NMAT)
xs = Vector{Vector{Float64}}(undef, NMAT)
ys = Vector{Vector{Float64}}(undef, NMAT)

@inline function fill_one!(blocks, xs, ys, k)
    A = Matrix{Float64}(undef, MROWS, NCOLS)
    fill!(A, 1.0 / k)
    blocks[k] = A
    xs[k] = fill(1.0, NCOLS)
    ys[k] = zeros(MROWS)
    return nothing
end

if ARM == "d"
    # chunk-affine first-touch: identical static partition to the sweep loop
    Threads.@threads :static for k in 1:NMAT
        fill_one!(blocks, xs, ys, k)
    end
else
    for k in 1:NMAT
        fill_one!(blocks, xs, ys, k)
    end
end

# show page placement (not just infer it), after fill / before sweep
try
    open(joinpath(OUT, "numastat_arm$(ARM).txt"), "w") do io
        run(pipeline(`numastat -p $(getpid())`; stdout=io, stderr=io))
    end
catch err
    @warn "numastat capture failed" err
end
try
    cp("/proc/self/numa_maps", joinpath(OUT, "numa_maps_arm$(ARM).txt"); force=true)
catch err
    @warn "numa_maps capture failed" err
end

function sweep!(ys, blocks, xs)
    if threaded
        Threads.@threads :static for k in eachindex(blocks)
            mul!(ys[k], blocks[k], xs[k])
        end
    else
        for k in eachindex(blocks)
            mul!(ys[k], blocks[k], xs[k])
        end
    end
    return nothing
end

for _ in 1:NWARM
    sweep!(ys, blocks, xs)
end

times = Float64[]
for _ in 1:NPASS
    t0 = time_ns()
    sweep!(ys, blocks, xs)
    push!(times, (time_ns() - t0) * 1e-9)
end

bytes = NMAT * (MROWS * NCOLS + MROWS + NCOLS) * 8
med = median(times)
gbps = bytes / med / 1e9

open(joinpath(OUT, "result_arm$(ARM).csv"), "w") do io
    println(io, "arm,nthreads,nmat,mrows,ncols,bytes_per_pass,nwarm,npass,median_s,min_s,max_s,gbps_median")
    println(io, join((ARM, nt, NMAT, MROWS, NCOLS, bytes, NWARM, NPASS,
                      med, minimum(times), maximum(times), gbps), ","))
end
open(joinpath(OUT, "passes_arm$(ARM).csv"), "w") do io
    println(io, "pass,seconds")
    for (i, t) in enumerate(times)
        println(io, "$i,$t")
    end
end
println("arm=$ARM nthreads=$nt bytes/pass=$bytes median=$(med) s -> $(round(gbps; digits=1)) GB/s")
