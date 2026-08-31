#=
Benchmark: sparse (working tree) vs dense (git HEAD a950790) cell list in
FLOWVPM.merge_particles!.

The dense implementation allocated two Int vectors of length
(extent/cell_size)^3 + 1 over the candidate bounding box, which OOM'd rotor
hover runs whenever a runaway particle stretched the box (job p018_L2_visc,
step ~1042, FLOWVPM_merging.jl:454). This script measures both methods on
identical data:

  Part A1: compact synthetic fields (the common case; regression check)
  Part A2: compact cluster + one outlier at growing distance (the pathology)
  Part B:  the real pre-OOM rotor states from data/p018_L2_visc_forensics/,
           including the reconstructed post-advection state the fatal merge saw

Dense runs predicted to allocate more than GUARD_BYTES are skipped and the
prediction 16*(n_cells+1) bytes is reported instead.

Run:  julia --project -t 1 benchmark/merge_sparse_vs_dense.jl
Outputs: benchmark/results/merge_sparse_vs_dense.csv, *_summary.md
=#

import FLOWVPM
using ReadVTK
using Printf
using Random

const RESULTS_DIR = joinpath(@__DIR__, "results")
const CSV_PATH = joinpath(RESULTS_DIR, "merge_sparse_vs_dense.csv")
const SUMMARY_PATH = joinpath(RESULTS_DIR, "merge_sparse_vs_dense_summary.md")
const FORENSICS_DIR = joinpath(dirname(@__DIR__), "data", "p018_L2_visc_forensics")

const GUARD_BYTES = 4 * 2^30

# ------------------------------------------------------------------ dense (old)
# Verbatim dense implementation from FLOWVPM.jl git HEAD a950790
# (src/FLOWVPM_merging.jl), namespaced so old and new run in one session on
# identical data. Shared helpers (_uf_union!, _merge_clusters_aggressive!) are
# imported from FLOWVPM and were not touched by the sparse rewrite.
module DenseMerge

import FLOWVPM: ParticleField, get_np, get_static, SIGMA_INDEX, GAMMA_INDEX,
                _uf_union!, _merge_clusters_aggressive!

function _build_cell_list!(
    sorted_indices::Vector{Int},
    offsets::Vector{Int},
    counts::Vector{Int},
    keys::Vector{Int},
    candidate_indices::Vector{Int},
    pfield::ParticleField,
    cell_size::Real,
    origin,
    Nx::Int,
    Ny::Int,
    Nz::Int,
)
    n_cells = Nx * Ny * Nz

    fill!(offsets, 0)

    for i in candidate_indices
        ix = clamp(floor(Int, (pfield.particles[1, i] - origin[1]) / cell_size), 0, Nx - 1)
        iy = clamp(floor(Int, (pfield.particles[2, i] - origin[2]) / cell_size), 0, Ny - 1)
        iz = clamp(floor(Int, (pfield.particles[3, i] - origin[3]) / cell_size), 0, Nz - 1)
        key = ix + iy * Nx + iz * Nx * Ny
        keys[i] = key
        offsets[key + 2] += 1
    end

    for i in 2:(n_cells + 1)
        offsets[i] += offsets[i - 1]
    end

    copyto!(counts, 1, offsets, 1, n_cells + 1)
    for i in candidate_indices
        key = keys[i] + 1
        counts[key] += 1
        sorted_indices[counts[key]] = i
    end

    return nothing
end

function merge_particles!(
    pfield::ParticleField;
    r_merge::Real=0.5,
    r_hash::Real=-1.0,
    sigma_relative::Bool=true,
    max_sigma_ratio::Real=2.0,
    skip_static::Bool=true,
    verbose::Bool=false,
    gamma_align_cos::Real=-1.0,
    on_representative::Union{Nothing,Function}=nothing,
)
    np = get_np(pfield)
    np <= 1 && return 0
    r_merge <= 0 && return 0
    max_sigma_ratio < 1 && return 0

    ws = pfield.merging_workspace
    candidate_indices = ws.candidate_indices
    empty!(candidate_indices)
    sizehint!(candidate_indices, np)

    xmin = typemax(eltype(pfield.particles))
    ymin = typemax(eltype(pfield.particles))
    zmin = typemax(eltype(pfield.particles))
    xmax = typemin(eltype(pfield.particles))
    ymax = typemin(eltype(pfield.particles))
    zmax = typemin(eltype(pfield.particles))
    sigma_sum = zero(eltype(pfield.particles))

    for i in 1:np
        if skip_static && get_static(pfield, i)
            continue
        end

        push!(candidate_indices, i)

        x = pfield.particles[1, i]
        y = pfield.particles[2, i]
        z = pfield.particles[3, i]
        sigma = pfield.particles[SIGMA_INDEX, i]

        xmin = min(xmin, x)
        ymin = min(ymin, y)
        zmin = min(zmin, z)
        xmax = max(xmax, x)
        ymax = max(ymax, y)
        zmax = max(zmax, z)
        sigma_sum += sigma
    end

    length(candidate_indices) <= 1 && return 0

    effective_r_hash = r_hash < 0.0 ? r_merge : r_hash
    mean_sigma = sigma_sum / length(candidate_indices)
    cell_size = sigma_relative ? effective_r_hash * mean_sigma : effective_r_hash
    if !(cell_size > 0)
        return 0
    end

    Nx = max(1, floor(Int, (xmax - xmin) / cell_size) + 1)
    Ny = max(1, floor(Int, (ymax - ymin) / cell_size) + 1)
    Nz = max(1, floor(Int, (zmax - zmin) / cell_size) + 1)
    n_cells = Nx * Ny * Nz

    sorted_indices = ws.sorted_indices
    resize!(sorted_indices, length(candidate_indices))

    offsets = ws.offsets
    resize!(offsets, n_cells + 1); fill!(offsets, 0)

    counts = ws.counts
    resize!(counts, n_cells + 1)

    keys = ws.keys
    resize!(keys, np)  # written before read for every candidate; no fill needed

    origin = (xmin, ymin, zmin)

    _build_cell_list!(sorted_indices, offsets, counts, keys, candidate_indices, pfield, cell_size, origin, Nx, Ny, Nz)

    parent = ws.parent
    rank = ws.rank
    resize!(parent, np)
    resize!(rank, np); fill!(rank, 0)
    @inbounds for i in 1:np
        parent[i] = i
    end

    for key in 0:(n_cells - 1)
        range_start = offsets[key + 1] + 1
        range_stop = offsets[key + 2]
        range_start > range_stop && continue

        for a in range_start:range_stop
            ia = sorted_indices[a]
            xi = pfield.particles[1, ia]
            yi = pfield.particles[2, ia]
            zi = pfield.particles[3, ia]
            sigma_i = pfield.particles[SIGMA_INDEX, ia]

            for b in (a + 1):range_stop
                ib = sorted_indices[b]
                sigma_j = pfield.particles[SIGMA_INDEX, ib]
                sigma_min = min(sigma_i, sigma_j)
                sigma_max = max(sigma_i, sigma_j)
                sigma_min <= 0 && continue
                sigma_max / sigma_min > max_sigma_ratio && continue

                if gamma_align_cos > -1.0
                    gax = pfield.particles[GAMMA_INDEX.start,     ia]
                    gay = pfield.particles[GAMMA_INDEX.start + 1, ia]
                    gaz = pfield.particles[GAMMA_INDEX.start + 2, ia]
                    gbx = pfield.particles[GAMMA_INDEX.start,     ib]
                    gby = pfield.particles[GAMMA_INDEX.start + 1, ib]
                    gbz = pfield.particles[GAMMA_INDEX.start + 2, ib]
                    ma2 = gax*gax + gay*gay + gaz*gaz
                    mb2 = gbx*gbx + gby*gby + gbz*gbz
                    if ma2 > 0 && mb2 > 0
                        cosθ = (gax*gbx + gay*gby + gaz*gbz) / sqrt(ma2 * mb2)
                        cosθ < gamma_align_cos && continue
                    end
                end

                dx = pfield.particles[1, ib] - xi
                dy = pfield.particles[2, ib] - yi
                dz = pfield.particles[3, ib] - zi
                dist2 = dx * dx + dy * dy + dz * dz
                r_pair = sigma_relative ? r_merge * sigma_min : r_merge
                dist2 < r_pair * r_pair && _uf_union!(parent, rank, ia, ib)
            end
        end
    end

    n_removed = _merge_clusters_aggressive!(pfield, ws; on_representative=on_representative)

    if verbose && n_removed > 0
        println("Merged $(length(candidate_indices)) candidate particles into $(length(candidate_indices) - n_removed) particles")
    end

    return n_removed
end

end # module DenseMerge

# --------------------------------------------------------------------- harness

"Prototype particle data from which fresh fields are built per repetition."
struct Proto
    X::Matrix{Float64}        # 3 x N
    Gamma::Matrix{Float64}    # 3 x N
    sigma::Vector{Float64}
    vol::Vector{Float64}
    circulation::Vector{Float64}
end

function build_field(proto::Proto)
    n = size(proto.X, 2)
    pfield = FLOWVPM.ParticleField(n)
    for i in 1:n
        FLOWVPM.add_particle(pfield,
            view(proto.X, :, i), view(proto.Gamma, :, i), proto.sigma[i];
            vol=proto.vol[i], circulation=proto.circulation[i])
    end
    return pfield
end

"Dense-grid cell count and predicted allocation for the two n_cells-sized buffers."
function dense_prediction(proto::Proto, cell_size)
    xmin = minimum(view(proto.X, 1, :)); xmax = maximum(view(proto.X, 1, :))
    ymin = minimum(view(proto.X, 2, :)); ymax = maximum(view(proto.X, 2, :))
    zmin = minimum(view(proto.X, 3, :)); zmax = maximum(view(proto.X, 3, :))
    Nx = max(1, floor(Int, (xmax - xmin) / cell_size) + 1)
    Ny = max(1, floor(Int, (ymax - ymin) / cell_size) + 1)
    Nz = max(1, floor(Int, (zmax - zmin) / cell_size) + 1)
    n_cells = big(Nx) * big(Ny) * big(Nz)
    return n_cells, 16 * (n_cells + 1)
end

"""
Best-of-`reps` measurement of `merge_fn` on fresh copies of `proto`.
Returns (time_s, bytes, n_removed, np_after). Each repetition rebuilds the
field, so the merging workspace is cold: reported bytes include the O(N)
workspace growth both methods share, plus whatever each method's cell list
costs on top.
"""
function measure(merge_fn, proto::Proto; reps=5, merge_kwargs...)
    best_t = Inf
    bytes = 0
    n_removed = -1
    np_after = -1
    for _ in 1:reps
        pfield = build_field(proto)
        GC.gc()
        stats = @timed merge_fn(pfield; merge_kwargs...)
        if stats.time < best_t
            best_t = stats.time
        end
        bytes = stats.bytes  # same every rep (cold workspace each time)
        n_removed = stats.value
        np_after = FLOWVPM.get_np(pfield)
    end
    return best_t, bytes, n_removed, np_after
end

const ROWS = NamedTuple[]

function record!(case, method, n, extent, n_cells, time_s, bytes, n_removed, np_after; note="")
    push!(ROWS, (; case, method, n, extent, n_cells, time_s, bytes, n_removed, np_after, note))
    tstr = time_s === missing ? "      --" : @sprintf("%8.4f", time_s)
    bstr = bytes === missing ? "        --" : _fmt_bytes(bytes)
    @printf("  %-28s %-8s N=%-8d ext=%9.3g cells=%-12.3g t=%s alloc=%10s removed=%s %s\n",
            case, method, n, extent, Float64(n_cells), tstr, bstr,
            n_removed === missing ? "--" : string(n_removed), note)
end

function _fmt_bytes(b)
    b < 1024^2 && return @sprintf("%.1f KiB", b / 1024)
    b < 1024^3 && return @sprintf("%.1f MiB", b / 1024^2)
    b < 1024^4 && return @sprintf("%.1f GiB", b / 1024^3)
    return @sprintf("%.1f TiB", Float64(b) / 1024^4)
end

"Run both methods on `proto`, guarding the dense one, and cross-check outcomes."
function bench_case!(case, proto::Proto; cell_size, reps=5, merge_kwargs...)
    n = size(proto.X, 2)
    n_cells, pred_bytes = dense_prediction(proto, cell_size)
    xspan = maximum(proto.X) - minimum(proto.X)

    t_new, b_new, rem_new, np_new =
        measure(FLOWVPM.merge_particles!, proto; reps, merge_kwargs...)
    record!(case, "sparse", n, xspan, n_cells, t_new, b_new, rem_new, np_new)

    if pred_bytes <= GUARD_BYTES
        t_old, b_old, rem_old, np_old =
            measure(DenseMerge.merge_particles!, proto; reps, merge_kwargs...)
        record!(case, "dense", n, xspan, n_cells, t_old, b_old, rem_old, np_old)
        if rem_old != rem_new || np_old != np_new
            @warn "METHOD MISMATCH" case rem_old rem_new np_old np_new
        else
            println("    outcomes match (removed=$(rem_new), np=$(np_new))")
        end
    else
        record!(case, "dense", n, xspan, n_cells, missing, pred_bytes, missing, missing;
                note="SKIPPED (predicted alloc exceeds guard)")
    end
    return nothing
end

# ------------------------------------------------------- Part A1: compact field

function compact_proto(n; cell_size, rng)
    extent = cbrt(n) * cell_size            # ~1 particle per occupied cell
    X = rand(rng, 3, n) .* extent
    Gamma = randn(rng, 3, n) .* 0.1
    sigma = fill(cell_size, n)
    return Proto(X, Gamma, sigma, ones(n), ones(n))
end

function part_a1!()
    println("\n== Part A1: compact fields (common case) ==")
    rng = Xoshiro(1234)
    h = 1.0
    for n in (10^3, 10^4, 10^5)
        proto = compact_proto(n; cell_size=h, rng)
        bench_case!("compact_N$(n)", proto; cell_size=h,
                    r_merge=0.5h, r_hash=h, sigma_relative=false, reps=5)
    end
end

# ------------------------------------------------------ Part A2: spread scaling

function part_a2!()
    println("\n== Part A2: compact cluster + one outlier (pathological case) ==")
    rng = Xoshiro(5678)
    h = 1.0
    n = 10^4
    base = compact_proto(n; cell_size=h, rng)
    for d_over_h in (10^2, 10^3, 10^4, 10^5, 10^6)
        X = copy(base.X)
        X[:, end] .= d_over_h * h
        proto = Proto(X, base.Gamma, base.sigma, base.vol, base.circulation)
        bench_case!("outlier_d$(d_over_h)h", proto; cell_size=h,
                    r_merge=0.5h, r_hash=h, sigma_relative=false, reps=3)
    end
end

# --------------------------------------------- Part B: rotor hover pre-OOM steps

"Load a forensics wake-particle VTP into a Proto (+ velocity for advection)."
function load_corpse(step)
    path = joinpath(FORENSICS_DIR, "p018_L2_visc_wake1_particles.$(step).vtp")
    vtk = VTKFile(path)
    X = get_points(vtk)
    pd = get_point_data(vtk)
    Gamma = get_data(pd["gamma"])
    sigma = vec(get_data(pd["sigma"]))
    vol = vec(get_data(pd["vol"]))
    circulation = vec(get_data(pd["circulation"]))
    velocity = get_data(pd["velocity"])
    return Proto(Matrix{Float64}(X), Matrix{Float64}(Gamma), Vector{Float64}(sigma),
                 Vector{Float64}(vol), Vector{Float64}(circulation)),
           Matrix{Float64}(velocity)
end

function part_b!()
    println("\n== Part B: rotor hover pre-OOM steps (job p018_L2_visc) ==")
    # Exact settings of the failing run: L2 rung, R = 0.119 m,
    # r_merge = 0.0027R, r_hash = 0.02R, absolute cell size.
    R = 0.119
    r_merge = 0.0027 * R
    r_hash = 0.02 * R
    dt = 3.086e-4
    kwargs = (; r_merge, r_hash, sigma_relative=false, max_sigma_ratio=2.0)

    local proto1041, vel1041
    for step in (1032, 1036, 1040, 1041)
        proto, vel = load_corpse(step)
        bench_case!("rotor_step$(step)", proto; cell_size=r_hash, reps=3, kwargs...)
        if step == 1041
            proto1041, vel1041 = proto, vel
        end
    end

    # The fatal state: step 1041 advected one Euler step by its own stored
    # velocity — the post-advection positions the step-1042 merge actually saw.
    Xadv = proto1041.X .+ dt .* vel1041
    fatal = Proto(Xadv, proto1041.Gamma, proto1041.sigma,
                  proto1041.vol, proto1041.circulation)
    bench_case!("rotor_step1042_fatal", fatal; cell_size=r_hash, reps=3, kwargs...)
end

# -------------------------------------------------------------------- reporting

function write_outputs()
    mkpath(RESULTS_DIR)
    open(CSV_PATH, "w") do io
        println(io, "case,method,n_particles,extent,n_cells_dense,time_s,bytes,n_removed,np_after,note")
        for r in ROWS
            println(io, join([r.case, r.method, r.n, r.extent, r.n_cells,
                              r.time_s === missing ? "" : r.time_s,
                              r.bytes === missing ? "" : r.bytes,
                              r.n_removed === missing ? "" : r.n_removed,
                              r.np_after === missing ? "" : r.np_after,
                              r.note], ","))
        end
    end

    open(SUMMARY_PATH, "w") do io
        println(io, "# merge_particles! sparse (working tree) vs dense (HEAD a950790) cell list\n")
        println(io, "Generated by `benchmark/merge_sparse_vs_dense.jl` on $(gethostname()), ",
                    "Julia $(VERSION), single thread. Time = best of reps on fresh fields ",
                    "(cold workspace); bytes = total allocation of one cold call. ",
                    "Dense rows marked SKIPPED report the predicted 16(n_cells+1)-byte ",
                    "allocation instead of running.\n")
        println(io, "| case | method | N | extent | dense n_cells | time [s] | alloc | removed |")
        println(io, "|---|---|---|---|---|---|---|---|")
        for r in ROWS
            t = r.time_s === missing ? "—" : @sprintf("%.4f", r.time_s)
            b = r.bytes === missing ? "—" : _fmt_bytes(r.bytes)
            b = r.note == "" ? b : b * " (predicted; skipped)"
            rem = r.n_removed === missing ? "—" : string(r.n_removed)
            println(io, "| $(r.case) | $(r.method) | $(r.n) | $(@sprintf("%.3g", r.extent)) | ",
                        "$(@sprintf("%.3g", Float64(r.n_cells))) | $t | $b | $rem |")
        end
    end
    println("\nWrote $(CSV_PATH)\nWrote $(SUMMARY_PATH)")
end

# ------------------------------------------------------------------------- main

function main()
    println("Warmup (compile both paths)...")
    rng = Xoshiro(1)
    warm = compact_proto(200; cell_size=1.0, rng)
    measure(FLOWVPM.merge_particles!, warm; reps=1, r_merge=0.5, r_hash=1.0, sigma_relative=false)
    measure(DenseMerge.merge_particles!, warm; reps=1, r_merge=0.5, r_hash=1.0, sigma_relative=false)

    part_a1!()
    part_a2!()
    part_b!()
    write_outputs()
end

main()
