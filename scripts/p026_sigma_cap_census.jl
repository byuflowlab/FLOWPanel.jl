# BRAINSTORM 026 Phase 0 (W6): census of sigma-cap crossings by attribution.
#
# Offline pass over an existing run's particle VTP series (no simulation).
# For each adjacent step pair, particles are identity-matched by nearest
# neighbor in position (grid hash; particles move ~U*dt per step, far less
# than the inter-particle spacing in a resolved wake). Each matched particle
# crossing sigma_cap between the two steps is classified by comparing its
# realized Δσ² against the exact per-step viscous budget 2ν·dt (uniform for
# the Euler CoreSpreading scheme, FLOWVPM_viscous.jl): a crossing whose Δσ²
# fits inside (1+tol)·2ν·dt is viscous-dominated; anything larger is
# rVPM-compression-dominated. This gates Mechanism A (isotropic viscous
# split) per the design doc §10 — see
# BRAINSTORM/026_sigma_growth_particle_splitting/particle_splitting_design.md.
#
# Usage:
#   julia --project scripts/p026_sigma_cap_census.jl <particles_dir> \
#         --sigma-cap=<m> [--nu=<m²/s>] [--dt=<s>] [--tol=0.5] \
#         [--match-frac=0.5] [--first=<idx>] [--last=<idx>]
#
#   <particles_dir> holds <name>_particles.<idx>.vtp files.
#   --dt defaults to the median timestep spacing in the sibling .pvd if
#   present, else it must be given. --nu defaults to 0 (inviscid run: every
#   crossing is rVPM by definition; the census then reports growth-rate
#   percentiles only). --match-frac is the NN-match cutoff as a fraction of
#   the candidate particle's sigma.
#
# Output: a markdown summary to stdout (counts, dominance split, top
# offenders, matching diagnostics).

import ReadVTK
using Printf
using Statistics

function parse_args(args)
    isempty(args) && error("usage: p026_sigma_cap_census.jl <particles_dir> --sigma-cap=<m> [...]")
    dir = args[1]
    opts = Dict{String,Float64}("tol" => 0.5, "match-frac" => 0.5, "nu" => 0.0,
                                "first" => -1.0, "last" => -1.0, "dt" => NaN,
                                "sigma-cap" => NaN)
    for a in args[2:end]
        m = match(r"^--([a-z-]+)=(.*)$", a)
        m === nothing && error("unrecognized argument: $a")
        haskey(opts, m.captures[1]) || error("unknown option --$(m.captures[1])")
        opts[m.captures[1]] = parse(Float64, m.captures[2])
    end
    isnan(opts["sigma-cap"]) && error("--sigma-cap is required")
    return dir, opts
end

function vtp_series(dir)
    files = filter(f -> endswith(f, ".vtp"), readdir(dir))
    parse_idx(f) = parse(Int, match(r"\.(\d+)\.vtp$", f).captures[1])
    idxs = sort!(parse_idx.(files))
    name = replace(files[1], r"\.\d+\.vtp$" => "")
    return name, idxs
end

function load_step(dir, name, idx)
    vtk = ReadVTK.VTKFile(joinpath(dir, "$name.$idx.vtp"))
    X = Float64.(ReadVTK.get_points(vtk))
    pd = ReadVTK.get_point_data(vtk)
    sigma = Float64.(ReadVTK.get_data(pd["sigma"]))
    return X, sigma
end

# dt from the sibling .pvd timestep attributes, if available
function infer_dt(dir)
    pvd = dir * ".pvd"
    isfile(pvd) || return NaN
    ts = Float64[]
    for line in eachline(pvd)
        m = match(r"timestep=\"([0-9.eE+-]+)\"", line)
        m === nothing || push!(ts, parse(Float64, m.captures[1]))
    end
    length(ts) < 2 && return NaN
    sort!(ts)
    return median(diff(ts))
end

# Nearest-neighbor match from step-b particles back to step-a particles via a
# uniform grid hash. Returns match index per b-particle (0 = unmatched).
function match_particles(Xa, sig_a, Xb, sig_b, match_frac)
    na = size(Xa, 2); nb = size(Xb, 2)
    matches = zeros(Int, nb)
    na == 0 && return matches
    cell = max(median(sig_a), 1e-12)
    inv_cell = 1 / cell
    key(x, y, z) = (floor(Int, x * inv_cell), floor(Int, y * inv_cell),
                    floor(Int, z * inv_cell))
    grid = Dict{NTuple{3,Int},Vector{Int}}()
    for i in 1:na
        push!(get!(grid, key(Xa[1,i], Xa[2,i], Xa[3,i]), Int[]), i)
    end
    for j in 1:nb
        r_max = match_frac * sig_b[j]
        best = 0; best_d2 = r_max^2
        kx, ky, kz = key(Xb[1,j], Xb[2,j], Xb[3,j])
        reach = max(1, ceil(Int, r_max * inv_cell))
        for dx in -reach:reach, dy in -reach:reach, dz in -reach:reach
            bucket = get(grid, (kx+dx, ky+dy, kz+dz), nothing)
            bucket === nothing && continue
            for i in bucket
                d2 = (Xa[1,i]-Xb[1,j])^2 + (Xa[2,i]-Xb[2,j])^2 +
                     (Xa[3,i]-Xb[3,j])^2
                if d2 < best_d2
                    best_d2 = d2; best = i
                end
            end
        end
        matches[j] = best
    end
    return matches
end

function main(args)
    dir, opts = parse_args(args)
    name, idxs = vtp_series(dir)
    dt = isnan(opts["dt"]) ? infer_dt(joinpath(dirname(dir), basename(dir))) : opts["dt"]
    isnan(dt) && error("--dt not given and no .pvd found to infer it")
    nu = opts["nu"]; cap = opts["sigma-cap"]; tol = opts["tol"]
    lo = opts["first"] < 0 ? idxs[1] : Int(opts["first"])
    hi = opts["last"] < 0 ? idxs[end] : Int(opts["last"])
    steps = [i for i in idxs if lo <= i <= hi]
    visc_budget = 2 * nu * dt

    n_cross_visc = 0; n_cross_rvpm = 0
    n_matched_total = 0; n_particles_total = 0
    offenders = Tuple{Float64,Int,Float64,Float64}[]  # (dsigma2, step, σ_a, σ_b)
    growth_rates = Float64[]  # positive Δσ² per matched particle-step

    Xa, sig_a = load_step(dir, name, steps[1])
    for k in 2:length(steps)
        Xb, sig_b = load_step(dir, name, steps[k])
        matches = match_particles(Xa, sig_a, Xb, sig_b, opts["match-frac"])
        n_matched = count(!iszero, matches)
        n_matched_total += n_matched
        n_particles_total += length(matches)
        nsub = steps[k] - steps[k-1]  # VTPs may be every N steps
        for j in eachindex(matches)
            i = matches[j]
            i == 0 && continue
            ds2 = sig_b[j]^2 - sig_a[i]^2
            ds2 > 0 && push!(growth_rates, ds2 / nsub)
            if sig_a[i] < cap <= sig_b[j]
                if ds2 <= (1 + tol) * visc_budget * nsub && nu > 0
                    n_cross_visc += 1
                else
                    n_cross_rvpm += 1
                end
                push!(offenders, (ds2 / nsub, steps[k], sig_a[i], sig_b[j]))
            end
        end
        Xa, sig_a = Xb, sig_b
    end

    sort!(offenders; by=first, rev=true)
    match_rate = n_matched_total / max(n_particles_total, 1)

    println("# 026 W6 sigma-cap crossing census")
    println()
    println("- series: `$dir` steps $(steps[1])–$(steps[end]) ($(length(steps)) VTPs)")
    @printf("- dt = %.3e s, nu = %.3e m²/s, sigma_cap = %.4e m, tol = %.2f\n",
            dt, nu, cap, tol)
    @printf("- viscous per-step Δσ² budget 2ν·dt = %.3e m²\n", visc_budget)
    @printf("- NN match rate: %.1f%% (%d / %d particle-steps)\n",
            100match_rate, n_matched_total, n_particles_total)
    println()
    println("| crossings | viscous-dominated | rVPM-dominated |")
    println("|---|---|---|")
    println("| $(n_cross_visc + n_cross_rvpm) | $n_cross_visc | $n_cross_rvpm |")
    println()
    if !isempty(growth_rates)
        q = quantile(growth_rates, [0.5, 0.9, 0.99, 1.0])
        @printf("Positive Δσ²/step percentiles: p50 = %.3e, p90 = %.3e, p99 = %.3e, max = %.3e m²\n",
                q...)
        nu > 0 && @printf("(max / viscous budget = %.1f×)\n", q[4] / visc_budget)
        println()
    end
    if !isempty(offenders)
        println("Top crossings by Δσ²/step:")
        println()
        println("| Δσ²/step (m²) | step | σ before (m) | σ after (m) |")
        println("|---|---|---|---|")
        for (ds2, s, sa, sb) in offenders[1:min(10, end)]
            @printf("| %.3e | %d | %.4e | %.4e |\n", ds2, s, sa, sb)
        end
    else
        println("No sigma_cap crossings found in this window.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
