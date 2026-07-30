# Item 005 (Phase 2) oscillation analyzer.
#
# Reads a stable-wake run's per-step force history, isolates the post-withdrawal
# settle window, DETRENDS the rising mean with a low-order polynomial, peak-picks
# the residual, and reports period consistency + a per-cycle detrended-amplitude
# trend (growing / steady / decaying) -- the E1/E2/E3 verdict signal.
#
# Usage:
#   julia --project examples/analyze_stable_wake_oscillation.jl <run_dir_or_csv> [settle_revs] [rpm]
#
# <run_dir_or_csv> : a data/<run> directory (auto-finds monitors/*force_system1.csv,
#                    falling back to <run>_CT_vs_rev.csv) OR a direct CSV path.
# settle_revs      : analyze only the final this-many revolutions (default: all).
# rpm              : rev = time*rpm/60 when reading a force CSV (default 6000).
#
# Prints a compact report and, when several runs are analyzed by a driver, the
# key scalars (period mean/std, detrended peak-to-peak, growth-per-cycle) on a
# single tagged line for easy grep/collection.

using LinearAlgebra
using Printf

# ---- CSV loading ------------------------------------------------------------

"Return (rev::Vector, CT::Vector) from a run dir or a CSV path."
function load_ct_history(path::AbstractString; rpm::Float64=6000.0, axial_dimension::Int=1)
    csv = path
    if isdir(path)
        run = basename(rstrip(path, '/'))
        mondir = joinpath(path, "monitors")
        force_csv = nothing
        if isdir(mondir)
            cands = filter(f -> endswith(f, "force_system1.csv") && occursin("force", f),
                           readdir(mondir))
            isempty(cands) || (force_csv = joinpath(mondir, first(sort(cands))))
        end
        if force_csv !== nothing
            csv = force_csv
        else
            ct_csv = joinpath(path, "$(run)_CT_vs_rev.csv")
            isfile(ct_csv) || error("No force CSV in $mondir and no $ct_csv")
            csv = ct_csv
        end
    end

    lines = readlines(csv)
    isempty(lines) && error("empty CSV: $csv")
    header = split(strip(lines[1]), ',')
    colidx = Dict(strip(h) => i for (i, h) in enumerate(header))

    rev = Float64[]
    ct  = Float64[]
    if haskey(colidx, "CFx")            # per-step force monitor CSV
        ti = colidx["time"]; ci = colidx["CFx"]
        for ln in lines[2:end]
            isempty(strip(ln)) && continue
            f = split(strip(ln), ',')
            t = parse(Float64, f[ti])
            cfx = parse(Float64, f[ci])
            push!(rev, t * rpm / 60)
            push!(ct, -cfx)             # CT = -CFx (axial_dimension=1)
        end
    elseif haskey(colidx, "revolution") # end-of-run CT_vs_rev CSV (CT stored as -force)
        ri = colidx["revolution"]; bi = colidx["CT_bernoulli"]
        for ln in lines[2:end]
            isempty(strip(ln)) && continue
            f = split(strip(ln), ',')
            push!(rev, parse(Float64, f[ri]))
            push!(ct, parse(Float64, f[bi]))   # already CT = -CFx in this file
        end
    else
        error("unrecognized CSV columns in $csv: $(header)")
    end
    return rev, ct, csv
end

# ---- detrend + peak pick ----------------------------------------------------

"Least-squares polynomial fit (degree d) evaluated at x; returns the trend vector."
function polytrend(x::Vector{Float64}, y::Vector{Float64}, d::Int)
    n = length(x)
    d = min(d, n - 1)
    d < 1 && return fill(sum(y) / n, n)
    # center/scale x for conditioning
    x0 = sum(x) / n
    s  = maximum(abs.(x .- x0)); s = s == 0 ? 1.0 : s
    u  = (x .- x0) ./ s
    V  = hcat((u .^ k for k in 0:d)...)
    c  = V \ y
    return V * c
end

"Indices of local maxima / minima of v (strict interior extrema)."
function local_extrema(v::Vector{Float64})
    maxi = Int[]; mini = Int[]
    for i in 2:length(v)-1
        if v[i] > v[i-1] && v[i] >= v[i+1]
            push!(maxi, i)
        elseif v[i] < v[i-1] && v[i] <= v[i+1]
            push!(mini, i)
        end
    end
    return maxi, mini
end

function analyze(path::AbstractString; settle_revs::Float64=Inf, rpm::Float64=6000.0,
                 detrend_degree::Int=3, tag::AbstractString="")
    rev, ct, csv = load_ct_history(path; rpm=rpm)
    finite = isfinite.(ct)
    rev, ct = rev[finite], ct[finite]
    n = length(ct)
    n < 6 && (println("  [$(tag)] too few points ($n) to analyze"); return)

    # settle window = final settle_revs revolutions
    if isfinite(settle_revs)
        cutoff = rev[end] - settle_revs
        keep = rev .>= cutoff
        rev, ct = rev[keep], ct[keep]
    end
    n = length(ct)
    n < 6 && (println("  [$(tag)] settle window too short ($n pts)"); return)

    trend = polytrend(rev, ct, detrend_degree)
    resid = ct .- trend

    maxi, mini = local_extrema(resid)
    peak_revs = rev[maxi]
    periods = length(peak_revs) >= 2 ? diff(peak_revs) : Float64[]
    period_mean = isempty(periods) ? NaN : sum(periods) / length(periods)
    period_std  = length(periods) >= 2 ?
        sqrt(sum((periods .- period_mean).^2) / (length(periods) - 1)) : NaN

    # per-cycle peak-to-peak from consecutive max/min pairs (detrended)
    amps = Float64[]; amp_revs = Float64[]
    allext = sort(vcat([(i, true) for i in maxi], [(i, false) for i in mini]), by=first)
    for k in 1:length(allext)-1
        i1, ismax1 = allext[k]; i2, ismax2 = allext[k+1]
        ismax1 == ismax2 && continue
        push!(amps, abs(resid[i1] - resid[i2]))
        push!(amp_revs, (rev[i1] + rev[i2]) / 2)
    end

    # amplitude trend: linear fit amp vs cycle index -> growth per cycle
    growth = NaN
    if length(amps) >= 3
        idx = collect(1.0:length(amps))
        g = polytrend(idx, amps, 1)
        growth = (g[end] - g[1]) / (length(amps) - 1)
    end

    raw_ptp = maximum(ct) - minimum(ct)
    det_ptp = isempty(amps) ? (maximum(resid) - minimum(resid)) : maximum(amps)
    ncycles = length(peak_revs) >= 1 ? (rev[end] - rev[1]) / (isnan(period_mean) ? Inf : period_mean) : 0.0

    verdict = isnan(growth) ? "indeterminate" :
              growth >  0.02 * det_ptp ? "GROWING" :
              growth < -0.02 * det_ptp ? "DECAYING" : "STEADY (limit cycle)"

    println("\n=== Oscillation analysis: $(tag == "" ? csv : tag) ===")
    @printf("  window: rev %.2f -> %.2f  (%d pts), settle_revs=%s\n",
            rev[1], rev[end], n, isfinite(settle_revs) ? string(settle_revs) : "all")
    @printf("  mean CT (window)         = %.5f\n", sum(ct)/length(ct))
    @printf("  raw peak-to-peak         = %.5f\n", raw_ptp)
    @printf("  detrended peak-to-peak   = %.5f  (deg-%d trend removed)\n", det_ptp, detrend_degree)
    @printf("  detected peaks/cycles    = %d peaks, ~%.1f cycles\n", length(peak_revs), ncycles)
    @printf("  period  mean=%.3f rev  std=%.3f rev  (n=%d intervals)\n",
            period_mean, period_std, length(periods))
    if !isempty(amps)
        @printf("  per-cycle detrended ptp  = [%s]\n",
                join((@sprintf("%.4f", a) for a in amps), ", "))
    end
    @printf("  amplitude growth/cycle   = %.5f  => %s\n", growth, verdict)
    @printf("RESULT %s period_mean=%.3f period_std=%.3f det_ptp=%.5f growth=%.5f ncycles=%.1f verdict=%s\n",
            tag == "" ? "run" : tag, period_mean, period_std, det_ptp, growth, ncycles,
            replace(verdict, ' ' => '_'))
    return
end

# ---- CLI --------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("usage: analyze_stable_wake_oscillation.jl <run_dir_or_csv> [settle_revs] [rpm]")
    path = ARGS[1]
    settle_revs = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : Inf
    rpm = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 6000.0
    tag = basename(rstrip(path, '/'))
    analyze(path; settle_revs=settle_revs, rpm=rpm, tag=tag)
end
