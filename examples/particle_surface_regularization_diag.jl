# particle_surface_regularization_diag.jl
#
# BRAINSTORM 018 Phase 12 C1 (spec: phase_12_spatial_rigor.md): offline
# measurement of the particle->surface velocity deficit caused by kernel
# regularization. Loads a saved settled wake state (particle VTP + body VTU at
# the same step), evaluates the wake-particle-induced velocity at every body
# panel centroid twice -- the run's regularized Gaussian-erf kernel at each
# particle's actual sigma, and the singular Biot-Savart limit (g == 1) -- and
# reports the relative deficit per panel against its nearest-particle d/sigma.
#
# Deliverables (written to OUT):
#   per_panel.csv   one row per body panel: centroid, r/R, min d/sigma,
#                   dominant-contributor d/sigma, |u_sing|, deficit metrics
#   binned.csv      log-binned median/p95 deficit vs min d/sigma + local slope
#   summary.txt     admissible thresholds d/sigma* (median and p95 at 0.5%),
#                   fitted log-log slope vs the 017 M4 static prior (~ -3.4)
#   deficit_vs_dsigma.png
#
# The single-particle deficit is EXACTLY 1 - g(d/sigma) (selftest gate); the
# many-particle field measures how much cancellation/summation moves the
# surface-seen deficit off that curve.
#
# Usage:
#   julia --project -t 6 examples/particle_surface_regularization_diag.jl
#   SELFTEST=1 julia --project examples/particle_surface_regularization_diag.jl
#
# ENV knobs: STATE_DIR (default ~/p018_L1_ov3_paraview), RUN (p018_L1_ov3),
# STEP (719), OUT (data/particle_surface_regularization_diag), MAKE_PLOT (true).

import ReadVTK
import FLOWVPM
import Random

const DO_PLOT = get(ENV, "SELFTEST", "0") != "1" &&
                get(ENV, "MAKE_PLOT", "true") == "true"
if DO_PLOT
    import PythonPlot as plt
end
import Statistics: median, quantile
import Printf: @sprintf
using Base.Threads: @threads, nthreads

# Regularization function g(rho): fraction of a particle's vorticity-induced
# swirl recovered at scaled radius rho = r/sigma; g -> 1 is the singular limit.
# Must match the run's velocity kernel (FLOWVPM default Gaussian-erf).
const G = FLOWVPM.g_gauserf

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

"Read particle positions/strengths/core sizes from a saved wake VTP."
function load_particles(vtp_path::AbstractString)
    isfile(vtp_path) || error("Particle VTP not found: $(vtp_path)")
    vtk = ReadVTK.VTKFile(vtp_path)
    np = vtk.n_points
    pd = ReadVTK.get_point_data(vtk)
    X = Matrix{Float64}(ReadVTK.get_points(vtk))                        # 3 x np
    gamma = Matrix{Float64}(reshape(ReadVTK.get_data(pd["gamma"]), 3, np))
    sigma = Vector{Float64}(vec(ReadVTK.get_data(pd["sigma"])))
    keep = [!iszero(gamma[1, i]) || !iszero(gamma[2, i]) || !iszero(gamma[3, i])
            for i in 1:np]
    return (; X = X[:, keep], gamma = gamma[:, keep], sigma = sigma[keep],
            np = count(keep), ndropped = np - count(keep))
end

"Read body panel centroids and unit normals from a saved body VTU."
function load_body(vtu_path::AbstractString)
    isfile(vtu_path) || error("Body VTU not found: $(vtu_path)")
    vtk = ReadVTK.VTKFile(vtu_path)
    pts = Matrix{Float64}(ReadVTK.get_points(vtk))                      # 3 x npts
    cells = ReadVTK.to_meshcells(ReadVTK.get_cells(vtk))
    ncell = length(cells)
    centroids = zeros(3, ncell)
    for (j, cell) in enumerate(cells)
        for i in cell.connectivity
            @views centroids[:, j] .+= pts[:, i]
        end
        centroids[:, j] ./= length(cell.connectivity)
    end
    cd = ReadVTK.get_cell_data(vtk)
    normals = Matrix{Float64}(reshape(ReadVTK.get_data(cd["normals"]), 3, ncell))
    for j in 1:ncell
        n = @views sqrt(sum(abs2, normals[:, j]))
        n > 0 && (@views normals[:, j] ./= n)
    end
    R = maximum(sqrt(pts[1, i]^2 + pts[2, i]^2) for i in axes(pts, 2))
    return (; centroids, normals, ncell, R)
end

# ---------------------------------------------------------------------------
# Core evaluation
# ---------------------------------------------------------------------------

"""
Evaluate at each target the regularized and singular particle-induced
velocities. Returns per-target u_sing (3 x n), du = u_reg - u_sing (3 x n),
rho_min = min_p |x - X_p|/sigma_p, and rho_dom = rho of the particle
contributing the largest |deficit| at that target.
"""
function evaluate_deficit(targets::AbstractMatrix, X::AbstractMatrix,
                          gamma::AbstractMatrix, sigma::AbstractVector)
    n = size(targets, 2)
    np = size(X, 2)
    u_sing = zeros(3, n)
    du = zeros(3, n)
    rho_min = fill(Inf, n)
    rho_dom = fill(NaN, n)
    @threads for j in 1:n
        xj = (targets[1, j], targets[2, j], targets[3, j])
        us1 = us2 = us3 = 0.0
        dd1 = dd2 = dd3 = 0.0
        rmin = Inf
        dommag = -1.0
        domrho = NaN
        @inbounds for p in 1:np
            r1 = xj[1] - X[1, p]; r2 = xj[2] - X[2, p]; r3 = xj[3] - X[3, p]
            r2sum = r1*r1 + r2*r2 + r3*r3
            r = sqrt(r2sum)
            r < eps() && continue        # target on top of a particle
            rho = r / sigma[p]
            rho < rmin && (rmin = rho)
            g1 = G(rho)
            inv4pir3 = 1.0 / (4pi * r2sum * r)
            # cross(r, Gamma_p)
            c1 = r2*gamma[3, p] - r3*gamma[2, p]
            c2 = r3*gamma[1, p] - r1*gamma[3, p]
            c3 = r1*gamma[2, p] - r2*gamma[1, p]
            us1 -= inv4pir3 * c1; us2 -= inv4pir3 * c2; us3 -= inv4pir3 * c3
            w = (g1 - 1.0) * inv4pir3
            dd1 -= w * c1; dd2 -= w * c2; dd3 -= w * c3
            mag = abs(w) * sqrt(c1*c1 + c2*c2 + c3*c3)
            mag > dommag && (dommag = mag; domrho = rho)
        end
        u_sing[1, j] = us1; u_sing[2, j] = us2; u_sing[3, j] = us3
        du[1, j] = dd1; du[2, j] = dd2; du[3, j] = dd3
        rho_min[j] = rmin
        rho_dom[j] = domrho
    end
    return (; u_sing, du, rho_min, rho_dom)
end

# ---------------------------------------------------------------------------
# Selftest: single-particle deficit must equal 1 - g(rho) exactly
# ---------------------------------------------------------------------------

function selftest()
    ok = true
    # Gate 1: g saturates to 1 far away
    ok &= abs(G(50.0) - 1.0) < 1e-12
    # Gate 2: single particle, exact deficit ratio and singular magnitude
    sigma = 0.137
    gamma = reshape([0.3, -0.2, 0.9], 3, 1)
    X = zeros(3, 1)
    for rho in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
        d = rho * sigma
        tgt = reshape([d, 0.0, 0.0], 3, 1)
        r = evaluate_deficit(tgt, X, gamma, [sigma])
        rel = sqrt(sum(abs2, r.du)) / sqrt(sum(abs2, r.u_sing))
        expect = 1.0 - G(rho)
        e1 = abs(rel - expect) / expect
        # |u_sing| = |r x Gamma| / (4 pi d^3) with r = (d,0,0)
        umag = sqrt(gamma[2]^2 + gamma[3]^2) * d / (4pi * d^3)
        e2 = abs(sqrt(sum(abs2, r.u_sing)) - umag) / umag
        e3 = abs(r.rho_min[1] - rho) / rho
        okk = e1 < 1e-12 && e2 < 1e-12 && e3 < 1e-12
        println(@sprintf("  rho=%5.2f  rel=%.6e  1-g=%.6e  err=%.2e  |u_s| err=%.2e  %s",
                         rho, rel, expect, e1, e2, okk ? "PASS" : "FAIL"))
        ok &= okk
    end
    # Gate 3: superposition -- 40 random particles vs naive re-sum
    rng = Random.MersenneTwister(4)
    Xr = randn(rng, 3, 40); Gr = randn(rng, 3, 40)
    sr = 0.05 .+ 0.2 .* rand(rng, 40)
    tgts = 3.0 .* randn(rng, 3, 7)
    r = evaluate_deficit(tgts, Xr, Gr, sr)
    maxerr = 0.0
    for j in 1:7
        us = zeros(3); ur = zeros(3)
        for p in 1:40
            rv = tgts[:, j] .- Xr[:, p]
            rr = sqrt(sum(abs2, rv))
            c = -[rv[2]*Gr[3,p]-rv[3]*Gr[2,p], rv[3]*Gr[1,p]-rv[1]*Gr[3,p],
                  rv[1]*Gr[2,p]-rv[2]*Gr[1,p]] ./ (4pi*rr^3)
            us .+= c
            ur .+= G(rr/sr[p]) .* c
        end
        maxerr = max(maxerr,
                     maximum(abs.(r.u_sing[:, j] .- us)) / maximum(abs.(us)),
                     maximum(abs.((r.du[:, j] .+ us) .- ur)) / maximum(abs.(ur)))
    end
    println(@sprintf("  superposition max rel err = %.2e  %s", maxerr,
                     maxerr < 1e-12 ? "PASS" : "FAIL"))
    ok &= maxerr < 1e-12
    println(ok ? "SELFTEST PASS" : "SELFTEST FAIL")
    return ok
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    state_dir = get(ENV, "STATE_DIR", joinpath(homedir(), "p018_L1_ov3_paraview"))
    run_name  = get(ENV, "RUN", "p018_L1_ov3")
    step      = parse(Int, get(ENV, "STEP", "719"))
    out       = get(ENV, "OUT", "data/particle_surface_regularization_diag")
    mkpath(out)

    vtp = joinpath(state_dir, "$(run_name)_wake1_particles",
                   "$(run_name)_wake1_particles.$(step).vtp")
    vtu = joinpath(state_dir, "$(run_name)_body1", "$(run_name)_body1.$(step).vtu")

    part = load_particles(vtp)
    body = load_body(vtu)
    # Campaign fixed operating point (every 018 run): 5400 RPM.
    Vtip = 5400 * 2pi / 60 * body.R

    println("state: $(run_name) step $(step)  particles=$(part.np) " *
            "(dropped $(part.ndropped) zero-strength)  panels=$(body.ncell)")
    println("R=$(round(body.R, digits=4)) m  Vtip=$(round(Vtip, digits=2)) m/s  " *
            "sigma: min=$(round(minimum(part.sigma), sigdigits=3)) " *
            "med=$(round(median(part.sigma), sigdigits=3)) " *
            "max=$(round(maximum(part.sigma), sigdigits=3))  threads=$(nthreads())")

    t = @elapsed res = evaluate_deficit(body.centroids, part.X, part.gamma, part.sigma)
    println(@sprintf("pairwise evaluation: %.1f s (%.2e pair evals)", t,
                     Float64(part.np) * body.ncell))

    # Per-panel metrics
    n = body.ncell
    usmag = [sqrt(sum(abs2, @view res.u_sing[:, j])) for j in 1:n]
    dumag = [sqrt(sum(abs2, @view res.du[:, j])) for j in 1:n]
    dun   = [abs(sum(res.du[i, j] * body.normals[i, j] for i in 1:3)) for j in 1:n]
    rel   = dumag ./ usmag
    reln  = dun ./ usmag
    rtip  = [sqrt(body.centroids[1, j]^2 + body.centroids[2, j]^2) / body.R for j in 1:n]

    open(joinpath(out, "per_panel.csv"), "w") do io
        println(io, "panel,x,y,z,r_over_R,rho_min,rho_dom,u_sing,rel_deficit,rel_deficit_normal,du_over_Vtip")
        for j in 1:n
            println(io, @sprintf("%d,%.6e,%.6e,%.6e,%.4f,%.4f,%.4f,%.6e,%.6e,%.6e,%.6e",
                                 j, body.centroids[1, j], body.centroids[2, j],
                                 body.centroids[3, j], rtip[j], res.rho_min[j],
                                 res.rho_dom[j], usmag[j], rel[j], reln[j],
                                 dumag[j] / Vtip))
        end
    end

    # Log-binned medians/p95 vs min d/sigma
    lo, hi = extrema(res.rho_min)
    edges = exp.(range(log(lo * 0.999), log(hi * 1.001), length = 25))
    rows = NamedTuple[]
    for b in 1:length(edges)-1
        sel = [j for j in 1:n if edges[b] <= res.rho_min[j] < edges[b+1]]
        isempty(sel) && continue
        push!(rows, (; rho = sqrt(edges[b] * edges[b+1]),
                     med = median(rel[sel]), p95 = quantile(rel[sel], 0.95),
                     med_n = median(reln[sel]), count = length(sel)))
    end
    open(joinpath(out, "binned.csv"), "w") do io
        println(io, "rho_center,median_rel_deficit,p95_rel_deficit,median_rel_deficit_normal,count")
        for r in rows
            println(io, @sprintf("%.4f,%.6e,%.6e,%.6e,%d", r.rho, r.med, r.p95,
                                 r.med_n, r.count))
        end
    end

    # Thresholds: smallest rho above which every bin stays under 0.5%
    function threshold(vals)
        thr = NaN
        for k in eachindex(rows)
            if all(vals[m] < 0.005 for m in k:length(rows))
                thr = rows[k].rho
                break
            end
        end
        return thr
    end
    thr_med = threshold([r.med for r in rows])
    thr_p95 = threshold([r.p95 for r in rows])

    # Log-log slope of the median curve (per-bin local + global LSQ), for
    # comparison with the 017 M4 static prior ~ (d/sigma)^-3.4. The
    # single-particle law 1-g decays super-polynomially, so the slope steepens
    # with rho; report the global fit AND the slope near the threshold.
    fitsel = [k for k in eachindex(rows) if rows[k].med > 0 && rows[k].count >= 10]
    lx = [log(rows[k].rho) for k in fitsel]
    ly = [log(rows[k].med) for k in fitsel]
    mx = sum(lx) / length(lx); my = sum(ly) / length(ly)
    slope = sum((lx .- mx) .* (ly .- my)) / sum(abs2, lx .- mx)

    summary = """
    particle_surface_regularization_diag -- $(run_name) step $(step)
    particles=$(part.np)  panels=$(body.ncell)  R=$(body.R)  Vtip=$(Vtip)
    min d/sigma over panels: $(round(lo, digits=3));  max: $(round(hi, digits=3))
    global rel deficit: median=$(round(median(rel), sigdigits=3))  p95=$(round(quantile(rel, 0.95), sigdigits=3))  max=$(round(maximum(rel), sigdigits=3))
    admissible threshold d/sigma* (deficit < 0.5% of local |u_sing|):
      by binned MEDIAN: $(round(thr_med, digits=2))
      by binned P95:    $(round(thr_p95, digits=2))
    global log-log slope of median deficit vs d/sigma: $(round(slope, digits=2))
      (017 M4 static prior: ~ -3.4; single-particle 1-g decays faster than any
       power, so a shallower global slope means many-particle summation, not
       the kernel, controls the tail)
    """
    write(joinpath(out, "summary.txt"), summary)
    println(summary)

    if DO_PLOT
        fig, ax = plt.subplots(figsize = (6.4, 4.2))
        rho = [r.rho for r in rows]
        ax.loglog(rho, [r.med for r in rows]; color = "#2a78d6", lw = 2,
                  label = "median over panels")
        ax.loglog(rho, [r.p95 for r in rows]; color = "#eb6834", lw = 2,
                  label = "95th percentile")
        rr = exp.(range(log(minimum(rho)), log(maximum(rho)), length = 200))
        ax.loglog(rr, 1 .- G.(rr); color = "#1baf7a", lw = 2, ls = "--",
                  label = "single particle 1−g(d/σ)")
        ax.axhline(0.005; color = "0.5", lw = 1, ls = ":")
        isfinite(thr_p95) && ax.axvline(thr_p95; color = "0.5", lw = 1, ls = ":")
        ax.set_xlabel("min particle clearance d/σ")
        ax.set_ylabel("relative velocity deficit |Δu|/|u_sing|")
        ax.set_title("Particle→surface regularization deficit ($(run_name) step $(step))")
        ax.legend(frameon = false)
        ax.grid(true, which = "both", alpha = 0.2, lw = 0.5)
        fig.tight_layout()
        fig.savefig(joinpath(out, "deficit_vs_dsigma.png"), dpi = 160)
        println("figure: $(joinpath(out, "deficit_vs_dsigma.png"))")
    end
end

if get(ENV, "SELFTEST", "0") == "1"
    selftest() || exit(1)
else
    main()
end
