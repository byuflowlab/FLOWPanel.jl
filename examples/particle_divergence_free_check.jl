# Numerical check: is the vortex-particle-induced velocity field divergence-free?
#
# BRAINSTORM item 009. We differentiate the regularized (Gaussian-erf) induced
# velocity kernel with ForwardDiff to form J_ab = du_a/dx_b and report tr(J),
# which is exactly div(u). Theory (see 009): each radial particle kernel gives a
# divergence-free velocity for any overlap, so tr(J) should sit at Float64
# round-off everywhere -- at particle centers, in the overlap zone, and far away.
#
# Run: julia --project=test examples/particle_divergence_free_check.jl

import ForwardDiff
using LinearAlgebra
using Random
using Printf

# --- Gaussian-erf regularization, matching FLOWVPM's default kernel -----------
# erf: prefer the exact SpecialFunctions.erf (has ForwardDiff rules); otherwise
# fall back to an AD-clean series + asymptotic form (~1e-15 over our sample range).
const HAVE_SF = let
    try
        @eval import SpecialFunctions
        true
    catch
        false
    end
end

if HAVE_SF
    _erf(x) = SpecialFunctions.erf(x)
else
    function _erf(x::T) where {T}
        z = abs(x)
        if z < 4
            # Maclaurin: erf(x) = 2/sqrt(pi) * sum_{n>=0} (-1)^n x^{2n+1}/(n!(2n+1))
            s = z; term = z; n = 0
            while abs(term) > 1e-18 * max(abs(s), one(T)) && n < 200
                n += 1
                term *= -z^2 / n
                s += term / (2n + 1)
            end
            e = 2 / sqrt(T(pi)) * s
        else
            e = one(T) - exp(-z^2) / (z * sqrt(T(pi)))  # asymptotic erfc
        end
        return x < 0 ? -e : e
    end
end

# g(rho) = erf(rho/sqrt2) - sqrt(2/pi) * rho * exp(-rho^2/2)
g_gauserf(rho) = _erf(rho / sqrt(2)) - sqrt(2 / pi) * rho * exp(-rho^2 / 2)

# phi(rho) = g(rho)/rho^3 -- the kernel's radial profile. It is a SMOOTH, even
# function of rho (finite at rho=0), but written as g/rho^3 it is 0/0 there, and
# rho=sqrt(r.r) is non-differentiable at r=0. To keep ForwardDiff clean at the
# particle center, express phi as a function of rho2 = rho^2 and use the Taylor
# series for small rho2:  phi = sqrt(2/pi)(1/3 - rho2/10 + rho2^2/56 - ...).
function phi_gauserf(rho2)
    if rho2 < 1e-4
        return sqrt(2 / pi) * (1 / 3 - rho2 / 10 + rho2^2 / 56)
    else
        rho = sqrt(rho2)
        return g_gauserf(rho) / (rho2 * rho)   # g/rho^3
    end
end

# --- Induced velocity from a set of particles ---------------------------------
# u(x) = -1/(4pi) sum_p g(|x-Xp|/sigma_p)/|x-Xp|^3 * (x-Xp) x Gamma_p
function induced_velocity(x, Xs, Gammas, sigmas)
    u = zeros(eltype(x), 3)
    @inbounds for p in eachindex(sigmas)
        r = x .- Xs[p]
        rho2 = (r[1]^2 + r[2]^2 + r[3]^2) / sigmas[p]^2
        f = phi_gauserf(rho2) / sigmas[p]^3            # = g(rho)/|r|^3
        cross = (r[2] * Gammas[p][3] - r[3] * Gammas[p][2],
                 r[3] * Gammas[p][1] - r[1] * Gammas[p][3],
                 r[1] * Gammas[p][2] - r[2] * Gammas[p][1])
        u = u .- (1 / (4pi)) * f .* collect(cross)
    end
    return u
end

# Deliberately-broken control: anisotropic (non-radial) core. Should NOT be
# divergence-free, proving the diagnostic has teeth.
function induced_velocity_broken(x, Xs, Gammas, sigmas; scale=(1.0, 2.0, 0.5))
    u = zeros(eltype(x), 3)
    @inbounds for p in eachindex(sigmas)
        r = x .- Xs[p]
        # anisotropic radius -> non-radial mollifier: K is no longer f(|r|)*r, so
        # curl(K) != 0 and the field is no longer divergence-free.
        rho2 = ((r[1] / scale[1])^2 + (r[2] / scale[2])^2 + (r[3] / scale[3])^2) / sigmas[p]^2
        f = phi_gauserf(rho2) / sigmas[p]^3
        cross = (r[2] * Gammas[p][3] - r[3] * Gammas[p][2],
                 r[3] * Gammas[p][1] - r[1] * Gammas[p][3],
                 r[1] * Gammas[p][2] - r[2] * Gammas[p][1])
        u = u .- (1 / (4pi)) * f .* collect(cross)
    end
    return u
end

# --- Divergence and curl via ForwardDiff --------------------------------------
function div_curl(velfun, x)
    J = ForwardDiff.jacobian(velfun, x)          # J[a,b] = du_a/dx_b
    divu = J[1, 1] + J[2, 2] + J[3, 3]           # tr(J) = div u
    curl = (J[3, 2] - J[2, 3], J[1, 3] - J[3, 1], J[2, 1] - J[1, 2])
    return divu, norm(J), sqrt(curl[1]^2 + curl[2]^2 + curl[3]^2)
end

# --- Test driver --------------------------------------------------------------
function run_regime(name, d; npts=400, broken=false, rng=Random.MersenneTwister(1))
    sigma = 1.0
    Xs = [(-d / 2, 0.0, 0.0), (d / 2, 0.0, 0.0)]
    Gammas = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.5)]
    sigmas = [sigma, sigma]

    velfun = broken ?
        (x -> induced_velocity_broken(x, Xs, Gammas, sigmas)) :
        (x -> induced_velocity(x, Xs, Gammas, sigmas))

    # Evaluation points: the two centers + random cloud spanning cores..far field.
    pts = Vector{Vector{Float64}}()
    push!(pts, collect(Xs[1]))
    push!(pts, collect(Xs[2]))
    L = max(d, 2 * sigma)
    for _ in 1:npts
        push!(pts, [(-1.5L + 3L * rand(rng)) for _ in 1:3])
    end

    max_absdiv = 0.0
    max_reldiv = 0.0
    min_curl = Inf
    sum_curl = 0.0
    for x in pts
        divu, Jn, curl = div_curl(velfun, x)
        reldiv = abs(divu) / (Jn + eps())
        max_absdiv = max(max_absdiv, abs(divu))
        max_reldiv = max(max_reldiv, reldiv)
        min_curl = min(min_curl, curl)
        sum_curl += curl
    end
    return (; name, d, npts=length(pts), max_absdiv, max_reldiv,
            min_curl, mean_curl=sum_curl / length(pts))
end

function main()
    println("Particle-induced velocity divergence check (ForwardDiff through kernel)")
    println("erf source: ", HAVE_SF ? "SpecialFunctions.erf" : "self-contained series")
    println()

    regimes = [("overlap a lot  (d=0.25 sigma)", 0.25),
               ("overlap a little (d=2 sigma)", 2.0),
               ("far apart        (d=20 sigma)", 20.0)]

    TOL = 1e-10   # Float64 round-off floor for tr(J)/||J|| (broken control is O(1))
    allpass = true
    rows = String[]

    @printf("%-32s %5s  %12s  %12s  %10s\n",
            "regime", "npts", "max|divu|", "max rel div", "mean|curl|")
    println("-"^80)
    for (nm, d) in regimes
        r = run_regime(nm, d)
        pass = r.max_reldiv < TOL && r.mean_curl > 1e-6
        allpass &= pass
        @printf("%-32s %5d  %12.3e  %12.3e  %10.3e  %s\n",
                nm, r.npts, r.max_absdiv, r.max_reldiv, r.mean_curl,
                pass ? "PASS" : "FAIL")
        push!(rows, @sprintf("%s,%g,%d,%.6e,%.6e,%.6e,PASS=%s",
                             nm, d, r.npts, r.max_absdiv, r.max_reldiv, r.mean_curl, pass))
    end

    println()
    println("Broken control (anisotropic core -> should FAIL divergence test):")
    rb = run_regime("BROKEN anisotropic (d=2 sigma)", 2.0; broken=true)
    broken_ok = rb.max_reldiv > 1e-6   # we WANT this to be non-tiny
    @printf("%-32s %5d  %12.3e  %12.3e  %10.3e  %s\n",
            rb.name, rb.npts, rb.max_absdiv, rb.max_reldiv, rb.mean_curl,
            broken_ok ? "div != 0 (good: metric has teeth)" : "UNEXPECTEDLY ZERO")
    push!(rows, @sprintf("%s,%g,%d,%.6e,%.6e,%.6e,brokenOK=%s",
                         rb.name, 2.0, rb.npts, rb.max_absdiv, rb.max_reldiv, rb.mean_curl, broken_ok))

    # Write CSV
    outdir = joinpath(@__DIR__, "..", "data", "particle_divergence_free")
    mkpath(outdir)
    open(joinpath(outdir, "divergence_report.csv"), "w") do io
        println(io, "regime,d_over_sigma,npts,max_absdiv,max_reldiv,mean_curl,flag")
        for row in rows
            println(io, row)
        end
    end

    println()
    final = allpass && broken_ok
    println(final ?
        "RESULT: PASS -- field is divergence-free to round-off for all overlaps, and the metric is non-vacuous." :
        "RESULT: FAIL -- investigate (see BRAINSTORM/009).")
    println("Wrote ", joinpath(outdir, "divergence_report.csv"))
    return final
end

main()
