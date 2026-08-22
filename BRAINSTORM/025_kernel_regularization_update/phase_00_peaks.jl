#!/usr/bin/env julia
# BRAINSTORM/025 Phase 0: peak |u| and |du/dh| of the three vortex-filament
# regularization families under the MATCHED-MATCHING-DISTANCE convention
# (item ruling 3): each family's core parameter is chosen so the distance at
# which its transverse profile matches the singular 1/h profile within
# relative tolerance `tol` equals a common Δr (normalized to 1).
#
# Profiles in units of Γ/(2π), transverse coordinate h (infinite straight
# filament; the 2D reduction of the segment kernel in
# src/FLOWPanel_elements_fmm.jl — see phase_00_peak_comparison.md):
#   singular            u = 1/h
#   Vatistas n=2        u = h / sqrt(h^4 + rc^4)          (as built)
#   compact-support     u = h / (h^2 + δ(h)), δ = (h-rc)^2 for h<rc else 0
#                       (the `regularize` family used by source/doublet)
#   Gaussian            u = (1 - exp(-h^2/(2σ^2))) / h    (Lamb–Oseen)
#
# Matching distances (relative error vs 1/h ≤ tol):
#   Vatistas   d = rc (2/tol)^(1/4)      [error ≈ ½(rc/h)^4]
#   compact    d = rc                    [exact beyond rc]
#   Gaussian   d = σ sqrt(2 ln(1/tol))   [error = exp(-h²/2σ²)]
#
# Run: julia --project -t 1 BRAINSTORM/025_kernel_regularization_update/phase_00_peaks.jl

using Printf

u_vat(h, rc) = h / sqrt(h^4 + rc^4)
u_cpt(h, rc) = h < rc ? h / (h^2 + (h - rc)^2) : 1 / h
u_gau(h, σ) = h < 1e-8*σ ? h / (2σ^2) : (1 - exp(-h^2 / (2σ^2))) / h

function peaks(u, param; hmax=3.0)
    hs = range(1e-6, hmax, length=2_000_001)
    umax, hstar = findmax(h -> u(h, param), hs)
    us = [u(h, param) for h in hs]
    dh = step(hs)
    # |du/dh| is the RADIAL-derivative metric (finding 4: labeled as such,
    # not as the full gradient norm)
    dumax = maximum(abs, diff(us)) / dh
    # operator norm of the transverse gradient of the infinite-filament
    # field: the nonzero entries are du/dh (radial) and u/h (curvature), so
    # the pointwise operator norm is max(|du/dh|, u/h); report its max over h
    duds = abs.(diff(us)) ./ dh
    uoverh = [us[i] / hs[i] for i in 1:length(hs)-1]
    opmax = maximum(max.(duds, uoverh))
    return umax, hs[hstar], dumax, opmax
end

println("Matched-MATCHING-DISTANCE convention (Δr normalized to 1):")
println("| tol | family | core param | u_max·Δr | at h/Δr | max|du/dh|·Δr² | max‖∇u‖_op·Δr² |")
println("|---|---|---|---:|---:|---:|---:|")
for tol in (1e-4, 1e-5, 1e-6)
    rc_v = 1 / (2 / tol)^0.25
    rc_c = 1.0
    σ_g = 1 / sqrt(2 * log(1 / tol))
    for (name, u, p) in (("Vatistas n=2", u_vat, rc_v),
                         ("compact", u_cpt, rc_c),
                         ("Gaussian", u_gau, σ_g))
        umax, hstar, dumax, opmax = peaks(u, p)
        @printf("| %.0e | %s | %.4g | %.3f | %.3f | %.2f | %.2f |\n",
            tol, name, p, umax, hstar, dumax, opmax)
    end
end

println("\nMatched-CORE-SIZE convention (rc normalized to 1; governs per Ryan 2026-08-20):")
println("| family | u_max·rc | at h/rc | max|du/dh|·rc² | max‖∇u‖_op·rc² |")
println("|---|---:|---:|---:|---:|")
for (name, u) in (("Vatistas n=2", u_vat), ("compact", u_cpt), ("Gaussian", u_gau))
    umax, hstar, dumax, opmax = peaks(u, 1.0)
    @printf("| %s | %.3f | %.3f | %.2f | %.2f |\n", name, umax, hstar, dumax, opmax)
end

# far-field sanity: relative error at h = Δr = 1 must be ≤ tol
println("\nFar-field check (rel err vs 1/h at h = Δr):")
for tol in (1e-4, 1e-5, 1e-6)
    rc_v = 1 / (2 / tol)^0.25
    σ_g = 1 / sqrt(2 * log(1 / tol))
    @printf("tol %.0e: vatistas %.2e  compact %.2e  gaussian %.2e\n", tol,
        abs(u_vat(1.0, rc_v) - 1), abs(u_cpt(1.0, 1.0) - 1), abs(u_gau(1.0, σ_g) - 1))
end
