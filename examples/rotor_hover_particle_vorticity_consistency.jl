#=##############################################################################
# DESCRIPTION (T4 of plans/20260709_pressure_monitor_reliability.md)
#   At one saved step (default 1200) of the settled rotor-hover cycle, compare
#   three candidate "particle-induced vorticity at the body control points":
#     (i)   FMM extra_outputs=3 induced vorticity, particle sources only
#           (this is exactly what feeds the lamb-vector PressureLaplace RHS,
#           minus the bound-sheet kappa added separately) — evaluated with both
#           the FastMultipoleBackend used in production and a DirectBackend to
#           separate FMM approximation error from kernel identity questions;
#     (ii)  direct basis sum  omega(x) = sum_p Gamma_p * zeta_sigma(|x-x_p|)
#           (the naive "sum of regularized vorticity blobs");
#     (iii) curl of the particle-induced velocity from the FMM Hessian
#           (antisymmetric part of body.velocity_gradient).
#   For a regularized particle with vector potential psi = G_sigma * Gamma,
#     curl(u) = zeta_sigma*Gamma + Hess(G_sigma)*Gamma,
#   so (iii) - (ii) should equal the analytic Hess(G)Gamma correction
#   (the Leray-projection discrepancy established in
#   examples/particle_vorticity_curl_vs_basis_check.jl). Which of (ii)/(iii)
#   the extra-output (i) matches is the empirical question this script answers.
#
#   Outputs: printed norm table + per-radial-bin RMS profile CSV under
#   data/rotor_hover_pressure_reliability/.
#
# Usage: julia -t auto --project examples/rotor_hover_particle_vorticity_consistency.jl
=###############################################################################

import FLOWPanel as pnl
import FLOWPanel.FLOWVPM as vpm
using Printf
using Statistics

include(joinpath(@__DIR__, "particle_gauserf_hessian_helpers.jl"))
using .ParticleGauserfHessian: zeta_sigma, hess_green_times_gamma_components

# ---------------- configuration ----------------------------------------------

const DATA_PATH = get(ENV, "T4_PATH", joinpath("data", "rotor_hover_relaxfilter2p0_ws"))
const RUN_NAME  = get(ENV, "T4_NAME", "rotor_hover_relaxfilter2p0_ws")
const STEP      = parse(Int, get(ENV, "T4_STEP", "1200"))
const OUTDIR    = joinpath("data", "rotor_hover_pressure_reliability")
const ROTOR_R   = 0.119
const NBINS     = 24

mkpath(OUTDIR)

# ---------------- load saved state (no recompute, no monitors) ----------------

println("Loading step $(STEP) of $(RUN_NAME) ...")
result = pnl.replay(DATA_PATH, RUN_NAME; steps=[STEP], recompute=())
body = result.systems[1]
wakes = result.wakes

wake_sources = pnl._collect_wake_sources(wakes)
particle_sources = Tuple(s for s in wake_sources if s isa vpm.ParticleField)
@assert length(particle_sources) == 1 "expected exactly one particle field, got $(length(particle_sources))"
pf = particle_sources[1]
np = pf.np
println("np = $(np) particles, kernel = $(typeof(pf.kernel))")

cps = body.controlpoints
ncp = size(cps, 2)
println("ncp = $(ncp) control points")

X     = copy(view(pf.particles, vpm.X_INDEX, 1:np))       # 3 x np
GAMMA = copy(view(pf.particles, vpm.GAMMA_INDEX, 1:np))   # 3 x np
SIGMA = copy(vec(view(pf.particles, vpm.SIGMA_INDEX, 1:np)))

# ---------------- (i) FMM extra-output induced vorticity ----------------------

function extra_output_omega(body, sources, backend)
    body.induced_vorticity .= 0
    pnl.influence!((body,), sources, backend; precalc=true,
        scalar_potential=false, velocity=false,
        velocity_gradient=(false,),
        extra_outputs=(3,))
    return copy(body.induced_vorticity)
end

println("\n(i) extra_outputs=3, FastMultipoleBackend ...")
t = @elapsed omega_fmm = extra_output_omega(body, particle_sources, pnl.FastMultipoleBackend())
@printf("    %.1f s\n", t)

println("(i') extra_outputs=3, DirectBackend ...")
t = @elapsed omega_direct_kernel = extra_output_omega(body, particle_sources, pnl.DirectBackend())
@printf("    %.1f s\n", t)

# ---------------- (ii) direct basis sum  sum_p Gamma_p zeta_sigma -------------

println("(ii) direct basis sum ...")
omega_basis = zeros(3, ncp)
t = @elapsed Threads.@threads for i in 1:ncp
    xi, yi, zi = cps[1, i], cps[2, i], cps[3, i]
    s1 = s2 = s3 = 0.0
    @inbounds for p in 1:np
        rx = xi - X[1, p]
        ry = yi - X[2, p]
        rz = zi - X[3, p]
        sig = SIGMA[p]
        rho2 = (rx * rx + ry * ry + rz * rz) / (sig * sig)
        rho2 > 60 && continue
        z = zeta_sigma(rho2, sig)
        s1 += z * GAMMA[1, p]
        s2 += z * GAMMA[2, p]
        s3 += z * GAMMA[3, p]
    end
    omega_basis[1, i] = s1
    omega_basis[2, i] = s2
    omega_basis[3, i] = s3
end
@printf("    %.1f s\n", t)

# ---------------- (iii) curl from FMM Hessian ---------------------------------

function hessian_curl_omega(body, sources, backend)
    body.velocity_gradient .= 0
    pnl.influence!((body,), sources, backend; precalc=true,
        scalar_potential=false, velocity=false,
        velocity_gradient=(true,),
        extra_outputs=(0,))
    G = body.velocity_gradient
    out = zeros(3, size(G, 3))
    @inbounds for i in axes(G, 3)
        out[1, i] = G[3, 2, i] - G[2, 3, i]
        out[2, i] = G[1, 3, i] - G[3, 1, i]
        out[3, i] = G[2, 1, i] - G[1, 2, i]
    end
    return out
end

println("(iii) Hessian curl, FastMultipoleBackend ...")
t = @elapsed omega_curl = hessian_curl_omega(body, particle_sources, pnl.FastMultipoleBackend())
@printf("    %.1f s\n", t)

println("(iii') Hessian curl, DirectBackend ...")
t = @elapsed omega_curl_direct = hessian_curl_omega(body, particle_sources, pnl.DirectBackend())
@printf("    %.1f s\n", t)

# ---------------- analytic Hess(G)Gamma correction -----------------------------
# curl(u_sigma) = zeta*Gamma + Hess(G_sigma)*Gamma per particle, so
# (iii) - (ii) should equal corr = sum_p Hess(G_sigma_p)(x - x_p) * Gamma_p
# (up to the sign convention baked into the helper; both signs reported).

println("analytic Hess(G)Gamma correction ...")
corr = zeros(3, ncp)
t = @elapsed Threads.@threads for i in 1:ncp
    xi, yi, zi = cps[1, i], cps[2, i], cps[3, i]
    c1 = c2 = c3 = 0.0
    @inbounds for p in 1:np
        rx = xi - X[1, p]
        ry = yi - X[2, p]
        rz = zi - X[3, p]
        q1, q2, q3 = hess_green_times_gamma_components(
            rx, ry, rz, SIGMA[p], GAMMA[1, p], GAMMA[2, p], GAMMA[3, p])
        c1 += q1
        c2 += q2
        c3 += q3
    end
    corr[1, i] = c1
    corr[2, i] = c2
    corr[3, i] = c3
end
@printf("    %.1f s\n", t)

# ---------------- norms + verdicts ---------------------------------------------

colnorms(A) = vec(sqrt.(sum(abs2, A; dims=1)))
rms(v) = sqrt(mean(abs2, v))

fields = [
    ("(i)   omega_fmm (extra_outputs, FMM)", omega_fmm),
    ("(i')  omega_direct_kernel (extra_outputs, direct)", omega_direct_kernel),
    ("(ii)  omega_basis (sum Gamma zeta)", omega_basis),
    ("(iii) omega_curl (Hessian antisym, FMM)", omega_curl),
    ("(iii')omega_curl_direct (Hessian antisym, direct)", omega_curl_direct),
    ("corr  analytic Hess(G)Gamma", corr),
]

println("\n=== field magnitudes (|omega| over $(ncp) control points) ===")
@printf("%-52s %12s %12s %12s\n", "field", "rms", "max", "median")
for (label, A) in fields
    m = colnorms(A)
    @printf("%-52s %12.4e %12.4e %12.4e\n", label, rms(m), maximum(m), median(m))
end

diffs = [
    ("(i)-(i')   FMM approximation error", omega_fmm .- omega_direct_kernel),
    ("(iii)-(iii') FMM Hessian approx error", omega_curl .- omega_curl_direct),
    ("(i')-(ii)  extra_output vs basis sum", omega_direct_kernel .- omega_basis),
    ("(i')-(iii') extra_output vs true curl", omega_direct_kernel .- omega_curl_direct),
    ("(iii')-(ii) curl vs basis (Leray gap)", omega_curl_direct .- omega_basis),
    ("(iii')-(ii)-corr", omega_curl_direct .- omega_basis .- corr),
    ("(iii')-(ii)+corr", omega_curl_direct .- omega_basis .+ corr),
    ("(i')-(ii)-corr", omega_direct_kernel .- omega_basis .- corr),
    ("(i')-(ii)+corr", omega_direct_kernel .- omega_basis .+ corr),
]

omega_ref = rms(colnorms(omega_basis))
println("\n=== differences (rms, relative to rms|omega_basis| = $(@sprintf("%.4e", omega_ref))) ===")
@printf("%-42s %12s %12s %10s\n", "difference", "rms", "max", "rel_rms")
for (label, D) in diffs
    m = colnorms(D)
    @printf("%-42s %12.4e %12.4e %10.4f\n", label, rms(m), maximum(m), rms(m) / omega_ref)
end

# ---------------- radial profile CSV -------------------------------------------
# rotor rotates about x; radial coordinate r = sqrt(y^2 + z^2)

r = [sqrt(cps[2, i]^2 + cps[3, i]^2) for i in 1:ncp]
edges = range(0, 1.05 * ROTOR_R; length=NBINS + 1)

profile_fields = [
    ("omega_basis", omega_basis),
    ("omega_extra_direct", omega_direct_kernel),
    ("omega_curl_direct", omega_curl_direct),
    ("d_extra_basis", omega_direct_kernel .- omega_basis),
    ("d_extra_curl", omega_direct_kernel .- omega_curl_direct),
    ("d_curl_basis", omega_curl_direct .- omega_basis),
    ("d_fmm_vs_direct", omega_fmm .- omega_direct_kernel),
]

csv_path = joinpath(OUTDIR, "T4_particle_vorticity_radial_profile_step$(STEP).csv")
open(csv_path, "w") do io
    println(io, "r_mid,count," * join(("rms_" * l for (l, _) in profile_fields), ","))
    for b in 1:NBINS
        idx = [i for i in 1:ncp if edges[b] <= r[i] < edges[b + 1]]
        isempty(idx) && continue
        vals = [rms(colnorms(A[:, idx])) for (_, A) in profile_fields]
        println(io, join(vcat([string((edges[b] + edges[b + 1]) / 2), string(length(idx))],
                              [@sprintf("%.6e", v) for v in vals]), ","))
    end
end
println("\nWrote $(csv_path)")
