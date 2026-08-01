#=##############################################################################
# code_audit Task 1b — Fix-candidate quantification: the user's alternative
# affine-Kutta formulation (docs/wake_solve_schemes.md, "Suggested alternative
# formulation for semiinfinite_wake=false"), quantified on the same frozen
# settled marched state as task1b_numerics.jl.
#
# Reuses ALL machinery from task1b_numerics.jl via include() — the base script
# rebuilds the settled state, assembles G_B/G_A/K_int, probes phi_tr, and
# computes the A/A'/B solves and the e4a residual. This script then:
#
#  1. Three-way + variants solve comparison on the identical frozen state:
#       A   (semiinfinite)                     — from base script
#       B   (production velocity-through-σ)    — marched μ*; also c=0 sanity
#       A'  (explicit potential in RHS)        — from base script
#       ALT (G_B μ̃ = −S(σ0+σ) + W c, γ = Cμ̃ − c) with
#           c = exact (C·φ_tr), trapezoid, Simpson
#  2. Analytic-prediction test: γ(c-exact) − Δμ_A' =?= −C·G_B⁻¹·e4a
#  3. Doc check 2 in isolation: gauge-optimized discrete Green identity
#     residual (direct e4a measurement, no solve).
#
# Outputs: code_audit/results/task1b/{task1b_altform_summary.txt,
#          task1b_altform_spanwise.csv}  (base-script outputs also rewritten)
# Run:  julia --project -t auto code_audit/scripts/task1b_alternative_formulation.jl
#       TASK1B_SMOKE=1 for a coarse/fast smoke run.
=###############################################################################

include(joinpath(@__DIR__, "task1b_numerics.jl"))

println("\n" * "="^78)
println("Task 1b extension — alternative affine-Kutta formulation")

alt_io = IOBuffer()
macro alogf(fmt, args...)
    esc(quote
        @printf($fmt, $(args...))
        @printf(alt_io, $fmt, $(args...))
    end)
end

@alogf("Task 1b alternative-formulation quantification — %s run\n", SMOKE ? "SMOKE" : "FULL")

# ----------------------------------------------------------------- [E1] setup
@assert length(wing.shedding) == 1 "script assumes a single shedding group"
@assert length(unique(upper)) == nst && isempty(intersect(Set(upper), Set(lower))) "W·c construction requires unique upper panels disjoint from lower"

"Free-wake velocity at arbitrary positions (direct backend)."
function probe_velocity(positions, sources, backend)
    probes = FastMultipole.ProbeSystem(length(positions), Float64)
    for (i, x) in enumerate(positions)
        probes.position[i] = SVector{3}(x[1], x[2], x[3])
        probes.gradient[i] = zero(SVector{3,Float64})
        probes.scalar_potential[i] = 0.0
    end
    pnl.influence!((probes,), sources, backend;
        precalc=false, scalar_potential=false, velocity=true)
    return [probes.gradient[i] for i in eachindex(positions)]
end

xus = [SVector{3}(wing.controlpoints[:, upper[i]]) for i in 1:nst]
xls = [SVector{3}(wing.controlpoints[:, lower[i]]) for i in 1:nst]
xms = [0.5 .* (xus[i] .+ xls[i]) for i in 1:nst]

# --------------------------- [E2] geometry check: straight chord inside body?
# Interior test at the chord midpoint: it must lie on the interior side of both
# TE panels' planes (outward normals). Doc warns straight chord validity is not
# automatic for nonconvex bodies.
s_int = zeros(nst, 2)
for i in 1:nst
    nu = SVector{3}(wing.normals[:, upper[i]])
    nl = SVector{3}(wing.normals[:, lower[i]])
    s_int[i, 1] = dot(xms[i] .- xus[i], nu)
    s_int[i, 2] = dot(xms[i] .- xls[i], nl)
end
@alogf("\n[E2] straight-chord interiority (midpoint side of TE panel planes; interior <= 0):\n")
@alogf("     max over stations: upper-side %.3e, lower-side %.3e; |d| midspan = %.4g m\n",
    maximum(s_int[:, 1]), maximum(s_int[:, 2]), norm(xus[imid] - xls[imid]))
if maximum(s_int) > 0
    @alogf("     WARNING: some chord midpoints lie outside a TE panel plane — straight chord suspect there.\n")
else
    @alogf("     All chord midpoints interior: straight centroid-to-centroid chord valid.\n")
end

# ------------------------------- [E3] Kutta-trace vector c: exact / trap / Simpson
c_exact = copy(Cphi)                                   # C * phi_tr (ground truth)
u_up  = probe_velocity(xus, wake_srcs, backend_dir)
u_lo  = probe_velocity(xls, wake_srcs, backend_dir)
u_mid = probe_velocity(xms, wake_srcs, backend_dir)

dvec = [xus[i] .- xls[i] for i in 1:nst]
c_trap_raw = [0.5 * dot(u_up[i] .+ u_lo[i], dvec[i]) for i in 1:nst]
c_simp_raw = [dot(u_lo[i] .+ 4 .* u_mid[i] .+ u_up[i], dvec[i]) / 6 for i in 1:nst]

# sign convention: c_e = ∫ u_f·dl from lower to upper must reproduce q_up − q_low
# with the code's potential/velocity sign pairing. Detect empirically.
sgn = cor(c_trap_raw, c_exact) >= 0 ? 1.0 : -1.0
c_trap = sgn .* c_trap_raw
c_simp = sgn .* c_simp_raw
@alogf("\n[E3] Kutta-trace vector c per shedding edge (nst=%d):\n", nst)
@alogf("     sign pairing u_f vs probed potential: %s (cor(c_trap_raw, C·phi_tr) = %+.4f)\n",
    sgn > 0 ? "u_f = +grad(phi)" : "u_f = -grad(phi) — FLIPPED", cor(c_trap_raw, c_exact))
@alogf("     station means: c_exact=%.6g  c_trap=%.6g  c_simp=%.6g\n",
    mean(c_exact), mean(c_trap), mean(c_simp))
relerr(a, b) = maximum(abs.(a .- b)) / maximum(abs.(b))
@alogf("     quadrature error vs exact: trap max-rel %.3e (rms %.3e), Simpson max-rel %.3e (rms %.3e)\n",
    relerr(c_trap, c_exact), rms(c_trap .- c_exact) / rms(c_exact),
    relerr(c_simp, c_exact), rms(c_simp .- c_exact) / rms(c_exact))
@alogf("     mid-span: c_exact=%.6g c_trap=%.6g c_simp=%.6g (dmu_A=%.6g)\n",
    c_exact[imid], c_trap[imid], c_simp[imid], dmu_A[imid])

# --------------------------------------------- [E4] W·c and the ALT solves
# W c = (G_B − K_int)·v with C·v = c  (v[upper_i] = c_i, else 0): the attached-
# wake coupling columns KwC_B factor as W∘C by construction.
function Wc_of(cvec)
    v = zeros(n)
    for i in 1:nst
        v[upper[i]] = cvec[i]
    end
    @assert maximum(abs.(dmu(v) .- cvec)) < 1e-14 * max(1.0, maximum(abs.(cvec)))
    return KwC_B * v
end

rhs_base = -phi_sStar    # = −S(σ0+σ), production RHS (free wake through σ)
function alt_solve(cvec)
    mu_t = luGB \ (rhs_base .+ Wc_of(cvec))
    return dmu(mu_t) .- cvec        # affine Kutta: γ = C μ̃ − c
end

gamma_exact = alt_solve(c_exact)
gamma_trap  = alt_solve(c_trap)
gamma_simp  = alt_solve(c_simp)
gamma_c0    = alt_solve(zeros(nst))   # sanity: must reproduce scheme B (μ*)

ratio(g)  = sum(g) / sum(dmu_A)
ratiom(g) = g[imid] / dmu_A[imid]
ratio_B   = ratio(dmu_star)
recov(g)  = (ratio(g) - ratio_B) / (1 - ratio_B)   # fraction of settled deficit recovered

@alogf("\n[E4] Three-way + variants, identical frozen geometry+wake (ratios vs scheme A):\n")
@alogf("     %-34s %10s %10s %14s\n", "scheme", "sum-ratio", "mid-span", "deficit-recov")
@alogf("     %-34s %10.4f %10.4f %14s\n", "A  (semiinfinite, frozen geom)", 1.0, 1.0, "--")
@alogf("     %-34s %10.4f %10.4f %14.1f%%\n", "B  (marched, velocity-through-σ)",
    ratio_B, ratiom(dmu_star), 0.0)
@alogf("     %-34s %10.4f %10.4f %14.1f%%\n", "A' (explicit potential, one-shot)",
    ratio(dmu_Ap), ratiom(dmu_Ap), 100 * recov(dmu_Ap))
@alogf("     %-34s %10.4f %10.4f %14.1f%%\n", "ALT c-exact (affine Kutta)",
    ratio(gamma_exact), ratiom(gamma_exact), 100 * recov(gamma_exact))
@alogf("     %-34s %10.4f %10.4f %14.1f%%\n", "ALT c-trapezoid",
    ratio(gamma_trap), ratiom(gamma_trap), 100 * recov(gamma_trap))
@alogf("     %-34s %10.4f %10.4f %14.1f%%\n", "ALT c-Simpson",
    ratio(gamma_simp), ratiom(gamma_simp), 100 * recov(gamma_simp))
@alogf("     %-34s %10.4f %10.4f %14s\n", "ALT c=0 (sanity -> B)",
    ratio(gamma_c0), ratiom(gamma_c0), "--")
@alogf("     c=0 sanity: max|γ(c=0) − Δμ*| = %.3e (must be ~solver eps)\n",
    maximum(abs.(gamma_c0 .- dmu_star)))
@alogf("     implied CL (linear scaling of settled CL=%.6g with sum γ): ", CL_settled)
@alogf("A'=%.5g  ALT-exact=%.5g  ALT-trap=%.5g\n",
    CL_settled * sum(dmu_Ap) / sum(dmu_star),
    CL_settled * sum(gamma_exact) / sum(dmu_star),
    CL_settled * sum(gamma_trap) / sum(dmu_star))

# ---------------------- [E5] analytic prediction: ALT(c-exact) − A' = −C·G_B⁻¹·e4a
pred_e4a = -Cimg(luGB \ e4a)
diff_meas = gamma_exact .- dmu_Ap
@alogf("\n[E5] Analytic prediction test (orchestrator): γ(c-exact) − Δμ_A' =?= −C·G_B⁻¹·e4a\n")
@alogf("     station means: measured %.6g, predicted %.6g; max|meas−pred| = %.3e (rel %.3e)\n",
    mean(diff_meas), mean(pred_e4a),
    maximum(abs.(diff_meas .- pred_e4a)),
    maximum(abs.(diff_meas .- pred_e4a)) / max(maximum(abs.(pred_e4a)), 1e-300))
@alogf("     retained-e4a term as fraction of measured deficit D1: mean = %+.1f%% (G_A-route E2b was %+.1f%%)\n",
    100 * mean(-pred_e4a) / mean(D1), 100 * mean(E2b) / mean(D1))
@alogf("     P5-elimination check: γ(c-exact) removes trace-shift exactly by construction;\n")
@alogf("     measured recovery ALT-exact vs A': %.1f%% vs %.1f%% of deficit (gap = retained e4a).\n",
    100 * recov(gamma_exact), 100 * recov(dmu_Ap))

# ---------------------- [E6] doc check 2 in isolation: gauge-optimized Green identity
# e4a = φ_σ,fw − (I−K_int)φ_tr from the base script. Constant-gauge freedom of
# φ_tr changes e4a by −α(1 − K_int·1); optimize α (should be tiny on a closed body).
g1 = ones(n) .- Kint * ones(n)
@alogf("\n[E6] Doc check 2 (discrete Green identity, no Kutta columns, direct):\n")
@alogf("     ||K_int·1 − 1||_inf = %.3e (closed-body constant mode)\n", maximum(abs.(g1)))
if norm(g1) / sqrt(n) > 1e-10
    alpha = dot(e4a, g1) / dot(g1, g1)
    e4a_gauged = e4a .- alpha .* g1
    @alogf("     ||e4a||/||φ_σ,fw||: raw = %.4g, gauge-optimized (α=%.3g) = %.4g\n",
        norm(e4a) / norm(phi_sfw), alpha, norm(e4a_gauged) / norm(phi_sfw))
else
    @alogf("     ||e4a||/||φ_σ,fw|| = %.4g (gauge-invariant: (I−K_int)·1 = 0 to machine precision)\n",
        norm(e4a) / norm(phi_sfw))
end
@alogf("     -> irreducible P1/P2 sampling residual; ALT formulation retains it, A' removes it.\n")

# --------------------------------------------------------------------- outputs
open(joinpath(RESULTS, "task1b_altform_spanwise.csv"), "w") do io
    println(io, "eta,y,dmu_A,dmu_marched,dmu_Aprime,gamma_alt_exact,gamma_alt_trap,gamma_alt_simp,c_exact,c_trap,c_simp,pred_negCGBinv_e4a,measured_altexact_minus_Aprime")
    for i in 1:nst
        println(io, join((eta[i], ysta[i], dmu_A[i], dmu_star[i], dmu_Ap[i],
            gamma_exact[i], gamma_trap[i], gamma_simp[i],
            c_exact[i], c_trap[i], c_simp[i], pred_e4a[i], diff_meas[i]), ","))
    end
end
open(joinpath(RESULTS, "task1b_altform_summary.txt"), "w") do io
    write(io, String(take!(alt_io)))
end
println("\nWrote alternative-formulation outputs under $(RESULTS)")
