# Numerical check: curl-of-J particle vorticity vs scalar zeta-basis vorticity.
#
# BRAINSTORM item 010. A single Gaussian-erf vortex particle is enough to show
# that curl(u) is the Leray-projected vorticity, not the scalar Gaussian
# zeta-basis field: at the particle center, curl(u) = (2/3) zeta(0) Gamma.
#
# Run: julia --project=test examples/particle_vorticity_curl_vs_basis_check.jl

import ForwardDiff
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "particle_gauserf_hessian_helpers.jl"))
using .ParticleGauserfHessian:
    HAVE_SF, ZETA0, zeta_sigma, green_profile, phi_gauserf,
    hess_green_times_gamma

function green_sigma(x, X, sigma)
    r = x .- X
    rho2 = dot(r, r) / sigma^2
    return green_profile(rho2) / (4pi * sigma)
end

function induced_velocity(x, X, Gamma, sigma)
    r = x .- X
    rho2 = dot(r, r) / sigma^2
    f = phi_gauserf(rho2) / sigma^3
    cross = (r[2] * Gamma[3] - r[3] * Gamma[2],
             r[3] * Gamma[1] - r[1] * Gamma[3],
             r[1] * Gamma[2] - r[2] * Gamma[1])
    return -(1 / (4pi)) * f .* collect(cross)
end

function curl_from_velocity(x, X, Gamma, sigma)
    J = ForwardDiff.jacobian(y -> induced_velocity(y, X, Gamma, sigma), x)
    return [J[3, 2] - J[2, 3], J[1, 3] - J[3, 1], J[2, 1] - J[1, 2]]
end

function omega_zeta(x, X, Gamma, sigma)
    r = x .- X
    rho2 = dot(r, r) / sigma^2
    return zeta_sigma(rho2, sigma) .* collect(Gamma)
end

function delta_analytic(x, X, Gamma, sigma)
    H = ForwardDiff.hessian(y -> green_sigma(y, X, sigma), x)
    return H * collect(Gamma)
end

# M1 diagnostic: the true solenoidal vorticity omega_J = omega_zeta +
# sum_i Hess(G_sigma_i)(x-X_i) Gamma_i, carrying the strengths Gamma unchanged.
# Brute-force O(N^2) loop (matches particle_divergence_free_check.jl); the
# production path replaces this with the existing hard-cutoff neighbor list /
# FMM Hessian sum (see examples/particle_overlap_residual.jl).
function omega_J_corrected(x, Xs, Gammas, sigmas)
    wj = zeros(eltype(x), 3)
    @inbounds for p in eachindex(sigmas)
        r = x .- Xs[p]
        rho2 = (r[1]^2 + r[2]^2 + r[3]^2) / sigmas[p]^2
        wj = wj .+ zeta_sigma(rho2, sigmas[p]) .* collect(Gammas[p])
        wj = wj .+ hess_green_times_gamma(r, sigmas[p], Gammas[p])
    end
    return wj
end

function omegaj_divergence(x, X, Gamma, sigma)
    G = ForwardDiff.jacobian(y -> curl_from_velocity(y, X, Gamma, sigma), x)
    return tr(G), norm(G)
end

# --- M2: effective strength Gamma_eff best-fitting omega_J -------------------
# Carry an effective strength whose scalar zeta-basis sum best-fits the
# solenoidal omega_J (the M1 field). Collocate at the particle centers:
#   M_{ki} = zeta_{sigma_i}(X_k - X_i)               (scalar blocks),
#   (H Gamma)_k = sum_i Hess(G_sigma_i)(X_k - X_i) Gamma_i,
# then Gamma_eff solves  M Gamma_eff = (M + H) Gamma, i.e.
#   Gamma_eff = Gamma + M^{-1} H Gamma.
# These dense brute-force assemblies are the offline analogue of FLOWVPM's
# rbf_conjugategradient core-reset (FLOWVPM_viscous.jl:314); production swaps in
# the neighbor-list / FMM sums from examples/particle_overlap_residual.jl.

# M_{ki} = zeta_{sigma_i}(X_k - X_i) (note source sigma_i; symmetric iff all
# sigma are equal, in which case M is also SPD and CG applies).
function overlap_matrix(Xs, sigmas)
    n = length(sigmas)
    M = zeros(n, n)
    @inbounds for k in 1:n, i in 1:n
        r = Xs[k] .- Xs[i]
        rho2 = (r[1]^2 + r[2]^2 + r[3]^2) / sigmas[i]^2
        M[k, i] = zeta_sigma(rho2, sigmas[i])
    end
    return M
end

# omega_J at the particle centers as a 3 x n matrix: (M + H) Gamma collocated.
function omega_J_nodes(Xs, Gammas, sigmas)
    n = length(sigmas)
    M = overlap_matrix(Xs, sigmas)
    Gmat = reduce(hcat, collect.(Gammas))       # 3 x n, column i = Gamma_i
    W = Gmat * transpose(M)                      # 3 x n, (M Gamma)_k = omega_zeta(X_k)
    @inbounds for k in 1:n
        acc = zeros(3)
        for i in 1:n
            acc = acc .+ hess_green_times_gamma(Xs[k] .- Xs[i], sigmas[i], Gammas[i])
        end
        W[:, k] = W[:, k] .+ acc                 # add (H Gamma)_k
    end
    return W, M
end

# Solve M Gamma_eff = (M+H) Gamma by a dense direct solve (robust for any sigma).
function gamma_eff(Xs, Gammas, sigmas)
    W, M = omega_J_nodes(Xs, Gammas, sigmas)     # 3 x n RHS = omega_J at centers
    Geff_T = M \ transpose(W)                    # n x 3
    return permutedims(Geff_T), M, W
end

# CG variant mirroring rbf_conjugategradient; valid when M is SPD (equal sigma).
function gamma_eff_cg(Xs, Gammas, sigmas; tol = 1e-13, itmax = 5000)
    W, M = omega_J_nodes(Xs, Gammas, sigmas)
    n = length(sigmas)
    Geff = zeros(3, n)
    for a in 1:3
        b = W[a, :]
        x = copy(b)
        r = b .- M * x
        p = copy(r)
        rs = dot(r, r)
        for _ in 1:itmax
            sqrt(rs) <= tol && break
            Ap = M * p
            alpha = rs / dot(p, Ap)
            x .+= alpha .* p
            r .-= alpha .* Ap
            rsnew = dot(r, r)
            p = r .+ (rsnew / rs) .* p
            rs = rsnew
        end
        Geff[a, :] = x
    end
    return Geff
end

# Field reconstructed from a strength set: sum_i zeta_{sigma_i}(x - X_i) S_i.
function zeta_field(x, Xs, Strengths, sigmas)
    w = zeros(eltype(x), 3)
    @inbounds for i in eachindex(sigmas)
        r = x .- Xs[i]
        rho2 = (r[1]^2 + r[2]^2 + r[3]^2) / sigmas[i]^2
        w = w .+ zeta_sigma(rho2, sigmas[i]) .* Strengths[i]
    end
    return w
end

# M4 projected-Gaussian basis field:
#   K_sigma Gamma = zeta_sigma Gamma + Hess(G_sigma) Gamma.
# This is the Leray-projected blob, i.e. the represented vorticity omega_J.
function projected_gaussian_field(x, Xs, Gammas, sigmas)
    return omega_J_corrected(x, Xs, Gammas, sigmas)
end

# Local matrix-valued divergence-free Gaussian/RBF kernel:
#   Kloc = Hess(phi) - Delta(phi) I, with phi = zeta_sigma.
# It is exactly divergence-free columnwise and local/L1, hence its whole-space
# integral must vanish by the M12 Fourier no-go theorem.
function local_divfree_gaussian_field(x, X, Gamma, sigma)
    r = x .- X
    r2 = dot(r, r)
    phi = zeta_sigma(r2 / sigma^2, sigma)
    yg = dot(r, Gamma)
    return phi .* ((2 / sigma^2 - r2 / sigma^4) .* Gamma .+
                   (yg / sigma^4) .* r)
end

# --- M6: projected-Gaussian overlap evolution (block-kernel collocation) -----
# Promote the M4 projected blob to the carried model: omega = sum_i K_{sigma_i} Gamma_i,
# K_sigma = zeta_sigma I + Hess(G_sigma) (= omega_J). The strength equation is the
# block-kernel collocation analogue of item 008:
#   sum_i K_{sigma_i}(X_k - X_i) Gdot_i = J(X_k) omega(X_k) - convection - viscous,
# i.e. a 3np x 3np system with 3x3 blocks K_{sigma_i}(X_k - X_i) instead of the
# scalar M_{ki} I. See the M6 section of BRAINSTORM/010.

# The 3x3 kernel block K_sigma(y) = zeta_sigma(|y|^2/sigma^2) I + Hess(G_sigma)(y).
# Its columns are hess_green_times_gamma applied to the basis vectors (the Hessian
# is linear in Gamma), plus the isotropic zeta term.
function kblock(y, sigma)
    rho2 = (y[1]^2 + y[2]^2 + y[3]^2) / sigma^2
    z = zeta_sigma(rho2, sigma)
    H = hcat(hess_green_times_gamma(y, sigma, [1.0, 0.0, 0.0]),
             hess_green_times_gamma(y, sigma, [0.0, 1.0, 0.0]),
             hess_green_times_gamma(y, sigma, [0.0, 0.0, 1.0]))
    return z * Matrix(I, 3, 3) + H
end

# Assemble the dense 3n x 3n block collocation operator with blocks
# 𝒦[k,i] = K_{sigma_i}(X_k - X_i) (note the SOURCE sigma_i, as in overlap_matrix).
function assemble_block_K(Xs, sigmas)
    n = length(sigmas)
    K = zeros(3n, 3n)
    @inbounds for k in 1:n, i in 1:n
        K[3k-2:3k, 3i-2:3i] = kblock(Xs[k] .- Xs[i], sigmas[i])
    end
    return K
end

function check_point(x, X, Gamma, sigma)
    wz = omega_zeta(x, X, Gamma, sigma)
    wj = curl_from_velocity(x, X, Gamma, sigma)
    delta = wj - wz
    da = delta_analytic(x, X, Gamma, sigma)
    divwj, gradwj_norm = omegaj_divergence(x, X, Gamma, sigma)

    delta_relerr = norm(delta - da) / max(norm(da), eps())
    div_rel = abs(divwj) / max(gradwj_norm, eps())
    discrepancy = norm(delta) / max(norm(wz), eps())
    two_thirds_err = norm(wj - (2 / 3) .* wz) / max(norm(wz), eps())
    return (; wz, wj, delta, da, divwj, gradwj_norm, delta_relerr, div_rel,
            discrepancy, two_thirds_err)
end

function main()
    sigma = 1.0
    X = [0.0, 0.0, 0.0]
    Gamma = [1.0, 0.2, 0.0]
    direction = normalize([1.0, 2.0, -0.5])
    radii = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0]

    println("Particle vorticity: curl-of-J vs zeta-basis check")
    println("erf source: ", HAVE_SF ? "SpecialFunctions.erf" : "self-contained series")
    println("sigma = $sigma, Gamma = $Gamma")
    println()

    @printf("%8s  %14s  %14s  %14s  %14s\n",
            "r/sigma", "||d||/||wz||", "delta relerr", "rel div wJ", "2/3 center err")
    println("-"^78)

    rows = String[]
    max_delta_relerr = 0.0
    max_div_rel = 0.0
    center_discrepancy = NaN
    center_two_thirds_err = NaN

    for radius in radii
        x = X .+ sigma * radius .* direction
        r = check_point(x, X, Gamma, sigma)
        max_delta_relerr = max(max_delta_relerr, r.delta_relerr)
        max_div_rel = max(max_div_rel, r.div_rel)
        if radius == 0.0
            center_discrepancy = r.discrepancy
            center_two_thirds_err = r.two_thirds_err
        end
        @printf("%8.3f  %14.6e  %14.6e  %14.6e  %14.6e\n",
                radius, r.discrepancy, r.delta_relerr, r.div_rel, r.two_thirds_err)
        push!(rows, @sprintf("%.6g,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e",
                             radius, r.discrepancy, r.delta_relerr, r.div_rel,
                             r.two_thirds_err, norm(r.wz), norm(r.wj)))
    end

    self_ok = abs(center_discrepancy - 1 / 3) < 1e-9 && center_two_thirds_err < 1e-9
    structure_ok = max_delta_relerr < 1e-9
    solenoidal_ok = max_div_rel < 1e-8
    nonvacuous_ok = all(norm(check_point(X .+ sigma * radius .* direction, X, Gamma, sigma).wz) > 0 for radius in radii[1:5])
    final = self_ok && structure_ok && solenoidal_ok && nonvacuous_ok

    outdir = joinpath(@__DIR__, "..", "data", "particle_vorticity_curl_vs_basis")
    mkpath(outdir)
    outfile = joinpath(outdir, "curl_vs_basis_report.csv")
    open(outfile, "w") do io
        println(io, "r_over_sigma,discrepancy_rel,delta_analytic_relerr,omegaJ_div_rel,two_thirds_center_err,norm_omega_zeta,norm_omega_J")
        for row in rows
            println(io, row)
        end
    end

    println()
    @printf("center discrepancy: %.16e (target 1/3 = %.16e)\n", center_discrepancy, 1 / 3)
    @printf("max delta analytic relative error: %.3e\n", max_delta_relerr)
    @printf("max relative div(omega_J): %.3e\n", max_div_rel)
    println(final ?
        "RESULT: PASS -- curl-of-J is the projected vorticity, with the 2/3 self-term at the center." :
        "RESULT: FAIL -- investigate the kernel/sign/Taylor implementation.")
    println("Wrote ", outfile)
    return final
end

# Opt-in closed-form self-check (M1): assert the analytic
# hess_green_times_gamma matches the ForwardDiff delta_analytic along the radial
# line. Not part of main()'s default gates. Run with M1_SELFTEST=1.
function m1_selftest()
    sigma = 1.0
    X = [0.0, 0.0, 0.0]
    Gamma = [1.0, 0.2, 0.0]
    direction = normalize([1.0, 2.0, -0.5])
    radii = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0]
    maxerr = 0.0
    for radius in radii
        x = X .+ sigma * radius .* direction
        analytic = hess_green_times_gamma(x .- X, sigma, Gamma)
        ad = delta_analytic(x, X, Gamma, sigma)
        err = norm(analytic - ad) / max(norm(ad), eps())
        maxerr = max(maxerr, err)
        @printf("  r/sigma=%6.3f  closed-form vs ForwardDiff relerr = %.3e\n", radius, err)
    end
    ok = maxerr < 1e-9
    @printf("M1 self-test max relerr = %.3e -- %s\n", maxerr, ok ? "PASS" : "FAIL")
    return ok
end

# Opt-in M2 projection self-check: Gamma_eff = Gamma + M^{-1} H Gamma. Not part
# of main()'s default gates. Run with M2_SELFTEST=1.
function m2_selftest()
    ok = true

    # (1) Single particle: the self-only projection collapses to Gamma_eff=(2/3)Gamma.
    Xs1 = [[0.0, 0.0, 0.0]]
    G1 = [[1.0, 0.2, 0.0]]
    s1 = [1.0]
    Geff1, _, _ = gamma_eff(Xs1, G1, s1)
    err_23 = norm(Geff1[:, 1] - (2 / 3) .* G1[1]) / norm(G1[1])
    @printf("  (1) single particle: ||Geff - (2/3)Gamma||/||Gamma|| = %.3e\n", err_23)
    ok &= err_23 < 1e-12

    # (2) Two overlapping particles (equal sigma => M is SPD): dense solve and CG
    #     agree, and the linear residual ||M Geff - omega_J_nodes|| is roundoff.
    Xs2 = [[0.0, 0.0, 0.0], [0.8, 0.3, -0.2]]
    G2 = [[1.0, 0.2, 0.0], [-0.4, 0.7, 0.5]]
    s2 = [1.0, 1.0]
    Geff2, M2, W2 = gamma_eff(Xs2, G2, s2)
    Geff2_cg = gamma_eff_cg(Xs2, G2, s2)
    cg_err = norm(Geff2 - Geff2_cg) / norm(Geff2)
    lin_res = norm(M2 * transpose(Geff2) - transpose(W2)) / norm(W2)
    @printf("  (2) two particles: CG vs dense = %.3e, linear residual = %.3e\n", cg_err, lin_res)
    ok &= cg_err < 1e-9 && lin_res < 1e-10

    # (3) Best-fit improvement: off-node, the Gamma_eff field is closer to omega_J
    #     than the raw omega_zeta field, but a nonzero floor remains (a scalar
    #     blob is never exactly solenoidal).
    Geff_cols = [Geff2[:, i] for i in 1:length(s2)]
    probes = [[0.4, 0.15, -0.1], [0.2, -0.3, 0.25], [1.1, 0.6, 0.1],
              [-0.3, 0.2, 0.4], [0.6, 0.0, -0.5]]
    res_raw = 0.0
    res_eff = 0.0
    den = 0.0
    for x in probes
        wj = omega_J_corrected(x, Xs2, G2, s2)
        wraw = zeta_field(x, Xs2, G2, s2)
        weff = zeta_field(x, Xs2, Geff_cols, s2)
        res_raw += norm(wraw - wj)^2
        res_eff += norm(weff - wj)^2
        den += norm(wj)^2
    end
    res_raw = sqrt(res_raw / den)
    res_eff = sqrt(res_eff / den)
    @printf("  (3) off-node field vs omega_J: raw Gamma = %.3e, Gamma_eff = %.3e (floor)\n",
            res_raw, res_eff)
    ok &= res_eff < res_raw

    # (4) Divergence reduction: the whole point of M2 is to best-fit a solenoidal
    #     field, so div(Gamma_eff field) should be smaller than div(raw field) --
    #     but nonzero (a scalar Gaussian blob is never exactly divergence-free).
    divfield(x, S) = tr(ForwardDiff.jacobian(y -> zeta_field(y, Xs2, S, s2), x))
    fieldnorm(x, S) = norm(ForwardDiff.jacobian(y -> zeta_field(y, Xs2, S, s2), x))
    div_raw = 0.0
    div_eff = 0.0
    gden = 0.0
    for x in probes
        div_raw += divfield(x, G2)^2
        div_eff += divfield(x, Geff_cols)^2
        gden += fieldnorm(x, G2)^2
    end
    div_raw = sqrt(div_raw / gden)
    div_eff = sqrt(div_eff / gden)
    @printf("  (4) off-node rel divergence: raw Gamma = %.3e, Gamma_eff = %.3e (floor)\n",
            div_raw, div_eff)
    ok &= div_eff < div_raw

    @printf("M2 self-test -- %s\n", ok ? "PASS" : "FAIL")
    return ok
end

# Opt-in M4 representation self-check: K_sigma = zeta_sigma I + Hess(G_sigma)
# is the Leray-projected Gaussian blob. This verifies representation facts only;
# it does not derive or wire a projected-basis strength update. Run with
# M4_SELFTEST=1.
function m4_selftest()
    ok = true
    sigma = 1.0
    X = [0.0, 0.0, 0.0]
    Gamma = [1.0, 0.2, 0.0]
    Xs = [X]
    Gammas = [Gamma]
    sigmas = [sigma]

    # (1) Center factor: K_sigma(0) Gamma = (2/3) zeta_sigma(0) Gamma.
    center = projected_gaussian_field(X, Xs, Gammas, sigmas)
    target_center = (2 / 3) * zeta_sigma(0.0, sigma) .* Gamma
    center_err = norm(center - target_center) / norm(target_center)
    @printf("  (1) K_sigma center factor relerr = %.3e\n", center_err)
    ok &= center_err < 1e-12

    # (2) K_sigma Gamma equals curl(u) from ForwardDiff through the velocity.
    probes = [[0.0, 0.0, 0.0], [0.125, 0.25, -0.0625], [0.5, -0.25, 0.2],
              [1.0, 0.2, -0.4], [2.0, -0.5, 0.1], [4.0, 1.0, -0.5]]
    max_curl_err = 0.0
    for x in probes
        wk = projected_gaussian_field(x, Xs, Gammas, sigmas)
        wj = curl_from_velocity(x, X, Gamma, sigma)
        max_curl_err = max(max_curl_err, norm(wk - wj) / max(norm(wj), eps()))
    end
    @printf("  (2) K_sigma Gamma vs curl(u) max relerr = %.3e\n", max_curl_err)
    ok &= max_curl_err < 1e-9

    # (3) Solenoidality: div(K_sigma Gamma) should be roundoff-small.
    max_div_rel = 0.0
    for x in probes
        J = ForwardDiff.jacobian(y -> projected_gaussian_field(y, Xs, Gammas, sigmas), x)
        max_div_rel = max(max_div_rel, abs(tr(J)) / max(norm(J), eps()))
    end
    @printf("  (3) relative div(K_sigma Gamma) max = %.3e\n", max_div_rel)
    ok &= max_div_rel < 1e-8

    # (4) Whole-space integral caveat: the Hessian tail is algebraic, and the
    # naive symmetric-volume integral of the pointwise K_sigma field tends toward
    # (2/3)Gamma, not Gamma. This catches the important M4 nuance that the
    # projected blob has a nonlocal k=0 / far-field convention caveat.
    n = 25
    L = 6.0 * sigma
    dx = 2L / n
    integral = zeros(3)
    for i in 1:n, j in 1:n, k in 1:n
        x = [-L + (i - 0.5) * dx, -L + (j - 0.5) * dx, -L + (k - 0.5) * dx]
        integral .+= projected_gaussian_field(x, Xs, Gammas, sigmas) .* dx^3
    end
    circulation_err = norm(integral - (2 / 3) .* Gamma) / norm(Gamma)
    @printf("  (4) finite-cube integral vs (2/3)Gamma relerr = %.3e\n", circulation_err)
    ok &= circulation_err < 3e-3

    @printf("M4 self-test -- %s\n", ok ? "PASS" : "FAIL")
    return ok
end

# Opt-in M6 projected-basis strength-evolution self-check. Verifies the three
# decisive M6 facts: the no-overlap reduction Gdot=J Gamma, a consistent two-particle
# block solve (with conditioning vs the scalar M), and that the core-spreading term
# d/dsigma(K_sigma Gamma) is finite/AD-consistent. Run with M6_SELFTEST=1.
function m6_selftest()
    ok = true
    sigma = 1.0
    Jrand = [0.30 -0.20 0.10; 0.05 -0.40 0.25; -0.15 0.20 0.35]

    # (1) No-overlap reduction. K(0) = (2/3) zeta(0) I, so the collocation equation
    #     K(0) Gdot = J (K(0) Gamma) collapses to Gdot = J Gamma to round-off.
    z0 = zeta_sigma(0.0, sigma)
    K0 = kblock([0.0, 0.0, 0.0], sigma)
    target0 = (2 / 3) * z0 * Matrix(I, 3, 3)
    center_err = norm(K0 - target0) / norm(target0)
    @printf("  (1) K(0) vs (2/3)zeta(0) I relerr = %.3e\n", center_err)
    ok &= center_err < 1e-12

    Gamma = [1.0, 0.2, -0.1]
    omega0 = K0 * Gamma                      # = (2/3) zeta(0) Gamma
    Gdot0 = K0 \ (Jrand * omega0)            # classic stretching (w.grad)u = J w
    red_err = norm(Gdot0 - Jrand * Gamma) / norm(Jrand * Gamma)
    @printf("  (1) no-overlap Gdot vs J Gamma relerr = %.3e\n", red_err)
    ok &= red_err < 1e-12

    # (2) Two-particle block solve. Assemble the 6x6 block operator, build the RHS
    #     from the analytic stretching J_k omega(X_k) with omega = sum_i K_i Gamma_i,
    #     solve, and report block-system conditioning vs the scalar overlap M.
    Xs = [[0.0, 0.0, 0.0], [0.8, 0.3, -0.2]]
    Gammas = [[1.0, 0.2, 0.0], [-0.4, 0.7, 0.5]]
    sigmas = [1.0, 1.0]
    Js = [Jrand, [-0.10 0.20 0.05; 0.30 0.15 -0.20; 0.10 -0.05 0.25]]

    # sanity: kblock * Gamma matches the single-particle projected field.
    kb_err = 0.0
    for k in 1:2
        kb_err = max(kb_err, norm(kblock(Xs[1] .- Xs[k], sigmas[k]) * Gammas[k] -
                                  projected_gaussian_field(Xs[1], [Xs[k]], [Gammas[k]], [sigmas[k]])))
    end
    ok &= kb_err < 1e-12

    Kblk = assemble_block_K(Xs, sigmas)
    b = zeros(6)
    for k in 1:2
        wk = projected_gaussian_field(Xs[k], Xs, Gammas, sigmas)
        b[3k-2:3k] = Js[k] * wk
    end
    Gdot = Kblk \ b
    res = norm(Kblk * Gdot - b) / norm(b)
    M = overlap_matrix(Xs, sigmas)
    @printf("  (2) block solve residual = %.3e, cond(K_block) = %.3e, cond(M_scalar) = %.3e\n",
            res, cond(Kblk), cond(M))
    ok &= res < 1e-10

    # Well-separated limit: off-diagonal blocks decay and Gdot_k -> J_k Gamma_k.
    # The decay is the ALGEBRAIC 1/r^3 Hessian tail (M4/M6 nonlocality), not the
    # Gaussian's exponential decay, so the residual shrinks like (sigma/d)^3 rather
    # than hitting machine zero. Check that scaling across two separations.
    function sep_residual(D)
        Xsf = [[0.0, 0.0, 0.0], [D, 0.0, 0.0]]
        Kf = assemble_block_K(Xsf, sigmas)
        bf = zeros(6)
        for k in 1:2
            wk = projected_gaussian_field(Xsf[k], Xsf, Gammas, sigmas)
            bf[3k-2:3k] = Js[k] * wk
        end
        Gdotf = Kf \ bf
        e = 0.0
        for k in 1:2
            e = max(e, norm(Gdotf[3k-2:3k] - Js[k] * Gammas[k]) / norm(Js[k] * Gammas[k]))
        end
        return e
    end
    e20, e40 = sep_residual(20.0), sep_residual(40.0)
    ratio = e20 / max(e40, eps())
    @printf("  (2) well-separated Gdot vs J Gamma relerr: d=20s %.3e, d=40s %.3e (ratio %.2f ~ 8 for 1/r^3 tail)\n",
            e20, e40, ratio)
    ok &= e40 < e20 && 5.0 < ratio < 11.0   # ~2^3 = 8 algebraic decay

    # (3) Core-spreading term d/dsigma(K_sigma Gamma). Verify it is finite and
    #     AD-consistent with a finite difference; it is the sigma-derivative basis
    #     drift that, with sigma^2 += 2 nu dt and Gamma fixed, realizes nu grad^2 omega
    #     (so it cancels rather than forcing Gdot -- inviscidly sigma_dot = 0).
    xprobe = [0.4, -0.2, 0.15]
    X0 = [[0.0, 0.0, 0.0]]
    dK_ad = ForwardDiff.derivative(
        s -> projected_gaussian_field(xprobe, X0, [Gamma], [s]), sigma)
    h = 1e-6
    fd = (projected_gaussian_field(xprobe, X0, [Gamma], [sigma + h]) .-
          projected_gaussian_field(xprobe, X0, [Gamma], [sigma - h])) ./ (2h)
    cs_err = norm(dK_ad - fd) / max(norm(fd), eps())
    @printf("  (3) d/dsigma(K Gamma) AD vs FD relerr = %.3e (term present; cancels viscously)\n", cs_err)
    ok &= cs_err < 1e-5

    @printf("M6 self-test -- %s\n", ok ? "PASS" : "FAIL")
    return ok
end

# Opt-in M12 theory self-check. This is a diagnostic companion to the theorem,
# not a production model: it illustrates that local/L1 exactly divergence-free
# kernels have zero volume moment, while physical circulation is a loop/surface
# Stokes quantity rather than the full-volume vector integral of omega.
function m12_selftest()
    ok = true
    sigma = 1.0
    X = [0.0, 0.0, 0.0]
    Gamma = [0.0, 0.0, 1.0]
    Xs = [X]
    Gammas = [Gamma]
    sigmas = [sigma]

    # (1) Volume moments: scalar zeta has moment Gamma; the local divergence-free
    #     matrix Gaussian has zero moment; M4's nonlocal projected blob has the
    #     symmetric-volume 2/3 Gamma convention caveat.
    n = 25
    L = 6.0 * sigma
    dx = 2L / n
    scalar_int = zeros(3)
    local_int = zeros(3)
    m4_int = zeros(3)
    for i in 1:n, j in 1:n, k in 1:n
        x = [-L + (i - 0.5) * dx, -L + (j - 0.5) * dx, -L + (k - 0.5) * dx]
        dV = dx^3
        scalar_int .+= omega_zeta(x, X, Gamma, sigma) .* dV
        local_int .+= local_divfree_gaussian_field(x, X, Gamma, sigma) .* dV
        m4_int .+= projected_gaussian_field(x, Xs, Gammas, sigmas) .* dV
    end
    scalar_err = norm(scalar_int - Gamma) / norm(Gamma)
    local_norm = norm(local_int) / norm(Gamma)
    m4_err = norm(m4_int - (2 / 3) .* Gamma) / norm(Gamma)
    @printf("  (1a) scalar zeta volume moment relerr to Gamma = %.3e\n", scalar_err)
    @printf("  (1b) local div-free Gaussian volume moment / ||Gamma|| = %.3e\n", local_norm)
    @printf("  (1c) M4 symmetric volume moment relerr to (2/3)Gamma = %.3e\n", m4_err)
    ok &= scalar_err < 1e-8 && local_norm < 3e-3 && m4_err < 3e-3

    # (2) Stokes check on a finite circular contour in the plane normal to Gamma:
    #     loop circulation of u equals vorticity flux of omega_J through the disk.
    #     This is the physical circulation object; it is not the full-volume
    #     vector moment tested above.
    R = 2.0 * sigma
    ntheta = 720
    nr = 160
    loop = 0.0
    for q in 1:ntheta
        theta = 2pi * (q - 0.5) / ntheta
        x = [R * cos(theta), R * sin(theta), 0.0]
        tangent = [-sin(theta), cos(theta), 0.0]
        loop += dot(induced_velocity(x, X, Gamma, sigma), tangent) * R * (2pi / ntheta)
    end
    flux = 0.0
    dr = R / nr
    dtheta = 2pi / ntheta
    for ir in 1:nr, it in 1:ntheta
        r = (ir - 0.5) * dr
        theta = (it - 0.5) * dtheta
        x = [r * cos(theta), r * sin(theta), 0.0]
        omega = projected_gaussian_field(x, Xs, Gammas, sigmas)
        flux += omega[3] * r * dr * dtheta
    end
    stokes_relerr = abs(loop - flux) / max(abs(loop), abs(flux), eps())
    @printf("  (2) finite-loop circulation vs Stokes flux relerr = %.3e (loop=%.6e, flux=%.6e)\n",
            stokes_relerr, loop, flux)
    ok &= stokes_relerr < 5e-4

    @printf("M12 self-test -- %s\n", ok ? "PASS" : "FAIL")
    return ok
end

main()

if get(ENV, "M1_SELFTEST", "0") != "0"
    println()
    println("M1 closed-form Hess(G_sigma) self-check:")
    m1_selftest()
end

if get(ENV, "M2_SELFTEST", "0") != "0"
    println()
    println("M2 effective-strength projection self-check:")
    m2_selftest()
end

if get(ENV, "M4_SELFTEST", "0") != "0"
    println()
    println("M4 projected-Gaussian representation self-check:")
    m4_selftest()
end

if get(ENV, "M6_SELFTEST", "0") != "0"
    println()
    println("M6 projected-basis strength-evolution self-check:")
    m6_selftest()
end

if get(ENV, "M12_SELFTEST", "0") != "0"
    println()
    println("M12 no-go theorem / circulation-moment self-check:")
    m12_selftest()
end
