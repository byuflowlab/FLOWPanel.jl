# =============================================================================
# Driver: M2 effective-strength diagnostic (BRAINSTORM item 010)
#
# Offline, read-only over SAVED particle VTP states. Quantifies whether
# Gamma_eff = Gamma + M^{-1} H Gamma is dominated by the universal 2/3 self
# correction or by neighbor overlap. This does not modify any particle update.
#
# Usage:
#   julia --project examples/particle_m2_effective_strength_diag.jl [RUN_NAME] [STEPS]
#   ENV overrides: RUN_NAME, STEPS (e.g. "340:359" or "359"), CUTOFF_FACTOR,
#                  NSTATES, TARGET_NSAMPLE, LOCAL_SOLVE_NSAMPLE,
#                  LOCAL_SOLVE_MAX_NEIGHBORS, RNG_SEED.
#
# Defaults to RUN_NAME=rotor_hover_pressure_comparison, last NSTATES=20 steps.
# =============================================================================

include(joinpath(@__DIR__, "particle_overlap_residual.jl"))
include(joinpath(@__DIR__, "particle_gauserf_hessian_helpers.jl"))
const POR = ParticleOverlapResidual

import LinearAlgebra: dot, norm
import Printf: @printf, @sprintf
import Random
import Statistics: mean
using .ParticleGauserfHessian: HAVE_SF, hess_green_times_gamma_components

relnorm(A, B; eps=1e-30) = norm(A) / (norm(B) + eps)
fieldnorm(A) = norm(A)

function cos_field(A, B)
    den = norm(A) * norm(B)
    return den > 0 ? dot(vec(A), vec(B)) / den : NaN
end

function _median(v)
    s = sort(collect(v)); n = length(s)
    return n == 0 ? NaN : (iseven(n) ? 0.5 * (s[n ÷ 2] + s[n ÷ 2 + 1]) : s[(n + 1) ÷ 2])
end

function selected_steps(part_dir, run_name, steps_arg, nstates)
    avail = Int[]
    for f in readdir(part_dir)
        m = match(Regex("$(run_name)_wake1_particles\\.(\\d+)\\.vtp\$"), f)
        m === nothing || push!(avail, parse(Int, m.captures[1]))
    end
    sort!(avail)
    isempty(avail) && error("No particle VTP files in $(part_dir)")

    steps =
        if !isempty(steps_arg)
            if occursin(':', steps_arg)
                lo, hi = parse.(Int, split(steps_arg, ':'))
                filter(s -> lo <= s <= hi, avail)
            else
                [parse(Int, steps_arg)]
            end
        else
            avail[max(1, end - nstates + 1):end]
        end
    isempty(steps) && error("No steps selected (requested '$(steps_arg)' from $(length(avail)) available).")
    return steps, avail
end

function select_indices(np::Int, nsample::Int, seed::Int)
    (nsample <= 0 || nsample >= np) && return collect(1:np)
    rng = Random.MersenneTwister(seed)
    idx = Random.randperm(rng, np)[1:nsample]
    sort!(idx)
    return idx
end

function collect_neighbors(g::POR.CSRGrid, X, k, cutoff2)
    kx, ky, kz = g.key[k]
    xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
    ids = Int[]
    d2s = Float64[]
    @inbounds for dz in -1:1, dy in -1:1, dx in -1:1
        c = get(g.cellid, (kx + dx, ky + dy, kz + dz), 0)
        c == 0 && continue
        for p in (g.start[c] + 1):g.start[c + 1]
            i = g.perm[p]
            rx = X[1,i] - xk1; ry = X[2,i] - xk2; rz = X[3,i] - xk3
            d2 = rx * rx + ry * ry + rz * rz
            d2 <= cutoff2 || continue
            push!(ids, i)
            push!(d2s, d2)
        end
    end
    return ids, d2s
end

function local_target_and_field(x, X, sigma, loc, coeffs, target_gamma)
    target = zeros(Float64, 3)
    field = zeros(Float64, 3)
    @inbounds for (j, i) in pairs(loc)
        rx = x[1] - X[1,i]; ry = x[2] - X[2,i]; rz = x[3] - X[3,i]
        d2 = rx * rx + ry * ry + rz * rz
        z = POR.zeta_sigma(d2, sigma[i])
        field[1] += z * coeffs[1,j]; field[2] += z * coeffs[2,j]; field[3] += z * coeffs[3,j]
        target[1] += z * target_gamma[1,j]
        target[2] += z * target_gamma[2,j]
        target[3] += z * target_gamma[3,j]
        q1, q2, q3 = hess_green_times_gamma_components(
            rx, ry, rz, sigma[i], target_gamma[1,j], target_gamma[2,j], target_gamma[3,j])
        target[1] += q1; target[2] += q2; target[3] += q3
    end
    return target, field
end

function local_m2_sample(state, g::POR.CSRGrid, cutoff2;
                         nsample::Int=0, max_neighbors::Int=160,
                         rng_seed::Int=1)
    nsample <= 0 && return (; local_n=0, local_rel_delta=NaN, local_res_raw=NaN,
                            local_res_self=NaN, local_res_jac=NaN,
                            local_solve_check=NaN, local_probe_n=0,
                            local_probe_res_raw=NaN, local_probe_res_self=NaN,
                            local_probe_res_jac=NaN, local_probe_res_solve=NaN)

    X, sigma, gamma, np = state.X, state.sigma, state.gamma, state.np
    centers = select_indices(np, nsample, rng_seed)
    raw2 = 0.0; self2 = 0.0; jac2 = 0.0; solve2 = 0.0; den2 = 0.0; delta2 = 0.0; base2 = 0.0
    praw2 = 0.0; pself2 = 0.0; pjac2 = 0.0; psolve2 = 0.0; pden2 = 0.0
    probe_used = 0
    used = 0

    for k in centers
        ids, d2s = collect_neighbors(g, X, k, cutoff2)
        isempty(ids) && continue
        order = sortperm(d2s)
        take = order[1:min(length(order), max_neighbors)]
        loc = ids[take]
        if !(k in loc)
            loc[end] = k
        end
        center = findfirst(==(k), loc)
        center === nothing && continue
        nloc = length(loc)

        M = zeros(Float64, nloc, nloc)
        W = zeros(Float64, 3, nloc)
        G = gamma[:, loc]
        @inbounds for a in 1:nloc, b in 1:nloc
            ib = loc[b]
            rx = X[1,loc[b]] - X[1,loc[a]]
            ry = X[2,loc[b]] - X[2,loc[a]]
            rz = X[3,loc[b]] - X[3,loc[a]]
            d2 = rx * rx + ry * ry + rz * rz
            M[a,b] = POR.zeta_sigma(d2, sigma[ib])
        end
        W .= G * transpose(M)
        @inbounds for a in 1:nloc
            ka = loc[a]
            h1 = 0.0; h2 = 0.0; h3 = 0.0
            for ib in loc
                rx = X[1,ka] - X[1,ib]
                ry = X[2,ka] - X[2,ib]
                rz = X[3,ka] - X[3,ib]
                q1, q2, q3 = hess_green_times_gamma_components(
                    rx, ry, rz, sigma[ib], gamma[1,ib], gamma[2,ib], gamma[3,ib])
                h1 += q1; h2 += q2; h3 += q3
            end
            W[1,a] += h1; W[2,a] += h2; W[3,a] += h3
        end

        Geff = permutedims(M \ transpose(W))
        wc = W[:, center]
        raw_field = vec(M[center, :]' * transpose(G))
        self_coeffs = (2 / 3) .* G
        jac_coeffs = copy(G)
        @inbounds for a in 1:nloc
            za = POR.zeta_sigma(0.0, sigma[loc[a]])
            wz_a = vec(M[a, :]' * transpose(G))
            hj_a = W[:, a] .- wz_a
            jac_coeffs[:, a] .= G[:, a] .+ hj_a ./ za
        end
        self_field = vec(M[center, :]' * transpose(self_coeffs))
        jac_field = vec(M[center, :]' * transpose(jac_coeffs))
        solve_field = vec(M[center, :]' * transpose(Geff))

        raw2 += norm(raw_field - wc)^2
        self2 += norm(self_field - wc)^2
        jac2 += norm(jac_field - wc)^2
        solve2 += norm(solve_field - wc)^2
        den2 += norm(wc)^2
        delta2 += norm(Geff[:, center] - G[:, center])^2
        base2 += norm(G[:, center])^2
        used += 1

        probe_order = [j for j in order if ids[j] != k]
        for j in probe_order[1:min(length(probe_order), 4)]
            ip = ids[j]
            probe = 0.5 .* (X[:, k] .+ X[:, ip])
            wt, wr = local_target_and_field(probe, X, sigma, loc, G, G)
            _, ws = local_target_and_field(probe, X, sigma, loc, self_coeffs, G)
            _, wj = local_target_and_field(probe, X, sigma, loc, jac_coeffs, G)
            _, wsolve = local_target_and_field(probe, X, sigma, loc, Geff, G)
            praw2 += norm(wr - wt)^2
            pself2 += norm(ws - wt)^2
            pjac2 += norm(wj - wt)^2
            psolve2 += norm(wsolve - wt)^2
            pden2 += norm(wt)^2
            probe_used += 1
        end
    end

    return (; local_n=used,
            local_rel_delta=sqrt(delta2 / max(base2, eps())),
            local_res_raw=sqrt(raw2 / max(den2, eps())),
            local_res_self=sqrt(self2 / max(den2, eps())),
            local_res_jac=sqrt(jac2 / max(den2, eps())),
            local_solve_check=sqrt(solve2 / max(den2, eps())),
            local_probe_n=probe_used,
            local_probe_res_raw=sqrt(praw2 / max(pden2, eps())),
            local_probe_res_self=sqrt(pself2 / max(pden2, eps())),
            local_probe_res_jac=sqrt(pjac2 / max(pden2, eps())),
            local_probe_res_solve=sqrt(psolve2 / max(pden2, eps())))
end

function analyze_m2(state; cutoff_factor::Float64=4.0,
                    target_nsample::Int=500,
                    local_solve_nsample::Int=0,
                    local_solve_max_neighbors::Int=160,
                    rng_seed::Int=1)
    np = state.np
    np == 0 && return nothing
    X, sigma, gamma = state.X, state.sigma, state.gamma
    sig_max = maximum(sigma); sig_min = minimum(sigma)
    cutoff = cutoff_factor * sig_max
    cutoff2 = cutoff * cutoff

    t0 = time()
    g = POR.build_csr(X, np, cutoff)
    centers = select_indices(np, target_nsample, rng_seed)
    nt = length(centers)
    omega = zeros(Float64, 3, nt)
    H = zeros(Float64, 3, nt)
    rowsum = zeros(Float64, nt)
    nbr = zeros(Int, nt)
    diag = Vector{Float64}(undef, nt)

    @inbounds for (icol, k) in pairs(centers)
        kx, ky, kz = g.key[k]
        xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
        o1 = 0.0; o2 = 0.0; o3 = 0.0
        h1 = 0.0; h2 = 0.0; h3 = 0.0
        rs = 0.0; cnt = 0
        for dz in -1:1, dy in -1:1, dx in -1:1
            c = get(g.cellid, (kx + dx, ky + dy, kz + dz), 0)
            c == 0 && continue
            for p in (g.start[c] + 1):g.start[c + 1]
                i = g.perm[p]
                rx_src = X[1,i] - xk1; ry_src = X[2,i] - xk2; rz_src = X[3,i] - xk3
                d2 = rx_src * rx_src + ry_src * ry_src + rz_src * rz_src
                d2 <= cutoff2 || continue
                si = sigma[i]
                z = POR.zeta_sigma(d2, si)
                o1 += z * gamma[1,i]; o2 += z * gamma[2,i]; o3 += z * gamma[3,i]
                q1, q2, q3 = hess_green_times_gamma_components(
                    -rx_src, -ry_src, -rz_src, si, gamma[1,i], gamma[2,i], gamma[3,i])
                h1 += q1; h2 += q2; h3 += q3
                rs += z; cnt += 1
            end
        end
        omega[1,icol] = o1; omega[2,icol] = o2; omega[3,icol] = o3
        H[1,icol] = h1; H[2,icol] = h2; H[3,icol] = h3
        rowsum[icol] = rs; nbr[icol] = cnt
        diag[icol] = POR.zeta_sigma(0.0, sigma[k])
    end

    Hself = similar(H)
    Gamma_centers = zeros(Float64, 3, nt)
    @inbounds for (icol, k) in pairs(centers)
        Gamma_centers[1,icol] = gamma[1,k]
        Gamma_centers[2,icol] = gamma[2,k]
        Gamma_centers[3,icol] = gamma[3,k]
        Hself[1,icol] = -(1 / 3) * diag[icol] * gamma[1,k]
        Hself[2,icol] = -(1 / 3) * diag[icol] * gamma[2,k]
        Hself[3,icol] = -(1 / 3) * diag[icol] * gamma[3,k]
    end
    Hnbr = H .- Hself
    target = omega .+ H

    omega_self = (2 / 3) .* omega
    residual_jac = NaN
    gamma_eff_jac_delta = NaN

    omega_curl_all = POR.curl_from_J(state)
    omega_curl = omega_curl_all[:, centers]
    domin = diag ./ max.(rowsum .- diag, eps())
    localres = local_m2_sample(state, g, cutoff2;
        nsample=local_solve_nsample,
        max_neighbors=local_solve_max_neighbors,
        rng_seed=rng_seed + 10_000)

    elapsed = time() - t0
    return (; np, ntarget=nt, cutoff, sig_min, sig_max,
            basis_curl_relerr=POR.relnorm(omega .- omega_curl, omega_curl),
            h_total_over_omega=fieldnorm(H) / max(fieldnorm(omega), eps()),
            h_self_over_omega=fieldnorm(Hself) / max(fieldnorm(omega), eps()),
            h_neighbor_over_self=fieldnorm(Hnbr) / max(fieldnorm(Hself), eps()),
            h_neighbor_over_total=fieldnorm(Hnbr) / max(fieldnorm(H), eps()),
            cos_h_total_self=cos_field(H, Hself),
            cos_h_total_gamma=cos_field(H, Gamma_centers),
            target_vs_curl_relerr=POR.relnorm(target .- omega_curl, omega_curl),
            gamma_eff_jac_delta,
            residual_raw=POR.relnorm(omega .- target, target),
            residual_self=POR.relnorm(omega_self .- target, target),
            residual_jac,
            domin_min=minimum(domin), domin_med=_median(domin),
            nbr_mean=sum(nbr) / nt, nbr_max=maximum(nbr),
            localres..., elapsed)
end

function synthetic_gates()
    println("\n[gate 1] isolated particle gives H_neighbor=0 and Jacobi/self=2/3")
    X = reshape([0.0, 0.0, 0.0], 3, 1)
    gamma = reshape([1.0, 0.2, -0.1], 3, 1)
    sigma = [1.0]
    z0 = POR.zeta_sigma(0.0, 1.0)
    vel = zeros(3, 1); vort = zeros(3, 1); J = zeros(9, 1)
    st = (; X, gamma, sigma, velocity=vel, vorticity=vort, J, np=1)
    r = analyze_m2(st; cutoff_factor=4.0, target_nsample=0)
    @printf("    h_neighbor/self=%.3e residual_self=%.3e residual_jac=%.3e\n",
            r.h_neighbor_over_self, r.residual_self, r.residual_jac)
    @assert r.h_neighbor_over_self < 1e-12
    @assert r.residual_self < 1e-12
    @assert abs(r.h_self_over_omega - 1/3) < 1e-12
    @assert abs(z0 - POR.ZETA0) < 1e-14
    println("    PASS")

    println("\n[gate 2] two-particle local solve improves over raw/self/Jacobi")
    X2 = [0.0 0.8; 0.0 0.3; 0.0 -0.2]
    G2 = [1.0 -0.4; 0.2 0.7; 0.0 0.5]
    s2 = [1.0, 1.0]
    st2 = (; X=X2, gamma=G2, sigma=s2, velocity=zeros(3, 2), vorticity=zeros(3, 2),
             J=zeros(9, 2), np=2)
    r2 = analyze_m2(st2; cutoff_factor=4.0, target_nsample=0,
                    local_solve_nsample=2, local_solve_max_neighbors=2)
    @printf("    local residuals raw=%.3e self=%.3e jac=%.3e solve=%.3e\n",
            r2.local_res_raw, r2.local_res_self, r2.local_res_jac, r2.local_solve_check)
    @assert r2.local_n == 2
    @assert r2.local_solve_check < 1e-12
    @assert r2.local_solve_check < r2.local_res_jac
    println("    PASS")
end

run_name = length(ARGS) >= 1 ? ARGS[1] : get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
steps_arg = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "STEPS", "")
cutoff_factor = parse(Float64, get(ENV, "CUTOFF_FACTOR", "4.0"))
nstates = parse(Int, get(ENV, "NSTATES", "20"))
local_solve_nsample = parse(Int, get(ENV, "LOCAL_SOLVE_NSAMPLE", "0"))
local_solve_max_neighbors = parse(Int, get(ENV, "LOCAL_SOLVE_MAX_NEIGHBORS", "160"))
target_nsample = parse(Int, get(ENV, "TARGET_NSAMPLE", "500"))
rng_seed = parse(Int, get(ENV, "RNG_SEED", "1"))

data_root = joinpath(@__DIR__, "..", "data")
part_dir = joinpath(data_root, run_name, "$(run_name)_wake1_particles")
isdir(part_dir) || error("Particle dir not found: $(part_dir)")
vtp_path(idx) = joinpath(part_dir, "$(run_name)_wake1_particles.$(idx).vtp")
steps, avail = selected_steps(part_dir, run_name, steps_arg, nstates)

println("="^78)
println("M2 effective-strength diagnostic (item 010)")
println("  run                 = $(run_name)")
println("  steps               = $(first(steps))...$(last(steps)) ($(length(steps)) states; $(length(avail)) available)")
println("  cutoff_fac          = $(cutoff_factor)")
println("  target_nsample      = $(target_nsample <= 0 ? "all" : string(target_nsample))")
println("  local_solve_nsample = $(local_solve_nsample)")
println("  local_max_neighbors = $(local_solve_max_neighbors)")
println("  rng_seed            = $(rng_seed)")
println("  erf source           = ", HAVE_SF ? "SpecialFunctions.erf" : "self-contained series")
println("="^78)

synthetic_gates()

csv_path = joinpath(data_root, run_name, "particle_m2_effective_strength_diag.csv")
header = ["step","np","ntarget","cutoff","sig_min","sig_max","basis_curl_relerr",
          "target_vs_curl_relerr",
          "h_total_over_omega","h_self_over_omega","h_neighbor_over_self",
          "h_neighbor_over_total","cos_h_total_self","cos_h_total_gamma",
          "gamma_eff_jac_delta","residual_raw","residual_self","residual_jac",
          "domin_min","domin_med","nbr_mean","nbr_max",
          "local_n","local_rel_delta","local_res_raw","local_res_self","local_res_jac",
          "local_solve_check","local_probe_n","local_probe_res_raw",
          "local_probe_res_self","local_probe_res_jac","local_probe_res_solve",
          "elapsed_s"]

rows = NamedTuple[]
println("\n[diagnostic] per-state M2 metrics")
@printf("  %5s %7s %9s %9s %9s %9s %9s %9s %9s %9s\n",
        "step", "ntarg", "basisCJ", "targetCJ", "H/omega", "Hnbr/Hs",
        "raw", "self", "jac", "nbr")
for s in steps
    st = POR.load_state(vtp_path(s))
    res = analyze_m2(st; cutoff_factor, target_nsample,
                     local_solve_nsample, local_solve_max_neighbors, rng_seed)
    res === nothing && (println("  step $s: 0 particles, skipped"); continue)
    push!(rows, (; step=s, res...))
    @printf("  %5d %7d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.0f\n",
            s, res.ntarget, res.basis_curl_relerr, res.target_vs_curl_relerr,
            res.h_total_over_omega,
            res.h_neighbor_over_self, res.residual_raw, res.residual_self,
            res.residual_jac, res.nbr_mean)
end

open(csv_path, "w") do io
    println(io, join(header, ","))
    for r in rows
        println(io, join((r.step, r.np, r.ntarget, r.cutoff, r.sig_min, r.sig_max,
                          r.basis_curl_relerr, r.target_vs_curl_relerr,
                          r.h_total_over_omega, r.h_self_over_omega,
                          r.h_neighbor_over_self, r.h_neighbor_over_total,
                          r.cos_h_total_self, r.cos_h_total_gamma,
                          r.gamma_eff_jac_delta, r.residual_raw,
                          r.residual_self, r.residual_jac,
                          r.domin_min, r.domin_med, r.nbr_mean, r.nbr_max,
                          r.local_n, r.local_rel_delta, r.local_res_raw,
                          r.local_res_self, r.local_res_jac, r.local_solve_check,
                          r.local_probe_n, r.local_probe_res_raw,
                          r.local_probe_res_self, r.local_probe_res_jac,
                          r.local_probe_res_solve,
                          r.elapsed), ","))
    end
end

if !isempty(rows)
    mean_field(f) = mean(getfield.(rows, f))
    println("\n" * "="^78)
    println("SUMMARY over $(length(rows)) states")
    @printf("  basis omega vs curl-of-J rel-err : mean %.3e\n", mean_field(:basis_curl_relerr))
    @printf("  target (omega+H) vs curl-of-J    : mean %.3e\n", mean_field(:target_vs_curl_relerr))
    @printf("  ||H Gamma|| / ||M Gamma||        : mean %.3e\n", mean_field(:h_total_over_omega))
    @printf("  ||H_neighbor|| / ||H_self||      : mean %.3e\n", mean_field(:h_neighbor_over_self))
    @printf("  cos(H_total, H_self)             : mean %.3f\n", mean_field(:cos_h_total_self))
    if all(r -> isnan(r.residual_jac), rows)
        println("  global Jacobi residual           : not evaluated in sampled-target mode")
        @printf("  residual raw / self              : %.3e / %.3e\n",
                mean_field(:residual_raw), mean_field(:residual_self))
    else
        @printf("  Jacobi ||Gamma_eff-Gamma||/||G|| : mean %.3e\n", mean_field(:gamma_eff_jac_delta))
        @printf("  residual raw / self / Jacobi     : %.3e / %.3e / %.3e\n",
                mean_field(:residual_raw), mean_field(:residual_self), mean_field(:residual_jac))
    end
    if any(r -> r.local_n > 0, rows)
        @printf("  local row residual raw/self/Jac/check: %.3e / %.3e / %.3e / %.3e\n",
                mean_field(:local_res_raw), mean_field(:local_res_self),
                mean_field(:local_res_jac), mean_field(:local_solve_check))
        @printf("  local probe residual raw/self/Jac/solve: %.3e / %.3e / %.3e / %.3e\n",
                mean_field(:local_probe_res_raw), mean_field(:local_probe_res_self),
                mean_field(:local_probe_res_jac), mean_field(:local_probe_res_solve))
    end
    @printf("  diag-dominance median            : mean %.3e\n", mean_field(:domin_med))
    @printf("  neighbors mean / max             : %.1f / %.0f\n", mean_field(:nbr_mean), mean_field(:nbr_max))
    println("  CSV: $(csv_path)")
    println("="^78)
end
