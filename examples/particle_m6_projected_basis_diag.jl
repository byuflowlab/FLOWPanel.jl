# =============================================================================
# Driver: M6 projected-basis strength-evolution diagnostic (BRAINSTORM item 010)
#
# Offline, read-only over SAVED particle VTP states. M6 promotes the M4 Leray-
# projected blob to the carried model and poses the strength-RATE equation in the
# matrix kernel K_sigma = zeta_sigma I + Hess(G_sigma):
#
#   sum_i K_{sigma_i}(X_k - X_i) Gdot_i = J(X_k) omega(X_k) - convection - viscous,
#   omega(x) = sum_i K_{sigma_i}(x - X_i) Gamma_i   (= curl(u), the M1 field).
#
# This is the block-kernel (3x3 blocks) analogue of item 008's scalar M Gdot = b.
# The diagnostic characterizes, on a settled rotor-hover state, three rate
# estimates and the block-operator conditioning:
#   * self   : Gdot_k = J_k Gamma_k                 (isolated VPM stretch)
#   * jac    : Gdot_k = K(0)^{-1} b_k = (3/(2 zeta(0))) J_k omega(X_k)
#                                                    (block-diagonal / block-Jacobi)
#   * solve  : center row of the local block solve  (overlap-aware)
#
# The convective basis-motion term -sum_i [(u_k-u_i).grad K_i] Gamma_i is the
# matrix analogue of POR's scalar conv term; it needs the 3rd-derivative tensor
# grad Hess(G) and is DEFERRED here (this is the "reduced" RHS, exactly as item
# 008 staged reduced-vs-sampled). It does NOT change any particle update.
#
# Usage:
#   julia --project examples/particle_m6_projected_basis_diag.jl [RUN_NAME] [STEPS]
#   ENV overrides: RUN_NAME, STEPS ("340:359" or "359"), CUTOFF_FACTOR, NSTATES,
#                  TARGET_NSAMPLE, BLOCK_SOLVE_NSAMPLE, BLOCK_MAX_NEIGHBORS, RNG_SEED.
# Defaults to RUN_NAME=rotor_hover_pressure_comparison, last NSTATES=20 steps.
# =============================================================================

include(joinpath(@__DIR__, "particle_overlap_residual.jl"))
include(joinpath(@__DIR__, "particle_gauserf_hessian_helpers.jl"))
const POR = ParticleOverlapResidual

import LinearAlgebra: dot, norm, cond, I
import Printf: @printf
import Random
import Statistics: mean
using .ParticleGauserfHessian: HAVE_SF, hess_green_times_gamma_components

# --- shared helpers ----------------------------------------------------------
function _median(v)
    s = sort(collect(v)); n = length(s)
    return n == 0 ? NaN : (iseven(n) ? 0.5 * (s[n ÷ 2] + s[n ÷ 2 + 1]) : s[(n + 1) ÷ 2])
end

# physical stretching (w.grad)u = M w with M[a,b]=du_a/dx_b (the _euler classic branch)
@inline s_phys(J, w) = (J[1]*w[1] + J[4]*w[2] + J[7]*w[3],
                        J[2]*w[1] + J[5]*w[2] + J[8]*w[3],
                        J[3]*w[1] + J[6]*w[2] + J[9]*w[3])

# 3x3 kernel block K_sigma(r) = zeta_sigma(|r|^2/sigma^2) I + Hess(G_sigma)(r).
function kblock33(rx, ry, rz, sigma)
    d2 = rx*rx + ry*ry + rz*rz
    z = POR.zeta_sigma(d2, sigma)
    K = zeros(3, 3)
    c1 = hess_green_times_gamma_components(rx, ry, rz, sigma, 1.0, 0.0, 0.0)
    c2 = hess_green_times_gamma_components(rx, ry, rz, sigma, 0.0, 1.0, 0.0)
    c3 = hess_green_times_gamma_components(rx, ry, rz, sigma, 0.0, 0.0, 1.0)
    @inbounds for a in 1:3
        K[a, 1] = c1[a]; K[a, 2] = c2[a]; K[a, 3] = c3[a]
        K[a, a] += z
    end
    return K
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
    ids = Int[]; d2s = Float64[]
    @inbounds for dz in -1:1, dy in -1:1, dx in -1:1
        c = get(g.cellid, (kx + dx, ky + dy, kz + dz), 0)
        c == 0 && continue
        for p in (g.start[c] + 1):g.start[c + 1]
            i = g.perm[p]
            rx = X[1,i] - xk1; ry = X[2,i] - xk2; rz = X[3,i] - xk3
            d2 = rx*rx + ry*ry + rz*rz
            d2 <= cutoff2 || continue
            push!(ids, i); push!(d2s, d2)
        end
    end
    return ids, d2s
end

# --- local block solve over a truncated neighbor cloud -----------------------
# Returns rate-estimate residuals (block solve vs self vs block-Jacobi), the
# block-solve linear residual (solve check), block & scalar conditioning, and the
# block diagonal-dominance proxy.
function local_m6_sample(state, g::POR.CSRGrid, cutoff2;
                         nsample::Int=0, max_neighbors::Int=160, rng_seed::Int=1)
    nan = NaN
    nsample <= 0 && return (; m6_local_n=0, gdot_solve_vs_self=nan,
                            gdot_jac_vs_self=nan, gdot_solve_vs_jac=nan,
                            block_solve_res=nan, block_cond=nan, scalar_cond=nan,
                            block_domin_med=nan, omega_vs_isolated=nan)

    X, sigma, gamma, J, np = state.X, state.sigma, state.gamma, state.J, state.np
    centers = select_indices(np, nsample, rng_seed)
    n_solve2 = 0.0; n_jac2 = 0.0; n_sj2 = 0.0; den_self2 = 0.0; den_solve2 = 0.0
    res_acc = 0.0; omega_iso2 = 0.0; omega_den2 = 0.0
    bconds = Float64[]; sconds = Float64[]; bdom = Float64[]
    used = 0

    for k in centers
        ids, d2s = collect_neighbors(g, X, k, cutoff2)
        isempty(ids) && continue
        order = sortperm(d2s)
        loc = ids[order[1:min(length(order), max_neighbors)]]
        (k in loc) || (loc[end] = k)
        center = findfirst(==(k), loc)
        center === nothing && continue
        nloc = length(loc)

        # block operator Kb (3nloc x 3nloc), block[a,b] = K_{sigma_b}(X_a - X_b),
        # and scalar overlap Ms (nloc x nloc) for a conditioning comparison.
        Kb = zeros(3nloc, 3nloc)
        Ms = zeros(nloc, nloc)
        @inbounds for a in 1:nloc, b in 1:nloc
            ia = loc[a]; ib = loc[b]
            rx = X[1,ia] - X[1,ib]; ry = X[2,ia] - X[2,ib]; rz = X[3,ia] - X[3,ib]
            blk = kblock33(rx, ry, rz, sigma[ib])
            Kb[3a-2:3a, 3b-2:3b] = blk
            Ms[a, b] = POR.zeta_sigma(rx*rx + ry*ry + rz*rz, sigma[ib])
        end

        # stacked Gamma, represented omega = Kb * Gamma, reduced RHS b_a = J_a . omega_a.
        Gst = zeros(3nloc)
        @inbounds for a in 1:nloc
            Gst[3a-2] = gamma[1,loc[a]]; Gst[3a-1] = gamma[2,loc[a]]; Gst[3a] = gamma[3,loc[a]]
        end
        omega = Kb * Gst
        bst = zeros(3nloc)
        @inbounds for a in 1:nloc
            Ja = view(J, :, loc[a])
            wa = (omega[3a-2], omega[3a-1], omega[3a])
            s = s_phys(Ja, wa)
            bst[3a-2] = s[1]; bst[3a-1] = s[2]; bst[3a] = s[3]
        end

        Gdot = Kb \ bst
        res_acc += norm(Kb * Gdot - bst) / max(norm(bst), eps())

        # center-row rate estimates
        Jc = view(J, :, k); Gc = view(gamma, :, k)
        self = collect(s_phys(Jc, Gc))
        z0 = POR.zeta_sigma(0.0, sigma[k])
        bc = bst[3center-2:3center]
        jac = (3 / (2 * z0)) .* bc                  # K(0)^{-1} b = block-Jacobi
        solve = Gdot[3center-2:3center]

        n_solve2 += norm(solve - self)^2
        n_jac2 += norm(jac - self)^2
        n_sj2 += norm(solve - jac)^2
        den_self2 += norm(self)^2
        den_solve2 += norm(solve)^2

        # how non-isolated is the represented vorticity here: omega(X_k) vs K(0)Gamma_k
        wc = omega[3center-2:3center]
        omega_iso2 += norm(wc - (2/3) * z0 .* collect(Gc))^2
        omega_den2 += norm(wc)^2

        push!(bconds, cond(Kb)); push!(sconds, cond(Ms))
        # block diagonal dominance: ||K(0)|| / sum_{b!=center} ||block(center,b)||
        diagn = norm(Kb[3center-2:3center, 3center-2:3center])
        offsum = 0.0
        for b in 1:nloc
            b == center && continue
            offsum += norm(Kb[3center-2:3center, 3b-2:3b])
        end
        push!(bdom, diagn / max(offsum, eps()))
        used += 1
    end

    used == 0 && return (; m6_local_n=0, gdot_solve_vs_self=nan,
                         gdot_jac_vs_self=nan, gdot_solve_vs_jac=nan,
                         block_solve_res=nan, block_cond=nan, scalar_cond=nan,
                         block_domin_med=nan, omega_vs_isolated=nan)
    return (; m6_local_n=used,
            gdot_solve_vs_self=sqrt(n_solve2 / max(den_self2, eps())),
            gdot_jac_vs_self=sqrt(n_jac2 / max(den_self2, eps())),
            gdot_solve_vs_jac=sqrt(n_sj2 / max(den_solve2, eps())),
            block_solve_res=res_acc / used,
            block_cond=_median(bconds), scalar_cond=_median(sconds),
            block_domin_med=_median(bdom),
            omega_vs_isolated=sqrt(omega_iso2 / max(omega_den2, eps())))
end

# --- per-state full-cloud (sampled-target) metrics ---------------------------
function analyze_m6(state; cutoff_factor::Float64=4.0, target_nsample::Int=500,
                    block_solve_nsample::Int=0, block_max_neighbors::Int=160,
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
    omega = zeros(3, nt); H = zeros(3, nt)
    nbr = zeros(Int, nt); rowsum = zeros(nt)

    @inbounds for (icol, k) in pairs(centers)
        kx, ky, kz = g.key[k]
        xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
        o1=0.0; o2=0.0; o3=0.0; h1=0.0; h2=0.0; h3=0.0; rs=0.0; cnt=0
        for dz in -1:1, dy in -1:1, dx in -1:1
            c = get(g.cellid, (kx+dx, ky+dy, kz+dz), 0)
            c == 0 && continue
            for p in (g.start[c]+1):g.start[c+1]
                i = g.perm[p]
                rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
                d2 = rx*rx + ry*ry + rz*rz
                d2 <= cutoff2 || continue
                si = sigma[i]
                z = POR.zeta_sigma(d2, si)
                o1 += z*gamma[1,i]; o2 += z*gamma[2,i]; o3 += z*gamma[3,i]
                q1, q2, q3 = hess_green_times_gamma_components(
                    -rx, -ry, -rz, si, gamma[1,i], gamma[2,i], gamma[3,i])
                h1 += q1; h2 += q2; h3 += q3; rs += z; cnt += 1
            end
        end
        omega[1,icol]=o1; omega[2,icol]=o2; omega[3,icol]=o3
        H[1,icol]=h1; H[2,icol]=h2; H[3,icol]=h3
        nbr[icol]=cnt; rowsum[icol]=rs
    end

    omega_curl = POR.curl_from_J(state)[:, centers]
    target = omega .+ H                                  # represented vorticity = M1 field
    diag = [POR.zeta_sigma(0.0, sigma[k]) for k in centers]
    domin = diag ./ max.(rowsum .- diag, eps())

    loc = local_m6_sample(state, g, cutoff2; nsample=block_solve_nsample,
                          max_neighbors=block_max_neighbors, rng_seed=rng_seed+10_000)
    elapsed = time() - t0
    return (; np, ntarget=nt, cutoff, sig_min, sig_max,
            basis_curl_relerr=POR.relnorm(omega .- omega_curl, omega_curl),
            target_vs_curl_relerr=POR.relnorm(target .- omega_curl, omega_curl),
            h_over_omega=norm(H) / max(norm(omega), eps()),
            domin_med=_median(domin), nbr_mean=sum(nbr)/nt, nbr_max=maximum(nbr),
            loc..., elapsed)
end

# --- synthetic gates ---------------------------------------------------------
function synthetic_gates()
    println("\n[gate 1] isolated particle: block solve recovers Gdot = J Gamma")
    X = reshape([0.0, 0.0, 0.0], 3, 1)
    gamma = reshape([1.0, 0.2, -0.1], 3, 1)
    sigma = [1.0]
    Jvec = [0.3, 0.05, -0.15, -0.2, -0.4, 0.2, 0.1, 0.25, 0.35]   # M[a,b]=du_a/dx_b
    J = reshape(Jvec, 9, 1)
    st = (; X, gamma, sigma, velocity=zeros(3,1), vorticity=zeros(3,1), J, np=1)
    r = analyze_m6(st; cutoff_factor=4.0, target_nsample=0,
                   block_solve_nsample=1, block_max_neighbors=1)
    @printf("    solve_vs_self=%.3e block_solve_res=%.3e block_cond=%.3e\n",
            r.gdot_solve_vs_self, r.block_solve_res, r.block_cond)
    @assert r.m6_local_n == 1
    @assert r.gdot_solve_vs_self < 1e-12      # no overlap => solve == self
    @assert r.block_solve_res < 1e-10
    @assert abs(r.block_cond - 1.0) < 1e-9
    println("    PASS")

    println("\n[gate 2] two particles: local block solve fits its rows (roundoff)")
    X2 = [0.0 0.8; 0.0 0.3; 0.0 -0.2]
    G2 = [1.0 -0.4; 0.2 0.7; 0.0 0.5]
    s2 = [1.0, 1.0]
    J2 = hcat(Jvec, [-0.1, 0.3, 0.1, 0.2, 0.15, -0.05, 0.05, -0.2, 0.25])
    st2 = (; X=X2, gamma=G2, sigma=s2, velocity=zeros(3,2), vorticity=zeros(3,2),
             J=J2, np=2)
    r2 = analyze_m6(st2; cutoff_factor=4.0, target_nsample=0,
                    block_solve_nsample=2, block_max_neighbors=2)
    @printf("    block_solve_res=%.3e block_cond=%.3e scalar_cond=%.3e solve_vs_self=%.3e\n",
            r2.block_solve_res, r2.block_cond, r2.scalar_cond, r2.gdot_solve_vs_self)
    @assert r2.m6_local_n == 2
    @assert r2.block_solve_res < 1e-10
    @assert r2.gdot_solve_vs_self > 0      # overlap actually changes the rate
    println("    PASS")
end

# --- driver ------------------------------------------------------------------
run_name = length(ARGS) >= 1 ? ARGS[1] : get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
steps_arg = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "STEPS", "")
cutoff_factor = parse(Float64, get(ENV, "CUTOFF_FACTOR", "4.0"))
nstates = parse(Int, get(ENV, "NSTATES", "20"))
target_nsample = parse(Int, get(ENV, "TARGET_NSAMPLE", "500"))
block_solve_nsample = parse(Int, get(ENV, "BLOCK_SOLVE_NSAMPLE", "20"))
block_max_neighbors = parse(Int, get(ENV, "BLOCK_MAX_NEIGHBORS", "160"))
rng_seed = parse(Int, get(ENV, "RNG_SEED", "1"))

data_root = joinpath(@__DIR__, "..", "data")
part_dir = joinpath(data_root, run_name, "$(run_name)_wake1_particles")
isdir(part_dir) || error("Particle dir not found: $(part_dir)")
vtp_path(idx) = joinpath(part_dir, "$(run_name)_wake1_particles.$(idx).vtp")

avail = Int[]
for f in readdir(part_dir)
    m = match(Regex("$(run_name)_wake1_particles\\.(\\d+)\\.vtp\$"), f)
    m === nothing || push!(avail, parse(Int, m.captures[1]))
end
sort!(avail)
isempty(avail) && error("No particle VTP files in $(part_dir)")
steps = isempty(steps_arg) ? avail[max(1, end-nstates+1):end] :
        (occursin(':', steps_arg) ?
            (lh = parse.(Int, split(steps_arg, ':')); filter(s -> lh[1] <= s <= lh[2], avail)) :
            [parse(Int, steps_arg)])
isempty(steps) && error("No steps selected from $(length(avail)) available.")

println("="^78)
println("M6 projected-basis strength-evolution diagnostic (item 010)")
println("  run                 = $(run_name)")
println("  steps               = $(first(steps))...$(last(steps)) ($(length(steps)) states; $(length(avail)) available)")
println("  cutoff_fac          = $(cutoff_factor)")
println("  target_nsample      = $(target_nsample <= 0 ? "all" : string(target_nsample))")
println("  block_solve_nsample = $(block_solve_nsample)")
println("  block_max_neighbors = $(block_max_neighbors)")
println("  rng_seed            = $(rng_seed)")
println("  erf source          = ", HAVE_SF ? "SpecialFunctions.erf" : "self-contained series")
println("="^78)

synthetic_gates()

csv_path = joinpath(data_root, run_name, "particle_m6_projected_basis_diag.csv")
header = ["step","np","ntarget","cutoff","sig_min","sig_max","basis_curl_relerr",
          "target_vs_curl_relerr","h_over_omega","domin_med","nbr_mean","nbr_max",
          "m6_local_n","gdot_solve_vs_self","gdot_jac_vs_self","gdot_solve_vs_jac",
          "block_solve_res","block_cond","scalar_cond","block_domin_med",
          "omega_vs_isolated","elapsed_s"]

rows = NamedTuple[]
println("\n[diagnostic] per-state M6 metrics")
@printf("  %5s %7s %9s %9s %9s %9s %9s %9s %9s\n",
        "step", "ntarg", "basisCJ", "tgtCJ", "solve/self", "jac/self",
        "blkcond", "blkdom", "nbr")
for s in steps
    st = POR.load_state(vtp_path(s))
    res = analyze_m6(st; cutoff_factor, target_nsample,
                     block_solve_nsample, block_max_neighbors, rng_seed)
    res === nothing && (println("  step $s: 0 particles, skipped"); continue)
    push!(rows, (; step=s, res...))
    @printf("  %5d %7d %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.0f\n",
            s, res.ntarget, res.basis_curl_relerr, res.target_vs_curl_relerr,
            res.gdot_solve_vs_self, res.gdot_jac_vs_self, res.block_cond,
            res.block_domin_med, res.nbr_mean)
end

open(csv_path, "w") do io
    println(io, join(header, ","))
    for r in rows
        println(io, join((r.step, r.np, r.ntarget, r.cutoff, r.sig_min, r.sig_max,
                          r.basis_curl_relerr, r.target_vs_curl_relerr, r.h_over_omega,
                          r.domin_med, r.nbr_mean, r.nbr_max, r.m6_local_n,
                          r.gdot_solve_vs_self, r.gdot_jac_vs_self, r.gdot_solve_vs_jac,
                          r.block_solve_res, r.block_cond, r.scalar_cond,
                          r.block_domin_med, r.omega_vs_isolated, r.elapsed), ","))
    end
end

if !isempty(rows)
    mf(f) = mean(getfield.(rows, f))
    println("\n" * "="^78)
    println("SUMMARY over $(length(rows)) states")
    @printf("  basis omega vs curl-of-J         : mean %.3e\n", mf(:basis_curl_relerr))
    @printf("  represented (omega+H) vs curl-J  : mean %.3e\n", mf(:target_vs_curl_relerr))
    @printf("  ||H||/||omega||                  : mean %.3e\n", mf(:h_over_omega))
    @printf("  Gdot solve-vs-self / jac-vs-self : %.3e / %.3e\n",
            mf(:gdot_solve_vs_self), mf(:gdot_jac_vs_self))
    @printf("  Gdot solve-vs-jac                : %.3e\n", mf(:gdot_solve_vs_jac))
    @printf("  block solve residual (solve check): %.3e\n", mf(:block_solve_res))
    @printf("  block cond / scalar cond (median): %.3e / %.3e\n",
            mf(:block_cond), mf(:scalar_cond))
    @printf("  block diag-dominance (median)    : %.3e\n", mf(:block_domin_med))
    @printf("  represented omega vs isolated    : %.3e\n", mf(:omega_vs_isolated))
    @printf("  scalar diag-dominance median     : %.3e\n", mf(:domin_med))
    @printf("  neighbors mean / max             : %.1f / %.0f\n", mf(:nbr_mean), mf(:nbr_max))
    println("  CSV: $(csv_path)")
    println("="^78)
end
