# =============================================================================
# Driver: Pedrizzetti vs overlap-residual efficacy (BRAINSTORM item 011)
#
# Offline, read-only over SAVED particle VTP states. Quantifies how much of the
# 008/M2/M6 overlap-aware strength correction the already-shipped Pedrizzetti
# relaxation captures on its own, gating whether an explicit overlap solve is
# worth its cost. Does NOT modify any particle update.
#
# Verified prereq (FLOWVPM.jl dev v4.1.0, src/FLOWVPM_relaxation.jl:41):
# relax_pedrizzetti rotates Gamma toward the curl-of-J  omega_J = (J6-J8, J7-J3,
# J2-J4) = curl(u), magnitude-preserving:
#   Gamma_i <- (1-f) Gamma_i + f ||Gamma_i|| omega_hat_J(X_i).
# corrected-Pedrizzetti adds a /sqrt(b2) renorm (restores some magnitude).
#
# Three framings (see BRAINSTORM/011):
#   A. Rate-residual: does adding the Pedrizzetti effective rate reduce the 008
#      residual ||M Gdot - b|| relative to the self-only stretch J Gamma?
#   B. Alignment decomposition (parameter-free): angle ∠(Gamma, omega_J) and the
#      ∥(magnitude, NOT addressable by a norm-preserving rotation) / ⊥(direction,
#      addressable) split of the overlap correction. Primary deliverable.
#   C. Local-solve cross-check: does the Pedrizzetti rotation point toward the M2
#      local overlap target Gamma_eff?  (local-solve reference, per M2 caveat.)
#
# Usage:
#   julia --project examples/particle_pedrizzetti_overlap_diag.jl [RUN_NAME] [STEPS]
#   ENV: RUN_NAME, STEPS ("340:359"/"359"), CUTOFF_FACTOR, NSTATES, TARGET_NSAMPLE,
#        RLXF, NSTEPS_RELAX, DT, RPM, NT, LOCAL_SOLVE_NSAMPLE, LOCAL_MAX_NEIGHBORS,
#        RNG_SEED.
# =============================================================================

include(joinpath(@__DIR__, "particle_overlap_residual.jl"))
include(joinpath(@__DIR__, "particle_gauserf_hessian_helpers.jl"))
const POR = ParticleOverlapResidual

import LinearAlgebra: dot, norm
import Printf: @printf
import Random
import Statistics: mean, quantile
using .ParticleGauserfHessian: HAVE_SF, hess_green_times_gamma_components

# --- helpers -----------------------------------------------------------------
function _median(v)
    s = sort(collect(v)); n = length(s)
    return n == 0 ? NaN : (iseven(n) ? 0.5 * (s[n ÷ 2] + s[n ÷ 2 + 1]) : s[(n + 1) ÷ 2])
end

@inline s_phys(J, w) = (J[1]*w[1] + J[4]*w[2] + J[7]*w[3],
                        J[2]*w[1] + J[5]*w[2] + J[8]*w[3],
                        J[3]*w[1] + J[6]*w[2] + J[9]*w[3])

# One Pedrizzetti increment direction (f=1 step): ||G|| Ghat_omegaJ - G.
# Returns the increment vector and the alignment cosine cos(G, omega_J).
@inline function pedrizzetti_increment(G, wJ)
    nG = sqrt(G[1]^2 + G[2]^2 + G[3]^2)
    nw = sqrt(wJ[1]^2 + wJ[2]^2 + wJ[3]^2)
    (nG == 0 || nw == 0) && return (0.0, 0.0, 0.0), NaN
    inc = (nG*wJ[1]/nw - G[1], nG*wJ[2]/nw - G[2], nG*wJ[3]/nw - G[3])
    cosang = (G[1]*wJ[1] + G[2]*wJ[2] + G[3]*wJ[3]) / (nG*nw)
    return inc, clamp(cosang, -1.0, 1.0)
end

# corrected-Pedrizzetti post-rotation Gamma (rlxf step), per FLOWVPM_relaxation.jl.
@inline function corrected_pedrizzetti_gamma(G, wJ, rlxf)
    nG = sqrt(G[1]^2 + G[2]^2 + G[3]^2)
    nw = sqrt(wJ[1]^2 + wJ[2]^2 + wJ[3]^2)
    (nG == 0 || nw == 0) && return (G[1], G[2], G[3])
    cosang = (G[1]*wJ[1] + G[2]*wJ[2] + G[3]*wJ[3]) / (nG*nw)
    b2 = 1 - 2*(1-rlxf)*rlxf*(1 - cosang)
    g1 = (1-rlxf)*G[1] + rlxf*nG*wJ[1]/nw
    g2 = (1-rlxf)*G[2] + rlxf*nG*wJ[2]/nw
    g3 = (1-rlxf)*G[3] + rlxf*nG*wJ[3]/nw
    s = b2 > 0 ? sqrt(b2) : 1.0
    return (g1/s, g2/s, g3/s)
end

function select_indices(np::Int, nsample::Int, seed::Int)
    (nsample <= 0 || nsample >= np) && return collect(1:np)
    rng = Random.MersenneTwister(seed)
    idx = Random.randperm(rng, np)[1:nsample]; sort!(idx); return idx
end

function collect_neighbors(g::POR.CSRGrid, X, k, cutoff2)
    kx, ky, kz = g.key[k]; xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
    ids = Int[]; d2s = Float64[]
    @inbounds for dz in -1:1, dy in -1:1, dx in -1:1
        c = get(g.cellid, (kx+dx, ky+dy, kz+dz), 0); c == 0 && continue
        for p in (g.start[c]+1):g.start[c+1]
            i = g.perm[p]
            rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
            d2 = rx*rx+ry*ry+rz*rz; d2 <= cutoff2 || continue
            push!(ids, i); push!(d2s, d2)
        end
    end
    return ids, d2s
end

# --- framing C: local M2 overlap target Gamma_eff at sampled centers ----------
function local_gamma_eff_alignment(state, g::POR.CSRGrid, cutoff2;
                                   nsample::Int=0, max_neighbors::Int=160, rng_seed::Int=1)
    nan = NaN
    nsample <= 0 && return (; c_n=0, c_cos_dpedr_deff=nan, c_frac_captured=nan,
                            c_rel_pedr=nan, c_rel_raw=nan)
    X, sigma, gamma, J, np = state.X, state.sigma, state.gamma, state.J, state.np
    centers = select_indices(np, nsample, rng_seed)
    cossum = 0.0; fracsum = 0.0; relp2 = 0.0; relr2 = 0.0; den2 = 0.0; used = 0
    for k in centers
        ids, d2s = collect_neighbors(g, X, k, cutoff2); isempty(ids) && continue
        order = sortperm(d2s); loc = ids[order[1:min(length(order), max_neighbors)]]
        (k in loc) || (loc[end] = k); center = findfirst(==(k), loc)
        center === nothing && continue; nloc = length(loc)
        M = zeros(nloc, nloc); W = zeros(3, nloc); G = gamma[:, loc]
        @inbounds for a in 1:nloc, b in 1:nloc
            ib = loc[b]
            rx = X[1,loc[a]]-X[1,ib]; ry = X[2,loc[a]]-X[2,ib]; rz = X[3,loc[a]]-X[3,ib]
            M[a,b] = POR.zeta_sigma(rx*rx+ry*ry+rz*rz, sigma[ib])
        end
        W .= G * transpose(M)
        @inbounds for a in 1:nloc
            ka = loc[a]; h1=0.0; h2=0.0; h3=0.0
            for ib in loc
                rx = X[1,ka]-X[1,ib]; ry = X[2,ka]-X[2,ib]; rz = X[3,ka]-X[3,ib]
                q1,q2,q3 = hess_green_times_gamma_components(rx,ry,rz,sigma[ib],
                    gamma[1,ib],gamma[2,ib],gamma[3,ib])
                h1+=q1; h2+=q2; h3+=q3
            end
            W[1,a]+=h1; W[2,a]+=h2; W[3,a]+=h3
        end
        Geff = permutedims(M \ transpose(W))
        Gc = G[:, center]; wJ = POR.curl_from_J(state)[:, k]
        d_eff = Geff[:, center] .- Gc
        inc, _ = pedrizzetti_increment(Gc, wJ); d_pedr = collect(inc)
        nde = norm(d_eff); ndp = norm(d_pedr)
        if nde > 0 && ndp > 0
            cossum += dot(d_eff, d_pedr)/(nde*ndp)
            fracsum += dot(d_eff, d_pedr)/(nde*ndp)   # signed fraction of d_eff along d_pedr direction
        end
        # how close does ONE rlxf=1 rotation land vs raw, relative to Geff
        Gpedr = Gc .+ d_pedr
        relp2 += norm(Gpedr .- Geff[:,center])^2; relr2 += norm(Gc .- Geff[:,center])^2
        den2 += norm(Geff[:,center])^2; used += 1
    end
    used == 0 && return (; c_n=0, c_cos_dpedr_deff=nan, c_frac_captured=nan,
                         c_rel_pedr=nan, c_rel_raw=nan)
    return (; c_n=used, c_cos_dpedr_deff=cossum/used, c_frac_captured=fracsum/used,
            c_rel_pedr=sqrt(relp2/max(den2,eps())), c_rel_raw=sqrt(relr2/max(den2,eps())))
end

# --- effective overlap (honest, replaces the misleading 4*sig_max cutoff count) -
# The raw neighbor count is taken within cutoff = 4*sig_max (the GLOBAL max core
# size), which spans ~1.7 R here and engulfs ~half the wake. This samples the TRUE
# local overlap instead: nearest-neighbor spacing h vs the LOCAL sigma, the count
# within 2*sigma_local, and the kernel-weighted effective neighbor count
# (off-diagonal zeta row-sum / diagonal) -- the number that actually drives the
# overlap operator.
function eff_overlap_sample(state, g::POR.CSRGrid, cutoff2; nsample::Int=500, rng_seed::Int=1)
    X, sigma, np = state.X, state.sigma, state.np
    idx = select_indices(np, nsample, rng_seed)
    dnn = Float64[]; soh = Float64[]; n2 = Int[]; ncut = Int[]; neff = Float64[]
    for k in idx
        ids, d2s = collect_neighbors(g, X, k, cutoff2)
        sk = sigma[k]; diag = POR.zeta_sigma(0.0, sk); two_sk2 = (2*sk)^2
        mind2 = Inf; rs = 0.0; cnt2 = 0
        @inbounds for (j, i) in enumerate(ids)
            d2 = d2s[j]
            i != k && d2 < mind2 && (mind2 = d2)
            rs += POR.zeta_sigma(d2, sigma[i])
            d2 <= two_sk2 && (cnt2 += 1)
        end
        isfinite(mind2) || continue
        h = sqrt(mind2)
        push!(dnn, h); push!(soh, sk/h)
        push!(n2, cnt2 - 1)              # exclude self
        push!(ncut, length(ids) - 1)
        push!(neff, (rs - diag)/max(diag, eps()))
    end
    return (; eff_n=length(dnn), h_local_med=_median(dnn),
            sigma_over_h_med=_median(soh), n_within_2sig_med=_median(n2),
            n_within_cutoff_med=_median(ncut), n_eff_kernel_med=_median(neff))
end

# --- per-state analysis ------------------------------------------------------
function analyze_pedrizzetti(state; cutoff_factor::Float64=4.0, rlxf::Float64=0.3,
                             rate_factor::Float64=1.0, target_nsample::Int=0,
                             local_nsample::Int=20, local_max_neighbors::Int=160,
                             eff_nsample::Int=500, rng_seed::Int=1)
    np = state.np; np == 0 && return nothing
    X, sigma, gamma, J = state.X, state.sigma, state.gamma, state.J
    sig_max = maximum(sigma); cutoff = cutoff_factor*sig_max; cutoff2 = cutoff*cutoff
    t0 = time()
    wJ = POR.curl_from_J(state)                       # = omega_J, what Pedrizzetti reads
    z0(i) = POR.zeta_sigma(0.0, sigma[i])

    # ---- Framing B: per-particle alignment + ∥/⊥ overlap-correction split ----
    angs = Float64[]; magfrac = Float64[]; wts = Float64[]
    bidx = select_indices(np, target_nsample, rng_seed)
    @inbounds for i in bidx
        G = view(gamma, :, i); w = view(wJ, :, i)
        nG = norm(G); nw = norm(w)
        (nG == 0 || nw == 0) && continue
        cosang = clamp(dot(G, w)/(nG*nw), -1.0, 1.0)
        push!(angs, acos(cosang)); push!(wts, nG)
        # overlap correction = omega_J - isolated self-term (2/3)zeta(0) G  (∥ G)
        self = (2/3)*z0(i)
        corr = (w[1]-self*G[1], w[2]-self*G[2], w[3]-self*G[3])
        ncorr = sqrt(corr[1]^2+corr[2]^2+corr[3]^2)
        if ncorr > 0
            ghat = (G[1]/nG, G[2]/nG, G[3]/nG)
            par = corr[1]*ghat[1] + corr[2]*ghat[2] + corr[3]*ghat[3]   # ∥ component (magnitude)
            push!(magfrac, abs(par)/ncorr)
        end
    end
    wsum = sum(wts)
    ang_med = _median(angs); ang_mean = isempty(angs) ? NaN : mean(angs)
    ang_p90 = isempty(angs) ? NaN : quantile(angs, 0.9)
    ang_wmean = wsum > 0 ? sum(angs .* wts)/wsum : NaN
    magfrac_med = _median(magfrac); magfrac_mean = isempty(magfrac) ? NaN : mean(magfrac)

    # ---- Framing A: rate-residual on the full cloud (one fused neighbor pass) -
    # candidate rates (3xnp):  self = J G ;  pedr = self + rate*(||G|| what_J - G)
    Vself = zeros(3, np); Vpedr = zeros(3, np); Vpedrc = zeros(3, np)
    @inbounds for i in 1:np
        Ji = view(J, :, i); Gi = view(gamma, :, i)
        sp = s_phys(Ji, Gi); Vself[1,i],Vself[2,i],Vself[3,i] = sp
        inc, _ = pedrizzetti_increment(Gi, view(wJ, :, i))
        Vpedr[1,i] = sp[1] + rate_factor*inc[1]
        Vpedr[2,i] = sp[2] + rate_factor*inc[2]
        Vpedr[3,i] = sp[3] + rate_factor*inc[3]
        gc = corrected_pedrizzetti_gamma(Gi, view(wJ, :, i), rlxf)
        Vpedrc[1,i] = sp[1] + (rate_factor/rlxf)*(gc[1]-Gi[1])    # /rlxf: rate_factor already carries rlxf
        Vpedrc[2,i] = sp[2] + (rate_factor/rlxf)*(gc[2]-Gi[2])
        Vpedrc[3,i] = sp[3] + (rate_factor/rlxf)*(gc[3]-Gi[3])
    end
    g = POR.build_csr(X, np, cutoff)
    eff = eff_overlap_sample(state, g, cutoff2; nsample=eff_nsample, rng_seed=rng_seed+20_000)
    # omega = M G ; conv ; M*Vself ; M*Vpedr ; M*Vpedrc  in one pass
    omega = zeros(3,np); conv = zeros(3,np)
    Ms = zeros(3,np); Mp = zeros(3,np); Mpc = zeros(3,np)
    U = state.velocity
    @inbounds for k in 1:np
        kx,ky,kz = g.key[k]; xk1,xk2,xk3 = X[1,k],X[2,k],X[3,k]
        uk1,uk2,uk3 = U[1,k],U[2,k],U[3,k]
        o1=0.0;o2=0.0;o3=0.0; c1=0.0;c2=0.0;c3=0.0
        s1=0.0;s2=0.0;s3=0.0; p1=0.0;p2=0.0;p3=0.0; q1=0.0;q2=0.0;q3=0.0
        for dz in -1:1, dy in -1:1, dx in -1:1
            cc = get(g.cellid, (kx+dx,ky+dy,kz+dz), 0); cc == 0 && continue
            for pp in (g.start[cc]+1):g.start[cc+1]
                i = g.perm[pp]
                rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
                d2 = rx*rx+ry*ry+rz*rz; d2 > cutoff2 && continue
                si = sigma[i]; si2 = si*si
                z = POR.ZETA0*exp(-d2/(2*si2))/(si*si2)
                o1+=z*gamma[1,i]; o2+=z*gamma[2,i]; o3+=z*gamma[3,i]
                s1+=z*Vself[1,i]; s2+=z*Vself[2,i]; s3+=z*Vself[3,i]
                p1+=z*Vpedr[1,i]; p2+=z*Vpedr[2,i]; p3+=z*Vpedr[3,i]
                q1+=z*Vpedrc[1,i]; q2+=z*Vpedrc[2,i]; q3+=z*Vpedrc[3,i]
                if i != k
                    w = -z*((uk1-U[1,i])*(-rx)+(uk2-U[2,i])*(-ry)+(uk3-U[3,i])*(-rz))/si2
                    c1+=gamma[1,i]*w; c2+=gamma[2,i]*w; c3+=gamma[3,i]*w
                end
            end
        end
        omega[1,k]=o1;omega[2,k]=o2;omega[3,k]=o3; conv[1,k]=c1;conv[2,k]=c2;conv[3,k]=c3
        Ms[1,k]=s1;Ms[2,k]=s2;Ms[3,k]=s3; Mp[1,k]=p1;Mp[2,k]=p2;Mp[3,k]=p3
        Mpc[1,k]=q1;Mpc[2,k]=q2;Mpc[3,k]=q3
    end
    # physical RHS b = J*omega - conv  (reduced = no conv)
    b_red = similar(omega)
    @inbounds for k in 1:np
        sp = s_phys(view(J,:,k), view(omega,:,k))
        b_red[1,k]=sp[1]; b_red[2,k]=sp[2]; b_red[3,k]=sp[3]
    end
    b_smp = b_red .- conv
    r_self_red = POR.relnorm(Ms .- b_red, b_red)
    r_pedr_red = POR.relnorm(Mp .- b_red, b_red)
    r_pedrc_red = POR.relnorm(Mpc .- b_red, b_red)
    r_self_smp = POR.relnorm(Ms .- b_smp, b_smp)
    r_pedr_smp = POR.relnorm(Mp .- b_smp, b_smp)

    # ---- Framing C: alignment vs local M2 Gamma_eff --------------------------
    cres = local_gamma_eff_alignment(state, g, cutoff2; nsample=local_nsample,
        max_neighbors=local_max_neighbors, rng_seed=rng_seed+10_000)

    return (; np, n_angle=length(angs),
            ang_med, ang_mean, ang_p90, ang_wmean, magfrac_med, magfrac_mean,
            r_self_red, r_pedr_red, r_pedrc_red, ratio_red=r_pedr_red/max(r_self_red,eps()),
            ratio_corr=r_pedrc_red/max(r_self_red,eps()),
            r_self_smp, r_pedr_smp, ratio_smp=r_pedr_smp/max(r_self_smp,eps()),
            eff..., cres..., elapsed=time()-t0)
end

# --- synthetic gates ---------------------------------------------------------
function synthetic_gates()
    println("\n[gate 1] isolated particle, omega_J ∥ Gamma: Pedrizzetti idle")
    G = [1.0, 0.2, -0.1]
    # build J so curl(J) = (J6-J8, J7-J3, J2-J4) = G  (parallel)
    J = zeros(9); J[2]=G[3]; J[6]=G[1]; J[7]=G[2]
    st = (; X=reshape([0.0,0.0,0.0],3,1), gamma=reshape(G,3,1), sigma=[1.0],
            velocity=zeros(3,1), vorticity=zeros(3,1), J=reshape(J,9,1), np=1)
    r = analyze_pedrizzetti(st; cutoff_factor=4.0, target_nsample=0, local_nsample=0)
    @printf("    angle_med=%.3e r_self_red=%.3e r_pedr_red=%.3e (both ~0: self exact, Pedrizzetti idle)\n",
            r.ang_med, r.r_self_red, r.r_pedr_red)
    @assert r.ang_med < 1e-6          # acos near 1 has a ~1e-8 roundoff floor
    @assert r.r_self_red < 1e-9       # isolated: self-only stretch exactly satisfies M Gdot = b
    @assert r.r_pedr_red < 1e-9       # Pedrizzetti increment is zero (omega_J ∥ Gamma)
    println("    PASS")

    println("\n[gate 2] two particles, non-parallel omega_J: angle>0, Pedrizzetti active")
    X2 = [0.0 0.8; 0.0 0.3; 0.0 -0.2]; G2 = [1.0 -0.4; 0.2 0.7; 0.0 0.5]
    J2 = hcat([0.0,0.0,0.0,0.0,0.0,0.3,0.1,0.0,0.0],
              [0.0,0.2,0.0,0.0,0.0,-0.1,0.4,0.0,0.0])
    st2 = (; X=X2, gamma=G2, sigma=[1.0,1.0], velocity=zeros(3,2), vorticity=zeros(3,2),
             J=J2, np=2)
    r2 = analyze_pedrizzetti(st2; cutoff_factor=4.0, target_nsample=0,
                             local_nsample=2, local_max_neighbors=2)
    @printf("    angle_med=%.3e ratio_red=%.4f magfrac_med=%.3f c_cos=%.3f\n",
            r2.ang_med, r2.ratio_red, r2.magfrac_med, r2.c_cos_dpedr_deff)
    @assert r2.ang_med > 0
    @assert r2.c_n == 2
    println("    PASS")
end

# --- driver ------------------------------------------------------------------
run_name = length(ARGS) >= 1 ? ARGS[1] : get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
steps_arg = length(ARGS) >= 2 ? ARGS[2] : get(ENV, "STEPS", "")
cutoff_factor = parse(Float64, get(ENV, "CUTOFF_FACTOR", "4.0"))
nstates = parse(Int, get(ENV, "NSTATES", "20"))
target_nsample = parse(Int, get(ENV, "TARGET_NSAMPLE", "0"))      # 0 = all particles for framing B
local_nsample = parse(Int, get(ENV, "LOCAL_SOLVE_NSAMPLE", "20"))
local_max_neighbors = parse(Int, get(ENV, "LOCAL_MAX_NEIGHBORS", "160"))
eff_nsample = parse(Int, get(ENV, "EFF_NSAMPLE", "500"))
rng_seed = parse(Int, get(ENV, "RNG_SEED", "1"))
rlxf = parse(Float64, get(ENV, "RLXF", "0.3"))
nsteps_relax = parse(Int, get(ENV, "NSTEPS_RELAX", "1"))
RPM = parse(Float64, get(ENV, "RPM", "6000"))
NT = parse(Int, get(ENV, "NT", "36"))
dt = haskey(ENV, "DT") ? parse(Float64, ENV["DT"]) : 60/RPM/NT
rate_factor = rlxf / (nsteps_relax * dt)

data_root = joinpath(@__DIR__, "..", "data")
part_dir = joinpath(data_root, run_name, "$(run_name)_wake1_particles")
isdir(part_dir) || error("Particle dir not found: $(part_dir)")
vtp_path(idx) = joinpath(part_dir, "$(run_name)_wake1_particles.$(idx).vtp")
avail = Int[]
for f in readdir(part_dir)
    m = match(Regex("$(run_name)_wake1_particles\\.(\\d+)\\.vtp\$"), f)
    m === nothing || push!(avail, parse(Int, m.captures[1]))
end
sort!(avail); isempty(avail) && error("No particle VTP files in $(part_dir)")
steps = isempty(steps_arg) ? avail[max(1,end-nstates+1):end] :
        (occursin(':', steps_arg) ?
            (lh = parse.(Int, split(steps_arg, ':')); filter(s -> lh[1] <= s <= lh[2], avail)) :
            [parse(Int, steps_arg)])
isempty(steps) && error("No steps selected from $(length(avail)) available.")

println("="^78)
println("Pedrizzetti vs overlap-residual diagnostic (item 011)")
println("  run            = $(run_name)")
println("  steps          = $(first(steps))...$(last(steps)) ($(length(steps)) states; $(length(avail)) available)")
println("  cutoff_fac     = $(cutoff_factor)")
println("  rlxf           = $(rlxf)  nsteps_relax = $(nsteps_relax)  dt = $(dt)")
println("  rate_factor    = $(rate_factor)  (= rlxf/(nsteps_relax*dt))")
println("  target_nsample = $(target_nsample <= 0 ? "all" : string(target_nsample)) (framing B)")
println("  local_nsample  = $(local_nsample) (framing C, max_nbr=$(local_max_neighbors))")
println("  erf source     = ", HAVE_SF ? "SpecialFunctions.erf" : "self-contained series")
println("="^78)

synthetic_gates()

csv_path = joinpath(data_root, run_name, "particle_pedrizzetti_overlap_diag.csv")
header = ["step","np","n_angle","ang_med","ang_mean","ang_p90","ang_wmean",
          "magfrac_med","magfrac_mean","r_self_red","r_pedr_red","r_pedrc_red",
          "ratio_red","ratio_corr","r_self_smp","r_pedr_smp","ratio_smp",
          "eff_n","h_local_med","sigma_over_h_med","n_within_2sig_med",
          "n_within_cutoff_med","n_eff_kernel_med",
          "c_n","c_cos_dpedr_deff","c_frac_captured","c_rel_pedr","c_rel_raw","elapsed_s"]
rows = NamedTuple[]
println("\n[diagnostic] per-state metrics (angles in rad)")
@printf("  %5s %8s %8s %9s %9s %9s %9s %8s\n",
        "step","ang_med","ang_p90","magf_med","r_self","r_pedr","ratio","c_cos")
for s in steps
    st = POR.load_state(vtp_path(s))
    res = analyze_pedrizzetti(st; cutoff_factor, rlxf, rate_factor, target_nsample,
                              local_nsample, local_max_neighbors, eff_nsample, rng_seed)
    res === nothing && (println("  step $s: 0 particles"); continue)
    push!(rows, (; step=s, res...))
    @printf("  %5d %8.4f %8.4f %9.3f %9.3e %9.3e %9.3f %8.3f\n",
            s, res.ang_med, res.ang_p90, res.magfrac_med,
            res.r_self_red, res.r_pedr_red, res.ratio_red, res.c_cos_dpedr_deff)
end

open(csv_path, "w") do io
    println(io, join(header, ","))
    for r in rows
        println(io, join((r.step, r.np, r.n_angle, r.ang_med, r.ang_mean, r.ang_p90,
            r.ang_wmean, r.magfrac_med, r.magfrac_mean, r.r_self_red, r.r_pedr_red,
            r.r_pedrc_red, r.ratio_red, r.ratio_corr, r.r_self_smp, r.r_pedr_smp,
            r.ratio_smp, r.eff_n, r.h_local_med, r.sigma_over_h_med,
            r.n_within_2sig_med, r.n_within_cutoff_med, r.n_eff_kernel_med,
            r.c_n, r.c_cos_dpedr_deff, r.c_frac_captured, r.c_rel_pedr,
            r.c_rel_raw, r.elapsed), ","))
    end
end

if !isempty(rows)
    mf(f) = mean(getfield.(rows, f))
    println("\n" * "="^78)
    println("SUMMARY over $(length(rows)) states")
    @printf("  [B] angle ∠(Gamma,omega_J) rad : median %.4f  mean %.4f  p90 %.4f  Γ-wt mean %.4f\n",
            mf(:ang_med), mf(:ang_mean), mf(:ang_p90), mf(:ang_wmean))
    @printf("      (median in deg = %.2f)\n", rad2deg(mf(:ang_med)))
    @printf("  [B] overlap-corr magnitude(∥) fraction : median %.3f  mean %.3f (Pedrizzetti CANNOT reach this)\n",
            mf(:magfrac_med), mf(:magfrac_mean))
    @printf("  [A] residual self / pedr / corrected   : %.3e / %.3e / %.3e\n",
            mf(:r_self_red), mf(:r_pedr_red), mf(:r_pedrc_red))
    @printf("  [A] residual ratio pedr/self (red,smp) : %.3f / %.3f   corrected/self %.3f\n",
            mf(:ratio_red), mf(:ratio_smp), mf(:ratio_corr))
    @printf("  [overlap] TRUE local σ/h median %.2f  (nearest-nbr spacing h_local %.4f)\n",
            mf(:sigma_over_h_med), mf(:h_local_med))
    @printf("  [overlap] neighbors: within 2σ_local %.0f | within 4σ_max cutoff %.0f | kernel-effective %.0f\n",
            mf(:n_within_2sig_med), mf(:n_within_cutoff_med), mf(:n_eff_kernel_med))
    @printf("  [C] cos(d_pedr, d_eff) vs local Γ_eff  : %.3f   frac captured %.3f\n",
            mf(:c_cos_dpedr_deff), mf(:c_frac_captured))
    @printf("  [C] ‖Γ_pedr-Γ_eff‖/‖Γ_eff‖ vs raw      : %.3f / %.3f\n",
            mf(:c_rel_pedr), mf(:c_rel_raw))
    println("  CSV: $(csv_path)")
    println("="^78)
end
