# p018 sigma-collapse forensics (Track F) -- direct-strain companion to
# scripts/p018_sigma_forensics.py. Computes the particle-induced velocity
# gradient J at selected targets with the same Gaussian-erf regularized VPM
# kernel FLOWVPM uses (g_gauserf via FLOWVPM.g_dgdr_gauserf), then the axial
# strain along each target's own strength axis, Z = Ghat . S . Ghat with
# S = (J + J')/2. This is the raw f=0,g=1 MM4 quantity (the model applies a
# further 1/5 factor; we report raw Z).
#
# Gates:
#   (a) healthy-719 aged column (ax/R > 1.3): median dt*Z should be ~ +1e-3
#   (b) single-particle analytic J vs central finite differences, rel err < 1e-6
# Feedback test: recompute Z at the collapsing targets EXCLUDING thin sources
# (sigma/sigma_shed < 0.25); a large drop => thin cores strain each other.
#
# NOTE the particle self-term: for any radial kernel the self-induced J at the
# particle center is purely antisymmetric ((c/sigma^3)[Gamma x]), so skipping
# the self interaction is EXACT for the symmetric part S and hence for Z.
#
# Run: JULIA_NUM_THREADS=6 julia --project scripts/p018_sigma_forensics_strain.jl

import ReadVTK
import FLOWVPM
import Statistics: median, quantile
import Random
import Printf: @sprintf
using Base.Threads: @threads, nthreads
using LinearAlgebra: norm, cross

const REPO = dirname(@__DIR__)
const OUTDIR = joinpath(REPO, "data", "p018_L2_visc_forensics")
const CORPSE_VTP = joinpath(OUTDIR, "p018_L2_visc_wake1_particles.1041.vtp")
const HEALTHY_VTP = expanduser(
    "~/p018_L1_ov3_paraview/p018_L1_ov3_wake1_particles/p018_L1_ov3_wake1_particles.719.vtp")

const R = 0.119
const DT = 60.0 / (5400 * 36)               # both runs
const SIGMA_SHED_CORPSE = 2pi * R * 2.88 / (36 * 26)
const SIGMA_SHED_HEALTHY = 4.45e-3
const THIN = 0.25                           # sigma/sigma_shed collapse threshold

# ---------------------------------------------------------------------------
function load_particles(path)
    vtk = ReadVTK.VTKFile(path)
    np = vtk.n_points
    pd = ReadVTK.get_point_data(vtk)
    X = Matrix{Float64}(ReadVTK.get_points(vtk))
    G = Matrix{Float64}(reshape(ReadVTK.get_data(pd["gamma"]), 3, np))
    sig = Vector{Float64}(vec(ReadVTK.get_data(pd["sigma"])))
    J = Matrix{Float64}(reshape(ReadVTK.get_data(pd["velocity_gradient"]), 9, np))
    keep = [norm(view(G, :, i)) > 0 for i in 1:np]
    (; X=X[:, keep], G=G[:, keep], sig=sig[keep], Jrec=J[:, keep], n=count(keep))
end

"""Velocity and Jacobian at point x induced by one regularized particle:
u = (1/4pi) * g(rho)/r^3 * (Gamma x r),  r = x - xp, rho = r/sigma.
J_ij = du_i/dx_j = (1/4pi)*[ qp*(r_j/r)*(Gamma x r)_i + q*(Gamma x e_j)_i ]
with q = g/r^3, qp = dgdr/(sigma*r^3) - 3g/r^4."""
@inline function uJ_single!(J::AbstractMatrix, x1, x2, x3, p1, p2, p3,
                            G1, G2, G3, sigma)
    r1 = x1 - p1; r2 = x2 - p2; r3 = x3 - p3
    r2sum = r1*r1 + r2*r2 + r3*r3
    r = sqrt(r2sum)
    r < 1e-14 && return (0.0, 0.0, 0.0)
    g, dgdr = FLOWVPM.g_dgdr_gauserf(r / sigma)
    q = g / (r2sum * r)
    qp = dgdr / (sigma * r2sum * r) - 3g / (r2sum * r2sum)
    c1 = G2*r3 - G3*r2; c2 = G3*r1 - G1*r3; c3 = G1*r2 - G2*r1  # Gamma x r
    f = 1 / (4pi)
    @inbounds for j in 1:3
        rj = (j == 1 ? r1 : j == 2 ? r2 : r3)
        a = f * qp * rj / r
        J[1, j] += a * c1
        J[2, j] += a * c2
        J[3, j] += a * c3
    end
    # q * (Gamma x e_j)_i = f*q*eps_{iaj}*G_a terms (purely antisymmetric)
    fq = f * q
    J[2, 1] += fq * G3;  J[3, 1] -= fq * G2
    J[1, 2] -= fq * G3;  J[3, 2] += fq * G1
    J[1, 3] += fq * G2;  J[2, 3] -= fq * G1
    return (f * q * c1, f * q * c2, f * q * c3)
end

"J (3x3xN) at targets Xt from all sources (X,G,sig) with source mask `use`."
function compute_J(Xt, X, G, sig, use::AbstractVector{Bool})
    n = size(Xt, 2); np = size(X, 2)
    Jout = zeros(3, 3, n)
    @threads for t in 1:n
        Jt = view(Jout, :, :, t)
        x1, x2, x3 = Xt[1, t], Xt[2, t], Xt[3, t]
        @inbounds for p in 1:np
            use[p] || continue
            uJ_single!(Jt, x1, x2, x3, X[1, p], X[2, p], X[3, p],
                       G[1, p], G[2, p], G[3, p], sig[p])
        end
    end
    Jout
end

"Z = Ghat . S . Ghat per target."
function axial_strain(J, Gt)
    n = size(J, 3)
    Z = zeros(n)
    for t in 1:n
        g = view(Gt, :, t); gh = g / max(norm(g), 1e-300)
        S = (J[:, :, t] + J[:, :, t]') / 2
        Z[t] = gh' * S * gh
    end
    Z
end

stats(z) = (median(z), quantile(z, 0.95), quantile(z, 0.05))

# ---------------------------------------------------------------------------
function gate_b()
    xp = [0.01, -0.02, 0.005]; Gp = [0.3, -0.1, 0.2]; sp = 2.3e-3
    xt = xp .+ [1.7e-3, -0.9e-3, 2.2e-3]
    X = reshape(xp, 3, 1); G = reshape(Gp, 3, 1)
    Jana = compute_J(reshape(xt, 3, 1), X, G, [sp], [true])[:, :, 1]
    u_at(x) = begin
        Jtmp = zeros(3, 3)
        u = uJ_single!(Jtmp, x[1], x[2], x[3], xp[1], xp[2], xp[3],
                       Gp[1], Gp[2], Gp[3], sp)
        collect(u)
    end
    h = 1e-7
    Jfd = zeros(3, 3)
    for j in 1:3
        e = zeros(3); e[j] = h
        Jfd[:, j] = (u_at(xt .+ e) .- u_at(xt .- e)) ./ (2h)
    end
    relerr = norm(Jana .- Jfd) / norm(Jfd)
    println("GATE(b): analytic vs FD Jacobian rel err = ", @sprintf("%.3e", relerr),
            "  -> ", relerr < 1e-6 ? "PASS" : "FAIL")
    relerr
end

function axcoord(X)  # axial/R with wake positive
    ax = X[1, :] ./ R
    median(ax) < 0 ? -ax : ax
end

function main()
    println("threads = ", nthreads())
    gate_b()

    # ---- gate (a): healthy aged column -----------------------------------
    H = load_particles(HEALTHY_VTP)
    axh = axcoord(H.X)
    rng = Random.MersenneTwister(0)
    aged = findall(axh .> 1.3)
    sub = aged[Random.randperm(rng, length(aged))[1:min(3000, length(aged))]]
    t = @elapsed Jh = compute_J(H.X[:, sub], H.X, H.G, H.sig, trues(H.n))
    Zh = axial_strain(Jh, H.G[:, sub])
    mh, p95h, p5h = stats(DT .* Zh)
    println(@sprintf("GATE(a): healthy719 aged column (ax/R>1.3, N=%d of %d, %.0f s): dt*Z med %+.3e  p95 %+.3e  p5 %+.3e  -> %s (expect ~ +1e-3, order of magnitude)",
                     length(sub), length(aged), t, mh, p95h, p5h,
                     (3e-4 < mh < 3e-3) ? "PASS" : (1e-4 < mh < 1e-2) ? "PASS(loose)" : "FAIL"))
    # convention cross-check vs recorded J
    Zrec = axial_strain(reshape(H.Jrec[:, sub], 3, 3, length(sub)), H.G[:, sub])
    cor_ = sum((Zh .- median(Zh)) .* (Zrec .- median(Zrec))) /
           (norm(Zh .- median(Zh)) * norm(Zrec .- median(Zrec)))
    println(@sprintf("  recomputed-vs-recorded Z correlation (healthy sub): %+.3f (positive => sign convention consistent; recorded J also holds panel/SFS contributions)", cor_))

    # ---- Analysis 4: corpse 1041 -----------------------------------------
    C = load_particles(CORPSE_VTP)
    srel = C.sig ./ SIGMA_SHED_CORPSE
    axc = axcoord(C.X)
    thin = srel .< THIN
    tgt_thin = findall(thin)
    healthy_pool = findall((srel .> 0.5) .& (axc .> 0.0))
    tgt_ref = healthy_pool[Random.randperm(rng, length(healthy_pool))[1:3000]]
    println("\ncorpse 1041: N=", C.n, "  thin(<%0.25 shed)=", length(tgt_thin))

    Jt_all = compute_J(C.X[:, tgt_thin], C.X, C.G, C.sig, trues(C.n))
    Zt_all = axial_strain(Jt_all, C.G[:, tgt_thin])
    Jt_exc = compute_J(C.X[:, tgt_thin], C.X, C.G, C.sig, .!thin)
    Zt_exc = axial_strain(Jt_exc, C.G[:, tgt_thin])
    Jr_all = compute_J(C.X[:, tgt_ref], C.X, C.G, C.sig, trues(C.n))
    Zr_all = axial_strain(Jr_all, C.G[:, tgt_ref])

    println("\n== Analysis 4: direct strain (raw Z = Ghat.S.Ghat; MM4 uses Z/5) ==")
    for (nm, Z) in [("thin targets, ALL sources", Zt_all),
                    ("thin targets, thin sources EXCLUDED", Zt_exc),
                    ("healthy-region ref targets, ALL sources", Zr_all)]
        m, p95, p5 = stats(DT .* Z)
        println(@sprintf("%-40s N=%5d  dt*Z: med %+.3e  p95 %+.3e  p5 %+.3e  med|dtZ| %.3e",
                         nm, length(Z), m, p95, p5, median(abs.(DT .* Z))))
    end
    ratio = median(abs.(Zt_exc)) / median(abs.(Zt_all))
    contraction_all = median(DT .* Zt_all)
    contraction_exc = median(DT .* Zt_exc)
    println(@sprintf("\nfeedback signature: med|Z| ratio (excl/all) = %.3f ; med dt*Z %+.3e -> %+.3e when thin sources removed",
                     ratio, contraction_all, contraction_exc))
    println(ratio < 0.5 ?
        "  => feedback CONFIRMED: thin cores dominate the strain on thin cores" :
        "  => feedback NOT confirmed: ambient (thick-core) strain dominates")

    open(joinpath(OUTDIR, "strain_corpse1041.csv"), "w") do io
        println(io, "index,ax_over_R,sig_rel,dtZ_all,dtZ_excl_thin")
        for (k, i) in enumerate(tgt_thin)
            println(io, @sprintf("%d,%.5f,%.5f,%.6e,%.6e", i, axc[i], srel[i],
                                 DT * Zt_all[k], DT * Zt_exc[k]))
        end
    end
    open(joinpath(OUTDIR, "strain_healthy719_gatea.csv"), "w") do io
        println(io, "index,ax_over_R,dtZ")
        for (k, i) in enumerate(sub)
            println(io, @sprintf("%d,%.5f,%.6e", i, axh[i], DT * Zh[k]))
        end
    end
    println("\nWrote strain_corpse1041.csv, strain_healthy719_gatea.csv in ", OUTDIR)
end

main()
