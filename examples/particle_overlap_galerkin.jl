# =============================================================================
# Galerkin overlap helpers for BRAINSTORM item 008.
#
# This module is intentionally script-local: it reuses the validated VTP loader,
# kernel normalization, CSR grid, and FLOWVPM stretching conventions from
# `particle_overlap_residual.jl`, then adds the Gaussian Gram operator
# M_G[k,i] = ∫ζ_k ζ_i dx = ζ_{sqrt(σ_k²+σ_i²)}(X_k-X_i).
# =============================================================================
module ParticleOverlapGalerkin

include(joinpath(@__DIR__, "particle_overlap_residual.jl"))
const POR = ParticleOverlapResidual

import LinearAlgebra: Symmetric, cholesky, dot, eigvals, norm, svdvals
import Serialization
import Statistics: mean, quantile

const ZETA0 = POR.ZETA0
const DEFAULT_R = 0.119
const EPSN = 1e-30

relnorm(A, B; eps=EPSN) = norm(A) / (norm(B) + eps)

median(v) = begin
    s = sort(collect(v)); n = length(s)
    n == 0 ? NaN : (iseven(n) ? 0.5 * (s[n ÷ 2] + s[n ÷ 2 + 1]) : s[(n + 1) ÷ 2])
end

function particle_path(run_name::AbstractString, step::Integer)
    stem = "$(run_name)_wake1_particles"
    return joinpath("data", run_name, stem, "$(stem).$(step).vtp")
end

@inline gram_zeta(d2::Float64, sk::Float64, si::Float64) =
    POR.zeta_sigma(d2, sqrt(sk*sk + si*si))

@inline gram_diag(si::Float64) = ZETA0 / (2sqrt(2.0) * si^3)

function isolated_rate(state; transposed::Bool=true, include_sfs::Bool=false)
    _, V = POR.isolated_rates(state; transposed)
    if include_sfs
        out = copy(V)
        @inbounds for i in 1:state.np
            fac = state.C[i] * state.sigma[i]^3 / ZETA0
            out[1,i] -= fac * state.SFS[1,i]
            out[2,i] -= fac * state.SFS[2,i]
            out[3,i] -= fac * state.SFS[3,i]
        end
        return out
    end
    return V
end

function sampled_rhs_centers(state; cutoff_factor::Float64=4.0, transposed::Bool=true,
                             include_sfs::Bool=false)
    np = state.np
    X, sigma, gamma, U, Jm = state.X, state.sigma, state.gamma, state.velocity, state.J
    cutoff = cutoff_factor * maximum(sigma)
    c2 = cutoff * cutoff
    g = POR.build_csr(X, np, cutoff)
    omega = zeros(Float64, 3, np)
    conv = zeros(Float64, 3, np)
    @inbounds for k in 1:np
        kx, ky, kz = g.key[k]
        xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
        uk1, uk2, uk3 = U[1,k], U[2,k], U[3,k]
        o1=0.0; o2=0.0; o3=0.0; c1=0.0; c2v=0.0; c3=0.0
        for dz in -1:1, dy in -1:1, dx in -1:1
            cid = get(g.cellid, (kx+dx, ky+dy, kz+dz), 0)
            cid == 0 && continue
            for p in (g.start[cid]+1):g.start[cid+1]
                i = g.perm[p]
                rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
                d2 = rx*rx + ry*ry + rz*rz
                d2 > c2 && continue
                si = sigma[i]; si2 = si*si
                z = ZETA0 * exp(-d2/(2*si2)) / (si*si2)
                o1 += z*gamma[1,i]; o2 += z*gamma[2,i]; o3 += z*gamma[3,i]
                if i != k
                    w = -z * ((uk1-U[1,i])*(-rx) + (uk2-U[2,i])*(-ry) + (uk3-U[3,i])*(-rz)) / si2
                    c1 += gamma[1,i]*w; c2v += gamma[2,i]*w; c3 += gamma[3,i]*w
                end
            end
        end
        omega[1,k]=o1; omega[2,k]=o2; omega[3,k]=o3
        conv[1,k]=c1; conv[2,k]=c2v; conv[3,k]=c3
    end
    rhs = zeros(Float64, 3, np)
    @inbounds for k in 1:np
        s = POR.s_fvpm(view(Jm, :, k), view(omega, :, k); transposed)
        rhs[1,k] = s[1] - conv[1,k]
        rhs[2,k] = s[2] - conv[2,k]
        rhs[3,k] = s[3] - conv[3,k]
        if include_sfs
            fac = state.C[k] * state.sigma[k]^3 / ZETA0
            rhs[1,k] -= fac * state.SFS[1,k]
            rhs[2,k] -= fac * state.SFS[2,k]
            rhs[3,k] -= fac * state.SFS[3,k]
        end
    end
    return rhs
end

function galerkin_rhs_midpoint(state, rhs_centers; cutoff_factor::Float64=4.0)
    np = state.np
    X, sigma = state.X, state.sigma
    vol = quadrature_weights(state; cutoff_factor)
    cutoff = cutoff_factor * maximum(sigma)
    c2 = cutoff * cutoff
    g = POR.build_csr(X, np, cutoff)
    B = zeros(Float64, 3, np)
    @inbounds for k in 1:np
        sk = sigma[k]; sk2 = sk*sk
        kx, ky, kz = g.key[k]
        xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
        b1=0.0; b2=0.0; b3=0.0
        for dz in -1:1, dy in -1:1, dx in -1:1
            cid = get(g.cellid, (kx+dx, ky+dy, kz+dz), 0)
            cid == 0 && continue
            for p in (g.start[cid]+1):g.start[cid+1]
                i = g.perm[p]
                rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
                d2 = rx*rx + ry*ry + rz*rz
                d2 > c2 && continue
                z = ZETA0 * exp(-d2/(2*sk2)) / (sk*sk2)
                w = vol[i] * z
                b1 += w*rhs_centers[1,i]; b2 += w*rhs_centers[2,i]; b3 += w*rhs_centers[3,i]
            end
        end
        B[1,k]=b1; B[2,k]=b2; B[3,k]=b3
    end
    return B
end

function quadrature_weights(state; cutoff_factor::Float64=4.0)
    if any(!iszero, state.vol)
        return state.vol
    end
    np = state.np
    X, sigma = state.X, state.sigma
    cutoff = cutoff_factor * maximum(sigma)
    c2 = cutoff * cutoff
    g = POR.build_csr(X, np, cutoff)
    w = zeros(Float64, np)
    @inbounds for k in 1:np
        sk = sigma[k]; sk2 = sk*sk
        kx, ky, kz = g.key[k]
        xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
        rowsum = 0.0
        for dz in -1:1, dy in -1:1, dx in -1:1
            cid = get(g.cellid, (kx+dx, ky+dy, kz+dz), 0)
            cid == 0 && continue
            for p in (g.start[cid]+1):g.start[cid+1]
                i = g.perm[p]
                rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
                d2 = rx*rx + ry*ry + rz*rz
                d2 > c2 && continue
                rowsum += ZETA0 * exp(-d2/(2*sk2)) / (sk*sk2)
            end
        end
        w[k] = rowsum > 0 ? 1.0 / rowsum : 0.0
    end
    return w
end

function build_rhs(state; cutoff_factor::Float64=4.0, transposed::Bool=true,
                   include_sfs::Bool=false)
    centers = sampled_rhs_centers(state; cutoff_factor, transposed, include_sfs)
    B = galerkin_rhs_midpoint(state, centers; cutoff_factor)
    return (; centers, B)
end

function galerkin_matvec!(out, state, V; cutoff_factor::Float64=4.0, lambda::Float64=0.0,
                          indices::Union{Nothing,Vector{Int}}=nothing,
                          sources::Union{Nothing,Vector{Int}}=nothing)
    X, sigma = state.X, state.sigma
    rows = indices === nothing ? collect(1:state.np) : indices
    cols = sources === nothing ? rows : sources
    fill!(out, 0.0)
    sigmax = maximum(@view sigma[cols])
    cutoff = cutoff_factor * sqrt(2.0) * sigmax
    c2 = cutoff * cutoff
    # Dense row/column subsets are small in Part 1; full-scale uses all rows.
    if length(cols) == state.np && length(rows) == state.np
        g = POR.build_csr(X, state.np, cutoff)
        @inbounds for k in rows
            kx, ky, kz = g.key[k]
            xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
            sk = sigma[k]
            a1=0.0; a2=0.0; a3=0.0
            for dz in -1:1, dy in -1:1, dx in -1:1
                cid = get(g.cellid, (kx+dx, ky+dy, kz+dz), 0)
                cid == 0 && continue
                for p in (g.start[cid]+1):g.start[cid+1]
                    i = g.perm[p]
                    rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
                    d2 = rx*rx + ry*ry + rz*rz
                    d2 > c2 && continue
                    z = gram_zeta(d2, sk, sigma[i])
                    a1 += z*V[1,i]; a2 += z*V[2,i]; a3 += z*V[3,i]
                end
            end
            if lambda != 0.0
                d = lambda * gram_diag(sk)
                a1 += d*V[1,k]; a2 += d*V[2,k]; a3 += d*V[3,k]
            end
            out[1,k]=a1; out[2,k]=a2; out[3,k]=a3
        end
    else
        colset = Set(cols)
        @inbounds for (rr, k) in enumerate(rows)
            xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
            sk = sigma[k]
            a1=0.0; a2=0.0; a3=0.0
            for i in cols
                rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
                d2 = rx*rx + ry*ry + rz*rz
                d2 > c2 && continue
                z = gram_zeta(d2, sk, sigma[i])
                a1 += z*V[1,i]; a2 += z*V[2,i]; a3 += z*V[3,i]
            end
            if lambda != 0.0 && k in colset
                d = lambda * gram_diag(sk)
                a1 += d*V[1,k]; a2 += d*V[2,k]; a3 += d*V[3,k]
            end
            out[1,k]=a1; out[2,k]=a2; out[3,k]=a3
        end
    end
    return out
end

function dense_gram(state, idx::Vector{Int}; cutoff_factor::Float64=Inf, lambda::Float64=0.0,
                    median_sigma::Bool=false)
    n = length(idx)
    X, sigma = state.X, state.sigma
    sig = median_sigma ? fill(median(sigma[idx]), n) : sigma[idx]
    M = Matrix{Float64}(undef, n, n)
    c2 = isfinite(cutoff_factor) ? (cutoff_factor * sqrt(2.0) * maximum(sig))^2 : Inf
    @inbounds for a in 1:n
        ia = idx[a]
        M[a,a] = gram_diag(sig[a]) * (1 + lambda)
        for b in (a+1):n
            ib = idx[b]
            rx = X[1,ib]-X[1,ia]; ry = X[2,ib]-X[2,ia]; rz = X[3,ib]-X[3,ia]
            d2 = rx*rx + ry*ry + rz*rz
            z = d2 <= c2 ? gram_zeta(d2, sig[a], sig[b]) : 0.0
            M[a,b] = z; M[b,a] = z
        end
    end
    return M
end

function dense_collocation(state, idx::Vector{Int})
    n = length(idx)
    X, sigma = state.X, state.sigma
    M = Matrix{Float64}(undef, n, n)
    @inbounds for a in 1:n, b in 1:n
        ia = idx[a]; ib = idx[b]
        rx = X[1,ib]-X[1,ia]; ry = X[2,ib]-X[2,ia]; rz = X[3,ib]-X[3,ia]
        M[a,b] = POR.zeta_sigma(rx*rx + ry*ry + rz*rz, sigma[ib])
    end
    return M
end

function jacobi_pcg(state, B; cutoff_factor::Float64=4.0, lambda::Float64=0.0,
                    tol::Float64=1e-6, maxiter::Int=150, x0=nothing)
    np = state.np
    X = x0 === nothing ? zeros(Float64, 3, np) : copy(x0)
    AX = zeros(Float64, 3, np)
    galerkin_matvec!(AX, state, X; cutoff_factor, lambda)
    R = B .- AX
    Minv = [1.0 / ((1 + lambda) * gram_diag(state.sigma[i])) for i in 1:np]
    Z = similar(R); P = similar(R); AP = similar(R)
    @inbounds for i in 1:np
        Z[:,i] .= Minv[i] .* R[:,i]
        P[:,i] .= Z[:,i]
    end
    rz = sum(R .* Z)
    bnorm = norm(B) + EPSN
    hist = Float64[norm(R) / bnorm]
    monotone = true
    for it in 1:maxiter
        galerkin_matvec!(AP, state, P; cutoff_factor, lambda)
        denom = sum(P .* AP)
        denom <= 0 && return (; x=X, history=hist, iterations=it-1, converged=false, monotone=false, breakdown=:nonpositive_curvature)
        alpha = rz / denom
        @. X = X + alpha * P
        @. R = R - alpha * AP
        push!(hist, norm(R) / bnorm)
        monotone &= hist[end] <= hist[end-1] * (1 + 1e-10)
        hist[end] <= tol && return (; x=X, history=hist, iterations=it, converged=true, monotone, breakdown=:none)
        @inbounds for i in 1:np
            Z[:,i] .= Minv[i] .* R[:,i]
        end
        rz_new = sum(R .* Z)
        beta = rz_new / rz
        @. P = Z + beta * P
        rz = rz_new
    end
    return (; x=X, history=hist, iterations=maxiter, converged=false, monotone, breakdown=:maxiter)
end

function residual_galerkin(state, V, B; cutoff_factor::Float64=4.0, lambda::Float64=0.0)
    MV = zeros(Float64, 3, state.np)
    galerkin_matvec!(MV, state, V; cutoff_factor, lambda)
    return relnorm(MV .- B, B)
end

function dense_metrics(M)
    symerr = norm(M .- transpose(M)) / (norm(M) + EPSN)
    ev = eigvals(Symmetric(M))
    return (; symerr, mineig=minimum(ev), maxeig=maximum(ev),
            cond=maximum(ev) / minimum(ev), eig=ev)
end

function select_region(state; region::AbstractString="tipvortex", np_target::Int=1000,
                       R::Float64=DEFAULT_R, axial_dim::Int=1)
    X, gamma = state.X, state.gamma
    gnorm = vec(sqrt.(sum(abs2, gamma; dims=1)))
    axial = vec(X[axial_dim, :])
    seed = if region == "oldwake"
        target = quantile(axial, 0.90)
        argmin(abs.(axial .- target))
    else
        candidates = findall(x -> 0.0 <= x <= R, axial)
        isempty(candidates) && (candidates = collect(1:state.np))
        candidates[argmax(gnorm[candidates])]
    end
    d = vec(sqrt.(sum(abs2, X .- X[:, seed]; dims=1)))
    ord = sortperm(d)
    interior = sort(ord[1:min(np_target, length(ord))])
    r_sel = maximum(d[interior])
    return seed, interior, r_sel, d
end

function buffered_indices(state, interior::Vector{Int}, d_seed, r_sel;
                          cutoff_factor::Float64=4.0, tau::Float64=1e-8)
    sigloc = state.sigma[interior]
    sigma_ring = cutoff_factor * sqrt(2.0) * maximum(sigloc)
    all_ring = findall(x -> x <= r_sel + sigma_ring, d_seed)
    if length(setdiff(all_ring, interior)) <= length(interior)
        return sort(all_ring), "sigma_ring", sigma_ring
    end
    mind = minimum(gram_diag.(sigloc))
    seff = sqrt(2.0) * maximum(sigloc)
    arg = tau * mind * seff^3 / ZETA0
    r_eff = arg > 0 && arg < 1 ? sqrt(-2 * seff^2 * log(arg)) : sigma_ring
    eff_ring = findall(x -> x <= r_sel + r_eff, d_seed)
    return sort(eff_ring), "effective_tau", r_eff
end

function write_csv(path, rows::Vector{<:NamedTuple})
    mkpath(dirname(path))
    isempty(rows) && return path
    names = propertynames(rows[1])
    open(path, "w") do io
        println(io, join(names, ","))
        for row in rows
            println(io, join((getproperty(row, n) for n in names), ","))
        end
    end
    return path
end

function save_reference(path::AbstractString; kwargs...)
    mkpath(dirname(path))
    open(path, "w") do io
        Serialization.serialize(io, NamedTuple(kwargs))
    end
    return path
end

function load_reference(path::AbstractString)
    open(path, "r") do io
        return Serialization.deserialize(io)
    end
end

end # module
