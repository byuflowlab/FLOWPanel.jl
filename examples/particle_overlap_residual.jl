# =============================================================================
# Particle-overlap vorticity-evolution residual diagnostic (BRAINSTORM item 008)
#
# Reusable, side-effect-free helpers for the *collocation* mat-vec residual
# diagnostic described in `BRAINSTORM/008_particle_overlap_vorticity_evolution.md`.
#
# The diagnostic poses the filtered-vorticity strength-update system  M Γ̇ = b
# and measures, on a SAVED particle state, how far the isolated single-particle
# strength-rate estimates sit from the overlap-projected RHS. It does NOT change
# any particle update; it answers the go/no-go question of whether an overlap
# solve would materially alter the live FLOWVPM behavior.
#
# Trust of saved fields (see item md / plan): inputs are the saved primitives
# (X, Γ, σ) plus the saved velocity-gradient J (the SAME J the live `_euler`
# consumed, evaluated at exactly the saved (X,Γ,σ); with the stable-wake run's
# flags it is the clean particle-only induced gradient WITH the self/overlap
# term). Saved `velocity` enters the convective term; saved `vorticity` is used
# ONLY as an independent cross-check / kernel-normalization calibration — the
# basis vorticity ω = Σ_j Γ_j ζ_j is recomputed here.
#
# Conventions (verified against FLOWVPM_timeintegration.jl `_euler`, lines 66-76,
# and FLOWVPM_particlefield.jl J/Γ index layout):
#   * The stored 9-vector J reshapes column-major to M[a,b] = ∂u_a/∂x_b.
#   * Physical stretching (ω·∇)u = M·ω is the `_euler` "classic" branch formula.
#   * The live FLOWVPM update (default `transposed=true`) uses the "transposed"
#     branch formula. These two are kept strictly separate per the item.
#   * get_W1 = J[6]-J[8] = ∂u₃/∂x₂-∂u₂/∂x₃ confirms the layout.
# =============================================================================
module ParticleOverlapResidual

import ReadVTK
import FLOWVPM
import LinearAlgebra: norm

# Particle regularization kernel. The run's velocity kernel is FLOWVPM's default
# Gaussian-error-function kernel; ζ(r)=const1·exp(-r²/2). Calibration against the
# saved `vorticity` field confirms this choice and its normalization.
const ZETA  = FLOWVPM.zeta_gauserf      # ζ(r), scalar of the *scaled* radius r=|d|/σ
const ZETA0 = FLOWVPM.zeta_gauserf(0.0) # ζ(0) = const1 = 1/(2π)^{3/2}

"""
    zeta_sigma(d2, sigma)

3-D regularized kernel value ζ_σ(d⃗) = ζ(|d⃗|/σ)/σ³, given the *squared* distance
`d2 = |d⃗|²` and core size `sigma`. (Gaussian form; matches `ZETA`.)
"""
@inline function zeta_sigma(d2::Float64, sigma::Float64)
    s2 = sigma * sigma
    return ZETA0 * exp(-d2 / (2 * s2)) / (sigma * s2)   # const1·exp(-ρ²/2)/σ³
end

# ---------------------------------------------------------------------------
# Loading a saved particle VTP state (mirrors src/FLOWPanel_warmstart.jl:229)
# ---------------------------------------------------------------------------
"""
    load_state(vtp_path) -> NamedTuple

Read one particle VTP file. Returns `(; X, gamma, sigma, velocity, vorticity, J,
vol, C, SFS, np)` with `X, gamma, velocity, vorticity, SFS` as 3×np,
`sigma`, `vol`, `C` as np-vectors, `J` as 9×np (column-major
M[a,b]=∂u_a/∂x_b per particle). Optional fields default to zeros when absent.
"""
function load_state(vtp_path::AbstractString)
    isfile(vtp_path) || error("Particle VTP not found: $(vtp_path)")
    vtk = ReadVTK.VTKFile(vtp_path)
    np = vtk.n_points
    pd = ReadVTK.get_point_data(vtk)
    required = ("gamma", "sigma", "velocity", "vorticity", "velocity_gradient")
    miss = filter(f -> !(f in keys(pd)), required)
    isempty(miss) || error("VTP missing field(s): $(join(miss, ", "))")

    X = Matrix{Float64}(ReadVTK.get_points(vtk))                  # 3×np
    gamma     = Matrix{Float64}(reshape(ReadVTK.get_data(pd["gamma"]), 3, np))
    sigma     = Vector{Float64}(vec(ReadVTK.get_data(pd["sigma"])))
    velocity  = Matrix{Float64}(reshape(ReadVTK.get_data(pd["velocity"]), 3, np))
    vorticity = Matrix{Float64}(reshape(ReadVTK.get_data(pd["vorticity"]), 3, np))
    J         = Matrix{Float64}(reshape(ReadVTK.get_data(pd["velocity_gradient"]), 9, np))
    vol = "vol" in keys(pd) ? Vector{Float64}(vec(ReadVTK.get_data(pd["vol"]))) : zeros(Float64, np)
    C = "C" in keys(pd) ? Vector{Float64}(vec(ReadVTK.get_data(pd["C"]))) : zeros(Float64, np)
    SFS = "SFS" in keys(pd) ? Matrix{Float64}(reshape(ReadVTK.get_data(pd["SFS"]), 3, np)) : zeros(Float64, 3, np)
    return (; X, gamma, sigma, velocity, vorticity, J, vol, C, SFS, np)
end

# ---------------------------------------------------------------------------
# Per-particle stretching operators (act on a single particle's 9-vector J)
# ---------------------------------------------------------------------------
"""
    s_phys(J, w) -> 3-tuple

Physical stretching (w·∇)u = M·w with stored layout M[a,b]=∂u_a/∂x_b
(the `_euler` "classic" branch formula). `J` is the 9-vector, `w` a 3-vector.
"""
@inline function s_phys(J, w)
    return (J[1]*w[1] + J[4]*w[2] + J[7]*w[3],
            J[2]*w[1] + J[5]*w[2] + J[8]*w[3],
            J[3]*w[1] + J[6]*w[2] + J[9]*w[3])
end

"""
    s_fvpm(J, G; transposed=true) -> 3-tuple

The exact FLOWVPM `_euler` stretching branch. With `transposed=true` (FLOWVPM
default) this is the "Γ⋅∇'U" transposed branch; with `false` it falls back to
the classic branch (= `s_phys`).
"""
@inline function s_fvpm(J, G; transposed::Bool=true)
    if transposed
        return (J[1]*G[1] + J[2]*G[2] + J[3]*G[3],
                J[4]*G[1] + J[5]*G[2] + J[6]*G[3],
                J[7]*G[1] + J[8]*G[2] + J[9]*G[3])
    else
        return s_phys(J, G)
    end
end

"""
    curl_from_J(state) -> 3×np

Vorticity from the stored velocity gradient, `ω = (J₆−J₈, J₇−J₃, J₂−J₄)` per
particle (FLOWVPM `get_W` layout). In regularized VPM, `∇×u_regularized` equals
the basis vorticity `Σ_j Γ_j ζ_j`, so this is the calibration reference for the
kernel normalization (the saved `vorticity` point-field is empty in some runs).
"""
function curl_from_J(state)
    np = state.np
    W = Matrix{Float64}(undef, 3, np)
    @inbounds for i in 1:np
        J = view(state.J, :, i)
        W[1,i] = J[6] - J[8]
        W[2,i] = J[7] - J[3]
        W[3,i] = J[2] - J[4]
    end
    return W
end

"""
    isolated_rates(state; transposed=true) -> (Vphys, Vfvpm)

Per-particle isolated single-particle strength-rate estimates (3×np each):
`Vphys[:,i] = (Γ_i·∇)u = J_i Γ_i` (physical) and
`Vfvpm[:,i] = S_FLOWVPM(J_i, Γ_i)` (live FLOWVPM branch).
"""
function isolated_rates(state; transposed::Bool=true)
    np = state.np
    Vphys = Matrix{Float64}(undef, 3, np)
    Vfvpm = Matrix{Float64}(undef, 3, np)
    @inbounds for i in 1:np
        Ji = view(state.J, :, i)
        Gi = view(state.gamma, :, i)
        sp = s_phys(Ji, Gi)
        sf = s_fvpm(Ji, Gi; transposed)
        Vphys[1,i], Vphys[2,i], Vphys[3,i] = sp
        Vfvpm[1,i], Vfvpm[2,i], Vfvpm[3,i] = sf
    end
    return Vphys, Vfvpm
end

# ---------------------------------------------------------------------------
# Neighbor search: uniform hash-grid cell list (no extra dependency)
# ---------------------------------------------------------------------------
struct CellGrid
    cell::Float64                                  # cell size = cutoff radius
    inv::Float64
    buckets::Dict{NTuple{3,Int}, Vector{Int}}
end

@inline _cellidx(X, i, inv) = (floor(Int, X[1,i]*inv), floor(Int, X[2,i]*inv), floor(Int, X[3,i]*inv))

"""
    build_grid(X, np, cutoff) -> CellGrid

Hash-grid with cell size = `cutoff`, so neighbors within `cutoff` lie in the
27 cells around a target's cell.
"""
function build_grid(X, np, cutoff::Float64)
    inv = 1.0 / cutoff
    buckets = Dict{NTuple{3,Int}, Vector{Int}}()
    @inbounds for i in 1:np
        key = _cellidx(X, i, inv)
        push!(get!(() -> Int[], buckets, key), i)
    end
    return CellGrid(cutoff, inv, buckets)
end

"""
    foreach_neighbor(f, grid, X, k)

Call `f(i, d2)` for every particle `i` whose squared distance `d2` to target `k`
is ≤ cutoff² (includes `i==k`, d2=0).
"""
@inline function foreach_neighbor(f, grid::CellGrid, X, k::Int)
    c2 = grid.cell * grid.cell
    ck = _cellidx(X, k, grid.inv)
    xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
    @inbounds for dz in -1:1, dy in -1:1, dx in -1:1
        bucket = get(grid.buckets, (ck[1]+dx, ck[2]+dy, ck[3]+dz), nothing)
        bucket === nothing && continue
        for i in bucket
            d1 = X[1,i]-xk1; d2v = X[2,i]-xk2; d3 = X[3,i]-xk3
            dd = d1*d1 + d2v*d2v + d3*d3
            dd <= c2 && f(i, dd)
        end
    end
end

# ---------------------------------------------------------------------------
# Overlap mat-vec and convective-basis term
# ---------------------------------------------------------------------------
"""
    overlap_matvec(grid, state, V) -> 3×np

`(M V)[:,k] = Σ_i ζ_{σ_i}(X_k - X_i) V[:,i]` over neighbors of each target k.
With `V = Γ` this is the basis vorticity ω(X_k).
"""
function overlap_matvec(grid::CellGrid, state, V::AbstractMatrix)
    X, sigma, np = state.X, state.sigma, state.np
    out = zeros(Float64, 3, np)
    @inbounds for k in 1:np
        a1 = 0.0; a2 = 0.0; a3 = 0.0
        foreach_neighbor(grid, X, k) do i, d2
            z = zeta_sigma(d2, sigma[i])
            a1 += z*V[1,i]; a2 += z*V[2,i]; a3 += z*V[3,i]
        end
        out[1,k] = a1; out[2,k] = a2; out[3,k] = a3
    end
    return out
end

# NOTE: the per-row sum (conditioning proxy) and the convective-basis term are
# computed inline in the fused `analyze` pass below — the single source of truth
# for that math. `overlap_matvec` (above) + the brute-force reference (below) are
# retained only for the validation gates.

# ---------------------------------------------------------------------------
# Brute-force O(N²) reference (validation only; no truncation)
# ---------------------------------------------------------------------------
function overlap_matvec_bruteforce(state, V::AbstractMatrix)
    X, sigma, np = state.X, state.sigma, state.np
    out = zeros(Float64, 3, np)
    @inbounds for k in 1:np
        for i in 1:np
            dx = X[1,i]-X[1,k]; dy = X[2,i]-X[2,k]; dz = X[3,i]-X[3,k]
            z = zeta_sigma(dx*dx+dy*dy+dz*dz, sigma[i])
            out[1,k] += z*V[1,i]; out[2,k] += z*V[2,i]; out[3,k] += z*V[3,i]
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Flat CSR cell list (counting-sort; no Dict/closure in the hot loop)
# ---------------------------------------------------------------------------
struct CSRGrid
    inv::Float64                              # 1/cell size
    cid::Vector{Int}                          # per-particle contiguous cell id
    key::Vector{NTuple{3,Int}}               # per-particle (cx,cy,cz)
    cellid::Dict{NTuple{3,Int},Int}          # (cx,cy,cz) -> cell id
    start::Vector{Int}                        # CSR offsets (length ncell+1)
    perm::Vector{Int}                         # particle ids grouped by cell
end

function build_csr(X, np, cutoff::Float64)
    inv = 1.0 / cutoff
    key = Vector{NTuple{3,Int}}(undef, np)
    cellid = Dict{NTuple{3,Int},Int}()
    cid = Vector{Int}(undef, np)
    @inbounds for i in 1:np
        k = (floor(Int, X[1,i]*inv), floor(Int, X[2,i]*inv), floor(Int, X[3,i]*inv))
        key[i] = k
        cid[i] = get!(cellid, k, length(cellid) + 1)
    end
    ncell = length(cellid)
    start = zeros(Int, ncell + 1)
    @inbounds for i in 1:np; start[cid[i]+1] += 1; end
    @inbounds for c in 1:ncell; start[c+1] += start[c]; end
    perm = Vector{Int}(undef, np)
    fill = copy(start)
    @inbounds for i in 1:np
        c = cid[i]; fill[c] += 1; perm[fill[c]] = i
    end
    return CSRGrid(inv, cid, key, cellid, start, perm)
end

# ---------------------------------------------------------------------------
# Residual assembly
# ---------------------------------------------------------------------------
relnorm(R, B; eps=1e-30) = norm(R) / (norm(B) + eps)

"""
    analyze(state; cutoff_factor=4.0, transposed=true) -> NamedTuple

Full per-state diagnostic (single fused neighbor pass). Returns the four relative
field-space residuals (reduced/sampled × physical/fvpm), the basis-ω-vs-curl-of-J
calibration error, conditioning/neighbor proxies, and timing. `cutoff_factor`
sets the truncation radius `r_cut = cutoff_factor · max σ`.
"""
function analyze(state; cutoff_factor::Float64=4.0, transposed::Bool=true)
    np = state.np
    np == 0 && return nothing
    X, sigma, gamma, U, Jm = state.X, state.sigma, state.gamma, state.velocity, state.J
    sig_max = maximum(sigma); sig_min = minimum(sigma)
    cutoff = cutoff_factor * sig_max
    c2 = cutoff * cutoff

    t0 = time()
    g = build_csr(X, np, cutoff)

    # isolated single-particle rate estimates (O(np))
    Vphys, Vfvpm = isolated_rates(state; transposed)

    # ---- single fused neighbor pass: ω=MΓ, M·Vphys, M·Vfvpm, rowsum, count, conv ----
    omega  = zeros(Float64, 3, np)
    Mvp    = zeros(Float64, 3, np)
    Mvf    = zeros(Float64, 3, np)
    conv   = zeros(Float64, 3, np)
    rowsum = zeros(Float64, np)
    nbr    = zeros(Int, np)
    @inbounds for k in 1:np
        kx, ky, kz = g.key[k]
        xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
        uk1, uk2, uk3 = U[1,k], U[2,k], U[3,k]
        o1=0.0; o2=0.0; o3=0.0; p1=0.0; p2=0.0; p3=0.0; f1=0.0; f2=0.0; f3=0.0
        cc1=0.0; cc2=0.0; cc3=0.0; rs=0.0; cnt=0
        for dz in -1:1, dy in -1:1, dx in -1:1
            c = get(g.cellid, (kx+dx, ky+dy, kz+dz), 0)
            c == 0 && continue
            @inbounds for p in (g.start[c]+1):g.start[c+1]
                i = g.perm[p]
                rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
                d2 = rx*rx + ry*ry + rz*rz
                d2 > c2 && continue
                si = sigma[i]; si2 = si*si
                z = ZETA0 * exp(-d2/(2*si2)) / (si*si2)
                o1 += z*gamma[1,i]; o2 += z*gamma[2,i]; o3 += z*gamma[3,i]
                p1 += z*Vphys[1,i]; p2 += z*Vphys[2,i]; p3 += z*Vphys[3,i]
                f1 += z*Vfvpm[1,i]; f2 += z*Vfvpm[2,i]; f3 += z*Vfvpm[3,i]
                rs += z; cnt += 1
                if i != k
                    # conv_k += Γ_i (u_k-u_i)·∇ζ_i ;  ∇ζ_i(X_k) = -z (X_k-X_i)/σ_i²
                    # note d⃗ = X_k - X_i = -(rx,ry,rz)
                    w = -z * ((uk1-U[1,i])*(-rx) + (uk2-U[2,i])*(-ry) + (uk3-U[3,i])*(-rz)) / si2
                    cc1 += gamma[1,i]*w; cc2 += gamma[2,i]*w; cc3 += gamma[3,i]*w
                end
            end
        end
        omega[1,k]=o1; omega[2,k]=o2; omega[3,k]=o3
        Mvp[1,k]=p1; Mvp[2,k]=p2; Mvp[3,k]=p3
        Mvf[1,k]=f1; Mvf[2,k]=f2; Mvf[3,k]=f3
        conv[1,k]=cc1; conv[2,k]=cc2; conv[3,k]=cc3
        rowsum[k]=rs; nbr[k]=cnt
    end

    # Calibration: in regularized VPM, ∇×u = Σ Γ_j ζ_j, so compare basis ω to
    # curl-of-J (always available). The saved `vorticity` point-field is empty in
    # some runs; report its norm so that case is visible rather than silent.
    omega_curl = curl_from_J(state)
    basis_curl_relerr  = relnorm(omega .- omega_curl, omega_curl)
    saved_vort_norm    = norm(state.vorticity)
    basis_savedvort_relerr = saved_vort_norm > 0 ? relnorm(omega .- state.vorticity, state.vorticity) : NaN

    # reduced RHS: J_k ω(X_k)
    b_red_phys = similar(omega); b_red_fvpm = similar(omega)
    @inbounds for k in 1:np
        Jk = view(Jm, :, k); wk = view(omega, :, k)
        sp = s_phys(Jk, wk); sf = s_fvpm(Jk, wk; transposed)
        b_red_phys[1,k], b_red_phys[2,k], b_red_phys[3,k] = sp
        b_red_fvpm[1,k], b_red_fvpm[2,k], b_red_fvpm[3,k] = sf
    end

    # sampled RHS: reduced minus convective-basis correction
    b_smp_phys = b_red_phys .- conv
    b_smp_fvpm = b_red_fvpm .- conv

    r_red_phys = relnorm(Mvp .- b_red_phys, b_red_phys)
    r_red_fvpm = relnorm(Mvf .- b_red_fvpm, b_red_fvpm)
    r_smp_phys = relnorm(Mvp .- b_smp_phys, b_smp_phys)
    r_smp_fvpm = relnorm(Mvf .- b_smp_fvpm, b_smp_fvpm)

    # conditioning proxy: diagonal dominance per row
    diag = [zeta_sigma(0.0, sigma[k]) for k in 1:np]
    domin = diag ./ max.(rowsum .- diag, eps())
    elapsed = time() - t0

    return (; np,
            cutoff, sig_min, sig_max,
            basis_curl_relerr, basis_savedvort_relerr, saved_vort_norm,
            conv_rel = relnorm(conv, b_red_phys),        # size of convective correction
            r_red_phys, r_red_fvpm, r_smp_phys, r_smp_fvpm,
            domin_min=minimum(domin), domin_med=_median(domin),
            nbr_mean=sum(nbr)/np, nbr_max=maximum(nbr),
            elapsed)
end

function _median(v)
    s = sort(v); n = length(s)
    return iseven(n) ? 0.5*(s[n÷2] + s[n÷2+1]) : s[(n+1)÷2]
end

end # module
