#!/usr/bin/env julia
# Wake-age-resolved Galerkin residual attribution for item 008.
#
# Attributes the full-scale Galerkin residual R = M_G Γ̇⁽⁰⁾ − b_G of the isolated
# single-particle update across wake-age bands (axial-position proxy) and σ bands
# (merge-history proxy). Row-restricted residual with the FULL-cloud operator:
# this is an attribution of the global residual, not a subsystem solve.
#
# Hypothesis under test (see the 2026-07-02 sections of the item-008 brainstorm
# doc): if the residual is concentrated in the old/merged/distorted wake and
# small in the young, freshly-shed wake, it is distortion-driven (classical
# quadrature-consistency broken), favoring remeshing/filament remedies and
# young-wake-restricted solves; if flat across age, it is intrinsic to the
# update rule at this overlap density.

include(joinpath(@__DIR__, "particle_overlap_galerkin.jl"))
const POG = ParticleOverlapGalerkin

import Printf: @printf
import Statistics: mean, median, quantile
import LinearAlgebra: norm
import Serialization

function parse_steps(s)
    if occursin(":", s)
        a, b = parse.(Int, split(s, ":"))
        return collect(a:b)
    end
    return parse.(Int, split(s, ","))
end

getbool(k, d) = parse(Bool, get(ENV, k, string(d)))
getfloat(k, d) = parse(Float64, get(ENV, k, string(d)))
getint(k, d) = parse(Int, get(ENV, k, string(d)))

run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
steps = parse_steps(get(ENV, "STEPS", "350"))
cutoff_factor = getfloat("CUTOFF_FACTOR", 4.0)
include_sfs = getbool("INCLUDE_SFS", false)
nbands = getint("NBANDS", 10)
axial_dim = getint("AXIAL_DIM", 1)
per_particle = getbool("PER_PARTICLE", false)

outdir = joinpath("data", "overlap_age_diag")

# Per-particle nearest-neighbor distance and near-duplicate flag (one CSR sweep).
function nn_sweep(state; cutoff_factor::Float64=4.0)
    np = state.np
    X, sigma = state.X, state.sigma
    cutoff = cutoff_factor * maximum(sigma)
    c2 = cutoff * cutoff
    g = POG.POR.build_csr(X, np, cutoff)
    hnn = fill(Inf, np)
    neardup = falses(np)
    Threads.@threads for k in 1:np
        kx, ky, kz = g.key[k]
        xk1, xk2, xk3 = X[1,k], X[2,k], X[3,k]
        sk = sigma[k]
        best = Inf
        dup = false
        for dz in -1:1, dy in -1:1, dx in -1:1
            cid = get(g.cellid, (kx+dx, ky+dy, kz+dz), 0)
            cid == 0 && continue
            @inbounds for p in (g.start[cid]+1):g.start[cid+1]
                i = g.perm[p]
                i == k && continue
                rx = X[1,i]-xk1; ry = X[2,i]-xk2; rz = X[3,i]-xk3
                d2 = rx*rx + ry*ry + rz*rz
                d2 > c2 && continue
                d2 < best && (best = d2)
                thr = 0.25 * min(sk, sigma[i])
                d2 < thr*thr && (dup = true)
            end
        end
        hnn[k] = sqrt(best)
        neardup[k] = dup
    end
    return hnn, neardup
end

# Equal-count quantile bands of `key`; returns vector of index vectors (low→high).
function quantile_bands(key::Vector{Float64}, nbands::Int)
    ord = sortperm(key)
    n = length(key)
    return [ord[(1 + div((b-1)*n, nbands)):div(b*n, nbands)] for b in 1:nbands]
end

function band_rows(bandkey_name, bands, key, state, R, B, gnorm, hnn, neardup;
                   run_name, step)
    # Per-particle row norms; robust per-row relative residual guarded by a
    # floor at 1e-6 × median row-RHS scale (heterogeneous-σ states can have
    # 1/σ³ outliers that dominate all energy-weighted norms).
    normR = vec(sqrt.(sum(abs2, R; dims=1)))
    normB = vec(sqrt.(sum(abs2, B; dims=1)))
    floorB = 1e-6 * median(normB)
    rrel = normR ./ max.(normB, floorB)
    rows = NamedTuple[]
    for (b, S) in enumerate(bands)
        RS = @view R[:, S]
        BS = @view B[:, S]
        overlap = state.sigma[S] ./ hnn[S]
        push!(rows, (; run_name, step, bandkey=bandkey_name, band=b,
            count=length(S),
            key_lo=minimum(key[S]), key_hi=maximum(key[S]),
            r_G=norm(RS) / (norm(BS) + POG.EPSN),
            rrel_med=median(rrel[S]), rrel_q90=quantile(rrel[S], 0.9),
            share_B2=sum(abs2, BS), share_R2=sum(abs2, RS),   # normalized below
            sigma_mean=mean(state.sigma[S]), sigma_med=median(state.sigma[S]),
            gamma_med=median(gnorm[S]),
            h_med=median(hnn[S]), overlap_med=median(overlap),
            neardup_frac=count(neardup[S]) / length(S)))
    end
    # normalize energy shares
    tB = sum(r.share_B2 for r in rows); tR = sum(r.share_R2 for r in rows)
    return [merge(r, (; share_B2=r.share_B2/tB, share_R2=r.share_R2/tR)) for r in rows]
end

function print_band_table(rows)
    @printf("%-6s %-4s %7s  [%9s, %9s]  %8s %8s %8s  %7s %7s  %8s %8s %7s\n",
            "key", "band", "count", "key_lo", "key_hi", "r_G", "rrelmed", "rrelq90",
            "shr_B2", "shr_R2", "sig_med", "ovlp_med", "dupfrac")
    for r in rows
        @printf("%-6s %-4d %7d  [%9.4f, %9.4f]  %8.4f %8.4f %8.4f  %7.4f %7.4f  %8.4f %8.2f %7.3f\n",
                r.bandkey, r.band, r.count, r.key_lo, r.key_hi, r.r_G,
                r.rrel_med, r.rrel_q90,
                r.share_B2, r.share_R2, r.sigma_med, r.overlap_med, r.neardup_frac)
    end
end

# EXTRA_J: path to a Serialization dump with fields (np, J_body, J_wrow)
# (written by particle_body_strain_audit.jl). When set, the omitted body and
# panel-wake-row gradients are ADDED to the saved particle-only J, so both
# sides of the overlap equation use the total strain. Output CSVs get a
# "_jtotal" suffix.
extra_j_path = get(ENV, "EXTRA_J", "")
suffix = isempty(extra_j_path) ? "" : "_jtotal"

for step in steps
    state = POG.POR.load_state(POG.particle_path(run_name, step))
    @printf("Loaded step %d np=%d include_sfs=%s cutoff_factor=%.1f\n",
            step, state.np, include_sfs, cutoff_factor)
    if !isempty(extra_j_path)
        extra = open(Serialization.deserialize, extra_j_path, "r")
        @assert extra.np == state.np "EXTRA_J np=$(extra.np) != state np=$(state.np)"
        state.J .+= extra.J_body .+ extra.J_wrow
        println("Added body + wake-row gradients from $(extra_j_path)")
    end

    rhs = POG.build_rhs(state; cutoff_factor, include_sfs)
    B = rhs.B
    gdot0 = POG.isolated_rate(state; include_sfs)
    MV = zeros(Float64, 3, state.np)
    POG.galerkin_matvec!(MV, state, gdot0; cutoff_factor)
    R = MV .- B

    # Whole-cloud check: must reproduce the recorded full-scale r_G(Γ̇⁰)
    # (0.4504575 for step 350, SFS off, cutoff 4, λ=0).
    rg_global = norm(R) / (norm(B) + POG.EPSN)
    @printf("Global r_G(Γ̇⁰) = %.7f\n", rg_global)

    hnn, neardup = nn_sweep(state; cutoff_factor)
    gnorm = vec(sqrt.(sum(abs2, state.gamma; dims=1)))
    axial = Vector{Float64}(vec(state.X[axial_dim, :]))
    sig = Vector{Float64}(state.sigma)

    rows_ax = band_rows("axial", quantile_bands(axial, nbands), axial,
                        state, R, B, gnorm, hnn, neardup; run_name, step)
    rows_sg = band_rows("sigma", quantile_bands(sig, nbands), sig,
                        state, R, B, gnorm, hnn, neardup; run_name, step)

    # Partition sanity: band energies must sum to the totals.
    for (rows, key) in ((rows_ax, axial), (rows_sg, sig))
        n_sum = sum(r.count for r in rows)
        e_B = abs(sum(r.share_B2 for r in rows) - 1.0)
        e_R = abs(sum(r.share_R2 for r in rows) - 1.0)
        @assert n_sum == state.np "band counts sum to $(n_sum) != np=$(state.np)"
        @assert e_B < 1e-12 && e_R < 1e-12 "band energy shares do not sum to 1"
    end

    println("\nAxial bands (age proxy; larger axial = older wake; NOTE: proxy only — " *
            "rollup makes age non-monotone near the rotor plane):")
    print_band_table(rows_ax)
    println("\nSigma bands (merge-history proxy):")
    print_band_table(rows_sg)

    POG.write_csv(joinpath(outdir, "$(run_name)_step$(step)_axial$(suffix).csv"), rows_ax)
    POG.write_csv(joinpath(outdir, "$(run_name)_step$(step)_sigma$(suffix).csv"), rows_sg)
    @printf("\nWrote %s\n", joinpath(outdir, "$(run_name)_step$(step)_{axial,sigma}$(suffix).csv"))

    if per_particle
        pp = [(; k, x1=state.X[1,k], x2=state.X[2,k], x3=state.X[3,k],
               axial=axial[k], sigma=sig[k], hnn=hnn[k], neardup=Int(neardup[k]),
               normR=norm(@view R[:,k]), normB=norm(@view B[:,k]))
              for k in 1:state.np]
        POG.write_csv(joinpath(outdir, "$(run_name)_step$(step)_particles.csv"), pp)
        @printf("Wrote %s\n", joinpath(outdir, "$(run_name)_step$(step)_particles.csv"))
    end
end
