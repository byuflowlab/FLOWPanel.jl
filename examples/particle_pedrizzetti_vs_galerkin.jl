#!/usr/bin/env julia
# Pedrizzetti replay/rate comparison against the Galerkin residual.

include(joinpath(@__DIR__, "particle_overlap_galerkin.jl"))
const POG = ParticleOverlapGalerkin

import Printf: @printf
import LinearAlgebra: dot, norm

parse_grid(s) = parse.(Float64, split(s, ","))
parse_steps(s) = occursin(":", s) ? collect(parse(Int, split(s, ":")[1]):parse(Int, split(s, ":")[2])) : parse.(Int, split(s, ","))
getfloat(k, d) = parse(Float64, get(ENV, k, string(d)))

@inline function omega_from_J(J)
    return (J[6]-J[8], J[7]-J[3], J[2]-J[4])
end

function euler_stretch_sfs(state, dt)
    Gp = copy(state.gamma)
    @inbounds for i in 1:state.np
        J = view(state.J, :, i)
        G1 = Gp[1,i]; G2 = Gp[2,i]; G3 = Gp[3,i]
        G1 += dt * (J[1]*G1 + J[2]*G2 + J[3]*G3)
        G2 += dt * (J[4]*G1 + J[5]*G2 + J[6]*G3)
        G3 += dt * (J[7]*G1 + J[8]*G2 + J[9]*G3)
        fac = dt * state.C[i] * state.sigma[i]^3 / POG.ZETA0
        Gp[1,i] = G1 - fac * state.SFS[1,i]
        Gp[2,i] = G2 - fac * state.SFS[2,i]
        Gp[3,i] = G3 - fac * state.SFS[3,i]
    end
    return Gp
end

function relax_gamma!(G, state, rlxf; corrected::Bool)
    @inbounds for i in 1:state.np
        J = view(state.J, :, i)
        w = omega_from_J(J)
        nw = sqrt(w[1]^2 + w[2]^2 + w[3]^2)
        iszero(nw) && continue
        g0 = (G[1,i], G[2,i], G[3,i])
        ng = sqrt(g0[1]^2 + g0[2]^2 + g0[3]^2)
        iszero(ng) && continue
        cosang = (g0[1]*w[1] + g0[2]*w[2] + g0[3]*w[3]) / (ng*nw)
        G[1,i] = (1-rlxf)*g0[1] + rlxf*ng*w[1]/nw
        G[2,i] = (1-rlxf)*g0[2] + rlxf*ng*w[2]/nw
        G[3,i] = (1-rlxf)*g0[3] + rlxf*ng*w[3]/nw
        if corrected
            b2 = 1 - 2*(1-rlxf)*rlxf*(1 - cosang)
            s = sqrt(b2)
            G[1,i] /= s; G[2,i] /= s; G[3,i] /= s
        end
    end
    return G
end

function effective_rate(state, dt, rlxf; corrected::Bool)
    Gp = euler_stretch_sfs(state, dt)
    rlxf == 0 && return (Gp .- state.gamma) ./ dt
    relax_gamma!(Gp, state, rlxf; corrected)
    return (Gp .- state.gamma) ./ dt
end

function optimal_scaled_residual(state, base, cand, B; cutoff_factor, lambda)
    Mbase = zeros(Float64, 3, state.np)
    Mdel = zeros(Float64, 3, state.np)
    delta = cand .- base
    POG.galerkin_matvec!(Mbase, state, base; cutoff_factor, lambda)
    POG.galerkin_matvec!(Mdel, state, delta; cutoff_factor, lambda)
    rem = B .- Mbase
    denom = sum(Mdel .* Mdel)
    s = denom > 0 ? sum(Mdel .* rem) / denom : 0.0
    r = norm(Mbase .+ s .* Mdel .- B) / (norm(B) + POG.EPSN)
    return s, r
end

run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
steps = parse_steps(get(ENV, "STEPS", "350"))
rlxfs = parse_grid(get(ENV, "RLXF_GRID", "0,0.06,0.3,0.6,1.0"))
cutoff_factor = getfloat("CUTOFF_FACTOR", 4.0)
lambda = getfloat("LAMBDA", 0.0)
dt = getfloat("DT", 1/3600)

rows = NamedTuple[]
for step in steps
    state = POG.POR.load_state(POG.particle_path(run_name, step))
    rhs = POG.build_rhs(state; cutoff_factor, include_sfs=true)
    base = POG.isolated_rate(state; include_sfs=true)
    rg_base = POG.residual_galerkin(state, base, rhs.B; cutoff_factor, lambda)
    @printf("Step %d base SFS-on r_G=%.6e\n", step, rg_base)
    for variant in ("vanilla", "corrected")
        corrected = variant == "corrected"
        for rlxf in rlxfs
            rate = effective_rate(state, dt, rlxf; corrected)
            rg = POG.residual_galerkin(state, rate, rhs.B; cutoff_factor, lambda)
            delta = rate .- base
            sopt, rg_opt = optimal_scaled_residual(state, base, rate, rhs.B; cutoff_factor, lambda)
            corr = norm(delta) / (norm(base) + POG.EPSN)
            push!(rows, (; run_name, step, variant, rlxf, cutoff_factor, lambda,
                rg_base, rg, sopt, rg_opt, correction_size=corr))
            @printf("  %-9s rlxf=%5.2f r_G=%.6e s*=%.3e r_G(s*)=%.6e corr=%.3e\n",
                    variant, rlxf, rg, sopt, rg_opt, corr)
        end
    end
end

out = joinpath("data", "pedrizzetti_vs_galerkin", "$(run_name)_summary.csv")
POG.write_csv(out, rows)
@printf("Wrote %s\n", out)

