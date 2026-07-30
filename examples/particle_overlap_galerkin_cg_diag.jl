#!/usr/bin/env julia
# Matrix-free Galerkin PCG diagnostic for item 008.

include(joinpath(@__DIR__, "particle_overlap_galerkin.jl"))
const POG = ParticleOverlapGalerkin

import Printf: @printf
import Statistics: mean
import LinearAlgebra: norm

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
lambda = getfloat("LAMBDA", 0.0)
include_sfs = getbool("INCLUDE_SFS", false)
maxiter = getint("MAXITER", 150)
tol = getfloat("TOL", 1e-6)

rows = NamedTuple[]
for step in steps
    state = POG.POR.load_state(POG.particle_path(run_name, step))
    @printf("Loaded step %d np=%d include_sfs=%s lambda=%.3e\n", step, state.np, include_sfs, lambda)
    rhs = POG.build_rhs(state; cutoff_factor, include_sfs)
    gdot0 = POG.isolated_rate(state; include_sfs)
    rg0 = POG.residual_galerkin(state, gdot0, rhs.B; cutoff_factor, lambda)
    @printf("Initial r_G=%.6e\n", rg0)
    t0 = time()
    pcg = POG.jacobi_pcg(state, rhs.B; cutoff_factor, lambda, tol, maxiter, x0=gdot0)
    elapsed = time() - t0
    rg = POG.residual_galerkin(state, pcg.x, rhs.B; cutoff_factor, lambda)
    corr = norm(pcg.x .- gdot0) / (norm(gdot0) + POG.EPSN)
    @printf("PCG iter=%d converged=%s monotone=%s breakdown=%s final_hist=%.6e r_G=%.6e corr=%.6e elapsed=%.1fs\n",
            pcg.iterations, pcg.converged, pcg.monotone, pcg.breakdown,
            pcg.history[end], rg, corr, elapsed)
    hist_path = joinpath("data", "overlap_galerkin_cg",
        "$(run_name)_step$(step)_sfs$(include_sfs)_lambda$(lambda)_history.csv")
    mkpath(dirname(hist_path))
    open(hist_path, "w") do io
        println(io, "iter,relres")
        for (i, r) in enumerate(pcg.history)
            println(io, "$(i-1),$r")
        end
    end
    push!(rows, (; run_name, step, np=state.np, cutoff_factor, lambda, include_sfs,
        maxiter, tol, rg_initial=rg0, rg_final=rg, iterations=pcg.iterations,
        converged=pcg.converged, monotone=pcg.monotone, breakdown=pcg.breakdown,
        correction_size=corr, elapsed))
end

out = joinpath("data", "overlap_galerkin_cg", "$(run_name)_summary.csv")
POG.write_csv(out, rows)
@printf("Wrote %s\n", out)

