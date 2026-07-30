#!/usr/bin/env julia
# Dense subsystem conditioning and buffered-interior solve for item 008.

include(joinpath(@__DIR__, "particle_overlap_galerkin.jl"))
const POG = ParticleOverlapGalerkin

import LinearAlgebra: Symmetric, cholesky, Diagonal, diag, eigvals, norm, svdvals
import Printf: @printf
import Statistics: quantile

function getenv_parse(T, key, default)
    parse(T, get(ENV, key, string(default)))
end

run_name = get(ENV, "RUN_NAME", "rotor_hover_pressure_comparison")
step = getenv_parse(Int, "STEP", 350)
region = get(ENV, "REGION", "tipvortex")
np_target = getenv_parse(Int, "NP_TARGET", 1000)
cutoff_factor = getenv_parse(Float64, "CUTOFF_FACTOR", 4.0)
lambda_default = getenv_parse(Float64, "LAMBDA", 0.0)
buffer_max_factor = getenv_parse(Float64, "BUFFER_MAX_FACTOR", 1.0)

path = POG.particle_path(run_name, step)
state = POG.POR.load_state(path)
@printf("Loaded %s: np=%d\n", path, state.np)

seed, I, r_sel, d_seed = POG.select_region(state; region, np_target)
IB, buffer_kind, buffer_width = POG.buffered_indices(state, I, d_seed, r_sel; cutoff_factor)
Bset = setdiff(IB, I)
if length(Bset) > round(Int, buffer_max_factor * length(I))
    outside = setdiff(sortperm(d_seed), I)
    keep = outside[1:round(Int, buffer_max_factor * length(I))]
    Bset = sort(keep)
    IB = sort(vcat(I, Bset))
    buffer_kind = "$(buffer_kind)_capped"
end
pos = Dict(id => j for (j, id) in enumerate(IB))
irows = [pos[i] for i in I]
brows = [pos[i] for i in Bset]
@printf("Region=%s seed=%d interior=%d buffer=%d kind=%s width=%.6g\n",
        region, seed, length(I), length(Bset), buffer_kind, buffer_width)

rhs = POG.build_rhs(state; cutoff_factor, include_sfs=false)
gdot0 = POG.isolated_rate(state; include_sfs=false)

Mfull = POG.dense_gram(state, IB; cutoff_factor=Inf, lambda=0.0)
base = POG.dense_metrics(Mfull)
Dfull = Diagonal([POG.gram_diag(state.sigma[i]) for i in IB])
Mtr = POG.dense_gram(state, IB; cutoff_factor=cutoff_factor, lambda=0.0)
trunc_error = norm(Mtr .- Mfull) / (norm(Mfull) + POG.EPSN)
tr = POG.dense_metrics(Mtr)
Mmed = POG.dense_gram(state, IB; cutoff_factor=Inf, lambda=0.0, median_sigma=true)
medm = POG.dense_metrics(Mmed)

colloc_cond = NaN
if length(IB) <= 3000
    C = POG.dense_collocation(state, IB)
    s = svdvals(C)
    colloc_cond = maximum(s) / minimum(s)
else
    @printf("Skipping dense collocation SVD at n=%d (reported as NaN).\n", length(IB))
end

Pinvhalf = Diagonal(1.0 ./ sqrt.(diag(Dfull)))
Pmat = Pinvhalf * Mfull * Pinvhalf
pev = eigvals(Symmetric(Pmat))
pq = [quantile(pev, q) for q in 0.0:0.1:1.0]

MII = Mfull[irows, irows]
BI = rhs.B[:, I]
rhs_eff = copy(BI)
if !isempty(Bset)
    MIB = Mfull[irows, brows]
    rhs_eff .-= permutedims(MIB * transpose(gdot0[:, Bset]))
end
diagI = [POG.gram_diag(state.sigma[i]) for i in I]
evI = eigvals(Symmetric(MII))
solve_lambda = lambda_default
if minimum(evI) <= 0 && solve_lambda == 0.0
    solve_lambda = -minimum(evI) / minimum(diagI) + 1e-12
end
Msolve = MII + solve_lambda * Diagonal(diagI)
chol = cholesky(Symmetric(Msolve))
x0 = gdot0[:, I]
xstar = permutedims(chol \ transpose(rhs_eff))
rg_sub = norm(permutedims(Msolve * transpose(xstar)) .- rhs_eff) / (norm(rhs_eff) + POG.EPSN)
corr_size = norm(xstar .- x0) / (norm(x0) + POG.EPSN)

lambda_rows = NamedTuple[]
xbase = xstar
for lam in (0.0, 1e-8, 1e-6, 1e-4, 1e-2)
    Mlam = MII + lam * Diagonal([POG.gram_diag(state.sigma[i]) for i in I])
    ev = eigvals(Symmetric(Mlam))
    xl = if minimum(ev) > 0
        permutedims(cholesky(Symmetric(Mlam)) \ transpose(rhs_eff))
    else
        fill(NaN, size(xbase))
    end
    push!(lambda_rows, (; run_name, step, region, np_target, n_interior=length(I),
        n_buffer=length(Bset), buffer_kind, lambda=lam, mineig=minimum(ev),
        cond=maximum(ev)/minimum(ev), solution_drift=norm(xl .- xbase)/(norm(xbase)+POG.EPSN)))
end

rows = NamedTuple[]
for cf in (Inf, 6.0, 4.0, 3.0, 2.0)
    M = POG.dense_gram(state, IB; cutoff_factor=cf, lambda=0.0)
    m = POG.dense_metrics(M)
    terr = isinf(cf) ? 0.0 : norm(M .- Mfull) / (norm(Mfull) + POG.EPSN)
    push!(rows, (; run_name, step, region, np_target, n_total=length(IB),
        n_interior=length(I), n_buffer=length(Bset), buffer_kind, buffer_width,
        buffer_max_factor,
        cutoff_factor=cf, lambda=lambda_default, sig_min=minimum(state.sigma[IB]),
        sig_med=POG.median(state.sigma[IB]), sig_max=maximum(state.sigma[IB]),
        symerr=m.symerr, mineig=m.mineig, cond=m.cond, colloc_cond,
        trunc_error=terr, median_sigma_cond=medm.cond,
        jacobi_cond=maximum(pev)/minimum(pev), jacobi_q0=pq[1], jacobi_q1=pq[2],
        jacobi_q5=pq[6], jacobi_q9=pq[10], jacobi_q10=pq[11],
        interior_mineig=minimum(evI), solve_lambda,
        buffered_rg=rg_sub, correction_size=corr_size))
end

outdir = joinpath("data", "overlap_subsystem_cond")
csv = joinpath(outdir, "$(run_name)_step$(step)_$(region)_np$(np_target).csv")
POG.write_csv(csv, rows)
POG.write_csv(joinpath(outdir, "$(run_name)_step$(step)_$(region)_np$(np_target)_lambda.csv"), lambda_rows)

ref_path = joinpath(outdir, "$(run_name)_step$(step)_$(region)_np$(np_target)_reference.bin")
POG.save_reference(ref_path; run_name, step, region, cutoff_factor, lambda=lambda_default,
    solve_lambda, interior=I, buffer=Bset, seed, xstar, rhs_eff, rg_sub, correction_size=corr_size)

@printf("SPD base: symerr=%.3e mineig=%.3e cond=%.3e\n", base.symerr, base.mineig, base.cond)
@printf("Truncated cf=%.2f: mineig=%.3e cond=%.3e trunc_error=%.3e\n",
        cutoff_factor, tr.mineig, tr.cond, trunc_error)
@printf("Jacobi-preconditioned cond=%.3e q=[%.3e %.3e %.3e]\n",
        maximum(pev)/minimum(pev), pq[1], pq[6], pq[11])
@printf("Buffered solve: interior_mineig=%.3e solve_lambda=%.3e r_G=%.3e correction_size=%.3e\n",
        minimum(evI), solve_lambda, rg_sub, corr_size)
@printf("Wrote %s\n", csv)
@printf("Saved reference %s\n", ref_path)
