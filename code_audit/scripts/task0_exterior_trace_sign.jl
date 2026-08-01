#=
Task 0 — Exterior-trace sign check (code_audit/plan.md)

Verifies that the `phi -= mu` conversion at
src/FLOWPanel_simulate_monitors.jl:391 (inside `_pressure_bernoulli_phi_dot!`)
produces the EXTERIOR trace of the scalar potential.

Two checks:

(A) Kernel-level (single panel, exact): for several individual body panels,
    evaluate the panel's own scalar potential via `pnl.induced` at the control
    point (self-limit path) and at cp ± eps*n_hat. Establishes which side of
    the panel the self-limit (+mu/2 by _self_limit) corresponds to, and that
    exterior = interior - mu under FLOWPanel's conventions.

(B) Full-body (monitor recipe): steady-solve the pitching_wing.jl wing at
    fixed alpha, probe total scalar potential (all panels + wake + uinf.x)
    at cp (self-limit) and at cp ± eps*ell*n_hat via the same
    `influence!((probes,), sources, backend)` call the monitor uses, and
    compare the recipe phi_surface = phi(cp) + uinf.x - mu against the
    off-surface exterior probe, the interior probe, and the raw interior value.

Run:  julia --project code_audit/scripts/task0_exterior_trace_sign.jl
=#

import FLOWPanel as pnl
import FLOWPanel.FastMultipole as FastMultipole
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: norm, dot
import Printf: @printf
import Statistics: mean, median

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(REPO, "examples", "pitching_wing.jl"))

# ---------------------------------------------------------------- build/solve
c = 1.0 * FT_TO_M
b = 4.0 * c
alpha_deg = 5.0
U = 30.0
rho = 1.225

body = build_pitching_wing_body(c, b;
    n_span=7, n_airfoil=31, n_endcap=5, endcap=:round, semiinfinite_wake=true)
Uinf = _uinf_from_alpha(U, alpha_deg)
set_wake_Das!(body, Uinf)   # defined in examples/pitching_wing.jl

backend = pnl.FastMultipoleBackend(expansion_order=10, multipole_acceptance=0.4,
                                   leaf_size=100)  # same default as PressureBernoulli
pressure = pnl.PressureBernoulli(rho; file=false)

pnl.steady!(body, pnl.ReferenceFrame(body), Uinf;
    body_solvers=pnl.Backslash(body), backend, monitors=(pressure,),
    path=nothing, name="task0", verbose=false)

n = body.ncells
i_mu = pnl.get_Gammai(body)
@assert pnl.has_grad_mu(body)
mu = [body.strength[p, i_mu] for p in 1:n]
println("wing: ncells = $n, alpha = $(alpha_deg) deg, |mu| range = ",
        extrema(abs.(mu)))

# normals: confirm outward (dot with cp - body center)
center = vec(mean(body.controlpoints, dims=2))
outward = [dot(SVector{3}(body.normals[:, p]), SVector{3}(body.controlpoints[:, p]) - center)
           for p in 1:n]
println("normals outward check: fraction with n.(cp-center) > 0 = ",
        count(>(0), outward) / n, "  (thin wing: near-TE/root panels may dip)")

# characteristic length per panel
ell = [pnl.characteristiclength_sqrtarea(body.nodes, view(body.cells, :, p)) for p in 1:n]
eps_rel = 1e-3

# which panels shed (skip for kernel-level check; their _induced_wake adds wake phi)
shed_panels = Set{Int}()
for shed in body.shedding
    shed isa AbstractMatrix || continue
    for j in axes(shed, 2), r in (1, 4)
        pa = shed[r, j]
        pa > 0 && push!(shed_panels, pa)
    end
end

# ----------------------------------------------------- (A) kernel-level check
println("\n=== (A) kernel-level single-panel check (pnl.induced) ===")
println("self-limit doc claims: phi(cp) = interior = +mu/2 (+ source PV);")
println("exterior should be interior - mu.  Testing d± = phi(cp ± eps*n) - phi(cp):")
ds = FastMultipole.DerivativesSwitch(true, false, false)
test_panels = [p for p in round.(Int, range(1, n; length=12)) if !(p in shed_panels)]
ratios_plus = Float64[]; ratios_minus = Float64[]
for p in test_panels
    cp = SVector{3}(body.controlpoints[:, p])
    nh = SVector{3}(body.normals[:, p])
    h = eps_rel * ell[p]
    phi0, _, _ = pnl.induced(cp, body, p, ds)
    phip, _, _ = pnl.induced(cp + h * nh, body, p, ds)
    phim, _, _ = pnl.induced(cp - h * nh, body, p, ds)
    push!(ratios_plus, (phip - phi0) / mu[p])
    push!(ratios_minus, (phim - phi0) / mu[p])
    @printf("  panel %4d: mu=%+ .4e  (phi(+)-phi(cp))/mu=%+ .5f  (phi(-)-phi(cp))/mu=%+ .5f\n",
            p, mu[p], ratios_plus[end], ratios_minus[end])
end
@printf("  median: (phi(+)-phi(cp))/mu = %+.5f   (phi(-)-phi(cp))/mu = %+.5f\n",
        median(ratios_plus), median(ratios_minus))
println("  interpretation: value ≈ -1 on the +n (exterior) side and ≈ 0 on the -n side")
println("  ⇒ self-limit = interior trace and exterior = interior - mu (recipe `phi -= mu` correct).")

# ------------------------------------------- (B) full-body monitor-path check
println("\n=== (B) full-body probe check (monitor influence! path) ===")
scalar_sources = (body,)   # steady semiinfinite-wake body carries its own wake

function probe_phi(positions::Vector{SVector{3,Float64}})
    probes = FastMultipole.ProbeSystem(length(positions), Float64)
    for (i, x) in enumerate(positions)
        probes.position[i] = x
        probes.scalar_potential[i] = 0.0
    end
    pnl.influence!((probes,), scalar_sources, backend;
        precalc=false, scalar_potential=true, gradient=false, hessian=(false,))
    return [probes.scalar_potential[i] + dot(Uinf, positions[i]) for i in eachindex(positions)]
end

cps = [SVector{3}(body.controlpoints[:, p]) for p in 1:n]
nhs = [SVector{3}(body.normals[:, p]) for p in 1:n]
phi_cp = probe_phi(cps)                                          # interior (self-limit) + uinf.x
phi_out = probe_phi([cps[p] + eps_rel * ell[p] * nhs[p] for p in 1:n])  # exterior side probe
phi_in = probe_phi([cps[p] - eps_rel * ell[p] * nhs[p] for p in 1:n])   # interior side probe
phi_recipe = phi_cp .- mu                                         # monitors:384-391 recipe

err_ext = phi_recipe .- phi_out    # should be O(eps)
err_int = phi_recipe .- phi_in     # should be ≈ -mu
raw_vs_int = phi_cp .- phi_in      # raw self-limit vs interior-side probe, should be O(eps)

mu_scale = median(abs.(mu))
stat(v) = (median(abs.(v)), maximum(abs.(v)))

m1, x1 = stat(err_ext); m2, x2 = stat(err_int); m3, x3 = stat(raw_vs_int)
m4, x4 = stat(err_ext ./ mu); m5, x5 = stat((err_int .+ mu) ./ mu)
@printf("median|mu| = %.4e,  offset = %.1e * sqrt(area)\n", mu_scale, eps_rel)
@printf("recipe - exterior probe : median|.| = %.3e  max|.| = %.3e   (per-mu: median %.3e, max %.3e)\n",
        m1, x1, m4, x4)
@printf("recipe - interior probe : median|.| = %.3e  max|.| = %.3e   (should be ~|mu|)\n", m2, x2)
@printf("recipe - interior + mu  : per-mu median %.3e, max %.3e     (≈0 confirms interior differs by full mu)\n",
        m5, x5)
@printf("raw phi(cp) - interior probe: median|.| = %.3e  max|.| = %.3e  (self-limit = interior trace)\n",
        m3, x3)

pass = m4 < 1e-2 && m5 < 1e-2 && (m3 / mu_scale) < 1e-2
println(pass ?
    "\nVERDICT: SIGN CORRECT — `phi -= mu` at monitors:391 yields the exterior trace." :
    "\nVERDICT: CHECK FAILED — investigate; recipe does not match exterior probe.")
