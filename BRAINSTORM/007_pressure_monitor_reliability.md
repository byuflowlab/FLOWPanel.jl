# 007 — Reliability of Bernoulli vs Laplace(lamb) vs Laplace(Du/Dt) surface-pressure monitors

**Goal.** The rotor-hover limit cycle (dji9443, RPM 6000, revs 30–39 of
`data/rotor_hover_relaxfilter2p0_ws[_ext]`) gives three cycle-mean CTs from three pressure
monitors — Laplace material-derivative 0.0713 ± 0.0015, Bernoulli (steady) 0.0685 ± 0.0013,
Laplace lamb-vector 0.0665 ± 0.0012 — against experiment 0.072. The ~4.8e-3 cross-monitor
spread exceeds the Bernoulli deficit, so before attributing the gap to physics we need a
principled basis for trusting one monitor over the others. This document records (a) the
theory analysis of what each monitor actually computes, (b) four numerical tests that each
isolate one discrepancy mechanism, and (c) a ranked reliability verdict.

## Theory (verified against code, 2026-07-09)

### Code conventions

- `kinematic_velocity!` (`src/FLOWPanel_frames.jl:338-353`) computes rigid-body point velocity
  in global axes and *subtracts* it from `body.velocity` while accumulating it into
  `body.velocity_kinematic`. Hence `body.velocity` = fluid velocity **relative to the moving
  panel** in inertial axes; inertial fluid velocity u = `body.velocity + body.velocity_kinematic`.
- `PressureBernoulli` (pressure assembly in `calcfield_P!`, `src/FLOWPanel_postprocess.jl:934-983`):
  p = ½ρ(U∞² − |body.velocity|²), plus −ρ·phi_dot only when `unsteady=true`. The phi_dot path
  (`_pressure_bernoulli_phi_dot!`, monitors.jl) probes the scalar potential of body+wake-panel
  sources (particles excluded — vector-potential only) via FastMultipole `influence!`, forms the
  panel-following FD Dφ/Dt = (φ_new−φ_old)/dt (history zero-init ⇒ first-call transient), then
  converts to Eulerian ∂φ/∂t by subtracting `velocity_kinematic·(body.velocity+velocity_kinematic)`.
- `PressureLaplace` (monitors.jl): FV surface Poisson Lp = b per body, CG (Krylov), gauge pinned
  at `reference_panel`. Both RHS forms build the convective term from **two-point edge finite
  differences of `body.velocity`** — the Hessian-based "material acceleration"
  (`_pressure_material_acceleration!`) is computed but **never used** in the solved RHS
  (`_pressure_rhs_from_acceleration!` is dead code). The ∂u/∂t edge term (velocity_dot midpoint
  FD) is **identical** between forms.
  - `:material_derivative` (`_pressure_rhs_from_edge_material_derivative!`): per edge,
    convective = `urel·Δu` with `urel` = edge-mean *tangent-projected* velocity but Δu =
    **raw un-projected** `body.velocity` difference.
  - `:lamb_vector` (`_pressure_rhs_from_lamb_vector!`): convective = Δ(½|u_t|²) + (ω×u_t)·r
    with ω = edge-mean of `body.induced_vorticity`.
- `body.induced_vorticity` = FMM `extra_outputs=3` induced vorticity (includes particle wake's
  regularized vorticity at control points; panel elements are potential flows so their off-sheet
  curl ≈ 0) **plus** bound sheet vorticity κ = n×∇ₛμ added by `_add_bound_surface_vorticity!`
  (`src/FLOWPanel_simulate.jl:57-94`, κ from `compute_mu_gradient!`).

### Bernoulli theory result

Decomposing exact unsteady Bernoulli at the moving panel, ∂φ/∂t = Dφ/Dt|panel − u_kin·∇φ; for
flow steady in the rotating frame Dφ/Dt|panel = 0 and p_exact = p_code + ½ρ|Ω×r|². The missing
rothalpy term is purely radial ⇒ **zero axial force on any closed body** (buoyancy). So for
perfectly periodic hover, steady Bernoulli **thrust** is exact; the only neglected
thrust-relevant piece is the rotating-frame unsteadiness −ρ Dφ/Dt|panel from the wake limit
cycle: cycle-mean ~0 (φ periodic at blade-fixed points), instantaneous
~ρ·(0.05·Γ)/T_passage ≈ 5–10 Pa vs ½ρ(Ωr)² ≈ 1.9 kPa at 0.75R, i.e. ~0.3–0.5% instantaneous,
≪ that in mean. T3 measures this directly.

### Lamb-vs-matderiv theory

Same continuum identity (u·∇)u = ∇(|u|²/2) + ω×u, same Laplacian/gauge/∂u/∂t term; the
discrepancy is entirely the convective RHS, from three mechanisms:

1. **Vorticity injection** — matderiv uses no vorticity model; lamb injects extra-output ω
   (particles' regularized vorticity at the surface — nonzero where Gaussian cores overlap the
   blade, e.g. tip vortex under the following blade) + κ. Sub-issues: (a) true ω = 0 at control
   points of a potential flow; (b) prior repo work
   (`examples/particle_vorticity_curl_vs_basis_check.jl`) established curl(u_σ) =
   *Leray-projected* vorticity ≠ Σ Γ_p ζ_σ, differing by Hess(G_σ)Γ terms (2/3·ζ(0)Γ at particle
   centers) — "kernel vorticity inconsistent with divergence-free velocity"; which of these
   `extra_outputs=3` returns is empirically testable (T4); (c) full-κ at an off-sheet control
   point (vs κ/2 vs 0) is a modeling ambiguity matderiv never faces.
2. **Projection asymmetry** — matderiv's raw Δu picks up curvature/normal-rotation terms;
   lamb tangent-projects everything.
3. **Plain O(h) discretization asymmetry** between the two edge stencils, even with perfect ω.

## Test matrix

| # | Test | Isolates |
|---|---|---|
| T4 | Step 1200: at 7288 control points compare (i) FMM extra-output ω from particle sources only, (ii) direct Σ Γ_p ζ_σ basis sum, (iii) curl from FMM Hessian antisym part; check (i)−(ii) vs analytic Leray/Hess(G_σ)Γ correction | What ω actually feeds lamb; div-free consistency |
| T1 | Both Laplace forms on wake-free analytic steady bodies (sphere + Weber swept wing), mesh sweep | Mechanisms 2+3 alone (ω ≈ 0 wake-free) |
| T2 | Rotor replay steps 1152:1223, lamb-ingredient swap variants in one pass: (a) ω−κ, (b) κ only, (c) Hessian-curl ω (+κ), (d) matderiv w/ projected Δu, + baselines | Attributes the 0.0665↔0.0713 CT gap term-by-term |
| T3 | Same pass: `PressureBernoulli(unsteady=true)` + steady twin | Actual dφ/dt CT contribution; thrust-neutrality prediction |

Outputs under `data/rotor_hover_pressure_reliability/`.

## Results

### T4 (2026-07-09): what vorticity actually feeds the lamb form — resolved

Script: `examples/rotor_hover_particle_vorticity_consistency.jl`, step 1200 of
`rotor_hover_relaxfilter2p0_ws` (39 689 particles, gauserf kernel, 7288 control points).
Radial-profile CSV: `data/rotor_hover_pressure_reliability/T4_particle_vorticity_radial_profile_step1200.csv`.

Field magnitudes (|ω| over control points, s⁻¹):

| field | rms | max | median |
|---|---|---|---|
| (i) extra_outputs=3, FMM backend | 93.22 | 308.6 | 88.5 |
| (i') extra_outputs=3, direct backend | 93.20 | 308.3 | 88.4 |
| (ii) basis sum Σ Γ_p ζ_σ | 88.87 | 539.6 | 78.3 |
| (iii'/iii) Hessian-curl, direct/FMM | 93.20 | 308.3 | 88.4 |
| analytic Σ Hess(G_σ)Γ correction | 49.86 | 349.2 | 18.8 |

Differences (rms, relative to rms|ω_basis| = 88.87):

| difference | rms | rel |
|---|---|---|
| (i')−(iii') extra-output vs true curl | 6.2e-13 | 0.0000 |
| (iii')−(ii)−corr (Leray gap vs analytic) | 7.5e-10 | 0.0000 |
| (iii')−(ii) curl vs basis sum | 49.86 | **0.561** |
| (i)−(i') FMM approximation error | 3.16 | 0.036 |
| (iii)−(iii') FMM Hessian approx error | 1.7e-5 | 0.0000 |

**Verdicts.**
1. `extra_outputs=3` returns **exactly curl(u_σ)** — the Leray-projected
   (divergence-free-consistent) vorticity — machine-identical to the antisymmetric part of
   the velocity Hessian. It is *not* the naive blob sum Σ Γ ζ_σ. The user's
   div-free-consistency concern is settled in the code's favor: the lamb RHS ingests a
   vorticity that is consistent with the velocity field it is paired with. T2 variant (c)
   (Hessian-curl swap) is therefore predicted to be a no-op vs baseline (up to FMM error).
2. The curl-vs-basis distinction is *large* here (56% rms), fully explained by the analytic
   Hess(G)Γ term — the tip-vortex cores overlapping the blades put substantial regularized
   vorticity (~90 s⁻¹ rms, ~310 s⁻¹ max) on the surface. Injecting this into the lamb RHS is
   a genuine physical-modeling choice, and it is the *whole* ω the lamb form uses (mechanism
   1 of the theory section is live, and T2's ω-swap variants quantify its CT impact).
3. FMM approximation error on the extra-output vorticity is ~3.6% rms of the field —
   non-negligible but an order below the modeling questions.

### T1 (2026-07-09): form asymmetry on wake-free analytic bodies — mechanisms 2+3 are negligible

Script: `examples/pressure_form_asymmetry_check.jl`; CSV:
`data/rotor_hover_pressure_reliability/T1_form_asymmetry_sweep.csv`.
Deviation from plan: `examples/helper_functions.jl` (Weber RAE101 wing) needs GeometricTools,
which is not installed in this workspace; used an equivalent analytic loft instead — 45°-swept
untapered NACA0012, AR=5, AOA=4.2°, open tips (CL converges 0.248→0.253→0.255, sane vs Weber
0.238 given the different airfoil). Sphere is a lat-long triangulation with exact Cp. Both
solved with `steady!` + `Backslash` + `DirectBackend`, both Laplace monitors in one pass,
`grad_mu_options=(; basis=:tri)` (the steady! quad default has empty TE stencils on a sharp
analytic TE).

| case | ncells | rms\|ω\| fed to lamb | (p_mat−p_lamb) L2/q | e_mat L2/q vs exact | dCL |
|---|---|---|---|---|---|
| sphere (4,16) | 128 | 0 | 9.8e-16 | 8.8e-2 | — |
| sphere (8,32) | 512 | 0 | 1.3e-15 | 4.5e-2 | — |
| sphere (16,64) | 2048 | 0 | 1.9e-15 | 2.9e-2 | — |
| wing (4,16) | 1536 | 27.6 | 1.6e-3 | — | +4e-5 |
| wing (8,32) | 6144 | 44.9 | 4.2e-3 | — | −6e-5 |
| wing (16,64) | 24576 | 71.3 | 4.2e-3 | — | −4e-5 |

**Verdicts.**
1. **The two RHS forms are algebraically identical when impermeability holds exactly.**
   ½(u_t,i+u_t,j)·(u_j−u_i) = ½(|u_t,j|²−|u_t,i|²) whenever the raw Δu has no normal
   content; the residual difference is `urel·(u_n,j n_j − u_n,i n_i)`. On the sphere
   (direct solve, u·n ≈ 1e-14) the forms agree to 1e-15·q — mechanisms 2+3 (projection and
   stencil asymmetry) are *not* independent error sources, they collapse to the
   normal-residual term.
2. On the lifting wing the lamb form additionally ingests κ (rms 27–71 s⁻¹, a genuine sheet
   quantity growing with mesh sharpness); the resulting pressure-field difference is local
   (L∞ up to 0.1q near the TE) but integrates away: |ΔCL| ≤ 6e-5, i.e. **0.02% of CL** —
   three orders below the rotor's 7% CT spread.
3. Sphere absolute error of both forms converges under refinement (0.088→0.045→0.029 ·q),
   satisfying the plan's verification item.
4. **Implication for the rotor:** the 0.0665↔0.0713 CT gap cannot come from the convective
   discretization asymmetry itself; the live suspects are (a) the particle-wake vorticity
   injected into the lamb RHS (T4 showed ~90 s⁻¹ rms on the surface) and (b) the u·n solver
   residual entering matderiv's raw Δu (the rotor solve is iterative, not machine-exact).
   T2 decomposes both.

### T2 gate (2026-07-09): replay κ-basis bug found and fixed; lamb's deficit traced to quad-κ

Gate (steps 1152:1154, `examples/rotor_hover_pressure_reliability_replay.jl`,
`RELIABILITY_MODE=gate`): recomputed matderiv reproduced the saved CT to 5e-12 and Bernoulli
to 0, but lamb was off by 3.9e-2 relative. Diagnosis: `simulate!` accumulates the bound-sheet
κ into `induced_vorticity` with the **quad** μ-gradient basis (`_steady_aerodynamics!`
normalizes `grad_mu_options` with `default_basis=:quad`), while replay's
`_recompute_replay_fields!` called the bare `_add_bound_surface_vorticity!`, whose default is
**tri** — a genuine replay↔simulate inconsistency. Swapping quad-κ in reproduces the saved
lamb CT to 2e-5. **Fixed in src**: `replay()` now takes `grad_mu_options` normalized with
`default_basis=:quad`, mirroring `simulate!` (`src/FLOWPanel_replay.jl`); `PressureLaplace`
gained `kappa_basis::Symbol=:quad` so the `lamb_vorticity` decomposition variants subtract
the same κ that was added. Gate now PASSES (worst dev 2e-5).

The gate numbers already resolve the attribution (3-step means):

| variant | CT | reading |
|---|---|---|
| lamb baseline (quad κ) | 0.06446 | = original saved lamb |
| lamb, κ removed (`:no_bound`, wake-ω only) | 0.06711 | quad-κ injection = **−2.65e-3** |
| lamb, wake-ω removed (`:bound_only`) | 0.06452 | wake-ω injection = only −6e-5 |
| lamb with tri-κ | 0.06707 | tri-κ nearly inert; the CT hit is specific to quad-κ |
| lamb, Hessian-curl ω (+quad κ) | 0.06413 | FMM extra-output vs exact curl: −3e-4 (T4-consistent) |
| matderiv (raw Δu) | 0.06983 | |
| matderiv, Δu projected off edge-normal | −1.898 | **pathological**: at sharp-TE edge pairs n_i ≈ −n_j, the edge-averaged normal degenerates and the projection removes real loading content under a huge conormal weight — not a usable variant, but it shows the raw-Δu form's TE edges are where the forms genuinely differ |
| Bernoulli steady | 0.06723 | |
| Bernoulli unsteady | 0.12495 | **broken for particle wakes**: `_pressure_bernoulli_phi_dot!` converts Dφ/Dt\|panel to ∂φ/∂t by subtracting u_kin·∇φ but approximates ∇φ with the *total* velocity (`body.velocity + velocity_kinematic`), which in hover is dominated by the particles' solenoidal (vector-potential) field that the probed φ does not contain — a systematic +0.058 CT, not the ≲1e-4 rotating-frame unsteadiness T3 set out to measure |

### T2+T3 full pass (2026-07-09): 72 steps (1152:1223, revs 32–34), gate PASS over all steps (worst 4.5e-5)

CSV: `data/rotor_hover_pressure_reliability/T2_T3_variant_CT_vs_step.csv`. Cycle means over
the 2-rev window (note: the window mean differs from the 10-rev means quoted in §goal — the
limit cycle has a ~9-rev period — so compare variants *within* this table, not against the
10-rev numbers):

| variant | mean CT | std | reading |
|---|---|---|---|
| laplace_lamb (baseline, quad κ) | 0.06570 | 0.00093 | = original monitor |
| laplace_lamb_nobound (wake-ω only) | 0.06748 | 0.00095 | quad-κ injection = **−1.78e-3** |
| laplace_lamb_boundonly (quad κ only) | 0.06578 | 0.00094 | wake-ω injection = **−8e-5** |
| laplace_lamb_trikappa | 0.06746 | 0.00094 | tri-κ inert (−2e-5); the deficit is quad-κ-specific |
| laplace_lamb_hesscurl | 0.06539 | 0.00093 | FMM extra-output vs exact curl: −3.1e-4 |
| laplace_matderiv | 0.07043 | 0.00130 | +2.9e-3 above ω-free lamb |
| laplace_matderiv_proj | −1.245 | 0.570 | pathological (degenerate TE edge normals), discard |
| bernoulli (steady) | 0.06760 | 0.00096 | **agrees with ω-free lamb Laplace to 1.2e-4** |
| bernoulli_unsteady | 0.12603 | 0.00192 | broken for particle wakes (see gate note) |

Supporting diagnostic — the only structural difference between matderiv and lamb (T1) is the
normal-velocity content of the raw edge Δu. On the loaded step-1200 surface velocity:
median |u·n|/|u| = 1.3% (healthy bulk), but q95 = 34%, q99 = 93%, with a tail of hub-root
panels (r/R ≈ 0.065) carrying unphysical |u| up to 1.2e4 m/s (tip speed is 75 m/s) and
|u·n| up to 4.6e3 m/s, plus tip/TE panels at ~350 m/s. The +2.9e-3 that matderiv adds over
the ω-free lamb form and Bernoulli rides on this tail — an artifact of the
kerneloffset-smoothed post-solve velocity evaluation, not resolved physics.

## Ranked reliability verdict

1. **Bernoulli (steady) — most reliable for cycle-mean hover thrust.**
   Theory: for flow periodic in the rotating frame, the term it omits (½ρ|Ω×r|² rothalpy) is
   purely radial and exactly thrust-neutral on a closed body; the wake-limit-cycle
   Dφ/Dt|panel is ~0.3–0.5% instantaneous and ~0 in cycle mean. Numerics: the independent
   FV-Poisson (lamb, ω-free) monitor agrees with it to 1.2e-4 CT — two very different
   discretizations landing together. (T3's direct measurement failed for implementation
   reasons — see below — but the theory argument stands on its own.)
2. **Laplace material-derivative — plausible but inflated.** Algebraically identical to the
   lamb form under exact impermeability (T1); its +2.9e-3 excess is exactly the u·n content
   of the raw edge Δu, dominated by a small tail of near-singular hub-root/tip panels. Its
   10-rev mean 0.0713 sitting closest to experiment 0.072 is therefore **partly fortuitous**
   and should not be read as the physical answer closing the gap.
3. **Laplace lamb-vector (as configured) — least defensible.** Its deficit relative to its
   own ω-free variant is entirely the quad-basis bound-sheet κ injection (−1.8e-3), which is
   discretization-fragile (tri-basis κ: −2e-5) and conceptually ambiguous (full κ at an
   off-sheet control point vs κ/2 vs 0 — mechanism 1c). The wake-particle ω it ingests is
   div-free-consistent (T4) but CT-inert (−8e-5).
4. **Bernoulli (unsteady) — do not use with particle wakes** until
   `_pressure_bernoulli_phi_dot!` is fixed: its u_kin·∇φ conversion approximates ∇φ with the
   total velocity, which includes the particles' solenoidal field absent from the probed φ
   (+0.058 CT systematic here).

**Implication for the reported CT.** The defensible monitor estimate for the settled 2.0R
cycle is the steady-Bernoulli / ω-free-Laplace pair: 10-rev mean **0.0685 ± 0.0013** vs
experiment 0.072. The remaining −0.0035 (−5%) is a genuine physics/discretization gap
(mesh, wake model, relaxation filter) — the monitor spread no longer brackets it, and the
earlier reading that "Laplace Du/Dt nearly closes the gap" should be retired.

**Code changes that came out of this study** (all tested, 294 postprocess + 89 replay unit
tests green):
- `replay()` now accumulates κ with the same quad-basis default as `simulate!` and exposes
  `grad_mu_options` (bug fix — was tri, silently shifting replayed lamb CT by 4%).
- `PressureLaplace` gained `lamb_vorticity ∈ (:induced, :no_bound, :bound_only,
  :hessian_curl)`, `kappa_basis ∈ (:quad, :tri)`, `project_edge_du::Bool` (diagnostic).
- `_bound_surface_vorticity!` / `_add_bound_surface_vorticity_into!` /
  `_subtract_bound_surface_vorticity!` helpers (κ isolated from accumulation).
- Replay now refreshes kinematic sidecars when an unsteady `PressureBernoulli` runs without
  a velocity recompute (latent wrong-phi_dot bug).
