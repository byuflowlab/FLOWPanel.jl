# Pedrizzetti Efficacy at Satisfying the Overlap Residual

## Context

Items 008 / 010 (M2, M6) pose an **overlap-aware strength-update** for the
particle wake: the consistent strength rate `\dot\Gamma` (or projected coefficient
`\Gamma_{\mathrm{eff}}`) solves an overlap system `M\dot\Gamma=b` (008, scalar
basis) or its block-kernel form (M6). Both require an expensive per-step CG /
block solve, and 010 showed that solve is ill-conditioned in the settled
rotor-hover wake (scalar `M` diagonal-dominance `\approx0.006`, block operator
even worse). Before investing in any such solve, we should answer a cheaper
gating question:

> **Does the relaxation FLOWVPM already runs (Pedrizzetti) approximately satisfy
> the overlap-aware equation on its own — well enough that an explicit overlap
> solve adds little?**

This is motivated by a verified fact (item 010, 2026-06-25 M12 section): Pedrizzetti
relaxation (`FLOWVPM.jl` dev v4.1.0, `src/FLOWVPM_relaxation.jl:41`) rotates each
particle strength toward the **curl-of-J** `\omega_J=\nabla\times u` (not the scalar
basis `\omega_\zeta`), magnitude-preserving:
`\Gamma_i\leftarrow(1-f)\Gamma_i+f\|\Gamma_i\|\hat\omega_J(X_i)`. Item 010 also
established that the misalignment Pedrizzetti corrects is **entirely overlap-induced**
(the `K_\sigma` self-term is parallel to `\Gamma_i`, so only the neighbor sum tilts
`\omega_J(X_i)` away from `\Gamma_i`). So Pedrizzetti is *already an overlap
correction* — just a crude, direction-only, per-particle one. The open quantitative
question is how much of the full overlap correction it captures.

Relation to the mission (rotor-hover `C_T`): if Pedrizzetti already does the
overlap job, the M2/M6 CG solve is unnecessary and the wake-modeling effort should
go elsewhere (item 012, or 004/005/006). If it does not, an overlap solve (or a
better relaxation) is a real lever. Either way this is a cheap, decisive gate.

## Proposed Path

Offline, read-only over saved particle VTP states; reuse the item-008
machinery (`examples/particle_overlap_residual.jl`: VTP loader, CSR neighbor
grid, `curl_from_J`, overlap mat-vec, convective term) and the M2 projection
(`examples/particle_m2_effective_strength_diag.jl`: analytic Hessian, `\Gamma_{\mathrm{eff}}`).

Two complementary framings (do both — they cross-check):

1. **Rate-residual framing (the 008 equation).** The consistent rate satisfies
   `M\dot\Gamma=b`. Form the *effective* per-step rate Pedrizzetti imposes,
   `\dot\Gamma_{\mathrm{pedr}}=\dot\Gamma_{\mathrm{stretch}}+\tfrac{f}{\Delta t_{\mathrm{relax}}}(\|\Gamma\|\hat\omega_J-\Gamma)`,
   and measure the relative residual `r=\|M\dot\Gamma-b\|/\|b\|` for:
   - `\dot\Gamma_{\mathrm{self}}` (self-only stretch `J\Gamma`, no relaxation),
   - `\dot\Gamma_{\mathrm{pedr}}` (stretch + Pedrizzetti increment).
   Report the residual-reduction fraction `r_{\mathrm{pedr}}/r_{\mathrm{self}}`.

2. **Projection-alignment framing (M2 `\Gamma_{\mathrm{eff}}`).** The overlap solve's
   target is `\Gamma_{\mathrm{eff}}=\Gamma+M^{-1}H\Gamma`. Apply one Pedrizzetti
   rotation to the saved `\Gamma` offline to get `\Gamma_{\mathrm{pedr}}`, and measure
   whether it moves toward `\Gamma_{\mathrm{eff}}`:
   `\|\Gamma_{\mathrm{pedr}}-\Gamma_{\mathrm{eff}}\|/\|\Gamma-\Gamma_{\mathrm{eff}}\|<1`?
   Decompose into the **direction** part (the angle `\angle(\Gamma_i,\omega_J(X_i))`,
   which Pedrizzetti targets) and the **magnitude** part (which it cannot, being
   norm-preserving) — 010 predicts Pedrizzetti captures direction but misses the
   neighbor-magnitude redistribution (neighbor Hessian `\approx3.1\times` self,
   `\cos\approx0.34`).

Both are fully offline on the existing `rotor_hover_pressure_comparison` settled
states (340:359). The per-particle `\angle(\Gamma_i,\omega_J(X_i))` distribution
(the deferred metric flagged in 010's M12 section) drops out of framing (2) and
quantifies Pedrizzetti's actual workload.

A stronger optional test: a small **A/B simulation** with relaxation on vs off (or
`nsteps_relax`/`rlxf` varied), comparing the overlap residual *trajectory* and
`C_T` — to see whether residual reduction tracks a real `C_T` difference.

## Acceptance Criteria

- A quantified residual-reduction fraction `r_{\mathrm{pedr}}/r_{\mathrm{self}}` on
  settled rotor-hover states, with the direction-vs-magnitude split and the
  `\angle(\Gamma_i,\omega_J(X_i))` distribution.
- A go/no-go verdict: does Pedrizzetti satisfy the overlap residual well enough
  (e.g. residual reduced to within the saved-`J`/FMM noise floor) that the explicit
  M2/M6 CG solve is **not** worth its cost — or does a measurable gap remain that an
  overlap solve (or improved relaxation) could close?
- Offline diagnostic script (extend the 008/M2 diagnostics; no live FLOWVPM change),
  with synthetic gates (isolated particle ⟹ Pedrizzetti idle, residual unchanged).

## Caveats

- Pedrizzetti acts on `\Gamma` as a periodic rotation every `nsteps_relax`, not as a
  continuous `\dot\Gamma`; forming the *effective* rate requires the relaxation
  cadence and factor used by the run. State assumptions explicitly.
- Saved states already reflect whatever relaxation the run used; the offline
  "apply-one-Pedrizzetti-step" test (framing 2) is the clean before/after, but a true
  on-vs-off comparison needs paired runs.
- Direction-only/magnitude-preserving by construction: do not expect Pedrizzetti to
  capture the magnitude part of `\Gamma_{\mathrm{eff}}`. The interesting number is how
  much of the *total* overlap correction is direction (capturable) vs magnitude (not).
- This is a gating study for the solve track, not a wake-physics result by itself;
  keep it offline and cheap.

## 2026-06-25 Results — Pedrizzetti captures direction, not the dominant magnitude

Implemented `examples/particle_pedrizzetti_overlap_diag.jl` (offline, read-only;
reuses the item-008 loader/CSR grid/`curl_from_J` and the M1/M2 analytic Hessian).
Synthetic gates pass: isolated particle with `ω_J∥Γ` ⟹ angle ≈ 0, self-only stretch
exactly satisfies `MΓ̇=b` (`r_self=0`), Pedrizzetti idle (`r_pedr=0`); two-particle
cloud ⟹ nonzero angle, Pedrizzetti active. Default run config (`relax_pedrizzetti`,
`rlxf=0.3`, `nsteps_relax=1`, `dt=60/RPM/NT`, RPM=6000, NT=36 ⟹ `rate_factor=1080`).

Summary over 20 settled states (`rotor_hover_pressure_comparison` 340:359,
`data/rotor_hover_pressure_comparison/particle_pedrizzetti_overlap_diag.csv`):

```text
[B] angle ∠(Γ, ω_J)            : median 0.0188 rad = 1.08°  (Γ-wt mean 4.0°, p90 8.4°)
[B] overlap-corr ∥(magnitude) fraction : median 1.000  mean 0.977   <- Pedrizzetti CANNOT reach
[A] residual self / pedr / corrected   : 0.315 / 2.09 / 2.09   (ratio pedr/self red 6.63, smp 3.96)
[C] cos(d_pedr, d_eff) vs local Γ_eff  : 0.029   frac captured 0.029
[C] ‖Γ_pedr − Γ_eff‖/‖Γ_eff‖ vs raw    : 1.000 / 1.000
[overlap] TRUE local σ/h median 4.19  (nearest-nbr h_local 0.0047 m)
[overlap] neighbors: within 2σ_local 394 | within 4σ_max cutoff 21430 | kernel-effective 181
```

**Verdict — Pedrizzetti does NOT substitute for the overlap solve.** Three
independent framings agree:

- **[B, primary] The overlap correction is ~98% magnitude, which Pedrizzetti
  structurally cannot touch.** At the settled (relaxation-on) state `Γ` is already
  nearly aligned with `ω_J` (median misalignment **1.1°**) — relaxation's *direction*
  job is essentially converged — so the leftover overlap correction
  `ω_J − (2/3)ζ(0)Γ` is almost entirely **parallel to `Γ`** (∥-fraction median 1.000,
  mean 0.977). A norm-preserving rotation cannot change magnitude, so it cannot
  address this. (Small angle and ∥-fraction≈1 are the same fact: once `Γ∥ω_J`, the
  correction is pure magnitude.)
- **[C] Pedrizzetti's rotation is nearly orthogonal to the overlap-solve
  correction.** Against the M2 local-cloud target, `cos(d_pedr, d_eff) ≈ 0.03` and one
  rotation moves `Γ` no closer to `Γ_eff` (`‖Γ_pedr−Γ_eff‖/‖Γ_eff‖ = 1.000`, same as
  raw). The two operations do different things.
- **[A] The marginal Pedrizzetti rate does not reduce the 008 residual** (ratio
  pedr/self ≈ 6.6 — it increases it). Rate-assumption-sensitive (the `1080` factor),
  so weighted below B/C, but corroborating: the increment direction is not the
  overlap-rate correction. corrected-Pedrizzetti is no better (its `/√b²` magnitude
  fix is `≈1` here because `θ` is tiny).

**Go/no-go:** Pedrizzetti is a converged *direction* fix and leaves the dominant
*magnitude* part of the overlap correction untouched — so an overlap-aware
**magnitude** correction (M2/M6, or a magnitude-carrying relaxation) is where the
unaddressed content lives. But item 010/M6 showed that magnitude solve is
**ill-conditioned/intractable** in this regime. Net across 008+010+011: the overlap
correction is real and dominantly magnitude; neither Pedrizzetti nor an explicit
solve resolves it here ⟹ the lever is **overlap reduction (item 012)**.

**Overlap re-characterization (corrects earlier "~21k neighbors").** The raw
neighbor count is taken within `cutoff = 4·σ_max`, and `σ_max≈0.052 m` is the
*global* max core size (`σ_max/σ_min≈15×`), so `cutoff≈0.21 m ≈ 1.7·R` engulfs ~half
the wake — the `~21k` count is a cutoff artifact, not overlap. Measured directly:
the **TRUE local `σ/h ≈ 4.2`** (nearest-neighbor spacing vs local `σ`), with **~394
particles within `2σ_local`** and a **kernel-weighted effective overlap ≈ 181**
(diag-dominance `≈0.007`). So the overlap is genuinely heavy (~4× nominal `σ/h≈1`,
~order-100 effective neighbors vs textbook tens), but ~100× smaller than the `21k`
cutoff count. Added `eff_overlap_sample` to the diagnostic for these metrics.

**Caveats.** Saved states are relaxation-on ⟹ angles are the residual at the relaxed
fixed point (marginal, not cumulative); framing A's rate factor is a modeling
assumption (swept via `RLXF`/`NSTEPS_RELAX`/`DT`); the local `Γ_eff` is a
local-solve reference (the global solve is ill-conditioned, per M2). The stronger
live relaxation-on/off A/B + `C_T` comparison remains an optional follow-up.
