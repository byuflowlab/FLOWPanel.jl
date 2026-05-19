# PressureLaplace logic review checklist

Source: `src/FLOWPanel_simulate_monitors.jl`

## Correctness bugs (change the answer)

- [x] **`velocity_dot` is not `∂u/∂t` — it is the rate following a moving panel.** ✅ Fixed.
  Resolution: `_pressure_velocity_dot!` and `_pressure_store_negative_velocity!` now buffer `u_inertial = tangent(body.velocity) + body.velocity_kinematic`, so the finite difference yields the panel-following rate `d/dt[u(x_p,t)]` with solver-tolerance normal leakage removed before time differencing. `_pressure_material_acceleration!` uses the analytic gradient of inertial velocity and the tangent-projected body-relative slip velocity for the convective term. Net result: `acceleration = d/dt[u(x_p,t)] + (u_rel·∇)u_inertial = Du/Dt` on the impermeable surface. The stale monitor test that expected raw Cartesian storage was updated to assert the tangent-projected buffer contract.

- [x] **RHS edge flux uses a difference, not a midpoint average.** ✅ Fixed.
  Resolution: `_pressure_rhs!` uses the midpoint-average flux `flux = ℓ · 0.5·(a_i + a_j)·ν_ij` and accumulates `b[i] += ρ flux`, `b[j] -= ρ flux`, matching the SPD co-normal FV Laplacian. The earlier difference form would apply a graph derivative to acceleration rather than represent the FV flux.

- [x] **Surface gradient is taken in the inertial frame even when the body rotates.** ✅ Resolved as a side-effect of the issue-1 fix (no independent bug).
  Rationale: the gradient path now uses the analytic FastMultipole Hessian for `∇u_induced` plus `[Ω]_×` for the rigid kinematic velocity. The moving-panel finite difference supplies `∂u/∂t + (U_body·∇)u`, and the body-relative slip supplies the remaining convective term, so the pieces telescope to `Du/Dt`.

- [x] **Convective term assumes `u·n = 0` exactly.** ✅ Fixed alongside issue 1.
  Resolution: the convective term now uses `u_rel = body.velocity`, for which `u_rel·n = 0` holds in the continuous limit on an impermeable surface. The explicit tangent projection is retained to scrub residual numerical normal component from the solver.

## Numerical / robustness issues

- [ ] **First-call transient is `u/dt`, not zero.**
  `velocity_dot` starts at zero, so the first call evaluates `(u_new − 0)/dt`. Docstring acknowledges this, but for a rotor at operating speed it can dominate the first-step pressure. Mitigation (preload `-(tangent(body.velocity) + body.velocity_kinematic)`) is not wired into `simulate!`.

- [x] **Tangent projection applied after the temporal term, not before.** ✅ Fixed.
  Resolution: the three sites that write `u_inertial`-like data (`_pressure_velocity_dot!`, `_pressure_store_negative_velocity!`, and the `u_inertial` fill at the top of `_pressure_material_acceleration!`) subtract the normal component of `body.velocity` per panel before adding `body.velocity_kinematic`. The RHS also projects acceleration per edge with the averaged edge normal used by the co-normal metric.

- [x] **Least-squares gradient ill-conditioned on coarse / TE-adjacent panels.** ✅ Superseded.
  Resolution: the least-squares surface-gradient path was removed from `PressureLaplace`. The convective term now uses the analytic FastMultipole velocity gradient accumulated in `body.velocity_gradient` plus the rigid-body kinematic Jacobian `[Ω]_×`, so gradient accuracy is no longer controlled by panel-neighbor stencil conditioning.

- [ ] **`topology_changed` check is too loose.**
  Line 362 only checks `body.ncells` and `size(body.velocity)`. If `body.cells` is reallocated with the same `ncells` (e.g. remesh that preserves count), `m.edges` and `m.L` go stale silently when `cache=true`. The geometry signature catches metric changes but not connectivity-with-same-counts changes.

## Non-issues / confirm-only

- [ ] **Gauge column shift conditional on edge endpoints (lines 553–558).** Correct for edges touching `reference_panel`; other off-diagonal entries were zero anyway. OK.

- [ ] **`te_info` only used when `hasproperty(body, :shedding_full)`.** For a `NonLiftingBody` the stencil exclusion is a silent no-op. Intended behavior; flagging only to confirm.

## Discretization accuracy issues (new — found in 2026-05-14 review against rotor mesh)

- [x] **FV Laplacian is non-K-orthogonal on triangulated surfaces.** ✅ Fixed with panel-centered co-normal FV metric.
  `_assemble_pressure_laplacian` no longer uses TPFA weights `w = ℓ_ij / d_ij` with `d_ij` the 3D centroid-to-centroid distance. It now uses the shared-edge co-normal `ν_ij` built from the edge tangent and averaged panel normal, oriented from panel `i` to panel `j`, with `w = ℓ_ij (ν_ij·r_ij) / |r_ij|²`. This keeps pressure unknowns panel-centered and makes the operator metric consistent with the RHS flux. The standard vertex cotangent Laplace-Beltrami was deliberately deferred because it is not a drop-in replacement for cell-centered `body.P`.

- [x] **RHS edge flux uses centroid-to-centroid direction, not the in-plane edge normal.** ✅ Fixed.
  `_pressure_rhs!` now uses the same oriented edge co-normal `ν_ij` as the Laplacian instead of `ê_ij = (cp_j - cp_i)/d`. The midpoint flux remains conservative: `b[i] += ρ f`, `b[j] -= ρ f`.

- [x] **LS gradient still rank-deficient on high-AR triangles with full 1-ring.** ✅ Superseded.
  The LS gradient is no longer used by `PressureLaplace`; see the analytic FastMultipole Hessian resolution above.

- [x] **Tangent projection at edges uses panel-i normal for both sides.** ✅ Fixed.
  `_pressure_rhs!` now projects each side's acceleration with the same averaged edge normal used to form `ν_ij`, then dots the midpoint value with that co-normal.

- [ ] **Single-point Dirichlet gauge is the noisiest possible.**
  Pinning `p[reference_panel] = reference_pressure` is fine for SPD but the pressure error grows with geodesic distance from the pin and the choice of pin biases the global field. Replacing with **zero-mean pressure** `Σ p_i A_i = 0` via either a Lagrange multiplier (saddle-point system) or rank-1 deflation (project the constant null space out of `b` each call and out of CG iterates) gives a uniform error distribution. Symmetric; CG-compatible via projection.

- [ ] **BDF1 temporal derivative.**
  `(u_new - u_old)/dt` is O(dt). For a rotor where ∂u/∂t carries real signal (blade-passage tones, etc.), a BDF2 stencil `(3u_new - 4u_old + u_older)/(2dt)` would help. Requires an additional `velocity_older` buffer and a first-two-step warm-up that falls back to BDF1.

## Highest-impact fixes

1. ~~**Convective frame.**~~ ✅ Done — see issue-1 resolution above.
2. ~~**RHS edge flux.**~~ ✅ Done — see issue-2 resolution above.
3. ~~**Surface gradient frame.**~~ ✅ Resolved by issue-1 fix; not an independent bug.

All three correctness bugs (issues 1, 2, 3) are now fixed. Issue 4 (convective term assumed `u·n = 0`) was also resolved alongside issue 1. Remaining items in the checklist are numerical-robustness / confirm-only.

## Next priorities (2026-05-14 review)

Ranked by expected impact on the rotor case:

1. ~~**Analytic ∇u via FastMultipole Hessian**~~ ✅ Done — replaces the LS gradient with a mesh-independent quantity.
2. ~~**Cotangent Laplace–Beltrami + edge-normal RHS**~~ — superseded by panel-centered co-normal FV operator/RHS; vertex cotangent remains a separate future discretization if pressure unknowns move off cells.
3. **Zero-mean gauge** — deferred; uniform error distribution.
4. ~~**Per-edge averaged normal in tangent projection**~~ ✅ Done as part of co-normal RHS.
5. **BDF2 temporal derivative** — deferred; only if unsteady signal matters at the dt being used.

## Capped-wing discrepancy follow-up (2026-05-16)

Context: `examples/simple_wing_capped_pressure_comparison.jl` shows
`PressureBernoulli` and `KuttaJoukowskiForce` in close agreement, while
`PressureLaplace` predicts a much larger force. One reproduced run on
`wing_ar4_naca0016_5.msh` gave:

- `PressureBernoulli`: `[-0.0193, ~0, 0.3916]`
- `PressureLaplace`: `[-0.3790, 0.0190, 12.7209]`
- `KuttaJoukowskiForce`: `[-0.0456, -0.00026, 0.3882]`

### Raw Hessian finite-difference check

Test performed: solve the capped Dirichlet wing, then compare the raw
panel-induced velocity Jacobian against central finite differences of the same
raw panel-induced probe velocity. This deliberately excludes freestream and the
`compute_mu_gradient!` doublet-jump correction, so the comparison targets only
the FastMultipole kernel Hessian / storage path.

Findings:

- At body control points, finite differences disagree with
  `body.velocity_gradient` at O(1) relative error. Median relative error over a
  25-panel sample was about `1.0` for `h = 1e-3, 3e-4, 1e-4, 3e-5, 1e-5`.
- Sign and transpose variants did not explain the mismatch. At `h = 1e-5`, the
  median relative errors were all still about `1.0` for `J`, `J'`, `-J`, and
  `-J'`.
- `body.velocity_gradient[:,:,p]` exactly matches `ProbeSystem.hessian` at the
  same control point, so this does not look like a body-buffer packing or
  row/column storage bug.
- Off the surface, the analytic probe Hessian matches finite differences well:
  - normal offset `1e-2`: median relative error `1.19e-5`
  - normal offset `3e-3`: median relative error `1.46e-5`
  - normal offset `1e-3`: median relative error `1.60e-5`
  - closer to the surface, worst-case errors grow, consistent with near-panel
    singular/discontinuous evaluation rather than a generic Hessian bug.

Interpretation: the raw kernel Hessian appears correct away from the singular
surface and is stored consistently. The control-point finite-difference failure
is not, by itself, evidence that the Hessian implementation is algebraically
wrong.

### Current strongest hypothesis

`PressureLaplace` is using `body.velocity_gradient` as the spatial gradient of
the velocity field that appears in the pressure-Poisson RHS, but
`calcfield_U!` builds the surface velocity used by `PressureBernoulli` from more
than the raw panel-induced velocity:

- raw panel/wake velocity from `influence!`
- the lifting-body surface doublet-gradient correction from `compute_mu_gradient!`
- freestream

The current Hessian path only represents the raw induced field accumulated by
`influence!(..., velocity_gradient=true)` plus any rigid-body kinematic
Jacobian. It does not obviously include the derivative of the
`compute_mu_gradient!` correction. That makes the Laplace RHS internally
inconsistent with the velocity magnitude used by Bernoulli and with the
surface velocity that Kutta-Joukowski effectively validates.

### Suggested next steps

1. **Manufactured-pressure residual.** Continue using Bernoulli pressure as a
   diagnostic pressure field:
   `r_B = L * (P_B .- P_B[ref]) - b_L`. A reproduced run gave
   `norm(r_B) / norm(b_L) ≈ 0.56`, so the main disagreement is in the assembled
   RHS/operator consistency, not only the final force integration.

2. **Restore the old LS-gradient RHS as an A/B test.** Temporarily reintroduce
   the pre-Hessian `_pressure_surface_gradient!` path from commit `1789542` and
   compare:
   - `b_hessian`
   - `b_surface_ls`
   - `L * (P_B .- P_B[ref])`

   If `b_surface_ls` is much closer to Bernoulli, the regression is the
   Hessian substitution as a model for the surface velocity gradient, not the
   sparse Poisson solve.

3. **Disable the doublet-gradient correction as a consistency test.** Run the
   pressure comparison with `calcfield_U!(...; doublet_gradient=false)` and
   recompute the Laplace RHS from the matching raw induced velocity/Hessian. If
   agreement improves substantially, the missing derivative of the
   doublet-gradient correction is confirmed as the leading issue.

4. **Compare edge-local pressure gradients.** For every shared edge, compare
   the Bernoulli pressure jump against the Laplace RHS acceleration flux:
   `(P_B[j] - P_B[i])` versus `-rho * a_edge * ds`. A global sign error should
   show up coherently; TE/cap-local spikes would point toward surface-gradient
   discontinuity handling.

5. **Decide the intended gradient definition.** Either:
   - use a surface-gradient reconstruction of the final surface velocity field
     in `PressureLaplace`, including TE one-sided treatment, or
   - derive and implement the missing spatial derivative of the
     `compute_mu_gradient!` contribution so the analytic Hessian path matches
     the velocity field used for pressure.

### Diagnostic log

- 2026-05-16: Starting with the manufactured-pressure residual diagnostic.
  Goal: quantify whether the Bernoulli pressure field satisfies the current
  co-normal Laplace operator and RHS before changing the gradient model.

- 2026-05-16: Manufactured-pressure residual result on
  `wing_ar4_naca0016_5.msh` with the current Hessian RHS:
  - `npanels = 1000`, `nedges = 1500`
  - `||b_L|| = 9.2798e6`
  - `||L * (P_B - P_B[ref])|| = 1.0848e7`
  - `||L * (P_B - P_B[ref]) - b_L|| / ||b_L|| = 0.5591`
  - cosine alignment between `L * P_B` and `b_L`: `0.8785`
  - absolute residual quantiles: median `405`, p90 `1.76e4`, p99 `5.37e5`

  The largest residual panels are concentrated near the trailing edge / span-tip
  region: top offenders include panels around `x ≈ 1.97`, `|y| ≈ 4.00`,
  `|z| ≈ 0.0027`. This is consistent with a surface-gradient/discontinuity
  issue rather than a uniform sparse-solve failure.

- 2026-05-16: Added `debug_pressure_laplace_ab.jl` to make the RHS A/B
  repeatable. It reconstructs the old pre-Hessian surface-LS convective term
  in a diagnostic-only path and compares it against the current Hessian RHS.

- 2026-05-16: Surface-LS RHS A/B result:
  - `||L * P_B - b_hessian|| / ||b_hessian|| = 0.5591`
  - `||L * P_B - b_surface_ls|| / ||b_surface_ls|| = 1.8306`
  - cosine alignment with `L * P_B`: current Hessian `0.8785`, surface-LS
    `0.4707`
  - `||b_surface_ls - b_hessian|| / ||b_hessian|| = 1.1038`
  - direct-solve force from current Hessian RHS:
    `[-0.4197, ~0, 12.5352]`
  - direct-solve force from surface-LS RHS:
    `[1.9522, ~0, 19.4428]`

  This does **not** support reverting to the old surface-LS gradient as the
  fix. The current Hessian RHS is substantially closer to the Bernoulli
  manufactured pressure than the old LS reconstruction on this capped-wing
  case. The remaining mismatch is still large and still localized near
  trailing-edge / tip geometry.

- 2026-05-16: `doublet_gradient=false` consistency test result:
  - `||b_L|| = 3.5641e6`
  - `||L * P_B|| = 1.5342e6`
  - `||L * P_B - b_L|| / ||b_L|| = 0.6859`
  - cosine alignment between `L * P_B` and `b_L`: `0.8304`
  - raw Bernoulli force: `[0.0158, ~0, 0.1711]`
  - raw Laplace/Hessian force: `[-0.4398, ~0, 4.4781]`

  Disabling the surface doublet-gradient correction does **not** improve the
  Bernoulli-vs-Laplace consistency. This weakens the hypothesis that the
  missing derivative of `compute_mu_gradient!` is the sole driver. It may still
  contribute to the physical mismatch, but the current residual also exists in
  a raw induced-velocity-only comparison.

- 2026-05-16: Edge-local pressure-flux comparison using the current full
  velocity field:
  - edge count: `1500`
  - compare each edge's Bernoulli pressure flux
    `w * (P_B[i] - P_B[j])` against the Laplace acceleration flux
    `rho * ell * 0.5 * (a_i + a_j) · nu_ij`
  - relative edge-flux norm error: `0.5063`
  - cosine alignment: `0.9450`
  - same-sign fraction: `0.7873`
  - absolute edge-flux difference quantiles: median `254`, p90 `6.13e3`,
    p99 `2.79e5`, max `2.03e6`

  This does not look like a simple global sign error. The dominant edge
  mismatches are again concentrated near `x ≈ 1.9745`, `|y| ≈ 4.00`, i.e. the
  trailing-edge / span-tip corner region. Several of the worst local pressure
  fluxes are larger than the acceleration flux by factors of roughly `2.4` to
  `2.6`; one smaller-sign edge has a ratio around `7.2`.

  Current implication: prioritize geometry/topology handling around the capped
  wing's trailing-edge/tip corner before changing the global Poisson solve.

- 2026-05-16: Localized to surface velocity, not the Laplace solve. Visualizing
  `body.velocity · body.normals` on the capped Dirichlet wing shows tangency
  residuals up to ~10 m/s (vs |Vinf|=56 m/s) concentrated on high-aspect-ratio
  sliver triangles along the cap–wing seam. The seam slivers have AR 20+ and
  alternating-sign normal residual — the classic signature of a poorly
  conditioned local reconstruction on near-degenerate panels, not a uniform bug.

- 2026-05-16: Kerneloffset sweep on the same mesh (Dirichlet, `:backslash`):
  | k_off | med rel \|u·n\| | p99 rel \|u·n\| | hi-AR med rel | rel resid \|L·P_B−b_L\|/\|b_L\| | F_L_z |
  | --- | --- | --- | --- | --- | --- |
  | 1e-2 | 0.0225 | 0.757 | 0.073 | 0.526 | 4.35e5 |
  | 3e-3 | 0.0225 | 0.753 | 0.032 | 0.622 | 4.10e5 |
  | 1e-3 | 0.0223 | 0.753 | 0.031 | 0.559 | 4.10e5 |
  | 3e-4 | 0.0221 | 0.753 | 0.030 | 0.561 | 4.10e5 |
  | 1e-4 | 0.0221 | 0.753 | 0.030 | 0.561 | 4.10e5 |

  Global tangency and force mismatch are essentially invariant across 2 orders
  of magnitude of kerneloffset. Hi-AR median improves ~2× from 1e-2 to 3e-4
  then plateaus. → Regularization mismatch between BIE matrix and velocity
  evaluation is **not** the leading driver.

- 2026-05-16: Switching the BC to Neumann (with one panel held out for
  rank-fix) drops `|u·n|` on the same sliver panels from O(10) m/s to
  O(1e-4) m/s — five orders of magnitude — confirming the BIE solve itself is
  fine. The residual pattern still shows a faint correlation with the sliver
  geometry at the 1e-4 level. Bernoulli vs Laplace agreement under Neumann
  pending follow-up.

- 2026-05-16: Mechanism hypothesis. The Dirichlet pipeline reconstructs the
  surface velocity as `U_∞ + u_induced + ∇_s μ`, where `∇_s μ` comes from
  `compute_mu_gradient!` (`FLOWPanel_postprocess.jl:95`): a 1-ring LS gradient
  of the doublet strength, constrained to the panel's tangent plane by a 1e4
  normal-projection penalty. On a sliver triangle the 1-ring neighbors line up
  nearly collinearly inside the tangent plane → the 2D Gram matrix becomes
  ill-conditioned → one tangent component of `∇_s μ` is unreliable → the
  reconstructed surface velocity picks up a spurious normal component. This
  fits all observations: (a) localized to slivers, (b) sign-flipping along the
  seam, (c) insensitive to kerneloffset, (d) fixed by switching to a BC that
  doesn't depend on `∇_s μ` for tangency. `PressureLaplace`'s convective term
  `(u·∇)u` lives off this corrupted tangent component, which is why Bernoulli
  (`|u|²`-dominated) and KJ (edge-circulation) survive while Laplace does not.

  Suggested surgical fix on the post-process side: detect AR > threshold or
  stencil tangent-plane κ > threshold and fall back to a more robust
  reconstruction (e.g. inherit `∇_s μ` from a smoother neighbor, or use a
  Green–Gauss flux around the dual cell). Long-term path: remesh capped
  geometries to avoid slivers, or document `PressureLaplace` as Neumann-only.

- 2026-05-16: Default flipped — `calcfield_U!` now defaults to
  `gradient_robust=false` (the AR mask + LS-fallback path is opt-in).
  Rationale: on the refined NACA 0012 capped wing the mask is degenerate
  (≤11 panels flagged at threshold=3, zero at threshold≥7) and on the
  coarse wing the AR-sweep showed the knob saturates well short of
  agreement (best CzL ≈ 5.0, still ~13× CzB), with catastrophic
  rel_resid spikes at intermediate thresholds (2.09 at threshold=15).
  No callers in src/ or test/ set `gradient_robust` explicitly, so this
  is purely a default flip — the path remains available via
  `gradient_robust=true` for diagnostics.

- 2026-05-16: Bernoulli vs Laplace under Neumann — the pending follow-up.
  Harness: `debug_bc_comparison.jl`. Mesh: `naca0012_nc101_nw26.msh` (refined).
  Three cases: pure Dirichlet, Dirichlet+AR=3 mask, pure Neumann (kernel
  = `pnl.VortexRing`, one mid-mesh panel dropped to break watertightness
  and resolve rank deficiency per `simple_wing_capped.jl:194`).
  `velocity_dot` zeroed throughout (single-shot harness).

  | case            | CzB    | CzL    | CzKJ    | BL_rel | BKJ_rel | max\|u·n\|/\|u\| | rel_resid |
  |-----------------|--------|--------|---------|--------|---------|------------------|-----------|
  | Dirichlet       | 0.3515 | 0.3056 | 0.2488  | 0.135  | 0.309   | 9.20e-1          | 0.940     |
  | Dirichlet+AR=3  | 0.3515 | 0.3058 | 0.2488  | 0.135  | 0.309   | 9.20e-1          | 1.109     |
  | Neumann         | 0.4490 | 1.1879 | -0.0059 | 1.645  | 1.311   | 4.81e-4          | 1.992     |

  Observations:
  - Neumann BIE actually enforces impermeability — `max|u·n|/|u|` drops
    by 3+ orders of magnitude (9.2e-1 → 4.8e-4). The dropped-panel +
    pure-`VortexRing` formulation is the right way to test the
    hypothesis on this geometry; the earlier `Union{Source,VortexRing}`
    attempt left `max|u·n|/|u|` at 0.92 (effectively still Dirichlet)
    because the coupled formulation routes surface velocity through
    `compute_mu_gradient!` the same way.
  - **The sliver-∇_s µ hypothesis is refuted as the sole driver.**
    Cleaning tangency made Laplace *worse*, not better: CzL went
    0.31 → 1.19, BL_rel went 0.135 → 1.65, rel_resid went 0.94 → 1.99.
    If contaminated `∇_s µ` were the dominant cause of the
    Bernoulli↔Laplace gap, removing the contamination should have
    narrowed it.
  - CzB itself shifted 0.35 → 0.45 between BCs (~28%), meaning the
    "trusted" Bernoulli baseline is not invariant to the surface
    velocity reconstruction path. The 13% Dirichlet BL_rel is partly a
    coincidence of where two independently-imperfect quantities happen
    to land.
  - CzKJ ≈ 0 under Neumann is a harness artifact —
    `KuttaJoukowskiForce` extracts circulation `γ` from the
    Source+Doublet representation; the pure-`VortexRing` body stores
    circulation in a different field, so the probe reads ~0. Not a
    physics result.
  - AR=3 row is uninformative on this mesh (only 11 panels flagged, no
    movement). Consistent with the threshold sweep — the AR mask is a
    coarse-mesh-only knob.

  Implication: the remaining gap is in the FV pressure operator / RHS
  itself, not in surface-velocity reconstruction. The most likely sub-
  candidates (from prior log entries) are (i) TE/tip-corner stencil
  handling in `_assemble_pressure_laplacian` (edge-flux comparison
  already localized worst residuals there), and (ii) the gauge / shed-
  edge treatment. Next diagnostic should target one of these directly
  rather than continuing on the velocity-reconstruction side.

- 2026-05-16: Open-wing data point. Resurrected `examples/sweptwing.jl`
  (Weber & Brebner 45° swept RAE 101 wing, AR=5, AOA=4.2°, single-loft
  `simplewing` body of 2880 triangular `VortexRing` panels, no tip caps,
  Backslash + DirectBackend) and added both Bernoulli and Laplace force
  integrations side by side. Experimental reference: `CLs_web[2] = 0.238`.

  ```
  CL_bernoulli = 0.2422   (error 1.77%)
  CL_laplace   = 220.5    (error 92527%)
  CLexp        = 0.238
  ```

  - Bernoulli is well within experimental tolerance — confirms the wing
    build, BIE solve, and Bernoulli post-process path are all correct on
    an open-wing geometry.
  - Laplace blows up by ~900× on a geometry with **no caps and no
    sliver triangles**. This is incompatible with the previously-leading
    "sliver `∇_s µ` contamination drives the gap" hypothesis (which was
    already refuted as the sole driver by the BC-comparison entry
    above; this is the second independent refutation).
  - The remaining most-likely failure modes are now wake-side: the open
    wing's free tip edges + shed TE feed into `body.shedding_full`, and
    the FV stencil in `_assemble_pressure_laplacian` treats free
    boundaries as zero-flux. The tip-vortex–induced velocity gradient
    near the free edges is also highly singular and likely contaminates
    `body.velocity_gradient` regardless of mesh quality.
  - Next diagnostic candidate: probe the per-panel Laplace pressure
    field on the sweptwing body and locate where it diverges. If the
    divergence is concentrated near the tip / TE free edges (likely),
    the FV operator's free-edge boundary condition is the next fix
    target. If it's distributed, suspect the Hessian convective term
    `(u·∇)u` near the wake.

- 2026-05-16: Sweptwing Laplace diagnostic — `gradient_mode` A/B + spatial
  localization. Harness: `debug_sweptwing_laplace.jl`. Same open swept
  wing (2880 panels, no caps, AOA=4.2°).

  | mode                | CL_L     | median \|P_L\| | max \|P_L\| | rel_resid |
  |---------------------|----------|----------------|-------------|-----------|
  | bernoulli (ref)     | 0.242    | 73             | 2.7e3       | —         |
  | :raw_hessian        | 220.5    | 1.5e5          | 2.5e6       | 0.999     |
  | :surface_velocity   | -10.8    | 3.7e3          | 6.3e4       | 0.879     |

  Spatial concentration of |P_L| (mean by region):

  | mode               | tip 10% | root 10% | mid 80% |
  |--------------------|---------|----------|---------|
  | :raw_hessian       | 1.2e5   | 1.5e5    | 2.9e5   |
  | :surface_velocity  | 5.6e3   | 4.7e3    | 5.1e3   |

  Top-K |P_L| panel locations:
  - :raw_hessian — clustered in **bulk** of wing (x≈0.7–1.1,
    |y|≈0.2–0.6), magnitudes ~2e6. NOT at edges.
  - :surface_velocity — clustered at **tip** edge (|y|≈1.0–1.24) and
    near TE (x≈1.4–1.7), magnitudes ~3e4–6e4.

  Interpretation:

  1. `:raw_hessian` is structurally wrong on lifting bodies. The FMM
     `body.velocity_gradient` does not include the spatial derivative
     of the `compute_mu_gradient!` (∇ₛµ) correction that
     `calcfield_U!` folds into `body.velocity`. The `(u·∇)u` term is
     thus computed with mismatched velocity and gradient, producing a
     spurious diffuse acceleration over the *entire wing surface* —
     median |P_L| is 2000× the Bernoulli median, with the largest
     panels in the bulk, not at edges. Consistent with the open-wing
     CL_L = 220 finding from the prior entry.
  2. `:surface_velocity` removes the bulk divergence (median |P_L|
     drops 40×, CL_L magnitude drops 20×, top-K location shifts from
     bulk to tip edge). This confirms the missing ∇(∇ₛµ) was the
     dominant error in `:raw_hessian` mode.
  3. The remaining `:surface_velocity` error (CL_L = −10.8 vs
     CL_exp = 0.238, factor ~45 over-prediction with sign flipped) is
     **concentrated at the wingtip**. The wing has no tip cap; the FV
     `_assemble_pressure_laplacian` therefore sees free edges at the
     tip whose flux is implicitly treated as zero — physically wrong
     (vorticity sheds off the tip into the free wake). Next likely
     root cause: the FV stencil's free-edge boundary condition. The
     `:surface_velocity` `rel_resid` of 0.879 (vs 0.999 for
     `:raw_hessian`) shows the operator is closer to admitting the
     Bernoulli pressure but still substantially inconsistent.

  Recommended next steps in priority order:
  - Audit `_assemble_pressure_laplacian` for how it treats panels with
    edges that have no neighbor (free tip / open boundary). If those
    edges are simply skipped (Neumann homogeneous), substitute a
    physically motivated free-edge BC (e.g. extrapolated tangential
    acceleration, or a Dirichlet from far-field).
  - For lifting Dirichlet bodies, document `:surface_velocity` as the
    recommended mode (and consider flipping the default for those
    body types). The `:raw_hessian` default is correct only when
    `body.velocity` *is* the raw influence-evaluated velocity, i.e.
    non-lifting / no `∇ₛµ` correction path.
  - VTK files written to `sweptwing_debug/{wing_bernoulli,
    wing_laplace_raw, wing_laplace_sv}.vtu` for ParaView inspection
    of the spatial pressure field.

- 2026-05-16: Sweptwing pressure-Poisson consistency check — new
  repeatable diagnostic script:
  `examples/sweptwing_pressure_poisson_check.jl`. This supersedes the
  prior "free-edge boundary condition is likely next" interpretation for
  the swept-wing case.

  Same open swept wing, same one-shot zero-history `velocity_dot`
  initialization as `examples/sweptwing.jl`:

  | quantity | value |
  | --- | --- |
  | panels | 2880 |
  | boundary panels | 96 |
  | interior panels | 2784 |
  | Bernoulli force | CL = 0.24221, CD = -0.00109 |
  | `surface_velocity` monitor force | CL = -10.76308, CD = -0.77737 |
  | `raw_hessian` monitor force | CL = 220.45288, CD = 21.35303 |

  Pressure ranges:

  | field | min | max | median |
  | --- | --- | --- | --- |
  | Bernoulli | -2.715e3 | 4.869e2 | -2.305e1 |
  | Laplace `surface_velocity` | -6.297e4 | 2.653e4 | 3.138e3 |
  | Laplace `raw_hessian` | -2.541e6 | 2.171e6 | -4.099e4 |

  Manufactured-pressure residuals using Bernoulli pressure as the
  diagnostic field:

  | RHS mode | `||L*pB||` | `||rhs||` | `||L*pB-rhs||/||rhs||` | direct CL |
  | --- | --- | --- | --- | --- |
  | `surface_velocity` | 3.0105e4 | 7.9216e4 | 0.8788 | -10.889 |
  | `raw_hessian` | 3.0105e4 | 3.9449e5 | 0.9989 | 220.262 |

  Boundary split for `surface_velocity`:

  | panel set | relative residual |
  | --- | --- |
  | boundary panels | 0.762 |
  | interior panels | 0.915 |

  This rules out a boundary-only explanation. The residual is at least as
  large on fully interior panels, so skipped open edges / homogeneous
  Neumann treatment are not the leading failure in the current swept-wing
  mismatch.

  Sparse solve check:

  | solve | relative residual |
  | --- | --- |
  | CG monitor result | 4.96e-3 |
  | direct `L \ rhs` | 3.74e-15 |

  The direct solve gives the same bad force, so this is not a Krylov
  convergence problem. It is the RHS that is inconsistent with the
  pressure field under the current discrete operator.

  Convective variant check with the `surface_velocity` gradient:

  | variant | `||L*pB-rhs||/||rhs||` |
  | --- | --- |
  | `G*u` | 0.8788 |
  | `transpose(G)*u` | 0.8686 |

  Transposing the reconstructed gradient marginally improves the norm but
  does not change the diagnosis.

  Mimetic Bernoulli edge-flux experiment:

  - Build an RHS directly from Bernoulli pressure jumps with the same
    co-normal weights used by `L`.
  - Solve the same sparse system.

  Results:

  | check | value |
  | --- | --- |
  | `||L*pB-b_mimetic||/||b_mimetic||` | 1.49e-16 |
  | `||p_mimetic-pB||/||pB||` | 1.11e-13 |
  | mimetic force | CL = 0.24221, CD = -0.00109 |

  Interpretation:

  - The panel-centered co-normal pressure operator can reproduce the
    Bernoulli pressure field and force to roundoff when it receives a
    flux RHS that is mimetic with its own weights.
  - Force integration is also consistent: the mimetic solve produces the
    Bernoulli CL/CD exactly.
  - Gauge is not the issue; the constant-pressure force contribution is
    negligible on this mesh and the mimetic shifted pressure reproduces
    the same force.
  - The remaining failure is the acceleration-to-edge-flux RHS:
    `_pressure_material_acceleration!` plus `_pressure_rhs!` does not
    produce the discrete divergence that corresponds to the Bernoulli
    pressure field.

  Updated recommendation:

  1. **Do not start by changing the Laplacian, gauge, free-edge boundary
     handling, force integration, or CG settings.** The mimetic RHS
     experiment shows those pieces can produce the right answer on the
     swept wing.
  2. **Replace the current pointwise-acceleration RHS with a mimetic
     edge-based RHS.** The immediate target should be an edge formula
     that uses the same co-normal weights as `_assemble_pressure_laplacian`
     and is algebraically consistent with the discrete Bernoulli pressure
     jump for steady potential flow. For the steady term, this likely
     means computing the edge-normal convective pressure flux from a
     kinetic-energy jump,
     `p_i - p_j = 0.5*rho*(|u_j|^2 - |u_i|^2)`, instead of reconstructing
     `(u_t · ∇_s)u` at panel centers and midpoint-averaging it to the edge.
     The unsteady term should receive the same treatment using a discrete
     potential-time derivative jump when available.
  3. **Add a regression test before the fix.** A small manufactured
     pressure/velocity case should assert that the PressureLaplace RHS
     satisfies `norm(L*p_B - b)/norm(b)` to a tight tolerance for a
     steady potential flow where Bernoulli is known. The swept-wing
     diagnostic can remain an example-level integration check because it
     is too expensive for unit tests.
  4. **Keep `:raw_hessian` out of the fix path for lifting bodies.** It is
     still structurally inconsistent with the final surface velocity on
     lifting bodies because `body.velocity_gradient` does not include the
     derivative of the `∇_s μ` correction folded into `body.velocity`.
     The repair should first make `:surface_velocity` mimetic, then decide
     whether `:raw_hessian` should be documented as non-lifting only or
     upgraded with the missing correction derivative.
