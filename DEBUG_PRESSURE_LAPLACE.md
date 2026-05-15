# PressureLaplace logic review checklist

Source: `src/FLOWPanel_simulate_monitors.jl`

## Correctness bugs (change the answer)

- [x] **`velocity_dot` is not `∂u/∂t` — it is the rate following a moving panel.** ✅ Fixed.
  Resolution: `_pressure_velocity_dot!` and `_pressure_store_negative_velocity!` now buffer `u_inertial = body.velocity + body.velocity_kinematic`, so the finite difference yields the panel-following rate `d/dt[u(x_p,t)]`. `_pressure_material_acceleration!` builds the surface gradient of `u_inertial` and uses the body-relative slip velocity `u_rel = body.velocity` (which is `u − U_body` by construction; `src/FLOWPanel_frames.jl:344-346`) for the convective term. Net result: `acceleration = d/dt[u(x_p,t)] + (u_rel·∇)u_inertial = Du/Dt`. A new `u_inertial::Vector{Matrix{Float64}}` field on the `PressureLaplace` struct holds the scratch buffer. Verified via `test/runtests_unit_postprocess.jl` (47/47 pass, including the static-body Bernoulli comparison which exercises the identity branch when `velocity_kinematic == 0`).

- [x] **RHS edge flux uses a difference, not a midpoint average.** ✅ Fixed.
  Resolution: `_pressure_rhs!` now uses the midpoint-average flux `flux = ℓ · 0.5·(a_i + a_j)·ê_ij` and accumulates `b[i] += ρ flux`, `b[j] -= ρ flux`, matching the SPD `-∇²` FV Laplacian (outward-from-i = +ê_ij, outward-from-j = -ê_ij). All 47 unit tests in `test/runtests_unit_postprocess.jl` continue to pass; the constant-field comparison is insensitive to this fix (both old and new forms give `b ≡ 0` under constancy), so the meaningful end-to-end discriminator is deferred to a non-trivial unsteady case run after issue 3 is also addressed.

- [x] **Surface gradient is taken in the inertial frame even when the body rotates.** ✅ Resolved as a side-effect of the issue-1 fix (no independent bug).
  Rationale: the gradient routine computes a single-snapshot least-squares spatial derivative using current-time `body.controlpoints` and `body.normals`, with `scalar = u_inertial` (post issue-1). That is exactly `∇u_inertial` at the moving control point's current position, which is what `Du/Dt = ∂u/∂t + (u·∇)u` requires. The "rotating stencil" concern was unfounded — the stencil is rebuilt each call from current geometry, so no time-difference is taken across the rotating stencil. The apparent frame-spin contamination came entirely from the issue-1 book-keeping mismatch: with `velocity_dot` interpreted as `∂u/∂t` and the convective velocity in the body-relative frame, the `(U_body·∇)u_inertial` correction was missing. After issue 1, the panel-following rate `d/dt[u_inertial(x_p,t)] = ∂u/∂t + (U_body·∇)u_inertial` plus `(u_rel·∇)u_inertial` telescopes to `∂u/∂t + (u·∇)u = Du/Dt` exactly. No gradient-routine change required.

- [x] **Convective term assumes `u·n = 0` exactly.** ✅ Fixed alongside issue 1.
  Resolution: the convective term now uses `u_rel = body.velocity`, for which `u_rel·n = 0` holds in the continuous limit on an impermeable surface. The explicit tangent projection is retained to scrub residual numerical normal component from the solver.

## Numerical / robustness issues

- [ ] **First-call transient is `u/dt`, not zero.**
  `velocity_dot` starts at zero, so the first call evaluates `(u_new − 0)/dt`. Docstring acknowledges this, but for a rotor at operating speed it can dominate the first-step pressure. Mitigation (preload `-body.velocity`) is not wired into `simulate!`.

- [x] **Tangent projection applied after the temporal term, not before.** ✅ Fixed.
  Resolution: the three sites that write `u_inertial`-like data (`_pressure_velocity_dot!`, `_pressure_store_negative_velocity!`, and the `u_inertial` fill at the top of `_pressure_material_acceleration!`) now subtract the normal component of `body.velocity` per panel before adding `body.velocity_kinematic`. The solver-tolerance normal residual no longer propagates into `velocity_dot` (panel-following time derivative) or into the surface gradient input. No new buffers — projection is inlined into the existing per-panel write loops at no measurable cost. The post-sum tangent projection in `_pressure_rhs!` is retained as cheap insurance for the rotating-frame residual (time-difference of tangent vectors at two different tangent planes can still have a small normal-new component). Verified that the postprocess unit test pass/error counts match the pre-change baseline (24 pass / 4 errored — all 4 errors are a pre-existing unrelated `make_octa_source_body` import issue).

- [x] **Least-squares gradient ill-conditioned on coarse / TE-adjacent panels.** ✅ Fixed.
  Resolution: `_pressure_surface_gradient!` now expands the 1-ring face-neighbor stencil to a 2-ring (neighbors-of-neighbors) whenever the panel is TE-adjacent **or** the 1-ring has fewer than 3 entries (boundary/corner panels). This mirrors the existing 2-ring pattern in `compute_mu_gradient!` (`src/FLOWPanel_postprocess.jl:146-160`) and preserves TE-edge exclusion across the second ring. Interior panels with a full 1-ring of 3 are unchanged (bit-identical path). All 47 unit tests in `test/runtests_unit_postprocess.jl` continue to pass.

- [ ] **`topology_changed` check is too loose.**
  Line 362 only checks `body.ncells` and `size(body.velocity)`. If `body.cells` is reallocated with the same `ncells` (e.g. remesh that preserves count), `m.edges` and `m.L` go stale silently when `cache=true`. The geometry signature catches metric changes but not connectivity-with-same-counts changes.

## Non-issues / confirm-only

- [ ] **Gauge column shift conditional on edge endpoints (lines 553–558).** Correct for edges touching `reference_panel`; other off-diagonal entries were zero anyway. OK.

- [ ] **`te_info` only used when `hasproperty(body, :shedding_full)`.** For a `NonLiftingBody` the stencil exclusion is a silent no-op. Intended behavior; flagging only to confirm.

## Highest-impact fixes

1. ~~**Convective frame.**~~ ✅ Done — see issue-1 resolution above.
2. ~~**RHS edge flux.**~~ ✅ Done — see issue-2 resolution above.
3. ~~**Surface gradient frame.**~~ ✅ Resolved by issue-1 fix; not an independent bug.

All three correctness bugs (issues 1, 2, 3) are now fixed. Issue 4 (convective term assumed `u·n = 0`) was also resolved alongside issue 1. Remaining items in the checklist are numerical-robustness / confirm-only.
