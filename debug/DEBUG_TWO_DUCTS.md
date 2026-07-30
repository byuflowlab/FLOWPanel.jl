# Debug Log: `examples/two_ducts_vpm.jl`

## Goal
Achieve flow tangency (RMS of U·n) within 5% of freestream on each duct.

## Setup
- Body 1: `RigidWakeBody{Union{ConstantSource, VortexRing}, 2}` (lifting, with `PanelWake`)
- Body 2: `NonLiftingBody{ConstantSource}` (non-lifting, no wake)
- Body 2 offset in z by 1.5 diameters
- AOA = 5°, Vinf = 30 m/s
- 1440 panels per body, NDIVS_theta=20, n_rfl=6

---

## Identified Issues (Initial Analysis)

| # | Issue | Severity | Location |
|---|-------|----------|----------|
| 1 | Missing wake/external potential in Dirichlet RHS | High | `FLOWPanel_solver.jl:447-485` |
| 2 | Multi-body cross-influence missing `scalar_potential` for Dirichlet | High | `FLOWPanel_solver.jl:533-538` |
| 3 | `Das` magnitude scaled by 0.005 (too small) | Medium | `two_ducts_vpm.jl:97` |
| 4 | `CPoffset=1e-12` too small for Dirichlet interior CPs | Medium | `two_ducts_vpm.jl:54,73` |

User requested fixing issues 3 and 4 first.

---

## Test 1: Fix Das magnitude + CPoffset=1e-4 (31 steps)

**Changes:**
- `Das` offset: removed `* 0.005`, now uses unit direction (matching `two_ducts.jl`)
- `CPoffset`: changed from `1e-12` to `1e-4` for both ducts
- n_steps = 31

**Results:**
```
Body 1: RMS(U·n) = 1.0079,  max|U·n| = 5.1068,  RMS/Vinf = 3.36%    ✓ (<5%)
Body 2: RMS(U·n) = 0.2131,  max|U·n| = 0.8749,  RMS/Vinf = 0.71%    ✓ (<5%)

Body 1 relative force error: 53.65%
Body 2 relative force error: 1036.84%
Body 2 max |Cp|: 2.78
```

**Flow tangency: PASSED.** Force errors very large.

---

## Test 2: CPoffset=1e-10 (user requested), 11 steps

**Changes:**
- `CPoffset`: changed to `1e-10` (user preference)
- n_steps reduced to 11 (user change for faster iteration)

**Results:**
```
Body 1: RMS(U·n) = 0.9751,  max|U·n| = 4.8753,  RMS/Vinf = 3.25%    ✓
Body 2: RMS(U·n) = 0.2113,  max|U·n| = 0.8722,  RMS/Vinf = 0.70%    ✓

Body 1 relative force error: 53.51%
Body 2 relative force error: 1049.89%
Body 2 max |Cp|: 2.85
```

**Flow tangency: PASSED.** CPoffset change from 1e-4 to 1e-10 made negligible difference.

---

## Test 3: DirectBackend (isolate FMM accuracy), 11 steps

**Changes:**
- Replaced `FastMultipoleBackend` with `DirectBackend()`
- Everything else same as Test 2

**Results:**
```
Body 1: RMS(U·n) = 0.9831,  max|U·n| = 4.9421,  RMS/Vinf = 3.28%    ✓
Body 2: RMS(U·n) = 0.2145,  max|U·n| = 0.88,    RMS/Vinf = 0.71%    ✓

Body 1 relative force error: 53.97%
Body 2 relative force error: 1046.86%
Body 2 max |Cp|: 2.83
```

**Conclusion: FMM vs Direct gives identical results.** FMM is NOT the source of force errors.

---

## Investigation: "Scalar potential for vector potential" Warnings

These warnings come from `FastMultipole.fmm!` at line 798. Traced the source:

- **NOT** from `simulate!` lines 150/162 — those use `scalar_potential=false`
- **FROM** `BackslashDirichlet.solve!` line 467: `influence!(self, self, backend; scalar_potential=true, velocity=false)`
- At that point, VortexRing strength is zero (set at line 464), so the warning is benign — zero-strength VortexRing contributes zero potential.
- Confirmed by Direct vs FMM comparison: identical results.

---

## Investigation: Dirichlet Solver Logic

Verified `BackslashDirichlet.solve!` for `RigidWakeBody` (`FLOWPanel_solver.jl:447-485`) against the correct Morino formulation:

| Step | Expected | Code | Status |
|------|----------|------|--------|
| 1. Velocity at CPs = freestream + wake + other bodies | `self.velocity` accumulated in simulate! and multi-body iteration | Lines 449, 531-539 | ✓ |
| 2. σ = -U·n | Loop over d=1:3, accumulate `-velocity[d,:] .* normals[d,:]` | Lines 458-464 | ✓ |
| 3. φ_source at interior CPs | Zero potential, `influence!(self, self, backend; scalar_potential=true)` | Lines 466-467 | ✓ |
| 4. G_doublet * μ = -φ_source | `ldiv!(strength[:,2], Glu, rhs)` | Line 480 | ✓ |
| 5. Restore velocity/potential | Saved/restored around solve | Lines 449-450, 483-484 | ✓ |

**Conclusion: Solver logic is correct.** Issues 1 and 2 from initial analysis were wrong — the Morino formulation does NOT require external potential in the Dirichlet RHS. The source strength σ = -U·n already captures external velocity, and the doublet equation is purely self-referential.

---

## Investigation: dphidt / Unsteady Bernoulli

`simulate!` uses `scalar_potential=false` at both influence calls (lines 150, 162). This means:
- `body.potential` is never updated during the loop (stays at 0)
- `dphidt = (-potential_old + potential_new) / dt = 0` always
- The unsteady Bernoulli correction `Cp -= 2*dphidt/Uref^2` is always zero
- Forces are computed purely from `Cp = 1 - (U/Uref)^2`

This may be intentional (velocity-based Cp only) or a missing feature (no unsteady pressure correction).

---

## Investigation: Force Reference Values

The reference in `two_ducts_vpm.jl`:
```julia
F1_ref = [1.5305, -109.619, 263.507]   # "BackslashDirichlet steady-state"
F2_ref = [-1.3068, 1.3759, 9.7267]
```

**Note:** `two_ducts.jl` currently uses `NonLiftingBody` for BOTH bodies (no lifting, no wake). The reference label says "BackslashDirichlet steady-state" which implies a lifting body1 — this may be from an older version of `two_ducts.jl` or a separate steady-state run.

The current unsteady sim with 11 steps (~0.07s, ~2 convective times) is far from steady state. The wake is only ~2m long on a ~1m chord duct.

---

## Current State of `two_ducts_vpm.jl`

Changes from original:
- `Das`: unit direction (removed `* 0.005`)
- `CPoffset`: `1e-10` (was `1e-12`)
- `n_steps`: 51 (was 31, then 11, now 51 — pending test)
- Backend: `FastMultipoleBackend` (restored after DirectBackend test)
- Added flow tangency diagnostic at end of file

---

## Open Questions

1. **Convergence**: Will forces approach the reference with more time steps (50-100+)?
2. **dphidt**: Should `simulate!` compute `scalar_potential=true` to enable the unsteady Bernoulli correction?
3. **Reference validity**: Are the F1_ref/F2_ref values from a lifting-body steady-state solve, or from the current non-lifting `two_ducts.jl`?
