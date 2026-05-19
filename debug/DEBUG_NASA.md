# Debug Log: Flow Potential Residuals ~0.6

## Problem
`flow_potential_residuals` returns ~0.6 instead of ~0 after `BackslashCoupled` solve in `test/runtests_unit_solver.jl`.

## Investigation Log

### Step 1: Reproduce baseline (single body, nasa_surface_spaced.msh)
- Flow Tangency Residuals: [3302.27]
- Flow Potential Residuals: [0.6019]

### Step 2: Check influence! dispatch
- Method at `solver.jl:904` (`influence!(::Tuple, ::Tuple, backend=DirectBackend(); optargs...)`) was identified as a suspect due to undefined `precalc` variable and potential kwargs conflict.
- **Finding**: The method IS matched by `which()` but NOT actually called. Julia's default-argument trampoline redirects the 2-arg call to a 3-arg call `influence!(targets, sources, DirectBackend(); ...)`, which then dispatches to the more specific `fmm.jl:93` method (`influence!(::Tuple, ::Tuple, ::DirectBackend; ...)`).
- **Conclusion**: The solver.jl `influence!` override is dead code for DirectBackend calls. NOT the cause of the bug. (Note: the `precalc` reference at line 906 is an undefined variable bug that should be fixed separately, but it's never reached.)

### Step 3: Linear system verification
- Added `G_copy * sol - rhs` check using a copy of G before LU factorization.
- **Finding**: Linear residual is ~1e-14 (machine precision). The LU solve is CORRECT.
- G diagnostics: size=113x113, cond~3447, diagonal in [-0.52, 0.55]

### Step 4: CP offset sensitivity analysis (single body, CPoffset=1e-14)
Tested `flow_potential_residuals` at various CP offsets (all interior/negative):
| Offset | Residual |
|--------|----------|
| -1e-14 (same as solve) | 6.7e-30 (essentially 0) |
| -1e-12 | 0.6019 |
| -1e-10 | 0.6019 |
| -1e-6  | 0.6018 |
| -1e-2  | 0.5465 |

**Key finding**: The residual is 0 ONLY at the exact solve offset (-1e-14). Any other interior offset gives ~0.6.

### Step 5: CPoffset sensitivity — single body
Tested with CPoffset=1e-3 and CPoffset=1e-2 in `generate_body`:
- CPoffset=1e-3: residual at -1e-10 = 4.7e-6 (vs 0.6 with 1e-14)
- CPoffset=1e-2: residual at -1e-10 = 7.4e-4

**Conclusion**: CPoffset=1e-14 is too small — CPs are at the doublet singularity. A reasonable value (1e-3 to 1e-2) works much better for single body.

### Step 6: Multi-body test — solution blow-up
With 2 bodies (nasa_wing.msh + nasa_surface_spaced_repaired.msh) and CPoffset=1e-14:
- G condition: 310,872
- Body 1 doublet strengths: ALL ≈ -5.75e15 (blow-up!)
- Linear residual: 11.37 (solve FAILED)

With CPoffset=1e-2:
- G condition: 746,096 (even worse!)
- Body 1 doublets: ALL ≈ -9.7e15 (still blow-up!)

### Step 7: G matrix structure analysis — ROOT CAUSE FOUND
With CPoffset=1e-2 and `flip_normals=true`:
- **G diagonal blocks: ALL POSITIVE (+0.46 to +0.53)** — should be ~-0.5 for interior CPs
- **G row sums: exactly 0** — confirms rank deficiency (CPs are EXTERIOR, not interior)
- **SVD: 2 exact zero singular values** — rank-2 null space
- **Cause**: CPs are on the EXTERIOR of the body

### Step 8: Normal direction analysis
The issue chain:
1. `generate_body` creates body with `flip_normals=true` → normals point **INWARD**
2. Initial `CPoffset = 1e-2` → CPs at `centroid + 0.01*L*inward_normal` → CPs are **INSIDE** (correct)
3. **BUG**: Solver does `b.CPoffset = -abs(b.CPoffset) = -0.01` for Dirichlet
4. Now CPs at `centroid - 0.01*L*inward_normal` → CPs are **OUTSIDE** (wrong!)
5. G_ii ≈ +0.5 (exterior), row sums ≈ 0 → rank deficiency → solution blow-up

**The solver's CPoffset sign logic assumes normals point outward, but `flip_normals=true` makes them point inward.**

### Step 9: Fix verification
Changed `flip_normals=true` → `flip_normals=false` in `generate_body` (both constructor calls):
- G row sums: ~1.0 (no rank deficiency)
- No zero singular values (smallest: 8.2e-4)
- Linear residual: 1.9e-13 (machine precision)
- Post-solve potential at solve CPs: rms ~1e-15
- **Flow Potential Residuals (off=-1e-10): [0.0017, 0.0026]** (was [0.6] / [1.15, 0.49])

## Root Cause
Two issues in `test/test_helpers.jl:generate_body`:

1. **`flip_normals=true` (PRIMARY)**: Flips normals to point inward. The solver's CPoffset logic (`-abs(CPoffset)` for Dirichlet) assumes outward normals, so CPs end up on the exterior instead of interior. This causes the G matrix to be rank-deficient and the solution to blow up for multi-body systems.

2. **`CPoffset=1e-14` (SECONDARY)**: Places CPs at the doublet singularity on the panel surface. Even with correct normal direction, this causes the self-influence evaluation to be numerically unstable, and the potential is only correct at the exact solve CP locations.

## Fix Applied
In `test/test_helpers.jl:generate_body`:
- Changed `CPoffset = 1e-14` → `CPoffset = 1e-2`
- Changed `flip_normals=true` → `flip_normals=false` (both constructor calls)

## Results After Fix
- Flow Potential Residuals: [0.0017, 0.0026] (previously [1.15, 0.49] or [0.6])
- Flow Tangency Residuals: [88.4, 3394] (still large — separate issue, likely related to mesh quality or wake modeling)
