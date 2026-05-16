# Dirichlet Rotor Hover Verification

## Goal

Verify the Dirichlet rotor-hover solve sequence:

1. Compute freestream plus wake-induced velocity on the body.
2. Choose source panel strengths.
3. Compute source-panel induced potential on interior-offset control points.
4. Solve doublet/vortex-ring strengths so the interior perturbation potential is zero.
5. Compute body-on-all influence on the exterior/output control points.

## Initial Code-Path Verification

- `simulate!` resets body and wake fields, applies freestream velocity, applies kinematics, snaps the wake trailing edge, then applies wake-source influence to the body before calling `solve!`.
- Dirichlet `solve!` calls `set_strengths!`, zeroes `body.potential`, computes body-on-body scalar potential, then calls `_solve!`.
- `set_strengths!` sets the fixed source column to `sigma = -U dot n` and clears the solved second column before the linear solve.
- `_G!` uses interior control points for `DBC=true` and stores scalar potential influence in the matrix.
- `_solve!` solves into strength column 2 for `Union{ConstantSource,VortexRing}` Dirichlet bodies.

## Initial Issue Found

The single-body Dirichlet `Backslash` path used `body.potential` as workspace for the interior source-induced potential and did not restore the pre-solve external potential afterward. This means step 5 could accumulate body-on-all potential onto stale interior solve workspace after control points were moved back outside.

The coupled solver already saves and restores external potential, so the single-body direct solver should do the same.

## Implemented Fix

- `solve!(body, solver)` for Dirichlet bodies now saves the incoming external potential before using `body.potential` as solve workspace.
- The solve still uses `body.potential` as the interior source-potential workspace for RHS construction.
- The `finally` block restores both `CPoffset` and the saved external potential, so post-solve body-on-all influence starts from the correct pre-solve external state.
- Unit tests now assert that single-body Dirichlet Backslash solves restore `body.potential` for both lifting and nonlifting source+doublet/source+vortex-style paths.

An initial implementation used the `Backslash.phi_ext` field directly, but the same Dirichlet `solve!` method also handles `KrylovSolver`. The final implementation uses a local saved potential vector so all Dirichlet solver types keep the same state-restoration behavior.

## Unit Verification

Command:

```julia
julia --project=. -e 'include("test/test_helpers.jl"); include("test/runtests_unit_solver.jl")'
```

Result:

- `Solvers`: 78 passed, 0 failed.
- The Dirichlet Backslash tests now check that `body.potential` is restored after solve.
- Existing Krylov Dirichlet tests also pass, confirming the generic restore logic does not require a Backslash-specific field.

## Rotor First-Step Verification

The rotor diagnostic reproduced the first simulation step without VTK output and checked each requested stage.

### Step 1: Freestream plus wake-induced velocity

- Dirichlet formulation confirmed: `DBC=true`.
- Interior control point offset confirmed: sample `(controlpoint - centroid) dot normal = -3.2457400531497063e-10`.
- After freestream: velocity norm extrema `(0.0001, 0.0001)`.
- After kinematics: velocity norm extrema `(5.004722211721764, 74.78155674072701)`.
- After wake influence: velocity norm extrema unchanged for the first empty-wake step, `(5.004722211721764, 74.78155674072701)`, and potential extrema `(0.0, 0.0)`.

### Step 2: Source panel strengths

- Expected source strengths were computed independently as `sigma = -U dot n`.
- Maximum absolute error between stored source strengths and expected values: `0.0`.
- Source strength extrema: `(-39.82088956717775, 68.3493361903501)`.

### Step 3: Interior source-panel potential

- The Dirichlet solve uses the interior-offset control points and computes scalar potential before RHS construction.
- The solve-time potential workspace is no longer allowed to leak out of `solve!`.

### Step 4: Doublet/vortex-ring strength solve

- Solved strength column is finite.
- Solved strength extrema: `(-0.21763883947306595, 0.14776901868042996)`.
- Recomputed interior perturbation potential residual extrema: `(-2.6987267346781844e-9, 2.7168715992968725e-10)`.
- Maximum absolute residual: `2.6987267346781844e-9`.
- External velocity restoration after solve: `0.0`.
- External potential restoration after solve: `0.0`.

### Step 5: Body-on-all influence

- Exterior control point offset confirmed after solve: sample `(controlpoint - centroid) dot normal = 3.2457400531497063e-10`.
- Body-on-all influence then produced body velocity norm extrema `(0.9549090547019256, 90.19568634974647)`.
- Body potential extrema after body-on-all influence: `(-0.14776899009904343, 0.21763883057948485)`.

## Current Assessment

The requested Dirichlet sequence is now verified for the first rotor step. The main defect found during verification was stale interior solve potential leaking into the post-solve body-on-all influence path; the fix restores external potential before step 5.
