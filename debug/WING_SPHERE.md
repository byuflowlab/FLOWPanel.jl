# Wing+Sphere Example Status

## What was accomplished

Created `examples/wing_sphere.jl` — a multi-body unsteady simulation with:

1. **Watertight NACA 0012 wing** (AR=2, chord=1m, span=2m) as `RigidWakeBody{Union{ConstantSource, VortexRing}}`
   - Built programmatically with structured mesh: 30 chordwise × 25 spanwise stations
   - Rounded end caps via `sin²(t)` taper scaling over 30% of semispan
   - Cap holes closed with center-node fan triangles at each tip (2 center nodes + 60 fan triangles)
   - 752 nodes, 1500 panels, 16 shedding edges along main wing trailing edge
   - Shedding matrix constructed manually with correct nia/nib local indices

2. **Sphere** (R=0.5m) as `NonLiftingBody{ConstantSource}`, centered at [2.5, 0, 0] (downstream)
   - Standard parametric grid approach from sphere.jl example
   - 840 nodes, 1600 panels

3. **Simulation infrastructure**
   - `PanelParticleWake` for wing (nwakerows=3, max_particles=20000)
   - Tuple-based `simulate!` with `(wing, sphere)` and `(wake, nothing)`
   - Static `ReferenceFrame` with `dependent_index=[1,2]`
   - 21 time steps, dt = chord/Vinf/5

4. **Flow tangency verification**
   - Sphere: **PASS** — RMS(U·n)/Vinf = 0.0%, max|U·n| ≈ 0.0
   - Wing: **FAIL** — RMS(U·n)/Vinf ≈ 6.0%, max|U·n| ≈ 10.75

## What remains to be fixed

### Primary issue: Wing flow tangency (~6% RMS)

The wing has poor flow tangency with max|U·n| = 10.75 (at Vinf=30), concentrated at specific panels. Root cause is NOT:
- FMM approximation (tested with DirectBackend — same result)
- Normal direction (verified all 1500 normals point outward)
- Solver choice (Backslash for AbstractLiftingBody already dispatches to BackslashDirichlet)

Likely causes to investigate:
1. **Shedding edge correctness** — The manually-constructed shedding matrix may have incorrect local node indices. The nia/nib assignments need verification against the actual cell vertex ordering. A single wrong index would cause the wake strength computation to be wrong, leading to poor BC enforcement near the TE.

2. **Cap panel quality** — The heavily-tapered cap panels (scale down to 0.02) create very elongated triangles that may cause conditioning issues in the influence matrix. The fan triangles at the tip center connect to a single node, creating many triangles sharing a vertex.

3. **Neighbor matrix** — The raw `(nodes, cells, shedding)` constructor defaults to `neighbor = zeros(Int, 3, ncells)` (no neighbor info). The `compute_mu_gradient!` function in the Dirichlet solver uses neighbor info for the doublet gradient computation. Missing neighbors would cause incorrect surface velocity recovery.

### Suggested next steps

1. **Verify shedding**: Print out the physical nodes at each shedding edge and verify they form a continuous chain along the TE. Check that upper/lower panel assignments are correct.

2. **Test with steady-state solve first**: Before the unsteady simulation, do a steady solve (`solve!` with semi-infinite wake) and check tangency with `calcfield_U!` + `apply_freestream!`. This eliminates wake dynamics as a variable.

3. **Populate neighbor matrix**: Build the neighbor connectivity from the cell data so the doublet gradient computation works correctly.

4. **Consider using `calc_shedding`**: If the mesh can be wrapped in a `GridTriangleSurface`, use the existing `calc_shedding` function instead of manual construction.

### Other notes

- Backend is currently `DirectBackend()` (slow, ~120s for 21 steps). Switch back to `FastMultipoleBackend` once tangency is resolved.
- The sphere works perfectly — zero tangency error confirms the simulation infrastructure and multi-body coupling are correct.
