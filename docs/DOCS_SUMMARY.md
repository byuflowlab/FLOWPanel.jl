# FLOWPanel.jl Documentation — Current State Summary

## Build System

- **Tool:** [Documenter.jl](https://documenter.juliadocs.org/), configured in `docs/make.jl`
- **Output:** HTML site deployed to GitHub Pages at `flow.byu.edu/FLOWPanel.jl/dev/`
- **Auto-generation:** `make.jl` calls two pre-build scripts before building:
  - `src/generate_index.jl` — generates `src/index.md` from the project `README.md`
  - `src/generate_examples.jl` — generates all files under `src/examples/` from Julia source files in `FLOWPanel/examples/`
- **DO NOT manually edit** `src/index.md` or anything in `src/examples/` — they are overwritten on every build
- Sidebar disabled; collapse level 1; favicon at `assets/favicon.ico`
- API Reference section exists but is commented out in `make.jl` (lines 84–87)

---

## Page Structure (42 pages)

### Intro (`src/index.md`)
*Auto-generated from project README via `generate_index.jl`.*

- FLOWPanel logo and tagline
- Capabilities list (structured mesh generation, unstructured import, multi-element, flexible solver, fast computation)
- Current limitations: explicitly states "No fast multipole method" (now outdated — FMM has been added)
- Tested on Julia v1.10
- Installation instructions: `] add FLOWPanel`
- Link to selected publication (Alvarez et al. 2023 AIAA AVIATION)
- About/contributors/license
- Sponsors banner

---

### Potential Flow (`src/potentialflow.md`)
Theory reference (327 lines). Covers:

- Helmholtz decomposition of velocity fields
- Laplace equation for incompressible irrotational flow
- Green's second identity derivation
- Boundary integral equation (BIE) formulation
- Source and doublet boundary integral expressions
- Internal Dirichlet and Neumann boundary condition formulations
- No code — purely mathematical derivations with LaTeX equations

---

### Elements

#### Panel Definition (`src/elements/paneldefinition.md`)
- Panel local coordinate system (ξ, η, ζ axes)
- Basis vectors derived from panel nodes
- Diagram of panel geometry and normal orientation
- Note: GeometricTools defines normals with the right-hand rule, opposite to Katz & Plotkin convention

#### Constant Source (`src/elements/constantsource.md`, 131 lines)
- Formulation of constant-strength source panel (planar polygon)
- Velocity kernel (analytic near-field, far-field cutoff)
- Potential kernel
- Notes on singularity handling: kernel offset `δ`, far-field threshold `ε`

#### Constant Doublet (`src/elements/constantdoublet.md`, 139 lines)
- Equivalence of constant doublet panel and vortex ring
- Velocity kernel via Biot-Savart (vortex ring interpretation)
- Potential kernel
- Same `δ`/`ε` parameters as source

#### Semi-Infinite Doublet (`src/elements/semiinfdoublet.md`, 196 lines)
- Horseshoe vortex formulation (two semi-infinite trailing legs + bound vortex)
- Velocity kernel (Biot-Savart for finite and semi-infinite segments)
- Used for lifting bodies with trailing wake in steady flow

#### Non-Planar Semi-Infinite Doublet (`src/elements/semiinfnonplanardoublet.md`)
- Extension of horseshoe vortex to non-planar trailing edge geometry
- Handles swept/dihedral trailing edges

#### Constant Vortex Sheet (`src/elements/constantvortexsheet.md`, 387 lines)
- Vortex sheet strength distribution over a planar panel
- Biot-Savart integration (analytic result with trigonometric/logarithmic terms)
- Both velocity and potential kernels
- Optimization notes: avoids repeated trig evaluations

---

### Geometry Engine

#### Grid Generation Overview (`src/geometry/gridgeneration.md`)
- High-level description of the three mesh generation strategies: lofting, body of revolution, transformations
- Points to GeometricTools.jl as the underlying engine

#### Lofting (`src/geometry/gridgeneration-loft.md`, 531 lines — largest doc)
- `gt.surface_loft` API: lofting a surface through a series of cross-section curves
- Detailed index notation for node/cell ordering in the resulting structured grid
- Step-by-step examples: simple wing, multi-section wing
- Chordwise and spanwise panel distribution parameters
- Diagrams of coordinate system indexing

#### Path Lofting (`src/geometry/gridgeneration-pathloft.md`)
- Extruding a cross-section along a 3D path
- Useful for non-planar geometries (bent tubes, curved surfaces)

#### Body of Revolution (`src/geometry/gridgeneration-rev.md`, 230 lines)
- `gt.surface_revolution` API: rotating a profile curve around an axis
- Parameters: profile coordinates, angular discretization, angle range
- Example: centerbody/hub geometry
- Output: structured azimuthal grid

#### Space Transformations (`src/geometry/gridgeneration-transf.md`, 132 lines)
- `gt.transform!`, `gt.lintransform!` for rotation, translation, scaling
- Composing multiple transformations
- Use case: positioning sub-geometries within an assembly

#### Triangulation (`src/geometry/gridgeneration-triang.md`, 277 lines)
- `gt.GridTriangleSurface(grid, dimsplit)`: converting quadrilateral grids to triangles
- `dimsplit` parameter controls which diagonal is used
- Impact on panel normal orientation and solver quality

#### GeometricTools Basics (`src/geometry/basics.md`)
- Overview of the GeometricTools.jl package and its role in FLOWPanel
- Grid types: `gt.Grid`, `gt.GridTriangleSurface`

#### Grid Operations (`src/geometry/basics-grid.md`)
- Reading/writing grids
- Node and cell access patterns
- Subsetting and splitting grids

#### Transformations (`src/geometry/basics-transformations.md`)
- Detailed transformation API reference

#### Looped Grid (`src/geometry/basics-loopedgrid.md`)
- Closed (periodic) structured grids for axisymmetric bodies
- Azimuthal continuity handling

#### Surface Grid (`src/geometry/basics-surfacegrid.md`, 237 lines)
- `gt.GridTriangleSurface` in depth: accessing nodes, cells, normals, areas
- Field storage on surface grids

#### Panel Gradient (`src/geometry/panel-gradient.md`)
- Computing field gradients over panels
- Finite-difference and analytic gradient options

---

### Tutorials (23 pages — all auto-generated from `examples/*.jl`)

All example pages follow the same pattern:
1. Visualization image
2. Brief description
3. Full Julia code block (truncated at a `break_flag` comment)
4. Result plots
5. CL/CD/Cm comparison tables
6. Tip block showing how to run via `include(joinpath(pnl.examples_path, ...))`

#### Swept Wing (3 pages)

**4.2° Angle of Attack** (`examples/sweptwing-4p2aoa.md`, 224 lines)
- 45° swept-back wing, rigid wake, vortex ring elements
- Geometry via `gt.surface_loft` with NACA 0012 sections
- Solver: `solve_backslash!` (old API — backslash solve)
- Post-processing: chordwise pressure distribution, spanwise loading, CL/CD/Cm
- Comparison to experimental data

**AOA Sweep** (`examples/sweptwing-aoasweep.md`)
- Sweeps AOA from -5° to 15°
- Plots CL, CD, Cm vs AOA; validates against experiment

**Solver Benchmark** (`examples/sweptwing-solver.md`, 242 lines)
- Compares: backslash (`\`), LU decomposition (`solve_ludiv!`), GMRES (`solve_gmres!`)
- Docs reference `FLOWPanel.solve_backslash!`, `FLOWPanel.solve_ludiv!`, `FLOWPanel.solve_gmres!`
- Timing and accuracy comparison
- Note: this uses the old solver API; new API uses solver type objects (`BackslashNeumann`, etc.)

#### Centerbody (3 pages)

**Source Elements** (`examples/centerbody-source.md`)
- Axisymmetric hub/centerbody using constant source panels
- Geometry via `gt.surface_revolution`
- Validates surface pressure coefficient vs analytic solution

**Slice** (`examples/centerbody-slice.md`)
- Post-processing: extracting and plotting field slices through the body
- Uses `calcfield_U!` and related post-processing functions

**Vortex Ring** (`examples/centerbody-vortexring.md`)
- Same centerbody but with vortex ring elements instead of source panels
- Comparison of source vs vortex ring formulations

#### Duct (3 pages)

**AOA Sweep** (`examples/duct-aoasweep.md`, 226 lines)
- Fan duct geometry, swept over angles of attack
- Demonstrates non-lifting body with doublet elements
- Validates against experimental data

**Fluid Domain** (`examples/duct-fluiddomain.md`)
- Generating a volumetric field grid around the duct
- Computing velocity field at off-body points
- Visualization in ParaView

**Least Squares** (`examples/duct-leastsquares.md`, 204 lines)
- Demonstrates least-squares solver for watertight mesh with redundant equations
- Uses `solve_leastsquares!` (old API)

#### Blended Wing Body (5 pages)

**CAD Model** (`examples/blendedwingbody-cad.md`)
- Importing a CAD geometry via SolidWorks → `.stl` → Gmsh → `.msh`
- Mesh processing with Gmsh Python API

**Unstructured Meshing** (`examples/blendedwingbody-gmsh.md`)
- Full Gmsh workflow: importing geometry, setting mesh size, generating surface mesh
- Exporting `.msh` file for use in FLOWPanel

**Trailing Edge** (`examples/blendedwingbody-TE.md`)
- Defining shedding edges on unstructured meshes
- Wake shedding geometry for lifting body analysis

**Aerodynamic Analysis** (`examples/blendedwingbody-aero.md`, 272 lines)
- Importing the Gmsh mesh via `FLOWPanel.read_gmsh` / `Meshes.jl`
- Full solve and post-processing: CL, CD, Cm, surface pressure
- Validation results

**GPU vs CPU** (`examples/blendedwingbody-gpucpu.md`, 236 lines)
- Optional CUDA backend for GPU acceleration
- Performance comparison: GPU vs CPU, timing breakdown
- Notes on memory requirements

#### Cessna 210 (4 pages)

**OpenVSP Geometry** (`examples/cessna-openvsp.md`)
- Exporting geometry from OpenVSP as DegenGeom/CompGeom files
- Converting to FLOWPanel body

**Trailing Edge** (`examples/cessna-TE.md`)
- Defining trailing edges from OpenVSP geometry for wake shedding
- Handling multi-surface configurations (wing + stabilizers)

**Aerodynamic Analysis** (`examples/cessna-aero.md`, 277 lines)
- Full Cessna 210 (wing + horizontal stabilizer + fuselage) analysis
- Multi-body solve
- CL, CD, Cm comparison to published data

**VSP Geometry Processing** (`examples/cessna-vspgeom.md`)
- Detailed walkthrough of OpenVSP geometry file format
- Node/cell extraction and coordinate frame handling

#### Lifting Line (4 pages)

**Weber Wing 4.2° AOA** (`examples/ll-weber-4p2aoa.md`, 245 lines)
- Coupling FLOWPanel with a lifting-line wake model
- Weber wing geometry (experimental validation case)
- `liftingline` module usage

**Weber Wing AOA Sweep** (`examples/ll-weber-aoasweep.md`)
- AOA sweep using the coupled panel + lifting-line model
- CL, CD vs AOA compared to experiment

**Stability Derivatives** (`examples/ll-stabilityderivatives.md`, 266 lines)
- Computing dCL/dα, dCm/dα using finite differences over the panel+LL system
- Demonstrates differentiation-friendly solve interface

**A50K27 Wing** (`examples/ll-a50k27.md`, 208 lines)
- A50K27 airfoil, validation against experimental data
- Demonstrates custom airfoil section in lifting-line model

---

## Resources Directory

- `data/` — airfoil coordinate files (RAE 101, NACA series, Selig), CL/CD comparison data (`.md` tables), sphere experimental pressure data
- `data/images` symlink → `../src/assets/images/`

---

## Assets

- `src/assets/images/` — 50+ PNG/SVG diagrams: panel geometry sketches, grid coordinate diagrams, result plots (CL/CD curves, pressure distributions), ParaView screenshots
- `src/assets/favicon.ico`

---

## Known Issues / Staleness

- `index.md` states "No fast multipole method" — this is now false (FMM backend added on `fastmultipole` branch)
- Solver benchmark (`sweptwing-solver.md`) documents `solve_backslash!`, `solve_ludiv!`, `solve_gmres!` — these are the old function-based API; new API uses solver type objects
- No documentation for new solver types: `BackslashNeumann`, `BackslashDirichlet`, `KrylovSolver`, `FGSSolver`
- No documentation for `FastMultipoleBackend` or `DirectBackend`
- No Quick Start page
- No dedicated Solvers reference page
- No documentation for matrix-based mesh input (nodes + cells arrays)
- API Reference section is commented out in `make.jl`
