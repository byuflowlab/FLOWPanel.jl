# Linear-Strength Panels in FLOWPanel

> **Status:** Scaffolding only. Section bodies are placeholders; derivations and citations to be filled in iteratively.

## 1. Motivation and Scope

- Why move beyond constant-strength panels: improved accuracy, smoother surface fields (pressure, velocity), fewer panels for equivalent accuracy, better behavior near edges and stagnation regions.
- Hard requirements driving this design:
  1. **FMM-accelerated** — must compose with the `FastMultipole` backend in `src/FLOWPanel_fmm.jl` and `src/FLOWPanel_elements_fmm.jl`.
  2. **Fully determined system** — no auxiliary constraints (e.g. `Σσ = 0`) needed to close the linear system.
  3. **Mathematically rigorous** — derivations must be explicit and verifiable, not heuristic.
- Out of scope (for now): higher-order (quadratic+) panels, curved/parametric panels, panel deformation across time steps.

## 2. Literature Review

*TODO — populate with full citations.*

- Classical references: Hess & Smith linear-source formulations; Morino's linear-doublet integral equation; Maskew (VSAERO); Katz & Plotkin Ch. 10–11.
- Node-based vs. panel-based parameterizations (Vassberg; Johnson).
- Closed-body indeterminacy in linear-doublet formulations, and equivalence of a C0-continuous linear-doublet sheet to a constant vortex sheet on each panel.
- FMM with linearly varying source/dipole densities: Greengard & Rokhlin extensions; high-order density handling in KIFMM / bbFMM / black-box FMM.
- Prior implementations to study: PANAIR (constant-source / linear-doublet hybrid), VSAERO, PMARC.

## 3. Parameterization

### 3.1 Node-valued strengths
- DOFs live at mesh vertices: `σ_v` (linear source), `μ_v` (linear doublet).
- Within each triangle, strength interpolates linearly via barycentric basis functions `N_1, N_2, N_3`.

### 3.2 Why this is fully determined
- DOF count = `n_nodes`; boundary-condition count = `n_nodes` (if collocated at nodes) or `n_cells` (if collocated at panels).
- Argument for why no global constraint like `Σσ = 0` is required — to be filled in once collocation choice (§6) is fixed.

### 3.3 Closed-body identity (linear doublet ↔ constant vortex sheet)
- Claim: on a closed surface with a C0-continuous linear-doublet distribution, the induced field equals that of a constant-strength vortex sheet on each panel, with sheet strength `γ = n̂ × ∇_s μ`.
- Proof deferred to §4.

### 3.4 Control-point placement
- Options: node-collocated; panel-centroid; dual-mesh centroids.
- Trade-offs: conditioning of the influence matrix; FMM near-field singularity handling (nodes lie *on* the surface); compatibility with existing `CPoffset * characteristiclength` machinery.

## 4. Kernel Derivations

*TODO — placeholder.*

### 4.1 Linear source panel
- Potential and velocity integrals over a flat triangle with `σ(ξ,η) = σ_1 N_1 + σ_2 N_2 + σ_3 N_3`.
- Reduction to constant-strength integrals plus first-moment integrals `∫ ξ/r dA`, `∫ η/r dA`.
- Edge-line-integral form (extension of Hess–Smith) for numerical robustness near the panel.

### 4.2 Linear doublet panel
- Reduction to an equivalent vortex-line representation along the panel edges via Stokes' theorem applied to `∇μ × n̂`.

### 4.3 Closed-surface reduction
- Sum of edge contributions across shared edges: interior edges cancel between adjacent panels with consistent vertex values, leaving only boundary / Kutta edges.
- Connects directly to the existing `VortexRing` element machinery.

### 4.4 Self-influence and near-field
- Closed-form or quadrature-based treatment of the singular `1/r` and `1/r²` kernels on the source triangle itself.

## 5. Multipole Coefficients for FMM

*TODO — placeholder.*

### 5.1 Source multipole expansion
- Multipole moments of a linearly varying `σ` about an expansion center `x_c`.
- Decompose as `(constant-σ moment) + (first-moment correction)`; both reduce to closed-form polynomial integrals over the triangle.

### 5.2 Doublet multipole expansion
- Linear `μ` ⇒ constant vortex strength on each edge (from §4.2/§4.3) ⇒ reuse `FastMultipole`'s existing vortex-segment expansion path.

### 5.3 Local expansion / target side
- Unchanged: targets see the same harmonic potential regardless of source-side parameterization.

### 5.4 Required interface additions
- New methods/trait specializations in `src/FLOWPanel_elements_fmm.jl`.
- Possible extensions to `ProbeSystem` if vertex-collocated control points are used.

## 6. Solver Formulation

*TODO — placeholder.*

### 6.1 Unknown vector
- Length `n_nodes` rather than `n_cells`.

### 6.2 Collocation choices
- **Node-collocated Neumann:** enforce `n̂·U = 0` at each node, using an area-weighted normal averaged from incident panels.
- **Panel-collocated Neumann:** keep current control-point machinery; system becomes rectangular (`n_cells × n_nodes`) — either least-squares or dual-mesh remapping needed. Discuss square-system conditions.
- **Dirichlet (interior `φ = 0`) variant.**

### 6.3 Square-system / Kutta closure
- DOF-balance argument: how Kutta conditions on lifting bodies absorb residual DOFs without ad-hoc constraints.

### 6.4 Compatibility with existing solvers (`src/FLOWPanel_solver.jl`)
- `BackslashNeumann` / `BackslashDirichlet`: matrix shape and assembly loop change.
- `KrylovSolver`: matrix-free path; needs node↔panel scatter/gather around `influence!`.
- `FGSSolver`: node-wise Gauss-Seidel; revisit relaxation parameter `rlx`.

## 7. Code-Change Roadmap

### 7.1 New element types
- `LinearSource`, `LinearDoublet` mirroring the constant variants in `src/FLOWPanel_elements.jl`.

### 7.2 Body storage changes
- `strength` is currently `n_cells × N`. For linear element types it must be node-indexed.
- Introduce a trait, e.g. `panel_dof_location(::Type{E}) = :cell | :node`, and either a separate `node_strength` field on `AbstractBody` or a generalized `strength` storage keyed by DOF location.

### 7.3 FMM integration
- New `_Uind!` / `_phi!` methods in `src/FLOWPanel_elements_fmm.jl` for the new element types.
- Extend interaction-list handling if needed for vertex-collocated targets.

### 7.4 Assembly / solve
- Update `_solve_*` paths in `src/FLOWPanel_solver.jl` and `src/FLOWPanel.jl` to scatter node strengths into panel evaluations and gather panel residuals into node equations.

### 7.5 Post-processing
- `src/FLOWPanel_postprocess.jl` currently assumes per-cell strengths. Pressure / force integration needs either a node→cell projection or a direct node-based form.

### 7.6 Wake coupling
- `src/FLOWPanel_wake.jl`: rigid-wake shedding edges must match the node-doublet representation along the trailing edge.
- Decide whether `RigidWakeBody` needs a separate "linear wake" type or can keep constant-strength wake panels.

### 7.7 Tests
- Add `test/runtests_unit_linearpanels.jl`.
- Analytical cases: sphere (source-only), ellipsoid, possibly Karman–Trefftz wing.
- Convergence-order study at multiple refinements (expect ~2nd order vs. constant panels' ~1st).

## 8. Open Questions

- Best control-point strategy under FMM when nodes lie on the surface (singularity / near-field).
- Interaction of node-based unknowns with the `ImplicitAD` differentiation path used by `solve_ludiv!`.
- Whether the rigid wake should also be promoted to linear doublets or remain constant-strength.

## 9. References

*TODO — bibliography to be populated.*
