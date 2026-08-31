# Item 015 Phase 4 — pinpointing the small-Das Dirichlet/Neumann divergence

**Written 2026-07-30 as a clear-context handoff.** Everything a fresh agent needs
is below with file:line references; no exploration should be necessary. Read
this file, then `BRAINSTORM/015_pressure_continuity_kutta_condition.md` (hub) and
the Phase 4 doc's last two sections only if you need the campaign narrative.

**Constraints, non-negotiable:**
- **Single-threaded** (`julia --project -t 1`). User directive. Also makes runs
  deterministic and dodges a pre-existing multithreaded FLOWVPM zero-particle
  bug in `runtests_unit_simulate.jl`.
- **No `src/` changes.** Driver/example edits only.
- **Do not touch the gate record.** `data/kutta_v1_attribution/ssw_sanity_gate.csv`,
  the 5% pass criterion, and `KUTTAV1_GATE_REFERENCE_BASIS` stay as they are.
  V1 stays blocked; reopening it is a separate decision.
- **T4 is GATED on explicit user approval — do not run it.** Run T1/T2/T3,
  report, and stop. See "T4 (GATED)" below.

---

## 1. Where this sits

Item 015 Phase 4 V1 is blocked at a Dirichlet-vs-Neumann formulation sanity gate
that fails at a 25.4% settled-CL gap (criterion ≤5%). A prior battery
(`KUTTAV1_STAGE=gate_diagnosis`, results in
`data/kutta_v1_attribution/gate_diagnosis_*.csv`) established:

- flat tip caps close only 3.9% of the gap → **eliminated**;
- the Dirichlet `G` kerneloffset is inert over 1e-8…1e-4 → **eliminated**;
- mesh refinement 480→1920 cells halves the unsteady gap and collapses the
  `:tri`/`:tri_robust` spread from 7.2pp to 0.35pp;
- the gap is monotone in η and approaches the semi-infinite anchor at large η.

That last point was reported as "η is the dominant lever". **The user correctly
objected that this framing is nearly tautological** — `Das` is the finite
stand-in for the semi-infinite row, so agreement at large `Das` is expected. The
real question is the opposite one.

## 2. The finding to explain

From `gate_diagnosis_summary.csv`, each body against **its own** semi-infinite
anchor, `|CL_KJ(η)| / |CL_KJ_steady|` (anchors: Neumann −0.399434,
Dirichlet-capped −0.387432):

| η | Neumann | Dirichlet |
|---:|---:|---:|
| 0.0625 | 0.934 | **0.767** |
| 0.125 | 0.962 | 0.823 |
| 0.25 *(gate)* | 0.975 | 0.876 |
| 0.5 | 0.982 | 0.920 |
| 1.0 | 0.986 | 0.950 |
| 2.0 | 0.989 | 0.970 |
| 4.0 | 0.991 | 0.982 |

**Both** bodies lose bound circulation as Das→0; the Dirichlet body loses ~3.5×
more (23.3% vs 6.6% at η=0.0625). This is a *differential sensitivity to the
near-TE wake*, and that asymmetry is what must be explained.

## 3. Established by code reading — do NOT re-derive

- **`Das` is freestream-parallel and frozen.**
  `_set_ssw_Das!(wing, config.eta*dt*full_uinf)`,
  `examples/suddenly_started_wing.jl:488` (setter at `:446-451`).
  `set_Das_eta_freestream=NaN` makes `initialize_Das!` return immediately
  (`src/FLOWPanel_simulate.jl:367-369`) and disables the per-step refresh
  (`:964-977`); `Das` is then only rotated with the body
  (`rotate_Das!`, `src/FLOWPanel_frames.jl:178-185`). **Not a suspect.**
- **Control points are identical for both formulations.** `calc_controlpoints!`
  (`src/FLOWPanel_abstractbody.jl:772-786`) puts the CP at the triangle centroid
  with **zero offset**, Neumann and Dirichlet alike. The "exterior for Neumann,
  interior for Dirichlet" line in the `Backslash` docstring
  (`src/FLOWPanel_solver.jl:208-211`) is **stale** — see
  `docs/exact_controlpoints.md`. Self-influence uses side-aware closed-form
  limits selected by `DBC` (`_self_limit`,
  `src/FLOWPanel_elements_fmm.jl:144-180`).
  *Consequence:* with identical geometry and identical wake strengths, the
  wake-induced velocity **at the control points is identical by construction**.
  The discrepancy cannot live in wake→CP evaluation; it must live in what each
  formulation does with that velocity.
- **The attached `Das` strip is a real influencing panel, in the LHS.** Two
  triangles per shedding panel, `(v1,v2,v1+Da)` and `(v1+Da,v2,v2+Db)` over
  adjacent TE **nodes** (`src/FLOWPanel_elements_fmm.jl:1043-1078`), emitted with
  the panel's own doublet strength (`:1036`). Upper/lower reversed vertex order
  (`shedding_full`, `src/FLOWPanel_liftingbody.jl:186-213`) makes the net strip
  circulation γ = μ_u − μ_l, so **G = G_body + W_attached·C** with `C` implicit.
  With `semiinfinite_wake=true` the same `_G!` path emits a semi-infinite sheet
  instead (`:1043-1044`), i.e. the whole wake sits in the matrix.
- **Wake geometry is NOT currently identical between the two runs.** Runs use
  `shed_with_induced_velocity=true`, `freestream_convection=false`, so every row
  convects with its own body's induced velocity (`propagate!`,
  `src/FLOWPanel_wake.jl:521-543`). The two wakes roll up differently from step 1.
  `SSWConfig.freestream_convection` (field at
  `examples/suddenly_started_wing.jl:115`) removes this: with it true, every row
  is translated rigidly by `dt*freestream` (`src/FLOWPanel_wake.jl:523-529`).
- **Strip/row-1 strength jump exists structurally.** The strip carries the
  **current** bound circulation (part of the LHS); wake row 1 carries the
  **previous** step's (`shed_wake!`, `src/FLOWPanel_simulate.jl:1204-1212`).
  Jump = γ_{k+1} − γ_k.

## 4. Leading hypothesis (traced to one line)

For the **Dirichlet** body the shed wake reaches the solve through exactly one
line — `set_strengths!`, `src/FLOWPanel_solver.jl:1083-1091`:

```julia
body.strength[:, 1] .= 0.0        # σ
for d in 1:3
    body.strength[:, 2] .= view(body.velocity, d, :)
    body.strength[:, 2] .*= view(body.normals, d, :)
    body.strength[:, 1] .-= body.strength[:, 2]
end
body.strength[:, 2] .= 0.0        # μ
```

i.e. `σ = −u·n`, where `u` = `body.velocity` already carries the shed-wake
induced velocity accumulated by `_sa_wake_influence!`
(`src/FLOWPanel_simulate.jl:468-470`). Then `rhs = −S·σ`
(`src/FLOWPanel_solver.jl:714`; `solve!` wrapper `:236-255`).
`formulation_prewake!`/`solve_formulation!` for the default
`VelocityThroughSources` are a no-op and a bare `solve!`
(`src/FLOWPanel_formulation.jl:788`, `:799-804`).

**So the shed wake is represented by an equivalent source sheet derived from its
normal velocity, and its direct contribution to the interior-potential condition
is dropped.** The exact RHS term is `−φ_wake`; the code substitutes
`−S·(u_wake·n)`.

The **Neumann** body has no analogue: `calc_bc_noflowthrough!`
(`src/FLOWPanel_solver.jl:703`, defn `:65-80`) consumes `u·n` directly — that
*is* its boundary condition, unconverted. `set_strengths!` for NK=1 (`:1099-1101`)
merely zeroes, and Neumann `solve!` (`:257-271`) never calls it.

As Das→0 the shed sheet approaches the TE, its potential contribution grows and
is progressively worse represented by a normal-velocity proxy — **only on the
Dirichlet side**. That is exactly the observed asymmetry and its Das-dependence.

**The repair already exists:** `GreenReconstruction`
(`src/FLOWPanel_formulation.jl:66`) transmits the wake as a reconstructed
potential trace `q` — solve `(I−B)q = Sσ`, then `G·μE = −S·σ0 − q`. It requires
Dirichlet (DBC=true, NK=2) + `Backslash` + `semiinfinite_wake=false`
(`_validate_formulation_common`, `:343-361`) — all satisfied by the gate cells.
Our A/jump pair **is** the legacy Kutta pair (`_is_legacy_kutta`,
`src/FLOWPanel_kutta.jl:436-438`), so `_validate_kutta_configuration` — which
would reject non-default formulations at `src/FLOWPanel_kutta.jl:489` — is
skipped entirely (`src/FLOWPanel_simulate.jl:878-880`). **`TraceCorrected` is
deprecated** and was rejected as a lever by DJI Gate A0 (~0.74%); do not use it.

## 5. Implementation

All edits in `examples/kutta_v1_attribution.jl` unless noted. Reuse the existing
machinery; do not write a parallel harness.

### 5.1 Existing machinery (current line numbers)

| symbol | line | role |
|---|---:|---|
| `_ssw_combo_run(config, combo; on_failure, grad_mu_options, settle_window, case_name)` | 230 | one unsteady run; calls `pnl.simulate!` |
| `run_ssw_combo(config, combo; optargs...)` | 277 | wrapper with pressure-failure fallback |
| `KUTTAV1_GATE_COMMON` | 313 | `(; AR=6.0, alpha_deg=5.0, n_span=12, n_airfoil=21, t_end_star=8.0, dt_star=0.125, backend_kind=:direct, save_vtk=false, verbose=false)` |
| `KUTTAV1_GATE_SETTLE` | 316 | `(7.0, 8.0)` |
| `KUTTAV1_GATE_BASES` | 327 | `(name=:tri, options=(; basis=:tri))`, `(name=:tri_robust, options=(; basis=:tri, tri_robust=true))` |
| `kj_lift_coefficient(config, body)` | 344 | `CL = 2∫Γdy/(U·S)`, Voronoi weights |
| `spanwise_agreement(ref, test)` | 376 | sign flips, Pearson, `scale_ratio` |
| `_run_gate_cell(cell, basis; battery, tag)` | 409 | **the cell runner — extend this** |
| `formulation_gap(n, d)` | 483 | the gap decomposition |
| `run_gate(output)` | 504 | the gate; **do not modify** |
| `KUTTAV1_DIAG_PAIR` | 627 | `neumann_uncapped` / `dirichlet_capped` specs |
| `_diag_basis(name)` | 640 | basis lookup |
| `_diag_push!(store, base, overrides, basis; battery, tag)` | 643 | run + file a cell |
| `run_gate_diagnosis(output; batteries)` | 662 | battery driver with env selection |

`_run_gate_cell` already has a `get_or(name, default)` override block reading
optional `eta`, `kerneloffset_over_c`, `n_span`, `n_airfoil` off the `cell`
NamedTuple, records them in `knobs`, and catches a throwing cell into a NaN row
with `status="failed:…"`. It also flags non-finite results as
`status="nonfinite"` (added because `kerneloffset=0` NaNs silently on the
Neumann body rather than throwing).

`SSWConfig` fields you need: `eta` (`:100`), `caps` (`:110`),
`freestream_convection` (`:115`), `kerneloffset_over_c` (`:130`).

### 5.2 Extend the cell spec

Add two optional `cell` fields to `_run_gate_cell`'s `get_or` block:

- `formulation` → default `pnl.VelocityThroughSources()`. Forward to
  `pnl.simulate!` via a new `formulation` kwarg on `_ssw_combo_run` (line 230)
  and `run_ssw_combo` (line 277). `simulate!` already accepts it
  (`src/FLOWPanel_simulate.jl:860`).
- `freestream_convection` → default `false`; it is an `SSWConfig` field, so it
  goes into the `overrides` NamedTuple alongside `eta` etc.

Record both in the row schema (`knobs`). For `formulation`, store a short label
(`"VelocityThroughSources"` / `"GreenReconstruction"`), not the struct — CSV.

**Guard, do not silently coerce:** `GreenReconstruction` is Dirichlet-only.
Assert `cell.bodytype === :dirichlet` when a non-default formulation is given,
and error otherwise with a clear message.

### 5.3 T1 — formulation axis (decisive for the observed CL gap)

η ∈ {0.0625, 0.125, 0.25, 0.5, 1.0} × three arms:

| arm | bodytype | caps | formulation |
|---|---|---|---|
| neumann | `:neumann` | `:none` | default (forced; Neumann rejects the others) |
| dirichlet_vts | `:dirichlet` | `:flat` | `VelocityThroughSources()` |
| dirichlet_green | `:dirichlet` | `:flat` | `GreenReconstruction()` (all defaults: `gauge=:area_mean`, `recompute_interval=1`, `green_solver=nothing`) |

`:tri` basis only — KJ is basis-independent (verified in the prior battery) and
this is about the solve. 15 runs.

**Primary metric:** circulation recovery `|CL_kj_unsteady| / |CL_kj_steady|` per
arm, per η. If `dirichlet_green` collapses onto the Neumann curve, the cause is
pinpointed. If it does not move, the hypothesis is dead and T2/T3 carry it.

### 5.4 T2 — geometry control

η ∈ {0.0625, 0.25, 1.0} × {neumann, dirichlet_vts} × `freestream_convection=true`.
6 runs. This makes the sheet straight and identical between formulations, and
analytically one uniform sheet — the regime in which the cancelling-interior-edge
argument documented at `examples/suddenly_started_wing.jl:757-766` says `Das`
should be nearly inert. Residual η-sensitivity here is **not** explainable by
wake shape, which removes the rollup confound from §3.

### 5.5 T3 — strip/row-1 strength jump

A monitor closure (monitors are called at `src/FLOWPanel_simulate.jl:1038-1047`,
**before** `propagate!`/`shed_wake!` at `:1098-1127` — the correct evaluation
point) recording per shedding edge per step:

```julia
pnl._get_wakestrength_Gamma(body, i, i_surf) - wake.strength[i_surf][1, 1, i]
```

(`_get_wakestrength_Gamma` / `_get_wakestrength_mu` at
`src/FLOWPanel_liftingbody.jl:1597-1615`; `PanelWake.strength` indexed
`[1, i_row, i_shed]`, `src/FLOWPanel_wake.jl:96-99`.) Report the settled jump
normalized by γ, per arm per η. Adds no runtime; fold into the T1/T2 runs by
appending the monitor to the tuple built in `prepare_suddenly_started_wing`
(`examples/suddenly_started_wing.jl:481`) — or, cleaner, accept an extra-monitor
argument in `_ssw_combo_run` so the SSW example is untouched.

### 5.6 Control assertions (cheap, once)

- `Das` ∥ U∞ to 1e-12 and unchanged between step 0 and the final step.
- Neumann and Dirichlet control points coincide to floating point on the shared
  mesh (guards the §3 zero-offset deduction).

### 5.7 New stage

Add `KUTTAV1_STAGE=gate_diagnosis2` → `run_gate_diagnosis2(output)` in
`main_kuttav1`'s dispatch and the file header. Write
`gate_diagnosis2_summary.csv`, `gate_diagnosis2_gaps.csv`,
`gate_diagnosis2_recovery.csv` (η × arm × recovery ratio),
`gate_diagnosis2_strip_jump.csv`. **No entry in `KUTTAV1_STAGE_GATES`** — this
stage gates nothing and is gated by nothing.

## 6. Run commands

```bash
cd /Users/ryan/Dropbox/research/projects/FLOWPanel.jl
KUTTAV1_STAGE=gate_diagnosis2 KUTTAV1_FORCE=true julia --project -t 1 \
    examples/kutta_v1_attribution.jl
```

Cost: ~11 min single-threaded (15 + 6 runs at 480 cells, ~30 s each; the
`GreenReconstruction` arm adds a one-time bordered LU per run).

## 7. Verification

1. **Harness self-checks, before interpreting anything:**
   - the η=0.25 `dirichlet_vts` and `neumann` cells reproduce the corresponding
     `gate_diagnosis_summary.csv` rows **bitwise** (`settled_CL`,
     `CL_kj_unsteady`, `steady_semiinfinite_CL`);
   - `CL_kj_steady` per body is unchanged across every new cell — the
     semi-infinite anchor cannot depend on η, formulation, or convection mode;
   - the two §5.6 control assertions pass.
2. Re-run `KUTTAV1_STAGE=gate` and confirm `ssw_sanity_gate.csv` and
   `gate_localization.csv` are **byte-identical** (md5 before/after).
3. `runtests_unit_kutta.jl` (658) + `runtests_unit_kutta_routeb.jl` (62) green,
   single-threaded, run from `test/`:
   ```bash
   cd test && julia --project=.. -t 1 -e 'include("runtests_unit_kutta.jl"); include("runtests_unit_kutta_routeb.jl")'
   ```
   Note `runtests_unit_simulate.jl` cannot be included standalone (needs
   `runtests.jl`'s preamble: `using Test, WriteVTK, LinearAlgebra`, `import
   FLOWPanel as pnl`, `include("test_helpers.jl")`), and errors at `-t 6` on a
   pre-existing FLOWVPM zero-particle bug unrelated to item 015.

## 8. Reporting

Append a "2026-07-30 — gate diagnosis II: small-Das attribution" section to
`BRAINSTORM/015_pressure_continuity_kutta_condition/phase_04_testing_verification_validation.md`
with the recovery curves, the T1 verdict, the T2 geometry control, and the T3
jump table. Answer explicitly the four questions the user raised, **including
the two §3 answers that needed no experiment** (Das direction/freezing; the
zero-offset control-point deduction).

**If T1 confirms**, amend the prior section's "η is the dominant lever" headline
and the `BRAINSTORM/INDEX.md` 015 Outcome cell: the Dirichlet arm's
η-sensitivity would be substantially an artifact of the default wake→body
transfer scheme rather than physics. That amendment is part of this work.

**If T1 disconfirms** (Green ≈ VTS), say so plainly, record the hypothesis as
eliminated, and hand the diagnosis to T2/T3 plus the gated T4.

Quote only numbers present in the new CSVs or simple arithmetic on their
published columns.

---

## 9. T4 (GATED — DO NOT RUN WITHOUT EXPLICIT USER APPROVAL)

User decision, 2026-07-30: *"gate T4 on user approval — I'm thinking it might
not be worth doing, but I want to see the results of the other tests before I
decide."* Run T1/T2/T3, report, and **stop**.

Recorded in full so it can be executed later without re-exploration. For context
when deciding: **T4a costs ~1 min and needs no solve at all**, and is the most
direct test of the §4 mechanism — T1 measures the mechanism's *effect on CL*,
T4a measures the *operator error itself*.

### T4a — local operator error (no solve, no truncation tail)

Both quantities refer to the same shed rows, so nothing is extrapolated:

- **Δ_code** = `−S·(u_wake·n)`. Get `u_wake` as the velocity delta across the
  wake influence pass — snapshot `body.velocity` before/after, the pattern in
  `formulation_prewake!`/`_split_sigma!` (`src/FLOWPanel_formulation.jl:760-790`)
  — then map σ→potential with `_source_potential!` (`:719-736`, save/restore
  wrapped).
- **Δ_exact** = `φ_wake`: `body.potential .= 0` then
  `influence!((body,), pnl._collect_wake_sources((wake,)), pnl.DirectBackend(); scalar_potential=true)`.

Report `‖Δ_code − Δ_exact‖∞ / ‖Δ_exact‖∞` and the L2 analogue vs η, for
`VelocityThroughSources` and `GreenReconstruction`. Prediction if §4 holds:
relative error is substantial and **grows as Das→0**, and shrinks sharply under
`GreenReconstruction`.

### T4b — global closure identity (the user's original formulation)

Same mesh, two bodies: `body_S` (`semiinfinite_wake=true`, whole wake in the
operator) and `body_U` (`semiinfinite_wake=false`, only the strip in the
operator). Prescribe the shed rows to the converged `γ_S` on straight geometry
at spacing `dt·U`, then check

    (G_U − G_S)·μ_S + Σ_k W_k·γ_S + W_tail·γ_S = 0

Implementation notes, all verified:

- **`Das` must be a UNIT vector for `body_S`** — `_phi_semiinfinite` errors if
  `|Da| ≠ 1` within `2eps()` (`src/FLOWPanel_elements_fmm.jl:1226-1229`). This is
  why `_ssw_steady_cl` sets `Das = full_uinf/‖full_uinf‖`
  (`examples/suddenly_started_wing.jl:577`). **You cannot sweep `Das` magnitude
  in a semi-infinite solve.**
- **Get `G` from your own `_G!` call**, not from `Backslash`: the constructor
  destroys `G` in place via `lu!` (`src/FLOWPanel_solver.jl:230`). Use
  `pnl._G!(Gcopy, body, body; kerneloffset=body.kerneloffset_panel, update_geometry=false)`.
  Optional finer split: `_assemble_B!` (`src/FLOWPanel_formulation.jl:469-480`,
  sets `body.suppress_attached_wake[]=true` → pure body operator) and
  `_assemble_W!` (`:486-517`), with `C` from `_shedding_edge_map`/
  `_apply_kutta_map!` (`:227-249`).
- **`W_tail` exactly**: `induced_semiinfinite` takes explicit vertices, no
  throwaway body needed (`src/FLOWPanel_elements_fmm.jl:1203`). Pass
  `v1 = TE1 + Da + L`, `v2 = TE2 + Db + L` and the unit direction. Vertex
  ordering and `Das` column indices must match `_induced_wake` (`:1012-1025`);
  `Das` is **vertex-based**, one column per TE node, indices in
  `shedding_full` rows 5/6 (`src/FLOWPanel_liftingbody.jl:181-211`).
- **Driving one solve with a prescribed wake:** use
  `pnl._steady_aerodynamics!(systems, (body,), (wake,), frames, uinf, body_solvers; backend_wake=..., backend_solve=..., backend_system=..., update_trailing_edges=false)`.
  **NOT `steady!`** — it hardcodes `wakes_tuple = Tuple(nothing for _ in systems_tuple)`
  (`src/FLOWPanel_simulate.jl:735`) and would silently apply **no wake at all**,
  producing a clean and meaningless result. With `update_trailing_edges=false`
  (the default, `:439`) nothing mutates `wake.nodes`/`wake.strength`:
  `update_TE!` is only called when `simulate!` passes `true` (`:982`), and
  `shed_wake!`/`propagate!` live in the stepping loop (`:1105-1127`), not in
  `_steady_aerodynamics!`.
- **Set `wake.nwakes[]`** or rows beyond it are invisible to `influence!`
  (`get_n_bodies(::PanelWake)`, `src/FLOWPanel_wake.jl:371`).
- **`DirectBackend` throughout**, and leave `wake_correction_active[] == false`.
  Two reasons: the buffer variant of `_induced_wake`
  (`src/FLOWPanel_elements_fmm.jl:1087-1157`) never applies `wake_strength_shift`,
  and `extra_farfield` (`src/FLOWPanel_fmm.jl:72-78`) is an FMM-only accuracy
  accommodation for semi-infinite sources.
- **Residual read-out**: no dedicated BC-residual routine exists. Cleanest is
  `body.potential .= 0; influence!(body, body, backend; scalar_potential=true)`
  then compare against `Sσ` with `_solve!`'s sign convention; see
  `_collect_coupled_operator!` (`src/FLOWPanel_solver.jl:536-547`) for the
  Dirichlet (`body.potential`) vs Neumann (`sum(velocity .* normals; dims=1)`)
  read-out. A `KrylovSolver(body; backend=DirectBackend())` is also callable as
  a bare matvec (`:368-401`) but clobbers `body.strength`/`velocity`/`potential`.

**Sign conventions are the one real implementation risk.** Anchor them by
requiring the identity to close for a configuration where it must — the Neumann
body — **before** reading any Dirichlet number. Expected outcome if §4 holds:
Neumann residual at round-off, Dirichlet residual far above it and growing as
Das→0.

### T5 — dropped, recorded so it is not re-derived

Strip aspect ratio vs wake-row spacing, via matched pairs of constant η·n_span.
The strip's cross dimension is the spanwise TE **node** spacing Δy = b/n_span
(`src/FLOWPanel_elements_fmm.jl:1043-1078`), so its aspect ratio Δy/|Das| is
governed by `n_span` — **not** by the chordwise panel size, which never enters
the strip geometry. At η=0.0625, n_span=12 the strip is a sliver of AR ≈ 64.
Dropped because it needs a mesh change, and `n_span` simultaneously sets the
spanwise resolution of the bound circulation and the KJ integral, so a
non-collapse could not be cleanly attributed. T4 tests the same suspicion on a
single mesh. This is the fallback if T4 comes back clean.
