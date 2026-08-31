# Item 015 Phase 4 — remediate V1 review findings and unblock the wing gate

## Context

The 2026-07-30 three-agent clear-context review of the V1 attempt found: (1) the
Dirichlet-vs-Neumann sanity-gate diagnosis is unreproducible (only 1 of ~10
numbers has a script/CSV behind it) and internally contradictory (Γ agrees ~3%
yet KJ CL differs 2.5×; Neumann's own KJ vs pressure CL disagree 3.3×), so the
doc's localization overclaims; (2) the scripted gate samples an unsettled
transient (t*∈[1.1,2]) and records no strength-level comparison, so it cannot
separate solve from velocity-reconstruction effects (`basis=:tri` on an all-tri
mesh is the repo's known-divergent configuration; prior sweptwing mirror work
traced a CL gap to exactly this); (3) `examples/sweptwing.jl`'s default body was
silently flipped to Dirichlet with no toggle and no formulation token in output
paths; (4) diamond-table blemishes (mixed Δp̂ normalizations across cells,
unnormalized "fixture CL", unflagged Route-B 24.5%-of-scale correction);
(5) doc/log gaps (hub log stale, "3–4%" is ~3.0%, V0-vs-V1 simulate-suite
contradiction unexplained, ssw0 resume can mask a failed gate).

**Decisive new fact (user directive):** the Dirichlet gate bodies are OPEN-TIP —
`suddenly_started_wing.jl:183` ("Triangular open-tip extrusion"), built
`watertight=false` (:300); `simplewing` (`examples/helper_functions.jl:203-247`)
lofts sections with no caps either. The interior-Dirichlet Green identity
assumes a closed surface (the solver warns for the Neumann analogue but is
silent for non-watertight Dirichlet, `src/FLOWPanel_solver.jl:261-265`). The
DJI campaign's convention was capped/Dirichlet vs uncapped/Neumann, and the
user reports **flat caps match Neumann no-caps best**. So the leading suspect
for the 14–26% gap is the missing caps, ahead of `:tri` reconstruction. The fix
plan makes flat caps the Dirichlet configuration and rebuilds the gate as a
scripted, discriminating experiment.

## Workstream A — flat tip caps for Dirichlet wing bodies (the physics fix)

1. **Generic utility** `add_flat_tip_caps(nodes, cells)` in
   `examples/helper_functions.jl` (usable by both SSW and sweptwing): detect
   the two open boundary loops (reuse the open-edge machinery behind
   `pnl.iswatertight` — it already identifies non-manifold/open edges), then
   close each loop with a centroid-fan triangulation (add one centroid node per
   tip, fan to all loop nodes). Centroid-fan preserves the SSW mesh's z→−z
   in-plane symmetry and y→−y tip-to-tip symmetry automatically when the
   contour nodes are symmetric (they are, mesh revision 2). Return updated
   (nodes, cells).
2. **SSW**: add `SSWConfig` field `caps::Symbol = :none` with `:flat`
   supported; **default `:flat` whenever `bodytype=:dirichlet`** (constructor
   promotion), `:none` for `:neumann` (legacy unchanged). Apply caps in the
   mesh step BEFORE body construction; keep shedding derived from the
   constructed probe body (invariant). Pass `watertight =
   pnl.iswatertight(...)` (should now be true for capped). Include `caps` in
   `_ssw_case_tag`. Extend `assert_ssw_mesh_symmetry` to accept capped meshes
   (cap faces must satisfy the same face-set reflection tests).
3. **sweptwing**: the Dirichlet variant gets the same utility applied to the
   `simplewing` output (only in the kutta/V1 path — see Workstream C).
4. **Unit checks** (in the driver or a small testset): capped SSW mesh is
   watertight; symmetry assertion passes; shedding pairing count unchanged vs
   uncapped; Neumann uncapped results bitwise unchanged vs pre-edit.

## Workstream B — rebuild the sanity gate as a scripted, discriminating experiment

Rewrite the gate stage of `examples/kutta_v1_attribution.jl` (`run_ssw_sanity`)
so every number the diagnosis needs is produced by the script and lands in CSV:

- **Settled comparison**: `t_end_star=8.0` (settle window mean over t*∈[7,8]),
  not 2.0. Keep dt/mesh as before (AR=6, α=5°, n_span=12, n_airfoil=21).
- **Cells** (each a CSV row): Neumann-uncapped (legacy oracle),
  Dirichlet-flat-capped (new configuration), Dirichlet-uncapped (retained once,
  to quantify how much of the old 14–26% gap the caps close — this is the
  direct test of the flat-caps-match-Neumann observation). Each × grad_mu
  basis ∈ {`:tri`, `:tri_robust`} (quad basis N/A on an all-tri mesh) to
  separate reconstruction sensitivity from the solve.
- **Recorded per cell**: settled CL, full CL history file, semi-infinite steady
  CL (legacy default kutta pair — validator not involved), spanwise TE
  Δstrength (Δγ or Δμ) distribution CSV (the driver already computes these in
  `run_ssw_combo`; stop discarding them in the gate path), and a **scripted KJ
  diagnostic**. Validate the KJ extractor first on the Neumann steady case
  against its own pressure CL and the lifting-line anchor
  CL_LL = 2πα·AR/(AR+2) ≈ 0.41; if it can't reproduce Neumann pressure CL to
  ~10%, report the extractor broken and exclude KJ from gate evidence
  (resolving the review's 3.3× internal contradiction one way or the other).
- **Pre-registered gate criterion**: settled-CL relative gap
  |Dir(capped) − Neu|/Neu ≤ 5% AND spanwise TE Δstrength distributions agree
  qualitatively (no sign flips, shape match). Hard-stop semantics unchanged.
- Runs are local, ≤6 threads; ~6 unsteady t*=8 runs at the gate mesh — minutes
  each; if any exceeds ~20 min, stop and route to HPC per `agent_policies/HPC.md`.

## Workstream C — sweptwing default restoration

- Restore `examples/sweptwing.jl:76` default to the legacy
  `RigidWakeBody{pnl.VortexRing, 1, Float64, false}`; add an opt-in switch
  (env `SWEPTWING_BODYTYPE=neumann|dirichlet`) selecting the Dirichlet
  source+doublet type (+ flat caps via Workstream A when dirichlet).
- Add a formulation token to `run_name`/`save_path`/`grid_tag` (:66-72) and to
  the mirrored-pressure cache load paths, so Neumann-era artifacts (incl. the
  tracked `sweptwing_Cps_mirrored.png` inputs) can never mix with Dirichlet
  outputs.
- The V1 sweptwing stage in the driver drives it with `dirichlet` + the
  finite-wake settings the kutta validator needs (`semiinfinite_wake=false`,
  explicit finite Das) — configured from the driver, leaving the classic
  example's standalone behavior at legacy defaults.

## Workstream D — driver metric fixes

In `examples/kutta_v1_attribution.jl`:
- **Unify diamond Δp̂ normalization**: one fixed scale (freestream q = ½ρ‖U∞‖²)
  for the pre AND post columns of every cell; record the kutta frozen
  `pressure_scale` as its own column. Drop the `_fill_pre_dp` cross-run
  planting or label the column honestly (`dp_hat_pre_from_jump_run`).
- **Iteration means**: exclude `:startup_jump` records from
  `mean_outer_iterations`/`mean_body_solves`; report startup steps separately.
- **Resume masking**: on `ssw0` resume, re-read `ssw0_gate.csv` and re-throw if
  `pass=false` (a checkpoint must not present a failed stage as complete);
  parse gate CSVs by header name, not last-column position.
- Remove dead `KUTTAV1_RHO`.
- Rerun the diamond stage after the normalization fix (cheap) so
  `diamond_summary.csv` is single-metric.

## Workstream E — documentation remediation

- `BRAINSTORM/015_pressure_continuity_kutta_condition/phase_04_testing_verification_validation.md`:
  mark the 2026-07-29 diagnosis numbers (t*=8 row, semi-infinite row, KJ CLs,
  Γ ranges) as unscripted and superseded; soften the localization sentence to
  "diagnostics mutually inconsistent; localization not established"; fix
  "3–4%" → ~3.0%; gloss "fixture CL" as unnormalized (area-dimensional) lift;
  state (and then note as fixed) the Δp̂ normalization mix; flag B/pressure's
  24.5%-of-scale correction and 4.2× worse jump baseline as a tracked
  relocated-sensitivity signal for V2/V3; reconcile the V0-green vs V1-not-green
  `runtests_unit_simulate.jl` statements (threading difference); add the new
  gate design (caps rationale: flat caps match Neumann no-caps best; open-tip
  Dirichlet violates the closed-surface Green identity) and all commands.
- Hub `BRAINSTORM/015_pressure_continuity_kutta_condition.md`: add progress-log
  entries (V1 attempted/blocked; review; remediation), update "Current phase".
- `BRAINSTORM/INDEX.md` 015 Outcome cell: update after the gate rerun (concise).

## Workstream F — rerun and conditional continuation

1. Rerun diamond (fixed metrics) and the new gate battery.
2. **Gate passes** (capped Dirichlet within 5%): proceed with the original V1
   sequence unchanged — `ssw0` zero-lift matrix → `ssw_alpha` four-cell matrix
   → sweptwing steady A-cells — using the capped Dirichlet configuration
   everywhere, then document the attribution tables in phase_04 + INDEX.
3. **Gate still fails**: document the capped result (incl. the uncapped
   comparison row and basis discriminator), keep V1 blocked, and record the
   surviving suspects (Dirichlet G kerneloffset regularization at
   `src/FLOWPanel_solver.jl:228` vs the DJI unregularized-by-design convention;
   `:tri` reconstruction) as the next gate-diagnosis targets. A negative result
   is a valid outcome.

## Out of scope (riders, do not fix here)

- SSW example test `n_span=1` breakage vs the even-`n_span` guard — item 014.
- Silent solver acceptance of non-watertight Dirichlet bodies (a warning
  symmetric to the Neumann one at `src/FLOWPanel_solver.jl:261-265` would be
  nice) — note it in the doc, propose as follow-up, don't change src/ in V1.
- `TraceCorrected` deprecation follow-ups.

## Verification

- `test/runtests_unit_kutta.jl` (658) + `runtests_unit_kutta_routeb.jl` (62)
  green after all edits (no src/ changes expected; driver/example edits only,
  plus `helper_functions.jl`).
- Cap utility checks (Workstream A.4), including Neumann-path bitwise
  no-change.
- Gate battery CSVs regenerate every number quoted in the updated doc
  (reproducibility is itself a review criterion next time).
- Keep ≤6 local threads throughout.
