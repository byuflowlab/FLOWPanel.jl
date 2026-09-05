# Static surface-vorticity conversion evidence

Simulation-free `DirectBackend` comparison for BRAINSTORM item 016 Phase 3
Stage 4b. The prescribed open sheet uses the fixed smooth maps

```julia
x(s) = 3expm1(1.1s)/expm1(1.1)
y(q) = 1.5sinh(0.8q)/sinh(0.8)
z(s,q) = 0.10sin(pi*s)cos(pi*q/2)
```

and a smooth non-affine strength field. The representative grid has 6
streamwise rows and 12 spanwise panels, `sigma=0.06`, and `overlap=1.3`.
The pure state is unchanged; the hybrid is formed by one conversion followed
by analytical removal of the outgoing source row (no propagation, shedding,
or aerodynamic solve). Probes cover both sheet sides at
`d/sigma = 0.25, 0.5, 1, 2, 4, 8`; `probe_errors.csv` is the compact handoff
family excerpt. The unit test repeats handoff, interior, root, and tip families.

> **Caveat added 2026-08-03 (Stage 4c finding R4).** `generate.jl` never sets
> `PanelWake.overflowed[]`, and
> `FastMultipole.get_n_bodies(::FilamentWrapper{<:PanelWake})`
> (`src/FLOWPanel_wake.jl:2641`) returns **zero** unless it is set. Every field in
> `probe_errors.csv` was therefore computed with the retained final filament
> emitting nothing — for both the pure and the hybrid state, and for all three
> strategies. Since the retained filament is the one `attribution`-dependent
> source, that table cannot have measured the half-jump-filament effect it was
> read as measuring. `ledgers.csv` and `refinement.csv` are unaffected (they read
> particle strengths and per-panel `kappa`, not fields). Use
> `generate_attribution.jl` and its CSVs instead; that harness sets `overflowed`
> and verifies the filament is live before drawing any conclusion.

`ledgers.csv` reports the exact deposited-particle vector ledgers
`sum(Gamma)` and `sum(x cross Gamma)`. For smooth conversion, the independently
assembled outgoing-filament circulation is included and agrees with the
deposited circulation to roundoff. `refinement.csv` records the final three
levels of the fixed smooth stretched-grid map; the centered `:split` rule is
second order, while `:upstream` is first order.

Decision: retain `:upstream`. Although `:split` has the expected refinement
order and is more accurate for the representative coarse gradient, it does not
lower handoff RMS at every standoff and its maximum error regresses at multiple
standoffs. The conjunctive default-change rule therefore fails before the
order/coarse advantages can change the default.

Reproduce the direct probe and ledger rows from the repository root with
`JULIA_NUM_THREADS=6 julia --project data/surface_vorticity_conversion_static/generate.jl`.

## Stage 4c — discriminating attribution criterion (2026-08-03)

`generate_attribution.jl` remediates review finding R3: Stage 4b did not select
`:upstream`, it failed to displace it. Every rule, threshold, and disposition was
committed to
`BRAINSTORM/016_surface_vorticity_particle_shedding/phase_03_implementation.md`
("Stage 4c — pre-registration") **before** this harness was first run.

Same prescribed sheet, maps, and strength field as `generate.jl`, with three
corrections: `overflowed[]` is set so the retained filament is a live source
(R4); the conversion runs on its steady-state branch with the incoming filament
already carrying `alpha`; and the reference is a refined sheet, never the
edge-jump sheet being replaced.

| File | Contents |
| --- | --- |
| `attribution_collapse.csv` / `_order.csv` | T0: alpha spread vs row extent `h_row`, and the fitted order `p` in `S ~ h_row^p` |
| `attribution_reference_convergence.csv` | T1 gate: the refined reference must be converged well below the alpha signal |
| `attribution_accuracy.csv` | T1: paired per-probe error of each mode against the Richardson-extrapolated reference |
| `attribution_hazard.csv` / `_verdict.csv` | T2: panel-induced stretching at deposited particles, by class, over 18 replicate conditions |
| `attribution_kappa_alpha.csv` | Per-panel error-minimizing, consistency, and second-order alpha, with kappa error at 0 / 0.5 / 1 |

**Headline results.**

- **T0 — the alpha spread does not collapse in the near field.** Fitted order
  rises monotonically with standoff: `p = -0.46` at `d/sigma = 0.25` through
  `+0.11` at 8 to `+0.91` at 64. The alpha difference is a dipole layer of
  thickness `h_row`; inside it (`d < h_row`, which covers every standoff Stage 4b
  probed) it does not decay, and outside it decays at roughly first order.
  **`alpha` is a near-field modelling choice and a far-field truncation
  artifact** — so its spread must be carried in the CT error budget, since the
  body sits in the near field.
- **T1 — conclusive, and regime-dependent with a clean crossover.** Reference
  order 1.99, Richardson residual `1.13e-4` against a paired signal of
  `7.40e-3` (65x margin). `:upstream` has lower error at **all 12 probes** at
  every standoff `d/sigma <= 2`; `:split` has lower error at **all 12** at every
  standoff `d/sigma >= 4`. The aggregate (56.9% favouring `:split`) hides a
  perfectly systematic structure — exactly the aggregate-only reporting failure
  Stage 4b was criticised for.
- **T2 — no veto.** Paired p95 stretching differences (`:split` minus
  `:upstream`) are positive but unresolved at the pre-registered 2-sigma
  threshold: interior `+2.38 +/- 1.94`, root/tip `+0.72 +/- 0.52`, perimeter
  `+4.02 +/- 3.58`.
- **alpha\* — `:upstream` is the best admissible weight for pointwise kappa.**
  The error-minimizing alpha is **1.45** (outside `[0,1]`), the consistency
  weight is 1.20, and the second-order weight on this graded mesh is 0.23 — not
  0.5. Relative kappa error is 0.468 / 0.309 / 0.150 at alpha = 0 / 0.5 / 1, so
  `:upstream` is 2x more accurate than `:split` here. Stage 4b's "second order"
  holds only in the refinement limit, where grading tends to 1 and both weights
  tend to 0.5.

**Verdict: `:upstream` retained, now on evidence.** The near field is the
regime the item designated in advance (Phase 1 sec. 5.7; Ryan's 2026-08-01
Stage 4b charge names "the near field this item exists to clean up"), and there
`:upstream` wins unanimously and is corroborated independently by the kappa
result. T4 (the pre-registered Phase 4 CT comparison) remains the mission-level
confirmation.
