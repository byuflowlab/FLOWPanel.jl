# Phase 3 — Implementation

**Status:** APPROVED 2026-08-03 (explicit user approval, following the
2026-08-03 clear-context re-review recommendation).

**Prerequisite:** approved Phase 2 architecture — satisfied 2026-08-01

**Phase approval:** [x]

Implement the architecture **as amended by Phase 2 §13**. In particular the
ghost is the pre-shift final active row read in place (Option B): `PanelWake`
storage, `capacity` derivation, and every source/probe/filament/VTK view are
unchanged, and only the new row 1 is staged, only when `nwakerows == 1`.
Conversion metadata belongs in `src/FLOWPanel_replay.jl`, and conversion
decoding must not reuse the `NoShed()` unknown-tag fallback.

## Staged execution

| Stage | Content | Status |
| --- | --- | --- |
| 0 | Gate bookkeeping; legacy golden-reference characterization test | **DONE** |
| 1 | Strategy axis and dispatch; no behavior change | **DONE** |
| 2 | Reconstruction core (tangent basis, rank-aware gradient, bilinear quadrature) | **DONE** |
| 3 | Transaction, deposition, root/tip closure, handoff | **DONE (divergence form; R1/R2 guards verified)** |
| 4b | Static sheet-vs-particle equivalence evidence (no timestepping) | **DONE** |
| 4 | Persistence: metadata, replay, warm start | **DONE** |
| 5 | Diagnostics and close-out | **DONE** |
| 4c | Discriminating attribution criterion (R3 remediation; pre-registered) | **DONE** — `:upstream` retained ON EVIDENCE; R4 raised against Stage 4b |

Stage 4b is user-requested acceptance evidence: with no time marching, convert
a prescribed sheet once and compare the deposited particles against the sheet
they replace — induced velocity on an off-body probe shell versus standoff,
far-field circulation and impulse moments, legacy-versus-smooth near-field
difference on the same sheet, and internal-edge cancellation across two
successive conversions.

**Added 2026-08-01 (Ryan):** Stage 4b must also settle the streamwise
`attribution` default by measuring `:split` against `:upstream` on the same
sheet — (a) near-field probe error, i.e. whether the half-jump filament
`:split` retains at the handoff degrades the near field this item exists to
clean up, and (b) resolution sensitivity, i.e. whether `:split`'s second-order
advantage actually materializes on a realistically stretched wake. Both rules
conserve circulation exactly, so the decision rests entirely on this evidence;
the current `:upstream` default is provisional.

## Stage 4c — pre-registration (written 2026-08-03, BEFORE the harness was run)

Stage 4c remediates R3. Stage 4b did not select `:upstream`; it failed to
displace it. This section fixes the knobs, metrics, thresholds, and the
disposition for **every** outcome including the null, and is committed before any
Stage 4c number exists. R3 arose because a decision rule appeared after the data;
a rule written afterwards can only ever return the incumbent. If any part of this
section is changed once results are in, the change and its reason must be logged
explicitly and the affected conclusion demoted.

**Instrument correction that shapes the design.** Stage 5's
`diagnose_nearfield` (`_surface_vorticity_nearfield`) is invoked at *preflight*
(`src/FLOWPanel_wake.jl:2509`) — before commit and before the row shift, with the
outgoing row still present and the retained filament still in its
`attribution`-independent legacy form. It is faithful to Phase 2 §8.3 as written
but is **structurally blind to `alpha`**, so it cannot be the hazard instrument.
T2 builds its own post-commit probe. Conversely
`data/surface_vorticity_conversion_static/generate.jl` already emulates the
post-shift state correctly (`pw.nwakes[] = N-1` after conversion makes
`_final_filament_strength` return the alpha-weighted retained filament); that
emulation is reused verbatim, because it is the state fresh particles live in.

### T0 — is `alpha` a truncation knob or a modelling choice? (runs first)

The centered-difference argument (`alpha = 1/2`) and the artificial-partition
argument (`alpha = 1`) conflict *only* because the outgoing row has finite
streamwise extent. As `h_row -> 0` the two faces coincide and the spread must
vanish. Which criterion is even appropriate depends on whether it does.

Fixed physical sheet and fixed `sigma`; sweep `nrows` over {4, 6, 8, 12, 16}
(`prescribed_sheet` uses `s = (r-1)/nrows` and assigns strength from centroid
*position*, so the sheet and the physical `mu` field are invariant under the
sweep). Probe positions are fixed in **physical space** across the sweep, never
tied to node indices. For each `nrows` x `alpha` in {0, 1/2, 1}: convert, emulate
post-shift, evaluate the hybrid field at the fixed probes. Report the spread
`S(h_row) = max over alpha pairs of ||u(alpha) - u(alpha')||`, raw and
normalized, and fit the observed order `p` in `S ~ h_row^p`. Report `h_row` as a
physical length alongside `sigma/h_row`, the real similarity parameter (rotor:
`sigma ~ 1.5` local chords, `h_row = dt*U`), which is what connects this to
018's dt ladder.

| Outcome | Disposition |
| --- | --- |
| `p` in [0.8, 1.2] and `S` extrapolated to production `h_row` below the Phase-4 CT tolerance | `alpha` is a discretization knob. Select on accuracy (T1), subject to T2's veto. Carry `S` in the error budget. |
| `p ~ 0` (no collapse) | `alpha` is a genuine model parameter. Carry the measured spread as model-form uncertainty on every reported CT, exactly as item 014 does for `Das`. **This is a result, not a deferral.** |
| Intermediate `p`, or `S` above tolerance at production `h_row` | Report as measured; T1/T2 still run, but the selection is provisional and the spread is carried. |

### T1 — accuracy against a refined reference (the selector)

The Stage 4b metric compared the hybrid against `u_sheet`, the representation
being replaced; it therefore scores fidelity to the edge-jump artefact this item
exists to remove. Replace it with a converged object.

- **Reference:** the same physical sheet at `nrows_fine`, all rows present. A
  constant-`mu` panel's field is exact, so the panel sheet is a
  piecewise-constant approximation to the continuous doublet distribution and
  converges in `nrows`. Convergence is **demonstrated**, at
  `nrows_fine` in {32, 64, 128}, and the residual reference error must be well
  below the alpha differences it is asked to resolve. If it is not, T1 is
  declared inconclusive and T3 fires.
- **Test object:** hybrid at production-like `nrows_coarse` — retained coarse
  rows plus alpha-converted particles, post-shift emulation.
- **Metric:** `E(alpha)` per probe, reported raw **and** normalized separately,
  so the numerator/denominator confound of the Stage 4b metric is visible rather
  than baked in.
- **Statistic:** the *paired* difference `E(:split) - E(:upstream)` per probe.
  The coarse-panel error is common-mode and cancels in the pairing; RMS-of-ratios
  was the wrong statistic.
- **Probes:** families {handoff, interior, root, tip}, standoffs `d/sigma` in
  {0.25, 0.5, 1, 2, 4, 8}, both sides, and at least 8 spanwise stations (Stage 4b
  used 2, giving M = 4 with no scatter estimate). Report median, p95, and max of
  the paired difference, not only its mean.

### T2 — the hazard, judged on the hazard (the veto)

Phase 1 §5.7's concern is a concentrated line coincident with fresh particles,
feeding stretching and relaxation — not sheet similarity.

Post-commit probe: after conversion and post-shift emulation, evaluate the
panel-induced velocity gradient `J` at each *deposited particle position*,
sourcing only `get_sources(pw)` in the post-shift state (retained rows plus the
alpha-weighted retained filament). Metric: stretching magnitude
`||J * Gamma_p|| / ||Gamma_p||` per particle, reported separately for interior,
root/tip, and perimeter classes — never pooled, per §8.3.

**Replicate set, clarified 2026-08-03 before any Stage 4c run** (the draft said
"8 span stations", which is a T1 *probe* concept and cannot apply to a metric
evaluated at every deposited particle): the replicates are `sigma` in
{0.5, 1, 2}*`sigma_0` x `nrows` in {6, 12} x `nspan` in {8, 12, 16} = 18
conditions. This resolves an ambiguous phrase; it does not touch the threshold
or the decision rule.

**Rule (lexicographic, not conjunctive; user selection 2026-08-03 — "a resolved
increase disqualifies"):**

> `:split` is **rejected** iff the paired difference in p95 stretching
> (`:split` minus `:upstream`) is positive **and** its magnitude exceeds **twice
> the standard deviation of that paired difference across the replicate
> conditions**. Otherwise `:split` is **adopted**, and T1's accuracy result is
> uncontested.

T2 is the veto; T1 is the selector. A null in T2 is a decision, not a stalemate.

### T3 — reference-free fallback (conditional)

Runs **only if** T1's reference fails its convergence demonstration, or the alpha
differences fall inside the reference error. Freeze the wake geometry, run K
consecutive conversions, and measure where the alpha fields differ. Since each
face is deposited exactly once regardless of alpha, the difference must be
confined to a one-row boundary layer at the partition; report its amplitude and
spatial extent (distance from the partition at which `|u(alpha1) - u(alpha2)|`
falls below a stated fraction of the local `|u|`). No truth object required.

### T4 — Phase 4 null, pre-registered here and recorded there

Recorded in `phase_04_testing_verification_validation.md`; not run (Phase 4 is
blocked).

### `alpha*` — diagnostic only (user selection 2026-08-03)

**Definition corrected 2026-08-03, before any Stage 4c run.** The plan proposed
`alpha* = h_{N-1}/(h_{N-1}+h_N)`. Deriving it properly shows that is not the
second-order weight. With `a = (h_{N-1}+h_N)/2` and `b = (h_N+h_{N+1})/2` the
upstream/downstream centroid separations, Taylor expansion of

```
kappa_streamwise ~ [ alpha*(muhat_A - muhat_G) + (1-alpha)*(muhat_G - muhat_D) ] / h_N
```

gives leading term `-[alpha*a + (1-alpha)*b] mu' / h_N` and second term
`[alpha*a^2 - (1-alpha)*b^2] mu'' / (2 h_N)`. So:

- cancelling `mu''` requires `alpha = b^2/(a^2+b^2)`;
- *consistency* (leading term equal to `-mu'`) requires `alpha*a + (1-alpha)*b = h_N`.

Both hold simultaneously only when `a = b = h_N`, i.e. a uniform mesh — where
every `alpha` is consistent and `alpha = 1/2` is second order, exactly matching
Stage 4b's measured 1.98/1.99 (a smooth map's grading ratio tends to 1 under
refinement, so the inconsistency vanishes in the limit). **On a fixed graded
mesh no `alpha` restores second order**: the `h_N` denominator is implicated as
well as the numerator weight. That is a substantive theory result and belongs in
the record regardless of what the harness measures.

Stage 4c therefore reports, per panel: the **error-minimizing** `alpha` (closed
form — `kappa(alpha)` is affine in `alpha`, so `||kappa(alpha) - kappa_exact||^2`
is a quadratic with an exact minimizer), the theoretical `b^2/(a^2+b^2)`, the
local spacing ratio, and the kappa error at `alpha` in {0, 1/2, 1, optimal}. The
affinity itself is verified by checking that the measured `:split` kappa equals
the midpoint of the `:downstream` and `:upstream` kappas.

**No fourth `attribution` mode is added to the source**, so the Phase 3 diff, its
validation, and its metadata/fingerprint surface are untouched. A material win is
the case for implementing one in a later phase, not a reason to widen this one.

## Prior-phase handoff

Phase 1 supplies the approved theory only after its gate is passed. Phase 2
must then provide the exact typed strategy API, state ownership,
perimeter/handoff ledger, transactional conversion protocol, diagnostics,
failure behavior, metadata/replay/warm-start representation, and file/test
map.

Read the hub and the approved
[Phase 2 architecture](phase_02_architecture.md) before beginning. Do not
infer unresolved architecture directly from the theory.

## Scope

Implement only the approved architecture around the existing
`PanelParticleWake` conversion seam. The likely source touchpoints are
`PanelWake` storage/source indexing, `PanelParticleWake`,
`_convert_to_particles!`, the current conversion-before-row-shift call in
`shed_wake!`, metadata/replay reconstruction, and warm-start state.

The existing edge-jump strategy remains the exact constructor default and
must retain current `method_trailing`/`method_unsteady` behavior. The
surface-vorticity strategy adds rank-aware gradient reconstruction,
physical-area quadrature, true-perimeter closure, handoff bookkeeping,
capacity preflight, and diagnostics. Its conversion uses a staged post-shift
view: the old final active row becomes a non-source ghost only after the new
active row needed by the stencil is constructed.

## Deliverables

- Approved concretely typed conversion strategies and constructor option.
- Explicit active-row capacity, always-allocated ghost geometry, and source
  views that ignore the ghost.
- Gradient, bilinear-quadrature, perimeter, and handoff implementation.
- A preflight/commit conversion transaction: a partially deposited outgoing
  row is never observable.
- Versioned metadata plus backward-compatible replay and warm-start behavior.
- Deterministic validation, diagnostics, and failure paths.
- Focused unit tests for the new strategy and byte-for-behavior legacy
  regression where practical.
- A progress log with exact verification commands and results.

## Required focused tests

- Constructor default and explicit legacy selection are equivalent.
- Existing legacy locations, strengths, counts, line policies, metadata,
  replay, and warm starts remain unchanged.
- Constant and affine planar reconstructions, rotated sign covariance,
  resolution floor, warped/nonuniform geometry, and a full-rank one-sided
  `nwakerows == 1` ghost stencil.
- Active source counts for `nwakerows=1,2,3`; the ghost never contributes
  panel influence.
- No gradient dependence on `Das`, body strengths, virtual TE panels, or an
  outside-zero value.
- Internal-edge cancellation, true perimeter retention, and single-owned
  panel/particle handoff.
- Root/tip versus interior panel-induced velocity-gradient diagnostics.
- Insufficient capacity, invalid geometry, and invalid configuration leave
  both the particle field and panel-row state unchanged.
- New metadata round-trips and old manifests select the legacy default.

## Acceptance gate and progress log

Phase 3 is complete only when the scoped diff matches the approved
architecture, focused tests pass, legacy behavior remains unchanged, and the
user explicitly approves the implementation. Phase 4 remains blocked until
that approval.

- **2026-08-03 — CLEAR-CONTEXT RE-REVIEW (post-remediation): approval
  RECOMMENDED; user checkbox still required.** Independently verified, fresh
  context: **(a) R1/R2 guards** exist at `_validate_conversion_filaments`
  (`src/FLOWPanel_wake.jl:866-877`), called from the `PanelParticleWake`
  constructor before `PanelWake`/particle-field allocation (:1791), legacy
  no-op method untouched; replay rejects `unsteady_filament=false` smooth
  metadata (`runtests_unit_replay.jl` invalid-manifest case, `@test_throws
  ArgumentError`). **(b) The external transfer ledger is genuinely
  source-derived now** — `S(w)` returns zero unless a `FilamentWrapper` is
  present in `get_sources(pw)` before consulting `_final_filament_strength`
  (`runtests_unit_wake.jl` "transfer matches the wake's own filament ledger").
  **(c) Stage 4c pre-registration compliance spot-checked in the harness**:
  `overflowed[] = true` is set (R4 fixed for its own measurements), the
  reference is Richardson from nrows 256/512, statistics are paired with
  median/p95/max, and the `alpha`-affinity check is present. **(d) Every
  Stage 4c number cross-checked against the checked-in CSVs holds**: T0 orders
  (−0.457 at `d/sigma`=0.25), T2 paired p95 `+2.383±1.940` / `+0.720±0.519` /
  `+4.021±3.577` with `split_rejected=false` at 18 replicates each, `alpha*`
  interior ≈1.46–1.47 with consistency weight 1.20 and second-order weight
  0.229, kappa errors 0.468/0.309/0.150, reference rungs 8.66e-3 → 1.68e-3 →
  4.23e-4 (order 1.99). **(e) The near-field designation is genuinely
  pre-dated** — Ryan's 2026-08-01 Stage 4b charge in this document names "the
  near field this item exists to clean up", and Phase 1 §5.7 records the
  handoff/root-tip hazard; both precede any Stage 4c number, so retaining
  `:upstream` on the near-field branch is not post-hoc. **(f) Suites re-run
  fresh at `JULIA_NUM_THREADS=6`**: wake 583/583, replay 125/125, warm start
  12/12 + 14/14, simulate 81/81, zero `Fail`/`Error`, exit 0.

  *One precision note, non-blocking (existing entries left untouched for
  record integrity):* the T1 "12 of 12 probes" unanimity refers to the
  **handoff (`aft`) family** — 12 of the 32 probes per standoff — which is the
  regime the charge designates. Root/tip agree in sign (0/2 at ≤1, 2/2 at ≥4;
  the crossover standoff `d/sigma=2` is mixed at 1/2), while `:control` at
  mid-chord favours `:split` at *all* standoffs including ≤2 (16/16), which is
  what locates the effect at the partition. Pooled across all families the
  near field is 16/32, so the unanimity claim must be read family-scoped, as
  the `:control` sentence already implies. No conclusion changes.

- **2026-08-03 — STAGE 4c COMPLETE: `:upstream` retained ON EVIDENCE; R3
  closed; new finding R4 raised against the Stage 4b evidence.** Harness
  `data/surface_vorticity_conversion_static/generate_attribution.jl`,
  simulation-free, `DirectBackend`, no source change. All rules were committed
  above before the first run; the two amendments made (T2's replicate set,
  the `alpha*` definition) were both made **before** any Stage 4c number existed
  and are logged in place.

  **R4 (new, against existing checked-in evidence) — Stage 4b's probe table was
  computed with the retained filament emitting nothing.**
  `FastMultipole.get_n_bodies(::FilamentWrapper{<:PanelWake})`
  (`src/FLOWPanel_wake.jl:2641`) returns **zero** unless `overflowed[]` is set,
  and only `shed_wake!` ever sets it. `generate.jl` builds its fixture
  statically and never does, so every field in `probe_errors.csv` — pure and
  hybrid, all three strategies — omitted the retained final filament. That
  filament is the *only* `attribution`-dependent source in the comparison, so
  the table cannot have measured the half-jump-filament effect it was read as
  measuring. `ledgers.csv`/`refinement.csv` are unaffected (particle strengths
  and per-panel `kappa`, not fields). Caveat added to the data README.
  This was caught by the pre-registered instrument self-check, which required
  the T2 probe to be demonstrably `alpha`-sensitive and measured exactly
  `0.000e+00` until `overflowed[]` was set. **Related trap worth carrying:** a
  ledger that asks only whether a `FilamentWrapper` is *present* repeats R2 in a
  subtler form — the wrapper can be present and emit zero bodies. `Stage 4c`'s
  ledger requires `get_n_bodies > 0`.

  **T0 — `alpha` does not collapse in the near field.** Fitted order in
  `S ~ h_row^p` rises monotonically with standoff: `-0.457` at `d/sigma = 0.25`,
  `-0.191` at 4, `+0.113` at 8, `+0.512` at 16, `+0.912` at 64. The `alpha`
  difference is a dipole layer of thickness `h_row`; inside it (`d < h_row`,
  which covers **every** standoff Stage 4b probed) it does not decay, and
  outside it it decays at roughly first order. Per the pre-registered
  disposition this is the third row (report as measured, carry the spread) in
  the near field and the first row in the far field. Operative statement:
  **`alpha` is a near-field modelling choice and a far-field truncation
  artifact**, so its spread belongs in the CT error budget — the body sits in
  the near field.

  **T1 — CONCLUSIVE, and regime-dependent with a clean crossover.** The
  reference gate passed only after Richardson extrapolation from `nrows`
  256/512 (observed order 1.99, residual `1.132e-4` against a paired signal of
  `7.396e-3`, a 65x margin). The first run used `gold = u(nrows=128)`, whose own
  discretization error (`8.7e-3` relative change over the preceding rung)
  exceeded the paired signal, so the gate correctly declared it **inconclusive**
  and its sign carried no weight.

  **Correction, entered when this was checked rather than asserted.** An earlier
  draft of this entry said the unconverged reference "gave the opposite sign"
  and that the gate caught a result that would have been published backwards.
  That is wrong. Two things changed between the runs — the reference *and* the
  probe set, which gained `d/sigma = 16, 32, 64`. Recomputing with the converged
  reference but restricted to the original standoff range gives a weighted mean
  paired difference of **`+1.313e-3` over 96 probes**, the same sign as the first
  run's `+1.377e-3`. **The sign flip is attributable to the added far-field
  standoffs, not to the reference fix**, and the first run's sign in fact agreed
  with the final near-field verdict on its own probe range. What the gate did do,
  correctly, is refuse to let a number whose uncertainty exceeded it count as
  evidence either way. Disaggregated by standoff, the outcome is unanimous on
  both sides of a crossover at `d/sigma ~ 3`: `:upstream` lower error at **12
  of 12** probes at each of `d/sigma = 0.25, 0.5, 1, 2`; `:split` lower at
  **12 of 12** at each of `4, 8, 16, 32, 64`. The aggregate statistic (56.9%
  favouring `:split`, mean paired `-1.207e-3`) hides that structure completely
  — the same aggregate-only failure R3 charged against Stage 4b. The
  `:control` family at mid-chord tracks the far-field branch, which locates the
  effect at the **partition**, not at the sheet.

  **T2 — no veto.** Paired p95 stretching differences (`:split` minus
  `:upstream`) over the 18 replicate conditions are positive but unresolved at
  the pre-registered 2-sigma threshold: interior `+2.383 +/- 1.940` (ratio
  1.23), root/tip `+0.720 +/- 0.519` (1.39), perimeter `+4.021 +/- 3.577`
  (1.12). So the §5.7 hazard is **not** resolved by this instrument, and
  `:split` is not vetoed.

  **`alpha*` — `:upstream` is the best admissible weight for pointwise
  `kappa`.** `kappa(alpha)` was verified affine to `1.7e-16`, so the
  error-minimizing `alpha` is exact: **1.45** (interior median, range
  [1.42, 1.47]) — *outside* `[0,1]`. The derived consistency weight is 1.20 and
  the second-order weight on this graded mesh is 0.23, **not** 0.5. Relative
  `kappa` error is 0.468 / 0.309 / 0.150 at `alpha` = 0 / 0.5 / 1, i.e.
  `:upstream` is 2x more accurate than `:split` here. This confirms the
  pre-registered theory result numerically: no admissible `alpha` is consistent
  on a graded mesh, and Stage 4b's "`:split` is second order" holds only in the
  refinement limit, where the grading ratio tends to 1 and both weights tend to
  1/2.

  **T3 did not fire** — its trigger was T1 being inconclusive, and T1 is
  conclusive.

  **Verdict.** The mechanical aggregate of my T1 rule points to `:split`; the
  disaggregated result points to `:upstream` in the near field and `:split` in
  the far field. Privileging the near field is **not** a post-hoc choice: Phase 1
  §5.7 and Ryan's 2026-08-01 Stage 4b charge ("whether the half-jump filament
  `:split` retains at the handoff degrades **the near field this item exists to
  clean up**") both designate it in advance. In that regime `:upstream` wins
  unanimously, and is independently corroborated by the `kappa` result. So
  **`:upstream` is retained, and the wording "retained pending a discriminating
  measurement" is now superseded: it is retained on evidence.** Two honest
  caveats stay attached: my Stage 4c pre-registration specified probes across
  standoffs but never said how to aggregate across regimes, which is a real gap
  in it; and T4 (the pre-registered Phase 4 CT comparison) remains the
  mission-level confirmation, now with a concrete prior — the near/far crossover
  means CT, a near-field-weighted integral, should favour `:upstream`.

  *Verification (six threads maximum).* Harness reruns deterministically and
  reproduces every CSV. No source file was modified, and the suites are
  unchanged from the R1/R2 entry below: wake **583/583**, replay **125/125**,
  warm start **12/12 + 14/14**, simulate **81/81**.

- **2026-08-03 — R1/R2 REMEDIATED; technical completion restored, approval
  still pending.** `PanelParticleWake` now rejects
  `SurfaceVorticityConversion` with either `unsteady_filament=false` or
  `include_final_filament=false` before allocating the panel wake or particle
  field. The legacy conversion still accepts and forwards both options. No
  handoff arithmetic, deposition rule, conversion metadata schema, or legacy
  line-policy representation changed.

  The external transfer ledger now asks `get_sources` whether a
  `FilamentWrapper` is actually emitted before including its circulation. New
  source-level tests independently reproduce both premises: the
  `unsteady_filament=false` aft face is exactly cancelled, while omitting the
  wrapper leaves the face uncancelled after a full-weight handoff. Constructor
  tests cover both invalid smooth options separately and together, and replay
  rejects smooth metadata that disables the required unsteady filament.

  Verification with `JULIA_NUM_THREADS=6`: wake **583/583**, replay
  **125/125**, warm start **12/12 + 14/14**, and simulate **81/81**. The full
  `test/runtests.jl` run reached Analytical Validation, exited 0, and contained
  zero `Fail`/`Error` entries. The 33-test Stage 0 legacy golden reference
  remains exact. R3 remains a non-blocking evidence caveat: `:upstream` is
  retained pending a discriminating measurement, not selected by the existing
  Stage 4b measurement. Phase 4 remains blocked on explicit Phase 3 approval.

- **2026-08-03 — CLEAR-CONTEXT REVIEW: approval WITHHELD.** The scoped diff was
  checked against the approved Phase 2 architecture (incl. §13 Option B and the
  §14 divergence-form amendment) and against the live source; the divergence
  ledger was re-derived independently; every verification number in the log was
  re-run.

  *What holds.* The face ledger is correct and self-consistent as derived from
  the package's contiguous ring orientation: the upstream face nets
  `(muhat_A - muhat_G)(v4-v1)`, the downstream face
  `(muhat_G - muhat_D)(v3-v2)`, each interior spanwise face is split `0.5/0.5`
  between the two panels that flank it and reassembles whole, and the root/tip
  lines `+muhat_1`/`-muhat_ncols` are exactly the two streamwise faces with no
  neighbour. The post-conversion complement is genuinely carried by the retained
  filament: `shed_wake!` copies `strength[N] -> strength[N+1]` on the shift, so
  `_final_filament_strength = -(alpha*strength[N] + (1-alpha)*strength[N+1])`
  leaves exactly `(1-alpha)(muhat_A - muhat_G)` on the panel side — verified
  against the filament's own `node(N+1,j) -> node(N+1,j+1)` orientation. At the
  first conversion `strength[N+1]` is still zero (the shift that writes it has
  not yet run), so `beta=1` deposits the physical starting vortex, confirming
  the 2026-08-01 fix. Option B is implemented as selected: no index map, source
  count, probe view, or VTK writer was touched. Persistence is strict where §9.1
  and §13.3 require it — `_deserialize_conversion` never reaches the `NoShed()`
  fallback, absent-block-selects-legacy is the only tolerated omission, and
  smooth metadata fabricates no line policies. **Verification independently
  reproduced:** `runtests_unit_wake.jl` 575/575, `runtests_unit_replay.jl`
  124/124, `runtests_unit_warmstart.jl` 12+14 = 26/26.

  *R1 (blocking) — the startup deposit is wrong when `unsteady_filament=false`.*
  `beta = handoff_before ? 1-alpha : 1` assumes the pre-handoff retained filament
  is `-strength[N+1]`, i.e. that the aft face carries an uncancelled starting
  vortex. That is only true for `unsteady_filament=true`. `PanelParticleWake`
  forwards `unsteady_filament` to `PanelWake` and the option is explicitly
  supported and unit-tested (`runtests_unit_wake.jl`, "PanelParticleWake
  forwards unsteady filament option"); with `unsteady_filament=false` the
  filament is `-strength[N]`, so the aft face carries **zero** circulation and
  the first smooth conversion injects a starting vortex that never existed.
  Measured on the standard fixture: circulation present on the aft face
  `[0, 0, 0]`, deposited `[0, -0.300, 0]` — the exact mirror image of the bug
  fixed on 2026-08-01, and invisible to the tests for the same reason (every
  smooth test runs at the default `unsteady_filament=true`; the startup testset
  even asserts `norm(startup) > 0.1`, which the broken configuration fails to
  reach because it is never constructed). Fix: either derive `beta` from the
  filament actually in force, or reject the combination at construction.

  *R2 (blocking, lower severity) — `include_final_filament=false` silently
  double-counts.* `get_sources(::PanelWake)` omits the `FilamentWrapper`
  entirely in that mode, so `_final_filament_strength` is never consulted and
  the handoff's cancelling/complement filament does not exist as a source. With
  `alpha=1` the deposited upstream face is then also still carried by the
  retained ring — every conversion, not once. The constructor accepts the
  combination (verified). The external ledger test cannot see it because `S(w)`
  calls `_final_filament_strength` directly rather than asking `get_sources`
  what the wake actually radiates. Fix: reject the combination, or make the
  ledger test source-derived. (Compounding, pre-existing and out of scope:
  `include_final_filament` is not serialized in the `PanelParticleWake` manifest
  branch at all, so it silently round-trips to `true`.)

  Both are one-line guards in the spirit of §2.2, which already raises
  `ArgumentError` when a legacy line policy is supplied alongside the smooth
  strategy. The precedent for refusing a silently-ignored option exists; these
  two are worse than ignored.

  *R3 (non-blocking, record integrity) — `:upstream` is retained by
  DEFAULT-PRESERVATION, not by measurement; the record should say so.* Three
  outcomes were possible in Stage 4b: measurement selects `:upstream`;
  measurement shows the two equivalent within resolved uncertainty; or
  measurement is inconclusive and the incumbent survives untested. What happened
  is the third. **No measurement in the Stage 4b set favours `:upstream`.**

  Of the two measurements Ryan's 2026-08-01 request named, only one returned a
  resolved sign, and it went the other way: `:split` is second order where
  `:upstream` is first (observed orders 1.982/1.991 vs 0.982/0.991, stable over
  three refinement levels; ~80x lower relative kappa error at the coarse rung —
  `refinement.csv`). The near-field axis returned mixed results and, on
  inspection, cannot resolve a sign at all (see the metric critique below).
  The verdict was produced by the *conjunction* — "`:split` must win on both" —
  and that rule was written after the data was in. A post-hoc conjunctive rule
  applied to mixed evidence can only ever return the incumbent; it has no power
  to displace anything. The operative cause of `:upstream` surviving is
  therefore that it was already the default.

  That disposition is defensible — conservatism is right when one displacing
  axis is clean and the other is muddy, and `:upstream` has an independent
  argument in its favour (it leaves nothing at the artificial panel/particle
  partition, which is what this item exists to remove). The objection is to the
  **record**, which currently reads as if the question were closed by evidence
  ("Stage 4b settles it"; "the conjunctive rule therefore fails"). It is not
  closed; it is deferred, and the one axis that answered cleanly answered
  against the retained default. This item has already been bitten by the same
  failure mode twice — the self-referential `expected_total` residual that
  *looked* like evidence, and (item 014) the `eta = 0.2` default that was never
  selected by measurement but was read by every later number as if it had been.
  Correct wording: **"`:upstream` retained pending a discriminating
  measurement."**

  *Why the near-field axis cannot resolve a sign.* The Stage 4b probe metric is

  ```
  E_rms(d) = rms_i || u_hybrid,i - u_sheet,i ||  /  rms_i || u_sheet,i ||
  E_max(d) = max_i || u_hybrid,i - u_sheet,i ||  /  max_i || u_sheet,i ||
  ```

  over the handoff family at each standoff. Five problems, in decreasing order
  of severity:

  1. **The reference is the representation being replaced, not truth.** It
     scores fidelity to the edge-jump sheet — the very near-field artefact item
     016 was commissioned to remove. `:legacy` scoring best at `d/sigma=0.25`
     (1.701 vs 1.812/1.761) is close to tautological: legacy reproduces the
     sheet's filaments. For the `:split`/`:upstream` pair the bias still runs one
     way, since `:split` retains a half-jump filament at the handoff and is thus
     the more sheet-like of the two — its near-in win is partly the metric.
  2. **M = 4.** The handoff family is two spanwise stations x two sides, on one
     grid at one `sigma`, with no repeat and no scatter estimate. The reported
     gaps are 1-3%.
  3. **Normalization drives the shape.** The curve peaks at `d/sigma = 1` and
     falls both ways; the raw deviation almost certainly decays monotonically in
     `d`, and the near-in suppression is the denominator (sheet speed diverging
     as the probe approaches the sheet), not hybrid fidelity. Numerator and
     denominator were never reported separately, so "wins near-in, loses
     mid-field" may be a statement about where the sheet is fast.
  4. **The two columns aggregate differently.** `E_rms` is RMS-over-RMS; `E_max`
     is max-of-deviations over max-of-magnitudes, not the max of pointwise
     relative error. They can rank differently for purely bookkeeping reasons.
  5. **The level is 0.8-3.3 for every strategy at every standoff.** The hybrid
     differs from the sheet by of order its own magnitude throughout the probed
     range, so this ranks schemes without establishing that any of them
     reproduces the sheet near-field.

  Amend the Stage 4b write-up accordingly, and record the open question rather
  than a closed one. A discriminating criterion still needs to be designed; the
  Stage 5 `diagnose_nearfield` panel-induced velocity-gradient distributions at
  the staged candidates are the closest existing instrument to the actual
  Phase 1 sec. 5.7 hazard (`:split`'s retained half-jump filament sitting exactly
  where fresh particles are deposited) and need no new machinery.

  Remediation: close R1 and R2 with guards plus tests that construct the
  offending configurations, amend the Stage 4b write-up per R3, then re-request
  approval. Nothing in the mathematics, the architecture conformance, or the
  persistence contract needs to change.

- **2026-08-03 — Stages 4b, 4, and 5 technically complete; approval gate
  reached.** No Phase 4 aerodynamic simulation was run.

  *Stage 4b.* Added a simulation-free `DirectBackend` comparison of a prescribed
  open smooth stretched sheet against retained panel rows/final filament plus
  one analytically removed outgoing row's particles. Both sheet sides are
  probed at `d/sigma = {0.25, 0.5, 1, 2, 4, 8}` in handoff, interior, root, and
  tip families. The compact handoff table, exact vector-circulation and impulse
  ledgers, and final-three-level refinement table live in
  `data/surface_vorticity_conversion_static/`. Legacy uses matched
  `SigmaOverlap(sigma, overlap)`. Existing external filament-ledger tests cover
  two real adjacent conversions and prove their shared face is transferred
  once across startup and steady-state branches. The centered `:split` rule
  retains its measured second-order refinement (1.982, 1.991) and coarse
  gradient advantage, but it does **not** lower handoff RMS at every standoff
  and its maximum error regresses at several standoffs. The conjunctive rule
  therefore fails: **`:upstream` remains the default**.

  *Stage 4.* `PanelParticleWake` metadata now carries strict conversion schema
  version 1. Legacy records type only and retains its historical line-policy
  fields; smooth records every constructor parameter and emits no fabricated
  line policies. Total omission selects the historical legacy default; present
  unknown/malformed schemas are hard errors. Each saved step records snapshot
  phase, active/overflow/live-row state, handoff flag/weight, conversion count,
  step and live-step identity, conversion fingerprint, and the non-source
  terminal strength needed to restore a weighted final filament without
  changing VTK layout. Replay restores and validates this state after every VTK
  load. Warm start restores it before replaying the saved step's skipped shed.
  Missing, ambiguous, mismatched, or malformed smooth state throws typed
  `WakeContinuationStateError`; historical legacy omission keeps its old
  defaults. A smooth interrupted run crossing a conversion boundary is bitwise
  identical to the uninterrupted run in panel state and active particles.

  *Stage 5.* `diagnose_nearfield=true` now evaluates a fresh direct probe system
  during preflight, sourcing only active panel-wake panels and the retained
  filament. Immutable count/min/max/mean/RMS/median/p95 summaries are stored
  separately for interior area, root/tip area, and perimeter line candidates;
  disabled conversions store `nothing`. Independent direct recomputation agrees
  to roundoff, a pre-existing particle with enormous strength has no effect,
  and non-finite diagnostic preflight remains transactional.

  *Verification (six threads maximum).* `runtests_unit_wake.jl` **575/575**;
  `runtests_unit_replay.jl` **124/124**; `runtests_unit_warmstart.jl` **26/26**;
  `runtests_unit_simulate.jl` all listed testsets pass; full `runtests.jl` runs
  through analytical validation with exit 0 and zero `Fail`/`Error` log entries.
  Stage 0's 33 exact legacy-golden assertions remain unchanged and green.

  Phase approval remains **unchecked**. Phase 4 remains blocked pending explicit
  user approval.

- **2026-08-01 — Stage 3 follow-up: startup-edge leak fixed, attribution option
  added.** Ryan's review of the divergence form raised two points; both checked
  out and the first was a live bug.

  *Bug: the starting vortex was being deleted.* The scheme deposits the whole
  upstream face and none of the downstream face, because the retained filament
  cancels the latter — but only once `particle_handoff_active` is set. On the
  **first** conversion the flag is still false, the filament is in legacy form,
  and the aft face carries an uncancelled `muhat_G - muhat_D`. Nothing deposited
  it and the row was then discarded. Measured on the standard fixture: an
  uncancelled `[0.0, -0.300, 0.0]` present and lost, with later conversions clean
  (`[0.0, 0.0, 0.0]`). **The existing conservation tests could not see it**,
  because they compare the transaction against its own `expected_total`, built by
  the same code that omitted the face — a test gap as much as a code bug. Fixed
  by depositing the aft face **whole** when `handoff_before == false`: before the
  first handoff it is the sheet's physical trailing boundary, not an interface
  with particles, exactly like the root/tip closures.

  *Option: streamwise attribution.* `SurfaceVorticityConversion(...;
  attribution=:upstream|:downstream|:split)` sets `alpha`, the fraction of the
  upstream face a conversion deposits (`1`, `0`, `1/2`). Every inter-row face is
  still deposited exactly once, in two instalments by the two conversions that
  flank it, so all three conserve exactly. `PanelWake` gained
  `particle_handoff_weight`, and `_final_filament_strength` became
  `-(alpha*strength[i_row] + (1-alpha)*strength[i_row+1])` when the handoff is
  active — an expression that reproduces both pre-existing forms at its
  endpoints, so `:downstream` is bit-identical to the legacy filament.
  `muhat_D` needed no new storage (`strength` already carries `nwakerows+1` rows;
  legacy reads the same value as `Gamma_tm1`). Default stays `:upstream`.
  Recorded as Phase 1 §9b.4 and the Phase 2 §14 amendment.

  *Verification.* Measured against an exact analytic gradient on a smoothly
  graded mesh, `:split` is **second order** (1.96, 1.99) where `:upstream` and
  `:downstream` are first (1.04, 1.02), and 40x more accurate at coarse
  resolution — the split makes the streamwise difference centered, so the leading
  stretching error cancels. This is *not* the Alternative B rejected in Phase 1
  §5.6 (that was retaining the full filament *and* depositing the complete area
  field). The cost is fidelity at the handoff: `:upstream` leaves nothing at the
  artificial panel/particle partition, `:split` leaves a half-jump filament there
  — exactly where fresh particles sit (§5.7).

  *Tests.* The decisive addition is a conservation check **independent of the
  deposition code**: every panel ring's perimeter vector sums to zero, so the
  wake's whole field-relevant content reduces to the retained filament,
  `S = sum_j filament_j * edge_j`, and the particles' gain must equal `S`'s loss
  across a real `shed_wake!` call. Verified to round-off for every attribution x
  {open, wrapping} x `nwakerows in {1,2,3}` over two consecutive conversions, so
  both the startup and steady-state branches are covered. (The driver must
  re-place node row 1 between steps, since `shed_wake!` leaves it a stale
  duplicate until the next `update_TE!`; only node row 1 moves and a complete
  ring cannot affect `S`.) Also added: the startup edge is deposited and
  `startup_edge_deposited` is true only on the first conversion; the retained
  filament matches `alpha` for each mode, with `:downstream` bit-identical to
  legacy; `:split` is second order and `:upstream` first, measured on a graded
  mesh with a strictly **linear** field (a curved one adds an opposite-signed
  truncation term that cancels against the stretching error and destroys the
  measured order — this cost two iterations to find); and `attribution`
  validation. Three pre-existing assertions were corrected rather than
  preserved, because the startup fix legitimately changes them: the `kappa`
  identity and the kappa-rule convergence test now convert on the steady-state
  branch, and the legacy-vs-smooth difference is now the attributed upstream
  share alone (both paths deposit the aft face at startup).

  `julia --project -e 'include("test/runtests_unit_wake.jl")'` → **549/549 pass**
  in `Free Wakes`, Stage 0 legacy golden reference still bit-exact;
  `runtests_unit_{simulate,replay,warmstart}.jl` and the full `runtests.jl` clean.

  *Deferred to Stage 4b (item 3 of Ryan's request):* measure `:split` against
  `:upstream` on the same sheet — near-field probe error (does the retained
  half-jump filament at the handoff degrade the near field this item exists to
  clean up?) and resolution sensitivity (does second order show up on a real
  stretched wake?). The default may change on that evidence.

- **2026-08-01 — Stage 3 complete: transaction, deposition, handoff.**

  *Source.* `src/FLOWPanel_wake.jl`: `PanelWake` gained
  `particle_handoff_active::Array{Bool,0}` (false by default, so legacy
  bookkeeping is untouched) and `_final_filament_strength` now selects the
  `unsteady_filament=false` strength row once the flag is set, cancelling the
  current final active row's own downstream edge per Phase 1 eq. (8). Added the
  `SurfaceParticleClass` enum, `SurfaceVorticityWorkspace{TF}` (reusable staging
  buffers, pre-sized from the wake geometry at construction and `empty!`-ed —
  not reallocated — per transaction), the `SurfaceVorticityPanelRecord`/
  `SurfaceVorticityConversionRecord`/`SurfaceVorticityDiagnostics` records with
  the mandatory §8.1/§8.2 fields and a cumulative ledger, and the typed
  `PanelParticleCapacityError` and `WakeConversionStateError`.
  `_convert_to_particles!(wake, ::SurfaceVorticityConversion, system)` implements
  the six-step preflight/commit protocol of §7. `src/FLOWPanel_simulate.jl`
  threads `system` into the conversion call; the legacy method ignores it.

  *Ghost sourcing is Option B*, exactly as §13.1 selected: the outgoing ghost is
  the pre-shift final active row `N` read in place, and no index map, source
  count, probe view, or VTK writer was touched.

  **Deposition is divergence form (Phase 1 §9b / Phase 2 §14, approved by Ryan
  2026-08-01).** The first cut of Stage 3 deposited
  `kappa = -n x grad(muhat)` using the Stage 2 centroid reconstruction, and
  needed two patches to work: an estimated upstream centroid at
  `nwakerows == 1` (the upstream row's nodes do not exist at conversion time),
  and a root/tip line strength extrapolated to the sheet edge (otherwise the
  full-area reconstruction plus a cell-value boundary line double-counts by one
  spanwise cell, measured: net `(0.150, -0.300, 0)` instead of
  `(0, -0.300, 0)`, accumulating every step).

  Both patches were symptoms of using a gradient reconstruction where the
  invariant that matters is circulation conservation. They are **deleted, not
  replaced**. The conversion now redistributes the stored ring assembly's own
  filaments — smear the internal partitions into area vorticity, keep the
  physical boundaries as lines:

  ```
  kappa_j = V_j / A_j,   Gamma_p = kappa_j * dA_p
  V_j = H_j + 0.5*(muhat_j - muhat_{j-1})*(v2-v1) + 0.5*(muhat_{j+1} - muhat_j)*(v3-v4)
  ```

  `A_j` is the *sum* of the subcell areas, which makes `sum(kappa_j*dA_p) == V_j`
  exactly; `kappa_j` is not projected onto the local normal (that would break the
  sum). Root/tip lines are the legacy filaments verbatim, `+muhat_1` and
  `-muhat_ncols`. Only the upstream *strength* is ever read — never the upstream
  panel's geometry — so `nwakerows == 1` (the production rotor configuration)
  needs nothing invented, and `system` is threaded in solely for
  `_get_wakestrength_mu`.

  Verified before adoption that the two rules converge to the same value: on a
  uniform mesh they are identical to 2e-15; on a smoothly graded one both are
  first-order accurate against the exact `-n x grad(muhat)` and their difference
  vanishes at first order, scaling as `(r-1)/2` in the row-to-row extent ratio.
  The trade is recorded honestly in Phase 1 §9b.3: at fixed resolution the
  reconstruction is the better *pointwise* estimator, while the divergence form
  has exactly zero conservation error at every resolution where the
  reconstruction leaks `O(r-1)` of circulation on every conversion. A per-step
  leak biases thrust cumulatively; a local kappa error does not.

  The Stage 2 reconstruction is retained and still evaluated per panel as a
  diagnostic; `kappa_difference` is a direct grid-nonuniformity measure and is
  pinned to converge under refinement.

  *Other implementation decisions.*
  - The relative Stokes residual uses an L1 accumulation of every signed
    contribution (`sum_j |H_j| + |face jumps|` plus the line magnitudes) as its
    denominator, so it stays meaningful when the ledger itself sums to zero — a
    wrapping chain or a constant field, where the norm of the total is pure
    round-off.
  - The `0.5/0.5` spanwise split moves *where* vorticity sits, never how much:
    both neighbours convert in the same transaction. The upstream face goes whole
    to the ghost because the upstream panel survives.
  - The diagnostic reconstruction invents no geometry either: at `N == 1` its
    streamwise leg is a zero displacement (reported unobservable), and a
    single-column wake with no legs at all gets a rank-0 placeholder.
  - `diagnose_nearfield=true` is rejected at `PanelParticleWake` construction
    with an `ArgumentError` naming Stage 5, rather than silently recording
    nothing.
  - Wrap detection reuses the legacy `5*eps()` node test verbatim so both
    strategies classify topology identically.

  *Tests.* New testset "Surface-vorticity conversion transaction (BRAINSTORM
  016)" in `test/runtests_unit_wake.jl` (271 assertions), covering: **exact
  circulation conservation on arbitrary geometry** — a stretched, warped,
  non-affine sheet as well as the uniform affine one, open and wrapping, at
  `nwakerows = 1, 2, 3`, checked against an independent recomputation of the
  filament content from stored strengths and nodes; `kappa` as the face content
  per unit area, equal to `-n x grad(muhat)` on the interior of a uniform affine
  fixture and to *half* the spanwise part on root/tip panels (the other half
  being the boundary line — which is exactly what closes the total);
  `Gamma_p = kappa*dA` and scalar circulation `|Gamma_p|/sqrt(dA)`; a single
  shared `sigma`; **convergence of the two kappa rules** — monotone in the
  stretch ratio, roughly halving per halving of excess stretch, with conservation
  exact at every level; root/tip lines equal to the legacy filaments exactly, and
  absent on a wrapping chain; **no upstream geometry at `nwakerows == 1`**
  (translating the body and tripling `Das` leaves every particle bitwise
  unchanged, while the shed strength does move the result) plus the
  `WakeConversionStateError` refusal when the body is absent; the legacy-vs-smooth
  net-circulation difference as the exact upstream-minus-downstream face identity;
  constant-strength elision; the handoff flag flipping and retargeting
  `_final_filament_strength`, with a legacy wake never setting it; repeated
  conversions accumulating the cumulative ledger without rewriting a committed
  record; transactional failure for capacity, geometry, and state with panel *and*
  particle state bitwise unchanged, plus the configuration rejection; every
  mandated §8.1/§8.2 diagnostic field; and `sigma`-driven subdivision refinement
  leaving the face content exactly invariant (area and `kappa` converge at `rtol`,
  since the 2x2 Gauss area of a warped panel is only second-order accurate).
  `make_conversion_fixture` gained `conversion`, `strength_fun`, and `node_fun`
  keywords (defaults unchanged, so the Stage 0 golden reference is untouched).

  Verification:
  `julia --project -e 'include("test/runtests_unit_wake.jl")'` → **468/468 pass**
  in `Free Wakes`, including the 33 Stage 0 legacy golden-reference assertions
  still bit-exact. `runtests_unit_{simulate,replay,warmstart}.jl` all pass.
  `julia --project -e 'include("test/runtests.jl")'` now runs **every** included
  file through to `runtests_analytical.jl` with zero `Fail`/`Error` lines and
  exit 0.

  **Stage R is moot.** The `data/pitching_wing_exp` blocker documented in the
  Stage 0–2 log and in the Phase 3 plan no longer exists: both
  `test/runtests_unit_pitching_wing_exp.jl` and its `include` line are absent
  from the working tree (removed in commit `905f7be`). The two-part gate
  workaround is therefore retired — run `test/runtests.jl` directly. The
  "judge by log content, not exit status" caution still stands generally,
  since `julia` exits 0 on a load error.

  *Not in this stage:* persistence of `particle_handoff_active`/
  `conversion_count`/conversion identity (Stage 4), the near-field
  velocity-gradient diagnostic (Stage 5), and any aerodynamic validation
  (Phase 4, still blocked).

- **2026-08-01 — Stages 1 and 2 complete.**

  *Stage 1 — strategy axis, no behavior change.* Added
  `AbstractPanelParticleConversion`, `LegacyEdgeJumpConversion` (the exact
  default), and `SurfaceVorticityConversion{TF}` with §2.1 validation, next to
  the `WakeSheddingMethod` hierarchy in `src/FLOWPanel_wake.jl`. Added the
  `DefaultWakeSheddingMethod` sentinel so the constructor can distinguish
  "caller said nothing" from "caller explicitly asked for the legacy default":
  `_resolve_line_policy` turns it into a fresh `OverlapPPS(1.3, 2)` for the
  legacy strategy, keeps an explicitly supplied policy, and raises
  `ArgumentError` if either line policy is supplied alongside the smooth
  strategy (§2.2). The sentinel itself throws if it ever reaches
  `_shed_particles!`. `PanelParticleWake` gained `conversion`,
  `conversion_workspace`, `conversion_diagnostics`, and `conversion_count`
  (three new type parameters); legacy instances store `nothing` for workspace
  and diagnostics. `_convert_to_particles!(wake)` now dispatches on
  `wake.conversion`, and the legacy method is the previous body verbatim.

  *Stage 2 — reconstruction core.* Pure, simulation-free functions per §5.2–§5.4:
  `_deterministic_tangent_basis`, `SurfaceGradientResult`,
  `_reconstruct_surface_gradient` (scaled `2x2` SVD with the rank-2/1/0 policy;
  the minimum-norm branch falls out of the same truncated-SVD sum),
  `_surface_vorticity` (the `kappa = -n x grad(muhat)` sign), the bilinear panel
  family (`_bilinear_position`, `_bilinear_derivatives`, `_bilinear_normal`),
  `_subdivision_counts`, `_subcell_area` (`2x2` Gauss on the Jacobian norm),
  `_validate_wake_panel`, and the typed `WakeGeometryError`.

  Verification (`julia --project -e 'include("test/runtests_unit_wake.jl")'`):
  **197/197 pass** in `Free Wakes` — 33 legacy golden reference (still bit-exact
  after the dispatch refactor, which is the point of Stage 0), 23 strategy axis,
  51 reconstruction core. The reconstruction tests cover constant and affine
  planar fields (streamwise, spanwise, combined), the exact package sign and its
  antisymmetry, rigid-rotation covariance of both the gradient and the derived
  surface vorticity, rank-1 minimum-norm and rank-0 behavior with observable
  flags, rejection of a vanishing metric scale and non-finite stencils, exact
  area partitioning and quadrature convergence on a warped panel, the
  subdivision floor and stretched-panel counts, and rejection of coincident,
  collapsed, bow-tie, and non-finite panels.

  Full-suite gate: `test/runtests.jl` passes lines 7–20 (all unit files) with no
  failures. **Pre-existing, unrelated blocker:**
  `test/runtests_unit_pitching_wing_exp.jl:4` aborts on a missing data asset
  (`data/pitching_wing_exp/load.jl`; the directory is absent from the repo and
  is not a pending deletion), which stops `runtests.jl` at line 21. Lines 22–27
  were therefore run separately and all pass. Note also that `julia` exits 0
  here even on such a load error, so the gate must be judged from log content,
  not exit status.

  *Design notes verified for Stage 3 (Option B mechanics).* Per-step ordering is
  `update_TE!` (`FLOWPanel_simulate.jl:460`, placing wake node row 1 at
  `body TE + Das`) → solve → `propagate!` → `propagate_kinematics!` →
  `shed_wake!` (:1125-1127). At conversion time the outgoing ghost is active row
  `N`: node rows `N`, `N+1` and strength row `N`. For `N >= 2` its streamwise
  upstream neighbour is existing active row `N-1`, already in storage — nothing
  is staged. For `N == 1` the upstream neighbour is the row the shift is about
  to create; `shed_wake!(::PanelWake, ...)` (:1190-1212) does **not** build new
  row-1 nodes, it copies row 1 into row 2 and writes only
  `strength[1,1,:] = _get_wakestrength_mu(...)`. So the staged upstream sample
  is: position = centroid of node row 1 (the current `TE + Das` line, i.e. the
  ghost's own upstream edge), strength = `_get_wakestrength_mu(system, ...)`.
  That is a genuine one-sided two-panel stencil with two distinct strengths at
  two distinct positions, and it introduces no new `Das` dependence beyond the
  geometry that already defines the ghost. **Consequence:** `system` must be
  threaded into the conversion, so Stage 3 should give the smooth strategy the
  entry point `_convert_to_particles!(wake, conversion, system)` and pass
  `system` from `shed_wake!(::PanelParticleWake, ...)` (:1224-1236); the legacy
  method can ignore it.

- **2026-08-01 — Stage 0 complete.** Phase 2 approval recorded across the hub,
  the phase document, and `BRAINSTORM/INDEX.md`, with the §13 amendment
  (ghost-storage Option B; the two §10/§9.1 corrections). Added
  `make_conversion_fixture` to `test/test_helpers.jl` — a deterministic
  `PanelParticleWake` whose panel-wake nodes and strengths are written directly
  so `_convert_to_particles!` can be exercised **with no time stepping** — and
  captured its output as
  `test/data/legacy_wake_conversion_reference.jl` for `wraps ∈ {false, true}`
  × `nwakerows ∈ {1, 2, 3}`. The new testset "Legacy edge-jump conversion
  golden reference (BRAINSTORM 016)" in `test/runtests_unit_wake.jl` compares
  positions, vector strengths, smoothing widths, the vestigial scalar
  circulation, and particle counts by **exact equality**, and asserts that the
  wrapping chain deposits fewer particles than the open chain (no root/tip
  closure). This closes the pre-existing coverage gap in which nothing pinned
  `_convert_to_particles!` directly. No source code has been modified yet.

  Command and result:
  `julia --project -e 'include("test/runtests_unit_wake.jl")'` →
  **123/123 pass** (`Free Wakes`), of which 33 are the new reference testset.

- **2026-07-29 — Phase seeded.** Blocked pending explicit Phase 2 approval;
  no implementation has started.
