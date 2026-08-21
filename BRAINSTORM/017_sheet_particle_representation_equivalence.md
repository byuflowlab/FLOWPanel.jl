# Sheet/Particle Representation Equivalence: When Do Particles Reproduce the Panel Wake They Replace?

## The question

Under what **quantitative conditions** does the particle representation of a wake region
reproduce the panel-sheet representation it replaced, to a stated tolerance? And does a
**convergent, physically consistent `PanelParticleWake` prescription** follow from that
map?

The starting axes are: particle core size σ/c, conversion distance from the body d/c,
particles per shed segment (PPS), overlap, and shedding scheme. **This list is explicitly
open-ended (Ryan, 2026-07-30): other conditions may be identified and explored as the study
proceeds, and the map extended accordingly.** One already identified is the **shedding
*policy*** — which quantity is held fixed as others refine:

- `SigmaPPS(σ, pps)` fixes σ and the particle count per segment, so particle spacing h is
  tied to the *segment length* — under σ-refinement at fixed PPS the overlap ratio σ/h
  **falls**;
- `SigmaOverlap(σ, overlap)` fixes σ and the overlap ratio, deriving PPS so h ≈ σ/overlap —
  spacing refines *with* σ;
- `OverlapPPS(overlap, pps)` (the legacy default) ties σ to the segment length, giving
  non-uniform σ across particle families (see below).

**VPM convergence theory constrains this choice**: classical results (Beale–Majda-type)
require the inter-particle spacing to decrease **faster** than σ (h/σ → 0, e.g. h ~ σ^q
with q > 1) for provable convergence of the regularized particle method. At fixed PPS,
σ-refinement moves *away* from the convergent regime; fixed overlap is the marginal case;
theory prefers overlap *growing* under refinement. Every refinement ladder in this study
must therefore state its (h, σ) path, and `SigmaOverlap` (or an explicit h ~ σ^q schedule)
is the theoretically defensible refinement policy — `SigmaPPS` remains useful for
*isolating* mechanisms (σ varied at fixed positions), not for claiming convergence. Note
the tension with item 012: heavy overlap is exactly the regime that breaks conditioning
and relaxation there — the admissibility map must expose this trade-off rather than hide
it.

Tolerances, fixed in advance (Ryan, 2026-07-30):

- **|ΔCL| ≤ 0.5%** vs the pure-sheet control (settled state), AND
- **per-station spanwise loading / circulation: max |ΔΓ(y)| / max|Γ| ≤ 1%** — integrated CL
  can hide compensating spanwise errors, so it is not sufficient alone. The spanwise TE
  μ-jump (Γ_TE(y)) is the quick-access metric.

**Wing first, rotor later.** The study runs on a steady/settled rectangular wing until the
admissibility map and mechanism are established; only then transfer to the rotor. The rotor
adds physics the wing lacks — the blade chops into its own wake — which is expected to
disturb any wing-derived prediction and is deliberately deferred (see Phase D).

## Why it matters

This item inherits the surviving core of item 014. There, `Das`/η was settled (physical
offset length on a log plateau, <1% model-form uncertainty), but the residual finding is
much larger and unowned until now:

- Rotor CT rises monotonically with the fraction of near wake represented as panel sheet vs
  particles, spanning **0.034–0.110 (3.3×)** across plausible discretizations — with
  experiment (0.072) *inside* the range. 006's 0.0713 cannot be called a prediction until
  this is resolved.
- The `nwakerows` (sheet-buffer) ladder is **non-monotone** (1/4/8/16/36 rows →
  0.07133 / 0.07431 / 0.07506 / 0.07304 / 0.07049): sheet extent is **not a convergent
  axis** at current σ, so no buffer prescription can be adopted (Ryan's 2026-07-30 ruling
  in 014).
- The wing sheet/particle split ladder (014, 2026-07-29) localized the sensitivity to
  **σ/chord ~ 1**, with a sign flip: chord-scale-σ particles near the body over-induce,
  farther away they under-induce. At σ/c = 0.08 the split is irrelevant (<0.1% over a 32×
  buffer range); at σ/c = 1.5 (rotor-like) the spread is 3.8% on a fixed-incidence wing —
  amplified on the rotor by inflow→loading feedback.

## What is established (imported from 014/012, verified in source)

- **The two wake representations are consistent when geometry matches**: the straight-sheet
  panel wake converges to the semi-infinite answer to 0.014% (014 verdict, 2026-07-29). The
  problem is therefore not a representation *inconsistency* — it is the particle
  approximation of the sheet at finite σ and finite quadrature.
- **`OverlapPPS` gives non-uniform σ across particle families.** σ = segment_length ·
  overlap / pps per converted segment (`_shed_particles!`,
  `src/FLOWPanel_wake.jl:659-696`): trailing particles use the streamwise row spacing
  (≈U·dt), unsteady particles use the spanwise panel width (≈b/n_span). A controlled σ
  sweep must use a **fixed-σ method** — `SigmaPPS(sigma, pps)` for mechanism isolation,
  `SigmaOverlap(sigma, overlap)` for convergence-claiming ladders (see the policy
  discussion under "The question"). Both are exported and accepted by the
  `method_trailing`/`method_unsteady` slots of `PanelParticleWake`
  (`src/FLOWPanel_wake.jl:659-696`, `:1049`).
- **Current shedding seeds trailing particles exactly along the streamwise vortex filament
  of the converted row** (`_convert_to_particles!`, `src/FLOWPanel_wake.jl:1276-1319`
  treats each constant-strength wake panel as a vortex ring; strength jumps become line
  particles laid on the panel edges). Near-duplicate/inline placement is structural (014 /
  overlap memory), and corrected-Pedrizzetti relaxation acts on these near-aligned
  particles every step.
- **Conversion happens only at the single oldest row** (`nwakes[]`), from inside
  `shed_wake!` (`src/FLOWPanel_simulate.jl:1223`, gated on buffer full). There is no
  offline conversion API — but `_convert_to_particles!` is a pure function of wake state,
  so a harness-side one-shot row-range converter is small.
- **Item 012 owns the σ/h (overlap) regime**: measured σ/h ≈ 4.2 in the settled rotor wake
  (~4× nominal), which is the same under-resolution seen here as σ ≈ 1.5 local chords at
  shedding. 017 supplies the quantitative wing oracle 012 lacked.
- **FMM confound**: `src/FLOWPanel_liftingbody.jl:996` inflates the FMM panel radius by
  |Das|; every rotor η run used `FastMultipoleBackend`. All 017 wing work uses
  `DirectBackend`; the rotor discriminator is Phase D / T4.

## Candidate mechanisms to discriminate — the study's targets

- **M1 — kernel smoothing.** The regularized particle kernel smears the sheet's induced
  field at scale σ; error at evaluation distance d expected ~O((σ/d)²), sign-changing with
  d (over-induction near, under-induction far — matches the 014 wing table's sign flip).
- **M2 — quadrature/lumping.** Discrete vortons approximate a continuous sheet; error
  controlled by PPS (particle spacing along the segment) independently of σ.
- **M3a — filament-alignment × relaxation.** Particles seeded on/near the preceding row's
  filament line, then acted on by corrected-Pedrizzetti; a *dynamic* error source invisible
  to any static field comparison. (This is the mechanism item 016's area-quadrature
  shedding is designed to remove.) **Caveat (Ryan, 2026-07-30): relaxation is not entirely
  artificial** — it makes particle dynamics more consistent with the governing equations —
  so a relaxation-OFF arm is NOT a clean discriminator and is omitted for now. If M3
  attribution becomes load-bearing, a future investigation should design a discriminator
  that respects relaxation's physical role (e.g. seeding-geometry perturbation at fixed
  relaxation, or the 016 shedding change itself as the intervention).
- **M3b — trajectory/rollup divergence.** The particle wake *geometry* evolves differently
  from the sheet wake under mutual induction — present even with ideal strengths and
  relaxation untouched. Routes to a different remedy than M3a (σ management / convection,
  item 012, vs shedding, item 016), so the two must not be lumped. Discriminator: the
  matched-geometry arm — reuse 014's `freestream_convection` mode for the particle wake if
  it extends to `PanelParticleWake` cheaply, so sheet and particle wakes share identical
  geometry and only strength dynamics differ.
- **M4 — core–body overlap.** σ ≈ chord means a newly shed particle's core engulfs the
  generating blade section; direct smoothing of the body's near field.

M1/M2/M4 are static (measurable offline); M3a/M3b are dynamic (appear only in time-marched
runs). The program below separates them by construction: Phase A measures the static
mechanisms in isolation; Phase B measures the total in dynamics; the difference is the
combined M3a+M3b budget (separated only as far as the matched-geometry arm allows).

## Study program

### Phase 0 — harness hardening (prerequisite; `examples/` only, no `src/` default changes)

In `examples/suddenly_started_wing.jl`:

1. `SSWConfig` fields: `shed_method::Symbol = :overlap_pps` (`:overlap_pps | :sigma_pps |
   :sigma_overlap`) and `sigma_over_c::Float64 = NaN` (required by the sigma methods). The
   particle branch of `prepare_suddenly_started_wing` (:431-447) builds
   `method_trailing`/`method_unsteady` from these (same method in both slots).
2. ENV hooks in `_ssw_config_from_env` (:966): `SSW_WAKE_MODEL`, `SSW_PANEL_ROWS`,
   `SSW_SHED_METHOD`, `SSW_SIGMA_OVER_C`, `SSW_PPS`, `SSW_OVERLAP`, `SSW_WAKE_CORE`
   (→ `wake_core_over_c`), `SSW_MAX_PARTICLES`.
3. Case-tag tokens in `_ssw_case_tag` (:518): wake core (when ≠ 1e-3) and shedding method
   (σ/c + pps for `SigmaPPS`), so no two distinct configs share a tag. **This closes the
   silent-resume defect class that bit 014 twice** (mesh revision; wake core).
4. Spanwise outputs per case at the settled tail: `gamma_te.csv` (span station, Γ_TE from
   the TE doublet-strength jump) and per-station loading (reuse the
   `SpanwiseLoadingMonitor` pattern, `code_audit/scripts/task1_steady_wake_consistency.jl:127`).

### Phase A — offline static-conversion field diagnostic (T1)

New driver `examples/ssw_representation_probe.jl`. Local, `DirectBackend`, ≤6 threads,
minutes per case.

1. Settle a **pure-panel** wing wake to t* = 20 (`wake_model=:panel`, all-steps buffer),
   in **two arms, both required**: `freestream_convection=true` (straight sheet — exactly
   reproducible geometry, cleanest scaling fits) and the rolled-up default. The rolled-up
   arm is not optional decoration: Phase B's dynamic runs roll up, so **the pre-registered
   Phase B prediction must be sourced from the rolled-up arm** — predicting a rolled-up
   ladder from straight-sheet fits would build a geometry mismatch into the closure test.
   Tip-vortex note (Ryan, 2026-07-30): a constant-doublet wake panel is exactly a vortex
   ring on its edges, so even the straight sheet carries a genuinely *concentrated* line
   filament at its tip edge (strength = the tip panel column's μ) — what the straight arm
   lacks is the *accumulation and geometry* of roll-up (inboard filaments spiraling into a
   contracted core carrying a growing fraction of total Γ). So the straight arm tests
   particle-vs-line-filament locally; only the rolled-up arm tests particle-cloud-vs-
   accumulated-core, the regime the rotor lives in. Settledness gated by drift+ripple
   (`settled_stats` pattern, `task1_steady_wake_consistency.jl:205`).
2. **One-shot row-range converter** (harness-side, adapted from `_convert_to_particles!`):
   rows r1..r2 → trailing/unsteady/tip particle sets in a fresh `FLOWVPM.ParticleField`
   via `SigmaPPS(σ, pps)`. No time stepping. **Conservation checks**: (0th moment)
   per-row particle Γ ledger equals the sheet row's strength jumps; (1st moment) per-row
   impulse Σ x×Γ vs the sheet row's — Γ-exact lumping still moves the impulse through
   position quantization, and the far field is governed by it, so this line yields M2's
   leading-order coefficient for free.
2a. **Frozen-wake re-solve (the quantitative A→B bridge).** For each conversion case,
   re-solve the body with rows r1..r2's sheet influence replaced by the particle set's
   influence in the RHS (same G matrix / factorization, new RHS — one linear solve per
   case). This yields **ΔCL and ΔΓ_TE(y) in exactly Phase B's metrics**, with the
   one-pass solve feedback (field error at CPs → Γ change) included. Field-error maps
   (step 3) alone are incommensurate with Phase B's outputs; this step is what makes the
   pre-registered prediction quantitative rather than qualitative.
3. Compare **region-only induced fields** — u and ∇u (the relaxation input) from sheet rows
   r1..r2 vs their particle replacement — at (a) body control points, (b) a Cartesian probe
   plane through the mid-span chord line, (c) a spanwise line just aft of the TE. Probes:
   `FastMultipole.ProbeSystem` + `influence!` (template
   `src/FLOWPanel_simulate_monitors_fieldprobe.jl:171-200`; sheet/particle entity split as
   in `_fieldprobe_entities`, same file :222-228).
4. Axes, varied **independently**: d/c (which rows), σ/c ∈ {0.05…1.5}, PPS at fixed σ
   (isolates M2), and re-emission at identical positions/strengths with σ varied alone
   (isolates M1). Metrics: RMS/max relative error at body CPs, error-vs-d curves, fitted
   scaling exponents.
5. **Trivial-limit validation before any science**: σ→small with PPS→large must reproduce
   the sheet field to <0.1% at all probes. This limit is *exact*, not approximate: a
   constant-doublet panel ≡ vortex ring on its edges, so line particles on those edges
   with σ→0, PPS→∞ converge to the sheet field identically — a persistent near-field
   discrepancy in this gate is a harness bug, not physics.

Deliverable: an **e(σ/c, d/c, PPS) error map** with per-mechanism scalings (fitted in the
transferable variables of Phase B's mapping deliverable, i.e. σ/d and h/σ, not raw PPS),
and a **written prediction of the Phase B split-ladder table** (pre-registered, sourced
from the rolled-up arm via the step-2a frozen-wake re-solve).

**Discussion deliverable — which moments should shedding target?** (Ryan, 2026-07-30; not
current work.) For 017's equivalence question the reference is necessarily the
piecewise-constant sheet (that is what the particles replace, and the step-2 ledger closes
against it). But for *method accuracy* the fork is real: the sheet's own moments already
differ from the underlying continuous vorticity distribution at O(h²), so a shedding rule
matching the panel model's moments exactly inherits the panel discretization error
faithfully, while one targeting the *continuous* distribution's moments could partially
cancel it (superconvergent-quadrature analogy). The two targets differ at the same order
as M2's lumping error, so it is testable within Phase A's machinery: compute both
discrepancies (particles-vs-sheet and sheet-vs-continuous, the latter from a spanwise
μ-refinement reference) and report which dominates, with a recommendation. Open to either
conclusion; feeds 016's quadrature design if the continuous target wins.

### Phase B — dynamic split ladder under controlled σ (T2)

Extend `examples/ssw_sheet_particle_split.jl` (ENV: `SSWSPS_METHOD`, `SSWSPS_SIGMAS`,
`SSWSPS_PPS`): `SigmaPPS` with σ/c ∈ {0.08, 0.15, 0.3, 0.6, 1.2} × buffer ∈
{0.25, 1, 4} c × PPS ∈ {1, 2, 4} — anchor corners plus a dense σ axis, not the full
factorial (~20–25 runs; **moderate** effort, local `DirectBackend`, ≤6 threads, minutes
each). Metrics per case: tail CL vs pure-panel control, **spanwise loading Δ**,
**Γ_TE(y) Δ**, particle count, a passive relaxation-magnitude diagnostic (how much
corrected-Pedrizzetti moves strengths per step — an M3a *indicator* only; no
relaxation-OFF arm, per the M3a caveat, and production keeps stock relaxation per the
standing policy), and the **startup CL(t) transient** (Wagner-like rise) scored
sheet-vs-particle — free with each run, and it is the only place the steady wing exercises
the **unsteady (spanwise) particle family**, whose settled-tail strengths are ≈0.

Deliverable: the **admissibility contour** — the (σ/c, buffer, PPS) region satisfying
|ΔCL| ≤ 0.5% AND per-station |ΔΓ| ≤ 1%. The output table must flag, cell by cell, which
cells are **admissibility-eligible** (on a theory-consistent (h, σ) path) vs
**mechanism-isolation only** (e.g. σ/c = 1.2 at pps = 1 sits far from every defensible
refinement policy and may not be quoted as a convergence point).

**Shedding-policy arm**: repeat the σ axis under `SigmaOverlap(σ, overlap)` at fixed
overlap (and, if warranted, an h ~ σ^q schedule with q > 1) alongside the `SigmaPPS` runs.
`SigmaPPS` isolates mechanisms; `SigmaOverlap` is the policy under which convergence claims
are theoretically defensible (h/σ must not grow under refinement — see "The question").
If the two policies disagree at matched σ, that disagreement is itself a finding (spacing
error M2 vs smoothing error M1) and the admissibility map must be stated per policy.

Also in Phase B:

- **Close 014's wake-core gap**: rerun the F2 core probes (1e-4c / 1e-2c) with the new
  `SSW_WAKE_CORE` hook + tag token on the rev-2 mesh (minutes).
- **Mechanism closure test**: compare the measured ladder against Phase A's pre-registered
  prediction (rolled-up arm, frozen-wake re-solve — same metrics, matched geometry).
  *Match* ⇒ the static mechanisms (M1/M2/M4) explain the split sensitivity — named and
  closed. *Mismatch* ⇒ the dynamic route (M3a/M3b) is implicated; the matched-geometry
  particle arm (if wired) then splits geometry divergence (M3b → 012) from
  strength-dynamics effects (M3a → **strengthens the recommendation to finish item 016**).
- **Transferable-variables mapping (deliverable, Ryan 2026-07-30)**: a written mapping
  from simulation parameters to the physics variables the map is stated in — σ/c (and
  σ/d), **h/σ** per particle family (streamwise h = U·dt/pps; spanwise h = b/n_span;
  rotor: h = Δψ·r/pps at station r, local chord for σ/c), d/σ. This absorbs the hidden
  **dt axis** (at fixed PPS, streamwise h is set by dt, so any PPS-labelled map silently
  depends on dt) and makes Phase D a *test* of the wing map rather than a re-fit —
  pre-register the rotor lookup before any Phase D run.

### Phase C — settledness rigor (cross-cutting)

Every quoted number from a settled tail with explicit drift+ripple gates; wake-extent-in-
chords diagnostics recorded per case (`task1_steady_wake_consistency.jl:157-181`);
`assert_ssw_mesh_symmetry` at sweep start; even `n_span`.

### Phase D — rotor transfer (DEFERRED; staged here, executed only after the wing map holds)

Pre-registered for a later agent:

- **T3 — σ-refined `nwakerows` ladder** {1, 4, 16} at wing-map-admissible σ.
  *Pre-registered prediction:* if M1/M2 govern, the ladder becomes flat or
  monotone-saturating (the non-monotonicity comes from the handoff moving through the
  near-over-induce / far-under-induce sign-flip region at chord-scale σ). First harvest
  the in-flight 014 jobs — 12955430 (`p2e_nrows72_das1p0`), 12943696 (`p2e_sigF_nofilt`,
  σ ≈ 0.37 local chords), 12950996 (`p2e_nt72_das2p0_ov6`) — they are prior points on
  these axes.
- **T4 — `DirectBackend` discriminator** (cheap rotor-side hygiene; may run in parallel
  with the wing study at Ryan's discretion): two η-ladder points on `DirectBackend` to
  close the FMM |Das|-radius mechanism.
- **T5 — cost model**: particle count and wall time vs σ/c on the rotor
  (σ = Δψ·r·overlap/pps gives the count scaling in closed form; one calibration point
  suffices). Verdict: σ-refinement affordable, or **sheet-until-harmless** (buffer length
  read off Phase A's e(d) curve — admissible only if T3 shows the nrows axis convergent at
  adequate σ), or route to **item 016** (smooth shedding) and/or **item 012** (σ
  management).
- **Not covered by the wing** (named transfer caveats): (a) the rotor blade chops into its
  own wake; (b) the **unsteady (spanwise) particle family** is only weakly constrained by
  the steady wing (settled strengths ≈0; the Phase B startup transient is its sole wing
  test), while the rotor's curved rows load it continuously. Any wing-derived prescription
  must be re-verified under both; this is the bridge to cross after the map holds, not
  before.

## Acceptance

017 is answered when there is:

1. a **quantitative admissibility map** — which conditions (starting axes σ/c, d/c, PPS,
   shedding policy; extended with any further conditions the study identifies) meet the
   0.5% CL and 1% per-station Γ tolerances on the wing, with each convergence claim stating
   its (h, σ) refinement path consistent with VPM theory (h/σ non-increasing);
2. **mechanism attribution**: the split error assigned to named mechanisms (M1/M2/M3a/
   M3b/M4) with measured scaling laws that *predict* the dynamic split-ladder table
   (Phase A step-2a frozen-wake re-solve → Phase B closure, matched rolled-up geometry),
   plus the **transferable-variables mapping** and the **moment-target discussion**
   (particles-vs-sheet vs sheet-vs-continuous) as written deliverables; and
3. a **prescription** — σ/c bound, sheet-buffer length, or scheme change — under which the
   wing's split sensitivity vanishes at rotor-like parameters, together with its rotor cost
   estimate. A negative outcome is acceptable and valuable: a demonstration that no
   affordable legacy-shedding configuration meets tolerance, with the cost curve, in which
   case the item's recommendation is to complete **item 016** (surface-vorticity shedding,
   fully designed, awaiting approval — NOT to be implemented under this item) and/or 012's
   overlap management.

Rotor transfer (Phase D) is a follow-on gate on top of 1–3, not a precondition for calling
the wing study complete.

## Cross-references

- `BRAINSTORM/014_first_wake_row_offset_selection.md` — origin: the η/`Das` resolution, the
  0.034–0.110 CT range, the non-monotone nwakerows ladder, the wing split table, Ryan's
  2026-07-30 ruling. Rotor run ledger:
  `logs/dji_convergence_20260722/phase_02e_unsteady_ct.md`.
- `BRAINSTORM/012_robust_resolution_overlap_management.md` — owns σ/h overlap management;
  017's admissibility map is the wing oracle for it.
- `BRAINSTORM/016_surface_vorticity_particle_shedding.md` — the designed alternative to
  filament-aligned shedding (M3 lever). Future recourse only; referenced, not implemented
  here.
- `examples/suddenly_started_wing.jl`, `examples/ssw_sheet_particle_split.jl` — the
  platforms; `code_audit/scripts/task1_steady_wake_consistency.jl` — settledness/extent
  diagnostics to borrow.
- Plan of record for staging: `/Users/ryan/.claude/plans/recursive-splashing-hopcroft.md`
  (session-local; the content above is self-contained).

## Status

**TECHNICALLY COMPLETE 2026-07-30.** Phases 0–C executed; see `### Technical outcome`
near the end of this file for the results and the transferable prescription. Phase D
(rotor transfer) remains deferred. The staging text below is retained as the record of
the plan as approved.

**Staged 2026-07-30.** No phase executed yet. Phase 0 (harness hardening) is the first
action; Phases A/B are local and cheap; Phase D is deferred until the wing map holds.

**2026-07-30 clear-context design review (approved by Ryan, folded in above).** Design
upheld; source-code claims spot-checked and accurate. Changes: (1) Phase A gains the
**frozen-wake re-solve** (step 2a) so the pre-registered prediction is in Phase B's own
metrics with solve feedback included; (2) M3 split into **M3a** (relaxation-coupled) and
**M3b** (trajectory/rollup divergence) — relaxation-OFF discriminator REJECTED by Ryan
(relaxation is partly physical; future-investigation note instead), matched-geometry arm
kept as the M3b probe; (3) rolled-up arm promoted to the required prediction source, with
the tip-vortex clarification (straight sheet has a concentrated tip *filament*; roll-up
adds accumulation + geometry); (4) map to be stated in transferable variables (σ/d, h/σ
per family) with a written simulation-parameter→physics-variable mapping deliverable,
absorbing the hidden dt axis; (5) impulse (1st-moment) ledger added, plus the open
**moment-target question** (match the sheet's moments or the continuous distribution's?)
as a discussion deliverable; (6) startup-transient scoring added as the unsteady-family
wing test, and that gap named as a Phase D transfer caveat; (7) effort estimates,
trivial-limit rationale, and per-cell admissibility-eligibility flagging added.

## 2026-07-30 implementation and execution ledger

Phase 0 and the Phase A/B drivers are implemented without changing production
`src/` defaults or item 016:

- `examples/suddenly_started_wing.jl`: fixed-σ shedding methods, ENV hooks,
  collision-proof wake tags (full policy, core, capacity), TE-circulation and
  spanwise-loading artifacts, CL/Γ settledness gates, and a passive
  corrected-Pedrizzetti recording wrapper.
- `examples/ssw_representation_probe.jl`: complete-ring row-range conversion,
  DirectBackend velocity/gradient probes, shell and cumulative-tail maps,
  final-three-state serialization, circulation/impulse ledgers, and three-state
  frozen BDF2 re-solves. Rolled-up predictions are immutable at
  `data/ssw_representation_probe/phase_b_prediction.csv`.
- `examples/ssw_sheet_particle_split.jl`: the 23-cell SigmaPPS set,
  fixed-overlap SigmaOverlap ladder, tail extension through t*=40, transferable
  variables, eligibility/admissibility flags, relaxation diagnostics,
  checkpoint resume, and a Phase-A-prediction hard gate.
- `examples/run_ssw_representation_study_hpc.slurm.sh`: six-thread, 64-GB,
  12-hour checkpointed launcher. Isolated cluster staging is
  `/home/rander39/projects/FLOWPanel017`.

Focused verification:

```text
test/runtests_example_suddenly_started_wing.jl: 100/100 pass
test/runtests_unit_wake.jl: 90/90 pass
test/runtests_unit_simulate.jl: all testsets pass
test/runtests_unit_postprocess.jl: 406/406 pass
```

The simulate suite's standalone invocation used one Julia thread and
`using LinearAlgebra: cross`: FLOWVPM's zero-particle direct-SFS test constructs
a zero-step threaded range above one thread, and the standalone file assumes
`cross` came from the broad runner. Neither issue is caused by 017.

### Phase A gate and first results

The real AR=6, n_span=24, t*=20 trivial-limit gate was run in both straight and
rolled-up arms. Both `(σ/c,PPS)=(0.01,32)` and `(0.005,64)` pass every probe
family. Worst normalized velocity error is `3.23e-6`; worst normalized
velocity-gradient error is `2.14e-4`, below 0.1%. Vector-circulation and impulse
ledgers close at O(1e-15). The aft-span line is offset 0.15c normal to the
sheet: the original on-edge line was singular and tested finite-core
regularization rather than converter correctness.

The full shell+cumulative field table has 1512 rows in
`data/ssw_representation_probe/field_metrics.csv`. A representative
rotor-bridge cell (rolled-up cumulative tail, σ/c=0.3, PPS=2, d/c=0.25) has
body-control-point max errors `|Δu|/U=0.1016` and
`|Δ∇u|/(U/c)=0.7966`; its 30,336 particles retain ledger errors below
`4e-14`.

The midpoint SigmaPPS rule preserves total vector circulation and first impulse
exactly for each straight subsegment, independent of PPS. Thus the proposed
"impulse quantization" coefficient is identically zero; the leading lumping
defect is a second/higher moment. The particles-vs-sheet versus
sheet-vs-refined-continuous comparison remains required before recommending a
moment target.

### Pre-registration, smoke closure, and core checks

The immutable rolled-up frozen table contains 168 cells. At
`(σ/c,PPS,d/c)=(0.3,2,0.25)`, frozen replacement predicts
`ΔCL=-18.356%` and max-station `ΔΓ=19.631%`; the settled dynamic smoke measures
`-0.2580%` and `0.2602%`. The residuals exceed the 0.25/0.5-percentage-point
closure gates by orders of magnitude. This pre-registers the dynamic M3 branch
for the full ladder; matched geometry is required before final attribution.

Corrected-mesh core smoke checks at the same cell:

| wake core/c | ΔCL | max ΔΓ |
|---:|---:|---:|
| 1e-4 | -0.1806% | 0.1875% |
| 1e-3 | -0.2580% | 0.2602% |
| 1e-2 | -0.4554% | 0.4488% |

All CL/Γ histories pass the 0.1% drift and 0.25% ripple gates. These are
mechanism-isolation SigmaPPS cells, not admissibility claims.

### Transfer variables and rotor cost projection

The output table maps parameters to `σ/d`, `d/σ`, trailing
`h/σ=dt*/(PPS·σ/c)`, and unsteady `h/σ=(AR/n_span)/(σ/c)`. For a rotor,
trailing `h=Δψ r/PPS`, unsteady spacing is the blade-span segment length, and
all σ ratios use local chord.

The existing rotor calibration is 41,778 particles at settled step 700 with
PPS=2 and σ/c≈1.5, with NT=36/depth-4 wall time near 9 hours. Holding retained
wake volume fixed gives the analytic legacy-shedding projection

$$
N_p \approx \frac{62{,}667}{\sigma/c}, \qquad
t_{\rm wall} \approx \frac{13.5\ {\rm h}}{\sigma/c}.
$$

Thus σ/c=0.3 projects ≈209k particles / 45 h; 0.15 projects ≈418k / 90 h;
0.08 projects ≈783k / 169 h. The finest point exceeds the configured 500k
capacity, and the latter two are not practical production runs. This is an
analytic Phase-D projection, not a rotor run.

### Final execution record

The first cluster attempt (`12960755`) stopped before Phase B because Julia
1.11 cannot deserialize the locally generated Julia-1.12 snapshots. The
launcher now forces host-native snapshot regeneration. The corrected full job
`12960766` completed all 38 particle cells plus the panel control in 3:05:54
on six threads. Moment-refinement job `12960933` completed in 1:14:49 and
matched-geometry job `12977811` completed in 0:12:57. Rotor runs were not
performed. A provenance-only recovery job (`12978237`, 2:50:55) persisted
relaxation and measured-extent diagnostics that the initial resumable
artifacts omitted; its 38 CL, Γ, and loading results are bit-for-bit identical
to `12960766`.

Commands used:

```text
JULIA_NUM_THREADS=6 julia --project=. -t 6 test/runtests_example_suddenly_started_wing.jl
JULIA_NUM_THREADS=6 julia --project=. -t 6 test/runtests_unit_wake.jl
JULIA_NUM_THREADS=1 julia --project=. -t 1 -e 'using LinearAlgebra: cross; include("test/runtests_unit_simulate.jl")'
JULIA_NUM_THREADS=6 julia --project=. -t 6 test/runtests_unit_postprocess.jl

SSW_NO_PLOT=true SSWRP_TRIVIAL_ONLY=true JULIA_NUM_THREADS=6 \
  julia --project=. -t 6 examples/ssw_representation_probe.jl
sbatch examples/run_ssw_representation_study_hpc.slurm.sh
sbatch examples/run_ssw_moment_refinement_hpc.slurm.sh
sbatch examples/run_ssw_matched_geometry_hpc.slurm.sh
python3 examples/analyze_ssw_representation_study.py
```

Final focused example result is 100/100; wake is 90/90 and postprocess is
406/406. The standalone simulate suite passes all testsets with the invocation
above. `git diff --check` passes. The real Phase-A trivial gate and the full
dynamic ladder use at most six threads.

### Phase-A scaling and frozen prediction

The asymptotic rolled-up body-probe fit for `σ/d≤0.3`, using PPS=8 to suppress
lumping error, is

```text
M1 velocity max = 1.2978e-4 (σ/d)^2.359, R²=0.586
M1 gradient max = 1.7466e-3 (σ/d)^3.020, R²=0.576
M2 velocity excess vs PPS=8 = 7.0662e-6 (h_trailing/σ)^1.127, R²=0.471
M4 body velocity max = 4.1357e-3 (d/σ)^-3.391, R²=0.858
```

M1 is consistent with an approximately quadratic regularized-kernel velocity
error in its asymptotic range; the gradient carries one additional inverse
length. M2 is weak and poorly collapsed, agreeing with the exact first-moment
ledger: PPS affects higher moments, not vector circulation or impulse. M4 is
the clearest field correlation. These are diagnostic fits, not universal
orders; the rolled-up tip core and mixed trailing/spanwise segment families
limit the one-variable collapse.

Both immutable frozen tables were written before their particle ladders:
`phase_b_prediction.csv` (SigmaPPS, 168 cells) and
`phase_b_prediction_sigma_overlap.csv` (15 cells). Static closure succeeds in
18/38 observed cells but fails globally. The worst dynamic cell
(`σ/c=1.2`, PPS=1, `d/c=0.25`) is predicted at `ΔCL=-35.602%`,
`ΔΓ=37.873%`, but observes `+2.421%`, `2.753%`.

### Dynamic admissibility map

All 38 cells satisfy the final-one-chord settledness gates at `t*=20`; no
5-chord extension was needed. Across the ladder the largest CL drift/ripple
is 0.0362%/0.0362%, and the largest TE-Γ drift/ripple is
0.0425%/0.0425%, comfortably below 0.1%/0.25%.

Final measured downstream extents span 19.929–20.214c. The passive
corrected-Pedrizzetti recorder sees per-update mean relative strength changes
of 1.062–1.600% across particle cases and per-particle maxima of
39.38–43.93%; every sample is obtained through a wrapper that calls the stock
update unchanged. The complete counts and values persist in each
`case_diagnostics.csv` and in `split_summary.csv`.

The SigmaPPS cells are mechanism-isolation-only. Under the pre-registered
fixed-overlap policy, 14/15 eligible cells are admissible:

| σ/c | d/c=0.25: ΔCL / ΔΓ | d/c=1: ΔCL / ΔΓ | d/c=4: ΔCL / ΔΓ |
|---:|---:|---:|---:|
| 0.08 | +0.0458% / 0.0742% | +0.0248% / 0.0300% | +0.0153% / 0.0156% |
| 0.15 | -0.1620% / 0.1645% | -0.1358% / 0.1355% | -0.0588% / 0.0634% |
| 0.30 | -0.2556% / 0.2577% | -0.2790% / 0.2740% | -0.1140% / 0.1192% |
| 0.60 | +0.2857% / 0.6012% | -0.4070% / 0.4117% | -0.1705% / 0.1794% |
| 1.20 | **+2.4214% / 2.7528%** | +0.0239% / 0.1822% | -0.3219% / 0.3443% |

The bold cell is the only failure. For every buffer, both `|ΔCL|` and `ΔΓ`
decrease from σ/c=0.3→0.15→0.08. The conditional `h∝σ^1.5`
growing-overlap path is therefore not triggered. Fixed overlap remains the
marginal (`h/σ≈constant`), admissibility-eligible policy specified in the
design; it is not a claim of the stronger provable `h/σ→0` limit.

The plotted map and complete tables are:

- `data/ssw_sheet_particle_split/analysis/admissibility_map.png`
- `data/ssw_sheet_particle_split/analysis/admissibility_table.csv`
- `data/ssw_sheet_particle_split/analysis/scaling_and_closure.png`
- `data/ssw_sheet_particle_split/analysis/prediction_observation_residuals.csv`

### Mechanism budget and startup transient

The matched-geometry controller preserves the normal relaxed Euler step and
then replaces every pre-existing particle position by its
freestream-convected position. It leaves newly shed particles and every
strength/radius unchanged. Results:

| case | default ΔCL | matched ΔCL | frozen ΔCL | default−matched (M3b) | matched−frozen (M3a bound) |
|---|---:|---:|---:|---:|---:|
| finest admissible, σ/c=.08, d/c=4 | +0.0153% | +0.0267% | -0.0005% | -0.0114 pp | +0.0272 pp |
| worst CL/Γ, σ/c=1.2, d/c=.25 | +2.4214% | +2.4663% | -35.6016% | -0.0449 pp | +38.0679 pp |

The Γ budgets are respectively `-0.0105/+0.0260 pp` and
`-0.0315/-35.0889 pp`. Thus trajectory/roll-up divergence M3b is small. The
large nonclosing budget is strength evolution at matched geometry—M3a in the
study's bounded sense, including relaxation-coupled strength dynamics—not a
static kernel or quadrature error. M1/M4 determine instantaneous field error
and the near-body failure boundary; M2 and M3b are secondary here.

The worst startup difference through `t*=4` is `max|ΔCL|=0.3083` for the
same `σ/c=1.2,d/c=0.25` cell. The history-level table and comparison plot are
`startup_transient_metrics.csv` and `startup_transient.png`. This confirms
that the unsteady/spanwise family is most stressed in the same inadmissible
near-body, chord-scale-core regime.

### Moment target

Midpoint SigmaPPS preserves each straight subsegment's vector circulation and
`x×Γ` moment algebraically, so particle-vs-sheet first-moment error is zero
(roundoff in the explicit ledgers). The full rolled-up sheet impulse changes
0.4349% from `n_span=24→48`, exceeding the 0.25% trigger. Against the required
`n_span=96` reference, the 24- and 48-strip errors are 0.3030% and 0.2031%.
The sheet-vs-refined-sheet error is therefore more than twice the
particle-vs-sheet error by an unambiguous margin. For accuracy-oriented future
shedding (item 016), target the refined/continuous distribution's moment
rather than merely reproducing the coarse panel sheet; for representation
equivalence, retain the exact sheet ledger as a non-negotiable conservation
check.

### Transferable prescription and negative/positive boundaries

The final table reports nominal `σ/c`, `σ/d`, `d/σ`, and the actual
ceiling-derived nominal segment spacings:

- SigmaPPS trailing `h/σ=dt*/(PPS·σ/c)` and unsteady
  `h/σ=(AR/n_span)/(PPS·σ/c)`;
- SigmaOverlap uses `h=L/ceil(overlap·L/σ)` separately for
  `L=Udt` and `L=b/n_span`;
- rotor trailing `h=Δψr/PPS`, unsteady `h=Δr/PPS`, with every σ normalized by
  local chord.

**Wing prescription:** with legacy shedding and overlap 1.3, use
`σ/c≤0.6` for any tested handoff distance `d/c∈[0.25,4]`. If chord-scale
`σ/c≈1.2` is unavoidable, retain at least one chord of panel sheet before
conversion. The shorter handoff is decisively inadmissible. This is a wing
map, not a rotor validation: blade–wake re-encounter and continuously loaded
curved spanwise rows remain Phase-D caveats.

The rotor projection from the existing 41,778-particle/9-hour calibration is:

| σ/c | projected particles | projected wall time |
|---:|---:|---:|
| 0.60 | 104k | 22.5 h |
| 0.30 | 209k | 45 h |
| 0.15 | 418k | 90 h |
| 0.08 | 783k | 169 h |

Thus the wing-safe `σ/c≤0.6` bound is analytically affordable but expensive;
the robust fine end is not affordable under the 500k capacity. The
transferable next action is a deferred rotor test at σ/c≈0.6 with fixed
overlap, or sheet-until-`d/c≥1` at chord-scale σ. If neither survives rotor
re-encounter, finish item 016's smooth surface-vorticity shedding and/or item
012's overlap management. Phase D remains deferred as required.

### Technical outcome

Phases 0–C are technically complete: the admissibility map, static scaling,
dynamic mechanism attribution, moment discussion, transferable-variable map,
startup comparison, and rotor cost projection are all present. The outcome is
mixed but actionable: legacy particles are equivalent on the wing over a
broad fixed-overlap region, but the near-body chord-scale-core cell fails, and
static frozen replacement does not predict its evolved strengths. Clear-
context and user approval remain deliberately unchecked.

## 2026-07-31 — T4 (FMM |Das|-radius discriminator) EXECUTED by BRAINSTORM/018

Rotor A/B on 40_40 (NT=6, velocity, SFS off both arms, sigma_overlap 2.0/4, matched physical Das via η·36/NT): FMM-vs-DirectBackend CT delta −0.0033% at Das=0.20c and −0.0019% at Das=1.64c — η-variation of the backend delta 0.0014% ≪ 0.1% gate ⇒ **the FMM panel-radius |Das| inflation does not contaminate the rotor Das axis**. T4 closed. Data: `data/p018_t4_discriminator/`. (Found en route: FLOWVPM SFS direct-path small-N bug, `Estr_direct_multithreaded` zero-step thread range when np < nthreads.)
