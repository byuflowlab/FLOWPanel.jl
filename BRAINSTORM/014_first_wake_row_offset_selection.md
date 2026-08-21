# First-Wake-Row Offset (`Das` / η): On What Principle Should It Be Chosen?

## The question

`Das` is the offset of the first wake row from the trailing edge, set as
`Das = max(η·dt·|V_te|, min_displacement)`. **On what principle should η be
selected, such that the resulting `CT` is a *converged prediction* rather than a
*tuned* number?**

This is an open question, not a task with a known method. It is posed here
because η turned out to be the single largest discretization sensitivity found
anywhere in the rotor-hover study, and because **every mechanism proposed to
explain that sensitivity has been refuted** (see below). The empirical behaviour
is well characterised; the principle for choosing η is not.

## Current status (as of 2026-07-31)

**Partly answered, still open.** η is the *wrong* parameter — the effect collapses onto the
physical length `Das` = η·dt·U, and dt-converged lift grows only logarithmically in `Das`
(+0.205%/doubling). So the draft rule is: pick `Das` as a physical length on the log plateau
(≈0.25–1.5 local chords) and carry ±~1% model-form uncertainty. **But that does not explain
the rotor**, whose floor-free sensitivity is **+37%** — ~77× the log law — so a
rotor-specific amplifier still exists. The wake-representation-inconsistency candidate was
eliminated; the surviving lead is the **sheet/particle split at rotor σ/chord**, spun off to
item **017** and now carried by item **018**'s clearance track.

Entry points for the latest state, in order:
- `## Ruling 2026-07-30 (Ryan): N=36 rows is NOT adopted…` — the standing ruling on `nwakerows`.
- `## 2026-07-30 — Sheet/particle question spun off to item 017`.
- `## 2026-07-31 — In-flight jobs harvested (by BRAINSTORM/018 Phase 0)` — final job disposition.

No jobs remain in flight for this item.

## Why it matters

On the DJI9443 production mesh at 5400 RPM (BRAINSTORM/006, 2026-07-28), the
10-rev cycle-mean `CT` depends on η far more strongly than on anything else
tested:

| lever | effect on CT |
|---|---|
| **`Das` η: 0.2 → 1.0** | **+16.0%** (0.06148 → 0.07133) |
| mesh refinement 7,288 → 36,752 panels | < cycle scatter |
| `GreenReconstruction` vs velocity | +0.9% (< scatter) |
| truncation depth 4R → 6R | ≤ 0.22% (null) |
| relaxation filter depth 0.5R vs 1.0R | +0.6% |
| far-field relaxation strength `rlxf` | null |

η=0.2 was hard-coded in the driver from the beginning and inherited by **every**
prior rotor number. Changing it moved `CT` from ~15% below experiment into
006's 0.068–0.072 acceptance band. **If η is a free parameter, then 0.0713 is a
tuned result and the validation claim is void.** If η is determinable on
principle, 0.0713 is a prediction. Nothing currently distinguishes these.

The sensitivity is not rotor-specific: the pitching-wing oracle shows **+13% over
a 4× `das` range** (`code_audit/log.md`, Task 1b), so this is a property of the
unsteady panel-wake scheme, not of hover.

## What is established (verified in source and by measurement)

**Geometry** — the near wake is partitioned between two representations:

| element | extent | source |
|---|---|---|
| body's attached doublet panel | TE → TE+`Das` | `src/FLOWPanel_elements_fmm.jl:1044-1066` |
| `PanelWake` row 1 | at TE+`Das` | `src/FLOWPanel_simulate.jl:1021` (`update_TE!`) |
| `PanelWake` row 2 | one convection step further | — |
| particles | shed from rows ≥ 1 | `src/FLOWPanel_wake.jl:1251,1256` |

`FLOWPanel_wake.jl` contains **no** reference to `Das`. The partition **tiles
space cleanly at any η** — panel on `[TE, TE+Das]`, particles from `TE+Das`
onward, abutting exactly, with no gap and no overlap. So η is a genuine modelling
parameter: it controls *how much of the near wake is a rigid flat sheet versus
particles*, and nothing in the scheme's internal consistency selects a value.

**`Das` is frozen at its `t=0` magnitude.** `initialize_Das!` runs once at `t=0`,
where the driver's `spinup_fraction(0) = SPINUP_START_FRACTION = 0.4`, and
`rotate_Das!` only rotates the stored vector — it never rescales. So in steps of
TE travel **at operating RPM**, `η_eff = 0.4·η`:

| nominal η | η_eff | `Das` at 0.75R | local chords | `Das`/σ |
|---:|---:|---:|---:|---:|
| 0.2 | 0.08 | 1.25 mm | 0.08 | 0.053 |
| 0.5 | 0.20 | 3.12 mm | 0.20 | 0.133 |
| 1.0 | 0.40 | 6.23 mm | 0.41 | 0.267 |
| 4.0 | 1.60 | 24.9 mm | 1.63 | 1.067 |

**Whether freezing `Das` at the spin-up value is intended or incidental is itself
open** — it means the offset does not track the operating condition it was
derived from.

**Empirical CT-vs-η (production mesh, 5400 RPM, NT=36, 4R, unfiltered,
10-rev cycle-means):** 0.06148 (η=0.2) → 0.07133 (1.0) → 0.07190 (4.0). Flat for
η ≥ 1 (+0.8%, inside cycle scatter). η=0.5 pending (job 12920967); its transient
suggests ~77% of the jump is already realised there, so the knee is below 0.5.

**Particle core size, measured** (settled step 700, 41,778 particles): σ tracks
`σ/r = Δψ·OVERLAP/PPS` to ~10% across the span. At 0.75R, **σ ≈ 1.5 local
chords** — newly shed particles engulf the blade section that produced them.
`Das/σ = 0.4·η/(OVERLAP/PPS) = 0.267·η`, exactly radius-independent.

## Mechanisms proposed and REFUTED — do not re-tread

1. **"Larger `Das` pushes the near wake out of influence, so CT asymptotes."**
   *Refuted:* the body carries an attached doublet panel on `[TE, TE+Das]` at all
   times. Nothing is removed from the near field; the region is always
   represented.
2. **"η<1 seeds particles inside the attached panel's span (double-counting)."**
   *Refuted:* `update_TE!` places `PanelWake` row 1 exactly at TE+`Das`, and
   particles are shed from rows ≥ 1. They abut the panel precisely. No overlap at
   any η, so there is no inconsistency to fix.
3. **"`Das` < σ makes the first particle's core straddle the TE."**
   *Refuted by measurement:* if this were controlling, CT should improve until
   `Das/σ ≳ 1` (η ≈ 4). CT is already flat from η=1, where the core still exceeds
   the offset **3.7×**.

**Net: the +16% is currently an empirical fact with no surviving explanation.**

## What would settle it — candidate directions

- **Temporal convergence at fixed η** (in flight: job 12921559, NT=72, η=1.0).
  If CT is dt-invariant at η=1, η is at least not a dt-artifact. *Caveat:* at
  fixed η, halving dt halves `Das` **and** σ as well, so three things refine
  together; `Das/σ` is invariant under it. A null result is strong; a non-null
  result cannot be attributed without a follow-up that varies them separately.
- **Separate the partition from the resolution.** Vary `Das` at fixed σ and fixed
  dt (and vice versa) to establish which of the two the +16% actually responds
  to. Currently they are confounded through `dt` in both `Das` and σ.
- **Ask what the attached panel is *for*.** It is rigid, flat and
  non-convecting over its length, so on physical grounds it should be as short as
  possible — the opposite of what CT prefers. That tension is unexplained and may
  be the most informative thread.
- **Un-freeze `Das`.** Re-derive it from the *current* `|V_te|` each step rather
  than from `t=0`, and see whether the η sensitivity persists. This also removes
  the 0.4 spin-up factor as a hidden variable.
- **Verification against a case with a known answer.** The pitching-wing oracle
  (`code_audit/log.md`) has reference data and shows the same sensitivity; if a
  principle for η is proposed, it must produce the right answer *there* too. This
  is the strongest available discriminator and does not require rotor runs.
- **Analytic rung.** Wagner / Theodorsen (BRAINSTORM/013 item 3) have exact
  shed-wake solutions; a 2D or sudden-start test would expose whether a defensible
  η follows from matching the analytic near-wake circulation distribution.

## Acceptance

014 is answered when there is a **stated rule for η that is justified other than
by the CT it produces**, that is verified on at least one case with an
independent reference answer (oracle wing or analytic), and that is either
dt-invariant or explicitly dt-dependent in a derived way. A negative outcome is
also acceptable and valuable: a demonstration that η is *irreducibly* a model
parameter, in which case its ±16% must be carried as a **model-form uncertainty**
on every `CT` this project reports — including 006's 0.0713.

## Cross-references

- `BRAINSTORM/006_rotor_hover_converge_stable_near_reference.md` — where the
  sensitivity was found; 2026-07-28 sections.
- `logs/dji_convergence_20260722/phase_02e_unsteady_ct.md` — full run tables,
  measurements, and the dated retractions of mechanisms 1–3.
- `BRAINSTORM/012_robust_resolution_overlap_management.md` — σ ≈ 1.5 chords at
  shedding is a sharper statement of the under-resolution 012 owns.
- `BRAINSTORM/013_unsteady_wake_validation_cases.md` — items 3 (Wagner/
  Theodorsen) and 1/2 (heaving/pitching foils) are the verification rungs.
- `code_audit/log.md` Task 1b — pitching-wing oracle, +13% over a 4× `das` range.
- Arc-vs-tangent `Das` construction (`DAS_KINEMATIC_ARC`, default true) is a
  *separate* question — that is about the offset's direction, not its length.
  A/B in flight as job 12921071.

## Added 2026-07-28 (Ryan): decouple `Das` from the panel-wake buffer, and question whether `Das` is needed at all

### Proposal 1 — smaller `Das`, more wake rows. Do they converge to the same answer?

There is already a **free doublet-sheet wake row** between the `Das` row and the
particles. With `nwakerows = N` (once `nwakes[]` saturates): row 1 sits at
TE+`Das` and is **rigid** — re-placed every step by `update_TE!`, it never
convects — while rows 2…N+1 are **free** and convect with the flow, and
particles are shed from rows N → N+1 (`_convert_to_particles!` uses
`nodes[:, nwakes, ·]` and `nodes[:, nwakes+1, ·]`).

So the particle handoff sits at roughly `Das + (N−1)·(one step of TE travel)`.

**This is the crux of the whole item.** At the `N = 1` used everywhere to date,
*the length of the rigid attached sheet* and *the distance to the particle
handoff* are *the same number*. They cannot be told apart, which is precisely why
the +16% has resisted explanation. Raising N separates them:

- If CT depends only on the **total** handoff distance `Das + (N−1)·travel`, then
  `Das` is merely a way to buy panel-wake extent — and extra **free** rows are
  the better way to buy it, since they convect while the `Das` row cannot. η then
  has no special status and should be set as small as the Kutta enforcement
  permits, with N carrying the extent.
- If CT depends specifically on **`Das`** at fixed total handoff distance, the
  effect is about the rigid attached sheet / Kutta geometry, and η is a genuine
  model parameter.

*Concrete test:* hold `Das + (N−1)·travel` fixed and trade one against the other
— e.g. (η=1.0, N=1) vs (η_eff→small, N=2) vs (η small, N=4) — and separately vary
N at fixed `Das`. This is the **decisive experiment for 014** and supersedes the
vaguer "separate the partition from the resolution" bullet above.

*Implementation status:* `PanelWake`/`PanelParticleWake` already accept
`nwakerows` (default 100 in the constructor; a commented-out `nwakerows=12`
sits at `examples/rotor_hover_pressure_comparison.jl:322`). But the driver
**hard-codes `nwakerows=1`** at line 335 — a third never-varied hard-coded
discretization choice, alongside η=0.2 and the `min_displacement` floor below.
It needs an ENV hook, and line 764 writes a **literal** `nwakerows = 1` into the
case metadata, which would silently lie once the value is varied. Fix both
before running the sweep.

### Proposal 2 — reframe: `Das` exists to enforce the Kutta condition

`Das` is not a physical modelling choice in origin. **It exists so the shed wake
has a finite geometry at the trailing edge against which the Kutta condition can
be enforced** — the attached doublet panel spanning TE → TE+`Das` carries the TE
circulation (`get_strength_doublet` + `wake_strength_shift`,
`src/FLOWPanel_elements_fmm.jl:1035-1038`). As `Das → 0` that panel degenerates,
so there is a *numerical* floor on `Das` that has nothing to do with the physics
of the wake.

That suggests the real question may not be "what η?" but **"is there a Kutta
enforcement that does not require a finite first-row offset at all?"** — e.g.
iterative pressure-equality Kutta, or a formulation that matches Γ_TE directly
without a finite attached panel. If such a scheme converges without `Das`, the
+16% sensitivity dissolves rather than being resolved, and the tuned-vs-predicted
worry disappears with it.

The repo already carries alternative formulation machinery
(`src/FLOWPanel_formulation.jl`: `GreenReconstruction`, `TraceCorrected`) and a
`correct_kuttacondition` option, so this is not a from-scratch direction.
Item 001 measured Kutta-condition sensitivity as negligible (Bernoulli Δ0.4%) —
but that tested the *correction*, not the *existence of a finite `Das`*, so it
does not close this off.

### New confound found while checking the above: the `min_displacement` floor

`_accumulate_Das!` (`src/FLOWPanel_simulate.jl:10`) applies
`displacement_length = max(|η·dt|·|V_te|, min_displacement)`, with
`min_displacement = DAS_MIN_DISPLACEMENT_R · R` defaulting to **0.01R** and
**never overridden by the launcher**. The kinematic term scales as `r`; the floor
does not. With `Das_kin/r = 0.4·η·Δψ = 0.0698·η` and floor `0.01R`:

| η | floor active below r/R | % of shedding span (0.111–1.0) clamped |
|---:|---:|---:|
| **0.2** | **0.716** | **68.1%** |
| 0.5 | 0.286 | 19.7% |
| 1.0 | 0.143 | 3.6% |
| 4.0 | 0.036 | 0% |

**At η=0.2 — the value hard-coded from the beginning and used for every prior
rotor number — `Das` is set by the floor, not by η, over ~68% of the blade**, and
is radially *constant* (0.01R) there instead of proportional to `r`. The η ladder
is therefore **not a clean single-parameter sweep**: 0.2 → 1.0 simultaneously
scales `Das` and releases the floor, changing its *spanwise distribution*.

This is a code-level fact. Whether it *explains* the +16% is **untested**.
Unlike the three refuted mechanisms above, it is a concrete difference in what
was actually computed, and it survives inspection so far.

### RECOMMENDED NEXT RUN — the floor test (highest priority in this item)

**Run this before any further η rungs, and before the `nwakerows` sweep.** It is
one job, it costs the same as any other Phase 2e case (~10 h), and it can
reframe 006's headline result.

*Configuration:* the completed unfiltered η=0.2 baseline, changing **one** thing
— the floor:

| setting | value |
|---|---|
| `DAS_ETA_KINEMATIC` | 0.2 (unchanged) |
| **`DAS_MIN_DISPLACEMENT_R`** | **1e-6** (default 0.01; floor effectively off) |
| everything else | as `p2e_vel_nt36_d4`: production mesh, 5400 RPM, velocity, NT=36, truncation 4R, `rlxf` 0.3, no relaxation filter, `SETTLE_REVS=12` (719 steps) |

*Comparison:* the 10-rev cycle-mean over revs 10–19 against the two existing
points — η=0.2 with floor **0.06148 ± 0.00152**, η=1.0 **0.07133 ± 0.00159**.

*Interpretation, fixed in advance:*
- **CT rises toward 0.0713** ⇒ the "η sensitivity" is substantially a **floor**
  sensitivity. The real parameter is `min_displacement`, not η, and the whole η
  ladder needs re-reading — including 006's claim that η is the magnitude lever.
- **CT stays near 0.0615** ⇒ the floor is irrelevant, η itself is the lever, and
  the η ladder stands as measured.
- **CT lands between** ⇒ both contribute, and the ladder must be repeated with
  the floor disabled throughout before any η can be called converged.

*Optional converse, only if the first result is ambiguous:* η=1.0 with a
deliberately **larger** floor, which should then reproduce the low CT.

*Note for whoever runs it:* this is recorded here as a recommendation only. It
has **no launcher case tag and no entry in the Phase 2e job ledger**
(`logs/dji_convergence_20260722/phase_02e_unsteady_ct.md`), so it will not be
picked up by an agent working from that log — it must be staged deliberately.

## Added 2026-07-28 (agent): theory rung — what η *should* be, derived from this scheme's discrete geometry

The discrete near wake is a piecewise-constant doublet distribution: the
attached panel on [TE, TE+`Das`] carries the *current* bound circulation Γ(t)
(`get_strength_doublet`), and `PanelWake` row k carries Γ(t − k·dt), starting at
TE+`Das`. Equivalently (constant-doublet ≡ edge vortex filaments), the junction
between the attached panel and row 1 carries a concentrated vortex of strength
ΔΓ = Γ(t) − Γ(t−dt) — the circulation shed *this step* — located at η·V·dt
behind the TE.

The continuum it approximates is the retarded-time doublet sheet
μ(x) = Γ(t − x/V): the shed increment ΔΓ is *distributed* over [TE, TE+V·dt]
(uniformly, for Γ̇ ≈ const over the step). The discrete scheme replaces that
smooth ramp with a step at x = η·V·dt. Three standard matching rules give three
constants:

| rule | η | order of the moment error |
|---|---:|---|
| left endpoint (jump at TE) | → 0 | first |
| **centroid / midpoint** (jump at the shed sheet's centroid) | **1/2** | **second** |
| right endpoint / retarded-position match (row k at x = k·V·dt) | 1 | first |

Katz & Plotkin's familiar "0.2–0.3 of the TE travel" arises in a *different*
convention (their newest wake panel spans [TE, c·V·dt] and carries the new wake
unknown), so it does not transfer to this scheme; the analogous derivation here
selects **η = 1/2** as the accuracy-optimal constant, with **any fixed η
first-order consistent** (`Das = η·V·dt → 0` as dt → 0, so the scheme converges
to the same continuum limit for every η, including η = 4).

**Predictions this makes (testable on the Wagner/Jones oracle, minutes):**

1. **P1 (consistency):** at fixed η, the computed indicial curve converges as
   dt → 0, and different η converge to the *same* limit — the η-spread across
   the matrix must → 0 with dt. If it does not, the scheme has an
   η-inconsistency (Branch B) and the Kutta-without-`Das` reframing takes over.
2. **P2 (superconvergence):** finite-dt error is smallest near η ≈ 1/2 and
   grows roughly linearly in |η − 1/2|·dt; large η (2, 4) shows a downstream
   *lag* signature (the whole wake sheet sits (η−1)·V·dt too far back, since
   row 1 is pinned at TE+`Das` and convection preserves the offset).
3. **P3 (rotor contrast):** a pure truncation error cannot *saturate* in η at
   fixed dt. The rotor CT's flatness for η ≥ 1 (η_eff ≥ 0.4) therefore does NOT
   look like this mechanism — it looks like a near-field threshold (floor
   release, or particle-handoff proximity), which is exactly what the staged
   runs (`p2e_das0p2_nofloor`, `p2e_nrows2_dassmall`, `p2e_nrows4_das1p0`,
   `p2e_das1p0_refresh`) separate.

Corollary if P1–P2 hold on the oracle: the **stated rule** 014's acceptance asks
for is *"η = 1/2 on accuracy grounds (centroid placement), any η ≥ small floor
admissible with dt-refinement, `Das` re-derived from the current TE velocity
each step (`DAS_REFRESH`), floor ≈ 0"* — justified independently of the CT it
produces — and the rotor's η ladder must then be re-read as a near-field
artifact question (floor / handoff), not an η-selection question.

## RESULTS 2026-07-29 (re-run on the corrected mesh): finite-AR Wagner oracle η×dt matrix — η is the wrong parameter; the physical offset length governs, logarithmically

Oracle-first program on `examples/suddenly_started_wing.jl` at finite AR (AR=100
is broken — its steady semi-infinite solve returns `CLsteady` ~1e5; at AR=6 the
steady solve is sane and R.T. Jones' AR=6 elliptic indicial function, NACA Rept.
681, is the external band). Matrix: η ∈ {0.125…4} × dt* ∈ {1/4…1/32},
rectangular NACA0012, n_span=24, n_airfoil=21, `DirectBackend`, pure panel wake
(no particles, no floor, no frozen `Das`), metrics masked to t* ≥ 1
(impulsive-start ring). Machinery `examples/ssw_eta_convergence.jl`; data in
`data/ssw_eta_convergence/` (`summary.csv`, `analysis.csv`, `das_analysis.csv`,
plots). `CLsteady = 0.38005447`; Jones tail at t*=7 is 0.998258.

These numbers **supersede an earlier run of the same matrix** that was computed
before the two topological mesh defects (see the wake-representation section
below) were fixed: the triangulation was invariant under neither reflection, and
the spanwise defect put ~10% asymmetry into shed circulation. Every finding below
was re-derived from scratch on the corrected mesh (`SSW_MESH_REVISION = 2`, now
carried in the case tag so this cannot silently recur). **All four findings
stand**; the numbers moved by <1% and one loose statement of F2 is corrected.

Tail CL/CL_steady (t* ∈ [6,7]):

| η\dt* | 1/4 | 1/8 | 1/16 | 1/32 |
|---:|---:|---:|---:|---:|
| 0.125 | 0.965144 | 0.964846 | 0.962248 | 0.958938 |
| 0.25 | 0.976255 | 0.975032 | 0.972862 | 0.970404 |
| 0.5 | 0.982375 | 0.980879 | 0.978903 | 0.976629 |
| 1 | 0.986208 | 0.984551 | 0.982711 | 0.980624 |
| 2 | 0.989089 | 0.987307 | 0.985456 | 0.983509 |
| 4 | 0.991515 | 0.989688 | 0.987782 | 0.985857 |

**F1 — the collapse variable is the physical offset `Das` = η·dt·U, not η.**
Along fixed-`Das` diagonals dt-refinement converges cleanly: successive-increment
ratios **0.39–0.47** across all five diagonals, giving Richardson limits

| `Das` | c/64 | c/32 | c/16 | c/8 | c/4 |
|---|---:|---:|---:|---:|---:|
| dt→0 limit | 0.97997 | 0.98200 | 0.98412 | 0.98618 | 0.98815 |

Along fixed-η rows there is **no dt→0 limit**: each halving subtracts a
near-constant amount and the increments do *not* contract (successive ratios
1.01–1.13 at η ≥ 1, i.e. ≈1). The sharpest form of this — and the direct
refutation of the theory rung's P1 — is that **the η-spread in the tail is
dt-invariant**: 0.991515 − 0.965144 = **2.637%** at dt*=1/4 versus
0.985857 − 0.958938 = **2.692%** at dt*=1/32. If η were ordinary truncation
error that spread would collapse under dt-refinement. It does not move.
**η is not a convergent discretization parameter of this scheme.**

**F2 — lift grows logarithmically in `Das`, +0.205%/doubling.** Fitting the
dt-converged (Richardson) limits against log₂`Das` gives slope
**+0.2053%/doubling with RMS residual 0.0035%** over c/64–c/4 — an excellent log
fit (`das_analysis.csv`, `log_law_richardson`). No saturation is seen at the
large-`Das` end: probes at fixed dt*=1/4 continue to rise through `Das` = 2c
(0.993749) and 4c (0.995893), at +0.223% and +0.214% per doubling.

*Correction to the previous statement of F2.* It quoted a single slope
"+0.23%/doubling, unsaturated 0.004c→4c" read off a fixed-dt row. That conflates
two dependencies: along a row of constant dt, η varies too, so the dt-error
varies with it. The per-doubling increment along the dt*=1/4 row actually runs
**1.11% → 0.61% → 0.38% → 0.29% → 0.24% → 0.223% → 0.214%** as `Das` goes
0.03c → 4c, and a naive log fit to that row has RMS residual 0.303% against a
0.396% slope — nearly meaningless. The logarithmic law is a statement about the
**dt-converged** limit, where it is clean; the steepening at small `Das` on a
fixed-dt row is F1's dt-error reappearing, not a steeper log law.

**F3 — mechanism reading (unchanged, and now independently corroborated).** The
rigid span [TE, TE+`Das`] (+ pinned row 1) holds the near wake flat; free rows
droop with the local downwash. Larger `Das` ⇒ flatter effective near wake ⇒ more
lift, approaching the flat-/semi-infinite-wake normalization logarithmically.
Neither limit is privileged inside the scheme. The wake-representation study
(below) supports this directly and independently: at converged wake length a
**straight** (freestream-convected) sheet gives **+0.96%** more lift than the
rolled-up wake — i.e. flattening the wake does raise CL, by an amount comparable
to several `Das` doublings. Jones (0.998258 at t*=7) sits above all finite-`Das`
values, mildly favouring larger `Das`, but elliptic-vs-rectangular and the t*=7
window prevent using it to pick the constant.

**F4 — the rotor's +37% is NOT the intrinsic log law.** With `η_eff = 0.4η`, the
rotor's η 0.2 → 1.0 moves `Das` at 0.75R from 0.08 to 0.41 local chords —
log₂(5.125) = **2.36 doublings**, so the log law predicts **+0.48%**. The
measured floor-free rotor change is **+36.8%** (0.05215 → 0.07133): about **77×**
the oracle slope. The rotor's flat top does match, though — η 1.0 → 4.0 is 2
doublings, predicted +0.41% against a measured +0.8% (inside cycle scatter). So
the ladder's plateau is the log law, and the collapse below η ≈ 0.5 needs a
rotor-specific amplifier: the `min_displacement` floor, the `nwakerows=1`
particle handoff, or the σ ≈ 1.5-chord overlap.

**Mesh sensitivity at the check corner** (η=0.25, dt*=1/8): n_span 24→48 changes
the curve by 0.435% max / 0.280% L2; n_airfoil 21→41 by 2.085% max / 0.874% L2
(`CLsteady` 0.38005 → 0.37918, −0.23%). Trends are resolved; absolute levels are
good to roughly ±1%, not better. Note n_span=24 is *not* itself mesh-converged in
absolute CL (steady 0.42122 / 0.38005 / 0.36491 at n_span 12/24/48) — but every
matrix entry is normalized by its own `CLsteady`, so the differential findings
are insulated from that.

**F5 — pre-registered prediction for NT=72 η=1.0 (job 12921559), restated.**
Frozen `Das` halves in physical units, i.e. one doubling down the log law, so CT
should drop only ≈0.2% (to ~0.0712). Materially unchanged by the re-run: the
slope moved +0.23% → +0.205%/doubling. A fall toward 0.0615 instead means the
rotor's amplified regime starts near `Das` ≈ 0.2 local chords, and 006's
magnitude claim needs the amplifier resolved first.

### Revised answer to this item's question (draft rule, pending rotor confirmation)

1. **Parameterize the first-row offset as a physical length** (fraction of local
   chord / TE travel at a fixed reference step), not as η per current step. The
   accidental freezing of `Das` at t=0 made the rotor ladder a fixed-length
   family — which is why it "converged"; `DAS_REFRESH` re-ties `Das` to dt and
   must NOT be used inside dt-refinement studies (scale η ∝ NT instead, as 006
   study D already prescribes).
2. **Choose the offset on the log plateau above the scheme's near-field scales**
   (wake core, panel size, particle σ) — in practice ≈0.25–1.5 local chords. The
   frozen η=1.0 rotor setting (≈0.4 chords at 0.75R) qualifies.
3. **Carry the log slope as the honest model-form uncertainty:** ≈0.2%/doubling
   in CL within the plateau — i.e. well under ±1% across any defensible choice,
   NOT ±37%. The ±37% belongs to the sub-plateau amplified regime, which is a
   configuration error to avoid (and whose rotor mechanism R1/R2 will name), not
   an uncertainty to carry.
4. Open for full acceptance: rotor attribution (R1/R2/R3 + 12921559 vs the F5
   prediction).

### Addendum — long-window quasi-steady run

t*=20, n_span=48, dt*=1/16 (`data/ssw_eta_convergence/long_fine/`), tails over
t* ∈ [19,20] against the Jones asymptote 1.000000:

| `Das` | c/64 | c/16 | c/4 |
|---|---:|---:|---:|
| tail | 0.992681 | 0.996590 | 0.997892 |
| deficit vs Jones | 0.73% | 0.34% | 0.21% |

Two refinements to F2/F3. (i) Most of the t*=7 deficit against Jones was **window
truncation, not scheme error** — settled, the scheme reaches within 0.2–0.7% of
the analytic asymptote, larger `Das` closest. (ii) The **quasi-steady** `Das`
sensitivity is weaker *and saturating*: ≈0.195%/doubling from c/64→c/16 but only
≈0.065%/doubling from c/16→c/4 (0.130%/doubling averaged over the four
doublings), against ≈0.205%/doubling in the t*=7 transient. So part of the
transient slope is a `Das`-dependent *rate* of approach rather than a difference
in the settled state, and in the settled state the plateau flattens further.
Plateau model-form uncertainty in the settled state is therefore well under 1%
across any defensible `Das`. The rotor's +37% remains un-log-like (F4 unchanged).

**Known gap — wake-core robustness NOT re-run.** F2's core probes (1e-4c / 1e-2c)
were not repeated, and cannot be as the code stands: `wake_core_over_c` has no
ENV hook and is **absent from `_ssw_case_tag`**, so two runs differing only in
core size share a tag and the second silently resumes the first. That is the same
silent-resume defect as the missing mesh revision. So the earlier finding that
the log law is not a core artifact (core 1e-3c → 1e-4c changed the tail <0.02% at
`Das` = c/32 and c/256; core 1e-2c ≈ `Das` lowered it −0.5%) **still rests on the
pre-fix mesh and is unconfirmed.** Closing it needs an `SSW_WAKE_CORE` env hook
plus a non-default core token in the tag first; the runs themselves are minutes.
See `data/ssw_eta_convergence/README.md`, which also records which case
directories in that tree are current (`*_m2`) versus dead (31 pre-fix dirs).

## Results (2026-07-28, evening): floor test and `nwakerows` test both report — first non-refuted signal

Ryan ran both recommended experiments directly. **Transient values (revs 5–10,
`SETTLE_REVS=12`)**; settled cycle-means pending. Full tables in
`logs/dji_convergence_20260722/phase_02e_unsteady_ct.md` (21:20 entry).

### 1. The `min_displacement` floor confound is RESOLVED — η attribution survives and grows

`p2e_das0p2_nofloor` (η=0.2, `das_min_R=1e-6`) vs the η=0.2 baseline, matched revs:

| rev | floor ON (0.01R) | floor OFF |
|---:|---:|---:|
| 5 | 0.06952 | 0.05786 |
| 10 | 0.06284 | 0.05411 |

Removing the floor **LOWERS** CT ~14% — the opposite of the failure mode this
test was designed to catch. The floor was **propping the η=0.2 end up**, so it
*masked* part of the η sensitivity rather than manufacturing it.

- The "the +16% is really a floor artifact" branch is **dead**.
- Floor-free, the ladder runs ~0.053 → 0.0713: the true η sensitivity is
  **~+30%, not +16%**. Every uncertainty statement in this item should use the
  larger figure.
- The clean (floor-free) η ladder is now the one worth completing.

### 2. `nwakerows=2` with tiny `Das`: CT collapses — free rows do NOT substitute for `Das`

`p2e_nrows2_dassmall` (η=0.1, floor off, **nwakerows=2**) gives **0.0366 / 0.0350
/ 0.0332** at revs 5 / 7 / 10 — about **half** the baseline, and less than half
the η=1.0 value.

The decisive comparison is against total handoff distance:

| case | rigid `Das` | free rows | total distance to particles | CT (rev 10) |
|---|---:|---:|---:|---:|
| η=1.0, nrows=1 | 0.40 steps | 0 | **0.40 steps** | 0.0737 |
| η=0.1, nrows=2 | 0.04 steps | 1 | **~1.04 steps** | 0.0332 |

The second has the **larger** total handoff distance and by far the **lower** CT.

⇒ **CT does not depend on the total distance to the particle handoff.** A free
convecting row is not interchangeable with a longer rigid `Das`. Proposal 1's
first branch ("`Das` is merely a way to buy panel-wake extent, and free rows are
the better way to buy it") is **refuted**. Its second branch — that CT depends
specifically on `Das` — is what the data support.

### 3. Where this leaves the item

The floor test and the `nwakerows` test are mutually consistent: **CT collapses
as `Das` shrinks, regardless of what is downstream of it.** After three refuted
mechanisms this is the first coherent surviving signal, and it moves the centre
of gravity of 014 from Proposal 1 to **Proposal 2**:

> `Das` exists to give the Kutta condition a finite wake geometry. As `Das → 0`
> the attached doublet panel degenerates. The η sensitivity may therefore be a
> **Kutta-enforcement** sensitivity, not a wake-representation one.

That reframes the acceptance criterion: a rule for η may not exist *because η is
the wrong question* — the productive direction is a Kutta enforcement whose
answer does not depend on a finite first-row offset.

**Recommended next steps, revised:**
1. Let the two runs settle and take 10-rev cycle-means before quoting numbers.
2. **Diagnose the circulation directly, not just CT** — compare Γ_TE (and the
   spanwise Γ distribution) across η at fixed everything else. If Γ_TE itself
   scales with `Das`, the Kutta condition is not being enforced consistently and
   that is the whole story. The `monitor03_bound_circulation` output already
   recorded per-run makes this nearly free — **do this before any more CT runs.**
3. Only then revisit whether a different Kutta enforcement removes the
   dependence.

**Note on the arc construction (separate question, now closed):** the arc-vs-
tangent A/B at η=1.0 agrees within 0.5% at all matched revs, so the offset's
*direction* is immaterial here. Only its *length* matters.

## Γ_TE diagnostic (2026-07-28): the η sensitivity acts through the BOUND CIRCULATION

Read-only analysis of `monitor03_bound_circulation` already written by each run;
no new compute. Settled window revs 10–19 (steps ≥ 360), span-averaged |Γ_TE|:

| η | mean Γ_TE | Γ ratio | CT | CT ratio |
|---:|---:|---:|---:|---:|
| 0.2 | 0.18540 | 1.000 | 0.06148 | 1.000 |
| 0.5 | 0.21138 | 1.140 | 0.06942 | 1.129 |
| 1.0 | 0.22193 | 1.197 | 0.07133 | 1.160 |
| 4.0 | 0.22550 | 1.216 | 0.07190 | 1.169 |

**Γ_TE tracks CT in near-lockstep.** The blade's own circulation changes with
`Das`; this is not the wake acting differently downstream on an unchanged
loading.

Spanwise breakdown — **the effect is concentrated INBOARD**:

| η | Γ(0.25R) | Γ(0.50R) | Γ(0.75R) | Γ(0.90R) |
|---:|---:|---:|---:|---:|
| 0.2 | 0.13185 | 0.24126 | 0.26344 | 0.15863 |
| 0.5 | 0.14760 | 0.28359 | 0.29601 | 0.17180 |
| 1.0 | 0.17894 | 0.29082 | 0.29730 | 0.17296 |
| 4.0 | 0.19576 | 0.29031 | 0.29081 | 0.16917 |

At 0.75R, Γ saturates by η=0.5 and slightly *declines* at η=4. At 0.25R it rises
monotonically, **+48%** from η=0.2 to η=4. Inboard is exactly where `Das/chord`
is smallest (0.03 chord at 0.25R vs 0.20 at 0.75R for η=0.2) — i.e. where the
attached panel is closest to degenerate. Proximity, not wake extent, is what
correlates with the effect.

**What this establishes, and what it does NOT.** It rules out any explanation in
which the wake acts downstream while blade loading stays put. It does **not**
uniquely establish the Kutta hypothesis: a change in wake-induced downwash at the
blade would also change effective incidence → Γ → CT. Kutta-enforcement and
induced-velocity routes are still confounded.

### RECOMMENDED NEXT — steady rigid-wake `Das` sweep (cheap, decisive)

Supersedes the CT-run recommendations above in priority. In a **steady** solve
with a rigid wake there is no wake evolution, so varying `Das` changes **only the
attached panel geometry**:

- **Γ still moves ~20% ⇒ the effect is Kutta/geometry**, unsteady wake dynamics
  are irrelevant to it, and 014 collapses to a Kutta-enforcement question
  (Proposal 2) that can be pursued entirely in cheap steady solves.
- **Γ is insensitive ⇒ the route is wake-induced velocity**, Proposal 2 is wrong,
  and the question returns to near-wake representation.

Cost is seconds per case, not ~10 h — item 001 already has the steady
machinery (`examples/rotor_hover_force_method_audit.jl`), and the same sweep can
be run on the pitching-wing oracle where an independent reference answer exists.
**Do this before any further unsteady η runs.**

## RECOMMENDED: steady prescribed-wake `Das` sweep — the cheap discriminator (with Ryan's caveat, verified)

### The question it answers

The Γ_TE diagnostic showed the η effect acts through the **bound circulation**,
but left two routes confounded:

- **(a) Kutta / paneling** — the attached panel's *length* changes how the Kutta
  condition is enforced, so Γ changes for purely numerical reasons.
- **(b) Induced velocity** — a different wake changes the downwash at the blade,
  changing effective incidence and hence Γ for a *physical* reason.

Both produce "Γ moves with `Das`", so the unsteady runs cannot separate them.

### How the test answers it

In a **steady solve with a prescribed (rigid) wake**, there is no wake evolution
and no particle handoff. If the wake *sheet geometry* is held fixed and `Das`
only sets **where the first panel boundary falls along that same sheet**, then
`Das` is a **pure paneling parameter** — the continuous solution cannot depend on
it. So:

- **Γ moves with `Das` ⇒ route (a).** The dependence is a discretization/Kutta
  artifact, its size measures how far the near-wake paneling is from converged,
  and 014 becomes a Kutta-enforcement question answerable entirely in **cheap
  steady solves** (seconds per case, no 10-h jobs, no particles).
- **Γ is insensitive ⇒ route (b).** The unsteady sensitivity then lives in the
  panel→particle transition rather than in the Kutta condition, Proposal 2 is
  wrong, and the question returns to near-wake representation.

This is the strongest available discriminator per unit cost, and it removes
particles, wake convection and the limit cycle from the problem in one step.

### Ryan's caveat — VERIFIED, and stronger than stated

**The obvious implementation does not work.** With `semiinfinite_wake=true`,
`Das` is not merely magnitude-insensitive — `_phi_semiinfinite`
(`src/FLOWPanel_elements_fmm.jl:1226`) **throws**:

```julia
if abs(d1*d1 + d2*d2 + d3*d3 - 1) > 2*eps()
    error("Found non-unitary semi-infinite direction" ...)
```

so any `|Das| ≠ 1` is a hard error. A `Das` sweep in semi-infinite mode is
impossible by construction.

**This is itself a finding for 014.** The struct comment reads
`Das::Vector{Array{TF,2}}  # Unitary direction of rigid wake (vertex-based)`
(`src/FLOWPanel_liftingbody.jl:68`) — **`Das` was designed as a unit direction
vector**, and its use as a finite offset *length* (`η·dt·|V_te|`) is a later
overload of the same field by the unsteady panel-wake path. One field now
carries two incompatible meanings: a *direction* in the semi-infinite steady
mode, and a *direction + length* in the finite unsteady mode. **That may be the
root of why no principle for η exists — the length was never a modelled quantity
in the original design, it arrived as a side effect.** This should be checked
against the git history of `Das` before drawing conclusions.

### Implementation — Ryan's workaround

Use `semiinfinite_wake=false` (so `|Das|` is live) and approximate the far wake
with **prescribed rigid panel-wake rows extending far enough downstream** to
stand in for the semi-infinite sheet. Requirements:

1. **Hold the far wake fixed while varying `Das`.** `update_TE!` places wake row
   1 at TE+`Das`, so naively changing `Das` translates the whole wake. Either
   re-place the downstream rows so the far-field geometry is identical across
   cases, or make the wake long enough (several R) that a shift of the
   `Das` range (≈6–25 mm) is negligible in the far field. State which was done.
2. **Converge the wake length** first — sweep total wake extent until Γ stops
   moving, so the `Das` result is not contaminated by a truncated wake.
3. **Report Γ, not just CT** — Γ_TE spanwise, since the unsteady effect was
   concentrated inboard where `Das`/chord is smallest.

### Cleanest variant — do this one first

Run it on a **steady, non-rotating wing at fixed incidence with a straight rigid
wake**. There, `Das` is *unambiguously* a pure paneling parameter: the wake sheet
is the same straight sheet regardless of where the first panel boundary sits, so
the exact solution is provably `Das`-independent and **any** Γ variation is
numerical. That gives a clean, interpretable, near-zero-cost answer before
touching the rotor's helical wake, and it can be checked against the
pitching-wing oracle (`code_audit/log.md`) where an independent reference exists.

**Priority: run this before any further unsteady η runs.**

## Settled results (2026-07-29): floor-free η sensitivity is +37%; `nwakerows` negative confirmed

Both experiments `COMPLETED`. 10-rev cycle-means, revs 10–19:

| case | `Das` (steps of TE travel) | rows | CT | std | rel | p-p |
|---|---:|---:|---:|---:|---:|---:|
| baseline η=0.2, floor ON | 0.08 (floor-clamped inboard) | 1 | 0.06148 | 0.00152 | 2.47% | 7.79% |
| **η=0.2, floor OFF** | 0.08 | 1 | **0.05215** | 0.00078 | 1.49% | 5.10% |
| η=1.0 | 0.40 | 1 | 0.07133 | 0.00159 | 2.23% | 7.15% |
| η=4.0 | 1.60 | 1 | 0.07190 | 0.00101 | 1.40% | 4.19% |
| **η=0.1, nrows=2, floor OFF** | 0.04 + 1 free row ≈ 1.04 total | 2 | **0.03362** | 0.00046 | 1.35% | 3.71% |

**1. The true η sensitivity is +36.8%, not +16%.** Floor-free, η=0.2 settles at
**0.05215 ± 0.00078** vs η=1.0's 0.07133. The `min_displacement = 0.01R` floor
was inflating the low end by **15.2%** and thereby masking **more than half** the
real dependence. Every uncertainty statement in this item and in 006 should use
**±37%**, not ±16%. The floor-free ladder is:

| `Das` [steps of TE travel] | 0.04 (2 rows) | 0.08 | 0.40 | 1.60 |
|---|---:|---:|---:|---:|
| CT | 0.03362 | 0.05215 | 0.07133 | 0.07190 |

CT climbs steeply from small `Das` and saturates by ≈0.4 steps.

**2. The `nwakerows` negative is confirmed at settled state.** η=0.1 with 2 rows
has the **larger** total handoff distance (≈1.04 steps vs 0.40 for η=1.0,
nrows=1) and **less than half** the CT (0.03362 vs 0.07133). Free convecting rows
do not substitute for a longer rigid `Das`. Comparing at equal row count is even
starker: `Das` 0.08 → 0.04 steps drops CT 0.05215 → 0.03362 (−36%) *despite* the
extra row. **`Das` length, not wake extent, is the controlling variable.**

**3. Ripple note — do not over-read.** The low-CT cases show lower *relative*
ripple (1.35–1.49% vs 2.2–2.7% mid-ladder), but absolute scatter largely tracks
CT magnitude (std 0.00046 → 0.00159 as CT rises), so much of that is
proportional. The genuine outlier remains **η=4.0: high CT (0.0719) with low
relative ripple (1.40%)**. No configuration yet meets 006's 0.5% / 2% gates.

**Consequence for this item's acceptance.** A ±37% dependence on a parameter with
no principled value is a very large model-form uncertainty — larger than the
entire gap this project set out to explain (panel 0.0506 vs experiment 0.072 is
+42%). Resolving 014 is therefore not a refinement of 006's result; it
substantially determines whether that result means anything.

## Un-freezing `Das`: first attempt INVALID (2026-07-29) — the question is still open

The "un-freeze `Das`" candidate direction above was attempted
(`p2e_das1p0_refresh`, job 12927924) and **cancelled as invalid**; no conclusion
about it can be drawn.

**Why it failed.** `Backslash` caches an LU factorization built once at
construction, and `G` includes the attached wake panel, so **`G` depends on
`Das`** — stated in `src/FLOWPanel_formulation.jl:357`. The driver skips
`initialize_Das!` when `DAS_REFRESH=true` (to avoid `+=` double-accumulation),
but the solver is constructed *after* that skip, so **`G` was factorized with
`Das = 0`**. Bound circulation collapsed **21.7×** (Γ 0.0105 vs 0.2279) and CT
fell to 0.0076 vs 0.0713.

**Constraint this places on the whole item.** Any experiment that varies `Das`
*during* a run is incompatible with a cached-operator solver unless `G` is
rebuilt. This is a real obstacle to the un-freeze test specifically, and it is
worth noting that it does **not** affect the frozen-`Das` ladder (η=0.2/0.5/1.0/
4.0), where `Das` is fixed before solver construction and only rotated
thereafter — those results stand.

**Re-attempt options:** retain the pre-init before solver construction (accepting
a slightly stale `G`), rebuild `G` when `Das` changes, or use a non-caching
solver. Note that during spin-up `Das` varies with ω by construction, so option 1
is only defensible if judged on the post-spin-up window.

**This also strengthens the case for the steady prescribed-wake sweep** (section
above): in a steady solve `Das` is fixed before the operator is built, so the
caching problem does not arise at all.

## HPC JOB STATUS + NEXT PLAN (2026-07-29) — supersedes earlier "in flight" statements

**All earlier "in flight" claims in this file are STALE.** Authoritative state:

| Job | Case | State | CT (10-rev, revs 10–19) |
|---|---|---|---|
| 12925945 | `p2e_das0p2_nofloor` | COMPLETED | 0.05215 ± 0.00078 |
| 12925947 | `p2e_nrows2_dassmall` (η=0.1, rows=2) | COMPLETED | 0.03362 ± 0.00046 |
| 12927923 | `p2e_nrows4_das1p0` (η=1.0, rows=4) | COMPLETED | **0.07431 ± 0.00096** |
| 12921071 | `p2e_das1p0_arc` | COMPLETED | 0.07108 (arc immaterial, <0.5%) |
| 12927924 | `p2e_das1p0_refresh` | **CANCELLED — INVALID** | see the refresh section above |
| **12921559** | `p2e_nt72_das1p0` | **RUNNING** (21 h, step ~1174/1438) | pending — F5's pre-registered prediction |

### `nrows4` result — free rows are NOT inert (corrects a pre-registered prediction)

At fixed `Das` (η=1.0), going `nwakerows` 1 → 4 gives **0.07133 → 0.07431
(+4.2%)**, with the lowest ripple of any high-CT case (1.30% / 4.52%). The
pre-registered "CT unchanged ⇒ free rows inert" branch is **wrong**.

This is a wake-**extent** effect (a longer panel sheet), not a subdivision
effect, and should saturate as rows → many. Scale comparison: free rows +4.2%
vs `Das` **+37%**. So `Das` dominates by ~an order of magnitude, but the
`nwakerows` axis is live and Proposal 1 is narrowed, not closed.

Note both `nrows4` (0.07431) and the η=1.0 + 0.5R filter combination (0.07625)
sit **above** 006's 0.068–0.072 band.

### Checking 12921559 after a context clear

Harvest command and the two gotchas (steps/rev = **72**, and
`scripts/p2e_status.sh`'s `reason=nan_inf` is a **false positive** on completed
runs) are written out in the "HPC JOB STATUS" section of
`BRAINSTORM/006_rotor_hover_converge_stable_near_reference.md`. Always confirm
the terminal state with `sacct -j 12921559` first — a running job's prefix is not
evidence.

### The next experiment is planned and ready to execute

**→ `plans/20260729_wake_representation_consistency.md`** — self-contained, no
exploration needed.

It leads with a **representation-consistency check** rather than another `Das`
sweep, because a steady uniform-Γ `Das` sweep is a *theorem*: constant-strength
doublet panels have exactly cancelling interior edges, so subdivision of a
uniform sheet is analytically irrelevant and the sweep would measure only
round-off. The open question it exposes instead is whether the code's **two wake
representations agree at all** — which the entire η ladder presumes.

Plan in brief: add a `PanelWake` mode that convects **all** rows with the
freestream only (today `shed_with_induced_velocity=false` does this for row 1
only), so the panel wake's geometry is *exactly* the semi-infinite wake's; then
check convergence to the semi-infinite answer as wake length grows, reusing
`examples/suddenly_started_wing.jl` which already has the impulsive start, the
`PanelWake`, an `SSW_ETA` `Das` knob and a unit-`Das` semi-infinite reference
(`_ssw_steady_cl`). Non-convergence would be a **direct candidate for the +37%**.

**Results of that plan are to be recorded HERE (014) as the headline verdict**,
with the run ledger in `logs/dji_convergence_20260722/phase_02e_unsteady_ct.md`
and the outcome cell in `BRAINSTORM/INDEX.md`. 006 only needs updating if the
outcome changes the standing of its 0.0713 magnitude claim.

## VERDICT 2026-07-29: the two wake representations ARE consistent — this candidate for the +37% is ELIMINATED

Executed `plans/20260729_wake_representation_consistency.md`. **Result: negative
for the hypothesis, clean and unambiguous.**

`PanelWake` gained a `freestream_convection` mode (`src/FLOWPanel_wake.jl`) in
which *every* row convects with the freestream only, so the sheet stays straight
along `U∞` — geometrically *identical* to the semi-infinite wake, truncated at
length `L`. Driver: `SSW_MODE=wake_consistency` in
`examples/suddenly_started_wing.jl`; data in
`data/suddenly_started_wing/wake_consistency.csv` (+ `.png`). AR=6, n_span=12,
n_airfoil=21, dt*=1/8, `DirectBackend` (mandatory — FMM inflates the panel
radius by `|Das|`, see below), `include_final_filament=false`, η=0.25.
`cl_semiinf = 0.42122` from `_ssw_steady_cl` (unit-`Das` semi-infinite solve).

| `L/c` | rows | `cl` straight sheet | rel err | `cl` rolled-up (control) | rel err |
|---:|---:|---:|---:|---:|---:|
| 1 | 8 | 0.33912 | −19.49% | 0.33459 | −20.57% |
| 2 | 16 | 0.37425 | −11.15% | 0.36993 | −12.18% |
| 4 | 32 | 0.40212 | −4.534% | 0.39801 | −5.509% |
| 8 | 64 | 0.41604 | −1.231% | 0.41207 | −2.173% |
| 16 | 128 | 0.42011 | −0.263% | 0.41616 | −1.202% |
| 32 | 256 | 0.42098 | −0.0577% | 0.41701 | −0.999% |
| 64 | 512 | 0.42116 | **−0.0136%** | 0.41719 | −0.956% |

**W1 — the straight-sheet panel wake converges to the semi-infinite answer.**
Monotone, and to **0.014%** at `L/c=64`. §5's criterion (<1% by `L/c≈32`) is met
by `L/c≈8`. The successive-error ratio over the last three doublings is
0.213 / 0.220 / 0.236, i.e. an **asymptotic order ≈ 2.2 in wake length**.
(The `_ssw_convergence_order` log-log fit reports +1.81, but that fit is over the
whole range and is dragged down by the shallow small-`L` region; the asymptotic
ratio is the meaningful number.) Second order is the *expected* physics: the
truncated sheet's missing far wake is closed by the starting-vortex loop at
distance `L`, whose induced velocity at the wing falls off as `1/L²`.

**W2 — the rolled-up control plateaus at −0.96%, as it should.** It does *not*
converge to the semi-infinite value, and must not: a rolled-up wake is a
genuinely different (drooped) geometry. That residual is a measurement of the
physical rollup effect on CL at this configuration, and it also confirms the
straight-sheet result is not an artifact of the comparison machinery — the same
machinery gives a non-zero answer where a non-zero answer is correct.

**W3 — not AR-specific.** Repeated at AR=3 (`data/suddenly_started_wing_ar3/`):
−0.376% / −0.084% / −0.020% at `L/c` = 8 / 16 / 32, ratios 0.22 / 0.24, order
≈2.1. Same convergence, same rate.

**Consequences.**
1. The plan's leading hypothesis is **dead**. The code's two wake
   representations agree to 0.014% when given identical geometry, so a
   representation inconsistency is **not** a candidate for the +37%. The §6
   diagnostic ladder (streamwise Γ uniformity, attached-strip-vs-row-1 strength,
   final-filament treatment) was gated on non-convergence and is **not
   triggered** — do not spend on it.
2. The whole η ladder's presumption that the representations agree is
   **verified**, not merely assumed.
3. The live follow-on is §7's "if it converges" branch — and specifically its
   *second* bullet, which is **not** gated on the outcome: sweep `Das` at fixed
   *large* `nwakerows`, so particles sit far downstream and the near wake is all
   panels at near-uniform Γ. Collapse of the sensitivity there ⇒ the +37% is the
   **panel→particle representation change**, routing to `BRAINSTORM/012`.
   `p2e_nrows4_das1p0` (η=1.0, nrows=4, CT 0.07431) is one point; it needs
   **η=0.2 at nrows=4** to form the pair.

### Two mesh defects found and fixed en route — they affect earlier oracle numbers

Both are in `suddenly_started_wing_mesh`, and both are *topological*: the node
cloud was symmetric, the triangulation was not.

- **Chordwise / trailing edge.** Every quad was split along the same diagonal,
  so the upper- and lower-surface TE panels were not reflections. At a given
  spanwise station the lower TE control point sat at x=0.99184, z=−0.00117 while
  the upper sat at x=0.98369, z=+0.00233 — **0.8% chord apart with double the
  |z|**. The Kutta condition was being enforced on an asymmetric pair.
- **Spanwise (found by Ryan's question "is the wing geometry mirrored about the
  x-z plane?" — it was not).** Under y→−y that same diagonal maps onto the
  *other* diagonal of the image quad, so **all 480 cells violated spanwise
  symmetry**. Measured consequence: **~10% asymmetry in shed circulation on a
  geometrically symmetric wing, which did NOT refine away** — 1.01e-1 / 1.13e-1 /
  1.10e-1 at n_span = 12 / 24 / 48. A fixed topological bias, not truncation
  error.

Fixed by keying the splitting diagonal on the **XOR** of the chordwise and
spanwise half-tests, so the split flips under each reflection. Spanwise Γ
asymmetry is now **1.3e-14** at all three resolutions and the TE pairs are
exactly mirrored. Requires even `n_span` (now enforced; an odd centre column
would have to be its own mirror image). Guard: `assert_ssw_mesh_symmetry`,
called at the top of `wake_consistency` mode — it passes on the corrected mesh
and flags all 480 cells of the old one.

**Effect on earlier results in THIS file — discharged 2026-07-29.**
`examples/ssw_eta_convergence.jl` includes `suddenly_started_wing.jl` and calls
`run_suddenly_started_wing`, so the η×dt matrix (F1–F4), the `Das` log law and
the long-window addendum had all been computed on the spanwise-asymmetric mesh.
**They have now been fully re-run on the corrected mesh and the results section
above replaced.** All four findings stand; the numbers moved by <1% and one loose
statement of F2 was corrected. The case tag now carries `SSW_MESH_REVISION`, so a
future mesh change cannot silently resume histories from the previous geometry —
which is precisely how this defect propagated unnoticed.

Also note steady CL at AR=6 moved with the mesh fix: n_span 12 / 24 / 48 now
give 0.42122 / 0.38005 / 0.36491 (was 0.39280 / 0.38019 / — ). The n_span=12
value in particular is coarse; the sequence is still converging, so absolute
levels from n_span=12 should be treated as indicative.

### Fourth candidate mechanism, untested — record before it is forgotten

`src/FLOWPanel_liftingbody.jl:996` inflates the FMM panel radius by `|Das|`
(`buffer[4,i] += max(|Da|,|Db|)`), so the multipole acceptance criterion is an
**explicit function of `Das`**. **Every unsteady rotor η run used
`FastMultipoleBackend`**, so FMM accuracy varied along the ladder. This is
independent of the three refuted mechanisms and untested. Cheapest check: rerun
two η ladder points with `DirectBackend`.

*Status 2026-07-29:* still untested on the rotor. The wake-consistency study
above was run entirely on `DirectBackend` precisely to keep this out of the
measurement, so it neither confirms nor refutes this mechanism. It remains the
cheapest open discriminator.

## RESULTS 2026-07-29: rotor discriminators overturn two hypotheses — the sheet/particle split governs CT

Full numbers and table in `logs/dji_convergence_20260722/phase_02e_unsteady_ct.md`
(2026-07-29 entry). Headlines: **floor-off η=0.2 gives CT = 0.0522** (the 0.01R
floor was *raising* `Das` inboard, so the ladder's low end was floor-propped and
the true η 0.2→1.0 sensitivity is +37%); **tiny-`Das` + free row gives 0.0336**
(the "only total handoff distance matters" hypothesis is dead — Proposal 1's
optimistic branch is closed); **η=1.0 + 4 rows gives 0.0743** (+4.2%, panel-
sheet extent is its own lever, above the band). Combined with 001's semi-
infinite-sheet 0.110: **CT rises monotonically with the fraction of near wake
represented as panel sheet vs particles, spanning 0.034–0.110 across plausible
discretizations, with experiment (0.072) inside the range.** The η=1, 1-row
0.0713 therefore cannot yet be called a prediction: the required demonstration
is a joint plateau in (`Das`, `nwakerows`) — e.g. does the nrows ladder
asymptote, and to what — or a principled near-wake treatment (Proposal 2 /
possibly a panel-sheet buffer long enough that particles only carry the far
wake, cf. truncation-depth null at 4R). NT=72 (12921559) pending as the
dt-refinement check with F5's 0.0708–0.0712 pre-registration. Oracle matrix is
being re-run on the corrected rev-2 mesh (rev-1 triangulation carried ~10% shed
asymmetry) to confirm the log-law findings.

## RESULTS 2026-07-29 (wing sheet/particle ladder): particles ≈ sheets at wing σ; σ/chord is the bridge variable to the rotor

Ryan's requested experiment — does a simple rectangular wing's circulation
depend on how much of the wake is panel sheet vs particles? Machinery:
`wake_model=:particle` mode in `examples/suddenly_started_wing.jl` (final-edge
filament ACTIVE, OverlapPPS unsteady shedding so the starting vortex survives
conversion) + sweep `examples/ssw_sheet_particle_split.jl`. AR=6, η=1,
dt*=1/8, t*=20, rev-2 mesh; data `data/ssw_sheet_particle_split/`
(`split_summary.csv`, `sheet_particle_split.png`). Two findings:

1. **At wing-scale σ (σ/c=0.08) the split does NOT matter**: buffer 0.25c→8c
   all within +0.06%…+0.01% of the pure-panel control. The rotor's sensitivity
   is NOT a generic particle deficit — pre-registered "flat" branch fires.
2. **The split matters exactly when σ ~ chord** (tail ΔCL vs pure panel;
   σ/c = 0.0625·overlap):

   | σ/c | buffer 0.25c | 1c | 4c |
   |---:|---:|---:|---:|
   | 0.08 | +0.06% | +0.04% | +0.03% |
   | 0.19 | −0.20% | −0.18% | −0.08% |
   | 0.38 | −0.20% | −0.32% | −0.12% |
   | 0.75 | **+0.77%** | −0.36% | −0.17% |
   | 1.5 (rotor-like) | **+3.32%** | +0.42% | −0.44% |

   Non-monotone with a sign flip: chord-scale-σ particles NEAR the body
   over-induce (smoothed core overlaps body/near-wake); farther away they
   under-induce. Buffer-length spread grows <0.1% → 3.8% as σ/c goes
   0.08 → 1.5. On the rotor (σ ≈ 1.5 local chords AND inflow→loading feedback
   that a fixed-α wing lacks) this is amplified into the observed
   0.034–0.110 CT range.

**Verdict for the convergent methodology:** the sheet/particle split is a
σ/chord problem, not a `Das` problem. `Das`/η is settled (physical length on
the log plateau, ±≲1%); the remaining requirement is **σ/chord ≪ 1 near the
body, or sheet representation until the wake is far enough that σ-smoothing is
harmless** — item 012's territory, now with a quantitative wing oracle for it.
Rotor-side: σ-refined unfiltered case `p2e_sigF_nofilt` submitted (12943696,
σ_med≈0.048R ≈ 0.37 local chords, 3× better than stock 1.13); sheet-extent
asymptote ladder `p2e_nrows{8,16,36}_das1p0` in queue (12943025–27).

**Bug fixed en route:** `_shed_particles!` admitted zero-strength particles;
|Γ|=0 NaNs corrected-Pedrizzetti relaxation and poisons the field (exact zeros
arise on symmetric wings and impulsive starts; the rotor never hits exact
zeros). Guard in `src/FLOWPanel_wake.jl`; wake unit tests pass.

## 2026-07-29 — Legacy-default riders flagged by the item-015 review

A clear-context integration review (run for item 015's Phase 3 approval)
confirmed two uncommitted item-014 changes in this worktree change
default-kwargs legacy trajectories and need their own attribution/regression
evidence before any commit that claims "legacy defaults unchanged":

1. **`set_Das_kinematic_arc::Bool=true` is a behavior-changing new default**
   (`src/FLOWPanel_simulate.jl` and `simulate_warmstart!`): kinematic `Das`
   accumulation now routes through `accumulate_Das_arc!` instead of
   `_accumulate_Das!`. The two agree only to first order in |ω|τ and are
   identical only for purely translating bodies, so every default-kwargs
   rotating-body run (i.e. every rotor case) gets a different first-wake-row
   offset and trajectory. Either default it to `false` or document it as an
   intentional 014 default change with regression evidence.
2. **Zero-Γ `_shed_particles!` early return** (`src/FLOWPanel_wake.jl`,
   the NaN fix above): legitimate, but it changes `PanelParticleWake` legacy
   particle counts/indices/VTK output for runs that previously created
   zero-strength particles (symmetric wings, impulsive starts). Should be
   attributed here, not silently bundled with unrelated work.

## Ruling 2026-07-30 (Ryan): N=36 rows is NOT adopted; isolate the sheet/particle variability first

The `nwakerows` ladder (1 → 0.07133, 4 → 0.07431, 8 → 0.07506, 16 → 0.07304,
36 → 0.07049) is **non-monotone — not a convergent metric** — so the N=36
configuration's in-band, low-ripple result is *not* adopted as the answer.
Continue investigating before adopting any sheet-buffer prescription.

**Directed next study: isolate the CAUSE of the sheet/particle variability with
a targeted experiment — a straight wing in steady conditions.** The natural
platform is the existing `examples/ssw_sheet_particle_split.jl` harness (or
`code_audit/scripts/task1_steady_wake_consistency.jl`): run to a settled steady
state at rotor-like σ/chord and vary the conversion point, σ, and core
treatments *while measuring more than CL* — e.g. compare the particle region's
induced velocity at the body/near-wake against the pure-sheet control field
directly, so the over-induction (near) vs under-induction (far) mechanism seen
in the 2026-07-29 wing ladder is localized to a term (kernel smoothing of a
sheet's self-influence, missing sheet coherence, shed-strength lumping, core
overlap with the body) rather than inferred from integrated lift. The study is
convergent only when a scheme makes the wing's split-sensitivity vanish at
rotor σ/chord; that scheme then transfers to the rotor.

### In-flight jobs for the next agent (harvest: 10-rev cycle-means, revs 10–19, from `data/<case>/<case>_CT_per_rev.csv` on the cluster)

| job | case | question |
|---|---|---|
| **12955430** | `p2e_nrows72_das1p0` | two revolutions of sheet — does CT continue below N=36's 0.0705 (no asymptote) or agree with it? |
| **12943696** | `p2e_sigF_nofilt` | σ-refined (OVERLAP=4, PPS=8, σ≈0.37 local chords) unfiltered — does σ refinement move CT at N=1, and is it stable? (006 study C) |
| **12950996** | `p2e_nt72_das2p0_ov6` | NT=72 with `Das` AND σ pinned at NT=36 physical values — isolates dt; ~0.0713 ⇒ dt converged, ~0.0734 ⇒ NT=36 under-resolved in time |

## 2026-07-30 — Sheet/particle question spun off to item 017

The directed study from the 2026-07-30 ruling above is now staged as its own item:
**`BRAINSTORM/017_sheet_particle_representation_equivalence.md`**. It owns the
sheet/particle representation question from here on — the quantitative admissibility map
(σ/c, conversion distance, PPS to tolerance |ΔCL| ≤ 0.5% and per-station |ΔΓ(y)| ≤ 1% on a
steady rectangular wing), mechanism attribution (kernel smoothing / quadrature lumping /
filament-alignment×relaxation / core–body overlap), and the eventual rotor transfer
(including the σ-refined `nwakerows` ladder, the still-open FMM |Das|-radius
`DirectBackend` discriminator, and harvest of the in-flight jobs tabled above). 014 stays
open only for its own residue: the `Das`/η rule (answered in draft; see "Revised answer")
and any remaining rotor confirmation. New sheet/particle results are recorded in 017, not
here.

## 2026-07-31 — In-flight jobs harvested (by BRAINSTORM/018 Phase 0)

- 12955430 `p2e_nrows72_das1p0`: CT̄(10–19) = 0.06931 ± 0.00154 — ladder continues falling past N=36; non-convergence of the legacy-σ nwakerows axis confirmed.
- 12943696 `p2e_sigF_nofilt`: STABLE until a genuine 32G memory ceiling (MaxRSS 33.4 GB) at step 403/719; CF ≈ 0.066 smooth — refined σ does not destabilize the unfiltered wake.
- 12950996 `p2e_nt72_das2p0_ov6` (COMPLETED, harvested 2026-07-31): CT̄(10–19) = 0.06852 ± 0.00157. The matched-physical-units dt isolation lands **−3.9% BELOW** the NT=36 partner (0.07133) — so the naive NT=72's 0.07337 (+2.9%, outside the log-law pre-registration) was indeed confound-dominated: halving Das pushed CT up while true dt refinement pulls it down. F5's pre-registered suspicion vindicated; dt is a live convergence axis (018 Phase 3, full ladder incl. NT=144).

Successor campaign: BRAINSTORM/018 (physical-Das ladder at sigma_overlap shedding, d≥4σ nwakerows criterion instead of a ladder).
