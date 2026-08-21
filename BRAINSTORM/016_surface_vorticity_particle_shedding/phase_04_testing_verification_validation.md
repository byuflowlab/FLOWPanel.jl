# Phase 4 — Testing, Verification, and Validation

**Status:** All planned measurements COMPLETE 2026-08-04 (Phase 3 approved
2026-08-03). §1 and §2 closed; §3 case 1 (wing) PASSED; §3 case 2 (rotor):
T4 attribution null PASSED, smooth-vs-legacy equivalence FAILED on the
unconfounded pair (see §2026-08-04 (c)). Two Ryan-approved attribution
discriminators in flight since 2026-08-05 (larger-d/σ 13051758, equal-h
13051763 — see §2026-08-05). Awaiting clear-context review and user approval.

**Prerequisite:** approved Phase 3 implementation — satisfied 2026-08-03

**Phase approval:** [ ]

## Prior-phase handoff

Phase 1 defines the mathematical identities and acceptance scenarios. Phase 2
must define the approved API, diagnostics, and state contracts. Phase 3 must
provide the implementation and focused mechanical-test record. Read all
approved prior phases before beginning.

## 1. Mathematical verification

Use deterministic synthetic sheets independent of aerodynamic validation:

- Constant \(\mu\): zero distributed particles, with circulation represented
  only by the retained true perimeter closure.
- Affine \(\mu\) on planar and rigidly rotated sheets: exact reconstructed
  gradient, correct package sign, and exact quadrature-integrated particle
  strength.
- Nonuniform planar and warped quad grids: measured convergence under panel
  refinement and particle-spacing refinement.
- Telescoping/discrete-Stokes ledger: internal panel boundaries cancel;
  physical perimeter and panel/particle handoff circulation are each counted
  exactly once.
- Resolution floor: for a nonzero full-rank field, one particle lies at the
  bilinear panel centroid when both subdivision counts are one.
- `nwakerows == 1`: exact one-sided streamwise reconstruction from the newly
  shed active row and outgoing non-source ghost, with no artificial
  outside-zero or `Das` sample.
- Cancelled-edge and retained-filament handoff alternatives: signed
  circulation/flux balance for one row, multiple spanwise panels, and
  repeated conversions.

Record tolerances, observed orders, signed vector balances, geometry scales,
and rank thresholds.

## 2. Mechanical regression

Prove that the new strategy axis is opt-in:

- compare omitted strategy, explicit legacy strategy, and the pre-change
  reference;
- require unchanged particle locations, strengths, counts,
  `method_trailing`/`method_unsteady` behavior, metadata, replay, and warm
  starts;
- round-trip smooth-strategy metadata and continuation state;
- verify active source counts remain \(N\) despite \(N+1\) stored panel
  strengths and \(N+2\) node rows;
- verify the ghost never contributes panel influence;
- exercise validation and failure paths; and
- force insufficient capacity and verify that neither particles nor panel
  rows are partially committed.

Run the narrow wake, simulate, replay, and warm-start suites selected by the
approved Phase 3 test map before broader regression.

## 3. Aerodynamic validation

Compare legacy edge-jump and smooth surface-vorticity strategies with
identical geometry, wake, particle, solver, timestep, maintenance, viscous,
SFS, and monitor settings. Only the conversion strategy may differ.

### Case order

1. **Suddenly started wing.** Establish signs, startup circulation,
   lift history, induced velocity, and convergence in a compact transient.
2. **Rotor hover.** Evaluate the production-relevant helical wake only after
   the wing case and all mathematical/mechanical gates pass.

### Quantities and refinement axes

For both cases compare:

- induced velocity at fixed probes and relevant surfaces;
- integrated circulation and impulse;
- particle count and spatial distribution;
- numerical stability and conservation drift; and
- \(C_L\) for the wing or \(C_T\) for the rotor.

Refine panel-row count, spanwise panels, timestep, \(\sigma\), overlap, and
root/tip line-particle spacing. Report cost and particle count alongside
accuracy so apparent smoothing caused only by coarser effective resolution
is visible. Separately report panel-induced velocity-gradient norms at new
root/tip and interior particle locations. Non-finite values, instability, or
nonconvergent growth trigger reconsideration of a conservative diffuse
root/tip closure. Negative results, lack of convergence, or sensitivity
merely transferred to another parameter are valid conclusions.

## Acceptance gate and progress log

Phase 4 and Item 016 are complete only after mathematical verification,
mechanical regression, and both ordered aerodynamic validations satisfy the
approved criteria and the user explicitly approves the conclusions. Only then
may the item-level columns in `BRAINSTORM/INDEX.md` be checked.

- **2026-07-29 — Phase seeded.** Blocked pending explicit Phase 3 approval;
  no tests or simulations have started.

## Carried in from Phase 3 (2026-08-01)

- **Streamwise attribution A/B.** `SurfaceVorticityConversion` accepts
  `attribution = :upstream` (default) | `:downstream` | `:split`. All three
  conserve circulation exactly; `:split` is second-order accurate where the
  others are first, but retains a half-jump filament at the panel/particle
  handoff. See Phase 1 §9b.4.

  **Corrected 2026-08-03 (review finding R3).** Stage 4b did *not* settle the
  default. It failed to displace it: no Stage 4b measurement favours
  `:upstream`; the only axis with a resolved sign (refinement order) favours
  `:split`; and the verdict came from a conjunctive rule written after the data
  was in, which can only ever return the incumbent. The near-field axis cannot
  resolve a sign at all — it scores fidelity to the edge-jump sheet being
  replaced, over M = 4 points on one grid, with normalization that drives the
  curve shape. `:upstream` is **retained pending a discriminating measurement**.
  Stage 4c (pre-registered in the Phase 3 document) is that measurement.

- **T4 — Phase 4 attribution null, PRE-REGISTERED 2026-08-03 before any Phase 4
  run.** When the ordered aerodynamic validation reaches the rotor:

  > If `|Delta CT|` between `alpha = 1` (`:upstream`) and `alpha = 1/2`
  > (`:split`) is below the item-018 campaign scatter — baseline B0 = 0.07170
  > with CI [0.07164, 0.07173], i.e. a half-width of about 5e-5 in CT (~0.07%) —
  > then attribution is documented as a **free parameter of measured size** and
  > carried at that magnitude in the CT error budget.

  The null is framed as a positive finding with a number attached, in the form
  item 014 produced for `Das`. Both arms must be run at identical geometry,
  wake, particle, solver, timestep, maintenance, viscous, SFS, and monitor
  settings; only `attribution` may differ. If Stage 4c's T0 finds no `h_row`
  collapse, this comparison is mandatory rather than confirmatory.
- **The startup edge is physical.** Before the first handoff the outgoing row's
  aft face is the sheet's trailing boundary (the starting vortex) and is
  deposited whole. Any validation that compares a run's total shed circulation
  against theory must account for it being present, not suppressed. A deliberate
  suppression switch, if wanted for rotor hover (cf. item 005), is a separate
  option and does not exist yet.
- **Conservation must be checked externally.** A residual measured against the
  transaction's own `expected_total` is self-referential and cannot detect an
  omitted face; that is how the startup leak survived the Stage 3 gate. Use the
  filament ledger: the wake's field-relevant content reduces to the retained
  final filament, and the particles' gain must equal its loss.

## 2026-08-03 — driver wiring + readiness bundle (prepared for Ryan's Phase 3 approval)

Requested by Ryan in the BRAINSTORM/018 session: "let's try the shedding
technique developed in 016 … let me know so I can verify that it's ready to
use." **Nothing here starts Phase 4 validation** — no aerodynamic case has been
run, and this section does not claim any physics result.

### What was wired

`examples/rotor_hover_pressure_comparison.jl` gained an opt-in switch:

| ENV | default | meaning |
| --- | --- | --- |
| `CONVERSION` | `legacy` | `legacy` = historical edge-jump (unchanged); `smooth` = `SurfaceVorticityConversion` |
| `CONVERSION_SIGMA` | tip σ = `2πR/NT · OVERLAP/PPS` | smoothing width; defaults to the σ the `sigma_overlap` policy already uses, so `smooth` sheds at the ladder's σ |
| `CONVERSION_OVERLAP` | `1.3` | sampling overlap (h = σ/overlap) |
| `ATTRIBUTION` | `upstream` | `:upstream` / `:downstream` / `:split` |
| `CONVERSION_DIAGNOSE` | `false` | `diagnose_nearfield` |

The two strategies need mutually exclusive wake options, so the differing
kwargs are built once and splatted: `legacy` passes
`method_trailing`/`method_unsteady=NoShed()`/`unsteady_filament=false`;
`smooth` passes `conversion=…`/`unsteady_filament=true` and passes **neither**
method (supplying either alongside the smooth strategy is a configuration error
the constructor rejects, and `unsteady_filament=false` is rejected by the R1
guard — the driver was previously configured in *both* of those rejected ways,
so `conversion=` could not simply be added).

`conversion` (and σ/overlap/attribution when smooth) is now written to
`_case_metadata.toml` and the diagnostics banner, so a run's conversion is
auditable after the fact — the same lesson as the BRAINSTORM/018 η leak.

### Evidence collected

- **Legacy path unchanged.** `runtests_unit_wake.jl` **583/583**, including the
  Stage-0 golden reference (33 exact-equality assertions against
  `test/data/legacy_wake_conversion_reference.jl`). `runtests_unit_replay.jl`
  125/125, `runtests_unit_warmstart.jl` and `runtests_unit_simulate.jl` clean.
  Zero Fail/Error anywhere (judged by log content — `julia` exits 0 on a load
  error).
- **The smooth path time-marches.** Coarse smoke (`RHPC_MESH=40_40`, NT=6,
  ~2.3 revs, `NWAKEROWS=1`, σ_conv = 0.1246) completed for both modes with
  `all finite = true`.
- **It is not inert.** Same smoke, matched step 15: legacy `CT_bernoulli`
  0.081756 vs smooth 0.083010 (**+1.5%**), diverging smoothly from step ~10.
  This is a *does-it-run* check on an unsettled coarse case and is **not**
  evidence about CT.

### What is still NOT true, and must be said plainly

1. **No aerodynamic validation exists.** All of this phase's ladder —
   mathematical verification on synthetic sheets, the suddenly-started wing,
   then rotor hover — is unstarted. A rotor CT from this path today would not
   be trustworthy.
2. **No end-to-end time-marched *test*.** `runtests_unit_simulate.jl` contains
   no reference to the conversion; the smoke above is a manual run, not a
   regression test.
3. **R4 is still open** (evidence integrity, not a code defect): Stage 4b's
   `probe_errors.csv` was computed with the retained filament emitting nothing,
   because `get_n_bodies(::FilamentWrapper{<:PanelWake})` returns zero unless
   `overflowed[]` is set, so that table cannot have measured the
   `attribution`-dependent effect it was read as measuring. `:upstream` is
   retained on Stage 4c's evidence, not on that table.
4. **`include_final_filament` is not serialized** in the `PanelParticleWake`
   manifest branch (`FLOWPanel_replay.jl:341-360`), so it silently round-trips
   to `true`. Benign for the smooth strategy (which requires `true`), but a
   latent replay hazard.
5. **Phase 3 approval is still unchecked**, and Phase 4 is formally blocked on
   it. Using this on the rotor IS Phase 4 work.

### Suggested first use, if approved

016's own ladder says wing before rotor. The cheapest honest first step is the
suddenly-started wing (`examples/ssw_*`), where BRAINSTORM/014/017 already
provide a converged sheet-vs-particle reference, **not** a rotor CT — a rotor
number would arrive with no way to tell a conversion error from the campaign's
open σ and dt errors.

## 2026-08-03 — §1/§2 evidence record: coverage map + gap-closure tests

Phase 3 was approved and Phase 4 opened this session. §1 and §2 were closed by
(a) mapping every spec bullet to the testset that proves it and (b) writing new
tests only for the bullets no existing testset covered. All testsets below live
in `test/runtests_unit_wake.jl` unless noted; suites after the additions:
**wake 605/605** (was 583), **simulate 88/88** (was 81), replay 125/125,
warm start 26/26, full `test/runtests.jl` exit 0 with zero `Fail`/`Error`.

### §1 mathematical verification — coverage map

| Spec bullet | Evidence (testset) | Status |
| --- | --- | --- |
| Constant μ → zero distributed particles, perimeter-only circulation | "constant strength deposits nothing" (zero field: np==0, elided; nonzero constant: zero area strength, root/tip ±0.7, residual <1e-14) | pre-existing |
| Affine μ: exact gradient, package sign, rotation covariance | "constant and affine planar fields", "package sign convention", "rigid-rotation covariance" (unit level) | pre-existing |
| Affine μ: exact quadrature-integrated particle strength END TO END | **NEW** "affine field deposits the exact integrated strength end to end" — per-panel particle sums from the particle field equal `deposited_strength` and the interior panel equals κ_exact·A to 1e-12 | new |
| Nonuniform planar + warped grids: measured convergence under panel refinement | **NEW** "kappa converges on nonuniform and warped sheets" — smooth non-affine μ (sin), graded rows with ratio→1, vs analytic pointwise κ (warped case uses \|dμ/ds\| on z=f(x)). Measured orders: nonuniform planar **1.62, 1.39**; warped **1.59, 1.45** (band asserted 0.6–2.2, errors monotone; coarsest warped rel-err 0.47 declining). Complements the pre-existing kappa-rule gap order (:485, first order) and split-attribution order (0.7–1.4 vs 1.6–2.5) | new + pre-existing |
| Particle-spacing refinement | "subdivision follows sigma / overlap" — deposited totals are *invariant* (rtol 1e-7 area / exact face content), which is the specified behavior: spacing refines sampling, not content | pre-existing |
| Telescoping/discrete-Stokes ledger | "circulation is conserved exactly on arbitrary geometry" (vs independent `expected_transfer`, <1e-12), "transfer matches the wake's own filament ledger" (external, `get_sources`-derived), "repeated conversions keep the ledger", "only true root and tip lines survive" (wrap → n_perimeter==0) | pre-existing |
| Resolution floor: ONE particle AT the bilinear centroid | **NEW** "resolution floor: one particle at the bilinear centroid" — (1,1) counts AND position == mean of the four corners to 1e-13, one per panel | new (position clause; counts pre-existing) |
| `nwakerows == 1` one-sided reconstruction, no Das/outside-zero sample | "no upstream geometry is read at nwakerows == 1" — body translation + 3×Das bit-identical; upstream strength matters; absent body throws | pre-existing |
| Handoff alternatives balance (one row, multi-span, repeated) | "the startup edge is deposited, not deleted", "the retained filament matches the attribution" (`:downstream` bit-identical legacy), "legacy and smooth differ only by the handoff share", ledger sweep attribution × nwakerows × wrap | pre-existing |

### §2 mechanical regression — coverage map

| Spec bullet | Evidence | Status |
| --- | --- | --- |
| Omitted vs explicit legacy vs pre-change reference | "Legacy edge-jump conversion golden reference" (33 exact equalities vs pre-016 capture), "Panel-particle conversion strategy axis" | pre-existing |
| Unchanged locations/strengths/counts/method_*/metadata/replay/warm starts | golden reference + strategy axis; legacy metadata round trip (`runtests_unit_replay.jl` "strict conversion metadata and continuation round trip"); legacy warm start (`runtests_unit_warmstart.jl`) | pre-existing |
| Smooth metadata + continuation round trip | replay round-trip testset (all six fields, invalid `unsteady_filament=false` manifest rejected); warm start across a conversion boundary | pre-existing |
| Active source count N despite the extra stored rows | **NEW** "active source count stays N; the stale row has no influence" — `get_n_bodies(pw) == N·ncols` with N+1 stored strength rows, unchanged across a conversion. **Note:** the bullet's "N+2 node rows" is Option A phrasing; under the approved Option B (Phase 2 §13.1) there is no ghost storage — nodes/strengths are N+1 rows and the ghost is the pre-shift final active row, itself a legitimate source | new (Option B mapping recorded) |
| Ghost never contributes panel influence | same NEW testset: what must be invisible under Option B is the stale strength row N+1 — after the smooth handoff (`:upstream`) retargets the filament, perturbing row N+1 leaves probe velocities bit-identical while perturbing an active row moves them (control). Also pins R4 at source level: `get_n_bodies(FilamentWrapper) == 0` until `overflowed[]`, `== ncols` after | new |
| Validation and failure paths | constructor validation, line-policy conflicts, sentinel, R1/R2 guards ("smooth conversion requires its final-filament sources"), attribution validation, geometry rejection | pre-existing |
| Transactional capacity failure | "failures are transactional" (capacity, NaN geometry, premature nwakes — full state snapshots) | pre-existing |
| End-to-end time-marched smooth test (driver-wiring caveat 2 above) | **NEW** `runtests_unit_simulate.jl` "simulate! smooth surface-vorticity conversion smoke" — real `simulate!` loop over 6 steps with buffer overflow: conversions fire, particles finite, external ledger totals close to <1e-10. Uses `freestream_convection=true` (the noop test solver leaves O(1) strengths whose induced velocity folds the crude synthetic wake — a fixture constraint, not a code limitation) | new — closes caveat 2 |

§1/§2 are **complete for this session's scope**. §3 case 1 (SSW A/B) runs via
the new `examples/ssw_conversion_ab.jl`; results recorded separately below.

## 2026-08-03 — §3 case 1 COMPLETE: suddenly-started-wing legacy-vs-smooth A/B

Driver `examples/ssw_conversion_ab.jl`; wiring in
`examples/suddenly_started_wing.jl` (`conversion=:smooth` opt-in; legacy case
tags verified byte-identical to pre-change, smooth appends `_svc_…`).
Setup: AR=6, n_span=24, n_airfoil=21, η=1.0, panel_rows=8, t*∈(0,20],
`wake_model=:particle`, FMM backend (per Ryan; the dt*=0.25 pair was
cross-checked against DirectBackend and reproduces ΔCL to 4 digits). Within
every pair ONLY the conversion strategy differs; the smooth σ default resolves
to the legacy unsteady-row σ = `pps_overlap·dt*/pps_n` so effective resolution
is matched. The pre-registered interpretation block sits in the driver header
and was committed before the first run. Data:
`data/ssw_conversion_ab/conversion_ab_summary.csv` (+ per-case artifacts,
`conversion_diagnostics.csv`, `probe_comparison.csv`, `conversion_ab.png`).

| label | dt* | σ_conv/c | ΔCL % (smooth−legacy) | Γ_TE err % | loading err % | both settled | ledger residual |
| --- | --- | --- | --- | --- | --- | --- | --- |
| primary | 0.25 | 0.1625 | +0.035 | 0.043 | 0.038 | yes | 1.4e−16 |
| primary | 0.125 | 0.08125 | −0.120 | 0.115 | 0.115 | yes | 3.5e−17 |
| primary | 0.0625 | 0.040625 | −0.090 | 0.091 | 0.094 | yes | 2.0e−17 |
| sigma ×0.5 | 0.125 | 0.040625 | −0.044 | 0.043 | 0.041 | yes | 1.3e−17 |
| sigma ×2 | 0.125 | 0.1625 | −0.176 | 0.170 | 0.172 | yes | 2.4e−17 |
| overlap 1.0 | 0.125 | 0.08125 | −0.121 | 0.116 | 0.116 | yes | 1.9e−16 |
| overlap 1.6 | 0.125 | 0.08125 | −0.120 | 0.114 | 0.114 | yes | 1.4e−16 |
| attribution :split | 0.125 | 0.08125 | −0.120 | 0.115 | 0.115 | yes | 2.1e−17 |

**Findings (per the pre-registered framing):**

1. **Stability/settledness: PASS.** Every arm of every pair settles under
   `ssw_settled_stats`; smooth CL ripple ≤ legacy's at dt*=0.25 (0.014% vs
   0.020%); all quantities finite; no growth.
2. **Conservation: exact.** The external filament-ledger residual is ≤2e−16 in
   all 8 smooth runs (72–1104 conversions each) — the §1 round-off-exact
   invariant survives real time marching.
3. **Magnitude: smooth ≡ legacy to ~0.1%.** |ΔCL| ≤ 0.18% everywhere, an order
   below the 0.5% admissibility band item 017 used for shedding variants; Γ_TE
   and loading profile errors track ΔCL 1:1. The dt* trend does not grow toward
   dt→0 (−0.120% → −0.090%).
4. **The residual delta is σ-linked, not scheme-intrinsic:** at fixed dt* it
   scales ~linearly with conversion σ (−0.044/−0.120/−0.176% at σ/c
   0.041/0.081/0.163) and is inert to overlap (sampling density) — consistent
   with a kernel-smoothing-scale effect, the effect §3 said to expose, not a
   circulation error (ledger exact). Reported per-σ accordingly.
5. **Attribution is a null on the wing:** :split vs :upstream differ by
   2e−4 pp in ΔCL — three orders below the σ effect. First datapoint for the
   pre-registered rotor T4 null.
6. **Induced velocity at fixed probes** (`run_ssw_conversion_probes`, matched
   t*=5/dt*=0.25 pair, one DirectBackend end-state evaluation over panels +
   filament + particles at 18 points): attached near-wake (x/c=1.5) agrees to
   ~1% of legacy RMS; deltas rise toward the starting-vortex wake front
   (max 14% at x/c=5, y≈0.6·semispan) where velocities are ~3× smaller and the
   two deposits legitimately differ at particle-discretization scale
   (absolute deltas ≤8e−3·U∞). Snapshot, not per-step history.

**§3 case 1 verdict:** the smooth strategy is aerodynamically equivalent to
legacy on the wing at matched resolution (deltas ~0.1%, σ-attributable), with
exactly conserved transfer and no stability penalty. The wing gate for
proceeding to the rotor (§ case order) is satisfied. Rotor hover (incl. T4)
remains deferred to a follow-up session.

## 2026-08-04 — §3 case 2 rotor hover A/B: SUBMITTED (no results yet)

Ryan asked for the pair to be launched. **Nothing below is a physics result** —
this section records the launch and its pre-flight verification only.

| arm | case tag | job | conversion |
| --- | --- | --- | --- |
| A (control) | `p018_conv_legacy` | 13036850 | `legacy` edge-jump |
| B (test) | `p018_conv_smooth` | 13036851 | `SurfaceVorticityConversion` |

Both `--time=48:00:00`, 64 cores / 64G, started 2026-08-04 10:08:07.

**Carrier = the item-018 B0 baseline** (`p018_*` block of
`examples/run_dji9443_hover_ct_hpc.slurm.sh`), so the two arms share every knob
by construction rather than by matching two hand-written configs: mesh
`45_185_ct4`, velocity formulation, RPM 5400, NT=36, depth 4R, rlxf 0.3,
η=1.0 (Das\* = 0.41c), nwakerows=4, overlap 2.0, pps 4, merge_r 0.0120,
settle 12, filter off, `das_arc=false`, inviscid. The new arms set **only**
`CONVERSION`, by explicit assignment rather than `${VAR:-default}` — the A/B's
defining knob must never be winnable by an inherited environment value (the
`p018_nt72` η leak, 2026-08-03).

σ is matched, not merely similar: the smooth arm's `CONVERSION_SIGMA` is left
unset, so the driver resolves it to `tip_sigma_default = 2πR/NT · OVERLAP/PPS`
— the same σ the `sigma_overlap` policy sheds at — exactly as in the wing A/B.
`ATTRIBUTION` is `:upstream` in both arms; the pre-registered T4 null
(`:upstream` vs `:split`) is a **separate** later pair and is deliberately not
folded into this one.

The legacy arm is run rather than reusing B0's stored 0.07170: `src/` has moved
since B0 (TE-selection hardening, CoreSpreading integration fix), and although
both were verified bit-identical for stock settings, an A/B whose control was
produced by a different code state is not a controlled A/B. Reproducing B0 is
therefore itself a check on those two changes.

**Pre-flight (per `018/ops_reference.md`):**

- Cluster `src/` md5-verified against local — **all 24 files identical**, as are
  the driver and (after deploy) the launcher. Note the cluster is at
  `5615ada` + scp'd working files, so the hash check is the only real
  provenance; the batched comparison first showed a phantom missing entry
  (documented MOTD gotcha) and was redone with a grep sentinel plus an
  individual check.
- `logs/slurm` present; no pre-existing `data/p018_conv_*` (no wipe); disk
  130G of the 200G budget, 1.9T free.
- Launcher banner extended with a `conversion:` line so the A/B's defining knob
  is auditable within seconds of start, not only from the much-later
  `_case_metadata.toml`.
- **Post-submission knob verification PASSED** — both banners read identically
  on mesh/NT/η/overlap/pps/merge_r/rlxf/nwakerows/das_arc/visc and differ only
  in `conversion:legacy` vs `conversion:smooth`.

Particle-count blowup was estimated before submission rather than discovered:
σ = 2πR/36 · 2.0/4 ≈ 1.04e−2 m, h = σ/1.3 ≈ 8e−3 m against a 0.889R shed span
and a Das ≈ 0.41c streamwise extent gives O(10) particles per blade per step,
comparable to legacy's σ-spaced edge shedding — far from the 500k cap.

### 2026-08-04 (b) — two PRIOR completed smooth rotor runs found; T4 ANSWERED

These are the arms submitted by the **2026-08-03 session 2** (jobs 13035886 =
`p016_smooth_up`, 13035887 = `p016_smooth_split`; see `log.md`, which records
the submission, the verified banners, and a harvest TODO). This session
rediscovered them on disk while inventorying, having read INDEX + this phase
doc but not `log.md` — an avoidable detour; **`log.md` is the item's session
narrative and must be read before concluding anything about run inventory.**
Both COMPLETE (720 steps ≈ 20 revs, 8.3–8.6 h, finished 2026-08-04 07:04 /
07:22). A full `_case_metadata.toml` diff shows
`p016_smooth_up` is identical to `p018_b0` in **every** field except the
conversion block, so these are the T4 pair on the B0 carrier at
σ_conv = 0.0103847 (= the `sigma_overlap` σ), η=1.0, nwakerows=4.

Reduced with the campaign's own `scripts/p018_analyze.py` over revs 10–20
(these are 20-rev runs; the campaign's 15-rev settling criterion is therefore
**not** met — all three report `M1 settled CHECK (< 15 revs)`).

| run | conversion | attribution | CT̄ (revs 10–20) | 95% CI |
| --- | --- | --- | --- | --- |
| `p018_b0` | legacy | — | 0.071866 | [0.071676, 0.072070] |
| `p016_smooth_up` | smooth | `:upstream` | 0.072214 | [0.072051, 0.072408] |
| `p016_smooth_split` | smooth | `:split` | 0.072220 | [0.072068, 0.072396] |

**T4 attribution null: PASSES on BOTH observables.** ΔCT(`:split`−`:upstream`)
= 5.7e−6 = **0.008%**, ~9× below the pre-registered threshold (≈5e−5, 0.07%).
Crucially the M2 circulation check — the observable ruling 9 exists to catch —
agrees rather than hiding a redistribution: **ε_Γ max 0.028% / RMS 0.013%**,
36× below the 1% criterion. Ledgers exactly conserved in both
(`residual_norm` 7.5e−18 `:split`, 2.7e−17 `:upstream`, over 715 conversions
and ~115k particles). So attribution is a **free parameter of measured size
~6e−6 in CT and ~3e−4 in Γ̄** on the rotor, to be carried at that magnitude in
the error budget. Note what this does *not* say: it does not favour
`:upstream`. It says the choice is immaterial at this resolution, so the
refinement-order argument (favouring `:split`, Phase 1 §9b.4) remains the only
axis with a resolved sign.

**Smooth vs legacy does NOT reproduce the wing's equivalence.** Against
`p018_b0`: CT̄ **+0.48%** (3.5e−4, about 2× the per-run CI half-width — resolved
but only just, and the two CIs abut at 0.072051 vs 0.072070), and
**ε_Γ max 1.049% / RMS 0.540% ⇒ M2 FAILS**, a hair over the 1% criterion.
On the wing the two strategies agreed to ~0.1% in CL *and* in Γ_TE; on the
rotor the circulation distribution moves ~1%. This is consistent with
014/017's conclusion that **σ/chord is the bridge variable** — rotor σ ≈ 1.5
local chords versus the wing A/B's σ/c ≈ 0.04–0.16 — i.e. the wing gate passing
was necessary but did not transfer, exactly the risk the case ordering was
meant to expose.

**Γ̄(r/R) plotted** — `scripts/p016_plot_conversion_gamma.py` (reuses the
018 loaders; categorical palette + line style because the two smooth arms
nearly coincide, which IS the T4 result). Figure and table:
`data/p016_conversion_gamma/gamma_conversion.{png,csv}`.

The shape matters, and it is **not** the σ-axis pattern:

- The delta is **one-signed and inboard-peaked**, not a redistribution. It
  rises from +0.06% at the tip to **+1.34% at r/R ≈ 0.18**, decaying
  monotonically outboard through +1.05% at the band edge (0.31), +0.44% at
  0.56, +0.16% at 0.94. Smooth adds circulation everywhere, most strongly
  inboard.
- Consequently **CT̄ and Γ̄ agree here**: the band integral moves **+0.506%**
  against ΔC̄T +0.48%. Contrast the σ-axis, where ∫Γ̄dr moved only +0.117%
  while ε_Γ hit 8.78% — that was the case CT̄ was blind to. This one CT̄ sees.
- **ε_Γ = 1.05% is clipped by the metric band.** The maximum sits at
  r/R ≈ 0.18, *outside* 0.3 ≤ r/R ≤ 0.95, so M2 as defined understates this
  failure mode; the true worst is 1.34%. Worth noting whenever a failure's
  worst excursion lands at a band edge.
- **The inboard-peaked signature matches 018's `d/σ` clearance family** —
  N=1 and small-`Das` both fail monotonically inboard, because `d = N·Das ∝ r`
  makes clearance intrinsically worst at small radius. A conversion that
  changes where circulation is deposited relative to the surface would be
  expected to show exactly this radial signature. Hypothesis, not a
  conclusion: it predicts the delta should shrink at larger `d/σ`, which the
  in-flight `Das`×N matrix could test.

**Against the pre-registered A/B band** recorded in `log.md`
(|ΔC̄T| ≤ 0.5%, ε_Γ,max ≤ 1%): CT **passes at 0.48%** and ε_Γ **fails at
1.049%** — both marginal, on opposite sides. The verdict therefore rests on
the observable the campaign already ruled is the binding one (018 ruling 9,
M2), so the honest reading is a **marginal FAIL pending the unconfounded
re-measurement** below, not a pass.

**Confounded, and this is why the in-flight pair matters:** `p018_b0`'s
metadata has no `conversion` field at all, proving it predates the 016 driver
wiring, so smooth-vs-legacy above straddles a `src` change. Jobs 13036850 /
13036851 (previous section) are the same-code-state pair that removes the
confound; the ε_Γ 1.049% must be re-measured against `p018_conv_legacy` before
it is treated as a property of the conversion rather than of the code state.
Ryan's ruling (2026-08-04): keep both arms running for exactly this reason.

Harvested locally to `data/p016_smooth_{up,split}/` (CT series, per-rev,
metadata, conversion ledger, bound-circulation monitor).

**Open until harvest:** settled CT̄ per arm over the campaign's 10-rev window,
Γ̄(r/R) via `scripts/p018_plot_gamma.py` (the σ-axis lesson from 018 applies —
CT̄ alone can read as agreement while the circulation distribution disagrees, so
ε_Γ is a required second observable here, not optional), the smooth arm's
external filament-ledger residual from `*_conversion_diagnostics.csv`, particle
count and cost, and whether the legacy arm reproduces B0 = 0.07170.

### 2026-08-04 (c) — unconfounded pair harvested: §3 case 2 verdict is FAIL, and it is a property of the conversion

Jobs 13036850 (`p018_conv_legacy`, 12h03) and 13036851 (`p018_conv_smooth`,
9h31) both COMPLETED exit 0, 720/720 steps. Pre-trust checks passed: the two
`*_case_metadata.toml` differ **only** in derived CT statistics (outputs) —
every knob identical, `conversion` the sole input difference, confirming the
submission-time banner audit. Harvested to `data/p018_conv_{legacy,smooth}/`
(CT series, per-rev, metadata, bound-circulation + force monitors, smooth
conversion ledger, both slurm logs). Reduction via `scripts/p018_analyze.py`,
matched windows revs 10–20.

| observable | legacy 13036850 | smooth 13036851 | Δ (smo−leg) | criterion | verdict |
| --- | --- | --- | --- | --- | --- |
| C̄T (revs 10–20) | 0.071847 [0.071659, 0.072063] | 0.072210 [0.072050, 0.072391] | **+0.505%** | ≤0.5% | **FAIL (marginal)** |
| ε_Γ max, 0.3≤r/R≤0.95 | — | — | **1.080%** (at r/R 0.313) | ≤1% | **FAIL** |
| ε_Γ RMS, in-band | — | — | 0.567% | (reported) | — |
| ε_Γ max, unclipped | — | — | 1.344% at r/R 0.178 | (reported) | — |
| ledger residual | n/a (legacy emits none) | 1.650e−17 (715 conversions, 115,051 particles) | — | exact | PASS |

**The confound is resolved against the conversion.** The unconfounded pair
reproduces the confounded numbers almost exactly (ΔC̄T +0.48%→+0.505%, ε_Γ
1.049%→1.080%): the ~1% circulation delta is a property of the smooth
conversion at rotor σ/chord, not of the code state. Simultaneously the
control closes the code-state question: `p018_conv_legacy` reproduces
`p018_b0` to **−0.026%** on the same window (0.071847 vs 0.071866, well
inside both CIs) — the TE-selection hardening + CoreSpreading integration fix
are CT-neutral, retroactively validating the campaign baseline.

**Robustness.** Window shift to revs 12–20 moves the verdict *further* into
FAIL (ΔC̄T +0.675%, ε_Γ max 1.570%) — legacy is still drifting down (5-rev
block drift 0.501%, monotone) while smooth is flat (0.100%), so the 20-rev
settling caveat (`M1 settled CHECK`, <15 settled revs) cannot rescue the
comparison; if anything the FAIL is conservative. Both arms carry the caveat.

**Structure.** Γ̄(r/R) delta is **one-signed** (smooth > legacy at every
station with |δ|>0.2%), inboard-peaked (max at r/R≈0.18, decaying to ~0 at
the tip), span integral +0.772% — the same radial signature as 018's small-`d`
family, consistent with the pre-registered d/σ-clearance hypothesis (§ above):
the conversion changes where circulation is deposited relative to the surface,
and clearance is intrinsically worst inboard. Prediction unchanged: the delta
should shrink at larger `d/σ` (testable against the 018 `Das`×N matrix).
Figure + per-station CSV: `data/p016_conversion_gamma/gamma_conv_pair.{png,csv}`
(delta computed with `p018_analyze.py`'s own `gamma_profile`/`_interp`, so it
reproduces `m2` bit-for-bit; in-band numbers match the `m2` subcommand).

**§3 case 2 verdict: smooth-vs-legacy equivalence FAILS at rotor σ/chord**
(σ_conv = 0.0104 ≈ 1.5 local chords) — marginally on C̄T, clearly on the
binding observable ε_Γ, and now unconfounded. The wing's ~0.1% equivalence
(§3 case 1) did not transfer, which is the outcome the case ordering was
designed to detect (014/017 σ/chord bridge). This is a valid negative result
per §3's pre-registration; it does not block the smooth conversion as an
*opt-in* feature (ledger exact, runs stable, cost comparable — legacy's longer
walltime is node variation, same step count), but it does mean smooth is NOT a
drop-in replacement at rotor resolution, and 016's aerodynamic equivalence
claim is scoped to wing-like σ/chord until a larger-d/σ rotor test passes.


### 2026-08-05 — particle count and cost (closing the last open harvest item)

Deposition rate: the smooth ledger records **115,051 particles over 715
conversion events ≈ 161 particles/step**. Legacy emits no ledger, so its
counts come from wake-population VTK. `p018_conv_legacy`'s own VTK was
deleted (Ryan-approved, ledger 2026-08-05), so the legacy column below is
**`p018_dasc0p41`** — still writing on the cluster and banner-verified to be
the exact B0-carrier legacy configuration at the same code state (das_chord
0.41, η 1.0, nwakerows 4, pps 4, overlap 2.0, merge_r 0.0120,
conversion:legacy). Smooth column is `p016_smooth_up` (full VTK preserved on
Ryan's Mac, `~/p016_smooth_vtk_paraview/`), same carrier, `:upstream`.

| step (rev) | legacy wake particles | smooth wake particles | legacy/smooth |
| --- | --- | --- | --- |
| 36 (1) | 4,598 | 3,500 | 1.31 |
| 72 (2) | 12,465 | 9,188 | 1.36 |
| 180 (5) | 36,760 | 26,790 | 1.37 |
| 360 (10) | 70,159 | 50,950 | 1.38 |
| 540 (15) | not yet reached | 59,350 | — |
| 719 (20) | not yet reached | 58,663 (peak 60,249 @ 655) | — |

**Smooth runs the wake with ~27% fewer particles at matched steps** (ratio
stable at ~1.37 once the column develops); smooth plateaus at ~59–60k under
the 4R truncation + merging. Consistent with the smooth arm's shorter
walltime (9h31 vs legacy 12h03), so that difference is likely not purely node
variation. Caveats: the legacy column is a same-code-state proxy run, not job
13036850 itself; its 540/720 rows can be filled when `p018_dasc0p41`
(job 13050925) finishes. Cost-per-particle-step is therefore comparable or better for
smooth — the equivalence FAIL above is a physics/accuracy scoping, not a cost
penalty.

**Larger-d/σ discriminator SUBMITTED (Ryan-approved 2026-08-05):** new
launcher arm `p018_conv_smooth_dasc0p82` (unconditional
`DAS_CHORD_FRACTION=0.82` + `CONVERSION=smooth`; launcher deployed
md5-verified; cluster `src/` untouched ⇒ same code state as the conv pair) as
job **13051758**, pairing with the in-flight legacy `p018_dasc0p82`
(13050926). Both banners verified: they differ only in `conversion:`. This
doubles d/σ at unchanged σ. Pre-registered readout: ε_Γ(smooth vs legacy) at
Das=0.82c vs the 1.080% at 0.41c — shrinks ⇒ clearance hypothesis confirmed
(FAIL is operating-point-specific; smooth adoptable at converged
configurations); persists ⇒ conversion-operator difference. Secondary
readout: whichever arm moves less between Das 0.41c→0.82c was closer to the
converged limit at small d/σ.

**Equal-h discriminator SUBMITTED (Ryan-approved 2026-08-05, prompted by an
outside-agent observation):** the conv A/B was not sampling-density matched —
legacy sheds at h = σ/2.0 (`OVERLAP=2.0`) while the smooth quadrature samples
at h = σ/1.3 (`CONVERSION_OVERLAP` default), a ~1.54× mismatch consistent
with the measured ~1.37 population ratio, and exactly the "apparent smoothing
from coarser effective resolution" confound §3 warns about. The wing A/B
passed at ~0.1% with this same mismatch, so h can only matter where σ/chord
is large. New arm `p018_conv_smooth_h2p0` (unconditional `CONVERSION=smooth`
+ `CONVERSION_OVERLAP=2.0`; σ, Das, carrier all matched) submitted as job
**13051763**, pairing with the completed `p018_conv_legacy` (13036850).
Readout: ε_Γ collapses from 1.080% ⇒ the FAIL is a sampling-density artifact
(and smooth's particle savings were partly artifact); unchanged ⇒ h is inert
at rotor scale and the d/σ / conversion-operator candidates remain. The two
discriminators are orthogonal (conv_overlap constant within the Das pair; Das
constant within the h pair), so their four joint outcomes cleanly separate
clearance, sampling, and operator explanations.
