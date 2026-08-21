# 016 Session Log

One dated entry per working session, **newest first**. This file is the append
target for all session narrative on this item: approval events, review findings,
bugs found and fixed, and measurement outcomes. The item file
(`../016_surface_vorticity_particle_shedding.md`) carries only objective, phase
gates, invariants, and current status; per-phase derivation and evidence stay in
that phase's document in this folder.

## Dated entries

- **2026-08-05 — particle-count/cost measurement recorded; larger-d/σ smooth
  arm SUBMITTED (Ryan-approved).** Cost item closed in phase_04 §2026-08-05:
  smooth deposits 161 particles/step (ledger) and runs the wake with **~27%
  fewer particles than legacy** at matched steps (ratio ~1.37 stable; legacy
  column = knob-verified same-code-state proxy `p018_dasc0p41`, since
  13036850's VTK was deleted — Ryan-approved deletion of all four conv/smooth
  run VTKs logged in 018 ledger; `p016_smooth_{up,split}` survive byte-exact
  on Ryan's Mac). Follow-up experiment: Ryan approved the pre-registered
  larger-d/σ discriminator — new launcher arm `p018_conv_smooth_dasc0p82`
  (DAS_CHORD_FRACTION=0.82 + CONVERSION=smooth, unconditional exports per
  ops rule; launcher deployed md5-verified, cluster `src/` untouched so the
  run shares the conv pair's code state) submitted as job **13051758**
  (24 h / 64c / 64G). Legacy partner = in-flight `p018_dasc0p82` (13050926).
  This doubles d/σ vs the conv A/B at unchanged σ: delta shrinks ⇒ clearance
  hypothesis confirmed and the FAIL is operating-point-specific; delta
  persists ⇒ conversion-operator difference (routes to internals, not 012).
  Whichever arm moves less toward the common large-d/σ limit was the more
  accurate one at small d/σ. Banner verified: `das_chord:0.82
  conversion:smooth`, carrier stock; legacy partner banner also verified —
  the pair differs only in `conversion:`. **Second discriminator (same day,
  Ryan-approved after an outside-agent observation): equal-h arm
  `p018_conv_smooth_h2p0` (job 13051763)** — the conv A/B compared legacy
  shedding at h=σ/2.0 against smooth quadrature at h=σ/1.3 (the
  `CONVERSION_OVERLAP` default), a ~1.54× sampling-density mismatch that is a
  candidate explanation for the ε_Γ FAIL (and for smooth's ~27% particle
  savings). This arm sets `CONVERSION_OVERLAP=2.0` so smooth samples at the
  same h legacy sheds at; pairs with the COMPLETED `p018_conv_legacy`
  (13036850); only sampling density differs from `p018_conv_smooth`. Readout:
  ε_Γ collapses from 1.080% ⇒ h-mismatch explains the FAIL (and the particle
  savings were partly artifact); unchanged ⇒ h inert at rotor scale, d/σ /
  operator explanations remain. Note the two discriminators are orthogonal:
  conv_overlap is constant within the Das pair, Das constant within the h
  pair. Wing context: the wing A/B achieved ~0.1% equivalence WITH the same
  h mismatch, so h can only matter where σ/chord is large.

- **2026-08-04 (b) — Phase 4 session 4: unconfounded pair harvested; §3 case 2
  = FAIL, attributed to the conversion.** Both arms COMPLETED (13036850 legacy
  12h03 / 13036851 smooth 9h31, 720 steps, exit 0); metadata knob-diff clean
  (only `conversion` differs). Matched windows revs 10–20
  (`scripts/p018_analyze.py`): legacy C̄T 0.071847 [0.071659, 0.072063], smooth
  0.072210 [0.072050, 0.072391] ⇒ **ΔC̄T +0.505% (marginal FAIL vs ≤0.5%)**;
  **ε_Γ max 1.080% in-band FAILS ≤1%** (RMS 0.567%; unclipped worst 1.344% at
  r/R 0.178). The pair reproduces the confounded numbers (+0.48%/1.049%) ⇒
  the delta is a **property of the smooth conversion at rotor σ/chord**, not
  the code state. Control: legacy reproduces `p018_b0` to −0.026% ⇒
  TE-selection hardening + CoreSpreading fix are CT-neutral (campaign baseline
  retroactively validated). Robust: revs 12–20 window makes it worse (+0.675%,
  1.570%) — legacy still drifting (block drift 0.50% monotone), smooth flat
  (0.10%); settled caveat stands but cannot rescue. Delta one-signed,
  inboard-peaked, span integral +0.772% — d/σ-clearance signature; prediction
  (shrinks at larger d/σ) testable on 018's Das×N matrix. Smooth ledger exact
  (1.65e−17; 715 conversions / 115,051 particles). Figure/CSV
  `data/p016_conversion_gamma/gamma_conv_pair.{png,csv}`. Ops: p016_smooth_up
  + p016_smooth_split fully copied byte-verified to Ryan's Mac
  (`~/p016_smooth_vtk_paraview/`, 7,919 files / ~8.4G each) ahead of disk
  pressure (home crossed 190G of the 200G budget); cluster-side VTK deletion
  requires Ryan (agent-blocked) — one-liner in session transcript. Evidence:
  Phase 4 doc §"2026-08-04 (c)".

- **2026-08-04 — Phase 4 session 3: T4 ANSWERED; smooth-vs-legacy marginal FAIL
  (confounded); unconfounded pair in flight.** Harvested the session-2 arms
  (13035886 `p016_smooth_up`, 13035887 `p016_smooth_split`; both COMPLETE, 720
  steps, 8.3–8.6 h) and reduced them with the campaign's own
  `scripts/p018_analyze.py` over revs 10–20. **M1:** legacy `p018_b0` 0.071866
  [0.071676, 0.072070]; smooth `:upstream` 0.072214 [0.072051, 0.072408];
  smooth `:split` 0.072220 [0.072068, 0.072396]. All three report
  `M1 settled CHECK (< 15 revs)` — these are 20-rev runs, so the campaign
  settling criterion is NOT met and every number here inherits that caveat.
  **T4 PASSES on both observables:** ΔCT(`:split`−`:upstream`) = 5.7e−6
  (0.008%), ~9× under the pre-registered 5e−5 threshold, and — the part that
  matters under 018 ruling 9 — M2 agrees rather than hiding a redistribution,
  ε_Γ max 0.028% / RMS 0.013%, 36× under the 1% criterion. Ledgers exactly
  conserved (residual_norm 7.5e−18 `:split`, 2.7e−17 `:upstream`; 715
  conversions, ~115k particles). Attribution is therefore a **free parameter of
  measured size ~6e−6 CT / ~3e−4 Γ̄**, to be carried at that magnitude in the
  error budget. It does NOT vindicate `:upstream`: it makes the choice
  immaterial at this resolution, leaving refinement order (favouring `:split`,
  Phase 1 §9b.4) the only axis with a resolved sign — so the Stage 4b/R3
  question "why is `:upstream` the default?" is answered with "it does not
  matter", not with evidence for it. **Smooth vs legacy: ΔC̄T +0.48% (passes
  the ≤0.5% band) but ε_Γ max 1.049% / RMS 0.540% (FAILS ≤1%)** — marginal on
  both, opposite sides, and the binding observable is the failing one. The
  wing's ~0.1% equivalence did NOT transfer to the rotor, consistent with
  014/017's σ/chord bridge variable (rotor σ ≈ 1.5 local chords vs the wing
  A/B's σ/c 0.04–0.16) — i.e. the case ordering worked as designed.
  **Confound:** `p018_b0` has no `conversion` field in its metadata, proving it
  predates the 016 driver wiring, so this comparison straddles a `src` change.
  Ryan's session-2 decision had been "B0 reused as the legacy arm (no fresh
  legacy run)"; **reversed by Ryan 2026-08-04** on that evidence — a fresh
  same-code-state pair was submitted on new launcher arms `p018_conv_legacy` /
  `p018_conv_smooth` (jobs **13036850** / **13036851**, ~9 h, banners verified
  identical on every knob except `conversion`; launcher banner gained a
  `conversion:` line so the defining knob is auditable in seconds). ε_Γ 1.049%
  must be re-measured against `p018_conv_legacy` before it is attributed to the
  conversion rather than to the code state. **Process:** this session
  rediscovered the session-2 runs by inventorying cluster disk, having read
  INDEX + the Phase 4 doc but not `log.md` — read `log.md` first when the
  question is "what has been run". Also: 31 GB of VTK reclaimed on the cluster
  (`p018_L1_ov3` copied byte-verified to Ryan's Mac first); disk watchdog armed;
  detail in `018/ledger.md`. Evidence: Phase 4 doc §"2026-08-04 (b)".

- **2026-08-03 — Phase 4 session 2: rotor A/B + T4 arms SUBMITTED (in flight).**
  Pre-submission work: (1) conversion-ledger emission added to
  `examples/rotor_hover_pressure_comparison.jl` (sibling
  `<run>_conversion_diagnostics.csv` written after `simulate!`, smooth runs
  only — legacy emits nothing and stays bit-identical; on a restart-chained run
  the ledger covers only the final segment). Validated locally: coarse smoke
  (40_40, NT=6) smooth residual_norm 2.7e−17 with 52 conversions / 3848
  particles, legacy no file, both all-finite; wake + simulate unit suites 0
  Fail/Error. (2) Cluster sync restored: driver +
  `src/FLOWPanel_{metadata,replay,simulate,wake,warmstart}.jl` deployed to
  `/home/rander39/projects/FLOWPanel.jl` and md5-verified (launcher already
  matched). Ryan's decisions this session: B0 reused as the legacy arm (no
  fresh legacy run), agent submits over the socket, ledger instrumentation
  added first, the 2 reserved 018 slots used now. Submitted on the `p018_b0`
  carrier: job **13035886** = `p016_smooth_up` (`ATTRIBUTION=upstream`), job
  **13035887** = `p016_smooth_split` (`ATTRIBUTION=split`); distinct
  `P018_RUN_NAME` per arm (launcher wipes `data/$RUN_NAME`). Knobs verified
  from both log banners per the 018 mandatory rule: exact B0 carrier
  (45_185_ct4 captess4 mesh, NT=36, das_eta 1.0, nwakerows 4, overlap 2.0,
  pps 4, merge_r 0.0120, depth 4R, rlxf 0.3, filter off, 18.5 revs / 720
  steps), `CONVERSION=smooth`, per-arm attribution correct, and
  `CONVERSION_SIGMA=0.0103847` = 2πR/NT·OVERLAP/PPS (the ladder σ ⇒ matched-σ
  A/B). ~24 h wall each. **Harvest session TODO:** health check (log content,
  not exit codes; OOM-in-merge = divergence), harvest CT CSVs + metadata +
  ledger + bound-circulation monitor CSVs, then M1 (`scripts/p018_analyze.py`
  vs B0 and up-vs-split), M2 Γ̄(r/R) (`p018_analyze.py` + `p018_plot_gamma.py`),
  ledger residual, particle-count/cost table; A/B verdict vs the 018 band
  (|ΔC̄T|≤0.5%, ε_Γ,max≤1%), T4 verdict vs the pre-registered 5e−5 CT
  threshold; then Phase 4 doc §3-case-2 evidence section + item/INDEX updates.
  Plan: `~/.claude/plans/finish-item-016-per-quiet-penguin.md`.

- **2026-08-03 — Phase 4 session 1 COMPLETE: §1/§2 closed, §3 wing A/B PASSED;
  rotor remains.** §1/§2: coverage map recorded in the Phase 4 doc; five
  gap-closure tests added (ghost/stale-row non-influence + source count, which
  also pins R4's `overflowed[]` gate at source level; resolution-floor centroid
  placement; κ convergence on nonuniform planar and warped sheets, measured
  orders 1.62/1.39 and 1.59/1.45; end-to-end affine-μ exact deposition; a
  time-marched smooth smoke in the simulate suite — closing the "no end-to-end
  test" caveat). Wake 605/605, simulate 88/88, replay 125/125, warm start
  26/26, full suite exit 0 / zero Fail-Error. §3 case 1: `conversion=:smooth`
  threaded into `examples/suddenly_started_wing.jl` (legacy tags byte-identical;
  smooth appends `_svc_…`; ledger totals → `conversion_diagnostics.csv`) and a
  new pre-registered driver `examples/ssw_conversion_ab.jl` ran an 8-pair
  matrix (dt* ladder × strategy, σ and overlap arms, attribution arm) plus an
  end-state induced-velocity probe pair. **Result: smooth ≡ legacy to ~0.1% in
  CL/Γ/loading, every pair settled, ledger residual ≤2e−16; residual delta is
  σ-linked (linear in conversion σ, inert to overlap); attribution :split vs
  :upstream is a 2e−4 pp null on the wing (first T4 datapoint); probes agree
  ~1% near-field, degrading only at the low-velocity starting-vortex front.**
  Backend note (Ryan mid-session): driver default switched from the
  template-inherited DirectBackend to FMM; the dt*=0.25 pair reproduces ΔCL to
  4 digits across backends. Wing gate for the rotor is satisfied; rotor A/B +
  pre-registered T4 deferred to a follow-up session. Detail in the Phase 4
  doc's two 2026-08-03 evidence sections; data in `data/ssw_conversion_ab/`.

- **2026-08-03 — Phase 3 APPROVED; Phase 4 opened.** Ryan explicitly approved
  Phase 3 in-session, following the clear-context re-review recommendation.
  Session scope for Phase 4: §1 mathematical verification + §2 mechanical
  regression (record the existing coverage map, close the identified gaps:
  ghost non-influence/active-source-count, resolution-floor centroid placement,
  panel-refinement convergence orders on nonuniform/warped grids, end-to-end
  affine-μ conservation through `_convert_to_particles!`, smooth smoke in the
  simulate suite) and §3 case 1, a suddenly-started-wing legacy-vs-smooth A/B
  (new `examples/ssw_conversion_ab.jl`; `conversion` option threaded into
  `examples/suddenly_started_wing.jl` with an append-only case-tag token).
  §3 case 2 (rotor hover, incl. pre-registered T4) is deferred to a follow-up
  session per the phase's case ordering. Item-level INDEX columns remain
  unchecked until all of Phase 4 completes.

- **2026-08-03 — Clear-context re-review of Phase 3 (post-remediation):
  approval RECOMMENDED; user checkbox pending.** Fresh-context verification of
  the R1/R2 guards, the source-derived (`get_sources`) transfer ledger, the
  replay rejection path, Stage 4c pre-registration compliance (incl.
  `overflowed[]` set, Richardson 256/512 reference, paired statistics,
  `alpha`-affinity check), and every Stage 4c number against the checked-in
  CSVs — all hold. The near-field designation pre-dates the data (Phase 1
  §5.7 + Ryan's 2026-08-01 charge), so `:upstream`-on-evidence stands. Suites
  re-run fresh: wake 583/583, replay 125/125, warm start 26/26, simulate
  81/81, exit 0, no `Fail`/`Error`. One non-blocking precision note logged in
  the Phase 3 doc: T1's "12/12" unanimity is scoped to the handoff family
  (pooled near-field is 16/32 because mid-chord `:control` tracks the
  far-field branch, as disclosed). No blockers found.

- **2026-08-03 — Stage 4c complete: R3 CLOSED, `:upstream` retained ON
  EVIDENCE; new finding R4 against the Stage 4b evidence.** A pre-registered,
  simulation-free criterion (`generate_attribution.jl`; all rules committed to
  the Phase 3 document before the first run) replaced Stage 4b's
  compare-to-the-artifact metric with a refined reference. **T1** is conclusive
  (Richardson reference, order 1.99, residual 65x below the signal) and
  regime-dependent with a clean crossover at `d/sigma ~ 3`: `:upstream` lower
  error at 12/12 probes for every standoff `<= 2`, `:split` at 12/12 for every
  standoff `>= 4`. **T2** finds no resolved stretching hazard, so no veto.
  **T0** shows the `alpha` spread does **not** collapse in the near field
  (fitted order `-0.46` at `d/sigma = 0.25` rising to `+0.91` at 64) — `alpha` is
  a near-field modelling choice and a far-field truncation artifact, and its
  spread belongs in the CT error budget. **`alpha*`** shows the error-minimizing
  weight is 1.45, *outside* `[0,1]`, so `:upstream` is the best admissible
  choice for pointwise `kappa` (2x better than `:split`), and confirms the
  derived theory result that no admissible `alpha` is consistent on a graded
  mesh. Since the near field is the regime designated in advance (Phase 1 §5.7;
  Ryan's 2026-08-01 charge), **`:upstream` is retained on evidence, not by
  default-preservation**. **R4:** Stage 4b's `probe_errors.csv` was computed with
  the retained filament emitting nothing —
  `get_n_bodies(::FilamentWrapper{<:PanelWake})` is zero unless `overflowed[]` is
  set, and a static fixture never sets it — so that table could not have measured
  the `attribution`-dependent effect it was read as measuring; caveat added to
  the data README. No source file changed; suites unchanged.

- **2026-08-03 — Phase 3 review blockers remediated; technical completion
  restored.** Smooth construction now rejects `unsteady_filament=false` and
  `include_final_filament=false`, the two configurations that invalidate its
  startup and handoff ledgers. Legacy behavior remains unchanged. The external
  transfer test now derives the retained filament from `get_sources`, with
  direct source-level regressions for both failure modes; replay rejects an
  invalid smooth/unsteady manifest. Verification: wake 583/583, replay 125/125,
  warm start 26/26, simulate 81/81, and the full suite through Analytical
  Validation with exit 0 and no `Fail`/`Error` entries. Approval remains
  unchecked; R3's attribution-evidence caveat remains, and Phase 4 is blocked.

- **2026-08-03 — Phase 3 clear-context review: approval WITHHELD.** The ledger,
  the Option B integration, the strict persistence contract, and all verification
  numbers (wake 575/575, replay 124/124, warm start 26/26) were independently
  re-derived or re-run and hold. Two unguarded constructor combinations break
  circulation conservation silently: **R1** `unsteady_filament=false` (a
  supported, unit-tested `PanelParticleWake` option) makes the first conversion
  deposit a starting vortex the sheet does not carry — measured, present `[0,0,0]`
  vs deposited `[0,-0.300,0]`, the mirror image of the 2026-08-01 bug; **R2**
  `include_final_filament=false` removes the filament from `get_sources`, so the
  handoff's cancelling filament is not a source and the deposited face is
  double-counted every step. Both escaped the suite because every smooth test
  uses the defaults, and the external ledger test reads
  `_final_filament_strength` directly instead of asking `get_sources` what the
  wake radiates. **R3** (non-blocking, record integrity): `:upstream` is retained
  by **default-preservation, not by measurement** — no Stage 4b measurement
  favours it; the one axis with a resolved sign (refinement order, 1.98/1.99 vs
  0.98/0.99) favours `:split`; the verdict came from a conjunctive rule written
  after the data was in, which can only ever return the incumbent. The near-field
  axis cannot resolve a sign at all: it scores fidelity to the edge-jump sheet
  being replaced (so it rewards the artefact this item removes), over M=4 points
  on one grid, with normalization that drives the curve shape and a deviation
  level of 0.8-3.3 for every strategy. Correct wording is "retained pending a
  discriminating measurement." Detail and remediation in the Phase 3 log.

- **2026-08-01 — Stage 3 follow-up: startup-edge leak fixed, attribution option
  added.** Ryan's review of the divergence form found a live bug and a free
  parameter.

  **Bug:** the scheme deposits the whole upstream face and none of the
  downstream face, because the retained filament cancels the latter — but only
  once the handoff flag is set. On the *first* conversion the flag is still
  false, so the aft face carried an uncancelled `muhat_G - muhat_D` (the
  starting vortex) that nothing deposited before the row was discarded.
  Measured: `[0.0, -0.300, 0.0]` lost, once. **The existing conservation tests
  could not see it** — they compare the transaction against its own
  `expected_total`, built by the same code that omitted the face. Fixed by
  depositing the aft face whole before the first handoff, where it is the
  sheet's physical trailing boundary rather than an interface with particles.

  **Option:** `attribution = :upstream | :downstream | :split` sets `alpha`, the
  fraction of the upstream face each conversion deposits. Every inter-row face
  is still deposited exactly once across the two conversions that flank it, so
  all three conserve exactly; the retained filament becomes
  `-(alpha*strength[N] + (1-alpha)*strength[N+1])`, reproducing both
  pre-existing forms at its endpoints (`:downstream` is bit-identical to
  legacy). `:split` is **second order** where the others are first (measured
  1.96/1.99 vs 1.04/1.02, 40x better at coarse resolution), because it makes the
  streamwise difference centered — but it retains a half-jump filament at the
  handoff, exactly where fresh particles sit (Phase 1 §5.7), whereas `:upstream`
  leaves nothing at that artificial partition. Default stays `:upstream`;
  **Stage 4b settles it on near-field evidence.** Recorded as Phase 1 §9b.4 and
  the Phase 2 §14 amendment. 549/549 in `Free Wakes`; full gate clean.

- **2026-08-01 — Stage 3 landed; deposition is divergence form.** The smooth
  conversion now deposits: preflight/commit transaction, area particles on a
  bilinear subdivision, true root/tip closure, `particle_handoff_active`
  retargeting `_final_filament_strength`, typed capacity/geometry/state failures,
  and the §8.1/§8.2 diagnostic records. 468/468 in `Free Wakes`; the Stage 0
  legacy golden reference is still bit-exact.

  The first cut used `kappa = -n x grad(muhat)` from the Stage 2 centroid
  reconstruction and needed two patches to work — an estimated upstream centroid
  at `nwakerows == 1`, and a root/tip line extrapolated to the sheet edge to stop
  a one-cell double-count. Review with Ryan established both were symptoms of
  using a gradient reconstruction where the invariant that matters is
  **circulation conservation**, and both were deleted rather than patched. The
  conversion now redistributes the stored ring assembly's own filaments: smear
  the internal partitions into area vorticity, keep the physical boundaries as
  lines. That is exactly conservative on any mesh, requires only the upstream
  *strength* and never its geometry (so `nwakerows == 1` — the production rotor
  configuration — needs nothing invented), and reuses the legacy root/tip
  filaments verbatim. Recorded as Phase 1 §9b and Phase 2 §14.

  Verified before adoption that the two rules converge to the same value:
  identical to 2e-15 on a uniform mesh, both first-order accurate on a graded
  one, difference vanishing at first order and scaling as `(r-1)/2` in the
  row-to-row extent ratio. The trade is recorded rather than glossed — at fixed
  resolution the reconstruction is the better pointwise estimator, while the
  divergence form has exactly zero conservation error where the reconstruction
  leaks `O(r-1)` of circulation every conversion. A per-step leak biases thrust
  cumulatively; a local kappa error does not.

- **2026-08-01 — Phase 2 APPROVED; Phase 3 unblocked and started.**
  Clear-context review checked the architecture against the live source and
  recorded one simplification plus two corrections as Phase 2 §13: the
  persistent ghost row and its ~20-site visibility exclusion are unnecessary
  (conversion runs before the row shift, so the ghost is the pre-shift final
  active row read in place — Option B selected, Option A recorded verbatim);
  wake serializers live in `src/FLOWPanel_replay.jl`, not
  `src/FLOWPanel_metadata.jl`; and the existing unknown-tag `NoShed()` fallback
  in `_deserialize_wake_shedding` must not be reused for conversion decoding.
  User explicitly approved. Phase 3 began with a golden-reference
  characterization test pinning the legacy conversion bit-exactly before any
  refactor (`test/data/legacy_wake_conversion_reference.jl`,
  `make_conversion_fixture` in `test/test_helpers.jl`). Phase 4 remains blocked.

- **2026-07-29 — Phase 2 technically complete; awaiting approval.** The
  architecture fixes the exact legacy-default/smooth-opt-in API, explicit
  active capacity and invisible ghost storage, rank-aware reconstruction,
  physical-area deposition, true root/tip closure, cancelled-edge handoff,
  transactional capacity/failure semantics, concrete diagnostics, and
  metadata/replay/warm-start contracts. The Phase 2 approval checkbox remains
  unchecked and Phase 3 remains blocked pending explicit user approval.

- **2026-07-29 — Phase 1 APPROVED; Phase 2 unblocked.** Clear-context review
  independently re-derived the handoff/row ledgers (7)–(11) and verified all
  seven code-anchored claims (sign convention, vertex order, final-filament
  cancellation, legacy handoff deposit, N+1/N+1 storage, zero-Γ guard, body
  surface-vorticity orientation) at exact source lines; the Alternative B
  rejection is forced by equation (10). User explicitly approved. Non-blocking
  polish notes recorded in the Phase 1 log (sign-chain code anchors, commit
  the uncommitted zero-Γ guard before Phase 3, duplicate (3)/(3a) label).
  Phase 2 may begin; Phases 3–4 remain blocked on their upstream approvals.

- **2026-07-29 — Phase 1 handoff ledger completed; technically complete.**
  The live package orientation gives
  $\boldsymbol H=(\hat\mu_A-\hat\mu_G)
  (\boldsymbol r_R-\boldsymbol r_L)$. The planar-affine streamwise area
  integral equals this vector, internal spanwise edges telescope, and
  repeated conversions leave no artificial row-boundary filament.
  Cancelled-edge Alternative A is selected; retained-filament Alternative B
  is rejected. Phase 1 approval remains unchecked and Phase 2 remains
  blocked pending explicit user approval.
- **2026-07-29 — Item opened and Phase 1 technically completed.** The sign
  convention, finite-difference reconstruction, rank-deficiency policy,
  perimeter/handoff ledger, surface quadrature, limitations, and acceptance
  scenarios are recorded. Phase 2 remains blocked pending explicit user
  approval of Phase 1.
- **2026-07-29 — Phase 1 reopened after user review.** `PanelWake` will be
  planned with an always-allocated non-source ghost panel so
  `nwakerows == 1` has a genuine two-panel one-sided stencil. The ghost is the
  outgoing panel being converted, not a previously converted passive marker.
  Cancelled-edge and retained-filament handoff formulations must now be
  derived and compared before Phase 1 can be technically complete.
