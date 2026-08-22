# 025: Kernel Regularization Update (Compact-Support vs Gaussian, Selectable, New Default)

**Opened:** 2026-08-20 (Ryan directive, staged by agent; spun out of 023/024
findings on radius inflation)
**Current phase:** Phases 0–3 COMPLETE. Default = Gaussian (Ryan convention ruling 2026-08-20); 018 production adopted Gaussian AND the tuned FMM knobs (Ryan rulings 2026-08-21, implemented in the launcher).
**Item-level approvals:**
- [ ] Technical Completion
- [ ] Clear-Context Approval
- [ ] User Approval After AI Discussion

## Current status and next actions

> **RESET BRIEF (e) — 2026-08-21 (Phase 3 DONE; both adoption rulings
> implemented).** Three warm-start continuations of the frozen 018 carrier
> `p018_cs_f1_l3p4` from step 1034, 117 steps each on the mature 181k-particle
> wake, scored from the monitor CSVs
> (`phase_03_018_compatibility.md`): **(1) family null** — gaussian vs vatistas
> (13247863 vs 13247862) differ by **+0.023% in CT̄**, half the parent run's own
> rev-to-rev drift and ~3% of one rev's peak-to-peak; **(2) knob null** —
> gaussian at the tuned knobs vs at production knobs (13290979 vs 13247863)
> differ by **+0.0038%**, worst single step 0.053% of CT, so the certified
> 1.8e-6/1.2e-6 field errors do NOT accumulate over 117 steps; **(3) cost** —
> 160.2 (vat, production knobs) → 140.2 (gau, production knobs) → **71.7 s/step**
> (gau, tuned), i.e. 2.23×, against the parent's own 178.1 s. Wake health
> identical across all three arms (n_particles to 0.1%, same min σ, same max
> Γ/σ², no tripwire). Bonus conditioning evidence: live-wake `max_u` 31.2 (vat)
> vs 24.3 m/s (gau), matching the Phase 0 matched-r_c peak ranking.
> **Implemented:** launcher default `FLOWPANEL_FILAMENT_REG` vatistas→gaussian
> and `FMM_{BODY,WAKE}_*` defaults 17/0.7/109 and 16/0.6/38 (018 launcher only —
> the driver defaults, hence every other study, are untouched); driver gained
> per-pass FMM knob env vars (shared `FMM_*` still the fallback) and now records
> both triples + the family in `case_metadata.toml`; `FLOWPanel.__init__` prints
> the pinned family. Deployed to the cluster md5-verified.
> **Reporting gotcha:** the driver's headline convergence block is meaningless
> for warm-start continuations (it zero-fills restored steps, so it reports
> CYCLE-MEAN 0.0229 ±145%, CONVERGED=false) — score continuations from per-rev
> blocks 30–32 or the monitor CSVs.
> **Deploy gotcha:** do NOT scp local `src/FLOWPanel.jl` to the cluster — the
> local copy includes `gpu_influence`/`gpu_wake`, which do not exist there
> (cost one 90 s job failure, 13290944). Patch the cluster's own copy.
> **Still Ryan's:** all commits (nothing is committed), the ∇φ ≠ u
> potential/velocity unification, 024's census.

> **RESET BRIEF (d) — 2026-08-20 (review fixes, all six findings applied).**
> (1) Gaussian `radius_inflation` is now GRADIENT-AWARE: solves
> `e^(-z)(1+2z)=tol` (fixed point), Δr/rc = 4.99/5.47/5.90 at tol
> 1e-4/1e-5/1e-6 — supersedes the velocity-only `√(2 ln 1/tol)` (≈5.3 at
> 1e-6) quoted in brief (c); Vatistas rule deliberately unchanged
> (legacy-pinned; gradient error ≤1.25·tol documented). (2) `∂c/∂x = −[s]×`
> sign fixed in phase_01 (implementation was already correct). (3) phase_01
> gains the shipped source/doublet derivations (`regularize` enters ONLY the
> doublet-velocity denominator) and the VortexRing potential/velocity
> family-mismatch KNOWN LIMITATION (also in Deferred — open Ryan decision).
> (4) phase_00 metrics relabeled radial |du/dh| + operator norm added (they
> coincide; matched-rc op-norm maxima G 0.50 / V 1.00 / C 2.55 — ranking
> unchanged). (5) matched-Δr bias direction corrected (conservative Vatistas
> radius PENALIZES Vatistas). (6) superseded ruling/gate/contingency
> annotated. Tests: regularization suite 261/261 and fmm suite fully green
> under BOTH the Gaussian default and the vatistas env pin.

> **RESET BRIEF (c) — 2026-08-20 (Ryan convention ruling; default = GAUSSIAN).**
> Ryan's review of Phase 0 exposed that the matched-matching-distance
> convention back-solves rc per family, while the operationally pinned
> quantity is the CORE SIZE (kerneloffset). At MATCHED rc the peaks are
> tol-independent and the decision FLIPS: u_max·rc = 0.45 (Gaussian) vs 0.71
> (Vatistas) vs 1.21 (compact); max|du/dh|·rc² = 0.50 vs 1.00 vs 2.55 —
> Gaussian lowest on both, so per ruling 1 **GAUSSIAN is now the default**
> (`FILAMENT_REGULARIZATION = Ref(GaussianRegularization)`), Δr = rc·√(2 ln
> 1/tol) ≈ 5.3·rc at tol 1e-6 (vs Vatistas 37.6·rc). Regularization suite
> 221/221 under the new default. Benchmarks re-running under Gaussian
> (mature-state A/B vs the Vatistas baselines) per Ryan 2026-08-20.

> **RESET BRIEF (b) — 2026-08-20 (Phases 0–2 executed, subagent).**
> **Phase 0 DECISION: compact-support wins** (peaks under matched matching
> distance, units Γ/2π·Δr: u_max 1.207 vs Gaussian 1.94–2.37 vs Vatistas
> 8.4–26.6; gradient·Δr² 2.55 vs 9.2–13.8 vs 141–1414; compact peaks are
> tol-independent and C¹ at the support boundary) — Gaussian is NOT lower, so
> per ruling 1 + contingency the default is now COMPACT-SUPPORT; Gaussian and
> Vatistas selectable. Docs: `phase_00_peak_comparison.md` (+`phase_00_peaks.jl`),
> `phase_01_theory.md` (all velocity/gradient formulae as (D, κ) pairs bound
> to code; per-family radius_inflation rules).
> **Phase 2 implemented:** `@enum FilamentRegularization` +
> `FILAMENT_REGULARIZATION` Ref + `set_filament_regularization!` (Symbol or
> enum) in `src/FLOWPanel_elements_fmm.jl`; per-family branches in
> `_bound_vortex_velocity`/`_bound_vortex_gradient` (Vatistas path
> byte-equivalent math); family-aware `radius_inflation(VortexRing,...)`;
> env pin `FLOWPANEL_FILAMENT_REG` read in `__init__`
> (`src/FLOWPanel.jl`) so frozen drivers can select without code changes.
> **Tests:** new `test/runtests_unit_regularization.jl` (221 pass: closed-form
> velocity per family, FD-verified gradients, far-field-at-Δr, compact
> exactness beyond rc, C¹ continuity, guards, plumbing) wired into
> `runtests.jl`; fmm/kernel_gradient/kutta suites updated where fixtures
> encoded Vatistas behavior (family-pinned, intent preserved); fmm,
> kernel_gradient, liftingbody, wake, simulate, kutta, kutta_routeb,
> analytical all PASS under the compact default, and the fmm suite also
> passes under the vatistas pin (legacy exactly reproduced). Ran on local
> julia 1.12.5.
> **Known/pre-existing:** `examples/sweptwing.jl` smoke fails in
> GeometricTools loft PLOTTING (PythonCall IndexError) before any solve —
> identical under the vatistas pin, unrelated to this change.
> **NOT done (Ryan gates):** Phase 3 (018 compatibility + warm-start A/B on
> the production carrier), any HPC run, any commit; item-level checkboxes
> untouched. NOTE: the compact default changes FMM radius inflation
> Δr 37.6·rc → rc — 023's cost numbers at Vatistas no longer describe the
> new default; re-measure via the 023 harness after Ryan's review.

> **RESET BRIEF (a) — 2026-08-20 (staging).** Motivation from 023: the
> VortexRing kernel's Vatistas n=2 regularization has an algebraic far-field
> tail (rel. error ≈ ½(rc/h)⁴), so its FMM radius inflation is
> Δr = rc·(2/tol)^{1/4} = 37.6·rc at the default tol 1e-6 — this single term
> forces the body pass ~99% into near-field direct evaluation (measured:
> 45.6 of 45.7 s) and sets its ~36 s/step cost floor. Both alternatives kill
> the tail: the codebase's existing compactly-supported `regularize` (exactly
> singular beyond rc ⇒ Δr = rc, tol-independent) and a Gaussian
> regularization (exponentially small tail ⇒ Δr ≈ rc·√(2·ln(1/tol)) ≈ 5.6·rc
> at 1e-6). **Ryan's decision rule: compare the two candidates' maximum
> velocity and/or velocity gradient (the regularization's job is conditioning
> — lower peak wins); if Gaussian is lower, Gaussian becomes the codebase
> default.** Either way the regularization becomes SELECTABLE so we can A/B
> and roll back; Vatistas stays available as the legacy option.

In-flight anchor (023): jobs 13246526/27 measure the body-pass cost with the
Vatistas Δr shrunk via `FMM_RADIUS_TOL` 1e-5/1e-4 (Δr 0.021/0.0119 m) — an
empirical proxy for how much step time the regularization change releases.

## Objective and scope

Replace the panel-kernel regularization default with the
better-conditioned of {existing compact-support, Gaussian}, implemented as a
**selectable regularization family** (legacy Vatistas retained), with full
derivations, direct/FMM-aligned implementation, and validation. Out of scope:
the FLOWVPM particle kernel (already Gaussian gauserf); the 024 far-field
kernel swap (separate item; its thresholds must be re-derived if this item
changes the default — noted there); changing `kerneloffset` magnitudes.

## Standing rulings (binding)

1. **(Ryan 2026-08-20)** Decision rule: lower max |u| and/or max |∇u| between
   compact-support and Gaussian picks the default; if Gaussian is lower, the
   codebase default CHANGES to Gaussian.
2. **(Ryan 2026-08-20)** Regularization must be user-selectable
   (compare/rollback path); the winner is the new default.
3. ~~Fair-comparison convention to be fixed in Phase 0 before any comparison is
   scored: candidates are compared at MATCHED far-field matching distance
   (same Δr at the same tol), not merely matched rc — peaks are meaningless
   across families otherwise. Convention itself is a Ryan sign-off.~~
   **SUPERSEDED (Ryan 2026-08-20, RESET BRIEF (c)):** the governing
   convention is MATCHED CORE SIZE (kerneloffset is the operationally pinned
   quantity); the matched-Δr comparison is retained in `phase_00` as a
   secondary framing only.
4. Behavioral-alignment constraint (WORKFLOW.md): FMM and direct backends
   must stay aligned for every selectable family; `radius_inflation` rules
   (src/FLOWPanel_elements_fmm.jl:982-991) must be extended per family in the
   same change.

## Phase gates

| Phase | Deliverable | Status | File |
|---|---|---|---|
| 0 | Peak comparison: max \|u\|, max \|∇u\| vs r/rc for compact-support, Gaussian, (Vatistas as reference), at the matched-Δr convention; default decision per ruling 1 | **DONE 2026-08-20 — ~~compact wins~~ superseded → GAUSSIAN wins under the matched-rc ruling (RESET BRIEF (c) + phase_00 addendum)** | `025_kernel_regularization_update/phase_00_peak_comparison.md` |
| 1 | Theory: derive ALL formulae per family — velocity, velocity gradient (and scalar potential where the kernel carries it) for the regularized bound-vortex segment / vortex ring and the source/doublet `regularize` path; matching `radius_inflation(family, rc, tol)` derivations; rendered per Math Syntax rules | **DONE 2026-08-20** | `025_kernel_regularization_update/phase_01_theory.md` |
| 2 | Implementation + testing: selectable family (construction-time knob + env plumbing), new default per Phase 0; unit kernel tests vs Phase 1 closed forms; direct/FMM alignment tests; wing + rotor smokes; solve-conditioning check (iteration counts); mature-state cost/error A/B via `benchmark/p023_fmm_tune.jl` machinery | **DONE 2026-08-20 (local tests; mature-state A/B deferred to Phase 3/Ryan)** | `025_kernel_regularization_update/phase_02_implementation.md` |
| 3 | 018 compatibility task: verify the 018 drivers/launcher against the new default — explicit pinning knob so the frozen RHPC driver's behavior is REPRODUCIBLE both ways (env override only, RHPC stays frozen), one warm-start A/B on the production carrier (CT + wall_s), ledger entry in 018 | **DONE 2026-08-21 — family null +0.023%, knob null +0.0038%, 178→71.7 s/step; both adoption rulings implemented in the 018 launcher** | `025_kernel_regularization_update/phase_03_018_compatibility.md` |

## Contingency chain

- *(historical — fired under the superseded matched-Δr convention, then
  overturned by the matched-rc ruling of RESET BRIEF (c); retained for the
  record)* Phase 0 finds compact-support has the lower peaks → compact-support becomes
  the default for ALL offset-regularized panel kernels (incl. VortexRing);
  the Gaussian machinery is still implemented as a selectable option only if
  its marginal cost is small — else dropped with a note.
- Gaussian's per-pair transcendental cost (exp/erf) measurably erodes the
  near-field savings in the Phase 2 A/B → re-score default on end-to-end
  step cost at fixed achieved error, not on peaks alone (return to Ryan).
- Conditioning regression (solve iterations up, or Kutta/CT shift beyond
  noise in the A/B) → keep Vatistas default, ship selectable families only,
  route findings to 018/021.

## Deferred

- **VortexRing potential/velocity regularization mismatch (open Ryan
  decision, 2026-08-20 review finding 3b):** under the Gaussian default the
  element's velocity is Gaussian-regularized while its scalar-potential
  branch (equivalent constant-doublet panel) remains compact-regularized, so
  ∇φ ≠ u inside the core region for that element. Unification requires a
  Gaussian-regularized doublet-potential derivation. Details:
  `phase_01_theory.md` "Scalar potential — KNOWN LIMITATION".

- 024's far-field-swap thresholds re-derivation under the new default
  (owned by 024, gated on this item's outcome).
- Particle (FLOWVPM) kernel changes — already Gaussian; out of scope.
- kerneloffset magnitude policy (separate 023 lever).
