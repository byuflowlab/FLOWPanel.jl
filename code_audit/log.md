# Unsteady-case code audit — progress log

Objective: (1) identify why the code fails on unsteady validation cases, (2) fix it.
Instructions and full task specs: `code_audit/plan.md`. Protocol: the master agent does
NOT do tasks itself — it launches one fresh-context subagent per task (sequentially).
Each subagent reads `code_audit/plan.md` + this log, does exactly one task, checks the
box, writes findings under the task, and returns. The master agent then launches the
next task's subagent; no context clearing between tasks is needed.

## Tasks

- [x] **Task 0 — Exterior-trace sign check** (monitors:391 `phi -= μ`; analytic + off-surface probe comparison)
- [x] **Task 1 — Steady wake-model consistency** (semiinfinite vs panel vs particle, spanwise lift, simple wing)
- [x] **Task 1b — Diagnose semiinfinite vs finite-wake solve discrepancy** (USER-STAGED PRIORITY; derive/audit the two formulations, verify by clean agent, hunt sign/orientation culprits, manual linear-system check; all later tasks paused)
- [x] **Task 2 — Unsteady panel vs particle wake**
- [x] **Task 3a — Sphere added-mass gate** (analytic, no wake; make permanent CI test)
- [ ] **Task 3b — Wagner/Theodorsen gates** (panel AND particle wake; make permanent CI tests)
- [ ] **Task 4 — Re-evaluate, conclude, plan the fix** (review findings below, name the broken rung + file:line, append fix tasks here)

## Findings

### Task 0

**Verdict: SIGN CORRECT.** The `phi -= μ` conversion at
`src/FLOWPanel_simulate_monitors.jl:391` yields the *exterior* scalar-potential trace.
No code changed.

Setup: pitching_wing.jl wing (n_span=7, n_airfoil=31, n_endcap=5, 644 cells),
steady solve at α=5°, `semiinfinite_wake=true`, monitor-default
`FastMultipoleBackend(10, 0.4, 100)`. Normals verified geometrically outward
(100% of panels have n̂·(cp−center)>0). Control points are exact centroids
(`calc_controlpoints!`, `src/FLOWPanel_abstractbody.jl:772-786`), so the
`_self_limit` path (`src/FLOWPanel_elements_fmm.jl:165-189`) fires at probes.

- **(A) Kernel level** (per-panel `pnl.induced` at cp and cp ± 1e-3·√A·n̂,
  Union{ConstantSource,ConstantDoublet} kernel): median
  (φ(cp+εn̂)−φ(cp))/μ = **−0.998**, (φ(cp−εn̂)−φ(cp))/μ = **−0.0018** across 10
  sample panels. So the self-limit value (+μ/2 per `_self_limit`) is the trace on
  the −n̂ (interior) side, and exterior = interior − μ, exactly as the docstring at
  elements_fmm:126-154 claims.
- **(B) Monitor path** (same `influence!((probes,), (body,), backend; scalar_potential=true)`
  call as `_pressure_bernoulli_phi_dot!`, φ = probe + uinf·x): recipe φ(cp)+uinf·x−μ vs
  exterior probe at cp+1e-3√A·n̂: median|Δ| = **2.15e-05**, max 2.98e-04 (O(offset);
  median|μ| = 0.429, per-μ median 5.6e-05; per-μ max 0.27 is a tiny-μ panel,
  |μ|min = 2e-4, absolute error still O(offset)). Against the interior probe the recipe
  differs by the full μ (median|Δ| = 0.429 ≈ median|μ|; (Δ+μ)/μ median 2.8e-04). Raw
  φ(cp) matches the interior-side probe to median 1.1e-04.

Consequence for the audit: the exterior-trace conversion feeding every φ̇ history is
sound; the unsteady-Bernoulli suspect list narrows to the BE/BDF2 ladder, the ALE
w·∇φ term, wake-model consistency, and unsteady shedding (Tasks 1-3).

Script: `code_audit/scripts/task0_exterior_trace_sign.jl` (runs in ~2 min locally;
prints verdict).

### Task 1

**Verdict: b→a FAILS (~11%, far outside the ~1% gate); c→b PASSES (0.02%).**
The unsteady march (either wake model) settles to a steady state ~11% below the
trusted semiinfinite steady solve; the panel→particle conversion is clean. This
reproduces the reported "marched-to-steady disagrees with steady solver" symptom
with **steady** PressureBernoulli and a scalar-potential panel wake — so it is
independent of the φ̇ machinery. The lever is shedding/first-wake-row (Das)
handling, not wake truncation and not the particle path.

Setup: pitching_wing.jl wing (n_span=15, n_airfoil=31, n_endcap=5, 1124 cells),
α=4.5° fixed, U=100.6 m/s, steady PressureBernoulli(correct_kuttacondition=false)
in all runs, Backslash + FastMultipoleBackend(8, 0.4, 40), dt=0.5 chords/step,
80 steps (wake extent 39.9 chords in both marched runs; particle run: 5-chord
panel buffer + particles to 39.9 chords, np=4340).

- Total CL: a (semiinfinite, steady) = **0.293859**; b (panel wake, settled) =
  **0.262051** (−10.82%); c (particle wake, settled) = **0.261992** (−10.84%;
  c vs b = −0.022%).
- Max sectional deviation for |η|≤0.8 (normalized by max cl_a): b vs a =
  **11.2%**, c vs b = **0.02%**. The deficit is nearly uniform across the span
  (global circulation reduction, not a tip/truncation signature).
- Settledness: CL(step) rises Wagner-like and asymptotes cleanly; last-10-step
  drift 6.2e-5 (panel) / 7.4e-5 (particle) relative. Truncation ruled out
  analytically (missing wake beyond 40 chords is a ≪1% downwash effect) and
  empirically (b at 40 panel-chords ≡ c at 5 panel-chords + particles).
- Attribution follow-up (`task1_followup_attribution.jl`):
  1. Orientation confound ruled out: semiinfinite steady solve on the *rotated*
     body with horizontal U∞ gives CL = 0.29385887, identical to case a.
  2. **Das-magnitude sensitivity is the smoking gun**: panel-wake march, 40
     steps — das=0.025c → CL 0.24296 (−7.2%), das=0.05c → 0.26183 (baseline),
     das=0.10c → 0.27530 (+5.1%). The settled steady state depends strongly on
     the first-wake-row length and moves *toward* the semiinfinite answer as
     das grows. A converged steady limit must not depend on this discretization
     parameter at O(10%).
  3. Deficit is in the *solve*, not pressure integration: TE shed circulation
     Δμ (steady vs marched, 40 steps) ratio = 0.912 mid-span, 0.916 summed.
- Implicated code (for Task 4): first-row wake geometry/kinematics — `update_TE!`
  re-pins row-1 to TE+Das each step (`src/FLOWPanel_simulate.jl:891-921`); Das is
  set once and never refreshed; `velocity_te` contains only freestream+kinematic,
  induced TE velocity never added (`src/FLOWPanel_abstractliftingbody.jl:43-47`,
  `src/FLOWPanel_frames.jl:267`); and/or how the das-long row-1 buffer enters the
  Kutta/solve influence vs the semiinfinite sheet in the steady path.
- Note for Task 2/3b: any unsteady validation using these wake models inherits
  this ~11% quasi-steady offset; per plan, split static offset vs dynamic
  increment when comparing to experiment.

Scripts: `code_audit/scripts/task1_steady_wake_consistency.jl` (~3 min),
`code_audit/scripts/task1_followup_attribution.jl` (~4 min).
Outputs: `code_audit/results/task1/task1_spanwise.csv`,
`task1_cl_history.csv`, `task1_overlay.png`, `task1_cl_history.png`.

### Task 1b

IN PROGRESS (2026-07-16).
- Step 1 done: no existing side-by-side derivation of the two wake-coupling schemes found.
  `docs/dirichlet_potential_theory.md` proves the key ingredient (wake-velocity-through-
  sources reproduces the wake interior potential up to a constant, verified by
  `test/dirichlet_potential_test.jl`) but does not compare the schemes. Created
  `docs/wake_solve_schemes.md`: derives scheme A (semiinfinite attached panel in matrix)
  vs scheme B (das-length attached panel in matrix + free wake through σ only; confirmed
  every wake→body `influence!` branch in `src/FLOWPanel_simulate.jl:358-444` passes
  `scalar_potential=false`), states equivalence preconditions P1–P4, and ranks P3 (handoff
  filament cancellation at TE+Das between the attached panel's downstream edge, strength
  Δμ(t), and free row-1's upstream edge, strength Δμ(t−dt), opposite orientation — any
  sign/geometry mismatch leaves a spurious ~Δμ/(4π·das) normal velocity at TE control
  points) as the top suspect for the das-sensitive, span-uniform ~9% deficit. Document
  includes the manual linear-system cross-substitution experiment design (§5).
- Step 2 done: clean-context reviewer verified `docs/wake_solve_schemes.md`, corrected the
  Green-identity ingredient (eq. 2), and added precondition P5 (Kutta coupling reads the
  wake-trace-shifted doublets) as top suspect.

#### Numerics (steps 3–4 — MECHANISM CONFIRMED)

**Verdict: the deficit is caused by the velocity-through-sources wake coupling itself
(scheme B's formulation), not by any geometry/sign/handoff bug.** Exact decomposition of
the settled deficit at the fixed point (identity D1 = −C·G_A⁻¹(T1+T2), closure 6e-13):

- **95.0% of the deficit is T2 = φ_tr − φ_σ,fw** — the error of representing the
  free-wake *potential* through wake-velocity-sampled sources instead of including it in
  the Dirichlet RHS. Split (exact, closure 1e-15):
  - **P5 (trace shift read by Kutta coupling), 57.6%**: E2a = −C·G_A⁻¹(K_int·φ_tr).
    Even a perfect single-layer representation shifts the solved doublets by the wake
    trace (μ_B = μ_phys + φ_tr-ish), and the Kutta/attached-wake coupling reads that
    shift. NOTE: the naive local estimate −C·φ_tr is only 4.1% of the deficit — the
    contamination is dominated by the *distributed chordwise structure* of K_int·φ_tr
    amplified through G_A⁻¹, not by the TE-pair jump alone. Correlation with the measured
    deficit is high (cor(D1, −Cφ_tr) = 0.996, everything is span-uniform) but magnitude
    requires the full eq.-(3)/G⁻¹ route.
  - **P1/P2 (sampled-σ consistency), 37.4%**: E2b = C·G_A⁻¹·e4a, where
    e4a = φ_σ,fw − (I−K_int)φ_tr is the eq.-(2) residual. e4a is only 1.5% of φ_σ,fw in
    norm (‖e4a‖/‖φ_σ,fw‖ = 0.0155, distributed, not TE-concentrated) but is amplified
    ~25× by the Kutta-coupled solve.
- **5.0% is T1** (composite sheet vs semi-infinite sheet at matched strengths): of which
  4.8% is *physical* wake relaxation (the marched wake sags up to 0.36 m ≈ 1.2c below
  the rigid plane over 40 chords) and 0.3% is truncation+row-discretization+core_size.
  **P3 (handoff filament/sign/geometry bug) is REFUTED**: the off-body composite-sheet
  check (das panel + free wake vs semi-infinite, μ* strengths) agrees to 0.6–1.7%, and
  the §5.1 sanity holds to 2e-15.

Key measurements (full Task-1 config, marched state reproduced to CL = 0.26205059 vs
Task 1's 0.262051; frozen-geometry scheme-A solve reproduces sum-Δμ ratio 0.9168 and
mid-span 0.9131 vs Task 1's 0.916/0.912):
- Cross-substitution residual r_A = G_A μ* + φ_σ,A is *distributed* (rms non-TE 0.029 >
  TE rows 0.0053) — consistent with the trace-shift/gauge contamination of μ*, not a
  TE-localized defect.
- One-shot A′ solve (das geometry + free-wake potential in RHS, marched wake state):
  sum(Δμ_A′)/sum(Δμ_A) = 0.9415 vs marched 0.9168 — recovers ~30% one-shot; at the A′
  *fixed point* (wake re-fed with stronger Δμ each step) the remaining error is T1-type
  only (~5% of the current deficit ⇒ ~0.5% in CL).
- Das-sensitivity is explained: both E2a and E2b scale with the near-TE wake standoff
  (φ_tr structure over the TE panels), monotone in das, span-uniform, identical for
  panel/particle wakes — matching every Task 1 signature.

**Fix direction (for Task 4): include the free-wake scalar potential in the Dirichlet
RHS (scheme A′) instead of (or in addition to removing) the u_fw contribution to σ.**
Concretely: in `_steady_aerodynamics!` (`src/FLOWPanel_simulate.jl:358-444`) probe the
wake with `scalar_potential=true` at body control points, carry φ_fw into the solve RHS
(`solve!` Dirichlet wrapper currently zeroes `body.potential`,
`src/FLOWPanel_solver.jl:227-246`), and *remove* the wake-induced velocity from
`set_strengths!`'s σ (keep it for post-solve surface velocity/monitors). This
eliminates T2 identically (both the P5 trace shift and the P1 sampling error); panel
wakes support scalar potential (PanelWake rows are doublet panels); particle wakes need
the panel-row buffer's potential plus a treatment for the particle far wake (particles
carry no scalar potential — the buffer-dominated φ_tr observed here suggests the near
rows carry most of it, but this needs a Task-4 quantification; alternatively keep
velocity-through-σ for the particle portion only).

Scripts: `code_audit/scripts/task1b_numerics.jl` (~5 min local; TASK1B_SMOKE=1 for a
fast smoke). Outputs: `code_audit/results/task1b/{task1b_summary.txt, task1b_run.log,
task1b_spanwise.csv, task1b_residual_cells.csv, task1b_spanwise.png}`.

#### Fix-candidate quantification (2026-07-17)

Three-way one-shot solve on the identical frozen settled state (Task-1 case-b march
reproduced; scheme A = semiinfinite steady on frozen geometry): production B vs A′
(free-wake potential in Dirichlet RHS) vs Ryan's affine-Kutta alternative from
`docs/wake_solve_schemes.md` (ALT: G_B μ̃ = −S(σ₀+σ) + W c, γ = Cμ̃ − c), with the
Kutta-trace vector c evaluated three ways (exact C·φ_tr; trapezoid and Simpson line
integrals of the free-wake velocity along the lower→upper TE control-point chord).

**Verdict: A′ dominates. The affine-Kutta alternative removes the P5 trace shift
exactly (identity verified to 3e-12) but retains the e4a sampling residual, as
predicted — one-shot deficit recovery 18.0% (ALT) vs 29.7% (A′). Quadrature choice
for c is a non-issue.**

- Sum-γ ratios vs scheme A (mid-span in parens): B = 0.9168 (0.9131);
  A′ = 0.9415 (0.9384), recovery 29.7%; ALT c-exact = ALT c-trap = ALT c-Simpson =
  0.9318 (0.9284), recovery 18.0%. ALT with c=0 reproduces B to 0 (exact sanity).
- Analytic identity confirmed: γ(ALT,c-exact) − γ_A′ = −C·G_B⁻¹·e4a to max abs err
  1.9e-13 (rel 3.4e-12). The retained-e4a term is 11.7% of the measured deficit via
  the one-shot G_B route (the fixed-point G_A-route figure is E2b = 37.4%).
- Fixed-point projection (from the exact Task-1b decomposition): A′ retains only T1
  → ~95% recovery (~0.5% CL error); ALT retains T1 + E2b → ~58% recovery (~4–5% CL
  still missing). One-shot implied CLs: A′ 0.26911, ALT 0.26633 (target A: 0.29386).
- c quadrature: trapezoid rel err 1.0e-4, Simpson 2.9e-9 vs exact — all three give
  identical γ to 4 digits. Straight centroid-to-centroid TE chord verified interior
  at all stations (|d| mid-span 4.0e-4 m). Sign pairing empirically u_f = +∇φ
  (cor(c_trap, C·φ_tr) = +1.0000).
- Discrete Green identity in isolation (doc check 2, gauge-clean since
  (I−K_int)·1 = 0 to 2.3e-13): ‖e4a‖/‖φ_σ,fw‖ = 0.0155 — an irreducible P1/P2
  point-sampling residual that any velocity-through-σ representation of the wake
  (including ALT) inherits; only the explicit-potential RHS (A′) removes it.

Consequence for Task 4: **implement scheme A′** — wake scalar potential in the
Dirichlet RHS, wake velocity removed from σ (kept for post-solve monitors). The
affine-Kutta route is not a substitute at current panel resolution; its appeal
(velocity-only wake evaluation, works for particle far wakes with no scalar
potential) could still make it a fallback for the particle-far-wake portion, but it
caps recovery at ~58% of the deficit here. Open Task-4 item stands: how much of
φ_tr the panel-row buffer carries for `PanelParticleWake` (buffer-dominated per the
Task-1b probes, needs quantification).

Script: `code_audit/scripts/task1b_alternative_formulation.jl` (includes
task1b_numerics.jl and extends it; ~6 min full, TASK1B_SMOKE=1 for smoke).
Outputs: `code_audit/results/task1b/{task1b_altform_summary.txt,
task1b_altform_spanwise.csv}`.

### Task 2

**Verdict: PASS — the unsteady particle path is not the failure mechanism.** With a
genuinely particle-populated wake, panel vs particle agree to **+0.35% in cl amplitude
and −0.17° in phase** at every station (gates were ~10% / a few deg), consistent with
Task 1's steady c→b result (0.02%). Combined with Task 1b, the unsteady failure
suspect list stays on the solve-side wake coupling (scheme A′ fix) and the φ̇ ladder
(Tasks 3a/3b), not on panel→particle conversion or particle influence.

Setup: pitching_wing.jl wing (n_span=15, n_airfoil=31, n_endcap=5, 1124 cells), the
validation case's own forcing α = 3.94° + 1.99°·sin(2π·4.01t) (k = ωc/2U = 0.0382),
3 cycles, dt = 0.5 chords/step (164.7 steps/cycle, 495 steps), wake trimmed at
2 spans (downstream boundary 8.75c aft of pivot), unsteady PressureBernoulli,
Backslash + FastMultipoleBackend(8, 0.4, 40). Per-station cycle-3 least-squares
sinusoid fits (fit residuals ~1e-3 relative; cycle-2→3 amplitude drift ≤0.004% ⇒
washed out).

- **Gotcha found (affects any particle-wake pitching-wing run using example
  defaults):** `prepare_pitching_wing` defaults
  `panel_wake_rows = wake_length/das = 160` rows ≈ 80 chords of panel buffer —
  the entire trimmed wake stays panels, and every converted particle is born ~80c
  downstream, far outside the 8.75c `FrameBox` cull, and is deleted on the next
  step's maintenance (conversion in `shed_wake!` runs *after* the cull inside
  `propagate!`, `src/FLOWPanel_wake.jl:1157-1184`; hence the residual np=62 = one
  fresh row). The first Task-2 run therefore compared a panel wake with itself
  (agreement 0.009% — a good buffer-invariance sanity but not the particle test).
- Corrected run: `panel_wake_rows=10` (5-chord buffer, matching Task 1c), so the
  5→8.75c band is real particles influencing the body (np settles ≈372).
  Cycle-3 fits vs the stored panel run: dA = +0.347…+0.359% (span-uniform),
  dphase = −0.170…−0.175°, mean-cl shift +0.35%. Total CL amplitude:
  panel 0.114607, particle 0.115014.
- Interpretation note: this delta bundles three tiny effects — particle influence
  itself, the missing particle scalar potential in φ̇ (`allow_partial=true`;
  particles carry no φ), and the buffer-length difference (160 vs 10 rows). At
  k=0.038 their sum is 0.36%/0.17°; no need to decompose further. The Task-4 open
  item (how much φ_tr the buffer carries for PanelParticleWake) stands, but Task 2
  shows the practical error at 5-chord buffers is sub-half-percent at this k.
- Both models inherit the Task-1 quasi-steady deficit identically (mean CL ≈ 0.233
  at ᾱ=3.94° vs the semiinfinite steady level); per plan this is the static offset,
  already attributed in Task 1b — not re-counted here.

Scripts: `code_audit/scripts/task2_unsteady_panel_vs_particle.jl` (base pair,
~18 min) and `task2b_particle_short_buffer.jl` (short-buffer particle rerun,
~9 min; reuses stored panel histories). Outputs under `code_audit/results/task2/`:
`task2_station_fits.csv`, `task2_histories.csv`, `task2_cl_history.png`,
`task2_amp_phase.png`, `task2b_station_fits.csv`, `task2b_histories.csv`,
`run.log`, `task2b_run.log`.

### Task 3a

**Verdict: PASS — the φ̇ machinery (φ probe, exterior-trace handling, BE/BDF2
ladder, ALE w·∇φ term, inertial-trace kinetic term) is clean on the analytic
sphere added-mass problem. No code changed.** Combined with Task 0 (trace sign)
and Tasks 1/1b/2, the unsteady-failure suspect list now excludes the entire
unsteady-Bernoulli pressure path for wakeless bodies; remaining levers are the
scheme-A′ wake-coupling fix (Task 1b) and shedding/wake history (Task 3b).

Setup: closed UV-sphere `NonLiftingBody{ConstantSource}` (R=1, ρ=1.225),
prescribed smooth-ramp translation U(t) = ½Umax(1−cos(πt/T)) x̂ (Umax=1, T=1 s)
via frames/maneuver, uinf=0, no wake, unsteady `PressureBernoulli` +
`ForceMonitor(NoNormalization)`, Backslash + FastMultipoleBackend(8, 0.4, 40).
Exact: p−p∞ = (ρR/2)U̇cosθ + (ρU²/8)(9cos²θ−5); F = −½ρVol·U̇.

- **Per-panel pressure** (peak-U̇ step, non-pole panels, normalized by
  max|P_exact|): 960 panels/dt=0.02 → median 1.6%, p90 3.7%, max 4.5%;
  2208 panels/dt=0.01 → median 1.0%, p90 2.4% (converging). Mesh-pole cap
  panels are *better* (max 0.3%), no polar pathology.
- **Force** (steps ≥ 3, ramp, |U̇| ≥ 20% peak): |Fx−F_exact|/|F_exact| =
  4.9% (528 panels) → 3.7% (960) → 2.7% (2208): pure spatial-discretization
  constant, converging to 0. Signature that the time ladder is clean: the
  ratio Fx/F_exact is flat in time (spread ~0.4% across the whole ramp,
  BE-startup step included after step 3). Note the discrete added mass
  *overshoots* continuum (ratio 1.037 at 960 panels) even though
  V_mesh/V_exact = 0.984 — a volume correction goes the wrong way; don't
  "fix" the residual that way.
- Lateral forces O(1e-13) of peak; post-ramp (U̇=0) force is exactly 0 to
  double precision. First-step (BE) and BDF2 steps both inside gates.
- Permanent CI gate: `test/runtests_unit_added_mass.jl` (12 tests, ~15 s;
  528 panels, dt=0.04; gates: p90 pressure < 7%, force < 7%, Fx/F_exact
  time-spread < 2%, lateral < 1e-8·F_peak, post-ramp < 2%·F_peak),
  registered in `test/runtests.jl` and in the `agent_policies/TESTING.md`
  matrix.

Script: `code_audit/scripts/task3a_sphere_added_mass.jl` (~1 min at defaults;
env `TASK3A_NTHETA/TASK3A_NPHI/TASK3A_DT` for refinement). Outputs:
`code_audit/results/task3a/{task3a_run.log, task3a_force_history.csv,
task3a_panel_pressure_peak.csv}` (960-panel default config).

### Task 3b

IN PROGRESS (2026-07-17).

### Task 4 — Conclusions & next steps

(pending)
