# Unsteady-case code audit: plan

## Objective (read this first)

**The entire objective of this effort is: (1) identify why FLOWPanel.jl is failing on
unsteady validation cases, and (2) fix it.** Everything below is in service of those two
goals. Steady validations pass; unsteady ones fail with wrong force magnitudes, and an
unsteady run marched to steady state disagrees with the steady solver.

## How to work this plan (protocol for every agent)

- Progress is tracked in `code_audit/log.md` (checkboxed tasks + per-task findings).
  **Orchestration**: a master agent runs the plan by launching **one fresh-context
  subagent per task**, sequentially (each task's subagent may depend on prior findings
  in log.md). The master agent does not execute tasks itself and does not need to clear
  its own context between tasks.
- **Per-task subagent protocol**: read `code_audit/log.md`, pick the first unchecked
  task, mark it in progress, and do exactly that one task. When complete: check it off
  in log.md, write a concise findings entry under that task (verdict, key numbers,
  file:line of anything implicated, paths to plots/CSVs/scripts produced), then return
  a summary to the master agent. Do not start the next task.
- Keep runs small (coarse meshes, tens of steps): these are correctness gates, not
  convergence studies. Prototype scripts in `code_audit/scripts/` (create it); final CI
  gates go in `test/`.
- **Run authorization & HPC**: the user has granted permission to launch scripts
  autonomously for this audit. If a simulation would take longer than ~20 minutes on this
  Mac, run it on HPC instead, using the repository at `/home/rander39/projects/FLOWPanel.jl`
  (its Manifest already points to the needed dependencies). Read `agent_policies/HPC.md`
  before any HPC/Slurm work.
- Required reads before touching code: `CLAUDE.md` (routes to `agent_policies/WORKFLOW.md`
  for src/ edits, `agent_policies/TESTING.md` for test selection).

## Background (self-contained; do not re-explore)

FLOWPanel.jl, branch `fastmultipole`. 3D panel method; unsteady time-marching via
`simulate!` with panel or panel+particle wakes.

**What fails / what works (from the user):**
- Rotor in hover is effectively steady in the body frame and is run with **steady**
  Bernoulli — do NOT chase unsteady/ALE/vector-potential theories on the rotor.
- The live failure: **unsteady Bernoulli on the pitching wing with a panel wake** (panel
  wake has a genuine scalar potential, so particle-solenoidal-contamination explanations
  do not apply).
- Static polars match experiment to α≈6°, so the steady solve + steady pressure are sound
  at low α. Suspects narrow to: the φ̇ machinery (exterior-trace sign, BE/BDF2 ladder,
  ALE term), wake-model consistency, and unsteady shedding.
- Ruled out by the user / earlier audits — do not re-investigate: steady-limit ALE algebra
  on the rotor; particle-wake ∂A/∂t on the rotor; first-step phi_dot=0; rotor reliability
  replay; backend & Kutta condition (BRAINSTORM/001); `semiinfinite_wake` misconfiguration
  (both failing examples set it correctly).

**Code map (verified file:line, 2026-07-16):**
- Unsteady Bernoulli: `src/FLOWPanel_simulate_monitors.jl` — monitor entry :206-237
  (steady mode → body-relative tangential trace via
  `_pressure_fill_relative_surface_velocity!` :232; unsteady mode → inertial trace
  `_pressure_fill_inertial_surface_velocity!` :231 plus phi_dot :222).
  `_pressure_bernoulli_phi_dot!` :349-425; `_variable_bdf2` :427-433. Final assembly
  `calcfield_P!` `src/FLOWPanel_postprocess.jl:944-993`: `P = ½ρ(U∞²−|U|²) − ρφ̇`.
- **Sign-check target**: monitors :384-391. φ = probed scalar potential + `uinf·x`; then
  `has_grad_mu(body) && (phi -= body.strength[p, get_Gammai(body)])` converts the
  canonical *interior* limit (self panel contributes +μ/2, `_self_limit` at
  `src/FLOWPanel_elements_fmm.jl:165-177`) to the *exterior* trace. FLOWPanel's
  stored-strength convention is opposite the usual exterior-minus-interior convention
  (see comment at monitors :2057-2058). A wrong sign here corrupts every φ̇.
- Frames: `kinematic_velocity!` `src/FLOWPanel_frames.jl:311-357`. After it runs:
  `body.velocity` = body-relative fluid velocity (induced − kinematic);
  `body.velocity_kinematic` = inertial rigid-body velocity w; inertial fluid velocity =
  sum of the two. The ALE conversion in phi_dot is `∂φ/∂t = Dφ/Dt − w·∇φ` with
  `∇φ = body.velocity + w − excluded_gradient` (monitors :402-410; `excluded_sources` are
  vector-potential-only wake sources, :376-382).
- Time loop: `src/FLOWPanel_simulate.jl:734-881` — per step: maneuver! → solve
  (`_steady_aerodynamics!` :755) → optional Γ under-relaxation :770-779 → monitors :795 →
  save → `propagate!` wake :853 → `propagate_kinematics!` :863 → `shed_wake!` :873.
  Shed row-1 strength = post-solve Δμ = μ_upper − μ_lower (`_get_wakestrength_mu`,
  `src/FLOWPanel_liftingbody.jl:1565-1573`; note deliberate sign flip for swapped node
  order, simulate.jl:947-954). Row-1 upstream nodes re-pinned to TE+Das every step
  (`update_TE!` simulate.jl:891-921). `Das` is set once by `initialize_Das!`
  (simulate.jl:281-315) or `set_wake_Das!`, never updated in the loop; `velocity_te`
  contains only freestream+kinematic (`src/FLOWPanel_abstractliftingbody.jl:43-47`,
  `src/FLOWPanel_frames.jl:267`) — induced TE velocity is never added (watch-list item).
- Wake models: `PanelWake` `src/FLOWPanel_wake.jl:76-157`; `PanelParticleWake` :895-998 —
  oldest panel row converts to particles in `_convert_to_particles!` :1230-1269 (trailing
  filaments from spanwise ΔΓ; shed particles from temporal ΔΓ, `Γ − Γ_tm1`, :1253-1256).
  Particle convection: explicit Euler `FLOWVPM._euler` (wake.jl:1157-1184).
- Wing under test: reuse the wing of `examples/pitching_wing.jl` — steady body built at
  :840 (`semiinfinite_wake=true`, `set_wake_Das!(body, Uinf)` :843); unsteady body at
  :949 (`semiinfinite_wake=false`), `set_wake_Das!(wing, Uinf(t0); magnitude=das)` :961.
  Wake model switch: `wake_model=:panel` or `:particle`.
- Spanwise loading: `calcfield_sectionalforce!` `src/FLOWPanel_postprocess.jl` (~:1146);
  `ForceMonitor` requires a pressure monitor earlier in the monitor tuple
  (`monitor_requires`/`audit_monitors`, monitors :16-53).
- Existing data possibly reusable for Task 2: `examples/pitching_wing_convergence.jl`
  outputs in `data/pitching_wing_convergence/` (panel wake by default; particle wake via
  env `PITCHCONV_WAKE_MODEL`).
- Forces: `ForceMonitor` integrates `:P` only — all unsteadiness must arrive through φ̇.
  `KuttaJoukowskiForce`/`SurfaceVorticityForce` are quasi-steady (no φ̇) and cannot match
  unsteady loads. `PressureLaplace` is unvalidated for unsteady; only its
  `:edge_difference` gradient mode is usable (`corrected_hessian`/`surface_velocity` are
  known-pathological, see `data/pitching_wing_pressure_comparison/*_summary.csv`).

## Tasks

### Task 0 — Exterior-trace sign check (cheap, decisive)
Verify that `phi -= μ` at monitors:391 yields the *exterior* φ trace under FLOWPanel's
sign conventions.
1. Analytic: trace the doublet scalar-potential kernel sign
   (`src/FLOWPanel_elements_fmm.jl`) against the K&P-vs-GeometricTools normal convention
   (CLAUDE.md notes element kernels carry explicit sign flips), and confirm interior limit
   +μ/2 ⇒ exterior = interior − μ (not +μ).
2. Numerical (no time-marching): steady-solve any lifting body (the pitching_wing.jl wing
   at fixed α), probe φ at points offset slightly outward along +n̂ (offset ~1e-3 × panel
   size) using the standard off-body probe, and compare to the on-surface φ from the
   :384-391 recipe. Agreement to O(offset) ⇒ sign correct; matching the interior value or
   off by a full μ ⇒ sign bug found (then fix it and re-run the check).
Deliverable: verdict paragraph + comparison numbers in log.md; script in
`code_audit/scripts/`.

### Task 1 — Steady wake-model consistency (spanwise lift)
Same wing as `examples/pitching_wing.jl`, fixed α in the validated range (4–5°), steady
conditions, **steady** PressureBernoulli everywhere so only the wake model varies:
  a. `semiinfinite_wake=true`, steady solve — trusted baseline;
  b. unsteady march, `wake_model=:panel`, `semiinfinite_wake=false`, constant Uinf/α, run
     until forces settle and the wake spans many chords;
  c. same as (b) with `wake_model=:particle`.
Compare spanwise sectional lift + total CL. b→a isolates finite-wake truncation +
shedding/Das; c→b isolates panel→particle conversion. Acceptance: settled b and c within
~1% of a in sectional lift away from tips. Deliverable: overlay plot + CSV + script;
deviations recorded in log.md.

### Task 1b — Diagnose the semiinfinite vs finite-wake solve discrepancy (USER-STAGED, PRIORITY)
Staged by the user 2026-07-16 after Task 1's finding (marched wake settles ~11% low; deficit
is in the solved circulation, Δμ ratio 0.912; strong Das-length sensitivity). **All tasks
beyond those already complete are paused until Task 1b is done** (Task 2's run may finish
in the background; its analysis waits).

User's understanding to verify: the `semiinfinite_wake=true` and `=false` solvers produce
different source and doublet strengths because in the former the wake-induced *potential*
is included in the doublet influence matrix, while in the latter the wake-induced potential
is not explicitly included — instead the wake-induced *velocity* is used to affect the
source strengths.

Steps:
1. Search `docs/` and `theory/` for an existing document deriving the two formulations
   (semiinfinite/attached-wake-in-matrix vs finite-wake/velocity-through-sources). If it
   exists, audit it for correctness. If not, create it (derive both schemes from the
   Dirichlet BIE as implemented: `solve!` paths in `src/FLOWPanel_solver.jl`,
   `_induced_wake` in `src/FLOWPanel_elements_fmm.jl:1004-1074` (direct) and :1076+
   (buffer), wake-on-body RHS assembly in `src/FLOWPanel_simulate.jl:358-444`).
2. Have the derivation verified for correctness by a NEW agent with clean context before
   relying on it.
3. With the verified derivation in hand, hypothesize why the solved circulation differs
   ~9% between the schemes; hunt the code for culprits, specifically including:
   - a sign error in how the transition (attached/das-length row-1) panel enters the
     influence matrix;
   - the shed wake row carrying an upside-down orientation (normal/winding flip);
   - anything else the derivation exposes (e.g. missing wake-potential term in the
     finite-wake RHS, source-strength contamination).
4. Decisive numerical check (user-suggested): take the settled unsteady solution and
   manually build the full linear systems for BOTH schemes (semiinfinite=true and =false)
   and verify each solution satisfies its own equations and quantify which terms differ —
   i.e., substitute the marched solution into the semiinfinite system's residual and vice
   versa, and compare row-by-row (wake columns, RHS terms) for a mid-span control point.
Deliverables: the (audited or new) theory document; a log.md entry naming the mechanism
with file:line; the manual-linear-system script in code_audit/scripts/.

### Task 2 — Unsteady panel vs particle wake (spanwise lift)
Same wing, gentle sinusoidal pitching (amplitude ~1–2°, reduced frequency within the
pitching_wing_convergence range), identical dt/discretization, unsteady PressureBernoulli.
**First** check whether `data/pitching_wing_convergence/` already holds runs for both wake
models (env `PITCHCONV_WAKE_MODEL`) — compare stored outputs before launching anything.
Compare per-station lift amplitude and phase, panel vs particle. Disagreement with Task 1
passing isolates the *unsteady* particle path (temporal-ΔΓ shed particles, particle
influence during the cycle). Deliverable: amplitude/phase table + plots; log.md entry.

### Task 3a — Sphere added-mass gate (analytic, no wake)
`NonLiftingBody` sphere, prescribed translation U(t) (smooth ramp) via frames/maneuver,
unsteady PressureBernoulli + ForceMonitor. No shedding ⇒ isolates φ probe, trace sign,
BE/BDF2 ladder, ALE term, inertial-trace kinetic term. Exact solution:
φ_ext = −½U(t)R³cosθ/r², surface pressure closed-form, total force F = −½ρVol·dU/dt.
Gate: per-panel P within a few % away from the poles; force error < ~2–5% once BDF2 is
active (step ≥ 2). Build as a script first, then distill into a `test/` gate
(pattern-match `test/runtests_unit_*.jl`, register in `test/runtests.jl`).

### Task 3b — Wagner/Theodorsen gates (panel AND particle wake)
High-AR flat/symmetric wing (reuse Step-1 wing or the pitching_wing.jl generator). Two
sub-gates, each run with `wake_model=:panel` and `wake_model=:particle`:
- Step change in α → CL(s)/CL_steady vs Wagner (Jones: ψ(s) ≈ 1 − 0.165e^{−0.0455s} −
  0.335e^{−0.3s}, s in semichords). Gate: within ~5–10% for s ≳ 1; CL(0⁺) ≈ 0.5·CL_steady.
- Gentle sinusoidal pitch → amplitude ratio + phase vs Theodorsen C(k). Gate: amplitude
  within ~10%, phase within a few degrees.
Normalize by the steady solver's own CL_steady (not 2π) so the gate tests the unsteady
*increment*, not the lifting-line limit. Make both permanent CI gates in `test/`.
Failure-space split: 3a fails ⇒ φ̇ machinery/sign; 3a passes + 3b(panel) fails ⇒
shedding/wake history; 3b(panel) passes + 3b(particle) fails ⇒ particle conversion.

### Task 4 — Re-evaluate, conclude, and plan the fix
Review all log.md entries from Tasks 0–3b. Deliverables:
1. A conclusion: which rung of the ladder breaks first, and the implicated code path
   (file:line).
2. A recommended fix (or, if already fixed during a task, confirmation that the gates now
   pass and the pitching-wing validation improves).
3. **Append new checkboxed tasks to log.md** for implementing/verifying the fix — the
   objective is not met until the unsteady failure is understood AND fixed, with the
   Task-3 CI gates passing and the pitching-wing case re-validated.
Approved follow-on diagnostics if the ladder is ambiguous: term-by-term pressure
decomposition fields (kinetic, −ρDφ/Dt, +ρw·∇φ) on a failing case; panel-wake-only rotor
sanity run; pitching-wing-vs-experiment split into static offset + dynamic increment;
`Das` magnitude sensitivity sweep (×0.5/×2; if sensitive, per-step Das refresh from
current TE velocity).
