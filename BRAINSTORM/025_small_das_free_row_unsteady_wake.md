# Small-Das + Free-Row Near Wake: Restoring Fully Unsteady Prediction

**Opened:** 2026-08-20 (Ryan-directed, from the 018 near-wake modeling-scope
discussion; ledger §2026-08-20 "RYAN NOTE (modeling-scope limitation)" +
its 0.25c-provenance amendment are the source record).

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

**Status: STAGED ONLY — do NOT implement.** Ryan 2026-08-20: "stage item
..., but don't implement now." Everything below exists so a fresh agent can
execute later without exploration. Related but distinct: BRAINSTORM/024
(N=0 convert-at-shed) goes the *opposite* direction (fewer free rows,
steady-model cleanup); both can coexist because they serve different
regimes.

## The problem this item solves (self-contained)

Current production architecture (018 campaign): a **prescribed rigid Das
row** (TE → TE+Das; the implicit Kutta carrier inside the linear solve) is
sized to enforce the particle-σ separation criterion — Das = λσ ≈ 1 local
chord under the co-scaling law (σ = 0.313c, λ = 3.4, C1 clearance band
d/σ ≈ 2.6–3.4). One free wake-panel row (nwakerows = 1) rides beyond it,
then particles.

**Ryan's recorded insight (2026-08-20): this is STEADY-STATE modeling
only.** Prescribing the first ~chord of wake replaces the region that
governs unsteady lift (Wagner/Theodorsen deficiency — freshly shed
vorticity convecting away from the TE) with a quasi-steady element; F1b
additionally bakes in a frozen settled-wake downwash table. No sub-transit
timescale response is possible: BVI, gusts, pitching, flutter are out of
scope by construction.

**The alternative (this item):** keep Das SMALL — just large enough to
enforce the Kutta condition in a converged sense — and let **N free
wake-panel rows** carry the near wake time-accurately (convecting with
kinematics + induced velocity) until the free sheet spans the 3–4σ
particle clearance; convert to particles only beyond it. That is the
classical time-accurate scheme with a particle far wake.

**Why it isn't done today:** past simulations blew up due to the free wake
panels (Ryan, 2026-08-20 — the instability is the blocker, not the
concept). Campaign-adjacent history: the 018 Phase-5 legacy-σ nwakerows
ladder was non-monotone and rejected (see
`018_dji9443_hover_convergence_campaign/phase_05_shedding_distance.md` and
the N=8 spot check `p018_nrows8`); the executing agent should read that
phase file for the recorded failure modes before designing screens.
Mitigations that postdate the blow-ups and change the odds: the 019
stability screen/tripwire protocol (min_sigma, max_u, |Γ|/σ², σ≤0;
absolutes under station-σ), viscous CoreSpreading particle handoff, the
merge policy, and the panel-wake kernel `core_size` regularization.

## Critical evidence status (do not re-derive; do not over-trust)

- **The ~0.25c "attachment band edge" is NOT a validated Kutta floor**
  (ledger amendment 2026-08-20). It is the empirical lower edge of 014's
  wing log-plateau (dt-converged CL grows +0.205%/doubling of Das), and its
  position was substantially set by particle-σ clearance at the measured
  σ/c (0.25c ≈ 3σ on that wing). 014's sheet/particle split sweep
  (`examples/ssw_sheet_particle_split.jl`, data
  `data/ssw_sheet_particle_split/`) shows CL moves <0.1% for sheet-buffer
  0.25c→8c once particles are kept away.
- **Therefore the Kutta-only lower bound for the rigid row is UNTESTED and
  plausibly well below 0.25c** — likely set by TE panel size / wake-core
  scale. **Phase 1 of this item MUST measure it**: a rigid-row-length
  mini-sweep with the free-row sheet holding clearance fixed (the first
  clean isolation of the Kutta requirement from clearance).
- σ-side stability boundaries (018/019): uniform σ < 0.030–0.035R ignited
  every viscous run; 0.04R survives; station-σ law screens live in
  phase_16. Free-ROW instability is a separate mode — screens must watch
  both.

## Code map (verified 2026-08-20 on branch `fastmultipole`)

All machinery exists; this item is configuration + stabilization, not new
architecture. Anchors (same map as BRAINSTORM/024 §Code map, which see for
line-level detail):

- `src/FLOWPanel_wake.jl` — `PanelWake` (struct ~:98; constructor sizes
  storage from `nwakerows`, default 100 standalone / 3 via
  `PanelParticleWake`): free rows already convect with local flow;
  knobs `shed_with_induced_velocity` (first-row convection),
  `unsteady_filament`, `include_final_filament`, `core_size` (panel-kernel
  regularization, default 1e-3), `freestream_convection` (=true forces
  straight-sheet; diagnostic only). `PanelParticleWake` :1814+ composes it
  with the particle field; conversion (`_convert_to_particles!` :2137+)
  converts the OLDEST row each step once the buffer is full, i.e. handoff
  distance ≈ Das + N·(local travel per step). `StationSigmaOverlap`
  (:731–778) is node-column-indexed and geometry-agnostic — works
  unchanged for any N.
- `src/FLOWPanel_simulate.jl` — `update_TE!` :1380 (pins newest row to
  TE+Das), `shed_wake!` :1416 (PanelWake shift+shed) / :1455
  (PanelParticleWake convert-then-shift).
- Driver `examples/rotor_hover_pressure_comparison.jl` — `NWAKEROWS` env
  (:163), Das knobs (`DAS_SIGMA_LAMBDA`/`SIGMA_CHORD_FRACTION` co-scaling;
  `DAS_UNIFORM_DSIGMA`; `DAS_ETA_KINEMATIC`). The launcher's uniform-D arm
  already encodes N-aware clearance bookkeeping
  (`Das_j = D·σ − (N−1)·travel_j`) — prior art for "total buffer = target,
  split between rigid row and free rows".
- Wing-scale fixtures for cheap screens: `examples/ssw_*` family (AR=6
  swept wing, 014's η/dt matrix and split sweep) — reuse these before any
  rotor run.
- **Kutta-cache constraint:** Das stays FROZEN per run (Route A attachment
  operator is assembled once — `src/FLOWPanel_kutta.jl` ~:528–535; per-step
  Das changes are rejected). Small-fixed Das is fine; do NOT couple Das to
  dt or per-step state (DAS_REFRESH is banned inside dt studies — 014).

## Staged plan (phases; each gated on Ryan)

1. **Kutta-floor mini-sweep (wing, cheap, local-scale):** fixed free-row
   sheet holding handoff at 3–4σ; sweep rigid Das ∈ {0.02, 0.05, 0.1,
   0.25, 0.5}c at fixed dt and mesh; measure CL/Γ̄ convergence and solve
   conditioning. Output: the actual Kutta-converged floor. (No src/
   changes expected — knobs exist.)
2. **Free-row stability screens (the real risk):** N ladder (e.g. 2, 4, 8)
   at small Das on the wing, then the coarse 40_40 rotor, 019 tripwire
   protocol, ≤4 threads locally or standard HPC screens. Levers if
   unstable: panel `core_size`, N, conversion cadence. Record failure
   modes against the phase_05 history.
3. **Rotor steady A/B:** small-Das+N-rows vs the production prescribed-Das
   carrier at matched settings — steady CT/Γ̄ should agree (it's the same
   converged state); disagreement localizes model-form error.
4. **Unsteady demonstration:** a case with known unsteady content
   (pitching wing from item 013's validation set, or prescribed RPM
   transient on the rotor) where the prescribed-Das model must fail and
   this one should not.

## Standing constraints

- STAGED ONLY until Ryan lifts the hold. When executing: 018-campaign
  rules apply (≤4 local threads; HPC per `agent_policies/HPC.md`; cluster
  Manifest stays Julia-1.11.7; deploy-skew quarantine on src/ —
  cluster-base + intended hunks, md5 first; judge from monitors CSVs,
  exit 0 ≠ health; sigma-ratio monitor columns are NaN under station-σ —
  absolutes only; notebook writes need Ryan).
- Dt-dependence is PROPER here (time-accurate free rows) — do not apply
  018's dt-independence requirement to the free sheet; DO keep the rigid
  Das dt-independent (Kutta cache + 014's oracle finding).
- Coordinate with BRAINSTORM/024: if 024's N=0 mode has landed, its
  convert-at-shed path and tests touch the same functions — rebase this
  item's reading on the then-current `shed_wake!`.

## Exit criteria

- [ ] Kutta floor measured and recorded (Phase 1)
- [ ] Free-row stability envelope mapped; blow-up mechanism named (Phase 2)
- [ ] Rotor steady equivalence A/B (Phase 3)
- [ ] One unsteady validation case demonstrating capability beyond the
      prescribed-Das model (Phase 4)
- [ ] Ryan review at each phase gate

## Log

- 2026-08-20 — Item staged on Ryan's order (stage, don't implement). All
  code anchors verified against the working tree same day; evidence
  status (0.25c provenance correction) recorded from ledger amendment.
