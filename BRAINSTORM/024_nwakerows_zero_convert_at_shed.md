# N=0 Convert-at-Shed: Eliminate the Free Wake-Panel Row

**Opened:** 2026-08-20 (Ryan-directed, from the 018 near-wake anatomy
discussion; see 018 ledger §2026-08-20 "RYAN NOTE (parked...)" — this item
promotes that parked cleanup to implementation).

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

**Status:** implementation delegated 2026-08-20. Not urgent; expected to be
result-neutral. Do NOT deploy to the cluster or wire into any live 018 run
without Ryan.

## What and why (self-contained — no exploration needed)

The unsteady near-wake behind a `RigidWakeBody` in a `PanelParticleWake`
simulation has three elements, in order downstream:

1. **Rigid Das row** — a prescribed doublet row TE → TE+Das (Das is a per-TE-
   node-column 3-vector). It is the *implicit Kutta carrier*: its influence is
   inside the linear solve (semi-infinite horseshoe in dense assembly, FMM
   buffer contract of 8 extra rows/panel, both Kutta routes). **This item does
   NOT touch it.**
2. **Free wake-panel row(s)** — `nwakerows` rows of free doublet panels that
   convect with the local flow (kinematics + induced velocity). The newest row
   is pinned to the rigid row's end by `update_TE!`. Production 018 runs use
   `nwakerows = 1`.
3. **Particles** — each step, once the row buffer is full, the *oldest* free
   row is converted to particles (`_convert_to_particles!`).

Consequence of N=1: first particles materialize at ≈ Das + one step of TE
travel from the body. Under the 018 co-scaling law (σ = 0.313·c, Das = λσ,
λ=3.4 production) that is ~6.6σ at NT=36, 0.75R — beyond the ~3–4σ clearance
that Das alone already provides (Phase-12 C1 band d/σ ≈ 2.6–3.4).

**Ryan's direction (2026-08-20):** since Das alone satisfies particle-σ
clearance, the free row is redundant buffering; he can see it causing problems
later and wants an **N=0 mode: convert the just-shed row to particles in the
same step it is created** (particles appear at the Das line; no free sheet
survives into the next solve). Expected result-neutral; the point is
simplification/optionality, not a fix. The 018 NT ladder
(13245449–454) doubles as the sensitivity read: the free-row extent shrinks
∝Δt (≈3.2σ → 0.8σ at 0.75R over NT 36→144), so NT-flat CT̄ under F1b is
evidence the row's extent is benign.

## Code map (verified 2026-08-20 on branch `fastmultipole`; all anchors checked)

- `src/FLOWPanel_wake.jl`
  - `PanelWake` struct ~:98–190; constructor sizes storage from nwakerows:
    `nodes[i_surf] :: (3, nwakerows+1, ncols+1)` (node rows),
    `strength[i_surf] :: (dim, nwakerows+1, ncols)` (:161, :165). Fields
    `nwakes::Array{Int,0}` (count of live rows) and `overflowed::Array{Bool,0}`
    (:103, :174–175, set true when the buffer first fills).
  - Docstring :76–96: `shed_with_induced_velocity` controls first-row
    convection; `unsteady_filament` / `include_final_filament` control the
    final-edge filament that represents/cancels the shifted-out row —
    **filament bookkeeping is the subtle part of this change** (see
    `_final_filament_strength`, and `write_vtk(::FilamentWrapper)` :603+ which
    early-returns while `!overflowed[]`).
  - `PanelParticleWake` struct :1814, constructor :1831+ (`nwakerows=3`
    default, passes through to `PanelWake` :1861).
  - `_convert_to_particles!` dispatch :2137–2150; legacy edge-jump method
    :2152; smooth `SurfaceVorticityConversion` :2304/:2328. Conversion is
    geometry-agnostic — reads actual node coordinates only; a comment at
    :2419 notes the smooth stencil is what already makes `nwakerows == 1`
    work (single-row wake completes its streamwise stencil from the body).
    `StationSigmaOverlap` is indexed by node column only (:731–778).
- `src/FLOWPanel_simulate.jl`
  - Time loop calls: `update_TE!` at :678 (pins newest-row nodes to TE+Das;
    PanelWake method :1380, PanelParticleWake forwards :1412) and
    `shed_wake!` at :1357.
  - `shed_wake!(::PanelWake, system)` :1416–1451: shifts rows back one slot,
    writes the new first-row strengths from the just-solved μ jump
    (`_get_wakestrength_mu`), sets `overflowed[]` when about to overflow,
    increments `nwakes[]` (capped at storage−1).
  - `shed_wake!(::PanelParticleWake, system)` :1455–1469: `buffer_full =
    nwakes[] >= n_rows-1`; **if full, converts the oldest row FIRST, then
    shifts + sheds** — i.e. conversion happens at row age = nwakerows steps.
- Driver: `examples/rotor_hover_pressure_comparison.jl` reads
  `NWAKEROWS` env (:163) → `PanelParticleWake(..., nwakerows=...)`; banner
  prints `nwakerows:`; metadata records it (grep `nwakerows` in the driver
  and `src/FLOWPanel_metadata.jl` for the exact keys).

## Design (Route decision left to implementer, but constraints are fixed)

Two candidate routes; pick whichever satisfies ALL invariants below with the
least surface area, and record the choice + rationale in this item's Log:

- **Route I (preferred if clean): true `nwakerows = 0`.** Storage collapses to
  one node row (the TE+Das line) and zero panel rows. In
  `shed_wake!(::PanelParticleWake)`: compute the would-be first row's
  strengths from the solved μ, place particles directly (reusing the
  conversion machinery on a virtual/staging row), never grow `nwakes`.
  Watch: every consumer that assumes ≥1 panel row (VTK writers, filament
  writer's `i_row = nwakes[]`, influence evaluation of `PanelWake`, warmstart
  save/load of wake state, replay manifest round trip).
- **Route II (fallback): age-0 conversion with N=1 storage.** Keep
  `nwakerows=1` arrays; after `shed_wake!(pw, ...)` writes the fresh row,
  immediately convert it and reset `nwakes[] = 0` so the next solve sees no
  free sheet. Fewer storage-shape changes, but the filament/overflow
  bookkeeping must be re-derived for the empty-buffer steady state, and the
  "convert BEFORE shift" order in the current code inverts.

**Invariants (hard requirements):**
1. **Bit-identity for nwakerows ≥ 1** — zero behavior change for every
   existing configuration. The 016 golden-reference testset in the wake suite
   is the gate (it pinned bit-identity for the station-method dispatch the
   same way; extend, do not weaken).
2. **Circulation conservation at N=0**: total shed circulation per step
   (trailing + unsteady families) matches the N=1 case's conversion output on
   the same μ history — same particle count per shed or a disclosed
   difference, no lost filament. The final-edge filament semantics
   (`unsteady_filament`, `include_final_filament`) must remain consistent:
   with no surviving sheet, whatever the final filament used to
   represent/cancel must be carried by particles or provably zero.
3. **Kutta carrier untouched**: no change to Das, `update_TE!` pinning of the
   rigid row, or any solve-side influence of the rigid row.
4. **Both conversion strategies** (`LegacyEdgeJumpConversion`,
   `SurfaceVorticityConversion`) and both shedding-method families work at
   N=0, or unsupported combinations error out loudly at construction.
5. **Warmstart + replay round-trip** at N=0 (metadata records nwakerows=0;
   loading an N=0 run reproduces state; replay manifest handles the empty
   row buffer).
6. Driver: `NWAKEROWS=0` accepted and validated; banner + metadata truthful.

**Physics watch-items to note in the Log (not blockers):** fresh regularized
blobs now sit adjacent to the rigid row's singular end filament (possible
solve noise — if a quick coarse smoke shows Kutta/CT noise vs N=1, record it
and stop; do not tune); sheet-vs-particle representation of the youngest
circulation is 017's measured territory (small in tested configs).

## Test plan (extend existing suites; do not fork)

- Wake suite (`test/` — the suite that currently reports 652 tests, incl. the
  016 golden-reference and station testsets): add an N=0 testset — storage
  shapes, one full shed cycle, circulation bookkeeping vs an N=1 reference on
  identical prescribed μ, filament strength accounting, both conversion
  strategies. Keep the golden-reference N≥1 tests green and untouched.
- Simulate suite: an N=0 smoke on the existing small fixture (wing or coarse
  rotor) — runs, sheds, no NaN/σ≤0, particle count per step matches
  expectation. **R4 gotcha (016): static fixtures emit NO filament unless
  `overflowed[]` is set** — at N=0 decide and document what `overflowed[]`
  means (recommend: set true from the first shed so downstream
  filament/VTK guards behave).
- Replay suite: metadata + warmstart round-trip with nwakerows=0.
- Optional local gate smoke (only if cheap): coarse `40_40` rotor, NT=36,
  ≤4 threads, `PARTICLE_SHEDDING=sigma_overlap`, a few revs; compare CT trace
  N=0 vs N=1 qualitatively. NOTE there is no REQUIRED_REVS knob — a full
  default run is ~2–3 h; banner + first-revs capture is enough.

## Standing constraints (018-campaign rules apply)

- **Local implementation only. NO cluster deploy, no sbatch, nothing under
  `data/` committed.** Cluster deploy of `src/` is quarantined (files carry
  other-session diffs; cluster base + intended hunks only, md5 before
  overwrite) and is a separate, Ryan-gated step.
- ≤4 threads for anything local; no long local runs without asking.
- `src/FLOWPanel_wake.jl` and `src/FLOWPanel_simulate.jl` carry uncommitted
  work from live sessions (018 F1b arc code, station-σ law, 020/021/022
  session hunks). **Edit surgically; never revert or reformat surrounding
  code; run the full affected testsets before declaring done.**
- Notebook: never write it. This item's Log below is the record.
- If an invariant cannot be met without touching the rigid-row/Kutta side,
  STOP and write the blocker in the Log for Ryan instead of proceeding.

## Exit criteria

- [ ] N=0 implemented (route recorded), all invariants demonstrated by tests
- [ ] Full wake + simulate + replay suites green, N≥1 golden tests untouched
- [ ] Driver knob + metadata + banner support
- [ ] Log entry: route, filament-semantics decision, test counts, any
      physics watch-item observations
- [ ] Ryan review; cluster deploy decision (separate, not this item)

## Log

- 2026-08-20 — Item opened and staged (018 session, Ryan: "stage another
  BRAINSTORM item to allow for N=0... it's relatively simple to clean up");
  code map verified against the working tree; implementation delegated to a
  subagent under the constraints above.

- 2026-08-20 — **Implemented (Route II: age-0 conversion on N=1 storage).**
  Route I (true zero-row storage) was rejected during design, before any code:
  the converted row needs real bilinear geometry (upstream edge = current
  TE+Das line, downstream edge = the previous line convected one step) and a
  strength-history slot for the unsteady face — i.e. exactly the N=1 arrays.
  A "virtual staging row" would have re-created that storage under another
  name while forcing every consumer (FMM views, filament writer, VTK,
  warmstart) onto new code paths. Route II keeps the existing arrays and adds
  one marker: `PanelWake.convert_at_shed::Bool`.

  **Mechanism** (`shed_wake!(::PanelParticleWake)` branch, FLOWPanel_simulate.jl):
  copy strength row 1 → row 2 (the shift is empty at `nwakes[]==0`), call the
  stock `shed_wake!(pw)` (nodes row 2 ← this step's convected row-1 line,
  fresh row-1 strengths from the solved μ), `update_TE!` to pin row 1 to the
  current TE+Das line, set `overflowed[] = true`, `_convert_to_particles!`
  on the now-full one-row buffer, then reset `nwakes[] = 0`. Particles appear
  at the Das line; no free sheet ever enters a solve. Kutta Route B live-row
  reservation errors loudly in this branch (unsupported).

  **Filament / overflowed semantics** (the subtle part):
  - `overflowed[]` at N=0 means "at least one shed has happened"; set at every
    shed (per the R4 gotcha recommendation), so filament sources and VTK
    guards behave. Warmstart VTK load infers it as `idx >= 1` (metadata
    continuation restore stays authoritative).
  - With `nwakes[]==0`, `i_row = nwakes[] = 0` naturally puts the retained
    final filament on node row 1 — the TE+Das line — where it must cancel the
    rigid Das row's exposed downstream ring edge (net edge during the next
    solve = Δμ, exactly mirroring the N≥1 near-field structure). Strength:
    legacy `unsteady_filament=true` branch already reads
    `-strength[1, i_row+1, j] = -μ_last` (in bounds, correct); the handoff
    branch got an `alpha == 0` short-circuit that is bit-identical for
    i_row ≥ 1 and the only reachable handoff case at i_row = 0.
  - **Attribution constraint:** with conversion in the same solve that
    produced the row, the upstream (rigid-row) face jump is identically zero
    and the per-step unsteady circulation is the downstream face.
    `:downstream` (α=0, β=1) deposits it every shed and closes the ledger
    exactly; `:upstream` would silently drop ALL unsteady circulation (β=0
    after the first commit and a permanently-zero handoff face) and `:split`
    would strand a half-jump on a filament row that no longer exists
    (BoundsError by construction). Constructor therefore requires
    `attribution=:downstream` for the smooth conversion at N=0, and both
    `unsteady_filament=true` + `include_final_filament=true` for either
    strategy (loud ArgumentErrors). Legacy conversion needs no change: its
    unsteady filament along node row 2 carries Γ_k − Γ_{k−1} via the
    row-2 strength history.

  **Files touched:**
  - `src/FLOWPanel_wake.jl` — `convert_at_shed` field + docstrings; keyword
    (default false) and `nwakerows >= 1` guard in `PanelWake`;
    `_logical_nwakerows`; `alpha==0` short-circuit in
    `_final_filament_strength`; `write_vtk(::PanelWake)` single-node-row
    grid for convert-at-shed (nodes/velocity, no cells) so warmstart/replay
    keep the live row-1 line; `PanelParticleWake` accepts `nwakerows=0`
    (maps to N=1 storage + marker) with the validations above.
  - `src/FLOWPanel_simulate.jl` — convert-at-shed branch in
    `shed_wake!(::PanelParticleWake)` (existing path untouched below it).
  - `src/FLOWPanel_replay.jl` — manifests record `_logical_nwakerows` (0);
    reconstruction needs no change (`nwakerows=0` flows through the
    constructor). Continuation restore already round-trips: at N=0 the
    `terminal_strength` row (`nwakes[]+1 = 1`) IS the filament-carrier row.
  - `src/FLOWPanel_warmstart.jl` — `_load_panel_wake_vtk!`: skip the CellData
    lookup for single-node-row grids (`dim1 > 1` guard, unreachable change
    for N≥1 files); convert-at-shed `overflowed[]` inference.
  - `examples/rotor_hover_pressure_comparison.jl` — `NWAKEROWS=0` accepted
    with validation (requires CONVERSION=smooth + ATTRIBUTION=downstream;
    the driver's legacy mode hardcodes `unsteady_filament=false`);
    `nwakerows_extent = max(N,1)` in the Das curvature-cap/θ and
    uniform-dσ/handoff-clearance formulas (at N=0: d_front = |Das|, particles
    at the Das line); banner tags "(convert-at-shed)"; metadata records
    `nwakerows = 0` truthfully (unchanged line).
  - Tests: new testsets only — wake suite "Convert-at-shed nwakerows=0"
    (construction validation; one-real-shed particle output proven EQUAL to a
    matched N=1 `_convert_to_particles!` for both strategies; second-shed
    strength history; filament strength/count at nwakes=0; VTK
    write/warmstart-load round trip), simulate suite N=0 smoke (legacy +
    smooth, ledger closure, conversion_count = one per step), replay suite
    manifest + continuation round trip, warmstart suite full-vs-warmstart
    bitwise consistency at N=0.

  **Test counts** (all local, ≤4 threads; N≥1 golden testsets untouched):
  - wake:      652/652 before → 707/707 after (+55)
  - simulate:  160/160 before → 177/177 after (+17)
  - replay:    125/125 before → 142/142 after (+17)
  - warmstart:   26/26 before →  40/40 after (+14)

  **Semantics notes / physics watch-items (not blockers):**
  - The N=0 row is one convection step "younger" than what N=1 converts (its
    upstream edge is the just-pinned Das line, pinned with the
    post-kinematics Das of the *current* step rather than the next solve's
    possibly-refreshed Das; identical for steady RPM without DAS_REFRESH).
  - Between sheds the retained filament sits on the row-1 line frozen for the
    one solve that uses it (at N≥1 the final row's nodes convect); the
    interface advances at every shed.
  - Fresh regularized blobs now sit adjacent to the rigid row's end filament
    at the Das line — the 024 solve-noise watch-item. Not measured here: the
    optional coarse rotor smoke was SKIPPED per the item's preference (the
    suites demonstrate every invariant; a CT-level A/B belongs to the 018
    ladder, cluster-side, Ryan-gated).
  - Kutta Route B (TEAnchoredAttachment) + N=0 is unsupported and errors at
    the first shed.

  All four suites first-run green after the change (exit 0); a per-testset
  diff of the wake suite's verbose table against the pre-change baseline
  shows every pre-existing testset (incl. the 016 golden references and the
  408-test conversion-transaction set) at identical pass counts — the only
  delta is the new 024 testset. The driver was syntax-checked
  (`Meta.parseall`); its N=0 gate logic is exercised at the constructor
  level by the wake-suite validation tests. The optional coarse `40_40`
  rotor smoke was SKIPPED (item marks it optional and prefers skipping when
  the suites demonstrate the invariants — they do). No commits made; all
  changes left in the working tree. Exit-criteria checkboxes left for Ryan.

- 2026-08-20 evening — **Cluster A/B AUTHORIZED BY RYAN ("launch another
  NT=36 with λ=4 for F1b to test it") and LAUNCHED.** Driver amendment:
  legacy conversion permitted at N=0 (unsteady_filament=(N==0) in the
  legacy kwargs; the smooth-only gate was driver plumbing, not physics —
  first local smoke failed loudly on exactly that flag, then ran clean
  50+ steps at 40_40/csarc-λ4.8/N=0). Deployed under quarantine proof
  (diffs 024-only; md5s in 018 ledger §2026-08-20 evening; backups
  .deploy_backups/pre024_20260820/). **Job 13246032 =
  p018_csarc_n0_l4p8** (csarc_l4p8 + NWAKEROWS=0), 48 h/64G, SETTLE 22.
  Disclosed A/B deltas vs N=1: {no free row} + {unsteady_filament
  false→true, forced}. Smoke VTK for Ryan:
  ~/p018_das_inspect/smoke_csarc_n0_l4p8/.
- 2026-08-20 late — **N=0 NT ladder launched (Ryan order):** 13246048/49
  = csarc_n0_nt{72,144}_l4p8 (exact-rate rlxf, NT·pps=432), joining
  13246032 (NT36). Scoring lens: N=0 removes the free-row extent term
  entirely, so N=0 NT-flatness isolates pps/dt model-def from row
  extent, and N1-vs-N0 at each NT isolates the row itself.
