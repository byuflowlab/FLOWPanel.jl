# 026 reset prompt — adaptive elongation shipped; verify, then commit sequencing (2026-09-11)

You are picking up BRAINSTORM 026 (resolution-preserving particle splitting)
in `/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (panel solver) +
`/Users/ryan/Dropbox/research/projects/FLOWVPM.jl` (VPM core). Read
`CLAUDE.md` and the routing policies it names before touching code.

## Documentation map (Ryan's requested format — respect it)

- **`splitting_theory.md`** (this directory) — THE theory document: a single
  living draft, always the latest design, updated IN PLACE with **no change
  history / no dated amendments**. If you change the design, rewrite the
  relevant section so the doc reads as one current draft.
- **`particle_splitting_design.md`** — the append-only decision/history
  record (§19 fractional gating, §20 separation-criteria audit, §20.7 the
  2026-09-09 rulings). Append dated sections; never rewrite.

## Where the item stands (2026-09-11, ALL UNCOMMITTED on both live checkouts)

Two uncommitted layers now stack on the live checkouts:

1. The §19 fractional-gating redesign (2026-09-08): per-mechanism fractions
   `f_visc/f_comp/f_elong` vs each particle's `sigma_0`, attempted
   (pre-clamp) Δσ² accumulators `dvisc`/`drvpm` (signed net), σ bounds as
   clamps, sigma_guard permanent (commit 8 CANCELLED).
2. The 2026-09-09 separation-criteria work (audit = design doc §20; rulings
   = §20.7; theory = `splitting_theory.md`):
   - **Adaptive in-line elongation kernel** in
     `FLOWVPM.jl/src/FLOWVPM_resolution_split.jl`: `_elongate_plan` gives
     `m = clamp(round(λ_att·σ₀/σ_c), 2, elongate_m_max)` and spacing
     `s = λ_att·(σ₀/Φ_t)/m`, with `λ_att = σ₀²/(σ₀² + drvpm)` (σ²L = const
     is EXACT for the rVPM channel — proof in theory doc §2);
     `_split_elongate_line!` places centroid-symmetric children,
     `Γ_p/m`, `σ_c = σ_p`, circulation unchanged. Saturation (m_ideal >
     m_max or attempted total collapse) falls back to spacing `σ_p/Φ_t`
     (retain overlap, under-tile). Properties: m·σ_c³ = σ₀³ (volume rule =
     merge inverse), exact composition, no per-particle L₀ state needed.
   - New `ResolutionSplitOpts` fields: `elongate_overlap` (Φ_t; **NaN =
     FLOWVPM default = legacy fixed pair2** at `elongate_offset_ratio`) and
     `elongate_m_max` (default 4). `split_particles!` return gained
     `n_children_elongate`.
   - **circ bookkeeping fix**: tri3/tetra4 children carry `circ/3`/`circ/4`
     (lengthwise bundle division); elongation keeps circ (crosswise cut).
   - **Compression stays m=3** (Ryan: m=2 too anisotropic); docstring +
     theory §5 note it is count-matched to f_comp ≈ √3−1 ≈ 0.73 and that
     higher m needs a more complicated child shape (future).
   - Driver (`examples/rotor_hover_pressure_comparison.jl`):
     `WAKE_SPLIT_ELONGATE_OVERLAP` (defaults to the shedding `OVERLAP`, so
     the DRIVER default is adaptive; NaN selects legacy pair2 as the
     comparison arm), `WAKE_SPLIT_ELONGATE_M_MAX` (default 4), and
     `MERGE_OVERLAP` (default NaN = OFF → production keeps the absolute
     0.02R merge radius; when set, flips MergeParticles to
     `sigma_relative=true, r = r_hash = 1/Φ_merge` — merging's native
     overlap criterion, no FLOWVPM change). Provenance dump prints
     `merge_overlap`.
   - Angular-impulse claim CORRECTED everywhere: in-line children are NOT
     exact by ± symmetry; error ≤ a²|Γ|/3 with a = (m−1)s/2, zero when
     axis ∥ Γ (t2 bound).

Key physics facts (derived + code-pinned, see theory doc §1–§4):
- Overlap convention Φ ≡ σ/h (`FLOWPanel_wake.jl:730,2113`); campaign 2.75.
- Production merging (MERGE_R_FACTOR=0.02, sigma_relative=false) is an
  ABSOLUTE radius 0.02R = 0.525σ* — it exceeds shed spacing σ*/Φ = 0.0139R,
  so merging (not OVERLAP) sets the wake's effective spacing floor, and
  split→merge cycles are a σ-PUMP (merge σ = ∛Σσ³ → ×m^(1/3) per cycle).
  Overlap-gated merging with Φ_merge = 3.5 > Φ_t = 2.75 kills this by
  construction; it is OFF by default pending an A/B because it merges
  strictly less than production (counts/cost rise).

## Verification status

GREEN: FLOWVPM `test/runtests_resolution_split.jl` 1011/1011 (includes new
t10 adaptive-elongation testset — plan math, saturation, conservation +
in-line geometry + angular bound, EXACT composition, capacity, legacy
fallback, opts validation — and t11 circulation shares; note `NTzero`-style
literals now include `n_children_elongate=0`). FLOWPanel
`runtests_unit_wake.jl` 730/730, `runtests_unit_replay.jl` 148/148.
Known unrelated pre-existing failure: `test/runtests_unit_warmstart.jl`
FIRST testset (immutable `WarmstartNoopSolver` as WeakKeyDict key) — not
ours.

INCOMPLETE: the sane-fraction driver smoke
(`WAKE_SPLIT_STRETCH=true WAKE_SPLIT_FRAC_COMPRESS=0.5
WAKE_SPLIT_FRAC_ELONGATE=0.3 SIGMA_CEIL=1.0 NREVS=0.2 RUN_MONITORS=false
julia --project -t 4 examples/rotor_hover_pressure_comparison.jl`) was
killed at step 294/467 (~61 s) — but everything observed was healthy:
adaptive kernel firing (`elongate=1 (children=3)` per event, exactly the
predicted m=3 at f_elong=0.3), compress events firing, zero capacity skips,
zero errors. Log (may be gone after reboot):
`/private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/8387cdae-00c1-4884-bee7-fc5eb0b48bab/scratchpad/smoke_split_adaptive.log`.

## Your tasks

1. **Finish verification** (local, ≤ 4 threads):
   a. Re-run the sane-fraction smoke above to completion; confirm sane
      counters and a clean exit.
   b. Run a second smoke adding `MERGE_OVERLAP=3.5` to exercise the
      overlap-gated merge knob end to end (expect: more particles retained
      than the default-merge smoke, no split/merge churn — split children
      at Φ_t are strictly below the merge threshold by construction).
   c. Optionally one legacy-arm smoke with
      `WAKE_SPLIT_ELONGATE_OVERLAP=NaN` to confirm the comparison arm still
      runs (2 children per elongate event).
2. **Commit sequencing — coordinate with Ryan (AskUserQuestion) before
   committing anything.** The uncommitted stack spans both repos: §19
   fractional gating + exposure removal, then the 2026-09-09 layer (theory
   doc, §20/§20.7, adaptive kernel + tests, circ fix, driver knobs).
   Propose a commit breakdown (likely: FLOWVPM kernel+tests / FLOWPanel
   driver+docs, §19 layer first) and let Ryan rule.
3. **Then the item returns to the §19 pending queue**: commit-7 campaign
   arms (§8.4 cap030/cap018, §9 s020v matrix) are Ryan-gated with fractions
   to re-derive in fraction space (f_comp near the §5 adequacy table;
   f_elong / D4 re-derivation) — now also choose per-arm
   `WAKE_SPLIT_ELONGATE_OVERLAP` / `MERGE_OVERLAP` (an A/B of
   MERGE_OVERLAP=3.5 vs absolute is the natural first discriminator, given
   the σ-pump finding). Campaign launches: worktrees + tags per the global
   campaign policy, Ryan-gated.

## Open questions on the books (unchanged)

- D2 (tri3 ring radius 0.6 vs moment-match 1.155) and D3 (tetra4 geometry)
  — both gated on the §3a kernel-fit study; tetra4 has zero production
  events.
- Higher-m compression child shapes (ring+center / two rings) — noted as
  future work in theory §5, deliberately not pursued.
- GPU gap: device-resident integrator twins still don't maintain
  dvisc/drvpm → splitting remains host-mirror-only.

## Ground rules (carry-over)

Local runs ≤ 4 threads. Commit only when Ryan says so. Campaign launches
Ryan-gated. Theory doc is rewritten in place; design doc is append-only.
NaN-disabled knob defaults; VTK/warm-start untouched unless per-particle
state is added (none was — the adaptive kernel needs no new state).
