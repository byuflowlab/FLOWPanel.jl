# 022 decision rules

Adopted from 018 (`018_dji9443_hover_convergence_campaign/decision_rules.md`)
with ground-specific additions. Pre-registered before any production run.

## M1 — windowed CT mean (primary)

- CT series = negated axial Bernoulli force channel (thrust convention),
  identical to 018.
- Window: the final complete revolutions of self-sustained hover (window
  entirely after freestream withdrawal, `window_in_hover=true`), ≥ 10 revs.
- Headline number = cycle-mean ± cycle-std over the window (limit-cycle-safe),
  as computed by the driver's per-rev block code (fork of RHPC Phase-2e code).
- IGE runs may settle slower (near-ground recirculation): judge window entry
  by per-rev block stationarity, not a fixed rev count. If per-rev means are
  still trending monotonically at end-of-run, the run is NOT harvestable —
  chain or extend before quoting M1.

## Headline comparison (Phase 2)

- Report CT_IGE, CT_OGE, and the ratio CT_IGE/CT_OGE per mesh rung, with
  cycle-std error bars propagated in quadrature on the ratio.
- Sanity anchor (direction + rough magnitude only, NOT an acceptance
  criterion): momentum theory / Cheeseman–Bennett at h/R=1 ⇒ ratio ≈ 1.07 at
  fixed power. A ratio < 1 or > ~1.2 demands diagnosis before reporting.

## Ground-convergence acceptance (Phase 3)

- One axis at a time off the first-light carrier (disc radius, panel length,
  truncation radius).
- ACCEPT a rung as converged when doubling the knob (halving for panel
  length) moves the M1 cycle-mean by **< 0.5%** AND the change is within the
  quadrature-combined cycle-std of the two runs. Adjustable only at a Ryan
  gate.

## Health gates (every run, before any M1 is quoted)

- Startup banner verified: mesh key/file, shedding edge counts (no
  circumferential edges), RPM/ρ/R, ground knobs, GS settings. Exit 0 ≠ health.
- GS outer loop: `gs_nonconverged = 0` in metadata; `gs_iters_max` well below
  `GS_MAX_OUTER` (default cap 50). Any non-converged solve invalidates the run.
- Ground tangency: RMS(U·n) on the ground stays bounded and small relative to
  tip speed (report `ground_tangency_rms_max / (ΩR)`; smoke establishes the
  expected scale). A step-change or growth trend = investigate.
- Below-ground census: report count and Σ|Γ| trajectories; leave-be is
  acceptable while Σ|Γ| below ground stays a small fraction of the total wake
  and CT shows no correlated drift.
- Wake health / sigma columns: judge absolutes, not ratios (018 lesson:
  sigma-ratio columns can be NaN).

## M2 (secondary, diagnostic)

- Γ̄(r/R) spanwise circulation profile from BoundCirculationMonitor: overlay
  IGE vs OGE. No acceptance threshold in 022; used to attribute WHERE the
  ground effect loads the blade (expected: inboard-loaded shift + tip-vortex
  standoff change).
