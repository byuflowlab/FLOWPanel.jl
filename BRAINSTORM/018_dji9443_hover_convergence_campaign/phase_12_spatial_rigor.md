# Phase 12 — Spatial-convergence rigor: h at fixed σ, σ-axis framing, d/σ clearance

Staged 2026-08-05 (Ryan-approved scope: h ladder at B0 σ with h/σ = 0.5, 0.25,
0.125; d/σ = offline kernel diagnostic + dynamic N spot-check at σ*; runs
submitted same day). Read `ops_reference.md` and `decision_rules.md` first.

## Objective

Make the campaign's spatial-convergence claim rigorous by separating three
axes that Phase 4 currently entangles:

- **(A) Particle spacing h at fixed σ** — never scored as a convergence axis
  anywhere in the campaign (`error_budget.md` had no h row before this phase).
- **(B) σ resolution** — owned by Phase 4, but its ladder co-refines h faster
  than σ (h ~ σ^q, q ≈ 1.22), so its σ claim is *conditional on* (A).
- **(C) Particle–surface kernel-regularization clearance d/σ** — the concern
  that d must be large enough that particle-on-surface induced velocity is not
  materially reduced by kernel regularization. Dynamically characterized at B0
  σ only (Phase 5); no direct kernel-level measurement exists yet.

Context for (A): the shedding machinery ties the knobs as
`h_tip = 2πR/(NT·P_PER_STEP)` and `σ_tip = OVERLAP·h_tip`
(`SigmaOverlap`, `src/FLOWPanel_wake.jl` ~684; driver σ line ~437), so
**h/σ = 1/OVERLAP exactly at the tip**, and raising `OVERLAP` and
`P_PER_STEP` by the same factor refines h at fixed σ. The only prior
fixed-σ-ish datapoint, `p018_L1_ov3` (13029923, h/σ 0.42→0.33 at L1 σ), was
scored solely on the settle-time hypothesis (refuted) and its σ drifted −1.8%;
no M1/M2 pair was ever computed. Temporal convergence, for the record, is
also NOT yet shown (Phase 3: NT=36 fails M1 by +0.52–0.60% at NT=72 at both
Das values; Richardson closure waits on `p018_nt144`) — Phase 3 owns that
axis and Phase 12 does not touch it.

## Part A — h ladder at fixed σ

Carrier = B0 exactly (NT=36, η=1.0 ⇒ Das=0.41c at 0.75R, N=4, σ/c=0.68,
velocity formulation, 4R truncation, stock relaxation, no filter). σ/c check:
1.363·OVERLAP/PPS = 0.68 on every rung. `MERGE_R_FACTOR` stays at B0's
0.0120 — it holds r_merge/σ and σ is unchanged.

| tag | OVERLAP | P_PER_STEP | h/σ | h/c | mem/time | status |
|---|---|---|---|---|---|---|
| `p018_b0` | 2.0 | 4 | 0.5 | 0.34 | done | reference: CT̄ = 0.07170, CI [0.07164, 0.07173], revs 15–29 |
| `p018_h0p25` | 4.0 | 8 | 0.25 | 0.17 | 64G / 48 h | submitted 2026-08-05 |
| `p018_h0p125` | 8.0 | 16 | 0.125 | 0.085 | 128G / 48 h | submitted 2026-08-05 |

Run plan:

- **Cold starts, `P018_SETTLE_REVS=22`** (~30 revs total, ~1080 steps) so a
  ≥15-rev M1 window (revs ~15–30) exists without chaining, matching B0's
  window. *Warm start was the plan but is unavailable*: the 2026-08-05 VTK
  sweep left B0 with no restart set (predates the retention amendment); the
  only surviving warm source is `p018_L1_ov3`@719 at σ/c 0.292, and seeding a
  fixed-σ study's wake with different-σ particles would contaminate its
  premise. Cold is clean and precedented (B0 was cold).
- Cost prior: L1 (PPS 11) ran ~45 h for 719 steps at 64G. `p018_h0p25`
  (PPS 8, ~2× B0 particles) should fit 1080 steps in 48 h or close;
  `p018_h0p125` (PPS 16, ~4× B0) likely needs a restart-chained `_s2`
  (recipe in `ops_reference.md`; last restartable step = last VTU index,
  NOT the CSV's last step). 128G per the σ-L2 memory ruling — particle
  counts are L2-like (~180k).
- Ops rules applied: new unconditional case arms (no `--export` knob
  overrides), banner verification after submission, one case per job,
  ≤10 active study jobs (these two take the count to exactly 10),
  VTK deletion after harvest keeping ONE complete restart set per run
  (disk was 174G of the 200G budget at submission — sweep promptly).

Caveats to record at harvest:

- h/σ = 1/OVERLAP is **tip-exact only**: inboard, `p_per_step =
  ceil(OVERLAP·dist/σ)` floors at small dist, so actual local h is finer than
  nominal. Log actual shed counts per station if cheap (metadata/monitors).
- σ/h = 8 at the finest rung is deep in the heavy-overlap regime flagged by
  008/011/012 (conditioning, merge workload). Record merge-event counts per
  rung — particle density ×4 at fixed r_merge/σ may change merge activity —
  but do NOT gate the verdict on them.
- Watch `max_particles=500_000` (driver hard-code): ~176k saturated expected
  at h0p125, margin OK, but verify the count actually saturates.

**Pre-registered decision rule** (written before any data): each adjacent
rung pair is scored on M1 (|ΔCT̄| ≤ 0.5%, ≥15-rev matched windows,
`scripts/p018_analyze.py m1`) AND M2 (ε_Γ,max ≤ 1% over 0.3 ≤ r/R ≤ 0.95,
`m2`), per ruling 9.

- PASS at 0.5→0.25 ⇒ h is converged at the production ratio; budget terms =
  the observed deltas; Phase 4's ladder becomes retroactively interpretable
  as a σ statement; the 0.125 rung serves as confirmation/Richardson.
- FAIL ⇒ h is a live error term: every Phase-4 Δ decomposes into σ and h
  parts, and σ* selection must move to fixed-h/σ rungs at the converged
  ratio. **Flag for Ryan before proceeding** — this reshapes Phase 4.

**Pre-registered prediction** (falsifiable): deltas small and monotone.
Discriminator: if the L1 3-lobe Γ̄ redistribution (dip at r/R ≈ 0.76)
appears here at *fixed σ*, the Phase-4 attribution of that mode to σ
(under-resolved tip vortex → 012) is wrong and it re-attributes to
h/overlap.

## Part B — σ-axis framing (no new runs)

Phase 4's ladder varies (σ, h) jointly: OVERLAP rises 2.0 → 2.4 → 2.88 down
the ladder, i.e. h ~ σ^q with q ≈ 1.22, deliberately refining h faster than
σ. That design is fine for detecting σ-axis error but cannot by itself
attribute a delta to σ rather than h. Phase 12A supplies the missing
h-invariance leg at the B0 rung. The campaign's σ-convergence claim is
therefore the conjunction: **(Phase 4 ladder verdict) ∧ (Phase 12A h-pass)**.
If 12A fails, Phase 4's deltas are joint (σ,h) statements and the σ* choice
must be revisited (see Part A fail branch). A dated pointer to this phase has
been added to `phase_04_sigma_ladder.md`'s log.

Same-day addition (Ryan): Phase 4 also gained **Ladder B — σ at FIXED h**
(OVERLAP alone at PPS 4, σ/c 0.68 → 0.477 → 0.34; range-limited by the
overlap condition h/σ ≤ 1, so not a σ→0 ladder). Together {12A (h|σ),
Ladder B (σ|h), Ladder A (joint)} form a one-factor decomposition around B0;
the shared 3-lobe discriminator is written in phase_04 §2026-08-05.

## Part C — d/σ clearance

### C0. What is already settled (do not rerun)

- d ≈ N·Das at η=1 (handoff distance); at the B0 carrier d/σ = 0.603·N·η.
- Dynamic N ladder at B0 σ: d/σ = 0.60 FAILS (−0.75% CT, ε_Γ 2.77%),
  1.21 PASSES (+0.01%, 0.35%), 2.41 = reference (Phase 5).
- Das×N cross matrix (2026-08-04): **NO d/σ collapse** — Das owns the Γ̄
  lobes; clearance penalty bites only at d/σ ≲ 1; N is a secondary modifier.
- Along the Phase-4 σ ladder at fixed physical Das, d/σ *improves*
  (0.60 → 1.38 → 2.72 in N=1 units), so fine σ is the safe side.

### C1. Offline kernel-regularization diagnostic (Mac, no cluster cost)

New script `examples/particle_surface_regularization_diag.jl` (to be written
by the executing agent; this spec is the contract): load a retained settled
wake state (a complete restart set per ruling 10 — currently
`p018_L1_ov3`@719 on the cluster; any settled state is acceptable since this
is a kernel-level measurement, but record which) plus the body; evaluate the
wake-particle-induced velocity at all body control points twice — the
regularized Gaussian-erf kernel at the particles' actual σ vs the singular
Biot–Savart limit — and report the relative deficit per control point against
its nearest-particle d/σ (plus a kernel-weighted variant). Reuse the existing
erf-kernel code paths; `examples/particle_divergence_free_check.jl` (009) has
the Taylor-branch handling near ρ→0 to crib from.

Deliverables: deficit-vs-d/σ curve; fitted exponent compared against the 017
M4 static prior (~(d/σ)^−3.4); a quantitative admissible threshold **d/σ***
defined as deficit < 0.5% of local |u|. This measures the regularization
deficit directly instead of inferring it from CT deltas.

### C2. Dynamic N spot-check at σ* — DEFERRED, gated on σ* selection (post-L2)

One N-doubling run at the final converged settings; pre-registered null:
|ΔCT̄| ≤ 0.5% AND ε_Γ,max ≤ 1%. Choose the pair with ruling 12 in mind
(prefer N=2 over N=1 when it permits significantly fewer particles per
step). This closes the loop: C1 predicts the threshold, C2 verifies it
dynamically at the operating point that matters.

## Exit criteria

- A: both rung pairs scored on M1+M2 with matched windows; verdict + budget
  rows (CT and Γ tables) filled; pass/fail branch executed.
- B: framing paragraph (above) + phase_04 cross-reference in place. DONE at
  staging.
- C1: d/σ* threshold published with the fitted curve.
- C2: executed at σ* (or explicitly re-gated if σ* slips).

## Log

### 2026-08-05 — staged and submitted

- Phase created; launcher arms `p018_h0p25` / `p018_h0p125` added
  (unconditional exports per ops rule 2); budget rows added; phase-gates row
  added to the item file.
- Cold-start deviation from the session plan recorded (B0 restart set gone;
  L1_ov3 is the only surviving warm source and is at the wrong σ).
- Submitted and **banner-verified** (mandatory ops rule), all cold with
  `P018_SETTLE_REVS=22`, η=1.0, NT=36, N=4, merge_r as specced, no filter:
  - `p018_h0p25` = job **13051772** (banner: overlap 4.0 / pps 8 /
    merge_r 0.0120 / settle 22) — 64G / 48 h.
  - `p018_h0p125` = job **13051773** (banner: overlap 8.0 / pps 16 /
    merge_r 0.0120 / settle 22) — 128G / 48 h.
  - Companion Phase-4 Ladder B rungs (see phase_04 §2026-08-05):
    `p018_ov1p4` = **13051774** (overlap 1.4 / pps 4 / merge_r 0.0084),
    `p018_ov1p0` = **13051775** (overlap 1.0 / pps 4 / merge_r 0.0060).
  - Study slots at 10/10 after these four; disk at 174G/200G at submission —
    the harvesting agent should sweep VTK promptly (keep one complete restart
    set per run).

### 2026-08-05 (b) — C1 offline kernel diagnostic EXECUTED

`examples/particle_surface_regularization_diag.jl` written and run on
`p018_L1_ov3` step 719 (Mac archive; 209,577 particles × 36,752 panels,
7.7e9 pair evals, 34 s at 6 threads). Selftest gates all PASS at machine
precision (single-particle deficit ≡ 1−g(d/σ) to 1e−16; superposition vs
naive re-sum 8e−16). Outputs in `data/particle_surface_regularization_diag/`
(per_panel.csv, binned.csv, summary.txt, deficit_vs_dsigma.png).

**Results:**

- **Admissible threshold d/σ\* (deficit < 0.5% of local wake-induced |u|):
  ≈ 2.6 by binned median, ≈ 3.4 by binned p95.** The binned curves hug the
  single-particle law 1−g(d/σ) from below (many-particle cancellation only
  helps), so 1−g is a usable conservative closed-form bound.
- Global log-log slope of the median deficit −3.29 — **consistent with the
  017 M4 static prior (−3.4)**; the tail is actually super-polynomial
  (Gaussian), so any power-law fit is range-dependent.
- Context for the campaign's criteria: B0's shedding handoff d/σ = 2.41
  (N=4) sits AT the median threshold and below the p95 one, yet the
  DYNAMIC tests (phase_05) show N=2 (d/σ=1.2) already passes M1+M2 to
  0.35%. Reading: the kernel deficit is a per-evaluation worst-case bound;
  the integrated observables (CT̄, Γ̄) tolerate ~5× more raw deficit because
  the deficit acts on the near-field of freshly shed rows whose influence
  partially cancels. Both numbers are now measured rather than argued.
- Caveat: min-d/σ per panel spans 0.011–6.6 over the whole body — the small
  end is settled wake hugging the blade surface (not the shedding handoff),
  where the relative-to-singular metric is the wrong physical yardstick
  (|u_sing| is unphysically large there); per_panel.csv carries |Δu|/V_tip
  for that regime. The curve itself is kernel-level and state-agnostic; the
  panel POPULATION along it is L1_ov3-specific (σ/c 0.292).

C1 exit criterion MET (threshold published with fit). C2 remains gated on
σ\*.

### 2026-08-05 (late) — disk pressure resolved, restart sets protected

The "sweep VTK promptly" note above is discharged: an automated, restart-set-aware
sweeper (`scripts/p018_vtk_sweeper.sh`) freed 66 GB from six CLOSED runs, taking
home from 178 G to 112 G with all four spatial-rigor runs untouched. All four are
on the protect list (`vtk_protect_list.txt`) precisely because their restart sets
are the `_s2` chain sources if the 48 h walltime cuts them short — remove a run
from that list once it has finished its 22 revs or has had its `_s2` launched.
Full manifest and the retention rules in `ledger.md` §"2026-08-05 (late)".

**VTK deletion log 2026-08-12 ~23:40 MDT (200G budget sweep, protect-list
entries removed by Ryan):** h-ladder runs swept of VTK — p018_h0p25
25544MB freed (kept restartable steps 1070–1079), p018_h0p125 20476MB
freed (kept 845–854; FLAGGED-for-Ryan run — metadata.toml + monitors
verified intact post-sweep, no CT CSV ever existed for the timeout
segment, force-monitor fallback preserved; `p018_h0p125_s2` untouched,
keeps 1070–1079). CSVs/TOML/monitors untouched throughout.
