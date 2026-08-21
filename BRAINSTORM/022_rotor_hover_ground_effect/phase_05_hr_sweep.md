# Phase 5 — h/R sweep vs experiment and theory

**Objective:** T/T_OGE vs h/R at 0.5, 1.5, 2.0 (h/R = 1.0 reused from Phase 2;
T_OGE denominator from Phase 1 at the **same mesh rung**), compared against the
experimental table below and the Cheeseman–Bennett curve. Validation phase —
capped at h/R = 2.0 per Ryan (effect is gone in experiment by 2.0).

## Gate

Blocked until the Phase-2 h/R = 1 case is harvested and healthy (banner + GS
converged + tangency bounded + M1 window per `decision_rules.md`) AND its
IGE/OGE ratio is judged sane against the ≈1.07 momentum-theory anchor. If
Phase 3 has selected a converged ground prescription by then, inherit it;
otherwise run on the Phase-2 carrier (4R disc / 0.15R panel / 3R truncation)
and say so in `ledger.md`.

**Gate status 2026-08-19: NOT MET, and the coarse route to it is gone.**
p022_ige_coarse (13207682) blew up at rev 17.6 and was cancelled; it has no
harvestable M1 window. The gate now rests on the fine rung p022_ige_fine
(13207681, still running) plus p022_oge_fine (13207679). Two consequences the
executing agent must weigh before submitting anything:

1. The blow-up is an unexplained IGE instability at the *same* h/R the sweep
   brackets. Running h/R = 0.5 — a stronger-recirculation case — before
   understanding it risks burning three 48 h jobs on the same failure mode.
   Check what p022_ige_fine does around its own rev 17–18 first.
2. If Phase 4's cull ablation is pulled forward in response (see
   `phase_02_ige_first_light.md`), the particle policy may change from `none`,
   which would change the Phase 5 carrier. Do not submit the sweep on a carrier
   that Phase 4 is about to invalidate.

## Reference data

Experiment (Ryan-provided 2026-08-19; source citation TBD — ask Ryan at
write-up). Full table recorded for reference; the sweep covers h/R ≤ 2.0.

| $h/R$ | 0.5 | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 | 4.0 | 5.5 | 7.0 |
|---|---|---|---|---|---|---|---|---|---|
| $T/T_{OGE}$ | $1.200 \pm 0.009$ | $1.078 \pm 0.008$ | $1.010 \pm 0.014$ | $1.004 \pm 0.004$ | $0.996 \pm 0.009$ | $0.999 \pm 0.023$ | $1.009 \pm 0.008$ | $0.998 \pm 0.012$ | $1.001 \pm 0.012$ |

Cheeseman–Bennett (fixed power):

$$
T/T_{OGE} = \frac{1}{1 - (R/4h)^2}
$$

giving 1.333 (h/R = 0.5), 1.067 (1.0), 1.029 (1.5), 1.016 (2.0). CB is known
to overpredict at small h — the experiment's 1.200 at h/R = 0.5 vs CB's 1.333
is expected, not a red flag.

## Cases (pre-registered)

One row per run; coarse rung 56_57, **48 h** walltime, all carrier knobs
identical to Phase 2 except the two below. Case arms are already written into
`examples/run_rotor_ground_effect_hpc.slurm.sh` (not submitted).

`TRUNCATION_DEPTH_R = GROUND_H_R + 3` holds the below-ground allowance fixed
across the sweep. Depth is referenced to the **rotor**, not the ground: the
truncation cylinder starts 0.5R upstream of the rotor plane and extends
`TRUNCATION_DEPTH_R`·R downstream, and the rotor plane itself sits at
axial/R = −0.0425, so

$$
\text{allowance} = (\texttt{TRUNCATION\_DEPTH\_R} - 0.5) - (h/R - 0.0425)
= \texttt{TRUNCATION\_DEPTH\_R} - h/R - 0.4575 .
$$

With depth = h/R + 3 this is **2.5425 R for every case — identical to the
Phase-2 h/R = 1 carrier**, which is the actual matching requirement (the round
number "3R" was the intent; 2.5425R is what the carrier actually uses).
Verified directly in the local smokes below.

| case | GROUND_H_R | TRUNCATION_DEPTH_R | ground at axial/R | note |
|---|---|---|---|---|
| p022_hr05 | 0.5 | 3.5 | 0.4575 | strongest effect; watch wall jet vs 4R disc (VTK) |
| p022_hr15 | 1.5 | 4.5 | 1.4575 | |
| p022_hr20 | 2.0 | 5.0 | 1.9575 | wake reaches ground later → judge settle per-rev blocks, not a fixed count |

The `TRUNC_RADIUS_R ≤ GROUND_RADIUS_R` constraint (Phase 3) carries over
unchanged (3.0 ≤ 4.0; the driver's warning stayed silent in both smokes).

### Why 48 h and not 24 h

Measured from the healthy OGE coarse run 13207680: ~87 s/step averaged over its
first 670 steps, ~110 s/step recently ⇒ the 1007-step schedule needs **25–31 h**.
The 24 h coarse walltime used in Phases 1–2 is undersized (13207680 will wall
around rev 25 of 28, short of the ≥10-rev settled window). Warm-start with
`GROUND_ENABLE=true` errors out by design
(`examples/rotor_hover_ground_effect.jl:1323`, untested multi-body restart
state), so **a walled run cannot be chained — only re-run cold**. That makes
the first submission's walltime load-bearing. h/R = 0.5 should be budgeted at
the top of the range or beyond: stronger recirculation retains more particles.

### Local pre-flight (2026-08-19, discharges ruling 3 for the new heights)

40_40 / NT=6 / 18-step smokes at 4 threads, `BERNOULLI_ONLY=true`:

- **h/R = 0.5** — exit 0. Ground built at axial/R = 0.4575 (4752 panels, 2463
  nodes), truncation depth 3.5R, allowance 2.5425R. **GS: 0 non-converged
  solves.** No `TRUNC_RADIUS_R` warning. Tangency RMS(U·n) 1.6e-4 → 4.3e-3 m/s
  over the first 10 steps.
- **h/R = 2.0** — exit 0. Ground built at axial/R = 1.9575, truncation depth
  5.0R, same 4752-panel disc, allowance 2.5425R. **GS: 0 non-converged
  solves.** No `TRUNC_RADIUS_R` warning. Tangency RMS(U·n) 6.6e-6 → 1.3e-3 m/s
  over the first 10 steps — an order of magnitude below the h/R = 0.5 case at
  the same step, as expected from the larger standoff.
- Gotcha for anyone rerunning these: the ground body requires
  `BERNOULLI_ONLY=true`. `PressureLaplace` is constructed for a single body and
  throws `ArgumentError: PressureLaplace was constructed for 1 bodies, got 2
  bodies` on the two-body IGE path. Production already sets this in the
  launcher's shared block; a hand-rolled smoke that omits it fails at the first
  monitor call, not at setup.

**Mesh-rung policy (pre-registered):** sweep on the coarse rung. If the
Phase-2 coarse vs fine IGE/OGE **ratio** at h/R = 1 agrees within combined
cycle-std, coarse ratios stand for the sweep; otherwise promote h/R = 0.5
(largest signal) to the fine rung as a confirm run.

## Preliminary read and resolvability (2026-08-19)

**Preliminary matched read — NOT an M1 harvest.** From the pre-blow-up portion
of p022_ige_coarse against p022_oge_coarse over revs 7–17 (post-withdrawal but
inside the settling transient; `scripts/p022_harvest.py ratio`):

| run | CT̄ | cycle-std | SE of mean |
|---|---|---|---|
| p022_oge_coarse | 0.075516 | 0.003093 (4.10%) | 0.000978 (1.30%) |
| p022_ige_coarse | 0.078888 | 0.005376 (6.82%) | 0.001700 (2.16%) |

Ratio = **1.045** ± 0.083 (cycle-std propagated) / ± 0.026 (SE propagated).
Direction is right and it is consistent with both the anchor (1.067) and the
experiment (1.078 ± 0.008) — but it does not discriminate between them.

**Resolvability.** This is the load-bearing finding for Phase 5. The per-rev
scatter is currently 4–7%; even at a settled ~1.4% per-rev std (018's converged
scale) over a 10-rev window the ratio standard error is ~0.6%. Against that
floor:

| h/R | experimental signal (T/T_OGE − 1) | vs ~0.6% ratio SE | verdict |
|---|---|---|---|
| 0.5 | +20.0% | ~33σ | clearly resolvable |
| 1.0 | +7.8% | ~13σ | resolvable |
| 1.5 | +1.0% | ~1.7σ | marginal |
| 2.0 | +0.4% | ~0.7σ | **below the floor** |

Consequences, pre-registered so the harvest is not read as more than it is:

- h/R = 0.5 and 1.0 carry the quantitative validation.
- h/R = 1.5 can support a directional statement at best; do not quote a
  discrepancy there as evidence of a modelling error.
- h/R = 2.0 can only support "consistent with no ground effect". That is still
  worth reporting — it is exactly what the experiment says (1.004 ± 0.004) —
  but it is a null result, not agreement to 0.4%.
- Reaching a genuine 0.4% resolution would need ~280 revs at the current
  scatter. Not feasible; do not attempt it by lengthening runs alone.

**OPEN QUESTION — raise with Ryan after the runs land (his instruction,
2026-08-19: record now, design later).** How to get a smoother CT signal so the
small-effect end of the sweep becomes resolvable. Candidate directions, none
selected: longer settled windows; paired-seed / common-mode cancellation
between the IGE and OGE members of a pair; harvesting the ratio per-rev rather
than as a ratio of independently-averaged means (correlated transients would
partly cancel); or accepting the null-result framing at h/R ≥ 1.5.

## Comparison table (fill at harvest)

Sim uncertainty = cycle-std of numerator and denominator propagated in
quadrature.

| h/R | CT_IGE (sim) | ratio sim | experiment | CB theory | sim − exp |
|---|---|---|---|---|---|
| 0.5 | — | — | 1.200 ± 0.009 | 1.333 | — |
| 1.0 | (Phase 2) | — | 1.078 ± 0.008 | 1.067 | — |
| 1.5 | — | — | 1.010 ± 0.014 | 1.029 | — |
| 2.0 | — | — | 1.004 ± 0.004 | 1.016 | — |

## Decision (pre-registered)

This is a validation phase, not a knob gate — report the curve. Flag any point
where |sim − exp| exceeds the experiment error bar plus the propagated sim
std. Qualitative pass: a monotonically decaying sim curve that lands inside
the experimental bars by h/R = 2.

## Notes

- h/R = 0.5 may develop stronger recirculation/fountain flow: if the tangency
  monitor or GS outer-iteration count degrades vs the h/R = 1 case, log it and
  consult the item-file contingency chain before touching knobs.
- Below-ground particle policy stays `none` (ruling 2) unless Phase 4 has
  ruled otherwise by run time.

## Tooling (built 2026-08-19, ready to use)

- `scripts/p022_harvest.py` — M1 cycle statistics, matched IGE/OGE ratios, and
  the figure CSVs. Subcommands `m1`, `ratio`, `figdata`. It exists because
  `scripts/p018_analyze.py` pins RPM 5400 in its force-monitor fallback with no
  CLI pass-through, and that fallback is the load-bearing path for 022 (the
  end-of-run `_CT_vs_rev.csv` is never written for a walled or cancelled run).
  Guards worth knowing: it refuses to average |CT| > 1 steps and reports them,
  it warns when the requested window is only partly covered, and it labels any
  window starting before rev 18 as preliminary.
- `figures/fig_hr_sweep.tex` + `figures/fig_hr_sweep/{experiment,theory_cb,sim}.csv`
  — standalone TikZ, pdflatex-verified. Regenerate the CSVs with
  `python3 scripts/p022_harvest.py figdata --point <h/R> <ige_run> <oge_run> <lo> <hi>`
  (repeat `--point` per case), then recompile from the `figures/` dir. Currently
  carries the single preliminary h/R = 1 point, annotated as such.
- Launcher arms `p022_hr05` / `p022_hr15` / `p022_hr20` in
  `examples/run_rotor_ground_effect_hpc.slurm.sh`, with the 48 h sbatch lines in
  the header usage block. **Not deployed to the cluster and not submitted.**

## Exit criteria

Comparison table filled + ratio-vs-h/R figure (sim, experiment with error
bars, CB curve) recorded in `ledger.md` and the item file.

## Log

- 2026-08-19 (staging): phase staged per Ryan request; experiment table
  received; gated on Phase-2 h/R = 1 harvest (jobs 13207681/82 in flight).
- 2026-08-19 (pre-work, gate still closed): tooling above built and verified;
  walltime corrected 24 h → 48 h from measured step rate; local h/R = 0.5 and
  2.0 geometry smokes passed (GS 0 non-converged, allowance invariant confirmed
  at 2.5425R); preliminary matched read and the resolvability analysis
  recorded. p022_ige_coarse blow-up documented here and in
  `phase_02_ige_first_light.md` — it removes the coarse route to the gate.
