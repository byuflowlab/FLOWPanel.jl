# Phase 13 — Panel→particle handoff clearance (`d/σ`), and its coupling to σ and dt

Status: **S0c MEASURED 2026-08-05.** Zero HPC cost. Tooling landed; findings below are arithmetic over the
production mesh geometry and are not subject to run-to-run scatter.

## 1. What the handoff distance is, and why it was never controlled

The distance from the trailing edge at which the panel wake hands off to particles is **not a knob**. It is
emergent from `Das`, `nwakerows`, and `dt`. Geometry (verified in `src/FLOWPanel_simulate.jl:1212-1300` and
`src/FLOWPanel_wake.jl:2079-2119`): node row 1 is rigidly re-pinned at `TE + Das` every step; rows 2..N convect
freely and are separated by one step of TE travel (the `Das` cancels between consecutive rows); particles are
deposited on the segment spanning node rows `N → N+1`. Hence, per shedding station `j`:

$$ d_{\text{front},j} = |Das_j| + (N-1)\,|\boldsymbol{\omega}\times\mathbf{r}_j|\,\Delta t $$

There is no length-based or σ-based conversion criterion anywhere in `src/`; the trigger is the integer
buffer-fullness test `pw.nwakes[] >= n_rows-1` (`FLOWPanel_simulate.jl:1290`).

**Tooling (S0c).** `examples/rotor_hover_pressure_comparison.jl` now prints the per-station `d/σ` profile from
the run's own geometry and optionally dumps it (`RHPC_HANDOFF_CSV`), with `RHPC_HANDOFF_PROFILE_ONLY=true` for
a seconds-long query. It is placed **before** `Backslash(rotor)`, which eagerly assembles and LU-factors a
dense `ncells²` matrix (~10 GB at `45_185_ct4`) — so the profile costs mesh load only. Reuses the driver's own
station indexing (the same `cells[shedding[3,j], shedding[1,j]]` nib/nia convention as `station_chords`), so it
cannot drift from what the solver actually does.

## 2. The measurement: B0 violates the C1 clearance bound inboard by ~4×

Production carrier (`45_185_ct4`, NT=36, RPM 5400, η=1.0 kinematic `Das`, N=4, OVERLAP 2.0 / PPS 4):

```
sigma = 0.010385 m  (sigma/R = 0.08727, one-step tip travel = 0.020769 m)
shedding1: d/sigma min 0.713 at r/R 0.100; max 6.787 at r/R 0.998; median 3.74
shedding2: identical  (clean y -> -y symmetry check)
```

Against Phase 12 C1's measured admissible `d/σ* ≈ 2.6` (median) / `3.4` (p95): the tip is comfortable at 6.79,
but **the inboard stations sit at 0.713 — roughly a quarter of the admissible clearance.** The blade's
innermost shed stations are therefore evaluating a particle field whose singularities are well inside the
kernel-deficit radius.

This is a single-number explanation for a pattern the campaign has hit repeatedly and attributed separately
each time: the **inboard-localized** Γ̄ errors of the `N=1` case (−2.8% at r/R 0.31), of small `Das`, and of the
σ ladder. All three are small-`d` cases, and `d/σ` is worst inboard by construction.

## 3. Root cause of the span dependence: chord and travel run in opposite directions

Measured on the production mesh:

| | root (r/R 0.10) | tip (r/R 1.00) | ratio |
| --- | --- | --- | --- |
| local chord / R | 0.1680 | 0.0723 | 0.43× |
| one-step travel / R | 0.0174 | 0.1742 | 10.0× |
| travel per chord | 0.104 | 2.41 | **23×** |

Each additional free wake row adds **0.10 local chords at the root but 2.41 at the tip**. So `d` cannot be made
uniform in *either* normalization by choosing a scalar `N`: `d/σ` rises steeply outboard, while `d/c` rises
even more steeply. This is the same structural pathology the campaign already diagnosed for `Das` itself
("inboard needs η≥2.8, tip needs η≤2.3, windows disjoint") and solved with chord-proportional `Das` — but the
free-row travel term is ∝ r and is *not* fixed by that change.

## 4. Two-sided window: the clearance bound binds at the root, the sheet bound binds at the tip

Lower bound: `d/σ ≥ 3.4` (C1 clearance) — binds inboard. Upper bound: total sheet length in local chords,
`d/c` — binds outboard, and is the mechanism behind the campaign's measured "large `Das` ⇒ **outboard**
deficit (rigid sheet vs rolled-up tip)". Sweeping κ (`DAS_CHORD_FRACTION`), `N`, and σ:

| σ/R | κ | N | min `d/σ` (root) | tip `d/c` | clears? |
| --- | --- | --- | --- | --- | --- |
| 0.0873 (B0) | η=1.0 | 4 | 0.713 | — | no |
| 0.0873 | 0.41 | 4 | 1.388 | 7.64 | no |
| 0.0873 | 1.00 | 8 | 3.322 | 17.86 | no |
| 0.050 | 1.00 | 2 | 3.709 | 3.41 | yes |
| 0.030 | 0.41 | 4 | 4.038 | 7.64 | yes |
| **0.020** | **0.41** | **2** | **4.316** | **2.82** | **yes** |
| 0.020 | 0.41 | 4 | 6.058 | 7.64 | yes |

Two structural results:

1. **At B0's σ the window is empty.** Reaching `d/σ ≥ 3.4` at the root requires N=18, which drives the tip to
   `d/σ = 34.7` and `d/c ≈ 40` — i.e. replacing essentially the whole near wake with sheet. No `N` fixes B0's σ.
2. **σ = 0.02R with κ = 0.41 and N = 2 is the best available configuration** — it clears the root bound
   (4.32) at the *smallest achievable tip sheet* (2.82 chords) of any clearing option. `N=1` does not work:
   with no free rows `d = κ·c`, whose minimum is at the **tip** (chord is smallest there), giving `d/σ = 1.48`.
   So N=2 is the minimum viable row count, which is also consistent with ruling 12 (prefer N=2 over N=1).

Caveat, stated rather than assumed: the *upper* bound is not quantified. 014/017's admissible band (0.25–1.5
local chords) was measured for the **rigid row-1 `Das`**, whereas rows 2..N convect freely and the
`freestream_convection` study showed the panel wake converges to the semi-infinite reference. So the free-row
sheet is a better wake model than a rigid one and the 1.5c bound should not be applied to `d` directly. If it
*did* apply, the window would be empty at every σ tested (min achievable tip `d/c` = 2.82). **Quantifying the
upper bound is the open follow-up.**

## 4b. Per-station `Das` makes `d/σ` exactly uniform — and it favours SMALL N

§3–§4 are conditioned on the two *implemented* `Das` span laws (η·travel ∝ r, or κ·c), neither of which
compensates the travel term. Ryan's correction: with `Das` free per station, uniformity is a construction,

$$ |Das_j| = D\,\sigma - (N-1)\,|\boldsymbol{\omega}\times\mathbf{r}_j|\,\Delta t $$

and the machinery already exists — `_set_Das_station_lengths!` (`FLOWPanel_simulate.jl:28-52`) accepts an
arbitrary per-station length array, and the driver already builds one for chord mode
(`rotor_hover_pressure_comparison.jl:596-610`). **Driver-level change; no `src/` work.**

Feasibility is set at the tip, since `Das_j` must stay positive and inside 014/017's 0.25–1.5 local-chord band:
`D ≥ [(N−1)·travel_tip + 0.25·c_tip]/σ`. With `travel_tip = 8.7σ` at σ = 0.02R, each added free row forces `D`
up by ~8.7 and pushes the *root* `Das` out of band:

| σ/R | N | required uniform `d/σ` | `Das`/c (tip → root) | tip `d/c` | verdict |
| --- | --- | --- | --- | --- | --- |
| **0.020** | **1** | **3.40** | **0.94 → 0.40** | **0.94** | **in band** |
| 0.020 | 2 | 9.68 | 0.27 → 1.05 | 2.68 | in band |
| 0.020 | 3 | 18.32 | 0.25 → 1.97 | 5.07 | `Das` band violated at root |
| 0.030 | 1 | 3.40 | 1.41 → 0.61 | 1.41 | in band |
| 0.030 | 2 | 6.45 | 0.27 → 1.05 | 2.68 | in band |
| 0.050 | 1 | 3.40 | 2.35 → 1.01 | 2.35 | violated at tip |
| 0.0873 (B0) | 1, 2, 3 | — | — | — | **violated at every N** |

Three consequences.

1. **The per-station construction rescues small N, not large N.** Uniformity is free at N=1 and expensive by
   N=3.
2. **The recommended configuration changes** from §4's "σ=0.02R, κ=0.41, N=2" to **σ = 0.02R, N = 1,
   `Das` = 3.4σ span-uniform in absolute length**: exactly uniform clearance at the C1 bound, `Das`/c inside
   band at both ends, and total sheet **under one local chord everywhere** (vs 2.68 at N=2, 7.64 at the
   current N=4).
3. **It removes the premise of ruling 12.** That ruling preferred N=2 over N=1 because `N` bought `d/σ`
   clearance a finer σ would otherwise pay for; here `Das` buys it directly. This also reframes the recorded
   N=1 failure (−0.75% CT, ε_Γ 2.77%), which was measured at `d/σ = 0.60` — confounded by clearance, not a
   verdict on N=1.

**Open caveat.** At N=1 there are no freely convecting rows: row 1 is rigidly re-pinned at TE+`Das` each step,
so the whole handoff region is rigid sheet — the representation the campaign found least faithful outboard
(`nwakerows` 1→4 = +4.2%). So N=1-uniform wins on clearance and sheet length but is untested on rigidity;
**N=2-uniform is the hedge, and the two form a clean A/B** at matched, uniform `d/σ`. This is the same
unquantified upper bound as §4, now sharpened into a concrete experiment.

### Why B0's σ fails, and how strongly

The two constraints are in **different units**: clearance is a multiple of σ (`d ≥ 3.4σ`), the `Das` band a
multiple of *local chord*. At B0, `3.4σ = 35.3 mm = 4.10 tip-chords = 1.77 root-chords` — already outside the
0.25–1.5c band before any row is added. Adding rows moves `d` from `Das` into travel, but travel is worth
**2.41 chords per row at the tip and 0.104 at the root**. Admissible window on `d`:

| N | window on `d` | as `D = d/σ` at B0 σ | verdict |
| --- | --- | --- | --- |
| 1 | [0.00788, 0.01291] m | [0.76, 1.24] | non-empty, caps below 3.4 |
| 2 | [0.02304, 0.03207] m | [2.22, 3.09] | non-empty, misses 3.4 by **10%** |
| 3 | [0.04361, 0.03414] m | lo > hi | **empty** |
| 4 | [0.06434, 0.03621] m | lo > hi | **empty** |

Two distinct failures. At N ≤ 2 the band is satisfiable but the largest permitted `d` is short of `3.4σ`. At
N ≥ 3 the window **inverts**: the tip's lower requirement `d ≥ 0.25c_tip + (N−1)·travel_tip` exceeds the root's
upper allowance `d ≤ 1.5c_root + (N−1)·travel_root`, crossing over at
`(N−1)(travel_tip − travel_root) > 1.5c_root − 0.25c_tip` ⇒ N ≥ 3. Beyond two rows the tip demands `d` grow
and the root demands it not.

**Strength of the claim.** "No solution at any N" is structurally true only for N ≥ 3; at N = 1–2 it is a
quantitative shortfall. And the N=2 shortfall is **only 10%**, measured against the 1.5c cap — the softest
number here, taken from 014/017's *rigid row-1* study and explicitly unquantified for convecting rows. A 10%
margin does not survive that uncertainty. **Honest statement: B0's σ is disfavoured, not excluded; σ = 0.02R
clears comfortably (0.94 tip-chords) while B0 sits at the edge of a criterion not yet pinned down.** This makes
the §6.4 upper-bound item load-bearing for the σ recommendation rather than a side quest.

## 5. The dt pin works exactly, and un-confounds the existing NT=72 result

Holding `Das` dt-independent (chord-fraction mode, or η ∝ NT) and scaling `N` so the *whole profile* is
invariant:

| NT | Δt | N now | min `d/σ` now | **N pinned** | min `d/σ` pinned | tip `d/σ` pinned |
| --- | --- | --- | --- | --- | --- | --- |
| 36 | 1× | 4 | 0.713 | **4** | 0.713 | 6.786 |
| 72 | ½× | 4 | **0.414** | **7** | 0.713 | 6.786 |
| 144 | ¼× | 4 | **0.264** | **13** | 0.713 | 6.786 |

The pin is exact — min *and* tip match to three decimals across all three rungs, because `Das` is dt-fixed and
the travel term scales as `(N−1)·Δt`. Without it, refining dt at fixed `N` **degrades** inboard clearance by
1.7× (NT=72) and 2.7× (NT=144).

**Consequence for the campaign's dt axis:** `p018_nt72_eta2` measured +0.52% against B0 and FAILED M1. That
rung ran at N=4, i.e. at min `d/σ = 0.414` versus B0's 0.713 — so it is **not a clean dt refinement**; it
confounds temporal truncation with a 1.7× loss of handoff clearance, in a regime already 4× below admissible.
The {36, 72, 144} Richardson close-out cannot be interpreted until the rungs are re-run pinned. At N=2 base
(the §4 recommendation) the pinned ladder is N = 2 / 3 / 5.

## 5b. Under the §4b construction the dt pin is automatic, and at N=1 it is exact

`Das_j = D·σ − (N−1)·travel_j` with `travel ∝ Δt` means `d_front,j = D·σ` **identically, at every station and
every timestep** — `Das` absorbs the dt dependence by construction. The §5 integer schedule (N = 4/7/13) is
then unnecessary: `D` and `N` are specified, and the pin maintains itself under refinement.

At **N = 1** it is exact and trivial: `d = Das = D·σ`, with no travel term at all, so the handoff distance
becomes a **pure geometric parameter — span-uniform and dt-independent**. That is precisely the invariant Ryan
asked for ("the front of the transition wake row at a constant distance from the trailing edge"), achieved
directly rather than by scheduling `N` against `Δt`.

Note the direction of the feasibility constraint under refinement: at fixed `D`, finer `Δt` makes
`Das_j = D·σ − (N−1)·travel_j` *grow*, so the dt ladder becomes easier, not harder, as it refines. The binding
case is always the coarsest timestep.

## 6. Actions

1. **Production config**: **σ = 0.02R, N = 1, `Das` = 3.4σ span-uniform** (§4b) — uniform clearance at the C1
   bound, `Das`/c in band at both ends, sheet under one chord everywhere, and a dt-independent handoff by
   construction. Gated on the σ-collapse mitigation (S0d), since 0.02R is the rung that diverged twice.
2. **N=1 vs N=2 at matched uniform `d/σ`** — the A/B that settles whether the rigid row-1 representation at
   N=1 costs anything once clearance is controlled. This is the experiment that decides §4b's recommendation.
3. **dt ladder**: if the η/κ `Das` laws are retained, pin via `NWAKEROWS` per rung (N = 4/7/13, launcher-only,
   `run_dji9443_hover_ct_hpc.slurm.sh:387-388`) and rerun `p018_nt72` to see whether its +0.52% survives. If
   the §4b construction is adopted instead, the pin is automatic and this reduces to re-running the ladder.
4. **Open**: quantify the upper (sheet-fidelity) bound on `d/c` for *convecting* rows — it is what decides
   between N=1 and N=2 above, and until it exists the ranking rests on the lower bound plus minimum sheet
   length, not a verified two-sided criterion.
5. Retro-label completed cases with their `d/σ` profiles so past inboard findings can be re-attributed.

Data: `dsigma_b0.csv`, `chords.csv` (regenerate with `RHPC_HANDOFF_PROFILE_ONLY=true`; κ=1.0 in chord mode
makes the `das` column equal the local chord, which is how the chord table above was obtained).

## 7. Log

### 2026-08-05 (~16:30 MDT) — Ryan decision: §4b pair APPROVED at σ=0.03R; fixed-κ discriminator ON HOLD

Ryan redirected the next submissions from the fixed-κ clearance discriminator to the §4b uniform-d_front
construction directly ("fixed trailing-edge-to-front-of-the-transition-row distance"). Carrier: NT=36 B0
kinematics, **σ = 0.03R** (stability-safe: between L1's σ/c 0.297, which ran clean, and L2's 0.151, which
diverged twice; σ=0.02R stays gated on the S0d σ-collapse mitigation). σ recipe at NT=36: OVERLAP 2.4 /
PPS 14 ⇒ σ/R 0.0299, σ/c(0.75R) 0.234; `MERGE_R_FACTOR` = 0.138·(σ/R) ≈ 0.00413. Both arms inviscid, cold,
`P018_SETTLE_REVS=12` (20 revs, matched 10–19 windows per κ-ladder precedent; `_s2` if marginal), 64G/48 h.

| tag | N | D | `Das` (span-uniform d_front = D·σ) | `Das`/c_local range (measured) | tip sheet d/c |
| --- | --- | --- | --- | --- | --- |
| `p018_ufront_n1` | 1 | 3.40 | 0.1017R uniform | 0.384 – 1.407 | 1.41 |
| `p018_ufront_n2` | 2 | 6.50 | `6.5σ − travel_j` (0.1771R root → 0.0203R tip) | 0.260 – 1.054 | 2.69 |

Measured by the driver's offline profile (2026-08-05, `RHPC_HANDOFF_PROFILE_ONLY=true`, mesh 45_185_ct4,
OVERLAP 2.4 / PPS 14 ⇒ σ/R 0.02992): **d/σ exactly uniform — min = max = median = D on both sheddings, both
arms**; both `Das`/c ranges inside the 0.25–1.5 band (the min sits at the blade's max chord, outboard of the
root — the driver's `station_chords` numbers supersede the §4b table's root/tip-endpoint estimates).
Implementation: `DAS_UNIFORM_DSIGMA` env in `examples/rotor_hover_pressure_comparison.jl` (errors on
infeasible `Das_j ≤ 0`, mutually exclusive with `DAS_CHORD_FRACTION`/`DAS_REFRESH`, metadata field
`das_uniform_dsigma`), launcher arms `p018_ufront_n1/n2` (unconditional exports) + `das_uniform:` banner
field, unit testset "uniform-d_front Das law identity" in `test/runtests_unit_simulate.jl` (suite green).

D=6.5 is N=2's minimal in-band value (feasibility floor D ≥ 6.43 from `0.25·c_tip`); a matched-D pair is
infeasible, so the A/B compares two *admissible* uniform-clearance configs — both at or above the C1 bound
(kernel deficit ≤ 0.5% everywhere in both), so residual differences attribute to the **sheet representation**
(N=1 fully rigid at 1.41 tip-chords vs N=2 with a free row and 2.69 tip-chords), which is exactly §6 action 2.

**Pre-registered readouts (written before any data):**

1. Clean A/B = `n1` vs `n2` (same σ, same uniform-clearance law): M1+M2 on matched windows. PASS ⇒ the
   rigid-row-1 representation costs nothing once clearance is controlled ⇒ §4b's N=1 recommendation stands
   (cheapest, sheet < 1.5c everywhere). FAIL ⇒ rigidity (or sheet length) is a live error term; N=2-uniform
   becomes the recommended law and the d/c upper bound (§6 action 4) becomes the binding open item.
2. Falsifiable clearance prediction: vs B0's Γ̄(r/R), the **inboard-localized deficit family (N=1 −2.8% at
   r/R 0.31, small-Das, σ-ladder lobes) should be ABSENT** in both arms (uniform d/σ ≥ 3.4 at every station).
   If an inboard deficit persists at uniform clearance, clearance is NOT its cause and §2's reinterpretation
   falls.
3. Caveat pre-committed: vs-B0 comparisons confound the σ change (0.68c → 0.234c) with the Das-law change —
   they are context, not verdicts. The verdict-bearing comparison is (1); the σ axis stays owned by Phase 4/12.

Submission gated on the mandatory src md5-sync (this session's uncommitted `WakeHealthMonitor` + driver
changes). Implementation is driver-level per §4b (station array → `set_Das_station_lengths`).

### 2026-08-05 (~18:30 MDT) — pair SUBMITTED

Gates before submission, all green: unit suite (incl. new "uniform-d_front Das
law identity" testset), NT=6 40_40 coarse smoke (N=2 arm, runs clean, finite
CT, known all-NaN `CT_kj` column only), offline exact-uniform d/σ profile both
arms, deploy of the 4 changed files md5-verified against the cluster
(simulate_monitors [WakeHealthMonitor], simulate [wall clock/flush], driver,
launcher — only those 4 differed; the rest of `src/` matched).

- `p018_ufront_n1` = job **13057253** (48 h / 64G / SETTLE=12)
- `p018_ufront_n2` = job **13057254** (48 h / 64G / SETTLE=12)

64G per the memory ruling (128G is authorized only for σ-L2 and finer;
σ/c 0.234 is coarser — expected ~224k final particles, between L1's 176k@64G-ok
and L2's 347k@64G-OOM; the WakeHealthMonitor tripwire is active in these runs,
the first campaign runs to carry it). Banner verification pending (mandatory);
cost estimate ~38 h for 20 revs. Score on matched 10–19 windows; note the
formulation-independence of the κ response (phase_02 §2026-08-05 (d)) now
makes this pair the primary Das-axis discriminator.
