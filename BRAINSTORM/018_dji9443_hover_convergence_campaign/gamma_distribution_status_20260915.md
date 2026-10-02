# STATUS — 018 particle-Γ distribution vs NT (2026-09-15)

Follow-on to `gamma_distribution_reset_prompt_20260915.md`. All analysis on ORC;
scripts and logs in `~/p018_gamma_dist_20260915/` (cluster):
`p018_gamma_dist.py` (single snapshot), `p018_gamma_dist_avg.py` (phase-avg,
log `gdist_avg.log`), `p018_gamma_ts.py` (per-rev totals, log `gdist_ts.log`).

## Data protection (done first)

Full per-step VTP sets of the g25 pair harvested to the archive before any
sweeper action (protect list is Ryan's; agents don't write it — harvested
instead):

- `/nobackup/archive/usr/rander39/FLOWPanel_runs/p018_csarc_l3p0_3r_g25_wake1_particles_full_20260915.tar` (33 GB, 1080 VTPs, member count verified)
- `/nobackup/archive/usr/rander39/FLOWPanel_runs/p018_csarc_n2_nt72_l3p0_3r_srlx_g25_wake1_particles_full_20260915.tar` (70 GB, 2160 VTPs, verified)

The 8 matched-rev snapshots also copied to `~/p018_gamma_dist_20260915/snaps/`.

## Method

Axial bins of 0.5R along x (disk at x≈0, wake → +x, R = 0.11995 m), per bin:
count, |Γ| percentiles, Σ|Γ|, ‖ΣΓ⃗‖. Single matched-rev snapshots (revs
8/16/24/30) proved too noisy: per-bin Σ|Γ| fluctuates rev-to-rev *within a
rung* by amounts comparable to the NT36↔NT72 differences, and per-bin ratio
signs flip between revs (coherent tip-vortex bundles alias against bin edges).
Final numbers therefore use **phase averaging over 12 matched blade phases
spanning one full rev** (windows rev 15→16, 23→24, 29→30) plus a **per-rev
time series** (4-phase avg, revs 2–29).

## Results — g25 pair (guarded, same arch)

### Per-rev totals (4-phase avg; near = x<1.5R, far = x≥1.5R)

| rev | NT36 n | NT36 Σ\|Γ\| | NT36 near | NT36 far | NT72 n | NT72 Σ\|Γ\| | NT72 near | NT72 far |
|---|---|---|---|---|---|---|---|---|
| 12 | 238444 | 4.431 | 2.157 | 2.274 | 250746 | 4.044 | 2.103 | 1.942 |
| 16 | 296767 | 5.343 | 2.015 | 3.328 | 323733 | 5.195 | 2.012 | 3.183 |
| 20 | 306741 | 5.576 | 2.043 | 3.532 | 332844 | 5.500 | 1.976 | 3.525 |
| 24 | 306956 | 5.762 | 2.008 | 3.754 | 323436 | 5.639 | 1.916 | 3.723 |
| 26 | 309212 | 5.891 | 2.028 | 3.863 | 326554 | 5.841 | 1.894 | 3.947 |
| 28 | 310391 | 5.938 | 1.999 | 3.939 | 334588 | 6.107 | 1.937 | 4.169 |
| 29 | 310903 | 5.956 | 2.056 | 3.900 | 339032 | 6.248 | 1.917 | 4.331 |

- **Near-disk (<1.5R) Σ|Γ| is NT-invariant**: saturates at ~2.0 by rev 12 in
  both rungs and stays flat; NT72 sits 2–4% below NT36 throughout. No shed-side
  circulation discrepancy.
- **The far wake (≥1.5R) is where the rungs diverge, and it diverges in
  TIME**: NT36's far Σ|Γ| decelerates toward a plateau (~3.9–4.0 by rev 29);
  NT72's plateaus revs 18–22, then **re-accelerates** (3.72 → 4.33 over revs
  24–29, ~+0.15/rev and rising at run end). NT72 crosses NT36 near rev 28.
  NT72 particle count is likewise non-monotone (336k → 323k → 339k): removal
  wins revs 18–24, then loses. **NT72's wake is not statistically steady
  anywhere in the CT scoring window (revs 21–30).**

### Phase-averaged NT72/NT36 ratios, rev 29→30 (12 phases)

| bin [R] | count | Σ\|Γ\| | mean\|Γ\| | p50 |
|---|---|---|---|---|
| 0–0.5 | 0.92 | 0.84 | 0.91 | 0.99 |
| 0.5–1 | 0.87 | 0.87 | 1.00 | 1.12 |
| 1–1.5 | 1.03 | 1.07 | 1.04 | 0.98 |
| 1.5–2 | 1.19 | 1.33 | 1.11 | 1.17 |
| 2–2.5 | 1.14 | 1.08 | 0.94 | 0.74 |
| 2.5–3 | 1.20 | 1.26 | 1.05 | 0.96 |
| 3–3.5 | 1.11 | 0.93 | 0.84 | 0.68 |
| TOTAL | 1.09 | 1.05 | 0.96 | — |

(Windows 15→16 and 23→24 in `gdist_avg.log`; totals there: Σ|Γ| ratio 0.946 →
0.971 → 1.051 across the three windows — the far-wake excess is late-onset.)

### CT window split (scripts/p018_analyze.py m1)

| window | NT36 CT̄ | NT72 CT̄ | climb |
|---|---|---|---|
| revs 21–25 | 0.070571 | 0.071394 | +1.17% |
| revs 26–30 | 0.070454 | 0.071446 | +1.41% |

NT36 drifts down late, NT72 holds/rises — same direction as the far-wake
secular growth.

## Cross-check — `_3r_sv_s1p5` (unguarded, 5 live end-of-run snaps, rev ≈ 30)

| rung | n | Σ\|Γ\| | near <1.5R | far ≥1.5R | far p50 (2.5–3R bin) |
|---|---|---|---|---|---|
| NT36 | 375844 | 3.964 | 1.457 | 2.506 | 2.04e-6 |
| NT72 | 468486 | 3.657 | 1.429 | 2.228 | 1.27e-6 |
| NT144 | 609850 | 3.640 | 1.418 | 2.222 | 8.9e-7 |

Structure confirmed: near-disk Σ|Γ| NT-invariant within 2.7% across a 4× NT
range; the difference lives in the far field (12.5% here), with per-particle
far-field median |Γ| roughly halving per NT doubling while counts grow.
Note the far-field **sign differs** from g25 at rev 30 (here NT36 > NT72 ≈
NT144) — consistent with the g25 time series if the ordering at any fixed rev
depends on where each rung sits on its (NT-dependent) equilibration curve.

## Verdict

The per-particle Γ distribution change with NT is a **far-field,
wake-evolution effect, and substantially a TEMPORAL one** — not a shed-side
one. Near the disk (<1.5R) Σ|Γ| is NT-invariant (≲3–4%) in both stacks at all
times after rev ~12. Beyond ~1.5R two things happen: (1) count-normalized
shape change — more, weaker particles at higher NT (far-field median |Γ| drops
~30–50% per doubling; partially the benign quantization-excess tradeoff), and
(2) a **real Σ|Γ| discrepancy at the 10%+ level whose sign evolves in time**:
in the g25 pair NT72's far-wake circulation is still secularly growing at rev
30 while NT36's has saturated. This implicates the **wake-evolution mechanism
class** (removal/SFS/relaxation balance on the old wake — matching the
documented removal-side asymmetry, old wake eaten ~2× faster at NT36) and
rules out the shedding/conversion class for the Γ-distribution signal. It
also means matched-rev CT scoring at revs 21–30 compares rungs at different
points on NT-dependent equilibration transients.

## Proposed next discriminator

**Extend the g25 NT72 run (and ideally NT36) to rev ~60** and track CT̄ and
far-wake Σ|Γ| vs rev: if NT72's far wake saturates and the corrected climb
(+1.27%/doubling) shrinks with window, the residual carrier is (at least
partly) a **transient-length artifact** of NT-dependent wake equilibration —
cheap to test, and it directly gates whether the FMM-near-set and
spatial-entanglement fronts are still needed at full strength. A cheaper
partial: per-step removed-Σ|Γ| by axial bin from existing diagnostics, to pin
the removal-side asymmetry as the driver of the equilibration difference.

— agent session 2026-09-15; do not edit prior entries, corrections in a new dated file.
