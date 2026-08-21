# 018 Error Budget (assembled in Phase 9)

Final claim form: $\bar{C_T} = X \pm \mathrm{CI_{stat}} \pm \Sigma$(systematic
terms below), $\bar\Gamma(r/R)$ likewise. Convergence and agreement with
experiment are separate claims; $C_{T,\mathrm{exp}} = 0.072$ is compared only
after this table is complete.

| # | term | how measured | value | status |
| --- | --- | --- | --- | --- |
| 1 | Statistical (finite window) | M1 bootstrap CI on the Phase 9 run | | |
| 2 | Das model-form | plateau slope from Phase 2 ladder (wing law ~0.2%/doubling) | | |
| 3 | dt residual | NT\* → 2·NT\* delta (Phase 3), Richardson if 3 rungs | | |
| 4 | σ extrapolation | $|\bar{C_T}(\sigma^*) - C_{T\infty}|$ from the Phase 4 fit; interpretation conditioned on term 16 (Ladder A co-refines h) and cross-checked by Ladder B (σ at fixed h, phase_04 §2026-08-05) | | |
| 5 | Shedding distance / d/σ clearance | N spot-check at σ\* (Phase 12 C2 pre-registered null, supersedes the N=8-at-B0 plan) + offline kernel-deficit threshold d/σ\* (Phase 12 C1) | C1 measured 2026-08-05: d/σ\* ≈ 2.6 (median) / 3.4 (p95); slope −3.3 ≈ 017 prior; 1−g(d/σ) is a conservative closed-form bound | C1 DONE; C2 gated on σ\* |
| 6 | Truncation depth | 6R vs 4R delta (Phase 6) + prior null ≤0.22% | | |
| 7 | Formulation (Green Δ) | Phase 7 ΔCT̄, ΔΓ̄ | | |
| 8 | Merging | Phase 8 ΔCT̄ (target ≤0.25%) | | |
| 9 | Relaxation | ≈ −0.005 CT, monotone in rlxf, NOT converged (006 ladder; stock rlxf per ruling) — carried as one-sided systematic. **Ladder VERDICT 2026-08-13: the downward slope is NOT MEASURABLE at this operating point — BOTH reduced rungs IGNITED from a healthy warm state** (rlxf 0.15 job 13157881: max_u breakout + min_sr collapse at step 1734–1735 = rev 48.2, +3.2 revs after the rev-45 handoff, peak max_u 47.5k m/s, min_sigma pinned at the 9.41e-5 viscous floor, scancel'd; rlxf 0.075 job 13157882: breakout at step 1684–1698 = rev 46.8–47.2, +~2 revs, peak max_u 1.12e6 m/s, wake self-annihilated 217k→11k particles, exit 0 with garbage forces — dose-response: quartered rlxf ignites ~1.2 revs sooner than halved). Stock rlxf=0.3 is LOAD-BEARING for stability at σ=0.04R/N=1 viscous (relaxation acts as stabilizer against Γ-side ignition at the viscous floor, consistent with the σ-blow-up mechanism). How to carry the term is a RYAN DECISION: reframe as model definition; probe upward rlxf>0.3; or joint σ–rlxf probe at a more stable point. Post-mortems: ledger §2026-08-13 | ≈ −0.005 one-sided, slope unmeasurable downward; upward [0.3,0.6] CT-flat within burst noise (|Δ|≲1.7% non-monotone), Γ̄ insensitive (M2 ≤1%) | downward ladder CLOSED 2026-08-13 (both reduced rungs ignite); **upward pair harvested 2026-08-14 (phase_09): stable, stability rises with rlxf (tail min_sr 0.149/0.153/0.201), slope UNRESOLVED at 13-rev windows ⇒ model-definition reframe directly supported; carrying mode = Ryan decision** |
| 10 | Pressure monitor | PressureBernoulli(unsteady=false) moving-body regression flag; same-state force-method spread 16.3% (001) disclosed | disclosure | |
| 11 | Dirichlet tangency | thrust contamination ½ρ(U·n)² ≈ 0.8–1% of CT at n=185 (2d App. G) | disclosure | |
| 12 | Radial truncation | 1.5R hard-coded, not converged | disclosure | |
| 13 | Mesh (prior steady evidence) | steady refinement moved CT < cycle scatter (006); cite 2c/2d | disclosure | |
| 14 | Mesh spanwise (unsteady) | 45 → 60 [→ 80] rung delta at final settings (Phase 10), CT and Γ | | |
| 15 | Mesh chordwise (unsteady) | 145 → 185 → 249 rung deltas at final settings (Phase 11), CT and Γ; interpret n=249 with the tangency systematic (term 11) | | |
| 16 | Particle spacing h (fixed σ) | Phase 12A rung-pair deltas (h/σ 0.5 → 0.25 → 0.125 at B0 σ/c 0.68); also conditions the interpretation of term 4 (Phase 4 co-refines h) | | OPEN — rungs submitted 2026-08-05 |

## Γ̄(r/R) budget terms (added 2026-08-03 session c)

$\bar\Gamma(r/R)$ is half the deliverable (ruling 9) but the table above is
CT-only; these rows carry the circulation-distribution uncertainty per axis.
Metric: $\varepsilon_\Gamma$ = max and RMS of $|\Delta\bar\Gamma| /
\max_r|\bar\Gamma|$ over $0.3 \le r/R \le 0.95$ (`p018_analyze.py m2`),
successive rungs on matched windows, PLUS the error *signature* (which lobe),
because CT̄ is structurally near-blind to redistribution modes (σ axis:
ε_Γ 8.78% at ΔCT̄ +1.34%, band integral +0.117%).

| # | term | how measured | value so far | status |
| --- | --- | --- | --- | --- |
| Γ1 | Statistical | per-section cycle-mean scatter over the M1 window (Phase 9 run) | | |
| Γ2 | Das model-form | successive-rung ε_Γ + signature (small Das: inboard deficit −8.9% @ r/R 0.31; large Das: outboard deficit −6.6% @ 0.84) | 8.96 / 3.11 / 5.32 % max across the 4-rung ladder — no passing pair | OPEN; re-test top rungs at σ\* |
| Γ3 | dt residual | ε_Γ NT→2NT (broad uniform shift — only axis where CT̄ and Γ̄ agree) | 0.979% max / 0.772% RMS (36→72 @ Das 0.205c) | OPEN (nt72_eta2/nt144 pending) |
| Γ4 | σ | ε_Γ per rung + 3-lobe amplitude fit vs σ (falsifiable L2 prediction: same dip at r/R≈0.76 deepens) | 8.776% max / 4.692% RMS (L0→L1) — **binding term of the campaign** | OPEN; blocked on viscous L2 (see ledger 2026-08-03c) |
| Γ5 | Shedding distance N | ε_Γ of N-rung vs N=4 (inboard-localized; collapses with Γ2 onto d/σ per the matrix test) | N=1: 2.77% (FAIL); N=2: 0.35% (PASS) | Das×N matrix in flight (13035910–13) |
| Γ6 | Viscous model | Γ̄ delta of the B0 CoreSpreading A/B (`p018_b0_visc`) — erratum: viscous silently inert in all runs to date | | OPEN (disclose + measure, Ryan 2026-08-03) |
| Γ7 | Merging / truncation / formulation | ε_Γ from Phases 6–8 nulls (merging target ≤0.5%) | | |
| Γ8 | Mesh (span/chord) | ε_Γ across Phase 10/11 rungs on a common r/R grid | | |
| Γ9 | Particle spacing h (fixed σ) | ε_Γ of Phase 12A rung pairs at fixed σ/c 0.68; discriminator: if the L1 3-lobe shape appears at fixed σ, Γ4's attribution to σ re-attributes to h/overlap | | OPEN — rungs submitted 2026-08-05 |

Also report (not budget terms): FMM/Direct discriminator outcome (Phase 2 side
job); min_displacement floor clamp fraction at Das\*; startup-schedule
insensitivity (staged recipe fixed across all cases — deltas are differential).
CT term 4 (σ) and Γ4 must quote the σ-instability root cause + viscous erratum
(ledger 2026-08-03 session c): all pre-fix runs are functionally inviscid-wake;
deltas between them remain internally consistent.
