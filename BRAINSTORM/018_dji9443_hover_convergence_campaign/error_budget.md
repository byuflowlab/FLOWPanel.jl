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
| 4 | σ extrapolation | $|\bar{C_T}(\sigma^*) - C_{T\infty}|$ from the Phase 4 fit | | |
| 5 | Shedding distance | N=8 vs N=4 spot-check delta (Phase 5) | | |
| 6 | Truncation depth | 6R vs 4R delta (Phase 6) + prior null ≤0.22% | | |
| 7 | Formulation (Green Δ) | Phase 7 ΔCT̄, ΔΓ̄ | | |
| 8 | Merging | Phase 8 ΔCT̄ (target ≤0.25%) | | |
| 9 | Relaxation | ≈ −0.005 CT, monotone in rlxf, NOT converged (006 ladder; stock rlxf per ruling) — carried as one-sided systematic | ≈ −0.005 | fixed by ruling |
| 10 | Pressure monitor | PressureBernoulli(unsteady=false) moving-body regression flag; same-state force-method spread 16.3% (001) disclosed | disclosure | |
| 11 | Dirichlet tangency | thrust contamination ½ρ(U·n)² ≈ 0.8–1% of CT at n=185 (2d App. G) | disclosure | |
| 12 | Radial truncation | 1.5R hard-coded, not converged | disclosure | |
| 13 | Mesh (prior steady evidence) | steady refinement moved CT < cycle scatter (006); cite 2c/2d | disclosure | |
| 14 | Mesh spanwise (unsteady) | 45 → 60 [→ 80] rung delta at final settings (Phase 10), CT and Γ | | |
| 15 | Mesh chordwise (unsteady) | 145 → 185 → 249 rung deltas at final settings (Phase 11), CT and Γ; interpret n=249 with the tangency systematic (term 11) | | |

Also report (not budget terms): FMM/Direct discriminator outcome (Phase 2 side
job); min_displacement floor clamp fraction at Das\*; startup-schedule
insensitivity (staged recipe fixed across all cases — deltas are differential).
