# Phase 2d Log — Tip-Cap Circulation Discrepancy Investigation

Plan: [Phase 2d — Tip-Cap Circulation Discrepancy](../../plans/dji_convergence_20260722/phase_02d_tip_cap_discrepancy.md)

Dashboard: [DJI convergence progress](../../plans/dji_convergence_20260722/dji_convergence_20260722.md)

## Current snapshot

Status: **Phase 2d COMPLETE and APPROVED by Ryan (2026-07-25). Appendix G added
2026-07-25** (Ryan's follow-up question on flow tangency; local-only, no HPC). Full
evidence and verification:
`data/dji_convergence_20260722/phase_02d_tip_cap_discrepancy/phase_02d_report.md`.

**Appendix G summary.** A Morino/Dirichlet solve never imposes tangency (σ = −U_ext·n; the
solver enforces the interior potential), so it was measured for the first time here. Two
enabling facts: after `steady!` `body.velocity` is already the exterior surface limit
(`src/FLOWPanel_simulate.jl:522-538`), and `kinematic_velocity!` does `Vcp -= dv`
(`src/FLOWPanel_frames.jl:350`) so the field is blade-relative and `U·n = 0` is exactly
right. Residuals were recomputed from the archived `_panels.csv` strengths **without
re-solving** — O(n) memory, which put the whole six-rung ladder (incl. n=201/249, whose
dense G is 12.8/19.6 GB) back on the laptop. Ordering gate 0.0 error on every case;
faithfulness gate vs a real solve at n=97 = 2.4e-6 (FMM-vs-direct truncation).

**The decisive control is Neumann.** On matched tip-capped root-open meshes, through the
identical evaluation path, the Neumann solution satisfies U·n = 0 to rms **4.8e-7·ΩR**
(tip-cap slivers included, 2.1e-6) while Dirichlet at the same rung sits at rms
**3.11e-1·ΩR** — five orders worse. That proves (a) the evaluation machinery, including
near-singular sliver pairs and the half-jump, is correct, and (b) the Dirichlet residual is
a genuine property of the Morino solution, not of post-processing or of small panels.

Findings: bulk tangency ~1–3% of local blade speed at n=145; residual concentrated at the
**tip-cap apex slivers** (top 1% of panels carry 93% of Σ(U·n)²; worst at r/R=0.998,
charlen 1.4e-4·R, |U| = 991 m/s). Mid-span (r/R 0.3–0.7) is best; TE panels are NOT a hot
spot. ∇μ reconstruction **exonerated** (tri/tri_robust/quad identical to 4 digits;
half-jump rms ~4e-7, six orders down). Under the study's chordwise-only refinement the
residual **grows** (inboard rms 2.49e-2→6.84e-2, ×2.75); the dominant cause is **panel
aspect ratio** — median inboard AR 7.0→17.8, residual scales ~linearly with AR, and within
the AR 8–16 band growth is only +29% vs +175% for the aggregate. Dir–Neu velocity gap
monotone-decreasing at all 7 stations on all 4 rungs, converging at least as fast as the
∫Γ gap (−41% at r/R=0.5 vs −19% for ∫Γ over 97→185); Dirichlet self-convergence ~first
order. Thrust contamination (exactly ΔP = ½ρ(U·n)², reference-free) 0.84% of this solve's
CT at n=145 rising to 3.24% at n=249, but small in absolute terms (ΔCT ~1e-4) because the
offending panels hold 0.009% of wetted area. **Production prescription unaffected**
(Dirichlet G is built from the unregularized potential kernel).

**A regularization "fix" was tested and REJECTED** (recorded in G2b because it is a trap).
Raising `kerneloffset_panel` to R·1e-3 cuts the Dirichlet max residual 80× and even removes
the divergence — but the same change applied to the Neumann control **destroys** its
tangency, rms 4.97e-7 → 4.88e-2. The offset consistent with the operator that produced the
solution is the only faithful one; the Dirichlet "improvement" is coincidental
cancellation, not accuracy. No kernel-offset change is recommended.

Actionable recommendation: **control panel aspect ratio — refine spanwise alongside
chordwise.** Also, for Phases 2e/5, prefer velocity-based convergence checks over Γ-based
ones: a converged ∫Γ does not certify a converged velocity field, and the tangency residual
and canonical-probe ΔU are cheap (no re-solve).

**Root cause: the round tip cap with OpenVSP CapUMinTess=3** (the 20260723 `capped`
recipe). It alone produces the fragile ±10% non-monotonic capped/Dirichlet tip
circulation ({185,201} overshoot artifact). Eliminated with direct evidence: extraction,
conditioning (H3), regularization (H5 — architecturally impossible: the Dirichlet G is
built from the UNREGULARIZED potential kernel), wake attachment (anomaly survives every
Das length and semi-infinite wake; the matched-wake design is gap-neutral), first-order
cap geometry, and spanwise sampling. Round-ct2, round-ct4, and flat caps agree ≤0.5% at
every rung — ct3 is the outlier. The earlier "{145,249} trustworthy" grouping was
coincidence inside ct3's scatter.

**Fix (outcome (i), Ryan prefers round): round tip cap CapUMinTess=4** —
`scripts/generate_dji9443_mesh.sh <span> <chord> flat round 4`. **F verified on a
six-rung ladder (97→249):** tip Γ@r/R=0.942 = 0.1282→0.1299 strictly monotone, climbing
toward Neumann from below, no overshoot; all integral gap bands monotone-decreasing; Dir
∫Γ changes ≤0.22%/rung from n=121 (numerically converged per the study's 1% criterion).
**Matched-tip Neumann references (Ryan's idea: root-open + tip-capped, `none round 4`)**
close the geometry share of the tip gap (Neu tip Γ 0.143→0.138); remaining Dir–Neu gap
(2.1% full-∫Γ at n=249, monotone) is the slow Phase 2b chordwise mechanism + the
still-unmatched root cap. **Production prescription: flat root + round-ct4 tip, 45 span,
n_airfoil ≥ 145 (185 preferred); do NOT use the 20260723 round-ct3 `capped` meshes for
Dirichlet.**

D (simple-case oracle) not needed — cause was isolated causally on the DJI meshes (E2,
approved by Ryan). No `src/` modifications were made (unit tests not required).

Goal: isolate with high confidence the cause of the anomalous, non-monotonic **tip**
(r/R ≳ 0.9) bound circulation in the Phase 2c capped/Dirichlet DJI solves, and decide
whether it is a fixable mesh/geometry issue, a fixed-regularization artifact, or an inherent
Morino/numerical tip sensitivity.

**End goal (Ryan, 2026-07-24): Dirichlet is the production formulation** for the DJI rotor
(Neumann is the reference oracle only), with a Green wake-potential reconstruction run
planned later (Phase 4) to close much of the experimental gap. So the deliverable is the
**mesh/settings prescription under which Dirichlet tip circulation is trustworthy** — not
just attribution, and "use interior stations only" is a last resort.

**Scope (re-scoped after review, 2026-07-24): A + C are the required core; B is cheap —
keep it; D and E are conditional** (run only if A+C leave the cause ambiguous, and reassess
with Ryan before starting them). New capability: `scripts/generate_dji9443_mesh.sh` lets
the agent generate DJI mesh variants (independent root/tip cap modes, `--cap-tess`) locally
— see the plan's tool section; this also documents the cap recipe (capped = flat root +
round tip, CapUMinTess=3).

**Success criterion (Ryan, 2026-07-24): fix-and-verify, not diagnose-and-report.** If a
study identifies an issue, the phase must also implement the solution — **targeted code
changes are in scope** (run the unit tests named in the plan's Verification), but **a
major code overhaul requires Ryan's permission before pursuing it**. The phase exits only via
**step F**: fix applied, ladder re-run (new meshes via the script as needed), Dirichlet tip
physical (approaching Neumann from below) and monotonic, and the Dir–Neu gap ≤1% (or
demonstrably converging with the required resolution/settings prescribed).

### Observation carried from Phase 2c

Capped/Dirichlet DJI tip Γ_TE is non-monotonic in n_airfoil and partly non-physical:
145→0.138 (below Neumann, expected), 185→0.153 & 201→0.154 (**overshoot** Neumann ≈0.143),
249→0.139 (below). Uncapped/Neumann tip is identical across all four meshes. The tip caps
were generated by an **identical recipe** for all meshes (per Ryan). This tip artifact made
the Phase 2c outboard gap falsely look "converged" (~0.28%) at 185/201; the honest,
interior-governed gap is ~2.2% and still decreasing.

### Hypotheses already ELIMINATED (read-only, during Phase 2d planning — do not redo)

grad_mu basis (bit-identical Γ :quad vs tri_robust); TE control-point misalignment (paired
upper/lower TE centroids chordwise-aligned to ~4e-6 m, same in 145c/185c); TE pairing /
single-sided tip (all 44 paired); watertightness (all watertight); cap normal winding flip
(all outward); degenerate/sliver caps (none; aspect ≤5); two cap recipes (identical recipe).
⇒ The tip TE **extraction** geometry is essentially identical between a low (145) and high
(185) mesh, so the ~10% tip-Γ difference is a genuine solved-field difference near the tip,
capped/Dirichlet-only, invisible to first-order geometry metrics.

### Open hypotheses (see plan for tests)

- H1 — cap↔surface seam / sharp-tip-TE-corner interaction.
- H2 — which quantity carries the overshoot (μ_upper / μ_lower / cap doublet).
- H3 — inherent Morino / (I−B) tip sensitivity or conditioning.
- H4 — which group is correct (185/201 overshoot ⇒ probably 145/249 are trustworthy).
- H5 — **fixed regularization radius vs shrinking / high-aspect-ratio tip panels**
  (Ryan-added). `kerneloffset_panel = R*1e-10`, `kernelcutoff = R*1e-13` are fixed, NOT
  panel-relative; at the thin tip, upper/lower panels are close and possibly high-AR, so a
  fixed `core_size` could act non-monotonically. Test via a kerneloffset/kernelcutoff sweep
  on 185c.

### Handoff for the executing agent

- Run **A → B → C**, then **stop and reassess before D/E** (they are conditional on A+C
  being inconclusive). Phase A (instrumented solve) is the decisive first step: it confirms
  field-vs-extraction (`circulation_slice` vs `circulation_te`) and decomposes
  μ_upper/μ_lower/cap. If C confirms H5, verify a panel-relative regularization restores
  tip monotonicity across the ladder — that setting becomes the production Dirichlet
  configuration.
- **The phase does not end at diagnosis.** Whatever the cause, implement the fix (mesh
  regeneration via the script, targeted `src/` code change + unit tests, or settings) and
  complete **step F**: re-run the ladder and verify a converged Dirichlet solution per the
  plan's F criteria, before writing the decision report. If the fix would require a major
  code overhaul, stop and ask Ryan's permission before pursuing it.
- Control meshes (if E is reached, or a quick variant helps earlier): generate locally with
  `scripts/generate_dji9443_mesh.sh` (usage in its header; e.g. `45 185 flat flat` for a
  flat-tip-cap control), then `scp` to the cluster for large solves.
- Compute note: 145c fits locally (~6 GB); **185c (~10 GB dense) is borderline on the 17 GB
  laptop — prefer HPC** (BYU `orc`; see `agent_policies/HPC.md`, the Phase 2c HPC handoff in
  its log, and `examples/run_dji9443_mesh_convergence_hpc.slurm.sh`). Cluster repo:
  `/home/rander39/projects/FLOWPanel.jl`, branch `fastmultipole`; keep cluster src in sync
  with local (a stale-src monitor bug in Phase 2c produced degenerate `r_over_R` — verify
  `src/FLOWPanel_simulate_monitors.jl` and `src/FLOWPanel_liftingbody.jl` md5 match local).
- Stop at the end of Phase 2d for Ryan's approval before Phase 3.

## Working records

| Date | Command or file | Purpose/result |
|---|---|---|
| 2026-07-24 | `examples/dji9443_mesh_convergence.jl` (+modes `tipdiag`, `tipgeom`; PHASE2D_* env) | Phase 2d instrumentation added (additive; default modes unchanged). `tipdiag`: full strength dump w/ region tags, μ_upper/μ_lower per TE section, loop-integral ∮V·dl circulation at offset stations (the monitor's crossing-sum slice is topologically zero on watertight bodies — replaced), 1-norm RCOND via LAPACK gecon, optional exact dense residual (PHASE2D_KEEP_A=1). `tipgeom`: experiment B geometry diff. |
| 2026-07-24 | `PHASE2C_MODE=tipgeom` → `data/.../phase_02d_tip_cap_discrepancy/geom/` | **B COMPLETE: no geometric invariant splits {145,249} from {185,201}.** Corner cap-fan counts split {145,185}=2/2 vs {201,249}=1/1 (wrong grouping); cap AR (med 1.6→2.0), seam valence (med 3), nn distances all vary smoothly with n_airfoil. Tip surface AR is high and growing (med 13→22, max 62→107). Weakens H1; strengthens H3/H5. |
| 2026-07-24 | `PHASE2C_MODE=tipdiag PHASE2C_CASE=dji45_145c` (local) | A(145c) instrumented solve: rcond=1.12e-5 (cond₁≈9e4, benign), Γ_TE matches Phase 2c exactly; slice confirmed identically zero (watertight topology, not station alignment). Loop-integral rerun in progress. |
| 2026-07-24 | `examples/run_dji9443_tipdiag_hpc.slurm.sh` | Phase 2d HPC launcher drafted: A-base runs (185c w/ exact residual, 185u, 201c, 249c) + C sweep on 185c (koff_panel 0/1e-7/1e-5/1e-4, kcutoff 0/1e-7). |
| 2026-07-24 | Slurm job 12893604 (COMPLETED, 32 min) → raw/ CSVs retrieved | **A+C round 1 results.** (1) **H3 ELIMINATED**: rcond smooth & benign across ladder (1.12e-5→3.15e-6, no 185/201 spike); 185c dense residual 5.7e-15. (2) **H5 ELIMINATED — architecturally**: entire koff/kcutoff sweep bit-identical incl. koff=0 and koff=R·1e-3; local test at R·1e-2 also bit-identical; root cause: `compute_source_dipole` applies `reg_term` only to velocity/gradient — the scalar potential (−μ·tan_term) is UNREGULARIZED, so kerneloffset can never touch the Dirichlet (potential-based) G. (3) **H2**: 185 tip overshoot carried SYMMETRICALLY by μ_upper AND μ_lower (both grow ~9–12% at r/R 0.93–0.97); 145 vs 249 agree to ~0.5%. (4) **H4 confirmed**: fixed-point near-field probes — 145c vs 249c tip fields agree to 0.24%, 185c deviates 4.4% at tip while 0.08% mid-span; {145,249} = trustworthy group, anomaly is a genuine localized solved-field difference. (5) μ deviation is TE-CONCENTRATED (mid-chord ~2–3%, TE ~12%) → wake-attachment (Kutta) machinery at the tip implicated. |
| 2026-07-24 | `wakestrip_check.jl` (scratchpad) | Tip CP↔wake-strip distances shrink smoothly with n_airfoil (1.3e-5→7.6e-6 m); no {185,201} grouping from proximity alone. Note: the doublet-potential kernel JUMPS across the finite attached strip (~4 mm at tip, unregularized tan_term) → side-coincidence sensitivity plausible; testing causally via Das/semiinf perturbations. |
| 2026-07-24 | Slurm job 12893855 (round 2) + local 145c controls | Causal wake-attachment test: DAS_ETA ∈ {0.1,0.4,1.0} on 185c, semiinf variants on 145c/185c/201c/249c/185u, eta controls on 145c/249c. New tipdiag env knobs PHASE2D_DAS_ETA / PHASE2D_DAS_MIN / PHASE2D_SEMIINF (semiinf normalizes Das to unit vectors — semi-infinite kernels require unitary directions). |
| 2026-07-25 | round-2 results (job 12893855 COMPLETED, 0 failures) | **Wake-attachment mechanism ELIMINATED.** Tip Γ level is strongly Das-length-dependent everywhere (sec3 145c: 0.113@η0.1 → 0.142@0.2 → 0.168@0.4 → 0.224 semiinf), but (a) Neumann moves in LOCKSTEP (Dir/Neu tip ratio 0.985–0.987 at every η and semiinf at 145 — the matched-wake design is gap-neutral), and (b) the {185,201} anomaly PERSISTS at every wake setting incl. semi-infinite (+10–11% over {145,249}, which agree to 0.7%; Neumann 185u≡145u to 0.2%). The defect is intrinsic to the capped/Dirichlet discretization of the 185/201 meshes. |
| 2026-07-25 | `plane_proximity.jl` (scratchpad) + cap-area check | tan_term plane-jump hazard ruled out: min CP-to-other-panel-plane distance (in-footprint) shrinks smoothly 4.3e-6→2.4e-6 m, none <1e-6, no grouping. Total tip-cap areas identical to 0.02% (8.772–8.774e-6 m²) — cap SURFACE identical; only its tessellation differs. |
| 2026-07-25 | E2 approved by Ryan → generated `dji9443_20260725_45_{145,185}_{flatcaps,capped_captess2,capped_captess4}.msh`; jobs 12893980 (failed: TE walk hit 0.1R bbox before the controls' deeper end node — driver `trace_shedding` now passes end_node=nothing when the detected inner node is <0.1R) → 12893982 (COMPLETED) + local 145 controls | **E2 CAP-TREATMENT MATRIX (tip Γ @r/R≈0.945): round-ct3 is the fragile outlier; flat/ct2/ct4 form a stable family.** ct3: 145→0.1421, 185→0.1578, 249→0.1427 (±10% scatter). flat: 145→0.1292, 185→0.1296; ct2: 0.1289/0.1293; ct4(185): 0.1294 — three treatments × two n_airfoil within 0.5%, ≈10% below open-tip Neumann (0.1443/0.1444). The earlier "{145,249} trustworthy" grouping was coincidence within ct3 fragility. **Root cause named: round CapUMinTess=3 tip cap.** |
| 2026-07-25 | interpolated full-span flat-vs-base comparison | **CONFOUND: cap settings change VSP's shared U-tessellation** — flat/ct2 meshes have 40 blade stations vs 44 (study), and differ from base across the WHOLE span (+11% @r/R=0.3, −4% @0.6), though flat-145≡flat-185 to 0.1% everywhere. Inboard flat sits ABOVE Neumann (+7% @0.3). Station-count/clustering confound must be removed before prescribing the flat family: solving a 49-span flat-145 (tip spacing ≈ study's) to separate spanwise-sampling from cap-recipe effects. |
| 2026-07-25 | 49-span flat-145 solve + validate_case on all controls | **Spanwise sampling exonerated**: flat-49 (44 stations) ≡ flat-45 (40 stations) to 0.1% at every station. All controls pass full topology validation (paired 40/40 or 44/44, monotone radii, symmetric blades). Root cap counts: flat kept root identical (1136) yet shifted inboard like ct2 (which changed root to 568) — root cap exonerated. Local Γ_TE differences between families partially integrate out (flat's inboard wave ±7% → only ~+0.9% inboard integral). |
| 2026-07-25 | flat-121 solve + integral-gap table (all cases, matched bins [0.125,0.975]) | **F DIRECTION SETTLED — stable-cap ladder satisfies the study criteria.** Flat 121→145→185: tip Γ monotone climbing, gaps monotone decreasing. Cap-treatment-insensitive (flat≡ct2≡ct4 ≤0.5%). Round-ct3 recipe confirmed fragile. Extending flat ladder: 97 local, 201/249 HPC job 12894030. |
| 2026-07-25 | Ryan: prefers ROUND caps; suggested matched tip-capped Neumann refs | Round-ct2/ct4 are as stable as flat (only ct3 is defective) → production family switched to **round CapUMinTess=4**. Generated ct4 ladder 97–249 + root-open/tip-capped Neumann refs (`none round 4`, `none flat`); solved 97/121/145 + Neu-145 locally, 201/249 + Neu-185 on HPC job 12894044 (flat 201/249 completed on 12894030 as corroboration). Tip-capping the Neumann body lowers its tip Γ 0.143→0.1382–0.1386 (cap-shape-insensitive), isolating geometry vs formulation in the tip gap. |
| 2026-07-25 | `phase_02d_report.md` + final six-rung ct4 ladder | **F VERIFIED.** ct4 tip Γ@0.942: 0.1282/0.1286/0.1290/0.1294/0.1295/0.1299 (97→249) vs matched Neu 0.138 — strictly monotone from below, no overshoot. Integral gaps vs matched Neu all monotone-decreasing (full 2.90→2.11%, inboard 2.03→1.43%, outboard 4.62→3.46%); Dir ∫Γ ≤0.22%/rung from n=121. Decision report written; phase closed pending Ryan's approval. |
| 2026-07-24 | `plans/.../phase_02d_tip_cap_discrepancy.md`, this log | Phase 2d staged: context, eliminated hypotheses (from read-only pre-work), open hypotheses incl. H5 (regularization), experiments A–E, decision rule. No experiments run. |
| 2026-07-24 | plan + log + dashboard | Phase 2d reviewed & re-scoped: A+C core, B kept, D/E conditional; decision rule reweighted (Dirichlet = production path); `scripts/generate_dji9443_mesh.sh` documented as an agent tool. |
| 2026-07-24 | plan + log + dashboard | Success criterion strengthened (Ryan): fix-and-verify, not diagnose-and-report. New required step F — implement the fix (targeted code changes in scope; a major code overhaul needs Ryan's permission first) and re-run the ladder to a demonstrated converged Dirichlet solution + settings. |
| 2026-07-25 | `examples/dji9443_mesh_convergence.jl` (+modes `tangency`, `gprobe_points`, `gprobe`) | Appendix G instrumentation (additive; no `src/` changes). `tangency`: rebuild body, load solved strengths from `<label>_panels.csv` behind an ordering gate, replay the `_steady_aerodynamics!` velocity assembly with NO solve (O(n) memory), dump per-panel residual + decomposition + per-route variants + full-surface VTK. `gprobe*`: canonical mesh-independent probe points (Phase 2d's `run_fieldprobe` rings are bbox-derived and differ ~9e-6 m between formulations, contaminating any Dir–Neu ΔU). |
| 2026-07-25 | `PHASE2C_MODE=tangency` × ct4 ladder {97,121,145,185,201,249} | G1. Ordering gate **0.0** on every rung; linearity ≤5.3e-11. |U·n|/ΩR rms 2.39e-1→5.00e-1, max 14.2→44.4 — **growing** with refinement. Concentrated: top 1% of panels carry 93% of Σ(U·n)²; worst are tip-cap apex slivers at r/R=0.998 (charlen 1.4e-4·R, |U| 991 m/s). Mid-span (r/R 0.3–0.7) rms ~1.3e-2. TE panels NOT a hot spot. Median |U| = 0.50·ΩR — bulk field correct. |
| 2026-07-25 | scratchpad `gate_faithful.jl` (n=97, real `steady!` + dense G) | **G0.4 faithfulness gate PASSED**: reconstruction vs genuine solve = 2.4e-6 rel (FMM-vs-direct truncation), strengths bit-identical, residual stats agree to 5 digits. The large residual is real, not a reconstruction artifact. |
| 2026-07-25 | scratchpad `tanalyze.jl`, `tanalyze2.jl`, `axischeck.jl` | G1/G2/G4 tables. ∇μ **exonerated** (tri/tri_robust/quad identical to 4 sig figs; half-jump rms ~4e-7). AR stratification: within-bin residual nearly flat, but median inboard AR 7.0→17.8 across the ladder ⇒ chordwise-only refinement inflates AR and AR drives the residual. Thrust axis verified (F_y=F_z=0 exactly). |
| 2026-07-25 | scratchpad `koffsweep.jl` (n=145, strengths held fixed) | Velocity-side regularization sweep koff/R ∈ {1e-10,1e-6,1e-4,1e-3,1e-2,3e-2}: apparent optimum at R·1e-3, max residual 21.9→0.275 (80×), rms 0.311→0.024. Phase 2d item 3 eliminated regularization for the SOLVE (potential kernel unregularized) — but velocity is exactly where `reg_term` acts. **Later REJECTED — see the Neumann control row below.** |
| 2026-07-25 | `koffsweep.jl` × full ladder at koff=R·1e-3 | At the raised offset the divergence disappears (rms flat 2.377e-2→2.444e-2; inboard rms decreasing −7.4% vs +175% at default). Looked like a fix; **it is not** — see next rows. |
| 2026-07-25 | `PHASE2C_MODE=tangency` on `dji45_{97,145}u_tipround4` | **Neumann control — the decisive result.** Same evaluation path: rms **4.78e-7 / 4.97e-7**, max ~1e-5, tip-cap rms 2.1e-6. The machinery is exact on an exactly-tangent solution ⇒ the Dirichlet residual (rms 3.11e-1 at the matched rung) is genuine, and the default offset is the faithful one. |
| 2026-07-25 | `koffsweep.jl` on `dji45_145u_tipround4` (robustness check of the proposed fix) | **Fix REFUTED.** Raising koff on the Neumann control destroys tangency: rms 4.965e-7 (koff=R·1e-10, consistent with the solve) → 7.56e-2 (1e-4) → 4.88e-2 (1e-3) → 1.06e-1 (1e-2). Moving off the solve-consistent offset biases the velocity field; the Dirichlet "improvement" was coincidental cancellation. Recommendation withdrawn before delivery; retained in the report as G2b so it is not re-proposed. |
| 2026-07-25 | `scripts/generate_dji9443_mesh.sh 45 {97,121} none round 4`; `PHASE2C_MODE=tipdiag` on both | G3 ladder extension: the two missing Neumann partners generated and solved locally (rcond 4.6e-7 / 3.7e-7), giving a 4-rung Dir–Neu velocity ladder instead of a 2-point difference. |
| 2026-07-25 | `PHASE2C_MODE=gprobe_points` + `gprobe` × 10 cases | G3. Dir–Neu |ΔU|/|U| monotone-decreasing at all 7 stations on all 4 rungs (r/R=0.5: 1.63%→0.97%; tip: 13.4%→10.2%), converging ≥ as fast as the ∫Γ gap (−41%/−24% vs −19% over 97→185). Dirichlet self-convergence ~first order in 1/n. |
| 2026-07-25 | scratchpad `dct_koff.jl` (n=145, 249 × koff {1e-10, 1e-3}) | G4. ΔP = ½ρ(U·n)² exactly (Bernoulli reference cancels in the blade-relative frame). At the solve-consistent offset: \|ΔCT/CT\| = **0.837%** (n=145) → **3.245%** (n=249), growing with the aspect-ratio trend; ΔCT stays ~1e-4 in absolute CT because the pathological panels are area-negligible (\|U\|>5× local = 0.0091% of wetted area). The raised-offset values (0.023%/0.025%) belong to the rejected G2b configuration and are not a result. |
| 2026-07-25 | report + plan + log + dashboard | Appendix G written; production prescription confirmed unaffected; open recommendation (raise the Dirichlet reconstruction kernel offset in `src/`) left for Ryan's decision. |

## Dated entries

### 2026-07-25 — Appendix G: flow-tangency residual (Ryan follow-up, local-only)

- Ryan asked, after Phase 2d closed: how well is flow tangency satisfied on the "converged"
  mesh, where is the residual largest/smallest, why if large, and do the Dirichlet and
  Neumann **velocity fields** agree and converge as well as the circulation does?
- Motivation the appendix makes explicit: all Phase 2b–2d evidence is circulation-based,
  and Γ is a TE-panel μ difference, whereas forces come from the reconstructed surface
  velocity. A Morino solve never imposes tangency, so it was simply unmeasured.
- Method note worth keeping: because `run_tipdiag` archived per-panel strengths, the entire
  ladder could be re-analyzed with **no re-solve** at O(n) memory. That is what made
  n=201/249 (dense G = 12.8/19.6 GB) tractable on the laptop under the local-only
  constraint. Two gates (ordering, faithfulness) make the shortcut defensible.
- Course correction during the work: Phase 2d's `_fieldprobe.csv` pairs are **not** valid
  for a cross-formulation ΔU comparison (bbox-derived rings differ ~9e-6 m between capped
  and root-open meshes, enough to move U comparably to the effect being measured). Replaced
  with a canonical point set.
- Headline: Neumann satisfies tangency to ~5e-7 through the identical evaluation path,
  Dirichlet to ~3e-1 — the residual is a real property of the Morino solution, and the
  growth under refinement is driven mainly by rising panel aspect ratio from
  chordwise-only refinement. Thrust contamination 0.84%→3.24% across the ladder but only
  ~1e-4 in absolute CT. The Phase 2d production prescription stands.
- **Process note worth keeping.** A kernel-offset change looked like a clean 80× fix and
  was written into the report as a recommendation. Running the same sweep on the Neumann
  control — whose true residual is known to be ~0 — showed it *destroys* tangency by five
  orders. The recommendation was retracted before delivery. Generalizable lesson: when a
  regularization knob "improves" a diagnostic, verify it against a case whose correct
  answer is known independently; otherwise an inconsistency with the solve operator reads
  as an accuracy gain.

### 2026-07-24 — Phase 2d staged (documents only)

- Ryan flagged the Phase 2c 145→185 tip jump and asked whether grad_mu changed across it
  (it did — 81–145 used `:quad`, 185–249 `:tri,tri_robust`), but re-solving dji45_145 with
  tri_robust gave bit-identical Γ ⇒ grad_mu ruled out. Per-station inspection localized the
  jump to r/R≳0.9 and showed 185/201 Dirichlet **overshoot** the flat Neumann tip while
  145/249 sit below it — a tip-cap artifact, non-monotonic in n_airfoil, capped-only.
- Planning-time read-only mesh analysis eliminated the obvious geometry causes (see snapshot).
  Ryan confirmed the caps use an identical recipe, ruling out two-recipes and sharpening the
  puzzle. Ryan added H5 (fixed regularization radius vs shrinking/high-AR tip panels).
- Wrote the Phase 2d plan + this log + dashboard entry. Execution handed to a fresh agent.

### 2026-07-24 — Review, end-goal clarification, and re-scope

- Reviewed Phase 2d against the study objective. Verdict: necessary (the tip artifact
  already produced one false convergence claim in Phase 2c, and Dir > Neu is non-physical
  and potentially a real code defect via H5), but the full A–E ladder was over-scoped.
- Ryan clarified the end goal: **Dirichlet is the production formulation** for the DJI
  rotor; this study verifies Dirichlet converges to the right answer and determines the
  mesh resolution/settings required. A later Green wake-potential-reconstruction run
  (Phase 4) is expected to close much of the validation gap. This drops any "fall back to
  Neumann" option and reweights the decision rule toward the actionable fixes (i)/(ii),
  with (iii) as a last resort that must still name a trustworthy Dirichlet configuration.
- Re-scoped: **A + C core, B kept (cheap), D/E conditional** on A+C being ambiguous
  (reassess with Ryan first).
- Documented `scripts/generate_dji9443_mesh.sh` (dashboard Tools section + plan): agents
  can generate DJI mesh variants with independent root/tip cap modes and `--cap-tess`.
  This resolves the "confirm the cap recipe" open item (capped = flat root + **round tip**,
  CapUMinTess=3) and upgrades E2 from Ryan-generated to agent-generated controls.
- Ryan then strengthened the success criterion: if a study identifies an issue, the phase
  must also implement the solution — targeted code changes in scope (a major code overhaul
  requires his permission first). Added required **step F** (fix
  + ladder re-run to a demonstrated converged Dirichlet solution and settings) as the
  phase's exit condition; the decision rule's three outcomes are now fix paths that must
  each end implemented and verified, with (iii) allowed only if the cause is proven
  inherent and unfixable, and even then requiring a converged configuration within the
  trustworthy mesh family.
