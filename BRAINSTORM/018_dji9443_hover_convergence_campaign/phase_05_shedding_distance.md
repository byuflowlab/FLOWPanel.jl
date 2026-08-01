# Phase 5 — Shedding-Distance Verification (no ladder)

**Objective:** verify the criterion-based `NWAKEROWS` choice. Ruling 4: no
nwakerows ladder — choose N high enough that particle-on-body influence isn't
stunted by smoothing-radius overlap; if all other parameters converge, N
shouldn't matter much.

## Criterion (documented for publication)

Particle-on-body velocity error scales ~ (d/σ)^−3.4 (017 M4 static fit),
d = handoff distance = Das + (N−1)·(local TE travel/step). With N=4 at NT=36:
d ≈ 4.5c at the tip, d/σ ≥ ~3.4 even at 0.5R for the coarsest rung
(σ/c=0.68), and d/σ grows ∝ 1/σ down the ladder. The rejected legacy ladder
(non-monotone at σ≈1.5c, nrows 1→72) is cited as motivation for the
criterion-based choice, not as data.

## Case

| tag | knobs | time |
| --- | --- | --- |
| `p018_nrows8` | B0 carrier + NWAKEROWS=8 (d ≈ 9.9c tip) | 24 h (+extension if unsettled) |

Compare vs `p018_b0` (N=4) on matched windows.

## Decision

|ΔCT̄| ≤ 0.5% and ε_Γ ≤ 1% ⇒ criterion stands; delta → error-budget term 5.
This also documents 017's T3 prediction (flat/monotone-saturating behavior at
admissible σ) with a two-point check. **If it fails:** do NOT add more rows —
that non-monotone axis is the known dead end; this is a 016 trigger
(pre-authorized; see Phase 4 contingency (d)) and the finding routes to
BRAINSTORM/014's open residue.

## Exit criteria

Spot-check delta in the ledger; criterion text finalized for the report.

## Reframing (Ryan ruling 11, 2026-07-31)

nwakerows is NOT a convergence axis — "convergence of nwakerows doesn't even
make sense"; the quantities to converge are particle shedding (σ, spacing)
and Das. N is a modeling choice whose *sensitivity* is measured and
disclosed. Phase 5 therefore becomes an N-sensitivity probe centered on the
small-N side (the physically minimal sheet), motivated by Phase 2's finding
that d = N·Das confounds the Das ladder at fixed N:

| tag | knobs | job | status |
| --- | --- | --- | --- |
| `p018_nrows1` | B0 carrier + NWAKEROWS=1 (d = Das = 0.6σ) | 13011982 | submitted 2026-07-31 |
| `p018_nrows2` | B0 carrier + NWAKEROWS=2 (d ≈ 1.2σ) | 13011983 | submitted 2026-07-31 |
| `p018_nrows8` | B0 carrier + NWAKEROWS=8 (d ≈ 4.9σ) | — | as originally planned, when a job slot frees |

Compare all against `p018_b0` (N=4) on matched windows (M1 + M2). The
N∈{1,2,4,8} spread — not a limit claim — is what enters the error budget /
disclosure text. If the spread is large, that strengthens the case that the
sheet→particle handoff representation (016) is where the real convergence
lever lives; a re-run of the Das ladder at the chosen N may follow.

## STAGED — Das × N × shedding interaction study (Ryan directive 2026-07-31)

**Ryan's framing:** particles must be small enough that the body does not lie
inside the heavily regularized σ region. Higher Das pushes particles away from
the body, as do more sheet rows — but nwakerows is not physically meaningful
and should be minimal (1–2 preferred). Explore the interaction and seek
*convergent criteria*.

**Dimensionless groups.** The interaction reduces to three ratios plus one
choice:

$$d/\sigma \;(\text{body clearance}), \qquad \mathrm{Das}/c \;(\text{physical
attachment length}), \qquad h/\sigma = 1/\mathrm{OVERLAP}, \qquad N$$

with handoff distance $d \approx \mathrm{Das} + (N-1)\cdot(\text{TE
travel/step})$ (row lengths ≈ Das at η=1). Prior for the clearance threshold:
017 M4's static fit, particle-on-body velocity error $\sim (d/\sigma)^{-3.4}$.

**Key structural point (answers "should Das be larger?"):** clearance is
$\mathrm{Das}/\sigma$, and σ is *already* a campaign convergence axis while
Das is pinned physical. At fixed Das = 0.41c, the clearance improves
automatically down the σ ladder: $\mathrm{Das}/\sigma$ = **0.60** (σ/c=0.68)
→ **1.38** (L1, 0.297) → **2.72** (L2, 0.151). So enlarging Das is the
*coarse-σ workaround* — it buys clearance but moves the Kutta attachment
length (the 014 log law says dt-converged loading grows ~0.2%/doubling, and
the rotor ladder is not even monotone at this σ) — whereas refining σ buys
clearance without touching the physical model. The convergent limit is
**σ→0 at fixed physical Das on the 014 plateau (0.25–1.5c), N minimal**; a
larger Das should only be adopted if the coarse-σ data demand it AND the
Das-dependence at fine σ is shown flat (Phase 2's deferred re-test).

**Staged experiments (in order; respect the 6-job cap):**

1. **[in flight]** `p018_nrows1`/`p018_nrows2` vs `p018_b0`: N ∈ {1,2,4} at
   Das/σ=0.60/1.2/2.4 total-clearance. Large spread here = clearance-starved
   at coarse σ (expected); its *sign and magnitude* calibrate the error model.
2. **Coarse-σ Das×N cross matrix** (cheap, ~13 h each): N ∈ {1,2} ×
   η ∈ {2.0, 4.0} (Das 0.82/1.64c), i.e. up to 4 runs as slots free
   (`P018_RUN_NAME=p018_nrows<N>_das<η>`, carrier + NWAKEROWS + DAS_ETA
   overrides). **Discriminator:** plot all {N, Das} runs (incl. the Phase 2
   ladder) against $d/\sigma$. If CT̄ and Γ̄(r/R) collapse onto a single
   curve in $d/\sigma$, the Phase-2 non-monotonicity was clearance physics,
   Das itself is (nearly) inert above the 014 floor, and the criterion is a
   clearance threshold. If they do NOT collapse, Das carries independent
   physics (attachment length) and must be converged separately at adequate
   clearance.
3. **σ-refinement at minimal N:** rerun N=1 (and one N=2 confirm) at L1
   (Das/σ=1.38) and — if the pilot gate passes, warm — at L2 (2.72), fixed
   Das=0.41c. **Convergence criterion sought:** N-insensitivity
   (|ΔCT̄| ≤ 0.5%, ε_Γ ≤ 1% between N=1/2/4) emerging once $d/\sigma$
   exceeds a measured threshold; report the threshold vs the 017
   $(d/\sigma)^{-3.4}$ prior.
4. **Das plateau re-test at fine σ** (Phase 2's deferred branch, now at
   minimal N): top-two Das points at L1, N=1. Pass = doubling Δ ≤ 0.5%
   consistent with the ~0.2%/doubling law ⇒ quote Das model-form term from
   the fine-σ slope, not the coarse one.
5. **Acceptance for the production configuration:** N ∈ {1,2} (per Ryan's
   preference, smallest N passing 3's criterion), Das on the 014 plateau with
   fine-σ plateau evidence, σ\* from Phase 4's fit, and the N/Das/clearance
   sensitivities all ≤ tolerance at σ\* — each quoted in the error budget.
   If no clearance threshold produces N-insensitivity even at L2 ⇒ the
   sheet→particle handoff representation itself is the blocker ⇒ **016
   trigger** (pre-authorized, ruling 7).

## STAGED — hub-radius mesh variant runs (Ryan directive 2026-07-31)

**Implemented:** `scripts/generate_dji9443_mesh.sh --hub-r-over-r X`
(contingency (c) machinery, now real). Mechanics: OpenVSP clamps each XSec
RadiusFrac against its neighbors, and both *crowding* the inner sections and
*cutting the fat transitional root* produce broken lofts — the working recipe
is: CutXSec the interior sections inside the new hub, keep the stock root
section, move it out to X (`SetParmValUpdate`, outboard-first). Validation
matrix (45_185 flat/round ct4): **0.15 VALID** (watertight, TE auto-seeds OK,
36,752 elems = production count, true root r/R 0.136 from loft tilt);
0.18 watertight but TE-chain detection fails; 0.20/0.25 non-watertight
(fat-root loft folds mid-blade, worsening with hub size) — both rejected and
deleted. **Usable hub ceiling with this recipe ≈ 0.15R.**
Production mesh: `examples/data/dji9443_20260731_45_185_capped_captess4_hub0p15.msh`
(md5 `9d748d8c693242f37ac7b4663850b348`, verified on cluster). Blade root
r/R 0.136 vs stock 0.010; shedding root clip 0.1 becomes inert (root is
outboard of it) — the wake now sheds nowhere near the center of rotation.

**Three staged runs** (submit as slots free; b0 carrier + `RHPC_MESH_FILE`;
TE seeds auto-detect; 24 h each; pairings are stock-hub twins already run
or in flight):

```
sbatch --job-name=fp-018-hub-b0 --time=24:00:00 \
  --export=ALL,P018_RUN_NAME=p018_hub0p15_b0,RHPC_MESH_FILE=dji9443_20260731_45_185_capped_captess4_hub0p15.msh \
  examples/run_dji9443_hover_ct_hpc.slurm.sh p018_b0
sbatch --job-name=fp-018-hub-nrows1 --time=24:00:00 \
  --export=ALL,P018_RUN_NAME=p018_hub0p15_nrows1,NWAKEROWS=1,RHPC_MESH_FILE=dji9443_20260731_45_185_capped_captess4_hub0p15.msh \
  examples/run_dji9443_hover_ct_hpc.slurm.sh p018_b0
sbatch --job-name=fp-018-hub-nrows2 --time=24:00:00 \
  --export=ALL,P018_RUN_NAME=p018_hub0p15_nrows2,NWAKEROWS=2,RHPC_MESH_FILE=dji9443_20260731_45_185_capped_captess4_hub0p15.msh \
  examples/run_dji9443_hover_ct_hpc.slurm.sh p018_b0
```

**Why these three:** `hub×b0` isolates the pure hub-size effect at the
campaign baseline (vs `p018_b0`); `hub×nrows1` and `hub×nrows2` pair with the
in-flight stock `p018_nrows1/2` to answer the decision question — **does
removing near-axis shedding shrink the N-sensitivity enough to make Ryan's
preferred N=1–2 viable?** If ΔCT̄(N=1↔4) collapses under the hub mesh but not
the stock one, near-root particles were the driver; if the N-sensitivity is
unchanged, the clearance problem is outboard (σ-driven) and the hub is a null
(record ≤0.3% and move on). Caveats to carry: the hub mesh removes real blade
area inboard of 0.136R (CT̄ is NOT directly comparable to experiment — only
hub-vs-stock *deltas* and Γ̄(r/R) overlays outboard of 0.2R are meaningful),
and M2 comparisons must interpolate to the common r/R grid (mesh section
stations differ).

## Log

- 2026-07-31 — Ruling 11 recorded; `p018_nrows1` (13011982) and `p018_nrows2`
  (13011983) submitted at the b0 carrier (η=1.0, σ/c=0.68), 24 h each.
  Study-job count at the 6-cap with nt72/L1/L1_warm/b0_s2.
- 2026-07-31 — Das×N×shedding interaction study staged (section above) per
  Ryan directive; step 2's cross matrix queues as slots free; steps 3–4 gate
  on Phase 4's L1/L1_warm results.
- 2026-07-31 — `--hub-r-over-r` implemented + validated (section above);
  hub0p15 mesh generated, shipped, md5-verified; 3 hub runs staged (submit
  when nrows1/nrows2/b0_s2 free slots). 0.18/0.20/0.25 variants rejected
  (loft defects) and deleted.
- 2026-07-31 — `p018_hub0p15_b0` first attempt (13015838) **FAILED at step 0**:
  `basis=:quad ... agglomerate 1094, stencil_size=0` — the hub mesh's tilted
  root slivers (~126 panels of area ~1e-4·median, r/R≈0.169) fail the
  μ-gradient fold test against ALL neighbors ⇒ isolated agglomerate ⇒ empty
  stencil ⇒ cond=Inf throw. The stock mesh carries the same *count* of
  slivers (136 at r/R 0.057) whose normals happen to stay admissible, which
  is why it never trips. **Fix (src):** `_quad_mu_diff_gradient!` now gives
  an *isolated* agglomerate (empty stencil — no growth depth can help) a zero
  μ gradient with a `@warn` instead of aborting; non-isolated non-finite
  reconstructions still throw. Justification: zero-area slivers contribute
  negligibly to force; basis stays :quad everywhere ⇒ hub-vs-stock
  comparability preserved. Gate: full runtests + formulation_test rerun, then
  md5 redeploy + resubmit.
- 2026-08-01 — **`p018_nrows1` (13011982) COMPLETED (11.1 h) — THE N-AXIS IS
  NOW BRACKETED AT σ/c=0.68.** CT̄(10–19) = **0.07133 ± 0.00014, −0.75% vs b0**,
  ε_Γ,max 2.77% / RMS 1.33% ⇒ **FAILS both M1 (≤0.5%) and M2 (≤1%)**.
  Combined with nrows2's clean pass, the N-sensitivity table at the campaign σ:

  | N | handoff d | CT̄ | ΔCT̄ vs b0 | ε_Γ,max | verdict |
  | --- | --- | --- | --- | --- | --- |
  | 1 | 0.6σ | 0.07133 | −0.75% | 2.77% | FAIL both |
  | 2 | 1.2σ | 0.07188 | +0.01% | 0.35% | **PASS both** |
  | 4 | 2.4σ | 0.07187 | — (ref) | — | reference |

  **Reading:** exactly the clearance behavior Ryan's framing predicts — the
  body sits inside the regularized σ region only at the smallest handoff, and
  the effect switches on sharply between 0.6σ and 1.2σ (monotone in d,
  right sign: too-close particles UNDER-predict thrust). Qualitatively
  consistent with 017's (d/σ)^−3.4 error prior. **Practical answer to
  "N=1 or 2 if possible": N=2 IS admissible at this σ — adopt N=2 as the
  provisional production choice** (Ryan's preference satisfied without a
  criterion violation), pending the σ-refinement check: since clearance is
  d/σ, N=1 should become admissible at finer σ (at L1, N=1 gives
  d/σ=1.38 > the bracketed threshold) — test at L1 before finalizing.
  Note the N=1 value coincides numerically with the legacy `p2e_das1p0`
  0.07133; unrelated (different shedding policy), do not read as corroboration.
  Hub twin `p018_hub0p15_nrows1` submitted (13016742) to test whether the N=1
  failure is driven by near-axis particles.
- 2026-08-01 — **`p018_nrows2` (13011983) COMPLETED (11.6 h): CT̄(10–19) =
  0.07188 ± 0.00017 — NULL vs b0's matched window (0.07187, Δ +0.01%).**
  N=4→2 (handoff 2.4σ→1.2σ) is inert at σ/c=0.68: first evidence that the
  N-sensitivity Ryan's clearance concern predicts, if present, lives below
  1.2σ (nrows1, d=0.6σ, still running decides). Hub twin `p018_hub0p15_nrows2`
  submitted (13016670). M2 ε_Γ vs b0 pending with the batch comparison.
- 2026-07-31 — **Gate closed and hub run resubmitted (13016062).**
  formulation_test 10/10; full suite green after fixing three UNRELATED
  pre-existing latent failures the deeper-running suite unmasked:
  (1) postprocess unit test updated to the new isolated-agglomerate contract
  (sentinel-accumulation form — compute_mu_gradient! accumulates);
  (2) FLOWVPM `Estr_direct_multithreaded` empty-field guard (the documented
  "small-N SFS direct bug" is precisely np=0: `divrem(0,T)` → zero-step
  range; one-line fix in FLOWVPM.jl — the T4 SFS_OFF workaround is now
  unnecessary for future backend A/Bs); (3) missing `cross` import in
  `test/runtests_unit_simulate.jl`. **Exception: `data/pitching_wing_exp/`
  fixture is MISSING** (present in session-start git status, untracked, gone
  from disk + cluster + Trash; no session command touched it) — its unit test
  cannot run; recovery = Ryan restores via Dropbox web deleted-files. All
  other suite files pass (unit + example + analytical, run individually past
  the missing-fixture stopper). Deployed with md5 verify: FLOWPanel
  `src/FLOWPanel_postprocess.jl`, FLOWVPM `src/FLOWVPM_subfilterscale_models.jl`
  (cluster `~/projects/FLOWVPM.jl`); full src/ sweep matches.

- 2026-08-01 — **`p018_hub0p15_b0` (13016062) COMPLETED 13:05:13.** Status
  bookkeeping only: the run finished on the cluster and **its output has not
  been retrieved, read, or interpreted** (recorded per Ryan's instruction to log
  completion without touching results). The "action on completion" for this job
  is therefore still outstanding — hub-vs-stock delta at N=4, remembering that
  the hub mesh removes blade area inboard of r/R 0.136, so only DELTAS and
  Γ̄(r/R) outboard of 0.2R are comparable, never absolute CT against experiment.
  Its two hub twins (`nrows2` 13016670, `nrows1` 13016742) were still running at
  the time of this entry, and the payoff question — whether the hub mesh rescues
  N=1 — needs 13016742, so no hub conclusion can be drawn from this job alone.
