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

> **RETRACTED 2026-08-03 — every claim in the paragraph above is wrong.**
> "0.15 VALID" was validated on watertightness + element count, which a
> 5.5×-oversized root chord passes cleanly; the mesh was in fact broken and
> produced CT = −0.474 at the first timestep. The "usable hub ceiling ≈ 0.15R"
> and the 0.18/0.20/0.25 rejections were **artifacts of the same generator
> bug**, not a property of hub size: with the bug fixed, hubs up to **0.40R**
> generate and pass the new gate. See "Root cause and fix" below.**
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

## 2026-08-03 — all three hub runs are VOID (mesh defect, not physics)

`p018_hub0p15_b0` (13016062), `p018_hub0p15_nrows2` (13016670) and
`p018_hub0p15_nrows1` (13016742) all completed and all produced **negative
thrust**:

| run | N | CT̄ (revs 15–19) |
| --- | --- | --- |
| `p018_hub0p15_b0` | 4 | **−0.2989** |
| `p018_hub0p15_nrows2` | 2 | **−0.3125** |
| `p018_hub0p15_nrows1` | 1 | **−0.2655** |

**These are not physical hub effects and must never be read as deltas.** The
diagnosis is unambiguous and does not require the wake at all:

- `CT_bernoulli = −0.474` at **step 1** of `p018_hub0p15_b0`, versus **+0.070**
  at step 1 of stock `p018_b0`. At step 1 no wake exists, so the defect is in
  the body's own geometry/pressure integration, not in shedding, clearance, or
  particle physics.
- Magnitude rules out the obvious explanation: a global normal/winding flip
  would give ≈ −0.07, not −0.47. Something is wrong with the surface itself,
  ~6.7× in magnitude and inverted in sign.
- Only 16 `basis=:quad found an isolated agglomerate` warnings appear in the
  run, so the root slivers the `FLOWPanel_postprocess.jl` guard was added for
  are **not** a sufficient explanation on their own.
- The histories are smooth and settled (per-rev ptp ~0.002–0.02), so this is
  not a divergence or an OOM-in-merge case — the solver is happily converging
  on a broken body.

**Consequence for the campaign:** the hub test's payoff — *does removing
near-axis shedding rescue N=1?* — is **unanswerable** until the mesh is fixed.
Per Ryan (2026-08-03), the mesh is to be diagnosed but **not** re-run; no
cluster slots are committed to this branch.

Staged for Ryan's ParaView inspection at `~/hub0p15_paraview/`:
`p018_hub0p15_b0_body1.0.vtu` (step 1, CT −0.474),
`p018_hub0p15_b0_body1.700.vtu` (settled), `p018_b0_body1.0.vtu` (known-good
stock control at the same step), and both `.msh` files.

Note for whoever picks up the diagnosis: the hub `.vtu` (3.8 MB) and `.msh`
(1.6 MB) are the **same size as stock**, so the hub recipe appears to have
preserved topology and moved nodes rather than removing panels — consistent
with the generator's "cut interior XSecs inside the hub, then move the stock
root out" description, and a natural way to end up with folded or inverted
cells at the relocated root. Check in this order: (1) per-cell normal
orientation and total signed volume vs stock; (2) count of degenerate /
near-zero-area cells; (3) whether the blade itself was rescaled (which is the
only candidate that explains a magnitude change a sign flip cannot);
(4) `scripts/generate_dji9443_mesh.sh --hub-r-over-r`. Remember the standing
`CLAUDE.md` invariant: shedding must come from the *constructed* body's
re-wound `.cells`, never raw mesh cells — a hub mesh that changes winding is
exactly where that bites.

### Revised N preference (Ryan, 2026-08-03)

Superseding the earlier "minimal N" framing: **prefer N=2 over N=1 if N=2
permits a significantly smaller number of particles shed per step.**

Rationale, and why it fits what this phase measured: N sets the handoff
distance d ≈ N·Das, and the failure mode here is a *clearance* one in d/σ
(N=1 at d=0.6σ FAILS: −0.75% CT, ε_Γ 2.77%; N=2 at d=1.2σ PASSES: +0.01%,
0.35%). Buying clearance with N rather than with σ permits a **coarser
admissible σ**, and since σ = 2πR·OVERLAP/(NT·PPS), coarser σ means a smaller
`P_PER_STEP` — directly fewer particles and a much cheaper run (cost model
≈ 9.0 h/(σ/c), ≈ 52.4k particles/(σ/c)).

This **re-frames the deferred N=1-at-L1 re-test** (Phase 4): the question is no
longer "is minimal N attainable?" but "does N=1 buy anything N=2 does not,
given N=2 lets us shed fewer particles?" If N=1 needs a finer σ (larger PPS)
to reach equal accuracy, N=2 wins outright. It refines ruling 11's direction
("try N=1 and N=2"); N remains a modeling choice whose sensitivity is measured
and disclosed, not driven to a limit.

### VTK retention (ruling 10), 2026-08-03

Deleted VTK (kept CSV/TOML/monitors) for `p018_hub0p15_nrows1` and
`p018_hub0p15_nrows2` — both void and redundant with `p018_hub0p15_b0`, whose
VTK is **retained** as the diagnostic copy while the mesh defect is
investigated. Local harvest of all three completed first.

Deviation from the "≤3 runs' VTK" cap, stated rather than glossed: several
harvested runs (`p018_b0`, `p018_das*`, `p018_L1*`, `p018_nrows*`) still hold
VTK. Cluster disk is at 9% of 2.0 T, and the σ-axis conclusions these runs
underpin have just FAILED M2 (ε_Γ 8.8%) and are likely to be re-examined, so
deleting them now would trade a real diagnostic asset for space that is not
scarce. Revisit once Phase 4 closes.

## 2026-08-03 (b) — N=1's failure is INBOARD-localized; the clearance premise is confirmed

Figure: `data/p018_gamma_ladders/gamma_nrows.png`.

ΔΓ̄ vs B0 (N=4), matched revs 10–19, normalized by max_r|Γ̄_B0|:

- **N=1 (d = 0.6σ):** a monotone **inboard deficit** — **−2.8% at r/R 0.31**,
  decaying through −1% at r/R ≈ 0.58, crossing zero near 0.65, and slightly
  *positive* (+0.4%) outboard. Band integral −0.809%.
- **N=2 (d = 1.2σ):** flat within ±0.35% everywhere. Band integral +0.077%.

**This confirms the premise the hub test was built to probe — independently of
the hub mesh, which turned out to be defective.** The N=1 error is not spread
over the blade; it is concentrated exactly where the near-axis hypothesis said
it would be.

### Why inboard, mechanistically

`Das = η·Δt·U_local` with `U_local = Ωr`, so the handoff distance
**d = N·Das ∝ r**. The clearance ratio d/σ is therefore *intrinsically smallest
inboard* — a single global N is a **radially varying clearance**, tightest at
the root. The `min_displacement` floor (0.01R) is active in that same inboard
region, which is where it was already known to bite hardest. So N=1 fails
inboard first not because of anything special about the hub, but because that
is where d/σ runs out first.

### Cross-axis corroboration (the strong form of the argument)

The Das ladder's small rung shows the **same inboard signature**: Das = 0.205c
gives −8.9% at r/R 0.31 decaying outboard (phase_02, 2026-08-03). N=1 and
Das = 0.205c are *both* small-`d` cases — one halves N, the other halves Das,
and `d = N·Das`. That two independent ways of shrinking the same product
produce the same spatially-localized error is direct evidence that **d/σ is the
governing group**, which is exactly what the staged Das×N cross matrix was
designed to test. That matrix is now much better motivated, and it has a sharp
prediction to check: the two ladders' Γ̄ deficits should **collapse onto one
curve when plotted against d/σ**.

### Consequence for the N decision (and ruling 12)

Because the clearance deficit is inboard and the inboard region carries little
thrust (Γ̄ is small there and the moment arm is short), N=1's −0.75% CT̄ penalty
is the *integrated shadow* of a locally much larger (−2.8%) circulation error.
Under ruling 12 the question is whether N=2 lets us shed fewer particles; the
answer this figure adds is that **N=1's deficit is structural and radially
systematic, not scatter** — so it will not wash out at finer σ merely by
averaging. It may still clear once d/σ rises (at L1, d/σ = 1.38), and that
remains the test worth running.

## 2026-08-03 (c) — hub mesh: root cause found and generator FIXED

Ryan's ParaView read ("the mesh looks crazy near the root — as if the chord
explodes and twist jumps") is exactly right, and is now measured, explained,
and fixed.

### Measurement

New tool `scripts/p018_mesh_profile.py` extracts chord(r/R) and twist(r/R) from
the OpenVSP loft rings in a `.msh`. Stock vs the broken hub mesh:

| ring | stock | hub0p15 (broken) |
| --- | --- | --- |
| 0 (root) | r/R 0.0100, chord 0.5400, twist 13.19° | r/R 0.1493, **chord 0.5400, twist 13.19°** — the r/R=0.01 section, moved but unchanged |
| 3 | r/R 0.0217, chord 0.5922 | r/R **0.1455** — *inboard of the root* — **chord 7.2546** |
| 4 | r/R 0.0442, chord 0.6599 | r/R 0.1695, **chord 7.1787** |

Max chord anywhere on the stock blade is **1.3276**, so 7.25 is a **5.5×
balloon**. Panel areas corroborate: max/median panel area is **70.8** on the
broken mesh vs **3.6** on stock.

### Root cause

`scripts/generate_dji9443_mesh.sh` moved the root XSec's `RadiusFrac` out to
`hub_r` but left the blade's **distributions** untouched. Those distributions
are PCurves keyed on r/R starting at the stock root 0.01, and `Update()` only
**clamps the first knot's t** to the new root — it keeps that knot's *value*
and leaves every other knot where it was. So the root inherited the r/R=0.01
chord/twist, and knots at r/R 0.10 and 0.14 ended up **inboard of the root**,
folding the loft; the spanwise skin then overshot wildly bridging the leftover
step. `CT = −0.474` at step 1 follows directly.

Critically, **it is not just chord and twist.** A PropGeom carries nine
distributions (chord, twist, rake, skew, sweep, thickness, CLi, axial and
tangential offsets), all keyed from 0.01. Fixing only chord+twist moved the
root section to the right values but still left a **2.35 chord balloon** in the
first spanwise segment, because thickness and the offset curves were still
folded. Only trimming *all nine* removed it.

### Fix (implemented)

Correct order of operations, now in the generator:

1. **Split every distribution at `hub_r`** while its domain still covers it,
   and keep only knots at/outboard of `hub_r`. The split point carries
   OpenVSP's own evaluation of the original curve at `hub_r`, so the retained
   blade is the stock blade restricted to r ≥ hub_r — **truncation with a
   re-anchored endpoint, not deletion**. (This is the answer to Ryan's
   question about dropping inboard data: the inboard *knots* go, but no
   physical information is lost, because the new root knot carries the original
   curve's value there. Splitting a CEDIT curve at `hub_r` makes that point a
   segment endpoint, so the retained count stays 3n+1 as CEDIT requires.
   `PCurveDeletePt` is a no-op on the leading point, hence the rebuild via
   `SetPCurve`.)
2. Cut the XSecs inside the hub; copy the first surviving outboard section's
   **airfoil** onto the root (`CopyXSec`/`PasteXSec` — replaces shape only,
   leaves `RadiusFrac`), since the r/R=0.01 shank shape describes a region the
   hub has removed; then move the root out to `hub_r`.
3. **Only then** install the trimmed distributions, so the clamp in (1) is a
   no-op — the first knot already equals the root station.

### Validation gate (the one that was missing)

`scripts/p018_mesh_profile.py check` — watertightness and element count cannot
see this class of defect, and demonstrably did not. The gate asserts: no ring
folds inboard of the root; no chord exceeds the stock maximum; **outboard of
the hub the blade reproduces stock** (it must differ only by removed inboard
area, or a hub-vs-stock delta is uninterpretable); the first interior ring
matches the stock distribution; and the root is not still the stock r/R=0.01
section. It **fails the old mesh on 3 checks and passes the new one on all 7**.

One subtlety encoded in the gate: compare *interior* rings, not the root cap.
A cap ring measures the true 2D chord (stock's own cap reads exactly
PCurve × R = 0.5400) while interior rings measure a 3D LE–TE distance ~4%
higher because of sweep — comparing a cap ring against interpolated interior
rings fails by that systematic offset alone.

### Result

New mesh `examples/data/dji9443_20260803_45_185_capped_captess4_hub0p15.msh`,
36,752 elements (production count):

| quantity | stock | hub OLD (broken) | hub NEW (fixed) |
| --- | --- | --- | --- |
| root section chord / twist | 0.5400 / 13.19° @ r/R 0.010 | **0.5400 / 13.19° @ 0.149** | **1.1029 / 20.15° @ 0.150** |
| max chord | 1.3276 | **7.2546** | **1.3276** |
| max/median panel area | 3.6 | **70.8** | **3.6** |
| sliver panels (<1e-3 median) | 136 | 126 | **64** |
| first interior ring vs stock | — | 7.1787 vs 1.2161 | **1.1820 vs 1.1858** |

The fixed mesh is better conditioned than stock on sliver count, which also
means the `FLOWPanel_postprocess.jl` isolated-agglomerate guard (added for this
mesh's root slivers) should rarely fire on it. **Keep that guard** — it is
correct on its own terms — but note it suppressed the only automatic signal
that the geometry was wrong, so its silence must not be read as health.

### Hub ceiling re-measured — the old limit was an artifact

With the fix, `--hub-r-over-r` **0.20, 0.25, 0.30 and 0.40 all pass the gate**
(previously 0.18 failed TE-chain detection and 0.20/0.25 were non-watertight).
The "usable hub ceiling ≈ 0.15R" is retracted; hub size is no longer the
binding constraint.

### Not done, deliberately

No hub case has been resubmitted (Ryan's ruling: fix and record, do not
revive). **The mesh is validated geometrically but has not been run**, so the
CT-sign confirmation — the thing that exposed the bug — is still outstanding
and needs a cluster job whenever the hub branch is revived. The old broken mesh
is retained on disk as the gate's negative fixture; it must never be run again.

### Fixed-mesh hub run submitted 2026-08-03 (job 13031568, `p018_hubfix_b0`)

Ryan authorized one hub run on the fixed mesh. Deployed and md5-verified on the
cluster: mesh `dji9443_20260803_45_185_capped_captess4_hub0p15.msh`
(`14857a4d08c8fc1e654c29c0370ef549`) and the driver
(`b4ac5ef1109479b07fbbc0831fbd108a`, which differs from the previous cluster
copy only by the BRAINSTORM/016 `CONVERSION` switch, defaulting to `legacy`).

Submitted as the **B0 carrier with only the mesh swapped**, so it is directly
comparable to `p018_b0` = 0.07170. Banner-verified knobs: NT=36, η=1.0,
OVERLAP 2.0 / PPS 4, merge_r 0.0120, rlxf 0.3, **nwakerows=4** — identical to
B0. New run name `p018_hubfix_b0` deliberately avoids the void
`p018_hub0p15_b0` directory.

**Why N=4 and not the N=1 "payoff" run.** With a single run, hub@N=1 is
uninterpretable: its delta would mix the hub effect with the N effect, and the
hub reference at N=4 does not exist yet. hub@N=4 (a) is the fixed mesh's first
test *in the solver* — the outstanding item from the geometry fix, since the
gate proves geometry but not physics — and (b) establishes the pure hub effect
that every later hub×N cell is measured against. The N=1 payoff needs its pair
and should be submitted next; one slot remains under the 6-job cap.

**First check on this run is a sign test, not a convergence test.** The broken
mesh gave `CT_bernoulli = −0.474` at step 1 (`CF_x = +0.47`; thrust ≈ −CF_x).
The fixed mesh must give CF_x ≈ −0.07 like every healthy run. Interpretation
caveat (unchanged): the hub removes blade area inboard of r/R ≈ 0.15, so
compare **deltas** and Γ̄(r/R) outboard of 0.2R — never absolute CT against
experiment. Note the M2 band (0.3 ≤ r/R ≤ 0.95) lies entirely outboard of the
hub cut, so band Γ̄ comparisons are meaningful.

#### Result (job 13031568): geometry fix VALIDATED in the solver, but the run DIVERGES

**The sign test passes decisively.** Step 0 `CFx = −0.08030`, i.e. thrust
**+0.080** — right sign, right magnitude class (stock B0 is +0.070 at step 1).
The broken mesh gave `CFx = +0.474`. **The geometry fix is confirmed physically,
not merely geometrically**; the `CT = −0.474` pathology is gone.

**But the run blew up within 17 steps** and was OOM-killed at 10:28 elapsed
(MaxRSS 66.9 GB vs 64 G request):

| step | 0 | 2 | 4 | 8 | 12 | 16 | 17 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CFx | −0.0803 | −0.1547 | −1.604 | −23.02 | −118.8 | −1725.7 | −3263.3 |

Roughly geometric growth (≈×2/step) from step 2 on. Per the standing diagnostic
rule this is **divergence, not a memory ceiling** — `|CF| = 3263` is recorded
*before* the kill, so the blow-up is proven rather than inferred, and the OOM is
the diverged cloud wrecking the merge spatial hash. **Do not resubmit with more
memory.**

Note the broken mesh did *not* diverge (it ran 20 revs at CT ≈ −0.30), but that
is not evidence the fix caused instability: the broken mesh's loading was
garbage throughout, so it never developed a physical root vortex to be unstable
about.

**Leading hypothesis, and it connects to today's Γ̄ finding.** The hub truncates
the blade at r/R = 0.15 where the stock distribution still carries substantial
circulation (chord 1.10, Γ̄ ≈ 0.22), so a **strong root vortex is shed right at
the hub edge**. Two things make that the worst possible place for it:

1. `shedding_r_over_R = 0.1` is now **inert** — the blade root (0.15) is
   outboard of the clip — so the innermost, full-strength strip sheds.
2. `Das = η·Δt·U_local` with `U_local = Ωr`, so at r/R = 0.15 the local speed is
   small, `Das` is floor-clamped by `min_displacement`, and **d/σ is at its
   worst exactly there** — the same radially-varying-clearance mechanism that
   this phase measured as N=1's inboard deficit.

So the hub configuration concentrates a strong shed vortex at the radius where
the clearance criterion is weakest. That is a coherent, testable explanation
rather than a generic "it diverged".

**Cheapest discriminators, none submitted** (Ryan authorized one run):

- `SHEDDING_R_OVER_R=0.2` — clip the innermost strip so the hub edge does not
  shed. This is contingency (b), an existing knob, and directly tests the
  root-vortex hypothesis.
- Raise `NWAKEROWS` (more clearance) or lower `DAS_MIN_DISPLACEMENT_R`.
- A smaller hub (0.20/0.25/0.30/0.40 now all generate cleanly) would *increase*
  root Γ; a hub *smaller* than 0.15 would decrease it, and is the gentler probe.

**Status: the mesh generator fix is done and validated; the hub branch itself is
now blocked on a shedding/clearance problem at the hub edge, not on geometry.**

#### ROOT CAUSE of the divergence (Ryan, from ParaView): the wake was shed from the ROOT CAP

Ryan's read of the downloaded files: *"Wake is being shed from the cap at the
hub, when it should only be shed from the trailing edge."* Confirmed
quantitatively by a geometry-only probe (no solve, no time stepping) that walks
`calc_shedding_from_seed` and bins every shed edge by radius:

| case | shed edges | r/R range | edges AT the blade root (= cap) |
| --- | --- | --- | --- |
| stock, cutoff 0.1 | **40** | 0.111 → 0.987 | **0** |
| hub 0.15, cutoff 0.1 (**the failed run**) | **136** | **0.1497** | **92** |
| hub 0.15, cutoff 0.2 | 41 | 0.217 → 0.988 | 0 |
| hub 0.15, auto-bumped 0.1997 (**the fix**) | **41** | 0.217 → 0.988 | **0** |

**68% of the shedding chain (92 of 136 edges) was the root-cap perimeter**, with
the innermost edge sitting exactly at the blade root radius. My earlier
"strong root vortex at weak d/σ clearance" hypothesis was wrong — the vortex was
not merely badly resolved, it was being shed off a **cap surface that has no
trailing edge at all**.

**Mechanism — a latent trap, not a one-off config slip.** `shedding_r_over_R`
(hard-coded 0.1) silently doubled as the guard that stops the TE chain from
wrapping onto the root cap: `calc_shedding_from_seed` walks the trailing edge
until the bbox cuts it off, and on the stock blade (root at r/R 0.0095) that
cutoff always did the job. **Any mesh whose root moves outboard of the cutoff
leaves the guard inert and the chain wraps the whole cap perimeter with NO
error** — the same class of silent-failure trap as the re-wound-cells invariant
in `CLAUDE.md`. This applies to *every* `--hub-r-over-r` mesh, so it would have
bitten the whole hub branch, not just this run.

Note the earlier BROKEN hub mesh also shed from its cap, but its cap was half
the chord (0.54 vs 1.10) and its loading was garbage, so it produced a wrong-sign
steady state instead of a blow-up. Both pathologies traced to the same root
region; only the fixed geometry made the cap shedding energetic enough to
diverge.

**Fix (in `examples/rotor_hover_pressure_comparison.jl`):**

1. `shedding_r_over_R` is now a real ENV knob, `SHEDDING_R_OVER_R` (it was
   hard-coded before — the campaign docs' claim that it was "an existing knob"
   for contingency (b) was **wrong**), plus `SHEDDING_ROOT_MARGIN` (default 0.05).
2. The cutoff auto-tracks the blade root:
   `effective = max(requested, blade_root_r_over_R + margin)`, with a warning
   when it fires. On stock this is a **no-op** (0.0095 + 0.05 < 0.1).
3. A **hard guard** counts shed edges at the blade root radius and `error`s if
   any survive. Cap shedding is silent and catastrophic, so it must abort, not
   warn.

**Verification:** stock 40_40 coarse smoke is **bit-identical** to the pre-change
run (`CT_bernoulli = 0.08175569303996208` at step 15), no auto-bump fires, guard
passes. On the hub mesh the auto-bumped cutoff (0.1997) yields 41 shed edges,
0 cap edges — matching stock's 40-edge TE topology.

**Status: the hub branch's two blockers are both now fixed** (mesh geometry, and
cap shedding). The fixed-mesh run has NOT been resubmitted. Before any new
submission the driver must be redeployed and md5-verified — the cluster copy
predates this guard.

#### CORRECTION 2026-08-03: the auto-bumped cutoff was the WRONG fix; anchor the trace instead

Ryan, inspecting the guarded run: *"there are two more trailing edge pairs, right
at the hub, that are not shedding."* Correct — the probe agrees: at cutoff 0.1
the hub mesh has **44 genuine (non-cap) TE edges**, and raising the cutoff to
`root + 0.05` kept only **41**. The guard traded cap shedding for silently
discarding real trailing edge.

**The underlying mistake was selecting TE edges by WHERE they are instead of by
WHAT they are.** `shedding_r_over_R` was doing three unrelated jobs at once:
separating the two blades, clipping the root as a modeling choice, and — only by
coincidence — stopping `trace_trailing_edge` before it wrapped onto the cap. Any
mesh that moves the blade root breaks that coincidence, and it fails **silently
in both directions**: too small ⇒ shed from the cap and diverge; too large ⇒
quietly drop trailing edge.

**The intrinsic criterion already existed and was already validated.**
`examples/dji9443_trailing_edge.jl` finds the TE by what it is — edges that are
predominantly radial (`|Δy|/len ≥ 0.7`) AND sharp (`1 − n̂₁·n̂₂ ≥ 0.2`) —
assembles them into a chain, requires a simple open path, and cross-validates by
tracing `outer → inner`. Cap edges run circumferentially, so they never enter the
chain. It returns `[outer, second_outer, inner]`; **the driver used only the
first two and threw `inner` away**, and `calc_shedding_from_seed` has accepted
`end_node` all along (`src/FLOWPanel_liftingbody.jl:1277`).

**Resolution — separate the two jobs, because they conflict.** The modeling clip
deliberately truncates the stock TE (whose true root node is at r/R 0.0095,
inboard of the 0.1 clip), so a walk *confined* by the clip can never reach
`end_node` — attempting that errors. So:

1. **Trace** the full TE anchored end-to-end (`end_node = te_indices[3]`), with
   the bbox reduced to blade separation only (cutoff 0). Cap-free by
   construction, and a wrong turn now **errors** instead of failing silently.
2. **Clip** afterwards as an explicit modeling choice, using the identical
   edge-midpoint criterion the bbox used — so the retained set, and therefore
   stock behavior, is unchanged.

Also replaced the regression guard: it now tests edge **direction**
(circumferential ⇒ cap) rather than edge **radius**, because the chain now
legitimately reaches the blade root and a radius test could not tell a root TE
edge from a cap edge.

`examples/dji9443_trailing_edge.jl`: the chain filter's `inner_radius ≤ 0.2·radius`
rejected every hub ≥ 0.2R. Replaced with a span test
(`outer_radius − inner_radius ≥ 0.5·radius`) that discriminates equally well
without assuming where the root is.

**Verification (all passed):**

| case | traced → kept | innermost r/R | circumferential (cap) edges |
| --- | --- | --- | --- |
| stock, clip 0.1, `end_node` | 44 → **40** | 0.111 | **0** |
| hub 0.15, clip 0.1, `end_node` | 44 → **44** | **0.1593** | **0** |
| hub 0.15, clip 0.1, **no** `end_node` | 136 → 136 | 0.1497 | **92** |

The third row isolates the cause: `end_node`, not the clip, is what excludes the
cap. Stock keeps exactly its production 40 edges at r/R 0.111.

5-step local runs (coarse 30×81; production 45×185 needs ~10 GB and the local box
has 16 GB): **stock is bit-identical** to the pre-change run
(`CT_bernoulli` 0.09042433795069622 → 0.1079982255879601), confirming campaign
comparability is untouched. Hub now sheds 29 edges from r/R 0.164 (was 0.223),
stays bounded, and sits ~2% below stock instead of ~4% — the expected direction
from recovering trailing edge. Tests: DJI TE detection **23/23**, RigidWakeBody
**42/42**. Detector now succeeds on 0.20/0.25/0.30/0.40 hub meshes.

**Not done:** `trace_trailing_edge` itself is unhardened (Ryan's scope call) — the
general fix, for cases with no `end_node`, would be to promote its `side_score`
from a sort key to a hard gate, or require the edge to run spanwise. And the
stock root clip is NOT re-baselined: shedding the full stock TE (r/R 0.0095
instead of 0.111) is more physically consistent but would move every 018 number.

## 2026-08-04 — Das×N cross matrix harvested: NO d/σ COLLAPSE; Das dominates

The pre-registered discriminator (STAGED step 2) is answered. Jobs 13035910–13
(`p018_nrows{1,2}_das{2p0,4p0}`, banner/metadata-verified: B0 carrier, NT=36,
σ/c=0.68, OVERLAP 2.0/PPS 4) harvested; matched windows revs 10–19 (all
settled by rev ~9). Tool: `scripts/p018_dasN_matrix.py` (reuses
`p018_analyze.py` loaders); figure
`data/p018_gamma_ladders/gamma_dasN_matrix.png`; rows in `ledger.md`.

Full 10-point field vs B0 (d/σ = 0.603·N·η at this σ):

| d/σ | cells (N·η) | CT̄ spread | ε_Γ,max vs B0 spread |
| --- | --- | --- | --- |
| 1.21 | N2·η1 / N4·η0.5 / N1·η2 | 0.07188 / 0.07006 / 0.07147 (**2.6%**) | 0.35 / 8.88 / 0.99 % |
| 2.41 | N4·η1 (ref) / N2·η2 / N1·η4 | 0.07187 / 0.07185 / 0.07016 (**2.4%**) | 0 / 1.01 / 5.90 % |
| 4.82 | N4·η2 / N2·η4 | 0.07230 / 0.07039 (**2.7%**) | 3.11 / 6.20 % |

**Verdict: the {N, Das} runs do NOT collapse onto d/σ — the branch "Das
carries independent physics (attachment length) and must be converged
separately at adequate clearance" is the one selected.** Sharper structure:

1. **The large-Das OUTBOARD deficit is a function of Das alone.** At
   Das=1.64c the r/R≈0.84 deficit is −5.90/−6.15/−6.61% for N=1/2/4 — the
   *first-row length*, not the total sheet extent N·Das, misrepresents the
   rolled-up tip vortex.
2. **The small-Das INBOARD deficit is also Das-dominated**: das0p5 (−8.88%
   at r/R 0.31) sits at the same d/σ=1.21 where N2·η1 is clean to 0.35% —
   so it is not explained by total-handoff clearance either.
3. **N is a secondary modifier with two distinct signatures:** the true
   clearance penalty exists only at d/σ ≲ 1 (N=1·η1, d=0.6σ, −0.75%), and a
   mild sheet-extent effect appears at Das ≥ 0.82c (N=2 vs N=4: −0.62%/
   −0.63% CT, ε_Γ 2.2/2.1% — **fails M1+M2**, so the "N=2 ≡ N=4 null" of
   2026-08-01 holds ONLY at Das\*=0.41c, where it is clean at 0.35%).
4. Consequences: (a) the production choice **N=2 at Das=0.41c stands on its
   own two-point null**; (b) **Das=0.82c's "live competitor" status is
   weakened** — it carries its own N-sensitivity; (c) the Das axis still
   needs its σ\*-re-test (Phase 2 no-plateau branch) — nothing here converges
   Das, it identifies which variable owns the error; (d) the earlier
   cross-axis argument (phase_05 2026-08-03 (b)) that "d/σ is the governing
   group" is CORRECTED: the shared inboard *location* of the small-d errors
   was real, but the magnitudes do not scale with d/σ.

### `p018_hub0p25_b0` (13032932) harvested — inboard region is aerodynamically minor

CT̄(14–19) = 0.068022 ± 0.00025 (deltas only — blade area removed inboard of
0.25R). Γ̄ vs B0 (matched 14–19): a steep **local** deficit at the new root
(−16% at r/R 0.29 → 0 by r/R 0.41 — the truncated blade's own root
unloading), and **|Δ| ≤ 1.2% everywhere outboard of r/R 0.45** (raw ε_Γ,max
10.7% is entirely that root-adjacent lobe; RMS 2.4%).

**Pre-registered question answered: the inboard region does NOT shape
outboard loading.** Deleting everything inboard of 0.25R leaves the outboard
Γ̄ at the ~1% level ⇒ the σ-axis 3-lobe error (dip at r/R 0.76) is not
root-driven, the tip-vortex reading stands, and inboard clearance stays a
footnote in the budget rather than a first-order term.

### Larger-hub HPC case submitted 2026-08-03: `p018_hub0p25_b0` (job 13032932)

Ryan asked for a single HPC case with a larger hub. **0.25R at N=4 on the B0
carrier** — only the mesh differs from B0 (0.07170), so the delta is the pure
hub effect.

**Why 0.25R and not 0.30–0.40** (all of which now generate cleanly): the M2 band
is 0.3 ≤ r/R ≤ 0.95, so a 0.30R hub would start exactly at the band's inner
edge, leaving no comparison margin. 0.25R removes substantially more disk than
0.15R while keeping the **entire** Γ̄ band on blade, so Γ̄ compares directly and
completely against B0.

**Why N=4 and not N=1:** with one run, hub@N=1 mixes the hub effect with the N
effect and has no same-hub reference to subtract.

**Pre-registered prediction** (from the 2026-08-03 Γ̄ analysis): d = N·Das ∝ r, so
clearance is intrinsically worst inboard — exactly the region a hub deletes. If
the inboard region is aerodynamically minor and the σ error really lives
outboard (tip vortex), CT̄ should fall modestly while **ε_Γ over 0.3–0.95 stays
small**. A *large* outboard ε_Γ would instead show the root region shaping
outboard loading, promoting inboard clearance from a footnote to a first-order
budget term.

**Pre-flight (all verified before submission):**

- Mesh gate **7/7**: max chord 1.3266 vs stock 1.3276; outboard blade matches
  stock to 0.68% chord / 0.31° twist; root section 1.3266 / 19.02° (not the
  stale-bug signature 0.5400 / 13.19°).
- Shedding topology, local probe **and** confirmed in the run's own log:
  **44 TE edges traced → 44 kept**, innermost at r/R **0.258**, **0
  circumferential (cap) edges**; the 0.1 clip is inert as intended.
- Deploy md5-verified for all three changed files (mesh `02e6ded8…`, driver
  `d19c514a…`, TE detector `0f8dce88…`).

**Early result: healthy.** `CFx = −0.0665` at step 0 (thrust **+0.067**, correct
sign — the broken 0.15R mesh gave +0.474 there), smooth and bounded through
step 348 (~rev 9.7) at +0.065. ~87 steps/h ⇒ ~12 h to a full 29-rev history.

Harvest as usual (M1 over ≥15 settled revs, M2 vs B0 on 0.3–0.95). **Interpret
deltas only** — the hub removes blade area, so absolute CT is not comparable to
experiment.
