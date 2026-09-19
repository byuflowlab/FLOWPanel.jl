# FGS production dagteam status — gate-1 closed 28/28, FLOWPanel plumbing done (2026-09-19d)

Continues `fgs_acceleration_reset_prompt_20260919d.md` (TASK 2 production
implementation). Spec = `fgs_acceleration_recommendation_20260918.md`; measured
basis = `fgs_acceleration_status_20260919c.md` (champion dagteam + F32full +
interleave 0-3 @ t16 → 2.10×).

## 1. Multi-system oracle: RESOLVED — fixture class, not a code bug

The 2 failing gate-1 tests (multi-system strength-recovery oracle → NaN for
lex AND dagteam) are closed. Controls run locally (-t4, worktree
`/private/tmp/fastmultipole-p021-fgs-accel-20260918`):

| Control | Result |
|---|---|
| union-600-as-ONE-system, leaf40 (reset-prompt discriminator) | DIVERGES (NaN) — identical to 2-system |
| union-600 leaf100, start at truth | dev = 0.0 (VACUOUS — solve exits at iter 0 without moving) |
| 2-system leaf100, start at truth (lex + dagteam) | dev = 0.0 (same vacuity) |
| 2-system / union leaf100, ZERO start | dev = 3.1e219 — diverges identically for lex, dagteam, and single-system union |
| fixed point, all-direct MAC 0.01, 4 sweeps | dev = 4e8 (amplification ~5e5 per sweep) |
| fixed point, all-direct MAC 0.01, EXACTLY 1 sweep | dev = 3.8e-8 (lex) / 4.0e-8 (dagteam) |

**Diagnosis:** the 1/r first-kind gravitational potential operator has
block-GS spectral radius ≫ 1 for ANY multi-leaf partition (leaf self-blocks
have zero diagonal — self term excluded at r²=0). A multi-sweep recovery
oracle can never converge on this fixture class; the ONLY convergent recovery
test in the existing suite (`solve_test.jl` "full solve") uses
`leaf_size=n_bodies`, i.e. a single-leaf dense direct solve. Divergence says
nothing about multi-system plumbing — the multi-system fill fix from
`29a55bf4` stands (all-direct control: one sweep from truth moves 3.8e-8).

**Fix (committed `f4d6b671`):** replaced the recovery oracle with a
**one-sweep fixed-point oracle**: all-direct (MAC=0.01, m2l empty asserted),
`max_iterations=1, inner_iterations=1, tolerance=0.0` (forces exactly one
sweep; `reverse_pass=false` default), true strengths must remain a fixed
point to ≤1e-6. Sensitivity is proven inside the test by negative controls:
a 1% scaling of the nonself operator (lex: `nonself_matrices.data`; dagteam:
all `Lmat`/`Umat` blocks — `Lmat[1]` alone is an empty root block, gotcha)
must trip the oracle, and does at dev ≈ 1.6e4 for both schedules.

**Gate-1 now 28/28 at -t4** (was 24/26; +2 fixed-point, +2 negative
controls, −2 old recovery).

## 2. Regression suites (plan step 2) — all green

- `test/fgs_rowpar_gate1_test.jl`: 15/15.
- `test/solve_test.jl` (full FGS portion, runtests preamble): every testset
  passes (562121 + 3066 + 100 + 20000 + 9 + 9 + 157).

## 3. FLOWPanel plumbing (plan step 3) — done, smoke-verified

Live FLOWPanel checkout (uncommitted, this thread's edits only):

- `src/FLOWPanel_solver.jl`: `FGSSolver` gains a `dagteam_precision::Symbol`
  field + kwarg (`:f64|:f32conv|:f32full`, default `:f64`), forwarded to
  `FastMultipole.FastGaussSeidel` **only when `sweep_order === :dagteam`** so
  the live env (FastMultipole `c18e4b46`, pre-dagteam) keeps working for all
  other sweep orders. `FGSPreconditioner` passes both through.
- `src/FLOWPanel_metadata.jl`: `dagteam_precision` recorded in both solver
  metadata dicts.
- `benchmark/fgs_cold_common.jl`: `sweep_order="dagteam"` accepted;
  `dagteam_precision` config key validated (`f64|f32conv|f32full`, requires
  `sweep_order=dagteam`, default f64) and forwarded by `cold_make` — same
  pattern as v22 `chunks`.

Smoke (scratch env deving FLOWPanel + FLOWVPM + the p021 worktree; sphere
source body, P4/MAC0.4/leaf20, tol 1e-10): dagteam f64 rel-dev vs lex
3.3e-15; f32conv 1.6e-9; f32full 3.4e-7. All three construct through
`FGSSolver` and solve. (F32 accuracy certification remains the independent
evaluator's job at R4 — gate NOT relaxed.)

Live-env compatibility verified: with the guard, a `:lexicographic`
`FGSSolver` constructs and solves in the LIVE FLOWPanel env against the old
FastMultipole `c18e4b46` — other threads using the live checkouts are
unaffected by the plumbing. Cold-harness validation unit-checked all four
ways (dagteam+f32full accepted; `dagteam_precision` without dagteam
rejected; `f16` rejected; dagteam default-precision accepted).

## Worktree state

FastMultipole `p021-fgs-accel-20260918` @ **`f4d6b671`** (oracle fix, on top
of `29a55bf4`). FLOWPanel live checkout carries the plumbing edits
uncommitted (mixed with other threads' state — commit needs Ryan's routing
call or a clean stage of just these 3 files).

## Numerical gate (job 13773687): PASSED 2026-09-19 — f32full certified at R4

Ryan approved launching both jobs; provenance =
`fgs_acceleration_provenance_20260919d.md` (pins, tags, decision rules,
submission record). The qos=test calibrate-only job completed clean
(COMPLETED marker, empty .err): staircase-calibrated **dagteam + f32full at
tolerance 3.43e-7 (27 iters, repeat delta bitwise 0.0)** and the colored
baseline twin at 5.22e-7 (26 iters — matches v21's accepted 26), every
accepted solve passing the independent evaluator at BC rel-L2 ≤ 1e-6, at the
champion placement (interleave 0-3 / cpunodebind 0-3, j16). **The champion
precision rung holds; no fallback needed.** Local M2 pre-submit gate
independently certified the same config (dagteam tol 8.9e-8, evaluator
~1.0e-6 accepted) before dying in its confirmation solve (local memory,
non-blocking). Run dir
`data/p021-cold-20260910/dagteam-numgate-13773687/` on orc, harvest owed.

## End-to-end A/B (job 13773689): PASSED 2026-09-19 — 2.26× vs the accepted baseline

COMPLETED 3:01:47 on m12-2-11, empty .err, COMPLETED marker present; all
controls green (incl. standalone gate-1 28/28 on cluster); every trial
evaluator-certified (≤1e-6) with repeat agreement ≤1e-8 enforced; iterations
stable (colored 26 / dagteam 27 in every arm); cross-order solution
agreement 2.0e-6 (informational). 40 trials per order per arm. Evidence
harvested to `fgs_r4_followup_evidence_20260914/dagteam-v23-13773689/` and
`dagteam-numgate-13773687/`.

| Arm | colored median s | dagteam median s | within-arm × |
|---|---|---|---|
| **j16-interleave (decision)** | 7.243 | **4.475** | **1.62** |
| j16-native (v21 anchor) | 9.434 | 5.951 | 1.59 |
| j64-interleave | 10.481 | 6.223 | 1.68 |

Verdict against the spec gates (baseline = accepted colored@j16 =
**10.116 s**, v21):

- **dagteam f32full @ j16-interleave = 4.475 s = 2.26× — beats the ≤5.058 s
  design target** (required ≤6.744 s passed with margin); measured better
  than the 4.81 s projection.
- The fixed within-arm decision rule (≥1.5) passes at **1.62×**: interleave
  placement alone accelerates the INCUMBENT colored to 7.243 s (a free
  1.40× from a launcher change — adopt for the colored fallback regardless).
- Anchor: j16-native colored 9.434 s vs v21's 10.116 s (6.8% faster, same
  26 iterations; node variation m12-2-11 — within acceptable drift, flagged
  not disqualifying).

**Recommendation: PROMOTE dagteam + f32full + interleave 0-3 @ j16 as the
production R4 FGS configuration** (both spec gates passed: numerical +
end-to-end); adopt interleave placement for the colored fallback rung too.

## PROMOTED (Ryan, 2026-09-19)

Ryan approved promotion the same day. Executed:

- Live FastMultipole production lineage (`../FastMultipole`, branch
  `flowpanel-20260817`) fast-forwarded `c18e4b46` → **`f4d6b671`** (clean ff;
  additive: dagteam executor + gate tests + replay harness + the
  multi-system fill fix; lex/colored paths behaviorally unchanged). The
  FLOWPanel live Manifest now loads dagteam-capable FastMultipole.
- Calibrated champion committed as **`benchmark/retained_r4_champion.toml`**
  (dagteam/f32full, tolerance 3.4309419310610173e-7, schema-validated;
  placement `numactl --interleave=0-3 --cpunodebind=0-3` @ j16/BLAS 1
  documented in the file header; colored fallback tolerance
  5.223162079893313e-7 recorded there).
- Live-env promotion smoke (sphere, dagteam f64/f32conv/f32full vs lex
  through FGSSolver against the merged live FastMultipole): result in the
  reset prompt / next status.

Origin pushes of the merged lineage remain Ryan-gated (standing ledger).

## Next (Ryan-gated)

1. **Numerical gate** (plan step 4): independent evaluator BC rel-L2 ≤ 1e-6
   on R4 decides :f32full vs fallback rung. HPC submission awaits Ryan.
2. **End-to-end A/B** (plan step 5): interleaved vs colored@j16 champion,
   ≥1.5× accepted throughput, full campaign ceremony (tag pins, worktrees,
   provenance before submission).
3. Ryan-pending ledger unchanged: origin pushes (incl. this branch), 4
   notebook entries, WeakKeyDict/warmstart fix.
