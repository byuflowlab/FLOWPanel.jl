# FGS acceleration: consolidated summary (closed at promotion, 2026-09-19)

One-page story of the BRAINSTORM 021 FGS-acceleration campaign: what was
tried, what won, why, and where the evidence lives. The prepared cold R4
body solve (58,192-panel DJI 9443 mesh, retained P8/MAC0.4/leaf100/inner3)
went from **10.116 s → 4.475 s (2.26×)**, promoted to production 2026-09-19.

## Result

| Configuration | R4 solve (median) | vs accepted baseline |
|---|---|---|
| lexicographic @ j64 (pre-campaign reference) | 10.95 s | 0.92× |
| colored @ j16 (v21, accepted baseline) | 10.116 s | 1.00× |
| colored @ j16 + interleave 0-3 (v23 side finding) | 7.243 s | 1.40× |
| **dagteam + f32full + interleave 0-3 @ j16 (PROMOTED)** | **4.475 s** | **2.26×** |

Champion = `benchmark/retained_r4_champion.toml` (calibrated tolerance
3.43e-7); the `numactl --interleave=0-3 --cpunodebind=0-3`, 16-thread/BLAS-1
placement is part of the champion (without it the executor collapses to
~1.09×). Fallback ladder: colored+interleave (7.24 s) → lexicographic.

## Why it is faster (three multiplicative levers)

1. **Memory placement (~1.4× alone).** The nonself coefficient cache
   (2.86 GB streamed per sweep, 81 sweeps/solve) was serially first-touched
   onto one NUMA node; interleaving pages across socket 0's four controllers
   raises useful bandwidth for every schedule (NUMA diagnostic:
   29.4 → 74 → 164 GB/s synthetic ladder).
2. **Full-F32 sweep precision (~1.5–1.9× on the stream).** Coefficients,
   sweep state, and leaf LU in Float32 halve the streamed bytes; the earlier
   convert-on-load variant was limited by the convert kernel itself (gate 2c:
   serial 27.8 F64 / 34.0 F32conv / 52.2 F32full GB/s). Residual and outer
   bookkeeping stay F64; accuracy is certified externally (below).
3. **Dagteam split dual-layout executor (+53% over row-parallel at matched
   precision).** The nonself operator is split into a target-major lower
   triangle (readiness-counter pull-DAG over the 48,167 directed lower
   edges — leaves solve as soon as their predecessors publish, no full-team
   handoff per leaf) and a source-major upper triangle computed as backward
   filler and reduced at the sweep boundary. Mathematically the exact
   lexicographic block-GS iterate (gate-1: 28/28 vs lex, ~1e-14
   sweep-level), bitwise-deterministic per run at any thread count.

The measured stream at the champion is 128.2 GB/s F64-equiv (gate 2d),
against a ~3.0 s non-sweep remainder — which now dominates the 4.475 s
solve, so further large gains must come from iteration count or remainder,
not bandwidth (spec §disposition).

## What was tried and rejected

- **Chunked hybrid GS/Jacobi (v22)**: faster sweeps (38.8 active threads)
  but iteration inflation 27→44; loses end-to-end at every thread count.
- **Row-parallel F64 (gate 2b)**: best 54.8 GB/s → 1.40×, below the 1.5×
  minimum everywhere.
- **F32 convert-on-load as the terminal precision**: certified GO at 1.73×
  (rev b) but dominated by f32full once the convert cost was isolated.
- **Socket-membind placement**: control collapsed to 36.8 GB/s — interleave
  is load-bearing.

## Accuracy of the winning configuration (certified, gate never relaxed)

Accuracy is judged by the campaign's **independent evaluator** — a certified
FMM evaluation of the boundary condition, separate from the solver's
internal residual (which can flatter a rounded F32 operator). The standing
acceptance is **authoritative BC relative L2 ≤ 1e-6**, identical to what
every prior champion passed:

- Numerical gate (job 13773687): the f32full staircase calibrated a positive
  tolerance (3.43e-7) with every accepted solve evaluator-certified at
  ≤ 1e-6 (measured ~1.0e-6 at the staircase edge, 27 iterations, repeat
  delta bitwise 0.0). Independently reproduced on macOS/Apple-BLAS in the
  local pre-submit gate.
- End-to-end A/B (job 13773689): all 240 accepted trials (3 arms × 2 orders
  × 40) re-passed the evaluator per trial, with repeat agreement ≤ 1e-8 and
  dagteam-vs-colored solution agreement 2.0e-6 (both orders certified
  at 1e-6, so cross-order agreement at that level is expected).

So yes — the promoted configuration meets the same 1e-6 independent
accuracy standard as the F64 baselines, verified per solve, not assumed
from the internal residual.

## Where the details live (chronological)

- v21 colored A/B: `fgs_r4_followup_evidence_20260914/colored-v21-13738665/`
- v22 chunked A/B (rejected): `.../chunked-v22-13749231/`
- NUMA placement diagnostic: `numa_placement_findings_20260918.md`
- Spec (recommendation + gates): `fgs_acceleration_recommendation_20260918.md`
- Replay gates 2a–2d: `fgs_acceleration_status_20260919{,b,c}.md` +
  `.../replay-gate2{,c,d}-*/`
- Production implementation + gate-1 + plumbing + v23 gates + PROMOTED:
  `fgs_acceleration_status_20260919d.md`
- v23 pins/jobs/decision rules: `fgs_acceleration_provenance_20260919d.md`;
  evidence `.../dagteam-{numgate-13773687,v23-13773689}/`
- Next-agent handoff: `fgs_acceleration_reset_prompt_20260919f.md`
