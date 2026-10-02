# FGS acceleration status — gates 2c (precision) + 2d (schedule) harvested (2026-09-19c)

Continues `fgs_acceleration_status_20260919b.md` (rev-b verdict) per
`fgs_acceleration_reset_prompt_20260919c.md`. Spec =
`fgs_acceleration_recommendation_20260918.md`. Provenance of the two jobs =
`fgs_acceleration_provenance_20260919c.md`.

## Run record

| Job | Gate | State | Elapsed | Node | End (2026-09-19) |
|---|---|---|---|---|---|
| 13773580 | 2c full-F32 precision ladder | COMPLETED 0:0 | 3:37 | m12-2-13 | 10:39:33 |
| 13773581 | 2d dagteam split executor | COMPLETED 0:0 | 3:35 | m12-2-13 | 10:43:08 |

Both `.err` empty; both logs end `== done ==`. (An earlier report that the
jobs "failed" was incorrect — judged by outputs per policy, both are clean.)
Evidence harvested to
`fgs_r4_followup_evidence_20260914/replay-gate2c-13773580/` and
`replay-gate2d-13773581/` (slurm logs + full results dirs). Numastat
snapshots captured non-empty this time (poll-and-keep-last fix verified).
Worktree/commit: `/home/rander39/wt-p021-fgs-gate2` @ FastMultipole
`e904e763`, CENSUS/EDGES overrides md5-verified (see provenance file).

Model (unchanged): T ≈ 3.0 s remainder + 81-sweep stream; baseline colored
@ j16 = 10.116 s; minimum ≤ 6.744 s (stream ≤ 3.744 s ⇔ 61.9 GB/s
F64-equiv); design target ≤ 5.058 s (stream ≤ 2.058 s ⇔ 112.7 GB/s).
All BW below is useful F64-equiv from min-of-12-sweeps.

## Gate 2c — precision ladder (rowpar unless noted)

Serial (node-0 bind, t1):

| Precision | s/sweep (min) | BW GB/s | 81-sweep stream |
|---|---|---|---|
| F64 | 0.1029 | 27.81 | 8.337 s |
| F32conv | 0.0843 | 33.98 | 6.825 s |
| **F32full** | **0.0548** | **52.22** | **4.441 s** |

Rowpar ladder:

| Config | BW GB/s | Stream s | T ≈ | × |
|---|---|---|---|---|
| F32full owner t4 | 76.20 | 3.043 | 6.04 | 1.67 |
| F32full owner t8 | 65.11 | 3.562 | 6.56 | 1.54 |
| **F32full owner t16** | **90.21** | **2.571** | **5.57** | **1.82** |
| F32full owner t32 | 83.49 | 2.777 | 5.78 | 1.75 |
| F32full interleave t16 | 88.37 | 2.624 | 5.62 | 1.80 |
| F32full interleave t32 | 84.37 | 2.749 | 5.75 | 1.76 |
| F32conv owner t16 | 67.54 | 3.433 | 6.43 | 1.57 |
| F32conv interleave t16 | 78.05 | 2.971 | 5.97 | 1.69 |
| F32conv interleave t32 | 72.39 | 3.203 | 6.20 | 1.63 |

Handoff kill switch: PASS again — 1068 handoffs/sweep, avg 3.132 µs
(budget < 5.8 µs).

**Precision verdict: F32full wins at every matched placement** (owner t16
90.2 vs 67.5; interleave t16 88.4 vs 78.1). Rev-b's anomaly is closed: the
convert kernel was the limiter; with F32 storage end-to-end the serial
stream alone moves 52.2 GB/s (vs F64's 27.8, ≈ the ideal bytes-halving).
Within rowpar, F32full alone does NOT reach the 112.7 GB/s design-target
rate (best 90.2).

## Gate 2d — dagteam split executor (interleave 0-3 unless noted)

Structural invariants matched the M2 smoke exactly: 1068 lower tasks,
3 roots, 48,167 edges, 1065 backward tasks, lower/upper bytes
1,513,294,744 / 1,349,555,288.

| Config | BW GB/s | Stream s | T ≈ | × |
|---|---|---|---|---|
| dagteam F64 t4 | 65.03 | 3.566 | 6.57 | 1.54 |
| dagteam F64 t8 | 69.89 | 3.318 | 6.32 | 1.60 |
| dagteam F64 t16 | 69.47 | 3.338 | 6.34 | 1.60 |
| dagteam F64 t32 | 60.27 | 3.847 | 6.85 | 1.48 |
| dagteam F64 t16 (sock membind ctrl) | 36.78 | 6.305 | 9.31 | 1.09 |
| dagteam F32conv t16 | 114.28 | 2.029 | 5.03 | 2.01 |
| **dagteam F32full t16** | **128.15** | **1.810** | **4.81** | **2.10** |
| dagteam F32conv t32 | 94.92 | 2.443 | 5.44 | 1.86 |
| dagteam F32full t32 | 114.35 | 2.028 | 5.03 | 2.01 |
| rowpar F64 t16 (in-job ref) | 50.81 | 4.564 | 7.56 | 1.34 |
| rowpar F32conv t16 (in-job ref) | 77.56 | 2.990 | 5.99 | 1.69 |
| rowpar F32full t16 (in-job ref) | 83.49 | 2.778 | 5.78 | 1.75 |

**Schedule verdict: dagteam wins best-vs-best and at every matched
precision/team/placement**, with margin:

- Matched t16/interleave: F64 69.5 vs 50.8 (+37%), F32conv 114.3 vs 77.6
  (+47%), F32full 128.2 vs 83.5 (+53%).
- Total-time terms at the champion: T ≈ 4.81 vs 5.78 s (−17%); dagteam is
  the only schedule that reaches the 2× design target, and even pure-F64
  dagteam clears the 1.5× minimum (rowpar F64 never did).
- The priced dag sim (1.76 s / ~132 GB/s), which rev b refused to promote
  on simulation alone, is now confirmed by real measurement (1.810 s /
  128.2 GB/s, within 3%).

Sanity: speedups vs matched-precision serial stay under the 2.849
work/span bound everywhere (max 2.45× at F32full t16). Socket-membind
control collapses to 36.8 GB/s → interleave placement is load-bearing.
Best team size is t16 (one thread per node-0..3 CCX group ×4); t32 is
uniformly worse.

## Picks and projection

- **Precision pick: F32full** (storage + state/accumulate/LU), fallback
  ladder unchanged: F32full → F32conv (certified rev-b GO, T ≈ 5.84 s) →
  selective per-block F64. **Adoption remains CONDITIONAL on the
  independent evaluator at BC rel-L2 ≤ 1e-6** — F32full accuracy has NOT
  yet been certified (only one decade of headroom above F32 eps 1.2e-7;
  the internal residual can flatter the rounded operator over 81 sweeps).
  The gate is not relaxed.
- **Schedule pick: PROMOTE the split (dagteam)** — the win is real,
  measured, and large enough to justify the rebuild under the spec rule
  (+53% stream BW over rowpar at the champion configuration; only
  schedule reaching 2×).
- **Champion configuration: dagteam + F32full + interleave 0-3 + t16 →
  projected T ≈ 3.0 + 1.81 = 4.81 s = 2.10× vs baseline 10.116 s** —
  beats the 2× design target (median-based projection 4.97 s, still ≥2×).
  If F32full fails the evaluator: dagteam + F32conv t16 → T ≈ 5.03 s =
  2.01× (still at target); dagteam F64 t8 → 6.32 s = 1.60× (above
  minimum) as the floor.

## Go/no-go

**GO for production implementation** (TASK 2 of the 20260919b prompt) on
the dagteam + F32full champion, with the accuracy gate as the first
production hurdle: gate-1 correctness harness vs the real implementation →
independent evaluator at 1e-6 (decides F32full vs fallback rung) →
end-to-end interleaved A/B vs unchanged champion at ≥1.5× accepted
throughput, full campaign ceremony. Two-system `FastGaussSeidel`
construction BoundsError at `c18e4b46` must be fixed or declared before
multi-rotor use.

**Awaiting Ryan's go before touching the production solver path.**
