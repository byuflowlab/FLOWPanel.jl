# 026 de-risk trio harvest (2026-09-15)

Runs 13691080/81/82 (tag `campaign/p026-derisk-20260914`), all COMPLETED,
gate_rc=0, gpu_gemv=all / cpu_gemv=0 / nan_lines=0. Run dirs MOVED to
`orc:~/projects/FLOWPanel.jl/data/` with symlinks left in the campaign
worktree per provenance data policy. Banners verified in all three arms:
f_visc=0.587 f_comp=0.73 f_elong=0.3, SIGMA_FLOOR_FRAC=0.1.

## Headline discrepancy: NREVS env was NOT honored

All three drivers ran **nrevs=8 physical + 1 spinup = 9 revs**, ignoring
the submitted `NREVS=20` (cap030) / `NREVS=12` (exp pair). cap030 banner:
NT=144 → 1296 steps = rev 9.0 exactly; exp arms: 36 steps/rev → 323 steps
= rev ~9.0. (The 09-15 reset prompt's "≈64.75 steps/rev" inference was
wrong — steps/rev matched the case defs; it's the rev count that didn't
take.) Root cause not yet traced (env plumbing vs case-def override) —
must be fixed/understood before the extension and before wave 2.

## Acceptance scoring (block in p026_derisk_20260914_provenance.md)

| criterion | verdict | evidence |
|---|---|---|
| completion or guard-death w/ interpretable telemetry | **PASS** | 3/3 COMPLETED, rc=0, clean gates |
| cap030 no cliff + 6–8 s/step class | **PASS within window / coverage INCOMPLETE** | mean 7.54 s/step (median 7.29, late rise to ~12 at 212k particles); NO adequacy stress observed — but run ended at rev 9.0 and historical cliff sits at rev ≈15.3 (step ~2200 @ NT144), so the cliff window was never reached |
| split/skip counters sane | **PASS** | zero capacity/mech_disabled skips, no NaN, monotonic sane growth (cap030: 12,575 elongate / 108 compress / **0 viscous**; particles 0→211,882) |
| exp pair separates merge policies | **PASS** (with caveat) | measurable separation, see below; but NO per-event merge telemetry exists in logs/CSVs — churn attribution is indirect |
| floor 0.1 doesn't destabilize healthy phase | **PASS** (no counter-evidence) | early phase smooth in all arms (CT smooth, 0 guard trips, 0–2 splits pre-step-107); floor-clamp events not instrumented per-step, so engagement count unmeasurable |

**Overall: PASS** → per Ryan's 09-14 ruling the remaining 11 arms keep
SIGMA_FLOOR_FRAC=0.1. Wave-2 submission itself remains Ryan-gated.

## cap030 (13691080)

- NT=144, 9 revs, 1296 steps, 7.62 s/step wall-average, H200.
- No FMM-adequacy instrumentation in this run's output (adequacy/sigma_max
  greps: 0 hits; wake-health `min_sigma_ratio` all-NaN under this config).
  `Current sigma growth` 1.00→1.06 vs Critical 1e9 (unrelated metric).
- max_gamma_over_sigma2: ~13 (s143) → peak 50.6 (s1007) → 35–47 at end.
- CT_bernoulli: 0.0765/0.0720/0.0790/0.0750 at 25/50/75/100%; rev8 mean
  0.07756 (ptp 0.00582), rev9 0.07522 (ptp 0.00533), convergence window
  flagged but not converged.
- **viscous=0 across all 901 split-event lines** despite viscous=true.
  Not inert globally — exp arms fire viscous late (424/1,995 events) — so
  reading: threshold not crossed by rev 9 in the cap030 regime. Watch item.

## Merge A/B: exp_split (A) vs exp_split_mo35 (B), 323 steps each

- **2.7× retention drop does NOT reproduce at campaign σ≈0.0381R** — as
  the Task-1 regime caveat predicted (smoke B at σ=0.26R doesn't transfer).
  B/A particle count: 0.95–0.92 through step 287, then crossover, final
  **1.122** (A 443,256 vs B 497,137).
- Split totals A vs B: viscous 424 vs **1,995** (4.7×), compress 268 vs
  187, elongate 117,048 vs **161,129** (1.38×), all concentrated late.
  Indirectly consistent with merge-churn freeing crowding room; direct
  merge telemetry ABSENT (no per-step merge counts or merged-pair σ stats
  logged — instrumentation gap to fix before wave-2 merge A/Bs if churn
  attribution is wanted).
- No σ-pump signature: B's min_sigma / p1_sigma_ratio trend LOWER than A
  late (0.103 vs 0.112; 0.136 vs 0.245 at end). No mean/max σ column
  exists, so mid/upper-percentile pump not ruled out.
- CT: differences within noise; final-2-rev cycle means A 0.075856 ±2.41%
  vs B 0.076393 ±0.77% (B tighter). Neither converged per Phase-2e
  criterion.
- Cost: merge gate free — median 5.65 vs 5.67 s/step (+1.6% total wall).
- Gotcha: `MERGE_OVERLAP=3.5` never appears in B's stdout banners — only
  in `*_case_metadata.toml` (`merge_overlap = 3.5`; A shows `nan`). Don't
  gate on log greps for it.

## Proposed cap030 extension (Ryan-gated)

Cliff question is unanswered. Propose a chained restart from the retained
final state (step 1295/rev 9) through rev ≈18 (cliff rev 15.3 + ~2.7 rev
margin): +~1300 steps at 10–14 s/step (count still growing) ≈ **4–6 h
H200**, same tag/worktree (code unchanged). Prereq: trace why NREVS was
ignored and set the rev target by whatever knob the driver actually honors;
verify RESTART_STEP path against the retained checkpoint. Optionally add
FMM-adequacy logging first (none exists in this driver's output), else the
cliff must again be inferred from s/step + γ/σ² trends.

## Telemetry gaps surfaced (for wave-2 hardening, cheap)

1. No per-event merge logging (counts, merged-pair σ). 2. No σ mean/max in
wake-health CSV. 3. No σ-floor clamp event counter. 4. No FMM-adequacy
metric in stdout. 5. NREVS env silently ignored.

## Hygiene facts (Task B)

- Smoke run dirs live PHYSICALLY inside `orc:~/wt026gpu/FLOWPanel.jl/data/`
  (NOT the shared root, contra the 09-15 reset prompt): scr_p026cpuv_split
  1.3G, scr_p026gpuv_split 1.3G, scr_p026gpuv_splitmerge 768M, logs 593K.
  Removing wt026gpu deletes them → Ryan decides: relocate to shared root
  first, archive, or discard (they were verification-only; parity results
  are recorded in the provenance file + commits).
- `data/p026_restart_gpu40_s950` (3.2M remnant) left for the archiver.
- A/B harvest intermediates (logs+CSV copies) in local scratchpad
  `.../scratchpad/p026ab/` — disposable.

Harvest method: two harvester subagents over ssh (read-only), ad-hoc
grep/python; sources = slurm-fp-p026dr-*.{out,err} + monitor02/04 +
CT_vs_rev/CT_per_rev CSVs + case_metadata.toml.
