# 018 GPU — merge-off λ3.0 ladder + 3R ladder results (2026-08-29, context-reset handoff)

Written for the next agent. Read this FIRST; do not re-derive anything in §1–§3.
Companions: `gpu_3r_ladder_status_20260827.md` (cliff root cause + band-aid, all
still true), `gpu_launcher_handoff_20260827.md` (launcher/silo/env, open items),
`phase_17_nprop_nt_ladder.md` (ladder design + the 2026-08-29 findings sections
— READ ITS LAST THREE SECTIONS, they carry the analysis this ladder tests).

## 1. 3R ladder: COMPLETE, 9/9 clean (2026-08-28/29)

All arms finished full step count, gate clean, no cliff recurrence, no restarts,
np max 433k ≪ 1.5M cap. Three NT36 arms were hopped from eng to mgh (cancel +
resubmit from the gh200 silo — eng-pending jobs canNOT be `scontrol update`d to
mgh: wrong arch/silo/gres/script-arg; use the resubmit pattern in §4).

| arm | job | CT̄ (10-rev) | scatter | Phase-2e conv | vs 1.5R prior |
|---|---|---|---|---|---|
| NT36 λ2.0 | 13507288 | 0.0706819 | ±0.059% | yes | — |
| NT36 λ3.0 | 13507289 | 0.0707752 | ±0.076% | yes | −0.13% ✓ |
| NT36 λ4.0 | 13507754 | 0.0709842 | ±0.069% | yes | −0.074% ✓ |
| NT72 λ2.0 | 13504010 | 0.0710983 | ±0.115% | yes | — |
| NT72 λ3.0 | 13504011 | 0.0718444 | ±0.109% | no | −0.16% ✓ |
| NT72 λ4.0 | 13504012 | 0.0725274 | ±0.10% | yes | — |
| NT144 λ2.0 | 13504007 | 0.0731787 | ±0.522% | no | — |
| NT144 λ3.0 | 13504008 | 0.073829 | ±1.31% | no | — |
| NT144 λ4.0 | 13504009 | 0.073532 | ±0.36% | no | — |

Headline: **λ converges (flat at fixed NT); NT does not** (+1.5–2%/doubling,
and all NT144 arms fail the Phase-2e scatter tolerance).

## 2. Root-cause work on the NT climb (2026-08-29, in phase_17 doc — summary)

1. **Gross shedding per rev is nearly NT-invariant** (early np slopes 14.7k /
   14.7k / 17.2k; mature young-cohort 20.2k / 21.4k / 24.5k per rev). NT144's
   +17–21% excess is the `max(1, ceil(overlap·dist/σ))` quantization in
   `SigmaOverlap` (`FLOWPanel_wake.jl:734`): the inner 18/41 stations are
   floor-bound at 1 particle/step at NT144. Excess is inboard, low-Γ.
2. **The mature-population difference (229k/325k/433k at λ3.0) is
   age-progressive REMOVAL, ~2× stronger at NT36**: np trajectories identical
   to rev ~10; NT36 peaks (244k @ rev 14) then shrinks while truncation is not
   yet engaged; axial density ratio NT144/NT36 grows 1.3 → 2.6 with wake age.
   Only merging is active mid-wake (MERGE_PARTICLES=true, MERGE_R_FACTOR=0.0055
   identical across arms) → **hypothesis: the NT ladder is secretly a
   merge-error ladder** (NT36 reads low CT̄ because merging eats its old wake).
3. Useful trick: newly shed particles in any VTP = σ exactly in the 41-station
   set (σ_j = 0.313·c_j; core spreading re-stamps all older σ) AND x < 0.02.
4. Estimated no-merge np: ~440k(NT36/72) / ~520k(NT144) at rev 30; saturation
   ~520–690k. Under the 1.5M cap; expect ~1.5–2× step cost at NT144.

## 3. LIVE: merge-off λ3.0 ladder (submitted 2026-08-29 ~23:10)

Env: 3R ladder env + `MERGE_PARTICLES=false`; run names `<case>_3r_nm`.
Full submit lines recoverable: `sacct -j <job> -X -P -o SubmitLine%600`.

| job | arm | where | state at handoff | wall |
|---|---|---|---|---|
| 13513063 | n4_nt144_l3p0 | mgh-1-2 (gh200 silo) | RUNNING from 23:15 | 24 h |
| 13513064 | n2_nt72_l3p0 | eng (h200 silo) | PENDING est. 08-30 07:07 | 14 h |
| 13513065 | l3p0 (NT36) | eng (h200 silo) | PENDING est. 08-30 22:00 | 8 h |

FIRST CHECK: confirm merge is actually off — `MERGE_PARTICLES=false` is in the
sbatch --export (verified via sacct) and the driver wires it as
`every = merge_particles ? 1 : 0` (driver line ~705), but the run banner was
not yet printed at handoff. Verify via the run's `.metadata.toml`
(`merge_particles = false`) or np/rev ≈ 14.7k with no mid-run decline.

Interpretation when done: if CT̄(NT) flattens vs §1's λ3.0 column
(0.0708/0.0718/0.0738), the NT climb was merge error, not integration error.
Also check: np vs the no-merge estimates in §2.4; step time (NT144 may run
~2× slower late); Phase-2e scatter at NT144 (does it calm down?).

## 4. Recipes / quirks (unchanged but restated)

- eng→mgh hop: `scancel <job>`, then from `~/FLOWPanel-018-gpu-gh200`:
  `sbatch --job-name=... --partition=mgh --gres=gpu:gh200:1 --constraint=arm
  --qos=normal --no-requeue --cpus-per-task=72 --mem=192G --time=08:00:00
  --output=logs/slurm/slurm-%x-%j.out --error=... --export=ALL,<same env>,
  P018_RUN_NAME=<name> examples/run_dji9443_hover_ct_gpu.slurm.sh gh200 <case>`
- ssh orc needs `bash -lc` for slurm (or `bash -ls` with a heredoc).
- `@info` silently suppressed in runs — instrument with `println(stderr, …)`.
- Walltime hit (exit 124) = valid partial → chain `_s2` via RESTART_STEP/NAME/PATH.
- Overflow of MAX_PARTICLES: raise it and relaunch (Ryan pre-approved).
- rsync with `--checksum`. sacct is authoritative.
- ParaView endpoint VTPs (λ3.0, first-10 + last-10):
  local `FLOWPanel.jl/data/p018_l3p0_3r_endpoints/{nt36,nt72,nt144}/`.
  Retention: runs keep only last ~1,500 VTPs — true steps 0–9 exist only for
  NT36 arms.

## 5. Pending / guardrails (ask Ryan first)

- Notebook + ledger entries for the WHOLE arc (cliff → band-aid → 3R ladder →
  λ/NT findings → merge-off ladder): approved-pending-draft, still undrafted.
- Commits: band-aid + launcher changes still uncommitted in local FLOWVPM.jl /
  FLOWPanel.jl (mirrored md5-identical to both silos).
- 026 particle-splitting implementation: designed, not implemented, on hold.
- Merge-count instrumentation in FLOWVPM (would confirm §2.2 directly): NOT
  done — touches FLOWVPM source, needs Ryan's go-ahead.
- Do NOT delete 1.5R-era data/p018_* dirs.
- Next session plan per Ryan: brainstorm other approaches to the NT problem.
