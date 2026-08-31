# 018 GPU — 3R ladder relaunch status (2026-08-27, end of session)

Written for the next agent at a context reset. Companions:
`gpu_nt144_cliff_findings_20260827.md` (RESOLVED — see §1),
`gpu_launcher_handoff_20260827.md` (launcher/silo/env; open items 1–7),
`../026_sigma_growth_particle_splitting/particle_splitting_design.md` (real fix design).

## 1. Cliff RESOLVED — root cause + band-aid (do not re-derive)

The NT144 cliff was NOT the counting sort. Root cause, confirmed closed-loop:
a single runaway particle grows sigma via the rVPM compression term
`σ ← σ − dt·σ·Z` (FLOWVPM_timeintegration.jl; 052c guard capped only shrink),
crossed the radix-FMM ell=3 adequacy limit — computed 0.0388975 vs measured
sigma_max 0.0388918 (step 2252, fast) → 0.0389060 (step 2253, slow) — which
made `_radix_sigma_outgrown!` (FLOWVPM_fmm_radix.jl:643) rebuild the cache at
ell=2 (64 cells / 267k particles, quasi-dense near field), 6.5→52 s/step,
×8 = exactly one ell decrement. Verified by: instrumented silo run (printed
ell=2; NOTE `@info` is silently suppressed in these runs — use
`println(stderr,…)`), VTP sigma extraction, and the fix: clamping 7 particles
restored 7.8–9.2 s/step at the same steps. Diagnostic instrumentation was
REMOVED from the gh200 silo afterwards (`.pre_sortdiag` backup files remain
in ~/FastMultipole-018-gpu-gh200/src/ — can be deleted).

Band-aid deployed (all mirrored md5-identical to both silos + orc main tree):
- `FLOWVPM.jl/src/FLOWVPM_timeintegration.jl` — sigma_guard `:ceil` key;
  CPU + GPU-broadcast paths `clamp(new_sig, sfloor, sceil)`.
- `FLOWPanel.jl/src/FLOWPanel_warmstart.jl` — restored particles also clamped
  at load (env `SIGMA_CEIL`); REQUIRED because the FMM cache is built at the
  first evaluation and nothing ever re-deepens a shallow cache.
- `FLOWPanel.jl/examples/rotor_hover_pressure_comparison.jl` — env knobs:
  `SIGMA_CEIL` (m, default Inf), `TRUNCATION_RADIUS_R` (default 1.5),
  `MAX_PARTICLES` (default 500000).
- `FLOWPanel.jl/examples/run_dji9443_hover_ct_hpc.slurm.sh` — three new
  `l2p0` case blocks (λ=2.0 replaces λ=2.4 in the ladder, Ryan 2026-08-27).

## 2. Live 9-arm Phase-17 ladder (submitted 2026-08-27 ~13:30)

All: `SIGMA_CEIL=0.030, TRUNCATION_RADIUS_R=3.0, MAX_PARTICLES=1500000,
P018_SETTLE_REVS=22`; run names `<case>_3r`; truncation DEPTH unchanged (4R).
Fresh from step 0. 1.5R-era data dirs preserved (do not delete).

| job | case (λ, NT) | where | state at handoff | wall |
|---|---|---|---|---|
| 13504007 | n5_nt144_l2p0 (2.0, 144) | mgh-1-1 | RUNNING, 2.2 s/step @135, banner radius=3.0R λ2.0 ✓ | 24 h |
| 13504008 | n4_nt144_l3p0 | mgh | PENDING (2nd GH200 node busy) | 20 h |
| 13504009 | n4_nt144_l4p0 | eng-1-1 | RUNNING, banner 3.0R ✓ | 24 h |
| 13504010–12 | n3_nt72_l2p0 / n2_nt72_l3p0 / n2_nt72_l4p0 | eng | PENDING | 14 h |
| 13504013–15 | n2_l2p0 / l3p0 / l4p0 (NT36) | eng | PENDING | 6 h |

Walltimes deliberately tight for backfill; the launcher's internal timeout
(wall − 10 min) + gate makes a timeout (exit 124) a restartable partial —
chain `_s2` restarts via `RESTART_STEP/NAME/PATH` (recipe + prior submit
lines: `sacct -j <job> -X -P -o SubmitLine%600`).

Cancelled this session (all superseded): 13502383/84 (NT72 1.5R), 13502385/86
+ 13503612/13 (NT144 pre-3R), 13503937 (2R run, data deleted), 13502567 (H200
arch test, moot), 13503398/13502670/13502720 (diagnostics), 13490700 (Phase-17
**CPU** NT144 λ2.4, obsolete; CPU NT72 sibling died NODE_FAIL earlier).

## 3. Watch items for the running ladder

1. np vs the 1.5M cap as NT144 wakes mature at 3R (wake-health monitor logs
   np per step; overflow errors loudly at add_particle — remedy: raise
   MAX_PARTICLES, relaunch).
2. Cliff must not recur: sigma capped at 0.030 < ell=3 limit 0.0389 (limits
   rise with box size L, so headroom only grows). If s/step jumps ~8×,
   suspect the same geometry mechanism — instrument with println, not @info.
3. Device memory at MAX_PARTICLES=1.5M: FMM cache ~3× the ~16 GB used at
   500k — fine on 96 GB, but unproven at maturity.
4. CT̄ per rev vs prior arms (±0.36 % scatter band; gate + CT tables in run
   logs); λ2.0 arms have no prior reference.

## 4. Pending / not done

- Ledger + notebook entries for the whole arc (cliff diagnosis → 026 design →
  band-aid → 3R/λ2.0 relaunch): NOT written, needs Ryan's approval.
- 026 particle-splitting implementation (design doc ready; band-aid to be
  removed when it lands).
- Launcher-handoff open items 1–7 remain (incl. metadata-at-start fix, gate
  proof, silo resync discipline — silos DO now match local for the files in
  §1, but 052c-era src drift noted in translate_batched_cuda.jl rect-panel
  section is still local-only).
- All of §1's changes are UNCOMMITTED in local FLOWVPM.jl / FLOWPanel.jl
  (plus the pre-existing uncommitted launcher work); commit needs Ryan.
