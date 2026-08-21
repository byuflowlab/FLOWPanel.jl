# Pause manifest — 2026-08-15 23:50 MDT (Ryan-directed weekend pause)

All rander39 jobs on orc scancel'd at ~23:52 MDT 2026-08-15 (queue
verified EMPTY). `scontrol suspend` is permission-denied for regular
users (tested job 13184024); kills were Ryan-authorized. **These jobs
did NOT fail — do not post-mortem them as ignitions.** Resume scheduled
00:00 MDT Monday 2026-08-17 (in-session cron 4a0fcf9d; if that died, a
fresh agent executes this manifest manually).

Restart steps are post-kill, verified present in ALL three output
families (body1, wake1_particles, wake1_filaments) AND registered in
body1.pvd. Monitors + CT CSVs harvested to local `data/<run>/` at pause.

## 018 runs — warm-chain resume (ops_reference "Restart chaining")

Template (from cluster repo top level):

```
sbatch --job-name=fp-018-<tag>_rs1 --time=<WALL> --mem=64G \
  --export=ALL,P018_RUN_NAME=<run>_rs1,P018_SETTLE_REVS=22,RESTART_STEP=<STEP>,RESTART_NAME=<run>,RESTART_PATH=data/<run> \
  examples/run_dji9443_hover_ct_hpc.slurm.sh <case_tag>
```

| run (= case_tag) | killed job | RESTART_STEP | rev at kill | target | WALL | notes |
| --- | --- | --- | --- | --- | --- | --- |
| `p018_cs_l2p4` | 13178762 | **834** | 23.2 | 30 revs (1080 steps) | 24:00:00 | Phase-16 ladder bottom |
| `p018_cs_l3p4` | 13178763 | **778** | 21.6 | 30 revs | 24:00:00 | ladder mid / P3 A/B arm |
| `p018_cs_l4p8` | 13178764 | **713** | 19.8 | 30 revs | 36:00:00 | ladder top; gos2 was SUBSIDING at kill (6707@677 → 3468@713) — hot phase passed, not ignition |
| `p018_upin_nt72_rlxf0p3` | 13183888 | **252** | 3.5 | 30 revs (2160 steps) | 72:00:00 | model-def arm; will TIMEOUT short of 30 like all NT=72 |
| `p018_upin_nt72_rlxf0p16334` | 13183889 | **194** | 2.7 | 30 revs | 72:00:00 | model-def arm |
| `p018_upin_nt144_rlxf0p3` | 13183998 | **268** | 1.86 | 72 h segment (orig plan ~rev 14.5) | 72:00:00 | model-def arm; **ASK RYAN before any chain beyond this segment** |

All six: SETTLE_REVS=22 (total-rev target unchanged — the driver runs to
the step count implied by SETTLE regardless of restart point). Banner-
verify each within minutes of start (Phase-16 arms: sigma_chord:0.313
sigma_floor:0 das_lambda:2.4/3.4/4.8; upin arms: NT:72|144 pps:6|3
rlxf:0.3|0.16334 das_uniform:3.4; all: overlap:2.75 merge_r:0.0055
nwakerows:1 visc:true settle:22, mesh 45_185_ct4) AND verify the
`simulate_warmstart!: resuming from step <STEP>` line. scancel on any
mismatch. Launcher md5 at pause: 9a5b2f2cb459d75a1b5986be56493f02
(local==cluster).

## Non-018 job (other session's work)

| job | name | script | resume |
| --- | --- | --- | --- |
| 13184015 | `fm041aH` | saved to `~/fm041aH_batch_13184015.sh.pausecopy` (also /tmp copy) | resubmit from WorkDir `/home/rander39`: `cd ~ && sbatch fm041aH_batch_13184015.sh.pausecopy` — 12 h wall; restarts FROM SCRATCH (checkpointing unknown; it belongs to the FastMultipole/021 line — if that session is active Monday, defer to it) |

## Progress/health at kill (for continuity of the watch)

- Ladder wake-health at kill (step, n_p, max_u, min_sigma, gos2):
  l2p4 833: 183752, 23.0, 8.41e-4, 232; l3p4 777: 197415, 15.1,
  7.86e-4, 441; l4p8 713: 197048, 18.2, 5.47e-4, 3468 (declining).
  All healthy; sigma-ratio columns NaN under station-σ law — judge
  absolutes (calibration in phase_16 §Log).
- upin arms all early/healthy (revs 1.9–3.5, gos2 3–6).
- After resume, re-arm the state+tripwire monitor on the six `_rs1`
  jobs (absolutes: max_u>500 or min_sigma<4e-4).
