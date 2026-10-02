# Reset prompt — NEW item 032 (omit shed locations) + 026 slate watch (2026-09-18c)

You are picking up work in
`/Users/ryan/Dropbox/research/projects/FLOWPanel.jl` (branch
**fastmultipole**) + sibling `/Users/ryan/Dropbox/research/projects/FLOWVPM.jl`
(branch **flowpanel**). Read `CLAUDE.md` and the policies it names
(`agent_policies/WORKFLOW.md` before editing source; `HPC.md` before
ORC). This file supersedes `rerunslate_reset_prompt_20260918b.md` (its
state section is still accurate background; its watch tasks are folded
in below). Durable 026 record: `rerunslate_provenance_20260918.md`.

## PRIMARY TASK — new BRAINSTORM item 032: omit shed locations from a panel object

Ryan's directive (2026-09-18, verbatim intent): *make a new item to add
a feature to omit shed locations from a panel object, and then try doing
that at the root-most shed locations. I suspect the strong root particle
strength is causing the instability here, and it isn't even strictly
physical.*

Motivating evidence (026 §22.3 rerun slate, in the provenance file):

- A1 (`scr_p026s9r2_explg_fs`) died at step 328/467 via the euler_exp
  substep guard (max-over-particles dt·|L| > 2048). Offending particle
  LOCALIZED from the step-327 VTP: **depth x = −0.37R (above the rotor
  plane, thrust side), radial r = 0.344R** — the fountain-flow
  recirculation region. Driver/victim structure: a σ-at-floor particle
  carrying |Γ|=5.4e-3 (Γ/σ² = 9.1e4, field max) imposes the gradient on
  near-zero-Γ neighbors 1.5σ away. Global field healthy (Ryan verified
  in ParaView: fountain surge near root, "not fully unstable (yet)").
- A2 (floor-only) died at step 274 — the EXACT step of its wave-2 twin
  ⇒ this ignition channel is independent of the (now fixed) merge
  σ-pump. Hypothesis: root-shed circulation feeding the fountain region
  is the source, and root shedding at the blade root is itself of
  dubious physicality (root cutout / hub interference in reality).

What to do:

1. Create `BRAINSTORM/032_omit_shed_locations.md` (+ add to
   `BRAINSTORM/INDEX.md`) laying out: motivation (above), design
   options, and a validation plan. Keep 026 cross-references.
2. Design + implement the feature: a user-facing way to EXCLUDE
   selected trailing-edge/shedding locations when building a
   `RigidWakeBody`'s shedding (e.g. a predicate or index mask applied
   at `calc_shedding_from_seed` level or a post-filter on the shedding
   matrix before body construction). Constraints:
   - **CRITICAL invariant** (CLAUDE.md): shedding must be computed from
     the CONSTRUCTED body's cells (`ensure_winding=true` re-winds in
     place) — build `noshedding` body, run `calc_shedding_from_seed` on
     ITS nodes/cells, filter, rebuild. The omission feature must not
     break this flow, and omitted edges must not silently reappear.
   - Bound-circulation bookkeeping: omitting a shed edge changes where
     the bound-vortex return path closes; check the Kutta/jump closure
     and `BoundCirculationMonitor` still make sense (flag, don't hide,
     any conservation consequence in the item file).
   - Driver plumbing: env knob in
     `examples/rotor_hover_pressure_comparison.jl` +
     `examples/run_p018_screen_hpc.slurm.sh` (e.g.
     `SHED_OMIT_ROOT_FRACTION=<r/R below which TE shedding is omitted>`
     or an explicit index list) so it's usable in screen cases; print it
     in the banner and metadata TOML.
   - Unit tests (shedding count/edges with and without omission; wing
     regression unchanged when knob off) per `TESTING.md` matrix.
3. Local smoke, then A/B on HPC (Ryan-gated submission): rerun
   `scr_p026s9_explg_fs`-class case with root-most shed locations
   omitted (start with the innermost station(s) ~ r < 0.15–0.2R; the
   DJI blade root geometry defines what "root-most" means — inspect the
   mesh/TE seed). Success signal: the fountain-region Γ concentration
   and the ~274–328-step guard trips disappear (or move out), CT
   changes stay small and explainable. Campaign rules apply (worktrees,
   tags, provenance) if it graduates beyond a smoke.

## SECONDARY — 026 slate watch (carried from 20260918b)

Slate scoreboard (jobs on ORC; judge by outputs, sacct unreliable;
tracebacks are in `logs/slurm/slurm-fp-052-scr-gpu-<jobid>.err`, the
.out ends at the GATE line):

| arm | job | state at handoff |
|---|---|---|
| A1 explg_fs | 13758586 (m13h) | DIED 328 — localized Γ-ignition (see provenance; VTK steps 278–327 local at `~/scr_p026s9r2_explg_fs_last50steps/`) |
| A2 explg_floor | 13763819 (eng) | DIED 274 (dt·\|L\|=5251.8) — exact wave-2 twin step |
| A3 ctrllg_fs | 13763820 (eng) | RUNNING (banner PASS 9/9); ctrl-channel probe, watch ~260–290 window |
| A4 explg_fs_cap3 | 13763821 (eng) | PENDING (eng serializes, 64 CPU/job); k=3 tripwire — banner must show SIGMA_CEIL=0.0071 and clamp=[…,0.0071]; key question: does the σ-cap arrest the ignition? Prediction given the close-pair mechanism: the cap does NOT bind the Γ-driver (its σ is at the FLOOR, not the ceiling) — watch for a same-class death |
| A5 explg_split | 13763822 (eng) | PENDING; wave-2 twin died @280 |

- Re-establish a watch (predecessor's Monitor died with its session):
  poll `ssh orc 'bash -lc "sacct -n -X -j 13763820,13763821,13763822 -o JobID,State"'`,
  anchor parsing on the job-id (MOTD noise). Banner-check A4/A5 on
  start; autopsy terminal arms from the .err; record everything in the
  provenance Banner/Outcomes tables.
- When all terminal: harvest (delegate `harvester` — MOVE run dirs from
  campaign worktree `data/` to `orc:~/projects/FLOWPanel.jl/data/` +
  symlink back), slate verdict, ledger line, launch `hpc-storage`
  (A1's worktree run dir is 11 G). Note 032's A/B may want to reuse
  these dirs for comparison — coordinate.

## Owed / parked (carried)

- Notebook entry for the 026 arc (Ryan "not yet" ×2) — offer only.
- orc branch-divergence ruling; github pushes of both branches+tags —
  Ryan-gated.
- scr_p026gpuv_split archive retry (hpc-storage) once quiet ≥24 h.
- NT144 cap-ladder rung parked until slate verdict.
- 021 silo cleanup (021 package §7).

## Ground rules

Local ≤4 threads. `ssh orc` needs live ControlMaster socket (`! ssh orc
echo ok` if 2FA blocks); `bash -lc` for Slurm commands; MOTD
contaminates output. Read `BYU_ORC_AGENTS.md` before ORC execution.
Commits, HPC submissions, and notebook writes are Ryan-gated. Campaign
pins: tag `campaign/p026-rerunslate-20260918` (FLOWPanel 6e52628,
FLOWVPM 8d4a3b4, FMM ac7230a6), worktrees at
`orc:~/campaigns/p026-derisk-20260914/`.
