# Phase 6 — Multi-rotor GPU acceptance (1/2/4 rotors, OGE + IGE, 10 revs ≤ 2 h)

## Authoritative acceptance contract — 2026-08-27

This section supersedes the older staging text below.

- Production formulation: `VelocityThroughSources`, with frozen Dirichlet
  sources, scalar-potential body coupling for Dirichlet targets, velocity
  coupling for Neumann targets, and hard block-GS convergence. Hybrid is an
  independent experimental comparison owned by FastMultipole item 052e and
  does not block Phase 6.
- Six independent cases: 1/2/4 rotors x OGE/IGE, fine `45_185_ct4` mesh,
  Gaussian filament regularization, IGE h/R=1.5 and truncation depth 4.5R.
- Run length: 54 spinup steps plus 360 acceptance steps at NT=36, for 414
  total steps. The authoritative total case elapsed time must be no more than
  7200 seconds.
- Performance promotion gate: early/mid/mature particle-population probes and
  a nonlinear fit against particles, panels, and GS sweeps must project no
  more than 6480 seconds with at least 20% device-memory reserve. A one-rev
  linear extrapolation is not admissible.
- Every block solve must remain finite and converge within 50 sweeps at
  normalized tolerance 1e-8. Residual scales are ΩR² for rotor potential and
  ΩR for ground velocity. Required CUDA routes must be nonzero, with zero
  prohibited or unclassified fallbacks.
- Acceptance is six separate submissions after rectangular CUDA parity,
  reduced FLOWPanel route parity, and short 4-rotor OGE/IGE shedding smokes.
  The carrier defaults to smoke and requires explicit acceptance confirmation.

**Historical objective (superseded above):** Run the 022 driver (`examples/rotor_hover_ground_effect.jl`,
`NROTORS` ∈ {1,2,4}) for **10 revolutions per case** (360 steps at NT = 36)
with the **major cost steps on GPU** (wake FMM, panel passes, solve), each
case finishing in **≤ 2 h GPU walltime** (per-case target, Ryan 2026-08-25;
"if possible" — see Decision). GPU support comes from the FastMultipole
extension item **052b**
(`../../../FastMultipole/MATRIX_OPERATOR_REFACTOR/052b-impl-multirotor-ige-gpu.md`),
which extends 052 (single-rotor H200 pipeline) and 052a (GH200 verdict:
discrete-memory path). The existing CPU `p022m_*` arms
(`scripts/p022m_submit_*.sh`, launcher
`examples/run_rotor_multi_ground_effect_hpc.slurm.sh`) serve as the
CPU-parity reference.

**Blocked on:**
1. **052b technical completion** (FastMultipole side; see its Phase gates).
2. **022 Phase-3 particle-policy verdict** — the `p022m` launcher carries a
   placeholder `GROUND_DAMP_BAND_R=0`; IGE arms must not be submitted until
   the production `GROUND_DAMP_BAND_R` / `GROUND_PARTICLE_POLICY` values are
   ruled.
3. **Ruling 7 (single particle field, Ryan 2026-08-25):** all rotors shed
   into ONE shared pfield. The driver currently builds one
   `PanelParticleWake` per rotor
   (`examples/rotor_hover_ground_effect.jl:815`; wakes tuple at
   :1258/:1262) — must be updated first (expected minimal change; one pfield
   referenced by all rotor wakes, or one merged wake). NROTORS=1 must remain
   bit-identical to the legacy path.

## Cases

Six acceptance arms; fine mesh 45_185_ct4 only (coarse 56_57 banned for IGE
— blew up at rev 17.6, RESET BRIEF (c)). Carrier knobs = 022 fixed operating
point (Ω = 628.32 rad/s, σ = 0.0400 R, OVERLAP = 2.75, P_PER_STEP = 12,
d/σ = 3.4, `FLOWPANEL_FILAMENT_REG=vatistas` pinned per the 025 hazard).
IGE arms: `IGE_H_R=1.0`, `IGE_TRUNC_DEPTH=4.0` (h/R + 3), ground disc per
the fixed operating point; damping knobs from the Phase-3 verdict.
Multi-rotor: `ROTOR_SPACING_R=2.7` (driver default), `ROTOR_DIRECTIONS`
default.

| tag | mesh | knobs | job | status |
|---|---|---|---|---|
| p022g_1r_oge | 45_185_ct4 | NROTORS=1, GROUND_ENABLE=false, NREVS=10, GPU | — | staged |
| p022g_2r_oge | 45_185_ct4 | NROTORS=2, GROUND_ENABLE=false, NREVS=10, GPU | — | staged |
| p022g_4r_oge | 45_185_ct4 | NROTORS=4, GROUND_ENABLE=false, NREVS=10, GPU | — | staged |
| p022g_1r_ige | 45_185_ct4 | NROTORS=1, GROUND_ENABLE=true, IGE_H_R=1.0, NREVS=10, GPU | — | staged |
| p022g_2r_ige | 45_185_ct4 | NROTORS=2, GROUND_ENABLE=true, IGE_H_R=1.0, NREVS=10, GPU | — | staged |
| p022g_4r_ige | 45_185_ct4 | NROTORS=4, GROUND_ENABLE=true, IGE_H_R=1.0, NREVS=10, GPU | — | staged |

Staging per the 052 stage convention, combined into one sbatch chain per the
FastMultipole cluster-jobs guidance (H200 queue waits are long):
- **Probe stage:** 1 rev per shape with per-pass GPU timers (wake FMM,
  panel pass 1/2, solve) → linear extrapolation to 10 revs vs the 2 h
  budget. The **4-rotor IGE** probe is the binding constraint
  (≤ 20 s/step average required; single-rotor 052 prediction is
  3.1 s/step, and the rectangular particles→panels pass scales
  ~panels×particles).
- **Acceptance stage:** the 10-rev runs, gated on the probe extrapolation.

## Decision (pre-registered)

A case PASSES iff:
1. 10 revolutions complete in **≤ 2 h GPU walltime**. If the probe
   extrapolates past 2 h for a case (expected risk: p022g_4r_ige), do NOT
   silently relax — report measured s/step, the dominant pass, and the
   scaling bottleneck; the 041k escape-hatch levers (far-field singular
   switch, fast-math, register blocking; 1.6–1.7× F32 on all-pairs) are the
   pre-identified mitigation, applied only as a logged ruling candidate.
2. GPU CT trace matches the CPU `p022m` reference arm for the same shape
   within the parity tolerance defined in 052b.
3. No blow-up through rev 10 (the known rev-17.6 ignition is out of range at
   10 revs, but log CF_x / ground-tangency health regardless).

## Exit criteria

- All six cases pass the pre-registered decision.
- Per-pass GPU timing table recorded (per case: wake FMM, panel passes,
  solve, other; s/step and share).
- IGE/OGE CT ratio reported per rotor count (1, 2, 4) with per-rev scatter,
  in `ledger.md`.
- Offer Ryan a notebook entry at this milestone (ask header + verbosity
  first, per journal rules).

## Task checklist (Ryan-ruled 2026-08-26; plan `~/.claude/plans/work-on-brainstorm-022-multi-rotor-gentle-frog.md`)

Amendments ruled 2026-08-26 superseding the text above: IGE arms at
**h/R = 1.5** (not 1.0; TRUNC_DEPTH = 4.5), filament reg **gaussian**
(018 production carrier, not the vatistas pin), mesh 45_185_ct4 unchanged.
GPU placement ruling: at each submission pick the GPU via the
**slurm-availability skill** — shortest predicted wait among **H200, H100,
GH200, B200** that likely won't preempt; when using H200 prefer
**`qos=eng`** unless the scheduler predicts a shorter wait on `qos=gpu`.
GH200 remains the primary target when waits are comparable.

- [x] **(DONE 2026-08-26, uncommitted)** Shared pfield (Ruling 7): N `PanelParticleWake`s
      reference ONE pfield; identity-dedupe in `_collect_wake_probes`/
      `_collect_wake_sources`, `apply_freestream!`, `propagate!` loop, VTK;
      driver builds rotors 2..N with `pfield=` rotor 1's,
      `max_particles = N*500_000`; damp-band `propagate!` override and
      `ground_diagnostics_monitor` dedupe. Gate: 1r_oge smoke CSVs
      byte-identical pre/post; 2r/4r np conserved, below-ground not N-folded.
- [x] **(DONE 2026-08-26: 13484022–25 RUNNING; 4r_ige ride+rescue ruling)** CPU parity reference: deploy via
      `scripts/deploy_022m.sh` (md5-verify; frozen tree must NOT receive the
      shared-pfield src edits) and submit the four `p022m_{2r,4r}_{oge,ige}`
      arms; verify banners (nrotors, ground, reg:gaussian, FMM knobs),
      first-hour s/step vs walltime; launch hpc-storage (VTK writers live).
- [ ] Commit the 022m scripts + shared-pfield work on `fastmultipole`.
- [x] **(DONE 2026-08-26, uncommitted; setup-only verified)** GPU knob wiring in the driver (pressure_comparison :650-682 pattern:
      `VPM_ARRAYTYPE=cuarray` + `FLOWPANEL_GPU_INFLUENCE=cuda`, radix
      lifecycle, pinned FMM p=4/ncrit=50/theta=0.4 on the GPU branch only;
      loud eligibility errors, `GPU_ALLOW_FALLBACK=1` escape; no GPU-S).
      Setup-only parses green 2026-08-26 for 4r IGE (`BERNOULLI_ONLY=true`
      required by the driver contract; the carrier exports it) and 1r OGE,
      after fixing a load-killing ParseError at driver :2040 (escaped
      quotes inside `$()`).
- [x] **(DONE 2026-08-26, uncommitted; not deployed, no jobs)** GPU launcher `examples/run_rotor_multi_ground_effect_gpu.slurm.sh` +
      `scripts/p022g_submit_*.sh` (arch modes gh200/h200 per the placement
      ruling; probe stage NREVS=2 with `FLOWPANEL_STEP_TIMERS`/
      `FLOWPANEL_GPU_TIMERS`, accept stage NREVS=10; julia 1.11.7 pinned).
      `bash -n` + all six `P022G_SETUP_ONLY=1` arms verified.
- [ ] Deploy to dedicated cluster tree `~/FLOWPanel-022g`
      (`scripts/deploy_022g.sh`, rsync --checksum) + ARM env
      `fp022env-gh200` offline build (052a Phase A recipe) if GH200 is the
      pick.
- [x] Verify 022 Phase-3 damp-band verdict (jobs 13246557/58) → RULED
      2026-08-26: both arms FAILED (PARTICLE OVERFLOW at 500k, ~rev 15–18,
      force blow-up, no CSVs; ppdamp showed 32,202 below-ground particles
      mid-blow-up). Phase-3 verdict void; Ryan applied the breach ruling —
      IGE arms run `GROUND_DAMP_BAND_R=0.1`.
- [ ] Submit six-shape GPU probe chain (placement per ruling); harvest
      per-pass timings; extrapolate vs the 2 h budget.
- [ ] Acceptance runs (NREVS=10) with probe-derived walltimes; 4r_ige over
      20 s/step ⇒ report dominant pass + escape-hatch ruling candidate.
- [ ] Ledger: IGE/OGE CT ratio per rotor count; offer notebook entry.

## Log

- 2026-08-26 (session 3) — **Warm-start blocker RESOLVED; rescue path for
  13484025 is unblocked.** The multi-rotor/ground warm-start NaN was
  root-caused (missing per-step freestream/velocity_te staging before
  `simulate_warmstart!`'s end-of-step replay — geometric, not backend) and
  fixed (`FLOWPanel_warmstart.jl` section 5.0). `runtests_unit_warmstart.jl`
  153/153 incl. the multi-rotor IGE testset (113/113 at 1e-10/1e-11 vs the
  uninterrupted arm); simulate 199/199, solver 412/412 (two test-code fixes
  documented in 052b). Driver restart guards for NROTORS>1/ground stay
  removed (now test-covered). Setup-only parses green for 4r IGE + 1r OGE
  after fixing a load-killing driver ParseError at :2040. Diff-overlap
  audit clean; `git diff --check` clean. All uncommitted; no deploys, no
  jobs.
- 2026-08-26 (session 2, evening) — **Shared pfield (Ruling 7) implemented
  and gate-verified; CPU arms live.**
  - src: `PanelParticleWake(...; pfield=)` shares one field across wakes;
    identity-dedupe in `_collect_wake_probes`/`_collect_wake_sources`,
    `apply_freestream!` (`include_pfield`), `propagate!`
    (`propagate_pfield`), and VTK (`include_pfield` — one particle cloud per
    step, no duplicate ParaView files). Driver: shared allocation
    `N×500_000` (`MAX_PARTICLES` override), damp-band override + ground
    diagnostics dedupe.
  - Gates: 1r_oge **bit-identical** to pre-change code (every CSV and VTP;
    caveat: the wake-health CSV wall_s column and a pre-existing final-step
    C-array run-order sensitivity — a pre-change rerun matches the
    post-change run byte-for-byte, so the discrepancy is baseline
    nondeterminism, not the change). 2r_oge: ΔCT ≤ 0.05%, per-rotor
    symmetry preserved, np 2674→2677 (+0.1%, union-merge pairing), single
    `wake1_particles` output. 4r_oge/2r_ige post-change smokes deferred —
    a parallel session (052b Phase A.1 warm-start) is actively editing the
    tree and shares `data/smoke_*`.
  - CPU arms deployed (cluster src reconciled: remote was an exact 3e5aa08
    snapshot; four files synced to local HEAD) and submitted:
    13484022 2r_oge (~39 s/step, ~11 h/48 h), 13484023 2r_ige (~149,
    ~42 h/72 h), 13484024 4r_oge (~71, ~20 h/72 h), 13484025 4r_ige
    (~286 s/step, **~80 h vs 72 h wall**). Banners verified incl.
    damp_band 0.1R on IGE (breach ruling applied; Phase-3 verdict VOID —
    both arms overflowed 500k ~rev 15–18 and blew up, ppdamp with 32k
    below-ground particles mid-blow-up).
  - **Ruling (Ryan): ride 13484025 to the wall + warm-start rescue** — it
    writes restart VTK every step; ground warm-start validation (052b
    Phase A.1 follow-on; OGE dedup landed, ground still blocked at driver
    :1630) must land within ~3 days to continue the lost tail from the
    last restart set. Do NOT cancel it.
- 2026-08-25 — Phase staged (Ryan request): extend FastMultipole 052/052a
  for multi-rotor + IGE GPU support; success = 10 revs for 1/2/4 rotors,
  OGE and IGE, ≤ 2 h walltime per case, major cost steps on GPU.
  FastMultipole item 052b created and registered in its START_HERE.md.
  Ruling 7 recorded (single shared particle field; per-rotor pfield
  construction at driver :815 must be updated, minimal change expected).
  No jobs submitted; no code changed.
