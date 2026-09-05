# RESET BRIEF — smooth-conversion / SFS / σ NT ladders (2026-09-02)

Supersedes `gpu_relaunched_ladders_status_20260831.md` as the entry point for
the NT-convergence investigation. Read that file's §5 guardrails — they all
still apply (sacct authoritative, rsync `--checksum`, `ssh orc` needs
`bash -lc` + live ControlMaster socket, NODE_FAIL restarts are fresh,
exit 124 → `_s2` chain, don't touch `fp-il-*`/`fp052*` jobs).

## 0. Where the investigation stands

The 2026-08-31 ladders all landed clean (exit 0:0). Verdict, recorded in the
notebook (`~/Dropbox/research/notebooks/journals/20260901.md`, `# 20260902`,
"NT Not Converging"): **CT-vs-NT anti-converges** (per-doubling deltas GROW:
+1.3–1.6% at 36→72, +2.2–3.1% at 72→144) and every mechanism tested is a
null — merge-off, fixed rlxf=0.3, and traditional-Pedrizzetti all preserve the
climb. λ converges, NT does not (phase_17_nprop_nt_ladder.md).

| Ladder (λ3.0, 3R) | NT36 | NT72 | NT144 |
|---|---|---|---|
| Reference (corrected-Pedrizzetti, merge on, progressive rlxf) | 0.07078 | 0.07184 | 0.07383 |
| Merge OFF (`_nm`) | 0.07081 | 0.07193 | 0.07413 |
| Fixed rlxf=0.3 | 0.07078 | 0.07171 | 0.07330 |
| Traditional Pedrizzetti (`tp`) | 0.07171 | 0.07273 | 0.07459 |
| `tp` + rlxf=0.3 | — | 0.07242 | 0.07415 |

Key structural facts established this session:

- The `l3p0` Das policy (`DAS_SIGMA_LAMBDA=3.0` on `SIGMA_CHORD_FRACTION=0.313`
  station σ, arc-placed) is **dt-independent**, and the nprop N=1/2/4 row
  scaling holds `d_front = |Das| + (N−1)·travel` constant. First-row placement
  and buffer extent are NOT the lever.
- Legacy branch (N≥1) runs `method_unsteady=NoShed()` +
  `unsteady_filament=false`. VERIFIED in `src/FLOWPanel_wake.jl:2859-2875`:
  the final filament then carries exactly −μ(last row) on that row's aft-edge
  nodes → **net-zero circulation at the sheet/particle interface by
  construction**. The unsteady-circulation carrier is entirely absent in the
  legacy ladders; `SurfaceVorticityConversion` (smooth) switches it on
  (requires `unsteady_filament=true`, constructor-enforced).
- SFS in all ladders so far: canonical `DynamicSFS(Estr_fmm,
  pseudo3level_beforeUJ, pseudo3level_positive_afterUJ; alpha=0.999,
  clippings=(clipping_backscatter,), maxC=1.0, rlxf=0.005)`
  (driver lines ~489-510). Constant across rungs. Note SFS rlxf is
  per-STEP (one of the few per-step-scaled quantities left).
- Smooth conversion enforces particle overlap AT DEPOSITION in both sheet
  directions: bilinear subdivision at h = σ/CONVERSION_OVERLAP
  (`_subdivision_counts`, wake.jl:1190; perimeter `p_per_step =
  ceil(overlap·dist/σ)`, wake.jl:2415). No post-shed enforcement (same as
  legacy).

## 1. Ryan's orders — FOUR new NT ladders (12 GPU jobs), in this sequence

All identical to the reference `l3p0` construction (same case tags, merge on,
3R, corrected-Pedrizzetti progressive rlxf, settle 22) except as noted.
CONVERSION_SIGMA = 0.0381R with R=0.119 m → 0.0045339 m (σ* from BRAINSTORM
019, Ryan-selected; pinned constant across rungs — the driver default is
NT-dependent, never leave it unset).

| # | Arm | Extra `--export` vars | Run-name suffix |
|---|---|---|---|
| 1 | smooth conversion, matched overlap | `CONVERSION=smooth,CONVERSION_SIGMA=0.0045339,CONVERSION_OVERLAP=2.75` | `_3r_sv` |
| 2 | smooth, leaner overlap | `CONVERSION=smooth,CONVERSION_SIGMA=0.0045339,CONVERSION_OVERLAP=2.0` | `_3r_sv_h2p0` |
| 3 | legacy shedding, no-backscatter SFS | `SFS_THREELEVEL=true` | `_3r_sfs3nb` |
| 4 | smooth, 1.5× σ | `CONVERSION=smooth,CONVERSION_SIGMA=0.0068009,CONVERSION_OVERLAP=2.75` | `_3r_sv_s1p5` |

Ladder 3 = `SFS_Cd_threelevel_nobackscatter` (Ryan chose it over the
projection variant). Ladder 4 isolates particle-smoothing width; Das law
untouched. ATTRIBUTION stays default (`upstream`; only N=0 requires
downstream). Launch order 1→2→3→4 ("after these", "then").

## 2. Submission recipe (recovered verbatim from job 13542774 via sacct)

Common: submit from the silo top level; case tags per rung
`p018_csarc_l3p0` (NT36) / `p018_csarc_n2_nt72_l3p0` (NT72) /
`p018_csarc_n4_nt144_l3p0` (NT144); common exports
`SIGMA_CEIL=0.030,TRUNCATION_RADIUS_R=3.0,MAX_PARTICLES=1500000,P018_SETTLE_REVS=22`;
run name = `<case-with-nt>_<suffix>` e.g. `p018_csarc_n2_nt72_l3p0_3r_sv`.
Walltimes: NT36 8h, NT72 14h, NT144 24h (padded over tp ladder actuals
2.5/5.4/13h — smooth cost unmeasured). Template (m13h NT72, ladder 1):

    cd ~/FLOWPanel-018-gpu-h200 && sbatch --job-name=fp-018gpu-n2_nt72_l3p0-3r-sv \
      --partition=m13h --qos=gpu --constraint=intel --gres=gpu:h200:1 --no-requeue \
      --cpus-per-task=64 --mem=192G --time=14:00:00 \
      --output=logs/slurm/slurm-%x-%j.out --error=logs/slurm/slurm-%x-%j.err \
      --export=ALL,SIGMA_CEIL=0.030,TRUNCATION_RADIUS_R=3.0,MAX_PARTICLES=1500000,P018_SETTLE_REVS=22,CONVERSION=smooth,CONVERSION_SIGMA=0.0045339,CONVERSION_OVERLAP=2.75,P018_RUN_NAME=p018_csarc_n2_nt72_l3p0_3r_sv \
      examples/run_dji9443_hover_ct_gpu.slurm.sh h200 p018_csarc_n2_nt72_l3p0

Pool variants:
- **eng**: `--partition=eng --qos=eng --gres=gpu:h200:1` (drop `--constraint`),
  same h200 silo/arch. Single node — our eng rungs pipeline behind each other.
- **mgh**: `--partition=mgh --qos=normal --constraint=arm --gres=gpu:gh200:1
  --cpus-per-task=72`, arch arg `gh200`, submit from
  `~/FLOWPanel-018-gpu-gh200`. mgh maxtime 1-00:00:00.

**Node spread (Ryan's order + ETA probe 2026-09-02 23:10):** eng and mgh start
IMMEDIATELY (eng preempts standby jobs — we hold preemption privilege there;
both mgh GH200s idle); m13h ~4 h wait. Assign per ladder: NT36 → mgh,
NT72 → eng, NT144 → m13h (rotate if capacity shifts; re-probe with
`sbatch --test-only`). **mgh gets NO NT144 rungs** (mgh-1-2 NODE_FAILed 3
NT144 attempts on 08-31; mgh-1-1 rebooted after a 61.6 GB GPU-memory leak,
usable but watch it; leak signature: "insufficient free GPU memory" with
nvidia-smi showing GB used and no processes).

## 3. BLOCKING PRECONDITIONS — do these before any submission

1. **Driver edit is LOCAL-ONLY and uncommitted** (matches campaign precedent:
   mirror to silos, no commit without Ryan). The edit
   (`examples/rotor_hover_pressure_comparison.jl` ~line 575): SIGMA_CHORD_FRACTION
   no longer errors under CONVERSION=smooth — station σ then drives Das
   placement only, `method_trailing` untouched, and it prints "Sigma chord
   mode with CONVERSION=smooth: ...". Mirror this ONE file to BOTH silos
   (`~/FLOWPanel-018-gpu-h200/examples/`, `~/FLOWPanel-018-gpu-gh200/examples/`)
   with `rsync --checksum`, verify `md5sum` both sides. Pre-edit silo diff vs
   local was ONLY this edit + the already-committed NaN-gate fmt fix, so
   whole-file mirror is safe.
2. **Verify the second smoke run** before believing the config works:
   `p018_smoke_sv2_local` (7 steps, NT36, smooth @ 2.75, CONVERSION_DIAGNOSE,
   SAVE_VTK=true, monitors off; log
   `/private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/49793aa5-d9b6-457a-bc97-805657d2f15f/scratchpad/smoke_sv2.log`,
   was at step 1/6 ~189 s/step at reset; data lands in repo `data/p018_smoke_sv2_local/`).
   The FIRST smoke exited 0 with finite forces and the construction accepted the
   combo, BUT 7 "fmm! ... sources or targets are empty" warnings landed at
   steps 4–5 (the LAST steps) — **it is unconfirmed that the smooth conversion
   deposited ANY particles**. Check: wake particle VTPs exist with
   NumberOfPoints > 0 by the final steps, and
   `p018_smoke_sv2_local_conversion_diagnostics.csv` shows nonzero deposits
   with small residuals. If zero particles → STOP, debug the conversion
   trigger (`_convert_to_particles!` requires nwakes[] == capacity;
   wake.jl:2530-2535) before burning GPU time.
3. **GPU-path risk, disclosed to Ryan**: SurfaceVorticityConversion has never
   run against the device-resident CuArray pfield (`VPM_ARRAYTYPE=cuarray`).
   Local smokes are CPU. If a device overload is missing the NT36 rung fails
   minutes after first shed — treat ladder-1 NT36 (mgh) as the canary before
   relying on the other rungs' walltimes.
4. After submitting jobs that write VTK, kick `hpc-storage` per HPC.md.

## 4. Verification & readout (per arm, unchanged recipes)

- On first RUNNING: fresh `*.metadata.toml` (mtime AFTER job start — stale
  same-name dirs from Aug runs sit in `~/projects/FLOWPanel.jl/data/`), and
  the log must echo `CONVERSION=smooth, CONVERSION_SIGMA=..., OVERLAP=...`
  (or `SFS` choice for ladder 3) plus "Sigma chord mode with CONVERSION=smooth".
  Note: jobs cd into the silo but run data lands under
  `~/projects/FLOWPanel.jl/data/` (launcher redirect; arc table
  `data/p018_cs_l3p4_rs1_te_downwash_te.csv` resolves there — now also copied
  into the LOCAL repo data/ for the smokes).
- Readout: `<run>_CT_per_rev.csv`, mean CT_mean over
  `in_convergence_window==true`; mid-run
  `monitors/*_monitor02_force_system1.csv` per-rev blocks.
- Decision: any arm whose NT36→72→144 climb flattens vs the reference column
  fingers that arm's mechanism (interface carrier / overlap / backscatter /
  smoothing width). If ALL still climb → back to
  `phase_17_nprop_nt_ladder.md` candidate list with four more nulls.

## 5. Local smoke env (for reproducing/debugging)

    NREVS=0.2 NT=36 RPM=5400 RHPC_MESH=45_185_ct4 FREESTREAM_RAMP_REVS=0.05 \
    FREESTREAM_HOLD_REVS=0.0 FREESTREAM_WITHDRAW_REVS=0.05 SETTLE_REVS=0.05 \
    NWAKEROWS=1 OVERLAP=2.75 P_PER_STEP=12 MERGE_R_FACTOR=0.0055 \
    PARTICLE_SHEDDING=sigma_overlap SIGMA_CHORD_FRACTION=0.313 SIGMA_FLOOR_R=0 \
    DAS_SIGMA_LAMBDA=3.0 DAS_ARC_PLACED=true DAS_ARC_HELIX_SOURCE=steady \
    DAS_ARC_TABLE=data/p018_cs_l3p4_rs1_te_downwash_te.csv \
    CORE_SPREADING_ACTIVE=true WAKE_CORE_BETA=1e9 CONVERSION=smooth \
    CONVERSION_SIGMA=0.0045339 CONVERSION_OVERLAP=2.75 CONVERSION_DIAGNOSE=true \
    TRUNCATION_RADIUS_R=3.0 SIGMA_CEIL=0.030 RUN_MONITORS=false \
    P018_RUN_NAME=p018_smoke_sv2_local SAVE_VTK=true \
    julia --project=. -t 4 examples/rotor_hover_pressure_comparison.jl

(≤4 threads locally. First smoke `smoke_sv.log` in the same scratchpad dir.
Smoke data dirs `data/p018_smoke_sv*_local/` are disposable — delete when done.)

## 6. Open items owed to Ryan

- Notebook: "NT Not Converging" table is WRITTEN (20260902 header, checkbox
  unticked). The new-campaign notebook entry is NOT yet proposed — ask before
  writing.
- Driver-edit commit decision (this repo + silos) pending Ryan.
- mgh-1-2 triple-NODE_FAIL still unticketed (rc.byu.edu/ticket).
- The 08-31 ladder results are harvested from log tails only; formal
  CT_per_rev-window scoring + λ-gate pass over those arms still undone.
