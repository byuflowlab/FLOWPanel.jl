# RESET PROMPT — 018 particle-Γ distribution vs NT (2026-09-15, v2)

You are a clean-context agent working BRAINSTORM item 018 (DJI-9443 hover
convergence campaign; skim the header + standing rulings of
`../018_dji9443_hover_convergence_campaign.md`, then `ops_reference.md` and
`decision_rules.md`). Never tune to CT_exp; ask Ryan before writing anything to
the lab notebook or `ledger.md` (both stale since 2026-09-08).

## Where the campaign stands (2026-09-15 session — do NOT re-derive)

- NT anti-converges: unguarded reference climb +1.51%/doubling
  (NT36 0.070775 → NT72 0.071844 → NT144 0.07383). All single-mechanism
  discriminators are nulls or worse (see `nt_discriminators_handoff_20260906.md`
  §17 and `expguard_reset_prompt_20260908.md`).
- **NEW (this session): exact-rate `SFS_RLXF` is NOT the carrier.** A clean
  guarded same-arch (gh200) pair, tag `campaign/p018-rlxfscaled-20260915`,
  provenance `rlxfscaled_provenance_20260915.md`:
  - job 13704962, run `p018_csarc_l3p0_3r_g25` (NT36, rlxf 0.005):
    CT̄(21–30) = 0.070519 [0.070478, 0.070577]
  - job 13704963, run `p018_csarc_n2_nt72_l3p0_3r_srlx_g25` (NT72,
    rlxf 0.0025031 = exact rate): CT̄ = 0.071417 [0.071306, 0.071548]
  - **Corrected climb +1.27%/doubling** — only ~16% below the reference climb;
    the srlx-inferred "39% of the climb" is superseded (it was a cross-stack
    comparison). Residual carrier unknown.
- **NEW (this session): gross sheds/rev is NOT constant across NT**:
  30,993 / 34,222 / 38,357 per rev at NT36/72/144 on `_3r_sv_s1p5`
  (+10.4%/doubling; `max(1,ceil)` station quantization). Birth σ IS constant
  (<1%: 6.92/6.96/6.95e-3 m). Audit script: `~/p018_sigma_shed_audit.py` on
  the cluster.
- Untested front runners for the residual climb: FMM near-set adequacy bias
  (Ryan hold since 09-05, `nt_discriminators_handoff_20260906.md` §2) and
  NT-axis spatial entanglement (shed-count excess + removal-side asymmetry,
  `phase_17_nprop_nt_ladder.md`).

## Your task (Ryan, 2026-09-15)

Check the **per-particle Γ (strength) distribution in the wake at different
NT** and determine WHERE it changes: near the rotor disk, in the far field,
everywhere — or something else (e.g. a shape change at fixed location, or a
count-vs-strength tradeoff).

Ladder: **use the new g25 pair** (runs above — same stack, guarded, both
COMPLETE, and their FULL per-step VTP sets, 1080/2160 files, were still live at
2026-09-15; protect or harvest before the VTK sweeper culls to newest-36 —
see `vtk_protect_list.txt` and `scripts/p018_vtk_sweeper.sh`). Cross-check any
headline finding on `_3r_sv_s1p5` NT36/72 (older clean pair, unguarded; only
last ~5 steps live, earlier steps in the archive tarball named in each run
dir's `ARCHIVED.txt`).

Method sketch (adapt as needed):

1. Matched revolutions (rev = step/NT): e.g. revs 8, 16, 24, 30 →
   steps 288/576/864/1079 (NT36) and 576/1152/1728/2159 (NT72).
2. For each snapshot: read `Points`, `gamma`, `sigma`; bin by axial distance
   x below the disk (disk at x ≈ 0, wake convects to +x, extent ~3.5R,
   R = 0.11995 m; 0.5R bins is a good start). Do NOT choose the axis by max
   extent — radial spread exceeds axial.
3. Per bin and per rung, report: particle count, per-particle |Γ| percentiles
   (p5/p25/p50/p75/p95) and mean, Σ|Γ|, and ‖Σ Γ⃗‖. Optionally Γ-weighted σ.
4. Compare NT36 vs NT72 per bin. IMPORTANT: NT72 sheds ~10% more particles
   per rev (quantization excess) — distinguish "more, weaker particles with
   the same Σ|Γ|" (benign discretization) from "different Σ|Γ|" (real
   circulation discrepancy). Report both raw and count-normalized views.

Starting observation to verify and extend (single snapshot, rev ≈ 30,
`_3r_sv_s1p5`, Σ|Γ| in 0.5R bins [0–0.5R … 3–3.5R]):

| rung | 0–.5 | .5–1 | 1–1.5 | 1.5–2 | 2–2.5 | 2.5–3 | 3–3.5 |
|---|---|---|---|---|---|---|---|
| NT36  | 0.461 | 0.442 | 0.391 | 0.314 | 0.612 | 0.666 | 0.914 |
| NT72  | 0.459 | 0.410 | 0.402 | 0.391 | 0.527 | 0.518 | 0.791 |

Near-disk agrees ≲1%; beyond ~2R differs 10–20% — hints at a wake-evolution
(far-field) difference rather than a shed-side one, but it is ONE snapshot
with no per-particle distribution breakdown. Interpretation guide:
- differs immediately after shed → shedding/conversion mechanism
  (phase_17 candidates 1 & 4, untested);
- consistent near disk, diverges downstream → wake-evolution mechanism
  (SFS/relaxation/removal; matches the documented removal-side asymmetry:
  old wake eaten ~2× faster at NT36);
- everywhere / shape change at fixed Σ → discretization-density effect
  (ties to the shed-count quantization excess).

## Data and tooling (ORC cluster, `ssh orc`, user rander39)

- Run dirs: `~/projects/FLOWPanel.jl/data/<run_name>/`; particle VTPs in
  `<run>_wake1_particles/<run>_wake1_particles.<step>.vtp`.
- Reader: `~/vtp_C.py` (`read_vtp(path, {"Points","gamma","sigma"})`,
  `run_files`, `step_of`) — VTPs are raw-appended VTK XML; meshio CANNOT read
  them. Example use: `~/p018_sigma_shed_audit.py`.
- Remote python needs a login shell for numpy: `ssh orc 'bash -lc "python3 …"'`;
  same for slurm. The login banner's ANSI codes glue to the first stdout line —
  start remote commands with a bare `echo` and never anchor greps to line
  starts on the first line.
- Print summaries/tables only; never cat VTP/CSV bytes into context. Archive
  tarballs live on `/nobackup/archive/usr/rander39/FLOWPanel_runs/…` — extract
  selected members to a scratch dir, never into the run dir.
- Scorer (if you need CTs): `python3 scripts/p018_analyze.py m1 --revs 21 30
  <run_name>` from `~/projects/FLOWPanel.jl`.

## Deliverable

Per-rung, per-rev, per-bin table (counts, |Γ| percentiles, Σ|Γ|, vector sum);
a one-paragraph verdict: WHERE the Γ distribution changes with NT (near disk /
far field / everywhere / other), raw vs count-normalized; what mechanism class
that implicates; proposed next discriminator. Record results in a dated status
file in this directory; offer (don't write) a notebook entry and ledger rows —
notebook owes entries for: Cd transient, rlxf derivation, Ladder C forensics,
expguard arc, and the 2026-09-15 session (audit + rlxf-scaled pair result).
