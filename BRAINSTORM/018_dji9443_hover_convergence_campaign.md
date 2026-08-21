# Publishable DJI 9443 Hover Convergence Campaign

**Opened:** 2026-07-30

**Current phase:** Phases 0–2 CLOSED; 3/4/5 executing; 16 (chord–σ co-scaling Das law) OPENED 2026-08-14 (see status below)

**Item-level approvals:** Technical [ ]; clear-context [ ]; user [ ]

## Current status and next actions (as of 2026-08-06 ~10:30 MDT, PREPPED FOR CONTEXT RESET)

> **RESET BRIEF (b) — read this first; then (a) below; then `ledger.md`
> §2026-08-06 entries. No in-session monitors survive a reset — poll
> `squeue` directly.**
>
> **RESET BRIEF 2026-08-18 (l, evening) — READ THIS FIRST; supersedes (k)
> ((k)'s standing rules bind; the (k) pause-window note below is HISTORY —
> pause+resume both executed cleanly, see pause_manifest_20260815.md).
> Work from disk; monitors die with the session — re-arm sacct polls
> (sentinel-echo; MOTD eats piped lines). Session record: STATUS REPORT
> items 22–27 + ledger §2026-08-17/18 + phase_16 §Log + F1b plan
> `~/.claude/plans/stage-f1b-by-writing-elegant-stonebraker.md`.**
>
> **NINE jobs in flight/pending:**
>
> | jobs | what | action on landing |
> | --- | --- | --- |
> | 13193513/514/515 | F1 trio `p018_cs_f1_l{2p4,3p4,4p8}` (β=0.6 cap), ETA ~Wed eve | **banners NEVER verified (resubmits after the Manifest incident) — verify from logs FIRST** (das_beta:0.6, sigma_chord 0.313, λ per case); then score: doubling-under-cap (M1 raw+quiet, M2), Γ̄ overlay vs uncapped `p018_cs_*` (did the +47% inboard excess collapse?), l3p4 cap-only A/B |
> | 13185008/09/10 | model-def arms nt72@0.3 / nt72@0.16334 / nt144@0.3, walls Wed ~24:00 | nt72@0.3 may COMPLETE (~2160 steps); others TIMEOUT ~rev 28.7 / ~rev 19.1. Score stability vs nt72@0.15 marks (precursor rev 14.5 / death 18.6) + matched windows vs NT36 carrier, labeled MODEL-DEF A/Bs. **nt144 landing ⇒ ASK RYAN about `_s2` chain (standing instruction — never chain unprompted)** |
> | 13206336/337 | `p018_scr_csarc_l{2p4,3p4}` F1b screens (TE table), 24 h, PENDING behind 021 traffic | banner-verify on start (`das_arc:true arc_src:steady arc_table:p018_cs_l3p4_rs1_te_downwash_te.csv` + F1b ACTIVE line, max-|u| 0.102 tip); tail-judge ABSOLUTES (sigma-ratio cols NaN); **PASS ⇒ submit TE ladder `p018_csarc_l{2p4,3p4,4p8}` 48 h/64G SETTLE 22 (pre-authorized per the approved F1b plan §Cases); IGNITE ⇒ stop, bring to Ryan (no pre-registered floor for the arc law)** |
> | 13206338 | `p018_csarc_mid_l3p4` (Ryan's midpoint-table sanity A/B, 30 revs, 48 h) | score m1+m2 vs `p018_csarc_l3p4` when both exist — "behavior doesn't change much" is the pass criterion (his words) |
>
> **F1b state (implemented + launched this session):** law =
> endpoint-on-arc via arc-length integration (`_set_Das_station_arc!`,
> src/FLOWPanel_simulate.jl; kwarg `set_Das_station_drifts`); driver knobs
> `DAS_ARC_PLACED/DAS_ARC_HELIX_SOURCE(steady|kinematic)/DAS_ARC_TABLE`;
> unit testset "endpoint-on-arc Das" green; Route A multi-segment REJECTED
> by Ryan; live/converge mode DEFERRED (Kutta operator cache); clamp OFF
> (Ryan pick — the cs table's tip upwash is REAL physics; clamp
> generalization beyond axial deferred on his go-ahead). **Steady tables**
> (probe monitor `examples/p018_te_downwash_probe_replay.jl` under
> pnl.replay, WAKE-only sources — body bound field EXCLUDED, 0.5–0.9 tip
> at a TE point sample): 4 tables local+cluster
> `data/p018_{cs_l3p4_rs1,ufront_n1_s040_visc_s3}_te_downwash_{te,mid}.csv`.
> **Ryan picks 2026-08-18: TE location + cs_l3p4_rs1 source + clamp OFF.**
> Probe findings: mid-span downwash 0.07 tip; `_s3` (rev 58) has inboard
> FOUNTAIN FLOW that rev-30 cs lacks (both locations = real maturity
> effect); locations agree 3–14% mid-span.
>
> **Ops:** deployed md5s local==cluster: driver `00bda73d…`, launcher
> `f2ec5789…`, simulate.jl `a398bb71…` (deploy was PROVEN cluster-base+
> F1b-hunks-only; simulate.jl carries other-session hunks in simulate!/
> shed_wake! — re-verify before any future whole-file deploy), probe
> script `75c2c499…`. **Cluster Manifest must stay Julia-1.11.7-resolved**
> (2026-08-18 incident: 021 re-resolved under 1.12.6, killed the first F1
> trio at load; snapshot-restored + ILUZero kept; their copy =
> Manifest.toml.p021_julia112_20260817; a 021 resubmit under 1.12 can
> re-break it). Disk unchecked since ~114G — re-check with 9 jobs.
> Queue carries heavy 021 traffic (p1-*, some 3-day walls) — expect
> csarc queueing delay. Standing: knobs from banner, judge from monitors
> CSV, exit 0 ≠ health, notebook needs Ryan, ≤4 local threads (Ryan
> killed one smoke over resource pressure — coordinate local runs).**
>
> **⚠ PAUSE WINDOW (HISTORY — executed cleanly 08-15→08-17): ALL rander39
> HPC jobs are scancel'd at 23:50 MDT 08-15 (suspend is
> permission-denied on orc; kills authorized) and warm-restarted 00:00
> MDT Mon 08-17 from `pause_manifest_20260815.md` (this directory) per
> ops_reference chain recipe. Both are in-session crons (78a94d26 /
> 4a0fcf9d) — if the session died in between, DO THE RESUME MANUALLY
> from the manifest: jobs found dead between those times were killed ON
> PURPOSE, not by ignition; do not post-mortem them as failures.**
>
> **RESET BRIEF 2026-08-15 (k, night) — READ THIS FIRST; supersedes (j)
> ((j)'s standing rules + (h)/(g) history still bind). Work from disk;
> no monitors survive a reset — re-arm sacct polls (ssh orc, `bash -lc`,
> sentinel-echo first; MOTD eats first piped lines). Full session record:
> "## STATUS REPORT for Ryan (2026-08-14 evening)" items 1–21 at the
> bottom of this file + ledger §2026-08-14/15 + phase_03 §2026-08-15 +
> phase_15 (c) + phase_16 §Log.**
>
> **SIX jobs in flight (submission-time banners verified except where
> noted):**
>
> | job | case | note / action on landing |
> | --- | --- | --- |
> | 13178762/63/64 | `p018_cs_l{2p4,3p4,4p8}` Phase-16 ladder, 48 h, land ~08-16 evening | wake-health TAIL first (sigma-ratio cols are NaN under station-σ — judge min_sigma/max_u ABSOLUTES vs carrier calibration in phase_16 §Log), then score per phase_16 pre-registration + §Post-ladder rule (M1 raw AND quiet-limit, M2, P1 doubling, P3 A/B vs `ufront_n1_s040_visc`); **l4p8 was running hot** (\|Γ\|/σ² 6.7e3 @ rev 18.8, rising; if it ignites that's a θ-correlated P2 signature ⇒ F1 candidate, NOT a screen bust) |
> | 13183888 | `p018_upin_nt72_rlxf0p3` (MODEL-DEF change, stock per-step rate), 72 h | stability read vs nt72 marks (precursor rev 14.5 / catastrophe 18.6); if it survives to TIMEOUT (~rev 24–30), matched windows vs NT=36 carrier, labeled model-def A/B |
> | 13183889 | `p018_upin_nt72_rlxf0p16334` (MODEL-DEF, exact continuous-rate 1−√(1−0.3)) | same; tests whether +9% dose over ignited 0.15 changes outcome |
> | 13183998 | `p018_upin_nt144_rlxf0p3` (MODEL-DEF, NT=144 @ stock rate), 72 h ⇒ ~rev 14–15 | **banner NOT yet verified — verify on start** (NT:144 pps:3 rlxf:0.3 das_uniform:3.4); **when it finishes: ASK RYAN whether to chain `_s2` (his standing instruction — likely yes if 13183888 looks good); never chain unprompted** |
>
> **Verdicts already logged this session (do not re-derive):** Phase-16
> screens PASS ⇒ ladder submitted; trunc6 H2 NOT killed (drift
> depth-dependent, level −1.2..−2.7%, fill-confounded) and its 6R `_s2`
> IGNITED rev 34.5 ⇒ **3.5R truncation deletion is LOAD-BEARING for
> stability**; upin_nt72 (rlxf 0.15) IGNITED rev 19.6 ⇒ dt ladder STOPPED
> per pre-registration (Ryan's model-def arms above are his directed
> follow-up, not a ladder resumption); nt144 harvested healthy @rev 14.45
> ⇒ **burst attribution 3-point verdict: fluctuation mode NUMERICAL**
> (within-rev std ×0.12 over NT 36→144; rev-14 burst absent at 144;
> caveats: pre-settled window, CT̄ NON-monotone 0.07513/0.07490/0.07775).
>
> **RYAN'S decisions (stage, don't enact):** M1 observable (attribution
> verdict supports quiet-limit; his call), 13183998 chain (ASK), Phase-16
> fallback selection if ladder fails, plan-16 items (exact-rate/spanwise
> — spanwise blocked on verifying the "80-station contraction" claim vs
> phase_10), notebook writes. **Adopted+enacted (binding): minimal-run
> plan items 12–15** (STATUS REPORT): Phase-16 close-out protocol, single
> selected-rung extension to 45–60 revs if unsettled at 30, do-not-run
> list (no low-rlxf reruns, no more upward rlxf, no old-nt144 chain, no
> new truncation/merge ladders, new tags for model changes), acceptance
> language (settings blind to CT_exp).
>
> **Ops:** local harvests current through nt144 (`data/p018_scr_cs_*`,
> `p018_trunc6_*`, `p018_upin_nt72`, `p018_upin_nt144`). Launcher md5
> 9a5b2f2cb459d75a1b5986be56493f02 local==cluster (arms through
> `p018_upin_nt144_rlxf0p3`). Disk ~166–180G of 200G, 6 jobs ⇒ sweep per
> HPC.md (never newest-10 of a live run). Deploy-skew quarantine still
> in force: NEVER whole-file rsync `src/` from this worktree
> (`_wake/_replay/_liftingbody/_metadata/_solver/_abstractbody/
> _elements_fmm` carry 021-session diffs; reconstruct = cluster base +
> intended hunks; md5 cluster copies BEFORE overwriting).**
>
> **RESET BRIEF 2026-08-14 (j) — SCHEDULED ~19:00 MDT RESUME (Ryan absent;
> work from disk, not conversation memory). Ryan's standing instruction:
> make best-judgment decisions and LOG every one — append each decision +
> result to a dated "## STATUS REPORT for Ryan (2026-08-14 evening)"
> section at the bottom of THIS file as you go; that section is what he
> reads on return. The HPC ControlMaster ssh socket should still be live
> (`ssh orc`, login shell `bash -lc` for slurm; MOTD eats first piped
> lines — sentinel-echo first).**
>
> Duties, in priority order:
>
> 1. **Phase 16 screens** (`p018_scr_cs_l2p4/l3p4`, resubmitted after the
>    13170749/50 deploy-skew UndefVar failure — CURRENT job IDs are in
>    `phase_16_chord_sigma_coscaling.md` §Log and the ledger's 2026-08-14
>    entries). They land ~mid-day; judge the wake-health TAIL over all ~8
>    revs (min_sr trend, max_u, σ<0, |Γ|/σ²), plus merge counts (tip-merge
>    aggressiveness flag) — never exit-0. **PASS ⇒ submit the 30-rev ladder
>    `p018_cs_l2p4/l3p4/l4p8` (48 h/64G, P018_SETTLE_REVS=22) — Ryan
>    PRE-AUTHORIZED. IGNITE ⇒ pre-registered contingency: re-screen at
>    `SIGMA_FLOOR_R=0.030` then 0.035 — also pre-authorized; record the
>    crossing station.** Banner-verify every submission (sigma_chord 0.313,
>    das_lambda, θ table, sigma_overlap, visc:true); scancel on mismatch.
> 2. **Brief (h) jobs**: harvest + score whatever landed of 13157751
>    (trunc6, H2), 13157833/834 (upin nt72/nt144 — per phase_03
>    §PRE-REGISTRATION + §STAGED burst attribution), 13159912/3 (rlxf-up
>    pair, stability screen first then matched ~50–60 m1+m2 vs `_s3`).
>    nt144 TIMEOUT ⇒ chain `_s2` per ops_reference.
> 3. **Decisions that are RYAN'S, not yours** — prepare, stage, and
>    recommend in the STATUS REPORT, but do NOT enact: the M1-observable
>    choice (raw vs quiet-limit; gated on upin burst attribution), any
>    Phase-16 fallback selection (F1/F2/F3), notebook entries (never write
>    the notebook without approval), anything destructive to prior runs.
> 4. Standing rules: knobs from the log banner; judge from
>    `data/<run>/monitors/*.csv`; disk ≤200G per `agent_policies/HPC.md`
>    (VTK sweeps, protect list `vtk_protect_list.txt`, never newest-10 of a
>    live run); deploy skew warning — cluster src DELIBERATELY lags local
>    021-session work in `FLOWPanel.jl/_liftingbody/_metadata/_solver/
>    _abstractbody/_elements_fmm .jl`: NEVER rsync whole dirty-worktree src
>    files; ship cluster-compatible patches only (see ledger 2026-08-14
>    incident).
>
> **RESET BRIEF 2026-08-14 (i) — PHASE 16 OPENED (Ryan-directed): chord–σ
> co-scaling, a new Das law.** After reviewing fig_das ("none of these look
> great"), Ryan directed a new phase: σ_j = s*·c_j and |Das|_j = λ·σ_j, so
> Das/c AND Das/σ are span-uniform simultaneously — the unique one-parameter
> family satisfying all three admissibility bands (chord 0.25–1.5c, clearance
> d/σ 2.6–3.4, curvature), whose pairwise conflicts explain the η/κ/uniform-D
> failures. Full motivation, λ-window derivation (λ ∈ [2.6, 4.8] at
> s* = 0.313), pre-registrations P1–P3, curvature diagnostic θ_j = N·Das_j/r_j,
> stability screen (tip σ 0.0226R < the 0.035R ignition threshold — screen
> GATES the ladder; `SIGMA_FLOOR_R` contingency), and ordered fallbacks
> (F1 curvature cap, F2 Green, F3 016/017 equivalence) →
> **`phase_16_chord_sigma_coscaling.md`**. Cases: screens
> `p018_scr_cs_l{2p4,3p4}` then ladder `p018_cs_l{2p4,3p4,4p8}` (30 revs,
> viscous N=1 carrier; l3p4 is a pure co-scaling A/B vs `ufront_n1_s040_visc`).
> Ryan pre-authorized the HPC submissions. Also this session: 13157532
> `p018_etadas_n4_s040_visc` landed+scored (N=1→4 under η-Das = −2.4% CT,
> ε_Γ 4.8% — N knob dominates the offset attribution; ledger); ufront series
> added to fig_das. Brief (h) below still governs the other five in-flight
> jobs' watch+harvest duties.**
>
> **RESET BRIEF 2026-08-13 (h, evening) — READ THIS FIRST; supersedes (g)
> below (keep (g) for its standing rules and job-action detail). Full
> detail: ledger §2026-08-13 entries + phase_03 §STAGED + phase_15
> §2026-08-13. The previous session's watch-and-harvest SUBAGENT and all
> monitors are DEAD — re-arm `sacct` polls (sentinel-echo before piped
> ssh; MOTD eats first lines) and resume the watch+harvest role. Judge
> runs from `data/<run>/monitors/*.csv`, never stdout; wake-health TAIL
> check before trusting any exit 0.**
>
> **SIX jobs in flight (all banner-VERIFIED at submission):**
>
> | job | case | action on landing |
> | --- | --- | --- |
> | 13157532 | `p018_etadas_n4_s040_visc` (diag OFF) | ~landing 08-13 night: score matched 15–29 vs n1_s040 AND etadas_n1 — completes the 2×2 {N}×{Das-law} square for the carrier's +2–3% offset (an unattributed INTERACTION) |
> | 13157751 | `p018_trunc6_n1_s040_visc` (5.5R) | H2: drift rate vs 3.5R production on matched windows; depth-dependent ⇒ truncation artifact, else H2 killed |
> | 13157833 | `p018_upin_nt72` | pinned dt rung (score per phase_03 §PRE-REGISTRATION) + H1/H3 inventory read + **BURST ATTRIBUTION (Ryan-directed, phase_03 §STAGED — flag its verdict prominently: it gates Ryan's M1-observable decision)** |
> | 13157834 | `p018_upin_nt144` | TIMEOUT ~rev 18 ⇒ chain `_s2` per ops_reference; then Richardson + burst attribution (nt144 rlxf 0.075 is dt-equivalent to stock — see phase_03 caveat) |
> | 13159912/3 | `p018_rlxf0p45/0p6_n1_s040_visc` (warm `_s2`@1619, SETTLE 50) | budget term 9 UPWARD slope: stability screen first, then m1+m2 matched ~50–60 vs `_s3` (= stock 0.3, harvested); compare magnitude vs 006's downward slope (~+0.0011/halving); both healthy ~10 revs past handoff at last check |
>
> **Landed & scored 2026-08-13 (all in ledger + phase files):**
> - **rlxf downward CLOSED**: 0.15 ignited rev 48.2, 0.075 rev ~47 (Γ-ignition
>   at viscous floor, dose-response) ⇒ stock 0.3 LOAD-BEARING at σ=0.04R/N=1;
>   error_budget row 9 updated; upward pair = Ryan's chosen follow-up.
> - **Mergeoff: H4b KILLED** (Σ|Γ| slope +3.10 vs +3.00 %/rev unchanged;
>   particle sink real ~0.8%/rev with no Σ|Γ| footprint; merge ΔCT +0.22%
>   3-rev = term-8 preview).
> - **`_s3` (60 revs): raw drift STEEPENS** (+3.1% monotone 45–57, CT̄ 0.0768)
>   BUT **quiet-limit stationary at 0.072975** (unchanged from _s2-only fit
>   0.072979); rev-45 seam CLEAN; Σ|Γ| now RISING +0.45%/rev in 45–58.
>   ⇒ raw M1 will not settle on this carrier; the stationary observable is
>   the quiet-epoch baseline ≈0.0730 (+1.3% vs exp). **PENDING RYAN
>   DECISION (after burst attribution): M1 observable = raw mean vs
>   quiet-limit.**
> - **Ruling 9 REAFFIRMED (Ryan): M2/Γ̄ stays co-equal — CT-only descope
>   declined.**
>
> **New assets this session:** living convergence figures —
> `~/Dropbox/research/notebooks/img/20260813_018_convergence/` (6 × fig_*.tex
> + CSV dirs + build.sh, gs 300dpi PNGs) generated by
> `scripts/p018_convergence_figs.py` (edit its `AXES` dict when a rung lands,
> rerun, `./build.sh`; notebook refs stable). Notebook `journals/20260803.md`
> § 20260813: CT-convergence status entry + 08-07→08-12 backlog entry +
> figures section (file at 967 lines ⇒ next substantial entry starts a NEW
> journal file). Ledger/error_budget/phase files current through the _s3
> harvest.
>
> **Ops:** disk grows ~15–20G/h with six jobs (sweep at ≥175G, budget 200G,
> `p018_vtk_sweeper.sh`; protect list now clear of the four big closed runs).
> Local uncommitted: launcher rlxf arms (deployed, md5 3ccab7d3),
> `scripts/p018_convergence_figs.py`. **Do NOT push local `src/` to the
> cluster** — local diffs are other-session (021) work; cluster src is the
> carrier-consistent state. Standing rules in (g) below still bind (no local
> jobs without asking; notebook writes need Ryan approval; walltime ≤72 h).**
>
> **RESET BRIEF 2026-08-12 (g, night) — superseded by (h); standing rules
> and per-job detail still apply; supersedes (f)
> and everything below. Full detail: ledger §2026-08-12 (many entries).
> YOUR JOB: watch the SIX in-flight runs and harvest each as it lands.
> No in-session watchers survive a reset — re-arm `sacct` polls for every
> job below (sentinel-echo before piped ssh output; MOTD eats first
> lines). Wake-health TAIL check before trusting any exit 0; judge from
> `data/<run>/monitors/*.csv`, never stdout. Disk watchdog NOT armed —
> sweep VTK manually per retention rules (keep last-10 steps + protect
> list).**
>
> **In-flight jobs (all banner-VERIFIED) and the action each needs:**
>
> | job | case | ETA | action on completion |
> | --- | --- | --- | --- |
> | 13157753 | `p018_mergeoff_n1_s040_visc` (MERGE_PARTICLES=false, warm `_s2`@1619, ~+5 revs, diag ON) | ~night 08-12 | **FIRST re-verify merge-off took**: n_p in its wake_health CSV must exceed `_s3`'s at matched steps (an armed checker died with the session). Then phase_15 (f): per-band Σ\|Γ\| slope (inventory CSV) + CT vs `_s3` over the matched ~5 revs. Σ\|Γ\| decline stops ⇒ merging owns the drift; unchanged ⇒ H4b killed |
> | 13157752 | `p018_ufront_n1_s040_visc_s3` (to ~60 revs, warm `_s2`@1619, SETTLE 50, diag ON) | ~08-14 | phase_15 (d)+(e): H4 trend (drift continue/arrest/reverse through rev 60); seam burst check at the rev-45 stitch (does within-rev std collapse like rev 30?); extend burst/quiet stats + banded inventory series; it is ALSO the mergeoff control |
> | 13157751 | `p018_trunc6_n1_s040_visc` (depth 6R = 5.5R downstream, cold ~30 revs, diag ON) | ~08-14 | phase_15 (c): drift rate vs `n1_s040_visc` (3.5R) on matched windows ⇒ H2 verdict (depth-dependent ⇒ H2; independent ⇒ killed) |
> | 13157532 | `p018_etadas_n4_s040_visc` (N=4 η-Das @σ0.04R visc, **diag OFF** for A/B purity) | ~late 08-13 | score matched 15–29 (`p018_analyze.py` m1+m2 + `p018_burst_stats.py`) vs BOTH `n1_s040_visc` and `etadas_n1` — completes the 2×2 {N=1,4}×{η,uniform-D} interaction square for the +2–3% carrier offset |
> | 13157833 | `p018_upin_nt72` (NT=72/PPS=6/rlxf 0.15, 72 h, diag ON) | ~08-15 | per-rung STABILITY SCREEN first (ignite ⇒ STOP ladder, bring 019 regime map to Ryan); then matched 15–29 M1+M2 per phase_03 §PRE-REGISTRATION; it is the H1/H3 MEASUREMENT run — analyze banded Σ\|Γ\|/ΣΓ⃗/σ trends (stationary ΣΓ⃗ under declining Σ\|Γ\| ⇒ merge cancellation) |
> | 13157834 | `p018_upin_nt144` (NT=144/PPS=3/rlxf 0.075, 72 h, diag ON) | TIMEOUT ~rev 18 ~08-15 | chain `_s2` (72 h allowed) from last VTU step per `ops_reference.md` recipe (VTU last = CSV last − 1); then Richardson {36: existing carrier, 72, 144} matched 15–29 |
> | ~~13157881~~ | `p018_rlxf0p15_n1_s040_visc` | **IGNITED rev 48.2** (+3.2 revs post-handoff), scancel'd 08-13 | downward ladder CLOSED — see error_budget row 9 + ledger §2026-08-13 |
> | ~~13157882~~ | `p018_rlxf0p075_n1_s040_visc` | **IGNITED rev ~47**, exit 0 with garbage (wake self-annihilated 217k→11k) | ↑ dose-response: quartered ignites ~1.2 revs sooner than halved; stock rlxf=0.3 is LOAD-BEARING at σ=0.04R/N=1 |
> | 13159912 | `p018_rlxf0p45_n1_s040_visc` (RELAX_RLXF=0.45, warm `_s2`@1619, +15 revs matched to `_s3`) | ~08-14 | UPWARD rlxf pair (Ryan 2026-08-13, after both downward rungs ignited): m1+m2 matched ~50–60 vs `_s3` (=0.3); upward CT-vs-rlxf slope, judged against the 006 legacy-carrier downward slope (~+0.0011/halving) → carry term 9 as measured local sensitivity |
> | 13159913 | `p018_rlxf0p6_n1_s040_visc` (RELAX_RLXF=0.6, otherwise ↑) | ~08-14 | ↑ same |
>
> **NOTE (Ryan, 2026-08-13): when the upin rungs (13157833/834) land, run
> the BURST ATTRIBUTION check** — the pre-registered numerical-vs-physical
> discriminator (band 0.09–0.17, episode timescale ~5–10 revs, CT per-rev
> std across NT 36/72/144: a numerical mode moves with Δt, a physical wake
> mode does not). Its verdict decides whether bursts belong in the M1 mean
> at all — i.e. whether the carrier's settledness claim rests on the raw
> mean or the stationary quiet-limit ≈0.0730 (`_s3` harvest, ledger
> §2026-08-13). Record the verdict in the ledger + fold into the min_sr
> attribution.
>
> **NEW (Ryan, 2026-08-12 evening): ruling 9 REAFFIRMED after a CT-only
> descope was proposed and declined — M2 (Γ̄(r/R)) stays co-equal. The rlxf
> ladder above was Ryan-directed to convert budget term 9 (≈ −0.005 CT
> one-sided, the largest named term) into a measured slope; decide the
> carrying rule when results land. Queue = 8 study jobs. Review context:
> `~/.claude/plans/review-progress-on-018-jiggly-sky.md`.**
>
> **Standing rules (Ryan, this session):** NO local jobs without asking;
> walltime ≤72 h allowed (ops_reference §"Walltime allowance"); notebook
> 08-07→08-12 ON HOLD (do not write; re-offer at a milestone);
> GeometricTools port deferred until AFTER 018 (see "Deferred" section);
> cluster `test/runtests.jl` carries a documented ssw-test exclusion.
>
> **Landed & scored this session (ledger §2026-08-12):** etadas A/B NULL
> in CT (−0.26%, Γ̄ redistribution 4.05% max) ⇒ offset = interaction, dt
> ladder pins at N=1 uniform-D rung 1 = `n1_s040_visc(+_s2)`; p15
> instrumentation gate PASSED take-8 (bit-identity PASS; WakeInventoryMonitor
> = `monitor05_wake_inventory` CSV, wake_health gains
> p1_sigma_ratio/argmin columns); burst/quiet decomposition
> (`scripts/p018_burst_stats.py`, REPORT-ONLY, s* pooled per carrier
> family only): carrier quiet-series drift −0.6 to −1.0%±1.5% ⇒ quiet
> baseline consistent with stationary; merge-cancellation (H4b) = leading
> Σ\|Γ\|-decline candidate; under-damping-near-σ_stab DISFAVORED; σ stays
> 0.04R. dt-ladder pre-registration APPROVED as drafted (phase_03
> §PRE-REGISTRATION).
>
> **Tools:** `p018_analyze.py m1 RUN --revs LO HI [--stitch SEG --restart-step S]`,
> `m2 A B --revs-a .. --revs-b ..`; `p018_burst_stats.py RUN[+SEG@STEP] ...
> [--lo --hi]`. Validated anchors: b0(+_s2@719) 15–29 = 0.071702;
> n1_s040(+_s2@1079) r≈0.67, quiet limit ≈0.0730.
>
> - **etadas A/B LANDED + SCORED: Das-law ΔCT ~NULL** (−0.26%, matched
>   15–29; M2 shows Γ̄ redistribution 4.05% max at unchanged CT) ⇒ carrier
>   offset = INTERACTION; **dt ladder pins at N=1 uniform-D, rung 1 =
>   `n1_s040_visc(+_s2)`**. Ryan's named N-knob follow-up
>   `p018_etadas_n4_s040_visc` = job **13157532** RUNNING (banner
>   verified; diagnostics OFF for A/B purity; ETA late 2026-08-13; score
>   matched 15–29 vs n1_s040 AND etadas — completes the 2×2
>   {N}×{Das-law} square).
> - **Phase-15 EXECUTING** (plan discussion HELD; decisions in ledger
>   "PLAN DISCUSSION"): WakeInventoryMonitor + attribution IMPLEMENTED,
>   deployed, unit-green; **gate job take-8 = 13157614 in flight**
>   (takes 1–7 died on env/skew — julia 1.11.7 pin, six stale test
>   files, missing fixtures/examples, and the PRE-EXISTING
>   GeometricTools break: see ops_reference "Pre-submission gate" crib).
>   **Four conditionals authorized, submit AFTER gate passes**: trunc6 +
>   mergeoff (launcher arms exist), `_s3` ~60 revs (chain recipe), seam
>   check (rides `_s3`; phase_15 (c)–(f) protocols).
> - **Pinned dt-ladder PRE-REGISTRATION DRAFTED** (phase_03 §"PRE-
>   REGISTRATION (2026-08-12)") — NT 36/72/144 @ PPS 12/6/3, arms
>   `p018_upin_nt72/144` NOT yet created; **awaits Ryan sign-off**.
> - **Drift reframing (burst_stats, ledger):** carrier quiet-series drift
>   −0.6 to −1.0% ±1.5% (consistent with stationary quiet baseline);
>   burst rectification = viscous-N=1 family property; merge-cancellation
>   (H4b) elevated to leading Σ|Γ|-decline candidate (ΣΓ⃗ columns + the
>   merge-off A/B discriminate).
> - **Standing session rules:** NO local jobs without asking Ryan (his
>   2026-08-12 direction); GeometricTools port deferred until AFTER 018
>   (see "Deferred" section); notebook 08-07→08-12 still UNWRITTEN.
>
> **RESET BRIEF 2026-08-12 (e) — superseded by (f); kept for history.
> Full detail: ledger §2026-08-12 entries (three). No in-session
> monitors survive a reset — poll `squeue` directly; re-arm a watchdog
> (sweeper at 100G budget 200G, report only >1G freed; cluster protect
> list authoritative md5 f1de3b12; disk idles ~122G at retention floors,
> so a bare >100G trigger fires every poll — gate notification on freed
> MB).**
>
> **ORDERED NEXT ACTIONS FOR A FRESH AGENT:**
>
> 1. **Check 13135245 `p018_etadas_n1_s040_visc`** (Das-LAW A/B:
>    η-kinematic 1.0 instead of uniform D=3.4σ, all else = n1_s040;
>    banner VERIFIED; was at step 936/1079 at 08:56 MDT 08-12, ~135
>    s/step ⇒ ETA late afternoon 08-12 — likely LANDED by the time you
>    read this). Check the wake-health tail before trusting exit 0, then
>    score ΔCT/M2 vs `p018_ufront_n1_s040_visc` on the MATCHED 15–29
>    window (both 30-rev runs share the settling transient; do NOT use
>    the 45-rev extension's late window for the A/B). ~Null ⇒
>    interaction; Ryan's named follow-up = N=4 η-Das @0.04R viscous.
>    The verdict picks the Das law — and hence which existing run is
>    rung 1 of the dt ladder (below).
> 2. **DISCUSS PLANS WITH RYAN BEFORE IMPLEMENTING ANYTHING** (his
>    explicit direction 2026-08-12). Two staged work items await that
>    discussion: (i) the Phase-15 drift-source study —
>    `phase_15_drift_source.md` (H1–H4, monitor design Ryan-approved at
>    the design level: 5 z/R bands + radial split, σ quantiles, no merge
>    stats in v1, instrumentation GATES the dt ladder); (ii) the pinned
>    dt-ladder pre-registration (windows, gates, per-rung stability
>    screens, cost) — draft it AFTER etadas lands, show Ryan before any
>    submission.
> 3. **Notebook entries for 08-07→08-12 are UNWRITTEN** (ledger-only).
>    Ryan talked through items 1–2 of a proposed 3-item entry this
>    session (Phase-3 close; σ-boundary + arrest + offset revision) but
>    has NOT approved any notebook text — re-offer per notebook policy,
>    including the min_sr/correlation/inventory findings.
>
> **STAGED DECISION (Ryan, 2026-08-12): pinned dt ladder at the §4b
> uniform-d_front N=1 construction — DO NOT design/submit until the etadas
> A/B (13135245) lands and prices the Das-law knob (it may change the
> carrier's Das law; any dt-independent law pins automatically at N=1).**
> **DECIDED (Ryan, 2026-08-12): schedule = NT 36/72/144 with PPS 12/6/3**
> (NT·PPS = 432 pinned ⇒ σ = 0.04R AND h pinned at OVERLAP 2.75 — the
> exact certified carrier operating point; cost ACCEPTED). Rationale
> trail: PPS must scale as 1/NT to pin σ and h; Ryan's PPS floor of 3 at
> NT=144 keeps overlap in the WAKE (row placement itself is
> PPS-independent at N=1 uniform-D; the constraint is streamwise sheet
> sampling — Ladder B priced overlap 1.0 at 7.86% max Γ̄ error).
> **The NT=36 rung already exists**: n1_s040_visc(+_s2, 45 revs) if
> uniform-D survives the etadas A/B, or the etadas run itself if η-Das
> wins — only NT=72/PPS=6 and NT=144/PPS=3 are new submissions.
> Also agreed: stability screen per rung (σ_stab and the viscous floor
> move with Δt; screens showed refinement tightens the ignition
> boundary). Prerequisite: a settledness standard on this carrier (1.1%
> monotone drift at 45 revs; ~0.3–0.5% deltas are unreadable without it —
> consider 019's drift-arrest criterion). Pre-registration (windows,
> gates, screen protocol, cost) to be drafted for Ryan AFTER etadas
> lands; no submission before he sees it.
>
> **STAGED (Ryan 2026-08-12): Phase-15 drift-source diagnosis study —
> full design in `phase_15_drift_source.md`** (H1 far-wake Σ|Γ|
> equilibration / H2 truncation artifact / H3 viscous σ-distribution /
> H4 undersampled ≳30-rev mode, each with kill conditions; H4b
> merge-rate hypothesis deferred, revisit if H1–H4 die). Ryan design
> approvals recorded there: 5 z/R bands + radial split; σ QUANTILES
> included; NO merge stats in v1; **instrumentation gates the dt
> ladder** (one bundled monitor src change before the NT=72/144
> submissions, so ladder rungs carry the diagnostics free).
> **Implementation awaits a plan discussion with Ryan — do not start.**
>
> **NOTE (post-dt-ladder analysis task): numerical-vs-physical
> discriminator for the min_sr / CT aperiodic fluctuation** — after the
> pinned dt ladder lands, compare the fluctuation's band (0.09–0.17),
> episode timescale (~5–10 revs), and CT per-rev std across NT 36/72/144:
> a numerical mode should move with Δt, a physical wake mode should not.
> Comes free with the ladder; record the verdict in the ledger and fold
> into the min_sr attribution (per-rev min_sr–CT correlation is NULL at
> 45 revs — detrended r≈+0.08, no consistent lag structure across arms).
>
> **LANDED 2026-08-11/12 (both scored, ledger §2026-08-12):**
>
> - **Phase 3 CLOSED (13134726 `nt144_s3`):** matched 10–19 CT̄(144) =
>   0.072023 [0.071952, 0.072085]; Richardson {36: 0.071866, 72: 0.072243,
>   144: 0.072023} **NON-MONOTONE — no observed order**; dt residual
>   bounded ≤0.5% p-t-p; order needs a clearance-pinned ladder (phase_13
>   §5). Close-out table in phase_03. Seam-gate gotcha: cut
>   monitor-reconstructed segments at S+1, not S (csv=monitor+1).
> - **13135200 `n1_s040_visc_s2` (45 revs): arrest HOLDS** (global min_sr
>   0.0976 @ rev ~31, monotone recovery to 0.153 @ 45); M1 late window
>   30–44 CT̄ 0.073632, drift 1.107% monotone — **still CHECK**. The 15–29
>   window (0.0753) was transient-inflated (hump revs 25–29 → drop → slow
>   rise); **ufront offset shrinks to +2.7% vs η-Das / +2.3% vs experiment**
>   — much of the "+4–8% offset" was settling transient.
>
> **Standing verdicts (all 2026-08-11 ledger):** σ-blow-up A/B: viscous
> survivors n1_s040 (FIRST arrest, min_sr 0.099→0.122), n1_s050, n2_s050;
> n2_s035 ignited rev 23.3, n2_s040 rev 26 ⇒ **viscous N=2 LESS stable than
> N=1** (reverses inviscid ordering); survival boundary N1 ≥0.04R, N2
> ≥0.05R. **No surviving arm converged**: CT̄ 0.0753–0.0772 all drifting
> monotone (M1 CHECK), cross-M2s ~5% FAIL — the ufront carrier sits +4–8%
> above the η-Das family AND experiment; **this offset is the dominant open
> discrepancy**. h-ladder: 0.5→0.25 PASS (ΔCT −0.02%, ε 0.22%, 3-lobe
> ABSENT); **0.25→0.125 FAILS both gates** (+2.1% CT, ε 5.1%, smooth inboard
> tilt NOT 3-lobe; seam 0.006% PASS, no particle-cap saturation) ⇒ h-claim
> rests on the first pair alone (σ/h=8 deep-overlap implicated); **FLAGGED,
> awaiting Ryan.** σ\*: 0.035R certification FAILED (ignition rev 27.6);
> provisional σ\* = 0.04R N=1 pending Ryan's survival-vs-arrest standard.
>
> **Decisions Ryan owes / next-work menu (do NOT start without direction
> unless listed as pre-authorized):** (1) h-ladder fail-branch handling;
> (2) extend `_s2` the surviving arms for settled windows — n1_s040 first
> (arrested σ\* candidate); (3) investigate the ufront +4–8% CT offset
> (confounds: N 4→1, uniform-D law, viscosity, σ — an A/B isolating one
> knob on the B0 carrier would decompose it); (4) trim-before-merge
> reordering still unapproved (OOM symptom gone via FLOWVPM sparse fix;
> ordering question open); (5) notebook entries for the 08-07→08-11
> results are UNWRITTEN (all results ledger-only; needs Ryan approval per
> notebook policy). Item 019 runs in a CONCURRENT session (P4 grid,
> max_dtZ column) — coordinate only via the cluster protect list,
> append-only; job cap is 20 (their auto-submitter may still check 10).
>
> **UPDATE 2026-08-06 ~17:00 MDT (superseded by the block above; kept for
> history):** σ\*=0.035R DECIDED
> (n1_visc ignited); `scr_ufdt_nt144` ignited step ~1104 and was scancel'd
> (whole ufdt family dead at ≈7–8 revs regardless of dt — ledger §13:20);
> ufdt corpses + 019 entries now on the VTK protect list (cluster copy is
> authoritative — md5 f1de3b12). **NEW: Ryan-directed viscous σ-ladder A/B
> RUNNING, banners verified (ledger §12:15/§12:35): 13061047
> `p018_ufront_n2_s035_visc` (N2 D5.6), 13061048 `n1_s040` (N1 D3.4),
> 13061049 `n2_s040` (N2 D4.9), 13061050 `n1_s050` (N1 D3.4), 13061051
> `n2_s050` (N2 D3.9); all viscous, SETTLE 22, ~ETA evening 2026-08-07.**
> Job cap now 20. Projections (made ~16:45 06 Aug): h0p25 photo-finish
> ~02:00 MDT 07 Aug; h0p125 TIMEOUT ~02:00 at ~step 850 ⇒ chain `_s2`;
> nt144_s2 TIMEOUT ~12:30 MDT 07 Aug at ~step 3130 ⇒ chain `_s3`;
> s035_visc completes ~noon 07 Aug (σ\* certification run); n2_visc +
> the five A/B arms complete ~18:00–22:00 MDT 07 Aug. 019 is ACTIVE in a
> CONCURRENT session (its own screens + auto-submitter; coordinate protect
> list via the cluster copy, append-only).
>
> **In-flight jobs and the action each needs on completion:**
>
> | job | case | action on completion |
> | --- | --- | --- |
> | 13058534 | `p018_ufront_n1_visc` (σ 0.03R viscous, N=1 D=3.4) | **DEAD 2026-08-06 — IGNITED** (max_u →67k m/s, \|Γ\|/σ² 1.4e9, σ pinned at viscous floor; merge OOM = symptom). **Tripwire fired ⇒ σ\* = 0.035R decided** (ledger §11:20). All-viscous A/B at 0.03R lost its N=1 arm |
> | 13060144 | `p018_ufront_n2_visc` (N=2 D=6.5 viscous) | banner VERIFIED (nwakerows:2, das_uniform:6.5, visc:true). **On watch: shares σ=0.03R that ignited its N=1 twin at step ~490** (healthy at step 166, min_sr 0.48). If it survives, it's an N=2-stabilizes datapoint; scoring per phase_13 §7 with the lost-arm caveat |
> | 13058988 | `p018_ufront_s035_visc` (σ 0.0349R, D=3.0, hedge) | **PROMOTED to production candidate by the σ\* decision.** Healthy (step 431, min_sr 0.316). Score vs B0 + as the likely σ₀ probe for item 019's demonstration |
> | 13057254 | `p018_ufront_n2` (inviscid N=2) | **DEAD 2026-08-06 — IGNITED at step ~667/774** (σ<0 mode, min_sr −2.2, max_u 12.8k), OOM 0:125. Time-to-ignition recorded: N=2 delays (667 vs 501) but does not prevent. Corpse VTK swept, restart set kept |
> | 13054739 | `p018_nt144_s2` | projected TIMEOUT ~80 steps short of 3600 — chain `_s3` from last VTU step (recipe `ops_reference.md`), then Phase 3 Richardson {0.071866, 0.072243, CT(144)} matched 10–19, WITH the phase_13 §5 clearance-confound caveat |
> | 13051772/3 | `p018_h0p25` / `p018_h0p125` | Δx closing leg: score vs B0 (M1+M2, windows 15–29); pre-registered: the 3-lobe Γ̄ mode should be ABSENT at fixed σ (it appeared ∝ σ-step in Ladder B). FAIL ⇒ flag Ryan before touching σ\* logic |
> | 13058048 | `scr_ufdt_nt144` | screen; stability-probe readout only (calibration gate FAILED — no aero deltas from screens) |
>
> **Also pending:** ~~stability-probe readouts of landed screens~~ DONE
> (ledger §11:20): `nt72` healthy; `c145` bounded; `s60` marginal 0.075;
> **`c249` and `s80` ended with σ<0 despite exit 0** — mesh refinement
> tightens the inviscid boundary; exit code ≠ stability verdict.
> ~~sweep screen-batch VTK~~ DONE (18/19 pre-swept by watchdog; freed
> `scr_ufdt_nt144` 1.2G live + `ufront_n2` corpse 1.2G; data/ 88G).
> **Watchdog did NOT survive the reset — disk sweeping is manual.**
> **trim-before-merge reordering** is the indicated cheap fix
> (`FLOWPanel_wake.jl:1636–1646`) — NOT yet implemented, needs Ryan
> approval as a src change (note: a FLOWVPM sparse-cell-list merge fix
> landed 2026-08-06 — see notebook # 20260806 — which removes the OOM
> symptom; the ordering question may still matter for deleting runaways).
>
> **NEW ITEM STAGED: `BRAINSTORM/019_sigma_selection_procedure.md`** —
> formalize + demonstrate the general σ-selection procedure (initialize
> max(2σ_eq, σ_stab=√(Γ_vΔt/2π)) → tripwire-probe ΔtZ_max → iterate),
> regime-map parameter sweep showcasing blow-up, publishable, reuse-first.
> Not started; several in-flight runs double as its grid cells.
>
> **Documents current through this morning:** `sigma_blowup_mechanism.md`
> (COMPLETE, confidence table final incl. σ_stab initializer);
> notebook "Core Size Impact on Viscous VPM Stability" (synthesized entry,
> regime flow chart + generator script in img/20260805_sigma_blowup/,
> procedure-form recommendation); phase_14 (screen harvest + failed
> calibration gate); ledger through §~09:50.

> **HEADLINE UPDATE 2026-08-06 (a) — σ-blow-up ROOT-CAUSED and documented;
> screen batch harvested; viscous production candidates in flight.**
>
> - **Root cause complete** (`sigma_blowup_mechanism.md`, confidence table
>   final; notebook entry "σ-Blow-Up Root Cause" written): the discrete map
>   σ←σ(1−ΔtZ) [+ viscous √(σ²+2νΔt)] is provably stable iff ΔtZ<2 with
>   fixed point ≈σ_eq; ignition is a **Γ blow-up at the viscous σ floor**
>   (corpse-measured, floor 9.41e−5 m = predicted √(2νΔt) to 2%); the
>   thin-core feedback hypothesis is REFUTED (ambient strain dominates);
>   merge-OOM = bounding-box symptom (723R/step runaway; merge runs BEFORE
>   trims — reordering is the indicated cheap fix); collapsing set is
>   unreachable by the current merge trigger (0.12%).
> - **Screen batch (phase_14): 8 complete, 6 IGNITED** — five with NEGATIVE
>   σ observed live (tripwire), the code-audit prediction realized. Inviscid
>   screen-condition boundary = σ/R ≈ 0.030. **Calibration gate FAIL** ⇒
>   screens are stability/plumbing/cost probes only; 8-rev aero deltas
>   untrusted (N-axis wrong-sign null + false positive). dt screen ladder
>   lost to ignition.
> - **`p018_L1_visc` harvested: viscosity delta at L1 σ NULL in Γ̄**
>   (ε 0.156%, ΔCT −0.38%); σ-axis 3-lobe is viscosity-independent.
>   Γ6(L1)=0.16%.
> - **In flight (production candidates + remainders):** `p018_ufront_n1_visc`
>   13058534 (σ 0.03R viscous, min_sr 0.381 — floor working);
>   `p018_ufront_s035_visc` 13058988 (0.035R, D=3.0, Ryan's pre-authorized
>   hedge — trigger was met); inviscid ufront pair (n1 contracting, min_sr
>   0.145 — watch); nt144_s2 (tight walltime); h0p25/h0p125; screens
>   nt72/mesh_c249/mesh_s60/mesh_s80/ufdt_nt144.
> - **Recommendation on record**: production = viscous uniform-d_front at
>   σ=0.03R (hedge 0.035R); trim-before-merge reordering + tripwire as
>   safeguards; ΔtZ guard held as insurance. σ=0.02R pursuit closed unless
>   physics demands it.

> **HEADLINE UPDATE 2026-08-05 (f, ~16:20 MDT session) — review + redirection.**
>
> - Fresh-context review of the campaign found **no methodological errors**;
>   two record corrections logged (`ledger.md` §~16:20): the watchdog's
>   "nt144_s2 has no output" reports were wrong — **`p018_nt144_s2` is warm
>   (opens step 2026) and healthy at ~2125/3600**; and **`p018_L1_visc`
>   warm-start is CONFIRMED** (`resuming from step 719`) — the FIRST-ACTION
>   check below is discharged.
> - **Ryan: fixed-κ clearance discriminator ON HOLD.** Next new submissions =
>   the **uniform-d_front pair** (phase_13 §4b): σ=0.03R, `p018_ufront_n1`
>   (N=1, D=3.4) + `p018_ufront_n2` (N=2, minimal in-band D≈6.5). This also
>   runs §6 action 2 (rigid-row-1 vs free-row A/B). Spec + pre-registration in
>   phase_13; **src md5-sync is still mandatory before submission.**
> - **All four landings harvested and scored this session (ledger §~17:00 +
>   §~18:00):** (a) **Ladder B COMPLETE** — σ at fixed h: CT deltas shrink
>   (+0.78% → +0.40%) but M2 fails every pair; the 3-lobe Γ̄ mode appears at
>   FIXED h with amplitude ∝ σ-step (dip −4.1%/−7.9% at r/R 0.74/0.76) ⇒
>   **the mode is a σ phenomenon** (h-ladder null = closing leg, in flight).
>   (b) **Green κ pair verdict NEGATIVE** — green κ 0.41→0.82 = +1.56% vs
>   velocity +1.80%: κ response is formulation-independent, Green does NOT
>   become production; green−velocity Δ small (−0.32%/−0.56%, M2 PASS both) =
>   early Phase-7 Green-Δ points. Supports the clearance reinterpretation.
> - **Phase 14 SCREEN BATCH SUBMITTED (Ryan: slot cap waived): 19 short
>   8-rev jobs (13058031–049)** — calibration ladder (5, gates whether short
>   windows predict settled verdicts), ufront D×N map (5), first CLEAN dt
>   ladder under the ufront law (3), σ ladder incl. the σ_eq 0.02R rung (2),
>   mesh ladder chord 145/249 + span 60/80 (4). Runner
>   `examples/run_p018_screen_hpc.slurm.sh`, scorer
>   `scripts/p018_screen_score.py`; pre-registrations + evidence-class rule
>   (screens never satisfy M1/M2) in **phase_14_screen.md**. First harvest
>   action = the calibration gate, before reading any campaign ladder.
> - **Uniform-d_front pair SUBMITTED** (phase_13 §7): `p018_ufront_n1` =
>   13057253, `p018_ufront_n2` = 13057254 (σ=0.03R, 48 h/64G/SETTLE=12; gates:
>   unit tests, NT=6 smoke, exact-uniform offline d/σ profile, md5-verified
>   deploy of the 4 changed files incl. WakeHealthMonitor — these are the
>   first campaign runs carrying the tripwire). Banner verification in
>   progress. Queue = 6 study jobs.

> **HEADLINE UPDATE 2026-08-05 (e) — the offline clearance diagnostic was run
> (it was the stated "Next" below). It REINTERPRETS the κ-ladder failure.
> Read this first. Detail: `phase_13_handoff_clearance.md`.**
>
> - **The panel→particle handoff distance is measured, and B0 violates the C1
>   clearance bound inboard by ~4×.** `d_front,j = |Das_j| + (N−1)|ω×r_j|Δt`;
>   on the production carrier the profile is **min `d/σ` = 0.713** (r/R 0.100),
>   tip 6.787, median 3.74, against C1's admissible 2.6–3.4. One number
>   explaining the campaign's repeated **inboard-localized** Γ̄ errors (N=1,
>   small `Das`, σ ladder) — all are small-`d` cases.
> - **THE κ LADDER WAS ALSO A CLEARANCE LADDER, entirely inside the
>   inadmissible region** — so its "failure to converge" is not a statement
>   about chord-proportional `Das`:
>
>   | κ | min `d/σ` | harvested CT̄ (10–19) |
>   |---|---|---|
>   | 0.25 | 1.080 | 0.071571 |
>   | 0.41 | 1.388 | 0.072154 |
>   | 0.82 | 2.178 | 0.073455 |
>
>   CT rises monotonically with clearance at a near-constant **≈2.4% per unit
>   `d/σ`** (2.63 and 2.28 across the two intervals). Growing deltas are what
>   climbing out of a clearance deficit looks like, not a divergent parameter.
>   **Caveat: κ and `d/σ` are collinear along this ladder**, so it cannot by
>   itself separate them; the discriminator is a run that moves `d/σ` at FIXED
>   κ (vary `N`, or σ). **Testable prediction: reaching `d/σ` = 3.4 from 2.178
>   costs a further ≈ +2.9% CT.**
> - **Consequence for the in-flight Green fallback (13051802/3):** it was
>   authorized on the reading that κ diverges. That reading is now in question.
>   The Green pair is still worth having, but **do not treat it as the response
>   to a κ divergence** until a fixed-κ clearance run has separated the two.
> - **`d/σ` is now auditable from every run's own banner** (`RHPC_HANDOFF_CSV`,
>   `RHPC_HANDOFF_PROFILE_ONLY=true` for a seconds-long query; placed before
>   `Backslash`, which LU-factors a dense ~10 GB matrix).
> - **Wake tripwire shipped** (`WakeHealthMonitor`, default ON, bit-identity
>   verified, suite green): per-step `n_particles`, `max_u`, `min_sigma_ratio`,
>   `max_gamma_over_sigma2`, `wall_s`. `min_sigma_ratio` starts at exactly 1.0
>   and the rVPM contraction is visible within 4 steps of the first particle —
>   the signal that ran unobserved for ~1000 steps in `p018_L2`. See
>   `ops_reference.md` §"Wake tripwire".
> - **Recommended production config** (arithmetic, not yet simulated):
>   **σ = 0.02R, N = 1, `Das` = 3.4σ span-uniform absolute** — uniform
>   clearance at the C1 bound, `Das`/c in the 014/017 band at both ends, sheet
>   under one local chord everywhere, and a **dt-independent handoff by
>   construction** (`Das_j = D·σ − (N−1)·travel_j` absorbs Δt). Note this is a
>   THIRD `Das` law, neither η nor κ. Gated on the σ-collapse mitigation.
> - **Open, and load-bearing:** the *upper* (sheet-fidelity) bound on `d/c` is
>   unquantified for freely-convecting rows. It decides N=1 vs N=2 and it is
>   what makes B0's σ "disfavoured" rather than "excluded" (B0 misses by only
>   10% at N=2, against the softest number in the analysis).
>
> **Next actions, in order:** (1) fixed-κ clearance discriminator to separate
> κ from `d/σ` (needs a slot; blocks the Green interpretation) — **ON HOLD per
> Ryan 2026-08-05 (f), superseded as next submission by the uniform-d_front
> pair**; (2) quantify
> the `d/c` upper bound (decides N=1 vs N=2, and whether B0's σ is genuinely
> excluded or only disfavoured); (3) S0d σ-collapse localization (gates
> σ=0.02R; Ryan: **locate before mitigating** — merging / adaptive refinement
> are candidates alongside a floor); (4) S0b offline Δx convergence — with the
> time-marched ladder run **as well**, per Ryan, since frozen-state resampling
> cannot see different `h` growing different wake structures.
>
> **⚠ BEFORE SUBMITTING ANYTHING — this session's code is UNCOMMITTED and NOT
> on the cluster:**
> - `src/FLOWPanel_simulate_monitors.jl` — `WakeHealthMonitor`, appended at EOF
> - `src/FLOWPanel_simulate.jl` — per-step wall clock + `flush(stdout)` (~1008)
> - `examples/rotor_hover_pressure_comparison.jl` — `d/σ` profile block (before
>   `Backslash`), and the `wake_health` wiring after the `monitors = …` tuple
>
> `ops_reference.md` mandates md5 verification of `src/` against the cluster
> before submission, and a stale-`src` mismatch has already produced one
> degenerate run. **Sync first, or the job runs without the tripwire.**
>
> **Local assets (no HPC needed for S0d / S0b):**
> `~/p018_L1_ov3_paraview/p018_L1_ov3_wake1_particles/…719.vtp` is present on
> the Mac. Three of S0d's four columns (`σ/σ_shed`, `σ/h_local`, `|Γ|/σ²`) need
> only `X`/`Γ`/`σ` from that file — only `dt·Z` needs a gradient pass.
> Regenerate the `d/σ` and chord tables any time with
> `RHPC_HANDOFF_PROFILE_ONLY=true` (seconds); in chord mode with κ=1.0 the
> `das` column equals the local chord, which is how the chord table was built.

> **HEADLINE UPDATE 2026-08-05 (d) — five jobs departed and were harvested;
> read this before the older text below.**
>
> - **κ ladder FAILED (the campaign's main fork).** Matched 10–19 windows:
>   κ=0.25 → 0.071571, 0.41 → 0.072154 (**+0.81%**), 0.82 → 0.073455
>   (**+1.80%**, ε_Γ 5.22%). Successive deltas *grow* ⇒ chord-proportional Das
>   does not converge. **The pre-authorized Green fallback is already in
>   flight** (13051802/3) — no submission was needed. Detail: phase_02
>   §2026-08-05, `ledger.md` §2026-08-05 (b).
> - **Viscous erratum is NULL at production σ**: `p018_b0_visc` ΔCT̄ = −0.078%,
>   ε_Γ = 0.19% ⇒ M1+M2 PASS. Budget row Γ6 = 0.19%.
> - **`p018_L2_visc` FAILED (diverged, not OOM)** — CoreSpreading did not
>   rescue the σ-L2 rung. Per Ryan: record and **defer the σ decision** to
>   Phase 12A + Ladder B; do NOT trigger the 016 contingency yet.
> - **`p018_nt144` chained**: `p018_nt144_s2` = job **13054739**, restart step
>   2025, 25 revs total, banner verified. Early partial window (revs 10–13,
>   from the monitor fallback) reads 0.072088 — preview only.
> - **`scripts/p018_analyze.py` now falls back to the force monitor** when a
>   run has no `CT_vs_rev.csv` (walltime-killed segments); validated
>   bit-identical on `p018_dasc0p41` (ΔCT = 0.0).
> - **Next**: the free offline clearance diagnostic (phase_02 §2026-08-05,
>   "Code finding") before any new Das ladder; then the Green κ verdict.

A fresh agent should read: this file → `ops_reference.md` → `decision_rules.md`
→ the phase file being executed. No repo exploration should be needed.
`ledger.md` holds every harvested number with provenance. Analysis is
scripted: **`scripts/p018_analyze.py`** (`m1` / `m2` / `windows`, with
`--stitch` for restart segments) implements `decision_rules.md` directly and
is validated by reproducing `p018_b0` = 0.071701 over revs 15–29. Use it
rather than hand-computing windows.

**Recent decisions a fresh agent must know (2026-08-04/05):**

- **Slot cap raised to 10** (Ryan 2026-08-04, supersedes the 6-cap wherever
  stated; unrelated `fm0*` and 016 `fp-016-*` jobs do not count).
- **Root cause of the Phase-2 no-plateau Das ladder found** (phase_02
  §2026-08-04 (b)): Das = η·Δt·Ωr is ∝ r, but the admissible band (014) is
  0.25–1.5 **local chords** — inboard on-plateau needs η ≥ 2.8, tip ≤ 1.5c
  needs η ≤ 2.3, windows disjoint ⇒ no single η can work. **Fix implemented,
  tested, deployed: `DAS_CHORD_FRACTION` = κ ⇒ |Das|_j = κ·c_local per
  station** (dt-independent; freeze factor and floor irrelevant). The κ
  ladder is IN FLIGHT (below) with pre-registered predictions in phase_02.
- **Standing Ryan authorization (2026-08-04):** if the κ ladder fails, the
  agent may execute the pre-registered fallback WITHOUT further approval —
  a **GreenReconstruction Das pair** (`RHPC_FORMULATION=green`, κ 0.41 +
  0.82, B0 carrier); decision rule + rationale (015 Battery II) in phase_02.
- **σ-axis root cause stands** (ledger §2026-08-03 session c): inviscid rVPM
  has no σ→0 limit (core contraction σ̇=−Zσ unbounded; CoreSpreading was
  silently inert campaign-wide; viscous equilibrium σ_eq=√(ν/Z̄)≈2.1 mm ≈
  L2's shed σ). The viscous A/B + L2 rescue are IN FLIGHT (below). Full
  story + math written to the lab notebook (`journals/20260803.md`
  §20260805) with a **TBD table that MUST be filled when the two viscous
  jobs land**.
- Das×N matrix harvested: **no d/σ collapse — Das itself owns both Γ̄
  lobes**; N=2≡N=4 null holds ONLY at 0.41c (phase_05 §2026-08-04).
- Hub 0.25R harvested: inboard aerodynamically minor (outboard Γ̄ ≤1.2%) —
  σ error is tip-driven, not root-driven.
- Phase 3: `p018_nt72_eta2` harvested 2026-08-05 — **+0.52% matched vs B0,
  FAILS M1** at production Das; dt error confirmed broad/uniform at two Das
  values; NT=144 closes the Richardson triple.
- Diagnosis synthesis + forward plan:
  `~/.claude/plans/018-per-brainstorm-index-md-has-lucky-spark.md`.
- **Phase 12 staged AND submitted (Ryan-approved 2026-08-05):** spatial rigor —
  h ladder at FIXED σ (h/σ 0.5→0.25→0.125 at B0 σ; Phase 4's ladder co-refines
  h faster than σ, so its σ claim is conditional on this axis), offline d/σ
  kernel-deficit diagnostic (Mac, script spec in phase_12), N spot-check at σ\*
  (deferred). Jobs 13051772/13051773 (table below). Note: temporal convergence
  is still NOT shown (NT=36 fails M1; nt144 pending) — Phase 3 owns it.
- **Phase 4 Ladder B added AND submitted (Ryan 2026-08-05):** σ at FIXED h
  (OVERLAP alone at PPS 4: σ/c 0.68 → 0.477 → 0.34). Range-limited by the
  overlap condition (h/σ ≤ 1) so NOT a σ→0 ladder — it is the one-factor
  complement: {Ladder B (σ|h), Phase 12A (h|σ)} decompose Ladder A's joint
  deltas. Jobs 13051774/13051775 (table below). Study slots now 10/10.
- **Phase 12 C1 MEASURED (2026-08-05, Mac, no cluster cost):** offline
  regularized-vs-singular kernel diagnostic
  (`examples/particle_surface_regularization_diag.jl`, selftests machine-
  precision) on `p018_L1_ov3`@719 → **admissible clearance d/σ\* ≈ 2.6
  (binned median) / 3.4 (p95)** at the 0.5% deficit criterion; slope −3.29 ≈
  017 prior; single-particle 1−g(d/σ) is a conservative closed-form bound.
  Detail + campaign reading (kernel deficit vs dynamic tolerance) in
  phase_12 §2026-08-05 (b); budget term 5 value filled; C2 gated on σ\*.
- **Disk/retention (Ryan 2026-08-05): a SEPARATE cleanup agent (not this
  session's) will keep the last 10 steps of every run.** Warmstart verified
  compatible with that policy (loader reads only `<run>_body1.pvd` + the
  four per-step files at the chosen index — `src/FLOWPanel_warmstart.jl:343,
  352`; prune by highest index so types stay aligned). **Sacrosanct:
  `p018_L1_ov3`@719 set (cluster + sole Mac backup `~/p018_L1_ov3_paraview/`)
  until L2_visc and chains are done.** Never prune a live run's newest steps.
- Notebook: `journals/20260803.md` §20260805 now has TWO open checklist
  items — the viscous A/B TBD table AND the spatial-rigor verdicts (fill
  when 13051772–75 land). Checkboxes are Ryan-ticked only.

### Where the campaign stands

**Phases 0, 1, 2 CLOSED (Phase 2 reopened as 2b — the κ ladder — after the
root-cause finding). Everything submitted through 2026-08-04 morning is
harvested; THIRTEEN 018 jobs are running (table below) — **Ryan explicitly
authorized the three 2026-08-05 (c) submissions (L1_visc + Green pair)
beyond the 10-cap**; treat 10 as the cap again once the queue drains, and
any further submission (`_s2` chains) waits for a completion.**

| quantity | value | source |
| --- | --- | --- |
| **Baseline B0** (η=1.0, σ/c=0.68, N=4, NT=36) | **CT̄ = 0.07170**, CI [0.07164, 0.07173], revs 15–29, M1 PASS | phase_01 |
| Das ladder {0.205, 0.41, 0.82, 1.64}c | {0.07006, 0.07187, 0.07230, 0.07084} — top rung non-monotone ⇒ **Das\* = 0.41c PROVISIONAL**; root cause = spanwise Das/c_local wedge (see above) ⇒ κ ladder in flight | phase_02 |
| dt at matched Das (0.205c) | NT 36→72: 0.07006 → **0.070473, +0.595% — FAILS M1**; ε_Γ 0.98% (marginal PASS) | phase_03 |
| **dt at production Das (0.41c)** (NEW 2026-08-05) | NT 36→72: **0.072243**, +0.52% matched (10–19) / +0.76% vs settled B0 — **FAILS M1**; ε_Γ 0.865% PASS, shift uniformly positive (no lobes) ⇒ dt error robust across Das; **NT=144 (running) closes the Richardson triple** | phase_03 |
| **σ rung L1** (σ/c 0.297) (NEW) | **CT̄ = 0.07267**, revs 17–31, window-insensitive to 0.02%; **Δ vs B0 = +1.34%** | phase_04 |
| **σ axis M2** (NEW, binding) | **ε_Γ = 8.78% max / 4.69% RMS — FAIL.** CT̄ looks nearly converged, Γ̄(r/R) is not ⇒ **L2 mandatory; σ budget term set by M2, not CT̄** | phase_04 |
| **Warm-pilot gate** (NEW) | cold 0.072411 vs warm 0.072280 on matched revs 22–31 = **−0.18% ⇒ PASS**; ε_Γ 0.29% ⇒ **cross-σ warm starts validated, L2 runs warm** | phase_04 |
| N-sensitivity (ruling 11) | N=1 (d=0.6σ) FAIL (−0.75%, ε_Γ 2.77%); N=2 (1.2σ) PASS (+0.01%, 0.35%) ⇒ N=2 provisional | phase_05 |
| **Das×N matrix** (NEW 2026-08-04) | **NO d/σ collapse — Das itself owns both Γ̄ lobes** (outboard deficit at 1.64c is N-independent; inboard deficit at 0.205c occurs at a clean d/σ). N=2≡N=4 null holds ONLY at Das\*=0.41c; 0.82c competitor weakened (own N-sensitivity −0.6%/2.2%) | phase_05 |
| **Hub 0.25R** (NEW 2026-08-04) | Γ̄ vs B0: −16% local root lobe (gone by r/R 0.41), **≤1.2% outboard of 0.45** ⇒ inboard region aerodynamically minor; σ-axis error NOT root-driven | phase_05 |
| **Hub-mesh branch** | **all three 0.15R runs VOID** — static mesh defect (fixed); hubfix diverged via cap shedding (fixed); 0.25R run healthy (above) | phase_05 |
| **L2 viscous rescue** (NEW 2026-08-04) | first attempt 13036477 FAILED (retention sweep deleted all cluster warm sources); step-719 set restored from Ryan's Mac archive; `p018_b0_visc` (13047290) + `p018_L2_visc` (13047296, warm from L1_ov3@719, 128G) RUNNING, banners verified | phase_04 |
| **d/σ clearance threshold (C1)** (NEW 2026-08-05) | kernel-level admissible **d/σ\* ≈ 2.6 (median) / 3.4 (p95)** at 0.5% deficit; slope −3.29 ≈ 017 prior; 1−g(d/σ) = conservative bound; dynamic observables tolerate ~5× more (N=2 at d/σ=1.2 passes) | phase_12 §2026-08-05 (b) |
| Campaign scatter | per-rev std ±0.0002–0.0003 (≈5× tighter than legacy) | ledger |

Cost model: **t ≈ 9.0 h/(σ/c)**, **N_p ≈ 52.4k/(σ/c)**.

**Cross-cutting insight (do not lose):** Das, σ and N are NOT independent axes
— they share the dimensionless clearance group **d/σ** (d ≈ N·Das). Raising
OVERLAP without raising PPS raises σ and undoes clearance. Any "converge one
axis at a time" plan must account for this coupling.

**Second cross-cutting insight (NEW 2026-08-03):** the two converged-so-far
axes push CT in **opposite directions** — dt refinement +0.60%, σ refinement
+1.34% on CT̄ but with a 8.8% circulation error. Do not let partial
cancellation between axes flatter the final number; the error budget must
carry them separately.

### Jobs on the cluster (snapshot 2026-08-05 ~02:25 MDT session (c) — re-check `squeue` first)

**~~FIRST ACTION~~ DONE 2026-08-05 ~16:20: `p018_L1_visc` warm-start CONFIRMED
(`simulate_warmstart!: resuming from step 719`).** Original instruction kept
for provenance: verify `p018_L1_visc` (13051801) actually
warm-started (its banner cannot show `RESTART_*`; submission at 02:12 was
too fresh for the marker, which appears a few minutes in — L2_visc's showed
at log line 24):
`ssh orc 'bash -lc "echo S; grep -i simulate_warmstart /home/rander39/projects/FLOWPanel.jl/logs/slurm/slurm-fp-018-L1visc-13051801.out"'`
— expect `simulate_warmstart!: resuming from step 719 (file count 720)`. If
it is absent once CT lines exist, the run started COLD: its SETTLE=32 plan
is then wrong (window math assumed resume at rev ~20) — kill, fix, resubmit.

**Gotcha that bit twice this session: the ssh MOTD can EAT THE FIRST LINE of
piped `squeue`/`md5sum` output** — a job can look "gone" when it is fine.
Echo a sentinel first (`echo SENTINEL; squeue ...`) or query the job id
directly before concluding anything departed. `p018_nt144` was at step ~1806
(of ~2448) and healthy at snapshot.

Check `squeue`/`sacct`; re-open access with Ryan's `! ssh orc -fN` if 2FA
blocks. **No in-session monitors survive a context clear — poll directly.**
Slot cap 10 (`fm0*` and `fp-016-*` jobs belong to other work streams).
Logs: `logs/slurm/slurm-<jobname>-<jobid>.out` in the cluster repo; stdout is
block-buffered, so judge a run by its incrementally-written
`data/<run>/monitors/*force_system1.csv`, never by log silence.

**All ten were banner-verified at submission** (mandatory rule,
`ops_reference.md`). All are on the B0 carrier (NT=36, σ/c=0.68, N=4, η=1.0)
unless a knob below says otherwise.

| job | case | phase | what it is / knobs | action on completion |
| --- | --- | --- | --- | --- |
| ~~13050924~~ **HARVESTED** | `p018_dasc0p25` | 2b | **chord-Das ladder** κ=0.25 (\|Das\|=0.25·c_local), 24 h | Harvest all three together. Judge against phase_02 §2026-08-04 (b) pre-registered predictions: successive-κ M1+M2, lobe structure vs uniform shift; κ=0.41-vs-B0 isolates the spanwise redistribution (same Das at 0.75R). **PASS ⇒ Das axis convergent, Das\* → κ\***. **FAIL ⇒ execute the pre-authorized Green fallback** (phase_02, Ryan 2026-08-04) — submit `RHPC_FORMULATION=green` κ=0.41 + 0.82 without asking |
| ~~13050925~~ **HARVESTED** | `p018_dasc0p41` | 2b | κ=0.41 (matches B0's Das at 0.75R exactly), 24 h | ↑ same |
| ~~13050926~~ **HARVESTED** | `p018_dasc0p82` | 2b | κ=0.82, 24 h | ↑ same |
| ~~13047290~~ **HARVESTED (null)** | `p018_b0_visc` | 4/Γ6 | B0 carrier + `CORE_SPREADING_ACTIVE=true` (the CoreSpreading-erratum A/B), cold, 24 h | M1+M2 vs B0 ⇒ error-budget rows Γ6 (+ CT term); **fill the TBD row in the notebook table** (`journals/20260803.md` §20260805) |
| ~~13047296~~ **FAILED (diverged)** | `p018_L2_visc` | 4 | σ-ladder L2 rescue: OV 2.88/PPS 26 ⇒ σ/c 0.151, **viscous ON**, warm from `p018_L1_ov3` step 719 (restored from Ryan's Mac), **128G**, SETTLE=50, 48 h | First question: does it survive past the inviscid ignition age (inviscid L2 died deterministically at step ~1251)? That is the direct test of the σ-blow-up mechanism. Then M2 vs L1 (pre-registered: r/R≈0.76 dip deepens, same 3-lobe shape). Expect to chain `p018_L2_visc_s2` (restart from last VTU step) before an M1 window exists. **Fill the notebook TBD row.** If it blows up again or fails M2 ⇒ 016 contingency (pre-authorized, ruling 7) |
| 13029922 → **chained as 13054739 `p018_nt144_s2`** | `p018_nt144` | 3 | NT=144, η=4.0 ⇒ Das 0.41c + σ/c 0.68 pinned, rlxf 0.075, ~17 revs, 48 h (was ~36 h in at snapshot) | Harvest; **needs a restart-chained `_s2` segment** for a ≥15-rev M1 window (recipe in `ops_reference.md`; restart from last VTU step). Then Phase 3 close-out: Richardson over {0.071866, 0.072243, CT(144)} matched 10–19 windows at Das 0.41c, report observed order → budget terms 3/Γ3 |
| 13051772 | `p018_h0p25` | 12 | **h ladder at FIXED σ** rung 2: OVERLAP 4.0 / PPS 8 ⇒ h/σ=0.25 at unchanged σ/c 0.68, merge_r 0.0120, cold, SETTLE=22 (~30 revs), 48 h | Score vs B0 on M1+M2 (matched ≥15-rev windows, `p018_analyze.py`) per phase_12 pre-registration → budget rows 16/Γ9. **FAIL ⇒ flag Ryan before touching Phase 4's σ\* logic** (its ladder co-refines h). Check the Γ̄ 3-lobe discriminator either way |
| 13051773 | `p018_h0p125` | 12 | h ladder rung 3: OVERLAP 8.0 / PPS 16 ⇒ h/σ=0.125, σ/c 0.68, merge_r 0.0120, cold, SETTLE=22, **128G**, 48 h | Confirmation/Richardson rung for the h axis; expect a `_s2` chain if walltime cuts ~1080 steps short (~4× B0 particle count). Same scoring as above vs the 0.25 rung |
| ~~13051774~~ **HARVESTED** | `p018_ov1p4` | 4B | Ladder B rung 2: σ/c 0.477 at fixed h | **+0.78% FAIL M1, ε_Γ 4.12% FAIL M2, 3-lobe at half amplitude** — phase_04 §2026-08-05 (e) |
| ~~13051775~~ **HARVESTED** | `p018_ov1p0` | 4B | Ladder B rung 3: σ/c 0.34, h/σ=1.0 (did NOT diverge) | **+1.18%/ε_Γ 7.86% FAIL both; 3-lobe ∝ σ-step ⇒ mode is σ, not h** — phase_04 §2026-08-05 (d,e) |
| 13057253 | `p018_ufront_n1` | 13 | **uniform-d_front pair** (phase_13 §7): σ=0.03R (OV 2.4/PPS 14), N=1, D=3.4 ⇒ d_front=3.4σ span-uniform, Das/c 0.38–1.41, cold, SETTLE=12, 64G, 48 h | Score n1-vs-n2 (M1+M2, matched 10–19) per phase_13 §7 pre-registration: rigid-vs-free-row A/B at admissible uniform clearance; falsifiable — inboard Γ̄ deficits must be ABSENT vs B0 or the clearance attribution falls. vs-B0 reads are context only (σ confounded) |
| 13057254 | `p018_ufront_n2` | 13 | ufront N=2 arm: D=6.5 (minimal in-band), Das_j=6.5σ−travel_j, Das/c 0.26–1.05, otherwise ↑ | ↑ same |
| 13051801 | `p018_L1_visc` | 4 | **Viscous L1** — the missing middle rung so the viscous σ ladder's L1→L2 M2 pair compares SAME physics (existing L1 is inviscid): OV 2.4/PPS 11, `CORE_SPREADING_ACTIVE=true`, **warm from `p018_L1_ov3`@719** (same source+σ caveat as L2_visc), SETTLE=32, 64G, 48 h | M1+M2 vs `p018_b0_visc` (viscous pair) AND vs inviscid L1 (isolates the viscosity delta at L1 σ → Γ6); becomes the L2_visc comparison partner. Expect `_s2` if walltime short |
| ~~13051802~~ **HARVESTED** | `p018_green_dasc0p41` | 2b/7 | Green κ=0.41 arm | **0.071921; Δ vs velocity −0.32%, M2 PASS** — phase_02 §2026-08-05 (c) |
| ~~13051803~~ **HARVESTED** | `p018_green_dasc0p82` | 2b/7 | Green κ=0.82 arm | **0.073046; green successive κ Δ +1.56% ⇒ fallback rule NEGATIVE, κ response formulation-independent, Green not adopted** — phase_02 §2026-08-05 (d) |

Also running but NOT 018's: `fp-016-conv-smooth-*` (016 discriminators,
another session) and `fm0*`.

**Warm-start caveat for `p018_L2_visc` (disclose with its result):** source is
the ov3 arm (OVERLAP 3.0, σ/c 0.292) rather than L1, and the σ jump is
0.292→0.151 (ratio 1.9 vs the validated pilot's 2.3). Also its restart set
was rebuilt from the Mac archive after the retention sweep deleted every
cluster warm source (ledger §2026-08-04; retention rule amended in
`ops_reference.md` — never delete the last complete restart set).

**Staged but NOT submitted:** the Green fallback pair (only on κ-ladder
FAIL); `p018_nt144_s2`; `p018_L2_visc_s2`; `p018_nrows8` spot-check; Phases
6–8, 10, 11 final-settings cases; Phase 9 last.

### Next actions, in priority order

1. **Harvest as jobs land** (`scripts/p018_analyze.py`; append to `ledger.md`
   + owning phase file; delete VTK per amended ruling 10 — keep one complete
   restart set per run).
2. **κ-ladder verdict** (the campaign's main fork): pass ⇒ Das axis closed at
   κ\*; fail ⇒ submit the pre-authorized Green pair (phase_02 fallback).
3. **Viscous A/B + L2 rescue**: judge, then **update the notebook TBD table**
   (`~/Dropbox/research/notebooks/journals/20260803.md`, §20260805 — Ryan
   asked for this explicitly). Chain `p018_L2_visc_s2` as needed.
4. **Phase 3 close-out** when `p018_nt144`(+`_s2`) lands: 3-rung Richardson
   at Das 0.41c (36: 0.071866, 72: 0.072243 on matched 10–19 windows).
5. Phase 4 σ fit only when viscous L2 has a settled window; judged on M2.
6. Phases 6–8, 10, 11 gated on final settings (which now await κ\*/NT\*/σ\*);
   Phase 9 last.

### Deferred until AFTER 018 completes (Ryan, 2026-08-12)

- **Port `examples/helper_functions.jl` off GeometricTools** (do NOT restore
  the dep — commit 7d52be0 removed it deliberately). Until then
  `runtests_example_suddenly_started_wing.jl` cannot load in the current
  project env and stays excluded (commented include) in the CLUSTER
  `test/runtests.jl`; ops_reference "Pre-submission gate" carries the note.
  Re-include the test when the port lands.

### Open items needing Ryan

- **`data/pitching_wing_exp/` fixture is MISSING** from disk (untracked; not
  in git, absent on cluster/Trash). `test/runtests_unit_pitching_wing_exp.jl`
  cannot run until restored via **Dropbox web → Deleted files**. Everything
  else in the suite passes.
- 25 GB ParaView set for `p018_L1_ov3` is at Ryan's `~/p018_L1_ov3_paraview/`
  (also serves as the off-cluster warm-start archive — do not delete: its
  step-719 set is the only backup of the σ-ladder warm source).

### Uncommitted source changes (all test-gated)

- `src/FLOWPanel_postprocess.jl` — isolated agglomerate (empty stencil) gets a
  zero μ gradient + warning instead of an `ArgumentError` abort. Non-isolated
  failures still throw.
- `~/Dropbox/research/projects/FLOWVPM.jl/src/FLOWVPM_subfilterscale_models.jl`
  — `Estr_direct_multithreaded` guards np=0 (the documented "small-N SFS direct
  bug"; the T4-era `SFS_OFF` workaround is now unnecessary).
- `test/runtests_unit_postprocess.jl`, `test/runtests_unit_simulate.jl`.
- `scripts/generate_dji9443_mesh.sh` — `--hub-r-over-r` (**produces a defective
  mesh — see phase_05**).
- **NEW 2026-08-03:** `examples/run_dji9443_hover_ct_hpc.slurm.sh` gains a
  `p018_L1_ov3` case arm (deployed, md5 `b2a8a2d48028c521662dd88253ef7f2f`).
- **NEW 2026-08-03 (b):** three analysis/geometry tools added —
  `scripts/p018_analyze.py` (M1/M2, validated against B0),
  `scripts/p018_plot_gamma.py` (Γ̄ ladders, reproduces M2's ε_Γ exactly),
  `scripts/p018_mesh_profile.py` (chord/twist mesh gate).
  `scripts/generate_dji9443_mesh.sh` hub recipe FIXED (see phase_05); new mesh
  `examples/data/dji9443_20260803_45_185_capped_captess4_hub0p15.msh`.
- **NEW 2026-08-04 (chord-proportional Das, all deployed + md5-verified,
  full src/ sweep matched):** `src/FLOWPanel_simulate.jl`
  (`set_Das_station_lengths` kwarg in `initialize_Das!` +
  `_set_Das_station_lengths!`; Das/velocity_te are VERTEX-based — station j =
  edge j's nib node, station n+1 = last nia), driver `DAS_CHORD_FRACTION` +
  `station_chords`, launcher `p018_dasc*` arms + `das_chord:` banner field,
  unit testset in `test/runtests_unit_simulate.jl` (7/7). Gates run:
  formulation_test 10/10, chord + stock 40_40 smokes finite. NOTE: the 016
  session also edits the launcher/driver concurrently — **always md5-verify
  local vs cluster immediately before any submission** rather than trusting
  this note.

**Standing diagnostic rule (Ryan, re-affirmed 2026-07-31): an OutOfMemoryError
thrown inside `merge_particles!` typically means the simulation blew up** —
the diverged particle cloud wrecks the merge spatial hash; memory is the
symptom, not the cause. Confirm via |CT| > 1 in the log tail (buffered stdout
truncates on SIGKILL — absence of insane CT in the log is NOT evidence of
boundedness); do not resubmit with more memory. The one recorded counterexample
is a run that was genuinely memory-capped: `p2e_sigF_nofilt` (32G request,
MaxRSS 33.4 GB, smooth CF tail) — distinguishable because CT stayed sane and
RSS approached the request.

**Standing submission rule (NEW 2026-08-03): verify every job's knobs from its
own log banner immediately after submission.** The launcher's override
precedence differs between case arms in both directions and neither warns —
this cost a 36.9 h run. Full statement in `ops_reference.md`.

## Objective and scope

Deliver a publishable convergence study of the DJI 9443 hover rotor (5400 RPM,
$C_{T,\mathrm{exp}} = 0.072$): a high-confidence converged thrust coefficient
and bound-circulation distribution $\Gamma(r/R)$, with **each discretization
axis converged separately** and an explicit error budget. Context: CT currently
spans 0.034–0.110 across plausible discretizations, so no single number is yet
a prediction. This item runs an HPC campaign on existing machinery
(`examples/rotor_hover_pressure_comparison.jl` +
`examples/run_dji9443_hover_ct_hpc.slurm.sh`, `p018_*` case tags) and reuses
recorded Phase-2e/006/014 results as data points.

All operational detail lives in the subdirectory
`018_dji9443_hover_convergence_campaign/`. A clean-context agent should read,
in order: this file → `ops_reference.md` → `decision_rules.md` → the phase file
it is executing. `ledger.md` is the single running results table.

## Standing rulings (binding on every phase)

1. **RPM = 5400 for all convergence jobs** (Ryan 2026-07-30).
2. **Pedrizzetti relaxation FILTER stays OFF in every run** — do not explore it
   in this study; recommending it in the final report is allowed if evidence
   suggests it would help (Ryan 2026-07-30). Stock corrected-Pedrizzetti
   relaxation itself stays ON with stock `rlxf` (Ryan 2026-07-28); relaxation is
   carried as an error-budget term (≈ −0.005 CT, monotone in rlxf, unconverged).
3. **Do not tune to the target**: $C_{T,\mathrm{exp}}$ is compared only after
   the error budget closes; it never selects settings (Ryan 2026-07-28).
4. **No nwakerows ladder** — the legacy-σ ladder was non-monotone and rejected
   (Ryan 2026-07-30, BRAINSTORM/014). Instead: choose `NWAKEROWS` high enough
   that particle-on-body influence isn't stunted by smoothing-radius overlap
   (handoff distance $d \gtrsim 4\sigma$; N=4 here), plus one N=8 spot-check.
5. **Memory: flat 64G, 64 cores per job, ≤6 concurrent** (Ryan 2026-07-30).
   **Amended 2026-08-04 (Ryan): concurrency cap raised to ≤10.**
   **Amended 2026-08-03 (Ryan): σ-ladder L2 and finer may request 128G.**
   `p018_L2` was OOM-killed at 64G while running cleanly (bounded CT), at
   ~347k particles; this is the same escalation that produced the 64G ruling
   itself after `p2e_sigF_nofilt` died at 32G. **Corollary for triage: judge
   |CT| from the force-monitor CSV, not the stdout log** (block-buffering
   truncates the log on SIGKILL), and **do not** use `MaxRSS ≈ ReqMem` as a
   memory-ceiling test — L2's MaxRSS was 35% of its request because `sacct`
   sampling missed the spike.
6. **Convergence metric: long averaging window; the key requirement is that the
   mean settles** (Ryan 2026-07-30). The old 5-rev gate is retired with a
   limit-cycle defense; see `decision_rules.md`.
7. **BRAINSTORM/016 (surface-vorticity shedding) is PRE-AUTHORIZED** as the
   contingency if particle properties cannot be converged with existing knobs
   (Ryan 2026-07-30). Try existing machinery first.
8. `DAS_REFRESH` stays **false** everywhere (`Backslash` LU caching bug);
   `DAS_KINEMATIC_ARC=false` (tangent) for comparability with the recorded
   ladders; PressureBernoulli(unsteady=false) caveat disclosed in publication.
9. **SFS (DynamicSFS) and viscous (CoreSpreading) particle models stay ON in
   every campaign run** (Ryan 2026-07-31; driver defaults already comply —
   `SFS_OFF` is permitted only in backend-A/B diagnostics like T4, never in a
   ladder rung). **Mesh spanwise and chordwise refinement are convergence axes
   in their own right** (Phases 10/11, respectively at fixed chord and fixed
   span). **The convergence observable is CT̄ AND the spanwise circulation
   distribution Γ̄(r/R), determined from the trailing-edge μ jump**
   (BoundCirculationMonitor `circulation_te`; slice estimator as cross-check)
   — every axis must pass BOTH (Ryan 2026-07-31).
10. **HPC disk retention (Ryan 2026-07-31):** after a run is harvested and
    deemed no longer useful, **delete its VTK files, keeping the CSV
    summaries** (`*_CT_vs_rev.csv`, `*_CT_per_rev.csv`, `*.toml`,
    `monitors/*.csv`). **Keep VTK for at most 3 runs at a time** (prefer the
    actively-running and most decision-relevant). If ParaView inspection of a
    run would help answer a convergence-behavior question, **ask Ryan to
    examine it** rather than hoarding histories.
11. **nwakerows is not a convergence axis (Ryan 2026-07-31, superseding the
    framing of ruling 4):** "convergence of nwakerows doesn't even make sense"
    — the quantities to converge are **particle shedding (σ, spacing) and
    Das**. Direction: try N=1 (and N=2) at campaign settings; treat N as a
    modeling choice whose sensitivity is measured and disclosed, not driven
    to a limit. (Ruling 4's d≳4σ criterion and the N=8 spot-check stand as
    sensitivity probes, not as a convergence claim.)

12. **N preference revised (Ryan 2026-08-03), superseding "minimal N":**
    **prefer N=2 over N=1 if N=2 permits a significantly smaller number of
    particles shed per step.** N sets the handoff distance d ≈ N·Das and the
    measured failure mode is clearance in d/σ, so buying clearance with N
    rather than σ permits a coarser admissible σ — and since
    σ = 2πR·OVERLAP/(NT·PPS), coarser σ means smaller `P_PER_STEP`, i.e.
    directly fewer particles and a much cheaper run. Re-frames the deferred
    N=1-at-L1 test: the question is whether N=1 buys anything N=2 does not,
    given N=2 lets us shed fewer particles.

## Fixed operating point and unit conversions

Mesh `45_185_ct4` (36,752 panels), $R = 0.119$ m, $\rho = 1.179$,
$c(0.75R) \approx 0.128R$; staged startup (SPINUP 1.5 revs @ 0.4 start
fraction, freestream peak 5 m/s, ramp/hold/withdraw 1.0/1.5/4 revs, SETTLE ≥ 12);
velocity formulation unless stated; `PARTICLE_SHEDDING=sigma_overlap`;
truncation 4R (radius 1.5R hard-coded — disclosed, not converged);
`NWAKEROWS=4`. Conversions (freeze factor 0.4 included; tangent Das):

$$\mathrm{Das}/c(0.75R) = 0.41\,\eta\,(36/\mathrm{NT}), \qquad
\sigma/R = \frac{2\pi}{\mathrm{NT}}\cdot\frac{\mathrm{OVERLAP}}{\mathrm{PPS}}, \qquad
\sigma/c(0.75R) = 1.363\,\frac{36}{\mathrm{NT}}\cdot\frac{\mathrm{OVERLAP}}{\mathrm{PPS}}$$

$h = \sigma/\mathrm{OVERLAP}$; classical convergence needs $h$ refined faster
than $\sigma$ ($\sigma^q/h$ const, $q>1$) ⇒ OVERLAP increases down the σ
ladder. `MERGE_R_FACTOR` $= 0.138\,(\sigma/R)$ holds $r_\mathrm{merge}/\sigma$
at the stock ratio. **Pinning σ across NT** (at fixed OVERLAP): halve
`P_PER_STEP` when NT doubles (σ ∝ OVERLAP/(NT·PPS)); the `SigmaOverlap` policy
then sheds the correct per-step count on its own.

## Phase gates

| Phase | Deliverable | Status |
| --- | --- | --- |
| [0](018_dji9443_hover_convergence_campaign/phase_00_harvest.md) | Harvest 3 in-flight cluster jobs; seed ledger; deploy scaffolding | **CLOSED 2026-07-31** — all 3 harvested; nt72match ⇒ Phase 3 full ladder incl. NT=144 |
| [1](018_dji9443_hover_convergence_campaign/phase_01_baseline_settling.md) | Baseline B0 (η=1.0, σ/c=0.68) + settling demo + restart validation + cost recalibration | **CLOSED 2026-07-31** — B0 = 0.07170 CI [0.07164,0.07173] (revs 15–29, M1 PASS); restart chaining validated (seam within natural variability); cost model recalibrated |
| [2](018_dji9443_hover_convergence_campaign/phase_02_das_length.md) | Physical Das ladder {0.205, 0.41, 0.82, 1.64}·c → Das\* = plateau lower edge; FMM/Direct discriminator | **RESOLVED 2026-07-31 (no-plateau branch): Das\* = 0.41c PROVISIONAL** — no pair passes M1+M2; top rung non-monotone (−2.02%); handoff confound flagged ⇒ Phase 5 decision-critical; T4 PASSED |
| [3](018_dji9443_hover_convergence_campaign/phase_03_timestep.md) | NT ladder at pinned physical Das AND σ (36 → 72 [→ 144]) → NT\* | **Matched-Das rung MEASURED 2026-08-03: NT 36→72 at Das 0.205c = +0.595%, FAILS M1 ⇒ dt not converged.** Corrected `p018_nt72_eta2` (13029921) + `p018_nt144` (13029922) RUNNING at Das=0.41c |
| [4](018_dji9443_hover_convergence_campaign/phase_04_sigma_ladder.md) | σ ladder 0.68 → 0.297 → 0.151 c with q ≈ 1.2 → production σ + extrapolation term; **+ Ladder B (2026-08-05): σ at fixed h, OVERLAP 2.0 → 1.4 → 1.0 at PPS 4 (range-limited one-factor complement)** | **L1 SETTLED: 0.07267 (+1.34% vs B0) but M2 FAILS at ε_Γ 8.78% ⇒ L2 mandatory. Warm-pilot gate PASSES (−0.18%) ⇒ L2 runs warm.** Viscous rescue in flight; Ladder B rungs 13051774/13051775 RUNNING |
| [5](018_dji9443_hover_convergence_campaign/phase_05_shedding_distance.md) | N-sensitivity + **staged Das×N×σ interaction study seeking convergent clearance criteria** (ruling 11; Ryan directive 2026-07-31 — see phase file's STAGED section) | N=1 FAIL / N=2 PASS recorded. **Hub branch VOID 2026-08-03 (mesh defect: CT −0.474 at step 1) — parked pending diagnosis, no slots committed.** Coarse-σ Das×N matrix still unrun |
| [6](018_dji9443_hover_convergence_campaign/phase_06_truncation.md) | 6R vs 4R confirmation at final settings | OPEN |
| [7](018_dji9443_hover_convergence_campaign/phase_07_green_delta.md) | GreenReconstruction ΔCT and ΔΓ(r/R) at converged settings | OPEN |
| [8](018_dji9443_hover_convergence_campaign/phase_08_merging.md) | Merging null demonstration + speedup report | OPEN |
| [9](018_dji9443_hover_convergence_campaign/phase_09_final_error_budget.md) | ≥30-settled-rev production run; CT̄ ± CI, Γ̄(r/R) ± CI; closed error budget | OPEN (last) |
| [10](018_dji9443_hover_convergence_campaign/phase_10_mesh_spanwise.md) | Mesh spanwise ladder 45 → 60 [→ 80] at n_airfoil=185, final wake settings (needs local mesh generation) | OPEN — staged 2026-07-31 |
| [11](018_dji9443_hover_convergence_campaign/phase_11_mesh_chordwise.md) | Mesh chordwise ladder 145 → 185 → 249 at span 45, final wake settings (meshes exist) | OPEN — staged 2026-07-31 |
| [12](018_dji9443_hover_convergence_campaign/phase_12_spatial_rigor.md) | Spatial rigor: h ladder at FIXED σ (h/σ 0.5 → 0.25 → 0.125 at B0 σ, M1+M2) + offline d/σ kernel-deficit threshold + N spot-check at σ\* + σ-axis conditionality note | OPEN — staged 2026-08-05; h rungs SUBMITTED (cold, SETTLE=22); C2 gated on σ\* |

Guessed parameter values are starting points — executing agents may adapt
based on results, recording deviations in the phase file and `ledger.md`.
Dependencies: P0 → all; P2 needs the B0 carrier definition (can run alongside
P1); P3 needs P0+P2; P4 needs P2; P5 needs P1; P6/P7/P8/P10/P11 need P4
(final settings); P9 last, after 10/11 fill their budget terms.
Respect the 6-job cap throughout.

## Contingency chain (Phase 4 failure path, in order)

(a) raise OVERLAP; (b) `SHEDDING_R_OVER_R` root-clip test (existing knob —
cheap proxy for a larger hub: near-hub particles barely move; 40_40 hub
exclusion measured ≤ +2e-4 CT); (c) hub-radius mesh variant — **IMPLEMENTED 2026-07-31**
(`scripts/generate_dji9443_mesh.sh --hub-r-over-r X`; recipe cuts interior
XSecs inside the hub then moves the stock root out; validated ceiling ≈ 0.15R,
larger hubs fold the loft non-watertight; hub0p15 mesh md5-verified on
cluster; 3 comparison runs staged in phase_05);
(d) **016 implementation (pre-authorized)** per its phase_03 architecture doc.

## Error budget (assembled in Phase 9)

Skeleton in `018_dji9443_hover_convergence_campaign/error_budget.md`: Das
model-form (~0.2%/doubling) · dt residual · σ extrapolation · N spot-check
spread · truncation null · Green Δ · merge Δ · relaxation (≈ −0.005 CT,
unconverged) · PressureBernoulli(unsteady=false) disclosure · Dirichlet
tangency contamination (2d App. G, ~0.8–1% of CT at n=185) · radial
truncation 1.5R fixed.

## Decision log

- 2026-08-04 — Ryan: slot cap raised to 10. Chord-proportional Das
  (`DAS_CHORD_FRACTION`) implemented and κ ladder submitted (Ryan: "yes, try
  it"). **Standing authorization (Ryan): if the κ ladder does not converge the
  Das axis, the executing agent may try the next most likely approach without
  further approval — pre-registered in phase_02 as the GreenReconstruction
  Das pair (015 Battery II rationale).**

- 2026-07-30 — Item opened; scaffolding written; rulings 1–8 recorded (planning
  session with Ryan). Correction vs the planning doc: pinning σ across NT
  requires **halving** `P_PER_STEP` as NT doubles, not doubling it
  ($\sigma = 2\pi R\,\mathrm{OVERLAP}/(\mathrm{NT}\cdot\mathrm{PPS})$).
- 2026-07-31 — Ruling 9 added (Ryan): mesh spanwise + chordwise refinement
  staged as Phases 10/11; SFS+viscous stay ON in all runs; Γ̄(r/R) from the
  TE μ jump is a co-equal convergence observable with CT̄ on every axis.
  Launcher gained `p018_chord145/chord249/span60/span80`; span meshes need
  local generation (Phase 10 prerequisite).

## STATUS REPORT for Ryan (2026-08-14 evening)

Scheduled-resume session (brief (j)), working from disk. Every decision +
result below, in order; details in `phase_16_chord_sigma_coscaling.md` §Log
and ledger §2026-08-14.

1. **Phase 16 screens (13170886 l2p4 / 13170887 l3p4): both COMPLETED
   (~9.5/10 h, 414 steps = 11.47 revs each) and both PASS the stability
   gate** — decision: PASS, per the pre-registered tripwire judged on the
   full-tail evidence:
   - No ignition signature at any step: max_u ≤ 41→24 (l2p4) / 53→33 (l3p4)
     m/s (ignitions hit 1e4–1e6); min_sigma ends 1.87e-3/1.79e-3 m — 19×
     above the 9.4e-5 viscous floor (floor-pinning is the ignition mark);
     no σ ≤ 0 anywhere.
   - Calibration vs the healthy uniform carrier `ufront_n1_s040_visc` at
     MATCHED age (cold revs 0–12, which went on to survive 60 revs):
     carrier rev-11 max|Γ|/σ² = 1.29e3 and min_sigma 1.09e-3; the screens
     show |Γ|/σ² 72 (l2p4, already declining) and 297 (l3p4, rising but
     decelerating) — 4–18× BELOW the healthy reference; min_sigma higher
     (healthier) than carrier. The min_sigma drift below the smallest shed
     σ is present in the uniform carrier too ⇒ normal code behavior, not a
     station-law artifact.
   - n_p budget tracks the carrier's trajectory and saturation shape
     (merge+truncation outflow normal); no tip-fusion signature (end n_p
     202.9k/215.2k vs carrier 213.2k at rev 11).
   - **FLAG (instrumentation gap, non-blocking):** `min_sigma_ratio` AND
     `p1_sigma_ratio` are NaN for all 414 steps in both screens — the
     wake-health monitor's reference σ is NaN under the station-σ law.
     Tripwire judged from min_sigma/max_u absolutes + carrier calibration
     instead. Worth a small monitor fix before any Phase-16 run we intend
     to watch in-flight by sigma-ratio; NOT redeployed tonight (deploy risk
     vs pre-authorized ladder).
   - Screen CT (disclosed, NOT judged — screens end inside the withdraw
     transient): l2p4 rev-11 block 0.0651, l3p4 0.0863 and rising.
2. **Decision: submitted the pre-authorized 30-rev ladder** (screens-pass
   condition met): **13178762 `p018_cs_l2p4`, 13178763 `p018_cs_l3p4`,
   13178764 `p018_cs_l4p8`**, 48 h/64G, `--export=ALL,P018_SETTLE_REVS=22`,
   submitted ~21:15 MDT, PENDING at submission. Banner verification on
   start: see below (this session).
3. **Decision: swept both screens' VTK** (harvested to local
   `data/p018_scr_cs_l{2p4,3p4}/` first — monitors + CT CSVs + metadata;
   sweeper kept newest-10 restartable steps): freed 19.2G; cluster data/
   123G vs 200G budget.
4. Brief-(h) jobs at session start: rlxf-up pair + etadas_n4 already
   landed+scored (prior session, ledger); 13157833/834 (upin nt72/nt144)
   still RUNNING, ~26 h left, land ~08-16 — burst attribution runs then
   (next session; your M1-observable decision stays gated on it).
5. **Ladder banners VERIFIED running** (~23:45 MDT): 13178762/63/64 show
   das_lambda 2.4/3.4/4.8, sigma_chord 0.313, sigma_floor 0, settle 22,
   carrier line identical to screens, θ_max 1.26/1.79/2.53 rad, cold
   start. No mismatch, nothing cancelled. ETA ~08-16.
6. **trunc6 (13157751) TIMEOUT at 48 h — but at step 1062/1080 (rev
   29.5), so the matched window survived. H2 first read: NOT KILLED —
   drift is depth-dependent.** Matched 15–29: CT̄ 0.074411 and drift
   −3.56% monotone DOWN vs carrier 0.075298 / +2.08% UP; late 22–29
   level effect −2.7%; M2 ε_max 1.87%. Disclosed confound: the 5.5R
   domain was still filling through the whole window (n_p rising
   6–9k/rev vs carrier saturation ~rev 12), so the DOWN drift mixes fill
   transient with depth response. Clean claims tonight: the carrier's
   upward drift is not depth-independent, and deeper truncation lowers
   CT. Detail in phase_15 (c).
7. **Judgment call (mine, logged for your review): submitted
   `p018_trunc6_n1_s040_visc_s2` = 13179250** (warm restart @1061,
   SETTLE 36 ⇒ extends to rev 44, 72 h/64G) to give H2 a settled 30–44
   window vs the carrier's `_s2` on the same revs. Rationale: the
   authorized conditional (c) as-run cannot separate fill transient from
   depth response, truncation is a named error-budget term, and the
   chain recipe is the campaign's routine pattern. It is one 72 h job —
   scancel 13179250 if you disagree; the cold-run harvest and tonight's
   read stand either way. Will banner+warmstart-verify when it starts.

### 2026-08-15 (evening) — two ignitions; dt ladder stopped

8. **trunc6_s2 (13179250) IGNITED at rev ~34.5** (exit 0 with garbage —
   caught by wake-health, not exit code). Not a seam artifact: the cold
   5.5R run had already left the carrier's health envelope BEFORE the
   handoff (|Γ|/σ² 5.3e3 vs carrier 1.6e3 at rev 29). **Finding: the
   3.5R truncation deletion is load-bearing for stability** — it removes
   the aged wake tail before Γ-ignition matures (joins stock rlxf and
   σ=0.04R in the 019 stabilizer set). H2's settled-window
   quantification is blocked at 6R; the cold-run read (item 6) is the
   H2 result of record. **Your options (staged, not enacted): (a) a
   4.5R arm, (b) accept the cold-run read and carry truncation as a
   one-sided level term −1.2..−2.7% (fill-confounded).** Corpse swept.
9. **upin_nt72 (13157833) FAILED at rev 19.6 — Γ-ignition** (min_sr
   collapse to 0.014 at rev 15, transient arrest, re-ignition; min_sigma
   pinned at the NT=72 viscous floor; wake self-annihilated; forces blew
   up). Since its rlxf 0.15 is dt-equivalent to stock, this reads on
   **dt refinement itself** at matched physical relaxation, σ, and
   clearance — refinement tightens the ignition boundary, as 014/019
   anticipated. **Pre-registered screen rule applied: dt ladder
   STOPPED.** nt144 (healthy at rev 14) is being allowed to hit its
   72 h wall (~rev 14.6 tonight) and will be harvested, but I am NOT
   submitting its `_s2` chain — the chain instruction assumed a passing
   screen, and the stop rule says bring you the 019 regime map before
   resubmitting anything. Richardson triple is dead at these knobs.
10. **Burst attribution is now PARTIAL-ONLY** (no rung reaches the
   matched 15–29 window): on shared clean revs 7–13, within-rev CT std
   NT36→NT72 = 0.00303→0.00199 (×0.66 — moves with Δt, leans
   NUMERICAL), but both rungs burst co-timed at rev 14 with similar
   size (that episode looks PHYSICAL); CT̄ differs −0.31%. nt144's
   harvest completes the 3-point trend tonight. **Your M1-observable
   decision remains gated on this; the discriminator now only covers
   the pre-settled era, which weakens whichever verdict emerges —
   factor that into how much weight the raw mean vs quiet-limit
   choice puts on it.**
11. Ops: swept trunc6_s2/upin_nt72/trunc6-cold corpses (~62G) after a
   226G disk peak; data/ at 166G. Phase-16 ladder healthy at revs 17–19
   (|Γ|/σ² 202/338/1800, within the carrier envelope; l4p8 warmest —
   under watch with absolute-value tripwires).

### 2026-08-15 — minimal-run plan: adopted parts ENACTED (Ryan-directed)

Ryan reviewed an external agent's "Minimal Additional-Run Plan" and my
assessment, then directed: "do the adopt-as-is parts now." Adopted and
now in force (the contested/pending parts are NOT adopted — see below):

12. **ADOPTED — Phase-16 close-out protocol** (= the existing
   pre-registration, reaffirmed): score the running λ ladder on matched
   windows with M1 (BOTH raw-mean and quiet-limit), M2 radial loading,
   burst/quiet decomposition, and wake-health gating; P1 doubling test
   λ 2.4→4.8; P3 A/B (cs_l3p4 vs ufront_n1_s040_visc). **No
   intermediate-λ rungs will be added to tune CT.** On failure, route
   by the registered signature (inboard/θ-correlated ⇒ F1 curvature
   cap; small-λ handoff ⇒ F2 Green pair); no scalar-ladder
   densification unless genuinely ambiguous.
13. **ADOPTED — selected-rung extension rule:** if the winning Phase-16
   rung passes settledness at 30 revs, no extension; otherwise
   warm-chain THAT ONE RUNG to 45–60 revs (`_s2` recipe), and it
   becomes the primary hover-validation case. If burst rectification
   persists, report raw mean AND the pre-registered quiet-limit
   side-by-side. (Execution follows the ladder harvest, ~08-16.)
14. **ADOPTED — do-not-run list (standing):** no rlxf 0.15/0.075
   reruns; no additional upward-rlxf rungs (0.3/0.45/0.6 already
   CT-flat within burst noise); no nt144 `_s2` chain (screen-fail stop
   reaffirmed); no new truncation-depth or merge ladders for the
   CT-validation objective (reinforced by the trunc6_s2 ignition); new
   case tags for any model-definition change — never overwrite or
   re-tag historical runs.
15. **ADOPTED — acceptance/validation language:** numerical settings
   selected blind to experimental CT; validation-grade runs require
   clean wake health + the settledness gate; report CT with CI,
   averaging window, convention, and Δ vs CT_exp = 0.072; if temporal
   or spanwise checks fail, the claim downgrades to "single-point CT
   agreement with unresolved numerical uncertainty," not converged
   validation.
16. **NOT adopted (pending your ruling):** plan items 3–5 — the
   exact-rate NT=72 run (rlxf 0.16334; my assessment: outcome nearly
   predetermined, 9% dose vs factor-2 sensitivity), NT=72 @ stock 0.3
   (my recommended swap as the only viable temporal probe;
   model-definition change), and the spanwise 45→60 rung (blocked on
   verifying the plan's "80-station contraction problems" claim
   against phase_10 — it conflicts with the recorded 40/56/80
   non-monotonicity where 40 & 80 AGREE). Any dt-axis submission also
   still requires your 019-regime-map sign-off per the stop rule.
17. Watch note (~19:30 MDT): **cs_l4p8 |Γ|/σ² rising fast** — 1.8e3 →
   6.7e3 over revs ~17→18.8 (above the carrier envelope; min_sigma
   5.3e-4 and max_u 19 still healthy). Consistent with the λ=4.8
   curvature-edge rung being the stressed one (θ_max 2.53 rad). If it
   ignites, that is a P2-relevant failure signature (θ-correlated ⇒
   F1), not a screen bust — the λ=2.4/3.4 screens carried the
   stability gate. Under tripwire watch.
18. **nt144 harvested (TIMEOUT rev 14.45, healthy to the end) and the
   burst attribution is now 3-point: the fluctuation mode is
   NUMERICAL.** Within-rev CT std on shared revs 7–13: 0.00303 →
   0.00199 → 0.00035 (NT 36/72/144) — monotone, ×0.12 overall — and
   the rev-14 co-timed burst vanishes at NT=144. Caveats: pre-settled
   window only; CT̄ is NON-monotone in NT (0.07513/0.07490/0.07775,
   +3.5% at 144), so the mean-CT dt axis is NOT resolved by this.
   Staged implication for your M1 call: bursts-as-numerical-
   rectification supports the quiet-limit observable.
19. **Your mid-session direction enacted:** launched both NT=72 arms
   as explicitly-labeled MODEL-DEFINITION changes (not dt-convergence
   rungs): **13183888 `p018_upin_nt72_rlxf0p3`** (stock per-step rate)
   and **13183889 `p018_upin_nt72_rlxf0p16334`** (exact
   continuous-rate equivalent). New launcher arms + case tags per the
   adopted provenance rule; deployed md5-verified (4135f810);
   72 h/64G, SETTLE 22 (30 revs), cold starts; banner verification on
   start. Scoring: stability vs nt72's rev-14.5/18.6 marks first, then
   matched windows vs the NT=36 carrier, always labeled model-def A/Bs.
20. Adopted-plan items 12–15 above are ENACTED (recorded in phase_16
   §Post-ladder rule, phase_03, and the ledger); items in 16 remain
   yours except the two NT=72 arms you just directed (now launched).
21a. **Weekend pause/resume (your direction) both executed on
   schedule:** all 7 jobs killed 23:52 Sat (manifest
   `pause_manifest_20260815.md`; queue verified empty; cs_l4p8's hot
   phase SUBSIDED before the kill — no ignition), and resumed 00:01
   Mon as warm chains **13185005–010** (`*_rs1`, restart steps
   834/778/713/252/194/268, SETTLE 22) + **13185011 fm041aH**
   (other-session job, from scratch). Banner + warmstart verification
   on start; monitor re-armed. nt144 arm still carries the
   ask-before-chain note.
22. **Phase-16 interim (l2p4 + l3p4 landed 08-17, healthy to rev 30,
   seams clean): co-scaling COLLAPSES the NT=36 burst mode.** Per-rev
   CT std ×7–10 smaller than the uniform carrier; **λ2.4 = the
   campaign's first raw-M1 PASS** (drift 0.18% non-monotone); λ3.4
   CHECK (0.52%, non-monotone). CT̄ 15–29: 0.071074 / 0.072442 ⇒
   λ-slope 2.4→3.4 = +1.92% (P1's 2.4→4.8 doubling verdict waits on
   l4p8, ETA Mon evening). **M2 FAILs so far**: rung pair ε_max 5.97%;
   P3 A/B vs uniform carrier 6.45% (bigger than the predicted 1–4%
   redistribution). Full P1/P2 verdicts + Γ̄ lobe analysis on l4p8's
   landing; no fallback selection without you.
27. **TE-probe 2×2 landed; your picks (TE location, cs_l3p4 source,
   clamp OFF) enacted; csarc screens + your midpoint sanity arm
   SUBMITTED: 13206336/337 (screens, 24 h) + 13206338 (mid A/B,
   30 revs).** Key probe findings: mid-span wake-only downwash ≈ 0.07
   tip, locations agree 3–14% mid-span; **the rev-58 `_s3` wake shows
   inboard fountain flow the rev-30 cs wake lacks** (both locations —
   real wake-maturity physics); the cs table's tip upwash is physical,
   so the clamp was NOT free (your question — answered before launch).
   Launch-blocking fixes: driver table reader (comment lines, row
   order). TE ladder gated on the screens. You killed the local smoke
   at rev 3.1 (resource pressure) — record clean; screens carry the
   stability gate.
25. **F1 trio failed at load (env, not code) and was resubmitted as
   13193513/514/515.** Cause: the cluster Manifest.toml was
   re-resolved under Julia 1.12.6 (Aug 17 evening, the 021 session's
   window); our 1.11.7 jobs hard-fail on 1.12 manifests. Fixed by
   snapshot-restore of the 1.11.7 manifest + conservative resolve
   (021's ILUZero kept; their 1.12 manifest saved aside).
   **Cross-session flag: a 021 resubmit under Julia 1.12 can re-break
   this — worth a word with that session's workflow.**
26. **F1b (endpoint-on-arc Das) implemented per the approved plan and
   your in-flight rulings** (Route A rejected; arc-length-integrated
   endpoint, no floors; general frame-walked construction; live mode
   deferred on the Kutta operator cache; axial-only clamp with
   generalization deferred on your go-ahead). src kernel + driver
   knobs deployed (simulate.jl proven cluster-base+my-hunks before
   deploy); 14/14 new unit tests pass, legacy tests untouched;
   kinematic-mode smoke running locally. **Probe measurement finding:
   a point probe at the TE is dominated by the body's own bound-sheet
   field (0.5–0.9 tip speed!) — the table therefore uses WAKE-only
   induced + Uinf drift, with all-sources kept as diagnostic
   columns.** The 2×2 probe job (cs_l3p4_rs1 / _s3 ends × te / mid)
   is 13206092; when it lands you get the comparison and pick the
   table before any csarc screens go up.
24. **F1 LAUNCHED on your order (Tue): 13193493/494/495 =
   `p018_cs_f1_l{2p4,3p4,4p8}`**, β = 0.6 rad via new driver knob
   `DAS_CURVATURE_BETA` (driver-only change, no src/). Gate smoke ran
   9 revs clean at λ4.8+cap; banner-capture verified cap ACTIVE (binds
   15/36 coarse stations, θ_max = 0.6 exactly). Deployed md5-verified;
   48 h/64G, SETTLE 22. Disclosed tradeoff: the cap pushes inboard
   clearance below the C1 band (λ_eff ≈ 1.9 at r/R 0.27) — intrinsic
   to adding the curvature band; watching M2 inboard. Scoring plan on
   landing: doubling-under-cap (M1 raw+quiet, M2), Γ̄ overlay vs
   uncapped rungs, and l3p4 cap-only A/B.
23. **PHASE-16 LADDER COMPLETE + SCORED (Mon evening). P1 FAIL, P2
   selects the curvature fallback, P3 shows the σ-law is the big
   lever.** Doubling 2.4→4.8: ΔCT̄ +9.63% (raw AND quiet agree), M2
   31.4%. The λ4.8 deviation is an INBOARD, θ-correlated, λ-monotone
   Γ̄ EXCESS (+46.7% of Γmax at r/R 0.27 → ~0 by 0.76; the
   pre-registration predicted a deficit — sign flipped, selector
   unchanged). λ3.4-vs-uniform A/B: ΔCT −3.79%, M2 6.45%, and
   co-scaling collapses the burst mode (λ2.4 = first raw-M1 PASS).
   **Per the pre-registered routing this selects F1 (curvature-capped
   |Das|_j = min(λ·s*·c_j, β·r_j), evidence β ≈ cap near θ 0.5–0.7
   rad) — STAGED FOR YOUR CONFIRMATION, nothing submitted.** Full
   numbers in phase_16 §Log; Γ̄ tables reproducible from local
   harvests.
21. **Third model-def arm (your direction): 13183998
   `p018_upin_nt144_rlxf0p3`** (NT=144/PPS=3, stock per-step rlxf 0.3;
   72 h reaches ~rev 14–15). NT=72 arms' banners VERIFIED
   (rlxf:0.3/0.16334, pps:6, das_uniform:3.4, cold). **STANDING
   NOTE-TO-ASK (your instruction): when 13183998 finishes, ask you
   whether to chain `_s2` — likely yes if the NT=72 rlxf0p3 arm
   (13183888) looks good. No chain without your answer.**

## STATUS REPORT for Ryan (2026-08-20)

Fresh-agent resume off brief (l). All nine in-flight jobs are accounted
for; details in ledger §2026-08-20 and phase_16 §Log 2026-08-20.

28. **F1 trio landed (48 h wall, revs ~28.5 of 30) and is scored — the
   curvature cap WORKS on the failure it targeted.** Banners verified
   post-hoc first (cap ACTIVE, θ_max 0.6, binds 10/13/17 of 41). On
   matched 15–28 (uncapped re-scored on the same window): CT̄
   0.070539/0.070605/0.069877 (λ2.4/3.4/4.8) ⇒ doubling −0.94% vs
   +9.68% uncapped; λ4.8 burstiness gone (std ×10 down); Γ̄ inboard
   excess −0.6..−0.9% of Γmax (was +34..+54). M2 doubling 5.63% —
   improved ×5.6 but above the 1% gate; the residual is an OUTBOARD
   lobe (r/R 0.78–0.87, ~−5.6%) that the uncapped ladder also has —
   cap-untouched, a separate λ-sensitivity in the tip-roll-up region.
   Cap-only A/B at λ3.4: ΔCT̄ −2.56%, Γ̄ change confined inboard of
   0.4. **No chains submitted; nothing further launched for F1.**
29. **F1b screens PASS ⇒ TE ladder submitted (your pre-authorization):
   13243083/084/085 `p018_csarc_l{2p4,3p4,4p8}`, 48 h/64G, SETTLE 22,
   banners verified running** (steady TE table cs_l3p4_rs1, no cap).
   Screens were clean on the pre-registered absolutes (zero σ≤0, max_u
   ≤23, min_sigma 1.78e-3, gos2 below carrier calibration). Lands
   ~08-22 morning. csarc_mid sanity arm lands tonight; its A/B vs
   `p018_csarc_l3p4` runs when both exist.
30. **nt72 model-def arms: both SURVIVED to revs 28.2/26.6 (walls) —
   the earlier nt72 ignition was the low-relaxation dose, not NT.**
   Matched 15–26 vs the NT36 carrier (0.075155): stock 0.3 → 0.077198
   (+2.72%) but bursty/drifting (std 3.2e-3, 6.6% monotone block
   drift); exact-continuous 0.16334 → 0.076717 (+2.08%), well-behaved
   (std 1.2e-3, drift 0.48%). M2 vs carrier 4.34%/2.93%. Dose within
   NT=72: −0.62%. **Interpretation and any further NT work are yours;
   nothing chained.**
31. **nt144@0.3 (13185010) DIED at rev 13.07 — Bus error (node fault,
   not physics: health calm at death).** Newest-10 restart set through
   step 1881 is retained. **YOUR CALL (folds into your standing nt144
   ask): (a) warm-resubmit from 1881 to finish the 30-rev arm
   (~20 h+), (b) drop it (the burst-attribution verdict already stands
   on the old nt144 run), or (c) drop and spend the slot elsewhere.
   My recommendation: (a) only if you still want the NT=144 model-def
   CT̄ point; the stability question it was probing is already
   answered by the nt72 arms.**
32. **Ops:** VTK sweeps freed 46.9G (data 165→120G of 200G); harvests
   local for all seven landed runs; TE-ladder take-1 submission died
   harmlessly to a quoting bug and was resubmitted clean (ledger);
   Manifest untouched at 1.11.7; sacct monitor re-armed on the four
   live jobs; no notebook writes.
