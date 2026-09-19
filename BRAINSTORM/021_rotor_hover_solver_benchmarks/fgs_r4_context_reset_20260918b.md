# v22 chunked A/B: context-reset handoff (2026-09-18)

Supersedes `fgs_r4_context_reset_20260918.md` — its task is DONE through §8:
the approved plan `fgs_chunked_hybrid_plan_20260918.md` (this directory) is
IMPLEMENTED, locally GATED, and DEPLOYED. **Job 13749231 is running on
m12-3-31** (submitted 2026-09-17 ~23:19 UTC, 12 h limit, started 23:21).
State at handoff: all cluster controls PASS (incl. FMM chain with
`fgs_chunked_test.jl`), j64 chunked calibration CERTIFIED (tolerance
5.348662427942506e-7, 44 iterations vs lex 27, confirmation repeat delta
exactly 0.0 — matches the local gate crossing bit-for-near-bit; ARM↔zen3
libm explains the last digits), j1 trials arm in progress. Expect ~10.5 h
total (v21 precedent) → done ~2026-09-18 morning UTC.

**Your task: execute plan §9 + §5 (watch → harvest → analysis → conditional
coloring revert → wrap).** Read the plan §§5–12 IN FULL first; §0–§4 are
background (implementation is finished — do not re-implement or re-gate).

## What was built (all committed; nothing dirty except pre-existing docs)

- FastMultipole `flowpanel-20260817` @ **`c18e4b46`** (= tag
  `campaign/p021-r4-chunked-fm-20260918-v22`): `sweep_order=:chunked`,
  kwarg `chunks::Int=64`; `build_chunk_map` (cost `m_i*n_i + n_i^2`,
  contiguous prefix-sum cuts, j-invariant); per-leaf intra/cross scatter
  partition — a segment is intra ONLY if all rows it writes lie in the
  source's own chunk (multi-leaf target branches straddling chunks are
  always deferred); `gs_sweep!` `:chunked` branch =
  `Threads.@threads :static` over chunks (solve_leaf! →
  compute_nonself_products! → scatter intra) + serial ascending deferred
  cross scatter after the barrier. Unsplit scatter path untouched
  (lexicographic bit-identity preserved by construction).
  `test/fgs_chunked_test.jl` (1723 cases: map validity/purity, nchunks=1 ≡
  lex bitwise, historical-loop lex regression guard, threaded ≡ serial
  chunk-major reference bitwise, cross-process -t1/-t4 bitwise via
  `test/fgs_chunked_threadcheck.jl`, repeat delta 0, finite histories).
- FLOWPanel `fastmultipole` @ **`67d570f`** (= tag
  `campaign/p021-r4-chunked-source-20260918-v22`): FGSSolver/
  FGSPreconditioner `chunks` pass-through + both metadata emitters;
  `fgs_cold_common.jl` — sweep_order ∈ (lexicographic,colored,chunked),
  `chunks` key legal only for chunked (cold_make defaults 64), screen axis
  now proposes chunked; v22 trio `benchmark/fgs_r4_chunked_ab.jl`,
  `benchmark/run_r4_chunked_ab.slurm.sh`,
  `test/runtests_r4_chunked_ab_driver.jl` (near-clones of v21;
  `CHUNKED_CONFIG` env; ab_check_pair allows key-set delta exactly
  {chunks}); solver suite gained a :chunked plumbing case.

Local §7 gate: ALL PASS (fresh 1.11.8 campaign-style env WITH
Meshes+StaticArrays; driver unit test also under 1.12.4; every
launcher-invoked script; driver end-to-end all three AB_MODEs at R4
j4/BLAS1 — local calibrate 5.348662428097527e-7 / 44 iters / delta 0.0;
local j4 medians lex 15.82 s, chunked 19.26 s — gate evidence only, NOT
rankings). Deployment per §8: worktrees + env + `pins.toml` at
`/home/rander39/campaigns/p021-r4-chunked-20260918-v22/` (FLOWVPM pin
unchanged `05c658f7` at `/home/rander39/campaigns/p021-cold-20260910-v1/`),
tags pushed to CLUSTER CLONES ONLY. Full record:
`fgs_r4_followup_evidence_20260914/v22-deployment/submission-provenance-13749231.md`.

## Remaining work (in order)

1. **Watch 13749231** via `hpc-monitor` / a Monitor loop. Judge by OUTPUTS
   (`COMPLETED` sentinel + per-stage `status.toml`), never sacct exit
   status; ≥300 s sacct spacing. Run dir:
   `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/chunked-v22-13749231/`
   (stages: calibrate/, j{1,4,16,64}-trials/, j64-activity/, COMPLETED).
   GOTCHAS (both bit this session): a non-login ssh shell has NO slurm
   binaries on PATH — use `/apps/slurm/latest/bin/squeue` absolute (empty
   output otherwise mimics "job gone"); a login shell (`bash -lc`) injects
   an ANSI banner into stdout. `ssh orc` needs the live ControlMaster
   socket — auth failure = STOP, never MFA. On failure: harvest evidence
   FIRST to `fgs_r4_followup_evidence_20260914/chunked-v22-13749231-FAILED/`,
   reproduce locally, postmortem (v17–v21 pattern; wall-limit contingency
   plan §11).
2. **Harvest** (delegate to `harvester`): full run dir →
   `fgs_r4_followup_evidence_20260914/chunked-v22-13749231/`,
   SHA256-verified against remote (v21 pattern: remote.sha256 manifest).
3. **Analysis** `chunked-v22-13749231/analysis/ab_summary.md` mirroring
   `colored-v21-13738665/analysis/ab_summary.md` (read it as the template):
   gates table (all §6 gates incl. repeat delta exactly 0, arm-invariant
   iterations per order, j1 direct-vs-FMM warmup crosscheck), per-arm
   lex/chunked medians + IQRs, comparison vs BOTH yardsticks (lex@j64
   10.96–11.20 s, colored@j16 **10.116 s**, cited from v21), iteration
   counts (expect chunked 44 / lex 27; ranking metric = total time to
   accepted accuracy absorbs it), j64 activity attribution (did the
   nearfield_update span shrink, at what avg active threads — colored
   engaged 42.65 but did NOT shrink; report vs the ~3.5–4.5 s Amdahl-style
   expectation).
4. **§5 coloring revert — verdict REVERT, executed LAST and ONLY IF**
   chunked's best certified operating point beats colored@j16 10.116 s.
   If it loses: KEEP coloring and report — that is a finding for Ryan, not
   a failure. Revert contents (plan §5, new commits only, v21 tags stay):
   FM — remove `:colored` branch, `color_leaves` + helpers, struct fields
   `leaf_colors`/`leaves_by_color` (keep `sweep_order`), `:colored` from
   validation, `fgs_coloring_test.jl` + registration (verify `417489d5` is
   coloring-unrelated via `git show --stat` first); FLOWPanel — drop
   "colored" from `fgs_cold_common.jl` validation, remove the v21 trio,
   update the solver comment. Full suites of BOTH repos green after. Also
   update the v22 launcher's FMM control chain if the coloring test is
   removed (it currently includes `fgs_coloring_test.jl`) — new commit,
   NOT a tag move; the pinned v22 tree stays as launched.
5. **Wrap:** update `../021_rotor_hover_solver_benchmarks.md`
   `## Current status` + decision log and `log.md` (NOT yet updated for
   v22 — implementation/submission still needs its decision-log entry);
   update memory `project_021_solver_benchmarks.md` (+ MEMORY.md hook —
   already current through submission); OFFER the notebook entry (now 3
   owed: v21 A/B, diagnostics ladder, v22) — never write without Ryan's
   approval.

## Facts you'd otherwise re-derive

- Chunked selected config (cluster-certified): retained R4 lex config
  (`benchmark/retained_r4_diagnostics.toml`, P=8 MAC=0.4 leaf=100 inner=3
  rlx=1.0 tol=3.479128881193055e-7) + sweep_order=chunked, chunks=64,
  tolerance=5.348662427942506e-7. On-cluster copy:
  `<run>/calibrate/results/chunked_selected.toml`.
- Trials: 8 alternating batches × 10 reps per arm (odd=lex, even=chunked),
  j∈{1,4,16,64}; per-trial repeat gate ≤1e-8; cross-order rel-L2
  informational (~3e-6 locally).
- Storage: pre-submit hpc-storage cycle archived 15 runs, /home 725.9→
  **359.1 G (UNDER the 400 G cap)**; ledger line appended to the 018
  ledger. 3 RECENT runs (208 GiB) still await Ryan's
  `--include-recent --only` approval (p018_csarc_l3p0_3r_g25_s2,
  p018_csarc_n2_nt72_l3p0_3r_srlx_g25_s2,
  p018_csarc_n2_nt72_l3p0_3r_srlx_nv_g25).
- Uncommitted-but-owed FLOWPanel docs (commit at wrap together with the
  021 item/log updates): the v22 provenance file, this handoff, the 018
  ledger line. Pre-existing dirty BRAINSTORM/docs files from other
  campaigns — leave them alone.
- Local scratch (gate env/runs) is session-temporary; all gate evidence
  that matters is recorded in the provenance file.

## Ryan-pending (do not act)

- Origin pushes: merged branches + v21 AND v22 campaign tags.
- Notebook entries (3 owed — offer at wrap).
- RECENT VTK archiving approval (208 GiB across 3 runs).
- Any chunked follow-up experiments if the staircase/result disappoints
  (§11: more/fewer chunks, under-relaxation = NEW experiments, Ryan's call).
