# FGS acceleration implementation — status 2026-09-18 (gates 1+2 built)

Implements the first moves of `fgs_acceleration_reset_prompt_20260918.md`
(spec: `fgs_acceleration_recommendation_20260918.md`). No production solver
code touched; all work is on a dedicated FastMultipole branch.

## Provenance

- FastMultipole dev worktree: `/private/tmp/fastmultipole-p021-fgs-accel-20260918`,
  branch `p021-fgs-accel-20260918` off `c18e4b46` (the Manifest-pinned HEAD;
  live checkout dirty only in `MATRIX_OPERATOR_REFACTOR/` docs, `src/` clean).
- Commit: `b0946c36` — gate-1 harness + gate-2 replay benchmark + HPC driver.
- All spec code anchors re-verified at `c18e4b46` (`gs_sweep!` now at
  `solve.jl:1214`; all others exact).

## Gate 1 — correctness model (BUILT, PASSING locally)

`test/fgs_rowpar_gate1_test.jl` (standalone or via runtests): a *shadow
executor* of the source-major row-parallel schedule runs against real
`FastGaussSeidel` state — production leaf solves kept; per leaf the nonself
product+scatter is split into contiguous row tiles executed by
`Threads.@threads`, preserving `+= old` / `-= new` as two per-row ops on the
production buffers.

| Check | Result |
|---|---|
| ntiles=1 vs production `:lexicographic` `gs_sweep!` | **bitwise identical** (required) |
| ntiles=2,3,4 at 4 threads, 3 sweeps, nonzero starts | **bitwise identical** (worst rel dev 0.0) |
| rigidly transformed solver (`transform_solver!` fixture) | **bitwise identical** |
| multiple source systems | **fixture cannot construct** — see finding |

Sub-range BLAS `dgemv` on Apple's stack is bit-identical to the full call;
this must be re-checked on the HPC OpenBLAS build (different kernels per m),
else certified via accuracy per the spec's trap list.

**Finding:** a two-system `FastGaussSeidel((sysA,sysB),(sysA,sysB))`
construction fails at `c18e4b46` with
`BoundsError ... reshape(view(::Vector{Float64},1:39824), 524, 76) at index [523:526, 40:76]`
in the nonself matrix fill — pre-existing, unrelated to this branch (no
existing test covers multi-system FGS). The R4 benchmark is single-system, so
gates proceed; the implementation-side gate must revisit (or the limitation
must be declared) before multi-system use.

## Gate 2 — real-shape sequence replay (BUILT; decision numbers await zen3)

`benchmark/fgs_sequence_replay.jl`: reconstructs the exact R4 sequence from
the saved census/edges (both git-tracked in FLOWPanel). Structural
verification against the spec, all exact: 1,068 leaves; per-source
`m_j = Σ n_i` over dependents for every leaf; 2,862,850,032 coefficient
bytes/sweep; 48,167 lower edges; unit longest path 279; byte-weighted
work/span 2.849. Modes:

- `serial` — lex loop, 1-thread BLAS (baseline; reproduces the 29.4 GB/s arm)
- `rowpar` — persistent spin-wait team (GC-safe), coordinator solves cached
  LU + publishes, workers do contiguous row tiles with private scratch and
  disjoint-segment scatter; adaptive small-block serial policy
  (`--small-bytes`); `--first-touch serial|owner`; `--f32` = Float32 storage
  with Float64 accumulate (convert-on-load kernel, strengths stay F64)
- `handoff` — same machinery, zero-work payload → µs per leaf handoff
- `dag` — pull-DAG audit: structural spans + measured-cost recurrence C_i
  (in-process per-leaf `ldiv!` timings) + critical-path list-schedule
  simulation with backward-filler drain model, priced by `--dag-bandwidth`
  / `--dag-handoff-us`

### Local smoke (Apple M2, ≤4 threads — indicative ONLY, not decision-grade)

| Arm | Result |
|---|---|
| handoff, team=4 | **2.52 µs avg** (budget < 5.8 µs) → 86,508 handoffs ≈ 0.22 s/solve |
| serial F64, 1 core | 18.3 GB/s useful → 12.7 s stream/solve |
| serial F32-storage, 1 core | 4.82 s stream/solve (≈2.6× vs serial F64) |
| rowpar F64 t2 (P-cores, owner touch) | 32.4 GB/s → 7.2 s stream |
| rowpar F64 t4 | SLOWER than serial — M2 E-cores drag the balanced tiles; uniform-core HPC nodes unaffected |
| dag @ B=29.4 GB/s, h=2.5 µs | span-limited ≥ W=8: 1.74–2.13 s lower-stream/solve; backward fully absorbed as filler at W≥8 |

"Useful BW" is always reported in F64-equivalent bytes so F32 arms compare
directly; projected 81-sweep stream times are the unambiguous numbers.

The M2 t4 regression is the expected asymmetric-core artifact (t2 on P-cores
beats serial), not a schedule failure; the handoff budget — the spec's kill
switch — passes with 2.3× margin locally.

### HPC driver (prepared, NOT submitted — Ryan-gated)

`benchmark/fgs_sequence_replay_orc.slurm.sh`: zen3 node, exclusive,
`qos=test` (<1 h), socket-0 `numactl` pinning, matrix = serial baselines ×
handoff ladder × rowpar {F64,F32}×{serial,owner touch}×{4,8,16,32} threads ×
interleave control, then `dag` priced with that node's measured B and h.
Census paths default to `~/projects/FLOWPanel.jl/BRAINSTORM/...` (tracked;
confirm the evidence commit is pushed/pulled on orc before submitting).
Open question for Ryan: does this diagnostic merit full campaign ceremony
(tagged worktree pins) or is the dev branch + recorded SHA sufficient for a
go/no-go benchmark whose numbers will be re-established in the end-to-end A/B?

## Next

1. Ryan: approve zen3 gate-2 submission (and campaign-vs-diagnostic call).
2. From its numbers: go/no-go on source-major rowpar (handoff < 5.8 µs and
   projected ≤ 6.744 s with margin), and the source-major vs pull-DAG pick.
3. Only after go: implement the persistent team + mixed-storage container
   change in FastMultipole proper (gate-1 harness then runs against the real
   implementation, plus the multi-system question).
