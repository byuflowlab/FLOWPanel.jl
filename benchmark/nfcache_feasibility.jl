#=##############################################################################
BRAINSTORM 021 W0 — near-field influence-matrix cache FEASIBILITY MAP.

Sweeps FastMultipole.estimate_nearfield_cache over a leaf ladder at one rung,
WITHOUT building any cache and WITHOUT solving. Answers: how large is the
cache, and how long would it take to build, as a function of leaf_size_source?

Why this exists (phase_15_caching_and_objective_prompt.md §1, Ryan 2026-08-24):
an uncapped cached tune walks to a DENSE operator — measured at R1, the cached
optimum landed at leaf 275-342 with 687-862 MB against a dense 8N^2 of 514 MB
(1.3-1.7x dense), with cache_capped=false, i.e. the cost optimum genuinely
wanted it. If assembly is free, matrix-ful always beats matrix-free. So
max_bytes is not a safety valve — it silently DETERMINES the published knob.
Ryan's ruling: the memory budget becomes a swept AXIS in Phase 2, and this map
is what picks the budget ladder.

Key pruning fact this map establishes per rung: budgets at or above the rung's
dense 8N^2 all admit the same (dense) optimum and collapse to a single
"unconstrained" point, so only budgets BELOW dense are distinct.

`bytes` is exact size-pass arithmetic — deterministic, hardware-independent,
and therefore VALID LOCALLY. `est_build_time` scales one sampled kernel
evaluation and IS a timing, so it is indicative only off-cluster (Ryan's Mac
scatters 22-39% at fixed knobs); SAMPLE=0 reports NaN for it and makes the
whole run deterministic.

Run:
  RUNG=R1 SKIP_B=1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    julia --project -t 4 benchmark/nfcache_feasibility.jl

Env:
  RUNG      required, R1..R7
  LEAVES    colon/comma-separated leaf ladder (default 9:15:21:30:45:70:100:
            150:225:340:500:750:1100)
  P, MAC    expansion order / multipole acceptance held fixed across the sweep
            (default: the rung's TUNED triple)
  SAMPLE    1 (default) to estimate build time, 0 for bytes only

SKIP_B=1 is expected: this script never solves, so the O(N^2) frozen-b
assembly (tens of minutes at R7) is pure waste.
Appends nfcache_feasibility.csv under benchmark/results/phase2/<mode>[/<rung>].
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

import FLOWPanel.FastMultipole as FM

# TUNED only carries rungs that have been tuned; this map does not need tuned
# knobs. `bytes` depends on the direct list, i.e. on MAC (and leaf), and is
# INDEPENDENT of expansion_order, so a fixed reference MAC across all rungs is
# what makes the map comparable rung-to-rung anyway.
const REF_P, REF_MAC = 17, 0.5
p_t, mac_t = haskey(TUNED, rung) ? TUNED[rung][1:2] : (REF_P, REF_MAC)
haskey(TUNED, rung) || println("no TUNED[$rung]; using reference " *
    "p=$REF_P mac=$REF_MAC (bytes are independent of p)")
p_sweep = parse(Int, get(ENV, "P", string(p_t)))
mac_sweep = parse(Float64, get(ENV, "MAC", string(mac_t)))
sample = get(ENV, "SAMPLE", "1") == "1"

const DEFAULT_LEAVES = "9:15:21:30:45:70:100:150:225:340:500:750:1100"
leaves = sort(unique(parse(Int, s) for s in
                     split(get(ENV, "LEAVES", DEFAULT_LEAVES), r"[:,]")))

outdir = joinpath(@__DIR__, "results", "phase2", banner.threading_mode)
get(ENV, "PER_RUNG_DIR", "0") == "1" && (outdir = joinpath(outdir, rung))
mkpath(outdir)
csv_path = joinpath(outdir, "nfcache_feasibility.csv")
header = "rung,mesh_file,n_panels,expansion_order,multipole_acceptance," *
    "leaf_size,cache_bytes,dense_bytes,bytes_over_dense,est_build_time," *
    "n_blocks,total_probe_pairs,t_plan_build,sampled,threading_mode," *
    "julia_threads,blas_threads,commit,fm_commit,date,notes"
fresh = !isfile(csv_path) || filesize(csv_path) == 0
if !fresh
    first(eachline(csv_path)) == header || error(
        "$csv_path has a pre-W0 schema; move it aside before appending")
end
io = open(csv_path, "a")
fresh && println(io, header)

n = rotor.ncells
dense_bytes = 8 * n^2   # the backslash_ldiv endpoint: a full dense G
println("\n--- nfcache feasibility: $rung, $n panels, p=$p_sweep mac=$mac_sweep ---")
println("dense 8N^2 = $(round(dense_bytes / 1024^3; digits=3)) GiB\n")

for leaf in leaves
    # The plan mirrors the Dirichlet operator FLOWPanel builds in
    # _apply_dirichlet_G! (scalar potential only, shrink=true) — the cache
    # blocks are sized by the output rows that operator requests, so a
    # different derivative switch would give different bytes.
    t_plan = @elapsed plan = FM.FmmPlan((rotor,), (rotor,);
        scalar_potential=true, gradient=false, hessian=false,
        expansion_order=p_sweep, multipole_acceptance=mac_sweep,
        leaf_size_source=leaf, shrink=true)
    est = FM.estimate_nearfield_cache(plan.target_tree, plan.source_tree,
        plan.direct_list, plan.derivatives_switches, (rotor,); sample)
    ratio = est.bytes / dense_bytes
    println("  leaf $(lpad(leaf, 5)): $(lpad(round(est.bytes / 1024^3; digits=3), 8)) GiB " *
            "($(lpad(round(ratio; digits=3), 7))x dense), $(est.n_blocks) blocks, " *
            "est build $(round(est.est_build_time; digits=2)) s")
    println(io, join(_csv_cell.([rung, msh_name, n, p_sweep, mac_sweep, leaf,
        est.bytes, dense_bytes, ratio, est.est_build_time, est.n_blocks,
        est.total_probe_pairs, t_plan, sample, banner.threading_mode,
        banner.julia_threads, banner.blas_threads, banner.commit,
        banner.fm_commit, time_string(),
        "estimate only, nothing built, nothing solved"]), ","))
    flush(io)
    plan = nothing; GC.gc()
end
close(io)
println("\nwrote $csv_path")
