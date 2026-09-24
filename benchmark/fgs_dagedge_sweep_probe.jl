#=##############################################################################
BRAINSTORM 021 L-shortening item #1: local raw-sweep timing sanity probe for
sweep_order=:dagedge vs :dagteam (NOT a benchmark — laptop thread counts
cannot show the j64 DAG-width starvation the edge split targets; this probe
only detects pathological per-task/scheduler overhead before any Ryan-gated
HPC campaign, and cross-checks executor agreement on a real rotor operator).

Times dagteam_initialize! + N x inner_sweeps! directly on a synthetic rhs
(no FMM, no evaluator), so the measured quantity is exactly the nearfield
sweep the gate-0 L analysis bounds.

Usage (local, ≤4 threads per house rules):
  RUNG=R2 SKIP_B=1 THREADING_MODE=multi EXPECT_JULIA_THREADS=4 \
    BENCH_BLAS_THREADS=8 DAGTEAM_PRECISION=f32full NSWEEPS=50 \
    julia --project=benchmark -t 4 benchmark/fgs_dagedge_sweep_probe.jl
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

using Printf
import FastMultipole

precision = Symbol(get(ENV, "DAGTEAM_PRECISION", "f32full"))
nsweeps = parse(Int, get(ENV, "NSWEEPS", "50"))
theta = parse(Int, get(ENV, "DAGEDGE_THETA", "4096"))

make(order) = pnl.FGSSolver(rotor; expansion_order=8, multipole_acceptance=0.4,
    leaf_size=100, inner_iterations=2, max_iterations=100, tolerance=1e-6,
    rlx=1.0, shrink=true, recenter=false, reverse_pass=false,
    cache_leaf_lu=true, sweep_order=order, dagteam_precision=precision,
    dagedge_theta=theta, verbose=false)

println("Constructing solvers (RUNG=$(rung), precision=$precision, theta=$theta)...")
t0 = time_ns()
s_team = make(:dagteam)
@printf("  dagteam constructed in %.1f s\n", (time_ns() - t0) / 1e9)
t0 = time_ns()
s_edge = make(:dagedge)
@printf("  dagedge constructed in %.1f s\n", (time_ns() - t0) / 1e9)

pe = s_edge.fgs.dagteam
base = pe.base
n_leaves = length(base.preds)
n_edge = sum(count(t -> t.kind == 0x01, lst) for lst in pe.lists)
n_small = sum(count(t -> t.kind == 0x02, lst) for lst in pe.lists)
n_back = sum(count(t -> t.kind == 0x04, lst) for lst in pe.lists)
@printf("plan: %d leaves, %d lower edges; dagedge tasks/sweep = %d (big %d, small-agg %d, back %d, finalizes %d)\n",
    n_leaves, sum(length, base.preds), pe.ntasks, n_edge, n_small, n_back, n_leaves)
@printf("sim (bytes): makespan %.1f MB, edge-L %.1f MB, workers %d\n",
    pe.sim_makespan / 1e6, pe.sim_edge_L / 1e6, length(pe.lists))

function run_sweeps!(solver, nsweeps)
    fgs = solver.fgs
    nstr = length(fgs.strengths)
    b = sin.(1.0 .* (1:nstr))
    fgs.strengths .= 0
    rhs = fgs.self_matrices.rhs
    rhs .= b
    zero_ff = zero(fgs.extra_right_hand_side)
    FastMultipole.dagteam_initialize!(rhs, fgs.dagteam, fgs.strengths)
    inner! = fgs.dagteam isa FastMultipole.DagEdgePlan ?
        FastMultipole.dagedge_inner_sweeps! : FastMultipole.dagteam_inner_sweeps!
    inner!(rhs, zero_ff, fgs.strengths, fgs.dagteam, 2)   # compile + warm
    fgs.strengths .= 0
    rhs .= b
    FastMultipole.dagteam_initialize!(rhs, fgs.dagteam, fgs.strengths)
    t0 = time_ns()
    inner!(rhs, zero_ff, fgs.strengths, fgs.dagteam, nsweeps)
    dt = (time_ns() - t0) / 1e9
    return dt, copy(fgs.strengths)
end

dt_team, x_team = run_sweeps!(s_team, nsweeps)
dt_edge, x_edge = run_sweeps!(s_edge, nsweeps)
using LinearAlgebra: norm
dev = norm(x_team .- x_edge, Inf) / max(norm(x_team, Inf), eps())

@printf("\ndagteam: %.4f s total, %.3f ms/sweep\n", dt_team, 1e3 * dt_team / nsweeps)
@printf("dagedge: %.4f s total, %.3f ms/sweep\n", dt_edge, 1e3 * dt_edge / nsweeps)
@printf("ratio dagteam/dagedge = %.3f (j=%d; L-bound only binds at high j)\n",
    dt_team / dt_edge, Threads.nthreads())
@printf("rel dev dagedge vs dagteam after %d sweeps: %.2e\n", nsweeps, dev)
