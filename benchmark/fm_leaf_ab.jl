#=##############################################################################
021 diagnostic (2026-08-24): is the new FastMultipole slower per apply, or did
stock tune_fmm simply choose a worse leaf_size?

R1 multi under FM d8258a7d came out ~2.5x slower per FMM apply than under
a9b734a at IDENTICAL iteration counts and BC error, with p/mac unchanged (17/0.5)
but apply leaf 9 -> 60. Two mutually exclusive explanations:
  (1) leaf 60 is a bad choice and the new FM is fast at leaf 9;
  (2) the new FM is genuinely ~2.5x slower per apply and the tuner enlarged leaf
      to compensate.
This sweeps leaf at fixed p/mac on the CURRENT FM and reports cost. If leaf 9 is
fast, it is (1). If every leaf is slow, it is (2) — a FastMultipole regression.
Writes nothing to the campaign results tree.
=###############################################################################
import FLOWPanel as pnl
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

const P, MAC = 17, 0.5
const LEAVES = (9, 21, 40, 60, 71)
target = 1e-6

println("\n=== leaf sweep at p=$P mac=$MAC — rung $rung, $(rotor.ncells) panels ===")
println("fm_commit = $(banner.fm_commit)")
println(rpad("leaf", 8), rpad("t_solve_min [s]", 18), rpad("niter", 8), "bc_rel_l2")
for leaf in LEAVES
    backend = pnl.FastMultipoleBackend(; expansion_order=P,
        multipole_acceptance=MAC, leaf_size=leaf)
    solver = pnl.KrylovSolver(rotor; method=:gmres, itmax=500, atol=1e-14,
        rtol=target, memory=50, backend=backend)
    # Use phase1_case's reset_cold! rather than a local partial reset. The
    # Dirichlet rhs is `-rotor.potential` (assemble_rhs! contract) and a solve
    # CLOBBERS potential, so a reset that restores only strengths and velocity
    # leaves every solve after the first one solving against a stale b. That is
    # why this script's bc_rel_l2 column previously read ~1.0-1.1 (2026-08-24);
    # the timings were unaffected (same operator, same niter to within 1) but
    # the accuracy column was meaningless.
    reset_cold!(); pnl._solve!(rotor, solver)              # warmup / JIT
    t = minimum(begin reset_cold!(); @elapsed pnl._solve!(rotor, solver) end for _ in 1:3)
    x = copy(vec(rotor.strength[:, solution_column]))
    e = bc_error!(rotor, x; rms_b, target_rel=target, safety=1.0,
                  max_expansion_order=20, multipole_acceptance=MAC,
                  leaf_size=leaf, backend=:fmm)
    println(rpad(leaf, 8), rpad(round(t; digits=4), 18), rpad(solver.niter, 8), e.rel_l2)
    solver = nothing; GC.gc()
end
