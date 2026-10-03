#=##############################################################################
034 Phase 0 deliverable 3: MEASURE whether the LHS (Backslash dense G + LU) is
rebuilt/refactored per step under the pitching maneuver.

Code hypothesis (to be confirmed here, file:line evidence in phase_00_harness.md):
- `Backslash(body)` assembles G and factorizes it ONCE (lu! aliases G's memory,
  src/FLOWPanel_solver.jl:489-497); per-step `_solve!` only refreshes the RHS
  unless `update_G=true` (src/FLOWPanel_solver.jl:1714-1736), and nothing in the
  simulate! path ever passes update_G=true.
- Per step, simulate! rotates the body nodes + Das (propagate_kinematics!,
  FLOWPanel_simulate.jl:1495) and mirrors the rigid delta into persistent solver
  state via transform_body_solvers! (:1505) — a documented no-op for Backslash
  because the dense Dirichlet operator is rotation-invariant
  (src/FLOWPanel_solver.jl:1248-1249).

Measurements on ladder rung 1 (n_span=7, n_airfoil=89, n_endcap=5 → 1920 cells),
~20 steps of the pitching maneuver, panel wake, default FMM backend:
1. per-step max|G - G_at_construction| and objectid(solver.Glu) — zero deviation
   and a stable objectid prove no rebuild/refactor occurs;
2. rotation-invariance: raw G re-assembled from the END-of-run (rotated)
   geometry vs the raw G assembled at t=0 — max-abs and relative-Frobenius diff;
3. counterfactual per-step rebuild cost on this rung: min-of-k G assembly time
   and LU time, plus the constructor's own recorded g_assembly_s / lu_s.

Run (local smoke; BENCH_BLAS_THREADS=8 per the macOS OpenBLAS note below):
      THREADING_MODE=single EXPECT_JULIA_THREADS=1 BENCH_BLAS_THREADS=8 \
          julia --project -t 1 benchmark/p034_phase0_lhs_rebuild.jl
Judge from the CSVs in BRAINSTORM/034_pitching_wing_solver_benchmarks/data/,
never stdout.
=###############################################################################

import FLOWPanel as pnl
import LinearAlgebra
import LinearAlgebra: lu!, norm
using Printf

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "..", "examples", "pitching_wing.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "BRAINSTORM",
    "034_pitching_wing_solver_benchmarks", "data", "phase0_lhs_rebuild")
mkpath(OUTDIR)

# Local-smoke quirk (macOS OpenBLAS 0.3.31, 2026-10-02): after any GEMM,
# BLAS.get_num_threads() reports the hardware count (8) regardless of what was
# set, and a 2000x2000 GEMM times identically at 1 vs 4 threads — the pin is
# neither observable nor effective on this build. assert_and_banner's strict
# single-mode pin is therefore unsatisfiable here; launch local smokes with
# BENCH_BLAS_THREADS=8 so the assert passes and the banner records the truth.
# Published (HPC/Linux) runs keep the strict pin. Smoke only — never publish.
banner = assert_and_banner()
open(joinpath(OUTDIR, "banner.txt"), "w") do io
    println(io, banner.text)
end

# ladder rung 1 (provisional, Phase 0)
sim = prepare_pitching_wing(;
    n_span=7, n_airfoil=89, n_endcap=5,
    n_cycles=0.12,                      # ~20 steps at c_per_dt = 0.5
    include_static_polar=false,
    save_vtk=false,
)
wing, solver = sim.wing, sim.solver
n = wing.ncells
n == 1920 || error("rung-1 cell count drifted: expected 1920, got $n")
nsteps = length(sim.t_range)

# raw G at t=0 geometry (the constructor's G is destroyed in place by lu!)
G0_raw = zeros(Float64, n, n)
pnl._G!(G0_raw, wing, wing; core_size=wing.core_size_panel, update_geometry=false)

# factor storage + factorization object as constructed
Gfact0 = copy(solver.G)
glu_id0 = objectid(solver.Glu)

# instrument via the maneuver callback: called once per step before the solve
step_dev = Float64[]
glu_stable = Bool[]
maneuver_inner! = sim.maneuver!
function instrumented_maneuver!(frames, systems, wakes, t)
    push!(step_dev, maximum(abs, solver.G .- Gfact0))
    push!(glu_stable, objectid(solver.Glu) == glu_id0)
    return maneuver_inner!(frames, systems, wakes, t)
end

t_sim = @elapsed pnl.simulate!((wing,), (sim.wake,), sim.frames,
    instrumented_maneuver!, sim.Uinf, sim.t_range;
    body_solvers=(solver,),
    backend=sim.backend,
    monitors=sim.monitors,
    path=nothing,
    name="p034_phase0_lhs",
    set_Das_eta_freestream=NaN,
    verbose=false,
)

# final check after the last step's kinematics
final_dev = maximum(abs, solver.G .- Gfact0)
final_glu_stable = objectid(solver.Glu) == glu_id0

open(joinpath(OUTDIR, "per_step.csv"), "w") do io
    println(io, "step,g_maxabs_dev_from_construction,glu_objectid_stable")
    for (i, (dev, stable)) in enumerate(zip(step_dev, glu_stable))
        println(io, "$(i-1),$dev,$stable")
    end
    println(io, "post_final,$final_dev,$final_glu_stable")
end

# rotation-invariance: re-assemble raw G from the END-of-run rotated geometry
pnl.calc_normals!(wing)
pnl.calc_controlpoints!(wing)
G1_raw = zeros(Float64, n, n)
t_assembly, _ = min_of_k(k=5) do
    fill!(G1_raw, 0.0)
    pnl._G!(G1_raw, wing, wing; core_size=wing.core_size_panel,
        update_geometry=false)
end
inv_maxabs = maximum(abs, G1_raw .- G0_raw)
inv_relfro = norm(G1_raw .- G0_raw) / norm(G0_raw)

# counterfactual LU cost (factorize a copy so G1_raw survives the k reps)
G_scratch = similar(G1_raw)
t_lu, _ = min_of_k(k=5) do
    G_scratch .= G1_raw
    lu!(G_scratch)
end

ctor = get(pnl._BACKSLASH_CONSTRUCTION_TIMINGS, solver, nothing)

open(joinpath(OUTDIR, "summary.csv"), "w") do io
    println(io, "key,value")
    println(io, "n_panels,$n")
    println(io, "n_steps,$nsteps")
    println(io, "steps_instrumented,$(length(step_dev))")
    println(io, "max_step_dev,$(isempty(step_dev) ? NaN : maximum(step_dev))")
    println(io, "final_dev,$final_dev")
    println(io, "all_glu_stable,$(all(glu_stable) && final_glu_stable)")
    println(io, "rotation_invariance_maxabs,$inv_maxabs")
    println(io, "rotation_invariance_relfro,$inv_relfro")
    println(io, "t_assembly_min_s,$t_assembly")
    println(io, "t_lu_min_s,$t_lu")
    println(io, "ctor_g_assembly_s,$(isnothing(ctor) ? NaN : ctor.g_assembly_s)")
    println(io, "ctor_lu_s,$(isnothing(ctor) ? NaN : ctor.lu_s)")
    println(io, "t_simulate_total_s,$t_sim")
    # denominator = nsteps-1 timesteps marched (t_range has nsteps points; the
    # maneuver callback also fires at t=0, hence steps_instrumented = nsteps)
    println(io, "t_simulate_per_step_s,$(t_sim / max(nsteps - 1, 1))")
end

println("Wrote $(joinpath(OUTDIR, "per_step.csv")) and summary.csv")
