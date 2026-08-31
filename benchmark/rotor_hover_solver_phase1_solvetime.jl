#=##############################################################################
BRAINSTORM 021 Phase 1 — clean per-solve timings for the Krylov roster
configs at the frozen/PROPOSED settings (companion to fgs_check.jl, which
already timed fgs/fgmres_fgs the same way; Backslash timings live in
runs.csv). Same frozen-b contract: cold solve via pnl._solve! with the
cached direct source-potential assembly restored (rhs = -potential), so
t_solve excludes the N^2 source assembly and monitor overhead — comparable
across configs. Setup (preconditioner/solver ctor) timed separately.

NOTE: local 4-thread numbers, indicative only — published timings are HPC
(ruling 5). Krylov applies rebuild FMM trees per iteration (Phase 2b lever),
included as production behavior.

Run: RUNG=R1 CONFIGS=krylov_gmres,krylov_jacobi,krylov_ilu \
     EXPECT_JULIA_THREADS=4 THREADING_MODE=multi julia --project -t 4 \
     benchmark/rotor_hover_solver_phase1_solvetime.jl
Appends solvetime.csv under benchmark/results/phase1/<mode>/.
=###############################################################################

import FLOWPanel as pnl
import FastMultipole
include(joinpath(pnl.examples_path, "dji9443_trailing_edge.jl"))
include(joinpath(@__DIR__, "common.jl"))

using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm, cross

banner = assert_and_banner()
results_dir = joinpath(@__DIR__, "results", "phase1", banner.threading_mode)
mkpath(results_dir)

const LADDER = Dict(
    "R1" => ("dji9443_20260813_23_73_capped_captess4.msh",   8016),
    "R2" => ("dji9443_20260813_33_105_capped_captess4.msh", 15760),
    "R3" => ("dji9443_20260813_45_145_capped_captess4.msh", 28752),
)
const TUNED = Dict(
    "R1" => (17, 0.5, 21),
    "R2" => (17, 0.5, 12),
    "R3" => (18, 0.5, 18),
)

rung = get(ENV, "RUNG", "")
haskey(LADDER, rung) || error("RUNG must be one of $(sort(collect(keys(LADDER)))); got \"$rung\"")
msh_name, n_expected = LADDER[rung]
k_reps = parse(Int, get(ENV, "K_REPS", "1"))

# ---- Dirichlet case, mirroring the Phase 1 driver (keep in sync) ----
magVinf = 0.0001; rho = 1.179; RPM = 6000; R = 0.119
Vinf = magVinf * [1.0, 0.0, 0.0]
dt = 60 / RPM / 36
core_size_panel = R * 1e-10
core_size_targets = 1e-3
kernelcutoff = R * 1e-13
radial_dimension = 2
shedding_r_over_R = 0.1

msh_file = joinpath(pnl.examples_path, "data", msh_name)
te_indices_1, te_indices_2 = find_dji9443_trailing_edge_indices(msh_file; watertight=true)
msh = pnl.read_gmsh(msh_file)
nodes0, cells0 = pnl.meshes2nodes_cells(msh)
nodes0 .*= R / maximum(nodes0[radial_dimension, :])

kernel = Union{pnl.ConstantSource, pnl.VortexRing}
body_kwargs = (; core_size=core_size_panel, core_size_panel,
                 core_size_targets, kernelcutoff,
                 semiinfinite_wake=false, watertight=true, DBC=true)
rotor0 = pnl.RigidWakeBody{kernel}(nodes0, cells0, pnl.noshedding; body_kwargs...)

function make_shedding_bbox(nodes, seed_nodes, radial_dimension, R, cutoff_r_over_R)
    radial_midpoint = sum(nodes[radial_dimension, seed_nodes]) / length(seed_nodes)
    radial_sign = sign(radial_midpoint)
    lower = [minimum(nodes[i, :]) for i in 1:size(nodes, 1)]
    upper = [maximum(nodes[i, :]) for i in 1:size(nodes, 1)]
    padding = max(sqrt(eps(eltype(nodes))) * R, R * 1e-6)
    lower .-= padding; upper .+= padding
    if radial_sign > 0
        lower[radial_dimension] = cutoff_r_over_R * R - padding
    else
        upper[radial_dimension] = -cutoff_r_over_R * R + padding
    end
    return (pnl.SVector{3}(lower...), pnl.SVector{3}(upper...))
end
function clip_shedding_root(nodes, shedding, cells, radial_dimension, R, clip_r_over_R)
    keep = Int[]
    for j in axes(shedding, 2)
        p, nia, nib = shedding[1, j], shedding[2, j], shedding[3, j]
        na, nb = cells[nia, p], cells[nib, p]
        mid = (nodes[radial_dimension, na] + nodes[radial_dimension, nb]) / 2
        abs(mid) / R >= clip_r_over_R && push!(keep, j)
    end
    return shedding[:, keep]
end
sheddings = map((te_indices_1, te_indices_2)) do te
    bbox = make_shedding_bbox(rotor0.nodes, te[1:2], radial_dimension, R, 0.0)
    full = pnl.calc_shedding_from_seed(rotor0.nodes, rotor0.cells, te[1], te[2];
        bbox, end_node=te[3], normal_jump_tol=0.2, max_turn_angle=pi/3)
    clip_shedding_root(rotor0.nodes, full, rotor0.cells, radial_dimension, R, shedding_r_over_R)
end
rotor = pnl.RigidWakeBody{kernel}(rotor0.nodes, rotor0.cells, collect(sheddings);
            ensure_winding=true, body_kwargs...)
rotor.ncells == n_expected ||
    error("Rung $rung: expected $n_expected panels, got $(rotor.ncells)")

frames = pnl.ReferenceFrame(rotor;
    origin=SVector{3}(0.0,0.0,0.0), v=SVector{3}(0.0,0.0,0.0),
    ω_axis=SVector{3}(-1.0,0.0,0.0), ω=2*pi*RPM/60,
    R=SMatrix{3,3}(1.0,0,0, 0,1.0,0, 0,0,1.0),
    name="vehicle", child_index=Int[], dependent_index=[1])
pnl.initialize_Das!((rotor,), frames, t -> Vinf, 0.0, dt; set_Das_eta_kinematic=0.2)
ω_vec = SVector{3}(-2*pi*RPM/60, 0.0, 0.0)
pnl.calc_normals!(rotor); pnl.calc_controlpoints!(rotor)
for i in 1:rotor.ncells
    rr = SVector{3}(rotor.controlpoints[:, i]...)
    rotor.velocity[:, i] .= Vinf .- cross(ω_vec, rr)
end
frozen_velocity = copy(rotor.velocity)

pnl.set_strengths!(rotor)
rotor.potential .= 0
pnl.influence!(rotor, rotor, pnl.DirectBackend(); scalar_potential=true, velocity=false)
potential_frozen = copy(rotor.potential)
b = zeros(rotor.ncells)
pnl.assemble_rhs!(b, rotor)
rms_b = norm(b) / sqrt(length(b))
println("$rung Dirichlet case: $(rotor.ncells) panels, rms(b) = $rms_b\n")

solution_column = 2
function reset_cold!()
    pnl.set_strengths!(rotor)
    rotor.strength[:, solution_column] .= 0
    rotor.velocity .= frozen_velocity
    rotor.potential .= potential_frozen
    return nothing
end

p_t, mac_t, leaf_t = TUNED[rung]
backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)
krylov_kw = (; itmax=500, atol=1e-14, rtol=1e-6, memory=50, backend=backend_apply)

rotor.core_size = rotor.core_size_panel   # ctor-time trees at panel offset
configs = [
    # dense reference per-solve cost: ctor = G assembly + LU (setup column);
    # Dirichlet _solve! with update_G=false is rhs-copy + triangular ldiv only
    ("backslash_ldiv", () -> pnl.Backslash(rotor)),
    ("krylov_gmres",  () -> pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...)),
    ("krylov_jacobi", () -> pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...,
                          preconditioner=pnl.FastMultipole.JacobiPreconditioner((rotor,); cell_size=R/4))),
    ("krylov_ilu",    () -> pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...,
                          preconditioner=pnl.ILUPreconditioner(rotor; leaf_size=10, multipole_acceptance=1.0,
                              max_pattern_entries=8192 * rotor.ncells))),
]
if haskey(ENV, "CONFIGS")
    want = split(ENV["CONFIGS"], ",")
    configs = [c for c in configs if c[1] in want]
end

header = "rung,mesh_file,n_panels,config,t_setup_ctor,t_solve_min,k_reps,niter," *
    "rms_residual_direct,rel_rms_direct,rms_b,radius_tol,threading_mode," *
    "julia_threads,commit,fm_commit,date,notes"
csv_path = joinpath(results_dir, "solvetime.csv")
fresh = !isfile(csv_path) || filesize(csv_path) == 0
io = open(csv_path, "a")
fresh && println(io, header)

r = zeros(rotor.ncells)
for (name, make_solver) in configs
    println("--- $name ---")
    t_setup = @elapsed solver = make_solver()
    t_solve, _ = min_of_k(() -> pnl._solve!(rotor, solver);
                          k=k_reps, warmup=1, setup! = reset_cold!)
    x = vec(rotor.strength[:, solution_column])
    rms_d, _ = pnl.true_residual!(r, rotor, x, b; backend=pnl.DirectBackend())
    niter = solver isa pnl.KrylovSolver || solver isa pnl.FGSSolver ? solver.niter : -1
    println("  setup=$(round(t_setup; digits=2))s solve=$(round(t_solve; digits=4))s " *
            "niter=$niter rel_rms=$(rms_d/rms_b)")
    println(io, join(_csv_cell.([rung, msh_name, rotor.ncells, name, t_setup,
        t_solve, k_reps, niter, rms_d, rms_d / rms_b, rms_b,
        pnl.FMM_RADIUS_TOL[], banner.threading_mode, banner.julia_threads,
        banner.commit, banner.fm_commit, time_string(),
        "frozen settings; t excludes N^2 source assembly (restored); " *
        "apply backend tuned($p_t,$mac_t,$leaf_t); tree rebuild per apply included"]), ","))
    flush(io)
    solver = nothing; GC.gc()
end
close(io)
println("\n$rung solvetime done. Rows appended to $csv_path")
