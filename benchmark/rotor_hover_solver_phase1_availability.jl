#=##############################################################################
BRAINSTORM 021 Phase 1 (B5) — roster availability re-proof on the DIRICHLET case.

Phase 0 proved the seven roster configs on the Neumann (uncapped) case only; the
campaign case is now Dirichlet per Ryan's 2026-08-13 ruling (RHPC formulation).
This script builds the R1 Dirichlet case exactly as the Phase 1 driver does and
attempts each roster config once at smoke-scale settings (rtol 1e-4, itmax 400):
construct → one cold solve → direct true residual vs the frozen Dirichlet b.

Output: PASS/FAIL table on stdout + availability.csv (append) in
benchmark/results/phase1/<mode>/. A FAIL is a finding for Ryan (config cannot
operate on the DBC system), not something to hack around.

Run:
  EXPECT_JULIA_THREADS=4 THREADING_MODE=multi julia --project -t 4 \
      benchmark/rotor_hover_solver_phase1_availability.jl
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

# ---- R1 Dirichlet case, mirroring the Phase 1 driver's construction ----
# (the driver is a script-style file, so its case section cannot be included;
# keep this block in sync with it)
magVinf = 0.0001; rho = 1.179; RPM = 6000; R = 0.119
Vinf = magVinf * [1.0, 0.0, 0.0]
dt = 60 / RPM / 36
core_size_panel = R * 1e-10
core_size_targets = 1e-3
kernelcutoff = R * 1e-13
radial_dimension = 2
shedding_r_over_R = 0.1

msh_name = "dji9443_20260813_23_73_capped_captess4.msh"
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

# frozen Dirichlet b (direct source-potential assembly)
pnl.set_strengths!(rotor)
rotor.potential .= 0
pnl.influence!(rotor, rotor, pnl.DirectBackend(); scalar_potential=true, velocity=false)
b = zeros(rotor.ncells)
pnl.assemble_rhs!(b, rotor)
rms_b = norm(b) / sqrt(length(b))
println("R1 Dirichlet case: $(rotor.ncells) panels, rms(b) = $rms_b\n")

# ---- roster attempts ----
solution_column = 2
backend_sim = pnl.FastMultipoleBackend()
krylov_kw = (; method=:gmres, itmax=400, atol=1e-6, rtol=1e-4, backend=backend_sim)
# NOTE the potential zero: coupled solvers (BackslashCoupled/KrylovCoupled)
# treat entry body.potential as an EXTERNAL potential and add it to their
# internally assembled source potential. Leaving the b-capture's assembly in
# rotor.potential doubled the coupled Dirichlet RHS (x = 2*x_ref, rel_rms
# exactly 1.0 — the 2026-08-13/14 "backslash_coupled x~0" rows measured this
# script bug, not a solver defect). Single-body paths zero it themselves.
reset!() = (rotor.strength[:, solution_column] .= 0; rotor.velocity .= frozen_velocity;
            rotor.potential .= 0; nothing)

configs = [
    ("backslash_coupled", () -> begin
        solver = pnl.BackslashCoupled((rotor,))
        pnl.solve!((rotor,), solver)
        solver
    end),
    ("backslash_iterative", () -> begin
        solver = pnl.Backslash(rotor)
        pnl.solve!((rotor,), (solver,))
        solver
    end),
    ("krylov_gmres", () -> begin
        solver = pnl.KrylovSolver(rotor; krylov_kw...)
        pnl.solve!(rotor, solver)
        solver
    end),
    ("krylov_jacobi", () -> begin
        P = pnl.FastMultipole.JacobiPreconditioner((rotor,); cell_size=R/4)
        solver = pnl.KrylovSolver(rotor; krylov_kw..., preconditioner=P)
        pnl.solve!(rotor, solver)
        solver
    end),
    ("fgs", () -> begin
        solver = pnl.FGSSolver(rotor; max_iterations=500, tolerance=1e-6, rlx=1.0,
            expansion_order=10, multipole_acceptance=0.4, leaf_size=150,
            shrink=true, recenter=false, inner_iterations=20,
            reverse_pass=false, verbose=false)
        pnl.solve!(rotor, solver)
        solver
    end),
    ("fgmres_fgs", () -> begin
        P = pnl.FGSPreconditioner(rotor; sweeps=1, inner_iterations=2,
            expansion_order=10, multipole_acceptance=0.4, leaf_size=150)
        solver = pnl.KrylovSolver(rotor; krylov_kw..., method=:fgmres, preconditioner=P)
        pnl.solve!(rotor, solver)
        solver
    end),
    ("krylov_ilu", () -> begin
        P = pnl.ILUPreconditioner(rotor; leaf_size=10, multipole_acceptance=1.0)
        solver = pnl.KrylovSolver(rotor; krylov_kw..., preconditioner=P)
        pnl.solve!(rotor, solver)
        solver
    end),
]

io = open(joinpath(results_dir, "availability.csv"), "a")
filesize(joinpath(results_dir, "availability.csv")) == 0 &&
    println(io, "config,status,rms_residual,rel_rms,error,commit,date")
r = zeros(rotor.ncells)
for (name, attempt) in configs
    reset!()
    status, rms_v, err_msg = try
        attempt()
        x = vec(rotor.strength[:, solution_column])
        rms_v, _ = pnl.true_residual!(r, rotor, x, b; backend=pnl.DirectBackend())
        ("PASS", rms_v, "")
    catch err
        ("FAIL", NaN, sprint(showerror, err)[1:min(end, 200)])
    end
    rel = rms_v / rms_b
    println("$(rpad(name, 20)) $status  rel_rms=$(rel)  $(err_msg)")
    println(io, join(_csv_cell.([name, status, rms_v, rel,
        replace(err_msg, "\n" => " "), banner.commit, time_string()]), ","))
end
close(io)
println("\nAvailability results appended to $(joinpath(results_dir, "availability.csv"))")
