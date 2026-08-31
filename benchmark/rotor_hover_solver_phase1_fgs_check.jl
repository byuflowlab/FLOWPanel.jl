#=##############################################################################
BRAINSTORM 021 Phase 1 — PENDING CHECK (Ryan 2026-08-14): does tune_fmm's
output transfer to FGS configs now that FastGaussSeidel's trees replay the
(shallower) SOURCE topology?

tune_fmm tunes through the standard two-tree fmm! path (independent target
tree); FastGaussSeidel replays the source topology on the target side, so the
tuner's error model / knob search may not describe the FGS far-field. This
script answers empirically, per rung (R1-R3), by running the two FGS-family
roster configs at the PROPOSED matched-residual target (rel 1e-6) under three
knob candidates and judging every run by its DIRECT-evaluated true residual
against the frozen Dirichlet b:

  prod    p=10, MAC=0.4, leaf=150   (Phase 0 / availability production knobs)
  tuned   per-rung tune_fmm output  (tune.csv 2026-08-14; leaf is small)
  hybrid  tuned p + MAC, leaf=150   (FGS leaf doubles as the dense
                                     near-field block size, which tune_fmm
                                     knows nothing about)

PASS for a candidate = internal convergence AND direct rel_rms <= target
(no new floor from the source-topology target tree). Iteration counts + wall
time decide the knob selection (logged as a judgment in
phase_01_consistency.md, PROPOSED pending Ryan).

Case construction mirrors rotor_hover_solver_phase1_availability.jl (keep in
sync). The Dirichlet source potential is assembled once with the DIRECT
backend and restored before every solve (pure function of geometry+BC), so
t_solve excludes the N^2 source assembly — noted in the CSV.

Run (one process per rung):
  RUNG=R1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi julia --project -t 4 \
      benchmark/rotor_hover_solver_phase1_fgs_check.jl

ENV: RUNG (R1-R3), TARGET_REL (default 1e-6), MAXIT (FGS outer cap, default
300), K_REPS (timed cold solves, default 1), FGMRES_ITMAX (default 500),
FGMRES_MEMORY (default 50). Appends fgs_check.csv under
benchmark/results/phase1/<mode>/.
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
# tune_fmm output per rung (tune.csv rows 2026-08-14, variant="tuned")
const TUNED = Dict(
    "R1" => (17, 0.5, 21),
    "R2" => (17, 0.5, 12),
    "R3" => (18, 0.5, 18),
)
const PROD = (10, 0.4, 150)

rung = get(ENV, "RUNG", "")
haskey(LADDER, rung) || error("RUNG must be one of $(sort(collect(keys(LADDER)))); got \"$rung\"")
msh_name, n_expected = LADDER[rung]
target_rel = parse(Float64, get(ENV, "TARGET_REL", "1e-6"))
maxit_fgs = parse(Int, get(ENV, "MAXIT", "300"))
k_reps = parse(Int, get(ENV, "K_REPS", "1"))
fgmres_itmax = parse(Int, get(ENV, "FGMRES_ITMAX", "500"))
fgmres_memory = parse(Int, get(ENV, "FGMRES_MEMORY", "50"))

# ---- Dirichlet case, mirroring the Phase 1 driver / availability script ----
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

# frozen Dirichlet b (direct source-potential assembly) + cached potential:
# rhs = -potential is a pure function of geometry + BC, so the assembly is
# restored (not recomputed) before every solve
pnl.set_strengths!(rotor)
rotor.potential .= 0
pnl.influence!(rotor, rotor, pnl.DirectBackend(); scalar_potential=true, velocity=false)
potential_frozen = copy(rotor.potential)
b = zeros(rotor.ncells)
pnl.assemble_rhs!(b, rotor)
rms_b = norm(b) / sqrt(length(b))
tol_abs = target_rel * rms_b   # FGS tolerance is ABSOLUTE max-abs; max >= rms
println("$rung Dirichlet case: $(rotor.ncells) panels, rms(b) = $rms_b, " *
        "target rel = $target_rel (FGS abs tol $tol_abs)\n")

solution_column = 2
"Restore the frozen pre-solve state (cold solve; FGS seeds from body.strength)."
function reset_cold!()
    pnl.set_strengths!(rotor)                       # col 1 = BC sources
    rotor.strength[:, solution_column] .= 0
    rotor.velocity .= frozen_velocity
    rotor.potential .= potential_frozen             # rhs = -potential contract
    return nothing
end

p_t, mac_t, leaf_t = TUNED[rung]
p_p, mac_p, leaf_p = PROD
candidates = [
    ("prod",   p_p, mac_p, leaf_p),
    ("tuned",  p_t, mac_t, leaf_t),
    ("hybrid", p_t, mac_t, leaf_p),
]

fgs_kwargs = (; max_iterations=maxit_fgs, tolerance=tol_abs, rlx=1.0,
                shrink=true, recenter=false, inner_iterations=20,
                reverse_pass=false, verbose=false)

header = "rung,mesh_file,n_panels,config,expansion_order,multipole_acceptance," *
    "leaf_size,tolerance_setting,target_rel,iterations,converged_internal," *
    "internal_final_residual,t_solve_min,k_reps,rms_residual_direct," *
    "max_residual_direct,rel_rms_direct,meets_target,rms_b,radius_tol," *
    "threading_mode,julia_threads,commit,fm_commit,date,notes"
csv_path = joinpath(results_dir, "fgs_check.csv")
fresh = !isfile(csv_path) || filesize(csv_path) == 0
io = open(csv_path, "a")
fresh && println(io, header)

r = zeros(rotor.ncells)
function report!(config, p, mac, leaf, tol_setting, iters, conv, internal_final,
                 t_solve, notes)
    x = vec(rotor.strength[:, solution_column])
    rms_d, rmax_d = pnl.true_residual!(r, rotor, x, b; backend=pnl.DirectBackend())
    rel = rms_d / rms_b
    meets = rel <= target_rel
    println("  -> iters=$iters conv=$conv internal_final=$internal_final " *
            "t=$(round(t_solve; digits=2))s direct rel_rms=$rel " *
            (meets ? "MEETS TARGET" : "MISSES TARGET"))
    println(io, join(_csv_cell.([rung, msh_name, rotor.ncells, config, p, mac,
        leaf, tol_setting, target_rel, iters, conv, internal_final, t_solve,
        k_reps, rms_d, rmax_d, rel, meets, rms_b, pnl.FMM_RADIUS_TOL[],
        banner.threading_mode, banner.julia_threads, banner.commit,
        banner.fm_commit, time_string(), notes]), ","))
    flush(io)
    return nothing
end

# ---- FGSSolver candidates ----
for (name, p, mac, leaf) in candidates
    config = "fgs_$name"
    println("--- $config (p=$p, MAC=$mac, leaf=$leaf) ---")
    solver = pnl.FGSSolver(rotor; expansion_order=p, multipole_acceptance=mac,
                           leaf_size=leaf, fgs_kwargs...)
    hist = pnl.ConvergenceHistory(:fgs_maxabs)
    reset_cold!()
    pnl._solve_history!(rotor, solver, hist)        # untimed history solve
    iters = length(hist)
    internal_final = isempty(hist.residual_internal) ? NaN : hist.residual_internal[end]
    conv = internal_final <= tol_abs
    t_solve, _ = min_of_k(() -> pnl._solve!(rotor, solver);
                          k=k_reps, warmup=0, setup! = reset_cold!)
    # the timed solve overwrote strengths from the same cold start; residual
    # below reflects the timed solve's converged state
    report!(config, p, mac, leaf, "abs_maxabs=$tol_abs", iters, conv,
            internal_final, t_solve,
            "FGSSolver inner=20 rlx=1.0 shrink=true maxit=$maxit_fgs; " *
            "t excludes N^2 source-potential assembly (restored, not rebuilt)")
    solver = nothing; GC.gc()
end

# ---- FGMRES + FGSPreconditioner candidates (backend apply = tuned knobs) ----
backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)
for (name, p, mac, leaf) in candidates
    config = "fgmres_fgs_$name"
    println("--- $config (precond p=$p, MAC=$mac, leaf=$leaf; apply=tuned) ---")
    P = pnl.FGSPreconditioner(rotor; sweeps=1, inner_iterations=2,
        expansion_order=p, multipole_acceptance=mac, leaf_size=leaf)
    solver = pnl.KrylovSolver(rotor; method=:fgmres, itmax=fgmres_itmax,
        atol=1e-14, rtol=target_rel, memory=fgmres_memory,
        backend=backend_apply, preconditioner=P, record_history=true)
    reset_cold!()
    pnl._solve!(rotor, solver)                      # untimed history solve
    iters = solver.niter
    conv = solver.solved
    internal_final = isempty(solver.history.residual_internal) ? NaN :
                     solver.history.residual_internal[end]
    solver.record_history = false
    t_solve, _ = min_of_k(() -> pnl._solve!(rotor, solver);
                          k=k_reps, warmup=0, setup! = reset_cold!)
    report!(config, p, mac, leaf, "rtol=$target_rel;atol=1e-14", iters, conv,
            internal_final, t_solve,
            "fgmres memory=$fgmres_memory itmax=$fgmres_itmax sweeps=1 inner=2; " *
            "apply backend p=$p_t/MAC=$mac_t/leaf=$leaf_t (tuned); " *
            "t excludes N^2 source-potential assembly (restored, not rebuilt)")
    P = nothing; solver = nothing; GC.gc()
end

close(io)
println("\n$rung FGS check done. Rows appended to $csv_path")
