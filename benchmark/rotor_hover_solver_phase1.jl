#=##############################################################################
BRAINSTORM 021 Phase 1 — consistency driver: no-wake frozen-RHS solves on the
mesh ladder (ledger.md).

Case definition (Ryan ruling 2026-08-13): matches
examples/rotor_hover_pressure_comparison.jl (RHPC), the production rotor
pipeline — Dirichlet formulation on capped/watertight captess4 meshes,
kernel Union{ConstantSource, VortexRing}, DBC=true, core_size_panel=R·1e-10,
core_size_targets=1e-3, kernelcutoff=R·1e-13, rho=1.179, RPM=6000; DJI
conventions x axial / y radial (scale by dim 2, rotation about −x, thrust
= −force_x). Shedding follows RHPC's recipe: full end_node-anchored TE trace
(bbox only separates the blades), explicit root clip at r/R=0.1, and the
circumferential (cap-wrap) guard — all traced on the *constructed* noshedding
body per the CLAUDE.md invariant.

Stages (per rung, one process per rung — dense G at R3 is ~6.6 GB):
  1. Backslash reference with separately timed assembly / factorization /
     solve (the matrix-ful cost check that sets the ladder's cost ceiling);
  2. FMM residual floor (ruling 11c): the reference solution's true residual
     through the FastMultipole operator over an expansion-order sweep
     (floor.csv). The frozen RHS b is captured Dirichlet-correctly: sources
     set from the BC, source potential assembled with the DIRECT backend,
     b = −potential at solver invocation (see runtests_unit_instrumentation.jl
     and the assemble_rhs! docstring); the Dirichlet solution is strength
     column 2.
  3. TUNE=1: FastMultipole.tune_fmm per rung at the campaign tolerance, then a
     one-at-a-time perturbation check (MAC ±0.1, p ±1, leaf ×/÷2) verifying
     near-optimality (tune.csv).

Run:
  RUNG=R1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi julia --project -t 4 \
      benchmark/rotor_hover_solver_phase1.jl

ENV knobs: RUNG (required: R1..R7), K_REPS (default 3 for local dev; published
runs use >= 5 per ruling 5), FMM_ORDERS (default "4,6,8,10,12,14"),
FMM_ACCEPTANCE (default 0.4), FMM_LEAF_SIZE (default 20), TUNE (0/1),
TUNE_TARGET_REL (default 1e-6), RADIUS_TOL (sets pnl.FMM_RADIUS_TOL[]; "Inf"
reproduces pre-2026-08-13 geometric radii for A/B), MESH_OVERRIDE (diagnostics
only: same-N different mesh file), HARDWARE_TAG.

CSVs append (never truncate — Phase 0 smoke-harness gotcha) under
benchmark/results/phase1/<threading_mode>/: runs.csv (reference solves),
floor.csv (floor sweep), tune.csv (tuning + perturbations).
=###############################################################################

import FLOWPanel as pnl
import FastMultipole
include(joinpath(pnl.examples_path, "dji9443_trailing_edge.jl"))
include(joinpath(@__DIR__, "common.jl"))

using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm, cross, lu!

################################################################################
# Banner + output layout
################################################################################

banner = assert_and_banner()

results_dir = joinpath(@__DIR__, "results", "phase1", banner.threading_mode)
mkpath(results_dir)
write(joinpath(results_dir, "banner.txt"), banner.text * "\n")

"Open `path` for appending, writing `header` first iff the file is new."
function open_append(path::AbstractString, header::AbstractString)
    fresh = !isfile(path) || filesize(path) == 0
    io = open(path, "a")
    fresh && println(io, header)
    return io
end

################################################################################
# Frozen mesh ladder (ledger.md; regenerate via scripts/generate_dji9443_mesh.sh
# with `U W capped --cap-tess 4`)
################################################################################

const LADDER = Dict(
    "R1" => ("dji9443_20260813_23_73_capped_captess4.msh",     8016),
    "R2" => ("dji9443_20260813_33_105_capped_captess4.msh",   15760),
    "R3" => ("dji9443_20260813_45_145_capped_captess4.msh",   28752),
    "R4" => ("dji9443_20260813_65_209_capped_captess4.msh",   58192),
    "R5" => ("dji9443_20260813_89_289_capped_captess4.msh",  108240),
    "R6" => ("dji9443_20260813_125_409_capped_captess4.msh", 212108),
    "R7" => ("dji9443_20260813_177_577_capped_captess4.msh", 419276),
)

rung = get(ENV, "RUNG", "")
haskey(LADDER, rung) || error("RUNG must be one of $(sort(collect(keys(LADDER)))); got \"$rung\"")
msh_name, n_expected = LADDER[rung]
# diagnostics only: run a rung's protocol on a different mesh file (same N)
msh_name = get(ENV, "MESH_OVERRIDE", msh_name)
msh_file = joinpath(pnl.examples_path, "data", msh_name)

################################################################################
# Case parameters (RHPC / production rotor pipeline)
################################################################################

magVinf = 0.0001
rho     = 1.179
RPM     = 6000
R       = 0.119
Vinf    = magVinf * [1.0, 0.0, 0.0]        # axial (x) freestream, effectively hover
nt      = 36
dt      = 60 / RPM / nt
k_reps  = parse(Int, get(ENV, "K_REPS", "3"))   # local dev; published runs k >= 5 (ruling 5)

core_size_panel   = R * 1e-10
core_size_targets = 1e-3
kernelcutoff         = R * 1e-13
init_Das_eta_kinematic = 0.2
shedding_r_over_R    = 0.1
radial_dimension     = 2

fmm_orders = parse.(Int, split(get(ENV, "FMM_ORDERS", "4,6,8,10,12,14"), ","))
fmm_acceptance = parse(Float64, get(ENV, "FMM_ACCEPTANCE", "0.4"))
fmm_leaf_size = parse(Int, get(ENV, "FMM_LEAF_SIZE", "20"))

if haskey(ENV, "RADIUS_TOL")
    pnl.FMM_RADIUS_TOL[] = parse(Float64, ENV["RADIUS_TOL"])
end
radius_tol = pnl.FMM_RADIUS_TOL[]
println("FMM radius-inflation tolerance: $radius_tol")

################################################################################
# Case construction (RHPC shedding recipe on the constructed body)
################################################################################

println("Rung $rung: $msh_file")
te_indices_1, te_indices_2 = find_dji9443_trailing_edge_indices(msh_file; watertight=true)

msh = pnl.read_gmsh(msh_file)
nodes0, cells0 = pnl.meshes2nodes_cells(msh)
nodes0 .*= R / maximum(nodes0[radial_dimension, :])   # y is radial for the DJI9443 meshes

kernel = Union{pnl.ConstantSource, pnl.VortexRing}
DBC = true
body_kwargs = (; core_size=core_size_panel, core_size_panel,
                 core_size_targets, kernelcutoff,
                 semiinfinite_wake=false, watertight=true, DBC)
rotor0 = pnl.RigidWakeBody{kernel}(nodes0, cells0, pnl.noshedding; body_kwargs...)

"Half-space bbox separating one blade (per RHPC / rotor_hover_convergence.jl)."
function make_shedding_bbox(nodes, seed_nodes, radial_dimension, R, cutoff_r_over_R)
    radial_midpoint = sum(nodes[radial_dimension, seed_nodes]) / length(seed_nodes)
    radial_sign = sign(radial_midpoint)
    radial_sign == 0 && error("Seed edge lies on the rotor axis; cannot determine shedding side")
    lower = [minimum(nodes[i, :]) for i in 1:size(nodes, 1)]
    upper = [maximum(nodes[i, :]) for i in 1:size(nodes, 1)]
    padding = max(sqrt(eps(eltype(nodes))) * R, R * 1e-6)
    lower .-= padding
    upper .+= padding
    radial_cutoff = cutoff_r_over_R * R
    if radial_sign > 0
        lower[radial_dimension] = radial_cutoff - padding
    else
        upper[radial_dimension] = -radial_cutoff + padding
    end
    return (pnl.SVector{3}(lower...), pnl.SVector{3}(upper...))
end

"RHPC modeling root clip: drop shed edges whose midpoint radius is inside the clip."
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

"RHPC cap-wrap guard: count shed edges running circumferentially (radial ratio < tol)."
function count_circumferential_shedding(nodes, shedding, cells, radial_dimension; tol=0.3)
    n = 0
    for j in axes(shedding, 2)
        p, nia, nib = shedding[1, j], shedding[2, j], shedding[3, j]
        na, nb = cells[nia, p], cells[nib, p]
        delta = nodes[:, na] - nodes[:, nb]
        len = sqrt(sum(abs2, delta))
        len > 0 && abs(delta[radial_dimension]) / len < tol && (n += 1)
    end
    return n
end

sheddings = map((te_indices_1, te_indices_2)) do te
    bbox = make_shedding_bbox(rotor0.nodes, te[1:2], radial_dimension, R, 0.0)
    full = pnl.calc_shedding_from_seed(rotor0.nodes, rotor0.cells, te[1], te[2];
        bbox, end_node=te[3], normal_jump_tol=0.2, max_turn_angle=pi/3)
    clipped = clip_shedding_root(rotor0.nodes, full, rotor0.cells, radial_dimension, R,
        shedding_r_over_R)
    nbad = count_circumferential_shedding(rotor0.nodes, clipped, rotor0.cells, radial_dimension)
    nbad == 0 || error("shedding includes $nbad circumferential edges — TE chain wrapped onto a cap")
    println("  TE traced: $(size(full, 2)) edges -> $(size(clipped, 2)) after root clip at r/R $(shedding_r_over_R)")
    clipped
end

rotor = pnl.RigidWakeBody{kernel}(rotor0.nodes, rotor0.cells, collect(sheddings);
            ensure_winding=true, body_kwargs...)
rotor.ncells == n_expected ||
    error("Rung $rung: expected $n_expected panels, got $(rotor.ncells)")

frames = pnl.ReferenceFrame(rotor;
    origin = SVector{3}(0.0, 0.0, 0.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(-1.0, 0.0, 0.0),
    ω = 2*pi * RPM/60,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = [1])
pnl.initialize_Das!((rotor,), frames, t -> Vinf, 0.0, dt;
    set_Das_eta_kinematic=init_Das_eta_kinematic)

# Frozen rotor-like RHS: apparent velocity Vinf - ω×r at each control point.
ω_vec = SVector{3}(-2*pi * RPM/60, 0.0, 0.0)   # ω·ω_axis
pnl.calc_normals!(rotor)
pnl.calc_controlpoints!(rotor)
for i in 1:rotor.ncells
    r = SVector{3}(rotor.controlpoints[1, i], rotor.controlpoints[2, i], rotor.controlpoints[3, i])
    rotor.velocity[:, i] .= Vinf .- cross(ω_vec, r)
end
frozen_velocity = copy(rotor.velocity)

################################################################################
# Frozen Dirichlet RHS + Backslash reference (timed assembly/factorize/solve)
################################################################################

# Dirichlet RHS capture (assemble_rhs! contract): sources from the BC, source
# potential assembled with the DIRECT backend, b = -potential at solver
# invocation. All configs are judged against this one frozen b.
pnl.set_strengths!(rotor)
rotor.potential .= 0
t_rhs = @elapsed pnl.influence!(rotor, rotor, pnl.DirectBackend();
    scalar_potential=true, velocity=false)
b = zeros(rotor.ncells)
pnl.assemble_rhs!(b, rotor)
rms_b = norm(b) / sqrt(length(b))
println("RHS: rms(b) = $rms_b (direct source-potential assembly $t_rhs s)")

println("\n--- Backslash reference (matrix-ful cost check) ---")
G = zeros(rotor.ncells, rotor.ncells)
t_assembly = @elapsed pnl._G!(G, rotor, rotor; core_size=rotor.core_size_panel)
t_factorize = @elapsed lu!(G)
println("assembly $t_assembly s, factorization $t_factorize s " *
        "(dense G = $(round(8 * rotor.ncells^2 / 2^30; digits=2)) GiB)")
G = nothing; GC.gc()   # timing copy freed before the solver builds its own G

# Dirichlet solution lives in strength column 2 (col 1 = BC-set sources)
solution_column = 2
solver = pnl.Backslash(rotor)

# Reference solution against exactly the frozen b: rotor.potential still holds
# the direct source-potential assembly captured above, and _solve! consumes it
# (the full solve! path would reassemble internally).
pnl._solve!(rotor, solver)
x_ref = copy(vec(rotor.strength[:, solution_column]))
r = zeros(rotor.ncells)
rms_ref, rmax_ref = pnl.true_residual!(r, rotor, x_ref, b; backend=pnl.DirectBackend())
println("reference solution: direct true residual rms = $rms_ref, max = $rmax_ref")

# Timed benchmark solves through the production tuple path (DirectBackend
# default). NOTE: for Dirichlet bodies each solve! reassembles the source
# potential internally, so t_solve_min includes that N² assembly, not just the
# triangular solves — recorded as such in the notes.
reset_cold! = () -> (rotor.strength[:, solution_column] .= 0;
                     rotor.velocity .= frozen_velocity; nothing)
solve_once! = () -> pnl.solve!((rotor,), (solver,))
reset_cold!()
t_cold_first = @elapsed solve_once!()
t_solve_min, _ = min_of_k(solve_once!; k=k_reps, warmup=1, setup! = reset_cold!)
reset_cold!()
alloc_solve = @allocated solve_once!()
println("timed solve (incl. Dirichlet source-potential assembly): $t_solve_min s (min of $k_reps)")

io_runs = open_append(joinpath(results_dir, "runs.csv"), runs_csv_header())
write_run_row!(io_runs;
    run_id="backslash_ref_$(rung)_$(banner.threading_mode)",
    phase="phase1_costcheck", solver_config="backslash_ref",
    mesh_file=msh_name, n_panels=rotor.ncells,
    threading_mode=banner.threading_mode,
    julia_threads=banner.julia_threads, blas_threads=banner.blas_threads,
    t_assembly, t_factorize, t_rhs, t_solve_min, k_reps,
    rms_residual=rms_ref, max_residual=rmax_ref,
    mem_state_bytes=solver_state_bytes(solver), alloc_solve_bytes=alloc_solve,
    commit=banner.commit, fm_commit=banner.fm_commit,
    julia_version=banner.julia_version, hardware_tag=banner.hardware_tag,
    filament_reg=banner.filament_reg,
    solver_settings=settings_string(pnl._solver_metadata_dict(solver)),
    backend_settings="",
    notes="phase1 Dirichlet (RHPC case); tuple block-GS path (single inner solve " *
          "on a single body since the 021 Phase-3 one-body early break, incl. " *
          "direct source-potential assembly); " *
          "cold_first $(t_cold_first)s; rms_b $(rms_b); radius_tol $(radius_tol)")
close(io_runs)

################################################################################
# FMM residual floor (ruling 11c): reference solution through the FMM operator
################################################################################

floor_backend(p, mac, leaf) = pnl.FastMultipoleBackend(;
    expansion_order=p, multipole_acceptance=mac, leaf_size=leaf)

"Evaluate the floor (reference residual through `backend`); returns (rms, max, t)."
function floor_point(backend)
    t_eval = @elapsed rms_p, rmax_p = pnl.true_residual!(r, rotor, x_ref, b; backend)
    return rms_p, rmax_p, t_eval
end

println("\n--- FMM residual floor sweep (MAC=$fmm_acceptance, leaf=$fmm_leaf_size) ---")
floor_header = "rung,mesh_file,n_panels,expansion_order,multipole_acceptance," *
    "leaf_size,rms_b,rms_residual,max_residual,rel_rms,t_eval,radius_tol," *
    "threading_mode,julia_threads,commit,fm_commit,date"
io_floor = open_append(joinpath(results_dir, "floor.csv"), floor_header)

# direct-backend row (expansion_order empty): the reference residual itself
println(io_floor, join(_csv_cell.([rung, msh_name, rotor.ncells, "", "", "",
    rms_b, rms_ref, rmax_ref, rms_ref / rms_b, "", radius_tol,
    banner.threading_mode, banner.julia_threads, banner.commit,
    banner.fm_commit, time_string()]), ","))

for p in fmm_orders
    rms_p, rmax_p, t_eval = floor_point(floor_backend(p, fmm_acceptance, fmm_leaf_size))
    rel = rms_p / rms_b
    println("p = $(lpad(p, 2)): rms = $rms_p  (rel $rel)  max = $rmax_p  " *
            "[$(round(t_eval; digits=2)) s incl. tree build]")
    println(io_floor, join(_csv_cell.([rung, msh_name, rotor.ncells, p,
        fmm_acceptance, fmm_leaf_size, rms_b, rms_p, rmax_p, rel, t_eval,
        radius_tol, banner.threading_mode, banner.julia_threads, banner.commit,
        banner.fm_commit, time_string()]), ","))
end
close(io_floor)

################################################################################
# TUNE=1: tune_fmm + one-at-a-time perturbation double-check (tune.csv)
################################################################################

if get(ENV, "TUNE", "0") == "1"
    target_rel = parse(Float64, get(ENV, "TUNE_TARGET_REL", "1e-6"))
    # Dirichlet apply is potential-based; 0.1 safety factor keeps the FMM error
    # subdominant to the matched-residual target
    error_tolerance = FastMultipole.PowerAbsolutePotential(0.1 * target_rel * rms_b)
    println("\n--- tune_fmm (PowerAbsolutePotential($(0.1 * target_rel * rms_b))) ---")

    # tune with the reference solution loaded (tuning depends on strengths)
    rotor.strength[:, solution_column] .= x_ref
    t_tune = @elapsed tuned, _cache = FastMultipole.tune_fmm((rotor,), (rotor,);
        error_tolerance, multipole_acceptances=range(0.3, 0.7, step=0.1),
        max_expansion_order=20,
        scalar_potential=true, gradient=false, hessian=false, verbose=true)
    println("tuned: $tuned ($(round(t_tune; digits=1)) s); cache discarded " *
            "(layout-bound; backend rebuilds trees per call)")

    tune_header = "rung,mesh_file,n_panels,variant,expansion_order," *
        "multipole_acceptance,leaf_size,t_tune,t_eval,rms_residual,rel_rms," *
        "dt_rel_vs_tuned,meets_target,radius_tol,threading_mode,julia_threads," *
        "commit,fm_commit,date"
    io_tune = open_append(joinpath(results_dir, "tune.csv"), tune_header)

    p0 = tuned.expansion_order
    mac0 = tuned.multipole_acceptance
    leaf0 = tuned.leaf_size_source[1]   # SVector, one entry per source system
    variants = [
        ("tuned",       p0,     mac0,                        leaf0),
        ("mac_minus",   p0,     max(mac0 - 0.1, 0.3),        leaf0),
        ("mac_plus",    p0,     min(mac0 + 0.1, 0.8),        leaf0),
        ("p_minus",     p0 - 1, mac0,                        leaf0),
        ("p_plus",      p0 + 1, mac0,                        leaf0),
        ("leaf_half",   p0,     mac0,                        max(leaf0 ÷ 2, 1)),
        ("leaf_double", p0,     mac0,                        2 * leaf0),
    ]
    t_tuned_eval = Ref(NaN)
    for (variant, p, mac, leaf) in variants
        backend = floor_backend(p, mac, leaf)
        # min-of-k timing (tree build included every call, as in production use)
        t_eval, _ = min_of_k(() -> pnl.true_residual!(r, rotor, x_ref, b; backend);
                             k=k_reps, warmup=1)
        rms_v, _, _ = floor_point(backend)
        rel = rms_v / rms_b
        meets = rel <= target_rel
        variant == "tuned" && (t_tuned_eval[] = t_eval)
        dt_rel = t_eval / t_tuned_eval[] - 1
        println("  $(rpad(variant, 12)) p=$(lpad(p,2)) mac=$mac leaf=$(lpad(leaf,4)): " *
                "t=$(round(t_eval; digits=3))s ($(round(100dt_rel; digits=1))% vs tuned), " *
                "rel_rms=$rel $(meets ? "" : "[misses target]")")
        println(io_tune, join(_csv_cell.([rung, msh_name, rotor.ncells, variant,
            p, mac, leaf, variant == "tuned" ? t_tune : "", t_eval, rms_v, rel,
            dt_rel, meets, radius_tol, banner.threading_mode,
            banner.julia_threads, banner.commit, banner.fm_commit,
            time_string()]), ","))
    end
    close(io_tune)

    println("tune verdict: inspect tune.csv — tuned point is near-optimal if no " *
            "target-meeting perturbation is >10-15% faster")
end

println("\nRung $rung done. CSVs appended under $results_dir")
