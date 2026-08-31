#=##############################################################################
BRAINSTORM 021 Phase 1 — agreement tables vs Backslash at matched residual
(exit criterion 1): strength-difference norms + CT from the three pressure
monitors (PressureLaplace+Force, PressureBernoulli+Force, KuttaJoukowski),
per rung x config, through the production `steady!` path (skeleton =
examples/rotor_hover_convergence.jl, monitors per
examples/rotor_hover_pressure_comparison.jl).

DESCOPED 2026-08-24 (Ryan): Phase 1's remaining purpose is exactly this file
— solve the ladder, up to and including R7, and verify that every solver
agrees on the SOLUTION. Tuning moved to Phase 2
(rotor_hover_solver_phase2_tune.jl), which is what publishes benchmarks; wall
clock recorded here is incidental and must not be quoted as a benchmark. Set
NFCACHE=<GiB> to give the Krylov configs a near-field cache purely for speed
— it changes no reported metric (the cached and kernel near-fields agree to
rtol 1e-12) and is skipped if it does not fit the budget.

ERROR METRIC IS THE BOUNDARY CONDITION, NOT A REFERENCE (Ryan 2026-08-24).
Each config is judged by `bc_rel_l2 = ||Ax-b||/rms_b` on that config's OWN
solution — reference-free and identically defined at every rung. Evaluated
through the FMM under an explicit error tolerance of `0.1 * target_rel * rms_b`
(an order of magnitude tighter than the criterion under test) and reported
alongside the FMM's own `bc_certified` flag; `RESIDUAL_BACKEND=direct` gives
the exact O(N^2) evaluation, which is only tractable at the small rungs. Cross-config CONSISTENCY is then measured
against the ENSEMBLE (agreement_spread.csv), not against a nominated run.

This is what makes R6/R7 tractable at all: a dense `backslash_ref` would be
8N^2 = 335 GiB at R6 and 1.31 TiB at R7. Backslash is still worth running
where it FITS (25 GiB at R4) as an independent non-iterative cross-check, but
it is no longer required, no longer privileged, and config ORDER no longer
matters. REFERENCE=<config> optionally restores the old reference-relative
columns; leave it unset.

Design (2026-08-14 autonomous session; judgments logged in
phase_01_consistency.md):
- (SUPERSEDED 2026-08-24, see above) Reference = `backslash_ref` run through
  the SAME steady! path; agreement metrics were vs its strengths/CTs. The
  primary metric is now the reference-free BC residual, and consistency is
  measured against the ensemble.
- All configs share one post-processing/monitor backend (`backend_system` =
  the per-rung TUNED apply knobs) so CT deltas reflect only the solved
  strengths, and share backend_solve=DirectBackend for the Dirichlet
  source-potential RHS assembly (exact rhs for every config; matrix-free
  operators still run through each solver's own internal backend).
- Iterative-config settings = the frozen/PROPOSED Phase 1 values
  (decision_rules.md): Krylov rtol=1e-6/atol=1e-14/itmax=500/memory=50;
  FGS p=10/MAC0.4/leaf150, inner=20, rlx=1.0, tol=1e-6*rms_b(rung) abs.
- Each config also gets a direct-evaluated true residual vs the frozen b of
  the phase1 driver contract (matched-residual verification).
- Monitors are freshly constructed per config; steady Bernoulli
  (unsteady=false, correct_kuttacondition=false) mirrors
  rotor_hover_convergence.jl — its known moving-body validity caveat is
  acceptable here because the metric is the CROSS-CONFIG DELTA at identical
  recovery, not absolute CT.

Run (one process per rung):
  RUNG=R1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi julia --project -t 4 \
      benchmark/rotor_hover_solver_phase1_agreement.jl

ENV: RUNG (R1-R3), CONFIGS (comma list to subset; default all).
Appends agreement.csv under benchmark/results/phase1/<mode>/.
=###############################################################################

import FLOWPanel as pnl
import FastMultipole
include(joinpath(pnl.examples_path, "dji9443_trailing_edge.jl"))
include(joinpath(@__DIR__, "common.jl"))

using FLOWPanel.FastMultipole.StaticArrays
using LinearAlgebra: norm, cross

banner = assert_and_banner()
# PER_RUNG_DIR=1 (HPC): isolate each rung's CSVs in results/phase1/<mode>/<rung>/
# — concurrent per-rung jobs on different nodes append over NFS, and
# non-atomic appends CLOBBER shared files (observed 2026-08-18: R1's
# staircase rows destroyed by R2–R4 writers). This driver holds agreement.csv
# open for hours and writes at the end, so it is the same hazard exactly.
# Merge per-rung dirs with p021_merge.jl; local sequential runs keep the flat
# layout. (Pattern copied verbatim from phase1_case.jl.)
_rung_sub = get(ENV, "PER_RUNG_DIR", "0") == "1" ? get(ENV, "RUNG", "") : ""
results_dir = joinpath(@__DIR__, "results", "phase1", banner.threading_mode,
                       _rung_sub)
mkpath(results_dir)

# Full frozen ladder (matches phase1_case.jl). Extended past R3 on 2026-08-24:
# the descope makes R7 this driver's headline target.
const LADDER = Dict(
    "R1" => ("dji9443_20260813_23_73_capped_captess4.msh",    8016),
    "R2" => ("dji9443_20260813_33_105_capped_captess4.msh",  15760),
    "R3" => ("dji9443_20260813_45_145_capped_captess4.msh",  28752),
    "R4" => ("dji9443_20260813_65_209_capped_captess4.msh",  58192),
    "R5" => ("dji9443_20260813_89_289_capped_captess4.msh", 108240),
    "R6" => ("dji9443_20260813_125_409_capped_captess4.msh", 212108),
    "R7" => ("dji9443_20260813_177_577_capped_captess4.msh", 419276),
)
# Apply knobs for the SHARED post-processing/monitor backend. Since the
# 2026-08-24 descope (below) these only have to be ACCURATE, not fast: this
# driver measures cross-config deltas at identical recovery, and every config
# sees the same backend_system, so a suboptimal leaf costs wall clock and
# changes no reported number. Untuned rungs therefore fall back to a reference
# triple instead of failing.
const REF_KNOBS = (17, 0.5, 21)
const TUNED = Dict(   # per-rung tuned apply knobs (tune.csv 2026-08-14)
    "R1" => (17, 0.5, 21),
    "R2" => (17, 0.5, 12),
    "R3" => (18, 0.5, 18),
)

# The criterion under test. Iterative configs solve to this, and the BC metric
# is evaluated an order of magnitude tighter (safety=0.1 below).
const target_rel = parse(Float64, get(ENV, "TARGET_REL", "1e-6"))

rung = get(ENV, "RUNG", "")
haskey(LADDER, rung) || error("RUNG must be one of $(sort(collect(keys(LADDER)))); got \"$rung\"")
msh_name, n_expected = LADDER[rung]

# ---- Dirichlet case, mirroring the Phase 1 driver (keep in sync) ----
magVinf = 0.0001; rho = 1.179; RPM = 6000; R = 0.119
Vinf = magVinf * [1.0, 0.0, 0.0]
dt = 60 / RPM / 36
core_size_panel = R * 1e-10
core_size_targets = 1e-3
kernelcutoff = R * 1e-13
radial_dimension = 2
shedding_r_over_R = 0.1
axial_dimension = 1

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

omega = 2 * pi * RPM / 60
frames = pnl.ReferenceFrame(rotor;
    origin=SVector{3}(0.0,0.0,0.0), v=SVector{3}(0.0,0.0,0.0),
    ω_axis=SVector{3}(-1.0,0.0,0.0), ω=omega,
    R=SMatrix{3,3}(1.0,0,0, 0,1.0,0, 0,0,1.0),
    name="vehicle", child_index=Int[], dependent_index=[1])
pnl.initialize_Das!((rotor,), frames, t -> Vinf, 0.0, dt; set_Das_eta_kinematic=0.2)

# frozen b for matched-residual verification (phase1 driver contract)
ω_vec = SVector{3}(-omega, 0.0, 0.0)
pnl.calc_normals!(rotor); pnl.calc_controlpoints!(rotor)
for i in 1:rotor.ncells
    rr = SVector{3}(rotor.controlpoints[:, i]...)
    rotor.velocity[:, i] .= Vinf .- cross(ω_vec, rr)
end
pnl.set_strengths!(rotor)
rotor.potential .= 0
pnl.influence!(rotor, rotor, pnl.DirectBackend(); scalar_potential=true, velocity=false)
# Frozen kinematic state. `steady!` runs post-processing monitors that
# OVERWRITE rotor.velocity (and leave core_size at the TARGETS offset), and
# `bc_error!` recomputes the BC source strengths from body.velocity — so
# without restoring these first, every BC number would be computed against the
# wrong sigma. Observed before this fix: bc_rel_l2 = 0.985 for a converged
# krylov_gmres, i.e. the metric read "BC completely unsatisfied".
frozen_velocity = copy(rotor.velocity)
b = zeros(rotor.ncells)
pnl.assemble_rhs!(b, rotor)
rms_b = norm(b) / sqrt(length(b))
rotor.potential .= 0
println("$rung Dirichlet case: $(rotor.ncells) panels, rms(b) = $rms_b\n")

p_t, mac_t, leaf_t = get(TUNED, rung, REF_KNOBS)
haskey(TUNED, rung) || println("no TUNED[$rung]; using reference apply knobs " *
    "$REF_KNOBS for backend_system (accuracy matters here, speed does not)")
backend_system = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)
fgs_tol = target_rel * rms_b

# NFCACHE=<GiB>: opportunistic near-field cache on the Krylov configs (Ryan
# 2026-08-24, "cache phase 1 if it improves speed, but it's not important to
# the outcome"). Purely a speed lever — the cached and kernel near-fields agree
# to rtol 1e-12, so no reported agreement metric moves. Skipped, with a notice,
# when the cache does not fit: see benchmark/nfcache_feasibility.jl for what a
# rung actually needs (measured: 0.20 GiB at R1 but 78 GiB at R7 even at the
# smallest leaf, so R6/R7 need a large budget or none at all).
nfcache_gib = parse(Float64, get(ENV, "NFCACHE", "0"))
nfcache_kw = if nfcache_gib > 0
    est = let plan = pnl.FastMultipole.FmmPlan((rotor,), (rotor,);
                scalar_potential=true, gradient=false, hessian=false,
                expansion_order=p_t, multipole_acceptance=mac_t,
                leaf_size_source=leaf_t, shrink=true)
        e = pnl.FastMultipole.estimate_nearfield_cache(plan.target_tree,
            plan.source_tree, plan.direct_list, plan.derivatives_switches,
            (rotor,); sample=false)
        plan = nothing; GC.gc(); e
    end
    budget = round(Int, nfcache_gib * 1024^3)
    if est.bytes > budget
        println("NFCACHE=$nfcache_gib GiB requested but $rung needs " *
                "$(round(est.bytes/1024^3; digits=2)) GiB at these knobs — " *
                "running UNCACHED (this changes no reported metric)")
        (;)
    else
        println("near-field cache ON for the Krylov configs: " *
                "$(round(est.bytes/1024^3; digits=2)) GiB of $nfcache_gib GiB")
        (; cache_nearfield=true, persistent_plan=true,
           nearfield_cache_max_bytes=budget)
    end
else
    (;)
end

krylov_kw = (; itmax=500, atol=1e-14, rtol=target_rel, memory=50,
             backend=backend_system, nfcache_kw...)

# ---- exit criterion: the BOUNDARY CONDITION, not a reference ---------------
# Ryan 2026-08-24: "we should be measuring error based on the boundary
# condition - not a reference. That way, it doesn't matter."
#
# He is right, and the instrument already exists: `true_residual!` below
# evaluates ||Ax - b|| through the DirectBackend, i.e. the EXACT Dirichlet BC
# residual of each config's own solution, with no reference solution anywhere.
# `rel_rms_direct = rms_d / rms_b` is therefore the primary metric, and it is
# defined identically at R1 and at R7.
#
# This dissolves the blocker that a dense `backslash_ref` is 335 GiB at R6 and
# 1.31 TiB at R7: no config needs to be exact for the others to be checked.
# A Backslash config is still worth running where it FITS, as an independent
# non-iterative cross-check, but it is no longer required and no longer
# privileged, and the config order no longer matters.
#
# Cross-config CONSISTENCY is then also reference-free: every config's solution
# is persisted and the spread is computed against the ENSEMBLE (mean and
# pairwise) in a summary pass at the end, instead of against one nominated run.
# See agreement_spread.csv.
#
# Caveat to keep in mind when reading the spread: at matched BC residual the
# expected solution difference still scales with the operator's condition
# number, so a small residual bounds the solution difference only up to
# cond(A). The spread table is what actually demonstrates consistency; the
# residual demonstrates correctness.
const REFERENCE = get(ENV, "REFERENCE", "")   # "" = no privileged reference
# :fmm (default, certified 10x tighter than target_rel) or :direct (exact,
# O(N^2) — only tractable at the small rungs, use as a cross-check)
const residual_backend = Symbol(get(ENV, "RESIDUAL_BACKEND", "fmm"))
residual_backend in (:fmm, :direct) ||
    error("RESIDUAL_BACKEND must be fmm or direct; got $residual_backend")

configs = [
    ("backslash_ref",     () -> pnl.Backslash(rotor)),
    ("backslash_coupled", () -> pnl.BackslashCoupled((rotor,))),
    ("krylov_gmres",      () -> pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...)),
    ("krylov_jacobi",     () -> pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...,
                                preconditioner=pnl.FastMultipole.JacobiPreconditioner((rotor,); cell_size=R/4))),
    # max_pattern_entries raised from the 512N default: the Barba direct list
    # at leaf 10 / MAC 1.0 grows with refinement on this sliver-heavy mesh.
    # Measured pattern density (2026-08-22, geometry-only probe over the frozen
    # ladder): R1 487/row, R6 2060/row, R7 3691/row — i.e. ~sqrt(N) per row,
    # ~N^1.5 total. The former 2048N limit was set at R3 and undershot R6 by
    # 0.6% (2059.5 vs 2048), which killed p1-table-R6-multi after 4 h. 8192N
    # clears R7 with ~2.2x headroom; it is a guard only, nothing is allocated
    # from it. NOTE the superlinear growth: ILU pattern memory is ~23 GiB of
    # entries at R7 and will not scale past the ladder.
    ("krylov_ilu",        () -> pnl.KrylovSolver(rotor; method=:gmres, krylov_kw...,
                                preconditioner=pnl.ILUPreconditioner(rotor; leaf_size=10, multipole_acceptance=1.0,
                                    max_pattern_entries=8192 * rotor.ncells))),
    ("fgs",               () -> pnl.FGSSolver(rotor; max_iterations=300, tolerance=fgs_tol,
                                rlx=1.0, expansion_order=10, multipole_acceptance=0.4,
                                leaf_size=150, shrink=true, recenter=false,
                                inner_iterations=20, reverse_pass=false, verbose=false)),
    ("fgmres_fgs",        () -> pnl.KrylovSolver(rotor; method=:fgmres, krylov_kw...,
                                preconditioner=pnl.FGSPreconditioner(rotor; sweeps=1,
                                    inner_iterations=2, expansion_order=10,
                                    multipole_acceptance=0.4, leaf_size=150))),
]
if haskey(ENV, "CONFIGS")
    # split on [:,] like phase2.jl: the wrappers promote the colon form
    # (a comma in --export=ALL,CONFIGS=a,b breaks sbatch's own parsing), and
    # a "," -only split would silently select NOTHING from a colon list
    want = split(ENV["CONFIGS"], r"[:,]")
    configs = [c for c in configs if c[1] in want]
end

header = "rung,mesh_file,n_panels,config,t_steady,niter,rms_residual," *
    "bc_rel_l2,bc_certified,residual_backend," *
    "x_relL2_vs_ref,x_relmax_vs_ref,CT_laplace,CT_bernoulli," *
    "CT_kj,dCT_laplace_pct,dCT_bernoulli_pct,dCT_kj_pct,rms_b,radius_tol," *
    "threading_mode,julia_threads,commit,fm_commit,filament_reg,date,notes"
const N_COLS = length(split(header, ","))
csv_path = joinpath(results_dir, "agreement.csv")
fresh = !isfile(csv_path) || filesize(csv_path) == 0
if !fresh
    first(eachline(csv_path)) == header || error(
        "$csv_path has a pre-2026-08-24 schema (the BC-metric columns are " *
        "new and the reference-relative ones are now optional); move it aside")
end
io = open(csv_path, "a")
fresh && println(io, header)

"Emit one CSV row, asserting it matches the header width."
function emit_row!(cells)
    length(cells) == N_COLS || error("agreement.csv row has $(length(cells)) " *
        "cells but the header has $N_COLS")
    println(io, join(_csv_cell.(cells), ","))
    flush(io)
end

make_monitors() = begin
    plap = pnl.PressureLaplace(rotor, rho;
        acceleration_form=:material_derivative, gradient_mode=:edge_difference,
        verbose=false)
    flap = pnl.ForceMonitor(1, 1; i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2R, 1),
        correct_kuttacondition=false, verbose=false)
    pbern = pnl.PressureBernoulli(rho; unsteady=false,
        correct_kuttacondition=false, backend=backend_system)
    fbern = pnl.ForceMonitor(1, 1; i_frame=1,
        normalization=pnl.RotorNormalization(rho, 2R, 1),
        correct_kuttacondition=false, verbose=false)
    kj = pnl.KuttaJoukowskiForce(rotor, 1, 1; rho, backend=backend_system,
        i_frame=1, normalization=pnl.RotorNormalization(rho, 2R, 1), verbose=false)
    (plap, flap, pbern, fbern, kj), (flap, fbern, kj)
end

r = zeros(rotor.ncells)
x_ref = Ref{Vector{Float64}}()
ct_ref = Ref{NTuple{3,Float64}}()

# Reference persistence: lets a memory-tight rung run heavy configs in their
# own process without re-running backslash_ref (16 GB local: two live dense
# Gs at R3 get the process killed). backslash_ref runs SAVE; a CONFIGS subset
# without backslash_ref LOADS.
xref_bin = joinpath(results_dir, "agreement_xref_$rung.bin")
ctref_csv = joinpath(results_dir, "agreement_ctref_$rung.csv")
if REFERENCE != "" && !any(c -> c[1] == REFERENCE, configs)
    isfile(xref_bin) || error("REFERENCE=$REFERENCE was nominated but that " *
        "config is not in CONFIGS and no saved reference exists at $xref_bin")
    x_ref[] = reinterpret(Float64, read(xref_bin)) |> collect
    length(x_ref[]) == rotor.ncells || error("saved x_ref length mismatch")
    ct_ref[] = Tuple(parse.(Float64, split(readline(ctref_csv), ",")))
    println("Loaded saved reference: ||x_ref|| = $(norm(x_ref[])), CT_ref = $(ct_ref[])")
end

for (name, make_solver) in configs
    println("--- $name ---")
    rotor.strength .= 0                      # cold start (FGS seeding gotcha)
    rotor.potential .= 0
    # restore the solve-pass core_size: the previous config's steady!
    # post-processing leaves the TARGETS offset (1e-3) active, and ctor-time
    # trees (ILUPreconditioner, FastGaussSeidel) would otherwise see
    # radius-inflated panels (ILU direct-list pattern guard trips at R1)
    rotor.core_size = rotor.core_size_panel
    solver = make_solver()
    monitors, (flap, fbern, kj) = make_monitors()
    t_steady = try
        @elapsed pnl.steady!((rotor,), frames, Vinf;
            body_solvers=(solver isa pnl.BackslashCoupled ? solver : (solver,)),
            backend=backend_system,
            backend_solve=pnl.DirectBackend(),
            monitors,
            path=nothing, verbose=false)
    catch err
        msg = sprint(showerror, err)[1:min(end, 200)]
        println("  FAILED: $msg")
        # blanks for t_steady..dCT_kj (columns 5-18), sized so the row cannot
        # drift out of step with the header
        emit_row!([rung, msh_name, rotor.ncells, name, fill("", 14)..., rms_b,
            pnl.FMM_RADIUS_TOL[], banner.threading_mode, banner.julia_threads,
            banner.commit, banner.fm_commit, banner.filament_reg, time_string(),
            "FAILED: " * replace(msg, "\n" => " ")])
        solver = nothing; monitors = nothing; GC.gc()
        continue
    end
    x = copy(vec(rotor.strength[:, 2]))
    # BC residual — FMM by default (Ryan 2026-08-24: DirectBackend is O(N^2),
    # 1.8e11 pair evaluations at R7, i.e. prohibitive). `bc_error!` evaluates
    # the BC through the FMM under an EXPLICIT error tolerance of
    # `safety * target_rel * rms_b` with safety=0.1 — an order of magnitude
    # tighter than the 1e-6 criterion being tested — and returns
    # `error_success`, the FMM's own certification that it met that tolerance.
    # An uncertified row is therefore visible rather than silently wrong.
    # RESIDUAL_BACKEND=direct restores the exact O(N^2) evaluation; worth doing
    # at R1-R2 as a cross-check that the FMM metric matches the exact one (see
    # also rotor_hover_solver_phase1_bcerror.jl, which validated this evaluator).
    local bc_rel, bc_cert, rms_d
    if residual_backend === :direct
        rotor.velocity .= frozen_velocity
        rotor.core_size = rotor.core_size_panel
        rms_d, _ = pnl.true_residual!(r, rotor, x, b; backend=pnl.DirectBackend())
        bc_rel = rms_d / rms_b
        bc_cert = true
    else
        rotor.velocity .= frozen_velocity
        rotor.core_size = rotor.core_size_panel
        e_bc = bc_error!(rotor, x; rms_b, target_rel, safety=0.1,
                         max_expansion_order=20, multipole_acceptance=mac_t,
                         leaf_size=leaf_t, backend=:fmm)
        bc_rel = e_bc.rel_l2
        bc_cert = e_bc.error_success
        rms_d = bc_rel * rms_b
        bc_cert || @warn "$name: FMM BC evaluation did NOT certify at " *
            "epsilon=$(e_bc.epsilon_requested) — bc_rel_l2=$bc_rel is an " *
            "upper bound of unknown tightness"
    end
    niter = solver isa pnl.KrylovSolver || solver isa pnl.FGSSolver ? solver.niter : -1
    CT_lap  = -flap.force[axial_dimension, 1]
    CT_bern = -fbern.force[axial_dimension, 1]
    CT_kj   = -kj.force[axial_dimension, 1]
    if name == REFERENCE
        x_ref[] = x
        ct_ref[] = (CT_lap, CT_bern, CT_kj)
        write(xref_bin, x)
        write(ctref_csv, join((CT_lap, CT_bern, CT_kj), ",") * "\n")
    end
    # Persist this config's solution so the reference-free ensemble spread can
    # be computed at the end regardless of which configs ran, in which order,
    # or in how many separate processes.
    write(joinpath(results_dir, "agreement_x_$(rung)_$(name).bin"), x)
    write(joinpath(results_dir, "agreement_ct_$(rung)_$(name).csv"),
          join((CT_lap, CT_bern, CT_kj), ",") * "\n")
    # Reference-relative columns are now OPTIONAL: blank unless a reference was
    # explicitly nominated with REFERENCE=<config> and has already run.
    relL2, relmax = "", ""
    dcts = ("", "", "")
    if isassigned(x_ref)
        dx = x .- x_ref[]
        relL2 = norm(dx) / norm(x_ref[])
        relmax = maximum(abs, dx) / maximum(abs, x_ref[])
        dcts = ((CT_lap, CT_bern, CT_kj) .- ct_ref[]) ./ ct_ref[] .* 100
    end
    println("  t=$(round(t_steady; digits=1))s niter=$niter " *
            "bc_rel_l2=$bc_rel (certified=$bc_cert, $residual_backend) " *
            "relL2=$relL2 relmax=$relmax CT(lap,bern,kj)=($CT_lap, $CT_bern, $CT_kj) " *
            "dCT%=$(dcts)")
    emit_row!([rung, msh_name, rotor.ncells, name, t_steady,
        niter, rms_d, bc_rel, bc_cert, residual_backend,
        relL2, relmax, CT_lap, CT_bern, CT_kj,
        dcts[1], dcts[2], dcts[3], rms_b, pnl.FMM_RADIUS_TOL[],
        banner.threading_mode, banner.julia_threads, banner.commit,
        banner.fm_commit, banner.filament_reg, time_string(),
        "steady! path; backend_system=tuned($p_t,$mac_t,$leaf_t); backend_solve=Direct; " *
        "monitors fresh per config; steady Bernoulli caveat accepted for deltas"])
    solver = nothing; monitors = nothing; GC.gc()
end

close(io)

# ---- reference-free ensemble spread ----------------------------------------
# Loads every per-config solution persisted for this rung (this process or an
# earlier one) and reports how far apart the solvers actually are, with no
# config privileged. Two views:
#   dev_from_mean  — each config vs the ensemble mean
#   max_pairwise   — the worst pair, which is the honest "do they agree?" number
# Both are relative L2. Written to agreement_spread.csv, one row per config
# plus one summary row.
spread_csv = joinpath(results_dir, "agreement_spread.csv")
xs = Pair{String,Vector{Float64}}[]
for f in sort(readdir(results_dir))
    m = match(Regex("^agreement_x_" * rung * "_(.+)\\.bin\$"), f)
    m === nothing && continue
    v = collect(reinterpret(Float64, read(joinpath(results_dir, f))))
    length(v) == rotor.ncells || (@warn "skipping $f: length mismatch"; continue)
    push!(xs, String(m.captures[1]) => v)
end
if length(xs) < 2
    println("\nensemble spread: only $(length(xs)) solution(s) for $rung — " *
            "need >= 2; run more configs (they may be separate processes)")
else
    names_ = first.(xs); vs = last.(xs)
    xbar = sum(vs) ./ length(vs)
    dev = [norm(v .- xbar) / norm(xbar) for v in vs]
    # in a function so the loop accumulators are not top-level soft-scope locals
    function _max_pairwise(names_, vs)
        best, arg = 0.0, ("", "")
        for i in eachindex(vs), j in (i+1):length(vs)
            d = norm(vs[i] .- vs[j]) / max(norm(vs[i]), norm(vs[j]))
            if d > best
                best = d
                arg = (names_[i], names_[j])
            end
        end
        return best, arg
    end
    maxpair, argpair = _max_pairwise(names_, vs)
    sh = "rung,n_panels,config,dev_from_mean_relL2,n_configs,max_pairwise_relL2," *
         "max_pairwise_configs,commit,fm_commit,date,notes"
    fresh_s = !isfile(spread_csv) || filesize(spread_csv) == 0
    ios = open(spread_csv, "a"); fresh_s && println(ios, sh)
    println("\n=== $rung ensemble spread over $(length(xs)) configs " *
            "(reference-free) ===")
    for (nm, d) in zip(names_, dev)
        println("  $(rpad(nm, 22)) dev from mean relL2 = $d")
        println(ios, join(_csv_cell.([rung, rotor.ncells, nm, d, length(xs),
            "", "", banner.commit, banner.fm_commit, time_string(),
            "reference-free; ensemble = $(join(names_, "|"))"]), ","))
    end
    println("  MAX PAIRWISE relL2 = $maxpair  ($(argpair[1]) vs $(argpair[2]))")
    println(ios, join(_csv_cell.([rung, rotor.ncells, "ALL", "", length(xs),
        maxpair, "$(argpair[1])|$(argpair[2])", banner.commit, banner.fm_commit,
        time_string(), "summary row; max over all config pairs"]), ","))
    close(ios)
    println("wrote $spread_csv")
end

println("\n$rung agreement done. Rows appended to $csv_path")
