#=##############################################################################
BRAINSTORM 021 — rigid_motion_tree_reuse_item.md sketch step 1: FGS unsteady
staleness discriminator.

Hypothesis under test (item's CORRECTNESS FLAG): FastGaussSeidel builds trees
+ interaction lists + dense near-field matrices ONCE at construction; under
unsteady rigid rotation the far-field expansions are formed about t=0 branch
centers, so solver error should grow with accumulated rotation angle while
KrylovSolver (fresh trees per apply) stays converged.

Method: one short unsteady run on the FGS trajectory (production
solve_formulation! with the campaign FGS solver at the Phase-1 frozen knobs).
Immediately after each production FGS solve, the SAME state is re-solved with
a freshly constructed KrylovSolver (:gmres, rtol=1e-9; trees rebuilt per
apply, so its solution reflects the CURRENT geometry). The per-step relative
L2 divergence of the solved strength column (mu for Dirichlet) is recorded
against rotation angle; the FGS strengths are restored afterward so the
march continues on the pure-FGS trajectory (the reference solve is
side-effect-free: the outer Dirichlet solve! reassembles body.potential from
scratch and restores it, and _solve! only writes body.strength).

Expected if the bug is real: solver-tolerance-class agreement at step 1
(zero accumulated rotation), growth with angle after.

Run (local, ≤4 threads):
  RUNG=R1 N_STEPS=8 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    julia --project -t 4 benchmark/rigid_motion_fgs_staleness.jl
Output: BRAINSTORM/021_rotor_hover_solver_benchmarks/rigid_motion_tree_reuse/
        fgs_staleness.csv (this diagnostic lives with the item's docs, NOT
        under benchmark/results/ — HPC chains are in flight there).
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
import FastMultipole

banner = assert_and_banner()

const LADDER = Dict(
    "R1" => ("dji9443_20260813_23_73_capped_captess4.msh",    8016),
    "R2" => ("dji9443_20260813_33_105_capped_captess4.msh",  15760),
    "R3" => ("dji9443_20260813_45_145_capped_captess4.msh",  28752),
)
rung = get(ENV, "RUNG", "R1")
haskey(LADDER, rung) || error("RUNG must be one of $(sort(collect(keys(LADDER)))); got \"$rung\"")
msh_ladder, n_expected = LADDER[rung]

knobs_dir = joinpath(@__DIR__, "results", "phase1",
                     get(ENV, "KNOBS_MODE", banner.threading_mode))
outdir = normpath(joinpath(@__DIR__, "..", "BRAINSTORM",
    "021_rotor_hover_solver_benchmarks", "rigid_motion_tree_reuse"))
mkpath(outdir)

target_rel = 1e-6
include(joinpath(@__DIR__, "phase1_knobs.jl"))

# --- RHPC setup-only include (pattern of rotor_hover_solver_unsteady.jl) ----
setdefault!(k, v) = haskey(ENV, k) || (ENV[k] = v)
ENV["RHPC_SETUP_ONLY"] = "1"
ENV["RHPC_MESH_FILE"] = msh_ladder
setdefault!("RUN_NAME", "rigid_motion_fgs_staleness_$(rung)")
setdefault!("RUN_MONITORS", "false")
setdefault!("SAVE_VTK", "false")
setdefault!("NT", "36")

include(joinpath(pnl.examples_path, "rotor_hover_pressure_comparison.jl"))

rotor.ncells == n_expected ||
    error("Rung $rung: expected $n_expected panels, got $(rotor.ncells)")

n_steps = parse(Int, get(ENV, "N_STEPS", "8"))
# DIAG=1 relaxes the growth-vs-angle premise guards for short diagnostic
# runs (e.g. the FGS_TOL_SCALE plateau check)
diag_mode = get(ENV, "DIAG", "0") == "1"
frozen_diag = get(ENV, "FROZEN_DIAG", "0") == "1"
min_steps = diag_mode ? 2 : 6
n_steps >= min_steps || error("premise: need >=$min_steps steps; got $n_steps")
t_range = range(0.0, step=dt, length=n_steps + 1)

# --- production FGS solver at the Phase-1 frozen knobs (as the unsteady driver)
winner = stage3_winner()
winner === nothing && error("no fgsprecond winner for $rung")
sc = staircase_for(winner.p, winner.mac, winner.leaf, winner.inner)
i_cross = findfirst(t -> t[4] <= target_rel, sc)
i_cross === nothing && error("winner staircase never crosses 1e-6")
tol_abs = margin_tol(sc, i_cross)
# FGS_TOL_SCALE: diagnostic knob for the wake-on plateau hypothesis — if the
# post-fix ~2e-3 divergence floor is the FGS stopping tolerance, it should
# scale with this factor
tol_abs *= parse(Float64, get(ENV, "FGS_TOL_SCALE", "1.0"))
fgs_solver = pnl.FGSSolver(rotor;
    expansion_order=winner.p, multipole_acceptance=winner.mac,
    leaf_size=winner.leaf, inner_iterations=winner.inner,
    max_iterations=300, tolerance=tol_abs, rlx=1.0, shrink=true,
    recenter=false, reverse_pass=false, verbose=false)
body_solvers = (fgs_solver,)
println("FGS knobs: p=$(winner.p) mac=$(winner.mac) leaf=$(winner.leaf) " *
        "inner=$(winner.inner) tol_abs=$tol_abs")

# premise guards (campaign convention: non-vacuous discriminator)
fgs_solver isa pnl.FGSSolver || error("premise: production solver is not FGSSolver")
pnl.has_dirichlet_bc(rotor) ||
    error("premise: expected the Dirichlet campaign body (DBC=true)")
size(rotor.strength, 2) == 2 ||
    error("premise: expected 2 strength columns (sigma, mu); got $(size(rotor.strength, 2))")

# reference apply backend at the Phase-1 tuned shared Krylov apply knobs
tuned = Dict("R1" => (17, 0.5, 21), "R2" => (17, 0.5, 12), "R3" => (18, 0.5, 18))
parsed = read_rows(joinpath(knobs_dir, "tune.csv"))
if parsed !== nothing
    cols, rows = parsed
    for c in rows
        length(c) >= length(cols) || continue
        (c[cols["rung"]] == rung && c[cols["variant"]] == "tuned") || continue
        tuned[rung] = (parse(Int, c[cols["expansion_order"]]),
                       parse(Float64, c[cols["multipole_acceptance"]]),
                       parse(Int, c[cols["leaf_size"]]))
    end
end
p_t, mac_t, leaf_t = tuned[rung]
backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)

# --- discriminator instrumentation -------------------------------------------
"Wraps the production formulation; after each FGS solve, re-solves the same
state with a FRESH KrylovSolver (current-geometry reference) and records the
strength divergence, then restores the FGS strengths so the march continues
on the pure-FGS trajectory."
struct StalenessProbe{F} <: pnl.AbstractSolveFormulation
    inner::F
    i_steps::Vector{Int}
    div_mu_rel_l2::Vector{Float64}   # ||mu_fgs - mu_ref||2 / ||mu_ref||2
    div_mu_maxabs::Vector{Float64}
    div_sigma_maxabs::Vector{Float64} # null check: sigma is BC-set, expect 0
    ref_niter::Vector{Int}
    ref_solved::Vector{Bool}
    mu_ref_norm::Vector{Float64}
    dense_steps::Vector{Int}
    div_fgs_dense::Vector{Float64}
    div_ref_dense::Vector{Float64}
    residual_fgs::Vector{Float64}
    residual_ref::Vector{Float64}
    residual_dense::Vector{Float64}
end
pnl.initialize_formulation(f::StalenessProbe, args...) =
    pnl.initialize_formulation(f.inner, args...)
pnl.formulation_prewake!(f::StalenessProbe, state, systems_tuple) =
    pnl.formulation_prewake!(f.inner, state, systems_tuple)

"Diagnostic-only application of the cached FGS operator to one frozen vector."
function frozen_fgs_components(body, solver, x)
    fgs = solver.fgs
    source_buffers = fgs.source_tree.buffers
    target_buffers = fgs.target_tree.buffers
    saved_strength = copy(body.strength)
    saved_source_buffers = copy.(source_buffers)
    saved_target_buffers = copy.(target_buffers)
    saved_strengths = copy(fgs.strengths)
    saved_nonself_rhs = copy(fgs.nonself_matrices.rhs)
    saved_old = copy(fgs.old_influence_storage)
    saved_extra = copy(fgs.extra_right_hand_side)
    saved_influences = copy.(fgs.influences_per_system)
    switches = FastMultipole.DerivativesSwitch(true, false, false, (body,))
    try
        body.strength[:, 1] .= 0
        body.strength[:, 2] .= x
        FastMultipole.system_to_buffer!(source_buffers, (body,),
            fgs.source_tree.sort_index_list)
        FastMultipole.update_by_leaf!(fgs.strengths, fgs.strengths_by_leaf,
            (body,), source_buffers, fgs.source_tree)
        FastMultipole.update_by_leaf!(source_buffers, (body,), fgs.strengths,
            fgs.strengths_by_leaf, fgs.source_tree)
        all(iszero, view(source_buffers[1], 5, :)) ||
            error("premise: sigma was not zeroed in FGS mu-operator pass")

        self = zeros(eltype(x), length(x))
        for i_leaf in eachindex(fgs.self_matrices.sizes)
            mat, _ = FastMultipole.get_matrix_vector(fgs.self_matrices, i_leaf)
            LinearAlgebra.mul!(view(self, fgs.strengths_by_leaf[i_leaf]), mat,
                view(fgs.strengths, fgs.strengths_by_leaf[i_leaf]))
        end

        nonself_rhs = zeros(eltype(x), length(x))
        fgs.nonself_matrices.rhs .= 0
        fgs.old_influence_storage .= 0
        FastMultipole.update_nonself_influence!(nonself_rhs, fgs.strengths,
            fgs.nonself_matrices, fgs.old_influence_storage, fgs.source_tree,
            fgs.target_tree, fgs.strengths_by_leaf, fgs.index_map,
            fgs.direct_list, fgs.targets_by_branch)
        nonself = -nonself_rhs

        empty_direct = Vector{Tuple{Int,Int}}()
        FastMultipole.reset!(target_buffers)
        FastMultipole.fmm!((body,), fgs.target_tree, (body,), fgs.source_tree,
            fgs.source_tree.leaf_size, fgs.m2l_list, empty_direct, switches,
            fgs.interaction_list_method; expansion_order=fgs.source_tree.expansion_order,
            error_tolerance=nothing, upward_pass=true, horizontal_pass=true,
            downward_pass=true, reset_target_tree=true, reset_source_tree=true,
            tune=false, update_target_systems=false,
            multipole_acceptance=fgs.multipole_acceptance,
            extra_farfield=fgs.extra_farfield)
        far_rhs = zeros(eltype(x), length(x))
        FastMultipole.influence!(far_rhs, fgs.influences_per_system,
            target_buffers, (body,), source_buffers, fgs.source_tree, switches)
        far = -far_rhs

        # Independent analytic near-field action over exactly the FGS direct
        # pairs (self included), used to arbitrate the cached dense blocks.
        FastMultipole.reset!(target_buffers)
        FastMultipole.fmm!((body,), fgs.target_tree, (body,), fgs.source_tree,
            fgs.source_tree.leaf_size, Vector{Tuple{Int,Int}}(),
            fgs.full_direct_list, switches, fgs.interaction_list_method;
            expansion_order=fgs.source_tree.expansion_order,
            error_tolerance=nothing, upward_pass=false, horizontal_pass=false,
            downward_pass=false, reset_target_tree=true, reset_source_tree=false,
            tune=false, update_target_systems=false,
            multipole_acceptance=fgs.multipole_acceptance,
            extra_farfield=false)
        near_rhs = zeros(eltype(x), length(x))
        FastMultipole.influence!(near_rhs, fgs.influences_per_system,
            target_buffers, (body,), source_buffers, fgs.source_tree, switches)
        near_direct = -near_rhs

        sortidx = fgs.target_tree.sort_index_list[1]
        leaf_body_order = Int[]
        for i_branch in fgs.source_tree.leaf_index
            append!(leaf_body_order,
                sortidx[fgs.source_tree.branches[i_branch].bodies_index[1]])
        end
        length(leaf_body_order) == length(x) ||
            error("premise: incomplete FGS leaf-to-body map")
        to_body(v) = (u = similar(v); u[leaf_body_order] .= v; u)
        return (; self=to_body(self), nonself=to_body(nonself),
                far=to_body(far), near_direct=to_body(near_direct))
    finally
        body.strength .= saved_strength
        for (dst, src) in zip(source_buffers, saved_source_buffers); dst .= src; end
        for (dst, src) in zip(target_buffers, saved_target_buffers); dst .= src; end
        fgs.strengths .= saved_strengths
        fgs.nonself_matrices.rhs .= saved_nonself_rhs
        fgs.old_influence_storage .= saved_old
        fgs.extra_right_hand_side .= saved_extra
        for (dst, src) in zip(fgs.influences_per_system, saved_influences); dst .= src; end
    end
end

function pnl.solve_formulation!(f::StalenessProbe, state, systems,
        systems_tuple, wakes_tuple, body_solvers;
        backend_solve, backend_wake, i_step::Int=0, kwargs...)
    # production FGS solve
    out = pnl.solve_formulation!(f.inner, state, systems, systems_tuple,
        wakes_tuple, body_solvers; backend_solve, backend_wake, i_step,
        kwargs...)
    body = systems_tuple[1]
    strength_fgs = copy(body.strength)
    potential_fgs = copy(body.potential)
    velocity_fgs = copy(body.velocity)

    # fresh-per-step reference: trees rebuilt per apply => current geometry
    ref = pnl.KrylovSolver(body; method=:gmres, itmax=500, atol=1e-14,
                           rtol=1e-9, memory=100, backend=backend_apply)
    pnl.solve!(body, ref; backend=backend_solve)

    strength_ref = copy(body.strength)
    mu_ref = view(strength_ref, :, 2)
    mu_fgs = view(strength_fgs, :, 2)
    push!(f.i_steps, i_step)
    push!(f.div_mu_rel_l2,
          sqrt(sum(abs2, mu_fgs .- mu_ref)) / sqrt(sum(abs2, mu_ref)))
    push!(f.div_mu_maxabs, maximum(abs.(mu_fgs .- mu_ref)))
    push!(f.div_sigma_maxabs,
          maximum(abs.(view(strength_fgs, :, 1) .- body.strength[:, 1])))
    push!(f.ref_niter, ref.niter)
    push!(f.ref_solved, ref.solved)
    push!(f.mu_ref_norm, sqrt(sum(abs2, mu_ref)))

    # Optional one-shot arbitration on the first state carrying free wake
    # particles. This is a third solve of the identical frozen state; the full
    # production FGS strength matrix is restored below before time marching
    # resumes.
    wake_np = sum(hasproperty(w, :pfield) ? w.pfield.np : 0 for w in wakes_tuple)
    if frozen_diag && isempty(f.dense_steps) && wake_np > 0
        !isempty(body_solvers[1].fgs.m2l_list) || error("premise: empty FGS M2L list")
        !isempty(body_solvers[1].fgs.direct_list) || error("premise: empty FGS direct list")
        body.core_size == body.core_size_panel || error(
            "premise: active core_size $(body.core_size) != panel offset $(body.core_size_panel)")
        LinearAlgebra.norm(view(strength_fgs, :, 1)) > 0 ||
            error("premise: frozen sigma is zero")

        dense = pnl.Backslash(body)
        pnl.solve!(body, dense; backend=backend_solve)
        strength_dense = copy(body.strength)
        mu_dense = view(strength_dense, :, 2)
        b = copy(dense.rhs)
        LinearAlgebra.norm(b) > 0 || error("premise: frozen dense RHS is zero")
        all(isfinite, mu_dense) || error("premise: dense solution is nonfinite")

        rel(x, y) = LinearAlgebra.norm(x - y) / LinearAlgebra.norm(y)
        residual = zeros(eltype(b), length(b))
        pnl.true_residual!(residual, body, mu_fgs, b; backend=pnl.DirectBackend())
        rf = LinearAlgebra.norm(residual) / LinearAlgebra.norm(b)
        pnl.true_residual!(residual, body, mu_ref, b; backend=pnl.DirectBackend())
        rk = LinearAlgebra.norm(residual) / LinearAlgebra.norm(b)
        pnl.true_residual!(residual, body, mu_dense, b; backend=pnl.DirectBackend())
        rd = LinearAlgebra.norm(residual) / LinearAlgebra.norm(b)

        xprobe = [sin(0.013 * i) + 0.2 * cos(0.029 * i) for i in 1:body.ncells]
        components = frozen_fgs_components(body, body_solvers[1], xprobe)
        direct_action = zeros(eltype(b), length(b))
        pnl.apply_G!(direct_action, body, xprobe, pnl.DirectBackend())
        near_cached = components.self + components.nonself
        far_direct = direct_action - components.near_direct
        component_rel(a, ref) = LinearAlgebra.norm(a - ref) /
            max(LinearAlgebra.norm(ref), eps(eltype(ref)))

        # Recreate the FGS external RHS from the exact dense RHS and compare
        # after applying the solver's tree permutation.
        body.potential .= -b
        switches = FastMultipole.DerivativesSwitch(true, false, false, (body,))
        FastMultipole.target_influence_to_buffer!(body_solvers[1].fgs.target_tree.buffers,
            (body,), switches, body_solvers[1].fgs.target_tree.sort_index_list)
        rhs_sorted = zeros(eltype(b), length(b))
        FastMultipole.influence!(rhs_sorted, body_solvers[1].fgs.influences_per_system,
            body_solvers[1].fgs.target_tree.buffers, (body,),
            body_solvers[1].fgs.source_tree.buffers,
            body_solvers[1].fgs.source_tree, switches)
        fgs = body_solvers[1].fgs
        sortidx = fgs.target_tree.sort_index_list[1]
        leaf_body_order = Int[]
        for i_branch in fgs.source_tree.leaf_index
            append!(leaf_body_order,
                sortidx[fgs.source_tree.branches[i_branch].bodies_index[1]])
        end
        rhs_body = similar(rhs_sorted); rhs_body[leaf_body_order] .= rhs_sorted
        rhs_rel = component_rel(rhs_body, b)

        push!(f.dense_steps, i_step)
        push!(f.div_fgs_dense, rel(mu_fgs, mu_dense))
        push!(f.div_ref_dense, rel(mu_ref, mu_dense))
        push!(f.residual_fgs, rf)
        push!(f.residual_ref, rk)
        push!(f.residual_dense, rd)
        println("\n=== frozen wake-on dense arbitration (step $i_step, np=$wake_np) ===")
        println("relative solution error: FGS=$(f.div_fgs_dense[end]) " *
                "fresh Krylov=$(f.div_ref_dense[end]) dense=0")
        println("true relative residual: FGS=$rf fresh Krylov=$rk dense=$rd")
        println("FGS RHS relative error: $rhs_rel")
        println("FGS cached-near vs analytic-near: " *
                "$(component_rel(near_cached, components.near_direct))")
        println("FGS far vs direct-complement: " *
                "$(component_rel(components.far, far_direct))")
        println("FGS total vs direct operator: " *
                "$(component_rel(near_cached + components.far, direct_action))")
    end

    # continue the march on the pure-FGS trajectory
    body.strength .= strength_fgs
    body.potential .= potential_fgs
    body.velocity .= velocity_fgs
    return out
end

probe = StalenessProbe(formulation, Int[], Float64[], Float64[], Float64[],
                       Int[], Bool[], Float64[], Int[], Float64[], Float64[],
                       Float64[], Float64[], Float64[])

# --- time march --------------------------------------------------------------
sim_wall = @elapsed pnl.simulate!(systems, wakes, frames, maneuver!,
    Uinf, t_range;
    set_Das_eta_kinematic=set_Das_refresh ? init_Das_eta_kinematic : NaN,
    set_Das_min_kinematic_displacement,
    set_Das_kinematic_arc,
    set_Das_refresh,
    monitors,
    formulation=probe,
    body_solvers, backend, backend_wake,
    wakerow_no_hessian_to_particles,
    body_hessian_to_particles,
    body_gradient_core_size,
    body_on_wake,
    panel_wake_on_particles,
    particle_hessian_self,
    particle_relax,
    bound_strength_rlx,
    verbose=true,
    path=save_path, name=run_name,
)

# --- post guards + CSV -------------------------------------------------------
# simulate! also solves once at t=0 with i_step=0 before stepping, so
# n_steps+1 solves is the expected count
nrec = length(probe.i_steps)
nrec == n_steps + 1 || @warn "recorded $nrec solves for $n_steps steps (+1 initial)"
nrec >= min_steps + 1 || error("premise: fewer than $(min_steps + 1) recorded solves ($nrec)")
all(probe.ref_solved) || error("premise: reference Krylov solve did not " *
    "converge at steps $(findall(!, probe.ref_solved)) — reference invalid")
all(>(0), probe.mu_ref_norm) || error("premise: zero reference mu norm " *
    "(vacuous divergence)")
all(>(0), probe.ref_niter) || error("premise: reference solver reports 0 " *
    "iterations")
angles = [360.0 * RPM / 60.0 * (i * dt) for i in probe.i_steps]
diag_mode || angles[end] >= 30.0 || error("premise: final rotation angle " *
    "$(angles[end])° < 30° — not enough accumulated rotation to discriminate")
maximum(probe.div_sigma_maxabs) == 0.0 ||
    @warn "sigma columns differ between FGS and reference (max " *
          "$(maximum(probe.div_sigma_maxabs))) — BC assembly not deterministic?"

csv_path = joinpath(outdir, "fgs_staleness.csv")
fresh = !isfile(csv_path) || filesize(csv_path) == 0
open(csv_path, "a") do io
    fresh && println(io, "rung,config,step,t,angle_deg,div_mu_rel_l2," *
        "div_mu_maxabs,div_sigma_maxabs,mu_ref_norm,ref_niter,ref_rtol," *
        "fgs_knobs,apply_p,apply_mac,apply_leaf,nt,n_steps,threading_mode," *
        "julia_threads,commit,fm_commit,filament_reg,date")
    arm = get(ENV, "ARM", "fgs_vs_fresh_krylov")
    for (k, i) in enumerate(probe.i_steps)
        println(io, join(_csv_cell.([rung, arm, i,
            i * dt, angles[k], probe.div_mu_rel_l2[k],
            probe.div_mu_maxabs[k], probe.div_sigma_maxabs[k],
            probe.mu_ref_norm[k], probe.ref_niter[k], 1e-9,
            "p=$(winner.p);mac=$(winner.mac);leaf=$(winner.leaf);" *
            "inner=$(winner.inner);tol_abs=$tol_abs",
            p_t, mac_t, leaf_t, parse(Int, ENV["NT"]), n_steps,
            banner.threading_mode, banner.julia_threads, banner.commit,
            banner.fm_commit, banner.filament_reg, time_string()]), ","))
    end
end

println("\n=== FGS staleness discriminator ($rung, $(n_steps) steps, " *
        "sim wall $(round(sim_wall; digits=1)) s) ===")
println("step  angle°   div_mu_rel_l2   div_mu_maxabs   ref_niter")
for (k, i) in enumerate(probe.i_steps)
    println(rpad(i, 6), rpad(round(angles[k]; digits=1), 9),
            rpad(@sprintf("%.3e", probe.div_mu_rel_l2[k]), 16),
            rpad(@sprintf("%.3e", probe.div_mu_maxabs[k]), 16),
            probe.ref_niter[k])
end
growth = probe.div_mu_rel_l2[end] / max(probe.div_mu_rel_l2[1], eps())
println("divergence growth factor step $(probe.i_steps[1])→$(probe.i_steps[end]): " *
        "$(round(growth; digits=1))x")
println("Rows appended to $csv_path")
