#=##############################################################################
BRAINSTORM 021 Phase 2b — near-field influence-matrix cache A/B
(nearfield_matrix_cache_plan.md §5).

Measures, at the Phase-1 frozen per-rung apply knobs:
  1. Isolated near-field: cached nearfield_matvec! vs kernel
     nearfield_multithread!, on the SAME FmmPlan buffers (min-of-k), plus
     cache build time, bytes, and the break-even apply count
     build_time / (t_kernel - t_cached).
  2. End-to-end: cold krylov_gmres solve with cache_nearfield off/on
     (reset_cold! before every rep; min-of-k), niter, certified BC via
     bc_error! for both, and the solution shift between them.

Run (local smoke; published rows are HPC-only per ruling 5):
  RUNG=R1 EXPECT_JULIA_THREADS=4 THREADING_MODE=multi CACHE_B=1 \
    julia --project -t 4 benchmark/rotor_hover_solver_phase2b_nearfield_cache.jl
Appends nearfield_cache_ab.csv under benchmark/results/phase2/<mode>.
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

target_rel = 1e-6
safety = 0.1
include(joinpath(@__DIR__, "phase1_knobs.jl"))

p_t, mac_t, leaf_t = TUNED[rung]
bc_fmm(x) = bc_error!(rotor, x; rms_b, target_rel, safety,
                      max_expansion_order=20, multipole_acceptance=mac_t,
                      leaf_size=leaf_t, backend=:fmm)

outdir = joinpath(@__DIR__, "results", "phase2", banner.threading_mode)
# PER_RUNG_DIR=1 (HPC): per-rung subdirs — concurrent per-rung jobs clobber
# shared CSVs via non-atomic NFS appends (Phase-1 lesson); phase2.jl applies
# the same rule, so tune_cached.csv lands where its consumer reads it
get(ENV, "PER_RUNG_DIR", "0") == "1" && (outdir = joinpath(outdir, rung))
mkpath(outdir)
csv_path = joinpath(outdir, "nearfield_cache_ab.csv")
header = "rung,mesh_file,n_panels,row_kind,t_nearfield_kernel,t_nearfield_cached," *
    "nearfield_speedup,t_cache_build,cache_bytes,break_even_applies," *
    "t_solve_min,k_reps,niter,bc_rel_l2_certified,bc_certified," *
    "solution_shift_rel_l2,apply_p,apply_mac,apply_leaf,threading_mode," *
    "julia_threads,blas_threads,commit,fm_commit,date,hardware_tag,notes"
fresh = !isfile(csv_path) || filesize(csv_path) == 0
io = open(csv_path, "a")
fresh && println(io, header)

# ---- Resume support (2026-08-24) --------------------------------------------
# Jobs run on the preemptible `standby` QOS with --requeue, so a preempted job
# restarts from scratch and would otherwise APPEND duplicates of whatever it had
# already written. Here the row identity is (rung, row_kind) — row_kind is
# column 4, there is no config column. The two existing manual gates (SKIP_AB,
# TUNE_CACHED) already carve the script into exactly the right units, so resume
# just drives them from what has landed instead of from the environment.
landed = Set{String}()
if !fresh
    for line in Iterators.drop(eachline(csv_path), 1)
        isempty(strip(line)) && continue
        c = split(line, ",")
        length(c) >= 4 && String(c[1]) == rung && push!(landed, String(c[4]))
    end
end
const AB_ROWS = ("nearfield_isolated", "solve_krylov_gmres_cache_off",
                 "solve_krylov_gmres_cache_on")
ab_done = all(k -> k in landed, AB_ROWS)
tune_done = "tune_cached" in landed
if !isempty(landed)
    println("RESUME: $rung already has " * join(sort(collect(landed)), ", ") *
            " — A/B " * (ab_done ? "SKIPPED" : "will run") *
            ", cached-tune " * (tune_done ? "SKIPPED" : "will run"))
end

hardware_tag = get(ENV, "HARDWARE_TAG", "local-smoke")
k_iso = parse(Int, get(ENV, "K_ISO", "7"))
k_solve = parse(Int, get(ENV, "K_SOLVE", "3"))
n_threads = Threads.nthreads()

# Near-field cache caps, shared by every section below. FastMultipole's 4 GiB
# default is BELOW what the upper ladder rungs need (measured R4: 4.47 GiB for
# 207,344 blocks — a 12% overshoot that failed job 13242665 outright), and the
# cache build is SERIAL, so the time cap is wall-clock rather than a per-thread
# budget (measured R1: 19.5 s at 4 threads vs 20.0 s at 1). Both are raised
# here and overridable per rung at sbatch time.
max_bytes_t = parse(Int, get(ENV, "TUNE_MAX_BYTES", string(32 * 1024^3)))
max_time_t = parse(Float64, get(ENV, "TUNE_MAX_BUILD_TIME", "1800.0"))
println("nearfield cache caps: $(round(max_bytes_t / 1024^3; digits=2)) GiB / " *
        "$(max_time_t) s build")

backend = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)
n = rotor.ncells
scratch = zeros(n)
x_probe = sin.(0.7 .* (1:n)) .+ 0.1

# SKIP_AB=1 runs only the TUNE_CACHED section (sections 1-2 already recorded)
if get(ENV, "SKIP_AB", "0") != "1" && !ab_done

# ---- 1. isolated near-field A/B on one plan --------------------------------
slot = Ref{Any}(nothing)
pnl._apply_dirichlet_G!(rotor, x_probe, backend, scratch; plan_slot=slot)
plan = slot[][1]
tt, st = plan.target_tree, plan.source_tree
@assert length(plan.direct_list) > 0

t_kernel = minimum(@elapsed FastMultipole.nearfield_multithread!(
        tt.buffers, tt.branches, (rotor,), st.buffers, st.branches,
        plan.derivatives_switches, plan.direct_list,
        plan.interaction_list_method, n_threads)
    for _ in 1:k_iso)

cache = FastMultipole.build_nearfield_cache!(plan, (rotor,), (rotor,);
    max_bytes=max_bytes_t, max_build_time=max_time_t)
t_cached = minimum(@elapsed FastMultipole.nearfield_matvec!(
        tt.buffers, cache, st.buffers; n_threads)
    for _ in 1:k_iso)

# exactness spot check on the full operator apply (V0 at rung scale)
slot_nc = Ref{Any}(nothing)
y_ref = copy(pnl._apply_dirichlet_G!(rotor, x_probe, backend, scratch; plan_slot=slot_nc))
slot_c = Ref{Any}(nothing)
y_cached = copy(pnl._apply_dirichlet_G!(rotor, x_probe, backend, scratch;
    plan_slot=slot_c, cache_nearfield=true,
    nearfield_cache_max_bytes=max_bytes_t, nearfield_cache_max_build_time=max_time_t))
apply_rel_err = norm(y_cached .- y_ref) / norm(y_ref)
println("apply exactness: rel_l2(cached - kernel) = $apply_rel_err")
@assert apply_rel_err < 1e-12

speedup = t_kernel / t_cached
break_even = cache.build_time / max(t_kernel - t_cached, eps())
println("nearfield kernel $t_kernel s, cached $t_cached s, speedup $(round(speedup; digits=2))x")
println("cache build $(cache.build_time) s, $(cache.bytes) bytes, " *
        "break-even $(round(break_even; digits=1)) applies")

println(io, join([rung, msh_name, n, "nearfield_isolated",
    t_kernel, t_cached, speedup, cache.build_time, cache.bytes, break_even,
    "", k_iso, "", "", "", "", p_t, mac_t, leaf_t, banner.threading_mode,
    n_threads, banner.blas_threads, banner.commit, banner.fm_commit,
    time_string(), hardware_tag,
    "isolated nearfield min-of-k on shared FmmPlan; apply_rel_err=$apply_rel_err"], ","))

# ---- 2. end-to-end cold krylov_gmres A/B -----------------------------------
results = Dict{Bool,Any}()
for cache_nearfield in (false, true)
    solver = pnl.KrylovSolver(rotor; method=:gmres, itmax=200, atol=1e-10,
        rtol=1e-8, memory=200, backend, cache_tree=true, cache_nearfield,
        nearfield_cache_max_bytes=max_bytes_t,
        nearfield_cache_max_build_time=max_time_t)
    t_best = Inf
    local niter
    # Untimed warm-up. Without it the FIRST arm of the loop absorbs all the
    # JIT for the solve path, which biases this A/B in favour of whichever arm
    # runs second — here the CACHED one, i.e. exactly the direction that would
    # flatter the feature under test (found 2026-08-24; the previously reported
    # 5.5x cold-solve speedup is therefore an upper bound). Every other solver
    # measurement in the campaign gets a warm-up via adaptive_min_of_k.
    reset_cold!()
    pnl.solve!(rotor, solver)
    for _ in 1:k_solve
        reset_cold!()
        t = @elapsed pnl.solve!(rotor, solver)
        t_best = min(t_best, t)
        niter = solver.niter
        solver.solved || error("krylov_gmres (cache_nearfield=$cache_nearfield) did not converge")
    end
    xsol = copy(rotor.strength[:, solution_column])
    bc = bc_fmm(xsol)
    results[cache_nearfield] = (; t_best, niter, bc_rel=bc.rel_l2,
                                  certified=bc.error_success, xsol)
    println("krylov_gmres cache_nearfield=$cache_nearfield: t=$t_best s, " *
            "niter=$niter, bc_rel_l2=$(bc.rel_l2) certified=$(bc.error_success)")
end

shift = norm(results[true].xsol .- results[false].xsol) / norm(results[false].xsol)
println("solution shift (on vs off): rel_l2 = $shift")

for cache_nearfield in (false, true)
    r = results[cache_nearfield]
    println(io, join([rung, msh_name, n,
        "solve_krylov_gmres_cache_" * (cache_nearfield ? "on" : "off"),
        "", "", "", cache_nearfield ? cache.build_time : "",
        cache_nearfield ? cache.bytes : "", "",
        r.t_best, k_solve, r.niter, r.bc_rel, r.certified,
        cache_nearfield ? shift : 0.0,
        p_t, mac_t, leaf_t, banner.threading_mode, n_threads,
        banner.blas_threads, banner.commit, banner.fm_commit, time_string(),
        hardware_tag,
        "cold min-of-k (ruling 12); cache built fresh inside each solve's first apply"], ","))
end
end   # SKIP_AB
close(io)

# ---- 3. cached-economics tune_fmm smoke (TUNE_CACHED=1) --------------------
# Re-tunes the apply knobs with the near-field evaluated through per-trial
# throwaway caches (build cost excluded; leaf changes stop at the
# max_bytes/max_build_time caps). Local smoke only — published knobs come
# from the HPC tune stage.
if get(ENV, "TUNE_CACHED", "0") == "1" && !tune_done
    epsilon = 0.1 * target_rel * rms_b
    tune_macs = [parse(Float64, s) for s in
                 split(get(ENV, "TUNE_MACS", "0.5"), r"[:,]")]   # colon-separated at sbatch --export
    println("\n--- cached-economics tune_fmm (macs=$tune_macs, caps: " *
            "$(max_bytes_t) B / $(max_time_t) s build) ---")
    t_tune = @elapsed tuned, _, tinfo = FastMultipole.tune_fmm((rotor,), (rotor,);
        error_tolerance=FastMultipole.PowerAbsolutePotential(epsilon),
        scalar_potential=true, gradient=false, hessian=false,
        leaf_size_source=leaf_t, multipole_acceptances=tune_macs,
        max_expansion_order=20, shrink=true,
        tune_nearfield_cache=true,
        nearfield_cache_max_bytes=max_bytes_t,
        nearfield_cache_max_build_time=max_time_t)
    leaf_tuned = maximum(tuned.leaf_size_source)
    println("cached-economics tuned: p=$(tuned.expansion_order) " *
            "mac=$(tuned.multipole_acceptance) leaf=$leaf_tuned " *
            "(kernel-tuned frozen: p=$p_t mac=$mac_t leaf=$leaf_t); " *
            "cache_capped=$(tinfo.cache_capped); t_tune=$t_tune s")
    # cache cost at the tuned knobs
    plan_t = FastMultipole.FmmPlan((rotor,), (rotor,);
        scalar_potential=true, gradient=false, hessian=false,
        expansion_order=tuned.expansion_order,
        multipole_acceptance=tuned.multipole_acceptance,
        leaf_size_source=tuned.leaf_size_source, shrink=true)
    est_t = FastMultipole.estimate_nearfield_cache(plan_t.target_tree,
        plan_t.source_tree, plan_t.direct_list, plan_t.derivatives_switches,
        (rotor,))
    println("cache at tuned knobs: $(est_t.bytes) bytes, est build " *
            "$(round(est_t.est_build_time; digits=2)) s")
    io2 = open(csv_path, "a")
    println(io2, join([rung, msh_name, n, "tune_cached",
        "", "", "", est_t.est_build_time, est_t.bytes, "",
        t_tune, 1, "", "", "", "",
        tuned.expansion_order, tuned.multipole_acceptance, leaf_tuned,
        banner.threading_mode, n_threads, banner.blas_threads, banner.commit,
        banner.fm_commit, time_string(), hardware_tag,
        "tune_fmm tune_nearfield_cache=true; macs=$(tune_macs); " *
        "cache_capped=$(tinfo.cache_capped); kernel-tuned frozen p$p_t/mac$mac_t/leaf$leaf_t"], ","))
    close(io2)

    # machine-readable knob row for the phase2 nfcache configs (csv-first,
    # latest row per rung wins — mirrors the Phase-1 tune.csv convention)
    knobs_csv = joinpath(outdir, "tune_cached.csv")
    knobs_header = "rung,expansion_order,multipole_acceptance,leaf_size," *
        "cache_capped,cache_bytes,est_build_time,commit,fm_commit,date"
    fresh_knobs = !isfile(knobs_csv) || filesize(knobs_csv) == 0
    io3 = open(knobs_csv, "a")
    fresh_knobs && println(io3, knobs_header)
    println(io3, join([rung, tuned.expansion_order,
        tuned.multipole_acceptance, leaf_tuned, tinfo.cache_capped,
        est_t.bytes, est_t.est_build_time, banner.commit, banner.fm_commit,
        time_string()], ","))
    close(io3)
    println("WROTE $knobs_csv")
end

println("WROTE $csv_path")
