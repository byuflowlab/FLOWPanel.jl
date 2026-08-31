#=##############################################################################
BRAINSTORM 021 Phase 2a — bottleneck attribution: statistical profile +
allocation profile of one isolated cold solve per roster config, at the
Phase-1 frozen knobs. Feeds Stage 2b's "whatever profiling reveals" lever
list; ranks the dominant per-step cost per solver.

Outputs per config under benchmark/results/phase2/profile/<mode>/<rung>/:
  profile_<config>_flat.txt   Profile.print(format=:flat, sortedby=:count),
                              C frames off, mincount-filtered
  profile_<config>_tree.txt   Profile.print(format=:tree, maxdepth capped)
  allocs_<config>.txt         Profile.Allocs top sites by total bytes
  profile_summary.csv         one row per config: samples, top frame,
                              alloc_total_bytes, top alloc site

Solver constructions mirror rotor_hover_solver_phase2.jl (keep in sync — the
componentized-timing contract there prevents sharing one constructor).

Run (dev target = R2/R3 local 4T; the publishable attribution reruns
mid-ladder on the pinned HPC node, both modes):
  RUNG=R2 CONFIGS=fgs:krylov_ilu EXPECT_JULIA_THREADS=4 THREADING_MODE=multi \
    CACHE_B=1 julia --project -t 4 benchmark/rotor_hover_solver_phase2_profile.jl
=###############################################################################

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "phase1_case.jl"))

using Profile
import Profile.Allocs

target_rel = 1e-6
include(joinpath(@__DIR__, "phase1_knobs.jl"))

p_t, mac_t, leaf_t = TUNED[rung]

want = split(get(ENV, "CONFIGS",
    "backslash_ldiv,krylov_gmres,krylov_jacobi,krylov_ilu,fgs,fgmres_fgs"),
    r"[:,]")

outdir = joinpath(@__DIR__, "results", "phase2", "profile",
                  banner.threading_mode, rung)
mkpath(outdir)
write(joinpath(outdir, "banner.txt"), banner.text * "\n")

summary_path = joinpath(outdir, "profile_summary.csv")
fresh = !isfile(summary_path) || filesize(summary_path) == 0
sio = open(summary_path, "a")
fresh && println(sio, "rung,config,n_samples,top_flat_frame,alloc_total_bytes," *
    "n_alloc_records,top_alloc_site,threading_mode,julia_threads,commit," *
    "fm_commit,date")

# ---- constructions (mirror rotor_hover_solver_phase2.jl; keep in sync) -----
backend_apply = pnl.FastMultipoleBackend(; expansion_order=p_t,
    multipole_acceptance=mac_t, leaf_size=leaf_t)
krylov_kw = (; itmax=500, atol=1e-14, rtol=target_rel, memory=50,
             backend=backend_apply)

function make_solver(config)
    config == "backslash_ldiv" && return pnl.Backslash(rotor)
    config == "krylov_gmres" && return pnl.KrylovSolver(rotor; method=:gmres,
        krylov_kw...)
    config == "krylov_jacobi" && return pnl.KrylovSolver(rotor; method=:gmres,
        krylov_kw..., preconditioner=pnl.FastMultipole.JacobiPreconditioner(
            (rotor,); cell_size=R/4))
    config == "krylov_ilu" && return pnl.KrylovSolver(rotor; method=:gmres,
        krylov_kw..., preconditioner=pnl.ILUPreconditioner(rotor;
            leaf_size=10, multipole_acceptance=1.0,
            max_pattern_entries=8192 * rotor.ncells))
    if config in ("fgs", "fgmres_fgs")
        winner = stage3_winner()
        winner === nothing && error("no fgsprecond winner for $rung")
        sc = staircase_for(winner.p, winner.mac, winner.leaf, winner.inner)
        i_cross = findfirst(t -> t[4] <= target_rel, sc)
        i_cross === nothing && error("winner staircase never crosses 1e-6")
        tol_abs = margin_tol(sc, i_cross)
        config == "fgs" && return pnl.FGSSolver(rotor;
            expansion_order=winner.p, multipole_acceptance=winner.mac,
            leaf_size=winner.leaf, inner_iterations=winner.inner,
            max_iterations=300, tolerance=tol_abs, rlx=1.0, shrink=true,
            recenter=false, reverse_pass=false, verbose=false)
        P = pnl.FGSPreconditioner(rotor; sweeps=winner.sweeps,
            inner_iterations=winner.inner, rlx=1.0, expansion_order=winner.p,
            multipole_acceptance=winner.mac, leaf_size=winner.leaf,
            shrink=true, recenter=false)
        return pnl.KrylovSolver(rotor; method=:fgmres, krylov_kw...,
            preconditioner=P)
    end
    error("Unknown config $(repr(config))")
end

"First profile line that resolves to FLOWPanel/FastMultipole code (skip
Base/task frames) — a crude but stable 'top frame' for the summary row."
function top_own_frame(lines)
    for l in lines
        (occursin("FLOWPanel", l) || occursin("FastMultipole", l)) && return strip(l)
    end
    return length(lines) >= 2 ? strip(lines[2]) : ""
end

for config in
        ("backslash_ldiv", "krylov_gmres", "krylov_jacobi", "krylov_ilu",
         "fgs", "fgmres_fgs")
    config in want || continue
    println("--- profiling $config ---")
    solver = make_solver(config)
    reset_cold!(); pnl._solve!(rotor, solver)      # compile + warm

    # statistical profile of one cold solve
    Profile.clear()
    reset_cold!()
    Profile.@profile pnl._solve!(rotor, solver)
    n_samples = length(Profile.fetch())
    flat_path = joinpath(outdir, "profile_$(config)_flat.txt")
    open(flat_path, "w") do io
        Profile.print(io; format=:flat, sortedby=:count, C=false, mincount=5)
    end
    open(joinpath(outdir, "profile_$(config)_tree.txt"), "w") do io
        Profile.print(io; format=:tree, C=false, mincount=10, maxdepth=40)
    end
    top_frame = top_own_frame(readlines(flat_path))

    # allocation profile of one cold solve (sampled — full tracking is too
    # slow at N^2 scales; rate covers the big sites)
    reset_cold!()
    Allocs.clear()
    Allocs.@profile sample_rate=0.01 pnl._solve!(rotor, solver)
    aresults = Allocs.fetch()
    by_site = Dict{String,Int}()
    total_bytes = 0
    for a in aresults.allocs
        total_bytes += a.size
        site = isempty(a.stacktrace) ? "unknown" : string(a.stacktrace[1])
        by_site[site] = get(by_site, site, 0) + a.size
    end
    sites = sort(collect(by_site); by=x -> -x[2])
    open(joinpath(outdir, "allocs_$(config).txt"), "w") do io
        println(io, "# sampled allocation profile (sample_rate=0.01), one cold solve")
        println(io, "# total sampled bytes = $total_bytes over $(length(aresults.allocs)) records")
        for (site, bytes) in sites[1:min(end, 30)]
            println(io, "$bytes\t$site")
        end
    end
    top_site = isempty(sites) ? "" : sites[1][1]

    println(sio, join(_csv_cell.([rung, config, n_samples, top_frame,
        total_bytes, length(aresults.allocs), top_site,
        banner.threading_mode, banner.julia_threads, banner.commit,
        banner.fm_commit, time_string()]), ","))
    flush(sio)
    println("  $n_samples samples; top: $top_frame")
    solver = nothing; GC.gc()
end

close(sio)
println("\n$rung profiles written to $outdir")
