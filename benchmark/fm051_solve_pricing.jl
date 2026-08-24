# Task 051 Stage 3 (Job C): price the production p018 solve on CPU.
#
# Question (051-plan Stage 3): does the CPU solve fit the revised U-only
# headroom of ~0.78 s/step (3.3 s budget − B′ stack 0.373 + 0.124 + 2.02 s from
# Job A/B measurements)?  p018's body solver is `Backslash` (dense, factored
# once at RHPC setup — the rotor is rigid in its rotating frame), so pricing is
# the WARM per-step `solve_formulation!` wall (RHS build + backsolve), not a
# Krylov niter × matvec product; job 13310123's 54.9 s solve wall was a cold
# first call (JIT + first-touch) and is re-measured here as the cold point.
#
# Reuses the parity harness's full-mode p018 warm start verbatim by including
# benchmark/fm051_pass_parity.jl (its main() is guarded by PROGRAM_FILE).
# Seam stays :off throughout — this is a pure CPU measurement; no GPU needed.
#
#   julia --project=$HOME/fm051env -t 64 benchmark/fm051_solve_pricing.jl
#
# Output: benchmark/results/fm051_solve_pricing.csv + a PRICE/VERDICT block.

ENV["FM051_MODE"] = "full"
ENV["FM051_SEAM"] = "host"   # harness include requires host|cuda; the seam is
                             # never armed here (set_gpu_influence!(:off) below)

module FM051Harness
    Base.include(@__MODULE__, joinpath(dirname(@__FILE__), "fm051_pass_parity.jl"))
end

using Printf
import Profile
const H = FM051Harness
const pnl = H.pnl

const N_WARM = parse(Int, get(ENV, "FM051_PRICE_REPEATS", "5"))
const HEADROOM_S = parse(Float64, get(ENV, "FM051_PRICE_HEADROOM", "0.783"))
# Stage-3 fix pricing: opt in the constant source-potential matrix S on the
# body's Backslash and re-price the warm solve on the dense-gemv path.
const PRICE_S = get(ENV, "FM051_PRICE_S", "1") == "1"

function main()
    pnl.set_gpu_influence!(:off)
    cfg, systems_tuple, wakes_tuple, parity_i_step = H.build_full_cfg()
    println("\nFM051 solve pricing -- p018 full config, i_step=$(parity_i_step), " *
            "threads=$(Threads.nthreads())")
    for (k, s) in enumerate(cfg.body_solvers)
        note = hasproperty(s, :niter) ? " (niter tracked)" : ""
        println("  body solver $k: $(typeof(s))$note")
    end

    pnl._set_core_sizes!(systems_tuple, :core_size_panel)
    solve_once!() = pnl.solve_formulation!(cfg.formulation, cfg.formulation_state,
        cfg.systems, systems_tuple, wakes_tuple, cfg.body_solvers;
        backend_solve=cfg.backend_solve, backend_wake=cfg.backend_wake,
        i_step=parity_i_step)

    t_cold = @elapsed solve_once!()
    @printf("  cold solve (JIT + first-touch): %.3f s\n", t_cold)

    warm = Float64[]
    for r in 1:N_WARM
        t = @elapsed solve_once!()
        push!(warm, t)
        for (k, s) in enumerate(cfg.body_solvers)
            hasproperty(s, :niter) && @printf("    solver %d niter: %d\n", k, s.niter)
        end
        @printf("  warm solve %d/%d: %.3f s\n", r, N_WARM, t)
    end
    med = sort(warm)[cld(length(warm), 2)]
    fits = med <= HEADROOM_S

    # ---- decomposition: where does the warm solve spend its time? ----------
    # The Backslash inner solve is RHS + ldiv! (ms at 36k), so a ~7 s warm wall
    # must live upstream (dispatch layers, Kutta/TE, influence assembly, ...).
    println("\nDECOMPOSITION")
    println("  systems type:       ", typeof(cfg.systems))
    println("  systems_tuple types: ", join((string(typeof(s)) for s in systems_tuple), " | "))
    if length(systems_tuple) == length(cfg.body_solvers)
        for (k, (body, s)) in enumerate(zip(systems_tuple, cfg.body_solvers))
            t1 = @elapsed pnl.solve!(body, s; backend=cfg.backend_solve)
            t2 = @elapsed pnl.solve!(body, s; backend=cfg.backend_solve)
            @printf("  inner solve!(body%d): %.3f s (repeat %.3f s)\n", k, t1, t2)
        end
    end
    Profile.clear()
    Profile.@profile solve_once!()
    println("  flat profile of one warm solve (top frames by count):")
    Profile.print(IOContext(stdout, :displaysize => (80, 200));
        format=:flat, sortedby=:count, mincount=max(50, div(length(Profile.fetch()), 50)))

    # ---- Stage-3 fix: source-potential matrix S (dense gemv) --------------
    # Assemble S once on the Backslash (opt-in, ~ncells^2 * 8 B), verify
    # S*σ against the fmm-influence! potential on the production state, then
    # re-price the warm solve on the gemv path (the solve! seam picks S up
    # automatically once it is attached).
    warm_S = Float64[]
    med_S = NaN
    t_assemble = NaN
    equiv_abs = NaN
    equiv_rel = NaN
    t_gemv = NaN
    warm_S_tuned = Float64[]
    med_S_tuned = NaN
    gemv_sweep = Tuple{String,Float64}[]
    t_gemv_blocked = NaN
    blas_best_nt = 0
    if PRICE_S
        println("\nSTAGE-3 FIX: source-potential matrix S")
        body = systems_tuple[1]
        solver = pnl._single_body_solver(cfg.body_solvers)
        if !(solver isa pnl.Backslash)
            println("  SKIP: body solver is $(typeof(solver)), not Backslash")
        else
            t_assemble = @elapsed pnl.assemble_source_potential!(solver, body)
            @printf("  S assembly (one-time, %d^2): %.1f s, %.2f GB\n",
                body.ncells, t_assemble, sizeof(solver.S)/1e9)

            # equivalence on the production state: fmm influence! vs S*σ
            pnl.set_strengths!(body)
            sigma = copy(body.strength[:, 1])
            potential_old = copy(body.potential)
            body.potential .= 0
            S_backup = solver.S
            solver.S = nothing   # force the fmm path through the same seam
            pnl._source_influence!(body, solver, cfg.backend_solve)
            phi_fmm = copy(body.potential)
            solver.S = S_backup
            body.potential .= potential_old
            phi_gemv = similar(phi_fmm)
            t_gemv = @elapsed pnl.LA.mul!(phi_gemv, solver.S, sigma)
            t_gemv = min(t_gemv, @elapsed pnl.LA.mul!(phi_gemv, solver.S, sigma))
            equiv_abs = maximum(abs, phi_gemv .- phi_fmm)
            equiv_rel = equiv_abs / max(maximum(abs, phi_fmm), eps())
            # 13391706 postmortem: exact-0 equivalence is GENUINE — the body
            # self-influence via FastMultipoleBackend routes 100% direct (no
            # farfield accepted; probed bitwise fmm==direct locally), so the
            # gemv replaces an exactly-dense evaluation. Norms below guard
            # against a vacuous zero (σ = 0 would also print 0/0).
            @printf("  gemv-vs-fmm potential: max abs %.3e, rel %.3e (0 is genuine; see comment)\n",
                equiv_abs, equiv_rel)
            @printf("  norms: |sigma| %.3e, |phi_fmm| %.3e, |phi_gemv| %.3e%s\n",
                maximum(abs, sigma), maximum(abs, phi_fmm), maximum(abs, phi_gemv),
                maximum(abs, sigma) == 0 ? "  [VACUOUS: sigma == 0]" : "")
            @printf("  bare gemv: %.3f s\n", t_gemv)

            for r in 1:N_WARM
                t = @elapsed solve_once!()
                push!(warm_S, t)
                @printf("  warm solve with S %d/%d: %.3f s\n", r, N_WARM, t)
            end
            med_S = sort(warm_S)[cld(length(warm_S), 2)]

            # ---- gemv threading diagnostics (13391706 follow-up) ----------
            # That run's bare gemv was 0.453 s = ~24 GB/s effective, roughly
            # single-core DRAM bandwidth: suspected OpenBLAS dgemv
            # under-threading (BLAS threads are set independently of
            # --threads) and/or NUMA first-touch mismatch. Sweep BLAS thread
            # counts, try a Julia-threads row-blocked gemv, then re-time the
            # warm solves at the best BLAS setting (the solve! seam calls
            # plain mul!, so set_num_threads is all it takes).
            println("\n  GEMV DIAGNOSTICS")
            blas0 = pnl.LA.BLAS.get_num_threads()
            @printf("    BLAS threads at start: %d (Julia threads: %d)\n",
                blas0, Threads.nthreads())
            y = similar(phi_gemv)
            gb = sizeof(solver.S)/1e9
            best_t = Inf
            blas_best_nt = blas0
            for nt in (1, 8, 16, 32, 64)
                pnl.LA.BLAS.set_num_threads(nt)
                t = minimum([(@elapsed pnl.LA.mul!(y, solver.S, sigma)) for _ in 1:3])
                push!(gemv_sweep, ("blas$(nt)", t))
                @printf("    dgemv BLAS threads=%2d: %.3f s  (%.0f GB/s)\n", nt, t, gb/t)
                if t < best_t
                    best_t = t
                    blas_best_nt = nt
                end
            end
            pnl.LA.BLAS.set_num_threads(1)
            nrows = size(solver.S, 1)
            nblk = Threads.nthreads()
            blocked_gemv! = (yv, S, x) -> begin
                Threads.@threads :static for b in 1:nblk
                    lo = div((b-1)*nrows, nblk) + 1
                    hi = div(b*nrows, nblk)
                    pnl.LA.mul!(view(yv, lo:hi), view(S, lo:hi, :), x)
                end
                yv
            end
            t_gemv_blocked = minimum([(@elapsed blocked_gemv!(y, solver.S, sigma)) for _ in 1:3])
            @printf("    row-blocked gemv (%d Julia threads, BLAS=1): %.3f s  (%.0f GB/s), max dev vs dgemv %.1e\n",
                nblk, t_gemv_blocked, gb/t_gemv_blocked,
                maximum(abs, y .- phi_gemv))

            pnl.LA.BLAS.set_num_threads(blas_best_nt)
            @printf("    re-timing warm solves at BLAS threads=%d\n", blas_best_nt)
            for r in 1:N_WARM
                t = @elapsed solve_once!()
                push!(warm_S_tuned, t)
                @printf("  warm solve with S (tuned BLAS) %d/%d: %.3f s\n", r, N_WARM, t)
            end
            med_S_tuned = sort(warm_S_tuned)[cld(length(warm_S_tuned), 2)]
            pnl.LA.BLAS.set_num_threads(blas0)
        end
    end

    println("\nPRICE")
    @printf("  warm solve median (fmm path): %.3f s   (min %.3f, max %.3f over %d)\n",
        med, minimum(warm), maximum(warm), N_WARM)
    if !isempty(warm_S)
        @printf("  warm solve median (S gemv path): %.3f s   (min %.3f, max %.3f over %d)\n",
            med_S, minimum(warm_S), maximum(warm_S), length(warm_S))
        fits = med_S <= HEADROOM_S
    end
    if !isempty(warm_S_tuned)
        @printf("  warm solve median (S, tuned BLAS=%d): %.3f s   (min %.3f, max %.3f over %d)\n",
            blas_best_nt, med_S_tuned, minimum(warm_S_tuned), maximum(warm_S_tuned),
            length(warm_S_tuned))
        fits = min(med_S, med_S_tuned) <= HEADROOM_S
    end
    @printf("  U-only headroom:   %.3f s   (3.3 − [0.373 + 0.124 + 2.02])\n", HEADROOM_S)
    println("VERDICT: ", fits ?
        "CPU solve FITS the headroom -- device-resident matvec stays a 052+ option" :
        "CPU solve DOES NOT fit -- escalate per 051-plan Stage 3 (device matvec)")

    mkpath(H.RESULTS_DIR)
    path = joinpath(H.RESULTS_DIR, "fm051_solve_pricing.csv")
    open(path, "w") do io
        println(io, "key,value")
        println(io, "i_step,", parity_i_step)
        println(io, "threads,", Threads.nthreads())
        println(io, "solver,", join((string(typeof(s)) for s in cfg.body_solvers), ";"))
        println(io, "cold_s,", t_cold)
        for (r, t) in enumerate(warm)
            println(io, "warm_$(r)_s,", t)
        end
        println(io, "warm_median_s,", med)
        for (r, t) in enumerate(warm_S)
            println(io, "warm_S_$(r)_s,", t)
        end
        println(io, "warm_S_median_s,", med_S)
        for (r, t) in enumerate(warm_S_tuned)
            println(io, "warm_S_tuned_$(r)_s,", t)
        end
        println(io, "warm_S_tuned_median_s,", med_S_tuned)
        println(io, "blas_best_nt,", blas_best_nt)
        for (label, t) in gemv_sweep
            println(io, "gemv_$(label)_s,", t)
        end
        println(io, "gemv_blocked_s,", t_gemv_blocked)
        println(io, "S_assemble_s,", t_assemble)
        println(io, "S_gemv_s,", t_gemv)
        println(io, "S_equiv_max_abs,", equiv_abs)
        println(io, "S_equiv_max_rel,", equiv_rel)
        println(io, "headroom_s,", HEADROOM_S)
        println(io, "fits_headroom,", fits)
    end
    println("Wrote ", path)
    return fits
end

if abspath(PROGRAM_FILE) == (@__FILE__) || isempty(PROGRAM_FILE)
    main()          # exit 0 either way: the price is a measurement, not a gate
end
