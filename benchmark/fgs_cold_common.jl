# Shared implementation for the frozen cold-solve timing and profiling drivers.
# Definitions only: cold_initialize! owns fixture loading and filesystem effects.
using TOML, SHA, Statistics, Serialization, Profile

const COLD_TARGET = 1e-6
const COLD_SEEDS = Dict("R1" => (6, 0.3, 50, 10, 17),
                        "R2" => (8, 0.4, 100, 10, 17),
                        "R3" => (6, 0.3, 100, 5, 16))

function cold_seed(rung, kind)
    p, mac, leaf, inner, kp = COLD_SEEDS[rung]
    common = Dict{String,Any}("kind" => kind, "rung" => rung,
        "P" => kind == "fgs" ? p : kp,
        "MAC" => kind == "fgs" ? mac : 0.65,
        "leaf" => kind == "fgs" ? leaf : 6)
    if kind == "fgs"
        merge!(common, Dict("inner" => inner, "max_iterations" => 300,
            "tolerance" => 0.0, "sweep_order" => "lexicographic",
            "cache_leaf_lu" => true, "rlx" => 1.0))
    elseif kind == "krylov_ilu"
        merge!(common, Dict("memory" => 50, "itmax" => 500,
            "atol" => 1e-14, "rtol" => 1e-6, "ilu_leaf" => 10,
            "ilu_MAC" => 1.0, "pattern_entries_per_panel" => 8192,
            "equilibrate" => false, "diagonal_shift" => 0.0,
            "cache_tree" => true, "persistent_plan" => true,
            "cache_nearfield" => false))
    else
        error("Unknown CONFIGS entry: $kind")
    end
    common
end

cold_write_toml(path, data) = open(io -> TOML.print(io, data; sorted=true), path, "w")
cold_id(c) = bytes2hex(sha256(sprint(io -> TOML.print(io, c; sorted=true))))[1:16]

function cold_configs(rung, kinds, stage; file="")
    if !isempty(file)
        doc = TOML.parsefile(file)
        configs = haskey(doc, "configs") ? doc["configs"] : [doc]
        all(c -> c["rung"] == rung, configs) || error("CONFIG_FILE rung mismatch")
        return filter(c -> c["kind"] in kinds, configs)
    end
    stage == "verify" && error("STAGE=verify requires CONFIG_FILE")
    configs = Dict{String,Any}[]
    for kind in kinds
        seed = cold_seed(rung, kind)
        push!(configs, seed)
        stage == "screen" || continue
        # Bounded one-factor-at-a-time neighbors; do not conflate thread effects.
        axes = kind == "fgs" ? [
            "P" => [seed["P"]-2, seed["P"]+2],
            "MAC" => [seed["MAC"]-0.1, seed["MAC"]+0.1],
            "leaf" => [max(10, seed["leaf"]÷2), 2seed["leaf"]],
            "inner" => [max(1, seed["inner"]÷2), 2seed["inner"]],
            "sweep_order" => ["colored"]] : [
            "P" => [seed["P"]-1, seed["P"]+1],
            "MAC" => [0.6, 0.7], "leaf" => [4, 8], "memory" => [100],
            "cache_nearfield" => [true]]
        for (key, values) in axes, value in values
            c = copy(seed); c[key] = value; push!(configs, c)
        end
        if kind == "krylov_ilu"
            for leaf in (5,10,20), mac in (0.8,1.0)
                c = copy(seed); c["ilu_leaf"] = leaf; c["ilu_MAC"] = mac
                push!(configs, c)
            end
            # Uncached arm is explicitly diagnostic; it cannot win selection.
            c = copy(seed); c["cache_tree"] = false; c["persistent_plan"] = false
            c["diagnostic"] = true; push!(configs, c)
        end
    end
    unique(cold_id, configs)
end

function cold_initialize!()
    stage = get(ENV, "STAGE", "baseline")
    stage in ("baseline", "screen", "verify") || error("Invalid STAGE")
    rung = get(ENV, "RUNG", "")
    haskey(COLD_SEEDS, rung) || error("Cold investigation initially supports R1–R3")
    for key in ("OUTDIR", "BENCH_CASE_ROOT")
        haskey(ENV, key) && isabspath(ENV[key]) || error("$key must be an explicit absolute path")
    end
    get(ENV, "CACHE_B", "0") == "0" || error("Cold investigation requires CACHE_B=0")
    get(ENV, "SKIP_B", "0") == "0" || error("Cold investigation requires SKIP_B=0")
    get(ENV, "FLOWPANEL_FILAMENT_REG", "") == "linegauss" || error("Set FLOWPANEL_FILAMENT_REG=linegauss")
    kinds = split(get(ENV, "CONFIGS", "fgs:krylov_ilu"), r"[:,]")
    all(k -> k in ("fgs", "krylov_ilu"), kinds) || error("Invalid CONFIGS")
    configs = cold_configs(rung, kinds, stage; file=get(ENV, "CONFIG_FILE", ""))
    isempty(configs) && error("No matching configurations")
    out = ENV["OUTDIR"]
    # mkdir is the single-writer claim; never append to an existing generation.
    mkpath(dirname(out)); mkdir(out)
    case = ENV["BENCH_CASE_ROOT"]
    isdir(case) && !isempty(readdir(case)) && error("BENCH_CASE_ROOT must be fresh/empty")
    mkpath(case)
    t0 = time_ns()
    Base.include(@__MODULE__, joinpath(@__DIR__, "common.jl"))
    Base.include(@__MODULE__, joinpath(@__DIR__, "phase1_case.jl"))
    Base.invokelatest(cold_assert_threads)
    global cold_fixture_seconds = (time_ns()-t0)/1e9
    global cold_memory_limit = parse(Float64, get(ENV, "MEMORY_GIB", "500"))*2.0^30
    isfinite(cold_memory_limit) && cold_memory_limit > 0 || error("Invalid MEMORY_GIB")
    Base.invokelatest(cold_provenance, out, stage)
    return configs, out, stage
end

"Reject thread drift after fixture loading and at each measurement boundary."
function cold_assert_threads()
    Threads.nthreads() == banner.julia_threads || error("Julia thread count drifted")
    actual = LinearAlgebra.BLAS.get_num_threads()
    actual == banner.blas_threads || error(
        "BLAS thread count drifted from $(banner.blas_threads) to $actual; " *
        "investigate the reset before benchmarking")
    return nothing
end

function cold_provenance(out, stage)
    packages = Dict{String,Any}()
    for (name, mod) in (("FLOWPanel", pnl), ("FastMultipole", pnl.FastMultipole), ("FLOWVPM", pnl.FLOWVPM))
        path = pkgdir(mod)
        packages[name] = Dict("path" => path, "sha" => readchomp(`git -C $path rev-parse HEAD`),
            "status" => readchomp(`git -C $path status --porcelain`),
            "tags" => split(readchomp(`git -C $path tag --points-at HEAD`), '\n'))
    end
    affinity = Sys.islinux() ? read("/proc/self/status", String) : "unavailable on this OS"
    cold_write_toml(joinpath(out, "provenance.toml"), Dict(
        "schema" => 1, "stage" => stage, "packages" => packages,
        "julia_version" => string(VERSION), "julia_threads" => Threads.nthreads(),
        "blas_threads" => LinearAlgebra.BLAS.get_num_threads(),
        "blas" => sprint(show, LinearAlgebra.BLAS.get_config()),
        "cpu" => sprint(show, Sys.cpu_info()), "affinity" => affinity,
        "hostname" => gethostname(), "date" => time_string(),
        "filament_reg" => string(pnl.FILAMENT_REGULARIZATION[]),
        "fixture_sha256" => bytes2hex(sha256(read(msh_file))),
        "rhs_sha256" => bytes2hex(sha256(reinterpret(UInt8, b))),
        "fixture_seconds_including_load" => cold_fixture_seconds,
        "direct_rhs_seconds" => t_rhs_assembly, "memory_ceiling_bytes" => cold_memory_limit,
        "timing_scope" => "frozen _solve!; reset and BC diagnostics excluded; no formulation subsolves"))
end

function cold_make(c; history=false)
    if c["kind"] == "fgs"
        return pnl.FGSSolver(rotor; expansion_order=c["P"], multipole_acceptance=c["MAC"],
            leaf_size=c["leaf"], inner_iterations=c["inner"], max_iterations=c["max_iterations"],
            tolerance=c["tolerance"], rlx=c["rlx"], shrink=true, recenter=false,
            reverse_pass=false, cache_leaf_lu=c["cache_leaf_lu"],
            sweep_order=Symbol(c["sweep_order"]), verbose=false,
            project_solution=false, solution_history_length=0)
    end
    P = pnl.ILUPreconditioner(rotor; leaf_size=c["ilu_leaf"],
        multipole_acceptance=c["ilu_MAC"], max_pattern_entries=c["pattern_entries_per_panel"]*rotor.ncells,
        equilibrate=c["equilibrate"], diagonal_shift=c["diagonal_shift"])
    backend = pnl.FastMultipoleBackend(; expansion_order=c["P"],
        multipole_acceptance=c["MAC"], leaf_size=c["leaf"])
    pnl.KrylovSolver(rotor; method=:gmres, backend, preconditioner=P,
        itmax=c["itmax"], atol=c["atol"], rtol=c["rtol"], memory=c["memory"],
        warmstart=false, record_history=history, cache_tree=c["cache_tree"],
        cache_nearfield=c["cache_nearfield"], persistent_plan=c["persistent_plan"],
        nearfield_cache_max_bytes=floor(Int, cold_memory_limit),
        nearfield_cache_donor=Ref{Any}(nothing))
end

function cold_reset!(solver)
    reset_cold!()
    rotor.core_size = rotor.core_size_panel
    solver.niter = 0; solver.solved = false
    pnl.begin_step_solution!(solver)
    if solver isa pnl.KrylovSolver
        solver.have_x_prev = false; solver.x_history_nsaved = 0
        fill!(solver.workspace.x, 0); fill!(solver.rhs, 0)
        pnl.reset!(solver.history)
    else
        solver.solution_history_nsaved = 0
    end
    @assert all(iszero, view(rotor.strength, :, 2))
end

function cold_validate(x; crosscheck=false)
    rotor.velocity .= frozen_velocity
    phi = zeros(rotor.ncells)
    e = bc_error!(rotor, x; rms_b, target_rel=COLD_TARGET, safety=0.1,
        max_expansion_order=20, multipole_acceptance=0.5, leaf_size=20, phi_out=phi)
    direct_rel = NaN; delta = NaN; direct_seconds = 0.0
    accepted = e.error_success && e.rel_l2 <= COLD_TARGET
    if crosscheck || (!e.error_success && rung in ("R1", "R2"))
        direct_phi = similar(phi)
        d = bc_error!(rotor, x; rms_b, backend=:direct, phi_out=direct_phi)
        direct_rel = d.rel_l2; direct_seconds = d.t_eval
        delta = norm(phi-direct_phi)/sqrt(rotor.ncells)/rms_b
        # A direct fallback is authoritative, but never relabel FMM as certified.
        accepted = d.rel_l2 <= COLD_TARGET && (!e.error_success || delta <= 0.1COLD_TARGET)
    end
    (; e..., direct_rel, evaluator_delta=delta, direct_seconds, accepted)
end

function cold_csv(path, rows)
    isempty(rows) && return
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ','))
        for row in rows
            println(io, join((_csv_cell(getproperty(row,k)) for k in names), ','))
        end
    end
end

# Diagnostic staircase only: snapshot overhead is excluded from all timing metrics.
function cold_calibrate!(c, dir)
    if c["kind"] == "fgs"
        c["tolerance"] > 0 && return
        probe = copy(c); probe["tolerance"] = 0.01COLD_TARGET*rms_b
        reset_cold!(); solver = cold_make(probe)
        snapshots = Vector{Float64}[]; residuals = Float64[]; iterations = Int[]
        callback = (iteration, residual) -> begin
            pnl.FastMultipole.buffer_to_system_strength!((rotor,), solver.fgs.source_tree)
            push!(snapshots, copy(rotor.strength[:,2]))
            push!(residuals, residual); push!(iterations, iteration)
            nothing
        end
        cold_reset!(solver); pnl._solve!(rotor, solver; callback)
        final = copy(rotor.strength[:,2])
        if isempty(snapshots) || final != snapshots[end]
            push!(snapshots, final); push!(residuals, NaN)
            push!(iterations, isempty(iterations) ? 1 : iterations[end]+1)
        end
        rows = NamedTuple[]; crossing = nothing
        for i in eachindex(snapshots)
            e = cold_validate(snapshots[i])
            push!(rows, (; iteration=iterations[i], internal_residual=residuals[i], e...))
            cold_csv(joinpath(dir, "calibration.csv"), rows)
            if crossing === nothing && e.accepted && isfinite(residuals[i]) && residuals[i] > 0
                crossing = i
            elseif crossing !== nothing && e.accepted && 0 < residuals[i] < residuals[crossing]
                c["tolerance"] = sqrt(residuals[i]*residuals[crossing])
                return
            end
        end
        error("FGS staircase has no certified crossing with a decreasing successor; capped candidate")
    else
        rows = NamedTuple[]
        for _ in 1:5
            reset_cold!(); solver = cold_make(c); cold_reset!(solver)
            pnl._solve!(rotor, solver)
            e = cold_validate(copy(rotor.strength[:,2]))
            push!(rows, (; rtol=c["rtol"], iterations=solver.niter, solved=solver.solved, e...))
            cold_csv(joinpath(dir, "calibration.csv"), rows)
            solver.solved && e.accepted && return
            c["rtol"] *= 0.1
        end
        error("ILU-GMRES calibration failed after five tolerances")
    end
end

function cold_trial(solver; setup_seconds=0.0, setup_bytes=0, setup_gc=0.0, crosscheck=false)
    cold_assert_threads()
    cold_reset!(solver)
    timed = @timed pnl._solve!(rotor, solver)
    cold_assert_threads()
    x = copy(rotor.strength[:,2])
    # summarysize traverses the tuple once, deduplicating body/solver references.
    retained = Base.summarysize((rotor, solver))
    e = cold_validate(x; crosscheck)
    row = (; setup_seconds, solve_seconds=timed.time,
        total_seconds=setup_seconds+timed.time, allocated_bytes=setup_bytes+timed.bytes,
        gc_seconds=setup_gc+timed.gctime, retained_bytes=retained,
        process_peak_rss_bytes=Sys.maxrss(), iterations=solver.niter,
        inner_sweeps=solver isa pnl.FGSSolver ? solver.niter*solver.inner_iterations : -1,
        fmm_passes=solver isa pnl.FGSSolver ? solver.niter+Int(solver.solved) : -1,
        formulation_subsolves=0, solved=solver.solved, e...,
        eligible=solver.solved && e.accepted && retained <= cold_memory_limit)
    return row, x
end

function cold_benchmark(c, dir)
    cold_calibrate!(c, dir)
    cold_write_toml(joinpath(dir, "config.toml"), c)
    # Compile constructors, lazy plan/cache builds, and solves before sampling.
    reset_cold!(); solver = cold_make(c)
    cold_trial(solver)
    warm, reference = cold_trial(solver; crosscheck=rung in ("R1", "R2"))
    warm.eligible || error("Warmup failed accuracy/status/memory gate")
    rows = NamedTuple[]
    for mode in ("prepared", "fresh")
        k = warm.solve_seconds < 60 ? 5 : warm.solve_seconds < 600 ? 3 : 2
        # Fresh setup may dominate; use an excluded, compiled fresh trial to choose k.
        if mode == "fresh"
            solver = nothing; GC.gc(); reset_cold!()
            setup = @timed cold_make(c)
            solver = setup.value
            firstrow, _ = cold_trial(solver; setup_seconds=setup.time)
            k = firstrow.total_seconds < 60 ? 5 : firstrow.total_seconds < 600 ? 3 : 2
        end
        for trial in 1:k
            setup_seconds = 0.0; setup_bytes = 0; setup_gc = 0.0
            if mode == "fresh"
                solver = nothing; GC.gc(); reset_cold!()
                setup = @timed cold_make(c)
                solver = setup.value
                setup_seconds=setup.time; setup_bytes=setup.bytes; setup_gc=setup.gctime
                setup = nothing
            end
            row, x = cold_trial(solver; setup_seconds, setup_bytes, setup_gc)
            agreement = norm(x-reference)/max(norm(reference), eps())
            push!(rows, (; mode, trial, row..., fresh_prepared_relative_delta=agreement,
                repeatable=agreement <= 1e-8))
            cold_csv(joinpath(dir, "trials.csv"), rows)
        end
    end
    summaries = map(("prepared", "fresh")) do mode
        rr = filter(r -> r.mode == mode, rows)
        times = [r.total_seconds for r in rr]
        (; mode, minimum_seconds=minimum(times), median_seconds=median(times),
            maximum_seconds=maximum(times), spread_seconds=maximum(times)-minimum(times),
            repetitions=length(rr), eligible=all(r -> r.eligible && r.repeatable, rr))
    end
    cold_csv(joinpath(dir, "summary.csv"), summaries)
    # Independent recorder run; enabling history cannot inflate timing trials.
    cold_convergence(c, dir)
    return summaries
end

function cold_convergence(c, dir)
    reset_cold!(); solver = cold_make(c; history=true); cold_reset!(solver)
    rows = NamedTuple[]
    if c["kind"] == "fgs"
        callback = (iteration, residual) -> (push!(rows, (; iteration, residual)); nothing)
        pnl._solve!(rotor, solver; callback)
    else
        pnl._solve!(rotor, solver)
        h = solver.history
        for i in eachindex(h.iter)
            push!(rows, (; iteration=h.iter[i], residual=h.residual_internal[i]))
        end
    end
    cold_csv(joinpath(dir, "convergence.csv"), rows)
end

function cold_profile(c, dir)
    cold_calibrate!(c, dir)
    cold_write_toml(joinpath(dir, "config.toml"), c)
    reset_cold!(); solver = cold_make(c)
    row, _ = cold_trial(solver; crosscheck=rung in ("R1", "R2"))
    row.eligible || error("Profile candidate failed validation")
    cold_reset!(solver); Profile.clear()
    Profile.@profile pnl._solve!(rotor, solver)
    serialize(joinpath(dir, "cpu_profile.jls"), Profile.retrieve())
    for format in (:flat, :tree)
        open(joinpath(dir, "cpu_$(format).txt"), "w") do io
            Profile.print(io; format, C=true, groupby=[:thread, :task])
        end
    end
    cold_reset!(solver); Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate=0.01 pnl._solve!(rotor, solver)
    allocations = Profile.Allocs.fetch()
    serialize(joinpath(dir, "allocation_profile.jls"), allocations)
    open(joinpath(dir, "allocations.txt"), "w") do io
        println(io, "sample_rate=0.01; sampled_bytes=", sum(a.size for a in allocations.allocs; init=0))
        println(io, "Actual bytes/GC/RSS are from the separate unprofiled trial in components.csv.")
        for a in sort(allocations.allocs; by=a -> a.size, rev=true)[1:min(end,100)]
            println(io, a.size, '\t', a.stacktrace)
        end
    end
    # Profiling run has no performance-winner selection.
    cold_csv(joinpath(dir, "components.csv"), [row])
end

function cold_main(; profile=false)
    configs, out, stage = cold_initialize!()
    Base.invokelatest(cold_run, configs, out, stage; profile)
end

function cold_run(configs, out, stage; profile=false)
    selected = Dict{String,Any}[]; best = Dict{String,Float64}()
    failures = 0
    for config in configs
        c = copy(config)
        dir = joinpath(out, rung, c["kind"]*"_"*cold_id(c),
            "j$(Threads.nthreads())_b$(LinearAlgebra.BLAS.get_num_threads())")
        mkpath(dir)
        cold_write_toml(joinpath(dir, "requested_config.toml"), c)
        try
            println("Cold investigation: ", c)
            if profile
                cold_profile(c, dir)
            else
                summary = cold_benchmark(c, dir)
                prepared = first(summary)
                if all(s -> s.eligible, summary) && !get(c, "diagnostic", false) &&
                        prepared.minimum_seconds < get(best, c["kind"], Inf)
                    best[c["kind"]] = prepared.minimum_seconds
                    filter!(x -> x["kind"] != c["kind"], selected); push!(selected, copy(c))
                end
            end
            cold_write_toml(joinpath(dir, "status.toml"), Dict("status" => "completed"))
        catch err
            err isa InterruptException && rethrow()
            failures += 1
            cold_write_toml(joinpath(dir, "status.toml"), Dict("status" => "failed", "error" => sprint(showerror, err, catch_backtrace())))
            @error "Recorded failed candidate" dir exception=(err, catch_backtrace())
        end
        GC.gc()
    end
    cold_write_toml(joinpath(out, "selected.toml"), Dict("configs" => selected))
    !profile && isempty(selected) && error("No eligible candidates; inspect status.toml and trials.csv")
    stage == "verify" && failures > 0 && error("Verification has failed candidates")
    return nothing
end
