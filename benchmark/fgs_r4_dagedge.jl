#!/usr/bin/env julia
# 021 :dagedge HPC benchmark driver (fgs_dagedge_design_20260924.md;
# reset prompt fgs_dagedge_benchmark_reset_prompt_20260924.md). One invocation
# = ONE fresh process = one (j, placement, idle) point running a BATCH of
# executor arms ("cold" = zero-initial-guess solves, NOT fresh process per arm
# — Ryan 2026-09-23). The Slurm launcher (run_r4_fgs_dagedge.slurm.sh)
# sequences (j, block) processes.
#
# Env (beyond the standard cold-harness set consumed by fgs_cold_common.jl):
#   DAGTEAM_CONFIG   — absolute path to this j's calibrated dagteam_selected.toml
#   DAGEDGE_ARMS     — comma list of arm specs, each "dagteam" or
#                      "dagedge:<theta>" (theta in bytes, 0 = split all big
#                      edges), run in list order within this process
#   DAGEDGE_DIAG     — 1 enables the FastMultipole diagnostics dict on every
#                      solve (separate matched runs; NEVER pool instrumented
#                      and uninstrumented rows — profile runs are not
#                      performance trials)
#   DAGEDGE_SOLVES   — timed prepared solves per arm after compile+warmup
#                      (default 5)
#   DAGEDGE_BLOCK    — fresh-process block id (int; recorded per row)
#   DAGEDGE_LABEL    — variant label (primary | ladder | theta ...)
#   DAGTEAM_WORKERS  — sweep-team size cap (0 = all threads, the default)
#   DAGTEAM_IDLE     — spin | backoff (dagteam idle policy; dagedge
#                      dependency-wait pause policy)
include(joinpath(@__DIR__, "fgs_cold_common.jl"))
using Statistics

const DAGEDGE_DIAG_KEYS = (:total_ns, :initialization_ns, :fmm_ns,
    :influence_mapping_ns, :residual_ns, :leaf_solve_ns, :nonself_product_ns,
    :scatter_ns, :remaining_iteration_ns, :final_update_ns,
    :dagteam_spawn_ns, :dagteam_join_ns, :dagteam_wait_ns, :dagteam_reduce_ns,
    :dagteam_busy_lower_ns, :dagteam_busy_back_ns, :dagteam_lockmgmt_ns,
    :dagteam_idle_ns, :dagteam_busy_max_ns, :dagteam_busy_min_ns,
    :dagteam_empty_pops, :dagteam_n_lower, :dagteam_n_back, :dagteam_team_size,
    :outer_count, :sweep_count, :leaf_visit_count)

# -1 sentinels keep one row schema across instrumented and uninstrumented runs
dagedge_diag_columns(diagd) = NamedTuple{map(k -> Symbol(:diag_, k), DAGEDGE_DIAG_KEYS)}(
    diagd === nothing ? map(_ -> -1, DAGEDGE_DIAG_KEYS) :
                        map(k -> Int(diagd[k]), DAGEDGE_DIAG_KEYS))

# Free schedule metadata captured at construction (reset prompt requirement):
# plan.ntasks, big/small/rootfin/back static-task counts, and the build-time
# list-scheduling simulation's makespan and byte-weighted critical path.
function dagedge_plan_columns(solver)
    plan = solver.fgs.dagteam
    if plan isa pnl.FastMultipole.DagEdgePlan
        kinds = [t.kind for w in plan.lists for t in w]
        (; plan_ntasks=plan.ntasks,
           plan_n_big=count(==(0x01), kinds), plan_n_small=count(==(0x02), kinds),
           plan_n_rootfin=count(==(0x03), kinds), plan_n_back=count(==(0x04), kinds),
           plan_workers=length(plan.lists),
           sim_makespan_bytes=plan.sim_makespan, sim_edge_L_bytes=plan.sim_edge_L)
    else
        (; plan_ntasks=plan === nothing ? -1 : plan.ntasks,
           plan_n_big=-1, plan_n_small=-1, plan_n_rootfin=-1, plan_n_back=-1,
           plan_workers=-1, sim_makespan_bytes=-1.0, sim_edge_L_bytes=-1.0)
    end
end

# Page-placement verification: per-NUMA-node resident pages after warmup.
function dagedge_numa_snapshot(out)
    Sys.islinux() || return
    pages = Dict{Int,Int}()
    for line in eachline("/proc/self/numa_maps")
        for token in split(line)
            m = match(r"^N(\d+)=(\d+)$", token)
            m === nothing && continue
            node = parse(Int, m.captures[1])
            pages[node] = get(pages, node, 0) + parse(Int, m.captures[2])
        end
    end
    rows = [(; node, pages=pages[node]) for node in sort(collect(keys(pages)))]
    cold_csv(joinpath(out, "numa_pages.csv"), rows)
    return nothing
end

function dagedge_parse_arms(spec)
    arms = NamedTuple[]
    for (k, entry) in enumerate(split(spec, ','))
        parts = split(strip(entry), ':')
        order = String(parts[1])
        if order == "dagteam"
            length(parts) == 1 || error("dagteam arm takes no theta")
            push!(arms, (; arm=k, sweep_order="dagteam", theta=-1))
        elseif order == "dagedge"
            length(parts) == 2 || error("dagedge arm spec must be dagedge:<theta>")
            theta = parse(Int, parts[2])
            theta >= 0 || error("dagedge theta must be >= 0")
            push!(arms, (; arm=k, sweep_order="dagedge", theta))
        else
            error("Arm spec must be dagteam or dagedge:<theta> (got $entry)")
        end
    end
    isempty(arms) && error("DAGEDGE_ARMS is empty")
    return arms
end

function dagedge_main()
    configs, out, stage = cold_initialize!()
    Base.invokelatest(run_dagedge, out)
end

function run_dagedge(out)
    file = get(ENV, "DAGTEAM_CONFIG", "")
    isabspath(file) && isfile(file) || error("DAGTEAM_CONFIG must be an existing absolute file")
    base = cold_check_config(TOML.parsefile(file); selected=true)
    base["kind"] == "fgs" && base["sweep_order"] == "dagteam" || error("Expected a calibrated dagteam configuration")
    arms = dagedge_parse_arms(get(ENV, "DAGEDGE_ARMS", ""))
    diag = get(ENV, "DAGEDGE_DIAG", "0") == "1"
    nsolves = parse(Int, get(ENV, "DAGEDGE_SOLVES", "5"))
    nsolves >= 1 || error("DAGEDGE_SOLVES must be positive")
    block = parse(Int, get(ENV, "DAGEDGE_BLOCK", "0"))
    label = get(ENV, "DAGEDGE_LABEL", "primary")
    workers = parse(Int, get(ENV, "DAGTEAM_WORKERS", "0"))
    workers >= 0 || error("DAGTEAM_WORKERS must be >= 0")
    idle = get(ENV, "DAGTEAM_IDLE", "spin")
    idle in ("spin", "backoff") || error("DAGTEAM_IDLE must be spin or backoff")

    # fixed workload (plan C): 27 outer x 3 inner, stopping disabled; solves
    # are not convergence-gated but must run exactly 27 outer iterations,
    # agree to the repeat gate, and pass the independent evaluator.
    base["inner"] == 3 || error("Fixed-work arm expects the champion inner=3 (got $(base["inner"]))")
    base["max_iterations"] = 27
    base["tolerance"] = 0.0
    workers != 0 && (base["dagteam_workers"] = workers)
    idle != "spin" && (base["dagteam_idle"] = idle)
    cp(file, joinpath(out, "input_dagteam_config.toml"))

    diagd = diag ? Dict{Symbol,UInt64}() : nothing
    expected_iterations = 27

    gate = function (row, agreement, arm, phase)
        ok_accuracy = row.accepted && row.authoritative_evaluator == "certified_fmm"
        ok_work = row.iterations == expected_iterations
        cold_require(ok_accuracy && ok_work && isfinite(agreement) && agreement <= 1e-8,
            "$label arm$arm $phase failed the accuracy/work/repeat gate " *
            "(accepted=$(row.accepted), evaluator=$(row.authoritative_evaluator), " *
            "iterations=$(row.iterations), agreement=$agreement)")
        return nothing
    end

    rows = NamedTuple[]
    summaries = Dict{String,Any}[]
    cross_reference = nothing   # first arm's solution: cross-executor agreement
    for spec in arms
        c = copy(base)
        c["sweep_order"] = spec.sweep_order
        if spec.sweep_order == "dagedge"
            c["dagedge_theta"] = spec.theta
        end
        cold_write_toml(joinpath(out, "config_arm$(spec.arm).toml"), c)

        reset_cold!()
        solver = cold_make(c)
        plan_cols = dagedge_plan_columns(solver)
        stamp = (row, phase, solve, agreement, cross) -> (; block, label,
            spec.arm, sweep_order=spec.sweep_order, theta=spec.theta,
            dagteam_workers=workers, dagteam_idle=idle, instrumented=diag,
            phase, solve, relative_solution_delta=agreement,
            cross_executor_delta=cross, row..., plan_cols...,
            dagedge_diag_columns(diagd)...)

        compile_row, reference = cold_trial(solver; diagnostics=diagd)
        gate(compile_row, 0.0, spec.arm, "compile")
        cross_reference === nothing && (cross_reference = reference)
        cross = norm(reference - cross_reference) / max(norm(cross_reference), eps())
        # dagteam and dagedge iterates are mathematically equal, not bitwise:
        # at f32full the reduce-order rounding paths differ (R2 smoke: 3.4e-7),
        # so this is a divergence tripwire, not the accuracy gate — the binding
        # gate is the independent evaluator on every solve.
        cross <= 1e-5 || error("$label arm$(spec.arm) cross-executor disagreement ($cross)")
        push!(rows, stamp(compile_row, "compile", 0, 0.0, cross))
        warm, warm_x = cold_trial(solver; diagnostics=diagd)
        warm_agreement = norm(warm_x - reference) / max(norm(reference), eps())
        gate(warm, warm_agreement, spec.arm, "warmup")
        push!(rows, stamp(warm, "warmup", 0, warm_agreement, cross))
        cold_csv(joinpath(out, "dagedge_solves.csv"), rows)
        spec.arm == 1 && dagedge_numa_snapshot(out)

        for solve in 1:nsolves
            row, x = cold_trial(solver; diagnostics=diagd)
            agreement = norm(x - reference) / max(norm(reference), eps())
            gate(row, agreement, spec.arm, "trial $solve")
            push!(rows, stamp(row, "trial", solve, agreement, cross))
            cold_csv(joinpath(out, "dagedge_solves.csv"), rows)
        end

        # residual history: one extra recorded solve outside the timing set;
        # both executors are deterministic so it replays the trials
        history = NamedTuple[]
        cold_reset!(solver)
        pnl._solve!(rotor, solver;
            callback=(iteration, residual) -> (push!(history, (; iteration, residual)); nothing))
        cold_assert_threads()
        all(r -> isfinite(r.residual), history) || error("Nonfinite residual in recorded history")
        cold_csv(joinpath(out, "residual_history_arm$(spec.arm).csv"), history)
        x_hist = copy(rotor.strength[:, 2])
        hist_agreement = norm(x_hist - reference) / max(norm(reference), eps())
        hist_agreement <= 1e-8 || error("History solve disagrees with trials ($hist_agreement)")

        trials = filter(r -> r.phase == "trial" && r.arm == spec.arm, rows)
        times = [r.solve_seconds for r in trials]
        push!(summaries, Dict{String,Any}(
            "arm" => spec.arm, "sweep_order" => spec.sweep_order,
            "theta" => spec.theta, "label" => label, "block" => block,
            "instrumented" => diag, "dagteam_workers" => workers,
            "dagteam_idle" => idle,
            "julia_threads" => Threads.nthreads(),
            "solves" => length(trials),
            "median_solve_seconds" => median(times),
            "min_solve_seconds" => minimum(times),
            "max_solve_seconds" => maximum(times),
            "cross_executor_delta" => cross,
            "plan_ntasks" => plan_cols.plan_ntasks,
            "plan_n_big" => plan_cols.plan_n_big,
            "plan_n_small" => plan_cols.plan_n_small,
            "plan_n_rootfin" => plan_cols.plan_n_rootfin,
            "plan_n_back" => plan_cols.plan_n_back,
            "plan_workers" => plan_cols.plan_workers,
            "sim_makespan_bytes" => plan_cols.sim_makespan_bytes,
            "sim_edge_L_bytes" => plan_cols.sim_edge_L_bytes,
            "history_iterations" => length(history),
            "history_agreement" => hist_agreement))
        cold_write_toml(joinpath(out, "dagedge_summary.toml"), Dict("arms" => summaries))
        solver = nothing
        GC.gc()
    end
    cold_write_toml(joinpath(out, "status.toml"), Dict("status" => "completed"))
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && dagedge_main()
