#!/usr/bin/env julia
# 021 FGS scalability diagnostic — Stage 1 driver (plan C,
# fgs_scalability_diagnostic_plan_20260921c.md). One invocation = ONE fresh
# process = one (arm, thread count, block, variant) configuration; the Slurm
# launcher (run_r4_fgs_stage1.slurm.sh) sequences blocks and A/Bs.
#
# Env (beyond the standard cold-harness set consumed by fgs_cold_common.jl):
#   DAGTEAM_CONFIG   — absolute path to this j's calibrated dagteam_selected.toml
#   STAGE1_ARM       — fixed    : fixed work (max_iterations=27, inner=3,
#                                 tolerance=0.0); solves are NOT convergence-
#                                 gated (solved=false by construction) but must
#                                 run exactly 27 outer iterations, agree to the
#                                 1e-8 repeat gate, and record BC certification
#                      accepted : the calibrated production configuration,
#                                 gated exactly like ab_trials
#   STAGE1_DIAG      — 1 enables the FastMultipole coarse-phase diagnostics
#                      dict on every solve (separate matched runs; never mix
#                      instrumented and uninstrumented rows in one analysis)
#   STAGE1_SOLVES    — timed prepared solves after compile+warmup (default 5)
#   STAGE1_BLOCK     — fresh-process block id (int; recorded per row)
#   STAGE1_LABEL     — variant label (ladder | placement-<name> | cap<N> ...)
#   DAGTEAM_WORKERS  — sweep-team size cap (0 = all threads, the default)
include(joinpath(@__DIR__, "fgs_cold_common.jl"))
using Statistics

const STAGE1_DIAG_KEYS = (:total_ns, :initialization_ns, :fmm_ns,
    :influence_mapping_ns, :residual_ns, :leaf_solve_ns, :nonself_product_ns,
    :scatter_ns, :remaining_iteration_ns, :final_update_ns,
    :dagteam_spawn_ns, :dagteam_join_ns, :dagteam_wait_ns, :dagteam_reduce_ns,
    # Stage-2 per-worker drain-loop aggregates (zero on pre-Stage-2 pins)
    :dagteam_busy_lower_ns, :dagteam_busy_back_ns, :dagteam_lockmgmt_ns,
    :dagteam_idle_ns, :dagteam_busy_max_ns, :dagteam_busy_min_ns,
    :dagteam_empty_pops, :dagteam_n_lower, :dagteam_n_back, :dagteam_team_size,
    :outer_count, :sweep_count, :leaf_visit_count)

# -1 sentinels keep one row schema across instrumented and uninstrumented runs
stage1_diag_columns(diagd) = NamedTuple{map(k -> Symbol(:diag_, k), STAGE1_DIAG_KEYS)}(
    diagd === nothing ? map(_ -> -1, STAGE1_DIAG_KEYS) :
                        map(k -> Int(diagd[k]), STAGE1_DIAG_KEYS))

# Page-placement verification (measurement contract): per-NUMA-node resident
# pages of this process after warmup, from /proc/self/numa_maps.
function stage1_numa_snapshot(out)
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

function stage1_main()
    configs, out, stage = cold_initialize!()
    Base.invokelatest(run_stage1, out)
end

function run_stage1(out)
    file = get(ENV, "DAGTEAM_CONFIG", "")
    isabspath(file) && isfile(file) || error("DAGTEAM_CONFIG must be an existing absolute file")
    c = cold_check_config(TOML.parsefile(file); selected=true)
    c["kind"] == "fgs" && c["sweep_order"] == "dagteam" || error("Expected a calibrated dagteam configuration")
    arm = get(ENV, "STAGE1_ARM", "fixed")
    arm in ("fixed", "accepted") || error("STAGE1_ARM must be fixed or accepted")
    diag = get(ENV, "STAGE1_DIAG", "0") == "1"
    nsolves = parse(Int, get(ENV, "STAGE1_SOLVES", "5"))
    nsolves >= 1 || error("STAGE1_SOLVES must be positive")
    block = parse(Int, get(ENV, "STAGE1_BLOCK", "0"))
    label = get(ENV, "STAGE1_LABEL", "ladder")
    workers = parse(Int, get(ENV, "DAGTEAM_WORKERS", "0"))
    workers >= 0 || error("DAGTEAM_WORKERS must be >= 0")
    idle = get(ENV, "DAGTEAM_IDLE", "spin")
    idle in ("spin", "backoff") || error("DAGTEAM_IDLE must be spin or backoff")

    if arm == "fixed"
        # plan-C fixed workload: 27 outer x 3 inner, stopping disabled
        c["inner"] == 3 || error("Fixed-work arm expects the champion inner=3 (got $(c["inner"]))")
        c["max_iterations"] = 27
        c["tolerance"] = 0.0
    end
    workers != 0 && (c["dagteam_workers"] = workers)
    idle != "spin" && (c["dagteam_idle"] = idle)
    cp(file, joinpath(out, "input_dagteam_config.toml"))
    cold_write_toml(joinpath(out, "config.toml"), c)

    diagd = diag ? Dict{Symbol,UInt64}() : nothing
    expected_iterations = arm == "fixed" ? 27 : -1

    gate = function (row, agreement, phase)
        ok_accuracy = row.accepted && row.authoritative_evaluator == "certified_fmm"
        ok_work = arm == "fixed" ? row.iterations == expected_iterations :
                                   row.solved && row.eligible
        cold_require(ok_accuracy && ok_work && isfinite(agreement) && agreement <= 1e-8,
            "$label $arm $phase failed the accuracy/work/repeat gate " *
            "(accepted=$(row.accepted), evaluator=$(row.authoritative_evaluator), " *
            "iterations=$(row.iterations), agreement=$agreement)")
        return nothing
    end

    reset_cold!()
    solver = cold_make(c)
    rows = NamedTuple[]
    stamp = (row, phase, solve, agreement) -> (; block, label, arm,
        dagteam_workers=workers, dagteam_idle=idle, instrumented=diag, phase, solve,
        relative_solution_delta=agreement, row..., stage1_diag_columns(diagd)...)

    compile_row, reference = cold_trial(solver; diagnostics=diagd)
    gate(compile_row, 0.0, "compile")
    push!(rows, stamp(compile_row, "compile", 0, 0.0))
    warm, warm_x = cold_trial(solver; diagnostics=diagd)
    warm_agreement = norm(warm_x - reference) / max(norm(reference), eps())
    gate(warm, warm_agreement, "warmup")
    push!(rows, stamp(warm, "warmup", 0, warm_agreement))
    cold_csv(joinpath(out, "stage1_solves.csv"), rows)
    stage1_numa_snapshot(out)

    for solve in 1:nsolves
        row, x = cold_trial(solver; diagnostics=diagd)
        agreement = norm(x - reference) / max(norm(reference), eps())
        gate(row, agreement, "trial $solve")
        push!(rows, stamp(row, "trial", solve, agreement))
        cold_csv(joinpath(out, "stage1_solves.csv"), rows)
    end

    # residual history (measurement contract): one extra recorded solve,
    # outside the timing set; dagteam is deterministic so it replays the trials
    history = NamedTuple[]
    cold_reset!(solver)
    pnl._solve!(rotor, solver;
        callback=(iteration, residual) -> (push!(history, (; iteration, residual)); nothing))
    cold_assert_threads()
    all(r -> isfinite(r.residual), history) || error("Nonfinite residual in recorded history")
    cold_csv(joinpath(out, "residual_history.csv"), history)
    x_hist = copy(rotor.strength[:, 2])
    hist_agreement = norm(x_hist - reference) / max(norm(reference), eps())
    hist_agreement <= 1e-8 || error("History solve disagrees with trials ($hist_agreement)")

    trials = filter(r -> r.phase == "trial", rows)
    times = [r.solve_seconds for r in trials]
    summary = Dict{String,Any}(
        "arm" => arm, "label" => label, "block" => block,
        "instrumented" => diag, "dagteam_workers" => workers,
        "dagteam_idle" => idle,
        "julia_threads" => Threads.nthreads(),
        "solves" => length(trials),
        "median_solve_seconds" => median(times),
        "min_solve_seconds" => minimum(times),
        "max_solve_seconds" => maximum(times),
        "iterations" => unique(r.iterations for r in trials),
        "history_iterations" => length(history),
        "history_agreement" => hist_agreement)
    cold_write_toml(joinpath(out, "stage1_summary.toml"), summary)
    cold_write_toml(joinpath(out, "status.toml"), Dict("status" => "completed"))
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && stage1_main()
