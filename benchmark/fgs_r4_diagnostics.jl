#!/usr/bin/env julia

include(joinpath(@__DIR__, "fgs_cold_common.jl"))
using Statistics
using Profile

function diag_seconds(d, key)
    Float64(get(d, key, 0)) * 1e-9
end

function write_census(solver, out)
    fgs = solver.fgs
    T = eltype(fgs.nonself_matrices.data)
    rows = NamedTuple[]
    for i in eachindex(fgs.nonself_matrices.sizes)
        m, n = fgs.nonself_matrices.sizes[i]
        push!(rows, (; leaf=i, m, n, matrix_elements=m*n,
            matrix_bytes=m*n*sizeof(T), scatter_entries=m,
            target_interactions=length(fgs.index_map[i]),
            source_strengths=length(fgs.strengths_by_leaf[i])))
    end
    cold_csv(joinpath(out, "gemv_census.csv"), rows)
    color_sizes = [length(x) for x in fgs.leaves_by_color]
    cold_write_toml(joinpath(out, "gemv_census.toml"), Dict(
        "element_type" => string(T),
        "element_bytes" => sizeof(T),
        "leaf_count" => length(rows),
        "nonempty_leaf_count" => count(r -> r.m > 0 && r.n > 0, rows),
        "matrix_elements" => sum(r.matrix_elements for r in rows),
        "matrix_bytes" => sum(r.matrix_bytes for r in rows),
        "rhs_bytes" => sizeof(fgs.nonself_matrices.rhs),
        "direct_interactions" => length(fgs.direct_list),
        "scatter_entries" => sum(r.scatter_entries for r in rows),
        "sweep_order" => string(fgs.sweep_order),
        "cached_leaf_lu" => fgs.leaf_lu_cache !== nothing,
        "color_count" => length(color_sizes),
        "color_sizes" => color_sizes))
end

function thread_snapshot()
    Sys.islinux() || return Dict{Int,Tuple{Int,Int}}()
    result = Dict{Int,Tuple{Int,Int}}()
    for path in readdir("/proc/self/task"; join=true)
        tid = tryparse(Int, basename(path)); tid === nothing && continue
        fields = split(read(joinpath(path, "stat"), String))
        result[tid] = (parse(Int, fields[14]) + parse(Int, fields[15]),
                       parse(Int, fields[39]))
    end
    result
end

function timed_row(solver, batch, trial, instrumented, reference, activity)
    diagnostics = instrumented ? Dict{Symbol,UInt64}() : nothing
    before = thread_snapshot()
    row, x = cold_trial(solver; diagnostics)
    after = thread_snapshot()
    for tid in sort!(collect(union(keys(before), keys(after))))
        b = get(before, tid, (0, -1)); a = get(after, tid, (0, -1))
        push!(activity, (; batch, trial, instrumented, tid,
            cpu_ticks=a[1]-b[1], cpu_before=b[2], cpu_after=a[2]))
    end
    agreement = norm(x-reference) / max(norm(reference), eps())
    cold_require(row.eligible && agreement <= 1e-8,
        "Diagnostic timing trial failed correctness/agreement gate")
    d = diagnostics === nothing ? Dict{Symbol,UInt64}() : diagnostics
    stage_sum = sum(get(d, k, 0) for k in
        (:initialization_ns, :fmm_ns, :influence_mapping_ns, :residual_ns,
         :leaf_solve_ns, :nonself_product_ns, :scatter_ns,
         :remaining_iteration_ns, :final_update_ns))
    return (; batch, trial, instrumented, solve_seconds=row.solve_seconds,
        allocated_bytes=row.allocated_bytes, gc_seconds=row.gc_seconds,
        iterations=row.iterations, estimated_inner_sweeps=row.estimated_inner_sweeps,
        estimated_fmm_passes=row.estimated_fmm_passes,
        relative_solution_delta=agreement, solved=row.solved,
        accepted=row.accepted, authoritative_evaluator=row.authoritative_evaluator,
        authoritative_rel_l2=row.authoritative_rel_l2,
        direct_fmm_delta=row.evaluator_delta,
        total_stage_seconds=diag_seconds(d, :total_ns),
        initialization_seconds=diag_seconds(d, :initialization_ns),
        fmm_seconds=diag_seconds(d, :fmm_ns),
        influence_mapping_seconds=diag_seconds(d, :influence_mapping_ns),
        residual_seconds=diag_seconds(d, :residual_ns),
        leaf_solve_seconds=diag_seconds(d, :leaf_solve_ns),
        nonself_product_seconds=diag_seconds(d, :nonself_product_ns),
        scatter_seconds=diag_seconds(d, :scatter_ns),
        remaining_iteration_seconds=diag_seconds(d, :remaining_iteration_ns),
        final_update_seconds=diag_seconds(d, :final_update_ns),
        exclusive_sum_seconds=Float64(stage_sum)*1e-9,
        unaccounted_seconds=diag_seconds(d, :total_ns)-Float64(stage_sum)*1e-9,
        outer_count=Int(get(d, :outer_count, 0)),
        sweep_count=Int(get(d, :sweep_count, 0)),
        leaf_visit_count=Int(get(d, :leaf_visit_count, 0)))
end

function history_control(solver, out, reference)
    records = Vector{NamedTuple}()
    histories = Vector{Vector{Float64}}()
    solutions = Vector{Vector{Float64}}()
    for instrumented in (false, true)
        residuals = Float64[]
        diagnostics = instrumented ? Dict{Symbol,UInt64}() : nothing
        cold_reset!(solver)
        pnl._solve!(rotor, solver; diagnostics,
            callback=(iteration, residual) -> push!(residuals, residual))
        x = copy(rotor.strength[:,2])
        e = cold_validate(x; crosscheck=true)
        agreement = norm(x-reference)/max(norm(reference), eps())
        cold_require(solver.solved && e.accepted && agreement <= 1e-8,
            "Instrumentation history control failed acceptance")
        push!(histories, residuals); push!(solutions, x)
        push!(records, (; instrumented, iterations=solver.niter,
            history_length=length(residuals), solution_delta=agreement,
            accepted=e.accepted, authoritative_evaluator=e.authoritative_evaluator,
            authoritative_rel_l2=e.authoritative_rel_l2,
            direct_rel_l2=e.direct_rel_l2, direct_fmm_delta=e.evaluator_delta))
    end
    same_history = histories[1] == histories[2]
    rel = norm(solutions[1]-solutions[2])/max(norm(solutions[1]), eps())
    cold_require(same_history && rel <= 1e-8,
        "Instrumentation changed convergence history or solution")
    cold_csv(joinpath(out, "instrumentation_equivalence.csv"),
        [(; records[1]..., counterpart_history_identical=same_history,
            counterpart_solution_rel_l2=rel),
         (; records[2]..., counterpart_history_identical=same_history,
            counterpart_solution_rel_l2=rel)])
end

function main()
    configs, out, stage = cold_initialize!()
    length(configs) == 1 || error("R4 diagnostics requires exactly one retained config")
    c = copy(only(configs))
    c["kind"] == "fgs" || error("R4 diagnostics requires FGS")
    cold_calibrate_selected!(c, out; selected=true)
    cold_write_toml(joinpath(out, "config.toml"), c)
    reset_cold!(); solver = cold_make(c)
    compile_row, _ = cold_trial(solver)
    cold_require(compile_row.eligible, "Compilation solve failed")
    warm, reference = cold_trial(solver; crosscheck=true)
    cold_require(warm.eligible, "Warmup failed")
    write_census(solver, out)
    history_control(solver, out, reference)

    rows = NamedTuple[]
    activity = NamedTuple[]
    reps = parse(Int, get(ENV, "COLD_DIAG_REPS", "10"))
    reps >= 10 || error("COLD_DIAG_REPS must be at least 10")
    for (batch, instrumented) in enumerate((false, true, false, true))
        for trial in 1:reps
            push!(rows, timed_row(solver, batch, trial, instrumented, reference, activity))
            cold_csv(joinpath(out, "diagnostic_trials.csv"), rows)
            cold_csv(joinpath(out, "thread_activity.csv"), activity)
        end
    end
    cold_reset!(solver)
    Profile.init(n=10^7, delay=0.001)
    Profile.clear()
    Profile.@profile pnl._solve!(rotor, solver)
    for format in (:flat, :tree)
        open(joinpath(out, "cpu_thread_complete_$(format).txt"), "w") do io
            Profile.print(io; format, C=true, groupby=[:thread, :task])
        end
    end
    cold_check_solution(solver, out, "cpu_thread_complete_validation"; reference)
    cold_write_toml(joinpath(out, "status.toml"), Dict("status"=>"completed"))
end

main()
