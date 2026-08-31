#=##############################################################################
Shared utilities for the BRAINSTORM 021 solver benchmark harnesses.

Implements the campaign's standing rulings (021 control doc):
- ruling 5/7: hand timing via `time_ns()`, min-of-k (k ≥ 5) after warmup,
  setup components timed separately from the per-step solve;
- ruling 6: thread-mode assertion + log banner (Julia and BLAS thread counts
  asserted at startup, written to the banner, and recorded in every CSV row);
- ruling 8: memory via `Base.summarysize` (resident solver state, provisional
  per decision_rules.md) and `@allocated` (per-solve allocation);
- CSV schema exactly as specified in
  BRAINSTORM/021_rotor_hover_solver_benchmarks/decision_rules.md.

Judge runs from these CSVs, never stdout (018 lesson).
=###############################################################################

import LinearAlgebra
import FLOWPanel as pnl
using Printf

################################################################################
# Thread-mode assertion + banner (ruling 6)
################################################################################

"""
    assert_and_banner(io=stdout) -> NamedTuple

Assert `Threads.nthreads() == ENV["EXPECT_JULIA_THREADS"]`, pin BLAS to one
thread in single mode (`ENV["THREADING_MODE"] == "single"`; BLAS defaults to
multithreaded even under `julia -t 1`, which would silently advantage the
dense-LU solvers), record the BLAS thread count in multi mode, and write the
log banner. Returns the banner facts for CSV rows, including the full banner
`text` so harnesses can persist it next to their CSVs (`banner.txt`).
"""
function assert_and_banner(io::IO=stdout)
    threading_mode = get(ENV, "THREADING_MODE", "single")
    threading_mode in ("single", "multi") ||
        error("THREADING_MODE must be \"single\" or \"multi\"; got \"$threading_mode\"")

    expect = parse(Int, get(ENV, "EXPECT_JULIA_THREADS", "1"))
    julia_threads = Threads.nthreads()
    julia_threads == expect || error(
        "Thread assertion failed: Threads.nthreads() = $julia_threads but " *
        "EXPECT_JULIA_THREADS = $expect. Launch with `julia -t $expect`.")

    # BLAS threads are pinned EXPLICITLY in BOTH modes (2026-08-24). Relying on
    # OPENBLAS_NUM_THREADS is fragile: OpenBLAS reads it at load time and
    # otherwise defaults to the machine's core count. On the campaign's
    # exclusive 128-core zen3 node that is 128 BLAS threads against Julia's 64
    # — a 2x oversubscription that would silently inflate every dense timing
    # (backslash_ldiv, and the Backslash reference solves), which is exactly
    # what Phase 1 is trying to measure. Asserting makes a mismatch a hard
    # failure instead of a value that is merely recorded and later believed.
    expected_blas = threading_mode == "single" ? 1 : julia_threads
    LinearAlgebra.BLAS.set_num_threads(expected_blas)
    blas_threads = LinearAlgebra.BLAS.get_num_threads()
    blas_threads == expected_blas || error(
        "BLAS thread pinning failed: got $blas_threads threads, expected " *
        "$expected_blas in $threading_mode mode")

    commit = _git_describe(@__DIR__)
    fm_commit = _git_describe(joinpath(@__DIR__, "..", "..", "FastMultipole"))
    hardware_tag = get(ENV, "HARDWARE_TAG", gethostname())
    julia_version = string(VERSION)

    # Filament regularization family (2026-08-22). The family is a GLOBAL
    # production default: RECORD it, never change it here. It was unrecorded
    # until now, which left several 2026-08-20/21 wake-carrying arms with an
    # indeterminate family (BRAINSTORM/025 landed selectability mid-day on
    # 08-20, so date alone cannot resolve them). Pattern lifted from
    # benchmark/fm051_pass_parity.jl:922.
    filament_reg = string(pnl.FILAMENT_REGULARIZATION[])

    text = join([
        "#"^70,
        "# 021 solver benchmark banner",
        "#   threading_mode = $threading_mode",
        "#   julia_threads  = $julia_threads",
        "#   blas_threads   = $blas_threads",
        "#   commit         = $commit",
        "#   fm_commit      = $fm_commit",
        "#   julia_version  = $julia_version",
        "#   hardware_tag   = $hardware_tag",
        "#   filament_reg   = $filament_reg",
        "#   date           = $(time_string())",
        "#"^70,
    ], "\n")
    println(io, text)

    return (; threading_mode, julia_threads, blas_threads,
              commit, fm_commit, julia_version, hardware_tag, filament_reg, text)
end

"Commit hash + dirty marker for the git repository containing `dir`."
function _git_describe(dir::AbstractString)
    sha = try
        readchomp(`git -C $dir rev-parse HEAD`)
    catch
        return "unknown"
    end
    dirty = try
        isempty(readchomp(`git -C $dir status --porcelain`)) ? "" : "-dirty"
    catch
        ""
    end
    return sha * dirty
end

time_string() = Libc.strftime("%Y-%m-%d %H:%M:%S", time())

################################################################################
# Timing (rulings 5, 7) and memory (ruling 8)
################################################################################

"""
    min_of_k(f; k=5, warmup=1, setup! = nothing) -> (t_min_seconds, result)

Run `f()` `warmup` times (compilation excluded), then `k` more times, timing
each with `time_ns()`. Returns the minimum time in seconds and the result of
the last timed call.

`setup!()` (untimed) runs before **every** call including warmups — use it to
reset solver-visible state (e.g. zero `body.strength`) so every rep is a
genuine cold solve. Iterative solvers that seed from current state (FGS reads
`body.strength` as its initial iterate) otherwise measure warm no-ops.
"""
function min_of_k(f; k::Int=5, warmup::Int=1, setup! = nothing)
    result = nothing
    for _ in 1:warmup
        setup! === nothing || setup!()
        result = f()
    end
    t_min = typemax(UInt64)
    for _ in 1:k
        setup! === nothing || setup!()
        t0 = time_ns()
        result = f()
        t_min = min(t_min, time_ns() - t0)
    end
    return t_min / 1e9, result
end

"""
    solver_state_bytes(solver) -> Int

Resident solver-state size in bytes with a **comparable object boundary**
across solver types (finalized ruling-8 method, 2026-08-12): total
`Base.summarysize(solver)` minus the summarysize of any `AbstractBody` the
solver references (directly, through nested operators/preconditioners, or as
a bodies tuple). `Base.summarysize` deduplicates shared references, so aliased
storage (e.g. `Backslash.G` ≡ `Glu.factors`) is counted once, and each
referenced body is counted exactly once in the total before subtraction.

In: matrices, factorizations, FMM trees/buffers, Krylov workspaces, scratch,
histories. Out: the body (mesh/strength/field arrays), which is not solver
state.
"""
function solver_state_bytes(solver)
    bodies = Base.IdSet{Any}()
    _collect_bodies!(bodies, solver, 0)
    total = Base.summarysize(solver)
    for body in bodies
        total -= Base.summarysize(body)
    end
    return total
end

function _collect_bodies!(found, obj, depth)
    depth > 4 && return nothing
    if obj isa pnl.AbstractBody
        push!(found, obj)
        return nothing
    end
    # don't walk into leaves or closures (any body captured by a closure is
    # also referenced by a direct field on these solver types)
    (obj isa Number || obj isa AbstractArray || obj isa AbstractString ||
     obj isa Symbol || obj isa Function || obj === nothing) && return nothing
    if obj isa Tuple
        for el in obj
            _collect_bodies!(found, el, depth + 1)
        end
        return nothing
    end
    isstructtype(typeof(obj)) || return nothing
    for fn in fieldnames(typeof(obj))
        isdefined(obj, fn) && _collect_bodies!(found, getfield(obj, fn), depth + 1)
    end
    return nothing
end

################################################################################
# CSV schema (decision_rules.md, finalized Phase 0 W4)
################################################################################

const RUNS_CSV_COLUMNS = (
    :run_id, :phase, :solver_config, :mesh_file, :n_panels, :threading_mode,
    :julia_threads, :blas_threads, :t_assembly, :t_factorize, :t_tree,
    :t_precond, :t_rhs, :t_solve_min, :k_reps, :iterations, :rms_residual,
    :max_residual, :mem_state_bytes, :alloc_solve_bytes, :commit, :fm_commit,
    :julia_version, :hardware_tag, :filament_reg, :solver_settings,
    :backend_settings, :notes,
)

runs_csv_header() = join(String.(collect(RUNS_CSV_COLUMNS)), ",")

"""
    settings_string(dict) -> String

Flatten a (possibly nested) metadata Dict — e.g. from
`pnl._solver_metadata_dict`/`pnl._backend_metadata_dict` — into a compact,
CSV-safe `k=v;k.sub=v;…` string with sorted keys.
"""
function settings_string(dict::AbstractDict; prefix::String="")
    parts = String[]
    for key in sort!(collect(keys(dict)))
        v = dict[key]
        name = isempty(prefix) ? String(key) : "$prefix.$key"
        if v isa AbstractDict
            sub = settings_string(v; prefix=name)
            isempty(sub) || push!(parts, sub)
        else
            push!(parts, "$name=$v")
        end
    end
    return join(parts, ";")
end

_csv_cell(x::Nothing) = ""
_csv_cell(x::AbstractFloat) = @sprintf("%.9g", x)
_csv_cell(x::AbstractString) = occursin(r"[,\"\n]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : x
_csv_cell(x) = string(x)

"""
    write_run_row!(io; kwargs...)

Write one `runs.csv` row. Keywords must be a subset of `RUNS_CSV_COLUMNS`;
missing columns are left empty (`nothing` also maps to empty — use it for
non-applicable cost components, e.g. `t_factorize` for matrix-free solvers).
"""
function write_run_row!(io::IO; kwargs...)
    unknown = setdiff(keys(kwargs), RUNS_CSV_COLUMNS)
    isempty(unknown) || error("Unknown runs.csv columns: $unknown")
    row = join((_csv_cell(get(kwargs, c, nothing)) for c in RUNS_CSV_COLUMNS), ",")
    println(io, row)
    return nothing
end

"""
    write_history_csv(path, history; residual_true=nothing)

Write a per-iteration convergence-history sidecar
(`iter, t_wall, residual_internal, residual_true`), `t_wall` in seconds since
the history's `reset!`. `residual_true` is optional (same length as the
history) — the internal metric stays labeled by `history.metric` in the
header comment.
"""
function write_history_csv(path::AbstractString, history;
                           residual_true=nothing)
    open(path, "w") do io
        println(io, "# metric = $(history.metric)")
        println(io, "iter,t_wall,residual_internal,residual_true")
        t_wall = (history.t_ns .- history.t0_ns) ./ 1e9
        for i in eachindex(history.iter)
            rt = residual_true === nothing ? "" : _csv_cell(residual_true[i])
            println(io, "$(history.iter[i]),$(_csv_cell(t_wall[i]))," *
                        "$(_csv_cell(history.residual_internal[i])),$rt")
        end
    end
    return nothing
end

################################################################################
# Certified BC-error evaluator (021 Phase 1, Ryan ruling 2026-08-15 #1)
################################################################################

"""
    bc_error!(body, x; rms_b, target_rel=1e-6, safety=0.1,
              max_expansion_order=20, multipole_acceptance=0.5, leaf_size=20,
              backend=:fmm) -> NamedTuple

Primary campaign metric (Ryan 2026-08-15): the boundary-condition error of a
candidate Dirichlet solution `x` — the interior perturbation-potential
residual at the control points — as a relative L2 against `rms_b`
(= RMS of the source potential φ_σ, since b = −φ_σ).

One influence pass with BOTH strength columns loaded: σ re-derived from the
frozen BC via `set_strengths!` (contract: `body.velocity` must hold the
apparent velocity at the control points on entry) and column 2 = `x`. The
control-point potential then IS the BC residual φ_σ + G_μ x = G x − b, so no
reference solution and no dense G are needed at any rung.

`backend=:fmm` (production) evaluates with the FMM under a dynamic
`PowerAbsolutePotential(safety × target_rel × rms_b)` error tolerance
(`max_expansion_order` = the dynamic-P cap) so evaluation error inflates the
metric ≤ `safety` at the target scale; `error_success=false` in the return
means some M2L pair could not certify its bound at the P cap — treat the
metric as UNCERTIFIED. `backend=:direct` is the O(N²) validation path
(exact near-field kernel; `error_success` trivially true).

Frozen shared BC definition (do not vary): perturbation gauge (φ∞ excluded —
whatever `body.potential` held on entry is ignored and restored), σ from the
apparent velocity, panel core_size (forced during evaluation). Body
state (`strength`, `velocity`, `potential`, `core_size`) is fully
restored — the evaluation is side-effect-free.

Pass `phi_out` (length `body.ncells`) to also receive the per-panel BC residual
φ itself, copied out before the body state is restored. Callers that want order
statistics (min/quartiles/max) rather than just the two summary norms need this;
it costs one copy and nothing else. Added 2026-08-25 for the Phase 3 per-step
BC check.

Returns `(; rel_l2, rel_max, error_success, t_eval, epsilon_requested)` where
`rel_l2 = RMS(φ)/rms_b` and `rel_max = max|φ|/rms_b`.
"""
function bc_error!(body::pnl.AbstractBody, x::AbstractVector; rms_b::Real,
                   target_rel::Real=1e-6, safety::Real=0.1,
                   max_expansion_order::Int=20, multipole_acceptance::Real=0.5,
                   leaf_size::Int=20, backend::Symbol=:fmm,
                   phi_out::Union{Nothing,AbstractVector}=nothing)
    pnl.has_dirichlet_bc(body) ||
        error("bc_error! implements the Dirichlet BC metric; got a Neumann body")
    strength_old = copy(body.strength)
    velocity_old = copy(body.velocity)
    potential_old = copy(body.potential)
    core_size_old = body.core_size
    epsilon_requested = safety * target_rel * rms_b
    local rel_l2, rel_max, error_success, t_eval
    try
        body.core_size = body.core_size_panel
        pnl.set_strengths!(body)              # col 1 = σ from BC, col 2 zeroed
        body.strength[:, 2] .= x
        body.potential .= 0
        if backend === :fmm
            t_eval = @elapsed outputs = pnl.FastMultipole.fmm!((body,), (body,);
                expansion_order=max_expansion_order,
                multipole_acceptance=Float64(multipole_acceptance),
                leaf_size_source=Int(leaf_size),
                error_tolerance=pnl.FastMultipole.PowerAbsolutePotential(epsilon_requested),
                scalar_potential=true, gradient=false, hessian=false,
                shrink=true)
            error_success = outputs[end]      # fmm! returns error_success last
        elseif backend === :direct
            t_eval = @elapsed pnl.influence!(body, body, pnl.DirectBackend();
                scalar_potential=true, velocity=false)
            error_success = true
        else
            error("bc_error!: unknown backend $backend (use :fmm or :direct)")
        end
        rel_l2 = (LinearAlgebra.norm(body.potential) / sqrt(body.ncells)) / rms_b
        rel_max = maximum(abs, body.potential) / rms_b
        phi_out === nothing || (phi_out .= body.potential)
    finally
        body.strength .= strength_old
        body.velocity .= velocity_old
        body.potential .= potential_old
        body.core_size = core_size_old
    end
    return (; rel_l2, rel_max, error_success, t_eval, epsilon_requested)
end

"""
    validate_runs_csv(path) -> Int

Re-parse a `runs.csv` and assert schema conformance (column names/count and
finite required numeric fields). Returns the number of data rows.
"""
function validate_runs_csv(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("Empty runs.csv at $path")
    header = split(lines[1], ",")
    expected = String.(collect(RUNS_CSV_COLUMNS))
    header == expected || error(
        "runs.csv header mismatch at $path:\n  got      $header\n  expected $expected")
    required = (:run_id, :solver_config, :mesh_file, :n_panels, :threading_mode,
                :julia_threads, :blas_threads, :t_solve_min, :rms_residual,
                :max_residual, :commit, :fm_commit, :julia_version,
                :filament_reg, :solver_settings)
    idx = Dict(Symbol(c) => i for (i, c) in enumerate(expected))
    nrows = 0
    for (li, line) in enumerate(lines[2:end])
        isempty(strip(line)) && continue
        cells = split(line, ",")  # harness cells contain no quoted commas except notes (last column)
        length(cells) >= length(expected) || error("Row $(li+1) has $(length(cells)) cells")
        for c in required
            isempty(cells[idx[c]]) && error("Row $(li+1): required column $c is empty")
        end
        for c in (:t_solve_min, :rms_residual, :max_residual)
            v = tryparse(Float64, cells[idx[c]])
            (v === nothing || !isfinite(v)) && error("Row $(li+1): column $c = '$(cells[idx[c]])' is not finite")
        end
        nrows += 1
    end
    return nrows
end
