# Staged Dirichlet wake-coupling diagnostic.
#
# Task 1 supplies the semi-infinite reference. Task 2 exercises the existing
# finite-wake velocity-through-sources route. Task 3 evaluates the prescribed
# finite PanelWake scalar potential directly. Neither changes production code.

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: dot, norm, ldiv!, svd, svdvals
import LinearAlgebra.BLAS
import Printf: @printf, @sprintf
import TOML
import SHA
import Statistics: mean

include(joinpath(@__DIR__, "..", "..", "examples", "pitching_wing.jl"))

const DIRICHLET_DATA_PATH = joinpath(@__DIR__, "..", "..", "data", "dirichlet_solve")
const DIRICHLET_COMPARISON_CSV = joinpath(DIRICHLET_DATA_PATH, "comparison.csv")
const DEFAULT_ALPHA_DEG = 3.94

function _cli_options(args)
    task = nothing
    alpha_deg = DEFAULT_ALPHA_DEG
    output_dir = DIRICHLET_DATA_PATH
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--task"
            i < length(args) || throw(ArgumentError("--task requires a value"))
            task = lowercase(args[i + 1])
            i += 2
        elseif startswith(arg, "--task=")
            task = lowercase(split(arg, "="; limit=2)[2])
            i += 1
        elseif arg == "--alpha-deg"
            i < length(args) || throw(ArgumentError("--alpha-deg requires a value"))
            alpha_deg = parse(Float64, args[i + 1])
            i += 2
        elseif startswith(arg, "--alpha-deg=")
            alpha_deg = parse(Float64, split(arg, "="; limit=2)[2])
            i += 1
        elseif arg == "--output-dir"
            i < length(args) || throw(ArgumentError("--output-dir requires a value"))
            output_dir = args[i + 1]
            i += 2
        elseif startswith(arg, "--output-dir=")
            output_dir = split(arg, "="; limit=2)[2]
            i += 1
        elseif !startswith(arg, "--") && isnothing(task)
            task = lowercase(arg)
            i += 1
        else
            throw(ArgumentError("unknown argument '$arg'"))
        end
    end
    isfinite(alpha_deg) || throw(ArgumentError("--alpha-deg must be finite"))
    return (; task=something(task, "task1"), alpha_deg, output_dir)
end

function _build_kutta_map(body)
    nedge = sum(size(s, 2) for s in body.shedding)
    C = zeros(Float64, nedge, body.ncells)
    row = 0
    for shedding in body.shedding, edge in axes(shedding, 2)
        row += 1
        upper = shedding[1, edge]
        lower = shedding[4, edge]
        1 <= upper <= body.ncells || error("invalid upper shedding-panel index $upper")
        1 <= lower <= body.ncells || error(
            "Task 1 requires paired shedding edges; edge $row has lower index $lower")
        C[row, upper] = 1.0
        C[row, lower] = -1.0
    end
    return C
end

function _common_exterior_probes(c, b)
    points = SVector{3,Float64}[]
    for x_c in (-0.5, 0.5, 1.5, 3.0), y_b in (0.0, 0.35), z_c in (-0.5, 0.5)
        push!(points, SVector{3,Float64}(x_c * c, y_b * b, z_c * c))
    end
    return points
end

function _direct_exterior_velocity(body, Uinf, points)
    switch = pnl.FastMultipole.DerivativesSwitch(false, true, false)
    velocities = zeros(Float64, 3, length(points))
    for (ip, target) in pairs(points)
        velocity = SVector{3,Float64}(Uinf)
        for source in 1:body.ncells
            _, induced_velocity, _ = pnl.induced(target, body, source, switch;
                kerneloffset=body.kerneloffset_targets)
            velocity += induced_velocity
        end
        velocities[:, ip] .= velocity
    end
    all(isfinite, velocities) || error("non-finite exterior velocity in Task 1")
    return velocities
end

function _write_probe_csv(path, points, velocities)
    open(path, "w") do io
        println(io, "probe,x_m,y_m,z_m,u_m_s,v_m_s,w_m_s")
        for i in eachindex(points)
            p = points[i]
            u = view(velocities, :, i)
            println(io, join((i, p[1], p[2], p[3], u[1], u[2], u[3]), ','))
        end
    end
end

function _write_comparison_csv(path, row)
    header = propertynames(row)
    open(path, "w") do io
        println(io, join(header, ','))
        println(io, join((getproperty(row, key) for key in header), ','))
    end
end

const _COMPARISON_COLUMNS = (
    :task, :case, :formulation, :status, :frozen_wake_step,
    :wake_length_c, :wake_rows, :ncells, :nshedding_edges,
    :linear_residual_abs, :linear_residual_rel, :kutta_constant_error,
    :gamma_min, :gamma_max, :gamma_mean, :gamma_l2,
    :transition_free_mismatch_l2, :transition_free_mismatch_rel,
    :CL, :CL_minus_task1, :relative_CL_difference,
    :relative_exterior_velocity_difference,
    :CL_minus_task2, :relative_CL_difference_task2,
    :relative_exterior_velocity_difference_task2,
    :q_f_l2, :Cq_f_l2, :fixed_point_iterations,
    :oracle_linear_residual_rel, :relative_oracle_mu_difference,
    :relative_oracle_CL_difference,
    # Task 4/5 (formulation validation) additive columns.  Old rows leave these
    # blank; no existing column is dropped or reinterpreted.
    :CL_minus_task3_single, :CL_minus_task3_iterated,
    :signed_post_correction_diff, :f_semiinf, :f_oracle,
    :q_l2, :green_q_vs_task3_trace_rel,
    :green_bordered_residual_rel, :gauge_row_abs,
    :explicit_potential_residual_rel, :blas_threads,
    # Task 5 line-integral trace diagnostics.
    :estimator, :quadrature, :interior_path, :path_depth,
    :c_error_rel, :c_error_l2, :c_error_max, :c_l2, :c_oracle_l2,
    :quadrature_spread_l2, :path_spread_l2,
    :circulation_consistency_abs, :circulation_consistency_rel,
    :state_init_elapsed_s, :eligible_solve_elapsed_s, :oracle_elapsed_s,
    :total_elapsed_s,
)

_csv_string(x) = x === missing ? "" : string(x)

function _read_comparison_rows(path)
    isfile(path) || return Dict{String,Dict{Symbol,String}}()
    lines = readlines(path)
    isempty(lines) && return Dict{String,Dict{Symbol,String}}()
    header = Symbol.(split(first(lines), ','))
    rows = Dict{String,Dict{Symbol,String}}()
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        values = split(line, ','; keepempty=true)
        row = Dict{Symbol,String}(zip(header, values))
        formulation = get(row, :formulation, "")
        case = get(row, :case, formulation == "semi_infinite_attached_wake" ? "baseline" : "")
        rows[formulation * "|" * case] = row
    end
    return rows
end

function _upsert_comparison(path, row)
    rows = _read_comparison_rows(path)
    incoming = Dict{Symbol,String}(
        key => _csv_string(hasproperty(row, key) ? getproperty(row, key) : missing)
        for key in _COMPARISON_COLUMNS)
    key = incoming[:formulation] * "|" * incoming[:case]
    rows[key] = incoming
    open(path, "w") do io
        println(io, join(_COMPARISON_COLUMNS, ','))
        ordered = sort(collect(values(rows)); by=r -> (
            tryparse(Int, get(r, :task, "")) === nothing ? typemax(Int) : parse(Int, r[:task]),
            get(r, :case, "")))
        for stored in ordered
            println(io, join((get(stored, column, "") for column in _COMPARISON_COLUMNS), ','))
        end
    end
    return path
end

function _write_task1_configuration(path; c, b, U, Uinf, dt, alpha_deg)
    config = Dict(
        "task" => 1,
        "formulation" => "semi_infinite_attached_wake",
        "geometry" => Dict(
            "airfoil" => "NACA 0015",
            "chord_m" => c,
            "span_m" => b,
            "aspect_ratio" => b / c,
            "n_airfoil" => 161,
            "n_span" => 13,
            "n_endcap" => 9,
            "watertight" => true,
        ),
        "flow" => Dict(
            "Uinf_ft_s" => 330.2,
            "Uinf_m_s" => U,
            "alpha_deg" => alpha_deg,
            "Uinf_vector_m_s" => collect(Uinf),
            "rho_kg_m3" => 1.225,
        ),
        "time_scale" => Dict("c_per_dt" => 0.5, "dt_s" => dt),
        "free_wake" => Dict("present" => false, "frozen_step" => -1),
        "backend" => "DirectBackend",
    )
    open(path, "w") do io
        TOML.print(io, config; sorted=true)
    end
end

function _lu_matrix_product(factorization, x)
    # Julia's LU satisfies L*U == A[p, :].  Undo that row permutation after
    # applying the factors so this is the original assembled operator times x.
    permuted = factorization.L * (factorization.U * x)
    product = similar(permuted)
    product[factorization.p] .= permuted
    return product
end

function run_task1(; path=DIRICHLET_DATA_PATH, alpha_deg=DEFAULT_ALPHA_DEG)
    mkpath(path)

    c = 1.0 * FT_TO_M
    b = 4.0 * c
    U = 330.2 * FT_TO_M
    rho = 1.225
    Uinf = _uinf_from_alpha(U, alpha_deg)
    dt = c / U * 0.5
    Sref = b * c

    body = build_pitching_wing_body(c, b;
        n_span=13, n_airfoil=161, n_endcap=9, endcap=:round,
        semiinfinite_wake=true)
    set_wake_Das!(body, Uinf)

    watertight, _ = pnl.iswatertight(body.nodes, body.cells)
    watertight || error("Task 1 body is not watertight")
    C = _build_kutta_map(body)
    kutta_constant_error = norm(C * ones(body.ncells), Inf)
    kutta_constant_error == 0.0 || error(
        "paired-edge Kutta map does not annihilate constants: $kutta_constant_error")

    direct = pnl.DirectBackend()
    solver = pnl.Backslash(body)
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false,
        backend=direct)
    force = pnl.ForceMonitor(1, 1;
        i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c),
        correct_kuttacondition=false,
        verbose=false,
    )
    run_name = "dirichlet_task1_semiinfinite"
    pnl.steady!(body, pnl.ReferenceFrame(body), Uinf;
        body_solvers=solver,
        backend=direct,
        backend_solve=direct,
        backend_system=direct,
        monitors=(pressure, force),
        path=path,
        name=run_name,
        dt=dt,
        clean_files=true,
        verbose=true,
    )

    mu = copy(view(body.strength, :, 2))
    gamma = C * mu
    residual = _lu_matrix_product(solver.Glu, mu) - solver.rhs
    residual_abs = norm(residual)
    residual_rel = residual_abs / max(norm(solver.rhs), eps(Float64))
    all(isfinite, mu) || error("non-finite Task 1 doublet solution")
    isfinite(residual_rel) || error("non-finite Task 1 linear residual")
    residual_rel <= 1e-10 || error("Task 1 relative residual $residual_rel exceeds 1e-10")

    Lhat, _, _ = _lift_drag_span_directions(Uinf)
    CL = dot(SVector{3,Float64}(force.force[:, 1]), Lhat)
    isfinite(CL) || error("non-finite Task 1 steady CL")

    probes = _common_exterior_probes(c, b)
    probe_velocity = _direct_exterior_velocity(body, Uinf, probes)
    _write_probe_csv(joinpath(path, "task1_exterior_probes.csv"), probes, probe_velocity)

    comparison = (;
        task=1,
        case="baseline",
        formulation="semi_infinite_attached_wake",
        status="reference",
        frozen_wake_step=-1,
        wake_length_c=Inf,
        wake_rows=0,
        ncells=body.ncells,
        nshedding_edges=size(C, 1),
        linear_residual_abs=residual_abs,
        linear_residual_rel=residual_rel,
        kutta_constant_error,
        gamma_min=minimum(gamma),
        gamma_max=maximum(gamma),
        gamma_mean=sum(gamma) / length(gamma),
        gamma_l2=norm(gamma),
        transition_free_mismatch_l2=NaN,
        transition_free_mismatch_rel=NaN,
        CL,
        CL_minus_task1=0.0,
        relative_CL_difference=0.0,
        relative_exterior_velocity_difference=NaN,
    )
    _upsert_comparison(joinpath(path, "comparison.csv"), comparison)
    _write_comparison_csv(joinpath(path, "task1_convergence.csv"), (;
        task=1, solve_index=1, linear_residual_abs=residual_abs,
        linear_residual_rel=residual_rel))
    _write_task1_configuration(joinpath(path, "task1.config.toml");
        c, b, U, Uinf, dt, alpha_deg)
    open(joinpath(path, "task1.metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "theory_file" => "docs/wake_solve_schemes.md",
            "theory_section" => "Discrete diagnostic construction: Task 1 semi-infinite baseline",
            "run_name" => run_name,
            "vtk_collection" => run_name * "_body1.pvd",
            "comparison_csv" => "comparison.csv",
            "probe_csv" => "task1_exterior_probes.csv",
            "alpha_deg" => alpha_deg,
        ); sorted=true)
    end

    @printf("Task 1 semi-infinite baseline\n")
    @printf("  panels / shedding edges: %d / %d\n", body.ncells, size(C, 1))
    @printf("  relative linear residual: %.8e\n", residual_rel)
    @printf("  transition circulation: min %.8e, max %.8e, L2 %.8e\n",
        minimum(gamma), maximum(gamma), norm(gamma))
    @printf("  steady CL: %.8g\n", CL)
    println("  exterior comparison: N/A (baseline probe field saved)")
    println("  output: $(abspath(path))")
    return (; body, solver, C, gamma, CL, residual_rel, probes, probe_velocity)
end

function _task1_grid_case(n_airfoil, n_span, n_endcap;
        U=330.2 * FT_TO_M, alpha_deg=3.94,
        c=1.0 * FT_TO_M, b=4.0 * c, rho=1.225,
        path=nothing, name="dirichlet_task1_grid_case")
    Uinf = _uinf_from_alpha(U, alpha_deg)
    Sref = b * c

    body = build_pitching_wing_body(c, b;
        n_span, n_airfoil, n_endcap, endcap=:round, semiinfinite_wake=true)
    set_wake_Das!(body, Uinf)
    C = _build_kutta_map(body)
    norm(C * ones(body.ncells), Inf) == 0.0 || error("invalid paired Kutta map")

    direct = pnl.DirectBackend()
    solver = pnl.Backslash(body)
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false,
        backend=direct)
    force = pnl.ForceMonitor(1, 1;
        i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c),
        correct_kuttacondition=false,
        verbose=false,
    )
    elapsed = @elapsed pnl.steady!(body, pnl.ReferenceFrame(body), Uinf;
        body_solvers=solver,
        backend=direct,
        backend_solve=direct,
        backend_system=direct,
        monitors=(pressure, force),
        path,
        name,
        dt=c / U * 0.5,
        verbose=false,
    )
    mu = copy(view(body.strength, :, 2))
    residual_abs = norm(_lu_matrix_product(solver.Glu, mu) - solver.rhs)
    residual_rel = residual_abs / max(norm(solver.rhs), eps(Float64))
    Lhat, _, _ = _lift_drag_span_directions(Uinf)
    CL = dot(SVector{3,Float64}(force.force[:, 1]), Lhat)
    all(isfinite, (CL, residual_abs, residual_rel)) || error(
        "non-finite Task 1 grid-convergence result")
    return (; n_airfoil, n_span, n_endcap, ncells=body.ncells,
        nshedding_edges=size(C, 1), CL, residual_abs, residual_rel,
        solve_elapsed_s=elapsed)
end

function run_task1_grid_convergence(; path=DIRICHLET_DATA_PATH)
    mkpath(path)
    grids = (
        (n_airfoil=81, n_span=7, n_endcap=5),
        (n_airfoil=121, n_span=10, n_endcap=7),
        (n_airfoil=161, n_span=13, n_endcap=9),
        (n_airfoil=201, n_span=16, n_endcap=11),
    )
    rows = NamedTuple[]
    for grid in grids
        @printf("Task 1 grid: n_airfoil=%d, n_span=%d, n_endcap=%d\n",
            grid.n_airfoil, grid.n_span, grid.n_endcap)
        row = _task1_grid_case(grid.n_airfoil, grid.n_span, grid.n_endcap)
        push!(rows, row)
        @printf("  panels=%d, CL=%.10f, residual=%.3e, solve=%.2f s\n",
            row.ncells, row.CL, row.residual_rel, row.solve_elapsed_s)
        GC.gc()
    end
    csv_path = joinpath(path, "task1_grid_convergence.csv")
    header = propertynames(first(rows))
    open(csv_path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join((getproperty(row, key) for key in header), ','))
        end
    end
    println("Task 1 grid convergence written to $(abspath(csv_path))")
    return rows
end

_relative_difference(value, reference) =
    (value - reference) / max(abs(reference), eps(Float64))

function _array_hash(A)
    bytes = reinterpret(UInt8, vec(copy(A)))
    return bytes2hex(SHA.sha256(bytes))
end

function _state_snapshot(body, wake)
    return (;
        body_nodes=copy(body.nodes),
        body_cells=copy(body.cells),
        body_Das=[copy(Das) for Das in body.Das],
        body_strength=copy(body.strength),
        wake_nodes=[copy(nodes) for nodes in wake.nodes],
        wake_strength=[copy(strength) for strength in wake.strength],
        wake_rows=wake.nwakes[],
        wake_options=(;
            core_size=wake.core_size,
            shed_with_induced_velocity=wake.shed_with_induced_velocity,
            unsteady_filament=wake.unsteady_filament,
            include_final_filament=wake.include_final_filament,
        ),
    )
end

function _assert_frozen_state(snapshot, body, wake; label)
    isequal(body.nodes, snapshot.body_nodes) || error("$label changed body nodes")
    isequal(body.cells, snapshot.body_cells) || error("$label changed body connectivity")
    all(isequal(a, b) for (a, b) in zip(body.Das, snapshot.body_Das)) ||
        error("$label changed attached transition geometry")
    wake.nwakes[] == snapshot.wake_rows || error("$label changed wake row count")
    all(isequal(a, b) for (a, b) in zip(wake.nodes, snapshot.wake_nodes)) ||
        error("$label changed prescribed wake nodes")
    all(isequal(a, b) for (a, b) in zip(wake.strength, snapshot.wake_strength)) ||
        error("$label changed prescribed wake strengths")
    options = (;
        core_size=wake.core_size,
        shed_with_induced_velocity=wake.shed_with_induced_velocity,
        unsteady_filament=wake.unsteady_filament,
        include_final_filament=wake.include_final_filament,
    )
    isequal(options, snapshot.wake_options) || error("$label changed wake options")
    size(body.strength, 2) == 2 || error("$label has unexpected body strength fields")
    return true
end

function _snapshot_checksums(snapshot)
    return Dict(
        "body_nodes_sha256" => _array_hash(snapshot.body_nodes),
        "body_cells_sha256" => _array_hash(snapshot.body_cells),
        "body_Das_sha256" => [_array_hash(A) for A in snapshot.body_Das],
        "body_strength_before_sha256" => _array_hash(snapshot.body_strength),
        "wake_nodes_sha256" => [_array_hash(A) for A in snapshot.wake_nodes],
        "wake_strength_sha256" => [_array_hash(A) for A in snapshot.wake_strength],
        "wake_rows" => snapshot.wake_rows,
        "wake_options" => Dict(string(k) => v for (k, v) in pairs(snapshot.wake_options)),
    )
end

function _write_rows_csv(path, rows)
    isempty(rows) && error("cannot write empty CSV $path")
    header = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            propertynames(row) == header || error("inconsistent CSV row schema for $path")
            println(io, join((getproperty(row, key) for key in header), ','))
        end
    end
    return path
end

function _finite_body(c, b, Uinf; das_length_c=0.5)
    body = build_pitching_wing_body(c, b;
        n_span=13, n_airfoil=161, n_endcap=9, endcap=:round,
        semiinfinite_wake=false)
    set_wake_Das!(body, Uinf; magnitude=das_length_c * c)
    return body
end

function _flat_wake(body, gamma, length_c, c; free_row_length_c=0.5, core_size=1e-3)
    nfree = round(Int, 2length_c - 1)
    nfree >= 1 || error("flat composite wake needs at least one free row")
    wake = pnl.PanelWake(body; nwakerows=nfree,
        include_final_filament=false, core_size)
    pnl.update_TE!(wake, body)
    Das = view(body.Das[1], :, 1)
    direction = Das ./ norm(Das) .* (free_row_length_c * c)
    for nodes in wake.nodes
        first_row = copy(view(nodes, :, 1, :))
        for row in 1:nfree+1
            view(nodes, :, row, :) .= first_row .+ (row - 1) .* direction
        end
    end
    for strength in wake.strength
        for row in 1:nfree
            view(strength, 1, row, :) .= gamma
        end
    end
    wake.nwakes[] = nfree
    return wake
end

function _exterior_velocity(body, wake, Uinf, points, backend)
    probes = pnl.FastMultipole.ProbeSystem(length(points), Float64)
    for i in eachindex(points)
        probes.position[i] = points[i]
        probes.gradient[i] = zero(SVector{3,Float64})
    end
    sources = isnothing(wake) ? (body,) : (body, pnl.get_sources(wake)...)
    pnl.influence!(probes, sources, backend;
        scalar_potential=false, velocity=true)
    velocities = zeros(Float64, 3, length(points))
    for i in eachindex(points)
        velocities[:, i] .= Uinf .+ probes.gradient[i]
    end
    all(isfinite, velocities) || error("non-finite finite-wake exterior velocity")
    return velocities
end

function _monitor_lift!(body, wake, frames, Uinf, rho, Sref, c, dt, backend;
        i_step=0, nt=1)
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false,
        backend, file=false)
    force = pnl.ForceMonitor(nt, 1;
        i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c),
        correct_kuttacondition=false,
        verbose=false,
        file=false,
    )
    context = pnl.MonitorContext()
    pnl._run_monitor!(pressure, context, (body,), (wake,), frames, Uinf,
        i_step, dt, i_step * dt)
    pnl._run_monitor!(force, context, (body,), (wake,), frames, Uinf,
        i_step, dt, i_step * dt)
    Lhat, _, _ = _lift_drag_span_directions(Uinf)
    return dot(SVector{3,Float64}(force.force[:, i_step + 1]), Lhat), pressure, force
end

function _frozen_solve(master_body, master_wake, solver, C, Uinf, rho, Sref, c,
        dt, probes, baseline_probe; label, backend=pnl.DirectBackend())
    body = deepcopy(master_body)
    wake = deepcopy(master_wake)
    snapshot = _state_snapshot(body, wake)
    frames = pnl.ReferenceFrame(body)
    pnl._steady_aerodynamics!(body, (body,), (wake,), frames, Uinf, solver;
        backend_wake=backend,
        backend_solve=backend,
        backend_system=backend,
        update_trailing_edges=false)
    _assert_frozen_state(snapshot, body, wake; label)

    mu = copy(view(body.strength, :, 2))
    gamma = C * mu
    all(isfinite, body.strength) || error("$label produced non-finite body strengths")
    residual = _lu_matrix_product(solver.Glu, mu) - solver.rhs
    residual_abs = norm(residual)
    residual_rel = residual_abs / max(norm(solver.rhs), eps(Float64))
    residual_rel <= 1e-10 || error("$label relative residual $residual_rel exceeds 1e-10")

    CL, pressure, force = _monitor_lift!(body, wake, frames, Uinf, rho, Sref, c,
        dt, backend)
    isfinite(CL) || error("$label produced non-finite CL")
    first_free = vec(copy(view(wake.strength[1], 1, 1, :)))
    mismatch = gamma - first_free
    mismatch_l2 = norm(mismatch)
    mismatch_rel = mismatch_l2 / max(norm(first_free), eps(Float64))
    velocity = _exterior_velocity(body, wake, Uinf, probes, backend)
    probe_rel = norm(velocity - baseline_probe) /
        max(norm(baseline_probe), eps(Float64))
    all(isfinite, (residual_abs, residual_rel, mismatch_l2, mismatch_rel, probe_rel)) ||
        error("$label produced a non-finite diagnostic")
    return (;
        body, wake, snapshot, pressure, force, gamma, CL,
        residual_abs, residual_rel, mismatch_l2, mismatch_rel,
        probe_velocity=velocity, probe_rel,
    )
end

function _task2_baseline(path; alpha_deg=DEFAULT_ALPHA_DEG)
    baseline = run_task1(; path, alpha_deg)
    if alpha_deg == DEFAULT_ALPHA_DEG
        isapprox(baseline.CL, 0.2747643938; atol=5e-11, rtol=0) || error(
            "Task 1 baseline changed: CL=$(baseline.CL)")
    end
    return baseline
end

function _finite_solver_template(c, b, Uinf; das_length_c=0.5)
    body = _finite_body(c, b, Uinf; das_length_c)
    return body, pnl.Backslash(body), _build_kutta_map(body)
end

function run_task2_flat(; path=DIRICHLET_DATA_PATH, baseline=nothing,
        finite_setup=nothing, das_length_c=0.5, free_row_length_c=0.5,
        alpha_deg=DEFAULT_ALPHA_DEG,
        case="flat", artifact_stem="task2_flat",
        formulation="finite_flat_velocity_through_sources")
    mkpath(path)
    c = 1.0 * FT_TO_M
    b = 4.0 * c
    U = 330.2 * FT_TO_M
    Uinf = _uinf_from_alpha(U, alpha_deg)
    rho = 1.225
    dt = 0.5c / U
    Sref = b * c
    baseline = isnothing(baseline) ? _task2_baseline(path; alpha_deg) : baseline
    template, solver, C = isnothing(finite_setup) ?
        _finite_solver_template(c, b, Uinf; das_length_c) : finite_setup
    actual_das_length_c = norm(view(template.Das[1], :, 1)) / c
    isapprox(actual_das_length_c, das_length_c; atol=10eps(Float64), rtol=0) ||
        error("finite solver template has Das/c=$actual_das_length_c, expected $das_length_c")
    norm(C * ones(template.ncells), Inf) == 0.0 || error("invalid finite Kutta map")

    lengths = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
    rows = NamedTuple[]
    last_result = nothing
    checksums = Dict{String,Any}()
    for length_c in lengths
        master_body = deepcopy(template)
        master_body.strength .= 0.0
        master_wake = _flat_wake(master_body, baseline.gamma, length_c, c;
            free_row_length_c)
        result = _frozen_solve(master_body, master_wake, solver, C, Uinf, rho,
            Sref, c, dt, baseline.probes, baseline.probe_velocity;
            label="Task 2 $case nominal L/c=$length_c")
        cl_delta = result.CL - baseline.CL
        free_wake_length_c = result.wake.nwakes[] * free_row_length_c
        total_wake_length_c = das_length_c + free_wake_length_c
        push!(rows, (;
            nominal_length_c=length_c,
            das_length_c,
            free_row_length_c,
            free_wake_length_c,
            total_wake_length_c,
            wake_rows=result.wake.nwakes[],
            linear_residual_abs=result.residual_abs,
            linear_residual_rel=result.residual_rel,
            CL=result.CL,
            CL_minus_task1=cl_delta,
            relative_CL_difference=_relative_difference(result.CL, baseline.CL),
            transition_gamma_min=minimum(result.gamma),
            transition_gamma_max=maximum(result.gamma),
            transition_gamma_l2=norm(result.gamma),
            transition_free_mismatch_l2=result.mismatch_l2,
            transition_free_mismatch_rel=result.mismatch_rel,
            relative_exterior_velocity_difference=result.probe_rel,
            frozen_state_preserved=true,
        ))
        checksums["L_c_$(Int(length_c))"] = _snapshot_checksums(result.snapshot)
        @printf("Task 2 %s nominal L/c=%g (total %.3g): rows=%d, CL=%.10f, delta=%+.6e, residual=%.3e\n",
            case, length_c, total_wake_length_c, result.wake.nwakes[], result.CL,
            cl_delta, result.residual_rel)
        last_result = result
    end
    tail_relative_change = abs(rows[end].CL - rows[end-1].CL) /
        max(abs(rows[end].CL), abs(rows[end-1].CL), eps(Float64))
    flat_converged = tail_relative_change <= 1e-3
    status = flat_converged ? "converged" : "nonconverged_at_64c"
    _write_rows_csv(joinpath(path, artifact_stem * "_length_convergence.csv"), rows)
    open(joinpath(path, artifact_stem * "_frozen_state.toml"), "w") do io
        TOML.print(io, Dict(
            "all_preserved" => true,
            "das_length_c" => das_length_c,
            "free_row_length_c" => free_row_length_c,
            "tail_relative_CL_change" => tail_relative_change,
            "status" => status,
            "cases" => checksums,
        ); sorted=true)
    end
    terminal_time = rows[end].total_wake_length_c * c / U
    pnl.write_vtk(joinpath(path, "dirichlet_" * artifact_stem * "_terminal_body"),
        last_result.body, 0, terminal_time; monitors=(last_result.pressure, last_result.force),
        overwrite=true)
    pnl.write_vtk(joinpath(path, "dirichlet_" * artifact_stem * "_terminal_wake"),
        last_result.wake, 0, terminal_time; overwrite=true)

    final = rows[end]
    _upsert_comparison(joinpath(path, "comparison.csv"), (;
        task=2, case, formulation,
        status, frozen_wake_step=-1,
        wake_length_c=final.total_wake_length_c,
        wake_rows=final.wake_rows, ncells=template.ncells,
        nshedding_edges=size(C, 1), linear_residual_abs=final.linear_residual_abs,
        linear_residual_rel=final.linear_residual_rel, kutta_constant_error=0.0,
        gamma_min=final.transition_gamma_min, gamma_max=final.transition_gamma_max,
        gamma_mean=sum(last_result.gamma) / length(last_result.gamma),
        gamma_l2=final.transition_gamma_l2,
        transition_free_mismatch_l2=final.transition_free_mismatch_l2,
        transition_free_mismatch_rel=final.transition_free_mismatch_rel,
        CL=final.CL, CL_minus_task1=final.CL_minus_task1,
        relative_CL_difference=final.relative_CL_difference,
        relative_exterior_velocity_difference=final.relative_exterior_velocity_difference,
    ))
    _write_task2_available_manifests(path; alpha_deg)
    return (; rows, result=last_result, converged=flat_converged,
        status, tail_relative_change, finite_setup=(template, solver, C), baseline)
end

function run_task2_flat_das005(; path=DIRICHLET_DATA_PATH, baseline=nothing,
        alpha_deg=DEFAULT_ALPHA_DEG)
    return run_task2_flat(; path, baseline, alpha_deg, das_length_c=0.05,
        free_row_length_c=0.5, case="flat_das005",
        artifact_stem="task2_flat_das005",
        formulation="finite_flat_das005_velocity_through_sources")
end

function _settling_metrics(times_c, CL; window_c=10.0)
    last_time = last(times_c)
    idx = findall(t -> t >= last_time - window_c - 10eps(last_time), times_c)
    length(idx) >= 2 || return (pp=Inf, drift=Inf, n=length(idx))
    tw = times_c[idx]
    yw = CL[idx]
    scale = max(abs(mean(yw)), eps(Float64))
    pp = (maximum(yw) - minimum(yw)) / scale
    centered_t = tw .- mean(tw)
    slope = dot(centered_t, yw .- mean(yw)) / max(dot(centered_t, centered_t), eps(Float64))
    drift = abs(slope) * (last(tw) - first(tw)) / scale
    return (; pp, drift, n=length(idx))
end

function run_task2_march(; path=DIRICHLET_DATA_PATH, baseline=nothing,
        finite_setup=nothing, alpha_deg=DEFAULT_ALPHA_DEG,
        das_length_c=0.05)
    mkpath(path)
    c = 1.0 * FT_TO_M
    b = 4.0 * c
    U = 330.2 * FT_TO_M
    Uinf = _uinf_from_alpha(U, alpha_deg)
    rho = 1.225
    dt = 0.5c / U
    Sref = b * c
    baseline = isnothing(baseline) ? _task2_baseline(path; alpha_deg) : baseline
    template, solver, C = isnothing(finite_setup) ?
        _finite_solver_template(c, b, Uinf; das_length_c) : finite_setup
    actual_das_length_c = norm(view(template.Das[1], :, 1)) / c
    isapprox(actual_das_length_c, das_length_c; atol=10eps(Float64), rtol=0) ||
        error("finite march template has Das/c=$actual_das_length_c, expected $das_length_c")
    body = deepcopy(template)
    body.strength .= 0.0
    wake = pnl.PanelWake(body; nwakerows=160)
    frames = pnl.ReferenceFrame(body)
    backend = pnl.FastMultipoleBackend(
        expansion_order=10, multipole_acceptance=0.4, leaf_size=100)
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false,
        backend, file=false)
    force = pnl.ForceMonitor(161, 1;
        i_frame=-1, normalization=pnl.WingNormalization(rho, Sref, c),
        correct_kuttacondition=false, verbose=false, file=false)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)

    rows = NamedTuple[]
    times_c = Float64[]
    CL_history = Float64[]
    selected_step = -1
    selected_metrics = (pp=Inf, drift=Inf, n=0)
    for i_step in 0:160
        time_c = 0.5i_step
        # Production performs this immediately before wake influence.  Make
        # that prescribed geometry update the boundary of the immutable solve.
        pnl.update_TE!(wake, body)
        solve_body = deepcopy(body)
        solve_wake = deepcopy(wake)
        snapshot = _state_snapshot(solve_body, solve_wake)
        pnl._steady_aerodynamics!(solve_body, (solve_body,), (solve_wake,), frames,
            Uinf, solver;
            backend_wake=backend, backend_solve=backend, backend_system=backend,
            update_trailing_edges=false)
        _assert_frozen_state(snapshot, solve_body, solve_wake;
            label="Task 2 march step $i_step")
        body = solve_body
        wake = solve_wake
        context = pnl.MonitorContext()
        pnl._run_monitor!(pressure, context, (body,), (wake,), frames, Uinf,
            i_step, dt, i_step * dt)
        pnl._run_monitor!(force, context, (body,), (wake,), frames, Uinf,
            i_step, dt, i_step * dt)
        Lhat, _, _ = _lift_drag_span_directions(Uinf)
        CL = dot(SVector{3,Float64}(force.force[:, i_step + 1]), Lhat)
        mu = copy(view(body.strength, :, 2))
        gamma = C * mu
        residual = _lu_matrix_product(solver.Glu, mu) - solver.rhs
        residual_abs = norm(residual)
        residual_rel = residual_abs / max(norm(solver.rhs), eps(Float64))
        all(isfinite, (CL, residual_abs, residual_rel)) ||
            error("non-finite Task 2 march result at step $i_step")
        residual_rel <= 1e-10 || error("Task 2 march residual exceeded tolerance")
        push!(times_c, time_c)
        push!(CL_history, CL)
        metrics = time_c >= 10 ? _settling_metrics(times_c, CL_history) :
            (pp=Inf, drift=Inf, n=length(times_c))
        checkpoint = time_c in (40.0, 60.0, 80.0)
        settled = checkpoint && metrics.pp <= 0.005 && metrics.drift <= 0.005
        first_free = wake.nwakes[] > 0 ? vec(view(wake.strength[1], 1, 1, :)) :
            fill(NaN, length(gamma))
        mismatch_l2 = wake.nwakes[] > 0 ? norm(gamma - first_free) : NaN
        mismatch_rel = wake.nwakes[] > 0 ? mismatch_l2 / max(norm(first_free), eps(Float64)) : NaN
        push!(rows, (;
            step=i_step, time_c, wake_rows=wake.nwakes[], CL,
            CL_minus_task1=CL - baseline.CL,
            relative_CL_difference=_relative_difference(CL, baseline.CL),
            linear_residual_abs=residual_abs,
            linear_residual_rel=residual_rel,
            transition_gamma_l2=norm(gamma),
            transition_free_mismatch_l2=mismatch_l2,
            transition_free_mismatch_rel=mismatch_rel,
            settling_window_samples=metrics.n,
            relative_CL_peak_to_peak=metrics.pp,
            normalized_linear_drift=metrics.drift,
            checkpoint, settled,
            frozen_during_solve=true,
        ))
        if checkpoint
            @printf("Task 2 march tU/c=%g: rows=%d, CL=%.10f, pp=%.4g, drift=%.4g, settled=%s\n",
                time_c, wake.nwakes[], CL, metrics.pp, metrics.drift, string(settled))
        end
        if settled || i_step == 160
            selected_step = i_step
            selected_metrics = metrics
            break
        end
        pnl.propagate!(wake, dt; step=i_step, frames)
        pnl.propagate_kinematics!((body,), frames, dt)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        pnl.shed_wake!(wake, body)
    end

    settled = selected_metrics.pp <= 0.005 && selected_metrics.drift <= 0.005
    status = settled ? "settled" : "nonconverged_at_80c"
    frozen = _frozen_solve(body, wake, solver, C, Uinf, rho, Sref, c, dt,
        baseline.probes, baseline.probe_velocity;
        label="Task 2 production frozen step $selected_step", backend)
    _write_rows_csv(joinpath(path, "task2_march_convergence.csv"), rows)
    open(joinpath(path, "task2_march_frozen_state.toml"), "w") do io
        TOML.print(io, Dict(
            "all_marching_solves_preserved_geometry_and_wake" => true,
            "frozen_reporting_solve_preserved" => true,
            "selected_step" => selected_step,
            "selected_time_c" => 0.5selected_step,
            "status" => status,
            "relative_CL_peak_to_peak" => selected_metrics.pp,
            "normalized_linear_drift" => selected_metrics.drift,
            "snapshot" => _snapshot_checksums(frozen.snapshot),
        ); sorted=true)
    end
    terminal_time = 0.5selected_step * c / U
    pnl.write_vtk(joinpath(path, "dirichlet_task2_march_terminal_body"),
        frozen.body, selected_step, terminal_time;
        monitors=(frozen.pressure, frozen.force), overwrite=true)
    pnl.write_vtk(joinpath(path, "dirichlet_task2_march_terminal_wake"),
        frozen.wake, selected_step, terminal_time; overwrite=true)

    _upsert_comparison(joinpath(path, "comparison.csv"), (;
        task=2, case="march", formulation="finite_rolled_velocity_through_sources",
        status, frozen_wake_step=selected_step, wake_length_c=0.5selected_step,
        wake_rows=frozen.wake.nwakes[], ncells=body.ncells,
        nshedding_edges=size(C, 1), linear_residual_abs=frozen.residual_abs,
        linear_residual_rel=frozen.residual_rel, kutta_constant_error=0.0,
        gamma_min=minimum(frozen.gamma), gamma_max=maximum(frozen.gamma),
        gamma_mean=sum(frozen.gamma) / length(frozen.gamma), gamma_l2=norm(frozen.gamma),
        transition_free_mismatch_l2=frozen.mismatch_l2,
        transition_free_mismatch_rel=frozen.mismatch_rel,
        CL=frozen.CL, CL_minus_task1=frozen.CL - baseline.CL,
        relative_CL_difference=_relative_difference(frozen.CL, baseline.CL),
        relative_exterior_velocity_difference=frozen.probe_rel,
    ))
    _write_task2_available_manifests(path; alpha_deg,
        march_das_length_c=das_length_c)
    return (; rows, result=frozen, settled, status, selected_step,
        metrics=selected_metrics, finite_setup=(template, solver, C), baseline)
end

function _write_task2_available_manifests(path; alpha_deg=DEFAULT_ALPHA_DEG,
        march_das_length_c=0.05)
    stored = _read_comparison_rows(joinpath(path, "comparison.csv"))
    flat = get(stored, "finite_flat_velocity_through_sources|flat", nothing)
    flat_das005 = get(stored,
        "finite_flat_das005_velocity_through_sources|flat_das005", nothing)
    march = get(stored, "finite_rolled_velocity_through_sources|march", nothing)
    config = Dict(
        "task" => 2,
        "formulation" => "finite_velocity_through_sources",
        "geometry" => Dict(
            "airfoil" => "NACA 0015", "chord_m" => 1.0 * FT_TO_M,
            "span_m" => 4.0 * FT_TO_M, "aspect_ratio" => 4.0,
            "n_airfoil" => 161, "n_span" => 13, "n_endcap" => 9,
        ),
        "flow" => Dict("Uinf_ft_s" => 330.2, "alpha_deg" => alpha_deg,
            "rho_kg_m3" => 1.225, "constant_incidence" => true,
            "pitching" => false),
        "transition" => Dict(
            "flat_primary_length_c" => 0.5,
            "march_length_c" => march_das_length_c,
            "c_per_dt" => 0.5,
        ),
        "flat" => Dict(
            "lengths_c" => [1, 2, 4, 8, 16, 32, 64],
            "das_length_c" => 0.5,
            "free_row_length_c" => 0.5,
            "prescribed_strength" => "fresh_task1_gamma",
            "tail_relative_tolerance" => 1e-3,
            "status" => isnothing(flat) ? "not_run" : get(flat, :status, ""),
        ),
        "flat_das005" => Dict(
            "nominal_lengths_c" => [1, 2, 4, 8, 16, 32, 64],
            "das_length_c" => 0.05,
            "free_row_length_c" => 0.5,
            "free_wake_lengths_c" => [0.5, 1.5, 3.5, 7.5, 15.5, 31.5, 63.5],
            "total_wake_lengths_c" => [0.55, 1.55, 3.55, 7.55, 15.55, 31.55, 63.55],
            "prescribed_strength" => "fresh_task1_gamma",
            "tail_relative_tolerance" => 1e-3,
            "status" => isnothing(flat_das005) ? "not_run" :
                get(flat_das005, :status, ""),
        ),
        "march" => Dict(
            "checkpoints_c" => [40, 60, 80], "settling_window_c" => 10,
            "das_length_c" => march_das_length_c,
            "relative_peak_to_peak_tolerance" => 0.005,
            "normalized_drift_tolerance" => 0.005,
            "default_induced_velocity_convection" => true,
            "status" => isnothing(march) ? "not_run" : get(march, :status, ""),
        ),
        "backend" => Dict(
            "flat" => "DirectBackend",
            "march" => "FastMultipoleBackend(order=10, theta=0.4, leaf_size=100)",
            "residual" => "direct assembled matrix product",
        ),
        "julia_threads_required" => 2,
    )
    open(joinpath(path, "task2.config.toml"), "w") do io
        TOML.print(io, config; sorted=true)
    end
    metadata = Dict{String,Any}(
        "theory_file" => "docs/wake_solve_schemes.md",
        "theory_section" => "Discrete diagnostic construction: Task 2 finite velocity-through-sources wake",
        "comparison_csv" => "comparison.csv",
        "flat_convergence_csv" => "task2_flat_length_convergence.csv",
        "flat_das005_convergence_csv" => "task2_flat_das005_length_convergence.csv",
        "march_convergence_csv" => "task2_march_convergence.csv",
        "flat_frozen_state" => "task2_flat_frozen_state.toml",
        "flat_das005_frozen_state" => "task2_flat_das005_frozen_state.toml",
        "march_frozen_state" => "task2_march_frozen_state.toml",
        "terminal_vtk" => [
            "dirichlet_task2_flat_terminal_body.pvd",
            "dirichlet_task2_flat_terminal_wake.pvd",
            "dirichlet_task2_flat_das005_terminal_body.pvd",
            "dirichlet_task2_flat_das005_terminal_wake.pvd",
            "dirichlet_task2_march_terminal_body.pvd",
            "dirichlet_task2_march_terminal_wake.pvd",
        ],
        "alpha_deg" => alpha_deg,
        "march_das_length_c" => march_das_length_c,
    )
    !isnothing(flat) && (metadata["flat_terminal_CL"] = parse(Float64, flat[:CL]))
    !isnothing(flat_das005) &&
        (metadata["flat_das005_terminal_CL"] = parse(Float64, flat_das005[:CL]))
    !isnothing(march) && (metadata["march_terminal_CL"] = parse(Float64, march[:CL]))
    baseline = get(stored, "semi_infinite_attached_wake|baseline", nothing)
    !isnothing(baseline) && (metadata["task1_CL"] = parse(Float64, baseline[:CL]))
    open(joinpath(path, "task2.metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return nothing
end

function run_task2(; path=DIRICHLET_DATA_PATH, alpha_deg=DEFAULT_ALPHA_DEG)
    baseline = _task2_baseline(path; alpha_deg)
    flat_setup = _finite_solver_template(1.0 * FT_TO_M, 4.0 * FT_TO_M,
        _uinf_from_alpha(330.2 * FT_TO_M, alpha_deg))
    flat = run_task2_flat(; path, baseline, finite_setup=flat_setup, alpha_deg)
    flat_das005 = run_task2_flat_das005(; path, baseline, alpha_deg)
    march = run_task2_march(; path, baseline, alpha_deg, das_length_c=0.05)
    _write_task2_available_manifests(path; alpha_deg)
    return (; baseline, flat, flat_das005, march)
end

# -----------------------------------------------------------------------------
# Task 3: direct prescribed-wake potential, diagnostic only
# -----------------------------------------------------------------------------

function _task3_parameters(alpha_deg)
    c = 1.0 * FT_TO_M
    b = 4.0 * c
    U = 330.2 * FT_TO_M
    return (;
        c, b, U, Uinf=_uinf_from_alpha(U, alpha_deg), rho=1.225,
        dt=0.5c / U, Sref=b * c,
    )
end

function _task3_assert_panel_sources(wake; label)
    wake isa pnl.PanelWake || error("$label requires a PanelWake")
    wake.nwakes[] > 0 || error("$label has no active wake-panel rows")
    pnl.FastMultipole.has_vector_potential(wake) && error(
        "$label PanelWake unexpectedly exposes only vector potential")
    final_filament = pnl.FilamentWrapper(wake)
    nfinal = pnl.FastMultipole.get_n_bodies(final_filament)
    nfinal == 0 || error("$label has $nfinal active final filaments")
    sources = (wake,)
    all(!pnl.FastMultipole.has_vector_potential(source) for source in sources) ||
        error("$label contains a vector-potential-only active source")
    return (; panel_sources=length(sources), active_final_filaments=nfinal,
        include_final_filament_option=wake.include_final_filament)
end

function _task3_assert_matching_seed(reference, candidate; label)
    isequal(reference.body_nodes, candidate.body_nodes) ||
        error("$label body nodes do not match Task 2")
    isequal(reference.body_cells, candidate.body_cells) ||
        error("$label body cells do not match Task 2")
    all(isequal(a, b) for (a, b) in zip(reference.body_Das, candidate.body_Das)) ||
        error("$label attached geometry does not match Task 2")
    reference.wake_rows == candidate.wake_rows ||
        error("$label wake row count does not match Task 2")
    isequal(reference.wake_options, candidate.wake_options) ||
        error("$label wake options do not match Task 2")
    all(isequal(a, b) for (a, b) in zip(reference.wake_nodes, candidate.wake_nodes)) ||
        error("$label wake nodes do not match Task 2")
    all(isequal(a, b) for (a, b) in zip(reference.wake_strength,
        candidate.wake_strength)) || error("$label initial wake strengths do not match Task 2")
    _array_hash(reference.body_strength) == _array_hash(candidate.body_strength) ||
        error("$label initial body-strength hash does not match Task 2")
    return true
end

function _task3_sigma0!(body, Uinf)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    sigma0 = pnl.calc_bc_noflowthrough(repeat(reshape(collect(Uinf), 3, 1),
        1, body.ncells), body.normals)
    body.strength[:, 1] .= sigma0
    body.strength[:, 2] .= 0.0
    return sigma0
end

function _task3_accumulated_potential(body, sources)
    saved = copy(body.potential)
    try
        fill!(body.potential, 0.0)
        pnl.influence!((body,), sources, pnl.DirectBackend();
            scalar_potential=true, velocity=false)
        q = copy(body.potential)
        all(isfinite, q) || error("non-finite Task 3 scalar potential")
        return q
    finally
        body.potential .= saved
    end
end

function _task3_surface_velocity!(body, wake, Uinf)
    saved_potential = copy(body.potential)
    try
        pnl.calcfield_U!(body, Uinf, (wake,); backend=pnl.DirectBackend())
    finally
        body.potential .= saved_potential
    end
    all(isfinite, body.velocity) || error("non-finite Task 3 surface velocity")
    return body.velocity
end

function _task3_panel_exterior_velocity(body, wake, Uinf, points)
    probes = pnl.FastMultipole.ProbeSystem(length(points), Float64)
    for i in eachindex(points)
        probes.position[i] = points[i]
        probes.gradient[i] = zero(SVector{3,Float64})
    end
    pnl.influence!((probes,), (body, wake), pnl.DirectBackend();
        scalar_potential=false, velocity=true)
    velocities = zeros(Float64, 3, length(points))
    for i in eachindex(points)
        velocities[:, i] .= Uinf .+ probes.gradient[i]
    end
    all(isfinite, velocities) || error("non-finite Task 3 exterior velocity")
    return velocities
end

function _task3_direct_solve(master_body, master_wake, solver, C, Uinf, rho,
        Sref, c, dt, probes, task1_probe, task2_probe; label)
    body = deepcopy(master_body)
    wake = deepcopy(master_wake)
    source_assertions = _task3_assert_panel_sources(wake; label)
    snapshot = _state_snapshot(body, wake)
    sigma0 = _task3_sigma0!(body, Uinf)
    source_potential = _task3_accumulated_potential(body, (body,))
    q_f = _task3_accumulated_potential(body, (wake,))
    rhs = -source_potential - q_f
    solver.rhs .= rhs
    ldiv!(view(body.strength, :, 2), solver.Glu, solver.rhs)
    _assert_frozen_state(snapshot, body, wake; label)
    isequal(view(body.strength, :, 1), sigma0) || error("$label changed sigma0")

    mu = copy(view(body.strength, :, 2))
    gamma = C * mu
    residual = _lu_matrix_product(solver.Glu, mu) - rhs
    residual_abs = norm(residual)
    residual_rel = residual_abs / max(norm(rhs), eps(Float64))
    all(isfinite, mu) && all(isfinite, gamma) && all(isfinite, q_f) &&
        all(isfinite, (residual_abs, residual_rel)) ||
        error("$label produced a non-finite solve result")
    residual_rel <= 1e-10 || error(
        "$label relative linear residual $residual_rel exceeds 1e-10")

    _task3_surface_velocity!(body, wake, Uinf)
    frames = pnl.ReferenceFrame(body)
    CL, pressure, force = _monitor_lift!(body, wake, frames, Uinf, rho, Sref, c,
        dt, pnl.DirectBackend())
    first_free = vec(copy(view(wake.strength[1], 1, 1, :)))
    mismatch_l2 = norm(gamma - first_free)
    mismatch_rel = mismatch_l2 / max(norm(first_free), eps(Float64))
    velocity = _task3_panel_exterior_velocity(body, wake, Uinf, probes)
    probe_rel_task1 = norm(velocity - task1_probe) /
        max(norm(task1_probe), eps(Float64))
    probe_rel_task2 = norm(velocity - task2_probe) /
        max(norm(task2_probe), eps(Float64))
    Cq_f = C * q_f
    all(isfinite, Cq_f) && all(isfinite, (CL, mismatch_l2, mismatch_rel,
        probe_rel_task1, probe_rel_task2)) ||
        error("$label produced a non-finite diagnostic")
    _assert_frozen_state(snapshot, body, wake; label="$label complete diagnostic")
    isequal(view(body.strength, :, 1), sigma0) ||
        error("$label changed sigma0 during post-processing")
    return (;
        body, wake, snapshot, source_assertions, sigma0, source_potential, q_f,
        Cq_f, rhs, pressure, force, mu, gamma, CL, residual_abs, residual_rel,
        mismatch_l2, mismatch_rel, probe_velocity=velocity,
        probe_rel_task1, probe_rel_task2,
    )
end

function _task3_active_strengths(wake)
    rows = wake.nwakes[]
    blocks = Matrix{Float64}[]
    for strength in wake.strength
        size(strength, 1) == 1 || error("Task 3 requires scalar wake strengths")
        push!(blocks, Matrix(view(strength, 1, 1:rows, :)))
    end
    return hcat(blocks...)
end

function _task3_set_active_strengths!(wake, eta)
    rows = wake.nwakes[]
    offset = 0
    for strength in wake.strength
        nedge = size(strength, 3)
        size(eta, ndims(eta)) >= offset + nedge ||
            error("invalid Task 3 strip strength")
        for edge in 1:nedge
            target = view(strength, 1, 1:rows, edge)
            if eta isa AbstractMatrix
                size(eta, 1) == rows || error("invalid Task 3 wake-row count")
                target .= view(eta, :, offset + edge)
            else
                target .= eta[offset + edge]
            end
        end
        offset += nedge
    end
    offset == size(eta, ndims(eta)) || error("unused Task 3 strip strengths")
    return wake
end

function _task3_inactive_strength_hashes(wake)
    rows = wake.nwakes[]
    return [_array_hash(copy(view(strength, :, rows+1:size(strength, 2), :)))
        for strength in wake.strength]
end

function _task3_strip_potential_matrix(body, wake, C)
    saved_strength = [copy(strength) for strength in wake.strength]
    P = zeros(Float64, body.ncells, size(C, 1))
    try
        zeros_eta = zeros(size(C, 1))
        _task3_set_active_strengths!(wake, zeros_eta)
        for edge in axes(C, 1)
            eta = copy(zeros_eta)
            eta[edge] = 1.0
            _task3_set_active_strengths!(wake, eta)
            P[:, edge] .= _task3_accumulated_potential(body, (wake,))
            _task3_set_active_strengths!(wake, zeros_eta)
        end
    finally
        for (strength, saved) in zip(wake.strength, saved_strength)
            strength .= saved
        end
    end
    all(isequal(a, b) for (a, b) in zip(wake.strength, saved_strength)) ||
        error("strip-potential assembly changed wake strengths")
    all(isfinite, P) || error("non-finite strip-potential matrix")
    return P
end

function _task3_fixed_geometry_iteration(master_body, master_wake, solver, C,
        Uinf, rho, Sref, c, dt, probes, task1_probe, task2_probe; label,
        max_iterations=200)
    body = deepcopy(master_body)
    wake = deepcopy(master_wake)
    geometry_snapshot = _state_snapshot(body, wake)
    inactive_hashes = _task3_inactive_strength_hashes(wake)
    rows = NamedTuple[]
    omega = 1.0
    previous_defect = Inf
    previous_CL = NaN
    insufficient_decrease = 0
    consecutive = 0
    converged = false
    result = nothing
    for iteration in 1:max_iterations
        before = _state_snapshot(body, wake)
        active_before = _task3_active_strengths(wake)
        result = _task3_direct_solve(body, wake, solver, C, Uinf, rho, Sref, c,
            dt, probes, task1_probe, task2_probe; label="$label iteration $iteration")
        _assert_frozen_state(before, result.body, result.wake;
            label="$label inner solve $iteration")
        defect = maximum(abs.(active_before .- reshape(result.gamma, 1, :)))
        lift_change = iteration == 1 ? Inf : abs(result.CL - previous_CL) /
            max(abs(result.CL), abs(previous_CL), eps(Float64))
        if defect <= 1e-8 && lift_change <= 1e-8
            consecutive += 1
        else
            consecutive = 0
        end
        push!(rows, (;
            iteration, omega, max_rowwise_strength_defect=defect,
            relative_lift_change=lift_change, consecutive_converged=consecutive,
            linear_residual_abs=result.residual_abs,
            linear_residual_rel=result.residual_rel, CL=result.CL,
            transition_gamma_min=minimum(result.gamma),
            transition_gamma_max=maximum(result.gamma),
            transition_gamma_l2=norm(result.gamma),
            transition_free_mismatch_l2=result.mismatch_l2,
            transition_free_mismatch_rel=result.mismatch_rel,
            q_f_l2=norm(result.q_f), Cq_f_l2=norm(result.Cq_f),
            frozen_inner_solve=true,
        ))
        @printf("%s iteration %d: omega=%.4g, defect=%.3e, dCL=%.3e, streak=%d\n",
            label, iteration, omega, defect, lift_change, consecutive)
        flush(stdout)
        if consecutive >= 3
            body, wake = result.body, result.wake
            converged = true
            break
        end
        insufficient_decrease = isfinite(previous_defect) &&
            defect > 0.95previous_defect ? insufficient_decrease + 1 : 0
        if insufficient_decrease >= 2
            omega = max(omega / 2, 1 / 16)
            insufficient_decrease = 0
        end
        eta = (1 - omega) .* active_before .+
            omega .* reshape(result.gamma, 1, :)
        body, wake = result.body, result.wake
        _task3_set_active_strengths!(wake, eta)
        isequal(body.nodes, geometry_snapshot.body_nodes) ||
            error("$label changed body geometry between iterations")
        all(isequal(a, b) for (a, b) in zip(wake.nodes,
            geometry_snapshot.wake_nodes)) || error("$label changed wake geometry")
        wake.nwakes[] == geometry_snapshot.wake_rows || error("$label changed row count")
        isequal(_task3_inactive_strength_hashes(wake), inactive_hashes) ||
            error("$label changed inactive wake-strength storage")
        previous_defect = defect
        previous_CL = result.CL
    end
    status = converged ? "converged" : "nonconverged_at_200"

    P = _task3_strip_potential_matrix(result.body, result.wake, C)
    rhs0 = -result.source_potential
    # Woodbury evaluates the stated augmented oracle without materializing a
    # second dense 6688-by-6688 matrix. G_delta is already LU-factored.
    x0 = solver.Glu \ rhs0
    Z = solver.Glu \ P
    small = C * Z
    for i in axes(small, 1)
        small[i, i] += 1.0
    end
    mu_oracle = x0 - Z * (small \ (C * x0))
    oracle_residual = _lu_matrix_product(solver.Glu, mu_oracle) +
        P * (C * mu_oracle) - rhs0
    oracle_residual_abs = norm(oracle_residual)
    oracle_residual_rel = oracle_residual_abs / max(norm(rhs0), eps(Float64))
    isfinite(oracle_residual_rel) || error("$label oracle residual is non-finite")
    oracle_residual_rel <= 1e-10 || error(
        "$label oracle residual $oracle_residual_rel exceeds 1e-10")
    oracle_body = deepcopy(result.body)
    oracle_wake = deepcopy(result.wake)
    oracle_body.strength[:, 1] .= result.sigma0
    oracle_body.strength[:, 2] .= mu_oracle
    _task3_set_active_strengths!(oracle_wake, C * mu_oracle)
    _task3_surface_velocity!(oracle_body, oracle_wake, Uinf)
    oracle_CL, oracle_pressure, oracle_force = _monitor_lift!(oracle_body,
        oracle_wake, pnl.ReferenceFrame(oracle_body), Uinf, rho, Sref, c, dt,
        pnl.DirectBackend())
    relative_oracle_mu_difference = norm(result.mu - mu_oracle) /
        max(norm(mu_oracle), eps(Float64))
    relative_oracle_CL_difference = abs(result.CL - oracle_CL) /
        max(abs(oracle_CL), eps(Float64))
    if converged
        relative_oracle_mu_difference <= 1e-8 || error(
            "$label iterated/oracle mu difference exceeds 1e-8")
        relative_oracle_CL_difference <= 1e-8 || error(
            "$label iterated/oracle lift difference exceeds 1e-8")
    end
    oracle = (; P, body=oracle_body, wake=oracle_wake, pressure=oracle_pressure,
        force=oracle_force, mu=mu_oracle, CL=oracle_CL, oracle_residual_abs,
        oracle_residual_rel, relative_oracle_mu_difference,
        relative_oracle_CL_difference)
    return (; rows, result, converged, status, omega, oracle,
        inactive_hashes, geometry_snapshot)
end

function _task3_write_terminal_vtk(path, artifact_stem, single, iterated, time,
        i_step=0)
    for (suffix, result) in (("single_shot", single), ("iterated", iterated))
        pnl.write_vtk(joinpath(path, "dirichlet_" * artifact_stem * "_" * suffix *
            "_body"), result.body, i_step, time;
            monitors=(result.pressure, result.force), overwrite=true)
        pnl.write_vtk(joinpath(path, "dirichlet_" * artifact_stem * "_" * suffix *
            "_wake"), result.wake, i_step, time; overwrite=true)
    end
end

function _task3_comparison!(path, case, formulation, single, iteration, baseline,
        task2_result; frozen_step=-1, wake_length_c)
    _upsert_comparison(joinpath(path, "comparison.csv"), (;
        task=3, case, formulation, status=iteration.status,
        frozen_wake_step=frozen_step, wake_length_c,
        wake_rows=single.wake.nwakes[], ncells=single.body.ncells,
        nshedding_edges=length(single.gamma),
        linear_residual_abs=single.residual_abs,
        linear_residual_rel=single.residual_rel, kutta_constant_error=0.0,
        gamma_min=minimum(single.gamma), gamma_max=maximum(single.gamma),
        gamma_mean=sum(single.gamma) / length(single.gamma),
        gamma_l2=norm(single.gamma),
        transition_free_mismatch_l2=single.mismatch_l2,
        transition_free_mismatch_rel=single.mismatch_rel,
        CL=single.CL, CL_minus_task1=single.CL - baseline.CL,
        relative_CL_difference=_relative_difference(single.CL, baseline.CL),
        relative_exterior_velocity_difference=single.probe_rel_task1,
        CL_minus_task2=single.CL - task2_result.CL,
        relative_CL_difference_task2=_relative_difference(single.CL,
            task2_result.CL),
        relative_exterior_velocity_difference_task2=single.probe_rel_task2,
        q_f_l2=norm(single.q_f), Cq_f_l2=norm(single.Cq_f),
        fixed_point_iterations=length(iteration.rows),
        oracle_linear_residual_rel=iteration.oracle.oracle_residual_rel,
        relative_oracle_mu_difference=iteration.oracle.relative_oracle_mu_difference,
        relative_oracle_CL_difference=iteration.oracle.relative_oracle_CL_difference,
    ))
end

function _task3_flat_row(length_c, das_length_c, free_row_length_c, task2, task3,
        baseline)
    free_length = task3.wake.nwakes[] * free_row_length_c
    return (;
        nominal_length_c=length_c, das_length_c, free_row_length_c,
        free_wake_length_c=free_length,
        total_wake_length_c=das_length_c + free_length,
        wake_rows=task3.wake.nwakes[], q_f_l2=norm(task3.q_f),
        Cq_f_l2=norm(task3.Cq_f), linear_residual_abs=task3.residual_abs,
        linear_residual_rel=task3.residual_rel, CL=task3.CL,
        CL_minus_task1=task3.CL - baseline.CL,
        relative_CL_difference=_relative_difference(task3.CL, baseline.CL),
        CL_minus_task2=task3.CL - task2.CL,
        relative_CL_difference_task2=_relative_difference(task3.CL, task2.CL),
        transition_gamma_min=minimum(task3.gamma),
        transition_gamma_max=maximum(task3.gamma),
        transition_gamma_l2=norm(task3.gamma),
        transition_free_mismatch_l2=task3.mismatch_l2,
        transition_free_mismatch_rel=task3.mismatch_rel,
        relative_exterior_velocity_difference_task1=task3.probe_rel_task1,
        relative_exterior_velocity_difference_task2=task3.probe_rel_task2,
        task2_seed_matched=true, frozen_state_preserved=true,
    )
end

function run_task3_flat(; path=DIRICHLET_DATA_PATH, baseline=nothing,
        alpha_deg=DEFAULT_ALPHA_DEG, das_length_c=0.5,
        free_row_length_c=0.5, case="flat", artifact_stem="task3_flat",
        formulation="finite_flat_direct_fixed_wake_potential")
    mkpath(path)
    p = _task3_parameters(alpha_deg)
    baseline = isnothing(baseline) ? _task2_baseline(path; alpha_deg) : baseline
    template, solver, C = _finite_solver_template(p.c, p.b, p.Uinf; das_length_c)
    norm(C * ones(template.ncells), Inf) == 0.0 || error("invalid Task 3 Kutta map")
    rows = NamedTuple[]
    checksums = Dict{String,Any}()
    terminal = nothing
    terminal_task2 = nothing
    for length_c in (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
        master_body = deepcopy(template)
        master_body.strength .= 0.0
        master_wake = _flat_wake(master_body, baseline.gamma, length_c, p.c;
            free_row_length_c)
        task2 = _frozen_solve(master_body, master_wake, solver, C, p.Uinf, p.rho,
            p.Sref, p.c, p.dt, baseline.probes, baseline.probe_velocity;
            label="Task 3 seed Task 2 $case L/c=$length_c")
        candidate_body, candidate_wake = deepcopy(master_body), deepcopy(master_wake)
        _task3_assert_matching_seed(_state_snapshot(master_body, master_wake),
            _state_snapshot(candidate_body, candidate_wake);
            label="Task 3 $case L/c=$length_c")
        task3 = _task3_direct_solve(candidate_body, candidate_wake, solver, C,
            p.Uinf, p.rho, p.Sref, p.c, p.dt, baseline.probes,
            baseline.probe_velocity, task2.probe_velocity;
            label="Task 3 $case L/c=$length_c")
        push!(rows, _task3_flat_row(length_c, das_length_c, free_row_length_c,
            task2, task3, baseline))
        checksums["L_c_$(Int(length_c))"] = Dict(
            "task2_seed" => _snapshot_checksums(_state_snapshot(master_body,
                master_wake)),
            "task3_single_shot" => _snapshot_checksums(task3.snapshot),
            "active_final_filaments" => task3.source_assertions.active_final_filaments,
            "panel_sources_only" => true,
        )
        @printf("Task 3 %s L/c=%g: CL=%.10f, Task2 delta=%+.6e, residual=%.3e\n",
            case, length_c, task3.CL, task3.CL - task2.CL, task3.residual_rel)
        terminal, terminal_task2 = task3, task2
    end
    _write_rows_csv(joinpath(path, artifact_stem * "_single_shot.csv"), rows)
    iteration = _task3_fixed_geometry_iteration(terminal.body, terminal.wake,
        solver, C, p.Uinf, p.rho, p.Sref, p.c, p.dt, baseline.probes,
        baseline.probe_velocity, terminal_task2.probe_velocity;
        label="Task 3 $case terminal")
    _write_rows_csv(joinpath(path, artifact_stem * "_iteration_history.csv"),
        iteration.rows)
    terminal_time = rows[end].total_wake_length_c * p.c / p.U
    _task3_write_terminal_vtk(path, artifact_stem, terminal, iteration.result,
        terminal_time)
    checksums["terminal_iteration"] = Dict(
        "status" => iteration.status,
        "iterations" => length(iteration.rows),
        "inactive_strength_hashes" => iteration.inactive_hashes,
        "oracle_residual_rel" => iteration.oracle.oracle_residual_rel,
        "relative_oracle_mu_difference" =>
            iteration.oracle.relative_oracle_mu_difference,
        "relative_oracle_CL_difference" =>
            iteration.oracle.relative_oracle_CL_difference,
    )
    _task3_comparison!(path, case, formulation, terminal, iteration, baseline,
        terminal_task2; wake_length_c=rows[end].total_wake_length_c)
    _write_task3_available_manifests(path; alpha_deg,
        invariants_update=Dict(artifact_stem => checksums))
    return (; rows, result=terminal, iteration, task2=terminal_task2, baseline,
        finite_setup=(template, solver, C))
end

function run_task3_flat_das005(; path=DIRICHLET_DATA_PATH, baseline=nothing,
        alpha_deg=DEFAULT_ALPHA_DEG)
    return run_task3_flat(; path, baseline, alpha_deg, das_length_c=0.05,
        free_row_length_c=0.5, case="flat_das005",
        artifact_stem="task3_flat_das005",
        formulation="finite_flat_das005_direct_fixed_wake_potential")
end

function run_task3_march(; path=DIRICHLET_DATA_PATH, baseline=nothing,
        alpha_deg=DEFAULT_ALPHA_DEG)
    mkpath(path)
    p = _task3_parameters(alpha_deg)
    baseline = isnothing(baseline) ? _task2_baseline(path; alpha_deg) : baseline
    task2_march = run_task2_march(; path, baseline, alpha_deg,
        das_length_c=0.05)
    task2 = task2_march.result
    template, solver, C = task2_march.finite_setup
    norm(C * ones(template.ncells), Inf) == 0.0 || error("invalid Task 3 Kutta map")
    master_body, master_wake = deepcopy(task2.body), deepcopy(task2.wake)
    reference = _state_snapshot(task2.body, task2.wake)
    candidate = _state_snapshot(master_body, master_wake)
    _task3_assert_matching_seed(reference, candidate;
        label="Task 3 march step $(task2_march.selected_step)")
    single = _task3_direct_solve(master_body, master_wake, solver, C, p.Uinf,
        p.rho, p.Sref, p.c, p.dt, baseline.probes, baseline.probe_velocity,
        task2.probe_velocity;
        label="Task 3 march step $(task2_march.selected_step)")
    row = (;
        selected_step=task2_march.selected_step,
        selected_time_c=0.5task2_march.selected_step,
        wake_rows=single.wake.nwakes[], q_f_l2=norm(single.q_f),
        Cq_f_l2=norm(single.Cq_f), linear_residual_abs=single.residual_abs,
        linear_residual_rel=single.residual_rel, CL=single.CL,
        CL_minus_task1=single.CL - baseline.CL,
        relative_CL_difference=_relative_difference(single.CL, baseline.CL),
        CL_minus_task2=single.CL - task2.CL,
        relative_CL_difference_task2=_relative_difference(single.CL, task2.CL),
        transition_gamma_min=minimum(single.gamma),
        transition_gamma_max=maximum(single.gamma),
        transition_gamma_l2=norm(single.gamma),
        transition_free_mismatch_l2=single.mismatch_l2,
        transition_free_mismatch_rel=single.mismatch_rel,
        relative_exterior_velocity_difference_task1=single.probe_rel_task1,
        relative_exterior_velocity_difference_task2=single.probe_rel_task2,
        task2_seed_matched=true, frozen_state_preserved=true,
    )
    _write_rows_csv(joinpath(path, "task3_march_single_shot.csv"), [row])
    iteration = _task3_fixed_geometry_iteration(single.body, single.wake, solver,
        C, p.Uinf, p.rho, p.Sref, p.c, p.dt, baseline.probes,
        baseline.probe_velocity, task2.probe_velocity; label="Task 3 march terminal")
    _write_rows_csv(joinpath(path, "task3_march_iteration_history.csv"),
        iteration.rows)
    time = 0.5task2_march.selected_step * p.c / p.U
    _task3_write_terminal_vtk(path, "task3_march", single, iteration.result,
        time, task2_march.selected_step)
    invariants = Dict(
        "task2_seed" => _snapshot_checksums(reference),
        "task3_single_shot" => _snapshot_checksums(single.snapshot),
        "selected_step" => task2_march.selected_step,
        "task2_status" => task2_march.status,
        "active_final_filaments" => single.source_assertions.active_final_filaments,
        "panel_sources_only" => true,
        "terminal_iteration" => Dict(
            "status" => iteration.status, "iterations" => length(iteration.rows),
            "inactive_strength_hashes" => iteration.inactive_hashes,
            "oracle_residual_rel" => iteration.oracle.oracle_residual_rel,
            "relative_oracle_mu_difference" =>
                iteration.oracle.relative_oracle_mu_difference,
            "relative_oracle_CL_difference" =>
                iteration.oracle.relative_oracle_CL_difference,
        ),
    )
    _task3_comparison!(path, "march",
        "finite_rolled_direct_fixed_wake_potential", single, iteration, baseline,
        task2; frozen_step=task2_march.selected_step,
        wake_length_c=0.5task2_march.selected_step)
    _write_task3_available_manifests(path; alpha_deg,
        invariants_update=Dict("task3_march" => invariants))
    return (; result=single, iteration, task2=task2_march, baseline,
        finite_setup=(template, solver, C))
end

function _read_toml_or_empty(path)
    isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()
end

function _write_task3_available_manifests(path; alpha_deg=DEFAULT_ALPHA_DEG,
        invariants_update=Dict{String,Any}())
    invariant_path = joinpath(path, "task3_invariants.toml")
    invariants = _read_toml_or_empty(invariant_path)
    merge!(invariants, invariants_update)
    invariants["all_single_shot_and_inner_solves_preserved"] = true
    invariants["C_times_one_inf"] = 0.0
    open(invariant_path, "w") do io
        TOML.print(io, invariants; sorted=true)
    end
    stored = _read_comparison_rows(joinpath(path, "comparison.csv"))
    forms = (
        ("flat", "finite_flat_direct_fixed_wake_potential|flat"),
        ("flat_das005", "finite_flat_das005_direct_fixed_wake_potential|flat_das005"),
        ("march", "finite_rolled_direct_fixed_wake_potential|march"),
    )
    statuses = Dict(name => (haskey(stored, key) ? stored[key][:status] : "not_run")
        for (name, key) in forms)
    config = Dict(
        "task" => 3, "formulation" => "direct_fixed_wake_potential",
        "geometry" => Dict("airfoil" => "NACA 0015", "chord_m" => 1FT_TO_M,
            "span_m" => 4FT_TO_M, "aspect_ratio" => 4.0,
            "n_airfoil" => 161, "n_span" => 13, "n_endcap" => 9),
        "flow" => Dict("Uinf_ft_s" => 330.2, "alpha_deg" => alpha_deg,
            "rho_kg_m3" => 1.225, "constant_incidence" => true),
        "single_shot" => Dict("flat_lengths_c" => [1,2,4,8,16,32,64],
            "flat_das_length_c" => 0.5, "short_das_length_c" => 0.05,
            "free_row_length_c" => 0.5, "statuses" => statuses),
        "iteration" => Dict("terminal_cases_only" => true,
            "initial_omega" => 1.0, "minimum_omega" => 0.0625,
            "defect_tolerance" => 1e-8, "relative_lift_tolerance" => 1e-8,
            "required_consecutive" => 3, "maximum_iterations" => 200,
            "poor_decrease_ratio" => 0.95, "poor_decrease_count" => 2),
        "linear_residual_tolerance" => 1e-10,
        "backend" => "DirectBackend", "julia_threads_required" => 2,
    )
    open(joinpath(path, "task3.config.toml"), "w") do io
        TOML.print(io, config; sorted=true)
    end
    terminal_vtk = String[]
    for stem in ("task3_flat", "task3_flat_das005", "task3_march"),
            route in ("single_shot", "iterated"), kind in ("body", "wake")
        push!(terminal_vtk, "dirichlet_$(stem)_$(route)_$(kind).pvd")
    end
    metadata = Dict(
        "theory_file" => "docs/wake_solve_schemes.md",
        "theory_section" => "Discrete diagnostic construction: Task 3 direct fixed-wake potential",
        "comparison_csv" => "comparison.csv",
        "single_shot_csv" => ["task3_flat_single_shot.csv",
            "task3_flat_das005_single_shot.csv", "task3_march_single_shot.csv"],
        "iteration_history_csv" => ["task3_flat_iteration_history.csv",
            "task3_flat_das005_iteration_history.csv",
            "task3_march_iteration_history.csv"],
        "invariants_toml" => "task3_invariants.toml",
        "terminal_vtk" => terminal_vtk, "alpha_deg" => alpha_deg,
    )
    open(joinpath(path, "task3.metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return nothing
end

function run_task3(; path=DIRICHLET_DATA_PATH, alpha_deg=DEFAULT_ALPHA_DEG)
    baseline = _task2_baseline(path; alpha_deg)
    flat = run_task3_flat(; path, baseline, alpha_deg)
    flat_das005 = run_task3_flat_das005(; path, baseline, alpha_deg)
    march = run_task3_march(; path, baseline, alpha_deg)
    _write_task3_available_manifests(path; alpha_deg)
    return (; baseline, flat, flat_das005, march)
end

# -----------------------------------------------------------------------------
# Task 4: Green-reconstruction formulation (velocity-only, explicit-potential)
#
# Runs `pnl.GreenReconstruction(gauge=:area_mean, green_solver=nothing)` on the
# same frozen wake states as Tasks 2/3, through the production
# `_steady_aerodynamics!` formulation hook.  The dense bordered-LU of the
# body-only Green operator (I−B) with the area-weighted gauge row is the
# authoritative full-resolution reference.  (I−B) is body-mesh-only (the
# attached wake is suppressed during its assembly), so a single Green state is
# reused across every wake state and both angles; only the mutable trace buffers
# are refreshed per solve.
# -----------------------------------------------------------------------------

const TASK4_FORMULATION = "GreenReconstruction(gauge=:area_mean, green_solver=nothing, recompute_interval=1)"

"Pin BLAS threads to a fixed, machine-capped value so recorded wall times are
comparable across cases; return the value actually set."
function _task4_pin_blas()
    nb = min(8, Sys.CPU_THREADS)
    BLAS.set_num_threads(nb)
    return BLAS.get_num_threads()
end

"Build the one-time body-only Green state (dense (I−B) assembly + bordered LU
with the area gauge row).  Body-mesh-only, hence reusable for every Task 4
case.  Validated once through the production `initialize_formulation` path."
function _task4_shared_green(template, solver; backend=pnl.DirectBackend())
    body = deepcopy(template)
    body.strength .= 0.0
    wake = _flat_wake(body, zeros(size(_build_kutta_map(body), 1)), 1.0,
        1.0 * FT_TO_M)
    formulation = pnl.GreenReconstruction(gauge=:area_mean, green_solver=nothing)
    state = pnl.initialize_formulation(formulation, (body,), (wake,), solver,
        backend, backend)
    return state.green
end

"Fresh Green-reconstruction runtime state that shares the immutable bordered-LU
`green` but owns fresh mutable trace buffers (`last_recompute=-1` forces a
recompute on first use, so no `q` leaks across cases)."
function _fresh_green_recon_state(body, green)
    TF = eltype(body.strength)
    N = body.ncells
    return pnl.GreenReconstructionState{TF, typeof(green)}(green,
        zeros(TF, 3, N), zeros(TF, N), zeros(TF, N), zeros(TF, N),
        zeros(TF, N), Ref(-1))
end

"Area-weighted zero-mean of a body-centroid field (matches the (I−B) gauge)."
function _remove_area_mean(x, areas)
    return x .- dot(areas, x) / sum(areas)
end

"""
One frozen Green-reconstruction solve on a copied body/wake.  Returns the solved
copies, the physical circulation γ=C·μE, formulation-specific residuals, and the
exterior/oracle diagnostics.  Gates: both formulation residuals ≤1e-10 relative,
finite fields, and geometry/wake immutability.
"""
function _task4_green_solve(master_body, master_wake, solver, green, C, Uinf, rho,
        Sref, c, dt, probes, task1_probe, task2_probe; label, i_step=0,
        compute_trace=true, backend=pnl.DirectBackend())
    body = deepcopy(master_body)
    wake = deepcopy(master_wake)
    snapshot = _state_snapshot(body, wake)
    state = _fresh_green_recon_state(body, green)
    formulation = pnl.GreenReconstruction(gauge=:area_mean, green_solver=nothing)
    frames = pnl.ReferenceFrame(body)

    pnl._steady_aerodynamics!(body, (body,), (wake,), frames, Uinf, solver;
        backend_wake=backend, backend_solve=backend, backend_system=backend,
        update_trailing_edges=false, formulation, formulation_state=state,
        i_step)
    _assert_frozen_state(snapshot, body, wake; label)

    mu = copy(view(body.strength, :, 2))          # explicit-potential μE
    sigma0_solved = copy(view(body.strength, :, 1))
    gamma = C * mu
    all(isfinite, body.strength) || error("$label produced non-finite strengths")

    # Formulation-specific residuals.  The dense bordered solve leaves
    # green.sol_b = [q; λ] and green.rhs_b = [Sσ; 0]; K·sol_b − rhs_b therefore
    # has first block (I−B)q + λa − Sσ and last entry the gauge row aᵀq.
    N = length(state.green.q)
    bordered = _lu_matrix_product(state.green.fact, state.green.sol_b) .-
        state.green.rhs_b
    green_res_abs = norm(view(bordered, 1:N))
    gauge_row_abs = abs(bordered[N + 1])
    green_res_rel = green_res_abs /
        max(norm(view(state.green.rhs_b, 1:N)), eps(Float64))
    # Explicit-potential residual G·μE + S·σ0 + q = G·μE − solver.rhs
    exp_residual = _lu_matrix_product(solver.Glu, mu) .- solver.rhs
    exp_res_abs = norm(exp_residual)
    exp_res_rel = exp_res_abs / max(norm(solver.rhs), eps(Float64))
    q = copy(state.green.q)
    all(isfinite, (green_res_rel, gauge_row_abs, exp_res_rel)) ||
        error("$label produced a non-finite formulation residual")
    green_res_rel <= 1e-10 || error(
        "$label Green bordered residual $green_res_rel exceeds 1e-10")
    exp_res_rel <= 1e-10 || error(
        "$label explicit-potential residual $exp_res_rel exceeds 1e-10")

    # Lift via steady PressureBernoulli + ForceMonitor (physical mode restored:
    # the formulation cleared any wake correction, γ = C·μE directly).
    CL, pressure, force = _monitor_lift!(body, wake, frames, Uinf, rho, Sref, c,
        dt, backend)
    isfinite(CL) || error("$label produced non-finite CL")

    first_free = vec(copy(view(wake.strength[1], 1, 1, :)))
    mismatch_l2 = norm(gamma - first_free)
    mismatch_rel = mismatch_l2 / max(norm(first_free), eps(Float64))

    velocity = _exterior_velocity(body, wake, Uinf, probes, backend)
    probe_rel_task1 = norm(velocity - task1_probe) /
        max(norm(task1_probe), eps(Float64))
    probe_rel_task2 = norm(velocity - task2_probe) /
        max(norm(task2_probe), eps(Float64))

    # Diagnostic (not a gate, and not part of the eligible solve path): the
    # direct free-wake trace q_f (the Task 3 oracle) at body centroids, compared
    # to the Green-reconstructed q after removing the area-weighted mean from
    # each so the comparison is gauge-independent.  Computed only for reported
    # results; inner iteration solves skip it (compute_trace=false).
    green_q_vs_trace_rel = NaN
    if compute_trace
        q_f = _task3_accumulated_potential(body, (wake,))
        areas = pnl.calc_areas(body.nodes, body.cells)
        dq = _remove_area_mean(q, areas) .- _remove_area_mean(q_f, areas)
        green_q_vs_trace_rel = norm(dq) /
            max(norm(_remove_area_mean(q_f, areas)), eps(Float64))
        isfinite(green_q_vs_trace_rel) ||
            error("$label produced a non-finite trace diagnostic")
    end
    _assert_frozen_state(snapshot, body, wake; label="$label complete diagnostic")

    all(isfinite, (mismatch_l2, mismatch_rel, probe_rel_task1, probe_rel_task2,
        norm(q))) ||
        error("$label produced a non-finite diagnostic")
    return (;
        body, wake, snapshot, pressure, force, mu, sigma0=sigma0_solved, gamma,
        CL, q, q_l2=norm(q), green_res_abs, green_res_rel, gauge_row_abs,
        exp_res_abs, exp_res_rel, residual_abs=exp_res_abs,
        residual_rel=exp_res_rel, mismatch_l2, mismatch_rel,
        probe_velocity=velocity, probe_rel_task1, probe_rel_task2,
        green_q_vs_trace_rel,
    )
end

"""
Fixed-geometry active-strength projection for Green reconstruction, mirroring
`_task3_fixed_geometry_iteration`: body and wake geometry stay frozen, only the
active free-wake strengths advance toward the transition circulation γ=C·μE with
adaptive relaxation.  Converged when the max active-row strength defect and the
relative lift change both stay ≤1e-8 for three consecutive iterations.
"""
function _task4_fixed_geometry_iteration(master_body, master_wake, solver, green,
        C, Uinf, rho, Sref, c, dt, probes, task1_probe, task2_probe; label,
        i_step0=0, max_iterations=200, backend=pnl.DirectBackend())
    body = deepcopy(master_body)
    wake = deepcopy(master_wake)
    geometry_snapshot = _state_snapshot(body, wake)
    inactive_hashes = _task3_inactive_strength_hashes(wake)
    rows = NamedTuple[]
    omega = 1.0
    previous_defect = Inf
    previous_CL = NaN
    insufficient_decrease = 0
    consecutive = 0
    converged = false
    result = nothing
    for iteration in 1:max_iterations
        before = _state_snapshot(body, wake)
        active_before = _task3_active_strengths(wake)
        result = _task4_green_solve(body, wake, solver, green, C, Uinf, rho,
            Sref, c, dt, probes, task1_probe, task2_probe;
            label="$label iteration $iteration", i_step=i_step0 + iteration,
            compute_trace=false, backend)
        _assert_frozen_state(before, result.body, result.wake;
            label="$label inner solve $iteration")
        defect = maximum(abs.(active_before .- reshape(result.gamma, 1, :)))
        lift_change = iteration == 1 ? Inf : abs(result.CL - previous_CL) /
            max(abs(result.CL), abs(previous_CL), eps(Float64))
        if defect <= 1e-8 && lift_change <= 1e-8
            consecutive += 1
        else
            consecutive = 0
        end
        push!(rows, (;
            iteration, omega, max_rowwise_strength_defect=defect,
            relative_lift_change=lift_change, consecutive_converged=consecutive,
            green_bordered_residual_rel=result.green_res_rel,
            gauge_row_abs=result.gauge_row_abs,
            explicit_potential_residual_rel=result.exp_res_rel,
            CL=result.CL, transition_gamma_min=minimum(result.gamma),
            transition_gamma_max=maximum(result.gamma),
            transition_gamma_l2=norm(result.gamma),
            transition_free_mismatch_l2=result.mismatch_l2,
            transition_free_mismatch_rel=result.mismatch_rel,
            q_l2=result.q_l2, frozen_inner_solve=true,
        ))
        @printf("%s iteration %d: omega=%.4g, defect=%.3e, dCL=%.3e, streak=%d\n",
            label, iteration, omega, defect, lift_change, consecutive)
        flush(stdout)
        if consecutive >= 3
            body, wake = result.body, result.wake
            converged = true
            break
        end
        insufficient_decrease = isfinite(previous_defect) &&
            defect > 0.95previous_defect ? insufficient_decrease + 1 : 0
        if insufficient_decrease >= 2
            omega = max(omega / 2, 1 / 16)
            insufficient_decrease = 0
        end
        eta = (1 - omega) .* active_before .+
            omega .* reshape(result.gamma, 1, :)
        body, wake = result.body, result.wake
        _task3_set_active_strengths!(wake, eta)
        isequal(body.nodes, geometry_snapshot.body_nodes) ||
            error("$label changed body geometry between iterations")
        all(isequal(a, b) for (a, b) in zip(wake.nodes,
            geometry_snapshot.wake_nodes)) || error("$label changed wake geometry")
        wake.nwakes[] == geometry_snapshot.wake_rows ||
            error("$label changed row count")
        isequal(_task3_inactive_strength_hashes(wake), inactive_hashes) ||
            error("$label changed inactive wake-strength storage")
        previous_defect = defect
        previous_CL = result.CL
    end
    status = converged ? "converged" : "nonconverged_at_$max_iterations"
    # One clean diagnostic solve on the converged active strengths so the
    # reported result carries the trace diagnostic (deterministic; reproduces
    # the last inner solve's μE/γ/CL).
    result = _task4_green_solve(body, wake, solver, green, C, Uinf, rho, Sref, c,
        dt, probes, task1_probe, task2_probe;
        label="$label converged diagnostic", i_step=i_step0 + length(rows) + 1,
        compute_trace=true, backend)
    return (; rows, result, converged, status, omega, inactive_hashes,
        geometry_snapshot)
end

"Two-tier closure fractions (both use absolute value; a value >1 is overshoot)."
function _task4_closures(cl_method, cl_semiinf, cl_old, cl_oracle)
    denom_semiinf = max(abs(cl_old - cl_semiinf), eps(Float64))
    denom_oracle = max(abs(cl_old - cl_oracle), eps(Float64))
    f_semiinf = 1 - abs(cl_method - cl_semiinf) / denom_semiinf
    f_oracle = 1 - abs(cl_method - cl_oracle) / denom_oracle
    return (; f_semiinf, f_oracle,
        signed_post_correction_diff=cl_method - cl_oracle)
end

"Upsert a Task 4 comparison row (unique formulation|case key, never overwrites
Tasks 1-3)."
function _task4_comparison_row!(path, case, formulation_id, status, frozen_step,
        wake_length_c, res, baseline_CL, task2_CL, t3_single_CL, t3_iter_CL,
        cl_old, cl_oracle, blas_threads)
    closures = _task4_closures(res.CL, baseline_CL, cl_old, cl_oracle)
    _upsert_comparison(joinpath(path, "comparison.csv"), (;
        task=4, case, formulation=formulation_id, status,
        frozen_wake_step=frozen_step, wake_length_c,
        wake_rows=res.wake.nwakes[], ncells=res.body.ncells,
        nshedding_edges=length(res.gamma), kutta_constant_error=0.0,
        linear_residual_abs=res.exp_res_abs, linear_residual_rel=res.exp_res_rel,
        gamma_min=minimum(res.gamma), gamma_max=maximum(res.gamma),
        gamma_mean=sum(res.gamma) / length(res.gamma), gamma_l2=norm(res.gamma),
        transition_free_mismatch_l2=res.mismatch_l2,
        transition_free_mismatch_rel=res.mismatch_rel,
        CL=res.CL, CL_minus_task1=res.CL - baseline_CL,
        relative_CL_difference=_relative_difference(res.CL, baseline_CL),
        relative_exterior_velocity_difference=res.probe_rel_task1,
        CL_minus_task2=res.CL - task2_CL,
        relative_CL_difference_task2=_relative_difference(res.CL, task2_CL),
        relative_exterior_velocity_difference_task2=res.probe_rel_task2,
        CL_minus_task3_single=res.CL - t3_single_CL,
        CL_minus_task3_iterated=res.CL - t3_iter_CL,
        signed_post_correction_diff=closures.signed_post_correction_diff,
        f_semiinf=closures.f_semiinf, f_oracle=closures.f_oracle,
        q_l2=res.q_l2, green_q_vs_task3_trace_rel=res.green_q_vs_trace_rel,
        green_bordered_residual_rel=res.green_res_rel,
        gauge_row_abs=res.gauge_row_abs,
        explicit_potential_residual_rel=res.exp_res_rel,
        blas_threads,
    ))
    return closures
end

"Run both Task 4 routes (single-shot + fixed-geometry iteration) plus the
VelocityThroughSources+Direct control on one frozen wake state.  `t3` is the
matching Task 3 return providing the frozen seed and both oracle CLs."
function _run_task4_case(path, case, artifact_stem, t3, task2, green, p,
        blas_threads; frozen_step=-1, backend=pnl.DirectBackend())
    template, solver, C = t3.finite_setup
    baseline = t3.baseline
    # `task2` is the frozen VelocityThroughSources result (fields .body, .wake,
    # .CL, .probe_velocity); callers normalize flat vs march before passing it.
    seed_body, seed_wake = task2.body, task2.wake
    wake_length_c = case == "march" ? 0.5frozen_step :
        (norm(view(template.Das[1], :, 1)) / p.c + task2.wake.nwakes[] * 0.5)
    cl_old = task2.CL
    t3_single_CL = t3.result.CL
    t3_iter_CL = t3.iteration.result.CL

    # single-shot: solve on the Task 2 active wake strengths
    single = _task4_green_solve(seed_body, seed_wake, solver, green, C, p.Uinf,
        p.rho, p.Sref, p.c, p.dt, baseline.probes, baseline.probe_velocity,
        task2.probe_velocity; label="Task 4 $case single-shot", i_step=0, backend)
    single_cl = _task4_comparison_row!(path, case,
        "green_reconstruction_single_shot", "single_shot", frozen_step,
        wake_length_c, single, baseline.CL, task2.CL, t3_single_CL, t3_iter_CL,
        cl_old, t3_single_CL, blas_threads)

    # fixed-geometry iteration
    iteration = _task4_fixed_geometry_iteration(seed_body, seed_wake, solver,
        green, C, p.Uinf, p.rho, p.Sref, p.c, p.dt, baseline.probes,
        baseline.probe_velocity, task2.probe_velocity;
        label="Task 4 $case iterated", backend)
    iter_cl = _task4_comparison_row!(path, case,
        "green_reconstruction_iterated", iteration.status, frozen_step,
        wake_length_c, iteration.result, baseline.CL, task2.CL, t3_single_CL,
        t3_iter_CL, cl_old, t3_iter_CL, blas_threads)

    # VelocityThroughSources + Direct control (backend-effect isolation)
    control = _frozen_solve(seed_body, seed_wake, solver, C, p.Uinf, p.rho,
        p.Sref, p.c, p.dt, baseline.probes, baseline.probe_velocity;
        label="Task 4 $case VTS-direct control", backend)
    _upsert_comparison(joinpath(path, "comparison.csv"), (;
        task=4, case, formulation="velocity_through_sources_direct_control",
        status="control", frozen_wake_step=frozen_step, wake_length_c,
        wake_rows=control.wake.nwakes[], ncells=control.body.ncells,
        nshedding_edges=length(control.gamma), kutta_constant_error=0.0,
        linear_residual_abs=control.residual_abs,
        linear_residual_rel=control.residual_rel,
        gamma_min=minimum(control.gamma), gamma_max=maximum(control.gamma),
        gamma_mean=sum(control.gamma) / length(control.gamma),
        gamma_l2=norm(control.gamma),
        transition_free_mismatch_l2=control.mismatch_l2,
        transition_free_mismatch_rel=control.mismatch_rel,
        CL=control.CL, CL_minus_task1=control.CL - baseline.CL,
        relative_CL_difference=_relative_difference(control.CL, baseline.CL),
        relative_exterior_velocity_difference=control.probe_rel,
        CL_minus_task2=control.CL - task2.CL,
        relative_CL_difference_task2=_relative_difference(control.CL, task2.CL),
        blas_threads,
    ))

    # per-case artifacts
    single_row = (; route="single_shot", CL=single.CL,
        CL_minus_task1=single.CL - baseline.CL,
        relative_CL_difference=_relative_difference(single.CL, baseline.CL),
        CL_minus_task2=single.CL - task2.CL, CL_minus_task3_single=single.CL - t3_single_CL,
        f_semiinf=single_cl.f_semiinf, f_oracle=single_cl.f_oracle,
        signed_post_correction_diff=single_cl.signed_post_correction_diff,
        green_bordered_residual_rel=single.green_res_rel,
        gauge_row_abs=single.gauge_row_abs,
        explicit_potential_residual_rel=single.exp_res_rel,
        q_l2=single.q_l2, green_q_vs_task3_trace_rel=single.green_q_vs_trace_rel,
        transition_free_mismatch_rel=single.mismatch_rel,
        relative_exterior_velocity_difference_task1=single.probe_rel_task1,
        relative_exterior_velocity_difference_task2=single.probe_rel_task2,
        fixed_point_iterations=0)
    iter_res = iteration.result
    iter_row = (; route="iterated", CL=iter_res.CL,
        CL_minus_task1=iter_res.CL - baseline.CL,
        relative_CL_difference=_relative_difference(iter_res.CL, baseline.CL),
        CL_minus_task2=iter_res.CL - task2.CL,
        CL_minus_task3_single=iter_res.CL - t3_iter_CL,
        f_semiinf=iter_cl.f_semiinf, f_oracle=iter_cl.f_oracle,
        signed_post_correction_diff=iter_cl.signed_post_correction_diff,
        green_bordered_residual_rel=iter_res.green_res_rel,
        gauge_row_abs=iter_res.gauge_row_abs,
        explicit_potential_residual_rel=iter_res.exp_res_rel,
        q_l2=iter_res.q_l2, green_q_vs_task3_trace_rel=iter_res.green_q_vs_trace_rel,
        transition_free_mismatch_rel=iter_res.mismatch_rel,
        relative_exterior_velocity_difference_task1=iter_res.probe_rel_task1,
        relative_exterior_velocity_difference_task2=iter_res.probe_rel_task2,
        fixed_point_iterations=length(iteration.rows))
    _write_rows_csv(joinpath(path, artifact_stem * "_routes.csv"),
        [single_row, iter_row])
    _write_rows_csv(joinpath(path, artifact_stem * "_iteration_history.csv"),
        iteration.rows)
    terminal_time = wake_length_c * p.c / p.U
    _task3_write_terminal_vtk(path, artifact_stem, single, iteration.result,
        terminal_time, max(frozen_step, 0))

    worst_res = maximum((single.green_res_rel, single.exp_res_rel,
        iter_res.green_res_rel, iter_res.exp_res_rel))
    @printf("Task 4 %s: single CL=%.10f (f_oracle=%.4f), iterated CL=%.10f (f_oracle=%.4f), worst residual=%.2e\n",
        case, single.CL, single_cl.f_oracle, iter_res.CL, iter_cl.f_oracle,
        worst_res)
    return (; single, iteration, control, single_closures=single_cl,
        iter_closures=iter_cl, worst_res, wake_length_c, frozen_step)
end

function run_task4_flat(; path=DIRICHLET_DATA_PATH, baseline=nothing,
        alpha_deg=DEFAULT_ALPHA_DEG, green=nothing, blas_threads=nothing)
    mkpath(path)
    blas_threads = isnothing(blas_threads) ? _task4_pin_blas() : blas_threads
    p = _task3_parameters(alpha_deg)
    t3 = run_task3_flat(; path, baseline, alpha_deg)
    green = isnothing(green) ? _task4_shared_green(t3.finite_setup[1],
        t3.finite_setup[2]) : green
    case = _run_task4_case(path, "flat", "task4_flat", t3, t3.task2, green, p,
        blas_threads)
    _write_task4_manifests(path; alpha_deg, blas_threads)
    return (; t3, case, green, blas_threads)
end

function run_task4_flat_das005(; path=DIRICHLET_DATA_PATH, baseline=nothing,
        alpha_deg=DEFAULT_ALPHA_DEG, green=nothing, blas_threads=nothing)
    mkpath(path)
    blas_threads = isnothing(blas_threads) ? _task4_pin_blas() : blas_threads
    p = _task3_parameters(alpha_deg)
    t3 = run_task3_flat_das005(; path, baseline, alpha_deg)
    green = isnothing(green) ? _task4_shared_green(t3.finite_setup[1],
        t3.finite_setup[2]) : green
    case = _run_task4_case(path, "flat_das005", "task4_flat_das005", t3,
        t3.task2, green, p, blas_threads)
    _write_task4_manifests(path; alpha_deg, blas_threads)
    return (; t3, case, green, blas_threads)
end

function run_task4_march(; path=DIRICHLET_DATA_PATH, baseline=nothing,
        alpha_deg=DEFAULT_ALPHA_DEG, green=nothing, blas_threads=nothing)
    mkpath(path)
    blas_threads = isnothing(blas_threads) ? _task4_pin_blas() : blas_threads
    p = _task3_parameters(alpha_deg)
    t3 = run_task3_march(; path, baseline, alpha_deg)
    green = isnothing(green) ? _task4_shared_green(t3.finite_setup[1],
        t3.finite_setup[2]) : green
    frozen_step = t3.task2.selected_step
    case = _run_task4_case(path, "march", "task4_march", t3, t3.task2.result,
        green, p, blas_threads; frozen_step)
    _write_task4_manifests(path; alpha_deg, blas_threads)
    return (; t3, case, green, blas_threads)
end

function _write_task4_manifests(path; alpha_deg=DEFAULT_ALPHA_DEG,
        blas_threads=BLAS.get_num_threads())
    stored = _read_comparison_rows(joinpath(path, "comparison.csv"))
    cases = ("flat", "flat_das005", "march")
    routes = ("green_reconstruction_single_shot", "green_reconstruction_iterated",
        "velocity_through_sources_direct_control")
    statuses = Dict("$route|$case" =>
        (haskey(stored, "$route|$case") ? stored["$route|$case"][:status] : "not_run")
        for route in routes, case in cases)
    config = Dict(
        "task" => 4, "formulation" => "green_reconstruction",
        "formulation_settings" => TASK4_FORMULATION,
        "geometry" => Dict("airfoil" => "NACA 0015", "chord_m" => 1FT_TO_M,
            "span_m" => 4FT_TO_M, "aspect_ratio" => 4.0,
            "n_airfoil" => 161, "n_span" => 13, "n_endcap" => 9),
        "flow" => Dict("Uinf_ft_s" => 330.2, "alpha_deg" => alpha_deg,
            "rho_kg_m3" => 1.225),
        "wake_states" => Dict("flat" => "Da=0.5c, 64c",
            "flat_das005" => "Da=0.05c, 63.55c",
            "march" => "settled rolled Da=0.05c"),
        "gauge" => "area_mean", "green_solver" => "dense bordered LU (nothing)",
        "green_state_shared_across_cases" => true,
        "iteration" => Dict("defect_tolerance" => 1e-8,
            "relative_lift_tolerance" => 1e-8, "required_consecutive" => 3,
            "maximum_iterations" => 200, "initial_omega" => 1.0,
            "minimum_omega" => 0.0625),
        "linear_residual_tolerance" => 1e-10,
        "velocity_only_constraint" =>
            "no scalar_potential requested from free-wake sources at body CPs",
        "backend" => Dict("physical_system" => "DirectBackend",
            "state_reconstruction_march" =>
                "FastMultipoleBackend(order=10, theta=0.4, leaf_size=100)"),
        "blas_threads" => blas_threads, "statuses" => statuses)
    open(joinpath(path, "task4.config.toml"), "w") do io
        TOML.print(io, config; sorted=true)
    end
    terminal_vtk = String[]
    for stem in ("task4_flat", "task4_flat_das005", "task4_march"),
            route in ("single_shot", "iterated"), kind in ("body", "wake")
        push!(terminal_vtk, "dirichlet_$(stem)_$(route)_$(kind).pvd")
    end
    invariants = Dict(
        "all_frozen_solves_preserved_geometry_and_wake" => true,
        "green_bordered_and_explicit_residuals_within_1e-10" => true,
        "green_operator_is_body_mesh_only" => true,
        "velocity_only_no_free_wake_scalar_potential" => true,
        "blas_threads" => blas_threads)
    open(joinpath(path, "task4_invariants.toml"), "w") do io
        TOML.print(io, invariants; sorted=true)
    end
    metadata = Dict(
        "theory_file" => "docs/wake_solve_schemes.md",
        "theory_section" => "Suggested alternative formulation for semiinfinite_wake=false",
        "comparison_csv" => "comparison.csv",
        "routes_csv" => ["task4_flat_routes.csv",
            "task4_flat_das005_routes.csv", "task4_march_routes.csv"],
        "iteration_history_csv" => ["task4_flat_iteration_history.csv",
            "task4_flat_das005_iteration_history.csv",
            "task4_march_iteration_history.csv"],
        "invariants_toml" => "task4_invariants.toml",
        "terminal_vtk" => terminal_vtk, "alpha_deg" => alpha_deg)
    open(joinpath(path, "task4.metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return nothing
end

function run_task4(; path=DIRICHLET_DATA_PATH, alpha_deg=DEFAULT_ALPHA_DEG)
    blas_threads = _task4_pin_blas()
    baseline = _task2_baseline(path; alpha_deg)
    # Build the shared body-only Green state once from a das=0.5c template; it is
    # body-mesh-only, so it serves every wake state and both angles.
    template, solver, _ = _finite_solver_template(1.0 * FT_TO_M, 4.0 * FT_TO_M,
        _uinf_from_alpha(330.2 * FT_TO_M, alpha_deg))
    green = _task4_shared_green(template, solver)
    flat = run_task4_flat(; path, baseline, alpha_deg, green, blas_threads)
    flat_das005 = run_task4_flat_das005(; path, baseline, alpha_deg, green,
        blas_threads)
    march = run_task4_march(; path, baseline, alpha_deg, green, blas_threads)
    _write_task4_manifests(path; alpha_deg, blas_threads)
    return (; baseline, flat, flat_das005, march, green, blas_threads)
end

################################################################################
# Green's-third-identity residual diagnostic (frozen finite-wake states).
#
# Additive, read-only-input diagnostic requested as a Task-4 follow-up.  It
# reconstructs each frozen free-wake from saved VTS (no re-marching) and measures
# how well the discrete Green identity  Sσ = (I−B)q  holds when q is evaluated
# DIRECTLY as the wake-induced potential trace on the body (not reconstructed),
# plus reports reconstructed-vs-direct q.  The residual depends only on the body
# mesh and the free-wake state, so it is shared by VelocityThroughSources and
# GreenReconstruction on the single-shot (seed) wake and diverges only on the
# per-method iterated wakes.  One template body + one shared bordered-LU serve
# every state and both angles (both are body-mesh-only; Da is irrelevant here).
################################################################################

"Structural wake options per frozen state (rolled carries a final filament; the
flat convergence states are strictly panel-only).  Matches the saved states so
the loaded PanelWake reproduces the Task-2/3/4 source set exactly."
_green_id_wake_opts(state::Symbol; core_size=0.001) = state === :rolled ?
    (; core_size, include_final_filament=true,
       shed_with_induced_velocity=true, unsteady_filament=true) :
    (; core_size, include_final_filament=false,
       shed_with_induced_velocity=true, unsteady_filament=true)

"Build a PanelWake on `body` matching a frozen state's structure and load its
saved nodes/strengths from `dir/stem/stem.<i_surf>.<idx>.vts`."
function _green_id_load_wake(body, dir, stem, idx, state::Symbol; nwakerows=200,
        core_size=0.001)
    opts = _green_id_wake_opts(state; core_size)
    wake = pnl.PanelWake(body; nwakerows, opts...)
    pnl.update_TE!(wake, body)
    pnl._load_panel_wake_vtk!(wake, dir, stem, idx)
    all(all(isfinite, s) for s in wake.strength) ||
        error("green-identity: non-finite loaded wake strength ($stem)")
    all(all(isfinite, n) for n in wake.nodes) ||
        error("green-identity: non-finite loaded wake nodes ($stem)")
    return wake
end

"Shared reducer: given the raw direct trace and the wake-only source strengths,
form Sσ, the identity residual r = Sσ − (I−B)q_direct (area-mean removed from
q_direct first), the reconstruction q_recon from the same bordered LU, and the
reconstructed-vs-direct q.  Uses `_lu_matrix_product` to apply the bordered
operator K forward (its top block is (I−B))."
function _green_identity_metrics(green, areas, q_direct_raw, sigma, body; label,
        backend=pnl.DirectBackend())
    N = body.ncells
    q_direct = _remove_area_mean(q_direct_raw, areas)
    Ssigma = zeros(Float64, N)
    saved_str = copy(body.strength); saved_pot = copy(body.potential)
    try
        pnl._source_potential!(Ssigma, body, sigma, backend)      # S·σ
    finally
        body.strength .= saved_str; body.potential .= saved_pot
    end
    Ssigma_l2 = norm(Ssigma)
    # (I−B)·q_direct = top block of K·[q_direct; 0]
    ImB_q = view(_lu_matrix_product(green.fact, vcat(q_direct, 0.0)), 1:N)
    resid_rel = norm(Ssigma .- ImB_q) / max(Ssigma_l2, eps(Float64))
    # reconstruction from the same bordered LU (mutates green buffers; serial use).
    # The gauge-fixed solve satisfies (I−B)q + λa = Sσ, aᵀq = 0, so the by-construction
    # check is the FULL bordered residual ‖K·[q;λ] − [Sσ;0]‖ (~machine eps); the
    # top-block-only quantity Sσ − (I−B)q equals λa (the finite-resolution Green
    # compatibility defect) and is a physical diagnostic, not a solve error.
    q_recon = copy(pnl._green_solve_q!(green, Ssigma))
    bordered = _lu_matrix_product(green.fact, green.sol_b) .- green.rhs_b
    recon_solve_rel = norm(bordered) / max(norm(green.rhs_b), eps(Float64))
    recon_solve_rel <= 1e-8 ||
        error("$label bordered reconstruction residual $recon_solve_rel exceeds 1e-8")
    compat_defect_rel = norm(Ssigma .-
        view(_lu_matrix_product(green.fact, vcat(q_recon, 0.0)), 1:N)) /
        max(Ssigma_l2, eps(Float64))
    q_recon_vs_direct_rel = norm(_remove_area_mean(q_recon, areas) .- q_direct) /
        max(norm(q_direct), eps(Float64))
    all(isfinite, (Ssigma_l2, resid_rel, q_recon_vs_direct_rel, recon_solve_rel)) ||
        error("$label produced a non-finite diagnostic")
    return (; Ssigma_l2, q_direct_l2=norm(q_direct),
        green_identity_residual_rel=resid_rel, q_recon_vs_direct_rel,
        compat_defect_rel, recon_solve_rel)
end

"Panel-wake case: q_direct = wake scalar-potential trace (Task-3 helper); σ from
wake-only velocity at control points.  Frozen-state immutability asserted."
function _green_identity_panel_case(body, wake, green, areas; label,
        backend=pnl.DirectBackend())
    N = body.ncells
    snapshot = _state_snapshot(body, wake)
    q_direct_raw = _task3_accumulated_potential(body, (wake,))   # wake-only φ
    sigma = zeros(Float64, N); sigma0 = zeros(Float64, N)
    saved_vel = copy(body.velocity)
    try
        body.velocity .= 0.0
        pnl.calc_normals!(body); pnl.calc_controlpoints!(body)
        pnl.influence!((body,), (wake,), backend;
            scalar_potential=false, velocity=true)              # wake-only u_f
        pnl._split_sigma!(sigma, sigma0, body, zeros(3, N))     # σ = −u_f·n
    finally
        body.velocity .= saved_vel
    end
    all(isfinite, sigma) || error("$label produced non-finite σ")
    m = _green_identity_metrics(green, areas, q_direct_raw, sigma, body; label,
        backend)
    _assert_frozen_state(snapshot, body, wake; label)
    return merge(m, (; sigma=copy(sigma), q_direct_raw=copy(q_direct_raw)))
end

"Semi-infinite reference: isolate the attached sheet's induced potential/velocity
via the suppress_attached_wake on/off difference (technique from `_assemble_W!`).
The sheet touches the trailing edge, so its support is NOT disjoint from the
body — the identity residual here is an interpretive reference, not a clean check."
function _green_identity_semiinf(path, alpha_deg, green, areas; label,
        backend=pnl.DirectBackend())
    baseline = run_task1(; path, alpha_deg)      # single steady solve, no march
    body = deepcopy(baseline.body)
    return _green_identity_semiinf_core(body, green, areas; label, backend)
end

"Isolate the attached-sheet trace/velocity on an already-solved semi-infinite
body (`body.strength[:,2]=μ`) and evaluate the Green identity.  Extracted so the
refinement ladder can feed bodies solved at each resolution."
function _green_identity_semiinf_core(body, green, areas; label,
        backend=pnl.DirectBackend(), velocity_core_size=nothing)
    N = body.ncells
    pnl.calc_normals!(body); pnl.calc_controlpoints!(body)
    was = body.suppress_attached_wake[]
    saved_pot = copy(body.potential); saved_vel = copy(body.velocity)
    orig_offset = body.kerneloffset          # regularized velocity core for the sheet
    local q_direct_raw, sigma
    try
        # scalar potential: full (attached wake on) minus body-only (wake off).
        # Always the original active kernel — the potential formula is analytic
        # and does not carry the velocity core.
        fill!(body.potential, 0.0); body.suppress_attached_wake[] = false
        pnl.influence!((body,), (body,), backend;
            scalar_potential=true, velocity=false)
        phi_full = copy(body.potential)
        fill!(body.potential, 0.0); body.suppress_attached_wake[] = true
        pnl.influence!((body,), (body,), backend;
            scalar_potential=true, velocity=false)
        phi_body = copy(body.potential)
        q_direct_raw = phi_full .- phi_body
        # velocity: same on/off difference → wake-only u_f.  For the unregularized
        # probe only, drop the sheet's velocity kernel offset around both calls.
        try
            velocity_core_size !== nothing && (body.kerneloffset = velocity_core_size)
            body.velocity .= 0.0; body.suppress_attached_wake[] = false
            pnl.influence!((body,), (body,), backend;
                scalar_potential=false, velocity=true)
            u_full = copy(body.velocity)
            body.velocity .= 0.0; body.suppress_attached_wake[] = true
            pnl.influence!((body,), (body,), backend;
                scalar_potential=false, velocity=true)
            u_body = copy(body.velocity)
            body.velocity .= u_full .- u_body
            sigma = zeros(Float64, N); sigma0 = zeros(Float64, N)
            pnl._split_sigma!(sigma, sigma0, body, zeros(3, N))
        finally
            body.kerneloffset = orig_offset
        end
    finally
        body.suppress_attached_wake[] = was
        body.potential .= saved_pot; body.velocity .= saved_vel
    end
    all(isfinite, q_direct_raw) && all(isfinite, sigma) ||
        error("$label produced a non-finite semi-infinite field")
    m = _green_identity_metrics(green, areas, q_direct_raw, sigma, body; label,
        backend)
    return merge(m, (; sigma=copy(sigma), q_direct_raw=copy(q_direct_raw),
        velocity_core_size_used=orig_offset))
end

"Read a Task-4 `green_q_vs_task3_trace_rel` value for a case/route as an
independent cross-check of the reconstructed-vs-direct q (identical math)."
function _green_id_task4_xcheck(path, case, route_kind)
    csv = joinpath(path, "task4_$(case)_routes.csv")
    isfile(csv) || return NaN
    lines = readlines(csv)
    length(lines) >= 2 || return NaN
    header = split(first(lines), ',')
    ri = findfirst(==("route"), header)
    ti = findfirst(==("green_q_vs_task3_trace_rel"), header)
    (isnothing(ri) || isnothing(ti)) && return NaN
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        vals = split(line, ','; keepempty=true)
        length(vals) >= max(ri, ti) || continue
        vals[ri] == route_kind || continue
        v = strip(vals[ti])
        return isempty(v) ? NaN : parse(Float64, v)
    end
    return NaN
end

function run_green_identity_residual(; path=DIRICHLET_DATA_PATH,
        alpha_deg=DEFAULT_ALPHA_DEG)
    mkpath(path)
    blas_threads = _task4_pin_blas()
    p = _task3_parameters(alpha_deg)
    # One body-mesh-only template + shared bordered-LU Green state (Da-independent,
    # angle-independent); reused for every wake state below.
    template, solver, _ = _finite_solver_template(p.c, p.b, p.Uinf)
    green = _task4_shared_green(template, solver)
    areas = pnl._panel_areas(template)
    body = deepcopy(template)
    pnl.calc_normals!(body); pnl.calc_controlpoints!(body)

    # constant-mode sanity: B·1≈1, so (I−B)·1≈0; ‖(I−B)·1‖∞ is the finite-
    # resolution defect (how close the constant is to an exact null vector)
    ones_top = view(_lu_matrix_product(green.fact,
        vcat(ones(Float64, body.ncells), 0.0)), 1:body.ncells)
    const_defect = norm(ones_top, Inf)

    # (wake_state, route label, collection prefix, collection suffix, t4 route)
    routes = (
        (:seed, "shared_single_shot", "dirichlet_task2_", "_terminal_wake",
            "single_shot"),
        (:green_reconstruction, "green_reconstruction_iterated",
            "dirichlet_task4_", "_iterated_wake", "iterated"),
        (:oracle, "task3_oracle_iterated", "dirichlet_task3_",
            "_iterated_wake", ""),
    )
    # (case name, wake-structure state, saved VTS index)
    states = (("flat", :flat, 0), ("flat_das005", :flat, 0), ("march", :rolled, 80))

    rows = NamedTuple[]
    for (wake_state, route_label, prefix, suffix, t4route) in routes
        for (case, structure, idx) in states
            stem = "$(prefix)$(case)$(suffix)"
            label = "green-identity $(alpha_deg)° $wake_state/$case"
            wake = _green_id_load_wake(body, path, stem, idx, structure)
            m = _green_identity_panel_case(body, wake, green, areas; label)
            xref = wake_state === :seed || wake_state === :green_reconstruction ?
                _green_id_task4_xcheck(path, case, t4route) : NaN
            xdelta = isnan(xref) ? "" :
                abs(m.q_recon_vs_direct_rel - xref)
            push!(rows, (; angle=alpha_deg, wake_state=String(wake_state),
                route=route_label, source_method="panel_wake",
                n_panels=body.ncells,
                Ssigma_l2=m.Ssigma_l2, q_direct_l2=m.q_direct_l2,
                green_identity_residual_rel=m.green_identity_residual_rel,
                q_recon_vs_direct_rel=m.q_recon_vs_direct_rel,
                green_compat_defect_rel=m.compat_defect_rel,
                recon_solve_rel=m.recon_solve_rel,
                task4_trace_xcheck_delta=xdelta, notes=case))
        end
    end

    # semi-infinite reference (one row): interpretive — disjoint-support violated
    ms = _green_identity_semiinf(path, alpha_deg, green, areas;
        label="green-identity $(alpha_deg)° semiinf")
    push!(rows, (; angle=alpha_deg, wake_state="semiinf",
        route="task1_semiinfinite", source_method="semiinf_attached",
        n_panels=body.ncells, Ssigma_l2=ms.Ssigma_l2,
        q_direct_l2=ms.q_direct_l2,
        green_identity_residual_rel=ms.green_identity_residual_rel,
        q_recon_vs_direct_rel=ms.q_recon_vs_direct_rel,
        green_compat_defect_rel=ms.compat_defect_rel,
        recon_solve_rel=ms.recon_solve_rel, task4_trace_xcheck_delta="",
        notes="semiinf_disjoint_support_TE"))

    csv_path = joinpath(path, "green_identity_residual.csv")
    _write_rows_csv(csv_path, rows)

    @printf("\nGreen's-third-identity residual  (alpha = %g deg, N = %d, backend = Direct)\n",
        alpha_deg, body.ncells)
    @printf("  (I-B)*1 constant-mode defect |.|inf = %.3e ; BLAS threads = %d\n",
        const_defect, blas_threads)
    @printf("  %-22s %-13s %11s %13s %11s %11s %10s\n",
        "wake_state", "case", "qdir_l2",
        "resid_rel", "qrec_vs_dir", "compat_def", "t4_xdel")
    for r in rows
        xd = r.task4_trace_xcheck_delta == "" ? "   --" :
            @sprintf("%.1e", r.task4_trace_xcheck_delta)
        @printf("  %-22s %-13s %11.3e %13.4e %11.4e %11.3e %10s\n",
            r.wake_state, r.notes, r.q_direct_l2,
            r.green_identity_residual_rel, r.q_recon_vs_direct_rel,
            r.green_compat_defect_rel, xd)
    end
    println("  wrote $(abspath(csv_path))")
    println("  note: semiinf row is interpretive (attached sheet touches the TE; " *
        "Green disjoint-support assumption violated).")
    return (; rows, csv_path, const_defect, blas_threads, alpha_deg)
end

# -----------------------------------------------------------------------------
# Follow-up: unregularized wake-velocity kernel Green-identity comparison.
#
# Paired probe.  For each frozen state we build two otherwise-identical wakes and
# evaluate the Green identity with the doublet-panel VELOCITY core on (core_size
# = 0.001) and off (core_size = 0.0).  The analytic doublet POTENTIAL never uses
# the velocity core, so q_direct must be invariant (gate 2).  DirectBackend only.
# -----------------------------------------------------------------------------

"Assemble one wide paired row from a regularized (velocity core 0.001) and an
unregularized (velocity core 0.0) Green-identity result.  Only the velocity core
should differ; enforce finiteness (gate 5) and analytic-potential invariance
(gate 2).  Ratios and percent change are unregularized-over-regularized, so a
negative percent change means the singular velocity improves the residual."
function _regularization_paired_metrics(mreg, munreg, reg_core, unreg_core; label,
        potential_gate=1e-13)
    r_reg = mreg.green_identity_residual_rel
    r_unreg = munreg.green_identity_residual_rel
    ratio = r_unreg / max(r_reg, eps(Float64))
    pct = 100 * (ratio - 1)
    sig_reg = mreg.sigma; sig_unreg = munreg.sigma
    sigma_rel = norm(sig_unreg .- sig_reg) / max(norm(sig_reg), eps(Float64))
    qd_reg = mreg.q_direct_raw; qd_unreg = munreg.q_direct_raw
    qdir_rel = norm(qd_unreg .- qd_reg) / max(norm(qd_reg), eps(Float64))
    trace_ratio = munreg.q_recon_vs_direct_rel /
        max(mreg.q_recon_vs_direct_rel, eps(Float64))
    # gate 5: every sampled field and reported metric must be finite
    (all(isfinite, sig_reg) && all(isfinite, sig_unreg) &&
        all(isfinite, qd_reg) && all(isfinite, qd_unreg)) ||
        error("$label produced a non-finite paired field (core=$unreg_core unsafe)")
    all(isfinite, (r_reg, r_unreg, ratio, pct, sigma_rel, qdir_rel, trace_ratio)) ||
        error("$label produced a non-finite paired metric")
    # gate 2: the analytic potential does not use the velocity core → invariant
    qdir_rel <= potential_gate || error(
        "$label potential changed with velocity core: q_direct rel=$qdir_rel > $potential_gate")
    return (;
        regularized_core_size=reg_core, unregularized_core_size=unreg_core,
        regularized_Ssigma_l2=mreg.Ssigma_l2,
        unregularized_Ssigma_l2=munreg.Ssigma_l2,
        regularized_residual_rel=r_reg, unregularized_residual_rel=r_unreg,
        residual_ratio_unreg_over_reg=ratio, residual_percent_change=pct,
        regularized_q_recon_vs_direct_rel=mreg.q_recon_vs_direct_rel,
        unregularized_q_recon_vs_direct_rel=munreg.q_recon_vs_direct_rel,
        trace_gap_ratio_unreg_over_reg=trace_ratio,
        sigma_relative_change=sigma_rel, q_direct_relative_change=qdir_rel,
        regularized_compat_defect_rel=mreg.compat_defect_rel,
        unregularized_compat_defect_rel=munreg.compat_defect_rel,
        regularized_recon_solve_rel=mreg.recon_solve_rel,
        unregularized_recon_solve_rel=munreg.recon_solve_rel)
end

"Gate 6 cross-pair audit: paired wakes must be identical apart from core_size."
function _assert_wake_pair_equal(wa, wb; label)
    wa.nwakes[] == wb.nwakes[] || error("$label paired wake row count differs")
    all(isequal(a, b) for (a, b) in zip(wa.nodes, wb.nodes)) ||
        error("$label paired wake nodes differ")
    all(isequal(a, b) for (a, b) in zip(wa.strength, wb.strength)) ||
        error("$label paired wake strengths differ")
    (wa.include_final_filament == wb.include_final_filament &&
        wa.shed_with_induced_velocity == wb.shed_with_induced_velocity &&
        wa.unsteady_filament == wb.unsteady_filament) ||
        error("$label paired wake options differ beyond core_size")
    return true
end

"Read established regularized residuals keyed by (wake_state, notes/case) for the
baseline-reproduction gate (gate 3).  Empty if the artifact is absent."
function _green_baseline_residuals(path)
    out = Dict{Tuple{String,String},Float64}()
    for row in _read_simple_csv(joinpath(path, "green_identity_residual.csv"))
        haskey(row, "green_identity_residual_rel") || continue
        v = tryparse(Float64, row["green_identity_residual_rel"])
        isnothing(v) && continue
        out[(get(row, "wake_state", ""), get(row, "notes", ""))] = v
    end
    return out
end

"Gate 3: the regularized branch must reproduce the established residual to
rtol=1e-10, atol=1e-13.  Warn (do not fail) when no baseline artifact exists."
function _check_regularized_baseline(baseline, key, value; label)
    if isempty(baseline)
        return nothing
    elseif haskey(baseline, key)
        isapprox(value, baseline[key]; rtol=1e-10, atol=1e-13) || error(
            "$label regularized residual $value does not reproduce baseline " *
            "$(baseline[key]) for $key")
    else
        @warn "$label: no baseline residual row for $key (gate 3 skipped for it)"
    end
    return nothing
end

function run_green_identity_regularization(; path=DIRICHLET_DATA_PATH,
        alpha_deg=DEFAULT_ALPHA_DEG)
    mkpath(path)
    blas_threads = _task4_pin_blas()
    p = _task3_parameters(alpha_deg)
    template, solver, _ = _finite_solver_template(p.c, p.b, p.Uinf)
    green = _task4_shared_green(template, solver)
    areas = pnl._panel_areas(template)
    body = deepcopy(template)
    pnl.calc_normals!(body); pnl.calc_controlpoints!(body)
    baseline = _green_baseline_residuals(path)

    routes = (
        (:seed, "shared_single_shot", "dirichlet_task2_", "_terminal_wake"),
        (:green_reconstruction, "green_reconstruction_iterated",
            "dirichlet_task4_", "_iterated_wake"),
        (:oracle, "task3_oracle_iterated", "dirichlet_task3_", "_iterated_wake"),
    )
    states = (("flat", :flat, 0), ("flat_das005", :flat, 0), ("march", :rolled, 80))

    rows = NamedTuple[]
    for (wake_state, route_label, prefix, suffix) in routes
        for (case, structure, idx) in states
            stem = "$(prefix)$(case)$(suffix)"
            label = "green-reg $(alpha_deg)° $wake_state/$case"
            wreg = _green_id_load_wake(body, path, stem, idx, structure;
                core_size=0.001)
            wunreg = _green_id_load_wake(body, path, stem, idx, structure;
                core_size=0.0)
            _assert_wake_pair_equal(wreg, wunreg; label)
            mreg = _green_identity_panel_case(body, wreg, green, areas; label)
            munreg = _green_identity_panel_case(body, wunreg, green, areas; label)
            pair = _regularization_paired_metrics(mreg, munreg, 0.001, 0.0; label)
            _check_regularized_baseline(baseline, (String(wake_state), case),
                pair.regularized_residual_rel; label)
            push!(rows, merge((; angle=alpha_deg, wake_state=String(wake_state),
                route=route_label, case, n_panels=body.ncells), pair))
        end
    end

    # semi-infinite attached sheet: no PanelWake core; the velocity core is the
    # body's active kerneloffset.  Regularized keeps it; unregularized zeroes it
    # around the velocity calls only.  Solve once, evaluate both.
    label = "green-reg $(alpha_deg)° semiinf"
    t1 = run_task1(; path, alpha_deg)
    msi_reg = _green_identity_semiinf_core(deepcopy(t1.body), green, areas;
        label, velocity_core_size=nothing)
    msi_unreg = _green_identity_semiinf_core(deepcopy(t1.body), green, areas;
        label, velocity_core_size=0.0)
    pair = _regularization_paired_metrics(msi_reg, msi_unreg,
        msi_reg.velocity_core_size_used, 0.0; label)
    _check_regularized_baseline(baseline, ("semiinf", "semiinf_disjoint_support_TE"),
        pair.regularized_residual_rel; label)
    push!(rows, merge((; angle=alpha_deg, wake_state="semiinf",
        route="task1_semiinfinite", case="semiinf_disjoint_support_TE",
        n_panels=body.ncells), pair))

    csv_path = joinpath(path, "green_identity_regularization.csv")
    _write_rows_csv(csv_path, rows)

    @printf("\nGreen identity: regularized vs unregularized velocity kernel  (alpha = %g deg, N = %d)\n",
        alpha_deg, body.ncells)
    @printf("  %-22s %-13s %13s %13s %9s %11s\n",
        "wake_state", "case", "resid_reg", "resid_unreg", "ratio", "sig_change")
    for r in rows
        @printf("  %-22s %-13s %13.4e %13.4e %9.4f %11.3e\n",
            r.wake_state, r.case, r.regularized_residual_rel,
            r.unregularized_residual_rel, r.residual_ratio_unreg_over_reg,
            r.sigma_relative_change)
    end
    println("  wrote $(abspath(csv_path))")
    println("  (ratio<1 ⇒ singular velocity reduces the residual; " *
        "gate 2 verified q_direct invariant to the velocity core)")
    return (; rows, csv_path, blas_threads, alpha_deg)
end

# -----------------------------------------------------------------------------
# Probe A: mesh-refinement of the Green third-identity residual
# Probe B: conditioning / near-null structure of (I−B)
#
# Both operate on the body-mesh-only Green operator at each resolution of the
# Task-1 grid ladder.  Read-only; no time-marching; DirectBackend throughout.
# -----------------------------------------------------------------------------

const _GREEN_REFINEMENT_GRIDS = (
    (n_airfoil=81, n_span=7, n_endcap=5),
    (n_airfoil=121, n_span=10, n_endcap=7),
    (n_airfoil=161, n_span=13, n_endcap=9),   # production, ~6,688 panels
    (n_airfoil=201, n_span=16, n_endcap=11),  # ~10,360 panels (run last)
)

"Build a finite RigidWakeBody at `grid` resolution with the attached transition
of length `das_length_c·c`; geometry refreshed."
function _refinement_finite_body(c, b, Uinf, grid; das_length_c=0.5)
    body = build_pitching_wing_body(c, b;
        n_span=grid.n_span, n_airfoil=grid.n_airfoil, n_endcap=grid.n_endcap,
        endcap=:round, semiinfinite_wake=false)
    set_wake_Das!(body, Uinf; magnitude=das_length_c * c)
    pnl.calc_normals!(body); pnl.calc_controlpoints!(body)
    return body
end

"Build and steady-solve a semi-infinite body at `grid` resolution; returns the
solved body plus its transition circulation γ=C·μ (the deterministic prescribed
strength for the flat states at the same resolution)."
function _refinement_semiinf_solve(c, b, Uinf, grid, rho)
    U = norm(Uinf)
    Sref = b * c
    body = build_pitching_wing_body(c, b;
        n_span=grid.n_span, n_airfoil=grid.n_airfoil, n_endcap=grid.n_endcap,
        endcap=:round, semiinfinite_wake=true)
    set_wake_Das!(body, Uinf)
    C = _build_kutta_map(body)
    norm(C * ones(body.ncells), Inf) == 0.0 || error("invalid refinement Kutta map")
    direct = pnl.DirectBackend()
    solver = pnl.Backslash(body)
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false, backend=direct)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c),
        correct_kuttacondition=false, verbose=false)
    pnl.steady!(body, pnl.ReferenceFrame(body), Uinf;
        body_solvers=solver, backend=direct, backend_solve=direct,
        backend_system=direct, monitors=(pressure, force),
        path=nothing, name="dirichlet_green_refinement_semiinf",
        dt=c / U * 0.5, verbose=false)
    gamma = C * copy(view(body.strength, :, 2))
    all(isfinite, gamma) || error("non-finite refinement semi-infinite γ")
    return (; body, gamma, C)
end

_refinement_row(grid, n_panels, h_proxy, angle, state, wake_length_c, m, sigma_min) = (;
    n_airfoil=grid.n_airfoil, n_span=grid.n_span, n_endcap=grid.n_endcap,
    n_panels, h_proxy, angle, state, wake_length_c,
    green_identity_residual_rel=m.green_identity_residual_rel,
    q_recon_vs_direct_rel=m.q_recon_vs_direct_rel,
    green_compat_defect_rel=m.compat_defect_rel,
    sigma_min, recon_solve_rel=m.recon_solve_rel)

"Least-squares slope of log(y) vs log(h) (the observed convergence rate)."
function _loglog_slope(h, y)
    (length(h) >= 2 && all(>(0), y)) || return NaN
    lx = log.(h); ly = log.(y)
    mx = sum(lx) / length(lx); my = sum(ly) / length(ly)
    denom = sum(abs2, lx .- mx)
    denom > 0 || return NaN
    return dot(lx .- mx, ly .- my) / denom
end

function run_green_identity_refinement(; path=DIRICHLET_DATA_PATH,
        alpha_degs=(DEFAULT_ALPHA_DEG, 45.0), max_panels=Inf)
    mkpath(path)
    blas_threads = _task4_pin_blas()
    c = 1.0 * FT_TO_M; b = 4.0 * c; U = 330.2 * FT_TO_M; rho = 1.225
    wake_length_c = 32.0        # fixed free-wake length across resolutions
    rows = NamedTuple[]

    for grid in _GREEN_REFINEMENT_GRIDS
        # Green LU is body-mesh-only (attached wake suppressed in B, Da- and
        # angle-independent): build once per resolution from a das=0.5 template.
        template = _refinement_finite_body(c, b,
            _uinf_from_alpha(U, first(alpha_degs)), grid; das_length_c=0.5)
        n_panels = template.ncells
        if n_panels > max_panels
            @printf("[refinement] skip grid n_panels=%d > max_panels=%g\n",
                n_panels, max_panels)
            template = nothing; GC.gc()
            continue
        end
        solver = pnl.Backslash(template)
        green = _task4_shared_green(template, solver)
        areas = pnl._panel_areas(template)
        # sigma_min proxy = ‖(I−B)·1‖₂/‖1‖₂ (how near the constant is to a null vector)
        ones_top = view(_lu_matrix_product(green.fact,
            vcat(ones(Float64, n_panels), 0.0)), 1:n_panels)
        sigma_min = norm(ones_top) / sqrt(n_panels)
        h_proxy = 1 / sqrt(n_panels)
        @printf("\n[refinement] grid (%d,%d,%d): n_panels=%d, sigma_min=%.3e, BLAS=%d\n",
            grid.n_airfoil, grid.n_span, grid.n_endcap, n_panels, sigma_min,
            blas_threads)

        for angle in alpha_degs
            Uinf = _uinf_from_alpha(U, angle)
            semi = _refinement_semiinf_solve(c, b, Uinf, grid, rho)
            gamma = semi.gamma

            f05 = _refinement_finite_body(c, b, Uinf, grid; das_length_c=0.5)
            w05 = _flat_wake(f05, gamma, wake_length_c, c)
            m05 = _green_identity_panel_case(f05, w05, green, areas;
                label="refine $(angle)° flat_da0.5 N=$n_panels")
            push!(rows, _refinement_row(grid, n_panels, h_proxy, angle,
                "flat_da0.5c", wake_length_c, m05, sigma_min))

            f005 = _refinement_finite_body(c, b, Uinf, grid; das_length_c=0.05)
            w005 = _flat_wake(f005, gamma, wake_length_c, c)
            m005 = _green_identity_panel_case(f005, w005, green, areas;
                label="refine $(angle)° flat_da0.05 N=$n_panels")
            push!(rows, _refinement_row(grid, n_panels, h_proxy, angle,
                "flat_da0.05c", wake_length_c, m005, sigma_min))

            msi = _green_identity_semiinf_core(semi.body, green, areas;
                label="refine $(angle)° semiinf N=$n_panels")
            push!(rows, _refinement_row(grid, n_panels, h_proxy, angle,
                "semiinf", Inf, msi, sigma_min))

            @printf("  angle=%5.2f°  resid_rel: flat0.5c=%.4e  flat0.05c=%.4e  semiinf=%.4e\n",
                angle, m05.green_identity_residual_rel,
                m005.green_identity_residual_rel, msi.green_identity_residual_rel)
            @printf("               qrec_vs_dir: flat0.5c=%.4e  flat0.05c=%.4e  semiinf=%.4e\n",
                m05.q_recon_vs_direct_rel, m005.q_recon_vs_direct_rel,
                msi.q_recon_vs_direct_rel)
        end
        template = nothing; solver = nothing; green = nothing; areas = nothing
        GC.gc()
    end

    csv_path = joinpath(path, "green_identity_refinement.csv")
    _write_rows_csv(csv_path, rows)

    println("\nGreen identity residual — mesh refinement (log-log slope of resid_rel vs h=1/√N)")
    @printf("  %-12s %-7s %10s %10s %10s\n", "state", "angle",
        "slope_r", "slope_qrec", "slope_compat")
    for state in ("flat_da0.5c", "flat_da0.05c", "semiinf"), angle in alpha_degs
        sel = filter(r -> r.state == state && r.angle == angle, rows)
        length(sel) >= 2 || continue
        h = [r.h_proxy for r in sel]
        sr = _loglog_slope(h, [r.green_identity_residual_rel for r in sel])
        sq = _loglog_slope(h, [r.q_recon_vs_direct_rel for r in sel])
        sc = _loglog_slope(h, [r.green_compat_defect_rel for r in sel])
        @printf("  %-12s %5.2f° %10.3f %10.3f %10.3f\n", state, angle, sr, sq, sc)
    end
    println("  (slope>0 ⇒ residual decays under refinement; slope≈0 ⇒ plateau)")
    println("  wrote $(abspath(csv_path))")
    return (; rows, csv_path, blas_threads)
end

"Established refinement residuals keyed by (state, angle, n_panels) for gate 3."
function _green_baseline_refinement(path)
    out = Dict{Tuple{String,Float64,Int},Float64}()
    for row in _read_simple_csv(joinpath(path, "green_identity_refinement.csv"))
        haskey(row, "green_identity_residual_rel") || continue
        st = get(row, "state", ""); a = tryparse(Float64, get(row, "angle", ""))
        n = tryparse(Int, get(row, "n_panels", ""))
        v = tryparse(Float64, row["green_identity_residual_rel"])
        (isnothing(a) || isnothing(n) || isnothing(v)) && continue
        out[(st, a, n)] = v
    end
    return out
end

function run_green_identity_regularization_refinement(; path=DIRICHLET_DATA_PATH,
        alpha_degs=(DEFAULT_ALPHA_DEG, 45.0), max_panels=Inf)
    mkpath(path)
    blas_threads = _task4_pin_blas()
    c = 1.0 * FT_TO_M; b = 4.0 * c; U = 330.2 * FT_TO_M; rho = 1.225
    wake_length_c = 32.0
    baseline = _green_baseline_refinement(path)
    rows = NamedTuple[]

    for grid in _GREEN_REFINEMENT_GRIDS
        template = _refinement_finite_body(c, b,
            _uinf_from_alpha(U, first(alpha_degs)), grid; das_length_c=0.5)
        n_panels = template.ncells
        if n_panels > max_panels
            @printf("[reg-refinement] skip grid n_panels=%d > max_panels=%g\n",
                n_panels, max_panels)
            template = nothing; GC.gc()
            continue
        end
        solver = pnl.Backslash(template)
        green = _task4_shared_green(template, solver)
        areas = pnl._panel_areas(template)
        ones_top = view(_lu_matrix_product(green.fact,
            vcat(ones(Float64, n_panels), 0.0)), 1:n_panels)
        sigma_min = norm(ones_top) / sqrt(n_panels)
        h_proxy = 1 / sqrt(n_panels)
        # gate 7: the finest (201,16,11) grid has a known operator defect
        mesh_valid_for_fit = grid.n_airfoil != 201
        @printf("\n[reg-refinement] grid (%d,%d,%d): n_panels=%d, sigma_min=%.3e, valid_for_fit=%s\n",
            grid.n_airfoil, grid.n_span, grid.n_endcap, n_panels, sigma_min,
            mesh_valid_for_fit)

        for angle in alpha_degs
            Uinf = _uinf_from_alpha(U, angle)
            semi = _refinement_semiinf_solve(c, b, Uinf, grid, rho)
            gamma = semi.gamma

            for (state, das_length_c) in (("flat_da0.5c", 0.5), ("flat_da0.05c", 0.05))
                fbody = _refinement_finite_body(c, b, Uinf, grid; das_length_c)
                wreg = _flat_wake(fbody, gamma, wake_length_c, c; core_size=1e-3)
                wunreg = _flat_wake(fbody, gamma, wake_length_c, c; core_size=0.0)
                lbl = "reg-refine $(angle)° $state N=$n_panels"
                _assert_wake_pair_equal(wreg, wunreg; label=lbl)
                mreg = _green_identity_panel_case(fbody, wreg, green, areas; label=lbl)
                munreg = _green_identity_panel_case(fbody, wunreg, green, areas; label=lbl)
                pair = _regularization_paired_metrics(mreg, munreg, 1e-3, 0.0; label=lbl)
                _check_regularized_baseline_refine(baseline, (state, angle, n_panels),
                    pair.regularized_residual_rel; label=lbl)
                push!(rows, merge((;
                    n_airfoil=grid.n_airfoil, n_span=grid.n_span,
                    n_endcap=grid.n_endcap, n_panels, h_proxy, angle, state,
                    wake_length_c, mesh_valid_for_fit), pair,
                    (; sigma_min)))
            end

            lbl = "reg-refine $(angle)° semiinf N=$n_panels"
            msi_reg = _green_identity_semiinf_core(deepcopy(semi.body), green, areas;
                label=lbl, velocity_core_size=nothing)
            msi_unreg = _green_identity_semiinf_core(deepcopy(semi.body), green, areas;
                label=lbl, velocity_core_size=0.0)
            pair = _regularization_paired_metrics(msi_reg, msi_unreg,
                msi_reg.velocity_core_size_used, 0.0; label=lbl)
            _check_regularized_baseline_refine(baseline, ("semiinf", angle, n_panels),
                pair.regularized_residual_rel; label=lbl)
            push!(rows, merge((;
                n_airfoil=grid.n_airfoil, n_span=grid.n_span,
                n_endcap=grid.n_endcap, n_panels, h_proxy, angle, state="semiinf",
                wake_length_c=Inf, mesh_valid_for_fit), pair, (; sigma_min)))

            @printf("  angle=%5.2f°  resid ratio (unreg/reg): flat0.5c=%.4f  flat0.05c=%.4f  semiinf=%.4f\n",
                angle,
                rows[end-2].residual_ratio_unreg_over_reg,
                rows[end-1].residual_ratio_unreg_over_reg,
                rows[end].residual_ratio_unreg_over_reg)
        end
        template = nothing; solver = nothing; green = nothing; areas = nothing
        GC.gc()
    end

    csv_path = joinpath(path, "green_identity_regularization_refinement.csv")
    _write_rows_csv(csv_path, rows)

    # log-log slopes vs h=1/√N for regularized and unregularized residuals and
    # trace gaps, first three valid grids only (gate 7).
    println("\nRegularization refinement — log-log slope vs h (valid grids only)")
    @printf("  %-12s %-7s %10s %10s %10s %10s\n", "state", "angle",
        "slope_r_reg", "slope_r_unr", "slope_q_reg", "slope_q_unr")
    for state in ("flat_da0.5c", "flat_da0.05c", "semiinf"),
            angle in alpha_degs
        sel = filter(r -> r.state == state && r.angle == angle &&
            r.mesh_valid_for_fit, rows)
        length(sel) >= 2 || continue
        h = [r.h_proxy for r in sel]
        srr = _loglog_slope(h, [r.regularized_residual_rel for r in sel])
        sru = _loglog_slope(h, [r.unregularized_residual_rel for r in sel])
        sqr = _loglog_slope(h, [r.regularized_q_recon_vs_direct_rel for r in sel])
        squ = _loglog_slope(h, [r.unregularized_q_recon_vs_direct_rel for r in sel])
        @printf("  %-12s %5.2f° %10.3f %10.3f %10.3f %10.3f\n",
            state, angle, srr, sru, sqr, squ)
    end
    println("\n  Paired residual ratio (unreg/reg) at each resolution:")
    for state in ("flat_da0.5c", "flat_da0.05c", "semiinf"),
            angle in alpha_degs
        sel = filter(r -> r.state == state && r.angle == angle, rows)
        isempty(sel) && continue
        rs = join((@sprintf("N=%d:%.3f", r.n_panels,
            r.residual_ratio_unreg_over_reg) for r in sel), "  ")
        @printf("  %-12s %5.2f°  %s\n", state, angle, rs)
    end
    println("  wrote $(abspath(csv_path))")
    return (; rows, csv_path, blas_threads)
end

"Gate 3 for refinement, keyed by (state, angle, n_panels)."
function _check_regularized_baseline_refine(baseline, key, value; label)
    if isempty(baseline)
        return nothing
    elseif haskey(baseline, key)
        isapprox(value, baseline[key]; rtol=1e-10, atol=1e-13) || error(
            "$label regularized residual $value does not reproduce baseline " *
            "$(baseline[key]) for $key")
    else
        @warn "$label: no baseline refinement row for $key (gate 3 skipped for it)"
    end
    return nothing
end

"Right-singular-vector character of the smallest ~10 modes of (I−B) at
production resolution: projection on the constant and trailing-edge localization."
function _green_smallmode_rows(F, body)
    N = body.ncells
    edges = pnl._shedding_edge_map(body)
    te_mask = falses(N)
    for u in edges.upper; u > 0 && (te_mask[u] = true); end
    for l in edges.lower; l > 0 && (te_mask[l] = true); end
    ones_n = ones(Float64, N) ./ sqrt(N)
    rows = NamedTuple[]
    for k in 1:min(10, N)
        col = N - k + 1                     # k-th smallest singular value
        v = view(F.V, :, col)
        te = sum(abs2, view(v, te_mask)) / max(sum(abs2, v), eps(Float64))
        push!(rows, (; index=k, sigma=F.S[col],
            const_projection=abs(dot(v, ones_n)),
            te_localization_fraction=te))
    end
    return rows
end

"Parse a flat header/row CSV into a Vector{Dict{String,String}} (generic)."
function _read_simple_csv(path)
    isfile(path) || return NamedTuple[]
    lines = readlines(path)
    length(lines) >= 2 || return Dict{String,String}[]
    header = String.(split(first(lines), ','))
    out = Dict{String,String}[]
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        vals = split(line, ','; keepempty=true)
        push!(out, Dict(zip(header, String.(vals))))
    end
    return out
end

function run_green_identity_conditioning(; path=DIRICHLET_DATA_PATH, max_panels=Inf)
    mkpath(path)
    blas_threads = _task4_pin_blas()
    c = 1.0 * FT_TO_M; b = 4.0 * c; U = 330.2 * FT_TO_M
    Uinf = _uinf_from_alpha(U, DEFAULT_ALPHA_DEG)
    rows = NamedTuple[]
    sigma_2_production = NaN
    production_N = 0

    for grid in _GREEN_REFINEMENT_GRIDS
        body = _refinement_finite_body(c, b, Uinf, grid; das_length_c=0.5)
        N = body.ncells
        if N > max_panels
            @printf("[conditioning] skip n_panels=%d > max_panels=%g\n", N, max_panels)
            body = nothing; GC.gc()
            continue
        end
        B = Matrix{Float64}(undef, N, N)
        pnl._assemble_B!(B, body)
        ImB = -B
        for i in 1:N; ImB[i, i] += 1.0; end
        B = nothing

        is_production = grid.n_airfoil == 161
        local F
        if is_production
            F = svd(ImB)
            s = F.S
        else
            s = svdvals(ImB)
        end
        smax = s[1]
        small = [s[N - k + 1] for k in 1:min(10, N)]   # smallest first
        while length(small) < 10; push!(small, NaN); end
        sigma_min = small[1]; sigma_2 = small[2]
        cond_gf = smax / sigma_2
        const_check = norm(ImB * ones(Float64, N)) / sqrt(N)
        @printf("\n[conditioning] n_panels=%d: sigma_max=%.4e sigma_min=%.4e sigma_2=%.4e\n",
            N, smax, sigma_min, sigma_2)
        @printf("               cond_gaugefixed=%.4e  #(<1e-3σmax)=%d  #(<1e-6σmax)=%d  const_check=%.4e\n",
            cond_gf, count(<(1e-3 * smax), s), count(<(1e-6 * smax), s), const_check)

        push!(rows, (; n_panels=N, sigma_max=smax, sigma_min=sigma_min,
            sigma_2=small[2], sigma_3=small[3], sigma_4=small[4], sigma_5=small[5],
            sigma_6=small[6], sigma_7=small[7], sigma_8=small[8], sigma_9=small[9],
            sigma_10=small[10],
            cond_gaugefixed=cond_gf,
            n_below_1e_3=count(<(1e-3 * smax), s),
            n_below_1e_6=count(<(1e-6 * smax), s),
            const_mode_check=const_check))

        if is_production
            sigma_2_production = sigma_2
            production_N = N
            smrows = _green_smallmode_rows(F, body)
            _write_rows_csv(joinpath(path, "green_identity_smallmodes.csv"), smrows)
            println("  smallest-mode character (production N=$N):")
            @printf("    %5s %12s %14s %18s\n",
                "index", "sigma", "const_proj", "te_localization")
            for r in smrows
                @printf("    %5d %12.4e %14.4e %18.4e\n",
                    r.index, r.sigma, r.const_projection, r.te_localization_fraction)
            end
            F = nothing
        end
        ImB = nothing; s = nothing; body = nothing; GC.gc()
    end

    csv_path = joinpath(path, "green_identity_conditioning.csv")
    _write_rows_csv(csv_path, rows)

    # Synthesis: predicted trace gap ‖r‖/σ₂ vs observed q_recon_vs_direct.
    if isfinite(sigma_2_production)
        println("\nPredicted trace gap  ‖r‖/σ₂  vs observed q_recon_vs_direct  (σ₂ at N=$production_N)")
        @printf("  %-8s %-14s %12s %12s %10s\n",
            "angle", "state", "predicted", "observed", "ratio")
        for (label, resid_csv) in (
                ("3.94", joinpath(path, "green_identity_residual.csv")),
                ("45", joinpath(path, "alpha45", "green_identity_residual.csv")))
            for row in _read_simple_csv(resid_csv)
                haskey(row, "green_identity_residual_rel") || continue
                haskey(row, "Ssigma_l2") && haskey(row, "q_direct_l2") || continue
                resid = tryparse(Float64, row["green_identity_residual_rel"])
                ssig = tryparse(Float64, row["Ssigma_l2"])
                qdir = tryparse(Float64, row["q_direct_l2"])
                obs = tryparse(Float64, get(row, "q_recon_vs_direct_rel", ""))
                (isnothing(resid) || isnothing(ssig) || isnothing(qdir) ||
                    isnothing(obs) || qdir == 0) && continue
                predicted = (resid * ssig / sigma_2_production) / qdir
                state = get(row, "notes", get(row, "wake_state", "?"))
                @printf("  %-8s %-14s %12.4e %12.4e %10.3f\n",
                    label, state, predicted, obs, predicted / max(obs, eps()))
            end
        end
        println("  (ratio≈1 ⇒ trace gap = discretization residual × (I−B) conditioning)")
    end
    println("  wrote $(abspath(csv_path))")
    return (; rows, csv_path, sigma_2_production, production_N, blas_threads)
end

# -----------------------------------------------------------------------------
# Task 5: velocity-only line-integral trace correction
# -----------------------------------------------------------------------------

const TASK5_VARIANTS = Dict(
    "trap" => (quadrature=:trapezoid, interior_path=:straight, path_depth=1.0),
    "simpson" => (quadrature=:simpson, interior_path=:straight, path_depth=1.0),
    "deformed" => (quadrature=:simpson, interior_path=:deformed, path_depth=1.0),
)

"Task-5 mutation snapshot. Body strengths and the affine correction channel are
checked separately because those are the only fields an eligible solve may alter."
function _task5_snapshot(body, wake)
    return (;
        body_nodes=copy(body.nodes), body_cells=copy(body.cells),
        body_shedding=[copy(x) for x in body.shedding],
        body_Das=[copy(x) for x in body.Das],
        body_strength=copy(body.strength),
        body_correction=[copy(x) for x in body.wake_strength_correction],
        body_shift=copy(body.wake_strength_shift),
        correction_active=body.wake_correction_active[],
        wake_nodes=[copy(x) for x in wake.nodes],
        wake_strength=[copy(x) for x in wake.strength],
        wake_node_shapes=[size(x) for x in wake.nodes],
        wake_strength_shapes=[size(x) for x in wake.strength],
        wake_rows=wake.nwakes[],
        wake_options=(; core_size=wake.core_size,
            shed_with_induced_velocity=wake.shed_with_induced_velocity,
            unsteady_filament=wake.unsteady_filament,
            include_final_filament=wake.include_final_filament),
    )
end

function _task5_assert_solve_mutations(snapshot, body, wake, state; label)
    isequal(body.nodes, snapshot.body_nodes) || error("$label changed body nodes")
    isequal(body.cells, snapshot.body_cells) || error("$label changed body cells")
    all(isequal(a, b) for (a, b) in zip(body.shedding, snapshot.body_shedding)) ||
        error("$label changed shedding maps")
    all(isequal(a, b) for (a, b) in zip(body.Das, snapshot.body_Das)) ||
        error("$label changed attached-wake geometry")
    [size(x) for x in wake.nodes] == snapshot.wake_node_shapes ||
        error("$label changed wake-node storage shape")
    [size(x) for x in wake.strength] == snapshot.wake_strength_shapes ||
        error("$label changed wake-strength storage shape")
    wake.nwakes[] == snapshot.wake_rows || error("$label changed wake row count")
    all(isequal(a, b) for (a, b) in zip(wake.nodes, snapshot.wake_nodes)) ||
        error("$label changed wake geometry")
    all(isequal(a, b) for (a, b) in zip(wake.strength, snapshot.wake_strength)) ||
        error("$label changed active or inactive wake strengths")
    options = (; core_size=wake.core_size,
        shed_with_induced_velocity=wake.shed_with_induced_velocity,
        unsteady_filament=wake.unsteady_filament,
        include_final_filament=wake.include_final_filament)
    isequal(options, snapshot.wake_options) || error("$label changed wake options")
    body.wake_correction_active[] || error("$label did not activate physical correction")
    expected = deepcopy(body)
    pnl.clear_wake_correction!(expected)
    pnl.set_wake_correction!(expected, state.c)
    all(isequal(a, b) for (a, b) in zip(body.wake_strength_correction,
        expected.wake_strength_correction)) || error("$label changed an unexpected c channel")
    isequal(body.wake_strength_shift, expected.wake_strength_shift) ||
        error("$label affine panel shift is inconsistent with c")
    return true
end

function _task5_snapshot_dict(snapshot)
    Dict(
        "body_nodes_sha256" => _array_hash(snapshot.body_nodes),
        "body_cells_sha256" => _array_hash(snapshot.body_cells),
        "body_shedding_sha256" => [_array_hash(x) for x in snapshot.body_shedding],
        "body_Das_sha256" => [_array_hash(x) for x in snapshot.body_Das],
        "body_strength_before_sha256" => _array_hash(snapshot.body_strength),
        "body_correction_before_sha256" => [_array_hash(x) for x in snapshot.body_correction],
        "body_shift_before_sha256" => _array_hash(snapshot.body_shift),
        "correction_active_before" => snapshot.correction_active,
        "wake_nodes_sha256" => [_array_hash(x) for x in snapshot.wake_nodes],
        "wake_strength_sha256" => [_array_hash(x) for x in snapshot.wake_strength],
        "wake_node_shapes" => [collect(x) for x in snapshot.wake_node_shapes],
        "wake_strength_shapes" => [collect(x) for x in snapshot.wake_strength_shapes],
        "wake_rows" => snapshot.wake_rows,
        "wake_options" => Dict(string(k) => v for (k, v) in pairs(snapshot.wake_options)),
    )
end

function _task5_formulation(spec; path_depth=spec.path_depth)
    pnl.TraceCorrected(; estimator=:line_integral, quadrature=spec.quadrature,
        interior_path=spec.interior_path, path_depth=path_depth,
        recompute_interval=1)
end

function _task5_physical_gamma(body)
    [pnl._get_wakestrength_Gamma(body, i, s)
        for (s, shed) in enumerate(body.shedding) for i in axes(shed, 2)]
end

"One eligible velocity-only solve, followed by physical diagnostics and only
then the direct scalar-potential oracle. A fresh TraceCorrectedState is created
for every call."
function _task5_solve(master_body, master_wake, solver, C, formulation, Uinf,
        rho, Sref, chord, dt, probes, task1_probe, task2_probe; label,
        i_step=0, compute_oracle=true, backend=pnl.DirectBackend())
    body, wake = deepcopy(master_body), deepcopy(master_wake)
    snapshot = _task5_snapshot(body, wake)
    state = nothing
    init_elapsed = @elapsed state = pnl.initialize_formulation(formulation,
        (body,), (wake,), solver, backend, backend)
    frames = pnl.ReferenceFrame(body)
    solve_elapsed = @elapsed pnl._steady_aerodynamics!(body, (body,), (wake,),
        frames, Uinf, solver; backend_wake=backend, backend_solve=backend,
        backend_system=backend, update_trailing_edges=false, formulation,
        formulation_state=state, i_step)
    _task5_assert_solve_mutations(snapshot, body, wake, state; label)

    mu = copy(view(body.strength, :, 2))
    sigma = copy(view(body.strength, :, 1))
    c_est = copy(state.c)
    gamma = _task5_physical_gamma(body)
    Ssigma = zeros(eltype(mu), length(mu))
    pnl._source_potential!(Ssigma, body, sigma, backend)
    residual = _lu_matrix_product(solver.Glu, mu) + Ssigma - state.W * c_est
    rhsnorm = max(norm(-Ssigma + state.W * c_est), eps(Float64))
    residual_abs, residual_rel = norm(residual), norm(residual) / rhsnorm
    gamma_error = gamma - (C * mu - c_est)
    gamma_consistency_abs = norm(gamma_error)
    gamma_consistency_rel = gamma_consistency_abs /
        max(norm(C * mu - c_est), eps(Float64))
    residual_rel <= 1e-10 || error("$label formulation residual $residual_rel exceeds 1e-10")
    gamma_consistency_rel <= 1e-10 || error(
        "$label circulation consistency $gamma_consistency_rel exceeds 1e-10")

    # Every downstream physical diagnostic runs while the affine correction is active.
    CL, pressure, force = _monitor_lift!(body, wake, frames, Uinf, rho, Sref,
        chord, dt, backend)
    velocity = _exterior_velocity(body, wake, Uinf, probes, backend)
    probe_rel_task1 = norm(velocity - task1_probe) / max(norm(task1_probe), eps(Float64))
    probe_rel_task2 = norm(velocity - task2_probe) / max(norm(task2_probe), eps(Float64))
    first_free = vec(copy(view(wake.strength[1], 1, 1, :)))
    mismatch_l2 = norm(gamma - first_free)
    mismatch_rel = mismatch_l2 / max(norm(first_free), eps(Float64))

    # Diagnostic oracle: deliberately outside the eligible solve and physical
    # velocity/pressure/probe path, on the exact same frozen geometry/strengths.
    oracle_elapsed = 0.0
    q_f = Float64[]
    c_oracle = Float64[]
    if compute_oracle
        oracle_elapsed = @elapsed begin
            q_f = _task3_accumulated_potential(body, (wake,))
            c_oracle = C * q_f
        end
    end
    c_error_l2 = compute_oracle ? norm(c_est - c_oracle) : 0.0
    c_error_rel = compute_oracle ? c_error_l2 / max(norm(c_oracle), eps(Float64)) : 0.0
    c_error_max = compute_oracle ? maximum(abs, c_est - c_oracle) : 0.0
    c_oracle_l2 = compute_oracle ? norm(c_oracle) : 0.0
    all(isfinite, body.strength) && all(isfinite, gamma) && all(isfinite, c_est) &&
        all(isfinite, (CL, residual_abs, residual_rel, gamma_consistency_abs,
            gamma_consistency_rel, mismatch_l2, mismatch_rel, probe_rel_task1,
            probe_rel_task2, init_elapsed, solve_elapsed, oracle_elapsed,
            c_error_l2, c_error_rel, c_error_max, c_oracle_l2)) ||
        error("$label produced a non-finite field")
    _task5_assert_solve_mutations(snapshot, body, wake, state;
        label="$label complete diagnostic")
    return (; body, wake, snapshot, state, pressure, force, mu, sigma, gamma,
        c=c_est, q_f, c_oracle, CL, residual_abs, residual_rel,
        gamma_consistency_abs, gamma_consistency_rel, mismatch_l2, mismatch_rel,
        probe_velocity=velocity, probe_rel_task1, probe_rel_task2,
        c_error_l2, c_error_rel, c_error_max, c_l2=norm(c_est), c_oracle_l2,
        init_elapsed, solve_elapsed, oracle_elapsed,
        total_elapsed=init_elapsed + solve_elapsed + oracle_elapsed)
end

function _task5_iteration(master_body, master_wake, solver, C, formulation,
        Uinf, rho, Sref, chord, dt, probes, task1_probe, task2_probe; label,
        max_iterations=200, backend=pnl.DirectBackend())
    body, wake = deepcopy(master_body), deepcopy(master_wake)
    geometry = _task5_snapshot(body, wake)
    inactive_hashes = _task3_inactive_strength_hashes(wake)
    rows = NamedTuple[]
    omega, previous_defect, previous_CL = 1.0, Inf, 0.0
    insufficient_decrease = consecutive = 0
    converged = false
    for iteration in 1:max_iterations
        active_before = _task3_active_strengths(wake)
        result = _task5_solve(body, wake, solver, C, formulation, Uinf, rho,
            Sref, chord, dt, probes, task1_probe, task2_probe;
            label="$label iteration $iteration", i_step=iteration,
            compute_oracle=false, backend)
        defect = maximum(abs.(active_before .- reshape(result.gamma, 1, :)))
        lift_change = iteration == 1 ? 0.0 : abs(result.CL - previous_CL) /
            max(abs(result.CL), abs(previous_CL), eps(Float64))
        lift_change_applicable = iteration > 1
        if defect <= 1e-8 && lift_change_applicable && lift_change <= 1e-8
            consecutive += 1
        else
            consecutive = 0
        end
        push!(rows, (; iteration, omega, max_rowwise_strength_defect=defect,
            relative_lift_change=lift_change, lift_change_applicable,
            consecutive_converged=consecutive,
            formulation_residual_rel=result.residual_rel,
            circulation_consistency_rel=result.gamma_consistency_rel,
            CL=result.CL, transition_gamma_min=minimum(result.gamma),
            transition_gamma_max=maximum(result.gamma),
            transition_gamma_l2=norm(result.gamma),
            transition_free_mismatch_l2=result.mismatch_l2,
            transition_free_mismatch_rel=result.mismatch_rel,
            c_l2=result.c_l2, eligible_solve_elapsed_s=result.solve_elapsed,
            frozen_inner_solve=true))
        @printf("%s iteration %d: omega=%.4g, defect=%.3e, dCL=%.3e, streak=%d\n",
            label, iteration, omega, defect, lift_change, consecutive)
        flush(stdout)
        body, wake = result.body, result.wake
        if consecutive >= 3
            converged = true
            break
        end
        insufficient_decrease = isfinite(previous_defect) &&
            defect > 0.95previous_defect ? insufficient_decrease + 1 : 0
        if insufficient_decrease >= 2
            omega = max(omega / 2, 1 / 16)
            insufficient_decrease = 0
        end
        eta = (1 - omega) .* active_before .+ omega .* reshape(result.gamma, 1, :)
        _task3_set_active_strengths!(wake, eta)
        isequal(body.nodes, geometry.body_nodes) || error("$label changed body geometry")
        all(isequal(a, b) for (a, b) in zip(wake.nodes, geometry.wake_nodes)) ||
            error("$label changed wake geometry")
        isequal(_task3_inactive_strength_hashes(wake), inactive_hashes) ||
            error("$label changed inactive wake strengths")
        previous_defect, previous_CL = defect, result.CL
    end
    status = converged ? "converged" : "nonconverged_at_$max_iterations"
    result = _task5_solve(body, wake, solver, C, formulation, Uinf, rho, Sref,
        chord, dt, probes, task1_probe, task2_probe;
        label="$label terminal diagnostic", i_step=length(rows)+1,
        compute_oracle=true, backend)
    return (; rows, result, converged, status, omega, inactive_hashes, geometry)
end

function _task5_route_row(route, result, iterations, baseline_CL, task2_CL,
        task3_CL, closures, spec, path_depth)
    (; route, status=route == "single" ? "single_shot" : iterations.status,
        quadrature=string(spec.quadrature), interior_path=string(spec.interior_path),
        path_depth, CL=result.CL, CL_minus_task1=result.CL-baseline_CL,
        relative_CL_difference=_relative_difference(result.CL, baseline_CL),
        CL_minus_task2=result.CL-task2_CL, CL_minus_task3=result.CL-task3_CL,
        f_semiinf=closures.f_semiinf, f_oracle=closures.f_oracle,
        signed_post_correction_diff=closures.signed_post_correction_diff,
        c_error_rel=result.c_error_rel, c_error_l2=result.c_error_l2,
        c_error_max=result.c_error_max, c_l2=result.c_l2,
        c_oracle_l2=result.c_oracle_l2,
        formulation_residual_rel=result.residual_rel,
        circulation_consistency_rel=result.gamma_consistency_rel,
        transition_free_mismatch_rel=result.mismatch_rel,
        relative_exterior_velocity_difference_task1=result.probe_rel_task1,
        relative_exterior_velocity_difference_task2=result.probe_rel_task2,
        fixed_point_iterations=route == "single" ? 0 : length(iterations.rows),
        state_init_elapsed_s=result.init_elapsed,
        eligible_solve_elapsed_s=result.solve_elapsed,
        oracle_elapsed_s=result.oracle_elapsed, total_elapsed_s=result.total_elapsed)
end

function _task5_comparison!(path, case, variant, route, status, result, spec,
        path_depth, frozen_step, wake_length_c, baseline_CL, task2_CL,
        task3_single_CL, task3_iter_CL, iterations, blas_threads)
    task3_CL = route == "single" ? task3_single_CL : task3_iter_CL
    closures = _task4_closures(result.CL, baseline_CL, task2_CL, task3_CL)
    formulation_id = "trace_line_$(variant)_$(route)"
    _upsert_comparison(joinpath(path, "comparison.csv"), (;
        task=5, case, formulation=formulation_id, status,
        frozen_wake_step=frozen_step, wake_length_c,
        wake_rows=result.wake.nwakes[], ncells=result.body.ncells,
        nshedding_edges=length(result.gamma),
        linear_residual_abs=result.residual_abs,
        linear_residual_rel=result.residual_rel, kutta_constant_error=0.0,
        gamma_min=minimum(result.gamma), gamma_max=maximum(result.gamma),
        gamma_mean=sum(result.gamma)/length(result.gamma), gamma_l2=norm(result.gamma),
        transition_free_mismatch_l2=result.mismatch_l2,
        transition_free_mismatch_rel=result.mismatch_rel,
        CL=result.CL, CL_minus_task1=result.CL-baseline_CL,
        relative_CL_difference=_relative_difference(result.CL, baseline_CL),
        relative_exterior_velocity_difference=result.probe_rel_task1,
        CL_minus_task2=result.CL-task2_CL,
        relative_CL_difference_task2=_relative_difference(result.CL, task2_CL),
        relative_exterior_velocity_difference_task2=result.probe_rel_task2,
        fixed_point_iterations=route == "single" ? 0 : length(iterations.rows),
        CL_minus_task3_single=result.CL-task3_single_CL,
        CL_minus_task3_iterated=result.CL-task3_iter_CL,
        signed_post_correction_diff=closures.signed_post_correction_diff,
        f_semiinf=closures.f_semiinf, f_oracle=closures.f_oracle,
        explicit_potential_residual_rel=result.residual_rel, blas_threads,
        estimator="line_integral", quadrature=string(spec.quadrature),
        interior_path=string(spec.interior_path), path_depth,
        c_error_rel=result.c_error_rel, c_error_l2=result.c_error_l2,
        c_error_max=result.c_error_max, c_l2=result.c_l2,
        c_oracle_l2=result.c_oracle_l2,
        circulation_consistency_abs=result.gamma_consistency_abs,
        circulation_consistency_rel=result.gamma_consistency_rel,
        state_init_elapsed_s=result.init_elapsed,
        eligible_solve_elapsed_s=result.solve_elapsed,
        oracle_elapsed_s=result.oracle_elapsed, total_elapsed_s=result.total_elapsed))
    return closures
end

function _task5_write_c_csv(path, result)
    rows = [(; edge=i, c=result.c[i], c_oracle=result.c_oracle[i],
        error=result.c[i]-result.c_oracle[i]) for i in eachindex(result.c)]
    _write_rows_csv(path, rows)
end

function _run_task5_route(path, case, variant, route, t3, p;
        path_depth=TASK5_VARIANTS[variant].path_depth,
        artifact_variant=variant, backend=pnl.DirectBackend(),
        blas_threads=BLAS.get_num_threads())
    spec = TASK5_VARIANTS[variant]
    formulation = _task5_formulation(spec; path_depth)
    template, solver, C = t3.finite_setup
    task2 = case == "march" ? t3.task2.result : t3.task2
    seed_body, seed_wake = task2.body, task2.wake
    baseline = t3.baseline
    frozen_step = case == "march" ? t3.task2.selected_step : -1
    wake_length_c = case == "march" ? 0.5frozen_step :
        norm(view(template.Das[1], :, 1))/p.c + task2.wake.nwakes[]*0.5
    stem = "task5_$(artifact_variant)_$(case)_$(route)"
    iterations = nothing
    if route == "single"
        result = _task5_solve(seed_body, seed_wake, solver, C, formulation,
            p.Uinf, p.rho, p.Sref, p.c, p.dt, baseline.probes,
            baseline.probe_velocity, task2.probe_velocity;
            label="Task 5 $variant $case single", backend)
        status = "single_shot"
        iterations = (; rows=NamedTuple[], status)
    else
        iterations = _task5_iteration(seed_body, seed_wake, solver, C,
            formulation, p.Uinf, p.rho, p.Sref, p.c, p.dt, baseline.probes,
            baseline.probe_velocity, task2.probe_velocity;
            label="Task 5 $variant $case iterated", backend)
        result, status = iterations.result, iterations.status
    end
    task3_single_CL, task3_iter_CL = t3.result.CL, t3.iteration.result.CL
    closures = _task5_comparison!(path, case, artifact_variant, route, status, result,
        spec, path_depth, frozen_step, wake_length_c, baseline.CL, task2.CL,
        task3_single_CL, task3_iter_CL, iterations, blas_threads)
    task3_CL = route == "single" ? task3_single_CL : task3_iter_CL
    route_row = _task5_route_row(route, result, iterations, baseline.CL,
        task2.CL, task3_CL, closures, spec, path_depth)
    _write_rows_csv(joinpath(path, stem * "_route.csv"), [route_row])
    history = route == "single" ? [(; iteration=0, omega=1.0,
        max_rowwise_strength_defect=0.0, relative_lift_change=0.0,
        lift_change_applicable=false, consecutive_converged=0,
        formulation_residual_rel=result.residual_rel,
        circulation_consistency_rel=result.gamma_consistency_rel, CL=result.CL,
        transition_gamma_min=minimum(result.gamma),
        transition_gamma_max=maximum(result.gamma),
        transition_gamma_l2=norm(result.gamma),
        transition_free_mismatch_l2=result.mismatch_l2,
        transition_free_mismatch_rel=result.mismatch_rel, c_l2=result.c_l2,
        eligible_solve_elapsed_s=result.solve_elapsed, frozen_inner_solve=true)] :
        iterations.rows
    _write_rows_csv(joinpath(path, stem * "_iteration.csv"), history)
    _task5_write_c_csv(joinpath(path, stem * "_trace.csv"), result)
    open(joinpath(path, stem * "_frozen_state.toml"), "w") do io
        TOML.print(io, Dict("snapshot" => _task5_snapshot_dict(result.snapshot),
            "solved_body_strength_sha256" => _array_hash(result.body.strength),
            "active_wake_strength_sha256" => _array_hash(_task3_active_strengths(result.wake)),
            "inactive_wake_strength_sha256" => _task3_inactive_strength_hashes(result.wake),
            "c_sha256" => _array_hash(result.c),
            "correction_shift_sha256" => _array_hash(result.body.wake_strength_shift),
            "all_mutation_gates_passed" => true); sorted=true)
    end
    terminal_time = wake_length_c * p.c / p.U
    pnl.write_vtk(joinpath(path, "dirichlet_$(stem)_body"), result.body,
        max(frozen_step, 0), terminal_time;
        monitors=(result.pressure, result.force), overwrite=true)
    pnl.write_vtk(joinpath(path, "dirichlet_$(stem)_wake"), result.wake,
        max(frozen_step, 0), terminal_time; overwrite=true)
    @printf("Task 5 %s %s %s: c_rel=%.4e, CL=%.10f, f_oracle=%.4f, residual=%.2e\n",
        variant, case, route, result.c_error_rel, result.CL, closures.f_oracle,
        result.residual_rel)
    return (; result, iterations, closures, route_row, frozen_step, wake_length_c)
end

_task5_parse_float(row, key) = parse(Float64, row[key])

function _task5_stored_baseline(path, alpha_deg)
    rows = _read_comparison_rows(joinpath(path, "comparison.csv"))
    key = "semi_infinite_attached_wake|baseline"
    haskey(rows, key) || return _task2_baseline(path; alpha_deg)
    probe_rows = _read_simple_csv(joinpath(path, "task1_exterior_probes.csv"))
    isempty(probe_rows) && return _task2_baseline(path; alpha_deg)
    probes = [SVector{3,Float64}(parse(Float64, r["x_m"]),
        parse(Float64, r["y_m"]), parse(Float64, r["z_m"])) for r in probe_rows]
    velocity = hcat(([parse(Float64, r["u_m_s"]), parse(Float64, r["v_m_s"]),
        parse(Float64, r["w_m_s"])] for r in probe_rows)...)
    return (; CL=_task5_parse_float(rows[key], :CL), probes,
        probe_velocity=velocity)
end

function _task5_last_history_CL(path, case)
    rows = _read_simple_csv(joinpath(path, "task3_$(case)_iteration_history.csv"))
    isempty(rows) && error("missing stored Task 3 iteration history for $case")
    return parse(Float64, last(rows)["CL"])
end

"Load an already-audited Task-2 terminal VTK wake and reproduce its physical
control solve. This is the same read-only reconstruction route used by the
Green-identity audit and avoids re-marching/re-projecting Tasks 1-3."
function _task5_load_stored_seed(path, case, alpha_deg)
    comparison = _read_comparison_rows(joinpath(path, "comparison.csv"))
    task2_form = case == "flat" ? "finite_flat_velocity_through_sources" :
        case == "flat_das005" ? "finite_flat_das005_velocity_through_sources" :
        case == "march" ? "finite_rolled_velocity_through_sources" :
        error("unknown Task 5 wake state $case")
    task3_form = case == "flat" ? "finite_flat_direct_fixed_wake_potential" :
        case == "flat_das005" ? "finite_flat_das005_direct_fixed_wake_potential" :
        "finite_rolled_direct_fixed_wake_potential"
    t2key, t3key = "$task2_form|$case", "$task3_form|$case"
    haskey(comparison, t2key) || error("missing stored Task 2 comparison row $t2key")
    haskey(comparison, t3key) || error("missing stored Task 3 comparison row $t3key")
    p = _task3_parameters(alpha_deg)
    das = case == "flat" ? 0.5 : 0.05
    template, solver, C = _finite_solver_template(p.c, p.b, p.Uinf;
        das_length_c=das)
    idx = case == "march" ? parse(Int, comparison[t2key][:frozen_wake_step]) : 0
    structure = case == "march" ? :rolled : :flat
    stem = "dirichlet_task2_$(case)_terminal_wake"
    wake = _green_id_load_wake(template, path, stem, idx, structure;
        nwakerows=200)
    baseline = _task5_stored_baseline(path, alpha_deg)
    task2 = _frozen_solve(template, wake, solver, C, p.Uinf, p.rho, p.Sref,
        p.c, p.dt, baseline.probes, baseline.probe_velocity;
        label="Task 5 stored $case Task 2 control", backend=pnl.DirectBackend())
    stored_t2_CL = _task5_parse_float(comparison[t2key], :CL)
    isapprox(task2.CL, stored_t2_CL; rtol=1e-8, atol=1e-10) || error(
        "Task 5 loaded $case control CL $(task2.CL) != stored $stored_t2_CL")
    t3_single = (; CL=_task5_parse_float(comparison[t3key], :CL))
    t3_iteration = (; result=(; CL=_task5_last_history_CL(path, case)))
    task2_field = case == "march" ? (; result=task2, selected_step=idx) : task2
    return (; result=t3_single, iteration=t3_iteration, task2=task2_field,
        baseline, finite_setup=(template, solver, C), seed_source="stored_terminal_vtk")
end

function _task5_build_seed(path, case, alpha_deg)
    stem = joinpath(path, "dirichlet_task2_$(case)_terminal_wake")
    if isdir(stem) && isfile(joinpath(path, "comparison.csv"))
        return _task5_load_stored_seed(path, case, alpha_deg)
    end
    # Reproducible fallback for a fresh output directory with no Task 1-3 data.
    baseline = _task2_baseline(path; alpha_deg)
    case == "flat" && return run_task3_flat(; path, baseline, alpha_deg)
    case == "flat_das005" && return run_task3_flat_das005(; path, baseline, alpha_deg)
    case == "march" && return run_task3_march(; path, baseline, alpha_deg)
    error("unknown Task 5 wake state $case")
end

function _write_task5_manifests(path; alpha_deg=DEFAULT_ALPHA_DEG,
        blas_threads=BLAS.get_num_threads())
    mkpath(path)
    config = Dict("task" => 5, "formulation" => "TraceCorrected",
        "estimator" => "line_integral", "alpha_deg" => alpha_deg,
        "variants" => Dict(k => Dict("quadrature" => string(v.quadrature),
            "interior_path" => string(v.interior_path),
            "path_depth" => v.path_depth) for (k, v) in TASK5_VARIANTS),
        "backend" => Dict("wake" => "DirectBackend", "solve" => "DirectBackend",
            "physical_system" => "DirectBackend"),
        "iteration" => Dict("defect_tolerance" => 1e-8,
            "relative_lift_tolerance" => 1e-8, "required_consecutive" => 3,
            "maximum_iterations" => 200, "minimum_relaxation" => 0.0625),
        "residual_tolerance" => 1e-10, "blas_threads" => blas_threads,
        "oracle" => "direct wake scalar potential after eligible solve only")
    open(joinpath(path, "task5.config.toml"), "w") do io
        TOML.print(io, config; sorted=true)
    end
    artifacts = sort(filter(x -> startswith(x, "task5_") ||
        startswith(x, "dirichlet_task5_"), readdir(path)))
    metadata = Dict("comparison_csv" => "comparison.csv",
        "path_depth_csv" => "task5_line_integral_path_depth.csv",
        "artifacts" => artifacts, "alpha_deg" => alpha_deg,
        "theory_audit" => "no correction required",
        "eligible_path" => "wake velocity only; scalar-potential oracle post-solve")
    open(joinpath(path, "task5.metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    invariant_files = filter(x -> startswith(x, "task5_") &&
        endswith(x, "_frozen_state.toml"), readdir(path))
    invariants = Dict("all_recorded_fields_finite" => true,
        "formulation_residuals_within_1e-10" => true,
        "circulation_consistency_within_1e-10" => true,
        "fresh_state_per_independent_solve" => true,
        "frozen_geometry_and_storage_exact" => true,
        "only_expected_affine_c_channel_changed" => true,
        "no_free_wake_scalar_potential_during_eligible_solve" => true,
        "snapshot_files" => sort(invariant_files))
    open(joinpath(path, "task5_invariants.toml"), "w") do io
        TOML.print(io, invariants; sorted=true)
    end
    return nothing
end

function _task5_spread_rows(results, alpha_deg)
    rows = NamedTuple[]
    for case in ("flat", "flat_das005", "march"), route in ("single", "iterated")
        all(haskey(results, (v, case, route)) for v in keys(TASK5_VARIANTS)) || continue
        trap = results[("trap", case, route)].result.c
        simp = results[("simpson", case, route)].result.c
        deform = results[("deformed", case, route)].result.c
        push!(rows, (; alpha_deg, case, route,
            quadrature_spread_l2=norm(trap-simp),
            quadrature_spread_rel=norm(trap-simp)/max(norm(simp), eps(Float64)),
            path_spread_l2=norm(simp-deform),
            path_spread_rel=norm(simp-deform)/max(norm(deform), eps(Float64))))
    end
    rows
end

function _task5_stored_trace(path, variant, case, route)
    rows = _read_simple_csv(joinpath(path,
        "task5_$(variant)_$(case)_$(route)_trace.csv"))
    isempty(rows) && return nothing
    return [parse(Float64, r["c"]) for r in rows]
end

function _task5_rewrite_comparison(path, rows)
    open(path, "w") do io
        println(io, join(_COMPARISON_COLUMNS, ','))
        ordered = sort(collect(values(rows)); by=r -> (
            tryparse(Int, get(r, :task, "")) === nothing ? typemax(Int) : parse(Int, r[:task]),
            get(r, :case, ""), get(r, :formulation, "")))
        for stored in ordered
            println(io, join((get(stored, column, "") for column in _COMPARISON_COLUMNS), ','))
        end
    end
end

"Finalize spread diagnostics from independently checkpointed selector results."
function finalize_task5(path; alpha_deg=DEFAULT_ALPHA_DEG)
    rows = NamedTuple[]
    comparison_path = joinpath(path, "comparison.csv")
    comparison = _read_comparison_rows(comparison_path)
    for case in ("flat", "flat_das005", "march"), route in ("single", "iterated")
        trap = _task5_stored_trace(path, "trap", case, route)
        simp = _task5_stored_trace(path, "simpson", case, route)
        deform = _task5_stored_trace(path, "deformed", case, route)
        any(isnothing, (trap, simp, deform)) && continue
        ql2 = norm(trap-simp)
        pl2 = norm(simp-deform)
        push!(rows, (; alpha_deg, case, route, quadrature_spread_l2=ql2,
            quadrature_spread_rel=ql2/max(norm(simp), eps(Float64)),
            path_spread_l2=pl2,
            path_spread_rel=pl2/max(norm(deform), eps(Float64))))
        for variant in ("trap", "simpson", "deformed")
            key = "trace_line_$(variant)_$(route)|$case"
            haskey(comparison, key) || continue
            comparison[key][:quadrature_spread_l2] = string(ql2)
            comparison[key][:path_spread_l2] = string(pl2)
        end
    end
    isempty(rows) || _write_rows_csv(joinpath(path,
        "task5_line_integral_spreads.csv"), rows)
    _task5_rewrite_comparison(comparison_path, comparison)
    _write_task5_manifests(path; alpha_deg)
    return rows
end

function run_task5_audit(; path=DIRICHLET_DATA_PATH,
        alpha45_path=joinpath(path, "alpha45"))
    finalize_task5(path; alpha_deg=3.94)
    finalize_task5(alpha45_path; alpha_deg=45.0)
    roots = ((3.94, path), (45.0, alpha45_path))
    ranking = NamedTuple[]
    all_rows = NamedTuple[]
    for variant in ("trap", "simpson", "deformed")
        errors = Float64[]
        for (alpha, root) in roots, route in ("single", "iterated")
            stored = _read_comparison_rows(joinpath(root, "comparison.csv"))
            key = "trace_line_$(variant)_$(route)|march"
            haskey(stored, key) || error("Task 5 audit missing $alpha° $key")
            row = stored[key]
            err = parse(Float64, row[:c_error_rel])
            residual = parse(Float64, row[:linear_residual_rel])
            circulation = parse(Float64, row[:circulation_consistency_rel])
            all(isfinite, (err, residual, circulation)) || error(
                "Task 5 audit found non-finite metrics in $alpha° $key")
            residual <= 1e-10 || error("Task 5 audit residual gate failed in $alpha° $key")
            circulation <= 1e-10 || error(
                "Task 5 audit circulation gate failed in $alpha° $key")
            push!(errors, err)
            push!(all_rows, (; variant, alpha_deg=alpha, route,
                c_error_rel=err, formulation_residual_rel=residual,
                circulation_consistency_rel=circulation,
                CL=parse(Float64, row[:CL]), f_oracle=parse(Float64, row[:f_oracle])))
        end
        push!(ranking, (; variant, worst_primary_rolled_c_error_rel=maximum(errors),
            mean_primary_rolled_c_error_rel=mean(errors)))
    end
    sort!(ranking; by=r -> r.worst_primary_rolled_c_error_rel)
    ranked = [(; rank=i, r...) for (i, r) in enumerate(ranking)]
    _write_rows_csv(joinpath(path, "task5_line_integral_primary_cells.csv"), all_rows)
    _write_rows_csv(joinpath(path, "task5_line_integral_ranking.csv"), ranked)
    return (; ranking=ranked, primary_cells=all_rows)
end

function run_task5(; path=DIRICHLET_DATA_PATH, alpha_deg=DEFAULT_ALPHA_DEG)
    mkpath(path)
    blas_threads = _task4_pin_blas()
    p = _task3_parameters(alpha_deg)
    results = Dict{Tuple{String,String,String},Any}()
    for case in ("flat", "flat_das005", "march")
        t3 = _task5_build_seed(path, case, alpha_deg)
        for variant in ("trap", "simpson", "deformed"), route in ("single", "iterated")
            results[(variant, case, route)] = _run_task5_route(path, case,
                variant, route, t3, p; blas_threads)
        end
    end
    spreads = _task5_spread_rows(results, alpha_deg)
    _write_rows_csv(joinpath(path, "task5_line_integral_spreads.csv"), spreads)
    _write_task5_manifests(path; alpha_deg, blas_threads)
    return (; results, spreads, blas_threads)
end

function run_task5_selector(selector; path=DIRICHLET_DATA_PATH,
        alpha_deg=DEFAULT_ALPHA_DEG)
    m = match(r"^task5-(trap|simpson|deformed)-(flat|flat-das005|march)-(single|iterated)$",
        replace(selector, '_' => '-'))
    isnothing(m) && error("invalid Task 5 selector $selector")
    variant, case_raw, route = m.captures
    case = replace(case_raw, '-' => '_')
    blas_threads = _task4_pin_blas()
    t3 = _task5_build_seed(path, case, alpha_deg)
    result = _run_task5_route(path, case, variant, route, t3,
        _task3_parameters(alpha_deg); blas_threads)
    _write_task5_manifests(path; alpha_deg, blas_threads)
    return result
end

function run_task5_depth_march(; path=DIRICHLET_DATA_PATH,
        alpha_deg=DEFAULT_ALPHA_DEG)
    isapprox(alpha_deg, 3.94; atol=1e-12) || error(
        "Task 5 path-depth sensitivity is specified only for the rolled 3.94° state")
    mkpath(path)
    blas_threads = _task4_pin_blas()
    p = _task3_parameters(alpha_deg)
    t3 = _task5_build_seed(path, "march", alpha_deg)
    rows = NamedTuple[]
    for depth in (0.5, 1.0, 2.0), route in ("single", "iterated")
        if depth == 1.0
            stored = _read_comparison_rows(joinpath(path, "comparison.csv"))
            key = "trace_line_deformed_$(route)|march"
            if haskey(stored, key)
                sr = stored[key]
                push!(rows, (; alpha_deg, path_depth=depth, route,
                    c_error_rel=parse(Float64, sr[:c_error_rel]),
                    c_error_l2=parse(Float64, sr[:c_error_l2]),
                    c_error_max=parse(Float64, sr[:c_error_max]),
                    c_l2=parse(Float64, sr[:c_l2]),
                    c_oracle_l2=parse(Float64, sr[:c_oracle_l2]),
                    CL=parse(Float64, sr[:CL]),
                    f_oracle=parse(Float64, sr[:f_oracle]),
                    formulation_residual_rel=parse(Float64, sr[:linear_residual_rel]),
                    circulation_consistency_rel=parse(Float64,
                        sr[:circulation_consistency_rel]),
                    fixed_point_iterations=parse(Int, sr[:fixed_point_iterations])))
                continue
            end
        end
        depth_tag = replace(@sprintf("%.1f", depth), "." => "p")
        out = _run_task5_route(path, "march", "deformed", route, t3, p;
            path_depth=depth, artifact_variant="deformed_depth$(depth_tag)",
            blas_threads)
        r = out.result
        push!(rows, (; alpha_deg, path_depth=depth, route,
            c_error_rel=r.c_error_rel, c_error_l2=r.c_error_l2,
            c_error_max=r.c_error_max, c_l2=r.c_l2,
            c_oracle_l2=r.c_oracle_l2, CL=r.CL,
            f_oracle=out.closures.f_oracle,
            formulation_residual_rel=r.residual_rel,
            circulation_consistency_rel=r.gamma_consistency_rel,
            fixed_point_iterations=route == "single" ? 0 : length(out.iterations.rows)))
    end
    _write_rows_csv(joinpath(path, "task5_line_integral_path_depth.csv"), rows)
    _write_task5_manifests(path; alpha_deg, blas_threads)
    return rows
end

function main(args=ARGS)
    options = _cli_options(args)
    task, alpha_deg, path = options.task, options.alpha_deg, options.output_dir
    task in ("1", "task1", "baseline", "semiinfinite") &&
        return run_task1(; path, alpha_deg)
    task in ("task1-grid", "task1_grid", "grid") &&
        return run_task1_grid_convergence(; path)
    task in ("task2-flat", "task2_flat", "flat") &&
        return run_task2_flat(; path, alpha_deg)
    task in ("task2-flat-das005", "task2_flat_das005", "flat-das005",
        "flat_das005") && return run_task2_flat_das005(; path, alpha_deg)
    task in ("task2-march", "task2_march", "march") &&
        return run_task2_march(; path, alpha_deg)
    task in ("2", "task2") && return run_task2(; path, alpha_deg)
    task in ("task3-flat", "task3_flat") &&
        return run_task3_flat(; path, alpha_deg)
    task in ("task3-flat-das005", "task3_flat_das005") &&
        return run_task3_flat_das005(; path, alpha_deg)
    task in ("task3-march", "task3_march") &&
        return run_task3_march(; path, alpha_deg)
    task in ("3", "task3") && return run_task3(; path, alpha_deg)
    task in ("task4-flat", "task4_flat") &&
        return run_task4_flat(; path, alpha_deg)
    task in ("task4-flat-das005", "task4_flat_das005") &&
        return run_task4_flat_das005(; path, alpha_deg)
    task in ("task4-march", "task4_march") &&
        return run_task4_march(; path, alpha_deg)
    task in ("4", "task4") && return run_task4(; path, alpha_deg)
    task in ("5", "task5") && return run_task5(; path, alpha_deg)
    task in ("task5-depth-march", "task5_depth_march") &&
        return run_task5_depth_march(; path, alpha_deg)
    task in ("task5-finalize", "task5_finalize") &&
        return finalize_task5(path; alpha_deg)
    task in ("task5-audit", "task5_audit") && return run_task5_audit(; path)
    startswith(replace(task, '_' => '-'), "task5-") &&
        return run_task5_selector(task; path, alpha_deg)
    task in ("green-identity", "green_identity", "green") &&
        return run_green_identity_residual(; path, alpha_deg)
    task in ("green-refinement", "green_refinement") &&
        return run_green_identity_refinement(; path)
    task in ("green-conditioning", "green_conditioning") &&
        return run_green_identity_conditioning(; path)
    task in ("green-regularization", "green_regularization") &&
        return run_green_identity_regularization(; path, alpha_deg)
    task in ("green-regularization-refinement", "green_regularization_refinement") &&
        return run_green_identity_regularization_refinement(; path)
    throw(ArgumentError("unknown task selector '$task'"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
