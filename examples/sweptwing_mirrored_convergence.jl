#=##############################################################################
# DESCRIPTION
#   Convergence sweep for the symmetric mirrored swept-wing discretization.
#
#   This intentionally does one steady solve per row: solve strengths, compute
#   Bernoulli pressure, integrate force. No Laplace solve, plots, VTK, or rotor
#   cases are run.
#
#   Environment knobs:
#     FLOWPANEL_SWEEP_CASES="8:40,10:48,12:64,14:72"
#       Case format is n_rfl:n_span. For chordwise-only convergence at fixed
#       spanwise sections, use values like "4:40,8:40,12:40,16:40".
#     FLOWPANEL_SWEEP_KERNEL_OFFSETS="1e-8,1e-10,1e-12"
#     FLOWPANEL_SWEEP_KERNEL_CASE="8:40"
#     FLOWPANEL_SWEEP_SOLVER="backslash" or "krylov"
#     FLOWPANEL_SWEEP_SOURCE_HALVES="positive,negative"
#     FLOWPANEL_SWEEP_MAX_PANELS=45000
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))

import LinearAlgebra
import CSV
import DataFrames: DataFrame
using Printf

const LA = LinearAlgebra

function parse_case_list(s::AbstractString)
    out = Tuple{Int, Int}[]
    for item in split(s, ",")
        isempty(strip(item)) && continue
        parts = split(strip(item), ":")
        length(parts) == 2 || error("Expected n_rfl:n_span case, got $(item)")
        push!(out, (parse(Int, parts[1]), parse(Int, parts[2])))
    end
    isempty(out) && error("No cases parsed from $(s)")
    return out
end

function parse_float_list(s::AbstractString)
    out = Float64[]
    for item in split(s, ",")
        isempty(strip(item)) && continue
        push!(out, parse(Float64, strip(item)))
    end
    isempty(out) && error("No kernel offsets parsed from $(s)")
    return out
end

function parse_symbol_list(s::AbstractString)
    out = Symbol[]
    for item in split(s, ",")
        isempty(strip(item)) && continue
        push!(out, Symbol(lowercase(strip(item))))
    end
    isempty(out) && error("No source halves parsed from $(s)")
    return out
end

function case_string(cases)
    return join(("$(n_rfl):$(n_span)" for (n_rfl, n_span) in cases), ",")
end

panel_count(n_rfl::Int, n_span::Int) = 24 * n_rfl * n_span

function filter_cases_by_panel_limit(cases, max_panels::Int)
    kept = Tuple{Int, Int}[]
    skipped = Tuple{Int, Int, Int}[]
    for (n_rfl, n_span) in cases
        panels = panel_count(n_rfl, n_span)
        if panels <= max_panels
            push!(kept, (n_rfl, n_span))
        else
            push!(skipped, (n_rfl, n_span, panels))
        end
    end
    isempty(kept) && error("All requested cases exceed max panel limit $(max_panels).")
    return kept, skipped
end

function sweptwing_constants()
    AOA = 4.2
    magVinf = 30.0
    Vinf = magVinf * [cosd(AOA), 0.0, sind(AOA)]
    rho = 1.225
    b = 98 * 0.0254
    ar = 5.0
    return (; AOA, magVinf, Vinf, rho, b, ar,
        tr=1.0, twist_root=0.0, twist_tip=0.0, lambda=45.0, gamma=0.0,
        airfoil="airfoil-rae101.csv",
        airfoil_path=joinpath(pnl.examples_path, "data"),
        Sref=b^2 / ar,
        c_ref=b / ar)
end

function discretization(n_rfl::Int, n_span::Int)
    rfl = [(0.25, n_rfl, 10.0, false),
           (0.50, n_rfl, 1.0, true),
           (0.25, n_rfl, 0.1, false)]
    span = [(1.0, n_span, 1.0, true)]
    return rfl, span
end

function mirrored_negative_half(c, n_rfl::Int, n_span::Int, core_size::Real)
    rfl, span = discretization(n_rfl, n_span)
    bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
    half = simplewing(c.b, c.ar, c.tr, c.twist_root, c.twist_tip, c.lambda, c.gamma;
        bodytype=bodytype,
        bodyoptargs=(; core_size),
        airfoil_root=c.airfoil,
        airfoil_tip=c.airfoil,
        airfoil_path=c.airfoil_path,
        rfl_NDIVS=rfl,
        span_NDIVS=span,
        delim=",",
        b_low=-1.0,
        b_up=0.0,
        verify_spline=false,
        verify_rflspline=false)

    half_nodes = half.nodes
    half_cells = half.cells
    mirror_index = Vector{Int}(undef, size(half_nodes, 2))
    nodes = copy(half_nodes)
    for ni in axes(half_nodes, 2)
        if abs(half_nodes[2, ni]) <= 100eps(Float64)
            mirror_index[ni] = ni
        else
            nodes = hcat(nodes, [half_nodes[1, ni], -half_nodes[2, ni], half_nodes[3, ni]])
            mirror_index[ni] = size(nodes, 2)
        end
    end

    half_centers_y = [sum(half_nodes[2, half_cells[:, ci]]) / 3 for ci in axes(half_cells, 2)]
    neg_order = sort(collect(axes(half_cells, 2)); by=ci -> half_centers_y[ci])
    pos_order = sort(collect(axes(half_cells, 2)); by=ci -> -half_centers_y[ci])
    cells = Matrix{Int}(undef, 3, 2 * size(half_cells, 2))
    out_ci = 0
    for ci in neg_order
        out_ci += 1
        cells[:, out_ci] .= half_cells[:, ci]
    end
    for ci in pos_order
        out_ci += 1
        cells[:, out_ci] .= reverse(mirror_index[half_cells[:, ci]])
    end

    te_nodes = Int[]
    for col in eachcol(half.shedding[1])
        pi, nia, nib = col[1], col[2], col[3]
        push!(te_nodes, half_cells[nia, pi])
        push!(te_nodes, half_cells[nib, pi])
    end
    full_te_nodes = unique(vcat(te_nodes, mirror_index[te_nodes]))
    sort!(full_te_nodes; by=ni -> nodes[2, ni])
    shedding = pnl.calc_shedding(nodes, cells, full_te_nodes, zeros(eltype(nodes), 3, 0))
    watertight, _ = pnl.iswatertight(nodes, cells)
    return bodytype(nodes, cells, [shedding]; watertight, ensure_winding=false,
        core_size)
end

function build_body(c, n_rfl::Int, n_span::Int; core_size::Real=1e-10,
                    source_half::Symbol=:positive)
    if source_half == :positive
        rfl, span = discretization(n_rfl, n_span)
        bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
        return simplewing_mirrored(c.b, c.ar, c.tr, c.twist_root, c.twist_tip,
            c.lambda, c.gamma;
            bodytype=bodytype,
            bodyoptargs=(; core_size),
            airfoil_root=c.airfoil,
            airfoil_tip=c.airfoil,
            airfoil_path=c.airfoil_path,
            rfl_NDIVS=rfl,
            span_NDIVS=span,
            delim=",",
            verify_spline=false,
            verify_rflspline=false)
    elseif source_half == :negative
        return mirrored_negative_half(c, n_rfl, n_span, core_size)
    else
        error("source_half must be :positive or :negative; got $(source_half)")
    end
end

function apply_wake_direction!(body, Vinf)
    wake_direction = reshape(Vinf ./ LA.norm(Vinf), :, 1)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
    end
    return body
end

function solver_for(body, solver_name::Symbol)
    if solver_name == :backslash
        return pnl.Backslash(body)
    elseif solver_name == :krylov
        return pnl.KrylovSolver(body; backend=pnl.FastMultipoleBackend(),
            method=:gmres, atol=1e-9, rtol=1e-9, itmax=200)
    else
        error("Unknown solver $(solver_name); use :backslash or :krylov")
    end
end

function solve_case(n_rfl::Int, n_span::Int; core_size::Real=1e-10,
                    solver_name::Symbol=:backslash, source_half::Symbol=:positive)
    c = sweptwing_constants()
    body = build_body(c, n_rfl, n_span; core_size, source_half)
    apply_wake_direction!(body, c.Vinf)

    Dhat = c.Vinf / LA.norm(c.Vinf)
    Lhat = LA.cross(Dhat, [0.0, 1.0, 0.0])
    pressure = pnl.PressureBernoulli(c.rho; correct_kuttacondition=false)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(c.rho, c.Sref, c.c_ref),
        correct_kuttacondition=false)

    elapsed = @elapsed pnl.steady!(body, pnl.ReferenceFrame(body), c.Vinf;
        body_solvers=solver_for(body, solver_name),
        backend=pnl.FastMultipoleBackend(),
        monitors=(pressure, force),
        path=nothing,
        verbose=false)

    y = view(body.controlpoints, 2, :)
    lift = [LA.dot(view(force.distributed_force, :, i), Lhat) for i in 1:body.ncells]
    areas = pnl.calc_areas(body)
    hmin = sqrt(minimum(areas))
    CL = LA.dot(force.force[:, 1], Lhat)
    CD = LA.dot(force.force[:, 1], Dhat)
    return (; n_rfl, n_span, panels=body.ncells, shedding=size(body.shedding[1], 2),
        core_size, offset_over_hmin=core_size / hmin, source_half,
        CL, CD, CL_error_pct=100 * (CL - 0.238) / 0.238,
        lift_positive=sum(lift[y .> 0]), lift_negative=sum(lift[y .< 0]),
        elapsed)
end

function print_header()
    @printf("%5s %6s %8s %6s %11s %11s %10s %11s %11s %11s %11s %8s\n",
        "nRFL", "nSpan", "panels", "shed", "kernel", "k/hmin",
        "CL", "CLerr%", "CD", "L y>0", "L y<0", "time_s")
end

function print_result(r)
    @printf("%5d %6d %8d %6d %11.1e %11.3e %+10.6f %+11.3f %+11.6f %+11.4f %+11.4f %8.2f\n",
        r.n_rfl, r.n_span, r.panels, r.shedding, r.core_size,
        r.offset_over_hmin, r.CL, r.CL_error_pct, r.CD, r.lift_positive,
        r.lift_negative, r.elapsed)
end

function main()
    cases = parse_case_list(get(ENV, "FLOWPANEL_SWEEP_CASES",
        "4:16,8:40,10:48,12:64,14:72"))
    max_panels = parse(Int, get(ENV, "FLOWPANEL_SWEEP_MAX_PANELS", "45000"))
    cases, skipped = filter_cases_by_panel_limit(cases, max_panels)
    offsets = parse_float_list(get(ENV, "FLOWPANEL_SWEEP_KERNEL_OFFSETS",
        "1e-8,1e-10,1e-12"))
    kernel_case = only(parse_case_list(get(ENV, "FLOWPANEL_SWEEP_KERNEL_CASE",
        "8:40")))
    solver_name = Symbol(lowercase(get(ENV, "FLOWPANEL_SWEEP_SOLVER", "backslash")))
    source_halves = parse_symbol_list(get(ENV, "FLOWPANEL_SWEEP_SOURCE_HALVES",
        get(ENV, "FLOWPANEL_SWEEP_SOURCE_HALF", "positive,negative")))
    source_half_string = join(string.(source_halves), ",")
    main_core_size = parse(Float64, get(ENV, "FLOWPANEL_SWEEP_KERNEL_OFFSET", "1e-10"))

    println("# Swept-wing mirrored convergence")
    println("# cases=$(case_string(cases))")
    if !isempty(skipped)
        println("# skipped cases above max_panels=$(max_panels): " *
            join(("$(nr):$(ns) ($(p) panels)" for (nr, ns, p) in skipped), ", "))
    end
    println("# solver=$(solver_name), source_halves=$(source_half_string), main core_size=$(main_core_size)")
    println("# max_panels=$(max_panels)")
    println("# CLexp=0.238")
    csv_path = get(ENV, "FLOWPANEL_SWEEP_CSV", "")
    csv_rows = NamedTuple[]
    for source_half in source_halves
        println("\n# Source half: $(source_half)")
        print_header()
        for (n_rfl, n_span) in cases
            r = solve_case(n_rfl, n_span; core_size=main_core_size,
                solver_name, source_half)
            print_result(r)
            push!(csv_rows, r)
            flush(stdout)
        end
    end
    if !isempty(csv_path)
        mkpath(dirname(csv_path))
        CSV.write(csv_path, DataFrame(csv_rows))
        println("\n# Wrote $(length(csv_rows)) rows to $(csv_path)")
    end

    if panel_count(kernel_case...) > max_panels
        println("\n# Kernel-offset sensitivity skipped: $(kernel_case[1]):$(kernel_case[2]) exceeds max_panels=$(max_panels)")
    else
        for source_half in source_halves
            println("\n# Kernel-offset sensitivity at $(kernel_case[1]):$(kernel_case[2]), source_half=$(source_half)")
            print_header()
            for core_size in offsets
                print_result(solve_case(kernel_case[1], kernel_case[2]; core_size,
                    solver_name, source_half))
                flush(stdout)
            end
        end
    end
end

main()
