## Procedural capped Dirichlet wing: PressureBernoulli versus PressureLaplace.
##
## The geometry follows pitching_wing.jl: a closed NACA 00xx extrusion with
## rounded, watertight tips.  Angle of attack is carried by the freestream so
## that every refinement level represents exactly the same surface.

import FLOWPanel as pnl
using LinearAlgebra: cross, norm
using Printf: @printf

function pressure_comparison_naca00xx(n::Integer=41; thickness::Real=0.15)
    n >= 21 || throw(ArgumentError("n_airfoil must be at least 21; got $(n)"))
    0 < thickness < 1 || throw(ArgumentError(
        "thickness must be a chord fraction in (0, 1); got $(thickness)"))

    n_half = max(11, cld(n, 2))
    beta = range(0.0, pi, length=n_half)
    x = 0.5 .* (1 .- cos.(beta))
    yt = 5 * thickness .* (
        0.2969 .* sqrt.(x) .- 0.1260 .* x .- 0.3516 .* x.^2 .+
        0.2843 .* x.^3 .- 0.1036 .* x.^4
    )
    lower = hcat(reverse(x), -reverse(yt))
    upper = hcat(x[2:end-1], yt[2:end-1])
    return vcat(lower, upper)
end

function pressure_comparison_wing_mesh(chord, span;
        thickness::Real=0.15, n_span::Integer=4, n_airfoil::Integer=41,
        n_endcap::Integer=5)
    chord > 0 || throw(ArgumentError("chord must be positive; got $(chord)"))
    span > 0 || throw(ArgumentError("span must be positive; got $(span)"))
    n_span >= 1 || throw(ArgumentError("n_span must be positive; got $(n_span)"))
    n_endcap >= 3 || throw(ArgumentError("n_endcap must be at least 3; got $(n_endcap)"))

    contour = pressure_comparison_naca00xx(n_airfoil; thickness)
    n_sec = size(contour, 1)
    n_chord = cld(n_sec, 2) + 1
    y = collect(range(-span / 2, stop=span / 2, length=n_span + 1))
    node_index(i, j) = i + (j - 1) * n_sec
    lower_index(k) = n_chord - k + 1
    upper_index(k) = k == 1 || k == n_chord ? lower_index(k) : n_chord + k - 1

    n_main_nodes = n_sec * length(y)
    n_cap_nodes = 2 * (n_chord - 2) * (n_endcap - 2)
    nodes = zeros(Float64, 3, n_main_nodes + n_cap_nodes)
    for (j, yj) in enumerate(y), i in 1:n_sec
        nodes[:, node_index(i, j)] .= (chord * contour[i, 1], yj,
                                       chord * contour[i, 2])
    end

    n_main_cells = 2 * n_sec * n_span
    n_cap_cells = 4 * (n_chord - 2) * (n_endcap - 1)
    cells = zeros(Int, 3, n_main_cells + n_cap_cells)
    k = 0
    for j in 1:n_span, i in 1:n_sec
        ip = i == n_sec ? 1 : i + 1
        n11, n21 = node_index(i, j), node_index(ip, j)
        n12, n22 = node_index(i, j + 1), node_index(ip, j + 1)
        if i < n_chord
            k += 1; cells[:, k] .= (n11, n21, n12)
            k += 1; cells[:, k] .= (n21, n22, n12)
        else
            k += 1; cells[:, k] .= (n11, n21, n22)
            k += 1; cells[:, k] .= (n11, n22, n12)
        end
    end

    next_node = n_main_nodes
    theta = range(-pi / 2, stop=pi / 2, length=n_endcap)
    for (jtip, side) in ((1, -1.0), (length(y), 1.0))
        rings = Vector{Vector{Int}}(undef, n_chord)
        rings[1] = [node_index(lower_index(1), jtip)]
        rings[end] = [node_index(lower_index(n_chord), jtip)]
        for ix in 2:n_chord-1
            lower = node_index(lower_index(ix), jtip)
            upper = node_index(upper_index(ix), jtip)
            x = nodes[1, lower]
            zmid = (nodes[3, lower] + nodes[3, upper]) / 2
            radius = (nodes[3, upper] - nodes[3, lower]) / 2
            ring = Vector{Int}(undef, n_endcap)
            ring[1], ring[end] = lower, upper
            for itheta in 2:n_endcap-1
                next_node += 1
                ring[itheta] = next_node
                nodes[:, next_node] .= (x,
                    side * (span / 2 + radius * cos(theta[itheta])),
                    zmid + radius * sin(theta[itheta]))
            end
            rings[ix] = ring
        end

        add_cap_cell!(tri) = begin
            k += 1
            cells[:, k] .= side > 0 ? tri : reverse(tri)
        end
        for q in 1:n_endcap-1
            add_cap_cell!((rings[1][1], rings[2][q], rings[2][q + 1]))
        end
        for ix in 2:n_chord-2, q in 1:n_endcap-1
            a, d = rings[ix], rings[ix + 1]
            if theta[q + 1] <= 0
                add_cap_cell!((a[q], d[q], a[q + 1]))
                add_cap_cell!((d[q], d[q + 1], a[q + 1]))
            else
                add_cap_cell!((a[q], d[q], d[q + 1]))
                add_cap_cell!((a[q], d[q + 1], a[q + 1]))
            end
        end
        for q in 1:n_endcap-1
            add_cap_cell!((rings[end - 1][q], rings[end][1], rings[end - 1][q + 1]))
        end
    end

    next_node == size(nodes, 2) || error("internal end-cap node count mismatch")
    k == size(cells, 2) || error("internal end-cap cell count mismatch")
    pnl.iswatertight(nodes, cells)[1] || error("procedural wing mesh is not watertight")
    area2_min = minimum(norm(cross(
        nodes[:, cells[2, i]] - nodes[:, cells[1, i]],
        nodes[:, cells[3, i]] - nodes[:, cells[1, i]])) for i in axes(cells, 2))
    area2_min > 100 * eps(Float64) * chord^2 || error(
        "procedural wing mesh contains a degenerate triangle")
    return nodes, cells
end

function pressure_comparison_shedding(nodes, cells, chord)
    x_te = float(chord)
    tol = max(100 * eps(Float64) * max(abs(x_te), 1.0), 1e-8 * abs(x_te))
    te_nodes = findall(i -> isapprox(nodes[1, i], x_te; atol=tol, rtol=0),
                       axes(nodes, 2))
    length(te_nodes) >= 2 || error("procedural wing has fewer than two trailing-edge nodes")
    sort!(te_nodes, by=i -> nodes[2, i])
    lower = [x_te - tol, minimum(nodes[2, te_nodes]) - tol,
             minimum(nodes[3, te_nodes]) - tol]
    upper = [x_te + tol, maximum(nodes[2, te_nodes]) + tol,
             maximum(nodes[3, te_nodes]) + tol]
    return pnl.calc_shedding_from_seed(nodes, cells, te_nodes[1], te_nodes[2];
        end_node=te_nodes[end], bbox=(lower, upper), normal_jump_tol=1.0,
        max_turn_angle=pi/2)
end

function set_pressure_comparison_wake!(body, direction)
    dhat = collect(direction) ./ norm(direction)
    for Das in body.Das, j in axes(Das, 2)
        Das[:, j] .= dhat
    end
    return body
end

function _resolve_pressure_comparison_span(chord::Real, aspect_ratio::Real,
                                           span::Union{Nothing,Real})
    chord > 0 || throw(ArgumentError("chord must be positive; got $(chord)"))
    aspect_ratio > 0 || throw(ArgumentError(
        "aspect_ratio must be positive; got $(aspect_ratio)"))
    resolved_span = isnothing(span) ? chord * aspect_ratio : float(span)
    resolved_span > 0 || throw(ArgumentError("span must be positive; got $(span)"))
    return resolved_span, resolved_span / chord
end

function build_pressure_comparison_wing(; chord::Real=1.0,
        aspect_ratio::Real=4.0, span::Union{Nothing,Real}=nothing,
        thickness::Real=0.15, n_span::Integer=4, n_airfoil::Integer=41,
        n_endcap::Integer=5, kernel_offset::Real=1e-6)
    bodytype = pnl.RigidWakeBody{
        Union{pnl.ConstantSource, pnl.ConstantDoublet}, 2, Float64, true}
    options = (; kerneloffset=kernel_offset * chord, kernelcutoff=1e-12 * chord,
        semiinfinite_wake=true, watertight=true)
    resolved_span, _ = _resolve_pressure_comparison_span(chord, aspect_ratio, span)
    nodes, cells = pressure_comparison_wing_mesh(chord, resolved_span;
        thickness, n_span, n_airfoil, n_endcap)

    # The first construction establishes final winding.  Shedding must be
    # detected from those cells because it stores cell-local edge indices.
    base = bodytype(nodes, cells, zeros(Int, 6, 0); options...)
    shedding = pressure_comparison_shedding(base.nodes, base.cells, chord)
    return bodytype(copy(base.nodes), copy(base.cells), [shedding]; options...)
end

function solve_pressure_comparison_wing!(body, Vinf;
        rho::Real=1.225, chord::Real=1.0, span::Real=4.0,
        backend=pnl.DirectBackend(), gradient_mode::Symbol=:edge_difference,
        acceleration_form::Symbol=:material_derivative)
    set_pressure_comparison_wake!(body, Vinf)
    normalization = pnl.WingNormalization(rho, chord * span, chord)
    bernoulli = pnl.PressureBernoulli(rho; unsteady=false,
        correct_kuttacondition=false, backend)
    force_bernoulli = pnl.ForceMonitor(1, 1; i_frame=-1, normalization,
        correct_kuttacondition=false, verbose=false)
    laplace = pnl.PressureLaplace((body,), rho; reference_panel=1,
        reference_pressure=0.0, gradient_mode, acceleration_form,
        atol=0.0, rtol=1e-14, itmax=10_000, verbose=false)
    force_laplace = pnl.ForceMonitor(1, 1; i_frame=-1, normalization,
        correct_kuttacondition=false, verbose=false)
    monitors = (bernoulli, force_bernoulli, laplace, force_laplace)

    solver = pnl.Backslash(body)
    # Supply the real monitor stack to steady! so its dependency audit enables
    # Hessians and exterior induced-vorticity accumulation before the solve.
    # This exercises the same production path used by simulations.
    runtime = @elapsed pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
        body_solvers=solver, backend, monitors, path=nothing, verbose=false)
    normal_leakage = vec(sum(body.velocity .* body.normals; dims=1))
    Gcheck = zeros(Float64, body.ncells, body.ncells)
    pnl._G!(Gcheck, body, body; kerneloffset=body.kerneloffset_panel,
        update_geometry=false)
    linear_residual = Gcheck * view(body.strength, :, 2) - solver.rhs
    solver_residual = norm(linear_residual) / max(norm(solver.rhs), eps())
    return (; bernoulli, laplace, force_bernoulli, force_laplace,
        normal_leakage, solver_residual, runtime)
end

function compare_pressure_models(; chord::Real=1.0, aspect_ratio::Real=4.0,
        span::Union{Nothing,Real}=nothing,
        thickness::Real=0.15, n_span::Integer=4, n_airfoil::Integer=41,
        n_endcap::Integer=5, alpha_deg::Real=5.88, speed::Real=56.0,
        rho::Real=1.225, backend=pnl.DirectBackend(),
        gradient_mode::Symbol=:edge_difference,
        acceleration_form::Symbol=:material_derivative)
    # Other gradient modes remain available for explicit diagnostics. Only the
    # edge-difference/material-derivative defaults are used for acceptance.
    resolved_span, nominal_aspect_ratio =
        _resolve_pressure_comparison_span(chord, aspect_ratio, span)
    body = build_pressure_comparison_wing(; chord, aspect_ratio,
        span=resolved_span, thickness,
        n_span, n_airfoil, n_endcap)
    Vinf = speed .* [cosd(alpha_deg), 0.0, sind(alpha_deg)]
    monitors = solve_pressure_comparison_wing!(body, Vinf;
        rho, chord, span=resolved_span, backend, gradient_mode, acceleration_form)

    p_b_raw = monitors.bernoulli.pressure[1]
    p_l_raw = monitors.laplace.p[1]
    areas = pnl.calc_areas(body)
    gauge_shift = sum(areas .* (p_b_raw .- p_l_raw)) / sum(areas)
    p_l = p_l_raw .+ gauge_shift
    pressure_error = p_l .- p_b_raw
    pressure_l2_error = sqrt(sum(areas .* pressure_error.^2) /
        max(sum(areas .* p_b_raw.^2), eps()))
    pressure_linf_error = maximum(abs, pressure_error) /
        max(maximum(abs, p_b_raw), eps())
    cf_b = vec(monitors.force_bernoulli.force[:, 1])
    cf_l = vec(monitors.force_laplace.force[:, 1])
    force_rel_error = norm(cf_l .- cf_b) / max(norm(cf_b), eps())

    return (; body, Vinf, pressure_bernoulli=copy(p_b_raw),
        pressure_laplace=copy(p_l), gauge_shift,
        pressure_rel_error=pressure_l2_error, pressure_l2_error,
        pressure_linf_error,
        force_rel_error, bernoulli_force=cf_b, laplace_force=cf_l,
        npanels=body.ncells, n_span, n_airfoil, n_endcap, gradient_mode,
        acceleration_form, chord=float(chord), span=resolved_span,
        aspect_ratio=nominal_aspect_ratio, reference_area=chord * resolved_span,
        raw_normal_leakage=monitors.normal_leakage,
        solver_residual=monitors.solver_residual, runtime=monitors.runtime)
end

const PRESSURE_COMPARISON_LEVELS = (
    (name="coarse", n_span=2, n_airfoil=21, n_endcap=3),
    (name="medium", n_span=4, n_airfoil=41, n_endcap=5),
    (name="fine",   n_span=8, n_airfoil=81, n_endcap=7),
)

function run_pressure_refinement(; levels=PRESSURE_COMPARISON_LEVELS, kwargs...)
    results = NamedTuple[]
    for level in levels
        result = compare_pressure_models(; n_span=level.n_span,
            n_airfoil=level.n_airfoil, n_endcap=level.n_endcap, kwargs...)
        push!(results, merge((; name=level.name), result))
    end
    return results
end

function print_pressure_result(result)
    @printf("%-7s panels=%5d  AR=%5.2f  L2(p)=%10.3e  Linf(p)=%10.3e  rel(F)=%10.3e  bc=%9.2e  t=%7.3fs\n",
        result.name, result.npanels, result.aspect_ratio,
        result.pressure_l2_error, result.pressure_linf_error,
        result.force_rel_error, result.solver_residual, result.runtime)
    @printf("         C_B = [% .8e, % .8e, % .8e]\n", result.bernoulli_force...)
    @printf("         C_L = [% .8e, % .8e, % .8e]\n", result.laplace_force...)
end

function main()
    println("\n#===== PROCEDURAL CAPPED-WING PRESSURE COMPARISON =====#")
    println("PressureLaplace: edge_difference + material_derivative")
    results = run_pressure_refinement()
    foreach(print_pressure_result, results)
    all(diff([r.pressure_l2_error for r in results]) .<= 1e-8) ||
        error("pressure error increased under refinement")
    all(diff([r.force_rel_error for r in results]) .<= 1e-8) ||
        error("force-vector error increased under refinement")
    fine = last(results)
    fine.pressure_rel_error < 3e-2 || error("fine pressure error exceeds 3e-2")
    fine.force_rel_error < 2e-3 || error("fine force error exceeds 2e-3")
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
