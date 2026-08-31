#=##############################################################################
# DESCRIPTION
#   Rectangular NACA 0015 wing pitching sinusoidally about the quarter chord.
#
#   The script first runs a static lift polar for alpha in 0:10 deg, then runs
#   the unsteady pitching case. Geometry and flow inputs follow the requested
#   English-unit case, but all solver inputs are converted to SI.
#
# AUTHORSHIP
#   * Created   : Jul 2026
#   * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import FLOWPanel: _run_monitor!, monitor_provides, monitor_requires
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: cross, dot, norm
import Printf: @printf, @sprintf
import Statistics: mean
import TOML
using LaTeXStrings

const DEFAULT_PITCHING_WING_PATH = joinpath("data", "pitching_wing")
const PITCHING_WING_NAME = "pitching_wing"
const PITCHING_WING_CONFIG_FILE = PITCHING_WING_NAME * ".config.toml"
const PITCHING_WING_PRESSURE_STATE_FILE = PITCHING_WING_NAME * ".pressure_state.csv"
const FT_TO_M = 0.3048

function _resolve_pitching_wing_dimensions(c_ft::Real, aspect_ratio::Real,
        span_ft::Union{Nothing,Real}, semispan_ft::Union{Nothing,Real})
    c_ft > 0 || throw(ArgumentError("c_ft must be positive; got $(c_ft)"))
    aspect_ratio > 0 || throw(ArgumentError(
        "aspect_ratio must be positive; got $(aspect_ratio)"))
    !isnothing(span_ft) && !isnothing(semispan_ft) && throw(ArgumentError(
        "Specify at most one of span_ft and the legacy semispan_ft keyword."))
    resolved_span_ft = if !isnothing(span_ft)
        float(span_ft)
    elseif !isnothing(semispan_ft)
        2 * float(semispan_ft)
    else
        float(aspect_ratio) * float(c_ft)
    end
    resolved_span_ft > 0 || throw(ArgumentError(
        "resolved span must be positive; got $(resolved_span_ft) ft"))
    return (; c=c_ft * FT_TO_M, b=resolved_span_ft * FT_TO_M,
        span_ft=resolved_span_ft, semispan_ft=resolved_span_ft / 2,
        aspect_ratio=resolved_span_ft / c_ft)
end

function naca00xx_contour(n::Integer=121; thickness::Real=0.15)
    n >= 21 || throw(ArgumentError("naca00xx_contour requires at least 21 points"))
    0 < thickness < 1 || throw(ArgumentError("thickness must be a chord fraction in (0, 1); got $(thickness)"))

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

function pitching_wing_mesh(c, b; n_span::Integer=13, n_airfoil::Integer=161,
                            n_endcap::Integer=9, endcap::Symbol=:round,
                            caps::Bool=true, thickness::Real=0.15)
    endcap == :round || throw(ArgumentError(
        "pitching_wing_mesh currently implements endcap=:round only; got $(endcap)."))
    n_span >= 1 || throw(ArgumentError("n_span must be positive; got $(n_span)"))
    n_endcap >= 3 || throw(ArgumentError(
        "n_endcap requires at least 3 semicircle points; got $(n_endcap)"))
    contour = naca00xx_contour(n_airfoil; thickness=thickness)
    n_sec = size(contour, 1)
    n_chord = cld(n_sec, 2) + 1
    y = collect(range(-b / 2, stop=b / 2, length=n_span + 1))

    node_index(i, j) = i + (j - 1) * n_sec
    lower_index(k) = n_chord - k + 1
    upper_index(k) = k == 1 || k == n_chord ? lower_index(k) : n_chord + k - 1

    n_main_nodes = n_sec * length(y)
    n_cap_nodes = caps ? 2 * (n_chord - 2) * (n_endcap - 2) : 0
    nodes = zeros(Float64, 3, n_main_nodes + n_cap_nodes)
    for (j, yj) in enumerate(y), i in 1:n_sec
        x = c * contour[i, 1]
        z = c * contour[i, 2]
        nodes[:, node_index(i, j)] .= (x, yj, z)
    end

    n_main_cells = 2 * n_sec * n_span
    n_cap_cells = caps ? 4 * (n_chord - 2) * (n_endcap - 1) : 0
    cells = zeros(Int, 3, n_main_cells + n_cap_cells)
    k = 0
    for j in 1:n_span
        for i in 1:n_sec
            ip = i == n_sec ? 1 : i + 1
            n11 = node_index(i, j)
            n21 = node_index(ip, j)
            n12 = node_index(i, j + 1)
            n22 = node_index(ip, j + 1)
            if i < n_chord
                # Reflect the upper-surface diagonal across z=0 on every lower-
                # surface quad.  Besides making the complete surface mesh
                # symmetric, this keeps the two trailing-edge off-edge vertices
                # at the same span station.
                k += 1
                cells[:, k] .= (n11, n21, n12)
                k += 1
                cells[:, k] .= (n21, n22, n12)
            else
                k += 1
                cells[:, k] .= (n11, n21, n22)
                k += 1
                cells[:, k] .= (n11, n22, n12)
            end
        end
    end

    next_node = n_main_nodes
    theta = range(-pi / 2, stop=pi / 2, length=n_endcap)
    caps && for (jtip, side) in ((1, -1.0), (length(y), 1.0))
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
            ring[1] = lower
            ring[end] = upper
            for itheta in 2:n_endcap-1
                next_node += 1
                ring[itheta] = next_node
                nodes[:, next_node] .= (
                    x,
                    side * (b / 2 + radius * cos(theta[itheta])),
                    zmid + radius * sin(theta[itheta]),
                )
            end
            rings[ix] = ring
        end

        add_cap_cell!(tri) = begin
            k += 1
            cells[:, k] .= side > 0 ? tri : reverse(tri)
        end

        # Rounded leading-edge corner collapses the first semicircle into the LE point.
        for q in 1:n_endcap-1
            add_cap_cell!((rings[1][1], rings[2][q], rings[2][q + 1]))
        end

        # Semicircular y-z sections cover the regular chordwise portion.  Use
        # reflected diagonals below and above z=0 so each cap triangle has a
        # partner under z -> -z.  The default odd n_endcap places a ring node
        # on z=0; for an even count the straddling interval uses the upper
        # convention because no diagonal of that quad is itself reflection-
        # symmetric without adding a node.
        for ix in 2:n_chord-2, q in 1:n_endcap-1
            a = rings[ix]
            d = rings[ix + 1]
            if theta[q + 1] <= 0
                add_cap_cell!((a[q], d[q], a[q + 1]))
                add_cap_cell!((d[q], d[q + 1], a[q + 1]))
            else
                add_cap_cell!((a[q], d[q], d[q + 1]))
                add_cap_cell!((a[q], d[q + 1], a[q + 1]))
            end
        end

        # The trailing-edge cone collapses the last semicircle into the closed TE point.
        for q in 1:n_endcap-1
            add_cap_cell!((rings[end - 1][q], rings[end][1], rings[end - 1][q + 1]))
        end
    end

    next_node == size(nodes, 2) || error("internal end-cap node count mismatch")
    k == size(cells, 2) || error("internal end-cap cell count mismatch")
    if caps
        watertight, _ = pnl.iswatertight(nodes, cells)
        watertight || error(
            "pitching-wing end-cap construction did not produce a watertight mesh")
    end
    area2_min = minimum(norm(cross(
        nodes[:, cells[2, i]] - nodes[:, cells[1, i]],
        nodes[:, cells[3, i]] - nodes[:, cells[1, i]],
    )) for i in axes(cells, 2))
    area2_min > 100 * eps(Float64) * c^2 || error(
        "pitching-wing mesh contains a degenerate triangle (twice-area=$(area2_min))")

    return nodes, cells
end

function calc_pitching_wing_shedding(nodes, cells, c)
    x_te = float(c)
    te_tol = max(100 * eps(Float64) * max(abs(x_te), 1.0), 1e-8 * abs(x_te))
    te_nodes = findall(i -> isapprox(nodes[1, i], x_te; atol=te_tol, rtol=0), axes(nodes, 2))
    length(te_nodes) >= 2 || throw(ArgumentError(
        "pitching wing mesh requires at least two trailing-edge nodes at x=$(x_te); found $(length(te_nodes))"))

    sort!(te_nodes, by = i -> nodes[2, i])
    lower = [
        x_te - te_tol,
        minimum(nodes[2, te_nodes]) - te_tol,
        minimum(nodes[3, te_nodes]) - te_tol,
    ]
    upper = [
        x_te + te_tol,
        maximum(nodes[2, te_nodes]) + te_tol,
        maximum(nodes[3, te_nodes]) + te_tol,
    ]

    return pnl.calc_shedding_from_seed(nodes, cells, te_nodes[1], te_nodes[2];
        end_node=te_nodes[end],
        bbox=(lower, upper),
        normal_jump_tol=1.0,
        max_turn_angle=pi/2,
    )
end

function set_wake_Das!(body, direction; magnitude::Real=1.0)
    dhat = SVector{3, Float64}(direction[1], direction[2], direction[3])
    dhat = dhat / sqrt(sum(abs2, dhat))
    for Das in body.Das
        for j in axes(Das, 2)
            Das[:, j] .= magnitude .* dhat
        end
    end
    return body
end

function build_pitching_wing_body(c, b; n_span::Integer=13, n_airfoil::Integer=161,
        n_endcap::Integer=9, endcap::Symbol=:round, semiinfinite_wake::Bool=false,
        thickness::Real=0.15)
    bodytype = pnl.RigidWakeBody{
        Union{pnl.ConstantSource, pnl.ConstantDoublet}, 2, Float64, true}
    bodyoptargs = (;
        core_size=1e-6 * c,
        kernelcutoff=1e-12 * c,
        semiinfinite_wake,
        watertight=true,
    )

    nodes, cells = pitching_wing_mesh(c, b; n_span, n_airfoil, n_endcap, endcap,
        caps=true, thickness)
    base = bodytype(nodes, cells, zeros(Int, 6, 0); bodyoptargs...)
    shedding = calc_pitching_wing_shedding(base.nodes, base.cells, c)
    return bodytype(copy(base.nodes), copy(base.cells), [shedding]; bodyoptargs...)
end

function pitching_wing_frame(body, pivot, alpha_init)
    Rinit = pnl.Rodrigues(SVector{3}(0.0, 1.0, 0.0), alpha_init)
    pnl.rotate_translate!(body, pivot, Rinit, SVector{3}(0.0, 0.0, 0.0))
    pnl.rotate_Das!(body, Rinit)

    return pnl.ReferenceFrame(body;
        origin = pivot,
        v = SVector{3}(0.0, 0.0, 0.0),
        ω_axis = SVector{3}(0.0, 1.0, 0.0),
        ω = 0.0,
        R = Rinit,
        name = "pitching_wing",
        child_index = Int[],
        dependent_index = [1],
    )
end

_uinf_from_alpha(U, alpha_deg) =
    SVector{3}(U * cosd(alpha_deg), 0.0, U * sind(alpha_deg))

function _lift_drag_span_directions(Uvec)
    Dhat = SVector{3}(Uvec[1], Uvec[2], Uvec[3])
    Dhat = Dhat / sqrt(sum(abs2, Dhat))
    Shat = SVector{3}(0.0, 1.0, 0.0)
    Lhat = cross(Dhat, Shat)
    Lhat = Lhat / sqrt(sum(abs2, Lhat))
    return Lhat, Dhat, Shat
end

function _interp_section_values(bin_eta, values, section_eta)
    order = sortperm(bin_eta)
    x = collect(bin_eta[order])
    y = collect(values[order])
    out = similar(collect(float.(section_eta)))
    for (k, eta) in pairs(section_eta)
        eta <= first(x) && (out[k] = first(y); continue)
        eta >= last(x) && (out[k] = last(y); continue)
        j = searchsortedlast(x, eta)
        t = (eta - x[j]) / (x[j + 1] - x[j])
        out[k] = (1 - t) * y[j] + t * y[j + 1]
    end
    return out
end

function _section_interpolation_bins(bin_eta, section_eta)
    x = collect(bin_eta[sortperm(bin_eta)])
    return map(section_eta) do eta
        eta <= first(x) && return (first(x), first(x))
        eta >= last(x) && return (last(x), last(x))
        j = searchsortedlast(x, eta)
        return (x[j], x[j + 1])
    end
end

function _print_section_interpolation_bins(io::IO, bin_eta, section_eta)
    for (eta, (eta_lower, eta_upper)) in
            zip(section_eta, _section_interpolation_bins(bin_eta, section_eta))
        @printf(io, "  requested η=%.4f, interpolation bins η=[%.4f, %.4f]\n",
            eta, eta_lower, eta_upper)
    end
    return nothing
end

function _section_labels(section_eta)
    return ["cl_eta_" * replace(@sprintf("%.4f", eta), "." => "p") for eta in section_eta]
end

function _write_static_polar_csv(path, section_eta, rows)
    mkpath(path)
    filename = joinpath(path, "pitching_wing_static_polar.csv")
    labels = _section_labels(section_eta)
    open(filename, "w") do io
        println(io, join(vcat(["alpha_deg", "CL"], labels), ","))
        for row in rows
            vals = Any[row.alpha_deg, row.CL]
            append!(vals, row.section_cl)
            println(io, join(vals, ","))
        end
    end
    return filename
end

function plot_pitching_wing_static_polar(rows, section_eta; path=DEFAULT_PITCHING_WING_PATH)
    plt = Core.eval(@__MODULE__, :(begin
        import PythonPlot as pitching_wing_plt
        pitching_wing_plt
    end))

    return Base.invokelatest(_plot_pitching_wing_static_polar, plt, rows, collect(section_eta), path)
end

function _plot_pitching_wing_static_polar(plt, rows, section_eta, path)
    mkpath(path)
    alpha = [row.alpha_deg for row in rows]

    section_path = joinpath(path, "pitching_wing_section_cl_polar.png")
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for (i, eta) in pairs(section_eta)
        ax.plot(alpha, [row.section_cl[i] for row in rows], "-o";
            color="C$(i - 1)",
            linewidth=1.5,
            markersize=3,
            label=L"$\eta = %$(round(eta, sigdigits=4))$",
        )
    end
    ax.set_xlabel(L"$\alpha$ (deg)")
    ax.set_ylabel(L"c_\ell")
    ax.grid(true, alpha=0.35)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(section_path, dpi=170)
    plt.pyplot.close(fig)

    CL_path = joinpath(path, "pitching_wing_CL_polar.png")
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(alpha, [row.CL for row in rows], "-o";
        color="C0",
        linewidth=1.5,
        markersize=3,
        label=L"C_L",
    )
    ax.set_xlabel(L"$\alpha$ (deg)")
    ax.set_ylabel(L"C_L")
    ax.grid(true, alpha=0.35)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(CL_path, dpi=170)
    plt.pyplot.close(fig)

    return (; section_cl=section_path, CL=CL_path)
end

function _write_unsteady_csv(path, t_range, period, alpha_history, CL_history, section_eta, section_history)
    mkpath(path)
    filename = joinpath(path, "pitching_wing_unsteady.csv")
    labels = _section_labels(section_eta)
    open(filename, "w") do io
        println(io, join(vcat(["time", "t_over_T", "alpha_deg", "CL"], labels), ","))
        for i in eachindex(t_range)
            vals = Any[t_range[i], t_range[i] / period, alpha_history[i], CL_history[i]]
            append!(vals, section_history[:, i])
            println(io, join(vals, ","))
        end
    end
    return filename
end

_pitching_wing_config_path(path) = joinpath(path, PITCHING_WING_CONFIG_FILE)
_pitching_wing_csv_path(path) = joinpath(path, PITCHING_WING_NAME * "_unsteady.csv")

function _pitching_wing_config(setup)
    return Dict{String, Any}(
        "format_version" => 2,
        "geometry" => Dict(
            "semispan_ft" => float(setup.semispan_ft),
            "span_ft" => float(setup.span_ft),
            "c_ft" => float(setup.c_ft),
            "aspect_ratio" => float(setup.aspect_ratio),
            "n_span" => Int(setup.n_span),
            "n_airfoil" => Int(setup.n_airfoil),
            "n_endcap" => Int(setup.n_endcap),
            "endcap" => String(setup.endcap),
        ),
        "kinematics" => Dict(
            "v_inf_ft_s" => float(setup.v_inf_ft_s),
            "rho" => float(setup.rho),
            "alpha_mean_deg" => float(setup.alpha_mean_deg),
            "alpha_amp_deg" => float(setup.alpha_amp_deg),
            "frequency_hz" => float(setup.frequency_hz),
            "pivot_chord_fraction" => float(setup.pivot_chord_fraction),
        ),
        "time" => Dict(
            "c_per_dt" => float(setup.c_per_dt),
            "dt" => float(setup.dt),
        ),
        "wake" => Dict(
            "model" => String(setup.wake_model),
            "das_chord_fraction" => float(setup.das_chord_fraction),
            "das" => float(setup.das),
            "panel_wake_rows" => Int(setup.panel_wake_rows),
            "max_particles" => Int(setup.max_particles),
            "wake_length_spans" => isnothing(setup.wake_length_spans) ?
                "disabled" : float(setup.wake_length_spans),
            "wake_length_m" => isnothing(setup.wake_length) ?
                "disabled" : float(setup.wake_length),
            "downstream_boundary_from_pivot_m" => isnothing(setup.wake_downstream_boundary) ?
                "disabled" : float(setup.wake_downstream_boundary),
            "method_trailing" => "OverlapPPS(1.3,2)",
            "method_unsteady" => "OverlapPPS(1.3,2)",
        ),
        "sections" => Dict(
            "eta" => collect(float.(setup.section_eta)),
            "n_section_bins" => Int(setup.n_section_bins),
        ),
    )
end

function _clear_pitching_wing_output!(path)
    isdir(path) && rm(path; recursive=true, force=true)
    mkpath(path)
    return path
end

function _write_pitching_wing_config(path, config)
    mkpath(path)
    filename = _pitching_wing_config_path(path)
    open(filename, "w") do io
        TOML.print(io, config; sorted=true)
    end
    return filename
end

function _config_mismatches(saved, requested; prefix="")
    mismatches = String[]
    for key in sort!(collect(keys(requested)))
        location = isempty(prefix) ? key : prefix * "." * key
        if !haskey(saved, key)
            push!(mismatches, "$(location) is missing")
        elseif requested[key] isa AbstractDict
            saved[key] isa AbstractDict || (push!(mismatches, "$(location) has a different type"); continue)
            append!(mismatches, _config_mismatches(saved[key], requested[key]; prefix=location))
        else
            a, b = saved[key], requested[key]
            compatible = if a isa Real && b isa Real
                isapprox(float(a), float(b); atol=100eps(Float64) * max(abs(float(b)), 1.0), rtol=1e-12)
            else
                a == b
            end
            compatible || push!(mismatches, "$(location): saved=$(repr(a)), requested=$(repr(b))")
        end
    end
    return mismatches
end

function _validate_pitching_wing_config(path, requested)
    filename = _pitching_wing_config_path(path)
    isfile(filename) || throw(ArgumentError(
        "Cannot restart legacy pitching-wing output without $(PITCHING_WING_CONFIG_FILE). " *
        "Start fresh with check_existing=false."))
    saved = TOML.parsefile(filename)
    mismatches = _config_mismatches(saved, requested)
    isempty(mismatches) || throw(ArgumentError(
        "Cannot restart pitching-wing output because its configuration is incompatible:\n  " *
        join(mismatches, "\n  ") * "\nStart fresh with check_existing=false or use the saved configuration."))
    return saved
end

function _pitching_wing_existing_state(path)
    files = (
        metadata=joinpath(path, PITCHING_WING_NAME * ".metadata.toml"),
        body_pvd=joinpath(path, PITCHING_WING_NAME * "_body1.pvd"),
        wake_pvd=joinpath(path, PITCHING_WING_NAME * "_wake1.pvd"),
        csv=_pitching_wing_csv_path(path),
        config=_pitching_wing_config_path(path),
        pressure_state=joinpath(path, PITCHING_WING_PRESSURE_STATE_FILE),
    )
    present = (; (key => isfile(value) for (key, value) in pairs(files))...)
    return (; files, present, any=any(values(present)), restartable=all(values(present)))
end

function _prompt_pitching_wing_restart(input::IO, output::IO)
    while true
        print(output, "Restart from existing data? (y/n) ")
        flush(output)
        eof(input) && throw(ArgumentError(
            "EOF while choosing pitching-wing restart mode. Pass check_existing=false " *
            "to clear and run fresh, or provide interactive input to restart."))
        answer = lowercase(strip(readline(input)))
        answer in ("y", "yes") && return true
        answer in ("n", "no") && return false
        println(output, "Please enter y or n.")
    end
end

_choose_pitching_wing_restart(check_existing, state, input, output) =
    check_existing && state.any ? _prompt_pitching_wing_restart(input, output) : false

function _latest_pitching_wing_step(path)
    state = _pitching_wing_existing_state(path)
    state.restartable || throw(ArgumentError(
        "Existing pitching-wing output is incomplete and cannot be restarted " *
        "(required: metadata, body/wake PVD, histories, and config). Start fresh with check_existing=false."))
    data = TOML.parsefile(state.files.metadata)
    steps = Int[Int(step["i_step"]) for step in get(data, "step", Any[])]
    isempty(steps) && throw(ArgumentError("Pitching-wing metadata contains no saved steps; start fresh with check_existing=false."))
    _, body_steps = pnl._read_pvd_steps(state.files.body_pvd)
    _, wake_steps = pnl._read_pvd_steps(state.files.wake_pvd)
    common = intersect(steps, body_steps, wake_steps)
    isempty(common) && throw(ArgumentError("Pitching-wing metadata and VTK manifests have no common saved step."))
    return maximum(common)
end

function _load_pitching_wing_history(path, t_range, period, section_eta, restart_step)
    filename = _pitching_wing_csv_path(path)
    lines = readlines(filename)
    isempty(lines) && throw(ArgumentError("Pitching-wing CSV is empty: $(filename)"))
    expected = vcat(["time", "t_over_T", "alpha_deg", "CL"], _section_labels(section_eta))
    header = split(strip(lines[1]), ',')
    header == expected || throw(ArgumentError(
        "Pitching-wing CSV header does not match the requested sections. " *
        "Expected $(join(expected, ",")); found $(join(header, ","))."))
    nrows = restart_step + 1
    length(lines) - 1 >= nrows || throw(ArgumentError(
        "Pitching-wing CSV has $(length(lines)-1) samples but restart step $(restart_step) requires $(nrows)."))
    values = Matrix{Float64}(undef, length(expected), nrows)
    for i in 1:nrows
        fields = split(strip(lines[i + 1]), ',')
        length(fields) == length(expected) || throw(ArgumentError("Malformed pitching-wing CSV row $(i + 1)."))
        try
            values[:, i] .= parse.(Float64, fields)
        catch err
            throw(ArgumentError("Malformed numeric value in pitching-wing CSV row $(i + 1): $(sprint(showerror, err))"))
        end
        isapprox(values[1, i], t_range[i]; atol=100eps(Float64) * max(abs(t_range[i]), 1.0), rtol=1e-10) ||
            throw(ArgumentError("Pitching-wing CSV time at sample $(i) is incompatible with the requested timestep."))
        isapprox(values[2, i], t_range[i] / period; atol=1e-10, rtol=1e-10) ||
            throw(ArgumentError("Pitching-wing CSV t_over_T at sample $(i) is inconsistent with time."))
    end
    return (;
        time=collect(values[1, :]),
        alpha=collect(values[3, :]),
        CL=collect(values[4, :]),
        section=copy(values[5:end, :]),
    )
end

function _write_pitching_wing_pressure_state(path, monitor, i_step)
    history = only(monitor.potential_history)
    filename = joinpath(path, PITCHING_WING_PRESSURE_STATE_FILE)
    open(filename, "w") do io
        println(io, "i_step,potential_history")
        for value in history
            println(io, "$(i_step),$(value)")
        end
    end
    return filename
end

function _seed_pitching_wing_pressure_history!(monitor, body, path, restart_step, t_range)
    pnl._load_body_vtk!(body, path, PITCHING_WING_NAME * "_body1", restart_step)
    pnl._pressure_bernoulli_ensure_storage!(monitor, (body,))
    filename = joinpath(path, PITCHING_WING_PRESSURE_STATE_FILE)
    lines = readlines(filename)
    length(lines) == body.ncells + 1 || throw(ArgumentError(
        "Pitching-wing pressure state has $(length(lines)-1) panels; expected $(body.ncells)."))
    strip(lines[1]) == "i_step,potential_history" || throw(ArgumentError(
        "Malformed pitching-wing pressure-state header in $(filename)."))
    for (i, line) in enumerate(lines[2:end])
        fields = split(strip(line), ',')
        length(fields) == 2 || throw(ArgumentError("Malformed pitching-wing pressure-state row $(i+1)."))
        saved_step = parse(Int, fields[1])
        saved_step == restart_step || throw(ArgumentError(
            "Pitching-wing pressure state is for step $(saved_step), but VTK restart step is $(restart_step)."))
        monitor.potential_history[1][i] = parse(Float64, fields[2])
    end
    # Restart safely with one valid sample: the next call uses backward Euler,
    # then BDF2 resumes. Use the known interval after the restart sample rather
    # than the truncated run's final-step fallback dt.
    monitor.history_count[1] = 1
    monitor.previous_dt[1] = t_range[restart_step + 2] - t_range[restart_step + 1]
    monitor.older_dt[1] = NaN
    monitor.last_step[1] = restart_step
    return monitor
end

function plot_pitching_wing_convergence(t_range, period, CL_history, section_eta,
                                         section_history; path=DEFAULT_PITCHING_WING_PATH)
    plt = Core.eval(@__MODULE__, :(begin
        import PythonPlot as pitching_wing_convergence_plt
        pitching_wing_convergence_plt
    end))
    return Base.invokelatest(_plot_pitching_wing_convergence, plt, collect(t_range) ./ period,
        collect(CL_history), collect(section_eta), section_history, path)
end

function _pitching_wing_convergence_data(t_range, period, CL_history, section_eta, section_history)
    size(section_history, 1) == length(section_eta) || throw(DimensionMismatch(
        "section history has $(size(section_history, 1)) rows for $(length(section_eta)) requested eta values"))
    size(section_history, 2) == length(t_range) || throw(DimensionMismatch(
        "section history has $(size(section_history, 2)) samples for $(length(t_range)) times"))
    length(CL_history) == length(t_range) || throw(DimensionMismatch(
        "CL history has $(length(CL_history)) samples for $(length(t_range)) times"))
    return (;
        cycles=collect(t_range) ./ period,
        CL=collect(CL_history),
        sections=[(; eta=float(eta), cl=collect(section_history[i, :]))
                  for (i, eta) in pairs(section_eta)],
    )
end

function _plot_pitching_wing_convergence(plt, cycles, CL_history, section_eta, section_history, path)
    data = _pitching_wing_convergence_data(cycles, 1.0, CL_history, section_eta, section_history)
    mkpath(path)
    CL_path = joinpath(path, "pitching_wing_CL_vs_cycle.png")
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(data.cycles, data.CL; color="C0", linewidth=1.3)
    ax.set_xlabel(L"t/T")
    ax.set_ylabel(L"C_L")
    ax.grid(true, alpha=0.35)
    fig.tight_layout()
    fig.savefig(CL_path, dpi=170)
    plt.pyplot.close(fig)

    section_path = joinpath(path, "pitching_wing_section_cl_vs_cycle.png")
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for (i, series) in pairs(data.sections)
        ax.plot(data.cycles, series.cl; color="C$(i - 1)", linewidth=1.3,
            label=L"$\eta = %$(round(series.eta, sigdigits=4))$")
    end
    ax.set_xlabel(L"t/T")
    ax.set_ylabel(L"c_\ell")
    ax.grid(true, alpha=0.35)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(section_path, dpi=170)
    plt.pyplot.close(fig)
    return (; CL=CL_path, section_cl=section_path)
end

# Split a run of α values into monotonic sweeps (mirrors
# `_pitching_wing_exp_sweep_ranges` in data/pitching_wing_exp/load.jl so the
# example matches the experimental hysteresis plots without coupling to it).
function _alpha_sweep_ranges(alpha)
    n = length(alpha)
    n <= 1 && return [1:max(n, 1)]
    ranges = UnitRange{Int}[]
    start_idx = 1
    last_sign = 0
    for i in 2:n
        step = alpha[i] - alpha[i - 1]
        sign_step = step > 0 ? 1 : step < 0 ? -1 : 0
        sign_step == 0 && continue
        if last_sign != 0 && sign_step != last_sign
            push!(ranges, start_idx:(i - 1))
            start_idx = i - 1
        end
        last_sign = sign_step
    end
    push!(ranges, start_idx:n)
    return ranges
end

# Indices of the last `n_loops` full cycles, always excluding the first
# (startup) cycle: matched cycles superimpose only once the limit cycle has
# settled, so a clean loop is itself evidence of convergence.
function _pitching_wing_last_cycles_window(t_range, period, n_loops)
    n_loops >= 1 || throw(ArgumentError("n_loops must be positive; got $(n_loops)"))
    t_end = last(t_range)
    tol = 1e-9 * max(abs(t_end), 1.0)
    lower = max(period, t_end - n_loops * period) - tol
    return findall(t -> t >= lower, t_range)
end

# α-hysteresis loops for the last `n_loops` full cycles (startup cycle always
# excluded), styled like the experimental curves in
# data/pitching_wing_exp/load.jl: one color per η, solid = increasing α /
# dashed = decreasing α. The color/η legend lives in the left panel and the
# black solid/dashed direction key in the right panel. Overlaying the final two
# cycles (the default) lets a coincident loop double as a convergence check.
function plot_pitching_wing_hysteresis(t_range, period, alpha_history, CL_history,
        section_eta, section_history; path=DEFAULT_PITCHING_WING_PATH, n_loops::Integer=2)
    plt = Core.eval(@__MODULE__, :(begin
        import PythonPlot as pitching_wing_hysteresis_plt
        pitching_wing_hysteresis_plt
    end))
    return Base.invokelatest(_plot_pitching_wing_hysteresis, plt, collect(t_range),
        float(period), collect(alpha_history), collect(CL_history),
        collect(section_eta), section_history, path, Int(n_loops))
end

function _plot_pitching_wing_hysteresis(plt, t_range, period, alpha_history,
        CL_history, section_eta, section_history, path, n_loops)
    mkpath(path)
    window = _pitching_wing_last_cycles_window(t_range, period, n_loops)
    alpha = alpha_history[window]
    CL = CL_history[window]
    sweeps = _alpha_sweep_ranges(alpha)

    hysteresis_path = joinpath(path, "pitching_wing_hysteresis.png")
    fig, axs = plt.subplots(1, 2; figsize=(11.0, 4.6), squeeze=false)
    ax_section = axs[0, 0]
    ax_CL = axs[0, 1]

    for (i, eta) in pairs(section_eta)
        cl = section_history[i, window]
        for range in sweeps
            increasing = length(range) < 2 || alpha[last(range)] >= alpha[first(range)]
            ax_section.plot(alpha[range], cl[range];
                color="C$(i - 1)", linestyle=increasing ? "-" : "--",
                marker="o", markersize=3, label=nothing)
        end
    end
    ax_section.set_xlabel(L"$\alpha$ (deg)")
    ax_section.set_ylabel(L"c_\ell")
    ax_section.grid(true, alpha=0.35)
    ax_section.legend(handles=[plt.matplotlib.lines.Line2D([], [];
        color="C$(i - 1)", marker="o", markersize=3,
        label=L"$\eta = %$(round(eta, sigdigits=4))$")
        for (i, eta) in pairs(section_eta)], fontsize=9)

    for range in sweeps
        increasing = length(range) < 2 || alpha[last(range)] >= alpha[first(range)]
        ax_CL.plot(alpha[range], CL[range];
            color="black", linestyle=increasing ? "-" : "--",
            marker="o", markersize=3, label=nothing)
    end
    ax_CL.set_xlabel(L"$\alpha$ (deg)")
    ax_CL.set_ylabel(L"C_L")
    ax_CL.grid(true, alpha=0.35)
    ax_CL.legend(handles=[
        plt.matplotlib.lines.Line2D([], []; color="black", linestyle="-",
            label=L"increasing $\alpha$"),
        plt.matplotlib.lines.Line2D([], []; color="black", linestyle="--",
            label=L"decreasing $\alpha$"),
    ], fontsize=9)

    fig.tight_layout()
    fig.savefig(hysteresis_path, dpi=170)
    plt.pyplot.close(fig)
    return hysteresis_path
end

mutable struct SectionLiftHistoryMonitor{TF, TS} <: pnl.AbstractMonitor
    spanwise::TS
    section_eta::Vector{TF}
    semispan::TF
    chord::TF
    rho::TF
    cl::Matrix{TF}
end

monitor_requires(::SectionLiftHistoryMonitor) = (:sectional_F,)
monitor_provides(::SectionLiftHistoryMonitor) = ()

function SectionLiftHistoryMonitor(spanwise, nt::Integer, section_eta, semispan, chord, rho)
    eta = collect(float.(section_eta))
    return SectionLiftHistoryMonitor(spanwise, eta, float(semispan), float(chord), float(rho),
        fill(NaN, length(eta), nt))
end

function _run_monitor!(m::SectionLiftHistoryMonitor, ctx::pnl.MonitorContext, systems, wakes,
                       frames::AbstractVector{<:pnl.ReferenceFrame}, uinf,
                       i_step::Int, dt::Real, t=nothing)
    speed2 = uinf[1]^2 + uinf[2]^2 + uinf[3]^2
    qinf = 0.5 * m.rho * speed2
    bin_eta = abs.(m.spanwise.bin_center) ./ m.semispan
    cl_bins = m.spanwise.load_components[1, :] ./ (qinf * m.chord)
    m.cl[:, i_step + 1] .= _interp_section_values(bin_eta, cl_bins, m.section_eta)
    return nothing
end

function _finite_summary(values)
    finite = filter(isfinite, collect(values))
    isempty(finite) && return "all NaN/Inf"
    return "min=$(round(minimum(finite), sigdigits=5)), max=$(round(maximum(finite), sigdigits=5)), last=$(round(finite[end], sigdigits=5))"
end

function run_pitching_wing_static_polar(;
        semispan_ft::Union{Nothing,Real}=nothing,
        span_ft::Union{Nothing,Real}=nothing,
        aspect_ratio::Real=4.0,
        c_ft::Real=1.0,
        v_inf_ft_s::Real=330.2,
        rho::Real=1.225,
        static_alpha_deg=0.0:1.0:10.0,
        section_eta=[0.25, 0.5, 0.75],
        n_span::Integer=13,
        n_airfoil::Integer=161,
        n_endcap::Integer=9,
        n_section_bins::Integer=max(2 * n_span, 20),
        endcap::Symbol=:round,
        path=DEFAULT_PITCHING_WING_PATH,
        save_vtk::Bool=true,
        plot::Bool=false,
        plot_path=path,
        backend=pnl.FastMultipoleBackend(expansion_order=8, multipole_acceptance=0.4, leaf_size=40),
    )
    dims = _resolve_pitching_wing_dimensions(c_ft, aspect_ratio, span_ft, semispan_ft)
    c, b = dims.c, dims.b
    U = v_inf_ft_s * FT_TO_M
    Sref = b * c
    rows = NamedTuple[]
    alpha_values = collect(static_alpha_deg)
    vtk_path = save_vtk ? path : nothing

    for (i_alpha, alpha_deg) in enumerate(alpha_values)
        body = build_pitching_wing_body(c, b;
            n_span, n_airfoil, n_endcap, endcap, semiinfinite_wake=true)
        Uinf = _uinf_from_alpha(U, alpha_deg)
        Lhat, Dhat, Shat = _lift_drag_span_directions(Uinf)
        set_wake_Das!(body, Uinf)

        pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
        force = pnl.ForceMonitor(length(alpha_values), 1;
            i_frame=-1,
            normalization=pnl.WingNormalization(rho, Sref, c),
            correct_kuttacondition=false,
            verbose=false,
        )
        spanwise = pnl.SpanwiseLoadingMonitor(n_section_bins, 1;
            components=(lift=Lhat, drag=Dhat),
            span_axis=Shat,
            per_length=true,
            normalization=pnl.NoSectionalNormalization(),
            file=false,
            binning=:span_overlap,
        )

        pnl.steady!(body, pnl.ReferenceFrame(body), Uinf;
            body_solvers=pnl.Backslash(body),
            backend,
            monitors=(pressure, force, spanwise),
            path=vtk_path,
            name="pitching_wing_static",
            i_run=i_alpha,
            verbose=false,
        )

        qinf = 0.5 * rho * U^2
        bin_eta = abs.(spanwise.bin_center) ./ (b / 2)
        cl_bins = spanwise.load_components[1, :] ./ (qinf * c)
        i_alpha == firstindex(alpha_values) &&
            _print_section_interpolation_bins(stdout, bin_eta, section_eta)
        section_cl = _interp_section_values(bin_eta, cl_bins, section_eta)
        CL = dot(SVector{3}(force.force[:, i_alpha]), Lhat)
        push!(rows, (; alpha_deg=float(alpha_deg), CL, section_cl))
        @printf("  static alpha = %6.2f deg: CL = %.8g\n", alpha_deg, CL)
    end

    csv_path = _write_static_polar_csv(path, section_eta, rows)
    println("  Wrote static polar CSV: $(csv_path)")
    save_vtk && println("  Wrote static VTK output under: $(path)")
    plot_paths = if plot
        paths = plot_pitching_wing_static_polar(rows, section_eta; path=plot_path)
        println("  Wrote static polar plots: $(paths.section_cl), $(paths.CL)")
        paths
    else
        ()
    end
    return (; rows, csv_path, plot_paths)
end

function prepare_pitching_wing(;
        semispan_ft::Union{Nothing,Real}=nothing,
        span_ft::Union{Nothing,Real}=nothing,
        aspect_ratio::Real=4.0,
        c_ft::Real=1.0,
        v_inf_ft_s::Real=330.2,
        rho::Real=1.225,
        alpha_mean_deg::Real=3.94,
        alpha_amp_deg::Real=1.99,
        frequency_hz::Real=4.01,
        n_cycles::Real=3,
        c_per_dt::Real=0.5,
        das_chord_fraction::Real=0.05,
        pivot_chord_fraction::Real=0.25,
        section_eta=[0.25, 0.5, 0.75],
        static_alpha_deg=0.0:1.0:10.0,
        n_span::Integer=13,
        n_airfoil::Integer=161,
        n_endcap::Integer=9,
        n_section_bins::Integer=max(2 * n_span, 20),
        panel_wake_rows::Union{Nothing,Integer}=nothing,
        max_particles::Integer=20000,
        wake_length_spans::Union{Nothing,Real}=2.0,
        wake_model::Symbol=:panel,
        endcap::Symbol=:round,
        save_vtk::Bool=true,
        path=DEFAULT_PITCHING_WING_PATH,
        include_static_polar::Bool=true,
        plot_static_polar::Bool=false,
        backend=pnl.FastMultipoleBackend(expansion_order=8, multipole_acceptance=0.4, leaf_size=40),
    )
    wake_model in (:panel, :particle) || throw(ArgumentError(
        "wake_model must be :panel or :particle; got $(wake_model)"))
    dims = _resolve_pitching_wing_dimensions(c_ft, aspect_ratio, span_ft, semispan_ft)
    c, b = dims.c, dims.b
    if !isnothing(wake_length_spans)
        isfinite(wake_length_spans) && wake_length_spans > 0 || throw(ArgumentError(
            "wake_length_spans must be positive and finite, or nothing to disable trimming; " *
            "got $(repr(wake_length_spans))."))
    end
    wake_length = isnothing(wake_length_spans) ? nothing : wake_length_spans * b
    wake_downstream_boundary = isnothing(wake_length) ? nothing :
        c * (1 - pivot_chord_fraction) + wake_length
    U = v_inf_ft_s * FT_TO_M
    Sref = b * c
    omega = 2 * pi * frequency_hz
    period = 1 / frequency_hz
    dt = c / U * c_per_dt
    t_range = range(0.0, stop=n_cycles * period, step=dt)
    length(t_range) >= 2 || throw(ArgumentError(
        "pitching wing simulation requires at least two time samples; got $(length(t_range))."))

    pivot = SVector{3}(pivot_chord_fraction * c, 0.0, 0.0)
    wing = build_pitching_wing_body(c, b;
        n_span, n_airfoil, n_endcap, endcap, semiinfinite_wake=false)

    alpha_mean = deg2rad(alpha_mean_deg)
    alpha_amp = deg2rad(alpha_amp_deg)
    frames = pitching_wing_frame(wing, pivot, alpha_mean)
    Uinf(t) = SVector{3}(U, 0.0, 0.0)
    das = das_chord_fraction * c
    resolved_panel_wake_rows = isnothing(panel_wake_rows) ?
        (isnothing(wake_length) ? 100 : max(1, ceil(Int, wake_length / das))) :
        Int(panel_wake_rows)
    resolved_panel_wake_rows >= 1 || throw(ArgumentError(
        "panel_wake_rows must be positive; got $(panel_wake_rows)"))
    set_wake_Das!(wing, Uinf(first(t_range)); magnitude=das)
    Lhat = SVector{3}(0.0, 0.0, 1.0)
    Dhat = SVector{3}(1.0, 0.0, 0.0)
    Shat = SVector{3}(0.0, 1.0, 0.0)

    function pitching_maneuver!(frames, systems, wakes, t)
        frame = frames[1]
        this_omega = alpha_amp * omega * cos(omega * t)
        frames[1] = typeof(frame)(
            frame.x,
            frame.v,
            frame.ω_axis,
            this_omega,
            frame.R,
            frame.Rp2g,
            frame.name,
            frame.parent_index,
            frame.child_index,
            frame.dependent_index,
        )
        return nothing
    end

    wake = if wake_model == :panel
        # A ConstantDoublet PanelWake has a scalar potential, making unsteady
        # Bernoulli a complete irrotational reference for pressure validation.
        pnl.PanelWake(wing; nwakerows=resolved_panel_wake_rows,
            include_final_filament=false)
    else
        particle_maintenance = isnothing(wake_downstream_boundary) ?
            pnl.ParticleMaintenance() :
            pnl.ParticleMaintenance(pnl.FrameBox(1,
                SVector(-Inf, -Inf, -Inf),
                SVector(wake_downstream_boundary, Inf, Inf)))
        pnl.PanelParticleWake(wing;
            nwakerows=resolved_panel_wake_rows,
            max_particles=max_particles,
            method_trailing=pnl.OverlapPPS(1.3, 2),
            method_unsteady=pnl.OverlapPPS(1.3, 2),
            particle_maintenance,
        )
    end
    solver = pnl.Backslash(wing)
    normalization = pnl.WingNormalization(rho, Sref, c)

    pressure_monitor = pnl.PressureBernoulli(rho; unsteady=true,
        allow_partial=(wake_model == :particle))
    force_monitor = pnl.ForceMonitor(length(t_range), 1;
        normalization,
        i_frame=-1,
        verbose=false,
    )
    spanwise_monitor = pnl.SpanwiseLoadingMonitor(n_section_bins, 1;
        components=(lift=Lhat, drag=Dhat),
        span_axis=Shat,
        per_length=true,
        normalization=pnl.NoSectionalNormalization(),
        file=false,
        binning=:span_overlap,
    )
    section_monitor = SectionLiftHistoryMonitor(spanwise_monitor, length(t_range),
        section_eta, b / 2, c, rho)
    monitors = (pressure_monitor, force_monitor, spanwise_monitor, section_monitor)

    return (;
        wing,
        wake,
        frames,
        maneuver! = pitching_maneuver!,
        Uinf,
        t_range,
        solver,
        backend,
        monitors,
        pressure_monitor,
        force_monitor,
        spanwise_monitor,
        section_monitor,
        static_polar = include_static_polar,
        setup = (;
            semispan_ft=dims.semispan_ft,
            span_ft=dims.span_ft,
            aspect_ratio=dims.aspect_ratio,
            c_ft,
            v_inf_ft_s,
            rho,
            alpha_mean_deg,
            alpha_amp_deg,
            frequency_hz,
            n_cycles,
            c_per_dt,
            das_chord_fraction,
            das,
            pivot_chord_fraction,
            section_eta,
            static_alpha_deg,
            n_span,
            n_airfoil,
            n_endcap,
            n_section_bins,
            panel_wake_rows=resolved_panel_wake_rows,
            max_particles,
            wake_model,
            wake_length_spans,
            wake_length,
            wake_downstream_boundary,
            endcap,
            c,
            b,
            Sref,
            reference_area=Sref,
            U,
            omega,
            period,
            dt,
            pivot,
            path = save_vtk ? path : nothing,
            csv_path = path,
            plot_static_polar,
        ),
    )
end

function run_pitching_wing(; wake_length_spans::Union{Nothing,Real}=2.0, kwargs...)
    return _run_pitching_wing(; wake_length_spans, kwargs...)
end

function _run_pitching_wing(;
        check_existing::Bool=true,
        plot_convergence::Bool=true,
        restart_input::IO=stdin,
        restart_output::IO=stdout,
        kwargs...)
    sim = prepare_pitching_wing(; kwargs...)
    (; wing, wake, frames, maneuver!, Uinf, t_range, solver, backend, monitors,
       pressure_monitor, force_monitor, section_monitor, setup) = sim

    state = _pitching_wing_existing_state(setup.csv_path)
    restart = false
    restart_step = nothing
    prior = nothing
    requested_config = _pitching_wing_config(setup)
    restart = _choose_pitching_wing_restart(check_existing, state,
        restart_input, restart_output)
    if restart
        _validate_pitching_wing_config(setup.csv_path, requested_config)
        isnothing(setup.path) && throw(ArgumentError(
            "Cannot restart when save_vtk=false; rerun fresh with check_existing=false."))
        restart_step = _latest_pitching_wing_step(setup.csv_path)
        restart_step + 1 < length(t_range) || throw(ArgumentError(
            "Requested n_cycles=$(setup.n_cycles) does not extend beyond saved step $(restart_step). " *
            "Increase n_cycles to continue the run."))
        prior = _load_pitching_wing_history(setup.csv_path, t_range, setup.period,
            setup.section_eta, restart_step)
        force_monitor.force[3, 1:restart_step+1] .= prior.CL
        section_monitor.cl[:, 1:restart_step+1] .= prior.section
        _seed_pitching_wing_pressure_history!(pressure_monitor, wing, setup.path, restart_step, t_range)
    else
        _clear_pitching_wing_output!(setup.csv_path)
        _write_pitching_wing_config(setup.csv_path, requested_config)
    end

    println("Pitching rectangular NACA 0015 wing")
    @printf("  c = %.6g m (%.4g ft), b = %.6g m (%.4g ft), AR = %.4g, S = %.6g m^2\n",
        setup.c, setup.c_ft, setup.b, 2 * setup.semispan_ft, setup.b / setup.c, setup.Sref)
    @printf("  U = %.6g m/s (%.4g ft/s), rho = %.6g kg/m^3\n",
        setup.U, setup.v_inf_ft_s, setup.rho)
    @printf("  alpha(t) = %.4g + %.4g sin(2*pi*%.4g*t) deg, pivot = %.4g c\n",
        setup.alpha_mean_deg, setup.alpha_amp_deg, setup.frequency_hz, setup.pivot_chord_fraction)
    @printf("  period = %.6g s, dt = %.6g s, steps = %d, steps/cycle ~= %.1f\n",
        setup.period, setup.dt, length(t_range), setup.period / setup.dt)
    println("  panels = $(wing.ncells), nodes = $(wing.nnodes), shedding edges = $(wing.nsheddings)")

    static = sim.static_polar ? run_pitching_wing_static_polar(;
        semispan_ft=setup.semispan_ft,
        c_ft=setup.c_ft,
        v_inf_ft_s=setup.v_inf_ft_s,
        rho=setup.rho,
        section_eta=setup.section_eta,
        static_alpha_deg=setup.static_alpha_deg,
        n_span=setup.n_span,
        n_airfoil=setup.n_airfoil,
        n_endcap=setup.n_endcap,
        n_section_bins=setup.n_section_bins,
        endcap=setup.endcap,
        path=setup.csv_path,
        save_vtk=!isnothing(setup.path),
        plot=setup.plot_static_polar,
        plot_path=setup.csv_path,
        backend,
    ) : nothing

    if restart
        pnl.simulate_warmstart!((wing,), (wake,), frames, maneuver!, Uinf, t_range;
            body_solvers=(solver,), backend, monitors, path=setup.path,
            name=PITCHING_WING_NAME, restart_path=setup.path,
            restart_name=PITCHING_WING_NAME, restart_step,
            set_Das_eta_freestream=NaN, clean_files=false, verbose=true)
    else
        pnl.simulate!((wing,), (wake,), frames, maneuver!, Uinf, t_range;
            body_solvers=(solver,),
            backend,
            monitors,
            path=setup.path,
            name=PITCHING_WING_NAME,
            set_Das_eta_freestream=NaN,
            verbose=true,
        )
    end

    bin_eta = abs.(sim.spanwise_monitor.bin_center) ./ (setup.b / 2)
    _print_section_interpolation_bins(stdout, bin_eta, setup.section_eta)

    alpha_history = [setup.alpha_mean_deg + setup.alpha_amp_deg * sin(setup.omega * t) for t in t_range]
    CL_history = force_monitor.force[3, :]
    section_history = section_monitor.cl

    csv_dir = setup.csv_path
    csv_path = _write_unsteady_csv(csv_dir, t_range, setup.period, alpha_history,
        CL_history, setup.section_eta, section_history)
    _write_pitching_wing_pressure_state(csv_dir, pressure_monitor, length(t_range) - 1)
    plot_paths = plot_convergence ? plot_pitching_wing_convergence(t_range, setup.period,
        CL_history, setup.section_eta, section_history; path=csv_dir) : ()
    hysteresis_path = plot_convergence ? plot_pitching_wing_hysteresis(t_range,
        setup.period, alpha_history, CL_history, setup.section_eta, section_history;
        path=csv_dir) : nothing

    post_transient = findall(t -> t >= setup.period, collect(t_range))
    CL_mean = isempty(post_transient) ? NaN : mean(CL_history[post_transient])
    CL_pp = isempty(post_transient) ? NaN :
        maximum(CL_history[post_transient]) - minimum(CL_history[post_transient])

    println("\nPitching-wing force summary")
    @printf("  CL mean after first cycle = %.8g\n", CL_mean)
    @printf("  CL peak-to-peak after first cycle = %.8g\n", CL_pp)
    println("  CL history: $(_finite_summary(CL_history))")
    println("  Wrote unsteady CSV: $(csv_path)")
    if plot_convergence
        println("  Wrote convergence plots: $(plot_paths.CL), $(plot_paths.section_cl)")
        println("  Wrote hysteresis plot: $(hysteresis_path)")
    end
    if !isnothing(setup.path)
        println("  Wrote VTK output under: $(setup.path)")
    end

    return merge(sim, (; static, alpha_history, CL_history, section_history,
                        CL_mean, CL_pp, csv_path, convergence_plot_paths=plot_paths,
                        hysteresis_plot_path=hysteresis_path,
                        restarted=restart, restart_step))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_pitching_wing()
end
