import FLOWPanel as pnl
using DelimitedFiles
using LinearAlgebra
using Printf
using PythonPlot
using Statistics

const plt = PythonPlot.matplotlib.pyplot

const SCRIPT_DIR = @__DIR__
const INPUT_MESH = joinpath(SCRIPT_DIR, "dji9443_medium.msh")
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "processed")
const RUN_NAME = "dji9443_medium"

# Fill these after inspecting `processed/dji9443_medium_raw.pvd` in Paraview.
const plane_idx = Int[338480, 346889, 355054, 317232, 314981, 322126, 326943] .+ 1
const radial_idx = 5894 + 1

const NSLICES = 15
const SLICE_MIN_POINTS = 400
const SLICE_DISTANCE_ATOL = 0.0
const N_RFL_SECS = 100
const USE_COSINE_RFL_BINS = true
const HUB_TOP_FRACTION = 0.75
const CST_NCOEFFS = 8
const CST_RECONSTRUCTION_POINTS = 201
const SHARP_TE_BLEND_START = 0.8

cst_output_path(output_dir=OUTPUT_DIR) = joinpath(output_dir, RUN_NAME * "_cst.csv")

function section_review_paths(index::Int; output_dir=OUTPUT_DIR)
    directory = joinpath(output_dir, @sprintf("section_%02d_airfoil_review", index))
    prefix = @sprintf("section_%02d", index)
    return (; directory,
            raw=joinpath(directory, prefix * "_raw.dat"),
            airfoil=joinpath(directory, prefix * ".dat"),
            leading_edge=joinpath(directory, prefix * "_leading_edge_row_guide.png"),
            raw_selection=joinpath(directory, prefix * "_checkpoint_1_raw_selection.png"),
            normalization=joinpath(
                directory, prefix * "_checkpoint_1_normalization_and_canonicalization.png"),
            blunt=joinpath(directory, prefix * "_checkpoint_2_blunt_cst.png"),
            taper=joinpath(directory, prefix * "_checkpoint_3_taper.png"),
            final=joinpath(directory, prefix * "_checkpoint_4_final_cst.png"))
end

"""
    MeshSlice 344, 21, 343

2-D blade section at radial coordinate `r`.

`nodes` is a 2xN matrix in raw section coordinates before airfoil
normalization. Row 1 is axial coordinate `x`; row 2 is thickness coordinate
`z` in the aligned rotor frame.
"""
struct MeshSlice
    r::Float64
    nodes::Matrix{Float64}
end

struct AnchoredMeshSlice
    r::Float64
    nodes::Matrix{Float64}
    upper_te::NTuple{2,Float64}
    le_idx::Int
    lower_te::NTuple{2,Float64}
    excluded_indices::Vector{Int}
    blend_start::Float64
end

"""
    Airfoil

Chord-normalized airfoil extracted from a `MeshSlice`.

Coordinates are stored as upper and lower surfaces in leading-edge to
trailing-edge order.
"""
struct Airfoil
    r::Float64
    twist::Float64
    chord::Float64
    axial::Float64
    tangential::Float64
    x_u::Vector{Float64}
    y_u::Vector{Float64}
    x_l::Vector{Float64}
    y_l::Vector{Float64}
end

normalize_vec(v) = v / norm(v)

function require_indices(nodes)
    if isempty(plane_idx)
        error("""
        Fill `plane_idx` near the top of this file.

        A raw mesh preview has been written to:
          $(joinpath(OUTPUT_DIR, RUN_NAME * "_raw.pvd"))

        Pick several node indices that lie in the rotor plane, set
        `const plane_idx = [...]`, and rerun this script.
        """)
    end
    if radial_idx <= 0
        error("""
        Fill `radial_idx` near the top of this file.

        Pick one node on the desired positive radial axis in the raw/aligned
        mesh preview, set `const radial_idx = <index>`, and rerun this script.
        """)
    end
    all(1 .<= plane_idx .<= size(nodes, 2)) ||
        error("`plane_idx` contains indices outside 1:$(size(nodes, 2)).")
    1 <= radial_idx <= size(nodes, 2) ||
        error("`radial_idx` must be in 1:$(size(nodes, 2)).")
end

function load_body(mesh_path)
    msh = pnl.read_gmsh(mesh_path)
    nodes, cells = pnl.meshes2nodes_cells(msh)
    body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(nodes), copy(cells); watertight=false)
    return nodes, cells, body
end

function write_body_vtk(pathstem, nodes, cells; overwrite=true)
    body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(nodes), copy(cells); watertight=false)
    return pnl.write_vtk(pathstem, body, 0, 0.0; overwrite)
end

function fit_plane(nodes, idx)
    pts = nodes[:, idx]
    point = vec(mean(pts; dims=2))
    centered = pts .- reshape(point, :, 1)
    U, _, _ = svd(centered)
    normal = normalize_vec(U[:, end])
    return point, normal
end

function alignment_from_manual_points(nodes)
    plane_point, plane_normal = fit_plane(nodes, plane_idx)
    signed_distance = vec((nodes .- reshape(plane_point, :, 1))' * plane_normal)

    max_positive = maximum(signed_distance)
    max_negative = abs(minimum(signed_distance))
    top_sign = max_positive <= max_negative ? 1.0 : -1.0
    top_distance = min(max_positive, max_negative)
    top_distance > 0.0 || error("Could not identify a nonzero hub-top distance.")

    hub_mask = signed_distance .* top_sign .>= HUB_TOP_FRACTION * top_distance
    any(hub_mask) || error("No points found on the top hub side with the configured threshold.")

    axis_point = vec(mean(nodes[:, hub_mask]; dims=2))
    origin = axis_point - dot(axis_point - plane_point, plane_normal) * plane_normal

    xhat = normalize_vec(-top_sign * plane_normal)
    radial_vector = nodes[:, radial_idx] - origin
    yhat = radial_vector - dot(radial_vector, xhat) * xhat
    norm(yhat) > 0.0 || error("`radial_idx` lies on the rotation axis; choose an off-axis point.")
    yhat = normalize_vec(yhat)
    zhat = normalize_vec(cross(xhat, yhat))
    yhat = normalize_vec(cross(zhat, xhat))

    basis = hcat(xhat, yhat, zhat)
    aligned_nodes = basis' * (nodes .- reshape(origin, :, 1))
    return aligned_nodes, origin, basis
end

function estimate_slice_width(rnodes)
    r_sorted = sort(unique(rnodes))
    length(r_sorted) < 2 && return 0.0
    gaps = diff(r_sorted)
    positive_gaps = gaps[gaps .> eps(Float64)]
    isempty(positive_gaps) ? 0.0 : median(positive_gaps)
end

function generate_slices(nodes; nslices=NSLICES, min_points=SLICE_MIN_POINTS,
                         distance_atol=SLICE_DISTANCE_ATOL)
    min_points > 0 || error("`min_points` must be positive.")
    distance_atol >= 0.0 || error("`distance_atol` must be nonnegative.")

    rnodes = nodes[2, :]
    positive_mask = rnodes .>= 0.0
    any(positive_mask) || error("No nodes found on the positive-y blade side.")

    positive_rnodes = rnodes[positive_mask]
    rmin = minimum(positive_rnodes)
    rmax = maximum(positive_rnodes)
    stations = range(rmin, rmax; length=nslices + 1)
    centers = [(stations[i] + stations[i + 1]) / 2 for i in 1:nslices]
    candidate_indices = findall(positive_mask)

    slices = MeshSlice[]
    for r in centers
        distances = abs.(rnodes[candidate_indices] .- r)
        order = sortperm(distances)
        cutoff = distances[order[min(min_points, length(order))]]
        section_indices = candidate_indices[distances .<= cutoff + distance_atol]
        section_nodes = nodes[:, section_indices]
        section_2d = vcat(-section_nodes[3, :]', -section_nodes[1, :]')
        push!(slices, MeshSlice(float(r), Matrix(section_2d)))
    end

    return slices
end

function transform_points(points, origin, twist, chord)
    c = cos(-twist)
    s = sin(-twist)
    x_shift = points[1, :] .- origin[1]
    y_shift = points[2, :] .- origin[2]
    x = (c .* x_shift .- s .* y_shift) ./ chord
    y = (s .* x_shift .+ c .* y_shift) ./ chord
    return x, y
end

function normalized_section(ms::AnchoredMeshSlice)
    upper_raw = collect(ms.upper_te)
    lower_raw = collect(ms.lower_te)
    upper_raw != lower_raw || error("Upper and lower trailing-edge coordinates must be distinct.")
    te_point = 0.5 .* (upper_raw .+ lower_raw)
    le_point = ms.nodes[:, ms.le_idx]
    chord_vec = te_point - le_point
    chord = norm(chord_vec)
    chord > 0.0 || error("Slice at r=$(ms.r) has zero chord.")

    twist = atan(chord_vec[2], chord_vec[1])
    x_all, y_all = transform_points(ms.nodes, le_point, twist, chord)
    anchor_points = hcat(upper_raw, lower_raw)
    anchor_x, anchor_y = transform_points(anchor_points, le_point, twist, chord)
    anchor_y[1] > anchor_y[2] || error("Upper trailing-edge coordinate must transform above the lower coordinate.")
    h_te = abs(anchor_y[1] - anchor_y[2]) / 2
    h_te > eps(Float64) || error("Trailing-edge thickness is degenerate after chord-normal projection.")

    excluded = falses(size(ms.nodes, 2))
    excluded[ms.excluded_indices] .= true
    retained_before_filter = .!excluded
    chordwise = retained_before_filter .& (x_all .>= 0.0) .& (x_all .<= 1.0)
    any(chordwise) || error("Slice at r=$(ms.r) has no retained points in 0 <= x/c <= 1.")
    retained_indices = findall(chordwise)
    ms.le_idx in retained_indices || error("Selected leading-edge row was discarded by the chordwise filter.")
    dropped_indices = findall(retained_before_filter .& .!chordwise)
    likely_missed = findall(retained_before_filter .&
                            (x_all .> max(anchor_x[1], anchor_x[2])))
    canonical_x = [1.0, 1.0]
    canonical_y = [h_te, -h_te]
    displacements = [hypot(anchor_x[i] - canonical_x[i],
                           anchor_y[i] - canonical_y[i]) for i in 1:2]

    return (; twist, chord, te_point, le_point, x_all, y_all,
            x=x_all[retained_indices], y=y_all[retained_indices], retained_indices,
            dropped_indices, likely_missed, anchor_x, anchor_y,
            canonical_x, canonical_y, x_offsets=anchor_x .- 1.0,
            displacements, h_te)
end

function rotated_section(ms::AnchoredMeshSlice)
    normalized = normalized_section(ms)
    return normalized.twist, normalized.chord, normalized.x, normalized.y,
           normalized.retained_indices, normalized.anchor_x, normalized.anchor_y,
           normalized.canonical_x, normalized.canonical_y
end

function add_unique!(indices, new_indices)
    for idx in new_indices
        idx in indices || push!(indices, idx)
    end
    return indices
end

function chordwise_bin_edges(x; n_sections=N_RFL_SECS, cosine=USE_COSINE_RFL_BINS)
    xmin = minimum(x)
    xmax = maximum(x)
    if cosine
        θ = range(0.0, π; length=n_sections + 1)
        spacing = @. 0.5 * (1.0 - cos(θ))
        return xmin .+ (xmax - xmin) .* spacing
    else
        return range(xmin, xmax; length=n_sections + 1)
    end
end

function sectioned_selig_indices(x, y; n_sections=N_RFL_SECS,
                                 cosine_bins=USE_COSINE_RFL_BINS,
                                 required_upper_indices=Int[],
                                 required_lower_indices=Int[],
                                 upper_te_pos::Int,
                                 lower_te_pos::Int)
    n_sections > 0 || error("`n_sections` must be positive.")

    upper_indices = Int[]
    lower_indices = Int[]
    x_edges = chordwise_bin_edges(x; n_sections=n_sections, cosine=cosine_bins)
    previous_upper_y = y[upper_te_pos]
    previous_lower_y = y[lower_te_pos]

    for i in n_sections:-1:1
        xlo = x_edges[i]
        xhi = x_edges[i + 1]
        in_section = i == n_sections ? findall((x .>= xlo) .& (x .<= xhi)) :
                                        findall((x .>= xlo) .& (x .< xhi))
        isempty(in_section) && continue
        upper = Int[]
        lower = Int[]
        for idx in in_section
            if abs(y[idx] - previous_upper_y) <= abs(y[idx] - previous_lower_y)
                push!(upper, idx)
            else
                push!(lower, idx)
            end
        end
        isempty(upper) || push!(upper_indices, upper[argmax(y[upper])])
        isempty(lower) || push!(lower_indices, lower[argmin(y[lower])])
        isempty(upper) || (previous_upper_y = mean(y[upper]))
        isempty(lower) || (previous_lower_y = mean(y[lower]))
    end

    add_unique!(upper_indices, required_upper_indices)
    add_unique!(lower_indices, required_lower_indices)

    upper_order = sortperm(x[upper_indices]; rev=true)
    lower_order = sortperm(x[lower_indices])
    upper_te_to_le = upper_indices[upper_order]
    lower_le_to_te = lower_indices[lower_order]
    return vcat(upper_te_to_le, lower_le_to_te)
end

function required_selig_indices(original_indices, le_idx, upper_te_pos, lower_te_pos)
    le = findfirst(==(le_idx), original_indices)
    le !== nothing || error("Selected LE row was discarded by the chordwise filter.")
    return [upper_te_pos, le], [le, lower_te_pos]
end

function ordered_surface_samples(ms::AnchoredMeshSlice, normalized=normalized_section(ms))
    x = vcat(normalized.x, normalized.canonical_x)
    y = vcat(normalized.y, normalized.canonical_y)
    upper_te_pos = length(x) - 1
    lower_te_pos = length(x)
    required_upper, required_lower = required_selig_indices(
        normalized.retained_indices, ms.le_idx, upper_te_pos, lower_te_pos)
    selig_indices = sectioned_selig_indices(x, y;
                                            required_upper_indices=required_upper,
                                            required_lower_indices=required_lower,
                                            upper_te_pos=upper_te_pos,
                                            lower_te_pos=lower_te_pos)
    xselig = x[selig_indices]
    yselig = y[selig_indices]
    length(xselig) >= 4 || error("Slice at r=$(ms.r) did not produce enough ordered airfoil points.")
    le_idx = argmin(xselig)
    xu = collect(reverse(xselig[1:le_idx]))
    yu = collect(reverse(yselig[1:le_idx]))
    xl = collect(xselig[le_idx:end])
    yl = collect(yselig[le_idx:end])
    xu[1] = xl[1] = 0.0
    yu[1] = yl[1] = 0.0
    xu[end] = xl[end] = 1.0
    yu[end] = normalized.h_te
    yl[end] = -normalized.h_te
    return (; xu, yu, xl, yl, selig_indices)
end

function selig_coordinates(af::Airfoil)
    lower_start = hypot(af.x_u[1] - af.x_l[1], af.y_u[1] - af.y_l[1]) <=
                  sqrt(eps(Float64)) ? 2 : 1
    x = vcat(reverse(af.x_u), af.x_l[lower_start:end])
    y = vcat(reverse(af.y_u), af.y_l[lower_start:end])
    return x, y
end

function write_dat(path, ms::MeshSlice)
    open(path, "w") do io
        println(io, "# axial thickness")
        for i in axes(ms.nodes, 2)
            @printf(io, "%.16e %.16e\n", ms.nodes[1, i], ms.nodes[2, i])
        end
    end
end

function write_xy_dat(path, x, y)
    open(path, "w") do io
        println(io, "# x_over_c y_over_c")
        for i in eachindex(x)
            @printf(io, "%.16e %.16e\n", x[i], y[i])
        end
    end
end

function write_indexed_dat(path, ms::MeshSlice)
    open(path, "w") do io
        println(io, "# row axial thickness")
        for i in axes(ms.nodes, 2)
            @printf(io, "%d %.16e %.16e\n", i, ms.nodes[1, i], ms.nodes[2, i])
        end
    end
end

function read_dat(path, r)
    raw = readdlm(path, comments=true, comment_char='#')
    size(raw, 2) >= 2 || error("Expected at least two columns in $path.")
    return MeshSlice(r, Matrix{Float64}(raw[:, 1:2]'))
end

function parse_coordinate_pair(line::AbstractString; default=nothing)
    answer = strip(line)
    isempty(answer) && default !== nothing && return default
    parts = split(replace(answer, ',' => ' '))
    length(parts) == 2 || error("Enter exactly two coordinates separated by whitespace or a comma.")
    values = try
        parse.(Float64, parts)
    catch
        error("Both trailing-edge coordinates must be finite numbers.")
    end
    all(isfinite, values) || error("Both trailing-edge coordinates must be finite numbers.")
    return (values[1], values[2])
end

function parse_exclusion_rows(line::AbstractString, nrows::Int;
                              le_idx::Union{Nothing,Int}=nothing,
                              default=Int[], min_retained=4 * CST_NCOEFFS)
    answer = strip(line)
    rows = if isempty(answer)
        copy(default)
    else
        parsed = Int[]
        for token in split(answer, ',')
            token = strip(token)
            isempty(token) && error("Malformed exclusion list.")
            if occursin('-', token)
                parts = split(token, '-'; keepempty=true)
                length(parts) == 2 && all(!isempty, parts) ||
                    error("Malformed exclusion range '$token'.")
                first_row = tryparse(Int, parts[1])
                last_row = tryparse(Int, parts[2])
                first_row !== nothing && last_row !== nothing ||
                    error("Malformed exclusion range '$token'.")
                first_row <= last_row || error("Descending exclusion ranges are not allowed.")
                append!(parsed, first_row:last_row)
            else
                row = tryparse(Int, token)
                row !== nothing || error("Malformed exclusion row '$token'.")
                push!(parsed, row)
            end
        end
        sort!(unique!(parsed))
    end
    all(1 .<= rows .<= nrows) || error("Exclusion rows must lie in 1:$nrows.")
    le_idx !== nothing && le_idx in rows && error("The selected leading-edge row cannot be excluded.")
    nrows - length(rows) >= min_retained ||
        error("Exclusions leave fewer than $min_retained retained scan rows.")
    return rows
end

format_coordinate(p) = @sprintf("%.10g, %.10g", p[1], p[2])
format_rows(rows) = isempty(rows) ? "none" : join(rows, ",")

function prompt_slice_row(label, nrows; default=nothing)
    while true
        suffix = default === nothing ? "" : " [$(default)]"
        print("$(label) row index (1-$(nrows))$(suffix): ")
        answer = strip(readline(stdin))
        isempty(answer) && default !== nothing && return default
        try
            idx = parse(Int, answer)
            1 <= idx <= nrows && return idx
        catch
        end
        println("Please enter an integer row index from 1 to $(nrows).")
    end
end

function prompt_skip_section(index; input=stdin, output=stdout)
    while true
        print(output, "Skip section $(index)? [y/N]: ")
        answer = lowercase(strip(readline(input)))
        answer in ("y", "yes") && return true
        (isempty(answer) || answer in ("n", "no")) && return false
        println(output, "Please answer y or n.")
    end
end

function prompt_exclusion_rows(nrows, le_idx; default=Int[])
    while true
        suffix = isempty(default) ? " [blank = none]" : " [$(format_rows(default))]"
        print("Rounded trailing-edge rows to exclude$(suffix): ")
        try
            return parse_exclusion_rows(readline(stdin), nrows; le_idx, default)
        catch err
            println(sprint(showerror, err))
        end
    end
end

function prompt_coordinate(label; default=nothing)
    while true
        suffix = default === nothing ? "" : " [$(format_coordinate(default))]"
        print("$(label) raw (x,z)$(suffix): ")
        try
            return parse_coordinate_pair(readline(stdin); default)
        catch err
            println(sprint(showerror, err))
        end
    end
end

function inverse_transform_point(x, y, origin, twist, chord)
    c = cos(twist)
    s = sin(twist)
    return (origin[1] + chord * (c * x - s * y),
            origin[2] + chord * (s * x + c * y))
end

function suggest_te_anchors(ms::MeshSlice, le_idx::Int, excluded_indices::Vector{Int};
                            near_te_fraction=0.95)
    n = size(ms.nodes, 2)
    retained_indices = setdiff(collect(1:n), excluded_indices)
    le_point = ms.nodes[:, le_idx]
    distances = [norm(ms.nodes[:, i] - le_point) for i in retained_indices]
    farthest_idx = retained_indices[argmax(distances)]
    chord_vec = ms.nodes[:, farthest_idx] - le_point
    chord = norm(chord_vec)
    chord > 0.0 || return nothing
    twist = atan(chord_vec[2], chord_vec[1])
    retained_points = ms.nodes[:, retained_indices]
    xr, yr = transform_points(retained_points, le_point, twist, chord)
    x_all, _ = transform_points(ms.nodes, le_point, twist, chord)
    xstar = maximum(x_all)
    near = findall(xr .>= near_te_fraction)
    upper = near[yr[near] .> 0.0]
    lower = near[yr[near] .< 0.0]
    if length(upper) < 2 || length(lower) < 2
        return (; upper=nothing, lower=nothing, xstar, origin=copy(le_point), twist,
                chord, reason="fewer than two retained near-TE points on each surface")
    end
    line_value(indices) = begin
        A = hcat(ones(length(indices)), xr[indices])
        coefficients = A \ yr[indices]
        coefficients[1] + coefficients[2] * xstar
    end
    upper_raw = inverse_transform_point(xstar, line_value(upper), le_point, twist, chord)
    lower_raw = inverse_transform_point(xstar, line_value(lower), le_point, twist, chord)
    return (; upper=upper_raw, lower=lower_raw, xstar, origin=copy(le_point), twist,
            chord, reason=nothing)
end

function plot_le_row_selection(ms::MeshSlice, index::Int)
    fig, ax = PythonPlot.subplots()
    ax.plot(ms.nodes[1, :], ms.nodes[2, :], "o"; color="tab:blue", markersize=3)
    for i in axes(ms.nodes, 2)
        ax.text(ms.nodes[1, i], ms.nodes[2, i], string(i); fontsize=5)
    end
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("axial x")
    ax.set_ylabel("thickness z")
    ax.set_title(@sprintf("Section %02d: choose leading-edge row, r = %.6e", index, ms.r))
    ax.grid(true)
    fig.tight_layout()
    return fig
end

function show_interactive_figures(figures...; pause_seconds=0.2)
    for figure in figures
        figure.canvas.draw()
        figure.show()
    end
    plt.show(; block=false)
    plt.pause(pause_seconds)
    return nothing
end

function plot_raw_selection(ms, index, le_idx, excluded, suggestion)
    fig, ax = PythonPlot.subplots()
    retained = setdiff(collect(axes(ms.nodes, 2)), excluded)
    ax.plot(ms.nodes[1, :], ms.nodes[2, :], "o"; color="0.75", markersize=3,
            label="all scan rows")
    ax.plot(ms.nodes[1, retained], ms.nodes[2, retained], "o"; markersize=3,
            label="retained")
    isempty(excluded) || ax.plot(ms.nodes[1, excluded], ms.nodes[2, excluded], "o";
                                  color="tab:red", markersize=4, label="excluded rounded TE")
    for i in axes(ms.nodes, 2)
        ax.text(ms.nodes[1, i], ms.nodes[2, i], string(i); fontsize=5)
    end
    ax.plot(ms.nodes[1, le_idx], ms.nodes[2, le_idx], "*"; markersize=10,
            label="selected LE")
    if suggestion !== nothing
        p0 = inverse_transform_point(suggestion.xstar, -0.1, suggestion.origin,
                                     suggestion.twist, suggestion.chord)
        p1 = inverse_transform_point(suggestion.xstar, 0.1, suggestion.origin,
                                     suggestion.twist, suggestion.chord)
        ax.plot([p0[1], p1[1]], [p0[2], p1[2]], "--"; color="tab:purple",
                label="suggested TE plane")
        suggestion.upper === nothing || ax.plot(
            [suggestion.upper[1], suggestion.lower[1]],
            [suggestion.upper[2], suggestion.lower[2]], "s-"; color="tab:green",
            label="suggested square TE")
    end
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("axial x")
    ax.set_ylabel("thickness z")
    ax.set_title(@sprintf("Section %02d anchor suggestions, r = %.6e", index, ms.r))
    ax.grid(true)
    ax.legend()
    fig.tight_layout()
    return fig
end

function select_airfoil_anchors(ms::MeshSlice, index::Int; previous=nothing,
                                review_paths=section_review_paths(index))
    mkpath(review_paths.directory)
    nrows = size(ms.nodes, 2)
    le_default = previous === nothing ? nothing : previous.le_idx
    exclusion_default = previous === nothing ? Int[] : previous.excluded_indices
    le_figure = plot_le_row_selection(ms, index)
    le_figure.savefig(review_paths.leading_edge; dpi=180)
    println("Leading-edge row guide: $(review_paths.leading_edge)")
    show_interactive_figures(le_figure)
    if previous === nothing && prompt_skip_section(index)
        plt.close(le_figure)
        println("Skipping section $index.")
        return nothing
    end
    le_idx = prompt_slice_row("Leading edge", nrows; default=le_default)
    excluded = prompt_exclusion_rows(nrows, le_idx; default=exclusion_default)
    plt.close(le_figure)
    suggestion = suggest_te_anchors(ms, le_idx, excluded)
    if suggestion === nothing || suggestion.upper === nothing
        reason = suggestion === nothing ? "degenerate provisional chord" : suggestion.reason
        println("No automatic square-TE suggestion: $reason. Enter both coordinates manually.")
    end
    anchor_figure = plot_raw_selection(ms, index, le_idx, excluded, suggestion)
    show_interactive_figures(anchor_figure)
    upper_default = previous === nothing ?
        (suggestion === nothing ? nothing : suggestion.upper) : previous.upper_te
    lower_default = previous === nothing ?
        (suggestion === nothing ? nothing : suggestion.lower) : previous.lower_te
    blend_start = previous === nothing ? SHARP_TE_BLEND_START : previous.blend_start
    while true
        upper_te = prompt_coordinate("Upper square trailing edge"; default=upper_default)
        lower_te = prompt_coordinate("Lower square trailing edge"; default=lower_default)
        anchored = AnchoredMeshSlice(ms.r, ms.nodes, upper_te, le_idx, lower_te,
                                     excluded, blend_start)
        try
            normalized_section(anchored)
            plt.close(anchor_figure)
            return anchored, suggestion
        catch err
            println("Invalid trailing-edge selection: ", sprint(showerror, err))
            upper_default = upper_te
            lower_default = lower_te
        end
    end
end

function plot_selection_checkpoint(ms::AnchoredMeshSlice, normalized, suggestion, index)
    raw_fig, raw = PythonPlot.subplots(; figsize=(7, 5))
    normalization_fig, normax = PythonPlot.subplots(; figsize=(7, 5))
    retained_raw = setdiff(collect(axes(ms.nodes, 2)), ms.excluded_indices)
    raw.plot(ms.nodes[1, :], ms.nodes[2, :], "o"; color="0.8", markersize=3,
             label="all rows")
    raw.plot(ms.nodes[1, retained_raw], ms.nodes[2, retained_raw], "o"; markersize=3,
             label="retained")
    isempty(ms.excluded_indices) || raw.plot(
        ms.nodes[1, ms.excluded_indices], ms.nodes[2, ms.excluded_indices], "o";
        color="tab:red", label="excluded rounded TE")
    for i in axes(ms.nodes, 2)
        raw.text(ms.nodes[1, i], ms.nodes[2, i], string(i); fontsize=5)
    end
    upper = collect(ms.upper_te); lower = collect(ms.lower_te)
    raw.plot([upper[1], lower[1]], [upper[2], lower[2]], "s-"; color="tab:green",
             label="accepted square TE")
    raw.plot(normalized.te_point[1], normalized.te_point[2], "D"; label="TE midpoint")
    raw.plot([normalized.le_point[1], normalized.te_point[1]],
             [normalized.le_point[2], normalized.te_point[2]], "k--"; label="chord")
    if suggestion !== nothing
        p0 = inverse_transform_point(suggestion.xstar, -0.1, suggestion.origin,
                                     suggestion.twist, suggestion.chord)
        p1 = inverse_transform_point(suggestion.xstar, 0.1, suggestion.origin,
                                     suggestion.twist, suggestion.chord)
        raw.plot([p0[1], p1[1]], [p0[2], p1[2]], ":"; color="tab:purple")
    end

    normax.plot(normalized.x, normalized.y, "o"; markersize=3, label="retained normalized")
    isempty(normalized.dropped_indices) || normax.plot(
        normalized.x_all[normalized.dropped_indices], normalized.y_all[normalized.dropped_indices],
        "x"; color="0.65", label="dropped outside [0,1]")
    isempty(normalized.likely_missed) || normax.plot(
        normalized.x_all[normalized.likely_missed], normalized.y_all[normalized.likely_missed],
        "o"; markerfacecolor="none", markeredgecolor="tab:red", markersize=7,
        label="possible missed rounded TE")
    normax.plot(normalized.anchor_x, normalized.anchor_y, "s"; label="transformed anchors")
    normax.plot(normalized.canonical_x, normalized.canonical_y, "D";
                label="canonical anchors")
    for i in 1:2
        normax.plot([normalized.anchor_x[i], normalized.canonical_x[i]],
                    [normalized.anchor_y[i], normalized.canonical_y[i]], "--"; color="0.4")
    end
    normax.plot([0.0, 1.0], [0.0, 0.0], "k--"; linewidth=1, label="chord")
    annotation = @sprintf("TE thickness %.5g\nΔx upper/lower %.4g, %.4g\ncorrection %.4g, %.4g",
                          2 * normalized.h_te, normalized.x_offsets...,
                          normalized.displacements...)
    normax.text(0.02, 0.98, annotation; transform=normax.transAxes, va="top")
    raw.set_aspect("equal", adjustable="box")
    raw.grid(true)
    raw.legend(fontsize=7)
    raw.set_title("Raw selection")
    raw_fig.tight_layout()

    normax.set_aspect("equal", adjustable="box")
    normax.grid(true)
    normax.legend(fontsize=7)
    normax.set_title("Normalization and canonicalization")
    normalization_fig.tight_layout()
    return raw_fig, normalization_fig
end

function count_upstream_of_le(ms::AnchoredMeshSlice)
    normalized = normalized_section(ms)
    retained = setdiff(collect(axes(ms.nodes, 2)), ms.excluded_indices)
    return count(normalized.x_all[retained] .< 0.0)
end

function approve_selection(ms::MeshSlice, index::Int; previous=nothing,
                           review_paths=section_review_paths(index))
    mkpath(review_paths.directory)
    current = previous
    while true
        selection = select_airfoil_anchors(ms, index; previous=current, review_paths)
        selection === nothing && return nothing
        anchored, suggestion = selection
        normalized = normalized_section(anchored)
        if any(abs.(normalized.x_offsets) .> 0.02)
            println("Warning: at least one transformed TE anchor differs from x/c=1 by more than 0.02.")
        end
        n_upstream = count_upstream_of_le(anchored)
        if n_upstream > 0
            println("$n_upstream retained point(s) lie upstream of the selected leading edge.")
            while true
                print("[r]emove those points via the chordwise filter, or [c]hoose again? ")
                answer = lowercase(strip(readline(stdin)))
                answer in ("r", "remove") && break
                if answer in ("c", "choose")
                    current = anchored
                    anchored = nothing
                    break
                end
                println("Please answer r or c.")
            end
            anchored === nothing && continue
        end
        raw_fig, normalization_fig = plot_selection_checkpoint(
            anchored, normalized, suggestion, index)
        show_interactive_figures(raw_fig, normalization_fig)
        while true
            print("Approve selection and canonicalization for section $index? [y/n]: ")
            answer = lowercase(strip(readline(stdin)))
            if answer in ("y", "yes")
                raw_fig.savefig(review_paths.raw_selection; dpi=180)
                normalization_fig.savefig(review_paths.normalization; dpi=180)
                plt.close(raw_fig)
                plt.close(normalization_fig)
                return anchored, normalized
            elseif answer in ("n", "no")
                plt.close(raw_fig)
                plt.close(normalization_fig)
                current = anchored
                break
            else
                println("Please answer y or n.")
            end
        end
    end
end

bernstein(r, n, x) = binomial(n, r) .* x .^ r .* (1 .- x) .^ (n .- r)

function half_cst(coefficients, x, dz, leading_edge_weight; N1=0.5, N2=1.0)
    nb = length(coefficients) - 1
    C = @. x^N1 * (1.0 - x)^N2
    S = zero.(x)
    for (i, coefficient) in enumerate(coefficients)
        S .+= coefficient .* bernstein(i - 1, nb, x)
    end
    y = @. C * S + x * dz
    y .+= leading_edge_weight .* x .* (1.0 .- x) .^ (length(coefficients) + 0.5)
    return y
end

function cst_coordinates(p; N::Int=CST_RECONSTRUCTION_POINTS)
    np = Int((length(p) - 2) / 2)
    pu = p[1:np]
    pl = p[(np + 1):(2 * np)]
    leading_edge_weight = p[end - 1]
    dz = p[end]
    θ = range(0.0, π; length=N)
    xu = @. 0.5 * (1.0 - cos(θ))
    xl = copy(xu)
    yu = half_cst(pu, xu, dz / 2.0, leading_edge_weight)
    yl = half_cst(pl, xl, -dz / 2.0, leading_edge_weight)
    return vcat(reverse(xu), xl[2:end]), vcat(reverse(yu), yl[2:end])
end

function fit_cst_half(x, y, n_coeffs; dz_half=0.0)
    nb = n_coeffs - 1
    C = @. x^0.5 * (1.0 - x)^1.0
    A = zeros(length(x), n_coeffs)
    for j in 0:nb
        A[:, j + 1] .= C .* bernstein(j, nb, x)
    end
    rhs = y .- x .* dz_half
    mask = @. x > 1e-10 && x < 1.0 - 1e-10
    return A[mask, :] \ rhs[mask]
end

function fit_cst(xu, yu, xl, yl; n_coeffs=CST_NCOEFFS)
    dz = yu[end] - yl[end]
    w_upper = fit_cst_half(xu, yu, n_coeffs; dz_half=dz / 2.0)
    w_lower = fit_cst_half(xl, yl, n_coeffs; dz_half=-dz / 2.0)
    return [w_upper; w_lower; 0.0; dz]
end

function cosine_grid(N::Int=CST_RECONSTRUCTION_POINTS)
    N >= 2 || error("Cosine grid requires at least two points.")
    θ = range(0.0, π; length=N)
    x = @. 0.5 * (1.0 - cos(θ))
    x[1] = 0.0
    x[end] = 1.0
    return x
end

function cst_halves(params, x)
    np = Int((length(params) - 2) / 2)
    leading_edge_weight = params[end - 1]
    dz = params[end]
    yu = half_cst(params[1:np], x, dz / 2, leading_edge_weight)
    yl = half_cst(params[(np + 1):(2np)], x, -dz / 2, leading_edge_weight)
    return yu, yl
end

function fit_blunt_cst(ms::AnchoredMeshSlice, normalized=normalized_section(ms);
                       n_coeffs=CST_NCOEFFS, N=CST_RECONSTRUCTION_POINTS)
    samples = ordered_surface_samples(ms, normalized)
    n_upper_interior = count(@. samples.xu > 1e-10 && samples.xu < 1 - 1e-10)
    n_lower_interior = count(@. samples.xl > 1e-10 && samples.xl < 1 - 1e-10)
    minimum((n_upper_interior, n_lower_interior)) >= n_coeffs + 2 || error(
        "Each surface needs at least $(n_coeffs + 2) interior points after binning; " *
        "found $n_upper_interior upper and $n_lower_interior lower.")
    params = fit_cst(samples.xu, samples.yu, samples.xl, samples.yl; n_coeffs)
    params[end] = 2 * normalized.h_te
    all(isfinite, params) || error("Blunt CST fit produced non-finite coefficients.")
    x_common = cosine_grid(N)
    yu, yl = cst_halves(params, x_common)
    yu[end] = normalized.h_te
    yl[end] = -normalized.h_te
    camber = 0.5 .* (yu .+ yl)
    half_thickness = 0.5 .* (yu .- yl)
    isapprox(camber[end], 0.0; atol=100eps(Float64), rtol=0.0) ||
        error("Blunt CST mean camber does not end at zero.")
    isapprox(half_thickness[end], normalized.h_te; atol=100eps(Float64), rtol=1e-12) ||
        error("Blunt CST did not preserve the canonical trailing-edge half-thickness.")
    thickness_tol = 1000eps(Float64) * max(1.0, maximum(abs, half_thickness))
    minimum(half_thickness) >= -thickness_tol ||
        error("Blunt CST upper and lower surfaces cross.")
    half_thickness .= max.(half_thickness, 0.0)
    return (; samples, params, x_common, yu, yl, camber, half_thickness)
end

quadratic_thickness_weight(s) = 1 - s^2
quadratic_thickness_weight_derivative(s) = -2s

function thickness_multiplier(x, blend_start::Real)
    0.0 <= blend_start < 1.0 || error("`blend_start` must be in [0, 1).")
    return map(x) do xi
        xi <= blend_start && return 1.0
        s = (xi - blend_start) / (1 - blend_start)
        quadratic_thickness_weight(s)
    end
end

function taper_blunt_geometry(blunt, blend_start::Real)
    w = thickness_multiplier(blunt.x_common, blend_start)
    yu = blunt.camber .+ w .* blunt.half_thickness
    yl = blunt.camber .- w .* blunt.half_thickness
    forward = blunt.x_common .<= blend_start
    yu[forward] .= blunt.yu[forward]
    yl[forward] .= blunt.yl[forward]
    yu[end] = 0.0
    yl[end] = 0.0
    return (; x=blunt.x_common, yu, yl, w,
            camber=0.5 .* (yu .+ yl), half_thickness=0.5 .* (yu .- yl))
end

function fit_final_cst(sharp; n_coeffs=CST_NCOEFFS, N=CST_RECONSTRUCTION_POINTS)
    params = fit_cst(sharp.x, sharp.yu, sharp.x, sharp.yl; n_coeffs)
    params[end] == 0.0 || error("Final sharp CST must have exactly zero trailing-edge thickness.")
    all(isfinite, params) || error("Final CST fit produced non-finite coefficients.")
    x = cosine_grid(N)
    yu, yl = cst_halves(params, x)
    yu[end] = 0.0
    yl[end] = 0.0
    return (; params, x, yu, yl, camber=0.5 .* (yu .+ yl))
end

function Airfoil(ms::AnchoredMeshSlice)
    normalized = normalized_section(ms)
    blunt = fit_blunt_cst(ms, normalized)
    sharp = taper_blunt_geometry(blunt, ms.blend_start)
    chord_center = 0.5 .* (normalized.le_point .+ normalized.te_point)
    return Airfoil(ms.r, normalized.twist, normalized.chord,
                   -chord_center[2], -chord_center[1],
                   copy(sharp.x), copy(sharp.yu), copy(sharp.x), copy(sharp.yl))
end

residual_stats(residual) = (; rms=sqrt(mean(abs2, residual)), maxabs=maximum(abs, residual))

function reconstruction_diagnostics(blunt, sharp, final)
    upper_blunt_fit = cst_halves(blunt.params, blunt.samples.xu)[1]
    lower_blunt_fit = cst_halves(blunt.params, blunt.samples.xl)[2]
    blunt_upper_residual = blunt.samples.yu .- upper_blunt_fit
    blunt_lower_residual = blunt.samples.yl .- lower_blunt_fit
    final_upper_residual = sharp.yu .- cst_halves(final.params, sharp.x)[1]
    final_lower_residual = sharp.yl .- cst_halves(final.params, sharp.x)[2]
    camber_error = final.camber .- sharp.camber
    return (; blunt_upper_residual, blunt_lower_residual,
            blunt_upper_stats=residual_stats(blunt_upper_residual),
            blunt_lower_stats=residual_stats(blunt_lower_residual),
            final_upper_residual, final_lower_residual,
            final_upper_stats=residual_stats(final_upper_residual),
            final_lower_stats=residual_stats(final_lower_residual), camber_error)
end

function coords_to_cst(xaf, yaf; n_coeffs=CST_NCOEFFS)
    le_idx = argmin(xaf)
    xu = reverse(xaf[1:le_idx])
    yu = reverse(yaf[1:le_idx])
    xl = xaf[le_idx:end]
    yl = yaf[le_idx:end]
    return fit_cst(xu, yu, xl, yl; n_coeffs)
end

function airfoil_to_cst(af::Airfoil; n_coeffs=CST_NCOEFFS)
    x, y = selig_coordinates(af)
    return coords_to_cst(x, y; n_coeffs)
end

function plot_blunt_diagnostics(ms, normalized, blunt, diagnostics, index)
    fig, axs = PythonPlot.subplots(1, 2; figsize=(13, 5))
    ax, residual_ax = axs[0], axs[1]
    samples = blunt.samples
    ax.plot(samples.xu, samples.yu, "o"; markersize=3, label="upper scan samples")
    ax.plot(samples.xl, samples.yl, "o"; markersize=3, label="lower scan samples")
    ax.plot(blunt.x_common, blunt.yu, "-"; label="upper blunt CST")
    ax.plot(blunt.x_common, blunt.yl, "-"; label="lower blunt CST")
    ax.plot(normalized.canonical_x, normalized.canonical_y, "D";
            label="canonical TE anchors")
    isempty(ms.excluded_indices) || ax.plot(
        normalized.x_all[ms.excluded_indices], normalized.y_all[ms.excluded_indices], "x";
        color="0.65", label="excluded raw rows")
    ax.set_xlim(0.7, 1.02)
    local_y = vcat(blunt.yu[blunt.x_common .>= 0.7], blunt.yl[blunt.x_common .>= 0.7])
    margin = 0.1 * max(1e-6, maximum(local_y) - minimum(local_y))
    ax.set_ylim(minimum(local_y) - margin, maximum(local_y) + margin)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Blunt CST geometry (TE close-up)")
    ax.grid(true); ax.legend(fontsize=7)

    residual_ax.plot(samples.xu, diagnostics.blunt_upper_residual, "o-";
                     markersize=3, label="upper residual")
    residual_ax.plot(samples.xl, diagnostics.blunt_lower_residual, "o-";
                     markersize=3, label="lower residual")
    residual_ax.axhline(0.0; color="k", linewidth=1)
    text = @sprintf("upper RMS/max %.3g / %.3g\nlower RMS/max %.3g / %.3g",
                    diagnostics.blunt_upper_stats.rms, diagnostics.blunt_upper_stats.maxabs,
                    diagnostics.blunt_lower_stats.rms, diagnostics.blunt_lower_stats.maxabs)
    residual_ax.text(0.02, 0.98, text; transform=residual_ax.transAxes, va="top")
    residual_ax.set_title("Blunt-fit residuals")
    residual_ax.grid(true); residual_ax.legend(fontsize=7)
    fig.suptitle(@sprintf("Section %02d blunt CST, r = %.6e", index, ms.r))
    fig.tight_layout()
    return fig
end

function plot_taper_diagnostics(blunt, sharp, blend_start, index)
    fig, axs = PythonPlot.subplots(1, 2; figsize=(13, 5))
    ax, dist_ax = axs[0], axs[1]
    ax.plot(blunt.x_common, blunt.yu, "--"; label="blunt upper")
    ax.plot(blunt.x_common, blunt.yl, "--"; label="blunt lower")
    ax.plot(sharp.x, sharp.yu, "-"; label="sharp target upper")
    ax.plot(sharp.x, sharp.yl, "-"; label="sharp target lower")
    ax.plot(blunt.x_common, blunt.camber, "k-"; linewidth=2, label="blunt camber")
    ax.plot(sharp.x, sharp.camber, ":"; color="tab:orange", linewidth=2,
            label="sharp camber")
    ax.axvline(blend_start; color="0.4", linestyle="--")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Camber-preserving thickness taper")
    ax.grid(true); ax.legend(fontsize=7)

    camber_difference = sharp.camber .- blunt.camber
    dist_ax.plot(sharp.x, 2 .* blunt.half_thickness, label="blunt thickness")
    dist_ax.plot(sharp.x, 2 .* sharp.half_thickness, label="tapered thickness")
    dist_ax.plot(sharp.x, camber_difference, label="camber change")
    dist_ax.axhline(0.0; color="k", linewidth=1)
    forward = sharp.x .<= blend_start
    geometry_change = maximum(vcat(abs.(sharp.yu[forward] .- blunt.yu[forward]),
                                   abs.(sharp.yl[forward] .- blunt.yl[forward])))
    text = @sprintf("max |camber change| %.3g\nfinal TE thickness %.3g\nmax forward change %.3g",
                    maximum(abs, camber_difference), 2 * sharp.half_thickness[end], geometry_change)
    dist_ax.text(0.02, 0.98, text; transform=dist_ax.transAxes, va="top")
    dist_ax.set_title("Thickness and camber distributions")
    dist_ax.grid(true); dist_ax.legend(fontsize=7)
    fig.suptitle(@sprintf("Section %02d taper from x/c = %.4g", index, blend_start))
    fig.tight_layout()
    return fig
end

function plot_final_diagnostics(sharp, final, diagnostics, index)
    fig, axs = PythonPlot.subplots(2, 2; figsize=(13, 9))
    geometry, camber_ax, residual_ax, te_ax = axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]
    geometry.plot(sharp.x, sharp.yu, "--"; label="sharp upper target")
    geometry.plot(sharp.x, sharp.yl, "--"; label="sharp lower target")
    geometry.plot(final.x, final.yu, "-"; label="final upper CST")
    geometry.plot(final.x, final.yl, "-"; label="final lower CST")
    geometry.set_aspect("equal", adjustable="box")
    geometry.set_title("Whole-airfoil final fit")
    geometry.grid(true); geometry.legend(fontsize=7)

    camber_ax.plot(sharp.x, sharp.camber, label="sharp target camber")
    camber_ax.plot(final.x, final.camber, label="final CST camber")
    camber_ax.plot(sharp.x, diagnostics.camber_error, label="camber error")
    camber_ax.axhline(0.0; color="k", linewidth=1)
    camber_ax.set_title("Mean camber")
    camber_ax.grid(true); camber_ax.legend(fontsize=7)

    residual_ax.plot(sharp.x, diagnostics.final_upper_residual, label="upper residual")
    residual_ax.plot(sharp.x, diagnostics.final_lower_residual, label="lower residual")
    residual_ax.axhline(0.0; color="k", linewidth=1)
    text = @sprintf("upper RMS/max %.3g / %.3g\nlower RMS/max %.3g / %.3g",
                    diagnostics.final_upper_stats.rms, diagnostics.final_upper_stats.maxabs,
                    diagnostics.final_lower_stats.rms, diagnostics.final_lower_stats.maxabs)
    residual_ax.text(0.02, 0.98, text; transform=residual_ax.transAxes, va="top")
    residual_ax.set_title("Target-minus-fit residuals")
    residual_ax.grid(true); residual_ax.legend(fontsize=7)

    te_mask = sharp.x .>= 0.9
    te_ax.plot(sharp.x[te_mask], sharp.yu[te_mask], "--"; label="upper target")
    te_ax.plot(sharp.x[te_mask], sharp.yl[te_mask], "--"; label="lower target")
    te_ax.plot(final.x[te_mask], final.yu[te_mask], "-"; label="upper CST")
    te_ax.plot(final.x[te_mask], final.yl[te_mask], "-"; label="lower CST")
    te_ax.plot(1.0, 0.0, "D"; label="coincident TE")
    te_ax.set_title("Trailing-edge close-up")
    te_ax.grid(true); te_ax.legend(fontsize=7)
    fig.suptitle(@sprintf("Section %02d final sharp CST", index))
    fig.tight_layout()
    return fig
end

function prompt_blend_start(current)
    while true
        print("New blend_start [$(current)]: ")
        answer = strip(readline(stdin))
        isempty(answer) && return current
        value = tryparse(Float64, answer)
        value !== nothing && 0.0 <= value < 1.0 && return value
        println("Please enter a finite value in [0, 1).")
    end
end

function approve_processed_airfoil(ms::MeshSlice, index::Int)
    review_paths = section_review_paths(index)
    mkpath(review_paths.directory)
    write_indexed_dat(review_paths.raw, ms)
    println("Wrote raw indexed section coordinates to $(review_paths.raw)")
    previous = nothing
    while true
        selection = approve_selection(ms, index; previous, review_paths)
        selection === nothing && return nothing
        anchored, normalized = selection
        blunt = fit_blunt_cst(anchored, normalized)
        blend_start = anchored.blend_start
        while true
            sharp = taper_blunt_geometry(blunt, blend_start)
            final = fit_final_cst(sharp)
            diagnostics = reconstruction_diagnostics(blunt, sharp, final)
            fig_a = plot_blunt_diagnostics(anchored, normalized, blunt, diagnostics, index)
            fig_b = plot_taper_diagnostics(blunt, sharp, blend_start, index)
            fig_c = plot_final_diagnostics(sharp, final, diagnostics, index)
            fig_a.savefig(review_paths.blunt; dpi=180)
            fig_b.savefig(review_paths.taper; dpi=180)
            fig_c.savefig(review_paths.final; dpi=180)
            show_interactive_figures(fig_a, fig_b, fig_c)
            print("Approve reconstruction for section $index? [y/n/b = change blend start]: ")
            answer = lowercase(strip(readline(stdin)))
            for fig in (fig_a, fig_b, fig_c)
                plt.close(fig)
            end
            if answer in ("y", "yes")
                chord_center = 0.5 .* (normalized.le_point .+ normalized.te_point)
                af = Airfoil(anchored.r, normalized.twist, normalized.chord,
                             -chord_center[2], -chord_center[1],
                             copy(sharp.x), copy(sharp.yu), copy(sharp.x), copy(sharp.yl))
                xselig, yselig = selig_coordinates(af)
                write_xy_dat(review_paths.airfoil, xselig, yselig)
                return af, final.params
            elseif answer in ("b", "blend")
                blend_start = prompt_blend_start(blend_start)
            elseif answer in ("n", "no")
                previous = AnchoredMeshSlice(anchored.r, anchored.nodes, anchored.upper_te,
                    anchored.le_idx, anchored.lower_te, anchored.excluded_indices, blend_start)
                break
            else
                println("Please answer y, n, or b.")
            end
        end
    end
end

function record_processed_section!(airfoils, all_params, result)
    result === nothing && return false
    af, params = result
    push!(airfoils, af)
    push!(all_params, params)
    return true
end

function write_cst_table(path, airfoils, all_params)
    header = ["r", "twist_rad", "chord", "axial", "tangential"]
    append!(header, ["wu$(i)" for i in 1:CST_NCOEFFS])
    append!(header, ["wl$(i)" for i in 1:CST_NCOEFFS])
    append!(header, ["leading_edge_weight", "dz"])

    open(path, "w") do io
        println(io, join(header, ","))
        for (af, params) in zip(airfoils, all_params)
            row = vcat([af.r, af.twist, af.chord, af.axial, af.tangential], params)
            println(io, join((@sprintf("%.16e", value) for value in row), ","))
        end
    end
end

function main()
    mkpath(OUTPUT_DIR)
    nodes, cells, _ = load_body(INPUT_MESH)
    raw_stem = joinpath(OUTPUT_DIR, RUN_NAME * "_raw")
    println("Writing raw mesh preview: $(raw_stem).pvd")
    write_body_vtk(raw_stem, nodes, cells)

    require_indices(nodes)

    aligned_nodes, origin, basis = alignment_from_manual_points(nodes)
    aligned_stem = joinpath(OUTPUT_DIR, RUN_NAME * "_aligned")
    println("Writing aligned mesh preview: $(aligned_stem).pvd")
    write_body_vtk(aligned_stem, aligned_nodes, cells)
    println("Alignment origin: ", origin)
    println("Alignment basis columns are x/y/z axes in raw coordinates:")
    display(basis)

    slices = generate_slices(aligned_nodes)
    airfoils = Airfoil[]
    all_params = Vector{Float64}[]
    for (i, slice) in enumerate(slices)
        result = approve_processed_airfoil(slice, i)
        record_processed_section!(airfoils, all_params, result)
    end

    cst_path = cst_output_path()
    write_cst_table(cst_path, airfoils, all_params)
    println("Wrote CST table: $cst_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
