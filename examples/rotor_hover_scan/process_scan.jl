import FLOWPanel as pnl
import GeoIO
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

const NSLICES = 1
const SLICE_MIN_POINTS = 400
const SLICE_DISTANCE_ATOL = 0.0
const N_RFL_SECS = 100
const USE_COSINE_RFL_BINS = true
const HUB_TOP_FRACTION = 0.75
const CST_NCOEFFS = 8
const CST_RECONSTRUCTION_POINTS = 201
const SHARP_TE_BLEND_START = 0.8

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
    upper_te_idx::Int
    le_idx::Int
    lower_te_idx::Int
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
    msh = GeoIO.load(mesh_path).geometry
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

function rotated_section(ms::AnchoredMeshSlice)
    te_point = 0.5 .* (ms.nodes[:, ms.upper_te_idx] .+ ms.nodes[:, ms.lower_te_idx])
    le_point = ms.nodes[:, ms.le_idx]
    chord_vec = te_point - le_point
    chord = norm(chord_vec)
    chord > 0.0 || error("Slice at r=$(ms.r) has zero chord.")

    twist = atan(chord_vec[2], chord_vec[1])
    c = cos(-twist)
    s = sin(-twist)
    x_shift = ms.nodes[1, :] .- le_point[1]
    y_shift = ms.nodes[2, :] .- le_point[2]
    x = (c .* x_shift .- s .* y_shift) ./ chord
    y = (s .* x_shift .+ c .* y_shift) ./ chord

    x_le = x[ms.le_idx]
    x_te = 0.5 * (x[ms.upper_te_idx] + x[ms.lower_te_idx])
    x[ms.upper_te_idx] = x_te
    x[ms.lower_te_idx] = x_te
    chordwise_mask = (x .>= x_le) .& (x .<= x_te)
    any(chordwise_mask) || error("Slice at r=$(ms.r) has no points between the selected LE and TE bounds.")

    return twist, chord, x[chordwise_mask], y[chordwise_mask], findall(chordwise_mask)
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

function required_selig_indices(original_indices, ms::AnchoredMeshSlice)
    upper_te = findfirst(==(ms.upper_te_idx), original_indices)
    le = findfirst(==(ms.le_idx), original_indices)
    lower_te = findfirst(==(ms.lower_te_idx), original_indices)
    upper_te !== nothing || error("Selected upper TE row was discarded by the chordwise filter.")
    le !== nothing || error("Selected LE row was discarded by the chordwise filter.")
    lower_te !== nothing || error("Selected lower TE row was discarded by the chordwise filter.")
    return [upper_te, le], [le, lower_te], upper_te, lower_te
end

function selig_ordered_slice(ms::AnchoredMeshSlice)
    _, _, x, y, original_indices = rotated_section(ms)
    required_upper, required_lower, upper_te_pos, lower_te_pos =
        required_selig_indices(original_indices, ms)
    selig_indices = sectioned_selig_indices(x, y;
                                            required_upper_indices=required_upper,
                                            required_lower_indices=required_lower,
                                            upper_te_pos=upper_te_pos,
                                            lower_te_pos=lower_te_pos)
    return MeshSlice(ms.r, ms.nodes[:, original_indices[selig_indices]])
end

function Airfoil(ms::AnchoredMeshSlice)
    twist, chord, x, y, original_indices = rotated_section(ms)
    required_upper, required_lower, upper_te_pos, lower_te_pos =
        required_selig_indices(original_indices, ms)
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

    sharpen_trailing_edge!(xu, yu, xl, yl)

    te_point = 0.5 .* (ms.nodes[:, ms.upper_te_idx] .+ ms.nodes[:, ms.lower_te_idx])
    le_point = ms.nodes[:, ms.le_idx]
    chord_center = 0.5 .* (le_point .+ te_point)
    tangential = -chord_center[1]
    axial = -chord_center[2]

    return Airfoil(ms.r, twist, chord, axial, tangential, xu, yu, xl, yl)
end

function sharpen_trailing_edge!(xu, yu, xl, yl;
                                blend_start::Float64=SHARP_TE_BLEND_START)
    0.0 <= blend_start < 1.0 || error("`blend_start` must be in [0, 1).")
    y_te_mid = 0.5 * (yu[end] + yl[end])
    delta_upper = y_te_mid - yu[end]
    delta_lower = y_te_mid - yl[end]
    blend_range = 1.0 - blend_start
    for i in eachindex(xu)
        if xu[i] >= blend_start
            t = (xu[i] - blend_start) / blend_range
            yu[i] += t * delta_upper
        end
    end
    for i in eachindex(xl)
        if xl[i] >= blend_start
            t = (xl[i] - blend_start) / blend_range
            yl[i] += t * delta_lower
        end
    end
    return xu, yu, xl, yl
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

function prompt_slice_row(label, nrows)
    while true
        print("$(label) row index [1-$(nrows)]: ")
        answer = strip(readline(stdin))
        try
            idx = parse(Int, answer)
            1 <= idx <= nrows && return idx
        catch
        end
        println("Please enter an integer row index from 1 to $(nrows).")
    end
end

function select_airfoil_anchors(ms::MeshSlice, index::Int)
    raw_path = joinpath(OUTPUT_DIR, "airfoils", @sprintf("section_%02d_raw.dat", index))
    write_indexed_dat(raw_path, ms)

    fig, ax = PythonPlot.subplots()
    ax.plot(ms.nodes[1, :], ms.nodes[2, :], "o"; markersize=2.5)
    if size(ms.nodes, 2) <= 400
        for i in axes(ms.nodes, 2)
            ax.text(ms.nodes[1, i], ms.nodes[2, i], string(i); fontsize=6)
        end
    end
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_title(@sprintf("Raw section %02d, r = %.6e", index, ms.r))
    ax.grid(true)
    fig.tight_layout()
    PythonPlot.show()
    plt.close(fig)

    println("Wrote raw indexed section coordinates to $raw_path")
    println("Select row indices from that file for the airfoil anchors.")
    upper_te_idx = prompt_slice_row("Upper trailing edge", size(ms.nodes, 2))
    le_idx = prompt_slice_row("Leading edge", size(ms.nodes, 2))
    lower_te_idx = prompt_slice_row("Lower trailing edge", size(ms.nodes, 2))
    le_idx != upper_te_idx || error("Leading edge row must differ from upper trailing edge row.")
    le_idx != lower_te_idx || error("Leading edge row must differ from lower trailing edge row.")

    return AnchoredMeshSlice(ms.r, ms.nodes, upper_te_idx, le_idx, lower_te_idx)
end

function count_upstream_of_le(ms::AnchoredMeshSlice)
    te_point = 0.5 .* (ms.nodes[:, ms.upper_te_idx] .+ ms.nodes[:, ms.lower_te_idx])
    le_point = ms.nodes[:, ms.le_idx]
    chord_vec = te_point - le_point
    chord = norm(chord_vec)
    chord > 0.0 || error("Slice at r=$(ms.r) has zero chord.")
    twist = atan(chord_vec[2], chord_vec[1])
    c = cos(-twist)
    s = sin(-twist)
    x_shift = ms.nodes[1, :] .- le_point[1]
    y_shift = ms.nodes[2, :] .- le_point[2]
    x = (c .* x_shift .- s .* y_shift) ./ chord
    return count(x .< 0.0)
end

function select_airfoil_anchors_validated(ms::MeshSlice, index::Int)
    current = select_airfoil_anchors(ms, index)
    while true
        n_upstream = count_upstream_of_le(current)
        n_upstream == 0 && return current
        println("$(n_upstream) point(s) lie upstream of the leading edge (x < 0 after rotation).")
        println("CST fitting requires all points to have x >= 0.")
        print("[r]emove those points and continue, or [c]hoose a different leading edge? ")
        answer = lowercase(strip(readline(stdin)))
        if answer in ("r", "remove")
            return current
        elseif answer in ("c", "choose")
            current = select_airfoil_anchors(ms, index)
        else
            println("Please answer r or c.")
        end
    end
end

function approve_airfoil(ms::MeshSlice, index::Int)
    mkpath(joinpath(OUTPUT_DIR, "airfoils"))
    editable_path = joinpath(OUTPUT_DIR, "airfoils", @sprintf("section_%02d.dat", index))

    current_slice = select_airfoil_anchors_validated(ms, index)
    while true
        af = Airfoil(current_slice)
        xselig, yselig = selig_coordinates(af)
        fig, ax = PythonPlot.subplots()
        ax.plot(xselig, yselig,
                "o"; markersize=2.5, label="slice/chord")
        ax.plot(xselig, yselig, "-"; linewidth=1.5, label="airfoil")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x/c")
        ax.set_ylabel("y/c")
        ax.set_title(@sprintf("Section %02d, r = %.6e", index, current_slice.r))
        ax.grid(true)
        ax.legend()
        fig.tight_layout()
        PythonPlot.show()
        plt.close(fig)

        print("Approve section $(index)? [y/n]: ")
        answer = lowercase(strip(readline(stdin)))
        if answer in ("y", "yes")
            write_xy_dat(editable_path, xselig, yselig)
            return af
        elseif answer in ("n", "no")
            write_xy_dat(editable_path, xselig, yselig)
            println("Wrote editable coordinates to $editable_path")
            println("Edit the file, then press Enter to reload it.")
            readline(stdin)
            current_slice = select_airfoil_anchors_validated(read_dat(editable_path, current_slice.r), index)
        else
            println("Please answer y or n.")
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

function approve_cst(index, r, editable_path; n_coeffs=CST_NCOEFFS,
                     N=CST_RECONSTRUCTION_POINTS)
    while true
        raw = readdlm(editable_path, comments=true, comment_char='#')
        size(raw, 2) >= 2 || error("Expected at least two columns in $editable_path.")
        xaf = Float64.(raw[:, 1])
        yaf = Float64.(raw[:, 2])
        params = coords_to_cst(xaf, yaf; n_coeffs=n_coeffs)
        x_cst, y_cst = cst_coordinates(params; N=N)

        fig, ax = PythonPlot.subplots()
        ax.plot(xaf, yaf, "o"; markersize=2.5, label="airfoil")
        ax.plot(x_cst, y_cst, "-"; linewidth=1.5, label="CST fit")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x/c")
        ax.set_ylabel("y/c")
        ax.set_title(@sprintf("Section %02d CST fit, r = %.6e", index, r))
        ax.grid(true)
        ax.legend()
        fig.tight_layout()
        PythonPlot.show()
        plt.close(fig)

        print("Approve CST fit for section $(index)? [y/n]: ")
        answer = lowercase(strip(readline(stdin)))
        if answer in ("y", "yes")
            return params
        elseif answer in ("n", "no")
            println("Edit $editable_path, then press Enter to reload and refit.")
            readline(stdin)
        else
            println("Please answer y or n.")
        end
    end
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
        af = approve_airfoil(slice, i)
        push!(airfoils, af)
        editable_path = joinpath(OUTPUT_DIR, "airfoils", @sprintf("section_%02d.dat", i))
        push!(all_params, approve_cst(i, slice.r, editable_path))
    end

    cst_path = joinpath(OUTPUT_DIR, RUN_NAME * "_cst.csv")
    write_cst_table(cst_path, airfoils, all_params)
    println("Wrote CST table: $cst_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
