# Semi-infinite Dirichlet solve for the top-level quad-facet "STL" wing.
# The file uses four vertices per outer loop, so GeoIO's triangle STL loader
# cannot preserve it.  This diagnostic explicitly triangulates every loop.

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: cross, dot, norm
import Printf: @printf
import TOML

include(joinpath(@__DIR__, "dirichlet_solve.jl"))

function load_quad_facet_stl(path)
    node_index = Dict{NTuple{3,Float64},Int}()
    points = NTuple{3,Float64}[]
    triangles = NTuple{3,Int}[]
    loop_vertices = NTuple{3,Float64}[]

    get_node!(xyz) = get!(node_index, xyz) do
        push!(points, xyz)
        length(points)
    end

    for raw_line in eachline(path)
        line = strip(raw_line)
        if startswith(line, "outer loop")
            empty!(loop_vertices)
        elseif startswith(line, "vertex")
            values = split(line)
            length(values) == 4 || error("invalid vertex line: $raw_line")
            push!(loop_vertices, tuple(parse.(Float64, values[2:4])...))
        elseif startswith(line, "endloop")
            length(loop_vertices) == 4 || error(
                "expected four vertices per loop, found $(length(loop_vertices))")
            ids = get_node!.(loop_vertices)
            for local_triangle in ((1, 2, 3), (1, 3, 4))
                tri = (ids[local_triangle[1]], ids[local_triangle[2]],
                    ids[local_triangle[3]])
                a, b, c = points[tri[1]], points[tri[2]], points[tri[3]]
                ux, uy, uz = b[1] - a[1], b[2] - a[2], b[3] - a[3]
                vx, vy, vz = c[1] - a[1], c[2] - a[2], c[3] - a[3]
                cx, cy, cz = uy * vz - uz * vy, uz * vx - ux * vz,
                    ux * vy - uy * vx
                cx * cx + cy * cy + cz * cz > 1e-28 && push!(triangles, tri)
            end
        end
    end

    nodes = zeros(Float64, 3, length(points))
    for (i, p) in pairs(points)
        nodes[:, i] .= p
    end
    cells = zeros(Int, 3, length(triangles))
    for (i, t) in pairs(triangles)
        cells[:, i] .= t
    end
    return nodes, cells
end

function build_stl_wing(path)
    nodes, cells = load_quad_facet_stl(path)
    # The file uses LE x=0 and TE x=-1. Reflect x so the established diagnostic
    # convention is LE x=0, TE x=1 with freestream primarily in +x.
    nodes[1, :] .*= -1

    bodytype = pnl.RigidWakeBody{
        Union{pnl.ConstantSource,pnl.ConstantDoublet}, 2, Float64, true}
    options = (;
        kerneloffset=1e-6,
        kernelcutoff=1e-12,
        semiinfinite_wake=true,
        watertight=true,
    )
    base = bodytype(nodes, cells, zeros(Int, 6, 0); options...)
    watertight, _ = pnl.iswatertight(base.nodes, base.cells)
    watertight || error("explicitly triangulated STL is not watertight")

    x_te = maximum(base.nodes[1, :])
    tol = 1e-8
    te_nodes = findall(i -> isapprox(base.nodes[1, i], x_te; atol=tol, rtol=0) &&
        isapprox(base.nodes[3, i], 0.0; atol=tol, rtol=0), axes(base.nodes, 2))
    sort!(te_nodes, by=i -> base.nodes[2, i])
    length(te_nodes) >= 2 || error("could not identify trailing-edge chain")
    lower = [x_te - tol, minimum(base.nodes[2, te_nodes]) - tol, -tol]
    upper = [x_te + tol, maximum(base.nodes[2, te_nodes]) + tol, tol]
    shedding = pnl.calc_shedding_from_seed(base.nodes, base.cells,
        te_nodes[1], te_nodes[2]; end_node=te_nodes[end], bbox=(lower, upper),
        normal_jump_tol=1.0, max_turn_angle=pi / 2)
    body = bodytype(copy(base.nodes), copy(base.cells), [shedding]; options...)
    return body
end

function run_stl_wing_case(;
        mesh_path=joinpath(@__DIR__, "naca0015_wing.stl"),
        output_path=joinpath(
            @__DIR__, "..", "..", "data", "dirichlet_solve", "naca0015_wing_stl"))
    mkpath(output_path)
    U = 1.0
    alpha_deg = 12.0
    rho = 1.225
    Sref = 6.576
    cref = 1.0
    Uinf = _uinf_from_alpha(U, alpha_deg)

    body = build_stl_wing(mesh_path)
    set_wake_Das!(body, Uinf)
    C = _build_kutta_map(body)
    norm(C * ones(body.ncells), Inf) == 0 || error("invalid paired Kutta map")

    direct = pnl.DirectBackend()
    solver = pnl.Backslash(body)
    pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false,
        backend=direct)
    force = pnl.ForceMonitor(1, 1;
        i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, cref),
        correct_kuttacondition=false,
        verbose=false)
    name = "naca0015_wing_stl_semiinfinite"
    pnl.steady!(body, pnl.ReferenceFrame(body), Uinf;
        body_solvers=solver, backend=direct, backend_solve=direct,
        backend_system=direct, monitors=(pressure, force), path=output_path,
        name, dt=0.5, clean_files=true, verbose=true)

    mu = copy(view(body.strength, :, 2))
    residual_abs = norm(_lu_matrix_product(solver.Glu, mu) - solver.rhs)
    residual_rel = residual_abs / max(norm(solver.rhs), eps(Float64))
    Lhat, _, _ = _lift_drag_span_directions(Uinf)
    CL = dot(SVector{3,Float64}(force.force[:, 1]), Lhat)
    lift_N = 0.5 * rho * U^2 * Sref * CL

    open(joinpath(output_path, "result.csv"), "w") do io
        println(io, "ncells,nnodes,nshedding_edges,CL,dimensional_lift_N,linear_residual_abs,linear_residual_rel")
        println(io, join((body.ncells, body.nnodes, size(C, 1), CL, lift_N,
            residual_abs, residual_rel), ','))
    end
    open(joinpath(output_path, "config.toml"), "w") do io
        TOML.print(io, Dict(
            "mesh" => abspath(mesh_path),
            "mesh_import" => "explicit triangulation of 1680 four-vertex facet loops",
            "coordinate_transform" => "x := -x (LE 0, TE 1)",
            "S_ref_m2" => Sref,
            "b_ref_m" => 6.576,
            "c_ref_m" => cref,
            "U_ref_m_s" => U,
            "rho_ref_kg_m3" => rho,
            "alpha_deg" => alpha_deg,
            "x_cg_m" => -0.25,
            "p_ref_Pa" => 101325.0,
            "mu_ref_Pa_s" => 1.789e-5,
            "M_ref" => 0.0,
            "trailing_edge_angle_deg" => 150.0,
            "trailing_edge_interpretation" => "not applied; explicit TE topology",
            "watertight" => true,
        ); sorted=true)
    end

    @printf("STL wing: nodes=%d, panels=%d, shedding edges=%d\n",
        body.nnodes, body.ncells, size(C, 1))
    @printf("  CL=%.10f, lift=%.10f N, relative residual=%.3e\n",
        CL, lift_N, residual_rel)
    println("  output=$(abspath(output_path))")
    return (; body, solver, C, CL, lift_N, residual_abs, residual_rel)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_stl_wing_case()
end
