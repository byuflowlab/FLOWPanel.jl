import FLOWPanel: norm, dot, cross
import Meshes
import GeoIO
import FLOWPanel as pnl
import GeometricTools as gt
using LinearAlgebra: diag
using StaticArrays: SVector

# ============================================================
# TRAILING EDGE DETECTION
# ============================================================
function find_trailing_edge_nodes(nodes, cells; tol=1e-8)

    xs = nodes[1, :]
    ys = nodes[2, :]

    xmax = maximum(xs)

    candidates = findall(x -> isapprox(x, xmax; atol=tol), xs)

    return candidates[end], candidates[end-1]
end

# ============================================================
# BODY GENERATION
# ============================================================
function generate_body(
    meshfile::String,
    chord::Float64,
    span::Float64,
    bodytype::Type{<:pnl.RigidWakeBody};
    translate::NTuple{3,Float64} = (0.0, 0.0, 0.0),
    scaling::Float64 = 1.0,
    Vinf::AbstractVector{<:Real} = zeros(3),
    kerneloffset::Float64 = 1.0e-6
)

    magVinf = norm(Vinf)

    msh = GeoIO.load(meshfile).geometry
    msh = msh |> Meshes.Scale(scaling)
    msh = msh |> Meshes.Translate(translate)

    grid = pnl.gt.GridTriangleSurface(msh)

    CPoffset = 1e-6
    empty_shedding = zeros(Int, 6, 0)

    # build temporary body
    body = bodytype(
        grid,
        [empty_shedding];
        CPoffset,
        kerneloffset,
        flip_normals = false
    )

    nodes = body.nodes
    cells = body.cells

    # trailing edge
    firstnode, secondnode = find_trailing_edge_nodes(nodes, cells)

    shedding = pnl.calc_shedding_from_seed(
        nodes,
        cells,
        firstnode,
        secondnode
    )

    # rebuild with wake
    body = bodytype(
        nodes,
        cells,
        [shedding];
        CPoffset,
        kerneloffset,
        flip_normals = false,
        ensure_winding = false
    )

    # initialize wake direction
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(Vinf / magVinf, 1, size(body.Das[i], 2))
    end

    return body
end

# ============================================================
# POSTPROCESSING
# ============================================================
function postprocess!(bodies, Vinf, rho, chords, span, scaling::Float64=1.0)

    Dhat = Vinf / norm(Vinf)
    Shat = SVector(0.0, 1.0, 0.0)
    Lhat = cross(Dhat, Shat)

    Sref = sum(chord * span * scaling^2 for chord in chords)

    pnl.calcfield_U!(bodies, Vinf)
    pnl.calcfield_P!(bodies, norm(Vinf), rho)
    pnl.calcfield_F!(bodies)

    LDS = pnl.calcfield_LDS!(zeros(3,3), bodies, Lhat, Dhat, cross(Lhat, Dhat))

    nondim = 0.5 * rho * norm(Vinf)^2 * Sref

    CL = sign(dot(LDS[:,1], Lhat)) * norm(LDS[:,1]) / nondim
    CD = sign(dot(LDS[:,2], Dhat)) * norm(LDS[:,2]) / nondim

    return CL, CD
end

# ============================================================
# EXPERIMENT SETUP
# ============================================================
run_names = ["wing.msh", "surface.msh"]
files = [joinpath(pnl.examples_path, "wing_aileron", n) for n in run_names]

AOAs = [
-9.5879,
-5.9886,
-0.5368,
 3.1682,
 6.9262,
10.4725,
14.2306,
18.4650
]

m = 0.0254
magVinf = 117.3 * 12 * m
rho = 1.225

c_body1, c_body2 = 10.0, 2.0
b = 60.0
chords = [c_body1, c_body2]

trs = [
    (0.0, 0.0, 0.0),
    (9.8489*m, 0.0, -0.12*m)
]

kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}
bodytype = pnl.RigidWakeBody{kernel}

out_file = joinpath(pnl.examples_path, "wing_aileron", "timing_summary.csv")

# ============================================================
# SOLVER RUNNER (WRITES DIRECTLY)
# ============================================================
function run_solver(io, name, solver_builder, AOAs, experimental)

    CLs = Float64[]
    t_build_total = 0.0
    t_solve_total = 0.0
    nps_tot = 0
    experimental_AOA = experimental[:, 1]
    experimental_CL = experimental[:, 2]
    println("Running $name")

    sq_error = 0.0

    for (i, (AOA, expAOA)) in enumerate(zip(AOAs, experimental_AOA))

        Vinf = magVinf * [cosd(AOA), 0.0, sind(AOA)]

        bodies = tuple([
            generate_body(file, chord, b, bodytype;
                translate=tr,
                scaling=m,
                Vinf=Vinf
            )
            for (file, chord, tr) in zip(files, chords, trs)
        ]...)

        for body in bodies
            body.velocity .= 0.0
            pnl.apply_freestream!(body, Vinf)
        end

        solver = solver_builder(bodies)

        t_build, t_solve = pnl.solve!(bodies, solver)

        CL, CD = postprocess!(bodies, Vinf, rho, chords, b, m)

        push!(CLs, CL)
        t_build_total += t_build
        t_solve_total += t_solve

        if i == 1
            nps_tot = sum(b.ncells for b in bodies)
        end

        # RMS accumulation
        sq_error += (CL - experimental_CL[i])^2

        write(io,
            "$name,$nps_tot,$AOA,$CL,$CD,$t_build,$t_solve,$(t_build + t_solve)\n"
        )
    end

    rms = sqrt(sq_error / length(AOAs))

    write(io,
        "$name,SUMMARY,$rms,$nps_tot,$t_build_total,$t_solve_total,$(t_build_total + t_solve_total)\n"
    )

end

# ============================================================
# MAIN
# ============================================================
experimental = [
 -9.58790170132325  -0.34195250659630627
 -7.841209829867704 -0.20211081794194818
 -5.988657844990541 -0.05435356200527597
 -4.136105860113403  0.09340369393140424
 -2.3894139886577754 0.2358839050132029
 -0.536862003780719  0.38100263852242966
  1.3156899810964084 0.5261213720316649
  3.168241965973529  0.6686015831134591
  5.02079395085066   0.8084432717678132
  6.926275992438619  0.9535620052770526
  8.672967863894144  1.0907651715039615
 10.472589792060482  1.2279683377308745
 12.325141776937613  1.3546174142480258
 14.230623818525515  1.4759894459102951
 15.077504725897917  1.528759894459108
 16.400756143667294  1.3308707124010595
 18.465028355387517  1.2807387862796877
]


open(out_file, "w") do io

    write(io, "solver,nps,AOA,CL,CD,t_build,t_solve,total_time\n")

    run_solver(io, "BackslashCoupled",
        bodies -> pnl.BackslashCoupled(bodies), AOAs, experimental
    )

    run_solver(io, "BackslashIterative",
        bodies -> tuple(
            pnl.Backslash(bodies[1]),
            pnl.Backslash(bodies[2])
        ), AOAs, experimental
    )

    # run_solver(io, "KrylovSolver",
    #     bodies -> tuple(
    #         pnl.KrylovSolver(bodies[1]),
    #         pnl.KrylovSolver(bodies[2])
    #     )
    # )

    # run_solver(io, "KrylovSolver-FMM",
    #     bodies -> tuple(
    #         pnl.KrylovSolver(bodies[1]; backend=pnl.FastMultipoleBackend()),
    #         pnl.KrylovSolver(bodies[2]; backend=pnl.FastMultipoleBackend())
    #     )
    # )

    # run_solver(io, "FGSSolver",
    #     bodies -> tuple(
    #         pnl.FGSSolver(bodies[1]; leaf_size=10000),
    #         pnl.FGSSolver(bodies[2]; leaf_size=10000)
    #     ), AOAs
    # )
end

println("Saved results to: ", out_file)