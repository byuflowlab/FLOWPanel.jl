import FLOWPanel: norm, dot, cross
import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import Meshes
import GeoIO
import GeometricTools as gt
using LinearAlgebra: diag
using StaticArrays: SVector
using Rotations


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
    rotate::Float64 = 0.0,
    scaling::Float64 = 1.0,
    Vinf::AbstractVector{<:Real} = zeros(3),
    kerneloffset::Float64 = 1.0e-6
)

    magVinf = norm(Vinf)

    msh = GeoIO.load(meshfile).geometry
    msh = msh |> Meshes.Scale(scaling)
    msh = msh |> Meshes.Rotate(RotY(rotate * pi/180))
    msh = msh |> Meshes.Translate(translate)

    grid = gt.GridTriangleSurface(msh)

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
# run_names = ["wing_13_13.msh", "surface_13_13.msh"]
run_names = ["wing_20_20.msh", "surface_20_20.msh"]
# run_names = ["wing_40_40.msh", "surface_40_40.msh"]
# run_names = ["wing_60_60.msh", "surface_60_60.msh"]


files = [joinpath(pnl.examples_path, "wing_aileron", n) for n in run_names]

AOAs = [
# -9.58790170132325
-7.841209829867704
# -5.988657844990541
# -4.136105860113403
# -2.3894139886577754
-0.536862003780719
#  1.3156899810964084
#  3.168241965973529
#  5.02079395085066
 6.926275992438619
#  8.672967863894144
# 10.472589792060482
# 12.325141776937613
 14.230623818525515
#  15.077504725897917
# 16.400756143667294
]

# AOAs = [0.0]

m = 0.0254
magVinf = 117.3 * 12 * m
rho = 1.225

c_body1, c_body2 = 10.0, 2.0
b = 60.0
chords = [c_body1, c_body2]

trs = [
    (0.0, 0.0, 0.0),
    (9.8489*m, 0.0, -0.22*m)
]

rts = [0.0, 10.0]

kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}
bodytype = pnl.RigidWakeBody{kernel}

out_file = joinpath(pnl.examples_path, "wing_aileron", "krylov_solvers.csv")

# ============================================================
# SOLVER RUNNER (WRITES DIRECTLY)
# ============================================================
function run_solver(io, name, solver_builder, AOAs, experimental; paraview::Bool=false, solve_kwargs=NamedTuple())

    CLs = Float64[]
    nps_tot = 0
    experimental_AOA = experimental[:, 1]
    experimental_CL = experimental[:, 2]
    println("Running $name")

    sq_error = 0.0
    m = 0.0254
    tol = 1e-1

    #------------------- GENERATE END PLATES ----------------------------------------------
    # Ground-plate geometry depends only on (b, m, tol), not on AOA, so build
    # them once and reuse the same body objects across the whole AOA sweep.
    left_center = [12.0*0.5*m, -b/2 * m - tol, 0.0]
    left_normal = [0.0, 1.0, 0.0]
    left_radius = 12.0 * 0.1
    left_plate = pnl.FlatGround(left_center, left_normal, left_radius; panel_length=12.0*0.1*m)
    right_center = [12.0*0.5*m, b/2* m + tol, 0.0]
    right_normal = [0.0, -1.0, 0.0]
    right_radius = 12.0 * 0.1
    right_plate = pnl.FlatGround(right_center, right_normal, right_radius; panel_length=12.0*0.1*m)
    println(left_plate.ncells + right_plate.ncells)

    for (i, (AOA, expCL)) in enumerate(zip(AOAs, experimental_CL))

        Vinf = magVinf * [cosd(AOA), 0.0, sind(AOA)]

        # Wing/aileron bodies carry a rigid wake whose direction (Das) is set
        # from Vinf, so their G matrix genuinely depends on AOA and they must
        # be rebuilt every iteration.
        bodies = tuple([
            generate_body(file, chord, b, bodytype;
                translate=tr,
                scaling=m,
                Vinf=Vinf,
                rotate=rt
            )
            for (file, chord, tr, rt) in zip(files, chords, trs, rts)
        ]..., left_plate, right_plate)

        for body in bodies
            body.velocity .= 0.0
            pnl.apply_freestream!(body, Vinf)
        end

        @show p = sum(b.ncells for b in bodies)

        solver, t_build = solver_builder(bodies)

        t_build_solve, t_solve = pnl.solve!(bodies, solver; solve_kwargs...)
        t_build += t_build_solve
        @show t_build, t_solve

        CL, CD = postprocess!(bodies, Vinf, rho, chords, b, m)

        push!(CLs, CL)

        if i == 1
            nps_tot = sum(b.ncells for b in bodies)
        end

        if i == 1 && paraview
            filestr1 = pnl.write_vtk(joinpath("examples", "wing_val"), bodies[1], 0, 0.0)
            files1 = split(filestr1, ", ")
            pvd1 = first(filter(f -> endswith(f, ".pvd"), files1))

            filestr2 = pnl.write_vtk(joinpath("examples", "surface_val"), bodies[2], 0, 0.0)
            files2 = split(filestr2, ", ")
            pvd2 = first(filter(f -> endswith(f, ".pvd"), files2))

            filestr3 = pnl.write_vtk(joinpath("examples", "left_val"), bodies[3], 0, 0.0)
            files3 = split(filestr3, ", ")
            pvd3 = first(filter(f -> endswith(f, ".pvd"), files3))

            filestr4 = pnl.write_vtk(joinpath("examples", "right_val"), bodies[4], 0, 0.0)
            files4 = split(filestr4, ", ")
            pvd4 = first(filter(f -> endswith(f, ".pvd"), files4))

            run(`paraview $pvd1 $pvd2 $pvd3 $pvd4`, wait=false)
        end

        # RMS accumulation
        # sq_error += (CL - expCL)^2

        write(io,
            "$name,$nps_tot,$AOA,$CL,$CD,$t_build,$t_solve,$(t_build + t_solve)\n"
        )
    end

end

# ============================================================
# MAIN
# ============================================================
experimental = [
#  -9.58790170132325  -0.34195250659630627
 -7.841209829867704 -0.20211081794194818
#  -5.988657844990541 -0.05435356200527597
#  -4.136105860113403  0.09340369393140424
#  -2.3894139886577754 0.2358839050132029
 -0.536862003780719  0.38100263852242966
#   1.3156899810964084 0.5261213720316649
#   3.168241965973529  0.6686015831134591
#   5.02079395085066   0.8084432717678132
  6.926275992438619  0.9535620052770526
#   8.672967863894144  1.0907651715039615
#  10.472589792060482  1.2279683377308745
#  12.325141776937613  1.3546174142480258
 14.230623818525515  1.4759894459102951
#  15.077504725897917  1.528759894459108
#  16.400756143667294  1.3308707124010595
]

is_new_file = !isfile(out_file) || filesize(out_file) == 0

open(out_file, "a") do io
    if is_new_file
        write(io, "solver,nps,AOA,CL,CD,t_build,t_solve,total_time\n")
    end

    # BackslashCoupled(bodies) only allocates dummy placeholder storage (not
    # the real G matrix, which solve! assembles and factors itself, timed
    # internally), so no build time is attributable to the constructor here.
    # The solver object itself is safe to cache across the AOA sweep (rather
    # than reconstructing the dummy placeholder every AOA) now that solve!
    # explicitly zeros solver.rhs before building the RHS -- see the
    # fill!(solver.rhs, 0) fix in src/FLOWPanel_solver.jl.
    # coupled_solver = Ref{Any}(nothing)
    # run_solver(io, "BackslashCoupled",
    #     bodies -> begin
    #         if coupled_solver[] === nothing
    #             coupled_solver[] = pnl.BackslashCoupled(bodies)
    #         end
    #         (coupled_solver[], 0.0)
    #     end, AOAs, experimental
    # )

    # Unlike BackslashCoupled, Backslash(body) assembles and LU-factors G
    # eagerly in its constructor (solve! reuses that G unchanged, so its own
    # returned build time is 0 here) -- so the G-matrix build time for the
    # iterative solvers has to be measured around the constructor calls
    # themselves.
    #
    # The two ground-plate bodies (bodies[3], bodies[4]) are the same objects
    # every AOA iteration and their G matrix is AOA-invariant, so their
    # Backslash solvers are built once and cached here (their build time is
    # only counted on the AOA where they're actually built), while the wing
    # solvers (bodies[1], bodies[2]) are rebuilt every AOA since their G
    # depends on wake direction (Das).
    # plate_solvers = Ref{Any}(nothing)
    # run_solver(io, "BackslashIterative",
    #     bodies -> begin
    #         t_build = 0.0
    #         if plate_solvers[] === nothing
    #             t_build += @elapsed begin
    #                 plate_solvers[] = (pnl.Backslash(bodies[3]), pnl.Backslash(bodies[4]))
    #             end
    #         end
    #         wing_solvers = nothing
    #         t_build += @elapsed begin
    #             wing_solvers = (pnl.Backslash(bodies[1]), pnl.Backslash(bodies[2]))
    #         end
    #         (tuple(wing_solvers..., plate_solvers[]...), t_build)
    #     end, AOAs, experimental
    # )

    # KrylovSolver is a per-body solver, so it plugs into the same Gauss-Seidel
    # outer loop as BackslashIterative (solve!(bodies::Tuple, solvers::Tuple)),
    # which needs solver_builder to return (tuple_of_solvers, t_build) --
    # returning the bare tuple of solvers here made `solver, t_build =
    # solver_builder(bodies)` silently truncate to solver=KrylovSolver(bodies[1]),
    # t_build=KrylovSolver(bodies[2]), which then failed solve!(bodies, solver)
    # with "solve! not implemented for body of type Tuple{...}".
    # solver.backend (FastMultipoleBackend, KrylovSolver's own default) only
    # covers each body's own GMRES operator. The outer Gauss-Seidel loop
    # (solve!(bodies::Tuple, solvers::Tuple)) separately defaults its own
    # `backend` to DirectBackend() for the cross-body coupling evaluation and
    # for the per-body fixed-source RHS -- _solve! never sees or uses that
    # value for the GMRES operator itself, but it does reach those two other
    # spots. Pass FastMultipoleBackend() here explicitly so every evaluation
    # in the solve, not just the GMRES iterations, uses FMM.
    #
    # WARM-START (what/why/what it returns):
    #   What: each KrylovSolver persists its own `unabbreviated_strengths`
    #   (the last strength vector it solved for). _solve! (src/FLOWPanel_solver.jl)
    #   now always passes that stored vector to `Krylov.warm_start!` as GMRES's
    #   initial guess before calling `Krylov.krylov_solve!`. This is not a
    #   flag/kwarg -- it's unconditional solver behavior. There is nothing to
    #   toggle here; it falls out for free from these 4 solver objects being
    #   built once per AOA (below) and then reused across every outer
    #   Gauss-Seidel sweep within that AOA's solve!(bodies, solvers) call.
    #   Why: the outer Gauss-Seidel loop re-solves every body's GMRES system
    #   from scratch on each of up to `max_outer_iterations=50` sweeps as the
    #   bodies converge on each other's coupling influence. Before this change,
    #   every sweep started GMRES from an all-zero guess even though the
    #   previous sweep's answer is already a very good guess (the outer
    #   residual is shrinking sweep to sweep). Warm-starting should let GMRES
    #   converge in fewer iterations per sweep as the outer loop tightens,
    #   which is a pure win in wall-clock time if the diagnostics below confirm
    #   it -- with no accuracy cost, since it only changes the initial guess,
    #   not the converged solution. On the very first sweep of the very first
    #   AOA, `unabbreviated_strengths` is still all-zeros (the constructor's
    #   default), so behavior there is identical to no warm-start.
    #   What it returns: same CSV row shape as every other solver here
    #   (solver,nps,AOA,CL,CD,t_build,t_solve,total_time appended to
    #   krylov_solvers.csv), so it can be compared directly against the
    #   existing "KrylovSolver"/"KrylovCoupled" rows for both accuracy (CL vs
    #   experimental) and timing. It ALSO prints one diagnostic line per body
    #   per outer sweep to stdout (not the CSV) -- `||x0||` (norm of the
    #   warm-start guess), `niter`/`solved` (GMRES iterations and convergence
    #   flag), and `ts` (that solve's wall time) -- plus one line per outer
    #   sweep for the cross-body `influence!` cost (`t_influence`). These are
    #   the TEMPORARY diagnostics added to src/FLOWPanel_solver.jl for the
    #   warm-start investigation; remove them there once this is confirmed.
    try
        run_solver(io, "KrylovSolver",
            bodies -> (
                tuple(
                    pnl.KrylovSolver(bodies[1]),
                    pnl.KrylovSolver(bodies[2]),
                    pnl.KrylovSolver(bodies[3]),
                    pnl.KrylovSolver(bodies[4])
                ),
                0.0
            ), AOAs, experimental;
            solve_kwargs=(backend=pnl.FastMultipoleBackend(),)
        )
    catch e
        println("Error running KrylovSolver: ", e)
    end

    # KrylovCoupled solves the whole tuple of bodies as one coupled system, so
    # it goes through its own dedicated solve!(bodies::Tuple,
    # solver::KrylovCoupled) method rather than the Gauss-Seidel outer loop --
    # solver_builder just needs to return (solver, t_build) directly, no inner
    # tuple-of-solvers like KrylovSolver above.
    # try
    #     run_solver(io, "KrylovCoupled",
    #         bodies -> (pnl.KrylovCoupled(bodies), 0.0), AOAs, experimental
    #     )
    # catch e
    #     println("Error running KrylovCoupled: ", e)
    # end

    # TEMPORARY comparison: is FastMultipoleBackend's per-call tree-rebuild
    # overhead (expansion_order=10, rebuilt from scratch every single GMRES
    # iteration since panel strengths change each time) actually worth it at
    # this small panel count (nps~2572), or is plain O(N^2) direct summation
    # faster in absolute terms here? Same solver/tolerances as "KrylovCoupled"
    # above, only the backend differs.
    # try
    #     run_solver(io, "KrylovCoupled-Direct",
    #         bodies -> (pnl.KrylovCoupled(bodies; backend=pnl.DirectBackend()), 0.0), AOAs, experimental;
    #         solve_kwargs=(backend=pnl.DirectBackend(),)
    #     )
    # catch e
    #     println("Error running KrylovCoupled-Direct: ", e)
    # end

    # try
    #     run_solver(io, "KrylovCoupled-FMM",
    #         bodies -> (pnl.KrylovCoupled(bodies; backend=pnl.FastMultipoleBackend()), 0.0), AOAs, experimental;
    #         solve_kwargs=(backend=pnl.FastMultipoleBackend(),)
    #     )
    # catch e
    #     println("Error running KrylovCoupled-FMM: ", e)
    # end

    # try
    #     run_solver(io, "KrylovSolver-FMM",
    #         bodies -> tuple(
    #             pnl.KrylovSolver(bodies[1]; backend=pnl.FastMultipoleBackend()),
    #             pnl.KrylovSolver(bodies[2]; backend=pnl.FastMultipoleBackend()),
    #             pnl.KrylovSolver(bodies[3]; backend=pnl.FastMultipoleBackend()),
    #             pnl.KrylovSolver(bodies[4]; backend=pnl.FastMultipoleBackend()),
    #         ), AOAs, experimental
    #     )
    # catch e
    #     println("Error running KrylovSolver-FMM: ", e)
    # end

    # try
    #     run_solver(io, "FGSSolver",
    #         bodies -> tuple(
    #             pnl.FGSSolver(bodies[1]; leaf_size=10000),
    #             pnl.FGSSolver(bodies[2]; leaf_size=10000),
    #             pnl.FGSSolver(bodies[3]; leaf_size=10000),
    #             pnl.FGSSolver(bodies[4]; leaf_size=10000),
    #         ), AOAs, experimental
    #     )
    # catch e
    #     println("Error running FGSSolver: ", e)
    # end
end

println("Saved results to: ", out_file)
