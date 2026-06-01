import FLOWPanel as pnl
import GeoIO
using LinearAlgebra: norm

const DEFAULT_SIMPLE_WING_CAPPED_DIRICHLET_SOLVERS = (
    :backslash,
    :krylov,
    :fgs,
    :backslash_coupled,
    :krylov_coupled,
)

function build_simple_wing_capped_dirichlet(;
    meshfile=joinpath(pnl.examples_path, "data", "naca0012_nc133.msh"),
    AOA=5.88,
    magVinf=56.0,
    stretch=1.0,
    cp_offset=1e-10,
    kernel_offset=1e-3,
    trailing_edge_x=nothing,
    trailing_edge_y_start=nothing,
    trailing_edge_y_stop=nothing,
)
    msh = GeoIO.load(meshfile).geometry
    nodes, cells = pnl.meshes2nodes_cells(msh)
    nodes[2, :] .*= stretch

    xte = isnothing(trailing_edge_x) ? maximum(nodes[1, :]) : trailing_edge_x
    ystart = isnothing(trailing_edge_y_start) ? minimum(nodes[2, :]) : trailing_edge_y_start
    ystop = isnothing(trailing_edge_y_stop) ? maximum(nodes[2, :]) : trailing_edge_y_stop
    span = ystop - ystart

    trailingedge = zeros(eltype(nodes), 3, 10_000)
    trailingedge[1, :] .= xte
    trailingedge[2, :] .= range(ystart, stop=ystop, length=size(trailingedge, 2))

    Vinf = magVinf .* [cosd(AOA), 0.0, sind(AOA)]
    bodytype = pnl.RigidWakeBody{
        Union{pnl.ConstantSource, pnl.ConstantDoublet},
        2,
        Float64,
        true,
    }

    shedding = pnl.calc_shedding(nodes, cells, trailingedge; tolerance=0.001 * span)
    body = bodytype(
        nodes,
        cells,
        shedding;
        cp_outer=(cp_offset > 0),
        kerneloffset=kernel_offset,
        ensure_winding=true,
        semiinfinite_wake=true,
    )

    shedding = pnl.calc_shedding(body.nodes, body.cells, trailingedge; tolerance=0.001 * span)
    body = bodytype(
        body.nodes,
        body.cells,
        shedding;
        cp_outer=(cp_offset > 0),
        kerneloffset=kernel_offset,
        ensure_winding=true,
        semiinfinite_wake=true,
    )

    wake_direction = reshape(Vinf ./ norm(Vinf), :, 1)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
    end

    pnl.apply_freestream!(body, Vinf)

    return body, Vinf
end

function dirichlet_potential_max_residual(body; backend=pnl.DirectBackend(), cp_off=-1e-10)
    potential_ext = copy(body.potential)
    cp_outer_old = body.cp_outer

    try
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.potential .= 0.0
        pnl.influence!((body,), (body,), backend; scalar_potential=true, velocity=false)
        return maximum(abs.(body.potential .+ potential_ext))
    finally
        body.cp_outer = cp_outer_old
        body.potential .= potential_ext
    end
end

function make_simple_wing_dirichlet_solver(body, solver_name::Symbol; backend=pnl.DirectBackend())
    if solver_name == :backslash
        return pnl.Backslash(body)
    elseif solver_name == :krylov
        return pnl.KrylovSolver(body; backend, atol=1e-10, rtol=1e-10, itmax=200)
    elseif solver_name == :fgs
        return pnl.FGSSolver(body; leaf_size=10_000, tolerance=1e-10, max_iterations=200)
    elseif solver_name == :backslash_coupled
        return pnl.BackslashCoupled((body,))
    elseif solver_name == :krylov_coupled
        return pnl.KrylovCoupled((body,); backend, atol=1e-10, rtol=1e-10, itmax=200)
    end

    error("Unsupported solver name: $solver_name")
end

function solve_simple_wing_capped_dirichlet!(
    body,
    solver_name::Symbol;
    backend=pnl.DirectBackend(),
)
    solver = make_simple_wing_dirichlet_solver(body, solver_name; backend)

    if solver_name == :backslash_coupled
        pnl.solve!((body,), solver; backend, update_G=true)
    elseif solver_name == :krylov_coupled
        pnl.solve!((body,), solver; backend)
    elseif solver_name == :fgs
        pnl.solve!(body, solver)
    else
        pnl.solve!(body, solver; backend)
    end

    residual = dirichlet_potential_max_residual(body; backend)
    return solver, residual
end

function run_simple_wing_capped_dirichlet(;
    meshfile=joinpath(pnl.examples_path, "data", "naca0012_nc133.msh"),
    solvers=DEFAULT_SIMPLE_WING_CAPPED_DIRICHLET_SOLVERS,
    backend=pnl.DirectBackend(),
    build_kwargs...,
)
    reference_strength = nothing
    results = NamedTuple[]

    for solver_name in solvers
        body, Vinf = build_simple_wing_capped_dirichlet(; meshfile, build_kwargs...)
        solver, residual = solve_simple_wing_capped_dirichlet!(body, solver_name; backend)
        mu = copy(view(body.strength, :, 2))
        sigma = copy(view(body.strength, :, 1))

        if reference_strength === nothing
            reference_strength = mu
        end

        push!(results, (
            solver=solver_name,
            residual=residual,
            doublet_relerr=norm(mu .- reference_strength) / max(norm(reference_strength), eps()),
            source_relerr=norm(sigma .- vec(sum(body.velocity .* body.normals; dims=1)) .* -1.0) /
                          max(norm(sigma), eps()),
            npanels=body.ncells,
            solver_object=solver,
            Vinf=Vinf,
        ))
    end

    return results
end

function main()
    println()
    println("#===== SIMPLE WING CAPPED DIRICHLET =====#")
    println("Solving capped wing without endplates using the Dirichlet formulation.")

    results = run_simple_wing_capped_dirichlet()
    println("solver\t\tresidual\t\tdoublet_relerr")
    for result in results
        println("$(result.solver)\t$(result.residual)\t$(result.doublet_relerr)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
