#=
When calling solve!, assume that the wake-, kinematic-, and freestream-induced velocities have already been applied,
and that fast-access panels vector has been updated
=#

abstract type AbstractBoundaryCondition end

abstract type FlowTangency <: AbstractBoundaryCondition end

abstract type AbstractFormulation end

abstract type DirectNeumann <: AbstractFormulation end

abstract type Dirichlet <: AbstractFormulation end

abstract type Scheme{F<:AbstractFormulation, BC<:AbstractBoundaryCondition} end

abstract type AbstractSolver end

#####
##### auxilliary functions
#####

function update_influence_matrix!(influence_matrix,
                                    panels::AbstractPanels{K,<:Any,<:Any,<:Any},
                                    ::Type{Scheme{DirectNeumann, FlowTangency}};
                                    panel_indices=1:length(panels.panels),
                                    set_unit_strength=false
                                    ) where {K <: AbstractKernel{1}}
    # check matrix size
    @assert size(influence_matrix,1) == size(influence_matrix,2) ""*
        "influence matrix should be square"

    @assert size(influence_matrix,1) == length(panel_indices) ""*
        "influence matrix size $(size(influence_matrix,1)) inconsistent with"*
        " number of panels $(length(panels.panels))"

    # set panels to unit strength
    if set_unit_strength
        panels_2_grid_strength!(panels; panel_indices)
        set_unit_strength!(panels; panel_indices)
    end

    kernel = K()

    # update influence matrix
    for (j,i_source) in enumerate(panel_indices)
        source_panel = panels.panels[i_source]
        for (i,i_target) in enumerate(panel_indices)
            target_panel = panels.panels[i_target]
            _, v, _ = induced(target_panel.control_point, source_panel, kernel,
                                      DerivativesSwitch(false,true,false))
            influence_matrix[i, j] = dot(v, target_panel.normal)
        end
    end

    # restore panel strengths
    set_unit_strength && grid_2_panels_strength!(panels; panel_indices)
end


"""
# Arguments

* `right_hand_side::Vector{Float64}`: vector of length equal to the number of panels representing the right hand side of the boundary condition to be satified
* `panels::AbstractPanels`: the panels array object
* `scheme::Scheme`: the scheme used for setting the boundary condition

# Keyword Arguments

* `panel_indices::UnitRange{Int64}`: indices of the panels whose right-hand-side is to be updated (might not be the entire array)

"""
function update_right_hand_side!(right_hand_side, panels::AbstractPanels, ::Type{Scheme{DirectNeumann, FlowTangency}}, freestream=SVector{3,Float64}(0,0,0); panel_indices = 1:length(panels.panels), reset=true)
    # check RHS size
    @assert length(right_hand_side) == length(panels.panels) "length of RHS $(length(right_hand_side)) inconsistent with number of panels $(length(panels.panels))"

    # reset RHS
    reset && (right_hand_side[panel_indices] .= zero(eltype(right_hand_side)))

    # update right hand side
    for i in panel_indices
        panel = panels.panels[i]
        velocity = panels.velocity[i]
        right_hand_side[i] -= dot(freestream + velocity, panel.normal)
    end
end



#####
##### LU Decomposition
#####

struct LUDecomposition{TA,TF,S<:Scheme} <: AbstractSolver
    influence_matrix::TA
    right_hand_side::Vector{TF}
    strengths::Vector{TF}
end

LUDecomposition(influence_matrix, right_hand_side, strengths, scheme) =
LUDecomposition{typeof(influence_matrix), eltype(right_hand_side), scheme}(influence_matrix, right_hand_side, strengths)

function LUDecomposition(panels::AbstractPanels{K,TF,<:Any,<:Any}, scheme;
                            save_lu=false
                            ) where {TF, K<:AbstractKernel{1}}
    # preallocate memory
    n_panels = length(panels.panels)
    influence_matrix = zeros(TF,n_panels,n_panels)
    right_hand_side = zeros(TF,n_panels)

    # update influence matrix
    update_influence_matrix!(influence_matrix, panels, scheme; set_unit_strength=true)

    # get LU decomposition
    lu_decomposition = nothing#lu(influence_matrix)
    save_lu && (influence_matrix = lu!(influence_matrix))

    # initialize strengths
    strengths = zeros(TF,n_panels)

    return LUDecomposition(influence_matrix, right_hand_side, strengths, scheme)
end

function LUDecomposition_benchmark(panels::AbstractPanels{K,TF,<:Any,<:Any},
                                    scheme; save_lu=true
                                    ) where {TF, K<:AbstractKernel{1}}
    # preallocate memory
    t_alloc = @elapsed begin
        n_panels = length(panels.panels)
        influence_matrix = zeros(TF,n_panels,n_panels)
        right_hand_side = zeros(TF,n_panels)
    end

    # update influence matrix
    t_aic = @elapsed update_influence_matrix!(influence_matrix, panels, scheme)
    t_lu = @elapsed save_lu && (influence_matrix = lu!(influence_matrix))

    # initialize strengths
    t_alloc += @elapsed strengths = zeros(TF,n_panels)

    return LUDecomposition(influence_matrix, right_hand_side, strengths, scheme), t_aic, t_lu, t_alloc
end

@inline function get_strength(strengths::Vector, i, old_strength::SVector{1,<:Any})
    return SVector{1}(strengths[i])
end

@inline function get_strength(strengths::Vector, i, old_strength::SVector{2,<:Any})
    return SVector{2}(old_strength[1], strengths[i])
end

function solve!(panels::AbstractPanels{K,<:Any,<:Any,<:Any},
                solver::LUDecomposition{<:Any,<:Any,S}, dt=0.0
                ) where {S, K<:AbstractKernel{1}}
    # unpack
    # (; influence_matrix, right_hand_side, strengths) = solver
    influence_matrix = solver.influence_matrix
    right_hand_side = solver.right_hand_side
    strengths = solver.strengths

    # apply freestream/panel velocity
    update_right_hand_side!(right_hand_side, panels, S)

    # solver for strengths
    strengths .= influence_matrix \ right_hand_side

    # update panels
    for i in 1:length(panels.panels)
        # (; vertices, control_point, normal, strength, radius) = panels.panels[i]
        vertices = panels.panels[i].vertices
        control_point = panels.panels[i].control_point
        normal = panels.panels[i].normal
        strength = panels.panels[i].strength
        radius = panels.panels[i].radius
        # strength = get_strength(strengths, i, old_strength)

        panels.panels[i] = Panel(vertices, control_point, normal, strength, radius)
        panels.strengths[i] = strength
    end
end

#------- iterative solvers -------#

struct FastLinearOperator{TP,TT,TDL,TDS,scheme}
    panels::TP
    fmm_toggle::Bool
    reuse_tree::Bool
    save_residual::Bool
    tree::TT
    m2l_list::Vector{SVector{2,Int32}}
    direct_list::TDL
    derivatives_switch::TDS
    expansion_order::Int64
    leaf_size::SVector{1,Int64}
    multipole_threshold::Float64
end

Base.eltype(flo::FastLinearOperator) = eltype(flo.panels)

function LinearAlgebra.mul!(C, flo::FastLinearOperator, B, α, β)
    flo(C, B, α, β; expansion_order=flo.expansion_order, leaf_size_source=flo.leaf_size,
       multipole_threshold=flo.multipole_threshold)
end

Base.size(flo::FastLinearOperator) = (length(flo.panels.panels), length(flo.panels.panels))

struct IterativeSolver{TF,TS<:Krylov.KrylovSolver,TA<:Union{Matrix{TF},FastLinearOperator},S<:Scheme} <: AbstractSolver
    solver::TS
    A::TA
    right_hand_side::Vector{TF}
end

IterativeSolver(solver, A, right_hand_side, scheme) =
    IterativeSolver{eltype(right_hand_side), typeof(solver), typeof(A), scheme}(solver, A, right_hand_side)

function IterativeSolver(panels::AbstractPanels{K,TF,<:Any,<:Any}, scheme,
                            krylov_solver::Type{<:Krylov.KrylovSolver}=Krylov.GmresSolver
                            ) where {TF, K<:AbstractKernel{1}}

    # preallocate memory
    n_panels = length(panels.panels)
    influence_matrix = zeros(TF,n_panels,n_panels)
    right_hand_side = zeros(TF,n_panels)

    # update influence matrix
    update_influence_matrix!(influence_matrix, panels, scheme)

    # create krylov solver
    solver = krylov_solver(influence_matrix, right_hand_side)

    return IterativeSolver(solver, influence_matrix, right_hand_side, scheme)
end

function IterativeSolver_benchmark(panels::AbstractPanels{K,TF,<:Any,<:Any},
                                    scheme,
                                    krylov_solver::Type{<:Krylov.KrylovSolver}=Krylov.GmresSolver
                                    ) where {TF, K<:AbstractKernel{1}}

    # preallocate memory
    t_alloc = @elapsed begin
        n_panels = length(panels.panels)
        influence_matrix = zeros(TF,n_panels,n_panels)
        right_hand_side = zeros(TF,n_panels)
    end

    # update influence matrix
    t_aic = @elapsed update_influence_matrix!(influence_matrix, panels, scheme)

    # create krylov solver
    t_alloc += @elapsed solver = krylov_solver(influence_matrix, right_hand_side)

    return IterativeSolver(solver, influence_matrix, right_hand_side, scheme), t_aic, t_alloc
end

function get_fmm_objects(panels::AbstractPanels{<:Any,TF,<:Any,<:Any}, expansion_order, leaf_size, multipole_threshold, self_induced, reuse_tree) where TF
    switch = FastMultipole.DerivativesSwitch(false, true, false, panels)
    if reuse_tree
        # tree = FastMultipole.Tree(panels; expansion_order, leaf_size, shrink_recenter=true)
        # farfield, nearfield = true, true
        # m2l_list, direct_list = FastMultipole.build_interaction_lists(tree.branches, tree.branches, leaf_size, multipole_threshold, farfield, nearfield, self_induced)
        # direct_list = FastMultipole.InteractionList(direct_list, panels, tree, panels, tree, switch)
        @warn "reuse_tree was requested but is not currently supported; proceeding without"
    # else
    end
        tree = FastMultipole.Tree(panels; expansion_order, leaf_size, shrink_recenter=true)
        m2l_list = Vector{SVector{2,Int32}}(undef,0)
        direct_list = FastMultipole.InteractionList(Matrix{TF}[], TF[], TF[], m2l_list)
    # end

    return switch, tree, m2l_list, direct_list
end

function FastLinearOperator(panels::AbstractPanels{<:Any,TF,<:Any,<:Any}, scheme; fmm_toggle=true, reuse_tree=true, save_residual=false, expansion_order=4, leaf_size=SVector{1}(18), multipole_threshold=0.3) where TF

    self_induced = true
    switch, tree, m2l_list, direct_list = get_fmm_objects(panels, expansion_order, leaf_size, multipole_threshold, self_induced, reuse_tree)
    FastMultipole.unsort!((panels,), tree)

    return FastLinearOperator{typeof(panels), typeof(tree), typeof(direct_list), typeof(switch), scheme}(panels, fmm_toggle, reuse_tree, save_residual, tree, m2l_list, direct_list, switch, expansion_order, leaf_size, multipole_threshold)
end

function (flo::FastLinearOperator{<:AbstractPanels{K,<:Any,<:Any,<:Any}, <:Any, <:Any, <:Any, Scheme{DirectNeumann, FlowTangency}} where {K<:AbstractKernel})(C, B, α, β; fmm_args...)
    # unpack operator
    panels = flo.panels
    fmm_toggle = flo.fmm_toggle
    reuse_tree = flo.reuse_tree
    tree = flo.tree
    m2l_list = flo.m2l_list
    direct_list = flo.direct_list
    switch = flo.derivatives_switch
    expansion_order = flo.expansion_order
    leaf_size = flo.leaf_size
    multipole_threshold = flo.multipole_threshold

    # update panels with provided strengths
    for i in eachindex(panels.panels)
        # (; vertices, control_point, normal, strength, radius) = panels.panels[i]
        vertices = panels.panels[i].vertices
        control_point = panels.panels[i].control_point
        normal = panels.panels[i].normal
        old_strength = panels.panels[i].strength
        radius = panels.panels[i].radius
        strength = get_strength(B, i, old_strength)

        panels.panels[i] = Panel(vertices, control_point, normal, strength, radius)
        panels.strengths[i] = strength
    end

    # reset velocity
    @assert iszero(β) "FMM accelerated matrix-vector product not defined for nonzero beta"
    reset_potential_velocity!(panels)

    # solve N-body problem
    if fmm_toggle
        if reuse_tree && false
            FastMultipole.resort!((panels,), tree)
            fmm!((panels,), tree, (panels,), tree, m2l_list, direct_list, (switch,))
        else
            fmm!(panels; velocity_gradient=false, expansion_order, leaf_size_source=leaf_size, multipole_threshold, fmm_args...)
        end
    else
        direct!(panels; velocity_gradient=false)
    end

    # save normal component
    for i in eachindex(panels.panels)
        normal = panels.panels[i].normal
        velocity = panels.velocity[i]
        C[i] = dot(normal, velocity)
    end
end

function MatrixFreeSolver(panels::AbstractPanels{K,TF,<:Any,<:Any}, scheme,
                            krylov_solver::Type{<:Krylov.KrylovSolver}=Krylov.GmresSolver;
                            fmm_toggle=true, reuse_tree=true,
                            expansion_order=4, leaf_size=SVector{1}(18),
                            multipole_threshold=0.3
                            ) where {TF, K<:AbstractKernel{1}}

    # define fast linear operator functor (for use with FMM) (avoids closure)
    flo = FastLinearOperator(panels, scheme; fmm_toggle, reuse_tree, expansion_order, leaf_size, multipole_threshold)

    # define linear operator object for use with Krylov.jl
    n_panels = length(panels.panels)
    A = flo

    # construct ::KrylovSolver
    right_hand_side = zeros(TF,n_panels)
    solver = krylov_solver(A, right_hand_side)

    # construct solver
    return IterativeSolver(solver, A, right_hand_side, scheme)
end

function MatrixFreeSolver(panels::AbstractPanels{K,TF,<:Any,<:Any}, scheme,
                            krylov_solver::Type{<:Krylov.KrylovSolver}=Krylov.GmresSolver;
                            fmm_toggle=true, reuse_tree=true,
                            expansion_order=4, leaf_size=SVector{1}(18),
                            multipole_threshold=0.3
                            ) where {K<:AbstractKernel{2}, TF}

    # define fast linear operator functor (for use with FMM) (avoids closure)
    flo = FastLinearOperator(panels, scheme; fmm_toggle, reuse_tree, expansion_order, leaf_size, multipole_threshold)

    # define linear operator object for use with Krylov.jl
    n_panels = length(panels.panels)
    A = flo
    # A = LinearOperators.LinearOperator(TF, n_panels, n_panels, false, false, flo)

    # construct ::KrylovSolver
    right_hand_side = zeros(TF,n_panels)
    solver = krylov_solver(A, right_hand_side)

    # construct solver
    return IterativeSolver(solver, A, right_hand_side, scheme)
end

function MatrixFreeSolver_benchmark(panels::AbstractPanels{K,TF,<:Any,<:Any},
                                    scheme,
                                    krylov_solver::Type{<:Krylov.KrylovSolver}=Krylov.GmresSolver;
                                    fmm_toggle=true, reuse_tree=true,
                                    expansion_order=4, leaf_size=SVector{1}(18),
                                    multipole_threshold=0.3
                                    ) where {TF, K<:AbstractKernel{1}}

    t_alloc = @elapsed begin
        # define fast linear operator functor (for use with FMM) (avoids closure)
        flo = FastLinearOperator(panels, scheme; fmm_toggle, reuse_tree, expansion_order, leaf_size, multipole_threshold)

        # define linear operator object for use with Krylov.jl
        n_panels = length(panels.panels)
        # A = LinearOperators.LinearOperator(TF, n_panels, n_panels, false, false, flo)
        A = flo

        # construct ::KrylovSolver
        right_hand_side = zeros(TF,n_panels)
        solver = krylov_solver(A, right_hand_side)

    end

    # construct solver
    return IterativeSolver(solver, A, right_hand_side, scheme), t_alloc
end

function (solver::IterativeSolver{<:Any,<:Krylov.GmresSolver,<:Any,<:Any})(A, b; solver_kwargs...)
    Krylov.gmres!(solver.solver, A, b; solver_kwargs...)
end

function solve!(panels::AbstractPanels{K,<:Any,<:Any,<:Any},
                solver::IterativeSolver{<:Any,<:Any,<:Any,scheme},
                dt=0.0;
                verbose=true,
                tolerance=1e-6, max_iterations=100,
                solver_kwargs...
                ) where {scheme, K<:AbstractKernel{1}}
    # unpack
    # (; influence_matrix, right_hand_side, strengths) = solver
    A = solver.A
    right_hand_side = solver.right_hand_side

    # apply freestream/panel velocity
    update_right_hand_side!(right_hand_side, panels, scheme)

    # solver for strengths
    solver_output = solver(A, right_hand_side; atol=tolerance, itmax=max_iterations, solver_kwargs...)
    #x, stats = Krylov.gmres(influence_matrix, right_hand_side; atol=tolerance, itmax=max_iterations, solver_kwargs...)
    #=
                   memory=20, M=I, N=I, ldiv::Bool=false,
                   restart::Bool=false, reorthogonalization::Bool=false,
                   atol::T=√eps(T), rtol::T=√eps(T), itmax::Int=0,
                   timemax::Float64=Inf, verbose::Int=0, history::Bool=false,
                   callback=solver->false, iostream::IO=kstdout)
    =#
    verbose && println("Finished GMRES: \n\tstatus = $(solver.solver.stats.status)\n\tniter = $(solver.solver.stats.niter)\n\tinconsistent = $(solver.solver.stats.inconsistent)")

    strengths = solver.solver.x

    # update panels
    for i in 1:length(panels.panels)
        # (; vertices, control_point, normal, strength, radius) = panels.panels[i]
        vertices = panels.panels[i].vertices
        control_point = panels.panels[i].control_point
        normal = panels.panels[i].normal
        old_strength = panels.panels[i].strength
        radius = panels.panels[i].radius
        strength = get_strength(strengths, i, old_strength)

        panels.panels[i] = Panel(vertices, control_point, normal, strength, radius)
        panels.strengths[i] = strength
    end
end

"""
Assumes `freestream` has NOT been applied to `panels` yet, but any other induced velocity has.
"""
function solve!(panels::AbstractPanels{K,TF,<:Any,<:Any},
                solver::IterativeSolver{<:Any,<:Any,<:Any,scheme},
                freestream::AbstractVector,
                dt=0.0;
                verbose=true,
                tolerance=1e-6, max_iterations=100,
                solver_kwargs...
                ) where {scheme, K<:AbstractKernel{2}, TF}
    # unpack
    A = solver.A
    right_hand_side = solver.right_hand_side

    # solve for source panel strength
    for i_panel in eachindex(panels.panels)
        normal = panels.panels[i_panel].normal
        new_strength = dot(normal, freestream) * 4*pi

        # update panel strength
        panels.strengths[i_panel] = SVector{2}(new_strength, 0.0)
    end

    # update panels
    grid_2_panels_strength!(panels)

    # apply induced velocity
    if solver.A.fmm_toggle
        fmm!(panels; multipole_threshold=solver.A.multipole_threshold, leaf_size_source=solver.A.leaf_size, expansion_order=solver.A.expansion_order)
    else
        direct!(panels)
    end

    # apply freestream
    apply_freestream!(panels, freestream)

    # update right-hand-side
    update_right_hand_side!(right_hand_side, panels, scheme)

    # stash source panel strengths
    source_strengths = zeros(TF, length(panels.panels))
    for i in eachindex(source_strengths)
        source_strengths[i] = panels.strengths[i][1]
        panels.strengths[i] = SVector{2,TF}(0.0,0.0)
    end
    grid_2_panels_strength!(panels)

    # solver for strengths
    solver_output = solver(A, right_hand_side; atol=tolerance, itmax=max_iterations, solver_kwargs...)
    #x, stats = Krylov.gmres(influence_matrix, right_hand_side; atol=tolerance, itmax=max_iterations, solver_kwargs...)
    #=
                   memory=20, M=I, N=I, ldiv::Bool=false,
                   restart::Bool=false, reorthogonalization::Bool=false,
                   atol::T=√eps(T), rtol::T=√eps(T), itmax::Int=0,
                   timemax::Float64=Inf, verbose::Int=0, history::Bool=false,
                   callback=solver->false, iostream::IO=kstdout)
    =#
    verbose && println("Finished GMRES: \n\tstatus = $(solver.solver.stats.status)\n\tniter = $(solver.solver.stats.niter)\n\tinconsistent = $(solver.solver.stats.inconsistent)")

    dipole_strengths = solver.solver.x

    # update panels
    for i in 1:length(panels.panels)
        # (; vertices, control_point, normal, strength, radius) = panels.panels[i]
        vertices = panels.panels[i].vertices
        control_point = panels.panels[i].control_point
        normal = panels.panels[i].normal
        # strength = panels.panels[i].strength
        radius = panels.panels[i].radius
        strength = SVector{2,TF}(source_strengths[i], dipole_strengths[i])

        panels.panels[i] = Panel(vertices, control_point, normal, strength, radius)
        panels.strengths[i] = strength
    end
end

#------- Fast Gauss-Seidel Iterations -------#

struct FastGaussSeidel{TF,TP,TT,TDL,TDS,S} <: AbstractSolver
    panels::TP
    external_velocity::Vector{SVector{3,TF}}
    reuse_tree::Bool
    tree::TT
    m2l_list::Vector{SVector{2,Int32}}
    direct_list::TDL
    self_direct_list::Vector{SVector{2,Int32}}
    fmm_matrix_maps::Vector{Vector{Tuple{Int32,Int32,UnitRange{Int64}}}}
    derivatives_switch::TDS
    influence_matrices::Vector{LU{TF,Matrix{TF},Vector{Int32}}}
    strengths::Vector{TF}
    external_right_hand_side::Vector{TF}
    internal_right_hand_side::Vector{TF}
    expansion_order::Int
    multipole_threshold::Float64
    leaf_size::SVector{1,Int}
    convergence_history::Vector{Float64}
end

function FastGaussSeidel(panels::AbstractPanels{TK,TF,<:Any,<:Any}, scheme::Type{Scheme{DirectNeumann, FlowTangency}}; reuse_tree=false, expansion_order=5, multipole_threshold=0.5, leaf_size=SVector{1}(20)) where {TK,TF}
    external_velocity = zeros(SVector{3,TF}, length(panels.panels))

    # fmm objects
    println("fmm objects...")
    @time begin
    self_induced = false

    # note this function sorts the panels into the resulting tree
    switch, tree, m2l_list, direct_list = get_fmm_objects(panels, expansion_order, leaf_size, multipole_threshold, self_induced, reuse_tree)
    self_direct_list = Vector{SVector{2,Int32}}(undef,length(tree.leaf_index))
    for (i,i_leaf) in enumerate(tree.leaf_index)
        self_direct_list[i] = SVector{2,Int32}(Int32(i_leaf),Int32(i_leaf))
    end
    end # @time

    # compile fmm matrix maps
    println("fmm matrix maps...")
    @time begin

    if reuse_tree

        fmm_matrix_maps = Vector{Vector{Tuple{Int32,Int32,UnitRange{Int64}}}}(undef, length(tree.leaf_index))

        for (i_map,i_leaf) in enumerate(tree.leaf_index)
            fmm_matrix_map = Tuple{Int32,Int32,UnitRange{Int64}}[]

            for (i_target, j_source) in direct_list.direct_list
                if i_target == i_leaf # found a source
                    # determine row range of corresponding influence matrix
                    i_row_start = Int32(1)
                    n_rows = Int32(0)

                    # get row indices of the influence matrix corresponding to this target
                    for (ii_target, jj_source) in direct_list.direct_list
                        if jj_source == j_source # found this source
                            n_rows = Int32(3)*Int32(length(tree.branches[ii_target].bodies_index)) # hard-coded to only calculate velocity
                            if ii_target == i_target # no need to keep searching
                                break
                            else # otherwise, increment the row start to skip this target leaf
                                i_row_start += n_rows
                            end
                        end
                    end

                    # determine the index of the source leaf (to obtain its influence matrix)
                    j_matrix = Int32(-1)
                    for (source_matrix_index, source_leaf_index) in enumerate(tree.leaf_index)
                        if source_leaf_index == j_source
                            j_matrix = Int32(source_matrix_index)
                            break
                        end
                    end

                    # add fmm_matrix_map entry
                    push!(fmm_matrix_map, (j_source, j_matrix, i_row_start:i_row_start+n_rows-1))
                end
            end

            fmm_matrix_maps[i_map] = fmm_matrix_map
        end
    else
        fmm_matrix_maps = Vector{Vector{Tuple{Int32,Int32,UnitRange{Int64}}}}(undef, 0)
    end
    end # @time

    println("create influence matrices...")

    # create influence matrices to solve for each leaf's strengths
    panels_2_grid_strength!(panels) # store current strengths so they aren't lost
    set_unit_strength!(panels.panels)

    # influence_matrices = Vector{Matrix{TF}}(undef, length(tree.leaf_index))
    influence_matrices = Vector{LU{TF,Matrix{TF},Vector{Int32}}}(undef, length(tree.leaf_index))
    for (i_A, i_leaf) in enumerate(tree.leaf_index)
        leaf = tree.branches[i_leaf]
        n_panels = length(leaf.bodies_index)
        A = Matrix{TF}(undef, n_panels, n_panels)
        update_influence_matrix!(A, panels, scheme; panel_indices=leaf.bodies_index, set_unit_strength=true)
        influence_matrices[i_A] = lu!(A)
    end

    println("set up strengths...")

    # restore strengths
    grid_2_panels_strength!(panels)

    # unsort panels
    FastMultipole.unsort!((panels,), tree)

    # strengths
    strengths = zeros(TF,length(panels.panels))

    # RHS
    external_right_hand_side = zeros(TF,length(panels.panels))
    internal_right_hand_side = zeros(TF,length(panels.panels))

    # convergence history
    convergence_history = Float64[]

    println("return solver...")

    return FastGaussSeidel{TF,typeof(panels),typeof(tree),typeof(direct_list),typeof(switch),scheme}(panels, external_velocity, reuse_tree, tree, m2l_list, direct_list, self_direct_list, fmm_matrix_maps, switch, influence_matrices, strengths, external_right_hand_side, internal_right_hand_side, expansion_order, multipole_threshold, leaf_size, convergence_history)
end

function FastGaussSeidel_benchmark(panels::AbstractPanels{TK,TF,<:Any,<:Any}, scheme::Type{Scheme{DirectNeumann, FlowTangency}}; expansion_order=5, multipole_threshold=0.5, leaf_size=SVector{1}(20), reuse_tree=false) where {TK,TF}
    t_solve, t_fmm = 0.0, 0.0

    t_solve += @elapsed external_velocity = zeros(SVector{3,TF}, length(panels.panels))

    # fmm objects
    t_fmm += @elapsed begin
        switch, tree, m2l_list, direct_list = get_fmm_objects(panels, expansion_order, leaf_size, multipole_threshold, self_induced, reuse_tree)
        self_direct_list = Vector{SVector{2,Int32}}(undef,length(tree.leaf_index))
        for (i,i_leaf) in enumerate(tree.leaf_index)
            self_direct_list[i] = SVector{2,Int32}(Int32(i_leaf),Int32(i_leaf))
        end
    end

    # compile fmm matrix maps
    t_fmm += @elapsed begin

    fmm_matrix_maps = Vector{Vector{Tuple{Int32,Int32,UnitRange{Int64}}}}(undef, length(tree.leaf_index))

    if reuse_tree

        for (i_map,i_leaf) in enumerate(tree.leaf_index)
            fmm_matrix_map = Tuple{Int32,Int32,UnitRange{Int64}}[]

            for (i_target, j_source) in direct_list.direct_list
                if i_target == i_leaf # found a source
                    # determine row range of corresponding influence matrix
                    i_row_start = Int32(1)
                    n_rows = Int32(0)

                    # get row indices of the influence matrix corresponding to this target
                    for (ii_target, jj_source) in direct_list.direct_list
                        if jj_source == j_source # found this source
                            n_rows = Int32(3)*Int32(length(tree.branches[ii_target].bodies_index)) # hard-coded to only calculate velocity
                            if ii_target == i_target # no need to keep searching
                                break
                            else # otherwise, increment the row start to skip this target leaf
                                i_row_start += n_rows
                            end
                        end
                    end

                    # determine the index of the source leaf (to obtain its influence matrix)
                    j_matrix = Int32(-1)
                    for (source_matrix_index, source_leaf_index) in enumerate(tree.leaf_index)
                        if source_leaf_index == j_source
                            j_matrix = Int32(source_matrix_index)
                            break
                        end
                    end

                    # add fmm_matrix_map entry
                    push!(fmm_matrix_map, (j_source, j_matrix, i_row_start:i_row_start+n_rows-1))
                end
            end

            fmm_matrix_maps[i_map] = fmm_matrix_map
        end

    end

    end # @elapsed

    # create influence matrices to solve for each leaf's strengths
    t_solve += @elapsed panels_2_grid_strength!(panels) # store current strengths so they aren't lost
    t_solve += @elapsed set_unit_strength!(panels.panels)
    # influence_matrices = Vector{Matrix{TF}}(undef, length(tree.leaf_index))
    t_solve += @elapsed influence_matrices = Vector{LU{TF,Matrix{TF},Vector{Int32}}}(undef, length(tree.leaf_index))
    for (i_A, i_leaf) in enumerate(tree.leaf_index)
        leaf = tree.branches[i_leaf]
        n_panels = length(leaf.bodies_index)
        t_solve += @elapsed A = Matrix{TF}(undef, n_panels, n_panels)
        t_solve += @elapsed update_influence_matrix!(A, panels, scheme; panel_indices=leaf.bodies_index, set_unit_strength=true)
        influence_matrices[i_A] = lu!(A)
    end

    # restore strengths
    t_solve += @elapsed grid_2_panels_strength!(panels)

    # unsort panels
    t_fmm += @elapsed FastMultipole.unsort!((panels,), tree)

    # strengths
    t_solve += @elapsed strengths = zeros(TF,length(panels.panels))

    # RHS
    t_solve += @elapsed external_right_hand_side = zeros(TF,length(panels.panels))
    t_solve += @elapsed internal_right_hand_side = zeros(TF,length(panels.panels))

    # convergence history
    convergence_history = Float64[]

    return FastGaussSeidel{TF,typeof(panels),typeof(tree),typeof(direct_list),typeof(switch),scheme}(panels, external_velocity, reuse_tree, tree, m2l_list, direct_list, self_direct_list, fmm_matrix_maps, switch, influence_matrices, strengths, external_right_hand_side, internal_right_hand_side, expansion_order, multipole_threshold, leaf_size, convergence_history), t_solve, t_fmm
end

function solve!(panels::AbstractPanels{K,TF,<:Any,<:Any},
                solver::FastGaussSeidel{<:Any,<:Any,<:Any,<:Any,<:Any,scheme},
                dt=0.0;
                verbose=false,
                tolerance=1e-6, max_iterations=100,
                relaxation=1.4,
                error_style=:fixed_point,
                history=false
                ) where {TF, scheme, K<:AbstractKernel{1}}

    @assert panels === solver.panels "solver was created with a different AbstractPanels object"

    # unpack solver
    reuse_tree = solver.reuse_tree
    influence_matrices = solver.influence_matrices
    strengths = solver.strengths
    external_right_hand_side = solver.external_right_hand_side
    internal_right_hand_side = solver.internal_right_hand_side
    external_velocity = solver.external_velocity
    switch = solver.derivatives_switch
    expansion_order = solver.expansion_order
    multipole_threshold = solver.multipole_threshold
    leaf_size = solver.leaf_size
    convergence_history = solver.convergence_history
    history && sizehint!(convergence_history, max_iterations)

    # create/retrieve tree
    if reuse_tree
        tree = solver.tree
        FastMultipole.resort!((panels,), tree)
        m2l_list = solver.m2l_list
        direct_list = solver.direct_list
    else
        tree = FastMultipole.Tree(panels; expansion_order, leaf_size, shrink_recenter=true)
        farfield, nearfield, self_induced = true, true, false
        m2l_list, direct_list = FastMultipole.build_interaction_lists(tree.branches, tree.branches, leaf_size, tree.leaf_index, multipole_threshold, farfield, nearfield, self_induced, FastMultipole.UnequalSpheres(), expansion_order)
    end
    branches = tree.branches

    # save sorted external velocity
    for i in eachindex(panels.velocity)
        external_velocity[i] = panels.velocity[i]
    end

    # apply external velocity to external RHS, overwriting whatever was there previously
    update_right_hand_side!(external_right_hand_side, panels, scheme; reset=true)

    # convergence criterion
    ε = zero(TF)

    #--- begin Fast Gauss Seidel iterations ---#
    i_iter = 0 # want this available outside the loop
    verbose && println("\nBegin Fast Gauss-Seidel Iterations:")
    expansion_order_val = Val(expansion_order)
    for _ in 1:max_iterations
        i_iter += 1

        # copy external to internal RHS
        internal_right_hand_side .= external_right_hand_side

        # clear velocity
        reset_potential_velocity!(panels)

        # get farfield influence
        FastMultipole.fmm!(tree, panels, m2l_list, direct_list, switch, expansion_order_val; reset_tree=true, nearfield=false)#, unsort_bodies=false)

        # apply it to the RHS
        update_right_hand_side!(internal_right_hand_side, panels, scheme; reset=false)

        # clear velocity
        reset_potential_velocity!(panels)

        # iterate over leaves
        # for (A_solve, A_fmm, i_leaf) in zip(influence_matrices, direct_list.influence_matrices, tree.leaf_index)
        @assert length(tree.leaf_index) == length(influence_matrices)
        for (i_matrix, i_leaf) in enumerate(tree.leaf_index)

            # unpack
            A_solve = influence_matrices[i_matrix]
            leaf = branches[i_leaf]

            # get index
            panel_indices = leaf.bodies_index

            # apply all direct interactions to this leaf using the current strengths
            if solver.reuse_tree

                # unpack fmm objects
                strengths_fmm = direct_list.strengths
                influence_fmm = direct_list.influence

                # prepare just enough for this leaf
                n_influence = 3 * length(leaf.bodies_index) # hard-coded to only velocity
                this_influence = view(influence_fmm, 1:n_influence)
                this_influence .= zero(eltype(this_influence))

                # loop over all sources
                for (j_source, j_matrix, row_index) in solver.fmm_matrix_maps[i_matrix]

                    # update strengths
                    source_branch = branches[j_source]
                    n_strengths = FastMultipole.update_strengths!(strengths_fmm, source_branch, panels)
                    this_strength = view(strengths_fmm, 1:n_strengths)

                    # obtain influence
                    fmm_matrix = view(direct_list.influence_matrices[j_matrix], row_index, :)

                    # mul!(C, A, B, α, β) -> C; C = A B α + C β
                    mul!(this_influence, fmm_matrix, this_strength, 1, 1)
                end

                # apply influence to this leaf
                FastMultipole.update_influence!(panels, this_influence, 0, leaf.bodies_index, switch)

            else
                for (i_target, j_source) in direct_list
                    if i_target == i_leaf # found a source
                        FastMultipole._direct!(panels, leaf.bodies_index, switch, panels, branches[j_source].bodies_index)
                    end
                end
            end

            # update RHS
            update_right_hand_side!(internal_right_hand_side, panels, scheme; reset=false, panel_indices)

            # extract unknowns and RHS for this matrix
            s = view(strengths, panel_indices)
            rhs = view(internal_right_hand_side, panel_indices)

            # solve for strengths
            s .= A_solve \ rhs

            # update strengths
            vector_2_panels_strengths!(panels.panels, strengths; panel_indices, relaxation)

        end # iterate over leaves

        # update convergence criterion and update panels.strengths
        ε = zero(TF)
        inv_relaxation = 1/relaxation
        if error_style==:fixed_point # more performant
            for i in eachindex(panels.strengths)
                new_strength = panels.panels[i].strength
                old_strength = panels.strengths[i]
                ε = max(norm(new_strength-old_strength)*inv_relaxation, ε)
                panels.strengths[i] = new_strength
            end
        elseif error_style==:L2 # better error metric
            # have to compute the n-body problem again with the new strengths
            reset_potential_velocity!(panels)
            FastMultipole.fmm!(tree, panels, m2l_list, direct_list, switch; reset_tree=true)#, unsort_bodies=false)
            FastMultipole.fmm!(tree, panels, m2l_list, solver.self_direct_list, switch; reset_tree=false, upward_pass=false, horizontal_pass=false, downward_pass=false)#, unsort_bodies=false)
            # FastMultipole.direct!(panels)
            for i in eachindex(external_velocity)
                panels.velocity[i] += external_velocity[i]
            end

            # L2 norm error
            for (v,n) in zip(panels.velocity, panels.normals)
                ε += dot(v,n)^2
            end
            ε = sqrt(ε)
        else
            throw("error method $error_style not found")
        end
        verbose && println("\ti=$i_iter, epsilon = $(ε)")
        history && push!(convergence_history, Float64(ε))

        if ε <= tolerance
            break
        end

    end # iterations

    if ε <= tolerance
        verbose && println("Fast Gauss-Seidel Solver Successful:\n\titerations: $i_iter\n\terror: $ε")
    else
        verbose && println("Fast Guass-Seidel Solver Failed:\n\titerations: $i_iter\n\terror: $ε")
    end

    # restore external velocity (so it isn't lost)
    for i in eachindex(panels.velocity)
         panels.velocity[i] = external_velocity[i]
    end

    # unsort back to original order
    FastMultipole.unsort!((panels,), tree)

    # update grid strengths
    panels_2_grid_strength!(panels)

    return nothing
end


#------- other solvers -------#

struct SVDSolver
    epsilon::Float64
end

function (solver::SVDSolver)(A, b)

    # singular value decomp
    U, S, V = svd(A)

    # pseudo-inverse
    for i in eachindex(S)
        if abs(S[i]) > solver.epsilon
            S[i] = 1 / S[i]
        end
    end
    S_pseudo_inv = Diagonal(S)

    # solve
    x = V * S_pseudo_inv * U' * b
    return x
end

struct BlockGaussSeidel
    n_iter::Int
    n::Int
    rlx::Float64
end

function BlockGaussSeidel(; n_iter, n, rlx)
    return BlockGaussSeidel(n_iter, n, rlx)
end

function _bgs!(x0, A, b, n, rlx)
    N = length(x0)
    b2 = copy(b)
    for istart in 1:n:N
        iend = min(istart+n-1, N)
        bb = view(b2, istart:iend)
        # bb .= b[istart:iend]

        # subtract left influence
        if istart > 1
            dA = view(A, istart:iend, 1:istart-1)
            dx = view(x0, 1:istart-1)
            mul!(bb, dA, dx, -1.0, 1.0)
        end

        # subtract right influence
        if iend < N
            dA = view(A, istart:iend, iend+1:N)
            dx = view(x0, iend+1:N)
            mul!(bb, dA, dx, -1.0, 1.0)
        end

        # solve for x
        dA = view(A, istart:iend, istart:iend)
        dx = view(x0, istart:iend)
        # dx .*= (1-rlx)
        # dx .+= rlx .* (dA \ bb)
        dx .= dA \ bb
    end
end

function _bgs_mt!(x0, A, b, n, rlx)
    N = length(x0)
    b2 = copy(b)
    Threads.@threads for istart in 1:n:N
        iend = min(istart+n-1, N)
        bb = view(b2, istart:iend)
        # bb .= b[istart:iend]

        # subtract left influence
        if istart > 1
            dA = view(A, istart:iend, 1:istart-1)
            dx = view(x0, 1:istart-1)
            mul!(bb, dA, dx, -1.0, 1.0)
        end

        # subtract right influence
        if iend < N
            dA = view(A, istart:iend, iend+1:N)
            dx = view(x0, iend+1:N)
            mul!(bb, dA, dx, -1.0, 1.0)
        end

        # solve for x
        dA = view(A, istart:iend, istart:iend)
        dx = view(x0, istart:iend)
        # dx .*= (1-rlx)
        # dx .+= rlx .* (dA \ bb)
        dx .= dA \ bb
    end
end

function bgs_singlethread!(x0, A, b, n, rlx, n_iter)
    for _ in 1:n_iter
        _bgs!(x0, A, b, n, rlx)
    end
end

function bgs_multithread!(x0, A, b, n, rlx, n_iter)
    for _ in 1:n_iter
        _bgs_mt!(x0, A, b, n, rlx)
    end
end

function (bgs::BlockGaussSeidel)(A, b, x0=ones(length(b)))
    n_iter = bgs.n_iter
    n = bgs.n
    rlx = bgs.rlx
    if Threads.nthreads() > 1
        bgs_multithread!(x0, A, b, n, rlx, n_iter)
    else
        bgs_singlethread!(x0, A, b, n, rlx, n_iter)
    end
    return x0
end

struct GaussSeidel
    n_iter::Int
    atol::Float64
end

function _gs!(x0, A, b)
    N = length(x0)
    for i in 1:N
        x = b[i]
        for j in 1:i-1
            x -= A[i,j] * x0[j]
        end
        for j in i+1:N
            x -= A[i,j] * x0[j]
        end
        x0[i] = x / A[i,i]
    end
end

function _gs_mt!(x0, A, b, n_threads)
    N = length(x0)
    n_per_thread, rem = divrem(n_threads, N)
    rem > 0 && (n_per_thread += 1)
    i_starts = 1:n_per_thread:N

    Threads.@threads for i_start in i_starts
        for i in i_start:min(i_start-1+n_per_thread, N)
            x = b[i]
            for j in 1:i-1
                x -= A[i,j] * x0[j]
            end
            for j in i+1:N
                x -= A[i,j] * x0[j]
            end
            x0[i] = x / A[i,i]
        end
    end
end

function gs_singlethread!(x0, A, b, n_iter)
    for _ in 1:n_iter
        _gs!(x0, A, b)
    end
end

function gs_multithread!(x0, A, b, n_iter, n_threads)
    for i in 1:n_iter
        println("iteration $i")
        _gs_mt!(x0, A, b, n_threads)
    end
end

function (gs::GaussSeidel)(A, b, x0=ones(length(b)); n_threads = Threads.nthreads())
    if n_threads == 1
        gs_singlethread!(x0, A, b, gs.n_iter)
    else
        gs_multithread!(x0, A, b, gs.n_iter, n_threads)
    end
    return x0
end

struct LUSolver end

function (solver::LUSolver)(A, b, x=similar(b))
    # solve
    x .= A \ b
    return x
end

function solve_panels!(panels::AbstractPanels{TK,TF,NK,NS}, solver; A=nothing, watertight=false) where {TK,TF,NK,NS}
    # set up freestream
    AOA = 5.0 # deg
    magVinf = 10.0 # m/s
    Vinf = magVinf*SVector{3}(cos(AOA*pi/180), sin(AOA*pi/180), 0)

    # update panels
    for i in 1:length(panels.panels)
        panels.velocity[i] = Vinf
    end

    # set source strengths
    σs = zeros(TF, length(panels.panels))
    for i in 1:length(panels.panels)
        σ = -2 * dot(Vinf, panels.panels[i].normal)
        panels.strengths[i] = SVector{2}(σ, 0.0)
        σs[i] = σ
    end
    grid_2_panels_strength!(panels)

    # source-induced velocity
    expansion_order, leaf_size_source, multipole_threshold = 8, 18, 0.4
    fmm!(panels; expansion_order, leaf_size_source, multipole_threshold)

    # create rhs
    println("Creating RHS...")
    rhs = zeros(length(panels.panels)-watertight)
    for i in eachindex(rhs)
        v = panels.velocity[i]
        n = panels.panels[i].normal
        rhs[i] = -dot(v, n)
    end

    # initial guess
    # x0 = [panels.panels[i].strength[2] for i in eachindex(panels.panels)]
    # x0 = ones(length(panels.panels))

    # set panel strengths to unity
    for i in eachindex(panels.panels)
        (; vertices, control_point, normal, strength, radius) = panels.panels[i]
        panels.panels[i] = eltype(panels.panels)(vertices, control_point, normal, SVector{2}(0.0, 1.0), radius)
    end

    # create matrix
    println("Creating influence matrix...")
    if isnothing(A)
        A = zeros(length(panels.panels)-watertight, length(panels.panels)-watertight)
        Threads.@threads for i_source in 1:length(panels.panels)-watertight
            panel = panels.panels[i_source]
            for i_target in 1:length(panels.panels)-watertight
                x = panels.panels[i_target].control_point
                n = panels.panels[i_target].normal
                _, v, _ = induced(x, panel, TK(); sigma=1e-4)
                A[i_target, i_source] = dot(v, n)
            end
        end
    end

    println("solving...")
    μ = solver(A, rhs)
    # μ = A \ rhs
    println("done.")

    # update strengths
    for i in 1:length(panels.panels)-watertight
        old_strength = σs[i]
        panels.strengths[i] = SVector{2}(old_strength, μ[i])
    end
    if watertight
        panels.strengths[end] = SVector{2}(σs[end], 0.0)
    end
    grid_2_panels_strength!(panels)

    # check residual
    reset_potential_velocity!(panels)
    for i in eachindex(panels.panels)
        panels.velocity[i] = Vinf
    end

    fmm!(panels; expansion_order, leaf_size_source, multipole_threshold)

    resid = [dot(panels.panels[i].normal, panels.velocity[i]) for i in eachindex(panels.panels)]

    return resid, solver, A, rhs, μ
end



#=
#####
##### FMM-accelerated Jacobi iterations
#####

struct FastJacobi{NK,TF,NP,FO,S<:Scheme} <: AbstractSolver
    # containers
    # sort_index::Vector{Int}
    influence_matrices::Vector{Vector{TF}} # use vectors as they may be resized
    strength_history::Matrix{TF}
    previous_dt::MVector{NP,TF}
    internal_right_hand_side::Vector{TF}
    external_right_hand_side::Vector{TF}
    old_potential::Vector{TF}
    old_velocity::Vector{SVector{3,TF}}
    # fmm options
    expansion_order::Int
    n_per_branch::Int
    multipole_acceptance_criterion::Float64
    shrink_recenter::Bool
    ndivisions::Int
    # solver options
    max_iter::Int
    epsilon::Float64
end

# function apply_boundary_condition!(right_hand_side, panels::AbstractPanels{UniformSource(),<:Any,<:Any,<:Any}, scheme::FlowTangency)
#     for (i,(velocity,panel)) in enumerate(zip(panels.velocity, panels.panels))
#         right_hand_side[i] = -dot(velocity, panel.normal)
#     end
# end

function FastJacobi(panels::AbstractPanels{<:Any,TF,NK,<:Any}, scheme;
    fit_order=2, max_iter=200, epsilon=1e-4, n_previous_steps=3, expansion_order=5, n_per_branch=50, multipole_acceptance_criterion=0.4, ndivisions=10, shrink_recenter=true
) where {TF,NK}

    @assert fit_order < n_previous_steps "`fit_order` and `n_previous_steps` inconsistent: cannot create a best fit of order $(fit_order) with $n_previous_steps data points"

    # initialize memory
    influence_matrices = Vector{TF}[]
    n_unknowns = length(panels.panels)*NK
    strength_history = zeros(TF,n_unknowns,n_previous_steps)
    previous_dt = MVector{n_previous_steps}(t for t in range(0.0,stop=1.0,length=n_previous_steps))
    internal_right_hand_side = Vector{TF}(undef,n_unknowns)
    external_right_hand_side = Vector{TF}(undef,n_unknowns)
    old_potential = similar(vec(panels.potential))
    old_velocity = similar(vec(panels.velocity))

    # sort into octree
    tree = FastMultipole.Tree(panels; expansion_order, n_per_branch, shrink_recenter, ndivisions)

    # update inflence matrices
    # branch_index = get_branch_index(tree.levels_index, n_per_branch, max_n_per_matrix)
    # sort_index = get_sort_index(tree.branches, branch_index)
    update_influence_matrices!(influence_matrices, panels, tree)

    # resort panels
    FastMultipole.unsort!(panels, tree)

    return FastJacobi{NK,TF,n_previous_steps,fit_order,scheme}(influence_matrices, strength_history, previous_dt, internal_right_hand_side, external_right_hand_side, old_potential, old_velocity, expansion_order, n_per_branch, multipole_acceptance_criterion, shrink_recenter, ndivisions, max_iter, epsilon)
end

"""
`isnothing(dt)` indicates no previous solution exists, so Taylor series cannot be used to guess the next solution
"""
function solve!(panels::AbstractPanels{TK,<:Any,<:Any,<:Any}, solver::FastJacobi{NK,TF,<:Any,FO,S}, dt=nothing; update_influence_matrices=false, verbose=false) where {TK,NK,TF,FO,S}
    # unpack solver
    # (; influence_matrices, strength_history, previous_dt, external_right_hand_side, internal_right_hand_side, old_potential, old_velocity,
    #     expansion_order, n_per_branch, multipole_acceptance_criterion, shrink_recenter, ndivisions, max_iter, epsilon) = solver
    # sort_index = solver.sort_index
    influence_matrices = solver.influence_matrices
    strength_history = solver.strength_history
    previous_dt = solver.previous_dt
    external_right_hand_side = solver.external_right_hand_side
    internal_right_hand_side = solver.internal_right_hand_side
    old_potential = solver.old_potential
    old_velocity = solver.old_velocity
    expansion_order = solver.expansion_order
    n_per_branch = solver.n_per_branch
    multipole_acceptance_criterion = solver.multipole_acceptance_criterion
    shrink_recenter = solver.shrink_recenter
    ndivisions = solver.ndivisions
    max_iter = solver.max_iter
    epsilon = solver.epsilon

    # sort into octree
    tree = FastMultipole.Tree(panels; expansion_order, n_per_branch, shrink_recenter, ndivisions)

    # update inflence matrices
    update_influence_matrices && ( update_influence_matrices!(influence_matrices, panels, tree) )

    # save velocity
    for i in eachindex(old_potential)
        old_potential[i] = panels.potential[i]
        old_velocity[i] = panels.velocity[i]
    end

    # 1. initial guess of solved strengths based on previous strengths
    guess_strengths!(strength_history, previous_dt, Val(FO), dt) # strength_history[:,i] contains the solved strengths from i steps ago;
                                                                 # guess_strengths! stores its guess in the final column, as it will be discarded
    solved_strengths = view(strength_history,:,size(strength_history,2))
    vector_2_panels_strengths!(panels.panels, solved_strengths)

    # 2. form external RHS using wake- and freestream-induced velocity (noting that panels.velocity HAS been sorted)
    for (i,panel) in enumerate(panels.panels)
        velocity = panels.velocity[i]
        external_right_hand_side[i] = -dot(velocity, panel.normal)
    end

    # 3. block jacobi iterations
    error = epsilon * 2 # start with error greater than epsilon to ensure at least 1 iteration
    verbose && (println("Begin Fast Jacobi Iterations:"))
    i_iter = 1

    while error > epsilon && i_iter < max_iter
        verbose && (println("\titeration $i_iter: error = $error"))

        # get far-field influence by running the FMM
        for i in eachindex(panels.velocity)
            panels.velocity[i] = SVector{3,TF}(0,0,0)
            panels.potential[i] = zero(TF)
        end
        FastMultipole.fmm!(tree, panels; multipole_acceptance_criterion, self_induced=false, unsort_bodies=false)

        # move to right-hand side
        internal_right_hand_side .= zero(eltype(internal_right_hand_side))
        for (i,(panel,velocity)) in enumerate(zip(panels.panels, panels.velocity))
            internal_right_hand_side[i] = -dot(velocity, panel.normal)
        end

        # add external rhs
        internal_right_hand_side .+= external_right_hand_side

        # solve blocks
        for (branch,influence_matrix) in zip(view(tree.branches,tree.leaf_index), influence_matrices)
            this_index = branch.bodies_index
            strengths = view(solved_strengths, this_index)
            rhs = view(internal_right_hand_side, this_index)
            this_matrix = reshape(influence_matrix, length(rhs), length(rhs))
            strengths .= this_matrix \ rhs
        end

        # update strengths
        for i in eachindex(panels.panels)
            # (; vertices, control_point, normal, radius) = panels.panels[i]
            vertices = panels.panels[i].vertices
            control_point = panels.panels[i].control_point
            normal = panels.panels[i].normal
            radius = panels.panels[i].radius
            new_strength = SVector{1}(solved_strengths[i])
            panels.panels[i] = Panel(vertices, control_point, normal, new_strength, radius)
        end

        # evaluate error
        for i in eachindex(old_potential)
            panels.potential[i] = old_potential[i]
            panels.velocity[i] = old_velocity[i]
        end
        FastMultipole.fmm!(tree, panels; multipole_acceptance_criterion, self_induced=true, unsort_bodies=false)
        error = zero(error)
        for (panel, velocity) in zip(panels.panels, panels.velocity)
            error += dot(velocity, panel.normal)^2
        end
        error = sqrt(error/length(panels.panels))

        # tracking variables
        i_iter += 1
    end

    verbose && (println("\nFast Jacobi iterations finished: error = $error\n"))

    # update strength history
    # println("Update before")
    # @show strength_history
    update_strength_history!(strength_history, previous_dt, panels.panels, dt)
    # println("Update after")
    # @show strength_history

    FastMultipole.unsort!(panels, tree)

    # update panels
    for (i,panel) in enumerate(panels.panels)
        panels.strengths[i] = panel.strength
    end

    return nothing
end

function update_influence_matrices!(influence_matrices, panels::AbstractPanels{TK,TF,<:Any,<:Any}, tree::FastMultipole.Tree) where {TK,TF}
    # resize influence matrices
    leaf_index = tree.leaf_index
    n_leaves = length(leaf_index)

    # update size and number of existing matrices
    for i_matrix in eachindex(influence_matrices)
        leaf = tree.branches[leaf_index[i_matrix]]
        resize!(influence_matrices[i_matrix], length(leaf.bodies_index)^2)
    end
    for i_remaining in length(influence_matrices)+1:n_leaves
        leaf = tree.branches[leaf_index[i_remaining]]
        this_influence_matrix = zeros(TF,length(leaf.bodies_index)^2)
        push!(influence_matrices, this_influence_matrix)
    end
    resize!(influence_matrices, n_leaves)

    # update influence matrices
    set_unit_strength!(panels.panels)
    for (leaf, influence_matrix) in zip(view(tree.branches, leaf_index),influence_matrices)
        this_influence_matrix = reshape(influence_matrix, length(leaf.bodies_index), length(leaf.bodies_index))
        # for (i_source_panel,source_panel) in enumerate(( panels.panels[sort_index[i]] for i in leaf.bodies_index ))
        for (i_source_panel,source_panel) in enumerate(view(panels.panels,leaf.bodies_index))
            for (i_target_panel,target_panel) in enumerate(view(panels.panels,leaf.bodies_index))
            # for (i_target_panel,target_panel) in enumerate(( panels.panels[sort_index[i]] for i in leaf.bodies_index ))
                _, velocity, _ = induced(target_panel.control_point, source_panel, TK; toggle_potential=false, toggle_velocity=true, toggle_hessian=false)
                this_influence_matrix[i_target_panel,i_source_panel] = dot(velocity, target_panel.normal)
            end
        end
    end

    return nothing
end

# @inline function get_branch_index(levels_index, n_per_branch, max_n_per_matrix)
#     # choose level: n_per_branch * 8^i_level >= max_n_per_matrix
#     i_level = length(levels_index) - max(0,min(Int(ceil(log(8, max_n_per_matrix / n_per_branch))), length(levels_index)-1))
#     branch_index = levels_index[i_level]
#     return branch_index
# end

# function get_sort_index(branches, branch_index)
#     sort_index = zeros(Int,length(branches[1].bodies_index))
#     get_sort_index!(sort_index, branches, branch_index, 1)
#     return sort_index
# end

# function get_sort_index!(sort_index, branches, branch_index, i_start)
#     for branch in view(branches, branch_index)
#         if branch.n_branches == 0
#             sort_index[i_start:i_start+length(branch.bodies_index)-1] .= branch.bodies_index
#             i_start += length(branch.bodies_index)
#         else
#             get_sort_index!(sort_index, branches, branch.branch_index, i_start)
#         end
#     end
# end

function guess_strengths!(strength_history, previous_dt::MVector{NP,TF}, ::Val{FO}, dt) where {TF,NP,FO}
    @assert FO < NP

    # scale so the earliest time is unity
    earliest_t = previous_dt[end]
    previous_dt ./= earliest_t

    # create best fit coefficients matrix
    best_fit_coefficients = SMatrix{NP,FO+1,TF,NP*(FO+1)}(previous_dt[i]^j for i in 1:NP, j in 0:FO)
    At = transpose(best_fit_coefficients)
    AtAinv = inv(At * best_fit_coefficients)

    # create extrapolation polynomial
    dt_scaled = dt / earliest_t
    dt_predictor = transpose(SVector{FO+1,TF}(dt_scaled^p for p in 0:FO))

    # create predictions
    for i_unknown in axes(strength_history,1)
        # form RHS
        b = SVector{NP,TF}(strength_history[i_unknown,j] for j in 1:NP)
        Atb = At * b

        # solve for coefficients
        a = AtAinv * Atb

        # extrapolate strength and store
        strength_history[i_unknown,end] = dt_predictor * a
    end

    # unscale
    previous_dt .*= earliest_t

    return nothing
end

"""
Overload for a 0th order fit (guesses the next solution equal to the previous solution).
"""
function guess_strengths!(strength_history, previous_dt::MVector{NP,TF}, ::Val{0}, dt) where {NP,TF}
    @views strength_history[:,end] .= strength_history[:,1]
end

"""
If `isnothing(dt)`, do nothing.
"""
function guess_strengths!(strength_history, previous_dt::MVector{NP,TF}, ::Val{FO}, dt::Nothing) where {NP,TF,FO}
    @views strength_history[:,end] .= strength_history[:,1]
end

function guess_strengths!(strength_history, previous_dt::MVector{NP,TF}, ::Val{0}, dt::Nothing) where {NP,TF,FO}
    @views strength_history[:,end] .= strength_history[:,1]
end

"""
Initializes `solver.strength_history` to entirely consist of the current strengths of `panels.panels`.
"""
function initialize_strength_history!(solver::FastJacobi, panels::AbstractPanels)
    # (; strength_history, previous_dt) = solver
    strength_history = solver.strength_history
    previous_dt = solver.previous_dt
    update_strength_history!(strength_history, previous_dt, panels, nothing)
end

"""
Update `strength_history` using the current strengths of `panels.panels`, such that `strength_history[:,1]` contains the current solution, `strength_history[:,2]` contains the previous solution, and so on. Note that the final column of `strength_history` is discarded.
"""
function update_strength_history!(strength_history, previous_dt::MVector{NT,<:Any}, panels::AbstractVector{<:Panel}, dt) where {NT}
    # advance columns forward one step
    n_columns = size(strength_history,2)
    # println("Sherlock!")
    for i_column in n_columns-1:-1:1
        # @show strength_history
        @views strength_history[:,i_column+1] .= strength_history[:,i_column]
    end
    # @show strength_history
    # println("Holmes")

    # populate the first column with the current solution
    panels_2_vector_strengths!(view(strength_history,:,1), panels)

    # update previous_dt
    for i in 2:NT
        previous_dt[i] = previous_dt[i-1] - dt
    end
    @assert previous_dt[1] == 0.0 "previous_dt[1] should always equal 0.0 as it always references time relative to the first column of `strength_history`"
end

"""
Dispatch on `dt::Nothing` to populate the entire state history with the current strengths of `panels.panels`.
"""
function update_strength_history!(strength_history, previous_dt::MVector{NT,<:Any}, panels::AbstractVector{<:Panel}, dt::Nothing) where {NT}
    # populate the first column with the current solution
    panels_2_vector_strengths!(view(strength_history,:,1), panels)

    # copy to other columns
    n_columns = size(strength_history,2)
    for i_column in 2:n_columns
        @views strength_history[:,i_column] .= strength_history[:,i_column-1]
    end
end

#####
##### fast gauss-seidel iterations
#####
struct GaussSeidelMatrices{TF}
    matrix_A::Vector{TF}
    size_A::Int
    matrix_branch_list::Vector{Int}
    bodies_index_A::Vector{Int}
    rhs_branch_list::Vector{Tuple{Int,Int}}
end

struct FastGaussSeidel{NK,TF,NP,FO,S<:Scheme} <: AbstractSolver
    # containers
    # sort_index::Vector{Int}
    influence_matrices::Vector{GaussSeidelMatrices{TF}} # use vectors as they may be resized
    branch_checklist::Vector{Bool} # keep track of which direct interactions have already been accounted for
    strength_history::Matrix{TF}
    previous_dt::MVector{NP,TF}
    internal_right_hand_side::Vector{TF}
    external_right_hand_side::Vector{TF}
    old_potential::Vector{TF}
    old_velocity::Vector{SVector{3,TF}}
    # fmm options
    expansion_order::Int
    n_per_branch::Int
    multipole_acceptance_criterion::Float64
    shrink_recenter::Bool
    ndivisions::Int
    # solver options
    max_iter::Int
    epsilon::Float64
end

function FastGaussSeidel(panels::AbstractPanels{<:Any,TF,NK,<:Any}, scheme;
    fit_order=2, max_iter=200, epsilon=1e-4, n_previous_steps=3, expansion_order=5, n_per_branch=50, multipole_acceptance_criterion=0.4, ndivisions=10, shrink_recenter=true
) where {TF,NK}

    @assert fit_order < n_previous_steps "`fit_order` and `n_previous_steps` inconsistent: cannot create a best fit of order $(fit_order) with $n_previous_steps data points"

    # initialize memory
    influence_matrices = GaussSeidelMatrices{TF}[]
    n_unknowns = length(panels.panels)*NK
    strength_history = ones(TF,n_unknowns,n_previous_steps)
    previous_dt = MVector{n_previous_steps}(t for t in range(0.0,stop=1.0,length=n_previous_steps))
    internal_right_hand_side = Vector{TF}(undef,n_unknowns)
    external_right_hand_side = Vector{TF}(undef,n_unknowns)
    old_potential = similar(vec(panels.potential))
    old_velocity = similar(vec(panels.velocity))

    # sort into octree
    tree = FastMultipole.Tree(panels; expansion_order, n_per_branch, shrink_recenter, ndivisions)
    branch_checklist = fill(false, length(tree.branches)) # only need the leaves, but this makes it easier to navigate, and Bools are cheap
    m2l_list, direct_list = FastMultipole.build_interaction_lists(tree.branches, tree.branches, multipole_acceptance_criterion, true, true, true)
    @show m2l_list direct_list multipole_acceptance_criterion n_per_branch
    # update inflence matrices
    # branch_index = get_branch_index(tree.levels_index, n_per_branch, max_n_per_matrix)
    # sort_index = get_sort_index(tree.branches, branch_index)
    update_influence_matrices!(influence_matrices, panels, tree, multipole_acceptance_criterion, branch_checklist)

    # resort panels
    FastMultipole.unsort!(panels, tree)

    return FastGaussSeidel{NK,TF,n_previous_steps,fit_order,scheme}(influence_matrices, branch_checklist, strength_history, previous_dt, internal_right_hand_side, external_right_hand_side, old_potential, old_velocity, expansion_order, n_per_branch, multipole_acceptance_criterion, shrink_recenter, ndivisions, max_iter, epsilon), tree
end


"""
`isnothing(dt)` indicates no previous solution exists, so Taylor series cannot be used to guess the next solution
"""
function solve!(panels::AbstractPanels{TK,<:Any,<:Any,<:Any}, solver::FastGaussSeidel{NK,TF,<:Any,FO,S}, dt=nothing; update_influence_matrices=false, verbose=false) where {TK,NK,TF,FO,S}
    # unpack solver
    # (; influence_matrices, strength_history, previous_dt, external_right_hand_side, internal_right_hand_side, old_potential, old_velocity,
    #     expansion_order, n_per_branch, multipole_acceptance_criterion, shrink_recenter, ndivisions, max_iter, epsilon) = solver
    # sort_index = solver.sort_index
    influence_matrices = solver.influence_matrices
    branch_checklist = solver.branch_checklist
    strength_history = solver.strength_history
    previous_dt = solver.previous_dt
    external_right_hand_side = solver.external_right_hand_side
    internal_right_hand_side = solver.internal_right_hand_side
    old_potential = solver.old_potential
    old_velocity = solver.old_velocity
    expansion_order = solver.expansion_order
    n_per_branch = solver.n_per_branch
    multipole_acceptance_criterion = solver.multipole_acceptance_criterion
    shrink_recenter = solver.shrink_recenter
    ndivisions = solver.ndivisions
    max_iter = solver.max_iter
    epsilon = solver.epsilon

    # sort into octree
    tree = FastMultipole.Tree(panels; expansion_order, n_per_branch, shrink_recenter, ndivisions)
    m2l_list, direct_list = FastMultipole.build_interaction_lists(tree.branches, tree.branches, multipole_acceptance_criterion, true, true, true)

    # update inflence matrices
    update_influence_matrices && ( update_influence_matrices!(influence_matrices, panels, tree, multipole_acceptance_criterion, branch_checklist) )

    # save velocity/potential
    for i in eachindex(old_potential)
        old_potential[i] = panels.potential[i]
        old_velocity[i] = panels.velocity[i]
    end

    # 1. initial guess of solved strengths based on previous strengths
    guess_strengths!(strength_history, previous_dt, Val(FO), dt) # strength_history[:,i] contains the solved strengths from i steps ago;
                                                                 # guess_strengths! stores its guess in the final column, as it will be discarded
    solved_strengths = view(strength_history,:,size(strength_history,2))
    vector_2_panels_strengths!(panels.panels, solved_strengths)
    @show solved_strengths

    # 2. form external RHS using wake- and freestream-induced velocity (noting that panels.velocity HAS already been sorted)
    for (i,panel) in enumerate(panels.panels)
        velocity = panels.velocity[i]
        external_right_hand_side[i] = -dot(velocity, panel.normal)
    end

    # 3. block gauss-seidel iterations
    error = epsilon * 2 # start with error greater than epsilon to ensure at least 1 iteration
    verbose && (println("Begin Fast Gauss-Seidel Iterations:"))
    i_iter = 1

    # get far-field influence by running the FMM to prepare the first iteration
    # reset_potential_velocity!(panels)
    # FastMultipole.fmm!(tree, panels; multipole_acceptance_criterion, nearfield=false, unsort_bodies=false)

    while error > epsilon && i_iter < max_iter
        verbose && (println("\titeration $i_iter: error = $error"))

        # temporarily try this in the while loop
        reset_potential_velocity!(panels)
        # FastMultipole.fmm!(tree, panels; multipole_acceptance_criterion, nearfield=false, unsort_bodies=false)
        # direct for debugging
        for (i_branch_target,i_branch_source) in m2l_list
            this_target_index = tree.branches[i_branch_target].bodies_index
            this_source_index = tree.branches[i_branch_source].bodies_index
            FastMultipole._direct!(panels, this_target_index, panels, this_source_index)
        end

        # add farfield fmm influence and external influence to the RHS
        set_right_hand_side!(internal_right_hand_side, external_right_hand_side, panels)
        # i_iter == 1 && (internal_right_hand_side .=  [-0.4141721227865395,0.6907473543746981,-0.3980402370037941,-0.7738129213817879,-0.5033733566366838,0.9017307671895831,0.6907473543746981,-0.6240273847307353,-0.4872414708539383])

        # @show internal_right_hand_side external_right_hand_side

        # @show length(influence_matrices)
        # solve blocks
        for matrix_object in influence_matrices
            # unpack matrix object
            matrix_A = reshape(matrix_object.matrix_A, matrix_object.size_A, matrix_object.size_A)
            bodies_index_A = matrix_object.bodies_index_A
            rhs_branch_list = matrix_object.rhs_branch_list

            # add the influence of previously-solved panels to RHS
            @show rhs_branch_list
            for (i_target, i_source) in rhs_branch_list
                # extract branches
                target_index = tree.branches[i_target].bodies_index
                source_index = tree.branches[i_source].bodies_index

                # reset panels
                panels.potential[target_index] .= zero(TF)
                for i in target_index
                    panels.velocity[i] = SVector{3,TF}(0,0,0)
                end

                # perform direct calculation
                FastMultipole._direct!(panels, target_index, panels, source_index)

                # add to RHS
                for i_target in target_index
                    velocity = panels.velocity[i_target]
                    normal = panels.panels[i_target].normal
                    internal_right_hand_side[i_target] -= dot(velocity, normal)
                end
            end

            # solve matrix equation
            strengths = view(solved_strengths, bodies_index_A)
            rhs = view(internal_right_hand_side, bodies_index_A)
            @show rhs
            strengths .= matrix_A \ rhs
            @show strengths
            # update panels
            vector_2_panels_strengths!(view(panels.panels,bodies_index_A), strengths)
        end

        #####
        ##### evaluate error
        #####

        # # update the farfield influence of the panels
        # reset_potential_velocity!(panels)
        # FastMultipole.fmm!(tree, panels; multipole_acceptance_criterion, nearfield=false, unsort_bodies=false)

        # # add farfield fmm influence and external influence to the RHS, which is now our residual
        # set_right_hand_side!(internal_right_hand_side, external_right_hand_side, panels)

        # # add all direct calculations to the rhs, leaving nothing on the lhs (making it the residual)
        # for (i_matrix, matrix_object) in enumerate(influence_matrices)
        #     # unpack matrix object
        #     matrix_A = reshape(matrix_object.matrix_A, matrix_object.size_A, matrix_object.size_A)
        #     bodies_index_A = matrix_object.bodies_index_A
        #     matrix_B = reshape(matrix_object.matrix_B, matrix_object.size_B)
        #     bodies_index_B_source = matrix_object.bodies_index_B_source

        #     # add the influence of A panels to RHS
        #     strengths = view(solved_strengths, bodies_index_A)
        #     rhs = view(internal_right_hand_side, bodies_index_A)
        #     # println("Before")
        #     # @show rhs matrix_A strengths internal_right_hand_side
        #     mul!(rhs, matrix_A, strengths, -1, 1)
        #     # println("After")
        #     # @show rhs matrix_A strengths internal_right_hand_side
        #     # println()
        #     # add the influence of B panels to RHS
        #     if length(bodies_index_B_source) > 0
        #         strengths = view(solved_strengths, bodies_index_B_source)
        #         rhs = view(internal_right_hand_side, tree.branches[tree.leaf_index[i_matrix]].bodies_index)
        #         # println("B Before")
        #         # @show rhs matrix_B strengths
        #         mul!(rhs, matrix_B, strengths, -1, 1)
        #         # println("B after")
        #         # @show rhs matrix_B strengths
        #         # println()
        #     end
        # end

        # for debugging purposes, evaluate the entire FMM
        reset_potential_velocity!(panels)
        # FastMultipole.fmm!(tree, panels; multipole_acceptance_criterion, nearfield=true, unsort_bodies=false)
        FastMultipole.direct!(panels)
        set_right_hand_side!(internal_right_hand_side, external_right_hand_side, panels)

        # Linf norm error of the residual vector
        error = zero(error)
        for i in eachindex(internal_right_hand_side)
            # error += internal_right_hand_side[i]^2
            error = max(abs(internal_right_hand_side[i]),error)
        end
        # error = sqrt(error/length(panels.panels))

        # tracking variables
        i_iter += 1
    end

    verbose && (println("\nFast Gauss-Seidel iterations finished: error = $error\n"))

    # update strength history
    update_strength_history!(strength_history, previous_dt, panels.panels, dt)

    FastMultipole.unsort!(panels, tree)

    # update panels
    for (i,panel) in enumerate(panels.panels)
        panels.strengths[i] = panel.strength
    end

    # restore velocity
    for i in eachindex(old_potential)
        panels.potential[i] += old_potential[i]
        panels.velocity[i] += old_velocity[i]
    end

    return nothing
end

function set_right_hand_side!(internal_right_hand_side, external_right_hand_side, panels)
    # set equal to external rhs
    internal_right_hand_side .= external_right_hand_side

    # move current panel velocity to right-hand side
    for (i,(panel,velocity)) in enumerate(zip(panels.panels, panels.velocity))
        internal_right_hand_side[i] -= dot(velocity, panel.normal)
    end
end

function update_influence_matrices!(influence_matrices::Vector{<:GaussSeidelMatrices}, panels::AbstractPanels{TK,TF,<:Any,<:Any}, tree::FastMultipole.Tree, multipole_acceptance_criterion, branch_checklist) where {TK,TF}
    # obtain direct interaction list
    farfield, nearfield, self_induced = true, true, true
    m2l_list, direct_list = FastMultipole.build_interaction_lists(tree.branches, tree.branches, multipole_acceptance_criterion, farfield, nearfield, self_induced)
    # @show m2l_list direct_list multipole_acceptance_criterion tree.levels_index tree.leaf_index

    # unpack leaf index
    leaf_index = tree.leaf_index
    n_leaves = length(leaf_index)
    @show leaf_index tree.levels_index m2l_list direct_list multipole_acceptance_criterion

    # reset direct checklist
    branch_checklist .= false

    resize!(influence_matrices,0)

    @assert length(influence_matrices) == 0 "should start with zero matrices"

    # allocate influence matrices
    for (i_leaf,i_branch) in enumerate(leaf_index)
        if !branch_checklist[i_branch] # this leaf hasn't been included in an influence matrix yet
            matrix_branch_list = Int[]
            bodies_index_A = Int[]
            rhs_branch_list = Tuple{Int,Int}[]

            # find which branches in the interaction list haven't been included in a solve yet
            for (i_target,i_source) in direct_list
                if i_target == i_branch # in the interaction list
                    if !branch_checklist[i_source] # this source hasn't already been included in an influence matrix
                        println("filling matrix_branch_list:")
                        @show i_target i_source
                        println()
                        push!(matrix_branch_list, i_source)
                        branch_checklist[i_source] = true # mark it so it isn't repeated
                    end
                end
            end

            # which branches should be added to the RHS
            for i_branch in matrix_branch_list # loop over branches included in this solve
                for (i_target, i_source) in direct_list # find all branches in the interaction list that aren't included in this solve
                    if i_target == i_branch && !(i_source in matrix_branch_list)
                        push!(rhs_branch_list, (i_target,i_source))
                    end
                end
            end

            # create Gauss-Seidel influence matrix
            n_influencers = 0
            for leaf in view(tree.branches,matrix_branch_list)
                n_influencers += length(leaf.bodies_index)
            end
            resize!(bodies_index_A, n_influencers)
            i_start = 1
            for leaf in view(tree.branches,matrix_branch_list)
                bodies_index_A[i_start:i_start+length(leaf.bodies_index)-1] .= leaf.bodies_index
                i_start += length(leaf.bodies_index)
            end
            matrix_A = zeros(TF,n_influencers^2)
            push!(influence_matrices, GaussSeidelMatrices(matrix_A, n_influencers, matrix_branch_list, bodies_index_A, rhs_branch_list))
        end
    end

    # update influence matrices
    set_unit_strength!(panels.panels)
    for matrix_object in influence_matrices

        # unpack matrix object
        matrix_A = reshape(matrix_object.matrix_A, matrix_object.size_A, matrix_object.size_A)
        bodies_index_A = matrix_object.bodies_index_A

        # update matrix A
        for (i_source_matrix, i_source_panel) in enumerate(bodies_index_A)
            source_panel = panels.panels[i_source_panel]
            for (i_target_matrix, i_target_panel) in enumerate(bodies_index_A)
                target_panel = panels.panels[i_target_panel]
                _, velocity, _ = induced(target_panel.control_point, source_panel, TK; toggle_potential=false, toggle_velocity=true, toggle_hessian=false)
                matrix_A[i_target_matrix, i_source_matrix] = dot(velocity, target_panel.normal)
            end
        end
    end

    return nothing
end
=#
