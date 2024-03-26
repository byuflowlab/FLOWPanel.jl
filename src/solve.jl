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
##### LU Decomposition
#####

struct LUDecomposition{TF,S<:Scheme} <: AbstractSolver
    influence_matrix::Array{TF,2}
    right_hand_side::Vector{TF}
    strengths::Vector{TF}
end

LUDecomposition(influence_matrix, right_hand_side, strengths, scheme) =
LUDecomposition{eltype(influence_matrix), scheme}(influence_matrix, right_hand_side, strengths)

function update_influence_matrix!(influence_matrix, panels::AbstractPanels{ConstantSource(),<:Any,<:Any,<:Any}, ::Type{Scheme{DirectNeumann, FlowTangency}})
    # check matrix size
    @assert size(influence_matrix,1) == size(influence_matrix,2) "influence matrix should be square"
    @assert size(influence_matrix,1) == length(panels.panels) "influence matrix size $(size(influence_matrix,1)) inconsistent with number of panels $(length(panels.panels))"

    # update influence matrix
    for (i_source,source_panel) in enumerate(panels.panels)
        for (i_target,target_panel) in enumerate(panels.panels)
            _, v, _ = induced(target_panel.control_point, source_panel, ConstantSource(); toggle_potential=false, toggle_velocity=true, toggle_hessian=false)
            influence_matrix[i_target, i_source] = dot(v, target_panel.normal)
        end
    end
end

function update_influence_matrix!(influence_matrix, panels::AbstractPanels{ConstantNormalDoublet(),<:Any,<:Any,<:Any}, ::Type{Scheme{DirectNeumann, FlowTangency}})
    # check matrix size
    @assert size(influence_matrix,1) == size(influence_matrix,2) "influence matrix should be square"
    @assert size(influence_matrix,1) == length(panels.panels) "influence matrix size $(size(influence_matrix,1)) inconsistent with number of panels $(length(panels.panels))"

    # update influence matrix
    for (i_source,source_panel) in enumerate(panels.panels)
        for (i_target,target_panel) in enumerate(panels.panels)
            _, v, _ = induced(target_panel.control_point, source_panel, ConstantNormalDoublet(); toggle_potential=false, toggle_velocity=true, toggle_hessian=false)
            influence_matrix[i_target, i_source] = dot(v, target_panel.normal)
            # i_source == i_target && (influence_matrix[i_target,i_source] /= 2)
        end
    end
end

function update_right_hand_side!(right_hand_side, panels::AbstractPanels, ::Type{Scheme{DirectNeumann, FlowTangency}}, freestream=SVector{3,Float64}(0,0,0))
    # check RHS size
    @assert length(right_hand_side) == length(panels.panels) "length of RHS $(length(right_hand_side)) inconsistent with number of panels $(length(panels.panels))"

    # update right hand side
    for (i,(panel,velocity)) in enumerate(zip(panels.panels, panels.velocity))
        right_hand_side[i] = -dot(freestream + velocity, panel.normal)
    end
end

function LUDecomposition(panels::AbstractPanels{ConstantSource(),TF,<:Any,<:Any}, scheme) where {TF}
    # preallocate memory
    n_panels = length(panels.panels)
    influence_matrix = zeros(TF,n_panels,n_panels)
    right_hand_side = zeros(TF,n_panels)

    # update influence matrix
    update_influence_matrix!(influence_matrix, panels, scheme)

    # get LU decomposition
    lu_decomposition = nothing#lu(influence_matrix)

    # initialize strengths
    strengths = zeros(TF,n_panels)

    return LUDecomposition(influence_matrix, right_hand_side, strengths, scheme)
end

function LUDecomposition(panels::AbstractPanels{ConstantNormalDoublet(),TF,<:Any,<:Any}, scheme) where {TF}
    # preallocate memory, assuming strength of the first panel is zero
    n_panels = length(panels.panels)
    influence_matrix = zeros(TF,n_panels,n_panels)
    right_hand_side = zeros(TF,n_panels)

    # update influence matrix
    update_influence_matrix!(influence_matrix, panels, scheme)

    # get LU decomposition
    lu_decomposition = nothing#lu(influence_matrix)

    # initialize strengths
    strengths = zeros(TF,n_panels)

    return LUDecomposition(influence_matrix, right_hand_side, strengths, scheme)
end

@inline function get_strength(strengths::Vector, i)
    return SVector{1}(strengths[i])
end

function solve!(panels::AbstractPanels{ConstantSource(),<:Any,<:Any,<:Any}, solver::LUDecomposition{<:Any,S}, dt=0.0) where S
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
        # strength = panels.panels[i].strength
        radius = panels.panels[i].radius
        strength = get_strength(strengths, i)

        panels.panels[i] = Panel(vertices, control_point, normal, strength, radius)
        panels.strengths[i] = strength
    end
end

function solve!(panels::AbstractPanels{ConstantNormalDoublet(),<:Any,<:Any,<:Any}, solver::LUDecomposition{<:Any,S}, dt=nothing) where S
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
        # strength = panels.panels[i].strength
        radius = panels.panels[i].radius
        strength = get_strength(strengths, i)
        panels.panels[i] = Panel(vertices, control_point, normal, strength, radius)
        panels.strengths[i] = strength
    end
end

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

# function apply_boundary_condition!(right_hand_side, panels::AbstractPanels{ConstantSource(),<:Any,<:Any,<:Any}, scheme::FlowTangency)
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
