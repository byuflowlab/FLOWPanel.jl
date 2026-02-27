#=##############################################################################
# DESCRIPTION
    Definition of waked paneled body types (implementations of
    AbstractLiftingBody).

# AUTHORSHIP
  * Created by  : Ryan Anderson
  * Email       : Ry.M.Anderson@gmail.com
  * Date        : Jan 2026
  * License     : MIT License
=###############################################################################


################################################################################
# DOUBLET-WAKE BODY TYPE
################################################################################
"""
  `FreeWakeBody{E::AbstractElement, N}(grid::gt.GridTriangleSurface,
shedding::Matrix{Int})`

Lifting body that is solver using a combination of N panel elements and an unsteady panel wake. 
`grid` is the grid surface (paneled geometry).
`shedding[:, i]` contains the information of the i-th edge along which to shed
the wake, where `shedding[1, i]` is the linear index of the panel shedding the
 wake, and `shedding[2:3, i]` are the indices of the nodes in that panel that
 make the edge. Since the wake is typically shed at the edge between two panels,
`shedding[3, i]` is the index of the partner panel (use -1 if none) and
`shedding[4:5, i]` are the node indices in that panel that make the edge.
The user must ensure that both edges are coincident, and the strength of the
wake is equal to the difference between the strengths of both panels.

  **Properties**
  * `nnodes::Int`                       : Number of nodes
  * `ncells::Int`                       : Number of cells
  * `fields::Vector{String}`            : Available fields (solutions)
  * `Oaxis::Matrix`                     : Coordinate system of body w.r.t. global
  * `O::Vector`                         : Origin of body w.r.t. global
  * `ncellsTE::Int`                     : Number of cells along trailing edge
  * `nnodesTE::Int`                     : Number of nodes along trailing edge

"""
mutable struct FreeWakeBody{E, N, TF} <: AbstractLiftingBody{E, N, TF}

    # User inputs
    grid::gt.GridTriangleSurface        # Paneled geometry
    shedding::Vector{Array{Int, 2}}             # Indicates edges along which to
                                        # shed the wake
    # Properties
    shedding_full::Matrix{Int}          # Map from panel index to shedding edge index (-1 if none)
    nnodes::Int                         # Number of nodes
    ncells::Int                         # Number of cells
    cells::Matrix{Int}                  # Cell connectivity (each column is a cell)
    nsheddings::Int                     # Number of shedding edges
    fields::Array{String, 1}            # Available fields (solutions)
    Oaxis::Array{TF,2}                  # Coordinate system of original grid
    O::Array{TF,1}                      # Position of CS of original grid

    # Internal variables
    strength::Array{TF, 2}              # strength[i,j] is the stength of the i-th panel with the j-th element type
    potential::Vector{TF}
    velocity::Array{TF, 2}              # velocity induced at control points
    CPoffset::Float64                   # Control point offset in normal direction
    kerneloffset::Float64               # Kernel offset to avoid singularities
    kernelcutoff::Float64               # Kernel cutoff to avoid singularities
    characteristiclength::Function      # Characteristic length of each panel
    watertight::Bool                    # Whether the body is watertight or not

    # wake properties
    wake_strength::Vector{Array{TF,2}}  # vector of structured grids of strength of each wake panel
    wake_nodes::Vector{Array{TF,2}}     # vector of structured grids of nodes bounding each wake panel

end

function FreeWakeBody{E, N, TF}(
                                grid, shedding::Vector{Matrix{Int,2}};
                                nnodes=grid.nnodes, ncells=grid.ncells,
                                cells=grid2cells(grid),
                                nsheddings=size(shedding,2),
                                fields=Array{String,1}(),
                                Oaxis = Matrix{TF}(I(3)), O = zeros(TF, 3),
                                strength=zeros(TF, grid.ncells, N),
                                potential=zeros(TF, grid.ncells),
                                velocity=zeros(TF, 3, grid.ncells),
                                CPoffset=1e-14,
                                kerneloffset=1e-8,
                                kernelcutoff=1e-14,
                                characteristiclength=characteristiclength_unitary,
                                check_mesh=true, watertight=true,
                                wake_strength=nothing,
                                wake_nodes=nothing,
                                n_wake_rows=ones(length(shedding))
                            ) where {E, N, TF}

    @assert _checkTE(grid, shedding) "Got invalid trailing edge"

    # Automated sanity checks for Meshes.jl mesh
    if check_mesh && typeof(grid.orggrid) <: gt.Meshes.Mesh

        mesh = grid.orggrid
        watertight = gt.isclosed(mesh)
        if !watertight
            @warn "Mesh is not watertight; results might be inaccurate"
        end

        # Check that control points lay outside the geometry (based on
        # winding number)
        normals = _calc_normals(grid)
        controlpoints = _calc_controlpoints(grid, normals;
                                            off=CPoffset,
                                            characteristiclength=characteristiclength)

        (minw, maxw) = calc_minmax_winding(mesh, controlpoints)

        if abs(minw) > eps()^0.75 || abs(maxw) >= eps()^0.75

            if watertight
                @warn "Found winding numbers other than 0, which might indicate"*
                " that control points are inside the geometry; flipping the"*
                " sign of `CPoffset` is recommended; (minw, maxw) ="*
                " $((minw, maxw))"
            else
                @warn "Found winding numbers other than 0, which might indicate"*
                " that control points are inside the geometry; however,"*
                " geometry is not watertight, so it might be ok."*
                " (minw, maxw) = $((minw, maxw))"
            end
        end

    end

    # generate full shedding map
    shedding_full = zeros(Int, 2, ncells)
    shedding_full .= -1
    
    # upper shedding edge
    for i in 1:nsheddings
        shedding_full[:, shedding[1, i]] .= view(shedding, 2:3, i)
    end

    # lower shedding edge
    for i in 1:nsheddings
        shedding_full[:, shedding[4, i]] .= view(shedding, 5:6, i)
    end

    # wake strength
    n_wake_cols = [size(shedding[i], 2) for i in eachindex(shedding)]
    if isnothing(wake_strength)
        wake_strength = [zeros(TF, nw, ns) for (nw, ns) in zip(n_wake_rows, n_wake_cols)]
    end

    # wake panel nodes
    if isnothing(wake_nodes)
        wake_nodes = [zeros(TF, 3, nw+1, ns+1) for (nw, ns) in zip(n_wake_rows, n_wake_cols)]
    end

    return RigidWakeBody{E, N, TF}(
                    grid, shedding, shedding_full,
                    nnodes, ncells, cells,
                    nsheddings,
                    Das, Dbs,
                    fields,
                    Oaxis, O,
                    strength,
                    potential,
                    velocity,
                    CPoffset,
                    kerneloffset,
                    kernelcutoff,
                    characteristiclength,
                    watertight,
                    wake_nodes,
                    wake_strength,
                )
end


function (RigidWakeBody{E})(grid, shedding; optargs...) where {E}
    return RigidWakeBody{E, _count(E), eltype(grid._nodes)}(grid, shedding; optargs...)
end

function (RigidWakeBody{E, N})(grid, shedding; optargs...) where {E, N}
    return RigidWakeBody{E, N, eltype(grid._nodes)}(grid, shedding; optargs...)
end

function (RigidWakeBody{E})(grid; optargs...) where {E}
    return RigidWakeBody{E}(grid, [zeros(Int, 6, 0)]; optargs...)
end

function (RigidWakeBody{E, N})(grid; optargs...) where {E, N}
    return RigidWakeBody{E, N}(grid, [zeros(Int, 6, 0)]; optargs...)
end

solved_field_name(::RigidWakeBody{ConstantSource, <:Any}) = "sigma"
solved_field_name(::RigidWakeBody{ConstantDoublet, <:Any}) = "mu"
solved_field_name(::RigidWakeBody{VortexRing, <:Any}) = "gamma"
solved_field_name(::RigidWakeBody{Union{ConstantSource, ConstantDoublet}, <:Any}) = "mu"

function save(body::RigidWakeBody, args...;
                out_wake::Bool=true, debug::Bool=false,
                wake_len::Number=1.0,
                wake_panel::Bool=false,
                wake_suffix="_wake",
                optargs...)


    str = ""
    str *= save_base(body, args...; debug=debug, optargs...)

    # Output the wake
    if out_wake
        str *= _savewake(body, args...; len=wake_len, panel=wake_panel,
                            optargs..., suffix=wake_suffix)
    end

    return str
end

#--- solve protocol ---#

function solve2!(self::FreeWakeBody{Union{ConstantSource,ConstantDoublet}, 2, TF}, Uinfs::Matrix{<:Real}, solver::BackslashDirichlet; backend=DirectBackend(), optargs...) where TF
    # get normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals; off=1e-10)
    CPs_inside = _calc_controlpoints(self, normals; off=-1e-10)

    # get source strengths
    ∂ϕ∂n = vec(sum(Uinfs .* .-normals, dims=1)) # known boundary condition
    self.strength[:, 1] .= ∂ϕ∂n
    self.strength[:, 2] .= 0.0

    # add source-induced potential to RHS
    rhs = solver.rhs
    rhs .= 0.0
    _phi!(self, CPs_inside, rhs, backend; optargs...)

    rhs .*= -1.0 # move to RHS
    
    # influence matrix for ϕ
    G = solver.G
    G .= 0.0

    method = 1

    if method == 1

        #--- method 1: solve for interior perturbation potential = 0 ---#

        # solve for doublet strengths such that the interior perturbation potential vanishes everywhere
        # then, the resulting velocity outside will be equal to the source strength, which satifsied flow tangency
        _G_phi!(self, ConstantDoublet, G, CPs_inside, backend; kerneloffset=1.0e-8, include_wake=true)
        
    elseif method == 2
        
        #--- method 2: (Morino formulation) solve for exterior potential = μ ---#
        
        _G_phi!(self, ConstantDoublet, G, CPs, backend; kerneloffset=1.0e-8, include_wake=true)

        # Fredholm second-kind equation
        for i in 1:size(G,1)
            G[i, i] += 1.0 # jump in potential equals doublet strength
        end
    end

    # Solve system of equations for the potential
    μ = G \ rhs

    # set doublet strengths
    self.strength[:, 2] .= μ
    
    # calculate inside potential
    # phis_inside = zeros(TF, self.ncells)
    # _phi!(self, CPs_inside, phis_inside, backend; include_wake=!EXTRA_FARFIELD[1], optargs...)

    # @show maximum(abs.(phis_inside))

    # calculate outside potential
    # phis_outside = zeros(TF, self.ncells)
    # _phi!(self, CPs, phis_outside, backend; include_wake=!EXTRA_FARFIELD[1], optargs...)

    # verify that source strength is equal to -∂ϕ/∂n on the boundary
    # us_outside = zeros(TF, 3, self.ncells)
    # _Uind!(self, CPs, us_outside, backend; include_wake=!EXTRA_FARFIELD[1], optargs...)
    # us_inside = zeros(TF, 3, self.ncells)
    # _Uind!(self, CPs_inside, us_inside, backend; include_wake=!EXTRA_FARFIELD[1], optargs...)

    # n∇ϕ = zeros(TF, 3, self.ncells)
    # _Uind!(self, CPs, n∇ϕ, backend; optargs...)
    # n∇ϕ = sum(n∇ϕ .* normals, dims=1)

    # freestream potential
    # phi_freestream = vec(sum(Uinfs .* CPs, dims=1))

    # Save solution
    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "Da", "vector", collect(eachcol(self.Das)), "system")
    add_field(self, "Db", "vector", collect(eachcol(self.Dbs)), "system")
    add_field(self, "sigma", "scalar", view(self.strength, :, 1), "cell")
    add_field(self, "mu", "scalar", view(self.strength, :, 2), "cell")
    # add_field(self, "phi_source", "scalar", -rhs, "cell")
    # add_field(self, "phi_inside", "scalar", phis_inside, "cell")
    # add_field(self, "phi_outside", "scalar", phis_outside, "cell")
    # add_field(self, "us_outside", "scalar", vec(sum(us_outside .* normals, dims=1)), "cell")
    # add_field(self, "us_inside", "scalar", vec(sum(us_inside .* normals, dims=1)), "cell")
    # add_field(self, "delta_u_normal", "scalar", vec(sum((us_outside - us_inside) .* normals, dims=1)), "cell")
    # add_field(self, "phi_freestream", "scalar", phi_freestream, "cell")
    # add_field(self, "ngradphi", "scalar", vec(n∇ϕ), "cell")
    add_field(self, "normals", "vector", collect(eachcol(normals)), "cell")
end