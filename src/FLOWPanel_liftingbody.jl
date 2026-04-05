#=##############################################################################
# DESCRIPTION
    Definition of non-lifting paneled body types (implementations of
    AbstractLiftingBody).

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Sep 2018
  * License     : MIT License
=###############################################################################


################################################################################
# RIGID-WAKE BODY TYPE
################################################################################
"""
  `RigidWakeBody{E::AbstractElement, N}(grid::gt.GridTriangleSurface,
shedding::Matrix{Int})`

Lifting body that is solver using a combination of N panel elements and a steady
rigid wake. `grid` is the grid surface (paneled geometry).

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
mutable struct RigidWakeBody{E, N, TF, DBC} <: AbstractLiftingBody{E, N, TF, DBC}

    # User inputs
    nodes::Matrix{TF}                         # 3xnnodes matrix where nodes[:, i] is the position of the i-th node
    shedding::Vector{Array{Int, 2}}           # Indicates edges along which to
                                              # shed the wake
    # Properties
    vtk_cells::Vector{WriteVTK.MeshCell{WriteVTK.VTKCellTypes.VTKCellType, Vector{Int64}}}      # Vector of WriteVTK cells
    neighbor::Matrix{Int}                     # 3xncells matrix where neighbor[i, j] is the linear index of the cell neighboring the i-th edge of the j-th cell (or 0 if it's a boundary)
    shedding_full::Matrix{Int}                # Map from panel index to shedding edge index (0 if none)
    nnodes::Int                               # Number of nodes
    ncells::Int                               # Number of cells
    cells::Matrix{Int}                        # Cell connectivity (each column is a cell)
    nsheddings::Int                           # Number of shedding edges
    Das::Vector{Array{TF, 2}}                   # Unitary direction of rigid wake (vertex-based)
    Oaxis::Array{TF,2}                  # Coordinate system of original grid
    O::Array{TF,1}                      # Position of CS of original grid

    # Fields
    Cp::Vector{TF}
    F::Matrix{TF}

    # Internal variables
    strength::Array{TF, 2}              # strength[i,j] is the stength of the i-th panel with the j-th element type
    potential::Vector{TF}
    dphidt::Vector{TF}                   # Time derivative of perturbation potential (∂φ/∂t)
    velocity::Array{TF, 2}              # velocity induced at control points
    controlpoints::Matrix{TF}           # 3xncells control points
    normals::Matrix{TF}                 # 3xncells panel normals
    velocity_te::Vector{Matrix{TF}}     # velocity_te[i] is the velocity induced at the trailing edge of the i-th shedding edge
    CPoffset::Float64                         # Control point offset in normal direction
    kerneloffset::Float64                     # Kernel offset to avoid singularities
    kernelcutoff::Float64                     # Kernel cutoff to avoid singularities
    characteristiclength::Function            # Characteristic length of each panel
    watertight::Bool                         # Whether the body is watertight or not

    # wake properties
    semiinfinite_wake::Bool
end

function RigidWakeBody{E, N, TF, DBC}(
                                nodes::Matrix{TF}, cells::Matrix{Int}, shedding;
                                vtk_cells::Vector{<:WriteVTK.MeshCell}=[WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, cells[:, i]) for i in 1:size(cells, 2)],
                                neighbor::Matrix{Int}=zeros(Int, 3, size(cells, 2)),
                                nnodes=size(nodes, 2), ncells=size(cells, 2),
                                nsheddings=size(shedding,2),
                                Oaxis = Matrix{TF}(I(3)), O = zeros(TF, 3),
                                Cp=zeros(TF, size(cells, 2)),
                                F=zeros(TF, 3, size(cells, 2)),
                                Das::Vector{Matrix{TF}} = [zeros(TF, 3, size(s,2)+1) for s in shedding],
                                strength=zeros(TF, size(cells, 2), N),
                                potential=zeros(TF, size(cells, 2)),
                                dphidt=zeros(TF, size(cells, 2)),
                                velocity=zeros(TF, 3, size(cells, 2)),
                                controlpoints=zeros(TF, 3, ncells),
                                normals=zeros(TF, 3, ncells),
                                velocity_te=[zeros(TF, 3, size(s,2)+1) for s in shedding],
                                CPoffset=1e-14,
                                kerneloffset=1e-8,
                                kernelcutoff=1e-14,
                                characteristiclength=characteristiclength_unitary,
                                check_mesh=true, watertight=true,
                                semiinfinite_wake=true
                            ) where {E, N, TF, DBC}

    # for i_surf in eachindex(shedding)
    #     @assert _checkTE(nodes, cells, shedding[i_surf]) "Got invalid trailing edge"
    # end

    if check_mesh
        # we skip watertight checks for raw nodes and cells inputs, assumes true or false as flagged from the user.
    end

    # generate full shedding map
    # Rows: [1]=TE1 cell-local index, [2]=TE2 cell-local index,
    #        [3]=i_surf (which Das surface), [4]=edge index,
    #        [5]=Das column for TE1 vertex, [6]=Das column for TE2 vertex
    shedding_full = zeros(Int, 6, ncells)
    shedding_full .= -1

    # loop over shedding surfaces
    for i_surf in eachindex(shedding)

        # upper shedding edge
        this_shedding = shedding[i_surf]
        for i in axes(this_shedding, 2)
            idx = this_shedding[1,i]
            shedding_full[1:2, idx] .= view(this_shedding, 2:3, i)
            shedding_full[3, idx] = i_surf # which i_surf of Das to use
            shedding_full[4, idx] = i # which edge index of Das to use
            # Upper panel: TE1=nia uses Das[:,i+1], TE2=nib uses Das[:,i]
            shedding_full[5, idx] = i + 1  # Das column for TE1 (nia vertex)
            shedding_full[6, idx] = i      # Das column for TE2 (nib vertex)
        end

        # lower shedding edge
        for i in axes(this_shedding, 2)
            idx = this_shedding[4,i]
            if idx > 0
                shedding_full[1:2, idx] .= view(this_shedding, 5:6, i)
                shedding_full[3, idx] = i_surf
                shedding_full[4, idx] = i
                # Lower panel: TE edge is reversed, so TE1 is at the physical
                # nib position and TE2 is at the physical nia position
                shedding_full[5, idx] = i      # Das column for TE1 (physical nib)
                shedding_full[6, idx] = i + 1  # Das column for TE2 (physical nia)
            end
        end

    end

    return RigidWakeBody{E, N, TF, DBC}(
                    nodes, shedding, vtk_cells, neighbor, shedding_full,
                    nnodes, ncells, cells,
                    nsheddings,
                    Das,
                    Oaxis, O,
                    Cp, F,
                    strength,
                    potential,
                    dphidt,
                    velocity,
                    controlpoints,
                    normals,
                    velocity_te,
                    CPoffset,
                    kerneloffset,
                    kernelcutoff,
                    characteristiclength,
                    watertight,
                    semiinfinite_wake
                )
end

function RigidWakeBody{E, N, TF}(
                                nodes::Matrix{TF}, cells::Matrix{Int}, shedding;
                                DBC::Bool=true,
                                optargs...
                            ) where {E, N, TF}
    return RigidWakeBody{E, N, TF, DBC}(nodes, cells, shedding; optargs...)
end

_write_vtk_body_fields!(vtk, body::RigidWakeBody) =
    (vtk["dphidt", WriteVTK.VTKCellData()] = body.dphidt)

function RigidWakeBody{E, N, TF, DBC}(
                                grid::gt.GridTriangleSurface, shedding;
                                check_mesh=true, CPoffset=1e-14, characteristiclength=characteristiclength_unitary, watertight=true, optargs...
                            ) where {E, N, TF, DBC}
    
    nodes = grid._nodes
    cells = grid2cells(grid)
    
    # Extract neighbor info from grid
    neighbor = zeros(Int, 3, grid.ncells)
    
    ndivscellsc = Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells))
    linc = LinearIndices(ndivscellsc)
    
    for ci in 1:grid.ncells
        for ni in 1:3                   # Iterate over neighbors
            ncoor = gt.neighbor(grid, ni, ci; preserveEdge=true)
            if ncoor[1] != 0
                nlin = linc[ncoor...]
                neighbor[ni, ci] = nlin
            end
        end
    end

    vtk_cells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, cells[:, i]) for i in 1:grid.ncells]

    # Automated sanity checks for Meshes.jl mesh
    if check_mesh && typeof(grid.orggrid) <: gt.Meshes.Mesh

        mesh = grid.orggrid
        watertight = gt.isclosed(mesh)

        # # Check that topology is consistent with the solver
        # if watertight && E<:VortexRing && N==1
        #     @warn "Requested direct vortex ring solver on an closed mesh;"*
        #     " least-squares solver is recommended instead"*
        #     " (use `RigidWakeBody{VortexRing, 1}`)"
        # end

        # Check that control points lay outside the geometry (based on
        # winding number)
        normals = zeros(3, size(cells, 2))
        calc_normals!(nodes, cells, normals)

        controlpoints = zeros(3, size(cells, 2))
        calc_controlpoints!(nodes, cells, controlpoints, normals;
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

    return RigidWakeBody{E, N, TF, DBC}(
                    nodes, cells, shedding;
                    vtk_cells=vtk_cells, neighbor=neighbor, watertight=watertight, CPoffset=CPoffset, characteristiclength=characteristiclength, optargs...
                )
end

function RigidWakeBody{E, N, TF}(
                                grid::gt.GridTriangleSurface, shedding;
                                DBC::Bool=true,
                                optargs...
                            ) where {E, N, TF}
    return RigidWakeBody{E, N, TF, DBC}(grid, shedding; optargs...)
end

function RigidWakeBody{E, N, TF, DBC}(
                                mesh::VSPGeom.TriMesh, shedding;
                                optargs...
                            ) where {E, N, TF, DBC}
    nodes, cells = trimesh2cells(mesh)
    return RigidWakeBody{E, N, TF, DBC}(TF.(nodes), cells, shedding; optargs...)
end

function RigidWakeBody{E, N, TF}(
                                mesh::VSPGeom.TriMesh, shedding;
                                DBC::Bool=true,
                                optargs...
                            ) where {E, N, TF}
    return RigidWakeBody{E, N, TF, DBC}(mesh, shedding; optargs...)
end

function (RigidWakeBody{E})(grid::gt.GridTriangleSurface, shedding; DBC::Bool=true, optargs...) where {E}
    return RigidWakeBody{E, kernel_dim(E), eltype(grid._nodes), DBC}(grid, shedding; optargs...)
end

function (RigidWakeBody{E, N})(grid::gt.GridTriangleSurface, shedding; DBC::Bool=true, optargs...) where {E, N}
    return RigidWakeBody{E, N, eltype(grid._nodes), DBC}(grid, shedding; optargs...)
end

function (RigidWakeBody{E})(mesh::VSPGeom.TriMesh, shedding; DBC::Bool=true, optargs...) where {E}
    nodes, _ = trimesh2cells(mesh)
    return RigidWakeBody{E, kernel_dim(E), eltype(nodes), DBC}(mesh, shedding; optargs...)
end

function (RigidWakeBody{E, N})(mesh::VSPGeom.TriMesh, shedding; DBC::Bool=true, optargs...) where {E, N}
    nodes, _ = trimesh2cells(mesh)
    return RigidWakeBody{E, N, eltype(nodes), DBC}(mesh, shedding; optargs...)
end

function (RigidWakeBody{E})(grid::gt.GridTriangleSurface; optargs...) where {E}
    return RigidWakeBody{E}(grid, Vector{Array{Int, 2}}(); optargs...)
end

function (RigidWakeBody{E, N})(grid::gt.GridTriangleSurface; optargs...) where {E, N}
    return RigidWakeBody{E, N}(grid, Vector{Array{Int, 2}}(); optargs...)
end

function (RigidWakeBody{E})(mesh::VSPGeom.TriMesh; optargs...) where {E}
    return RigidWakeBody{E}(mesh, Vector{Array{Int, 2}}(); optargs...)
end

function (RigidWakeBody{E, N})(mesh::VSPGeom.TriMesh; optargs...) where {E, N}
    return RigidWakeBody{E, N}(mesh, Vector{Array{Int, 2}}(); optargs...)
end

function (RigidWakeBody{E})(nodes::Matrix{TF}, cells::Matrix{Int}, shedding; DBC::Bool=true, optargs...) where {E, TF}
    return RigidWakeBody{E, kernel_dim(E), TF, DBC}(nodes, cells, shedding; optargs...)
end

function (RigidWakeBody{E, N})(nodes::Matrix{TF}, cells::Matrix{Int}, shedding; DBC::Bool=true, optargs...) where {E, N, TF}
    return RigidWakeBody{E, N, TF, DBC}(nodes, cells, shedding; optargs...)
end

function (RigidWakeBody{E})(nodes::Matrix{TF}, cells::Matrix{Int}; optargs...) where {E, TF}
    return RigidWakeBody{E}(nodes, cells, Vector{Array{Int, 2}}(); optargs...)
end

function (RigidWakeBody{E, N})(nodes::Matrix{TF}, cells::Matrix{Int}; optargs...) where {E, N, TF}
    return RigidWakeBody{E, N}(nodes, cells, Vector{Array{Int, 2}}(); optargs...)
end

solved_field_name(::RigidWakeBody{ConstantSource, <:Any}) = "sigma"
solved_field_name(::RigidWakeBody{ConstantDoublet, <:Any}) = "mu"
solved_field_name(::RigidWakeBody{VortexRing, <:Any}) = "gamma"
solved_field_name(::RigidWakeBody{Union{ConstantSource, ConstantDoublet}, <:Any}) = "mu"
solved_field_name(::RigidWakeBody{Union{ConstantSource, VortexRing}, <:Any}) = "gamma"

wake_name(::RigidWakeBody{ConstantDoublet, <:Any}) = "mu"
wake_name(::RigidWakeBody{VortexRing, <:Any}) = "gamma"
wake_name(::RigidWakeBody{Union{ConstantSource, ConstantDoublet}, <:Any}) = "mu"
wake_name(::RigidWakeBody{Union{ConstantSource, VortexRing}, <:Any}) = "gamma"


################################################################################
# VORTEX RING SOLVER
################################################################################
function solve(self::RigidWakeBody{VortexRing, 1},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2};
                solver=solve_ludiv!, solver_optargs=(),
                optargs...
                ) where {T1, T2}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, size(self.Das[1], 2))
        error("Invalid Das;"*
              " expected size (3, $(size(self.Das[1], 2))), got $(size(Das))")
    end

    T = promote_type(T1, T2)

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # update Das
    @assert length(self.Das) == 1 "Matrix-valued Das input is only supported for a single shedding surface."
    self.Das[1] .= Das

    # Compute geometric matrix (left-hand-side influence matrix)
    G = zeros(T, self.ncells, self.ncells)
    _G_Uvortexring!(self, G, CPs, normals; optargs...)

    # Calculate boundary conditions (right-hand side of system of equations)
    RHS = calc_bc_noflowthrough(Uinfs, normals)

    # Solve system of equations
    Gamma = zeros(T, self.ncells)
    solver(Gamma, G, RHS; solver_optargs...)

    # Save solution
    self.strength[:, 1] .= Gamma

    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "Da", "vector", collect(eachcol(Das)), "system")
    add_field(self, "Gamma", "scalar", view(self.strength, :, 1), "cell")
end



function solve(self::RigidWakeBody{VortexRing, 2},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2};
                solver=solve_ludiv!, solver_optargs=(),
                elprescribe::AbstractArray{Tuple{Int, Float64}}=[(1, 0.0)],
                GPUArray=Array{promote_type(T1, T2)},
                optargs...
                ) where {T1, T2}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, size(self.Das[1], 2))
        error("Invalid Das;"*
              " expected size (3, $(size(self.Das[1], 2))), got $(size(Das))")
    end

    T = promote_type(T1, T2)

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix) and boundary
    # conditions (right-hand-side) converted into a least-squares problem
    println("BEGIN CONSTRUCTING G...")
    G, RHS = _G_U_RHS(self, Uinfs, CPs, normals, elprescribe;
                                                GPUArray=GPUArray,
                                                optargs...)
    println("FINISHED CONSTRUCTING G.")

    # Solve system of equations
    println("Solving system of equations...")
    Gamma = GPUArray(undef, self.ncells-length(elprescribe))
    solver(Gamma, G, RHS; solver_optargs...)
    println("Finished solving system of equations.")

    # Port solution back to CPU if solved in GPU
    if !(GPUArray <: Array)
        Gamma = Array{T}(Gamma)
    end

    # Save solution
    set_solution(self, nothing, Gamma, elprescribe, Uinfs)

end


function _G_Uvortexring!(self::RigidWakeBody,
                            G::Arr1, CPs::Arr2, normals::Arr3;
                            optargs...
                       ) where{ T1, Arr1<:AbstractArray{T1, 2},
                                T2, Arr2<:AbstractArray{T2, 2},
                                T3, Arr3<:AbstractArray{T3, 2}}

    N = self.ncells
    M = size(CPs, 2)

    Das = self.Das

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    elseif size(normals, 2)!=M
        error("normals matrix with invalid dimensions;"*
              " got $(size(normals)), expected (3, $M).")
    end


    # Build geometric matrix from panel contributions
    panels = 1:self.ncells
    chunks = collect(Iterators.partition(panels, max(length(panels) ÷ Threads.nthreads(), 3*Threads.nthreads())))

    Threads.@threads for chunk in chunks      # Distribute panel iteration among all CPU threads

        # for (pj, Gslice) in enumerate(eachcol(G))
        for pj in chunk                       # Iterate over panels

            panel = self.cells[:, pj]

            U_vortexring(
                              self.nodes,                        # All nodes
                              panel,                             # Indices of nodes that make this panel
                              1.0,                               # Unitary strength
                              CPs,                               # Targets
                              view(G, :, pj);                    # Velocity of j-th panel on every CP
                              # Gslice;
                              dot_with=normals,                  # Normal of every CP
                              offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                              cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                              optargs...
                             )

         end
    end


    # Add wake contributions
    sheddings = 1:self.nsheddings
    chunks = collect(Iterators.partition(sheddings, max(length(sheddings) ÷ Threads.nthreads(), 3*Threads.nthreads())))

    # for chunk in chunks        # Distribute wake panel iteration among all CPU threads
    Threads.@threads for chunk in chunks        # Distribute wake panel iteration among all CPU threads

        # Pre-allocate memory for panel calculation
        TE = zeros(Int, 2)

        # for (ei, (pi, nia, nib, pj, nja, njb)) in enumerate(eachcol(self.shedding))
        for ei in chunk                          # Iterate over wake-shedding panels

            pi, nia, nib, pj, nja, njb = view(self.shedding, :, ei)

            # Indicate nodes in the upper shedding edge
            TE[1] = self.cells[nia, pi]
            TE[2] = self.cells[nib, pi]
            da1, da2, da3 = Das[1, ei+1], Das[2, ei+1], Das[3, ei+1]
            db1, db2, db3 = Das[1, ei], Das[2, ei], Das[3, ei]

            U_semiinfinite_horseshoe(
                              self.nodes,                        # All nodes
                              TE,                                # Indices of nodes that make the shedding edge
                              da1, da2, da3,                     # Semi-infinite direction da
                              db1, db2, db3,                     # Semi-infinite direction db
                              1.0,                               # Unitary strength
                              CPs,                               # Targets
                              view(G, :, pi);                    # Velocity of upper wake panel on every CP
                              dot_with=normals,                  # Normal of every CP
                              offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                              cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                              optargs...
                             )

             if pj != -1
                 # Indicate nodes in the lower shedding edge
                 TE[1] = self.cells[nja, pj]
                 TE[2] = self.cells[njb, pj]

                 U_semiinfinite_horseshoe(
                                   self.nodes,                        # All nodes
                                   TE,                                # Indices of nodes that make the shedding edge
                                   db1, db2, db3,                     # Semi-infinite direction da (flipped in lower panel)
                                   da1, da2, da3,                     # Semi-infinite direction db
                                   1.0,                               # Unitary strength
                                   CPs,                               # Targets
                                   view(G, :, pj);                    # Velocity of lower wake panel on every CP
                                   dot_with=normals,                  # Normal of every CP
                                   offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                                   cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                                   optargs...
                                  )
             end
        end
    end
end

_get_Gdims(self::RigidWakeBody{VortexRing, N}) where {N} = (self.ncells, self.ncells)


################################################################################
# VORTEX RING + UNIFORM VORTEX SHEET SOLVER
################################################################################

function (RigidWakeBody{Union{VortexRing, UniformVortexSheet}})(grid, args...; optargs...)
    return RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 3}(grid, args...;
                                                    strength=zeros(grid.ncells, 3), optargs...)
end

function solve(self::RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 3},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2};
                solver=solve_ludiv!, solver_optargs=(),
                elprescribe_index::Int=1, elprescribe_value=0,
                weight_gammat=0, weight_gammao=1
                ) where {T1, T2}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings+1)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings+1)), got $(size(Das))")
    end

    T = promote_type(T1, T2)

    # store Das
    self.Das .= Das

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    # and boundary conditions (right-hand side of system of equations)
    G = zeros(T, self.ncells, self.ncells)
    RHS = zeros(T, self.ncells)

    _G_U_RHS!(self, G, RHS, Uinfs, CPs, normals,
                elprescribe_index, elprescribe_value,
                weight_gammat, weight_gammao)

    # Solve system of equations
    Gamma = zeros(T, self.ncells)
    solver(Gamma, G, RHS; solver_optargs...)

    # Save vortex ring circulations
    self.strength[:, 1] .= Gamma
    self.strength[elprescribe_index, 1] = elprescribe_value

    # Save vortex sheet strength
    gamma = Gamma[elprescribe_index]
    self.strength[:, 2] .= gamma*weight_gammat
    self.strength[1:2:end, 2] .*= -1
    self.strength[:, 3] .= gamma*weight_gammao
    self.strength[1:2:end, 3] .*= -1

    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "Da", "vector", collect(eachcol(Das)), "system")
    add_field(self, "Gamma", "scalar", view(self.strength, :, 1), "cell")

    tangents = _calc_tangents(self)
    obliques = _calc_obliques(self)
    aux = zip(eachcol(tangents), eachcol(obliques),
                view(self.strength, :, 2), view(self.strength, :, 3))
    gammas = [gammat*t + gammao*o for (t, o, gammat, gammao) in aux]
    add_field(self, "gamma", "vector", gammas, "cell")
end


function _G_U_RHS!(self::RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 3},
                    G, RHS, Uinfs, CPs, normals,
                    elprescribe_index::Int, elprescribe_value::Number,
                    weight_gammat::Number, weight_gammao::Number;
                    optargs...
                    )

    # Calculate normal velocity of freestream for boundary condition
    calc_bc_noflowthrough!(RHS, Uinfs, normals)

    # -------------- Influence of vortex rings -------------------------

    # Calculate influence of all vortex rings
    _G_Uvortexring!(self, G, CPs, normals; optargs...)

    # Move influence of prescribed vortex ring element to right-hand side
    for i in 1:length(RHS)
        RHS[i] -= elprescribe_value*G[i, elprescribe_index]
        G[i, elprescribe_index] = 0
    end

    # -------------- Influence of vortex sheet ------------------------

    # Influence of vortex sheet on each CP gets stored here
    Gslice = view(G, :, elprescribe_index)

    # Calculate influence of each panel on each CP
    for pj in 1:self.ncells         # Iterate over panels

        panel = self.cells[:, pj]

        s = pj%2==1 ? -1 : 1        # Alternate + and - strengths to get them all aligned

        U_constant_vortexsheet(
                          self.nodes,                        # All nodes
                          panel,                             # Indices of nodes that make this panel
                          s*weight_gammat,                   # Tangential strength
                          s*weight_gammao,                   # Oblique strength
                          CPs,                               # Targets
                          # view(G, :, pj);                  # Agglomerate velocity of j-th panel on every CP
                          Gslice;
                          dot_with=normals,                  # Normal of every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                          optargs...
                         )
    end

    return G, RHS
end

# function influence!(self::RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 3},
#                      target_body::AbstractBody, backend::DirectBackend;
#                      scalar_potential=false, velocity=false,
#                      velocity_gradient=false, optargs...)
#     # scalar_potential is a no-op for this body type
#     if velocity
#         _Uvortexring!(self, target_body.controlpoints, target_body.velocity, backend; stri=1, optargs...)
#         _Uconstantvortexsheet!(self, target_body.controlpoints, target_body.velocity; strti=2, stroi=3, optargs...)
#     end
#     return nothing
# end


function _Uconstantvortexsheet!(self::RigidWakeBody, targets, out;
                                                strti=2, stroi=3, optargs...)

    # Iterates over body panels
    for i in 1:self.ncells

        panel = self.cells[:, i]

        # Velocity of i-th panel on every target
        U_constant_vortexsheet(
                            self.nodes,                        # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, strti],           # Tangential strength
                            self.strength[i, stroi],           # Oblique strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                            optargs...
                         )
    end

end





################################################################################
# FASTMULTIPOLE BACKEND SUPPORT
################################################################################

FastMultipole.has_vector_potential(::AbstractBody{VortexRing, 1, <:Any}) = true

FastMultipole.has_vector_potential(::AbstractBody{Union{VortexRing, UniformVortexSheet}, 2, <:Any}) = true

has_semiinfinite_wake(self::RigidWakeBody) = self.semiinfinite_wake

function FastMultipole.body_to_multipole!(system::RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing}, 2, TF}, multipole_coefficients, buffer, center, bodies_index, harmonics, expansion_order) where TF
    # loop over bodies
    for i_body in bodies_index
        # relative body position
        # i1 = 5 + strength_dims(system)
        # x1 = SVector{3}(buffer[i1,i_body], buffer[i1+1,i_body], buffer[i1+2,i_body])
        # x2 = SVector{3}(buffer[i1+3,i_body], buffer[i1+4,i_body], buffer[i1+5,i_body])
        # x3 = SVector{3}(buffer[i1+6,i_body], buffer[i1+7,i_body], buffer[i1+8,i_body])
        x1 = FastMultipole.get_vertex(buffer, system, i_body, 1)
        x2 = FastMultipole.get_vertex(buffer, system, i_body, 2)
        x3 = FastMultipole.get_vertex(buffer, system, i_body, 3)
        x0 = x1 - center
        xu = x2 - x1
        xv = x3 - x1

        # get normal
        normal = cross(xu, xv)
        normal /= norm(normal)

        # get strength
        strength = FastMultipole.SVector{2,TF}(buffer[5, i_body], buffer[6, i_body])
        # strength = strength .* scale_strength

        # update values
        element = FastMultipole.Panel{FastMultipole.SourceDipole}
        FastMultipole.body_to_multipole_panel!(element, multipole_coefficients, harmonics, x0, xu, xv, normal, strength, expansion_order)

        # add wake panel
        if !system.semiinfinite_wake
            # which vertices connect to the wake
            idx1 = Int(buffer[end-7, i_body])
            if idx1 > 0
                # first vertex and wake
                Dax, Day, Daz = buffer[end-6, i_body], buffer[end-5, i_body], buffer[end-4, i_body]
                Da = FastMultipole.SVector{3}(Dax, Day, Daz)
                vs = (x1, x2, x3)
                v1 = vs[idx1]
                v1w = v1 + Da

                # second vertex and wake
                idx2 = Int(buffer[end-3, i_body])
                v2 = vs[idx2]
                Dbx, Dby, Dbz = buffer[end-2, i_body], buffer[end-1, i_body], buffer[end, i_body]
                Db = FastMultipole.SVector{3}(Dbx, Dby, Dbz)
                v2w = v2 + Db
                
                # body-to-multipole: first triangle
                x0 = v1 - center
                xu = v2 - v1
                xv = v1w - v1
                
                # get normal
                normal = cross(xu, xv)
                normal /= norm(normal)

                # doublet strength
                strength_doublet = FastMultipole.SVector{1,TF}(strength[2])

                # update values
                element = FastMultipole.Panel{FastMultipole.Dipole}
                FastMultipole.body_to_multipole_panel!(element, multipole_coefficients, harmonics, x0, xu, xv, normal, strength_doublet, expansion_order)

                # body-to-multipole: second triangle
                x0 = v1w - center
                xu = v2 - v1w
                xv = v2w - v1w

                # get normal
                normal = cross(xu, xv)
                normal /= norm(normal)

                # update values
                FastMultipole.body_to_multipole_panel!(element, multipole_coefficients, harmonics, x0, xu, xv, normal, strength_doublet, expansion_order)
            end
        end
    end
end

# multipole coefficients for vortex ring
function body_to_multipole_vl!(multipole_coefficients, harmonics, x1, x2, center, gamma, expansion_order)
    # delta
    xu = x2 - x1
    x0 = x1 - center

    # vector strength
    Gamma = gamma * xu / norm(xu)

    # update values
    FastMultipole.body_to_multipole_filament!(FastMultipole.Filament{FastMultipole.Vortex}, multipole_coefficients, harmonics, x0, xu, Gamma, expansion_order)
end

function FastMultipole.body_to_multipole!(system::AbstractBody{VortexRing, NK}, multipole_coefficients, source_buffer::Matrix{TF}, center, bodies_index, harmonics, expansion_order) where {TF, NK}
    # loop over bodies
    for i_source in bodies_index
       
        # get vertices
        v1x = source_buffer[5+NK, i_source]
        v1y = source_buffer[6+NK, i_source]
        v1z = source_buffer[7+NK, i_source]
        v1 = FastMultipole.StaticArrays.SVector{3,TF}(v1x, v1y, v1z)
        v2x = source_buffer[8+NK, i_source]
        v2y = source_buffer[9+NK, i_source]
        v2z = source_buffer[10+NK, i_source]
        v2 = FastMultipole.StaticArrays.SVector{3,TF}(v2x, v2y, v2z)
        v3x = source_buffer[11+NK, i_source]
        v3y = source_buffer[12+NK, i_source]
        v3z = source_buffer[13+NK, i_source]
        v3 = FastMultipole.StaticArrays.SVector{3,TF}(v3x, v3y, v3z)

        # extract strength from buffer
        gamma = source_buffer[5, i_source]

        # side 1
        body_to_multipole_vl!(multipole_coefficients, harmonics, v1, v2, center, gamma, expansion_order)

        # side 2
        body_to_multipole_vl!(multipole_coefficients, harmonics, v2, v3, center, gamma, expansion_order)

        # side 3
        body_to_multipole_vl!(multipole_coefficients, harmonics, v3, v1, center, gamma, expansion_order)

    end
end

# multipole coefficients for vortex ring + uniform vortex sheet
function FastMultipole.body_to_multipole!(::AbstractBody{Union{VortexRing, UniformVortexSheet}, 2, <:Any}, multipole_coefficients, source_buffer::Matrix{TF}, center, bodies_index, harmonics, expansion_order) where {TF}
    # loop over bodies
    for i_source in bodies_index
       
        # get vertices
        v1x = source_buffer[5+NK, i_source]
        v1y = source_buffer[6+NK, i_source]
        v1z = source_buffer[7+NK, i_source]
        v2x = source_buffer[8+NK, i_source]
        v2y = source_buffer[9+NK, i_source]
        v2z = source_buffer[10+NK, i_source]
        v3x = source_buffer[11+NK, i_source]
        v3y = source_buffer[12+NK, i_source]
        v3z = source_buffer[13+NK, i_source]

        # get normal (cross(v2-v1, v3-v1))
        nx = (v2y - v1y)*(v3z - v1z) - (v2z - v1z)*(v3y - v1y)
        ny = (v2z - v1z)*(v3x - v1x) - (v2x - v1x)*(v3z - v1z)
        nz = (v2x - v1x)*(v3y - v1y) - (v2y - v1y)*(v3x - v1x)
        normal = SVector{3,TF}(nx, ny, nz)

        # assumble vectors
        v1 = SVector{3,TF}(v1x, v1y, v1z)
        v2 = SVector{3,TF}(v2x, v2y, v2z)
        v3 = SVector{3,TF}(v3x, v3y, v3z)

        # extract strength from buffer
        gamma_ring = source_buffer[5, i_source]
        gamma_sheet = SVector{3,TF}(source_buffer[6, i_source], source_buffer[7, i_source], source_buffer[8, i_source])

        # vortex ring contribution
        x0 = v1 - center
        xu = v2 - v1
        xv = v3 - v1
        FastMultipole.body_to_multipole_panel!(FastMultipole.Panel{FastMultipole.Dipole}, multipole_coefficients, harmonics, x0, xu, xv, normal, gamma_ring, expansion_order)

        # vortex sheet contribution
        FastMultipole.body_to_multipole_panel!(FastMultipole.Panel{FastMultipole.Vortex}, multipole_coefficients, harmonics, x0, xu, xv, normal, gamma_sheet, expansion_order)
    end
end

#--- additional FastMultipole interface functions for RigidWakeBody ---#

"""
    additional_source_system_to_buffer!(buffer, i_buffer, system::RigidWakeBody, i_body, ilast)

Updates 8 rows to the buffer in column `i_buffer` corresponding to the `i_body`th body,
beginning with row `ilast+1`:

1. Which of the 3 vertices of this source sheds the first wake vertex (-1 if no wake is shed).
2. (If a wake is shed,) the \$x\$ delta from the first wake-shedding vertex to the next wake point
3. (If a wake is shed,) the \$y\$ delta from the first wake-shedding vertex to the next wake point
4. (If a wake is shed,) the \$z\$ delta from the first wake-shedding vertex to the next wake point
5. Which of the 3 vertices of this source sheds the second wake vertex (-1 if no wake is shed).
6. (If a wake is shed,) the \$x\$ delta from the second wake-shedding vertex to the next wake point
7. (If a wake is shed,) the \$y\$ delta from the second wake-shedding vertex to the next wake point
8. (If a wake is shed,) the \$z\$ delta from the second wake-shedding vertex to the next wake point

**Returns**

- `ilast+nrows` where `nrows` is the number of updated rows

"""
function additional_source_system_to_buffer!(buffer, i_buffer, system::RigidWakeBody, i_body, ilast)
    # index of first node of shedding edge (TE1)
    idx1 = system.shedding_full[1, i_body]
    buffer[ilast+1, i_buffer] = idx1

    # wake delta from the first node (TE1), using correct Das column
    update_radius = zero(eltype(buffer))
    if idx1 > 0
        i_surf = system.shedding_full[3, i_body]
        das_col_1 = system.shedding_full[5, i_body]
        Dax, Day, Daz = view(system.Das[i_surf], :, das_col_1)
        buffer[ilast+2:ilast+4, i_buffer] .= (Dax, Day, Daz)
        update_radius = sqrt(Dax*Dax + Day*Day + Daz*Daz)
    end

    # index of second node of shedding edge (TE2)
    idx2 = system.shedding_full[2, i_body]
    buffer[ilast+5, i_buffer] = system.shedding_full[2, i_body]

    # wake delta from the second node (TE2), using correct Das column
    if idx2 > 0
        i_surf = system.shedding_full[3, i_body]
        das_col_2 = system.shedding_full[6, i_body]
        Dax, Day, Daz = view(system.Das[i_surf], :, das_col_2)
        buffer[ilast+6:ilast+8, i_buffer] .= (Dax, Day, Daz)
        update_radius = max(update_radius, sqrt(Dax*Dax + Day*Day + Daz*Daz))
    end

    if !system.semiinfinite_wake
        buffer[4, i_buffer] += update_radius # conservative update to accomadate wake size
    end

    return ilast+8
end

additional_data_per_body(system::RigidWakeBody) = 8 # for the two node indices of the shedding edge

#--- interface functions for O(N) solver ---#

function FastMultipole.value_to_strength!(source_buffer, ::RigidWakeBody{<:Any,1,<:Any}, i_body, value)
    source_buffer[5, i_body] = value
end

function FastMultipole.value_to_strength!(source_buffer, ::RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing},2,<:Any}, i_body, value)
    source_buffer[5, i_body] = zero(value)
    source_buffer[6, i_body] = value
end

function FastMultipole.strength_to_value(strength, source_system::RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing}, 2, <:Any})
    return strength[2] # strength of constant doublet
end

function FastMultipole.buffer_to_system_strength!(system::RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing},2,<:Any}, i_body, source_buffer, i_buffer)
    system.strength[i_body, 1] = zero(eltype(source_buffer))
    system.strength[i_body, 2] = source_buffer[6, i_buffer]
end

function FastMultipole.influence!(influence, target_buffer, source_system::RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing},2,<:Any}, source_buffer)
    influence .= view(target_buffer, 4, :)
end


function FastMultipole.buffer_to_target_system!(target_system::RigidWakeBody, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    if PS
        phi = target_buffer[4, i_buffer]
        target_system.potential[i_target] += phi
    end

    if VS
        vx = target_buffer[5, i_buffer]
        vy = target_buffer[6, i_buffer]
        vz = target_buffer[7, i_buffer]
        target_system.velocity[1, i_target] += vx
        target_system.velocity[2, i_target] += vy
        target_system.velocity[3, i_target] += vz
    end
end

function FastMultipole.extra_farfield!(target_buffer, target_bodies_index, source_system::RigidWakeBody{<:Any,NK,<:Any}, source_buffer, source_bodies_index, switch::FastMultipole.DerivativesSwitch{PS,GS,HS}) where {NK,PS,GS,HS}

    # loop over targets
    for i_target in target_bodies_index
        
        # store influence to stack memory
        TF = eltype(target_buffer)
        ϕ = zero(TF)
        u = FastMultipole.StaticArrays.@SVector zeros(TF, 3)
        ∇u = FastMultipole.StaticArrays.@SMatrix zeros(TF, 3, 3)
        
        # target location                
        t = FastMultipole.StaticArrays.SVector{3,Float64}(target_buffer[1, i_target], target_buffer[2, i_target], target_buffer[3, i_target])

        # loop over sources
        for i_source in source_bodies_index
            idx1 = Int(source_buffer[end-7, i_source]) # index of first node of shedding edge
            if idx1 > 0

                # Da
                Dax, Day, Daz = source_buffer[end-6, i_source],
                                source_buffer[end-5, i_source],
                                source_buffer[end-4, i_source]

                # first vertex
                v1x = source_buffer[4+NK+1+(idx1-1)*3, i_source]
                v1y = source_buffer[4+NK+2+(idx1-1)*3, i_source]
                v1z = source_buffer[4+NK+3+(idx1-1)*3, i_source]

                # second vertex
                idx2 = Int(source_buffer[end-3, i_source]) # index of second node of shedding edge

                # Db
                Dbx, Dby, Dbz = source_buffer[end-2, i_source],
                                source_buffer[end-1, i_source],
                                source_buffer[end, i_source]

                # vertex
                v2x = source_buffer[4+NK+1+(idx2-1)*3, i_source]
                v2y = source_buffer[4+NK+2+(idx2-1)*3, i_source]
                v2z = source_buffer[4+NK+3+(idx2-1)*3, i_source]

                # strength
                strength = get_strength_doublet(source_system, source_buffer, i_source)

                # compute farfield influence
                this_ϕ, this_u, this_∇u = induced_semiinfinite(t, ConstantDoublet, v1x, v1y, v1z, v2x, v2y, v2z, Dax, Day, Daz, strength, switch; kerneloffset=source_system.kerneloffset)
                if PS
                    ϕ += this_ϕ
                end
                if GS
                    u += this_u
                end
                if HS
                    ∇u += this_∇u
                end
            end
        end

        # update target buffer
        if PS
            target_buffer[4, i_target] += ϕ
        end
        if GS
            target_buffer[5, i_target] += u[1]
            target_buffer[6, i_target] += u[2]
            target_buffer[7, i_target] += u[3]
        end
        if HS
            target_buffer[8, i_target] += ∇u[1, 1]
            target_buffer[9, i_target] += ∇u[1, 2]
            target_buffer[10, i_target] += ∇u[1, 3]
            target_buffer[11, i_target] += ∇u[2, 1]
            target_buffer[12, i_target] += ∇u[2, 2]
            target_buffer[13, i_target] += ∇u[2, 3]
            target_buffer[14, i_target] += ∇u[3, 1]
            target_buffer[15, i_target] += ∇u[3, 2]
            target_buffer[16, i_target] += ∇u[3, 3]
        end
    end
end

##### END OF FASTMULTIPOLE BACKEND SUPPORT #####################################


################################################################################
# COMMON FUNCTIONS
################################################################################
"""
    `calc_shedding(grid::GridTriangleSurface{gt.Meshes.SimpleMesh},
trailingedge::Matrix; tolerance=1e2*eps())`

Given an unstructured `grid` and a collection of points (line) `trailingedge`,
it finds the points in `grid` that are closer than `tolerance` to the line,
and automatically builds a `shedding` matrix that can be used to shed the wake
from this trailing edge.

Note: It is important that the points in `trailingedge` have been previously
    sorted to be contiguous to each other, otherwise the resulting `shedding`
    might have panels that are not contiguous to each other, fail to recognize
    panels that are at the trailing edge, or unphysically large trailing
    vortices.
"""
function calc_shedding(nodes, cells, trailingedge::Union{Matrix, Function};
                            periodic::Bool=false,
                            tolerance=1e2*eps(), debug=false
                            )
    # Identify the nodes that are on the TE line
    TEindices = gt.identifyedge(nodes, trailingedge; tolerance=tolerance)
    TEindices = [nodei for (nodei, pointi) in TEindices]

    return calc_shedding(nodes, cells, TEindices, trailingedge; periodic, tolerance, debug)
end

function calc_shedding(grid::gt.GridTriangleSurface{G}, trailingedge::Union{Matrix, Function};
                            periodic::Bool=false,
                            tolerance=1e2*eps(), debug=false
                            ) where {G<:gt.Meshes.SimpleMesh}
    # Identify the nodes that are on the TE line
    TEindices = gt.identifyedge(grid._nodes, trailingedge; tolerance=tolerance)
    TEindices = [nodei for (nodei, pointi) in TEindices]

    return calc_shedding(grid._nodes, grid2cells(grid), TEindices, trailingedge; periodic, tolerance, debug)
end

function calc_shedding(nodes::Matrix, cells::Matrix{Int}, TEindices, trailingedge::Union{Matrix, Function};
                            periodic::Bool=false,
                            tolerance=1e2*eps(), debug=false
                            )

    topology = gt.getHalfEdgeTopology(nodes, cells)
    connec = cells

    # Return if no TE nodes were identified
    if length(TEindices)==0
        return noshedding
    end

    # Append first node at end if expected to be periodic
    if periodic
        push!(TEindices, TEindices[1])
    end

    # All node pairs that could form an edge at the TE
    paircandidates = zip(view(TEindices, 1:length(TEindices)-1), view(TEindices, 2:length(TEindices)))

    # All node pairs that actually form an edge at the TE
    pairs = [pair for pair in paircandidates if haskey(topology.edge4pair, pair)]

    # Fetch all the first halfedge of each edges (node pairs) along the TE
    halfedges = [gt.Meshes.half4pair(topology, pair) for pair in pairs]

    # Build shedding matrix
    shedding = zeros(Int, 6, length(halfedges))

    for (ei, halfedge) in enumerate(halfedges)

        pair = pairs[ei]

        # pi is the panel associated with this half edge
        # pj is the panel associated with the other half

        # Case: Single-sided edge
        if isnothing(halfedge.elem) || isnothing(halfedge.half.elem)

            if isnothing(halfedge.half.elem)
                pi = halfedge.elem
            else
                pi = halfedge.half.elem
            end

            # Declare the other half as inexistent
            pj = -1

        # Case: Two-sided edge
        else

            # Identify which panel is "on top" and which "bottom" by matching
            # the order of the node pair
            inds1 = view(connec, :, halfedge.elem)
            inds2 = view(connec, :, halfedge.half.elem)

            if (
                    (inds1[1]==pair[1] && inds1[2]==pair[2])
                    ||
                    (inds1[2]==pair[1] && inds1[3]==pair[2])
                    ||
                    (inds1[3]==pair[1] && inds1[1]==pair[2])
                )

                pi = halfedge.elem
                pj = halfedge.half.elem

            elseif (
                    (inds2[1]==pair[1] && inds2[2]==pair[2])
                    ||
                    (inds2[2]==pair[1] && inds2[3]==pair[2])
                    ||
                    (inds2[3]==pair[1] && inds2[1]==pair[2])
                )

                pi = halfedge.half.elem
                pj = halfedge.elem

            else
                error("Logic error: Could not match panel to node pair")
            end

        end

        # Nodes of first half
        nia = findfirst(globindex -> globindex==pair[2], view(connec, :, pi))  # Local-index of the first node
        nib = findfirst(globindex -> globindex==pair[1], view(connec, :, pi))  # Local-index of the second node

        # Nodes of other half
        if pj != -1
            nja = findfirst(globindex -> globindex==pair[1], view(connec, :, pj))  # Local-index of the second node
            njb = findfirst(globindex -> globindex==pair[2], view(connec, :, pj))  # Local-index of the first node
        else
            nja = njb = -1
        end

        shedding[:, ei] .= (pi, nia, nib, pj, nja, njb)

    end

    if debug
        display("TEindices")
        display(TEindices)
        display("paircandidates")
        display(collect(paircandidates))
        display("pairs")
        display(pairs)
        display("halfedges")
        display(halfedges)
        display("shedding")
        display(shedding)

        points = [nodes[:, i] for i in TEindices]
        gt.generateVTK("TEindices", points; keep_points=true)

        lines = [[i-1, j-1] for (i, j) in paircandidates]
        gt.generateVTK("paircandidates", eachcol(nodes); lines=lines)

        lines = [[i-1, j-1] for (i, j) in pairs]
        gt.generateVTK("pairs", eachcol(nodes); lines=lines)
    end

    return shedding
end





##### INTERNAL FUNCTIONS  ######################################################
function _get_wakestrength_mu(self::RigidWakeBody, i, isurf=1; stri=1)
    # if self.use_wake_strength
    #     return self.wake_strength[i], zero(self.wake_strength[i])
    # else
        strength1 = self.strength[self.shedding[isurf][1, i], stri]
        strength2 = self.shedding[isurf][4, i] != -1 ? self.strength[self.shedding[isurf][4, i], stri] : 0.0
        return strength1, strength2
    # end
end
function _get_wakestrength_mu(self::RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing},2,<:Any}, i, isurf=1)
    # if self.use_wake_strength
    #     return self.wake_strength[i], zero(self.wake_strength[i])
    # else
        stri = 2
        strength1 = self.strength[self.shedding[isurf][1, i], stri]
        strength2 = self.shedding[isurf][4, i] != -1 ? self.strength[self.shedding[isurf][4, i], stri] : 0.0
        return strength1, strength2
    # end
end
function _get_wakestrength_Gamma(self::RigidWakeBody, i, isurf=1; stri=1)
    # if self.use_wake_strength
    #     return self.wake_strength[i]
    # else
        strength1 = self.strength[self.shedding[isurf][1, i], stri]
        strength2 = self.shedding[isurf][4, i] != -1 ? self.strength[self.shedding[isurf][4, i], stri] : 0.0
        return strength1 - strength2
    # end
end
function _get_wakestrength_Gamma(self::RigidWakeBody{<:Union{ConstantSource, ConstantDoublet, VortexRing},2,<:Any}, i, isurf=1)
    # if self.use_wake_strength
    #     return self.wake_strength[i]
    # else
        stri = 2
        strength1 = self.strength[self.shedding[isurf][1, i], stri]
        strength2 = self.shedding[isurf][4, i] != -1 ? self.strength[self.shedding[isurf][4, i], stri] : 0.0
        return strength1 - strength2
    # end
end

# RigidWakeBody hook: contributes panel strength
# to the VTK output produced by the generic write_vtk(name, body::AbstractBody, ...).

function _write_vtk_other_fields!(vtm, name, body::RigidWakeBody, idx)

    # set wake length for visualization (if semi-infinite wake is enabled)
    if body.semiinfinite_wake
        for Das in body.Das
            for i in axes(Das, 2)
                Das[:, i] .*= SEMIINFINITE_LENGTH[]
            end
        end
    end

    # save wake as doublet panels
    strength_label = wake_name(body)
    for i_surf in eachindex(body.shedding)
        shedding = body.shedding[i_surf]
        Das = body.Das[i_surf]

        n_wakes = size(shedding, 2)
        points = zeros(typeof(body.nodes[1]), 3, 4 * n_wakes)
        cells = Vector{WriteVTK.MeshCell{WriteVTK.VTKCellTypes.VTKCellType, FastMultipole.SVector{4,Int64}}}(undef, n_wakes)
        strengths = zeros(typeof(body.strength[1]), n_wakes)

        for i in 1:n_wakes
            pi = shedding[1, i]
            nia, nib = shedding[2, i], shedding[3, i]

            idx1 = body.cells[nia, pi]
            idx2 = body.cells[nib, pi]

            p1 = body.nodes[:, idx1]
            p2 = body.nodes[:, idx2]

            p3 = p2 + Das[:, i]
            p4 = p1 + Das[:, i+1]

            points[:, 1 + 4*(i-1)] = p1
            points[:, 2 + 4*(i-1)] = p2
            points[:, 3 + 4*(i-1)] = p3
            points[:, 4 + 4*(i-1)] = p4

            cells[i] = WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_QUAD, FastMultipole.SVector{4,Int64}(1+4*(i-1), 2+4*(i-1), 3+4*(i-1), 4+4*(i-1)))

            mu_upper, mu_lower = _get_wakestrength_mu(body, i, i_surf)
            strengths[i] = mu_upper - mu_lower
        end

        WriteVTK.vtk_grid(vtm, name * "_tw.$i_surf.$idx.vtu", points, cells) do vtk
            vtk[strength_label] = strengths
        end
    end

    # restore original values of Da if they were modified for visualization
    if body.semiinfinite_wake
        for Das in body.Das
            for i in axes(Das, 2)
                Das[:, i] ./= SEMIINFINITE_LENGTH[]
            end
        end
    end
end

#### END OF LIFTING BODY  ######################################################
