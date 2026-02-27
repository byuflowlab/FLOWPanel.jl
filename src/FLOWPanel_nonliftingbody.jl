#=##############################################################################
# DESCRIPTION
    Definition of non-lifting paneled body types (implementations of
    AbstractBody).

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Jun 2018
  * License     : MIT License
=###############################################################################


################################################################################
# NON-LIFTING BODY TYPE
################################################################################

"""
  `NonLiftingBody{E::AbstractElement, N}(grid::gt.GridTriangleSurface)`

Non-lifting body that is solved using a combination of N panel elements.
`grid` is the grid surface (paneled geometry).

  **Properties**
  * `nnodes::Int`                       : Number of nodes
  * `ncells::Int`                       : Number of cells
  * `fields::Vector{String}`            : Available fields (solutions)
  * `Oaxis::Matrix`                     : Coordinate system of body w.r.t. global
  * `O::Vector`                         : Origin of body w.r.t. global

"""
struct NonLiftingBody{E, N, TF} <: AbstractBody{E, N, TF}

    # User inputs
    grid::gt.GridTriangleSurface              # Paneled geometry
    
    # Properties
    nnodes::Int                               # Number of nodes
    ncells::Int                               # Number of cells
    cells::Matrix{Int}                        # Cell connectivity (each column is a cell)
    fields::Array{String, 1}                  # Available fields (solutions)
    Oaxis::Array{TF,2}                  # Coordinate system of original grid
    O::Array{TF,1}                      # Position of CS of original grid

    # Internal variables
    strength::Array{TF, 2}              # strength[i,j] is the stength of the i-th panel with the j-th element type
    velocity::Array{TF,2}               # Velocity at control points
    CPoffset::Float64                   # Control point offset in normal direction
    kerneloffset::Float64               # Kernel offset to avoid singularities
    kernelcutoff::Float64               # Kernel cutoff to avoid singularities
    characteristiclength::Function      # Characteristic length of each panel
    watertight::Bool                     # Whether the body is watertight or not
    inside_offset::Float64               # Offset to compute inside control points
end

function NonLiftingBody{E, N, TF}(
                grid;
                nnodes=grid.nnodes, ncells=grid.ncells,
                cells=grid2cells(grid),
                fields=Array{String,1}(),
                Oaxis=Array{TF,2}(1.0I, 3, 3), O=zeros(TF,3),
                strength=zeros(grid.ncells, N),
                velocity=zeros(3, grid.ncells),
                CPoffset=1e-14,
                kerneloffset=1e-8,
                kernelcutoff=1e-14,
                characteristiclength=characteristiclength_unitary,
                check_mesh=true, watertight=false,
                inside_offset=1e-6
              ) where {E, N, TF}
    # check if mesh is watertight
    if check_mesh && typeof(grid.orggrid) <: gt.Meshes.Mesh
        mesh = grid.orggrid
        watertight = gt.isclosed(mesh)
    end

    return NonLiftingBody{E, N, TF}(
                grid,
                nnodes, ncells, cells,
                fields,
                Oaxis, O,
                strength,
                velocity,
                CPoffset,
                kerneloffset,
                kernelcutoff,
                characteristiclength,
                watertight,
                inside_offset
              )
end

function (NonLiftingBody{E})(grid::gt.GridTriangleSurface; optargs...) where {E}
    return NonLiftingBody{E, _count(E), eltype(grid._nodes)}(grid; optargs...)
end

function save(body::NonLiftingBody, args...; optargs...)
    return save_base(body, args...; optargs...)
end

calc_elprescribe(::NonLiftingBody{ConstantSource, 1}) = Tuple{Int,Float64}[]
calc_elprescribe(body::NonLiftingBody{VortexRing, 1}) = body.watertight ? [(1, 0.0)] : Tuple{Int,Float64}[]
calc_elprescribe(body::NonLiftingBody{ConstantDoublet, 1}) = body.watertight ? [(1, 0.0)] : Tuple{Int,Float64}[]

solved_field_name(::NonLiftingBody{ConstantSource, 1}) = "sigma"
solved_field_name(::NonLiftingBody{ConstantDoublet, 1}) = "mu"
solved_field_name(::NonLiftingBody{VortexRing, 1}) = "gamma"
solved_field_name(::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2}) = "mu"

#### END OF NON-LIFTING BODY  ##################################################








################################################################################
# CONSTANT-SOURCE SOLVER
################################################################################
function solve(self::NonLiftingBody{ConstantSource, 1},
                                              Uinfs::AbstractArray{<:Number, 2})

    if size(Uinfs, 2) != self.ncells
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    end

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    G = zeros(self.ncells, self.ncells)
    _G_U!(self, G, CPs, normals)

    # Solve system of equations
    sigma = _solve(self, normals, G, Uinfs)

    # Save solution
    self.strength[:, 1] .= sigma

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "sigma", "scalar", view(self.strength, :, 1), "cell")

    rhs = .- vec(sum(normals .* Uinfs, dims=1))
    return G, rhs
end

function _solve(::NonLiftingBody{ConstantSource, 1}, normals, G, Uinfs)

    # Define right-hand side
    lambda = [-dot(Uinf, normal) for (Uinf, normal) in
                                        zip(eachcol(Uinfs), eachcol(normals))]

    # Solve the system of equations
    sigma = G\lambda

    return sigma
end

"""
Computes the geometric matrix (left-hand side matrix of the system of equation)
and stores it under `G`.

**ARGUMENTS**
  * `G::Array{T,2}`                     : Pre-allocated output memory.
  * `CPs::Array{T,2}`                   : Control points.
  * `normals::Array{T,2}`               : Normal associated to every CP.
"""
function _G_U!(self::NonLiftingBody{ConstantSource, 1},
                    G::Arr1, CPs::Arr2, normals::Arr3;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2},
                        T3, Arr3<:AbstractArray{T3, 2}}

    println("=====\nHERE,,,,,!!!\n=====")
    _G_U_constantsource!(self, G, CPs, normals; optargs...)
end

function _G_U_constantsource!(self, G, CPs, normals; optargs...)
    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    elseif size(normals, 2)!=M
        error("normals matrix with invalid dimensions;"*
              " got $(size(normals)), expected (3, $M).")
    end

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Build geometric matrix
    for (pj, Gslice) in enumerate(eachcol(G))

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        U_constant_source(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          # view(G, :, pj);                  # Velocity of j-th panel on every CP
                          Gslice;
                          dot_with=normals,                  # Normal of every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          optargs...
                         )
    end
end

function _Uind!(self::NonLiftingBody{ConstantSource, 1}, targets, out, backend::DirectBackend;
                                                                     optargs...)


    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Iterates over panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Velocity of i-th panel on every target
        U_constant_source(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 1],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            optargs...
                         )
    end
end

function _phi!(self::NonLiftingBody{ConstantSource, 1}, targets, out, backend::DirectBackend; optargs...)

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Iterates over panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Potential of i-th panel on every target
        phi_constant_source(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 1],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            optargs...
                         )
    end
end

_get_Gdims(self::NonLiftingBody{ConstantSource, 1}) = (self.ncells, self.ncells)



################################################################################
# CONSTANT-DOUBLET SOLVER
################################################################################
function solve(self::NonLiftingBody{ConstantDoublet, 1},
                                              Uinfs::AbstractArray{<:Number, 2})

    if size(Uinfs, 2) != self.ncells
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    end

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    G = zeros(self.ncells, self.ncells)
    _G_U!(self, G, CPs, normals)

    # Solve system of equations
    mu, lambda = _solve(self, normals, G, Uinfs)

    # Save solution
    self.strength[:, 1] .= mu

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "mu", "scalar", view(self.strength, :, 1), "cell")

    return G, lambda
end

function _solve(::NonLiftingBody{ConstantDoublet, 1}, normals, G, Uinfs)

    # Define right-hand side
    lambda = [-dot(Uinf, normal) for (Uinf, normal) in
                                        zip(eachcol(Uinfs), eachcol(normals))]

    # Solve the system of equations
    mu = G\lambda

    return mu, lambda
end

function _G_U!(self::AbstractBody{<:Any,<:Any,TF}, kernel, G, CPs, normals, backend::FastMultipoleBackend; kerneloffset=1.0e-3, include_wake=true, optargs...) where TF
    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    end

    # Build geometric matrix
    strength = FastMultipole.StaticArrays.SVector{1,TF}(1.0) # unit strength
    derivatives_switch=FastMultipole.DerivativesSwitch(false,true,false) # only potential
    for i_source in 1:N

        # get vertices
        v1i, v2i, v3i = self.cells[1,i_source], self.cells[2,i_source], self.cells[3,i_source]
        v1x, v1y, v1z = self.grid._nodes[1, v1i], self.grid._nodes[2, v1i], self.grid._nodes[3, v1i]
        v2x, v2y, v2z = self.grid._nodes[1, v2i], self.grid._nodes[2, v2i], self.grid._nodes[3, v2i]
        v3x, v3y, v3z = self.grid._nodes[1, v3i], self.grid._nodes[2, v3i], self.grid._nodes[3, v3i]
        
        for i_target in 1:M
            # get target
            tx, ty, tz = CPs[1, i_target], CPs[2, i_target], CPs[3, i_target]
            target = FastMultipole.StaticArrays.SVector{3,TF}(tx, ty, tz)

            # compute influence
            _, u, _ = induced(target, kernel, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, strength, derivatives_switch; kerneloffset)

            # update G
            G[i_target, i_source] = u[1] * normals[1, i_target] + u[2] * normals[2, i_target] + u[3] * normals[3, i_target]
        end
    end

    # add wake influence
    if include_wake
        _G_U_wake!(self, kernel, G, CPs, normals, backend; kerneloffset, optargs...)
    end
end

function _G_U_wake!(self::AbstractBody{<:Any,<:Any,TF}, kernel, G, CPs, normals, backend::FastMultipoleBackend; optargs...) where TF
    return nothing
end

"""
Computes the geometric matrix (left-hand side matrix of the system of equation)
and stores it under `G`.

**ARGUMENTS**
  * `G::Array{T,2}`                     : Pre-allocated output memory.
  * `CPs::Array{T,2}`                   : Control points.
  * `normals::Array{T,2}`               : Normal associated to every CP.
"""
function _G_U!(self::NonLiftingBody{<:Union{ConstantDoublet, VortexRing}, 1},
                    G::Arr1, CPs::Arr2, normals::Arr3;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2},
                        T3, Arr3<:AbstractArray{T3, 2}}

    _G_U_constantdoublet!(self, G, CPs, normals; optargs...)
end

function _G_U_constantdoublet!(self, G, CPs, normals; optargs...)
    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    elseif size(normals, 2)!=M
        error("normals matrix with invalid dimensions;"*
              " got $(size(normals)), expected (3, $M).")
    end

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Build geometric matrix
    for (pj, Gslice) in enumerate(eachcol(G))

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        U_constant_doublet(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          # view(G, :, pj);                  # Velocity of j-th panel on every CP
                          Gslice;
                          dot_with=normals,                  # Normal of every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                          optargs...
                         )
    end
end

function _G_phi_constantdoublet!(self, G, CPs, backend::DirectBackend; optargs...)
    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    end

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Build geometric matrix
    for (pj, Gslice) in enumerate(eachcol(G))

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        phi_constant_doublet(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            1.0,                               # Unitary strength
                            CPs,                           # Targets
                            Gslice;                               # Outputs
                            optargs...
                         )
    end
end

function _G_phi!(self::AbstractBody{<:Any,<:Any,TF}, kernel, G, CPs, backend::FastMultipoleBackend; kerneloffset=1.0e-8, include_wake=true, optargs...) where TF
    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    end

    # Build geometric matrix
    strength = FastMultipole.StaticArrays.SVector{1,TF}(1.0) # unit strength
    derivatives_switch=FastMultipole.DerivativesSwitch(true,false,false) # only potential
    for i_source in 1:N

        # get vertices
        v1i, v2i, v3i = self.cells[1,i_source], self.cells[2,i_source], self.cells[3,i_source]
        v1x, v1y, v1z = self.grid._nodes[1, v1i], self.grid._nodes[2, v1i], self.grid._nodes[3, v1i]
        v2x, v2y, v2z = self.grid._nodes[1, v2i], self.grid._nodes[2, v2i], self.grid._nodes[3, v2i]
        v3x, v3y, v3z = self.grid._nodes[1, v3i], self.grid._nodes[2, v3i], self.grid._nodes[3, v3i]
        
        for i_target in 1:M
            # get target
            tx, ty, tz = CPs[1, i_target], CPs[2, i_target], CPs[3, i_target]
            target = FastMultipole.StaticArrays.SVector{3,TF}(tx, ty, tz)

            # compute influence
            phi, _ = induced(target, kernel, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, strength, derivatives_switch; kerneloffset)

            # update G
            G[i_target, i_source] = phi
        end
    end

    # add wake influence
    if include_wake
        _G_phi_wake!(self, kernel, G, CPs, backend; derivatives_switch, kerneloffset, optargs...)
    end
end

function _G_phi_wake!(self::AbstractBody{<:Any,<:Any,TF}, kernel, G, CPs, backend::FastMultipoleBackend; kerneloffset=1.0e-3, optargs...) where TF
    return nothing
end

function _G_phi_constantsource!(self, G, CPs; optargs...)
    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    end

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Build geometric matrix
    for (pj, Gslice) in enumerate(eachcol(G))

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        phi_constant_source(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            1.0,                               # Unitary strength
                            CPs,                               # Targets
                            Gslice;                            # Outputs
                            optargs...
                         )
    end
end

function _Uind!(self::NonLiftingBody{<:Union{ConstantDoublet, VortexRing}, 1}, targets, out, backend::DirectBackend;
                                                                     optargs...)


    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Iterates over panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Velocity of i-th panel on every target
        U_constant_doublet(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 1],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                            optargs...
                         )
    end
end

function _phi!(self::NonLiftingBody{ConstantDoublet, 1},
                                                       targets, out, backend::DirectBackend; optargs...)

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Iterates over panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Potential of i-th panel on every target
        phi_constant_doublet(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 1],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            optargs...
                         )
    end
end

_get_Gdims(self::NonLiftingBody{ConstantDoublet, 1}) = (self.ncells, self.ncells)



################################################################################
# CONSTANT-SOURCE+DOUBLET SOLVER
################################################################################

"""
    freestream_potential(U, x)

Return the potential of a uniform freestream `Uinf` at point `x`, assuming
φ(0) = 0. Both `U` and `x` must be 3-element vectors.
"""
function freestream_potential(x, Uinf)
    return dot(Uinf, x)
end

function solve(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2, <:Any},
                                          Uinfs::AbstractArray{<:Number, 2}; solve_type=2, backend=DirectBackend(), optargs...)

    if solve_type == 1
        # --------- Pre-Calculations
        N = self.ncells

        if size(Uinfs, 2) != N
            error("Invalid Uinfs; expected size (3, $(N)), got $(size(Uinfs))")
        end

        # Compute normals and control points
        normals = _calc_normals(self)
        CPs = _calc_controlpoints(self, normals, off=1e-6)
        # CPs_inside = _calc_controlpoints(self, normals; off=-self.inside_offset)

        # --------- Compute geometric matrix
        G = zeros(N, N)
        _G_U_constantsource!(self, G, CPs, normals; optargs...)
        # _G_phi_U!(self, G, CPs, CPs_inside, normals)

        # --------- RHS
        rhs = zeros(N)
        rhs .= sum(Uinfs .* normals; dims=1) |> vec .*= -1.0

        # --------- solver for source strengths
        sigma = G\rhs
        self.strength[:, 1] .= sigma
        self.strength[:, 2] .= 0.0

        # --------- use source strengths to evaluate ϕ (doublet RHS)
        # phis = G * sigma # perturbation potential (due to sources)
        rhs .= 0.0
        _phi!(self, CPs, rhs, backend; optargs..., offset=self.kerneloffset)
        
        # --------- add freestream potential
        # phis .+= vec(sum(Uinfs .* CPs, dims=1))
        
        # --------- construct influence matrix for doublets
        G .= 0.0
        _G_phi_constantdoublet!(self, G, CPs; optargs...)
        # _G_U_constantdoublet!(self, G, CPs, normals; optargs...)
        
        # --------- Solve system of equations
        mu = G\rhs
        
        # --------- Store results
        # Save solution
        self.strength[:, 1] .= 0.0
        self.strength[:, 2] .= mu

        _solvedflag(self, true)
        add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
        add_field(self, "sigma", "scalar", view(self.strength, :, 1), "cell")
        add_field(self, "mu", "scalar", view(self.strength, :, 2), "cell")

        return G, rhs

    elseif solve_type == 2 # dirichlet bc with doublet panels only

        # --------- Pre-Calculations
        N = self.ncells

        if size(Uinfs, 2) != N
            error("Invalid Uinfs; expected size (3, $(N)), got $(size(Uinfs))")
        end

        # Compute normals and control points
        normals = _calc_normals(self)
        # CPs = _calc_controlpoints(self, normals, off=1e-6)
        CPs_inside = _calc_controlpoints(self, normals; off=0.0)# -self.inside_offset)

        # --------- Compute geometric matrix
        G = zeros(N, N)
        _G_phi!(self, ConstantDoublet, G, CPs_inside, backend; optargs...)
        for i in 1:N
            G[i,i] = 0.5
        end

        # --------- RHS
        # rhs = zeros(N)
        # rhs .= sum(Uinfs .* normals; dims=1) |> vec .*= -1.0
        rhs = -vec(sum(Uinfs .* CPs_inside; dims=1)) # dirichlet bc: constant interior potential = 0

        # --------- solver for doublet strengths
        mu = G\rhs
        self.strength[:, 1] .= 0.0
        self.strength[:, 2] .= mu
        
        # --------- Store results
        _solvedflag(self, true)
        add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
        add_field(self, "sigma", "scalar", view(self.strength, :, 1), "cell")
        add_field(self, "mu", "scalar", view(self.strength, :, 2), "cell")

        return G, rhs

    elseif solve_type == 3
        #--- Morino formulation ---#

        # get ∂ϕ∂n (boundary condition)
        normals = _calc_normals(self)
        ∂ϕ∂n = -vec(sum(Uinfs .* normals; dims=1))

        # compute source term in BIE
        self.strength[:, 1] .= ∂ϕ∂n
        self.strength[:, 2] .= 0.0

        # # source-influence
        CPs = _calc_controlpoints(self, normals; off=1e-8)
        ϕ_source = zeros(self.ncells)
        _phi!(self, CPs, ϕ_source, backend; optargs..., offset=self.kerneloffset)

        # compute geometric matrix for doublets
        N = self.ncells
        G = zeros(N, N)
        _G_phi!(self, ConstantDoublet, G, CPs, backend; optargs...)
        # G .*= -1.0 # move to LHS

        # add I term
        for i in 1:N
            G[i,i] += 1.0
        end

        # solve for doublet strengths
        μ = G \ ϕ_source

        # store strengths
        # self.strength[:, 1] .*= 2 # σ = 2 ∂ϕ/∂n
        self.strength[:, 2] .= μ # μ = 2 * solution of BIE on the boundary
        # NOTE: the factor of 2 is because the n-body solver is meant for the fluid domain,
        # and has an extra 1/2 factor

        # --------- Store results
        _solvedflag(self, true)
        add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
        add_field(self, "sigma", "scalar", view(self.strength, :, 1), "cell")
        add_field(self, "mu", "scalar", view(self.strength, :, 2), "cell")
        add_field(self, "normals", "vector", collect(eachcol(normals)), "cell")

        return G, ϕ_source

    end
end

# SIMULTANEOUS SOLVE FOR SOURCES AND DOUBLETS
# function solve(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2, <:Any},
#                                           Uinfs::AbstractArray{<:Number, 2}; _G_fun=_G_phi_U!)

#     # --------- Pre-Calculations
#     N = self.ncells

#     if size(Uinfs, 2) != N
#         error("Invalid Uinfs; expected size (3, $(N)), got $(size(Uinfs))")
#     end

#     # Compute normals and control points
#     normals = _calc_normals(self)
#     CPs = _calc_controlpoints(self, normals)
#     CPs_inside = _calc_controlpoints(self, normals; off=-self.inside_offset)

#     # --------- Compute geometric matrix
#     G = zeros(2N, 2N)
#     _G_fun(self, G, CPs, CPs_inside, normals)
#     # _G_phi_U!(self, G, CPs, CPs_inside, normals)

#     # --------- Define right-hand side
#     rhs = zeros(2N)
#     phis_inf = sum(Uinfs .* CPs_inside, dims=1)
#     rhs[1:N] = -phis_inf[:]
#     U_normal = sum(Uinfs .* normals, dims=1)
#     rhs[N+1:2N] .-= U_normal[:]

#     # --------- Solve system of equations
#     sigma_mu = G\rhs

#     # --------- Store results
#     # Save solution
#     self.strength[:, 1] .= sigma_mu[1:N]
#     self.strength[:, 2] .= sigma_mu[N+1:2N]

#     _solvedflag(self, true)
#     add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
#     add_field(self, "sigma", "scalar", view(self.strength, :, 1), "cell")
#     add_field(self, "mu", "scalar", view(self.strength, :, 2), "cell")

#     return G, rhs
# end

function _solve(::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2}, G, nphis)

    # Solve the system of equations
    mu = G\nphis

    return mu
end

"""
Computes the geometric matrix (left-hand side matrix of the system of equation)
and stores it under `G`.

**ARGUMENTS**
  * `G::Array{T,2}`                     : Pre-allocated output memory.
  * `CPs::Array{T,2}`                   : Control points.
  * `normals::Array{T,2}`               : Normal associated to every CP.
"""
function _G_phi!(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2},
                    G::Arr1, CPs::Arr2;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2}}

    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    elseif size(normals, 2)!=M
        error("normals matrix with invalid dimensions;"*
              " got $(size(normals)), expected (3, $M).")
    end

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Build geometric matrix
    for (pj, Gslice) in enumerate(eachcol(G))

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        phi_constant_doublet(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          Gslice;                            # Potential of j-th panel on every CP
                          optargs...
                         )
    end
end


"""
Computes the geometric matrix (left-hand side matrix of the system of equation)
and stores it under `G`.

**ARGUMENTS**
  * `G::Array{T,2}`                     : Pre-allocated output memory.
  * `CPs::Array{T,2}`                   : Control points.
  * `normals::Array{T,2}`               : Normal associated to every CP.
"""
function _G_phi_U!(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2},
                    G::Arr1, CPs::Arr2, CPs_inside::Arr2, normals::Arr3;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2},
                        T3, Arr3<:AbstractArray{T3, 2}}

    N = self.ncells * 2
    M = size(CPs, 2) * 2

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    elseif size(normals, 2)!=M>>1
        error("normals matrix with invalid dimensions;"*
              " got $(size(normals)), expected (3, $M).")
    end

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Build geometric matrix: source panels as sending
    for pj in 1:self.ncells

        # get source panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        # get slice of G
        Gslice = view(G, 1:self.ncells, pj)

        phi_constant_source(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            1.0,                               # Unitary strength
                            CPs_inside,                           # Targets
                            Gslice;                               # Outputs
                            offset=self.kerneloffset           # Offset of kernel to avoid singularities
                         )

        # get slice of G
        Gslice = view(G, self.ncells+1:2*self.ncells, pj)

        U_constant_source(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          # view(G, :, pj);                  # Velocity of j-th panel on every CP
                          Gslice;
                          dot_with=normals,                  # Normal of every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          optargs...
                         )
    end

    # Build geometric matrix: doublet panels as sending
    for pj in 1:self.ncells

        # get source panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        # get slice of G
        Gslice = view(G, 1:self.ncells, pj + self.ncells)

        phi_constant_doublet(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            1.0,                               # Unitary strength
                            CPs_inside,                        # Targets
                            Gslice                                # Result
                         )

        # get slice of G
        Gslice = view(G, self.ncells+1:2*self.ncells, pj + self.ncells)

        U_constant_doublet(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          # view(G, :, pj);                  # Velocity of j-th panel on every CP
                          Gslice;
                          dot_with=normals,                  # Normal of every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          optargs...
                         )
    end
end
function _G_phi_phi!(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2},
                    G::Arr1, CPs::Arr2, CPs_inside::Arr2, normals::Arr3;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2},
                        T3, Arr3<:AbstractArray{T3, 2}}

    N = self.ncells * 2
    M = size(CPs, 2) * 2

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    elseif size(normals, 2)!=M>>1
        error("normals matrix with invalid dimensions;"*
              " got $(size(normals)), expected (3, $M).")
    end

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Build geometric matrix: source panels as sending
    for pj in 1:self.ncells

        # get source panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        # get slice of G
        Gslice = view(G, 1:self.ncells, pj)

        phi_constant_source(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            1.0,                               # Unitary strength
                            CPs_inside,                           # Targets
                            Gslice;                               # Outputs
                            offset=self.kerneloffset           # Offset of kernel to avoid singularities
                         )

        # get slice of G
        Gslice = view(G, self.ncells+1:2*self.ncells, pj)

        U_constant_source(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          # view(G, :, pj);                  # Velocity of j-th panel on every CP
                          Gslice;
                          dot_with=normals,                  # Normal of every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          optargs...
                         )
    end

    # Build geometric matrix: doublet panels as sending
    for pj in 1:self.ncells

        # get source panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        # get slice of G
        Gslice = view(G, 1:self.ncells, pj + self.ncells)

        phi_constant_doublet(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            1.0,                               # Unitary strength
                            CPs_inside,                        # Targets
                            Gslice                                # Result
                         )

        # get slice of G
        Gslice = view(G, self.ncells+1:2*self.ncells, pj + self.ncells)

        U_constant_doublet(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          # view(G, :, pj);                  # Velocity of j-th panel on every CP
                          Gslice;
                          dot_with=normals,                  # Normal of every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          optargs...
                         )
    end
end

function _Uind!(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2}, targets, out, backend::DirectBackend;
                                                                     optargs...)


    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Iterates over panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Velocity of i-th panel on every target
        U_constant_source(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 1],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            optargs...
                         )

        # Velocity of i-th panel on every target
        U_constant_doublet(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 2],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                            optargs...
                         )
    end
end

function _phi!(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2},
                                                       targets, out, backend::DirectBackend; optargs...)

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Iterates over panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Potential of i-th panel on every target
        phi_constant_source(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 1],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            optargs...
                         )

        # Potential of i-th panel on every target
        phi_constant_doublet(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 2],               # Unitary strength
                            targets,                           # Targets
                            out                               # Outputs
                         )
    end
end

_get_Gdims(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2}) = (self.ncells, self.ncells)




################################################################################
# FASTMULTIPOLE BACKEND SUPPORT
################################################################################
# function _Uind!(self::NonLiftingBody, targets, out, backend::FastMultipoleBackend; optargs...)
#     # wrap targets in a probe system
#     TF = eltype(targets)
#     potential = Vector{TF}(undef, 0) # unused
#     hessian = Array{TF, 3}(undef, 0, 0, 0)  # unused
#     probe_system = FastMultipole.ProbeSystemArray(targets, potential, out, hessian)

#     # perform N-body calculation
#     FastMultipole.fmm!(probe_system, self; expansion_order=backend.expansion_order,
#                                         multipole_acceptance=backend.multipole_acceptance,
#                                         leaf_size_source=backend.leaf_size,
#                                         hessian=false,
#                                         gradient=true, 
#                                         scalar_potential=false)

#     return nothing
# end

# function _phi!(self::NonLiftingBody, targets, out, backend::FastMultipoleBackend; optargs...)
#     # wrap targets in a probe system
#     TF = eltype(targets)
#     velocity = Array{TF, 2}(undef, 0, 0)  # unused
#     hessian = Array{TF, 3}(undef, 0, 0, 0)  # unused
#     probe_system = FastMultipole.ProbeSystemArray(targets, out, velocity, hessian)

#     # perform N-body calculation
#     FastMultipole.fmm!(probe_system, self; expansion_order=backend.expansion_order,
#                                         multipole_acceptance=backend.multipole_acceptance,
#                                         leaf_size_source=backend.leaf_size,
#                                         hessian=false,
#                                         gradient=false, 
#                                         scalar_potential=true)

#     return nothing
# end

FastMultipole.has_vector_potential(::AbstractBody{ConstantSource, 1}) = false

FastMultipole.has_vector_potential(::AbstractBody{ConstantDoublet, 1}) = false

FastMultipole.has_vector_potential(::AbstractBody{Union{ConstantSource, ConstantDoublet}, 2, <:Any}) = false

FastMultipole.body_to_multipole!(system::AbstractBody{ConstantSource, 1, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Source}, system, args...)

FastMultipole.body_to_multipole!(system::AbstractBody{ConstantDoublet, 1, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Dipole}, system, args...; scale_strength=1.0)

FastMultipole.body_to_multipole!(system::AbstractBody{Union{ConstantSource,ConstantDoublet}, 2, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.SourceDipole}, system, args...; scale_strength=FastMultipole.StaticArrays.SVector(1.0, 1.0))

##### END OF FASTMULTIPOLE BACKEND SUPPORT #####################################

################################################################################
# ABSTRACT SOLVER INTERFACE
################################################################################

function solve2!(self::NonLiftingBody{TK,1,TFG}, Uinfs::Array{TFS, 2}, solver::AbstractMatrixfulSolver{false};
        backend=DirectBackend(),
        update_G::Bool=false,   # Whether to update the influence matrix G
        strength_name=get_strength_name(self),  # Name of the strength field to solve for
        optargs...              # Additional optional arguments to _G_U!
    ) where {TK, TFG, TFS<:Real}

    # formatting assertions
    @assert size(self.strength, 2) == 1 "AbstractBody{<:Any, 1} expected to have single strength per element; got size $(size(self.strength, 2))."

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)
    CPs_inside = _calc_controlpoints(self, normals; off=-1e-10)
    
    # update influence matrix (if requested)
    if update_G

        # Update geometric matrix (left-hand-side influence matrix)
        _G_U!(self, TK, solver.G, CPs, normals, backend; optargs...)
    end

    # generate RHS
    TF = promote_type(TFG, TFS)
    RHS = zeros(TF, self.ncells)
    calc_bc_noflowthrough!(RHS, Uinfs, normals)

    # Solve system of equations
    solution = zeros(TF, self.ncells)
    solve_matrix!(solution, solver.G, RHS, solver)
    
    # set solved flag
    _solvedflag(self, true)
    
    # Assign solution to body element strengths
    # _assign_elementstrengths!(self, solution)
    self.strength .= solution

    # verify that source strength is equal to -∂ϕ/∂n on the boundary
    us_outside = zeros(TF, 3, self.ncells)
    _Uind!(self, CPs, us_outside, backend; optargs...)
    us_inside = zeros(TF, 3, self.ncells)
    _Uind!(self, CPs_inside, us_inside, backend; optargs...)
    add_field(self, "us_inside", "vector", collect(eachcol(us_inside)), "cell")
    add_field(self, "delta_u_normal", "scalar", vec(sum((us_outside - us_inside) .* normals, dims=1)), "cell")

    # save solution fields
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, strength_name, "scalar", view(self.strength, :, 1), "cell")
    add_field(self, "normals", "vector", collect(eachcol(normals)), "cell")

    return nothing
end

function FastMultipole.value_to_strength!(source_buffer, ::NonLiftingBody, i_body, value, rlx)
    prev_value = source_buffer[5, i_body]
    source_buffer[5, i_body] = rlx*value + (1.0 - rlx)*prev_value
end

function FastMultipole.value_to_strength!(source_buffer, ::NonLiftingBody, i_body, value)
    source_buffer[5, i_body] = value
end

function FastMultipole.target_influence_to_buffer!(target_buffer, i_buffer, derivatives_switch, target_system::NonLiftingBody, i_target)
    vx, vy, vz = target_system.velocity[1, i_target], target_system.velocity[2, i_target], target_system.velocity[3, i_target]
    target_buffer[5, i_buffer] = vx
    target_buffer[6, i_buffer] = vy
    target_buffer[7, i_buffer] = vz
end

function FastMultipole.buffer_to_target_system!(target_system::NonLiftingBody, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    vx, vy, vz = target_buffer[5, i_buffer], target_buffer[6, i_buffer], target_buffer[7, i_buffer]
    target_system.velocity[1, i_target] = vx
    target_system.velocity[2, i_target] = vy
    target_system.velocity[3, i_target] = vz
end

##### END OF ABSTRACT SOLVER INTERFACE ##########################################


################################################################################
# COMMON FUNCTIONS
################################################################################
"""
  `generate_loft(bodytype::Type{B}, args...; bodyoptargs=(), optargs...)  where
{B<:NonLiftingBody}`

Generates a lofted non-lifting body. See documentation of
`GeometricTools.generate_loft` for a description of the arguments of this
function.
"""
function generate_loft(bodytype::Type{B}, args...; bodyoptargs=(),
                        dimsplit::Int64=2, optargs...) where {B<:NonLiftingBody}
    # Loft the surface geometry
    grid = gt.generate_loft(args...; optargs...)

    # Split the quadrialateral panels into triangles
    # dimsplit = 2              # Dimension along which to split
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)

    return bodytype(triang_grid; bodyoptargs...)
end

"""
  `generate_revolution_(bodytype::Type{B}, args...; bodyoptargs=(), optargs...)
where {B<:NonLiftingBody}`

Generates a non-lifting body of a body of revolution. See documentation of
`GeometricTools.surface_revolution` for a description of the arguments of this
function.
"""
function generate_revolution(bodytype::Type{B}, args...; bodyoptargs=(),
                             dimsplit::Int64=2, loop_dim::Int64=2,
                             optargs...)  where {B<:NonLiftingBody}
    # Revolve the geometry
    grid = gt.surface_revolution(args...; loop_dim=loop_dim, optargs...)

    # Split the quadrialateral panels into triangles
    # dimsplit = 2              # Dimension along which to split
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)

    return bodytype(triang_grid; bodyoptargs...)
end

##### END OF COMMON FUNTIONS ###################################################
