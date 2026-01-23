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
    CPoffset::Float64                   # Control point offset in normal direction
    kerneloffset::Float64               # Kernel offset to avoid singularities
    kernelcutoff::Float64               # Kernel cutoff to avoid singularities
    characteristiclength::Function      # Characteristic length of each panel
    watertight::Bool                     # Whether the body is watertight or not

end

function NonLiftingBody{E, N, TF}(
                grid;
                nnodes=grid.nnodes, ncells=grid.ncells,
                cells=grid2cells(grid),
                fields=Array{String,1}(),
                Oaxis=Array{TF,2}(1.0I, 3, 3), O=zeros(TF,3),
                strength=zeros(grid.ncells, N),
                CPoffset=1e-14,
                kerneloffset=1e-8,
                kernelcutoff=1e-14,
                characteristiclength=characteristiclength_unitary,
                check_mesh=true, watertight=false
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
                CPoffset,
                kerneloffset,
                kernelcutoff,
                characteristiclength,
                watertight
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
    mu = _solve(self, normals, G, Uinfs)

    # Save solution
    self.strength[:, 1] .= mu

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "mu", "scalar", view(self.strength, :, 1), "cell")
end

function _solve(::NonLiftingBody{ConstantDoublet, 1}, normals, G, Uinfs)

    # Define right-hand side
    lambda = [-dot(Uinf, normal) for (Uinf, normal) in
                                        zip(eachcol(Uinfs), eachcol(normals))]

    # Solve the system of equations
    mu = G\lambda

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
function _G_U!(self::NonLiftingBody{ConstantDoublet, 1},
                    G::Arr1, CPs::Arr2, normals::Arr3;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2},
                        T3, Arr3<:AbstractArray{T3, 2}}

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

function _Uind!(self::NonLiftingBody{ConstantDoublet, 1}, targets, out, backend::DirectBackend;
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
                                                       targets, out; optargs...)

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
function solve(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2},
                                          Uinfs::AbstractArray{<:Number, 2})

    # --------- Pre-Calculations
    N = self.ncells

    if size(Uinfs, 2) != N
        error("Invalid Uinfs; expected size (3, $(N)), got $(size(Uinfs))")
    end

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # --------- Compute strength of sources
    # We use U_z(x_cp) ≈ −σ/2 to set up σ = 2 U∞⋅n
    # sigma = [2*dot(n, Uinf) for (n, Uinf) in zip(eachcol(normals), eachcol(Uinfs))]
    # sigma = [-2*dot(n, Uinf) for (n, Uinf) in zip(eachcol(normals), eachcol(Uinfs))]

    sigma = [-dot(n, Uinf) for (n, Uinf) in zip(eachcol(normals), eachcol(Uinfs))]

    #=
    # Compute geometric matrix of source terms (left-hand-side influence matrix)
    Gsrc = zeros(N, N)
    _G_U!(self, Gsrc, CPs, normals)

    # Define right-hand side
    lambda = [-dot(Uinf, normal) for (Uinf, normal) in
    # lambda = [dot(Uinf, normal) for (Uinf, normal) in           #
                                        zip(eachcol(Uinfs), eachcol(normals))]

    # Solve the system of equations
    sigma = Gsrc\lambda
    =#


    # --------- Compute strength of doublets
    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Compute potential field induced by sources at each control point (right-hand side)
    nphis = deepcopy(sigma)
    nphis .= 0

    for (pj, sgm) in enumerate(sigma)

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        phi_constant_source(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          -sgm,                              # Source strength (flipped to move the potential to the RHS)
                          # sgm,
                          CPs,                               # Targets
                          nphis;                             # Potential of j-th panel on every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                         )
    end

    #=
    # Compute geometric matrix (left-hand-side influence matrix)
    G = zeros(self.ncells, self.ncells)
    _G_phi!(self, G, CPs)

    # Solve system of equations
    mu = _solve(self, G, nphis)
    =#

    mu = nphis


    #=
    # --------- Compute coupled strengths
    G = zeros(2*N, 2*N)

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    # Build geometric block of U
    for (pj, Gslice) in enumerate(eachcol(view(G, 1:N, 1:N)))

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
                         )
    end
    for (pj, Gslice) in enumerate(eachcol(view(G, 1:N, N+1:2*N)))

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
                         )
    end

    # Build geometric block of phi
    for (pj, Gslice) in enumerate(eachcol(view(G, N+1:2*N, 1:N)))

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        phi_constant_source(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          Gslice;                            # Potential of j-th panel on every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                         )
    end
    for (pj, Gslice) in enumerate(eachcol(view(G, N+1:2*N, N+1:2*N)))

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        phi_constant_doublet(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          1.0,                               # Unitary strength
                          CPs,                               # Targets
                          Gslice;                            # Potential of j-th panel on every CP
                         )
    end

    # Define right-hand side
    lambda = zeros(2*N)
    lambda[1:N] .= (-dot(Uinf, normal) for (Uinf, normal) in
                                        zip(eachcol(Uinfs), eachcol(normals)))

    # Solve the system of equations
    sigmamu = G\lambda
    sigma = sigmamu[1:N]
    mu = sigmamu[N+1:end]
    =#

    # --------- Store results
    # Save solution
    self.strength[:, 1] .= sigma
    self.strength[:, 2] .= mu

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "sigma", "scalar", view(self.strength, :, 1), "cell")
    add_field(self, "mu", "scalar", view(self.strength, :, 2), "cell")
end

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
function _G_U!(self::NonLiftingBody{Union{ConstantSource, ConstantDoublet}, 2},
                    G::Arr1, CPs::Arr2, normals::Arr3;
                    optargs...
               ) where{ T1, Arr1<:AbstractArray{T1, 2},
                        T2, Arr2<:AbstractArray{T2, 2},
                        T3, Arr3<:AbstractArray{T3, 2}}

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
                                                       targets, out; optargs...)

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
                            out;                               # Outputs
                            optargs...
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

FastMultipole.body_to_multipole!(system::AbstractBody{ConstantSource, 1, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Source}, system, args...)

FastMultipole.body_to_multipole!(system::AbstractBody{ConstantDoublet, 1, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.Dipole}, system, args...)

FastMultipole.body_to_multipole!(system::AbstractBody{Union{ConstantSource,ConstantDoublet}, 2, <:Any}, args...) =
    FastMultipole.body_to_multipole!(FastMultipole.Panel{FastMultipole.SourceDipole}, system, args...)

##### END OF FASTMULTIPOLE BACKEND SUPPORT #####################################

################################################################################
# ABSTRACT SOLVER INTERFACE
################################################################################

function solve2!(self::NonLiftingBody{<:Any,1,TFG}, Uinfs::Array{TFS, 2}, solver::AbstractMatrixfulSolver{false};
        update_G::Bool=false,   # Whether to update the influence matrix G
        strength_name=get_strength_name(self),  # Name of the strength field to solve for
        optargs...              # Additional optional arguments to _G_U!
    ) where {TFG, TFS<:Real}

    # formatting assertions
    @assert size(self.strength, 2) == 1 "AbstractBody{<:Any, 1} expected to have single strength per element; got size $(size(self.strength, 2))."

    # update influence matrix (if requested)
    normals = _calc_normals(self)
    if update_G
        # Compute normals and control points
        CPs = _calc_controlpoints(self, normals)

        # Update geometric matrix (left-hand-side influence matrix)
        _G_U!(self, solver.G, CPs, normals; optargs...)
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
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, strength_name, "scalar", view(self.strength, :, 1), "cell")

    return nothing
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
