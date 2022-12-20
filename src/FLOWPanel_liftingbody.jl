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
  * `nnodes::Int64`                     : Number of nodes
  * `ncells::Int64`                     : Number of cells
  * `fields::Array{String, 1}`          : Available fields (solutions)
  * `Oaxis::Array{T<:Real, 2} `         : Coordinate system of original grid
  * `O::Array{T<:Real,1} `              : Position of CS of original grid
  * `ncellsTE::Int64`                   : Number of cells along trailing edge
  * `nnodesTE::Int64`                   : Number of nodes along trailing edge

"""
struct RigidWakeBody{E, N} <: AbstractLiftingBody{E, N}

    # User inputs
    grid::gt.GridTriangleSurface              # Paneled geometry
    shedding::Array{Int, 2}                   # Indicates edges along which to
                                              # shed the wake
    # Properties
    nnodes::Int64                             # Number of nodes
    ncells::Int64                             # Number of cells
    nsheddings::Int64                         # Number of shedding edges
    fields::Array{String, 1}                  # Available fields (solutions)
    Oaxis::Array{<:Number,2}                  # Coordinate system of original grid
    O::Array{<:Number,1}                      # Position of CS of original grid

    # Internal variables
    strength::Array{<:Number, 2}              # strength[i,j] is the stength of the i-th panel with the j-th element type
    CPoffset::Float64                         # Control point offset in normal direction
    kerneloffset::Float64                     # Kernel offset to avoid singularities
    kernelcutoff::Float64                     # Kernel cutoff to avoid singularities
    characteristiclength::Function            # Characteristic length of each panel

    RigidWakeBody{E, N}(
                    grid, shedding;
                    nnodes=grid.nnodes, ncells=grid.ncells,
                    nsheddings=size(shedding,2),
                    fields=Array{String,1}(),
                    Oaxis=Array(1.0I, 3, 3), O=zeros(3),
                    strength=zeros(grid.ncells, N),
                    CPoffset=1e-14,
                    kerneloffset=1e-8,
                    kernelcutoff=1e-14,
                    characteristiclength=characteristiclength_unitary
                  ) where {E, N} = _checkTE(grid, shedding) ? new(
                    grid, shedding,
                    nnodes, ncells,
                    nsheddings,
                    fields,
                    Oaxis, O,
                    strength,
                    CPoffset,
                    kerneloffset,
                    kernelcutoff,
                    characteristiclength
                  ) : error("Got invalid trailing edge.")

end

function (RigidWakeBody{E})(grid, shedding; optargs...) where {E}
    return RigidWakeBody{E, _count(E)}(grid, shedding; optargs...)
end

function (RigidWakeBody{E})(grid; optargs...) where {E}
    return RigidWakeBody{E}(grid, zeros(Int, 6, 0); optargs...)
end



function save(body::RigidWakeBody, args...;
                out_wake::Bool=true, debug::Bool=false,
                wake_len::Number=1.0,
                wake_panel::Bool=true, optargs...)


    str = ""
    str *= save_base(body, args...; debug=debug, optargs...)

    # Output the wake
    if out_wake || debug
        str *= _savewake(body, args...; len=wake_len, panel=wake_panel, optargs...)
    end

    return str
end




################################################################################
# VORTEX RING SOLVER
################################################################################
function solve(self::RigidWakeBody{VortexRing, 1},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2},
                Dbs::AbstractMatrix{T3};
                solver=solve_ludiv!, solver_optargs=()
                ) where {T1, T2, T3}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings)), got $(size(Das))")
    elseif size(Dbs) != (3, self.nsheddings)
        error("Invalid Dbs;"*
              " expected size (3, $(self.nsheddings)), got $(size(Dbs))")
    end

    T = promote_type(T1, T2, T3)

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    G = zeros(T, self.ncells, self.ncells)
    _G_U!(self, G, CPs, normals, Das, Dbs)

    # Calculate boundary conditions (right-hand side of system of equations)
    RHS = calc_bc_noflowthrough(Uinfs, normals)

    # Solve system of equations
    Gamma = zeros(T, self.ncells)
    solver(Gamma, G, RHS; solver_optargs...)

    # Save solution
    self.strength[:, 1] .= Gamma

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "Da", "vector", collect(eachcol(Das)), "system")
    add_field(self, "Db", "vector", collect(eachcol(Dbs)), "system")
    add_field(self, "Gamma", "scalar", view(self.strength, :, 1), "cell")
end

"""
Computes the geometric matrix (left-hand side matrix of the system of equation)
and stores it under `G`.

**ARGUMENTS**
  * `G::Array{T,2}`                     : Pre-allocated output memory.
  * `CPs::Array{T,2}`                   : Control points.
  * `normals::Array{T,2}`               : Normal associated to every CP.
  * `Das::Array{T,2}`                   : Unitary direction of semi-infinite
                                            vortex at point `a` of each
                                            trailing edge panel.
  * `Dbs::Array{T,2}`                   : Unitary direction of semi-infinite
                                            vortex at point b of each
                                            trailing edge panel.
"""
function _G_U!(self::RigidWakeBody{VortexRing, 1},
                    G::Arr1, CPs::Arr2, normals::Arr3, Das, Dbs;
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

    # Build geometric matrix: Panels
    for (pj, Gslice) in enumerate(eachcol(G))

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        U_vortexring(
                          self.grid.orggrid.nodes,           # All nodes
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


    # Add wake contributions
    TE = zeros(Int, 2)
    for (ei, (pi, nia, nib, pj, nja, njb)) in enumerate(eachcol(self.shedding)) # Iterate over wake-shedding panels

        # Fetch nodes of upper wake panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pi, lin, ndivscells, cin)

        # Indicate nodes in the upper shedding edge
        TE[1] = panel[nia]
        TE[2] = panel[nib]
        da1, da2, da3 = Das[1, ei], Das[2, ei], Das[3, ei]
        db1, db2, db3 = Dbs[1, ei], Dbs[2, ei], Dbs[3, ei]

        U_semiinfinite_horseshoe(
                          self.grid.orggrid.nodes,           # All nodes
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
             # Fetch nodes of lower wake panel
             panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                 self.grid, pj, lin, ndivscells, cin)

             # Indicate nodes in the lower shedding edge
             TE[1] = panel[nja]
             TE[2] = panel[njb]

             U_semiinfinite_horseshoe(
                               self.grid.orggrid.nodes,           # All nodes
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

    #=
    # Back-diagonal correction to avoid matrix singularity in closed geometries.
    # Seee Eq. 2.19 in Lewis, R. (1991), "Vortex Element Methods for Fluid
    # Dynamic Analysis of Engineering Systems"

    println("PROTOTYPE BACK DIAGONAL CORRECTION")
    # TODO: Remove memory allocation associated with areas
    areas = calc_areas(self)

    for m in 1:size(G, 2)

        rowi = size(G, 1) + 1 - m

        val = 0
        for n in 1:size(G, 1)
            if n != rowi

                val += G[n, m] * areas[n]

            end
        end

        # G[rowi, m] -= val/areas[rowi]
        G[rowi, m] = -val/areas[rowi]

    end
    =#
end

function _Uind!(self::RigidWakeBody{VortexRing, 1}, targets, out; optargs...)


    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    Das, Dbs = get_field(self, "Da")["field_data"], get_field(self, "Db")["field_data"]

    # Iterates over body panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Velocity of i-th panel on every target
        U_vortexring(
                            self.grid.orggrid.nodes,           # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 1],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                            optargs...
                         )
    end

    # Add wake contribution
    TE = zeros(Int, 2)
    for (ei, (pi, nia, nib, pj, nja, njb)) in enumerate(eachcol(self.shedding)) # Iterate over wake-shedding panels

        strengthi, strengthj = _get_wakestrength_mu(self, ei)

        # Fetch nodes of upper wake panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pi, lin, ndivscells, cin)

        # Indicate nodes in the upper shedding edge
        TE[1] = panel[nia]
        TE[2] = panel[nib]
        da1, da2, da3 = Das[ei]
        db1, db2, db3 = Dbs[ei]

        U_semiinfinite_horseshoe(
                          self.grid.orggrid.nodes,           # All nodes
                          TE,                                # Indices of nodes that make the shedding edge
                          da1, da2, da3,                     # Semi-infinite direction da
                          db1, db2, db3,                     # Semi-infinite direction db
                          strengthi,                         # Strength
                          targets,                           # Targets
                          out;                               # Outputs
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                          optargs...
                         )

         if pj != -1
             # Fetch nodes of lower wake panel
             panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                 self.grid, pj, lin, ndivscells, cin)

             # Indicate nodes in the lower shedding edge
             TE[1] = panel[nja]
             TE[2] = panel[njb]

             U_semiinfinite_horseshoe(
                               self.grid.orggrid.nodes,           # All nodes
                               TE,                                # Indices of nodes that make the shedding edge
                               db1, db2, db3,                     # Semi-infinite direction da (flipped in lower panel)
                               da1, da2, da3,                     # Semi-infinite direction db
                               strengthj,                         # Unitary strength
                               targets,                           # Targets
                               out;                               # Outputs
                               offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                               cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                               optargs...
                              )
         end
    end
end


function _phi!(self::RigidWakeBody{VortexRing, 1}, targets, out; optargs...)

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    Das, Dbs = get_field(self, "Da")["field_data"], get_field(self, "Db")["field_data"]

    # Iterates over body panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Potential of i-th panel on every target
        phi_constant_doublet(
                            self.grid.orggrid.nodes,           # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 1],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            optargs...
                         )
    end

    # Add wake contribution
    TE = zeros(Int, 2)
    for (ei, (pi, nia, nib, pj, nja, njb)) in enumerate(eachcol(self.shedding)) # Iterate over wake-shedding panels

        strengthi, strengthj = _get_wakestrength_mu(self, ei)

        # Fetch nodes of upper wake panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pi, lin, ndivscells, cin)

        # Indicate nodes in the upper shedding edge
        TE[1] = panel[nia]
        TE[2] = panel[nib]
        da1, da2, da3 = Das[ei]
        db1, db2, db3 = Dbs[ei]

        phi_semiinfinite_doublet(
                          self.grid.orggrid.nodes,           # All nodes
                          TE,                                # Indices of nodes that make the shedding edge
                          da1, da2, da3,                     # Semi-infinite direction da
                          db1, db2, db3,                     # Semi-infinite direction db
                          strengthi,                         # Strength
                          targets,                           # Targets
                          out;                               # Outputs
                          optargs...
                         )

         if pj != -1
             # Fetch nodes of lower wake panel
             panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                 self.grid, pj, lin, ndivscells, cin)

             # Indicate nodes in the lower shedding edge
             TE[1] = panel[nja]
             TE[2] = panel[njb]

             phi_semiinfinite_doublet(
                               self.grid.orggrid.nodes,           # All nodes
                               TE,                                # Indices of nodes that make the shedding edge
                               db1, db2, db3,                     # Semi-infinite direction da (flipped in lower panel)
                               da1, da2, da3,                     # Semi-infinite direction db
                               strengthj,                         # Unitary strength
                               targets,                           # Targets
                               out;                               # Outputs
                               optargs...
                              )
         end
    end
end


##### INTERNAL FUNCTIONS  ######################################################
function _get_wakestrength_mu(self::RigidWakeBody, i)
    strength1 = self.strength[self.shedding[1, i], 1]
    strength2 = self.shedding[4, i] != -1 ? self.strength[self.shedding[4, i], 1] : 0
    return strength1, strength2
end
function _get_wakestrength_Gamma(self::RigidWakeBody, i)
    strength1 = self.strength[self.shedding[1, i], 1]
    strength2 = self.shedding[4, i] != -1 ? self.strength[self.shedding[4, i], 1] : 0
    return strength1 - strength2
end

"""
Outputs a vtk file with the wake.
"""
function _savewake(self::RigidWakeBody, filename::String;
                    len::Number=1.0, panel::Bool=true, suffix="_wake",optargs...)

    if check_field(self, "Da")==false
        error("Requested to save wake, but Da field was not found."*
              " Run solve(...) first")
    elseif check_field(self, "Db")==false
        error("Requested to save wake, but Db field was not found."*
              " Run solve(...) first")
    end


    nodes = self.grid.orggrid.nodes
    nedges = size(self.shedding, 2)
    points = zeros(3, 4*nedges)

    # Collect points along shedding edge
    tricoor, quadcoor, lin, ndivscells, cin = gt.generate_getcellt_args(self.grid)

    for (si, (pi, nia, nib, pj)) in enumerate(eachcol(self.shedding))

        # Convert node indices from panel-local to global
        pia = gt.get_cell_t(tricoor, quadcoor, self.grid, pi, nia, lin, ndivscells, cin)
        pib = gt.get_cell_t(tricoor, quadcoor, self.grid, pi, nib, lin, ndivscells, cin)

        Da = get_fieldval(self, "Da", si)
        Db = get_fieldval(self, "Db", si)

        for i in 1:3
            # Shedding edge
            points[i, si + 0*nedges] = nodes[i, pia]
            points[i, si + 1*nedges] = nodes[i, pib]

            # End points of semi-infinite wake
            points[i, si + 2*nedges] = points[i, si + 0*nedges] + len*Da[i]
            points[i, si + 3*nedges] = points[i, si + 1*nedges] + len*Db[i]
        end

    end

    if panel    # Case that wake is saved as doublet panels

        arrpoints = collect(eachcol(points))
        str = ""

        for si in 1:2   # si==1, upper side; si==2, lower side
            cells = Array{Int64,1}[]
            cell_data = []

            # Wake panels
            for i in 0:nedges-1
                cell = [i+2*nedges, i+0*nedges, i+1*nedges, i+3*nedges]
                push!(cells, si==1 ? cell : reverse(cell))
            end

            # Strength of wake along shedding edge
            strengths = [_get_wakestrength_mu(self, i)[si] for i in 1:nedges]

            push!(cell_data, Dict( "field_name"=> "mu",
                                    "field_type"=> "scalar",
                                    "field_data"=> strengths))

            for cell in cells; cell .+= 1; end
            normals = [[gt._calc_n1(points, cell), gt._calc_n2(points, cell), gt._calc_n3(points, cell)] for cell in cells]
            for cell in cells; cell .-= 1; end

            push!(cell_data, Dict( "field_name"=> "normal",
                                    "field_type"=> "vector",
                                    "field_data"=> normals))

            # Generate VTK
            str *= gt.generateVTK(filename*suffix*(si==1 ? "_up" : "_lo"),
                                    arrpoints; cells=cells,
                                    cell_data=cell_data, optargs...)
        end

        return str

    else        # Case that wake is saved as horseshoes

        cells = Array{Int64,1}[]
        cell_data = []

        # Lines that make the wake
        for i in 0:nedges-1
            line = [i+2*nedges, i+0*nedges, i+1*nedges, i+3*nedges]
            push!(cells, line)
        end

        # Strength of wake along shedding edge
        strengths = [_get_wakestrength_Gamma(self, i) for i in 1:nedges]

        push!(cell_data, Dict( "field_name"=> "Gamma",
                                "field_type"=> "scalar",
                                "field_data"=> strengths))

        points = collect(eachcol(points))

        # Generate VTK
        return gt.generateVTK(filename*suffix, points; cells=cells,
                                cell_data=cell_data,
                                override_cell_type=4, optargs...)

    end
end

_get_Gdims(self::RigidWakeBody{VortexRing, 1}) = (self.ncells, self.ncells)
#### END OF LIFTING BODY  ######################################################
