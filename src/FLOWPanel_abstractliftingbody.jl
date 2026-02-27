#=##############################################################################
# DESCRIPTION
    Lifting paneled body types definition.
# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Sep 2018
  * License     : MIT License
=###############################################################################

################################################################################
# ABSTRACT LIFTING BODY BODY TYPE
################################################################################
"""
  Implementations of AbstractLiftingBody are expected to have the following
  fields:
  * `U::Array{Int64,1}`                 : Indices of all panels along the upper
                                          side of the trailing edge.
  * `L::Array{Int64,1}`                 : Indices of all panels along the lower
                                          side of the trailing edge.
  * `ncellsTE::Int64`                   : Number of cells along trailing edge
  * `nnodesTE::Int64`                   : Number of nodes along trailing edge

  NOTE: U and L are assumed to have the same number of points.

  in addition to the same properties and functions expected from an AbstractBody
  implementation. The following functions also need to be implemented.

  ```julia

  # Impose boundary conditions to solve for element strengths
  function solve(self::AbstractLiftingBody, Uinfs::Array{<:Real, 2},
                                            D::Array{<:Real, 2}, args...)
    .
    .
    .
  end

  # Outputs a vtk file with the wake
  function _savewake(self::AbstractLiftingBody, filename::String;
                                  len::Real=1.0, upper::Bool=true, optargs...)
    .
    .
    .
  end
  ```
"""
abstract type AbstractLiftingBody{E, N, TF} <: AbstractBody{E, N, TF} end

"""
    `solve(body::AbstractBody, Uinfs::Array{<:Real, 2})`

Impose boundary conditions to solve for element strengths. `Uinds[:, i]` is the
velocity at the i-th control point used in the boundary condition.

`Das[:, i]` is the unitary vector pointing in the direction that the
semi-infinite wake is shed from the first node of the i-th shedding edge, while
`Dbs[:, i]` is for the second node.
NOTE: These directions are expected to point from the node out to infinite.
"""
function solve(self::AbstractLiftingBody, Uinfs::AbstractMatrix,
               Das::AbstractMatrix, Dbs::AbstractMatrix)
    error("solve(...) for body type $(typeof(self)) has not been implemented yet!")
end

##### COMMON FUNCTIONS  ########################################################
"""
  `generate_loft_liftbody(bodytype::Type{<:AbstractLiftingBody}, args...; optargs...)`
Generates a lofted lifting body of type `bodytype`. See documentation of
`GeometricTools.generate_loft` for a description of the arguments of this
function.
"""
function generate_loft_liftbody(bodytype::Type{B}, args...;
                                bodyoptargs=(), dimsplit::Int=2,
                                overwrite_shedding=nothing,
                                optargs...
                               ) where {B<:AbstractBody}
    # Lofts the surface geometry
    grid = gt.generate_loft(args...; optargs...)

    # Splits the quadrialateral panels into triangles
    # dimsplit = 2              # Dimension along which to split
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)

    ndivs = gt.get_ndivscells(triang_grid)              # Cells in each dimension
    U = [ Base._sub2ind(ndivs, ndivs[1]-1, i) for i in 1:ndivs[2] ] # Upper LE cells
    L = [ Base._sub2ind(ndivs, 2, i) for i in 1:ndivs[2] ]          # Lower LE cells

    if isnothing(overwrite_shedding)

        nedges = length(U)
        shedding = zeros(Int, 6, nedges)
        for (ei, (u, l)) in enumerate(zip(U, L))
            shedding[1, ei] = u
            shedding[2, ei] = 3
            shedding[3, ei] = 2

            shedding[4, ei] = l
            shedding[5, ei] = 3
            shedding[6, ei] = 2
        end

    else
        shedding = overwrite_shedding
    end

    @show bodytype
    if bodytype <: AbstractLiftingBody
        return bodytype(triang_grid, shedding; bodyoptargs...)
    elseif bodytype <: NonLiftingBody
        return bodytype(triang_grid; bodyoptargs...)
    else
        error("Body type $(bodytype) is not a lifting or non-lifting body.")
    end
end

"""
  `generate_revolution_liftbody(bodytype::Type{<:AbstractLiftingBody}, args...; optargs...)`
Generates a lifting body type `bodytype` of a body of revolution. See
documentation of `GeometricTools.surface_revolution` for a description of the
arguments of this function.

> **NOTE:** For a complete revolution generating a water-tight grid, the system
of equations solved in the body type `RigidWakeBody{VortexRing}` becomes
singular (this is because in the absence of an open edge, the solution
depends only in the difference between adjacent panels rather than the strengths
themselves), leading to vortex ring strengths that are astronomically large. To
avoid this situation, use `RigidWakeBody{VortexRing, 1}`, which indicates the
solver to center the solution around 0.
"""
function generate_revolution_liftbody(bodytype::Type{B}, args...;
                                                  bodyoptargs=(),
                                                  gridprocessing=nothing,
                                                  dimsplit::Int=1,
                                                  # loop_dim::Int=1,
                                                  loop_dim::Int=2,
                                                  axis_angle=270,
                                                  overwrite_shedding=nothing,
                                                  closed_contour=true,
                                                  optargs...
                                      ) where {B<:AbstractLiftingBody}
    # Revolves the geometry
    grid = gt.surface_revolution(args...; loop_dim=loop_dim,
                                            axis_angle=axis_angle, optargs...)

    # Intermediate processing of grid: rotate to align centerline with x-axis
    if gridprocessing==nothing
        Oaxis = gt.rotation_matrix2(0, 0, 90)
        O = zeros(3)
        gt.lintransform!(grid, Oaxis, O)

    # User-defined intermediate processing of grid
    else
        gridprocessing(grid)
    end

    # Splits the quadrialateral panels into triangles
    # dimsplit = 2              # Dimension along which to split
    triang_grid = gt.GridTriangleSurface(grid, dimsplit)

    if isnothing(overwrite_shedding)

        ndivs = gt.get_ndivscells(triang_grid)                 # Cells in each dimension
        U = [ Base._sub2ind(ndivs, ndivs[1]-1, i) for i in 1:ndivs[2] ] # Upper LE cells
        L = [ Base._sub2ind(ndivs, 2, i) for i in 1:ndivs[2] ]          # Lower LE cells

        nedges = length(U)
        shedding = zeros(Int, 6, nedges)
        for (ei, (u, l)) in enumerate(zip(U, L))
            shedding[1, ei] = u
            shedding[2, ei] = 3
            shedding[3, ei] = 2

            shedding[4, ei] = closed_contour ? l : -1
            shedding[5, ei] = 3
            shedding[6, ei] = 2
        end
    else
        shedding = overwrite_shedding
    end

    return bodytype(triang_grid, shedding; bodyoptargs...)
end

##### COMMON INTERNAL FUNCTIONS  ###############################################
"""
Checks correction definition of trailing edge
"""
function _checkTE(grid, shedding::Array{Int, 2}; tol=1e1*eps())

    nodes = grid._nodes

    # Correct number of inputs
    if size(shedding, 1) != 6
        return false
    end

    tricoor, quadcoor, lin, ndivscells, cin = gt.generate_getcellt_args(grid)

    # Check that node position along edge of each side are coincident
    for (pi, nia, nib, pj, nja, njb) in eachcol(shedding)

        if pj != -1

            # Convert node indices from panel-local to global
            pia = gt.get_cell_t(tricoor, quadcoor, grid, pi, nia, lin, ndivscells, cin)
            pib = gt.get_cell_t(tricoor, quadcoor, grid, pi, nib, lin, ndivscells, cin)
            pja = gt.get_cell_t(tricoor, quadcoor, grid, pj, nja, lin, ndivscells, cin)
            pjb = gt.get_cell_t(tricoor, quadcoor, grid, pj, njb, lin, ndivscells, cin)

            for i in 1:3
                if abs(nodes[i, pia] - nodes[i, pjb]) > tol
                    println("rabbit2\t$(abs(nodes[i, pia] - nodes[i, pjb]))")
                    return false
                elseif abs(nodes[i, pib] - nodes[i, pja]) > tol
                    println("rabbit3\t$(abs(nodes[i, pib] - nodes[i, pja]))")
                    return false
                end
            end

        end

    end

    return true
end

function _G_phi_wake!(self::AbstractLiftingBody{<:Any,<:Any,TF}, kernel, G, CPs, backend::FastMultipoleBackend; kerneloffset=1.0e-8, optargs...) where TF
    # Add wake contributions
    sheddings = 1:self.nsheddings
    chunks = collect(Iterators.partition(sheddings, max(length(sheddings) ÷ Threads.nthreads(), 3*Threads.nthreads())))
    Das, Dbs = self.Das, self.Dbs
    derivatives_switch = FastMultipole.DerivativesSwitch(true,false,false)

    # for chunk in chunks        # Distribute wake panel iteration among all CPU threads
    # Threads.@threads for chunk in chunks        # Distribute wake panel iteration among all CPU threads
    for chunk in chunks        # Distribute wake panel iteration among all CPU threads

        # for (ei, (pi, nia, nib, pj, nja, njb)) in enumerate(eachcol(self.shedding))
        for i_source in chunk                          # Iterate over wake-shedding panels

            pi, nia, nib, pj, nja, njb = view(self.shedding, :, i_source)

            # Fetch nodes of upper wake panel
            nodes_idx = (self.cells[1, pi], self.cells[2, pi], self.cells[3, pi])

            TE1 = nodes_idx[nia]
            TE2 = nodes_idx[nib]
            v1x = self.grid._nodes[1, TE1]
            v1y = self.grid._nodes[2, TE1]
            v1z = self.grid._nodes[3, TE1]
            v2x = self.grid._nodes[1, TE2]
            v2y = self.grid._nodes[2, TE2]
            v2z = self.grid._nodes[3, TE2]

            # direction of trailing semi-infinite wake
            da1, da2, da3 = Das[1, i_source], Das[2, i_source], Das[3, i_source]
            db1, db2, db3 = Dbs[1, i_source], Dbs[2, i_source], Dbs[3, i_source]
            @assert isapprox(da1, db1) && isapprox(da2, db2) && isapprox(da3, db3) "Inconsistent wake directions in _G_phi_wake!"

            for i_target in axes(CPs, 2)
                # get target
                tx, ty, tz = CPs[1, i_target], CPs[2, i_target], CPs[3, i_target]
                target = FastMultipole.StaticArrays.SVector{3,TF}(tx, ty, tz)

                # compute influence
                phi, _ = induced_semiinfinite(target, kernel, v1x, v1y, v1z, v2x, v2y, v2z, da1, da2, da3, 1.0, derivatives_switch; kerneloffset)

                # update G
                G[i_target, pi] += phi
            end

            # lower wake panel (if it exists)
            if pj != -1
                # Fetch nodes of lower wake panel
                nodes_idx = (self.cells[1, pj], self.cells[2, pj], self.cells[3, pj])

                TE1 = nodes_idx[nja]
                TE2 = nodes_idx[njb]
                v1x = self.grid._nodes[1, TE1]
                v1y = self.grid._nodes[2, TE1]
                v1z = self.grid._nodes[3, TE1]
                v2x = self.grid._nodes[1, TE2]
                v2y = self.grid._nodes[2, TE2]
                v2z = self.grid._nodes[3, TE2]

                for i_target in axes(CPs, 2)
                    # get target
                    tx, ty, tz = CPs[1, i_target], CPs[2, i_target], CPs[3, i_target]
                    target = FastMultipole.StaticArrays.SVector{3,TF}(tx, ty, tz)

                    # compute influence
                    phi, _ = induced_semiinfinite(target, kernel, v1x, v1y, v1z, v2x, v2y, v2z, da1, da2, da3, 1.0, derivatives_switch; kerneloffset)

                    # update G
                    G[i_target, pj] += phi
                end

            end
        end
    end
end

##### END OF ABSTRACT LIFTING BODY #############################################
