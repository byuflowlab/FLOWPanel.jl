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
    AbstractLiftingBody

Abstract supertype for lifting-surface body models with shedding-edge and
wake-direction data in addition to the fields required by [`AbstractBody`](@ref).
"""
abstract type AbstractLiftingBody{E, N, TF, DBC} <: AbstractBody{E, N, TF, DBC} end

"""
    solve(body::AbstractLiftingBody, Uinfs, Das)

Solve a lifting-body boundary-value problem using control-point velocities
`Uinfs` and wake-direction data `Das`.
"""
function solve(self::AbstractLiftingBody, Uinfs::AbstractMatrix,
               Das::AbstractMatrix)
    error("solve(...) for body type $(typeof(self)) has not been implemented yet!")
end

##### COMMON FUNCTIONS  ########################################################
"""
    generate_loft_liftbody(bodytype, args...; bodyoptargs=(), dimsplit=2, overwrite_shedding=nothing, optargs...)

Generate a lofted lifting body, triangulate it, and build the shedding-edge map
expected by lifting-body constructors.
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
    generate_revolution_liftbody(bodytype, args...; bodyoptargs=(), optargs...)

Generate a lifting body of revolution, optionally process the intermediate
surface grid, and construct the shedding-edge connectivity required by
rigid-wake solvers.
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

    return bodytype(triang_grid, [shedding]; bodyoptargs...)
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

function extra_reset!(body::AbstractLiftingBody)
    for vte in body.velocity_te
        vte .= 0.0
    end
end

function extra_apply_freestream!(body::AbstractLiftingBody, uinf)
    for i in eachindex(body.velocity_te)
        eachcol(body.velocity_te[i]) .+= Ref(uinf)
    end
end

function _G_phi_wake!(self::AbstractLiftingBody{<:Any,<:Any,TF}, kernel, G, CPs, backend::FastMultipoleBackend; kerneloffset=1.0e-8, optargs...) where TF
    # Add wake contributions
    sheddings = 1:self.nsheddings
    chunks = collect(Iterators.partition(sheddings, max(length(sheddings) ÷ Threads.nthreads(), 3*Threads.nthreads())))
    Das = self.Das
    derivatives_switch = FastMultipole.DerivativesSwitch(true,false,false)

    for isurf in eachindex(self.shedding)
        # for chunk in chunks        # Distribute wake panel iteration among all CPU threads
        # Threads.@threads for chunk in chunks        # Distribute wake panel iteration among all CPU threads
        for chunk in chunks        # Distribute wake panel iteration among all CPU threads

            # for (ei, (pi, nia, nib, pj, nja, njb)) in enumerate(eachcol(self.shedding))
            for i_source in chunk                          # Iterate over wake-shedding panels

                pi, nia, nib, pj, nja, njb = view(self.shedding[isurf], :, i_source)

                # Fetch nodes of upper wake panel
                nodes_idx = (self.cells[1, pi], self.cells[2, pi], self.cells[3, pi])

                TE1 = nodes_idx[nia]
                TE2 = nodes_idx[nib]
                v1x = self.nodes[1, TE1]
                v1y = self.nodes[2, TE1]
                v1z = self.nodes[3, TE1]
                v2x = self.nodes[1, TE2]
                v2y = self.nodes[2, TE2]
                v2z = self.nodes[3, TE2]

                # direction of trailing semi-infinite wake
                da1, da2, da3 = Das[isurf][1, i_source+1], Das[isurf][2, i_source+1], Das[isurf][3, i_source+1]
                db1, db2, db3 = Das[isurf][1, i_source], Das[isurf][2, i_source], Das[isurf][3, i_source]
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
                    v1x = self.nodes[1, TE1]
                    v1y = self.nodes[2, TE1]
                    v1z = self.nodes[3, TE1]
                    v2x = self.nodes[1, TE2]
                    v2y = self.nodes[2, TE2]
                    v2z = self.nodes[3, TE2]

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
end

##### END OF ABSTRACT LIFTING BODY #############################################
