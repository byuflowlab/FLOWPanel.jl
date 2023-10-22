#=##############################################################################
# DESCRIPTION
    Definition of methods for post-processing solver results.

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Oct 2022
  * License     : MIT License
=###############################################################################


################################################################################
# VELOCITY FIELDS
################################################################################
"""
    calcfield_U!(out::Matrix,
                 sourcebody::AbstractBody, targetbody::AbstractBody,
                 controlpoints::Matrix, Uinfs::Matrix; fieldname="U")

Calculate the velocity induced by `sourcebody` on `controlpoints` and save
it as a field of name `fieldname` under `targetbody`. The field includes the
freestream velocity `Uinfs`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_U!(out::Arr1, sourcebody::AbstractBody, targetbody::AbstractBody,
                        controlpoints::Arr2, Uinfs::Arr3;
                        fieldname="U", addfield=true, optargs...
                        ) where {   Arr1<:AbstractArray{<:Number,2},
                                    Arr2<:AbstractArray{<:Number,2},
                                    Arr3<:AbstractArray{<:Number,2}}

    # ERROR CASES
    if check_solved(sourcebody)==false
        error("Source body hasn't been solved yet."*
              " Please call `solve(...)` function first.")
    elseif size(controlpoints, 1)!=3 || size(controlpoints, 2)!=targetbody.ncells
        error("Invalid `controlpoints` matrix."*
              " Expected size $((3, targetbody.ncells)); got $(size(controlpoints)).")
    elseif size(Uinfs, 1)!=3 || size(Uinfs, 2)!=targetbody.ncells
        error("Invalid `Uinfs` matrix."*
              " Expected size $((3, targetbody.ncells)); got $(size(Uinfs)).")
    elseif size(out, 1)!=3 || size(out, 2)!=targetbody.ncells
        error("Invalid `out` matrix."*
              " Expected size $((3, targetbody.ncells)); got $(size(out)).")
    end

    # Add freestream
    out .+= Uinfs

    # Add induced velocity at each control point
    Uind!(sourcebody, controlpoints, out; optargs...)

    # Save field in body
    if addfield
        add_field(targetbody, fieldname, "vector", eachcol(out), "cell")
    end

    return out
end

"""
    calcfield_U!(out::Matrix,
                    sourcebody::AbstractBody, targetbody::AbstractBody;
                    offset=nothing, characteristiclength=nothing, optargs...)

Calculate the velocity induced by `sourcebody` on control points computed
using `offset` and `characteristiclength`, and save it as a field in
`targetbody`. The field includes the freestream velocity stored as field
`\"Uinf\"` in `targetbody`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_U!(out::Arr,
                        sourcebody::AbstractBody, targetbody::AbstractBody;
                        offset=nothing, characteristiclength=nothing,
                        optargs...
                        ) where {Arr<:AbstractArray{<:Number,2}}

    @assert check_field(targetbody, "Uinf") ""*
        "Target body doesn't have freestream field `\"Uinf\"`."*
        " Please call `add_field(targetbody, \"Uinf\", ...)` first."

    Uinfs = hcat(get_field(targetbody, "Uinf")["field_data"]...)

    # Optional arguments for calc_controlpoints
    cp_optargs = (off=offset, characteristiclength=characteristiclength)
    cp_optargs = ((key, val) for (key, val) in pairs(cp_optargs) if val!=nothing)

    # Calculate control points
    normals = calc_normals(targetbody)
    controlpoints = calc_controlpoints(targetbody, normals; cp_optargs...)

    # Calculate field on control points
    calcfield_U!(out, sourcebody, targetbody, controlpoints, Uinfs; optargs...)
end

"""
    calcfield_U(args...; optargs...)

Similar to [`calcfield_U!`](@ref) but without in-place calculation (`out` is not
needed).
"""
function calcfield_U(sourcebody, targetbody, args...; optargs...)
    out = zeros(3, targetbody.ncells)
    return calcfield_U!(out, sourcebody, targetbody, args...; optargs...)
end

"""
    calcfield_Uoff!(args...; optargs...) = calcfield_U(args...; optargs..., fieldname="Uoff")

See documentation of `calcfield_U!(...)`.
"""
calcfield_Uoff!(args...; optargs...) = calcfield_U!(args...; optargs..., fieldname="Uoff")
calcfield_Uoff(args...; optargs...) = calcfield_U(args...; optargs..., fieldname="Uoff")




"""
    calcfield_Ugradmu_cell!(out::Matrix, body::AbstractBody;
                            fieldname="Ugradmu")

Calculate the surface velocity on `body` due to changes in the constant
doublet strength using the Green-Gauss method and save it as a field of name `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).

TODO: Avoid the large gradient at the trailing edge recognizing the trailing
        edge and omiting the neighbor that should be the wake.
"""
function calcfield_Ugradmu_cell!(out::AbstractMatrix, body::AbstractBody,
                                areas::AbstractVector,
                                normals::AbstractMatrix,
                                controlpoints::AbstractMatrix;
                                fieldname="Ugradmu", addfield=true,
                                Gammai=1,
                                maxgrad=Inf,
                                )
    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==body.ncells ""*
        "Invalid `out` matrix."*
        " Expected size $((3, body.ncells)); got $(size(out))."
    @assert length(areas)==body.ncells ""*
        "Invalid `areas` vector."*
        " Expected length $(body.ncells); got $(length(areas))."
    @assert size(normals, 1)==3 && size(normals, 2)==body.ncells ""*
        "Invalid `normals` matrix."*
        " Expected size $((3, body.ncells)); got $(size(normals))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==body.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, body.ncells)); got $(size(controlpoints))."

    # Fetch data
    Gammas = view(body.strength, :, Gammai)
    nodes = body.grid.orggrid.nodes

    # Pre-allocate memory
    (tri_out, tricoor, quadcoor,
        quad_out, lin, ndivscells, cin) = gt.generate_getcellt_args!(body.grid)

    ndivscellsc = Tuple(collect( 1:(d != 0 ? d : 1) for d in body.grid._ndivscells))
    linc = LinearIndices(ndivscellsc)
    cinc = CartesianIndices(ndivscellsc)

    ncoor = ones(Int, 3)                # Stores coordinates of neighbor here

    # Iterate over cells
    for ci in 1:body.ncells             # Iterate over linear indexing
        ccoor = cinc[ci]                # Cartesian indexing of this cell

        # Fetch the cell
        panel = gt.get_cell_t!(tri_out, quadcoor, quad_out,
                        body.grid, collect(Tuple(ccoor)), lin, ndivscells)

        for ni in 1:3                   # Iterate over neighbors

            # Obtain coordinates of ni-th neighbor
            ncoor = gt.neighbor(body.grid, ni, ci; preserveEdge=true)

            if ncoor[1] != 0
                # Linear indexing of this neighbor
                nlin = linc[ncoor...]

                ei, ej = ni, ni%3 + 1

                # r = pj - pi
                r1 = nodes[1, tri_out[ej]] - nodes[1, tri_out[ei]]
                r2 = nodes[2, tri_out[ej]] - nodes[2, tri_out[ei]]
                r3 = nodes[3, tri_out[ej]] - nodes[3, tri_out[ei]]

                # d = r⨉n / |r⨉n| (normal to edge)
                d1 = r2*normals[3, ci] - r3*normals[2, ci]
                d2 = r3*normals[1, ci] - r1*normals[3, ci]
                d3 = r1*normals[2, ci] - r2*normals[1, ci]

                # # d = (cpj - cpi) / |cpj - cpi| (centroid to centroid)
                # d1 = controlpoints[1, nlin] - controlpoints[1, ci]
                # d2 = controlpoints[2, nlin] - controlpoints[2, ci]
                # d3 = controlpoints[3, nlin] - controlpoints[3, ci]

                dmag = sqrt(d1^2 + d2^2 + d3^2)
                d1 /= dmag
                d2 /= dmag
                d3 /= dmag

                # Use Green-Gauss method to compute gradient of circulation
                # where the interpolated gamma at each face (edge) is used

                # Compute vector from one edge vertex to cell-center
                vecMain = nodes[1:3, tri_out[ei]] - controlpoints[1:3, ci]
                vecNear = nodes[1:3, tri_out[ei]] - controlpoints[1:3, nlin]

                # Compute approx. distance of cell-center to edge
                # Common denominator has been cancelled out
                dMain = norm(cross(vecMain, [r1, r2, r3]))
                dNear = norm(cross(vecNear, [r1, r2, r3]))

                # r = [r1, r2, r3]
                # rhat = r/norm(r)
                # dMain = norm( vecMain - dot(vecMain, rhat)*rhat )
                # dNear = norm( vecNear - dot(vecNear, rhat)*rhat )

                # Compute inverse distance weighted interpolation factor
                f = dNear/(dMain + dNear)

                # Override interpolation factor to 0.5 for debugging
                # This is just averaging between gamma
                # f = 0.5

                # Compute face gamma
                faceGamma = f*Gammas[ci] + (1.0-f)*Gammas[nlin]

                # Invert direction of vector if normals point inward
                sgn = body.CPoffset==0 ? 1 : sign(body.CPoffset)

                # Add contribution from face gamma
                mag = faceGamma * sqrt(r1^2 + r2^2 + r3^2) / areas[ci]

                # if abs(mag) < maxgrad
                    out[1, ci] -= 0.5 * sgn * d1 * mag
                    out[2, ci] -= 0.5 * sgn * d2 * mag
                    out[3, ci] -= 0.5 * sgn * d3 * mag
                # end

            end
        end

    end

    # Quick and dirty fix to omit the high gradient at the trailing edge
    for ci in 1:body.ncells
        if norm(view(out, :, ci)) >= maxgrad
            out[:, ci] *= 0
        end
    end

    # Save field in body
    if addfield
        add_field(body, fieldname, "vector", eachcol(out), "cell")
    end

    return out
end

function calcfield_Ugradmu_cell!(out::AbstractMatrix, body::RigidWakeBody,
                                areas::AbstractVector,
                                normals::AbstractMatrix,
                                controlpoints::AbstractMatrix;
                                fieldname="Ugradmu", addfield=true,
                                Gammai=1,
                                maxgrad=Inf,
                                smoothPass=0, smoothRows=[0]
                                )
    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==body.ncells ""*
        "Invalid `out` matrix."*
        " Expected size $((3, body.ncells)); got $(size(out))."
    @assert length(areas)==body.ncells ""*
        "Invalid `areas` vector."*
        " Expected length $(body.ncells); got $(length(areas))."
    @assert size(normals, 1)==3 && size(normals, 2)==body.ncells ""*
        "Invalid `normals` matrix."*
        " Expected size $((3, body.ncells)); got $(size(normals))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==body.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, body.ncells)); got $(size(controlpoints))."

    # Fetch data
    Gammas = view(body.strength, :, Gammai)
    nodes = body.grid.orggrid.nodes

    # Pre-allocate memory
    (tri_out, tricoor, quadcoor,
        quad_out, lin, ndivscells, cin) = gt.generate_getcellt_args!(body.grid)

    ndivscellsc = Tuple(collect( 1:(d != 0 ? d : 1) for d in body.grid._ndivscells))
    linc = LinearIndices(ndivscellsc)
    cinc = CartesianIndices(ndivscellsc)

    ncoor = ones(Int, 3)                # Stores coordinates of neighbor here

    # Iterate over cells
    for ci in 1:body.ncells             # Iterate over linear indexing
        ccoor = cinc[ci]                # Cartesian indexing of this cell


        # Fetch the cell
        panel = gt.get_cell_t!(tri_out, quadcoor, quad_out,
                        body.grid, collect(Tuple(ccoor)), lin, ndivscells)

        for ni in 1:3                   # Iterate over neighbors

            # Obtain coordinates of ni-th neighbor
            ncoor = gt.neighbor(body.grid, ni, ci; preserveEdge=true)

            if ncoor[1] != 0
                # Linear indexing of this neighbor
                nlin = linc[ncoor...]

                ei, ej = ni, ni%3 + 1

                # r = pj - pi
                r1 = nodes[1, tri_out[ej]] - nodes[1, tri_out[ei]]
                r2 = nodes[2, tri_out[ej]] - nodes[2, tri_out[ei]]
                r3 = nodes[3, tri_out[ej]] - nodes[3, tri_out[ei]]

                # d = r⨉n / |r⨉n| (normal to edge)
                d1 = r2*normals[3, ci] - r3*normals[2, ci]
                d2 = r3*normals[1, ci] - r1*normals[3, ci]
                d3 = r1*normals[2, ci] - r2*normals[1, ci]

                # # d = (cpj - cpi) / |cpj - cpi| (centroid to centroid)
                # d1 = controlpoints[1, nlin] - controlpoints[1, ci]
                # d2 = controlpoints[2, nlin] - controlpoints[2, ci]
                # d3 = controlpoints[3, nlin] - controlpoints[3, ci]

                dmag = sqrt(d1^2 + d2^2 + d3^2)
                d1 /= dmag
                d2 /= dmag
                d3 /= dmag

                # Use Green-Gauss method to compute gradient of circulation
                # where the interpolated gamma at each face (edge) is used

                # Compute vector from one edge vertex to cell-center
                vecMain = nodes[1:3, tri_out[ei]] - controlpoints[1:3, ci]
                vecNear = nodes[1:3, tri_out[ei]] - controlpoints[1:3, nlin]

                # Compute approx. distance of cell-center to edge
                # Common denominator has been cancelled out
                dMain = norm(cross(vecMain, [r1, r2, r3]))
                dNear = norm(cross(vecNear, [r1, r2, r3]))

                # r = [r1, r2, r3]
                # rhat = r/norm(r)
                # dMain = norm( vecMain - dot(vecMain, rhat)*rhat )
                # dNear = norm( vecNear - dot(vecNear, rhat)*rhat )

                # Compute inverse distance weighted interpolation factor
                f = dNear/(dMain + dNear)

                # Override interpolation factor to 0.5 for debugging
                # This is just averaging between gamma
                # f = 0.5

                # Compute face gamma
                faceGamma = f*Gammas[ci] + (1.0-f)*Gammas[nlin]

                # Invert direction of vector if normals point inward
                sgn = body.CPoffset==0 ? 1 : sign(body.CPoffset)

                # Add contribution from face gamma
                mag = faceGamma * sqrt(r1^2 + r2^2 + r3^2) / areas[ci]

                # if abs(mag) < maxgrad
                    out[1, ci] -= 0.5 * sgn * d1 * mag
                    out[2, ci] -= 0.5 * sgn * d2 * mag
                    out[3, ci] -= 0.5 * sgn * d3 * mag
                # end
            end

        end

    end

    # Iterate over TE cells
    for (pi, nia, nib, pj, nja, njb) in eachcol(body.shedding)

        sides = pj!=-1 ? ((pi, nia, nib), (pj, nja, njb)) : ((pi, nia, nib),)


        # for (ci, ei, ej) in sides               # Iterate over both sides
            # # Identify neighbor index where the wake is
            # ni =    (ei==1 && ej==2) || (ei==2 && ej==1) ? 1 :
            #         (ei==2 && ej==3) || (ei==3 && ej==2) ? 2 :
            #         (ei==3 && ej==1) || (ei==1 && ej==3) ? 3 :
            #         error("Logic error: Invalid trailing edge!")

        for (ci, _, _) in sides               # Iterate over both sides

            for ni in 1:3                   # Iterate over neighbors

                ccoor = cinc[ci]                # Cartesian indexing of this cell

                # Obtain coordinates of ni-th neighbor
                ncoor = gt.neighbor(body.grid, ni, ci; preserveEdge=true)

                if ncoor[1] != 0
                    # Linear indexing of this neighbor
                    nlin = linc[ncoor...]

                    ei, ej = ni, ni%3 + 1

                    # Fetch the cell
                    panel = gt.get_cell_t!(tri_out, quadcoor, quad_out,
                                    body.grid, collect(Tuple(ccoor)), lin, ndivscells)

                    # Obtain coordinates of ni-th neighbor
                    ncoor = gt.neighbor(body.grid, ni, ci; preserveEdge=true)

                    if ncoor[1] != 0
                        # Linear indexing of this neighbor
                        nlin = linc[ncoor...]

                        # r = pj - pi
                        r1 = nodes[1, tri_out[ej]] - nodes[1, tri_out[ei]]
                        r2 = nodes[2, tri_out[ej]] - nodes[2, tri_out[ei]]
                        r3 = nodes[3, tri_out[ej]] - nodes[3, tri_out[ei]]

                        # d = r⨉n / |r⨉n| (normal to edge)
                        d1 = r2*normals[3, ci] - r3*normals[2, ci]
                        d2 = r3*normals[1, ci] - r1*normals[3, ci]
                        d3 = r1*normals[2, ci] - r2*normals[1, ci]

                        # # d = (cpj - cpi) / |cpj - cpi| (centroid to centroid)
                        # d1 = controlpoints[1, nlin] - controlpoints[1, ci]
                        # d2 = controlpoints[2, nlin] - controlpoints[2, ci]
                        # d3 = controlpoints[3, nlin] - controlpoints[3, ci]

                        dmag = sqrt(d1^2 + d2^2 + d3^2)
                        d1 /= dmag
                        d2 /= dmag
                        d3 /= dmag

                        # Use Green-Gauss method to compute gradient of circulation
                        # where the interpolated gamma at each face (edge) is used

                        # Compute vector from one edge vertex to cell-center
                        vecMain = nodes[1:3, tri_out[ei]] - controlpoints[1:3, ci]
                        vecNear = nodes[1:3, tri_out[ei]] - controlpoints[1:3, nlin]

                        # Compute approx. distance of cell-center to edge
                        # Common denominator has been cancelled out
                        dMain = norm(cross(vecMain, [r1, r2, r3]))
                        dNear = norm(cross(vecNear, [r1, r2, r3]))

                        # r = [r1, r2, r3]
                        # rhat = r/norm(r)
                        # dMain = norm( vecMain - dot(vecMain, rhat)*rhat )
                        # dNear = norm( vecNear - dot(vecNear, rhat)*rhat )

                        # Compute inverse distance weighted interpolation factor
                        f = dNear/(dMain + dNear)

                        # Override interpolation factor to 0.5 for debugging
                        # This is just averaging between gamma
                        # f = 0.5

                        # Compute face gamma
                        faceGamma = f*Gammas[ci] + (1.0-f)*Gammas[nlin]

                        # Invert direction of vector if normals point inward
                        sgn = body.CPoffset==0 ? 1 : sign(body.CPoffset)

                        # Add contribution from face gamma
                        mag = faceGamma * sqrt(r1^2 + r2^2 + r3^2) / areas[ci]

                        # Cancels out the neighbor where the wake is supposed to be
                        # if abs(mag) < maxgrad
                            out[1, ci] += 0.5 * sgn * d1 * mag
                            out[2, ci] += 0.5 * sgn * d2 * mag
                            out[3, ci] += 0.5 * sgn * d3 * mag
                        # end
                    end
                end
            end

        end

    end

    # Smoothen gradient of edge cells AFTER computation of all gradients
    if smoothRows[1] != 0 && smoothPass != 0
        if body.grid.orggrid.loop_dim == 2 && body.grid.dimsplit == 1
            for pass = 1:smoothPass
                for i in smoothRows, j in 1:body.grid.orggrid.NDIVS[2]
                    ci = linc[i, j, 1]

                    out[1:3, ci] .= 0.0
                    denom = 0

                    for ni in 1:3  # Iterate over neighbors

                        # Obtain coordinates of ni-th neighbor
                        ncoor = gt.neighbor(body.grid, ni, ci; preserveEdge=true)

                        if ncoor[1] != 0
                            denom += 1
                            # Linear indexing of this neighbor
                            nlin = linc[ncoor...]

                            out[1, ci] += out[1, nlin]
                            out[2, ci] += out[2, nlin]
                            out[3, ci] += out[3, nlin]
                        end
                    end
                    # Average of the gradient of neighboring cells
                    out[1:3, ci] = out[1:3, ci] ./ denom
                end
            end
        end
    end

    # Quick and dirty fix to omit the high gradient at the trailing edge
    for ci in 1:body.ncells
        if norm(view(out, :, ci)) >= maxgrad
            out[:, ci] *= 0
        end
    end

    # Save field in body
    if addfield
        add_field(body, fieldname, "vector", eachcol(out), "cell")
    end

    return out
end

function calcfield_Ugradmu_cell!(out::AbstractMatrix, mbody::MultiBody,
                                areas::AbstractVector,
                                normals::AbstractMatrix,
                                controlpoints::AbstractMatrix, args...;
                                fieldname="Ugradmu", addfield=true,
                                optargs...
                                )

    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==mbody.ncells ""*
        "Invalid `out` matrix."*
        " Expected size $((3, mbody.ncells)); got $(size(out))."
    @assert length(areas)==mbody.ncells ""*
        "Invalid `areas` vector."*
        " Expected length $(mbody.ncells); got $(length(areas))."
    @assert size(normals, 1)==3 && size(normals, 2)==mbody.ncells ""*
        "Invalid `normals` matrix."*
        " Expected size $((3, mbody.ncells)); got $(size(normals))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==mbody.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, mbody.ncells)); got $(size(controlpoints))."

    counter = 0

    for body in mbody.bodies

        offset = body.ncells
        thisout = view(out, 1:3, (1:offset) .+ counter)
        thisareas = view(areas, (1:offset) .+ counter)
        thisnormals = view(normals, 1:3, (1:offset) .+ counter)
        thiscontrolpoints = view(controlpoints, 1:3, (1:offset) .+ counter)

        calcfield_Ugradmu_cell!(thisout, body, thisareas, thisnormals,
                                thiscontrolpoints, args...;
                                fieldname=fieldname, addfield=addfield,
                                optargs...)
        counter += offset
    end

    if addfield && !(fieldname in mbody.fields)
        push!(mbody.fields, fieldname)
    end

    return out
end

function calcfield_Ugradmu_cell!(out::AbstractMatrix, body::AbstractBody; optargs...)

    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals)
    areas = calc_areas(body)

    return calcfield_Ugradmu_cell!(out, body, areas, normals, controlpoints; optargs...)
end

"""
    calcfield_Ugradmu_cell(body::AbstractBody; fieldname="Ugradmu")

Similar to [`calcfield_Ugradmu_cell!`](@ref) but without in-place calculation
(`out` is not needed).
"""
function calcfield_Ugradmu_cell(body::AbstractBody; optargs...)
    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals)
    areas = calc_areas(body)

    out = zeros(3, body.ncells)
    calcfield_Ugradmu_cell!(out, body, areas, normals, controlpoints; optargs...)
    return out
end

################################################################################
# GRADIENT COMPUTATION USING NODAL VALUES
################################################################################

function calcfield_Ugradmu!(out::AbstractMatrix, body::AbstractBody,
                                    areas::AbstractVector, normals::AbstractMatrix,
                                    controlpoints::AbstractMatrix;
                                    fieldname="Ugradmu", addfield=true, Gammai=1,
                                    sharpTE=false)

    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==body.ncells ""*
        "Invalid `out` matrix."*
        " Expected size $((3, body.ncells)); got $(size(out))."
    @assert length(areas)==body.ncells ""*
        "Invalid `areas` vector."*
        " Expected length $(body.ncells); got $(length(areas))."
    @assert size(normals, 1)==3 && size(normals, 2)==body.ncells ""*
        "Invalid `normals` matrix."*
        " Expected size $((3, body.ncells)); got $(size(normals))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==body.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, body.ncells)); got $(size(controlpoints))."

    # Fetch data
    nodes = body.grid.orggrid.nodes

    # This algorithm might be inaccurate if the body has multiple types of elements
    Gammas = view(body.strength, :, Gammai)

    # Compute nodal data for each node
    nodal_data = gt.get_nodal_data(body.grid, Gammas; areas=areas)

    # Pre-allocate required arrays
    A = Array{Float64}(undef, 3, 3)
    t0 = zeros(2)
    t1 = zeros(2)
    t2 = zeros(2)
    t3 = zeros(2)
    e1 = zeros(3)
    e2 = zeros(3)
    grad = zeros(3)

    # Use CPoffset as a flag to know if the normals are flipped into the
    # body. If that's the case, then it flips the sign of the nodal quantity
    # to have the effect of implicitly flipping the normals out
    gamma_sign = (-1)^(body.CPoffset<0)

    # Compute cell-based gradient for each cell
    for i = 1:prod(body.grid._ndivscells[1:2])
        # Convert cell vertices to a local x,y coordinate frame
        vtx = gt.get_cell(body.grid, i)
        gt.project_3d_2d!(t2, t3, e1, e2,
                          nodes[:, vtx[1]],
                          nodes[:, vtx[2]],
                          nodes[:, vtx[3]])

        # The (x, y) coordinate of t1 is always at origin
        t0 = @. (t1 + t2 + t3) / 3.0

        # Get gradient of plane formed by three scalar values at vertices
        gt.get_tri_gradient!(grad, t1, t2, t3, t0, e1, e2, A, gamma_sign*nodal_data[vtx])

        # Transform slopes back to global coordinate system
        out[:, i] = @. (grad[2]*e1 + grad[3]*e2)
    end

    # If it's a sharp TE, do not use contribution from the other side of the mesh
    # while converting from cell to nodal data
    # NOTE: This automatically calculates the TE with some assumptions which
    # might not generally hold true. Redo this part using the information in
    # body.sheddings (in the case of a lifting over) to iterate over the TE
    # instead
    if sharpTE && body.grid.orggrid.loop_dim == 1
        if body.grid.dimsplit == 1
            # Compute TE node indices
            lin_node = LinearIndices(body.grid._ndivsnodes)
            TE_idx = lin_node[1, :, 1]

            # Compute trailing cell indices that share vertices with TE nodes
            lin_cell = LinearIndices(body.grid._ndivscells[1:2])
            cells_U = vec(lin_cell[end-1:end, :])
            cells_L = vec(lin_cell[1:2, :])

            nodal_data_U, nodal_data_L = gt.get_nodal_data_TEcells(body.grid, Gammas,
                                                                   TE_idx,
                                                                   cells_U, cells_L;
                                                                   areas=areas)

            # Overwrite TE node values for cells on upper side
            nodal_data[TE_idx] .= nodal_data_U
            for i in cells_U
                # Convert cell vertices to a local x,y coordinate frame
                vtx = gt.get_cell(body.grid, i)
                gt.project_3d_2d!(t2, t3, e1, e2,
                                  nodes[:, vtx[1]],
                                  nodes[:, vtx[2]],
                                  nodes[:, vtx[3]])

                # The (x, y) coordinate of t1 is always at origin
                t0 = @. (t1 + t2 + t3) / 3.0

                # Get gradient of plane formed by three scalar values at vertices
                gt.get_tri_gradient!(grad, t1, t2, t3, t0, e1, e2, A, gamma_sign*nodal_data[vtx])

                # Transform slopes back to global coordinate system
                out[:, i] = @. (grad[2]*e1 + grad[3]*e2)
            end

            # Overwrite TE node values for cells on lower side
            nodal_data[TE_idx] .= nodal_data_L
            for i in cells_L
                # Convert cell vertices to a local x,y coordinate frame
                vtx = gt.get_cell(body.grid, i)
                gt.project_3d_2d!(t2, t3, e1, e2,
                                  nodes[:, vtx[1]],
                                  nodes[:, vtx[2]],
                                  nodes[:, vtx[3]])

                # The (x, y) coordinate of t1 is always at origin
                t0 = @. (t1 + t2 + t3) / 3.0

                # Get gradient of plane formed by three scalar values at vertices
                gt.get_tri_gradient!(grad, t1, t2, t3, t0, e1, e2, A, gamma_sign*nodal_data[vtx])

                # Transform slopes back to global coordinate system
                out[:, i] = @. (grad[2]*e1 + grad[3]*e2)
            end
        end
    end

    # Gamma / 2
    @. out *= -0.5

    # Save field in body
    if addfield
        add_field(body, "gamma_node", "scalar", nodal_data, "node")
        add_field(body, fieldname, "vector", eachcol(out), "cell")
    end
end

function calcfield_Ugradmu!(out::AbstractMatrix, mbody::MultiBody,
        areas::AbstractVector,
        normals::AbstractMatrix,
        controlpoints::AbstractMatrix, args...;
        fieldname="Ugradmu", addfield=true,
        optargs...
    )

    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==mbody.ncells ""*
    "Invalid `out` matrix."*
    " Expected size $((3, mbody.ncells)); got $(size(out))."
    @assert length(areas)==mbody.ncells ""*
        "Invalid `areas` vector."*
        " Expected length $(mbody.ncells); got $(length(areas))."
    @assert size(normals, 1)==3 && size(normals, 2)==mbody.ncells ""*
        "Invalid `normals` matrix."*
        " Expected size $((3, mbody.ncells)); got $(size(normals))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==mbody.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, mbody.ncells)); got $(size(controlpoints))."

    counter = 0

    for body in mbody.bodies

        offset = body.ncells
        thisout = view(out, 1:3, (1:offset) .+ counter)
        thisareas = view(areas, (1:offset) .+ counter)
        thisnormals = view(normals, 1:3, (1:offset) .+ counter)
        thiscontrolpoints = view(controlpoints, 1:3, (1:offset) .+ counter)

        calcfield_Ugradmu!(thisout, body, thisareas, thisnormals,
                                thiscontrolpoints, args...;
                                fieldname=fieldname, addfield=addfield,
                                optargs...)
        counter += offset
    end

    if addfield && !(fieldname in mbody.fields)
        push!(mbody.fields, fieldname)
    end

    return out
end

function calcfield_Ugradmu!(out::AbstractMatrix, body::AbstractBody; optargs...)
    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals)
    areas = calc_areas(body)

    calcfield_Ugradmu!(out, body, areas, normals, controlpoints; optargs...)
    return out
end

function calcfield_Ugradmu(body::AbstractBody; optargs...)
    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals)
    areas = calc_areas(body)

    out = zeros(3, body.ncells)
    calcfield_Ugradmu!(out, body, areas, normals, controlpoints; optargs...)
    return out
end

################################################################################
# PRESSURE FIELDS
################################################################################
"""
    calcfield_Cp!(out::Vector, body::AbstractBody, Us, Uref;
                            U_fieldname="U", fieldname="Cp")

Calculate the pressure coefficient
``C_p = 1 - \\left(\\frac{u}{U_\\mathrm{ref}}\\right)^2}``, where is the
velocity `Us` of each control point. The ``C_p`` is saved as a field named
`fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_Cp!(out::Arr1,
                        body::Union{NonLiftingBody, AbstractLiftingBody},
                        Us::Arr2, Uref::Number;
                        correct_kuttacondition=true,
                        fieldname="Cp", addfield=true
                        ) where {Arr1<:AbstractArray{<:Number,1},
                                 Arr2<:AbstractArray{<:Number,2}}

    # Calculate pressure coefficient
    for (i, U) in enumerate(eachcol(Us))
        out[i] += 1 - (norm(U)/Uref)^2
    end

    # Kutta-condition correction bringing the pressure on both sides of the TE
    # to be equal (average between upper and lower)
    if correct_kuttacondition && typeof(body) <: AbstractLiftingBody

        # Iterate over TE panels
        for (pi, nia, nib, pj, nja, njb) in eachcol(body.shedding)
            if pj != -1
                ave = (out[pi] + out[pi+1] + out[pj] + out[pj-1]) / 4
                out[pi] = ave
                out[pi+1] = ave
                out[pj] = ave
                out[pj-1] = ave
            else
                ave = (out[pi] + out[pi+1] ) / 2
                out[pi] = ave
                out[pi+1] = ave
            end
        end

    end

    # Save field in body
    if addfield
        add_field(body, fieldname, "scalar", out, "cell")
    end

    return out
end


function calcfield_Cp!(out::AbstractVector, mbody::MultiBody, Us::AbstractMatrix,
                        args...; addfield=true, fieldname="Cp", optargs...)


    # Error cases
    @assert length(out)==mbody.ncells ""*
        "Invalid `out` vector."*
        " Expected length $(mbody.ncells); got $(length(out))."
    @assert size(Us, 1)==3 && size(Us, 2)==mbody.ncells ""*
        "Invalid `Us` matrix."*
        " Expected size $((3, mbody.ncells)); got $(size(Us))."

    counter = 0

    for body in mbody.bodies

        offset = body.ncells
        thisout = view(out, (1:offset) .+ counter)
        thisUs = view(Us, 1:3, (1:offset) .+ counter)

        calcfield_Cp!(thisout, body, thisUs, args...;
                        fieldname=fieldname, addfield=addfield, optargs...)
        counter += offset
    end

    if addfield && !(fieldname in mbody.fields)
        push!(mbody.fields, fieldname)
    end

    return out
end

"""
    calcfield_Cp!(out::Vector, body::AbstractBody, Uref;
                            U_fieldname="U", fieldname="Cp")

Calculate the pressure coefficient
``C_p = 1 - \\left(\\frac{u}{U_\\mathrm{ref}}\\right)^2}``, where ``u`` is
the velocity field named `U_fieldname` under `body`. The ``C_p`` is saved
as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_Cp!(out, body, Uref; U_fieldname="U", optargs...)
    # Error case
    @assert check_field(body, U_fieldname) ""*
        "Field $(U_fieldname) not found;"*
       " Please run `calcfield_U(args...; fieldname=$(U_fieldname), optargs...)`"

    Us = hcat(get_field(body, U_fieldname)["field_data"]...)

    return calcfield_Cp!(out, body, Us, Uref; optargs...)
end

"""
    calcfield_Cp(args...; optargs...)

Similar to [`calcfield_Cp!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_Cp(body::AbstractBody, args...; optargs...) = calcfield_Cp!(zeros(body.ncells), body, args...; optargs...)




################################################################################
# FORCE FIELDS
################################################################################
"""
    calcfield_F!(out::Vector, body::AbstractBody,
                         areas::Vector, normals::Matrix, Us::Matrix,
                         Uinf::Number, rho::Number;
                         fieldname="F")

Calculate the force of each element
``F = - C_p \\frac{\\rho U_\\infty}{2} A \\hat{\\mathbf{n}}``, where ``C_p``is
calculated from the velocity `Us` at each control point, ``A`` is the area of
each element given in `areas`, and ``\\hat{\\mathbf{n}}`` is the normal of each
element given in `normals`. ``F`` is saved as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_F!(out::Arr0, body::Union{NonLiftingBody, AbstractLiftingBody},
                         areas::Arr1, normals::Arr2, Us::Arr3,
                         Uinf::Number, rho::Number;
                         correct_kuttacondition=true,
                         addfield=true, fieldname="F"
                         ) where {   Arr0<:AbstractArray{<:Number,2},
                                     Arr1<:AbstractArray{<:Number,1},
                                     Arr2<:AbstractArray{<:Number,2},
                                     Arr3<:AbstractArray{<:Number,2}}

    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==body.ncells ""*
        "Invalid `out` matrix."*
        " Expected size $((3, body.ncells)); got $(size(out))."
    @assert length(areas)==body.ncells ""*
        "Invalid `areas` vector."*
        " Expected length $(body.ncells); got $(length(areas))."
    @assert size(normals, 1)==3 && size(normals, 2)==body.ncells ""*
        "Invalid `normals` matrix."*
        " Expected size $((3, body.ncells)); got $(size(normals))."
    @assert size(Us, 1)==3 && size(Us, 2)==body.ncells ""*
        "Invalid `Us` matrix."*
        " Expected size $((3, body.ncells)); got $(size(Us))."

    # If F = -Cp * 0.5*ρ*u∞^2 * A * hat{n}, where Cp = 1 - (u/u∞)^2,
    # we can calculate F directly as F = 0.5*ρ*(u^2 - u∞^2) * A * hat{n}
    for (i, (U, area, normal)) in enumerate(zip(eachcol(Us), areas, eachcol(normals)))
        val = 0.5*rho*(norm(U)^2 - Uinf^2) * area
        out[1, i] += val*normal[1]
        out[2, i] += val*normal[2]
        out[3, i] += val*normal[3]
    end

    # Kutta-condition correction bringing the pressure on both sides of the TE
    # to be equal (average between upper and lower)
    # NOTE: This overwrites any previous force value instead of accumulating it
    if correct_kuttacondition && typeof(body) <: AbstractLiftingBody

        q = 0.5*rho*Uinf^2

        # Iterate over TE panels
        for (pi, nia, nib, pj, nja, njb) in eachcol(body.shedding)

            if pj != -1
                # Calculate average Cp, where Cp = 1 - (u/u∞)^2,
                aveCp = 1 - (   (norm(view(Us, :, pi))/Uinf)^2 +
                                (norm(view(Us, :, pi+1))/Uinf)^2 +
                                (norm(view(Us, :, pj))/Uinf)^2 +
                                (norm(view(Us, :, pj-1))/Uinf)^2
                            ) / 4

                # Convert Cp to force as F = -Cp * 0.5*ρ*u∞^2 * A * hat{n}
                out[1, pi] = -aveCp * q * areas[pi] * normals[1, pi]
                out[2, pi] = -aveCp * q * areas[pi] * normals[2, pi]
                out[3, pi] = -aveCp * q * areas[pi] * normals[3, pi]
                out[1, pi+1] = -aveCp * q * areas[pi+1] * normals[1, pi+1]
                out[2, pi+1] = -aveCp * q * areas[pi+1] * normals[2, pi+1]
                out[3, pi+1] = -aveCp * q * areas[pi+1] * normals[3, pi+1]
                out[1, pj] = -aveCp * q * areas[pj] * normals[1, pj]
                out[2, pj] = -aveCp * q * areas[pj] * normals[2, pj]
                out[3, pj] = -aveCp * q * areas[pj] * normals[3, pj]
                out[1, pj-1] = -aveCp * q * areas[pj-1] * normals[1, pj-1]
                out[2, pj-1] = -aveCp * q * areas[pj-1] * normals[2, pj-1]
                out[3, pj-1] = -aveCp * q * areas[pj-1] * normals[3, pj-1]

            else
                # Calculate average Cp, where Cp = 1 - (u/u∞)^2,
                aveCp = 1 - (   (norm(view(Us, :, pi))/Uinf)^2 +
                                (norm(view(Us, :, pi+1))/Uinf)^2
                            ) / 2

                # Convert Cp to force as F = -Cp * 0.5*ρ*u∞^2 * A * hat{n}
                out[1, pi] = -aveCp * q * areas[pi] * normals[1, pi]
                out[2, pi] = -aveCp * q * areas[pi] * normals[2, pi]
                out[3, pi] = -aveCp * q * areas[pi] * normals[3, pi]
                out[1, pi+1] = -aveCp * q * areas[pi+1] * normals[1, pi+1]
                out[2, pi+1] = -aveCp * q * areas[pi+1] * normals[2, pi+1]
                out[3, pi+1] = -aveCp * q * areas[pi+1] * normals[3, pi+1]

            end
        end

    end

    # Save field in body
    if addfield
        add_field(body, fieldname, "vector", eachcol(out), "cell")
    end

    return out
end


function calcfield_F!(out::AbstractMatrix, mbody::MultiBody,
                        areas::AbstractVector, normals::AbstractMatrix, Us::AbstractMatrix,
                        args...;
                        addfield=true, fieldname="F",
                        optargs...)


    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==mbody.ncells ""*
        "Invalid `out` matrix."*
        " Expected size $((3, mbody.ncells)); got $(size(out))."
    @assert length(areas)==mbody.ncells ""*
        "Invalid `areas` vector."*
        " Expected length $(mbody.ncells); got $(length(areas))."
    @assert size(normals, 1)==3 && size(normals, 2)==mbody.ncells ""*
        "Invalid `normals` matrix."*
        " Expected size $((3, mbody.ncells)); got $(size(normals))."
    @assert size(Us, 1)==3 && size(Us, 2)==mbody.ncells ""*
        "Invalid `Us` matrix."*
        " Expected size $((3, mbody.ncells)); got $(size(Us))."

    counter = 0

    for body in mbody.bodies

        offset = body.ncells
        thisout = view(out, 1:3, (1:offset) .+ counter)
        thisareas = view(areas, (1:offset) .+ counter)
        thisnormals = view(normals, 1:3, (1:offset) .+ counter)
        thisUs = view(Us, 1:3, (1:offset) .+ counter)

        calcfield_F!(thisout, body, thisareas, thisnormals,
                                thisUs, args...;
                                fieldname=fieldname, addfield=addfield,
                                optargs...)
        counter += offset
    end

    if addfield && !(fieldname in mbody.fields)
        push!(mbody.fields, fieldname)
    end

    return out
end

"""
    calcfield_F!(out::Vector, body::AbstractBody,
                            Uinf::Number, rho::Number;
                            U_fieldname="U", optargs...
                         )

Calculate the force of each element
``F = - C_p \\frac{\\rho U_\\infty}{2} A \\hat{\\mathbf{n}}``, where ``C_p``is
calculated from the velocity `Us` field `U_fieldname`, ``A`` is the area of
each element, and ``\\hat{\\mathbf{n}}`` is the normal of each element. ``F``
is saved as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_F!(out::Arr, body::AbstractBody,
                        Uinf::Number, rho::Number;
                        U_fieldname="U", optargs...
                     ) where {Arr<:AbstractArray{<:Number,2}}
    # Error cases
    @assert check_field(body, U_fieldname) ""*
        "Field $(U_fieldname) not found;"*
        " Please run `calcfield_U(args...; fieldname=$(U_fieldname), optargs...)`"

    Us = hcat(get_field(body, U_fieldname)["field_data"]...)
    areas = calc_areas(body)
    normals = calc_normals(body; flipbyCPoffset=true)

    return calcfield_F!(out, body, areas, normals, Us, Uinf, rho; optargs...)
end

"""
    calcfield_F(args...; optargs...)

Similar to [`calcfield_F!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_F(body::AbstractBody, args...; optargs...) = calcfield_F!(zeros(3, body.ncells), body, args...; optargs...)


"""
    calcfield_sectionalforce!(outf::Matrix, outpos::Vector,
                                        body::Union{NonLiftingBody, AbstractLiftingBody},
                                        controlpoints::Matrix, Fs::Matrix;
                                        dimspan=2, dimchord=1,
                                        spandirection=[0, 1, 0],
                                        fieldname="sectionalforce"
                                        )

Calculate the sectional force (a vectorial force per unit span) along the span.
This is calculated from the force `Fs` and the control points `controlpoints`
and saved as a field named `fieldname`.

The field is calculated in-place on `outf` while the spanwise position of each
section is stored under `outpos`.
"""
function calcfield_sectionalforce!(outf::Arr0, outpos::Arr1,
                                    body::Union{NonLiftingBody, AbstractLiftingBody},
                                    controlpoints::Arr2, Fs::Arr3;
                                    dimspan=2, dimchord=1,
                                    spandirection=[0, 1, 0],
                                    fieldname="sectionalforce", addfield=true
                                    ) where {   Arr0<:AbstractArray{<:Number,2},
                                                Arr1<:AbstractArray{<:Number,1},
                                                Arr2<:AbstractArray{<:Number,2},
                                                Arr3<:AbstractArray{<:Number,2}}



    lin, gdims = get_linearindex(body)      # LinearIndex and grid dimensions

    # Error cases
    @assert size(outf, 1)==3 && size(outf, 2)==gdims[dimspan] ""*
        "Invalid `outf` matrix."*
        " Expected size $((3, gdims[dimspan])); got $(size(outf))."
    @assert length(outpos)==gdims[dimspan] ""*
        "Invalid `outpos` matrix."*
        " Expected length $(gdims[dimspan]); got $(length(outpos))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==body.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, body.ncells)); got $(size(controlpoints))."
    @assert size(Fs, 1)==3 && size(Fs, 2)==body.ncells ""*
        "Invalid `Fs` matrix."*
        " Expected size $((3, body.ncells)); got $(size(Fs))."

    # Pre-allocate memory
    coor = ones(Int, 3)                     # Cartesian coordinates (indices)
    lincoors = zeros(Int, gdims[dimchord])  # Linear coordinate (index)
    outf .= 0

    # Integrate force in the chordwise direction along the span
    for j in 1:gdims[dimspan] # Iterate over span

        for i in 1:gdims[dimchord] # Iterate over chord

            coor[dimchord] = i
            coor[dimspan] = j
            lincoors[i] = lin[coor...]

            # Add force to this section
            outf[1, j] += Fs[1, lincoors[i]]
            outf[2, j] += Fs[2, lincoors[i]]
            outf[3, j] += Fs[3, lincoors[i]]

        end

        # Calculate span position of this section
        spanpos = mean(dot(spandirection, Xcp)
                        for Xcp in eachcol(view(controlpoints, :, lincoors)))
        outpos[j] = spanpos

    end

    # Convert force to be per unit span
    for j in 1:gdims[dimspan] # Iterate over span
        deltapos =  j==1 ?              outpos[j+1]-outpos[j] :
                    j==length(outpos) ? outpos[j]-outpos[j-1] :
                                        (outpos[j+1]-outpos[j-1])/2

        outf[:, j] /= abs(deltapos)
    end

    # Save field in body
    if addfield
        add_field(body, fieldname, "vector", eachcol(outf), "system")
        add_field(body, fieldname*"-pos", "vector", eachcol(outpos), "system")
    end

    return outf, outpos
end

"""
    calcfield_sectionalforce!(outFs::Matrix, outpos::Vector,
                                    body::Union{NonLiftingBody, AbstractLiftingBody};
                                    F_fieldname="F", optargs...
                                    )

Calculate the sectional force (a vectorial force per unit span) along the span.
This is calculated from the force field `F_fieldname` and saved as a field named
`fieldname`.

The field is calculated in-place on `outFs` while the spanwise position of each
section is stored under `outpos`.
"""
function calcfield_sectionalforce!(outFs::Arr0, outpos::Arr1,
                                    body::Union{NonLiftingBody, AbstractLiftingBody};
                                    F_fieldname="F",
                                    offset=nothing, characteristiclength=nothing,
                                    optargs...
                                    ) where {   Arr0<:AbstractArray{<:Number,2},
                                                Arr1<:AbstractArray{<:Number,1}}
    # Error cases
    @assert check_field(body, F_fieldname) ""*
        "Field $(F_fieldname) not found;"*
        " Please run `calcfield_F(args...; fieldname=$(F_fieldname), optargs...)`"

    Fs = hcat(get_field(body, F_fieldname)["field_data"]...)

    # Optional arguments for calc_controlpoints
    cp_optargs = (off=offset, characteristiclength=characteristiclength)
    cp_optargs = ((key, val) for (key, val) in pairs(cp_optargs) if val!=nothing)

    # Calculate control points
    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals; cp_optargs...)

    return calcfield_sectionalforce!(outFs, outpos, body,
                                            controlpoints, Fs; optargs...)
end


"""
    calcfield_sectionalforce(args...; optargs...)

Similar to [`calcfield_sectionalforce!`](@ref) but without in-place calculation
(`outFs` nor `outpos` are needed).
"""
function calcfield_sectionalforce(body::Union{NonLiftingBody, AbstractLiftingBody}, args...;
                                                        dimspan=2, optargs...)

    lin, gdims = get_linearindex(body)      # LinearIndex and grid dimensions

    outFs = zeros(3, gdims[dimspan])
    outpos = zeros(gdims[dimspan])

    return calcfield_sectionalforce!(outFs, outpos, body, args...;
                                                    dimspan=dimspan, optargs...)
end

"""
    calcfield_Ftot!(out::AbstractVector, body::AbstractBody,
                            Fs::AbstractMatrix; fieldname="Ftot")

Calculate the integrated force of this body, which is a three-dimensional vector.
This is calculated from the force of each element given in `Fs` and saved as a
field named `fieldname`.

The field is calculated in-place and added to `out`.
"""
function calcfield_Ftot!(out::AbstractVector, body::AbstractBody,
                            Fs::AbstractMatrix; fieldname="Ftot", addfield=true)

    # Error case
    @assert length(out)==3 ""*
        "Invalid `out` vector. Expected length $(3); got $(length(out))."

    for i in 1:3
        out[i] += sum(view(Fs, i, :))
    end

    # Save field in body
    if addfield
        add_field(body, fieldname, "vector", out, "system")
    end

    return out
end

"""
    calcfield_Ftot!(out::AbstractVector, body::AbstractBody;
                                    F_fieldname="F", optargs...)

Calculate the integrated force of this body, which is a three-dimensional vector.
This is calculated from the force field `F_fieldname` and saved as a field named
`fieldname`.

The field is calculated in-place and added to `out`.
"""
function calcfield_Ftot!(out, body; F_fieldname="F", optargs...)
    # Error case
    @assert check_field(body, F_fieldname) ""*
        "Field $(F_fieldname) not found;"*
        " Please run `calcfield_F(args...; fieldname=$(F_fieldname), optargs...)`"

    Fs = hcat(get_field(body, F_fieldname)["field_data"]...)

    return calcfield_Ftot!(out, body, Fs; optargs...)
end

"""
    calcfield_Ftot(body, args...; optargs...) = calcfield_Ftot!(zeros(3), body, args...; optargs...)

Similar to [`calcfield_Ftot!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_Ftot(body, args...; optargs...) = calcfield_Ftot!(zeros(3), body, args...; optargs...)

"""
    calcfield_LDS!(out::Matrix, body::AbstractBody, Fs::Matrix,
                    Lhat::Vector, Dhat::Vector, Shat::Vector)

Calculate the integrated force decomposed as lift, drag, and sideslip according
to the orthonormal basis `Lhat`, `Dhat`, `Shat`.
This is calculated from the force of each element given in `Fs`.
`out[:, 1]` is the lift vector and is saved as the field "L".
`out[:, 2]` is the drag vector and is saved as the field "D".
`out[:, 3]` is the sideslip vector and is saved as the field "S".

The field is calculated in-place on `out`.
"""
function calcfield_LDS!(out::AbstractMatrix, body::AbstractBody,
                        Fs::AbstractMatrix,
                        Lhat::AbstractVector, Dhat::AbstractVector,
                        Shat::AbstractVector;
                        addfield=true)
    # Error case
    @assert size(out, 1)==3 && size(out, 2)==3 ""*
        "Invalid `out` matrix. Expected size $((3, 3)); got $(size(out))."
    @assert abs(norm(Lhat) - 1) <= 2*eps() ""*
        "Lhat=$(Lhat) is not a unitary vector"
    @assert abs(norm(Dhat) - 1) <= 2*eps() ""*
        "Dhat=$(Dhat) is not a unitary vector"
    @assert abs(norm(Shat) - 1) <= 2*eps() ""*
        "Shat=$(Shat) is not a unitary vector"

    # Calculate Ftot (integrated force)
    for i in 1:3
        out[i, 3] += sum(view(Fs, i, :))
    end

    # Project Ftot in each direction
    out[:, 1] = Lhat
    out[:, 1] *= dot(view(out, :, 3), Lhat)
    out[:, 2] = Dhat
    out[:, 2] *= dot(view(out, :, 3), Dhat)
    aux = dot(view(out, :, 3), Shat)
    out[:, 3] = Shat
    out[:, 3] *= aux

    # Save field in body
    if addfield
        add_field(body, "L", "vector", view(out, :, 1), "system")
        add_field(body, "D", "vector", view(out, :, 2), "system")
        add_field(body, "S", "vector", view(out, :, 3), "system")
    end

    return out
end

"""
    calcfield_LDS!(out::Matrix, body::AbstractBody,
                    Lhat::Vector, Dhat::Vector, Shat::Vector; F_fieldname="F")

Calculate the integrated force decomposed as lift, drag, and sideslip according
to the orthonormal basis `Lhat`, `Dhat`, `Shat`.
This is calculated from the force field `F_fieldname`.
"""
function calcfield_LDS!(out, body, Lhat, Dhat, Shat; F_fieldname="F", optargs...)
    # Error case
    @assert check_field(body, F_fieldname) ""*
        "Field $(F_fieldname) not found;"*
        " Please run `calcfield_F(args...; fieldname=$(F_fieldname), optargs...)`"

    Fs = hcat(get_field(body, F_fieldname)["field_data"]...)

    return calcfield_LDS!(out, body, Fs, Lhat, Dhat, Shat; optargs...)
end

"""
    calcfield_LDS!(out, body, Lhat, Dhat; optargs...)

`Shat` is calculated automatically from `Lhat` and `Dhat`,
"""
function calcfield_LDS!(out, body, Lhat, Dhat; optargs...)
    return calcfield_LDS!(out, body, Lhat, Dhat, cross(Lhat, Dhat); optargs...)
end


"""
    calcfield_LDS(body, args...; optargs...) = calcfield_LDS!(zeros(3, 3), body, args...; optargs...)

Similar to [`calcfield_LDS!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_LDS(body, args...; optargs...) = calcfield_LDS!(zeros(3, 3), body, args...; optargs...)









################################################################################
# MOMENT FIELDS
################################################################################
"""
    calcfield_Mtot!(out::AbstractVector, body::AbstractBody,
                                Xac::AbstractVector, controlpoints::AbstractMatrix,
                                Fs::AbstractMatrix;
                                fieldname="Ftot", addfield=true)

Calculate the integrated moment of this body (which is a three-dimensional
vector) with respect to the aerodynamic center `Xac`.
This is calculated from the force and position of each element given in `Fs`
and `controlpoints`, respectively, and saved as a field named `fieldname`.

The field is calculated in-place and added to `out`.
"""
function calcfield_Mtot!(out::AbstractVector, body::AbstractBody,
                            Xac::AbstractVector, controlpoints::AbstractMatrix,
                            Fs::AbstractMatrix;
                            fieldname="Mtot", addfield=true)
    # Error case
    @assert length(out)==3 ""*
        "Invalid `out` vector. Expected length 3; got $(length(out))."
    @assert length(Xac)==3 ""*
        "Invalid `Xac` vector. Expected length 3; got $(length(Xac))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==body.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, body.ncells)); got $(size(controlpoints))."
    @assert size(Fs, 1)==3 && size(Fs, 2)==body.ncells ""*
        "Invalid `Fs` matrix."*
        " Expected size $((3, body.ncells)); got $(size(Fs))."

    # Calculate Mtot (integrated moment)
    for (X, F) in zip(eachcol(controlpoints), eachcol(Fs))
        out[1] += (X[2] - Xac[2])*F[3] - (X[3] - Xac[3])*F[2]
        out[2] += (X[3] - Xac[3])*F[1] - (X[1] - Xac[1])*F[3]
        out[3] += (X[1] - Xac[1])*F[2] - (X[2] - Xac[2])*F[1]
    end

    # Save field in body
    if addfield
        add_field(body, fieldname, "vector", out, "system")
    end

    return out
end

"""
    calcfield_Mtot!(out, body, Xac; F_fieldname="F",
                    offset=nothing, characteristiclength=nothing, optargs...)

Calculate the integrated moment of this body (which is a three-dimensional
vector) with respect to the aerodynamic center `Xac`.
This is calculated from the force field `F_fieldname` and saved as a field named
`fieldname`.

The field is calculated in-place and added to `out`.
"""
function calcfield_Mtot!(out, body, Xac; F_fieldname="F",
                            offset=nothing, characteristiclength=nothing,
                            optargs...)
    # Error case
    @assert check_field(body, F_fieldname) ""*
        "Field $(F_fieldname) not found;"*
        " Please run `calcfield_F(args...; fieldname=$(F_fieldname), optargs...)`"

    Fs = hcat(get_field(body, F_fieldname)["field_data"]...)

    # Optional arguments for calc_controlpoints
    cp_optargs = (off=offset, characteristiclength=characteristiclength)
    cp_optargs = ((key, val) for (key, val) in pairs(cp_optargs) if val!=nothing)

    # Calculate control points
    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals; cp_optargs...)

    return calcfield_Mtot!(out, body, Xac, controlpoints, Fs; optargs...)
end

"""
    calcfield_Mtot(body, args...; optargs...) = calcfield_Mtot!(zeros(3), body, args...; optargs...)

Similar to [`calcfield_Mtot!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_Mtot(body, args...; optargs...) = calcfield_Mtot!(zeros(3), body, args...; optargs...)



"""
    calcfield_lmn!(out::Matrix, body::AbstractBody,
                    Xac::AbstractVector, controlpoints::AbstractMatrix,
                    Fs::Matrix, lhat::Vector, mhat::Vector, nhat::Vector)

Calculate the integrated moment of this body with respect to the aerodynamic
center `Xac` and decompose it as rolling, pitching, and yawing moments according
to the orthonormal basis `lhat`, `mhat`, `nhat`, repsectively.
This is calculated from the force and position of each element given in `Fs`
and `controlpoints`, respectively.
`out[:, 1]` is the rolling moment vector and is saved as the field "Mroll".
`out[:, 2]` is the pitching moment vector and is saved as the field "Mpitch".
`out[:, 3]` is the yawing moment vector and is saved as the field "Myaw".

The field is calculated in-place on `out`.
"""
function calcfield_lmn!(out::AbstractMatrix, body::AbstractBody,
                        Xac::AbstractVector, controlpoints::AbstractMatrix,
                        Fs::AbstractMatrix,
                        lhat::AbstractVector, mhat::AbstractVector,
                        nhat::AbstractVector;
                        addfield=true)
    # Error case
    @assert size(out, 1)==3 && size(out, 2)==3 ""*
        "Invalid `out` matrix. Expected size $((3, 3)); got $(size(out))."
    @assert length(Xac)==3 ""*
        "Invalid `Xac` vector. Expected length 3; got $(length(Xac))."
    @assert size(controlpoints, 1)==3 && size(controlpoints, 2)==body.ncells ""*
        "Invalid `controlpoints` matrix."*
        " Expected size $((3, body.ncells)); got $(size(controlpoints))."
    @assert size(Fs, 1)==3 && size(Fs, 2)==body.ncells ""*
        "Invalid `Fs` matrix."*
        " Expected size $((3, body.ncells)); got $(size(Fs))."
    @assert abs(norm(lhat) - 1) <= 2*eps() ""*
        "lhat=$(lhat) is not a unitary vector"
    @assert abs(norm(mhat) - 1) <= 2*eps() ""*
        "mhat=$(mhat) is not a unitary vector"
    @assert abs(norm(nhat) - 1) <= 2*eps() ""*
        "nhat=$(nhat) is not a unitary vector"

    # Calculate Mtot (integrated moment)
    for (X, F) in zip(eachcol(controlpoints), eachcol(Fs))
        out[1, 3] += (X[2] - Xac[2])*F[3] - (X[3] - Xac[3])*F[2]
        out[2, 3] += (X[3] - Xac[3])*F[1] - (X[1] - Xac[1])*F[3]
        out[3, 3] += (X[1] - Xac[1])*F[2] - (X[2] - Xac[2])*F[1]
    end

    # Project Mtot in each direction
    out[:, 1] = lhat
    out[:, 1] *= dot(view(out, :, 3), lhat)
    out[:, 2] = mhat
    out[:, 2] *= dot(view(out, :, 3), mhat)
    aux = dot(view(out, :, 3), nhat)
    out[:, 3] = nhat
    out[:, 3] *= aux

    # Save field in body
    if addfield
        add_field(body, "Mroll", "vector", view(out, :, 1), "system")
        add_field(body, "Mpitch", "vector", view(out, :, 2), "system")
        add_field(body, "Myaw", "vector", view(out, :, 3), "system")
    end

    return out
end


"""
    calcfield_lmn!(out, body, Xac, lhat, mhat, nhat; F_fieldname="F",
                    offset=nothing, characteristiclength=nothing, optargs...)

Calculate the integrated moment of this body with respect to the aerodynamic
center `Xac` and decompose it as rolling, pitching, and yawing moments according
to the orthonormal basis `lhat`, `mhat`, `nhat`, repsectively.
This is calculated from the force field `F_fieldname`.

The field is calculated in-place on `out`.
"""
function calcfield_lmn!(out, body, Xac, lhat, mhat, nhat; F_fieldname="F",
                            offset=nothing, characteristiclength=nothing,
                            optargs...)
    # Error case
    @assert check_field(body, F_fieldname) ""*
        "Field $(F_fieldname) not found;"*
        " Please run `calcfield_F(args...; fieldname=$(F_fieldname), optargs...)`"

    Fs = hcat(get_field(body, F_fieldname)["field_data"]...)

    # Optional arguments for calc_controlpoints
    cp_optargs = (off=offset, characteristiclength=characteristiclength)
    cp_optargs = ((key, val) for (key, val) in pairs(cp_optargs) if val!=nothing)

    # Calculate control points
    normals = calc_normals(body)
    controlpoints = calc_controlpoints(body, normals; cp_optargs...)

    return calcfield_lmn!(out, body, Xac, controlpoints, Fs, lhat, mhat, nhat;
                                                                     optargs...)
end

"""
    calcfield_lmn!(out, body, Xac, lhat, mhat; optargs...)

`nhat` is calculated automatically from `lhat` and `mhat`,
"""
function calcfield_lmn!(out, body, Xac, lhat, mhat; optargs...)
    return calcfield_lmn!(out, body, Xac, lhat, mhat, cross(lhat, mhat); optargs...)
end


"""
    calcfield_lmn(body, args...; optargs...) = calcfield_lmn!(zeros(3, 3), body, args...; optargs...)

Similar to [`calcfield_lmn!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_lmn(body, args...; optargs...) = calcfield_lmn!(zeros(3, 3), body, args...; optargs...)
