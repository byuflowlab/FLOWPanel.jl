#=##############################################################################
# DESCRIPTION
    Definition of methods for post-processing solver results.

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Oct 2022
  * License     : MIT License
=###############################################################################

"""
    calcfield_U!(out::Matrix,
                 sourcebody::AbstractBody, targetbody::AbstractBody,
                 controlpoints::Matrix; fieldname="U")

Calculate the velocity induced by `sourcebody` on `controlpoints` and save
it as a field of name `fieldname` under `targetbody`. The field includes the
freestream velocity stored as field `\"Uinf\"` in `targetbody`.

The field is calculated in place on `out`.
"""
function calcfield_U!(out::Arr1, sourcebody::AbstractBody, targetbody::AbstractBody,
                        controlpoints::Arr2; fieldname="U"
                        ) where {   Arr1<:AbstractArray{<:Number,2},
                                    Arr2<:AbstractArray{<:Number,2}}

    # ERROR CASES
    if check_solved(sourcebody)==false
        error("Source body hasn't been solved yet."*
              " Please call `solve(...)` function first.")
    elseif check_field(targetbody, "Uinf")==false
        error("Target body doesn't have freestream field `\"Uinf\"`.")*
              " Please call `add_field(targetbody, \"Uinf\", ...)` first."
    elseif size(controlpoints, 1)!=3 || size(controlpoints, 2)!=targetbody.ncells
        error("Invalid `controlpoints` matrix."*
              " Expected size $((3, targetbody.ncells)); got $(size(controlpoints)).")
    elseif size(out, 1)!=3 || size(out, 2)!=targetbody.ncells
        error("Invalid `out` matrix."*
              " Expected size $((3, targetbody.ncells)); got $(size(out)).")
    end

    # Add freestream
    for (i, Uinf) in enumerate(get_field(targetbody, "Uinf")["field_data"])
        out[:, i] .= Uinf
    end

    # Add induced velocity at each control point
    Uind!(sourcebody, controlpoints, out)

    # Save field in body
    add_field(targetbody, fieldname, "vector", eachcol(out), "cell")

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

The field is calculated in place on `out`.
"""
function calcfield_U!(out::Arr, sourcebody::AbstractBody, targetbody::AbstractBody;
                        offset=nothing, characteristiclength=nothing, optargs...
                        ) where {Arr<:AbstractArray{<:Number,2}}

    # Optional arguments for _calc_controlpoints
    optargs = (off=offset, characteristiclength=characteristiclength)
    optargs = ((key, val) for (key, val) in optargs if val!=nothing)

    # Calculate control points
    normals = _calc_normals(targetbody)
    controlpoints = _calc_controlpoints(targetbody, normals; optargs...)

    # Calculate field on control points
    calcfield_U!(sourcebody, targetbody, controlpoints, out; optargs...)
end

"""
    calcfield_U(args...; optargs...)

Just like `calcfield_U` but without in-place calculation (`out` is not needed).
"""
function calcfield_U(sourcebody, targetbody, controlpoints; optargs...)
    out = similar(controlpoints)
    return calcfield_U!(out, sourcebody, targetbody, controlpoints; optargs...)
end
function calcfield_U(sourcebody, targetbody::AbstractBody; optargs...)
    out = zeros(3, targetbody.ncells)
    return calcfield_U!(out, sourcebody, targetbody; optargs...)
end

"""
    calcfield_Uoff!(args...; optargs...) = calcfield_U(args...; optargs..., fieldname="Uoff")

See documentation of `calcfield_U!(...)`.
"""
calcfield_Uoff!(args...; optargs...) = calcfield_U!(args...; optargs..., fieldname="Uoff")
calcfield_Uoff(args...; optargs...) = calcfield_U(args...; optargs..., fieldname="Uoff")

"""
    calcfield_Cp!(out::Vector, body::AbstractBody, Uref;
                            U_fieldname="U", fieldname="Cp")

Calculate the pressure coefficient
``C_p = 1 - \\left(\\frac{u}{U_\\mathrm{ref}}\\right)^2}``, where ``u`` is
the velocity field named `U_fieldname` under `body`. The ``C_p`` is saved
as a field named `fieldname`.

The field is calculated in place on `out`.
"""
function calcfield_Cp!(out::Vector, body::AbstractBody, Uref::Number;
                                                U_fieldname="U", fieldname="Cp")

    # Error cases
    @assert check_field(body, U_fieldname) "Field $(U_fieldname) not found;"*
                                           "Please run `calcfield_U(args...; fieldname=$(U_fieldname), optargs...)`"

    # Calculate pressure coefficient
    U = get_field(body, U_fieldname)["field_data"]
    map!(u -> 1 - (norm(u)/Uref)^2, out, U)

    # Save field in body
    add_field(body, fieldname, "scalar", out, "cell")

    return out
end

"""
    calcfield_Cp(args...; optargs...)

Just like `calcfield_Cp` but without in-place calculation (`out` is not needed).
"""
calcfield_Cp(body::AbstractBody, args...; optargs...) = calcfield_Cp!(zeros(body.ncells), body, args...; optargs...)

"""
    slicefield(body::AbstractBody, fieldname::String,
                position::Number, direction::Vector, row::Bool)

Return a slice of the field `fieldname` in `body` corresponding to the row or
column ("row" is the first dimension of the grid, "column" is the second
dimension) that is the closest to `position` calculated as the projection
of the average cell position in the direction `direction`.

> **Example:** For a wing with its span aligned along the y-axis, the pressure
along a slice of the wing at the spanwise position y=0.5 is obtained as
`slicefield(wing, "Cp", 0.5, [0, 1, 0], false)`.
"""
slicefield(body::AbstractBody, args...; optargs...) = slicefield(body, _calc_controlpoints(body), args...; optargs...)

"""
    slicefield(body::AbstractBody, controlpoints::Matrix,
                    fieldname::String,
                    position::Number, direction::Vector, row::Bool)

Same thing, but with the option of providing the control points as to save
wastefull memory allocation.
"""
function slicefield(body::AbstractBody, controlpoints::Arr,
                    fieldname::String,
                    position::Number, direction::Vector, row::Bool;
                    reduce=true
                    ) where {Arr<:AbstractArray{<:Number,2}}

    # Fetch field
    field = get_field(body, fieldname)["field_data"]

    # Find index of row or column slicing the field
    gdim = row ? 1 : 2                          # Dimension to slice
    islice, pos, errmin, lin = find_i(body, controlpoints,
                                            position, gdim, -1; xdir=direction)
    # Slice field
    ncell = get_ndivscells(body)[row ? 2 : 1]   # Number of cells in the slice
    indices = collect(row==1 ? lin[islice, j] : lin[j, islice] for j in 1:ncell)

    slice = field[indices]

    # Points along the slice
    slicepoints = controlpoints[:, indices]

    # Reduce the implicit double-column of the triangular grid into one column
    if reduce
        slice = [ (slice[2*(i-1) + 1]+slice[2*(i-1) + 2])/2 for i in 1:Int(length(slice)/2)]
        slicepoints = [ (slicepoints[i, 2*(j-1) + 1]+slicepoints[i, 2*(j-1) + 2])/2 for i in 1:size(slicepoints,1), j in 1:Int(size(slicepoints, 2)/2)]
    end

    return slicepoints, slice
end
