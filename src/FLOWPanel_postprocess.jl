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
    calcfield_U!(bodies; backend=DirectBackend(), reset=true, convolve_panels=true, doublet_gradient=true)
    calcfield_U!(body; optargs...)

Compute and store the control-point velocity field for one or more bodies.
"""
function calcfield_U!(bodies::Tuple, uinf=SVector{3,Float64}(0.0, 0.0, 0.0), wakes=Tuple(nothing for _ in bodies);
        backend::AbstractBackend=DirectBackend(),
        reset=true,
        convolve_panels=true,
        doublet_gradient=true
    )

    # reset velocity
    if reset
        for body in bodies
            body.velocity .= zero(eltype(body.velocity))
        end
    end

    # recalculate normals/control points on the target body
    for body in bodies
        calc_normals!(body)
        calc_controlpoints!(body; off=abs(body.CPoffset))
    end

    # apply freestream at control points
    for body in bodies
        apply_freestream!(body, uinf)
    end

    # add wake-induced velocity at control points
    wakes = Tuple(w for w in wakes if !isnothing(w))
    length(wakes) > 0 && influence!(bodies, wakes, backend; scalar_potential=false, velocity=true)

    # Add induced velocity at each control point
    convolve_panels && influence!(bodies, bodies, backend; scalar_potential=false, velocity=true)

    # add doublet gradient (if applicable)
    for body in bodies
        if has_grad_mu(body) && doublet_gradient
            compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
                body.cells,
                body.neighbor,
                view(body.strength, :, get_Gammai(body)),
                view(body.shedding_full, 1:2, :), scale=0.5)
            
            # alternatively, comment out the above function and:
            # targetbody.velocity .*= 2.0

            # alternatively, comment out the above and:
            # targetbody.velocity .= zero(eltype(targetbody.velocity))
            # compute_mu_gradient!(targetbody.velocity, targetbody.controlpoints, targetbody.normals,
                # targetbody.cells,
                # targetbody.neighbor,
                # view(targetbody.strength, :, get_Gammai(targetbody)),
                # view(targetbody.shedding_full, 1:2, :), 
                # scale=1.0)
        end
    end

    return nothing
end

calcfield_U!(targetbody; optargs...) = calcfield_U!((targetbody,); optargs...)

################################################################################
# GRADIENT COMPUTATION
################################################################################

"""
    compute_mu_gradient!(grad_mu, controlpoints, normals, cells, neighbors, mu, te_info; scale=0.5)

Compute the surface gradient of doublet or vortex-ring strength using a local
least-squares stencil, with one-sided handling at trailing-edge panels.
"""
function compute_mu_gradient!(grad_mu,
                            controlpoints::AbstractMatrix{Float64},
                            normals::AbstractMatrix{Float64},
                            cells::AbstractMatrix{Int},
                            neighbors::AbstractMatrix{Int},
                            mu::AbstractVector{Float64},
                            te_info::AbstractMatrix{Int};
                            scale=0.5
                        )

    _, M = size(cells)

    # Pre-allocate array for stencil to avoid allocations inside the loop
    stencil = Int[]
    sizehint!(stencil, 10)

    for i in 1:M
        empty!(stencil)

        # Check if current panel is at the Trailing Edge
        is_te = te_info[1, i] > 0 && te_info[2, i] > 0

        if is_te
            # Global vertex IDs of the TE edge
            te_v1 = cells[te_info[1, i], i]
            te_v2 = cells[te_info[2, i], i]
        else
            te_v1 = -1
            te_v2 = -1
        end

        # 1. Gather immediate valid neighbors
        for k in 1:3
            n_idx = neighbors[k, i]
            if n_idx <= 0
                continue
            end

            if is_te
                # Check if this neighbor shares BOTH trailing edge vertices
                has_v1 = (cells[1, n_idx] == te_v1) || (cells[2, n_idx] == te_v1) || (cells[3, n_idx] == te_v1)
                has_v2 = (cells[1, n_idx] == te_v2) || (cells[2, n_idx] == te_v2) || (cells[3, n_idx] == te_v2)

                # If it shares the TE edge, exclude it from the stencil (one-sided bias)
                if has_v1 && has_v2
                    continue
                end
            end
            push!(stencil, n_idx)
        end

        # 2. Expand stencil upstream if it's a TE panel to ensure robust Least Squares
        if is_te
            n_current = length(stencil)
            for s in 1:n_current
                s_idx = stencil[s]
                
                # Look at neighbors of the interior neighbors
                for k in 1:3
                    nn_idx = neighbors[k, s_idx]
                    if nn_idx > 0 && nn_idx != i && !(nn_idx in stencil)
                        push!(stencil, nn_idx)
                    end
                end
            end
        end

        # 3. Build normal equations for Least Squares: (A^T A) * grad = A^T b
        ATA = zeros(Float64, 3, 3)
        ATb = zeros(Float64, 3)
        mean_sq_dist = 0.0

        for j in stencil
            dx = controlpoints[1, j] - controlpoints[1, i]
            dy = controlpoints[2, j] - controlpoints[2, i]
            dz = controlpoints[3, j] - controlpoints[3, i]
            dmu = mu[j] - mu[i]

            ATA[1, 1] += dx * dx; ATA[1, 2] += dx * dy; ATA[1, 3] += dx * dz
            ATA[2, 1] += dy * dx; ATA[2, 2] += dy * dy; ATA[2, 3] += dy * dz
            ATA[3, 1] += dz * dx; ATA[3, 2] += dz * dy; ATA[3, 3] += dz * dz

            ATb[1] += dx * dmu
            ATb[2] += dy * dmu
            ATb[3] += dz * dmu

            mean_sq_dist += dx^2 + dy^2 + dz^2
        end

        mean_sq_dist = length(stencil) > 0 ? (mean_sq_dist / length(stencil)) : 1.0

        # 4. Constrain gradient to the panel surface (Penalty method)
        # This acts as an orthogonal constraint mapping it cleanly to the 3D plane
        nx = normals[1, i]
        ny = normals[2, i]
        nz = normals[3, i]
        
        # Scale penalty by geometry size (mean squared dist) for matrix conditioning
        penalty = 1e4 * mean_sq_dist
        ATA[1, 1] += penalty * nx * nx; ATA[1, 2] += penalty * nx * ny; ATA[1, 3] += penalty * nx * nz
        ATA[2, 1] += penalty * ny * nx; ATA[2, 2] += penalty * ny * ny; ATA[2, 3] += penalty * ny * nz
        ATA[3, 1] += penalty * nz * nx; ATA[3, 2] += penalty * nz * ny; ATA[3, 3] += penalty * nz * nz

        # 5. Add minor Tikhonov regularization in case of an exactly coplanar/rank-deficient system
        reg = 1e-10 * mean_sq_dist
        ATA[1, 1] += reg
        ATA[2, 2] += reg
        ATA[3, 3] += reg

        # 6. Solve the 3x3 local system
        g = ATA \ ATb
        grad_mu[1, i] -= g[1] * scale
        grad_mu[2, i] -= g[2] * scale
        grad_mu[3, i] -= g[3] * scale
    end

    return grad_mu
end

################################################################################
# PRESSURE FIELDS
################################################################################
"""
    calcfield_Cp!(out, body, Us, Uref; dphidt=nothing, correct_kuttacondition=true, clip=nothing)
    calcfield_Cp!(body, Uref; dphidt=nothing, optargs...)
    calcfield_Cp!(bodies, Uref; dphidt=..., correct_kuttacondition=..., optargs...)

Compute and store the surface pressure coefficient field.
"""
function calcfield_Cp!(out::Arr1,
                        body::Union{NonLiftingBody, AbstractLiftingBody},
                        Us::Arr2, Uref::Number;
                        dphidt::Union{Nothing, AbstractVector}=nothing,
                        correct_kuttacondition=true,
                        clip::Union{Nothing, Function}=nothing,
                        ) where {Arr1<:AbstractArray{<:Number,1},
                                 Arr2<:AbstractArray{<:Number,2}}

    # Calculate pressure coefficient
    for (i, U) in enumerate(eachcol(Us))
        out[i] += 1 - (norm(U)/Uref)^2
    end

    # Unsteady Bernoulli term: -2(∂φ/∂t) / V∞²
    if !isnothing(dphidt)
        inv_Uref2 = 2.0 / Uref^2
        for i in eachindex(out)
            out[i] -= inv_Uref2 * dphidt[i]
        end
    end

    # Kutta-condition correction bringing the pressure on both sides of the TE
    # to be equal (average between upper and lower)
    if correct_kuttacondition && typeof(body) <: AbstractLiftingBody

        # Iterate over TE panels
        for shedding in body.shedding
            for (pi, nia, nib, pj, nja, njb) in eachcol(shedding)
                if pj != -1
                    ave = (out[pi] + out[pj]) / 2
                    out[pi] = ave
                    out[pj] = ave
                end
            end
        end

    end

    # Clip values if requested
    if !isnothing(clip)
        for (i, Cp) in enumerate(out)
            out[i] = clip(Cp)
        end
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
calcfield_Cp!(body::AbstractBody, Uref; dphidt=nothing, optargs...) = calcfield_Cp!(body.Cp, body, body.velocity, Uref; dphidt, optargs...)

function calcfield_Cp!(bodies::Tuple, Uref; dphidt=fill(nothing, length(bodies)), correct_kuttacondition=fill(true, length(bodies)), optargs...) 
    for (i, body) in enumerate(bodies)
        calcfield_Cp!(body.Cp, body, body.velocity, Uref; dphidt=dphidt[i], correct_kuttacondition=correct_kuttacondition[i], optargs...)
    end
end

################################################################################
# FORCE FIELDS
################################################################################
"""
    calcfield_F!(out, body, areas, normals, Cps, Uinf, rho; correct_kuttacondition=true)
    calcfield_F!(body, Uinf, rho; correct_kuttacondition=true)
    calcfield_F!(bodies, Uinf, rho; correct_kuttacondition=...)

Compute and store distributed surface forces from the current pressure field.
"""
function calcfield_F!(out::Arr0, body::AbstractBody,
                         areas::Arr1, normals::Arr2, Cps::Arr3,
                         Uinf::Number, rho::Number;
                         correct_kuttacondition=true,
                         ) where {   Arr0<:AbstractArray{<:Number,2},
                                     Arr1<:AbstractArray{<:Number,1},
                                     Arr2<:AbstractArray{<:Number,2},
                                     Arr3<:AbstractArray{<:Number,1}}

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
    @assert length(Cps)==body.ncells ""*
        "Invalid `Cps` vector."*
        " Expected length $(body.ncells); got $(length(Cps))."

    # # If F = -Cp * 0.5*ρ*u∞^2 * A * hat{n}, where Cp = 1 - (u/u∞)^2,
    # # we can calculate F directly as F = 0.5*ρ*(u^2 - u∞^2) * A * hat{n}
    # for (i, (U, area, normal)) in enumerate(zip(eachcol(Us), areas, eachcol(normals)))
    #     val = 0.5*rho*(norm(U)^2 - Uinf^2) * area
    #     out[1, i] += val*normal[1]
    #     out[2, i] += val*normal[2]
    #     out[3, i] += val*normal[3]
    # end

    for (i, (Cp, area, normal)) in enumerate(zip(Cps, areas, eachcol(normals)))
        val = -0.5*rho*Uinf^2 * Cp * area
        out[1, i] += val*normal[1]
        out[2, i] += val*normal[2]
        out[3, i] += val*normal[3]
    end

    # Kutta-condition correction bringing the pressure on both sides of the TE
    # to be equal (average between upper and lower)
    # NOTE: This overwrites any previous force value instead of accumulating it
    if correct_kuttacondition && typeof(body) <: AbstractLiftingBody

        # if typeof(body.grid) <: gt.GridTriangleSurface{gt.Meshes.SimpleMesh}
        #     @warn "Kutta correction requested in calcfield_F, but"*
        #             " current implementation is wrong for unstructured meshes!"
        # end

        q = 0.5*rho*Uinf^2

        # Iterate over TE panels
        for shedding in body.shedding
            for (pi, nia, nib, pj, nja, njb) in eachcol(shedding)

                if pj != -1
                    # Average Cp across upper (pi) and lower (pj) TE panels
                    aveCp = (Cps[pi] + Cps[pj]) / 2

                    # Convert Cp to force as F = -Cp * 0.5*ρ*u∞^2 * A * hat{n}
                    out[1, pi] = -aveCp * q * areas[pi] * normals[1, pi]
                    out[2, pi] = -aveCp * q * areas[pi] * normals[2, pi]
                    out[3, pi] = -aveCp * q * areas[pi] * normals[3, pi]
                    out[1, pj] = -aveCp * q * areas[pj] * normals[1, pj]
                    out[2, pj] = -aveCp * q * areas[pj] * normals[2, pj]
                    out[3, pj] = -aveCp * q * areas[pj] * normals[3, pj]

                end
            end
        end
    end

    return out
end

calcfield_F!(body::AbstractBody, Uinf::Number, rho::Number; correct_kuttacondition=true) =
    calcfield_F!(body.F, body, calc_areas(body), body.normals, body.Cp, Uinf, rho; correct_kuttacondition)

function calcfield_F!(bodies::Tuple, Uinf::Number, rho::Number; correct_kuttacondition=fill(true, length(bodies)))
    for (i, body) in enumerate(bodies)
        calcfield_F!(body.F, body, calc_areas(body), body.normals, body.Cp, Uinf, rho; correct_kuttacondition=correct_kuttacondition[i])
    end
end

"""
    calcfield_sectionalforce!(outf, outpos, body, controlpoints, Fs; dimspan=2, dimchord=1, spandirection=[0, 1, 0], ...)
    calcfield_sectionalforce!(outf, outpos, body; optargs...)
    calcfield_sectionalforce(body; optargs...)

Integrate distributed panel forces into sectional loads along a chosen spanwise
direction.
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
    # if addfield
    #     add_field(body, fieldname, "vector", eachcol(outf), "system")
    #     add_field(body, fieldname*"-pos", "vector", eachcol(outpos), "system")
    # end

    return outf, outpos
end

"""
    calcfield_sectionalforce!(outFs::Matrix, outpos::Vector,
                                    body::Union{NonLiftingBody, AbstractLiftingBody};
                                    optargs...
                                    )

Calculate the sectional force (a vectorial force per unit span) along the span.
This is calculated from the force field stored in `body.F` and saved as a field
named `fieldname`.

The field is calculated in-place on `outFs` while the spanwise position of each
section is stored under `outpos`.
"""
function calcfield_sectionalforce!(outFs::Arr0, outpos::Arr1,
                                    body::Union{NonLiftingBody, AbstractLiftingBody};
                                    offset=nothing, characteristiclength=nothing,
                                    optargs...
                                    ) where {   Arr0<:AbstractArray{<:Number,2},
                                                Arr1<:AbstractArray{<:Number,1}}
    Fs = body.F

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
    calcfield_Ftot!(out, body, Fs; fieldname="Ftot")
    calcfield_Ftot!(out, body; optargs...)
    calcfield_Ftot(body, args...; optargs...)

Integrate the stored distributed force field into a total force vector.
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
    # if addfield
    #     add_field(body, fieldname, "vector", out, "system")
    # end

    return out
end

"""
    calcfield_Ftot!(out::AbstractVector, body::AbstractBody;
                                    optargs...)

Calculate the integrated force of this body, which is a three-dimensional vector.
This is calculated from the force field stored in `body.F` and saved as a field
named `fieldname`.

The field is calculated in-place and added to `out`.
"""
function calcfield_Ftot!(out, body; optargs...)
    return calcfield_Ftot!(out, body, body.F; optargs...)
end

"""
    calcfield_Ftot(body, args...; optargs...) = calcfield_Ftot!(zeros(3), body, args...; optargs...)

Similar to [`calcfield_Ftot!`](@ref) but without in-place calculation (`out` is
not needed).
"""
calcfield_Ftot(body, args...; optargs...) = calcfield_Ftot!(zeros(3), body, args...; optargs...)

"""
    calcfield_LDS!(out, body, Fs, Lhat, Dhat, Shat; addfield=true)
    calcfield_LDS!(out, body, Lhat, Dhat, Shat; optargs...)
    calcfield_LDS!(out, body, Lhat, Dhat; optargs...)
    calcfield_LDS(body, args...; optargs...)

Resolve the total force into lift, drag, and side-force components.
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
    # if addfield
    #     add_field(body, "L", "vector", view(out, :, 1), "system")
    #     add_field(body, "D", "vector", view(out, :, 2), "system")
    #     add_field(body, "S", "vector", view(out, :, 3), "system")
    # end

    return out
end

function calcfield_LDS!(out::AbstractMatrix,
                        bodies::Tuple,
                        Fs::Vector{<:AbstractMatrix},
                        Lhat::AbstractVector, Dhat::AbstractVector,
                        Shat::AbstractVector;
                        addfield=true)

    fill!(out[:,3], 0.0)

    for Fbody in Fs
        out[:,3] .+= sum(Fbody, dims=2)[:]
    end

    out[:,1] = Lhat .* dot(out[:,3], Lhat)
    out[:,2] = Dhat .* dot(out[:,3], Dhat)
    out[:,3] = Shat .* dot(out[:,3], Shat)

    return out
end

"""
    calcfield_LDS!(out::Matrix, body::AbstractBody,
                    Lhat::Vector, Dhat::Vector, Shat::Vector; optargs...)

Calculate the integrated force decomposed as lift, drag, and sideslip according
to the orthonormal basis `Lhat`, `Dhat`, `Shat`.
This is calculated from the force field stored in `body.F`.
"""
function calcfield_LDS!(out, body, Lhat, Dhat, Shat; optargs...)
    return calcfield_LDS!(out, body, body.F, Lhat, Dhat, Shat; optargs...)
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
    calcfield_Mtot!(out, body, Xac, controlpoints, Fs; fieldname="Mtot", addfield=true)
    calcfield_Mtot!(out, body, Xac; optargs...)
    calcfield_Mtot(body, args...; optargs...)

Integrate the stored distributed force field into a total moment about a
reference point.
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
    # if addfield
    #     add_field(body, fieldname, "vector", out, "system")
    # end

    return out
end

"""
    calcfield_Mtot!(out, body, Xac;
                    offset=nothing, characteristiclength=nothing, optargs...)

Calculate the integrated moment of this body (which is a three-dimensional
vector) with respect to the aerodynamic center `Xac`.
This is calculated from the force field stored in `body.F` and saved as a field
named `fieldname`.

The field is calculated in-place and added to `out`.
"""
function calcfield_Mtot!(out, body, Xac;
                            offset=nothing, characteristiclength=nothing,
                            optargs...)
    Fs = body.F

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
    calcfield_lmn!(out, body, Xac, controlpoints, Fs, lhat, mhat, nhat; addfield=true)
    calcfield_lmn!(out, body, Xac, lhat, mhat, nhat; optargs...)
    calcfield_lmn!(out, body, Xac, lhat, mhat; optargs...)
    calcfield_lmn(body, args...; optargs...)

Resolve the total moment into user-supplied rolling, pitching, and yaw axes.
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
    # if addfield
    #     add_field(body, "Mroll", "vector", view(out, :, 1), "system")
    #     add_field(body, "Mpitch", "vector", view(out, :, 2), "system")
    #     add_field(body, "Myaw", "vector", view(out, :, 3), "system")
    # end

    return out
end


"""
    calcfield_lmn!(out, body, Xac, lhat, mhat, nhat;
                    offset=nothing, characteristiclength=nothing, optargs...)

Calculate the integrated moment of this body with respect to the aerodynamic
center `Xac` and decompose it as rolling, pitching, and yawing moments according
to the orthonormal basis `lhat`, `mhat`, `nhat`, repsectively.
This is calculated from the force field stored in `body.F`.

The field is calculated in-place on `out`.
"""
function calcfield_lmn!(out, body, Xac, lhat, mhat, nhat;
                            offset=nothing, characteristiclength=nothing,
                            optargs...)
    Fs = body.F

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
