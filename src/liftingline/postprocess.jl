#=##############################################################################
# DESCRIPTION
    Methods for postprocessing solver results of non-linear lifting line method

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################

"""
Calculate forces and moments of this lifting line.

Before calling this function, the local velocity needs to be computed as 
follows:
```julia
ll.Us .= Uinfs
selfUind!(ll)
```
The local velocity is also automatically computed in `solve(ll, ...)`, so no
need to recompute it if `calc_forcesmoments` is called right after `solve`.
"""
function calc_forcesmoments(ll::LiftingLine, 
                            Uinfs::AbstractMatrix, 
                            Uinf_ref::AbstractVector, 
                            rho::Number; 

                            X0 = zeros(3),              # (m) center about which to calculate moments

                            use_Uind_for_force = true,  # Whether to use Uind as opposed to selfUind for force postprocessing
                                                        # (`true` for more accurate spanwise cd distribution, but worse integrated CD)

                            distributions = nothing,    # Give it an array and it will push spanwise distributions there

                            # Unit vectors
                            Dhat = Uinf_ref/norm(Uinf_ref), # Drag direction
                            Shat = [0, 1, 0],           # Span direction
                            Lhat = cross(Dhat, Shat),   # Lift direction

                            lhat = Dhat,                # Rolling direction
                            mhat = Shat,                # Pitching direction
                            nhat = Lhat,                # Yawing direction
                            )


    # NOTE: Coefficients must be evaluated using the velocity from 
    #       the effective horseshoes as shown below, which is automatically
    #       computed by the solver already, so these lines are commented out to
    #       avoid redundant computation
    # ll.Us .= Uinfs
    # selfUind!(ll)

    # Calculate stripwise coefficients
    calcfield_cl(ll)
    calcfield_cd(ll)
    calcfield_cm(ll)

    # Convert velocity to effective swept velocity
    # NOTE: Forces are most accurate with the velocity from the original horseshoes,
    #       as done in the conditional statement here
    if use_Uind_for_force
        ll.Us .= Uinfs
        Uind!(ll, ll.midpoints, ll.Us)
    end
    calc_UΛs!(ll, ll.Us)

    # Force per stripwise element integrating lift and drag coefficient
    calcfield_F(ll, rho)

    # Integrated force
    Ftot = calcfield_Ftot(ll)

    # Integrated force decomposed into lift and drag
    LDS = calcfield_LDS(ll, Lhat, Dhat, Shat)

    L = LDS[:, 1]
    D = LDS[:, 2]
    S = LDS[:, 3]

    # Loading distribution (force per unit span)
    if !isnothing(distributions)
        fs = calcfield_f(ll)

        lds = decompose(fs, Lhat, Dhat)

        ypos = (ll.ypositions[2:end] .+ ll.ypositions[1:end-1]) / 2
        l = lds[1, :]
        d = lds[2, :]
        s = lds[3, :]

        push!(distributions, (; spanposition=ypos, lift_distribution=l, 
                                        drag_distribution=d, side_distribution=s))
    end

    # Integrated moment
    Mtot = calcfield_Mtot(ll, X0, rho)

    # Moment decomposed into axes
    lmn = calcfield_lmn(ll, lhat, mhat, nhat)
    roll, pitch, yaw = collect(eachcol(lmn))

    # Outputs
    return (;   lift=L, drag=D, side=S, roll, pitch, yaw, 
                Ftot, Mtot,
                Dhat, Shat, Lhat, lhat, mhat, nhat)
end

"""
Same than `calc_forcesmoments` but converting forces and moments to 
non-dimensional coefficients.
"""
function calc_forcemoment_coefficients(ll::LiftingLine, 
                                            Uinfs::AbstractMatrix,              # Local inflow velocity
                                            Uinf_ref::AbstractVector,           # (m/s) reference freestream velocity vector
                                            rho::Number,                        # (kg/m^3) air density
                                            cref::Number,                       # (m) reference chord
                                            bref::Number;                       # (m) reference span
                                            Aref = cref*bref,                   # (m^2) reference area
                                            magUinf = norm(Uinf_ref),           # (m/s) reference freestream velocity magnitude
                                            distributions = nothing,            # Give it an array and it will push spanwise distributions there
                                            forcesmoments = nothing,            # Pre-calculated dimensional forces and moments
                                            optargs...)

    q = 0.5*rho*magUinf^2               # Dynamic pressure

    # Calculate dimensional forces and moment if not already provided
    if isnothing(forcesmoments)
        forcesmoments = calc_forcesmoments(ll, Uinfs, Uinf_ref, rho; 
                                                    distributions, optargs...)
    end

    # Fetch forces and moments
    (; lift, drag, side) = forcesmoments
    (; roll, pitch, yaw) = forcesmoments

    # Fetch unit vectors
    (; Dhat, Shat, Lhat) = forcesmoments
    (; lhat, mhat, nhat) = forcesmoments

    # Coefficients
    CL = sign(dot(lift, Lhat)) * norm(lift) / (q*Aref)
    CD = sign(dot(drag, Dhat)) * norm(drag) / (q*Aref)
    CY = sign(dot(side, Shat)) * norm(side) / (q*Aref)
    
    Cl = sign(dot(roll, lhat)) * norm(roll) / (q*Aref*cref)
    Cm = sign(dot(pitch, mhat)) * norm(pitch) / (q*Aref*cref)
    Cn = sign(dot(yaw, nhat)) * norm(yaw) / (q*Aref*cref)

    # Non-dimensional force and moment distributions
    if !isnothing(distributions)

        distrs = distributions[end]

        spanposition = distrs.spanposition
        cl = distrs.lift_distribution / (q*Aref/bref)
        cd = distrs.drag_distribution / (q*Aref/bref)
        cy = distrs.side_distribution / (q*Aref/bref)

        distributions[end] = (; spanposition, cl, cd, cy)
    end

    # Outputs
    return (;   CL, CD, CY, Cl, Cm, Cn,
                Dhat, Shat, Lhat, lhat, mhat, nhat,
                q, Aref, bref, cref)
end

################################################################################
# FORCE FIELDS
################################################################################
"""
    calcfield_Fkj!(out::Vector, body::AbstractBody,
                         areas::Vector, normals::Matrix, Cps::Vector,
                         Uinf::Number, rho::Number;
                         fieldname="Fkj")

Calculate the force of each element using the Kutta-Joukowski theorem,
``F = \\rho \\mathbf{u} \\times \\mathbf{l} \\Gamma``. 
``F`` is saved as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_Fkj!(out::AbstractMatrix, ll::LiftingLine, rho::Number;
                                                addfield=true, fieldname="Fkj")

    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==ll.nelements ""*
        "Invalid `out` matrix."*
        " Expected size $((3, ll.nelements)); got $(size(out))."

    for ei in 1:ll.nelements                # Iterate over horseshoes
        for (i, ip1) in (                   # Iterate over bound vortices
                            # (1, 2),           # A - Ap
                              (2, 3),           # B - A
                            # (3, 4)            # Bp - B
                        )

            dl1 = ll.horseshoes[1, ip1, ei] - ll.horseshoes[1, i, ei]
            dl2 = ll.horseshoes[2, ip1, ei] - ll.horseshoes[2, i, ei]
            dl3 = ll.horseshoes[3, ip1, ei] - ll.horseshoes[3, i, ei]

            out[1, ei] += rho * ll.Gammas[ei] * (ll.Us[2, ei]*dl3 - ll.Us[3, ei]*dl2)
            out[2, ei] += rho * ll.Gammas[ei] * (ll.Us[3, ei]*dl1 - ll.Us[1, ei]*dl3)
            out[3, ei] += rho * ll.Gammas[ei] * (ll.Us[1, ei]*dl2 - ll.Us[2, ei]*dl1)

        end
    end

    # Save field in lifting line
    if addfield
        add_field(ll, fieldname, "vector", eachcol(out), "cell")
    end

    return out
end


"""
    calcfield_Fkj(args...; optargs...)

Similar to [`calcfield_Fkj!`](@ref) but automatically pre-allocating `out` if it
hasn't been pre-allocated yet
"""
function calcfield_Fkj(ll::LiftingLine{R}, args...; fieldname="Fkj", optargs...) where {R}
    
    out = zeros(R, 3, ll.nelements)

    return calcfield_Fkj!(out, ll, args...; optargs...)

end


"""
    calcfield_fkj(args...; optargs...)

Similar to [`calcfield_fkj!`](@ref) but automatically pre-allocating `out` if it
hasn't been pre-allocated yet
"""
function calcfield_fkj(ll::LiftingLine{R}, args...; F_fieldname="Fkj", fieldname="fkj", optargs...) where {R}

    return calcfield_f!(ll, args...; F_fieldname, fieldname, optargs...)

end







"""
    calcfield_cl!(out::Vector, liftingline::LiftingLine; fieldname="cl")

Calculate effective lift coefficient cl on the swept section (C_𝐿Λ in Goates 
2022 nomenclature). ``cl`` is saved as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_cl!(out::AbstractVector, 
                        ll::LiftingLine;
                        addfield=true, fieldname="cl")

    # Error cases
    @assert length(out)==ll.nelements ""*
        "Invalid `out` vector."*
        " Expected size $((ll.nelements, )); got $(size(out))."

    for (ei, (element, aoa)) in enumerate(zip(ll.elements, ll.aoas))  # Iterate over stripwise elements

        sweep = calc_sweep(ll, ei)

        # Calculate swept sectional cl (C_𝐿Λ in Goates 2022, Eq. (28))
        clΛ = calc_sweptcl(element, sweep, aoa, view(ll.elements_settings, ei, :)...)

        out[ei] = clΛ

    end

    # Save field in lifting line
    if addfield
        add_field(ll, fieldname, "scalar", out, "cell")
    end

    return out
end

"""
    calcfield_cl(args...; optargs...)

Similar to [`calcfield_cl!`](@ref) but automatically pre-allocating `out` if it
hasn't been pre-allocated yet
"""
function calcfield_cl(ll::LiftingLine{R}, args...; fieldname="cl", optargs...) where {R}
    
    out = zeros(R, ll.nelements)

    return calcfield_cl!(out, ll, args...; optargs...)

end

"""
    calcfield_cd!(out::Vector, liftingline::LiftingLine; fieldname="cl")

Calculate drag coefficient cd. ``cd`` is saved as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).

**NOTE:** cd is most accurate when `calcfield_cd!` is called on the velocity 
    induced by the effective horseshoes as follows:

```julia
ll.Us .= Uinfs
selfUind!(ll, ll.Us)

calcfield_cd(ll)
```
"""
function calcfield_cd!(out::AbstractVector, 
                        ll::LiftingLine;
                        addfield=true, fieldname="cd")

    # Error cases
    @assert length(out)==ll.nelements ""*
        "Invalid `out` vector."*
        " Expected size $((ll.nelements, )); got $(size(out))."

    for (ei, (element, aoa)) in enumerate(zip(ll.elements, ll.aoas))  # Iterate over stripwise elements
        
        out[ei] = calc_cd(element, aoa, view(ll.elements_settings, ei, :)...)

        # NOTE: Goates 2022 JoA, Sec. V.E, recommends using the effective swept
        #       AOA, but we are getting too high of a cd. Hence, here we switch
        #       to the unswept AOA, assumming that Us is the unswept velocity.

        # # Calculate unswept AOA
        # aoa_unswept = calc_aoa(ll, ll.Us, ei)

        # out[ei] = calc_cd(element, aoa_unswept)

    end

    # Save field in lifting line
    if addfield
        add_field(ll, fieldname, "scalar", out, "cell")
    end

    return out
end

"""
    calcfield_cd(args...; optargs...)

Similar to [`calcfield_cd!`](@ref) but automatically pre-allocating `out` if it
hasn't been pre-allocated yet
"""
function calcfield_cd(ll::LiftingLine{R}, args...; fieldname="cd", optargs...) where {R}
    
    out = zeros(R, ll.nelements)

    return calcfield_cd!(out, ll, args...; optargs...)

end

"""
    calcfield_cm!(out::Vector, liftingline::LiftingLine; fieldname="cm")

Calculate effective pitching moment cm on the swept section (C_mc/4Λ in Goates 
2022 nomenclature). ``cm`` is saved as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_cm!(out::AbstractVector, 
                        ll::LiftingLine;
                        addfield=true, fieldname="cm")

    # Error cases
    @assert length(out)==ll.nelements ""*
        "Invalid `out` vector."*
        " Expected size $((ll.nelements, )); got $(size(out))."

    for (ei, (element, aoa)) in enumerate(zip(ll.elements, ll.aoas))  # Iterate over stripwise elements

        sweep = calc_sweep(ll, ei)

        # Calculate swept sectional cm (C_mc/4Λ in Goates 2022, Eq. (30))
        cmΛ = calc_sweptcm(element, sweep, aoa, view(ll.elements_settings, ei, :)...)

        out[ei] = cmΛ

    end

    # Save field in lifting line
    if addfield
        add_field(ll, fieldname, "scalar", out, "cell")
    end

    return out
end


"""
    calcfield_cm(args...; optargs...)

Similar to [`calcfield_cm!`](@ref) but automatically pre-allocating `out` if it
hasn't been pre-allocated yet
"""
function calcfield_cm(ll::LiftingLine{R}, args...; fieldname="cm", optargs...) where {R}
    
    out = zeros(R, ll.nelements)

    return calcfield_cm!(out, ll, args...; optargs...)

end

"""
    calcfield_F!(out::Vector, liftingLine::LiftingLine,
                         cls::Vector, cds::Vector, 
                         rho::Number;
                         fieldname="F")

Calculate loading distribution (force per unit span) using the Kutta-Joukowski 
theorem,
``f = \\rho \\mathbf{u} \\times \\hat{\\mathbf{y}} \\Gamma``. 
``f`` is saved as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).

NOTE: Make sure that the following has been run on the lifting line solution to
calculate the effective swept velocity before calling this function:

```julia
ll.Us .= Uinfs
Uind!(ll, ll.midpoints, ll.Us)
calc_UΛs!(ll, ll.Us)

calcfield_F(liftingline, rho)
```

This is equivalent to computing the velocity on the original horseshoes (not the
effective ones) which is then converted to effective velocity on the swept 
sections.
"""
function calcfield_F!(out::AbstractMatrix, ll::LiftingLine{R}, 
                        cls::AbstractVector, cds::AbstractVector, rho::Number;
                                                addfield=true, fieldname="f") where R

    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==ll.nelements ""*
        "Invalid `out` matrix."*
        " Expected size $((3, ll.nelements)); got $(size(out))."

    @assert length(cls)==ll.nelements ""*
        "Invalid `cls` vector."*
        " Expected size $((ll.nelements, )); got $(size(cls))."

    @assert length(cds)==ll.nelements ""*
        "Invalid `cds` vector."*
        " Expected size $((ll.nelements, )); got $(size(cds))."

    for ei in 1:ll.nelements                # Iterate over stripwise elements

        # Velocity
        U1 = ll.Us[1, ei]
        U2 = ll.Us[2, ei]
        U3 = ll.Us[3, ei]
        magU = sqrt(U1^2 + U2^2 + U3^2)

        # Lifting filament length
        dl1 = ll.horseshoes[1, 3, ei] - ll.horseshoes[1, 2, ei]
        dl2 = ll.horseshoes[2, 3, ei] - ll.horseshoes[2, 2, ei]
        dl3 = ll.horseshoes[3, 3, ei] - ll.horseshoes[3, 2, ei]

        # Span length of this element
        ds = abs(dl1*ll.spans[1, ei] + dl2*ll.spans[2, ei] + dl3*ll.spans[3, ei])

        # Dynamic pressure
        q = 0.5*rho*magU^2

        # Integrated lift and drag over this span section
        L = cls[ei] * q * ll.chords[ei] * ds
        D = cds[ei] * q * ll.chords[ei] * ds

        # Drag direction: direction of freestream
        hatD1 = U1
        hatD2 = U2
        hatD3 = U3

        aux = sqrt(hatD1^2 + hatD2^2 + hatD3^2)

        hatD1 /= aux
        hatD2 /= aux
        hatD3 /= aux

        # Lift direction: hatD x bound vortex
        hatL1 = hatD2*ll.lines[3, ei] - hatD3*ll.lines[2, ei]
        hatL2 = hatD3*ll.lines[1, ei] - hatD1*ll.lines[3, ei]
        hatL3 = hatD1*ll.lines[2, ei] - hatD2*ll.lines[1, ei]

        aux = sqrt(hatL1^2 + hatL2^2 + hatL3^2)

        hatL1 /= aux
        hatL2 /= aux
        hatL3 /= aux

        # Add lift and drag
        out[1, ei] += L*hatL1 + D*hatD1
        out[2, ei] += L*hatL2 + D*hatD2
        out[3, ei] += L*hatL3 + D*hatD3
        
    end

    # Save field in lifting line
    if addfield
        add_field(ll, fieldname, "vector", eachcol(out), "cell")
    end

    return out
end


"""
    calcfield_F(args...; optargs...)

Similar to [`calcfield_F!`](@ref) but automatically pre-allocating `out` if it
hasn't been pre-allocated yet
"""
function calcfield_F(ll::LiftingLine{R}, args...; 
                        cl_fieldname="cl", cd_fieldname="cd", 
                        fieldname="F", optargs...) where {R}

    # Error cases
    @assert check_field(ll, cl_fieldname) ""*
        "Field $(cl_fieldname) not found;"*
        " Please run `calcfield_cl(args...; fieldname=$(cl_fieldname), optargs...)`"

    @assert check_field(ll, cd_fieldname) ""*
        "Field $(cd_fieldname) not found;"*
        " Please run `calcfield_cd(args...; fieldname=$(cd_fieldname), optargs...)`"

    cls = get_field(ll, cl_fieldname)["field_data"]
    cds = get_field(ll, cd_fieldname)["field_data"]
    
    out = zeros(R, 3, ll.nelements)

    return calcfield_F!(out, ll, cls, cds, args...; fieldname, optargs...)

end


function calcfield_Ftot(ll::LiftingLine{R}, args...; optargs...) where R 
    return calcfield_Ftot!(zeros(R, 3), ll, args...; optargs...)
end

function calcfield_LDS(ll::LiftingLine{R}, args...; optargs...) where R 
    return calcfield_LDS!(zeros(R, 3, 3), ll, args...; optargs...)
end


"""
    calcfield_f!(out::Vector, liftingline::LiftingLine, Fs::Matrix;
                         fieldname="f")

Calculate loading distribution (force per unit span).

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_f!(out::AbstractMatrix, ll::LiftingLine, Fs::AbstractMatrix;
                                                addfield=true, fieldname="f")

    # Error cases
    @assert size(out, 1)==3 && size(out, 2)==ll.nelements ""*
        "Invalid `out` matrix."*
        " Expected size $((3, ll.nelements)); got $(size(out))."

    @assert size(Fs, 1)==3 && size(Fs, 2)==ll.nelements ""*
        "Invalid `Fs` matrix."*
        " Expected size $((3, ll.nelements)); got $(size(Fs))."

    for ei in 1:ll.nelements                # Iterate over stripwise elements

        # Lifting filament length
        dl1 = ll.horseshoes[1, 3, ei] - ll.horseshoes[1, 2, ei]
        dl2 = ll.horseshoes[2, 3, ei] - ll.horseshoes[2, 2, ei]
        dl3 = ll.horseshoes[3, 3, ei] - ll.horseshoes[3, 2, ei]

        # Span length of this element
        ds = abs(dl1*ll.spans[1, ei] + dl2*ll.spans[2, ei] + dl3*ll.spans[3, ei])

        for i in 1:3
            out[i, ei] += Fs[i, ei] / ds
        end

    end

    # Save field in lifting line
    if addfield
        add_field(ll, fieldname, "vector", eachcol(out), "cell")
    end

    return out
end


"""
    calcfield_f(args...; optargs...)

Similar to [`calcfield_f!`](@ref) but automatically pre-allocating `out` if it
hasn't been pre-allocated yet
"""
function calcfield_f(ll::LiftingLine{R}, args...; F_fieldname="F", fieldname="f", optargs...) where {R}

    # Error case
    @assert check_field(ll, F_fieldname) ""*
        "Field $(F_fieldname) not found;"*
        " Please run `calcfield_F(args...; fieldname=$(F_fieldname), optargs...)`"

    Fs = hcat(get_field(ll, F_fieldname)["field_data"]...)
    
    out = zeros(R, 3, ll.nelements)

    return calcfield_f!(out, ll, Fs, args...; fieldname, optargs...)

end









################################################################################
# MOMENT FIELDS
################################################################################
"""
    calcfield_Mtot!(out::AbstractVector, liftingline::LiftingLine,
                                X0::AbstractVector, midpoints::AbstractMatrix,
                                Fs::AbstractMatrix,
                                cmΛs::AbstractVector, chords::AbstractVector;
                                fieldname="Mtot", addfield=true)

Calculate the integrated moment of this lifting line (which is a 
three-dimensional vector) with respect to the center `X0`.
This is calculated from the force and position of each element given in `Fs` at
`midpoints`, respectively, in addition to the pitching moment of the stripwise
element, and saved as a field named `fieldname`.

The field is calculated in-place and added to `out`.

NOTE: Make sure that the following has been run on the lifting line solution to
calculate the effective swept velocity before calling this function:

```julia
ll.Us .= Uinfs
Uind!(ll, ll.midpoints, ll.Us)
calc_UΛs!(ll, ll.Us)

calcfield_Mtot(liftingline, X0, rho)
```

This is equivalent to computing the velocity on the original horseshoes (not the
effective ones) which is then converted to effective velocity on the swept 
sections.
"""
function calcfield_Mtot!(out::AbstractVector, ll::LiftingLine,
                            X0::AbstractVector, midpoints::AbstractMatrix,
                            Fs::AbstractMatrix, 
                            cmΛs::AbstractVector, chords::AbstractVector,
                            rho::Number;
                            fieldname="Mtot", addfield=true)
    # Error case
    @assert length(out)==3 ""*
        "Invalid `out` vector. Expected length 3; got $(length(out))."
    @assert length(X0)==3 ""*
        "Invalid `X0` vector. Expected length 3; got $(length(X0))."
    @assert size(midpoints, 1)==3 && size(midpoints, 2)==ll.nelements ""*
        "Invalid `midpoints` matrix."*
        " Expected size $((3, ll.nelements)); got $(size(midpoints))."
    @assert size(Fs, 1)==3 && size(Fs, 2)==ll.nelements ""*
        "Invalid `Fs` matrix."*
        " Expected size $((3, ll.nelements)); got $(size(Fs))."
    @assert length(cmΛs)==ll.nelements ""*
        "Invalid `cmΛs` vector. Expected length $(ll.nelements); got $(length(cmΛs))."
    @assert length(chords)==ll.nelements ""*
        "Invalid `chords` vector. Expected length $(ll.nelements); got $(length(chords))."

    # Calculate moment from force component
    for (X, F) in zip(eachcol(midpoints), eachcol(Fs))
        out[1] += (X[2] - X0[2])*F[3] - (X[3] - X0[3])*F[2]
        out[2] += (X[3] - X0[3])*F[1] - (X[1] - X0[1])*F[3]
        out[3] += (X[1] - X0[1])*F[2] - (X[2] - X0[2])*F[1]
    end

    # Calculate moment from stripwise pitching moment coefficient
    for ei in 1:ll.nelements

        sweep = calc_sweep(ll, ei)

        # Velocity
        U1 = ll.Us[1, ei]
        U2 = ll.Us[2, ei]
        U3 = ll.Us[3, ei]
        magU = sqrt(U1^2 + U2^2 + U3^2)

        # Lifting filament length
        dl1 = ll.horseshoes[1, 3, ei] - ll.horseshoes[1, 2, ei]
        dl2 = ll.horseshoes[2, 3, ei] - ll.horseshoes[2, 2, ei]
        dl3 = ll.horseshoes[3, 3, ei] - ll.horseshoes[3, 2, ei]

        # Span length of this element
        ds = abs(dl1*ll.spans[1, ei] + dl2*ll.spans[2, ei] + dl3*ll.spans[3, ei])

        # Dynamic pressure
        q = 0.5*rho*magU^2

        # Integrated pitching moment over this span section
        area = ll.chords[ei] * ds
        M = cmΛs[ei] * q * ll.chords[ei]*cosd(sweep) * area

        # Add pitching moment as aligned with the lifting line
        out[1] += M*ll.lines[1, ei]
        out[2] += M*ll.lines[2, ei]
        out[3] += M*ll.lines[3, ei]
    end


    # Save field in lifting line
    if addfield
        add_field(ll, fieldname, "vector", out, "system")
    end

    return out
end

"""
    calcfield_Mtot!(out, liftingline, X0; F_fieldname="F", cm_fieldname="cm",
                                                                    optargs...)

Calculate the integrated moment of this body (which is a three-dimensional
vector) with respect to the center `X0`.
This is calculated from the force field `F_fieldname` and sectional pitching 
moment coefficient `cm_fieldname` and saved as a field named `fieldname`.

The field is calculated in-place and added to `out`.
"""
function calcfield_Mtot!(out, ll::LiftingLine, X0, rho; 
                            F_fieldname="F", cm_fieldname="cm", optargs...)
    # Error case
    @assert check_field(ll, F_fieldname) ""*
        "Field $(F_fieldname) not found;"*
        " Please run `calcfield_F(args...; fieldname=$(F_fieldname), optargs...)`"

    @assert isnothing(cm_fieldname) || check_field(ll, cm_fieldname) ""*
        "Field $(cm_fieldname) not found;"*
        " Please run `calcfield_cm(args...; fieldname=$(cm_fieldname), optargs...)`"

    Fs = hcat(get_field(ll, F_fieldname)["field_data"]...)

    if isnothing(cm_fieldname)
        cmΛs = zeros(ll.nelements)
    else
        cmΛs = get_field(ll, cm_fieldname)["field_data"]
    end

    return calcfield_Mtot!(out, ll, X0, ll.midpoints, Fs, cmΛs, ll.chords, rho; optargs...)
end


function calcfield_Mtot(ll::LiftingLine{R}, args...; optargs...) where R 
    return calcfield_Mtot!(zeros(R, 3), ll, args...; optargs...)
end


"""
    calcfield_lmn!(out::Matrix, liftingline::LiftingLine,
                    Mtot:Vector, lhat::Vector, mhat::Vector, nhat::Vector)

Decompose the integrated moment `Mtot` as rolling, pitching, and yawing moments 
according to the orthonormal basis `lhat`, `mhat`, `nhat`, repsectively.

`out[:, 1]` is the rolling moment vector and is saved as the field "Mroll".
`out[:, 2]` is the pitching moment vector and is saved as the field "Mpitch".
`out[:, 3]` is the yawing moment vector and is saved as the field "Myaw".

The field is calculated in-place on `out`.
"""
function calcfield_lmn!(out::AbstractMatrix, ll::LiftingLine,
                        Mtot::AbstractVector,
                        lhat::AbstractVector, mhat::AbstractVector,
                        nhat::AbstractVector;
                        addfield=true)
    # Error case
    @assert size(out, 1)==3 && size(out, 2)==3 ""*
        "Invalid `out` matrix. Expected size $((3, 3)); got $(size(out))."
    @assert length(Mtot)==3 ""*
        "Invalid `Mtot` vector. Expected length 3; got $(length(Mtot))."
    @assert abs(norm(lhat) - 1) <= eps(10.0) ""*
        "lhat=$(lhat) is not a unitary vector"
    @assert abs(norm(mhat) - 1) <= eps(10.0) ""*
        "mhat=$(mhat) is not a unitary vector"
    @assert abs(norm(nhat) - 1) <= eps(10.0) ""*
        "nhat=$(nhat) is not a unitary vector"

    # Project Mtot in each direction
    Ml = dot(Mtot, lhat)
    Mm = dot(Mtot, mhat)
    Mn = dot(Mtot, nhat)
    
    for i in 1:3
        out[i, 1] += Ml * lhat[i]
        out[i, 2] += Mm * mhat[i]
        out[i, 3] += Mn * nhat[i]
    end

    # Save field in liftingline
    if addfield
        add_field(ll, "Mroll", "vector", view(out, :, 1), "system")
        add_field(ll, "Mpitch", "vector", view(out, :, 2), "system")
        add_field(ll, "Myaw", "vector", view(out, :, 3), "system")
    end

    return out
end


"""
    calcfield_lmn!(out, liftingline, lhat, mhat, nhat; Mtot_fieldname="Mtot",
                                                                    optargs...)

Decompose the integrated moment `Mtot` as rolling, pitching, and yawing moments 
according to the orthonormal basis `lhat`, `mhat`, `nhat`, repsectively.
This is calculated from the total moment field `Mtot_fieldname`.

The field is calculated in-place on `out`.
"""
function calcfield_lmn!(out, ll::LiftingLine, lhat, mhat, nhat; Mtot_fieldname="Mtot",
                            optargs...)
    # Error case
    @assert check_field(ll, Mtot_fieldname) ""*
        "Field $(Mtot_fieldname) not found;"*
        " Please run `calcfield_Mtot(args...; fieldname=$(Mtot_fieldname), optargs...)`"

    Mtot = get_field(ll, Mtot_fieldname)["field_data"]

    return calcfield_lmn!(out, ll, Mtot, lhat, mhat, nhat; optargs...)
end

"""
    calcfield_lmn!(out, liftingline, lhat, mhat; optargs...)

`nhat` is calculated automatically from `lhat` and `mhat`,
"""
function calcfield_lmn!(out, ll::LiftingLine, lhat, mhat; optargs...)
    return calcfield_lmn!(out, ll, lhat, mhat, cross(lhat, mhat); optargs...)
end

function calcfield_lmn(ll::LiftingLine{R}, args...; optargs...) where R 
    return calcfield_lmn!(zeros(R, 3, 3), ll, args...; optargs...)
end