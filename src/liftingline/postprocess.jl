#=##############################################################################
# DESCRIPTION
    Methods for postprocessing solver results of non-linear lifting line method

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################


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

    return calcfield_f(ll, args...; F_fieldname, fieldname, optargs...)

end







"""
    calcfield_cl!(out::Vector, liftingline::LiftingLine; fieldname="cl")

Calculate effective lift coefficient cl on the swept section (C_𝐿Λ in Goates 
2022 nomenclature). ``cl`` is saved as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_cl!(out::AbstractVector, 
                        ll::LiftingLine{<:Number, <:SimpleAirfoil};
                        addfield=true, fieldname="cl")

    # Error cases
    @assert length(out)==ll.nelements ""*
        "Invalid `out` vector."*
        " Expected size $((ll.nelements, )); got $(size(out))."

    for (ei, (element, aoa)) in enumerate(zip(ll.elements, ll.aoas))  # Iterate over stripwise elements

        sweep = calc_sweep(ll, ei)

        # Calculate swept sectional cl (C_𝐿Λ in Goates 2022, Eq. (28))
        clΛ = calc_sweptcl(element, sweep, aoa)

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
                        ll::LiftingLine{<:Number, <:SimpleAirfoil};
                        addfield=true, fieldname="cd")

    # Error cases
    @assert length(out)==ll.nelements ""*
        "Invalid `out` vector."*
        " Expected size $((ll.nelements, )); got $(size(out))."

    for (ei, (element, aoa)) in enumerate(zip(ll.elements, ll.aoas))  # Iterate over stripwise elements
        
        # out[ei] = calc_cd(element, aoa)

        # NOTE: Goates 2022 JoA, Sec. V.E, recommends using the effective swept
        #       AOA, but we are getting too high of a cd. Hence, here we switch
        #       to the unswept AOA, assumming that Us is the unswept velocity.

        # Calculate unswept AOA
        aoa_unswept = calc_aoa(ll, ll.Us, ei)

        out[ei] = calc_cd(element, aoa_unswept)

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
                        ll::LiftingLine{<:Number, <:SimpleAirfoil};
                        addfield=true, fieldname="cm")

    # Error cases
    @assert length(out)==ll.nelements ""*
        "Invalid `out` vector."*
        " Expected size $((ll.nelements, )); got $(size(out))."

    for (ei, (element, aoa)) in enumerate(zip(ll.elements, ll.aoas))  # Iterate over stripwise elements

        sweep = calc_sweep(ll, ei)

        # Calculate swept sectional cm (C_mc/4Λ in Goates 2022, Eq. (30))
        cmΛ = calc_sweptcm(element, sweep, aoa)

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