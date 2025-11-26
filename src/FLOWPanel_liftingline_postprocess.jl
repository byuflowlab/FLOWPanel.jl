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
                         fieldname="F")

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
    
    if fieldname in ll.fields
        out = get_field(ll, fieldname)["field_data"]
        out .= 0
    else
        out = zeros(R, 3, ll.nelements)
    end

    return calcfield_Fkj!(out, ll, args...; optargs...)

end




"""
    calcfield_fkj!(out::Vector, body::AbstractBody,
                         areas::Vector, normals::Matrix, Cps::Vector,
                         Uinf::Number, rho::Number;
                         fieldname="F")

Calculate loading distribution (force per unit span) using the Kutta-Joukowski 
theorem,
``f = \\rho \\mathbf{u} \\times \\hat{\\mathbf{y}} \\Gamma``. 
``f`` is saved as a field named `fieldname`.

The field is calculated in-place and added to `out` (hence, make sure that `out`
starts with all zeroes).
"""
function calcfield_fkj!(out::AbstractMatrix, ll::LiftingLine, rho::Number;
                                                addfield=true, fieldname="fkj")

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

            # Filament length
            dl1 = ll.horseshoes[1, ip1, ei] - ll.horseshoes[1, i, ei]
            dl2 = ll.horseshoes[2, ip1, ei] - ll.horseshoes[2, i, ei]
            dl3 = ll.horseshoes[3, ip1, ei] - ll.horseshoes[3, i, ei]

            # Velocity
            U1 = ll.Us[1, ei]
            U2 = ll.Us[2, ei]
            U3 = ll.Us[3, ei]
            magU = sqrt(U1^2 + U2^2 + U3^2)

            # Calculate tranversal length (counter-projection of U): dy = dl - (dl⋅hatU)hatU
            dldotU = dl1*U1 + dl2*U2 + dl3*U3
            dy1 = dl1 - dldotU*U1/magU^2
            dy2 = dl2 - dldotU*U2/magU^2
            dy3 = dl3 - dldotU*U3/magU^2
            magdy = sqrt(dy1^2 + dy2^2 + dy3^2)

            out[1, ei] += rho * ll.Gammas[ei] * (U2*dy3 - U3*dy2)/magdy
            out[2, ei] += rho * ll.Gammas[ei] * (U3*dy1 - U1*dy3)/magdy
            out[3, ei] += rho * ll.Gammas[ei] * (U1*dy2 - U2*dy1)/magdy

        end
    end

    # Save field in lifting line
    if addfield
        add_field(ll, fieldname, "vector", eachcol(out), "cell")
    end

    return out
end


"""
    calcfield_fkj(args...; optargs...)

Similar to [`calcfield_fkj!`](@ref) but automatically pre-allocating `out` if it
hasn't been pre-allocated yet
"""
function calcfield_fkj(ll::LiftingLine{R}, args...; fieldname="fkj", optargs...) where {R}
    
    if fieldname in ll.fields
        out = get_field(ll, fieldname)["field_data"]
        out .= 0
    else
        out = zeros(R, 3, ll.nelements)
    end

    return calcfield_fkj!(out, ll, args...; optargs...)

end