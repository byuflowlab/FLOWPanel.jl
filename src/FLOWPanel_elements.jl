#=##############################################################################
# DESCRIPTION
    Definition of panel elements

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Jul 2018
  * License     : MIT License
=###############################################################################

abstract type AbstractElement end
struct ConstantSource <: AbstractElement end
struct ConstantDoublet <: AbstractElement end
struct VortexRing <: AbstractElement end
struct ConstantVortexSheet <: AbstractElement end
struct UniformVortexSheet <: AbstractElement end

################################################################################
# SOURCE ELEMENTS
################################################################################
"""
Computes the velocity induced by a panel of vertices `nodes[:, panel]` and
constant strength source `strength` on the targets `targets`. It adds the
velocity at the i-th target to out[i].

Implementation of equations in Hess & Smith (1967).
"""
function U_constant_source(nodes::Arr1, panel,
                              strength::Number,
                              targets::Arr2, out::Arr3;
                              dot_with=nothing,
                              offset=1e-8
                          ) where{T1, Arr1<:AbstractArray{T1,2},
                                  T2, Arr2<:AbstractArray{T2,2},
                                  T3, Arr3<:AbstractArray{T3}}

    nt = size(targets, 2)                   # Number of targets
    no = dot_with!=nothing ? length(out) : size(out, 2) # Number of outputs
    nn = length(panel)                      # Number of nodes

    if no!=nt
        error("Invalid `out` argument. Expected size $(nt), got $(no).")
    end

    # Tangent, oblique, and normal vectors
    # t1, t2, t3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
    # o1, o2, o3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
    # n1, n2, n3 = gt._calc_n1(nodes, panel), gt._calc_n2(nodes, panel), gt._calc_n3(nodes, panel)

    # Here we flip the oblique and tangent vector since both Katz & Plotkin and
    # Hess & Smith define their panels opposite to the right-hand rule (while
    # GeometricTools defines the normal following the right-hand rule)
    t1, t2, t3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
    o1, o2, o3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
    n1, n2, n3 = -gt._calc_n1(nodes, panel), -gt._calc_n2(nodes, panel), -gt._calc_n3(nodes, panel)

    # Panel local coordinate system
    # NOTE: normal is pointing out of the body, which differs from Katz and Plotkin
    # @inbounds O = nodes[1]                         # Origin
    @inbounds Oi = panel[1]                         # Index of node that is the origin
    # xhat, yhat, zhat = t, o, n         # Unit vectors
    # Oaxis = hcat(xhat, yhat, zhat)'    # Transformation matrix

    # Iterate over targets
    for ti in 1:nt

        # Target position in panel coordinate system
        # X = Oaxis*(targets[:, ti]-O)
        @inbounds begin
            x = t1*(targets[1,ti]-nodes[1,Oi]) + t2*(targets[2,ti]-nodes[2,Oi]) + t3*(targets[3,ti]-nodes[3,Oi])
            y = o1*(targets[1,ti]-nodes[1,Oi]) + o2*(targets[2,ti]-nodes[2,Oi]) + o3*(targets[3,ti]-nodes[3,Oi])
            z = n1*(targets[1,ti]-nodes[1,Oi]) + n2*(targets[2,ti]-nodes[2,Oi]) + n3*(targets[3,ti]-nodes[3,Oi])
        end

        V1, V2, V3 = zero(T3), zero(T3), zero(T3)
        dtheta = 2*pi

        nR0::Int = 0

        @simd for i in 1:nn

            @inbounds begin
                pi, pj = panel[i], panel[i%nn + 1]

                # Convert nodes to panel coordinate system
                xi = t1*(nodes[1,pi]-nodes[1,Oi]) + t2*(nodes[2,pi]-nodes[2,Oi]) + t3*(nodes[3,pi]-nodes[3,Oi])
                yi = o1*(nodes[1,pi]-nodes[1,Oi]) + o2*(nodes[2,pi]-nodes[2,Oi]) + o3*(nodes[3,pi]-nodes[3,Oi])
                zi = n1*(nodes[1,pi]-nodes[1,Oi]) + n2*(nodes[2,pi]-nodes[2,Oi]) + n3*(nodes[3,pi]-nodes[3,Oi])
                xj = t1*(nodes[1,pj]-nodes[1,Oi]) + t2*(nodes[2,pj]-nodes[2,Oi]) + t3*(nodes[3,pj]-nodes[3,Oi])
                yj = o1*(nodes[1,pj]-nodes[1,Oi]) + o2*(nodes[2,pj]-nodes[2,Oi]) + o3*(nodes[3,pj]-nodes[3,Oi])
                zj = n1*(nodes[1,pj]-nodes[1,Oi]) + n2*(nodes[2,pj]-nodes[2,Oi]) + n3*(nodes[3,pj]-nodes[3,Oi])
            end

            dij = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
            ri = sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)
            rj = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2)

            Qij = log( (ri+rj+dij)/(ri+rj-dij + offset) )

            Sij = (yj-yi)/dij
            Cij = (xj-xi)/dij

            siji = (xi-x)*Cij + (yi-y)*Sij
            sijj = (xj-x)*Cij + (yj-y)*Sij
            Rij = (x-xi)*Sij - (y-yi)*Cij

            Jij = atan( Rij*abs(z)*( ri*sijj - rj*siji ) , ri*rj*Rij^2 + z^2*sijj*siji )

            V1 -= Sij*Qij
            V2 += Cij*Qij
            V3 -= Jij

            dtheta *= Rij>=0
            nR0 += Rij==0
            # nR0 += (abs(Rij) < offset)
        end

        V3 += dtheta
        V3 *= sign(z)        # Isn't this sign already accounted for in atan2?
        V3 *= !(nR0>1)       # Singularity fix of any z position aligned with node


        # NOTE: Katz and Plotkin's potential differs from Hess and Smith's by
        #       this factor
        V1 *= 1/(4*pi)
        V2 *= 1/(4*pi)
        V3 *= 1/(4*pi)

        if dot_with!=nothing
            @inbounds out[ti] += strength*(V1*t1 + V2*o1 + V3*n1)*dot_with[1,ti]
            @inbounds out[ti] += strength*(V1*t2 + V2*o2 + V3*n2)*dot_with[2,ti]
            @inbounds out[ti] += strength*(V1*t3 + V2*o3 + V3*n3)*dot_with[3,ti]
        else
            @inbounds out[1, ti] += strength*(V1*t1 + V2*o1 + V3*n1)
            @inbounds out[2, ti] += strength*(V1*t2 + V2*o2 + V3*n2)
            @inbounds out[3, ti] += strength*(V1*t3 + V2*o3 + V3*n3)
        end

    end
end



"""
Computes the potential induced by a panel of vertices `nodes[:, panel]` and
constant strength source `strength` on the targets `targets`. It adds the
potential at the i-th target to out[i].

Implementation of equations in Hess & Smith (1967), p. 53.
"""
function phi_constant_source(nodes::Arr1, panel,
                              strength::Number,
                              targets::Arr2, out::Arr3;
                              offset::Real=1e-8
                            ) where{T1, Arr1<:AbstractArray{T1,2},
                                    T2, Arr2<:AbstractArray{T2,2},
                                    T3, Arr3<:AbstractArray{T3}}

    nt = size(targets, 2)                   # Number of targets
    no = length(out)                        # Number of outputs
    nn = length(panel)                      # Number of nodes

    if no!=nt
        error("Invalid `out` argument. Expected size $(nt), got $(no).")
    end

    # Tangent, oblique, and normal vectors
    # t1, t2, t3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
    # o1, o2, o3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
    # n1, n2, n3 = gt._calc_n1(nodes, panel), gt._calc_n2(nodes, panel), gt._calc_n3(nodes, panel)

    # Here we flip the oblique and tangent vector since both Katz & Plotkin and
    # Hess & Smith define their panels opposite to the right-hand rule (while
    # GeometricTools defines the normal following the right-hand rule)
    t1, t2, t3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
    o1, o2, o3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
    n1, n2, n3 = -gt._calc_n1(nodes, panel), -gt._calc_n2(nodes, panel), -gt._calc_n3(nodes, panel)

    # Panel local coordinate system
    # NOTE: normal is pointing out of the body, which differs from Katz and Plotkin
    @inbounds Oi = panel[1]                         # Index of node that is the origin
    # xhat, yhat, zhat = t, o, n         # Unit vectors
    # Oaxis = hcat(xhat, yhat, zhat)'    # Transformation matrix

    # Convert nodes to panel coordinate system
    # Pnodes = [Oaxis*(node-O) for node in nodes]

    # Iterate over targets
    for ti in 1:nt


        # Target position in panel coordinate system
        # X = Oaxis*(targets[:, ti]-O)
        @inbounds begin
            x = t1*(targets[1,ti]-nodes[1,Oi]) + t2*(targets[2,ti]-nodes[2,Oi]) + t3*(targets[3,ti]-nodes[3,Oi])
            y = o1*(targets[1,ti]-nodes[1,Oi]) + o2*(targets[2,ti]-nodes[2,Oi]) + o3*(targets[3,ti]-nodes[3,Oi])
            z = n1*(targets[1,ti]-nodes[1,Oi]) + n2*(targets[2,ti]-nodes[2,Oi]) + n3*(targets[3,ti]-nodes[3,Oi])
        end

        dtheta = 2*pi

        # Iterate over nodes
        for i in 1:nn

            @inbounds begin
                pi, pj = panel[i], panel[i%nn + 1]

                # Convert nodes to panel coordinate system
                xi = t1*(nodes[1,pi]-nodes[1,Oi]) + t2*(nodes[2,pi]-nodes[2,Oi]) + t3*(nodes[3,pi]-nodes[3,Oi])
                yi = o1*(nodes[1,pi]-nodes[1,Oi]) + o2*(nodes[2,pi]-nodes[2,Oi]) + o3*(nodes[3,pi]-nodes[3,Oi])
                zi = n1*(nodes[1,pi]-nodes[1,Oi]) + n2*(nodes[2,pi]-nodes[2,Oi]) + n3*(nodes[3,pi]-nodes[3,Oi])
                xj = t1*(nodes[1,pj]-nodes[1,Oi]) + t2*(nodes[2,pj]-nodes[2,Oi]) + t3*(nodes[3,pj]-nodes[3,Oi])
                yj = o1*(nodes[1,pj]-nodes[1,Oi]) + o2*(nodes[2,pj]-nodes[2,Oi]) + o3*(nodes[3,pj]-nodes[3,Oi])
                zj = n1*(nodes[1,pj]-nodes[1,Oi]) + n2*(nodes[2,pj]-nodes[2,Oi]) + n3*(nodes[3,pj]-nodes[3,Oi])
            end

            dij = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
            ri = sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)
            rj = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2)

            Qij = log( (ri+rj+dij)/(ri+rj-dij + offset) )

            Sij = (yj-yi)/dij
            Cij = (xj-xi)/dij

            siji = (xi-x)*Cij + (yi-y)*Sij
            sijj = (xj-x)*Cij + (yj-y)*Sij
            Rij = (x-xi)*Sij - (y-yi)*Cij

            Jij = atan( Rij*abs(z)*( ri*sijj - rj*siji ) , ri*rj*Rij^2 + z^2*sijj*siji )

            @inbounds out[ti] -= strength/(4*π) * (Rij*Qij + abs(z)*Jij)

            dtheta *= Rij>=0
        end

        @inbounds out[ti] += strength/(4*π) * abs(z)*dtheta

    end
end

# """
# Returns the potential induced by a panel of vertices `nodes[:, panel]` and constant
# strength source `strength` on the targets `targets`. It adds the potential at
# the i-th target to out[i].
#
# Implementation of equations in Katz and Plotkin Sec. 10.4.1.
# NOTE: THOSE EQUATIONS ARE WRONG. THEY GIVE THE CORRECT POTENTIAL ONLY ABOVE
#       THE PANEL, AND NON-SENSE UNDERNEATH IT.
# """
# function phi_constant_source(nodes::Arr1, panel,
#                               strength::Number,
#                               targets::Arr2, out::Arr3;
#                               offset::Real=1e-8
#                             ) where{T1, Arr1<:AbstractArray{T1,2},
#                                     T2, Arr2<:AbstractArray{T2,2},
#                                     T3, Arr3<:AbstractArray{T3}}
#
#     nt = size(targets, 2)                   # Number of targets
#     no = length(out)                        # Number of outputs
#     nn = length(panel)                      # Number of nodes
#
#     if no!=nt
#         error("Invalid `out` argument. Expected size $(nt), got $(no).")
#     end
#
#     # Tangent, oblique, and normal vectors
#     # t1, t2, t3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
#     # o1, o2, o3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
#     # n1, n2, n3 = gt._calc_n1(nodes, panel), gt._calc_n2(nodes, panel), gt._calc_n3(nodes, panel)
#
#     # Here we flip the oblique and tangent vector since both Katz & Plotkin and
#     # Hess & Smith define their panels opposite to the right-hand rule (while
#     # GeometricTools defines the normal following the right-hand rule)
#     t1, t2, t3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
#     o1, o2, o3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
#     n1, n2, n3 = -gt._calc_n1(nodes, panel), -gt._calc_n2(nodes, panel), -gt._calc_n3(nodes, panel)
#
#     # Panel local coordinate system
#     # NOTE: normal is pointing out of the body, which differs from Katz and Plotkin
#     @inbounds Oi = panel[1]                         # Index of node that is the origin
#     # xhat, yhat, zhat = t, o, n         # Unit vectors
#     # Oaxis = hcat(xhat, yhat, zhat)'    # Transformation matrix
#
#     # Convert nodes to panel coordinate system
#     # Pnodes = [Oaxis*(node-O) for node in nodes]
#
#     # Iterate over nodes
#     for i in 1:nn
#
#         @inbounds begin
#             pi, pj = panel[i], panel[i%nn + 1]
#
#             # Convert nodes to panel coordinate system
#             xi = t1*(nodes[1,pi]-nodes[1,Oi]) + t2*(nodes[2,pi]-nodes[2,Oi]) + t3*(nodes[3,pi]-nodes[3,Oi])
#             yi = o1*(nodes[1,pi]-nodes[1,Oi]) + o2*(nodes[2,pi]-nodes[2,Oi]) + o3*(nodes[3,pi]-nodes[3,Oi])
#             zi = n1*(nodes[1,pi]-nodes[1,Oi]) + n2*(nodes[2,pi]-nodes[2,Oi]) + n3*(nodes[3,pi]-nodes[3,Oi])
#             xj = t1*(nodes[1,pj]-nodes[1,Oi]) + t2*(nodes[2,pj]-nodes[2,Oi]) + t3*(nodes[3,pj]-nodes[3,Oi])
#             yj = o1*(nodes[1,pj]-nodes[1,Oi]) + o2*(nodes[2,pj]-nodes[2,Oi]) + o3*(nodes[3,pj]-nodes[3,Oi])
#             zj = n1*(nodes[1,pj]-nodes[1,Oi]) + n2*(nodes[2,pj]-nodes[2,Oi]) + n3*(nodes[3,pj]-nodes[3,Oi])
#         end
#
#         # Iterate over targets
#         @simd for ti in 1:nt
#
#             # Target position in panel coordinate system
#             # X = Oaxis*(targets[:, ti]-O)
#             @inbounds begin
#                 x = t1*(targets[1,ti]-nodes[1,Oi]) + t2*(targets[2,ti]-nodes[2,Oi]) + t3*(targets[3,ti]-nodes[3,Oi])
#                 y = o1*(targets[1,ti]-nodes[1,Oi]) + o2*(targets[2,ti]-nodes[2,Oi]) + o3*(targets[3,ti]-nodes[3,Oi])
#                 z = n1*(targets[1,ti]-nodes[1,Oi]) + n2*(targets[2,ti]-nodes[2,Oi]) + n3*(targets[3,ti]-nodes[3,Oi])
#             end
#
#             dij = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
#             mij = (yj - yi)/(xj - xi)
#             ri = sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)
#             rj = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2)
#             ei = (x - xi)^2 + (z-zi)^2
#             ej = (x - xj)^2 + (z-zj)^2
#             hi = (x - xi)*(y - yi)
#             hj = (x - xj)*(y - yj)
#
#             Pij = (x - xi)*(yj - yi) - (y - yi)*(xj - xi)
#             Qij = log( (ri+rj+dij)/(ri+rj-dij + offset) )
#             # Rij = atan(mij*ei-hi, z*ri) - atan(mij*ej-hj, z*rj)
#             Rij = atan( (mij*ei-hi) / (z*ri) ) - atan( (mij*ej-hj) / (z*rj) )
#
#             @inbounds out[ti] -= strength/(4*π) * (Pij/dij * Qij - abs(z)*Rij)
#         end
#
#     end
# end


################################################################################
# VORTEX FILAMENT ELEMENTS
################################################################################
"""
Computes the velocity induced by a bound vortex with beginning and end points
`[pa1, pa2, pa3]` and `[pb1, pb2, pb3]`, respectively, and vortex strength
`strength` on the targets `targets`. It adds the velocity at the i-th target to
out[i].
"""
function U_boundvortex( pa1::Number, pa2::Number, pa3::Number,
                        pb1::Number, pb2::Number, pb3::Number,
                        strength::Number,
                        targets::Arr2, out::Arr3;
                        dot_with=nothing,
                        cutoff=1e-14, offset=1e-8,
                       ) where{   T2, Arr2<:AbstractArray{T2,2},
                                  T3, Arr3<:AbstractArray{T3}}

    nt = size(targets, 2)                   # Number of targets
    no = dot_with!=nothing ? length(out) : size(out, 2) # Number of outputs

    if no!=nt
        error("Invalid `out` argument. Expected size $(nt), got $(no).")
    end

    # rij = pj - pi
    rij1 = pb1 - pa1
    rij2 = pb2 - pa2
    rij3 = pb3 - pa3

    # Iterate over targets
    # @simd for ti in 1:nt
    Threads.@threads for ti in 1:nt

        @inbounds begin
            # ri = x - pi
            ri1 = targets[1, ti] - pa1
            ri2 = targets[2, ti] - pa2
            ri3 = targets[3, ti] - pa3

            # rj = x - pj
            rj1 = targets[1, ti] - pb1
            rj2 = targets[2, ti] - pb2
            rj3 = targets[3, ti] - pb3
        end

        # ri × rj
        rixrj1 = ri2*rj3 - ri3*rj2
        rixrj2 = ri3*rj1 - ri1*rj3
        rixrj3 = ri1*rj2 - ri2*rj1

        # ‖ ri × rj ‖^2
        dotrixrj = rixrj1^2 + rixrj2^2 + rixrj3^2

        # rij ⋅ (hat{ri} - hat{rj}), add core offset to avoid singularity
        normri = sqrt(ri1^2 + ri2^2 + ri3^2) + offset
        normrj = sqrt(rj1^2 + rj2^2 + rj3^2) + offset
        rijdothat = rij1*(ri1/normri - rj1/normrj) + rij2*(ri2/normri - rj2/normrj) + rij3*(ri3/normri - rj3/normrj)

        if dotrixrj > cutoff^2 # This makes the self induced velocity zero

            # aux = strength * rijdothat / (4*π*(sqrt(dotrixrj) + offset)^2) # NOTE: This is correct expression, but it adds an extra sqrt
            aux = strength * rijdothat / (4*π*(dotrixrj + offset^2))

            # NOTE: Negative sign is not needed since we defined rij = rj - ri
            if dot_with!=nothing
                @inbounds out[ti] += aux * ( rixrj1*dot_with[1,ti] + rixrj2*dot_with[2,ti] + rixrj3*dot_with[3,ti] )
            else
                @inbounds out[1, ti] += aux * rixrj1
                @inbounds out[2, ti] += aux * rixrj2
                @inbounds out[3, ti] += aux * rixrj3
            end
        end
    end

end



"""
`U_vortexring(nodes::Matrix, panel::Array{Int}, strength::Real,
                targets::Matrix, out::Matrix; dot_with=nothing,
                closed_ring=true,
                cutoff=1e-14, offset=1e-8)`

Computes the velocity induced by a vortex ring panel of vertices
`nodes[:, panel]` and vortex strength `strength` on the targets `targets`. It
adds the velocity at the i-th target to out[i].
"""
function U_vortexring(nodes::Arr1, panel, strength, targets, out;
                        dot_with=nothing,
                        # closed_ring::Bool=true,
                        cutoff=1e-14, offset=1e-8,
                        omit_wake=false
                     ) where{T1, Arr1<:AbstractArray{T1,2}}

    nn = length(panel)                      # Number of nodes

    # Iterate over nodes
    for i in 1:nn

        # if closed_ring || i != nn

            @inbounds begin
                pi, pj = panel[i], panel[i%nn + 1]
                pa1, pa2, pa3 = nodes[1, pi], nodes[2, pi], nodes[3, pi]
                pb1, pb2, pb3 = nodes[1, pj], nodes[2, pj], nodes[3, pj]
            end

            U_boundvortex(pa1, pa2, pa3, pb1, pb2, pb3, strength, targets, out;
                            dot_with=dot_with, cutoff=cutoff, offset=offset)
        # end

    end
end



"""
Computes the velocity induced by a semi-infinite vortex starting at point
`[p1, p2, p3]` going out in the unitary direction `[d1, d2, d3]` with vortex
strength `strength` on the targets `targets`. It adds the velocity at the i-th
target to out[i].

A negative strength makes it such that the vortex comes from infinite and
ends at the point.
"""
function U_semiinfinite_vortex( p1::Number, p2::Number, p3::Number,
                                d1::Number, d2::Number, d3::Number,
                                strength::Number,
                                targets::Arr2, out::Arr3;
                                dot_with=nothing,
                                cutoff=1e-14, offset=1e-8,
                               ) where{   T2, Arr2<:AbstractArray{T2,2},
                                          T3, Arr3<:AbstractArray{T3}}


    nt = size(targets, 2)                               # Number of targets
    no = dot_with!=nothing ? length(out) : size(out, 2) # Number of outputs

    if no!=nt
        error("Invalid `out` argument. Expected size $(nt), got $(no).")
    elseif abs(d1^2 + d2^2 + d3^2 - 1) > 2*eps()
        error("Found non-unitary semi-infinite direction"*
                " norm([d1, d2, d3]) = norm($([d1,d2,d3])) = $(norm((d1,d2,d3)))")
    end


    # Iterate over targets
    # @simd for ti in 1:nt
    Threads.@threads for ti in 1:nt

        # Split vortex into bound and semi-infinite sections
        # p0 = p + [(x-p)⋅d]d
        @inbounds xmpdotd = (targets[1, ti] - p1)*d1 + (targets[2, ti] - p2)*d2 + (targets[3, ti] - p3)*d3
        p01 = p1 + xmpdotd*d1
        p02 = p2 + xmpdotd*d2
        p03 = p3 + xmpdotd*d3

        # ----------------- Bound Vortex ---------------------------------------
        if (p01-p1)^2 + (p02-p2)^2 + (p03-p3)^2 > offset^2 # Check that there is a bound section

            # rij = pj - pi
            rij1 = p01 - p1
            rij2 = p02 - p2
            rij3 = p03 - p3

            @inbounds begin
                # ri = x - pi
                ri1 = targets[1, ti] - p1
                ri2 = targets[2, ti] - p2
                ri3 = targets[3, ti] - p3

                # rj = x - pj
                rj1 = targets[1, ti] - p01
                rj2 = targets[2, ti] - p02
                rj3 = targets[3, ti] - p03
            end

            # ri × rj
            rixrj1 = ri2*rj3 - ri3*rj2
            rixrj2 = ri3*rj1 - ri1*rj3
            rixrj3 = ri1*rj2 - ri2*rj1

            # ‖ ri × rj ‖^2
            dotrixrj = rixrj1^2 + rixrj2^2 + rixrj3^2

            # rij ⋅ (hat{ri} - hat{rj}), add core offset to avoid singularity
            normri = sqrt(ri1^2 + ri2^2 + ri3^2) + offset
            normrj = sqrt(rj1^2 + rj2^2 + rj3^2) + offset
            rijdothat = rij1*(ri1/normri - rj1/normrj) + rij2*(ri2/normri - rj2/normrj) + rij3*(ri3/normri - rj3/normrj)

            if dotrixrj > cutoff^2 # This makes the self induced velocity zero

                # aux = strength * rijdothat / (4*π*(sqrt(dotrixrj) + offset)^2) # NOTE: This is correct expression, but it adds an extra sqrt
                aux = strength * rijdothat / (4*π*(dotrixrj + offset^2))

                # NOTE: Negative sign is not needed since we defined rij = rj - ri
                if dot_with!=nothing
                    @inbounds out[ti] += aux * ( rixrj1*dot_with[1,ti] + rixrj2*dot_with[2,ti] + rixrj3*dot_with[3,ti] )
                else
                    @inbounds out[1, ti] += aux * rixrj1
                    @inbounds out[2, ti] += aux * rixrj2
                    @inbounds out[3, ti] += aux * rixrj3
                end
            end

        end

        # ----------------- Semi-Infinite Vortex -------------------------------
        # h = ‖x - p0‖
        # @inbounds h = sqrt( (targets[1, ti] - p01)^2 + (targets[2, ti] - p02)^2 + (targets[3, ti] - p03)^2 )
        @inbounds hsqr = (targets[1, ti] - p01)^2 + (targets[2, ti] - p02)^2 + (targets[3, ti] - p03)^2

        # hhat = (x - p0) / ‖x - p0‖
        # h1 = (targets[1, ti] - p01)/h
        # h2 = (targets[2, ti] - p02)/h
        # h3 = (targets[3, ti] - p03)/h
        h1 = (targets[1, ti] - p01)
        h2 = (targets[2, ti] - p02)
        h3 = (targets[3, ti] - p03)

        # nhat = dhat × hhat
        n1 = d2*h3 - d3*h2
        n2 = d3*h1 - d1*h3
        n3 = d1*h2 - d2*h1

        # if h > cutoff # This makes the self induced velocity zero
        if hsqr > cutoff^2

            # aux = strength / (4*π*(h + offset))
            # aux = strength / (4*π*(sqrt(h2) + offset)^2) # NOTE: This is correct expression, but it adds an extra sqrt
            aux = strength / (4*π*(hsqr + offset^2))

            if dot_with!=nothing
                @inbounds out[ti] += aux * (n1*dot_with[1,ti] + n2*dot_with[2,ti] + n3*dot_with[3,ti])
            else
                @inbounds out[1, ti] += aux * n1
                @inbounds out[2, ti] += aux * n2
                @inbounds out[3, ti] += aux * n3
            end

        end

    end

end



"""
Computes the velocity induced by a semi-infinite horseshoe of strength
`strength` on the targets `targets`. The semi-infinite horseshoe comes from ∞
to `nodes[:, TE[1]]` and goes out to ∞ from `nodes[:, TE[2]]`, with a bound
vortex going from `nodes[:, TE[1]]` to `nodes[:, TE[2]]`. The direction is of
the semi-infinite sections is given by `[d1, d2, d3]`. It adds the velocity at
the i-th target to out[i].
"""
function U_semiinfinite_horseshoe(nodes::Arr1,
                                    TE, d1::Number, d2::Number, d3::Number,
                                    strength::Number,
                                    targets, out;
                                    dot_with=nothing,
                                    cutoff=1e-14, offset=1e-8,
                                    omit_wake=false
                                  ) where{T1, Arr1<:AbstractArray{T1,2}}

    if !omit_wake
        # Semi-infinite vortex coming in (from ∞ to pi)
        U_semiinfinite_vortex(  nodes[1, TE[1]], nodes[2, TE[1]], nodes[3, TE[1]],
                                d1, d2, d3, -strength,
                                targets, out;
                                dot_with=dot_with,
                                cutoff=cutoff, offset=offset
                             )
    end

    # Bound vortex (from pi and pj)
    U_boundvortex(  nodes[1, TE[1]], nodes[2, TE[1]], nodes[3, TE[1]],
                    nodes[1, TE[2]], nodes[2, TE[2]], nodes[3, TE[2]],
                    strength,
                    targets, out;
                    dot_with=dot_with,
                    cutoff=cutoff, offset=offset
                 )

    if !omit_wake
        # Semi-infinite vortex going to (from pj to ∞)
        U_semiinfinite_vortex(  nodes[1, TE[2]], nodes[2, TE[2]], nodes[3, TE[2]],
                                d1, d2, d3, strength,
                                targets, out;
                                dot_with=dot_with,
                                cutoff=cutoff, offset=offset
                             )
    end
end

"""
Computes the velocity induced by a semi-infinite horseshoe of strength
`strength` on the targets `targets`. The semi-infinite horseshoe comes from ∞
to pa and goes out to ∞ from pb, with a bound vortex going from pa to pb. The
semi-infinite directions of pa and pb are da and db, respectively. It adds the
velocity at the i-th target to out[i].

DEFINITIONS
* pa = `nodes[:, TE[1]]`
* pb = `nodes[:, TE[2]]`
* da = `[da1, da2, da3]`
* db = `[db1, db2, db3]`
"""
function U_semiinfinite_horseshoe(nodes::Arr1,
                                    TE,
                                    da1::Number, da2::Number, da3::Number,
                                    db1::Number, db2::Number, db3::Number,
                                    strength::Number,
                                    targets, out;
                                    dot_with=nothing,
                                    cutoff=1e-14, offset=1e-8,
                                    omit_wake=false
                                  ) where{T1, Arr1<:AbstractArray{T1,2}}

    if !omit_wake
        # Semi-infinite vortex coming in (from ∞ to pi)
        U_semiinfinite_vortex(  nodes[1, TE[1]], nodes[2, TE[1]], nodes[3, TE[1]],
                                da1, da2, da3, -strength,
                                targets, out;
                                dot_with=dot_with,
                                cutoff=cutoff, offset=offset
                             )
    end

    # Bound vortex (from pi and pj)
    U_boundvortex(  nodes[1, TE[1]], nodes[2, TE[1]], nodes[3, TE[1]],
                    nodes[1, TE[2]], nodes[2, TE[2]], nodes[3, TE[2]],
                    strength,
                    targets, out;
                    dot_with=dot_with,
                    cutoff=cutoff, offset=offset
                 )
    if !omit_wake
        # Semi-infinite vortex going to (from pj to ∞)
        U_semiinfinite_vortex(  nodes[1, TE[2]], nodes[2, TE[2]], nodes[3, TE[2]],
                                db1, db2, db3, strength,
                                targets, out;
                                dot_with=dot_with,
                                cutoff=cutoff, offset=offset
                             )
    end
end

################################################################################
# DOUBLET ELEMENTS
################################################################################
"""
Computes the velocity induced by a panel of vertices `nodes[:, panel]` and
constant strength doublet `strength` on the targets `targets`. It adds the
velocity at the i-th target to out[i].
"""
const U_constant_doublet = U_vortexring
# U_constant_doublet(args...; optargs...) = U_vortexring(args...; optargs...) #<--- This produces memory allocation for some reason


"""
Computes the potential induced by a panel of vertices `nodes[:, panel]` and
constant strength doublet `strength` on the targets `targets`. It adds the
potential at the i-th target to out[i].

Implementation of equations in Hess & Smith (1967), using ϕ_µ = −µ*(n⋅∇ϕ_σ/σ).
"""
function phi_constant_doublet(nodes::Arr1, panel,
                              strength::Number,
                              targets::Arr2, out::Arr3
                          ) where{T1, Arr1<:AbstractArray{T1,2},
                                  T2, Arr2<:AbstractArray{T2,2},
                                  T3, Arr3<:AbstractArray{T3}}

    nt = size(targets, 2)                   # Number of targets
    no = length(out)                        # Number of outputs
    nn = length(panel)                      # Number of nodes

    if no!=nt
        error("Invalid `out` argument. Expected size $(nt), got $(no).")
    end

    # Tangent, oblique, and normal vectors
    # t1, t2, t3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
    # o1, o2, o3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
    # n1, n2, n3 = gt._calc_n1(nodes, panel), gt._calc_n2(nodes, panel), gt._calc_n3(nodes, panel)

    # Here we flip the oblique and tangent vector since both Katz & Plotkin and
    # Hess & Smith define their panels opposite to the right-hand rule (while
    # GeometricTools defines the normal following the right-hand rule)
    t1, t2, t3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
    o1, o2, o3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
    n1, n2, n3 = -gt._calc_n1(nodes, panel), -gt._calc_n2(nodes, panel), -gt._calc_n3(nodes, panel)

    # Panel local coordinate system
    # NOTE: normal is pointing out of the body, which differs from Katz and Plotkin
    # @inbounds O = nodes[1]                         # Origin
    @inbounds Oi = panel[1]                         # Index of node that is the origin
    # xhat, yhat, zhat = t, o, n         # Unit vectors
    # Oaxis = hcat(xhat, yhat, zhat)'    # Transformation matrix

    # Iterate over targets
    for ti in 1:nt

        # Target position in panel coordinate system
        # X = Oaxis*(targets[:, ti]-O)
        @inbounds begin
            x = t1*(targets[1,ti]-nodes[1,Oi]) + t2*(targets[2,ti]-nodes[2,Oi]) + t3*(targets[3,ti]-nodes[3,Oi])
            y = o1*(targets[1,ti]-nodes[1,Oi]) + o2*(targets[2,ti]-nodes[2,Oi]) + o3*(targets[3,ti]-nodes[3,Oi])
            z = n1*(targets[1,ti]-nodes[1,Oi]) + n2*(targets[2,ti]-nodes[2,Oi]) + n3*(targets[3,ti]-nodes[3,Oi])
        end

        V3 = zero(T3)
        dtheta = 2*pi
        nR0::Int = 0

        @simd for i in 1:nn

            @inbounds begin
                pi, pj = panel[i], panel[i%nn + 1]

                # Convert nodes to panel coordinate system
                xi = t1*(nodes[1,pi]-nodes[1,Oi]) + t2*(nodes[2,pi]-nodes[2,Oi]) + t3*(nodes[3,pi]-nodes[3,Oi])
                yi = o1*(nodes[1,pi]-nodes[1,Oi]) + o2*(nodes[2,pi]-nodes[2,Oi]) + o3*(nodes[3,pi]-nodes[3,Oi])
                zi = n1*(nodes[1,pi]-nodes[1,Oi]) + n2*(nodes[2,pi]-nodes[2,Oi]) + n3*(nodes[3,pi]-nodes[3,Oi])
                xj = t1*(nodes[1,pj]-nodes[1,Oi]) + t2*(nodes[2,pj]-nodes[2,Oi]) + t3*(nodes[3,pj]-nodes[3,Oi])
                yj = o1*(nodes[1,pj]-nodes[1,Oi]) + o2*(nodes[2,pj]-nodes[2,Oi]) + o3*(nodes[3,pj]-nodes[3,Oi])
                zj = n1*(nodes[1,pj]-nodes[1,Oi]) + n2*(nodes[2,pj]-nodes[2,Oi]) + n3*(nodes[3,pj]-nodes[3,Oi])
            end

            dij = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
            ri = sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)
            rj = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2)

            Sij = (yj-yi)/dij
            Cij = (xj-xi)/dij

            siji = (xi-x)*Cij + (yi-y)*Sij
            sijj = (xj-x)*Cij + (yj-y)*Sij
            Rij = (x-xi)*Sij - (y-yi)*Cij

            Jij = atan( Rij*abs(z)*( ri*sijj - rj*siji ) , ri*rj*Rij^2 + z^2*sijj*siji )

            V3 -= Jij

            dtheta *= Rij>=0
            nR0 += Rij==0
            # nR0 += (abs(Rij) < offset)
        end

        V3 += dtheta
        V3 *= sign(z)        # Isn't this sign already accounted for in atan2?
        V3 *= !(nR0>1)       # Singularity fix of any z position aligned with node


        # NOTE: Katz and Plotkin's potential differs from Hess and Smith's by
        #       this factor
        V3 *= 1/(4*pi)

        # out[ti] -= strength*V3        # No need for - since it was already accounted for using the negative normal
        out[ti] += strength*V3

    end
end


# """
# Returns the potential induced by a panel of vertices `nodes[:, panel]` and constant
# strength doublet `strength` on the targets `targets`. It adds the potential at
# the i-th target to out[i].
#
# Implementation of equations in Katz and Plotkin Sec. 10.4.2.
# # NOTE: THOSE EQUATIONS ARE WRONG. THEY GIVE THE CORRECT POTENTIAL ONLY BELOW
# #       THE PANEL, AND NON-SENSE ABOVE IT.
# """
# function phi_constant_doublet(nodes::Arr1, panel,
#                               strength::Number,
#                               targets::Arr2, out::Arr3
#                              ) where{T1, Arr1<:AbstractArray{T1,2},
#                                      T2, Arr2<:AbstractArray{T2,2},
#                                      T3, Arr3<:AbstractArray{T3}}
#
#     nt = size(targets, 2)                   # Number of targets
#     no = length(out)                        # Number of outputs
#     nn = length(panel)                      # Number of nodes
#
#     if no!=nt
#         error("Invalid `out` argument. Expected size $(nt), got $(no).")
#     end
#
#     # Tangent, oblique, and normal vectors
#     # t1, t2, t3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
#     # o1, o2, o3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
#     # n1, n2, n3 = gt._calc_n1(nodes, panel), gt._calc_n2(nodes, panel), gt._calc_n3(nodes, panel)
#
#     # Here we flip the oblique and tangent vector since both Katz & Plotkin and
#     # Hess & Smith define their panels opposite to the right-hand rule (while
#     # GeometricTools defines the normal following the right-hand rule)
#     t1, t2, t3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
#     o1, o2, o3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
#     n1, n2, n3 = -gt._calc_n1(nodes, panel), -gt._calc_n2(nodes, panel), -gt._calc_n3(nodes, panel)
#
#     # Panel local coordinate system
#     # NOTE: normal is pointing out of the body, which differs from Katz and Plotkin
#     @inbounds Oi = panel[1]                         # Index of node that is the origin
#     # xhat, yhat, zhat = t, o, n         # Unit vectors
#     # Oaxis = hcat(xhat, yhat, zhat)'    # Transformation matrix
#
#     # Convert nodes to panel coordinate system
#     # Pnodes = [Oaxis*(node-O) for node in nodes]
#
#     # Iterate over nodes
#     for i in 1:nn
#
#         @inbounds begin
#             pi, pj = panel[i], panel[i%nn + 1]
#
#             # Convert nodes to panel coordinate system
#             xi = t1*(nodes[1,pi]-nodes[1,Oi]) + t2*(nodes[2,pi]-nodes[2,Oi]) + t3*(nodes[3,pi]-nodes[3,Oi])
#             yi = o1*(nodes[1,pi]-nodes[1,Oi]) + o2*(nodes[2,pi]-nodes[2,Oi]) + o3*(nodes[3,pi]-nodes[3,Oi])
#             zi = n1*(nodes[1,pi]-nodes[1,Oi]) + n2*(nodes[2,pi]-nodes[2,Oi]) + n3*(nodes[3,pi]-nodes[3,Oi])
#             xj = t1*(nodes[1,pj]-nodes[1,Oi]) + t2*(nodes[2,pj]-nodes[2,Oi]) + t3*(nodes[3,pj]-nodes[3,Oi])
#             yj = o1*(nodes[1,pj]-nodes[1,Oi]) + o2*(nodes[2,pj]-nodes[2,Oi]) + o3*(nodes[3,pj]-nodes[3,Oi])
#             zj = n1*(nodes[1,pj]-nodes[1,Oi]) + n2*(nodes[2,pj]-nodes[2,Oi]) + n3*(nodes[3,pj]-nodes[3,Oi])
#         end
#
#
#         # Iterate over targets
#         @simd for ti in 1:nt
#
#             # Target position in panel coordinate system
#             # X = Oaxis*(targets[:, ti]-O)
#             @inbounds begin
#                 x = t1*(targets[1,ti]-nodes[1,Oi]) + t2*(targets[2,ti]-nodes[2,Oi]) + t3*(targets[3,ti]-nodes[3,Oi])
#                 y = o1*(targets[1,ti]-nodes[1,Oi]) + o2*(targets[2,ti]-nodes[2,Oi]) + o3*(targets[3,ti]-nodes[3,Oi])
#                 z = n1*(targets[1,ti]-nodes[1,Oi]) + n2*(targets[2,ti]-nodes[2,Oi]) + n3*(targets[3,ti]-nodes[3,Oi])
#             end
#
#             mij = (yj - yi)/(xj - xi)
#             ri = sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)
#             rj = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2)
#             ei = (x - xi)^2 + (z-zi)^2
#             ej = (x - xj)^2 + (z-zj)^2
#             hi = (x - xi)*(y - yi)
#             hj = (x - xj)*(y - yj)
#
#             @inbounds out[ti] += strength/(4*π) * (atan(mij*ei-hi, z*ri) - atan(mij*ej-hj, z*rj))
#         end
#
#     end
# end



################################################################################
# SEMI-INFINITE DOUBLET ELEMENTS
################################################################################
"""
Computes the velocity induced by a semi-infinite doublet panel of strength
`strength` on the targets `targets`. The semi-infinite panel starts at the
trailing edge between `nodes[:, TE[1]]` and `nodes[:, TE[2]]` and extends
infinitely in the direction `[d1, d2, d3]`. It adds the potential at
the i-th target to out[i].
"""
const U_semiinfinite_doublet = U_semiinfinite_horseshoe

"""
Computes the potential induced by a semi-infinite doublet panel of strength
`strength` on the targets `targets`. The semi-infinite panel starts at the
trailing edge between `nodes[:, TE[1]]` and `nodes[:, TE[2]]` and extends
infinitely in the direction `[d1, d2, d3]`. It adds the potential at
the i-th target to out[i].

Implementation of equations in Moran, J. (1984), Appendix F, p. 445.
"""
function phi_semiinfinite_doublet(nodes::Arr1,
                                    TE, d1::Number, d2::Number, d3::Number,
                                    strength::Number,
                                    targets::Arr2, out::Arr3
                                  ) where{T1, Arr1<:AbstractArray{T1,2},
                                          T2, Arr2<:AbstractArray{T2,2},
                                          T3, Arr3<:AbstractArray{T3}}

    nt = size(targets, 2)                   # Number of targets
    no = length(out)                        # Number of outputs

    if no!=nt
        error("Invalid `out` argument. Expected size $(nt), got $(no).")
    elseif abs(d1^2 + d2^2 + d3^2 - 1) > 2*eps()
        error("Found non-unitary semi-infinite direction"*
                " norm([d1, d2, d3]) = norm($([d1,d2,d3])) = $(norm((d1,d2,d3)))")
    end

    # Split panel into a bound panel and a semi-infinite panel
    @inbounds begin
        p_i1, p_i2, p_i3 = nodes[1,TE[1]], nodes[2,TE[1]], nodes[3,TE[1]]
        p_j1, p_j2, p_j3 = nodes[1,TE[2]], nodes[2,TE[2]], nodes[3,TE[2]]

        # p_a = p_i + [(p_j-p_i)⋅d]d
        pijdotd = (p_j1 - p_i1)*d1 + (p_j2 - p_i2)*d2 + (p_j3 - p_i3)*d3
        pa1 = p_i1 + pijdotd*d1
        pa2 = p_i2 + pijdotd*d2
        pa3 = p_i3 + pijdotd*d3

        # Panel local coordinate system
        x1, x2, x3 = d1, d2, d3
        nrmpbpa = sqrt( (p_j1-pa1)^2 + (p_j2-pa2)^2 + (p_j3-pa3)^2  )
        y1, y2, y3 = (p_j1-pa1)/nrmpbpa, (p_j2-pa2)/nrmpbpa, (p_j3-pa3)/nrmpbpa
        z1, z2, z3 = x2*y3-x3*y2, x3*y1-x1*y3, x1*y2-x2*y1

        O1, O2, O3 = p_i1, p_i2, p_i3                   # Origin
        # Oaxis = hcat(xhat, yhat, zhat)'               # Transformation matrix
    end


    # Iterate over targets
    for ti in 1:nt

        # Target position in panel coordinate system
        # X = Oaxis*(targets[:, ti]-O)
        @inbounds begin
            x = x1*(targets[1,ti]-O1) + x2*(targets[2,ti]-O2) + x3*(targets[3,ti]-O3)
            y = y1*(targets[1,ti]-O1) + y2*(targets[2,ti]-O2) + y3*(targets[3,ti]-O3)
            z = z1*(targets[1,ti]-O1) + z2*(targets[2,ti]-O2) + z3*(targets[3,ti]-O3)
        end

        # ------------ Potential of bound panel
        if abs(pijdotd) > 2*eps()               # <--- Avoids the case that there is no bound panel

            V3 = zero(T3)
            dtheta = 2*pi
            nR0::Int = 0

            @simd for i in 1:3 # iterate over p_a, p_i, and p_j

                if i==1         # Case (pi, pj) = (p_a, p_i)
                    pi1, pi2, pi3 = pa1, pa2, pa3
                    pj1, pj2, pj3 = p_i1, p_i2, p_i3
                elseif i==2     # Case (pi, pj) = (p_i, p_j)
                    pi1, pi2, pi3 = p_i1, p_i2, p_i3
                    pj1, pj2, pj3 = p_j1, p_j2, p_j3
                else            # Case (pi, pj) = (p_j, p_a)
                    pi1, pi2, pi3 = p_j1, p_j2, p_j3
                    pj1, pj2, pj3 = pa1, pa2, pa3
                end

                # Convert nodes to panel coordinate system
                xi = x1*(pi1-O1) + x2*(pi2-O2) + x3*(pi3-O3)
                yi = y1*(pi1-O1) + y2*(pi2-O2) + y3*(pi3-O3)
                zi = z1*(pi1-O1) + z2*(pi2-O2) + z3*(pi3-O3)
                xj = x1*(pj1-O1) + x2*(pj2-O2) + x3*(pj3-O3)
                yj = y1*(pj1-O1) + y2*(pj2-O2) + y3*(pj3-O3)
                zj = z1*(pj1-O1) + z2*(pj2-O2) + z3*(pj3-O3)

                dij = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
                ri = sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)
                rj = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2)

                Sij = (yj-yi)/dij
                Cij = (xj-xi)/dij

                siji = (xi-x)*Cij + (yi-y)*Sij
                sijj = (xj-x)*Cij + (yj-y)*Sij
                Rij = (x-xi)*Sij - (y-yi)*Cij

                Jij = atan( Rij*abs(z)*( ri*sijj - rj*siji ) , ri*rj*Rij^2 + z^2*sijj*siji )

                V3 -= Jij

                dtheta *= Rij>=0
                nR0 += Rij==0
                # nR0 += (abs(Rij) < offset)
            end

            V3 += dtheta
            V3 *= sign(z)        # Isn't this sign already accounted for in atan2?
            V3 *= !(nR0>1)       # Singularity fix of any z position aligned with node


            # NOTE: Katz and Plotkin's potential differs from Hess and Smith's by
            #       this factor
            V3 *= 1/(4*pi)

            # out[ti] -= strength*V3        # No need for - since it was already accounted for using the negative normal
            out[ti] += strength*V3
        end


        # ------------ Potential of semi-infinite panel
        val = zero(T3)

        # Convert nodes to panel coordinate system
        xa = x1*(pa1-O1) + x2*(pa2-O2) + x3*(pa3-O3)
        ya = y1*(pa1-O1) + y2*(pa2-O2) + y3*(pa3-O3)
        # za = z1*(pa1-O1) + z2*(pa2-O2) + z3*(pa3-O3)
        xb = x1*(p_j1-O1) + x2*(p_j2-O2) + x3*(p_j3-O3)
        yb = y1*(p_j1-O1) + y2*(p_j2-O2) + y3*(p_j3-O3)
        # zb = z1*(p_j1-O1) + z2*(p_j2-O2) + z3*(p_j3-O3)

        # TODO: What is the domain of evaluation of this atan function in the theory?
        val += atan((yb-y)/z) + atan( (yb-y)*(x-xa) / (z*sqrt((x-xa)^2 + (yb-y)^2 + z^2)) )
        val -= atan((ya-y)/z) + atan( (ya-y)*(x-xa) / (z*sqrt((x-xa)^2 + (ya-y)^2 + z^2)) )
        # val += atan(yb-y, z) + atan( (yb-y)*(x-xa), z*sqrt((x-xa)^2 + (yb-y)^2 + z^2) )
        # val -= atan(ya-y, z) + atan( (ya-y)*(x-xa), z*sqrt((x-xa)^2 + (ya-y)^2 + z^2) )

        # out[ti] -= strength/(4*pi) * val
        out[ti] += strength/(4*pi) * val # NOTE: For some reason I ended up having to flip this sign to match the potential of the finite doublet panel
    end
end

"""
Computes the potential induced by a semi-infinite horseshoe of strength
`strength` on the targets `targets`. The semi-infinite horseshoe comes from ∞
to pa and goes out to ∞ from pb, with a bound vortex going from pa to pb. The
semi-infinite directions of pa and pb are da and db, respectively. It adds the
velocity at the i-th target to out[i].

DEFINITIONS
* pa = `nodes[:, TE[1]]`
* pb = `nodes[:, TE[2]]`
* da = `[da1, da2, da3]`
* db = `[db1, db2, db3]`

Implementation of equations in Moran, J. (1984), Appendix F, p. 445, splitting
the non-planar semi-infinite panel into two planar sections.
"""
function phi_semiinfinite_doublet(nodes::Arr1,
                                    TE,
                                    da1::Number, da2::Number, da3::Number,
                                    db1::Number, db2::Number, db3::Number,
                                    strength::Number,
                                    targets::Arr2, out::Arr3;
                                    infinitelength::Real=1e12
                                  ) where{T1, Arr1<:AbstractArray{T1,2},
                                          T2, Arr2<:AbstractArray{T2,2},
                                          T3, Arr3<:AbstractArray{T3}}

    nt = size(targets, 2)                   # Number of targets
    no = length(out)                        # Number of outputs

    if no!=nt
        error("Invalid `out` argument. Expected size $(nt), got $(no).")
    elseif abs(da1^2 + da2^2 + da3^2 - 1) > 2*eps()
        error("Found non-unitary semi-infinite direction"*
              " norm([da1, da2, da3]) = norm($([da1,da2,da3])) = $(norm((da1,da2,da3)))")
    elseif abs(db1^2 + db2^2 + db3^2 - 1) > 2*eps()
        error("Found non-unitary semi-infinite direction"*
              " norm([db1, db2, db3]) = norm($([db1,db2,db3])) = $(norm((db1,db2,db3)))")
    elseif da1*db1 + da2*db2 + da3*db3 <= 0
        error("Invalid semi-infinite directions:"*
              " da=$([da1,da2,da3]) is opposite or orthogonal to db=$([db1,db2,db3])")
    end

    # Calculate potential as if panel was planar using da
    phi_semiinfinite_doublet(nodes, TE, da1, da2, da3, strength, targets, out)

    # Calculate triangular semi-infinite section between da and db going from pb
    if abs(da1-db1) > 2*eps() || abs(da2-db2) > 2*eps() || abs(da3-db3) > 2*eps()

        # Numerical approximation of semi-infinite panel with a large panel
        @inbounds begin
            pb1, pb2, pb3 = nodes[1, TE[2]], nodes[2, TE[2]], nodes[3, TE[2]]
        end
        pda1, pda2, pda3 = pb1+infinitelength*da1, pb2+infinitelength*da2, pb3+infinitelength*da3
        pdb1, pdb2, pdb3 = pb1+infinitelength*db1, pb2+infinitelength*db2, pb3+infinitelength*db3

        # Panel local coordinate system
        x1, x2, x3 = da1, da2, da3
        z1, z2, z3 = da2*db3-da3*db2, da3*db1-da1*db3, da1*db2-da2*db1
        nrmz = sqrt(z1^2 + z2^2 + z3^2)
        z1 /= nrmz
        z2 /= nrmz
        z3 /= nrmz
        y1, y2, y3 = z2*x3-z3*x2, z3*x1-z1*x3, z1*x2-z2*x1

        O1, O2, O3 = pb1, pb2, pb3                      # Origin
        # Oaxis = hcat(xhat, yhat, zhat)'               # Transformation matrix

        # Iterate over targets
        for ti in 1:nt

            # Target position in panel coordinate system
            # X = Oaxis*(targets[:, ti]-O)
            @inbounds begin
                x = x1*(targets[1,ti]-O1) + x2*(targets[2,ti]-O2) + x3*(targets[3,ti]-O3)
                y = y1*(targets[1,ti]-O1) + y2*(targets[2,ti]-O2) + y3*(targets[3,ti]-O3)
                z = z1*(targets[1,ti]-O1) + z2*(targets[2,ti]-O2) + z3*(targets[3,ti]-O3)
            end

            V3 = zero(T3)
            dtheta = 2*pi
            nR0::Int = 0

            @simd for i in 1:3 # iterate over pda, pdb, and pb

                # if i==1         # Case (pi, pj) = (pda, pdb)
                #     pi1, pi2, pi3 = pda1, pda2, pda3
                #     pj1, pj2, pj3 = pdb1, pdb2, pdb3
                # elseif i==2     # Case (pi, pj) = (pdb, pb)
                #     pi1, pi2, pi3 = pdb1, pdb2, pdb3
                #     pj1, pj2, pj3 = pb1, pb2, pb3
                # else            # Case (pi, pj) = (pb, pda)
                #     pi1, pi2, pi3 = pb1, pb2, pb3
                #     pj1, pj2, pj3 = pda1, pda2, pda3
                # end

                if i==1         # Case (pi, pj) = (pda, pb)
                    pi1, pi2, pi3 = pda1, pda2, pda3
                    pj1, pj2, pj3 = pb1, pb2, pb3
                elseif i==2     # Case (pi, pj) = (pb, pdb)
                    pi1, pi2, pi3 = pb1, pb2, pb3
                    pj1, pj2, pj3 = pdb1, pdb2, pdb3
                else            # Case (pi, pj) = (pdb, pda)
                    pi1, pi2, pi3 = pdb1, pdb2, pdb3
                    pj1, pj2, pj3 = pda1, pda2, pda3
                end

                # Convert nodes to panel coordinate system
                xi = x1*(pi1-O1) + x2*(pi2-O2) + x3*(pi3-O3)
                yi = y1*(pi1-O1) + y2*(pi2-O2) + y3*(pi3-O3)
                zi = z1*(pi1-O1) + z2*(pi2-O2) + z3*(pi3-O3)
                xj = x1*(pj1-O1) + x2*(pj2-O2) + x3*(pj3-O3)
                yj = y1*(pj1-O1) + y2*(pj2-O2) + y3*(pj3-O3)
                zj = z1*(pj1-O1) + z2*(pj2-O2) + z3*(pj3-O3)

                dij = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
                ri = sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)
                rj = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2)

                Sij = (yj-yi)/dij
                Cij = (xj-xi)/dij

                siji = (xi-x)*Cij + (yi-y)*Sij
                sijj = (xj-x)*Cij + (yj-y)*Sij
                Rij = (x-xi)*Sij - (y-yi)*Cij

                Jij = atan( Rij*abs(z)*( ri*sijj - rj*siji ) , ri*rj*Rij^2 + z^2*sijj*siji )

                V3 -= Jij

                dtheta *= Rij>=0
                nR0 += Rij==0
                # nR0 += (abs(Rij) < offset)
            end

            V3 += dtheta
            V3 *= sign(z)        # Isn't this sign already accounted for in atan2?
            V3 *= !(nR0>1)       # Singularity fix of any z position aligned with node


            # NOTE: Katz and Plotkin's potential differs from Hess and Smith's by
            #       this factor
            V3 *= 1/(4*pi)

            # out[ti] -= strength*V3        # No need for - since it was already accounted for using the negative normal
            out[ti] += strength*V3

        end

    end
end





################################################################################
# VORTEX SHEET ELEMENTS
################################################################################
"""
Computes the velocity induced by a panel of vertices `nodes[:, panel]` and
constant vortex sheet (with strength components `gammat` and `gammao`) on the
targets `targets`. It adds the velocity at the i-th target to out[i].

Implementation of equations in Pate's 2017 doctoral dissertation,
*A Surface Vorticity Method for Wake-Body Interactions*, Appendix A.
"""
function U_constant_vortexsheet(nodes::Arr1, panel,
                                gammat::Number, gammao::Number,
                                targets::Arr2, out::Arr3;
                                dot_with=nothing,
                                cutoff=1e-14, offset=1e-8
                              ) where{T1, Arr1<:AbstractArray{T1,2},
                                      T2, Arr2<:AbstractArray{T2,2},
                                      T3, Arr3<:AbstractArray{T3}}

    nt = size(targets, 2)                   # Number of targets
    no = dot_with!=nothing ? length(out) : size(out, 2) # Number of outputs
    nn = length(panel)                      # Number of nodes

    if no!=nt
        error("Invalid `out` argument. Expected size $(nt), got $(no).")
    end

    # @warn("Sheet thickness has been hardcoded!")

    #=
        TODO
        * [ ] Implement efficient H00, H10, and H01 formulas
    =#

    # Tangent, oblique, and normal vectors
    t1, t2, t3 = gt._calc_t1(nodes, panel), gt._calc_t2(nodes, panel), gt._calc_t3(nodes, panel)
    o1, o2, o3 = gt._calc_o1(nodes, panel), gt._calc_o2(nodes, panel), gt._calc_o3(nodes, panel)
    n1, n2, n3 = gt._calc_n1(nodes, panel), gt._calc_n2(nodes, panel), gt._calc_n3(nodes, panel)

    @inbounds p1i = panel[1]                    # Index of first vertex

    # V1, V2, V3 = zero(T3), zero(T3), zero(T3)

    # Iterate over targets
    for ti in 1:nt

        # V1 *= 0
        # V2 *= 0
        # V3 *= 0

        V1, V2, V3 = zero(T3), zero(T3), zero(T3)

        # Projection of target onto plane of the panel
        @inbounds begin
            z = n1*(targets[1,ti]-nodes[1,p1i]) + n2*(targets[2,ti]-nodes[2,p1i]) + n3*(targets[3,ti]-nodes[3,p1i])
            px1 = targets[1,ti] - z*n1
            px2 = targets[2,ti] - z*n2
            px3 = targets[3,ti] - z*n3
        end

        # Iterate over triangles of integration
        # @simd
        for Si in 1:nn


            @inbounds begin
                # Indices of first and second vertices of this triangle
                pi = panel[Si]
                pip1 = Si == nn ? panel[1] : panel[Si+1]

                # Vertices in the coordinate system of integration
                q11 = nodes[1, pi] - px1
                q12 = nodes[2, pi] - px2
                q13 = nodes[3, pi] - px3
                q21 = nodes[1, pip1] - px1
                q22 = nodes[2, pip1] - px2
                q23 = nodes[3, pip1] - px3
            end

            sqrtq2mq1 = sqrt((q21 - q11)^2 + (q22 - q12)^2 + (q23 - q13)^2)

            # Axes of coordinate system for integration
            c21 = (q21 - q11) / sqrtq2mq1
            c22 = (q22 - q12) / sqrtq2mq1
            c23 = (q23 - q13) / sqrtq2mq1
            c11 = c22*n3 - c23*n2
            c12 = c23*n1 - c21*n3
            c13 = c21*n2 - c22*n1

            # Decompose strength with coordinate system of integration
            a00 = gammat*(t1*c11 + t2*c12 + t3*c13) + gammao*(o1*c11 + o2*c12 + o3*c13)
            b00 = gammat*(t1*c21 + t2*c22 + t3*c23) + gammao*(o1*c21 + o2*c22 + o3*c23)

            # Limits of integration
            a = q11*c11 + q12*c12 + q13*c13
            l1 = q11*c21 + q12*c22 + q13*c23
            l2 = q21*c21 + q22*c22 + q23*c23

            # Regularized height
            h = sqrt(z^2 + offset^2)
            # h = sqrt(z^2 + (1.0*1e-2)^2)
            # h = sqrt(z^2 + (2.5*1e-3)^2)
            # h = sqrt(z^2 + (1.0*1e-8)^2)

            sqrtl1a = sqrt(l1^2 + a^2)
            sqrtl2a = sqrt(l2^2 + a^2)
            sqrtl1ah = sqrt(l1^2 + a^2 + h^2)
            sqrtl2ah = sqrt(l2^2 + a^2 + h^2)
            logh = log(h)

            # Integral terms
            # NOTE: What is the codomain of atan? Should this be atan or atan2?
            H00 =  1/h * atan(a*l2, a^2 + h^2 + h*sqrtl2ah)
            H00 -= 1/h * atan(a*l1, a^2 + h^2 + h*sqrtl1ah)
            H10 =  l2/sqrtl2a * log(sqrtl2ah + sqrtl2a) - log(l2 + sqrtl2ah) - l2*logh/sqrtl2a
            H10 -= l1/sqrtl1a * log(sqrtl1ah + sqrtl1a) - log(l1 + sqrtl1ah) - l1*logh/sqrtl1a
            H01 =  a/sqrtl2a * (logh - log(sqrtl2ah + sqrtl2a))
            H01 -= a/sqrtl1a * (logh - log(sqrtl1ah + sqrtl1a))

            # Avoid zero divided by zero when the projection lays on a vertex
            if abs(a)<=cutoff && (abs(l2)<=cutoff || abs(l1)<=cutoff)
                nothing
            else
                V1 += z*b00*H00*c11 - z*a00*H00*c21 + (b00*H10 - a00*H01)*n1
                V2 += z*b00*H00*c12 - z*a00*H00*c22 + (b00*H10 - a00*H01)*n2
                V3 += z*b00*H00*c13 - z*a00*H00*c23 + (b00*H10 - a00*H01)*n3
            end
        end

        V1 /= 4*pi
        V2 /= 4*pi
        V3 /= 4*pi

        if dot_with!=nothing
            @inbounds out[ti] += V1*dot_with[1,ti] + V2*dot_with[2,ti] + V3*dot_with[3,ti]
        else
            @inbounds out[1, ti] += V1
            @inbounds out[2, ti] += V2
            @inbounds out[3, ti] += V3
        end

    end

end
