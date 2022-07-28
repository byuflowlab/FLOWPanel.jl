#=##############################################################################
# DESCRIPTION
    Definition of panel elements

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jul 2018
  * License   : MIT License
=###############################################################################



const SMOOTH = 1e-2               # Smoothing radius of source panel
const SMOOTH2 = 2e-3              # Smoothing radius for vortex ring
# const SMOOTH2 = 1e-8              # Smoothing radius for vortex ring
const SMOOTH3 = SMOOTH2           # Smoothing radius for semi-infinite vortex
  # SMOOTH3 = 1e-16
const SMOOTH4 = 1e-5              # Smoothing radius for doublet panel
const SMOOTH5 = 1e-8              # Smoothing radius for vortex ring
# SMOOTH5 = 1e-10


################################################################################
# SOURCE ELEMENTS
################################################################################



"""
Returns the velocity induced by a panel of vertices `nodes` and constant
strength source `strength` on the targets `targets`. It adds the velocity at the
i-th target to out[i].

Implementation of equations in Katz and Plotkin Sec. 10.4.1.
"""
function Vconstant_source(nodes::Array{Arr1,1}, strength::RType,
                          targets::Array{Array{T2,1},1},
                          out::Array{Arr3,1};
                          # out;
                          dot_with=nothing
                          ) where{T1, Arr1<:AbstractArray{T1}, T2<:RType,
                                    T3, Arr3<:AbstractArray{T3}}
    if size(out)!=size(targets)
        error("Invalid `out` argument."*
                " Expected size $(size(targets)), got $(size(out)).")
    end

    nn = size(nodes, 1)                      # Number of nodes

    # Tangent, oblique, and normal vectors
    t1, t2, t3 = gt._calc_t1(nodes), gt._calc_t2(nodes), gt._calc_t3(nodes)
    o1, o2, o3 = gt._calc_o1(nodes), gt._calc_o2(nodes), gt._calc_o3(nodes)
    n1, n2, n3 = gt._calc_n1(nodes), gt._calc_n2(nodes), gt._calc_n3(nodes)

    # Panel local coordinate system
    # NOTE: normal is pointing out of the body, which differs from Katz and Plotkin
    O = nodes[1]                         # Origin
    # xhat, yhat, zhat = t, o, n         # Unit vectors
    # Oaxis = hcat(xhat, yhat, zhat)'    # Transformation matrix

    # Iterates over targets
    for ti in 1:size(targets, 1)

        # Target position in panel coordinate system
        # X = Oaxis*(targets[ti]-O)
        x = t1*(targets[ti][1]-O[1]) + t2*(targets[ti][2]-O[2]) + t3*(targets[ti][3]-O[3])
        y = o1*(targets[ti][1]-O[1]) + o2*(targets[ti][2]-O[2]) + o3*(targets[ti][3]-O[3])
        z = n1*(targets[ti][1]-O[1]) + n2*(targets[ti][2]-O[2]) + n3*(targets[ti][3]-O[3])

        V1, V2, V3 = zero(T3), zero(T3), zero(T3)
        dtheta = 2*pi

        nR0 = 0

        for i in 1:nn
            pi, pj = nodes[i], nodes[i%nn + 1]

            # Converts nodes to panel coordinate system
            xi = t1*(pi[1]-O[1]) + t2*(pi[2]-O[2]) + t3*(pi[3]-O[3])
            yi = o1*(pi[1]-O[1]) + o2*(pi[2]-O[2]) + o3*(pi[3]-O[3])
            zi = n1*(pi[1]-O[1]) + n2*(pi[2]-O[2]) + n3*(pi[3]-O[3])
            xj = t1*(pj[1]-O[1]) + t2*(pj[2]-O[2]) + t3*(pj[3]-O[3])
            yj = o1*(pj[1]-O[1]) + o2*(pj[2]-O[2]) + o3*(pj[3]-O[3])
            zj = n1*(pj[1]-O[1]) + n2*(pj[2]-O[2]) + n3*(pj[3]-O[3])

            dij = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
            ri = sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)
            rj = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2)

            Qij = log( (ri+rj+dij)/(ri+rj-dij + SMOOTH) )

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
            out[ti] += strength*(V1*t1 + V2*o1 + V3*n1)*dot_with[ti][1]
            out[ti] += strength*(V1*t2 + V2*o2 + V3*n2)*dot_with[ti][2]
            out[ti] += strength*(V1*t3 + V2*o3 + V3*n3)*dot_with[ti][3]
        else
            out[ti][1] += strength*(V1*t1 + V2*o1 + V3*n1)
            out[ti][2] += strength*(V1*t2 + V2*o2 + V3*n2)
            out[ti][3] += strength*(V1*t3 + V2*o3 + V3*n3)
        end

    end
end


"""
Returns the potential induced by a panel of vertices `nodes` and constant
strength source `strength` on the targets `targets`. It adds the potential at
the i-th target to out[i].

Implementation of equations in Katz and Plotkin Sec. 10.4.1.
"""
function phi_constant_source(nodes::Array{Arr1,1}, strength::RType,
                              targets::Array{Array{T2,1},1},
                              # out::Array{Array{T,1},1}
                              out
                             ) where{T1, Arr1<:AbstractArray{T1}, T2<:RType}
    if size(out)!=size(targets)
        error("Invalid `out` argument."*
              " Expected size $(size(targets)), got $(size(out)).")
    end

    nn = size(nodes, 1)                      # Number of nodes

    # Tangent, oblique, and normal vectors
    t1, t2, t3 = gt._calc_t1(nodes), gt._calc_t2(nodes), gt._calc_t3(nodes)
    o1, o2, o3 = gt._calc_o1(nodes), gt._calc_o2(nodes), gt._calc_o3(nodes)
    n1, n2, n3 = gt._calc_n1(nodes), gt._calc_n2(nodes), gt._calc_n3(nodes)

    # Panel local coordinate system
    # NOTE: normal is pointing out of the body, which differs from Katz and Plotkin
    O = nodes[1]                         # Origin
    # xhat, yhat, zhat = t, o, n         # Unit vectors
    # Oaxis = hcat(xhat, yhat, zhat)'    # Transformation matrix

    # Converts nodes to panel coordinate system
    # Pnodes = [Oaxis*(node-O) for node in nodes]

    # Iterates over targets
    for ti in 1:size(targets, 1)

        phi = 0

        # Target position in panel coordinate system
        # X = Oaxis*(targets[ti]-O)
        x = t1*(targets[ti][1]-O[1]) + t2*(targets[ti][2]-O[2]) + t3*(targets[ti][3]-O[3])
        y = o1*(targets[ti][1]-O[1]) + o2*(targets[ti][2]-O[2]) + o3*(targets[ti][3]-O[3])
        z = n1*(targets[ti][1]-O[1]) + n2*(targets[ti][2]-O[2]) + n3*(targets[ti][3]-O[3])


        for i in 1:nn
            pi, pj = nodes[i], nodes[i%nn + 1]

            # Converts nodes to panel coordinate system
            xi = t1*(pi[1]-O[1]) + t2*(pi[2]-O[2]) + t3*(pi[3]-O[3])
            yi = o1*(pi[1]-O[1]) + o2*(pi[2]-O[2]) + o3*(pi[3]-O[3])
            zi = n1*(pi[1]-O[1]) + n2*(pi[2]-O[2]) + n3*(pi[3]-O[3])
            xj = t1*(pj[1]-O[1]) + t2*(pj[2]-O[2]) + t3*(pj[3]-O[3])
            yj = o1*(pj[1]-O[1]) + o2*(pj[2]-O[2]) + o3*(pj[3]-O[3])
            zj = n1*(pj[1]-O[1]) + n2*(pj[2]-O[2]) + n3*(pj[3]-O[3])

            dij = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
            mij = (yj - yi)/(xj - xi)
            ri = sqrt((x-xi)^2 + (y-yi)^2 + (z-zi)^2)
            rj = sqrt((x-xj)^2 + (y-yj)^2 + (z-zj)^2)
            ei = (x - xi)^2 + (z-zi)^2
            ej = (x - xj)^2 + (z-zj)^2
            hi = (x - xi)*(y - yi)
            hj = (x - xj)*(y - yj)

            Pij = (x - xi)*(yj - yi) - (y - yi)*(xj - xi)
            Qij = log( (ri+rj+dij)/(ri+rj-dij + SMOOTH) )
            Rij = atan(mij*ei-hi, z*ri) - atan(mij*ej-hj, z*rj)

            phi += Pij/dij * Qij - abs(z)*Rij
        end

        phi *= strength

        # NOTE: Katz and Plotkin's potential differs from Hess and Smith's by
        #       this factor
        phi *= -1/(4*pi)

        out[ti] += phi
    end
end


################################################################################
# DOUBLET ELEMENTS
################################################################################

"""
Returns the velocity induced by a panel of vertices `nodes` and constant
strength doublet `strength` on the targets `targets`. It adds the velocity at
the i-th target to out[i].
"""
function Vconstant_doublet(nodes::Array{Array{T1,1},1}, strength::RType,
                          targets::Array{Array{T2,1},1},
                          out;
                          dot_with::Union{Array{Array{T3,1},1}, Nothing}=nothing
                          ) where{T1<:RType, T2<:RType, T3<:RType}
  if size(out)!=size(targets)
    error("Invalid `out` argument."*
          " Expected size $(size(targets)), got $(size(out)).")
  end

  nn = size(nodes, 1)                      # Number of nodes

  # Tangent, oblique, and normal vectors
  t, o, n = gt._calc_unitvectors(nodes)

  # Coordinate system defined by Hess & Smith
  unitxi, uniteta, unitz = o, t, -n        # Unit vectors
  O = nodes[1]                             # Origin
  Oaxis = hcat(unitxi, uniteta, unitz)'    # Transformation matrix

  # Converts nodes to H&S coordinate system
  HSnodes = [Oaxis*(node-O) for node in nodes]

  # Iterates over targets
  for ti in 1:size(targets, 1)
    HSX = Oaxis*(targets[ti]-O)
    V = Array{T2}(undef, 3)

    for i in 1:nn
      xi, xj = HSnodes[i], HSnodes[i%nn + 1]

      ri = norm(HSX-xi)
      rj = norm(HSX-xj)

      if (ri>SMOOTH4 && rj>SMOOTH4 &&
                        abs( ri*rj - dot(HSX-xi, HSX-xj) )>SMOOTH4*SMOOTH4)

        aux = (ri+rj)/( ri*rj*( ri*rj - dot(HSX-xi, HSX-xj) ) )
        V[1] -= HSX[3]*(xj[2]-xi[2])*aux
        V[2] += HSX[3]*(xj[1]-xi[1])*aux
        V[3] += ( (HSX[1]-xj[1])*(HSX[2]-xi[2]) -
                                        (HSX[1]-xi[1])*(HSX[2]-xj[2]) )*aux
      end

    end

    if dot_with!=nothing
      out[ti] += strength/(4*pi)*dot( V[1]*unitxi + V[2]*uniteta + V[3]*unitz,
                                                                  dot_with[ti])
    else
      out[ti] .+= strength/(4*pi)*(V[1]*unitxi + V[2]*uniteta + V[3]*unitz)
    end

  end
end


################################################################################
# VORTEX ELEMENTS
################################################################################
"""
Returns the velocity induced by a vortex ring panel of vertices `nodes` and
vortex strength `strength` on the targets `targets`. It adds the velocity at the
i-th target to out[i].
"""
function Vvortexring(nodes::Array{Arr1,1}, strength::RType,
                          targets::Array{Array{T2,1},1},
                          out;
                          dot_with::Union{Array{Array{T3,1},1}, Nothing}=nothing,
                          closed_ring::Bool=true
                          ) where{Arr1<:AbstractArray, T2<:RType, T3<:RType}
  if size(out)!=size(targets)
    error("Invalid `out` argument."*
          " Expected size $(size(targets)), got $(size(out)).")
  end

  nn = size(nodes, 1)                      # Number of nodes

  # Iterates over targets
  for ti in 1:size(targets, 1)
    V = zeros(3)

    for i in 1:(nn - 1*!closed_ring)
      p1, p2 = nodes[i], nodes[i%nn + 1]
      r1 = targets[ti] - p1
      r2 = targets[ti] - p2
      crossr1r2 = cross(r1,r2)
      # # This if statement avoids the singularity at the vortex line
      # if dot(crossr1r2,crossr1r2) > SMOOTH2*SMOOTH2
      #   V += crossr1r2/dot(crossr1r2,crossr1r2) * dot(
      #                                         (p1-p2), r1/norm(r1) - r2/norm(r2) )
      # end
        V += crossr1r2/(dot(crossr1r2,crossr1r2)+SMOOTH5) * dot(
                                              (p1-p2), r1/(norm(r1)+SMOOTH5) - r2/(norm(r2)+SMOOTH5) )
    end

    if dot_with!=nothing
      out[ti] -= strength/(4*pi)*dot(V, dot_with[ti])
    else
      out[ti] .-= strength/(4*pi)*V
    end

  end
end


"""
Returns the velocity induced by a semi-infinite vortex starting at point `p` in
the unitary direction `D` and vortex strength `strength` on the targets
`targets`. It adds the velocity at the i-th target to out[i].
"""
function Vsemiinfinitevortex(p::Array{T1,1}, D::Array{T2,1}, strength::RType,
                              targets::Array{Array{T3,1},1},
                              out;
                              dot_with::Union{Array{Array{T3,1},1}, Nothing}=nothing,
                              check::Bool=true
                              ) where{T1<:RType, T2<:RType, T3<:RType}
  # ERROR CASES
  if size(out)!=size(targets)
    error("Invalid `out` argument."*
          " Expected size $(size(targets)), got $(size(out)).")
  elseif check && abs(norm(D)-1.0)>1e-8
    error("Received non-unitary infinite direction D! (norm(D) = $(norm(D)))")
  end

  # Iterates over targets
  for ti in 1:size(targets, 1)

    p2 = p + dot(targets[ti]-p, D)*D
    h = norm(targets[ti]-p2)

    if h>SMOOTH3

      # Adds semi-infinite section
      if dot_with!=nothing
        out[ti] += dot(
                   strength / (4*pi*h) * cross(D, (targets[ti]-p2)/h ),
                                                                  dot_with[ti])
      else
        out[ti] .+= strength / (4*pi*h) * cross(D, (targets[ti]-p2)/h )
      end

      # Adds bound vortex section
      Vvortexring([p, p2], strength, targets[ti:ti], view(out, ti:ti);
                                          dot_with=dot_with!=nothing ? dot_with[ti:ti] : nothing,
                                          closed_ring=false)
    end

  end
end
