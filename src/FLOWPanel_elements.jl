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
"""
function Vconstant_source(nodes::Array{Arr1,1}, strength::RType,
                          targets::Array{Array{T2,1},1},
                          # out::Array{Array{T,1},1}
                          out;
                          dot_with=nothing
                          ) where{Arr1<:AbstractArray, T2<:RType}
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
    V = zeros(T2, 3)
    dtheta = 2*pi

    nR0 = 0

    for i in 1:nn
      xi, xj = HSnodes[i], HSnodes[i%nn + 1]

      dij = norm(xj-xi)
      ri = norm(HSX-xi)
      rj = norm(HSX-xj)

      #   println("ri,rj,dij=$ri,$rj,$dij")

      Qij = log( (ri+rj+dij)/(ri+rj-dij + SMOOTH) )

      Sij = (xj[2]-xi[2])/dij
      Cij = (xj[1]-xi[1])/dij

      siji = (xi[1]-HSX[1])*Cij + (xi[2]-HSX[2])*Sij
      sijj = (xj[1]-HSX[1])*Cij + (xj[2]-HSX[2])*Sij
      Rij = (HSX[1]-xi[1])*Sij - (HSX[2]-xi[2])*Cij

      Jij = atan( Rij*abs(HSX[3])*( ri*sijj - rj*siji ) ,
                   ri*rj*Rij^2 + HSX[3]^2*sijj*siji)

      #  println("Rij,ri,rj,siji,sijj,HSX=$Rij,$ri,$rj,$siji,$sijj,$HSX")
      #  println("Rij,ri,rj,siji,sijj,HSX,Sij,Cij=$Rij,$ri,$rj,$siji,$sijj,$HSX,$Sij,$Cij")

      V[1] -= Sij*Qij
      V[2] += Cij*Qij
      V[3] -= Jij

      dtheta *= Rij>=0
      nR0 += Rij==0
    end

    V[3] += dtheta
    V[3] *= sign(HSX[3])   # Isn't this sign already accounted for in atan2?
    V[3] *= !(nR0>1)       # Singularity fix of any z position aligned with node

    if dot_with!=nothing
      out[ti] += dot( strength*(V[1]*unitxi + V[2]*uniteta + V[3]*unitz),
                                                                  dot_with[ti])
    else
      out[ti] .+= strength*(V[1]*unitxi + V[2]*uniteta + V[3]*unitz)
    end

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
