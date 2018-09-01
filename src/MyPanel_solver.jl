"""
  Three-dimensional panel method solver. All these solvers asumes that every
  panel is a planar polygon of any number of vertices higher or equal to three.

  # AUTHORSHIP
    * Author    : Eduardo J. Alvarez
    * Email     : Edo.AlvarezR@gmail.com
    * Created   : Jul 2018
    * License   : AGPL-3.0
"""
module PanelSolver

# GeometricTools from https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
const gt = GeometricTools


const SMOOTH = 1e-2               # Smoothing radius of source panel
const SMOOTH2 = 2e-3              # Smoothing radius for vortex ring
# const SMOOTH2 = 1e-32              # Smoothing radius for vortex ring
const SMOOTH3 = SMOOTH2           # Smoothing radius for semi-infinite vortex



"""
Calculates and returns the geometric matrix of a collection of constant-source
panels on the given control points with their associated normals.

**ARGUMENTS**
  * `nodes::Array{T,2}`                 : All nodes in the collection of
                                          panels.
  * `panels::Array{Array{Int64,1},1}`   : Node connectivity data defining each
                                          panel, where `panels[i][j]` is the
                                          index in `nodes` of the j-th node in
                                          the i-th panel.
  * `CPs::Array{Array{Float64,1},1}`    : Control points.
  * `normals::Array{Array{Float64,1},1}`: Normal associated to every CP.
"""
function G_constant_source(nodes::Array{T,2},
                                panels::Array{Array{Int64,1},1},
                                CPs::Array{Array{T,1},1},
                                normals::Array{Array{T,1},1}) where{T<:Real}
  N = size(panels, 1)
  G = zeros(N, N)

  # Builds geometric matrix
  for j in 1:N # Iterates over columns (panels)
      Vconstant_source(
                        [nodes[:,ind] for ind in panels[j]], # Nodes in j-th panel
                        1.0,                               # Unitary strength,
                        CPs,                               # Targets
                        view(G, :, j);                     # Velocity of j-th
                                                           # panel on every CP
                        dot_with=normals                   # Normal of every CP
                      )
  end

  return G
end


"""
Returns the velocity induced by a panel of vertices `nodes` and constant
strength source `strength` on the targets `targets`. It adds the velocity at the
i-th target to out[i].
"""
function Vconstant_source(nodes::Array{Array{T,1},1}, strength::Real,
                          targets::Array{Array{T,1},1},
                          # out::Array{Array{T,1},1}
                          out;
                          dot_with=nothing
                          ) where{T<:Real}
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
    V = zeros(3)
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

      Jij = atan2( Rij*abs(HSX[3])*( ri*sijj - rj*siji ) ,
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
      out[ti] += strength*(V[1]*unitxi + V[2]*uniteta + V[3]*unitz)
    end

  end
end



"""
Calculates and returns the geometric matrix of a collection of vortex-ring
panels on the given control points with their associated normals. It models
the wake as a rigid steady wake.

**ARGUMENTS**
  * `nodes::Array{T,2}`                 : All nodes in the collection of
                                          panels.
  * `panels::Array{Array{Int64,1},1}`   : Node connectivity data defining each
                                          panel, where `panels[i][j]` is the
                                          index in `nodes` of the j-th node in
                                          the i-th panel.
  * `U::Array{Int64,1}`                 : Indices of all panels along the upper
                                          side of the trailing edge (from left
                                          to right).
  * `L::Array{Int64,1}`                 : Indices of all panels along the lower
                                          side of the trailing edge (from left
                                          to right).
  * `D::Array{Array{Float64,1},1}`      : Unitary direction of every
                                          semi-infinite vortex.
  * `CPs::Array{Array{Float64,1},1}`    : Control points.
  * `normals::Array{Array{Float64,1},1}`: Normal associated to every CP.
"""
function G_vortexring_rigid(nodes::Array{T,2},
                                panels::Array{Array{Int64,1},1},
                                CPs::Array{Array{T,1},1},
                                normals::Array{Array{T,1},1}) where{T<:Real}
  # N = size(panels, 1)
  # G = zeros(N, N)
  #
  # # Builds geometric matrix
  # for j in 1:N # Iterates over columns (panels)
  #     Vconstant_source(
  #                       [nodes[:,ind] for ind in panels[j]], # Nodes in j-th panel
  #                       1.0,                               # Unitary strength,
  #                       CPs,                               # Targets
  #                       view(G, :, j);                     # Velocity of j-th
  #                                                          # panel on every CP
  #                       dot_with=normals                   # Normal of every CP
  #                     )
  # end
  #
  # return G
end


"""
Calculates and returns the geometric matrix of a collection of vortex-ring
panels on the given control points with their associated normals.

**ARGUMENTS**
  * `nodes::Array{T,2}`                 : All nodes in the collection of
                                          panels.
  * `panels::Array{Array{Int64,1},1}`   : Node connectivity data defining each
                                          panel, where `panels[i][j]` is the
                                          index in `nodes` of the j-th node in
                                          the i-th panel.
  * `CPs::Array{Array{Float64,1},1}`    : Control points.
  * `normals::Array{Array{Float64,1},1}`: Normal associated to every CP.
"""
function G_vortexring(nodes::Array{T,2},
                                panels::Array{Array{Int64,1},1},
                                CPs::Array{Array{T,1},1},
                                normals::Array{Array{T,1},1}) where{T<:Real}
  N = size(panels, 1)
  G = zeros(N, N)

  # Builds geometric matrix
  for j in 1:N # Iterates over columns (panels)
      Vvortexring(
                        [nodes[:,ind] for ind in panels[j]], # Nodes in j-th panel
                        1.0,                               # Unitary strength,
                        CPs,                               # Targets
                        view(G, :, j);                     # Velocity of j-th
                                                           # panel on every CP
                        dot_with=normals                   # Normal of every CP
                      )
  end

  return G
end

"""
Returns the velocity induced by a vortex ring panel of vertices `nodes` and
vortex strength `strength` on the targets `targets`. It adds the velocity at the
i-th target to out[i].
"""
function Vvortexring(nodes::Array{Array{T,1},1}, strength::Real,
                          targets::Array{Array{T,1},1},
                          out;
                          dot_with=nothing, closed_ring::Bool=true
                          ) where{T<:Real}
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
      # This if statement avoids the singularity at the vortex line
      if dot(crossr1r2,crossr1r2) > SMOOTH2^2
        V += crossr1r2/dot(crossr1r2,crossr1r2) * dot(
                                              p1-p2, r1/norm(r1) - r2/norm(r2) )
      end
    end

    if dot_with!=nothing
      out[ti] += strength/(4*pi)*dot(V, dot_with[ti])
    else
      out[ti] += strength/(4*pi)*V
    end

  end
end


"""
Returns the velocity induced by a semi-infinite vortex starting at point `p` in
the unitary direction `D` and vortex strength `strength` on the targets
`targets`. It adds the velocity at the i-th target to out[i].
"""
function Vsemiinfinitevortex(p::Array{T,1}, D::Array{T,1}, strength::Real,
                              targets::Array{Array{T,1},1},
                              out;
                              dot_with=nothing, check=true
                              ) where{T<:Real}
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
        out[ti] += strength / (4*pi*h) * cross(D, (targets[ti]-p2)/h )
      end

      # Adds bound vortex section
      Vvortexring([p2, p], strength, targets[ti:ti], view(out, ti:ti);
                                          dot_with=dot_with, closed_ring=false)
    end
  end
end




end # END OF MODULE
