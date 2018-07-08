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


"""
  Returns the velocity induced by a
"""
function constant_source(nodes::Array{Array{T},1}, strength::Real,
                         targets::Array{Array{T},1}; out=nothing) where{T<:Real}
  if out==nothing
    _out = zeros(Real, 3, size(targets,1))
  elseif size(out)!=(3, size(targets, 1))
    error("Invalid `out` argument."*
          " Expected size (3,$(size(targets,1))), got $(size(out)).")
  else
    _out = out
  end

  nn = size(nodes, 1)                      # Number of nodes

  # Tangent, oblique, and normal vectors
  t, o, n = gt.get_unitvectors(nodes)

  # Coordinate system defined by Hess & Smith
  unitxi, unittheta, unitz = o, t, -n      # Unit vectors
  O = nodes[1]                             # Origin
  Oaxis = hcat(unitxi, unittheta, unitz)'  # Transformation matrix

  # Converts nodes to H&S coordinate system
  HSnodes = [Oaxis*(node-O) for node in nodes]

  # Iterates over targets
  for ti in 1:size(targets, 1)
    HSX = Oaxis*(targets[ti]-O)
    V = zeros(3)
    dtheta = 2*pi

    for i in 1:nn
      xi, xj = HSnodes[i], HSnodes[i%nn + 1]

      dij = norm(xj-xi)
      ri = norm(HSX-xi)
      rj = norm(HSX-xj)

      Sij = (xj[2]-xi[2])/dij
      Cij = (xj[1]-xi[1])/dij
      Qij = log( (ri+rj+dij)/(ri+rj-dij) )

      siji = (xi[1]-HSX[1])*Cij + (xi[2]-HSX[2])*Sij
      sijj = (xj[1]-HSX[1])*Cij + (xj[2]-HSX[2])*Sij
      Rij = (HSX[1]-xi[1])*Sij - (HSX[2]-xi[2])*Cij
      Jij = atan2( Rij*abs(HSX[3])*( ri*sijj - rj*sijj ) ,
                   ri*rj*Rij^2 + HSX[3]^2*sijj*siji)

      V[1] -= Sij*Qij
      V[2] += Cij*Qij
      V[3] -= Jij

      dtheta *= Rij>=0
    end

    V[3] += dtheta
    V[3] *= sign(HSX[3])

    _out[:, ti] += V
  end

  return _out
end









end # END OF MODULE
