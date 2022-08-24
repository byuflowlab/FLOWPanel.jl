#=##############################################################################
# DESCRIPTION
    Functions for generating geometry matrices (a.k.a., influence matrices).

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jul 2018
  * License   : MIT License
=###############################################################################





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
function G_constant_source(nodes::Array{T1,2},
                                panels::Array{Array{Int64,1},1},
                                CPs::Array{Array{T2,1},1},
                                normals::Array{Array{T3,1},1}
                            ) where{T1<:RType, T2<:RType, T3<:RType}
  N = size(panels, 1)
  # G = Array{Float64}(undef, N, N)
  G = zeros(N, N)

  # Builds geometric matrix
  for j in 1:N # Iterates over columns (panels)
      U_constant_source(
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
Calculates and returns the geometric matrix of a collection of constant-doublet
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
function G_constant_doublet(nodes::Array{T1,2},
                                panels::Array{Array{Int64,1},1},
                                CPs::Array{Array{T2,1},1},
                                normals::Array{Array{T3,1},1}
                            ) where{T1<:RType, T2<:RType, T3<:RType}
  N = size(panels, 1)
  # G = Array{Float64}(undef, N, N)
  G = zeros(N, N)

  # Builds geometric matrix
  for j in 1:N # Iterates over columns (panels)
      U_constant_doublet(
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
Calculates and returns the geometric matrix of a collection of vortex-ring
panels of a lifting body on the given control points with their associated
normals.

**ARGUMENTS**
  * `nodes::Array{T,2}`                 : All nodes in the collection of
                                          panels.
  * `panels::Array{Array{Int64,1},1}`   : Node connectivity data defining each
                                          panel, where `panels[i][j]` is the
                                          index in `nodes` of the j-th node in
                                          the i-th panel.
  * `U::Array{Int64,1}`                 : Indices of all panels along the upper
                                          side of the trailing edge.
  * `L::Array{Int64,1}`                 : Indices of all panels along the lower
                                          side of the trailing edge.
  * `CPs::Array{Array{Float64,1},1}`    : Control points.
  * `normals::Array{Array{Float64,1},1}`: Normal associated to every CP.
"""
function G_lifting_vortexring(nodes::Array{T1,2},
                                panels::Array{Array{Int64,1},1},
                                U::Array{Int64,1},
                                L::Array{Int64,1},
                                CPs::Array{Array{T2,1},1},
                                normals::Array{Array{T3,1},1}
                            ) where{T1<:RType, T2<:RType, T3<:RType}

  # Sorts trailing edge indices
  _U = sort(U)
  _L = sort(L)

  cur_u = 1                 # Index of current upper trailing edge cell
  cur_l = 1                 # Index of current lower trailing edge cell

  N = size(panels, 1)
  # G = Array{Float64}(undef, N, N)    # NOTE: This leads to unstable memory
                                       #  allocation when running this function
                                       #  multiple times. Prefer zeros(N, N)
                                       #  instead.
  G = zeros(N, N)

  # Builds geometric matrix --- Vortex rings
  for j in 1:N # Iterates over columns (panels)

    U_vortexring(# Nodes in j-th panel: puts first node last to match TE
                # [nodes[:,ind] for ind in vcat(panels[j][2:end], panels[j][1])],
                # Nodes in j-th panel: puts last node first to match TE
                [nodes[:,ind] for ind in vcat(panels[j][end], panels[j][1:end-1])],
                # Unitary strength,
                1.0,
                # Targets
                CPs,
                # Velocity of j-th panel on every CP
                view(G, :, j);
                # Normal of every CP
                dot_with=normals,
                # Checks for TE
                closed_ring= !( (cur_u<=size(_U,1) && j==_U[cur_u]) ||
                                (cur_l<=size(_L,1) && j==_L[cur_l]))
                )

    # if !prod(isnan.(view(G, :, j)).==false); println("j=$j\t$(findall(x->isnan(x), view(G, :, j)))"); end;

      if cur_u<=size(_U,1) && j==_U[cur_u]; cur_u+=1; end;
      if cur_l<=size(_L,1) && j==_L[cur_l]; cur_l+=1; end;
  end


  # println(!prod(isnan.(G)).==false)
  # println(length(filter(x-> isnan(x), G)))

  return G
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
function G_vortexring(nodes::Array{T1,2},
                                panels::Array{Array{Int64,1},1},
                                CPs::Array{Array{T2,1},1},
                                normals::Array{Array{T3,1},1}
                        ) where{T1<:RType, T2<:RType, T3<:RType}
  N = size(panels, 1)
  # G = Array{Float64}(undef, N, N)
  G = zeros(N, N)

  # Builds geometric matrix
  for j in 1:N # Iterates over columns (panels)
      U_vortexring(
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
