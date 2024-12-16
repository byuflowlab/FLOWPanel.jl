#=##############################################################################
# DESCRIPTION
    Definition of non-lifting paneled body types (implementations of
    AbstractLiftingBody).

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Sep 2018
  * License     : MIT License
=###############################################################################


################################################################################
# RIGID-WAKE BODY TYPE
################################################################################
"""
  `RigidWakeBody{E::AbstractElement, N}(grid::gt.GridTriangleSurface,
shedding::Matrix{Int})`

Lifting body that is solver using a combination of N panel elements and a steady
rigid wake. `grid` is the grid surface (paneled geometry).

`shedding[:, i]` contains the information of the i-th edge along which to shed
the wake, where `shedding[1, i]` is the linear index of the panel shedding the
 wake, and `shedding[2:3, i]` are the indices of the nodes in that panel that
 make the edge. Since the wake is typically shed at the edge between two panels,
`shedding[3, i]` is the index of the partner panel (use -1 if none) and
`shedding[4:5, i]` are the node indices in that panel that make the edge.
The user must ensure that both edges are coincident, and the strength of the
wake is equal to the difference between the strengths of both panels.

  **Properties**
  * `nnodes::Int`                       : Number of nodes
  * `ncells::Int`                       : Number of cells
  * `fields::Vector{String}`            : Available fields (solutions)
  * `Oaxis::Matrix`                     : Coordinate system of body w.r.t. global
  * `O::Vector`                         : Origin of body w.r.t. global
  * `ncellsTE::Int`                     : Number of cells along trailing edge
  * `nnodesTE::Int`                     : Number of nodes along trailing edge

"""
struct RigidWakeBody{E, N} <: AbstractLiftingBody{E, N}

    # User inputs
    grid::gt.GridTriangleSurface              # Paneled geometry
    shedding::Array{Int, 2}                   # Indicates edges along which to
                                              # shed the wake
    # Properties
    nnodes::Int                               # Number of nodes
    ncells::Int                               # Number of cells
    nsheddings::Int                           # Number of shedding edges
    fields::Array{String, 1}                  # Available fields (solutions)
    Oaxis::Array{<:Number,2}                  # Coordinate system of original grid
    O::Array{<:Number,1}                      # Position of CS of original grid

    # Internal variables
    strength::Array{<:Number, 2}              # strength[i,j] is the stength of the i-th panel with the j-th element type
    CPoffset::Float64                         # Control point offset in normal direction
    kerneloffset::Float64                     # Kernel offset to avoid singularities
    kernelcutoff::Float64                     # Kernel cutoff to avoid singularities
    characteristiclength::Function            # Characteristic length of each panel

    function RigidWakeBody{E, N}(
                                grid, shedding;
                                nnodes=grid.nnodes, ncells=grid.ncells,
                                nsheddings=size(shedding,2),
                                fields=Array{String,1}(),
                                Oaxis=Array(1.0I, 3, 3), O=zeros(3),
                                strength=zeros(grid.ncells, N),
                                CPoffset=1e-14,
                                kerneloffset=1e-8,
                                kernelcutoff=1e-14,
                                characteristiclength=characteristiclength_unitary,
                                check_mesh=true,
                                ) where {E, N}

        @assert _checkTE(grid, shedding) "Got invalid trailing edge"

        # Automated sanity checks for Meshes.jl mesh
        if check_mesh && typeof(grid.orggrid) <: gt.Meshes.Mesh

            mesh = grid.orggrid
            watertight = gt.isclosed(mesh)

            # Check that topology is consistent with the solver
            if watertight && E<:VortexRing && N==1
                @warn "Requested direct vortex ring solver on an closed mesh;"*
                " least-squares solver is recommended instead"*
                " (use `RigidWakeBody{VortexRing, 2}`)"
            end

            # Check that control points lay outside the geometry (based on
            # winding number)
            normals = _calc_normals(grid)
            controlpoints = _calc_controlpoints(grid, normals;
                                                off=CPoffset,
                                                characteristiclength=characteristiclength)

            (minw, maxw) = calc_minmax_winding(mesh, controlpoints)

            if abs(minw) > eps()^0.75 || abs(maxw) >= eps()^0.75

                if watertight
                    @warn "Found winding numbers other than 0, which might indicate"*
                    " that control points are inside the geometry; flipping the"*
                    " sign of `CPoffset` is recommended; (minw, maxw) ="*
                    " $((minw, maxw))"
                else
                    @warn "Found winding numbers other than 0, which might indicate"*
                    " that control points are inside the geometry; however,"*
                    " geometry is not watertight, so it might be ok."*
                    " (minw, maxw) = $((minw, maxw))"
                end
            end

        end

        return new(
                    grid, shedding,
                    nnodes, ncells,
                    nsheddings,
                    fields,
                    Oaxis, O,
                    strength,
                    CPoffset,
                    kerneloffset,
                    kernelcutoff,
                    characteristiclength
                    )

        end

end

function (RigidWakeBody{E})(grid, shedding; optargs...) where {E}
    return RigidWakeBody{E, _count(E)}(grid, shedding; optargs...)
end

function (RigidWakeBody{E})(grid; optargs...) where {E}
    return RigidWakeBody{E}(grid, zeros(Int, 6, 0); optargs...)
end



function save(body::RigidWakeBody, args...;
                out_wake::Bool=true, debug::Bool=false,
                wake_len::Number=1.0,
                wake_panel::Bool=false,
                wake_suffix="_wake",
                optargs...)


    str = ""
    str *= save_base(body, args...; debug=debug, optargs...)

    # Output the wake
    if out_wake
        str *= _savewake(body, args...; len=wake_len, panel=wake_panel,
                            optargs..., suffix=wake_suffix)
    end

    return str
end




################################################################################
# VORTEX RING SOLVER
################################################################################
function solve(self::RigidWakeBody{VortexRing, 1},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2},
                Dbs::AbstractMatrix{T3};
                solver=solve_ludiv!, solver_optargs=(),
                optargs...
                ) where {T1, T2, T3}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings)), got $(size(Das))")
    elseif size(Dbs) != (3, self.nsheddings)
        error("Invalid Dbs;"*
              " expected size (3, $(self.nsheddings)), got $(size(Dbs))")
    end

    T = promote_type(T1, T2, T3)

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    G = zeros(T, self.ncells, self.ncells)
    _G_U!(self, G, CPs, normals, Das, Dbs; optargs...)

    # Calculate boundary conditions (right-hand side of system of equations)
    RHS = calc_bc_noflowthrough(Uinfs, normals)

    # Solve system of equations
    Gamma = zeros(T, self.ncells)
    solver(Gamma, G, RHS; solver_optargs...)

    # Save solution
    self.strength[:, 1] .= Gamma

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "Da", "vector", collect(eachcol(Das)), "system")
    add_field(self, "Db", "vector", collect(eachcol(Dbs)), "system")
    add_field(self, "Gamma", "scalar", view(self.strength, :, 1), "cell")
end

"""
Computes the geometric matrix (left-hand side matrix of the system of equation)
and stores it under `G`.

**ARGUMENTS**
  * `G::Array{T,2}`                     : Pre-allocated output memory.
  * `CPs::Array{T,2}`                   : Control points.
  * `normals::Array{T,2}`               : Normal associated to every CP.
  * `Das::Array{T,2}`                   : Unitary direction of semi-infinite
                                            vortex at point `a` of each
                                            trailing edge panel.
  * `Dbs::Array{T,2}`                   : Unitary direction of semi-infinite
                                            vortex at point b of each
                                            trailing edge panel.
"""
_G_U!(self::RigidWakeBody{VortexRing, 1}, args...; optargs...) = _G_Uvortexring!(self, args..., optargs...)







function solve(self::RigidWakeBody{VortexRing, 2},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2},
                Dbs::AbstractMatrix{T3};
                solver=solve_ludiv!, solver_optargs=(),
                elprescribe::AbstractArray{Tuple{Int, Float64}}=[(1, 0.0)],
                GPUArray=Array{promote_type(T1, T2, T3)},
                optargs...
                ) where {T1, T2, T3}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings)), got $(size(Das))")
    elseif size(Dbs) != (3, self.nsheddings)
        error("Invalid Dbs;"*
              " expected size (3, $(self.nsheddings)), got $(size(Dbs))")
    end

    T = promote_type(T1, T2, T3)

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix) and boundary
    # conditions (right-hand-side) converted into a least-squares problem
    G, RHS = _G_U_RHS(self, Uinfs, CPs, normals, Das, Dbs, elprescribe;
                                                GPUArray=GPUArray,
                                                optargs...)

    # Solve system of equations
    Gamma = GPUArray(undef, self.ncells-length(elprescribe))
    solver(Gamma, G, RHS; solver_optargs...)

    # Port solution back to CPU if solved in GPU
    if !(GPUArray <: Array)
        Gamma = Array{T}(Gamma)
    end

    # Save solution
    set_solution(self, nothing, Gamma, elprescribe, Uinfs, Das, Dbs)

end

calc_elprescribe(::RigidWakeBody{VortexRing, 2}) = [(1, 0.0)]

function set_solution(self::RigidWakeBody{VortexRing, 2},
                        dummy, Gammals, elprescribe,
                        Uinfs, Das, Dbs)

    # Save vortex ring circulations: add Gamma and prescribed strengths
    prev_eli = 0
    if length(elprescribe)==0
        self.strength[:, 1] .= Gammals
    else
        for (i, (eli, elval)) in enumerate(elprescribe)
            self.strength[(prev_eli+1):(eli-1), 1] .= view(Gammals, (prev_eli+2-i):(eli-i))
            self.strength[eli, 1] = elval

            if i==length(elprescribe) && eli!=size(self.strength, 1)
                self.strength[eli+1:end, 1] .= view(Gammals, (eli-i+1):size(Gammals, 1))
            end

            prev_eli = eli
        end
    end

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "Da", "vector", collect(eachcol(Das)), "system")
    add_field(self, "Db", "vector", collect(eachcol(Dbs)), "system")
    add_field(self, "Gamma", "scalar", view(self.strength, :, 1), "cell")
end

function _G_U_RHS(self::RigidWakeBody{VortexRing, 2}, args...; optargs...)
    return _G_U_RHS_leastsquares(self, args...; optargs...)
end

function _G_U_RHS!(self::RigidWakeBody{VortexRing, 2}, args...; optargs...)
    return _G_U_RHS_leastsquares!(self, args...; optargs...)
end

function _G_U_RHS_leastsquares(self::AbstractBody,
                                Uinfs::AbstractMatrix{T1}, CPs, normals,
                                Das::AbstractMatrix{T2},
                                Dbs::AbstractMatrix{T3},
                                elprescribe::AbstractArray{Tuple{Int, T4}},
                                args...;
                                GPUArray=Array{promote_type(T1, T2, T3, T4)},
                                optargs...
                                ) where {T1, T2, T3, T4}

    T = promote_type(T1, T2, T3, T4)

    n = self.ncells
    npres = length(elprescribe)

    G = zeros(T, n, n)
    Gred = zeros(T, n, n-npres)
    tGred = zeros(T, n-npres, n)
    gpuGred = GPUArray(undef, size(Gred))
    Gls = GPUArray(undef, n-npres, n-npres)
    RHS = zeros(T, n)
    RHSls = GPUArray(undef, n-npres)

    _G_U_RHS_leastsquares!(self, G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
                Uinfs, CPs, normals, Das, Dbs,
                elprescribe,
                args...; optargs...)

    return Gls, RHSls
end

function _G_U_RHS_leastsquares!(self::AbstractBody,
                                G, Gred, tGred, gpuGred, Gls, RHS, RHSls,
                                Uinfs, CPs, normals,
                                Das, Dbs,
                                elprescribe::AbstractArray{Tuple{Int, T}};
                                optargs...
                                ) where {T<:Number}

    n = self.ncells
    npres = length(elprescribe)

    @assert size(G, 1)==n && size(G, 2)==n ""*
        "Invalid $(size(G, 1))x$(size(G, 2)) matrix G; expected $(n)x$(n)"
    @assert size(Gred, 1)==n && size(Gred, 2)==n-npres ""*
        "Invalid $(size(Gred, 1))x$(size(Gred, 2)) matrix Gred; expected $(n)x$(n-npres)"
    @assert size(tGred, 1)==n-npres && size(tGred, 2)==n ""*
        "Invalid $(size(tGred, 1))x$(size(tGred, 2)) matrix tGred; expected $(n-npres)x$(n)"
    @assert size(Gls, 1)==n-npres && size(Gls, 2)==n-npres ""*
        "Invalid $(size(Gls, 1))x$(size(Gls, 2)) matrix Gls; expected $(n-npres)x$(n-npres)"

    @assert length(RHS)==n "Invalid RHS length $(length(RHS)); expected $(n)"
    @assert length(RHSls)==n-npres "Invalid RHSls length $(length(RHSls)); expected $(n-npres)"

    # Sort prescribed elements by index
    sort!(elprescribe, by = x -> x[1])

    # Calculate normal velocity of freestream for boundary condition
    calc_bc_noflowthrough!(RHS, Uinfs, normals)

    # -------------- Influence of vortex rings -------------------------
    # Calculate influence of all vortex rings
    _G_Uvortexring!(self, G, CPs, normals, Das, Dbs; optargs...)

    # Move influence of prescribed vortex ring elements to right-hand side
    for (eli, elval) in elprescribe
        for i in 1:length(RHS)
            RHS[i] -= elval*G[i, eli]
            # G[i, eli] = 0
        end
    end

    # -------------- Least-squares problem ----------------------------
    # Gred = view(G, :, vcat(1:elprescribe_index-1, elprescribe_index+1:size(G, 2)))

    # Reduce G: copy G into Gred without the prescribed elements
    prev_eli = 0
    for (i, (eli, elval)) in enumerate(elprescribe)

        Gred[:, (prev_eli+2-i):(eli-i)] .= view(G, :, (prev_eli+1):(eli-1))

        if i==length(elprescribe) && eli!=size(G, 2)
            Gred[:, (eli-i+1):end] .= view(G, :, eli+1:size(G, 2))
        end

        prev_eli = eli
    end

    # Produce least-squares matrix
    if typeof(gpuGred) <: Array   # Case: CPU arrays

        # tGred = transpose(Gred)               # <- Very slow to multiply later on
        # tGred = collect(transpose(Gred))      # <- Much faster but allocating memory
        # tGred = permutedims(Gred)             # <- Ditto
        permutedims!(tGred, Gred, [2, 1])       # <- No memory allocation

        # RHSls = Gred'*RHS
        LA.mul!(RHSls, tGred, RHS)

        # Gls = Gred'*Gred
        LA.mul!(Gls, tGred, Gred)

    else                          # Case: GPU arrays

        copyto!(gpuGred, Gred)
        tGred = transpose(gpuGred)

        # RHSls = Gred'*RHS
        LA.mul!(RHSls, tGred, typeof(RHSls)(RHS))

        # Gls = Gred'*Gred
        LA.mul!(Gls, tGred, gpuGred)

    end

    return Gls, RHSls
end







function _G_Uvortexring!(self::RigidWakeBody,
                            G::Arr1, CPs::Arr2, normals::Arr3, Das, Dbs;
                            optargs...
                       ) where{ T1, Arr1<:AbstractArray{T1, 2},
                                T2, Arr2<:AbstractArray{T2, 2},
                                T3, Arr3<:AbstractArray{T3, 2}}

    N = self.ncells
    M = size(CPs, 2)

    if size(G, 1)!=M || size(G, 2)!=N
        error("Matrix G with invalid dimensions;"*
              " got $(size(G)), expected ($M, $N).")
    elseif size(normals, 2)!=M
        error("normals matrix with invalid dimensions;"*
              " got $(size(normals)), expected (3, $M).")
    end


    # Build geometric matrix from panel contributions
    panels = 1:self.ncells
    chunks = collect(Iterators.partition(panels, max(length(panels) ÷ Threads.nthreads(), 3*Threads.nthreads())))

    Threads.@threads for chunk in chunks      # Distribute panel iteration among all CPU threads

        # Pre-allocate memory for panel calculation
        tri_out, tricoor, quadcoor, quad_out, lin, ndivscells, cin = gt.generate_getcellt_args!(self.grid)

        # for (pj, Gslice) in enumerate(eachcol(G))
        for pj in chunk                       # Iterate over panels

            panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                self.grid, pj, lin, ndivscells, cin)

            U_vortexring(
                              self.grid._nodes,                  # All nodes
                              panel,                             # Indices of nodes that make this panel
                              1.0,                               # Unitary strength
                              CPs,                               # Targets
                              view(G, :, pj);                    # Velocity of j-th panel on every CP
                              # Gslice;
                              dot_with=normals,                  # Normal of every CP
                              offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                              cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                              optargs...
                             )

         end
    end


    # Add wake contributions
    sheddings = 1:self.nsheddings
    chunks = collect(Iterators.partition(sheddings, max(length(sheddings) ÷ Threads.nthreads(), 3*Threads.nthreads())))

    Threads.@threads for chunk in chunks        # Distribute wake panel iteration among all CPU threads

        # Pre-allocate memory for panel calculation
        TE = zeros(Int, 2)
        _tri_out, _tricoor, _quadcoor, _quad_out, _lin, _ndivscells, _cin = gt.generate_getcellt_args!(self.grid)

        # for (ei, (pi, nia, nib, pj, nja, njb)) in enumerate(eachcol(self.shedding))
        for ei in chunk                          # Iterate over wake-shedding panels

            pi, nia, nib, pj, nja, njb = view(self.shedding, :, ei)

            # Fetch nodes of upper wake panel
            panel = gt.get_cell_t!(_tri_out, _tricoor, _quadcoor, _quad_out,
                                                self.grid, pi, _lin, _ndivscells, _cin)

            # Indicate nodes in the upper shedding edge
            TE[1] = panel[nia]
            TE[2] = panel[nib]
            da1, da2, da3 = Das[1, ei], Das[2, ei], Das[3, ei]
            db1, db2, db3 = Dbs[1, ei], Dbs[2, ei], Dbs[3, ei]

            U_semiinfinite_horseshoe(
                              self.grid._nodes,                  # All nodes
                              TE,                                # Indices of nodes that make the shedding edge
                              da1, da2, da3,                     # Semi-infinite direction da
                              db1, db2, db3,                     # Semi-infinite direction db
                              1.0,                               # Unitary strength
                              CPs,                               # Targets
                              view(G, :, pi);                    # Velocity of upper wake panel on every CP
                              dot_with=normals,                  # Normal of every CP
                              offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                              cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                              optargs...
                             )

             if pj != -1
                 # Fetch nodes of lower wake panel
                 panel = gt.get_cell_t!(_tri_out, _tricoor, _quadcoor, _quad_out,
                                                     self.grid, pj, _lin, _ndivscells, _cin)

                 # Indicate nodes in the lower shedding edge
                 TE[1] = panel[nja]
                 TE[2] = panel[njb]

                 U_semiinfinite_horseshoe(
                                   self.grid._nodes,                  # All nodes
                                   TE,                                # Indices of nodes that make the shedding edge
                                   db1, db2, db3,                     # Semi-infinite direction da (flipped in lower panel)
                                   da1, da2, da3,                     # Semi-infinite direction db
                                   1.0,                               # Unitary strength
                                   CPs,                               # Targets
                                   view(G, :, pj);                    # Velocity of lower wake panel on every CP
                                   dot_with=normals,                  # Normal of every CP
                                   offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                                   cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                                   optargs...
                                  )
             end
        end

    end

end

_Uind!(self::RigidWakeBody{VortexRing, N}, args...; optargs...) where {N} = _Uvortexring!(self, args...; stri=1, optargs...)

function _Uvortexring!(self::RigidWakeBody, targets, out; stri=1, optargs...)


    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    Das, Dbs = get_field(self, "Da")["field_data"], get_field(self, "Db")["field_data"]

    # Iterates over body panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Velocity of i-th panel on every target
        U_vortexring(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, stri],            # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                            optargs...
                         )
    end

    # Add wake contribution
    TE = zeros(Int, 2)
    for (ei, (pi, nia, nib, pj, nja, njb)) in enumerate(eachcol(self.shedding)) # Iterate over wake-shedding panels

        strengthi, strengthj = _get_wakestrength_mu(self, ei; stri=stri)

        # Fetch nodes of upper wake panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pi, lin, ndivscells, cin)

        # Indicate nodes in the upper shedding edge
        TE[1] = panel[nia]
        TE[2] = panel[nib]
        da1, da2, da3 = Das[ei]
        db1, db2, db3 = Dbs[ei]

        U_semiinfinite_horseshoe(
                          self.grid._nodes,                  # All nodes
                          TE,                                # Indices of nodes that make the shedding edge
                          da1, da2, da3,                     # Semi-infinite direction da
                          db1, db2, db3,                     # Semi-infinite direction db
                          strengthi,                         # Strength
                          targets,                           # Targets
                          out;                               # Outputs
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                          optargs...
                         )

         if pj != -1
             # Fetch nodes of lower wake panel
             panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                 self.grid, pj, lin, ndivscells, cin)

             # Indicate nodes in the lower shedding edge
             TE[1] = panel[nja]
             TE[2] = panel[njb]

             U_semiinfinite_horseshoe(
                               self.grid._nodes,                  # All nodes
                               TE,                                # Indices of nodes that make the shedding edge
                               db1, db2, db3,                     # Semi-infinite direction da (flipped in lower panel)
                               da1, da2, da3,                     # Semi-infinite direction db
                               strengthj,                         # Unitary strength
                               targets,                           # Targets
                               out;                               # Outputs
                               offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                               cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                               optargs...
                              )
         end
    end
end


function _phi!(self::RigidWakeBody{VortexRing, N}, targets, out; optargs...) where {N}

    # Pre-allocate memory for panel calculation
    lin = LinearIndices(self.grid._ndivsnodes)
    ndivscells = vcat(self.grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in self.grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    Das, Dbs = get_field(self, "Da")["field_data"], get_field(self, "Db")["field_data"]

    # Iterates over body panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Potential of i-th panel on every target
        phi_constant_doublet(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, 1],               # Unitary strength
                            targets,                           # Targets
                            out;                               # Outputs
                            optargs...
                         )
    end

    # Add wake contribution
    TE = zeros(Int, 2)
    for (ei, (pi, nia, nib, pj, nja, njb)) in enumerate(eachcol(self.shedding)) # Iterate over wake-shedding panels

        strengthi, strengthj = _get_wakestrength_mu(self, ei)

        # Fetch nodes of upper wake panel
        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pi, lin, ndivscells, cin)

        # Indicate nodes in the upper shedding edge
        TE[1] = panel[nia]
        TE[2] = panel[nib]
        da1, da2, da3 = Das[ei]
        db1, db2, db3 = Dbs[ei]

        phi_semiinfinite_doublet(
                          self.grid._nodes,                  # All nodes
                          TE,                                # Indices of nodes that make the shedding edge
                          da1, da2, da3,                     # Semi-infinite direction da
                          db1, db2, db3,                     # Semi-infinite direction db
                          strengthi,                         # Strength
                          targets,                           # Targets
                          out;                               # Outputs
                          optargs...
                         )

         if pj != -1
             # Fetch nodes of lower wake panel
             panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                                 self.grid, pj, lin, ndivscells, cin)

             # Indicate nodes in the lower shedding edge
             TE[1] = panel[nja]
             TE[2] = panel[njb]

             phi_semiinfinite_doublet(
                               self.grid._nodes,                  # All nodes
                               TE,                                # Indices of nodes that make the shedding edge
                               db1, db2, db3,                     # Semi-infinite direction da (flipped in lower panel)
                               da1, da2, da3,                     # Semi-infinite direction db
                               strengthj,                         # Unitary strength
                               targets,                           # Targets
                               out;                               # Outputs
                               optargs...
                              )
         end
    end
end

_get_Gdims(self::RigidWakeBody{VortexRing, N}) where {N} = (self.ncells, self.ncells)


################################################################################
# VORTEX RING + UNIFORM VORTEX SHEET SOLVER
################################################################################

function (RigidWakeBody{Union{VortexRing, UniformVortexSheet}})(grid, args...; optargs...)
    return RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 2}(grid, args...;
                                                    strength=zeros(grid.ncells, 3), optargs...)
end

function solve(self::RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 2},
                Uinfs::AbstractMatrix{T1},
                Das::AbstractMatrix{T2},
                Dbs::AbstractMatrix{T3};
                solver=solve_ludiv!, solver_optargs=(),
                elprescribe_index::Int=1, elprescribe_value=0,
                weight_gammat=0, weight_gammao=1
                ) where {T1, T2, T3}

    if size(Uinfs) != (3, self.ncells)
        error("Invalid Uinfs;"*
              " expected size (3, $(self.ncells)), got $(size(Uinfs))")
    elseif size(Das) != (3, self.nsheddings)
        error("Invalid Das;"*
              " expected size (3, $(self.nsheddings)), got $(size(Das))")
    elseif size(Dbs) != (3, self.nsheddings)
        error("Invalid Dbs;"*
              " expected size (3, $(self.nsheddings)), got $(size(Dbs))")
    end

    T = promote_type(T1, T2, T3)

    # Compute normals and control points
    normals = _calc_normals(self)
    CPs = _calc_controlpoints(self, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    # and boundary conditions (right-hand side of system of equations)
    G = zeros(T, self.ncells, self.ncells)
    RHS = zeros(T, self.ncells)

    _G_U_RHS!(self, G, RHS, Uinfs, CPs, normals, Das, Dbs,
                elprescribe_index, elprescribe_value,
                weight_gammat, weight_gammao)

    # Solve system of equations
    Gamma = zeros(T, self.ncells)
    solver(Gamma, G, RHS; solver_optargs...)

    # Save vortex ring circulations
    self.strength[:, 1] .= Gamma
    self.strength[elprescribe_index, 1] = elprescribe_value

    # Save vortex sheet strength
    gamma = Gamma[elprescribe_index]
    self.strength[:, 2] .= gamma*weight_gammat
    self.strength[1:2:end, 2] .*= -1
    self.strength[:, 3] .= gamma*weight_gammao
    self.strength[1:2:end, 3] .*= -1

    _solvedflag(self, true)
    add_field(self, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(self, "Da", "vector", collect(eachcol(Das)), "system")
    add_field(self, "Db", "vector", collect(eachcol(Dbs)), "system")
    add_field(self, "Gamma", "scalar", view(self.strength, :, 1), "cell")

    tangents = _calc_tangents(self)
    obliques = _calc_obliques(self)
    aux = zip(eachcol(tangents), eachcol(obliques),
                view(self.strength, :, 2), view(self.strength, :, 3))
    gammas = [gammat*t + gammao*o for (t, o, gammat, gammao) in aux]
    add_field(self, "gamma", "vector", gammas, "cell")
end

function _G_U_RHS!(self::RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 2},
                    G, RHS, Uinfs, CPs, normals, Das, Dbs,
                    elprescribe_index::Int, elprescribe_value::Number,
                    weight_gammat::Number, weight_gammao::Number;
                    optargs...
                    )

    # Calculate normal velocity of freestream for boundary condition
    calc_bc_noflowthrough!(RHS, Uinfs, normals)

    # -------------- Influence of vortex rings -------------------------

    # Calculate influence of all vortex rings
    _G_Uvortexring!(self, G, CPs, normals, Das, Dbs; optargs...)

    # Move influence of prescribed vortex ring element to right-hand side
    for i in 1:length(RHS)
        RHS[i] -= elprescribe_value*G[i, elprescribe_index]
        G[i, elprescribe_index] = 0
    end

    # -------------- Influence of vortex sheet ------------------------

    # Pre-allocate memory for panel calculation
    (tri_out, tricoor, quadcoor, quad_out,
                lin, ndivscells, cin) = gt.generate_getcellt_args!(self.grid)

    # Influence of vortex sheet on each CP gets stored here
    Gslice = view(G, :, elprescribe_index)

    # Calculate influence of each panel on each CP
    for pj in 1:self.ncells         # Iterate over panels

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                            self.grid, pj, lin, ndivscells, cin)

        s = pj%2==1 ? -1 : 1        # Alternate + and - strengths to get them all aligned

        U_constant_vortexsheet(
                          self.grid._nodes,                  # All nodes
                          panel,                             # Indices of nodes that make this panel
                          s*weight_gammat,                   # Tangential strength
                          s*weight_gammao,                   # Oblique strength
                          CPs,                               # Targets
                          # view(G, :, pj);                  # Agglomerate velocity of j-th panel on every CP
                          Gslice;
                          dot_with=normals,                  # Normal of every CP
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                          optargs...
                         )
    end

    return G, RHS
end

function _Uind!(self::RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 2},
                                            targets, out; optargs...)

    # Velocity induced by vortex rings
    _Uvortexring!(self, targets, out; stri=1, optargs...)

    # Velocity induced by vortex sheets
    _Uconstantvortexsheet!(self, targets, out; strti=2, stroi=3, optargs...)
end


function _Uconstantvortexsheet!(self::RigidWakeBody, targets, out;
                                                strti=2, stroi=3, optargs...)

    # Pre-allocate memory for panel calculation
    (tri_out, tricoor, quadcoor, quad_out,
                lin, ndivscells, cin) = gt.generate_getcellt_args!(self.grid)

    # Iterates over body panels
    for i in 1:self.ncells

        panel = gt.get_cell_t!(tri_out, tricoor, quadcoor, quad_out,
                                             self.grid, i, lin, ndivscells, cin)

        # Velocity of i-th panel on every target
        U_constant_vortexsheet(
                            self.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            self.strength[i, strti],           # Tangential strength
                            self.strength[i, stroi],           # Oblique strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                            optargs...
                         )
    end

end


function _phi!(self::RigidWakeBody{Union{VortexRing, UniformVortexSheet}, 2}, args...; optargs...)
    nothing
end




################################################################################
# COMMON FUNCTIONS
################################################################################
"""
    `calc_shedding(grid::GridTriangleSurface{gt.Meshes.SimpleMesh},
trailingedge::Matrix; tolerance=1e2*eps())`

Given an unstructured `grid` and a collection of points (line) `trailingedge`,
it finds the points in `grid` that are closer than `tolerance` to the line,
and automatically builds a `shedding` matrix that can be used to shed the wake
from this trailing edge.

Note: It is important that the points in `trailingedge` have been previously
    sorted to be contiguous to each other, otherwise the resulting `shedding`
    might have panels that are not contiguous to each other, fail to recognize
    panels that are at the trailing edge, or unphysically large trailing
    vortices.
"""
function calc_shedding(grid::gt.GridTriangleSurface{G}, trailingedge::Union{Matrix, Function};
                            periodic::Bool=false,
                            tolerance=1e2*eps(), debug=false
                            ) where {G<:gt.Meshes.SimpleMesh}

    nodes = grid._nodes
    topology = grid._halfedgetopology
    connec = grid.orggrid.topology.connec

    # Identify the nodes that are on the TE line
    TEindices = gt.identifyedge(nodes, trailingedge; tolerance=tolerance)
    TEindices = [nodei for (nodei, pointi) in TEindices]

    # Return if no TE nodes were identified
    if length(TEindices)==0
        return noshedding
    end

    # Append first node at end if expected to be periodic
    if periodic
        push!(TEindices, TEindices[1])
    end

    # All node pairs that could form an edge at the TE
    paircandidates = zip(view(TEindices, 1:length(TEindices)-1), view(TEindices, 2:length(TEindices)))

    # All node pairs that actually form an edge at the TE
    pairs = [pair for pair in paircandidates if haskey(topology.edge4pair, pair)]

    # Fetch all the first halfedge of each edges (node pairs) along the TE
    halfedges = [gt.Meshes.half4pair(topology, pair) for pair in pairs]

    # Build shedding matrix
    shedding = zeros(Int, 6, length(halfedges))

    for (ei, halfedge) in enumerate(halfedges)

        pair = pairs[ei]

        # pi is the panel associated with this half edge
        # pj is the panel associated with the other half

        # Case: Single-sided edge
        if isnothing(halfedge.elem) || isnothing(halfedge.half.elem)

            if isnothing(halfedge.half.elem)
                pi = halfedge.elem
            else
                pi = halfedge.half.elem
            end

            # Declare the other half as inexistent
            pj = -1

        # Case: Two-sided edge
        else

            # Identify which panel is "on top" and which "bottom" by matching
            # the order of the node pair
            inds1 = connec[halfedge.elem].indices
            inds2 = connec[halfedge.half.elem].indices

            if (
                    (inds1[1]==pair[1] && inds1[2]==pair[2])
                    ||
                    (inds1[2]==pair[1] && inds1[3]==pair[2])
                    ||
                    (inds1[3]==pair[1] && inds1[1]==pair[2])
                )

                pi = halfedge.elem
                pj = halfedge.half.elem

            elseif (
                    (inds2[1]==pair[1] && inds2[2]==pair[2])
                    ||
                    (inds2[2]==pair[1] && inds2[3]==pair[2])
                    ||
                    (inds2[3]==pair[1] && inds2[1]==pair[2])
                )

                pi = halfedge.half.elem
                pj = halfedge.elem

            else
                error("Logic error: Could not match panel to node pair")
            end

        end

        # Nodes of first half
        nia = findfirst(globindex -> globindex==pair[2], connec[pi].indices)  # Local-index of the first node
        nib = findfirst(globindex -> globindex==pair[1], connec[pi].indices)  # Local-index of the second node

        # Nodes of other half
        if pj != -1
            nja = findfirst(globindex -> globindex==pair[1], connec[pj].indices)  # Local-index of the second node
            njb = findfirst(globindex -> globindex==pair[2], connec[pj].indices)  # Local-index of the first node
        else
            nja = njb = -1
        end

        shedding[:, ei] .= (pi, nia, nib, pj, nja, njb)

    end

    if debug
        display("TEindices")
        display(TEindices)
        display("paircandidates")
        display(collect(paircandidates))
        display("pairs")
        display(pairs)
        display("halfedges")
        display(halfedges)
        display("shedding")
        display(shedding)

        points = [nodes[:, i] for i in TEindices]
        gt.generateVTK("TEindices", points; keep_points=true)

        lines = [[i-1, j-1] for (i, j) in paircandidates]
        gt.generateVTK("paircandidates", eachcol(nodes); lines=lines)

        lines = [[i-1, j-1] for (i, j) in pairs]
        gt.generateVTK("pairs", eachcol(nodes); lines=lines)
    end

    return shedding
end





##### INTERNAL FUNCTIONS  ######################################################
function _get_wakestrength_mu(self::RigidWakeBody, i; stri=1)
    strength1 = self.strength[self.shedding[1, i], stri]
    strength2 = self.shedding[4, i] != -1 ? self.strength[self.shedding[4, i], stri] : 0
    return strength1, strength2
end
function _get_wakestrength_Gamma(self::RigidWakeBody, i; stri=1)
    strength1 = self.strength[self.shedding[1, i], stri]
    strength2 = self.shedding[4, i] != -1 ? self.strength[self.shedding[4, i], stri] : 0
    return strength1 - strength2
end

"""
Outputs a vtk file with the wake.
"""
function _savewake(self::RigidWakeBody, filename::String;
                    len::Number=1.0, panel::Bool=true, suffix="_wake", optargs...)

    if check_field(self, "Da")==false
        error("Requested to save wake, but Da field was not found."*
              " Run solve(...) first")
    elseif check_field(self, "Db")==false
        error("Requested to save wake, but Db field was not found."*
              " Run solve(...) first")
    end


    nodes = self.grid._nodes
    nedges = size(self.shedding, 2)
    points = zeros(3, 4*nedges)

    # Collect points along shedding edge
    tricoor, quadcoor, lin, ndivscells, cin = gt.generate_getcellt_args(self.grid)

    for (si, (pi, nia, nib, pj)) in enumerate(eachcol(self.shedding))

        # Convert node indices from panel-local to global
        pia = gt.get_cell_t(tricoor, quadcoor, self.grid, pi, nia, lin, ndivscells, cin)
        pib = gt.get_cell_t(tricoor, quadcoor, self.grid, pi, nib, lin, ndivscells, cin)

        Da = get_fieldval(self, "Da", si)
        Db = get_fieldval(self, "Db", si)

        for i in 1:3
            # Shedding edge
            points[i, si + 0*nedges] = nodes[i, pia]
            points[i, si + 1*nedges] = nodes[i, pib]

            # End points of semi-infinite wake
            points[i, si + 2*nedges] = points[i, si + 0*nedges] + len*Da[i]
            points[i, si + 3*nedges] = points[i, si + 1*nedges] + len*Db[i]
        end

    end

    if panel    # Case that wake is saved as doublet panels

        arrpoints = collect(eachcol(points))
        str = ""

        for si in 1:2   # si==1, upper side; si==2, lower side
            cells = Array{Int64,1}[]
            cell_data = []

            # Wake panels
            for i in 0:nedges-1
                cell = [i+2*nedges, i+0*nedges, i+1*nedges, i+3*nedges]
                push!(cells, si==1 ? cell : reverse(cell))
            end

            # Strength of wake along shedding edge
            strengths = [_get_wakestrength_mu(self, i)[si] for i in 1:nedges]

            push!(cell_data, Dict( "field_name"=> "mu",
                                    "field_type"=> "scalar",
                                    "field_data"=> strengths))

            for cell in cells; cell .+= 1; end
            normals = [[gt._calc_n1(points, cell), gt._calc_n2(points, cell), gt._calc_n3(points, cell)] for cell in cells]
            for cell in cells; cell .-= 1; end

            push!(cell_data, Dict( "field_name"=> "normal",
                                    "field_type"=> "vector",
                                    "field_data"=> normals))

            # Generate VTK
            str *= gt.generateVTK(filename*suffix*(si==1 ? "_up" : "_lo"),
                                    arrpoints; cells=cells,
                                    cell_data=cell_data, optargs...)
        end

        return str

    else        # Case that wake is saved as horseshoes

        cells = Array{Int64,1}[]
        cell_data = []

        # Lines that make the wake
        for i in 0:nedges-1
            line = [i+2*nedges, i+0*nedges, i+1*nedges, i+3*nedges]
            push!(cells, line)
        end

        # Strength of wake along shedding edge
        strengths = [_get_wakestrength_Gamma(self, i) for i in 1:nedges]

        push!(cell_data, Dict( "field_name"=> "Gamma",
                                "field_type"=> "scalar",
                                "field_data"=> strengths))

        points = collect(eachcol(points))

        # Generate VTK
        return gt.generateVTK(filename*suffix, points; cells=cells,
                                cell_data=cell_data,
                                override_cell_type=4, optargs...)

    end
end
#### END OF LIFTING BODY  ######################################################
