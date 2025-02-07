#=##############################################################################
# DESCRIPTION
    Definition of 2D body types

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Jan 2025
  * License     : MIT License
=###############################################################################

################################################################################
# 2D BODY TYPE
################################################################################
"""
Numerical observations:
* LinearSource2D bodies seem to be numerically ill conditioned on closed bodies.
    This is probably related to the fact that a linear source actually resembles
    a point vortex, so closing the body leads to the reduction in degrees of
    freedom that we see in a watertight body with point vortex elements.
    To avoid that, make sure that the body is open and use `check_mesh=false`.

* ConstantSource2D seems to be have a similar ill conditioning for some reason.
    I have observed that leaving the TE open or making the TE sharp and singular
    (meaning, the TE points very close to each other) lead to good results.

"""
struct Body2D{E<:AbstractElement2D, N} <: AbstractBody{E, N}

    # User inputs
    grid::gt.GridSurface2D                       # Paneled geometry

    # Properties
    nnodes::Int                               # Number of nodes
    ncells::Int                               # Number of cells
    fields::Array{String, 1}                  # Available fields (solutions)
    Oaxis::Array{<:Number,2}                  # Coordinate system of original grid
    O::Array{<:Number,1}                      # Position of CS of original grid

    # Internal variables
    strength::Matrix                          # strength[i,j] is the stength of the i-th panel with the j-th element type
    CPoffset::Float64                         # Control point offset in normal direction
    kerneloffset::Float64                     # Kernel offset to avoid singularities
    characteristiclength::Function            # Characteristic length of each panel


    function Body2D{E, N}(
                            grid;
                            nnodes=grid.nnodes, ncells=grid.ncells,
                            fields=Array{String,1}(),
                            Oaxis=Array(1.0I, 3, 3), O=zeros(3),
                            strength=zeros(grid.ncells, N),
                            CPoffset=1e-14,
                            kerneloffset=1e-8,
                            characteristiclength=characteristiclength_unitary,
                            check_mesh=true
                            ) where {E, N}

        if check_mesh
            check2D(E, N, grid)
        end

        return new(grid,
                    nnodes, ncells,
                    fields,
                    Oaxis, O,
                    strength,
                    CPoffset, kerneloffset, characteristiclength)
    end

    function Body2D{E}(args...; optargs...) where E
        return Body2D{E, multiplicity(E)}(args...; optargs...)
    end

    function Body2D{E}(args...; optargs...) where E <: Union{LinearSource2D, UniformVortex2D}
        return Body2D{E, 2}(args...; optargs...)
    end

end


function check2D(E::Type, N::Int, grid; debug=false)

    # Check that grid connectivity is contiguous
    if E <: Union{LinearSource2D, UniformVortex2D} && N==2

        prev_panel = nothing
        for pi in 1:grid.ncells

            panel = gt.get_cell(grid, pi)

            if pi != 1
                @assert panel[1] == prev_panel[2] ""*
                    "Body of element type $(E) is expected to have a contiguous"*
                    " connectivity; it found a discontinuity between panels "*
                    "$(pi-1) and $(pi) ($(prev_panel) and $(panel))"
            end

            prev_panel = panel
        end

        if debug; println("Passed!"); end;

    end

    # Check that grid is closed
    if E <: Union{LinearSource2D, UniformVortex2D} && N==2

        @assert gt.get_cell(grid, 1)[1] == gt.get_cell(grid, grid.ncells)[2] ""*
            "Body of element type $(E) is expected to be a closed contour"

        if debug; println("Passed!"); end;

    end

end


save(body::Body2D, args...; optargs...) = save_base(body, args...; optargs...)




################################################################################
# ORIGINAL HESS-SMITH SOLVER (CONSTANT SOURCE + UNIFORM VORTEX + KUTTA CONDITION)
################################################################################
# See https://github.com/EdoAlvarezR/MyPanel2D/blob/master/docs/MyPanel2D_documentation.ipynb

HessSmithBodies = Union{
                            Body2D{Union{ConstantSource2D, UniformVortex2D}, 2},
                            Body2D{Union{LinearSource2D, UniformVortex2D}, 2}
                        }

_get_Gdims(body::HessSmithBodies) = (body.ncells+1, body.ncells+1)

function solve(body::HessSmithBodies,
                Uinfs::AbstractMatrix{T};
                solver=solve_ludiv!, solver_optargs=(),
                kuttapanels=(1, body.ncells),
                GPUArray=Array{T},
                optargs...
                ) where {T<:Number}

    @assert size(Uinfs) == (3, body.ncells) ""*
        "Invalid Uinfs; expected size (3, $(body.ncells)), got $(size(Uinfs))"

    # Compute normals and control points
    normals = _calc_normals(body)
    tangents = _calc_tangents(body)
    CPs = _calc_controlpoints(body, normals)

    # Compute geometric matrix (left-hand-side influence matrix)
    G = zeros(T, _get_Gdims(body)...)
    _G_U!(body, G, CPs, normals, tangents, kuttapanels; optargs...)

    # Calculate boundary conditions (right-hand side of system of equations)
    RHS = calc_bc_noflowthrough(Uinfs, normals)

    # Add row of Kutta condition to RHS
    # NOTE: Here we make the RHS negative so we don't need to flip the sign of
    #       the LHS row
    push!(RHS, 0)
    for kpi in kuttapanels
        RHS[end] -= dot(tangents[:, kpi], Uinfs[:, kpi])
    end

    # Solve system of equations
    strengths = GPUArray(undef, body.ncells+1)
    solver(strengths, G, RHS; solver_optargs...)

    # Port solution back to CPU if solved in GPU
    if !(GPUArray <: Array)
        Gamma = Array{T}(Gamma)
    end

    # Save solution
    set_solution(body, strengths, Uinfs)

end


function set_solution(body::HessSmithBodies, strengths, Uinfs)

    @. body.strength[:, 1] = strengths[1:end-1]
    @. body.strength[:, 2] = strengths[end]

    _solvedflag(body, true)
    add_field(body, "Uinf", "vector", collect(eachcol(Uinfs)), "cell")
    add_field(body, "sigma", "scalar", view(body.strength, :, 1), "cell")
    add_field(body, "gamma", "scalar", view(body.strength, :, 2), "cell")
end


function _G_U!(body::Body2D{Union{ConstantSource2D, UniformVortex2D}, 2},
                G::AbstractMatrix, CPs::AbstractMatrix,
                normals::AbstractMatrix, tangents::AbstractMatrix,
                kuttapanels::Tuple{Int, Int}; optargs...
                )

    N = body.ncells+1
    M = size(CPs, 2)+1

    @assert size(G, 1)==M && size(G, 2)==N ""*
        "Matrix G with invalid dimensions; got $(size(G)), expected $((M, N))."

    @assert size(normals, 2)==size(CPs, 2) ""*
        "normals matrix with invalid dimensions; got $(size(normals)), expected $((2, M))."

    # Build geometric matrix from panel contributions
    panels = 1:body.ncells
    # chunks = collect(Iterators.partition(panels, max(length(panels) ÷ Threads.nthreads(), 3*Threads.nthreads())))

    # Threads.@threads for chunk in chunks      # Distribute panel iteration among all CPU threads

        # Pre-allocate memory for panel calculation
        out, coor, lin, ndivscells, cin = gt.generate_getcellt_args!(body.grid)

        # Impose no-flow-throug condition
        # for pj in chunk                       # Iterate over panels
        for pj in panels

            coor[1] = pj
            panel = gt.get_cell_t!(out, body.grid, coor, lin, ndivscells)

            # Add contribution of source panel on every target
            U_2D_constant_source(
                                body.grid._nodes,               # All nodes
                                panel,                          # Indices of nodes that make this panel
                                1.0,                            # Unitary strength
                                CPs,                            # Targets
                                view(G, 1:M-1, pj);             # Add velocity of j-th panel on every CP
                                dot_with=normals,               # Normal of every CP
                                offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                optargs...
                                )

            # Add contribution of vortex panel on every target
            U_2D_constant_vortex(
                                body.grid._nodes,               # All nodes
                                panel,                          # Indices of nodes that make this panel
                                1.0,                            # Unitary strength
                                CPs,                            # Targets
                                view(G, 1:M-1, N);              # Add velocity of j-th panel on every CP
                                dot_with=normals,               # Normal of every CP
                                offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                optargs...
                                )
         end

         # Impose Kutta condition
         # for pj in chunk                       # Iterate over panels
         for pj in panels

             coor[1] = pj
             panel = gt.get_cell_t!(out, body.grid, coor, lin, ndivscells)

             for kuttapanel in kuttapanels     # Iterate over Kutta panels

                 kCP = view(CPs, :, kuttapanel:kuttapanel)
                 ktangent = view(tangents, :, kuttapanel:kuttapanel)

                 # Add contribution of source panel on Kutta panel
                 U_2D_constant_source(
                                     body.grid._nodes,               # All nodes
                                     panel,                          # Indices of nodes that make this panel
                                     1.0,                            # Unitary strength
                                     kCP,                            # Kutta panel
                                     view(G, M, pj);                 # Add velocity of j-th panel on Kutta panel
                                     dot_with=ktangent,              # Tangent of Kutta panel
                                     offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                     optargs...
                                     )

                 # Add contribution of vortex panel on Kutta panel
                 U_2D_constant_vortex(
                                     body.grid._nodes,               # All nodes
                                     panel,                          # Indices of nodes that make this panel
                                     1.0,                            # Unitary strength
                                     kCP,                            # Kutta panel
                                     view(G, M, N);                  # Add velocity of j-th panel on Kutta panel
                                     dot_with=ktangent,              # Tangent of Kutta panel
                                     offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                     optargs...
                                     )
             end
          end

          # # Add self-induced component of surface velocity (-∇μ/2 = -γ/2) to
          # # the Kutta condition
          # for kuttapanel in kuttapanels     # Iterate over Kutta panels
          #     G[M, N] += -1/2
          #     # G[M, N] -= 1
          # end

    # end

end

function _Uind!(body::Body2D{Union{ConstantSource2D, UniformVortex2D}, 2},
                                                    targets, out; optargs...)

    # Pre-allocate memory for panel calculation
    pout, coor, lin, ndivscells, cin = gt.generate_getcellt_args!(body.grid)

    # Iterates over panels
    for i in 1:body.ncells

        coor[1] = i
        panel = gt.get_cell_t!(pout, body.grid, coor, lin, ndivscells)

        # Velocity of i-th panel on every target
        U_2D_constant_source(
                            body.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            body.strength[i, 1],               # Source strength
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=body.kerneloffset,          # Offset of kernel to avoid singularities
                            optargs...
                         )

         U_2D_constant_vortex(
                             body.grid._nodes,                  # All nodes
                             panel,                             # Indices of nodes that make this panel
                             body.strength[i, 2],               # Vortex strength
                             targets,                           # Targets
                             out;                               # Outputs
                             offset=body.kerneloffset,          # Offset of kernel to avoid singularities
                             optargs...
                          )
    end
end



################################################################################
# MODIFIED HESS-SMITH SOLVER (LINEAR SOURCE + UNIFORM VORTEX + KUTTA CONDITION)
################################################################################
function _G_U!(body::Body2D{Union{LinearSource2D, UniformVortex2D}, 2},
                G::AbstractMatrix, CPs::AbstractMatrix,
                normals::AbstractMatrix, tangents::AbstractMatrix,
                kuttapanels::Tuple{Int, Int}; optargs...
                )

    N = body.ncells+1
    M = size(CPs, 2)+1

    @assert size(G, 1)==M && size(G, 2)==N ""*
        "Matrix G with invalid dimensions; got $(size(G)), expected $((M, N))."

    @assert size(normals, 2)==size(CPs, 2) ""*
        "normals matrix with invalid dimensions; got $(size(normals)), expected $((2, M))."

    # Build geometric matrix from panel contributions
    panels = 1:body.ncells
    # chunks = collect(Iterators.partition(panels, max(length(panels) ÷ Threads.nthreads(), 3*Threads.nthreads())))

    # Threads.@threads for chunk in chunks      # Distribute panel iteration among all CPU threads

        # Pre-allocate memory for panel calculation
        out, coor, lin, ndivscells, cin = gt.generate_getcellt_args!(body.grid)

        # Impose no-flow-throug condition
        # for pj in chunk                       # Iterate over panels
        for pj in panels

            coor[1] = pj
            panel = gt.get_cell_t!(out, body.grid, coor, lin, ndivscells)

            # Add contribution of first source node on every target
            U_2D_linear_source(
                                body.grid._nodes,               # All nodes
                                panel,                          # Indices of nodes that make this panel
                                1.0,                            # Unitary strength for first node
                                0.0,                            # Omit contribution of second node
                                CPs,                            # Targets
                                view(G, 1:M-1, pj);             # Add velocity of first source node of j-th panel on every CP
                                dot_with=normals,               # Normal of every CP
                                offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                optargs...
                                )

            # Add contribution of second source node on every target
            U_2D_linear_source(
                                body.grid._nodes,               # All nodes
                                panel,                          # Indices of nodes that make this panel
                                0.0,                            # Omit contribution of first node
                                1.0,                            # Unitary strength for second node
                                CPs,                            # Targets
                                view(G, 1:M-1, pj==body.ncells ? 1 : pj+1); # Add velocity of second source node of j-th panel on every CP
                                dot_with=normals,               # Normal of every CP
                                offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                optargs...
                                )

            # Add contribution of vortex panel on every target
            U_2D_constant_vortex(
                                body.grid._nodes,               # All nodes
                                panel,                          # Indices of nodes that make this panel
                                1.0,                            # Unitary strength
                                CPs,                            # Targets
                                view(G, 1:M-1, N);              # Add velocity of j-th panel on every CP
                                dot_with=normals,               # Normal of every CP
                                offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                optargs...
                                )
         end

         # Impose Kutta condition
         # for pj in chunk                       # Iterate over panels
         for pj in panels

             coor[1] = pj
             panel = gt.get_cell_t!(out, body.grid, coor, lin, ndivscells)

             for kuttapanel in kuttapanels     # Iterate over Kutta panels

                 kCP = view(CPs, :, kuttapanel:kuttapanel)
                 ktangent = view(tangents, :, kuttapanel:kuttapanel)

                 # Add contribution of first source node on Kutta panel
                 U_2D_linear_source(
                                     body.grid._nodes,               # All nodes
                                     panel,                          # Indices of nodes that make this panel
                                     1.0,                            # Unitary strength for first node
                                     0.0,                            # Omit contribution of second node
                                     kCP,                            # Kutta panel
                                     view(G, M, pj);                 # Add velocity of first source node of j-th panel on Kutta panel
                                     dot_with=ktangent,              # Tangent of Kutta panel
                                     offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                     optargs...
                                     )

                 # Add contribution of second source node on Kutta panel
                 U_2D_linear_source(
                                     body.grid._nodes,               # All nodes
                                     panel,                          # Indices of nodes that make this panel
                                     0.0,                            # Omit contribution of first node
                                     1.0,                            # Unitary strength for second node
                                     kCP,                            # Kutta panel
                                     view(G, M, pj==body.ncells ? 1 : pj+1); # Add velocity of first source node of j-th panel on Kutta panel
                                     dot_with=ktangent,              # Tangent of Kutta panel
                                     offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                     optargs...
                                     )

                 # Add contribution of vortex panel on Kutta panel
                 U_2D_constant_vortex(
                                     body.grid._nodes,               # All nodes
                                     panel,                          # Indices of nodes that make this panel
                                     1.0,                            # Unitary strength
                                     kCP,                            # Kutta panel
                                     view(G, M, N);                  # Add velocity of j-th panel on Kutta panel
                                     dot_with=ktangent,              # Tangent of Kutta panel
                                     offset=body.kerneloffset,       # Offset of kernel to avoid singularities
                                     optargs...
                                     )
             end
          end

    # end

end

function _Uind!(body::Body2D{Union{LinearSource2D, UniformVortex2D}, 2},
                                                    targets, out; optargs...)

    # Pre-allocate memory for panel calculation
    pout, coor, lin, ndivscells, cin = gt.generate_getcellt_args!(body.grid)

    # Iterates over panels
    for i in 1:body.ncells

        coor[1] = i
        panel = gt.get_cell_t!(pout, body.grid, coor, lin, ndivscells)

        # Velocity of i-th panel on every target
        U_2D_linear_source(
                            body.grid._nodes,                  # All nodes
                            panel,                             # Indices of nodes that make this panel
                            body.strength[i, 1],               # Strength at first node
                            body.strength[i==body.ncells ? 1 : i+1, 1], # Strength at second node
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=body.kerneloffset,          # Offset of kernel to avoid singularities
                            optargs...
                         )

         U_2D_constant_vortex(
                             body.grid._nodes,                  # All nodes
                             panel,                             # Indices of nodes that make this panel
                             body.strength[i, 2],               # Vortex strength
                             targets,                           # Targets
                             out;                               # Outputs
                             offset=body.kerneloffset,          # Offset of kernel to avoid singularities
                             optargs...
                          )
    end
end

#### COMMON FUNCTIONS ###########################################################
function _phi!(body::Body2D, args...; optargs...)
    nothing
end
#### END OF 2D BODY  ###########################################################
