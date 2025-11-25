#=##############################################################################
# DESCRIPTION
    Non-linear lifting line method

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################

import LaTeXStrings: @L_str

################################################################################
# LIFTING LINE STRUCT
################################################################################
struct LiftingLine{ R<:Number, 
                    S<:StripwiseElement, 
                    VectorType<:AbstractVector{R}, 
                    MatrixType<:AbstractMatrix{R}, 
                    TensorType<:AbstractArray{R, 3},
                    LI<:LinearIndices}

    # Internal properties
    grid::gt.Grid                               # Flat-geometry grid
    linearindices::LI                           # Linear indices of grid.nodes where linearindices[i, j] 
                                                # is TE (j==1) or LE (j==2) of the i-th row of nodes

    ypositions::Vector{Float64}                 # Non-dimensional y-position of nodes, 2*y/b

    nelements::Int                              # Number of stripwise elements
    elements::Vector{S}                         # Stripwise elements

    # Pre-allocated memory for solver
    aerocenters::VectorType                     # Aerodynamic center of each stripwise element

    horseshoes::TensorType                      # Horseshoes nodes, where horseshoes[i, j, n] 
                                                # is the i-th coordinate of the j-th node in 
                                                # the n-th horseshoe

    Dinfs::TensorType                           # Direction of each semi-infinite vortex filament
                                                # (freestream direction at each TE), where
                                                # Dinfs[i, j, n] is the i-th coordinate of the 
                                                # j-th semi-infinite filament (1==a, 2==b) of the
                                                # n-th horeseshoe

    midpoints::MatrixType                       # Midpoint along bound vortex for probing the velocity
    controlpoints::MatrixType                   # Control point of each horseshoe
    normals::MatrixType                         # Normal of each horseshoe

    aoas::VectorType                            # (deg) angle of attack seen by each stripwise element
    Gammas::VectorType                          # Circulation of each horseshoe

    G::MatrixType                               # Geometry matrix in linear system of equations
    RHS::VectorType                             # Right-hand-side vector in linear system of equations

    # Solver settings
    kerneloffset::Float64                       # Kernel offset to avoid singularities
    kernelcutoff::Float64                       # Kernel cutoff to avoid singularities


    function LiftingLine{R}(
                            airfoil_distribution, 
                            args...;
                            element_optargs=(), 
                            initial_aerodynamic_centers=1/4,
                            initial_controlpoint_positions=3/4,
                            initial_Uinf=[1, 0, 0],
                            kerneloffset=1e-8,
                            kernelcutoff=1e-14,
                            arraytype::Type=Array,
                            optargs...
                            ) where {R<:Number}

        # Define concrete array types
        VectorType = arraytype{R, 1}
        MatrixType = arraytype{R, 2}
        TensorType = arraytype{R, 3}

        # ------------------ DISCRETIZE WING -----------------------------------
        (; b, ypositions, chords, twists, 
        sweeps, dihedrals, spanaxiss, 
        symmetric) = _discretize_wing_parameterization(args...; optargs...)

        # ------------------ PREALLOCATE GRID MEMORY ---------------------------
        P_min = [0, 0, 0]               # Lower boundary span, chord, dummy
        P_max = [1, 1, 0]               # Upper boundary span, chord, dummy

        NDIVS = [length(ypositions)-1, 1, 0]   # Divisions of span, chord, and dummy collapsed dimension (flat surface)

        # Generate parametric grid
        grid = gt.Grid(P_min, P_max, NDIVS)

        # Linear indices of grid
        linearindices = LinearIndices(grid._ndivsnodes)
        
        # ------------------ MORPH GRID INTO WING GEOMETRY ---------------------
        _morph_grid_wing!(grid, b, ypositions, chords, twists, sweeps, dihedrals, 
                                            spanaxiss, symmetric, linearindices; 
                                            center=true)

        # ------------------ GENERATE STRIPWISE ELEMENTS -----------------------
        ypositions_elements = (ypositions[2:end] + ypositions[1:end-1]) / 2

        elements = _generate_stripwise_elements(airfoil_distribution, 
                                                ypositions_elements; 
                                                element_optargs...)
        nelements = length(elements)

        # ------------------ PRE-ALLOCATE SOLVER MEMORY ------------------------
        aerocenters = VectorType(undef, nelements)
        horseshoes = TensorType(undef, 3, 4, nelements)
        Dinfs = TensorType(undef, 3, 2, nelements)
        midpoints = MatrixType(undef, 3, nelements)
        controlpoints = MatrixType(undef, 3, nelements)
        normals = MatrixType(undef, 3, nelements)
        aoas = VectorType(undef, nelements)
        Gammas = VectorType(undef, nelements)
        G = MatrixType(undef, nelements, nelements)
        RHS = VectorType(undef, nelements)

        aoas .= 0
        Gammas .= 0
        G .= 0
        RHS .= 0

        # ------------------ INITIALIZE SOLVER SETTINGS ------------------------
        aerocenters .= initial_aerodynamic_centers

        calc_horseshoes!(horseshoes, grid.nodes, linearindices, nelements,
                                                        aerocenters)

        calc_midpoints!(midpoints, horseshoes, nelements)

        calc_controlpoints!(controlpoints, grid.nodes, linearindices, nelements,
                                                 initial_controlpoint_positions)

        calc_normals!(normals, controlpoints, horseshoes, nelements)

        calc_Dinfs!(Dinfs, initial_Uinf, nelements)


        new{R,
            eltype(elements),
            VectorType, MatrixType, TensorType, 
            typeof(linearindices)}(
                                grid, linearindices,
                                ypositions, 
                                nelements, elements,
                                aerocenters, horseshoes, Dinfs, 
                                midpoints, controlpoints, normals,
                                aoas, Gammas,
                                G, RHS,
                                kerneloffset, kernelcutoff)

    end

    LiftingLine(args...; optargs...) = LiftingLine{Float64}(args...; optargs...)

end

"""
Morph the lifting-line wing geometry into a new geometry
"""
function remorph!(self::LiftingLine, args...; 
                    recenter=false, 
                    aerodynamic_centers=1/4, controlpoint_positions=3/4, 
                    Uinf=[1, 0, 0],
                    optargs...)

    # Discretize parameterization
    (; b, ypositions, chords, twists, 
    sweeps, dihedrals, spanaxiss, symmetric) = _discretize_wing_parameterization(args...; optargs...)

    # Morph existing wing geometry into the new geometry
    _morph_grid_wing!(self.grid, b, ypositions, chords, twists, sweeps, dihedrals, 
                            spanaxiss, symmetric, self.linearindices; center=recenter)

    # Update horseshoe geometries
    self.aerocenters .= aerocenters

    calc_horseshoes!(self)

    calc_midpoints!(self)

    calc_controlpoints!(self, controlpoint_positions)

    calc_normals!(self)

    calc_Dinfs!(self, Uinf)

    return nothing
end


function save(self::LiftingLine, filename::AbstractString; 
                    format="vtk", 
                    horseshoe_suffix="_horseshoe",
                    midpoint_suffix="_midp", 
                    controlpoint_suffix="_cp", 
                    planar_suffix="_planar", 
                    wake_length_factor=0.25,
                    optargs...)

    str = ""

    # ------------- OUTPUT HORSESHOES ------------------------------------------

    npoints = size(self.horseshoes, 2)                  # Number of points in a horseshoe
    nhorseshoes = size(self.horseshoes, 3)

    # Flatten the tensor of horseshoe point into a matrix of points
    horseshoes = self.horseshoes
    horseshoes = reshape(horseshoes, ( size(horseshoes, 1), npoints*nhorseshoes ))

    # Connect the horseshoe points (0-indexed for VTK format)
    lines = [ collect( (0:npoints-1) .+ (ei-1)*npoints ) for ei in 1:nhorseshoes]

    # Determine a characteristic length for the wake
    wake_length = norm( maximum(horseshoes; dims=2) - minimum(horseshoes; dims=2) )
    wake_length *= wake_length_factor

    # Add fictitious wake end points to the horseshoes
    horseshoes = hcat(horseshoes, zeros(3, 2*nhorseshoes))

    for ei in 1:nhorseshoes

        # TE points
        Ap = horseshoes[:, lines[ei][1]+1]
        Bp = horseshoes[:, lines[ei][end]+1]

        # Indices of wake end points
        ai = npoints*nhorseshoes + 2*(ei-1) + 1
        bi = npoints*nhorseshoes + 2*(ei-1) + 2

        # Add wake end points
        horseshoes[:, ai] .= Ap + wake_length*self.Dinfs[:, 1, ei]
        horseshoes[:, bi] .= Bp + wake_length*self.Dinfs[:, 2, ei]

        # Add points to the VTK line definition
        lines[ei] = vcat(ai-1, lines[ei], bi-1)

    end

    # Format the points as a vector of vectors 
    horseshoes = eachcol(horseshoes)

    horseshoes_data = [
                            Dict(   "field_name"  => "Gamma",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.Gammas)

                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.aoas)
                        ]

    str *= gt.generateVTK(filename*horseshoe_suffix, horseshoes; cells=lines,
                                cell_data=horseshoes_data,
                                override_cell_type=4, optargs...)

    # ------------- OUTPUT MIDPOINTS -------------------------------------------
    midpoints = eachcol(self.midpoints)

    midpoints_data = [
                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.aoas)
                        ]

    str *= gt.generateVTK(filename*midpoint_suffix, midpoints; 
                                point_data=midpoints_data, optargs...)

    # ------------- OUTPUT CONTROL POINTS --------------------------------------
    controlpoints = eachcol(self.controlpoints)

    controlpoints_data = [
                            Dict(   "field_name"  => "normal",
                                    "field_type"  => "vector",
                                    "field_data"  => eachcol(self.normals)),

                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.aoas)
                        ]

    str *= gt.generateVTK(filename*controlpoint_suffix, controlpoints; 
                                point_data=controlpoints_data, optargs...)

    #  ------------- OUTPUT PLANAR GEOMETRY ------------------------------------
    planar_data = [
                            Dict(   "field_name"  => "Gamma",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.Gammas)

                            Dict(   "field_name"  => "angleofattack",
                                    "field_type"  => "scalar",
                                    "field_data"  => self.aoas)
                        ]

    str *= gt.save(self.grid, filename*planar_suffix; 
                                cell_data=planar_data, format, optargs...)

    return str
end


function calc_horseshoes!(self::LiftingLine, args...; optargs...) 
    return calc_horseshoes!(self.horseshoes, 
                                self.grid.nodes, self.linearindices, 
                                self.nelements, 
                                self.aerocenters,
                                args...; optargs...) 
end

function calc_horseshoes!(horseshoes::AbstractArray,
                            nodes::AbstractMatrix, linearindices::LinearIndices, 
                            nelements::Int,
                            aerocenters::AbstractVector
                            )
    
    for ei in 1:nelements               # Iterate over horseshoes
        for i in 1:3                    # Iterate over coordinates

            TEa = nodes[i, linearindices[ei, 1]]    # TE at a-side
            TEb = nodes[i, linearindices[ei+1, 1]]  # TE at b-side
            LEa = nodes[i, linearindices[ei, 2]]    # LE at a-side
            LEb = nodes[i, linearindices[ei+1, 2]]  # LE at b-side

                                                    # Aerodynamic center at a and b sides
                                                    # NOTE: this assumes contiguous horseshoes
            ACa = ei==1 ?         aerocenters[ei] : (aerocenters[ei-1] + aerocenters[ei]  )/2
            ACb = ei==nelements ? aerocenters[ei] : (aerocenters[ei]   + aerocenters[ei+1])/2

            # Horseshoe are defined from Ap -> A -> B -> Bp
            horseshoes[i, 1, ei] = TEa                      # Ap point
            horseshoes[i, 2, ei] = LEa + ACa*(TEa-LEa)      # A point
            horseshoes[i, 3, ei] = LEb + ACb*(TEb-LEb)      # B point
            horseshoes[i, 4, ei] = TEb                      # Bp point

        end
    end

end

function calc_midpoints!(self::LiftingLine, args...; optargs...) 
    return calc_midpoints!(self.midpoints, 
                                self.horseshoes, 
                                self.nelements, 
                                args...; optargs...) 
end

function calc_midpoints!(midpoints::AbstractMatrix,
                            horseshoes::AbstractArray,
                            nelements::Int;
                            position::Number = 0.5
                            )
    for ei in 1:nelements               # Iterate over horseshoes
        for i in 1:3                    # Iterate over coordinates

            midpoints[i, ei] = position * (horseshoes[i, 2, ei] + horseshoes[i, 3, ei])

        end
    end

end

function calc_controlpoints!(self::LiftingLine, args...; optargs...) 
    return calc_controlpoints!(self.controlpoints, 
                                self.grid.nodes, self.linearindices, 
                                self.nelements, 
                                args...; optargs...) 
end

function calc_controlpoints!(controlpoints::AbstractMatrix,
                            nodes::AbstractMatrix, linearindices::LinearIndices, 
                            nelements::Int,
                            position::Number
                            )
    for ei in 1:nelements               # Iterate over horseshoes
        for i in 1:3                    # Iterate over coordinates

            TEa = nodes[i, linearindices[ei, 1]]    # TE at a-side
            TEb = nodes[i, linearindices[ei+1, 1]]    # TE at b-side
            LEa = nodes[i, linearindices[ei, 2]]    # LE at a-side
            LEb = nodes[i, linearindices[ei+1, 2]]    # LE at b-side

            controlpoints[i, ei] = LEa + position*(TEa-LEa)     # Chordwise position at a-side
            controlpoints[i, ei] += LEb + position*(TEb-LEb)    # Chordwise position at b-side
            controlpoints[i, ei] /= 2                           # Take the average of the two

        end
    end

end


function calc_normals!(self::LiftingLine) 
    return calc_normals!(self.normals, self.controlpoints, self.horseshoes,
                                                                self.nelements)
end

function calc_normals!(normals::AbstractMatrix, 
                            controlpoints::AbstractMatrix, 
                            horseshoes::AbstractArray,
                            nelements::Int
                            )
    for ei in 1:nelements               # Iterate over horseshoes

        # dX = CP - A
        dX1 = controlpoints[1, ei] - horseshoes[1, 2, ei]
        dX2 = controlpoints[2, ei] - horseshoes[2, 2, ei]
        dX3 = controlpoints[3, ei] - horseshoes[3, 2, ei]

        # dY = B - A
        dY1 = horseshoes[1, 3, ei] - horseshoes[1, 2, ei]
        dY2 = horseshoes[2, 3, ei] - horseshoes[2, 2, ei]
        dY3 = horseshoes[3, 3, ei] - horseshoes[3, 2, ei]

        # normal = dX × dY / ||dX × dY||
        normals[1, ei] = dX2*dY3 - dX3*dY2
        normals[2, ei] = dX3*dY1 - dX1*dY3
        normals[3, ei] = dX1*dY2 - dX2*dY1
        normals[:, ei] /= sqrt(normals[1, ei]^2 + normals[2, ei]^2 + normals[3, ei]^2)

    end

end


function calc_Dinfs!(self::LiftingLine, Uinfs) 
    return calc_Dinfs!(self.Dinfs, Uinfs, self.nelements)
end

function calc_Dinfs!(Dinfs::AbstractArray, Uinf::AbstractVector, nelements::Int)
    return calc_Dinfs!(Dinfs, repeat(Uinf, 1, nelements), nelements)
end

function calc_Dinfs!(Dinfs::AbstractArray, Uinfs::AbstractMatrix, nelements::Int)

    for ei in 1:nelements               # Iterate over horseshoes

        for i in 1:3                    # Iterate over coordinates

                                                    # Freestream at a and b sides
                                                    # NOTE: this assumes contiguous horseshoes
            Uinfa = ei==1 ?         Uinfs[i, ei] : (Uinfs[i, ei-1] + Uinfs[i, ei]  )/2
            Uinfb = ei==nelements ? Uinfs[i, ei] : (Uinfs[i, ei]   + Uinfs[i, ei+1])/2

            # Semi-infite direction at Ap (non-unit vector)
            Dinfs[i, 1, ei] = Uinfa

            # Semi-infite direction at Bp (non-unit vector)
            Dinfs[i, 2, ei] = Uinfb

        end

        # Normalize the unit vectors
        Dinfs[:, 1, ei] /= sqrt(Dinfs[1, 1, ei]^2 + Dinfs[2, 1, ei]^2 + Dinfs[3, 1, ei]^2)
        Dinfs[:, 2, ei] /= sqrt(Dinfs[1, 2, ei]^2 + Dinfs[2, 2, ei]^2 + Dinfs[3, 2, ei]^2)

    end

end

"""
    Uind!(self::LiftingLine, targets, out, args...; optargs...)

Returns the velocity induced by the wing on the targets `targets`, which is a
3xn matrix. It adds the velocity at the i-th target to `out[:, i]`.
"""
function Uind!(self::LiftingLine, targets::AbstractMatrix, out::AbstractMatrix; optargs...)

    # Add bound vortex contributions
    for ei in 1:self.nelements                              # Iterates over horseshoes

        # Velocity of i-th horseshoe on every target
        U_vortexring(
                            view(self.horseshoes, :, :, ei),   # All nodes
                            1:4,                               # Indices of nodes that make this panel (closed ring)
                            self.Gammas[ei],                   # Horseshoe circulation
                            targets,                           # Targets
                            out;                               # Outputs
                            offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                            cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                            optargs...
                         )
    end

    # Add wake contributions
    TE = [1, size(self.horseshoes, 2)]                      # Indices of TE nodes in each horseshoe
    
    for ei in 1:self.nelements                              # Iterate over semi-infinite segments

        da1 = self.Dinfs[1, 1, ei]
        da2 = self.Dinfs[2, 1, ei]
        da3 = self.Dinfs[3, 1, ei]
        db1 = self.Dinfs[1, 2, ei]
        db2 = self.Dinfs[2, 2, ei]
        db3 = self.Dinfs[3, 2, ei]

        U_semiinfinite_horseshoe(
                          view(self.horseshoes, :, :, ei),                  # All nodes
                          TE,                                # Indices of nodes that make the shedding edge
                          da1, da2, da3,                     # Semi-infinite direction da
                          db1, db2, db3,                     # Semi-infinite direction db
                          self.Gammas[ei],                   # Filament circulation
                          targets,                           # Targets
                          out;                               # Outputs
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                          optargs...
                         )

    end
end





##### SOLVER ###################################################################

function solve_linear(self::LiftingLine, Uinf::AbstractVector, 
                                            args...; optargs...) 
    solve_linear(self, repeat(Uinf, 1, self.nelements), args...; optargs...)
end

function solve_linear(self::LiftingLine, Uinfs::AbstractMatrix;
                        controlpoint_positions=3/4, 
                        solver=solve_ludiv!, solver_optargs=(),
                        addfields=true, raise_warn=false,
                        optargs...)

    # Update semi-infinite wake to align with freestream
    calc_Dinfs!(self, Uinfs)

    # Update control point locations
    calc_controlpoints!(self, controlpoint_positions)

    # Update normals since control points where updated
    calc_normals!(self)

    # Compute geometric matrix (left-hand-side influence matrix)
    self.G .= 0
    _G_U!(self; optargs...)

    # Calculate boundary conditions (right-hand side of system of equations)
    self.RHS .= 0
    calc_bc_noflowthrough!(self.RHS, Uinfs, self.normals)

    # Solve system of equations
    solver(self.Gammas, self.G, self.RHS; solver_optargs...)

    if addfields
        gt.add_field(self.grid, "Uinf", "vector", collect(eachcol(Uinfs)), "cell"; raise_warn)
        gt.add_field(self.grid, "Gamma", "scalar", self.Gammas, "cell"; raise_warn)
    end

    return nothing

end

function _G_U!(self::LiftingLine; optargs...)

    # Build geometric matrix from panel contributions
    elements = 1:self.nelements
    chunks = collect(Iterators.partition(elements, max(length(elements) ÷ Threads.nthreads(), 3*Threads.nthreads())))

    Threads.@threads for chunk in chunks      # Distribute panel iteration among all CPU threads

        for ei in chunk                       # Iterate over elements

            U_vortexring(
                              view(self.horseshoes, :, :, ei),   # All nodes
                              1:4,                               # Indices of nodes that make this panel (closed ring)
                              1.0,                               # Unitary strength
                              self.controlpoints,                # Targets
                              view(self.G, :, ei);               # Velocity of ei-th panel on every CP
                              dot_with=self.normals,             # Normal of every CP
                              offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                              cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                              optargs...
                             )

         end
    end


    # Add wake contributions
    TE = [1, size(self.horseshoes, 2)]          # Indices of TE nodes in each horseshoe

    Threads.@threads for chunk in chunks        # Distribute wake iteration among all CPU threads

        for ei in chunk                          # Iterate horseshoes

            da1 = self.Dinfs[1, 1, ei]
            da2 = self.Dinfs[2, 1, ei]
            da3 = self.Dinfs[3, 1, ei]
            db1 = self.Dinfs[1, 2, ei]
            db2 = self.Dinfs[2, 2, ei]
            db3 = self.Dinfs[3, 2, ei]

            U_semiinfinite_horseshoe(
                              view(self.horseshoes, :, :, ei),   # All nodes
                              TE,                                # Indices of nodes that make the shedding edge
                              da1, da2, da3,                     # Semi-infinite direction da
                              db1, db2, db3,                     # Semi-infinite direction db
                              1.0,                               # Unitary strength
                              self.controlpoints,                # Targets
                              view(self.G, :, ei);                    # Velocity of wake panel on every CP
                              dot_with=self.normals,                  # Normal of every CP
                              offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                              cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                              optargs...
                             )

        end

    end
end



















##### INTERNAL FUNCTIONS #######################################################
"Discretize a parametric wing into segments"
function _discretize_wing_parameterization(;
                        # -------- Geometry parameters -------------------------
                        b::R0 = 1.0,                                            # (m) wing span
                        chord_distribution::Matrix{R1} = [0 1; 1 1],            # Chord distribution (nondim y-position 2*y/b, nondim chord c/b)
                        twist_distribution::Matrix{R2} = [0 0; 1 0],            # Twist distribution (nondim y-position 2*y/b, twist (deg))
                        sweep_distribution::Matrix{R3} = [0 0; 1 0],            # Sweep distribution (nondim y-position 2*y/b, sweep (deg))
                        dihedral_distribution::Matrix{R4} = [0 0; 1 0],         # Dihedral distribution (nondim y-position 2*y/b, dihedral (deg))
                        spanaxis_distribution::Matrix{R5} = [0 0.25; 1 0.25],   # Span-axis distribution: chordwise point about which the wing is twisted, swept, and dihedralized (nondim y-position 2*y/b, nondim chord-position x/c)

                        # -------- Discretization parameters -------------------
                        nelements::Int = 40,                                    # Number of stripwise elements per semi-span (ignored if `discretization` is provided)
                        discretization = [(1.0, nelements, 10.0, true)],        # Multi-discretization of wing (seg length, ndivs, expansion, central)
                        ypos_lo::Number = 0.0,                                  # Lower bound of non-dimensional span to discretize
                        ypos_up::Number = 1.0,                                  # Upper bound of non-dimensional span to discretize
                        interpolation::Function = math.linear,                  # Interpolation scheme. Example: `FLOWMath.linear` or `FLOWMath.akima`
                        symmetric::Bool = true,                                  # Whether the wing is symmetric
                        plot_discretization::Bool = true

                        ) where {R0, R1, R2, R3, R4, R5}

    R = promote_type(R0, R1, R2, R3, R4, R5)

    # ------------------ DISCRETIZE WING ---------------------------------------

    # Discretize the wing span
    ypositions = gt.multidiscretize(identity, ypos_lo, ypos_up, discretization)

    # Convert ypos from Vector{Any} to Vector{Float} for type stability
    ypositions = Float64.(ypositions)

    # Do an Akima spline through the input distributions and probe them at the discretized stations
    chords = interpolation(chord_distribution[:, 1], chord_distribution[:, 2], ypositions)
    twists = interpolation(twist_distribution[:, 1], twist_distribution[:, 2], ypositions)
    sweeps = interpolation(sweep_distribution[:, 1], sweep_distribution[:, 2], ypositions)
    dihedrals = interpolation(dihedral_distribution[:, 1], dihedral_distribution[:, 2], ypositions)
    spanaxiss = interpolation(spanaxis_distribution[:, 1], spanaxis_distribution[:, 2], ypositions)

    # Mirror the wing if symmetric
    if symmetric
        @assert minimum(ypositions) >= 0 ""*
            "Invalid input y-positions! They must be positive in a symmetric wing; got $(ypositions)"

        rng = ypositions[1]==0 ? (2:length(ypositions)) : (1:length(ypositions))

        ypositions = vcat(-reverse(ypositions[rng]), ypositions)
        chords = vcat(reverse(chords[rng]), chords)
        twists = vcat(reverse(twists[rng]), twists)
        sweeps = vcat(-reverse(sweeps[rng]), sweeps)
        dihedrals = vcat(-reverse(dihedrals[rng]), dihedrals)
        spanaxiss = vcat(reverse(spanaxiss[rng]), spanaxiss)
    end


    # Plot input and discretized distributions for verification
    if plot_discretization

        fig = plt.figure(figsize = [7*2, 0.5*5*3]*7/9)
        axs = fig.subplots(3, 2)
        axs = permutedims(axs, (2, 1))

        stl_inp = "o"
        fmt_inp = (color="black", label="Input")
        stl_out = "-s"
        fmt_out = (color="steelblue", markersize=3, alpha=0.8, linewidth=1, label="Output")

        for (axi, (inp, out, ylabel)) in enumerate([
                                                    (chord_distribution, chords, L"Chord $c/b$")
                                                    (twist_distribution, twists, L"Twist ($^\circ$)")
                                                    (sweep_distribution, sweeps, L"Sweep ($^\circ$)")
                                                    (dihedral_distribution, dihedrals, L"Dihedral ($^\circ$)")
                                                    (spanaxis_distribution, spanaxiss, L"Span axis position $x/c$")
                                            ])

            ax = axs[axi]

            ax.plot(inp[:, 1], inp[:, 2], stl_inp; fmt_inp...)
            ax.plot(ypositions, out, stl_out; fmt_out...)

            ax.set_ylabel(ylabel)

        end

        for ax in axs
            
            ax.set_xlabel(L"Span position ($2y/b$)")
            
            ax.spines["top"].set_visible(false)
            ax.spines["right"].set_visible(false)

            ax.legend(loc="best", frameon=false, fontsize=8)
        end
            
        fig.tight_layout()
    end

    return (; b, ypositions, chords, twists, sweeps, dihedrals, spanaxiss, symmetric)
end

function _morph_grid_wing!(grid, b, ypositions, chords, twists, sweeps, dihedrals, 
                            spanaxiss, symmetric, linearindices; center=false)

    @assert grid._ndivsnodes[1] == length(ypositions) ""*
        "Invalid grid! Received grid of $(grid._ndivsnodes[1]) spanwise nodes,"*
        " and $(ypositions) y-positions"

    # ------------------ MORPH GRID INTO WING GEOMETRY -------------------------

    ypos_prev = ypositions[1]
    x0_prev, y0_prev, z0_prev = 0, ypos_prev*(b/2), 0
    sweep_prev, dihedral_prev = 0, 0

    for (i, (ypos, cob, twist, sweep, dihedral, xoc)) in enumerate(zip(
                                                                ypositions, chords, twists, sweeps, dihedrals, spanaxiss
                                                            ))
        # Dimensionalize parameters
        dy = (ypos - ypos_prev) * (b/2)
        chord = cob*b
        
        # Define global position of span axis for this section
        x0 = x0_prev + dy*tand(symmetric && ypos<=0 ? sweep_prev : sweep)
        y0 = y0_prev + dy
        z0 = z0_prev + dy*tand(symmetric && ypos<=0 ? dihedral_prev : dihedral)

        # Define global coordinate of leading edge
        xLE = x0 - xoc*chord*cosd(twist)
        yLE = y0
        zLE = z0 + xoc*chord*sind(twist)

        # Define global coordinate of trailing edge
        xTE = x0 + (1-xoc)*chord*cosd(twist)
        yTE = y0
        zTE = z0 - (1-xoc)*chord*sind(twist)

        # Write leading edge point into the grid nodes
        iLE = linearindices[i, 2]
        grid.nodes[1, iLE] = xLE
        grid.nodes[2, iLE] = yLE
        grid.nodes[3, iLE] = zLE
        
        # Write trailing edge point into the grid nodes
        iTE = linearindices[i, 1]
        grid.nodes[1, iTE] = xTE
        grid.nodes[2, iTE] = yTE
        grid.nodes[3, iTE] = zTE

        ypos_prev = ypos
        x0_prev = x0
        y0_prev = y0
        z0_prev = z0
        sweep_prev = sweep
        dihedral_prev = dihedral

    end

    # Center the wing nose on the origin
    if center
        xorigin = minimum(view(grid.nodes, 1, :))
        yorigin = (minimum(view(grid.nodes, 2, :)) + maximum(view(grid.nodes, 2, :))) / 2
        zorigin = grid.nodes[3, findmin(view(grid.nodes, 1, :))[2]]

        grid.nodes[1, :] .-= xorigin
        grid.nodes[2, :] .-= yorigin
        grid.nodes[3, :] .-= zorigin
    end

    return nothing

end

function _generate_stripwise_elements(airfoil_distribution, ypositions; 
                                        extrapolated=true, plot_polars=true, 
                                        optargs...)

    # Create baseline stripwise elements from polars
    airfoils = _read_polars(airfoil_distribution; optargs...)

    # Identify StripwiseElement types
    element_types = unique([typeof(airfoil) for (ypos, airfoil) in airfoils])
    element_types = Union{element_types...}

    # Extrapolate polars to +-180 deg
    if extrapolated
        airfoils_extrapolated = [(ypos, extrapolate(airfoil)) for (ypos, airfoil) in airfoils]
    else
        airfoils_extrapolated = airfoils
    end

    # Blend the elements along the span
    lo_i = 1                                        # Index of lower-bound element
    airfoils_to_blend = airfoils_extrapolated
    airfoils_blended = []
    elements = element_types[]

    for ypos in ypositions

        # Find upper-bound element
        up_i = findfirst(x -> ypos < x[1], airfoils_to_blend)

        # Case that there is no upper bound: default to last element
        if isnothing(up_i)
            lo_i = length(airfoils_to_blend)
            up_i = lo_i
        end

        # Fix lower-bound element if it found a new upper bound
        if up_i - lo_i >= 2
            lo_i = up_i - 1
        end

        # Fetch bounds
        ypos_lo, airfoil_lo = airfoils_to_blend[lo_i]
        ypos_up, airfoil_up = airfoils_to_blend[up_i]

        # Catch case that airfoil bounds are different element types
        # TODO: Define logic for blending disimilar StripwiseElements
        @assert typeof(airfoil_lo) <: typeof(airfoil_up) ""*
            "Requested to blend dissimilar StripwiseElement types "*
            "$(typeof(airfoil_lo)) and $(typeof(airfoil_up))"

        # Determine blending weight
        weight = (ypos - ypos_lo) / (ypos_up - ypos_lo)

        # Catch case that bounds are the same
        if isinf(weight) || isnan(weight)
            weight = 0.0
        end

        # Blend elements
        blended_airfoil = blend(airfoil_lo, airfoil_up, weight)

        push!(elements, blended_airfoil)
        plot_polars && push!(airfoils_blended, (ypos, blended_airfoil))
        
    end

    # Plot polars for verification
    if plot_polars
        _plot_polars(airfoils, airfoils_extrapolated, airfoils_blended)
    end

    return elements

end

function _read_polars(airfoil_distribution; optargs...)

    return [(ypos, element(polar_file; optargs...)) for (ypos, polar_file, element) in airfoil_distribution]

end

function _plot_polars(airfoils, airfoils_extrapolated, airfoils_blended)

    stl_org = ""
    stl_extrap = "-"
    stl_blnd = "-"

    fmt_org = (; marker=".", alpha=0.5)
    fmt_extrap = (; linewidth=1, alpha=0.5)
    fmt_blnd = (; linewidth=1, alpha=0.25)

    # Compare raw vs extrapolated
    fig1 = plt.figure(figsize = [7, 0.75*5*3]*7/9 )
    axs1 = fig1.subplots(3, 1)

    fig = fig1
    axs = axs1
    fig.suptitle("Extrapolation comparison")

    for (ypos, airfoil) in airfoils
        
        clr = plt.cm.gnuplot(ypos)
        stl = stl_org
        fmt = fmt_org

        axs[1].plot(airfoil.alpha, airfoil.cl, stl; label=L"$2y/b = $"*"$(ypos)", color=clr, fmt...)
        axs[2].plot(airfoil.alpha, airfoil.cd, stl; label=L"$2y/b = $"*"$(ypos)", color=clr, fmt...)
        axs[3].plot(airfoil.alpha, airfoil.cm, stl; label=L"$2y/b = $"*"$(ypos)", color=clr, fmt...)

    end

    for (ypos, airfoil) in airfoils_extrapolated
        
        clr = plt.cm.gnuplot(ypos)
        stl = stl_extrap
        fmt = fmt_extrap

        axs[1].plot(airfoil.alpha, airfoil.cl, stl; color=clr, fmt...)
        axs[2].plot(airfoil.alpha, airfoil.cd, stl; color=clr, fmt...)
        axs[3].plot(airfoil.alpha, airfoil.cm, stl; color=clr, fmt...)

    end

    ax = axs[1]
    ax.set_ylabel(L"Lift $c_\ell$")

    ax = axs[2]
    ax.set_ylabel(L"Drag $c_d$")

    ax = axs[3]
    ax.set_ylabel(L"Pitching moment $c_m$")

    for ax in axs
        ax.set_xlabel(L"Angle of attack ($^\circ$)")
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        ax.legend(loc="best", frameon=false, fontsize=8)
    end
        
    fig.tight_layout()

    # Compare blends
    fig2 = plt.figure(figsize = [7, 0.75*5*3]*7/9 )
    axs2 = fig2.subplots(3, 1)

    fig = fig2
    axs = axs2
    fig.suptitle("Blending comparison")

    for (ypos, airfoil) in airfoils_blended
        
        clr = plt.cm.gnuplot(ypos)
        stl = stl_blnd
        fmt = fmt_blnd

        axs[1].plot(airfoil.alpha, airfoil.cl, stl; label=L"$2y/b = $"*"$(ypos)", color=clr, fmt...)
        axs[2].plot(airfoil.alpha, airfoil.cd, stl; label=L"$2y/b = $"*"$(ypos)", color=clr, fmt...)
        axs[3].plot(airfoil.alpha, airfoil.cm, stl; label=L"$2y/b = $"*"$(ypos)", color=clr, fmt...)

    end

    ax = axs[1]
    ax.set_ylabel(L"Lift $c_\ell$")

    ax = axs[2]
    ax.set_ylabel(L"Drag $c_d$")

    ax = axs[3]
    ax.set_ylabel(L"Pitching moment $c_m$")

    for ax in axs
        ax.set_xlabel(L"Angle of attack ($^\circ$)")
        ax.spines["top"].set_visible(false)
        ax.spines["right"].set_visible(false)
        # ax.legend(loc="best", frameon=false, fontsize=8)
    end
        
    fig.tight_layout()

    return (fig1, axs1), (fig2, axs2)
end
##### END OF LIFTING LINE ######################################################