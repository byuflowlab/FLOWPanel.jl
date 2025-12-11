#=##############################################################################
# DESCRIPTION
    Non-linear lifting line method

# AUTHORSHIP
  * Created by  : Eduardo J. Alvarez
  * Email       : Edo.AlvarezR@gmail.com
  * Date        : Nov 2025
  * License     : MIT License
=###############################################################################


################################################################################
# LIFTING LINE STRUCT
################################################################################
struct LiftingLine{ R<:Number, 
                    S<:StripwiseElement, N,
                    VectorType<:AbstractVector{R}, 
                    MatrixType<:AbstractMatrix{R}, 
                    TensorType<:AbstractArray{R, 3},
                    TensorType2<:AbstractArray{R, 4},
                    LI<:LinearIndices} <: AbstractBody{S, N}

    # Internal properties
    grid::gt.Grid                               # Flat-geometry grid
    linearindices::LI                           # Linear indices of grid.nodes where linearindices[i, j] 
                                                # is TE (j==1) or LE (j==2) of the i-th row of nodes
    fields::Vector{String}                      # Available fields (solutions)

    ypositions::Vector{Float64}                 # Non-dimensional y-position of nodes, 2*y/b

    nelements::Int                              # Number of stripwise elements
    elements::Vector{S}                         # Stripwise elements

    deltasb::Float64                            # Blending distance, deltasb = 2*dy/b
    deltajoint::Float64                         # Joint distance, deltajoint = dx/c
    sigmafactor::Float64                        # Dragging line amplification factor
    sigmaexponent::Float64                      # Dragging line amplification exponent

    # Pre-allocated memory for solver
    aerocenters::VectorType                     # Aerodynamic center of each stripwise element
    strippositions::VectorType                  # Position of each stripwise element within the bound vortex (0==a, 1==b)

    horseshoes::TensorType                      # Horseshoes nodes, where horseshoes[i, j, n] 
                                                # is the i-th coordinate of the j-th node in 
                                                # the n-th horseshoe

    effective_horseshoes::TensorType2           # Horseshoes nodes of the effective lift line,
                                                # where effective_horseshoes[:, :, :, ei] is the
                                                # effective horseshoes seen by the ei-th stripwise element

    Dinfs::TensorType                           # Direction of each semi-infinite vortex filament
                                                # (freestream direction at each TE), where
                                                # Dinfs[i, j, n] is the i-th coordinate of the 
                                                # j-th semi-infinite filament (1==a, 2==b) of the
                                                # n-th horeseshoe

    midpoints::MatrixType                       # Midpoint along lifting line for probing the velocity
    controlpoints::MatrixType                   # Control point of each horseshoe
    tangents::MatrixType                        # Tangent of each horseshoe (direction of streamwise element). Orthogonal bases: span x normal = tangent
    spans::MatrixType                           # Spanwise unit vector of each horseshoe (direction of span at each streamwise element). Orthogonal bases: normal x tangent = span
    normals::MatrixType                         # Normal of each horseshoe (normal to streamwise element). Orthogonal bases: tangent x span = normal

    swepttangents::MatrixType                   # Unit vector in the direction of the effective swept chord (u_aΛ in Goates 2022 notation). line x sweptnormal = swepttangent
    lines::MatrixType                           # Unit vector in the direction of the lifting line (u_sΛ in Goates 2022 notation). sweptnormal x swepttangent = line
    sweptnormals::MatrixType                    # Unit vector in the normal to the lifting line and effective swep chord (-u_nΛ in Goates 2022 notation). swepttangent x line = sweptnormal

    auxtangents::MatrixType                     # Auxiliary memory for calculating tangents of effective lifting line

    aoas::VectorType                            # (deg) swept angle of attack seen by each stripwise element (𝛼_Λ in Goates 2022 notation)
    claeros::VectorType                         # Purely-aerodynamic sectional lift coefficient of each stripwise element (used for calculating Gamma)
    Gammas::VectorType                          # Circulation of each horseshoe (lifting line)
    sigmas::VectorType                          # Source strength of each horseshoe (dragging line)
    Us::MatrixType                              # Velocity at each midpoint
    chords::VectorType                          # Dimensional chord of each stripwise element

    G::MatrixType                               # Geometry matrix in linear system of equations (matrix in left-hand side)
    RHS::VectorType                             # Right-hand-side vector in linear system of equations

    residuals::VectorType                       # Non-linear solver residuals
    Geff::TensorType                            # Precomputed geometric matrix for evaluating the self-induced velocity by the effective horseshoes on each midpoint. 
                                                # Geff[mi, ei, i] is the i-th coordinate of the unitary-strength velocity induced by the ei-th effective horseshoe on the mi-th midpoint.
    # Solver settings
    kerneloffset::Float64                       # Kernel offset to avoid singularities
    kernelcutoff::Float64                       # Kernel cutoff to avoid singularities


    function LiftingLine{R}(
                            airfoil_distribution, 
                            args...;
                            element_optargs=(), 
                            aerodynamic_centers=1/4,
                            strip_positions=0.5,
                            controlpoint_position=3/4,
                            deltasb=1.0,
                            deltajoint=0.15,
                            sigmafactor=1.0,
                            sigmaexponent=1.0,
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
        TensorType2 = arraytype{R, 4}

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
        strippositions = VectorType(undef, nelements)
        horseshoes = TensorType(undef, 3, 4, nelements)
        effective_horseshoes = TensorType2(undef, 3, 4, nelements, nelements)
        Dinfs = TensorType(undef, 3, 2, nelements)

        midpoints = MatrixType(undef, 3, nelements)
        controlpoints = MatrixType(undef, 3, nelements)

        tangents = MatrixType(undef, 3, nelements)
        spans = MatrixType(undef, 3, nelements)
        normals = MatrixType(undef, 3, nelements)

        swepttangents = MatrixType(undef, 3, nelements)
        lines = MatrixType(undef, 3, nelements)
        sweptnormals = MatrixType(undef, 3, nelements)

        auxtangents = MatrixType(undef, 3, nelements)

        aoas = VectorType(undef, nelements)
        claeros = VectorType(undef, nelements)
        Gammas = VectorType(undef, nelements)
        Gammas = VectorType(undef, nelements)
        sigmas = VectorType(undef, nelements)
        Us = MatrixType(undef, 3, nelements)
        chords = VectorType(undef, nelements)
        G = MatrixType(undef, nelements, nelements)
        RHS = VectorType(undef, nelements)
        residuals = VectorType(undef, nelements)
        Geff = TensorType(undef, nelements, nelements, 3)

        aoas .= 0
        claeros .= 0
        Gammas .= 0
        sigmas .= 0
        Us .= 0
        G .= 0
        RHS .= 0
        residuals .= 0
        Geff .= 0

        # ------------------ INITIALIZE SOLVER SETTINGS ------------------------
        aerocenters .= aerodynamic_centers
        strippositions .= strip_positions

        calc_horseshoes!(horseshoes, grid.nodes, linearindices, nelements,
                                                        aerocenters)

        calc_midpoints!(midpoints, horseshoes, strippositions, nelements)

        calc_chords!(chords, grid.nodes, linearindices, strippositions, nelements)

        calc_controlpoints!(controlpoints, grid.nodes, linearindices, 
                                                strippositions, nelements,
                                                controlpoint_position)

        calc_tangents!(tangents, horseshoes, strippositions, nelements)
        calc_normals!(normals, tangents, horseshoes, nelements)
        calc_spans!(spans, normals, tangents, nelements)

        calc_lines!(lines, horseshoes, nelements)
        calc_swepttangents!(swepttangents, lines, tangents, nelements)
        calc_sweptnormals!(sweptnormals, swepttangents, lines, nelements)

        calc_Dinfs!(Dinfs, initial_Uinf, nelements)

        calc_effective_horseshoes!(effective_horseshoes, horseshoes, midpoints,
                                                tangents, spans, 
                                                ypositions, strippositions, 
                                                nelements; deltasb)

        jointerize!(effective_horseshoes, auxtangents,
                                                strippositions, nelements,
                                                normals, chords; deltajoint)

        calc_Geff!(Geff, effective_horseshoes, Dinfs,
                    midpoints,
                    swepttangents, lines, sweptnormals,
                    nelements;
                    offset=kerneloffset, cutoff=kernelcutoff)

        S = eltype(elements)
        new{R,
            S, _count(S),
            VectorType, MatrixType, TensorType, TensorType2,
            typeof(linearindices)}(
                                grid, linearindices, String[],
                                ypositions, 
                                nelements, elements,
                                deltasb, deltajoint, sigmafactor, sigmaexponent,
                                aerocenters, strippositions,
                                horseshoes, effective_horseshoes, Dinfs, 
                                midpoints, controlpoints, 
                                tangents, spans, normals,
                                swepttangents, lines, sweptnormals,
                                auxtangents,
                                aoas, claeros, Gammas, sigmas, Us, chords,
                                G, RHS, 
                                residuals, Geff,
                                kerneloffset, kernelcutoff)

    end

    LiftingLine(args...; optargs...) = LiftingLine{Float64}(args...; optargs...)

end

"""
Morph the lifting-line wing geometry into a new geometry
"""
function remorph!(self::LiftingLine, args...; 
                    recenter=false, 
                    aerodynamic_centers=1/4,
                    controlpoint_positions=3/4, 
                    Uinf=[1, 0, 0],
                    deltasb=self.deltasb,
                    deltajoint=self.deltajoint,
                    optargs...)

    # Discretize parameterization
    (; b, ypositions, chords, twists, 
    sweeps, dihedrals, spanaxiss, symmetric) = _discretize_wing_parameterization(args...; optargs...)

    # Morph existing wing geometry into the new geometry
    _morph_grid_wing!(self.grid, b, ypositions, chords, twists, sweeps, dihedrals, 
                            spanaxiss, symmetric, self.linearindices; center=recenter)

    # Update horseshoe geometries
    self.aerocenters .= aerodynamic_centers

    calc_horseshoes!(self)

    calc_midpoints!(self)
    calc_chords!(self)
    calc_controlpoints!(self, controlpoint_positions)

    calc_tangents!(self)
    calc_normals!(self)
    calc_spans!(self)

    calc_lines!(self)
    calc_swepttangents!(self)
    calc_sweptnormals!(self)

    calc_Dinfs!(self, Uinf)
    
    calc_effective_horseshoes!(self; deltasb)

    jointerize!(self; deltajoint)

    calc_Geff!(self)

    # Reset solution
    self.aoas .= 0
    self.claeros .= 0
    self.Gammas .= 0
    self.sigmas .= 0
    self.Us .= 0

    return nothing
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


function calc_effective_horseshoes!(self::LiftingLine, args...; optargs...) 
    return calc_effective_horseshoes!(self.effective_horseshoes,
                                            self.horseshoes, self.midpoints,
                                            self.tangents, self.spans,
                                            self.ypositions, 
                                            self.strippositions,
                                            self.nelements, args...; 
                                            deltasb=self.deltasb, optargs...) 
end

"""
Calculate Reid's effective lifting line of horseshoes as explained in Goates
2022, JoA, "Modern Implementation and Evaluation of Lifting-Line Theory
for Complex Geometries".
"""
function calc_effective_horseshoes!(effective_horseshoes::AbstractArray{<:Number, 4}, 
                                            horseshoes::AbstractArray{<:Number, 3},
                                            midpoints::AbstractMatrix,
                                            tangents::AbstractMatrix,
                                            spans::AbstractMatrix,
                                            ypositions::AbstractVector,
                                            strippositions::AbstractVector,
                                            nelements::Int;
                                            deltasb=1.0     # Blending distance
                                            )
    
    for ei in 1:nelements               # Iterate over stripwise elements

        # Lifting filament direction
        dl1 = horseshoes[1, 3, ei] - horseshoes[1, 2, ei]
        dl2 = horseshoes[2, 3, ei] - horseshoes[2, 2, ei]
        dl3 = horseshoes[3, 3, ei] - horseshoes[3, 2, ei]
        magdl = sqrt(dl1^2 + dl2^2 + dl3^2)

        dl1 /= magdl
        dl2 /= magdl
        dl3 /= magdl

        # TE direction
        dTE1 = horseshoes[1, 4, ei] - horseshoes[1, 1, ei]
        dTE2 = horseshoes[2, 4, ei] - horseshoes[2, 1, ei]
        dTE3 = horseshoes[3, 4, ei] - horseshoes[3, 1, ei]
        magdTE = sqrt(dTE1^2 + dTE2^2 + dTE3^2)

        dTE1 /= magdTE
        dTE2 /= magdTE
        dTE3 /= magdTE

        # TE midpoint
        midTE1 = horseshoes[1, 1, ei] + strippositions[ei]*(horseshoes[1, 4, ei] - horseshoes[1, 1, ei])
        midTE2 = horseshoes[2, 1, ei] + strippositions[ei]*(horseshoes[2, 4, ei] - horseshoes[2, 1, ei])
        midTE3 = horseshoes[3, 1, ei] + strippositions[ei]*(horseshoes[3, 4, ei] - horseshoes[3, 1, ei])

        # Estimate semi-span length to dimensionalize ypos
        semispan = magdl*abs( dl1*spans[1, ei] + dl2*spans[2, ei] + dl3*spans[3, ei] ) / (ypositions[ei+1] - ypositions[ei])

        sweep = calc_sweep(horseshoes, tangents, spans, ei)     # Sweep of this (original) horseshoe. NOTE: should this use ei or ni?
        sigma = (2 / (deltasb*cosd(sweep)) )^2   # Blending parameter

        # Non-dimensional y-position of stripwise element
        ypos_e = ypositions[ei] + strippositions[ei]*(ypositions[ei+1] - ypositions[ei])

        for ni in 1:nelements           # Iterate over effective horseshoes

            # Non-dimensional y-position of a and b sides
            ypos_a = ypositions[ni]
            ypos_b = ypositions[ni+1]
            w_a = exp(-sigma * (ypos_a - ypos_e)^2) # Blending weight on a-side
            w_b = exp(-sigma * (ypos_b - ypos_e)^2) # Blending weight on b-side

            for (i, dl, dTE, midTE) in zip(1:3, (dl1, dl2, dl3), (dTE1, dTE2, dTE3), (midTE1, midTE2, midTE3))
                
                # Effective Ap
                effective_horseshoes[i, 1, ni, ei] = (1-w_a)*horseshoes[i, 1, ni] + w_a*(midTE + dTE*semispan*(ypositions[ni] - ypos_e)/cosd(sweep))

                # Effective A
                effective_horseshoes[i, 2, ni, ei] = (1-w_a)*horseshoes[i, 2, ni] + w_a*(midpoints[i, ei] + dl*semispan*(ypositions[ni] - ypos_e)/cosd(sweep))

                # Effective B
                effective_horseshoes[i, 3, ni, ei] = (1-w_b)*horseshoes[i, 3, ni] + w_b*(midpoints[i, ei] + dl*semispan*(ypositions[ni+1] - ypos_e)/cosd(sweep))

                # Effective Bp
                effective_horseshoes[i, 4, ni, ei] = (1-w_b)*horseshoes[i, 4, ni] + w_b*(midTE + dTE*semispan*(ypositions[ni+1] - ypos_e)/cosd(sweep))

            end

        end
    end

end


"""
Convert Weissenger VLM horseshoes into Reid's joint horseshoes as explained in
Reid (2022), "A General Approach to Lifting-Line Theory, Applied to Wings With
Sweep", Sec. 2.2.3. 
"""
function jointerize!(self::LiftingLine; optargs...)

    jointerize!(self.effective_horseshoes, self.auxtangents, 
                    self.strippositions, self.nelements,
                    self.normals, self.chords; 
                    deltajoint=self.deltajoint, 
                    optargs...)
end

function jointerize!(effective_horseshoes::AbstractArray{<:Number, 4}, 
                        auxtangents::AbstractMatrix, 
                        strippositions::AbstractVector, nelements::Int,
                        args...; optargs...)

    for ei in 1:nelements                   # Iterate over stripwise elements

        # Fetch effective horseshoes seen by this element
        horseshoes = view(effective_horseshoes, :, :, :, ei)

        # Calculate tangent vectors of the effective lifting line
        tangents = auxtangents
        calc_tangents!(tangents, horseshoes, strippositions, nelements)

        # Modify the effective horseshoes to be joint
        jointerize!(horseshoes, tangents, args..., nelements; optargs...)
    end

end

function jointerize!(horseshoes::AbstractArray{R, 3}, tangents::AbstractMatrix, 
                        normals::AbstractMatrix, chords::AbstractVector,
                        nelements::Int;
                        deltajoint=0.15) where {R}

    if deltajoint < 0
        return
    end

    prev_tangent1::R = zero(R)
    prev_tangent2::R = zero(R)
    prev_tangent3::R = zero(R)

    this_tangent1::R = zero(R)
    this_tangent2::R = zero(R)
    this_tangent3::R = zero(R)

    prev_chord::R = zero(R)
    this_chord::R = zero(R)

    for ei in 0:nelements                   # Iterate over horseshoes

        if ei != nelements

            # Lifting filament direction of the next horseshoe
            dl1 = horseshoes[1, 3, ei+1] - horseshoes[1, 2, ei+1]
            dl2 = horseshoes[2, 3, ei+1] - horseshoes[2, 2, ei+1]
            dl3 = horseshoes[3, 3, ei+1] - horseshoes[3, 2, ei+1]
            magdl = sqrt(dl1^2 + dl2^2 + dl3^2)

            dl1 /= magdl
            dl2 /= magdl
            dl3 /= magdl

            # Calculate the tangent vector that is orthonormal to the filament
            next_tangent1 = dl2*normals[3, ei+1] - dl3*normals[2, ei+1]
            next_tangent2 = dl3*normals[1, ei+1] - dl1*normals[3, ei+1]
            next_tangent3 = dl1*normals[2, ei+1] - dl2*normals[1, ei+1]
            magtangent = sqrt(next_tangent1^2 + next_tangent2^2 + next_tangent3^2)

            next_tangent1 /= magtangent
            next_tangent2 /= magtangent
            next_tangent3 /= magtangent

            next_chord = chords[ei+1]

        else
            next_tangent1 = this_tangent1
            next_tangent2 = this_tangent2
            next_tangent3 = this_tangent3

            next_chord = this_chord
        end

        # Proceed to jointerize this horseshoe if it is not the initialization step
        if ei != 0

            # Tangent and chord on a-side
            tangenta1 = (prev_tangent1 + this_tangent1)/2
            tangenta2 = (prev_tangent2 + this_tangent2)/2
            tangenta3 = (prev_tangent3 + this_tangent3)/2
            chorda = (prev_chord + this_chord)/2

            # Tangent and chord on b-side
            tangentb1 = (this_tangent1 + next_tangent1)/2
            tangentb2 = (this_tangent2 + next_tangent2)/2
            tangentb3 = (this_tangent3 + next_tangent3)/2
            chordb = (this_chord + next_chord)/2

            # Override original Ap and Bp with joint Ap and joint Bp
            horseshoes[1, 1, ei] = horseshoes[1, 2, ei] + (deltajoint*chorda)*tangenta1
            horseshoes[2, 1, ei] = horseshoes[2, 2, ei] + (deltajoint*chorda)*tangenta2
            horseshoes[3, 1, ei] = horseshoes[3, 2, ei] + (deltajoint*chorda)*tangenta3

            horseshoes[1, 4, ei] = horseshoes[1, 3, ei] + (deltajoint*chordb)*tangentb1
            horseshoes[2, 4, ei] = horseshoes[2, 3, ei] + (deltajoint*chordb)*tangentb2
            horseshoes[3, 4, ei] = horseshoes[3, 3, ei] + (deltajoint*chordb)*tangentb3

        end

        if ei == 0
            this_tangent1 = next_tangent1
            this_tangent2 = next_tangent2
            this_tangent3 = next_tangent3

            this_chord = next_chord
        end

        # Shift tangents
        prev_tangent1 = this_tangent1
        prev_tangent2 = this_tangent2
        prev_tangent3 = this_tangent3

        this_tangent1 = next_tangent1
        this_tangent2 = next_tangent2
        this_tangent3 = next_tangent3

        prev_chord = this_chord
        this_chord = next_chord
        
    end
end

"Rotate joint segments to align with freestream rather than surface"
function align_joints_with_Uinfs!(ll::LiftingLine, Uinfs::AbstractMatrix)

    for ei in 1:ll.nelements                   # Iterate over stripwise elements

        # Fetch effective horseshoes seen by this element
        horseshoes = view(ll.effective_horseshoes, :, :, :, ei)

        # Modify the effective horseshoes to aling joints with freestream
        align_joints_with_Uinfs!(horseshoes, ll.chords, ll.lines, Uinfs, ll.nelements, ll.deltajoint)
    end

end

function align_joints_with_Uinfs!(horseshoes::AbstractArray{R, 3}, 
                                    chords::AbstractVector,
                                    lines::AbstractMatrix,
                                    Uinfs::AbstractMatrix, 
                                    nelements::Int,
                                    deltajoint::Number) where {R}

    if deltajoint < 0
        return
    end

    prev_DinfΛ1::R = zero(R)
    prev_DinfΛ2::R = zero(R)
    prev_DinfΛ3::R = zero(R)

    this_DinfΛ1::R = zero(R)
    this_DinfΛ2::R = zero(R)
    this_DinfΛ3::R = zero(R)

    prev_chord::R = zero(R)
    this_chord::R = zero(R)

    for ei in 0:nelements                   # Iterate over horseshoes

        if ei != nelements
            
            # Lifting filament direction
            dl1 = horseshoes[1, 3, ei+1] - horseshoes[1, 2, ei+1]
            dl2 = horseshoes[2, 3, ei+1] - horseshoes[2, 2, ei+1]
            dl3 = horseshoes[3, 3, ei+1] - horseshoes[3, 2, ei+1]
            magdl = sqrt(dl1^2 + dl2^2 + dl3^2)

            lines1 = dl1/magdl
            lines2 = dl2/magdl
            lines3 = dl3/magdl

            # Project the velocity onto the filament direction
            UsΛ = Uinfs[1, ei+1]*lines1

            # Substract filament-component from the velocity
            UinfΛ1 = Uinfs[1, ei+1] - UsΛ*lines1
            UinfΛ2 = Uinfs[2, ei+1] - UsΛ*lines2
            UinfΛ3 = Uinfs[3, ei+1] - UsΛ*lines3
            magUinfΛ = sqrt(UinfΛ1^2 + UinfΛ2^2 + UinfΛ3^2)

            # Counter-projected freestream direction
            next_DinfΛ1 = UinfΛ1 / magUinfΛ
            next_DinfΛ2 = UinfΛ2 / magUinfΛ
            next_DinfΛ3 = UinfΛ3 / magUinfΛ

            next_chord = chords[ei+1]

        else
            next_DinfΛ1 = this_DinfΛ1
            next_DinfΛ2 = this_DinfΛ2
            next_DinfΛ3 = this_DinfΛ3

            next_chord = this_chord
        end

        # Proceed to jointerize this horseshoe if it is not the initialization step
        if ei != 0

            # Tangent and chord on a-side
            DinfΛa1 = (prev_DinfΛ1 + this_DinfΛ1)/2
            DinfΛa2 = (prev_DinfΛ2 + this_DinfΛ2)/2
            DinfΛa3 = (prev_DinfΛ3 + this_DinfΛ3)/2
            chorda = (prev_chord + this_chord)/2

            # Tangent and chord on b-side
            DinfΛb1 = (this_DinfΛ1 + next_DinfΛ1)/2
            DinfΛb2 = (this_DinfΛ2 + next_DinfΛ2)/2
            DinfΛb3 = (this_DinfΛ3 + next_DinfΛ3)/2
            chordb = (this_chord + next_chord)/2

            # Override original Ap and Bp with joint Ap and joint Bp
            horseshoes[1, 1, ei] = horseshoes[1, 2, ei] + (deltajoint*chorda)*DinfΛa1
            horseshoes[2, 1, ei] = horseshoes[2, 2, ei] + (deltajoint*chorda)*DinfΛa2
            horseshoes[3, 1, ei] = horseshoes[3, 2, ei] + (deltajoint*chorda)*DinfΛa3

            horseshoes[1, 4, ei] = horseshoes[1, 3, ei] + (deltajoint*chordb)*DinfΛb1
            horseshoes[2, 4, ei] = horseshoes[2, 3, ei] + (deltajoint*chordb)*DinfΛb2
            horseshoes[3, 4, ei] = horseshoes[3, 3, ei] + (deltajoint*chordb)*DinfΛb3

        end

        if ei == 0
            this_DinfΛ1 = next_DinfΛ1
            this_DinfΛ2 = next_DinfΛ2
            this_DinfΛ3 = next_DinfΛ3

            this_chord = next_chord
        end

        # Shift DinfΛs
        prev_DinfΛ1 = this_DinfΛ1
        prev_DinfΛ2 = this_DinfΛ2
        prev_DinfΛ3 = this_DinfΛ3

        this_DinfΛ1 = next_DinfΛ1
        this_DinfΛ2 = next_DinfΛ2
        this_DinfΛ3 = next_DinfΛ3

        prev_chord = this_chord
        this_chord = next_chord
        
    end

end

# function align_joints_with_Uinfs!(ll::LiftingLine, Uinfs::AbstractMatrix)

#     for ei in 1:ll.nelements                   # Iterate over stripwise elements

#         # Fetch effective horseshoes seen by this element
#         horseshoes = view(ll.effective_horseshoes, :, :, :, ei)

#         # Calculate tangent vectors of the effective lifting line
#         tangents = ll.auxtangents
#         calc_tangents!(tangents, horseshoes, ll.strippositions, ll.nelements)

#         # Modify the effective horseshoes to aling joints with freestream
#         align_joints_with_Uinfs!(horseshoes, ll.sweptnormals, ll.swepttangents, ll.lines, Uinfs, ll.nelements)
#     end

# end

# function align_joints_with_Uinfs!(horseshoes::AbstractArray{R, 3}, 
#                                     normals::AbstractMatrix,
#                                     tangents::AbstractMatrix,
#                                     lines::AbstractMatrix,
#                                     Uinfs::AbstractMatrix, 
#                                     nelements::Int) where {R}

#     M = zeros(R, 3, 3)

#     for ei in 1:nelements

#         # Calculate freestream angle relative to the surface
#         aoa_uinf = calc_aoa(Uinfs, normals, tangents, ei)

#         for ji in (1, 4)
#             # Calculate current joint angle relative to the surface
#             aoa_joint = calc_aoa(horseshoes[1, ji, ei], horseshoes[2, ji, ei], horseshoes[3, ji, ei], 
#                                                         normals, tangents, ei)
            
#             # Rotation matrix
#             gt.axis_rotation!(M, view(lines, :, ei), aoa_joint - aoa_uinf)

#             # Rotate joint node
#             X1 = horseshoes[1, ji, ei]
#             X2 = horseshoes[2, ji, ei]
#             X3 = horseshoes[3, ji, ei]

#             X01 = horseshoes[1, ji==1 ? 2 : 3, ei]
#             X02 = horseshoes[2, ji==1 ? 2 : 3, ei]
#             X03 = horseshoes[3, ji==1 ? 2 : 3, ei]

#             new_X1 = M[1, 1]*(X1-X01) + M[1, 2]*(X2-X02) + M[1, 3]*(X3-X03) + X01
#             new_X2 = M[2, 1]*(X1-X01) + M[2, 2]*(X2-X02) + M[2, 3]*(X3-X03) + X02
#             new_X3 = M[3, 1]*(X1-X01) + M[3, 2]*(X2-X02) + M[3, 3]*(X3-X03) + X03

#             horseshoes[1, ji, ei] = new_X1
#             horseshoes[2, ji, ei] = new_X2
#             horseshoes[3, ji, ei] = new_X3
#         end
        
#     end

#     # Take the average of overlaping joints
#     for ei in 1:nelements-1
#         for i in 1:3
            
#             new_X = (horseshoes[i, 4, ei] + horseshoes[i, 1, ei+1])/2

#             horseshoes[i, 4, ei] = new_X
#             horseshoes[i, 1, ei+1] = new_X
#         end
#     end

# end


# "Rotate joint segments to align with freestream rather than surface"
# function align_joints_with_Uinfs!(ll::LiftingLine, Uinfs::AbstractMatrix)

#     for ei in 1:ll.nelements                   # Iterate over stripwise elements

#         # Fetch effective horseshoes seen by this element
#         horseshoes = view(ll.effective_horseshoes, :, :, :, ei)

#         # Calculate tangent vectors of the effective lifting line
#         tangents = ll.auxtangents
#         calc_tangents!(tangents, horseshoes, ll.strippositions, ll.nelements)

#         # Modify the effective horseshoes to aling joints with freestream
#         align_joints_with_Uinfs!(horseshoes, ll.normals, ll.tangents, ll.lines, Uinfs, ll.nelements)
#     end

# end

# function align_joints_with_Uinfs!(horseshoes::AbstractArray{R, 3}, 
#                                     normals::AbstractMatrix,
#                                     tangents::AbstractMatrix,
#                                     lines::AbstractMatrix,
#                                     Uinfs::AbstractMatrix, 
#                                     nelements::Int) where {R}

#     M = zeros(R, 3, 3)

#     for ei in 1:nelements

#         # Calculate freestream angle relative to the surface
#         aoa_uinf = calc_aoa(Uinfs, normals, tangents, ei)

#         for ji in (1, 4)
#             # Calculate current joint angle relative to the surface
#             aoa_joint = calc_aoa(horseshoes[1, ji, ei], horseshoes[2, ji, ei], horseshoes[3, ji, ei], 
#                                                         normals, tangents, ei)
            
#             # Rotation matrix
#             gt.axis_rotation!(M, view(lines, :, ei), aoa_joint - aoa_uinf)

#             # Rotate joint node
#             X1 = horseshoes[1, ji, ei]
#             X2 = horseshoes[2, ji, ei]
#             X3 = horseshoes[3, ji, ei]

#             X01 = horseshoes[1, ji==1 ? 2 : 3, ei]
#             X02 = horseshoes[2, ji==1 ? 2 : 3, ei]
#             X03 = horseshoes[3, ji==1 ? 2 : 3, ei]

#             new_X1 = M[1, 1]*(X1-X01) + M[1, 2]*(X2-X02) + M[1, 3]*(X3-X03) + X01
#             new_X2 = M[2, 1]*(X1-X01) + M[2, 2]*(X2-X02) + M[2, 3]*(X3-X03) + X02
#             new_X3 = M[3, 1]*(X1-X01) + M[3, 2]*(X2-X02) + M[3, 3]*(X3-X03) + X03

#             horseshoes[1, ji, ei] = new_X1
#             horseshoes[2, ji, ei] = new_X2
#             horseshoes[3, ji, ei] = new_X3
#         end
        
#     end

#     # Take the average of overlaping joints
#     for ei in 1:nelements-1
#         for i in 1:3
            
#             new_X = (horseshoes[i, 4, ei] + horseshoes[i, 1, ei+1])/2

#             horseshoes[i, 4, ei] = new_X
#             horseshoes[i, 1, ei+1] = new_X
#         end
#     end

# end

function calc_midpoints!(self::LiftingLine, args...; optargs...) 
    return calc_midpoints!(self.midpoints, 
                                self.horseshoes, 
                                self.strippositions,
                                self.nelements, 
                                args...; optargs...) 
end

function calc_midpoints!(midpoints::AbstractMatrix,
                            horseshoes::AbstractArray,
                            positions::AbstractVector,
                            nelements::Int,
                            )
    for ei in 1:nelements               # Iterate over horseshoes
        for i in 1:3                    # Iterate over coordinates

            midpoints[i, ei] = horseshoes[i, 2, ei] + positions[ei]*(horseshoes[i, 3, ei] - horseshoes[i, 2, ei])

        end
    end

end

function calc_chords!(self::LiftingLine, args...; optargs...) 
    return calc_chords!(self.chords, 
                                self.grid.nodes, self.linearindices,
                                self.strippositions,
                                self.nelements, 
                                args...; optargs...) 
end

function calc_chords!(chords::AbstractVector,
                            nodes::AbstractMatrix, linearindices::LinearIndices, 
                            positions::AbstractVector,
                            nelements::Int
                            )

    for ei in 1:nelements               # Iterate over horseshoes

        TEai = linearindices[ei, 1]             # Index of TE at a-side
        TEbi = linearindices[ei+1, 1]           # Index of TE at b-side
        LEai = linearindices[ei, 2]             # Index of LE at a-side
        LEbi = linearindices[ei+1, 2]           # Index of LE at b-side

        chorda = sqrt( (nodes[1, TEai]-nodes[1, LEai])^2 + (nodes[2, TEai]-nodes[2, LEai])^2 + (nodes[3, TEai]-nodes[3, LEai])^2 )
        chordb = sqrt( (nodes[1, TEbi]-nodes[1, LEbi])^2 + (nodes[2, TEbi]-nodes[2, LEbi])^2 + (nodes[3, TEbi]-nodes[3, LEbi])^2 )

        chord = chorda + positions[ei]*(chordb - chorda)

        chords[ei] = chord

    end

end

function calc_controlpoints!(self::LiftingLine, args...; optargs...) 
    return calc_controlpoints!(self.controlpoints, 
                                self.grid.nodes, self.linearindices, 
                                self.strippositions,
                                self.nelements, 
                                args...; optargs...) 
end

function calc_controlpoints!(controlpoints::AbstractMatrix,
                            nodes::AbstractMatrix, linearindices::LinearIndices, 
                            positionsy::AbstractVector,
                            nelements::Int,
                            positionx::Number
                            )
    for ei in 1:nelements               # Iterate over horseshoes
        for i in 1:3                    # Iterate over coordinates

            TEa = nodes[i, linearindices[ei, 1]]        # TE at a-side
            TEb = nodes[i, linearindices[ei+1, 1]]      # TE at b-side
            LEa = nodes[i, linearindices[ei, 2]]        # LE at a-side
            LEb = nodes[i, linearindices[ei+1, 2]]      # LE at b-side

            posx_a = LEa + positionx*(TEa-LEa)           # Chordwise position at a-side
            posx_b = LEb + positionx*(TEb-LEb)           # Chordwise position at b-side

            controlpoints[i, ei] = posx_a + positionsy[ei]*(posx_b - posx_a) # Blend a and b sides

        end
    end

end

function calc_tangents!(self::LiftingLine) 
    return calc_tangents!(self.tangents, self.horseshoes, 
                                    self.strippositions, self.nelements)
end

function calc_tangents!(tangents::AbstractMatrix, horseshoes::AbstractArray, 
                            positions::AbstractVector, nelements::Int)
    for ei in 1:nelements               # Iterate over horseshoes

        # dA = Ap - A
        dA1 = horseshoes[1, 1, ei] - horseshoes[1, 2, ei]
        dA2 = horseshoes[2, 1, ei] - horseshoes[2, 2, ei]
        dA3 = horseshoes[3, 1, ei] - horseshoes[3, 2, ei]

        # dB = Bp - B
        dB1 = horseshoes[1, 4, ei] - horseshoes[1, 3, ei]
        dB2 = horseshoes[2, 4, ei] - horseshoes[2, 3, ei]
        dB3 = horseshoes[3, 4, ei] - horseshoes[3, 3, ei]

        # tangent = normalized( dA + position*(dB - dA) )
        tangents[1, ei] = dA1 + positions[ei]*(dB1 - dA1)
        tangents[2, ei] = dA2 + positions[ei]*(dB2 - dA2)
        tangents[3, ei] = dA3 + positions[ei]*(dB3 - dA3)
        tangents[:, ei] /= sqrt(tangents[1, ei]^2 + tangents[2, ei]^2 + tangents[3, ei]^2)

    end

end

function calc_normals!(self::LiftingLine) 
    return calc_normals!(self.normals, self.tangents, self.horseshoes,
                                                                self.nelements)
end

function calc_normals!(normals::AbstractMatrix, 
                            tangents::AbstractMatrix, 
                            horseshoes::AbstractArray,
                            nelements::Int
                            )
    for ei in 1:nelements               # Iterate over horseshoes

        # # dX = CP - A
        # dX1 = controlpoints[1, ei] - horseshoes[1, 2, ei]
        # dX2 = controlpoints[2, ei] - horseshoes[2, 2, ei]
        # dX3 = controlpoints[3, ei] - horseshoes[3, 2, ei]

        # dX = tangent
        dX1 = tangents[1, ei]
        dX2 = tangents[2, ei]
        dX3 = tangents[3, ei]

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


function calc_spans!(self::LiftingLine) 
    return calc_spans!(self.spans, self.normals, self.tangents, self.nelements)
end

function calc_spans!(spans::AbstractMatrix,
                            normals::AbstractMatrix,
                            tangents::AbstractMatrix,
                            nelements::Int
                            )
    for ei in 1:nelements               # Iterate over horseshoes

        # span = normal x tangent
        spans[1, ei] = normals[2, ei]*tangents[3, ei] - normals[3, ei]*tangents[2, ei]
        spans[2, ei] = normals[3, ei]*tangents[1, ei] - normals[1, ei]*tangents[3, ei]
        spans[3, ei] = normals[1, ei]*tangents[2, ei] - normals[2, ei]*tangents[1, ei]
        spans[:, ei] /= sqrt(spans[1, ei]^2 + spans[2, ei]^2 + spans[3, ei]^2)

    end

end


function calc_lines!(self::LiftingLine)
    calc_lines!(self.lines, self.horseshoes, self.nelements)
end

function calc_lines!(lines::AbstractMatrix, 
                        horseshoes::AbstractArray, nelements::Int)
    for ei in 1:nelements

        # Lifting filament direction
        dl1 = horseshoes[1, 3, ei] - horseshoes[1, 2, ei]
        dl2 = horseshoes[2, 3, ei] - horseshoes[2, 2, ei]
        dl3 = horseshoes[3, 3, ei] - horseshoes[3, 2, ei]
        magdl = sqrt(dl1^2 + dl2^2 + dl3^2)

        lines[1, ei] = dl1/magdl
        lines[2, ei] = dl2/magdl
        lines[3, ei] = dl3/magdl

    end
end


function calc_swepttangents!(self::LiftingLine)
    calc_swepttangents!(self.swepttangents, self.lines, self.tangents, self.nelements)
end

function calc_swepttangents!(swepttangents::AbstractMatrix, 
                                lines::AbstractMatrix, tangents::AbstractMatrix, 
                                nelements::Int)
    for ei in 1:nelements

        # Project tangent onto lifting line
        tline = tangents[1, ei]*lines[1, ei] + tangents[2, ei]*lines[2, ei] + tangents[3, ei]*lines[3, ei]

        # Substract projected component
        t1 = tangents[1, ei] - tline*lines[1, ei]
        t2 = tangents[2, ei] - tline*lines[2, ei]
        t3 = tangents[3, ei] - tline*lines[3, ei]
        magt = sqrt(t1^2 + t2^2 + t3^2)

        # Save unit vector
        swepttangents[1, ei] = t1/magt
        swepttangents[2, ei] = t2/magt
        swepttangents[3, ei] = t3/magt

    end                   
end


function calc_sweptnormals!(self::LiftingLine)
    calc_sweptnormals(self.sweptnormals, self.swepttangents, self.lines, self.nelements)
end

function calc_sweptnormals!(sweptnormals::AbstractMatrix, 
                                swepttangents::AbstractMatrix, 
                                lines::AbstractMatrix,
                                nelements::Int)
    for ei in 1:nelements

        # n = t x l
        n1 = swepttangents[2, ei]*lines[3, ei] - swepttangents[3, ei]*lines[2, ei]
        n2 = swepttangents[3, ei]*lines[1, ei] - swepttangents[1, ei]*lines[3, ei]
        n3 = swepttangents[1, ei]*lines[2, ei] - swepttangents[2, ei]*lines[1, ei]
        magn = sqrt(n1^2 + n2^2 + n3^2)

        # Save unit vector
        sweptnormals[1, ei] = t=n1/magn
        sweptnormals[2, ei] = t=n2/magn
        sweptnormals[3, ei] = t=n3/magn

    end                   
end

function calc_sweep(self::LiftingLine, args...)
    return calc_sweep(self.horseshoes, self.tangents, self.spans, args...)
end

function calc_sweep(horseshoes::AbstractArray, 
                    tangents::AbstractMatrix, spans::AbstractMatrix, ei::Int)
    
    # Lifting filament
    dl1 = horseshoes[1, 3, ei] - horseshoes[1, 2, ei]
    dl2 = horseshoes[2, 3, ei] - horseshoes[2, 2, ei]
    dl3 = horseshoes[3, 3, ei] - horseshoes[3, 2, ei]

    # Projections on to span and tangent directions
    dy = dl1*spans[1, ei] + dl2*spans[2, ei] + dl3*spans[3, ei]
    dx = dl1*tangents[1, ei] + dl2*tangents[2, ei] + dl3*tangents[3, ei]

    return atand(dx, dy)
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
Precompute geometric matrix for evaluating the self-induced velocity by the 
effective horseshoes on each midpoint. 

Geff[mi, ei, i] is the i-th coordinate of the unitary-strength velocity induced 
by the ei-th effective horseshoe on the mi-th midpoint.

NOTE: this precomputation includes both the surface and wake induced velocity.
"""
function calc_Geff!(self::LiftingLine; optargs...)
    calc_Geff!(self.Geff, self.effective_horseshoes, 
                self.Dinfs,
                self.midpoints,
                self.swepttangents, self.lines, self.sweptnormals, 
                self.nelements; 
                offset=self.kerneloffset, 
                cutoff=self.kernelcutoff,
                optargs...)
end

function calc_Geff!(Geff::AbstractArray{<:Number, 3}, 
                    effective_horseshoes::AbstractArray{<:Number, 4}, 
                    Dinfs::AbstractArray{<:Number, 3},
                    midpoints::AbstractMatrix,
                    tangents::AbstractMatrix, 
                    spans::AbstractMatrix, 
                    normals::AbstractMatrix,
                    nelements::Int;
                    optargs...)

    # Velocity projected on the tangent unit vector
    calc_Geff!(view(Geff, :, :, 1), effective_horseshoes, Dinfs, 
                            midpoints, tangents, nelements; optargs...)

    # Velocity projected on the spanwise unit vector
    calc_Geff!(view(Geff, :, :, 2), effective_horseshoes, Dinfs, 
                            midpoints, spans, nelements; optargs...)

    # Velocity projected on the normal unit vector
    calc_Geff!(view(Geff, :, :, 3), effective_horseshoes, Dinfs, 
                            midpoints, normals, nelements; optargs...)

end

function calc_Geff!(Geff::AbstractMatrix, 
                    effective_horseshoes::AbstractArray{<:Number, 4}, 
                    Dinfs::AbstractArray{<:Number, 3},
                    midpoints::AbstractMatrix,
                    dot_with::AbstractMatrix, 
                    nelements::Int;
                    optargs...)

    Geff .= 0                                   # Erase previous values

    TE = [1, size(effective_horseshoes, 2)]     # Indices of TE nodes in each horseshoe

    # Build geometric matrix from panel contributions
    Threads.@threads for mi in 1:nelements                         # Iterate over midpoints (iterations distributed among all CPU threads)

        # Fetch this midpoints
        targets = view(midpoints, :, mi:mi)

        # Fetch effective horseshoes seen by this middle point
        horseshoes = view(effective_horseshoes, :, :, :, mi)

        for ei in 1:nelements               # Iterate over the effective horseshoes seen by this middle point

            # Add surface contribution
            U_vortexring(
                            view(horseshoes, :, :, ei),   # All nodes in this horseshoe
                            1:4,                          # Indices of nodes that make this horseshoe (closed ring)
                            1.0,                          # Unitary strength
                            targets,                      # Midpoint as the target
                            view(Geff, mi:mi, ei:ei);     # Velocity of ei-th horseshoe on the mi-th midpoint
                            dot_with=view(dot_with, :, mi:mi), # Dot the velocity by the orthonormal vector of this midpoint
                            optargs...
                            )
            
            # Add wake contribution
            da1 = Dinfs[1, 1, ei]
            da2 = Dinfs[2, 1, ei]
            da3 = Dinfs[3, 1, ei]
            db1 = Dinfs[1, 2, ei]
            db2 = Dinfs[2, 2, ei]
            db3 = Dinfs[3, 2, ei]

            U_semiinfinite_horseshoe(
                            view(horseshoes, :, :, ei),   # All nodes in this horseshoe
                            TE,                           # Indices of nodes that make the shedding edge
                            da1, da2, da3,                # Semi-infinite direction da
                            db1, db2, db3,                # Semi-infinite direction db
                            1.0,                          # Unitary strength
                            targets,                      # Midpoint as the target
                            view(Geff, mi:mi, ei:ei);     # Velocity of ei-th horseshoe on the mi-th midpoint
                            dot_with=view(dot_with, :, mi:mi), # Dot the velocity by the orthonormal vector of this midpoint
                            optargs...
                            )

        end

    end

end


"""
Self-induced velocity using the effective lifting line geometry for each 
stripwise element
"""
function selfUind!(self::LiftingLine; optargs...)
    selfUind!(self, self.Us; optargs...)
end
function selfUind!(self::LiftingLine, out; optargs...)
    selfUind!(self, self.Gammas, out; optargs...)
end

function selfUind!(self::LiftingLine, Gammas::AbstractVector, 
                    out::AbstractMatrix; precomputed=true, optargs...)
    if precomputed
        _selfUind_precomputed!(self, Gammas, out; optargs...)
    else
        _selfUind_lazy!(self, Gammas, out; optargs...)
    end
end

function _selfUind_lazy!(self::LiftingLine, Gammas::AbstractVector, 
                    out::AbstractMatrix; optargs...)

    for ei in 1:self.nelements                      # Iterate over stripwise elements

        targets = view(self.midpoints, :, ei:ei)
        this_out = view(out, :, ei:ei)
        effective_horseshoes = view(self.effective_horseshoes, :, :, :, ei)

        # Velocity of all the effective horseshoes on this stripwise element
        Uind!(self, Gammas, targets, this_out; 
                                    horseshoes=effective_horseshoes, optargs...)
    end

end

function _selfUind_precomputed!(self::LiftingLine, Gammas::AbstractVector, 
                                out::AbstractMatrix)

    # Convert current velocity in out from global to local coordinates
    for ei in 1:self.nelements

        # Global coordinates
        U1 = out[1, ei]
        U2 = out[2, ei]
        U3 = out[3, ei]

        # Local coordinates
        out[1, ei] = U1*self.swepttangents[1, ei] + U2*self.swepttangents[2, ei] + U3*self.swepttangents[3, ei]
        out[2, ei] = U1*self.lines[1, ei] + U2*self.lines[2, ei] + U3*self.lines[3, ei]
        out[3, ei] = U1*self.sweptnormals[1, ei] + U2*self.sweptnormals[2, ei] + U3*self.sweptnormals[3, ei]

    end

    # Calculate velocity by horseshoes in local coordinates, adding them to out
    # without overwritting
    LA.mul!(view(out, 1, :), view(self.Geff, :, :, 1), Gammas, 1, 1)
    LA.mul!(view(out, 2, :), view(self.Geff, :, :, 2), Gammas, 1, 1)
    LA.mul!(view(out, 3, :), view(self.Geff, :, :, 3), Gammas, 1, 1)

    # Re-project the velocity from local coordinates to global
    for ei in 1:self.nelements

        # Local coordinates
        U1 = out[1, ei]
        U2 = out[2, ei]
        U3 = out[3, ei]

        # Global coordinates
        for i in 1:3
            out[i, ei] = U1*self.swepttangents[i, ei] + U2*self.lines[i, ei] + U3*self.sweptnormals[i, ei]
        end

    end

end

"""
    Uind!(self::LiftingLine, targets, out, args...; optargs...)

Returns the velocity induced by the wing on the targets `targets`, which is a
3xn matrix. It adds the velocity at the i-th target to `out[:, i]`.
"""
function Uind!(self::LiftingLine, targets::AbstractMatrix, args...; optargs...)
    Uind!(self, self.Gammas, targets, args...; optargs...)
end

function Uind!(self::LiftingLine, Gammas::AbstractVector, 
                    targets::AbstractMatrix, out::AbstractMatrix; 
                    horseshoes::AbstractArray=self.horseshoes,
                    add_boundvortices=true,
                    optargs...)

    if add_boundvortices

        # Add bound vortex contributions
        for ei in 1:self.nelements                              # Iterate over horseshoes

            # Velocity of i-th horseshoe on every target
            U_vortexring(
                                view(horseshoes, :, :, ei),        # All nodes
                                1:4,                               # Indices of nodes that make this panel (closed ring)
                                Gammas[ei],                        # Horseshoe circulation
                                targets,                           # Targets
                                out;                               # Outputs
                                offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                                cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                                optargs...
                            )
        end

    end

    # Add wake contributions
    TE = [1, size(horseshoes, 2)]                           # Indices of TE nodes in each horseshoe
    
    for ei in 1:self.nelements                              # Iterate over semi-infinite segments

        da1 = self.Dinfs[1, 1, ei]
        da2 = self.Dinfs[2, 1, ei]
        da3 = self.Dinfs[3, 1, ei]
        db1 = self.Dinfs[1, 2, ei]
        db2 = self.Dinfs[2, 2, ei]
        db3 = self.Dinfs[3, 2, ei]

        U_semiinfinite_horseshoe(
                          view(horseshoes, :, :, ei),        # All nodes
                          TE,                                # Indices of nodes that make the shedding edge
                          da1, da2, da3,                     # Semi-infinite direction da
                          db1, db2, db3,                     # Semi-infinite direction db
                          Gammas[ei],                        # Filament circulation
                          targets,                           # Targets
                          out;                               # Outputs
                          offset=self.kerneloffset,          # Offset of kernel to avoid singularities
                          cutoff=self.kernelcutoff,          # Kernel cutoff to avoid singularities
                          optargs...
                         )

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
                                        extrapolatepolar=true, plot_polars=true, 
                                        optargs...)

    # Create baseline stripwise elements from polars
    airfoils = _read_polars(airfoil_distribution; optargs...)

    # Identify StripwiseElement types
    element_types = unique([typeof(airfoil) for (ypos, airfoil) in airfoils])
    element_types = Union{element_types...}

    # Extrapolate polars to +-180 deg
    if extrapolatepolar
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
        axs[1].plot([airfoil.alpha0], [calc_cl(airfoil, airfoil.alpha0)], "*"; color=clr, fmt...)
        axs[2].plot(airfoil.alpha, airfoil.cd, stl; label=L"$2y/b = $"*"$(ypos)", color=clr, fmt...)
        axs[3].plot(airfoil.alpha, airfoil.cm, stl; label=L"$2y/b = $"*"$(ypos)", color=clr, fmt...)

    end

    for (ypos, airfoil) in airfoils_extrapolated
        
        clr = plt.cm.gnuplot(ypos)
        stl = stl_extrap
        fmt = fmt_extrap

        axs[1].plot(airfoil.alpha, airfoil.cl, stl; color=clr, fmt...)
        axs[1].plot([airfoil.alpha0], [calc_cl(airfoil, airfoil.alpha0)], "*"; color=clr, fmt...)
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
        axs[1].plot([airfoil.alpha0], [calc_cl(airfoil, airfoil.alpha0)], "*"; color=clr, fmt...)
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