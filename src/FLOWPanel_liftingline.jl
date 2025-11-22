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
struct LiftingLine

    # User inputs
    grid::gt.Grid                                # Flat-geometry grid

    # Internal properties
    nelements::Int                               # Number of stripwise elements
    ypositions::Vector{Float64}                  # Non-dimensional y-position of nodes, 2*y/b

    function LiftingLine(args...; optargs...)

        # ------------------ DISCRETIZE WING -----------------------------------
        (; b, ypositions, chords, twists, 
        sweeps, dihedrals, spanaxiss, symmetric) = _discretize_wing_parameterization(args...; optargs...)

        # ------------------ PREALLOCATE GRID MEMORY ---------------------------
        P_min = [0, 0, 0]               # Lower boundary span, chord, dummy
        P_max = [1, 1, 0]               # Upper boundary span, chord, dummy

        NDIVS = [length(ypositions)-1, 1, 0]   # Divisions of span, chord, and dummy collapsed dimension (flat surface)

        # Generate parametric grid
        grid = gt.Grid(P_min, P_max, NDIVS)
        
        # ------------------ MORPH GRID INTO WING GEOMETRY ---------------------
        _morph_grid_wing!(grid, b, ypositions, chords, twists, sweeps, dihedrals, 
                                spanaxiss, symmetric; center=true)

        new(grid, grid._ndivsnodes[1]-1, ypositions)

    end

end

"""
Morph the lifting-line wing geometry into a new geometry
"""
function remorph!(self::LiftingLine, args...; recenter=false, optargs...)

    (; b, ypositions, chords, twists, 
    sweeps, dihedrals, spanaxiss, symmetric) = _discretize_wing_parameterization(args...; optargs...)

    _morph_grid_wing!(self.grid, b, ypositions, chords, twists, sweeps, dihedrals, 
                            spanaxiss, symmetric; center=recenter)
end

function save(self::LiftingLine, args...; format="vtk", optargs...)
    return gt.save(self.grid, args...; format, optargs...)
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
                            spanaxiss, symmetric; center=false)

    @assert grid._ndivsnodes[1] == length(ypositions) ""*
        "Invalid grid! Received grid of $(grid._ndivsnodes[1]) spanwise nodes,"*
        " and $(ypositions) y-positions"

    # ------------------ MORPH GRID INTO WING GEOMETRY -------------------------

    linearindices = LinearIndices(grid._ndivsnodes)

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
##### END OF LIFTING LINE ######################################################