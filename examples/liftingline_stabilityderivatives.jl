#=##############################################################################
# DESCRIPTION
    
    Example on how to calculate stability derivatives using automatic 
    differentiation. This example uses Weber's wing with 5deg dihedral, and
    the nonlinear lifting line method.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jan 2026
  * License   : MIT License
=###############################################################################

import ForwardDiff: Dual, Partials, value, partials

import FLOWPanel as pnl
import FLOWPanel: mean, norm, dot, cross

import PythonPlot as plt
import PythonPlot: @L_str
include(joinpath(pnl.examples_path, "plotformat.jl"))

import CSV
import DataFrames: DataFrame


run_name        = "ll-stabilityderivatives"     # Name of this run

save_path       = run_name                      # Where to save outputs
airfoil_path    = joinpath(pnl.examples_path, "data") # Where to find 2D polars

paraview        = true                          # Whether to visualize with Paraview


# ----------------- SIMULATION PARAMETERS --------------------------------------

# NOTE: Here we define a sideslip angle in addition to the angle of attack

alpha           = 4.2                           # (deg) angle of attack
beta            = 0.0                           # (deg) sideslip angle

# NOTE: In the dual numbers we will implicitely defined the first partial to be 
#       the derivative w.r.t. angle of attack and the second partial to be
#       the derivative w.r.t. sideslip angle

alpha           = Dual(alpha, Partials((1.0, 0.0))) # Convert angle of attack into dual number for automatic differentiation
beta            = Dual(beta,  Partials((0.0, 1.0))) # Convert sideslip angle into dual number for automatic differentiation
NumType         = typeof(alpha)                 # Number type for LiftingLine


magUinf         = 49.7                          # (m/s) freestream velocity magnitude
Uinf            = magUinf*pnl.direction(; alpha, beta) # Freestream vector

rho             = 1.225                         # (kg/m^3) air density


# ------------------ GEOMETRY PARAMETERS ---------------------------------------

# High-level parameters
b               = 98*0.0254                     # (m) wing span

# Discretization parameters
nelements       = 40                            # Number of stripwise elements per semi-span
discretization  = [                             # Multi-discretization of wing (seg length, ndivs, expansion, central)
                    (1.00,  nelements, 1/10, false),
                  ]
symmetric       = true                          # Whether the wing is symmetric

# Chord distribution (nondim y-position, nondim chord)
chord_distribution = [
#   2*y/b   c/b
    0.0     0.2
    1.0     0.2
]

# Twist distribution (nondim y-position, twist)
twist_distribution = [
#   2*y/b   twist (deg)
    0.0     0.0
    1.0     0.0
]

# Sweep distribution (nondim y-position, sweep)
sweep_distribution = [
#   2*y/b   sweep (deg)
    0.0     45.0
    1.0     45.0
]

# Dihedral distribution (nondim y-position, dihedral)
dihedral_distribution = [
#   2*y/b   dihedral (deg)
    0.0     5.0
    1.0     5.0
]

# Span-axis distribution: chordwise point about which the wing is twisted, 
# swept, and dihedralized (nondim y-position, nondim chord-position)
spanaxis_distribution = [
#   2*y/b   x/c
    0.0     0.25
    1.0     0.25
]

# Airfoil contour distribution (nondim y-position, polar, airfoil type)
airfoil_distribution = [
#    2*y/b  polar file                        airfoil type
    (0.00, "rae101-Re1p7e6-smooth180-2.csv",  pnl.SimpleAirfoil),
    (1.00, "rae101-Re1p7e6-smooth180-2.csv",  pnl.SimpleAirfoil)
]

element_optargs = (;    path = airfoil_path,
                        plot_polars = true,
                        extrapolatepolar = false,   # Whether to extrapolate the 2D polars to ±180 deg
                    )


# ------------------ SOLVER PARAMETERS -----------------------------------------
deltasb         = 1.0                           # Blending distance, deltasb = 2*dy/b
deltajoint      = 1.0                           # Joint distance, deltajoint = dx/c

sigmafactor     = 0.0                           # Dragging line amplification factor (set to -1.0 for rebust post-stall method)
sigmaexponent   = 4.0                           # Dragging line amplification exponent (no effects if `sigmafactor==0.0`)

                                                # Nonlinear solver
solver          = pnl.SimpleNonlinearSolve.SimpleDFSane()              # Indifferent to initial guess, but somewhat not robust
# solver        = pnl.SimpleNonlinearSolve.SimpleTrustRegion()         # Trust region needs a good initial guess, but it converges very reliably

solver_optargs  = (; 
                    abstol = 1e-9,  
                    maxiters = 800,
                    )

align_joints_with_Uinfs = false                 # Whether to align joint bound vortices with the freestream

use_Uind_for_force = true                       # Whether to use Uind as opposed to selfUind for force postprocessing
                                                # (`true` for more accurate spanwise cd distribution, but worse integrated CD)

# NOTE: We need to manually define the direction of lift, drag, and span, as 
#       well as roll, pitching, and yawing moment since we now have a sideslip
#       angle

Dhat            = pnl.direction(; alpha)        # Drag direction
Shat            = [0, 1, 0]                     # Span direction
Lhat            = cross(Dhat, Shat)             # Lift direction

X0              = [-0.25 * chord_distribution[1, 2]*b, 0, 0] # Center about which to calculate moments
lhat            = Dhat                          # Rolling direction
mhat            = Shat                          # Pitching direction
nhat            = Lhat                          # Yawing direction

cref            = chord_distribution[1, 2]*b    # (m) reference chord


# ------------------ GENERATE LIFTING LINE -------------------------------------

ll = pnl.LiftingLine{NumType}(
                                airfoil_distribution; 
                                b, chord_distribution, twist_distribution,
                                sweep_distribution, dihedral_distribution,
                                spanaxis_distribution,
                                discretization,
                                symmetric,
                                deltasb, deltajoint, sigmafactor, sigmaexponent,
                                element_optargs,
                                plot_discretization = true,
                                )

display(ll)


# ------------------ CALL NONLINEAR SOLVER -------------------------------------

# Freestream velocity at each stripwise element
Uinfs = repeat(Uinf, 1, ll.nelements)

# Run solver
result, solver_cache = pnl.solve(ll, Uinfs; 
                                    debug=true,             # `true` returns the residual rms
                                    aoas_initial_guess=alpha, 
                                    align_joints_with_Uinfs, 
                                    solver, solver_optargs)

# Check solver success
success = pnl.SimpleNonlinearSolve.SciMLBase.successful_retcode(result)
@show success

@show solver_cache[:fcalls]
# display(result)

# (Optional) plot residual convergence
fig = plt.figure(figsize=[7, 0.5*5]*7/9)
ax = fig.gca()

ax.plot(1:solver_cache[:fcalls], solver_cache[:residual_rms], "-k", linewidth=1)

ax.set_xlabel("Iteration")
ax.set_xlim([0, solver_cache[:fcalls]])

ax.set_ylabel("Residual RMS")
ax.set_yscale("log")

ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)

fig.tight_layout()


# ------------------ POSTPROCESSING --------------------------------------------

# Calculate force and moment coefficients
calcs = pnl.calc_forcemoment_coefficients(ll, Uinfs, Uinf, 
                                            rho, cref, b;
                                            X0, 
                                            Dhat, Shat, Lhat,
                                            lhat, mhat, nhat,
                                            use_Uind_for_force)

# Unpack calculations
(; CD, CY, CL) = calcs                      # Drag, side, and lift forces
(; Cl, Cm, Cn) = calcs                      # Roll, pitch, and yawing moment
(; q, Aref, bref, cref) = calcs             # Reference dynamic pressure, area, span, and chord

@show dCLdα = partials(CL)[1]
@show dCDdα = partials(CD)[1]
@show dCmdα = partials(Cm)[1]

@show dCLdβ = partials(CL)[2]
@show dCDdβ = partials(CD)[2]
@show dCmdβ = partials(Cm)[2]

@show CL = value(CL)
@show CD = value(CD)
@show Cm = value(Cm)



# ------------------ OUTPUT SOLUTION -------------------------------------------
if !isnothing(save_path)
    
    str = pnl.save(ll, run_name; path=save_path, debug=true) # Use `debug=true` to output the effective horseshoes

    if paraview
        run(`paraview --data=$(joinpath(save_path, str))`)
    end
    
end