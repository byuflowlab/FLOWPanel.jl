In this example we solve the aerodynamics of Weber's $45^\circ$ swept-back 
wing at an angle of attack of $4.2^\circ$ using a the lifting line solver
with a rigid wake model.

# $4.2^\circ$ Angle of Attack
```julia
#=##############################################################################
# DESCRIPTION
    45deg swept-back wing at an angle of attack of 4.2deg. This wing has an
    aspect ratio of 5.0, a RAE 101 airfoil section with 12% thickness, and no
    dihedral, twist, nor taper. This test case matches the experimental setup
    of Weber, J., and Brebner, G., “Low-Speed Tests on 45-deg Swept-Back Wings,
    Part I,” Tech. rep., 1951.

    This example uses the nonlinear lifting line method.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jan 2026
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import FLOWPanel: mean, norm, dot, cross

import PyPlot as plt



run_name        = "ll-weber"                    # Name of this run

save_path       = run_name                      # Where to save outputs
airfoil_path    = joinpath(pnl.examples_path, "data") # Where to find 2D polars

paraview        = false                         # Whether to visualize with Paraview


# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 4.2                           # (deg) angle of attack
magUinf         = 49.7                          # (m/s) freestream velocity
Uinf            = magUinf*[cosd(AOA), 0, sind(AOA)] # Freestream

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
    0.0     0.0
    1.0     0.0
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

Dhat            = Uinf/norm(Uinf)               # Drag direction
Shat            = [0, 1, 0]                     # Span direction
Lhat            = cross(Dhat, Shat)             # Lift direction

X0              = [0.0 * chord_distribution[1, 2]*b, 0, 0] # Center about which to calculate moments
lhat            = Dhat                          # Rolling direction
mhat            = Shat                          # Pitching direction
nhat            = Lhat                          # Yawing direction

cref            = chord_distribution[1, 2]*b    # (m) reference chord
nondim          = 0.5*rho*magUinf^2*b*cref      # Normalization factor


# ------------------ GENERATE LIFTING LINE -------------------------------------

ll = pnl.LiftingLine{Float64}(
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
                                    aoas_initial_guess=AOA, 
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

# NOTE: Coefficients must be evaluated using the velocity from 
#       the effective horseshoes as shown below, which is automatically
#       computed by the solver already, so these lines are commented out to
#       avoid redundant computation
# ll.Us .= Uinfs
# pnl.selfUind!(ll)

# Calculate stripwise coefficients
pnl.calcfield_cl(ll)
pnl.calcfield_cd(ll)
pnl.calcfield_cm(ll)

# Convert velocity to effective swept velocity
# NOTE: Forces are most accurate with the velocity from the original horseshoes,
#       as done in the conditional statement here
if use_Uind_for_force
    ll.Us .= Uinfs
    pnl.Uind!(ll, ll.midpoints, ll.Us)
end
pnl.calc_UΛs!(ll, ll.Us)

# Force per stripwise element integrating lift and drag coefficient
pnl.calcfield_F(ll, rho)

# Integrated force
Ftot = pnl.calcfield_Ftot(ll)

# Integrated force decomposed into lift and drag
LDS = pnl.calcfield_LDS(ll, Lhat, Dhat, Shat)

L = LDS[:, 1]
D = LDS[:, 2]

# Loading distribution (force per unit span)
fs = pnl.calcfield_f(ll)

lds = pnl.decompose(fs, Lhat, Dhat)

ypos = (ll.ypositions[2:end] .+ ll.ypositions[1:end-1]) / 2
l = lds[1, :]
d = lds[2, :]

# Integrated moment
Mtot = pnl.calcfield_Mtot(ll, X0, rho)

# Moment decomposed into axes
lmn = pnl.calcfield_lmn(ll, lhat, mhat, nhat)
roll, pitch, yaw = collect(eachcol(lmn))

# Coefficients
CL = sign(dot(L, Lhat)) * norm(L) / nondim
CD = sign(dot(D, Dhat)) * norm(D) / nondim

cl = l / (nondim/b)
cd = d / (nondim/b)

Cl = sign(dot(roll, lhat)) * norm(roll) / (nondim*cref)
Cm = sign(dot(pitch, mhat)) * norm(pitch) / (nondim*cref)
Cn = sign(dot(yaw, nhat)) * norm(yaw) / (nondim*cref)


# ------------------ OUTPUT SOLUTION -------------------------------------------
if !isnothing(save_path)
    
    str = pnl.save(ll, run_name; path=save_path, debug=true) # Use `debug=true` to output the effective horseshoes

    if paraview
        run(`paraview --data=$(joinpath(save_path, str))`)
    end
    
end


```
(see the complete example under
[examples/liftingline_weber.jl](https://github.com/byuflowlab/FLOWPanel.jl/blob/master/examples/liftingline_weber.jl)
to see how to postprocess the spanwise loading that is plotted below)

```@raw html
<center>
    <br><br><b>Spanwise loading distribution</b>
    <img src="../../assets/images/ll-weber-loading.png" alt="Pic here" style="width: 100%;"/>
</center>
```


|           | Experimental  | FLOWPanel Lifting Line    | Error | `VSPAERO 3D Panel` |
| --------: | :-----------: | :-----------------------: | :---- |  :----: |
| $C_L$   | 0.238         | 0.247    | 3.6% | *`0.257`* |
| $C_D$   | 0.005         | 0.0055    | 10.2% | *`0.0033`* |


!!! details "Tip"
    You can also automatically run this example and generate these plots
    with the following command:
    ```julia
    import FLOWPanel as pnl

    include(joinpath(pnl.examples_path, "liftingline_weber.jl"))
    ```

