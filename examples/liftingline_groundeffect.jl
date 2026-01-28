#=##############################################################################
# DESCRIPTION
    Lifting line in ground effect compared against experimental results 
    described in Kalman et al. (1970), "APPLICATION OF THE DOUBLET-LATTICE METHOD
    TO NONPLANAR CONFIGURATIONS IN SUBSONIC FLOW".


# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jan 2026
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import FLOWPanel: mean, norm, dot, cross

import PythonPlot as plt
import PythonPlot: @L_str, pyconvert
include(joinpath(pnl.examples_path, "plotformat.jl"))


run_name        = "ll-ground"                   # Name of this run
save_path       = run_name                      # Where to save outputs
save_outputs    = false                         # Whether to save outputs for docs or not

cache           = Dict()                        # Initialize model cache
calcs           = Dict()                        # Results will get stored here

# aspect_ratios   = [1, 2, 4]                     # Aspect ratios to sweep, b/c
# h_over_cs       = 0.18:0.02:1.2               # Ground distance ratios to sweep, h/c
# h_over_cs       = 0.3:0.02:1.2
# h_over_cs       = [0.3, 0.5, 1.0, 1.2]
# h_over_cs       = [Inf, 12, 8, 4, 2, 1]
aspect_ratios   = [8]
h_over_cs       = [Inf, 1.0, 0.75, 0.5, 0.25]

for aspect_ratio in aspect_ratios

    calcs[aspect_ratio] = Dict()

    for h_over_c in h_over_cs

        b = 3.0                                 # (m) wing span
        chord = b/aspect_ratio                  # (m) wing chord
        h = h_over_c * chord                    # (m) distance to ground plane

        results = pnl.run_liftingline(;

                stability_derivatives = !true,                   # Whether to calculate stability derivatives
                verbose         = true,
                # save_path       = "temps/",
                # paraview        = true,
                
                # ----------------- SIMULATION PARAMETERS --------------------------------------
                alpha           = 1.0,                          # (deg) angle of attack
                beta            = 0.0,                          # (deg) sideslip angle
                
                magUinf         = 50.0,                         # (m/s) freestream velocity magnitude
                ground_distance = h,                            # (m) distance to ground
                
                # ------------------ GEOMETRY PARAMETERS ---------------------------------------
                
                # High-level parameters
                b,                                              # (m) wing span
                
                # Discretization parameters
                nelements       = 40,                           # Number of stripwise elements per semi-span
                
                # Chord distribution (nondim y-position, nondim chord)
                chord_distribution = [
                #   2*y/b   c/b
                    0.0     chord/b
                    1.0     chord/b
                ],
                
                # Twist distribution (nondim y-position, twist)
                twist_distribution = [
                #   2*y/b   twist (deg)
                    0.0     0.0
                    1.0     0.0
                ],
                
                # Sweep distribution (nondim y-position, sweep)
                sweep_distribution = [
                #   2*y/b   sweep (deg)
                    0.0     0.0
                    1.0     0.0
                ],
                
                # Dihedral distribution (nondim y-position, dihedral)
                dihedral_distribution = [
                #   2*y/b   dihedral (deg)
                    0.0     0.0
                    1.0     0.0
                ],
                
                # Airfoil contour distribution (nondim y-position, polar, airfoil type)
                airfoil_distribution = [
                #    2*y/b  polar file                                airfoil type
                    (0.00, "n64_1_A612-Re0p5e6-neuralfoil180-3.csv",  pnl.SimpleAirfoil),
                    (1.00, "n64_1_A612-Re0p5e6-neuralfoil180-3.csv",  pnl.SimpleAirfoil)
                ],
                
                airfoil_path    = joinpath(pnl.examples_path, "data"), # Where to find 2D polars
                
                
                # ------------------ SOLVER PARAMETERS -----------------------------------------
                deltasb         = 1.0,                          # Blending distance, deltasb = 2*dy/b
                deltajoint      = 1.0,                          # Joint distance, deltajoint = dx/c
                
                sigmafactor     = 0.0,                          # Dragging line amplification factor (set to -1.0 for rebust post-stall method)
                sigmaexponent   = 4.0,                          # Dragging line amplification exponent (no effects if `sigmafactor==0.0`)
                
                # solver          = pnl.optimization_solver,      # Nonlinear solver
                solver = pnl.NonlinearSolve.SimpleBroyden(),
                
                align_joints_with_Uinfs = false,                # Whether to align joint bound vortices with the freestream
                
                use_Uind_for_force = true,                      # Whether to use Uind as opposed to selfUind for force postprocessing
                                                                # (`true` for more accurate spanwise cd distribution, but worse integrated CD)

                distributions = true,                          # Whether to output spanwise distributions

                cache                                           # Model cache
                
        )

        # Fetch results
        (;      CD, CY, CL, Cl, Cm, Cn,                       # Force and moment coefficients
        
                dCDdα, dCYdα, dCLdα, dCldα, dCmdα, dCndα,     # Stability derivatives
                dCDdβ, dCYdβ, dCLdβ, dCldβ, dCmdβ, dCndβ,
        
                q, Aref, bref, cref,                          # Non-dimensionalization parameters
        
                D, Y, L, l, m, n,                             # Dimensional forces and moments (scalars)
        
                Dhat, Shat, Lhat,                             # Direction of forces and moments (vectors)
                lhat, mhat, nhat,
        
                lift, drag, side, roll, pitch, yaw,           # Dimensional forces and moments (vectors)
                Ftot, Mtot,                                   # Total forces and moment vectors

                X0,                                           # Point about which the moments where calculated

                distributions,                                # Spanwise distributions
        
                ll,                                           # Lifting line object that was used to obtain these results
                ) = results


        polar = pnl.run_polarsweep(results.ll, 50.0, 1.225, results.X0, results.cref, results.bref; Aref=results.Aref, 
                                aoa_sweeps = [range(0, -20, step=-0.5), range(0.5, 20, step=0.5)],
                                solver=pnl.optimization_solver,
                                # solver=pnl.NonlinearSolve.SimpleBroyden()
                                )

        calcs[aspect_ratio][h_over_c] = merge(results, (; polar))

    end
end


# ----------------- COMPARISON TO EXPERIMENT -----------------------------------

stl_exp = "-o"
fmt_exp = (; label="Experimental", color="k")

stl_ll = ".-"
fmt_ll = (; label="FLOWPanel LL", color="steelblue", markersize=4, alpha=1.0)

fig = plt.figure(figsize=[7*1, 5*4*0.75]*2/3)
axs = fig.subplots(4, 1)
axs = pyconvert(Array, axs)

# axs = [axs[j, i] for i in 1:size(axs, 1), j in 1:size(axs, 2)]

for aspect_ratio in aspect_ratios

    CLs_ll = [calcs[aspect_ratio][h_over_c].success ? calcs[aspect_ratio][h_over_c].CL : NaN for h_over_c in h_over_cs]
    dCLdαs_ll = [calcs[aspect_ratio][h_over_c].success ? calcs[aspect_ratio][h_over_c].dCLdα*180/pi : NaN for h_over_c in h_over_cs]
    dCDdαs_ll = [calcs[aspect_ratio][h_over_c].success ? calcs[aspect_ratio][h_over_c].dCDdα*180/pi : NaN for h_over_c in h_over_cs]
    dCmdαs_ll = [calcs[aspect_ratio][h_over_c].success ? calcs[aspect_ratio][h_over_c].dCmdα*180/pi : NaN for h_over_c in h_over_cs]

    # dCLdα vs h/c
    ax = axs[1]
    # ax.plot(h_over_cs_exp, dCLdαs_exp, stl_exp; fmt_exp...)
    ax.plot(h_over_cs, dCLdαs_ll, stl_ll; fmt_ll...)
    ax.set_ylabel(L"$C_{L\alpha}$")
    ax.set_ylim([1.0, 8.0])

    # dCDdα vs h/c
    ax = axs[2]
    # ax.plot(h_over_cs_exp, dCDdαs_exp, stl_exp; fmt_exp...)
    ax.plot(h_over_cs, dCDdαs_ll, stl_ll; fmt_ll...)
    ax.set_ylabel(L"$C_{D\alpha}$")

    # dCmdα vs h/c
    ax = axs[3]
    # ax.plot(h_over_cs_exp, dCmdαs_exp, stl_exp; fmt_exp...)
    ax.plot(h_over_cs, dCmdαs_ll, stl_ll; fmt_ll...)
    ax.set_ylabel(L"$C_{m\alpha}$")
    ax.set_ylim([-0.4, 0.2])

    # CL vs h/c
    ax = axs[4]
    ax.plot(h_over_cs, CLs_ll, stl_ll; fmt_ll...)
    ax.set_ylabel(L"$C_L$")
end

for (axi, ax) in enumerate(axs)

    ax.set_xlim([0, 1.4])
    ax.set_xticks(0:0.2:1.4)
    ax.set_xlabel(L"Ground distance $h/c$")

    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)

    ax.legend(loc="best", frameon=false, fontsize=8, reverse=true)
end

fig.tight_layout()

# if save_outputs
#     fig.savefig(joinpath(fig_path, "$(run_name)-sweep-$(suffix).png"),
#                                                 dpi=300, transparent=true)
# end