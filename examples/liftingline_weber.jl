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
import PyPlot: @L_str
include(joinpath(pnl.examples_path, "plotformat.jl"))

import CSV
import DataFrames: DataFrame


run_name        = "ll-weber"                    # Name of this run

save_path       = run_name                      # Where to save outputs
airfoil_path    = joinpath(pnl.examples_path, "data") # Where to find 2D polars

paraview        = false                         # Whether to visualize with Paraview
save_outputs    = false                         # Whether to save outputs for docs or not


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


# ----------------- COMPARISON TO EXPERIMENTAL DATA ----------------------------

# Where to save figures (default to re-generating the figures that are used
# in the docs)
fig_path = joinpath(pnl.examples_path, "..", "docs", "resources", "images")
outdata_path = joinpath(pnl.examples_path, "..", "docs", "resources", "data")

# Weber's experimental loading distribution from Table 3
alphas_web = [2.1, 4.2, 6.3, 8.4, 10.5]
ypos_web = y2b_web = [0.0, 0.041, 0.082, 0.163, 0.245, 0.367, 0.510, 0.653, 0.898, 0.949]
cls_web = [
            0.118 0.121 0.126 0.129 0.129 0.129 0.131 0.125 NaN 0.087;      # AOA = 2.1deg
            0.235 0.241 0.248 0.253 0.251 0.251 0.251 0.246 0.192 0.171;    # AOA = 4.2deg
            0.351 0.358 0.367 0.374 0.375 0.373 0.377 0.365 NaN 0.256;      # AOA = 6.3deg
            0.466 0.476 0.483 0.494 0.494 0.493 0.493 0.48 NaN 0.34;        # AOA = 8.4deg
            0.577 0.589 0.597 0.607 0.611 0.605 0.599 0.587 0.415 0.401     # AOA = 10.5deg
         ]
cds_web = [
            0.044 0.014 0.007 0.002 0.0 -0.001 0.0 -0.001 -0.009 -0.01;     # AOA = 0deg
            0.047 0.016 0.01 0.004 0.002 0.001 0.002 0.001 NaN -0.009;
            0.059 0.025 0.016 0.009 0.007 0.006 0.006 0.004 -0.002 -0.007;
            0.078 0.039 0.028 0.017 0.015 0.012 0.012 0.009 NaN -0.005;
            0.104 0.058 0.045 0.029 0.026 0.023 0.022 0.016 NaN -0.001;
            0.138 0.084 0.065 0.044 0.041 0.037 0.035 0.026 0.009 0.004
         ]

# Add the zero-AOA drag that they say they substracted in the experiment
# cds_web = mapslices(x-> x .+ cds_web[1, :], cds_web[2:end, :]; dims=2)

# Integrated coefficients from Table 4B
CLs_web = [0.121, 0.238, 0.350, 0.456, 0.559]
CDs_web = [nothing, 0.005, 0.012, 0.022, 0.035]

# File with VSP solution
vsp_file = joinpath(pnl.examples_path, "data", "weber_vspaero.csv")
# color_vsp = "tab:orange"
color_vsp = "goldenrod"
alpha_vsp = 0.75

lbl_web = "Experimental"
lbl_vsp = "VSPAERO Panel"

lbl_ll = "FLOWPanel LL"
color_ll = "steelblue"

# --------- Integrated forces: lift and induced drag
CLexp = CLs_web[2]
CDexp = CDs_web[2]

@show CL
@show CD

@show CLexp
@show CDexp

println("CL error:\t$(round(abs(CL-CLexp)/CLexp*100, digits=2))%")
println("CD error:\t$(round(abs(CD-CDexp)/CDexp*100, digits=2))%")

if save_outputs
    str = """
    |           | Experimental  | FLOWPanel Lifting Line    | Error | `VSPAERO 3D Panel` |
    | --------: | :-----------: | :-----------------------: | :---- |  :----: |
    | \$C_L\$   | 0.238         | $(round(CL, digits=3))    | $(round(abs(CL-CLexp)/CLexp*100, digits=1))% | *`0.257`* |
    | \$C_D\$   | 0.005         | $(round(CD, digits=4))    | $(round(abs(CD-CDexp)/CDexp*100, digits=1))% | *`0.0033`* |
    """

    open(joinpath(outdata_path, run_name*"-CLCD.md"), "w") do f
        println(f, str)
    end
end


# --------- Spanwise loading

cls_vsp = [3.884799999999999920e-01,5.206999999999999823e-02,1.146999999999999964e-01,1.115700000000000025e-01,1.437400000000000067e-01,1.428300000000000125e-01,1.532600000000000073e-01,1.908199999999999896e-01,1.643699999999999883e-01,2.065800000000000136e-01,1.799199999999999966e-01,1.906599999999999961e-01,1.981499999999999928e-01,2.024299999999999988e-01,2.312200000000000089e-01,2.334199999999999886e-01,2.392199999999999882e-01,2.444699999999999929e-01,2.497900000000000120e-01,2.543799999999999950e-01,2.581600000000000006e-01,2.615899999999999892e-01,2.644699999999999829e-01,2.663400000000000212e-01,2.687999999999999834e-01,2.700699999999999767e-01,2.715299999999999936e-01,2.722600000000000020e-01,2.727100000000000080e-01,2.731500000000000039e-01,2.733099999999999974e-01,2.731600000000000139e-01,2.727100000000000080e-01,2.717899999999999761e-01,2.708800000000000097e-01,2.694500000000000228e-01,2.685199999999999809e-01,2.666000000000000036e-01,2.648699999999999943e-01,2.627200000000000091e-01,2.609699999999999798e-01,2.584899999999999975e-01,2.561399999999999788e-01,2.536499999999999866e-01,2.513900000000000023e-01,2.480900000000000050e-01,2.457900000000000085e-01,2.442300000000000026e-01,2.435299999999999965e-01,2.435499999999999887e-01,2.441399999999999959e-01,2.456399999999999972e-01,2.479200000000000015e-01,2.513199999999999878e-01,2.539399999999999991e-01,2.564500000000000113e-01,2.589400000000000035e-01,2.612999999999999767e-01,2.630100000000000215e-01,2.652800000000000158e-01,2.671700000000000186e-01,2.690299999999999914e-01,2.699500000000000233e-01,2.711799999999999766e-01,2.722200000000000175e-01,2.733700000000000019e-01,2.741700000000000248e-01,2.745299999999999963e-01,2.745699999999999807e-01,2.743800000000000128e-01,2.740699999999999803e-01,2.735699999999999799e-01,2.722300000000000275e-01,2.706500000000000017e-01,2.685600000000000209e-01,2.669900000000000051e-01,2.642800000000000149e-01,2.606299999999999728e-01,2.571800000000000197e-01,2.533299999999999996e-01,2.485300000000000009e-01,2.437300000000000022e-01,2.377299999999999969e-01,2.355700000000000016e-01,2.069000000000000006e-01,2.023099999999999898e-01,1.945399999999999907e-01,1.833100000000000007e-01,2.095499999999999863e-01,1.673099999999999865e-01,1.935300000000000076e-01,1.562300000000000078e-01,1.461800000000000044e-01,1.464999999999999913e-01,1.139799999999999980e-01,1.164099999999999996e-01,5.419000000000000206e-02,3.897700000000000053e-01]
cds_vsp = [1.577999999999999889e-02,-2.197999999999999954e-02,-7.799999999999999859e-04,-8.139999999999999666e-03,-4.579999999999999870e-03,-2.251999999999999835e-02,-7.500000000000000156e-04,-4.530000000000000172e-03,-3.010000000000000089e-03,-1.217999999999999985e-02,-8.840000000000000635e-03,-9.050000000000000752e-03,-1.330000000000000019e-03,-1.916999999999999954e-02,-6.020000000000000177e-03,-1.685000000000000039e-02,-7.100000000000000408e-03,-1.976999999999999938e-02,-7.749999999999999944e-03,-1.623000000000000137e-02,-7.539999999999999827e-03,-1.620999999999999872e-02,-8.099999999999999561e-03,-1.097000000000000058e-02,-9.820000000000000603e-03,-1.315000000000000023e-02,-7.380000000000000275e-03,-9.579999999999999974e-03,-5.130000000000000011e-03,-1.076000000000000047e-02,-4.259999999999999898e-03,-7.089999999999999948e-03,-2.899999999999999800e-03,-7.119999999999999593e-03,-8.000000000000000383e-04,-2.809999999999999998e-03,-1.600000000000000131e-04,-1.920000000000000049e-03,3.870000000000000176e-03,-5.599999999999999509e-04,1.699999999999999905e-03,-1.510000000000000057e-03,4.819999999999999632e-03,3.000000000000000076e-05,7.510000000000000182e-03,1.270000000000000078e-03,9.379999999999999449e-03,5.309999999999999616e-03,2.518000000000000099e-02,2.522999999999999896e-02,5.239999999999999866e-03,9.329999999999999752e-03,1.299999999999999940e-03,7.459999999999999618e-03,-8.000000000000000654e-05,4.740000000000000289e-03,-1.649999999999999991e-03,1.559999999999999972e-03,-7.399999999999999894e-04,3.680000000000000111e-03,-2.010000000000000068e-03,1.299999999999999886e-04,-2.599999999999999881e-03,-7.100000000000000191e-04,-7.139999999999999646e-03,-2.949999999999999931e-03,-7.040000000000000251e-03,-4.170000000000000095e-03,-1.059000000000000045e-02,-5.089999999999999906e-03,-9.679999999999999369e-03,-7.450000000000000025e-03,-1.316999999999999942e-02,-9.900000000000000813e-03,-1.111999999999999968e-02,-8.259999999999999981e-03,-1.632000000000000117e-02,-7.599999999999999985e-03,-1.634999999999999995e-02,-7.889999999999999444e-03,-1.965000000000000080e-02,-6.969999999999999633e-03,-1.692000000000000101e-02,-6.190000000000000190e-03,-1.922000000000000097e-02,-9.700000000000000505e-04,-8.649999999999999703e-03,-8.750000000000000833e-03,-1.280000000000000061e-02,-3.160000000000000048e-03,-4.959999999999999999e-03,-8.499999999999999526e-04,-2.248999999999999957e-02,-4.680000000000000132e-03,-7.919999999999999957e-03,-5.900000000000000296e-04,-2.031000000000000166e-02,1.408000000000000050e-02]
ypos_vsp = [1.241400000000000059e+00,1.234960000000000058e+00,1.227500000000000036e+00,1.219410000000000105e+00,1.210609999999999964e+00,1.201149999999999940e+00,1.190930000000000044e+00,1.179920000000000080e+00,1.168129999999999891e+00,1.155370000000000008e+00,1.141739999999999977e+00,1.127010000000000067e+00,1.111220000000000097e+00,1.094319999999999959e+00,1.076310000000000100e+00,1.057239999999999958e+00,1.036969999999999947e+00,1.015460000000000029e+00,9.927099999999999813e-01,9.686799999999999855e-01,9.433700000000000419e-01,9.167800000000000393e-01,8.889399999999999524e-01,8.598599999999999577e-01,8.295900000000000496e-01,7.981700000000000461e-01,7.656699999999999617e-01,7.321699999999999875e-01,6.977499999999999813e-01,6.625299999999999523e-01,6.266100000000000003e-01,5.901199999999999779e-01,5.532000000000000250e-01,5.159700000000000397e-01,4.786000000000000254e-01,4.412200000000000011e-01,4.039900000000000158e-01,3.670399999999999774e-01,3.305100000000000260e-01,2.945499999999999785e-01,2.592700000000000005e-01,2.247899999999999898e-01,1.912200000000000011e-01,1.586500000000000132e-01,1.271500000000000130e-01,9.679000000000000103e-02,6.761999999999999955e-02,3.969000000000000306e-02,1.302000000000000032e-02,-1.302000000000000032e-02,-3.969000000000000306e-02,-6.761999999999999955e-02,-9.679000000000000103e-02,-1.271500000000000130e-01,-1.586500000000000132e-01,-1.912200000000000011e-01,-2.247899999999999898e-01,-2.592700000000000005e-01,-2.945499999999999785e-01,-3.305100000000000260e-01,-3.670399999999999774e-01,-4.039900000000000158e-01,-4.412200000000000011e-01,-4.786000000000000254e-01,-5.159700000000000397e-01,-5.532000000000000250e-01,-5.901199999999999779e-01,-6.266100000000000003e-01,-6.625299999999999523e-01,-6.977499999999999813e-01,-7.321699999999999875e-01,-7.656699999999999617e-01,-7.981700000000000461e-01,-8.295900000000000496e-01,-8.598599999999999577e-01,-8.889399999999999524e-01,-9.167800000000000393e-01,-9.433700000000000419e-01,-9.686799999999999855e-01,-9.927099999999999813e-01,-1.015460000000000029e+00,-1.036969999999999947e+00,-1.057239999999999958e+00,-1.076310000000000100e+00,-1.094319999999999959e+00,-1.111220000000000097e+00,-1.127010000000000067e+00,-1.141739999999999977e+00,-1.155370000000000008e+00,-1.168129999999999891e+00,-1.179920000000000080e+00,-1.190930000000000044e+00,-1.201149999999999940e+00,-1.210609999999999964e+00,-1.219410000000000105e+00,-1.227500000000000036e+00,-1.234960000000000058e+00,-1.241400000000000059e+00] / (b/2)

cls = cl
cds = cd

fig = plt.figure(figsize = [7*2, 0.75*5]*2/3 )
axs = fig.subplots(1, 2)

ax = axs[1]

ax.plot(ypos_web, cls_web[2, :], ":ok", linewidth=0.5, label=lbl_web)
ax.plot(-ypos_web, cls_web[2, :], ":ok", linewidth=0.5)

ax.plot(ypos_vsp, cls_vsp, "-", label=lbl_vsp, color=color_vsp, alpha=alpha_vsp)

ax.plot(ypos, cls, "-", label=lbl_ll, color=color_ll, alpha=0.9)

ax.set_ylabel(L"Sectional lift $c_\ell$")
ax.set_ylim([0, 0.4])
ax.set_yticks(0:0.1:0.4)

ax = axs[2]

ax.plot(ypos_web, cds_web[3, :], ":ok", linewidth=0.5, label=lbl_web)
ax.plot(-ypos_web, cds_web[3, :], ":ok", linewidth=0.5)

ax.plot(ypos_vsp, cds_vsp, "-", label=lbl_vsp, color=color_vsp, linewidth=1, alpha=alpha_vsp)

ax.plot(ypos, cds, "-", label=lbl_ll, color=color_ll, alpha=0.9)

ax.set_ylabel(L"Sectional drag $c_d$")
ax.set_ylim([-0.02, 0.06])
ax.set_yticks(-0.02:0.02:0.06)


# ax = axs[3]
# ax.plot(ypos, ll.Gammas, "-", label=lbl_ll, color=color_ll, alpha=0.9)

# ax.set_ylabel(L"Circulation $\Gamma$ (m$^2$/s)")
# ax.set_ylim([0, 6])

# ax = axs[4]
# ax.plot(ypos, ll.aoas, "-", label=lbl_ll, color=color_ll, alpha=0.9)

# ax.set_ylabel(L"AOA ($^\circ$)")
# # ax.set_ylim([0, 6])
# ax.set_ylim([-AOA, 2.0*AOA])


# ax = axs[5]
# ax.plot(ypos, [pnl.calc_cl(element, aoa) for (element, aoa) in zip(ll.elements, ll.aoas)], 
#                                     "-", label=lbl_ll, color=color_ll, alpha=0.9)

# ax.set_ylabel(L"$c_\ell(\alpha)$")
# ax.set_ylim([0, 0.8])


for (axi, ax) in enumerate(axs)
    ax.set_xlabel(L"Span position $2y/b$")
    ax.set_xlim([-1, 1])
    ax.set_xticks(-1:0.5:1)
    
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    
    if axi==1
        ax.legend(loc="best", frameon=false, fontsize=10)
    end
end
    
fig.tight_layout()

if save_outputs
    fig.savefig(joinpath(fig_path, "$(run_name)-loading.png"),
                                                dpi=300, transparent=true)
end













# ----------------- AOA SWEEP --------------------------------------------------

# Sequence of sweeps to run
# NOTE: To help convergence and speed it up the sweep, we recommend starting
#       each sweep at 0deg AOA since the sweep steps over points using the last
#       converged solution as the initial guess
sweep1 = [0, 2.1, 4.2, 6.3, 8.4, 10.5]      # Same AOAs than Weber's experiment
sweep2 = range(0, -50, step=-0.5)           # Sweep from 0 into deep negative stall (-50deg)
sweep3 = range(0, 50, step=0.5)             # Sweep from 0 into deep positive stall (50deg)

distributions = []

@time wingpolar = pnl.run_polarsweep(ll,
                            magUinf, rho, X0, cref, b;
                            aoa_sweeps = (sweep1, sweep2, sweep3),
                            # sweepname = run_name,
                            # plots_path= save_outputs ? fig_path : nothing,
                            # extraplots_path = save_outputs ? fig_path : nothing,
                            output_distributions=distributions,
                            solver,
                            solver_optargs,
                            align_joints_with_Uinfs, 
                            use_Uind_for_force
                        )

# ----------------- COMPARISON TO EXPERIMENTAL SWEEP ---------------------------


data_machupx = Dict("alphas" => [-10.,  -9.,  -8.,  -7.,  -6.,  -5.,  -4.,  -3.,  -2.,  -1.,   0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,  11.],
 "CLs" => [-0.4939648598822694,
  -0.4444346397982473,
  -0.47480103840512605,
  -0.42414095499754806,
  -0.36992927522536434,
  -0.3094171238910714,
  -0.23930388354878754,
  -0.17808475251266476,
  -0.11925004497556596,
  -0.059966215576743254,
  1.366679602399278e-05,
  0.05995851045764972,
  0.11956872611797778,
  0.17599122860401487,
  0.2406293029383441,
  0.3082255551517388,
  0.37113242675690006,
  0.4239420066905544,
  0.47544976737499633,
  0.5227324043055217,
  0.516056449700583,
  0.5322438510926505],
 "CDs" => [0.040237800503509405,
  0.03079483494461626,
  0.043057351009387394,
  0.032320943360616756,
  0.022788709546785125,
  0.014510514319782659,
  0.00750172140813918,
  0.0021288228308744103,
  -0.0017126324862144803,
  -0.004064777146157755,
  -0.0048449943462536334,
  -0.0040583450674485004,
  -0.0017190871016840282,
  0.0020676803524148084,
  0.0073560514490787965,
  0.014431987045697457,
  0.02287357188287931,
  0.0322996947043889,
  0.0431106703848565,
  0.05481600105509516,
  0.04897820125626967,
  0.05028709695541397],
 "CMs" => [0.5986595473657705,
  0.5561746175221842,
  0.5787161415320228,
  0.5180836656498263,
  0.452908222110708,
  0.3790805216893495,
  0.29180254040878706,
  0.21651899289322413,
  0.14501656696731627,
  0.07298229825529091,
  -1.6735171985013218e-05,
  -0.07296517080052517,
  -0.14545667524810976,
  -0.2135657723524868,
  -0.293467565310629,
  -0.37741146258040253,
  -0.45421411960807184,
  -0.5178780874455109,
  -0.5796307969056694,
  -0.6357358892238342,
  -0.6338902900354635,
  -0.6445854052130995])

data_aerosandbox = Dict(
 "CL" => [-4.99490119e-01, -5.15191316e-01, -4.99853073e-01, -4.74567375e-01,
        -4.45511267e-01, -4.13997174e-01, -3.80073805e-01, -3.43332909e-01,
        -3.03114587e-01, -2.59215428e-01, -2.12935841e-01, -1.67995225e-01,
        -1.25635784e-01, -8.50368120e-02, -4.30390958e-02,  6.16549224e-19,
         4.30390958e-02,  8.50368120e-02,  1.25635784e-01,  1.67995225e-01,
         2.12935841e-01,  2.59215428e-01,  3.03114587e-01,  3.43332909e-01,
         3.80073805e-01,  4.13997174e-01,  4.45511267e-01,  4.74567375e-01,
         4.99853073e-01,  5.15191316e-01,  4.99490119e-01],
 "CD" => [0.07199318, 0.05604791, 0.03946754, 0.03479451, 0.03052084,
        0.02666502, 0.02315368, 0.01986456, 0.01673803, 0.01375567,
        0.01089447, 0.00823693, 0.0065706 , 0.005744  , 0.00522838,
        0.00503451, 0.00522838, 0.005744  , 0.0065706 , 0.00823693,
        0.01089447, 0.01375567, 0.01673803, 0.01986456, 0.02315368,
        0.02666502, 0.03052084, 0.03479451, 0.03946754, 0.05604791,
        0.07199318],
 "Cm" => [ 7.46878772e-01,  7.85000441e-01,  7.55916379e-01,  7.23611735e-01,
         6.84354474e-01,  6.40751176e-01,  5.92685083e-01,  5.38901924e-01,
         4.77669080e-01,  4.07964830e-01,  3.32027700e-01,  2.59081304e-01,
         1.93198450e-01,  1.31398862e-01,  6.67966657e-02, -2.01321212e-17,
        -6.67966657e-02, -1.31398862e-01, -1.93198450e-01, -2.59081304e-01,
        -3.32027700e-01, -4.07964830e-01, -4.77669080e-01, -5.38901924e-01,
        -5.92685083e-01, -6.40751176e-01, -6.84354474e-01, -7.23611735e-01,
        -7.55916379e-01, -7.85000441e-01, -7.46878772e-01],
 "alpha" => [-15, -14, -13, -12, -11, -10,  -9,  -8,  -7,  -6,  -5,  -4,  -3,
         -2,  -1,   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
         11,  12,  13,  14,  15])


# --------- Load distribution

function plot_distribution(distributions, sweep1; suffix="loading")

    fig = plt.figure(figsize=[7*2, 5*1*0.8]*2/3)
    axs = fig.subplots(1, 2)

    ypos = distributions[1].yposition
    aoas = [aoa for (aoa, cl) in zip(distributions[1].AOA, distributions[1].cl) if aoa in sweep1]
    cls = [cl for (aoa, cl) in zip(distributions[1].AOA, distributions[1].cl) if aoa in sweep1]
    cds = [cd for (aoa, cd) in zip(distributions[1].AOA, distributions[1].cd) if aoa in sweep1]

    for (axi, (ax, vals_exp)) in enumerate(zip(axs, [cls_web, cds_web[2:end, :]]))

        first = true

        for (AOA, cl, cd) in zip(aoas, cls, cds)

            rowi = findfirst(a -> a==AOA, alphas_web)

            if rowi != nothing && AOA in (axi==1 ? [2.1, 4.2, 6.3, 8.4] : [4.2, 6.3, 8.4])

                # Filter out NaNs
                ys = vals_exp[rowi, :]
                xs = [val for (vali, val) in enumerate(y2b_web) if !isnan(ys[vali])]
                ys = [val for (vali, val) in enumerate(ys) if !isnan(ys[vali])]

                # Plot experimental
                for f in [-1, 1]
                    ax.plot(f*xs, ys, "o--k",
                                label=("Experimental"^(f==1))^first,
                                linewidth=0.5, markersize=5, alpha=1.0)
                end

                # Plot FLOWPanel
                ax.plot(ypos, axi==1 ? cl : cd, "-", label=lbl_ll^first,
                                color=color_ll, markersize=8, linewidth=1)

                first = false
            end

        end

        xlims = [0, 1]
        ax.set_xlim(xlims)
        ax.set_xticks(xlims[1]:0.2:xlims[end])
        ax.set_xlabel(L"Span position $2y/b$")

        if axi==1
            ylims = [0, 0.6]
            ax.set_ylim(ylims)
            ax.set_yticks(ylims[1]:0.2:ylims[end])
            ax.set_ylabel(L"Sectional lift $c_\ell$")

            ax.legend(loc="best", frameon=false, fontsize=6)

            ax.annotate(L"\alpha=8.4^\circ", [0.38, 0.53], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)
            ax.annotate(L"\alpha=6.3^\circ", [0.38, 0.402], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)
            ax.annotate(L"\alpha=4.2^\circ", [0.38, 0.275], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)
            ax.annotate(L"\alpha=2.1^\circ", [0.38, 0.145], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)
        else
            ylims = [-0.04, 0.12]
            ax.set_ylim(ylims)
            ax.set_yticks(ylims[1]:0.04:ylims[end])
            ax.set_ylabel(L"Sectional drag $c_d$")


            ax.annotate(L"\alpha=8.4^\circ", [0.25, 0.030], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6, rotation=-10)
            ax.annotate(L"\alpha=4.2^\circ", [0.25, -0.005], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6, rotation=-5)

            ax.annotate(L"\alpha=6.3^\circ", [0.5, 0.035], xycoords="data", fontsize=ANNOT_SIZE, color="black", alpha=0.6)

            ax.annotate("", [0.4, 0.0145], xycoords="data",
                        xytext=[0.5, 0.035], textcoords="data",
                        arrowprops=Dict(:facecolor=>"black", :linewidth=>0, :alpha=>0.4,
                                        :shrink=>0, :width=>1.0, :headwidth=>5.0, :headlength=>7))
        end

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)
    end

    fig.tight_layout()

    if save_outputs
        fig.savefig(joinpath(fig_path, "$(run_name)-sweep-$(suffix).png"),
                                                    dpi=300, transparent=true)
    end
end

plot_distribution(distributions, sweep1; suffix="loading")



# --------- Integrated forces and moment: lift, drag, and pitching moment

stl_web = "-o"
fmt_web = (; label=lbl_web, color="k")

stl_asb = "-^"
fmt_asb = (; label="AeroSandbox LL", color="green", markersize=3, alpha=0.2)

stl_mux = "-s"
fmt_mux = (; label="MachUpX LL", color="orchid", markersize=3, alpha=0.3)

stl_vsp = "-o"
fmt_vsp = (; label=lbl_vsp, color="goldenrod", markersize=4, alpha=0.6)

stl_ll = ".-"
fmt_ll = (; label=lbl_ll, color=color_ll, markersize=4, alpha=1.0)

function plot_polars(wingpolar; suffix="CLCDCm", 
                        aoalims=[-17, 17], aoaticks=-15:5:15,
                        CLlims=[-0.75, 1.0], CLticks=-0.5:0.5:CLlims[end],
                        CDlims=[0, 0.06], CDticks=CDlims[1]:0.02:CDlims[end],
                        Cmlims=[-1.5, 1.5], Cmticks=Cmlims[1]:0.5:Cmlims[end]
                        )
    AOAs = wingpolar.AOA
    CLs = wingpolar.CL
    CDs = wingpolar.CD
    Cms = wingpolar.Cm

    # VSPAERO CL and CD
    data_vsp = CSV.read(vsp_file, DataFrame; skipto=397, limit=419-397+1)
    alphas_vsp = [val for val in data_vsp[1, 2:end]]
    CDi_vsp = [val for val in data_vsp[3, 2:end]]
    CDtot_vsp = [val for val in data_vsp[6, 2:end]]
    CL_vsp = [val for val in data_vsp[11, 2:end]]
    CMy_vsp = [val for val in data_vsp[16, 2:end]]

    fig = plt.figure(figsize=[7*2, 5*2*0.75]*2/3)
    axs = fig.subplots(2, 2)

    axs = [axs[j, i] for i in 1:size(axs, 1), j in 1:size(axs, 2)]

    # CL vs AOA
    ax = axs[1]
    ax.plot(alphas_web, CLs_web, stl_web; fmt_web...)
    ax.plot(alphas_vsp, CL_vsp, stl_vsp; fmt_vsp...)
    ax.plot(data_aerosandbox["alpha"], data_aerosandbox["CL"], stl_asb; fmt_asb...)
    ax.plot(data_machupx["alphas"], data_machupx["CLs"], stl_mux; fmt_mux...)
    ax.plot(AOAs, CLs, stl_ll; fmt_ll...)

    if !isnothing(CLlims)
        ax.set_ylim(CLlims)
    end
    if !isnothing(CLticks)
        ax.set_yticks(CLticks)
    end
    ax.set_ylabel(L"Lift coeff. $C_L$")

    # CD vs AOA
    ax = axs[2]
    ax.plot(alphas_web, CDs_web, stl_web; fmt_web...)
    ax.plot(alphas_vsp, CDtot_vsp, stl_vsp; fmt_vsp...)
    ax.plot(data_aerosandbox["alpha"], data_aerosandbox["CD"], stl_asb; fmt_asb...)
    # ax.plot(data_machupx["alphas"], data_machupx["CDs"], stl_mux; fmt_mux...)
    ax.plot(AOAs, CDs, stl_ll; fmt_ll...)

    if !isnothing(CDlims)
        ax.set_ylim(CDlims)
    end
    if !isnothing(CDticks)
        ax.set_yticks(CDticks)
    end
    ax.set_ylabel(L"Drag coeff. $C_D$")

    # CL vs CD
    ax = axs[3]
    ax.plot(CDs_web, CLs_web, stl_web; fmt_web...)
    ax.plot(CDtot_vsp, CL_vsp, stl_vsp; fmt_vsp...)
    ax.plot(data_aerosandbox["CD"], data_aerosandbox["CL"], stl_asb; fmt_asb...)
    # ax.plot(data_machupx["CDs"], data_machupx["CLs"], stl_mux; fmt_mux...)
    ax.plot(CDs, CLs, stl_ll; fmt_ll...)

    if !isnothing(CLlims)
        ax.set_ylim(CLlims)
    end
    if !isnothing(CLticks)
        ax.set_yticks(CLticks)
    end
    ax.set_ylabel(L"Lift coeff. $C_L$")
    if !isnothing(CDlims)
        ax.set_xlim(CDlims)
    end
    if !isnothing(CDticks)
        ax.set_xticks(CDticks)
    end
    ax.set_xlabel(L"Drag coeff. $C_D$")

    # Cm vs AOA
    ax = axs[4]
    ax.plot(alphas_vsp, CMy_vsp, stl_vsp; fmt_vsp...)
    ax.plot(data_aerosandbox["alpha"], data_aerosandbox["Cm"], stl_asb; fmt_asb...)
    ax.plot(data_machupx["alphas"], data_machupx["CMs"], stl_mux; fmt_mux...)
    ax.plot(AOAs, Cms, stl_ll; fmt_ll...)

    ax.set_ylim(Cmlims)
    ax.set_yticks(Cmticks)
    ax.set_ylabel(L"Pitching moment $C_m$")

    for (axi, ax) in enumerate(axs)

        if axi != 3
            ax.set_xlim(aoalims)
            ax.set_xticks(aoaticks)
            ax.set_xticklabels(["$val"*L"^\circ" for val in aoaticks])
            ax.set_xlabel(L"Angle of attack $\alpha$")
        end

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        ax.legend(loc="best", frameon=false, fontsize=8, reverse=true)
    end

    fig.tight_layout()

    if save_outputs
        fig.savefig(joinpath(fig_path, "$(run_name)-sweep-$(suffix).png"),
                                                    dpi=300, transparent=true)
    end

end

plot_polars(wingpolar; suffix="CLCDCm")

plot_polars(wingpolar; suffix="CLCDCm-zoomout", 
                        aoalims=[-40, 40], aoaticks=-40:20:40,
                        CDlims=[0, 0.4], CDticks=0:0.1:0.4,
                        )



# Proceed to do the following sweeps only if we are regenerating the docs
if save_outputs

    # ----------------- AOA SWEEP 2 ------------------------------------------------
    distributions = []

    @time wingpolar2 = pnl.run_polarsweep(ll,
                                magUinf, rho, X0, cref, b;
                                use_Uind_for_force = false,
                                aoa_sweeps = (sweep1, sweep2, sweep3),
                                # sweepname = run_name,
                                # plots_path= save_outputs ? fig_path : nothing,
                                # extraplots_path = save_outputs ? fig_path : nothing,
                                output_distributions=distributions,
                                solver,
                                solver_optargs,
                                align_joints_with_Uinfs, 
                            )

    # --------- Load distribution 2

    plot_distribution(distributions, sweep1; suffix="loading2")

    # --------- Integrated forces and moment: lift, drag, and pitching moment
    plot_polars(wingpolar2; suffix="CLCDCm2", 
                            aoalims=[-17, 17], aoaticks=-15:5:15,
                            CLlims=[-0.75, 1.0],
                            CDlims=[0, 0.06], 
                            Cmlims=[-1.5, 1.5],
                            )


    plot_polars(wingpolar2; suffix="CLCDCm-zoomout2", 
                            aoalims=[-40, 40], aoaticks=-40:20:40,
                            CDlims=[0, 0.4], CDticks=0:0.1:0.4,
                            )
end