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
plt.rc("font", family="STIXGeneral")            # Text font
plt.rc("mathtext", fontset="stix")                  # Math font
plt.rc("font", size=12)          # controls default text sizes
plt.rc("axes", titlesize=12)     # fontsize of the axes title
plt.rc("axes", labelsize=14)     # fontsize of the x and y labels
plt.rc("xtick", labelsize=12)    # fontsize of the tick labels
plt.rc("ytick", labelsize=12)    # fontsize of the tick labels
plt.rc("legend", fontsize=12)    # legend fontsize
plt.rc("figure", titlesize=18)   # fontsize of the figure title


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
#    2*y/b  polar file            airfoil type
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
                    maxiters = 400,
                    )

align_joints_with_Uinfs = false                 # Whether to align joint bound vortices with the freestream


use_Uind_for_force = true                       # Whether to use Uind as opposed to selfUind for force postprocessing

Dhat            = Uinf/norm(Uinf)               # Drag direction
Shat            = [0, 1, 0]                     # Span direction
Lhat            = cross(Dhat, Shat)             # Lift direction

X0              = [0.0 * chord_distribution[1, 2]*b, 0, 0] # Center about which to calculate moments
lhat            = Dhat                          # Rolling direction
mhat            = Shat                          # Pitching direction
nhat            = Lhat                          # Yawing direction

cref            = chord_distribution[1, 2]*b    # Reference chord
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

save_outputs = true                        # Whether to save outputs or not

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

ax.plot(ypos_web, cls_web[2, :], ":ok", label=lbl_web)
ax.plot(-ypos_web, cls_web[2, :], ":ok")

ax.plot(ypos_vsp, cls_vsp, "-", label=lbl_vsp, color=color_vsp, alpha=alpha_vsp)

ax.plot(ypos, cls, "-", label=lbl_ll, color=color_ll, alpha=0.9)

ax.set_ylabel(L"Sectional lift $c_\ell$")
ax.set_ylim([0, 0.4])
ax.set_yticks(0:0.1:0.4)

ax = axs[2]

ax.plot(ypos_web, cds_web[3, :], ":ok", label=lbl_web)
ax.plot(-ypos_web, cds_web[3, :], ":ok")

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
    ax.set_xlabel(L"Spanwise position $2y/b$")
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