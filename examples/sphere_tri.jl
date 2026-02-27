#=##############################################################################
# DESCRIPTION
    Example of inviscid flow around a sphere

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Aug 2022
  * License   : MIT License
=###############################################################################


import FLOWPanel as pnl

import PyPlot as plt
import CSV
import DataFrames: DataFrame
import PyPlot: @L_str
import GeoIO

save_path       = "temps/"                  # Where to save results
file_name       = "sphere00"                # Prefix of output files
paraview        = true                      # Whether to visualize results in paraview

# -------- Simulation Parameters -----------------------------------------------
nu              = 1.443e-5                  # (m^2/s) kinematic viscosity
Re              = 8800                      # Reynolds number V*d/nu
R               = 1                         # (m) radius of sphere
magVinf         = Re*nu/(2*R)               # (m/s) freestream velocity
magVinf = 1.0
Vinf            = magVinf*[1.0,0,0]         # (m/s) freestream vector

# Solver options
# kernel    = Union{pnl.ConstantSource} # Elements to use
kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}

# -------- Generate Geometry ---------------------------------------------------
# meshfile = "mysphere.msh"
meshfile = "sphere_r1.msh"
msh = GeoIO.load(meshfile).geometry

# Generates grid
grid = pnl.gt.GridTriangleSurface(msh)

# Creates non lifting body
body = pnl.NonLiftingBody{kernel}(grid)


# -------- Call Solver ------------------------------------------------------------
# Freestream at every control point
Vinfs = repeat(Vinf, 1, body.ncells)

# Solve body
@time pnl.solve(body, Vinfs; solve_type=3)

# -------- Postprocess -------------------------------------------------------------
# induced velocity
normals = pnl._calc_normals(body)
ctrlpts = pnl._calc_controlpoints(body, normals; off=1e-8)
U = zeros(3, size(ctrlpts, 2))
U = pnl.calcfield_U!(U, body, body, ctrlpts, Vinfs)

# contribution due to gradient of doublet strength
Ugradmu = pnl.calcfield_Ugradmu(body; off=1e-8, force_cellTE=false, sharpTE=false, Gammai=2)

# total velocity
pnl.addfields(body, "Ugradmu", "U")

# pressure/force coefficients
Cps = pnl.calcfield_Cp(body, magVinf)
Fs = pnl.calcfield_F(body, magVinf, 1.225) # rho


# -------- Visualization ----------------------------------------------------------
str = save_path*"/"

# Save vtk
str *= pnl.save(body, file_name; path=save_path, debug=true)

# Call Paraview
# if paraview
#     run(`paraview --data=$str`)
# end


# -------- Verification and Validation -----------------------------------------
# Compare predicted pressure coefficient to analytic and experimental results
fig = plt.figure(figsize=[7, 5]*2/3)
ax = fig.gca()

# Invicid potential flow
theo_theta = range(0, stop=pi, length=180)
theo_Cp = 1 .- 9 / 4 * (sin.(theo_theta)) .^ 2

ax.plot(theo_theta * 180 / pi, theo_Cp, "-k", label="Analytic (inviscid)")

# FLOWPanel (inviscid)
# dim_plot = 2                                # Plot along this coordinate dimension
# ncells = body.grid._ndivscells[dim_plot]    # Number of panels along this dimension
# dtheta = (P_max[dim_plot]-P_min[dim_plot])/ncells
# pnl_theta = [i * dtheta - dtheta / 2 for i in 1:ncells] # Angular position of each panel

# if dim_plot==1
#     dim2 = Int(floor(body.grid._ndivscells[2] / 2))
#     pnl_U = [pnl.get_fieldval(body, "U", [i, dim2]) for i in 1:ncells]
# else
#     dim1 = Int(floor(body.grid._ndivscells[1] / 2))
#     pnl_U1 = [pnl.get_fieldval(body, "U", [dim1-1, i]) for i in 1:ncells]
#     pnl_U2 = [pnl.get_fieldval(body, "U", [dim1, i]) for i in 1:ncells]

#     # Average the triangular panels in this dimension, otherwise the pressure
#     # is jumpy
#     pnl_U = (pnl_U1 .+ pnl_U2)/2
# end

# pnl_Cp = [1 - (pnl.norm(U) / magVinf)^2 for U in pnl_U]
# ax.plot(pnl_theta * 180 / pi, pnl_Cp, "or", alpha=0.8, label="FLOWPanel (inviscid)")


# # Experimental
# exp_theta = [-99.99939258190653, -93.24452374625109, -88.5321741771383,
#             -83.8143578451843, -79.68148513723852, -73.51976006985308,
#             -67.68399073687408, -63.918074484643704, -58.70961618769218,
#             -53.80524657378232, -48.6077218025132, -43.71565240499601,
#             -37.96871796818648, -34.242435746554804, -28.76542272502941,
#             -24.737785201776674, -18.944383280816965, -13.720891386052173,
#             -9.06457613606166, -3.8028169014084483, 1.7630310162863907,
#             7.645267833415602, 13.234349493185533, 17.67404426559358,
#             23.305493337382813, 28.353365475874085, 33.147716487604896,
#             36.765346797767734, 41.53783075813374, 47.19524695341863,
#             51.39918757829997, 56.19080520861016, 61.82498766181999,
#             66.2701491970692, 71.58520936942412, 77.75650127178167,
#             82.7169052048138, 87.67457575642538, 92.04183592118753,
#             97.29402832086862]
# exp_Cp = [-0.44314946281462353, -0.4639535325158495, -0.5011275198359972,
#             -0.5476253748908546, -0.5964693823317264, -0.6056338028169015,
#             -0.558862609619984, -0.48184199536843697, -0.3651569796135299,
#             -0.22983182111537137, -0.0944990698910444, 0.06180479101021241,
#             0.2600888348961694, 0.40470749022436514, 0.5633575035116358,
#             0.6939979499639347, 0.8130291181048556, 0.9040734975893094,
#             0.9624691545499413, 0.988246459891424, 0.9953836224896551,
#             0.962901939941536, 0.9304126646672488, 0.8582665806157701,
#             0.7535173303974795, 0.6440909608594968, 0.4670589575186972,
#             0.2969894840742573, 0.15725295167229825, 0.0082153297141343,
#             -0.16183895827796912, -0.3342090277514138, -0.4436202118370598,
#             -0.5250901636232483, -0.5902205686951896, -0.6157017577161077,
#             -0.5759462434987279, -0.5315287954139931, -0.48013363198056225,
#             -0.4380395581033367]
# ax.plot(exp_theta, exp_Cp, "^k",
#             alpha=0.5, markersize=3, label="Experimental (viscous)")


# fig.suptitle("Pressure distribution over sphere")
# ax.set_ylim([1.5, -1.5])
# ax.set_xlim([0, 180])
# ax.set_xticks(0:45:180)
# ax.set_xlabel(L"Angle $\theta$ ($^\circ$)")
# ax.set_ylabel(L"$C_p=\frac{p-p_\infty}{\frac{1}{2}\rho U_\infty^2}$")
# ax.spines["right"].set_visible(false)
# ax.spines["top"].set_visible(false)
# ax.legend(loc="best", frameon=false)

fig.tight_layout()
# fig.savefig(joinpath(save_path, file_name*"-Cp.png"), dpi=300)
