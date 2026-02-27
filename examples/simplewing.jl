#=##############################################################################
# DESCRIPTION
    45deg swept-back wing at an angle of attack of 4.2deg. This wing has an
    aspect ratio of 5.0, a RAE 101 airfoil section with 12% thickness, and no
    dihedral, twist, nor taper. This test case matches the experimental setup
    of Weber, J., and Brebner, G., “Low-Speed Tests on 45-deg Swept-Back Wings,
    Part I,” Tech. rep., 1951.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Dec 2022
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import PyPlot as plt

save_path       = "."                      # Where to save outputs
airfoil_path    = joinpath(pnl.examples_path, "data") # Where to find airfoil contours

paraview        = true                          # Whether to visualize with Paraview

# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 8.4                           # (deg) angle of attack
magVinf         = 30.0                          # (m/s) freestream velocity
Vinf            = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)] # Freestream

rho             = 1.225                         # (kg/m^3) air density

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
b               = 98*0.0254                     # (m) span length
ar              = 5.0                           # Aspect ratio b/c_tip
c = b / ar                                      # (m) root chord length
tr              = 1.0                           # Taper ratio c_tip/c_root
twist_root      = 0                             # (deg) twist at root
twist_tip       = 0                             # (deg) twist at tip
lambda          = 0                             # (deg) sweep
gamma           = 0                             # (deg) dihedral
airfoil         = "airfoil-rae101.csv"          # Airfoil contour file

# ----- Chordwise discretization

# NOTE: NDIVS is the number of divisions (panels) in each dimension. This an be
#       either an integer, or an array of tuples as shown below

n_rfl           = 8                             # Control number of chordwise panels
# n_rfl         = 16                            # <-- uncomment this for finer discretization

#           # 0 to 0.25 of the airfoil has `n_rfl` panels at a geometric expansion of 10 that is not central
NDIVS_rfl = [ (0.25, n_rfl,   10.0, false),
            # 0.25 to 0.75 of the airfoil has `n_rfl` panels evenly spaced
              (0.50, n_rfl,    1.0, true),
            # 0.75 to 1.00 of the airfoil has `n_rfl` panels at a geometric expansion of 0.1 that is not central
              (0.25, n_rfl, 1/10.0, false)]

# NOTE: A geometric expansion of 10 that is not central means that the last
#       panel is 10 times larger than the first panel. If central, the
#       middle panel is 10 times larger than the peripheral panels.

# ----- Spanwise discretization
n_span          = 15                            # Number of spanwise panels on each side of the wing
# n_span        = 60                            # <-- uncomment this for finer discretization

NDIVS_span_l    = [(1.0, n_span, 1.0, false)]  # Discretization of left side
NDIVS_span_r    = [(1.0, n_span, 1.0, false)]  # Discretization of right side


# ----------------- GENERATE BODY ----------------------------------------------
println("Generating wing...")

#= NOTE: Here we loft each side of the wing independently. One could also loft
        the entire wing at once from left tip to right tip, but the sweep of the
        wing would lead to an asymmetric discretization with the panels of left
        side side would have a higher aspect ratio than those of the right side.
        To avoid that, instead we loft the left side from left to right, then we
        loft to right side from right to left, and we combine left and right
        sides into a MultiBody that represents the wing.
=#

# Arguments for lofting the left side of the wing
bodyoptargs_l = (;
                    CPoffset=1e-14,                 # Offset control points slightly in the positive normal direction
                )

# Loft left side of the wing from left to right
# kernel = pnl.ConstantSource                     # Kernel type to use
kernel = pnl.ConstantDoublet                     # Kernel type to use
# kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}                     # Kernel type to use

run_name = ""
if kernel <: pnl.ConstantSource
    run_name = "simplewing_source"
elseif kernel <: pnl.ConstantDoublet
    run_name = "simplewing_doublet"
elseif kernel <: pnl.VortexRing
    run_name = "simplewing_vortexring"
elseif kernel == Union{pnl.ConstantSource, pnl.ConstantDoublet}
    run_name = "simplewing_source_doublet"
else
    error("Unknown kernel type!")
end

bodytype = pnl.NonLiftingBody{kernel}    # Elements and wake model
@time wing = pnl.simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                            bodytype=bodytype, bodyoptargs=bodyoptargs_l,
                            airfoil_root=airfoil, airfoil_tip=airfoil,
                            airfoil_path=airfoil_path,
                            rfl_NDIVS=NDIVS_rfl,
                            delim=",",
                            span_NDIVS=NDIVS_span_l,
                            b_low=-1.0, b_up=0.0,
                            verify_spline=false,
                            verify_rflspline=false,
                           )
@show typeof(wing)

if typeof(wing) <: Union{pnl.MultiBody, pnl.AbstractLiftingBody}
    wing.Das .= repeat(Vinf ./ magVinf, 1, wing.nsheddings)
    wing.Dbs .= repeat(Vinf ./ magVinf, 1, wing.nsheddings)
end

# Freestream at every control point
Uinfs = repeat(Vinf, 1, wing.ncells)

println("Number of panels:\t$(wing.ncells)")


# ----------------- CALL SOLVER ------------------------------------------------
println("Solving body...")

# Solve body (panel strengths) giving `Uinfs` as boundary conditions and
# `Das` and `Dbs` as trailing edge rigid wake direction

# # uncomment to use original (unabstracted) solver
# Das = repeat(Vinf ./ magVinf, 1, body.nsheddings)
# Dbs = repeat(Vinf ./ magVinf, 1, body.nsheddings)
# @time pnl.solve(body, Uinfs, Das, Dbs)

# solver = pnl.Backslash(wing; least_squares=false)
# # elprescribe = Tuple{Int,Float64}[] # [(1, 0.0)]   # Prescribe strength of first panel to be zero
# # solver = pnl.KrylovSolver(body;
# #             method=:gmres,
# #             itmax=20,
# #             atol=1e-4,
# #             rtol=1e-4,
# #             elprescribe,
# #             backend=pnl.FastMultipoleBackend(
# #                         expansion_order=7,
# #                         multipole_acceptance=0.4,
# #                         leaf_size=10
# #                     )
# #         )
# pnl.solve2!(wing, Uinfs, solver; elprescribe=Tuple{Int,Float64}[])

G, rhs = pnl.solve(wing, Uinfs)


# ----------------- POST PROCESSING --------------------------------------------
println("Post processing...")

# Calculate surface velocity induced by the wing on itself
# backend = pnl.FastMultipoleBackend(
#                                 expansion_order=7,
#                                 multipole_acceptance=0.4,
#                                 leaf_size=10
#                             )
backend = pnl.DirectBackend()
@time Us = pnl.calcfield_U(wing, wing; backend)

# NOTE: Since the boundary integral equation of the potential flow has a
#       discontinuity at the boundary, we need to add the gradient of the
#       doublet strength to get an accurate surface velocity

# Calculate surface velocity U_∇μ due to the gradient of the doublet strength
if kernel <: Union{pnl.ConstantDoublet, pnl.VortexRing}
    UDeltaGamma = pnl.calcfield_Ugradmu(wing; force_cellTE=false)
    
    # Add both velocities together
    pnl.addfields(wing, "Ugradmu", "U")
end

Us = pnl.get_field(wing, "U")["field_data"] # save for comparison/verification
Us_matrix = zeros(3, length(Us))
for i in 1:length(Us)
    Us_matrix[:, i] = Us[i]
end
Us = Us_matrix
Us_doublet = deepcopy(Us) # save for comparison/verification



# Calculate pressure coefficient
@time Cps = pnl.calcfield_Cp(wing, magVinf)

# Calculate the force of each panel
@time Fs = pnl.calcfield_F(wing, magVinf, rho)

Ftot = sum(Fs, dims=2)
@show Ftot

# check flow tangency
normals = pnl._calc_normals(wing)
Udotn = sum(Us .* normals, dims=1)
resid = maximum(abs.(Udotn))
println("Max flow tangency residual: $resid")

# ----------------- VISUALIZATION ----------------------------------------------
if paraview
    str = save_path*"/"

    # Save wing as a VTK
    str *= pnl.save(wing, run_name; path=save_path)

    # Call Paraview
    # run(`paraview --data=$(str)`)
end

function plot_Cp(wing, spanposs, b, c; clearme=true)

    # slice wing
    spandirection = [0, 1, 0]

    figname = "simplewing"
    fig = plt.figure(figname, figsize=(8,6))
    if clearme
        plt.clf()
        fig.add_subplot(111, xlabel="x/c", ylabel="Cp")
    end
    axs = fig.get_axes()

    stls = ["-", "--", "-.", ":"]

    normals = pnl._calc_normals(wing)
    controlpoints = pnl.calc_controlpoints(wing, normals)
    for (i, spanpos) in enumerate(spanposs)
        points, Cps = pnl.slicefield(wing, controlpoints, "Cp", spanpos*b/2, spandirection, false)
        axs[1].plot(points[1, :] ./ c, Cps, label="2y/b=$(round(spanpos, digits=2))", stls[i])
    end

    axs[1].legend()
end

plot_Cp(wing, [0.5, 0.8, 0.9, 1.0], b, c; clearme=true)