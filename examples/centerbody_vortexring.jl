bodytype = pnl.RigidWakeBody{pnl.VortexRing}    # Elements and wake model

body = pnl.generate_revolution_liftbody(bodytype, points, NDIVS_theta;
                                        # Loop the azimuthal dimension to close the surface
                                        loop_dim=2,
                                        # Rotate the axis of rotation to align with x-axis
                                        axis_angle=90,
                                        # Indicate that this body is open at the trailing edge
                                        closed_contour=false)


# ----------------- CALL SOLVER ------------------------------------------------
println("Solving body...")

# Freestream at every control point
Uinfs = repeat(Vinf, 1, body.ncells)

# Unitary direction of semi-infinite vortex at points `a` and `b` of each
# trailing edge panel
# NOTE: In this case they are empty arrays since there is no wake
Das = repeat(Vinf/magVinf, 1, body.nsheddings+1)

# Solve body (panel strengths) giving `Uinfs` as boundary conditions and
# `Das` as trailing edge rigid wake direction
body.velocity .= Uinfs
body.Das .= Das
solver = pnl.Backslash(body)
@time pnl.solve!(body, solver)


# ----------------- POST PROCESSING --------------------------------------------
println("Post processing...")
# NOTE: A thick body with only vortex ring elements leads to a surface velocity
#       that is inaccurate at the exact surface of the body, but that
#       approximates the physical solution away from the surface. For this
#       reason, we probe the velocity used to calculate Cp slightly away from
#       the body

# Calculate surface velocity on the body
Us = pnl.calcfield_U(body, body)

# NOTE: Since the boundary integral equation of the potential flow has a
#       discontinuity at the boundary, we need to add the gradient of the
#       doublet strength to get an accurate surface velocity

# Calculate surface velocity U_∇μ due to the gradient of the doublet strength
UDeltaGamma = pnl.calcfield_Ugradmu(body)

# Add both velocities together
pnl.addfields(body, "Ugradmu", "U")

# Calculate gauge pressure
Ps = pnl.calcfield_P(body, magVinf, rho)

# Calculate the force of each panel
Fs = pnl.calcfield_F(body)


# ----------------- COMPARISON TO EXPERIMENTAL DATA ----------------------------
include(joinpath(pnl.examples_path, "centerbody_postprocessing.jl"))

# Where to save figures (default to re-generating the figures that are used
# in the docs)
fig_path = joinpath(pnl.examples_path, "..", "docs", "resources", "images")
outdata_path = joinpath(pnl.examples_path, "..", "docs", "resources", "data")

if save_outputs
    fig.savefig(joinpath(fig_path, "$(run_name)-velocity-vortexring.png"),
                                                dpi=300, transparent=true)
end
