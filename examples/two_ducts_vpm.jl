#=##############################################################################
# DESCRIPTION
    Two-duct unsteady simulation using the tuple-based simulate! interface.
    Matches the geometry of two_ducts.jl: body1 is a lifting duct, body2 is a
    non-lifting duct offset in z by 1.5 diameters.

# AUTHORSHIP
  * Author    : Ryan Anderson
  * Created   : Apr 2026
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import CSV
import DataFrames: DataFrame
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: norm

# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 5                             # (deg) angle of attack
magVinf         = 30.0                          # (m/s) freestream velocity
rho             = 1.225                         # (kg/m^3) air density
Vinf            = magVinf * [cosd(AOA), 0.0, sind(AOA)]

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
filename        = joinpath(pnl.examples_path, "data", "naca662015.csv")
contour         = CSV.read(filename, DataFrame)

aspectratio     = 0.6                           # Duct trailing edge aspect ratio l/d
d               = 2 * 0.835                     # (m) duct diameter

# ----------------- SOLVER PARAMETERS ------------------------------------------
NDIVS_theta     = 20                            # Number of azimuthal panels
n_rfl           = 6                             # Chordwise panel control

NDIVS_rfl = [(0.25, n_rfl, 10.0, false),
             (0.50, n_rfl,  1.0, true),
             (0.25, n_rfl,  0.1, false)]

# ----------------- HELPER: GENERATE DUCTS ------------------------------------
function make_lifting_duct(; semiinfinite_wake=false)
    kernel   = Union{pnl.ConstantSource, pnl.VortexRing}
    bodytype = pnl.RigidWakeBody{kernel}

    xs, ys = pnl.gt.rediscretize_airfoil(contour[:, 1], contour[:, 2],
                                          NDIVS_rfl, NDIVS_rfl; verify_spline=false)
    ys[end] = ys[1]
    xs *= d * aspectratio
    ys *= d * aspectratio
    ys .+= d / 2
    points = hcat(xs, ys)

    return pnl.generate_revolution_liftbody(bodytype, points, NDIVS_theta;
                bodyoptargs=(CPoffset=1e-10, kerneloffset=1e-2,
                             kernelcutoff=1e-14,
                             characteristiclength=(args...)->d*aspectratio,
                             semiinfinite_wake))
end

function make_nonlifting_duct()
    kernel   = pnl.ConstantSource
    bodytype = pnl.NonLiftingBody{kernel}

    xs, ys = pnl.gt.rediscretize_airfoil(contour[:, 1], contour[:, 2],
                                          NDIVS_rfl, NDIVS_rfl; verify_spline=false)
    ys[end] = ys[1]
    xs *= d * aspectratio
    ys *= d * aspectratio
    ys .+= d / 2
    points = hcat(xs, ys)

    return pnl.generate_revolution_liftbody(bodytype, points, NDIVS_theta;
                bodyoptargs=(CPoffset=1e-10, kerneloffset=1e-2,
                             kernelcutoff=1e-14,
                             characteristiclength=(args...)->d*aspectratio))
end

# ----------------- GENERATE TWO DUCTS ----------------------------------------
body1 = make_lifting_duct()
# body2 = make_nonlifting_duct()
body2 = make_lifting_duct()

# Offset second duct in z by 1.5 diameters (matching two_ducts.jl)
pnl.rotate!(body2, 0, 0, 0; translation=[0.0, 0.0, 1.5*d])

println("Body 1 panels: $(body1.ncells)")
println("Body 2 panels: $(body2.ncells)")

# ----------------- STEADY-STATE REFERENCE SOLVE ------------------------------
# Solve the same two-duct geometry with semi-infinite wake steady-state method
# to provide a reference for the unsteady simulation comparison.
println("\nSolving steady-state reference...")

ref1 = make_lifting_duct(; semiinfinite_wake=true)
ref2 = make_lifting_duct(; semiinfinite_wake=true)
pnl.rotate!(ref2, 0, 0, 0; translation=[0.0, 0.0, 1.5*d])

# Set wake directions
ref1.Das[1] .= repeat(Vinf / magVinf, 1, size(ref1.Das[1], 2))
ref2.Das[1] .= repeat(Vinf / magVinf, 1, size(ref2.Das[1], 2))

# Apply freestream
pnl.apply_freestream!(ref1, Vinf)
pnl.apply_freestream!(ref2, Vinf)

# Solve with BackslashDirichlet (multi-body outer iteration)
ref_solver1 = pnl.BackslashDirichlet(ref1)
ref_solver2 = pnl.BackslashDirichlet(ref2)

@time pnl.solve!((ref1, ref2), (ref_solver1, ref_solver2);
    backend=fill(pnl.FastMultipoleBackend(), 2),
    max_outer_iterations=50,
    outer_tolerance=1e-8,
    verbose=true)

# Post-process reference solution
pnl.calcfield_U!((ref1, ref2), Vinf; backend=pnl.FastMultipoleBackend())
pnl.calcfield_Cp!((ref1, ref2), magVinf)
pnl.calcfield_F!((ref1, ref2), magVinf, rho)

F1_ref = sum(ref1.F, dims=2)
F2_ref = sum(ref2.F, dims=2)
println("Steady-state Body 1 force: ", F1_ref)
println("Steady-state Body 2 force: ", F2_ref)

# Flow tangency diagnostic (must use tuple form to include cross-body influence)
println("\n--- Flow Tangency ---")
pnl.calcfield_U!((ref1, ref2), Vinf; backend=pnl.FastMultipoleBackend())
for (i, body) in enumerate((ref1, ref2))
    Udotn = sum(body.velocity .* body.normals, dims=1)
    rms = sqrt(sum(Udotn .^ 2) / body.ncells)
    maxr = maximum(abs.(Udotn))
    println("  Body $i: RMS(U·n) = $(round(rms; digits=4)),  max|U·n| = $(round(maxr; digits=4)),  RMS/Vinf = $(round(rms/magVinf*100; digits=2))%")
end

# write VTK for reference solution
pnl.write_vtk("two_ducts_vpm/reference_solution_1", ref1, 0, 0.0)
pnl.write_vtk("two_ducts_vpm/reference_solution_2", ref2, 0, 0.0)

println()

# ----------------- SIMULATION SETUP -------------------------------------------
Uinf(t) = Vinf

l       = d * aspectratio
dt_val  = l / magVinf / 5                  # ~0.2 chord per step
n_steps = 21
t_range = range(0.0, step=dt_val, length=n_steps)

# Set wake directions for lifting body (unit direction, scaled by Das offset)
body1.Das[1] .= repeat(Vinf / magVinf * 0.05, 1, size(body1.Das[1], 2))
if typeof(body2) <: pnl.AbstractLiftingBody
    body2.Das[1] .= repeat(Vinf / magVinf * 0.05, 1, size(body2.Das[1], 2))
end

# Create wake for lifting body only; non-lifting body has no wake
wake1 = pnl.PanelWake(body1; nwakerows=100)
wake2 = typeof(body2) <: pnl.NonLiftingBody ? nothing : pnl.PanelWake(body2; nwakerows=100)

# Backend
backend = pnl.FastMultipoleBackend(;
    expansion_order=10,
    multipole_acceptance=0.4,
    leaf_size=20
)

# Solvers (one per body)
solver1 = pnl.Backslash(body1)
solver2 = pnl.Backslash(body2)

# Reference frames (static)
frames = pnl.ReferenceFrame(body1;
    origin = SVector{3}(0.0, 0.0, 0.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 1.0, 0.0),
    ω = 0.0,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = [1, 2]
)

# Maneuver (no-op for static frames)
maneuver!(frames, systems, wakes, t) = nothing

# ----------------- RUN SIMULATION ---------------------------------------------
systems = (body1, body2)
wakes   = (wake1, wake2)
body_solvers = (solver1, solver2)

println("\nBegin two-duct simulation ($(n_steps) steps)...")
@time pnl.simulate!(systems, wakes, frames, maneuver!, Uinf, t_range;
    body_solvers, backend, rho, verbose=true,
    path="two_ducts_vpm", name="two_ducts"
)

# ----------------- REPORT RESULTS ---------------------------------------------
println("\n=== RESULTS ===")
F1 = sum(body1.F, dims=2)
F2 = sum(body2.F, dims=2)
println("Body 1 total force: ", F1)
println("Body 2 total force: ", F2)
println("Body 1 max |Cp|: ", maximum(abs.(body1.Cp)))
println("Body 2 max |Cp|: ", maximum(abs.(body2.Cp)))

# Flow tangency diagnostic (must use tuple form to include cross-body influence)
println("\n--- Flow Tangency ---")
pnl.calcfield_U!(systems, Vinf, wakes; backend=pnl.FastMultipoleBackend())
for (i, body) in enumerate(systems)
    Udotn = sum(body.velocity .* body.normals, dims=1)
    rms = sqrt(sum(Udotn .^ 2) / body.ncells)
    maxr = maximum(abs.(Udotn))
    println("  Body $i: RMS(U·n) = $(round(rms; digits=4)),  max|U·n| = $(round(maxr; digits=4)),  RMS/Vinf = $(round(rms/magVinf*100; digits=2))%")
end

println("\n--- Comparison to steady-state ---")
for (i, (f, fref)) in enumerate(zip([F1, F2], [F1_ref, F2_ref]))
    fvec = [f[1], f[2], f[3]]
    frefvec = [fref[1], fref[2], fref[3]]
    err = norm(fvec - frefvec) / norm(frefvec)
    println("  Body $i relative force error: $(round(err*100; digits=2))%")
end
