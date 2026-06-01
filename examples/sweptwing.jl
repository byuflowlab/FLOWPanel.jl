#=##############################################################################
# DESCRIPTION
    45deg swept-back wing at an angle of attack of 4.2deg. This wing has an
    aspect ratio of 5.0, a RAE 101 airfoil section with 12% thickness, and no
    dihedral, twist, nor taper. This test case matches the experimental setup
    of Weber, J., and Brebner, G., "Low-Speed Tests on 45-deg Swept-Back Wings,
    Part I," Tech. rep., 1951.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Dec 2022
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using GeometricTools: PyPlot as plt
import LinearAlgebra
const _norm  = LinearAlgebra.norm
const _dot   = LinearAlgebra.dot
const _cross = LinearAlgebra.cross

run_name        = "sweptwing000"                # Name of this run

save_path       = joinpath("data", run_name)    # Where to save outputs
airfoil_path    = joinpath(pnl.examples_path, "data") # Where to find airfoil contours

paraview        = true                         # Whether to visualize with Paraview

# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 4.2                           # (deg) angle of attack
magVinf         = 30.0                          # (m/s) freestream velocity
Vinf            = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)] # Freestream

rho             = 1.225                         # (kg/m^3) air density

# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
b               = 98*0.0254                     # (m) span length
ar              = 5.0                           # Aspect ratio b/c_tip
tr              = 1.0                           # Taper ratio c_tip/c_root
twist_root      = 0                             # (deg) twist at root
twist_tip       = 0                             # (deg) twist at tip
lambda          = 45                            # (deg) sweep
gamma           = 0                             # (deg) dihedral
airfoil         = "airfoil-rae101.csv"          # Airfoil contour file

# ----- Chordwise discretization
n_rfl           = 8                             # Control number of chordwise panels
NDIVS_rfl = [ (0.25, n_rfl,   10.0, false),
              (0.50, n_rfl,    1.0, true),
              (0.25, n_rfl, 1/10.0, false)]

# ----- Spanwise discretization (full span, single loft)
# Uniform distribution: with `central=true, expansion=20` the *root* panels
# end up coarsest and the tips finest, which makes the inner Cp slice
# (2y/b ≈ 0.04) too under-resolved to plot cleanly. Uniform spacing keeps
# panel size constant across the span.
n_span_full     = 40                            # Number of spanwise panels across full span
NDIVS_span      = [(1.0, n_span_full, 1.0, true)]


# ----------------- GENERATE BODY ----------------------------------------------
println("Generating body...")

bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}

bodyoptargs = (; cp_outer=true)

@time body = simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
                        bodytype=bodytype, bodyoptargs=bodyoptargs,
                        airfoil_root=airfoil, airfoil_tip=airfoil,
                        airfoil_path=airfoil_path,
                        rfl_NDIVS=NDIVS_rfl,
                        delim=",",
                        span_NDIVS=NDIVS_span,
                        b_low=-1.0, b_up=1.0,
                       )
wake_direction = reshape(Vinf ./ magVinf, :, 1)
for i in eachindex(body.Das)
    body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
end

# Freestream at every control point
Uinfs = repeat(Vinf, 1, body.ncells)

println("Number of panels:\t$(body.ncells)")


# ----------------- CALL SOLVER ------------------------------------------------
println("Solving body...")

backend = pnl.DirectBackend()

pnl.apply_freestream!(body, Vinf)
solver = pnl.Backslash(body)
@time pnl.solve!(body, solver; backend)


# ----------------- POST PROCESSING: BERNOULLI ---------------------------------
println("Post processing (Bernoulli)...")

@time pnl.calcfield_U!(body, Vinf; backend)

# Bernoulli pressure (`calcfield_P!` evaluates 0.5 ρ (|U∞|² − |u|²); the
# ∇_s µ correction is already folded into body.velocity by calcfield_U!)
@time pnl.calcfield_P!(body, magVinf, rho)

# Per-panel force from Bernoulli pressure
@time pnl.calcfield_F!(body)


# ----------------- INTEGRATED FORCES: BERNOULLI -------------------------------
Dhat = Vinf/_norm(Vinf)        # Drag direction
Shat = [0, 1, 0]               # Span direction
Lhat = _cross(Dhat, Shat)      # Lift direction

LDS_bernoulli = pnl.calcfield_LDS(body, Lhat, Dhat)
L_b = LDS_bernoulli[:, 1]
D_b = LDS_bernoulli[:, 2]

nondim = 0.5*rho*magVinf^2*b^2/ar   # Normalization factor (Sref = b^2/ar)
CL_bernoulli = sign(_dot(L_b, Lhat)) * _norm(L_b) / nondim
CD_bernoulli = sign(_dot(D_b, Dhat)) * _norm(D_b) / nondim


# ----------------- POST PROCESSING: LAPLACE -----------------------------------
println("Post processing (Laplace)...")

body_l = deepcopy(body)
body_l.needs_velocity_gradient[] = true
body_l.velocity .= 0.0
body_l.velocity_gradient .= 0.0
pnl.calcfield_U!(body_l, Vinf; backend)
pnl.influence!((body_l,), (body_l,), backend;
    scalar_potential=false, velocity=false, velocity_gradient=true)

Sref = b^2 / ar
c_ref = b / ar
normalization = pnl.WingNormalization(rho, Sref, c_ref)

frames_l = pnl.ReferenceFrame(body_l)
pressure_laplace = pnl.PressureLaplace((body_l,), rho;
    reference_panel=1, reference_pressure=0.0, verbose=false,
    unsteady=false,
    gradient_mode=:surface_velocity)
force_laplace = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
    correct_kuttacondition=false, verbose=false)

pressure_laplace.velocity_dot[1] .= 0.0
@time pressure_laplace((body_l,), (nothing,), frames_l, Vinf, 0, 1.0)
@time force_laplace((body_l,), (nothing,), frames_l, Vinf, 0, 1.0)

# `WingNormalization` already divides by 0.5 ρ |Vinf|² Sref, so force_laplace.force
# is in coefficient form. Project onto Lhat/Dhat to recover CL/CD.
F_lap = force_laplace.force[:, 1]
CL_laplace = _dot(F_lap, Lhat)
CD_laplace = _dot(F_lap, Dhat)


# ----------------- VISUALIZATION ----------------------------------------------
if paraview
    mkpath(save_path)
    pnl.write_vtk(joinpath(save_path, run_name*"_bernoulli_AOA$(AOA)"),
                  body, 0, 0.0; overwrite=true)
    pnl.write_vtk(joinpath(save_path, run_name*"_laplace_AOA$(AOA)"),
                  body_l, 0, 0.0; overwrite=true)
    println("Wrote VTK files to $(save_path)/")
end


# ----------------- COMPARISON TO EXPERIMENTAL DATA ----------------------------
include(joinpath(pnl.examples_path, "sweptwing_postprocessing.jl"))

save_outputs = false

fig_path = joinpath(pnl.examples_path, "..", "docs", "resources", "images")
outdata_path = joinpath(pnl.examples_path, "..", "docs", "resources", "data")

# `plot_Cps` uses the unstructured-mesh-compatible `slice_scalarfield` primitive
# and works on this body. `plot_deltaCps` / `plot_loading` still depend on
# structured-grid helpers (`slicefield` / `calcfield_sectionalforce`) and remain
# disabled until those are ported.
make_plots_cps = true

# --------- Summary (printed before plotting so we get numbers even if plots fail)
CLexp = CLs_web[2]
CDexp = CDs_web[2]

println("\n#===== INTEGRATED CL/CD =====#")
@show CL_bernoulli CL_laplace CLexp
@show CD_bernoulli CD_laplace CDexp
println("Bernoulli vs Laplace CL diff: $(round(abs(CL_laplace-CL_bernoulli), sigdigits=4))")
println("Bernoulli vs Laplace CD diff: $(round(abs(CD_laplace-CD_bernoulli), sigdigits=4))")
println("Bernoulli CL error: $(round(abs(CL_bernoulli-CLexp)/CLexp*100, digits=2))%")
println("Laplace   CL error: $(round(abs(CL_laplace-CLexp)/CLexp*100, digits=2))%")

if make_plots_cps
    side = 1
    spanposs_cps = side*parse.(Float64, keys(weber_Cps["$AOA"]))[[2, 4, 5, 7]]
    # 45° sweep: LE x at spanwise position y is |y|*tan(λ).
    xLE_fn = y -> abs(y) * tan(lambda * pi / 180)
    fig1, axs = plot_Cps(body, spanposs_cps, b, rho, magVinf;
                                xscaling=ar/b, AOA=AOA,
                                xlims=[-0.1, 1.1], ylims=[1.0, -0.7], stl="-",
                                slicetol=0.013*b, xLE_fn=xLE_fn)
    fig1.tight_layout()
    fig1.savefig(joinpath(@__DIR__, "..", "sweptwing_Cps.png"), dpi=150)
    println("Saved Cp plot to sweptwing_Cps.png")
end
