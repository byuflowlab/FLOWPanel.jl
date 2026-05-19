## Hypothesis 1 from the analysis: the surface-LS gradient reconstruction
## itself is amplifying small velocity errors. Test by comparing:
##   - the analytic FMM Hessian (`body.velocity_gradient`) used by :raw_hessian
##   - the surface-LS reconstruction (`m.surface_velocity_gradient[1]`) used
##     by :surface_velocity
## against the expected physical scale `~ |Vinf| / c_ref`.
##
## If `:surface_velocity` is amplifying by ~50× relative to the Hessian or
## the physical scale, that explains the uniform 50× over-prediction of
## both median |P_L| and integrated CL_L.

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import LinearAlgebra
const _norm = LinearAlgebra.norm
using Statistics: quantile, mean, median
using Printf

# Wing assembly (copy of sweptwing.jl)
AOA = 4.2; magVinf = 30.0
Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]
rho = 1.225
b = 98*0.0254; ar = 5.0; tr = 1.0; lambda = 45
airfoil = "airfoil-rae101.csv"
airfoil_path = joinpath(pnl.examples_path, "data")
n_rfl = 8
NDIVS_rfl = [(0.25, n_rfl, 10.0, false),
             (0.50, n_rfl,  1.0, true),
             (0.25, n_rfl, 1/10.0, false)]
NDIVS_span = [(1.0, 30, 20.0, true)]
bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
bodyoptargs = (; CPoffset=1e-14)

c_ref = b / ar
println("Physical scale: |Vinf|/c_ref = $(magVinf/c_ref) /s")

body = simplewing(b, ar, tr, 0, 0, lambda, 0;
                  bodytype=bodytype, bodyoptargs=bodyoptargs,
                  airfoil_root=airfoil, airfoil_tip=airfoil,
                  airfoil_path=airfoil_path, rfl_NDIVS=NDIVS_rfl,
                  delim=",", span_NDIVS=NDIVS_span, b_low=-1.0, b_up=1.0)
wake_direction = reshape(Vinf ./ magVinf, :, 1)
for i in eachindex(body.Das)
    body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
end
println("npanels = $(body.ncells)")

backend = pnl.DirectBackend()
pnl.apply_freestream!(body, Vinf)
pnl.solve!(body, pnl.Backslash(body); backend)

# Run both modes; capture the gradient tensors
body_r = deepcopy(body)
body_r.needs_velocity_gradient[] = true
body_r.velocity .= 0.0; body_r.velocity_gradient .= 0.0
pnl.calcfield_U!(body_r, Vinf; backend)
pnl.influence!((body_r,), (body_r,), backend;
    scalar_potential=false, velocity=false, velocity_gradient=true)
G_raw = copy(body_r.velocity_gradient)   # 3 × 3 × ncells

body_s = deepcopy(body)
body_s.needs_velocity_gradient[] = true
body_s.velocity .= 0.0; body_s.velocity_gradient .= 0.0
pnl.calcfield_U!(body_s, Vinf; backend)
pnl.influence!((body_s,), (body_s,), backend;
    scalar_potential=false, velocity=false, velocity_gradient=true)
frames_s = pnl.ReferenceFrame(body_s)
pl_s = pnl.PressureLaplace((body_s,), rho;
    reference_panel=1, reference_pressure=0.0, verbose=false,
    gradient_mode=:surface_velocity)
fm_s = pnl.ForceMonitor(1, 1; i_frame=-1,
    normalization=pnl.WingNormalization(rho, b^2/ar, c_ref),
    correct_kuttacondition=false, verbose=false)
pl_s.velocity_dot[1] .= 0.0
pl_s((body_s,), (nothing,), frames_s, Vinf, 0, 1.0)
G_surf = copy(pl_s.surface_velocity_gradient[1])

# Frobenius norms per panel
function frob_per_panel(G)
    np = size(G, 3)
    return [sqrt(sum(G[:, :, p].^2)) for p in 1:np]
end

f_raw  = frob_per_panel(G_raw)
f_surf = frob_per_panel(G_surf)

println("\n#===== ‖∇u‖_F per-panel quantiles =====#")
@printf "%-22s  %-10s  %-10s  %-10s  %-10s\n" "source" "median" "p90" "p99" "max"
@printf "%-22s  %-10.3g  %-10.3g  %-10.3g  %-10.3g\n" "raw_hessian (FMM)"     median(f_raw)  quantile(f_raw, 0.9)  quantile(f_raw, 0.99)  maximum(f_raw)
@printf "%-22s  %-10.3g  %-10.3g  %-10.3g  %-10.3g\n" "surface_velocity (LS)" median(f_surf) quantile(f_surf, 0.9) quantile(f_surf, 0.99) maximum(f_surf)
@printf "%-22s  %-10.3g\n" "physical scale Vinf/c" magVinf/c_ref

ratio = median(f_surf) / median(f_raw)
@printf "\nmedian(LS) / median(Hessian) = %.2f\n" ratio
@printf "median(LS) / (Vinf/c)        = %.2f\n" median(f_surf) / (magVinf/c_ref)
@printf "median(Hessian) / (Vinf/c)   = %.2f\n" median(f_raw)  / (magVinf/c_ref)
