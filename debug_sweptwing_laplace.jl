## Localize the source of the PressureLaplace CL blow-up on the swept wing
## (open lifting body, no caps). Reuses the modernized sweptwing.jl wing
## assembly + Bernoulli/Laplace pipeline, then adds:
##
##   1. quantiles of |body.P| for Bernoulli vs Laplace
##   2. top-K panels by |P_L|, with control-point coords
##   3. manufactured-pressure residual ‖L (P_B − P_B[ref]) − b_L‖ / ‖b_L‖
##      (does the Bernoulli pressure even satisfy the assembled operator?)
##   4. VTK dump of both pressure fields for ParaView inspection
##
## Run both gradient_mode options for direct comparison.

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import LinearAlgebra
const _norm = LinearAlgebra.norm
const _dot = LinearAlgebra.dot
const _cross = LinearAlgebra.cross
using Statistics: quantile, mean
using Printf

# ----- geometry / freestream (copy of sweptwing.jl) -----
AOA = 4.2
magVinf = 30.0
Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]
rho = 1.225
b = 98*0.0254
ar = 5.0
tr = 1.0
lambda = 45
airfoil = "airfoil-rae101.csv"
airfoil_path = joinpath(pnl.examples_path, "data")
n_rfl = 8
NDIVS_rfl = [(0.25, n_rfl, 10.0, false),
             (0.50, n_rfl,  1.0, true),
             (0.25, n_rfl, 1/10.0, false)]
NDIVS_span = [(1.0, 30, 20.0, true)]

bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
bodyoptargs = (; CPoffset=1e-14)

println("Generating body...")
body = simplewing(b, ar, tr, 0, 0, lambda, 0;
                  bodytype=bodytype, bodyoptargs=bodyoptargs,
                  airfoil_root=airfoil, airfoil_tip=airfoil,
                  airfoil_path=airfoil_path, rfl_NDIVS=NDIVS_rfl,
                  delim=",", span_NDIVS=NDIVS_span, b_low=-1.0, b_up=1.0)
wake_direction = reshape(Vinf ./ magVinf, :, 1)
for i in eachindex(body.Das)
    body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
end
println("Number of panels: $(body.ncells)")

# ----- solve + Bernoulli (baseline) -----
backend = pnl.DirectBackend()
pnl.apply_freestream!(body, Vinf)
pnl.solve!(body, pnl.Backslash(body); backend)
pnl.calcfield_U!(body, Vinf; backend)
pnl.calcfield_P!(body, magVinf, rho)
pnl.calcfield_F!(body)

Dhat = Vinf/_norm(Vinf); Shat = [0,1,0]; Lhat = _cross(Dhat, Shat)
nondim = 0.5*rho*magVinf^2*b^2/ar
Sref = b^2/ar
c_ref = b/ar
LDS_b = pnl.calcfield_LDS(body, Lhat, Dhat)
CL_b = sign(_dot(LDS_b[:,1], Lhat)) * _norm(LDS_b[:,1]) / nondim
P_B = copy(body.P)

println("\n#===== BERNOULLI =====#")
@printf "CL_B = %.5f   |P_B|: median=%.3g  p90=%.3g  p99=%.3g  max=%.3g\n" CL_b quantile(abs.(P_B), 0.5) quantile(abs.(P_B), 0.9) quantile(abs.(P_B), 0.99) maximum(abs.(P_B))

function run_laplace(mode::Symbol)
    body_l = deepcopy(body)
    body_l.needs_velocity_gradient[] = true
    body_l.velocity .= 0.0
    body_l.velocity_gradient .= 0.0
    pnl.calcfield_U!(body_l, Vinf; backend)
    pnl.influence!((body_l,), (body_l,), backend;
        scalar_potential=false, velocity=false, velocity_gradient=true)
    frames = pnl.ReferenceFrame(body_l)
    normalization = pnl.WingNormalization(rho, Sref, c_ref)
    pl = pnl.PressureLaplace((body_l,), rho;
        reference_panel=1, reference_pressure=0.0, verbose=false,
        gradient_mode=mode)
    fm = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
        correct_kuttacondition=false, verbose=false)
    pl.velocity_dot[1] .= 0.0
    pl((body_l,), (nothing,), frames, Vinf, 0, 1.0)
    fm((body_l,), (nothing,), frames, Vinf, 0, 1.0)
    F = fm.force[:, 1]
    CL = _dot(F, Lhat)

    P_L = copy(body_l.P)
    L = pl.L[1]
    b_L = pl.b[1]
    P_B_shift = P_B .- P_B[1]
    rel_resid = _norm(L * P_B_shift .- b_L) / max(_norm(b_L), eps())

    return (mode=mode, body=body_l, CL=CL, P_L=P_L, L=L, b_L=b_L, rel_resid=rel_resid)
end

println("\n#===== LAPLACE: gradient_mode=:raw_hessian =====#")
r_raw = run_laplace(:raw_hessian)
@printf "CL_L = %.3f   |P_L|: median=%.3g  p90=%.3g  p99=%.3g  max=%.3g\n" r_raw.CL quantile(abs.(r_raw.P_L), 0.5) quantile(abs.(r_raw.P_L), 0.9) quantile(abs.(r_raw.P_L), 0.99) maximum(abs.(r_raw.P_L))
@printf "Bernoulli pressure satisfies (L P_B − b_L) / |b_L| = %.3e\n" r_raw.rel_resid

println("\n#===== LAPLACE: gradient_mode=:surface_velocity =====#")
r_sv = run_laplace(:surface_velocity)
@printf "CL_L = %.3f   |P_L|: median=%.3g  p90=%.3g  p99=%.3g  max=%.3g\n" r_sv.CL quantile(abs.(r_sv.P_L), 0.5) quantile(abs.(r_sv.P_L), 0.9) quantile(abs.(r_sv.P_L), 0.99) maximum(abs.(r_sv.P_L))
@printf "Bernoulli pressure satisfies (L P_B − b_L) / |b_L| = %.3e\n" r_sv.rel_resid

# ----- top-K offenders for each mode -----
function show_top_offenders(label, P, body; K=10)
    println("\n#===== TOP-$K |P| panels: $label =====#")
    idx = sortperm(abs.(P); rev=true)[1:K]
    @printf "%-5s  %-10s  %-10s  %-10s  %-12s\n" "panel" "x" "y" "z" "P"
    for p in idx
        cp = view(body.controlpoints, :, p)
        @printf "%-5d  %-10.4f  %-10.4f  %-10.4f  %-12.3e\n" p cp[1] cp[2] cp[3] P[p]
    end
end

show_top_offenders("Bernoulli", P_B, body)
show_top_offenders("Laplace (raw_hessian)",     r_raw.P_L, r_raw.body)
show_top_offenders("Laplace (surface_velocity)", r_sv.P_L,  r_sv.body)

# ----- VTK dumps for ParaView -----
mkpath("sweptwing_debug")
body_vtk = deepcopy(body)
body_vtk.P .= P_B
pnl.write_vtk(joinpath("sweptwing_debug", "wing_bernoulli"), body_vtk, 0, 0.0; overwrite=true)
pnl.write_vtk(joinpath("sweptwing_debug", "wing_laplace_raw"), r_raw.body, 0, 0.0; overwrite=true)
pnl.write_vtk(joinpath("sweptwing_debug", "wing_laplace_sv"),  r_sv.body,  0, 0.0; overwrite=true)
println("\nWrote VTK files to sweptwing_debug/")

# ----- y-position distribution of top offenders -----
println("\n#===== Edge-of-wing concentration check =====#")
y_max = maximum(abs.(body.controlpoints[2, :]))
for (label, P, body_x) in [("raw_hessian", r_raw.P_L, r_raw.body),
                           ("surface_velocity", r_sv.P_L, r_sv.body)]
    rel_y = abs.(body_x.controlpoints[2, :]) ./ y_max
    near_tip = rel_y .> 0.9
    near_root = rel_y .< 0.1
    @printf "  %-18s  |P| mean (tip 10%%) = %.3e   |P| mean (root 10%%) = %.3e   |P| mean (mid 80%%) = %.3e\n" label mean(abs.(P[near_tip])) mean(abs.(P[near_root])) mean(abs.(P[.!near_tip .& .!near_root]))
end
