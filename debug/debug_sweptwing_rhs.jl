## Hypothesis 2: a uniform scaling / sign factor in `_pressure_rhs!`
## (or its inputs) inflates the Laplace pressure by ~45×.
##
## Test logic: Bernoulli gives the physically-correct pressure field P_B.
## If the assembled L is correct and only b has a uniform scale error S,
## then L (P_B - P_B[ref]) and b should be related by `b ≈ S * L P_B`
## element-by-element. The distribution of the ratio tells us:
##   - tight cluster around a single S  → uniform scaling bug; S is the factor
##   - broad / multimodal               → not a simple scaling; bug is local

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import LinearAlgebra
const _norm = LinearAlgebra.norm
const _dot  = LinearAlgebra.dot
const _cross = LinearAlgebra.cross
using Statistics: quantile, median, mean, std
using Printf

# wing assembly (same as sweptwing.jl)
AOA = 4.2; magVinf = 30.0
Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]
rho = 1.225
b = 98*0.0254; ar = 5.0
airfoil = "airfoil-rae101.csv"
airfoil_path = joinpath(pnl.examples_path, "data")
n_rfl = 8
NDIVS_rfl = [(0.25, n_rfl, 10.0, false),
             (0.50, n_rfl,  1.0, true),
             (0.25, n_rfl, 1/10.0, false)]
NDIVS_span = [(1.0, 30, 20.0, true)]
bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
c_ref = b / ar
Sref = b^2 / ar

body = simplewing(b, ar, 1.0, 0, 0, 45, 0;
                  bodytype=bodytype, bodyoptargs=(; cp_outer=true),
                  airfoil_root=airfoil, airfoil_tip=airfoil,
                  airfoil_path=airfoil_path, rfl_NDIVS=NDIVS_rfl,
                  delim=",", span_NDIVS=NDIVS_span, b_low=-1.0, b_up=1.0)
wake_direction = reshape(Vinf ./ magVinf, :, 1)
for i in eachindex(body.Das)
    body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
end

backend = pnl.DirectBackend()
pnl.apply_freestream!(body, Vinf)
pnl.solve!(body, pnl.Backslash(body); backend)
pnl.calcfield_U!(body, Vinf; backend)
pnl.calcfield_P!(body, magVinf, rho)
pnl.calcfield_F!(body)
P_B = copy(body.P)

# Run PressureLaplace (:surface_velocity) and grab L, b
body_l = deepcopy(body)
body_l.needs_velocity_gradient[] = true
body_l.velocity .= 0.0; body_l.velocity_gradient .= 0.0
pnl.calcfield_U!(body_l, Vinf; backend)
pnl.influence!((body_l,), (body_l,), backend;
    scalar_potential=false, velocity=false, velocity_gradient=true)
frames_l = pnl.ReferenceFrame(body_l)
pl = pnl.PressureLaplace((body_l,), rho;
    reference_panel=1, reference_pressure=0.0, verbose=false,
    gradient_mode=:surface_velocity)
pl.velocity_dot[1] .= 0.0
pl((body_l,), (nothing,), frames_l, Vinf, 0, 1.0)

L = pl.L[1]
b_act = pl.b[1]
P_B_shift = P_B .- P_B[1]
LP_B = L * P_B_shift

# Ratio per panel, restricted to entries where LP_B is sufficiently nonzero
# (avoid division by ~0 contributing nonsense to the distribution)
mag_threshold = 0.01 * median(abs.(LP_B))
mask = abs.(LP_B) .> mag_threshold
n_used = sum(mask)
ratios = b_act[mask] ./ LP_B[mask]

println("\n#===== b_actual ./ (L P_B) per-panel =====#")
@printf "panels used (|LP_B| > 1%% of median): %d / %d\n" n_used length(b_act)
@printf "ratio median:   %+10.3f\n" median(ratios)
@printf "ratio mean:     %+10.3f\n" mean(ratios)
@printf "ratio std:      %+10.3f\n" std(ratios)
@printf "ratio p10:      %+10.3f\n" quantile(ratios, 0.10)
@printf "ratio p25:      %+10.3f\n" quantile(ratios, 0.25)
@printf "ratio p75:      %+10.3f\n" quantile(ratios, 0.75)
@printf "ratio p90:      %+10.3f\n" quantile(ratios, 0.90)
@printf "fraction with ratio in [median*0.5, median*1.5]: %.2f\n"  sum(abs.(ratios .- median(ratios)) .< 0.5*abs(median(ratios))) / n_used

# CL of a hypothetical "rescaled" pressure: shift b by 1/median(ratio), re-solve
shift = 1 / median(ratios)
b_rescaled = b_act .* shift
b_rescaled[pl.reference_panel] = pl.reference_pressure
P_rescaled = L \ b_rescaled
body_x = deepcopy(body_l)
body_x.P .= P_rescaled
pnl.calcfield_F!(body_x)
Dhat = Vinf/_norm(Vinf); Shat = [0,1,0]; Lhat = _cross(Dhat, Shat)
nondim = 0.5*rho*magVinf^2*b^2/ar
LDS = pnl.calcfield_LDS(body_x, Lhat, Dhat)
CL_rescaled = sign(_dot(LDS[:,1], Lhat)) * _norm(LDS[:,1]) / nondim

# Original Bernoulli and Laplace for reference
LDS_B = pnl.calcfield_LDS(body, Lhat, Dhat)
CL_B = sign(_dot(LDS_B[:,1], Lhat)) * _norm(LDS_B[:,1]) / nondim
pnl.calcfield_F!(body_l)
LDS_L = pnl.calcfield_LDS(body_l, Lhat, Dhat)
CL_L = sign(_dot(LDS_L[:,1], Lhat)) * _norm(LDS_L[:,1]) / nondim

println("\n#===== CL after uniform rescaling of b =====#")
@printf "CL_bernoulli         = %+8.4f\n" CL_B
@printf "CL_laplace (sv mode) = %+8.4f\n" CL_L
@printf "CL after b *= %+6.4f = %+8.4f\n" shift CL_rescaled
@printf "CL_exp               = %+8.4f\n" 0.238
