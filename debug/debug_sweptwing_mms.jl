## Manufactured-solution test of the L vs RHS internal consistency.
##
## Pick a smooth analytic `P_target(x)` on the wing. By Bernoulli:
##    a_target = -∇P_target / ρ
## Then run the SAME edge-loop machinery as `_pressure_rhs!` with this
## constructed acceleration field to produce `b_synth`. Independently
## compute `L * P_target` from the assembled operator.
##
## Internal consistency requires `b_synth ≈ L P_target` per panel. If
## they differ by a uniform scalar, that's the operator/RHS convention
## factor. If they differ pointwise, the operator and RHS edge loops
## are evaluating physically different things.
##
## Test cases:
##   P1(x) = x        (linear in x)         → ∇P = (1,0,0),    a = (-1/ρ,0,0)
##   P2(x) = x^2      (quadratic in x)      → ∇P = (2x,0,0),   ∇²P = 2
##   P3(x) = x*z      (bilinear)            → ∇P = (z,0,x),    ∇²P = 0

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))
using LinearAlgebra: norm, dot
using Statistics: quantile, median, mean
using Printf

# --- wing setup (same as sweptwing.jl) ---
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

body = simplewing(b, ar, 1.0, 0, 0, 45, 0;
                  bodytype=bodytype, bodyoptargs=(; CPoffset=1e-14),
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

# Stand up a PressureLaplace instance just to get L, edges, & ref panel.
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
edges = pl.edges[1]
ref_p = pl.reference_panel
ref_P = pl.reference_pressure
n = body_l.ncells

# Replicate the per-edge co-normal weight calculation (matches
# `_pressure_edge_conormal_weight` in src).
function edge_metric(body, edge_a, edge_b, i, j)
    r = body.controlpoints[:, j] .- body.controlpoints[:, i]
    rmag2 = sum(r .* r)
    e = body.nodes[:, edge_b] .- body.nodes[:, edge_a]
    ell = norm(e)
    t = e ./ ell
    n_avg = body.normals[:, i] .+ body.normals[:, j]
    n_avg ./= norm(n_avg)
    nu = [t[2]*n_avg[3] - t[3]*n_avg[2],
          t[3]*n_avg[1] - t[1]*n_avg[3],
          t[1]*n_avg[2] - t[2]*n_avg[1]]
    nu ./= norm(nu)
    if dot(nu, r) < 0
        nu .*= -1
    end
    w = ell * dot(nu, r) / rmag2
    return w, ell, nu, n_avg
end

# Reproduce `_pressure_rhs!` edge loop with arbitrary acceleration.
function build_b(body, edges, acceleration, rho, ref_p, ref_P)
    b = zeros(size(acceleration, 2))
    for k in axes(edges, 2)
        ea, eb, i, j = edges[1,k], edges[2,k], edges[3,k], edges[4,k]
        _, ell, nu, n_avg = edge_metric(body, ea, eb, i, j)
        ai = acceleration[:, i]
        aj = acceleration[:, j]
        # tangent-project both sides using the averaged edge normal
        ai_t = ai .- dot(ai, n_avg) .* n_avg
        aj_t = aj .- dot(aj, n_avg) .* n_avg
        flux = ell * 0.5 * (dot(ai_t, nu) + dot(aj_t, nu))
        b[i] += rho * flux
        b[j] -= rho * flux
    end
    b[ref_p] = ref_P
    return b
end

function test_field(label, P_target_fn, gradP_fn)
    P_t = [P_target_fn(view(body_l.controlpoints, :, p)) for p in 1:n]
    accel = zeros(3, n)
    for p in 1:n
        accel[:, p] .= -gradP_fn(view(body_l.controlpoints, :, p)) ./ rho
    end
    # Gauge-shift P_target so reference panel sits at ref_P (operator pins it)
    P_t_shifted = P_t .- (P_t[ref_p] - ref_P)
    LP = L * P_t_shifted
    b_s = build_b(body_l, edges, accel, rho, ref_p, ref_P)

    diff = LP .- b_s
    # exclude reference row (artificially set to ref_P / 1·P[ref])
    mask = trues(n); mask[ref_p] = false
    rel_err = norm(diff[mask]) / max(norm(b_s[mask]), eps())
    cos_align = dot(LP[mask], b_s[mask]) / (norm(LP[mask]) * norm(b_s[mask]))

    println("\n#===== Manufactured solution: $label =====#")
    @printf "‖LP_target‖   = %.3e   ‖b_synth‖ = %.3e\n" norm(LP[mask]) norm(b_s[mask])
    @printf "‖LP - b‖/‖b‖ = %.3e\n" rel_err
    @printf "cosine(LP, b) = %+.6f\n" cos_align
    # element-wise ratio for the magnitude diagnostic
    big = abs.(b_s[mask]) .> 0.01 * median(abs.(b_s[mask]))
    if any(big)
        ratios = LP[mask][big] ./ b_s[mask][big]
        @printf "ratio LP/b median=%+.4f mean=%+.4f p25=%+.4f p75=%+.4f\n" median(ratios) mean(ratios) quantile(ratios, 0.25) quantile(ratios, 0.75)
    end
end

test_field("P = x",          c -> c[1],         c -> [1.0, 0.0, 0.0])
test_field("P = x^2",        c -> c[1]^2,       c -> [2*c[1], 0.0, 0.0])
test_field("P = x*z",        c -> c[1]*c[3],    c -> [c[3], 0.0, c[1]])
test_field("P = x^2 + y^2",  c -> c[1]^2 + c[2]^2, c -> [2*c[1], 2*c[2], 0.0])
