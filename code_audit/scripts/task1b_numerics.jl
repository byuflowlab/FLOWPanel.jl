#=##############################################################################
# code_audit Task 1b — Numerical phase: decisive checks from
# docs/wake_solve_schemes.md §5, prioritizing the P5 experiment.
#
# Rebuilds the settled marched state with the exact Task 1 case-b settings
# (pitching_wing.jl wing, n_span=15, n_airfoil=31, n_endcap=5, alpha=4.5 deg,
# dt=0.5 chords/step, 80 steps, panel wake, das=0.05c), then on the frozen
# geometry+state:
#
#  A. P5 experiment: probe the free-wake SCALAR POTENTIAL phi_tr at the body
#     control points, form the TE-pair trace jump C*phi_tr per shedding
#     station, and compare its spanwise distribution/magnitude against the
#     measured Delta-mu deficit (steady minus marched).
#  B. Manual linear systems (§5 steps 1-4, reviewer-corrected):
#     1. G_B (das attached panel) assembled by _G!; verify G_B mu* = -phi_sigma*.
#     2. G_A (semi-infinite attached panel, unit Das) on the same geometry;
#        solve with sigma_A = -Uinf.n; compare Delta-mu against Task 1 (0.912).
#     3. Cross-substitution residual r_A = G_A mu* + phi_sigma_A; spanwise/
#        chordwise structure (TE rows vs distributed).
#     4. Corrected term split at mid-span TE control-point pair, plus:
#        (a) phi_sigma_fw =?= (I - K_int) phi_tr   (eq. 2 / P1-P2 check)
#        (b) off-body composite-sheet check: phi_att,das + phi_fw =?= phi_att,inf
#            with mu* strengths (P3/P4 check).
#     A' hybrid solve: mu_A' = G_B \ (-phi_sigma_A - phi_tr) — scheme B geometry
#     with the free wake entering through its potential instead of sigma.
#
# Outputs: code_audit/results/task1b/{task1b_spanwise.csv, task1b_residual_cells.csv,
#          task1b_summary.txt, task1b_spanwise.png}
# Run:  julia --project -t auto code_audit/scripts/task1b_numerics.jl
#       TASK1B_SMOKE=1 for a coarse/fast smoke run.
=###############################################################################

import FLOWPanel as pnl
import FLOWPanel.FastMultipole as FastMultipole
using FLOWPanel.FastMultipole.StaticArrays
import LinearAlgebra: dot, norm, lu
import Printf: @printf, @sprintf
import Statistics: mean, median, cor

include(joinpath(@__DIR__, "..", "..", "examples", "pitching_wing.jl"))

const RESULTS = joinpath(@__DIR__, "..", "results", "task1b")
mkpath(RESULTS)
const SMOKE = get(ENV, "TASK1B_SMOKE", "0") == "1"

summary_io = IOBuffer()
macro logf(fmt, args...)
    esc(quote
        @printf($fmt, $(args...))
        @printf(summary_io, $fmt, $(args...))
    end)
end

# ------------------------------------------------------------------ parameters
alpha_deg    = 4.5
c_ft         = 1.0
aspect_ratio = 4.0
dims  = _resolve_pitching_wing_dimensions(c_ft, aspect_ratio, nothing, nothing)
c, b  = dims.c, dims.b
U     = 330.2 * FT_TO_M
rho   = 1.225
Sref  = b * c
qinf  = 0.5 * rho * U^2

n_span, n_airfoil, n_endcap = SMOKE ? (5, 21, 5) : (15, 31, 5)
c_per_dt = 0.5
nsteps   = SMOKE ? 12 : 80
dt       = c / U * c_per_dt
t_range  = range(0.0, step=dt, length=nsteps + 1)
das      = 0.05 * c
pivot    = SVector{3}(0.25 * c, 0.0, 0.0)

backend_fmm = pnl.FastMultipoleBackend(expansion_order=8,
    multipole_acceptance=0.4, leaf_size=40)
backend_dir = pnl.DirectBackend()

@logf("Task 1b numerics — %s run\n", SMOKE ? "SMOKE" : "FULL")
@logf("  wing: c=%.4g m, b=%.4g m, n_span=%d n_airfoil=%d n_endcap=%d\n",
    c, b, n_span, n_airfoil, n_endcap)
@logf("  alpha=%.2f deg, U=%.4g m/s, dt=%.2g chords/step, nsteps=%d, das=%.3g m (0.05c)\n",
    alpha_deg, U, c_per_dt, nsteps, das)

# ------------------------------------------- rebuild the settled marched state
println("\n[1] Rebuilding settled marched state (Task 1 case b)...")
wing = build_pitching_wing_body(c, b;
    n_span, n_airfoil, n_endcap, semiinfinite_wake=false)
frames = pitching_wing_frame(wing, pivot, deg2rad(alpha_deg))
Uinf(t) = SVector{3}(U, 0.0, 0.0)
uinf = Uinf(0.0)
set_wake_Das!(wing, uinf; magnitude=das)

wake = pnl.PanelWake(wing; nwakerows=nsteps + 2, include_final_filament=false)
steady_maneuver!(frames, systems, wakes, t) = nothing

pressure = pnl.PressureBernoulli(rho; correct_kuttacondition=false)
force = pnl.ForceMonitor(length(t_range), 1;
    i_frame=-1, normalization=pnl.WingNormalization(rho, Sref, c),
    correct_kuttacondition=false, verbose=false)

tmarch = @elapsed pnl.simulate!((wing,), (wake,), frames, steady_maneuver!, Uinf, t_range;
    body_solvers=(pnl.Backslash(wing),),
    backend=backend_fmm,
    monitors=(pressure, force),
    path=nothing, name="task1b_panel",
    set_Das_eta_freestream=NaN,
    verbose=false)

CL_settled = force.force[3, end]
@logf("  marched %d steps in %.1f s; settled CL = %.8g (Task 1: 0.262051)\n",
    nsteps, tmarch, CL_settled)

n = wing.ncells
mu_star    = copy(wing.strength[:, 2])
sigma_star = copy(wing.strength[:, 1])

# freeze kerneloffset regime at the solve value (production uses kerneloffset_panel)
wing.kerneloffset = wing.kerneloffset_panel
@assert wing.kerneloffset_panel == wing.kerneloffset_targets  # single-knob body here

# TE pair bookkeeping (upper = shedding row 1, lower = row 4)
shed = wing.shedding[1]
nst = size(shed, 2)
upper = shed[1, :]
lower = shed[4, :]
@assert all(lower .> 0) "expected closed TE with lower panels everywhere"
ysta  = [wing.controlpoints[2, upper[i]] for i in 1:nst]
eta   = ysta ./ (b / 2)
dmu(mu) = [mu[upper[i]] - mu[lower[i]] for i in 1:nst]

# ------------------------------------------------------------------- utilities
"Interior source-layer potential at the body CPs — exact production RHS path."
function body_source_potential(sigma, backend)
    old_strength = copy(wing.strength)
    old_potential = copy(wing.potential)
    wing.strength[:, 1] .= sigma
    wing.strength[:, 2] .= 0.0
    wing.potential .= 0.0
    pnl.influence!(wing, wing, backend; scalar_potential=true, velocity=false)
    out = copy(wing.potential)
    wing.strength .= old_strength
    wing.potential .= old_potential
    return out
end

"Scalar potential of `sources` (tuple) at arbitrary positions."
function probe_potential(positions, sources, backend)
    probes = FastMultipole.ProbeSystem(length(positions), Float64)
    for (i, x) in enumerate(positions)
        probes.position[i] = SVector{3}(x[1], x[2], x[3])
        probes.scalar_potential[i] = 0.0
    end
    pnl.influence!((probes,), sources, backend;
        precalc=false, scalar_potential=true, velocity=false)
    return [probes.scalar_potential[i] for i in eachindex(positions)]
end

"Direct per-panel body influence (incl. attached wake) at an off-body point."
function body_phi_direct(x::SVector{3,Float64})
    ds = FastMultipole.DerivativesSwitch(true, false, false)
    phi = 0.0
    for i in 1:n
        p, _, _ = pnl.induced(x, wing, i, ds; kerneloffset=wing.kerneloffset_panel)
        phi += p
    end
    return phi
end

set_unit_Das!() = (for D in wing.Das, j in axes(D, 2); D[:, j] ./= norm(view(D, :, j)); end)
saved_Das = [copy(D) for D in wing.Das]
restore_Das!() = (for (D, S) in zip(wing.Das, saved_Das); D .= S; end)

# ------------------------------------------------- [2] free-wake potential trace
println("\n[2] Probing free-wake scalar potential at body control points...")
wake_srcs = pnl.get_sources(wake)
cps = [SVector{3}(wing.controlpoints[:, p]) for p in 1:n]
tprobe = @elapsed phi_tr = probe_potential(cps, wake_srcs, backend_dir)
phi_tr_fmm = probe_potential(cps, wake_srcs, backend_fmm)
@logf("  phi_tr probed (direct, %.1f s); |phi_tr| median=%.4g max=%.4g\n",
    tprobe, median(abs.(phi_tr)), maximum(abs.(phi_tr)))
@logf("  direct-vs-FMM probe agreement: max|diff| = %.3e (rel to max|phi_tr| = %.3e)\n",
    maximum(abs.(phi_tr .- phi_tr_fmm)), maximum(abs.(phi_tr .- phi_tr_fmm)) / maximum(abs.(phi_tr)))

Cphi = [phi_tr[upper[i]] - phi_tr[lower[i]] for i in 1:nst]   # C * phi_tr per station

# ----------------------------------------------------- [3] assemble the systems
println("\n[3] Assembling influence matrices (direct _G!)...")
GB = zeros(n, n)
t1 = @elapsed pnl._G!(GB, wing, wing; kerneloffset=wing.kerneloffset_panel)

# K_int: same nodes/cells, no shedding (ensure_winding=false — cells already wound)
clone = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}, 2, Float64, true}(
    copy(wing.nodes), copy(wing.cells), zeros(Int, 6, 0);
    kerneloffset=wing.kerneloffset_panel,
    kernelcutoff=wing.kernelcutoff,
    semiinfinite_wake=false, watertight=true,
    ensure_winding=false, check_mesh=false)
pnl.calc_normals!(clone); pnl.calc_controlpoints!(clone)
@assert maximum(abs.(clone.controlpoints .- wing.controlpoints)) < 1e-12
Kint = zeros(n, n)
t2 = @elapsed pnl._G!(Kint, wing, clone; kerneloffset=wing.kerneloffset_panel)

# G_A: semi-infinite attached wake on the frozen geometry, unit Das
wing.semiinfinite_wake = true
set_unit_Das!()
GA = zeros(n, n)
t3 = @elapsed pnl._G!(GA, wing, wing; kerneloffset=wing.kerneloffset_panel)
wing.semiinfinite_wake = false
restore_Das!()
@logf("  G_B %.1f s, K_int %.1f s, G_A %.1f s (n=%d)\n", t1, t2, t3, n)

KwC_B = GB .- Kint    # das attached-panel coupling columns
KwC_A = GA .- Kint    # semiinfinite attached-panel coupling columns

# RHS potentials
sigma_A = [-dot(uinf, SVector{3}(wing.normals[:, p])) for p in 1:n]
phi_sStar = body_source_potential(sigma_star, backend_fmm)   # production-path RHS of B
phi_sA    = body_source_potential(sigma_A, backend_fmm)
phi_sfw   = phi_sStar .- phi_sA                              # source potential of the u_fw part

# ------------------------------------------------------------ [4] solves/residuals
println("\n[4] Solving and cross-substituting...")
luGB = lu(GB); luGA = lu(GA)

# sanity: marched solution satisfies scheme B's system
rB = GB * mu_star .+ phi_sStar
@logf("\n(§5.1) scheme B sanity: ||G_B mu* + phi_sigma*|| / ||phi_sigma*|| = %.3e (max row %.3e)\n",
    norm(rB) / norm(phi_sStar), maximum(abs.(rB)) / maximum(abs.(phi_sStar)))

# scheme A solve on frozen geometry
mu_A = luGA \ (-phi_sA)
# scheme A' (hybrid): das geometry, free wake through POTENTIAL (not sigma)
mu_Ap = luGB \ (-(phi_sA .+ phi_tr))

dmu_star, dmu_A, dmu_Ap = dmu(mu_star), dmu(mu_A), dmu(mu_Ap)
@logf("(§5.2) scheme A on frozen geometry: sum(dmu*)/sum(dmu_A) = %.4f (Task 1: 0.916)\n",
    sum(dmu_star) / sum(dmu_A))
imid = argmin(abs.(ysta))
@logf("        mid-span dmu*/dmu_A = %.4f (Task 1: 0.912)\n", dmu_star[imid] / dmu_A[imid])

# cross-substitution residual of marched state in scheme A's system
rA = GA * mu_star .+ phi_sA
te_rows = sort(unique(vcat(upper, lower)))
other_rows = setdiff(1:n, te_rows)
rms(v) = sqrt(mean(abs2, v))
@logf("(§5.3) cross-substitution r_A = G_A mu* + phi_sigma_A:\n")
@logf("        rms(all)=%.4g rms(TE rows)=%.4g rms(non-TE)=%.4g max|r_A|=%.4g (at cell %d)\n",
    rms(rA), rms(rA[te_rows]), rms(rA[other_rows]), maximum(abs.(rA)), argmax(abs.(rA)))
@logf("        scale: max|phi_sigma_A| = %.4g\n", maximum(abs.(phi_sA)))

# ------------------------------------------------------ [5] P5 deficit comparison
println("\n[5] P5 comparison: deficit vs C*phi_tr ...")
D1 = dmu_A .- dmu_star          # measured deficit per station (steady - marched)
D2 = dmu_A .- dmu_Ap            # scheme-difference part (P3/P4; ~0 if composite sheet OK)
D3 = dmu_Ap .- dmu_star         # part explained by wake-potential-through-sigma (P5 route)
D4 = -Cphi                      # leading P5 prediction for D3
GinvKwC_phi = luGB \ (KwC_B * phi_tr)
D5 = -(Cphi .- [GinvKwC_phi[upper[i]] - GinvKwC_phi[lower[i]] for i in 1:nst])  # full eq-(3) P5 prediction

@logf("\n  station-mean values: D1(measured)=%.5g D2(A-A')=%.5g D3(A'-marched)=%.5g\n",
    mean(D1), mean(D2), mean(D3))
@logf("  P5 predictions:      D4(-C.phi_tr)=%.5g D5(eq.3 full)=%.5g\n", mean(D4), mean(D5))
@logf("  correlations: cor(D1,D4)=%.4f cor(D1,D3)=%.4f cor(D3,D5)=%.4f cor(D1,D2)=%.4f\n",
    cor(D1, D4), cor(D1, D3), cor(D3, D5), cor(D1, D2))
@logf("  magnitude ratios: mean(D3)/mean(D1)=%.4f mean(D4)/mean(D1)=%.4f mean(D5)/mean(D3)=%.4f mean(D2)/mean(D1)=%.4f\n",
    mean(D3) / mean(D1), mean(D4) / mean(D1), mean(D5) / mean(D3), mean(D2) / mean(D1))

# geometric P5 scale estimate at mid-span
xu = SVector{3}(wing.controlpoints[:, upper[imid]])
xl = SVector{3}(wing.controlpoints[:, lower[imid]])
@logf("  mid-span TE CP pair separation |xu-xl| = %.4g m (das = %.4g m)\n", norm(xu - xl), das)
@logf("  mid-span: D1=%.5g D2=%.5g D3=%.5g D4=%.5g D5=%.5g dmu*=%.5g\n",
    D1[imid], D2[imid], D3[imid], D4[imid], D5[imid], dmu_star[imid])

# --------------------------------------------- [6] term split + eq-(2) checks
println("\n[6] Term split and corrected §5.4 checks...")
att_A = KwC_A * mu_star    # semiinfinite attached-wake potential at CPs, strengths mu*
att_B = KwC_B * mu_star    # das attached-wake potential at CPs
for (lbl, row) in (("upper", upper[imid]), ("lower", lower[imid]))
    @logf("  mid-span TE %s CP (cell %d):\n", lbl, row)
    @logf("    (i)   semiinf attached phi   = %+.6g\n", att_A[row])
    @logf("    (ii)  das attached phi       = %+.6g\n", att_B[row])
    @logf("    (iii) free-wake phi_tr       = %+.6g\n", phi_tr[row])
    @logf("    (iv)  phi_sigma_fw           = %+.6g\n", phi_sfw[row])
    @logf("    (ii)+(iii)-(i) [naive P3 gauge-confounded on-body] = %+.6g\n",
        att_B[row] + phi_tr[row] - att_A[row])
end

# (a) eq-(2): phi_sigma_fw =?= (I - K_int) phi_tr over all CPs
e4a = phi_sfw .- (phi_tr .- Kint * phi_tr)
@logf("\n(§5.4a) ||phi_sigma_fw - (I-K_int)phi_tr|| / ||phi_sigma_fw|| = %.4g\n",
    norm(e4a) / norm(phi_sfw))
@logf("        max|err| = %.4g, max|phi_sigma_fw| = %.4g, cor = %.5f\n",
    maximum(abs.(e4a)), maximum(abs.(phi_sfw)),
    cor(phi_sfw, phi_tr .- Kint * phi_tr))
@logf("        TE-row rms(err)=%.4g vs non-TE rms(err)=%.4g\n",
    rms(e4a[te_rows]), rms(e4a[other_rows]))

# (b) off-body composite-sheet check with mu* strengths
println("\n(§5.4b) off-body composite-sheet check (P3/P4)...")
x_te_mid = 0.5 .* (xu .+ xl)
offpts = [x_te_mid + SVector{3}(ddx * c, 0.0, ddz * c)
          for ddx in (0.5, 1.0, 2.0) for ddz in (0.5, -0.5)]

old_strength = copy(wing.strength)
wing.strength[:, 1] .= 0.0
wing.strength[:, 2] .= mu_star
phiB_att = [body_phi_direct(x) for x in offpts]     # finite attached wake (das)
wing.semiinfinite_wake = true
set_unit_Das!()
phiA_att = [body_phi_direct(x) for x in offpts]     # semiinfinite attached wake
wing.semiinfinite_wake = false
restore_Das!()
wing.strength .= old_strength

phi_fw_off = probe_potential(offpts, wake_srcs, backend_dir)
@logf("  point (dx/c,dz/c)   phi_A,att     phi_B,att+phi_fw   rel.err\n")
for (k, x) in enumerate(offpts)
    comp = phiB_att[k] + phi_fw_off[k]
    @logf("  (%+.1f,%+.1f)      %+.6e   %+.6e   %.3e\n",
        (x[1] - x_te_mid[1]) / c, (x[3] - x_te_mid[3]) / c,
        phiA_att[k], comp, abs(comp - phiA_att[k]) / max(abs(phiA_att[k]), 1e-300))
end

# ------------------------------------------ [7] exact deficit decomposition
# Exact identity at the settled fixed point:
#   r_A = G_A mu* + phi_sigma_A = T1 + T2,   D1 = -C G_A^-1 (T1 + T2)
#   T1 = phi_rem(mu*) - phi_tr : rigid remainder sheet (TE+das -> inf, strengths
#        C mu*) minus actual free-wake potential  -> composite-sheet consistency
#        at MATCHED current strengths (geometry relaxation + truncation + Q3
#        regularization + any handoff bug)
#   T2 = phi_tr - phi_sigma_fw : wake-potential-through-sources error
#        = K_int phi_tr (P5 route) + e4a (P1/P2 sampling error)
println("\n[7] Exact deficit decomposition D1 = E1(T1) + E2(T2)...")
phi_rem = (KwC_A .- KwC_B) * mu_star
T1 = phi_rem .- phi_tr
T2 = phi_tr .- phi_sfw
@assert maximum(abs.((T1 .+ T2) .- rA)) < 1e-10 * max(1.0, maximum(abs.(rA)))

Cimg(v) = [v[upper[i]] - v[lower[i]] for i in 1:nst]   # C * vector
E1  = -Cimg(luGA \ T1)
E2  = -Cimg(luGA \ T2)
# T2 = phi_tr - phi_sfw = Kint*phi_tr - e4a  (by definition of e4a)
E2a = -Cimg(luGA \ (Kint * phi_tr))   # pure P5 (trace shift through Kutta coupling)
E2b = Cimg(luGA \ e4a)                # P1/P2 sigma-sampling error; E2 = E2a + E2b
@logf("\n  exact split of measured deficit (station means):\n")
@logf("    D1 = %.5g;  E1(T1, composite-sheet) = %.5g (%.1f%%);  E2(T2, potential-vs-sigma) = %.5g (%.1f%%)\n",
    mean(D1), mean(E1), 100 * mean(E1) / mean(D1), mean(E2), 100 * mean(E2) / mean(D1))
@logf("    E2 split: E2a(P5 trace/Kutta) = %.5g (%.1f%% of D1); E2b(P1/P2 sampling) = %.5g (%.1f%% of D1)\n",
    mean(E2a), 100 * mean(E2a) / mean(D1), mean(E2b), 100 * mean(E2b) / mean(D1))
@logf("    closure: max|E1+E2-D1| = %.3e, max|E2a+E2b-E2| = %.3e\n",
    maximum(abs.(E1 .+ E2 .- D1)), maximum(abs.(E2a .+ E2b .- E2)))
@logf("    A' vs A: sum(dmu_Ap)/sum(dmu_A) = %.4f (marched/A = %.4f)\n",
    sum(dmu_Ap) / sum(dmu_A), sum(dmu_star) / sum(dmu_A))

# T1 split: geometry relaxation vs truncation/discretization, via a synthetic
# FLAT wake (same strengths/connectivity, nodes on the rigid sheet at U*dt pitch)
swake = deepcopy(wake)
xhat = SVector{3}(1.0, 0.0, 0.0)
for i_surf in eachindex(swake.nodes)
    nodes = swake.nodes[i_surf]
    for s in axes(nodes, 3)
        base = SVector{3}(nodes[1, 1, s], nodes[2, 1, s], nodes[3, 1, s])
        for k in 2:swake.nwakes[]+1
            nodes[:, k, s] .= base .+ (k - 1) * U * dt .* xhat
        end
    end
end
phi_flat = probe_potential(cps, pnl.get_sources(swake), backend_dir)
T1_geom  = phi_flat .- phi_tr          # relaxed-vs-rigid wake geometry (physics)
T1_rest  = phi_rem .- phi_flat         # truncation + row discretization + core_size
E1_geom  = -Cimg(luGA \ T1_geom)
E1_rest  = -Cimg(luGA \ T1_rest)
@logf("  T1 split (station means of C-image through G_A^-1):\n")
@logf("    E1_geom(relaxed wake geometry) = %.5g (%.1f%% of D1)\n",
    mean(E1_geom), 100 * mean(E1_geom) / mean(D1))
@logf("    E1_rest(truncation+row-discretization+regularization) = %.5g (%.1f%% of D1)\n",
    mean(E1_rest), 100 * mean(E1_rest) / mean(D1))

# wake-sheet displacement off the rigid plane (context for T1_geom)
dz_all = Float64[]
for i_surf in eachindex(wake.nodes)
    nodes = wake.nodes[i_surf]
    for s in axes(nodes, 3), k in 2:wake.nwakes[]+1
        push!(dz_all, nodes[3, k, s] - nodes[3, 1, s])
    end
end
@logf("  relaxed wake z-displacement off rigid sheet: median=%.4g m, max|.|=%.4g m (c=%.4g m)\n",
    median(dz_all), maximum(abs.(dz_all)), c)

# ------------------------------------------------------------------------- CSVs
open(joinpath(RESULTS, "task1b_spanwise.csv"), "w") do io
    println(io, "eta,y,dmu_marched,dmu_A,dmu_Aprime,Cphi_tr,D1_measured_deficit,D2_A_minus_Aprime,D3_Aprime_minus_marched,D4_negCphi,D5_eq3,E1_composite,E2_pot_vs_sigma,E2a_P5,E2b_P1P2,E1_geom,E1_rest")
    for i in 1:nst
        println(io, join((eta[i], ysta[i], dmu_star[i], dmu_A[i], dmu_Ap[i],
            Cphi[i], D1[i], D2[i], D3[i], D4[i], D5[i],
            E1[i], E2[i], E2a[i], E2b[i], E1_geom[i], E1_rest[i]), ","))
    end
end
open(joinpath(RESULTS, "task1b_residual_cells.csv"), "w") do io
    println(io, "cell,x,y,z,rA,rB,phi_tr,phi_sigma_fw,err4a,is_te_upper,is_te_lower")
    us, ls = Set(upper), Set(lower)
    for p in 1:n
        println(io, join((p, wing.controlpoints[1, p], wing.controlpoints[2, p],
            wing.controlpoints[3, p], rA[p], rB[p], phi_tr[p], phi_sfw[p], e4a[p],
            Int(p in us), Int(p in ls)), ","))
    end
end

# ------------------------------------------------------------------------- plot
try
    import PythonPlot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
    axs[0].plot(eta, dmu_A, "-k", label="Δμ steady (A, frozen geom)")
    axs[0].plot(eta, dmu_star, "--", color="C0", label="Δμ marched (B)")
    axs[0].plot(eta, dmu_Ap, ":", color="C2", label="Δμ A′ (wake via potential)")
    axs[0].set_xlabel("η"); axs[0].set_ylabel("Δμ"); axs[0].legend(fontsize=8)
    axs[0].grid(true, alpha=0.35)
    axs[1].plot(eta, D1, "-k", label="D1 measured deficit (A−B)")
    axs[1].plot(eta, D3, "--", color="C2", label="D3 = A′−B")
    axs[1].plot(eta, D4, ":", color="C3", label="D4 = −C·φ_tr")
    axs[1].plot(eta, D5, "-.", color="C4", label="D5 = eq.(3) P5")
    axs[1].plot(eta, D2, "-", color="C7", alpha=0.6, label="D2 = A−A′ (P3/P4)")
    axs[1].plot(eta, E1, "-", color="C1", alpha=0.8, label="E1 (composite sheet, exact)")
    axs[1].plot(eta, E2, "--", color="C5", alpha=0.8, label="E2 (pot-vs-σ, exact)")
    axs[1].plot(eta, E1_geom, ":", color="C8", label="E1_geom (wake relaxation)")
    axs[1].set_xlabel("η"); axs[1].set_ylabel("Δμ deficit"); axs[1].legend(fontsize=8)
    axs[1].grid(true, alpha=0.35)
    fig.suptitle("Task 1b: Δμ deficit vs free-wake trace jump (P5)")
    fig.tight_layout()
    fig.savefig(joinpath(RESULTS, "task1b_spanwise.png"), dpi=170)
    plt.pyplot.close(fig)
catch err
    @warn "plotting failed" err
end

open(joinpath(RESULTS, "task1b_summary.txt"), "w") do io
    write(io, String(take!(summary_io)))
end
println("\nWrote outputs under $(RESULTS)")
