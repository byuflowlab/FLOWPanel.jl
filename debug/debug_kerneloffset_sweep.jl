## Kerneloffset sweep for capped Dirichlet wing.
## For each kerneloffset, build/solve the body, evaluate raw induced velocity +
## Hessian, run PressureBernoulli + PressureLaplace, and report:
##   - median/p99 of |u·n| and |u·n|/|u| over the surface
##   - same statistics restricted to high-AR panels (top decile by aspect ratio)
##   - relative residual ||L*(P_B - P_B[ref]) - b_L|| / ||b_L||
##   - integrated Bernoulli and Laplace forces (Z component, body frame)
##
## Hypothesis: if shrinking kerneloffset → 0 reduces |u·n| on sliver panels and
## brings Bernoulli/Laplace forces into agreement, regularization mismatch
## between the BIE matrix (using same kernel) and the velocity evaluation
## (regularized) is the leading driver of the localized error.

import FLOWPanel as pnl
using LinearAlgebra: norm, dot
using Statistics: median, quantile
using Printf

include(joinpath(pnl.examples_path, "simple_wing_capped_pressure_comparison.jl"))

function panel_aspect_ratios(body)
    AR = zeros(body.ncells)
    for p in 1:body.ncells
        ns = body.cells[:, p]
        v1 = body.nodes[:, ns[1]]
        v2 = body.nodes[:, ns[2]]
        v3 = body.nodes[:, ns[3]]
        e1 = norm(v2 .- v1)
        e2 = norm(v3 .- v2)
        e3 = norm(v1 .- v3)
        emin = min(e1, e2, e3)
        emax = max(e1, e2, e3)
        AR[p] = emin > 0 ? emax / emin : Inf
    end
    return AR
end

cross3(a, b) = [a[2]*b[3] - a[3]*b[2],
                a[3]*b[1] - a[1]*b[3],
                a[1]*b[2] - a[2]*b[1]]

function integrate_force(b)
    F = zeros(3)
    for p in 1:b.ncells
        ns = b.cells[:, p]
        v1 = b.nodes[:, ns[1]]; v2 = b.nodes[:, ns[2]]; v3 = b.nodes[:, ns[3]]
        area = 0.5 * norm(cross3(v2 .- v1, v3 .- v1))
        F .+= -b.P[p] * area .* b.normals[:, p]
    end
    return F
end

function tangency_stats(body)
    udotn = zeros(body.ncells)
    umag = zeros(body.ncells)
    for p in 1:body.ncells
        n = @view body.normals[:, p]
        u = @view body.velocity[:, p]
        udotn[p] = dot(u, n)
        umag[p] = norm(u)
    end
    rel = abs.(udotn) ./ max.(umag, eps())
    return abs.(udotn), rel, umag
end

function run_sweep(; meshfile=joinpath(pnl.examples_path, "data", "wing_ar4_naca0016_5.msh"),
                    AOA=5.88, magVinf=56.0, rho=1.225,
                    kernel_offsets=[1e-2, 3e-3, 1e-3, 3e-4, 1e-4],
                    backend=pnl.DirectBackend(),
                    dt=1.0)

    println()
    println("#===== KERNELOFFSET SWEEP =====#")
    @printf("mesh: %s\n", basename(meshfile))
    @printf("AOA=%.3f deg, |Vinf|=%.2f m/s, rho=%.3f\n", AOA, magVinf, rho)

    header = @sprintf("%-9s %-8s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s",
        "k_off", "npanels", "med|u.n|", "p99|u.n|", "med|u.n|/|u|", "p99|u.n|/|u|",
        "hiAR med rel", "rel_resid", "F_B_z", "F_L_z")
    println(header)
    println("-"^length(header))

    for k_off in kernel_offsets
        body = build_pressure_comparison_wing(; kernel_offset=k_off)
        Vinf = magVinf .* [cosd(AOA), 0.0, sind(AOA)]
        set_pressure_comparison_wake!(body, Vinf)
        pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
            body_solvers=pnl.Backslash(body), backend, verbose=false)

        body.needs_velocity_gradient[] = true
        body.velocity .= 0.0
        body.velocity_gradient .= 0.0
        pnl.calcfield_U!(body, Vinf; backend)
        pnl.influence!((body,), (body,), backend;
            scalar_potential=false, velocity=false, velocity_gradient=true)

        AR = panel_aspect_ratios(body)
        ar_thresh = quantile(AR, 0.9)
        hi_AR = AR .>= ar_thresh

        absudotn, relUN, umag = tangency_stats(body)
        med_u = median(absudotn)
        p99_u = quantile(absudotn, 0.99)
        med_r = median(relUN)
        p99_r = quantile(relUN, 0.99)
        med_r_hiAR = median(relUN[hi_AR])

        # Bernoulli pressure on a fresh copy
        body_b = deepcopy(body)
        frames_b = pnl.ReferenceFrame(body_b)
        pb_mon = pnl.PressureBernoulli(rho; unsteady=false,
            correct_kuttacondition=false, backend=backend)
        pb_mon((body_b,), (nothing,), frames_b, Vinf, 0, dt)

        # Laplace pressure on a fresh copy
        body_l = deepcopy(body)
        frames_l = pnl.ReferenceFrame(body_l)
        pl_mon = pnl.PressureLaplace((body_l,), rho;
            reference_panel=1, reference_pressure=0.0, verbose=false)
        # Preload velocity_dot to suppress first-call ∂u/∂t transient.
        # Use the same tangent + kinematic basis the monitor uses internally.
        for p in 1:body_l.ncells
            n = @view body_l.normals[:, p]
            u = @view body_l.velocity[:, p]
            tang = u .- dot(u, n) .* n
            pl_mon.velocity_dot[1][:, p] .= -(tang .+ body_l.velocity_kinematic[:, p])
        end
        pl_mon((body_l,), (nothing,), frames_l, Vinf, 0, dt)

        # Manufactured-pressure residual (Bernoulli pressure threaded through L)
        L = pl_mon.L[1]
        b_L = pl_mon.b[1]
        P_B_ref = body_b.P .- body_b.P[1]
        LP = L * P_B_ref
        resid_rel = norm(LP .- b_L) / max(norm(b_L), eps())

        F_B = integrate_force(body_b)
        F_L = integrate_force(body_l)

        @printf("%-9.2e %-8d %-12.3e %-12.3e %-12.3e %-12.3e %-12.3e %-12.3e %-12.3e %-12.3e\n",
            k_off, body.ncells, med_u, p99_u, med_r, p99_r, med_r_hiAR,
            resid_rel, F_B[3], F_L[3])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_sweep()
end
