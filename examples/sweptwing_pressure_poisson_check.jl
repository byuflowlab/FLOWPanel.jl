#=##############################################################################
# Diagnostic: swept-wing lift and surface-velocity reconstruction checks.
#
# This is intentionally example-local. It builds Weber/Brebner swept wings,
# runs the current monitor-owned pressure/force path, and prints compact tables
# for the observed CL underprediction without writing VTK, plots, or touching src/.
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))

import LinearAlgebra
using Printf
using Statistics: median, quantile

const LA = LinearAlgebra
const CLexp = 0.238

function sweptwing_constants()
    AOA = 4.2
    magVinf = 30.0
    Vinf = magVinf * [cosd(AOA), 0.0, sind(AOA)]
    rho = 1.225
    b = 98 * 0.0254
    ar = 5.0
    return (; AOA, magVinf, Vinf, rho, b, ar,
        tr=1.0, twist_root=0.0, twist_tip=0.0, lambda=45.0, gamma=0.0,
        airfoil="airfoil-rae101.csv", airfoil_path=joinpath(pnl.examples_path, "data"),
        Sref=b^2 / ar, c_ref=b / ar)
end

function discretization(n_rfl::Int, n_span::Int)
    rfl = [(0.25, n_rfl, 10.0, false),
           (0.50, n_rfl, 1.0, true),
           (0.25, n_rfl, 0.1, false)]
    span = [(1.0, n_span, 1.0, true)]
    return rfl, span
end

function build_sweptwing(; builder::Symbol=:mirrored, n_rfl::Int=4, n_span::Int=16,
                         bodyoptargs=(;))
    c = sweptwing_constants()
    rfl, span = discretization(n_rfl, n_span)
    bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
    kwargs = (;
        bodytype,
        bodyoptargs=bodyoptargs,
        airfoil_root=c.airfoil,
        airfoil_tip=c.airfoil,
        airfoil_path=c.airfoil_path,
        rfl_NDIVS=rfl,
        span_NDIVS=span,
        delim=",",
    )
    body = if builder == :mirrored
        simplewing_mirrored(c.b, c.ar, c.tr, c.twist_root, c.twist_tip,
                            c.lambda, c.gamma; kwargs...)
    elseif builder == :simple
        simplewing(c.b, c.ar, c.tr, c.twist_root, c.twist_tip,
                   c.lambda, c.gamma; kwargs..., b_low=-1.0, b_up=1.0)
    else
        throw(ArgumentError("builder must be :mirrored or :simple; got $(builder)."))
    end

    wake_direction = reshape(c.Vinf ./ c.magVinf, :, 1)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
    end
    return body, c
end

function force_coefficients_from_pressure(body, pressure, c, Lhat, Dhat; correct_kutta::Bool)
    F = zeros(3, body.ncells)
    pnl.calcfield_F!(F, body, pnl.calc_areas(body), body.normals, pressure;
                     correct_kuttacondition=correct_kutta)
    Ftot = pnl.calcfield_Ftot(body, F)
    qS = 0.5 * c.rho * c.magVinf^2 * c.Sref
    return LA.dot(Ftot, Lhat) / qS, LA.dot(Ftot, Dhat) / qS, F
end

function run_pressure_case(; builder::Symbol, n_rfl::Int, n_span::Int,
                           pressure::Symbol, gradient_mode::Symbol=:raw_hessian,
                           correct_kutta::Bool=false, bodyoptargs=(;))
    body, c = build_sweptwing(; builder, n_rfl, n_span, bodyoptargs)
    backend = pnl.DirectBackend()
    Dhat = c.Vinf / LA.norm(c.Vinf)
    Shat = [0.0, 1.0, 0.0]
    Lhat = LA.cross(Dhat, Shat)
    frames = pnl.ReferenceFrame(body)
    normalization = pnl.WingNormalization(c.rho, c.Sref, c.c_ref)

    pressure_monitor = pressure == :bernoulli ?
        pnl.PressureBernoulli(c.rho; correct_kuttacondition=correct_kutta) :
        pnl.PressureLaplace((body,), c.rho; reference_panel=1,
            reference_pressure=0.0, gradient_mode=gradient_mode, verbose=false)
    force_monitor = pnl.ForceMonitor(1, 1; i_frame=-1, normalization=normalization,
        correct_kuttacondition=correct_kutta, verbose=false)

    pnl.steady!(body, frames, c.Vinf;
        body_solvers=pnl.Backslash(body),
        backend=backend,
        monitors=(pressure_monitor, force_monitor),
        path=nothing,
        verbose=false)

    p = pressure == :bernoulli ? pressure_monitor.pressure[1] : pressure_monitor.p[1]
    CL_direct, CD_direct, F_direct =
        force_coefficients_from_pressure(body, p, c, Lhat, Dhat; correct_kutta)
    Fmon = force_monitor.force[:, 1]
    return (; builder, n_rfl, n_span, panels=body.ncells, pressure, gradient_mode,
        correct_kutta, body, constants=c, Lhat, Dhat, pressure_field=copy(p),
        CL=LA.dot(Fmon, Lhat), CD=LA.dot(Fmon, Dhat),
        CL_direct, CD_direct,
        force_monitor_diff=LA.norm(force_monitor.distributed_force .- F_direct) /
            max(LA.norm(F_direct), eps()))
end

function half_lift(case)
    y = view(case.body.controlpoints, 2, :)
    F = zeros(3, case.body.ncells)
    pnl.calcfield_F!(F, case.body, pnl.calc_areas(case.body), case.body.normals,
                     case.pressure_field; correct_kuttacondition=case.correct_kutta)
    lift = [LA.dot(view(F, :, i), case.Lhat) for i in 1:case.body.ncells]
    return (; positive=sum(lift[y .> 0]), negative=sum(lift[y .< 0]))
end

function mirror_pair_stats(body)
    y = view(body.controlpoints, 2, :)
    pos = findall(y .> 1e-9)
    neg = findall(y .< -1e-9)
    used = falses(length(neg))
    distances = Float64[]
    normal_dots = Float64[]

    for i in pos
        target = (body.controlpoints[1, i], -body.controlpoints[2, i],
                  body.controlpoints[3, i])
        best_k = 0
        best_d = Inf
        for (k, j) in enumerate(neg)
            used[k] && continue
            d = sqrt((body.controlpoints[1, j] - target[1])^2 +
                     (body.controlpoints[2, j] - target[2])^2 +
                     (body.controlpoints[3, j] - target[3])^2)
            if d < best_d
                best_d = d
                best_k = k
            end
        end
        best_k == 0 && continue
        j = neg[best_k]
        used[best_k] = true
        push!(distances, best_d)
        mirrored_normal = (body.normals[1, j], -body.normals[2, j], body.normals[3, j])
        push!(normal_dots, LA.dot(view(body.normals, :, i), mirrored_normal))
    end

    return (; n=length(distances),
        pair_median=median(distances),
        pair_max=maximum(distances),
        normal_dot_min=minimum(normal_dots),
        normal_dot_p01=quantile(normal_dots, 0.01),
        normal_dot_median=median(normal_dots))
end

function explicit_surface_velocity(solved_body, Vinf, backend;
                                   scale::Real=0.5,
                                   basis::Symbol=:quad,
                                   tri_robust::Bool=false)
    body = deepcopy(solved_body)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    fill!(body.velocity, 0.0)
    pnl.influence!((body,), (body,), backend; scalar_potential=false, velocity=true)

    pnl.compute_mu_gradient!(body.velocity, body.controlpoints, body.normals,
        body.cells, body.neighbor, view(body.strength, :, pnl.get_Gammai(body)),
        view(body.shedding_full, 1:2, :);
        scale=scale,
        nodes=body.nodes,
        grad_mu_options=basis === :tri ? (; basis, tri_robust) : (; basis))
    pnl.apply_freestream!(body, Vinf)
    return body.velocity
end

function current_surface_velocity(solved_body, Vinf, backend;
                                  basis::Symbol=:quad,
                                  tri_robust::Bool=false)
    body = deepcopy(solved_body)
    pnl.calcfield_U!(body, Vinf; backend,
        grad_mu_options=basis === :tri ? (; basis, tri_robust) : (; basis))
    return body.velocity
end

function pressure_force_from_velocity(body, velocity, c, Lhat, Dhat; correct_kutta::Bool=false)
    p = zeros(body.ncells)
    pnl.calcfield_P!(p, body, velocity, c.magVinf, c.rho, nothing;
                     correct_kuttacondition=correct_kutta)
    CL, CD, _ = force_coefficients_from_pressure(body, p, c, Lhat, Dhat; correct_kutta)
    return CL, CD
end

function normal_residual(body, velocity, magVinf)
    vals = [LA.dot(view(velocity, :, i), view(body.normals, :, i)) for i in 1:body.ncells]
    return sqrt(sum(abs2, vals) / length(vals)) / magVinf, maximum(abs, vals) / magVinf
end

function reconstruction_metrics(case; scale::Real=0.5,
                                basis::Symbol=:quad,
                                tri_robust::Bool=false)
    backend = pnl.DirectBackend()
    U_current = current_surface_velocity(case.body, case.constants.Vinf, backend;
        basis, tri_robust)
    U_explicit = explicit_surface_velocity(case.body, case.constants.Vinf, backend;
        scale, basis, tri_robust)
    d = vec(U_current .- U_explicit)
    mag = [LA.norm(view(U_current .- U_explicit, :, i)) for i in 1:case.body.ncells]
    nrm = LA.norm(U_current)
    normal_rms_current, normal_max_current =
        normal_residual(case.body, U_current, case.constants.magVinf)
    normal_rms_explicit, normal_max_explicit =
        normal_residual(case.body, U_explicit, case.constants.magVinf)
    CL_current, CD_current =
        pressure_force_from_velocity(case.body, U_current, case.constants, case.Lhat, case.Dhat)
    CL_explicit, CD_explicit =
        pressure_force_from_velocity(case.body, U_explicit, case.constants, case.Lhat, case.Dhat)
    return (; scale, basis, tri_robust,
        rel=LA.norm(d) / max(nrm, eps()),
        maxdiff=maximum(mag),
        p95=quantile(mag, 0.95),
        median=median(mag),
        normal_rms_current, normal_max_current,
        normal_rms_explicit, normal_max_explicit,
        CL_current, CD_current, CL_explicit, CD_explicit)
end

function print_case_table(rows)
    println("\n#===== Pressure/force sweep =====#")
    @printf("%-9s %3s %4s %6s %-9s %-16s %-5s %10s %10s %10s %10s %9s\n",
        "builder", "rfl", "span", "panels", "pressure", "grad_mode", "kutta",
        "CL", "CD", "CLerr%", "CLdir", "Frel")
    for r in rows
        @printf("%-9s %3d %4d %6d %-9s %-16s %-5s %+10.5f %+10.5f %10.2f %+10.5f %9.2e\n",
            string(r.builder), r.n_rfl, r.n_span, r.panels, string(r.pressure),
            string(r.gradient_mode), string(r.correct_kutta), r.CL, r.CD,
            100 * (r.CL - CLexp) / CLexp, r.CL_direct, r.force_monitor_diff)
    end
end

function print_reconstruction_table(rows)
    println("\n#===== Surface velocity reconstruction sweep =====#")
    @printf("%5s %-5s %-6s %11s %11s %11s %11s %11s %11s %10s %10s\n",
        "scale", "basis", "tri_rb", "relU", "max|dU|", "p95|dU|",
        "nRMS cur", "nRMS exp", "CL cur", "CL exp", "dCL")
    for r in rows
        @printf("%5.2f %-5s %-6s %11.3e %11.3e %11.3e %11.3e %11.3e %+11.5f %+10.5f %+10.5f\n",
            r.scale, string(r.basis), string(r.tri_robust), r.rel,
            r.maxdiff, r.p95, r.normal_rms_current, r.normal_rms_explicit,
            r.CL_current, r.CL_explicit, r.CL_explicit - r.CL_current)
    end
end

function print_winding_check()
    println("\n#===== Mirrored winding / matched-topology check =====#")
    rows = (
        ("simple nspan=32", run_pressure_case(; builder=:simple, n_rfl=4, n_span=32,
            pressure=:bernoulli, correct_kutta=false)),
        ("mirrored nspan=16", run_pressure_case(; builder=:mirrored, n_rfl=4, n_span=16,
            pressure=:bernoulli, correct_kutta=false)),
        ("mirrored ensure", run_pressure_case(; builder=:mirrored, n_rfl=4, n_span=16,
            pressure=:bernoulli, correct_kutta=false,
            bodyoptargs=(; ensure_winding=true))),
    )
    @printf("%-20s %6s %6s %10s %10s %10s %10s\n",
        "case", "panels", "shed", "CL", "CLerr%", "L y>0", "L y<0")
    for (label, row) in rows
        halves = half_lift(row)
        @printf("%-20s %6d %6d %+10.5f %10.2f %+10.3f %+10.3f\n",
            label, row.panels, size(row.body.shedding[1], 2), row.CL,
            100 * (row.CL - CLexp) / CLexp, halves.positive, halves.negative)
    end
    println("`ensure_winding=true` leaves the mirrored result unchanged; local panel winding is not the CL deficit source.")
    println("The matched simplewing case is not left/right symmetric, so its CL agreement is not a valid mirrored-wing reference.")

    println("\n#===== Simplewing mesh mirror-pair check =====#")
    @printf("%-20s %8s %12s %12s %12s %12s\n",
        "case", "pairs", "pair med", "pair max", "n dot min", "n dot p01")
    for (label, row) in rows[1:2]
        stats = mirror_pair_stats(row.body)
        @printf("%-20s %8d %12.4e %12.4e %12.4f %12.4f\n",
            label, stats.n, stats.pair_median, stats.pair_max,
            stats.normal_dot_min, stats.normal_dot_p01)
    end
    println("The full-span simplewing surface is symmetric as a loft, but its triangular panels are not mirror-paired.")
    println("That panel-level asymmetry is large enough to produce different pressure/circulation on the two halves.")
end

function main()
    full = lowercase(get(ENV, "FLOWPANEL_SWEPTWING_FULL_DIAG", "false")) in ("1", "true", "yes")
    resolutions = full ? [(4, 16), (6, 24), (8, 40)] : [(3, 10), (4, 16)]
    builders = (:mirrored, :simple)
    kutta_options = (false, true)
    pressure_variants = ((:bernoulli, :raw_hessian),
                         (:laplace, :raw_hessian),
                         (:laplace, :surface_velocity))

    cases = []
    for builder in builders, (n_rfl, n_span) in resolutions,
        correct_kutta in kutta_options, (pressure, gradient_mode) in pressure_variants
        push!(cases, run_pressure_case(; builder, n_rfl, n_span,
            pressure, gradient_mode, correct_kutta))
    end
    print_case_table(cases)

    baseline = first(r for r in cases if r.builder == :mirrored &&
        r.n_rfl == resolutions[end][1] && r.n_span == resolutions[end][2] &&
        r.pressure == :bernoulli && !r.correct_kutta)

    recon_rows = []
    for basis in (:quad, :tri), tri_robust in (false, true),
        scale in (0.0, 0.25, 0.5, 0.75, 1.0)
        basis === :quad && tri_robust && continue
        push!(recon_rows, reconstruction_metrics(baseline;
            scale, basis, tri_robust))
    end
    print_reconstruction_table(recon_rows)

    print_winding_check()

    c = baseline.constants
    println("\n#===== Normalization check =====#")
    @printf("Sref=b^2/AR=%.8f  c_ref=b/AR=%.8f  q=0.5*rho*V^2=%.8f\n",
        c.Sref, c.c_ref, 0.5 * c.rho * c.magVinf^2)
    @printf "rho cancels in CL when pressure and force normalization both use rho consistently.\n"
end

main()
