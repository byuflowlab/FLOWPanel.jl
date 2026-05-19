#=##############################################################################
# Diagnostic: swept-wing PressureLaplace pressure-Poisson consistency check.
#
# This script intentionally leaves sweptwing.jl untouched.  It builds the same
# swept wing, compares Bernoulli pressure against PressureLaplace, and prints
# compact operator/RHS evidence instead of plots or VTK files.
=###############################################################################

import FLOWPanel as pnl
include(joinpath(pnl.examples_path, "helper_functions.jl"))

import LinearAlgebra
import SparseArrays
using Printf
using Statistics: mean, median, quantile

const LA = LinearAlgebra
const _norm = LinearAlgebra.norm
const _dot = LinearAlgebra.dot
const _cross = LinearAlgebra.cross

function tangent_velocity(body)
    out = similar(body.velocity)
    @inbounds for p in 1:body.ncells
        n = view(body.normals, :, p)
        u = view(body.velocity, :, p)
        un = _dot(u, n)
        out[:, p] .= u .- un .* n .+ body.velocity_kinematic[:, p]
    end
    return out
end

function boundary_panel_mask(cells::AbstractMatrix{Int})
    edge_counts = Dict{Tuple{Int,Int}, Int}()
    @inbounds for p in axes(cells, 2)
        nnodes = size(cells, 1)
        for a in 1:nnodes
            n1 = cells[a, p]
            n2 = cells[a == nnodes ? 1 : a + 1, p]
            edge = n1 < n2 ? (n1, n2) : (n2, n1)
            edge_counts[edge] = get(edge_counts, edge, 0) + 1
        end
    end

    mask = falses(size(cells, 2))
    @inbounds for p in axes(cells, 2)
        nnodes = size(cells, 1)
        for a in 1:nnodes
            n1 = cells[a, p]
            n2 = cells[a == nnodes ? 1 : a + 1, p]
            edge = n1 < n2 ? (n1, n2) : (n2, n1)
            if edge_counts[edge] == 1
                mask[p] = true
                break
            end
        end
    end
    return mask
end

function force_coefficients_from_pressure(body, pressure, Lhat, Dhat, rho, magVinf, Sref)
    body_x = deepcopy(body)
    fill!(body_x.F, 0.0)
    body_x.P .= pressure
    pnl.calcfield_F!(body_x; correct_kuttacondition=false)
    LDS = pnl.calcfield_LDS(body_x, Lhat, Dhat)
    qS = 0.5 * rho * magVinf^2 * Sref
    CL = sign(_dot(LDS[:, 1], Lhat)) * _norm(LDS[:, 1]) / qS
    CD = sign(_dot(LDS[:, 2], Dhat)) * _norm(LDS[:, 2]) / qS
    return CL, CD
end

function pressure_range(p)
    return minimum(p), maximum(p), median(p), quantile(p, 0.01), quantile(p, 0.99)
end

function print_pressure_range(label, p)
    pmin, pmax, p50, p01, p99 = pressure_range(p)
    @printf "%-26s min=%+11.4e  max=%+11.4e  median=%+11.4e  p01=%+11.4e  p99=%+11.4e\n" label pmin pmax p50 p01 p99
end

function rhs_from_acceleration(pl, body, acceleration)
    b = zeros(body.ncells)
    @inbounds for k in axes(pl.edges[1], 2)
        edge_a, edge_b, i, j = pl.edges[1][:, k]
        w = pnl._pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]
        r1 = body.controlpoints[1, j] - body.controlpoints[1, i]
        r2 = body.controlpoints[2, j] - body.controlpoints[2, i]
        r3 = body.controlpoints[3, j] - body.controlpoints[3, i]
        aedge_dot_r = 0.5 * (
            (acceleration[1, i] + acceleration[1, j]) * r1 +
            (acceleration[2, i] + acceleration[2, j]) * r2 +
            (acceleration[3, i] + acceleration[3, j]) * r3)
        flux = pl.rho * w * aedge_dot_r
        b[i] += flux
        b[j] -= flux
    end
    b[pl.reference_panel] = pl.reference_pressure
    return b
end

function rhs_from_edge_material_derivative(pl, body, velocity_dot)
    b = zeros(body.ncells)
    @inbounds for k in axes(pl.edges[1], 2)
        edge_a, edge_b, i, j = pl.edges[1][:, k]
        w = pnl._pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]
        r1 = body.controlpoints[1, j] - body.controlpoints[1, i]
        r2 = body.controlpoints[2, j] - body.controlpoints[2, i]
        r3 = body.controlpoints[3, j] - body.controlpoints[3, i]
        udot_edge_dot_r = 0.5 * (
            (velocity_dot[1, i] + velocity_dot[1, j]) * r1 +
            (velocity_dot[2, i] + velocity_dot[2, j]) * r2 +
            (velocity_dot[3, i] + velocity_dot[3, j]) * r3)

        ni = view(body.normals, :, i)
        nj = view(body.normals, :, j)
        ui = view(body.velocity, :, i)
        uj = view(body.velocity, :, j)
        urel_i = ui .- _dot(ui, ni) .* ni
        urel_j = uj .- _dot(uj, nj) .* nj
        du = view(body.velocity, :, j) .- view(body.velocity, :, i)
        convective_edge_dot_r = _dot(0.5 .* (urel_i .+ urel_j), du)

        flux = pl.rho * w * (udot_edge_dot_r + convective_edge_dot_r)
        b[i] += flux
        b[j] -= flux
    end
    b[pl.reference_panel] = pl.reference_pressure
    return b
end

function rhs_from_conormal_acceleration(pl, body, acceleration)
    b = zeros(body.ncells)
    @inbounds for k in axes(pl.edges[1], 2)
        edge_a, edge_b, i, j = pl.edges[1][:, k]
        w, ell, nu1, nu2, nu3, n1, n2, n3 =
            pnl._pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)
        ai_n = acceleration[1, i] * n1 + acceleration[2, i] * n2 + acceleration[3, i] * n3
        aj_n = acceleration[1, j] * n1 + acceleration[2, j] * n2 + acceleration[3, j] * n3
        ai = (acceleration[1, i] - ai_n * n1) * nu1 +
             (acceleration[2, i] - ai_n * n2) * nu2 +
             (acceleration[3, i] - ai_n * n3) * nu3
        aj = (acceleration[1, j] - aj_n * n1) * nu1 +
             (acceleration[2, j] - aj_n * n2) * nu2 +
             (acceleration[3, j] - aj_n * n3) * nu3
        flux = ell * 0.5 * (ai + aj)
        b[i] += pl.rho * flux
        b[j] -= pl.rho * flux
        _ = w
    end
    b[pl.reference_panel] = pl.reference_pressure
    return b
end

function rhs_from_bernoulli_edge_flux(pl, body, p)
    b = zeros(body.ncells)
    @inbounds for k in axes(pl.edges[1], 2)
        edge_a, edge_b, i, j = pl.edges[1][:, k]
        w = pnl._pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]
        flux = w * (p[i] - p[j])
        b[i] += flux
        b[j] -= flux
    end
    b[pl.reference_panel] = p[pl.reference_panel]
    return b
end

function convective_acceleration(body, velocity_dot, G; transpose_gradient::Bool=false)
    acc = copy(velocity_dot)
    u_t = tangent_velocity(body)
    @inbounds for p in 1:body.ncells
        u1, u2, u3 = u_t[1, p], u_t[2, p], u_t[3, p]
        if transpose_gradient
            acc[1, p] += G[1, 1, p] * u1 + G[2, 1, p] * u2 + G[3, 1, p] * u3
            acc[2, p] += G[1, 2, p] * u1 + G[2, 2, p] * u2 + G[3, 2, p] * u3
            acc[3, p] += G[1, 3, p] * u1 + G[2, 3, p] * u2 + G[3, 3, p] * u3
        else
            acc[1, p] += G[1, 1, p] * u1 + G[1, 2, p] * u2 + G[1, 3, p] * u3
            acc[2, p] += G[2, 1, p] * u1 + G[2, 2, p] * u2 + G[2, 3, p] * u3
            acc[3, p] += G[3, 1, p] * u1 + G[3, 2, p] * u2 + G[3, 3, p] * u3
        end
    end
    return acc
end

function frobenius_per_panel(G)
    return [sqrt(sum(abs2, view(G, :, :, p))) for p in axes(G, 3)]
end

function build_sweptwing()
    AOA = 4.2
    magVinf = 30.0
    Vinf = magVinf * [cosd(AOA), 0.0, sind(AOA)]
    rho = 1.225

    b = 98 * 0.0254
    ar = 5.0
    tr = 1.0
    lambda = 45
    airfoil = "airfoil-rae101.csv"
    airfoil_path = joinpath(pnl.examples_path, "data")

    n_rfl = 8
    NDIVS_rfl = [(0.25, n_rfl, 10.0, false),
                 (0.50, n_rfl, 1.0, true),
                 (0.25, n_rfl, 1 / 10.0, false)]
    NDIVS_span = [(1.0, 30, 20.0, true)]

    bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
    body = simplewing(b, ar, tr, 0, 0, lambda, 0;
        bodytype=bodytype,
        bodyoptargs=(; CPoffset=1e-14),
        airfoil_root=airfoil,
        airfoil_tip=airfoil,
        airfoil_path=airfoil_path,
        rfl_NDIVS=NDIVS_rfl,
        delim=",",
        span_NDIVS=NDIVS_span,
        b_low=-1.0,
        b_up=1.0)

    wake_direction = reshape(Vinf ./ magVinf, :, 1)
    for i in eachindex(body.Das)
        body.Das[i] .= repeat(wake_direction, 1, size(body.Das[i], 2))
    end

    return body, Vinf, magVinf, rho, b^2 / ar, b / ar
end

function run_laplace_case(body, Vinf, rho, Sref, c_ref, mode::Symbol; unsteady::Bool=false)
    backend = pnl.DirectBackend()
    body_l = deepcopy(body)
    body_l.needs_velocity_gradient[] = true
    body_l.velocity .= 0.0
    body_l.velocity_gradient .= 0.0
    pnl.calcfield_U!(body_l, Vinf; backend)
    pnl.influence!((body_l,), (body_l,), backend;
        scalar_potential=false,
        velocity=false,
        velocity_gradient=true)

    pl = pnl.PressureLaplace((body_l,), rho;
        reference_panel=1,
        reference_pressure=0.0,
        verbose=false,
        gradient_mode=mode,
        unsteady=unsteady)
    fm = pnl.ForceMonitor(1, 1;
        i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c_ref),
        correct_kuttacondition=false,
        verbose=false)

    pl.velocity_dot[1] .= unsteady ? 0.0 : -tangent_velocity(body_l)
    pl((body_l,), (nothing,), pnl.ReferenceFrame(body_l), Vinf, 0, 1.0)
    fm((body_l,), (nothing,), pnl.ReferenceFrame(body_l), Vinf, 0, 1.0)

    return body_l, pl, fm.force[:, 1]
end

function summarize_operator(label, body_l, pl, p_b_shift, Lhat, Dhat, rho, magVinf, Sref, boundary_mask)
    L = pl.L[1]
    b_rhs = pl.b[1]
    p_direct = L \ b_rhs
    lp_b = L * p_b_shift
    residual = lp_b .- b_rhs
    cg_residual = L * pl.p[1] .- b_rhs
    direct_residual = L * p_direct .- b_rhs

    CL_direct, CD_direct = force_coefficients_from_pressure(body_l, p_direct, Lhat, Dhat, rho, magVinf, Sref)
    @printf "\n#===== %s =====#\n" label
    @printf "direct-solve force:       CL=%+11.5f  CD=%+11.5f\n" CL_direct CD_direct
    @printf "norms:                    ||L*pB||=%11.4e  ||rhs||=%11.4e  ||L*pB-rhs||=%11.4e  rel=%9.3e\n" _norm(lp_b) _norm(b_rhs) _norm(residual) (_norm(residual) / max(_norm(b_rhs), eps()))
    @printf "solver residuals:         CG rel=%9.3e  direct rel=%9.3e\n" (_norm(cg_residual) / max(_norm(b_rhs), eps())) (_norm(direct_residual) / max(_norm(b_rhs), eps()))
    @printf "boundary residual rel:    boundary=%9.3e  interior=%9.3e  panels(boundary/interior)=%d/%d\n" (
        _norm(residual[boundary_mask]) / max(_norm(b_rhs[boundary_mask]), eps())) (
        _norm(residual[.!boundary_mask]) / max(_norm(b_rhs[.!boundary_mask]), eps())) sum(boundary_mask) sum(.!boundary_mask)
    print_pressure_range("Laplace pressure", pl.p[1])
    return lp_b, residual
end

function main()
    println("Building swept-wing pressure-Poisson diagnostic...")
    body, Vinf, magVinf, rho, Sref, c_ref = build_sweptwing()
    backend = pnl.DirectBackend()
    pnl.apply_freestream!(body, Vinf)
    pnl.solve!(body, pnl.Backslash(body); backend)
    pnl.calcfield_U!(body, Vinf; backend)
    pnl.calcfield_P!(body, magVinf, rho; correct_kuttacondition=false)
    pnl.calcfield_F!(body; correct_kuttacondition=false)

    Dhat = Vinf / _norm(Vinf)
    Shat = [0.0, 1.0, 0.0]
    Lhat = _cross(Dhat, Shat)
    CL_b, CD_b = force_coefficients_from_pressure(body, body.P, Lhat, Dhat, rho, magVinf, Sref)
    p_b = copy(body.P)
    p_b_shift = p_b .- p_b[1]
    boundary_mask = boundary_panel_mask(body.cells)

    println("\n#===== setup =====#")
    @printf "panels=%d  boundary panels=%d  interior panels=%d\n" body.ncells sum(boundary_mask) sum(.!boundary_mask)
    @printf "Vinf=(%+.6f,%+.6f,%+.6f)  rho=%.4f  Sref=%.6f  cref=%.6f\n" Vinf[1] Vinf[2] Vinf[3] rho Sref c_ref
    @printf "PressureLaplace primary runs use unsteady=false; velocity_dot still updates after each call\n"

    println("\n#===== Bernoulli baseline, no Kutta pressure averaging =====#")
    @printf "Bernoulli force:          CL=%+11.5f  CD=%+11.5f\n" CL_b CD_b
    print_pressure_range("Bernoulli pressure", p_b)

    body_sv, pl_sv, F_sv = run_laplace_case(body, Vinf, rho, Sref, c_ref, :surface_velocity; unsteady=false)
    body_raw, pl_raw, F_raw = run_laplace_case(body, Vinf, rho, Sref, c_ref, :raw_hessian)
    body_sv_u, pl_sv_u, F_sv_u = run_laplace_case(body, Vinf, rho, Sref, c_ref, :surface_velocity; unsteady=true)

    @printf "\n#===== monitor forces =====#\n"
    @printf "%-24s CL=%+11.5f  CD=%+11.5f\n" "surface_velocity" _dot(F_sv, Lhat) _dot(F_sv, Dhat)
    @printf "%-24s CL=%+11.5f  CD=%+11.5f\n" "raw_hessian" _dot(F_raw, Lhat) _dot(F_raw, Dhat)
    @printf "%-24s CL=%+11.5f  CD=%+11.5f\n" "surface_velocity u" _dot(F_sv_u, Lhat) _dot(F_sv_u, Dhat)

    lp_b_sv, residual_sv = summarize_operator("surface_velocity RHS", body_sv, pl_sv, p_b_shift, Lhat, Dhat, rho, magVinf, Sref, boundary_mask)
    summarize_operator("raw_hessian RHS", body_raw, pl_raw, p_b_shift, Lhat, Dhat, rho, magVinf, Sref, boundary_mask)

    f_raw = frobenius_per_panel(body_raw.velocity_gradient)
    f_sv = frobenius_per_panel(pl_sv.surface_velocity_gradient[1])
    println("\n#===== gradient scale check =====#")
    @printf "%-22s median=%11.4e  p90=%11.4e  p99=%11.4e  max=%11.4e\n" "raw_hessian" median(f_raw) quantile(f_raw, 0.90) quantile(f_raw, 0.99) maximum(f_raw)
    @printf "%-22s median=%11.4e  p90=%11.4e  p99=%11.4e  max=%11.4e\n" "surface_velocity" median(f_sv) quantile(f_sv, 0.90) quantile(f_sv, 0.99) maximum(f_sv)
    @printf "%-22s median=%11.4e\n" "surface/raw ratio" (median(f_sv) / max(median(f_raw), eps()))

    b_edge = rhs_from_edge_material_derivative(pl_sv, body_sv, zeros(size(body_sv.velocity)))
    b_gu = rhs_from_acceleration(pl_sv, body_sv, convective_acceleration(body_sv, zeros(size(body_sv.velocity)), pl_sv.surface_velocity_gradient[1]; transpose_gradient=false))
    b_gtu = rhs_from_acceleration(pl_sv, body_sv, convective_acceleration(body_sv, zeros(size(body_sv.velocity)), pl_sv.surface_velocity_gradient[1]; transpose_gradient=true))
    b_conormal = rhs_from_conormal_acceleration(pl_sv, body_sv, pl_sv.acceleration[1])
    println("\n#===== acceleration RHS residuals, surface_velocity gradient only =====#")
    @printf "%-24s rel ||L*pB-rhs||/||rhs|| = %9.3e\n" "edge material" (_norm(lp_b_sv .- b_edge) / max(_norm(b_edge), eps()))
    @printf "%-24s rel ||L*pB-rhs||/||rhs|| = %9.3e\n" "G*u" (_norm(lp_b_sv .- b_gu) / max(_norm(b_gu), eps()))
    @printf "%-24s rel ||L*pB-rhs||/||rhs|| = %9.3e\n" "transpose(G)*u" (_norm(lp_b_sv .- b_gtu) / max(_norm(b_gtu), eps()))
    @printf "%-24s rel ||L*pB-rhs||/||rhs|| = %9.3e\n" "old conormal flux" (_norm(lp_b_sv .- b_conormal) / max(_norm(b_conormal), eps()))
    @printf "%-24s rel ||current rhs - edge|| = %9.3e\n" "current vs edge" (_norm(pl_sv.b[1] .- b_edge) / max(_norm(pl_sv.b[1]), eps()))

    b_mimetic = rhs_from_bernoulli_edge_flux(pl_sv, body_sv, p_b_shift)
    p_mimetic = pl_sv.L[1] \ b_mimetic
    CL_mim, CD_mim = force_coefficients_from_pressure(body_sv, p_mimetic, Lhat, Dhat, rho, magVinf, Sref)
    mimetic_residual = pl_sv.L[1] * p_b_shift .- b_mimetic
    println("\n#===== mimetic Bernoulli edge-flux RHS =====#")
    @printf "operator check:           ||L*pB-b_mimetic||/||b_mimetic|| = %9.3e\n" (_norm(mimetic_residual) / max(_norm(b_mimetic), eps()))
    @printf "direct solve vs pB:       ||p_mimetic-pB||/||pB|| = %9.3e\n" (_norm(p_mimetic .- p_b_shift) / max(_norm(p_b_shift), eps()))
    @printf "mimetic force:            CL=%+11.5f  CD=%+11.5f\n" CL_mim CD_mim

    println("\n#===== conclusion cue =====#")
    @printf "The current PressureLaplace RHS uses material acceleration with the same two-point edge weights as L.\n"
    @printf "The old conormal-flux and kinetic-pressure RHS values are diagnostic references only.\n"
    @printf "surface_velocity residual mean abs=%11.4e, median abs=%11.4e\n" mean(abs.(residual_sv)) median(abs.(residual_sv))
end

main()
