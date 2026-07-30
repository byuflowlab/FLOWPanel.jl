import FLOWPanel as pnl

using LinearAlgebra: dot, norm
using Statistics: median, quantile

include(joinpath(pnl.examples_path, "simple_wing_capped_pressure_comparison.jl"))

function surface_gradient_old!(grad_out, controlpoints, normals, cells,
                               neighbors, scalar, te_info,
                               ATA, ATb, stencil)
    n = size(cells, 2)
    @inbounds for i in 1:n
        empty!(stencil)

        is_te = te_info !== nothing && te_info[1, i] > 0 && te_info[2, i] > 0
        te_v1 = is_te ? cells[te_info[1, i], i] : -1
        te_v2 = is_te ? cells[te_info[2, i], i] : -1

        for k in 1:3
            j = neighbors[k, i]
            j <= 0 && continue
            if is_te
                has_v1 = cells[1, j] == te_v1 || cells[2, j] == te_v1 || cells[3, j] == te_v1
                has_v2 = cells[1, j] == te_v2 || cells[2, j] == te_v2 || cells[3, j] == te_v2
                has_v1 && has_v2 && continue
            end
            push!(stencil, j)
        end

        if is_te || length(stencil) < 3
            n_current = length(stencil)
            for s in 1:n_current
                s_idx = stencil[s]
                for k in 1:3
                    nn_idx = neighbors[k, s_idx]
                    nn_idx <= 0 && continue
                    nn_idx == i && continue
                    nn_idx in stencil && continue
                    if is_te
                        has_v1 = cells[1, nn_idx] == te_v1 || cells[2, nn_idx] == te_v1 || cells[3, nn_idx] == te_v1
                        has_v2 = cells[1, nn_idx] == te_v2 || cells[2, nn_idx] == te_v2 || cells[3, nn_idx] == te_v2
                        has_v1 && has_v2 && continue
                    end
                    push!(stencil, nn_idx)
                end
            end
        end

        fill!(ATA, 0.0)
        fill!(ATb, 0.0)
        mean_sq_dist = 0.0
        for j in stencil
            dx = controlpoints[1, j] - controlpoints[1, i]
            dy = controlpoints[2, j] - controlpoints[2, i]
            dz = controlpoints[3, j] - controlpoints[3, i]
            ds = scalar[j] - scalar[i]
            ATA[1, 1] += dx * dx; ATA[1, 2] += dx * dy; ATA[1, 3] += dx * dz
            ATA[2, 1] += dy * dx; ATA[2, 2] += dy * dy; ATA[2, 3] += dy * dz
            ATA[3, 1] += dz * dx; ATA[3, 2] += dz * dy; ATA[3, 3] += dz * dz
            ATb[1] += dx * ds
            ATb[2] += dy * ds
            ATb[3] += dz * ds
            mean_sq_dist += dx^2 + dy^2 + dz^2
        end
        mean_sq_dist = isempty(stencil) ? 1.0 : mean_sq_dist / length(stencil)

        nx, ny, nz = normals[1, i], normals[2, i], normals[3, i]
        penalty = 1e4 * mean_sq_dist
        ATA[1, 1] += penalty * nx * nx; ATA[1, 2] += penalty * nx * ny; ATA[1, 3] += penalty * nx * nz
        ATA[2, 1] += penalty * ny * nx; ATA[2, 2] += penalty * ny * ny; ATA[2, 3] += penalty * ny * nz
        ATA[3, 1] += penalty * nz * nx; ATA[3, 2] += penalty * nz * ny; ATA[3, 3] += penalty * nz * nz

        reg = 1e-10 * mean_sq_dist
        ATA[1, 1] += reg; ATA[2, 2] += reg; ATA[3, 3] += reg

        g = ATA \ ATb
        grad_out[1, i] = g[1]
        grad_out[2, i] = g[2]
        grad_out[3, i] = g[3]
    end
    return grad_out
end

function old_surface_ls_acceleration(body)
    n = body.ncells
    u_inertial = zeros(Float64, 3, n)
    acceleration = zeros(Float64, 3, n)
    grad = zeros(Float64, 3, 3, n)
    grad_comp = zeros(Float64, 3, n)
    ATA = zeros(Float64, 3, 3)
    ATb = zeros(Float64, 3)
    stencil = Int[]
    te_info = hasproperty(body, :shedding_full) ? view(body.shedding_full, 1:2, :) : nothing

    @inbounds for p in 1:n
        nx, ny, nz = body.normals[1, p], body.normals[2, p], body.normals[3, p]
        ux, uy, uz = body.velocity[1, p], body.velocity[2, p], body.velocity[3, p]
        un = ux * nx + uy * ny + uz * nz
        u_inertial[1, p] = ux - un * nx + body.velocity_kinematic[1, p]
        u_inertial[2, p] = uy - un * ny + body.velocity_kinematic[2, p]
        u_inertial[3, p] = uz - un * nz + body.velocity_kinematic[3, p]
    end

    for comp in 1:3
        fill!(grad_comp, 0.0)
        surface_gradient_old!(grad_comp, body.controlpoints, body.normals,
            body.cells, body.neighbor, view(u_inertial, comp, :), te_info,
            ATA, ATb, stencil)
        @inbounds for i in 1:n, k in 1:3
            grad[comp, k, i] = grad_comp[k, i]
        end
    end

    @inbounds for i in 1:n
        ur1, ur2, ur3 = body.velocity[1, i], body.velocity[2, i], body.velocity[3, i]
        nx, ny, nz = body.normals[1, i], body.normals[2, i], body.normals[3, i]
        un = ur1 * nx + ur2 * ny + ur3 * nz
        ut1 = ur1 - un * nx
        ut2 = ur2 - un * ny
        ut3 = ur3 - un * nz
        for comp in 1:3
            acceleration[comp, i] += ut1 * grad[comp, 1, i] +
                                     ut2 * grad[comp, 2, i] +
                                     ut3 * grad[comp, 3, i]
        end
    end

    return acceleration
end

function rhs_from_acceleration!(b, pressure_laplace, body, acceleration; rho)
    fill!(b, 0.0)
    edges = pressure_laplace.edges[1]
    reference_panel = pressure_laplace.reference_panel
    @inbounds for k in axes(edges, 2)
        edge_a, edge_b, i, j = edges[1, k], edges[2, k], edges[3, k], edges[4, k]
        _, ell, nu1, nu2, nu3, n1, n2, n3 =
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
        b[i] += rho * flux
        b[j] -= rho * flux
    end
    b[reference_panel] = pressure_laplace.reference_pressure
    return b
end

function force_coefficients(body, rho, Vinf)
    chord = maximum(body.nodes[1, :]) - minimum(body.nodes[1, :])
    span = maximum(body.nodes[2, :]) - minimum(body.nodes[2, :])
    normalization = pnl.WingNormalization(rho, chord * span, chord)
    monitor = pnl.ForceMonitor(1, 1; i_frame=-1, normalization,
        correct_kuttacondition=false, verbose=false)
    monitor((body,), (nothing,), pnl.ReferenceFrame(body), Vinf, 0, 1.0)
    return vec(monitor.force[:, 1])
end

function main()
    backend = pnl.DirectBackend()
    rho = 1.225
    body_ref = build_pressure_comparison_wing()
    Vinf = 56.0 .* [cosd(5.88), 0.0, sind(5.88)]
    set_pressure_comparison_wake!(body_ref, Vinf)
    pnl.steady!(body_ref, pnl.ReferenceFrame(body_ref), Vinf;
        body_solvers=pnl.Backslash(body_ref), backend, verbose=false)
    prepare_pressure_state!(body_ref, Vinf; backend)

    body_b = deepcopy(body_ref)
    body_h = deepcopy(body_ref)
    body_ls = deepcopy(body_ref)

    pressure_bernoulli = pnl.PressureBernoulli(rho; unsteady=false,
        correct_kuttacondition=false, backend=backend)
    pressure_bernoulli((body_b,), (nothing,), pnl.ReferenceFrame(body_b), Vinf, 0, 1.0)
    pB = body_b.P .- body_b.P[1]

    pressure_laplace = pnl.PressureLaplace((body_h,), rho;
        reference_panel=1, reference_pressure=0.0)
    pressure_laplace.velocity_dot[1] .= expected_negative_tangent_velocity(body_h)
    pnl._pressure_velocity_dot!(pressure_laplace, body_h, 1, 1.0)
    pnl._pressure_rhs!(pressure_laplace, body_h, 1)

    b_hessian = copy(pressure_laplace.b[1])
    LpB = pressure_laplace.L[1] * pB

    acceleration_ls = old_surface_ls_acceleration(body_ls)
    b_ls = similar(b_hessian)
    rhs_from_acceleration!(b_ls, pressure_laplace, body_ls, acceleration_ls; rho)

    p_hessian = pressure_laplace.L[1] \ b_hessian
    p_ls = pressure_laplace.L[1] \ b_ls
    body_h.P .= p_hessian
    body_ls.P .= p_ls

    r_h = LpB - b_hessian
    r_ls = LpB - b_ls

    println("npanels=", body_ref.ncells)
    println("rel residual LpB vs b_hessian = ", norm(r_h) / norm(b_hessian))
    println("rel residual LpB vs b_surface_ls = ", norm(r_ls) / norm(b_ls))
    println("cosine LpB,b_hessian = ", dot(LpB, b_hessian) / (norm(LpB) * norm(b_hessian)))
    println("cosine LpB,b_surface_ls = ", dot(LpB, b_ls) / (norm(LpB) * norm(b_ls)))
    println("||b_surface_ls - b_hessian|| / ||b_hessian|| = ", norm(b_ls - b_hessian) / norm(b_hessian))
    println("abs residual p99 hessian = ", quantile(abs.(r_h), 0.99))
    println("abs residual p99 surface_ls = ", quantile(abs.(r_ls), 0.99))
    println("force Bernoulli = ", force_coefficients(body_b, rho, Vinf))
    println("force Laplace hessian RHS = ", force_coefficients(body_h, rho, Vinf))
    println("force Laplace surface-LS RHS = ", force_coefficients(body_ls, rho, Vinf))
end

main()
