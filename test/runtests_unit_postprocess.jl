using Test
import FLOWPanel as pnl
using LinearAlgebra: dot
using SparseArrays

function make_planar_gradient_mesh()
    nodes = Float64[
        0 1 2 0 1 2 0 1 2;
        0 0 0 1 1 1 2 2 2;
        0 0 0 0 0 0 0 0 0;
    ]
    cells = Int[
        1 1 2 2 4 4 5 5;
        2 5 3 6 5 8 6 9;
        5 4 6 5 8 7 9 8;
    ]
    return nodes, cells
end

function make_postprocess_seeded_te_mesh()
    nodes = Float64[
        1 1 1 0 0 0 0 0;
        0 1 2 0 1 2 0 2;
        0 0 0 1 1 1 -1 -1;
    ]
    cells = Int[
        4 5 7 8;
        2 3 1 2;
        1 2 2 3;
    ]
    return nodes, cells
end

@testset verbose=true "Postprocess" begin
    @testset "PressureLaplace monitor" begin
        body = make_octa_source_body()
        monitor = pnl.PressureLaplace((body,), 1.2; reference_panel=1, reference_pressure=0.0)

        @test pnl.monitor_provides(monitor) == (:P,)
        @test pnl.audit_monitors((monitor, pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization()))) !== nothing
        @test_throws ArgumentError pnl.PressureLaplace(1.0)

        body.velocity .= 0.0
        body.velocity[1, :] .= 1.0
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.25)
        @test all(isapprox.(monitor.velocity_dot[1][1, :], -1.0; atol=1e-12))
        @test all(isapprox.(monitor.velocity_dot[1][2, :], 0.0; atol=1e-12))
        @test all(isapprox.(monitor.velocity_dot[1][3, :], 0.0; atol=1e-12))
        @test all(isfinite, body.P)

        body.velocity[1, :] .+= 0.5
        body.velocity[2, :] .-= 0.25
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 1, 0.25)
        @test all(isapprox.(monitor.velocity_dot[1][1, :], -1.5; atol=1e-12))
        @test all(isapprox.(monitor.velocity_dot[1][2, :], 0.25; atol=1e-12))
        @test all(isapprox.(monitor.velocity_dot[1][3, :], 0.0; atol=1e-12))
    end

    @testset "PressureLaplace sparse matrix and solve" begin
        body = make_octa_source_body()
        L = pnl._assemble_pressure_laplacian(body, 1)
        @test issparse(L)
        @test isapprox(Matrix(L), Matrix(L)'; atol=1e-12)
        @test L[1, 1] == 1.0
        @test all(L[1, 2:end] .== 0.0)
        @test all(L[2:end, 1] .== 0.0)

        L0 = pnl._assemble_pressure_laplacian(body, 0)
        @test all(isapprox.(vec(sum(L0; dims=2)), 0.0; atol=1e-12))

        p_exact = collect(range(0.0, 1.0; length=body.ncells))
        p_exact[1] = 0.0
        b = L * p_exact
        monitor = pnl.PressureLaplace((body,), 1.0)
        monitor.b[1] .= b
        pnl._pressure_solve!(monitor, 1)
        @test isapprox(monitor.p[1], p_exact; atol=1e-8)
    end

    @testset "PressureLaplace cache invalidation and force integration" begin
        body = make_octa_source_body()
        pressure = pnl.PressureLaplace((body,), 1.0; reference_panel=1)
        force = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())
        frames = pnl.ReferenceFrame(body)

        body.velocity[1, :] .= 0.2 .* (1:body.ncells)
        pressure((body,), (nothing,), frames, zeros(3), 0, 0.1)
        first_L = pressure.L[1]
        pressure((body,), (nothing,), frames, zeros(3), 1, 0.1)
        @test pressure.L[1] === first_L

        old = body.nodes[1, 1]
        body.nodes[1, 1] = old + 0.1
        pressure((body,), (nothing,), frames, zeros(3), 2, 0.1)
        @test pressure.L[1] !== first_L

        force((body,), (nothing,), frames, zeros(3), 0, 0.1)
        @test all(isfinite, body.P)
        @test all(isfinite, body.F)
        @test all(isfinite, force.force)
    end

    @testset "PressureLaplace Bernoulli constant-field comparison" begin
        body_b = make_octa_source_body()
        body_l = make_octa_source_body()
        body_b.velocity[1, :] .= 0.4
        body_b.velocity[2, :] .= -0.2
        body_l.velocity .= body_b.velocity

        laplace = pnl.PressureLaplace((body_l,), 1.0; reference_panel=1)
        laplace.velocity_dot[1] .= -body_l.velocity
        pnl.calcfield_P!(body_b.P, body_b, body_b.velocity, 1.0, 1.0, zeros(body_b.ncells);
            correct_kuttacondition=false)
        laplace((body_l,), (nothing,), pnl.ReferenceFrame(body_l), [1.0, 0.0, 0.0], 0, 0.1)

        p_b = body_b.P .- body_b.P[1]
        p_l = body_l.P .- body_l.P[1]
        @test isapprox(p_l, p_b; atol=1e-10)
    end

    @testset "compute_mu_gradient! interior recovery" begin
        nodes, cells = make_planar_gradient_mesh()
        body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells;
            watertight=false,
            ensure_winding=false)

        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body; off=0.0)

        mu = vec(body.controlpoints[1, :] .+ 2 .* body.controlpoints[2, :])
        grad_mu = zeros(3, body.ncells)
        te_info = zeros(Int, 2, body.ncells)

        pnl.compute_mu_gradient!(grad_mu, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info; scale=1.0)

        # Panels 1, 2, 4, 5, 7, and 8 have enough in-plane stencil support to
        # recover the exact constant gradient on this small mesh.
        exact_panels = (1, 2, 4, 5, 7, 8)
        expected = [-1.0, -2.0, 0.0]

        for i in exact_panels
            @test isapprox(grad_mu[:, i], expected; atol=1e-9)
        end

        for i in 1:body.ncells
            @test abs(dot(grad_mu[:, i], body.normals[:, i])) ≤ 1e-10
        end

        grad_half = zeros(3, body.ncells)
        pnl.compute_mu_gradient!(grad_half, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info; scale=0.5)

        for i in exact_panels
            @test isapprox(grad_half[:, i], 0.5 .* grad_mu[:, i]; atol=1e-10)
        end
    end

    @testset "compute_mu_gradient! trailing-edge stencil isolation" begin
        nodes, cells = make_postprocess_seeded_te_mesh()
        shedding = pnl.calc_shedding_from_seed(nodes, cells, 1, 2;
            bbox=([0.8, -0.1, -0.1], [1.1, 2.1, 0.1]),
            end_node=3)
        body = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.VortexRing}}(
            nodes, cells, [shedding];
            check_mesh=false,
            watertight=false)

        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body; off=0.0)

        upper_side = vec(body.controlpoints[3, :] .> 0)
        te_info = view(body.shedding_full, 1:2, :)

        mu_ref = vec(body.controlpoints[2, :])
        mu_perturbed = copy(mu_ref)
        mu_perturbed[.!upper_side] .+= 100.0

        grad_ref = zeros(3, body.ncells)
        grad_perturbed = zeros(3, body.ncells)
        pnl.compute_mu_gradient!(grad_ref, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu_ref, te_info; scale=1.0)
        pnl.compute_mu_gradient!(grad_perturbed, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu_perturbed, te_info; scale=1.0)

        @test isapprox(grad_ref[:, upper_side], grad_perturbed[:, upper_side]; atol=1e-12)
        @test isapprox(grad_ref[:, .!upper_side], grad_perturbed[:, .!upper_side]; atol=1e-12)

        grad_no_te_ref = zeros(3, body.ncells)
        grad_no_te_perturbed = zeros(3, body.ncells)
        no_te_info = zeros(Int, 2, body.ncells)
        pnl.compute_mu_gradient!(grad_no_te_ref, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu_ref, no_te_info; scale=1.0)
        pnl.compute_mu_gradient!(grad_no_te_perturbed, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu_perturbed, no_te_info; scale=1.0)

        @test maximum(abs.(grad_no_te_ref[:, upper_side] .- grad_no_te_perturbed[:, upper_side])) > 1.0
        @test maximum(abs.(grad_no_te_ref[:, .!upper_side] .- grad_no_te_perturbed[:, .!upper_side])) > 1.0
    end
end
