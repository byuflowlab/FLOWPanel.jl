using Test
import FLOWPanel as pnl
import FastMultipole
using LinearAlgebra: dot
using SparseArrays
using StaticArrays: SMatrix, SVector

if !isdefined(@__MODULE__, :make_octa_source_body)
    include("test_helpers.jl")
end

struct PostprocessVectorPotentialDummy end
FastMultipole.has_vector_potential(::PostprocessVectorPotentialDummy) = true

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

function expected_negative_tangent_velocity(body)
    expected = similar(body.velocity)
    for p in 1:body.ncells
        n = body.normals[:, p]
        u = body.velocity[:, p]
        expected[:, p] .= -(u .- dot(u, n) .* n .+ body.velocity_kinematic[:, p])
    end
    return expected
end

function make_skewed_two_panel_body()
    nodes = Float64[
        0.0 1.0 0.0 1.2;
        0.0 0.0 1.0 1.1;
        0.0 0.0 0.0 0.4;
    ]
    cells = Int[
        1 2;
        2 4;
        3 3;
    ]
    body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells;
        watertight=false,
        ensure_winding=false)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body; off=0.0)
    return body
end

@testset verbose=true "Postprocess" begin
    @testset "PressureBernoulli unsteady monitor-owned phi_dot" begin
        body = make_octa_source_body()
        uinf_old = [1.0, 0.0, 0.0]
        uinf_new = [1.2, -0.1, 0.05]
        dt = 0.2
        w = [0.3, -0.2, 0.1]
        for p in 1:body.ncells
            body.velocity_kinematic[:, p] .= w
            body.velocity[:, p] .= uinf_new .- w
        end

        monitor = pnl.PressureBernoulli(1.0; unsteady=true, backend=pnl.DirectBackend())
        pnl._pressure_bernoulli_ensure_storage!(monitor, (body,))
        for p in 1:body.ncells
            monitor.potential_history[1][p] = -dot(uinf_old, body.controlpoints[:, p])
        end

        phi_dot = pnl._pressure_bernoulli_phi_dot!(monitor, body, 1, (), uinf_new, dt)
        for p in 1:body.ncells
            phi_old = dot(uinf_old, body.controlpoints[:, p])
            phi_new = dot(uinf_new, body.controlpoints[:, p])
            expected = (phi_new - phi_old) / dt - dot(w, uinf_new)
            @test isapprox(phi_dot[p], expected; atol=1e-12)
            @test isapprox(monitor.potential_history[1][p], -phi_new; atol=1e-12)
        end

        scalar_sources = pnl._filter_scalar_potential_sources((body, PostprocessVectorPotentialDummy()))
        @test scalar_sources == (body,)

        steady_from_nothing = zeros(body.ncells)
        steady_from_zero = zeros(body.ncells)
        pnl.calcfield_P!(steady_from_nothing, body, body.velocity, 1.0, 1.0, nothing;
            correct_kuttacondition=false)
        pnl.calcfield_P!(steady_from_zero, body, body.velocity, 1.0, 1.0, zeros(body.ncells);
            correct_kuttacondition=false)
        @test steady_from_nothing == steady_from_zero
    end

    @testset "PressureLaplace monitor" begin
        body = make_octa_source_body()
        monitor = pnl.PressureLaplace((body,), 1.2; reference_panel=1, reference_pressure=0.0)

        @test pnl.monitor_provides(monitor) == (:P,)
        @test pnl.audit_monitors((monitor, pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization()))) !== nothing
        @test_throws ArgumentError pnl.PressureLaplace(1.0)

        body.velocity .= 0.0
        body.velocity[1, :] .= 1.0
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.25)
        @test isapprox(monitor.velocity_dot[1], expected_negative_tangent_velocity(body); atol=1e-12)
        @test all(isfinite, body.P)

        body.velocity[1, :] .+= 0.5
        body.velocity[2, :] .-= 0.25
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 1, 0.25)
        @test isapprox(monitor.velocity_dot[1], expected_negative_tangent_velocity(body); atol=1e-12)
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

    @testset "PressureLaplace co-normal metric" begin
        body = make_skewed_two_panel_body()
        edges = pnl._pressure_panel_edges(body)
        @test size(edges, 2) == 1
        edge_a, edge_b, i, j = edges[:, 1]
        w, ell, nu1, nu2, nu3, n1, n2, n3 =
            pnl._pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)
        r = body.controlpoints[:, j] .- body.controlpoints[:, i]
        d = sqrt(dot(r, r))
        @test w > 0.0
        @test !isapprox(w, ell / d; atol=1e-12)
        @test dot([nu1, nu2, nu3], r) > 0.0
        @test isapprox(dot([nu1, nu2, nu3], [n1, n2, n3]), 0.0; atol=1e-12)

        L0 = pnl._assemble_pressure_laplacian(body, 0)
        @test isapprox(Matrix(L0), [w -w; -w w]; atol=1e-12)
        L = pnl._assemble_pressure_laplacian(body, 1)
        @test isapprox(Matrix(L), [1.0 0.0; 0.0 w]; atol=1e-12)
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

    @testset "KuttaJoukowskiForce frame rotation" begin
        body = make_plate_vortex_body()
        body.strength[:, 1] .= [0.8, -0.35]
        frames = pnl.ReferenceFrame(body;
            ω_axis=SVector{3}(0.0, 0.0, 1.0),
            R=SMatrix{3,3}(0.0, -1.0, 0.0,
                           1.0,  0.0, 0.0,
                           0.0,  0.0, 1.0))
        uinf = zeros(3)

        global_monitor = pnl.KuttaJoukowskiForce(body, 1, 1;
            backend=pnl.DirectBackend(),
            i_frame=-1,
            normalization=pnl.NoNormalization())
        frame_monitor = pnl.KuttaJoukowskiForce(body, 1, 1;
            backend=pnl.DirectBackend(),
            i_frame=1,
            normalization=pnl.NoNormalization())

        global_monitor((body,), (nothing,), frames, uinf, 0, 0.1)
        frame_monitor((body,), (nothing,), frames, uinf, 0, 0.1)

        _, R_f2g = pnl.frame_global_transform(frames, 1)
        expected_frame_force = transpose(R_f2g) * SVector{3}(global_monitor.force[:, 1]...)
        @test isapprox(frame_monitor.force[:, 1], collect(expected_frame_force); atol=1e-10)
    end

    @testset "PressureLaplace Bernoulli constant-field comparison" begin
        body_b = make_octa_source_body()
        body_l = make_octa_source_body()
        body_b.velocity[1, :] .= 0.4
        body_b.velocity[2, :] .= -0.2
        body_l.velocity .= body_b.velocity

        laplace = pnl.PressureLaplace((body_l,), 1.0; reference_panel=1)
        laplace.velocity_dot[1] .= expected_negative_tangent_velocity(body_l)
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
