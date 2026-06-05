using Test
import FLOWPanel as pnl
import FastMultipole
using LinearAlgebra: cross, dot, norm
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
    pnl.calc_controlpoints!(body)
    return body
end

function make_skewed_two_panel_doublet_body()
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
    body = pnl.NonLiftingBody{pnl.ConstantDoublet}(nodes, cells;
        watertight=false,
        ensure_winding=false)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.strength[:, 1] .= [0.2, 1.0]
    return body
end

function make_spanwise_loading_body(ys)
    n = length(ys)
    nodes = zeros(Float64, 3, 3n)
    cells = zeros(Int, 3, n)
    for (i, y) in enumerate(ys)
        j = 3 * (i - 1)
        nodes[:, j + 1] .= [0.0, y, 0.0]
        nodes[:, j + 2] .= [1.0, y, 0.0]
        nodes[:, j + 3] .= [0.0, y, 1.0]
        cells[:, i] .= (j + 1, j + 2, j + 3)
    end
    body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells;
        watertight=false,
        ensure_winding=false)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    return body
end

function make_surface_vorticity_gradient_body()
    nodes, cells = make_planar_gradient_mesh()
    body = pnl.NonLiftingBody{pnl.ConstantDoublet}(nodes, cells;
        watertight=false,
        ensure_winding=false)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.strength[:, 1] .= vec(body.controlpoints[1, :] .+ 2 .* body.controlpoints[2, :])
    return body
end

function make_bound_circulation_side_body()
    nodes = Float64[
         1.0  1.0  2.0  -1.0 -1.0 -2.0;
         0.0  1.0  0.5   0.0  1.0  0.5;
         0.0  0.0  0.0   0.0  0.0  0.0;
    ]
    cells = Int[
        1 4;
        2 5;
        3 6;
    ]
    shedding = [reshape(Int[1, 1, 2, -1, -1, -1], 6, 1)]
    body = pnl.RigidWakeBody{pnl.VortexRing}(nodes, cells, shedding;
        check_mesh=false,
        watertight=false,
        ensure_winding=false)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.strength[:, 1] .= [2.0, 100.0]
    return body
end

function rotated_bound_body(body, R)
    nodes = Matrix(R * body.nodes)
    new_body = pnl.RigidWakeBody{pnl.VortexRing}(nodes, copy(body.cells), deepcopy(body.shedding);
        check_mesh=false,
        watertight=false,
        ensure_winding=false)
    pnl.calc_normals!(new_body)
    pnl.calc_controlpoints!(new_body)
    new_body.strength .= body.strength
    return new_body
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
        @test !monitor.unsteady
        @test monitor.acceleration_form == :material_derivative
        @test pnl.monitor_requires_body_hessian(monitor)
        @test pnl.audit_monitors((monitor, pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization()))) !== nothing
        @test_throws ArgumentError pnl.PressureLaplace(1.0)
        @test_throws ArgumentError pnl.PressureLaplace((body,), 1.0; gradient_mode=:unknown)
        @test_throws ArgumentError pnl.PressureLaplace((body,), 1.0; acceleration_form=:unknown)

        lamb = pnl.PressureLaplace((body,), 1.2; reference_panel=1,
            acceleration_form=:lamb_vector)
        @test lamb.acceleration_form == :lamb_vector
        @test !pnl.monitor_requires_body_hessian(lamb)
        @test pnl.monitor_requires_induced_vorticity(lamb)

        body.velocity .= 0.0
        body.velocity[1, :] .= 1.0
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.25)
        @test isapprox(monitor.velocity_dot[1], expected_negative_tangent_velocity(body); atol=1e-12)
        @test all(isfinite, monitor.p[1])

        body.velocity[1, :] .+= 0.5
        body.velocity[2, :] .-= 0.25
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 1, 0.25)
        @test isapprox(monitor.velocity_dot[1], expected_negative_tangent_velocity(body); atol=1e-12)
    end

    @testset "PressureLaplace surface velocity gradient mode" begin
        body = make_octa_source_body()
        monitor = pnl.PressureLaplace((body,), 1.0;
            reference_panel=1,
            gradient_mode=:surface_velocity)

        @test !pnl.monitor_requires_body_hessian(monitor)

        body.velocity .= 0.0
        body.velocity[1, :] .= 0.3
        body.velocity[2, :] .= -0.1
        body.velocity_gradient .= NaN
        monitor.velocity_dot[1] .= expected_negative_tangent_velocity(body)

        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.25)

        @test all(isfinite, monitor.p[1])
        @test all(isfinite, monitor.acceleration[1])
        @test all(isfinite, monitor.surface_velocity_gradient[1])

        lamb = pnl.PressureLaplace((body,), 1.0;
            reference_panel=1,
            gradient_mode=:surface_velocity,
            acceleration_form=:lamb_vector)
        body.velocity_gradient .= NaN
        lamb.velocity_dot[1] .= expected_negative_tangent_velocity(body)
        lamb((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.25)
        @test all(isfinite, lamb.p[1])
        @test all(isfinite, lamb.surface_velocity_gradient[1])
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

    @testset "PressureLaplace rebuild policy and force integration" begin
        body = make_octa_source_body()
        pressure = pnl.PressureLaplace((body,), 1.0; reference_panel=1)
        force = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())
        frames = pnl.ReferenceFrame(body)

        body.velocity[1, :] .= 0.2 .* (1:body.ncells)
        pressure((body,), (nothing,), frames, zeros(3), 0, 0.1)
        first_L = pressure.L[1]
        pressure((body,), (nothing,), frames, zeros(3), 1, 0.1)
        @test pressure.L[1] === first_L

        body.nodes[1, :] .+= 0.2
        body.nodes[2, :] .-= 0.1
        pressure((body,), (nothing,), frames, zeros(3), 2, 0.1)
        @test pressure.L[1] === first_L

        old = body.nodes[1, 1]
        body.nodes[1, 1] = old + 0.1
        pressure((body,), (nothing,), frames, zeros(3), 3, 0.1)
        @test pressure.L[1] === first_L

        rebuilding = pnl.PressureLaplace((body,), 1.0;
            reference_panel=1,
            rebuild_every_step=true)
        body.velocity[1, :] .= 0.2 .* (1:body.ncells)
        rebuilding((body,), (nothing,), frames, zeros(3), 4, 0.1)
        rebuild_first_L = rebuilding.L[1]
        rebuilding((body,), (nothing,), frames, zeros(3), 5, 0.1)
        @test rebuilding.L[1] !== rebuild_first_L

        ctx = pnl.MonitorContext()
        pnl.monitor_register!(ctx, :P, 1, rebuilding.p[1])
        pnl._run_monitor!(force, ctx, (body,), (nothing,), frames, zeros(3), 0, 0.1)
        @test all(isfinite, rebuilding.p[1])
        @test all(isfinite, force.distributed_force)
        @test all(isfinite, force.force)
    end

    @testset "SpanwiseLoadingMonitor binning and validation" begin
        body = make_spanwise_loading_body([0.125, 0.375, 0.625, 0.875])
        force = Float64[
            10.0 20.0 30.0 40.0;
             0.0  0.0  0.0  0.0;
             1.0  2.0  3.0  4.0;
        ]
        ctx = pnl.MonitorContext()
        pnl.monitor_register!(ctx, :F, 1, force)
        monitor = pnl.SpanwiseLoadingMonitor(2, 1;
            components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]))

        @test pnl.monitor_requires(monitor) == (:F,)
        @test_throws ArgumentError pnl.audit_monitors((monitor,))
        @test pnl.audit_monitors((pnl.SurfaceVorticityForce(make_plate_vortex_body(), 1, 1), monitor)) !== nothing
        @test isapprox(monitor.span_axis, SVector(0.0, 1.0, 0.0); atol=1e-12)
        @test_throws ArgumentError pnl.SpanwiseLoadingMonitor(2, 1;
            components=(lift=[0.0, 0.0, 2.0], drag=[1.0, 0.0, 0.0]))
        @test_throws ArgumentError pnl.SpanwiseLoadingMonitor(2, 1;
            components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]),
            span_axis=[0.0, 0.0, 1.0])

        pnl._run_monitor!(monitor, ctx, (body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.1)
        @test monitor.counts == [2, 2]
        @test monitor.panel_bin_id == [1, 1, 2, 2]
        @test isapprox(monitor.bin_center, [0.3125, 0.6875]; atol=1e-12)
        @test isapprox(monitor.load_components, [3.0 7.0; 30.0 70.0]; atol=1e-12)
        @test isapprox(monitor.force_components, monitor.load_components; atol=1e-12)

        explicit = pnl.SpanwiseLoadingMonitor(2, 1;
            components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]),
            span_axis=[0.0, 1.0, 0.0])
        @test isapprox(explicit.span_axis, monitor.span_axis; atol=1e-12)
    end

    @testset "SpanwiseLoadingMonitor frame, interpolation, normalization, CSV, VTK" begin
        R = SMatrix{3,3}(0.0, -1.0, 0.0,
                         1.0,  0.0, 0.0,
                         0.0,  0.0, 1.0)
        origin = SVector(2.0, -1.0, 0.5)
        local_body = make_spanwise_loading_body([0.0, 1.0])
        nodes_global = similar(local_body.nodes)
        for i in axes(local_body.nodes, 2)
            nodes_global[:, i] .= origin + R * SVector{3}(local_body.nodes[:, i])
        end
        body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes_global, copy(local_body.cells);
            watertight=false,
            ensure_winding=false)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        frames = pnl.ReferenceFrame(body; origin=origin, R=R)

        f_local = [0.0 0.0; 0.0 0.0; 2.0 8.0]
        f_global = Matrix(R * f_local)
        ctx = pnl.MonitorContext()
        pnl.monitor_register!(ctx, :F, 1, f_global)
        monitor = pnl.SpanwiseLoadingMonitor(3, 1; i_frame=1,
            components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]),
            per_length=true)
        pnl._run_monitor!(monitor, ctx, (body,), (nothing,), frames, zeros(3), 0, 0.1)
        @test monitor.counts == [1, 0, 1]
        @test monitor.panel_bin_id == [1, 3]
        @test isapprox(monitor.force_components[1, :], [2.0, 5.0, 8.0]; atol=1e-12)
        @test isapprox(monitor.load_components[1, :], [6.0, 15.0, 24.0]; atol=1e-12)

        fs = pnl.SpanwiseLoadingMonitor(1, 1;
            components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]),
            normalization=pnl.FreestreamSectionalNormalization(1.0, 2.0),
            per_length=true)
        ctx = pnl.MonitorContext()
        pnl.monitor_register!(ctx, :F, 1, [0.0 0.0; 0.0 0.0; 4.0 4.0])
        pnl._run_monitor!(fs, ctx, (local_body,), (nothing,), pnl.ReferenceFrame(local_body), [2.0, 0.0, 0.0], 0, 0.1)
        @test isapprox(fs.load_components[1, 1], 2.0; atol=1e-12)

        rotor = pnl.SpanwiseLoadingMonitor(1, 1;
            components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]),
            normalization=pnl.RotorSectionalNormalization(2.0, 2.0, 1))
        rframes = pnl.ReferenceFrame(local_body; ω=3.0)
        ctx = pnl.MonitorContext()
        pnl.monitor_register!(ctx, :F, 1, [0.0 0.0; 0.0 0.0; 36.0 36.0])
        pnl._run_monitor!(rotor, ctx, (local_body,), (nothing,), rframes, zeros(3), 0, 0.1)
        @test isapprox(rotor.load_components[1, 1], 2.0; atol=1e-12)

        csv = joinpath(mktempdir(), "span", "loading.csv")
        csvmon = pnl.SpanwiseLoadingMonitor(1, 1;
            components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]),
            csv_path=csv)
        ctx = pnl.MonitorContext()
        pnl.monitor_register!(ctx, :F, 1, [0.0 0.0; 0.0 0.0; 1.0 1.0])
        pnl._run_monitor!(csvmon, ctx, (local_body,), (nothing,), pnl.ReferenceFrame(local_body), zeros(3), 0, 0.1, 1.25)
        pnl._run_monitor!(csvmon, ctx, (local_body,), (nothing,), pnl.ReferenceFrame(local_body), zeros(3), 1, 0.1, 1.50)
        rows = readlines(csv)
        @test rows[1] == "step,time,bin,bin_center,bin_width,count,lift,drag"
        @test length(rows) == 3
        @test startswith(rows[2], "0,1.25,1,")
        @test startswith(rows[3], "1,1.5,1,")

        vtkdir = mktempdir()
        pnl.write_vtk(joinpath(vtkdir, "span_body"), local_body, 0, 0.0;
            monitors=(csvmon,), i_system=1, overwrite=true)
        vtk = pnl.ReadVTK.VTKFile(joinpath(vtkdir, "span_body", "span_body.0.vtu"))
        cell_data = pnl.ReadVTK.get_cell_data(vtk)
        @test "spanwise bin id" in keys(cell_data)
        @test vec(pnl.ReadVTK.get_data(cell_data["spanwise bin id"])) == [1, 1]
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

    @testset "KuttaJoukowskiForce relative translational velocity" begin
        body = make_plate_vortex_body()
        body.strength[:, 1] .= [0.8, -0.35]
        U = [1.2, -0.4, 0.25]
        W = [0.7, 0.3, -0.15]

        fixed_monitor = pnl.KuttaJoukowskiForce(body, 1, 1;
            backend=pnl.DirectBackend(),
            normalization=pnl.NoNormalization())
        moving_monitor = pnl.KuttaJoukowskiForce(body, 1, 1;
            backend=pnl.DirectBackend(),
            normalization=pnl.NoNormalization())

        fixed_monitor((body,), (nothing,), pnl.ReferenceFrame(body), U, 0, 0.1)
        moving_frames = pnl.ReferenceFrame(body; v=SVector{3}(W...))
        moving_monitor((body,), (nothing,), moving_frames, U .+ W, 0, 0.1)

        @test isapprox(moving_monitor.force[:, 1], fixed_monitor.force[:, 1]; atol=1e-10)
    end

    @testset "KuttaJoukowskiForce rotational edge kinematics" begin
        body = make_plate_vortex_body()
        origin = SVector{3}(0.25, -0.1, 0.05)
        v = SVector{3}(0.3, -0.2, 0.4)
        ω_axis = SVector{3}(0.0, 0.0, 1.0)
        ω = 2.5
        frames = pnl.ReferenceFrame(body;
            origin,
            v,
            ω_axis,
            ω)

        monitor = pnl.KuttaJoukowskiForce(body, 1, 1;
            backend=pnl.DirectBackend(),
            normalization=pnl.NoNormalization())
        monitor((body,), (nothing,), frames, zeros(3), 0, 0.1)

        k = 1
        a = monitor.edge_node_a[k]
        b = monitor.edge_node_b[k]
        midpoint = SVector{3}(
            0.5 * (body.nodes[1, a] + body.nodes[1, b]),
            0.5 * (body.nodes[2, a] + body.nodes[2, b]),
            0.5 * (body.nodes[3, a] + body.nodes[3, b]),
        )
        expected = v + cross(ω_axis * ω, midpoint - origin)

        @test isapprox(monitor.velocity_kinematic[:, k], collect(expected); atol=1e-12)
    end

    @testset "KuttaJoukowskiForce preserves kernel offsets" begin
        body = make_plate_vortex_body()
        other = make_octa_source_body()
        body.kerneloffset_panel = 1e-9
        body.kerneloffset_targets = 2e-3
        body.kerneloffset = 7e-4
        other.kerneloffset_panel = 3e-9
        other.kerneloffset_targets = 4e-3
        other.kerneloffset = 9e-4
        frames = pnl.ReferenceFrame(body)

        monitor = pnl.KuttaJoukowskiForce(body, 1, 1;
            backend=pnl.DirectBackend(),
            normalization=pnl.NoNormalization())
        monitor((body, other), (nothing, nothing), frames, zeros(3), 0, 0.1)

        @test body.kerneloffset == 7e-4
        @test other.kerneloffset == 9e-4
    end

    @testset "BoundCirculationMonitor constructor and TE jump" begin
        @test_throws ArgumentError pnl.BoundCirculationMonitor(
            make_octa_source_body(), 1, 1; i_frame=1, radial_dimension=2, R=1.0)

        body = make_plate_vortex_body()
        body.strength[:, 1] .= [0.3, 1.1]
        monitor = pnl.BoundCirculationMonitor(body, 1, 1;
            i_frame=1,
            radial_dimension=1,
            R=1.0,
            section_tol=0.6)

        @test pnl.monitor_requires(monitor) == ()
        @test pnl.monitor_provides(monitor) == ()
        @test size(monitor.r_over_R) == (1, 1)
        @test monitor.valid_section == trues(1, 1)

        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.1)

        @test isapprox(monitor.r_over_R[1, 1], 0.5; atol=1e-12)
        @test isapprox(monitor.circulation_te[1, 1, 1], 0.8; atol=1e-12)
    end

    @testset "BoundCirculationMonitor rotor-frame slicing" begin
        body = make_bound_circulation_side_body()
        monitor = pnl.BoundCirculationMonitor(body, 1, 1;
            i_frame=1,
            radial_dimension=2,
            R=1.0,
            section_tol=0.1)
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.1)

        @test isapprox(monitor.r_over_R[1, 1], 0.5; atol=1e-12)
        @test isapprox(monitor.circulation_te[1, 1, 1], 2.0; atol=1e-12)
        @test isapprox(monitor.circulation_slice[1, 1, 1], 2.0; atol=1e-12)

        narrow = pnl.BoundCirculationMonitor(body, 1, 1;
            i_frame=1,
            radial_dimension=2,
            R=1.0,
            section_tol=0.01)
        body.controlpoints[2, 1] = 0.55
        narrow((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.1)
        @test isapprox(narrow.circulation_slice[1, 1, 1], 0.0; atol=1e-12)
    end

    @testset "BoundCirculationMonitor rotor-frame invariance" begin
        body = make_bound_circulation_side_body()
        Rz = SMatrix{3,3}(0.0, -1.0, 0.0,
                          1.0,  0.0, 0.0,
                          0.0,  0.0, 1.0)
        rotated = rotated_bound_body(body, Rz)

        base_monitor = pnl.BoundCirculationMonitor(body, 1, 1;
            i_frame=1,
            radial_dimension=2,
            R=1.0,
            section_tol=0.1)
        rotated_monitor = pnl.BoundCirculationMonitor(rotated, 1, 1;
            i_frame=1,
            radial_dimension=2,
            R=1.0,
            section_tol=0.1)

        base_monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.1)
        frames = pnl.ReferenceFrame(rotated; R=Rz)
        rotated_monitor((rotated,), (nothing,), frames, zeros(3), 0, 0.1)

        @test isapprox(rotated_monitor.r_over_R, base_monitor.r_over_R; atol=1e-12)
        @test isapprox(rotated_monitor.circulation_te, base_monitor.circulation_te; atol=1e-12)
        @test isapprox(rotated_monitor.circulation_slice, base_monitor.circulation_slice; atol=1e-12)
    end

    @testset "SurfaceVorticityForce analytic panel force and integration" begin
        body = make_surface_vorticity_gradient_body()
        body.velocity[1, :] .= 3.0
        body.velocity[2, :] .= -0.5
        body.velocity[3, :] .= 0.25
        rho = 1.7
        monitor = pnl.SurfaceVorticityForce(body, 1, 1;
            rho,
            normalization=pnl.NoNormalization())

        @test pnl.monitor_provides(monitor) == (:F,)
        @test pnl.monitor_requires(monitor) == ()

        frames = pnl.ReferenceFrame(body)
        monitor((body,), (nothing,), frames, zeros(3), 0, 0.1)

        exact_panels = (1, 2, 4, 5, 7, 8)
        grad_mu = SVector{3}(1.0, 2.0, 0.0)
        V = SVector{3}(3.0, -0.5, 0.25)
        areas = pnl.calc_areas(body)
        for i in exact_panels
            n = SVector{3}(body.normals[:, i]...)
            kappa = -cross(n, grad_mu)
            expected = rho .* cross(V, kappa) .* areas[i]
            @test isapprox(monitor.distributed_force[:, i], collect(expected); atol=1e-9)
        end

        expected_force = pnl.calcfield_Ftot(body, monitor.distributed_force)
        expected_moment = pnl.calcfield_Mtot(body, zeros(3), body.controlpoints, monitor.distributed_force)
        @test isapprox(monitor.force[:, 1], expected_force; atol=1e-12)
        @test isapprox(monitor.moment[:, 1], expected_moment; atol=1e-12)
    end

    @testset "SurfaceVorticityForce frame rotation and zero velocity" begin
        body = make_surface_vorticity_gradient_body()
        body.velocity[1, :] .= 1.0
        body.velocity[3, :] .= -0.2
        frames = pnl.ReferenceFrame(body;
            R=SMatrix{3,3}(0.0, -1.0, 0.0,
                           1.0,  0.0, 0.0,
                           0.0,  0.0, 1.0))

        global_monitor = pnl.SurfaceVorticityForce(body, 1, 1;
            normalization=pnl.NoNormalization())
        frame_monitor = pnl.SurfaceVorticityForce(body, 1, 1;
            i_frame=1,
            normalization=pnl.NoNormalization())
        global_monitor((body,), (nothing,), frames, zeros(3), 0, 0.1)
        frame_monitor((body,), (nothing,), frames, zeros(3), 0, 0.1)

        _, R_f2g = pnl.frame_global_transform(frames, 1)
        expected_frame_force = transpose(R_f2g) * SVector{3}(global_monitor.force[:, 1]...)
        expected_frame_moment = transpose(R_f2g) * SVector{3}(global_monitor.moment[:, 1]...)
        @test isapprox(frame_monitor.force[:, 1], collect(expected_frame_force); atol=1e-12)
        @test isapprox(frame_monitor.moment[:, 1], collect(expected_frame_moment); atol=1e-12)

        zero_monitor = pnl.SurfaceVorticityForce(body, 1, 1;
            normalization=pnl.NoNormalization())
        body.velocity .= 0.0
        zero_monitor((body,), (nothing,), frames, zeros(3), 0, 0.1)
        @test isapprox(zero_monitor.distributed_force, zeros(size(zero_monitor.distributed_force)); atol=1e-12)
        @test isapprox(zero_monitor.force[:, 1], zeros(3); atol=1e-12)
    end

    @testset "Bound surface vorticity accumulation" begin
        body = make_surface_vorticity_gradient_body()
        seed = SVector{3}(0.3, -0.7, 1.1)
        for i in axes(body.induced_vorticity, 2)
            body.induced_vorticity[:, i] .= seed
        end

        pnl._add_bound_surface_vorticity!(body)

        exact_panels = (1, 2, 4, 5, 7, 8)
        grad_mu = SVector{3}(1.0, 2.0, 0.0)
        for i in exact_panels
            n = SVector{3}(body.normals[:, i]...)
            expected = collect(seed) .+ collect(cross(n, grad_mu))
            @test isapprox(body.induced_vorticity[:, i], expected; atol=1e-9)
        end
    end

    @testset "SurfaceVorticityForce constructor validation" begin
        body = make_octa_source_body()
        @test_throws ArgumentError pnl.SurfaceVorticityForce(body, 1, 1)
    end

    @testset "PressureLaplace Bernoulli constant-field comparison" begin
        body_b = make_octa_source_body()
        body_l = make_octa_source_body()
        body_b.velocity[1, :] .= 0.4
        body_b.velocity[2, :] .= -0.2
        body_l.velocity .= body_b.velocity

        laplace = pnl.PressureLaplace((body_l,), 1.0; reference_panel=1)
        laplace.velocity_dot[1] .= expected_negative_tangent_velocity(body_l)
        p_b_raw = zeros(body_b.ncells)
        pnl.calcfield_P!(p_b_raw, body_b, body_b.velocity, 1.0, 1.0, zeros(body_b.ncells);
            correct_kuttacondition=false)
        laplace((body_l,), (nothing,), pnl.ReferenceFrame(body_l), [1.0, 0.0, 0.0], 0, 0.1)

        p_b = p_b_raw .- p_b_raw[1]
        p_l = laplace.p[1] .- laplace.p[1][1]
        @test isapprox(p_l, p_b; atol=1e-10)
    end

    @testset "PressureLaplace acceleration RHS projection" begin
        body = make_skewed_two_panel_body()
        laplace = pnl.PressureLaplace((body,), 1.2; reference_panel=1)
        acceleration = zeros(3, body.ncells)
        p_exact = [0.0, 2.0]
        edge_a, edge_b, i, j = laplace.edges[1][:, 1]
        r = body.controlpoints[:, j] .- body.controlpoints[:, i]
        acceleration[:, i] .= ((p_exact[i] - p_exact[j]) / laplace.rho) .* r ./ dot(r, r)
        acceleration[:, j] .= acceleration[:, i]

        pnl._pressure_rhs_from_acceleration!(laplace.b[1], laplace, body, 1, acceleration)

        @test isapprox(laplace.b[1], laplace.L[1] * p_exact; atol=1e-12)
    end

    @testset "PressureLaplace edge material derivative RHS projection" begin
        body = make_skewed_two_panel_body()
        laplace = pnl.PressureLaplace((body,), 1.2; reference_panel=1)
        velocity_dot = zeros(3, body.ncells)
        p_exact = [0.0, 2.0]
        edge_a, edge_b, i, j = laplace.edges[1][:, 1]
        r = body.controlpoints[:, j] .- body.controlpoints[:, i]
        velocity_dot[:, i] .= ((p_exact[i] - p_exact[j]) / laplace.rho) .* r ./ dot(r, r)
        velocity_dot[:, j] .= velocity_dot[:, i]

        pnl._pressure_rhs_from_edge_material_derivative!(laplace.b[1], laplace, body, 1, nothing)
        @test isapprox(laplace.b[1], zeros(2); atol=1e-12)

        pnl._pressure_rhs_from_edge_material_derivative!(laplace.b[1], laplace, body, 1, velocity_dot)

        @test isapprox(laplace.b[1], laplace.L[1] * p_exact; atol=1e-12)
    end

    @testset "PressureLaplace Lamb vector RHS projection" begin
        body = make_skewed_two_panel_body()
        laplace = pnl.PressureLaplace((body,), 1.2; reference_panel=1,
            acceleration_form=:lamb_vector)
        p_exact = [0.0, -2.0]
        q = -p_exact ./ laplace.rho
        body.velocity .= 0.0
        edge_a, edge_b, i, j = laplace.edges[1][:, 1]
        r = body.controlpoints[:, j] .- body.controlpoints[:, i]
        for p in 1:body.ncells
            n = body.normals[:, p]
            t = r .- dot(r, n) .* n
            t ./= norm(t)
            body.velocity[:, p] .= sqrt(2.0 * q[p]) .* t
        end
        body.velocity_gradient .= 0.0
        body.angular_velocity .= 0.0

        pnl._pressure_rhs_from_lamb_vector!(laplace.b[1], laplace, body, 1, nothing)

        @test isapprox(laplace.b[1], laplace.L[1] * p_exact; atol=1e-12)

        body.velocity .= 0.0
        pnl._pressure_rhs_from_lamb_vector!(laplace.b[1], laplace, body, 1, nothing)
        @test isapprox(laplace.b[1], zeros(2); atol=1e-12)

        laplace.u_inertial[1] .= 0.0
        body.velocity .= 0.0
        for p in 1:body.ncells
            n = body.normals[:, p]
            t = r .- dot(r, n) .* n
            t ./= norm(t)
            body.velocity[:, p] .= t
        end
        body.velocity_gradient .= NaN
        body.induced_vorticity .= 0.0
        pnl._pressure_rhs_from_lamb_vector!(laplace.b[1], laplace, body, 1, nothing)
        @test isapprox(laplace.b[1], zeros(2); atol=1e-12)

        body.induced_vorticity[3, :] .= 1.0
        pnl._pressure_rhs_from_lamb_vector!(laplace.b[1], laplace, body, 1, nothing)
        @test !isapprox(laplace.b[1], zeros(2); atol=1e-12)

        doublet_body = make_skewed_two_panel_doublet_body()
        doublet_laplace = pnl.PressureLaplace((doublet_body,), 1.2;
            reference_panel=1,
            acceleration_form=:lamb_vector)
        doublet_laplace.u_inertial[1] .= 0.0
        doublet_body.velocity .= 0.0
        doublet_body.velocity[1, :] .= 0.7
        doublet_body.velocity[2, :] .= -0.2
        pnl._add_bound_surface_vorticity!(doublet_body)
        pnl._pressure_rhs_from_lamb_vector!(doublet_laplace.b[1], doublet_laplace,
            doublet_body, 1, nothing)
        surface_only_rhs = copy(doublet_laplace.b[1])
        @test !isapprox(surface_only_rhs, zeros(2); atol=1e-12)

        doublet_body.induced_vorticity[3, :] .+= 0.5
        pnl._pressure_rhs_from_lamb_vector!(doublet_laplace.b[1], doublet_laplace,
            doublet_body, 1, nothing)
        @test !isapprox(doublet_laplace.b[1], surface_only_rhs; atol=1e-12)
    end

    @testset "PressureLaplace unsteady toggle" begin
        body_steady = make_skewed_two_panel_body()
        body_unsteady = make_skewed_two_panel_body()
        for p in 1:body_steady.ncells
            body_steady.velocity[:, p] .= [0.2 + 0.1p, -0.05p, 0.03]
        end
        body_unsteady.velocity .= body_steady.velocity

        steady = pnl.PressureLaplace((body_steady,), 1.0; reference_panel=1)
        unsteady = pnl.PressureLaplace((body_unsteady,), 1.0; reference_panel=1,
            unsteady=true)
        steady.velocity_dot[1] .= 0.0
        unsteady.velocity_dot[1] .= 0.0

        steady((body_steady,), (nothing,), pnl.ReferenceFrame(body_steady), zeros(3), 0, 0.25)
        unsteady((body_unsteady,), (nothing,), pnl.ReferenceFrame(body_unsteady), zeros(3), 0, 0.25)

        @test isapprox(steady.velocity_dot[1], expected_negative_tangent_velocity(body_steady); atol=1e-12)
        @test isapprox(unsteady.velocity_dot[1], expected_negative_tangent_velocity(body_unsteady); atol=1e-12)
        @test !isapprox(steady.b[1], unsteady.b[1]; atol=1e-12)
    end

    @testset "compute_mu_gradient! interior recovery" begin
        nodes, cells = make_planar_gradient_mesh()
        body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells;
            watertight=false,
            ensure_winding=false)

        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)

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

    @testset "compute_mu_gradient! quad-pair consistency" begin
        nodes, cells = make_planar_gradient_mesh()
        body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells;
            watertight=false,
            ensure_winding=false)

        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)

        mu = [0.0, 2.0, -1.0, 5.0, 1.5, -3.0, 4.0, -2.5]
        te_info = zeros(Int, 2, body.ncells)
        grad_raw = zeros(3, body.ncells)
        grad_quad = zeros(3, body.ncells)

        pnl.compute_mu_gradient!(grad_raw, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info; scale=1.0)
        pnl.compute_mu_gradient!(grad_quad, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info; scale=1.0, nodes=body.nodes)

        pairs = ((1, 2), (3, 4), (5, 6), (7, 8))
        @test maximum(maximum(abs.(grad_raw[:, i] .- grad_raw[:, j])) for (i, j) in pairs) > 1e-3
        for (i, j) in pairs
            @test isapprox(grad_quad[:, i], grad_quad[:, j]; atol=1e-12)
        end
        for i in 1:body.ncells
            @test abs(dot(grad_quad[:, i], body.normals[:, i])) ≤ 1e-10
        end
    end

    @testset "compute_surface_velocity_gradient! interior recovery" begin
        nodes, cells = make_planar_gradient_mesh()
        body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells;
            watertight=false,
            ensure_winding=false)

        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)

        u = zeros(3, body.ncells)
        u[1, :] .= vec(body.controlpoints[1, :] .+ 2 .* body.controlpoints[2, :])
        u[2, :] .= vec(-3 .* body.controlpoints[1, :] .+ 0.5 .* body.controlpoints[2, :])
        u[3, :] .= 2.0
        grad_u = zeros(3, 3, body.ncells)
        te_info = zeros(Int, 2, body.ncells)

        pnl.compute_surface_velocity_gradient!(grad_u, u, body.controlpoints,
            body.normals, body.cells, body.neighbor, te_info)

        exact_panels = (1, 2, 4, 5, 7, 8)
        for i in exact_panels
            @test isapprox(grad_u[1, :, i], [1.0, 2.0, 0.0]; atol=1e-9)
            @test isapprox(grad_u[2, :, i], [-3.0, 0.5, 0.0]; atol=1e-9)
            @test isapprox(grad_u[3, :, i], [0.0, 0.0, 0.0]; atol=1e-9)
        end

        for i in 1:body.ncells, k in 1:3
            @test abs(dot(grad_u[k, :, i], body.normals[:, i])) ≤ 1e-10
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
        pnl.calc_controlpoints!(body)

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
