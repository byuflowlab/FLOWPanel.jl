using Test
import FLOWPanel as pnl
import FastMultipole
using LinearAlgebra: cross, dot, mul!, norm
using SparseArrays
using StaticArrays: SMatrix, SVector

if !isdefined(@__MODULE__, :compare_pressure_models)
    include(joinpath(@__DIR__, "..", "examples",
        "simple_wing_capped_pressure_comparison.jl"))
end

if !isdefined(@__MODULE__, :make_octa_source_body)
    include("test_helpers.jl")
end

struct PostprocessVectorPotentialDummy end
FastMultipole.has_vector_potential(::PostprocessVectorPotentialDummy) = true

struct PostprocessVectorWake <: pnl.AbstractFreeWake
    source::PostprocessVectorPotentialDummy
end
pnl.get_sources(w::PostprocessVectorWake) = (w.source,)
pnl._scalar_potential_sources(::PostprocessVectorWake) = ()

function make_bernoulli_limit_body(kind::Symbol)
    if kind == :doublet
        body = pnl.NonLiftingBody{pnl.ConstantDoublet}(copy(NODES_2TRI), copy(CELLS_2TRI);
            watertight=false, ensure_winding=false)
    elseif kind == :ring
        body = make_plate_vortex_body()
    elseif kind == :source_doublet
        body = pnl.NonLiftingBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}(
            copy(NODES_2TRI), copy(CELLS_2TRI);
            DBC=true, watertight=false, ensure_winding=false)
    elseif kind == :source_ring
        base = make_plate_vortex_body()
        body = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.VortexRing}}(
            copy(base.nodes), copy(base.cells), deepcopy(base.shedding);
            DBC=true, check_mesh=false, watertight=false, ensure_winding=false)
        for Da in body.Das
            Da .= repeat([1.0, 0.0, 0.0], 1, size(Da, 2))
        end
    else
        error("unknown Bernoulli limit body kind $(kind)")
    end
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.strength .= 0.0
    gammai = pnl.get_Gammai(body)
    body.strength[:, gammai] .= range(0.2, 0.8; length=body.ncells)
    return body
end

function raw_controlpoint_potential(body, backend)
    probes = FastMultipole.ProbeSystem(body.ncells, Float64)
    for p in 1:body.ncells
        probes.position[p] = SVector{3, Float64}(body.controlpoints[:, p])
        probes.scalar_potential[p] = 0.0
    end
    pnl.influence!((probes,), (body,), backend;
        precalc=false, scalar_potential=true, gradient=false, hessian=(false,))
    return copy(probes.scalar_potential)
end

function make_active_final_filament_wake(body)
    wake = pnl.PanelWake(body; nwakerows=2, include_final_filament=true)
    wake.nwakes[] = 1
    wake.overflowed[] = true
    for (nodes, strength) in zip(wake.nodes, wake.strength)
        nodes[:, 2, 1] .= (-0.5, -0.7, 0.6)
        nodes[:, 2, 2] .= (1.5, -0.7, 0.6)
        strength[1, 2, 1] = 0.9
    end
    return wake
end

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

function make_spanwise_overlap_body()
    nodes = Float64[
        0.0 1.0 0.0;
        0.25 0.75 0.25;
        0.0 0.0 1.0;
    ]
    cells = reshape(Int[1, 2, 3], 3, 1)
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
    @testset "Procedural capped-wing pressure comparison" begin
        chord, span = 1.0, 4.0
        nodes, cells = pressure_comparison_wing_mesh(chord, span;
            thickness=0.15, n_span=2, n_airfoil=21, n_endcap=3)
        @test pnl.iswatertight(nodes, cells) == (true, Int[])
        area2 = [norm(cross(
            nodes[:, cells[2, i]] - nodes[:, cells[1, i]],
            nodes[:, cells[3, i]] - nodes[:, cells[1, i]]))
            for i in axes(cells, 2)]
        @test minimum(area2) > 100 * eps(Float64) * chord^2

        body = build_pressure_comparison_wing(; chord, span,
            n_span=2, n_airfoil=21, n_endcap=3)
        @test body.ncells == 152
        @test size(body.shedding[1]) == (6, 2)
        te_y = Float64[]
        for col in axes(body.shedding[1], 2)
            panel = body.shedding[1][1, col]
            edge_nodes = body.cells[body.shedding[1][2:3, col], panel]
            @test all(body.nodes[1, edge_nodes] .≈ chord)
            append!(te_y, body.nodes[2, edge_nodes])
        end
        @test collect(extrema(te_y)) ≈ [-span / 2, span / 2]

        levels = (
            (name="coarse", n_span=2, n_airfoil=21, n_endcap=3),
            (name="fine", n_span=4, n_airfoil=41, n_endcap=5),
        )
        coarse, fine = run_pressure_refinement(; levels)
        @test fine.pressure_rel_error <= max(coarse.pressure_rel_error, 1e-8)
        @test fine.force_rel_error <= max(coarse.force_rel_error, 1e-7)
        @test fine.pressure_rel_error < 3e-2
        @test fine.force_rel_error < 2e-3
        @test fine.gradient_mode == :edge_difference
        @test fine.acceleration_form == :material_derivative

        default_levels = PRESSURE_COMPARISON_LEVELS
        @test [(l.n_span, l.n_airfoil, l.n_endcap) for l in default_levels] ==
              [(2, 21, 3), (4, 41, 5), (8, 81, 7)]
        @test [build_pressure_comparison_wing(; n_span=l.n_span,
                    n_airfoil=l.n_airfoil, n_endcap=l.n_endcap).ncells
               for l in default_levels] == [152, 624, 2216]
        overridden = compare_pressure_models(; chord=2.0, aspect_ratio=9.0,
            span=6.0, n_span=2, n_airfoil=21, n_endcap=3)
        @test overridden.aspect_ratio == 3.0
        @test overridden.reference_area == 12.0
    end

    @testset "PressureBernoulli unsteady monitor-owned phi_dot" begin
        body = make_octa_source_body()
        w = [0.3, -0.2, 0.1]
        monitor = pnl.PressureBernoulli(1.0; unsteady=true, backend=pnl.DirectBackend())
        pnl._pressure_bernoulli_ensure_storage!(monitor, (body,))
        us = ([1.0, 0.0, 0.0], [1.2, -0.1, 0.05], [1.5, 0.2, -0.1])
        dts = (0.2, 0.3, 0.4)
        phis = [dot(u, body.controlpoints[:, p]) for u in us, p in 1:body.ncells]
        for (step, u) in enumerate(us)
            body.velocity_kinematic .= w
            body.velocity .= u .- w
            pnl._pressure_bernoulli_phi_dot!(monitor, body, 1, (), (), u, step - 1, dts[step])
        end
        for p in 1:body.ncells
            expected_D = pnl._variable_bdf2(phis[3,p], phis[2,p], phis[1,p], dts[2], dts[1])
            @test isapprox(monitor.phi_dot[1][p], expected_D - dot(w, us[3]); atol=1e-12)
            @test isapprox(monitor.potential_history[1][p], phis[3,p]; atol=1e-12)
        end
        @test isempty(pnl.PressureBernoulli(1.0).potential_history)

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

    @testset "PressureBernoulli exterior potential limits" begin
        backends = (pnl.DirectBackend(),
            pnl.FastMultipoleBackend(expansion_order=12,
                multipole_acceptance=0.4, leaf_size=2))
        for backend in backends, kind in (:doublet, :ring, :source_doublet, :source_ring)
            body = make_bernoulli_limit_body(kind)
            monitor = pnl.PressureBernoulli(1.0; unsteady=true, backend)
            frames = pnl.ReferenceFrame(body)
            gammai = pnl.get_Gammai(body)

            raw1 = raw_controlpoint_potential(body, backend)
            monitor((body,), (nothing,), frames, zeros(3), 0, 0.2)
            exterior1 = raw1 .- body.strength[:, gammai]
            @test isapprox(monitor.potential_history[1], exterior1;
                atol=2e-9, rtol=2e-9)

            body.strength[:, gammai] .+= range(0.05, 0.15; length=body.ncells)
            raw2 = raw_controlpoint_potential(body, backend)
            monitor((body,), (nothing,), frames, zeros(3), 1, 0.2)
            exterior2 = raw2 .- body.strength[:, gammai]
            @test isapprox(monitor.potential_history[1], exterior2;
                atol=2e-9, rtol=2e-9)
            @test isapprox(monitor.phi_dot[1], (exterior2 .- exterior1) ./ 0.2;
                atol=2e-8, rtol=2e-8)
        end

        source = make_octa_source_body()
        raw = raw_controlpoint_potential(source, pnl.DirectBackend())
        monitor = pnl.PressureBernoulli(1.0; unsteady=true, backend=pnl.DirectBackend())
        monitor((source,), (nothing,), pnl.ReferenceFrame(source), zeros(3), 0, 0.1)
        @test isapprox(monitor.potential_history[1], raw; atol=1e-12)
    end


    @testset "PressureBernoulli steady relative projected trace" begin
        rho = 1.0
        uinf = [1.0, 0.0, 0.0]
        body = make_octa_source_body()
        body.velocity .= 0.0
        body.velocity_kinematic .= 0.0
        for p in 1:body.ncells
            body.velocity[:,p] .= [0.4, -0.2, 0.1] .+ 7.0 .* body.normals[:,p]
            body.velocity_kinematic[:,p] .= [0.3, 0.1, -0.2]
        end
        monitor = pnl.PressureBernoulli(rho; backend=pnl.DirectBackend(),
            correct_kuttacondition=false)
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), uinf, 0, 1.0)

        # Steady mode uses the tangential projection of the body-relative
        # velocity, with velocity_kinematic excluded.
        expected = similar(body.velocity)
        pnl._pressure_fill_relative_surface_velocity!(expected, body)
        @test monitor.inertial_velocity[1] == expected
        expected_pressure = [0.5 * rho * (norm(uinf)^2 - norm(expected[:, p])^2)
                             for p in 1:body.ncells]
        @test monitor.pressure[1] ≈ expected_pressure

        # Normal leakage is projected out: doubling the normal component of the
        # relative velocity leaves the steady pressure unchanged.
        for p in 1:body.ncells
            body.velocity[:,p] .+= 7.0 .* body.normals[:,p]
        end
        leak_monitor = pnl.PressureBernoulli(rho; backend=pnl.DirectBackend(),
            correct_kuttacondition=false)
        leak_monitor((body,), (nothing,), pnl.ReferenceFrame(body), uinf, 0, 1.0)
        @test leak_monitor.pressure[1] ≈ expected_pressure
        for p in 1:body.ncells
            body.velocity[:,p] .-= 7.0 .* body.normals[:,p]
        end

        # Changing velocity_kinematic must not change the steady pressure.
        body.velocity_kinematic .= 0.0
        for p in 1:body.ncells
            body.velocity_kinematic[:,p] .= [-5.0, 2.5, 4.0]
        end
        kin_monitor = pnl.PressureBernoulli(rho; backend=pnl.DirectBackend(),
            correct_kuttacondition=false)
        kin_monitor((body,), (nothing,), pnl.ReferenceFrame(body), uinf, 0, 1.0)
        @test kin_monitor.pressure[1] ≈ expected_pressure

        # Unsteady mode still uses the inertial trace,
        # tangent(body.velocity) + velocity_kinematic.
        expected_inertial = similar(body.velocity)
        pnl._pressure_fill_inertial_surface_velocity!(expected_inertial, body)
        @test expected_inertial ≈ expected .+ body.velocity_kinematic
        unsteady_monitor = pnl.PressureBernoulli(rho; unsteady=true,
            backend=pnl.DirectBackend(), correct_kuttacondition=false)
        unsteady_monitor((body,), (nothing,), pnl.ReferenceFrame(body), uinf, 0, 1.0)
        @test unsteady_monitor.inertial_velocity[1] == expected_inertial
    end

    @testset "PressureBernoulli steady rotor loading is first order" begin
        # Regression for the ef1fe1e steady kinetic-energy bug: on a moving
        # body, the inertial trace (relative + kinematic) reduces the surface
        # velocity to the small perturbation on both sides of a loaded pair,
        # cancelling the first-order loading. The relative form must retain it.
        rho = 1.2
        W = 30.0
        delta = 0.5
        body = make_octa_source_body()
        body.velocity .= 0.0
        body.velocity_kinematic .= 0.0
        tangents = similar(body.velocity)
        for p in 1:body.ncells
            n = body.normals[:, p]
            t = cross(n, [0.37, -0.61, 0.42])
            t = t / norm(t)
            tangents[:, p] .= t
            # Blade-frame picture: kinematic motion +W t, relative surface flow
            # -(W ± delta) t depending on the side of the pair.
            s = isodd(p) ? delta : -delta
            body.velocity[:, p] .= -(W + s) .* t
            body.velocity_kinematic[:, p] .= W .* t
        end
        monitor = pnl.PressureBernoulli(rho; backend=pnl.DirectBackend(),
            correct_kuttacondition=false)
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), [0.0, 0.0, 0.0], 0, 1.0)
        # First-order loading: P(-side) - P(+side) = 2 rho W delta.
        p_plus = monitor.pressure[1][1]   # panel with |u_rel| = W + delta
        p_minus = monitor.pressure[1][2]  # panel with |u_rel| = W - delta
        @test isapprox(p_minus - p_plus, 2 * rho * W * delta; rtol=1e-12)
        # The broken inertial steady form gives an exactly zero difference here.
        broken = zeros(body.ncells)
        inertial = similar(body.velocity)
        pnl._pressure_fill_inertial_surface_velocity!(inertial, body)
        pnl.calcfield_P!(broken, body, inertial, 0.0, rho, nothing;
            correct_kuttacondition=false)
        @test isapprox(broken[2] - broken[1], 0.0; atol=1e-10)
    end

    @testset "Unsteady history order and vector-source policy" begin
        # Variable-step BDF2 differentiates quadratics exactly.
        h, hp = 0.3, 0.2
        f(t) = 2t^2 - 3t + 1
        @test isapprox(pnl._variable_bdf2(f(1.0), f(1-h), f(1-h-hp), h, hp), 1.0; atol=1e-14)

        # Cubic truncation error is second order on equal-step refinement.
        cubic_error(h) = abs(pnl._variable_bdf2(1.0, (1-h)^3, (1-2h)^3, h, h) - 3.0)
        @test isapprox(cubic_error(0.1) / cubic_error(0.05), 4.0; rtol=1e-10)

        body = make_plate_vortex_body()
        wake = make_active_final_filament_wake(body)
        frames = pnl.ReferenceFrame(body)
        default_monitor = pnl.PressureBernoulli(1.0; unsteady=true,
            backend=pnl.DirectBackend())
        @test_throws ArgumentError default_monitor((body,), (wake,), frames, zeros(3), 0, 0.1)

        partial = pnl.PressureBernoulli(1.0; unsteady=true, allow_partial=true,
            backend=pnl.DirectBackend())
        @test_logs (:warn, r"allow_partial=true") partial((body,), (wake,), frames, zeros(3), 0, 0.1)
        @test_logs partial((body,), (wake,), frames, zeros(3), 1, 0.1)
        @test all(isfinite, partial.pressure[1])

        vector_body = pnl.RigidWakeBody{
            Union{pnl.VortexRing, pnl.UniformVortexSheet}, 2, Float64, true}(
            copy(body.nodes), copy(body.cells), deepcopy(body.shedding);
            check_mesh=false, watertight=false, ensure_winding=false)
        pnl.calc_normals!(vector_body); pnl.calc_controlpoints!(vector_body)
        rejected = pnl.PressureBernoulli(1.0; unsteady=true, allow_partial=true,
            backend=pnl.DirectBackend())
        @test_throws ArgumentError rejected((vector_body,), (nothing,),
            pnl.ReferenceFrame(vector_body), zeros(3), 0, 0.1)
    end


    @testset "PressureBernoulli scalar ALE split and Galilean translation" begin
        body = make_plate_vortex_body()
        wake = make_active_final_filament_wake(body)
        excluded = last(pnl.get_sources(wake))
        probes = FastMultipole.ProbeSystem(body.ncells, Float64)
        for p in 1:body.ncells
            probes.position[p] = SVector{3, Float64}(body.controlpoints[:, p])
            probes.gradient[p] = zero(eltype(probes.gradient))
        end
        pnl.influence!((probes,), (excluded,), pnl.DirectBackend();
            precalc=false, scalar_potential=false, gradient=true, hessian=(false,))
        excluded_velocity = reduce(hcat, probes.gradient)
        w = [0.35, -0.2, 0.15]
        retained = [0.6, 0.1, -0.25]
        body.velocity_kinematic .= w
        body.velocity .= retained .+ excluded_velocity .- w
        monitor = pnl.PressureBernoulli(1.0; unsteady=true, allow_partial=true,
            backend=pnl.DirectBackend())
        frames = pnl.ReferenceFrame(body)
        @test_logs (:warn, r"partial diagnostic") monitor((body,), (wake,), frames, zeros(3), 0, 0.1)
        @test_logs monitor((body,), (wake,), frames, zeros(3), 1, 0.1)
        @test isapprox(monitor.phi_dot[1], fill(-dot(w, retained), body.ncells);
            atol=2e-10, rtol=2e-10)
        total_surface = similar(body.velocity)
        pnl._pressure_fill_inertial_surface_velocity!(total_surface, body)
        kinetic_only = zeros(body.ncells)
        pnl.calcfield_P!(kinetic_only, body, total_surface, 0.0, 1.0, nothing)
        @test isapprox(monitor.pressure[1], kinetic_only .- monitor.phi_dot[1];
            atol=2e-10, rtol=2e-10)

        translating = make_octa_source_body()
        u = [0.4, -0.3, 0.2]
        translating.velocity_kinematic .= u
        translating.velocity .= 0.0
        galilean = pnl.PressureBernoulli(1.0; unsteady=true,
            backend=pnl.DirectBackend())
        galilean((translating,), (nothing,), pnl.ReferenceFrame(translating), u, 0, 0.25)
        translating.nodes .+= u .* 0.25
        pnl.calc_controlpoints!(translating)
        translating.velocity_kinematic .= u
        translating.velocity .= 0.0
        galilean((translating,), (nothing,), pnl.ReferenceFrame(translating), u, 1, 0.25)
        @test isapprox(galilean.phi_dot[1], zeros(translating.ncells); atol=2e-10)
        @test isapprox(galilean.pressure[1], zeros(translating.ncells); atol=2e-10)
    end

    @testset "Pressure Kutta averaging is opt-in" begin
        body = make_plate_vortex_body()
        U = zeros(3, body.ncells)
        U[1, :] .= (1.0, 2.0)
        default_pressure = zeros(body.ncells)
        averaged_pressure = zeros(body.ncells)
        pnl.calcfield_P!(default_pressure, body, U, 0.0, 1.0, nothing)
        pnl.calcfield_P!(averaged_pressure, body, U, 0.0, 1.0, nothing;
            correct_kuttacondition=true)
        @test default_pressure[1] != default_pressure[2]
        @test averaged_pressure[1] == averaged_pressure[2] ==
            (default_pressure[1] + default_pressure[2]) / 2
        @test !pnl.PressureBernoulli(1.0).correct_kuttacondition
    end

    @testset "PressureLaplace monitor" begin
        body = make_octa_source_body()
        monitor = pnl.PressureLaplace((body,), 1.2; reference_panel=1, reference_pressure=0.0)

        @test pnl.monitor_provides(monitor) == (:P,)
        @test !monitor.unsteady
        @test monitor.acceleration_form == :material_derivative
        @test monitor.gradient_mode == :edge_difference
        @test !pnl.monitor_requires_body_hessian(monitor)
        @test pnl.monitor_requires_induced_vorticity(monitor)
        @test pnl.audit_monitors((monitor, pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization()))) !== nothing
        @test_throws ArgumentError pnl.PressureLaplace(1.0)
        @test_throws ArgumentError pnl.PressureLaplace((body,), 1.0; gradient_mode=:unknown)
        @test_throws ArgumentError pnl.PressureLaplace((body,), 1.0; acceleration_form=:unknown)

        lamb = @test_logs (:warn, r"deprecated") pnl.PressureLaplace((body,), 1.2;
            reference_panel=1, gradient_mode=:edge_difference,
            acceleration_form=:lamb_vector)
        @test lamb.acceleration_form == :lamb_vector
        @test !pnl.monitor_requires_body_hessian(lamb)
        @test pnl.monitor_requires_induced_vorticity(lamb)

        body.velocity .= 0.0
        body.velocity[1, :] .= 1.0
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.25)
        @test isempty(monitor.velocity_dot)
        @test all(isfinite, monitor.p[1])

        body.velocity[1, :] .+= 0.5
        body.velocity[2, :] .-= 0.25
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 1, 0.25)
        @test isempty(monitor.velocity_dot)
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
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.25)

        @test all(isfinite, monitor.p[1])
        @test all(isfinite, monitor.acceleration[1])
        @test all(isfinite, monitor.surface_velocity_gradient[1])

        lamb = pnl.PressureLaplace((body,), 1.0;
            reference_panel=1,
            gradient_mode=:surface_velocity,
            acceleration_form=:lamb_vector)
        body.velocity_gradient .= NaN
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
        monitor = @test_logs (:warn, r"unsupported diagnostic") pnl.PressureLaplace((body,), 1.0;
            gradient_mode=:corrected_hessian)
        x = collect(range(-0.4, 0.7; length=body.ncells))
        @test monitor.pressure_operator[1].n == body.ncells
        @test monitor.pressure_operator[1].reference_panel == 1
        @test isapprox(monitor.pressure_operator[1].row_diagonal,
            [monitor.L[1][i, i] for i in 1:body.ncells]; atol=1e-12)
        y_op = zeros(body.ncells)
        mul!(y_op, monitor.pressure_operator[1], x)
        @test isapprox(y_op, monitor.L[1] * x; atol=1e-12)

        y_scaled = collect(range(0.2, 1.1; length=body.ncells))
        y_expected = 0.25 .* y_scaled .+ 1.7 .* (monitor.L[1] * x)
        mul!(y_scaled, monitor.pressure_operator[1], x, 1.7, 0.25)
        @test isapprox(y_scaled, y_expected; atol=1e-12)
        y_call = zeros(body.ncells)
        monitor.pressure_operator[1](y_call, x, 1.0, 0.0)
        @test isapprox(y_call, monitor.L[1] * x; atol=1e-12)

        monitor.b[1] .= b
        pnl._pressure_solve!(monitor, 1)
        @test isapprox(monitor.p[1], p_exact; atol=1e-8)
        @test monitor.workspace[1].stats.solved
        @test monitor.absolute_residual[1] ≤ 1e-8
        @test monitor.relative_residual[1] ≤ 1e-8

        failed = pnl.PressureLaplace((body,), 1.0;
            gradient_mode=:surface_velocity, itmax=1)
        failed.b[1] .= L * collect(range(0.0, 2.0; length=body.ncells))
        @test_logs (:warn, r"did not converge") pnl._pressure_solve!(failed, 1)
        @test !failed.workspace[1].stats.solved
        @test failed.convergence_warned
        @test isfinite(failed.absolute_residual[1])
        @test isfinite(failed.relative_residual[1])
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
        first_operator = pressure.pressure_operator[1]
        pressure((body,), (nothing,), frames, zeros(3), 1, 0.1)
        @test pressure.L[1] === first_L
        @test pressure.pressure_operator[1] === first_operator

        body.nodes[1, :] .+= 0.2
        body.nodes[2, :] .-= 0.1
        pressure((body,), (nothing,), frames, zeros(3), 2, 0.1)
        @test pressure.L[1] === first_L
        @test pressure.pressure_operator[1] === first_operator

        old = body.nodes[1, 1]
        body.nodes[1, 1] = old + 0.1
        pressure((body,), (nothing,), frames, zeros(3), 3, 0.1)
        @test pressure.L[1] === first_L
        @test pressure.pressure_operator[1] === first_operator

        rebuilding = pnl.PressureLaplace((body,), 1.0;
            reference_panel=1,
            rebuild_every_step=true)
        body.velocity[1, :] .= 0.2 .* (1:body.ncells)
        rebuilding((body,), (nothing,), frames, zeros(3), 4, 0.1)
        rebuild_first_L = rebuilding.L[1]
        rebuild_first_operator = rebuilding.pressure_operator[1]
        rebuilding((body,), (nothing,), frames, zeros(3), 5, 0.1)
        @test rebuilding.L[1] !== rebuild_first_L
        @test rebuilding.pressure_operator[1] !== rebuild_first_operator

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
        @test_throws ArgumentError pnl.SpanwiseLoadingMonitor(2, 1;
            components=(lift=[0.0, 0.0, 1.0],), binning=:area_clip)

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

        overlap_body = make_spanwise_overlap_body()
        overlap_force = [0.0; 0.0; 10.0;;]
        ctx = pnl.MonitorContext()
        pnl.monitor_register!(ctx, :F, 1, overlap_force)
        overlap = pnl.SpanwiseLoadingMonitor(2, 1;
            components=(lift=[0.0, 0.0, 1.0],),
            span_axis=[0.0, 1.0, 0.0],
            binning=:span_overlap)
        pnl._run_monitor!(overlap, ctx, (overlap_body,), (nothing,),
            pnl.ReferenceFrame(overlap_body), zeros(3), 0, 0.1)
        @test overlap.counts == [1, 1]
        @test overlap.panel_bin_id == [1]
        @test isapprox(overlap.bin_center, [0.375, 0.625]; atol=1e-12)
        @test isapprox(overlap.force_components[1, :], [5.0, 5.0]; atol=1e-12)
        @test isapprox(sum(overlap.force_components[1, :]), 10.0; atol=1e-12)

        degenerate = pnl.SpanwiseLoadingMonitor(2, 1;
            components=(lift=[0.0, 0.0, 1.0],),
            span_axis=[0.0, 1.0, 0.0],
            binning=:span_overlap)
        ctx = pnl.MonitorContext()
        pnl.monitor_register!(ctx, :F, 1, [0.0; 0.0; 3.0;;])
        one_panel = make_spanwise_loading_body([0.25])
        pnl._run_monitor!(degenerate, ctx, (one_panel,), (nothing,),
            pnl.ReferenceFrame(one_panel), zeros(3), 0, 0.1)
        @test degenerate.counts == [1, 0]
        @test degenerate.panel_bin_id == [1]
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

        csvmon = pnl.SpanwiseLoadingMonitor(1, 1;
            components=(lift=[0.0, 0.0, 1.0], drag=[1.0, 0.0, 0.0]))
        csvdir = joinpath(mktempdir(), "monitors")
        csv = joinpath(csvdir, "post_monitor01_spanwise_system1.csv")
        ctx = pnl.MonitorContext()
        pnl.monitor_register!(ctx, :F, 1, [0.0 0.0; 0.0 0.0; 1.0 1.0])
        pnl._run_monitor!(csvmon, ctx, (local_body,), (nothing,), pnl.ReferenceFrame(local_body), zeros(3), 0, 0.1, 1.25)
        pnl.write_monitor_csv!(csvmon, csvdir, "post", 1, ctx, (local_body,), 0, 0.1; overwrite=true)
        pnl._run_monitor!(csvmon, ctx, (local_body,), (nothing,), pnl.ReferenceFrame(local_body), zeros(3), 1, 0.1, 1.50)
        pnl.write_monitor_csv!(csvmon, csvdir, "post", 1, ctx, (local_body,), 1, 0.1)
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

    @testset "BoundCirculationMonitor shedding local-node indices" begin
        nodes, cells = make_seeded_te_mesh()
        final_cells = pnl.ensure_consistent_winding(nodes, cells; watertight=false)
        bbox = ([0.8, -0.1, -0.1], [1.1, 2.1, 0.1])
        shedding = pnl.calc_shedding_from_seed(
            nodes, final_cells, 1, 2; bbox, end_node=3)
        body = pnl.RigidWakeBody{pnl.VortexRing}(
            nodes, final_cells, shedding;
            check_mesh=false, watertight=false, ensure_winding=false)
        body.strength[:, 1] .= 1.0
        pnl.calc_controlpoints!(body)

        monitor = pnl.BoundCirculationMonitor(body, 1, 1;
            i_frame=1, radial_dimension=2, R=2.0)
        monitor((body,), (nothing,), pnl.ReferenceFrame(body), zeros(3), 0, 0.1)

        @test monitor.r_over_R[:, 1] ≈ [0.25, 0.75]
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

    @testset "Monitor CSV hooks" begin
        dir = joinpath(mktempdir(), "monitors")

        body = make_octa_source_body()
        pressure = zeros(body.ncells)
        ctx = pnl.MonitorContext()
        pnl.monitor_set_time!(ctx, 0.0)
        pnl.monitor_register!(ctx, :P, 1, pressure)
        force = pnl.ForceMonitor(2, 1; normalization=pnl.NoNormalization())
        pnl._run_monitor!(force, ctx, (body,), (nothing,), pnl.ReferenceFrame(body), [1.0, 0.0, 0.0], 0, 0.1)
        pnl.write_monitor_csv!(force, dir, "case", 1, ctx, (body,), 0, 0.1; overwrite=true)
        pnl.monitor_set_time!(ctx, 0.1)
        pnl._run_monitor!(force, ctx, (body,), (nothing,), pnl.ReferenceFrame(body), [1.0, 0.0, 0.0], 1, 0.1)
        pnl.write_monitor_csv!(force, dir, "case", 1, ctx, (body,), 1, 0.1)
        rows = readlines(joinpath(dir, "case_monitor01_force_system1.csv"))
        @test rows[1] == "step,time,CFx,CFy,CFz,CMx,CMy,CMz"
        @test length(rows) == 3

        doublet = make_surface_vorticity_gradient_body()
        doublet.velocity[1, :] .= 1.0
        surface = pnl.SurfaceVorticityForce(doublet, 1, 1; normalization=pnl.NoNormalization())
        surface((doublet,), (nothing,), pnl.ReferenceFrame(doublet), zeros(3), 0, 0.1)
        pnl.write_monitor_csv!(surface, dir, "case", 2, ctx, (doublet,), 0, 0.1; overwrite=true)
        rows = readlines(joinpath(dir, "case_monitor02_surface_vorticity_force_system1.csv"))
        @test rows[1] == "step,time,CFx,CFy,CFz,CMx,CMy,CMz"
        @test length(rows) == 2

        kj_body = make_plate_vortex_body()
        kj_body.strength[:, 1] .= [0.8, -0.35]
        kj = pnl.KuttaJoukowskiForce(kj_body, 1, 1;
            backend=pnl.DirectBackend(),
            normalization=pnl.NoNormalization())
        kj((kj_body,), (nothing,), pnl.ReferenceFrame(kj_body), zeros(3), 0, 0.1)
        pnl.write_monitor_csv!(kj, dir, "case", 3, ctx, (kj_body,), 0, 0.1; overwrite=true)
        rows = readlines(joinpath(dir, "case_monitor03_kutta_joukowski_force_system1.csv"))
        @test rows[1] == "step,time,CFx,CFy,CFz"
        @test length(rows) == 2

        bound = pnl.BoundCirculationMonitor(kj_body, 1, 1;
            i_frame=1,
            radial_dimension=1,
            R=1.0,
            section_tol=0.6)
        bound((kj_body,), (nothing,), pnl.ReferenceFrame(kj_body), zeros(3), 0, 0.1)
        pnl.write_monitor_csv!(bound, dir, "case", 4, ctx, (kj_body,), 0, 0.1; overwrite=true)
        rows = readlines(joinpath(dir, "case_monitor04_bound_circulation_system1.csv"))
        @test rows[1] == "step,time,blade,section,r_over_R,circulation_te,circulation_slice"
        @test length(rows) == 2

        lap_body = make_skewed_two_panel_body()
        laplace = pnl.PressureLaplace((lap_body,), 1.0; reference_panel=1)
        laplace((lap_body,), (nothing,), pnl.ReferenceFrame(lap_body), zeros(3), 0, 0.1)
        pnl.write_monitor_csv!(laplace, dir, "case", 5, ctx, (lap_body,), 0, 0.1; overwrite=true)
        rows = readlines(joinpath(dir, "case_monitor05_pressure_laplace_system1.csv"))
        @test rows[1] == "step,time,system,panels,rebuild,cg_iters,cg_solved,absolute_residual,relative_residual"
        @test occursin(",1,2,false,", rows[2])

        bernoulli = pnl.PressureBernoulli(1.0)
        pnl.write_monitor_csv!(bernoulli, dir, "case", 6, ctx, (body,), 0, 0.1; overwrite=true)
        @test !isfile(joinpath(dir, "case_monitor06_pressure_bernoulli_system1.csv"))

        suppressed = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization(), file=false)
        pnl.write_monitor_csv!(suppressed, dir, "case", 7, ctx, (body,), 0, 0.1; overwrite=true)
        @test !isfile(joinpath(dir, "case_monitor07_force_system1.csv"))
    end

    @testset "PressureLaplace Bernoulli constant-field comparison" begin
        body_b = make_octa_source_body()
        body_l = make_octa_source_body()
        body_b.velocity[1, :] .= 0.4
        body_b.velocity[2, :] .= -0.2
        body_l.velocity .= body_b.velocity

        laplace = pnl.PressureLaplace((body_l,), 1.0; reference_panel=1)
        p_b_raw = zeros(body_b.ncells)
        u_b = similar(body_b.velocity)
        pnl._pressure_fill_inertial_surface_velocity!(u_b, body_b)
        pnl.calcfield_P!(p_b_raw, body_b, u_b, 1.0, 1.0, zeros(body_b.ncells);
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

    @testset "PressureLaplace lamb_vorticity variants" begin
        body = make_skewed_two_panel_body()
        @test_throws ArgumentError pnl.PressureLaplace((body,), 1.2;
            lamb_vorticity=:unknown)
        @test_throws ArgumentError pnl.PressureLaplace((body,), 1.2;
            kappa_basis=:unknown)

        # :hessian_curl flips the Hessian requirement even for the lamb form.
        hess = pnl.PressureLaplace((body,), 1.2; acceleration_form=:lamb_vector,
            lamb_vorticity=:hessian_curl)
        @test pnl.monitor_requires_body_hessian(hess)
        @test pnl.monitor_requires_induced_vorticity(hess)
        base = pnl.PressureLaplace((body,), 1.2; acceleration_form=:lamb_vector)
        @test base.lamb_vorticity == :induced
        @test !pnl.monitor_requires_body_hessian(base)

        # Baseline :induced ingests body.induced_vorticity verbatim, so the RHS
        # is bit-identical to the pre-kwarg behavior.
        for p in 1:body.ncells
            body.velocity[:, p] .= [0.4 + 0.05p, -0.15, 0.1p]
        end
        body.induced_vorticity[1, :] .= 0.3
        body.induced_vorticity[3, :] .= -0.8
        pnl._pressure_rhs_from_lamb_vector!(base.b[1], base, body, 1, nothing)
        @test base.omega_used[1] == body.induced_vorticity

        # Source-only body: kappa == 0, so :no_bound == :induced and
        # :bound_only gives an omega-free RHS.
        nob = pnl.PressureLaplace((body,), 1.2; acceleration_form=:lamb_vector,
            lamb_vorticity=:no_bound)
        pnl._pressure_rhs_from_lamb_vector!(nob.b[1], nob, body, 1, nothing)
        @test isapprox(nob.b[1], base.b[1]; atol=1e-14)
        bnd = pnl.PressureLaplace((body,), 1.2; acceleration_form=:lamb_vector,
            lamb_vorticity=:bound_only)
        pnl._pressure_rhs_from_lamb_vector!(bnd.b[1], bnd, body, 1, nothing)
        @test all(iszero, bnd.omega_used[1])

        # Doublet body: kappa != 0 and the variants decompose additively:
        # omega(:no_bound) + omega(:bound_only) == omega(:induced).
        doublet_body = make_skewed_two_panel_doublet_body()
        doublet_body.velocity .= 0.0
        doublet_body.velocity[1, :] .= 0.7
        doublet_body.induced_vorticity .= 0.0
        doublet_body.induced_vorticity[2, :] .= 0.4     # fake wake-induced part
        pnl._add_bound_surface_vorticity!(doublet_body) # += kappa (tri basis)
        # kappa_basis must match the basis used to accumulate kappa above.
        variants = Dict(sym => pnl.PressureLaplace((doublet_body,), 1.2;
                acceleration_form=:lamb_vector, lamb_vorticity=sym,
                kappa_basis=:tri)
            for sym in (:induced, :no_bound, :bound_only))
        for m in values(variants)
            pnl._pressure_rhs_from_lamb_vector!(m.b[1], m, doublet_body, 1, nothing)
        end
        @test !all(iszero, variants[:bound_only].omega_used[1])
        @test isapprox(variants[:no_bound].omega_used[1] .+ variants[:bound_only].omega_used[1],
            variants[:induced].omega_used[1]; atol=1e-12)
        # and :no_bound recovers exactly the fake wake-induced part
        @test isapprox(variants[:no_bound].omega_used[1][2, :], fill(0.4, doublet_body.ncells);
            atol=1e-12)
        @test isapprox(variants[:no_bound].omega_used[1][[1, 3], :],
            zeros(2, doublet_body.ncells); atol=1e-12)

        # :hessian_curl reads the antisymmetric part of body.velocity_gradient
        # (+ kappa, zero for the source body).
        body.velocity_gradient .= 0.0
        for p in 1:body.ncells
            body.velocity_gradient[2, 1, p] = 0.5   # d(u2)/d(x1)
            body.velocity_gradient[1, 2, p] = -0.25 # d(u1)/d(x2)
        end
        hc = pnl.PressureLaplace((body,), 1.2; acceleration_form=:lamb_vector,
            lamb_vorticity=:hessian_curl)
        pnl._pressure_rhs_from_lamb_vector!(hc.b[1], hc, body, 1, nothing)
        @test isapprox(hc.omega_used[1][3, :], fill(0.75, body.ncells); atol=1e-14)
        @test isapprox(hc.omega_used[1][1:2, :], zeros(2, body.ncells); atol=1e-14)
    end

    @testset "PressureLaplace rotational edge correction" begin
        body = make_skewed_two_panel_body()
        edge = pnl.PressureLaplace((body,), 1.2; reference_panel=1,
            gradient_mode=:edge_difference)
        @test !pnl.monitor_requires_body_hessian(edge)
        @test pnl.monitor_requires_induced_vorticity(edge)
        body.velocity[1, :] .= 0.7
        body.velocity[2, :] .= -0.2
        pnl._pressure_fill_inertial_surface_velocity!(edge.u_inertial[1], body)
        body.induced_vorticity .= 0.0
        pnl._pressure_rhs_from_edge_material_derivative!(edge.b[1], edge, body, 1, nothing)
        symmetric_rhs = copy(edge.b[1])
        body.induced_vorticity[3, :] .= 1.0
        pnl._pressure_rhs_from_edge_material_derivative!(edge.b[1], edge, body, 1, nothing)
        @test !isapprox(edge.b[1], symmetric_rhs; atol=1e-12)

        # For a nonsymmetric manufactured gradient, the corrected contraction
        # recovers r'Gq; naive q'du alone generally recovers q'Gr instead.
        G = [0.2 1.1 0.0; -0.4 0.3 0.0; 0.0 0.0 0.0]
        q = [0.7, -0.2, 0.0]
        edge_a, edge_b, i, j = edge.edges[1][:,1]
        for p in 1:body.ncells
            edge.u_inertial[1][:,p] .= G * body.controlpoints[:,p]
            body.velocity[:,p] .= q
        end
        omega_manufactured =
            [G[3,2] - G[2,3], G[1,3] - G[3,1], G[2,1] - G[1,2]]
        body.induced_vorticity .= omega_manufactured
        pnl._pressure_rhs_from_edge_material_derivative!(edge.b[1], edge,
            body, 1, nothing)
        r = body.controlpoints[:,j] - body.controlpoints[:,i]
        qi = q - dot(q, body.normals[:,i]) * body.normals[:,i]
        qj = q - dot(q, body.normals[:,j]) * body.normals[:,j]
        qedge = 0.5 * (qi + qj)
        w = pnl._pressure_edge_conormal_weight(body, edge_a, edge_b, i, j)[1]
        corrected = edge.rho * w * dot(r, G * qedge)
        naive = edge.rho * w * dot(qedge, G * r)
        nonreference_value = i == edge.reference_panel ? -edge.b[1][j] : edge.b[1][i]
        @test nonreference_value ≈ corrected atol=1e-12
        @test !isapprox(corrected, naive; atol=1e-12)

        nodes, cells = make_planar_gradient_mesh()
        doublet = pnl.NonLiftingBody{pnl.ConstantDoublet}(nodes, cells;
            watertight=false, ensure_winding=false)
        pnl.calc_normals!(doublet); pnl.calc_controlpoints!(doublet)
        doublet.strength[:,1] .= range(-1.0, 1.0; length=doublet.ncells)
        wake_curl = repeat([0.2, -0.3, 0.4], 1, doublet.ncells)
        pnl._bound_surface_vorticity!(doublet.induced_vorticity, doublet;
            grad_mu_options=(; basis=:quad))
        doublet.induced_vorticity .+= wake_curl
        edge_doublet = pnl.PressureLaplace((doublet,), 1.0;
            gradient_mode=:edge_difference)
        omega = pnl._pressure_fill_edge_exterior_vorticity!(edge_doublet,
            doublet, 1)
        @test isapprox(omega, wake_curl; atol=1e-12)
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
        steady((body_steady,), (nothing,), pnl.ReferenceFrame(body_steady), zeros(3), 0, 0.25)
        unsteady((body_unsteady,), (nothing,), pnl.ReferenceFrame(body_unsteady), zeros(3), 0, 0.25)

        @test isempty(steady.velocity_dot)
        @test isapprox(unsteady.velocity_dot[1], zeros(size(body_unsteady.velocity)); atol=1e-12)
        @test isapprox(steady.b[1], unsteady.b[1]; atol=1e-12)

        body_unsteady.velocity[1, :] .+= 0.4
        unsteady((body_unsteady,), (nothing,), pnl.ReferenceFrame(body_unsteady), zeros(3), 1, 0.5)
        @test !isapprox(unsteady.velocity_dot[1], zeros(size(body_unsteady.velocity)); atol=1e-12)

        # Reusing a monitor from an earlier/nonconsecutive step safely reseeds.
        body_unsteady.velocity[2, :] .+= 0.7
        unsteady((body_unsteady,), (nothing,), pnl.ReferenceFrame(body_unsteady), zeros(3), 0, 0.25)
        @test isapprox(unsteady.velocity_dot[1], zeros(size(body_unsteady.velocity)); atol=1e-12)
    end

    @testset "PressureLaplace ALE rotation, final projection, and jump gradient" begin
        nodes, cells = make_planar_gradient_mesh()
        body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells;
            watertight=false, ensure_winding=false)
        pnl.calc_normals!(body); pnl.calc_controlpoints!(body)
        monitor = pnl.PressureLaplace((body,), 1.0;
            gradient_mode=:corrected_hessian)

        # Manufactured translating-grid ALE field: a = D_g u + Gq.
        translating = pnl.PressureLaplace((body,), 1.0; unsteady=true,
            gradient_mode=:corrected_hessian)
        Dg = [0.3, -0.4, 0.5]
        Gm = [0.2 -0.1 0.7; 0.6 0.3 -0.2; -0.4 0.8 0.1]
        q = [0.9, -0.25, 0.0]
        translating.velocity_dot[1] .= Dg
        body.velocity .= q
        body.velocity_kinematic .= [0.2, 0.1, -0.3]
        for p in 1:body.ncells
            body.velocity_gradient[:,:,p] .= Gm
        end
        pnl._pressure_fill_inertial_surface_velocity!(translating.u_inertial[1], body)
        pnl._pressure_material_acceleration!(translating, body, 1)
        expected_ale = Dg + Gm * q
        @test all(isapprox.(translating.acceleration[1], expected_ale; atol=1e-12))
        @test all(isapprox.(translating.tangential[1][1:2,:],
            expected_ale[1:2]; atol=1e-12))
        @test all(abs.(translating.tangential[1][3,:]) .< 1e-12)

        # Uniform inertial velocity on a rotating grid has zero spatial
        # gradient. body.angular_velocity must not be added to the Hessian.
        Ω = [0.0, 0.0, 2.0]
        uI = [0.4, -0.3, 0.0]
        body.angular_velocity .= Ω
        for p in 1:body.ncells
            w = cross(Ω, body.controlpoints[:,p])
            body.velocity_kinematic[:,p] .= w
            body.velocity[:,p] .= uI .- w
        end
        body.velocity_gradient .= 0.0
        pnl._pressure_fill_inertial_surface_velocity!(monitor.u_inertial[1], body)
        pnl._pressure_material_acceleration!(monitor, body, 1)
        @test isapprox(monitor.acceleration[1], zeros(size(body.velocity)); atol=1e-12)

        # A normal component is retained while forming G*q, then removed once
        # from the completed acceleration.
        body.velocity_kinematic .= 0.0
        body.velocity .= 0.0
        body.velocity[1,:] .= 1.0
        body.velocity_gradient .= 0.0
        body.velocity_gradient[:,1,:] .= reshape([1.0, 2.0, 3.0], 3, 1)
        pnl._pressure_fill_inertial_surface_velocity!(monitor.u_inertial[1], body)
        pnl._pressure_material_acceleration!(monitor, body, 1)
        @test all(abs.(monitor.tangential[1][3,:]) .< 1e-12)
        @test all(isapprox.(monitor.tangential[1][1,:], 1.0; atol=1e-12))
        @test all(isapprox.(monitor.tangential[1][2,:], 2.0; atol=1e-12))

        doublet = pnl.NonLiftingBody{pnl.ConstantDoublet}(nodes, cells;
            watertight=false, ensure_winding=false)
        pnl.calc_normals!(doublet); pnl.calc_controlpoints!(doublet)
        # Nonconstant paired-triangle μ data produces a nonuniform half-jump.
        doublet.strength[:,1] .= [0.0, 2.0, -1.0, 5.0, 1.5, -3.0, 4.0, -2.5]
        doublet.velocity .= 0.0; doublet.velocity[1,:] .= 0.6
        doublet.velocity_gradient .= 0.0
        jump_monitor = pnl.PressureLaplace((doublet,), 1.0;
            gradient_mode=:corrected_hessian)
        pnl._pressure_fill_inertial_surface_velocity!(jump_monitor.u_inertial[1], doublet)
        pnl._pressure_material_acceleration!(jump_monitor, doublet, 1)
        @test any(abs.(jump_monitor.jump_velocity_gradient[1]) .> 1e-10)
        @test any(abs.(jump_monitor.acceleration[1]) .> 1e-10)
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

    @testset "compute_mu_gradient! option normalization" begin
        nodes, cells = make_planar_gradient_mesh()
        body = pnl.NonLiftingBody{pnl.ConstantSource}(nodes, cells;
            watertight=false,
            ensure_winding=false)

        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)

        mu = [0.0, 2.0, -1.0, 5.0, 1.5, -3.0, 4.0, -2.5]
        te_info = zeros(Int, 2, body.ncells)

        @test pnl._normalize_grad_mu_options((; basis=:tri)).basis === :tri
        @test pnl._normalize_grad_mu_options((; basis=:quad)).basis === :quad
        @test pnl._normalize_grad_mu_options((; basis=:quad)).quad_grow === true
        @test pnl._normalize_grad_mu_options((; basis=:quad)).quad_grow_max_depth == 2
        @test pnl._normalize_grad_mu_options((;); default_basis=:quad).basis === :quad
        @test pnl._normalize_grad_mu_options((; tri_robust_ar_threshold=12.0)).tri_robust_ar_threshold == 12.0
        for old_key in (:mode, :robust, :robust_ar_threshold, :bfs_enabled,
                        :bfs_stop, :bfs_cond_max, :bfs_max_depth,
                        :bfs_target_healthy, :ar_threshold,
                        :gradient_quad_consistent)
            @test_throws ArgumentError pnl._normalize_grad_mu_options(NamedTuple{(old_key,)}((true,)))
        end
        @test_throws ArgumentError pnl._normalize_grad_mu_options((; basis=:quad, tri_robust=false))
        @test_throws ArgumentError pnl._normalize_grad_mu_options((; basis=:quad, tri_robust_ar_threshold=10.0))
        @test_throws ArgumentError pnl._normalize_grad_mu_options((; basis=:quad, tri_robust_max_depth=4))
        @test_throws ArgumentError pnl._normalize_grad_mu_options((; basis=:quad, tri_robust_target_healthy=6))
        @test_throws MethodError pnl.compute_mu_gradient!(zeros(3, body.ncells),
            body.controlpoints, body.normals, body.cells, body.neighbor, mu, te_info;
            scale=1.0, nodes=body.nodes, quad_mode=:mu_diff)
        @test_throws MethodError pnl.compute_mu_gradient!(zeros(3, body.ncells),
            body.controlpoints, body.normals, body.cells, body.neighbor, mu, te_info;
            scale=1.0, bad_panel_mask=falses(body.ncells))
        @test_throws ArgumentError pnl.compute_mu_gradient!(zeros(3, body.ncells),
            body.controlpoints, body.normals, body.cells, body.neighbor, mu, te_info;
            scale=1.0, grad_mu_options=(; basis=:quad))
        @test_throws ArgumentError pnl.compute_mu_gradient!(zeros(3, body.ncells),
            body.controlpoints, body.normals, body.cells, body.neighbor, mu, te_info;
            scale=1.0, grad_mu_options=(; basis=:tri, tri_robust=true))

        default_quad = zeros(3, body.ncells)
        canonical_tri = zeros(3, body.ncells)
        pnl.compute_mu_gradient!(default_quad, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info;
            scale=1.0, nodes=body.nodes)
        pnl.compute_mu_gradient!(canonical_tri, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info;
            scale=1.0, nodes=body.nodes, grad_mu_options=(; basis=:tri))
        @test !isapprox(default_quad, canonical_tri; atol=1e-12)

        canonical_mu_diff = zeros(3, body.ncells)
        pnl.compute_mu_gradient!(canonical_mu_diff, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info;
            scale=1.0, nodes=body.nodes,
            grad_mu_options=(; basis=:quad,
                quad_grow=true, quad_grow_stop=:depth,
                quad_grow_cond_max=12.0, quad_grow_max_depth=2))
        @test all(isfinite, canonical_mu_diff)

        canonical_mu_diff_nogrow = zeros(3, body.ncells)
        pnl.compute_mu_gradient!(canonical_mu_diff_nogrow, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info;
            scale=1.0, nodes=body.nodes,
            grad_mu_options=(; basis=:quad, quad_grow=false))
        @test all(isfinite, canonical_mu_diff_nogrow)

        bad = falses(body.ncells)
        bad[1] = true
        bfs_blocked = zeros(3, body.ncells)
        bfs_grown = zeros(3, body.ncells)
        pnl._compute_mu_gradient_masked!(bfs_blocked, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info;
            scale=1.0, bad_panel_mask=bad,
            grad_mu_options=(; basis=:tri, tri_robust_max_depth=0, tri_robust_target_healthy=1))
        pnl._compute_mu_gradient_masked!(bfs_grown, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info;
            scale=1.0, bad_panel_mask=bad,
            grad_mu_options=(; basis=:tri, tri_robust_max_depth=4, tri_robust_target_healthy=3))
        @test bfs_blocked[:, 1] ≈ zeros(3)
        @test norm(bfs_grown[:, 1]) > 0

        robust_manual = zeros(3, body.ncells)
        robust_public = zeros(3, body.ncells)
        ar_mask = pnl.panel_aspect_ratio_mask(body.nodes, body.cells; threshold=1.3)
        pnl._compute_mu_gradient_masked!(robust_manual, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info;
            scale=1.0,
            bad_panel_mask=any(ar_mask) ? ar_mask : nothing,
            nodes=body.nodes,
            grad_mu_options=(; basis=:tri, tri_robust_max_depth=4, tri_robust_target_healthy=3))
        pnl.compute_mu_gradient!(robust_public, body.controlpoints, body.normals,
            body.cells, body.neighbor, mu, te_info;
            scale=1.0,
            nodes=body.nodes,
            grad_mu_options=(; basis=:tri, tri_robust=true, tri_robust_ar_threshold=1.3,
                tri_robust_max_depth=4, tri_robust_target_healthy=3))
        @test robust_public ≈ robust_manual

        one_nodes = body.nodes[:, body.cells[:, 1]]
        one_cells = reshape([1, 2, 3], 3, 1)
        one_body = pnl.NonLiftingBody{pnl.ConstantSource}(one_nodes, one_cells;
            watertight=false,
            ensure_winding=false)
        pnl.calc_normals!(one_body)
        pnl.calc_controlpoints!(one_body)
        @test_throws ArgumentError pnl.compute_mu_gradient!(zeros(3, one_body.ncells),
            one_body.controlpoints, one_body.normals, one_body.cells,
            one_body.neighbor, [1.0], zeros(Int, 2, one_body.ncells);
            scale=1.0, nodes=one_body.nodes, grad_mu_options=(; basis=:quad))
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

        # Supplying nodes selects the paired-quad reconstruction. It remains
        # exact for constant vector gradients on highly anisotropic split quads.
        anisotropic_nodes = copy(nodes)
        anisotropic_nodes[1, :] .*= 100.0
        anisotropic = pnl.NonLiftingBody{pnl.ConstantSource}(anisotropic_nodes, cells;
            watertight=false, ensure_winding=false)
        pnl.calc_normals!(anisotropic); pnl.calc_controlpoints!(anisotropic)
        u_aniso = zeros(3, anisotropic.ncells)
        u_aniso[1, :] .= anisotropic.controlpoints[1, :] .+ 2 .* anisotropic.controlpoints[2, :]
        u_aniso[2, :] .= -3 .* anisotropic.controlpoints[1, :] .+ 0.5 .* anisotropic.controlpoints[2, :]
        grad_aniso = zeros(3, 3, anisotropic.ncells)
        pnl.compute_surface_velocity_gradient!(grad_aniso, u_aniso,
            anisotropic.controlpoints, anisotropic.normals, anisotropic.cells,
            anisotropic.neighbor, te_info; nodes=anisotropic.nodes)
        for i in 1:anisotropic.ncells
            @test isapprox(grad_aniso[1, :, i], [1.0, 2.0, 0.0]; atol=1e-8)
            @test isapprox(grad_aniso[2, :, i], [-3.0, 0.5, 0.0]; atol=1e-8)
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

        u_ref = zeros(3, body.ncells); u_ref[1, :] .= mu_ref
        u_perturbed = zeros(3, body.ncells); u_perturbed[1, :] .= mu_perturbed
        Gu_ref = zeros(3, 3, body.ncells); Gu_perturbed = similar(Gu_ref)
        pnl.compute_surface_velocity_gradient!(Gu_ref, u_ref, body.controlpoints,
            body.normals, body.cells, body.neighbor, te_info; nodes=body.nodes,
            grad_mu_options=(; basis=:tri))
        pnl.compute_surface_velocity_gradient!(Gu_perturbed, u_perturbed, body.controlpoints,
            body.normals, body.cells, body.neighbor, te_info; nodes=body.nodes,
            grad_mu_options=(; basis=:tri))
        @test isapprox(Gu_ref[:, :, upper_side], Gu_perturbed[:, :, upper_side]; atol=1e-12)
        @test isapprox(Gu_ref[:, :, .!upper_side], Gu_perturbed[:, :, .!upper_side]; atol=1e-12)
    end
end
