using Test
import FLOWPanel as pnl
import FLOWVPM
const ReadVTK = pnl.ReadVTK
const TOML = pnl.TOML

if !isdefined(@__MODULE__, :make_octa_source_body)
    include("test_helpers.jl")
end

function make_replay_surface_vorticity_body()
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
    body = pnl.NonLiftingBody{pnl.ConstantDoublet}(nodes, cells;
        watertight=false,
        ensure_winding=false)
    pnl.calc_normals!(body)
    pnl.calc_controlpoints!(body)
    body.strength[:, 1] .= vec(body.controlpoints[1, :] .+ 2 .* body.controlpoints[2, :])
    for i in 1:body.ncells
        body.velocity[:, i] .= (1.0 + 0.2i, -0.3 + 0.05i, 0.15)
    end
    return body
end

@testset "VTK replay" begin
    @testset "body VTK omits monitor-owned fields by default" begin
        body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_1TRI), copy(CELLS_1TRI);
            ensure_winding=false)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.strength[:, 1] .= 1.0

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)

        vtk = ReadVTK.VTKFile(joinpath(path, "run_body1", "run_body1.0.vtu"))
        keys_cell_data = keys(ReadVTK.get_cell_data(vtk))
        @test !("gauge pressure" in keys_cell_data)
        @test !("F" in keys_cell_data)
    end

    @testset "monitor VTK fields are configurable" begin
        body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_1TRI), copy(CELLS_1TRI);
            ensure_winding=false)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.strength[:, 1] .= 1.0
        body.velocity[1, :] .= 1.0
        pressure = pnl.PressureBernoulli(1.0)
        force = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())
        frames = pnl.ReferenceFrame(body)
        ctx = pnl.MonitorContext()
        pnl._run_monitor!(pressure, ctx, (body,), (nothing,), frames, [0.0, 0.0, 0.0], 0, 0.1)
        pnl._run_monitor!(force, ctx, (body,), (nothing,), frames, [0.0, 0.0, 0.0], 0, 0.1)

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0;
            monitors=(pressure, force), i_system=1, overwrite=true)

        vtk = ReadVTK.VTKFile(joinpath(path, "run_body1", "run_body1.0.vtu"))
        keys_cell_data = keys(ReadVTK.get_cell_data(vtk))
        @test "gauge pressure" in keys_cell_data
        @test "F" in keys_cell_data

        suppressed_pressure = pnl.PressureBernoulli(1.0; vtk_fields=())
        suppressed_force = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization(), vtk_fields=())
        ctx = pnl.MonitorContext()
        pnl._run_monitor!(suppressed_pressure, ctx, (body,), (nothing,), frames, [0.0, 0.0, 0.0], 0, 0.1)
        pnl._run_monitor!(suppressed_force, ctx, (body,), (nothing,), frames, [0.0, 0.0, 0.0], 0, 0.1)
        path2 = mktempdir()
        pnl.write_vtk(joinpath(path2, "run_body1"), body, 0, 0.0;
            monitors=(suppressed_pressure, suppressed_force), i_system=1, overwrite=true)
        vtk2 = ReadVTK.VTKFile(joinpath(path2, "run_body1", "run_body1.0.vtu"))
        keys_cell_data2 = keys(ReadVTK.get_cell_data(vtk2))
        @test !("gauge pressure" in keys_cell_data2)
        @test !("F" in keys_cell_data2)
    end

    @testset "monitor VTK field collisions are uniquely named" begin
        body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_1TRI), copy(CELLS_1TRI);
            ensure_winding=false)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.strength[:, 1] .= 1.0
        body.velocity[1, :] .= 1.0

        bernoulli = pnl.PressureBernoulli(1.0)
        laplace = pnl.PressureLaplace((body,), 1.0; reference_panel=1)
        frames = pnl.ReferenceFrame(body)
        ctx = pnl.MonitorContext()
        pnl._run_monitor!(bernoulli, ctx, (body,), (nothing,), frames, [0.0, 0.0, 0.0], 0, 0.1)
        laplace.p[1] .= bernoulli.pressure[1] .+ 1.0

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0;
            monitors=(bernoulli, laplace), i_system=1, overwrite=true)

        vtk = ReadVTK.VTKFile(joinpath(path, "run_body1", "run_body1.0.vtu"))
        cell_data = ReadVTK.get_cell_data(vtk)
        keys_cell_data = keys(cell_data)
        @test "gauge pressure" in keys_cell_data
        @test "gauge pressure (PressureLaplace #2)" in keys_cell_data
        @test ReadVTK.get_data(cell_data["gauge pressure"]) ≈ bernoulli.pressure[1]
        @test ReadVTK.get_data(cell_data["gauge pressure (PressureLaplace #2)"]) ≈ laplace.p[1]

        suppressed = pnl.PressureBernoulli(1.0; vtk_fields=())
        ctx = pnl.MonitorContext()
        pnl._run_monitor!(suppressed, ctx, (body,), (nothing,), frames, [0.0, 0.0, 0.0], 0, 0.1)
        path2 = mktempdir()
        pnl.write_vtk(joinpath(path2, "run_body1"), body, 0, 0.0;
            monitors=(suppressed, laplace), i_system=1, overwrite=true)

        vtk2 = ReadVTK.VTKFile(joinpath(path2, "run_body1", "run_body1.0.vtu"))
        keys_cell_data2 = keys(ReadVTK.get_cell_data(vtk2))
        @test "gauge pressure" in keys_cell_data2
        @test !("gauge pressure (PressureLaplace #2)" in keys_cell_data2)
    end

    @testset "force monitor VTK field collisions are uniquely named" begin
        body = make_replay_surface_vorticity_body()
        force = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())
        surface = pnl.SurfaceVorticityForce(body, 1, 1; normalization=pnl.NoNormalization())
        force.distributed_force = fill(1.0, 3, body.ncells)
        surface.distributed_force .= 2.0

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0;
            monitors=(force, surface), i_system=1, overwrite=true)

        vtk = ReadVTK.VTKFile(joinpath(path, "run_body1", "run_body1.0.vtu"))
        cell_data = ReadVTK.get_cell_data(vtk)
        keys_cell_data = keys(cell_data)
        @test "F" in keys_cell_data
        @test "F (SurfaceVorticityForce #2)" in keys_cell_data
        @test ReadVTK.get_data(cell_data["F"]) ≈ force.distributed_force
        @test ReadVTK.get_data(cell_data["F (SurfaceVorticityForce #2)"]) ≈ surface.distributed_force

        suppressed = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization(), vtk_fields=())
        suppressed.distributed_force = fill(3.0, 3, body.ncells)
        path2 = mktempdir()
        pnl.write_vtk(joinpath(path2, "run_body1"), body, 0, 0.0;
            monitors=(suppressed, surface), i_system=1, overwrite=true)

        vtk2 = ReadVTK.VTKFile(joinpath(path2, "run_body1", "run_body1.0.vtu"))
        keys_cell_data2 = keys(ReadVTK.get_cell_data(vtk2))
        @test "F" in keys_cell_data2
        @test !("F (SurfaceVorticityForce #2)" in keys_cell_data2)
    end

    @testset "recompute velocity from loaded strengths" begin
        body = make_octa_source_body()
        expected = deepcopy(body)
        expected.velocity .= 0
        pnl.influence!(expected, expected, pnl.DirectBackend(); velocity=true)

        body.velocity .= 0
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)

        result = pnl.replay(path, "run"; recompute=(:velocity,), backend=pnl.DirectBackend())
        @test result.systems[1].velocity ≈ expected.velocity
    end

    @testset "recompute velocity uses saved metadata freestream" begin
        body = make_octa_source_body()
        uinf = [2.0, -1.0, 0.5]
        expected = deepcopy(body)
        expected.velocity .= 0
        pnl.influence!(expected, expected, pnl.DirectBackend(); velocity=true)
        pnl.apply_freestream!(expected, uinf)

        body.velocity .= 0
        frames = pnl.ReferenceFrame(body)
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)
        pnl._append_metadata_step_toml(path, "run", frames, 0, 0.0; uinf)

        result = pnl.replay(path, "run"; recompute=(:velocity,), backend=pnl.DirectBackend(),
            Uinf=t -> [99.0, 99.0, 99.0])
        @test result.systems[1].velocity ≈ expected.velocity
    end

    @testset "force monitor uses loaded pressure" begin
        body = make_replay_surface_vorticity_body()
        pressure = pnl.PressureBernoulli(1.0)
        force_loaded = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())
        force_expected = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())
        frames = pnl.ReferenceFrame(body)
        ctx = pnl.MonitorContext()
        pressure.pressure = [collect(1.0:body.ncells)]
        pnl.monitor_register!(ctx, :P, 1, pressure.pressure[1])
        pnl._run_monitor!(force_expected, ctx, (body,), (nothing,), frames, [0.0, 0.0, 0.0], 0, 0.1)
        force_loaded.distributed_force = fill(99.0, 3, body.ncells)

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0;
            monitors=(pressure, force_loaded), i_system=1, overwrite=true)

        force = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())
        result = pnl.replay(path, "run"; monitors=(force,), recompute=(), backend=pnl.DirectBackend())
        @test result.monitors[1].distributed_force ≈ force_expected.distributed_force
        @test_throws ArgumentError pnl.replay(path, "run"; recompute=(:F,), backend=pnl.DirectBackend())
    end

    @testset "auxiliary recompute preserves loaded velocity" begin
        body = make_plate_vortex_body()
        pnl.calc_controlpoints!(body)
        for i in axes(body.velocity, 2)
            body.velocity[:, i] .= (1.0 + i, -0.5 * i, 0.25)
        end
        loaded_velocity = copy(body.velocity)
        frames = pnl.ReferenceFrame(body;
            ω_axis=pnl.FastMultipole.SVector{3}(0.0, 0.0, 1.0),
            ω=2.0)

        pnl._recompute_replay_fields!((body,), (nothing,), frames,
            pnl.FastMultipole.SVector{3,Float64}(0.0, 0.0, 0.0),
            Set([:induced_vorticity]), pnl.DirectBackend(), pnl.DirectBackend())

        @test body.velocity ≈ loaded_velocity
        @test any(!iszero, body.velocity_kinematic)
        @test body.angular_velocity ≈ [0.0, 0.0, 2.0]
    end

    @testset "monitor factory sees reconstructed bodies" begin
        body = make_plate_vortex_body()
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)

        called = Ref(false)
        factory = (systems, wakes, frames, t_range) -> begin
            called[] = true
            @test length(systems) == 1
            @test systems[1].ncells == body.ncells
            (pnl.KuttaJoukowskiForce(systems[1], length(t_range), 1; backend=pnl.DirectBackend()),)
        end

        result = pnl.replay(path, "run"; monitor_factory=factory, backend=pnl.DirectBackend())
        @test called[]
        @test length(result.monitors) == 1
    end

    @testset "surface vorticity replay comparison monitors" begin
        body = make_replay_surface_vorticity_body()
        body.kerneloffset_panel = 1e-9
        body.kerneloffset_targets = 1e-2
        body.kerneloffset = body.kerneloffset_panel
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)

        kj_only = Ref{Any}()
        surface_only = Ref{Any}()
        factory_diagnostic = (systems, wakes, frames, t_range) -> begin
            surface_only[] = pnl.SurfaceVorticityForce(systems[1], length(t_range), 1;
                normalization=pnl.NoNormalization())
            kj_only[] = pnl.KuttaJoukowskiForce(systems[1], length(t_range), 1;
                backend=pnl.DirectBackend(),
                normalization=pnl.NoNormalization())
            (surface_only[], kj_only[])
        end
        pnl.replay(path, "run"; monitor_factory=factory_diagnostic, backend=pnl.DirectBackend())

        kj_after_aux = Ref{Any}()
        surface_after_aux = Ref{Any}()
        bernoulli_force = Ref{Any}()
        laplace_lamb_force = Ref{Any}()
        factory_comparison = (systems, wakes, frames, t_range) -> begin
            bernoulli = pnl.PressureBernoulli(1.0; correct_kuttacondition=false)
            bernoulli_force[] = pnl.ForceMonitor(length(t_range), 1;
                normalization=pnl.NoNormalization(),
                correct_kuttacondition=false)
            surface_after_aux[] = pnl.SurfaceVorticityForce(systems[1], length(t_range), 1;
                normalization=pnl.NoNormalization())
            kj_after_aux[] = pnl.KuttaJoukowskiForce(systems[1], length(t_range), 1;
                backend=pnl.DirectBackend(),
                normalization=pnl.NoNormalization())
            laplace_lamb = pnl.PressureLaplace(systems[1], 1.0;
                acceleration_form=:lamb_vector,
                reference_panel=1)
            laplace_lamb_force[] = pnl.ForceMonitor(length(t_range), 1;
                normalization=pnl.NoNormalization(),
                correct_kuttacondition=false)
            (bernoulli, bernoulli_force[], surface_after_aux[], kj_after_aux[],
                laplace_lamb, laplace_lamb_force[])
        end
        pnl.replay(path, "run"; monitor_factory=factory_comparison, backend=pnl.DirectBackend())

        @test kj_only[].force[:, 1] ≈ kj_after_aux[].force[:, 1]
        @test surface_only[].force[:, 1] ≈ surface_after_aux[].force[:, 1]
        @test all(isfinite, surface_after_aux[].force[:, 1])
        @test all(isfinite, kj_after_aux[].force[:, 1])
        @test all(isfinite, bernoulli_force[].force[:, 1])
        @test all(isfinite, laplace_lamb_force[].force[:, 1])
        @test !isapprox(surface_after_aux[].force[:, 1], zeros(3); atol=1e-12)
        @test !isapprox(bernoulli_force[].force[:, 1], zeros(3); atol=1e-12)
    end

    @testset "missing rigid metadata errors and reconstruct can fill" begin
        body = make_plate_vortex_body()
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)
        open(joinpath(path, "run.replay.toml"), "w") do io
            println(io, "[[body]]")
            println(io, "i = 1")
            println(io, "kind = \"RigidWakeBody\"")
            println(io, "strength_names = [\"gamma\"]")
        end

        @test_throws ArgumentError pnl.replay(path, "run")

        reconstruct = (p, name, metadata) -> begin
            ((body,), (nothing,), pnl.ReferenceFrame(body))
        end
        result = pnl.replay(path, "run"; reconstruct)
        @test result.systems[1] isa pnl.RigidWakeBody
    end

    @testset "metadata replay reconstructs particle wake settings" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            nwakerows=2,
            max_particles=128,
            method_trailing=pnl.NoShed(),
            method_unsteady=pnl.SigmaOverlap(0.25, 4.0),
            particle_maintenance=pnl.ParticleMaintenance((
                pnl.MinGamma(1e-3),
                pnl.GlobalCylinder([0.0, -1.0, -1.0], [2.0, 0.0, 0.0], 1.5),
                pnl.MergeParticles(every=2, r=0.4, r_hash=0.3, sigma_relative=false, max_sigma_ratio=1.7, skip_static=false),
            )),
            viscous=FLOWVPM.CoreSpreading(1.5e-5, 0.01, FLOWVPM.zeta_fmm; beta=1.5),
            SFS=FLOWVPM.SFS_Cd_twolevel_nobackscatter,
            # non-default relaxation so we can verify it round-trips (the
            # PanelParticleWake default is correctedpedrizzetti)
            relaxation=FLOWVPM.relaxation_pedrizzetti,
        )
        frames = pnl.ReferenceFrame(body)
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)
        pnl.write_vtk(joinpath(path, "run_wake1"), wake, 0, 0.0; overwrite=true)
        pnl._write_metadata_toml(path, "run", (body,), (wake,), frames, [0.0, 0.1],
            (pnl.Backslash(body),), pnl.DirectBackend(), pnl.DirectBackend(), pnl.DirectBackend(), ();
            solver_options=(; particle_relax=false, body_on_wake=true, bound_strength_rlx=0.5))
        pnl._append_metadata_step_toml(path, "run", frames, 0, 0.0)

        # the run-affecting simulate! toggles are recorded under [simulation]
        meta = TOML.parsefile(joinpath(path, "run.metadata.toml"))
        @test meta["simulation"]["particle_relax"] == false
        @test meta["simulation"]["body_on_wake"] == true
        @test meta["simulation"]["bound_strength_rlx"] == 0.5
        # the relaxation scheme is recorded under the wake's pfield_optargs
        @test meta["wake"][1]["pfield_optargs"]["relaxation"]["type"] == "FLOWVPM.relax_pedrizzetti"

        result = pnl.replay(path, "run"; steps=0, recompute=())
        @test result.wakes[1] isa pnl.PanelParticleWake
        @test result.wakes[1].method_trailing isa pnl.NoShed
        @test result.wakes[1].method_unsteady isa pnl.SigmaOverlap
        @test result.wakes[1].particle_maintenance.trim_policies[2] isa pnl.GlobalCylinder
        @test result.wakes[1].pfield.SFS === FLOWVPM.SFS_Cd_twolevel_nobackscatter
        # relaxation scheme survives the metadata round-trip (not silently reset)
        @test result.wakes[1].pfield.relaxation.relax === FLOWVPM.relax_pedrizzetti
    end

    @testset "particle wake explicit field round trips" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; nwakerows=2, max_particles=8,
            viscous=FLOWVPM.CoreSpreading(1.5e-5, 0.01, FLOWVPM.zeta_fmm; beta=1.5),
            SFS=FLOWVPM.SFS_Cd_twolevel_nobackscatter)
        wake.panel_wake.nwakes[] = 1
        wake.pfield.np = 3
        for i in 1:wake.pfield.np
            wake.pfield.particles[FLOWVPM.X_INDEX, i] .= (10.0 + i, 20.0 + i, 30.0 + i)
            wake.pfield.particles[FLOWVPM.GAMMA_INDEX, i] .= (1.0 + i, 2.0 + i, 3.0 + i)
            wake.pfield.particles[FLOWVPM.SIGMA_INDEX, i] = 2.0 * i
            wake.pfield.particles[FLOWVPM.VOL_INDEX, i] = 3.0 * i
            wake.pfield.particles[FLOWVPM.CIRCULATION_INDEX, i] = 4.0 * i
            wake.pfield.particles[FLOWVPM.U_INDEX, i] .= (5.0 + i, 6.0 + i, 7.0 + i)
            wake.pfield.particles[FLOWVPM.VORTICITY_INDEX, i] .= (8.0 + i, 9.0 + i, 10.0 + i)
            wake.pfield.particles[FLOWVPM.C_INDEX, i] .= (11.0 + i, 12.0 + i, 13.0 + i)
            wake.pfield.particles[FLOWVPM.J_INDEX, i] .= collect(14.0 + 9 * (i - 1):14.0 + 9 * i - 1)
        end
        expected = copy(view(wake.pfield.particles, :, 1:wake.pfield.np))

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_wake1"), wake, 0, 0.0; overwrite=true)

        vtk = ReadVTK.VTKFile(joinpath(path, "run_wake1_particles", "run_wake1_particles.0.vtp"))
        point_data = ReadVTK.get_point_data(vtk)
        @test !("particle_state" in keys(point_data))
        @test "C" in keys(point_data)

        loaded = pnl.PanelParticleWake(body; nwakerows=2, max_particles=8,
            viscous=FLOWVPM.CoreSpreading(1.5e-5, 0.01, FLOWVPM.zeta_fmm; beta=1.5),
            SFS=FLOWVPM.SFS_Cd_twolevel_nobackscatter)
        loaded.pfield.np = 5
        loaded.pfield.particles[:, 1:5] .= -1
        pnl._load_panel_particle_wake_vtk!(loaded, path, "run_wake1", 0)

        @test loaded.pfield.np == wake.pfield.np
        @test view(loaded.pfield.particles, :, 1:loaded.pfield.np) ≈ expected
        @test all(iszero, loaded.pfield.particles[:, loaded.pfield.np + 1:end])
        @test loaded.pfield.viscous isa FLOWVPM.CoreSpreading
        @test loaded.pfield.SFS === FLOWVPM.SFS_Cd_twolevel_nobackscatter
    end

    @testset "unified metadata without requested frame state errors" begin
        body = make_plate_vortex_body()
        frames = pnl.ReferenceFrame(body; ω=2.0)
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)
        pnl._write_metadata_toml(path, "run", (body,), (nothing,), frames, [0.0],
            (pnl.Backslash(body),), pnl.DirectBackend(), pnl.DirectBackend(), pnl.DirectBackend(), ())

        @test_throws ArgumentError pnl.replay(path, "run"; steps=0, recompute=())
    end
end
