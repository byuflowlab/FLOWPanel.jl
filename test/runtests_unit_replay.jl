using Test
import FLOWPanel as pnl
import FLOWVPM
import TOML

if !isdefined(@__MODULE__, :make_octa_source_body)
    include("test_helpers.jl")
end

@testset "VTK replay" begin
    @testset "loaded pressure drives ForceMonitor" begin
        body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_1TRI), copy(CELLS_1TRI);
            ensure_winding=false)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.strength[:, 1] .= 1.0
        body.P .= 2.0
        pnl.calcfield_F!(body)

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)

        monitor = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())
        result = pnl.replay(path, "run"; monitors=(monitor,), recompute=())

        @test result.steps == [0]
        @test result.systems[1].P == body.P
        @test monitor.force[:, 1] ≈ pnl.calcfield_Ftot(body)
    end

    @testset "recompute force overwrites saved F" begin
        body = pnl.NonLiftingBody{pnl.ConstantSource}(copy(NODES_1TRI), copy(CELLS_1TRI);
            ensure_winding=false)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.strength[:, 1] .= 1.0
        body.P .= 3.0
        body.F .= 99.0

        expected = zeros(size(body.F))
        pnl.calcfield_F!(expected, body, pnl.calc_areas(body), body.normals, body.P)

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)

        result = pnl.replay(path, "run"; recompute=(:F,))
        @test result.systems[1].F ≈ expected
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
        body.P .= 1.0
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
                pnl.MergeParticles(every=2, r=0.4, r_hash=0.3, sigma_relative=false, max_sigma_ratio=1.7, skip_static=false),
            )),
            viscous=FLOWVPM.CoreSpreading(1.5e-5, 0.01, FLOWVPM.zeta_fmm; beta=1.5),
            SFS=FLOWVPM.SFS_Cd_twolevel_nobackscatter,
        )
        frames = pnl.ReferenceFrame(body)
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "run_body1"), body, 0, 0.0; overwrite=true)
        pnl.write_vtk(joinpath(path, "run_wake1"), wake, 0, 0.0; overwrite=true)
        pnl._write_metadata_toml(path, "run", (body,), (wake,), frames, [0.0, 0.1],
            (pnl.Backslash(body),), pnl.DirectBackend(), pnl.DirectBackend(), pnl.DirectBackend(), ())

        result = pnl.replay(path, "run"; steps=0, recompute=())
        @test result.wakes[1] isa pnl.PanelParticleWake
        @test result.wakes[1].method_trailing isa pnl.NoShed
        @test result.wakes[1].method_unsteady isa pnl.SigmaOverlap
        @test result.wakes[1].pfield.SFS === FLOWVPM.SFS_Cd_twolevel_nobackscatter
    end
end
