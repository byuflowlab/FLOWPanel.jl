using Test
import FLOWPanel as pnl
import FLOWVPM
const TOML = pnl.TOML
import FastMultipole
using LinearAlgebra: norm
using StaticArrays: SVector

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

@testset verbose=true "Free Wakes" begin
    @testset "SigmaOverlap chooses particle count from fixed sigma" begin
        pfield = FLOWVPM.ParticleField(10, Float64;
            fmm=FLOWVPM.FMM(autotune_reg_error=false))
        r1 = SVector(0.0, 0.0, 0.0)
        r2 = SVector(1.0, 0.0, 0.0)
        sigma = 0.2
        target_overlap = 0.7

        pnl._shed_particles!(pfield, r1, r2, 1.0, pnl.SigmaOverlap(sigma, target_overlap))

        @test pfield.np == 4
        @test all(pfield.particles[FLOWVPM.SIGMA_INDEX, 1:pfield.np] .== sigma)
        @test sigma * pfield.np / norm(r2 - r1) >= target_overlap
        @test pfield.particles[1, 1:pfield.np] ≈ [0.125, 0.375, 0.625, 0.875]
    end

    @testset "PanelWake first row shedding velocity option" begin
        body = make_plate_vortex_body()
        dt = 0.25
        uinf = [1.0, 2.0, 3.0]
        induced_first = [10.0, 20.0, 30.0]
        induced_downstream = [100.0, 200.0, 300.0]

        wake_total = pnl.PanelWake(body; nwakerows=2)
        wake_total.nwakes[] = 1
        pnl.apply_freestream!(wake_total, uinf)
        wake_total.velocity[1][:, 1, :] .+= induced_first
        wake_total.velocity[1][:, 2, :] .+= induced_downstream
        pnl.propagate!(wake_total, dt)

        @test wake_total.nodes[1][:, 1, 1] ≈ dt .* (uinf .+ induced_first)
        @test wake_total.nodes[1][:, 2, 1] ≈ dt .* (uinf .+ induced_downstream)

        wake_freestream = pnl.PanelWake(body; nwakerows=2,
            shed_with_induced_velocity=false)
        wake_freestream.nwakes[] = 1
        pnl.apply_freestream!(wake_freestream, uinf)
        wake_freestream.velocity[1][:, 1, :] .+= induced_first
        wake_freestream.velocity[1][:, 2, :] .+= induced_downstream
        pnl.propagate!(wake_freestream, dt)

        @test wake_freestream.nodes[1][:, 1, 1] ≈ dt .* uinf
        @test wake_freestream.nodes[1][:, 2, 1] ≈ dt .* (uinf .+ induced_downstream)
    end

    @testset "PanelWake optional final filament" begin
        body = make_plate_vortex_body()
        default_wake = pnl.PanelWake(body; nwakerows=2)
        panel_only = pnl.PanelWake(body; nwakerows=2,
            include_final_filament=false)
        @test length(pnl.get_sources(default_wake)) == 2
        @test pnl.get_sources(panel_only) == (panel_only,)
        @test default_wake.include_final_filament
        @test !panel_only.include_final_filament
    end

    @testset "PanelParticleWake forwards shedding velocity option" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; shed_with_induced_velocity=false)

        @test !wake.panel_wake.shed_with_induced_velocity
    end

    @testset "PanelParticleWake stores particle kernel offset" begin
        body = make_plate_vortex_body()
        default_wake = pnl.PanelParticleWake(body)
        offset_wake = pnl.PanelParticleWake(body; particle_kerneloffset=0.25)

        @test isnan(default_wake.particle_kerneloffset)
        @test offset_wake.particle_kerneloffset == 0.25
        @test body.kerneloffset_targets == 0.25
    end

    @testset "PanelParticleWake forwards SFS model" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; SFS=FLOWVPM.SFS_Cd_twolevel_nobackscatter)

        @test wake.pfield.SFS === FLOWVPM.SFS_Cd_twolevel_nobackscatter
        @test FLOWVPM.isLES(wake.pfield)
    end

    @testset "PanelParticleWake forwards viscous model" begin
        body = make_plate_vortex_body()
        viscous = FLOWVPM.CoreSpreading(1.5e-5, 0.01, FLOWVPM.zeta_fmm)
        wake = pnl.PanelParticleWake(body; viscous)

        @test wake.pfield.viscous === viscous
        @test FLOWVPM.iscorespreading(wake.pfield.viscous)
    end

    @testset "PanelWake final filament strength option" begin
        body = make_plate_vortex_body()

        function buffered_filament_strength(; unsteady_filament)
            wake = pnl.PanelWake(body; nwakerows=2, unsteady_filament)
            wake.nwakes[] = 2
            wake.overflowed[] = true
            wake.strength[1][1, 2, 1] = 11.0
            wake.strength[1][1, 3, 1] = 22.0
            wake.nodes[1][:, 3, 1] .= [0.0, 0.0, 0.0]
            wake.nodes[1][:, 3, 2] .= [1.0, 0.0, 0.0]

            buffer = zeros(12, 1)
            FastMultipole.source_system_to_buffer!(buffer, 1, pnl.FilamentWrapper(wake), 1)
            return buffer[5, 1]
        end

        @test buffered_filament_strength(unsteady_filament=true) == -22.0
        @test buffered_filament_strength(unsteady_filament=false) == -11.0
    end

    @testset "PanelWake steady final filament cancels last wake row influence" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelWake(body; nwakerows=1, core_size=1e-3, unsteady_filament=false)
        wake.nwakes[] = 1
        wake.overflowed[] = true

        Γ = 2.5
        wake.strength[1][1, 1, 1] = Γ
        wake.strength[1][1, 2, 1] = -3.0 # ignored when unsteady_filament=false

        v1 = SVector(0.0, 0.0, 0.0)
        v2 = SVector(1.0, 0.0, 0.0)
        v3 = SVector(1.0, 1.0, 0.0)
        v4 = SVector(0.0, 1.0, 0.0)
        wake.nodes[1][:, 1, 1] .= v1
        wake.nodes[1][:, 2, 1] .= v2
        wake.nodes[1][:, 2, 2] .= v3
        wake.nodes[1][:, 1, 2] .= v4

        switch = FastMultipole.DerivativesSwitch(false, true, true)
        panel_buffer = zeros(FastMultipole.data_per_body(wake), 1)
        FastMultipole.source_system_to_buffer!(panel_buffer, 1, wake, 1)

        filament = pnl.FilamentWrapper(wake)
        filament_buffer = zeros(FastMultipole.data_per_body(filament), 1)
        FastMultipole.source_system_to_buffer!(filament_buffer, 1, filament, 1)

        function filament_velocity(target, a, b, strength)
            return pnl._bound_vortex_velocity(target - a, target - b, true, wake.core_size) * strength
        end

        function filament_gradient(target, a, b, strength)
            return pnl._bound_vortex_gradient(a - target, b - target, true, wake.core_size) * strength
        end

        for target in (
            SVector(0.25, 0.35, 0.7),
            SVector(-0.4, 0.8, 0.5),
            SVector(1.0 + 1e-6, 0.45, 1e-6),
        )
            _, panel_velocity, panel_gradient = pnl.induced(target, wake, panel_buffer, 1, switch)
            final_velocity = filament_velocity(target, v2, v3, filament_buffer[5, 1])
            final_gradient = filament_gradient(target, v2, v3, filament_buffer[5, 1])
            combined_velocity = panel_velocity + final_velocity
            combined_gradient = panel_gradient + final_gradient

            expected_velocity =
                filament_velocity(target, v1, v2, Γ) +
                filament_velocity(target, v3, v4, Γ) +
                filament_velocity(target, v4, v1, Γ)
            expected_gradient =
                filament_gradient(target, v1, v2, Γ) +
                filament_gradient(target, v3, v4, Γ) +
                filament_gradient(target, v4, v1, Γ)

            @test combined_velocity ≈ expected_velocity atol=1e-12 rtol=1e-12
            @test combined_gradient ≈ expected_gradient atol=1e-10 rtol=1e-10
        end
    end

    @testset "PanelParticleWake forwards unsteady filament option" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; unsteady_filament=false)

        @test !wake.panel_wake.unsteady_filament
    end

    @testset "Kernel offset override regularizes body-to-particle influence" begin
        body = make_plate_vortex_body()
        body.kerneloffset = 1e-8
        wake = pnl.PanelParticleWake(body; particle_kerneloffset=0.1)
        FLOWVPM.add_particle(wake.pfield, [0.5, 1e-4, 1e-4], [0.0, 0.0, 0.0], 0.01)

        function particle_speed(kerneloffset)
            FLOWVPM._reset_particles(wake.pfield)
            old = body.kerneloffset
            body.kerneloffset = kerneloffset
            try
                pnl.influence!((wake.pfield,), (body,), pnl.DirectBackend();
                    velocity=true, velocity_gradient=(false,))
                return norm(view(wake.pfield.particles, FLOWVPM.U_INDEX, 1))
            finally
                body.kerneloffset = old
            end
        end

        baseline = particle_speed(body.kerneloffset_panel)
        regularized = particle_speed(body.kerneloffset_targets)

        @test regularized < baseline
    end

    @testset "Self panel conditioning leaves particle targets on target offset" begin
        body = make_plate_vortex_body()
        body.kerneloffset_panel = 1e-8
        body.kerneloffset_targets = 0.1
        body.kerneloffset = body.kerneloffset_targets
        wake = pnl.PanelParticleWake(body; particle_kerneloffset=body.kerneloffset_targets)
        FLOWVPM.add_particle(wake.pfield, [0.5, 1e-4, 1e-4], [0.0, 0.0, 0.0], 0.01)

        FLOWVPM._reset_particles(wake.pfield)
        body.velocity .= 0
        pnl.influence!((body, wake.pfield), (body,), pnl.DirectBackend();
            velocity=true,
            velocity_gradient=(false, false),
            direct_conditioning=pnl._self_panel_kerneloffset_conditioning())
        conditioned_velocity = copy(view(wake.pfield.particles, FLOWVPM.U_INDEX, 1))

        FLOWVPM._reset_particles(wake.pfield)
        body.kerneloffset = body.kerneloffset_targets
        pnl.influence!((wake.pfield,), (body,), pnl.DirectBackend();
            velocity=true, velocity_gradient=(false,))
        target_offset_velocity = copy(view(wake.pfield.particles, FLOWVPM.U_INDEX, 1))

        @test body.kerneloffset == body.kerneloffset_targets
        @test conditioned_velocity ≈ target_offset_velocity atol=1e-12 rtol=1e-12
    end

    @testset "ParticleMaintenance separates mixed policies" begin
        maintenance = pnl.ParticleMaintenance((
            pnl.MinGamma(1e-3),
            pnl.MergeParticles(every=2, r_hash=0.25),
            pnl.GlobalBox([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]),
        ))

        @test maintenance.trim_policies isa Tuple{<:pnl.MinGamma,<:pnl.GlobalBox}
        @test maintenance.functional_policies isa Tuple{<:pnl.MergeParticles}
        @test_throws ArgumentError pnl.ParticleMaintenance((pnl.MinGamma(1e-3), :not_a_policy))
    end

    @testset "PanelParticleWake particle maintenance trims by strength" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.MinGamma(1e-3),)))

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [1e-4, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [3.0, 0.0, 0.0], [5e-4, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 2
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 1)) == [2e-3, 0.0, 0.0]
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 2)) == [1e-3, 0.0, 0.0]
    end

    @testset "PanelParticleWake particle maintenance trims by strength range" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.MinGamma(1e-3), pnl.MaxGamma(3e-3))))

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [5e-4, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [4e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 1
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 1)) == [2e-3, 0.0, 0.0]
    end

    @testset "PanelParticleWake particle maintenance trims by global box" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.GlobalBox([0.0, -1.0, -1.0], [2.0, 1.0, 1.0]),)))

        FLOWVPM.add_particle(wake.pfield, [-1.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [3e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [3.0, 0.0, 0.0], [4e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 2
        @test sort([FLOWVPM.get_X(wake.pfield, i)[1] for i in 1:wake.pfield.np]) == [0.0, 2.0]
    end

    @testset "PanelParticleWake particle maintenance trims by global cylinder" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.GlobalCylinder([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0),)))

        FLOWVPM.add_particle(wake.pfield, [-0.1, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.5, 0.0], [3e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 1.1, 0.0], [4e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 1.0, 0.0], [5e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.1, 0.0, 0.0], [6e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        retained = [SVector{3}(FLOWVPM.get_X(wake.pfield, i)) for i in 1:wake.pfield.np]
        @test wake.pfield.np == 3
        @test SVector(0.0, 0.0, 0.0) in retained
        @test SVector(1.0, 0.5, 0.0) in retained
        @test SVector(2.0, 1.0, 0.0) in retained
        @test_throws ArgumentError pnl.GlobalCylinder([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0)
        @test_throws ArgumentError pnl.GlobalCylinder([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], -1.0)
    end

    @testset "PanelParticleWake particle maintenance trims by frame box" begin
        body = make_plate_vortex_body()
        frames = pnl.ReferenceFrame(body; origin=SVector(10.0, 0.0, 0.0))
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.FrameBox(1, [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]),)))

        FLOWVPM.add_particle(wake.pfield, [10.5, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [12.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)

        @test_throws ArgumentError pnl.propagate!(wake, 0.0; relax=false, step=1)

        pnl.propagate!(wake, 0.0; relax=false, step=1, frames)

        @test wake.pfield.np == 1
        @test collect(FLOWVPM.get_X(wake.pfield, 1)) == [10.5, 0.0, 0.0]
    end

    @testset "PanelParticleWake particle maintenance merges particles" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body;
            particle_maintenance=pnl.ParticleMaintenance((pnl.MergeParticles(every=1, r=0.5, r_hash=0.5),)))

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [0.1, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 1
    end

    @testset "RelaxationPlaneFilter validation and frame refresh" begin
        body = make_plate_vortex_body()

        @test_throws ArgumentError pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
        @test_throws ArgumentError pnl.RelaxationPlaneFilter([0.0, 0.0], [1.0, 0.0, 0.0])
        @test_throws ArgumentError pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [Inf, 0.0, 0.0])

        global_filter = pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        global_wake = pnl.PanelParticleWake(body;
            relaxation=FLOWVPM.Relaxation(FLOWVPM.relax_correctedpedrizzetti, 1, 0.3, global_filter))
        pnl.propagate!(global_wake, 0.0; relax=false, step=1)

        frame_filter = pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]; i_frame=1)
        frame_wake = pnl.PanelParticleWake(body;
            relaxation=FLOWVPM.Relaxation(FLOWVPM.relax_correctedpedrizzetti, 1, 0.3, frame_filter))
        @test_throws ArgumentError pnl.propagate!(frame_wake, 0.0; relax=false, step=1)

        bad_frame_filter = pnl.RelaxationPlaneFilter([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]; i_frame=99)
        bad_frame_wake = pnl.PanelParticleWake(body;
            relaxation=FLOWVPM.Relaxation(FLOWVPM.relax_correctedpedrizzetti, 1, 0.3, bad_frame_filter))
        frames = pnl.ReferenceFrame(body)
        @test_throws ArgumentError pnl.propagate!(bad_frame_wake, 0.0; relax=false, step=1, frames)
    end

    @testset "Unified wake metadata round trip" begin
        body = make_plate_vortex_body()
        relaxation = pnl.plane_filtered_relaxation(
            FLOWVPM.relaxation_correctedpedrizzetti,
            SVector(0.0, 0.0, -2.0),
            SVector(0.0, 0.0, -2.0);
            i_frame=1)
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
            relaxation,
        )
        frames = pnl.ReferenceFrame(body)
        path = mktempdir()
        pnl._write_metadata_toml(path, "run", (body,), (wake,), frames, [0.0, 0.1],
            (pnl.Backslash(body),), pnl.DirectBackend(), pnl.DirectBackend(), pnl.DirectBackend(), ())

        metadata = TOML.parsefile(joinpath(path, "run.metadata.toml"))
        @test metadata["wake"][1]["method_trailing"]["type"] == "NoShed"
        @test metadata["wake"][1]["pfield_optargs"]["viscous"]["type"] == "FLOWVPM.CoreSpreading"
        @test metadata["wake"][1]["pfield_optargs"]["SFS"]["type"] == "FLOWVPM.SFS_Cd_twolevel_nobackscatter"
        @test metadata["wake"][1]["particle_maintenance"]["trim_policies"][2]["type"] == "GlobalCylinder"
        @test metadata["wake"][1]["particle_maintenance"]["trim_policies"][2]["origin"] == [0.0, -1.0, -1.0]
        @test metadata["wake"][1]["particle_maintenance"]["trim_policies"][2]["extrude"] == [2.0, 0.0, 0.0]
        @test metadata["wake"][1]["particle_maintenance"]["trim_policies"][2]["radius"] == 1.5
        @test metadata["wake"][1]["pfield_optargs"]["relaxation"]["filter"]["type"] == "RelaxationPlaneFilter"
        @test metadata["wake"][1]["pfield_optargs"]["relaxation"]["filter"]["point"] == [0.0, 0.0, -2.0]
        @test metadata["wake"][1]["pfield_optargs"]["relaxation"]["filter"]["normal"] == [0.0, 0.0, -1.0]
        @test metadata["wake"][1]["pfield_optargs"]["relaxation"]["filter"]["i_frame"] == 1

        reconstructed = pnl._construct_wakes_from_manifest((body,), metadata)
        @test reconstructed[1] isa pnl.PanelParticleWake
        @test reconstructed[1].method_trailing isa pnl.NoShed
        @test reconstructed[1].method_unsteady isa pnl.SigmaOverlap
        @test reconstructed[1].particle_maintenance.trim_policies[2] isa pnl.GlobalCylinder
        @test reconstructed[1].particle_maintenance.functional_policies[1] isa pnl.MergeParticles
        @test reconstructed[1].pfield.SFS === FLOWVPM.SFS_Cd_twolevel_nobackscatter
        reconstructed_filter = reconstructed[1].pfield.relaxation.filter
        @test reconstructed_filter isa pnl.RelaxationPlaneFilter
        @test reconstructed_filter.point_local == SVector(0.0, 0.0, -2.0)
        @test reconstructed_filter.normal_local == SVector(0.0, 0.0, -1.0)
        @test reconstructed_filter.i_frame == 1

        step1_frames = pnl.ReferenceFrame(body; ω_axis=SVector(0.0, 1.0, 0.0), ω=1.0)
        step2_frames = pnl.ReferenceFrame(body; ω_axis=SVector(0.0, 1.0, 0.0), ω=2.0)
        replacement_frames = pnl.ReferenceFrame(body; ω_axis=SVector(0.0, 1.0, 0.0), ω=3.0)

        pnl._append_metadata_step_toml(path, "run", step1_frames, 1, 0.1)
        pnl._append_metadata_step_toml(path, "run", step2_frames, 2, 0.2)
        pnl._append_metadata_step_toml(path, "run", replacement_frames, 2, 0.25)

        metadata = TOML.parsefile(joinpath(path, "run.metadata.toml"))
        @test length(metadata["step"]) == 2
        @test pnl._metadata_step_frames(metadata, 1)[1]["omega"] == 1.0
        @test pnl._metadata_step_frames(metadata, 2)[1]["omega"] == 3.0
        @test pnl._metadata_step_frames(metadata, 3) === nothing
    end
end
