using Test
import FLOWPanel as pnl
import FLOWVPM

@testset verbose=true "Free Wakes" begin
    @testset "PanelParticleWake weak-particle trimming" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; gamma_trim_threshold=1e-3)

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [1e-4, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.0, 0.0], [1e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [2.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [3.0, 0.0, 0.0], [5e-4, 0.0, 0.0], 1.0)

        pnl.trim_weak_particles!(wake)

        @test wake.pfield.np == 2
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 1)) == [2e-3, 0.0, 0.0]
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 2)) == [1e-3, 0.0, 0.0]
    end

    @testset "PanelParticleWake propagate! trims weak particles" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; gamma_trim_threshold=1e-3)

        FLOWVPM.add_particle(wake.pfield, [0.0, 0.0, 0.0], [5e-4, 0.0, 0.0], 1.0)
        FLOWVPM.add_particle(wake.pfield, [1.0, 0.0, 0.0], [2e-3, 0.0, 0.0], 1.0)

        pnl.propagate!(wake, 0.0; relax=false, step=1)

        @test wake.pfield.np == 1
        @test collect(FLOWVPM.get_Gamma(wake.pfield, 1)) == [2e-3, 0.0, 0.0]
    end

    @testset "PanelParticleWake stores merge hash radius" begin
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; merge_r=0.5, merge_r_hash=0.25)

        @test wake.merge_r == 0.5
        @test wake.merge_r_hash == 0.25
    end
end
