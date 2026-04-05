using Test
import FLOWPanel as pnl
import GeometricTools as gt

@testset verbose=true "RigidWakeBody" begin
    @testset "construction" begin
        body = make_plate_vortex_body()
        @test body.nnodes == 4
        @test body.ncells == 2
        @test length(body.shedding) == 1
        @test size(body.Das[1]) == (3, size(body.shedding[1], 2) + 1)
        @test size(body.strength) == (body.ncells, 1)
    end

    @testset "shedding edge consistency" begin
        body = make_plate_vortex_body()
        @test size(body.shedding[1], 1) == 6
        @test size(body.shedding_full) == (6, body.ncells)
        te_panels = findall(!=(-1), vec(body.shedding_full[3, :]))
        @test !isempty(te_panels)
        @test all(body.shedding_full[4, te_panels] .> 0)
    end

    @testset "extra_reset!" begin
        body = make_plate_vortex_body()
        for vte in body.velocity_te
            vte .= 3.0
        end
        pnl.reset!(body)
        @test all(all(vte .== 0) for vte in body.velocity_te)
    end
end
