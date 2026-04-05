using Test
import FLOWPanel as pnl
import GeometricTools as gt
using Statistics: mean

@testset verbose=true "Analytical Validation" begin
    @testset "Sphere in uniform flow - Neumann" begin
        body = make_sphere_source_body(ntheta=12, nphi=24)
        solve_source_body!(body)

        cps = body.controlpoints
        cpvals = body.Cp
        xcoords = view(cps, 1, :)
        ycoords = view(cps, 2, :)
        zcoords = view(cps, 3, :)

        i_stag = argmax(xcoords)
        @test isapprox(cpvals[i_stag], 1.0; atol=0.05)

        eq_mask = abs.(xcoords) .< 0.05
        @test any(eq_mask)
        cp_equator = mean(cpvals[eq_mask])
        @test isapprox(cp_equator, -1.25; atol=0.2)

        q = 0.5 * pi
        ftot = pnl.calcfield_Ftot!(zeros(3), body, body.F)
        @test abs(ftot[1]) / q < 0.05
        @test abs(ftot[3]) / q < 0.05

        pos_y = cpvals[ycoords .> 0]
        neg_y = cpvals[ycoords .< 0]
        @test !isempty(pos_y)
        @test !isempty(neg_y)
        @test abs(mean(pos_y) - mean(neg_y)) < 0.02

        pos_z = cpvals[zcoords .> 0]
        neg_z = cpvals[zcoords .< 0]
        @test !isempty(pos_z)
        @test !isempty(neg_z)
        @test abs(mean(pos_z) - mean(neg_z)) < 0.02
    end
end
