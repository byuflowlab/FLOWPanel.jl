using Test
import FLOWPanel as pnl
import GeometricTools as gt
using Statistics: mean

@testset verbose=true "Analytical Validation" begin
    @testset "Sphere in uniform flow - Neumann" begin
        body = make_sphere_source_body(ntheta=12, nphi=24)
        solve_source_body!(body)

        cps = body.controlpoints
        pvals = body.P
        xcoords = view(cps, 1, :)
        ycoords = view(cps, 2, :)
        zcoords = view(cps, 3, :)

        # With rho=1, Uinf=1: q = 0.5, so P_stag = q*1.0 = 0.5, P_equator = q*(-1.25) = -0.625
        i_stag = argmax(xcoords)
        @test isapprox(pvals[i_stag], 0.5; atol=0.025)

        eq_mask = abs.(xcoords) .< 0.05
        @test any(eq_mask)
        p_equator = mean(pvals[eq_mask])
        @test isapprox(p_equator, -0.625; atol=0.1)

        q = 0.5 * pi
        ftot = pnl.calcfield_Ftot!(zeros(3), body, body.F)
        @test abs(ftot[1]) / q < 0.05
        @test abs(ftot[3]) / q < 0.05

        pos_y = pvals[ycoords .> 0]
        neg_y = pvals[ycoords .< 0]
        @test !isempty(pos_y)
        @test !isempty(neg_y)
        @test abs(mean(pos_y) - mean(neg_y)) < 0.01

        pos_z = pvals[zcoords .> 0]
        neg_z = pvals[zcoords .< 0]
        @test !isempty(pos_z)
        @test !isempty(neg_z)
        @test abs(mean(pos_z) - mean(neg_z)) < 0.01
    end
end
