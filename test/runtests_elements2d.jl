#=##############################################################################
# DESCRIPTION
    Unit tests of elements
=###############################################################################

using Test
import Printf: @printf
import FLOWPanel as pnl
import FLOWPanel: gt

try
    verbose
catch
    global verbose = true
end
v_lvl = 0


@testset verbose=verbose "2D Element Tests" begin

    offset = 1e-12

    # Nodes of planar panel
    nodes = [
                -0.5 0.0;
                0.5 0.0;
            ]
    nodes = collect(nodes')
    nnodes = size(nodes, 2)

    panel = [1, 2]


    panel_length = pnl.norm(nodes[:, panel[2]] - nodes[:, panel[1]])

    # Translate and rotate panel
    rotation = 27
    O = [1.0, 1.0, 0.0]

    Oaxis = gt.rotation_matrix2(0, 0, rotation)
    O2d, Oaxis2d = O[1:2], Oaxis[1:2, 1:2]
    nodes = Oaxis2d'*nodes + repeat(O2d, 1, nnodes)

    # Panel strengths
    gamma1 = -0.53
    gamma2 = 0.76

    normal = [gt._calc_n1_2D(nodes, panel), gt._calc_n2_2D(nodes, panel)]
    tangent = [gt._calc_t1_2D(nodes, panel), gt._calc_t2_2D(nodes, panel)]

    # Targets where to prove the velocity
    targets = hcat(
                    nodes[:, 1] + 1e-6*tangent,
                    nodes[:, 2] - 1e-5*tangent,
                    (nodes[:, 2] + nodes[:, 1])/2,
                    nodes[:, 1] + 0.123*(nodes[:, 2] - nodes[:, 1]),
                  )
    ntargets = size(targets, 2)

    # Array where to store velocities
    Us = zeros(3, ntargets)

    # Calculate velocities
    pnl.U_2D_linear_vortex(nodes, panel,
                                gamma1, gamma2,
                                targets, Us;
                                offset=offset
                              )

    @test abs(pnl.dot(Us[:, 3], tangent) - (gamma1 + gamma2)/4) <= offset
    @test abs(pnl.dot(Us[:, 3], normal) - (gamma2 - gamma1)/(2*pi)) <= offset
    @test abs(pnl.dot(Us[:, 1], tangent) - gamma1/2) <= 1e7*offset
    @test abs(pnl.dot(Us[:, 2], tangent) - gamma2/2) <= 1e7*offset

end
