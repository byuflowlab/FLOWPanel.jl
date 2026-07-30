using FLOWPanel
const pnl = FLOWPanel

# reduce number of shedding locations to make debugging easier
body.shedding = reshape(body.shedding[:,1], :, 1)
body.nsheddings = 1

include_wake = true

# recalculate G matrix
normals = pnl._calc_normals(body)
CPs = pnl._calc_controlpoints(body, normals)
G1 = deepcopy(solver.G)
G1 .= 0.0
backend = pnl.FastMultipoleBackend(leaf_size=1000000)
# pnl._G_phi!(body, pnl.ConstantSource, G1, CPs, backend; kerneloffset=1.0e-8, include_wake)
pnl._G_phi!(body, pnl.ConstantDoublet, G1, CPs, backend; kerneloffset=1.0e-8, include_wake)

function test_rigid_wake(body, include_wake, G1, CPs; ichoose=body.ncells-1, istart=1, imax = body.ncells)

    body.strength[:,1] .= zeros(body.ncells)
    test_rhs = zeros(body.ncells)
    tfmm, tmat = 0.0, 0.0

    # for i in istart:imax
    #     if i % 2 == 0
    #         @show i, tfmm/(i-1), tmat/(i-1)
    #     end
    #     test_rhs .= 0.0
    #     body.strength[:,2] .= 0.0
    #     body.strength[i,2] = 1.0
    #     # body.strength[:,2] .= rand(body.ncells)
    #     # body.strength[1,2] = 1.0
    #     # backend = pnl.DirectBackend()
    #     tfmm += @elapsed pnl._phi!(body, CPs, test_rhs, backend; kerneloffset=1.0e-8, include_wake)
    #     tmat += @elapsed debug_rhs = G1 * body.strength[:,2]
    #     resid = maximum(abs.(test_rhs .- debug_rhs))
    #     if resid > 1e-6
    #         @show i, resid
    #     end

    #     # _rigid_wake_phi!(body, targets, out; optargs...)
    # end

    # body.strength[:,2] .= 0.0 #rand(body.ncells)
    body.strength[:,2] = rand(body.ncells)
    pnl._phi!(body, CPs, test_rhs, backend; kerneloffset=1.0e-8, include_wake)
    debug_rhs = G1 * body.strength[:,2]
    resid = maximum(abs.(test_rhs .- debug_rhs))
    @show resid

    return resid, debug_rhs, test_rhs
end

resid, debug_rhs, test_rhs = test_rigid_wake(body, include_wake, G1, CPs; ichoose=2170, istart=body.ncells-2, imax = body.ncells)