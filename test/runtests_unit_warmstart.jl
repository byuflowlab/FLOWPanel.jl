using Test
import FLOWPanel as pnl

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

struct WarmstartNoopSolver <: pnl.AbstractSolver end
pnl.solve!(::pnl.AbstractBody, ::WarmstartNoopSolver; kwargs...) = nothing

@testset "simulate_warmstart! consistency (PanelParticleWake)" begin
    import FastMultipole

    function make_warmstart_body()
        base = make_plate_vortex_body()
        body = pnl.RigidWakeBody{Union{pnl.ConstantSource,pnl.VortexRing}}(
            copy(base.nodes),
            copy(base.cells),
            deepcopy(base.shedding);
            check_mesh=false,
            watertight=false,
        )
        for i in eachindex(body.Das)
            body.Das[i] .= repeat([1.0, 0.0, 0.0], 1, size(body.Das[i], 2))
        end
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.strength[:, 2] .= 1.0
        return body
    end

    function setup_warmstart_case()
        body = make_warmstart_body()
        wake = pnl.PanelParticleWake(body; nwakerows=2, max_particles=500)
        frames = pnl.ReferenceFrame(body;
            v=FastMultipole.SVector(0.01, 0.0, 0.0),
            name="vehicle",
        )
        return body, wake, frames
    end

    Uinf = t -> [1.0, 0.0, 0.0]
    maneuver = (frames, systems, wakes, t) -> nothing
    t_range = collect(range(0.0; step=0.05, length=10))

    body_A, wake_A, frames_A = setup_warmstart_case()
    path_A = mktempdir()
    pnl.simulate!((body_A,), (wake_A,), frames_A, maneuver, Uinf, t_range;
        body_solvers=(WarmstartNoopSolver(),),
        backend=pnl.DirectBackend(),
        name="run",
        path=path_A,
    )

    body_B, wake_B, frames_B = setup_warmstart_case()
    path_B = mktempdir()
    pnl.simulate!((body_B,), (wake_B,), frames_B, maneuver, Uinf, t_range[1:5];
        body_solvers=(WarmstartNoopSolver(),),
        backend=pnl.DirectBackend(),
        name="run",
        path=path_B,
    )

    body_C, wake_C, frames_C = setup_warmstart_case()
    pnl.simulate_warmstart!((body_C,), (wake_C,), frames_C, maneuver, Uinf, t_range;
        body_solvers=(WarmstartNoopSolver(),),
        backend=pnl.DirectBackend(),
        name="run",
        path=path_B,
    )

    @test body_A.nodes == body_C.nodes
    @test body_A.strength == body_C.strength
    @test body_A.potential == body_C.potential
    @test frames_A[1].x == frames_C[1].x
    @test frames_A[1].R == frames_C[1].R

    @test wake_A.panel_wake.nwakes[] == wake_C.panel_wake.nwakes[]
    for i in eachindex(wake_A.panel_wake.nodes)
        @test wake_A.panel_wake.nodes[i] == wake_C.panel_wake.nodes[i]
        @test wake_A.panel_wake.strength[i] == wake_C.panel_wake.strength[i]
    end

    @test wake_A.pfield.np == wake_C.pfield.np
    np = wake_A.pfield.np
    @test view(wake_A.pfield.particles, :, 1:np) == view(wake_C.pfield.particles, :, 1:np)
end
