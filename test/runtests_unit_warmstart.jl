using Test
import FLOWPanel as pnl
const TOML = pnl.TOML

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

struct WarmstartNoopSolver <: pnl.AbstractSolver end
pnl._solve!(::pnl.AbstractBody, ::WarmstartNoopSolver; kwargs...) = nothing

@testset "smooth conversion warm start crosses a conversion boundary" begin
    import FastMultipole

    function setup_smooth_case()
        base = make_plate_vortex_body()
        body = pnl.RigidWakeBody{Union{pnl.ConstantSource,pnl.VortexRing}}(
            copy(base.nodes), copy(base.cells), deepcopy(base.shedding);
            check_mesh=false, watertight=false)
        for i in eachindex(body.Das)
            body.Das[i] .= repeat([1.0, 0.0, 0.0], 1, size(body.Das[i],2))
        end
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.strength[:,2] .= 1.0
        conversion = pnl.SurfaceVorticityConversion(0.08; overlap=1.3,
            attribution=:split)
        wake = pnl.PanelParticleWake(body; nwakerows=2,
            max_particles=20000, conversion)
        frames = pnl.ReferenceFrame(body;
            v=FastMultipole.SVector(0.01,0.0,0.0), name="vehicle")
        return body,wake,frames
    end

    Uinf = t -> [1.0,0.0,0.0]
    maneuver = (frames,systems,wakes,t) -> nothing
    t_range = collect(range(0.0; step=0.05, length=8))
    opts = (body_solvers=(WarmstartNoopSolver(),),
        backend=pnl.DirectBackend(), grad_mu_options=(;basis=:tri), name="run")

    body_a,wake_a,frames_a = setup_smooth_case()
    path_a = mktempdir()
    pnl.simulate!((body_a,),(wake_a,),frames_a,maneuver,Uinf,t_range;
        opts...,path=path_a)

    body_b,wake_b,frames_b = setup_smooth_case()
    path_b = mktempdir()
    pnl.simulate!((body_b,),(wake_b,),frames_b,maneuver,Uinf,t_range[1:5];
        opts...,path=path_b)
    @test wake_b.conversion_count[] > 0

    body_c,wake_c,frames_c = setup_smooth_case()
    pnl.simulate_warmstart!((body_c,),(wake_c,),frames_c,maneuver,Uinf,t_range;
        opts...,path=path_b)
    @test body_a.nodes == body_c.nodes
    @test body_a.strength == body_c.strength
    @test wake_a.panel_wake.nodes == wake_c.panel_wake.nodes
    @test wake_a.panel_wake.strength == wake_c.panel_wake.strength
    @test wake_a.panel_wake.particle_handoff_active[] ==
          wake_c.panel_wake.particle_handoff_active[]
    @test wake_a.panel_wake.particle_handoff_weight[] ==
          wake_c.panel_wake.particle_handoff_weight[]
    @test wake_a.conversion_count[] == wake_c.conversion_count[]
    @test wake_a.pfield.np == wake_c.pfield.np
    @test view(wake_a.pfield.particles,:,1:wake_a.pfield.np) ==
          view(wake_c.pfield.particles,:,1:wake_c.pfield.np)

    # Missing and identity-mismatched smooth continuation state are typed hard
    # failures, before the saved step's end-of-step shedding can be replayed.
    metadata_file = joinpath(path_b,"run.metadata.toml")
    original = TOML.parsefile(metadata_file)
    missing = deepcopy(original)
    delete!(missing["step"][end],"wake_continuation")
    open(metadata_file,"w") do io
        TOML.print(io,missing)
    end
    body_d,wake_d,frames_d = setup_smooth_case()
    @test_throws pnl.WakeContinuationStateError pnl.simulate_warmstart!(
        (body_d,),(wake_d,),frames_d,maneuver,Uinf,t_range;
        opts...,path=path_b,restart_step=4)

    mismatch = deepcopy(original)
    mismatch["step"][end]["wake_continuation"][1]["conversion_fingerprint"] = "wrong"
    open(metadata_file,"w") do io
        TOML.print(io,mismatch)
    end
    body_e,wake_e,frames_e = setup_smooth_case()
    @test_throws pnl.WakeContinuationStateError pnl.simulate_warmstart!(
        (body_e,),(wake_e,),frames_e,maneuver,Uinf,t_range;
        opts...,path=path_b,restart_step=4)
end

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
    grad_mu_options = (; basis=:tri)

    body_A, wake_A, frames_A = setup_warmstart_case()
    path_A = mktempdir()
    pnl.simulate!((body_A,), (wake_A,), frames_A, maneuver, Uinf, t_range;
        body_solvers=(WarmstartNoopSolver(),),
        backend=pnl.DirectBackend(),
        grad_mu_options,
        name="run",
        path=path_A,
    )

    body_B, wake_B, frames_B = setup_warmstart_case()
    path_B = mktempdir()
    pnl.simulate!((body_B,), (wake_B,), frames_B, maneuver, Uinf, t_range[1:5];
        body_solvers=(WarmstartNoopSolver(),),
        backend=pnl.DirectBackend(),
        grad_mu_options,
        name="run",
        path=path_B,
    )

    body_C, wake_C, frames_C = setup_warmstart_case()
    pnl.simulate_warmstart!((body_C,), (wake_C,), frames_C, maneuver, Uinf, t_range;
        body_solvers=(WarmstartNoopSolver(),),
        backend=pnl.DirectBackend(),
        grad_mu_options,
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

    pnl._write_frame_state_toml(path_B, "run", frames_B, 4, t_range[5]; truncate=true)
    rm(joinpath(path_B, "run.metadata.toml"); force=true)

    body_D, wake_D, frames_D = setup_warmstart_case()
    pnl.simulate_warmstart!((body_D,), (wake_D,), frames_D, maneuver, Uinf, t_range;
        body_solvers=(WarmstartNoopSolver(),),
        backend=pnl.DirectBackend(),
        grad_mu_options,
        name="run",
        path=path_B,
        restart_step=4,
    )

    @test body_C.nodes == body_D.nodes
    @test body_C.strength == body_D.strength
    @test frames_C[1].x == frames_D[1].x
    @test wake_C.panel_wake.nwakes[] == wake_D.panel_wake.nwakes[]
end

@testset "convert-at-shed nwakerows=0 warm start (BRAINSTORM 024)" begin
    import FastMultipole

    function setup_n0_case()
        base = make_plate_vortex_body()
        body = pnl.RigidWakeBody{Union{pnl.ConstantSource,pnl.VortexRing}}(
            copy(base.nodes), copy(base.cells), deepcopy(base.shedding);
            check_mesh=false, watertight=false)
        for i in eachindex(body.Das)
            body.Das[i] .= repeat([1.0, 0.0, 0.0], 1, size(body.Das[i], 2))
        end
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.strength[:, 2] .= 1.0
        conversion = pnl.SurfaceVorticityConversion(0.08; overlap=1.3,
            attribution=:downstream)
        wake = pnl.PanelParticleWake(body; nwakerows=0,
            max_particles=20000, conversion)
        frames = pnl.ReferenceFrame(body;
            v=FastMultipole.SVector(0.01, 0.0, 0.0), name="vehicle")
        return body, wake, frames
    end

    Uinf = t -> [1.0, 0.0, 0.0]
    maneuver = (frames, systems, wakes, t) -> nothing
    t_range = collect(range(0.0; step=0.05, length=8))
    opts = (body_solvers=(WarmstartNoopSolver(),),
        backend=pnl.DirectBackend(), grad_mu_options=(; basis=:tri), name="run")

    body_a, wake_a, frames_a = setup_n0_case()
    path_a = mktempdir()
    pnl.simulate!((body_a,), (wake_a,), frames_a, maneuver, Uinf, t_range;
        opts..., path=path_a)
    @test wake_a.panel_wake.nwakes[] == 0
    @test wake_a.conversion_count[] == length(t_range) - 1

    body_b, wake_b, frames_b = setup_n0_case()
    path_b = mktempdir()
    pnl.simulate!((body_b,), (wake_b,), frames_b, maneuver, Uinf, t_range[1:5];
        opts..., path=path_b)
    @test wake_b.conversion_count[] > 0

    body_c, wake_c, frames_c = setup_n0_case()
    pnl.simulate_warmstart!((body_c,), (wake_c,), frames_c, maneuver, Uinf,
        t_range; opts..., path=path_b)
    @test body_a.nodes == body_c.nodes
    @test body_a.strength == body_c.strength
    @test wake_a.panel_wake.nwakes[] == wake_c.panel_wake.nwakes[] == 0
    @test wake_a.panel_wake.overflowed[] == wake_c.panel_wake.overflowed[] == true
    @test wake_a.panel_wake.nodes == wake_c.panel_wake.nodes
    @test wake_a.panel_wake.strength == wake_c.panel_wake.strength
    @test wake_a.panel_wake.particle_handoff_active[] ==
          wake_c.panel_wake.particle_handoff_active[]
    @test wake_a.panel_wake.particle_handoff_weight[] ==
          wake_c.panel_wake.particle_handoff_weight[]
    @test wake_a.conversion_count[] == wake_c.conversion_count[]
    @test wake_a.pfield.np == wake_c.pfield.np
    @test view(wake_a.pfield.particles, :, 1:wake_a.pfield.np) ==
          view(wake_c.pfield.particles, :, 1:wake_c.pfield.np)
end
