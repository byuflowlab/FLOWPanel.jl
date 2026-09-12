using Test
import FLOWPanel as pnl
const TOML = pnl.TOML

if !isdefined(@__MODULE__, :make_plate_vortex_body)
    include("test_helpers.jl")
end

struct WarmstartNoopSolver <: pnl.AbstractSolver end
pnl._solve!(::pnl.AbstractBody, ::WarmstartNoopSolver; kwargs...) = nothing

mutable struct WarmstartIGEAudit <: pnl.AbstractMonitor
    rotor_forces::Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}
    ground_strength::Dict{Int,Vector{Float64}}
    ground_tangency::Dict{Int,Tuple{Float64,Float64}}
    gs_status::Dict{Int,NamedTuple}
    force1::Any
    force2::Any
end

function (m::WarmstartIGEAudit)(systems, wakes, frames, uinf,
        i_step::Int, dt::Real)
    m.rotor_forces[i_step] =
        (copy(m.force1.force[:, i_step + 1]), copy(m.force2.force[:, i_step + 1]))
    ground = systems[3]
    m.ground_strength[i_step] = copy(ground.strength[:, 1])
    un = vec(sum(ground.velocity .* ground.normals; dims=1))
    m.ground_tangency[i_step] =
        (sqrt(sum(abs2, un) / length(un)), maximum(abs, un))
    m.gs_status[i_step] = pnl.block_gs_status(IGE_TEST_SOLVERS[])
    return nothing
end

const IGE_TEST_SOLVERS = Ref{Any}(())

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

@testset "multi-rotor IGE warm start: shared and legacy particle layouts" begin
    import FastMultipole

    function setup_ige_case(shared::Bool, ntime::Int)
        r1 = make_plate_vortex_body()
        r2 = make_plate_vortex_body()
        r2.nodes[1, :] .+= 2.0
        pnl.calc_normals!(r2)
        pnl.calc_controlpoints!(r2)
        ground = pnl.FlatGround([1.0, 0.5, -1.0], [0.0, 0.0, 1.0], 3.0;
            panel_length=1.5)

        w1 = pnl.PanelParticleWake(r1; nwakerows=100, max_particles=4000,
            freestream_convection=true, SFS=pnl.FLOWVPM.noSFS,
            relaxation=pnl.FLOWVPM.relaxation_none)
        w2 = pnl.PanelParticleWake(r2; nwakerows=100, max_particles=4000,
            freestream_convection=true,
            SFS=pnl.FLOWVPM.noSFS,
            relaxation=pnl.FLOWVPM.relaxation_none,
            pfield=shared ? w1.pfield : nothing)
        s1 = pnl.Backslash(r1)
        s2 = pnl.Backslash(r2; shared_operator=s1)
        sg = pnl.FlatGroundSolver(ground)
        solvers = (s1, s2, sg)
        IGE_TEST_SOLVERS[] = solvers

        pressure = pnl.PressureBernoulli(1.0; unsteady=false,
            correct_kuttacondition=false)
        f1 = pnl.ForceMonitor(ntime, 1; normalization=pnl.NoNormalization(),
            correct_kuttacondition=false)
        f2 = pnl.ForceMonitor(ntime, 2; normalization=pnl.NoNormalization(),
            correct_kuttacondition=false)
        audit = WarmstartIGEAudit(
            Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}(),
            Dict{Int,Vector{Float64}}(),
            Dict{Int,Tuple{Float64,Float64}}(),
            Dict{Int,NamedTuple}(), f1, f2)
        frames = pnl.ReferenceFrame(r1; name="vehicle")
        return (r1, r2, ground), (w1, w2, nothing), frames, solvers,
            (pressure, f1, f2, audit), audit
    end

    Uinf = t -> [0.2, 0.0, 0.1]
    t_range = collect(range(0.0; step=0.04, length=8))
    restart_step = 3
    maneuver = (frames, systems, wakes, t) -> nothing
    common = (backend=pnl.DirectBackend(),
        backend_wake=pnl.FastMultipoleBackend(expansion_order=4,
            multipole_acceptance=0.4, leaf_size=20),
        grad_mu_options=(; basis=:tri),
        name="ige", particle_hessian_self=false, particle_relax=false,
        panel_wake_on_particles=false)

    for shared in (true, false)
        # Isolated nonempty VTP round trip for both checkpoint layouts. Keep
        # these serialization probes out of the reduced physical continuation
        # so the restart comparison exercises only panel/ground dynamics.
        layout_systems, layout_wakes, _, _, _, _ =
            setup_ige_case(shared, length(t_range))
        pnl.FLOWVPM.add_particle(layout_wakes[1].pfield, (0.25, 0.25, 0.4),
            (0.0, 0.0, 0.01), 0.2)
        pnl.FLOWVPM.add_particle(layout_wakes[2].pfield, (2.25, 0.25, 0.4),
            (0.0, 0.0, -0.01), 0.2)
        for (body, wake) in zip(layout_systems[1:2], layout_wakes[1:2])
            pnl.update_TE!(wake, body)
            pnl.shed_wake!(wake, body)
        end
        layout_path = mktempdir()
        seen_layout = ()
        for (iw, wake) in enumerate(layout_wakes[1:2])
            repeat = any(p -> p === wake.pfield, seen_layout)
            pnl.write_vtk(joinpath(layout_path, "layout_wake$(iw)"), wake,
                restart_step, t_range[restart_step + 1];
                include_pfield=!repeat)
            repeat || (seen_layout = (seen_layout..., wake.pfield))
        end
        p1 = joinpath(layout_path, "layout_wake1_particles",
            "layout_wake1_particles.$(restart_step).vtp")
        p2 = joinpath(layout_path, "layout_wake2_particles",
            "layout_wake2_particles.$(restart_step).vtp")
        @test isfile(p1)
        @test isfile(p2) == !shared
        vtk1 = pnl.ReadVTK.VTKFile(p1)
        @test all(isfinite, pnl.ReadVTK.get_points(vtk1))
        pdata1 = pnl.ReadVTK.get_point_data(vtk1)
        for field in ("velocity", "gamma", "C", "SFS", "velocity_gradient")
            @test all(isfinite, pnl.ReadVTK.get_data(pdata1[field]))
        end
        _, probe_wakes, _, _, _, _ = setup_ige_case(shared, length(t_range))
        seen_probe = ()
        for (iw, wake) in enumerate(probe_wakes[1:2])
            repeat = any(p -> p === wake.pfield, seen_probe)
            pnl._load_panel_particle_wake_vtk!(wake, layout_path,
                "layout_wake$(iw)", restart_step; include_pfield=!repeat)
            repeat || (seen_probe = (seen_probe..., wake.pfield))
        end
        for (wl, wp) in zip(layout_wakes[1:2], probe_wakes[1:2])
            @test wl.pfield.np == wp.pfield.np
            @test view(wl.pfield.particles, 1:9, 1:wl.pfield.np) ≈
                  view(wp.pfield.particles, 1:9, 1:wp.pfield.np)
        end

        systems_a, wakes_a, frames_a, solvers_a, monitors_a, audit_a =
            setup_ige_case(shared, length(t_range))
        path_a = mktempdir()
        pnl.simulate!(systems_a, wakes_a, frames_a, maneuver, Uinf, t_range;
            body_solvers=solvers_a, monitors=monitors_a, path=path_a, common...)

        systems_b, wakes_b, frames_b, solvers_b, monitors_b, _ =
            setup_ige_case(shared, length(t_range))
        path_b = mktempdir()
        pnl.simulate!(systems_b, wakes_b, frames_b, maneuver, Uinf,
            t_range[1:restart_step + 1]; body_solvers=solvers_b,
            monitors=monitors_b, path=path_b, common...)

        systems_c, wakes_c, frames_c, solvers_c, monitors_c, audit_c =
            setup_ige_case(shared, length(t_range))
        pnl.simulate_warmstart!(systems_c, wakes_c, frames_c, maneuver, Uinf,
            t_range; restart_step, body_solvers=solvers_c,
            monitors=monitors_c, path=path_b, common...)

        continued = (restart_step + 1):(length(t_range) - 1)
        @test length(continued) >= 3
        for i in continued
            @test audit_c.rotor_forces[i][1] ≈ audit_a.rotor_forces[i][1] atol=1e-11 rtol=1e-10
            @test audit_c.rotor_forces[i][2] ≈ audit_a.rotor_forces[i][2] atol=1e-11 rtol=1e-10
            @test audit_c.ground_strength[i] ≈ audit_a.ground_strength[i] atol=1e-11 rtol=1e-10
            @test audit_c.ground_tangency[i][1] ≈ audit_a.ground_tangency[i][1] atol=1e-11 rtol=1e-10
            @test audit_c.ground_tangency[i][2] ≈ audit_a.ground_tangency[i][2] atol=1e-11 rtol=1e-10
            @test audit_c.gs_status[i] == audit_a.gs_status[i]
            @test audit_c.gs_status[i].converged
            @test audit_c.gs_status[i].iterations <= audit_c.gs_status[i].cap
        end

        for (wa, wc) in zip(wakes_a[1:2], wakes_c[1:2])
            @test wa.pfield.np == wc.pfield.np
            np = wa.pfield.np
            @test view(wa.pfield.particles, :, 1:np) ≈
                  view(wc.pfield.particles, :, 1:np) atol=1e-11 rtol=1e-10
            nw = wa.panel_wake.nwakes[]
            @test nw == wc.panel_wake.nwakes[]
            for isurf in eachindex(wa.panel_wake.nodes)
                @test view(wa.panel_wake.nodes[isurf], :, 1:nw + 1, :) ≈
                      view(wc.panel_wake.nodes[isurf], :, 1:nw + 1, :) atol=1e-11 rtol=1e-10
                @test view(wa.panel_wake.strength[isurf], :, 1:nw, :) ≈
                      view(wc.panel_wake.strength[isurf], :, 1:nw, :) atol=1e-11 rtol=1e-10
            end
        end
        if shared
            @test wakes_c[1].pfield === wakes_c[2].pfield
            @test wakes_c[1].pfield.np == wakes_a[1].pfield.np
        else
            @test wakes_c[1].pfield !== wakes_c[2].pfield
        end
    end
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

@testset "warm start ResolutionSplitState persistence (026 Phase 2)" begin
    vpm = pnl.FLOWVPM

    function setup_rsplit_wake()
        body = make_plate_vortex_body()
        wake = pnl.PanelParticleWake(body; nwakerows=2, max_particles=100,
            SFS=vpm.noSFS, relaxation=vpm.relaxation_none)
        pnl.update_TE!(wake, body)
        pnl.shed_wake!(wake, body)
        return wake
    end

    function seed_particles!(wake; np=4)
        for i in 1:np
            vpm.add_particle(wake.pfield, (0.3i, 0.1, 0.2),
                (0.01, 0.0, 0.02i), 0.1 + 0.01i)
        end
        return wake
    end

    @testset "A: full-block exact restore" begin
        wake = seed_particles!(setup_rsplit_wake())
        pf = wake.pfield
        np = pf.np
        rs = vpm.enable_resolution_split!(pf)
        for i in 1:np
            rs.sigma_0[i] = 0.05 + 0.01i
            rs.axis[1, i] = 0.1i; rs.axis[2, i] = -0.2i; rs.axis[3, i] = 0.3i
            rs.weight[i] = 0.4i
            rs.dvisc[i] = 1e-4 * i
            rs.drvpm[i] = 2e-4 * i
        end

        path = mktempdir()
        pnl.write_vtk(joinpath(path, "w1"), wake, 7, 0.35)

        wake2 = setup_rsplit_wake()
        @test wake2.pfield.resolution_split === nothing
        pnl._load_panel_particle_wake_vtk!(wake2, path, "w1", 7)
        pf2 = wake2.pfield
        rs2 = pf2.resolution_split
        @test pf2.np == np
        @test rs2 !== nothing
        @test rs2.sigma_0[1:np] == rs.sigma_0[1:np]
        @test rs2.axis[:, 1:np] == rs.axis[:, 1:np]
        @test rs2.weight[1:np] == rs.weight[1:np]
        @test rs2.dvisc[1:np] == rs.dvisc[1:np]
        @test rs2.drvpm[1:np] == rs.drvpm[1:np]
        # slots beyond np stay zero
        @test all(iszero, rs2.sigma_0[np+1:end])
        @test all(iszero, rs2.weight[np+1:end])
        @test all(iszero, rs2.axis[:, np+1:end])

        # fractional trigger fires identically after restore: give particle 1
        # a viscous credit past the f_visc=0.5 threshold on its restored σ₀
        opts = vpm.ResolutionSplitOpts(f_visc=0.5,
            enable_viscous_split=true, enable_stretch_split=true)
        rs2.dvisc[1] = 2.0 * rs2.sigma_0[1]^2   # > ((1.5)² − 1)·σ₀²
        counters = vpm.split_particles!(pf2, opts)
        @test counters.n_split_viscous == 1
        @test pf2.np == np + 3
    end

    @testset "B: absent block leaves state nothing (feature off)" begin
        wake = seed_particles!(setup_rsplit_wake())
        @test wake.pfield.resolution_split === nothing
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "w1"), wake, 3, 0.1)

        wake2 = setup_rsplit_wake()
        pnl._load_panel_particle_wake_vtk!(wake2, path, "w1", 3)
        @test wake2.pfield.np == wake.pfield.np
        @test wake2.pfield.resolution_split === nothing
    end

    @testset "C: partial block is a typed hard failure" begin
        wake = seed_particles!(setup_rsplit_wake())
        pf = wake.pfield
        np = pf.np
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "w1"), wake, 4, 0.2)
        vtp_dir = joinpath(path, "w1_particles")
        cells = [pnl.WriteVTK.MeshCell(pnl.WriteVTK.PolyData.Verts(), 1:np)]
        Xp = pf.particles[vpm.X_INDEX, 1:np]
        vtp = pnl.WriteVTK.vtk_grid(joinpath(vtp_dir, "w1_particles.4.vtp"),
            Xp, cells; compress=false)
        for (fname, idxs) in (("gamma", vpm.GAMMA_INDEX),
                ("sigma", vpm.SIGMA_INDEX), ("vol", vpm.VOL_INDEX),
                ("circulation", vpm.CIRCULATION_INDEX),
                ("velocity", vpm.U_INDEX), ("vorticity", vpm.VORTICITY_INDEX),
                ("C", vpm.C_INDEX), ("SFS", vpm.SFS_INDEX))
            vtp[fname, pnl.WriteVTK.VTKPointData()] = pf.particles[idxs, 1:np]
        end
        vtp["velocity_gradient", pnl.WriteVTK.VTKPointData()] =
            reshape(pf.particles[vpm.J_INDEX, 1:np], 3, 3, np)
        vtp["rsplit_sigma_0", pnl.WriteVTK.VTKPointData()] = fill(0.1, np)
        vtp["rsplit_weight", pnl.WriteVTK.VTKPointData()] = fill(0.2, np)
        pnl.WriteVTK.vtk_save(vtp)

        wake2 = setup_rsplit_wake()
        @test_throws ArgumentError pnl._load_panel_particle_wake_vtk!(
            wake2, path, "w1", 4)
    end

    @testset "E: legacy 6-field save (with rsplit_exposure) migrates" begin
        # Pre-fractional-gating checkpoints carry rsplit_exposure and
        # applied-Δσ² accumulators: sigma_0/axis/weight restore, exposure is
        # dropped, dvisc/drvpm are zeroed (semantics changed to attempted
        # Δσ², Ryan 2026-09-08) with a warning.
        wake = seed_particles!(setup_rsplit_wake())
        pf = wake.pfield
        np = pf.np
        rs = vpm.enable_resolution_split!(pf)
        for i in 1:np
            rs.sigma_0[i] = 0.05 + 0.01i
            rs.axis[1, i] = 0.1i; rs.axis[2, i] = -0.2i; rs.axis[3, i] = 0.3i
            rs.weight[i] = 0.4i
            rs.dvisc[i] = 1e-4 * i
            rs.drvpm[i] = 2e-4 * i
        end
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "w1"), wake, 5, 0.25)
        # graft the legacy exposure array onto the fresh save
        vtp_dir = joinpath(path, "w1_particles")
        cells = [pnl.WriteVTK.MeshCell(pnl.WriteVTK.PolyData.Verts(), 1:np)]
        Xp = pf.particles[vpm.X_INDEX, 1:np]
        vtp = pnl.WriteVTK.vtk_grid(joinpath(vtp_dir, "w1_particles.5.vtp"),
            Xp, cells; compress=false)
        for (fname, idxs) in (("gamma", vpm.GAMMA_INDEX),
                ("sigma", vpm.SIGMA_INDEX), ("vol", vpm.VOL_INDEX),
                ("circulation", vpm.CIRCULATION_INDEX),
                ("velocity", vpm.U_INDEX), ("vorticity", vpm.VORTICITY_INDEX),
                ("C", vpm.C_INDEX), ("SFS", vpm.SFS_INDEX))
            vtp[fname, pnl.WriteVTK.VTKPointData()] = pf.particles[idxs, 1:np]
        end
        vtp["velocity_gradient", pnl.WriteVTK.VTKPointData()] =
            reshape(pf.particles[vpm.J_INDEX, 1:np], 3, 3, np)
        vtp["rsplit_sigma_0", pnl.WriteVTK.VTKPointData()] = rs.sigma_0[1:np]
        vtp["rsplit_axis", pnl.WriteVTK.VTKPointData()] = rs.axis[:, 1:np]
        vtp["rsplit_weight", pnl.WriteVTK.VTKPointData()] = rs.weight[1:np]
        vtp["rsplit_exposure", pnl.WriteVTK.VTKPointData()] = fill(0.3, np)
        vtp["rsplit_dvisc", pnl.WriteVTK.VTKPointData()] = rs.dvisc[1:np]
        vtp["rsplit_drvpm", pnl.WriteVTK.VTKPointData()] = rs.drvpm[1:np]
        pnl.WriteVTK.vtk_save(vtp)

        wake2 = setup_rsplit_wake()
        @test_logs (:warn, r"legacy") match_mode=:any (
            pnl._load_panel_particle_wake_vtk!(wake2, path, "w1", 5))
        rs2 = wake2.pfield.resolution_split
        @test rs2 !== nothing
        @test rs2.sigma_0[1:np] == rs.sigma_0[1:np]
        @test rs2.axis[:, 1:np] == rs.axis[:, 1:np]
        @test rs2.weight[1:np] == rs.weight[1:np]
        @test all(iszero, rs2.dvisc[1:np])      # accumulators re-armed
        @test all(iszero, rs2.drvpm[1:np])
    end

    @testset "D: stale legacy split_* fields ignored silently (D-B)" begin
        # Checkpoints written before the legacy splitting path was removed
        # (2026-09-08) carry six split_* point-data arrays. Fabricate one and
        # verify it loads with no warning/error — the stale fields are inert
        # extra point data — while the rsplit_* block still restores exactly.
        wake = seed_particles!(setup_rsplit_wake())
        pf = wake.pfield
        np = pf.np
        path = mktempdir()
        pnl.write_vtk(joinpath(path, "w1"), wake, 5, 0.25)
        vtp_dir = joinpath(path, "w1_particles")
        cells = [pnl.WriteVTK.MeshCell(pnl.WriteVTK.PolyData.Verts(), 1:np)]
        Xp = pf.particles[vpm.X_INDEX, 1:np]
        vtp = pnl.WriteVTK.vtk_grid(joinpath(vtp_dir, "w1_particles.5.vtp"),
            Xp, cells; compress=false)
        for (fname, idxs) in (("gamma", vpm.GAMMA_INDEX),
                ("sigma", vpm.SIGMA_INDEX), ("vol", vpm.VOL_INDEX),
                ("circulation", vpm.CIRCULATION_INDEX),
                ("velocity", vpm.U_INDEX), ("vorticity", vpm.VORTICITY_INDEX),
                ("C", vpm.C_INDEX), ("SFS", vpm.SFS_INDEX))
            vtp[fname, pnl.WriteVTK.VTKPointData()] = pf.particles[idxs, 1:np]
        end
        vtp["velocity_gradient", pnl.WriteVTK.VTKPointData()] =
            reshape(pf.particles[vpm.J_INDEX, 1:np], 3, 3, np)
        # legacy Session-2-era block (exact field names/types the old writer
        # used; the rsplit_* block below is ALSO legacy-format — it carries
        # rsplit_exposure — so the loader takes the migration path)
        vtp["split_sigma_0", pnl.WriteVTK.VTKPointData()] = fill(0.3, np)
        vtp["split_H_chi", pnl.WriteVTK.VTKPointData()] = fill(0.1, np)
        vtp["split_hold", pnl.WriteVTK.VTKPointData()] = fill(Int32(2), np)
        vtp["split_cooldown", pnl.WriteVTK.VTKPointData()] = fill(Int32(1), np)
        vtp["split_dsigma2_visc", pnl.WriteVTK.VTKPointData()] = fill(1e-4, np)
        vtp["split_dsigma2_rvpm", pnl.WriteVTK.VTKPointData()] = fill(-2e-4, np)
        # plus a full rsplit_* block that must still restore
        vtp["rsplit_sigma_0", pnl.WriteVTK.VTKPointData()] = collect(0.05 .+ 0.01*(1:np))
        vtp["rsplit_axis", pnl.WriteVTK.VTKPointData()] = repeat([0.1, -0.2, 0.3], 1, np)
        vtp["rsplit_weight", pnl.WriteVTK.VTKPointData()] = fill(0.4, np)
        vtp["rsplit_exposure", pnl.WriteVTK.VTKPointData()] = fill(-0.05, np)
        vtp["rsplit_dvisc", pnl.WriteVTK.VTKPointData()] = fill(1e-4, np)
        vtp["rsplit_drvpm", pnl.WriteVTK.VTKPointData()] = fill(2e-4, np)
        pnl.WriteVTK.vtk_save(vtp)

        wake2 = setup_rsplit_wake()
        # stale split_* fields stay silently inert; the rsplit_exposure field
        # makes this a legacy save, so ONLY the migration warning fires
        @test_logs (:warn, r"legacy") match_mode=:any (
            pnl._load_panel_particle_wake_vtk!(wake2, path, "w1", 5))
        pf2 = wake2.pfield
        @test pf2.np == np
        rs2 = pf2.resolution_split
        @test rs2 !== nothing
        @test rs2.sigma_0[1:np] == collect(0.05 .+ 0.01*(1:np))
        @test all(rs2.weight[1:np] .== 0.4)
        @test all(iszero, rs2.dvisc[1:np])   # legacy accumulators re-armed
        @test all(iszero, rs2.drvpm[1:np])
    end
end
