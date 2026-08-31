using Test

const FLOWPANEL_ROOT = normpath(joinpath(@__DIR__, ".."))
const DRIVER = joinpath(FLOWPANEL_ROOT, "examples", "rotor_hover_ground_effect.jl")

function setup_probe(formulation, nrotors, ground; extra=Dict{String,String}())
    env = copy(ENV)
    merge!(env, Dict(
        "RHPC_SETUP_ONLY" => "true",
        "RHPC_FORMULATION" => formulation,
        "NROTORS" => string(nrotors),
        "GROUND_ENABLE" => string(ground),
        "BERNOULLI_ONLY" => "true",
        "RUN_MONITORS" => "true",
        "SAVE_VTK" => "false",
        "RHPC_MESH" => "40_40",
        "RHPC_BACKEND" => "direct",
        "PARTICLE_BODY_OVERLAP_ACTION" => "off",
        "FLOWPANEL_FILAMENT_REG" => "gaussian",
        "NT" => "36",
        "NREVS" => "10",
        "SPINUP_REVS" => "1.5",
        "FREESTREAM_RAMP_REVS" => "1.0",
        "FREESTREAM_HOLD_REVS" => "1.5",
        "FREESTREAM_WITHDRAW_REVS" => "4.0",
        "SETTLE_REVS" => "3.5",
        "GROUND_H_R" => "1.5",
        "GROUND_PARTICLE_POLICY" => "none",
        "GROUND_DAMP_BAND_R" => "0.1",
        "TRUNCATION_DEPTH_R" => "4.5",
        "PARTICLE_SHEDDING" => "sigma_overlap",
        "OVERLAP" => "2.75",
        "P_PER_STEP" => "12",
        "NWAKEROWS" => "1",
        "JULIA_NUM_THREADS" => "2",
    ))
    merge!(env, extra)
    cmd = setenv(`$(Base.julia_cmd()) --project=$FLOWPANEL_ROOT $DRIVER`, env)
    io = IOBuffer()
    process = run(pipeline(ignorestatus(cmd), stdout=io, stderr=io))
    return process.exitcode, String(take!(io))
end

@testset "rotor-hover production setup matrix (052b C)" begin
    # VTS is the 052b production formulation: cover every accepted rotor-count
    # x OGE/IGE shape. Each arm constructs independent simulation state.
    for nrotors in (1, 2, 4), ground in (false, true)
        code, output = setup_probe("velocity", nrotors, ground)
        @test code == 0
        @test occursin("RHPC_SETUP_OK formulation=velocity", output)
        @test occursin("nrotors=$nrotors", output)
        @test occursin("ground=$ground", output)
        @test occursin("require_outer_convergence=true", output)
        @test occursin("circulation_monitors=$nrotors", output)
        @test occursin("Total run length: 10.0 revs (414 steps)", output)
    end

    # Hybrid remains available as an independent 052e comparison, but is not
    # an acceptance dependency.
    for ground in (false, true)
        code, output = setup_probe("hybrid", 4, ground)
        @test code == 0
        @test occursin("RHPC_SETUP_OK formulation=hybrid", output)
        @test occursin("require_outer_convergence=true", output)
        @test occursin("circulation_monitors=4", output)
    end

    # Exact production handedness and explicitly validated owner classes:
    # +1 rotors 1/3 share owner 1; -1 rotors 2/4 share owner 2.
    for ground in (false, true)
        code, output = setup_probe("velocity", 4, ground; extra=Dict(
            "ROTOR_DIRECTIONS" => "1,-1,1,-1",
            "ROTOR_OPERATOR_OWNERS" => "1,2,1,2",
            "ROTOR_OPERATOR_CLASSES_VALIDATED" => "true"))
        @test code == 0
        @test occursin("directions=[1, -1, 1, -1]", output)
        @test occursin("operator_owners=[1, 2, 1, 2]", output)
        @test occursin("classes_validated=true", output)
    end

    # The older Green reconstruction remains a single-rotor diagnostic and is
    # rejected, rather than silently substituted, for multirotor production.
    code, output = setup_probe("green", 2, false)
    @test code != 0
    @test occursin("requires RHPC_FORMULATION=velocity or hybrid", output)
end
