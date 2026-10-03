#=##############################################################################
p033 R4 thread-scaling campaign (Ryan 2026-10-02): BRAINSTORM-018 wake-physics
defaults ported onto the 021 warm-start vehicle.

Every knob below is an existing env var of examples/rotor_hover_pressure_
comparison.jl (which rotor_hover_solver_unsteady.jl includes with
RHPC_SETUP_ONLY=1) — this file only sets DEFAULTS via _setdefault!, so an
explicit value on the submit line always wins. Activated by WAKE_ENV_018=1;
inert otherwise, so the historical 021 arms replay unchanged.

Source of the values: 018 ntladder r3 common env + the r4-plan operator
cleanup (split/merge off, exact-rate corrected-Pedrizzetti relaxation at
NT=36, SFS off as the fixed-Cs bracket — no D3.4 Cs value exists yet).
SIGMA_FLOOR_R=0 per the 018 case env (the r4 plan's 0.25 guard floor is part
of the unexecuted r4 ladder). FLOWPANEL_FILAMENT_REG=linegauss is operator-
level and is exported by the launchers, not here.

Include AFTER _setdefault! is defined (fgs_r4_warmstart_ab.jl).
=###############################################################################

if get(ENV, "WAKE_ENV_018", "0") == "1"
    const WAKE_ENV_018_DEFAULTS = [
        # shedding / discretization (PARTICLE_SHEDDING=sigma_overlap is the
        # p018_* dispatcher-level setting SIGMA_CHORD_FRACTION hard-requires)
        "PARTICLE_SHEDDING"          => "sigma_overlap",
        "OVERLAP"                    => "2.75",
        "P_PER_STEP"                 => "12",
        "NWAKEROWS"                  => "1",
        "DAS_ETA_KINEMATIC"          => "1.0",
        "TRUNCATION_DEPTH_R"         => "4",
        # sigma law (frozen 019/018 ruling: sigma_j = 0.313 c_j, lambda=3 arc)
        "SIGMA_CHORD_FRACTION"       => "0.313",
        "SIGMA_FLOOR_R"              => "0",
        "DAS_SIGMA_LAMBDA"           => "3.0",
        "DAS_ARC_PLACED"             => "true",
        "DAS_ARC_HELIX_SOURCE"       => "steady",
        "DAS_ARC_TABLE"              => joinpath("data",
                                        "p018_cs_l3p4_rs1_te_downwash_te.csv"),
        # core spreading / regularization
        "CORE_SPREADING_ACTIVE"      => "true",
        "WAKE_CORE_BETA"             => "1e9",
        # domain / population
        "TRUNCATION_RADIUS_R"        => "3.0",
        "PARTICLE_OMIT_ROOT_R_OVER_R" => "0.12",
        "MAX_PARTICLES"              => "1500000",
        # operator cleanup (018 r4 plan: ladder-confounding operators OFF)
        "WAKE_SPLIT_VISCOUS"         => "false",
        "WAKE_SPLIT_STRETCH"         => "false",
        "MERGE_PARTICLES"            => "false",
        # relaxation: corrected-Pedrizzetti at the exact rate r(NT)=
        # 1-(1-0.3)^(36/NT); 0.3 IS the rate at NT=36 — recompute if NT changes
        "RELAX_SCHEME"               => "correctedpedrizzetti",
        "RELAX_RLXF"                 => "0.3",
        # SFS bracket arm (fixed-coefficient ConstantSFS awaits a D3.4 Cs)
        "SFS_OFF"                    => "true",
    ]
    for (k, v) in WAKE_ENV_018_DEFAULTS
        _setdefault!(k, v)
    end

    # The Das arc table is load-bearing (its absence killed ntladder r1 ~1 min
    # into every arm) and content-pinned: fail fast on a missing or wrong file.
    import SHA as _P033_SHA
    let das = ENV["DAS_ARC_TABLE"],
        expected = "640ba059cf57d6456ded0bf65721326160d90987fa449d1b6cc276d69fe755bf"
        isfile(das) || error("WAKE_ENV_018: Das arc table not found at " *
            "$(abspath(das)) — deploy step missing (data symlink?)")
        digest = bytes2hex(open(_P033_SHA.sha256, das))
        digest == expected || error("WAKE_ENV_018: Das arc table sha256 " *
            "mismatch at $(abspath(das)): got $digest, expected $expected " *
            "(md5 08375291ed2b542ea946d09730e0b629)")
    end

    println("WAKE PHYSICS (018-ported, WAKE_ENV_018=1; defaults, env wins):")
    for (k, _) in WAKE_ENV_018_DEFAULTS
        println("  $k = $(ENV[k])")
    end
end
