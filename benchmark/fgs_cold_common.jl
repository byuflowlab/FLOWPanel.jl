# Shared implementation for the frozen cold-solve timing and profiling drivers.
# Definitions only: cold_initialize! owns fixture loading and filesystem effects.
using TOML, SHA, Statistics, Serialization, Profile

const COLD_TARGET = 1e-6
const COLD_SEEDS = Dict("R1" => (6, 0.3, 50, 10, 17),
                        "R2" => (8, 0.4, 100, 10, 17),
                        "R3" => (6, 0.3, 100, 5, 16),
                        # Provisional R2-derived starting point, NOT tuned R4 knobs.
                        "R4" => (8, 0.4, 100, 3, 17))

function cold_seed(rung, kind)
    p, mac, leaf, inner, kp = COLD_SEEDS[rung]
    common = Dict{String,Any}("kind" => kind, "rung" => rung,
        "P" => kind == "fgs" ? p : kp,
        "MAC" => kind == "fgs" ? mac : 0.65,
        "leaf" => kind == "fgs" ? leaf : 6)
    if kind == "fgs"
        merge!(common, Dict("inner" => inner, "max_iterations" => 300,
            "tolerance" => 0.0, "sweep_order" => "lexicographic",
            "cache_leaf_lu" => true, "rlx" => 1.0))
    elseif kind == "krylov_ilu"
        merge!(common, Dict("memory" => 50, "itmax" => 500,
            "atol" => 1e-14, "rtol" => 1e-6, "ilu_leaf" => 10,
            "ilu_MAC" => 1.0, "pattern_entries_per_panel" => 8192,
            "equilibrate" => false, "diagonal_shift" => 0.0,
            "cache_tree" => true, "persistent_plan" => true,
            "cache_nearfield" => false))
    else
        error("Unknown CONFIGS entry: $kind")
    end
    common
end

cold_write_toml(path, data) = open(io -> TOML.print(io, data; sorted=true), path, "w")
cold_id(c) = bytes2hex(sha256(sprint(io -> TOML.print(io, c; sorted=true))))[1:16]

# Pure validation: called before package loading and filesystem mutation.
function cold_check_config(c; selected=false)
    c isa AbstractDict || error("Configuration must be a table")
    get(c, "kind", nothing) in ("fgs", "krylov_ilu") || error("Invalid kind")
    haskey(COLD_SEEDS, get(c, "rung", nothing)) || error("Invalid rung")
    seed = cold_seed(c["rung"], c["kind"])
    Set(keys(c)) ⊆ union(Set(keys(seed)), Set(["diagnostic", "chunks", "dagteam_precision"])) || error("Unknown configuration field")
    Set(keys(seed)) ⊆ Set(keys(c)) || error("Missing configuration field")
    for (key, default) in seed
        value = c[key]
        valid = default isa Bool ? value isa Bool :
            default isa Integer ? value isa Integer && !(value isa Bool) :
            default isa Real ? value isa Real && !(value isa Bool) && isfinite(value) :
            value isa AbstractString
        valid || error("Invalid type/value for $key")
    end
    get(c, "diagnostic", false) isa Bool || error("diagnostic must be Boolean")
    for key in ("P", "leaf", "inner", "max_iterations", "memory", "itmax", "ilu_leaf", "pattern_entries_per_panel", "chunks")
        haskey(c, key) && !(0 < c[key] <= typemax(Int) ÷ 419276) && error("Invalid $key range")
    end
    0 < c["MAC"] <= 1 || error("Invalid MAC")
    if c["kind"] == "fgs"
        c["sweep_order"] in ("lexicographic", "colored", "chunked", "dagteam") || error("Invalid sweep_order")
        if c["sweep_order"] == "chunked"
            # cold_make applies the same 64 default; an explicit key is still
            # required to be a positive integer
            chunks = get(c, "chunks", 64)
            (chunks isa Integer && !(chunks isa Bool) && chunks >= 1) || error("Invalid chunks")
        else
            haskey(c, "chunks") && error("chunks requires sweep_order=chunked")
        end
        if c["sweep_order"] == "dagteam"
            # cold_make applies the same f64 default; an explicit key must name
            # a supported precision mode
            get(c, "dagteam_precision", "f64") in ("f64", "f32conv", "f32full") || error("Invalid dagteam_precision")
        else
            haskey(c, "dagteam_precision") && error("dagteam_precision requires sweep_order=dagteam")
        end
        0 < c["rlx"] < 2 || error("Invalid rlx")
        c["tolerance"] >= 0 || error("Invalid tolerance")
        selected && c["tolerance"] <= 0 && error("Selected FGS requires calibrated positive tolerance")
    else
        0 < c["ilu_MAC"] <= 1 || error("Invalid ilu_MAC")
        c["atol"] >= 0 && c["rtol"] > 0 || error("Invalid Krylov tolerance")
        c["diagonal_shift"] >= 0 || error("Invalid diagonal_shift")
        c["cache_nearfield"] && !c["cache_tree"] && error("cache_nearfield requires cache_tree")
        c["persistent_plan"] && !c["cache_tree"] && error("persistent_plan requires cache_tree")
    end
    return c
end

# Explicit one-factor-at-a-time screen roster, e.g. SCREEN_SET="inner:1,2,3,5"
# or "MAC:0.3,0.5;leaf:25,50,200". Values parse with the seed field's type so
# integer/Boolean axes cannot promote to Float64; the seed itself always runs.
function cold_screen_axes(spec, seed)
    axes = Pair{String,Vector}[]
    for entry in split(spec, ';')
        parts = split(entry, ':')
        length(parts) == 2 || error("SCREEN_SET entries must be key:v1,v2,...")
        key = String(parts[1])
        haskey(seed, key) || error("SCREEN_SET key $key not in $(seed["kind"]) seed")
        default = seed[key]
        values = map(split(parts[2], ',')) do v
            s = String(strip(v))
            default isa Bool ? parse(Bool, s) :
                default isa Integer ? parse(Int, s) :
                default isa Real ? parse(Float64, s) : s
        end
        push!(axes, key => values)
    end
    return Tuple(axes)
end

function cold_configs(rung, kinds, stage; file="")
    screen_set = get(ENV, "SCREEN_SET", "")
    base_file = get(ENV, "SCREEN_BASE_FILE", "")
    if !isempty(base_file)
        stage == "screen" || error("SCREEN_BASE_FILE requires STAGE=screen")
        isempty(file) || error("SCREEN_BASE_FILE conflicts with CONFIG_FILE")
        isempty(screen_set) && error("SCREEN_BASE_FILE requires explicit SCREEN_SET")
        isabspath(base_file) && isfile(base_file) || error("SCREEN_BASE_FILE must be an existing absolute file")
    end
    !isempty(screen_set) && stage != "screen" && error("SCREEN_SET requires STAGE=screen")
    !isempty(screen_set) && !isempty(file) && error("SCREEN_SET conflicts with CONFIG_FILE")
    if !isempty(file)
        doc = TOML.parsefile(file)
        haskey(doc, "configs") && Set(keys(doc)) != Set(["configs"]) && error("Unknown CONFIG_FILE field")
        configs = haskey(doc, "configs") ? doc["configs"] : [doc]
        configs isa AbstractVector && !isempty(configs) || error("configs must be a nonempty array")
        foreach(c -> cold_check_config(c; selected=true), configs)
        all(c -> c["rung"] == rung, configs) || error("CONFIG_FILE rung mismatch")
        return filter(c -> c["kind"] in kinds, configs)
    end
    stage == "verify" && error("STAGE=verify requires CONFIG_FILE")
    seeds = if isempty(base_file)
        [cold_seed(rung, kind) for kind in kinds]
    else
        withenv("SCREEN_BASE_FILE" => nothing, "SCREEN_SET" => nothing) do
            cold_configs(rung, kinds, "verify"; file=base_file)
        end
    end
    configs = Dict{String,Any}[]
    for original in seeds
        seed = copy(original)
        kind = seed["kind"]
        # A saved winner supplies settings, never a tolerance for its neighbors.
        !isempty(base_file) && kind == "fgs" && (seed["tolerance"] = 0.0)
        push!(configs, seed)
        stage == "screen" || continue
        # Bounded one-factor-at-a-time neighbors; do not conflate thread effects.
        # Tuples preserve each axis's element types (numeric vector promotion
        # would turn integer and Boolean Krylov settings into Float64 values).
        axes = !isempty(screen_set) ? cold_screen_axes(screen_set, seed) :
            kind == "fgs" ? (
            "P" => [seed["P"]-2, seed["P"]+2],
            "MAC" => [seed["MAC"]-0.1, seed["MAC"]+0.1],
            "leaf" => [max(10, seed["leaf"]÷2), 2seed["leaf"]],
            "inner" => [max(1, seed["inner"]÷2), 2seed["inner"]],
            "sweep_order" => ["chunked"]) : (
            "P" => [seed["P"]-1, seed["P"]+1],
            "MAC" => [0.6, 0.7], "leaf" => [4, 8], "memory" => [100],
            "cache_nearfield" => [true])
        for (key, values) in axes, value in values
            c = copy(seed); c[key] = value; push!(configs, c)
        end
        if kind == "krylov_ilu" && isempty(screen_set)
            for leaf in (5,10,20), mac in (0.8,1.0)
                c = copy(seed); c["ilu_leaf"] = leaf; c["ilu_MAC"] = mac
                push!(configs, c)
            end
            # Uncached arm is explicitly diagnostic; it cannot win selection.
            c = copy(seed); c["cache_tree"] = false; c["persistent_plan"] = false
            c["diagnostic"] = true; push!(configs, c)
        end
    end
    foreach(cold_check_config, configs)
    unique(cold_id, configs)
end

function cold_preflight(; profile=false)
    stage = get(ENV, "STAGE", "baseline")
    stage in ("baseline", "screen", "verify") || error("Invalid STAGE")
    rung = get(ENV, "RUNG", "")
    haskey(COLD_SEEDS, rung) || error("Cold investigation supports R1–R4")
    for key in ("OUTDIR", "BENCH_CASE_ROOT")
        haskey(ENV, key) && isabspath(ENV[key]) || error("$key must be an explicit absolute path")
    end
    out, case = normpath(ENV["OUTDIR"]), normpath(ENV["BENCH_CASE_ROOT"])
    # Resolve existing ancestors so symlink aliases cannot overlap.
    canonical(p) = ispath(p) ? realpath(p) : joinpath(canonical(dirname(p)), basename(p))
    out, case = canonical(out), canonical(case)
    (out == case || startswith(out, case*"/") || startswith(case, out*"/")) && error("Output and fixture paths overlap")
    (ispath(out) || islink(out)) && error("OUTDIR must be a new generation")
    (ispath(case) || islink(case)) && (!isdir(case) || !isempty(readdir(case))) && error("BENCH_CASE_ROOT must be fresh/empty")
    for path in (out, case)
        parent = dirname(path)
        while !ispath(parent) && !islink(parent)
            parent = dirname(parent)
        end
        isdir(parent) || error("Path ancestor is not a directory")
        # access(W_OK|X_OK), without creating a probe file.
        Sys.isunix() && ccall(:access, Cint, (Cstring, Cint), parent, 3) != 0 && error("Path parent is not writable")
    end
    for key in ("CACHE_B", "SKIP_B")
        get(ENV, key, "0") == "0" || error("Cold investigation requires $key=0")
    end
    get(ENV, "FLOWPANEL_FILAMENT_REG", "") == "linegauss" || error("Set FLOWPANEL_FILAMENT_REG=linegauss")
    mode = get(ENV, "THREADING_MODE", "single")
    mode in ("single", "multi") || error("Invalid THREADING_MODE")
    get(ENV,"KNOBS_MODE",mode) in ("single","multi") || error("Invalid KNOBS_MODE")
    get(ENV,"PER_RUNG_DIR","0") in ("0","1") || error("Invalid PER_RUNG_DIR")
    parse(Int,get(ENV,"K_REPS","1")) > 0 || error("Invalid K_REPS")
    get(ENV, "COLD_PREPARED_ONLY", "0") in ("0", "1") || error("Invalid COLD_PREPARED_ONLY")
    parse(Int, get(ENV, "COLD_PROFILE_REPS", "1")) > 0 || error("Invalid COLD_PROFILE_REPS")
    parse(Int, get(ENV, "COLD_MIN_REPS", "1")) > 0 || error("Invalid COLD_MIN_REPS")
    jt = parse(Int, get(ENV, "EXPECT_JULIA_THREADS", "1"))
    jt > 0 && Threads.nthreads() == jt || error("Julia thread request mismatch")
    bt = parse(Int, get(ENV, "BENCH_BLAS_THREADS", "1"))
    bt > 0 || error("Invalid BENCH_BLAS_THREADS")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "BLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")
        parse(Int, get(ENV, key, "0")) == bt || error("Set $key to requested BLAS count before startup")
    end
    memory = parse(Float64, get(ENV, "MEMORY_GIB", "500"))*2.0^30
    isfinite(memory) && 0 < memory <= 500*2.0^30 || error("MEMORY_GIB must be in (0,500]")
    kinds = split(get(ENV, "CONFIGS", "fgs:krylov_ilu"), r"[:,]")
    all(k -> k in ("fgs", "krylov_ilu"), kinds) || error("Invalid CONFIGS")
    file = get(ENV, "CONFIG_FILE", "")
    profile && isempty(file) && error("Profiling requires selected CONFIG_FILE")
    !isempty(file) && (!isabspath(file) || !isfile(file)) && error("CONFIG_FILE must be an existing absolute file")
    configs = cold_configs(rung, kinds, stage; file)
    isempty(configs) && error("No matching configurations")
    mesh = Dict("R1"=>"23_73", "R2"=>"33_105", "R3"=>"45_145", "R4"=>"65_209")[rung]
    isfile(joinpath(@__DIR__, "..", "examples", "data", "dji9443_20260813_$(mesh)_capped_captess4.msh")) || error("Missing rotor mesh")
    pins = get(ENV, "CAMPAIGN_PINS", "")
    if !isempty(pins)
        isabspath(pins) && isfile(pins) || error("CAMPAIGN_PINS must be an existing absolute file")
        doc = TOML.parsefile(pins)
        for name in ("FLOWPanel", "FastMultipole", "FLOWVPM")
            pin = doc["packages"][name]
            all(k -> get(pin,k,nothing) isa String, ("path","tag","sha")) || error("Invalid campaign pin")
            isabspath(pin["path"]) && isdir(pin["path"]) || error("Invalid pinned path")
        end
    end
    return (; configs, out, case, stage, memory, file)
end

function cold_initialize!(; profile=false)
    inputs = cold_preflight(; profile)
    (; configs, out, case, stage, memory, file) = inputs
    mkpath(dirname(out)); mkdir(out)
    mkpath(case)
    global cold_memory_limit = memory
    global cold_selected_file = file
    !isempty(file) && cp(file, joinpath(out, "input_config.toml"))
    t0 = time_ns()
    Base.include(@__MODULE__, joinpath(@__DIR__, "common.jl"))
    Base.invokelatest(cold_packages) # Provenance gate before fixture work.
    Base.include(@__MODULE__, joinpath(@__DIR__, "phase1_case.jl"))
    Base.invokelatest(cold_assert_threads)
    global cold_fixture_seconds = (time_ns()-t0)/1e9
    Base.invokelatest(cold_provenance, out, stage)
    return configs, out, stage
end

"Reject thread drift after fixture loading and at each measurement boundary."
function cold_assert_threads()
    Threads.nthreads() == banner.julia_threads || error("Julia thread count drifted")
    actual = LinearAlgebra.BLAS.get_num_threads()
    actual == banner.blas_threads || error(
        "BLAS thread count drifted from $(banner.blas_threads) to $actual; " *
        "investigate the reset before benchmarking")
    return nothing
end

function cold_packages()
    expected = isempty(get(ENV, "CAMPAIGN_PINS", "")) ? nothing :
        TOML.parsefile(ENV["CAMPAIGN_PINS"])["packages"]
    packages = Dict{String,Any}()
    for (name, mod) in (("FLOWPanel", pnl), ("FastMultipole", pnl.FastMultipole), ("FLOWVPM", pnl.FLOWVPM))
        path = pkgdir(mod)
        pin = expected === nothing ? nothing : expected[name]
        if pin !== nothing && get(pin, "deployment", "git_worktree") == "rsync"
            manifest = pin["content_manifest"]
            isabspath(manifest) && isfile(manifest) || error("Missing rsync content manifest for $name")
            bytes2hex(sha256(read(manifest))) == pin["content_manifest_sha256"] ||
                error("Rsync content manifest hash mismatch for $name")
            success(Cmd(`sha256sum --quiet -c $manifest`; dir=path)) ||
                error("Rsync deployed content mismatch for $name")
            packages[name] = Dict("path" => path, "sha" => pin["sha"],
                "status" => "rsync_content_verified", "tags" => [pin["tag"]],
                "deployment" => "rsync",
                "content_manifest" => manifest,
                "content_manifest_sha256" => pin["content_manifest_sha256"])
        else
            packages[name] = Dict("path" => path,
                "sha" => readchomp(`git -C $path rev-parse HEAD`),
                "status" => readchomp(`git -C $path status --porcelain`),
                "tags" => split(readchomp(`git -C $path tag --points-at HEAD`), '\n'),
                "deployment" => "git_worktree")
        end
    end
    if expected !== nothing
        for (name, facts) in packages
            pin = expected[name]
            path, tag = facts["path"], pin["tag"]
            realpath(path) == realpath(pin["path"]) || error("Loaded $name outside pinned deployment")
            if facts["deployment"] == "git_worktree"
                facts["sha"] == pin["sha"] && isempty(facts["status"]) || error("Dirty or wrong $name commit")
                readchomp(`git -C $path cat-file -t refs/tags/$tag`) == "tag" || error("Campaign tag must be annotated")
                revision = tag * "^{commit}"
                readchomp(`git -C $path rev-parse $revision`) == facts["sha"] || error("Execution tag mismatch")
                isfile(joinpath(path,".git")) || error("$name is not a git worktree")
            end
        end
    end
    return packages
end

function cold_provenance(out, stage)
    packages = cold_packages()
    project = Base.active_project()
    manifest = joinpath(dirname(project), "Manifest.toml")
    isfile(manifest) || error("Campaign requires Manifest.toml")
    cp(project, joinpath(out,"Project.toml"))
    cp(manifest, joinpath(out,"Manifest.toml"))
    pins = get(ENV,"CAMPAIGN_PINS", "")
    !isempty(pins) && cp(pins, joinpath(out,"campaign_pins.toml"))
    base_file = get(ENV, "SCREEN_BASE_FILE", "")
    !isempty(base_file) && cp(base_file, joinpath(out, "screen_bases.toml"))
    affinity = Sys.islinux() ? read("/proc/self/status", String) : "unavailable on this OS"
    cold_write_toml(joinpath(out, "provenance.toml"), Dict(
        "manifest_sha256" => bytes2hex(sha256(read(manifest))),
        "selected_config_sha256" => isempty(cold_selected_file) ? "generated" : bytes2hex(sha256(read(cold_selected_file))),
        "requested_julia_threads" => parse(Int,ENV["EXPECT_JULIA_THREADS"]),
        "requested_blas_threads" => parse(Int,ENV["BENCH_BLAS_THREADS"]),
        "thread_environment" => Dict(k=>get(ENV,k,"") for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","BLAS_NUM_THREADS")),
        "prepared_only" => get(ENV, "COLD_PREPARED_ONLY", "0") == "1",
        "profile_reps" => parse(Int, get(ENV, "COLD_PROFILE_REPS", "1")),
        "minimum_reps" => parse(Int, get(ENV, "COLD_MIN_REPS", "1")),
        "screen_set" => get(ENV, "SCREEN_SET", ""),
        "screen_base_sha256" => isempty(base_file) ? "rung_seed" : bytes2hex(sha256(read(base_file))),
        "schema" => 2, "stage" => stage, "packages" => packages,
        "julia_version" => string(VERSION), "julia_threads" => Threads.nthreads(),
        "blas_threads" => LinearAlgebra.BLAS.get_num_threads(),
        "blas" => sprint(show, LinearAlgebra.BLAS.get_config()),
        "cpu" => sprint(show, Sys.cpu_info()), "affinity" => affinity,
        "hostname" => gethostname(), "date" => time_string(),
        "filament_reg" => string(pnl.FILAMENT_REGULARIZATION[]),
        "fixture_sha256" => bytes2hex(sha256(read(msh_file))),
        "rhs_sha256" => bytes2hex(sha256(reinterpret(UInt8, b))),
        "fixture_seconds_including_load" => cold_fixture_seconds,
        "direct_rhs_seconds" => t_rhs_assembly, "memory_ceiling_bytes" => cold_memory_limit,
        "timing_scope" => "frozen _solve!; reset and BC diagnostics excluded; no formulation subsolves"))
end

function cold_make(c; history=false)
    if c["kind"] == "fgs"
        solver = pnl.FGSSolver(rotor; expansion_order=c["P"], multipole_acceptance=c["MAC"],
            leaf_size=c["leaf"], inner_iterations=c["inner"], max_iterations=c["max_iterations"],
            tolerance=c["tolerance"], rlx=c["rlx"], shrink=true, recenter=false,
            reverse_pass=false, cache_leaf_lu=c["cache_leaf_lu"],
            sweep_order=Symbol(c["sweep_order"]), chunks=get(c, "chunks", 64),
            dagteam_precision=Symbol(get(c, "dagteam_precision", "f64")), verbose=false,
            project_solution=false, solution_history_length=0)
        cold_assert_threads()
        return solver
    end
    P = pnl.ILUPreconditioner(rotor; leaf_size=c["ilu_leaf"],
        multipole_acceptance=c["ilu_MAC"], max_pattern_entries=c["pattern_entries_per_panel"]*rotor.ncells,
        equilibrate=c["equilibrate"], diagonal_shift=c["diagonal_shift"])
    backend = pnl.FastMultipoleBackend(; expansion_order=c["P"],
        multipole_acceptance=c["MAC"], leaf_size=c["leaf"])
    solver = pnl.KrylovSolver(rotor; method=:gmres, backend, preconditioner=P,
        itmax=c["itmax"], atol=c["atol"], rtol=c["rtol"], memory=c["memory"],
        warmstart=false, record_history=history, cache_tree=c["cache_tree"],
        cache_nearfield=c["cache_nearfield"], persistent_plan=c["persistent_plan"],
        nearfield_cache_max_bytes=floor(Int, cold_memory_limit),
        nearfield_cache_donor=Ref{Any}(nothing))
    @assert solver.kop.nearfield_cache_donor[] === nothing
    @assert solver.kop.plan_slot[] === nothing
    cold_assert_threads()
    return solver
end

function cold_reset!(solver)
    reset_cold!()
    rotor.core_size = rotor.core_size_panel
    solver.niter = 0; solver.solved = false
    pnl.begin_step_solution!(solver)
    if solver isa pnl.KrylovSolver
        solver.have_x_prev = false; solver.x_history_nsaved = 0
        fill!(solver.workspace.x, 0); fill!(solver.rhs, 0)
        fill!(solver.x_prev, 0); fill!(solver.x_history, 0); fill!(solver.x0_scratch, 0)
        pnl.reset!(solver.history)
        @assert all(iszero, solver.workspace.x) && !solver.have_x_prev
        @assert solver.x_history_nsaved == 0 && isempty(solver.history.iter)
    else
        solver.solution_history_nsaved = 0
        fill!(solver.solution_history, 0)
    end
    @assert all(iszero, view(rotor.strength, :, 2))
end

function cold_acceptance(fmm_rel, certified, direct_rel, delta, rung)
    direct = isfinite(direct_rel)
    evaluator = certified ? "certified_fmm" : direct && rung in ("R1", "R2", "R4") ? "direct_fallback" : "uncertified_fmm"
    authoritative_rel = evaluator == "direct_fallback" ? direct_rel : fmm_rel
    accepted = evaluator != "uncertified_fmm" && isfinite(authoritative_rel) && authoritative_rel <= COLD_TARGET
    if certified && direct
        accepted &= direct_rel <= COLD_TARGET && isfinite(delta) && delta <= 1e-7
    end
    return (; authoritative_evaluator=evaluator, authoritative_rel_l2=authoritative_rel, accepted)
end

function cold_validate(x; crosscheck=false)
    rotor.velocity .= frozen_velocity
    phi = zeros(rotor.ncells)
    e = bc_error!(rotor, x; rms_b, target_rel=COLD_TARGET, safety=0.1,
        max_expansion_order=20, multipole_acceptance=0.5, leaf_size=20, phi_out=phi)
    direct_rel = NaN; delta = NaN; direct_seconds = 0.0
    if crosscheck || (!e.error_success && rung in ("R1", "R2", "R4"))
        direct_phi = similar(phi)
        d = bc_error!(rotor, x; rms_b, backend=:direct, phi_out=direct_phi)
        direct_rel = d.rel_l2; direct_seconds = d.t_eval
        delta = norm(phi-direct_phi)/sqrt(rotor.ncells)/rms_b
    end
    acceptance = cold_acceptance(e.rel_l2, e.error_success, direct_rel, delta, rung)
    crosscheck && !isfinite(direct_rel) && error("Direct crosscheck returned nonfinite residual")
    cold_assert_threads()
    (; fmm_rel_l2=e.rel_l2, fmm_rel_max=e.rel_max, fmm_certified=e.error_success,
        fmm_seconds=e.t_eval, epsilon_requested=e.epsilon_requested,
        direct_rel_l2=direct_rel, evaluator_delta=delta, direct_seconds, acceptance...)
end

cold_require(ok, message) = ok ? nothing : error(message)
cold_calibrate_selected!(c, dir; selected=false) = selected ? cold_check_config(c; selected=true) : cold_calibrate!(c, dir)

function cold_csv(path, rows)
    isempty(rows) && return
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ','))
        for row in rows
            println(io, join((_csv_cell(getproperty(row,k)) for k in names), ','))
        end
    end
end

# Diagnostic staircase only: snapshot overhead is excluded from all timing metrics.
function cold_calibrate!(c, dir)
    if c["kind"] == "fgs"
        c["tolerance"] > 0 && return
        probe = copy(c); probe["tolerance"] = 0.01COLD_TARGET*rms_b
        reset_cold!(); solver = cold_make(probe)
        snapshots = Vector{Float64}[]; residuals = Float64[]; iterations = Int[]
        callback = (iteration, residual) -> begin
            pnl.FastMultipole.buffer_to_system_strength!((rotor,), solver.fgs.source_tree)
            push!(snapshots, copy(rotor.strength[:,2]))
            push!(residuals, residual); push!(iterations, iteration)
            nothing
        end
        cold_reset!(solver); pnl._solve!(rotor, solver; callback)
        cold_assert_threads()
        final = copy(rotor.strength[:,2])
        if isempty(snapshots) || final != snapshots[end]
            push!(snapshots, final); push!(residuals, NaN)
            push!(iterations, isempty(iterations) ? 1 : iterations[end]+1)
        end
        rows = NamedTuple[]; crossing = nothing
        for i in eachindex(snapshots)
            e = cold_validate(snapshots[i])
            push!(rows, (; iteration=iterations[i], internal_residual=residuals[i], e...))
            cold_csv(joinpath(dir, "calibration.csv"), rows)
            if crossing === nothing && e.accepted && isfinite(residuals[i]) && residuals[i] > 0
                crossing = i
            elseif crossing !== nothing && e.accepted && 0 < residuals[i] < residuals[crossing]
                c["tolerance"] = sqrt(residuals[i]*residuals[crossing])
                return
            end
        end
        error("FGS staircase has no certified crossing with a decreasing successor; capped candidate")
    else
        rows = NamedTuple[]
        for _ in 1:5
            reset_cold!(); solver = cold_make(c); cold_reset!(solver)
            pnl._solve!(rotor, solver)
            cold_assert_threads()
            e = cold_validate(copy(rotor.strength[:,2]))
            push!(rows, (; rtol=c["rtol"], iterations=solver.niter, solved=solver.solved, e...))
            cold_csv(joinpath(dir, "calibration.csv"), rows)
            solver.solved && e.accepted && return
            c["rtol"] *= 0.1
        end
        error("ILU-GMRES calibration failed after five tolerances")
    end
end

function cold_trial(solver; setup_seconds=0.0, setup_bytes=0, setup_gc=0.0,
                    crosscheck=false, diagnostics=nothing)
    cold_assert_threads()
    cold_reset!(solver)
    timed = @timed pnl._solve!(rotor, solver; diagnostics)
    cold_assert_threads()
    x = copy(rotor.strength[:,2])
    # summarysize traverses the tuple once, deduplicating body/solver references.
    retained = Base.summarysize((rotor, solver))
    e = cold_validate(x; crosscheck)
    row = (; setup_seconds, solve_seconds=timed.time,
        total_seconds=setup_seconds+timed.time, allocated_bytes=setup_bytes+timed.bytes,
        gc_seconds=setup_gc+timed.gctime, retained_bytes=retained,
        process_peak_rss_bytes=Sys.maxrss(), iterations=solver.niter,
        estimated_inner_sweeps=solver isa pnl.FGSSolver ? solver.niter*solver.inner_iterations : -1,
        estimated_fmm_passes=solver isa pnl.FGSSolver ? solver.niter+Int(solver.solved) : -1,
        work_count_kind="estimates; -1 means unavailable", formulation_subsolves=0, solved=solver.solved, e...,
        eligible=solver.solved && e.accepted && retained <= cold_memory_limit && Sys.maxrss() <= cold_memory_limit)
    return row, x
end

function cold_benchmark(c, dir; selected=false)
    cold_calibrate_selected!(c, dir; selected)
    cold_write_toml(joinpath(dir, "config.toml"), c)
    # Compile constructors, lazy plan/cache builds, and solves before sampling.
    reset_cold!(); solver = cold_make(c)
    compile_row, _ = cold_trial(solver)
    cold_csv(joinpath(dir, "compile_validation.csv"), [compile_row])
    cold_require(compile_row.eligible, "Compilation solve failed validation")
    warm, reference = cold_trial(solver; crosscheck=rung in ("R1", "R2"))
    cold_csv(joinpath(dir, "warmup.csv"), [warm])
    warm.eligible || error("Warmup failed accuracy/status/memory gate")
    rows = NamedTuple[]
    # COLD_PREPARED_ONLY=1: the optimization campaign excludes construction and
    # warm starts, so fresh-scope sampling (constructor-dominated) is skipped.
    modes = get(ENV, "COLD_PREPARED_ONLY", "0") == "1" ? ("prepared",) : ("prepared", "fresh")
    for mode in modes
        k = warm.solve_seconds < 60 ? 5 : warm.solve_seconds < 600 ? 3 : 2
        # Fresh setup may dominate; use an excluded, compiled fresh trial to choose k.
        if mode == "fresh"
            solver = nothing; GC.gc(); reset_cold!()
            setup = @timed cold_make(c)
            solver = setup.value
            setup_seconds, setup_bytes, setup_gc = setup.time, setup.bytes, setup.gctime
            setup = nothing
            firstrow, fresh_x = cold_trial(solver; setup_seconds, setup_bytes, setup_gc)
            cold_csv(joinpath(dir, "fresh_warmup.csv"), [firstrow])
            cold_require(firstrow.eligible && norm(fresh_x-reference)/max(norm(reference), eps()) <= 1e-8,
                "Excluded fresh solve failed validation/agreement")
            k = firstrow.total_seconds < 60 ? 5 : firstrow.total_seconds < 600 ? 3 : 2
        end
        k = max(k, parse(Int, get(ENV, "COLD_MIN_REPS", "1")))
        for trial in 1:k
            setup_seconds = 0.0; setup_bytes = 0; setup_gc = 0.0
            if mode == "fresh"
                solver = nothing; GC.gc(); reset_cold!()
                setup = @timed cold_make(c)
                solver = setup.value
                setup_seconds=setup.time; setup_bytes=setup.bytes; setup_gc=setup.gctime
                setup = nothing
            end
            row, x = cold_trial(solver; setup_seconds, setup_bytes, setup_gc)
            agreement = norm(x-reference)/max(norm(reference), eps())
            push!(rows, (; mode, trial, row..., fresh_prepared_relative_delta=agreement,
                repeatable=agreement <= 1e-8))
            cold_csv(joinpath(dir, "trials.csv"), rows)
            cold_require(row.eligible && agreement <= 1e-8, "Trial failed correctness, memory, or repeatability gate")
        end
    end
    summaries = map(modes) do mode
        rr = filter(r -> r.mode == mode, rows)
        times = [r.total_seconds for r in rr]
        (; mode, minimum_seconds=minimum(times), median_seconds=median(times),
            maximum_seconds=maximum(times), spread_seconds=maximum(times)-minimum(times),
            repetitions=length(rr), eligible=all(r -> r.eligible && r.repeatable, rr))
    end
    cold_csv(joinpath(dir, "summary.csv"), summaries)
    # Independent recorder run; enabling history cannot inflate timing trials.
    cold_convergence(c, dir)
    return summaries
end

function cold_convergence(c, dir)
    reset_cold!(); solver = cold_make(c; history=true); cold_reset!(solver)
    rows = NamedTuple[]
    if c["kind"] == "fgs"
        callback = (iteration, residual) -> (push!(rows, (; iteration, residual)); nothing)
        pnl._solve!(rotor, solver; callback)
    else
        pnl._solve!(rotor, solver)
        h = solver.history
        for i in eachindex(h.iter)
            push!(rows, (; iteration=h.iter[i], residual=h.residual_internal[i]))
        end
    end
    cold_csv(joinpath(dir, "convergence.csv"), rows)
    cold_check_solution(solver, dir, "convergence_validation")
end

function cold_check_solution(solver, dir, name; reference=nothing)
    cold_assert_threads()
    x = copy(rotor.strength[:,2])
    e = cold_validate(x; crosscheck=rung in ("R1", "R2"))
    retained = Base.summarysize((rotor, solver))
    rss = Sys.maxrss()
    agreement = reference === nothing ? 0.0 : norm(x-reference)/max(norm(reference), eps())
    eligible = solver.solved && e.accepted && retained <= cold_memory_limit && rss <= cold_memory_limit && agreement <= 1e-8
    cold_csv(joinpath(dir, name*".csv"), [(; solved=solver.solved, retained_bytes=retained,
        process_peak_rss_bytes=rss, relative_solution_delta=agreement, e..., eligible)])
    cold_require(eligible, "$name failed verification")
end

function cold_profile(c, dir; selected=true)
    cold_calibrate_selected!(c, dir; selected)
    cold_write_toml(joinpath(dir, "config.toml"), c)
    reset_cold!(); solver = cold_make(c)
    row, reference = cold_trial(solver; crosscheck=rung in ("R1", "R2"))
    cold_csv(joinpath(dir, "unprofiled_trial.csv"), [row])
    row.eligible || error("Profile candidate failed validation")
    # COLD_PROFILE_REPS accumulates samples over repeated prepared solves for
    # denser attribution; resets stay outside every recorded region and the
    # final solution is validated below.
    Profile.clear()
    for _ in 1:parse(Int, get(ENV, "COLD_PROFILE_REPS", "1"))
        cold_reset!(solver)
        Profile.@profile pnl._solve!(rotor, solver)
    end
    serialize(joinpath(dir, "cpu_profile.jls"), Profile.retrieve())
    for format in (:flat, :tree)
        open(joinpath(dir, "cpu_$(format).txt"), "w") do io
            Profile.print(io; format, C=true, groupby=[:thread, :task])
        end
    end
    cold_check_solution(solver, dir, "cpu_validation"; reference)
    cold_reset!(solver); Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate=0.01 pnl._solve!(rotor, solver)
    allocations = Profile.Allocs.fetch()
    serialize(joinpath(dir, "allocation_profile.jls"), allocations)
    open(joinpath(dir, "allocations.txt"), "w") do io
        println(io, "sample_rate=0.01; sampled_bytes=", sum(a.size for a in allocations.allocs; init=0))
        println(io, "Actual bytes/GC/RSS are from the separate unprofiled trial in unprofiled_trial.csv.")
        for a in sort(allocations.allocs; by=a -> a.size, rev=true)[1:min(end,100)]
            println(io, a.size, '\t', a.stacktrace)
        end
    end
    cold_check_solution(solver, dir, "allocation_validation"; reference)
end

# Constructor/reset integration checks on the actual frozen fixture, before timing.
function cold_smoke(configs, out)
    calibrated = Dict{String,Any}[]
    for config in configs
        c = copy(config)
        dir = joinpath(out,c["kind"]); mkdir(dir)
        try
            cold_calibrate_selected!(c, dir; selected=!isempty(cold_selected_file))
            cold_write_toml(joinpath(dir,"config.toml"),c)
            reset_cold!(); first_solver = cold_make(c)
            reset_cold!(); second_solver = cold_make(c)
            independent = if c["kind"] == "krylov_ilu"
                first_solver.kop.nearfield_cache_donor !== second_solver.kop.nearfield_cache_donor &&
                first_solver.kop.plan_slot !== second_solver.kop.plan_slot
            else
                first_solver.fgs !== second_solver.fgs && first_solver.fgs.source_tree !== second_solver.fgs.source_tree
            end
            cold_require(independent, "Constructors share solver cache state")
            rows = NamedTuple[]; reference = nothing
            for (i, solver) in enumerate((first_solver, first_solver, second_solver))
                # R4 uses certified evaluation, with direct fallback if inconclusive.
                row, x = cold_trial(solver; crosscheck=rung != "R4")
                reference === nothing && (reference = x)
                agreement = norm(x-reference)/max(norm(reference),eps())
                push!(rows,(; trial=i, independent, zero_reset_verified=true, relative_solution_delta=agreement, row...))
                cold_csv(joinpath(dir,"smoke.csv"),rows)
                cold_require(row.eligible && agreement <= 1e-8, "Smoke correctness/memory/agreement failure")
            end
            first_solver = nothing; second_solver = nothing
            push!(calibrated,c)
            cold_write_toml(joinpath(dir,"status.toml"),Dict("status"=>"completed"))
        catch err
            cold_write_toml(joinpath(dir,"status.toml"),Dict("status"=>"failed","error"=>sprint(showerror,err,catch_backtrace())))
            rethrow()
        end
    end
    cold_write_toml(joinpath(out,"selected.toml"),Dict("configs"=>calibrated))
end

function cold_main(; profile=false, smoke=false)
    configs, out, stage = cold_initialize!(; profile)
    if smoke
        Base.invokelatest(cold_smoke, configs, out)
    else
        Base.invokelatest(cold_run, configs, out, stage; profile)
    end
end

function cold_run(configs, out, stage; profile=false)
    selected = Dict{String,Any}[]; best = Dict{String,Float64}()
    failures = 0
    for config in configs
        c = copy(config)
        dir = joinpath(out, rung, c["kind"]*"_"*cold_id(c),
            "j$(Threads.nthreads())_b$(LinearAlgebra.BLAS.get_num_threads())")
        mkpath(dir)
        cold_write_toml(joinpath(dir, "requested_config.toml"), c)
        try
            println("Cold investigation: ", c)
            if profile
                cold_profile(c, dir; selected=!isempty(cold_selected_file))
            else
                summary = cold_benchmark(c, dir; selected=!isempty(cold_selected_file))
                prepared = first(summary)
                if all(s -> s.eligible, summary) && !get(c, "diagnostic", false) &&
                        prepared.minimum_seconds < get(best, c["kind"], Inf)
                    best[c["kind"]] = prepared.minimum_seconds
                    filter!(x -> x["kind"] != c["kind"], selected); push!(selected, copy(c))
                end
            end
            !isempty(cold_selected_file) && cold_require(c == config, "Selected configuration mutated")
            cold_write_toml(joinpath(dir, "status.toml"), Dict("status" => "completed"))
        catch err
            err isa InterruptException && rethrow()
            failures += 1
            cold_write_toml(joinpath(dir, "status.toml"), Dict("status" => "failed", "error" => sprint(showerror, err, catch_backtrace())))
            @error "Recorded failed candidate" dir exception=(err, catch_backtrace())
            # Generated screen candidates record failures and continue (e.g. a
            # roster point whose calibration cannot cross the accuracy gate);
            # baseline/verify and selected executions still stop at the first
            # failed arm.
            if stage == "screen" && isempty(cold_selected_file)
                GC.gc()
                continue
            end
            rethrow()
        end
        GC.gc()
    end
    !profile && isempty(selected) && error("No eligible candidates; inspect status.toml and trials.csv")
    cold_write_toml(joinpath(out, "selected.toml"), Dict("configs" => selected))
    stage == "verify" && failures > 0 && error("Verification has failed candidates")
    return nothing
end
