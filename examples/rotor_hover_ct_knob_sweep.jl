#!/usr/bin/env julia
# CT knob sweep: warm-started perturbation runs on the settled rotor-hover
# baseline (data/rotor_hover_pressure_comparison, RHPC_MESH=40_40, NT=36,
# corrected Pedrizzetti rlxf=0.3). Each scenario resumes from RESTART_STEP
# (default 350 ≈ rev 9.72, freestream fully withdrawn) and continues to
# 14 revs with ONE knob changed, so CT effects read as perturbations on the
# same settled wake.
#
# Knobs under test (see plans/20260702_convergence.md):
#   1. BODY_HESSIAN_TO_PARTICLES / PANEL_WAKE_HESSIAN_TO_PARTICLES — restore
#      surface-induced strain on particles (suspect: instability).
#   2. KERNELOFFSET_TARGETS — surface->particle smoothing radius (WAKE_CORE_SIZE
#      pinned to the baseline 1e-3 so shed-particle cores are NOT co-swept).
#   3. RELAX_RLXF halving ladder + RELAX_FILTER_DOWNSTREAM_R — relaxation as a
#      rate λ=rlxf/dt, extrapolated to λ→0 (artificial-dissipation treatment).
#
# ENV:
#   SCENARIOS  comma list (default "control,bodyhess"); "all" for everything
#   THREADS    julia + BLAS threads per scenario subprocess (default 36)
#   SMOKE      true => ~11-step continuation instead of ~154 (default false)
#   RESTART_STEP (default 350)
#
# Each scenario is an isolated subprocess (instability is an expected OUTCOME);
# output overwrites data/ct_knob_<name>/ in place; log at data/ct_knob_<name>.log.
# Summary (merged across invocations): data/ct_knob_sweep_summary.csv.

import Printf: @printf
import Statistics: mean, std

const RESTART_NAME = "rotor_hover_pressure_comparison"
const RESTART_PATH = joinpath("data", RESTART_NAME)
const RESTART_STEP = parse(Int, get(ENV, "RESTART_STEP", "350"))
const THREADS = get(ENV, "THREADS", "36")
const SMOKE = parse(Bool, get(ENV, "SMOKE", "false"))
const NT = 36
const RPM = 6000.0
const DT = 60 / RPM / NT
const RESTART_REV = RESTART_STEP / NT

# Baseline ENV shared by every scenario. Phases: freestream fully withdrawn by
# rev 9 (< restart rev 9.72), so the continuation is pure hover under any
# baseline phase split; SETTLE_REVS sets total length (2+3+4+settle revs).
base_env() = Dict(
    "RHPC_MESH" => "40_40",
    "NT" => string(NT),
    "RPM" => "6000",
    "FREESTREAM_RAMP_REVS" => "2",
    "FREESTREAM_HOLD_REVS" => "3",
    "FREESTREAM_WITHDRAW_REVS" => "4",
    "SETTLE_REVS" => SMOKE ? "1.05" : "5",   # 10.05 revs (smoke) / 14 revs
    "RESTART_STEP" => string(RESTART_STEP),
    "RESTART_NAME" => RESTART_NAME,
    "RESTART_PATH" => RESTART_PATH,
    "SAVE_VTK" => "true",
    "RUN_KJ" => "false",
    # BLAS threads: the example also calls LinearAlgebra.BLAS.set_num_threads
    # from BLAS_NUM_THREADS/OMP_NUM_THREADS at startup, so this is enforced
    # regardless of BLAS vendor or library load order.
    "OPENBLAS_NUM_THREADS" => THREADS,
    "OMP_NUM_THREADS" => THREADS,
    "BLAS_NUM_THREADS" => THREADS,
    "JULIA_NUM_THREADS" => THREADS,
)

# scenario name => ENV overrides (single knob each; control = pure continuation)
const SCENARIO_DEFS = [
    "control"          => Dict{String,String}(),
    "bodyhess"         => Dict("BODY_HESSIAN_TO_PARTICLES" => "true"),
    # split regularization: body->particle velocity at KERNELOFFSET_TARGETS
    # (1e-3), but its GRADIENT at a 4x larger offset — smooths the constant-
    # doublet-panel |∇U| bumpiness felt by near-blade particles (not strictly
    # physical; stability aid for bodyhess)
    "bodyhess_gradoff" => Dict("BODY_HESSIAN_TO_PARTICLES" => "true",
                               "BODY_GRADIENT_KERNELOFFSET" => "4e-3"),
    "wakerowhess"      => Dict("PANEL_WAKE_HESSIAN_TO_PARTICLES" => "true"),
    "koff_5e-4"        => Dict("KERNELOFFSET_TARGETS" => "5e-4", "WAKE_CORE_SIZE" => "1e-3"),
    "koff_2p5e-4"      => Dict("KERNELOFFSET_TARGETS" => "2.5e-4", "WAKE_CORE_SIZE" => "1e-3"),
    "rlxf_0p15"        => Dict("RELAX_RLXF" => "0.15"),
    "rlxf_0p075"       => Dict("RELAX_RLXF" => "0.075"),
    "relaxfilter_0p5R" => Dict("RELAX_FILTER_DOWNSTREAM_R" => "0.5"),
    # disable particle merging (baseline merges every step at r=0.02R): tests
    # whether merge-induced vorticity diffusion is suppressing CT
    "merge_off"        => Dict("MERGE_PARTICLES" => "false"),
]
const SCENARIO_RLXF = Dict("control" => 0.3, "rlxf_0p15" => 0.15, "rlxf_0p075" => 0.075)

function run_scenario(name, overrides)
    run_name = "ct_knob_$(name)"
    datadir = joinpath("data", run_name)
    logpath = joinpath("data", "$(run_name).log")
    rm(datadir; force=true, recursive=true)

    env = merge(copy(ENV), base_env(), overrides, Dict("RUN_NAME" => run_name))
    cmd = setenv(`julia --project=. -t $(THREADS) examples/rotor_hover_pressure_comparison.jl`, env)
    @printf("=== %s: launching (log: %s) ...\n", name, logpath)
    t0 = time()
    ok = success(pipeline(cmd; stdout=logpath, stderr=logpath))
    elapsed = time() - t0
    @printf("=== %s: %s after %.1f min\n", name, ok ? "finished" : "FAILED/UNSTABLE (see log)", elapsed / 60)
    return (; run_name, datadir, logpath, ok, elapsed)
end

# --- CT CSV parsing ---------------------------------------------------------
function parse_ct(datadir, run_name)
    path = joinpath(datadir, "$(run_name)_CT_vs_rev.csv")
    isfile(path) || return nothing
    header = String[]
    rows = Vector{Float64}[]
    for (i, line) in enumerate(eachline(path))
        if i == 1
            header = split(strip(line), ",")
        else
            push!(rows, [something(tryparse(Float64, x), NaN) for x in split(line, ",")])
        end
    end
    isempty(rows) && return nothing
    M = reduce(hcat, rows)  # ncol × nrow
    return (; header, M)
end

function ct_metrics(ct)
    # continuation rows only (monitor history is zero before the restart step)
    rev = ct.M[2, :]
    cont = findall(r -> r > RESTART_REV + 1e-9, rev)
    isempty(cont) && return nothing
    maxrev = maximum(rev[cont])
    lastrev = [k for k in cont if rev[k] > maxrev - 1.0]
    methods = ct.header[3:end]
    out = Dict{String,NamedTuple}()
    for (j, m) in enumerate(methods)
        v = ct.M[2+j, lastrev]
        vf = filter(isfinite, v)
        out[m] = (; mean=isempty(vf) ? NaN : mean(vf),
                  std=isempty(vf) ? NaN : (length(vf) > 1 ? std(vf) : 0.0),
                  ptp=isempty(vf) ? NaN : maximum(vf) - minimum(vf),
                  nan_frac=1 - length(vf) / max(length(v), 1))
    end
    return (; maxrev, out)
end

function np_stability(datadir, run_name)
    pdir = joinpath(datadir, "$(run_name)_wake1_particles")
    isdir(pdir) || return (; np_ok=false, size_last=0, size_max=0)
    files = filter(f -> endswith(f, ".vtp"), readdir(pdir))
    isempty(files) && return (; np_ok=false, size_last=0, size_max=0)
    steps = [parse(Int, split(f, ".")[end-1]) for f in files]
    sizes = [filesize(joinpath(pdir, f)) for f in files]
    ord = sortperm(steps)
    keep = [k for k in ord if steps[k] > RESTART_STEP]  # continuation only
    isempty(keep) && return (; np_ok=false, size_last=0, size_max=0)
    size_last = sizes[keep[end]]
    size_max = maximum(sizes[keep])
    return (; np_ok=size_last > 0.5 * size_max, size_last, size_max)
end

# --- main --------------------------------------------------------------------
selected = let s = get(ENV, "SCENARIOS", "control,bodyhess")
    lowercase(s) == "all" ? first.(SCENARIO_DEFS) : String.(strip.(split(s, ",")))
end
defs = Dict(SCENARIO_DEFS)
for s in selected
    haskey(defs, s) || error("Unknown scenario '$(s)'. Known: $(join(first.(SCENARIO_DEFS), ", "))")
end

results = NamedTuple[]
for name in selected
    r = run_scenario(name, defs[name])
    ct = parse_ct(r.datadir, r.run_name)
    met = ct === nothing ? nothing : ct_metrics(ct)
    nps = np_stability(r.datadir, r.run_name)
    ct_b = met === nothing ? NaN : met.out["CT_bernoulli"].mean
    ct_b_std = met === nothing ? NaN : met.out["CT_bernoulli"].std
    ct_lv = met === nothing ? NaN : met.out["CT_laplace_lamb"].mean
    stable = r.ok && met !== nothing && isfinite(ct_b) && nps.np_ok
    push!(results, (; scenario=name, status=r.ok ? "completed" : "crashed",
        stable, CT_bernoulli=ct_b, CT_bernoulli_std=ct_b_std,
        CT_laplace_lamb=ct_lv,
        maxrev=met === nothing ? NaN : met.maxrev,
        np_size_last=nps.size_last, np_size_max=nps.size_max,
        elapsed_min=r.elapsed / 60))
end

# merge with previous invocations' summary so the workstation can run scenarios
# incrementally
summary_path = joinpath("data", "ct_knob_sweep_summary.csv")
prev = Dict{String,String}()
if isfile(summary_path)
    lines = readlines(summary_path)
    for line in lines[2:end]
        prev[split(line, ",")[1]] = line
    end
end
names = propertynames(results[1])
for r in results
    prev[r.scenario] = join((getproperty(r, n) for n in names), ",")
end
open(summary_path, "w") do io
    println(io, join(names, ","))
    for k in sort(collect(keys(prev)))
        println(io, prev[k])
    end
end

println("\n=== Sweep summary (last-rev CT, continuation only) ===")
@printf("%-18s %-10s %-7s %12s %12s %12s %10s\n",
        "scenario", "status", "stable", "CT_bern", "CT_bern_std", "CT_lapl_lv", "min")
for r in results
    @printf("%-18s %-10s %-7s %12.5f %12.2e %12.5f %10.1f\n",
            r.scenario, r.status, string(r.stable), r.CT_bernoulli,
            r.CT_bernoulli_std, r.CT_laplace_lamb, r.elapsed_min)
end
ctrl = findfirst(r -> r.scenario == "control", results)
if ctrl !== nothing && isfinite(results[ctrl].CT_bernoulli)
    c = results[ctrl].CT_bernoulli
    println("\nΔCT vs control (Bernoulli):")
    for r in results
        r.scenario == "control" && continue
        isfinite(r.CT_bernoulli) || continue
        @printf("  %-18s %+.5f (%+.2f%%)\n", r.scenario, r.CT_bernoulli - c, 100 * (r.CT_bernoulli - c) / c)
    end
end

# Richardson extrapolation of CT(λ) to λ→0 over the stable rlxf ladder
ladder = [(SCENARIO_RLXF[r.scenario] / DT, r.CT_bernoulli)
          for r in results if haskey(SCENARIO_RLXF, r.scenario) && r.stable && isfinite(r.CT_bernoulli)]
if length(ladder) >= 3
    sort!(ladder; by=first)
    (l1, c1), (l2, c2), (l3, c3) = ladder[1:3]   # λ ascending; ~factor-2 ladder
    # CT(λ) = CT0 + a·λ^p with λ3/λ2 ≈ λ2/λ1 = 2  =>  p from CT differences
    r21 = c2 - c1; r32 = c3 - c2
    if r21 != 0 && sign(r32) == sign(r21)
        p = log(abs(r32 / r21)) / log(l3 / l2)
        ct0 = c1 - r21 / ((l2 / l1)^p - 1)
        @printf("\nRelaxation-rate extrapolation: λ = %s s⁻¹ → CT = %s\n",
                join(round.(first.(ladder[1:3]); digits=1), ", "),
                join(round.(last.(ladder[1:3]); sigdigits=5), ", "))
        @printf("  fitted order p = %.2f;  CT(λ→0) ≈ %.5f\n", p, ct0)
    else
        println("\nRelaxation ladder non-monotone; no λ→0 extrapolation (report raw values).")
    end
elseif !isempty(ladder)
    println("\nRelaxation ladder has $(length(ladder)) stable point(s); need 3 for λ→0 extrapolation.")
end
@printf("\nWrote %s\n", summary_path)
