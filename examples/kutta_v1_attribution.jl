#=##############################################################################
Item 015, Phase 4 V1: wake-attachment / Kutta-closure attribution matrix.

Run one stage at a time:

    KUTTAV1_STAGE=diamond julia --project -t 6 examples/kutta_v1_attribution.jl
    KUTTAV1_STAGE=gate julia --project -t 6 examples/kutta_v1_attribution.jl
    KUTTAV1_STAGE=gate_diagnosis julia --project -t 1 examples/kutta_v1_attribution.jl
    KUTTAV1_STAGE=gate_diagnosis2 julia --project -t 1 examples/kutta_v1_attribution.jl
    KUTTAV1_STAGE=gate_diagnosis2 KUTTAV1_DIAG2_BATTERIES=mesh julia --project -t 1 \
        examples/kutta_v1_attribution.jl
    KUTTAV1_STAGE=ssw0 julia --project -t 6 examples/kutta_v1_attribution.jl
    KUTTAV1_STAGE=ssw_alpha julia --project -t 6 examples/kutta_v1_attribution.jl
    KUTTAV1_STAGE=sweptwing julia --project -t 6 examples/kutta_v1_attribution.jl

Output defaults to data/kutta_v1_attribution. Existing completed stage summaries
are resume checkpoints; set KUTTAV1_FORCE=true to overwrite a stage.
=###############################################################################

import FLOWPanel as pnl
import LinearAlgebra as LA
import Printf: @printf, @sprintf
import Statistics: median, mean

include(joinpath(@__DIR__, "suddenly_started_wing.jl"))

const KUTTAV1_DEFAULT_OUTPUT = joinpath("data", "kutta_v1_attribution")
const KUTTAV1_COMBOS = (
    (label="A_jump", route=:A, pressure=false),
    (label="A_pressure", route=:A, pressure=true),
    (label="B_jump", route=:B, pressure=false),
    (label="B_pressure", route=:B, pressure=true),
)

_kuttav1_bool(name, default=false) =
    lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")

function _csv_value(x)
    if x isa AbstractString || x isa Symbol
        s = replace(string(x), '"' => "\"\"")
        return occursin(r"[,\"]", s) ? "\"$s\"" : s
    elseif x isa Bool
        return string(x)
    elseif x isa Real
        return @sprintf("%.16g", x)
    end
    return _csv_value(string(x))
end

function write_namedtuple_csv(path, rows)
    isempty(rows) && error("refusing to write empty CSV: $path")
    mkpath(dirname(path))
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for row in rows
            println(io, join((_csv_value(getproperty(row, n)) for n in names), ","))
        end
    end
    return path
end

function append_namedtuple_csv(path, rows)
    isempty(rows) && return path
    mkpath(dirname(path))
    exists = isfile(path)
    names = propertynames(first(rows))
    open(path, "a") do io
        exists || println(io, join(string.(names), ","))
        for row in rows
            println(io, join((_csv_value(getproperty(row, n)) for n in names), ","))
        end
    end
    return path
end

function combo_objects(combo, rho; on_failure=:error, pressure_scale=:auto,
        correction_scale=:auto)
    attachment = combo.route == :A ?
        pnl.RigidTransitionAttachment() : pnl.TEAnchoredAttachment()
    closure = if combo.pressure
        pnl.PressureContinuityKutta(pnl.SteadyBernoulliProvider(rho);
            pressure_tolerance=1e-6, correction_tolerance=1e-6,
            on_failure, store_diagnostics=true, pressure_scale, correction_scale)
    else
        pnl.JumpKutta()
    end
    return attachment, closure
end

function paired_te_metrics(body, rho, pressure_scale)
    residuals = Float64[]
    for shedding in body.shedding, i in axes(shedding, 2)
        upper, lower = shedding[1, i], shedding[4, i]
        uu = LA.norm(view(body.velocity, :, upper))
        ul = LA.norm(view(body.velocity, :, lower))
        push!(residuals, rho / 2 * (ul^2 - uu^2) / pressure_scale)
    end
    absres = abs.(residuals)
    return (; max_dp_hat=isempty(absres) ? NaN : maximum(absres),
        median_dp_hat=isempty(absres) ? NaN : median(absres),
        residuals)
end

function te_gamma_rows(case_name, combo_label, body)
    gamma_i = pnl.get_Gammai(body)
    rows = NamedTuple[]
    for shedding in body.shedding, i in axes(shedding, 2)
        upper, lower = shedding[1, i], shedding[4, i]
        y = 0.5 * (body.controlpoints[2, upper] + body.controlpoints[2, lower])
        gamma = body.strength[upper, gamma_i] - body.strength[lower, gamma_i]
        push!(rows, (; case=case_name, combo=combo_label, edge=i, y, gamma))
    end
    sort!(rows; by=r -> r.y)
    return rows
end

function gamma_symmetry(rows)
    n = length(rows)
    n == 0 && return NaN
    maximum(abs(rows[i].gamma - rows[n + 1 - i].gamma) for i in 1:n)
end

function diagnostic_metrics(closure)
    diags = pnl.kutta_diagnostics(closure)
    isempty(diags) && return (; c_inf_scaled=0.0, mean_outer_iterations=0.0,
        body_solves=0, mean_body_solves=0.0, startup_steps=0,
        startup_body_solves=0, backtracks=0, restarts=0,
        max_r_inf_scaled=NaN, status="jump", disposition="jump",
        pressure_scale=NaN)
    cscaled = Float64[]
    for d in diags
        (isempty(d.c) || !isfinite(d.correction_scale)) && continue
        append!(cscaled, (abs(c) / d.correction_scale for c in d.c if isfinite(c)))
    end
    residuals = [d.r_inf_scaled for d in diags if isfinite(d.r_inf_scaled)]
    # A `:startup_jump` step runs no nonlinear iteration by design (Route B's
    # deterministic cold start), so averaging it in dilutes the pressure-solve
    # cost by an amount that depends only on run length. Report it separately.
    startup = [d for d in diags if d.disposition === :startup_jump]
    solved = [d for d in diags if d.disposition !== :startup_jump]
    nsolved = max(length(solved), 1)
    return (; c_inf_scaled=isempty(cscaled) ? 0.0 : maximum(cscaled),
        mean_outer_iterations=sum(d.outer_iterations for d in solved; init=0) / nsolved,
        body_solves=sum(d.body_solves for d in diags),
        mean_body_solves=sum(d.body_solves for d in solved; init=0) / nsolved,
        startup_steps=length(startup),
        startup_body_solves=sum(d.body_solves for d in startup; init=0),
        backtracks=sum(d.backtracks for d in diags),
        restarts=sum(d.restarts for d in diags),
        max_r_inf_scaled=isempty(residuals) ? NaN : maximum(residuals),
        status=join(unique(string(d.status) for d in diags), "|"),
        disposition=join(unique(string(d.disposition) for d in diags), "|"),
        pressure_scale=diags[end].pressure_scale)
end

function make_diamond_body(; nspan=1, thick=0.06, das=0.3)
    ys = range(0, 1; length=nspan + 1)
    nodes = Float64[]
    for y in ys
        append!(nodes, (0.0, y, 0.0))
        append!(nodes, (0.5, y, thick))
        append!(nodes, (1.0, y, 0.0))
        append!(nodes, (0.5, y, -thick))
    end
    nodes = reshape(nodes, 3, :)
    idx(j, k) = (j - 1) * 4 + k
    cells = Int[]
    for j in 1:nspan
        le1, up1, te1, lo1 = idx(j, 1), idx(j, 2), idx(j, 3), idx(j, 4)
        le2, up2, te2, lo2 = idx(j + 1, 1), idx(j + 1, 2), idx(j + 1, 3), idx(j + 1, 4)
        append!(cells, (le1, up1, up2)); append!(cells, (le1, up2, le2))
        append!(cells, (up1, te1, te2)); append!(cells, (up1, te2, up2))
        append!(cells, (le1, le2, lo2)); append!(cells, (le1, lo2, lo1))
        append!(cells, (lo1, lo2, te2)); append!(cells, (lo1, te2, te1))
    end
    cells = reshape(cells, 3, :)
    shedding = pnl.calc_shedding_from_seed(nodes, cells, idx(1, 3), idx(2, 3))
    bodytype = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.VortexRing}}
    body = bodytype(nodes, cells, [shedding]; check_mesh=false, watertight=false,
        semiinfinite_wake=false)
    for Das in body.Das
        Das .= repeat([das, 0.0, 0.0], 1, size(Das, 2))
    end
    return body
end

function run_diamond_combo(combo)
    rho = 1.2
    body = make_diamond_body()
    wake = pnl.PanelWake(body; nwakerows=8)
    frames = pnl.ReferenceFrame(body; name="vehicle")
    solver = pnl.Backslash(body)
    attachment, closure = combo_objects(combo, rho)
    u = [1.0, 0.0, 0.3]
    t_range = collect(range(0.0; step=0.05, length=4))
    Dhat = u / LA.norm(u)
    Lhat = LA.cross(Dhat, [0.0, 1.0, 0.0])
    pressure = pnl.PressureBernoulli(rho; backend=pnl.DirectBackend(),
        correct_kuttacondition=false)
    force = pnl.ForceMonitor(length(t_range), 1; i_frame=-1,
        normalization=pnl.NoNormalization(), correct_kuttacondition=false,
        verbose=false)
    Uinf(t) = u
    maneuver!(frames, systems, wakes, t) = nothing
    elapsed = @elapsed pnl.simulate!((body,), (wake,), frames, maneuver!, Uinf,
        t_range; body_solvers=(solver,), backend=pnl.DirectBackend(), path=nothing,
        monitors=(pressure, force), grad_mu_options=(; basis=:tri),
        wake_attachment=attachment,
        kutta_closure=closure)
    dm = diagnostic_metrics(closure)
    # One fixed normalization for every cell and both columns: the kutta closure
    # freezes its own per-step `pressure_scale`, which differs cell to cell, so
    # using it here made the four Δp̂ columns mutually incomparable.
    q = 0.5 * rho * LA.dot(u, u)
    te = paired_te_metrics(body, rho, q)
    gamma = te_gamma_rows("diamond", combo.label, body)
    CL = LA.dot(force.force[:, end], Lhat) / q
    summary = (; case="diamond", combo=combo.label, route=combo.route,
        closure=combo.pressure ? "pressure" : "jump", CL,
        dp_hat_scale=q, kutta_pressure_scale=dm.pressure_scale,
        max_dp_hat_pre_from_jump_run=NaN, max_dp_hat_post=te.max_dp_hat,
        median_dp_hat_post=te.median_dp_hat, dm.c_inf_scaled,
        dm.mean_outer_iterations, dm.body_solves, dm.mean_body_solves,
        dm.startup_steps, dm.startup_body_solves,
        dm.backtracks, dm.restarts, dm.max_r_inf_scaled, dm.status,
        dm.disposition, elapsed_s=elapsed, gamma_symmetry=gamma_symmetry(gamma))
    return (; summary, gamma)
end

_ssw_default_case_name(config) = config.alpha_deg == 0 ? "ssw0" : "ssw_alpha"

function _ssw_combo_run(config, combo; on_failure=:error,
        grad_mu_options=SSW_GRAD_MU_OPTIONS, settle_window=nothing,
        case_name=_ssw_default_case_name(config),
        formulation=pnl.VelocityThroughSources(), extra_monitors=())
    sim = prepare_suddenly_started_wing(config)
    attachment, closure = combo_objects(combo, config.rho; on_failure,
        pressure_scale=0.5 * config.rho * config.U^2,
        correction_scale=config.U * config.c)
    # Extra monitors are appended, never inserted: `sim.monitors` has an
    # ordering contract (pressure -> force -> spanwise ...) that `audit_monitors`
    # enforces, and read-only diagnostics have nothing to provide upstream.
    monitors = isempty(extra_monitors) ? sim.monitors :
        (sim.monitors..., extra_monitors...)
    elapsed = @elapsed pnl.simulate!((sim.wing,), (sim.wake,), sim.frames,
        sim.maneuver!, sim.Uinf, sim.t_range;
        body_solvers=(sim.solver,), backend=sim.backend, monitors=monitors,
        path=nothing, set_Das_eta_freestream=NaN,
        grad_mu_options, formulation, wake_attachment=attachment,
        kutta_closure=closure, verbose=false)
    q = 0.5 * config.rho * config.U^2
    qS = q * config.AR * config.c^2
    CL = vec(sim.lift_hat' * sim.force.force) ./ qS
    CD = vec(sim.drag_hat' * sim.force.force) ./ qS
    t_star = collect(sim.t_range) .* config.U ./ config.c
    dm = diagnostic_metrics(closure)
    te = paired_te_metrics(sim.wing, config.rho, q)
    gamma = te_gamma_rows(case_name, combo.label, sim.wing)
    history = [(; case=case_name, combo=combo.label, step=i - 1,
        t_star=t_star[i], CL=CL[i], CD=CD[i]) for i in eachindex(sim.t_range)]
    if settle_window === nothing
        settled_start = max(2, length(CL) - max(3, round(Int, 1 / config.dt_star)) + 1)
        settled = settled_start:length(CL)
    else
        settled = [i for i in 2:length(CL)
            if settle_window[1] - 1e-9 <= t_star[i] <= settle_window[2] + 1e-9]
        isempty(settled) && error("settle window $settle_window contains no samples")
    end
    settled_CL = mean(CL[settled])
    summary = (; case=case_name,
        combo=combo.label, route=combo.route,
        closure=combo.pressure ? "pressure" : "jump", CL=settled_CL,
        settled_samples=length(settled),
        settled_CL_spread=maximum(CL[settled]) - minimum(CL[settled]),
        dp_hat_scale=q, kutta_pressure_scale=dm.pressure_scale,
        max_dp_hat_pre_from_jump_run=NaN, max_dp_hat_post=te.max_dp_hat,
        median_dp_hat_post=te.median_dp_hat, dm.c_inf_scaled,
        dm.mean_outer_iterations, dm.body_solves, dm.mean_body_solves,
        dm.startup_steps, dm.startup_body_solves,
        dm.backtracks, dm.restarts, dm.max_r_inf_scaled, dm.status,
        dm.disposition, elapsed_s=elapsed, gamma_symmetry=gamma_symmetry(gamma))
    return (; summary, gamma, history, body=sim.wing, CL, CD, t_star)
end

function run_ssw_combo(config, combo; optargs...)
    try
        return _ssw_combo_run(config, combo; optargs...)
    catch err
        combo.pressure || rethrow()
        @warn "pressure cell failed; rerunning with on_failure=:jump" combo=combo.label exception=err
        mapped = _ssw_combo_run(config, combo; on_failure=:jump, optargs...)
        summary = merge(mapped.summary,
            (; status="pressure_error_mapped|" * mapped.summary.status,
                disposition="pressure_error_mapped|" * mapped.summary.disposition))
        return merge(mapped, (; summary))
    end
end

#--- formulation sanity gate ---------------------------------------------------
#
# Every number this gate's diagnosis needs is produced here and lands in CSV.
# The 2026-07-29 version compared two runs at an unsettled t*<=2, recorded no
# strength-level comparison, and its supporting numbers were hand-collected;
# the resulting localization could not be reproduced or falsified.
#
# Discriminators, one CSV row per (formulation x grad_mu basis) cell:
#   * caps      — the leading suspect. `simplewing`/`suddenly_started_wing_mesh`
#                 loft OPEN tips, and the interior-Dirichlet Green identity
#                 assumes a closed surface (the solver warns for the Neumann
#                 analogue but is silent here). Flat caps are reported to match
#                 Neumann-no-caps best. Uncapped Dirichlet is retained as a row
#                 so the fraction of the old 14-26% gap the caps close is
#                 measured, not asserted.
#   * grad_mu basis — `:tri` on an all-triangle mesh is the repo's known
#                 divergent configuration (BRAINSTORM chord-divergence sweep),
#                 and prior sweptwing-mirror work traced a CL gap to exactly
#                 this velocity reconstruction. `:tri_robust` separates
#                 reconstruction sensitivity from the solve.
# Quad bases are N/A on an all-triangle mesh.

const KUTTAV1_GATE_COMMON = (; AR=6.0, alpha_deg=5.0, n_span=12, n_airfoil=21,
    t_end_star=8.0, dt_star=0.125, backend_kind=:direct,
    save_vtk=false, verbose=false)
const KUTTAV1_GATE_SETTLE = (7.0, 8.0)
const KUTTAV1_GATE_CELLS = (
    (label="neumann_uncapped", bodytype=:neumann, caps=:none),
    (label="dirichlet_capped", bodytype=:dirichlet, caps=:flat),
    (label="dirichlet_uncapped", bodytype=:dirichlet, caps=:none),
)
# `grad_mu` reconstruction variants. `:tri` is the legacy SSW setting
# (SSW_GRAD_MU_OPTIONS) and the repo's known-divergent configuration on an
# all-triangle mesh; `:tri_robust` is the same basis with the robust stencil
# (`tri_robust=true`), which is the direct discriminator. Quad bases are N/A
# here — the SSW mesh has no quads to recover.
const KUTTAV1_GATE_BASES = (
    (name=:tri, options=(; basis=:tri)),
    (name=:tri_robust, options=(; basis=:tri, tri_robust=true)),
)
const KUTTAV1_GATE_REFERENCE_BASIS = :tri   # the legacy SSW_GRAD_MU_OPTIONS basis

"""
    kj_lift_coefficient(config, body) -> CL

Kutta–Joukowski lift coefficient from the paired trailing-edge bound
circulation: `CL = 2 ∫Γ(y) dy / (U·S)`. Γ is the TE strength jump (Δγ for the
Neumann vortex-ring body, Δμ for the Dirichlet source+doublet body, both being
the bound circulation). Spanwise weights are the edges' Voronoi widths, with
the outer boundaries at the tips, which is exact for the uniform SSW mesh.

The sign follows the strength convention, so compare magnitudes.
"""
function kj_lift_coefficient(config, body)
    rows = te_gamma_rows("kj", "kj", body)      # sorted by y
    isempty(rows) && return NaN
    ys = [r.y for r in rows]
    b = config.AR * config.c
    n = length(ys)
    bounds = Vector{Float64}(undef, n + 1)
    bounds[1], bounds[n + 1] = -b / 2, b / 2
    for i in 2:n
        bounds[i] = 0.5 * (ys[i - 1] + ys[i])
    end
    circulation = sum(rows[i].gamma * (bounds[i + 1] - bounds[i]) for i in 1:n)
    return 2 * circulation / (config.U * config.AR * config.c^2)
end

"Pearson correlation; NaN when either input has no variance."
function _pearson(x, y)
    length(x) == length(y) || throw(ArgumentError("length mismatch"))
    length(x) >= 2 || return NaN
    xb, yb = mean(x), mean(y)
    sxy = sum((x .- xb) .* (y .- yb))
    sxx = sum(abs2, x .- xb)
    syy = sum(abs2, y .- yb)
    (sxx > 0 && syy > 0) || return NaN
    return sxy / sqrt(sxx * syy)
end

"""
Compare two spanwise TE Δstrength distributions defined on the same shedding
edges: count sign flips and measure shape agreement. Both formulations share the
mesh and the shedding pairing, so the edges correspond one-to-one.
"""
function spanwise_agreement(reference_rows, test_rows)
    length(reference_rows) == length(test_rows) ||
        return (; sign_flips=-1, shape_correlation=NaN, max_rel_diff=NaN,
            scale_ratio=NaN, edges=0, matched=false)
    g_ref = [r.gamma for r in reference_rows]
    g_test = [r.gamma for r in test_rows]
    flips = count(i -> g_ref[i] != 0 && g_test[i] != 0 &&
        sign(g_ref[i]) != sign(g_test[i]), eachindex(g_ref))
    scale = maximum(abs, g_ref)
    # Least-squares level factor: with a near-unit correlation, `scale_ratio`
    # is the single number separating "same spanwise shape, different
    # circulation level" from a genuine shape disagreement.
    denom = sum(abs2, g_ref)
    return (; sign_flips=flips, shape_correlation=_pearson(g_ref, g_test),
        max_rel_diff=scale > 0 ? maximum(abs.(g_test .- g_ref)) / scale : NaN,
        scale_ratio=denom > 0 ? sum(g_ref .* g_test) / denom : NaN,
        edges=length(g_ref), matched=true)
end

"""
Read-only monitor recording, at every step and for every shedding edge, the
attached-strip / first-wake-row strength jump. The strip carries the *current*
bound circulation (it is part of the LHS); wake row 1 carries the *previous*
step's, planted by `shed_wake!`. Monitors run before `propagate!`/`shed_wake!`,
so this samples the jump the solve actually saw.
"""
mutable struct StripJumpRecorder
    rows::Vector{NamedTuple}
end
StripJumpRecorder() = StripJumpRecorder(NamedTuple[])

function (rec::StripJumpRecorder)(systems, wakes, frames, uinf, i_step, dt)
    body = systems[1]
    wake = wakes[1]
    panel_wake = wake isa pnl.PanelParticleWake ? wake.panel_wake : wake
    for i_surf in eachindex(body.shedding)
        for i in axes(body.shedding[i_surf], 2)
            gamma = pnl._get_wakestrength_Gamma(body, i, i_surf)
            row1 = panel_wake.strength[i_surf][1, 1, i]
            push!(rec.rows, (; step=i_step, t=i_step * dt, surf=i_surf, edge=i,
                gamma, row1, jump=gamma - row1))
        end
    end
    return nothing
end

"""
Control assertion (plan §5.6): `Das` is freestream-parallel and frozen. Compared
against the value `prepare_suddenly_started_wing` installs, on the body as it
stands after the final step — `Das` is only ever rotated with the body, and this
body never rotates, so any drift is a bug.
"""
function _das_control(config, body)
    drag_hat, _ = _ssw_directions(config)
    dt = step(ssw_time_range(config))
    expected = config.eta * dt * config.U * collect(drag_hat)
    uhat = collect(drag_hat)
    par_err = 0.0
    drift = 0.0
    for Das in body.Das, j in axes(Das, 2)
        d = view(Das, :, j)
        nd = LA.norm(d)
        nd > 0 && (par_err = max(par_err, LA.norm(LA.cross(collect(d), uhat)) / nd))
        drift = max(drift, LA.norm(collect(d) .- expected) / max(LA.norm(expected), eps()))
    end
    return (; das_parallel_err=par_err, das_drift=drift)
end

"""
    _run_gate_cell(cell, basis; battery, tag) -> (; row, gamma, steady_gamma, history)

Run one configuration and emit the uniform battery row schema. `cell` carries
`label`, `bodytype`, `caps`, and may override any of `eta`,
`kerneloffset_over_c`, `n_span`, `n_airfoil`; anything absent falls back to
`KUTTAV1_GATE_COMMON`/`SSWConfig` defaults, so the gate's own call sites are
unaffected. Every knob is recorded in the row, so a battery row and a gate row
are directly comparable.

A cell that throws is recorded with NaN metrics and a `status` naming the
failure rather than being silently dropped or retried at a different setting —
`kerneloffset = 0` may be genuinely singular for the vortex-ring body.
"""
function _run_gate_cell(cell, basis; battery="gate", tag=nothing,
        record_strip_jump=false)
    combo = first(KUTTAV1_COMBOS)        # A/jump: no attachment or closure runtime
    get_or(name, default) = hasproperty(cell, name) ? getproperty(cell, name) : default
    formulation = get_or(:formulation, pnl.VelocityThroughSources())
    # Guard, never coerce: `GreenReconstruction` is Dirichlet-only (NK=2,
    # DBC=true) and `simulate!` would reject a Neumann body deeper in, after the
    # setup cost. Fail here, with the cell named.
    formulation isa pnl.VelocityThroughSources || cell.bodytype === :dirichlet ||
        error("formulation $(nameof(typeof(formulation))) requires " *
              "bodytype=:dirichlet; cell $(cell.label) is $(cell.bodytype)")
    overrides = (; eta=get_or(:eta, SSWConfig().eta),
        kerneloffset_over_c=get_or(:kerneloffset_over_c,
            SSWConfig().kerneloffset_over_c),
        n_span=get_or(:n_span, KUTTAV1_GATE_COMMON.n_span),
        n_airfoil=get_or(:n_airfoil, KUTTAV1_GATE_COMMON.n_airfoil),
        freestream_convection=get_or(:freestream_convection,
            SSWConfig().freestream_convection))
    config = SSWConfig(; KUTTAV1_GATE_COMMON..., bodytype=cell.bodytype,
        caps=cell.caps, overrides...)
    label = something(tag, "$(cell.label)_$(basis.name)")
    println("$battery: $label")
    knobs = (; cell=label, battery, bodytype=cell.bodytype, caps=cell.caps,
        basis=basis.name, config.eta, config.kerneloffset_over_c,
        config.n_span, config.n_airfoil,
        formulation=string(nameof(typeof(formulation))),
        config.freestream_convection)

    recorder = record_strip_jump ? StripJumpRecorder() : nothing
    extra_monitors = isnothing(recorder) ? () : (recorder,)
    local unsteady, steady
    try
        unsteady = run_ssw_combo(config, combo; grad_mu_options=basis.options,
            settle_window=KUTTAV1_GATE_SETTLE, case_name=label,
            formulation, extra_monitors)
        steady = _ssw_steady_cl(config, pnl.DirectBackend();
            grad_mu_options=basis.options)
    catch err
        @warn "battery cell failed; recorded as NaN, not retried" cell=label exception=err
        nan_row = (; knobs..., settled_CL=NaN,
            settle_lo=KUTTAV1_GATE_SETTLE[1], settle_hi=KUTTAV1_GATE_SETTLE[2],
            settled_samples=0, settled_CL_spread=NaN, CL_final=NaN,
            steady_semiinfinite_CL=NaN, steady_semiinfinite_CD=NaN,
            CL_kj_unsteady=NaN, CL_kj_steady=NaN, ncells=0, watertight=false,
            shedding_edges=0, max_dp_hat_post=NaN, median_dp_hat_post=NaN,
            gamma_symmetry=NaN, das_parallel_err=NaN, das_drift=NaN,
            status="failed:$(nameof(typeof(err)))",
            elapsed_s=NaN)
        return (; row=nan_row, gamma=NamedTuple[], steady_gamma=NamedTuple[],
            history=NamedTuple[], strip=NamedTuple[])
    end

    watertight, _ = pnl.iswatertight(unsteady.body.nodes, unsteady.body.cells)
    row = (; knobs...,
        settled_CL=unsteady.summary.CL,
        settle_lo=KUTTAV1_GATE_SETTLE[1], settle_hi=KUTTAV1_GATE_SETTLE[2],
        unsteady.summary.settled_samples, unsteady.summary.settled_CL_spread,
        CL_final=unsteady.CL[end],
        steady_semiinfinite_CL=steady.cl, steady_semiinfinite_CD=steady.cd,
        CL_kj_unsteady=kj_lift_coefficient(config, unsteady.body),
        CL_kj_steady=kj_lift_coefficient(config, steady.wing),
        ncells=unsteady.body.ncells, watertight,
        shedding_edges=sum(size(s, 2) for s in unsteady.body.shedding),
        unsteady.summary.max_dp_hat_post, unsteady.summary.median_dp_hat_post,
        unsteady.summary.gamma_symmetry,
        _das_control(config, unsteady.body)...,
        # A cell can fail without throwing: `kerneloffset = 0` makes the
        # vortex-ring near-wake self-influence singular and NaN propagates
        # silently through the solve. Non-finite results must say so.
        status=(isfinite(unsteady.summary.CL) &&
                isfinite(kj_lift_coefficient(config, unsteady.body)) &&
                isfinite(steady.cl)) ? "ok" : "nonfinite",
        elapsed_s=unsteady.summary.elapsed_s + steady.elapsed)
    steady_gamma = te_gamma_rows(label, "steady", steady.wing)
    strip = isnothing(recorder) ? NamedTuple[] :
        _settled_strip_jump(knobs, config, recorder)
    return (; row, gamma=unsteady.gamma, steady_gamma, history=unsteady.history,
        strip)
end

"""
Reduce a `StripJumpRecorder` to one row per shedding edge, averaged over the
settled window. Normalization is one number per cell — the settled-window
maximum |Γ| — so every edge of a cell is scaled identically and the near-tip
edges (where Γ→0) cannot manufacture a large ratio.
"""
function _settled_strip_jump(knobs, config, recorder)
    lo, hi = KUTTAV1_GATE_SETTLE
    sel = [r for r in recorder.rows
        if lo - 1e-9 <= r.t * config.U / config.c <= hi + 1e-9]
    isempty(sel) && return NamedTuple[]
    gscale = maximum(abs(r.gamma) for r in sel)
    rows = NamedTuple[]
    for key in sort(unique((r.surf, r.edge) for r in sel))
        rs = [r for r in sel if (r.surf, r.edge) == key]
        absjump = [abs(r.jump) for r in rs]
        push!(rows, (; knobs..., surf=key[1], edge=key[2], samples=length(rs),
            gamma_mean=mean(r.gamma for r in rs),
            jump_mean=mean(r.jump for r in rs),
            abs_jump_mean=mean(absjump), max_abs_jump=maximum(absjump),
            gamma_scale=gscale,
            mean_jump_over_gamma=gscale > 0 ? mean(absjump) / gscale : NaN,
            max_jump_over_gamma=gscale > 0 ? maximum(absjump) / gscale : NaN))
    end
    return rows
end

"""
    formulation_gap(neumann_row, dirichlet_row) -> NamedTuple

Split a Dirichlet-vs-Neumann gap into a circulation part and a
force-reconstruction part. `CL_kj_*` depends only on the TE strength jump, so it
is blind to `grad_mu` and to pressure integration; the `ForceMonitor` CL depends
on both. `*_pressure_over_kj` is each body's own pressure CL divided by its own
KJ CL — those need not be 1 (KJ is a span integral, not an identity with the
surface pressure integral), so it is the *difference* between the two ratios
that carries information.

This is the gate's own localization arithmetic, factored out so the gate and the
diagnosis battery cannot drift apart.
"""
function formulation_gap(n, d)
    reldiff(a, b) = abs(a - b) / max(abs(b), eps(Float64))
    return (; unsteady_settled_gap=reldiff(d.settled_CL, n.settled_CL),
        unsteady_kj_gap=reldiff(abs(d.CL_kj_unsteady), abs(n.CL_kj_unsteady)),
        steady_pressure_gap=reldiff(d.steady_semiinfinite_CL,
            n.steady_semiinfinite_CL),
        steady_kj_gap=reldiff(abs(d.CL_kj_steady), abs(n.CL_kj_steady)),
        neumann_pressure_over_kj=n.settled_CL / abs(n.CL_kj_unsteady),
        dirichlet_pressure_over_kj=d.settled_CL / abs(d.CL_kj_unsteady))
end

"""
    gap_factorization(n, d) -> NamedTuple

Split the settled-CL ratio into the two factors it is the *exact* product of:

    CL_d / CL_n  =  (|KJ_d| / |KJ_n|)  ×  ((CL_d/|KJ_d|) / (CL_n/|KJ_n|))
                     circulation             force reconstruction

The first factor is blind to `grad_mu` and to pressure integration; the second
is exactly the part that is not. `steady_reconstruction_ratio` is the same
reconstruction quantity built from the semi-infinite *steady* solve, which no
formulation kwarg can reach (`_ssw_steady_cl` runs `steady!` on a
`semiinfinite_wake=true` body, and `GreenReconstruction` is rejected there by
`_validate_formulation_common`). That makes it the comparator: an unsteady
reconstruction ratio that equals it is the steady discretization error and
nothing more; one that departs from it is something the marching wake added.

`factorization_residual` is an identity and must land at round-off — it is a
free harness check on this whole decomposition.
"""
function gap_factorization(n, d)
    per_kj(cl, kj) = cl / abs(kj)
    circ = abs(d.CL_kj_unsteady) / abs(n.CL_kj_unsteady)
    recon = per_kj(d.settled_CL, d.CL_kj_unsteady) /
        per_kj(n.settled_CL, n.CL_kj_unsteady)
    settled_ratio = d.settled_CL / n.settled_CL
    return (; circulation_ratio=circ, reconstruction_ratio=recon,
        steady_circulation_ratio=abs(d.CL_kj_steady) / abs(n.CL_kj_steady),
        steady_reconstruction_ratio=
            per_kj(d.steady_semiinfinite_CL, d.CL_kj_steady) /
            per_kj(n.steady_semiinfinite_CL, n.CL_kj_steady),
        settled_ratio, factorization_residual=circ * recon - settled_ratio)
end

"""
    run_gate(output) -> (; gate_row, rows)

The formulation sanity gate. Pre-registered criterion, fixed before the runs:
the settled-CL gap between flat-capped Dirichlet and uncapped Neumann (both at
the legacy `:tri` basis) is at most 5%, AND their spanwise TE Δstrength
distributions agree qualitatively — no sign flips, shape correlation ≥ 0.98.
The Kutta–Joukowski numbers are diagnostics only, and are admitted as evidence
only if the extractor first reproduces the Neumann steady pressure CL to 10%.
"""
function run_gate(output)
    cells = NamedTuple[]
    gamma = NamedTuple[]
    steady_gamma = NamedTuple[]
    history = NamedTuple[]
    for cell in KUTTAV1_GATE_CELLS, basis in KUTTAV1_GATE_BASES
        result = _run_gate_cell(cell, basis)
        push!(cells, result.row)
        append!(gamma, result.gamma)
        append!(steady_gamma, result.steady_gamma)
        append!(history, result.history)
    end
    write_namedtuple_csv(joinpath(output, "gate_summary.csv"), cells)
    write_namedtuple_csv(joinpath(output, "gate_cl_history.csv"), history)
    write_namedtuple_csv(joinpath(output, "gate_spanwise_gamma.csv"), gamma)
    write_namedtuple_csv(joinpath(output, "gate_spanwise_gamma_steady.csv"),
        steady_gamma)

    find(label) = only(r for r in cells
        if r.cell == "$(label)_$(KUTTAV1_GATE_REFERENCE_BASIS)")
    neu = find("neumann_uncapped")
    dir = find("dirichlet_capped")
    dir_open = find("dirichlet_uncapped")

    # KJ extractor validation on the Neumann steady case, against its own
    # pressure-integrated CL and the lifting-line anchor 2*pi*a*AR/(AR+2).
    alpha = deg2rad(KUTTAV1_GATE_COMMON.alpha_deg)
    AR = KUTTAV1_GATE_COMMON.AR
    cl_ll = 2pi * alpha * AR / (AR + 2)
    kj_err = abs(abs(neu.CL_kj_steady) - abs(neu.steady_semiinfinite_CL)) /
        max(abs(neu.steady_semiinfinite_CL), eps(Float64))
    kj_valid = isfinite(kj_err) && kj_err <= 0.10
    kj = (; check="kj_extractor_vs_neumann_steady_pressure",
        CL_kj_steady=neu.CL_kj_steady,
        CL_pressure_steady=neu.steady_semiinfinite_CL,
        relative_error=kj_err, CL_lifting_line=cl_ll, tolerance=0.10,
        valid=kj_valid,
        note=kj_valid ? "KJ admitted as gate evidence" :
            "KJ extractor cannot reproduce the Neumann steady pressure CL; " *
            "excluded from gate evidence")
    write_namedtuple_csv(joinpath(output, "gate_kj_validation.csv"), [kj])

    # Localization table, computed by `formulation_gap` so the gate and the
    # diagnosis battery can never drift apart.
    at(label, basis) = only(r for r in cells if r.cell == "$(label)_$(basis)")
    localization = NamedTuple[]
    for basis in KUTTAV1_GATE_BASES, capcase in ("dirichlet_capped", "dirichlet_uncapped")
        push!(localization, (; comparison="$(capcase)_vs_neumann_uncapped",
            basis=basis.name,
            formulation_gap(at("neumann_uncapped", basis.name),
                at(capcase, basis.name))...))
    end
    # The basis rows compare a cell against itself at the other reconstruction,
    # so the "Dirichlet" slot is the tri_robust run; the pressure_over_kj
    # columns are not meaningful in that pairing and are dropped to NaN.
    for c in KUTTAV1_GATE_CELLS
        g = formulation_gap(at(c.label, :tri), at(c.label, :tri_robust))
        push!(localization, (; comparison="$(c.label)_basis_tri_robust_vs_tri",
            basis=:both, g.unsteady_settled_gap, g.unsteady_kj_gap,
            g.steady_pressure_gap, g.steady_kj_gap,
            neumann_pressure_over_kj=NaN, dirichlet_pressure_over_kj=NaN))
    end
    write_namedtuple_csv(joinpath(output, "gate_localization.csv"), localization)

    gsel(label) = [r for r in gamma if r.case == "$(label)_$(KUTTAV1_GATE_REFERENCE_BASIS)"]
    agree = spanwise_agreement(gsel("neumann_uncapped"), gsel("dirichlet_capped"))
    agree_open = spanwise_agreement(gsel("neumann_uncapped"),
        gsel("dirichlet_uncapped"))

    gap(a, b) = abs(a - b) / max(abs(b), eps(Float64))
    capped_gap = gap(dir.settled_CL, neu.settled_CL)
    uncapped_gap = gap(dir_open.settled_CL, neu.settled_CL)
    basis_sensitivity = maximum(
        gap(only(r.settled_CL for r in cells if r.cell == "$(c.label)_tri_robust"),
            only(r.settled_CL for r in cells if r.cell == "$(c.label)_tri"))
        for c in KUTTAV1_GATE_CELLS)

    pass = isfinite(capped_gap) && capped_gap <= 0.05 &&
        agree.matched && agree.sign_flips == 0 &&
        isfinite(agree.shape_correlation) && agree.shape_correlation >= 0.98
    gate_row = (; gate="dirichlet_neumann_A_jump_settled_capped",
        basis=KUTTAV1_GATE_REFERENCE_BASIS,
        neumann_CL=neu.settled_CL, dirichlet_capped_CL=dir.settled_CL,
        dirichlet_uncapped_CL=dir_open.settled_CL,
        relative_gap=capped_gap, relative_gap_uncapped=uncapped_gap,
        gap_closed_by_caps=uncapped_gap > 0 ?
            1 - capped_gap / uncapped_gap : NaN,
        max_basis_sensitivity=basis_sensitivity,
        agree.sign_flips, agree.shape_correlation, agree.max_rel_diff,
        agree.scale_ratio,
        uncapped_shape_correlation=agree_open.shape_correlation,
        uncapped_scale_ratio=agree_open.scale_ratio,
        kj_extractor_valid=kj_valid, pass)
    write_namedtuple_csv(joinpath(output, "ssw_sanity_gate.csv"), [gate_row])

    @printf("gate: Neumann %.6f, Dirichlet capped %.6f (gap %.2f%%), Dirichlet uncapped %.6f (gap %.2f%%)\n",
        neu.settled_CL, dir.settled_CL, 100 * capped_gap, dir_open.settled_CL,
        100 * uncapped_gap)
    @printf("gate: spanwise sign flips %d, shape correlation %.5f, max basis sensitivity %.2f%%, KJ extractor valid: %s\n",
        agree.sign_flips, agree.shape_correlation, 100 * basis_sensitivity,
        kj_valid)
    return (; gate_row, rows=cells)
end

#--- gate diagnosis battery ----------------------------------------------------
#
# Reporting only. This stage does NOT touch `ssw_sanity_gate.csv`, the pass
# criterion, or `KUTTAV1_GATE_REFERENCE_BASIS`, and it gates nothing. It answers
# the three questions the gate itself left open, given that the gap decomposes
# into a circulation factor (KJ vs KJ, 12.8% unsteady but only 3.0% with a
# semi-infinite wake) and a force-reconstruction factor:
#
#   A  eta ladder  — is the circulation factor a `Das = eta*dt*U` artifact? It
#                    is absent in the semi-infinite steady solve and appears
#                    only under the marching PanelWake, which is what makes the
#                    near-wake offset the leading suspect.
#   B  mesh ladder — does the gap refine away, and is a 5% criterion attainable
#                    at all for two formulations on this geometry? The whole
#                    gate ran at one resolution.
#   C  kerneloffset ladder — does the Dirichlet G regularization matter?
#                    `SSWConfig` uses 1e-6*c while the DJI convention is an
#                    unregularized Dirichlet G.

const KUTTAV1_DIAG_PAIR = (
    (label="neumann_uncapped", bodytype=:neumann, caps=:none),
    (label="dirichlet_capped", bodytype=:dirichlet, caps=:flat),
)
# Extended past eta=1 to test whether the KJ gap flattens at the semi-infinite
# value (3.0%) or keeps falling; item 014 found rotor CT flat for eta >= 1.
const KUTTAV1_DIAG_ETAS = (0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0)
const KUTTAV1_DIAG_MESHES = ((n_span=12, n_airfoil=21), (n_span=24, n_airfoil=41))
# 0.0 is the DJI unregularized convention; it may be singular for the
# vortex-ring body, in which case the cell is recorded as failed, never retried
# at a different offset.
const KUTTAV1_DIAG_OFFSETS = (1e-4, 1e-6, 1e-8, 0.0)

_diag_basis(name) = only(b for b in KUTTAV1_GATE_BASES if b.name == name)

"Run one battery cell and file its row, spanwise gammas and history."
function _diag_push!(store, base, overrides, basis; battery, tag,
        record_strip_jump=false)
    cell = merge(base, overrides)
    result = _run_gate_cell(cell, basis; battery, tag, record_strip_jump)
    push!(store.rows, result.row)
    append!(store.gamma, result.gamma)
    append!(store.history, result.history)
    hasproperty(store, :strip) && append!(store.strip, result.strip)
    return result.row
end

"""
    run_gate_diagnosis(output; batteries)

`batteries` selects sub-batteries (default all three), overridable with
`KUTTAV1_DIAG_BATTERIES=eta,kerneloffset`. Single-threaded runs are slow enough
that the refined-mesh ladder deserves to be launched separately from the cheap
ladders, and a partial battery is still a complete answer to its own question.
Output files are suffixed when a subset is selected so a partial run cannot
overwrite a full one.
"""
function run_gate_diagnosis(output;
        batteries=Tuple(Symbol(strip(x)) for x in
            split(get(ENV, "KUTTAV1_DIAG_BATTERIES", "eta,mesh,kerneloffset"), ',')))
    for b in batteries
        b in (:eta, :mesh, :kerneloffset) ||
            error("unknown battery $b; use eta, mesh, kerneloffset")
    end
    suffix = length(batteries) == 3 ? "" : "_" * join(sort(collect(string.(batteries))), "-")
    store = (; rows=NamedTuple[], gamma=NamedTuple[], history=NamedTuple[])
    pairs = NamedTuple[]        # (battery, axis value, basis, neumann tag, dirichlet tag)

    # --- A: eta ladder -----------------------------------------------------
    for eta in (:eta in batteries ? KUTTAV1_DIAG_ETAS : ()), basis in KUTTAV1_GATE_BASES
        tags = Dict{Symbol, String}()
        for base in KUTTAV1_DIAG_PAIR
            tag = "eta$(_ssw_num_tag(eta))_$(base.label)_$(basis.name)"
            _diag_push!(store, base, (; eta), basis; battery="etaladder", tag)
            tags[base.bodytype] = tag
        end
        push!(pairs, (; battery="etaladder", axis="eta", axis_value=eta,
            basis=basis.name, neumann=tags[:neumann], dirichlet=tags[:dirichlet]))
    end

    # --- B: mesh ladder ----------------------------------------------------
    for mesh in (:mesh in batteries ? KUTTAV1_DIAG_MESHES : ()), basis in KUTTAV1_GATE_BASES
        tags = Dict{Symbol, String}()
        for base in KUTTAV1_DIAG_PAIR
            tag = "ns$(mesh.n_span)na$(mesh.n_airfoil)_$(base.label)_$(basis.name)"
            _diag_push!(store, base, mesh, basis; battery="meshladder", tag)
            tags[base.bodytype] = tag
        end
        push!(pairs, (; battery="meshladder", axis="ncells_proxy",
            axis_value=Float64(2 * (mesh.n_airfoil - 1) * mesh.n_span),
            basis=basis.name, neumann=tags[:neumann], dirichlet=tags[:dirichlet]))
    end

    # --- C: kerneloffset ladder (reference basis only; KJ is basis-blind) ---
    basis = _diag_basis(KUTTAV1_GATE_REFERENCE_BASIS)
    for offset in (:kerneloffset in batteries ? KUTTAV1_DIAG_OFFSETS : ())
        tags = Dict{Symbol, String}()
        for base in KUTTAV1_DIAG_PAIR
            tag = "ko$(_ssw_num_tag(offset))_$(base.label)_$(basis.name)"
            _diag_push!(store, base, (; kerneloffset_over_c=offset), basis;
                battery="kerneloffset", tag)
            tags[base.bodytype] = tag
        end
        push!(pairs, (; battery="kerneloffset", axis="kerneloffset_over_c",
            axis_value=offset, basis=basis.name,
            neumann=tags[:neumann], dirichlet=tags[:dirichlet]))
    end

    write_namedtuple_csv(joinpath(output, "gate_diagnosis_summary$(suffix).csv"), store.rows)
    write_namedtuple_csv(joinpath(output, "gate_diagnosis_cl_history$(suffix).csv"),
        store.history)
    write_namedtuple_csv(joinpath(output, "gate_diagnosis_spanwise_gamma$(suffix).csv"),
        store.gamma)

    byname = Dict(r.cell => r for r in store.rows)
    gsel(tag) = [r for r in store.gamma if r.case == tag]
    gaps = NamedTuple[]
    for p in pairs
        n, d = byname[p.neumann], byname[p.dirichlet]
        agree = spanwise_agreement(gsel(p.neumann), gsel(p.dirichlet))
        push!(gaps, (; p.battery, p.axis, p.axis_value, p.basis,
            neumann_CL=n.settled_CL, dirichlet_CL=d.settled_CL,
            formulation_gap(n, d)..., agree.shape_correlation, agree.scale_ratio,
            agree.sign_flips, ncells_neumann=n.ncells, ncells_dirichlet=d.ncells,
            status=n.status == "ok" && d.status == "ok" ? "ok" :
                "$(n.status)|$(d.status)"))
    end
    write_namedtuple_csv(joinpath(output, "gate_diagnosis_gaps$(suffix).csv"), gaps)

    for g in gaps
        @printf("%-13s %-18s %10.5g  basis=%-10s  settled gap %7.2f%%  KJ gap %7.2f%%  %s\n",
            g.battery, g.axis, g.axis_value, string(g.basis),
            100 * g.unsteady_settled_gap, 100 * g.unsteady_kj_gap, g.status)
    end
    return (; rows=store.rows, gaps)
end

#--- gate diagnosis II: small-Das attribution ----------------------------------
#
# Reporting only; gates nothing, is gated by nothing, and does not touch
# `ssw_sanity_gate.csv` or the 5% criterion.
#
# The first battery established that *both* formulations lose bound circulation
# as Das -> 0, the Dirichlet body ~3.5x more. Agreement at LARGE Das is
# tautological (Das is the finite stand-in for the semi-infinite row); the
# differential sensitivity at SMALL Das is the finding that needs a cause.
#
#   T1 formulation axis — for the Dirichlet body the shed wake reaches the solve
#          only through `set_strengths!`, as an equivalent source sheet
#          sigma = -u_wake.n; its direct contribution to the interior-potential
#          condition is dropped (the exact RHS term is -phi_wake). The Neumann
#          body has no analogue: u.n *is* its boundary condition, unconverted.
#          `GreenReconstruction` transmits the wake as a reconstructed potential
#          trace instead, so if that substitution is the cause, the
#          `dirichlet_green` arm must collapse onto the Neumann curve.
#   T2 geometry control — `freestream_convection=true` makes both wakes the same
#          straight uniform sheet, removing the rollup confound: the two runs
#          currently convect each row with their own body's induced velocity.
#   T3 strip/row-1 strength jump — the attached strip carries the current bound
#          circulation while wake row 1 carries the previous step's. Recorded on
#          every T1/T2 run at no extra cost.

const KUTTAV1_DIAG2_ARMS = (
    (arm="neumann", label="neumann_uncapped", bodytype=:neumann, caps=:none,
        formulation=pnl.VelocityThroughSources()),
    (arm="dirichlet_vts", label="dirichlet_capped", bodytype=:dirichlet,
        caps=:flat, formulation=pnl.VelocityThroughSources()),
    (arm="dirichlet_green", label="dirichlet_green", bodytype=:dirichlet,
        caps=:flat, formulation=pnl.GreenReconstruction()),
)
const KUTTAV1_DIAG2_ETAS = (0.0625, 0.125, 0.25, 0.5, 1.0)
const KUTTAV1_DIAG2_FC_ETAS = (0.0625, 0.25, 1.0)
# T4 (below) — the mesh ladder under each formulation. Three rungs, not two: a
# two-rung agreement is not convergence (the DJI Phase 2c lesson). `n_span` must
# be even and `n_airfoil` odd and >= 21 (`suddenly_started_wing_mesh` /
# `naca0012_contour`); the largest rung is ~4436 cells with caps, still under
# `SSWConfig.backslash_max_panels`, so the solver stays `Backslash` — which
# `GreenReconstruction` requires.
const KUTTAV1_DIAG2_MESHES = ((n_span=12, n_airfoil=21),
                              (n_span=24, n_airfoil=41),
                              (n_span=36, n_airfoil=61))
# The battery-I tags this stage must reproduce bitwise, per arm, at eta = 0.25.
const KUTTAV1_DIAG2_REPRO = Dict(
    "neumann" => "eta0p25_neumann_uncapped_tri",
    "dirichlet_vts" => "eta0p25_dirichlet_capped_tri")
# Battery-I meshladder tags, per arm, for the two rungs it also ran.
const KUTTAV1_DIAG2_MESH_REPRO = Dict(
    "neumann" => "neumann_uncapped_tri", "dirichlet_vts" => "dirichlet_capped_tri")
const KUTTAV1_DIAG2_REPRO_COLS = ("settled_CL", "CL_kj_unsteady",
    "steady_semiinfinite_CL")

# Rows arrive either as this run's NamedTuples or as strings parsed back out of a
# previously written CSV. `_csv_value` is the same %.16g formatter that wrote
# those files, so comparing formatted strings compares the stored numbers.
_repro_value(r::NamedTuple, col) = _csv_value(getproperty(r, Symbol(col)))
_repro_value(r::AbstractDict, col) = String(strip(r[col]))

_prior_rows(path) = isfile(path) ?
    Dict{String, Any}(r["cell"] => r for r in read_gate_csv(path)) :
    Dict{String, Any}()

"One bitwise-reproduction check row, or a failed row naming what was missing."
function _repro_check(newrow, oldrow, name, detail)
    (isnothing(newrow) || isnothing(oldrow)) &&
        return (; check=name, detail="$detail: row missing, not checked",
            value=NaN, pass=false)
    bad = [c for c in KUTTAV1_DIAG2_REPRO_COLS
        if _repro_value(newrow, c) != _repro_value(oldrow, c)]
    return (; check=name,
        detail=isempty(bad) ? detail : "$detail mismatched: " * join(bad, " "),
        value=Float64(length(bad)), pass=isempty(bad))
end

"""
Control assertion (plan §5.6): the Neumann and Dirichlet bodies place their
control points at the same points. `calc_controlpoints!` uses the triangle
centroid with zero offset for both, so with identical geometry the wake-induced
velocity at the control points is identical by construction — which is what
confines the formulation discrepancy to what each solve does with that velocity.
Compared on the uncapped mesh, the one both formulations can be built on.
"""
function _controlpoint_control(; n_span=KUTTAV1_GATE_COMMON.n_span,
        n_airfoil=KUTTAV1_GATE_COMMON.n_airfoil)
    base = (; KUTTAV1_GATE_COMMON..., n_span, n_airfoil, caps=:none)
    bn = build_suddenly_started_wing(SSWConfig(; base..., bodytype=:neumann))
    bd = build_suddenly_started_wing(SSWConfig(; base..., bodytype=:dirichlet))
    # Control points are populated by the solve, not the constructor; fill them
    # here or the comparison below is vacuously zero-vs-zero.
    pnl.calc_normals!(bn); pnl.calc_controlpoints!(bn)
    pnl.calc_normals!(bd); pnl.calc_controlpoints!(bd)
    same = size(bn.controlpoints) == size(bd.controlpoints) &&
        size(bn.nodes) == size(bd.nodes) && bn.cells == bd.cells
    cp_diff = same ? maximum(abs, bn.controlpoints .- bd.controlpoints) : NaN
    node_diff = same ? maximum(abs, bn.nodes .- bd.nodes) : NaN
    return [
        (; check="mesh_topology_identical",
            detail="ns$(n_span)_na$(n_airfoil)_capsnone",
            value=Float64(same), pass=same),
        (; check="max_abs_node_diff",
            detail="neumann_vs_dirichlet", value=node_diff,
            pass=node_diff == 0.0),
        (; check="max_abs_controlpoint_diff",
            detail="neumann_vs_dirichlet (scale $(maximum(abs, bn.controlpoints)))",
            value=cp_diff, pass=cp_diff == 0.0),
    ]
end

"""
    run_gate_diagnosis2(output; batteries)

T1 (formulation axis), T2 (freestream-convection geometry control), T3
(strip/row-1 strength jump, folded into every run) and T4 (mesh ladder under
each formulation). `:tri` basis only — the KJ circulation metric is basis-blind
(battery I: the `:tri`/`:tri_robust` spread collapses under refinement) and this
stage is about the solve, not the velocity reconstruction.

`batteries` selects sub-batteries (default all), overridable with
`KUTTAV1_DIAG2_BATTERIES=mesh`. The mesh ladder is an order of magnitude more
expensive than the other two and deserves to be launched on its own; output
files are suffixed when a subset is selected, so a partial run can never
overwrite a full one.
"""
function run_gate_diagnosis2(output;
        batteries=Tuple(Symbol(strip(x)) for x in
            split(get(ENV, "KUTTAV1_DIAG2_BATTERIES",
                "formulation,convection,mesh"), ',')))
    for b in batteries
        b in (:formulation, :convection, :mesh) ||
            error("unknown battery $b; use formulation, convection, mesh")
    end
    suffix = length(batteries) == 3 ? "" :
        "_" * join(sort(collect(string.(batteries))), "-")
    out(name) = joinpath(output, "gate_diagnosis2_$(name)$(suffix).csv")
    # Read every prior record BEFORE writing anything: a full run overwrites the
    # unsuffixed summary this function may want to reproduce against.
    prior_diag1 = _prior_rows(joinpath(output, "gate_diagnosis_summary.csv"))
    prior_diag2 = _prior_rows(joinpath(output, "gate_diagnosis2_summary.csv"))

    basis = _diag_basis(KUTTAV1_GATE_REFERENCE_BASIS)
    store = (; rows=NamedTuple[], gamma=NamedTuple[], history=NamedTuple[],
        strip=NamedTuple[])
    pairs = NamedTuple[]
    arm_of = Dict{String, String}()     # cell tag -> arm name

    # --- T1: formulation axis ---------------------------------------------
    for eta in (:formulation in batteries ? KUTTAV1_DIAG2_ETAS : ())
        tags = Dict{String, String}()
        for arm in KUTTAV1_DIAG2_ARMS
            tag = "eta$(_ssw_num_tag(eta))_$(arm.arm)"
            _diag_push!(store, arm, (; eta), basis; battery="formulation", tag,
                record_strip_jump=true)
            tags[arm.arm] = tag
            arm_of[tag] = arm.arm
        end
        for d in ("dirichlet_vts", "dirichlet_green")
            push!(pairs, (; battery="formulation", axis="eta", axis_value=eta,
                arm=d, neumann=tags["neumann"], dirichlet=tags[d]))
        end
    end

    # --- T2: geometry control ---------------------------------------------
    for eta in (:convection in batteries ? KUTTAV1_DIAG2_FC_ETAS : ())
        tags = Dict{String, String}()
        for arm in KUTTAV1_DIAG2_ARMS[1:2]
            tag = "fc_eta$(_ssw_num_tag(eta))_$(arm.arm)"
            _diag_push!(store, arm, (; eta, freestream_convection=true), basis;
                battery="freestream_convection", tag, record_strip_jump=true)
            tags[arm.arm] = tag
            arm_of[tag] = arm.arm
        end
        push!(pairs, (; battery="freestream_convection", axis="eta",
            axis_value=eta, arm="dirichlet_vts", neumann=tags["neumann"],
            dirichlet=tags["dirichlet_vts"]))
    end

    # --- T4: mesh ladder under each formulation ---------------------------
    # Tests whether the ~15% settled-CL floor that survives `GreenReconstruction`
    # is the force-reconstruction factor. eta is left at the SSWConfig default
    # (0.25, the gate value) so this is a spatial ladder only.
    for mesh in (:mesh in batteries ? KUTTAV1_DIAG2_MESHES : ())
        tags = Dict{String, String}()
        for arm in KUTTAV1_DIAG2_ARMS
            tag = "ns$(mesh.n_span)na$(mesh.n_airfoil)_$(arm.arm)"
            _diag_push!(store, arm, mesh, basis; battery="meshladder", tag,
                record_strip_jump=true)
            tags[arm.arm] = tag
            arm_of[tag] = arm.arm
        end
        for d in ("dirichlet_vts", "dirichlet_green")
            push!(pairs, (; battery="meshladder", axis="ncells_proxy",
                axis_value=Float64(2 * (mesh.n_airfoil - 1) * mesh.n_span),
                arm=d, neumann=tags["neumann"], dirichlet=tags[d]))
        end
    end
    isempty(store.rows) && error("no batteries selected; nothing to run")

    write_namedtuple_csv(out("summary"), store.rows)
    write_namedtuple_csv(out("cl_history"), store.history)
    write_namedtuple_csv(out("spanwise_gamma"), store.gamma)
    write_namedtuple_csv(out("strip_jump"), store.strip)

    byname = Dict(r.cell => r for r in store.rows)
    gsel(tag) = [r for r in store.gamma if r.case == tag]
    gaps = NamedTuple[]
    for p in pairs
        n, d = byname[p.neumann], byname[p.dirichlet]
        agree = spanwise_agreement(gsel(p.neumann), gsel(p.dirichlet))
        push!(gaps, (; p.battery, p.axis, p.axis_value, p.arm,
            neumann_CL=n.settled_CL, dirichlet_CL=d.settled_CL,
            formulation_gap(n, d)..., gap_factorization(n, d)...,
            agree.shape_correlation, agree.scale_ratio, agree.sign_flips,
            ncells_neumann=n.ncells, ncells_dirichlet=d.ncells,
            status=n.status == "ok" && d.status == "ok" ? "ok" :
                "$(n.status)|$(d.status)"))
    end
    write_namedtuple_csv(out("gaps"), gaps)

    # Primary metric: circulation recovery against each body's OWN semi-infinite
    # anchor, which is what makes the Dirichlet/Neumann comparison a statement
    # about small-Das sensitivity rather than about the two bodies' levels.
    recovery = NamedTuple[]
    for r in store.rows
        arm = arm_of[r.cell]
        push!(recovery, (; r.cell, r.battery, arm, r.eta,
            r.freestream_convection, r.formulation, r.bodytype, r.caps,
            r.n_span, r.n_airfoil, r.ncells,
            r.CL_kj_unsteady, r.CL_kj_steady, r.settled_CL,
            r.steady_semiinfinite_CL,
            kj_recovery=abs(r.CL_kj_unsteady) / abs(r.CL_kj_steady),
            pressure_recovery=abs(r.settled_CL) / abs(r.steady_semiinfinite_CL),
            # Each body's own pressure CL per unit of its own KJ CL. The
            # Dirichlet/Neumann ratio of THIS column is the reconstruction factor
            # in `gap_factorization`; published per cell so the factorization can
            # be re-derived from the recovery table alone.
            pressure_per_kj=r.settled_CL / abs(r.CL_kj_unsteady),
            steady_pressure_per_kj=r.steady_semiinfinite_CL / abs(r.CL_kj_steady),
            r.das_parallel_err, r.das_drift, r.status))
    end
    write_namedtuple_csv(out("recovery"), recovery)

    # --- harness self-checks ----------------------------------------------
    checks = _controlpoint_control()
    lookup(d, k) = get(d, k, nothing)

    # Bitwise reproduction of the battery-I cells this stage re-runs unchanged.
    for (arm, oldtag) in sort(collect(KUTTAV1_DIAG2_REPRO))
        haskey(byname, "eta0p25_$(arm)") || continue
        push!(checks, _repro_check(lookup(byname, "eta0p25_$(arm)"),
            lookup(prior_diag1, oldtag),
            "repro_$(arm)_vs_gate_diagnosis", oldtag))
    end
    # Mesh-ladder rungs against battery I's own mesh ladder, where they overlap.
    for mesh in KUTTAV1_DIAG2_MESHES, (arm, oldsuffix) in
            sort(collect(KUTTAV1_DIAG2_MESH_REPRO))
        stem = "ns$(mesh.n_span)na$(mesh.n_airfoil)"
        newtag, oldtag = "$(stem)_$(arm)", "$(stem)_$(oldsuffix)"
        # Rung 3 has no battery-I counterpart; that is expected, not a failure.
        (haskey(byname, newtag) && haskey(prior_diag1, oldtag)) || continue
        push!(checks, _repro_check(byname[newtag], prior_diag1[oldtag],
            "repro_mesh_$(stem)_$(arm)_vs_gate_diagnosis", oldtag))
    end
    # Mesh rung 1 is the same configuration as the T1 eta=0.25 cells; reproduce
    # them whether they were computed in this run or a previous one.
    for arm in KUTTAV1_DIAG2_ARMS
        stem = "ns$(KUTTAV1_DIAG2_MESHES[1].n_span)na$(KUTTAV1_DIAG2_MESHES[1].n_airfoil)"
        newtag, reftag = "$(stem)_$(arm.arm)", "eta0p25_$(arm.arm)"
        haskey(byname, newtag) || continue
        ref = haskey(byname, reftag) ? byname[reftag] : lookup(prior_diag2, reftag)
        push!(checks, _repro_check(byname[newtag], ref,
            "repro_mesh_rung1_$(arm.arm)_vs_formulation_battery", reftag))
    end

    # The semi-infinite anchor is a property of the body AND ITS MESH; it cannot
    # depend on eta, formulation, or convection mode. Grouping by (bodytype,
    # caps) alone would fail spuriously the moment the mesh ladder runs, because
    # `CL_kj_steady` legitimately refines (-0.399434 at 480 vs -0.396583 at 1920).
    for key in sort(unique((string(r.bodytype), string(r.caps), r.n_span,
            r.n_airfoil) for r in store.rows))
        vals = [r.CL_kj_steady for r in store.rows
            if (string(r.bodytype), string(r.caps), r.n_span, r.n_airfoil) == key &&
                isfinite(r.CL_kj_steady)]
        spread = isempty(vals) ? NaN : maximum(vals) - minimum(vals)
        push!(checks, (; check="steady_anchor_invariant",
            detail="$(key[1])_caps$(key[2])_ns$(key[3])na$(key[4]) over " *
                "$(length(vals)) cells, " *
                "CL_kj_steady=$(isempty(vals) ? NaN : first(vals))",
            value=spread, pass=!isempty(vals) && spread == 0.0))
    end
    # The circulation x reconstruction split is an identity, so its residual is a
    # free check on the whole decomposition this battery reports.
    let residuals = [abs(g.factorization_residual) for g in gaps
            if isfinite(g.factorization_residual)]
        worst = isempty(residuals) ? NaN : maximum(residuals)
        push!(checks, (; check="gap_factorization_identity",
            detail="max |circ*recon - settled ratio| over $(length(residuals)) pairs",
            value=worst, pass=!isempty(residuals) && worst <= 1e-12))
    end
    # Das control (plan §5.6), aggregated over every cell.
    for (name, col) in (("das_parallel_to_freestream", :das_parallel_err),
            ("das_frozen_over_run", :das_drift))
        worst = maximum(getproperty(r, col) for r in store.rows)
        push!(checks, (; check=name, detail="max over $(length(store.rows)) cells",
            value=worst, pass=isfinite(worst) && worst <= 1e-12))
    end
    write_namedtuple_csv(out("checks"), checks)

    println()
    @printf("%-28s %-22s %8s %7s %10s %10s\n", "cell", "formulation", "eta",
        "ncells", "kj_recov", "p_recov")
    for r in recovery
        @printf("%-28s %-22s %8.4f %7d %10.4f %10.4f\n", r.cell, r.formulation,
            r.eta, r.ncells, r.kj_recovery, r.pressure_recovery)
    end
    println()
    # `recon` next to `steady_recon` IS the test of the force-reconstruction
    # claim; printed adjacent so the comparison cannot be missed.
    @printf("%-22s %-13s %10s %-16s %9s %9s %9s %9s\n", "battery", "axis",
        "value", "arm", "settled%", "KJ%", "recon", "steady")
    for g in gaps
        @printf("%-22s %-13s %10.4g %-16s %9.2f %9.2f %9.4f %9.4f  %s\n",
            g.battery, g.axis, g.axis_value, g.arm,
            100 * g.unsteady_settled_gap, 100 * g.unsteady_kj_gap,
            g.reconstruction_ratio, g.steady_reconstruction_ratio, g.status)
    end
    println()
    for c in checks
        @printf("check %-40s %-12.5g pass=%-5s  %s\n", c.check, c.value,
            string(c.pass), c.detail)
    end
    if !isempty(store.strip)
        println()
        @printf("%-28s %10s %14s %14s\n", "cell", "gamma_max",
            "mean|jump|/g", "max|jump|/g")
        for key in sort(unique(r.cell for r in store.strip))
            rs = [r for r in store.strip if r.cell == key]
            @printf("%-28s %10.5f %14.5f %14.5f\n", key, first(rs).gamma_scale,
                mean(r.mean_jump_over_gamma for r in rs),
                maximum(r.max_jump_over_gamma for r in rs))
        end
    end
    return (; rows=store.rows, gaps, recovery, checks, strip=store.strip)
end

# The "pre-closure" residual is NOT measured inside the pressure run; it is the
# post-closure residual of the *jump* run on the same route, planted here as the
# closure's starting point. The column name says so.
function _fill_pre_dp(rows)
    jump_dp = Dict((r.route, r.max_dp_hat_post) for r in rows if r.closure == "jump")
    return [merge(r, (; max_dp_hat_pre_from_jump_run=get(jump_dp, r.route, NaN)))
        for r in rows]
end

function run_diamond(output)
    rows = NamedTuple[]
    gamma = NamedTuple[]
    for combo in KUTTAV1_COMBOS
        println("diamond: $(combo.label)")
        result = run_diamond_combo(combo)
        push!(rows, result.summary)
        append!(gamma, result.gamma)
    end
    rows = _fill_pre_dp(rows)
    write_namedtuple_csv(joinpath(output, "diamond_summary.csv"), rows)
    write_namedtuple_csv(joinpath(output, "diamond_spanwise_gamma.csv"), gamma)
    return rows
end

"""
Read a gate CSV and return its rows as `Dict{String,String}`. Gate records are
written by `write_namedtuple_csv` and contain no comma-bearing fields, so a
plain split is sufficient — but columns are addressed by *name*, never by
position, so adding a column can no longer silently redefine the pass flag.
"""
function read_gate_csv(path)
    isfile(path) || error("missing gate record $path")
    lines = filter(!isempty, strip.(readlines(path)))
    length(lines) >= 2 || error("gate record $path has no data rows")
    header = string.(split(lines[1], ','))
    return [Dict(zip(header, string.(split(line, ',')))) for line in lines[2:end]]
end

function assert_gate_passed(path)
    for row in read_gate_csv(path)
        haskey(row, "pass") || error("gate record $path has no `pass` column")
        lowercase(strip(row["pass"])) == "true" ||
            error("prerequisite gate did not pass: $path ($(row))")
    end
    return true
end

function run_ssw0(output)
    gate_path = joinpath(output, "ssw_sanity_gate.csv")
    isfile(gate_path) || run_gate(output)
    assert_gate_passed(gate_path)

    config = SSWConfig(AR=6.0, alpha_deg=0.0, n_span=12, n_airfoil=21,
        t_end_star=2.0, dt_star=0.125, backend_kind=:direct,
        bodytype=:dirichlet, save_vtk=false, verbose=false)
    rows = NamedTuple[]
    gamma = NamedTuple[]
    history = NamedTuple[]
    for combo in KUTTAV1_COMBOS
        println("ssw0: $(combo.label)")
        result = run_ssw_combo(config, combo)
        push!(rows, result.summary)
        append!(gamma, result.gamma)
        append!(history, result.history)
    end
    rows = _fill_pre_dp(rows)
    write_namedtuple_csv(joinpath(output, "ssw0_summary.csv"), rows)
    write_namedtuple_csv(joinpath(output, "ssw0_spanwise_gamma.csv"), gamma)
    write_namedtuple_csv(joinpath(output, "ssw0_cl_history.csv"), history)
    max_cl = maximum(abs(r.CL) for r in rows)
    max_dp = maximum(r.max_dp_hat_post for r in rows)
    max_c = maximum(r.c_inf_scaled for r in rows)
    max_sym = maximum(r.gamma_symmetry for r in rows)
    gate = (; gate="ssw0_symmetry", max_abs_CL=max_cl, max_dp_hat=max_dp,
        max_c_scaled=max_c, max_gamma_symmetry=max_sym,
        pass=max_cl <= 1e-10 && max_dp <= 1e-8 && max_c <= 1e-8 &&
            max_sym <= 1e-10)
    write_namedtuple_csv(joinpath(output, "ssw0_gate.csv"), [gate])
    gate.pass || error("SSW zero-lift gate failed: $gate")
    return rows
end

const KUTTAV1_STAGE_GATES = Dict(
    :gate => ("ssw_sanity_gate.csv",),
    :ssw0 => ("ssw_sanity_gate.csv", "ssw0_gate.csv"),
    :ssw_alpha => ("ssw_sanity_gate.csv", "ssw0_gate.csv"),
    :sweptwing => ("ssw_sanity_gate.csv", "ssw0_gate.csv"),
)

function require_ssw_gates(output, stage=:ssw_alpha)
    for filename in get(KUTTAV1_STAGE_GATES, stage, ())
        path = joinpath(output, filename)
        isfile(path) ||
            error("missing prerequisite $path; run KUTTAV1_STAGE=ssw0 first")
        assert_gate_passed(path)
    end
    return true
end

function run_ssw_alpha(output)
    require_ssw_gates(output)
    config = SSWConfig(AR=6.0, alpha_deg=5.0, n_span=12, n_airfoil=21,
        t_end_star=8.0, dt_star=0.125, backend_kind=:direct,
        bodytype=:dirichlet, save_vtk=false, verbose=false)
    rows = NamedTuple[]
    gamma = NamedTuple[]
    history = NamedTuple[]
    for combo in KUTTAV1_COMBOS
        println("ssw_alpha: $(combo.label)")
        result = run_ssw_combo(config, combo)
        push!(rows, result.summary)
        append!(gamma, result.gamma)
        append!(history, result.history)
    end
    rows = _fill_pre_dp(rows)
    write_namedtuple_csv(joinpath(output, "ssw_alpha_summary.csv"), rows)
    write_namedtuple_csv(joinpath(output, "ssw_alpha_spanwise_gamma.csv"), gamma)
    write_namedtuple_csv(joinpath(output, "ssw_alpha_cl_history.csv"), history)
    base = only(r for r in rows if r.combo == "A_jump")
    ap = only(r for r in rows if r.combo == "A_pressure")
    bj = only(r for r in rows if r.combo == "B_jump")
    bp = only(r for r in rows if r.combo == "B_pressure")
    attribution = [
        (; metric="CL", baseline=base.CL, closure=ap.CL - base.CL,
            geometry=bj.CL - base.CL,
            composition=bp.CL - base.CL,
            interaction=bp.CL - ap.CL - bj.CL + base.CL),
        (; metric="max_dp_hat_post", baseline=base.max_dp_hat_post,
            closure=ap.max_dp_hat_post - base.max_dp_hat_post,
            geometry=bj.max_dp_hat_post - base.max_dp_hat_post,
            composition=bp.max_dp_hat_post - base.max_dp_hat_post,
            interaction=bp.max_dp_hat_post - ap.max_dp_hat_post -
                bj.max_dp_hat_post + base.max_dp_hat_post),
    ]
    write_namedtuple_csv(joinpath(output, "ssw_alpha_attribution.csv"), attribution)
    return rows
end

function make_sweptwing_body()
    if !isdefined(@__MODULE__, :simplewing_mirrored)
        include(joinpath(@__DIR__, "helper_functions.jl"))
    end
    b = 98 * 0.0254
    bodytype = pnl.RigidWakeBody{Union{pnl.ConstantSource, pnl.ConstantDoublet}}
    # Interior Dirichlet assumes a closed surface; `simplewing` lofts open tips.
    body = simplewing_mirrored(b, 5.0, 1.0, 0, 0, 45, 0;
        bodytype, bodyoptargs=(; semiinfinite_wake=false), caps=:flat,
        airfoil_root="airfoil-rae101.csv", airfoil_tip="airfoil-rae101.csv",
        airfoil_path=joinpath(@__DIR__, "data"),
        rfl_NDIVS=[(0.25, 24, 10.0, false), (0.50, 24, 1.0, true),
            (0.25, 24, 0.1, false)],
        span_NDIVS=[(1.0, 24, 1.0, true)])
    Vinf = 30.0 .* [cosd(4.2), 0.0, sind(4.2)]
    for Das in body.Das
        Das .= repeat(Vinf ./ LA.norm(Vinf), 1, size(Das, 2))
    end
    return body, Vinf, b
end

function run_sweptwing_combo(combo)
    combo.route == :A || error("steady swept-wing V1 supports Route A only")
    rho = 1.225
    body, Vinf, b = make_sweptwing_body()
    Sref = b^2 / 5.0
    c_ref = b / 5.0
    Dhat = Vinf / LA.norm(Vinf)
    Shat = [0.0, 1.0, 0.0]
    Lhat = LA.cross(Dhat, Shat)
    pressure = pnl.PressureBernoulli(rho)
    force = pnl.ForceMonitor(1, 1; i_frame=-1,
        normalization=pnl.WingNormalization(rho, Sref, c_ref),
        correct_kuttacondition=false, verbose=false)
    spanwise = pnl.SpanwiseLoadingMonitor(24, 1;
        components=(lift=Lhat, drag=Dhat), span_axis=Shat, per_length=true,
        normalization=pnl.NoSectionalNormalization())
    attachment, closure = combo_objects(combo, rho)
    elapsed = @elapsed pnl.steady!(body, pnl.ReferenceFrame(body), Vinf;
        body_solvers=pnl.Backslash(body), backend=pnl.DirectBackend(),
        monitors=(pressure, force, spanwise), path=nothing,
        grad_mu_options=(; basis=:quad), wake_attachment=attachment,
        kutta_closure=closure, verbose=false)
    dm = diagnostic_metrics(closure)
    q = 0.5 * rho * LA.dot(Vinf, Vinf)
    te = paired_te_metrics(body, rho, q)
    gamma = te_gamma_rows("sweptwing", combo.label, body)
    lift_i = only(findall(==(:lift), spanwise.component_names))
    sectional = [(; case="sweptwing", combo=combo.label, bin=i,
        y2b=2 * spanwise.bin_center[i] / b,
        sectional_cl=spanwise.load_components[lift_i, i] / (q * c_ref))
        for i in eachindex(spanwise.bin_center)]
    summary = (; case="sweptwing", combo=combo.label, route=:A,
        closure=combo.pressure ? "pressure" : "jump",
        CL=LA.dot(force.force[:, 1], Lhat),
        dp_hat_scale=q, kutta_pressure_scale=dm.pressure_scale,
        max_dp_hat_pre_from_jump_run=NaN,
        max_dp_hat_post=te.max_dp_hat, median_dp_hat_post=te.median_dp_hat,
        dm.c_inf_scaled, dm.mean_outer_iterations, dm.body_solves,
        dm.mean_body_solves, dm.startup_steps, dm.startup_body_solves,
        dm.backtracks, dm.restarts, dm.max_r_inf_scaled,
        dm.status, dm.disposition, elapsed_s=elapsed,
        gamma_symmetry=gamma_symmetry(gamma))
    return (; summary, gamma, sectional)
end

function run_sweptwing(output)
    require_ssw_gates(output)
    rows = NamedTuple[]
    gamma = NamedTuple[]
    sectional = NamedTuple[]
    for combo in KUTTAV1_COMBOS[1:2]
        println("sweptwing: $(combo.label)")
        result = run_sweptwing_combo(combo)
        push!(rows, result.summary)
        append!(gamma, result.gamma)
        append!(sectional, result.sectional)
    end
    rows = _fill_pre_dp(rows)
    write_namedtuple_csv(joinpath(output, "sweptwing_summary.csv"), rows)
    write_namedtuple_csv(joinpath(output, "sweptwing_spanwise_gamma.csv"), gamma)
    write_namedtuple_csv(joinpath(output, "sweptwing_sectional_loading.csv"), sectional)
    write_namedtuple_csv(joinpath(output, "sweptwing_route_b_na.csv"),
        [(; combo="B_jump", status="N/A",
            reason="steady! rejects TEAnchoredAttachment by design"),
         (; combo="B_pressure", status="N/A",
            reason="steady! rejects TEAnchoredAttachment by design")])
    return rows
end

function main_kuttav1()
    Threads.nthreads() <= 6 ||
        error("Kutta V1 local runs are limited to at most 6 Julia threads")
    stage = Symbol(lowercase(get(ENV, "KUTTAV1_STAGE", "diamond")))
    output = get(ENV, "KUTTAV1_OUTPUT", KUTTAV1_DEFAULT_OUTPUT)
    mkpath(output)
    summary = joinpath(output, "$(stage)_summary.csv")
    if isfile(summary) && !_kuttav1_bool("KUTTAV1_FORCE", false)
        # A checkpoint records that the stage RAN, not that its gates passed.
        # Re-verify, so a resume can never present a failed hard gate as a
        # completed stage.
        require_ssw_gates(output, stage)
        println("Resume checkpoint exists: $summary")
        return summary
    end
    if stage == :diamond
        run_diamond(output)
    elseif stage == :gate
        run_gate(output)
    elseif stage == :gate_diagnosis
        run_gate_diagnosis(output)
    elseif stage == :gate_diagnosis2
        run_gate_diagnosis2(output)
    elseif stage == :ssw0
        run_ssw0(output)
    elseif stage == :ssw_alpha
        run_ssw_alpha(output)
    elseif stage == :sweptwing
        run_sweptwing(output)
    else
        error("KUTTAV1_STAGE must be diamond, gate, gate_diagnosis, " *
            "gate_diagnosis2, ssw0, ssw_alpha, or sweptwing")
    end
    println("Completed Kutta V1 stage $stage; output: $output")
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_kuttav1()
end
