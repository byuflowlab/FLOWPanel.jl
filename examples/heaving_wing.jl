#=##############################################################################
# DESCRIPTION
    Heaving and pitching rectangular NACA 0012 wing. This example mirrors the
    setup in VortexLattice.jl/test/heaving_wing.jl, using FLOWPanel's finite-
    thickness lifting-body and unsteady panel-wake machinery.

# AUTHORSHIP
  * Created   : Jun 2026
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
using FLOWPanel.FastMultipole.StaticArrays
import Printf: @printf
import Statistics: mean

const DEFAULT_HEAVING_WING_PATH = joinpath("data", "heaving_wing")

function naca0012_contour(n::Integer=121; thickness::Real=0.12)
    n >= 21 || throw(ArgumentError("naca0012_contour requires at least 21 points"))

    n_half = max(11, cld(n, 2))
    β = range(0.0, pi, length=n_half)
    x = 0.5 .* (1 .- cos.(β))
    yt = 5 * thickness .* (
        0.2969 .* sqrt.(x) .- 0.1260 .* x .- 0.3516 .* x.^2 .+
        0.2843 .* x.^3 .- 0.1036 .* x.^4
    )

    # Start at the trailing edge, proceed around the lower side to the leading
    # edge, then return along the upper side. The trailing edge is a single node.
    lower = hcat(reverse(x), -reverse(yt))
    upper = hcat(x[2:end-1], yt[2:end-1])
    return vcat(lower, upper)
end

function heaving_wing_mesh(c, b; n_span::Integer=13, n_airfoil::Integer=161)
    contour = naca0012_contour(n_airfoil)
    n_sec = size(contour, 1)
    y = collect(range(-b/2, stop=b/2, length=n_span + 1))

    node_index(i, j) = i + (j - 1) * n_sec
    nodes = zeros(Float64, 3, n_sec * length(y))
    for (j, yj) in enumerate(y), i in 1:n_sec
        nodes[:, node_index(i, j)] .= (c * contour[i, 1], yj, c * contour[i, 2])
    end

    cells = zeros(Int, 3, 2 * n_sec * n_span)
    k = 0
    for j in 1:n_span
        for i in 1:n_sec
            ip = i == n_sec ? 1 : i + 1
            n11 = node_index(i, j)
            n21 = node_index(ip, j)
            n12 = node_index(i, j + 1)
            n22 = node_index(ip, j + 1)
            k += 1
            cells[:, k] .= (n11, n21, n22)
            k += 1
            cells[:, k] .= (n11, n22, n12)
        end
    end

    shedding = calc_heaving_wing_shedding(nodes, cells, c)

    return nodes, cells, [shedding]
end

function calc_heaving_wing_shedding(nodes, cells, c)
    x_te = float(c)
    te_tol = max(100 * eps(Float64) * max(abs(x_te), 1.0), 1e-8 * abs(x_te))
    te_nodes = findall(i -> isapprox(nodes[1, i], x_te; atol=te_tol, rtol=0), axes(nodes, 2))
    length(te_nodes) >= 2 || throw(ArgumentError(
        "heaving wing mesh requires at least two trailing-edge nodes at x=$(x_te); " *
        "found $(length(te_nodes))"))

    sort!(te_nodes, by = i -> nodes[2, i])
    lower = [
        x_te - te_tol,
        minimum(nodes[2, te_nodes]) - te_tol,
        minimum(nodes[3, te_nodes]) - te_tol,
    ]
    upper = [
        x_te + te_tol,
        maximum(nodes[2, te_nodes]) + te_tol,
        maximum(nodes[3, te_nodes]) + te_tol,
    ]
    bbox = (lower, upper)

    return pnl.calc_shedding_from_seed(nodes, cells, te_nodes[1], te_nodes[2];
        end_node=te_nodes[end],
        bbox,
        normal_jump_tol=1.0,
        max_turn_angle=pi/2,
    )
end

function set_chordline_Das!(body, c)
    xhat = SVector{3}(1.0, 0.0, 0.0)
    min_positive = eps(float(c))

    for (shedding, Das) in zip(body.shedding, body.Das)
        for (j, col) in enumerate(eachcol(shedding))
            te_panel, _, nib, _, _, _ = col
            te_node = body.cells[nib, te_panel]

            best_length = float(c)
            for node in body.cells[:, te_panel]
                node == te_node && continue
                edge = SVector{3}(
                    body.nodes[1, te_node] - body.nodes[1, node],
                    body.nodes[2, te_node] - body.nodes[2, node],
                    body.nodes[3, te_node] - body.nodes[3, node],
                )
                chordwise_length = abs(edge[1])
                chordwise_length > min_positive && (best_length = min(best_length, chordwise_length))
            end

            Das[:, j] .= xhat .* best_length
        end

        te_panel, nia, _, _, _, _ = shedding[:, end]
        te_node = body.cells[nia, te_panel]
        best_length = float(c)
        for node in body.cells[:, te_panel]
            node == te_node && continue
            edge = SVector{3}(
                body.nodes[1, te_node] - body.nodes[1, node],
                body.nodes[2, te_node] - body.nodes[2, node],
                body.nodes[3, te_node] - body.nodes[3, node],
            )
            chordwise_length = abs(edge[1])
            chordwise_length > min_positive && (best_length = min(best_length, chordwise_length))
        end
        Das[:, end] .= xhat .* best_length
    end

    return body
end

function heaving_wing_frame(body, pivot, θ_init, h0, θ0, ω, ψ)
    Rinit = pnl.Rodrigues(SVector{3}(0.0, 1.0, 0.0), θ_init)
    pnl.rotate_translate!(body, pivot, Rinit, SVector{3}(0.0, 0.0, 0.0))
    pnl.rotate_Das!(body, Rinit)

    return pnl.ReferenceFrame(body;
        origin = pivot,
        v = SVector{3}(0.0, 0.0, h0 * ω),
        ω_axis = SVector{3}(0.0, 1.0, 0.0),
        ω = θ0 * ω * cos(ψ),
        R = Rinit,
        name = "wing",
        child_index = Int[],
        dependent_index = [1],
    )
end

function prepare_heaving_wing(;
        St::Real=0.3,
        heave_amplitude::Real=0.25,
        pitch_amplitude::Real=15*pi/180,
        phase::Real=90*pi/180,
        n_cycles::Real=5,
        c_per_dt::Real=0.02,
        b_pivot::Real=1/3,
        c::Real=3.81e-2,
        Re_c::Real=1100.0,
        nu::Real=1e-6,
        rho::Real=1.225,
        AR::Real=6,
        n_span::Integer=13,
        n_airfoil::Integer=31,
        panel_wake_rows::Integer=1,
        max_particles::Integer=20000,
        save_vtk::Bool=true,
        path=DEFAULT_HEAVING_WING_PATH,
        include_kutta_joukowski::Bool=false,
        backend=pnl.FastMultipoleBackend(expansion_order=8, multipole_acceptance=0.4, leaf_size=40),
    )
    U = Re_c * nu / c
    h0 = heave_amplitude * c
    A = 2 * h0
    θ0 = pitch_amplitude
    ψ = phase
    ω = St * U / A * 2*pi
    period = 2*pi / ω
    dt = c / U * c_per_dt
    t_range = range(0.0, stop=n_cycles * period, step=dt)
    length(t_range) >= 2 || throw(ArgumentError(
        "heaving wing simulation requires at least two time samples; got $(length(t_range)). " *
        "Increase n_cycles or decrease c_per_dt."))

    b = AR * c
    S0 = b * c
    pivot = SVector{3}(b_pivot * c, 0.0, 0.0)

    bodytype = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64, false}
    bodyoptargs = (;
        cp_outer=true,
        kerneloffset=1e-6 * c,
        kernelcutoff=1e-12 * c,
        semiinfinite_wake=false,
        watertight=false,
    )

    nodes, cells, shedding = heaving_wing_mesh(c, b; n_span, n_airfoil)
    wing = bodytype(nodes, cells, shedding; bodyoptargs...)
    set_chordline_Das!(wing, c)

    θ_init = θ0 * sin(ψ)
    frames = heaving_wing_frame(wing, pivot, θ_init, h0, θ0, ω, ψ)
    Uinf(t) = SVector{3}(U, 0.0, 0.0)

    function heaving_maneuver!(frames, systems, wakes, t)
        frame = frames[1]
        this_v = h0 * ω * cos(ω * t)
        this_ω = θ0 * ω * cos(ω * t + ψ)
        frames[1] = typeof(frame)(
            frame.x,
            SVector{3}(0.0, 0.0, this_v),
            frame.ω_axis,
            this_ω,
            frame.R,
            frame.Rp2g,
            frame.name,
            frame.parent_index,
            frame.child_index,
            frame.dependent_index,
        )
        return nothing
    end

    # wake = pnl.PanelWake(wing; nwakerows=panel_wake_rows)
    wake = pnl.PanelParticleWake(wing;
        nwakerows=panel_wake_rows,
        max_particles=max_particles,
        method_trailing=pnl.OverlapPPS(1.3, 2),
        method_unsteady=pnl.OverlapPPS(1.3, 2),
    )
    solver = pnl.Backslash(wing)
    normalization = pnl.WingNormalization(rho, S0, c)

    pressure_monitor = pnl.PressureBernoulli(rho; unsteady=true, allow_partial=true)
    force_monitor = pnl.ForceMonitor(length(t_range), 1;
        normalization,
        i_frame=-1,
        verbose=false,
    )

    monitors = if include_kutta_joukowski
        kj_monitor = pnl.KuttaJoukowskiForce(wing, length(t_range), 1;
            rho,
            backend,
            normalization,
            i_frame=-1,
            verbose=false,
        )
        (pressure_monitor, force_monitor, kj_monitor)
    else
        (pressure_monitor, force_monitor)
    end

    return (;
        wing,
        wake,
        frames,
        maneuver! = heaving_maneuver!,
        Uinf,
        t_range,
        solver,
        backend,
        monitors,
        pressure_monitor,
        force_monitor,
        setup = (;
            St,
            heave_amplitude,
            pitch_amplitude,
            phase,
            n_cycles,
            c_per_dt,
            b_pivot,
            c,
            Re_c,
            nu,
            rho,
            AR,
            b,
            S0,
            U,
            h0,
            A,
            ω,
            period,
            dt,
            n_span,
            n_airfoil,
            panel_wake_rows,
            max_particles,
            θ_init,
            path = save_vtk ? path : nothing,
        ),
    )
end

function _write_heaving_wing_forces_csv(path, t_range, period, CT_history, CL_history)
    mkpath(path)
    filename = joinpath(path, "heaving_wing_forces.csv")
    open(filename, "w") do io
        println(io, "time,t_over_T,CT,CL")
        for i in eachindex(t_range)
            @printf(io, "%.16e,%.16e,%.16e,%.16e\n",
                t_range[i], t_range[i] / period, CT_history[i], CL_history[i])
        end
    end
    return filename
end

function _finite_summary(values)
    finite = filter(isfinite, collect(values))
    isempty(finite) && return "all NaN/Inf"
    return "min=$(round(minimum(finite), sigdigits=5)), max=$(round(maximum(finite), sigdigits=5)), last=$(round(finite[end], sigdigits=5))"
end

function run_heaving_wing(; kwargs...)
    sim = prepare_heaving_wing(; kwargs...)
    (; wing, wake, frames, maneuver!, Uinf, t_range, solver, backend, monitors, force_monitor, setup) = sim

    println("Heaving/pitching NACA 0012 wing")
    @printf("  c = %.6g m, b = %.6g m, AR = %.4g, S = %.6g m^2\n", setup.c, setup.b, setup.AR, setup.S0)
    @printf("  Re_c = %.6g, nu = %.6g m^2/s, U = %.6g m/s, rho = %.6g kg/m^3\n",
        setup.Re_c, setup.nu, setup.U, setup.rho)
    @printf("  St = %.6g, h0/c = %.6g, theta0 = %.6g deg, phase = %.6g deg\n",
        setup.St, setup.heave_amplitude, setup.pitch_amplitude * 180/pi, setup.phase * 180/pi)
    @printf("  period = %.6g s, dt = %.6g s, steps = %d, steps/cycle ~= %.1f\n",
        setup.period, setup.dt, length(t_range), setup.period / setup.dt)
    println("  panels = $(wing.ncells), nodes = $(wing.nnodes), shedding edges = $(wing.nsheddings)")

    pnl.simulate!((wing,), (wake,), frames, maneuver!, Uinf, t_range;
        body_solvers=(solver,),
        backend,
        monitors,
        path=setup.path,
        name="heaving_wing",
        set_Das_eta_freestream=NaN,
        verbose=true,
    )

    CT_history = force_monitor.force[1, :]
    CL_history = force_monitor.force[3, :]
    post_transient = findall(t -> t >= setup.period, collect(t_range))

    if isempty(post_transient)
        CT_mean = NaN
        CL_mean = NaN
        CL_pp = NaN
    else
        CT_mean = mean(CT_history[post_transient])
        CL_mean = mean(CL_history[post_transient])
        CL_pp = maximum(CL_history[post_transient]) - minimum(CL_history[post_transient])
    end

    csv_dir = isnothing(setup.path) ? DEFAULT_HEAVING_WING_PATH : setup.path
    csv_path = _write_heaving_wing_forces_csv(csv_dir,
        t_range, setup.period, CT_history, CL_history)

    println("\nForce coefficient summary")
    @printf("  CT mean after first cycle = %.8g\n", CT_mean)
    @printf("  CL mean after first cycle = %.8g\n", CL_mean)
    @printf("  CL peak-to-peak after first cycle = %.8g\n", CL_pp)
    println("  CT history: $(_finite_summary(CT_history))")
    println("  CL history: $(_finite_summary(CL_history))")
    println("  Wrote force history CSV: $(csv_path)")
    if !isnothing(setup.path)
        println("  Wrote VTK output under: $(setup.path)")
    end

    return merge(sim, (; CT_history, CL_history, CT_mean, CL_mean, CL_pp, csv_path))
end

# if abspath(PROGRAM_FILE) == @__FILE__
#     run_heaving_wing()
# end
run_heaving_wing()
