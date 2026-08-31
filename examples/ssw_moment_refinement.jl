# Spanwise sheet-moment refinement diagnostic for BRAINSTORM 017.

if !isdefined(@__MODULE__, :SSWConfig)
    include(joinpath(@__DIR__, "suddenly_started_wing.jl"))
end

import Serialization

function ssw_sheet_moments(wake::pnl.PanelWake)
    circulation = zeros(3)
    impulse = zeros(3)
    for surface in eachindex(wake.nodes), row in 1:wake.nwakes[]
        nodes = wake.nodes[surface]
        strength = wake.strength[surface]
        for column in axes(strength, 3)
            gamma = strength[1, row, column]
            ring = (
                (SVector{3}(nodes[:, row, column]),
                    SVector{3}(nodes[:, row + 1, column])),
                (SVector{3}(nodes[:, row + 1, column]),
                    SVector{3}(nodes[:, row + 1, column + 1])),
                (SVector{3}(nodes[:, row + 1, column + 1]),
                    SVector{3}(nodes[:, row, column + 1])),
                (SVector{3}(nodes[:, row, column + 1]),
                    SVector{3}(nodes[:, row, column])),
            )
            for (r1, r2) in ring
                vector_gamma = gamma * (r2 - r1)
                circulation .+= vector_gamma
                impulse .+= LA.cross(0.5 * (r1 + r2), vector_gamma)
            end
        end
    end
    return (; circulation, impulse)
end

function _sswmr_run_state(n_span, output)
    state_path = joinpath(output, "rolled_nspan$(n_span).jls")
    isfile(state_path) && return Serialization.deserialize(state_path)
    config = SSWConfig(AR=6.0, n_span=n_span, n_airfoil=21,
        dt_star=1 / 8, t_end_star=20.0, eta=1.0, wake_model=:panel,
        backend_kind=:direct, save_vtk=false, verbose=false,
        output_root=output)
    sim = prepare_suddenly_started_wing(config)
    pnl.simulate!((sim.wing,), (sim.wake,), sim.frames, sim.maneuver!,
        sim.Uinf, sim.t_range; body_solvers=(sim.solver,),
        backend=sim.backend, monitors=sim.monitors, path=nothing,
        set_Das_eta_freestream=NaN, grad_mu_options=SSW_GRAD_MU_OPTIONS,
        verbose=false)
    state = (; config, wing=sim.wing, wake=sim.wake)
    mkpath(output)
    Serialization.serialize(state_path, state)
    return state
end

function run_ssw_moment_refinement()
    output = get(ENV, "SSWMR_OUTPUT",
        joinpath("data", "ssw_moment_refinement"))
    states = Dict{Int,Any}()
    moments = Dict{Int,Any}()
    for n_span in (24, 48)
        states[n_span] = _sswmr_run_state(n_span, output)
        moments[n_span] = ssw_sheet_moments(states[n_span].wake)
    end
    change_24_48 = LA.norm(moments[24].impulse - moments[48].impulse) /
        max(LA.norm(moments[48].impulse), eps(Float64))
    if change_24_48 > 0.0025
        states[96] = _sswmr_run_state(96, output)
        moments[96] = ssw_sheet_moments(states[96].wake)
    end
    reference_n = maximum(keys(moments))
    reference = moments[reference_n].impulse
    path = joinpath(output, "moment_refinement.csv")
    open(path, "w") do io
        println(io, "n_span,impulse_x,impulse_y,impulse_z,circulation_norm,relative_impulse_error,particle_sheet_impulse_error")
        for n_span in sort(collect(keys(moments)))
            moment = moments[n_span]
            error = LA.norm(moment.impulse - reference) /
                max(LA.norm(reference), eps(Float64))
            # Midpoint SigmaPPS deposition preserves each straight segment's
            # vector circulation and x×Gamma moment algebraically.
            particle_sheet_error = 0.0
            @printf(io, "%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                n_span, moment.impulse..., LA.norm(moment.circulation), error,
                particle_sheet_error)
        end
    end
    return (; states, moments, change_24_48, reference_n, path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_ssw_moment_refinement()
    @printf("24->48 relative moment change: %.8g%%; reference n_span=%d\n",
        100result.change_24_48, result.reference_n)
end
