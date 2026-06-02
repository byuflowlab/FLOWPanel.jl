#=##############################################################################
# DESCRIPTION
    Fan duct replicating the experiment reported by V. P. Hill (1978), "A
    Surface Vorticity Theory for Propeller Ducts and Turbofan Engine Cowls in
    Non-Axisymmetric Incompressible Flow." The same experiment is also
    discussed in Sections 4.4 and 6.3.1 of Lewis, R. (1991), "Vortex Element
    Methods for Fluid Dynamic Analysis of Engineering Systems."

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jan 2023
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import GeometricTools as gt
include(joinpath(pnl.examples_path, "helper_functions.jl"))
import CSV
import DataFrames: DataFrame

include(joinpath(pnl.examples_path, "duct_postprocessing_pyplot.jl"))

run_name        = "duct-hill00"             # Name of this run

save_path       = joinpath("data", run_name) # Where to save outputs
fluiddomain     = false                     # Whether to generate fluid domain
paraview        = true                      # Whether to visualize with Paraview
call_paraview   = false                     # Whether to call Paraview at the end

save_plots      = false                     # Whether to save plots or not
# Where to save plots (default to re-generating the figures that are used
# in the docs)
fig_path = joinpath(pnl.examples_path, "..", "docs", "resources", "images")


# ----------------- SIMULATION PARAMETERS --------------------------------------
AOAs            = [0, 5, 15]                # (deg) angles of attack to evaluate
# AOAs            = [5]                # (deg) angles of attack to evaluate
magVinf         = 30.0                      # (m/s) freestream velocity
rho             = 1.225                     # (kg/m^3) air density


# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
# Read duct contour (table in figure 7.4 of Lewis 1991)
filename        = joinpath(pnl.examples_path, "data", "naca662015.csv")
contour         = CSV.read(filename, DataFrame)

aspectratio     = 0.6                       # Duct trailing edge aspect ratio l/d
d               = 2*0.835                   # (m) duct diameter


# ----------------- SOLVER PARAMETERS ------------------------------------------
# Discretization
NDIVS_theta     = 80                        # Number of azimuthal panels
# NDIVS_theta     = 20                        # Number of azimuthal panels

# NOTE: NDIVS is the number of divisions (panels) in each dimension. This can be
#       either an integer, or an array of tuples as shown below

n_rfl           = 8                        # This controls the number of chordwise panels
# n_rfl           = 6                        # This controls the number of chordwise panels

NDIVS_rfl_up = [                            # Discretization of airfoil upper surface
            # 0 to 0.25 of the airfoil has `n_rfl` panels at a geometric expansion of 10 that is not central
                # (0.25, n_rfl >> 1,   1.0, false),
                (0.25, n_rfl,   10.0, false),
            # 0.25 to 0.75 of the airfoil has `n_rfl` panels evenly spaced
                (0.50, n_rfl,    1.0, true),
            # 0.75 to 1.00 of the airfoil has `n_rfl` panels at a geometric expansion of 0.1 that is not central
                # (0.25, n_rfl >> 1,    1.0, false)]
                (0.25, n_rfl,    0.1, false)]

NDIVS_rfl_lo = NDIVS_rfl_up                 # Discretization of airfoil lower surface

# NOTE: A geometric expansion of 10 that is not central means that the last
#       panel is 10 times larger than the first panel. If central, the
#       middle panel is 10 times larger than the peripheral panels.

# Solver: Vortex-ring least-squares
# bodytype        = pnl.RigidWakeBody{pnl.VortexRing, 1, Float64} # Elements and wake model
# bodytype        = pnl.RigidWakeBody{pnl.ConstantDoublet, 1, Float64} # Elements and wake model
# kernel = Union{pnl.ConstantSource, pnl.ConstantDoublet}
kernel = Union{pnl.ConstantSource, pnl.VortexRing}
# kernel = pnl.VortexRing
bodytype = pnl.RigidWakeBody{kernel} # Elements and wake model


# ----------------- GENERATE BODY ----------------------------------------------
# Re-discretize the contour of the body of revolution according to NDIVS
xs, ys = gt.rediscretize_airfoil(contour[:, 1], contour[:, 2],
                                        NDIVS_rfl_up, NDIVS_rfl_lo;
                                        verify_spline=false)

# Make sure that the contour is closed
ys[end] = ys[1]
ys .*= 1.0

# Scale contour by duct length
xs *= d*aspectratio
ys *= d*aspectratio

# Move contour to the radial position
ys .+= d/2

# Collect points that make the contour of the body of revolution
points = hcat(xs, ys)

# Generate body of revolution
body = generate_revolution_liftbody(bodytype, points, NDIVS_theta;
                                        bodyoptargs = (
                                                        kerneloffset=1e-2,
                                                        kernelcutoff=1e-14,
                                                        characteristiclength=(args...)->d*aspectratio,
                                                        semiinfinite_wake=false
                                            )
                                        )

# body.shedding[1] = body.shedding[1][:,1:1]
# body.nsheddings = 1
# body.shedding_full .= -1
# idx_shed = body.shedding[1, 1]
# body.shedding_full[1:2, idx_shed] .= view(body.shedding[1], 2:3, 1)
# body.shedding_full[3:4, idx_shed] .= 1
# idx_shed = body.shedding[1][4, 1]
# body.shedding_full[1:2, idx_shed] .= view(body.shedding[1], 5:6, 1)
# body.shedding_full[3:4, idx_shed] .= 1

println("Number of panels:\t$(body.ncells)")

vtks = save_path*"/"                        # String with VTK output files


# ----------------- CALL SOLVER ------------------------------------------------
global solver = nothing
# for (i, AOA) in enumerate(AOAs)             # Sweep over angle of attack
i = 2
AOA = AOAs[i]

    println("Solving body...")

    # Freestream vector
    Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]

    # Unitary direction of semi-infinite vortex at points `a` and `b` of each
    # trailing edge panel
    body.Das[1] .= repeat(Vinf/magVinf, 1, size(body.Das[1], 2))

    # Solve body (panel strengths) giving `Uinfs` as boundary conditions and
    # `Das` as trailing edge rigid wake direction
    # @time pnl.solve(body, Uinfs, Das)
    leaf_size = 100
    expansion_order = 10
    multipole_acceptance = 0.4
    backend = pnl.FastMultipoleBackend(;
                                    expansion_order,
                                    multipole_acceptance,
                                    leaf_size=20
                                )
    # backend = pnl.DirectBackend()
    # global solver = pnl.Backslash(body; least_squares=true)
    # solver = pnl.KrylovSolver(body;
    #     method=:gmres,
    #     itmax=20,
    #     atol=1e-4,
    #     rtol=1e-4,
    #     # elprescribe=Tuple{Int,Float64}[],   # No prescribed strengths
    #     backend=pnl.FastMultipoleBackend(
    #                 expansion_order=7,
    #                 multipole_acceptance=0.4,
    #                 leaf_size=10
    #             )
    # )
    # function test_solver(inner_iterations, reverse_pass)
    #     println("Initializing solver with inner_iterations=$inner_iterations, reverse_pass=$reverse_pass...")
        println("Initializaing solver...")
        # @time solver = pnl.FGSSolver(body;
        #     max_iterations=500,         # Maximum number of iterations
        #     tolerance=1.0e-6,            # Convergence tolerance
        #     rlx=1.0,                  # Relaxation factor
        #     expansion_order,
        #     multipole_acceptance,
        #     leaf_size=10000,
        #     shrink=true,
        #     recenter=false,
        #     inner_iterations=20,
        #     reverse_pass=false,
        #     verbose=false
        # )
        solver = pnl.Backslash(body)
        frames = pnl.ReferenceFrame(body)
        pressure_monitor = pnl.PressureBernoulli(rho; correct_kuttacondition=true)
        force_monitor = pnl.ForceMonitor(1, 1; normalization=pnl.NoNormalization())
        monitors = (pressure_monitor, force_monitor)

        name = typeof(solver) <: pnl.Backslash ? "duct_dirichlet" : "duct"
        name *= kernel == pnl.VortexRing ? "_vortexring" :
                kernel == Union{pnl.ConstantSource, pnl.ConstantDoublet} ? "_source_doublet" :
                ""

        println("\nSolving...")

        # profile
        # using Profile, PProf
    @time pnl.steady!(body, frames, Vinf;
        body_solvers=solver,
        backend=backend,
        monitors=monitors,
        i_run=1,
        dt=1.0,
        path=paraview ? save_path : nothing,
        name=name,
    )
    paraview && (global vtks *= joinpath(save_path, "$(name)_body1.pvd"))
        # @profile pnl.solve!(body, solver; backend)
        # Profile.clear()
        # @profile pnl.solve!(body, solver; backend)
        # pprof()

    # end

    # for reverse_pass in [true, false]
    #     for inner_iterations in [1, 2, 4, 8, 16, 32]
    #         test_solver(inner_iterations, reverse_pass)
    #     end
    # end

    # ----------------- POST PROCESSING ----------------------------------------
    println("\nPost processing...")

    # check boundary condition
    Us_tot = body.velocity
    Udotn = sum(Us_tot .* body.normals, dims=1)
    tangency_resid = maximum(abs.(Udotn))
    println("Max flow tangency residual: $tangency_resid")

    if pnl.has_dirichlet_bc(body)
        potential_old = copy(body.potential)
        pnl.calc_normals!(body)
        pnl.calc_controlpoints!(body)
        body.potential .= 0
        pnl.influence!(body, body, pnl.DirectBackend(); scalar_potential=true, velocity=false)
        resid = maximum(abs.(body.potential))
        body.potential .= potential_old
        println("Max interior potential residual: $resid")
    end

    # ----------------- COMPARISON TO EXPERIMENTAL DATA ------------------------
    # Plot surface pressure along slices of the duct
    fig, axs = plot_Cp(body, AOA, rho, magVinf; pressure=pressure_monitor.pressure[1])

    if save_plots
        fname = "$(run_name)-Cp-AOA$(ceil(Int, AOA)).png"
        fig.savefig(joinpath(fig_path, fname), dpi=300, transparent=true)
    end

    # ----------------- VISUALIZATION ------------------------------------------
    # Compute fluid domain and save as VTK
    if fluiddomain
        global vtks *= generate_fluiddomain(body, AOA, Vinf, d,
                                                aspectratio, save_path; num=i)
    end

# end

# Call Paraview
if paraview && call_paraview
    run(`paraview --data=$(vtks)`)
end
