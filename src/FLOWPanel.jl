"""
    FLOWPanel

Three-dimensional panel-method toolkit for non-lifting and lifting aerodynamic
surface models, wake models, solver backends, and post-processing utilities.
"""
module FLOWPanel

export  solve, save, influence!,
        ILUPreconditioner,
        get_ndivscells, get_ndivsnodes,
        get_cart2lin_cells, get_cart2lin_nodes,
        get_field, get_fieldval, add_field,
        iswatertight,
        ensure_consistent_winding, ensure_consistent_winding!,
        calc_normals!, calc_normals,
        calc_tangents!, calc_tangents,
        calc_obliques!, calc_obliques,
        calc_controlpoints!, calc_controlpoints,
        calc_areas!, calc_areas,
        calc_shedding, calc_shedding_from_seed, trace_trailing_edge,
        meshes2nodes_cells, read_gmsh,
        PressureBernoulli, PressureLaplace,
        JacobiPressurePreconditioner, NoPressurePreconditioner,
        IncompleteCholeskyPressurePreconditioner, AMGPressurePreconditioner,
        ForceMonitor, KuttaJoukowskiForce, SurfaceVorticityForce,
        BoundCirculationMonitor, SpanwiseLoadingMonitor, DragPolarMonitor,
        WingNormalization, NoNormalization, RotorNormalization,
        NoSectionalNormalization, FreestreamSectionalNormalization,
        RotorSectionalNormalization,
        steady!, simulate_warmstart!, initialize_Das!,
        AbstractSolveFormulation, VelocityThroughSources,
        GreenReconstruction, HybridWakePotential, TraceCorrected,
        DirectWakePotential,
        ParticleBodyOverlapPolicy, ParticleBodyOverlapReport,
        ParticleBodyOverlapError, particle_body_overlap,
        check_particle_body_overlap!,
        set_wake_correction!, clear_wake_correction!,
        AbstractWakeAttachment, RigidTransitionAttachment, TEAnchoredAttachment,
        AbstractKuttaClosure, JumpKutta, PressureContinuityKutta,
        AbstractKuttaPressureProvider, SteadyBernoulliProvider,
        kutta_diagnostics, KuttaDiagnostics, KuttaConvergenceError,
        replay, ReplayResult, migrate_metadata_toml

# ------------ GENERIC MODULES -------------------------------------------------
import LinearAlgebra as LA
import LinearAlgebra: I, lu!, ldiv!
import Krylov
import LinearOperators
import SparseArrays
import Requires: @require
# import SimpleNonlinearSolve
import FastMultipole
import ILUZero
using FastMultipole.StaticArrays: @SVector, SVector, SMatrix
import ReadVTK
using WriteVTK
import Meshes
import TOML

# ------------ FLOW LAB MODULES ------------------------------------------------
import ImplicitAD as IAD
import ImplicitAD: ForwardDiff as FD, ReverseDiff as RD
import FLOWMath as math
import FLOWVPM
import VSPGeom

# ------------ GLOBAL VARIABLES AND DATA STRUCTURES ----------------------------
const module_path = splitdir(@__FILE__)[1]      # Path to this module
                                                # Default path to data files
const def_data_path = joinpath(module_path, "..", "docs", "resources", "data")
                                            # Default path to airfoil data files
const def_rfl_path = joinpath(def_data_path, "airfoils")
                                                # Path to examples
const examples_path = joinpath(module_path, "..", "examples")

# 1/4pi
const ONE_OVER_4PI = 1.0 / (4.0 * pi)

# Discretization parameter type
const ndivstype = Union{Float64, AbstractVector, Nothing}

# Shedding matrix for a RigidWakeBody without shedding
const noshedding = zeros(Int, 6, 0)

const SEMIINFINITE_LENGTH = Array{Float64,0}(undef)
SEMIINFINITE_LENGTH[] = 10.0

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["elements", "fmm",
                    "abstractbody", "nonliftingbody",
                    "abstractliftingbody", "liftingbody",
                    "instrumentation", "solver",
                    "elements_fmm", "frames",
                    "liftingline",
                    "utils", "postprocess",
                    "wake", "gpu_influence", "gpu_wake", "particle_body_overlap",
                    "formulation", "kutta", "simulate_monitors", "simulate_monitors_fieldprobe", "metadata", "simulate", "warmstart",
                    "replay",
                    ]
  include("FLOWPanel_"*header_name*".jl")
end

const DEBUG = Array{Bool,0}(undef)
DEBUG[] = false

"""
    __init__()

Load optional runtime integrations, including plotting monitors when `PythonPlot`
is available.
"""
function __init__()

    # Conditionally load monitors if PythonPlot is available
    try
        @require PythonPlot="274fc56d-3b97-40fa-a1cd-1b4a50311bf9" begin

            import .PythonPlot as plt

            for header_name in ["monitor"]
              include("FLOWPanel_"*header_name*".jl")
            end

        end

    catch e
        @warn "PythonPlot is not available; monitors will not be loaded"
    end

    # BRAINSTORM 025: pin the filament regularization family from the
    # environment so frozen drivers can select it without code changes
    # (compact/gaussian/vatistas; default compact).
    if haskey(ENV, "FLOWPANEL_FILAMENT_REG")
        set_filament_regularization!(Symbol(lowercase(ENV["FLOWPANEL_FILAMENT_REG"])))
        # Log provenance: nothing else in a job log identifies the family, and
        # the two families are only distinguishable after the fact by the FMM
        # radius-inflation warning (37.6*rc vs 5.9*rc).
        println("FLOWPanel: filament regularization = $(FILAMENT_REGULARIZATION[]) ",
                "(pinned by FLOWPANEL_FILAMENT_REG)")
    end

end



end # END OF MODULE
