"""
    FLOWPanel

Three-dimensional panel-method toolkit for non-lifting and lifting aerodynamic
surface models, wake models, solver backends, and post-processing utilities.
"""
module FLOWPanel

export  solve, save, influence!,
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
        meshes2nodes_cells,
        PressureBernoulli, PressureLaplace,
        JacobiPressurePreconditioner, NoPressurePreconditioner,
        IncompleteCholeskyPressurePreconditioner, AMGPressurePreconditioner,
        ForceMonitor, KuttaJoukowskiForce,
        WingNormalization, NoNormalization, RotorNormalization,
        simulate_warmstart!, initialize_Das!

# ------------ GENERIC MODULES -------------------------------------------------
import LinearAlgebra as LA
import LinearAlgebra: I, lu!, ldiv!
import Krylov
import LinearOperators
import SparseArrays
import Requires: @require
# import SimpleNonlinearSolve
import FastMultipole
using FastMultipole.StaticArrays: @SVector, SVector, SMatrix
using WriteVTK
import Meshes

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
                    "solver",
                    "elements_fmm", "frames",
                    "liftingline",
                    "utils", "postprocess",
                    "wake", "simulate_monitors", "simulate", "warmstart",
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

end



end # END OF MODULE
