#=##############################################################################
# DESCRIPTION
    36deg swept-back wing with a cambered NACA 64(1)-612 airfoil from the NACA 
    RM A50K27 Flying Wing​ report. This wing has an aspect ratio of 15 and tapper 
    ratio of 0.5. The experiment was run at Re_c = 2e6 and Mach 0.27.

    This example uses the nonlinear lifting line method.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jan 2026
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import FLOWPanel: mean, norm, dot, cross

import PyPlot as plt
import PyPlot: @L_str
include(joinpath(pnl.examples_path, "plotformat.jl"))

import CSV
import DataFrames: DataFrame


run_name        = "ll-a50k27"                   # Name of this run

save_path       = run_name                      # Where to save outputs
airfoil_path    = joinpath(pnl.examples_path, "data") # Where to find 2D polars

paraview        = false                         # Whether to visualize with Paraview
save_outputs    = false                         # Whether to save outputs for docs or not


# ----------------- SIMULATION PARAMETERS --------------------------------------
magUinf         = 91.3                          # (m/s) freestream velocity
rho             = 1.225                         # (kg/m^3) air density


# ------------------ GEOMETRY PARAMETERS ---------------------------------------

# High-level parameters
b               = 2*1.5490                      # (m) wing span

# Discretization parameters
nelements       = 40                            # Number of stripwise elements per semi-span
discretization  = [                             # Multi-discretization of wing (seg length, ndivs, expansion, central)
                    (1.00,  nelements, 1/10, false),
                  ]
symmetric       = true                          # Whether the wing is symmetric

# Chord distribution (nondim y-position, nondim chord)
chord_distribution = [
#   2*y/b   c/b
    0.0     0.132989
    1.0     0.0664945
]

# Twist distribution (nondim y-position, twist)
twist_distribution = [
#   2*y/b   twist (deg)
    0.0     0.0
    1.0     0.0
]

# Sweep distribution (nondim y-position, sweep)
sweep_distribution = [
#   2*y/b   sweep (deg)
    0.0     36.2603
    1.0     36.2603
]

# Dihedral distribution (nondim y-position, dihedral)
dihedral_distribution = [
#   2*y/b   dihedral (deg)
    0.0     0.0
    1.0     0.0
]

# Span-axis distribution: chordwise point about which the wing is twisted, 
# swept, and dihedralized (nondim y-position, nondim chord-position)
spanaxis_distribution = [
#   2*y/b   x/c
    0.0     0.25
    1.0     0.25
]

# Airfoil contour distribution (nondim y-position, polar, airfoil type)
airfoil_distribution = [
#    2*y/b  polar file                                airfoil type
    (0.00, "n64_1_A612-Re0p5e6-neuralfoil180-3.csv",  pnl.SimpleAirfoil),
    (1.00, "n64_1_A612-Re0p5e6-neuralfoil180-3.csv",  pnl.SimpleAirfoil)
    # (0.00, "n64_1_A612-Re0p5e6-smooth180-4.csv",  pnl.SimpleAirfoil),
    # (1.00, "n64_1_A612-Re0p5e6-smooth180-4.csv",  pnl.SimpleAirfoil)
]

element_optargs = (;    path = airfoil_path,
                        plot_polars = true,
                        extrapolatepolar = false,   # Whether to extrapolate the 2D polars to ±180 deg
                    )


# ------------------ SOLVER PARAMETERS -----------------------------------------
deltasb         = 1.0                           # Blending distance, deltasb = 2*dy/b
deltajoint      = 1.0                           # Joint distance, deltajoint = dx/c

sigmafactor     = 0.0                           # Dragging line amplification factor (set to -1.0 for rebust post-stall method)
sigmaexponent   = 1.0                           # Dragging line amplification exponent (no effects if `sigmafactor==0.0`)

                                                # Nonlinear solver
solver          = pnl.SimpleNonlinearSolve.SimpleDFSane()              # Indifferent to initial guess, but somewhat not robust
# solver        = pnl.SimpleNonlinearSolve.SimpleTrustRegion()         # Trust region needs a good initial guess, but it converges very reliably

solver_optargs  = (; 
                    abstol = 1e-9,  
                    maxiters = 800,
                    )

align_joints_with_Uinfs = true                  # Whether to align joint bound vortices with the freestream
use_Uind_for_force = false                      # Whether to use Uind as opposed to selfUind for force postprocessing
                                                # (`true` for more accurate spanwise cd distribution, but worse integrated CD)

X0              = [0.595, 0, 0]                 # Center about which to calculate moments

cref            = 0.2472                        # (m) reference chord
Aref            = 0.957                         # (m^2) reference area
nondim          = 0.5*rho*magUinf^2*Aref        # Normalization factor


# ------------------ GENERATE LIFTING LINE -------------------------------------

ll = pnl.LiftingLine{Float64}(
                                airfoil_distribution; 
                                b, chord_distribution, twist_distribution,
                                sweep_distribution, dihedral_distribution,
                                spanaxis_distribution,
                                discretization,
                                symmetric,
                                deltasb, deltajoint, sigmafactor, sigmaexponent,
                                element_optargs,
                                plot_discretization = true,
                                )

display(ll)


# ------------------ OUTPUT GEOMETRY -------------------------------------------
if !isnothing(save_path)
    
    str = pnl.save(ll, run_name; path=save_path, debug=true) # Use `debug=true` to output the effective horseshoes

    if paraview
        run(`paraview --data=$(joinpath(save_path, str))`)
    end
    
end



# ----------------- AOA SWEEP --------------------------------------------------

# Sequence of sweeps to run
# NOTE: To help convergence and speed it up the sweep, we recommend starting
#       each sweep at 0deg AOA since the sweep steps over points using the last
#       converged solution as the initial guess
sweep_neg = range(0, -50, step=-0.5)        # Sweep from 0 into deep negative stall (-50deg)
sweep_pos = range(0, 50, step=0.5)          # Sweep from 0 into deep positive stall (50deg)

distributions = []

@time wingpolar = pnl.run_polarsweep(ll,
                            magUinf, rho, X0, cref, b; 
                            Aref,
                            aoa_sweeps = (sweep_neg, sweep_pos),
                            # sweepname = run_name,
                            # plots_path= save_outputs ? fig_path : nothing,
                            # extraplots_path = save_outputs ? fig_path : nothing,
                            output_distributions=distributions,
                            solver,
                            solver_optargs,
                            align_joints_with_Uinfs, 
                            use_Uind_for_force
                        )

# ----------------- COMPARISON TO EXPERIMENTAL DATA ---------------------------

# Where to save figures (default to re-generating the figures that are used
# in the docs)
fig_path = joinpath(pnl.examples_path, "..", "docs", "resources", "images")
outdata_path = joinpath(pnl.examples_path, "..", "docs", "resources", "data")

data_exp = Dict(
 "alpha_CL" => [-5.97426  , -4.97462  , -3.97515  , -2.92982  , -1.93035  ,
        -0.931063 ,  0.0685855,  1.06788  ,  2.06753  ,  3.06735  ,
         4.11125  ,  5.11054  ,  6.10858  ,  7.19352  ,  8.09879  ,
         9.18444  , 10.1802   , 11.1764   , 12.172    , 13.1674   ],
 "CL" => [-0.101053 , -0.0231579,  0.0526316,  0.136842 ,  0.212632 ,
         0.286316 ,  0.364211 ,  0.437895 ,  0.515789 ,  0.595789 ,
         0.663158 ,  0.736842 ,  0.795789 ,  0.814737 ,  0.844211 ,
         0.871579 ,  0.903158 ,  0.941053 ,  0.970526 ,  0.997895 ],
 "alpha_CD" => [-3.97515  , -2.92982  , -1.93035  , -0.931063 ,  0.0685855,
         1.06788  ,  2.06753  ,  3.06735  ,  4.11125  ,  5.11054  ,
         6.10858  ,  8.09879  ,  9.18444  , 10.1802   , 11.1764   ,
        12.172    , 13.1674   ],
 "CD" => [0.00834592, 0.0083899 , 0.00933814, 0.0102853 , 0.010778  ,
        0.0135382 , 0.0149386 , 0.0186075 , 0.022723  , 0.0272962 ,
        0.0314071 , 0.0386831 , 0.0436847 , 0.0495952 , 0.0573199 ,
        0.0654955 , 0.0745777 ],
 "polar_CD" => [0.00834592, 0.0083899 , 0.00933814, 0.0102853 , 0.010778  ,
        0.0135382 , 0.0149386 , 0.0186075 , 0.022723  , 0.0272962 ,
        0.0314071 , 0.0386831 , 0.0436847 , 0.0495952 , 0.0573199 ,
        0.0654955 , 0.0745777 , 0.0890969 , 0.103164  , 0.124479  ],
 "polar_CL" => [0.0512156, 0.133679 , 0.211938 , 0.288083 , 0.362101 , 0.438293 ,
        0.51445  , 0.594895 , 0.662665 , 0.738906 , 0.798218 , 0.842813 ,
        0.872548 , 0.906535 , 0.942684 , 0.974617 , 1.00657  , 1.03444  ,
        1.06442  , 1.08613  ],
 "alpha_Cm" => [-6.03774562, -5.01494306, -3.93954365, -2.9277364 , -1.89823568,
        -0.88378519,  1.12056138,  2.07016703,  3.0706195 ,  4.04682072,
         5.09912094,  6.09522133,  6.80776737,  7.26302621,  7.78516742,
         9.16147189, 10.17521789, 11.27594987, 12.26197843, 13.60765419],
 "Cm" => [-0.0997 , -0.0989 , -0.106  , -0.109  , -0.111  , -0.117  ,
        -0.127  , -0.132  , -0.135  , -0.142  , -0.14   , -0.132  ,
        -0.104  , -0.0988 , -0.0772 , -0.0588 , -0.0412 , -0.0213 ,
         0.00109,  0.0283 ]
)


alphas_vsp = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
CLs_vsp = [-0.132535789, -0.048719081, 0.03513067, 0.11897715, 0.202784821, 0.28652272, 0.37014917, 0.453635042, 0.536942698, 0.620034102, 0.702887187, 0.785460061, 0.867717773, 0.949625582, 1.031153268, 1.112271021]
CDs_vsp = [0.005040841, 0.004834691, 0.005176697, 0.00606653, 0.007503003, 0.009483998, 0.012006806, 0.015067818, 0.018662681, 0.022785756, 0.027433557, 0.032597544, 0.038270588, 0.044444463, 0.051109996, 0.058258825]
CMs_vsp = [-0.102929546, -0.106897688, -0.110862096, -0.114829482, -0.118805477, -0.122792914, -0.126795469, -0.130816755, -0.134863039, -0.138938458, -0.143070238, -0.147244647, -0.15145668, -0.155719011, -0.160031645, -0.164409958]


data_machupx = Dict(
  "alphas" => [ -20.0,   -18.0,   -17.0,   -14.0,   -13.0,   -12.0,   -11.0,   -10.0,   -8.0,   -7.0,   -6.0,   -5.0,   -4.0,   -3.0,   -2.0,   -1.0,   0.0,   1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,   10.0,   11.0,   12.0,   13.0,   14.0,   15.0,   16.0,   17.0,   18.0,   19.0,   20.0,   21.0,   22.0],  
  "CLs" => [-0.30003383558109326,  -0.2902681040145641,  -0.28610495445332534,  -0.27687304695108234,  -0.2749120028891531,  -0.2734109184885327,  -0.27240901926892036,  -0.27196447651320677,  -0.2673755224982744,  -0.18303990030065215,  -0.10183740808378725,  -0.02363792533124699,  0.05437467577504556,  0.13167438022307004,  0.20853532935968153,  0.28134457745489255,  0.35852821356970593,  0.4367245557511367,  0.5130567541179603,  0.5864953578755846,  0.6493400375973891,  0.6928164717758137,  0.7377429573675656,  0.7868395980290865,  0.8330113316368875,  0.8737841435173891,  0.9032960968357137,  0.9247053155750491,  0.857600343071871,  0.757979929170746,  0.8157990884362425,  0.7824552554903349,  0.7059002499406142,  0.7150745067248339,  0.6626079262914916,  0.6133467364380408,  0.5829416037848351,  0.5686701910144,  0.5755667192036986],
  "CDs" => [-0.10564349722643289,   -0.0823926267933612,   -0.07056990671319108,   -0.03439788828978355,   -0.022224930511816423,   -0.01002427047332819,   0.0022197133444258757,   0.014212348322147114,   0.023374706192339077,   0.011662760592502942,   0.0018851163129974388,   -0.00552793358016962,   -0.010782580828240034,   -0.013794138464614412,   -0.01447578680892808,   -0.012163377978522086,   -0.008087369278167358,   -0.0019857515410504046,   0.0063823769476897664,   0.016933729508596073,   0.02931486156322386,   0.0423988847557252,   0.05597009079198652,   0.07113056036549091,   0.08766278771684907,   0.10514094675569198,   0.12223182001600856,   0.13861915054905632,   0.12426478739779515,   0.09741488564584272,   0.10510526410014874,   0.08611859608219757,   0.05175824231820578,   0.04558778268635248,   0.05253172362403583,   0.03160724540922949,   0.03031946591029029,   0.010226533091251348,   -0.0044513871151397565],
  "CDis" => [0.08273385325937993,  0.07592844340062109,  0.07240188568036963,  0.06126934891875641,  0.05736603337255056,  0.05336912378051757,  0.049277607569836276,  0.0450880759373625,  0.03554995514738124,  0.021490657329537907,  0.010351835876668216,  0.001997392522781464,  -0.00391299949011832,  -0.007360987147308098,  -0.008391799747112516,  -0.006907331670120758,  -0.003246244183673866,  0.0028114805723793916,  0.011275619257535485,  0.022047367996362366,  0.034788940152968115,  0.0485551164496591,  0.06367728353456649,  0.08067157455254315,  0.09907346462592805,  0.11840998755214524,  0.13767998213199797,  0.15676336872091126,  0.15736640832628007,  0.14572241036656586,  0.16789918746531055,  0.1690022704413281,  0.16053044285034881,  0.17200687148219082,  0.17594513589352917,  0.17523696868212724,  0.1825416137557727,  0.1845021956891498,  0.1924372361464497],
  "CDvs" => [-0.18837735048581283,  -0.1583210701939823,  -0.14297179239356073,  -0.09566723720853997,  -0.07959096388436698,  -0.06339339425384576,  -0.0470578942254104,  -0.03087572761521538,  -0.012175248955042163,  -0.009827896737034967,  -0.008466719563670776,  -0.007525326102951085,  -0.006869581338121714,  -0.006433151317306313,  -0.006083987061815564,  -0.005256046308401327,  -0.004841125094493493,  -0.004797232113429796,  -0.004893242309845718,  -0.005113638487766292,  -0.005474078589744254,  -0.006156231693933893,  -0.007707192742579967,  -0.009541014187052237,  -0.011410676909078968,  -0.013269040796453252,  -0.01544816211598942,  -0.01814421817185493,  -0.03310162092848491,  -0.04830752472072315,  -0.06279392336516179,  -0.08288367435913051,  -0.10877220053214305,  -0.12641908879583832,  -0.12341341226949334,  -0.14362972327289775,  -0.15222214784548244,  -0.17427566259789845,  -0.19688862326158946],
  "CDcorrs" => [ 0.08273385,  0.07592844,  0.07240189,  0.06126935,  0.05736603,  0.05336912,  0.04927761,  0.04508808,  0.03554996,  0.02149066,  0.01035184,  0.00199739, -0.003913  , -0.00736099, -0.0083918 , -0.00690733, -0.00324624,  0.00281148,  0.01127562,  0.02204737,  0.03478894,  0.04855512,  0.06367728,  0.08067157,  0.09907346,  0.11840999,  0.13767998,  0.15676337,  0.15736641,  0.14572241,  0.16789919,  0.16900227,  0.16053044,  0.17200687,  0.17594514,  0.17523697,  0.18254161,  0.1845022 ,  0.19243724],
  "CMs" => [0.6358704591633543,  0.6054636116463847,  0.5917037414357292,  0.5575275928027122,  0.5485602486310102,  0.5404146832447573,  0.5332419707637415,  0.5271883281918367,  0.5125902369395189,  0.3323469735834919,  0.15696909460046493,  -0.013846475451795884,  -0.18477177181865756,  -0.3542716708459615,  -0.5228330408579934,  -0.6817114144583686,  -0.8521508434800384,  -1.02373634298832,  -1.1907361779698473,  -1.350835629815154,  -1.4812654732443986,  -1.5654614266502935,  -1.663505764329214,  -1.767964845934618,  -1.8612519851803104,  -1.9425194500253369,  -1.9959500355371413,  -2.0328359627175403,  -1.9351958266348908,  -1.6449069528956746,  -1.7413198238661531,  -1.7636818189826222,  -1.5913676348351178,  -1.588472952497899,  -1.4618720268610268,  -1.4005731695687769,  -1.2921961375068662,  -1.2995837911091752,  -1.3184911465730837]
 )

data_aerosandbox = Dict(
 "CL" => [-0.35422358, -0.33730102, -0.30719094, -0.28093264, -0.27510972,
        -0.33436453, -0.36700134, -0.3900839 , -0.38643296, -0.36618758,
        -0.31674076, -0.26238666, -0.2048178 , -0.14499605, -0.08368502,
        -0.02119285,  0.0423122 ,  0.10659571,  0.17042054,  0.23102446,
         0.29405869,  0.35785676,  0.42143149,  0.48376024,  0.5440495 ,
         0.60336228,  0.65490008,  0.6952417 ,  0.73568179,  0.77434743,
         0.80922296,  0.84157319,  0.87015225,  0.89309723,  0.90884231,
         0.9158099 ,  0.9114969 ,  0.88579377,  0.84636064,  0.78012818,
         0.70374976,  0.63993396,  0.60217408,  0.53360447,  0.51929955],
 "CD" => [0.24706945, 0.21904496, 0.19666945, 0.178282  , 0.15409159,
        0.11881611, 0.08491982, 0.05463231, 0.03666856, 0.02035821,
        0.01670138, 0.01388506, 0.01165345, 0.00992104, 0.00863956,
        0.00776366, 0.00727488, 0.00718515, 0.00734891, 0.0072467 ,
        0.00805719, 0.0094544 , 0.01121247, 0.01328725, 0.01563414,
        0.01829329, 0.02157101, 0.0255072 , 0.02918074, 0.03274226,
        0.03649641, 0.04041473, 0.0446092 , 0.04924513, 0.0545536 ,
        0.06107782, 0.06977351, 0.09049961, 0.10493688, 0.16128909,
        0.2281551 , 0.26871724, 0.32676672, 0.36345343, 0.41702237],
 "CDi" => [4.12748760e-03, 4.34637846e-03, 3.68465744e-03, 3.13765609e-03,
        3.08879031e-03, 6.21358790e-03, 6.29765822e-03, 6.24096081e-03,
        4.95168542e-03, 4.45274242e-03, 3.35352894e-03, 2.31069287e-03,
        1.41258677e-03, 7.09563642e-04, 2.36474940e-04, 1.50110766e-05,
        6.15401802e-05, 3.88666879e-04, 9.90968000e-04, 1.81992854e-03,
        2.95193907e-03, 4.37433785e-03, 6.06985393e-03, 7.99702359e-03,
        1.01080712e-02, 1.24304375e-02, 1.46548734e-02, 1.65291752e-02,
        1.84391062e-02, 2.03560785e-02, 2.22002625e-02, 2.39827477e-02,
        2.56223673e-02, 2.70068809e-02, 2.80474715e-02, 2.86740331e-02,
        2.88096523e-02, 2.56102997e-02, 2.74783112e-02, 2.92649107e-02,
        2.78855689e-02, 2.40495746e-02, 2.33628138e-02, 1.54362442e-02,
        1.53842345e-02],
 "CDp" => [0.24294196, 0.21469858, 0.19298479, 0.17514434, 0.1510028 ,
        0.11260252, 0.07862216, 0.04839135, 0.03171688, 0.01590547,
        0.01334785, 0.01157437, 0.01024086, 0.00921147, 0.00840309,
        0.00774865, 0.00721334, 0.00679648, 0.00635794, 0.00542677,
        0.00510525, 0.00508006, 0.00514262, 0.00529023, 0.00552607,
        0.00586286, 0.00691614, 0.00897803, 0.01074163, 0.01238618,
        0.01429615, 0.01643199, 0.01898684, 0.02223825, 0.02650613,
        0.03240379, 0.04096386, 0.06488931, 0.07745856, 0.13202418,
        0.20026953, 0.24466767, 0.30340391, 0.34801718, 0.40163813],
 "Cm" => [ 0.08195622,  0.03186318,  0.00354558, -0.01579885, -0.04720275,
        -0.1038634 , -0.16222441, -0.14897749, -0.12374406, -0.09615987,
        -0.0960067 , -0.09839135, -0.10234988, -0.10760644, -0.11381567,
        -0.12067485, -0.12791346, -0.13523181, -0.14200872, -0.14740805,
        -0.15531066, -0.16268802, -0.16966138, -0.17537821, -0.17972784,
        -0.18458661, -0.18224602, -0.17237972, -0.16630964, -0.15996378,
        -0.1531454 , -0.14527834, -0.1359593 , -0.12505893, -0.1124417 ,
        -0.09883905, -0.08552482, -0.07218549, -0.07353346, -0.09565867,
        -0.12051599, -0.13363426, -0.14857782, -0.23194069, -0.24902935],
 "alpha" => [-20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10,  -9,  -8,
         -7,  -6,  -5,  -4,  -3,  -2,  -1,   0,   1,   2,   3,   4,   5,
          6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,
         19,  20,  21,  22,  23,  24]
)


stl_exp = "-o"
fmt_exp = (; label="Experimental", color="k")

stl_asb = "-^"
fmt_asb = (; label="AeroSandbox LL", color="green", markersize=3, alpha=0.2)

stl_mux = "-s"
fmt_mux = (; label="MachUpX LL", color="orchid", markersize=3, alpha=0.3)

stl_vsp = "-o"
fmt_vsp = (; label="VSPAERO Panel", color="goldenrod", markersize=4, alpha=0.6)

stl_ll = ".-"
fmt_ll = (; label="FLOWPanel LL", color="steelblue", markersize=4, alpha=1.0)

function plot_polars(wingpolar; suffix="CLCDCm", 
                        aoalims=[-23, 25], aoaticks=-20:10:25,
                        CLlims=[-0.75, 1.5], CLticks=-0.5:0.5:CLlims[end],
                        CDlims=[0, 0.125], CDticks=CDlims[1]:0.025:CDlims[end],
                        Cmlims=[-0.3, 0.1], Cmticks=Cmlims[1]:0.1:Cmlims[end]
                        )
    AOAs = wingpolar.AOA
    CLs = wingpolar.CL
    CDs = wingpolar.CD
    Cms = wingpolar.Cm

    fig = plt.figure(figsize=[7*2, 5*2*0.75]*2/3)
    axs = fig.subplots(2, 2)

    axs = [axs[j, i] for i in 1:size(axs, 1), j in 1:size(axs, 2)]

    # CL vs AOA
    ax = axs[1]
    ax.plot(data_exp["alpha_CL"], data_exp["CL"], stl_exp; fmt_exp...)
    ax.plot(alphas_vsp, CLs_vsp, stl_vsp; fmt_vsp...)
    ax.plot(data_aerosandbox["alpha"], data_aerosandbox["CL"], stl_asb; fmt_asb...)
    ax.plot(data_machupx["alphas"], data_machupx["CLs"], stl_mux; fmt_mux...)
    ax.plot(AOAs, CLs, stl_ll; fmt_ll...)

    if !isnothing(CLlims)
        ax.set_ylim(CLlims)
    end
    if !isnothing(CLticks)
        ax.set_yticks(CLticks)
    end
    ax.set_ylabel(L"Lift coeff. $C_L$")

    # CD vs AOA
    ax = axs[2]
    ax.plot(data_exp["alpha_CD"], data_exp["CD"], stl_exp; fmt_exp...)
    ax.plot(alphas_vsp, CDs_vsp, stl_vsp; fmt_vsp...)
    ax.plot(data_aerosandbox["alpha"], data_aerosandbox["CD"], stl_asb; fmt_asb...)
    # ax.plot(data_machupx["alphas"], data_machupx["CDs"], stl_mux; fmt_mux...)
    ax.plot(AOAs, CDs, stl_ll; fmt_ll...)

    if !isnothing(CDlims)
        ax.set_ylim(CDlims)
    end
    if !isnothing(CDticks)
        ax.set_yticks(CDticks)
    end
    ax.set_ylabel(L"Drag coeff. $C_D$")

    # CL vs CD
    ax = axs[3]
    ax.plot(data_exp["polar_CD"], data_exp["polar_CL"], stl_exp; fmt_exp...)
    ax.plot(CDs_vsp, CLs_vsp, stl_vsp; fmt_vsp...)
    ax.plot(data_aerosandbox["CD"], data_aerosandbox["CL"], stl_asb; fmt_asb...)
    # ax.plot(data_machupx["CDs"], data_machupx["CLs"], stl_mux; fmt_mux...)
    ax.plot(CDs, CLs, stl_ll; fmt_ll...)

    if !isnothing(CLlims)
        ax.set_ylim(CLlims)
    end
    if !isnothing(CLticks)
        ax.set_yticks(CLticks)
    end
    ax.set_ylabel(L"Lift coeff. $C_L$")
    if !isnothing(CDlims)
        ax.set_xlim(CDlims)
    end
    if !isnothing(CDticks)
        ax.set_xticks(CDticks)
    end
    ax.set_xlabel(L"Drag coeff. $C_D$")

    # Cm vs AOA
    ax = axs[4]
    ax.plot(data_exp["alpha_Cm"], data_exp["Cm"], stl_exp; fmt_exp...)
    ax.plot(alphas_vsp, CMs_vsp, stl_vsp; fmt_vsp...)
    ax.plot(data_aerosandbox["alpha"], data_aerosandbox["Cm"], stl_asb; fmt_asb...)
    # ax.plot(data_machupx["alphas"], data_machupx["CMs"], stl_mux; fmt_mux...)
    ax.plot(AOAs, Cms, stl_ll; fmt_ll...)

    ax.set_ylim(Cmlims)
    ax.set_yticks(Cmticks)
    ax.set_ylabel(L"Pitching moment $C_m$")

    for (axi, ax) in enumerate(axs)

        if axi != 3
            ax.set_xlim(aoalims)
            ax.set_xticks(aoaticks)
            ax.set_xticklabels(["$val"*L"^\circ" for val in aoaticks])
            ax.set_xlabel(L"Angle of attack $\alpha$")
        end

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        ax.legend(loc="best", frameon=false, fontsize=8, reverse=true)
    end

    fig.tight_layout()

    if save_outputs
        fig.savefig(joinpath(fig_path, "$(run_name)-sweep-$(suffix).png"),
                                                    dpi=300, transparent=true)
    end

end

plot_polars(wingpolar; suffix="CLCDCm",
                        CLlims=[-0.5, 1.25], CLticks=-0.5:0.5:1.25,
                        Cmlims=[-0.2, 0.1], Cmticks=-0.2:0.1:0.1
                        )

plot_polars(wingpolar; suffix="CLCDCm-zoomout", 
                        aoalims=[-60, 60], aoaticks=-60:20:60,
                        CLlims=[-1.0, 1.5], CLticks=-1.0:0.5:1.5,
                        CDlims=[0, 0.4], CDticks=0:0.1:0.4,
                        Cmlims=[-0.5, 0.25], Cmticks=-0.5:0.25:0.25
                        )