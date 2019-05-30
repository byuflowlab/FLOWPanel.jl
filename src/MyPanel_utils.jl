#=##############################################################################
# DESCRIPTION
    Utilities.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Sep 2018
  * License   : AGPL-3.0
=###############################################################################


"""
`simplewing(b, ar, tr, twist_root, twist_tip, lambda, gamma;
bodytype=RigidWakeBody,
span_NDIVS="automatic", rfl_NDIVS="automatic",
airfoil_root="naca6412.dat", airfoil_tip="naca6412.dat",
airfoil_path=def_rfl_path)`

Generates a symmetric single-section wing.

**ARGUMENTS**
  * `b::Real`         : Span.
  * `ar::Real`        : Aspect ratio defined as b/c_tip.
  * `tr::Real`        : Taper ratio defined as c_tip/c_root.
  * `twist_root::Real`: (deg) twist of the root.
  * `twist_tip::Real` : (deg) twist of the tip.
  * `lambda::Real`    : (deg) sweep.
  * `gamma::Real`     : (deg) dihedral.

**OPTIONAL ARGUMENTS**
  * `bodytype::Type{LBodyTypes}`: Type of lifting body to generate.
  * `span_NDIVS::ndivstype`     : Spanwise divisions.
  * `rfl_NDIVS::ndivstype`    : Chordwise divisions.
  * `airfoil_root::String`      : File to root airfoil contour.
  * `airfoil_tip::String`       : File to tip airfoil contour.
  * `airfoil_path::String`      : Path to airfoil files.

NOTE: See gt.multidscretize for a description of arguments of type `ndivstype`.
NOTE2: In the current implementation, sweep and dihedral are done about the LE.
"""
function simplewing(b::RType, ar::RType, tr::RType, twist_root::RType,
                      twist_tip::RType, lambda::RType, gamma::RType;
                      bodytype::Type{LBodyTypes}=RigidWakeBody,
                      span_NDIVS::ndivstype=nothing,
                      rfl_NDIVS::ndivstype=nothing,
                      airfoil_root::String="naca6412.dat",
                      airfoil_tip::String="naca6412.dat",
                      airfoil_path::String=def_rfl_path,
                      spl_s::Real=0.0000001,
                      rflspl_s::Real=0.00000001,
                      verify_spline::Bool=true,
                      verify_rflspline::Bool=true,
                      opt_args...
                      )

  # ----------------- GEOMETRY DESCRIPTION -------------------------------------
  c_tip = b/ar                        # Tip chord
  c_root = c_tip/tr                   # Root chord
  semispan = b/2                      # (m) semi-span length

  y_tip = b/2
  x_tip = y_tip*tan(lambda*pi/180)
  z_tip = y_tip*tan(gamma*pi/180)

  chords = [0.00 c_root/semispan;     # (semi-span position, chord c/semib)
            1.00 c_tip/semispan]

  twists = [0.0 twist_root;           # (semi-span position, twist (deg))
            1.0 twist_tip]

  x_pos = [0.00 0;                    # (semi-span position, LE x-position x/semib)
           1.00 x_tip/semispan]

  z_pos = [0.00 0;                    # (semi-span position, LE z-position x/semib)
           1.00 z_tip/semispan]

  airfoil_files = [(0.0, airfoil_root), # (semi-span position, airfoil file)
                   (1.0, airfoil_tip)]

 # ----------------- DISCRETIZATION 0000----------------------------------------

  # Defines divisions
  if span_NDIVS==nothing
    b_NDIVS = [(1.0, 35, 20.0, true)]  # Span cell sections
  else
    b_NDIVS = span_NDIVS
  end

  if rfl_NDIVS==nothing
    urfl_NDIVS = [(0.25, 7,   10.0, false),   # Cells on upper side of airfoils
                  (0.50,  5,    1.0, true),
                  (0.25,  6, 1/10.0, false)]
  else
    urfl_NDIVS = rfl_NDIVS
  end

  lrfl_NDIVS = urfl_NDIVS             # Cells on lower side of airfoils


  # ----------------- LOFTING PARAMETERS ---------------------------------------
  b_low = -1.0                        # Lower bound of span lofting
  b_up = 1.0                          # Upper bound of span lofting
  symmetric = true                    # Lofting symmetric about b=0
  spl_k = 1                           # Spline order of distributions along span
  # spl_s = 0.0000001                 # Spline smoothing of distribution along span
  # rflspl_s = 0.00000001             # Spline smoothing of airfoil cross sections.
  # verify_spline = false             # Plots the splined distributions
  # verify_rflspline = true           # Plots the splined airfoil cross sections

  return generate_loft_liftbody(bodytype, airfoil_files, airfoil_path,
                                        urfl_NDIVS, lrfl_NDIVS,
                                        semispan, b_low, b_up, b_NDIVS,
                                        chords, twists, x_pos, z_pos;
                                        dimsplit=1,
                                        symmetric=symmetric,
                                        spl_k=spl_k, spl_s=spl_s,
                                        verify_spline=verify_spline,
                                        verify_rflspline=verify_rflspline,
                                        rflspl_s=rflspl_s,
                                        opt_args...
                                    )
end
