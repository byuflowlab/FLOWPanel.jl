#=##############################################################################
# DESCRIPTION
    Utilities.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jul 2018
  * License   : AGPL-3.0
=###############################################################################

"""
  Generates a lofted wings

  **Arguments**
  * `bscale::Float64`         : Semi-span scale.
  * `b_low::Flaot64`          : Scaled low bound of the span.
  * `b_up::Flaot64`           : Scaled upper bound of the span. To generate
                                a symmetric wing, give it b_low=-1, b_up=1; for
                                a semi-span, give it b_low=0, b_up=1. All
                                y/bscale value in the following arguments must
                                go from 0 to 1.
  * `NDIVS`                   : Number of division along the arc-length of every
                                cross section (airfoil) and along the span. Give
                                it an array with two values for a uniform
                                grid, or sections formatted for
                                GeometricTools.multidiscretize() for a complex
                                refining.
  * `chords::Array{Float64,2}`: Chord distribution along the span in the form
                                [(y/bscale, c/bscale)].
  * `twists::Array{Float64,2}`: Twist distribution along the span in the form
                                [(y/bscale, deg)].
  * `LE_x::Array{Float64,2}`  : x-position (chordwise) of leading edge along the
                                span in the form [(y/bscale, x/bscale)].
  * `LE_z::Array{Float64,2}`  : z-position (dihedral-wise) of leading edge along
                                the span in the form [(y/bscale, z/bscale)].
  * `airfoils`                : Airfoil cross sections along the span in the
                                form [(y/bscale, file_name)], where `file_name`
                                the String to the file that contains the points
                                of the section. Airfoil points must start
                                at the trailing edge, go around the top side
                                to the leading edge, and back to the trailing
                                edge around the bottom side.

  **Optional Arguments**
  * `tilt_z::{Array{Float64,2}}`            : Tilting about the z-axis of
                                              every span cross section in the
                                              form [(y/bscale, deg)].
  * `spline_k`, `spline_bc`, `spline_s`     : Spline parameters.

"""
function generate_wing(bscale::Real, b_low::Real, b_up::Real, NDIVS,
                        chords::Array{T,2}, twists::Array{T,2},
                        LE_x::Array{T,2}, LE_z::Array{T,2},
                        airfoils::Array{Tuple{Float64,String}, 1};
                        # MORE GEOMETRIC OPTIONS
                        tilt_z=nothing,
                        # SPLINE OPTIONS
                        spline_k::Int64=5, spline_bc::String="extrapolate",
                        spline_s::Real=0.001, verify_spline::Bool=true,
                        # AIRFOIL FILE OPTIONS
                        rfl_header_len::Int64=0, rfl_delim::String=",",
                        verify_rfl::Bool=true, zoom_factor::Real=1.0,
                        side_legend=(1.15, 1.1),
                        rfl_s::Real=0.00001, rfl_k="automatic",
                        # OUTPUT OPTIONS
                        save_path=nothing, paraview::Bool=true,
                        file_name::String="mywing"
                       ) where{T<:Real}

  if size(NDIVS,1)!=2
    error("Invalid divisions NDIVS."*
                " Expected two dimensions, got $(size(NDIVS,1)).")

  elseif b_low>=b_up
    error("Invalid span bounds (b_low>=b_up). b_low=$blow, b_up=$b_up")

  end

  # ----------------- PARAMETRIC GRID ------------------------------------------
  P_min = [0, b_low, 0]            # Lower boundary arclength, span, dummy
  P_max = [1, b_up, 0 ]            # Upper boundary arclength, span, dummy
  loop_dim = 1                     # Loops the arclength dimension

  # Adds dummy division
  multidivtype = Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}
  if typeof(NDIVS)==multidivtype
    _NDIVS = vcat(NDIVS, [[(1.0, 0, 0.0, false)]])
  elseif typeof(NDIVS)==Array{Int64, 1}
    _NDIVS = vcat(NDIVS, 0)
  else
    error("Invalid NDIVS type $(typeof(NDIVS))."*
          " Expected $(Array{Int64,1}) or $multidivtype.")
  end

  grid = gt.Grid(P_min, P_max, _NDIVS, loop_dim)


  # ----------------- GEOMETRY SPLINES -----------------------------------------
  # Splines all distributions for a smooth geometry
  _spl_chord = Dierckx.Spline1D(chords[:, 1], chords[:, 2];
                      k= size(chords,1)>=spline_k ? spline_k : size(chords,1)-1,
                                s=spline_s, bc=spline_bc)
  _spl_twist = Dierckx.Spline1D(twists[:, 1], twists[:, 2];
                      k= size(twists,1)>=spline_k ? spline_k : size(twists,1)-1,
                                s=spline_s, bc=spline_bc)
  _spl_LE_x = Dierckx.Spline1D(LE_x[:, 1], LE_x[:, 2];
                      k= size(LE_x,1)>=spline_k ? spline_k : size(LE_x,1)-1,
                                s=spline_s, bc=spline_bc)
  _spl_LE_z = Dierckx.Spline1D(LE_z[:, 1], LE_z[:, 2];
                      k= size(LE_z,1)>=spline_k ? spline_k : size(LE_z,1)-1,
                                s=spline_s, bc=spline_bc)
  if tilt_z!=nothing
    _spl_tlt_z = Dierckx.Spline1D(tilt_z[:, 1], tilt_z[:, 2];
                      k= size(tilt_z,1)>=spline_k ? spline_k : size(tilt_z,1)-1,
                                s=spline_s, bc=spline_bc)
  end

  # ----------------- SPLINE VERIFICATION --------------------------------------
  if verify_spline
    nnodesspan = gt.get_ndivsnodes(grid)[2]    # Number of nodes along span
    y_poss = [gt.get_node(grid, [1,i])[2] for i in 1:nnodesspan]  # Span positions

    fig = plt.figure("spl_verif", figsize=(7*2,5*1))

    plt.subplot(121)
    plt.plot(LE_x[:,1], LE_x[:,2], "og", label="Org LE x", alpha=0.5)
    plt.plot(LE_z[:,1], LE_z[:,2], "ob", label="Org LE z", alpha=0.5)
    plt.plot(y_poss, [_spl_LE_x(y) for y in y_poss], "--g", label="Spline LE x")
    plt.plot(y_poss, [_spl_LE_z(y) for y in y_poss], "--b", label="Spline LE z")
    plt.xlabel(plt.L"y/b_{scale}")
    plt.ylabel(plt.L"x/b_{scale}, z/b_{scale}")
    plt.grid(true, color="0.8", linestyle="--")
    plt.legend(loc="best")

    plt.subplot(122)
    p1 = plt.plot(twists[:,1], twists[:,2], "og", label="Org Twist", alpha=0.5)
    p2 = plt.plot(y_poss, [_spl_twist(y) for y in y_poss], "--g",
                                                          label="Spline twist")
    pextra = []
    if tilt_z!=nothing
      pextra1 = plt.plot(tilt_z[:,1], tilt_z[:,2], "or", label="Org tilt z",
                                                                      alpha=0.5)
      pextra2 = plt.plot(y_poss, [_spl_tlt_z(y) for y in y_poss], "--r",
                                                          label="Spline tilt z")
      pextra = vcat(pextra, [pextra1[1], pextra2[1]])
    end
    plt.ylabel("Twist (deg)")

    plt.grid(true, color="0.8", linestyle="--")
    plt.xlabel(plt.L"y/b_{scale}")

    plt.twinx()
    p3 = plt.plot(chords[:,1], chords[:,2], "ob", label="Org Chord", alpha=0.5)
    p4 = plt.plot(y_poss, [_spl_chord(y) for y in y_poss], "--b",
                                                          label="Spline chord")
    plt.ylabel(plt.L"c/b_{scale}")

    ps = vcat([p1[1], p2[1], p3[1], p4[1]], pextra)
    plt.legend(ps, [p[:get_label]() for p in ps], loc="best")
  end

  # ----------------- AIRFOIL PARAMETERIZATION ---------------------------------
  airfoil_funs = []

  for (pos, airfoil_file) in airfoils
    # Reads the original airfoil geometry from airfoiltools.com
    org_x, org_y = gt.readcontour(airfoil_file;
                                    header_len=rfl_header_len, delim=rfl_delim)

    # Separate upper and lower sides to make the contour injective in x
    upper, lower = gt.splitcontour(org_x, org_y)

    # Parameterize both sides independently
    fun_upper = gt.parameterize(upper[1], upper[2], zeros(upper[1]); inj_var=1,
                                                                      s=rfl_s)
    fun_lower = gt.parameterize(lower[1], lower[2], zeros(lower[1]); inj_var=1,
                                                                      s=rfl_s)

    push!(airfoil_funs, [pos, (fun_upper, fun_lower)])
  end

  # ----------------- AIRFOIL VERIFICATION -------------------------------------
  if verify_rfl

    # Arc-length positions along airfoils
    arcs = [ gt.get_node(grid, [i, 1])[1] for i in 1:gt.get_ndivsnodes(grid)[1]]

    clrs = "bgrycmk"
    clrs = clrs^Int(ceil((size(airfoils,1)/length(clrs))))

    for (i,(pos, airfoil_file)) in enumerate(airfoils)
      # Reads the original airfoil geometry from airfoiltools.com
      org_x, org_y = gt.readcontour(airfoil_file;
                                    header_len=rfl_header_len, delim=rfl_delim)

      gt.plot_airfoil(org_x, org_y; label="pos=$pos", style="."*clrs[i],
                                    zoom_factor=zoom_factor,
                                    side_legend=side_legend, alpha=0.25)

      # Parameterized airfoil
      points_up = [ airfoil_funs[i][2][1](arc) for arc in arcs]
      points_low = [ airfoil_funs[i][2][2](arc) for arc in arcs]
      points = vcat(reverse(points_up), points_low)
      x = [p[1] for p in points]
      y = [p[2] for p in points]

      gt.plot_airfoil(x, y; style="--"*clrs[i], zoom_factor=zoom_factor,
                            side_legend=side_legend)
    end

  end

  # ----------------- SURFACE GRID ---------------------------------------------
  # Auxiliary function for weighting values across span
  function calc_vals(span, array)

    # Finds bounding airfoil position
    val_in, val_out = nothing, array[1]
    for val in array[2:end]
        val_in = val_out
        val_out = val
        if val[1]>=abs(span); break; end
    end
    pos_in = val_in[1]
    val_in = val_in[2]
    pos_out = val_out[1]
    val_out = val_out[2]

    weight = (abs(span)-pos_in)/(pos_out-pos_in)

    return weight, val_in, val_out
  end

  # Space transformation function
  function my_space_transform(X)
    span = X[2]                     # y/bscale span position
    chord = _spl_chord(abs(span))   # c/bscale chord length
    twist = _spl_twist(abs(span))   # twist
    le_x = _spl_LE_x(abs(span))     # x/bscale LE position
    le_z = _spl_LE_z(abs(span))     # z/bscale LE position

    # Merges airfoil contours at this span position
    weight, rfl_in, rfl_out = calc_vals(span, airfoil_funs)
    fun_upper_in, fun_lower_in = rfl_in
    fun_upper_out, fun_lower_out = rfl_out

    # Arc-length on upper or lower side of airfoil
    if X[1]<0.5
        s = 1-2*X[1]
        fun_in = fun_lower_in # Goes over lower side first
        fun_out = fun_lower_out
    else
        s = 2*(X[1]-0.5)
        fun_in = fun_upper_in
        fun_out = fun_upper_out
    end

    # Point over airfoil contour
    point = weight*fun_out(s)+(1-weight)*fun_in(s)

    # Scales the airfoil contour by the normalized chord length
    point = chord*point

    # Applies twist to the airfoil point
    tlt_z = tilt_z!=nothing ?  _spl_tlt_z(abs(span)) : 0.0
    point = gt.rotation_matrix(-twist, -tlt_z, 0)*point

    # Places the point relative to LE and scales by span scale
    point = [point[1]+le_x, span+point[3], point[2]+le_z]*bscale


    return point
  end

  # Transforms the quasi-two dimensional parametric grid into the wing surface
  gt.transform!(grid, my_space_transform)

  # Splits the quadrialateral panels into triangles
  dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  if save_path!=nothing
    # Outputs a vtk file
    gt.save(triang_grid, file_name; path=save_path)

    if paraview
      # Calls paraview
      run(`paraview --data=$save_path$file_name.vtk`)
    end
  end

end

function example_wing(; N::Int64=2, save_path=joinpath(module_path,"../temps/"),
                        paraview=true)
  file_name = "paneledwing"

  # ----------------- READS AND PARAMETERIZES AIRFOILS -------------------------
  semispan = 10            # (m) semi-span length

  chords = [(0, 2.5),      # (semi-span position, chord length (m))
            (0.25, 2.0),
            (1, 1.0)]

  x_pos = [(0, 0),         # (semi-span position, leading edge x-position (m))
           (0.25, semispan/40 ),
           (1, semispan/8 )]

  z_pos = [(0, 0),         # (semi-span position, leading edge z-position (m))
           (0.25, semispan/100 ),
           (1, semispan/50 )]

  twist = [(0, 5),         # (semi-span position, twist (deg))
           (1, 0)]

  airfoils = [(0, "naca6412.dat"), # (semi-span position, airfoil geometry)
              (1,"naca6412.dat")]

  airfoil_funs = []

  for (pos, airfoil_file) in airfoils
      # Reads the original airfoil geometry from airfoiltools.com
      org_x, org_y = gt.readcontour(joinpath(data_path, airfoil_file); header_len=1)

      # Separate upper and lower sides to make the contour injective in x
      upper, lower = gt.splitcontour(org_x, org_y)

      # Parameterize both sides independently
      fun_upper = gt.parameterize(upper[1], upper[2], zeros(upper[1]); inj_var=1)
      fun_lower = gt.parameterize(lower[1], lower[2], zeros(lower[1]); inj_var=1)

      push!(airfoil_funs, [pos, (fun_upper, fun_lower)])
  end


  # ----------------- CREATES PANEL GRID ----------------------------------------
  P_min = [0, -1, 0]            # Lower boundaries arclength, span, dummy
  P_max = [1, 1, 0 ]            # Upper boundaries arclength, span, dummy
  NDIVS = [20, 60, 0]*N         # 20 arclength cells, 60 span cells, 0 dummy
  loop_dim = 1                  # Loop the arclength dimension

  grid = gt.Grid(P_min, P_max, NDIVS, loop_dim)

  gt.plot(grid; labelnodes=!true, labelcells=!true, labelndivs=true,
                  fontsize=8, fig_name="org", title_str="Original Grid")


  # Auxiliary function for weighting values across span
  function calc_vals(span, array)

      # Finds bounding airfoil position
      val_in, val_out = nothing, array[1]
      for val in array[2:end]
          val_in = val_out
          val_out = val
          if val[1]>=abs(span); break; end
      end
      pos_in = val_in[1]
      val_in = val_in[2]
      pos_out = val_out[1]
      val_out = val_out[2]

      weight = (abs(span)-pos_in)/(pos_out-pos_in)

      return weight, val_in, val_out
  end

  # Creates a space transformation function
  function my_space_transform(X)
      span = X[2]

      # Calculates chord
      weight, chord_in, chord_out = calc_vals(span, chords)
      chord = weight*chord_out+(1-weight)*chord_in

      # Calculates airfoil geometry
      weight, rfl_in, rfl_out = calc_vals(span, airfoil_funs)
      fun_upper_in, fun_lower_in = rfl_in
      fun_upper_out, fun_lower_out = rfl_out

      # Arc-length on upper or lower side of airfoil
      if X[1]<0.5
          s = 1-2*X[1]
          fun_in = fun_lower_in # Goes over lower side first
          fun_out = fun_lower_out
      else
          s = 2*(X[1]-0.5)
          fun_in = fun_upper_in
          fun_out = fun_upper_out
      end

      # Point over airfoil contour
      point =  weight*fun_out(s)+(1-weight)*fun_in(s)
      point = chord*point

      # Twist
      weight, twist_in, twist_out = calc_vals(span, twist)
      this_twist = weight*twist_out+(1-weight)*twist_in

      # Applies twist to the airfoil point
      point = gt.rotation_matrix(-this_twist, 0, 0)*point

      # Leading edge x-position
      weight, x_in, x_out = calc_vals(span, x_pos)
      le_x = weight*x_out+(1-weight)*x_in

      # Leading edge z-position
      weight, z_in, z_out = calc_vals(span, z_pos)
      le_z = weight*z_out+(1-weight)*z_in

      # Span position
      y = X[2]*semispan

      return [point[1]+le_x, y, point[2]+le_z]
  end

  # Transforms the quasi-two dimensional grid into the wing surface
  gt.transform!(grid, my_space_transform)


  lims = [-semispan, semispan]
  gt.plot(grid; labelnodes=!true, labelcells=!true, labelndivs=true,
                  fontsize=8, title_str="Transformed grid",
                  xlims=lims/2.*[0,1], ylims=lims, zlims=lims/10);

  # # Adds some dummy example fields
  # gt.add_field(grid, "node_index", "scalar", [i for i in 1:grid.nnodes], "node")
  # gt.add_field(grid, "cell_index", "scalar", [i for i in 1:grid.ncells], "cell")

  # Splits the quadrialateral panels into triangles
  dimsplit = 2              # Dimension along which to split
  triang_grid = gt.GridTriangleSurface(grid, dimsplit)

  # Adds some dummy example fields
  gt.add_field(triang_grid, "node_index", "scalar",
                      [i for i in 1:triang_grid.nnodes], "node")
  gt.add_field(triang_grid, "cell_index", "scalar",
                      [i for i in 1:triang_grid.ncells], "cell")
  gt.add_field(triang_grid, "normal", "vector",
                      [gt.get_normal(triang_grid, i)
                         for i in 1:triang_grid.ncells], "cell")
  gt.add_field(triang_grid, "tangent", "vector",
                      [gt.get_tangent(triang_grid, i)
                         for i in 1:triang_grid.ncells], "cell")

  if save_path!=nothing
      # Outputs a vtk file
      gt.save(triang_grid, file_name; path=save_path)

      if paraview
        # Calls paraview
        run(`paraview --data=$save_path$file_name.vtk`)
      end
  end

  return triang_grid
end
