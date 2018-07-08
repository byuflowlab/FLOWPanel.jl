#=##############################################################################
# DESCRIPTION
    Utilities.
# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jul 2018
  * License   : AGPL-3.0
=###############################################################################

function example_wing(; N::Int64=2)
  file_name = "paneledwing"
  paraview = true

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
      org_x, org_y = readcontour(joinpath(data_path, airfoil_file); header_len=1)

      # Separate upper and lower sides to make the contour injective in x
      upper, lower = splitcontour(org_x, org_y)

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
          fun_in = fun_upper_in
          fun_out = fun_upper_out
      else
          s = 2*(X[1]-0.5)
          fun_in = fun_lower_in
          fun_out = fun_lower_out
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
  dimsplit = 1              # Dimension along which to split
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


  if paraview
      # Outputs a vtk file
      gt.save(triang_grid, file_name)

      # Calls paraview
      run(`paraview --data=$file_name.vtk`)

      # Delete vtk file
      run(`rm -f $file_name.vtk`)
  end
end



"Receives a .dat file as pulled from airfoiltools.com containing the x and y
contour coordinates of an airfoil, and returns arrays x and y."
function readcontour(file_name; header_len=1)
  x, y = Float64[], Float64[]

  open(file_name) do f
    for (i,line) in enumerate(eachline(f))

      # Ignores header
      if i<=header_len
        nothing
      # Parses each line
      else
        this_x, this_y = split(line)
        push!(x, parse(Float64, this_x))
        push!(y, parse(Float64, this_y))
      end

    end
  end
  return x,y
end


"""
  Receives an airfoil contour and splits it up in upper and lower surfaces as
  divided by the chord line. It returns `(upper, lower)` with `upper=(x,y)` the
  points of the upper surface, ditto for `lower`. Both `upper` and `lower` are
  given in increasing order in x (i.e., from leading to trailing edge).
"""
function splitcontour(x,y)
  # ERROR CASES
  if !(x[1] in [0.0, 1.0])
    error("Invalid contour. x[1] must be either 0.0 or 1.0, got $(x[1]).")
  end

  # Flag indicating whether the contour start at the trailing or leading edge
  start_TE = x[1]==1.0

  # Find the opposite end of the contour
  end_i = -1
  for (i, xi) in enumerate(x)
    if i==1
      nothing
    # Case of starting from the trailing edge
    elseif start_TE && xi > x[i-1]
      end_i = i-1
      break
    # Case of leading edge
    elseif !start_TE  && xi < x[i-1]
      end_i = i-1
      break
    end
  end

  # ERROR CASE
  if end_i==-1
    error("Logic error! End of airfoil not found!")
  end

  # Splits them up
  x_sec1, y_sec1 = x[1:end_i], y[1:end_i]
  x_sec2, y_sec2 = x[end_i:end], y[end_i:end]

  # Sorts them from LE to TE
  if x_sec1[1] > 0.5; reverse!(x_sec1); reverse!(y_sec1); end;
  if x_sec2[1] > 0.5; reverse!(x_sec2); reverse!(y_sec2); end;

  # Determines upper and lower surfaces
  if mean(y_sec1) > mean(y_sec2)
    upper = [x_sec1, y_sec1]
    lower = [x_sec2, y_sec2]
  else
    upper = [x_sec2, y_sec2]
    lower = [x_sec1, y_sec1]
  end

  return upper, lower
end
