```@raw html
<center>
  <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna000.jpg" alt="Pic here" style="width: 100%;"/>
</center>
```

# Import Mesh and Solve
Here we import the mesh into FLOWPanel using
[Meshes.jl](https://juliageometry.github.io/MeshesDocs), identify the
trailing edge, and run the watertight solver.


!!! info "Default OpenVSP files"
    We have pre-generated and uploaded an OpenVSP mesh to this example so that
    you can run this section without needing to complete the previous
    sections. However, if you would like to use your own mesh, simply change
    `read_path`, `meshfile`, and `trailingedges` to point to
    your files.

```julia
#=##############################################################################
# DESCRIPTION
    Cessna 210 aircraft. The mesh used in this analysis was created using
    OpenVSP.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Apri 2024
  * License   : MIT License
=###############################################################################

import FLOWPanel as pnl
import FLOWPanel: norm, dot, cross

import Meshes
import GeoIO
import Rotations: RotX, RotY, RotZ

# import CUDA                               # Uncomment this to use GPU (if available)


run_name        = "cessna"                  # Name of this run

save_path       = run_name                  # Where to save outputs
paraview        = true                      # Whether to visualize with Paraview
read_path       = joinpath(pnl.examples_path, "data") # Where to read Gmsh files from

# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 4.0                       # (deg) freestream angle of attack
magVinf         = 180 * 0.514444            # (m/s) freestream velocity
rho             = 0.71461                   # (kg/m^3) air density at 17300 ft


# ----------------- GEOMETRY DESCRIPTION ---------------------------------------
meshfile        = joinpath(read_path, "cessna.msh")    # Gmsh file to read

offset          = [0, 0, 0]                 # Offset to center the mesh
rotation        = RotX(0*pi/180) * RotY(0*pi/180) * RotX(0*pi/180) # Rotation to align mesh
scaling         = 0.3048                    # Factor to scale mesh dimensions (ft -> m)

                                            # Gmsh files with trailing edges
trailingedges   = [ #   ( Gmsh file,                 span direction )
                        ("cessna-TE-leftwing.msh",      [0, 1, 0]),
                        ("cessna-TE-rightwing.msh",     [0, 1, 0]),
                        ("cessna-TE-leftelevator.msh",  [0, 1, 0]),
                        ("cessna-TE-rightelevator.msh", [0, 1, 0]),
                        ("cessna-TE-rudder.msh",        [0, 0, 1]),
                    ]

flip            = true                      # Whether to flip control points against the direction of normals
                                            # NOTE: use `flip=true` if the normals
                                            #       point inside the body

bref            = 11.2014                   # (m) reference span
cref            = 1.449324                  # (m) reference chord
Sref            = 16.2344                   # (m^2) reference area
Xac             = [2.815, 0, 0.4142]        # (m) aerodynamic center for
                                            # calculation of moments

# ----------------- SOLVER SETTINGS -------------------------------------------

# Solver: direct linear solver for open bodies
# bodytype = pnl.RigidWakeBody{pnl.VortexRing} # Wake model and element type

# Solver: least-squares solver for watertight bodies
bodytype        = pnl.RigidWakeBody{pnl.VortexRing, 2}


# ----------------- GENERATE BODY ----------------------------------------------
# Read Gmsh mesh
msh = GeoIO.load(meshfile)
msh = msh.geometry

# Transform the original mesh: Translate, rotate, and scale
msh = msh |> Meshes.Translate(offset...) |> Meshes.Rotate(rotation) |> Meshes.Scale(scaling)

# Wrap Meshes object into a Grid object from GeometricTools
grid = pnl.gt.GridTriangleSurface(msh)

# Read all trailing edges
sheddings = []

for (trailingedgefile, spandir) in trailingedges

    # Read Gmsh line of trailing edge
    TEmsh = GeoIO.load(joinpath(read_path, trailingedgefile))
    TEmsh = TEmsh.geometry

    # Apply the same transformations of the mesh to the trailing edge
    TEmsh = TEmsh |> Meshes.Translate(offset...) |> Meshes.Rotate(rotation) |> Meshes.Scale(scaling)

    # Convert TE Meshes object into a matrix of points used to identify the trailing edge
    trailingedge = pnl.gt.vertices2nodes(TEmsh.vertices)

    # Sort TE points from "left" to "right" according to span direction
    trailingedge = sortslices(trailingedge; dims=2, by = X -> pnl.dot(X, spandir))

    # Function for identifying if a point is close to a junction
    if trailingedgefile in ["cessna-TE-leftwing.msh", "cessna-TE-rightwing.msh"]

        criterion = X -> abs(X[2]) <= 0.67

    elseif trailingedgefile in ["cessna-TE-leftelevator.msh", "cessna-TE-rightelevator.msh"]

        criterion = X -> abs(X[2]) <= 0.23

    elseif trailingedgefile in ["cessna-TE-rudder.msh"]

        criterion = X -> X[3] <= 0.65

    else

        criterion = X -> false

    end

    # Filter out any points that are close to junctions
    tokeep = filter( i -> !criterion(trailingedge[:, i]), 1:size(trailingedge, 2) )
    trailingedge = trailingedge[:, tokeep]

    # Generate TE shedding matrix
    shedding = pnl.calc_shedding(grid, trailingedge; tolerance=0.001*bref)

    push!(sheddings, shedding)

end

# Combine all TE shedding matrices into one
shedding = hcat(sheddings...)

# Generate paneled body
body = bodytype(grid, shedding; CPoffset=(-1)^flip * 1e-14)

println("Number of panels:\t$(body.ncells)")


# ----------------- CALL SOLVER ------------------------------------------------
println("Solving body...")

# Freestream vector
Vinf = magVinf*[cos(AOA*pi/180), 0, sin(AOA*pi/180)]

# Freestream at every control point
Uinfs = repeat(Vinf, 1, body.ncells)

# Unitary direction of semi-infinite vortex at points `a` and `b` of each
# trailing edge panel
Das = repeat(Vinf/magVinf, 1, body.nsheddings)
Dbs = repeat(Vinf/magVinf, 1, body.nsheddings)

# Solve body (panel strengths) giving `Uinfs` as boundary conditions and
# `Das` and `Dbs` as trailing edge rigid wake direction
@time pnl.solve(body, Uinfs, Das, Dbs)

# Uncomment this to use GPU instead (if available)
# @time pnl.solve(body, Uinfs, Das, Dbs; GPUArray=CUDA.CuArray{Float32})

# ----------------- POST PROCESSING ----------------------------------------
println("Post processing...")

# Calculate surface velocity U on the body
Us = pnl.calcfield_U(body, body)

# NOTE: Since the boundary integral equation of the potential flow has a
#       discontinuity at the boundary, we need to add the gradient of the
#       doublet strength to get an accurate surface velocity

# Calculate surface velocity U_∇μ due to the gradient of the doublet strength
UDeltaGamma = pnl.calcfield_Ugradmu(body)
# UDeltaGamma = pnl.calcfield_Ugradmu(body; sharpTE=true, force_cellTE=false)

# Save this intermediate result in a separate field for debugging
pnl.add_field(body, "Uinfind", "vector", collect.(eachcol(Us)), "cell")

# Add both velocities together
pnl.addfields(body, "Ugradmu", "U")

# Calculate pressure coefficient (based on U + U_∇μ)
@time Cps = pnl.calcfield_Cp(body, magVinf)

# Calculate the force of each panel (based on Cp)
@time Fs = pnl.calcfield_F(body, magVinf, rho)

# Calculate total force of the vehicle decomposed as lift, drag, and sideslip
Dhat = Vinf/norm(Vinf)                      # Drag direction
Shat = [0, 1, 0]                            # Span direction
Lhat = cross(Dhat, Shat)                    # Lift direction

LDS = pnl.calcfield_LDS(body, Lhat, Dhat)
L, D, S = collect(eachcol(LDS))

# Force coefficients
nondim = 0.5*rho*magVinf^2*Sref             # Normalization factor
CL = sign(dot(L, Lhat)) * norm(L) / nondim
CD = sign(dot(D, Dhat)) * norm(D) / nondim

# Integrated moment decomposed into rolling, pitching, and yawing moments
lhat = Dhat                                 # Rolling direction
mhat = Shat                                 # Pitching direction
nhat = Lhat                                 # Yawing direction

lmn = pnl.calcfield_lmn(body, Xac, lhat, mhat, nhat)
roll, pitch, yaw = collect(eachcol(lmn))

# Moment coefficients
nondim = 0.5*rho*magVinf^2*Sref*cref # Normalization factor

Cl = sign(dot(roll, lhat)) * norm(roll) / nondim
Cm = sign(dot(pitch, mhat)) * norm(pitch) / nondim
Cn = sign(dot(yaw, nhat)) * norm(yaw) / nondim

@show L
@show D
@show CL
@show CD
@show Cl
@show Cm
@show Cn

# ----------------- VISUALIZATION ------------------------------------------

# Save body as VTK
vtks = save_path*"/"                        # String with VTK output files
vtks *= pnl.save(body, run_name; path=save_path, debug=false)

# NOTE: using `debug=true` outputs the control points and normal, but it takes
#       much longer since it also calculates the velocity at those points

# Call Paraview
if paraview
    run(`paraview --data=$(vtks)`)
end
```
```@raw html
<span style="font-size: 0.9em; color:gray;"><i>
    Number of panels: 24,000. <br>
    Run time: ~100 seconds on a Dell Precision 7760 laptop (no GPU). <br>
</i></span>
<br><br>
```

|                         | VSPAERO | FLOWPanel |
| ----------------------: | :-----: | :-------: |
| Lift ``C_L``            | 0.499   | 0.554     |
| Drag ``C_D``            | 0.0397  | 0.0496    |


```@raw html
<center>
  <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna002.jpg" alt="Pic here" style="width: 100%;"/>
  <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna005.jpg" alt="Pic here" style="width: 100%;"/>
</center>
<br><br>
<br><br>
```



!!! details "Tip"
    You can also automatically run this example
    with the following command:
    ```julia
    import FLOWPanel as pnl

    include(joinpath(pnl.examples_path, "cessna.jl"))
    ```


!!! details "Checklist for importing meshes"
    1. **Check whether normals point into the body:** Using the flag
        `debug=true` in `pnl.save(body, run_name; path=save_path, debug=true)`
        will output the control points of the body along with the associated
        normal vector of each panel.
            We recommend opening the body and control points in ParaView and
        visualizing the normals with the Glyph filter.
            Whenever the normals are pointing into the body, the user needs
        to flip the offset of the
        control points with `CPoffset=-1e-14` or any other negligibly small
        negative number. This won't flip the normals outwards, but it will flip
        the zero-potential domain from outwards back to inside the body
        (achieved by shifting the control points slightly into the body).
        If you pull up the solution in ParaView and realize that the surface
        velocity is much smaller than the freestream everywhere along the
        aircraft, that's an indication that the normals are point inwards
        and you need to set `CPoffset` to be negative.
    2. **Check that the trailing edge was correctly identified:**
        `pnl.save(body, run_name; path=save_path)` automatically outputes the
        wake.
            We recommend opening the body and wake in ParaView and visually
        inspecting that the wake runs along the trailing edge line that you
        defined under `trailingedge`.
            If not successful, increase the resolution of `trailingedge` and tighten
        the tolerance to something small like
        `pnl.calc_shedding(grid, trailingedge; tolerance=0.0001*span)`.
    3. **Choose the right solver for the geometry:**
        Use the least-squares solver with watertight bodies
        (`bodytype = pnl.RigidWakeBody{pnl.VortexRing, 2}`), and the direct
        linear solver with open bodies
        (`bodytype = pnl.RigidWakeBody{pnl.VortexRing}`). The least-squares
        solver runs much faster in GPU
        (`pnl.solve(body, Uinfs, Das, Dbs; GPUArray=CUDA.CuArray{Float32})`),
        but it comes at the price of sacrificing accuracy (single precision
        numbers as opposed to double).

!!! tip "Visualization"
    To help you practice in ParaView, we have uploaded the solution files
    of this simulation along with the ParaView state file (`.pvsm`) that
    we used to generate the visualizations shown above:
    [DOWNLOAD](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna.zip)

    To open in ParaView: `File → Load State → cessna.pvsm` then
    select "Search files under specified directory" and point it to the
    folder with the outputs of FLOWPanel.

