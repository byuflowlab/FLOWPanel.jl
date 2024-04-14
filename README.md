<img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/julianlogo-flowpanel06.png" alt="FLOWpanel logo" style="width:100%">


<p align="right">
  <span style="color:#7aa098;">
    <i>Three-dimensional panel method for low-speed aerodynamics</i>
  </span>
</p>

<p align="right">
  <a href="https://github.com/byuflowlab/FLOWPanel.jl">
    <img src="https://img.shields.io/badge/code-open%20source-brightgreen.svg">
  </a>

  <a href="https://flow.byu.edu/FLOWPanel.jl">
    <img src="https://img.shields.io/badge/docs-stable-blue.svg">
  </a>
</p>

<p><br></p>

### Capabilities

  > **Structured mesh generation:**
  > *Lofts*
  > *• Bodies of revolution*
  > *• User-defined space transformations*
  >
  > **Unstructured mesh import:**
  > *[OpenVSP](https://openvsp.org)*
  > *• CAD (e.g., SolidWorks) + [Gmsh](https://gmsh.info)*
  > *• [Meshes.jl](https://juliageometry.github.io/MeshesDocs)*
  >
  > **Multi Element:**
  > *Source-panel or doublet/vortex-ring for non-lifting bodies*
  > *• Doublet/vortex-ring for lifting bodies with rigid wake model*
  >
  > **Flexible Solver:**
  > *Direct LU decomposition*
  > *• Iterative Krylov (GMRES)*
  > *• Handling of both watertight and open meshes*
  > *• Multiple bodies*
  >
  > **Fast Computation:**
  > *GPU enabled*
  > *• CPU threaded*
  > *• Low memory allocation*
  > *• Type stable*
  > *• Adjoint solver through [ImplicitAD](https://github.com/byuflowlab/ImplicitAD.jl) for gradient-based optimization*
  > *• Support for both [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [ReverseDiff](https://github.com/JuliaDiff/ReverseDiff.jl)*
  >
  > **Current Limitations:**
  > *Inviscid (no viscous drag nor stall)*
  > *• Incompressible*
  > *• No fast multipole method (though a 30,000 panel problem can still be solved in as fast as 90 seconds on a laptop)*
  >
  > *Coded in [the Julia language](https://www.infoworld.com/article/3284380/what-is-julia-a-fresh-approach-to-numerical-computing.html) for Linux, MacOS, and Windows WSL.*


Tested in Julia v1.10

<p><br></p>

<p align="left">
  <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb009.jpg" alt="img" style="width:100%">
</p>

<p><br></p>

* *Documentation:* [flow.byu.edu/FLOWPanel.jl](https://flow.byu.edu/FLOWPanel.jl)
* *Code:* [github.com/byuflowlab/FLOWPanel.jl](https://github.com/byuflowlab/FLOWPanel.jl)

<p><br></p>

### Installation Instructions
1. Download and install [Julia](https://julialang.org/)
2. Type in the Julia REPL: `] add FLOWPanel`
3. *(optional)* For visualization, install [ParaView](https://www.paraview.org/) and make sure that [it is callable from the terminal](https://flow.byu.edu/FLOWUnsteady/installation/general/#paraview) typing `paraview`

Copy and paste any of [the examples](https://flow.byu.edu/FLOWPanel.jl/stable/examples/sweptwing-4p2aoa/) directly in the Julia REPL to run your first simulation

<p><br></p>

### Selected Publications
See the following publication for an in-depth dive into the theory and validation:

* E. J. Alvarez, C. Joseph, & A. Ning (2023), "Vortex Particle Method for Electric Ducted Fan in Non-Axisymmetric Flow," *AIAA AVIATION Forum*. [**[VIDEO]**](https://www.youtube.com/watch?v=BQpar3A0X-w&hd=1) [**[SLIDES]**](http://edoalvar2.groups.et.byu.net/public/FLOWUnsteady/alvarez_2023-SLIDES-VPM_for_EDF_in_Non_Axisymmetric_Flow.pdf) [**[PDF]**](https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=7676&context=facpub)

<p><br></p>

### About
FLOWPanel is an open-source project jointly led by [Whisper Aero](http://whisper.aero/) and the [FLOW Lab](http://flow.byu.edu/) at Brigham Young University.
All contributions are welcome.

  * Developers/contributors : [Eduardo J. Alvarez](https://www.edoalvarez.com/) (main) and [Cibin Joseph](https://github.com/cibinjoseph)
  * Email           : edo.alvarezr@gmail.com
  * Created         : August 2021
  * License         : MIT License

<p><br></p>

### Sponsors

<p align="right">
  <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/sponsors01.png" alt="img" style="width:90%">
  <br><br><br>
</p>

<p><br></p>


<p align="center">
  <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/cessna002.jpg" alt="img" style="width:75%">
  <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/zeroebwb000.jpg" alt="img" style="width:75%">
</p>

<p><br></p>

<p align="center">
  <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/duct-hill-aoa15-slice02-small.png" alt="img" style="width:45%"><img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/light/sphere01_2.gif" alt="img" style="width:45%">
</p>
