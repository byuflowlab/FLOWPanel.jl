![FLOWPanel logo](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/julianlogo-flowpanel06.png)

*Three-dimensional panel method for high-Reynolds aerodynamics*

[![](https://img.shields.io/badge/code-open%20source-brightgreen.svg)](https://github.com/byuflowlab/FLOWPanel.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](http://flow.byu.edu/FLOWPanel.jl/dev/)


### Features
* Structured mesh generation through space transformation, lofting, or body of revolution
* Source-panel or vortex-ring solver for non-lifting bodies
* Vortex-ring solver for lifting bodies with rigid wake
* Low-memory allocation
* Direct and iterative Krylov (GMRES) [solvers](http://flow.byu.edu/FLOWPanel.jl/dev/examples/sweptwing-solver/)
* Automatic differentiation for gradient-based optimization (supports [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl), [ReverseDiff](https://github.com/JuliaDiff/ReverseDiff.jl), and [ImplicitAD](https://github.com/byuflowlab/ImplicitAD.jl))

Developed and tested in Julia v1.6.

### Installation Instructions
1. Download and install [Julia](https://julialang.org/)
2. In the Julia REPL: `] add FLOWPanel`
3. (optional) For visualization, install [Paraview](https://www.paraview.org/) and make sure that [it is callable from the terminal](https://flow.byu.edu/FLOWUnsteady/tutorials/installation-instructions/#Paraview) typing `paraview`

Copy and paste any of [the examples](http://flow.byu.edu/FLOWPanel.jl/dev/examples/sweptwing-4p2aoa/) directly in the REPL to run your first simulation

### Sponsors

![sponsors](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/sponsors01.png)

### About
FLOWPanel is an open-source project jointly led by [Whisper Aero](http://whisper.aero/) and the [FLOW Lab](http://flow.byu.edu/) at Brigham Young University.
All contributions are welcome.

  * Main Developer  : [Eduardo J. Alvarez](https://edoalvarez.com/)
  * Email           : edo.alvarezr@gmail.com
  * Created         : August 2021
  * License         : MIT License


Source code: [LINK](https://github.com/byuflowlab/FLOWPanel.jl)

Documentation: [LINK](http://flow.byu.edu/FLOWPanel.jl/dev/)

![sphere](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/duct-hill-aoa15-slice02-small.png)
![sphere](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/light/sphere01_2.gif)
