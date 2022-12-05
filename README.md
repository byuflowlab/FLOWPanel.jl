![FLOWPanel logo](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/julianlogo-flowpanel06.png)

*Three-dimensional panel method for high-Reynolds aerodynamics*

[![](https://img.shields.io/badge/code-open%20source-brightgreen.svg)](https://github.com/byuflowlab/FLOWPanel.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](http://flow.byu.edu/FLOWPanel.jl/dev/)


### Features
* Structured mesh generation through space transformation, lofting, or body of revolution.
* Source-panel solver for non-lifting bodies.
* Vortex-ring solver for lifting bodies.
* Low-memory allocation.
* Automatic differentiation for gradient-based optimization: ForwardDiff.

Developed and tested in Julia v1.6.

### Installation Instructions
TODO: register this package in official Julia registry.

Dependencies
  * GeometricTools: `https://github.com/byuflowlab/GeometricTools.jl`

### Sponsors

![sponsors](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/sponsors01.png)

### To-Do's
* Finish documentation
* Swept-wing validation and example
* Implement actuator disk for ducted fan
* Test in gradient-based optimization.

### Future Work
* Currently, all solvers use a direct explicit-matrix inversion which works well for problems with less than 1000 panels. In order the scale the solver, it is recommended that an indirect matrix inversion method (e.g., GMRES) is implemented in future work along with an implicit matrix evaluation through FMM.

### About
FLOWPanel is an open-source project jointly led by the [FLOW Lab](http://flow.byu.edu/) at Brigham Young University and [Whisper Aero](http://whisper.aero/).
All contributions are welcome.

  * Main Developer  : [Eduardo J. Alvarez](https://edoalvarez.com/)
  * Email           : edo.alvarezr@gmail.com
  * Created         : August 2021
  * License         : MIT License

![sphere](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/light/sphere01_2.gif)
![box](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/light/box01_2.gif)
![hub](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/light/hub03_2.gif)
