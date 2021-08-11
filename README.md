<img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/sphere01_2.gif" alt="Vid here" style="width: 500px;"/>

![](http://edoalvar2.groups.et.byu.net/public/FLOWPanel/sphere01_2.gif)

<img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/sphere01_2.gif" alt="Vid here" style="width: 500px;"/>

![](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/sphere01_2.gif)

<img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/sphere01_2.gif" alt="Databay showcase gif" title="Databay showcase gif" width="500"/>

![](https://edoalvar2.groups.et.byu.net/public/FLOWPanel/sphere01_2.gif)


<img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/sphere01_2.gif?raw=true" alt="Databay showcase gif" title="Databay showcase gif" width="500"/>

# FLOWPanel
Three-dimensional panel method for high-Reynolds aerodynamics.
Developed and tested on Julia v1.6. See [documentation](https://nbviewer.jupyter.org/github/byuflowlab/FLOWPanel.jl/blob/master/docs/Documentation.ipynb) under [`docs/Documentation.ipynb`](https://nbviewer.jupyter.org/github/byuflowlab/FLOWPanel.jl/blob/master/docs/Documentation.ipynb) for theory, implementation, examples, and validation of this code.

### CURRENT CAPABILITIES
* Arbitrary mesh generation through space transformation, lofting, or body of revolution.
* Source-panel solver for non-lifting bodies.

### WORK IN PROGRESS
* Doublet-panel solver for non-lifting bodies has been implemented and verified, but has not been validated yet.
* Vortex-ring-panel solver for lifting bodies with rigid wake has been implemented and verified, but has not been validated yet.

### FUTURE WORK
* Currently, all solvers use a direct explicit-matrix inversion which works well for problems with less than 1000 panels. In order the scale the solver, it is recommended that an indirect matrix inversion method (e.g., GMRES) is implemented in future work along with an implicit matrix evaluation through FMM.

# Dependencies
  * GeometricTools: `https://github.com/byuflowlab/GeometricTools.jl`

# Authorship
  * Main Developer  : Eduardo J. Alvarez
  * Email           : Edo.AlvarezR@gmail.com
  * Created         : August 2021
  * License         : MIT License

<table style="width:100%">
  <tr>
    <td>
      <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/box00_2.gif" alt="Vid here" style="width: 500px;"/>
    </td>
    <td>
      <img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/box01_2.gif" alt="Vid here" style="width: 500px;"/>
    </td>
  </tr>
</table>

<img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/hub03_2.gif" alt="Vid here" style="width: 500px;"/>
<img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/wing05.gif" alt="Vid here" style="width: 500px;"/>
<img src="http://edoalvar2.groups.et.byu.net/public/FLOWPanel/prop07.gif" alt="Vid here" style="width: 500px;"/>
