# VSPGeom.jl

[VSPGeom.jl](https://flow.byu.edu/VSPGeom.jl) is a Julia package for reading
and processing the geometry of OpenVSP, which we can also use with
FLOWPanel.
For instance, in OpenVSP we can generate a CFD-quality surface mesh that
gets exported as an STL file.
The STL can then be read following [these instructions](https://flow.byu.edu/VSPGeom.jl/dev/howto/#STL-Files) and wrapped as
a [Meshes.jl](https://juliageometry.github.io/MeshesDocs) object that is
passed to FLOWPanel.

```@raw html
<center>
  <img src="https://edoalvar2.groups.et.byu.net/public/FLOWPanel/openvsp004.png" alt="Pic here" style="width: 100%;"/>
</center>
```

