# `GeometricTools.jl`

Some of the solvers implemented in FLOWPanel require flat panels and a structured surface grid that facilitates the declaration of edges along which to specify the Kutta condition.
Hence, all paneled bodies are required to be structured grids.

The package [`GeometricTools.jl`](https://github.com/byuflowlab/GeometricTools.jl) implements a variety of methods for the generation of grids and meshes, refinement, and space transformations, and is chosen as the geometric engine for FLOWPanel. Here below we show the low level manipulation of `GeometricTools.jl` for the sake of showing how to manually define complex geometries based on space transformations and localized refinement. The section [Methods](@ref) under Grid Generation shows predefined methods for creating common geometries, *e.g.*, lofted surfaces, surfaces of revolution, etc.
