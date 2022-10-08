# Methods

While the previous section outlines the basics on how to generate geometry manually using GeometricTools, this section shows how geometry can be automatically generated.

Three methods are reviewed:
* Surface lofting
* Surface of revolution
* Space transformation

Each of these methods return a grid of **quadrilateral panels that are not necessarily planar**. Because some element types require planar panels, we also review a method for generating a triangular grid out of any grid obtained from the previous methods.
