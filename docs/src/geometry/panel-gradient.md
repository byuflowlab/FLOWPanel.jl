# Panel Gradients

Computing the gradient of panel circulation is necessary when using Vortex ring elements to model a geometry. Unlike in structured meshes, gradient estimation in unstructured meshes is not straightforward due to the lack of consistent connectivity. Multiple methods have been proposed to compute gradients on unstructured meshes that utilize either the cell-centered or node-based circulation strengths [\[1\]](https://ntrs.nasa.gov/api/citations/20140011550/downloads/20140011550.pdf). A few of these are listed below:

1. Green-Gauss gradient method
    - Simple face averaging
    - Inverse distance weighted face interpolation
    - Weighted Least squares Face interpolation
    - Weighted Tri-Linear Face Interpolation
2. Least squares gradient method
3. Curvilinear gradient method 

Since cell-centered approaches have inherent difficulties handling boundaries and 'holes' in the geometry, a node-based approach is utlized here. This requires converting cell-centered panel strengths to node-based. An area-weighted averaging procedure enables this conversion. Once node-based strengths are obtained, the gradient is computed as the inclination of the plane passing through the scalar values at each cell's vertex. Steps involved in this procedure are described below.

Consider a single triangular panel element with the scalar values $\phi_1$, $\phi_2$, $\phi_3$ at its three vertices $P_1$, $P_2$, $P_3$. First, the three vertices are converted to a local two-dimensional frame.

```@raw html
<center>
  <img src="../../assets/images/geometry-panel-gradient.png" alt="panel gradient" style="width: 380px;"/>
</center>
<br><br>
```

A plane passing through $\phi_1$, $\phi_2$, $\phi_3$ satisfies the equations,

$\alpha + \beta (x_1 - x_0) + \gamma (y_1 - y_0) = \phi_1\\
 \alpha + \beta (x_2 - x_0) + \gamma (y_2 - y_0) = \phi_2\\
 \alpha + \beta (x_3 - x_0) + \gamma (y_3 - y_0) = \phi_3$

where $(x_0, y_0)$ is the cell center, and $\alpha$ the value of $\phi$ at the cell-center. $\beta$ and $\gamma$ are the gradients (or slope) in the $x$ and $y$ directions in the two-dimensional frame.

The unknowns $\alpha$, $\beta$, and $\gamma$ may be obtained by solving the system of linear equations given by,

$\begin{bmatrix}
1 & (x_1-x_0) & (y_1-y_0) \\
1 & (x_1-x_0) & (y_1-y_0) \\
1 & (x_1-x_0) & (y_1-y_0)
\end{bmatrix} \begin{bmatrix}
\alpha \\
\beta \\
\gamma
\end{bmatrix} = \begin{bmatrix}
\phi_1 \\
\phi_2 \\
\phi_3
\end{bmatrix}$

The slopes $\beta$ and $\gamma$ are converted back to three-dimensions using the appropriate basis vectors as shown below,

$\mathbf{\nabla} \phi = \beta \mathbf{e_x} + \gamma \mathbf{e_y}$

### References
[1] Sozer, E., Brehm, C., & Kiris, C. C. (2014). Gradient calculation methods on arbitrary polyhedral unstructured meshes for cell-centered CFD solvers. 52nd Aerospace Sciences Meeting, 1–24.
