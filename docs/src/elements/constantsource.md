# Constant-Strength Source

  (Adapted fron Hess, J. L., & Smith, A. M. O. (1967), *Calculation of potential flow about arbitrary bodies*)

  Given the following **planar** polygonal panel

```@raw html
<center>
  <img src="../../assets/images/panel-HS00.png" alt="Pic here" style="width: 400px;"/>
</center>
```

the potential at an arbitrary point $(x,\, y,\, z)$ is calculated in the panel's coordinate system with constant-strength source as

```math
  \begin{align*}
        \phi
    =
        -\frac{\sigma}{4\pi} \int\limits_{S} \frac{\mathrm{d}S}{r}
  ,\end{align*}
```
with $r=\sqrt{(x-\xi)^2 + (y-\eta)^2 + z^2}$.

> **NOTE:** The local coordinate system used by Hess and Smith (and also Katz and Plotkin) follows the opposite to the right-hand rule to define its normal. For this reason, the local coordinate system implemented in all $\phi$ and $\mathbf{u}$ functions in FLOWPanel define $\hat{\boldsymbol\xi} = \hat{\mathbf{o}}$, $\hat{\boldsymbol\eta} = \hat{\mathbf{t}}$, and $\hat{\mathbf{n}}_\mathrm{HS} = -\hat{\mathbf{n}}$.

Integrating over a planar element with $n$ vertices results in

```math
\begin{align*}
        \phi(x,y,z)
    =
        -\frac{\sigma}{4\pi}
        \sum\limits_{i,j\in A} \left[
            R_{ij} Q_{ij}
            + %
            \vert z \vert J_{ij}
        \right]
,\end{align*}
```

where $A = \{(1,2),\,\dots,\,(n-1, n),\,(n, 1) \}$ and

```math
\begin{align*}
        \bullet \quad & S_{i,j} = \frac{\eta_j - \eta_i}{d_{i,j}}
    \\
        \bullet \quad & C_{i,j} = \frac{\xi_j - \xi_i}{d_{i,j}}
    \\
        \bullet \quad & Q_{i,j} = \ln{\left(\frac{ r_i+r_j+d_{i,j} }{ r_i+r_j-d_{i,j} }\right)}
    \\
        \bullet \quad & J_{i,j} = \arctan{\left(\frac{
                                R_{i,j}\lvert z \rvert( r_i s_{i,j}^{(j)} - r_j s_{i,j}^{(i)})
                            }{
                                r_i r_j R_{i,j}^2 + z^2 s_{i,j}^{(j)} s_{i,j}^{(i)}
                            }\right)}
    \\
        \bullet \quad & s_{i,j}^{(k)} = (\xi_k - x)C_{i,j} + (\eta_k - y)S_{i,j}
    \\
        \bullet \quad & R_{i,j} = (x - \xi_i)S_{i,j} - (y - \eta_i)C_{i,j}
    \\
        \bullet \quad & d_{i,j} = \sqrt{(\xi_j-\xi_i)^2 + (\eta_j-\eta_i)^2}
    \\
        \bullet \quad & r_i = \sqrt{(x-\xi_i)^2 + (y-\eta_i)^2 + z^2}
.\end{align*}
```

> **NOTE:** The $\arctan$ defined in $J_{ij}$ is intended to be evaluated in the range $-\pi$ to $\pi$.

The velocity induced by the panel at $(x,y,z)$ is calculated as

```math
\begin{align*}
    \mathbf{u} = \nabla \phi
\end{align*}
```
which results in
```math
\begin{align*}
        \bullet \quad & U_x = \frac{\partial\varphi}{\partial x} = - \frac{\sigma}{4\pi} \sum\limits_{i=1}^n S_i Q_i
    \\
        \bullet \quad & U_y = \frac{\partial\varphi}{\partial y} = \frac{\sigma}{4\pi} \sum\limits_{i=1}^n C_i Q_i
    \\
        \bullet \quad & U_z = \frac{\partial\varphi}{\partial z} = \text{sgn}(z) \frac{\sigma}{4\pi} \left(  \Delta\theta - \sum\limits_{i=1}^n J_i \right)
,\end{align*}
```

where $\Delta\theta=2\pi$ if $(x,y,0)$ lies inside the quadrilateral, $\Delta\theta=0$ if not. Hess & Smith mentions that the point lies inside the quadrilateral iff all $R_{i}$ are positive.

> **NOTE:** The $U_z$ velocity poses a discontinuity at the surface since $\lim\limits_{z\rightarrow \pm 0} U_z (0, 0, z) = \pm \frac{\sigma}{2}$. Hess and Smith recommends setting $U_z (0, 0, 0) = + \frac{\sigma}{2}$. In FLOWPanel, however, we let $U_z (0, 0, 0) = 0$, but we also have shifted all control points slightly in the direction of $\hat{\mathbf{n}}$. Remembering that $\hat{\mathbf{n}}_\mathrm{HS} = -\hat{\mathbf{n}}$, the control points are thus shifted in the $-z$ direction, effectively obtaining $\boxed{U_z (\mathbf{x}_\mathrm{cp}) \approx - \frac{\sigma}{2}}$.

> **NOTE 2:** The offset of the control points is controlled through the properties `body.CPoffset` and `body.characteristiclength` of the Body type. Here, `body.CPoffset` is a small non-dimensional number $s$ and `body.characteristiclength` is a user-defined function of the form `(nodes, panel) -> l` that returns a characteristic length $\ell$ (which can be either computed based on the panel, or it can be set the same for all panels). Each control point $\mathbf{x}_\mathrm{cp}$ is then computed as $\mathbf{x}_\mathrm{cp} = \mathbf{x}_\mathrm{centroid} + s\ell\hat{\mathbf{n}}$.

> **NOTE 3:** By default, FLOWPanel sets the characteristic length to be the square root of the panel area, but **it is strongly recommended that the user provides their own characteristic length, and that this length is the same for all panels**.


**ASSUMPTIONS**

* The panel is a polygon of $n$ number of vertices with $n\ge3$.
* The polygon is planar, *i.e.,* all vertices lay on the same plane.
* Vectors betweens nodes 1 and 2 and nodes 1 and 3 are not collinear.
* **The polygon must be concave**.

**OBSERVATIONS**
* The term $Q_{i,j}$ makes this formulation singular at all vertices and edges of a panel; hence, FLOWPanel adds a small epsilon to the denominator of the log argument to avoid the singularity.

> **NOTE:** The small offset added to the denominator of $Q_{i,j}$ corresponds to `body.kerneloffset`.

The potential and velocity field of a source panel of unitary strength ($\sigma=1$) is shown below

```@raw html
<center>
  <table>
      <tr>
          <td>
              <img src="../../assets/images/panel-source-phi00.png" alt="Pic here" width="450px">
          </td>
          <td>
              <img src="../../assets/images/panel-source-u02.png" alt="Pic here" width="300px">
          </td>
      </tr>
  </table>
</center>
```

```@raw html
<center>
  <br>$\nabla \phi$ = $\mathbf{u}$ verification
  <img src="../../assets/images/panel-source-divphivsu00.png" alt="Pic here" style="width: 700px;"/>
</center>
```
