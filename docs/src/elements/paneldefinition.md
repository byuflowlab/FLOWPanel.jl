# Panel Definition

In order to compute and/or solve the boundary integral equation, we discretize the numerical boundaries of our domain through numerical elements, forming what is called a *boundary element method*. or BEM.
The geometry of our problem is discretized into panels. To each panel then we associate one or multiple types of elements described in this section.

Given a planar polygonal panel, we define the following unitary vectors following the right-hand rule as follows

```@raw html
<center>
  <img src="../../assets/images/panelcs00.png" alt="Pic here" width="400px">
</center>
```

```math
\begin{align*}
        \text{Tangent vector}\qquad
        &
        \hat{\mathbf{t}} = \frac{\mathbf{p}_2 - \mathbf{p}_1}{\Vert \mathbf{p}_2 - \mathbf{p}_1 \Vert}
    \\
        \text{Normal vector}\qquad
        &
        \hat{\mathbf{n}} = \frac{
                                    \left( \mathbf{p}_2 - \mathbf{p}_1 \right)
                                    \times
                                    \left( \mathbf{p}_3 - \mathbf{p}_1 \right)
                                }{\left\Vert
                                    \left( \mathbf{p}_2 - \mathbf{p}_1 \right)
                                    \times
                                    \left( \mathbf{p}_3 - \mathbf{p}_1 \right)
                                \right\Vert}
    \\
        \text{Oblique vector}\qquad
        &
        \hat{\mathbf{o}} = \hat{\mathbf{n}} \times \hat{\mathbf{t}}
\end{align*}
```

This defines the orthonormal basis $\left( \hat{\mathbf{t}},\, \hat{\mathbf{o}},\, \hat{\mathbf{n}} \right)$ which, along with the panel's centroid, define the panel's local coordinate system.
This basis follows $\hat{\mathbf{t}} \times \hat{\mathbf{o}} = \hat{\mathbf{n}}$, $\hat{\mathbf{o}} \times \hat{\mathbf{n}} = \hat{\mathbf{t}}$, and $\hat{\mathbf{n}} \times \hat{\mathbf{t}} = \hat{\mathbf{o}}$.
