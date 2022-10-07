# Constant-Strength Doublet (Vortex Ring)

Recall that the potential field induced by doublets is given by

```math
\begin{align*}
        \phi_\mu(\mathbf{x})
    = %
        \frac{f_{\tiny \partial V_b}(\mathbf{x})}{\pi }
        \int\limits_{\partial V_b}
            \mu \hat{\mathbf{n}} \cdot  \nabla \left(\frac{1}{r}\right)
        \,\mathrm{d}S
,\end{align*}
```

while the component induced by sources is

```math
\begin{align*}
        \phi_\sigma(\mathbf{x})
    = %
        -\frac{f_{\tiny \partial V_b}(\mathbf{x})}{\pi }
        \int\limits_{\partial V_b}
            \sigma \frac{1}{r}
        \,\mathrm{d}S
.\end{align*}
```

Assuming constant-strength panels, each strength comes out of the integrals and we rewrite the above equations as

```math
\begin{align*}
        \phi_\mu(\mathbf{x})
    = %
        \mu G_\mu (\mathbf{x})
    \quad , \quad
        G_\mu (\mathbf{x})
    \equiv
        \frac{f_{\tiny \partial V_b}(\mathbf{x})}{\pi }
        \int\limits_{\partial V_b}
             \hat{\mathbf{n}} \cdot  \nabla \left(\frac{1}{r}\right)
        \,\mathrm{d}S
,\end{align*}
```

and

```math
\begin{align*}
        \phi_\sigma(\mathbf{x})
    = %
        \sigma G_\sigma (\mathbf{x})
    \quad , \quad
        G_\sigma (\mathbf{x})
    \equiv
        -\frac{f_{\tiny \partial V_b}(\mathbf{x})}{\pi }
        \int\limits_{\partial V_b}
            \frac{1}{r}
        \,\mathrm{d}S
,\end{align*}
```

Notice that computation of the potential field induced by the doublet element is similar to the computation of the normal velocity induced by the source element:

```math
\begin{align*}
        \hat{\mathbf{n}} \cdot \nabla\phi_\sigma
    =
        -\sigma \frac{f_{\tiny \partial V_b}}{\pi }
        \int\limits_{\partial V_b}
            \hat{\mathbf{n}} \cdot \nabla \left(
                \frac{1}{r}
            \right)
        \,\mathrm{d}S
    =
        -\sigma G_\mu.
\end{align*}
```

Hence, we can simply reuse the computation of the source-induced velocity to calculate the potential of the doublet element as

```math
\begin{align*}
        \phi_\mu
    =
        -\mu \frac{\hat{\mathbf{n}} \cdot \mathbf{u}_\sigma}{\sigma}
,\end{align*}
```

where $\hat{\mathbf{n}} \cdot \mathbf{u}_\sigma = U_z$ with $U_z$ as defined in the previous section.

In order to calculate the velocity field induced by the doublet, instead of calculating $\nabla\phi_\mu$ we take advantage of the fact that the velocity induced by a constant-strength doublet panel is the same than the velocity induced by a vortex ring (see Katz and Plotkin Sec. 10.4.3),

```math
\begin{align*}
        \mathbf{u}_\mathrm{ring} \left( \mathbf{x} \right)
    =
        \frac{\Gamma}{4\pi}
        \sum\limits_{i,j\in A}
            \frac{\mathbf{r}_i \times \mathbf{r}_j}{ \Vert \mathbf{r}_i \times \mathbf{r}_j \Vert^2}
            \mathbf{r}_{ij} \cdot \left(
                \frac{\mathbf{r}_i}{r_i} - \frac{\mathbf{r}_j}{r_j}
            \right)
,\end{align*}
```

when $\Gamma = \mu$, and where $\mathbf{r}_{ij} = \mathbf{p}_j-\mathbf{p}_i$, $\mathbf{r}_i = \mathbf{x} - \mathbf{p}_i$, $\mathbf{r}_j = \mathbf{x} - \mathbf{p}_j$, $A = \{(1,2),\,\dots,\,(n-1, n),\,(n, 1) \}$ and $n$ the number of vertices that make the panel, each with position $\mathbf{p}_k$.

> **NOTE:** For the same reasons explained in the previous section, $\phi_\mu$ poses a discontinuity at the surface since $\lim\limits_{z\rightarrow \pm 0} \phi_\mu (0, 0, z) = -\frac{\mu}{\sigma} \lim\limits_{z\rightarrow \pm 0} U_z (0, 0, z) = \mp \frac{\mu}{2}$. In FLOWPanel, we let $\phi_\mu (0, 0, 0) = 0$, but we also have shifted all control points slightly in the direction of $\hat{\mathbf{n}}$. Remembering that $\hat{\mathbf{n}}_\mathrm{HS} = -\hat{\mathbf{n}}$, the control points are thus shifted in the $-z$ direction, effectively obtaining $\boxed{\phi_\mu (\mathbf{x}_\mathrm{cp}) \approx \frac{\mu}{2}}$.

**OBSERVATIONS**
* The velocity induced by a segment of the vortex ring shown above becomes singular when the denominator terms $r_i$, $r_j$, or $\Vert \mathbf{r}_i \times \mathbf{r}_j \Vert$  approach $0$. For this reason, FLOWPanel adds a small epsilon to each of these terms, while also defining a cutoff threshold for $\Vert \mathbf{r}_i \times \mathbf{r}_j \Vert$ under which the velocity induced becomes $0$ (thus, self-induced velocity is forced to be zero).

> **NOTE:** The small offset added to the denominator terms $r_i$, $r_j$, and $\Vert \mathbf{r}_i \times \mathbf{r}_j \Vert$ corresponds to `body.kerneloffset`, while the cutoff threshold for self-induced velocity corresponds to `body.kernelcutoff`.

The potential and velocity field of a doublet panel of unitary strength ($\mu=1$) is shown below

```@raw html
<center>
  <table>
      <tr>
          <th>
              <img src="../../assets/images/panel-doublet-phi02.png" alt="Pic here" width="450px">
          </th>
          <th>
              <img src="../../assets/images/panel-doublet-u01.png" alt="Pic here" width="300px">
          </th>
      </tr>
  </table>
</center>
```


```@raw html
<center>
  <br>$\nabla \phi$ = $\mathbf{u}$ verification
  <img src="../../assets/images/panel-doublet-divphivsu00.png" alt="Pic here" style="width: 900px;"/>
</center>
```
