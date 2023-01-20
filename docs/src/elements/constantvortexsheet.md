# Constant-Strength Vortex Sheet

(Adapted from David [Pate's 2017](https://smartech.gatech.edu/handle/1853/58265) doctoral dissertation *A Surface Vorticity Method for
Wake-Body Interactions*,
and [Pate and German (2018)](https://arc.aiaa.org/doi/abs/10.2514/1.J057120), *A Surface Vorticity Panel Method*)


The velocity induced by a vortex sheet at a target point $\mathbf{x}$ is the integral of the Biot-Savart law over the surface of the sheet:

```math
\begin{align*}
        \mathbf{u}(\mathbf{x})
    =
        \frac{1}{4\pi}
        \int\limits_S
            \frac{\mathbf{r} \times \boldsymbol\gamma}{\Vert \mathbf{r} \Vert^3}
        \,\mathrm{d} S'
,\end{align*}
```

where $\mathrm{r} = \mathbf{x}' - \mathbf{x}$ is the vector directed from the target to the differential $\mathrm{d}S'$.

In order to evaluate this integral, we follow the procedure put forth by [Pate 2017](https://smartech.gatech.edu/handle/1853/58265).
First, we define a local coordinate system at the projection of point $\mathbf{x}$ unto the plane of the triangle as shown below:

```@raw html
<center>
  <img src="../../assets/images/pate-figa1.png" alt="Pic here" style="width:50%;"/>
  <br>
  <span style="color:gray">(figure retrieved from <a href="https://smartech.gatech.edu/handle/1853/58265">Pate 2017</a>)</span>
  <br><br>
</center>
```

The projection point $\mathbf{p}_x$ (or $P_B$ in the figure above) is calculated as $\mathbf{p}_x = \mathbf{x} - \left[ (\mathbf{x} - \mathbf{p}_1)\cdot \hat{\mathbf{n}} \right] \hat{\mathbf{n}}$, and the orthonormal basis $(\mathbf{b}_1, \mathbf{b}_2, \mathbf{b}_3)$  is completely arbitrary except for requiring $\mathbf{b}_3 = \hat{\mathbf{n}}$.
We also define the vector directed from $\mathbf{p}_x$ to the differential as $\boldsymbol\rho \equiv \mathbf{x}' - \mathbf{p}_x$.
We then define the following three triangles $S_1, S_2, S_3$ that will be used for the integration:

```@raw html
<center>
  <img src="../../assets/images/pate-figa2.png" alt="Pic here" style="width:35%;"/>
  <br>
  <span style="color:gray">(figure retrieved from <a href="https://smartech.gatech.edu/handle/1853/58265">Pate 2017</a>)</span>
  <br><br>
</center>
```
where, $\mathbf{q}_i = \mathbf{p}_i - \mathbf{p}_x$.
$S_1$ is the triangle of vertices $(\mathbf{p}_x, \mathbf{p}_1, \mathbf{p}_2)$, $S_2$ is $(\mathbf{p}_x, \mathbf{p}_2, \mathbf{p}_3)$, and $S_3$ is $(\mathbf{p}_x, \mathbf{p}_3, \mathbf{p}_1)$.
More generally, the vertices of triangle $S_i$ in the global coordinate system are $(\mathbf{p}_x, \mathbf{p}_i, \mathbf{p}_{i+1})$, with $\mathbf{p}_{4} = \mathbf{p}_{1}$.
In the coordinate system with origin at $\mathbf{p}_x$ and basis $(\mathbf{b}_1, \mathbf{b}_2, \mathbf{b}_3)$, the vertices of triangle $S_i$ are

```math
\begin{align*}
    \boxed{
        (\mathbf{0}, \mathbf{q}_1, \mathbf{q}_{2})_i
    =
        (\mathbf{p}_x, \mathbf{p}_i - \mathbf{p}_x, \mathbf{p}_{i+1} - \mathbf{p}_x)
    }
\end{align*}
```

For each triangle $S_i$, Pate then proceeds to define the following orthonormal basis $(\mathbf{c}_1, \mathbf{c}_2, \mathbf{c}_3)$:


```@raw html
<center>
  <img src="../../assets/images/pate-figa3.png" alt="Pic here" style="width:35%;"/>
  <br>
  <span style="color:gray">(figure retrieved from <a href="https://smartech.gatech.edu/handle/1853/58265">Pate 2017</a>)</span>
  <br><br>
</center>
```

```math
\begin{align*}
            \mathbf{c}_3
        & \equiv
            \mathbf{b}_3 = \hat{\mathbf{n}}
    \\
            \mathbf{c}_2
        & \equiv
            \frac{\mathbf{q}_2 - \mathbf{q}_1}{\Vert \mathbf{q}_2 - \mathbf{q}_1 \Vert}
    \\
            \mathbf{c}_1
        & \equiv
            \mathbf{c}_2 \times \mathbf{c}_3
\end{align*}
```

The vector $\boldsymbol\rho$ is then decomposed in this basis as

```math
\begin{align*}
        \boldsymbol\rho
    & =
        u \mathbf{c}_1 + v \mathbf{c}_2
\end{align*}
```
where
```math
\begin{align*}
        u
    & \equiv
        \boldsymbol\rho \cdot \mathbf{c}_1
    \\
        v
    & \equiv
        \boldsymbol\rho \cdot \mathbf{c}_2
.\end{align*}
```

The full-distance vector $\mathbf{r}$ (pointing from $\mathbf{x}$ to the differential) is decomposed as

```math
\begin{align*}
        \mathbf{r}
    & =
        u \mathbf{c}_1 + v \mathbf{c}_2 - z \mathbf{c}_3
\end{align*}
```
where
```math
\begin{align*}
        z
    & \equiv
        \left( \mathbf{x} - \mathbf{p}_1 \right) \cdot \hat{\mathbf{n}}
,\end{align*}
```

while the constant strength of the vortex sheet is decomposed as

```math
\begin{align*}
        \boldsymbol\gamma
    & =
        a_{00} \mathbf{c}_1 + b_{00} \mathbf{c}_2
,\end{align*}
```
where
```math
\begin{align*}
        a_{00}
    & \equiv
        \boldsymbol\gamma \cdot \mathbf{c}_1
    \\
        b_{00}
    & \equiv
        \boldsymbol\gamma \cdot \mathbf{c}_2
.\end{align*}
```

The numerator in the integral of the Biot-Savart law then becomes
```math
\begin{align*}
        \mathbf{r} \times \boldsymbol\gamma
    & =
        z b_{00} \mathbf{c}_1 - z a_{00} \mathbf{c}_2 + (u b_{00} - v a_{00}) \mathbf{c}_3
.\end{align*}
```

In order to avoid the singularity when the denominator of the integral becomes zero, we add a small offset $\delta$ to the distance of the denominator defining $\boxed{r \equiv \sqrt{u^2 + v^2 + z^2 + \delta^2}}$.
We also add the same offset to the height of the target as $\boxed{h \equiv \sqrt{z^2 + \delta^2}}$, effectively defining a thickness for the sheet.

The integral over the original triange $S$ can then be expressed in terms of $S_1$, $S_2$, and $S_3$ as

```math
\begin{align*}
        \int\limits_S
            \frac{\mathbf{r} \times \boldsymbol\gamma}{r^3}
        \,\mathrm{d} S'
    =
        \sum\limits_i
        \int\limits_{S_i}
            \frac{\mathbf{r} \times \boldsymbol\gamma}{r^3}
        \,\mathrm{d} S_i'
,\end{align*}
```

The integral terms $H_{00}$, $H_{10}$, and $H_{01}$ are given in Pate 2017, Appendix A.4:
```math
\begin{align*}
    &
           H_{00}
       =
           \frac{1}{h} \tan^{-1}\left.\left(
               \frac{al}{ a^2 + h^2 + h\sqrt{l^2 + a^2 + h^2}}
           \right)\right\rvert_{l_1}^{l_2}
    \\ &
           H_{10}
       =
           \left.
           \left[
               \frac{l}{\sqrt{l^2 + a^2}} \ln\left( \sqrt{l^2 + a^2 + h^2} + \sqrt{l^2 + a^2} \right)
               -
               \ln\left( l + \sqrt{l^2 + a^2 + h^2} \right)
               -
               \frac{l \ln h }{\sqrt{l^2 + a^2}}
           \right]
           \right\rvert_{l_1}^{l_2}
    \\ &
           H_{01}
       =
           \left.
           \left\{
               \frac{a}{\sqrt{l^2 + a^2}}
               \left[
                   \ln h
                   -
                   \ln\left( \sqrt{l^2 + a^2 + h^2} + \sqrt{l^2 + a^2} \right)
               \right]
           \right\}
           \right\rvert_{l_1}^{l_2}
\end{align*}
```
and the integral over $S_i$ is calculated as

```math
\begin{align*}
            \int\limits_{S_i}
                \frac{\mathbf{r} \times \boldsymbol\gamma}{r^3}
            \,\mathrm{d} S_i'
        =
            z b_{00} \mathbf{c}_1
            \underbrace{
                \int\limits_{S_i}  
                    \frac{1}{r^3}
                \,\mathrm{d} S_i'
            }_{H_{00}}
        -
            z a_{00} \mathbf{c}_2
            \underbrace{
                \int\limits_{S_i}
                    \frac{1}{r^3}
                \,\mathrm{d} S_i'
            }_{H_{00}}
        +
            b_{00} \mathbf{c}_3
            \underbrace{
                \int\limits_{S_i}
                    \frac{u}{r^3}
                \,\mathrm{d} S_i'
            }_{H_{10}}
        -
            a_{00} \mathbf{c}_3
            \underbrace{
                \int\limits_{S_i}
                    \frac{v}{r^3}
                \,\mathrm{d} S_i'
            }_{H_{01}}
.\end{align*}
```

The integral terms $H_{00}$, $H_{10}$, and $H_{01}$ are given in Pate 2017, Appendix A.4:
```math
\begin{align*}
    &
           H_{00}
       =
           \frac{1}{h} \tan^{-1}\left.\left(
               \frac{al}{ a^2 + h^2 + h\sqrt{l^2 + a^2 + h^2}}
           \right)\right\rvert_{l_1}^{l_2}
    \\ &
           H_{10}
       =
           \left.
           \left[
               \frac{l}{\sqrt{l^2 + a^2}} \ln\left( \sqrt{l^2 + a^2 + h^2} + \sqrt{l^2 + a^2} \right)
               -
               \ln\left( l + \sqrt{l^2 + a^2 + h^2} \right)
               -
               \frac{l \ln h }{\sqrt{l^2 + a^2}}
           \right]
           \right\rvert_{l_1}^{l_2}
    \\ &
           H_{01}
       =
           \left.
           \left\{
               \frac{a}{\sqrt{l^2 + a^2}}
               \left[
                   \ln h
                   -
                   \ln\left( \sqrt{l^2 + a^2 + h^2} + \sqrt{l^2 + a^2} \right)
               \right]
           \right\}
           \right\rvert_{l_1}^{l_2}
\end{align*}
```
where
```math
\begin{align*}
    &
       a = \mathbf{q}_1 \cdot \mathbf{c}_1
   \\ &
       l_1 = \mathbf{q}_1 \cdot \mathbf{c}_2
   \\ &
       l_2 = \mathbf{q}_2 \cdot \mathbf{c}_2
   \\ &
       h = \sqrt{z^2 + \delta^2}
   \\ &
       z = \left( \mathbf{x} - \mathbf{p}_1 \right) \cdot \hat{\mathbf{n}}
\end{align*}
```
> **NOTE:** $H_{10}$ and $H_{01}$ becomes undefined when the projection of $\mathbf{x}$ onto the panel plane lays on one of the vertices, leading to the a division of zero by zero in $\frac{l}{\sqrt{l^2 + a^2}}$ and $\frac{a}{\sqrt{l^2 + a^2}}$. Both $H_{10}$ and $H_{01}$ actually converge to 0 in that situation, and can be manually prescibed to that value.

Finally, the velocity induced by the sheet is evaluated as
```math
\begin{align*}
    \boxed{
            \mathbf{u}(\mathbf{x})
        =
            \frac{1}{4\pi}
            \int\limits_S
                \frac{\mathbf{r} \times \boldsymbol\gamma}{\Vert \mathbf{r} \Vert^3}
            \,\mathrm{d} S'
        =
            \frac{1}{4\pi}
            \sum\limits_i
            \left[
                z b_{00} H_{00} \mathbf{c}_1
                -
                z a_{00} H_{00} \mathbf{c}_2
                +
                \left(b_{00} H_{10} - a_{00} H_{01}\right) \mathbf{c}_3
            \right]
    \,
    }
,\end{align*}
```

To speed up computation, we can use [trigonometric](https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Arctangent_addition_formula) and [logarithmic](https://en.wikipedia.org/wiki/Logarithm#Logarithmic_identities) identities to rewrite the formulas given by Pate as
```math
\begin{align*}
    &
           H_{00}
       =
           \frac{1}{h} \tan^{-1}\left(
               \frac{m_2 - m_1}{1 + m_1 m_2}
           \right)
           \qquad \text{ where } \quad
           m_i = \frac{al_i}{ a^2 + h^2 + h\sqrt{l_i^2 + a^2 + h^2}}
    \\ &
           H_{10}
       =
           \ln\left( \frac{m_2}{m_1} \right)
           \qquad \text{ where } \quad
           m_i = \frac{1}{l_i + \sqrt{l_i^2 + a^2 + h^2}}
               \left(
                   \frac{\sqrt{l_i^2 + a^2 + h^2} + \sqrt{l_i^2 + a^2}}{h}
               \right)^\frac{l_i}{\sqrt{l_i^2 + a^2}}
    \\ &
           H_{01}
       =
               -
               \ln \left(
                   \frac{m_2}{m_1}
               \right)
           \qquad \text{ where } \quad
           m_i = \left( \frac{\sqrt{l_i^2 + a^2 + h^2} + \sqrt{l_i^2 + a^2}}{h} \right)^\frac{a}{\sqrt{l_i^2 + a^2}}
,\end{align*}
```
thus reducing the number of $\log$ and $\tan^{-1}$ evaluations.


```@raw html
<center>
  <table>
      <tr>
          <td>
              <img src="../../assets/images/vortexsheet-viz-w02.png" alt="Pic here" width="100%">
          </td>
          <td>
              <img src="../../assets/images/vortexsheet-viz-u02.png" alt="Pic here" width="100%">
          </td>
      </tr>
  </table>
</center>
```


```@raw html
<center>
    <img src="../../assets/images/constantvortexsheet01-probe1.png" alt="Pic here" width="90%">
    <img src="../../assets/images/constantvortexsheet01-probe2.png" alt="Pic here" width="90%">
    <img src="../../assets/images/constantvortexsheet01-probe3.png" alt="Pic here" width="90%">
</center>
```
