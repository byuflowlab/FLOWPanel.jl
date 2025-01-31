# 2D Elements

## Linear-Strength Vortex Panel
*(Adapted from Katz and Plotkin's Low Speed Aerodynamics, Secs 10.3.3 and 11.4.2a)*

A linear-strenght vortex panel is composed of a constant strength panel
superimposed with a panel that ramps up its strength between 0 and a given value.
Hence, to define the linear-strength vortex panel first we define this two
components as follows.

### Constant Strength

Given a vortex element laying between points $\mathbf{p}_1$ and $\mathbf{p}_2$ and constant strength $\gamma_a$,
we wish to calculate the velocity induced by the element at an arbitrary point $\mathbf{p}_3$.

The induced velocity is given in Eqs. (10.39) and (10.40) of Katz & Plotkin as follows:

```math
\begin{align*}
        u
    & =
        \frac{\gamma_a}{2\pi} \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
    \\
        v
    & =
        -
        \frac{\gamma_a}{4\pi} \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
\end{align*}
```

```@raw html
<center>
  <br><i>2D vortex panel with $\gamma_1 = 1$ and $\gamma_2 = 1$</i>
  <table>
      <tr>
          <td>
              <img src="../../assets/images/vortexpanel2d-constant000.png" alt="Pic here" width="100%">
          </td>
          <td>
              <img src="../../assets/images/vortexpanel2d-constant001.png" alt="Pic here" width="100%">
          </td>
      </tr>
  </table>
</center>
```


```@raw html
<center>
    <img src="../../assets/images/vortexpanel-viz000-probe1.png" alt="Pic here" width="90%">
    <img src="../../assets/images/vortexpanel-viz000-probe2.png" alt="Pic here" width="90%">
</center>
```

### Simple Linear Distribution

Given a vortex element laying between points $\mathbf{p}_1$ and $\mathbf{p}_2$ and strength varying linearly from $0$ to $\gamma_b$ between $x_1$ and $x_2$ as given by
```math
\begin{align*}
    \gamma(x) = \gamma_b x
,\end{align*}
```
we wish to calculate the velocity induced by the element at an arbitrary point $\mathbf{p}_3$.

The induced velocity is given in Eqs. (10.72) and (10.73) of Katz & Plotkin as follows:

```math
\begin{align*}
        u
    & =
        -\frac{\gamma_b}{4\pi}
        \left[
            y \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
            -
            2x \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
        \right]
    \\
        v
    & =
        -\frac{\gamma_b}{2\pi}
        \left[
            \frac{x}{2} \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
            +
            (x_1-x_2)
            +
            y \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
        \right]
\end{align*}
```


```@raw html
<center>
  <br><i>2D vortex panel with $\gamma_1 = 0$ and $\gamma_2 = 1$</i>
  <table>
      <tr>
          <td>
              <img src="../../assets/images/vortexpanel2d-linear014.png" alt="Pic here" width="100%">
          </td>
          <td>
              <img src="../../assets/images/vortexpanel2d-linear015.png" alt="Pic here" width="100%">
          </td>
      </tr>
  </table>
</center>
```


```@raw html
<center>
    <img src="../../assets/images/vortexpanel-viz001-probe1.png" alt="Pic here" width="90%">
    <img src="../../assets/images/vortexpanel-viz001-probe2.png" alt="Pic here" width="90%">
    <img src="../../assets/images/vortexpanel-viz001-probe3.png" alt="Pic here" width="90%">
</center>
```


### Compound Linear Distribution

```@raw html
<center>
  <table>
      <tr>
          <td>
            <center>
                <img src="../../assets/images/2Dlinearvortex.png" alt="Pic here" width="30%">
              </center>
          </td>
      </tr>
  </table>
</center>
```

Given a vortex element laying between points $\mathbf{p}_1$ and $\mathbf{p}_2$ and strength varying linearly from $\gamma_1$ to $\gamma_2$ as,
```math
\begin{align*}
    \gamma(x) = \gamma_1 + (\gamma_2 - \gamma_1) \frac{x - x_1}{x_2 - x_1}
\end{align*}
```
we wish to calculate the velocity induced by the element at an arbitrary point $\mathbf{p}_3$.

Notice that the strength distribution can be reduced to a superposition of the two previous cases with
```math
\begin{align*}
        \gamma_a
    & =
        \gamma_1 - (\gamma_2 - \gamma_1) \frac{x_1}{x_2 - x_1}
    =
        \frac{\gamma_1 x_2 - \gamma_2 x_1}{x_2 - x_1}
    \\
        \gamma_b
    & =
        \frac{\gamma_2 - \gamma_1}{x_2 - x_1}
\end{align*}
```

Expressing $\mathrm{p}_1$, $\mathrm{p}_2$, and $\mathrm{p}_3$ in terms of the element's coordinate system shown above, the velocity induced at $\mathrm{p}_3$ by the constant-strength component is
```math
\begin{align*}
        u_\mathrm{const}
    & =
        \frac{1}{2\pi} \frac{\gamma_1 x_2 - \gamma_2 x_1}{x_2 - x_1}
        \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
    \\
        v_\mathrm{const}
    & =
        -
        \frac{1}{4\pi} \frac{\gamma_1 x_2 - \gamma_2 x_1}{x_2 - x_1}
        \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
\end{align*}
```
and the velocity induced by the linear-strength component is
```math
\begin{align*}
        u_\mathrm{lin}
    & =
        -\frac{1}{4\pi} \frac{\gamma_2 - \gamma_1}{x_2 - x_1}
        \left[
            y \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
            -
            2x \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
        \right]
    \\
        v_\mathrm{lin}
    & =
        -\frac{1}{2\pi} \frac{\gamma_2 - \gamma_1}{x_2 - x_1}
        \left[
            \frac{x}{2} \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
            +
            (x_1-x_2)
            +
            y \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
        \right]
.\end{align*}
```
The total induced velocity is given as the superposition
```math
\begin{align*}
        u
    & =
        u_\mathrm{const} + u_\mathrm{lin}
    \\
        v
    & =
        v_\mathrm{const} + v_\mathrm{lin}
\end{align*}
```
or
```math
\begin{align*}
        u
    & =
        \frac{1}{2\pi} \frac{\gamma_1 x_2 - \gamma_2 x_1}{x_2 - x_1}
        \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
        -
        \frac{1}{4\pi} \frac{\gamma_2 - \gamma_1}{x_2 - x_1}
        \left[
            y \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
            -
            2x \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
        \right]
    \\
        v
    & =
        -
        \frac{1}{4\pi} \frac{\gamma_1 x_2 - \gamma_2 x_1}{x_2 - x_1}
        \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
        -
        \frac{1}{2\pi} \frac{\gamma_2 - \gamma_1}{x_2 - x_1}
        \left[
            \frac{x}{2} \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
            +
            (x_1-x_2)
            +
            y \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
        \right]
\end{align*}
```
which simplifies to
```math
\begin{align*}
        u
    & =
        -
        \frac{y}{4\pi} \frac{\gamma_2 - \gamma_1}{x_2 - x_1}
        \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
        +
        \frac{1}{2\pi} \frac{\gamma_1 (x_2-x)+ \gamma_2 (x-x_1)}{x_2 - x_1}
        \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
    \\
        v
    & =
        -
        \frac{1}{4\pi}
        \frac{\gamma_1 (x_2 - x) + \gamma_2 (x - x_1)}{x_2 - x_1}
        \ln\frac{(x-x_1)^2 + y^2}{(x-x_2)^2 + y^2}
        -
        \frac{y}{2\pi} \frac{\gamma_2 - \gamma_1}{x_2 - x_1}
        \left( \tan^{-1}\frac{y}{x-x_2} - \tan^{-1}\frac{y}{x-x_1} \right)
        +
        \frac{\gamma_2 - \gamma_1}{2\pi}
\end{align*}
```

These equations are simplified after going back to the vector definition of points $\mathbf{p}_1$, $\mathbf{p}_2$, and $\mathbf{p}_3$, and defining $r_1 \equiv \Vert \mathbf{p}_3 - \mathbf{p}_1 \Vert$ and $r_2 \equiv \Vert \mathbf{p}_3 - \mathbf{p}_2 \Vert$:

```math
\begin{align*}
    \bullet \quad &
        u
    =
        -
        \frac{\Delta \gamma}{2\pi} \frac{y}{d}
        \ln\frac{r_1}{r_2}
        +
        \frac{1}{2\pi} \left( \gamma_1 + \Delta \gamma \frac{x-x_1}{d} \right)
        \left( \theta_2 - \theta_1 \right)
    \\
    \bullet \quad &
        v
    =
        -
        \frac{\Delta \gamma}{2\pi} \frac{y}{d}
        \left( \theta_2 - \theta_1 \right)
        -
        \frac{1}{2\pi}
        \left( \gamma_1 + \Delta \gamma \frac{x-x_1}{d} \right)
        \ln\frac{r_1}{r_2}
        +
        \frac{\Delta \gamma}{2\pi}
    \\
    \bullet \quad &
        \Delta \gamma
    =
        \gamma_2 - \gamma_1
    \\
    \bullet \quad &
        d
    =
        x_2 - x_1
    \\
    \bullet \quad &
        \theta_1
    =
        \tan^{-1}\frac{y}{x-x_1}
    \\
    \bullet \quad &
        \theta_2
    =
        \tan^{-1}\frac{y}{x-x_2}
\end{align*}
```

> **NOTES**
> * The vortex panel as defined here is such that a panel with $+\gamma$
>   strength is a vortex aligned with -z, which is counterintuitive. However,
>   $u$ at $y = +0 $ (slight offset in the positive normal direction)
>   is $+\frac{\gamma}{2}$, which is desirable.
> * There is a discontinuity on the panel at $y = \pm 0$ which requires a small offset in the implementation to avoid random behaviors with floating point precision.
> * There are singularities at $(x, y)=(x_1, 0)$ and $(x, y)=(x_2, 0)$, however, we see that $\lim\limits_{x\rightarrow_- x_1} u(x, \pm 0) = 0$, $\lim\limits_{x\rightarrow_+ x_1} u(x, \pm 0) = \pm \frac{\gamma_1}{2}$, $\lim\limits_{x\rightarrow_- x_2} u(x, \pm 0) = \pm \frac{\gamma_2}{2}$, $\lim\limits_{x\rightarrow_+ x_2} u(x, \pm 0) = 0$, which need to be harcoded with logic in the implementation.
> * Singularities can be regularized by defining the regularizing factor $\delta$ and $r_1 \equiv \Vert \mathbf{p}_3 - \mathbf{p}_1  + \delta \Vert$ and $r_2 \equiv \Vert \mathbf{p}_3 - \mathbf{p}_2 + \delta \Vert$

We now evaluate the velocity at $y = 0$ to obtain the velocity induced by the panel on itself.
Notice that $\theta_1$ and $\theta_2$ are the angles from the vertices to the probe, which, at $y = \pm 0$, they are $0$ outside of the panel and $\pm \pi$ on between the panel. Hence, we get
```math
\begin{align*}
        u(x, \pm 0)
    & =
        \pm \frac{\gamma_1}{2}
        \pm \frac{\Delta\gamma}{2}\frac{x - x_1}{d}
    \\
        v(x, 0)
    & =
        -
        \frac{1}{2\pi}
        \left( \gamma_1 + \Delta \gamma \frac{x-x_1}{d} \right)
        \ln\frac{r_1}{r_2}
        +
        \frac{\Delta \gamma}{2\pi}
\end{align*}
```
Furthermore, at the center of the panel ($x = \frac{x_1 + x_2}{2}$)
```math
\begin{align*}
        u(\frac{x_1 + x_2}{2}, \pm 0)
    & =
        \pm \frac{\gamma_1 + \gamma_2}{4}
    \\
        v(\frac{x_1 + x_2}{2}, 0)
    & =
        \frac{\gamma_2 - \gamma_1}{2\pi}
\end{align*}
```

```@raw html
<center>
  <br><i>2D vortex panel with $\gamma_1 = -1$ and $\gamma_2 = 1$</i>
  <table>
      <tr>
          <td>
              <img src="../../assets/images/vortexpanel2d-compound000.png" alt="Pic here" width="100%">
          </td>
          <td>
              <img src="../../assets/images/vortexpanel2d-compound001.png" alt="Pic here" width="100%">
          </td>
      </tr>
  </table>
</center>
```


```@raw html
<center>
    <img src="../../assets/images/vortexpanel-viz002-probe1.png" alt="Pic here" width="90%">
    <img src="../../assets/images/vortexpanel-viz002-probe2.png" alt="Pic here" width="90%">
</center>
```
