# Problem Formulation

(Adapted from Katz and Plotkin's *Low Speed Aerodynamics*, Sec 3.2)

By [Helmholtz' decomposition theorem](https://en.wikipedia.org/wiki/Helmholtz_decomposition), any velocity field $\mathbf{u}$ can be decomposed into a uniform component, an irrotational component, and a solenoidal component as

```math
\begin{align*}
        \mathbf{u}
    =
        \underbrace{\mathbf{u}_\infty}_\text{uniform}
        + \underbrace{\nabla\phi}_\text{irrotational}
        + \underbrace{\nabla\times\boldsymbol{\psi}}_\text{solenoidal}
,\end{align*}
```

where $\mathbf{u}_\infty(t)$ is the freestream, $\phi (\mathbf{x}, t)$ is a scalar-potential field, and $\boldsymbol\psi (\mathbf{x}, t)$ is a vector-potential field.

>**NOTE:** Even though here we show $\mathbf{u}_\infty$ as its own component, the freestream can be rolled into the potential field as a linearly-varying potential, giving a uniform velocity field.

Assuming incompressible flow, the continuity equation poses a Laplace equation for the scalar-potential as

```math
\begin{align*}
        \nabla \cdot \mathbf{u} = 0
    \quad\Rightarrow\quad
        \nabla^2 \phi = 0
.\end{align*}
```

Before we continue, let us derive a useful integral identity: [Green's second identity](https://en.wikipedia.org/wiki/Green%27s_identities).
Let $f$ and $g$ be two differentiable functions, we have

```math
\begin{align*}
    \nabla \cdot \left( f \nabla g - g \nabla f \right) = f \nabla^2 g - g \nabla^2 f
.\end{align*}
```

The divergence theorem then leads to

```math
\begin{align*}
        \int\limits_{\partial V}
            \left( f \nabla g - g \nabla f \right) \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    & =
        \int\limits_{V}
            \left( f \nabla^2 g - g \nabla^2 f \right)
        \,\mathrm{d}V
,\end{align*}
```
over any volume $V$. Here we have defined the normal $\hat{\mathbf{n}}$  as pointing outward from $V$.

Using Green's second identity, we now define $f(\mathbf{x}) \equiv \frac{1}{\Vert \mathbf{x} - \mathbf{x}_p \Vert^2} = \frac{1}{r}$ and $g(\mathbf{x}) \equiv \phi(\mathbf{x})$ where $r$ is the distance to a point $\mathbf{x}_p$ defining an arbitrary fixed center, and the identity gives

```math
\begin{align*}
        \int\limits_{\partial V}
            \left[ \frac{1}{r} \nabla \phi - \phi \nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    & =
        \int\limits_{V}
            \left[ \frac{1}{r} \cancel{\nabla^2 \phi}^{\,\,0} - \phi \cancel{\nabla^2 \left(\frac{1}{r}\right)}^{\,\,0} \right]
        \,\mathrm{d}V
\end{align*}
```

```math
\begin{align*}
    \Rightarrow \boxed{
        \int\limits_{\partial V}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    = %
        0
    }
.\end{align*}
```

The integrand becomes singular when evaluated at $\mathbf{x}_p$ since $\lim \limits_{r\rightarrow 0 }\frac{1}{r} = \infty$, which we will use to obtain some expressions depending on whether $\mathbf{x}_p$ lays inside, outside, or at the boundary of $V$.

In the case that $\mathbf{x}_p \cancel{\in} V$, the integral equation is automatically satisfied. In the case that $\mathbf{x}_p \in V$, we introduce a hole in $V$ in the form of a sphere of radius $\epsilon$ surrounding $\mathbf{x}_p$, and the integral equation is then automatically satisfied over the domain $V \backslash V_\epsilon$:

```math
\begin{align*}
    &
        \int\limits_{\partial (V \backslash V_\epsilon)}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    = %
        0
    \\
    \Leftrightarrow &
        \int\limits_{\partial V}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
        -
        \int\limits_{\partial V_\epsilon}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    = %
        0
.\end{align*}
```

Since $r$ is centered at $\mathbf{x}_p$ and the normal $\hat{\mathbf{n}}$ points radially inwards the sphere (since all normals point outside of the volume $V$), the second integral becomes

```math
\begin{align*}
        \int\limits_{\partial V_\epsilon}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    = %
        - \int\limits_{\partial V_\epsilon}
            \left( \frac{1}{r} \frac{\partial \phi}{\partial r} + \frac{\phi}{r^2} \right)
        \,\mathrm{d}S
.\end{align*}
```

In the limit that $\epsilon$ is infinitesimally small, the integral can be approximated as

```math
\begin{align*}
        \int\limits_{\partial V_\epsilon}
            \left( \frac{1}{r} \frac{\partial \phi}{\partial r} - \frac{\phi}{r^2} \right)
        \,\mathrm{d}S
    \approx
        - \left(
            \frac{\partial \phi}{\partial r}(\mathbf{x}_p) \lim\limits_{r\rightarrow 0} \frac{1}{r}
            +
            \phi (\mathbf{x}_p) \lim\limits_{r\rightarrow 0} \frac{1}{r^2}
        \right)
        \lim\limits_{\epsilon\rightarrow 0} 4 \pi \epsilon^2
,\end{align*}
```
where $4 \pi \epsilon^2$ is the surface area of the sphere.

Assuming $\phi(\mathbf{x}_p) \neq \pm\infty$ and $\frac{\partial \phi}{\partial r}(\mathbf{x}_p) \neq \pm\infty$, we have

```math
\begin{align*}
        \int\limits_{\partial V_\epsilon}
            \left( \frac{1}{r} \frac{\partial \phi}{\partial r} - \frac{\phi}{r^2} \right)
        \,\mathrm{d}S
    & \approx
            - 4 \pi \frac{\partial \phi}{\partial r}(\mathbf{x}_p) \cancel{\lim\limits_{\rho\rightarrow 0} \frac{\rho ^2}{\rho}}^{\,\,{0}}
            -
            4 \pi \phi(\mathbf{x}_p) \cancel{\lim\limits_{\rho\rightarrow 0} \frac{\rho^2}{\rho^2}}^{\,\,1}
.\end{align*}
```
Thus,

```math
\begin{align*}
        \int\limits_{\partial V_\epsilon}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    & =
        4 \pi \phi(\mathbf{x}_p)
.\end{align*}
```

Substituting this back into the original equation,

```math
\begin{align*}
        \int\limits_{\partial V}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
        -
        4 \pi \phi(\mathbf{x}_p)
    = %
        0
,\end{align*}
```
and since $\mathbf{x}_p$ is any arbitrary point inside $V$, we conclude

```math
\begin{align*}
    \boxed{
        \phi(\mathbf{x})
    = %
        \frac{1}{4 \pi }
        \int\limits_{\partial V}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    }
    , \quad \forall \mathbf{x} \in V
.\end{align*}
```

If $\mathbf{x}_p$ lays in the boundary of $V$ (*i.e.*, $\mathbf{x}_p \in \partial V$), we repeat the same derivation except that only half of the sphere $V_\epsilon$ is contained in the integration domain, and we conclude

```math
\begin{align*}
    \boxed{
        \phi(\mathbf{x})
    = %
        \frac{1}{2 \pi }
        \int\limits_{\partial V}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    }
    , \quad \forall \mathbf{x} \in \partial V
.\end{align*}
```

We merge the two equations to get

```math
\begin{align*}
    \therefore
    \boxed{
        \phi(\mathbf{x})
    = %
        \frac{f_{\tiny \partial V}(\mathbf{x})}{\pi }
        \int\limits_{\partial V}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    }
    , \quad \forall \mathbf{x} \in V
    , \quad \text{where }
    f_{\tiny \partial V}(\mathbf{x}) = \begin{cases}
                                          1/4  & \mathbf{x} \in V \backslash \partial V \\
                                          1/2  & \mathbf{x} \in \partial V
                                        \end{cases}
.\end{align*}
```

Now, consider the case that $V$ is infinitely large in all directions (with bound $\partial V_\infty$) while also having a hole (bound $\partial V_b$, which represents a solid body immersed in the domain), leading to $\partial V = \partial V_b \cup \partial V_\infty$.
In $V$, the potential is calculated as

```math
\begin{align*}
        \phi(\mathbf{x})
    = %
        \frac{f_{\tiny \partial V_b}(\mathbf{x})}{\pi }
        \int\limits_{\partial V_b}
            \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
        +
        \underbrace{
            \frac{1}{4\pi }
            \int\limits_{\partial V_\infty}
                \left[ \frac{1}{r} \nabla \phi - \phi\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
            \,\mathrm{d}S
        }_{
        \equiv \phi_\infty(\mathbf{x})
        }
.\end{align*}
```

When the isolated body is considered as a volume of its own, $V_b$, we can define an internal potential (*i.e.*, internal to the body) as $\phi_i$. When $\mathbf{x}_p$ is defined as being outside of $V_b$ (but still inside $V$), we have

```math
\begin{align*}
        \frac{f_{\tiny \partial V_b}(\mathbf{x})}{\pi }
        \int\limits_{\partial V_b}
            \left[ \frac{1}{r} \nabla \phi_i - \phi_i\nabla \left(\frac{1}{r}\right) \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
    = %
        0
,\end{align*}
```

from the Laplace equation.

Since the superposition of two solutions to the Laplace equation yields another valid solution, we can add this internal potential to the previous equation to obtain

```math
\begin{align*}
        \phi(\mathbf{x})
    = %
        \frac{f_{\tiny \partial V_b}(\mathbf{x})}{\pi }
        \int\limits_{\partial V_b}
            \left[
                \frac{1}{r} \nabla \left(\phi - \phi_i\right)
                -
                \left(\phi - \phi_i\right) \nabla \left(\frac{1}{r}\right)
            \right] \cdot \hat{\mathbf{n}}
        \,\mathrm{d}S
        +
        \phi_\infty(\mathbf{x})
,\end{align*}
```

where the negative sign accompanying $\phi_i$ comes from requiring $\hat{\mathbf{n}}$ to point outward from $V$.

This equation poses a boundary integral equation (BIE), with the potential $\phi$ being fully determined by its value at the boundaries. We can then define

```math
\begin{align*}
    &
    \boxed{-\sigma \equiv \frac{\partial \phi}{\partial n} - \frac{\partial \phi_i}{\partial n}}
    \\ &
    \boxed{-\mu \equiv \phi - \phi_i}
\end{align*}
```

resulting in the following BIE,

```math
\begin{align*}
    \boxed{
            \phi(\mathbf{x})
        = %
            -\frac{f_{\tiny \partial V_b}(\mathbf{x})}{\pi }
            \int\limits_{\partial V_b}
                \left[
                    \sigma \frac{1}{r}
                    -
                    \mu \hat{\mathbf{n}} \cdot  \nabla \left(\frac{1}{r}\right)
                \right]  
            \,\mathrm{d}S
            +
            \phi_\infty(\mathbf{x})
    }
    , \quad \text{where }
    f_{\tiny \partial V_b}(\mathbf{x}) = \begin{cases}
                                          1/4  & \mathbf{x} \in V \backslash \partial V_b \\
                                          1/2  & \mathbf{x} \in \partial V_b
                                        \end{cases}
,\end{align*}
```

where $\mu$ and $\sigma$ are unknowns that we will later determine imposing boundary conditions.
