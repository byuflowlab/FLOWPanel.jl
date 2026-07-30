# Surface Vorticity Force

`SurfaceVorticityForce` is a diagnostic force monitor that estimates the
surface vortex sheet from the panel strength field and integrates a
control-point Kutta-Joukowski force.

For a doublet or vortex-ring panel strength `mu`, the tangential surface
gradient defines the sheet-strength vector

```math
\kappa = -n \times \nabla_s \mu .
```

The leading minus sign follows FLOWPanel's stored strength convention. In many
panel-method references the doublet strength is defined as the exterior-minus-
interior potential jump; FLOWPanel's stored `mu`/`gamma` has the opposite sign,
so the physical sheet-strength vector uses `-∇s mu` under this convention.

The monitor reports the body-side force contribution

```math
dF = \rho \left(V_{cp} \times \kappa\right) dS ,
```

where `V_cp` is the already-stored relative velocity `body.velocity` at the panel control point.
It does not add `body.velocity_kinematic`.

```julia
monitor = SurfaceVorticityForce(body, nt, i_system;
    rho=1.225,
    i_frame=-1,
    normalization=NoNormalization(),
    correct_kuttacondition=true)
```

The monitor selects the `"gamma"` strength column when present, otherwise
`"mu"`. Bodies without either strength are rejected. Each call overwrites
`body.F`, so ordering matters when combining it with pressure-based force
monitors. It provides `:F` and has no `:P` dependency.
