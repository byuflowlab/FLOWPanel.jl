# Non-Planar Semi-Infinite Doublet (Vortex Horseshoe)

Suppose that we have a semi-infinite panel with a outgoing semi-infinite direction different than the incoming direction. That means that the semi-infinite panel is no longer planar. In this circumstances, the induced velocity $\mathbf{u}$ is computed just as explained in the previous section since it's simply a horseshoe. The computation of the potential field, however, needs a some adaptation.

From $\mathbf{p}_b$, we split the panel into two, and create two planar sections as shown below. The computation is then done on the $-\hat{\mathbf{d}}_a,\, \mathbf{p}_i,\, \mathbf{p}_j, +\hat{\mathbf{d}}_a$ section as explained before, while the $-\hat{\mathbf{d}}_a,\, \mathbf{p}_j, +\hat{\mathbf{d}}_b$ is approximated numerically with a large panel.

```@raw html
<center>
  <img src="../../assets/images/semiinfinite-nonplanar-doublet02.png" alt="Pic here" width="450px">
</center>
```

The potential and velocity field of a non-planar semi-infinte doublet panel (or non-planar vortex horseshoe) of unitary strength ($\mu=1$ or $\Gamma=1$) is shown below

```@raw html
<center>
  <table>
      <tr>
          <th>
              <img src="../../assets/images/panel-semiinfinite-nonplanar-doublet-phi00.png" alt="Pic here" width="450px">
          </th>
          <th>
              <img src="../../assets/images/panel-semiinfinite-nonplanar-doublet-u00.png" alt="Pic here" width="450px">
          </th>
      </tr>
  </table>
</center>
```

```@raw html
<center>
  <br>$\nabla \phi$ = $\mathbf{u}$ verification
  <img src="../../assets/images/panel-semiinfinitenonplanardoublet-divphivsu00.png" alt="Pic here" style="width: 700px;"/>
</center>
```
