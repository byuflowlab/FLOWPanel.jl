using PyPlot
using FLOWPanel

TK = FLOWPanel.ConstantDoublet
v1x, v1y, v1z = 0.0, 0.0, 0.0
v2x, v2y, v2z = 0.0, 1.0, 0.0
d1, d2, d3 = 1.0, 0.0, 0.0
strength = 1.0
kerneloffset = 0.0
zs = range(-1,stop=1,length=100)
targets = [[0.5, 0.5, z] for z in zs]
phis = [FLOWPanel._phi_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength; kerneloffset) for target in targets]
vs = [FLOWPanel._U_semiinfinite(target, TK, v1x, v1y, v1z, v2x, v2y, v2z, d1, d2, d3, strength; kerneloffset)[3] for target in targets]

fig = figure("check doublet horseshoe")
fig.clear()
fig.add_subplot(121, xlabel="z", ylabel="phi")
fig.add_subplot(122, xlabel="z", ylabel="Uz")
ax = fig.get_axes()[1]
ax.plot(zs, phis)
ax = fig.get_axes()[2]
ax.plot(zs, vs)
