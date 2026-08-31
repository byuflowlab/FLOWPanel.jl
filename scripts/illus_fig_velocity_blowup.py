#!/usr/bin/env python3
"""fig_velocity_blowup data: peak induced velocity and peak velocity gradient of
a single Gaussian-regularized vortex particle of UNIT strength, as a function of
its core size sigma.

For the Winckelmans/FLOWVPM Gaussian-erf kernel

    g(rho) = erf(rho/sqrt2) - sqrt(2/pi) rho exp(-rho^2/2),   rho = r/sigma
    u(x)   = (1/4pi) g(rho) (Gamma x r)/r^3

the field is self-similar in rho, so with |Gamma| = 1 the space-maxima are exact
power laws

    max|u|      = C_u / sigma^2
    max|grad u| = C_g / sigma^3

with C_u, C_g pure numbers computed here by maximizing over (rho, theta).
"""
import csv
import math
import os

import numpy as np
from scipy.special import erf

OUT = os.path.expanduser(
    "~/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations")
DATA = os.path.join(OUT, "fig_velocity_blowup")
os.makedirs(DATA, exist_ok=True)


def g(rho):
    """Gaussian-erf regularization function."""
    return erf(rho / math.sqrt(2.0)) - math.sqrt(2.0 / math.pi) * rho * np.exp(
        -0.5 * rho ** 2)


def dg(rho):
    """g'(rho) = sqrt(2/pi) rho^2 exp(-rho^2/2)."""
    return math.sqrt(2.0 / math.pi) * rho ** 2 * np.exp(-0.5 * rho ** 2)


def shape_constants_explicit(n_rho=4001, n_theta=181):
    """Reference build of C_u, C_g from the full tensor, for cross-checking.

    With u_i = (1/4pi) eps_ijk Gamma_j r_k K(r), K = g(r/sigma)/r^3,
        du_i/dx_l = (1/4pi) eps_ijk Gamma_j [delta_kl K + r_k r_l K'/r].
    Gamma = e_z, sigma = 1; the configuration is axisymmetric about Gamma, so
    the Frobenius norm depends only on (rho, theta) -- sample that plane.
    """
    rho = np.linspace(1e-6, 30.0, n_rho)
    C_u = np.max(g(rho) / rho ** 2) / (4.0 * math.pi)

    K = g(rho) / rho ** 3
    Kp = dg(rho) / rho ** 3 - 3.0 * g(rho) / rho ** 4  # dK/dr at sigma=1
    eps = np.zeros((3, 3, 3))
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        eps[i, j, k] = 1.0
        eps[i, k, j] = -1.0
    Gam = np.array([0.0, 0.0, 1.0])
    M1 = np.einsum("ijk,j->ik", eps, Gam)
    best = 0.0
    for th in np.linspace(0.0, math.pi, n_theta):
        rhat = np.array([math.sin(th), 0.0, math.cos(th)])
        A = np.cross(Gam, rhat)
        for m in range(len(rho)):
            r = rho[m]
            grad = M1 * K[m] + np.outer(A * r, rhat * r) * (Kp[m] / r)
            best = max(best, np.linalg.norm(grad))
    return C_u, best / (4.0 * math.pi)


def shape_constants():
    """Vectorized C_u, C_g (same math as shape_constants_explicit)."""
    rho = np.linspace(1e-6, 30.0, 200001)
    C_u = np.max(g(rho) / rho ** 2) / (4.0 * math.pi)

    K = g(rho) / rho ** 3
    Kp = dg(rho) / rho ** 3 - 3.0 * g(rho) / rho ** 4
    th = np.linspace(0.0, math.pi, 721)
    st, ct = np.sin(th), np.cos(th)
    # Gamma = e_z, rhat = (st, 0, ct)  =>  Gamma x rhat = (0, st, 0)
    # M1_ik = eps_ijk Gamma_j = [[0,1,0],[-1,0,0],[0,0,0]]
    # nonzero grad entries: (0,1) = K, (1,0) = -K + o10, (1,2) = o12
    RK = rho[:, None] * Kp[:, None]
    o10 = st[None, :] ** 2 * RK
    o12 = (st * ct)[None, :] * RK
    Kb = K[:, None]
    f2 = Kb ** 2 + (o10 - Kb) ** 2 + o12 ** 2
    return C_u, math.sqrt(np.max(f2)) / (4.0 * math.pi)


C_u, C_g = shape_constants()
C_u_ref, C_g_ref = shape_constants_explicit()
assert abs(C_u - C_u_ref) / C_u_ref < 1e-3, (C_u, C_u_ref)
assert abs(C_g - C_g_ref) / C_g_ref < 5e-3, (C_g, C_g_ref)
print(f"C_u = {C_u:.6e}   C_g = {C_g:.6e}")


def write(name, rows, header):
    with open(os.path.join(DATA, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


sig = np.logspace(-4, -1, 400)
write("u_unit.csv", [(s, C_u / s ** 2) for s in sig], ["sigma", "umax"])
write("gradu_unit.csv", [(s, C_g / s ** 3) for s in sig], ["sigma", "gradmax"])
write("constants.csv", [("C_u", C_u), ("C_g", C_g)], ["name", "value"])
