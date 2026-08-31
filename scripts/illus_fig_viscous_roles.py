#!/usr/bin/env python3
"""F2a/F2b fig_viscous_roles_{a,b} data for the sigma/VPM illustration
campaign (plans/sigma_vpm_illustrations_20260827). Shared data dir
fig_viscous_roles/ serves both figures.

F2a (broken, inviscid): strain thinning is unopposed. Continuous
sigma(t) = sigma0 exp(-Z t) stays positive; the explicit-Euler iterates
sigma_{n+1} = sigma_n (1 - dtZ) cross zero for dtZ > 1 and oscillate negative.

F2b (fixed, viscous): d(sigma^2)/dt = -2 Z sigma^2 + 2 nu  gives
sigma^2(t) = nu/Z + (sigma0^2 - nu/Z) exp(-2 Z t): every start converges to
the strain-diffusion equilibrium sigma_eq = sqrt(nu/Z).
"""
import csv
import math
import os

OUT = os.path.expanduser(
    "~/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations")
DATA = os.path.join(OUT, "fig_viscous_roles")
os.makedirs(DATA, exist_ok=True)


def write(name, rows, header):
    with open(os.path.join(DATA, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


# ---- F2a: inviscid, time in units of steps; sigma in units of sigma0 -------
NSTEP = 12
for tag, dtZ in [("mild", 0.5), ("severe", 1.5)]:
    rows_e, rows_c = [], []
    s = 1.0
    for n in range(NSTEP + 1):
        rows_e.append((n, s))
        s *= (1.0 - dtZ)
    for i in range(200):
        t = NSTEP * i / 199
        rows_c.append((t, math.exp(-dtZ * t)))
    write(f"inviscid_euler_{tag}.csv", rows_e, ["step", "sigma_rel"])
    write(f"inviscid_exact_{tag}.csv", rows_c, ["step", "sigma_rel"])

# ---- F2b: viscous, time in units of 1/(2Z); sigma in units of sigma_eq -----
for tag, s0 in [("above", 2.0), ("below", 0.4)]:
    rows = []
    for i in range(300):
        t = 4.0 * i / 299        # t in units of 1/(2Z)
        s2 = 1.0 + (s0 ** 2 - 1.0) * math.exp(-t)
        rows.append((t, math.sqrt(max(s2, 0.0))))
    write(f"viscous_{tag}.csv", rows, ["tau", "sigma_over_eq"])

# pure diffusion reference (no strain): sigma^2 = s0^2 + (nu-like) growth,
# in the same units: s2 = s0^2 + t/2  (arbitrary but monotone; labeled as such)
rows = [(4.0 * i / 299, math.sqrt(0.4 ** 2 + (4.0 * i / 299) / 2)) for i in range(300)]
write("viscous_diffusion_only.csv", rows, ["tau", "sigma_over_eq"])

print("wrote", DATA)
