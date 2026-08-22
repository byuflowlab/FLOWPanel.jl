#!/usr/bin/env python3
"""BRAINSTORM 020 Phase 2 Leg 1: the strain ceiling, measured vs physics.

Scores the pre-registered sec.7.1 thresholds (phase_01_theory.md):
  "stability-set" iff p >= 1.5 over the ladder AND Gamma_implied in
  [0.5, 2] Gamma_v for >= 3 consecutive rungs.
Primary numbers come from the 019 pipeline outputs (survivor-side direct-M
fit); a supplementary PRE-ONSET direct M is computed here for ignited runs
(019 P1 sanctions "the last pre-ignition M reading"), because whole-run max
dtZ of an ignited run is a divergence transient (019 wave-1 deviation note).
Onset = first step with max_u > 1000 m/s (tripwire convention).

Physics target: Z_phys = Gamma_v/(2 pi r_c^2), r_c/c in {0.01, 0.05}
(literature tip-vortex core bracket), c = c(0.75R) = 0.0873R/0.68 (018
sigma-ladder table anchor: sigma/R=0.0873 <-> 0.68c). Gap in decades at the
measured ceiling Z_res = M_boundary/dt.

Outputs: BRAINSTORM/020_sigma_aware_subgrid_closure/figures/fig_leg1/*.csv
"""

import csv
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "BRAINSTORM", "020_sigma_aware_subgrid_closure",
                   "figures", "fig_leg1")

DT = 60.0 / 5400.0 / 36.0
GV = 0.277922                      # m^2/s (019 pipeline gamma_v)
R = 0.1190                         # m
C075 = 0.0873 * R / 0.68           # chord at 0.75R from 018 ladder anchor
U_TRIP = 1000.0                    # m/s ignition onset tripwire

# Campaign E screens with the direct max_dtZ column (sigma/R from 019 margins)
RUNS = [
    ("scr_p019_s015",  0.01496, False, "ignited"),
    ("scr_p019_s015v", 0.01496, True,  "ignited"),
    ("scr_p019_s020v", 0.01995, True,  "ignited"),
    ("scr_p019_s025",  0.02493, False, "ignited"),
    ("scr_p019_s025v", 0.02493, True,  "ignited"),
    ("scr_p019_s030v", 0.02992, True,  "survivor"),
    ("scr_p019_sstab", 0.03117, True,  "survivor"),
    ("scr_p019_s038",  0.03808, False, "survivor"),
    ("scr_p019_s038v", 0.03808, True,  "survivor"),
]


def wake_health(run):
    path = os.path.join(ROOT, "data", run, "monitors",
                        f"{run}_monitor04_wake_health_system1.csv")
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def pre_onset_max_dtz(rows):
    onset = None
    for r in rows:
        try:
            if float(r["max_u"]) > U_TRIP:
                onset = int(float(r["step"]))
                break
        except ValueError:
            continue
    m = 0.0
    for r in rows:
        try:
            step = int(float(r["step"]))
            v = float(r["max_dtZ"])
        except (ValueError, KeyError):
            continue
        if onset is not None and step >= onset:
            break
        if math.isfinite(v):
            m = max(m, v)
    return onset, m


def fit_p(points):
    """log-log LS fit M = A sigma^-p over (sigma_over_R, M)."""
    xs = [math.log(s) for s, m in points]
    ys = [math.log(m) for s, m in points]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    return -slope


def main():
    os.makedirs(OUT, exist_ok=True)
    rows_out = []
    print("== Leg 1: margins (direct max_dtZ; ignited runs also PRE-ONSET) ==")
    print(f"{'run':22s} {'s/R':>8s} {'visc':>5s} {'outcome':>9s} "
          f"{'onset':>6s} {'M_pre_onset':>12s} {'M_whole_run':>12s} "
          f"{'G_impl_pre':>11s}")
    for run, sor, visc, outcome in RUNS:
        rows = wake_health(run)
        onset, m_pre = pre_onset_max_dtz(rows)
        m_all = max(float(r["max_dtZ"]) for r in rows
                    if r.get("max_dtZ") not in (None, "", "NaN"))
        sigma = sor * R
        g_impl = 2 * math.pi * sigma * sigma * m_pre / DT if m_pre > 0 else float("nan")
        rows_out.append(dict(run=run, sigma_over_R=sor, viscous=visc,
                             outcome=outcome, onset_step=onset,
                             M_pre_onset=m_pre, M_whole_run=m_all,
                             gamma_implied_pre=g_impl))
        print(f"{run:22s} {sor:8.5f} {str(visc):>5s} {outcome:>9s} "
              f"{str(onset) if onset else '-':>6s} {m_pre:12.4e} {m_all:12.4e} "
              f"{g_impl:11.4f}")

    with open(os.path.join(OUT, "leg1_margins.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        [w.writerow(r) for r in rows_out]

    # figure-facing subsets (fig_leg1.tex reads these directly)
    subsets = {
        "leg1_margins_visc_survivor.csv":
            [r for r in rows_out if r["viscous"] and r["outcome"] == "survivor"],
        "leg1_margins_visc_ignited.csv":
            [r for r in rows_out if r["viscous"] and r["outcome"] == "ignited"],
        "leg1_margins_inviscid.csv":
            [r for r in rows_out if not r["viscous"]],
    }
    for name, sub in subsets.items():
        with open(os.path.join(OUT, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
            w.writeheader()
            [w.writerow(r) for r in sub]

    # p-fits: survivors (pipeline-equivalent) and pre-onset all-rungs (suppl.)
    for label, pts in [
        ("viscous survivors (pipeline)",
         [(r["sigma_over_R"], r["M_pre_onset"]) for r in rows_out
          if r["viscous"] and r["outcome"] == "survivor"]),
        ("viscous all rungs, pre-onset M (supplementary)",
         [(r["sigma_over_R"], r["M_pre_onset"]) for r in rows_out
          if r["viscous"] and r["M_pre_onset"] > 0]),
        ("inviscid all rungs, pre-onset M (supplementary)",
         [(r["sigma_over_R"], r["M_pre_onset"]) for r in rows_out
          if not r["viscous"] and r["M_pre_onset"] > 0]),
    ]:
        if len(pts) >= 2:
            print(f"p-fit [{label}]: p = {fit_p(pts):.3f} over {len(pts)} pts")

    # Gamma_implied band check (pre-registered [0.5, 2] Gamma_v)
    lo, hi = 0.5 * GV, 2.0 * GV
    in_band = [(r["run"], r["gamma_implied_pre"]) for r in rows_out
               if r["viscous"] and lo <= r["gamma_implied_pre"] <= hi]
    print(f"\nGamma_v = {GV:.4f}; band [{lo:.4f}, {hi:.4f}] m^2/s; "
          f"viscous rungs in band (pre-onset): {len(in_band)} -> {in_band}")

    # Ceiling vs physics target
    m_boundary = max(r["M_pre_onset"] for r in rows_out
                     if r["viscous"] and r["outcome"] == "survivor")
    z_res = m_boundary / DT
    print(f"\n== ceiling vs physics ==")
    print(f"resolved ceiling (last viscous survivor): M = {m_boundary:.3f} "
          f"-> Z_res = {z_res:.0f} 1/s")
    gap_rows = []
    for rcfrac in (0.05, 0.01):
        rc = rcfrac * C075
        z_phys = GV / (2 * math.pi * rc * rc)
        dec = math.log10(z_phys / z_res)
        gap_rows.append(dict(rc_over_c=rcfrac, rc_m=rc, Z_phys=z_phys,
                             dtZ_phys=DT * z_phys, decades_gap=dec))
        print(f"r_c = {rcfrac:.2f}c = {rc*1e3:.3f} mm: Z_phys = {z_phys:.3e} 1/s "
              f"(dtZ = {DT*z_phys:.1f}), gap = {dec:.2f} decades")
    with open(os.path.join(OUT, "leg1_gap.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(gap_rows[0].keys()))
        w.writeheader()
        [w.writerow(r) for r in gap_rows]
    print(f"(chord anchor c(0.75R) = {C075*1e3:.2f} mm from 018 sigma-ladder "
          f"table: 0.0873R <-> 0.68c)")


if __name__ == "__main__":
    main()
