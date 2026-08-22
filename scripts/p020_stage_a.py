#!/usr/bin/env python3
"""BRAINSTORM 020 Phase 2 Leg 2 Stage A: offline re-integration of recorded
collapse data under the exact (exponential) local map vs the live forward
Euler map.

Pre-registration (phase_01_theory.md sec.7.2) assumed tracked per-particle
trajectories; none were retained (recorded deviation, item Log). The
reinterpretation, faithful to prediction P-2:

  A1. Pointwise one-step census over every recorded per-particle dtZ
      (corpse step-1041 strained subset, whole recorded field, healthy
      step-719 control): count spurious-flip exceedances of the Euler map
      and compare multiplier magnitudes against the exact map (which has
      none by construction).
  A2. Trace re-integration: infer per-step dtZ from recorded sigma traces
      (ignition core + argmin-sigma trace of the viscous corpse; min_sr
      trace of the inviscid ufront_n1 twin), re-integrate both maps from
      the first row, overlay. Exact map: no negative sigma, tracks the
      recorded contraction into the floor (fidelity failure persists).

Convention: forensics CSVs record RAW dtZ (f=0,g=1 stretching projection);
the live update uses Z_MM4 = raw/5. Spurious Euler thresholds in MM4 units:
transverse-Gamma flip 2/3, sigma flip 1, sigma geometric 2.

Outputs (deterministic, byte-identical on rerun):
  BRAINSTORM/020_sigma_aware_subgrid_closure/figures/fig_stage_a/*.csv
  stdout census + trace tables (evidence-pack source of record).
"""

import csv
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FOR = os.path.join(ROOT, "data", "p018_L2_visc_forensics")
OUT = os.path.join(ROOT, "BRAINSTORM", "020_sigma_aware_subgrid_closure",
                   "figures", "fig_stage_a")

DT = 3.0864197530864197e-04      # s (60/5400/36), 019 pipeline convention
NU = 1.433418e-05                # m^2/s, metadata wake_nu (authoritative)
RAW_TO_MM4 = 1.0 / 5.0           # live Z = raw/5 (f=0, g=1/5)

# MM4-unit thresholds (phase_01_theory.md sec.1.3 table)
TH_GPERP = 2.0 / 3.0             # transverse-Gamma sign flip
TH_SIG = 1.0                     # negative sigma (inviscid)
TH_SIG2 = 2.0                    # sigma geometric divergence


def read_csv(path):
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        return [row for row in r]


def census(rows, key, label, scale):
    """One-step Euler-vs-exact exceedance census over per-particle dtZ."""
    xs = []
    for row in rows:
        try:
            v = float(row[key]) * scale
        except ValueError:
            continue
        if math.isfinite(v):
            xs.append(v)
    n = len(xs)
    n_gperp = sum(1 for x in xs if x > TH_GPERP)
    n_sig = sum(1 for x in xs if x > TH_SIG)
    n_sig2 = sum(1 for x in xs if x > TH_SIG2)
    # worst Euler multiplier magnitudes vs exact (exact are e^{-3x}, e^{-x} <= e^{|3x|})
    worst = max(xs) if xs else float("nan")
    mult_gperp_euler = abs(1 - 3 * worst)
    mult_sig_euler = abs(1 - worst)
    return dict(label=label, n=n, max_dtZ_mm4=worst,
                n_flip_gperp=n_gperp, n_flip_sigma=n_sig, n_geom_sigma=n_sig2,
                worst_euler_mult_gperp=mult_gperp_euler,
                worst_euler_mult_sigma=mult_sig_euler,
                worst_exact_mult_gperp=math.exp(-3 * worst) if worst < 200 else 0.0,
                worst_exact_mult_sigma=math.exp(-worst) if worst < 200 else 0.0)


def infer_dtz_viscous(sig0, sig1):
    """Invert sigma1 = sqrt((sig0*(1-dtZ))^2 + 2 nu dt) for dtZ (Euler+CS)."""
    s2 = sig1 * sig1 - 2 * NU * DT
    if s2 < 0:
        s2 = 0.0
    return 1.0 - math.sqrt(s2) / sig0


def reintegrate(sig_start, dtzs, viscous, mode):
    """Forward re-integration of a sigma trace under euler|exp composed with
    CoreSpreading when viscous."""
    sig = sig_start
    out = [sig]
    for x in dtzs:
        if mode == "euler":
            sig = sig * (1.0 - x)
        else:
            sig = sig * math.exp(-x)
        if viscous:
            sig = math.sqrt(sig * sig + 2 * NU * DT)
        out.append(sig)
    return out


def main():
    os.makedirs(OUT, exist_ok=True)

    # ---------- A1: pointwise census ---------------------------------------
    tables = []
    tables.append(census(read_csv(os.path.join(FOR, "strain_corpse1041.csv")),
                         "dtZ_all", "corpse1041_collapsing_subset(864)", RAW_TO_MM4))
    tables.append(census(read_csv(os.path.join(FOR, "recordedJ_dtZ_corpse1041.csv")),
                         "dtZ", "corpse1041_whole_field(177k)", RAW_TO_MM4))
    tables.append(census(read_csv(os.path.join(FOR, "strain_healthy719_gatea.csv")),
                         "dtZ", "healthy719_aged_column_control", RAW_TO_MM4))
    tables.append(census(read_csv(os.path.join(FOR, "recordedJ_dtZ_healthy719.csv")),
                         "dtZ", "healthy719_whole_field_control", RAW_TO_MM4))

    with open(os.path.join(OUT, "census.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(tables[0].keys()))
        w.writeheader()
        for t in tables:
            w.writerow(t)

    print("== A1 census (MM4 units; raw forensics dtZ / 5) ==")
    print(f"{'population':40s} {'N':>7s} {'max dtZ':>9s} "
          f"{'>2/3':>6s} {'>1':>6s} {'>2':>6s}")
    for t in tables:
        print(f"{t['label']:40s} {t['n']:7d} {t['max_dtZ_mm4']:9.3f} "
              f"{t['n_flip_gperp']:6d} {t['n_flip_sigma']:6d} {t['n_geom_sigma']:6d}")
    print("exact map: 0 flips at any dtZ by construction "
          "(multipliers e^-3x, e^-x > 0)")

    # ---------- A2: trace re-integrations ----------------------------------
    # (a) viscous corpse, ignition-core trace (max Gos2 particle)
    core = read_csv(os.path.join(FOR, "ignition_core.csv"))
    sigs = [float(r["sigma"]) for r in core]
    steps = [int(float(r["step"])) for r in core]
    dtzs = [infer_dtz_viscous(sigs[i], sigs[i + 1]) for i in range(len(sigs) - 1)]
    eul = reintegrate(sigs[0], dtzs, True, "euler")
    exp = reintegrate(sigs[0], dtzs, True, "exp")
    with open(os.path.join(OUT, "trace_ignition_core.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "sigma_recorded", "dtZ_inferred_mm4",
                    "sigma_reint_euler", "sigma_reint_exp"])
        for i, s in enumerate(steps):
            w.writerow([s, f"{sigs[i]:.6e}",
                        f"{dtzs[i - 1]:.6e}" if i > 0 else "",
                        f"{eul[i]:.6e}", f"{exp[i]:.6e}"])

    # (b) viscous corpse, argmin-sigma trace (particle-switch caveat labeled)
    death = read_csv(os.path.join(FOR, "death_trajectory.csv"))
    dsig = [float(r["min_sigma"]) for r in death]
    dsteps = [int(float(r["step"])) for r in death]
    ddtz = [infer_dtz_viscous(dsig[i], dsig[i + 1]) for i in range(len(dsig) - 1)]
    deul = reintegrate(dsig[0], ddtz, True, "euler")
    dexp = reintegrate(dsig[0], ddtz, True, "exp")
    with open(os.path.join(OUT, "trace_argmin_sigma.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "min_sigma_recorded", "dtZ_inferred_mm4",
                    "sigma_reint_euler", "sigma_reint_exp"])
        for i, s in enumerate(dsteps):
            w.writerow([s, f"{dsig[i]:.6e}",
                        f"{ddtz[i - 1]:.6e}" if i > 0 else "",
                        f"{deul[i]:.6e}", f"{dexp[i]:.6e}"])

    # (c) inviscid ufront_n1 min_sigma_ratio trace through its death window
    wh = read_csv(os.path.join(
        ROOT, "data", "p018_ufront_n1", "monitors",
        "p018_ufront_n1_monitor04_wake_health_system1.csv"))
    rows = [(int(float(r["step"])), float(r["min_sigma_ratio"])) for r in wh
            if r["min_sigma_ratio"] not in ("", "NaN")]
    # death window: last 40 recorded steps (ends min_sr < 0)
    rows = rows[-40:]
    srs = [v for _, v in rows]
    usteps = [s for s, _ in rows]
    # inviscid Euler inversion: sr1 = sr0*(1-dtZ)  (valid through the flip)
    udtz = [1.0 - srs[i + 1] / srs[i] if srs[i] > 0 else float("nan")
            for i in range(len(srs) - 1)]
    ueul = [srs[0]]
    uexp = [srs[0]]
    for x in udtz:
        ueul.append(ueul[-1] * (1.0 - x) if math.isfinite(x) else float("nan"))
        uexp.append(uexp[-1] * math.exp(-min(x, 300.0)) if math.isfinite(x)
                    else uexp[-1])
    with open(os.path.join(OUT, "trace_ufront_n1_minsr.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "min_sr_recorded", "dtZ_inferred_mm4",
                    "minsr_reint_euler", "minsr_reint_exp"])
        for i, s in enumerate(usteps):
            w.writerow([s, f"{srs[i]:.6e}",
                        f"{udtz[i - 1]:.6e}" if i > 0 and math.isfinite(udtz[i - 1]) else "",
                        f"{ueul[i]:.6e}" if math.isfinite(ueul[i]) else "nan",
                        f"{uexp[i]:.6e}"])

    # ---------- verdict summary --------------------------------------------
    def neg_events(seq):
        return sum(1 for v in seq if v < 0)

    print("\n== A2 trace re-integration ==")
    print(f"{'trace':28s} {'euler reproduces?':>18s} {'neg-sig euler':>14s} "
          f"{'neg-sig exp':>12s} {'exp floor / recorded floor':>27s}")
    for name, rec, e, x in [
            ("ignition_core (viscous)", sigs, eul, exp),
            ("argmin_sigma (viscous)", dsig, deul, dexp),
            ("ufront_n1 min_sr (invisc)", srs, ueul, uexp)]:
        # Euler re-integration must reproduce the recorded trace (inversion
        # consistency); compare where recorded stayed positive.
        ok = all(abs(e[i] - rec[i]) <= 1e-9 * abs(rec[i]) + 1e-30
                 for i in range(len(rec)) if rec[i] > 0 and math.isfinite(e[i]))
        floor_ratio = min(x) / min(abs(v) for v in rec if v != 0)
        print(f"{name:28s} {str(ok):>18s} {neg_events(e):>14d} "
              f"{neg_events(x):>12d} {floor_ratio:27.3f}")

    print("\nPass criteria (sec.7.2 Stage A): exact map has zero negative-sigma "
          "events and zero\nmultiplier-magnitude excess; and still reaches <= 1.1x "
          "the recorded floor (fidelity\nfailure persists).")


if __name__ == "__main__":
    main()
