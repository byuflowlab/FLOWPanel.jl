#!/usr/bin/env python3
"""BRAINSTORM/022 harvest: CT cycle statistics, IGE/OGE ratios, figure data.

Implements the M1 metric defined in
BRAINSTORM/022_rotor_hover_ground_effect/decision_rules.md: CT = negated axial
Bernoulli force channel, cycle-mean +/- cycle-std over a settled rev window,
where cycle-std is the n-1 standard deviation *across per-rev block means*
(matching the driver's own definition in examples/rotor_hover_ground_effect.jl).

Why this is not scripts/p018_analyze.py: that script pins RPM_DEFAULT = 5400 in
its force-monitor fallback with no CLI pass-through, and 022 runs at RPM 6000.
The fallback is the load-bearing path here -- `<run>_CT_vs_rev.csv` is only
written after the time loop finishes, and 022's coarse runs are hitting the
walltime or being cancelled. Shared helpers are imported from p018_analyze
rather than duplicated.

Usage
  p022_harvest.py m1 RUN [--revs LO HI] [--rpm 6000]
  p022_harvest.py ratio IGE_RUN OGE_RUN [--revs LO HI] [--rpm 6000]
  p022_harvest.py figdata [--out DIR] [--rpm 6000]
                          [--point H_OVER_R IGE_RUN OGE_RUN LO HI] ...

Rev windows are inclusive-exclusive integer revolution indices [LO, HI), matching
the `rev_block` convention in the *_CT_per_rev.csv files.
"""

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p018_analyze import DATA, PLACEHOLDER_TOL, _mean, _std, per_rev  # noqa: E402

# 022 operating point (item file, standing ruling 1).
RPM_022 = 6000.0

# Freestream is fully withdrawn at ramp + hold + withdraw = 1 + 1.5 + 4 revs;
# an M1 window is only valid if it lies entirely after this.
HOVER_START_REV = 6.5

# SETTLE_REVS = 20 in the launcher's shared carrier block, so the intended M1
# window is revs 18-28. A window that merely clears HOVER_START_REV is still
# inside the settling transient and is reported as preliminary.
SETTLED_START_REV = 18

# |CT| above this is a solver blow-up, not a physical sample. Flagged loudly and
# excluded -- a single such step (p022_ige_coarse hit CT ~ 372 at step 636) would
# otherwise swamp a cycle mean.
BLOWUP_CT = 1.0

# Experiment, provided by Ryan 2026-08-19. Source citation TBD.
# h/R -> (T/T_OGE, +/- uncertainty)
EXPERIMENT = [
    (0.5, 1.200, 0.009),
    (1.0, 1.078, 0.008),
    (1.5, 1.010, 0.014),
    (2.0, 1.004, 0.004),
    (2.5, 0.996, 0.009),
    (3.0, 0.999, 0.023),
    (4.0, 1.009, 0.008),
    (5.5, 0.998, 0.012),
    (7.0, 1.001, 0.012),
]


def cheeseman_bennett(h_over_R):
    """Momentum-theory IGE thrust augmentation at fixed power, T/T_OGE."""
    return 1.0 / (1.0 - (1.0 / (4.0 * h_over_R)) ** 2)


def _ct_path(run):
    return os.path.join(DATA, run, f"{run}_CT_vs_rev.csv")


def _force_monitor_path(run):
    return os.path.join(DATA, run, "monitors",
                        f"{run}_monitor02_force_system1.csv")


def load_ct(run, rpm=RPM_022):
    """Return ([(step, revolution, CT)], [blowup rows]) for a 022 run.

    Prefers the end-of-run CT CSV; falls back to the incrementally-written force
    monitor (CT = -CFx, revolution = time*rpm/60, step_csv = step_monitor + 1),
    which is the only source for a run that did not reach the end of its loop.
    """
    rows, blowups = [], []

    if os.path.exists(_ct_path(run)):
        with open(_ct_path(run)) as fh:
            for rec in csv.DictReader(fh):
                try:
                    ct = float(rec["CT_bernoulli"])
                    step, rev = int(rec["step"]), float(rec["revolution"])
                except (TypeError, ValueError, KeyError):
                    continue
                if not math.isfinite(ct) or abs(ct) < PLACEHOLDER_TOL:
                    continue
                (blowups if abs(ct) > BLOWUP_CT else rows).append((step, rev, ct))
        return rows, blowups

    if not os.path.exists(_force_monitor_path(run)):
        raise SystemExit(f"no CT source for run {run!r} under {DATA}")

    sys.stderr.write(
        f"note: {run} has no CT_vs_rev.csv (run did not finish); "
        f"reconstructing CT from the force monitor at RPM {rpm:g}\n")
    with open(_force_monitor_path(run)) as fh:
        for rec in csv.DictReader(fh):
            try:
                ct = -float(rec["CFx"])
                t = float(rec["time"])
                step = int(rec["step"]) + 1
            except (TypeError, ValueError, KeyError):
                continue
            if not math.isfinite(ct) or abs(ct) < PLACEHOLDER_TOL:
                continue
            (blowups if abs(ct) > BLOWUP_CT else rows).append(
                (step, t * rpm / 60.0, ct))
    return rows, blowups


def cycle_stats(rows, lo, hi):
    """Cycle-mean +/- cycle-std over integer rev blocks in [lo, hi).

    cycle_std is the n-1 std across per-rev block means (the driver's
    definition), NOT across individual samples. Also returns the standard error
    of that mean, which is what an IGE/OGE ratio comparison actually needs.
    """
    bins = per_rev(rows)
    revs = sorted(r for r in bins if lo <= r < hi)
    if not revs:
        raise SystemExit(f"no samples in rev window [{lo}, {hi})")
    block_means = [_mean(bins[r]) for r in revs]
    mean = _mean(block_means)
    std = _std(block_means)
    return {
        "revs": revs,
        "block_means": block_means,
        "n_blocks": len(revs),
        "missing": sorted(set(range(lo, hi)) - set(revs)),
        "mean": mean,
        "std": std,
        "rel_std": std / abs(mean) if mean else float("nan"),
        "sem": std / math.sqrt(len(revs)) if len(revs) > 1 else float("nan"),
        "in_hover": lo >= HOVER_START_REV,
        "settled": lo >= SETTLED_START_REV,
    }


def _report_blowups(run, blowups):
    if not blowups:
        return
    sys.stderr.write(
        f"WARNING: {run} has {len(blowups)} step(s) with |CT| > {BLOWUP_CT} "
        f"(solver blow-up), excluded from all statistics:\n")
    for step, rev, ct in blowups[:5]:
        sys.stderr.write(f"    step {step} (rev {rev:.2f}): CT = {ct:.4g}\n")
    first = min(b[0] for b in blowups)
    sys.stderr.write(
        f"    first bad step is {first}; any window reaching it is NOT harvestable\n")


def _print_m1(run, st, blowups):
    print(f"{run}: revs [{st['revs'][0]}, {st['revs'][-1]}] "
          f"({st['n_blocks']} blocks)")
    for r, bm in zip(st["revs"], st["block_means"]):
        print(f"    rev {r:3d}: {bm:.6f}")
    print(f"  cycle-mean  = {st['mean']:.6f}")
    print(f"  cycle-std   = {st['std']:.6f}  ({100 * st['rel_std']:.2f}%)")
    print(f"  SE of mean  = {st['sem']:.6f}  ({100 * st['sem'] / abs(st['mean']):.2f}%)")
    print(f"  window_in_hover = {st['in_hover']} (hover starts at rev {HOVER_START_REV})")
    if not st["in_hover"]:
        print("  NOTE: window starts before freestream withdrawal completes -- "
              "PRELIMINARY, not an M1 harvest")
    elif not st["settled"]:
        print(f"  NOTE: window starts before rev {SETTLED_START_REV} "
              f"(SETTLE_REVS=20) -- inside the settling transient, PRELIMINARY")
    if st["missing"]:
        print(f"  WARNING: requested window is only partly covered -- "
              f"{len(st['missing'])} rev(s) absent from the data "
              f"({st['missing'][0]}..{st['missing'][-1]}); the run did not get "
              f"this far. Statistics are over {st['n_blocks']} block(s) only.")
    if st["n_blocks"] < 2:
        print("  WARNING: fewer than 2 rev blocks -- cycle-std is meaningless")
    if blowups:
        print(f"  NOTE: {len(blowups)} blow-up step(s) excluded; see stderr")


def cmd_m1(args):
    rows, blowups = load_ct(args.run, args.rpm)
    _report_blowups(args.run, blowups)
    st = cycle_stats(rows, *args.revs)
    _print_m1(args.run, st, blowups)


def ratio_stats(ige, oge, lo, hi, rpm=RPM_022):
    """Matched IGE/OGE thrust ratio with quadrature-propagated uncertainty."""
    out = {}
    for key, run in (("ige", ige), ("oge", oge)):
        rows, blowups = load_ct(run, rpm)
        _report_blowups(run, blowups)
        out[key] = cycle_stats(rows, lo, hi)
        out[key + "_blowups"] = blowups
    ratio = out["ige"]["mean"] / out["oge"]["mean"]
    rel = lambda k, f: out[k][f] / abs(out[k]["mean"])  # noqa: E731
    out["ratio"] = ratio
    out["ratio_std"] = ratio * math.hypot(rel("ige", "std"), rel("oge", "std"))
    out["ratio_sem"] = ratio * math.hypot(rel("ige", "sem"), rel("oge", "sem"))
    return out


def cmd_ratio(args):
    lo, hi = args.revs
    st = ratio_stats(args.ige, args.oge, lo, hi, args.rpm)
    for key, run in (("ige", args.ige), ("oge", args.oge)):
        _print_m1(run, st[key], st[key + "_blowups"])
        print()
    print(f"ratio IGE/OGE over revs [{lo}, {hi}) = {st['ratio']:.4f}")
    print(f"  +/- {st['ratio_std']:.4f}  (cycle-std propagated, per-rev spread)")
    print(f"  +/- {st['ratio_sem']:.4f}  (SE of mean propagated, uncertainty on the mean)")
    if not (st["ige"]["in_hover"] and st["oge"]["in_hover"]):
        print("  PRELIMINARY: window is not entirely post-freestream-withdrawal")


def _write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")


def cmd_figdata(args):
    os.makedirs(args.out, exist_ok=True)

    _write_csv(os.path.join(args.out, "experiment.csv"),
               ["h_over_R", "ratio", "eminus", "eplus"],
               [(h, v, e, e) for h, v, e in EXPERIMENT])

    # Fine grid for a smooth theory curve; CB diverges as h -> R/4 so start well
    # above it. 0.4R is already below the lowest experimental point.
    grid = [0.4 + 0.01 * i for i in range(int((2.2 - 0.4) / 0.01) + 1)]
    _write_csv(os.path.join(args.out, "theory_cb.csv"),
               ["h_over_R", "ratio"],
               [(f"{h:.2f}", f"{cheeseman_bennett(h):.6f}") for h in grid])

    sim_rows = []
    for h_over_R, ige, oge, lo, hi in args.point:
        st = ratio_stats(ige, oge, int(lo), int(hi), args.rpm)
        flag = "settled" if (st["ige"]["settled"] and st["oge"]["settled"]) \
            else "preliminary"
        sim_rows.append((h_over_R, f"{st['ige']['mean']:.6f}",
                         f"{st['ratio']:.4f}", f"{st['ratio_std']:.4f}",
                         f"{st['ratio_std']:.4f}", f"{lo}-{hi}", flag))
    _write_csv(os.path.join(args.out, "sim.csv"),
               ["h_over_R", "ct_ige", "ratio", "eminus", "eplus", "window", "status"],
               sim_rows)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rpm", type=float, default=RPM_022,
                   help="rotor speed for the force-monitor fallback (default 6000)")
    sub = p.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("m1", help="cycle statistics for one run")
    m.add_argument("run")
    m.add_argument("--revs", nargs=2, type=int, default=(18, 28),
                   metavar=("LO", "HI"))
    m.set_defaults(func=cmd_m1)

    r = sub.add_parser("ratio", help="matched IGE/OGE thrust ratio")
    r.add_argument("ige")
    r.add_argument("oge")
    r.add_argument("--revs", nargs=2, type=int, default=(18, 28),
                   metavar=("LO", "HI"))
    r.set_defaults(func=cmd_ratio)

    f = sub.add_parser("figdata", help="emit CSVs backing fig_hr_sweep.tex")
    f.add_argument("--out", default=os.path.join(
        os.path.dirname(DATA), "BRAINSTORM", "022_rotor_hover_ground_effect",
        "figures", "fig_hr_sweep"))
    f.add_argument("--point", nargs=5, action="append", default=[],
                   metavar=("H_OVER_R", "IGE_RUN", "OGE_RUN", "LO", "HI"),
                   help="add a simulation point; repeatable")
    f.set_defaults(func=cmd_figdata)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
