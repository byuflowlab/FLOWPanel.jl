#!/usr/bin/env python3
"""Burst/quiet-epoch decomposition of per-step CT (BRAINSTORM/018 phase 15).

REPORT-ONLY analysis layer: decision_rules.md is unchanged; M1 stays the
pre-registered gate. Motivated by the 2026-08-12 finding (ledger, "Phase-15
subtask (a)") that per-rev mean CT is burst-rectified: it tracks the
within-rev std s_k (r ~= 0.61), so the monotone-drift statistic conflates a
slowly-evolving quiet baseline with the episode timing of bursts.

Per arm and rev window, from the same CT series p018_analyze.py scores:
  - per-rev mean CT_k and within-rev std s_k (the burst-amplitude proxy);
  - quiet-epoch CT level two ways: (a) intercept of CT_k = a + b*s_k
    (extrapolation to the burst-free limit), (b) mean over quiet revs
    (s_k < s*), with s* pre-registered as the POOLED lower quartile of s_k
    across all arms passed in one invocation (override with --sstar);
  - burst incidence, episode count/durations;
  - drift of the quiet-rev series: linear slope with AR(1)-corrected CI,
    plus an M1-style first-half/second-half delta;
  - per-rev min(min_sigma_ratio) from the wake-health CSV (when present) and
    its detrended lag-0 correlation with s_k -- the under-damping check;
  - cross-arm summary table (drift-rate comparison across sigma).

Run syntax mirrors p018_analyze stitching, compactly:
  python3 scripts/p018_burst_stats.py RUN[+SEGMENT@RESTART_STEP] [RUN2 ...] \
      [--lo 10] [--hi 44] [--sstar X]

s* POOLING RULE (measured 2026-08-12): pool s* only across arms of ONE
carrier family. Families differ ~5-10x in per-rev scatter (B0 eta-Das N=4 vs
viscous N=1 ufront), so a cross-family pooled quartile classifies the noisier
family as all-burst and the quieter as all-quiet. For cross-family
comparisons use the regression intercept (threshold-free) column.

Validation (2026-08-12): p018_b0 revs 15-29 per-rev grand mean reproduces
0.071701 (p018_analyze M1); p018_ufront_n1_s040_visc(+_s2@1079) revs 10-44
reproduces r(mean, std) ~= 0.61 and quiet-limit intercept ~= 0.0730 from the
phase-15 subtask (a) subagent.
"""

import argparse
import csv
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p018_analyze import DATA, load_ct, stitch, _mean, _std  # noqa: E402


# ----------------------------------------------------------------------------
# loading

def parse_run_spec(spec):
    """'run', or 'parent+segment@step' -> (label, rows)."""
    if "+" in spec:
        parent, rest = spec.split("+", 1)
        segment, restart = rest.split("@", 1)
        return parent, stitch(parent, segment, int(restart))
    return spec, load_ct(spec)


def _wake_health_path(run):
    hits = sorted(glob.glob(os.path.join(
        DATA, run, "monitors", f"{run}_monitor*_wake_health_system1.csv")))
    return hits[-1] if hits else None


def load_min_sr(spec):
    """Per-step (step, min_sigma_ratio) honoring the same stitch as the CT."""
    def one(run, lo=None, hi=None):
        path = _wake_health_path(run)
        if path is None:
            return []
        out = []
        with open(path) as fh:
            for rec in csv.DictReader(fh):
                try:
                    step = int(rec["step"])
                    sr = float(rec["min_sigma_ratio"])
                except (KeyError, TypeError, ValueError):
                    continue
                if not math.isfinite(sr):
                    continue
                if lo is not None and step <= lo:
                    continue
                if hi is not None and step > hi:
                    continue
                out.append((step, sr))
        return out

    if "+" in spec:
        parent, rest = spec.split("+", 1)
        segment, restart = rest.split("@", 1)
        restart = int(restart)
        return sorted(one(parent, hi=restart) + one(segment, lo=restart))
    return one(spec)


# ----------------------------------------------------------------------------
# statistics helpers

def linreg(xs, ys):
    """OLS y = a + b x -> (a, b, r)."""
    n = len(xs)
    mx, my = _mean(xs), _mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0 or n < 3:
        return my, 0.0, float("nan")
    b = sxy / sxx
    return my - b * mx, b, sxy / math.sqrt(sxx * syy)


def detrend(xs, ys):
    a, b, _ = linreg(xs, ys)
    return [y - (a + b * x) for x, y in zip(xs, ys)]


def pearson(xs, ys):
    return linreg(xs, ys)[2]


def lag1_rho(vals):
    if len(vals) < 3:
        return 0.0
    return pearson(vals[:-1], vals[1:])


def trend_with_ar1_ci(revs, vals):
    """Slope per rev, with a 95% CI inflated by AR(1) effective sample size."""
    n = len(vals)
    a, b, _ = linreg(revs, vals)
    resid = [v - (a + b * r) for r, v in zip(revs, vals)]
    rho = max(-0.99, min(0.99, lag1_rho(resid)))
    n_eff = max(3.0, n * (1 - rho) / (1 + rho))
    mx = _mean(revs)
    sxx = sum((r - mx) ** 2 for r in revs)
    if sxx == 0 or n < 4:
        return b, float("nan"), rho
    s2 = sum(x * x for x in resid) / (n - 2)
    se = math.sqrt(s2 / sxx) * math.sqrt(n / n_eff)
    return b, 1.96 * se, rho


def episodes(flags):
    """Lengths of consecutive-True runs."""
    out, cur = [], 0
    for f in flags:
        if f:
            cur += 1
        elif cur:
            out.append(cur)
            cur = 0
    if cur:
        out.append(cur)
    return out


# ----------------------------------------------------------------------------
# per-arm decomposition

def per_rev_stats(rows, lo, hi):
    """[(rev, mean, std, n)] over integer revs in [lo, hi], full revs only."""
    bins = {}
    counts = {}
    for _, rev, ct in rows:
        k = int(math.floor(rev + 1e-9))
        bins.setdefault(k, []).append(ct)
    if not bins:
        return []
    full = max(len(v) for v in bins.values())
    out = []
    for k in sorted(bins):
        if not (lo <= k <= hi):
            continue
        v = bins[k]
        if len(v) < full:            # partial first/last rev: skip
            continue
        out.append((k, _mean(v), _std(v), len(v)))
    return out


def analyze_arm(spec, lo, hi):
    label, rows = parse_run_spec(spec)
    prs = per_rev_stats(rows, lo, hi)
    if len(prs) < 5:
        raise SystemExit(f"{label}: <5 full revs in [{lo}, {hi}]")
    revs = [p[0] for p in prs]
    means = [p[1] for p in prs]
    stds = [p[2] for p in prs]

    # per-rev min_sr aligned on rev via the CT rows' step->rev map
    step_rev = {r[0]: r[1] for r in rows}
    sr_bins = {}
    for step, sr in load_min_sr(spec):
        if step in step_rev:
            k = int(math.floor(step_rev[step] + 1e-9))
            sr_bins.setdefault(k, []).append(sr)
    min_sr = [min(sr_bins[k]) if k in sr_bins else float("nan") for k in revs]

    a, b, r = linreg(stds, means)
    return {
        "label": label, "revs": revs, "means": means, "stds": stds,
        "min_sr": min_sr, "intercept": a, "slope": b, "r_mean_std": r,
    }


def quiet_stats(arm, sstar):
    revs, means, stds = arm["revs"], arm["means"], arm["stds"]
    quiet = [s < sstar for s in stds]
    qrevs = [r for r, q in zip(revs, quiet) if q]
    qmeans = [m for m, q in zip(means, quiet) if q]
    out = {
        "n_quiet": len(qmeans),
        "incidence": 1.0 - len(qmeans) / len(means),
        "episodes": episodes([not q for q in quiet]),
        "quiet_mean": _mean(qmeans) if qmeans else float("nan"),
        "trend": (float("nan"), float("nan"), float("nan")),
        "half_delta_pct": float("nan"),
    }
    if len(qmeans) >= 5:
        out["trend"] = trend_with_ar1_ci(qrevs, qmeans)
        h = len(qmeans) // 2
        m1_, m2_ = _mean(qmeans[:h]), _mean(qmeans[h:])
        out["half_delta_pct"] = 100.0 * (m2_ - m1_) / m1_
    return out


# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("runs", nargs="+",
                    help="RUN or RUN+SEGMENT@RESTART_STEP")
    ap.add_argument("--lo", type=int, default=10)
    ap.add_argument("--hi", type=int, default=10 ** 6)
    ap.add_argument("--sstar", type=float, default=None,
                    help="quiet threshold on within-rev std "
                         "(default: pooled lower quartile across arms)")
    args = ap.parse_args()

    arms = [analyze_arm(spec, args.lo, args.hi) for spec in args.runs]

    if args.sstar is None:
        pooled = sorted(s for arm in arms for s in arm["stds"])
        h = (len(pooled) - 1) * 0.25
        i = int(h)
        sstar = pooled[i] + (h - i) * (pooled[min(i + 1, len(pooled) - 1)]
                                       - pooled[i])
        src = f"pooled p25 over {len(arms)} arm(s)"
    else:
        sstar, src = args.sstar, "user"
    print(f"s* (quiet threshold on within-rev std) = {sstar:.6f}  [{src}]")
    print(f"window: revs [{args.lo}, {args.hi}], full revs only\n")

    for arm in arms:
        q = quiet_stats(arm, sstar)
        grand = _mean(arm["means"])
        slope, ci, rho = q["trend"]
        drift_pct = (100.0 * slope * (arm["revs"][-1] - arm["revs"][0])
                     / q["quiet_mean"]) if q["n_quiet"] >= 5 else float("nan")
        ci_pct = (100.0 * ci * (arm["revs"][-1] - arm["revs"][0])
                  / q["quiet_mean"]) if q["n_quiet"] >= 5 else float("nan")

        # s_k vs per-rev min_sr, both detrended (under-damping check)
        have = [(r, s, m) for r, s, m in
                zip(arm["revs"], arm["stds"], arm["min_sr"])
                if math.isfinite(m)]
        if len(have) >= 5:
            hr = [h[0] for h in have]
            ds = detrend(hr, [h[1] for h in have])
            dm = detrend(hr, [h[2] for h in have])
            r_sr = pearson(ds, dm)
            sr_note = f"{r_sr:+.3f} (n={len(have)})"
        else:
            sr_note = "n/a (no wake-health CSV)"

        print(f"== {arm['label']}  (revs {arm['revs'][0]}-{arm['revs'][-1]}, "
              f"n={len(arm['revs'])})")
        print(f"   grand per-rev mean CT      : {grand:.6f}")
        print(f"   r(mean_k, std_k)           : {arm['r_mean_std']:+.3f}   "
              f"[burst rectification]")
        print(f"   quiet limit (regression a) : {arm['intercept']:.6f}   "
              f"(slope {arm['slope']:+.4f})")
        print(f"   quiet-rev mean (s_k < s*)  : {q['quiet_mean']:.6f}   "
              f"(n_quiet={q['n_quiet']})")
        print(f"   burst incidence            : {q['incidence']:.2f}   "
              f"episodes {q['episodes']}")
        print(f"   quiet-series drift         : {drift_pct:+.3f}% "
              f"± {ci_pct:.3f}% over window  (AR1 rho={rho:.2f}); "
              f"half-split {q['half_delta_pct']:+.3f}%")
        print(f"   detrended r(s_k, min_sr_k) : {sr_note}")
        print()

    if len(arms) > 1:
        print("cross-arm summary (drift-rate comparison; sigma from run names):")
        print(f"{'arm':42s} {'grand':>9s} {'quiet_reg':>9s} {'quiet_avg':>9s} "
              f"{'drift%':>8s} {'r(m,s)':>7s}")
        for arm in arms:
            q = quiet_stats(arm, sstar)
            slope, _, _ = q["trend"]
            drift_pct = (100.0 * slope * (arm["revs"][-1] - arm["revs"][0])
                         / q["quiet_mean"]) if q["n_quiet"] >= 5 \
                else float("nan")
            print(f"{arm['label']:42s} {_mean(arm['means']):9.6f} "
                  f"{arm['intercept']:9.6f} {q['quiet_mean']:9.6f} "
                  f"{drift_pct:8.3f} {arm['r_mean_std']:7.3f}")


if __name__ == "__main__":
    main()
