#!/usr/bin/env python3
"""Create the BRAINSTORM 017 tables, fits, and publication-ready diagnostics."""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path
import statistics

os.environ.setdefault("MPLCONFIGDIR", "/tmp/flowpanel-mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "data" / "ssw_sheet_particle_split"
PROBE = ROOT / "data" / "ssw_representation_probe"
OUTPUT = SPLIT / "analysis"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def f(row: dict[str, str], name: str) -> float:
    return float(row[name])


def linfit_log(x: list[float], y: list[float]) -> tuple[float, float, float, int]:
    pairs = [(math.log(a), math.log(b)) for a, b in zip(x, y)
             if a > 0 and b > 0 and math.isfinite(a) and math.isfinite(b)]
    if len(pairs) < 2:
        return math.nan, math.nan, math.nan, len(pairs)
    xs, ys = zip(*pairs)
    xm, ym = statistics.mean(xs), statistics.mean(ys)
    sxx = sum((a - xm) ** 2 for a in xs)
    slope = sum((a - xm) * (b - ym) for a, b in pairs) / sxx
    intercept = ym - slope * xm
    residual = sum((b - (intercept + slope * a)) ** 2 for a, b in pairs)
    total = sum((b - ym) ** 2 for b in ys)
    r2 = 1 - residual / total if total else 1.0
    return math.exp(intercept), slope, r2, len(pairs)


def predictions() -> tuple[dict[tuple, tuple], dict[tuple, tuple]]:
    pps = {}
    for row in read_csv(PROBE / "phase_b_prediction.csv"):
        key = (f(row, "sigma_over_c"), int(row["pps"]),
               f(row, "buffer_over_c"))
        pps[key] = (f(row, "predicted_delta_CL_percent"),
                    f(row, "predicted_gamma_error_percent"))
    overlap = {}
    for row in read_csv(PROBE / "phase_b_prediction_sigma_overlap.csv"):
        key = (f(row, "sigma_over_c"), f(row, "overlap"),
               f(row, "buffer_over_c"))
        overlap[key] = (f(row, "predicted_delta_CL_percent"),
                        f(row, "predicted_gamma_error_percent"))
    return pps, overlap


def closure_table(summary: list[dict[str, str]]) -> list[dict]:
    pps, overlap = predictions()
    rows = []
    for row in summary:
        method = row["method"]
        if method == "sigma_pps":
            prediction = pps.get((f(row, "sigma_over_c"), int(row["pps"]),
                                  f(row, "buffer_over_c")))
        elif method == "sigma_overlap":
            prediction = overlap.get((f(row, "sigma_over_c"),
                                      f(row, "overlap"),
                                      f(row, "buffer_over_c")))
        else:
            prediction = None
        if prediction is None:
            continue
        observed_cl = f(row, "delta_CL_percent")
        observed_gamma = f(row, "gamma_error_percent")
        residual_cl = observed_cl - prediction[0]
        residual_gamma = observed_gamma - prediction[1]
        rows.append({
            "method": method,
            "sigma_over_c": row["sigma_over_c"],
            "pps": row["pps"],
            "overlap": row["overlap"],
            "buffer_over_c": row["buffer_over_c"],
            "predicted_delta_CL_percent": prediction[0],
            "observed_delta_CL_percent": observed_cl,
            "CL_residual_percentage_points": residual_cl,
            "predicted_gamma_error_percent": prediction[1],
            "observed_gamma_error_percent": observed_gamma,
            "gamma_residual_percentage_points": residual_gamma,
            "static_closure": abs(residual_cl) <= 0.25
                              and abs(residual_gamma) <= 0.5,
        })
    return rows


def scaling_fits() -> tuple[list[dict], list[dict[str, str]]]:
    field = read_csv(PROBE / "field_metrics.csv")
    rolled_body = [r for r in field if r["arm"] == "rolledup"
                   and r["selection"] == "tail" and r["probe"] == "body"]
    fits = []
    high_pps = [r for r in rolled_body if int(r["pps"]) == 8]
    asymptotic = [r for r in high_pps
                  if f(r, "sigma_over_c") / f(r, "buffer_over_c") <= .3]
    x = [f(r, "sigma_over_c") / f(r, "buffer_over_c")
         for r in asymptotic]
    for metric in ("velocity_max", "gradient_max"):
        coefficient, exponent, r2, count = linfit_log(
            x, [f(r, metric) for r in asymptotic])
        fits.append({"mechanism": "M1_kernel_smoothing", "metric": metric,
                     "variable": "sigma_over_d", "coefficient": coefficient,
                     "exponent": exponent, "r_squared": r2, "n": count})
    x4 = [f(r, "buffer_over_c") / f(r, "sigma_over_c") for r in high_pps]
    coefficient, exponent, r2, count = linfit_log(
        x4, [f(r, "velocity_max") for r in high_pps])
    fits.append({"mechanism": "M4_core_body_overlap",
                 "metric": "body_velocity_max", "variable": "d_over_sigma",
                 "coefficient": coefficient, "exponent": exponent,
                 "r_squared": r2, "n": count})

    excess_x, excess_y = [], []
    reference = {(f(r, "sigma_over_c"), f(r, "buffer_over_c")):
                 f(r, "velocity_max") for r in high_pps}
    for row in rolled_body:
        pps = int(row["pps"])
        if pps == 8:
            continue
        key = (f(row, "sigma_over_c"), f(row, "buffer_over_c"))
        excess = abs(f(row, "velocity_max") - reference[key])
        excess_x.append(0.125 / (pps * f(row, "sigma_over_c")))
        excess_y.append(excess)
    coefficient, exponent, r2, count = linfit_log(excess_x, excess_y)
    fits.append({"mechanism": "M2_quadrature_lumping",
                 "metric": "velocity_excess_vs_pps8",
                 "variable": "h_trailing_over_sigma",
                 "coefficient": coefficient, "exponent": exponent,
                 "r_squared": r2, "n": count})
    return fits, rolled_body


def history(path: Path) -> dict[float, float]:
    return {f(row, "time_star"): f(row, "CL") for row in read_csv(path)}


def startup_table(summary: list[dict[str, str]]) -> tuple[list[dict], Path | None]:
    control_candidates = [p for p in SPLIT.glob("*/history.csv")
                          if "_pw" not in p.parent.name]
    if not control_candidates:
        return [], None
    control_path = control_candidates[0]
    control = history(control_path)
    rows = []
    for row in summary:
        case_path = SPLIT / row["tag"] / "history.csv"
        if not case_path.exists():
            continue
        case = history(case_path)
        times = sorted(t for t in control.keys() & case.keys() if 0 < t <= 4)
        errors = [case[t] - control[t] for t in times]
        rows.append({
            "method": row["method"], "sigma_over_c": row["sigma_over_c"],
            "pps": row["pps"], "overlap": row["overlap"],
            "buffer_over_c": row["buffer_over_c"],
            "startup_CL_RMS": math.sqrt(sum(e * e for e in errors) / len(errors)),
            "startup_CL_max_abs": max(map(abs, errors)),
            "startup_CL_peak_signed": max(errors, key=abs),
            "samples_through_t4": len(times), "tag": row["tag"],
        })
    return rows, control_path


def make_plots(summary: list[dict[str, str]], closure: list[dict],
               rolled_body: list[dict[str, str]], startup: list[dict],
               control_path: Path | None) -> None:
    methods = sorted({r["method"] for r in summary})
    fig, axes = plt.subplots(1, len(methods), figsize=(5 * len(methods), 4),
                             squeeze=False)
    color_max = max(1, max(f(r, "gamma_error_percent") for r in summary))
    color_artist = None
    for ax, method in zip(axes[0], methods):
        subset = [r for r in summary if r["method"] == method]
        for row in subset:
            admissible = row["admissible"].lower() == "true"
            eligible = row["eligible"].lower() == "true"
            color_artist = ax.scatter(
                f(row, "sigma_over_c"), f(row, "buffer_over_c"),
                c=f(row, "gamma_error_percent"), cmap="viridis",
                vmin=0, vmax=color_max,
                marker="o" if admissible else ("s" if eligible else "x"),
                s=70)
        ax.set(xscale="log", yscale="log", xlabel=r"$\sigma/c$",
               ylabel=r"$d/c$", title=method)
        ax.grid(True, which="both", alpha=.25)
    fig.suptitle("Admissibility map: circles pass, squares eligible/fail, x mechanism-only")
    color_ax = fig.add_axes([.92, .18, .018, .58])
    fig.colorbar(color_artist, cax=color_ax,
                 label="max station Γ error (%)")
    fig.subplots_adjust(left=.08, right=.88, bottom=.16, top=.78, wspace=.22)
    fig.savefig(OUTPUT / "admissibility_map.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for pps in sorted({int(r["pps"]) for r in rolled_body}):
        points = [r for r in rolled_body if int(r["pps"]) == pps]
        axes[0].loglog([f(r, "sigma_over_c") / f(r, "buffer_over_c")
                        for r in points],
                       [f(r, "velocity_max") for r in points], "o",
                       label=f"PPS={pps}")
    axes[0].set(xlabel=r"$\sigma/d$", ylabel="body max |Δu|/U",
                title="M1/M4 field scaling")
    axes[0].grid(True, which="both", alpha=.25)
    axes[0].legend()
    axes[1].scatter([abs(f(r, "CL_residual_percentage_points"))
                     for r in closure],
                    [abs(f(r, "gamma_residual_percentage_points"))
                     for r in closure], s=25)
    axes[1].axvline(.25, color="k", linestyle="--")
    axes[1].axhline(.5, color="k", linestyle="--")
    axes[1].set(xlabel="|CL frozen-dynamic residual| (points)",
                ylabel="|Γ frozen-dynamic residual| (points)",
                title="Static closure")
    axes[1].grid(True, alpha=.25)
    fig.tight_layout()
    fig.savefig(OUTPUT / "scaling_and_closure.png", dpi=180)
    plt.close(fig)

    if startup and control_path is not None:
        control = history(control_path)
        worst = max(startup, key=lambda r: r["startup_CL_max_abs"])
        case = history(SPLIT / worst["tag"] / "history.csv")
        times = sorted(control.keys() & case.keys())
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(times, [control[t] for t in times], label="panel control")
        ax.plot(times, [case[t] for t in times],
                label=f'{worst["method"]}, σ/c={worst["sigma_over_c"]}')
        ax.set(xlabel=r"$tU/c$", ylabel=r"$C_L$",
               title="Worst startup-transient discrepancy")
        ax.grid(True, alpha=.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUTPUT / "startup_transient.png", dpi=180)
        plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    summary = read_csv(SPLIT / "split_summary.csv")
    closure = closure_table(summary)
    fits, rolled_body = scaling_fits()
    startup, control_path = startup_table(summary)
    write_csv(OUTPUT / "admissibility_table.csv", summary)
    write_csv(OUTPUT / "prediction_observation_residuals.csv", closure)
    write_csv(OUTPUT / "scaling_fits.csv", fits)
    write_csv(OUTPUT / "startup_transient_metrics.csv", startup)

    matched_path = SPLIT / "matched_geometry" / "matched_geometry_summary.csv"
    mechanism = read_csv(matched_path) if matched_path.exists() else closure
    write_csv(OUTPUT / "mechanism_budget.csv", mechanism)
    make_plots(summary, closure, rolled_body, startup, control_path)

    admissible = [r for r in summary if r["admissible"].lower() == "true"]
    settled = [r for r in summary if r["settled"].lower() == "true"]
    closed = [r for r in closure if str(r["static_closure"]).lower() == "true"]
    with (OUTPUT / "study_summary.md").open("w") as stream:
        stream.write("# BRAINSTORM 017 analysis summary\n\n")
        stream.write(f"- Dynamic cells: {len(summary)}; settled: {len(settled)}; "
                     f"admissible: {len(admissible)}.\n")
        stream.write(f"- Frozen/dynamic closure: {len(closed)}/{len(closure)} cells "
                     "meet 0.25 CL and 0.5 Γ percentage-point tolerances.\n")
        if admissible:
            stream.write("- Admissible cells:\n")
            for row in admissible:
                stream.write(f"  - {row['method']}: σ/c={row['sigma_over_c']}, "
                             f"d/c={row['buffer_over_c']}, "
                             f"ΔCL={f(row, 'delta_CL_percent'):+.4f}%, "
                             f"ΔΓ={f(row, 'gamma_error_percent'):.4f}%.\n")
        stream.write("- Scaling fits (`error = coefficient × variable^exponent`):\n")
        for fit in fits:
            stream.write(f"  - {fit['mechanism']} {fit['metric']}: "
                         f"{fit['coefficient']:.5g} × {fit['variable']}^"
                         f"{fit['exponent']:.4g}, R²={fit['r_squared']:.4g}.\n")
        if startup:
            worst = max(startup, key=lambda r: r["startup_CL_max_abs"])
            stream.write(f"- Worst startup |ΔCL|={worst['startup_CL_max_abs']:.6g} "
                         f"for {worst['method']} σ/c={worst['sigma_over_c']}, "
                         f"d/c={worst['buffer_over_c']}.\n")
        stream.write("- Matched-geometry mechanism budget: "
                     + ("available.\n" if matched_path.exists()
                        else "pending controller runs.\n"))


if __name__ == "__main__":
    main()
