#!/usr/bin/env python3
"""Audit harvested P021 R4 FGS diagnostics (stdlib only).

Usage: audit_r4_diag.py RUN_DIRECTORY

This is a data audit, not provenance verification.  It accepts a completed
diagnostic root and reports arm completeness, correctness gates, timing
medians, stage-timer reconciliation, census files, and counter availability.
It never fills missing values or treats profile sample fractions as a wall-time
budget.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import statistics
import sys
import tomllib
from pathlib import Path

JOBS = (1, 4, 8, 16, 32, 64)
EXPECTED_BATCHES = (1, 2, 3, 4)
REPS = 10
TOL_BC = 1e-6
TOL_REPEAT = 1e-8
TOL_DIRECT = 1e-7
STAGES = (
    "initialization_seconds", "fmm_seconds", "influence_mapping_seconds",
    "residual_seconds", "leaf_solve_seconds", "nonself_product_seconds",
    "scatter_seconds", "remaining_iteration_seconds", "final_update_seconds",
)


def read_csv(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def fnum(row, key):
    """Return finite float or None; absent/blank/NaN is deliberately invalid."""
    value = row.get(key)
    if value is None or str(value).strip() == "":
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def truth(row, key):
    return str(row.get(key, "")).strip().lower() == "true"


def status_completed(path):
    try:
        with path.open("rb") as f:
            doc = tomllib.load(f)
        return doc.get("status") == "completed"
    except (OSError, tomllib.TOMLDecodeError):
        return False


def discover_arms(root):
    arms = {}
    for p in root.rglob("diagnostic_trials.csv"):
        match = re.search(r"(?:^|/)j(\d+)[_-]b1/results$", p.parent.as_posix())
        if match:
            j = int(match.group(1))
            arms[j] = p.parent
    return arms


def gate_rows(rows):
    failures = []
    for i, row in enumerate(rows, 1):
        if not truth(row, "solved"):
            failures.append(f"row {i}: solved != true")
        if not truth(row, "accepted"):
            failures.append(f"row {i}: accepted != true")
        if row.get("authoritative_evaluator") != "certified_fmm":
            failures.append(f"row {i}: authoritative_evaluator is not certified_fmm")
        bc = fnum(row, "authoritative_rel_l2")
        if bc is None or bc > TOL_BC:
            failures.append(f"row {i}: authoritative_rel_l2 missing/nonfinite/>1e-6")
        repeat = fnum(row, "relative_solution_delta")
        if repeat is None or repeat > TOL_REPEAT:
            failures.append(f"row {i}: relative_solution_delta missing/nonfinite/>1e-8")
        if row.get("direct_fmm_delta", "").strip().lower() not in ("", "nan"):
            direct = fnum(row, "direct_fmm_delta")
            if direct is None or direct > TOL_DIRECT:
                failures.append(f"row {i}: evaluated direct_fmm_delta missing/nonfinite/>1e-7")
        solve_time = fnum(row, "solve_seconds")
        if solve_time is None or solve_time <= 0:
            failures.append(f"row {i}: solve_seconds missing/nonfinite/nonpositive")
    return failures


def median_finite(rows, key):
    vals = [fnum(r, key) for r in rows]
    vals = [x for x in vals if x is not None]
    return statistics.median(vals) if vals else None


def fmt(x):
    return "missing" if x is None else f"{x:.6g}"


def written_precision(value):
    """Conservative absolute rounding error for a decimal CSV token."""
    s = str(value).strip().lower()
    if not s or s in ("nan", "inf", "+inf", "-inf"):
        return 0.0
    mantissa = s.split("e", 1)[0]
    places = len(mantissa.split(".", 1)[1]) if "." in mantissa else 0
    return 0.5 * 10.0 ** (-places)


def reconciliation_tolerance(row, keys):
    return max(1e-12, sum(written_precision(row.get(k, "")) for k in keys) + 1e-12)


def audit_arm(j, arm, root):
    failures = []
    trial_path = arm / "diagnostic_trials.csv"
    rows = read_csv(trial_path)
    batch_counts = {b: sum(str(r.get("batch", "")) == str(b) for r in rows) for b in EXPECTED_BATCHES}
    expected_modes = {1: "false", 2: "true", 3: "false", 4: "true"}
    pattern_ok = all([str(r.get("instrumented", "")).lower() for r in rows if str(r.get("batch", "")) == str(b)] == [expected_modes[b]] * REPS for b in EXPECTED_BATCHES)
    trial_sets_ok = all(sorted(int(r.get("trial", "-1")) for r in rows if str(r.get("batch", "")) == str(b) and str(r.get("trial", "")).isdigit()) == list(range(1, REPS + 1)) for b in EXPECTED_BATCHES)
    complete_rows = len(rows) == len(EXPECTED_BATCHES) * REPS and all(v == REPS for v in batch_counts.values()) and pattern_ok and trial_sets_ok
    if not complete_rows:
        failures.append(f"rows/batch or instrumented pattern invalid: counts={batch_counts}, expected modes=false,true,false,true, unique trials=1:10")
    failures.extend(gate_rows(rows))
    if not status_completed(arm / "status.toml"):
        failures.append("arm status.toml absent or status != completed")

    # Profile validation is an independent correctness check written by the
    # thread-complete profile section of fgs_r4_diagnostics.jl.
    pv = arm / "cpu_thread_complete_validation.csv"
    if not pv.is_file():
        failures.append("cpu_thread_complete_validation.csv absent")
    else:
        prow = read_csv(pv)
        if (len(prow) != 1 or not truth(prow[0], "eligible") or
                not truth(prow[0], "solved") or not truth(prow[0], "accepted") or
                prow[0].get("authoritative_evaluator") != "certified_fmm"):
            failures.append("profile validation absent or eligible != true")
        else:
            for key, limit in (("authoritative_rel_l2", TOL_BC), ("relative_solution_delta", TOL_REPEAT)):
                x = fnum(prow[0], key)
                if x is None or x > limit:
                    failures.append(f"profile validation {key} missing/nonfinite/>threshold")
    for name in ("cpu_thread_complete_flat.txt", "cpu_thread_complete_tree.txt"):
        if not (arm / name).is_file():
            failures.append(f"profile output absent: {name}")

    # Instrumentation control: exactly one uninstrumented and one instrumented
    # history, identical history and solution, with the same numerical gates.
    eq_path = arm / "instrumentation_equivalence.csv"
    eq_note = "missing"
    if not eq_path.is_file():
        failures.append("instrumentation_equivalence.csv absent")
    else:
        eq = read_csv(eq_path)
        modes = {str(r.get("instrumented", "")).lower() for r in eq}
        if len(eq) != 2 or modes != {"true", "false"}:
            failures.append("instrumentation equivalence must have exactly false/true rows")
        for r in eq:
            if (not truth(r, "accepted") or
                    r.get("authoritative_evaluator") != "certified_fmm"):
                failures.append("instrumentation equivalence numerical gate failed")
            x = fnum(r, "authoritative_rel_l2")
            if x is None or x > TOL_BC:
                failures.append("instrumentation equivalence authoritative_rel_l2 missing/nonfinite/>1e-6")
            # history_control invokes cold_validate(...; crosscheck=true): its
            # direct comparison is mandatory even though timed R4 rows may
            # legitimately leave direct_fmm_delta as NaN.
            for key, limit in (("direct_rel_l2", TOL_BC), ("direct_fmm_delta", TOL_DIRECT)):
                x = fnum(r, key)
                if x is None or x > limit:
                    failures.append(f"instrumentation equivalence {key} missing/nonfinite/>threshold")
            x = fnum(r, "solution_delta")
            counterpart = fnum(r, "counterpart_solution_rel_l2")
            if (x is None or x > TOL_REPEAT or counterpart is None or counterpart > TOL_REPEAT or
                    str(r.get("counterpart_history_identical", "")).lower() != "true"):
                failures.append("instrumentation changed history/solution or control is incomplete")
        if len(eq) == 2:
            eq_note = f"history_identical={str(eq[0].get('counterpart_history_identical')).lower()}"

    # Stage timers are independently measured.  Report only additive stage
    # totals and residual; no fractions are inferred from Profile samples.
    stage_rows = [r for r in rows if str(r.get("instrumented", "")).lower() == "true"]
    stage_reconcile = None
    wall_minus_stage = None
    if stage_rows:
        vals = []
        for r in stage_rows:
            total = fnum(r, "total_stage_seconds")
            excl = fnum(r, "exclusive_sum_seconds")
            unacc = fnum(r, "unaccounted_seconds")
            if total is None or excl is None or unacc is None:
                failures.append("stage timer row has absent/nonfinite total/exclusive/unaccounted")
            else:
                actual = [fnum(r, key) for key in STAGES]
                if any(x is None for x in actual):
                    failures.append("stage timer row has absent/nonfinite stage component")
                elif abs(sum(actual) - excl) > reconciliation_tolerance(r, STAGES + ("exclusive_sum_seconds",)):
                    failures.append("stage components do not reconcile to exclusive_sum_seconds")
                elif abs(total - excl - unacc) > reconciliation_tolerance(r, ("total_stage_seconds", "exclusive_sum_seconds", "unaccounted_seconds")):
                    failures.append("total_stage_seconds-exclusive_sum_seconds does not reconcile to unaccounted_seconds")
                vals.append((total, excl, unacc))
        if vals:
            stage_reconcile = tuple(statistics.median(v[i] for v in vals) for i in range(3))
            wall_minus_stage = median_finite(stage_rows, "solve_seconds") - stage_reconcile[0]
    else:
        failures.append("stage timer columns absent")

    # Census is availability only; values are not interpreted as provenance.
    census = [name for name in ("dependency_edges.csv", "gemv_census.csv", "gemv_census.toml") if (arm / name).is_file()]
    if len(census) != 3:
        failures.append("census incomplete (need dependency_edges.csv, gemv_census.csv, gemv_census.toml)")
    counter = "thread_activity.csv" if (arm / "thread_activity.csv").is_file() else "unavailable"
    root_counter_files = sorted(p.name for p in root.glob("hardware_counter_capability*") if p.is_file())
    uninst = [r for r in rows if str(r.get("instrumented", "")).lower() == "false"]
    inst = [r for r in rows if str(r.get("instrumented", "")).lower() == "true"]
    med_u = median_finite(uninst, "solve_seconds")
    med_i = median_finite(inst, "solve_seconds")
    overhead = None if med_u is None or med_i is None else med_i - med_u
    rel_overhead = None if med_u in (None, 0) else overhead / med_u
    return {
        "j": j, "arm": str(arm), "rows": len(rows), "batch_counts": batch_counts,
        "complete": not failures, "failures": failures, "median_uninstrumented": med_u,
        "median_instrumented": med_i, "overhead_seconds": overhead,
        "overhead_fraction": rel_overhead, "stage_reconcile": stage_reconcile,
        "wall_minus_stage": wall_minus_stage,
        "batch_medians": {b: median_finite([r for r in rows if str(r.get("batch", "")) == str(b)], "solve_seconds") for b in EXPECTED_BATCHES},
        "batch_spreads": {b: (max([fnum(r, "solve_seconds") for r in rows if str(r.get("batch", "")) == str(b) and fnum(r, "solve_seconds") is not None]) - min([fnum(r, "solve_seconds") for r in rows if str(r.get("batch", "")) == str(b) and fnum(r, "solve_seconds") is not None])) if any(fnum(r, "solve_seconds") is not None for r in rows if str(r.get("batch", "")) == str(b)) else None for b in EXPECTED_BATCHES},
        "eq": eq_note, "census": "+".join(census) if census else "unavailable",
        "counter": counter, "hardware_counter_capability": "+".join(root_counter_files) if root_counter_files else "unavailable",
    }


def print_report(root, results, root_ok):
    print(f"R4 diagnostic audit: {root}")
    print("data_audit_only=true (provenance is not verified)")
    print(f"root_COMPLETED={'PASS' if root_ok else 'FAIL'}")
    capability = sorted(p.name for p in root.glob("hardware_counter_capability*") if p.is_file())
    print("hardware_counter_capability=" + ("+".join(capability) if capability else "unavailable"))
    print("\n| arm | rows | batches | uninstrumented median s | instrumented median s | overhead s | overhead % | result |\n|---|---:|---|---:|---:|---:|---:|---|")
    for r in sorted(results, key=lambda x: x["j"]):
        pct = None if r["overhead_fraction"] is None else 100 * r["overhead_fraction"]
        batches = ",".join(f"{b}:{r['batch_counts'][b]}" for b in EXPECTED_BATCHES)
        print(f"| j{r['j']}/b1 | {r['rows']} | {batches} | {fmt(r['median_uninstrumented'])} | {fmt(r['median_instrumented'])} | {fmt(r['overhead_seconds'])} | {fmt(pct)} | {'PASS' if r['complete'] else 'FAIL'} |")
    print("\nChecks:")
    for r in sorted(results, key=lambda x: x["j"]):
        stage = r["stage_reconcile"]
        stage_text = "missing" if stage is None else "total/exclusive/unaccounted=" + "/".join(fmt(x) for x in stage) + f"; solve-total_stage={fmt(r['wall_minus_stage'])}"
        batches = "; ".join(f"b{b} median/spread={fmt(r.get('batch_medians', {}).get(b))}/{fmt(r.get('batch_spreads', {}).get(b))}" for b in EXPECTED_BATCHES)
        print(f"  j{r['j']}: profile + equivalence={r['eq']}; stage={stage_text}; {batches}; census={r['census']}; counters={r['counter']}")
        for failure in r["failures"]:
            print(f"    FAIL: {failure}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("run_directory", type=Path)
    args = ap.parse_args(argv)
    root = args.run_directory.resolve()
    if not root.is_dir():
        print(f"error: no run directory: {root}", file=sys.stderr)
        return 2
    root_ok = (root / "COMPLETED").is_file()
    arms = discover_arms(root)
    results = [audit_arm(j, arms[j], root) for j in JOBS if j in arms]
    for j in JOBS:
        if j not in arms:
            results.append({"j": j, "arm": "missing", "rows": 0, "batch_counts": {b: 0 for b in EXPECTED_BATCHES}, "complete": False, "failures": ["arm directory j%d-b1/results absent" % j], "median_uninstrumented": None, "median_instrumented": None, "overhead_seconds": None, "overhead_fraction": None, "stage_reconcile": None, "wall_minus_stage": None, "batch_medians": {}, "batch_spreads": {}, "eq": "missing", "census": "unavailable", "counter": "unavailable"})
    print_report(root, results, root_ok)
    return 0 if root_ok and all(r["complete"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
