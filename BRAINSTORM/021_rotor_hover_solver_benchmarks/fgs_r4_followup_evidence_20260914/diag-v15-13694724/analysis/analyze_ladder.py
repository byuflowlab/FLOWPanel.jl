#!/usr/bin/env python3
"""
Thread-scaling analysis for FGS solver diagnostic ladder.
Analyzes j1, j4, j8, j16, j32, j64 arms; each with 40 trials (batches 1-4, 10 each).
Generates markdown summary with scaling, stage medians, Amdahl attribution, and anomalies.
"""
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

def pct(xs, q):
    """Compute percentile q in [0,1]."""
    xs = sorted(xs)
    if not xs:
        return float("nan")
    x = (len(xs)-1)*q
    lo, hi = math.floor(x), math.ceil(x)
    return xs[lo] if lo == hi else xs[lo] + (xs[hi]-xs[lo])*(x-lo)

def vals(name, rows):
    """Extract float values for a column name."""
    return [float(r[name]) for r in rows if r[name] and r[name].strip()]

def fmt(x, precision=6):
    """Format number."""
    if math.isnan(x):
        return "NaN"
    if isinstance(x, int):
        return str(x)
    if x < 0.001 or x > 1e6:
        return f"{x:.3e}"
    return f"{x:.{precision}g}"

# Root directory
ROOT = Path("/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/BRAINSTORM/021_rotor_hover_solver_benchmarks/fgs_r4_followup_evidence_20260914/diag-v15-13694724")

# Arm specifications
arms = ["j1", "j4", "j8", "j16", "j32", "j64"]

# Load trial data for all arms
arm_data = {}
for arm in arms:
    results_dir = ROOT / f"{arm}-b1" / "results"
    if not results_dir.exists():
        print(f"WARNING: {results_dir} not found", flush=True)
        continue

    trials = list(csv.DictReader((results_dir / "diagnostic_trials.csv").open()))
    uninstrumented = [r for r in trials if r["instrumented"].lower() == "false"]
    instrumented = [r for r in trials if r["instrumented"].lower() == "true"]

    arm_data[arm] = {
        "uninstrumented": uninstrumented,
        "instrumented": instrumented,
        "results_dir": results_dir
    }

# === TABLE 1: SCALING ===
print("\n## Table 1: Scaling (Uninstrumented Median Wall Time)")
print("| Arm | Median (s) | Speedup vs j1 | Parallel Efficiency | Instrumentation Overhead % |")
print("|---|---:|---:|---:|---:|")

j1_median = None
scaling_data = {}

for arm in arms:
    if arm not in arm_data:
        continue

    uninst_times = vals("solve_seconds", arm_data[arm]["uninstrumented"])
    if not uninst_times:
        continue

    median = pct(uninst_times, 0.5)
    scaling_data[arm] = {"median": median}

    if arm == "j1":
        j1_median = median

    if j1_median is None:
        print(f"| {arm} | {fmt(median, 4)} | — | — | — |")
    else:
        threads = int(arm[1:])
        speedup = j1_median / median
        efficiency = speedup / threads

        # Instrumentation overhead from audit
        # j1: -0.24%, j4: 0.89%, j8: -0.59%, j16: -2.15%, j32: 0.38%, j64: -0.13%
        overhead_pcts = {"j1": -0.24098, "j4": 0.890416, "j8": -0.594008,
                         "j16": -2.148, "j32": 0.384687, "j64": -0.131182}
        overhead_pct = overhead_pcts.get(arm, float("nan"))

        print(f"| {arm} | {fmt(median, 4)} | {fmt(speedup, 4)} | {fmt(efficiency, 4)} | {fmt(overhead_pct, 3)} |")

# === TABLE 2: STAGE MEDIANS (INSTRUMENTED) ===
print("\n## Table 2: Stage Medians from Instrumented Profiles")
print("| Arm | Total (s) | FMM Total (s) | Nonself Prod. (s) | Initialize (s) | Scatter (s) | Leaf Solve (s) |")
print("|---|---:|---:|---:|---:|---:|---:|")

stage_data = {}

for arm in arms:
    if arm not in arm_data:
        continue

    inst = arm_data[arm]["instrumented"]
    if not inst:
        continue

    stage_names = ["total_stage_seconds", "fmm_seconds", "nonself_product_seconds",
                   "initialization_seconds", "scatter_seconds", "leaf_solve_seconds"]

    stages = {}
    for name in stage_names:
        v = vals(name, inst)
        stages[name] = pct(v, 0.5) if v else float("nan")

    stage_data[arm] = stages

    print(f"| {arm} | {fmt(stages['total_stage_seconds'], 4)} | " +
          f"{fmt(stages['fmm_seconds'], 4)} | {fmt(stages['nonself_product_seconds'], 4)} | " +
          f"{fmt(stages['initialization_seconds'], 4)} | {fmt(stages['scatter_seconds'], 4)} | " +
          f"{fmt(stages['leaf_solve_seconds'], 4)} |")

# === TABLE 3: STAGE SPEEDUPS (vs j1) ===
print("\n## Table 3: Stage Speedups vs j1")
if "j1" in stage_data:
    j1_stages = stage_data["j1"]
    print("| Arm | FMM Speedup | Nonself Prod. Speedup | Initialize Speedup | Scatter Speedup | Leaf Solve Speedup |")
    print("|---|---:|---:|---:|---:|---:|")

    for arm in arms:
        if arm not in stage_data:
            continue

        stages = stage_data[arm]
        fmm_su = j1_stages["fmm_seconds"] / stages["fmm_seconds"] if stages["fmm_seconds"] > 0 else float("nan")
        nonself_su = j1_stages["nonself_product_seconds"] / stages["nonself_product_seconds"] if stages["nonself_product_seconds"] > 0 else float("nan")
        init_su = j1_stages["initialization_seconds"] / stages["initialization_seconds"] if stages["initialization_seconds"] > 0 else float("nan")
        scatter_su = j1_stages["scatter_seconds"] / stages["scatter_seconds"] if stages["scatter_seconds"] > 0 else float("nan")
        leaf_su = j1_stages["leaf_solve_seconds"] / stages["leaf_solve_seconds"] if stages["leaf_solve_seconds"] > 0 else float("nan")

        print(f"| {arm} | {fmt(fmm_su, 4)} | {fmt(nonself_su, 4)} | {fmt(init_su, 4)} | {fmt(scatter_su, 4)} | {fmt(leaf_su, 4)} |")

# === TABLE 4: ITERATIONS CONSTANCY ===
print("\n## Table 4: Iterations and Inner Sweeps per Arm")
print("| Arm | Iterations (p50) | Estimated Inner Sweeps (p50) | FMM Passes (p50) |")
print("|---|---:|---:|---:|")

for arm in arms:
    if arm not in arm_data:
        continue

    all_trials = arm_data[arm]["uninstrumented"] + arm_data[arm]["instrumented"]

    iters = vals("iterations", all_trials)
    sweeps = vals("estimated_inner_sweeps", all_trials)
    passes = vals("estimated_fmm_passes", all_trials)

    if iters:
        p50_iters = pct(iters, 0.5)
    else:
        p50_iters = float("nan")

    if sweeps:
        p50_sweeps = pct(sweeps, 0.5)
    else:
        p50_sweeps = float("nan")

    if passes:
        p50_passes = pct(passes, 0.5)
    else:
        p50_passes = float("nan")

    # Check for variance
    if iters:
        iters_min, iters_max = min(iters), max(iters)
        iters_flag = "" if iters_min == iters_max else f" (var: {iters_min:.0f}–{iters_max:.0f})"
    else:
        iters_flag = ""

    print(f"| {arm} | {fmt(p50_iters, 0)}{iters_flag} | {fmt(p50_sweeps, 0)} | {fmt(p50_passes, 0)} |")

# === TABLE 5: THREAD ACTIVITY SUMMARY ===
print("\n## Table 5: Thread Activity Summary (from validation CSV)")
print("(Rows: unique thread IDs; cpu_ticks sum over all samples in thread_activity.csv)")
print("| Arm | Thread Count | Min CPU Ticks | p50 CPU Ticks | Max CPU Ticks |")
print("|---|---:|---:|---:|---:|")

for arm in arms:
    if arm not in arm_data:
        continue

    results_dir = arm_data[arm]["results_dir"]
    activity_csv = results_dir / "thread_activity.csv"

    if not activity_csv.exists():
        print(f"| {arm} | — | — | — | — |")
        continue

    activity = list(csv.DictReader(activity_csv.open()))
    by_tid = defaultdict(list)

    for r in activity:
        try:
            ticks = int(r["cpu_ticks"])
            by_tid[r["tid"]].append(ticks)
        except (ValueError, KeyError):
            pass

    all_ticks = []
    for tid_ticks in by_tid.values():
        all_ticks.extend(tid_ticks)

    if all_ticks:
        min_ticks = min(all_ticks)
        p50_ticks = pct(all_ticks, 0.5)
        max_ticks = max(all_ticks)
        thread_count = len(by_tid)
    else:
        thread_count = len(by_tid)
        min_ticks = max_ticks = p50_ticks = float("nan")

    print(f"| {arm} | {thread_count} | {fmt(min_ticks, 0)} | {fmt(p50_ticks, 0)} | {fmt(max_ticks, 0)} |")

# === AMDAHL / SCALING ANOMALIES ===
print("\n## Amdahl-Style Attribution")

if "j1" in stage_data and "j64" in stage_data:
    j1_total = stage_data["j1"]["total_stage_seconds"]
    j64_total = stage_data["j64"]["total_stage_seconds"]
    j64_speedup = j1_total / j64_total

    print(f"Overall j1→j64 speedup (instrumented): {fmt(j64_speedup, 4)}x")
    print(f"j1 total: {fmt(j1_total, 4)}s, j64 total: {fmt(j64_total, 4)}s")

    # Find stages that plateau
    print("\nStages with limited scaling (flat from j32→j64):")
    j1_stages = stage_data["j1"]
    j32_stages = stage_data.get("j32", {})
    j64_stages = stage_data["j64"]

    for stage_name in ["initialization_seconds", "scatter_seconds", "exclusive_sum_seconds"]:
        j1_val = j1_stages.get(stage_name, float("nan"))
        j32_val = j32_stages.get(stage_name, float("nan"))
        j64_val = j64_stages.get(stage_name, float("nan"))

        if not math.isnan(j1_val) and not math.isnan(j64_val):
            ratio_j32_j64 = j32_val / j64_val if j64_val > 0 else float("nan")
            if 0.9 < ratio_j32_j64 < 1.1:  # Flat within 10%
                pct_j64 = 100.0 * j64_val / j64_total if j64_total > 0 else 0
                print(f"  - {stage_name}: {fmt(j1_val, 4)}s (j1) → {fmt(j32_val, 4)}s (j32) → {fmt(j64_val, 4)}s (j64), " +
                      f"ratio j32/j64={fmt(ratio_j32_j64, 3)}, {fmt(pct_j64, 1)}% of j64 wall time")

print("\n## Anomaly Summary")

# Check iteration variance
print("Iterations variance:")
variance_found = False
for arm in arms:
    if arm not in arm_data:
        continue
    all_trials = arm_data[arm]["uninstrumented"] + arm_data[arm]["instrumented"]
    iters = vals("iterations", all_trials)
    if iters and (max(iters) != min(iters)):
        print(f"  - {arm}: iterations vary from {min(iters):.0f} to {max(iters):.0f}")
        variance_found = True

if not variance_found:
    print("  - None detected; all arms constant at expected values.")

print("\nProfile CSV availability and row counts:")
for arm in arms:
    if arm not in arm_data:
        continue
    results_dir = arm_data[arm]["results_dir"]
    for fname in ["cpu_thread_complete_flat.txt", "thread_activity.csv", "gemv_census.csv"]:
        fpath = results_dir / fname
        if fpath.exists():
            size_mb = fpath.stat().st_size / (1024*1024)
            print(f"  - {arm}: {fname} exists ({fmt(size_mb, 2)} MB)")

print("\n## Notes")
print("- Instrumentation overhead: from audit-full.txt medians.")
print("- Stage times from instrumented trials only (batches 2 & 4).")
print("- Thread activity: cpu_ticks from validation, not wall time.")
print("- Script: analysis/analyze_ladder.py")
