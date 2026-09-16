#!/usr/bin/env python3
import csv, math, re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent / "Dropbox"  # unused
# Pass the harvested results directory explicitly.
import sys
if len(sys.argv) != 2:
    raise SystemExit("usage: analyze_j1.py <.../diag-v15-13694724/j1-b1/results>")
d = Path(sys.argv[1]).resolve()
trials = list(csv.DictReader((d / "diagnostic_trials.csv").open()))
inst = [r for r in trials if r["instrumented"].lower() == "true"]
census = list(csv.DictReader((d / "gemv_census.csv").open()))
edges = list(csv.DictReader((d / "dependency_edges.csv").open()))
def vals(name, rows):
    return [float(r[name]) for r in rows]
def pct(xs, q):
    xs = sorted(xs)
    if not xs: return float("nan")
    x = (len(xs)-1)*q
    lo, hi = math.floor(x), math.ceil(x)
    return xs[lo] if lo == hi else xs[lo] + (xs[hi]-xs[lo])*(x-lo)
def fmt(x):
    return f"{x:.6g}" if isinstance(x, float) else str(x)

print("STAGE_MEDIANS_INSTRUMENTED")
total = pct(vals("total_stage_seconds", inst), .5)
for name in ["initialization_seconds","fmm_seconds","influence_mapping_seconds","residual_seconds",
             "leaf_solve_seconds","nonself_product_seconds","scatter_seconds",
             "remaining_iteration_seconds","final_update_seconds","exclusive_sum_seconds","unaccounted_seconds"]:
    m = pct(vals(name, inst), .5)
    print(f"{name}\tmedian={m:.9f}\tfrac_of_median_total={m/total:.6f}")
print(f"total_stage_seconds\tmedian={total:.9f}")
print("STAGE_FRACTION_CAVEAT\tFractions use instrumented diagnostic medians; never speed claims.")

print("CENSUS_SUMMARY")
for name in ["m","n","matrix_bytes","matrix_elements","scatter_entries","target_interactions","dependent_leaves","conflict_degree","potential_color"]:
    xs = [float(r[name]) for r in census]
    print(f"{name}\tmin={fmt(min(xs))}\tp50={fmt(pct(xs,.5))}\tp90={fmt(pct(xs,.9))}\tp99={fmt(pct(xs,.99))}\tmax={fmt(max(xs))}\ttotal={fmt(sum(xs))}")
print(f"leaf_count\t{len(census)}")
print("toml_expected_matrix_bytes\t2862850032")
print("toml_expected_matrix_elements\t357856254")
print("toml_expected_directed_dependency_count\t95390")
print("toml_expected_scatter_entries\t5431340")

# Independently cross-check exported directed dependency edges and colors.
adj = defaultdict(set)
for e in edges:
    s, t = int(e["source_leaf"]), int(e["dependent_leaf"])
    adj[s].add(t)
edge_set = {(int(e["source_leaf"]), int(e["dependent_leaf"])) for e in edges}
duplicate_count = len(edges) - len(edge_set)
reported_dep = {int(r["leaf"]): int(r["dependent_leaves"]) for r in census}
reported_deg = {int(r["leaf"]): int(r["conflict_degree"]) for r in census}
reported_color = {int(r["leaf"]): int(r["potential_color"]) for r in census}
und = defaultdict(set)
for s,t in edge_set:
    if s != t:
        und[s].add(t); und[t].add(s)
# Conflict graph includes self? Edges are directed dependencies; conflict degree is
# symmetrized neighbors, with self excluded.
degree_mismatch = [(k, len(und[k]), reported_deg[k]) for k in reported_deg if len(und[k]) != reported_deg[k]]
dep_mismatch = [(k, len(adj[k]), reported_dep[k]) for k in reported_dep if len(adj[k]) != reported_dep[k]]
color_pairs = []
for s,t in edge_set:
    if s != t and reported_color.get(s) == reported_color.get(t):
        color_pairs.append((s,t,reported_color[s]))
color_groups = Counter(reported_color.values())
print("DEPENDENCY_CROSSCHECK")
print(f"edge_rows={len(edges)} unique_edges={len(edge_set)} duplicates={duplicate_count}")
print(f"directed_outdegree_mismatches={len(dep_mismatch)}")
print(f"symmetrized_degree_mismatches={len(degree_mismatch)}")
print(f"same_color_conflicts={len(color_pairs)}")
print(f"color_count={len(color_groups)} min_group={min(color_groups.values())} p50_group={pct(list(color_groups.values()),.5):.6g} max_group={max(color_groups.values())}")
print("color_groups=" + ",".join(f"{k}:{v}" for k,v in sorted(color_groups.items())))

# Profile sample coverage and thread activity.
profile = (d / "cpu_thread_complete_flat.txt").read_text()
threads = re.findall(r"^Thread (\d+).*?Total snapshots: (\d+).*?Utilization: ([0-9.]+)%", profile, re.M)
print("PROFILE_COVERAGE")
for tid, snaps, util in threads:
    print(f"profile_thread={tid} snapshots={snaps} utilization_percent={util}")
v = list(csv.DictReader((d / "thread_activity.csv").open()))
by_tid = defaultdict(list)
for r in v: by_tid[r["tid"]].append(int(r["cpu_ticks"]))
for tid, rows in sorted(by_tid.items()):
    print(f"validation_tid={tid} rows={len(rows)} cpu_ticks_min={min(rows)} p50={pct(rows,.5):.6g} max={max(rows)}")
print("PROFILE_CAVEAT\tProfile samples are coverage only, not wall-time fractions.")

# Provenance key extraction, compact.
prov = (d / "provenance.toml").read_text()
for key in ["julia_version","hostname","julia_threads","requested_julia_threads",
            "blas_threads","requested_blas_threads","rhs_sha256","manifest_sha256",
            "prepared_only","timing_scope"]:
    m = re.search(rf"^{re.escape(key)}\s*=\s*(.*)$", prov, re.M)
    if m: print(f"PROV_{key}={m.group(1)}")
for section in ["FLOWPanel","FastMultipole","FLOWVPM"]:
    block = re.search(rf"\[packages\.{section}\](.*?)(?=\n\[|\Z)", prov, re.S)
    if block:
        print("PROV_PACKAGE_" + section)
        for key in ["path","tag","sha","content_manifest_sha256"]:
            m = re.search(rf"^{key}\s*=\s*(.*)$", block.group(1), re.M)
            if m: print(f"{key}={m.group(1)}")
