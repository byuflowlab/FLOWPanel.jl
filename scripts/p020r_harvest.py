#!/usr/bin/env python3
"""Harvest BRAINSTORM 020 Phase-2R wake-health and terminal sigma scores.

Run on the cluster, whose Python environment provides VTK:

    python3 scripts/p020r_harvest.py data/scr_p020r_geom_s020v 0.0023736

The script reads monitor CSVs as the source of record and the highest-numbered
particle VTP for the registered field-median sigma/sigma_shed score. It writes
``p020r_harvest.csv`` inside the run directory and prints the same one-row
record. A sustained margin event is defined reproducibly as at least 10
consecutive recorded steps with max_dtZ > 0.2.
"""

import csv
import glob
import math
import os
import re
import struct
import sys
import zlib


def numeric_step(path):
    match = re.search(r"\.(\d+)\.vtp$", path)
    return int(match.group(1)) if match else -1


def finite(value):
    return math.isfinite(value)


def first_sustained(rows, key, threshold, count=10):
    run = 0
    start = None
    for row in rows:
        value = float(row[key])
        if finite(value) and value > threshold:
            if run == 0:
                start = int(row["step"])
            run += 1
            if run >= count:
                return start
        else:
            run = 0
            start = None
    return None


def terminal_sigma(vtp_path, sigma_shed):
    try:
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(vtp_path)
        reader.Update()
        array = reader.GetOutput().GetPointData().GetArray("sigma")
        if array is None:
            raise RuntimeError(f"sigma array absent from {vtp_path}")
        sigma = vtk_to_numpy(array)
    except ModuleNotFoundError:
        # Minimal VTK XML appended-data reader for FLOWVPM's binary VTP output:
        # UInt64 headers, little endian, vtkZLibDataCompressor. Avoid adding a
        # heavyweight VTK dependency merely to extract one scalar array.
        with open(vtp_path, "rb") as stream:
            raw = stream.read()
        app_tag = raw.index(b"<AppendedData")
        underscore = raw.index(b"_", app_tag) + 1
        header = raw[:app_tag].decode("utf-8")
        match = re.search(
            r'<DataArray type="Float64" Name="sigma"[^>]*offset="(\d+)"',
            header,
        )
        if match is None:
            raise RuntimeError(f"sigma array absent from {vtp_path}")
        cursor = underscore + int(match.group(1))
        nblocks, block_size, last_size = struct.unpack_from("<QQQ", raw, cursor)
        cursor += 24
        sizes = struct.unpack_from(f"<{nblocks}Q", raw, cursor)
        cursor += 8 * nblocks
        chunks = []
        for size in sizes:
            chunks.append(zlib.decompress(raw[cursor:cursor + size]))
            cursor += size
        payload = b"".join(chunks)
        expected = (nblocks - 1) * block_size + last_size
        if len(payload) != expected or len(payload) % 8:
            raise RuntimeError(f"malformed sigma payload in {vtp_path}")
        sigma = struct.unpack(f"<{len(payload) // 8}d", payload)
    values = sorted(float(x) / sigma_shed for x in sigma)
    n = len(values)
    median = values[n // 2] if n % 2 else 0.5 * (values[n // 2 - 1] + values[n // 2])
    return n, values[0], values[max(0, math.ceil(0.01 * n) - 1)], median


def main():
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {sys.argv[0]} RUN_DIR SIGMA_SHED")
    run_dir = os.path.abspath(sys.argv[1])
    sigma_shed = float(sys.argv[2])
    tag = os.path.basename(run_dir.rstrip(os.sep))
    health = os.path.join(run_dir, "monitors", f"{tag}_monitor04_wake_health_system1.csv")
    with open(health, newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"no health rows in {health}")

    snapshots = glob.glob(os.path.join(
        run_dir, f"{tag}_wake1_particles", f"{tag}_wake1_particles.*.vtp"))
    if not snapshots:
        raise RuntimeError(f"no particle VTP snapshots under {run_dir}")
    terminal_vtp = max(snapshots, key=numeric_step)
    terminal_np, terminal_min_sr, terminal_p1_sr, terminal_median_sr = terminal_sigma(
        terminal_vtp, sigma_shed)

    def floats(key):
        return [float(row[key]) for row in rows
                if key in row and row[key] not in (None, "") and finite(float(row[key]))]
    p1_values = floats("p1_sigma_ratio")
    onset = next((int(row["step"]) for row in rows
                  if finite(float(row["max_u"])) and float(row["max_u"]) > 1000), None)
    summary = {
        "tag": tag,
        "last_step": int(rows[-1]["step"]),
        "n_health_rows": len(rows),
        "max_particles": max(int(float(row["n_particles"])) for row in rows),
        "terminal_particles_monitor": int(float(rows[-1]["n_particles"])),
        "terminal_particles_vtp": terminal_np,
        "negative_sigma_rows": sum(float(row["min_sigma"]) < 0 for row in rows
                                   if finite(float(row["min_sigma"]))),
        "min_sigma_ratio": min(floats("min_sigma_ratio")),
        "min_p1_sigma_ratio": min(p1_values) if p1_values else "",
        "max_u": max(floats("max_u")),
        "max_gamma_over_sigma2": max(floats("max_gamma_over_sigma2")),
        "max_dtZ": max(floats("max_dtZ")),
        "first_max_u_gt_1000_step": "" if onset is None else onset,
        "first_sustained_max_dtZ_gt_0p2_step": "" if (
            sustained := first_sustained(rows, "max_dtZ", 0.2)) is None else sustained,
        "terminal_vtp_step": numeric_step(terminal_vtp),
        "terminal_min_sigma_ratio": terminal_min_sr,
        "terminal_p1_sigma_ratio": terminal_p1_sr,
        "terminal_median_sigma_ratio": terminal_median_sr,
    }
    out = os.path.join(run_dir, "p020r_harvest.csv")
    with open(out, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary))
        writer.writeheader()
        writer.writerow(summary)
    history = os.path.join(run_dir, "p020r_harvest_history.csv")
    prior_steps = set()
    if os.path.exists(history):
        with open(history, newline="") as stream:
            prior_steps = {int(row["terminal_vtp_step"]) for row in csv.DictReader(stream)}
    if summary["terminal_vtp_step"] not in prior_steps:
        with open(history, "a", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(summary))
            if stream.tell() == 0:
                writer.writeheader()
            writer.writerow(summary)
    print(",".join(summary))
    print(",".join(str(summary[key]) for key in summary))
    print(f"wrote {out}")
    print(f"updated {history}")


if __name__ == "__main__":
    main()
