#!/usr/bin/env python3
"""Reduce particle VTPs to sigma/sigma_shed quantiles + histogram.
Usage: vtp_sigma_reduce.py <run> <sigma_shed_m> <step0> <step1> <outtag>
Reads data/<run>/<run>_wake1_particles/<run>_wake1_particles.<s>.vtp
(appended raw, zlib-compressed, UInt64 headers). Emits
p019_p4b_<outtag>_quantiles.csv and p019_p4b_<outtag>_hist.csv in cwd."""
import sys, re, zlib, struct, csv, math

def read_sigma(path):
    blob = open(path, "rb").read()
    m = re.search(rb'Name="sigma"[^>]*offset="(\d+)"', blob)
    off = int(m.group(1))
    start = blob.index(b'<AppendedData encoding="raw">')
    us = blob.index(b"_", start) + 1
    p = us + off
    nb, bs, lbs = struct.unpack_from("<QQQ", blob, p)
    sizes = struct.unpack_from("<%dQ" % nb, blob, p + 24)
    q = p + 24 + 8 * nb
    out = b""
    for s in sizes:
        out += zlib.decompress(blob[q:q + s]); q += s
    return struct.unpack("<%dd" % (len(out) // 8), out)

run, sshed, s0, s1, tag = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
NB, LO, HI = 60, 0.0, 1.5
hist = [0] * NB
qrows = []
allv = []
for s in range(s0, s1 + 1):
    v = [x / sshed for x in read_sigma(
        f"data/{run}/{run}_wake1_particles/{run}_wake1_particles.{s}.vtp")]
    v.sort()
    n = len(v)
    def q(p): return v[min(n - 1, int(p * n))]
    qrows.append([s, n, f"{v[0]:.5f}", f"{q(.01):.5f}", f"{q(.05):.5f}",
                  f"{q(.25):.5f}", f"{q(.50):.5f}", f"{q(.75):.5f}",
                  f"{q(.95):.5f}", f"{v[-1]:.5f}"])
    for x in v:
        b = int((x - LO) / (HI - LO) * NB)
        if 0 <= b < NB: hist[b] += 1
    allv.extend((v[0], q(.25), q(.5), q(.75)))
with open(f"p019_p4b_{tag}_quantiles.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["step","n","min","p1","p5","p25","median","p75","p95","max"])
    w.writerows(qrows)
tot = sum(hist)
with open(f"p019_p4b_{tag}_hist.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["bin_lo","bin_hi","count","frac"])
    for i, c in enumerate(hist):
        w.writerow([f"{LO+(HI-LO)*i/NB:.4f}", f"{LO+(HI-LO)*(i+1)/NB:.4f}",
                    c, f"{c/tot:.6e}"])
meds = sorted(float(r[6]) for r in qrows)
print(f"{tag}: {len(qrows)} steps, window-median-of-medians "
      f"{meds[len(meds)//2]:.4f}, pooled n {tot}")
