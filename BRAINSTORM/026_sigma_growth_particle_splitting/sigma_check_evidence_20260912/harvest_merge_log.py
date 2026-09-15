import csv, math, sys
from collections import defaultdict

SHED = 0.031154
PATH = sys.argv[1]
LO, HI = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (100, 10**9)

rows = []
with open(PATH) as f:
    r = csv.DictReader(f)
    for row in r:
        p = int(row["pass"])
        if LO <= p <= HI:
            rows.append((p, float(row["sigma_i"]), float(row["sigma_j"]),
                         float(row["dist"]), float(row["sigma0_i"]), float(row["sigma0_j"])))

if not rows:
    print("no events in window"); sys.exit()

npass = len(set(p for p,*_ in rows))
print(f"window passes {LO}..{max(p for p,*_ in rows)}: {len(rows)} merge events over {npass} passes = {len(rows)/npass:.1f}/step")

def cls(s0):
    if abs(s0-SHED) < 0.02*SHED: return "shed"
    if s0 < 0.9*SHED: return "split-child"
    if s0 > 1.1*SHED: return "merge-product"
    return "near-shed"

# member classification
comp = defaultdict(int); n_mem = 0
at_release = 0; fresh = 0
smin_list = []
for p, si, sj, d, s0i, s0j in rows:
    smin_list.append(min(si, sj))
    for s, s0 in ((si, s0i), (sj, s0j)):
        comp[cls(s0)] += 1; n_mem += 1
        if abs(s-SHED) < 0.001*SHED and abs(s0-SHED) < 0.001*SHED:
            at_release += 1
        if s0 > 0 and abs(s/s0 - 1) < 0.05:
            fresh += 1
print("member sigma_0 composition: " + "  ".join(f"{k}={v/n_mem:.3f}" for k,v in sorted(comp.items())))
print(f"at-release members (sigma==sigma0==shed to 0.1%): {at_release/n_mem:.3f}")
print(f"fresh members (|sigma/sigma0 - 1| < 5%): {fresh/n_mem:.3f}")

smin_list.sort()
q = lambda f: smin_list[int(f*(len(smin_list)-1))]
print("pair sigma_min/shed quantiles 5/25/50/75/95: " +
      " ".join(f"{q(f)/SHED:.2f}" for f in (0.05,0.25,0.5,0.75,0.95)))
edges = [0.25*k for k in range(15)]
hist = [0]*14
for s in smin_list:
    b = min(13, int(s/SHED/0.25))
    hist[b] += 1
print("hist sigma_min/shed @0.25 bins from 0:", hist)

# pair-type classification (both members)
pair_types = defaultdict(int)
for p, si, sj, d, s0i, s0j in rows:
    key = tuple(sorted((cls(s0i), cls(s0j))))
    pair_types[key] += 1
print("pair type counts:", dict(sorted(pair_types.items(), key=lambda kv:-kv[1])))

# per-100-pass merge rate evolution
byp = defaultdict(int)
for p,*_ in rows: byp[p] += 1
ps = sorted(byp)
print("merge rate: first pass", ps[0], "last", ps[-1])
for lo in range(ps[0]//50*50, ps[-1]+1, 50):
    w = [byp.get(p,0) for p in range(lo, lo+50)]
    if any(w): print(f"  passes {lo}-{lo+49}: {sum(w)/50:.1f}/step")
