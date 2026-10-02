#!/usr/bin/env python3
"""Track C prep (2026-09-28): interface + span profile on the real R4 lower DAG.

Session-prep measurement for the Ryan Track C conversation — NOT the official
C-T1 deliverable (that would be benchmark/fgs_triangular_interface_profile.jl
per fgs_exact_triangular_solver_plan_20260926.md, with R2 + larger meshes).

Inputs: the A-T2 exports in BRAINSTORM/033_atheory_20260926/
  fgs_dag_L_graph_R4.csv  (leaf, nof=leaf size, ptot=total predecessor size)
  fgs_dag_L_edges_R4.csv  (target, source) lower block edges

Method (mirrors Stage A of the plan):
  - contiguous leaf-range partitions K in {2,4,8,16,32,64}, cut by cumulative
    weight w_i = n_i*ptot_i + n_i^2;
  - interface C = source leaves with >=1 outgoing cross-subdomain lower edge;
  - structural fill of S = I + R_C B^{-1} E_C by reachability propagation of
    each incoming interface column block through the target subdomain's local
    lower DAG;
  - work proxy: gemv n_t*n_s per block edge + n_t^2 per leaf back-solve;
  - Track C sweep ceiling T_C >= max-subdomain local span (solve 1)
    + interface solve (serial total-work AND parallel span variants)
    + max-subdomain local span with E_C z gemvs folded in (solve 2).
"""
import csv
from collections import defaultdict

d = "/Users/ryan/Dropbox/research/projects/FLOWPanel.jl/BRAINSTORM/033_atheory_20260926"

leaves = {}
with open(f"{d}/fgs_dag_L_graph_R4.csv") as f:
    for row in csv.DictReader(f):
        leaves[int(row["leaf"])] = (int(row["nof"]), float(row["ptot"]))
ids = sorted(leaves)
n = {i: leaves[i][0] for i in ids}
ptot = {i: leaves[i][1] for i in ids}
N = sum(n.values())
L = len(ids)

preds = defaultdict(list)
succs = defaultdict(list)
E = 0
with open(f"{d}/fgs_dag_L_edges_R4.csv") as f:
    for row in csv.DictReader(f):
        t, s = int(row["target"]), int(row["source"])
        preds[t].append(s); succs[s].append(t); E += 1

pc = sorted(len(preds[i]) for i in ids)
print(f"leaves={L} unknowns={N} lower_edges={E}")
print(f"preds/leaf: mean={E/L:.1f} median={pc[L//2]} p90={pc[int(0.9*L)]} max={pc[-1]}")
print(f"leaf size: mean={N/L:.1f} min={min(n.values())} max={max(n.values())}")

depth = {}
for i in ids:  # sorted ids are topological (block lower triangular)
    depth[i] = 1 + max((depth[s] for s in preds[i]), default=0)
print(f"original lower-DAG depth (leaf levels): {max(depth.values())}")

w = [n[i]*ptot[i] + n[i]**2 for i in ids]
W = sum(w)
cum = []
c = 0.0
for x in w:
    c += x; cum.append(c)

def partition(K):
    sub = {}
    k = 0
    for idx, i in enumerate(ids):
        while k < K-1 and cum[idx] > (k+1)*W/K:
            k += 1
        sub[i] = k
    return sub

def span(pred_of, extra=None):
    sp = {}
    for t in ids:
        cost = n[t]**2 + sum(n[t]*n[s] for s in pred_of[t]) + (extra.get(t, 0) if extra else 0)
        sp[t] = cost + max((sp[s] for s in pred_of[t] if s in sp), default=0)
    return max(sp.values(), default=0)

W0 = sum(n[t]**2 + sum(n[t]*n[s] for s in preds[t]) for t in ids)
S0 = span(preds)
print(f"\nbaseline work proxy W0={W0:.3e}, weighted span S0={S0:.3e}, "
      f"infinite-thread ceiling W0/S0={W0/S0:.2f}x")

hdr = ("K  crossE ifcL  ifc_unk ifc_frac Sblocks S_MB_f64 Sdepth locdepth ctorcols "
       "imbal  localspan  S_work    S_span   TC_serS   TC_parS  vsS0par")
print("\n" + hdr)
for K in (2, 4, 8, 16, 32, 64):
    sub = partition(K)
    ifc = set(); crossE = 0
    cross_by = defaultdict(set)
    xtra = defaultdict(float)
    for t in ids:
        for s in preds[t]:
            if sub[s] != sub[t]:
                crossE += 1
                ifc.add(s)
                cross_by[(sub[t], s)].add(t)
                xtra[t] += n[t]*n[s]
    ifc_unk = sum(n[i] for i in ifc)
    local_succ = defaultdict(list)
    lpred = {t: [s for s in preds[t] if sub[s] == sub[t]] for t in ids}
    for s in ids:
        for t in succs[s]:
            if sub[s] == sub[t]:
                local_succ[s].append(t)
    Sadj = defaultdict(set)
    S_blocks = 0; S_bytes = 0
    for (dsub, c0), entries in cross_by.items():
        seen = set(); stack = list(entries)
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            stack.extend(v for v in local_succ[u] if v not in seen)
        rows = [r for r in seen if r in ifc]
        S_blocks += len(rows)
        S_bytes += sum(n[r]*n[c0] for r in rows) * 8
        for r in rows:
            if r != c0:
                Sadj[r].add(c0)
    sd = {}
    for r in sorted(ifc):
        sd[r] = 1 + max((sd[c0] for c0 in Sadj[r] if c0 in sd), default=0)
    ld = {}; locdepth = 0
    for i in ids:
        ld[i] = 1 + max((ld[s] for s in lpred[i]), default=0)
        locdepth = max(locdepth, ld[i])
    ls = span(lpred)
    ls2 = span(lpred, extra=xtra)
    S_work = sum(n[r]*n[c0] for r in Sadj for c0 in Sadj[r])
    ssp = {}
    for r in sorted(ifc):
        cost = sum(n[r]*n[c0] for c0 in Sadj[r])
        ssp[r] = cost + max((ssp[c0] for c0 in Sadj[r] if c0 in ssp), default=0)
    S_span = max(ssp.values(), default=0)
    TC_ser = ls + S_work + ls2
    TC_par = ls + S_span + ls2
    sub_unk = defaultdict(int)
    for i in ids:
        sub_unk[sub[i]] += n[i]
    imbal = max(sub_unk.values())/(N/K)
    print(f"{K:<3d}{crossE:7d}{len(ifc):5d}{ifc_unk:9d}{ifc_unk/N:9.3f}{S_blocks:8d}"
          f"{S_bytes/1e6:9.1f}{max(sd.values(), default=0):7d}{locdepth:9d}"
          f"{len(cross_by):9d}{imbal:7.2f}{ls:11.2e}{S_work:10.2e}{S_span:10.2e}"
          f"{TC_ser:10.2e}{TC_par:10.2e}{S0/TC_par:8.2f}x")
