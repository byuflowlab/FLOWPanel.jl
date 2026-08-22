#!/usr/bin/env python3
"""p018 sigma-collapse forensics (Track F).

Post-mortem of the diverged viscous rescue run p018_L2_visc (Slurm 13047296,
OOM in merge_particles! after step 1041) against the healthy warm-start parent
p018_L1_ov3 step 719.

Analyses (see BRAINSTORM 018):
  1. (axial/R, radial/R) maps of count, median/min sigma/sigma_shed, max |G|/s^2
     for corpse step 1041 and healthy step 719; census of the collapsing
     population (sigma/sigma_shed < {0.5, 0.25, 0.1}).
  2. Death trajectory over steps 1032..1041: global min sigma, argmin particle
     position/|Gamma|/NN-dist, field max |G|/s^2, recorded max |u|.
  3. Merge-candidacy audit of the collapsing set (sigma/sigma_shed<0.25) at
     step 1041 vs the absolute merge radius r_merge = 0.0027*R and vs 2*sigma.
  4. (baseline only) Z = Ghat.S.Ghat from the *recorded* velocity_gradient
     field; the from-scratch kernel recomputation + source-exclusion feedback
     test live in scripts/p018_sigma_forensics_strain.jl.

Reads VTK XML PolyData with zlib-compressed appended binary directly (no vtk
module needed). Writes CSVs and PNGs to data/p018_L2_visc_forensics/.

Usage: python3 scripts/p018_sigma_forensics.py
"""

import os
import re
import zlib
import struct
import numpy as np
from scipy.spatial import cKDTree

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(REPO, "data", "p018_L2_visc_forensics")
CORPSE_DIR = OUTDIR  # the ten retained VTPs were copied here
HEALTHY_VTP = os.path.expanduser(
    "~/p018_L1_ov3_paraview/p018_L1_ov3_wake1_particles/"
    "p018_L1_ov3_wake1_particles.719.vtp")

R = 0.119                      # rotor radius [m]
RPM = 5400.0
NT = 36
DT = 60.0 / (RPM * NT)         # 3.0864e-4 s (both runs)
SIGMA_SHED_CORPSE = 2 * np.pi * R * 2.88 / (NT * 26)   # 2.300e-3 m (L2_visc)
SIGMA_SHED_HEALTHY = 4.45e-3                           # m (L1_ov3, given)
R_MERGE = 0.0027 * R                                   # 3.213e-4 m absolute

STEPS = list(range(1032, 1042))


# ---------------------------------------------------------------------------
# Minimal VTK XML PolyData reader (zlib appended, header_type UInt64)
# ---------------------------------------------------------------------------

def read_vtp(path, fields=("gamma", "sigma", "velocity", "velocity_gradient")):
    with open(path, "rb") as f:
        raw = f.read()
    iapp = raw.index(b"<AppendedData")
    istart = raw.index(b"_", iapp) + 1
    xml = raw[:iapp].decode("utf-8", errors="replace")
    out = {}
    specs = {}
    for m in re.finditer(r'<DataArray[^>]*Name="([^"]+)"[^>]*offset="(\d+)"',
                         xml):
        tag = m.group(0)
        name = m.group(1)
        ncomp = int(re.search(r'NumberOfComponents="(\d+)"', tag).group(1))
        dtype = {"Float64": np.float64, "Float32": np.float32,
                 "Int64": np.int64}[re.search(r'type="(\w+)"', tag).group(1)]
        specs[name] = (int(m.group(2)), ncomp, dtype)

    def decode(offset):
        p = istart + offset
        nblk, blksz, lastsz = struct.unpack_from("<3Q", raw, p)
        sizes = struct.unpack_from("<%dQ" % nblk, raw, p + 24)
        p += 24 + 8 * nblk
        chunks = []
        for s in sizes:
            chunks.append(zlib.decompress(raw[p:p + s]))
            p += s
        return b"".join(chunks)

    for name in ("Points",) + tuple(fields):
        off, ncomp, dtype = specs[name]
        a = np.frombuffer(decode(off), dtype=dtype)
        out[name] = a.reshape(-1, ncomp) if ncomp > 1 else a
    return out


def load(path):
    d = read_vtp(path)
    X, G, s = d["Points"], d["gamma"], d["sigma"]
    live = np.linalg.norm(G, axis=1) > 0        # drop zero-strength placeholders
    return {"X": X[live], "G": G[live], "sigma": s[live],
            "u": d["velocity"][live],
            "J": d["velocity_gradient"][live],
            "ndrop": int((~live).sum())}


def cyl(X):
    """axial coord = x (thrust axis), radial = sqrt(y^2+z^2); returns in units
    of R with axial sign flipped if needed so the wake extends to positive
    values (downstream)."""
    ax = X[:, 0] / R
    rad = np.sqrt(X[:, 1] ** 2 + X[:, 2] ** 2) / R
    if np.median(ax) < 0:
        ax = -ax
    return ax, rad


def zstrain(G, J):
    """Z = Ghat . S . Ghat with S=(J+J^T)/2, J as recorded (row-major 3x3)."""
    Gh = G / np.maximum(np.linalg.norm(G, axis=1, keepdims=True), 1e-300)
    Jm = J.reshape(-1, 3, 3)
    S = 0.5 * (Jm + np.transpose(Jm, (0, 2, 1)))
    return np.einsum("ni,nij,nj->n", Gh, S, Gh)


def q(a, p):
    return np.quantile(a, p) if len(a) else np.nan


# ---------------------------------------------------------------------------
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    corpse = {st: load(os.path.join(
        CORPSE_DIR, "p018_L2_visc_wake1_particles.%d.vtp" % st))
        for st in STEPS}
    healthy = load(HEALTHY_VTP)

    # ---- gate (c): healthy min sigma/sigma_shed well above corpse ----------
    hmin = healthy["sigma"].min() / SIGMA_SHED_HEALTHY
    cmin = corpse[1041]["sigma"].min() / SIGMA_SHED_CORPSE
    habs, cabs = healthy["sigma"].min(), corpse[1041]["sigma"].min()
    print("GATE(c): healthy min sigma = %.4e m (rel %.4g) | corpse-1041 min "
          "sigma = %.4e m (rel %.4g); ratios abs %.2fx rel %.2fx -> %s"
          % (habs, hmin, cabs, cmin, habs / cabs, hmin / cmin,
             "PASS" if habs > 2 * cabs and hmin > 2 * cmin else "FAIL"))

    # ---- Analysis 1: maps --------------------------------------------------
    def make_map(d, sshed, label):
        ax, rad = cyl(d["X"])
        gmag = np.linalg.norm(d["G"], axis=1)
        srel = d["sigma"] / sshed
        gos2 = gmag / d["sigma"] ** 2
        axb = np.arange(-0.5, 4.51, 0.25)
        rdb = np.arange(0.0, 2.01, 0.1)
        ia = np.clip(np.digitize(ax, axb) - 1, 0, len(axb) - 2)
        ir = np.clip(np.digitize(rad, rdb) - 1, 0, len(rdb) - 2)
        rows = []
        for i in range(len(axb) - 1):
            for j in range(len(rdb) - 1):
                m = (ia == i) & (ir == j)
                n = int(m.sum())
                if n == 0:
                    continue
                rows.append((0.5 * (axb[i] + axb[i + 1]),
                             0.5 * (rdb[j] + rdb[j + 1]), n,
                             np.median(srel[m]), srel[m].min(),
                             gos2[m].max()))
        arr = np.array(rows)
        np.savetxt(os.path.join(OUTDIR, "map_%s.csv" % label), arr,
                   delimiter=",", header="axial_over_R,radial_over_R,count,"
                   "median_sig_rel,min_sig_rel,max_Gos2", comments="")
        return ax, rad, srel, gmag, gos2, arr

    axc, radc, srelc, gmc, gos2c, mapc = make_map(
        corpse[1041], SIGMA_SHED_CORPSE, "corpse1041")
    axh, radh, srelh, gmh, gos2h, maph = make_map(
        healthy, SIGMA_SHED_HEALTHY, "healthy719")

    print("\n== Analysis 1: collapsing-population census (corpse step 1041, "
          "N=%d live, %d dropped) ==" % (len(srelc), corpse[1041]["ndrop"]))
    print("thresh  count   frac      med_ax/R  med_rad/R  med|G|      "
          "p95|G|      med_sig/shed")
    for th in (0.5, 0.25, 0.1):
        m = srelc < th
        n = int(m.sum())
        print("%5.2f  %6d  %.5f   %8.3f  %8.3f   %.3e  %.3e  %.4f"
              % (th, n, n / len(srelc), q(axc[m], .5), q(radc[m], .5),
                 q(gmc[m], .5), q(gmc[m], .95), q(srelc[m], .5)))
    mh = srelh < 0.5
    print("healthy719 count below 0.5: %d (frac %.5f); healthy min sig/shed "
          "%.4g" % (mh.sum(), mh.mean(), srelh.min()))

    # where do the collapsers live (bin list, thresh 0.25)
    m25 = srelc < 0.25
    print("\ncollapser (sig/shed<0.25) location quantiles: ax/R p5/p50/p95 = "
          "%.2f/%.2f/%.2f ; rad/R = %.2f/%.2f/%.2f"
          % (q(axc[m25], .05), q(axc[m25], .5), q(axc[m25], .95),
             q(radc[m25], .05), q(radc[m25], .5), q(radc[m25], .95)))

    # figure: maps
    fig, axs = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    for k, (a, r, s, ttl) in enumerate([
            (axc, radc, srelc, "corpse 1041"), (axh, radh, srelh,
                                                "healthy 719")]):
        sc = axs[0, k].scatter(a, r, c=np.log10(s), s=.5, cmap="viridis",
                               vmin=-2, vmax=0.5)
        axs[0, k].set_title("%s: log10 sigma/shed" % ttl)
        plt.colorbar(sc, ax=axs[0, k])
    for k, (a, r, g, ttl) in enumerate([
            (axc, radc, gos2c, "corpse 1041"), (axh, radh, gos2h,
                                                "healthy 719")]):
        sc = axs[1, k].scatter(a, r, c=np.log10(np.maximum(g, 1e-12)), s=.5,
                               cmap="magma")
        axs[1, k].set_title("%s: log10 |G|/sigma^2" % ttl)
        plt.colorbar(sc, ax=axs[1, k])
    for a in axs.flat:
        a.set_xlabel("axial/R")
        a.set_ylabel("radial/R")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "maps_sigma_gos2.png"), dpi=140)
    plt.close(fig)

    # ---- Analysis 2: death trajectory -------------------------------------
    print("\n== Analysis 2: death trajectory 1032->1041 ==")
    print("step  N       min_sig[m]  sig/shed  ax/R   rad/R  |G|argmin   "
          "NNdist[m]   |G|/s^2argmin  max|G|/s^2  max|u|rec  p99.9|u|")
    rows = []
    for st in STEPS:
        d = corpse[st]
        s = d["sigma"]
        i = int(np.argmin(s))
        gm = np.linalg.norm(d["G"], axis=1)
        umag = np.linalg.norm(d["u"], axis=1)
        tree = cKDTree(d["X"])
        nnd = tree.query(d["X"][i], k=2)[0][1]
        gos2 = gm / s ** 2
        a1, r1 = cyl(d["X"][i:i + 1])
        rows.append((st, len(s), s[i], s[i] / SIGMA_SHED_CORPSE, a1[0], r1[0],
                     gm[i], nnd, gos2[i], gos2.max(), umag.max(),
                     np.quantile(umag, 0.999)))
        print("%d %7d  %.4e  %.5f  %5.2f  %5.2f  %.3e  %.4e  %.4e  %.4e  "
              "%8.1f  %8.1f" % rows[-1])
    arr = np.array(rows)

    # ignition core: per-step argmax |G|/sigma^2 particle
    print("\n-- ignition core: argmax |G|/sigma^2 per step --")
    print("step   Gos2max     sigma[m]    sig/shed  |G|        ax/R  rad/R  "
          "NNdist[m]   |u|rec")
    irows = []
    for st in STEPS:
        d = corpse[st]
        gm = np.linalg.norm(d["G"], axis=1)
        gos2 = gm / d["sigma"] ** 2
        i = int(np.argmax(gos2))
        tree = cKDTree(d["X"])
        nnd = tree.query(d["X"][i], k=2)[0][1]
        a1, r1 = cyl(d["X"][i:i + 1])
        urec = np.linalg.norm(d["u"][i])
        irows.append((st, gos2[i], d["sigma"][i],
                      d["sigma"][i] / SIGMA_SHED_CORPSE, gm[i], a1[0], r1[0],
                      nnd, urec))
        print("%d  %.4e  %.4e  %.5f  %.3e  %5.2f  %5.2f  %.4e  %8.1f"
              % irows[-1])
    np.savetxt(os.path.join(OUTDIR, "ignition_core.csv"), np.array(irows),
               delimiter=",", header="step,Gos2max,sigma,sig_rel,Gmag,"
               "ax_over_R,rad_over_R,nn_dist,u_rec", comments="")

    np.savetxt(os.path.join(OUTDIR, "death_trajectory.csv"), arr,
               delimiter=",", header="step,N,min_sigma,min_sig_rel,ax_over_R,"
               "rad_over_R,G_argmin,nn_dist,Gos2_argmin,max_Gos2,max_u,"
               "p999_u", comments="")

    fig, axs = plt.subplots(1, 3, figsize=(13, 3.6))
    axs[0].semilogy(arr[:, 0], arr[:, 2], "o-")
    axs[0].set_title("global min sigma [m]")
    axs[1].semilogy(arr[:, 0], arr[:, 9], "o-")
    axs[1].set_title("field max |G|/sigma^2")
    axs[2].semilogy(arr[:, 0], arr[:, 10], "o-", label="max")
    axs[2].semilogy(arr[:, 0], arr[:, 11], "s--", label="p99.9")
    axs[2].set_title("recorded |u| [m/s]")
    axs[2].legend()
    for a in axs:
        a.set_xlabel("step")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "death_trajectory.png"), dpi=140)
    plt.close(fig)

    # ---- Analysis 3: merge candidacy --------------------------------------
    d = corpse[1041]
    tree = cKDTree(d["X"])
    idx = np.where(m25)[0]
    dd, _ = tree.query(d["X"][idx], k=2)
    nnd = dd[:, 1]
    elig_abs = nnd < R_MERGE
    elig_2s = nnd < 2 * d["sigma"][idx]
    print("\n== Analysis 3: merge candidacy (sig/shed<0.25 set, N=%d, "
          "r_merge=%.4e m) ==" % (len(idx), R_MERGE))
    print("NN dist quantiles [m]: p5 %.3e  p25 %.3e  p50 %.3e  p75 %.3e  "
          "p95 %.3e" % tuple(q(nnd, p) for p in (.05, .25, .5, .75, .95)))
    print("eligible by NN<r_merge(abs): %d/%d = %.4f" %
          (elig_abs.sum(), len(idx), elig_abs.mean()))
    print("eligible by NN<2*sigma_local: %d/%d = %.4f" %
          (elig_2s.sum(), len(idx), elig_2s.mean()))
    np.savetxt(os.path.join(OUTDIR, "merge_candidacy_1041.csv"),
               np.column_stack([idx, d["sigma"][idx], nnd,
                                elig_abs, elig_2s]),
               delimiter=",", header="index,sigma,nn_dist,elig_abs,elig_2sig",
               comments="")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(np.log10(nnd), bins=60)
    ax.axvline(np.log10(R_MERGE), color="r", label="r_merge abs")
    ax.axvline(np.log10(2 * np.median(d["sigma"][idx])), color="g",
               label="2*med sigma")
    ax.set_xlabel("log10 NN dist [m] (collapsing set)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "merge_candidacy_hist.png"), dpi=140)
    plt.close(fig)

    # ---- Analysis 4 baseline: Z from recorded velocity_gradient ------------
    print("\n== Analysis 4 (recorded-J baseline; kernel recompute in "
          "p018_sigma_forensics_strain.jl) ==")
    for label, dat, sshed in (("corpse1041", corpse[1041],
                               SIGMA_SHED_CORPSE),
                              ("healthy719", healthy, SIGMA_SHED_HEALTHY)):
        Jnorm = np.abs(dat["J"]).max()
        if Jnorm == 0:
            print("%s: recorded velocity_gradient is all-zero -> recorded-J "
                  "baseline UNAVAILABLE" % label)
            continue
        Z = zstrain(dat["G"], dat["J"])
        srel = dat["sigma"] / sshed
        ax, rad = cyl(dat["X"])
        subsets = {"all": np.ones(len(Z), bool),
                   "aged column ax/R>1.3": ax > 1.3,
                   "collapsing sig/shed<0.25": srel < 0.25}
        for nm, m in subsets.items():
            if m.sum() == 0:
                continue
            dz = DT * Z[m]
            print("%s / %-26s N=%6d  dt*Z: med %+0.3e  p95 %+0.3e  "
                  "p5 %+0.3e" % (label, nm, m.sum(), q(dz, .5), q(dz, .95),
                                 q(dz, .05)))
        np.savetxt(os.path.join(OUTDIR, "recordedJ_dtZ_%s.csv" % label),
                   np.column_stack([ax, rad, srel, DT * Z]), delimiter=",",
                   header="ax_over_R,rad_over_R,sig_rel,dtZ", comments="")

    # ---- dump collapsing-subset targets for the Julia strain script --------
    sub = idx if len(idx) <= 3000 else np.random.default_rng(0).choice(
        idx, 3000, replace=False)
    np.savetxt(os.path.join(OUTDIR, "strain_targets_corpse1041.csv"),
               np.column_stack([sub, corpse[1041]["X"][sub],
                                corpse[1041]["G"][sub],
                                corpse[1041]["sigma"][sub]]),
               delimiter=",", header="index,x,y,z,Gx,Gy,Gz,sigma",
               comments="")
    print("\nWrote strain targets (N=%d) for Julia companion." % len(sub))


if __name__ == "__main__":
    main()
