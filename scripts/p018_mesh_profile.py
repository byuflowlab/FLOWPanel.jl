#!/usr/bin/env python3
"""Extract and compare chord(r/R) and twist(r/R) from DJI9443 gmsh meshes.

Why this exists: the hub-variant mesh
`dji9443_20260731_45_185_capped_captess4_hub0p15.msh` was accepted by a gate
that checked watertightness and element count. It passed both while carrying a
root chord 5.5x too large and the wrong twist, and produced CT = -0.474 at the
first timestep. Watertightness cannot see a geometry error; a chord/twist
profile can.

Use as the acceptance gate for any `--hub-r-over-r` mesh:

    p018_mesh_profile.py show  MESH
    p018_mesh_profile.py check HUB_MESH --stock STOCK_MESH --hub 0.15

`check` requires that, outboard of the hub, the hub blade reproduces the stock
blade (it must differ ONLY by removing inboard area -- otherwise a hub-vs-stock
delta is not interpretable), and that the root section matches the stock
distribution evaluated at the hub radius rather than the stock r/R=0.01
section.
"""

import argparse
import sys

import numpy as np

RING = 184          # nodes per OpenVSP loft ring in the *185-chordwise* family.
                    # = CHORDWISE - 1. Other resolutions MUST pass --ring, or the
                    # fixed-stride slicing silently mixes rings and the comparison
                    # is meaningless (e.g. an 81-chordwise mesh needs --ring 80).


def read_nodes(path):
    """Node coordinates in file order (gmsh 4.1 ASCII)."""
    lines = open(path).read().split("\n")
    i = lines.index("$Nodes")
    nblk = int(lines[i + 1].split()[0])
    j = i + 2
    pts = []
    for _ in range(nblk):
        n = int(lines[j].split()[3])
        j += 1 + n                      # block header + node tags
        for k in range(n):
            pts.append([float(x) for x in lines[j + k].split()])
        j += n
    return np.array(pts)


def profile(path, ring=RING):
    """[(r/R, chord, twist_deg)] for one blade, in file (radial) order.

    Chord is the maximum pairwise distance within a ring (the LE-TE line);
    twist is that line's angle out of the rotor disk plane. The rotor axis is
    x; the radial direction is y.
    """
    P = read_nodes(path)
    if len(P) % ring:
        print(f"warning: {len(P)} nodes is not a multiple of {ring}", file=sys.stderr)
    R = np.abs(P[:, 1]).max()
    out = []
    for i in range(len(P) // ring):
        r = P[i * ring:(i + 1) * ring]
        # Both blades share one node list, so one fixed-stride group straddles
        # the blade-1 tip / blade-2 root boundary. Such a group is not a loft
        # ring and its "chord" is ~R; require the whole group on one blade.
        if not np.all(r[:, 1] > 0):
            continue
        d = np.linalg.norm(r[:, None, :] - r[None, :, :], axis=-1)
        a, b = np.unravel_index(np.argmax(d), d.shape)
        v = r[b] - r[a]
        out.append((float(np.mean(r[:, 1]) / R), float(d[a, b]),
                    float(np.degrees(np.arctan2(abs(v[0]), abs(v[2]))))))
    return out


def interp(prof, x, col):
    r = np.array([p[0] for p in prof])
    y = np.array([p[col] for p in prof])
    o = np.argsort(r)
    return float(np.interp(x, r[o], y[o]))


def cmd_show(a):
    prof = profile(a.mesh, a.ring)
    print(f"{'ring':>4} {'r/R':>8} {'chord':>9} {'twist':>8}")
    for i, (r, c, t) in enumerate(prof):
        print(f"{i:4d} {r:8.4f} {c:9.4f} {t:8.2f}")
    print(f"\nmax chord {max(p[1] for p in prof):.4f}  rings {len(prof)}")


def cmd_check(a):
    hub, stock = profile(a.mesh, a.ring), profile(a.stock, a.ring)
    ok = True

    # 1. No ring may fold inboard of the root station.
    # Root-cap rings legitimately sit AT the root station, so compare with a
    # tolerance; a genuine fold (the pre-fix mesh had a ring 0.004 inboard of
    # the root) is far outside it.
    r0 = hub[0][0]
    folded = [(i, r) for i, (r, _, _) in enumerate(hub) if r < r0 - 1e-3]
    if folded:
        ok = False
        print(f"FAIL fold: {len(folded)} ring(s) inboard of the root r/R={r0:.4f}: "
              + ", ".join(f"#{i}@{r:.4f}" for i, r in folded[:5]))
    else:
        print(f"pass fold: no ring inboard of root r/R={r0:.4f}")

    # 2. No chord may exceed the stock blade's maximum (catches the balloon).
    cmax_stock = max(p[1] for p in stock)
    cmax_hub = max(p[1] for p in hub)
    if cmax_hub > cmax_stock * 1.02:
        ok = False
        print(f"FAIL chord bound: hub max chord {cmax_hub:.4f} exceeds stock "
              f"{cmax_stock:.4f} by {100*(cmax_hub/cmax_stock-1):.1f}%")
    else:
        print(f"pass chord bound: hub max {cmax_hub:.4f} vs stock {cmax_stock:.4f}")

    # 3. Outboard of the hub (plus margin) the hub blade must BE the stock blade.
    lo = a.hub + a.margin
    band = [(r, c, t) for r, c, t in hub if r >= lo]
    if not band:
        ok = False
        print(f"FAIL outboard: no rings outboard of r/R={lo:.3f}")
    else:
        dc = max(abs(c - interp(stock, r, 1)) / cmax_stock for r, c, _ in band)
        dt = max(abs(t - interp(stock, r, 2)) for r, _, t in band)
        c_ok, t_ok = dc <= a.chord_tol, dt <= a.twist_tol
        ok &= c_ok and t_ok
        print(f"{'pass' if c_ok else 'FAIL'} outboard chord: max dev "
              f"{100*dc:.3f}% of max chord (tol {100*a.chord_tol:.3f}%), "
              f"{len(band)} rings r/R>={lo:.3f}")
        print(f"{'pass' if t_ok else 'FAIL'} outboard twist: max dev "
              f"{dt:.3f} deg (tol {a.twist_tol:.3f})")

    # 4. The blade must start from the stock distribution evaluated at the hub,
    #    NOT from the stock r/R=0.01 section the pre-fix generator left there.
    #
    #    Compare the first INTERIOR ring, not the root cap rings: a cap ring
    #    measures the true 2D chord (stock's own cap at r/R=0.01 reads exactly
    #    0.5400 = PCurve x R), while interior rings measure a 3D LE-TE distance
    #    that runs ~4% high because of sweep/rake. Checking a cap ring against
    #    an interpolation of interior rings compares unlike quantities and
    #    fails by that systematic offset alone.
    r0 = hub[0][0]
    interior = [(r, c, t) for r, c, t in hub if r > r0 + 1e-3]
    if not interior:
        ok = False
        print("FAIL root: no interior ring outboard of the root cap")
    else:
        r_i, c_got, t_got = interior[0]
        c_exp, t_exp = interp(stock, r_i, 1), interp(stock, r_i, 2)
        c_ok = abs(c_got - c_exp) / cmax_stock <= a.root_tol
        t_ok = abs(t_got - t_exp) <= a.twist_tol * 2
        ok &= c_ok and t_ok
        print(f"{'pass' if c_ok else 'FAIL'} root chord: first interior ring at "
              f"r/R={r_i:.4f} is {c_got:.4f} vs stock {c_exp:.4f}")
        print(f"{'pass' if t_ok else 'FAIL'} root twist: {t_got:.2f} deg vs stock "
              f"{t_exp:.2f} deg")

    # 5. Explicit bug-signature check: the root must NOT still be carrying the
    #    stock r/R=0.01 section (chord 0.5400 / twist 13.19 at the hub radius).
    stale = (abs(hub[0][1] - stock[0][1]) / cmax_stock < 0.01
             and abs(hub[0][2] - stock[0][2]) < 0.5)
    ok &= not stale
    print(f"{'FAIL' if stale else 'pass'} root not stale: root section "
          f"{hub[0][1]:.4f}/{hub[0][2]:.2f}deg vs stock r/R={stock[0][0]:.4f} "
          f"section {stock[0][1]:.4f}/{stock[0][2]:.2f}deg")

    print("\n" + ("ALL CHECKS PASS" if ok else "GATE FAILED"))
    return 0 if ok else 1


def main():
    p = argparse.ArgumentParser()
    s = p.add_subparsers(dest="cmd", required=True)
    q = s.add_parser("show"); q.add_argument("mesh")
    q.add_argument("--ring", type=int, default=RING)
    q = s.add_parser("check")
    q.add_argument("mesh")
    q.add_argument("--stock", required=True)
    q.add_argument("--hub", type=float, required=True)
    q.add_argument("--margin", type=float, default=0.05)
    q.add_argument("--chord-tol", type=float, default=0.01)   # of max chord
    # 0.5 deg absorbs interpolation noise between two meshes whose radial
    # stations differ, while still catching the 6.9 deg root error that the
    # pre-fix generator produced.
    q.add_argument("--twist-tol", type=float, default=0.5)    # degrees
    q.add_argument("--root-tol", type=float, default=0.03)    # of max chord
    q.add_argument("--ring", type=int, default=RING,
                   help="nodes per loft ring = CHORDWISE-1 (184 for the 185 family)")
    a = p.parse_args()
    return cmd_show(a) if a.cmd == "show" else cmd_check(a)


if __name__ == "__main__":
    sys.exit(main())
