#!/usr/bin/env python3
"""
render_dump.py -- render a relaxed minimum from your own runs, in the same style.

    python3 render_dump.py dumps/vw1_14_c3.lammps --cfg cfg/vw1_14.lmp
    python3 render_dump.py dumps/def_06_c1.lammps --cfg cfg/def_06.lmp --R 6
    python3 render_dump.py dumps/vw1_14_c3.lammps --full          # whole cell

--cfg   the cfg file the run included: its vacated sites are drawn as red
        markers and the close-up is centred on the defect (positions follow
        the dump's box, so relaxation of the cell is accounted for)
--R     close-up radius (A); default 5.5
--out   png name (default: dump name + .png)

Needs: pip install ovito (+ a Vulkan driver; on a headless Linux box
apt install mesa-vulkan-drivers and export
VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/lvp_icd.json). On a laptop with a
GPU it works as is.  Uses renders.py from the same folder.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import renders as R  # noqa: E402


def read_dump(fn):
    lines = open(fn).read().splitlines()
    nat = int(lines[lines.index("ITEM: NUMBER OF ATOMS") + 1])
    j = [k for k, l in enumerate(lines) if l.startswith("ITEM: BOX BOUNDS")][0]
    lo, L = [], []
    for k in range(3):
        a, b = map(float, lines[j + 1 + k].split()[:2])
        lo.append(a)
        L.append(b - a)
    h = [k for k, l in enumerate(lines) if l.startswith("ITEM: ATOMS")][0]
    cols = lines[h].split()[2:]
    d = np.array([list(map(float, l.split())) for l in lines[h + 1:h + 1 + nat]])
    xyz = d[:, [cols.index("x"), cols.index("y"), cols.index("z")]] - np.array(lo)
    return xyz, d[:, cols.index("type")].astype(int), np.array(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--cfg")
    ap.add_argument("--R", type=float, default=5.5)
    ap.add_argument("--cell", nargs=3, type=int, default=[8, 5, 8])
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()
    P, T, L = read_dump(a.dump)
    out = a.out or os.path.splitext(os.path.basename(a.dump))[0] + ".png"
    if a.full:
        R.render(P, T, L, out, cam=(-1, -1.5, -0.9), ortho=False, size=(900, 720))
        print("wrote", out)
        return
    host = T != 3
    vac = []
    if a.cfg:
        v, add = R.parse_cfg(a.cfg)
        vac = [np.array(f) * L / np.array(a.cell) for f, _ in v]
    he = P[T == 3]
    pts = np.array(vac + list(he)) if (vac or len(he)) else P[:1]
    cen = pts[0] if not len(pts) else pts.mean(axis=0)
    extra = list(he) + vac
    et = [3] * len(he) + [4] * len(vac)
    R.closeup(P[host], T[host], L, cen, extra, et, out, R=a.R)
    print("wrote", out)


if __name__ == "__main__":
    main()
