#!/usr/bin/env python3
"""
rand_orientation_gif.py - animated GIF of the random surfaces, one frame per surface.

Each frame draws random surface s the way it sits in the ORIGINAL (unrotated)
lattice: the crystal block in its own frame, with the plane that the step-1
deck turns into the irradiated surface.  The plane's Miller label sits in the
top-left corner.  Optionally the right half shows what the He sees: the top
view of the same surface after rotation.

Angles come from rand_angles_lammps.csv (the LAMMPS random stream of the
decks), so frame s is surface s of the runs; pass --surfaces-dat to read them
from the step-1 log instead.

Edit the options below, or override them on the command line:
    python rand_orientation_gif.py                       # W, surfaces 1-50
    python rand_orientation_gif.py --mat WC --first 1 --last 20 --ms 800
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

import rand_geom as g
import rand_surface_views as v

# ============================ OPTIONS =====================================
MATERIAL = "W"            # "W" or "WC"
FRAME_FIRST = 1           # first random surface s
FRAME_LAST = 50           # last random surface s (inclusive)
FRAME_STEP = 1            # every n-th surface
MS_PER_FRAME = 1200       # display time per frame (ms); 1200 ms = 0.83 frames/s
HOLD_FIRST_MS = 2500      # longer pause on the first frame (ms)
REFERENCE_FRAMES = True   # start with the low-index planes: W (001),(011),(111); WC basal, prism, a-plane
SHOW_TOPVIEW = True       # right panel: top view of the rotated surface (what the He sees)
LOOP = 0                  # 0 = loop forever
DPI = 90                  # GIF resolution (90 dpi x 12 in = 1080 px wide)
FIGSIZE = (12.0, 5.6)     # inches; use (6.4, 5.6) with SHOW_TOPVIEW = False
SURFACES_DAT = None       # e.g. "rand_WC_surfaces.dat" to use the angles a run actually logged
OUT_DIR = "figs"
# ===========================================================================

# low-index references; small origin shifts keep lattice planes off the cut plane z = 0
# (a plane of atoms exactly at z = 0 is split by round-off; the decks' random origin avoids this)
REFS = {
    "W": [("(001) reference", (0, 0, 0, 0.05, 0.05, 0.05)), ("(011) reference", (0, 45, 0, 0.05, 0.05, 0.05)),
          ("(111) reference", (0, 54.7356, 45, 0.05, 0.05, 0.05))],
    "WC": [("basal reference, C-terminated", (0, 0, 0, 0.05, 0.05, 0.45)), ("prism reference", (0, 90, 0, 0.05, 0.05, 0.05)),
           ("a-plane reference", (0, 90, 90, 0.05, 0.05, 0.05))],
}


def render_frame(mat, title, angles, show_top):
    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    if show_top:
        ax3 = fig.add_axes([0.0, 0.0, 0.5, 0.86], projection="3d")
        ax2 = fig.add_axes([0.53, 0.04, 0.44, 0.80])
    else:
        ax3 = fig.add_axes([0.0, 0.0, 1.0, 0.86], projection="3d")
    n = v.draw_block(ax3, mat, angles, show_regular=True)
    lab, dev = g.plane_label(mat, n)
    # ---- label, top-left corner
    fig.text(0.012, 0.975, "%s  %s%s" % (mat, "" if dev < 0.5 else "≈", lab), fontsize=20, weight="bold",
             va="top", ha="left", color="#1d3f8f",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#1d3f8f", lw=1.2))
    fig.text(0.012, 0.865, "%s\nn = (%.2f, %.2f, %.2f)   β = %.0f°  γ = %.0f°"
             % (title, n[0], n[1], n[2], angles[1], angles[2]), fontsize=10.5, va="top", color="#333333")
    fig.text(0.25 if show_top else 0.5, 0.015, "unrotated crystal: blue plane = surface cut, orange = regular top face",
             ha="center", fontsize=9, color="#555555")
    if show_top:
        v.draw_topview(ax2, mat, angles, "top view after rotation (what the He sees)")
    fig.canvas.draw()
    img = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3])
    plt.close(fig)
    return img


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mat", default=MATERIAL, choices=["W", "WC"])
    ap.add_argument("--first", type=int, default=FRAME_FIRST)
    ap.add_argument("--last", type=int, default=FRAME_LAST)
    ap.add_argument("--step", type=int, default=FRAME_STEP)
    ap.add_argument("--ms", type=int, default=MS_PER_FRAME, help="ms per frame")
    ap.add_argument("--no-refs", action="store_true", help="skip the low-index reference frames")
    ap.add_argument("--no-top", action="store_true", help="3D block only")
    ap.add_argument("--surfaces-dat", default=SURFACES_DAT)
    ap.add_argument("--out", default=OUT_DIR)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    show_top = SHOW_TOPVIEW and not a.no_top

    frames, durs, src = [], [], set()
    if REFERENCE_FRAMES and not a.no_refs:
        for title, ang in REFS[a.mat]:
            frames.append(render_frame(a.mat, title, ang, show_top)); durs.append(a.ms)
    for s in range(a.first, a.last + 1, a.step):
        ang, sname = g.get_angles(s, a.surfaces_dat)
        src.add(sname)
        frames.append(render_frame(a.mat, "random surface %d" % s, ang, show_top)); durs.append(a.ms)
        print("frame s=%d (%s)" % (s, sname))
    durs[0] = max(durs[0], HOLD_FIRST_MS)
    # one shared 256-colour palette for every frame (no colour flicker between frames)
    w, h = frames[0].size
    sample = Image.new("RGB", (w, h * min(len(frames), 6)))
    for k, f in enumerate(frames[:: max(1, len(frames) // 6)][:6]):
        sample.paste(f, (0, k * h))
    ref = sample.quantize(colors=256, method=Image.MEDIANCUT)
    pal = [f.quantize(palette=ref, dither=Image.Dither.NONE) for f in frames]
    path = os.path.join(a.out, "rand_orientations_%s_s%d-%d.gif" % (a.mat, a.first, a.last))
    pal[0].save(path, save_all=True, append_images=pal[1:], duration=durs, loop=LOOP, optimize=True, disposal=2)
    print("wrote %s  (%d frames, %.1f MB; angles: %s)"
          % (path, len(pal), os.path.getsize(path) / 1e6, ", ".join(sorted(src))))


if __name__ == "__main__":
    main()
