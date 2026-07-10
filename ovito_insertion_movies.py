#!/usr/bin/env python3
"""
ovito_insertion_movies.py
=========================
GIFs + final stills + erosion analysis from the per-cycle dumps of the
continuous He-insertion campaign (`in_cont`). Each cycle writes one frame:

    dump.bulk-<i>.lammps   W + C only             -> "bulk" views
    dump.full-<i>.lammps   W + C + accumulated He  -> "bulk + He" views

OVITO loads a numbered sequence via a '*' wildcard (sorted numerically), so the
movies/analysis grow automatically as cycles finish.

SURFACE REGRESSION (atom-count method):
    regression = (atoms removed) / (N0 / H0) = dN * H0 / N0
where N0 and H0 are the initial solid-atom count and slab thickness, so N0/H0 is
the linear atom density (atoms per A of height across the cross-section).
Removing a fraction of the atoms thins the slab by that fraction of its
thickness. This is tied to material actually removed and is robust to thermal
surface roughness. The older geometric estimate (z-percentile drop) is kept as a
comparison column/curve.

Produces, alongside the dumps:
  * <base>_<view>.gif / <base>_<view>_final.png   movies + final stills
  * erosion.csv             cycle, atoms_removed, regression (atom & geometric),
                            surface_z, sputtered W, sputtered C
  * erosion_regression.png  regression vs. cycle (atom-count + geometric)
  * erosion_sputter.png     cumulative sputtered W & C vs. cycle

Colors:  W = gray,  C = blue,  He = purple.

----------------------------------------------------------------------------
RUN ON THE CLUSTER (headless), from the folder holding the dumps:
    conda activate ovito
    pip install pillow matplotlib
    export QT_QPA_PLATFORM=offscreen
    export PYTHONUNBUFFERED=1
    python ovito_insertion_movies.py [optional /path/to/dumps]
----------------------------------------------------------------------------
"""

import os
import re
import sys
import csv
import glob
import bisect
import tempfile
import shutil

import numpy as np
from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer, PythonViewportOverlay
from ovito.qt_compat import QtGui

if QtGui.QGuiApplication.instance() is None:
    _qt_app = QtGui.QGuiApplication(sys.argv)

try:
    from PIL import Image
except ImportError:
    sys.exit("Pillow is required to assemble GIFs.  Install it:  pip install pillow")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

# ===========================================================================
# CONFIG
# ===========================================================================
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()
INPUT_DIR = sys.argv[1] if len(sys.argv) > 1 else SCRIPT_DIR
OUTPUT_DIR = INPUT_DIR

GIF_SIZE = (900, 680)
STILL_SIZE = (1600, 1200)
FPS = 12
FRAME_STRIDE = 1
GIF_AO = False
ISO_DIR = (-1.0, -1.0, -0.6)

SHOW_COUNTER = True
SURFACE_PCT = 98.0                            # percentile z for the geometric estimate
WRITE_CSV = True
MAKE_PLOTS = True

SPUTTER_FILE = "sputter_summary.dat"
SPUTTER_COLS = (0, 1, 2)

TYPE_STYLE = {                                # type_id: (name, RGB 0..1, radius A)
    1: ("W", (0.55, 0.55, 0.57), 1.10),       # gray
    2: ("C", (0.20, 0.45, 0.85), 0.75),       # blue
    3: ("He", (0.62, 0.22, 0.80), 0.60),      # purple
}

JOBS = [
    {"label": "bulk (W + C)", "pattern": "dump.bulk-*.lammps", "base": "bulk",
     "views": {"iso": {"gif": True, "still": True},
               "side": {"gif": True, "still": True}}},
    {"label": "bulk + He", "pattern": "dump.full-*.lammps", "base": "full",
     "views": {"iso": {"gif": True, "still": True},
               "side": {"gif": True, "still": True}}},
]
# ===========================================================================


def cycle_of(path):
    m = re.search(r"-(\d+)\.lammps$", os.path.basename(path))
    return int(m.group(1)) if m else -1


def load_sputter(path, cols):
    if not os.path.exists(path):
        return None
    cc, wc, cco = cols
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            try:
                rows.append((int(float(parts[cc])), float(parts[wc]), float(parts[cco])))
            except (IndexError, ValueError):
                continue
    if not rows:
        return None
    rows.sort(key=lambda r: r[0])
    cyc_list, cumW, cumC = [], [], []
    sw = sc = 0.0
    last = None
    for cyc, w, c in rows:
        sw += w; sc += c
        if cyc == last:
            cumW[-1], cumC[-1] = sw, sc
        else:
            cyc_list.append(cyc); cumW.append(sw); cumC.append(sc); last = cyc
    return cyc_list, cumW, cumC


def cum_at(cum, cyc):
    if cum is None:
        return None
    cyc_list, cumW, cumC = cum
    idx = bisect.bisect_right(cyc_list, cyc) - 1
    if idx < 0:
        return (0, 0)
    return (int(round(cumW[idx])), int(round(cumC[idx])))


def set_style(pipeline):
    data = pipeline.compute(0)
    present = {t.id for t in data.particles.particle_types.types}
    ptypes = pipeline.source.data.particles_.particle_types_
    for tid, (name, color, radius) in TYPE_STYLE.items():
        if tid in present:
            t = ptypes.type_by_id_(tid)
            t.name, t.color, t.radius = name, color, radius


def solid_stats(pipeline, frame):
    """(n_solid, surface_z, z_min, z_max) over the W+C atoms of a frame."""
    data = pipeline.compute(frame)
    z = data.particles["Position"][:, 2]
    ptype = data.particles["Particle Type"][:]
    zs = z[(ptype == 1) | (ptype == 2)]
    if zs.size == 0:
        return 0, float("nan"), float("nan"), float("nan")
    return int(zs.size), float(np.percentile(zs, SURFACE_PCT)), float(zs.min()), float(zs.max())


# ---------------------------------------------------------------------------
# Erosion analysis (one pass): atom-count regression + geometric comparison
# ---------------------------------------------------------------------------
def pick_erosion_pattern():
    for pat in ("dump.bulk-*.lammps", "dump.full-*.lammps"):
        if glob.glob(os.path.join(INPUT_DIR, pat)):
            return pat
    return None


def analyze_erosion(cum):
    """Return (rows, reg_by_cycle).
    rows = [(cycle, atoms_removed, regression_atom, regression_geom, surface_z, W, C)].
    reg_by_cycle maps cycle -> atom-count regression (used for the GIF printout)."""
    pat = pick_erosion_pattern()
    if pat is None:
        print("  no dumps available for erosion analysis")
        return [], {}
    files = sorted(glob.glob(os.path.join(INPUT_DIR, pat)), key=cycle_of)
    cycles = [cycle_of(f) for f in files]
    pipeline = import_file(os.path.join(INPUT_DIR, pat))
    n = pipeline.num_frames

    N0, surf0, zmin0, zmax0 = solid_stats(pipeline, 0)
    H0 = zmax0 - zmin0
    lin_density = (N0 / H0) if (N0 and H0 > 0) else float("nan")   # atoms per A of height
    print(f"  baseline cycle {cycles[0]}: N0={N0} solid atoms, slab H0~{H0:.2f} A, "
          f"linear density ~{lin_density:.1f} atoms/A, surf0~{surf0:.2f} A")
    if not np.isfinite(lin_density):
        print("  WARNING: degenerate baseline; atom-count regression will be NaN.")

    rows, reg_by_cycle = [], {}
    for k in range(n):
        N, surf, _zmin, _zmax = solid_stats(pipeline, k)
        dN = N0 - N
        reg_atom = (dN / lin_density) if np.isfinite(lin_density) else float("nan")
        reg_geom = surf0 - surf
        ca = cum_at(cum, cycles[k])
        w = ca[0] if ca else ""
        c = ca[1] if ca else ""
        rows.append((cycles[k], dN, round(reg_atom, 3), round(reg_geom, 3),
                     round(surf, 3), w, c))
        reg_by_cycle[cycles[k]] = reg_atom
    last = rows[-1]
    print(f"  final cycle {last[0]}: atoms removed={last[1]}, "
          f"regression(atom-count)~{last[2]:+.2f} A, geometric~{last[3]:+.2f} A")
    return rows, reg_by_cycle


def write_csv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cycle", "atoms_removed", "regression_A", "regression_geom_A",
                    "surface_z_A", "sputtered_W", "sputtered_C"])
        w.writerows(rows)
    print("  wrote", path)


def make_plots(rows, outdir, has_sputter):
    if not HAVE_MPL:
        print("  matplotlib not installed - skipping plots (pip install matplotlib)")
        return
    cyc = [r[0] for r in rows]
    reg_atom = [r[2] for r in rows]
    reg_geom = [r[3] for r in rows]

    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.plot(cyc, reg_atom, "-", color="#b03030", lw=1.7, label="atom-count")
    ax.plot(cyc, reg_geom, "--", color="#888888", lw=1.4, label="geometric (z %ile)")
    ax.set_xlabel("insertion cycle")
    ax.set_ylabel("surface regression (\u00c5)")
    ax.set_title("Surface regression vs. cycle")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p1 = os.path.join(outdir, "erosion_regression.png")
    fig.savefig(p1, dpi=150); plt.close(fig)
    print("  wrote", p1)

    if has_sputter:
        W = [r[5] for r in rows]
        C = [r[6] for r in rows]
        fig, ax = plt.subplots(figsize=(7, 4.3))
        ax.plot(cyc, W, "-", color="#777777", lw=1.6, label="W")
        ax.plot(cyc, C, "-", color="#2858b4", lw=1.6, label="C")
        ax.set_xlabel("insertion cycle")
        ax.set_ylabel("cumulative sputtered atoms")
        ax.set_title("Cumulative sputtering vs. cycle")
        ax.grid(alpha=0.3); ax.legend()
        fig.tight_layout()
        p2 = os.path.join(outdir, "erosion_sputter.png")
        fig.savefig(p2, dpi=150); plt.close(fig)
        print("  wrote", p2)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def iso_viewport(size):
    vp = Viewport(type=Viewport.Type.Ortho)
    vp.camera_dir = ISO_DIR
    vp.zoom_all(size=size)
    vp.fov *= 1.10
    return vp


def side_viewport(size):
    vp = Viewport(type=Viewport.Type.Front)
    vp.zoom_all(size=size)
    vp.fov *= 1.05
    return vp


def make_viewport(view, size):
    return iso_viewport(size) if view == "iso" else side_viewport(size)


def make_renderer(ao):
    return TachyonRenderer(ambient_occlusion=ao, shadows=False,
                           direct_light_intensity=1.1)


def counter_overlay(frame_info):
    def render(args):
        info = frame_info.get(args.frame)
        if info is None:
            return
        cyc, cw, cc = info
        p = args.painter
        W, H = p.window().width(), p.window().height()
        fs = max(14, int(H * 0.028))
        font = p.font(); font.setPixelSize(fs); font.setBold(True); p.setFont(font)
        x, y, dy = int(0.03 * W), int(0.07 * H), int(fs * 1.25)
        p.setPen(QtGui.QPen(QtGui.QColor(30, 30, 30)))
        p.drawText(x, y, f"cycle {cyc}")
        if cw is not None:
            p.setPen(QtGui.QPen(QtGui.QColor(120, 120, 128)))
            p.drawText(x, y + dy, f"sputtered W: {cw}")
            p.setPen(QtGui.QPen(QtGui.QColor(40, 90, 180)))
            p.drawText(x, y + 2 * dy, f"sputtered C: {cc}")
    return render


def render_gif(pipeline, view, base, frame_info, cycles, reg_by_cycle):
    n = pipeline.num_frames
    vp = make_viewport(view, GIF_SIZE)
    if SHOW_COUNTER:
        vp.overlays.append(PythonViewportOverlay(function=counter_overlay(frame_info)))

    frames = list(range(0, n, FRAME_STRIDE))
    if frames[-1] != n - 1:
        frames.append(n - 1)
    tmp = tempfile.mkdtemp(prefix="ovito_frames_")
    paths = []
    for j, k in enumerate(frames):
        fp = os.path.join(tmp, f"f{j:05d}.png")
        vp.render_image(filename=fp, frame=k, size=GIF_SIZE,
                        background=(1, 1, 1), renderer=make_renderer(GIF_AO))
        paths.append(fp)
        reg = reg_by_cycle.get(cycles[k])
        if reg is None or not np.isfinite(reg):
            print(f"    [{view}] frame {j + 1}/{len(frames)}  cycle {cycles[k]}")
        else:
            print(f"    [{view}] frame {j + 1}/{len(frames)}  cycle {cycles[k]}: "
                  f"regression(atom-count) ~{reg:+.2f} A")

    out_gif = os.path.join(OUTPUT_DIR, f"{base}_{view}.gif")
    imgs = [Image.open(p).convert("RGB") for p in paths]
    imgs[0].save(out_gif, save_all=True, append_images=imgs[1:],
                 duration=max(1, int(1000 / FPS)), loop=0, optimize=True)
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"    wrote {out_gif}  ({len(frames)} frames)")


def render_still(pipeline, view, base, frame_info):
    n = pipeline.num_frames
    vp = make_viewport(view, STILL_SIZE)
    if SHOW_COUNTER:
        vp.overlays.append(PythonViewportOverlay(function=counter_overlay(frame_info)))
    out_png = os.path.join(OUTPUT_DIR, f"{base}_{view}_final.png")
    vp.render_image(filename=out_png, frame=n - 1, size=STILL_SIZE,
                    background=(1, 1, 1), renderer=make_renderer(True))
    print(f"    wrote {out_png}  (final cycle {frame_info.get(n - 1, ('?',))[0]})")


def make_job(job, cum, reg_by_cycle):
    pat = os.path.join(INPUT_DIR, job["pattern"])
    files = sorted(glob.glob(pat), key=cycle_of)
    if not files:
        print(f"  [skip] no files match {job['pattern']} in {INPUT_DIR}")
        return
    cycles = [cycle_of(f) for f in files]
    pipeline = import_file(pat)
    n = pipeline.num_frames
    print(f"  {job['label']}: {n} cycles  ({job['pattern']})")
    set_style(pipeline)
    pipeline.add_to_scene()

    frame_info = {}
    for k, cyc in enumerate(cycles):
        ca = cum_at(cum, cyc)
        frame_info[k] = (cyc, ca[0] if ca else None, ca[1] if ca else None)

    try:
        for view, out in job["views"].items():
            if out.get("gif"):
                print(f"  -- {view} GIF --")
                render_gif(pipeline, view, job["base"], frame_info, cycles, reg_by_cycle)
            if out.get("still"):
                print(f"  -- {view} final still --")
                render_still(pipeline, view, job["base"], frame_info)
    finally:
        pipeline.remove_from_scene()


def main():
    print(f"Reading dumps from: {INPUT_DIR}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    cum = load_sputter(os.path.join(INPUT_DIR, SPUTTER_FILE), SPUTTER_COLS)
    if cum is None:
        print(f"  NOTE: {SPUTTER_FILE} not found/empty - sputter counters/plot disabled.")

    print("=== erosion analysis (atom-count method) ===")
    rows, reg_by_cycle = analyze_erosion(cum)
    if rows:
        if WRITE_CSV:
            write_csv(rows, os.path.join(OUTPUT_DIR, "erosion.csv"))
        if MAKE_PLOTS:
            make_plots(rows, OUTPUT_DIR, has_sputter=(cum is not None))

    for job in JOBS:
        print(f"=== {job['label']} ===")
        make_job(job, cum, reg_by_cycle)

    print("Done. Output in:", os.path.abspath(OUTPUT_DIR))


if __name__ == "__main__":
    main()
