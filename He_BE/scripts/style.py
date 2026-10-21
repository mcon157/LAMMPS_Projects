"""Shared colours and drawing helpers for the schematics and GIFs.

Atom colours follow the render convention used throughout the project:
W gray, C blue, He yellow, vacancy red.  Line/accent colours are the same
fixed order as the paper figures.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch

W_C, C_C, HE_C, VAC_C = "#8c8c8c", "#2a78d6", "#f2c12e", "#d62f2f"
W_EDGE, C_EDGE, HE_EDGE = "#4d4d4d", "#174a87", "#9a7400"
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID, PAPER = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff"
SOFT = {"in": "#e8f0fb", "proc": "#fdf0e7", "ana": "#e6f6ef", "out": "#fdf6e1", "lim": "#fbecef"}
EDGE = {"in": S[0], "proc": S[1], "ana": S[2], "out": "#b07800", "lim": "#c0485f"}

plt.rcParams.update({
    "font.size": 10, "font.family": "DejaVu Sans", "axes.edgecolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "savefig.bbox": "tight",
    "savefig.facecolor": PAPER, "figure.facecolor": PAPER,
})

# ---------------------------------------------------------------- WC lattice
A, C_LAT = 2.91862, 2.81498
B = np.sqrt(3.0) * A
# orthogonal cell: W at z=0, C at z=c/2; the empty trigonal prisms (channels)
# sit at the same height as C
W_BASIS = [(0.0, 0.0), (0.5, 0.5)]
C_BASIS = [(0.0, 2 / 3), (0.5, 1 / 6)]
CH_BASIS = [(0.5, 5 / 6), (0.0, 1 / 3)]


def wc_sites(nx, ny, basis):
    out = []
    for i in range(-1, nx + 1):
        for j in range(-1, ny + 1):
            for fx, fy in basis:
                out.append(((i + fx) * A, (j + fy) * B))
    return np.array(out)


def atom(ax, xy, kind, r=None, z=3, alpha=1.0, lw=0.8):
    col, edge, rr = {"W": (W_C, W_EDGE, 0.62), "C": (C_C, C_EDGE, 0.38),
                     "He": (HE_C, HE_EDGE, 0.42)}[kind]
    ax.add_patch(Circle(xy, r or rr, fc=col, ec=edge, lw=lw, zorder=z, alpha=alpha))


def vacancy(ax, xy, r=0.62, z=4, lw=1.6):
    ax.add_patch(Circle(xy, r, fc="none", ec=VAC_C, lw=lw, ls=(0, (2.5, 1.5)), zorder=z))


def box(ax, x, y, w, h, text, kind="in", fs=9, weight="normal", ha="center",
        tc=INK, lw=1.2, alpha=1.0, pad=0.02):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=%g,rounding_size=0.02" % pad,
                       fc=SOFT[kind], ec=EDGE[kind], lw=lw, alpha=alpha, zorder=2)
    ax.add_patch(p)
    tx = x + w / 2 if ha == "center" else x + 0.012
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fs, color=tc,
            weight=weight, zorder=3, linespacing=1.35)
    return p


def arrow(ax, p0, p1, color=INK2, lw=1.4, ms=12, style="-|>", ls="-", z=1, rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, lw=lw,
                                 color=color, ls=ls, zorder=z,
                                 connectionstyle="arc3,rad=%g" % rad))


def legend_atoms(ax, x, y, dx=0.16, fs=8.5, items=("W", "C", "He", "vac"), r=0.018):
    """Atom legend in axes-fraction coordinates (ax must have aspect equal in data
    units, so draw markers with scatter in axes coords instead)."""
    sty = {"W": (W_C, W_EDGE, "-"), "C": (C_C, C_EDGE, "-"), "He": (HE_C, HE_EDGE, "-"),
           "vac": ("none", VAC_C, "--")}
    lab = {"W": "W", "C": "C", "He": "He", "vac": "vacancy"}
    for k, it in enumerate(items):
        fc, ec, ls = sty[it]
        ax.scatter([x + k * dx], [y], s=110, c=[fc] if fc != "none" else "none",
                   edgecolors=ec, linewidths=1.5 if it == "vac" else 0.8,
                   linestyle=ls, transform=ax.transAxes, zorder=10, clip_on=False)
        ax.text(x + k * dx + 0.035, y, lab[it], transform=ax.transAxes, va="center",
                fontsize=fs, color=INK)


def save(fig, path, dpi=200):
    fig.savefig(path + ".png", dpi=dpi)
    fig.savefig(path + ".pdf")
    plt.close(fig)
    print("wrote", path + ".{png,pdf}")


def fig_to_image(fig, dpi=100):
    """Render a figure to a PIL image (for GIF frames)."""
    from PIL import Image
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    return Image.fromarray(buf[:, :, :3].copy())


def write_gif(frames, path, ms=90, hold_last=1500, hold=None):
    """frames: list of PIL images.  hold: {index: ms} for extra-long frames."""
    from PIL import Image
    durs = [ms] * len(frames)
    for i, d in (hold or {}).items():
        durs[i] = d
    durs[-1] = max(durs[-1], hold_last)
    pal = [f.convert("P", palette=Image.ADAPTIVE, colors=128) for f in frames]
    pal[0].save(path, save_all=True, append_images=pal[1:], duration=durs, loop=0,
                optimize=True, disposal=2)
    print("wrote %s (%d frames)" % (path, len(frames)))
