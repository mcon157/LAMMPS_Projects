#!/usr/bin/env python3
"""
cim_replica_schematic.py - how the CIM replicas branch from the main run.

Draws the main CIM trajectory (cycles 0 -> end), the checkpoints written every
500 cycles, and N independent 100-cycle replicas fanning out of the 0, 500 and
1000-cycle checkpoints, each giving one yield sample; mean +- 95% CI per branch.

Usage:  python cim_replica_schematic.py [--nrep 10] [--ncyc 100] [--out cim_replica_branches.png]
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAIN_END = 1550          # cycles of the main run drawn
BRANCHES = (0, 500, 1000)
PER_CYCLE = 2.80e16      # He m^-2 per cycle (WC basal impact area)
C_MAIN, C_REP, C_CI, C_CK = "#4b4b4b", "#e2711d", "#7d0a12", "#2b7bba"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nrep", type=int, default=10)
    ap.add_argument("--ncyc", type=int, default=100)
    ap.add_argument("--out", default="cim_replica_branches.png")
    a = ap.parse_args()
    rng = np.random.default_rng(3)

    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    # main run
    ax.annotate("", xy=(MAIN_END + 60, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", lw=6, color="#bfbfbf", shrinkA=0, shrinkB=0))
    ax.text(MAIN_END + 70, 0, "main CIM run\n(one trajectory)", va="center", fontsize=10, color=C_MAIN)
    for ck in range(0, MAIN_END + 1, 500):
        used = ck in BRANCHES
        ax.plot(ck, 0, "D", ms=11, color=C_CK if used else "#bfbfbf", mec="k", mew=0.6, zorder=5)
        ax.text(ck, -0.17, "he_cycle_%d.restart" % ck if ck else "wche0.restart", ha="center", va="top",
                fontsize=8.5, color=C_CK if used else "#888888")
    # replicas
    for b in BRANCHES:
        ys = np.linspace(0.18, 1.0, a.nrep)
        for k, y in enumerate(ys):
            x1 = b + a.ncyc * 2.6            # stretched so the fan is readable
            xs = np.linspace(b, x1, 30)
            yy = y * np.sin(np.linspace(0, np.pi / 2, 30))
            taper = (xs - b) / (x1 - b)
            ax.plot(xs, yy + 0.008 * rng.standard_normal(30).cumsum() * taper, color=C_REP, lw=1.1, alpha=0.85)
        ax.text(b + a.ncyc * 2.6 + 12, 1.02, "rep 1…%d" % a.nrep, fontsize=8.5, color=C_REP, va="bottom")
        # per-branch yield estimate
        xm = b + a.ncyc * 2.6 + 55
        ax.errorbar(xm, 0.58, yerr=0.2, fmt="o", color=C_CI, capsize=5, ms=6, lw=1.6)
        ax.text(xm + 18, 0.58, r"$\bar Y \pm t\,s/\sqrt{N}$", fontsize=10, color=C_CI, va="center")
    ax.text(30, 1.32, "each replica: same restart, seed = 123 + 237 i + 100003·rep  →  %d cycles  →  "
            r"$Y_r$ = atoms removed / %d" % (a.ncyc, a.ncyc), fontsize=10, color="#333333")
    # fluence axis
    ax.set_xlim(-60, MAIN_END + 260)
    ax.set_ylim(-0.42, 1.45)
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.set_xticks(BRANCHES + (1500,))
    ax.set_xlabel("cycle of the main run  (replica fans drawn %.1f× longer than 100 cycles for clarity)" % 2.6,
                  fontsize=9.5)
    sec = ax.secondary_xaxis("top", functions=(lambda c: c * PER_CYCLE / 1e19, lambda f: f * 1e19 / PER_CYCLE))
    sec.set_xlabel(r"fluence ($10^{19}$ He m$^{-2}$, WC basal)", fontsize=9.5)
    sec.set_xticks([0, 1.4, 2.8, 4.2])
    for s_ in ("left", "right"):
        sec.spines[s_].set_visible(False)
    fig.tight_layout()
    fig.savefig(a.out, dpi=200)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
