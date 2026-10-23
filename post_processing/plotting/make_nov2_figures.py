#!/usr/bin/env python3
"""make_nov2_figures.py - rebuild every CIM comparison figure for an update folder.

Run from the folder that holds CIM/<surface>/<E>eV/erosion.csv, with this script
and the two plotting scripts in plotting/:

    python plotting/make_nov2_figures.py

It copies each erosion.csv into figures/csv/ under the names the plotting scripts
expect (erosion_WC_basal_C_100eV.csv, ...), adds surface_height_A to the W(001)
files (z_surf - z_surf at cycle 1), and writes
    figures/          full-size figures + erosion_summary_equal_fluence.csv
    figures/slide2/   6.0 x 3.9 in (two per slide)
    figures/slide3/   4.6 x 3.6 in (three per slide)
The WC figures share axes (LIMITS below) so panels sit side by side.
"""
import os
import shutil
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.getcwd()
CIM = os.path.join(ROOT, "CIM")
FIG = os.path.join(ROOT, "figures")
CSV = os.path.join(FIG, "csv")
CSV_W = os.path.join(FIG, "csv_w001")

SURF = {"basal-C": "WC_basal_C", "basal-W": "WC_basal_W", "prism": "WC_prism", "W001": "W001"}
ENERGIES = [50, 100, 200]
LIMITS = dict(xlim=("0", "41"), reg=("-8", "0.5"), sput=("0", "340"), w_height=("-0.3", "0.8"))
COMPARE_CYCLE = "1200"
SIZES = {"": None, "slide2": ("6.0", "3.9"), "slide3": ("4.6", "3.6")}


def run(*args):
    cmd = [sys.executable] + [str(a) for a in args]
    print(">", " ".join(cmd[1:]))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def main():
    os.makedirs(CSV, exist_ok=True)
    os.makedirs(CSV_W, exist_ok=True)
    for folder, name in SURF.items():
        for E in ENERGIES:
            src = os.path.join(CIM, folder, "%deV" % E, "erosion.csv")
            if not os.path.isfile(src):
                continue
            dst = os.path.join(CSV, "erosion_%s_%deV.csv" % (name, E))
            shutil.copyfile(src, dst)
            if name == "W001":
                d = pd.read_csv(src)
                d["surface_height_A"] = d["z_surf_A"] - d["z_surf_A"].iloc[0]
                d.to_csv(os.path.join(CSV_W, os.path.basename(dst)), index=False)
            print("  %s -> %s" % (os.path.relpath(src, ROOT), os.path.relpath(dst, ROOT)))

    en = os.path.join(HERE, "plot_erosion_energies.py")
    co = os.path.join(HERE, "plot_erosion_compare.py")
    for sub, size in SIZES.items():
        out = os.path.join(FIG, sub) if sub else FIG
        fs = ["--figsize", *size] if size else []
        run(en, "--dir", CSV, "--out", out, "--planes", "bC", "bW", "pr", "--x-axis", "ns",
            "--xlim", *LIMITS["xlim"], "--reg-ylim", *LIMITS["reg"],
            "--sput-ylim", *LIMITS["sput"], *fs)
        if sub != "slide2":
            run(en, "--dir", CSV_W, "--out", out, "--planes", "W001", "--x-axis", "ns",
                "--reg-col", "surface_height_A", "--reg-smooth", "25",
                "--reg-label", "Surface height change (Å)", "--reg-legend-loc", "upper left",
                "--reg-ylim", *LIMITS["w_height"],
                "--reg-title", "Surface height (swelling) - {plane}", *fs)
        for E in (100, 200):
            extra = ["--cluster-cols", "2"] if sub == "slide2" else []
            run(co, "--dir", CSV, "--out", out, "--energy", E, "--no-table", *fs, *extra)
    # order-deficit version of W(001), for the backup slide
    run(en, "--dir", CSV, "--out", FIG, "--planes", "W001", "--x-axis", "ns",
        "--mode", "separate", "--suffix", "_RD_noisy")
    extra = os.path.join(FIG, "erosion_sputtering_W001_RD_noisy.png")
    if os.path.isfile(extra):
        os.remove(extra)                      # identical to erosion_sputtering_W001.png
    run(co, "--dir", CSV, "--out", FIG, "--energy", 100, "--compare-cycle", COMPARE_CYCLE,
        "--cluster-cols", "2")
    print("done: figures in", FIG)


if __name__ == "__main__":
    main()
