#!/usr/bin/env python3
"""
analyze_cim_replicas.py - sputtering yield with error bars from CIM replicas.

Reads every <MAIN>/reps/bNNNN_rRR/ made by make_cim_replicas.sh, sums the
W and C atoms recorded in sputtering.dat over the cycles each replica ran
(cycle_time.dat), and reports per branch point.  Works for WC replicas
(sputtering.dat: cycle step time nW nC KE_W KE_C) and pure-W replicas
(cycle step time nW KE_W); the material is read from args.txt
(-var Eion = W deck) or forced with --mat.  Per branch:

    n replicas, mean yield, SD, SEM, 95% CI (Student t), dispersion index

The yield of replica r is Y_r = (atoms sputtered) / (cycles completed),
i.e. atoms per incident He.  Replicas are independent by construction
(different random streams from the same checkpoint), so the SEM is simply
SD/sqrt(n).  The dispersion index var(counts)/mean(counts) tells you
whether sputtering behaves like a Poisson process (~1) or comes in bursts
(>1, several atoms per event) - it sets how many replicas you need.

Optionally the main run is added as one more sample per window
(--main-csv erosion_WC_basal_C_100eV.csv, cumulative sputtered_C/W columns).

Outputs in --out:
    cim_replica_yields.csv   one row per replica
    cim_replica_summary.csv  one row per branch point
    cim_replica_yields.png   yield vs fluence at each branch point (C for WC, W for pure W)

Usage:
    python analyze_cim_replicas.py <MAIN>/reps --area 3567 --label "WC basal C-term, 100 eV"
    python analyze_cim_replicas.py <MAIN>/reps --area 3567 --main-csv erosion_WC_basal_C_100eV.csv \
        --sim 0.24 0.01            # SIM yield +- error, drawn as a band
"""
import argparse
import glob
import math
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T95 = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57, 6: 2.45, 7: 2.36, 8: 2.31,
       9: 2.26, 10: 2.23, 12: 2.18, 15: 2.13, 20: 2.09, 30: 2.04}


def t95(dof):
    if dof <= 0:
        return float("nan")
    keys = sorted(T95)
    for k in keys:
        if dof <= k:
            return T95[k]
    return 1.96


def read_numeric(path, ncol_min):
    rows = []
    if not os.path.isfile(path):
        return np.empty((0, ncol_min))
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split()
            try:
                rows.append([float(x) for x in parts[:ncol_min]])
            except ValueError:
                continue
    return np.array(rows) if rows else np.empty((0, ncol_min))


def material_of(d, forced=None):
    if forced:
        return forced
    a = os.path.join(d, "args.txt")
    if os.path.isfile(a) and "-var Eion" in open(a).read():
        return "W"
    return "WC"


def replica_row(d, forced_mat=None):
    m = re.search(r"b(\d+)_r(\d+)$", d.rstrip("/"))
    branch, rep = int(m.group(1)), int(m.group(2))
    mat = material_of(d, forced_mat)
    ct = read_numeric(os.path.join(d, "cycle_time.dat"), 2)       # cycle t_cycle ...
    # WC: cycle step time nW nC KE_W KE_C   |   W: cycle step time nW KE_W
    sp = read_numeric(os.path.join(d, "sputtering.dat"), 5 if mat == "WC" else 4)
    if ct.size == 0:
        return None
    first, last = branch + 1, int(ct[:, 0].max())
    done = int(np.sum((ct[:, 0] >= first) & (ct[:, 0] <= last)))
    sel = (sp[:, 0] >= first) & (sp[:, 0] <= last) if sp.size else np.array([], bool)
    nW = int(sp[sel, 3].sum()) if sp.size else 0
    nC = int(sp[sel, 4].sum()) if (sp.size and mat == "WC") else 0
    t_ok = bool(np.all(np.abs(ct[:, 1] - 20.0) < 0.01)) if ct.shape[1] > 1 else None
    return dict(material=mat, branch=branch, rep=rep, first_cycle=first, last_cycle=last, cycles=done,
                W=nW, C=nC, Y_W=nW / done if done else np.nan, Y_C=nC / done if done else np.nan,
                clock_20ps=t_ok, dir=d)


def main_rows(csv, branches, ncyc):
    d = pd.read_csv(csv).set_index("cycle")
    out = []
    for b in branches:
        lo, hi = b, b + ncyc                     # cumulative at b (0 -> 0) and at b + ncyc
        if hi > d.index.max():
            continue
        col = lambda name, c: d.loc[c, name] if name in d.columns else 0
        c0 = 0 if lo == 0 else col("sputtered_C", lo)
        w0 = 0 if lo == 0 else col("sputtered_W", lo)
        nC, nW = int(col("sputtered_C", hi) - c0), int(col("sputtered_W", hi) - w0)
        out.append(dict(branch=b, rep=0, first_cycle=lo + 1, last_cycle=hi, cycles=ncyc, W=nW, C=nC,
                        Y_W=nW / ncyc, Y_C=nC / ncyc, clock_20ps=None, dir="main run"))
    return out


def summarize(df, sp):
    """sp = the species whose yield gets the replica statistics (C for WC, W for pure W)."""
    rows = []
    for b, g in df.groupby("branch"):
        n = len(g)
        y = g["Y_" + sp].values
        sd = y.std(ddof=1) if n > 1 else float("nan")
        sem = sd / math.sqrt(n) if n > 1 else float("nan")
        disp = g[sp].var(ddof=1) / g[sp].mean() if n > 1 and g[sp].mean() > 0 else float("nan")
        totW, totN = int(g["W"].sum()), int(g["cycles"].sum())
        # W: usually 0-few events -> Poisson 95% upper limit on the pooled yield
        w_ul = (3.0 if totW == 0 else totW + 1.96 * math.sqrt(totW)) / totN
        rows.append({"branch": b, "n": n, "species": sp, "Y_mean": y.mean(), "Y_SD": sd, "Y_SEM": sem,
                     "Y_CI95": t95(n - 1) * sem if n > 1 else float("nan"),
                     "dispersion": disp, "W_total": totW, "impacts_total": totN,
                     "Y_W_pooled": totW / totN, "Y_W_95UL": w_ul,
                     "n_for_10pct_SEM": math.ceil((sd / y.mean() / 0.10) ** 2) if n > 1 and y.mean() > 0 else None})
    return pd.DataFrame(rows)


def plot(df, summ, area, label, sim, out, sp="C"):
    per = 1.0 / (area * 1e-20) / 1e19 if area else None
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    rng = np.random.default_rng(0)
    for b, g in df.groupby("branch"):
        ncyc = g["cycles"].median()
        x = (b + ncyc / 2) * per if per else b + ncyc / 2
        jit = (rng.random(len(g)) - 0.5) * (0.04 * (per or 1) * 100)
        main = g["rep"] == 0
        ax.scatter(x + jit[~main], g.loc[~main, "Y_" + sp], s=16, color="#e2711d", alpha=0.45, lw=0)
        if main.any():
            ax.scatter([x] * int(main.sum()), g.loc[main, "Y_" + sp], s=40, marker="D", facecolor="none",
                       edgecolor="k", label="main run" if b == df["branch"].min() else None)
    for _, r in summ.iterrows():
        ncyc = df.loc[df["branch"] == r["branch"], "cycles"].median()
        x = (r["branch"] + ncyc / 2) * per if per else r["branch"] + ncyc / 2
        ax.errorbar(x, r["Y_mean"], yerr=r["Y_CI95"], fmt="o", color="#7d0a12", ms=7,
                    capsize=5, lw=1.6, label="mean ± 95% CI" if r["branch"] == summ["branch"].min() else None)
    if sim:
        ax.axhspan(sim[0] - sim[1], sim[0] + sim[1], color="#2b7bba", alpha=0.18, label="SIM (pristine)")
        ax.axhline(sim[0], color="#2b7bba", lw=1.2)
    ax.set_xlabel(r"Fluence at window centre ($10^{19}$ He m$^{-2}$)" if per else "Cycle at window centre")
    ax.set_ylabel("%s sputtering yield (atoms / He)" % sp)
    ax.set_title("CIM replicas: %s yield vs fluence" % sp + (" - " + label if label else ""), fontsize=12)
    ax.grid(True, color="#b8b8b8", ls=":", lw=0.8, alpha=0.7)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9, loc="best")
    fig.tight_layout()
    p = os.path.join(out, "cim_replica_yields.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print("wrote", p)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reps", help="the reps/ directory made by make_cim_replicas.sh")
    ap.add_argument("--area", type=float, default=None, help="impact area in A^2 (fluence axis)")
    ap.add_argument("--label", default="")
    ap.add_argument("--main-csv", default=None, help="erosion CSV of the main run (adds it as rep 0)")
    ap.add_argument("--ncyc", type=int, default=100, help="cycles per replica window")
    ap.add_argument("--sim", nargs=2, type=float, default=None, metavar=("Y", "ERR"))
    ap.add_argument("--mat", choices=["WC", "W"], default=None, help="override the material read from args.txt")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()

    rows = [r for r in (replica_row(d, a.mat) for d in sorted(glob.glob(os.path.join(a.reps, "b*_r*")))) if r]
    if not rows:
        raise SystemExit("no finished replicas found under %s" % a.reps)
    df = pd.DataFrame(rows)
    short = df[df["cycles"] < a.ncyc]
    if len(short):
        print("WARNING: %d replica(s) have fewer than %d cycles (still running?)" % (len(short), a.ncyc))
    if (df["clock_20ps"] == False).any():                                      # noqa: E712
        print("WARNING: some replicas have cycles that are not 20.00 ps - check cycle_time.dat")
    if a.main_csv:
        df = pd.concat([df, pd.DataFrame(main_rows(a.main_csv, sorted(df["branch"].unique()), a.ncyc))],
                       ignore_index=True)
    os.makedirs(a.out, exist_ok=True)
    df.sort_values(["branch", "rep"]).to_csv(os.path.join(a.out, "cim_replica_yields.csv"), index=False)
    mats = set(df["material"].dropna())
    if len(mats) > 1:
        raise SystemExit("replicas of both W and WC under %s - analyse one campaign at a time" % a.reps)
    sp = "W" if mats == {"W"} else "C"
    summ = summarize(df, sp)
    summ.to_csv(os.path.join(a.out, "cim_replica_summary.csv"), index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30, "display.float_format", "{:.4g}".format):
        print(summ.to_string(index=False))
    plot(df, summ, a.area, a.label, a.sim, a.out, sp)


if __name__ == "__main__":
    main()
