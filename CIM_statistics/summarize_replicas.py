#!/usr/bin/env python3
"""summarize_replicas.py - replica statistics for every run folder.

    python3 summarize_replicas.py [runs/WC_basal_C_50eV ...] [--include-incomplete]

Reads the seed folders seed_b<branch>_s<seed>/ of each run folder, sums the
atoms in each seed's sputtering.dat over the cycles it ran, and writes into
the run folder:
    replica_yields.csv    per seed: cycles, W and C sputtered, yields, complete,
                          cycle-clock check
    replica_summary.csv   per save point: seeds, mean yield, SD, SEM, 95 % CI
                          (Student t), dispersion (var/mean; 1 = Poisson),
                          seeds for a 10 % SEM, pooled W yield + 95 % upper limit
    replica_yields.png    per-seed yields and mean +- 95 % CI vs fluence
plus runs/all_replica_summary.csv.  Statistics use the C yield for WC and the W
yield for pure W; only complete seeds count unless --include-incomplete.
Needs numpy; matplotlib optional.
"""
import argparse
import csv
import glob
import math
import os
import re

import numpy as np

# impact area (A^2) of the central 2/3 x 2/3 of each slab -> fluence axis
AREA = {"basal": 3567.0, "prism": 3420.0, "W": 3490.0}
T95 = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57, 6: 2.45, 7: 2.36, 8: 2.31, 9: 2.26,
       10: 2.23, 12: 2.18, 15: 2.13, 20: 2.09, 25: 2.06, 30: 2.04, 40: 2.02, 60: 2.00}


def t95(dof):
    for k in sorted(T95):
        if dof <= k:
            return T95[k]
    return 1.96


def numeric_rows(path, ncol):
    out = []
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if not p or p[0].startswith("#"):
                continue
            try:
                out.append([float(x) for x in p[:ncol]])
            except ValueError:
                pass
    return out


def run_info(run):
    name = os.path.basename(os.path.normpath(run))
    if os.path.isfile(os.path.join(run, "wche0.restart")):
        mat = "WC"
    elif os.path.isfile(os.path.join(run, "w0.restart")):
        mat = "W"
    else:
        a = glob.glob(os.path.join(run, "seed_b*_s*", "args.txt"))
        mat = "W" if a and "-var Eion" in open(a[0]).read() else "WC"
    m = re.search(r"_(\d+)eV$", name)
    energy = int(m.group(1)) if m else None
    area = AREA["W"] if mat == "W" else next((v for k, v in AREA.items() if k in name.lower()), None)
    return name, mat, energy, area


def seed_row(d, mat):
    m = re.search(r"seed_b(\d+)_s(\d+)$", os.path.normpath(d))
    branch, seed = int(m.group(1)), int(m.group(2))
    args = open(os.path.join(d, "args.txt")).read().split() if os.path.isfile(os.path.join(d, "args.txt")) else []
    av = {args[i + 1]: args[i + 2] for i in range(len(args) - 2) if args[i] == "-var"}
    first = int(av.get("istart", branch + 1))
    last_planned = int(av.get("nimplants", av.get("nend", first + 99)))
    rep = int(av.get("rep", branch + seed))
    ct = numeric_rows(os.path.join(d, "cycle_time.dat"), 2)
    cyc = sorted({int(r[0]) for r in ct if first <= r[0] <= last_planned})
    done = len(cyc)
    last = cyc[-1] if cyc else first - 1
    t_cycle = float(av.get("t_cycle", 20.0))
    clock_ok = all(abs(r[1] - t_cycle) < 0.01 for r in ct if first <= r[0] <= last_planned) if ct else None
    sp = numeric_rows(os.path.join(d, "sputtering.dat"), 5 if mat == "WC" else 4)
    sel = [r for r in sp if first <= r[0] <= last]
    nW = int(sum(r[3] for r in sel))
    nC = int(sum(r[4] for r in sel)) if mat == "WC" else 0
    return dict(branch=branch, seed=seed, rep=rep, first_cycle=first, last_cycle_planned=last_planned,
                cycles_done=done, complete=int(done >= last_planned - first + 1),
                W=nW, C=nC, Y_W=nW / done if done else float("nan"), Y_C=nC / done if done else float("nan"),
                clock_ok=("" if clock_ok is None else int(clock_ok)), folder=os.path.basename(os.path.normpath(d)))


def summarize(rows, mat, area, use_incomplete):
    sp = "W" if mat == "W" else "C"
    out = []
    for b in sorted({r["branch"] for r in rows}):
        g = [r for r in rows if r["branch"] == b and (use_incomplete or r["complete"]) and r["cycles_done"] > 0]
        n = len(g)
        if n == 0:
            out.append(dict(branch=b, n_seeds=0))
            continue
        y = np.array([r["Y_" + sp] for r in g])
        cnt = np.array([r[sp] for r in g], float)
        mean = y.mean()
        sd = y.std(ddof=1) if n > 1 else float("nan")
        sem = sd / math.sqrt(n) if n > 1 else float("nan")
        disp = cnt.var(ddof=1) / cnt.mean() if n > 1 and cnt.mean() > 0 else float("nan")
        totW = int(sum(r["W"] for r in g))
        totN = int(sum(r["cycles_done"] for r in g))
        ul = (3.0 if totW == 0 else totW + 1.96 * math.sqrt(totW)) / totN
        ncyc = np.median([r["cycles_done"] for r in g])
        centre = b + ncyc / 2.0
        out.append(dict(branch=b, n_seeds=n, n_incomplete=sum(1 for r in rows if r["branch"] == b and not r["complete"]),
                        species=sp, cycles_per_seed=int(ncyc), window_centre_cycle=centre,
                        fluence_centre_m2=(centre / (area * 1e-20)) if area else "",
                        Y_mean=mean, Y_SD=sd, Y_SEM=sem, Y_CI95=t95(n - 1) * sem if n > 1 else float("nan"),
                        dispersion=disp,
                        seeds_for_10pct_SEM=math.ceil((sd / mean / 0.10) ** 2) if n > 1 and mean > 0 else "",
                        W_total=totW, impacts_total=totN, Y_W_pooled=totW / totN, Y_W_95UL=ul))
    return out


def write_csv(path, rows, keys=None):
    if not rows:
        return
    keys = keys or list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def plot(run, name, rows, summ, mat, area):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    sp = "W" if mat == "W" else "C"
    xf = (lambda c: c / (area * 1e-20) / 1e19) if area else (lambda c: c)
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    rng = np.random.default_rng(0)
    for s in summ:
        if not s.get("n_seeds"):
            continue
        g = [r for r in rows if r["branch"] == s["branch"] and r["cycles_done"] > 0]
        x = xf(s["window_centre_cycle"])
        span = xf(s["window_centre_cycle"] + 20) - x
        ax.scatter(x + rng.uniform(-span, span, len(g)), [r["Y_" + sp] for r in g], s=18, color="#eb6834",
                   alpha=0.55, lw=0, label="one seed" if s is summ[0] else None)
        ax.errorbar(x, s["Y_mean"], yerr=s["Y_CI95"] if s["n_seeds"] > 1 else None, fmt="o", color="#222222",
                    ms=7, capsize=5, label="mean ± 95% CI" if s is summ[0] else None)
    ax.set_xlabel("fluence at window centre (10$^{19}$ He m$^{-2}$)" if area else "cycle at window centre")
    ax.set_ylabel("%s sputtering yield (atoms / He)" % sp)
    ax.set_title("%s: %s yield from seed replicas" % (name, sp), fontsize=12)
    ax.set_ylim(bottom=0)
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(run, "replica_yields.png"), dpi=200)
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*", default=["runs"],
                    help="run folders, or a folder that contains run folders (default: runs)")
    ap.add_argument("--include-incomplete", action="store_true")
    a = ap.parse_args()

    folders = []
    for p in a.runs:
        if glob.glob(os.path.join(p, "seed_b*_s*")):
            folders.append(p)
        else:
            folders += sorted(d for d in glob.glob(os.path.join(p, "*")) if glob.glob(os.path.join(d, "seed_b*_s*")))
    if not folders:
        raise SystemExit("no run folders with seed_b*_s* folders found in " + " ".join(a.runs))

    allrows = []
    for run in folders:
        name, mat, energy, area = run_info(run)
        rows = [seed_row(d, mat) for d in sorted(glob.glob(os.path.join(run, "seed_b*_s*"))) if os.path.isdir(d)]
        summ = summarize(rows, mat, area, a.include_incomplete)
        write_csv(os.path.join(run, "replica_yields.csv"), rows)
        write_csv(os.path.join(run, "replica_summary.csv"), summ)
        png = plot(run, name, rows, summ, mat, area)
        bad_clock = sum(1 for r in rows if r["clock_ok"] == 0)
        print("%s (%s, %s eV): %d seed folders, %d complete%s%s" % (
            name, mat, energy, len(rows), sum(r["complete"] for r in rows),
            "" if not bad_clock else ", WARNING: %d with a cycle length != t_cycle (20 ps)" % bad_clock,
            "" if png else " (no plot: matplotlib missing)"))
        for s in summ:
            if s.get("n_seeds"):
                print("   branch %5d: n=%2d  Y_%s = %.4g +- %.2g (95%% CI)  dispersion %.2f" % (
                    s["branch"], s["n_seeds"], s["species"], s["Y_mean"], s["Y_CI95"], s["dispersion"]))
            else:
                print("   branch %5d: no complete seeds yet" % s["branch"])
        for s in summ:
            allrows.append(dict(run=name, material=mat, energy_eV=energy, **s))
    top = os.path.dirname(os.path.normpath(folders[0])) or "."
    write_csv(os.path.join(top, "all_replica_summary.csv"), allrows)
    print("wrote replica_yields.csv / replica_summary.csv in each run folder and %s"
          % os.path.join(top, "all_replica_summary.csv"))


if __name__ == "__main__":
    main()
