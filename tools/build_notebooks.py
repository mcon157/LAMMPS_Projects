#!/usr/bin/env python3
"""Build the .ipynb notebooks from the .py scripts and run them, so the figures,
renders and GIFs are stored inside each notebook.

The .py files stay the source (they run headless on a cluster); every notebook
is regenerated from its script, so edit the .py and rebuild:

    python tools/build_notebooks.py                     # all notebooks
    python tools/build_notebooks.py rand_schematics     # just one
    python tools/build_notebooks.py --data paths.json   # point data-driven notebooks at your data
    python tools/build_notebooks.py --no-exec           # code only, no outputs

paths.json maps a notebook name to the settings used while it runs, e.g.
    {"ovito_insertion_movies": {"NB_RUN": "/path/to/basal_C/200eV"},
     "plot_erosion_energies":  {"NB_CSV": "/path/to/figures/csv"}}
The saved notebook keeps the generic settings; only its outputs come from that run.

Needs nbformat, nbclient, ipykernel, Pillow, plus what each script needs.
"""
import argparse
import ast
import copy
import json
import os
import re
import shutil
import sys
import tempfile

import nbformat
from nbclient import NotebookClient

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---------------------------------------------------------------------------
# One entry per notebook.
#   params   the settings cell (user-facing); ARGS becomes sys.argv.  Notebook
#            variables start with NB_ so the script's own globals cannot overwrite them.
#   execute  "all"  run every cell here (self-contained scripts)
#            "show" code cells are stored but not run (they need cluster data or
#                   OVITO); only the display cells run, on example outputs
#   show     (folder variable, [file names]) displayed at the end
#   extra    code run after the script (demo / plotting cells)
#   note     markdown shown above the outputs
# ---------------------------------------------------------------------------
NOTEBOOKS = {
    "rand_schematics": dict(
        script="SIM_random/rand_schematics.py", execute="all",
        params='NB_OUT = "nb_figs"\nARGS = ["--n", "50", "--out", NB_OUT]',
        show=("NB_OUT", ["rand_steps_WC.png", "rand_crosssection_WC.png", "rand_orientations.png",
                          "rand_example_surfaces.png"])),
    "rand_surface_views": dict(
        script="SIM_random/rand_surface_views.py", execute="all",
        params='NB_OUT = "nb_figs"\nARGS = ["--s", "7", "--out", NB_OUT]',
        show=("NB_OUT", ["rand_regular_vs_random.png"])),
    "rand_orientation_gif": dict(
        script="SIM_random/rand_orientation_gif.py", execute="all",
        params='NB_OUT = "nb_figs"\nARGS = ["--mat", "W", "--first", "1", "--last", "30", "--out", NB_OUT]',
        extra=['# the same for WC\nsys.argv = [SCRIPT, "--mat", "WC", "--first", "1", "--last", "30", "--out", NB_OUT]\nmain()'],
        show=("NB_OUT", ["*.gif"])),
    "analyze_rand": dict(
        script="SIM_random/analyze_rand.py", execute="show",
        params='NB_RUN = "."          # folder with rand_W_*eV.dat and rand_W_surfaces.dat\n'
               'ARGS = ["--mat", "W", "--dir", NB_RUN, "--out", NB_RUN]',
        note="Example output: yield of each random W surface on the inverse pole figure, 300 and 500 eV.",
        show=("NB_RUN", ["rand_W_300eV_ipf.png", "rand_W_500eV_ipf.png"])),
    "cim_replica_schematic": dict(
        script="CIM_statistics/cim_replica_schematic.py", execute="all",
        params='NB_OUT = "nb_figs"\nos.makedirs("nb_figs", exist_ok=True)\n'
               'ARGS = ["--nrep", "10", "--ncyc", "100", "--out", os.path.join(NB_OUT, "cim_replica_branches.png")]',
        show=("NB_OUT", ["cim_replica_branches.png"])),
    "summarize_replicas": dict(
        script="CIM_statistics/summarize_replicas.py", execute="show",
        params='NB_RUNS = "runs"          # the CIM_replicas/runs folder\nARGS = [NB_RUNS]',
        note="Run this on a CIM_replicas folder: it writes replica_yields.csv, replica_summary.csv and "
             "replica_yields.png into each run folder (no example data is included).",
        show=None),
    "elstop": dict(
        script="WC_electronic_stopping/elstop.py", execute="all",
        params='ARGS = ["basis", "--backend", "rustbca", "-o", "basis_nb.dat"]',
        extra=['''# stopping cross sections of the basis, and a 50/50 WC table built from it
import numpy as np, matplotlib.pyplot as plt
rows = [l.split() for l in open("basis_nb.dat") if l.strip() and not l.startswith("#")]
head = [l for l in open("basis_nb.dat") if l.startswith("# E_eV")][0].split()[1:]
data = np.array(rows, float)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for j, name in enumerate(head[1:], 1):
    proj, targ = name.split("@")
    ax[0 if targ == "W" else 1].loglog(data[:, 0], data[:, j], label=proj + " in " + targ)
for a, t in zip(ax, ("W", "C")):
    a.set_xlabel("E (eV)"); a.set_ylabel("sigma (eV A$^2$)"); a.set_title("target " + t); a.legend()
plt.tight_layout(); plt.show()''',
               '''sys.argv = [SCRIPT, "table", "--basis", "basis_nb.dat", "--types", "W", "C", "He",
            "--target", "W:0.5,C:0.5", "--rho", "15.63", "-o", "elstop_WC_nb.table"]
main()
print(open("elstop_WC_nb.table").read()[:600])''']),
    "ovito_insertion_movies": dict(
        script="post_processing/ovito_insertion_movies.py", execute="show",
        params='NB_RUN = "."          # run folder with dump.bulk-*/dump.full-* and sputter_summary.dat\n'
               'ARGS = [NB_RUN]       # add "--fast" for metrics and plots only, "--surface", "prism", ...',
        note="Example output: WC basal (0001), C-terminated, 200 eV He, up to cycle 1226 (final-frame "
             "stills; the full GIFs are too large to embed).",
        show=("NB_RUN", ["erosion_combined.png", "composition_profile.png", "he_profile.png",
                          "full_iso_final.png", "full_side_final.png", "planes_side_final.png", "depth_slice.png",
                          "full_iso.gif"])),
    "he_cluster_movies": dict(
        script="post_processing/he_cluster_movies.py", execute="show",
        params='NB_RUN = "."          # run folder with dump.full-*.lammps\n'
               'ARGS = ["--dir", NB_RUN, "--jobs", "4"]',
        note="Example output: WC basal (0001), C-terminated, 200 eV He.",
        show=("os.path.join(NB_RUN, 'he_clusters')", ["he_cluster_evolution.png", "he_depth_profile.png",
                                                       "he_clusters_final.png", "he_clusters.gif"])),
    "plot_erosion_energies": dict(
        script="post_processing/plotting/plot_erosion_energies.py", execute="all",
        params='NB_CSV = "."          # erosion_WC_basal_C_50eV.csv, ... (see PLANE_FILES)\n'
               'NB_OUT = "nb_figs"\nos.makedirs(NB_OUT, exist_ok=True)\n'
               'ARGS = ["--dir", NB_CSV, "--out", NB_OUT]',
        show=("NB_OUT", ["erosion_regression_bC.png", "erosion_sputtering_bC.png", "erosion_regression_bW.png",
                          "erosion_sputtering_bW.png", "erosion_regression_pr.png", "erosion_sputtering_pr.png",
                          "erosion_regression_W001.png", "erosion_sputtering_W001.png"])),
    "plot_erosion_compare": dict(
        script="post_processing/plotting/plot_erosion_compare.py", execute="all",
        params='NB_CSV = "."          # the same CSVs as plot_erosion_energies.py\n'
               'NB_OUT = "nb_figs"\nos.makedirs(NB_OUT, exist_ok=True)\n'
               'ARGS = ["--dir", NB_CSV, "--out", NB_OUT, "--energy", "100"]',
        show=("NB_OUT", ["surface_W_pct_100eV.png", "he_retention_100eV.png", "he_clusters_100eV.png",
                          "W001_surface_He.png"])),
    "schematics": dict(
        script="He_BE/scripts/schematics.py", execute="all",
        params='ARGS = []              # writes ../schematics/S1..S6 (.png + .pdf)',
        show=("os.path.join('..', 'schematics')", ["S1_model_schema.png", "S2_wc_structure.png",
                                                    "S3_cluster_family.png", "S4_binding_bookkeeping.png",
                                                    "S5_search_protocol.png", "S6_diagnostics.png"])),
    "gifs": dict(
        script="He_BE/scripts/gifs.py", execute="all",
        params='ARGS = []              # all four; or e.g. ["G2"]',
        note="Illustrative animations: the energies are made up.", show_kw="gif_kb=900",
        show=("os.path.join('..', 'gifs')", ["G1_workflow.gif", "G2_search.gif", "G3_vacancy_filling.gif",
                                              "G4_jw_vacancy_rule.gif"])),
    "renders": dict(
        script="He_BE/scripts/renders.py", execute="show",
        params='NB_RENDERS = os.path.join("..", "renders")   # where renders.py writes\n'
               'ARGS = []              # needs OVITO, and ../../he_be_v3, ../../wc_ef for R2-R4',
        note="Renders from the binding-energy runs (made with OVITO's Tachyon renderer).",
        show=("NB_RENDERS", ["R1_wc_supercell.png", "R2_v3_clusters.png", "R3_reference_phases.png",
                             "R4_point_defects.png"])),
}

SHOW_HELPER = r'''# display helper: embeds the outputs in this notebook (big GIFs and PNGs are shrunk first)
import glob, io
from IPython.display import Image, display
from PIL import Image as PILImage

def _shrink_png(path, width=1200):
    im = PILImage.open(path)
    if im.width > width:
        im = im.resize((width, round(im.height * width / im.width)), PILImage.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=True)
    return buf.getvalue()

def _shrink_gif(path, max_kb):
    src = PILImage.open(path)
    frames, durations = [], []
    try:
        while True:
            frames.append(src.convert("RGB"))
            durations.append(src.info.get("duration", 80))
            src.seek(src.tell() + 1)
    except EOFError:
        pass
    w0 = min(frames[0].width, 640)
    # smaller frames first; drop frames only if that is not enough
    for step, width in [(1, w0), (1, int(w0 * 0.8)), (1, int(w0 * 0.65)), (2, int(w0 * 0.8)),
                        (2, int(w0 * 0.65)), (4, int(w0 * 0.65)), (8, int(w0 * 0.5))]:
        sel = frames[::step]
        size = (width, round(frames[0].height * width / frames[0].width))
        sel = [f.resize(size, PILImage.LANCZOS).quantize(colors=128) for f in sel]
        buf = io.BytesIO()
        sel[0].save(buf, "GIF", save_all=True, append_images=sel[1:], loop=0,
                    duration=int(sum(durations[::step]) / len(sel)) * step, optimize=True)
        if buf.tell() <= max_kb * 1024:
            break
    return buf.getvalue()

def show(folder, names, png_kb=700, gif_kb=1500):
    for pat in names:
        for p in sorted(glob.glob(os.path.join(folder, pat))):
            data = open(p, "rb").read()
            if p.lower().endswith(".gif"):
                if len(data) > gif_kb * 1024:
                    data = _shrink_gif(p, gif_kb)
            elif len(data) > png_kb * 1024:
                data = _shrink_png(p)
            print(os.path.basename(p))
            display(Image(data=data, format="png"))   # GIF bytes still animate in the browser
'''


# ---------------------------------------------------------------------------
def split_script(src):
    """Script -> (docstring, [code cell texts])."""
    tree = ast.parse(src)
    lines = src.splitlines()
    doc = ast.get_docstring(tree, clean=True) or ""
    start = 0
    if tree.body and isinstance(tree.body[0], ast.Expr) and isinstance(getattr(tree.body[0], "value", None), ast.Constant) \
            and isinstance(tree.body[0].value.value, str):
        start = tree.body[0].end_lineno
    body = lines[start:]
    off = start
    cuts = {0}
    nodes = [n for n in tree.body if n.lineno > start]
    # the import block ends at the first assignment / def / class
    first_def = next((n for n in nodes if isinstance(n, (ast.Assign, ast.AnnAssign, ast.AugAssign,
                                                          ast.FunctionDef, ast.ClassDef))), None)
    if first_def is not None:
        cuts.add(_lead(lines, first_def.lineno - 1) - off)
    for i, l in enumerate(body):
        banner = re.match(r"^# ?[=\-]{8,}", l) or re.match(r"^# --- \w", l)
        if banner and (i == 0 or not re.match(r"^#", body[i - 1])):
            cuts.add(i)
    for n in nodes:
        if isinstance(n, ast.If) and "__name__" in ast.dump(n.test):
            cuts.add(_lead(lines, n.lineno - 1) - off)
    cuts = sorted(c for c in cuts if 0 <= c < len(body))
    cells = []
    for a, b in zip(cuts, cuts[1:] + [len(body)]):
        chunk = body[a:b]
        cells += _split_long(chunk, a + off, nodes, lines)
    out = []
    for c in cells:
        t = "\n".join(c).strip("\n")
        if not t.strip():
            continue
        comment_only = re.match(r"^(#[^\n]*\n?)+$", t + "\n")
        if out and (comment_only or (t.count("\n") < 6 and "__name__" not in t and
                                     out[-1].count("\n") < 120)):
            out[-1] += "\n\n" + t                  # no tiny cells
        else:
            out.append(t)
    return doc, out


def _lead(lines, i):
    """Index of the first comment/decorator line directly above line i."""
    while i > 0 and (lines[i - 1].startswith("#") or lines[i - 1].startswith("@")):
        i -= 1
    return i


def _split_long(chunk, base, nodes, lines, limit=160):
    if len(chunk) <= limit:
        return [chunk]
    starts = sorted({_lead(lines, n.lineno - 1) - base for n in nodes
                     if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and base <= n.lineno - 1 < base + len(chunk)})
    pieces, cur = [], 0
    for s in starts:
        if s - cur >= 40:
            pieces.append(chunk[cur:s])
            cur = s
    pieces.append(chunk[cur:])
    return [p for p in pieces if p]


def build(name, spec, overrides, run):
    path = os.path.join(ROOT, spec["script"])
    src = open(path).read()
    doc, cells = split_script(src)
    script = os.path.basename(path)
    nb = nbformat.v4.new_notebook()
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    first, _, rest = doc.partition("\n")
    md = "# %s\n\n*Generated from `%s` by `tools/build_notebooks.py`. Edit the script, then rebuild.*" % (script, script)
    if doc:
        md += "\n\n```text\n%s\n```" % doc.strip()
    nb.cells.append(nbformat.v4.new_markdown_cell(md))
    params = "# settings (the same as the script's command line)\nimport os, sys\nSCRIPT = %r\n%s\n" \
             "sys.argv = [SCRIPT] + ARGS\n__file__ = os.path.abspath(SCRIPT)" % (script, spec["params"])
    nb.cells.append(nbformat.v4.new_code_cell(params))
    skip = spec["execute"] == "show"
    if skip:
        nb.cells.append(nbformat.v4.new_markdown_cell(
            "The cells below are the script.  They need the run's dumps (and OVITO where it is imported), "
            "so they were not run when this notebook was built; the outputs further down come from a "
            "cluster run."))
    for c in cells:
        # in a notebook sys.exit() raises; call main() directly
        c = re.sub(r"^(\s*)sys\.exit\((main\([^)]*\))\)\s*$", r"\1\2", c, flags=re.M)
        cell = nbformat.v4.new_code_cell(c)
        if skip:
            cell.metadata["tags"] = ["skip-execution"]
        nb.cells.append(cell)
    for e in spec.get("extra", []):
        cell = nbformat.v4.new_code_cell(e)
        if skip:
            cell.metadata["tags"] = ["skip-execution"]
        nb.cells.append(cell)
    if spec.get("show"):
        folder, names = spec["show"]
        nb.cells.append(nbformat.v4.new_markdown_cell("## Outputs" + ("\n\n" + spec["note"] if spec.get("note") else "")))
        nb.cells.append(nbformat.v4.new_code_cell(SHOW_HELPER))
        kw = (", " + spec["show_kw"]) if spec.get("show_kw") else ""
        nb.cells.append(nbformat.v4.new_code_cell("show(%s, %r%s)" % (folder, names, kw)))
    elif spec.get("note"):
        nb.cells.append(nbformat.v4.new_markdown_cell(spec["note"]))
    out = os.path.splitext(path)[0] + ".ipynb"
    if run:
        execute(nb, path, overrides.get(name, {}))
    nbformat.write(nb, out)
    size = os.path.getsize(out)
    print("%-28s -> %s (%.1f MB)" % (name, os.path.relpath(out, ROOT), size / 1e6))


def execute(nb, path, over):
    """Run a copy of the notebook in a scratch copy of the repo; keep only the outputs."""
    tmp = tempfile.mkdtemp(prefix="nbbuild_")
    try:
        work = os.path.join(tmp, "repo")
        shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns(".git", "*.ipynb", "__pycache__"))
        cwd = os.path.join(work, os.path.dirname(os.path.relpath(path, ROOT)))
        run_nb = copy.deepcopy(nb)
        p = run_nb.cells[1]
        for k, v in over.items():
            p.source = re.sub(r"^%s = .*$" % re.escape(k), "%s = %r" % (k, v), p.source, flags=re.M)
        client = NotebookClient(run_nb, timeout=1800, kernel_name="python3",
                                resources={"metadata": {"path": cwd}})
        client.execute()
        hide = [work, cwd, tmp] + [str(v) for v in over.values()]
        for src, dst in zip(run_nb.cells, nb.cells):
            if dst.cell_type == "code" and dst is not nb.cells[1]:
                dst.outputs = [_clean(o, hide) for o in src.get("outputs", [])]
                dst.execution_count = src.get("execution_count")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _clean(o, hide):
    """Remove build-machine paths from printed text."""
    if o.get("output_type") == "stream":
        t = o["text"]
        for h in sorted(hide, key=len, reverse=True):
            if h:
                t = t.replace(h + os.sep, "").replace(h, ".")
        o["text"] = t
    return o


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("names", nargs="*", help="notebook names (default: all)")
    ap.add_argument("--data", help="JSON with per-notebook settings for the run")
    ap.add_argument("--no-exec", action="store_true")
    a = ap.parse_args()
    over = json.load(open(a.data)) if a.data else {}
    for name, spec in NOTEBOOKS.items():
        if a.names and name not in a.names:
            continue
        build(name, spec, over, not a.no_exec)


if __name__ == "__main__":
    main()
