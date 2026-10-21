# He_BE: He binding in WC (reference pack)

Documentation for the He solution, clustering and vacancy-trapping calculations in WC, with bcc W as the benchmark. Every reported energy is a 0 K conjugate-gradient minimum. MD is used only to search between minima, never to measure an energy.

The calculation decks and their analysis are **not** in this folder:
- `he_be_v3/` (`make_v3.py`, `analyze_v3.py`);
- `wc_ef/` (`in.phases`, `make_ef.py`).

The pack describes them. `renders.py` R2–R4 need them next to this folder.

| Path | Contents |
|---|---|
| `model_schema.yaml` | machine-readable schema: every input (value, source, sensitivity), the process, expected outputs with ranges, limitations |
| `sheet/quick_reference.pdf` (`.tex`) | 11-page landscape look-up sheet. See the page list below the table |
| `scripts/schematics.py` / `.ipynb` | S1–S6 schematics → `schematics/` (PNG 200 dpi + PDF) |
| `scripts/gifs.py` / `.ipynb` | G1–G4 illustrative animations → `gifs/` (about 1–2 min) |
| `scripts/renders.py` / `.ipynb` | R1–R4 ray-traced renders of the cells (OVITO) → `renders/` |
| `scripts/render_dump.py` | renders one of your relaxed minima in the same style |
| `scripts/style.py` | shared colours and drawing helpers |

The sheet covers:
1. the schema tables, expected outputs and limitations;
2. the workflow and files;
3. every LAMMPS command used and its knobs;
4. every equation in the analysis and diagnostics;
5. a glossary.

Compile it from `sheet/` with `latexmk -pdf quick_reference.tex`.

The notebooks show every figure, GIF and render, so the output folders are not stored here.

## Figures

| | Shows |
|---|---|
| S1 model schema | inputs → LAMMPS → analysis → outputs, with a limitations panel |
| S2 WC structure | top and side views: W, C, channel site, V<sub>W</sub>, V<sub>C</sub> (real a, c) |
| S3 cluster family | the eight v3 clusters, each a parent plus one site |
| S4 binding bookkeeping | which cells enter each binding energy |
| S5 search protocol | temperature ladder (CG → 300 / 614.6 / 900 K anneals → quench → CG), choice of the minimum |
| S6 diagnostics | intact / escaped / trap-mutated, and how each is detected |
| G1 workflow | the pipeline stage by stage |
| G2 search | 1D energy-landscape view of the anneal-ladder search |
| G3 vacancy filling | He added to a W vacancy one at a time until the last one escapes |
| G4 JW vacancy rule | Juslin & Wirth's "remove the highest-PE neighbour" rule for growing voids in W |
| R1 WC supercell | the 8×5×8 cell (1280 atoms): 3D, down c, down a |
| R2 v3 clusters | the eight vacancy clusters with 4 He |
| R3 reference phases | the six reference crystals of `in.phases` |
| R4 point defects | the formation-energy defect cells |

Colours: W gray, C blue, He yellow, vacancy red. Only S2 and the renders use real geometry. The energies in S5 and in the GIFs are made up and labelled as such.

## Run

```bash
cd scripts
python3 schematics.py && python3 gifs.py           # numpy, matplotlib, Pillow
pip install ovito && python3 renders.py            # R2-R4 need ../../he_be_v3 and ../../wc_ef (else it stops after R1)
python3 render_dump.py dumps/vw1_14_c3.lammps --cfg cfg/vw1_14.lmp       # one of your minima
python3 render_dump.py dumps/vw1_14_c3.lammps --full                     # whole cell
```

The renders need a Vulkan driver:
- **Laptop with a GPU.** Works as is.
- **Headless Linux.** Install `mesa-vulkan-drivers`, then `export VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/lvp_icd.json`.
