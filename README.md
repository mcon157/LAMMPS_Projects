# LAMMPS workflows: He on W and WC

LAMMPS input decks, Slurm job scripts and Python analysis for:
- **He sputtering and surface erosion** of tungsten (W) and tungsten carbide (WC), with:
  - a continuous-impact model (CIM, one He per cycle into an evolving slab);
  - a random-orientation single-impact model (SIM);
  - seed replicas for error bars;
- **electronic stopping** in the CIM runs;
- **surface binding energies** of W and WC;
- **He binding in WC**: a reference pack of schematics, renders and a quick-reference sheet.

Each deck is a single LAMMPS input that runs as-is inside its own run folder, with the potential files it needs copied next to it.

## Layout

| Folder | What it does | Main files |
|---|---|---|
| [`potentials/`](potentials) | interatomic potentials used by every deck | `WCH.tersoff`, `w_eam4.fs`, `WHe.table`, `CHe.table`, `HeHe.table` |
| [`WC_therm/`](WC_therm) | WC thermal expansion, giving the lattice constants at T | `in.step1` |
| [`WC_CIM/`](WC_CIM) | WC continuous-impact runs: build a slab (step 1), then He cycles (step 2) | `in_*_092926.step1`, `in_CIM_092926.step2`, `run.slurm` |
| [`WC_electronic_stopping/`](WC_electronic_stopping) | the CIM with composition-aware electronic stopping | `in_CIM_ESdyn.step2`, `elstop.py`, `build_lammps.sh`, `run.slurm` |
| [`CIM_statistics/`](CIM_statistics) | seed replicas branched from CIM save points (error bars) | `setup_replicas.sh`, `submit_all.sh`, `summarize_replicas.py` |
| [`SIM_random/`](SIM_random) | single impacts on randomly oriented hemispheres (polycrystal average) | `in_rand_*.step1/2`, `run_rand.slurm`, `analyze_rand.py` |
| [`W_SBE/`](W_SBE), [`WC_SBE/`](WC_SBE) | surface binding energy and kickout threshold | `in.step*` |
| [`post_processing/`](post_processing) | OVITO metrics and renders, He clusters, cross-run plots | `ovito_insertion_movies.py`, `he_cluster_movies.py`, `plotting/` |
| [`He_BE/`](He_BE) | He-in-WC binding-energy reference pack | `model_schema.yaml`, `scripts/`, `sheet/quick_reference.pdf` |
| [`tools/`](tools) | rebuilds the notebooks from the scripts | `build_notebooks.py` |

Every folder has its own README with the run commands, inputs, outputs and options.

## How the pieces fit

```
potentials/ ─┬─> WC_therm  (a, c at T) ──> lattice constants in every WC deck
             │
             ├─> WC_CIM step 1 (slab: wche0.restart + cycle-0 dumps)
             │      └─> step 2, chained 3 h jobs (run.slurm): one He per 20 ps cycle
             │            ├─> dumps, sputtering.dat, sputter_summary.dat, he_cycle_<N>.restart
             │            ├─> post_processing/ovito_insertion_movies.py  -> erosion.csv, figures, renders
             │            ├─> post_processing/he_cluster_movies.py       -> He cluster statistics, GIF
             │            └─> CIM_statistics: seeds from he_cycle_<N>.restart -> yield +- 95 % CI
             │      post_processing/plotting/: all energies / surfaces side by side
             │
             ├─> WC_electronic_stopping: the same CIM + fix electron/stopping (newer LAMMPS)
             ├─> SIM_random: random-orientation single impacts (W benchmark, then WC)
             └─> W_SBE / WC_SBE: surface binding energies
```

Conventions shared by all decks and scripts:

| | WC decks | W decks |
|---|---|---|
| atom types | W 1, C 2, He 3 | W 1, He 2 |
| substrate T | 614.6 K (0.2 T<sub>m</sub>) | 739 K |
| lattice | a = 2.9295 Å, c = 2.8242 Å at 614.6 K | 3.165 Å (EAM4) |

Further shared conventions:
- **Units.** `metal` units throughout.
- **CIM cycle.** 20 ps from one insertion to the next, so the flux is defined.
- **He launch.** At 30° from the surface normal.
- **Impact area.** The central 2/3 × 2/3 of the cell. Fluence = cycles / area.

## Requirements

**LAMMPS** (decks written for 29Oct2020; also run with 29Sep2021):
- `MANYBODY` (tersoff/zbl, eam/fs);
- core `pair table`, `fix dt/reset`;
- `EXTRA-FIX` for `fix electron/stopping`, which needs a newer LAMMPS. `WC_electronic_stopping/build_lammps.sh` builds one.

**Slurm** for the job scripts. Edit the lines marked `EDIT`: account, partition, cores, modules.

**Python 3.6+**:

| Package | Needed for |
|---|---|
| numpy, matplotlib | all analysis and plotting scripts |
| pandas | `analyze_rand.py`, `plot_erosion_compare.py`, `make_nov2_figures.py` (faster reading in `he_cluster_movies.py`) |
| scipy | `he_cluster_movies.py` |
| Pillow | GIFs |
| ovito (`pip install ovito`) | `ovito_insertion_movies.py`, `He_BE/scripts/renders.py`. Run headless with `QT_QPA_PLATFORM=offscreen`; `renders.py` also needs a Vulkan driver |
| jupyter, nbformat, nbclient | the notebooks and rebuilding them |
| *optional:* pysrim (+ SRIM) | `elstop.py basis --backend srim` |
| *optional:* LAMMPS Python module | `SIM_random/rand_geom.py` (otherwise it reads `rand_angles_lammps.csv`) |

`pip install -r requirements.txt` installs everything.

## Scripts and notebooks

Each analysis script has a notebook next to it (`script.py` → `script.ipynb`). The notebook holds the same code split into cells, and its figures, renders and GIFs are stored inside it, so GitHub shows them without running anything.
- **Run them.** Set the settings cell at the top (the same options as the command line), then *Run All*.
- **Where the shown outputs come from.**
  - Self-contained scripts (schematics, orientation GIFs, stopping tables) were run when the notebook was built.
  - Scripts that need simulation dumps (OVITO, He clusters) show outputs from a WC basal C-terminated 200 eV CIM run.
  - The erosion plots were built from the erosion CSVs of the W(001), WC basal and prism runs.
- **The `.py` is the source.** On a cluster, run the script, not the notebook. After editing a script, rebuild its notebook with `python tools/build_notebooks.py <name>`.

## Quick start (WC CIM on a Slurm cluster)

```bash
mkdir -p run/basal_C/100eV && cd run/basal_C/100eV
cp ../../../potentials/{WCH.tersoff,WHe.table,CHe.table,HeHe.table} .
cp ../../../WC_CIM/in_basal_092926.step1 ../../../WC_CIM/run.slurm .
cp ../../../WC_CIM/in_CIM_092926.step2 in.step2      # set 'variable e index 100' in it
lmp -in in_basal_092926.step1                        # slab -> wche0.restart (run once, e.g. as a job)
sbatch run.slurm                                     # chained jobs until NIMP cycles
python ../../../post_processing/ovito_insertion_movies.py --surface basal .
```
